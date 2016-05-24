

// Copyright 2016 UTK JICS AACE
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.


/*  $Id: blastn_app.cpp 461339 2015-03-09 18:07:37Z ivanov $
 * ===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
 *  the author's official duties as a United States Government employee and
 *  thus cannot be copyrighted.  This software/database is freely available
 *  to the public for use. The National Library of Medicine and the U.S.
 *  Government have not placed any restriction on its use or reproduction.
 *
 *  Although all reasonable efforts have been taken to ensure the accuracy
 *  and reliability of the software and data, the NLM and the U.S.
 *  Government do not and cannot warrant the performance or results that
 *  may be obtained by using this software or data. The NLM and the U.S.
 *  Government disclaim all warranties, express or implied, including
 *  warranties of performance, merchantability or fitness for any particular
 *  purpose.
 *
 *  Please cite the author in any work or product based on this material.
 *
 * ===========================================================================
 *
 * Authors:  Christiam Camacho
 *
 */

/** @file blastn_app.cpp
 * BLASTN command line application
 */

#ifndef SKIP_DOXYGEN_PROCESSING
static char const rcsid[] = 
	"$Id: blastn_app.cpp 461339 2015-03-09 18:07:37Z ivanov $";
#endif /* SKIP_DOXYGEN_PROCESSING */

#include <mpi.h>

#include <string>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <stdint.h>

#include <ncbi_pch.hpp>
#include <corelib/ncbiapp.hpp>
#include <algo/blast/api/local_blast.hpp>
#include <algo/blast/api/remote_blast.hpp>
#include <algo/blast/blastinput/blast_fasta_input.hpp>
#include <algo/blast/blastinput/blastn_args.hpp>
#include <algo/blast/api/objmgr_query_data.hpp>
#include <algo/blast/format/blast_format.hpp>
#include "blast_app_util.hpp"

#include <omp.h>

#ifndef SKIP_DOXYGEN_PROCESSING
USING_NCBI_SCOPE;
USING_SCOPE(blast);
USING_SCOPE(objects);
#endif

// HPC-BLAST::
// Define some global variables for MPI.

int g_rank, g_size;
int rep_rank, rep_group;

// Define globals for thread parameters
int num_thread_groups = 0;
int num_team_leaders = 0;
int num_searcher_threads = 0;
int num_threads_total = 0;
int num_ranks_per_group = 0;
int num_replication_groups = 0;
int num_ranks_total = 0;
std::string num_threads_str;
std::string db_str;
std::string query_str;
std::string out_str;


inline bool file_exists(const std::string& name)
{
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

struct membuf : std::streambuf
{
  membuf(char *buf, char *end) { this->setg(buf,buf,end); }
};

class CBlastnApp : public CNcbiApplication
{
public:
    /** @inheritDoc */
    CBlastnApp() {
        CRef<CVersion> version(new CVersion());
        version->SetVersionInfo(new CBlastVersion());
        SetFullVersion(version);
    }
private:
    /** @inheritDoc */
    virtual void Init();
    /** @inheritDoc */
    virtual int Run();

    /// This application's command line args
    CRef<CBlastnAppArgs> m_CmdLineArgs; 
};

void CBlastnApp::Init()
{
    // formulate command line arguments

    m_CmdLineArgs.Reset(new CBlastnAppArgs());

    // read the command line

    HideStdArgs(fHideLogfile | fHideConffile | fHideFullVersion | fHideXmlHelp | fHideDryRun);
    SetupArgDescriptions(m_CmdLineArgs->SetCommandLine());
}

int CBlastnApp::Run(void)
{
    int status = BLAST_EXIT_SUCCESS;

    bool working = true;

    // Declare the rank's buffer.
    char **buffer = 0;
    int *bsize = 0;

    buffer = (char**) malloc( num_thread_groups * sizeof(char*) );

    if ( buffer == 0 )
      return -1;

    for ( int i=0; i < num_thread_groups; ++i )
      buffer[i] = 0;

    bsize = (int*) malloc( num_thread_groups * sizeof(int) );

    if ( bsize == 0 )
      return -1;

    // Write out the global vars.
    std::cout<<"COMM_WORLD has "<<g_size<<" ranks."<<std::endl;
    std::cout<<"My rank is "<<g_rank<<std::endl;
    std::cout<<"My rep_group is "<<rep_group<<" and rep_rank is "<<rep_rank<<std::endl;
    std::cout<<"num_thread_groups = "<<num_thread_groups<<std::endl;
    std::cout<<"num_team_leaders = "<<num_team_leaders<<std::endl;
    std::cout<<"num_searcher_threads = "<<num_searcher_threads<<std::endl;
    std::cout<<"num_ranks_per_group = "<<num_ranks_per_group<<std::endl;
    std::cout<<"num_replication_groups = "<<num_replication_groups<<std::endl;

    for ( int i=0; i < num_thread_groups; ++i )
      {
	// Get the file name.
	std::stringstream qfile_name;
	qfile_name << query_str << "." <<rep_group<< "." <<i;

	// Open a stream to the file.
	std::ifstream ifs(qfile_name.str().c_str(), std::ifstream::binary );

	// Get the size of the file.
	ifs.seekg(0,ifs.end);
	bsize[i] = ifs.tellg();
	ifs.seekg(0,ifs.beg);  // Go back to the beginning.
	
	// Allocate the buffer.
	buffer[i] = (char*) malloc( bsize[i] * sizeof(char) );
	if ( buffer[i] == 0 )
	  return -1;

	// Read the file into the buffer.
	ifs.read( buffer[i], bsize[i] );

	// Close the stream.
	ifs.close();

	std::cout<<"Read in "<<bsize[i]<<" bytes from file "<<qfile_name.str()<<std::endl;
      }

    std::cout.flush();

#pragma omp parallel default(shared) num_threads( num_thread_groups*num_team_leaders )
    {
      int tid = omp_get_thread_num();

      // Compute the thread group I belong to.
      int thread_group_id = tid / num_team_leaders;
      int local_tid = tid % num_team_leaders;

      try {

        // Allow the fasta reader to complain on invalid sequence input
        SetDiagPostLevel(eDiag_Warning);
        SetDiagPostPrefix("hpc_blastn");

        /*** Get the BLAST options ***/
        const CArgs& args = GetArgs();

	CArgs my_args(GetArgs());

	// try to initialize a local copy.
	CRef<CBlastnAppArgs> my_CmdLineArgs;
	my_CmdLineArgs.Reset(new CBlastnAppArgs());

        CRef<CBlastOptionsHandle> opts_hndl;
        if(RecoverSearchStrategy(my_args, my_CmdLineArgs)){
        	opts_hndl.Reset(&*my_CmdLineArgs->SetOptionsForSavedStrategy(my_args));
        }
        else {
        	opts_hndl.Reset(&*my_CmdLineArgs->SetOptions(my_args));
        }
        const CBlastOptions& opt = opts_hndl->GetOptions();

        /*** Initialize the database/subject ***/
        CRef<CBlastDatabaseArgs> db_args(my_CmdLineArgs->GetBlastDatabaseArgs());
	
	CRef<CSearchDatabase> my_search_db = db_args->GetSearchDatabase();
	CRef<CSeqDB> my_seq_db = my_search_db->GetSeqDb();

	// Retreive the extents of the database from the CSeqDB object and compute the
	// actual extents threads in a given thread group will use.
	int db_num_seqs = my_seq_db->GetNumSeqs();

	int num_seqs_per_team_leader, my_oid_start, my_oid_stop;

	if ( db_num_seqs < num_team_leaders )
	  {
	    // Assign each available team leader at most a single database sequence.
	    if ( local_tid < db_num_seqs )
	      {
		my_oid_start = (local_tid + 1);
		my_oid_stop =  (local_tid + 2);
	      }
	    else
	      {
		my_oid_start = 0;
		my_oid_stop  = 0;
	      }
	  }
	else
	  {
	    num_seqs_per_team_leader = (db_num_seqs + num_team_leaders-1)/num_team_leaders;
	    my_oid_start =  local_tid      *num_seqs_per_team_leader;
	    my_oid_stop = ( local_tid + 1 )*num_seqs_per_team_leader;

	    if ( local_tid == (num_team_leaders-1) )
	      my_oid_stop = db_num_seqs;
	  }

	// Reset this thread's view into the database.
	my_seq_db->SetIterationRange(my_oid_start, my_oid_stop);
	
        CRef<CLocalDbAdapter> db_adapter;
        CRef<CScope> scope;
        InitializeSubject(db_args, opts_hndl, my_CmdLineArgs->ExecuteRemotely(),
                         db_adapter, scope);
        _ASSERT(db_adapter && scope);

	#pragma omp barrier

        /*** Get the query sequence(s) ***/
        CRef<CQueryOptionsArgs> query_opts = 
            my_CmdLineArgs->GetQueryOptionsArgs();
        SDataLoaderConfig dlconfig =
            InitializeQueryDataLoaderConfiguration(query_opts->QueryIsProtein(),
                                                   db_adapter);
        CBlastInputSourceConfig iconfig(dlconfig, query_opts->GetStrand(),
                                     query_opts->UseLowercaseMasks(),
                                     query_opts->GetParseDeflines(),
                                     query_opts->GetRange());

	membuf my_membuf( &(buffer[thread_group_id][0]), 
			  &(buffer[thread_group_id][ bsize[thread_group_id] - 1 ]) );

	auto_ptr<CNcbiIstream> my_new_input_stream;

	my_new_input_stream.reset(new CNcbiIstream(&my_membuf));

        if(IsIStreamEmpty(*my_new_input_stream)) {
           	ERR_POST(Warning << "Query is Empty!");
           	//return BLAST_EXIT_SUCCESS;
        }
        CBlastFastaInputSource fasta(*my_new_input_stream, iconfig);
        CBlastInput input(&fasta);

        // Initialize the megablast database index now so we can know whether an indexed search will be run.
        // This is only important for the reference in the report, but would be done anyway.
        if (opt.GetUseIndex() && !my_CmdLineArgs->ExecuteRemotely()) {
            CRef<CBlastOptions> my_options(&(opts_hndl->SetOptions()));
            CSetupFactory::InitializeMegablastDbIndex(my_options);
        }

	// Redirect thread's output file.
	string my_string (args[kArgOutput].AsString());
	string my_output;
	std::stringstream sstream;

	// Construct the thread's output file name.
	sstream << my_string << "." << rep_group << "." << rep_rank << "." << thread_group_id << "." << local_tid;
	my_output = sstream.str();
	
	auto_ptr<CNcbiOstream> my_new_output;
	my_new_output.reset(new CNcbiOfstream(my_output.c_str()));

        /*** Get the formatting options ***/
        CRef<CFormattingArgs> fmt_args(my_CmdLineArgs->GetFormattingArgs());
        CBlastFormat formatter(opt, *db_adapter,
                               fmt_args->GetFormattedOutputChoice(),
                               query_opts->GetParseDeflines(),
                               *my_new_output,
                               fmt_args->GetNumDescriptions(),
                               fmt_args->GetNumAlignments(),
                               *scope,
                               opt.GetMatrixName(),
                               fmt_args->ShowGis(),
                               fmt_args->DisplayHtmlOutput(),
                               opt.GetQueryGeneticCode(),
                               opt.GetDbGeneticCode(),
                               opt.GetSumStatisticsMode(),
                               my_CmdLineArgs->ExecuteRemotely(),
                               db_adapter->GetFilteringAlgorithm(),
                               fmt_args->GetCustomOutputFormatSpec(),
                               my_CmdLineArgs->GetTask() == "megablast",
                               opt.GetMBIndexLoaded());
                               
        formatter.SetQueryRange(query_opts->GetRange());
        formatter.SetLineLength(fmt_args->GetLineLength());
        if((fmt_args->GetFormattedOutputChoice() ==  CFormattingArgs::eXml2 ||
           fmt_args->GetFormattedOutputChoice() ==  CFormattingArgs::eJson)
           && args[kArgOutput].AsString() != "-")
        	formatter.SetBaseFile(args[kArgOutput].AsString());
        formatter.PrintProlog();

        /*** Process the input ***/
        CBatchSizeMixer mixer(SplitQuery_GetChunkSize(opt.GetProgram())-1000);
        int batch_size = my_CmdLineArgs->GetQueryBatchSize();
        if (batch_size) {
            input.SetBatchSize(batch_size);
        } else {
            Int8 total_len = formatter.GetDbTotalLength();
            if (total_len > 0) {
                /* the optimal hits per batch scales with total db size */
                mixer.SetTargetHits(total_len / 3000);
            }
            input.SetBatchSize(mixer.GetBatchSize());
        }

	#pragma omp single
	working = true;

	bool last_search = false;
	bool got_work = false;

	#pragma omp barrier

        for (; working; formatter.ResetScopeHistory()) {

	  my_seq_db->SetIterationRange(my_oid_start, my_oid_stop);
	  #pragma omp barrier

	  CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
	    
	  if ( !query_batch->Empty() )
	    got_work = true;
	  else
	    got_work = false;
          
          #pragma omp barrier

	  CRef<CSearchResultSet> results;

	  // No support for remote searches.
	  if ( got_work )  {
	    CRef<IQueryFactory> queries(new CObjMgr_QueryFactory(*query_batch));
	      
	    SaveSearchStrategy(my_args, my_CmdLineArgs, queries, opts_hndl);

	    CLocalBlast lcl_blast(queries, opts_hndl, db_adapter);
	    lcl_blast.SetNumberOfThreads(my_CmdLineArgs->GetNumThreads());
	    results = lcl_blast.Run();
	    if (!batch_size)
	      input.SetBatchSize(mixer.GetBatchSize(lcl_blast.GetNumExtensions()));
	  }

	  // Wait for all threads to finish. Then, reset all threads to 'see' the entire DB for reporting output.
          #pragma omp barrier
	  my_seq_db->SetIterationRange(0,db_num_seqs);
          #pragma omp barrier

	  // No support for archived format writes.
	  if ( got_work )  {
	    BlastFormatter_PreFetchSequenceData(*results, scope);
	    ITERATE(CSearchResultSet, result, *results)
	      {
		formatter.PrintOneResultSet(**result, query_batch);
	      }
	  }

	  #pragma omp single
	  working = false;

	  #pragma omp critical
	  {
	    if ( !input.End() )
	      working = true;
	  }
	  #pragma omp barrier

          #pragma omp flush(working)

	  #pragma omp barrier
	}

        formatter.PrintEpilog(opt);

        if (my_CmdLineArgs->ProduceDebugOutput()) {
            opts_hndl->GetOptions().DebugDumpText(NcbiCerr, "BLAST options", 1);
        }

    } CATCH_ALL(status)

	}   // end the parallel block

    return status;
}

#ifndef SKIP_DOXYGEN_PROCESSING
int main(int argc, char* argv[] /*, const char* envp[]*/)
{
  MPI_Init(&argc,&argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &g_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &g_size);

  int host_name_len;
  char host_name[MPI_MAX_PROCESSOR_NAME];
  
  MPI_Get_processor_name( &(host_name[0]), &host_name_len );

  double start_time, end_time, runtime, maxruntime;

  start_time = MPI_Wtime();

  // Redirect output to a file.
  char new_out_stream_name[100];
  char new_err_stream_name[100];

  sprintf(new_out_stream_name, "hpc_blast.stdout.%05d", g_rank);
  sprintf(new_err_stream_name, "hpc_blast.stderr.%05d", g_rank);

  freopen( new_out_stream_name, "w", stdout );
  freopen( new_err_stream_name, "w", stderr );
  
  omp_set_nested(1);

  // Root rank (rank 0) does some checking on file existence and parameter correctness.
  int state = 1;    // Variable to indicate the state of the job. Assume its good until proven bad.
  int i = 0;        // Loop counter over argc.
  int j = 0;         // Additional loop counter.
  
  int db_index = 0; // The index in argv where the database name appears.
  int out_index = 0; //The index in argv where the output name is.

  char *db_name_for_rank = 0;
  char *output_file_name = 0;

  // Array for packing information to send all data in a single call.
  int bcast_data[5];

  if ( g_rank == 0 )
    {
      // Tasks:
      // 1) Parse out db name, query name, and num_threads parameters from command line.
      // 2) Check the existence and read the job parameters file.
      // 3) Sanity check on the values.
      // 4) Verify the existence of all files implied by the job parameters.
      // 5) Broadcast the job parameters to all ranks, or tell them to quit.

      // Task 1 : Parse out database name, query filename, and number of threads for the search engine from the command line.
      for ( i=0; i < argc; ++i )
	{
	  if ( (strcmp("-db",argv[i]) == 0 ) && (i+1) < argc )
	    {
	      db_str.assign(argv[i+1]);
	      db_index = i+1;
	      break;
	    }
	}

      if ( i == argc )  // didn't find a db option.
	{
	  std::cerr<<"FATAL ERROR: '-db' was not found on the command line."<<std::endl;
	  std::cerr.flush();
	  state = -1;
	}

      for ( i=0; i < argc; ++i )
	{
	  if ( (strcmp("-query",argv[i]) == 0 ) && (i+1) < argc )
	    {
	      query_str.assign(argv[i+1]);
	      break;
	    }
	}

      if ( i == argc )  // didn't find a query option.
	{
	  std::cerr<<"FATAL ERROR: '-query' was not found on the command line."<<std::endl;
	  std::cerr.flush();
	  state = -1;
	}

      for ( i=0; i < argc; ++i )
	{
	  if ( (strcmp("-num_threads",argv[i]) == 0 ) && (i+1) < argc )
	    {
	      num_threads_str.assign(argv[i+1]);
	      num_searcher_threads = atoi( argv[i+1] );
	      break;
	    }
	}

      if ( i == argc )  // didn't find a num_threads option.
	{
	  // Silence warnings.
	  //std::cerr<<"WARNING: '-num_threads' was not found on the command line. Assuming a value of 1."<<std::endl;
	  //std::cerr.flush();
	  
	  num_threads_str.assign("1");
	  num_searcher_threads = 1;
	}

      for ( i=0; i < argc; ++i )
	{
	  if ( (strcmp("-out",argv[i]) == 0 ) && (i+1) < argc )
	    {
	      out_str.assign(argv[i+1]);
	      out_index = i+1;
	      break;
	    }
	}

      if ( i == argc ) // didn't find the output file option
	{
	  std::cerr<<"FATAL ERROR: '-out' was not specified on the command line."<<std::endl;
	  std::cerr.flush();
	  state=-1;
	}

      // Task 2 : Check the existence of and read in the job parameters file.
      if ( file_exists("job_params") )
	{
	  std::ifstream job_input("job_params");
	  job_input >> num_thread_groups >> num_team_leaders >> num_ranks_per_group >> num_replication_groups;
	  job_input.close();
	}
      else
	{
	  state = -1;
	  std::cerr<<"FATAL ERROR: job_params file not found!"<<std::endl;
	  std::cerr.flush();
	}

      // Task 3 : Sanity check on the values.
      if ( state == 1 )
	{
	  // Check num_searcher_threads.
	  if ( num_searcher_threads < 1 )
	    {
	      std::cerr<<"FATAL ERROR: num_searcher_threads must be positive."<<std::endl;
	      std::cerr<<"  num_searcher_threads = "<<num_searcher_threads<<std::endl;
	      std::cerr.flush();
	      state = -1;
	    }
	  else if ( num_searcher_threads > 100 ) // Max db partitions by the search engine.
	    {
	      std::cerr<<"FATAL ERROR: The number of searcher threads for the search engine exceeds the number of"<<std::endl;
	      std::cerr<<"             database partitions the engine supplies."<<std::endl;
	      std::cerr.flush();
	      state = -1;
	    }

	  // Check num_thread_groups.
	  if ( num_thread_groups < 1 )
	    {
	      std::cerr<<"FATAL ERROR: num_thread_groups must be positive."<<std::endl;
	      std::cerr<<"  num_thread_groups = "<<num_thread_groups<<std::endl;
	      std::cerr.flush();
	      state = -1;
	    }
	  else if ( num_thread_groups > 240 && num_thread_groups < 361 )
	    {
	      std::cerr<<"WARNING: num_thread_groups will oversubscribe resources, but not too high."<<std::endl;
	      std::cerr<<"  240 < (num_thread_groups="<<num_thread_groups<<") <= 360"<<std::endl;
	      std::cerr.flush();
	    }
	  else if ( num_thread_groups > 360 )
	    {
	      std::cerr<<"FATAL ERROR: num_thread_groups is too large."<<std::endl;
	      std::cerr<<"  360 < num_thread_groups="<<num_thread_groups<<std::endl;
	      std::cerr.flush();
	      state = -1;
	    }

	  // Check num_team_leaders.
	  if ( num_team_leaders < 1 )
	    {
	      std::cerr<<"FATAL ERROR: num_team_leaders must be positive."<<std::endl;
	      std::cerr<<"  num_team_leaders = "<<num_team_leaders<<std::endl;
	      std::cerr.flush();
	      state = -1;
	    }
	  else if ( num_team_leaders > 240 && num_team_leaders < 361 )
	    {
	      std::cerr<<"WARNING: num_team_leaders will oversubscribe resources, but not too high."<<std::endl;
	      std::cerr<<"  240 < (num_team_leaders="<<num_team_leaders<<") <= 360"<<std::endl;
	      std::cerr.flush();
	    }
	  else if ( num_team_leaders > 360 )
	    {
	      std::cerr<<"FATAL ERROR: num_team_leaders is too large."<<std::endl;
	      std::cerr<<"  360 < num_team_leaders="<<num_team_leaders<<std::endl;
	      std::cerr.flush();
	      state = -1;
	    }
	  
	  // Check the total number of threads used by a rank.
	  num_threads_total = num_thread_groups * num_team_leaders * num_searcher_threads ;
	  if ( num_threads_total < 1 )
	    {
	      std::cerr<<"FATAL ERROR: num_threads_total must be positive."<<std::endl;
	      std::cerr<<"  num_threads_total = "<<num_threads_total<<std::endl;
	      std::cerr.flush();
	      state = -1;
	    }
	  else if ( num_threads_total > 240 && num_threads_total < 361 )
	    {
	      std::cerr<<"WARNING: num_threads_total will oversubscribe resources, but not too high."<<std::endl;
	      std::cerr<<"  240 < (num_threads_total="<<num_threads_total<<") <= 360"<<std::endl;
	      std::cerr.flush();
	    }
	  else if ( num_threads_total > 360 )
	    {
	      std::cerr<<"FATAL ERROR: num_threads_total is too large."<<std::endl;
	      std::cerr<<"  360 < num_threads_total="<<num_threads_total<<std::endl;
	      std::cerr.flush();
	      state = -1;
	    }

	  // Test MPI values.
	  if ( num_replication_groups < 1 )
	    {
	      std::cerr<<"FATAL ERROR: num_replication_groups must be positive."<<std::endl;
	      std::cerr<<"  num_replication_groups = "<<num_replication_groups<<std::endl;
	      std::cerr.flush();
	      state = -1;
	    }

	  if ( num_ranks_per_group < 1 )
	    {
	      std::cerr<<"FATAL ERROR: num_ranks_per_group must be positive."<<std::endl;
	      std::cerr<<"  num_ranks_per_group = "<<num_ranks_per_group<<std::endl;
	      std::cerr.flush();
	      state = -1;
	    }

	  num_ranks_total = num_replication_groups * num_ranks_per_group;

	  if ( num_ranks_total != g_size )
	    {
	      std::cerr<<"FATAL ERROR: num_ranks_total requested by the job does not match number of ranks in COMM_WORLD."<<std::endl;
	      std::cerr<<"  num_ranks_total = "<<num_ranks_total<<std::endl;
	      std::cerr<<"  Size of MPI_COMM_WORLD = "<<g_size<<std::endl;
	      std::cerr.flush();
	      state = -1;
	    }

	} // End if valid state value for Task 3.

      // Task 4 : Verify existence of all files that the parameters imply should exist.
      //if ( state == 1 )
      if ( 0 )
	{
	  // First check database files.
	  std::string file_name;

	  if ( num_ranks_per_group == 1 )
	    {
	      file_name.assign(db_str);

	      file_name += ".nin";
	      if ( ! file_exists(file_name) )
		{
		  std::cerr<<"FATAL ERROR: The database file does not exist: "<<file_name<<std::endl;
		  std::cerr.flush();
		  state = -1;
		}

	      file_name.assign(db_str);
	      file_name += ".nsq";
	      if ( ! file_exists(file_name) )
		{
		  std::cerr<<"FATAL ERROR: The database file does not exist: "<<file_name<<std::endl;
		  std::cerr.flush();
		  state = -1;
		}

	      file_name.assign(db_str);
	      file_name += ".nhr";
	      if ( ! file_exists(file_name) )
		{
		  std::cerr<<"FATAL ERROR: The database file does not exist: "<<file_name<<std::endl;
		  std::cerr.flush();
		  state = -1;
		}
	    }
	  else  // multiple partitions.
	    {
	      for ( i=0; i < num_ranks_per_group; ++i )
		{
		  std::stringstream sstm;

		  sstm << db_str;

		  if ( i < 10 )
		    sstm << ".0" << i;
		  else
		    sstm << "." << i;

		  file_name.assign( sstm.str() + ".nin" );
		  if ( ! file_exists(file_name) )
		    {
		      std::cerr<<"FATAL ERROR: The database file does not exist: "<<file_name<<std::endl;
		      std::cerr.flush();
		      state = -1;
		    }

		  file_name.assign( sstm.str() + ".nsq" );
		  if ( ! file_exists(file_name) )
		    {
		      std::cerr<<"FATAL ERROR: The database file does not exist: "<<file_name<<std::endl;
		      std::cerr.flush();
		      state = -1;
		    }

		  file_name.assign( sstm.str() + ".nhr" );
		  if ( ! file_exists(file_name) )
		    {
		      std::cerr<<"FATAL ERROR: The database file does not exist: "<<file_name<<std::endl;
		      std::cerr.flush();
		      state = -1;
		    }
		}
	    }

	  // Now check query file. Assume we will always partition the queries ( or at least pack them so they will always be in the given file naming convention.

	  for ( i=0; i < num_replication_groups; ++i )
	    {
	      
	      for ( j=0; j < num_thread_groups; ++j )
		{
		  std::stringstream sstm;

		  sstm << query_str << "." << i << "." << j;

		  file_name = sstm.str();

		  if ( ! file_exists(file_name) )
		    {
		      std::cerr<<"FATAL ERROR: The query file does not exist: "<<file_name<<std::endl;
		      std::cerr.flush();
		      state = -1;
		    }
		}
	      
	    }
	  
	} // End if valid state value for Task 4.

    } // End if root rank.

  // Task 5 : Broadcast the state of the job to all ranks. If a valid job, broadcast parameters to all ranks.

  if ( g_rank == 0 ) 
    {
      std::cout<<"Root is about to broadcast the state variable: "<<state<<std::endl;
      std::cout.flush();
    }
  
  MPI_Bcast ( &state, 1, MPI_INT, 0, MPI_COMM_WORLD );

  if ( state == -1 )
    {
      MPI_Finalize();
      return 1;
    }
  else
    {
      // Root packs the data for broadcast.
      if ( g_rank == 0 )
	{
	  bcast_data[0] = num_thread_groups;
	  bcast_data[1] = num_team_leaders;
	  bcast_data[2] = num_searcher_threads;
	  bcast_data[3] = num_ranks_per_group;
	  bcast_data[4] = num_replication_groups;

	  std::cout<<"Sending bcast_data[0] = num_thread_groups      = "<<num_thread_groups<<std::endl;
	  std::cout<<"Sending bcast_data[0] = num_team_leaders       = "<<num_team_leaders<<std::endl;
	  std::cout<<"Sending bcast_data[0] = num_searcher_threads   = "<<num_searcher_threads<<std::endl;
	  std::cout<<"Sending bcast_data[0] = num_ranks_per_group    = "<<num_ranks_per_group<<std::endl;
	  std::cout<<"Sending bcast_data[0] = num_replication_groups = "<<num_replication_groups<<std::endl;
	  std::cout.flush();
	}

      MPI_Bcast ( bcast_data, 5, MPI_INT, 0, MPI_COMM_WORLD );

      if ( g_rank != 0 )
	{
	  // Unpack the sent data.
	  num_thread_groups      = bcast_data[0];
	  num_team_leaders       = bcast_data[1];
	  num_searcher_threads   = bcast_data[2];
	  num_ranks_per_group    = bcast_data[3];
	  num_replication_groups = bcast_data[4];

	  std::cout<<"Received bcast_data[0] = num_thread_groups      = "<<num_thread_groups<<std::endl;
	  std::cout<<"Received bcast_data[0] = num_team_leaders       = "<<num_team_leaders<<std::endl;
	  std::cout<<"Received bcast_data[0] = num_searcher_threads   = "<<num_searcher_threads<<std::endl;
	  std::cout<<"Received bcast_data[0] = num_ranks_per_group    = "<<num_ranks_per_group<<std::endl;
	  std::cout<<"Received bcast_data[0] = num_replication_groups = "<<num_replication_groups<<std::endl;
	  std::cout.flush();

	  // Also get the database name and query file name from the command line. We know it is there as implied by the state.
	  for ( i=0; i < argc; ++i )
	    {
	      if ( (strcmp("-db",argv[i]) == 0 ) && (i+1) < argc )
		{
		  db_str.assign(argv[i+1]);
		  db_index = i+1;
		  break;
		}
	    }

	  for ( i=0; i < argc; ++i )
	    {
	      if ( (strcmp("-query",argv[i]) == 0 ) && (i+1) < argc )
		{
		  query_str.assign(argv[i+1]);
		  break;
		}
	    }

	  for ( i=0; i < argc; ++i )
	    {
	      if ( (strcmp("-out",argv[i]) == 0 ) && (i+1) < argc )
		{
		  out_str.assign(argv[i+1]);
		  out_index= i+1;
		  break;
		}
	    }
	  
	}
    }

  MPI_Barrier( MPI_COMM_WORLD );

  // Establish the rank's ID in the replication group.
  rep_rank = g_rank % num_ranks_per_group;          // Which db partition to read.
  rep_group = g_rank / num_ranks_per_group;         // Which query partition to read.

  // If the database is partitioned, we will have to change the process' db file name so that it will open the appropriate file.
  if ( num_ranks_per_group > 1 )
    {
      db_name_for_rank = new char[ strlen(argv[db_index]) + 10 ];

      if ( rep_rank < 10 )
	sprintf(db_name_for_rank,"%s.0%d",argv[db_index],rep_rank);
      else
	sprintf(db_name_for_rank,"%s.%d",argv[db_index],rep_rank);

      argv[db_index] = db_name_for_rank;
    }

  // Allow for each rank to have its own output file. This is purely for the sake of the NCBI toolkit which will attempt to open it for writing later (which is overridden in the threaded section).
  if ( g_size > 1 )
    {
      output_file_name = new char[ strlen(argv[out_index]) + 10];

      sprintf(output_file_name,"%s.%d",argv[out_index],g_rank);

      argv[out_index] = output_file_name;
    }

  std::cout<<"Rank "<<g_rank<<" is running on "<<host_name<<std::endl;
  std::cout<<"Rank "<<g_rank<<" : my database is "<<argv[db_index]<<std::endl;

  // Start the BLAST search.

  int blast_return = CBlastnApp().AppMain(argc, argv, 0, eDS_Default, 0);

  end_time = MPI_Wtime();

  runtime = end_time - start_time;

  std::cout<<"Runtime is "<<runtime<<std::endl;
  std::cout.flush();

  // Catch all nodes at the end.
  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Reduce( &runtime, &maxruntime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );

  if ( g_rank == 0 )
    std::cout<<"Job runtime is "<<maxruntime<<std::endl;

  MPI_Finalize();

  if ( db_name_for_rank != 0 )
    delete [] db_name_for_rank;

  if ( output_file_name != 0 )
    delete [] output_file_name;

  return blast_return;
}
#endif /* SKIP_DOXYGEN_PROCESSING */
