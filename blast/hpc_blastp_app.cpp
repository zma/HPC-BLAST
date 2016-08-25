

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


/*  $Id: blastp_app.cpp 461340 2015-03-09 18:08:15Z ivanov $
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

/** @file blastp_app.cpp
 * BLASTP command line application
 */

#ifndef SKIP_DOXYGEN_PROCESSING
static char const rcsid[] =
	"$Id: blastp_app.cpp 461340 2015-03-09 18:08:15Z ivanov $";
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
#include <algo/blast/blastinput/blastp_args.hpp>
#include <algo/blast/api/objmgr_query_data.hpp>
#include <algo/blast/format/blast_format.hpp>
#include "blast_app_util.hpp"

#include "../../algo/blast/api/blast_aux_priv.hpp"
#include "../../algo/blast/api/blast_seqalign.hpp"
#include "../../algo/blast/core/blast_psi_priv.h"//for serialization
#include <objects/seq/Seq_annot.hpp>//for output of alignments ceb
#include <omp.h>
#include <fstream>

#ifndef SKIP_DOXYGEN_PROCESSING
USING_NCBI_SCOPE;
USING_SCOPE(blast);
USING_SCOPE(objects);
#endif

#define USE_RESTART 1
#define CLEAN_UP_OBJECTS 1
#define USE_SHALLOW_COPY 0

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

// Combine HSPResults within a thread group (across team leaders)
// Input is array of serialized HSPResults objects
BlastHSPResults* gbhrpfncCombineTeamLeaderBlastHSPResults  (char ** HSPResultsBuffers,
							    size_t * buffer_size,
							    int tid,
							    int num_team_leaders );

// Combine BlastHSPResults across team leaders
// Input is array of HSPResults objects
BlastHSPResults* gbhrpfncCombineTeamLeaderBlastHSPResults ( BlastHSPResults** abhrppBlastHSPResults,
							    int aintThreadID,
							    int aintNumberOfTeamLeaders );

// Output objects to file
void gvoifncDumpBlastHSPResults (BlastHSPResults *abhrpBlastHSPResults);

void gvoifncDumpBlastQueryInfo (BlastQueryInfo *abqipBlastQueryInfo);

void gvoifncDumpBlast_KarlinBlk (Blast_KarlinBlk **abkbppBlast_KarlinBlk, 
				 int aintNumberOfContexts);

void gvoifncDumpBlast_KarlinBlk1 (Blast_KarlinBlk *abkbpBlast_KarlinBlk);

void gvoifncDumpSBlastScoreMatrix(SBlastScoreMatrix *absmpSBlastScoreMatrix, 
				  int alphabet_size);

void gvoifncDumpBlastScoreBlk(BlastScoreBlk *absbpBlastScoreBlk);

//Combine objects across thread groups
BlastScoreBlk *gbsbpfncCombineThreadGroupBlastScoreBlk(int aintNumberOfThreadGroups, 
						       BlastScoreBlk **ScoreBlkGroup);

BlastHSPResults *gbhrpfncCombineThreadGroupBlastHSPResults(int aintNumberOfThreadGroups, 
							   BlastHSPResults **abhrppBlastHSPResults); 

BlastQueryInfo *gbqipfncCombineThreadGroupBlastQueryInfo(int aintNumberOfThreadGroups, 
							 BlastQueryInfo **abqippBlastQueryInfo); 

CRef<CBlastQueryVector> gbqvfncCombineReplicationGroupQueries(int aintNumberOfThreadGroups,
							      CRef<CBlastQueryVector>* abqvpQueryBatchGroup);					    

BlastHSPResults* gbhrpfncDuplicateBlastHSPResults( BlastHSPResults* HSPResultsArray);
BlastScoreBlk *gbsbpfncDuplicateBlastScoreBlk(BlastScoreBlk *ScoreBlkGroup);


CRef<CBlastQueryVector> gbqvfncCombineReplicationGroupQueries(int aintNumberOfThreadGroups,
							      char** buffer,
							      int gbl_bsize,
							      int* lintpBatchStartPosition,
							      int* lintpBatchEndPosition,
							      CBlastInputSourceConfig* iconfig,
							      int QueryBatchSize,
							      int gbl_query_count,
							      CRef<CScope> gbl_scope);

//Create search results set from objects
CRef<CSearchResultSet> gsrsfncCollectResults(BlastHSPResults* hsp_results,
					     BlastQueryInfo* m_QueryInfo,
					     BlastScoreBlk* m_ScoreBlk,
					     CRef<IQueryFactory> m_QueryFactory,
					     const CBlastOptions* m_Options,
					     CRef<IBlastSeqInfoSrc> m_SeqInfoSrc);


BlastHSPList* BlastHSPListDuplicate(const BlastHSPList *bhspl);







//==========================================


class CBlastpApp : public CNcbiApplication
{
public:
   /** @inheritDoc */
   CBlastpApp() {
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
   CRef<CBlastpAppArgs> m_CmdLineArgs;

   //ceb array of pointers for each threads objects necessary to build output 
   BlastHSPResults ** tHSPResults = NULL;
   BlastQueryInfo ** tQueryInfo = NULL;
   BlastScoreBlk ** tBlastScoreBlk = NULL;

   BlastHSPResults ** HSPResultsGroup = NULL;
   BlastQueryInfo ** QueryInfoGroup = NULL;
   BlastScoreBlk ** BlastScoreBlkGroup = NULL;
   CRef<CBlastQueryVector>* QueryBatchGroup = NULL;
};

//ceb function to test for existence of file
bool fexists(const char *filename)
{
  ifstream ifile(filename);
  return ifile.good();
};


//ceb This function searches a file for n occurances of a given phrase.
//    It returns the file position of the line following the 
//    nth phrase location. If n occurances of the phrase are not found,
//    the NULL value is returned.
long get_last_good_position(string filename, string phrase, int nqueries)
{
    std::ifstream infile(filename,std::fstream::in);
    if ( !infile.is_open() ) {
        cerr << "Error while opening file." << endl;
        exit( EXIT_FAILURE );
    }    

    int count=0;
    string line;   
    //printf("nqueries=%d\n",nqueries);
    //printf("searching output file for phrase: %s\n",phrase.c_str());
    while( !infile.eof() ) 
    {
       std::getline(infile,line);  
       //printf("current line: %s\n",line.c_str());
       if( strstr(line.c_str(), phrase.c_str()) != NULL)
       {
           count++;
           //printf("%dth phrase %s found\n",count,phrase.c_str());
           if (count==nqueries) 
           {
               long pos = infile.tellg();
               //printf("last completed report ended at pos %d\n",pos);
               infile.close();
               //return current position in the file
               return pos;
           }
       }
    }   
    //printf("%dth phrase %s found\n",count,phrase.c_str());
    infile.close();
    return NULL;
}




void CBlastpApp::Init()
{
    // formulate command line arguments

    m_CmdLineArgs.Reset(new CBlastpAppArgs());

    // read the command line

    HideStdArgs(fHideLogfile | fHideConffile | fHideFullVersion | fHideXmlHelp | fHideDryRun);
    SetupArgDescriptions(m_CmdLineArgs->SetCommandLine());
}


int CBlastpApp::Run(void)
{
    int status = BLAST_EXIT_SUCCESS;

    bool working = true;

    // Declare the rank's buffer.
    char **buffer = 0;
    int *bsize = 0;

    int n_batches_precomputed=0;


    buffer = (char**) malloc( num_thread_groups * sizeof(char*) );
    if ( buffer == 0 ){ return -1; }
    for ( int i=0; i < num_thread_groups; ++i ){ buffer[i] = 0; }

    bsize = (int*) malloc( num_thread_groups * sizeof(int) );
    if ( bsize == 0 ){ return -1; }

    // Query buffer end posititon for each query by thread group
    int **lintppThreadGroupQueryBufferEndPosition=0;
    lintppThreadGroupQueryBufferEndPosition=(int**)malloc(num_thread_groups*sizeof(int*));
    if (lintppThreadGroupQueryBufferEndPosition==0) {return -1;}
    
    // Query buffer batch start position by thread group
    int *lintpBatchStartPosition=0;
    lintpBatchStartPosition=(int*)malloc(num_thread_groups*sizeof(int));
    if (lintpBatchStartPosition==0) {return -1;}
    
    // Query buffer batch end position by thread group
    int *lintpBatchEndPosition=0;
    lintpBatchEndPosition=(int*)malloc(num_thread_groups*sizeof(int));
    if (lintpBatchEndPosition==0) {return -1;}


    // Write out the global vars.
    std::cout<<"COMM_WORLD has "<<g_size<<" ranks."<<std::endl;
    std::cout<<"My rank is "<<g_rank<<std::endl;
    std::cout<<"My rep_group is "<<rep_group<<" and rep_rank is "<<rep_rank<<std::endl;
    std::cout<<"num_thread_groups = "<<num_thread_groups<<std::endl;
    std::cout<<"num_team_leaders = "<<num_team_leaders<<std::endl;
    std::cout<<"num_searcher_threads = "<<num_searcher_threads<<std::endl;
    std::cout<<"num_ranks_per_group = "<<num_ranks_per_group<<std::endl;
    std::cout<<"num_replication_groups = "<<num_replication_groups<<std::endl;

    int gbl_bsize=0;//ceb

    //Load query file into buffer
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

       if ( buffer[i] == 0 ){ return -1; }

       // Read the file into the buffer.
       ifs.read( buffer[i], bsize[i] );

       // Close the stream.
       ifs.close();

       std::cout<<"Read in "<<bsize[i]<<" bytes from file "<<qfile_name.str()<<std::endl;

       //std::cout<<"\n init buffer["<<i<<"]:\n"<<buffer[i]<<std::endl<<std::endl; std::cout.flush();

       gbl_bsize+=bsize[i];


       //Set up variables for use in tracking postion of input buffers for output merging 
       // Find number of queries in this thread group
       int lintNumberOfThreadGroupQueries=0;
       for (int j=0; j<bsize[i]; ++j) {
	 //a query entry is denoted by the ">gi" sequence of characters
	 if (buffer[i][j]=='>') {  
           //a query entry can have multiple >gi entries prior to the listing of the sequence
	   // we need to wait till we reach a return character before we can count the query
	   while(buffer[i][j]!='\n' && j<bsize[i]){
	     ++j;
	   }
	   lintNumberOfThreadGroupQueries++; 
	 }
       }

       //std::cout<<"lintNumberOfThreadGroupQueries="
       //<<lintNumberOfThreadGroupQueries<<std::endl;std::cout.flush();

       // Allocate space for query buffer end position by thread group
       lintppThreadGroupQueryBufferEndPosition[i]=(int*)malloc(lintNumberOfThreadGroupQueries*sizeof(int));
       
       // Set query buffer end position by thread group
       int lintCurrentQuery=-1;
       for (int j=0; j<bsize[i]; ++j) {
	 //a query entry is denoted by the ">gi" sequence of characters
	 if (buffer[i][j]=='>'){
	   //the next query must be preceded by one or more return characters
	   
	   while(j<bsize[i])
	   {
	     ++j;
	     if(buffer[i][j]=='\n' && buffer[i][j+1]=='>'){break;}
	   }

	   lintCurrentQuery++;
	   if (j>1) {
	     lintppThreadGroupQueryBufferEndPosition[i][lintCurrentQuery]=j;
	     //std::cout<<"lintppThreadGroupQueryBufferEndPosition["<<i<<"]["<<lintCurrentQuery<<"]="<<lintppThreadGroupQueryBufferEndPosition[i][lintCurrentQuery-1]<<std::endl;std::cout.flush();
	   }
	 }

       }

       // Set final query buffer end position 
       lintppThreadGroupQueryBufferEndPosition[i][lintNumberOfThreadGroupQueries-1]=bsize[i]-1;
       //std::cout<<"lintppThreadGroupQueryBufferEndPosition["<<i<<"]["
       //<<lintNumberOfThreadGroupQueries-1<<"]="
       //<<lintppThreadGroupQueryBufferEndPosition[i][lintNumberOfThreadGroupQueries-1]
       //<<std::endl;std::cout.flush();
    }
    std::cout.flush();

    // Allocate the thread buffers for serialized data (only the first ptr, the threads will allocate for us)
    int number_of_threads = num_thread_groups*num_team_leaders ;

    // Pointers to thread objects used to consolidate output
    tHSPResults = new BlastHSPResults*[number_of_threads];
    tQueryInfo = new BlastQueryInfo*[number_of_threads];
    tBlastScoreBlk = new BlastScoreBlk*[number_of_threads];

    //variables to point to objects merged within each thread group
    HSPResultsGroup = new BlastHSPResults*[num_thread_groups];
    QueryInfoGroup = new BlastQueryInfo*[num_thread_groups];
    BlastScoreBlkGroup = new BlastScoreBlk*[num_thread_groups];
    QueryBatchGroup = new CRef<CBlastQueryVector>[num_thread_groups];



#pragma omp parallel default(shared) num_threads( num_thread_groups*num_team_leaders )
    { 

       int tid = omp_get_thread_num();
 
       // Compute the thread group I belong to.
       int thread_group_id = tid / num_team_leaders;
       int local_tid = tid % num_team_leaders;

       // Running total of thread group queries emitted from GetNextSeqBatch
       int lintCurrentNumberOfThreadGroupQueries=0;
       
       // Calculate the next query buffer start position from the current query buffer end position (add 1) 
       int lintNextBatchStartPosition=0;

       //ceb
       int gbl_query_count=0;

       try 
       {
          // Allow the fasta reader to complain on invalid sequence input
          SetDiagPostLevel(eDiag_Warning);
          SetDiagPostPrefix("hpc_blastp");

          /*** Get the BLAST options ***/
          const CArgs& args = GetArgs();

	  CArgs my_args(GetArgs());

	  // try to initialize a local copy.
	  CRef<CBlastpAppArgs> my_CmdLineArgs;
	  my_CmdLineArgs.Reset(new CBlastpAppArgs());

          CRef<CBlastOptionsHandle> opts_hndl;
          if(RecoverSearchStrategy(my_args, my_CmdLineArgs)) 
          {
             opts_hndl.Reset(&*my_CmdLineArgs->SetOptionsForSavedStrategy(my_args));
          }
          else 
          {
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
	     if ( local_tid == (num_team_leaders-1) ){  my_oid_stop = db_num_seqs; }
	  }

	  // Reset this thread's view into the database.
	  my_seq_db->SetIterationRange(my_oid_start, my_oid_stop);

          CRef<CLocalDbAdapter> db_adapter;
          CRef<CScope> lcl_scope;

          InitializeSubject(db_args, opts_hndl, my_CmdLineArgs->ExecuteRemotely(), db_adapter, lcl_scope);

          _ASSERT(db_adapter && lcl_scope);

	  //create a global scope object for thread 0
          CRef<CScope> gbl_scope;
          if(tid==0) InitializeSubject(db_args, opts_hndl, my_CmdLineArgs->ExecuteRemotely(),db_adapter, gbl_scope);//ceb

	  #pragma omp barrier


          /*** Get the query sequence(s) ***/
          CRef<CQueryOptionsArgs> query_opts = my_CmdLineArgs->GetQueryOptionsArgs();

          SDataLoaderConfig dlconfig =
             InitializeQueryDataLoaderConfiguration(query_opts->QueryIsProtein(),
                                                    db_adapter);

          CBlastInputSourceConfig iconfig(dlconfig, query_opts->GetStrand(),
                                          query_opts->UseLowercaseMasks(),
                                          query_opts->GetParseDeflines(),
                                          query_opts->GetRange());


	  //Set up thread local input stream
	  //buffer[i] contains full query input from each file[i]
	  membuf lcl_query_membuf( &(buffer[thread_group_id][0]),
				   &(buffer[thread_group_id][ bsize[thread_group_id] - 1 ]) );
	  auto_ptr<CNcbiIstream> lcl_query_input_stream;
	  lcl_query_input_stream.reset(new CNcbiIstream(&lcl_query_membuf));
          if(IsIStreamEmpty(*lcl_query_input_stream)){ERR_POST(Warning << "Query is Empty!");}
          CBlastFastaInputSource lcl_fasta(*lcl_query_input_stream, iconfig);
          CBlastInput lcl_input(&lcl_fasta, my_CmdLineArgs->GetQueryBatchSize());

	  // Redirect thread's output file.
	  string my_string (args[kArgOutput].AsString());
	  string my_output;
	  std::stringstream sstream;

	  // Construct the thread's output file name.
	  sstream << my_string << "." << rep_group << "." << rep_rank << "." << thread_group_id << "." << local_tid;
	  my_output = sstream.str();

	  auto_ptr<CNcbiOstream> my_new_output;

#if USE_RESTART //Restart code


	  //ceb Check to see if restart file exists for this output file (set of queries)
	  //    If restart file does not exist, then we default to no restart.
	  //    If we fail gracefully at the end, we will create a restart file
	  //    Note: restart files correspond to thread output, so there may be restart files for some threads
	  //    and not for others.
	  bool restarted=false;//ceb
	  int my_batches_precomputed=0;
	  int query_count=0;
	  string my_restart;
	  auto_ptr<CNcbiFstream> restart_stream;
	  if(tid==0)
	    {
	      //ceb create restart file name. Identical to output file name but with .restart appended
	      sstream << ".restart";
	      my_restart = sstream.str();

	      //ceb Open stream for restart file
	      restart_stream.reset(new CNcbiFstream(my_restart.c_str(),std::ios::in));
	      if(restart_stream->good())
		{  //ceb If file exists, then open input stream and read count
		  *restart_stream >> my_batches_precomputed >> query_count;
		  //printf("Opening restart file for reading: %d batches with %d queries precomputed\n",n_batches_precomputed,query_count);
		  restart_stream->close();//close this file. We may recreate it at the end
		  if(my_batches_precomputed>0) restarted=true;
		}



              //ceb Setting up output file
              if(restarted && fexists(my_output.c_str()))
		{ //ceb If file exists, then open output stream and begin write at end of file
		  //printf("Opening previous output file\n");
		  //ceb locate pointer location of last good position based on number of batch queries successfully written
		  string phrase = "Effective search space used:";//phrase indicates end of query output  
		  long pos = get_last_good_position(my_output.c_str(), phrase, query_count);
		  //printf("restart file pos=%d\n",pos);
		  //ceb open file as output and move write pointer to position
		  //File must be read/write in order to alter the position of the initial write
		  CNcbiOfstream* my_output_stream = new CNcbiOfstream(my_output.c_str(),std::fstream::in | std::fstream::out);
		  //move pointer
		  if(pos!=NULL) 
		    {
		      my_output_stream->seekp(pos,ios_base::beg);
		    }
		  //assign pointer to stream
		  my_new_output.reset(my_output_stream);
		}
	      else
		{  //ceb else create new file and open for writing
		  //printf("Opening new output file\n");
		  my_new_output.reset(new CNcbiOfstream(my_output.c_str()));
		  //ceb if a restart file exists but the output file does not, we have to start from scratch
		  restarted=false;
		  my_batches_precomputed=0;
		}
	    } 


	  //make sure all threads get the number of precomputed batches
          //#pragma omp single
	    //if(n_batches_precomputed < my_batches_precomputed)
	  //{n_batches_precomputed = my_batches_precomputed;}
          if (tid==0){n_batches_precomputed = my_batches_precomputed;}
          #pragma omp flush(n_batches_precomputed)
          #pragma omp barrier


#else //no restart
	  if (tid==0){ my_new_output.reset(new CNcbiOfstream(my_output.c_str()));}
#endif

          /*** Get the formatting options ***/
          CRef<CFormattingArgs> fmt_args(my_CmdLineArgs->GetFormattingArgs());
	  CBlastFormat* formatter;
          if(tid==0)
	  {
            formatter= new CBlastFormat(opt, 
					*db_adapter,
					fmt_args->GetFormattedOutputChoice(),
					query_opts->GetParseDeflines(),
					*my_new_output,
					fmt_args->GetNumDescriptions(),
					fmt_args->GetNumAlignments(),
					*gbl_scope,
					opt.GetMatrixName(),
					fmt_args->ShowGis(),
					fmt_args->DisplayHtmlOutput(),
					opt.GetQueryGeneticCode(),
					opt.GetDbGeneticCode(),
					opt.GetSumStatisticsMode(),
					my_CmdLineArgs->ExecuteRemotely(),
					db_adapter->GetFilteringAlgorithm(),
					fmt_args->GetCustomOutputFormatSpec());
	    
            formatter->SetQueryRange(query_opts->GetRange());
            formatter->SetLineLength(fmt_args->GetLineLength());
            if((fmt_args->GetFormattedOutputChoice() ==  CFormattingArgs::eXml2 ||
              fmt_args->GetFormattedOutputChoice() ==  CFormattingArgs::eJson)
              && args[kArgOutput].AsString() != "-")
        	{formatter->SetBaseFile(args[kArgOutput].AsString());}
#if USE_RESTART
            //ceb Skip printing the file prologue if we are restarting from a previous solution
            if(!restarted){formatter->PrintProlog();}
#else
            formatter->PrintProlog();
#endif
	  }

	  #pragma omp single
	  {working = true;}

	  bool last_search = false;
	  bool got_work = false;

	  #pragma omp barrier
	  //int cntr=0;//ceb

          /*** Process the input ***/
	  //#if USE_RESTART
	  int batch_count=0;
	  //#endif
	  //if(tid==0)std::cout<<"checkpt 0.1 \n";std::cout.flush();

          for (; working; ) 
          {
	    //if(tid==0)std::cout<<"checkpt 0.2 \n";std::cout.flush();

	    //clear memory from last scope object
	    if(tid==0){formatter->ResetScopeHistory();}
	    lcl_scope->ResetDataAndHistory();//ceb

	    my_seq_db->SetIterationRange(my_oid_start, my_oid_stop);
	    #pragma omp barrier

            CRef<CBlastQueryVector> query_batch(lcl_input.GetNextSeqBatch(*lcl_scope));

	    //#if USE_RESTART
	    //ceb increment completed batch counter for log file
	    batch_count++;//ceb
	    //#endif

	    //if(tid==0)std::cout<<"checkpt 1 rank="<<g_rank<<" batch_count= "<<batch_count<<"\n";    std::cout.flush();

	    // Get number of queries in batch
	    int lintNumberOfQueriesInBatch=query_batch->Size();

	    // Only team leader 0 needs to do this
	    if (local_tid==0) {// If first batch

	      //std::cout<<"checkpt 1.1 local_tid="<<local_tid<<" query_batch->Size()="<<query_batch->Size()<<"\n";std::cout.flush();

	      if (lintCurrentNumberOfThreadGroupQueries==0) {// Start at 0 
		lintpBatchStartPosition[thread_group_id] = 0;
	      } 
	      else{
		lintpBatchStartPosition[thread_group_id] = lintNextBatchStartPosition;
	      }
	      //if(tid==0)std::cout<<"checkpt 1.2 \n";std::cout.flush();
	      // Update running total 
	      lintCurrentNumberOfThreadGroupQueries += lintNumberOfQueriesInBatch;  
	      // Set query buffer end posiiton 
	      lintpBatchEndPosition[thread_group_id] = lintppThreadGroupQueryBufferEndPosition[thread_group_id][lintCurrentNumberOfThreadGroupQueries-1];
	      // Get ready for next iteraton 
	      lintNextBatchStartPosition = lintppThreadGroupQueryBufferEndPosition[thread_group_id][lintCurrentNumberOfThreadGroupQueries-1]+1;
	      //std::cout<<"checkpt 1.2.1 local_tid="<<local_tid<<" lintNextBatchStartPosition="<<lintNextBatchStartPosition<<"\n";std::cout.flush();
	    }

	    //if(tid==0)std::cout<<"checkpt 1.3 \n";std::cout.flush();

#pragma omp barrier
	    //reset flag
	    got_work=false;

#if USE_RESTART
	    // Check to see if this batch has been previously finished. 
	    // If so, the skip searching and writing by setting got_work to false for this batch.
	    if ( !query_batch->Empty() && (batch_count > n_batches_precomputed) ) 
	      {got_work = true;}
#else
	    if ( !query_batch->Empty() )
	      {got_work = true;}
#endif
	    //if(tid==0)std::cout<<"checkpt 1.4 \n";std::cout.flush();

	    #pragma omp barrier


	    CRef<CSearchResultSet> results;
            CLocalBlast* lcl_blast = NULL;
	    //if(tid==0)std::cout<<"checkpt 1.5 \n";std::cout.flush();

	    // No support for remote searches.
	    if ( got_work ){
	      //if(tid==0)std::cout<<"checkpt 1.6 \n";std::cout.flush();
	      
	       CRef<IQueryFactory> queries(new CObjMgr_QueryFactory(*query_batch));
	       SaveSearchStrategy(my_args, my_CmdLineArgs, queries, opts_hndl);

               //Create a local blast searcher for this thread
	       lcl_blast = new CLocalBlast(queries, opts_hndl, db_adapter);
	       lcl_blast->SetNumberOfThreads(my_CmdLineArgs->GetNumThreads());

	       //Perform search
	       //if(tid==0)std::cout<<"checkpt 1.7 searching\n";std::cout.flush();
	       results = lcl_blast->Run();
	       //if(tid==0)std::cout<<"checkpt 1.8 collecting results\n";std::cout.flush();

               //Return objects saved from local_blast object
               //Collect pieces from local_blast object used to create output
	       //if(tid==0)std::cout<<"checkpt 1.8 \n";std::cout.flush();

	       //Make duplicates of these objects so we can delete the local_blast object
               //tHSPResults[tid] =    lcl_blast->GetHSPResults();
	       tHSPResults[tid] = gbhrpfncDuplicateBlastHSPResults(lcl_blast->GetHSPResults());
               //tQueryInfo[tid] =     lcl_blast->GetQueryInfo();
	       tQueryInfo[tid] = BlastQueryInfoDup(lcl_blast->GetQueryInfo());
               //tBlastScoreBlk[tid] = lcl_blast->GetScoreBlk();
	       //this is currently a shallow copy
	       tBlastScoreBlk[tid] = gbsbpfncDuplicateBlastScoreBlk(lcl_blast->GetScoreBlk());

#if CLEAN_UP_OBJECTS
	       //if(tid==0)std::cout<<"free checkpt 0.1\n";    std::cout.flush();
	       delete(lcl_blast);
	       //if(tid==0)std::cout<<"free checkpt 0.2\n";    std::cout.flush();
#endif
	    }
	    else{//not got_work

	      //if(tid==0)std::cout<<"checkpt 2 rank="<<g_rank<<" got_work= "<<got_work<<"\n";std::cout.flush();

	       //Fill in non working thread contributions with some empty (but not NULL) structs
	      tHSPResults[tid] =    Blast_HSPResultsNew(0);
	      //results of BlastQueryInfoNew is not handled appropriately in combine function
	      //tQueryInfo[tid] =     BlastQueryInfoNew(eBlastTypeBlastp,0);//{0,0,0,NULL,0,NULL};
	      tQueryInfo[tid] =     &BlastQueryInfo{0,0,0,NULL,0,NULL};
	      tBlastScoreBlk[tid] = BlastScoreBlkNew(BLASTAA_SEQ_CODE,0);//type could be anything here
              query_batch->clear();
	    }

            #pragma omp barrier //ceb needed to merge results within thread group

	    // Merge the results in parallel within the thread group.
	    if ( local_tid == 0 ){       

	      //if(local_tid==0)std::cout<<"checkpt 3 thread_group_id="<<thread_group_id<<"\n";std::cout.flush();

	      //ceb All threads for this thread group need to be finished by this
	      // point so that we can combine results.
	      // so we need a barrier before entering this scope

	      //Concatenate and sort HSP results (by e-value) from threads in a thread group
	      HSPResultsGroup[thread_group_id] = gbhrpfncCombineTeamLeaderBlastHSPResults( tHSPResults,tid,num_team_leaders);
	      //if(local_tid==0)std::cout<<"checkpt 3.1 \n";std::cout.flush();

	      //These values should be the same inside a thread group so we only need to keep one thread
#if USE_SHALLOW_COPY
	      QueryInfoGroup[thread_group_id] = tQueryInfo[tid];
#else
	      QueryInfoGroup[thread_group_id] = BlastQueryInfoDup(tQueryInfo[tid]);
#endif
	      //if(local_tid==0)std::cout<<"checkpt 3.2 \n";std::cout.flush();
	      
	      //ceb ScoreBlk should be the same size but each corresponds to different database segment
	      //  so it's not certain how to merge these results for a thread group or what the ScoreBlk contains
	      BlastScoreBlkGroup[thread_group_id] = tBlastScoreBlk[tid];
	      //Create a duplicate of ScoreBlkGroup
	      //BlastScoreBlkGroup[thread_group_id] = gbsbpfncDuplicateBlastScoreBlk(tBlastScoreBlk[tid]);
	      //if(local_tid==0)std::cout<<"checkpt 3.3 \n";std::cout.flush();

	      
	      //ceb a vector for consolidating query batches for output
	      //we only need one entry for each thread group
	      QueryBatchGroup[thread_group_id] = query_batch;
	      //if(local_tid==0)std::cout<<"checkpt 3.4 \n";std::cout.flush();

//for(int i=0;i<query_batch->Size();++i)
//{
//std::cout<<"Gid: "<< thread_group_id <<"lcl_query_batch["<<i<<"].GetQueryId() ="<<(((query_batch.GetObject())[i].GetObject()).GetQueryId().GetObject()).AsFastaString()<<"\n"<<std::cout.flush();
//}

	    }

//if(local_tid==0)std::cout<<"checkpt 3.5 \n";std::cout.flush();
	    // Wait for all threads to finish. Then, reset all threads to 'see' the entire DB for reporting output.
            #pragma omp barrier
	       my_seq_db->SetIterationRange(0,db_num_seqs);
            #pragma omp barrier
//if(local_tid==0)std::cout<<"checkpt 3.6 \n";std::cout.flush();

#if CLEAN_UP_OBJECTS
	    //if(tid==0)std::cout<<"free checkpt 1\n";    std::cout.flush();
	    //free local thread objects
	    Blast_HSPResultsFree(tHSPResults[tid]);
	    //if(tid==0)std::cout<<"free checkpt 1.1\n";    std::cout.flush();

	    if(tQueryInfo[tid]->num_queries>0){BlastQueryInfoFree(tQueryInfo[tid]);}//this function does not handle empty objects
	    //if(tid==0)std::cout<<"free checkpt 1.2\n";    std::cout.flush();

	    //crashes here!!
	    //BlastScoreBlkFree(tBlastScoreBlk[tid]);
	    //delete doesn't crash code but free does??
	    //delete(tBlastScoreBlk[tid]);
	    //if(tid==0)std::cout<<"free checkpt 2\n";    std::cout.flush();
#endif

	    //we need to collect queries and store them in a big query_batch vector
            //This section should be called every time since at least one thread should have work here
	    if ( tid == 0 ){ //Merge our consolidated thread group data across thread groups and output to file 
	      //thread info will be merged across thread groups prior to creating the results object           
	      BlastHSPResults* HSPResultsCollection = NULL; 
	      BlastQueryInfo* QueryInfoCollection = NULL;
	      BlastScoreBlk* ScoreBlkCollection = NULL; 
	      CRef<CBlastQueryVector> QueryBatchCollection(new CBlastQueryVector);  
             
	      //std::cout<<"checkpt 4 \n";std::cout.flush();
	      HSPResultsCollection = gbhrpfncCombineThreadGroupBlastHSPResults(num_thread_groups,HSPResultsGroup);
	      //std::cout<<"checkpt 5 \n";std::cout.flush();
	      QueryInfoCollection = gbqipfncCombineThreadGroupBlastQueryInfo(num_thread_groups,QueryInfoGroup);
	      //std::cout<<"checkpt 6 \n";std::cout.flush();
	      ScoreBlkCollection = gbsbpfncCombineThreadGroupBlastScoreBlk(num_thread_groups,BlastScoreBlkGroup);
	      //std::cout<<"checkpt 7 \n";std::cout.flush();
	      QueryBatchCollection = gbqvfncCombineReplicationGroupQueries(num_thread_groups,
									   buffer,
									   gbl_bsize,
									   lintpBatchStartPosition,
									   lintpBatchEndPosition,
									   &iconfig,
									   my_CmdLineArgs->GetQueryBatchSize(),
									   gbl_query_count,
									   gbl_scope);
	       
	      //std::cout<<"checkpt 8 \n";std::cout.flush();
#if USE_RESTART
	       query_count += HSPResultsCollection->num_queries;//for restart
#endif
	       gbl_query_count += QueryBatchCollection->Size();
	       
#if CLEAN_UP_OBJECTS
	       //std::cout<<"free checkpt 3\n";    std::cout.flush();

	       //free objects merged within thread groups
	       for(int gid=0;gid<num_thread_groups;++gid)
	       {
		 Blast_HSPResultsFree(HSPResultsGroup[gid]);
		 if(QueryInfoGroup[gid]->num_queries>0){BlastQueryInfoFree(QueryInfoGroup[gid]);}
		 //BlastScoreBlkFree(BlastScoreBlkGroup[gid]);
		 //ok
		 //delete(BlastScoreBlkGroup[gid]);
	       }
	       //std::cout<<"free checkpt 4\n";    std::cout.flush();
#endif
	       if(HSPResultsCollection->num_queries>0)
	       {
		 //std::cout<<"checkpt 9 num_queries="<<HSPResultsCollection->num_queries<<"\n";std::cout.flush();
		 
		 CRef<IQueryFactory> queries(new CObjMgr_QueryFactory(*QueryBatchCollection));
		 //std::cout<<"checkpt 10 \n";std::cout.flush();
		 
		 //This SeqInfoSrc needs to be created and filled with appropriate info
		 CRef<IBlastSeqInfoSrc> SeqInfoSrc;
		 SeqInfoSrc.Reset(db_adapter->MakeSeqInfoSrc());
		 //std::cout<<"checkpt 11 \n";std::cout.flush();
		 
		 results = gsrsfncCollectResults(HSPResultsCollection, 
						 QueryInfoCollection,
						 ScoreBlkCollection,
						 queries, 
						 &opt,
						 SeqInfoSrc);
		 
		 //std::cout<<"checkpt 12 \n";std::cout.flush();
		 BlastFormatter_PreFetchSequenceData(*results, gbl_scope);
		 
		 //std::cout<<"checkpt 13 \n";std::cout.flush();
		 
		 //Loop over results
		 ITERATE(CSearchResultSet, result, *results)
		   {
		     formatter->PrintOneResultSet(**result, QueryBatchCollection);
		   }
		 //std::cout<<"checkpt 14 \n";std::cout.flush();
		 
		 //Need to make sure we don't accidentally update the number of completed batches before they actually finish writing
		 my_new_output->flush();
		 //std::cout<<"checkpt 15 \n";std::cout.flush();
	       }

#if CLEAN_UP_OBJECTS
	       //if(tid==0)std::cout<<"free checkpt 5\n";    std::cout.flush();
	       //free objects merged across thread groups
	       //crashes code
	       //Blast_HSPResultsFree(HSPResultsCollection);
	       //delete(HSPResultsCollection);
	       //ok
	       BlastQueryInfoFree(QueryInfoCollection);
	       //crashes code
	       //BlastScoreBlkFree(ScoreBlkCollection);
	       //if(tid==0)std::cout<<"free checkpt 6\n";    std::cout.flush();
#endif


#if USE_RESTART
	       // Record successful output in restart log file
	       {
		 restart_stream->open(my_restart.c_str(),std::fstream::trunc | std::fstream::out);
		 //printf("writing to restart_file= %s\n",my_restart.c_str());
		 //write total number of batches and queries processed so far
		 *restart_stream << batch_count << " " << query_count << '\n';
		 restart_stream->close();
	       }
#endif

	    }//end if tid==0

            #pragma omp single
	    {working = false;}

	    #pragma omp critical
	    {
	       if ( !lcl_input.End() ){ working = true;}
	    }
	    #pragma omp barrier

            #pragma omp flush(working)

            #pragma omp barrier
	    //if(tid==0){std::cout<<"checkpt 16\n";std::cout.flush();}

         }//end working loop over superqueries
	    
	  //if(tid==0){std::cout<<"checkpt 17 \n";std::cout.flush();}


	 if(tid==0) {formatter->PrintEpilog(opt);}

	 //if(tid==0){std::cout<<"checkpt 18\n";std::cout.flush();}

#if USE_RESTART
	 // Complete output successful.
	 // We can ensure next execution creates new file by setting restart counter
	 // values to zero.
	 if(tid==0)
	 {
	   restart_stream->open(my_restart.c_str(),std::fstream::trunc | std::fstream::out);
	   //printf("reset restart_file= %s\n",my_restart.c_str());
	   //write total number of batches and queries processed so far
	   *restart_stream << 0 << " " << 0 << '\n';
	   restart_stream->close();  
	 }
#endif
	 //if(tid==0){std::cout<<"checkpt 18\n";std::cout.flush();}
         if (my_CmdLineArgs->ProduceDebugOutput()) {
            opts_hndl->GetOptions().DebugDumpText(NcbiCerr, "BLAST options", 1);
	    //if(tid==0){std::cout<<"checkpt 19\n";std::cout.flush();}
         }
      } CATCH_ALL(status)
	  //if(tid==0){std::cout<<"checkpt 20\n";std::cout.flush();}
   }   // end the parallel block.
    //std::cout<<"checkpt 21\n";std::cout.flush();

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

	      file_name += ".pin";
	      if ( ! file_exists(file_name) )
		{
		  std::cerr<<"FATAL ERROR: The database file does not exist: "<<file_name<<std::endl;
		  std::cerr.flush();
		  state = -1;
		}

	      file_name.assign(db_str);
	      file_name += ".psq";
	      if ( ! file_exists(file_name) )
		{
		  std::cerr<<"FATAL ERROR: The database file does not exist: "<<file_name<<std::endl;
		  std::cerr.flush();
		  state = -1;
		}

	      file_name.assign(db_str);
	      file_name += ".phr";
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

		  file_name.assign( sstm.str() + ".pin" );
		  if ( ! file_exists(file_name) )
		    {
		      std::cerr<<"FATAL ERROR: The database file does not exist: "<<file_name<<std::endl;
		      std::cerr.flush();
		      state = -1;
		    }

		  file_name.assign( sstm.str() + ".psq" );
		  if ( ! file_exists(file_name) )
		    {
		      std::cerr<<"FATAL ERROR: The database file does not exist: "<<file_name<<std::endl;
		      std::cerr.flush();
		      state = -1;
		    }

		  file_name.assign( sstm.str() + ".phr" );
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

  int blast_return = CBlastpApp().AppMain(argc, argv, 0, eDS_Default, 0);

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


// Functions to combine objects across within thread groups ==================================

BlastHSPResults* gbhrpfncCombineTeamLeaderBlastHSPResults ( char ** HSPResultsBuffers,
							    size_t * buffer_size,
							    int tid,
							    int num_team_leaders )
{
  Int4 num_queries = 0;
  Int4 dummy_int4 = 0;  // dummy variable.
  double dummy_double = 0.;
  Boolean dummy_boolean;
  Boolean all_heapified = 0;

  Int4 total_hsplist_current = 0;
  Int4 total_hsplist_count = 0;
  Int4 total_hsplist_max=0;
  Int4 total_hsparray_allocated = 0;

  double worst_evalue_all=0.;
  Int4 low_score_all=0;

  // Need a buffer pointer for each different buffer.
  char **buff_ptr = NULL;

  Int4 *hsplist_counts = NULL;
  //ceb
  Int4 *total_hsplist_counts = NULL;

  buff_ptr = (char**)malloc(num_team_leaders * sizeof(char*));
  hsplist_counts = (Int4*)malloc(num_team_leaders * sizeof(char*));

  //std::cout<<std::endl;//ceb

  BlastHSPResults *hsp_results = NULL;    // The guy we're returning.

  // Before we unpack the buffers into the object we have to know how many of certain objects their are.
  // So, first we crawl across the buffers to count the stuff up.
  for ( int i=0; i < num_team_leaders; ++i )
    {
      //std::cout<< "chkpt -1 tid: " << tid <<" i: "<< i <<std::endl;std::cout.flush();
      buff_ptr[i] = HSPResultsBuffers[tid+i];
    }


  // Get the number of queries -- its the same across all threads in a thread group.
  memcpy( &num_queries , buff_ptr[0], sizeof(Int4) );
  //ceb Seems to fix the offset mismatch error when going in to query 0+1
  for ( int i=0; i < num_team_leaders; ++i )
    { buff_ptr[i] += sizeof(Int4); }//each thread will have num_queries, we only read it in once outside this loop and then skip here.


  //ceb
  //total_hsplist_counts = (Int4*)malloc(num_queries * sizeof(char*));
  //for(int i=0;i<num_queries;++i){total_hsplist_counts[i]=0;}


  //if(tid>=0){std::cout<<"tid: "<<tid<<" Found "<<num_queries<<" num_queries in the first buffer."<<std::endl;}

  hsp_results = Blast_HSPResultsNew(num_queries);

  // We have to move across the buffer in a very particular way (reverse of how its packed.)
  // See traceback_stage.cpp for details on how its packed.

  for ( int BlastHSPResults_query_index = 0; BlastHSPResults_query_index < hsp_results->num_queries ; ++BlastHSPResults_query_index )
    {    
      // Allocate the BlastHitList for this query.
      BlastHitList *bhl = Blast_HitListNew(1);

      //if(tid>=0){std::cout<<"tid: "<<tid<<" Query "<<BlastHSPResults_query_index<<std::endl;std::cout.flush();}

      // hsplist_count
      total_hsplist_count=0;//ceb Reset here for each query
      for ( int i=0; i < num_team_leaders; ++i )
	{
	  // read past num_queries
	  //buff_ptr[i] += sizeof(Int4);//each thread will have num_queries, we only read it in once outside this loop and then skip here.
          memcpy( &dummy_int4, buff_ptr[i], sizeof(Int4) );
	  buff_ptr[i] += sizeof(Int4);
	  
	  //if(tid>=0){std::cout<<"  "<<"Team leader "<<i<<" has hsplist_count "<<dummy_int4<<std::endl; std::cout.flush();}

          //ceb We want to merge hsplist data across multiple team leaders in a thread group for a given query.
	  total_hsplist_count += dummy_int4;

	  //total_hsplist_count[BlastHSPResults_query_index] += dummy_int4;//ceb

	  hsplist_counts[i] = dummy_int4;
	}

      //if(tid>=0){std::cout<<"tid: "<<tid<<" Query "<<BlastHSPResults_query_index<<" Total hsplist_count "<<total_hsplist_count<<std::endl; std::cout.flush();}

      bhl->hsplist_count = total_hsplist_count;

      //if(tid==0){std::cout<<"checkpt 1 "<<std::endl; std::cout.flush();}

      // hsplist_max
      for ( int i=0; i < num_team_leaders; ++i )
	{
	  memcpy( &dummy_int4, buff_ptr[i], sizeof(Int4) );
	  buff_ptr[i] += sizeof(Int4);

	  //if(tid>=0){std::cout<<"  "<<"Team leader "<<i<<" has hsplist_max "<<dummy_int4<<std::endl; std::cout.flush();}

	  if ( i == 0 )
	    {
	      //memcpy( &dummy_int4, buff_ptr[i], sizeof(Int4) );
	      total_hsplist_max = dummy_int4;
	    }
	}

      //if(tid>=0){std::cout<<"tid: "<<tid<<" pre mult total_hsplist_max: "<<total_hsplist_max<< std::endl; std::cout.flush();}

      while ( total_hsplist_max < total_hsplist_count )
	{total_hsplist_max *= 2;}

      //if(tid>=0){std::cout<<"tid: "<<tid<<" post mult Total hsplist max is "<<total_hsplist_max<<std::endl;std::cout.flush();}

      bhl->hsplist_max = total_hsplist_max;

      // worst evalue
      memcpy( &dummy_double, buff_ptr[0], sizeof(double));
      buff_ptr[0] += sizeof(double);
      worst_evalue_all = dummy_double;
      //if(tid>=0){std::cout<<"tid: "<<tid<<"checkpt 4 "<<std::endl; std::cout.flush();}
      //if(tid>=0){std::cout<<"tid: "<<tid<<"  "<<"Team leader "<<"0"<<" has worst evalue "<<dummy_double<<std::endl;std::cout.flush();}
      for ( int i=1; i < num_team_leaders; ++i )
	{
	  memcpy( &dummy_double, buff_ptr[i], sizeof(double) );
	  buff_ptr[i] += sizeof(double);
	  //if(tid>=0){std::cout<<"tid: "<<tid<<"  "<<"Team leader "<<i<<" has worst evalue "<<dummy_double<<std::endl;std::cout.flush();}
	  if( dummy_double > worst_evalue_all )
	    {worst_evalue_all = dummy_double;}
	}
      //if(tid>=0){std::cout<<"tid: "<<tid<<"  "<<"Worst evalue all is "<<worst_evalue_all<<std::endl;std::cout.flush();}
      bhl->worst_evalue = worst_evalue_all;


      // best bit(low) score
      memcpy( &dummy_int4, buff_ptr[0], sizeof(Int4));
      buff_ptr[0] += sizeof(Int4);
      low_score_all = dummy_int4;
      //if(tid>=0){std::cout<<"tid: "<<tid<<"checkpt 6 "<<std::endl; std::cout.flush();}
      //if(tid>=0){std::cout<<"  "<<"Team leader "<<"0"<<" has low score "<<dummy_int4<<std::endl;std::cout.flush();}
      for ( int i=1; i < num_team_leaders; ++i )
	{
	  memcpy( &dummy_int4, buff_ptr[i], sizeof(Int4));
	  buff_ptr[i] += sizeof(Int4);
	  //if(tid>=0){std::cout<<"tid: "<<tid<<"  "<<"Team leader "<<i<<" has low score "<<dummy_int4<<std::endl;std::cout.flush();}
	  if( dummy_int4 < low_score_all )
	    {low_score_all = dummy_int4;}
	}
      //if(tid>=0){std::cout<<"tid: "<<tid<<"  "<<"Low score all is "<<low_score_all<<std::endl;std::cout.flush();}
      bhl->low_score = low_score_all;

      //heapified boolean
      for ( int i=0; i < num_team_leaders; ++i )
	{
	  memcpy( &dummy_boolean, buff_ptr[i], sizeof(Boolean) );
	  buff_ptr[i] += sizeof(Boolean);
	  //if(tid>=0)std::cout<<"tid: "<<tid<<"  "<<"Team leader "<<i<<" has heapified "<<(int)dummy_boolean<<std::endl;std::cout.flush();
	  if ( dummy_boolean == 1 )
	    {all_heapified = dummy_boolean;}
	}
      // if any are heapified assume the merge can be.
      //if(tid>=0){std::cout<<"tid: "<<tid<<"  "<<"all heapified is "<<(int)all_heapified<<std::endl;std::cout.flush();}
      bhl->heapified = all_heapified;

      //hsplist_current
      for ( int i=0; i < num_team_leaders; ++i )
	{
	  memcpy( &dummy_int4, buff_ptr[i], sizeof(Int4) );
	  buff_ptr[i] += sizeof(Int4);
	  //if(tid>=0){std::cout<<"tid: "<<tid<<"  "<<"Team leader "<<i<<" has hsplist_current "<<dummy_int4<<std::endl;std::cout.flush();}
	}
      //if(tid>=0){std::cout<<"tid: "<<tid<<"checkpt 10 "<<std::endl; std::cout.flush();}


      // take the first one or the max. add 100 until it surpasses hsplist_count!
      while ( dummy_int4 < total_hsplist_count )
	{dummy_int4 += 100;}

      total_hsplist_current = dummy_int4;
      //if(tid>=0){std::cout<<"tid: "<<tid<<"  Total hsplist_current is "<<total_hsplist_current<<std::endl; std::cout.flush();}
      bhl->hsplist_current = total_hsplist_current;

      // Allocate the hspList to current long.
      bhl->hsplist_array = (BlastHSPList**)calloc( bhl->hsplist_current, sizeof(BlastHSPList*) );

      //if(tid>=0){std::cout<<"tid: "<<tid<<"checkpt 11 "<<std::endl; std::cout.flush();}

      Int4 global_hsplist_index = 0;
      // Loop over each team leader and process the hit lists.
      for ( int i=0; i < num_team_leaders; ++i )
	{
	  //if(tid>=0){std::cout<<"tid: "<<tid<<"checkpt 12 hsplist_counts["<<i<<"]: "<<hsplist_counts[i]<<std::endl; std::cout.flush();}

	  //for ( Int4 BlastHitList_hsplist_index; BlastHitList_hsplist_index < hsplist_counts[i]; ++BlastHitList_hsplist_index )
	  for ( int BlastHitList_hsplist_index=0; BlastHitList_hsplist_index < hsplist_counts[i]; ++BlastHitList_hsplist_index )
	    {
	      // Allocate a BlastHSPList ptr
	      BlastHSPList* bhspl = (BlastHSPList*) calloc(1,sizeof(BlastHSPList) );

	      memcpy( &(bhspl->oid) , buff_ptr[i], sizeof(Int4) );
	      buff_ptr[i] += sizeof(Int4);
	      //if(tid>=0){std::cout<<"tid: "<<tid<<" bhspl->oid: "<<bhspl->oid<<std::endl; std::cout.flush();}
	      
	      memcpy( &(bhspl->query_index) , buff_ptr[i], sizeof(Int4) );
	      buff_ptr[i] += sizeof(Int4);

	      memcpy( &(bhspl->hspcnt) , buff_ptr[i], sizeof(Int4) );
	      buff_ptr[i] += sizeof(Int4);

	      memcpy( &(bhspl->allocated) , buff_ptr[i], sizeof(Int4) );
	      buff_ptr[i] += sizeof(Int4);

	      memcpy( &(bhspl->hsp_max) , buff_ptr[i], sizeof(Int4) );
	      buff_ptr[i] += sizeof(Int4);
	      //if(tid>=0){std::cout<<"tid: "<<tid<<" bhspl->hsp_max: "<<bhspl->hsp_max<<std::endl; std::cout.flush();}

	      memcpy( &(bhspl->do_not_reallocate) , buff_ptr[i], sizeof(Boolean) );
	      buff_ptr[i] += sizeof(Boolean);

	      memcpy( &(bhspl->best_evalue) , buff_ptr[i], sizeof(double) );
	      buff_ptr[i] += sizeof(double);
	      //if(tid>=0){std::cout<<"tid: "<<tid<<" bhspl->best_evalue "<<bhspl->best_evalue<<std::endl; std::cout.flush();}

	      // Allocate the array for High Scoring segment pairs.
	      bhspl->hsp_array = (BlastHSP**) calloc( bhspl->allocated, sizeof(BlastHSP*) );
	      //if(tid>=0){std::cout<<"tid: "<<tid<<"bhspl->hspcnt "<<bhspl->hspcnt<<std::endl; std::cout.flush();}

	      // Loop over all HSPs between the query and the subject sequence.
	      for ( int BlastHSP_hsp_array_index=0; BlastHSP_hsp_array_index < bhspl->hspcnt; ++BlastHSP_hsp_array_index )
		{
		  BlastHSP* bhsp = (BlastHSP*) calloc (1, sizeof(BlastHSP) );

		  memcpy ( &(bhsp->score), buff_ptr[i] , sizeof(Int4) );
		  buff_ptr[i] += sizeof(Int4);
		  
		  memcpy ( &(bhsp->num_ident), buff_ptr[i] , sizeof(Int4) );
		  buff_ptr[i] += sizeof(Int4);
		
		  memcpy ( &(bhsp->bit_score), buff_ptr[i] , sizeof(double) );
		  buff_ptr[i] += sizeof(double);

		  memcpy ( &(bhsp->evalue), buff_ptr[i] , sizeof(double) );
		  buff_ptr[i] += sizeof(double);

		  memcpy ( &(bhsp->context), buff_ptr[i] , sizeof(Int4) );
		  buff_ptr[i] += sizeof(Int4);

		  memcpy ( &(bhsp->num), buff_ptr[i] , sizeof(Int4) );
		  buff_ptr[i] += sizeof(Int4);

		  memcpy ( &(bhsp->comp_adjustment_method), buff_ptr[i] , sizeof(Int2) );
		  buff_ptr[i] += sizeof(Int2);

		  memcpy ( &(bhsp->num_positives), buff_ptr[i] , sizeof(Int4) );
		  buff_ptr[i] += sizeof(Int4);

		  memcpy ( &(bhsp->query.frame), buff_ptr[i] , sizeof(Int2) );
		  buff_ptr[i] += sizeof(Int2);

		  memcpy ( &(bhsp->query.offset), buff_ptr[i] , sizeof(Int4) );
		  buff_ptr[i] += sizeof(Int4);

		  memcpy ( &(bhsp->query.end), buff_ptr[i] , sizeof(Int4) );
		  buff_ptr[i] += sizeof(Int4);

		  memcpy ( &(bhsp->query.gapped_start), buff_ptr[i] , sizeof(Int4) );
		  buff_ptr[i] += sizeof(Int4);

		  memcpy ( &(bhsp->subject.frame), buff_ptr[i] , sizeof(Int2) );
		  buff_ptr[i] += sizeof(Int2);

		  memcpy ( &(bhsp->subject.offset), buff_ptr[i] , sizeof(Int4) );
		  buff_ptr[i] += sizeof(Int4);

		  memcpy ( &(bhsp->subject.end), buff_ptr[i] , sizeof(Int4) );
		  buff_ptr[i] += sizeof(Int4);

		  memcpy ( &(bhsp->subject.gapped_start), buff_ptr[i] , sizeof(Int4) );
		  buff_ptr[i] += sizeof(Int4);
		  //if(tid>=0){std::cout<<"bhsp->subject.gapped_start "<<bhsp->subject.gapped_start<<std::endl; std::cout.flush();}

		  // The next value might be a setinel value.
		  Int4 sentinel;
                  
 
		  memcpy ( &sentinel , buff_ptr[i], sizeof(Int4) );

		  if ( sentinel )
		    {
		      SPHIHspInfo *my_pat_info = (SPHIHspInfo*) calloc( 1, sizeof(SPHIHspInfo) );

		      memcpy( &(my_pat_info->index) , buff_ptr[i], sizeof(Int4) );
		      buff_ptr[i] += sizeof(Int4);

		      memcpy( &(my_pat_info->length) , buff_ptr[i], sizeof(Int4) );
		      //if(tid>=0){std::cout<<"tid: "<<tid<<"my_pat_info->length "<<my_pat_info->length<<std::endl; std::cout.flush();}
		      buff_ptr[i] += sizeof(Int4);

		      bhsp->pat_info = my_pat_info ;
		    }
		  else
		    {
		      bhsp->pat_info = NULL;
		      buff_ptr[i] += sizeof(Int4);
		    }

		  // The next value might be a sentinel value on GapEditScript.
		  memcpy ( &sentinel , buff_ptr[i] , sizeof(Int4) );

		  if ( sentinel )
		    {
		      GapEditScript *my_gap_info = (GapEditScript*) calloc(1, sizeof(GapEditScript) );

		      memcpy( &(my_gap_info->size) , buff_ptr[i] , sizeof(Int4) );
		      //if(tid>=0){std::cout<<"tid: "<<tid<<" my_gap_info->size "<<my_gap_info->size<<std::endl; std::cout.flush();}

		      buff_ptr[i] += sizeof(Int4);

		      my_gap_info->op_type = (EGapAlignOpType*) calloc( my_gap_info->size , sizeof(EGapAlignOpType) );
		      memcpy( &(my_gap_info->op_type[0]), buff_ptr[i] , my_gap_info->size * sizeof(EGapAlignOpType) );
		      buff_ptr[i] += ( my_gap_info->size * sizeof(EGapAlignOpType) );

		      my_gap_info->num = (Int4*) calloc( my_gap_info->size , sizeof(Int4) );
		      memcpy( &(my_gap_info->num[0]), buff_ptr[i] , my_gap_info->size * sizeof(Int4) );
		      buff_ptr[i] += (my_gap_info->size * sizeof(Int4) );

		      bhsp->gap_info = my_gap_info ;
		    }
		  else
		    {
		      bhsp->gap_info = NULL; 
		      buff_ptr[i] += sizeof(Int4);
		    }

		  // Now attach this BlastHSP to the parent BlastHSPList.
		  bhspl->hsp_array[ BlastHSP_hsp_array_index ] = bhsp ;

		} // End of loop over the BlastHSP_hsp_array_index -> the high scoring pair data

              //if(tid==0){std::cout<<"tid: "<<tid<<" checkpt 19 "<<std::endl; std::cout.flush();}

	      // Now attach this BlastHSPList to the parent BlastHitList
	      bhl->hsplist_array[ global_hsplist_index ] = bhspl ;

	      ++global_hsplist_index;  // increment the global index
	      
	    }  // End loop over BlastHitList_hsplist_index

	} // End Loop over team leaders.

      //if(tid==0){std::cout<<"tid: "<<tid<<" checkpt 20 "<<std::endl; std::cout.flush();}

      // Now we can add this BlastHitList to the parent BlastHSPResults
      hsp_results->hitlist_array[ BlastHSPResults_query_index ] = bhl ;

  }

  //if(tid>=0){std::cout<<"tid: "<<tid<<" final hsp_results->num_queries: "<<hsp_results->num_queries <<std::endl; std::cout.flush();}
  // for ( int BlastHSPResults_query_index = 0; BlastHSPResults_query_index < hsp_results->num_queries ; ++BlastHSPResults_query_index )
  //  {
  //    if(tid>=0){std::cout<<"tid: "<<tid<<" Query "<<BlastHSPResults_query_index<<std::endl; std::cout.flush();}
  //    //if(tid>=0){std::cout<<"tid: "<<tid<<" final hsp_results->hitlist_array["<< BlastHSPResults_query_index <<"]->worst_evalue: "<< hsp_results->hitlist_array[BlastHSPResults_query_index]->worst_evalue <<std::endl; std::cout.flush();}
  //    if(tid>=0){std::cout<<"tid: "<<tid<<" final hsp_results->hitlist_array["<< BlastHSPResults_query_index <<"]->hsplist_count: "<< hsp_results->hitlist_array[BlastHSPResults_query_index]->hsplist_count <<std::endl; std::cout.flush();}
  //    //if(tid>=0){std::cout<<"tid: "<<tid<<" final hsp_results->hitlist_array["<< BlastHSPResults_query_index <<"]->hsplist_max: "<< hsp_results->hitlist_array[BlastHSPResults_query_index]->hsplist_max <<std::endl; std::cout.flush();}
  //    //if(tid>=0){std::cout<<"tid: "<<tid<<" final hsp_results->hitlist_array["<< BlastHSPResults_query_index <<"]->low_score: "<< hsp_results->hitlist_array[BlastHSPResults_query_index]->low_score <<std::endl; std::cout.flush();}
  //    //if(tid>=0){std::cout<<"tid: "<<tid<<" final hsp_results->hitlist_array["<< BlastHSPResults_query_index <<"]->heapified: "<< hsp_results->hitlist_array[BlastHSPResults_query_index]->heapified <<std::endl; std::cout.flush();}
  //    if(tid>=0){std::cout<<"tid: "<<tid<<" final hsp_results->hitlist_array["<< BlastHSPResults_query_index <<"]->hsplist_current: "<< hsp_results->hitlist_array[BlastHSPResults_query_index]->hsplist_current <<std::endl; std::cout.flush();}
  //  }



  //if(tid>=0){std::cout<<"tid: "<<tid<<" Pre SortByValue "<<std::endl; std::cout.flush();}
  // Now we can sort the HSPResults.
  Blast_HSPResultsSortByEvalue( hsp_results ) ;
  //if(tid>=0){std::cout<<"tid: "<<tid<<" Post SortByValue "<<std::endl; std::cout.flush();}


  return hsp_results;

}








BlastHSPResults* gbhrpfncCombineTeamLeaderBlastHSPResults ( BlastHSPResults** HSPResultsArray,
							    int tid,
							    int num_team_leaders )
{
  Int4 num_queries = 0;
  Int4 dummy_int4 = 0;  // dummy variable.
  double dummy_double = 0.;
  Boolean dummy_boolean;
  Boolean all_heapified = 0;

  Int4 total_hsplist_current = 0;
  Int4 total_hsplist_count = 0;
  Int4 total_hsplist_max=0;
  Int4 total_hsparray_allocated = 0;

  double worst_evalue_all=0.;
  Int4 low_score_all=0;

  // Need a buffer pointer for each different buffer.
  Int4 *hsplist_counts = NULL;
  BlastHSPResults **hsp_results =  new BlastHSPResults*[num_team_leaders];    // use this to point to a value in HSPResultsArray
  BlastHSPResults *merged_hsp_results = NULL;    // The guy we're returning.

  hsplist_counts = (Int4*)malloc(num_team_leaders * sizeof(char*));

  for ( int i=0; i < num_team_leaders; ++i )
  {
    hsp_results[i] = HSPResultsArray[tid+i];
  }

  // Get the number of queries -- its the same across all threads in a thread group.
  num_queries = (hsp_results[0])->num_queries;
  //std::cout<<"(hsp_results[0])->num_queries= "<<num_queries <<std::endl;std::cout.flush();
  merged_hsp_results = Blast_HSPResultsNew(num_queries);

  // We have to move across the buffer in a very particular way (reverse of how its packed.)
  // See traceback_stage.cpp for details on how its packed.

  for ( int BlastHSPResults_query_index = 0; BlastHSPResults_query_index < merged_hsp_results->num_queries ; ++BlastHSPResults_query_index )
    {    
      // Allocate the BlastHitList for this query.
      BlastHitList *merged_bhl = Blast_HitListNew(1);

      // hsplist_count
      total_hsplist_count=0;//ceb Reset here for each query
      for ( int i=0; i < num_team_leaders; ++i )
	{
          BlastHitList *bhl = (hsp_results[i])->hitlist_array[BlastHSPResults_query_index];
          //dummy_int4 = bhl->hsplist_count;
	  //total_hsplist_count += dummy_int4;
	  //hsplist_counts[i] = dummy_int4;
	  total_hsplist_count += bhl->hsplist_count;
	  hsplist_counts[i] = bhl->hsplist_count;
	}

      merged_bhl->hsplist_count = total_hsplist_count;

      // hsplist_max
      total_hsplist_max=0;
      for ( int i=0; i < num_team_leaders; ++i )
	{
          BlastHitList *bhl = (hsp_results[i])->hitlist_array[BlastHSPResults_query_index];
	  if ( i == 0 )
	    {
	      total_hsplist_max = bhl->hsplist_max;
	    }
	}

      while ( total_hsplist_max < total_hsplist_count )
	{total_hsplist_max *= 2;}

      merged_bhl->hsplist_max = total_hsplist_max;

      // worst evalue
      worst_evalue_all = ((hsp_results[0])->hitlist_array[BlastHSPResults_query_index])->worst_evalue;
      for ( int i=1; i < num_team_leaders; ++i )
	{
          BlastHitList *bhl = (hsp_results[i])->hitlist_array[BlastHSPResults_query_index];
          dummy_double = bhl->worst_evalue;
	  if( dummy_double > worst_evalue_all )
	    {worst_evalue_all = dummy_double;}
	}
      merged_bhl->worst_evalue = worst_evalue_all;


      // best bit(low) score
      low_score_all = ((hsp_results[0])->hitlist_array[BlastHSPResults_query_index])->low_score;
      for ( int i=1; i < num_team_leaders; ++i )
	{
          BlastHitList *bhl = (hsp_results[i])->hitlist_array[BlastHSPResults_query_index];
	  dummy_int4 = bhl->low_score;
	  if( dummy_int4 < low_score_all )
	    {low_score_all = dummy_int4;}
	}
      merged_bhl->low_score = low_score_all;

      //heapified boolean
      all_heapified=0;
      for ( int i=0; i < num_team_leaders; ++i )
	{
          BlastHitList *bhl = (hsp_results[i])->hitlist_array[BlastHSPResults_query_index];
	  dummy_boolean = bhl->heapified;
	  if ( dummy_boolean == 1 )
	    {all_heapified = dummy_boolean;}
	}
      // if any are heapified assume the merge can be.
      merged_bhl->heapified = all_heapified;

      //hsplist_current
      dummy_int4 = ((hsp_results[0])->hitlist_array[BlastHSPResults_query_index])->hsplist_current;

      // take the first one or the max. add 100 until it surpasses hsplist_count!
      while ( dummy_int4 < total_hsplist_count )
	{dummy_int4 += 100;}

      total_hsplist_current = dummy_int4;
      merged_bhl->hsplist_current = total_hsplist_current;

      // Allocate the hspList to current long.
      merged_bhl->hsplist_array = (BlastHSPList**)calloc( merged_bhl->hsplist_current, sizeof(BlastHSPList*) );

      Int4 global_hsplist_index = 0;
      // Loop over each team leader and process the hit lists.
      for ( int i=0; i < num_team_leaders; ++i )
      {
        BlastHitList *bhl = (hsp_results[i])->hitlist_array[BlastHSPResults_query_index];
	for ( int BlastHitList_hsplist_index=0; BlastHitList_hsplist_index < hsplist_counts[i]; ++BlastHitList_hsplist_index )
	    {
	      // Allocate a BlastHSPList ptr
	      BlastHSPList *bhspl = bhl->hsplist_array[BlastHitList_hsplist_index];

#if USE_SHALLOW_COPY
	      BlastHSPList* merged_bhspl = (BlastHSPList*) calloc(1,sizeof(BlastHSPList) );

	      merged_bhspl->oid = bhspl->oid;
	      merged_bhspl->query_index = bhspl->query_index;
	      merged_bhspl->hspcnt = bhspl->hspcnt;
	      merged_bhspl->allocated = bhspl->allocated;
	      merged_bhspl->hsp_max = bhspl->hsp_max;
	      merged_bhspl->do_not_reallocate = bhspl->do_not_reallocate;
              merged_bhspl->best_evalue = bhspl->best_evalue;

	      // Allocate the array for High Scoring segment pairs.
	      merged_bhspl->hsp_array = (BlastHSP**) calloc( merged_bhspl->allocated, sizeof(BlastHSP*) );

	      // Loop over all HSPs between the query and the subject sequence.
	      for ( int BlastHSP_hsp_array_index=0; BlastHSP_hsp_array_index < bhspl->hspcnt; ++BlastHSP_hsp_array_index )
		{
		  BlastHSP* merged_bhsp = (BlastHSP*) calloc (1, sizeof(BlastHSP) );
		  BlastHSP* bhsp = bhspl->hsp_array[BlastHSP_hsp_array_index];

		  merged_bhsp->score = bhsp->score;
		  merged_bhsp->num_ident = bhsp->num_ident;
		  merged_bhsp->bit_score = bhsp->bit_score;
		  merged_bhsp->evalue = bhsp->evalue;
		  merged_bhsp->context = bhsp->context;
		  merged_bhsp->num = bhsp->num;
		  merged_bhsp->comp_adjustment_method = bhsp->comp_adjustment_method;
		  merged_bhsp->num_positives = bhsp->num_positives;
		  merged_bhsp->query.frame = bhsp->query.frame;
		  merged_bhsp->query.offset = bhsp->query.offset;
		  merged_bhsp->query.end = bhsp->query.end;
		  merged_bhsp->query.gapped_start = bhsp->query.gapped_start;
		  merged_bhsp->subject.frame = bhsp->subject.frame;
		  merged_bhsp->subject.offset = bhsp->subject.offset;
		  merged_bhsp->subject.end = bhsp->subject.end;
		  merged_bhsp->subject.gapped_start = bhsp->subject.gapped_start;

		  SPHIHspInfo *my_pat_info = bhsp->pat_info;
		  SPHIHspInfo *merged_pat_info = NULL;

                  if(my_pat_info!=NULL)
		  {
		    merged_pat_info = (SPHIHspInfo*) calloc( 1, sizeof(SPHIHspInfo) );
		    merged_pat_info->index = my_pat_info->index;
		    merged_pat_info->length = my_pat_info->length;
		  }
		  merged_bhsp->pat_info = merged_pat_info ;
		  
		  GapEditScript* my_gap_info = bhsp->gap_info;
		  GapEditScript *merged_gap_info = NULL;

                  if(my_gap_info != NULL)
		  {
		    merged_gap_info = (GapEditScript*) calloc(1, sizeof(GapEditScript) );
		    merged_gap_info->size = my_gap_info->size;
		    merged_gap_info->op_type = (EGapAlignOpType*) calloc( merged_gap_info->size , sizeof(EGapAlignOpType) );
		    memcpy(&(merged_gap_info->op_type[0]) , &(my_gap_info->op_type[0]) , my_gap_info->size * sizeof(EGapAlignOpType) );
		    merged_gap_info->num = (Int4*) calloc( merged_gap_info->size , sizeof(Int4) );
		    memcpy( &(merged_gap_info->num[0]) , &(my_gap_info->num[0]) , my_gap_info->size * sizeof(Int4) );
		  }
		  merged_bhsp->gap_info = merged_gap_info;

		  // Now attach this BlastHSP to the parent BlastHSPList.
		  merged_bhspl->hsp_array[ BlastHSP_hsp_array_index ] = merged_bhsp ;

		} // End of loop over the BlastHSP_hsp_array_index -> the high scoring pair data


#else
	      // Now attach this BlastHSP to the parent BlastHSPList.
	      BlastHSPList* merged_bhspl = BlastHSPListDuplicate(bhspl);
#endif


	      // Now attach this BlastHSPList to the parent BlastHitList
	      merged_bhl->hsplist_array[ global_hsplist_index ] = merged_bhspl ;
	      ++global_hsplist_index;  // increment the global index
	      
	    }  // End loop over BlastHitList_hsplist_index

	} // End Loop over team leaders.


      // Now we can add this BlastHitList to the parent BlastHSPResults
      merged_hsp_results->hitlist_array[ BlastHSPResults_query_index ] = merged_bhl ;
  }

  // Now we can sort the HSPResults.
  Blast_HSPResultsSortByEvalue( merged_hsp_results ) ;

  return merged_hsp_results;
}













BlastHSPResults* gbhrpfncDuplicateBlastHSPResults ( BlastHSPResults* HSPResultsArray)
{
  int tid = 0;
  int num_team_leaders=1; 
  Int4 num_queries = 0;
  Int4 dummy_int4 = 0;  // dummy variable.
  double dummy_double = 0.;
  Boolean dummy_boolean;
  Boolean all_heapified = 0;

  Int4 total_hsplist_current = 0;
  Int4 total_hsplist_count = 0;
  Int4 total_hsplist_max=0;
  Int4 total_hsparray_allocated = 0;

  double worst_evalue_all=0.;
  Int4 low_score_all=0;

  // Need a buffer pointer for each different buffer.
  Int4 *hsplist_counts = NULL;
  BlastHSPResults **hsp_results =  new BlastHSPResults*[num_team_leaders];    // use this to point to a value in HSPResultsArray
  BlastHSPResults *merged_hsp_results = NULL;    // The guy we're returning.

  hsplist_counts = (Int4*)malloc(num_team_leaders * sizeof(char*));

  for ( int i=0; i < num_team_leaders; ++i )
  {
    hsp_results[i] = HSPResultsArray;
  }

  // Get the number of queries -- its the same across all threads in a thread group.
  num_queries = (hsp_results[0])->num_queries;
  //std::cout<<"(hsp_results[0])->num_queries= "<<num_queries <<std::endl;std::cout.flush();
  merged_hsp_results = Blast_HSPResultsNew(num_queries);

  // We have to move across the buffer in a very particular way (reverse of how its packed.)
  // See traceback_stage.cpp for details on how its packed.

  for ( int BlastHSPResults_query_index = 0; BlastHSPResults_query_index < merged_hsp_results->num_queries ; ++BlastHSPResults_query_index )
    {    
      // Allocate the BlastHitList for this query.
      BlastHitList *merged_bhl = Blast_HitListNew(1);

      // hsplist_count
      total_hsplist_count=0;//ceb Reset here for each query
      for ( int i=0; i < num_team_leaders; ++i )
	{
          BlastHitList *bhl = (hsp_results[i])->hitlist_array[BlastHSPResults_query_index];
          dummy_int4 = bhl->hsplist_count;
	  total_hsplist_count += dummy_int4;
	  hsplist_counts[i] = dummy_int4;
	}

      merged_bhl->hsplist_count = total_hsplist_count;

      // hsplist_max
      total_hsplist_max=0;
      for ( int i=0; i < num_team_leaders; ++i )
	{
          BlastHitList *bhl = (hsp_results[i])->hitlist_array[BlastHSPResults_query_index];
	  if ( i == 0 )
	    {
	      total_hsplist_max = bhl->hsplist_max;
	    }
	}

      while ( total_hsplist_max < total_hsplist_count )
	{total_hsplist_max *= 2;}

      merged_bhl->hsplist_max = total_hsplist_max;

      // worst evalue
      worst_evalue_all = ((hsp_results[0])->hitlist_array[BlastHSPResults_query_index])->worst_evalue;
      for ( int i=1; i < num_team_leaders; ++i )
	{
          BlastHitList *bhl = (hsp_results[i])->hitlist_array[BlastHSPResults_query_index];
          dummy_double = bhl->worst_evalue;
	  if( dummy_double > worst_evalue_all )
	    {worst_evalue_all = dummy_double;}
	}
      merged_bhl->worst_evalue = worst_evalue_all;


      // best bit(low) score
      low_score_all = ((hsp_results[0])->hitlist_array[BlastHSPResults_query_index])->low_score;
      for ( int i=1; i < num_team_leaders; ++i )
	{
          BlastHitList *bhl = (hsp_results[i])->hitlist_array[BlastHSPResults_query_index];
	  dummy_int4 = bhl->low_score;
	  if( dummy_int4 < low_score_all )
	    {low_score_all = dummy_int4;}
	}
      merged_bhl->low_score = low_score_all;

      //heapified boolean
      all_heapified=0;
      for ( int i=0; i < num_team_leaders; ++i )
	{
          BlastHitList *bhl = (hsp_results[i])->hitlist_array[BlastHSPResults_query_index];
	  dummy_boolean = bhl->heapified;
	  if ( dummy_boolean == 1 )
	    {all_heapified = dummy_boolean;}
	}
      // if any are heapified assume the merge can be.
      merged_bhl->heapified = all_heapified;

      //hsplist_current
      dummy_int4 = ((hsp_results[0])->hitlist_array[BlastHSPResults_query_index])->hsplist_current;

      // take the first one or the max. add 100 until it surpasses hsplist_count!
      while ( dummy_int4 < total_hsplist_count )
	{dummy_int4 += 100;}

      total_hsplist_current = dummy_int4;
      merged_bhl->hsplist_current = total_hsplist_current;

      // Allocate the hspList to current long.
      merged_bhl->hsplist_array = (BlastHSPList**)calloc( merged_bhl->hsplist_current, sizeof(BlastHSPList*) );

      Int4 global_hsplist_index = 0;
      // Loop over each team leader and process the hit lists.
      for ( int i=0; i < num_team_leaders; ++i )
      {
        BlastHitList *bhl = (hsp_results[i])->hitlist_array[BlastHSPResults_query_index];
	for ( int BlastHitList_hsplist_index=0; BlastHitList_hsplist_index < hsplist_counts[i]; ++BlastHitList_hsplist_index )
	    {
	      // Allocate a BlastHSPList ptr
	      BlastHSPList* merged_bhspl = (BlastHSPList*) calloc(1,sizeof(BlastHSPList) );
	      BlastHSPList *bhspl = bhl->hsplist_array[BlastHitList_hsplist_index];

	      merged_bhspl->oid = bhspl->oid;
	      merged_bhspl->query_index = bhspl->query_index;
	      merged_bhspl->hspcnt = bhspl->hspcnt;
	      merged_bhspl->allocated = bhspl->allocated;
	      merged_bhspl->hsp_max = bhspl->hsp_max;
	      merged_bhspl->do_not_reallocate = bhspl->do_not_reallocate;
              merged_bhspl->best_evalue = bhspl->best_evalue;

	      // Allocate the array for High Scoring segment pairs.
	      merged_bhspl->hsp_array = (BlastHSP**) calloc( merged_bhspl->allocated, sizeof(BlastHSP*) );

	      // Loop over all HSPs between the query and the subject sequence.
	      for ( int BlastHSP_hsp_array_index=0; BlastHSP_hsp_array_index < bhspl->hspcnt; ++BlastHSP_hsp_array_index )
		{
		  BlastHSP* merged_bhsp = (BlastHSP*) calloc (1, sizeof(BlastHSP) );
		  BlastHSP* bhsp = bhspl->hsp_array[BlastHSP_hsp_array_index];

		  merged_bhsp->score = bhsp->score;
		  merged_bhsp->num_ident = bhsp->num_ident;
		  merged_bhsp->bit_score = bhsp->bit_score;
		  merged_bhsp->evalue = bhsp->evalue;
		  merged_bhsp->context = bhsp->context;
		  merged_bhsp->num = bhsp->num;
		  merged_bhsp->comp_adjustment_method = bhsp->comp_adjustment_method;
		  merged_bhsp->num_positives = bhsp->num_positives;
		  merged_bhsp->query.frame = bhsp->query.frame;
		  merged_bhsp->query.offset = bhsp->query.offset;
		  merged_bhsp->query.end = bhsp->query.end;
		  merged_bhsp->query.gapped_start = bhsp->query.gapped_start;
		  merged_bhsp->subject.frame = bhsp->subject.frame;
		  merged_bhsp->subject.offset = bhsp->subject.offset;
		  merged_bhsp->subject.end = bhsp->subject.end;
		  merged_bhsp->subject.gapped_start = bhsp->subject.gapped_start;

		  SPHIHspInfo *my_pat_info = bhsp->pat_info;
		  SPHIHspInfo *merged_pat_info = NULL;

                  if(my_pat_info!=NULL)
		  {
		    merged_pat_info = (SPHIHspInfo*) calloc( 1, sizeof(SPHIHspInfo) );
		    merged_pat_info->index = my_pat_info->index;
		    merged_pat_info->length = my_pat_info->length;
		  }
		  merged_bhsp->pat_info = merged_pat_info ;
		  
		  GapEditScript* my_gap_info = bhsp->gap_info;
		  GapEditScript *merged_gap_info = NULL;

                  if(my_gap_info != NULL)
		  {
		    merged_gap_info = (GapEditScript*) calloc(1, sizeof(GapEditScript) );
		    merged_gap_info->size = my_gap_info->size;
		    merged_gap_info->op_type = (EGapAlignOpType*) calloc( merged_gap_info->size , sizeof(EGapAlignOpType) );
		    memcpy(&(merged_gap_info->op_type[0]) , &(my_gap_info->op_type[0]) , my_gap_info->size * sizeof(EGapAlignOpType) );
		    merged_gap_info->num = (Int4*) calloc( merged_gap_info->size , sizeof(Int4) );
		    memcpy( &(merged_gap_info->num[0]) , &(my_gap_info->num[0]) , my_gap_info->size * sizeof(Int4) );
		  }
		  merged_bhsp->gap_info = merged_gap_info;

		  // Now attach this BlastHSP to the parent BlastHSPList.
		  merged_bhspl->hsp_array[ BlastHSP_hsp_array_index ] = merged_bhsp ;

		} // End of loop over the BlastHSP_hsp_array_index -> the high scoring pair data

	      // Now attach this BlastHSPList to the parent BlastHitList
	      merged_bhl->hsplist_array[ global_hsplist_index ] = merged_bhspl ;
	      ++global_hsplist_index;  // increment the global index
	      
	    }  // End loop over BlastHitList_hsplist_index

	} // End Loop over team leaders.


      // Now we can add this BlastHitList to the parent BlastHSPResults
      merged_hsp_results->hitlist_array[ BlastHSPResults_query_index ] = merged_bhl ;
  }

  // Now we can sort the HSPResults.
  Blast_HSPResultsSortByEvalue( merged_hsp_results ) ;

  return merged_hsp_results;
}








BlastHSPList* BlastHSPListDuplicate(const BlastHSPList *bhspl)
{
  // Allocate a BlastHSPList ptr
#if 0
  BlastHSPList* new_bhspl = BlastHSPListDup(bhspl);

#else  
  BlastHSPList* new_bhspl = Blast_HSPListNew(bhspl->hsp_max);
  new_bhspl->oid = bhspl->oid;
  new_bhspl->query_index = bhspl->query_index;
  new_bhspl->hspcnt = bhspl->hspcnt;
  new_bhspl->allocated = bhspl->allocated;
  new_bhspl->hsp_max = bhspl->hsp_max;
  new_bhspl->do_not_reallocate = bhspl->do_not_reallocate;
  new_bhspl->best_evalue = bhspl->best_evalue;
 
  // Allocate the array for High Scoring segment pairs.
  //new_bhspl->hsp_array = (BlastHSP**) calloc( new_bhspl->allocated, sizeof(BlastHSP*) );
  //new_bhspl->hsp_array = Blast_HSPListNew(new_bhspl->hsp_max);
  
  // Loop over all HSPs between the query and the subject sequence.
  for ( int BlastHSP_hsp_array_index=0; BlastHSP_hsp_array_index < bhspl->hspcnt; ++BlastHSP_hsp_array_index )
    {
      //BlastHSP* new_bhsp = (BlastHSP*) calloc (1, sizeof(BlastHSP) );
      //BlastHSP* new_bhsp = Blast_HSPNew();
      BlastHSP* bhsp = bhspl->hsp_array[BlastHSP_hsp_array_index];
#if 1
      //works
      BlastHSP* new_bhsp = NULL;
      /* Do not pass the edit script, because we don't want to tranfer 
	 ownership. */
      Blast_HSPInit(bhsp->query.offset, bhsp->query.end, bhsp->subject.offset, 
		    bhsp->subject.end, bhsp->query.gapped_start, 
		    bhsp->subject.gapped_start, bhsp->context, 
		    bhsp->query.frame, bhsp->subject.frame, bhsp->score, 
		    NULL, &new_bhsp);
      new_bhsp->evalue = bhsp->evalue;
      new_bhsp->num = bhsp->num;
      new_bhsp->num_ident = bhsp->num_ident;
      new_bhsp->bit_score = bhsp->bit_score;
      new_bhsp->comp_adjustment_method = bhsp->comp_adjustment_method;
      if (bhsp->gap_info) {
	new_bhsp->gap_info = GapEditScriptDup(bhsp->gap_info);
      }
      
      if (bhsp->pat_info) {
	/* Copy this HSP's pattern data. */
	new_bhsp->pat_info = 
	  (SPHIHspInfo*) BlastMemDup(bhsp->pat_info, sizeof(SPHIHspInfo));
      }
      
#else
      BlastHSP* new_bhsp = Blast_HSPNew();
      
      new_bhsp->score = bhsp->score;
      new_bhsp->num_ident = bhsp->num_ident;
      new_bhsp->bit_score = bhsp->bit_score;
      new_bhsp->evalue = bhsp->evalue;
      new_bhsp->context = bhsp->context;
      new_bhsp->num = bhsp->num;
      new_bhsp->comp_adjustment_method = bhsp->comp_adjustment_method;
      new_bhsp->num_positives = bhsp->num_positives;
      new_bhsp->query.frame = bhsp->query.frame;
      new_bhsp->query.offset = bhsp->query.offset;
      new_bhsp->query.end = bhsp->query.end;
      new_bhsp->query.gapped_start = bhsp->query.gapped_start;
      new_bhsp->subject.frame = bhsp->subject.frame;
      new_bhsp->subject.offset = bhsp->subject.offset;
      new_bhsp->subject.end = bhsp->subject.end;
      new_bhsp->subject.gapped_start = bhsp->subject.gapped_start;
      
      SPHIHspInfo *my_pat_info = bhsp->pat_info;
      SPHIHspInfo *new_pat_info = NULL;
      
      if(my_pat_info!=NULL)
	{
	  new_pat_info = (SPHIHspInfo*) calloc( 1, sizeof(SPHIHspInfo) );
	  new_pat_info->index = my_pat_info->index;
	  new_pat_info->length = my_pat_info->length;
	}
      new_bhsp->pat_info = new_pat_info ;
      
      GapEditScript* my_gap_info = bhsp->gap_info;
      GapEditScript *new_gap_info = NULL;
      
      if(my_gap_info != NULL)
	{
	  new_gap_info = (GapEditScript*) calloc(1, sizeof(GapEditScript) );
	  new_gap_info->size = my_gap_info->size;
	  new_gap_info->op_type = (EGapAlignOpType*) calloc( new_gap_info->size , sizeof(EGapAlignOpType) );
	  memcpy(&(new_gap_info->op_type[0]) , &(my_gap_info->op_type[0]) , my_gap_info->size * sizeof(EGapAlignOpType) );
	  new_gap_info->num = (Int4*) calloc( new_gap_info->size , sizeof(Int4) );
	  memcpy( &(new_gap_info->num[0]) , &(my_gap_info->num[0]) , my_gap_info->size * sizeof(Int4) );
	}
      new_bhsp->gap_info = new_gap_info;
#endif

      // Now attach this BlastHSP to the parent BlastHSPList.
      new_bhspl->hsp_array[ BlastHSP_hsp_array_index ] = new_bhsp ;
      
    } // End of loop over the BlastHSP_hsp_array_index -> the high scoring pair data
#endif  

  return new_bhspl;  
}


//======================================================================




// Functions for serializing and unserializing the objects necessary to 
// create an output file
//========================================================

char* PackHSPResultsBuffer(BlastHSPResults* HSPResults,
                           size_t& HSPResultsBSize)
{
    /*==========================================================================
      Here is where I will throw in my test code to see if I can serialize the
      BLASTHSPResults structure that was filled from the traceback search into
      a char* buffer. After filling it in, I will then deserialize into a new
      BLASTHSPResults structure and see if I get the same results.
      ==========================================================================*/

    BlastHSPResults *hsp_results=HSPResults;

    char * buffer = NULL;
    size_t buffer_size = 0;
    size_t bsize = 0;
    char * buffer_ptr = NULL;

    // hsp_results will have come out of Blast_RunTracebackSearchWithInterrupt with
    // all the results. Now lets traverse the structures.

    buffer_size += sizeof( hsp_results->num_queries );
    bsize += sizeof( Int4 );

    for ( int BlastHSPResults_query_index = 0; BlastHSPResults_query_index < hsp_results->num_queries ; ++BlastHSPResults_query_index )
      {
	BlastHitList *bhl = hsp_results->hitlist_array[BlastHSPResults_query_index];

	bsize += sizeof( Int4 );
	bsize += sizeof( Int4 );
	bsize += sizeof( double );
	bsize += sizeof( Int4 );
	bsize += sizeof( Boolean );
	bsize += sizeof( Int4 );

	for ( int BlastHitList_hsplist_index = 0 ; BlastHitList_hsplist_index < bhl->hsplist_count; ++BlastHitList_hsplist_index )
	  {
	    BlastHSPList *bhspl = bhl->hsplist_array[BlastHitList_hsplist_index];

	    bsize += sizeof( Int4 );
	    bsize += sizeof( Int4 );
	    bsize += sizeof( Int4 );
	    bsize += sizeof( Int4 );
	    bsize += sizeof( Int4 );
	    bsize += sizeof( Boolean );
	    bsize += sizeof( double );

	    for ( int BlastHSP_hsp_array_index = 0; BlastHSP_hsp_array_index < bhspl->hspcnt; ++BlastHSP_hsp_array_index )
	      {
		BlastHSP* bhsp = bhspl->hsp_array[BlastHSP_hsp_array_index];

		bsize += sizeof( Int4 );
		bsize += sizeof( Int4 );
		bsize += sizeof( double );
		bsize += sizeof( double );
		bsize += sizeof( Int4 );
		bsize += sizeof( Int4 );
		bsize += sizeof( Int2 );
		bsize += sizeof( Int4 );

		BlastSeg *bs_q = &(bhsp->query);
		BlastSeg *bs_s = &(bhsp->subject);

		SPHIHspInfo * my_pat_info = bhsp->pat_info;
		GapEditScript* my_gap_info = bhsp->gap_info;

		bsize += sizeof( Int2 );
		bsize += sizeof( Int4 );
		bsize += sizeof( Int4 );
		bsize += sizeof( Int4 );

		bsize += sizeof( Int2 );
		bsize += sizeof( Int4 );
		bsize += sizeof( Int4 );
		bsize += sizeof( Int4 );

		if ( my_pat_info )
		  {
		    bsize += sizeof( Int4 );
		    bsize += sizeof( Int4 );
		  }
		else
		  {
		    bsize += sizeof( Int4 );
		  }

		if ( my_gap_info )
		  {
		    bsize += sizeof( Int4 );
		    bsize += (sizeof( Int4 ) * my_gap_info->size );
		    bsize += (sizeof( EGapAlignOpType ) * my_gap_info->size);
		  }
		else
		  {
                    bsize += sizeof( Int4 );
		  }
	      }
	  }
      }

    // Now that I have a feel on the size of the buffer and how to walk over all the fields, 
    //   I will serialize the structures into a buffer and then deserialize
    //   them into a new set of structures and verify that the results are correct.

    buffer = (char*) calloc( bsize , sizeof(char) );

    buffer_ptr = buffer;

    // First put in the number of queries.
    memcpy ( buffer_ptr, &(hsp_results->num_queries), sizeof(Int4) );
    // Increment the buffer pointer
    buffer_ptr += sizeof(Int4);

    //std::cout<<"pack hsp_results->num_queries: "<<hsp_results->num_queries<<std::endl;

    for ( int BlastHSPResults_query_index = 0; BlastHSPResults_query_index < hsp_results->num_queries ; ++BlastHSPResults_query_index )
      {
	//std::cout<<"BlastHSPResults_query_index "<<BlastHSPResults_query_index<<std::endl;

	BlastHitList *bhl = hsp_results->hitlist_array[BlastHSPResults_query_index];

	memcpy ( buffer_ptr, &(bhl->hsplist_count), sizeof(Int4) );
	buffer_ptr += sizeof(Int4);
	//std::cout<<"bhl->hsplist_count "<<bhl->hsplist_count<<std::endl;

	memcpy ( buffer_ptr, &(bhl->hsplist_max), sizeof(Int4) );
	buffer_ptr += sizeof(Int4);
	//std::cout<<"bhl->hsplist_max "<<bhl->hsplist_max<<std::endl;
	//std::cout<<"bhl->hsplist_max "<<(double)*(buffer_ptr)<<std::endl;

	memcpy ( buffer_ptr, &(bhl->worst_evalue), sizeof(double) );
	buffer_ptr += sizeof(double);
	//std::cout<<"bhl->worst_evalue "<<bhl->worst_evalue<<std::endl;
	//std::cout<<"bhl->worst_evalue "<<(double)*(buffer_ptr)<<std::endl;

	memcpy ( buffer_ptr, &(bhl->low_score), sizeof(Int4) );
	buffer_ptr += sizeof(Int4);

	memcpy ( buffer_ptr, &(bhl->heapified), sizeof(Boolean) );
	buffer_ptr += sizeof(Boolean);

	memcpy ( buffer_ptr, &(bhl->hsplist_current), sizeof(Int4) );
	buffer_ptr += sizeof(Int4);
	//std::cout<<"bhl->hsplist_current "<<bhl->hsplist_current<<std::endl;

	for ( int BlastHitList_hsplist_index = 0 ; BlastHitList_hsplist_index < bhl->hsplist_count; ++BlastHitList_hsplist_index )
	  {
	    BlastHSPList *bhspl = bhl->hsplist_array[BlastHitList_hsplist_index];

	    memcpy ( buffer_ptr, &(bhspl->oid), sizeof(Int4) );
	    buffer_ptr += sizeof(Int4);
	    //std::cout<<"bhspl->oid "<<bhspl->oid<<std::endl;

	    memcpy ( buffer_ptr, &(bhspl->query_index), sizeof(Int4) );
	    buffer_ptr += sizeof(Int4);
	    //std::cout<<"bhspl->query_index "<<bhspl->query_index<<std::endl;

	    memcpy ( buffer_ptr, &(bhspl->hspcnt), sizeof(Int4) );
	    buffer_ptr += sizeof(Int4);
	    //std::cout<<"bhspl->hspcnt "<<bhspl->hspcnt<<std::endl;

	    memcpy ( buffer_ptr, &(bhspl->allocated), sizeof(Int4) );
	    buffer_ptr += sizeof(Int4);
	    //std::cout<<"bhspl->allocated "<<bhspl->allocated<<std::endl;

	    memcpy ( buffer_ptr, &(bhspl->hsp_max), sizeof(Int4) );
	    buffer_ptr += sizeof(Int4);
	    //std::cout<<"bhspl->hsp_max "<<bhspl->hsp_max<<std::endl;

	    memcpy ( buffer_ptr, &(bhspl->do_not_reallocate), sizeof(Boolean) );
	    buffer_ptr += sizeof(Boolean);
	    //std::cout<<"bhspl->do_not_reallocate "<<bhspl->do_not_reallocate<<std::endl;

	    memcpy ( buffer_ptr, &(bhspl->best_evalue), sizeof(double) );
	    buffer_ptr += sizeof(double);
	    //std::cout<<"bhspl->best_evalue "<<bhspl->best_evalue<<std::endl;

	    //std::cout<<"bhspl->hspcnt "<<bhspl->hspcnt<<std::endl;
	    for ( int BlastHSP_hsp_array_index = 0; BlastHSP_hsp_array_index < bhspl->hspcnt; ++BlastHSP_hsp_array_index )
	      {
		BlastHSP* bhsp = bhspl->hsp_array[BlastHSP_hsp_array_index];

		memcpy ( buffer_ptr, &(bhsp->score), sizeof(Int4) );
		buffer_ptr += sizeof(Int4);

		memcpy ( buffer_ptr, &(bhsp->num_ident), sizeof(Int4) );
		buffer_ptr += sizeof(Int4);

		memcpy ( buffer_ptr, &(bhsp->bit_score), sizeof(double) );
		buffer_ptr += sizeof(double);

		memcpy ( buffer_ptr, &(bhsp->evalue), sizeof(double) );
		buffer_ptr += sizeof(double);

		memcpy ( buffer_ptr, &(bhsp->context), sizeof(Int4) );
		buffer_ptr += sizeof(Int4);

		memcpy ( buffer_ptr, &(bhsp->num), sizeof(Int4) );
		buffer_ptr += sizeof(Int4);

		memcpy ( buffer_ptr, &(bhsp->comp_adjustment_method), sizeof(Int2) );
		buffer_ptr += sizeof(Int2);

		memcpy ( buffer_ptr, &(bhsp->num_positives), sizeof(Int4) );
		buffer_ptr += sizeof(Int4);

		BlastSeg *bs_q = &(bhsp->query);
		BlastSeg *bs_s = &(bhsp->subject);

		SPHIHspInfo * my_pat_info = bhsp->pat_info;
		GapEditScript* my_gap_info = bhsp->gap_info;

		memcpy ( buffer_ptr, &(bs_q->frame), sizeof(Int2) );
		buffer_ptr += sizeof(Int2);

		memcpy ( buffer_ptr, &(bs_q->offset), sizeof(Int4) );
		buffer_ptr += sizeof(Int4);

		memcpy ( buffer_ptr, &(bs_q->end), sizeof(Int4) );
		buffer_ptr += sizeof(Int4);

		memcpy ( buffer_ptr, &(bs_q->gapped_start), sizeof(Int4) );
		buffer_ptr += sizeof(Int4);

		memcpy ( buffer_ptr, &(bs_s->frame), sizeof(Int2) );
		buffer_ptr += sizeof(Int2);

		memcpy ( buffer_ptr, &(bs_s->offset), sizeof(Int4) );
		buffer_ptr += sizeof(Int4);

		memcpy ( buffer_ptr, &(bs_s->end), sizeof(Int4) );
		buffer_ptr += sizeof(Int4);

		memcpy ( buffer_ptr, &(bs_s->gapped_start), sizeof(Int4) );
		//std::cout<<"bs_s->gapped_start "<<bs_s->gapped_start<<std::endl;
		buffer_ptr += sizeof(Int4);


		if ( my_pat_info )
		  {
		    memcpy ( buffer_ptr, &(my_pat_info->index), sizeof(Int4) );
		    buffer_ptr += sizeof(Int4);

		    memcpy ( buffer_ptr, &(my_pat_info->length), sizeof(Int4) );
		    //std::cout<<"my_pat_info->length "<<my_pat_info->length<<std::endl;
		    buffer_ptr += sizeof(Int4);
		  }
		else
		  {
		    Int4 ZERO = 0;

		    memcpy ( buffer_ptr, &(ZERO), sizeof(Int4) );
		    buffer_ptr += sizeof(Int4);
		  }

		if ( my_gap_info )
		  {
		    memcpy ( buffer_ptr, &(my_gap_info->size), sizeof(Int4) );
		    //std::cout<<"my_gap_info->size "<<my_gap_info->size<<std::endl;
		    buffer_ptr += sizeof(Int4);

		    memcpy ( buffer_ptr, (my_gap_info->op_type), my_gap_info->size * sizeof(EGapAlignOpType) );
		    buffer_ptr += ( my_gap_info->size * sizeof(EGapAlignOpType) );

		    memcpy ( buffer_ptr, (my_gap_info->num), my_gap_info->size * sizeof(Int4) );
		    buffer_ptr += ( my_gap_info->size * sizeof(Int4) );
		  }
		else
		  {
		    Int4 ZERO = 0;

		    memcpy ( buffer_ptr, &(ZERO), sizeof(Int4) );
		    buffer_ptr += sizeof(Int4);
		  }

	      }
	  }
      }

    // Send the buffer.
    //MPI_Send(&bsize, sizeof(size_t), MPI_BYTE, 0, 1, MPI_COMM_WORLD );

    //MPI_Send(buffer, int(bsize), MPI_BYTE, 0, 1, MPI_COMM_WORLD );

    // Set class members.
    //HSPResultsBuffer = buffer;
    HSPResultsBSize = bsize;
    return buffer;
    //std::cout<<"HSPResultsBSize "<<HSPResultsBSize<<std::endl;
}






BlastHSPResults*  UnpackHSPResultsBuffer(char * HSPResultsBuffer,
                                         size_t HSPResultsBSize)
{
    /*==========================================================================
      Here is where I will throw in my test code to see if I can serialize the
      BLASTHSPResults structure that was filled from the traceback search into
      a char* buffer. After filling it in, I will then deserialize into a new
      BLASTHSPResults structure and see if I get the same results.
      ==========================================================================*/

    BlastHSPResults *my_hsp_results(0);

    char * buffer = HSPResultsBuffer;
    //size_t buffer_size = 0;

    size_t bsize = HSPResultsBSize;
    char * buffer_ptr = NULL;
    buffer_ptr=buffer;

    // Recv the packed buffer.
    //MPI_Status status;
    //MPI_Recv(&bsize, sizeof(size_t), MPI_BYTE, 1, 1, MPI_COMM_WORLD, &status);
    //buffer = (char*) calloc( bsize, sizeof(char) );
    //MPI_Recv(buffer, int(bsize), MPI_BYTE, 1, 1, MPI_COMM_WORLD, &status);

    Int4 buf_num_queries;
    memcpy( &buf_num_queries , buffer_ptr, sizeof(Int4) );
    buffer_ptr += sizeof(Int4);

    my_hsp_results = Blast_HSPResultsNew( buf_num_queries );

    for ( int BlastHSPResults_query_index = 0; BlastHSPResults_query_index < my_hsp_results->num_queries ; ++BlastHSPResults_query_index )
      {
	//std::cout<<"BlastHSPResults_query_index "<<BlastHSPResults_query_index<<std::endl;
	BlastHitList *bhl = Blast_HitListNew(1);

	memcpy( &(bhl->hsplist_count), buffer_ptr, sizeof(Int4) );
	buffer_ptr += sizeof( Int4 );

	memcpy( &(bhl->hsplist_max), buffer_ptr, sizeof(Int4) );
	buffer_ptr += sizeof( Int4 );

	memcpy( &(bhl->worst_evalue), buffer_ptr, sizeof(double) );
	buffer_ptr += sizeof( double );

	memcpy( &(bhl->low_score), buffer_ptr, sizeof(Int4) );
	buffer_ptr += sizeof( Int4 );

	memcpy( &(bhl->heapified), buffer_ptr, sizeof(Boolean) );
	buffer_ptr += sizeof( Boolean );

	memcpy( &(bhl->hsplist_current), buffer_ptr, sizeof(Int4) );
	buffer_ptr += sizeof( Int4 );
	//std::cout<<"bhl->hsplist_current "<<bhl->hsplist_current<<std::endl;

	// Allocate the internal BlastHSPList array.
	bhl->hsplist_array = (BlastHSPList**) calloc( bhl->hsplist_current , sizeof(BlastHSPList*) );

	if ( bhl->hsplist_current < bhl->hsplist_count )
	  {
	    std::cout<<"For the blasthitlist the count was greater than the current value!"<<std::endl;
	  }

	for ( int BlastHitList_hsplist_index = 0 ; BlastHitList_hsplist_index < bhl->hsplist_count; ++BlastHitList_hsplist_index )
	  {
	    BlastHSPList* bhspl = (BlastHSPList*) calloc( 1, sizeof(BlastHSPList) );

	    memcpy( &(bhspl->oid) , buffer_ptr, sizeof(Int4) );
	    buffer_ptr += sizeof(Int4);

	    memcpy( &(bhspl->query_index) , buffer_ptr, sizeof(Int4) );
	    buffer_ptr += sizeof(Int4);

	    memcpy( &(bhspl->hspcnt) , buffer_ptr, sizeof(Int4) );
	    buffer_ptr += sizeof(Int4);

	    memcpy( &(bhspl->allocated) , buffer_ptr, sizeof(Int4) );
	    buffer_ptr += sizeof(Int4);

	    memcpy( &(bhspl->hsp_max) , buffer_ptr, sizeof(Int4) );
	    buffer_ptr += sizeof(Int4);

	    memcpy( &(bhspl->do_not_reallocate) , buffer_ptr, sizeof(Boolean) );
	    buffer_ptr += sizeof(Boolean);

	    memcpy( &(bhspl->best_evalue) , buffer_ptr, sizeof(double) );
	    buffer_ptr += sizeof(double);

	    bhspl->hsp_array = (BlastHSP**) calloc( bhspl->allocated, sizeof(BlastHSP*) );

	    if ( bhspl->allocated < bhspl->hspcnt )
	      std::cout<<"  for the blasthsplist the allocated was smaller than the hspcnt!"<<std::endl;

	    for ( int BlastHSP_hsp_array_index=0; BlastHSP_hsp_array_index < bhspl->hspcnt; ++BlastHSP_hsp_array_index )
	      {
		BlastHSP* bhsp = (BlastHSP*) calloc (1, sizeof(BlastHSP) );

		memcpy ( &(bhsp->score), buffer_ptr , sizeof(Int4) );
		buffer_ptr += sizeof(Int4);

		memcpy ( &(bhsp->num_ident), buffer_ptr , sizeof(Int4) );
		buffer_ptr += sizeof(Int4);

		memcpy ( &(bhsp->bit_score), buffer_ptr , sizeof(double) );
		buffer_ptr += sizeof(double);

		memcpy ( &(bhsp->evalue), buffer_ptr , sizeof(double) );
		buffer_ptr += sizeof(double);

		memcpy ( &(bhsp->context), buffer_ptr , sizeof(Int4) );
		buffer_ptr += sizeof(Int4);

		memcpy ( &(bhsp->num), buffer_ptr , sizeof(Int4) );
		buffer_ptr += sizeof(Int4);

		memcpy ( &(bhsp->comp_adjustment_method), buffer_ptr , sizeof(Int2) );
		buffer_ptr += sizeof(Int2);

		memcpy ( &(bhsp->num_positives), buffer_ptr , sizeof(Int4) );
		buffer_ptr += sizeof(Int4);

		memcpy ( &(bhsp->query.frame), buffer_ptr , sizeof(Int2) );
		buffer_ptr += sizeof(Int2);

		memcpy ( &(bhsp->query.offset), buffer_ptr , sizeof(Int4) );
		buffer_ptr += sizeof(Int4);

		memcpy ( &(bhsp->query.end), buffer_ptr , sizeof(Int4) );
		buffer_ptr += sizeof(Int4);

		memcpy ( &(bhsp->query.gapped_start), buffer_ptr , sizeof(Int4) );
		buffer_ptr += sizeof(Int4);

		memcpy ( &(bhsp->subject.frame), buffer_ptr , sizeof(Int2) );
		buffer_ptr += sizeof(Int2);

		memcpy ( &(bhsp->subject.offset), buffer_ptr , sizeof(Int4) );
		buffer_ptr += sizeof(Int4);

		memcpy ( &(bhsp->subject.end), buffer_ptr , sizeof(Int4) );
		buffer_ptr += sizeof(Int4);

		memcpy ( &(bhsp->subject.gapped_start), buffer_ptr , sizeof(Int4) );
		buffer_ptr += sizeof(Int4);

		// The next value might be a setinel value.
		Int4 sentinel;

		memcpy ( &sentinel , buffer_ptr, sizeof(Int4) );

		if ( sentinel )
		  {
		    SPHIHspInfo *my_pat_info = (SPHIHspInfo*) calloc( 1, sizeof(SPHIHspInfo) );
		    memcpy( &(my_pat_info->index) , buffer_ptr, sizeof(Int4) );
		    buffer_ptr += sizeof(Int4);
		    memcpy( &(my_pat_info->length) , buffer_ptr, sizeof(Int4) );
		    buffer_ptr += sizeof(Int4);
		    bhsp->pat_info = my_pat_info ;
		  }
		else
		  {
		    bhsp->pat_info = NULL;
		    buffer_ptr += sizeof(Int4);
		  }

		// The next value might be a sentinel value on GapEditScript.
		memcpy ( &sentinel , buffer_ptr , sizeof(Int4) );

		if ( sentinel )
		  {
		    GapEditScript *my_gap_info = (GapEditScript*) calloc(1, sizeof(GapEditScript) );

		    memcpy( &(my_gap_info->size) , buffer_ptr , sizeof(Int4) );
		    buffer_ptr += sizeof(Int4);

		    my_gap_info->op_type = (EGapAlignOpType*) calloc( my_gap_info->size , sizeof(EGapAlignOpType) );
		    memcpy( &(my_gap_info->op_type[0]), buffer_ptr , my_gap_info->size * sizeof(EGapAlignOpType) );
		    buffer_ptr += ( my_gap_info->size * sizeof(EGapAlignOpType) );

		    my_gap_info->num = (Int4*) calloc( my_gap_info->size , sizeof(Int4) );
		    memcpy( &(my_gap_info->num[0]), buffer_ptr , my_gap_info->size * sizeof(Int4) );
		    buffer_ptr += (my_gap_info->size * sizeof(Int4) );

		    bhsp->gap_info = my_gap_info ;
		  }
		else
		  {
		    bhsp->gap_info = NULL;

		    buffer_ptr += sizeof(Int4);
		  }
		// Now attach this BlastHSP to the parent BlastHSPList.
		bhspl->hsp_array[ BlastHSP_hsp_array_index ] = bhsp ;
	      }
	    // Now attach this BlastHSPList to the parent BlastHitList.
	    bhl->hsplist_array[ BlastHitList_hsplist_index ] = bhspl ;
	  }
	// Now attach this BlastHitList to the parent BlastHSPResults
	my_hsp_results->hitlist_array[ BlastHSPResults_query_index ] = bhl ;
      }
    // Switch between my results and the original results.
    return my_hsp_results;

    // Clean up the mess.
    //delete [] buffer;
}




char* PackScoreBlkBuffer(BlastScoreBlk* ScoreBlk,
                         size_t& ScoreBlkBufferSize)
{
    // ============================================================================
    // Lets dump the m_InternalData->BlastScoreBlk structure.

    BlastScoreBlk* sbp = ScoreBlk;

    char *buffer = NULL;
    char *buff_ptr = NULL;
    size_t bsize = 0;

    // Pack data first
    bsize += sizeof(Boolean);
    bsize += sizeof(Uint1);
    bsize += sizeof(Int2);
    bsize += sizeof(Int2);
    std::string sbp_name = sbp->name;
    int str_length = strlen(sbp_name.c_str());
    bsize += sizeof(int);   // store the length of the string for the matrix name.
    bsize += ( str_length * sizeof(char) );

    // Let's walk over the comments.
    ListNode *newnode = sbp->comments;

    // count up the number of comments first.
    int num_items = 0;

    if ( newnode )
      {
	num_items++;
	newnode = newnode->next;

	while ( newnode )
	  {
	    num_items++;
	    newnode = newnode->next;
	  }
      }

    newnode = sbp->comments;
    bsize += sizeof(int);   // store the number of comments first.

    if ( newnode )
      {
	bsize += sizeof( Uint1 );

	if ( newnode->ptr )
	  {
	    sbp_name = (char*)newnode->ptr;
	    str_length = strlen( sbp_name.c_str() );
	    bsize += sizeof(int);        // store the length of the string.
	    bsize += ( str_length * sizeof (char) );
	  }

	newnode = newnode->next;

	while ( newnode->next != NULL )
	  {
	    if ( newnode->ptr )
	      {
		sbp_name = (char*)newnode->ptr;
		bsize += sizeof( Uint1 );
		str_length = strlen( sbp_name.c_str() );
		bsize += sizeof(int);
		bsize += ( str_length * sizeof (char) );
	      }
	    newnode = newnode->next;
	  }
      }

    // store a sentinel value for presence of matrix
    bsize += sizeof(int);

    if ( sbp->matrix )
      {
	bsize += sizeof( size_t );
	bsize += sizeof( size_t );
	bsize += ( sbp->matrix->nrows * sbp->matrix->ncols * sizeof(int) ); // Be careful of overflow
	bsize += sizeof(double);
	bsize += ( sbp->alphabet_size * sizeof(double) );
      }

    // store a sentinel value for presence of psi_matrix
    bsize += sizeof(int);

    if ( sbp->psi_matrix )
      {
	bsize += sizeof( size_t );
	bsize += sizeof( size_t );
	bsize += ( sbp->psi_matrix->pssm->nrows * sbp->psi_matrix->pssm->ncols * sizeof(int) ); // Be careful of overflow
	bsize += sizeof(double);
	bsize += ( sbp->alphabet_size * sizeof( double) );
	bsize += ( sbp->psi_matrix->pssm->nrows * sbp->psi_matrix->pssm->ncols * sizeof(double) ); // Be careful of overflow
	bsize += ( 5 * sizeof(double) );
      }

    bsize += sizeof(Boolean);
    bsize += sizeof(Boolean);
    bsize += sizeof(Int4);
    bsize += sizeof(Int4);
    bsize += sizeof(Int4);
    bsize += sizeof(Int4);
    bsize += sizeof(double);
    bsize += sizeof(Boolean);

    // go ahead and count number of contexts here since this is where we will save it in the buffer.
    bsize += sizeof(Int4);

    // store a sentinel value for presence of sblast score block
    bsize += sizeof(int);

    if ( sbp->sfp )
      {
	for ( Int4 ncont = 0; ncont < sbp->number_of_contexts; ++ncont )
	  {
	    Blast_ScoreFreq* my_sfp = NULL;
	    my_sfp = sbp->sfp[ncont];

	    if ( my_sfp )
	      {
		bsize += ( 4 * sizeof(Int4) );
		bsize += sizeof(double);
		int score_range = my_sfp->score_max - my_sfp->score_min + 1;
		bsize += ( score_range * sizeof(double) );
	      }
	  }
      }

    // store a sentinel value for presence of karlin stats block
    // ** THIS IS AN ALIAS PTR **
    // 0 means points to NULL, 1 means points to kbp_std, 2 means points to kbp_psi
    bsize += sizeof(int);

    // store a sentinel value for presence of Karlin Stats Block for gapped alignments
    // ** THIS IS AN ALIAS PTR **
    // 0 means points to NULL, 1 means points to kbp_gap_std, 2 means points to kbp_gap_psi
    bsize += sizeof(int);

    // store a sentinel value for presence of Gumbel block
    bsize += sizeof(int);

    if ( sbp->gbp )
      {
	bsize += ( 11 * sizeof(double) );
	bsize += sizeof(Int8);
	bsize += sizeof(Boolean);
      }

    // store a sentinel value for presence of karlin stats block for ungapped alignment
    bsize += sizeof(int);

    if ( sbp->kbp_std )
      {
	for ( Int4 ncont = 0; ncont < sbp->number_of_contexts; ++ncont )
	  {
	    Blast_KarlinBlk* my_kbp_std = NULL;
	    my_kbp_std = sbp->kbp_std[ncont];

	    if ( my_kbp_std )
	      {
		bsize += ( 5 * sizeof(double) );
	      }
	  }
      }

    // store a sentinel value for presence of karlin stats block for position based alignments
    bsize += sizeof(int);

    if ( sbp->kbp_psi )
      {
	for ( Int4 ncont = 0; ncont < sbp->number_of_contexts; ++ncont )
	  {
	    Blast_KarlinBlk* my_kbp_psi = NULL;
	    my_kbp_psi = sbp->kbp_psi[ncont];

	    if ( my_kbp_psi )
	      {
		bsize += ( 5 * sizeof(double) );
	      }
	  }
      }

    // store a sentinel value for presence of std karlin stats block
    bsize += sizeof(int);

    if ( sbp->kbp_gap_std )
      {
	for ( Int4 ncont = 0; ncont < sbp->number_of_contexts; ++ncont )
	  {
	    Blast_KarlinBlk* my_kbp_gap_std = NULL;
	    my_kbp_gap_std = sbp->kbp_gap_std[ncont];

	    if ( my_kbp_gap_std )
	      {
		bsize += ( 5 * sizeof(double) );
	      }
	  }
      }

    // store a sentinel value for presence of position based std karlin stats block
    bsize += sizeof(int);

    if ( sbp->kbp_gap_psi )
      {
	for ( Int4 ncont = 0; ncont < sbp->number_of_contexts; ++ncont )
	  {
	    Blast_KarlinBlk* my_kbp_gap_psi = NULL;
	    my_kbp_gap_psi = sbp->kbp_gap_psi[ncont];

	    if ( my_kbp_gap_psi )
	      {
		bsize += ( 5 * sizeof(double) );
	      }
	  }
      }

    // store a sentinel value for presence of ideal karlin stats block
    bsize += sizeof(int);

    if ( sbp->kbp_ideal )
      {
	Blast_KarlinBlk* my_kbp_ideal = sbp->kbp_ideal;

	if ( my_kbp_ideal )
	  {
	    bsize += ( 5 * sizeof(double) );
	  }

      }

    // store a sentinel value for presence of ambiguous symbol resolution array
    bsize += sizeof(int);

    if ( sbp->ambiguous_res )
      {
	bsize += ( sbp->ambig_size * sizeof(Uint1) );
      }

    // THESE WILL BE SERIALIZED BEFORE THE ARRAY SO THAT WE CAN SIZE ON A DIFFERENT PROCESS.
    bsize += sizeof(Int2);
    bsize += sizeof(Int2);
    bsize += sizeof(Boolean);


    // Now serialize the data.

    buffer = (char*)calloc(bsize, sizeof(char));
    buff_ptr = buffer;

    memcpy( buff_ptr, &(sbp->protein_alphabet), sizeof(Boolean));
    buff_ptr += sizeof(Boolean);

    memcpy( buff_ptr, &(sbp->alphabet_code), sizeof(Uint1));
    buff_ptr += sizeof(Uint1);

    memcpy( buff_ptr, &(sbp->alphabet_size), sizeof(Int2));
    buff_ptr += sizeof(Int2);

    memcpy( buff_ptr, &(sbp->alphabet_start), sizeof(Int2));
    buff_ptr += sizeof(Int2);

    // Store the length of the string ahead of the string.
    str_length = strlen(sbp_name.c_str());
    memcpy( buff_ptr, &(str_length), sizeof(int));
    buff_ptr += sizeof(int);

    memcpy( buff_ptr, sbp->name, str_length*sizeof(char) );
    buff_ptr += (str_length * sizeof(char) );

    // Pack up the comments ListNode structure. We will assume that the only data stored are strings.
    newnode = sbp->comments;

    // First, we need to know how many items are in the linked list given by the ListNode structure.
    num_items = 0;

    if ( newnode )
      {
	num_items++;
	newnode = newnode->next;

	while ( newnode )
	  {
	    num_items++;
	    newnode = newnode->next;
	  }
      }
    memcpy( buff_ptr, &num_items, sizeof(int) );
    buff_ptr += sizeof(int);

    if ( newnode )
      {
	memcpy( buff_ptr, &(newnode->choice), sizeof(Uint1) );
	buff_ptr += sizeof(Uint1);

	// Store the string length ahead of the actual string.
	sbp_name = (char*)newnode->ptr;
	str_length = strlen( sbp_name.c_str() );
	bsize += ( sizeof(int) );
	memcpy( buff_ptr, newnode->ptr, str_length*sizeof(char) );
	buff_ptr += ( str_length * sizeof(char) );



	newnode = newnode->next;

	while ( newnode )
	  {
	    memcpy( buff_ptr, &(newnode->choice), sizeof(Uint1) );
	    buff_ptr += sizeof(Uint1);

	    // Store the string length ahead of the actual string.
	    sbp_name = (char*)newnode->ptr;
	    str_length = strlen( sbp_name.c_str() );
	    bsize += ( sizeof(int) );
	    memcpy( buff_ptr, newnode->ptr, str_length*sizeof(char) );
	    buff_ptr += ( str_length * sizeof(char) );
	    newnode = newnode->next;
	  }
      }

    // Now serialize the SBlastScoreMatrix.
    if ( sbp->matrix )
      {
	num_items = 1; // has matrix
	memcpy( buff_ptr, &num_items, sizeof(int) );
	buff_ptr += sizeof(int);

	memcpy( buff_ptr, &(sbp->matrix->ncols), sizeof(size_t) );
	buff_ptr += sizeof(size_t);

	memcpy( buff_ptr, &(sbp->matrix->nrows), sizeof(size_t) );
	buff_ptr += sizeof(size_t);

	for ( size_t this_row=0; this_row < sbp->matrix->nrows; ++this_row )
	  {
	    for ( size_t this_col=0; this_col < sbp->matrix->ncols; ++this_col )
	      {
		memcpy( buff_ptr, &(sbp->matrix->data[this_row][this_col]), sizeof(int) );
		buff_ptr += sizeof(int);
	      }
	  }

	for ( Int2 freqs_ind = 0; freqs_ind < sbp->alphabet_size ; ++freqs_ind )
	  {
	    memcpy( buff_ptr, &(sbp->matrix->freqs[freqs_ind]), sizeof(double) );
	    buff_ptr += sizeof(double);
	  }

	memcpy( buff_ptr, &(sbp->matrix->lambda), sizeof(double) );
	buff_ptr += sizeof(double);
      }
    else
      {
	num_items = 0; //no matrix.
	memcpy( buff_ptr, &num_items, sizeof(int) );
	buff_ptr += sizeof(int);
      }

    // now we store the SPsiBlastScoreMatrix -- first we'll save its SBlastScoreMatrix -- but only if its been allocated

    if ( sbp->psi_matrix )
      {
	num_items = 1; //have psi_matrix.
	memcpy( buff_ptr, &num_items, sizeof(int) );
	buff_ptr += sizeof(int);

	memcpy( buff_ptr, &(sbp->psi_matrix->pssm->ncols), sizeof(size_t) );
	buff_ptr += sizeof(size_t);

	memcpy( buff_ptr, &(sbp->psi_matrix->pssm->nrows), sizeof(size_t) );
	buff_ptr += sizeof(size_t);

	for ( size_t this_row=0; this_row < sbp->psi_matrix->pssm->nrows; ++this_row )
	  {
	    for ( size_t this_col=0; this_col < sbp->psi_matrix->pssm->ncols; ++this_col )
	      {
		memcpy( buff_ptr, &(sbp->psi_matrix->pssm->data[this_row][this_col]), sizeof(int) );
		buff_ptr += sizeof(int);
	      }
	  }

	for ( Int2 freqs_ind = 0; freqs_ind < sbp->alphabet_size ; ++freqs_ind )
	  {
	    memcpy( buff_ptr, &(sbp->psi_matrix->pssm->freqs[freqs_ind]), sizeof(double) );
	    buff_ptr += sizeof(double);
	  }

	memcpy( buff_ptr, &(sbp->psi_matrix->pssm->lambda), sizeof(double) );
	buff_ptr += sizeof(double);

	// now store the frequency ratios table
	for ( size_t this_row=0; this_row < sbp->psi_matrix->pssm->nrows; ++this_row )
	  {
	    for ( size_t this_col=0; this_col < sbp->psi_matrix->pssm->ncols; ++this_col )
	      {
		memcpy( buff_ptr, &(sbp->psi_matrix->freq_ratios[this_row][this_col]), sizeof(double) );
		buff_ptr += sizeof(double);
	      }
	  }

	// and lastly the associated Blast_KarlinBlk
	memcpy( buff_ptr, &(sbp->psi_matrix->kbp->Lambda), sizeof(double) );
	buff_ptr += sizeof(double);

	memcpy( buff_ptr, &(sbp->psi_matrix->kbp->K), sizeof(double) );
	buff_ptr += sizeof(double);

	memcpy( buff_ptr, &(sbp->psi_matrix->kbp->logK), sizeof(double) );
	buff_ptr += sizeof(double);

	memcpy( buff_ptr, &(sbp->psi_matrix->kbp->H), sizeof(double) );
	buff_ptr += sizeof(double);

	memcpy( buff_ptr, &(sbp->psi_matrix->kbp->paramC), sizeof(double) );
	buff_ptr += sizeof(double);
      }
    else
      {
	num_items = 0; //no matrix.
	memcpy( buff_ptr, &num_items, sizeof(int) );
	buff_ptr += sizeof(int);
      }

    // more POD members.
    memcpy( buff_ptr, &(sbp->matrix_only_scoring), sizeof(Boolean) );
    buff_ptr += sizeof(Boolean);

    memcpy( buff_ptr, &(sbp->complexity_adjusted_scoring), sizeof(Boolean) );
    buff_ptr += sizeof(Boolean);

    memcpy( buff_ptr, &(sbp->loscore), sizeof(Int4) );
    buff_ptr += sizeof(Int4);

    memcpy( buff_ptr, &(sbp->hiscore), sizeof(Int4) );
    buff_ptr += sizeof(Int4);

    memcpy( buff_ptr, &(sbp->penalty), sizeof(Int4) );
    buff_ptr += sizeof(Int4);

    memcpy( buff_ptr, &(sbp->reward), sizeof(Int4) );
    buff_ptr += sizeof(Int4);

    memcpy( buff_ptr, &(sbp->scale_factor), sizeof(double) );
    buff_ptr += sizeof(double);

    memcpy( buff_ptr, &(sbp->read_in_matrix), sizeof(Boolean) );
    buff_ptr += sizeof(Boolean);

    // Now we loop over the Blast_ScoreFreq array and write each one first.
    // --- Need the 'number_of_contexts' variable to parse out much of the structures below this point so record it now.

    memcpy( buff_ptr, &(sbp->number_of_contexts), sizeof(Int4) );
    buff_ptr += sizeof(Int4);

    if ( sbp->sfp )
      {
	num_items = 1; //has matrix.
	memcpy( buff_ptr, &num_items, sizeof(int) );
	buff_ptr += sizeof(int);

	Blast_ScoreFreq* my_sfp = NULL;
	for ( Int4 ncont=0; ncont < sbp->number_of_contexts; ++ncont )
	  {
	    my_sfp = sbp->sfp[ncont];

	    memcpy( buff_ptr, &(my_sfp->score_min), sizeof(Int4) );
	    buff_ptr += sizeof(Int4);

	    memcpy( buff_ptr, &(my_sfp->score_max), sizeof(Int4) );
	    buff_ptr += sizeof(Int4);

	    memcpy( buff_ptr, &(my_sfp->obs_min), sizeof(Int4) );
	    buff_ptr += sizeof(Int4);

	    memcpy( buff_ptr, &(my_sfp->obs_max), sizeof(Int4) );
	    buff_ptr += sizeof(Int4);

	    memcpy( buff_ptr, &(my_sfp->score_avg), sizeof(double) );
	    buff_ptr += sizeof(double);

	    int score_range = my_sfp->score_max - my_sfp->score_min + 1;
	    for ( int score_index = 0; score_index < score_range; ++score_index )
	      {
		memcpy( buff_ptr, &(my_sfp->sprob0[score_index]), sizeof(double) );
		buff_ptr += sizeof(double);
	      }
	  }

	// remember to shift sprob when I make the new structure.
      }
    else
      {
	num_items = 0; //no matrix.
	memcpy( buff_ptr, &num_items, sizeof(int) );
	buff_ptr += sizeof(int);
      }

    // serialize the kbp structure.
    if ( sbp->gbp )
      {
	num_items = 1; //has matrix.
	memcpy( buff_ptr, &num_items, sizeof(int) );
	buff_ptr += sizeof(int);

	Blast_GumbelBlk* my_gbp = NULL;
	my_gbp = sbp->gbp;

	memcpy( buff_ptr, &(my_gbp->Lambda), sizeof(double) );
	buff_ptr += sizeof(double);

	memcpy( buff_ptr, &(my_gbp->C), sizeof(double) );
	buff_ptr += sizeof(double);

	memcpy( buff_ptr, &(my_gbp->G), sizeof(double) );
	buff_ptr += sizeof(double);

	memcpy( buff_ptr, &(my_gbp->a), sizeof(double) );
	buff_ptr += sizeof(double);

	memcpy( buff_ptr, &(my_gbp->Alpha), sizeof(double) );
	buff_ptr += sizeof(double);

	memcpy( buff_ptr, &(my_gbp->Sigma), sizeof(double) );
	buff_ptr += sizeof(double);

	memcpy( buff_ptr, &(my_gbp->a_un), sizeof(double) );
	buff_ptr += sizeof(double);

	memcpy( buff_ptr, &(my_gbp->Alpha_un), sizeof(double) );
	buff_ptr += sizeof(double);

	memcpy( buff_ptr, &(my_gbp->b), sizeof(double) );
	buff_ptr += sizeof(double);

	memcpy( buff_ptr, &(my_gbp->Beta), sizeof(double) );
	buff_ptr += sizeof(double);

	memcpy( buff_ptr, &(my_gbp->Tau), sizeof(double) );
	buff_ptr += sizeof(double);

	memcpy( buff_ptr, &(my_gbp->db_length), sizeof(Int8) );
	buff_ptr += sizeof(Int8);

	memcpy( buff_ptr, &(my_gbp->filled), sizeof(Boolean) );
	buff_ptr += sizeof(Boolean);
      }
    else
      {
	num_items = 0; //no matrix.
	memcpy( buff_ptr, &num_items, sizeof(int) );
	buff_ptr += sizeof(int);
      }

    // serialize the kbp_std structure.
    if ( sbp->kbp_std)
      {
	num_items = 1; //has matrix.
	memcpy( buff_ptr, &num_items, sizeof(int) );
	buff_ptr += sizeof(int);

	for ( Int4 ncont = 0; ncont < sbp->number_of_contexts; ++ncont )
	  {
	    Blast_KarlinBlk* my_kbp_std = NULL;
	    my_kbp_std = sbp->kbp_std[ncont];

	    memcpy( buff_ptr, &(my_kbp_std->Lambda), sizeof(double) );
	    buff_ptr += sizeof(double);

	    memcpy( buff_ptr, &(my_kbp_std->K), sizeof(double) );
	    buff_ptr += sizeof(double);

	    memcpy( buff_ptr, &(my_kbp_std->logK), sizeof(double) );
	    buff_ptr += sizeof(double);

	    memcpy( buff_ptr, &(my_kbp_std->H), sizeof(double) );
	    buff_ptr += sizeof(double);

	    memcpy( buff_ptr, &(my_kbp_std->paramC), sizeof(double) );
	    buff_ptr += sizeof(double);
	  }
      }
    else
      {
	num_items = 0; //no matrix.
	memcpy( buff_ptr, &num_items, sizeof(int) );
	buff_ptr += sizeof(int);
      }

    // serialize the kbp_psi structure.
    if ( sbp->kbp_psi)
      {
	num_items = 1; //has matrix.
	memcpy( buff_ptr, &num_items, sizeof(int) );
	buff_ptr += sizeof(int);

	for ( Int4 ncont = 0; ncont < sbp->number_of_contexts; ++ncont )
	  {
	    Blast_KarlinBlk* my_kbp_psi = NULL;
	    my_kbp_psi = sbp->kbp_psi[ncont];

	    memcpy( buff_ptr, &(my_kbp_psi->Lambda), sizeof(double) );
	    buff_ptr += sizeof(double);

	    memcpy( buff_ptr, &(my_kbp_psi->K), sizeof(double) );
	    buff_ptr += sizeof(double);

	    memcpy( buff_ptr, &(my_kbp_psi->logK), sizeof(double) );
	    buff_ptr += sizeof(double);

	    memcpy( buff_ptr, &(my_kbp_psi->H), sizeof(double) );
	    buff_ptr += sizeof(double);

	    memcpy( buff_ptr, &(my_kbp_psi->paramC), sizeof(double) );
	    buff_ptr += sizeof(double);
	  }
      }
    else
      {
	num_items = 0; //no matrix.
	memcpy( buff_ptr, &num_items, sizeof(int) );
	buff_ptr += sizeof(int);
      }

    // serialize the kbp_gap_std structure.
    if ( sbp->kbp_gap_std)
      {
	num_items = 1; //has matrix.
	memcpy( buff_ptr, &num_items, sizeof(int) );
	buff_ptr += sizeof(int);

	for ( Int4 ncont = 0; ncont < sbp->number_of_contexts; ++ncont )
	  {
	    Blast_KarlinBlk* my_kbp_gap_std = NULL;
	    my_kbp_gap_std = sbp->kbp_gap_std[ncont];

	    memcpy( buff_ptr, &(my_kbp_gap_std->Lambda), sizeof(double) );
	    buff_ptr += sizeof(double);

	    memcpy( buff_ptr, &(my_kbp_gap_std->K), sizeof(double) );
	    buff_ptr += sizeof(double);

	    memcpy( buff_ptr, &(my_kbp_gap_std->logK), sizeof(double) );
	    buff_ptr += sizeof(double);

	    memcpy( buff_ptr, &(my_kbp_gap_std->H), sizeof(double) );
	    buff_ptr += sizeof(double);

	    memcpy( buff_ptr, &(my_kbp_gap_std->paramC), sizeof(double) );
	    buff_ptr += sizeof(double);
	  }
      }
    else
      {
	num_items = 0; //no matrix.
	memcpy( buff_ptr, &num_items, sizeof(int) );
	buff_ptr += sizeof(int);
      }

    // serialize the kbp_gap_psi structure.
    if ( sbp->kbp_gap_psi)
      {
	num_items = 1; //has matrix.
	memcpy( buff_ptr, &num_items, sizeof(int) );
	buff_ptr += sizeof(int);

	for ( Int4 ncont = 0; ncont < sbp->number_of_contexts; ++ncont )
	  {
	    Blast_KarlinBlk* my_kbp_gap_psi = NULL;
	    my_kbp_gap_psi = sbp->kbp_gap_psi[ncont];

	    memcpy( buff_ptr, &(my_kbp_gap_psi->Lambda), sizeof(double) );
	    buff_ptr += sizeof(double);

	    memcpy( buff_ptr, &(my_kbp_gap_psi->K), sizeof(double) );
	    buff_ptr += sizeof(double);

	    memcpy( buff_ptr, &(my_kbp_gap_psi->logK), sizeof(double) );
	    buff_ptr += sizeof(double);

	    memcpy( buff_ptr, &(my_kbp_gap_psi->H), sizeof(double) );
	    buff_ptr += sizeof(double);

	    memcpy( buff_ptr, &(my_kbp_gap_psi->paramC), sizeof(double) );
	    buff_ptr += sizeof(double);

	  }
      }
    else
      {
	num_items = 0; //no matrix.
	memcpy( buff_ptr, &num_items, sizeof(int) );
	buff_ptr += sizeof(int);
      }

    // serialize the kbp_ideal structure.
    if ( sbp->kbp_ideal)
      {
	num_items = 1; //has matrix.
	memcpy( buff_ptr, &num_items, sizeof(int) );
	buff_ptr += sizeof(int);

	Blast_KarlinBlk* my_kbp_ideal = NULL;
	my_kbp_ideal = sbp->kbp_ideal;

	memcpy( buff_ptr, &(my_kbp_ideal->Lambda), sizeof(double) );
	buff_ptr += sizeof(double);

	memcpy( buff_ptr, &(my_kbp_ideal->K), sizeof(double) );
	buff_ptr += sizeof(double);

	memcpy( buff_ptr, &(my_kbp_ideal->logK), sizeof(double) );
	buff_ptr += sizeof(double);

	memcpy( buff_ptr, &(my_kbp_ideal->H), sizeof(double) );
	buff_ptr += sizeof(double);

	memcpy( buff_ptr, &(my_kbp_ideal->paramC), sizeof(double) );
	buff_ptr += sizeof(double);
      }
    else
      {
	num_items = 0; //no matrix.
	memcpy( buff_ptr, &num_items, sizeof(int) );
	buff_ptr += sizeof(int);
      }

    if ( sbp->kbp)
      {
	if ( ( (void*)sbp->kbp ) == ( (void*)sbp->kbp_std ) )
	  {
	    num_items = 1; // points to kbp_std
	  }
	else
	  {
	    num_items = 2; // points to kbp_psi
	  }

	memcpy( buff_ptr, &num_items, sizeof(int) );
	buff_ptr += sizeof(int);
      }
    else
      {
	num_items = 0; //no matrix.
	memcpy( buff_ptr, &num_items, sizeof(int) );
	buff_ptr += sizeof(int);
      }

    if ( sbp->kbp_gap )
      {
	if ( ( (void*)sbp->kbp_gap ) == ( (void*)sbp->kbp_gap_std ) )
	  {
	    num_items = 1; // points to kbp_gap_std
	  }
	else
	  {
	    num_items = 2; // points to kbp_gap_psi
	  }

	memcpy( buff_ptr, &num_items, sizeof(int) );
	buff_ptr += sizeof(int);
      }
    else
      {
	num_items = 0; //no matrix.
	memcpy( buff_ptr, &num_items, sizeof(int) );
	buff_ptr += sizeof(int);
      }

    memcpy( buff_ptr, &(sbp->ambig_size), sizeof(Int2) );
    buff_ptr += sizeof(Int2);

    memcpy( buff_ptr, &(sbp->ambig_occupy), sizeof(Int2) );
    buff_ptr += sizeof(Int2);

    if ( sbp->ambiguous_res )
      {
	num_items = 1; //has matrix.
	memcpy( buff_ptr, &num_items, sizeof(int) );
	buff_ptr += sizeof(int);

	for ( int ambig_index = 0; ambig_index < sbp->ambig_size; ++ambig_index )
	  {
	    memcpy( buff_ptr, &(sbp->ambiguous_res[ambig_index]), sizeof(Uint1) );
	    buff_ptr += sizeof(Uint1);
	  }
      }
    else
      {
	num_items = 0; //no matrix.
	memcpy( buff_ptr, &num_items, sizeof(int) );
	buff_ptr += sizeof(int);
      }

    memcpy( buff_ptr, &(sbp->round_down), sizeof(Boolean) );
    buff_ptr += sizeof(Boolean);

    // Send the buffer.
    //MPI_Send(&bsize, sizeof(size_t), MPI_BYTE, 0, 1, MPI_COMM_WORLD );

    //MPI_Send(buffer, int(bsize), MPI_BYTE, 0, 1, MPI_COMM_WORLD );

    // Delete the buffer.
    //free(buffer);

    // Save the BlastScoreBlk to the object's data member.
    ScoreBlkBufferSize = bsize;
    
    return buffer;
}




BlastScoreBlk* UnpackScoreBlkBuffer(char * ScoreBlkBuffer,
                                    size_t ScoreBlkBufferSize)
{
  // ============================================================================
  // Lets recv the m_InternalData->BlastScoreBlk structure.

  char *buffer = ScoreBlkBuffer;
  char *buff_ptr = NULL;
  size_t bsize = ScoreBlkBufferSize;
  int str_length = 0;
  int num_items=0;

  // Recv the packed buffer.
  //MPI_Status status_mpi;
  //MPI_Recv(&bsize, sizeof(size_t), MPI_BYTE, 1, 1, MPI_COMM_WORLD, &status_mpi);
  //buffer = (char*) calloc( bsize, sizeof(char) );
  //MPI_Recv(buffer, int(bsize), MPI_BYTE, 1, 1, MPI_COMM_WORLD, &status_mpi);

  buff_ptr = buffer;

  // Now lets unpack the buffer into a new Blast_ScoreBlk
  BlastScoreBlk* my_sbp = (BlastScoreBlk*)calloc(1,sizeof(BlastScoreBlk));

  memcpy( &(my_sbp->protein_alphabet), buff_ptr, sizeof(Boolean) );
  buff_ptr += sizeof(Boolean);

  memcpy( &(my_sbp->alphabet_code), buff_ptr, sizeof(Uint1) );
  buff_ptr += sizeof(Uint1);

  memcpy( &(my_sbp->alphabet_size), buff_ptr, sizeof(Int2) );
  buff_ptr += sizeof(Int2);

  memcpy( &(my_sbp->alphabet_start), buff_ptr, sizeof(Int2) );
  buff_ptr += sizeof(Int2);

  memcpy( &(str_length), buff_ptr, sizeof(int) );
  buff_ptr += sizeof(int);

  // allocate the 'name' member.
  my_sbp->name = (char*)malloc( (str_length+1)*sizeof(char) );

  memcpy( my_sbp->name, buff_ptr, str_length*sizeof(char) );
  buff_ptr += (str_length * sizeof(char) );

  my_sbp->name[str_length] = '\0';

  // Get the number of items for the 'comments' member.
  memcpy( &num_items, buff_ptr, sizeof(int) );
  buff_ptr += sizeof(int);

  if ( num_items )
    {
      // Make the head of the linked list.
      ListNode* newnode = (ListNode*) calloc(1, sizeof(ListNode));

      memcpy( &(newnode->choice), buff_ptr, sizeof(Uint1) );
      buff_ptr += sizeof(Uint1);

      memcpy( &(str_length), buff_ptr, sizeof(int) );
      buff_ptr += sizeof(int);

      newnode->ptr = calloc( (str_length+1), sizeof(char) );

      // Make a string.
      char *my_string = (char*) malloc( (str_length+1)*sizeof(char) );

      memcpy( my_string, buff_ptr, (str_length*sizeof(char) ) );
      buff_ptr += ( str_length*sizeof(char) );

      my_string[str_length] = '\0';

      memcpy( newnode->ptr, my_string, ((str_length+1)*sizeof(char)) );

      newnode->next = NULL;

      free(my_string);

      // set the comments ptr to the newnode.
      my_sbp->comments = newnode;

      for ( int ncomms=1; ncomms < num_items; ++ncomms )
	{
	  // Make a new node.
	  newnode->next = (ListNode*) calloc(1, sizeof(ListNode));

	  // Advance our current spot in the list.
	  newnode = newnode->next;

	  memcpy( &(newnode->choice), buff_ptr, sizeof(Uint1) );
	  buff_ptr += sizeof(Uint1);

	  memcpy( &(str_length), buff_ptr, sizeof(int) );
	  buff_ptr += sizeof(int);

	  newnode->ptr = calloc( (str_length+1), sizeof(char) );

	  // Make a string.
	  my_string = (char*) malloc( (str_length+1)*sizeof(char) );

	  memcpy( my_string, buff_ptr, (str_length*sizeof(char) ) );
	  buff_ptr += ( str_length*sizeof(char) );

	  my_string[str_length] = '\0';

	  memcpy( newnode->ptr, my_string, ((str_length+1)*sizeof(char)) );

	  newnode->next = NULL;

	  free(my_string);

	}
    }
  else
    {
      my_sbp->comments = NULL;
    }

  // Did the BlastScoreBlk have a matrix.
  memcpy( &num_items, buff_ptr, sizeof(int) );
  buff_ptr += sizeof(int);

  if ( num_items )
    {
      size_t nr, nc;

      memcpy( &nc, buff_ptr, sizeof(size_t) );
      buff_ptr += sizeof(size_t);

      memcpy( &nr, buff_ptr, sizeof(size_t) );
      buff_ptr += sizeof(size_t);

      my_sbp->matrix = (SBlastScoreMatrix*) calloc( 1, sizeof(SBlastScoreMatrix) );

      my_sbp->matrix->data = (int**) _PSIAllocateMatrix( nc, nr, sizeof(int) );

      my_sbp->matrix->freqs = (double*) calloc( my_sbp->alphabet_size, sizeof(double) );

      for ( size_t this_row=0; this_row < nr; ++this_row )
	{
	  for ( size_t this_col=0; this_col < nc; ++this_col )
	    {
	      memcpy( &(my_sbp->matrix->data[this_row][this_col]), buff_ptr, sizeof(int) );
	      buff_ptr += sizeof(int);
	    }
	}

      for ( Int2 freqs_ind = 0; freqs_ind < my_sbp->alphabet_size ; ++freqs_ind )
	{
	  memcpy( &(my_sbp->matrix->freqs[freqs_ind]), buff_ptr, sizeof(double) );
	  buff_ptr += sizeof(double);
	}

      memcpy( &(my_sbp->matrix->lambda), buff_ptr, sizeof(double) );
      buff_ptr += sizeof(double);
    }
  else
    {
      my_sbp->matrix = NULL;
    }

  // Did the SPsiBlastScoreBlk have a matrix.
  memcpy( &num_items, buff_ptr, sizeof(int) );
  buff_ptr += sizeof(int);

  if ( num_items )
    {
      size_t nr, nc;
      memcpy( &nc, buff_ptr, sizeof(size_t) );
      buff_ptr += sizeof(size_t);

      memcpy( &nr, buff_ptr, sizeof(size_t) );
      buff_ptr += sizeof(size_t);

      my_sbp->psi_matrix = SPsiBlastScoreMatrixNew( nc );

      for ( size_t this_row=0; this_row < nr; ++this_row )
	{
	  for ( size_t this_col=0; this_col < nc; ++this_col )
	    {
	      memcpy( &(my_sbp->psi_matrix->pssm->data[this_row][this_col]), buff_ptr, sizeof(int) );
	      buff_ptr += sizeof(int);
	    }
	}

      for ( Int2 freqs_ind = 0; freqs_ind < my_sbp->alphabet_size ; ++freqs_ind )
	{
	  memcpy( &(my_sbp->psi_matrix->pssm->freqs[freqs_ind]), buff_ptr, sizeof(double) );
	  buff_ptr += sizeof(double);
	}

      memcpy( &(my_sbp->psi_matrix->pssm->lambda), buff_ptr, sizeof(double) );
      buff_ptr += sizeof(double);

      // now store the frequency ratios table
      for ( size_t this_row=0; this_row < nr; ++this_row )
	{
	  for ( size_t this_col=0; this_col < nc; ++this_col )
	    {
	      memcpy( &(my_sbp->psi_matrix->freq_ratios[this_row][this_col]), buff_ptr, sizeof(double) );
	      buff_ptr += sizeof(double);
	    }
	}

      // and lastly the associated Blast_KarlinBlk
      memcpy( &(my_sbp->psi_matrix->kbp->Lambda), buff_ptr, sizeof(double) );
      buff_ptr += sizeof(double);

      memcpy( &(my_sbp->psi_matrix->kbp->K), buff_ptr, sizeof(double) );
      buff_ptr += sizeof(double);

      memcpy( &(my_sbp->psi_matrix->kbp->logK), buff_ptr, sizeof(double) );
      buff_ptr += sizeof(double);

      memcpy( &(my_sbp->psi_matrix->kbp->H), buff_ptr, sizeof(double) );
      buff_ptr += sizeof(double);

      memcpy( &(my_sbp->psi_matrix->kbp->paramC), buff_ptr, sizeof(double) );
      buff_ptr += sizeof(double);

    }
  else
    {
      my_sbp->psi_matrix = NULL;
    }

  // more POD members.
  memcpy( &(my_sbp->matrix_only_scoring), buff_ptr, sizeof(Boolean) );
  buff_ptr += sizeof(Boolean);

  memcpy( &(my_sbp->complexity_adjusted_scoring), buff_ptr, sizeof(Boolean) );
  buff_ptr += sizeof(Boolean);

  memcpy( &(my_sbp->loscore), buff_ptr, sizeof(Int4) );
  buff_ptr += sizeof(Int4);

  memcpy( &(my_sbp->hiscore), buff_ptr, sizeof(Int4) );
  buff_ptr += sizeof(Int4);

  memcpy( &(my_sbp->penalty), buff_ptr, sizeof(Int4) );
  buff_ptr += sizeof(Int4);

  memcpy( &(my_sbp->reward), buff_ptr, sizeof(Int4) );
  buff_ptr += sizeof(Int4);

  memcpy( &(my_sbp->scale_factor), buff_ptr, sizeof(double) );
  buff_ptr += sizeof(double);

  memcpy( &(my_sbp->read_in_matrix), buff_ptr, sizeof(Boolean) );
  buff_ptr += sizeof(Boolean);

  // --- Need the 'number_of_contexts' variable to get the size of remaining arrays of ptr's

  memcpy( &(my_sbp->number_of_contexts), buff_ptr, sizeof(Int4) );
  buff_ptr += sizeof(Int4);

  memcpy( &num_items, buff_ptr, sizeof(int) );
  buff_ptr += sizeof(int);

  // Is there a Blast_ScoreFreq array
  if ( num_items )
    {
      my_sbp->sfp = (Blast_ScoreFreq**) calloc( my_sbp->number_of_contexts, sizeof(Blast_ScoreFreq*) );

      Blast_ScoreFreq* my_sfp = NULL;
      for ( Int4 ncont=0; ncont < my_sbp->number_of_contexts; ++ncont )
	{
	  my_sfp = (Blast_ScoreFreq*) calloc( 1, sizeof(Blast_ScoreFreq) );

	  memcpy( &(my_sfp->score_min), buff_ptr, sizeof(Int4) );
	  buff_ptr += sizeof(Int4);

	  memcpy( &(my_sfp->score_max), buff_ptr, sizeof(Int4) );
	  buff_ptr += sizeof(Int4);

	  memcpy( &(my_sfp->obs_min), buff_ptr, sizeof(Int4) );
	  buff_ptr += sizeof(Int4);

	  memcpy( &(my_sfp->obs_max), buff_ptr, sizeof(Int4) );
	  buff_ptr += sizeof(Int4);

	  memcpy( &(my_sfp->score_avg), buff_ptr, sizeof(double) );
	  buff_ptr += sizeof(double);

	  int score_range = my_sfp->score_max - my_sfp->score_min + 1;

	  my_sfp->sprob0 = (double*) calloc( score_range, sizeof(double) );

	  for ( int score_index = 0; score_index < score_range; ++score_index )
	    {
	      memcpy( &(my_sfp->sprob0[score_index]), buff_ptr, sizeof(double) );
	      buff_ptr += sizeof(double);
	    }

	  // This is to mirror the shift done in the NCBI code.
	  my_sfp->sprob = my_sfp->sprob0 - my_sfp->score_min ;

	  my_sbp->sfp[ncont] = my_sfp;
	  my_sfp = NULL;
	}
    }
  else
    {
      my_sbp->sfp = NULL;
    }

  // Karlin-Altschul parameters.
  memcpy( &num_items, buff_ptr, sizeof(int) );
  buff_ptr += sizeof(int);

  // Gumbel Block parameters.
  if ( num_items )
    {
      my_sbp->gbp = (Blast_GumbelBlk*) calloc( 1, sizeof(Blast_GumbelBlk) );

      memcpy( &(my_sbp->gbp->Lambda), buff_ptr, sizeof(double) );
      buff_ptr += sizeof(double);

      memcpy( &(my_sbp->gbp->C), buff_ptr, sizeof(double) );
      buff_ptr += sizeof(double);

      memcpy( &(my_sbp->gbp->G), buff_ptr, sizeof(double) );
      buff_ptr += sizeof(double);

      memcpy( &(my_sbp->gbp->a), buff_ptr, sizeof(double) );
      buff_ptr += sizeof(double);

      memcpy( &(my_sbp->gbp->Alpha), buff_ptr, sizeof(double) );
      buff_ptr += sizeof(double);

      memcpy( &(my_sbp->gbp->Sigma), buff_ptr, sizeof(double) );
      buff_ptr += sizeof(double);

      memcpy( &(my_sbp->gbp->a_un), buff_ptr, sizeof(double) );
      buff_ptr += sizeof(double);

      memcpy( &(my_sbp->gbp->Alpha_un), buff_ptr, sizeof(double) );
      buff_ptr += sizeof(double);

      memcpy( &(my_sbp->gbp->b), buff_ptr, sizeof(double) );
      buff_ptr += sizeof(double);

      memcpy( &(my_sbp->gbp->Beta), buff_ptr, sizeof(double) );
      buff_ptr += sizeof(double);

      memcpy( &(my_sbp->gbp->Tau), buff_ptr, sizeof(double) );
      buff_ptr += sizeof(double);

      memcpy( &(my_sbp->gbp->db_length), buff_ptr, sizeof(Int8) );
      buff_ptr += sizeof(Int8);

      memcpy( &(my_sbp->gbp->filled), buff_ptr, sizeof(Boolean) );
      buff_ptr += sizeof(Boolean);

    }
  else
    {
      my_sbp->gbp = NULL;
    }

  memcpy( &num_items, buff_ptr, sizeof(int) );
  buff_ptr += sizeof(int);

  // Karlin-Altschul for ungapped alignments.
  if ( num_items )
    {
      my_sbp->kbp_std = (Blast_KarlinBlk**) calloc(my_sbp->number_of_contexts, sizeof(Blast_KarlinBlk*));

      Blast_KarlinBlk* my_kbp = NULL;

      for ( Int4 ncont = 0; ncont < my_sbp->number_of_contexts; ++ncont )
	{
	  Blast_KarlinBlk* my_kbp = (Blast_KarlinBlk*) calloc( 1, sizeof(Blast_KarlinBlk) );

	  memcpy( &(my_kbp->Lambda), buff_ptr, sizeof(double) );
	  buff_ptr += sizeof(double);

	  memcpy( &(my_kbp->K), buff_ptr, sizeof(double) );
	  buff_ptr += sizeof(double);

	  memcpy( &(my_kbp->logK), buff_ptr, sizeof(double) );
	  buff_ptr += sizeof(double);

	  memcpy( &(my_kbp->H), buff_ptr, sizeof(double) );
	  buff_ptr += sizeof(double);

	  memcpy( &(my_kbp->paramC), buff_ptr, sizeof(double) );
	  buff_ptr += sizeof(double);

	  my_sbp->kbp_std[ncont] = my_kbp;

	  my_kbp = NULL;
	}
    }
  else
    {
      my_sbp->kbp_std = NULL;
    }

  memcpy( &num_items, buff_ptr, sizeof(int) );
  buff_ptr += sizeof(int);

  // Karlin-Altschul for position based alignments.
  if ( num_items )
    {
      my_sbp->kbp_psi = (Blast_KarlinBlk**) calloc(my_sbp->number_of_contexts, sizeof(Blast_KarlinBlk*));

      Blast_KarlinBlk* my_kbp = NULL;

      for ( Int4 ncont = 0; ncont < my_sbp->number_of_contexts; ++ncont )
	{
	  Blast_KarlinBlk* my_kbp = (Blast_KarlinBlk*) calloc( 1, sizeof(Blast_KarlinBlk) );

	  memcpy( &(my_kbp->Lambda), buff_ptr, sizeof(double) );
	  buff_ptr += sizeof(double);

	  memcpy( &(my_kbp->K), buff_ptr, sizeof(double) );
	  buff_ptr += sizeof(double);

	  memcpy( &(my_kbp->logK), buff_ptr, sizeof(double) );
	  buff_ptr += sizeof(double);

	  memcpy( &(my_kbp->H), buff_ptr, sizeof(double) );
	  buff_ptr += sizeof(double);

	  memcpy( &(my_kbp->paramC), buff_ptr, sizeof(double) );
	  buff_ptr += sizeof(double);

	  my_sbp->kbp_psi[ncont] = my_kbp;

	  my_kbp = NULL;
        }
    }
  else
    {
      my_sbp->kbp_psi = NULL;
    }

  memcpy( &num_items, buff_ptr, sizeof(int) );
  buff_ptr += sizeof(int);

  // Karlin-Altschul for std (not psi) alignments.
  if ( num_items )
    {
      my_sbp->kbp_gap_std = (Blast_KarlinBlk**) calloc(my_sbp->number_of_contexts, sizeof(Blast_KarlinBlk*));

      Blast_KarlinBlk* my_kbp = NULL;

      for ( Int4 ncont = 0; ncont < my_sbp->number_of_contexts; ++ncont )
	{
	  Blast_KarlinBlk* my_kbp = (Blast_KarlinBlk*) calloc( 1, sizeof(Blast_KarlinBlk) );

	  memcpy( &(my_kbp->Lambda), buff_ptr, sizeof(double) );
	  buff_ptr += sizeof(double);

	  memcpy( &(my_kbp->K), buff_ptr, sizeof(double) );
	  buff_ptr += sizeof(double);

	  memcpy( &(my_kbp->logK), buff_ptr, sizeof(double) );
	  buff_ptr += sizeof(double);

	  memcpy( &(my_kbp->H), buff_ptr, sizeof(double) );
	  buff_ptr += sizeof(double);

	  memcpy( &(my_kbp->paramC), buff_ptr, sizeof(double) );
	  buff_ptr += sizeof(double);

	  my_sbp->kbp_gap_std[ncont] = my_kbp;

	  my_kbp = NULL;
	}
    }
  else
    {
      my_sbp->kbp_gap_std = NULL;
    }

  memcpy( &num_items, buff_ptr, sizeof(int) );
  buff_ptr += sizeof(int);
  //std::cout<<"unpack KA gap psi num_items: "<<num_items <<std::endl;

  // Karlin-Altschul for psi alignments.
  if ( num_items )
    {
      my_sbp->kbp_gap_psi = (Blast_KarlinBlk**) calloc(my_sbp->number_of_contexts, sizeof(Blast_KarlinBlk*));
      //std::cout<<"unpack sbp->number_of_contexts: "<<my_sbp->number_of_contexts <<std::endl;

      Blast_KarlinBlk* my_kbp = NULL;

      for ( Int4 ncont = 0; ncont < my_sbp->number_of_contexts; ++ncont )
	{
	  Blast_KarlinBlk* my_kbp = (Blast_KarlinBlk*) calloc( 1, sizeof(Blast_KarlinBlk) );

	  memcpy( &(my_kbp->Lambda), buff_ptr, sizeof(double) );
	  buff_ptr += sizeof(double);

	  memcpy( &(my_kbp->K), buff_ptr, sizeof(double) );
	  buff_ptr += sizeof(double);

	  memcpy( &(my_kbp->logK), buff_ptr, sizeof(double) );
	  buff_ptr += sizeof(double);

	  memcpy( &(my_kbp->H), buff_ptr, sizeof(double) );
	  buff_ptr += sizeof(double);

	  memcpy( &(my_kbp->paramC), buff_ptr, sizeof(double) );
	  buff_ptr += sizeof(double);

	  my_sbp->kbp_gap_psi[ncont] = my_kbp;

	  my_kbp = NULL;
	  //std::cout<<"unpack sbp->kbp_gap_psi["<<ncont<<"]->Lambda: "<<my_sbp->kbp_gap_psi[ncont]->Lambda <<std::endl;
	  //std::cout<<"unpack sbp->kbp_gap_psi["<<ncont<<"]->paramC: "<<my_sbp->kbp_gap_psi[ncont]->paramC <<std::endl;

	}
    }
  else
    {
      my_sbp->kbp_gap_psi = NULL;
    }

  memcpy( &num_items, buff_ptr, sizeof(int) );
  buff_ptr += sizeof(int);

  // Karlin-Altschul with ideal values.
  if ( num_items )
    {
      my_sbp->kbp_ideal = (Blast_KarlinBlk*) calloc(my_sbp->number_of_contexts, sizeof(Blast_KarlinBlk));

      memcpy( &(my_sbp->kbp_ideal->Lambda), buff_ptr, sizeof(double) );
      buff_ptr += sizeof(double);

      memcpy( &(my_sbp->kbp_ideal->K), buff_ptr, sizeof(double) );
      buff_ptr += sizeof(double);

      memcpy( &(my_sbp->kbp_ideal->logK), buff_ptr, sizeof(double) );
      buff_ptr += sizeof(double);

      memcpy( &(my_sbp->kbp_ideal->H), buff_ptr, sizeof(double) );
      buff_ptr += sizeof(double);

      memcpy( &(my_sbp->kbp_ideal->paramC), buff_ptr, sizeof(double) );
      buff_ptr += sizeof(double);

    }
  else
    {
      my_sbp->kbp_ideal = NULL;
    }


  // Read the alias pointers.
  memcpy( &num_items, buff_ptr, sizeof(int) );
  buff_ptr += sizeof(int);

  if ( num_items == 0 )
    {
      my_sbp->kbp = NULL;
    }
  else if ( num_items == 1 )
    {
      my_sbp->kbp = my_sbp->kbp_std;
    }
  else if ( num_items == 2 )
    {
      my_sbp->kbp = my_sbp->kbp_std;
    }
  else
    {
      std::cerr<<"    GOT NUM ITEMS "<<num_items<<", so I don't know where to point sbp->kbp"<<std::endl;
    }

  memcpy( &num_items, buff_ptr, sizeof(int) );
  buff_ptr += sizeof(int);

  if ( num_items == 0 )
    {
      my_sbp->kbp_gap = NULL;
    }
  else if ( num_items == 1 )
    {
      my_sbp->kbp_gap = my_sbp->kbp_gap_std;
    }
  else if ( num_items == 2 )
    {
      my_sbp->kbp_gap = my_sbp->kbp_gap_std;
    }
  else
    {
      std::cerr<<"    GOT NUM ITEMS "<<num_items<<", so I don't know where to point sbp->kbp_gap"<<std::endl;
    }

  // Unpack the size of the ambiguous res array.
  memcpy( &(my_sbp->ambig_size), buff_ptr, sizeof(Int2) );
  buff_ptr += sizeof(Int2);

  // Unpack how much of it is used.
  memcpy( &(my_sbp->ambig_occupy), buff_ptr, sizeof(Int2) );
  buff_ptr += sizeof(Int2);

  my_sbp->ambiguous_res = (Uint1*)calloc( my_sbp->ambig_size, sizeof(Uint1) );

  // Grab the array sentinel value.
  memcpy( &num_items, buff_ptr, sizeof(int) );
  buff_ptr += sizeof(int);

  if ( num_items )
    {
      for ( int ambig_index=0; ambig_index < my_sbp->ambig_size; ambig_index++ )
	{
	  memcpy( &(my_sbp->ambiguous_res[ambig_index]), buff_ptr, sizeof(Uint1) );
	  buff_ptr += sizeof(Uint1);
	}
    }
  else
    {
      my_sbp->ambiguous_res = NULL;
    }

  memcpy( &(my_sbp->round_down), buff_ptr, sizeof(Boolean) );

  return my_sbp;
}




char* PackQueryInfoBuffer(BlastQueryInfo* QueryInfo,
                          size_t& QueryInfoBufferSize)
{
    // Dump the QueryInfo data structure.
    //BlastQueryInfo* bqi = m_InternalData->m_QueryInfo;
    BlastQueryInfo* bqi = QueryInfo;

    char* buffer = NULL;
    char* buff_ptr = NULL;
    size_t bsize = 0;

    // Determine size of buffer for data structure.
    bsize += sizeof(Int4);
    bsize += sizeof(Int4);
    bsize += sizeof(int);

    // Next is the BlastContextInfo pointer -- count up size of the structure and multiply be number of contexts.
    for ( Int4 ctx = bqi->first_context; ctx <= bqi->last_context; ++ctx )
      {
	bsize += sizeof(Int4);	
	bsize += sizeof(Int4);
	bsize += sizeof(Int8);
	bsize += sizeof(Int4);
	bsize += sizeof(Int4);
	bsize += sizeof(Int1);
	bsize += sizeof(Boolean);
      }
    bsize += sizeof(Uint4);


    buffer = (char*)calloc(bsize, sizeof(char));

    buff_ptr = buffer;

    // Serialize the data.
    memcpy( buff_ptr, &(bqi->first_context), sizeof(Int4));
    buff_ptr += sizeof(Int4);

    memcpy( buff_ptr, &(bqi->last_context), sizeof(Int4));
    buff_ptr += sizeof(Int4);

    memcpy( buff_ptr, &(bqi->num_queries), sizeof(int));
    buff_ptr += sizeof(int);

    //std::cout<<"pack bqi->last_context: "<<bqi->last_context<<std::endl;

    for ( Int4 ctx = bqi->first_context; ctx <= bqi->last_context; ++ctx )
      {
	memcpy( buff_ptr, &(bqi->contexts[ctx].query_offset), sizeof(Int4));
	buff_ptr += sizeof(Int4);

	memcpy( buff_ptr, &(bqi->contexts[ctx].query_length), sizeof(Int4));
	buff_ptr += sizeof(Int4);

	memcpy( buff_ptr, &(bqi->contexts[ctx].eff_searchsp), sizeof(Int8));
	buff_ptr += sizeof(Int8);

	memcpy( buff_ptr, &(bqi->contexts[ctx].length_adjustment), sizeof(Int4));
	buff_ptr += sizeof(Int4);

	memcpy( buff_ptr, &(bqi->contexts[ctx].query_index), sizeof(Int4));
	buff_ptr += sizeof(Int4);

	memcpy( buff_ptr, &(bqi->contexts[ctx].frame), sizeof(Int1));
	buff_ptr += sizeof(Int1);

	memcpy( buff_ptr, &(bqi->contexts[ctx].is_valid), sizeof(Boolean));
	buff_ptr += sizeof(Boolean);
      }
    
    memcpy( buff_ptr, &(bqi->max_length), sizeof(Uint4));
    //buff_ptr += sizeof(Uint4);

    // Save the QueryInfo
    QueryInfoBufferSize = bsize;
    return buffer;
}



BlastQueryInfo* UnpackQueryInfoBuffer(char * QueryInfoBuffer,
                                      size_t QueryInfoBufferSize)
{
    // Dump the QueryInfo data structure.
    BlastQueryInfo* bqi = new BlastQueryInfo();

    char* buffer = QueryInfoBuffer;
    char* buff_ptr = NULL;
    size_t bsize = QueryInfoBufferSize;

    buff_ptr = buffer;

    // Serialize the data.
    memcpy( &(bqi->first_context), buff_ptr, sizeof(Int4));
    buff_ptr += sizeof(Int4);

    memcpy( &(bqi->last_context), buff_ptr, sizeof(Int4));
    buff_ptr += sizeof(Int4);

    memcpy( &(bqi->num_queries), buff_ptr, sizeof(int));
    buff_ptr += sizeof(int);

    int num_contexts = bqi->last_context - bqi->first_context; 
    //std::cout<<"unpack num_contexts: "<<num_contexts<<std::endl;

    bqi->contexts = (BlastContextInfo*) calloc(num_contexts+1,sizeof(BlastContextInfo));

    //std::cout<<"unpack bqi->last_context: "<<bqi->last_context<<std::endl;

    for ( Int4 ctx = bqi->first_context; ctx <= bqi->last_context; ++ctx )
      {
	memcpy( &(bqi->contexts[ctx].query_offset), buff_ptr, sizeof(Int4));
	buff_ptr += sizeof(Int4);

	memcpy( &(bqi->contexts[ctx].query_length), buff_ptr, sizeof(Int4));
	buff_ptr += sizeof(Int4);

	memcpy( &(bqi->contexts[ctx].eff_searchsp), buff_ptr, sizeof(Int8));
	buff_ptr += sizeof(Int8);

	memcpy( &(bqi->contexts[ctx].length_adjustment), buff_ptr, sizeof(Int4));
	buff_ptr += sizeof(Int4);

	memcpy( &(bqi->contexts[ctx].query_index), buff_ptr, sizeof(Int4));
	buff_ptr += sizeof(Int4);

	memcpy( &(bqi->contexts[ctx].frame), buff_ptr, sizeof(Int1));
	buff_ptr += sizeof(Int1);

	memcpy( &(bqi->contexts[ctx].is_valid), buff_ptr, sizeof(Boolean));
	buff_ptr += sizeof(Boolean);
      }
    
    memcpy( &(bqi->max_length), buff_ptr, sizeof(Uint4));
    //buff_ptr += sizeof(Uint4);
    
    // Save the QueryInfo
    return bqi;
}

//========================================================





//Functions to output selected object to file ============================

void gvoifncDumpBlastHSPResults (BlastHSPResults *abhrpBlastHSPResults) {

  std::cout<<"num_queries "<<abhrpBlastHSPResults->num_queries<<std::endl;
  std::cout.flush();

  for (int i=0; i<(abhrpBlastHSPResults->num_queries); i++) {

    std::cout<<"  [BlastHitList "<<i<<"]"<<std::endl;
    std::cout.flush();

    std::cout<<"  hsplist_count "<<abhrpBlastHSPResults->hitlist_array[i]->hsplist_count<<std::endl;
    std::cout.flush();

    std::cout<<"  hsplist_max "<<abhrpBlastHSPResults->hitlist_array[i]->hsplist_max<<std::endl;
    std::cout.flush();

    std::cout<<"  worst_evalue "<<abhrpBlastHSPResults->hitlist_array[i]->worst_evalue<<std::endl;
    std::cout.flush();

    std::cout<<"  low_score "<<abhrpBlastHSPResults->hitlist_array[i]->low_score<<std::endl;
    std::cout.flush();

    std::cout<<"  heapified "<<abhrpBlastHSPResults->hitlist_array[i]->heapified<<std::endl;
    std::cout.flush();

    for (int j=0; j<(abhrpBlastHSPResults->hitlist_array[i]->hsplist_count); j++) {

      std::cout<<"    [BlastHSPList "<<j<<"]"<<std::endl;
      std::cout.flush();

      std::cout<<"    oid "<<abhrpBlastHSPResults->hitlist_array[i]->hsplist_array[j]->oid<<std::endl;
      std::cout.flush();

      std::cout<<"    query_index "<<abhrpBlastHSPResults->hitlist_array[i]->hsplist_array[j]->query_index<<std::endl;
      std::cout.flush();

      for (int k=0; k<(abhrpBlastHSPResults->hitlist_array[i]->hsplist_array[j]->hspcnt); k++) {

        std::cout<<"      [BlastHSP "<<k<<"]"<<std::endl;
        std::cout.flush();

        std::cout<<"      score "<<abhrpBlastHSPResults->hitlist_array[i]->hsplist_array[j]->hsp_array[k]->score<<std::endl;
        std::cout.flush();

        std::cout<<"      num_ident "<<abhrpBlastHSPResults->hitlist_array[i]->hsplist_array[j]->hsp_array[k]->num_ident<<std::endl;
        std::cout.flush();

        std::cout<<"      bit_score "<<abhrpBlastHSPResults->hitlist_array[i]->hsplist_array[j]->hsp_array[k]->bit_score<<std::endl;
        std::cout.flush();

        std::cout<<"      evalue "<<abhrpBlastHSPResults->hitlist_array[i]->hsplist_array[j]->hsp_array[k]->evalue<<std::endl;
        std::cout.flush();

        std::cout<<"        query frame "<<abhrpBlastHSPResults->hitlist_array[i]->hsplist_array[j]->hsp_array[k]->query.frame<<std::endl;
        std::cout.flush();

        std::cout<<"        query offset "<<abhrpBlastHSPResults->hitlist_array[i]->hsplist_array[j]->hsp_array[k]->query.offset<<std::endl;
        std::cout.flush();

        std::cout<<"        query end "<<abhrpBlastHSPResults->hitlist_array[i]->hsplist_array[j]->hsp_array[k]->query.end<<std::endl;
        std::cout.flush();

        std::cout<<"        query gapped_start "<<abhrpBlastHSPResults->hitlist_array[i]->hsplist_array[j]->hsp_array[k]->query.gapped_start<<std::endl;
        std::cout.flush();

        std::cout<<"        subject frame "<<abhrpBlastHSPResults->hitlist_array[i]->hsplist_array[j]->hsp_array[k]->subject.frame<<std::endl;
        std::cout.flush();

        std::cout<<"        subject offset "<<abhrpBlastHSPResults->hitlist_array[i]->hsplist_array[j]->hsp_array[k]->subject.offset<<std::endl;
        std::cout.flush();

        std::cout<<"        subject end "<<abhrpBlastHSPResults->hitlist_array[i]->hsplist_array[j]->hsp_array[k]->subject.end<<std::endl;
        std::cout.flush();

        std::cout<<"        gapped_start "<<abhrpBlastHSPResults->hitlist_array[i]->hsplist_array[j]->hsp_array[k]->subject.gapped_start<<std::endl;
        std::cout.flush();

        std::cout<<"      context "<<abhrpBlastHSPResults->hitlist_array[i]->hsplist_array[j]->hsp_array[k]->context<<std::endl;
        std::cout.flush();

        for (int l=0; l<(abhrpBlastHSPResults->hitlist_array[i]->hsplist_array[j]->hsp_array[k]->gap_info->size); l++) {

          std::cout<<"        op_type "<<abhrpBlastHSPResults->hitlist_array[i]->hsplist_array[j]->hsp_array[k]->gap_info->op_type[l]<<std::endl;
          std::cout.flush();

        }

        for (int l=0; l<(abhrpBlastHSPResults->hitlist_array[i]->hsplist_array[j]->hsp_array[k]->gap_info->size); l++) {

          std::cout<<"        num "<<abhrpBlastHSPResults->hitlist_array[i]->hsplist_array[j]->hsp_array[k]->gap_info->num[l]<<std::endl;
          std::cout.flush();

        }

        std::cout<<"      size "<<abhrpBlastHSPResults->hitlist_array[i]->hsplist_array[j]->hsp_array[k]->gap_info->size<<std::endl;
        std::cout.flush();

        std::cout<<"      num "<<abhrpBlastHSPResults->hitlist_array[i]->hsplist_array[j]->hsp_array[k]->num<<std::endl;
        std::cout.flush();

        std::cout<<"      comp_adjustment_method "<<abhrpBlastHSPResults->hitlist_array[i]->hsplist_array[j]->hsp_array[k]->comp_adjustment_method<<std::endl;
        std::cout.flush();

        if (abhrpBlastHSPResults->hitlist_array[i]->hsplist_array[j]->hsp_array[k]->pat_info) {

          std::cout<<"*** pat_info ***"<<std::endl;
          std::cout.flush();

        }

        std::cout<<"      num_positives "<<abhrpBlastHSPResults->hitlist_array[i]->hsplist_array[j]->hsp_array[k]->num_positives<<std::endl;
        std::cout.flush();

      }

      std::cout<<"    hspcnt "<<abhrpBlastHSPResults->hitlist_array[i]->hsplist_array[j]->hspcnt<<std::endl;
      std::cout.flush();

      std::cout<<"    allocated "<<abhrpBlastHSPResults->hitlist_array[i]->hsplist_array[j]->allocated<<std::endl;
      std::cout.flush();

      std::cout<<"    hsp_max "<<abhrpBlastHSPResults->hitlist_array[i]->hsplist_array[j]->hsp_max<<std::endl;
      std::cout.flush();

      std::cout<<"    do_not_reallocate "<<abhrpBlastHSPResults->hitlist_array[i]->hsplist_array[j]->do_not_reallocate<<std::endl;
      std::cout.flush();

      std::cout<<"    best_evalue "<<abhrpBlastHSPResults->hitlist_array[i]->hsplist_array[j]->best_evalue<<std::endl;
      std::cout.flush();

    }

    std::cout<<"  hsplist_current "<<abhrpBlastHSPResults->hitlist_array[i]->hsplist_current<<std::endl;
    std::cout.flush();

  }

}



void gvoifncDumpBlastQueryInfo (BlastQueryInfo *abqipBlastQueryInfo) {

  std::cout<<"first_context "<<abqipBlastQueryInfo->first_context<<std::endl;
  std::cout.flush();

  std::cout<<"last_context "<<abqipBlastQueryInfo->last_context<<std::endl;
  std::cout.flush();

  std::cout<<"num_queries "<<abqipBlastQueryInfo->num_queries<<std::endl;
  std::cout.flush();

  for (int i=abqipBlastQueryInfo->first_context; i<=(abqipBlastQueryInfo->last_context); i++) {

    std::cout<<"query_offset "<<abqipBlastQueryInfo->contexts[i].query_offset<<std::endl;
    std::cout.flush();

    std::cout<<"query_length "<<abqipBlastQueryInfo->contexts[i].query_length<<std::endl;
    std::cout.flush();

    std::cout<<"eff_searchsp "<<abqipBlastQueryInfo->contexts[i].eff_searchsp<<std::endl;
    std::cout.flush();

    std::cout<<"length_adjustment "<<abqipBlastQueryInfo->contexts[i].length_adjustment<<std::endl;
    std::cout.flush();

    std::cout<<"query_index "<<abqipBlastQueryInfo->contexts[i].query_index<<std::endl;
    std::cout.flush();

    std::cout<<"frame "<<abqipBlastQueryInfo->contexts[i].frame<<std::endl;
    std::cout.flush();

    std::cout<<"is_valid "<<abqipBlastQueryInfo->contexts[i].is_valid<<std::endl;
    std::cout.flush();

  }

  std::cout<<"max_length "<<abqipBlastQueryInfo->max_length<<std::endl;
  std::cout.flush();

  if (abqipBlastQueryInfo->pattern_info) {
    std::cout<<"*** pattern_info *** "<<std::endl;
    std::cout.flush();
  }

}


void gvoifncDumpBlast_KarlinBlk (Blast_KarlinBlk **abkbppBlast_KarlinBlk, int aintNumberOfContexts) {

  for (int i=0; i<aintNumberOfContexts; i++) {

    std::cout<<"  Blast_KarlinBlk] "<<i<<" "<<std::endl;
    std::cout.flush();

    std::cout<<"    Lambda "<<abkbppBlast_KarlinBlk[i]->Lambda<<std::endl;
    std::cout.flush();

    std::cout<<"    K "<<abkbppBlast_KarlinBlk[i]->K<<std::endl;
    std::cout.flush();

    std::cout<<"    logK "<<abkbppBlast_KarlinBlk[i]->logK<<std::endl;
    std::cout.flush();

    std::cout<<"    H "<<abkbppBlast_KarlinBlk[i]->H<<std::endl;
    std::cout.flush();

    std::cout<<"    paramC "<<abkbppBlast_KarlinBlk[i]->paramC<<std::endl;
    std::cout.flush();

  }

}

void gvoifncDumpBlast_KarlinBlk1 (Blast_KarlinBlk *abkbpBlast_KarlinBlk) {

  std::cout<<"  Lambda "<<abkbpBlast_KarlinBlk->Lambda<<std::endl;
  std::cout.flush();

  std::cout<<"  K "<<abkbpBlast_KarlinBlk->K<<std::endl;
  std::cout.flush();

  std::cout<<"  logK "<<abkbpBlast_KarlinBlk->logK<<std::endl;
  std::cout.flush();

  std::cout<<"  H "<<abkbpBlast_KarlinBlk->H<<std::endl;
  std::cout.flush();

  std::cout<<"  paramC "<<abkbpBlast_KarlinBlk->paramC<<std::endl;
  std::cout.flush();

}


void gvoifncDumpSBlastScoreMatrix(SBlastScoreMatrix *absmpSBlastScoreMatrix, int alphabet_size) {

  std::cout<<" read nc: "<<absmpSBlastScoreMatrix->ncols<<std::endl;
  std::cout.flush();

  std::cout<<" read nr: "<<absmpSBlastScoreMatrix->nrows<<std::endl;
  std::cout.flush();

  for (size_t i=0; i<absmpSBlastScoreMatrix->nrows; ++i) {
    for (size_t j=0; j<absmpSBlastScoreMatrix->ncols; ++j) {
      std::cout<<"    unpacked absmpSBlastScoreMatrix "<<i<<","<<j<<","<<absmpSBlastScoreMatrix->data[i][j]<<std::endl;
      std::cout.flush();
    }
  }
	    
  for (Int2 i=0; i<alphabet_size ; ++i) {
    std::cout<<"    unpacked freqs "<<i<<" "<<absmpSBlastScoreMatrix->data[i]<<std::endl;
    std::cout.flush();
  }

  std::cout<<"  unpacked lambda : "<<absmpSBlastScoreMatrix->lambda<<std::endl;
  std::cout.flush();

} 


void gvoifncDumpBlastScoreBlk(BlastScoreBlk *absbpBlastScoreBlk) {

  std::cout<<"unpacked protein alphabet"<<std::endl;
  std::cout<<"read : "<<(int)absbpBlastScoreBlk->protein_alphabet<<endl;
  std::cout.flush();

  std::cout<<"unpacked alphabet code"<<std::endl;
  std::cout<<"read : "<<(int)absbpBlastScoreBlk->alphabet_code<<endl;
  std::cout.flush();

  std::cout<<"unpacked alphabet size"<<std::endl;
  std::cout<<"read : "<<absbpBlastScoreBlk->alphabet_size<<endl;
  std::cout.flush();

  std::cout<<"unpacked alphabet start"<<std::endl;
  std::cout<<"read : "<<absbpBlastScoreBlk->alphabet_start<<endl;
  std::cout.flush();

  std::cout<<"unpacked name"<<std::endl;
  std::cout<<"read : "<<absbpBlastScoreBlk->name<<std::endl;
  std::cout.flush();

  ListNode *llndpNode=absbpBlastScoreBlk->comments;

  while (llndpNode) {

    std::cout<<"  choice"<<std::endl;
    std::cout<<"  read : "<<llndpNode->choice<<std::endl;
    std::cout.flush();

    std::cout<<"  ptr"<<std::endl;
    std::cout<<"  read : "<<llndpNode->ptr<<std::endl;
    std::cout.flush();

  }



  if (absbpBlastScoreBlk->matrix) {
    std::cout<<"  matrix "<<std::endl;
    std::cout.flush();

    gvoifncDumpSBlastScoreMatrix(absbpBlastScoreBlk->matrix,absbpBlastScoreBlk->alphabet_size);
  }

  if (absbpBlastScoreBlk->psi_matrix) {
    std::cout<<"  psi_matrix, pssm "<<std::endl;
    std::cout.flush();

    gvoifncDumpSBlastScoreMatrix(absbpBlastScoreBlk->psi_matrix->pssm,absbpBlastScoreBlk->alphabet_size);

    std::cout<<"  psi_matrix, freq_ratios "<<std::endl;
    std::cout.flush();

    for (size_t i=0; i<absbpBlastScoreBlk->psi_matrix->pssm->nrows; ++i) {
      for (size_t j=0; j<absbpBlastScoreBlk->psi_matrix->pssm->ncols; ++j) {
        std::cout<<"    freq_ratio["<<i<<"],["<<j<<"], "<<absbpBlastScoreBlk->psi_matrix->freq_ratios[i][j]<<std::endl;
        std::cout.flush();
      }
    }
	    
    std::cout<<"    Lambda "<<absbpBlastScoreBlk->psi_matrix->kbp->Lambda<<std::endl;
    std::cout.flush();

    std::cout<<"    K "<<absbpBlastScoreBlk->psi_matrix->kbp->K<<std::endl;
    std::cout.flush();

    std::cout<<"    logK "<<absbpBlastScoreBlk->psi_matrix->kbp->logK<<std::endl;
    std::cout.flush();

    std::cout<<"    H "<<absbpBlastScoreBlk->psi_matrix->kbp->H<<std::endl;
    std::cout.flush();

    std::cout<<"    paramC "<<absbpBlastScoreBlk->psi_matrix->kbp->paramC<<std::endl;
    std::cout.flush();

  }



  std::cout<<"matrix_only_scoring "<<absbpBlastScoreBlk->matrix_only_scoring<<std::endl;
  std::cout.flush();

  std::cout<<"complexity_adjusted_scoring "<<absbpBlastScoreBlk->complexity_adjusted_scoring<<std::endl;
  std::cout.flush();

  std::cout<<"loscore "<<absbpBlastScoreBlk->loscore<<std::endl;
  std::cout.flush();

  std::cout<<"hiscore "<<absbpBlastScoreBlk->hiscore<<std::endl;
  std::cout.flush();

  std::cout<<"penalty "<<absbpBlastScoreBlk->penalty<<std::endl;
  std::cout.flush();

  std::cout<<"reward "<<absbpBlastScoreBlk->reward<<std::endl;
  std::cout.flush();

  std::cout<<"scale_factor "<<absbpBlastScoreBlk->scale_factor<<std::endl;
  std::cout.flush();

  std::cout<<"read_in_matrix "<<absbpBlastScoreBlk->read_in_matrix<<std::endl;
  std::cout.flush();



  if (absbpBlastScoreBlk->sfp) {

    for (Int4 i=0; i<absbpBlastScoreBlk->number_of_contexts; ++i) {
    
      std::cout<<"[sfp] "<<i<<std::endl;
      std::cout.flush();

      std::cout<<"  score_min "<<absbpBlastScoreBlk->sfp[i]->score_min<<std::endl;
      std::cout.flush();

      std::cout<<"  score_max "<<absbpBlastScoreBlk->sfp[i]->score_max<<std::endl;
      std::cout.flush();

      std::cout<<"  obs_min "<<absbpBlastScoreBlk->sfp[i]->obs_min<<std::endl;
      std::cout.flush();

      std::cout<<"  obs_max "<<absbpBlastScoreBlk->sfp[i]->obs_max<<std::endl;
      std::cout.flush();

      std::cout<<"  score_avg "<<absbpBlastScoreBlk->sfp[i]->score_avg<<std::endl;
      std::cout.flush();

      for (int j=0; j< absbpBlastScoreBlk->sfp[i]->score_max; ++j ) {

        std::cout<<"    sprob0["<<j<<"]"<<absbpBlastScoreBlk->sfp[i]->sprob0[j]<<std::endl;
        std::cout.flush();

      }

      for ( int j = 0; j < absbpBlastScoreBlk->sfp[i]->score_max; ++j ) {

        std::cout<<"    sprob["<<j<<"]"<<absbpBlastScoreBlk->sfp[i]->sprob[j]<<std::endl;
        std::cout.flush();

      }

    }

  }




  if (absbpBlastScoreBlk->kbp) {

    std::cout<<"  kbp "<<std::endl;
    std::cout.flush();

    gvoifncDumpBlast_KarlinBlk (absbpBlastScoreBlk->kbp, absbpBlastScoreBlk->number_of_contexts);
 
  }

  if (absbpBlastScoreBlk->kbp_gap) {

    std::cout<<"  kbp_gap "<<std::endl;
    std::cout.flush();

    std::cout.flush();

    gvoifncDumpBlast_KarlinBlk (absbpBlastScoreBlk->kbp_gap, absbpBlastScoreBlk->number_of_contexts);
 
  }

  if (absbpBlastScoreBlk->gbp) {

    std::cout<<"  gbp "<<std::endl;
    std::cout.flush();

    std::cout<<"  Lambda "<<absbpBlastScoreBlk->gbp->Lambda<<std::endl;
    std::cout.flush();

    std::cout<<"  C "<<absbpBlastScoreBlk->gbp->C<<std::endl;
    std::cout.flush();

    std::cout<<"  G "<<absbpBlastScoreBlk->gbp->G<<std::endl;
    std::cout.flush();

    std::cout<<"  a "<<absbpBlastScoreBlk->gbp->a<<std::endl;
    std::cout.flush();

    std::cout<<"  Alpha "<<absbpBlastScoreBlk->gbp->Alpha<<std::endl;
    std::cout.flush();

    std::cout<<"  Sigma "<<absbpBlastScoreBlk->gbp->Sigma<<std::endl;
    std::cout.flush();

    std::cout<<"  a_un "<<absbpBlastScoreBlk->gbp->a_un<<std::endl;
    std::cout.flush();

    std::cout<<"  Alpha_un "<<absbpBlastScoreBlk->gbp->Alpha_un<<std::endl;
    std::cout.flush();

    std::cout<<"  b "<<absbpBlastScoreBlk->gbp->b<<std::endl;
    std::cout.flush();

    std::cout<<"  Beta "<<absbpBlastScoreBlk->gbp->Beta<<std::endl;
    std::cout.flush();

    std::cout<<"  Tau "<<absbpBlastScoreBlk->gbp->Tau<<std::endl;
    std::cout.flush();

    std::cout<<"  db_length "<<absbpBlastScoreBlk->gbp->db_length<<std::endl;
    std::cout.flush();

    std::cout<<"  filled "<<absbpBlastScoreBlk->gbp->filled<<std::endl;
    std::cout.flush();

  }

  if (absbpBlastScoreBlk->kbp_std) {

    std::cout<<"  kbp_std "<<std::endl;
    std::cout.flush();

    gvoifncDumpBlast_KarlinBlk (absbpBlastScoreBlk->kbp_std, absbpBlastScoreBlk->number_of_contexts);
 
  }

  if (absbpBlastScoreBlk->kbp_psi) {

    std::cout<<"  kbp_psi "<<std::endl;
    std::cout.flush();

    gvoifncDumpBlast_KarlinBlk (absbpBlastScoreBlk->kbp_psi, absbpBlastScoreBlk->number_of_contexts);
 
  }

  if (absbpBlastScoreBlk->kbp_gap_std) {

    std::cout<<"  kbp_gap_std "<<std::endl;
    std::cout.flush();

    gvoifncDumpBlast_KarlinBlk (absbpBlastScoreBlk->kbp_gap_std, absbpBlastScoreBlk->number_of_contexts);
 
  }

  if (absbpBlastScoreBlk->kbp_gap_psi) {

    std::cout<<"  kbp_gap_pst "<<std::endl;
    std::cout.flush();

    gvoifncDumpBlast_KarlinBlk (absbpBlastScoreBlk->kbp_gap_psi, absbpBlastScoreBlk->number_of_contexts);
 
  }

  if (absbpBlastScoreBlk->kbp_ideal) {

    std::cout<<"  kbp_ideal "<<std::endl;
    std::cout.flush();

    gvoifncDumpBlast_KarlinBlk1 (absbpBlastScoreBlk->kbp_ideal);
 
  }

  std::cout<<"number_of_contexts "<<absbpBlastScoreBlk->number_of_contexts<<std::endl;
  std::cout.flush();

  for ( int i=0; i<absbpBlastScoreBlk->ambig_size; ++i) {

    std::cout<<"  i "<<absbpBlastScoreBlk->ambiguous_res[i]<<std::endl;
    std::cout.flush();

  }

  std::cout<<"ambig_size "<<absbpBlastScoreBlk->ambig_size<<std::endl;
  std::cout.flush();

  std::cout<<"ambig_occupy "<<absbpBlastScoreBlk->ambig_occupy<<std::endl;
  std::cout.flush();

  std::cout<<"round_down "<<absbpBlastScoreBlk->round_down<<std::endl;
  std::cout.flush();

}



// Functions to combine objects across thread groups

BlastScoreBlk *gbsbpfncCombineThreadGroupBlastScoreBlk(int aintNumberOfThreadGroups, BlastScoreBlk **ScoreBlkGroup) {

  BlastScoreBlk* lbsbpBlastScoreBlk = (BlastScoreBlk*)calloc(1,sizeof(BlastScoreBlk));

  lbsbpBlastScoreBlk->protein_alphabet=ScoreBlkGroup[0]->protein_alphabet;

  lbsbpBlastScoreBlk->alphabet_code=ScoreBlkGroup[0]->alphabet_code;

  lbsbpBlastScoreBlk->alphabet_size=ScoreBlkGroup[0]->alphabet_size;

  lbsbpBlastScoreBlk->alphabet_start=ScoreBlkGroup[0]->alphabet_start;

  lbsbpBlastScoreBlk->name=ScoreBlkGroup[0]->name;

  lbsbpBlastScoreBlk->comments=NULL;

  ListNode *llndpCurrentNode=NULL;

  ListNode *llndpPreviousNode=NULL;


  for (int i=0; i<aintNumberOfThreadGroups; i++) {

    if (ScoreBlkGroup[i]->comments) {

      if (lbsbpBlastScoreBlk->comments) {

        llndpPreviousNode=lbsbpBlastScoreBlk->comments;

        llndpCurrentNode=lbsbpBlastScoreBlk->comments->next;

         while (llndpCurrentNode) {

           llndpPreviousNode=llndpCurrentNode;

           llndpCurrentNode=llndpCurrentNode->next;

         }

         llndpPreviousNode->next=ScoreBlkGroup[i]->comments;

       } else {

         lbsbpBlastScoreBlk->comments=ScoreBlkGroup[i]->comments;

      }

    }

  }


  lbsbpBlastScoreBlk->matrix = NULL;

  for (int i=0; i<aintNumberOfThreadGroups; i++) {

    if (ScoreBlkGroup[i]->matrix) {

      lbsbpBlastScoreBlk->matrix=ScoreBlkGroup[i]->matrix;

    }

  }


  lbsbpBlastScoreBlk->psi_matrix = NULL;

  for (int i=0; i<aintNumberOfThreadGroups; i++) {
    if (ScoreBlkGroup[i]->psi_matrix) {
      lbsbpBlastScoreBlk->psi_matrix=ScoreBlkGroup[i]->psi_matrix;
    }
  }


  lbsbpBlastScoreBlk->matrix_only_scoring=ScoreBlkGroup[0]->matrix_only_scoring;
  lbsbpBlastScoreBlk->complexity_adjusted_scoring=ScoreBlkGroup[0]->complexity_adjusted_scoring;
  lbsbpBlastScoreBlk->loscore=ScoreBlkGroup[0]->loscore;
  lbsbpBlastScoreBlk->hiscore=ScoreBlkGroup[0]->hiscore;
  lbsbpBlastScoreBlk->penalty=ScoreBlkGroup[0]->penalty;
  lbsbpBlastScoreBlk->reward=ScoreBlkGroup[0]->reward;
  lbsbpBlastScoreBlk->scale_factor=ScoreBlkGroup[0]->scale_factor;
  lbsbpBlastScoreBlk->read_in_matrix=ScoreBlkGroup[0]->read_in_matrix;
  


  int lintTotalNumberOfContexts=0;

  for (int i=0; i<aintNumberOfThreadGroups; i++) {
    lintTotalNumberOfContexts+=ScoreBlkGroup[i]->number_of_contexts;
    //lintTotalNumberOfContexts+=ScoreBlkGroup[0]->number_of_contexts;
  }

  lbsbpBlastScoreBlk->number_of_contexts=lintTotalNumberOfContexts;

  lbsbpBlastScoreBlk->sfp=NULL;

  if (lbsbpBlastScoreBlk->number_of_contexts>0) { 

    if (ScoreBlkGroup[0]->sfp) {

      lbsbpBlastScoreBlk->sfp = (Blast_ScoreFreq**) calloc( lbsbpBlastScoreBlk->number_of_contexts, sizeof(Blast_ScoreFreq*) );

      int lintGlobalContextNumber=-1;

      for (int i=0; i<aintNumberOfThreadGroups; i++) {

        for (Int4 j=0; j < ScoreBlkGroup[i]->number_of_contexts; ++j) {

          lintGlobalContextNumber++;

          lbsbpBlastScoreBlk->sfp[lintGlobalContextNumber]=ScoreBlkGroup[i]->sfp[j];

        }

      }

    }

  }



  lbsbpBlastScoreBlk->kbp=NULL;

  if (lintTotalNumberOfContexts>0) {
    if (ScoreBlkGroup[0]->kbp) {
      lbsbpBlastScoreBlk->kbp = (Blast_KarlinBlk**) calloc( lbsbpBlastScoreBlk->number_of_contexts, sizeof(Blast_KarlinBlk*) );
      int lintGlobalContextNumber=-1;

      for (int i=0; i<aintNumberOfThreadGroups; i++) {
        for (Int4 j=0; j < ScoreBlkGroup[i]->number_of_contexts; ++j) {
          lintGlobalContextNumber++;
          lbsbpBlastScoreBlk->kbp[lintGlobalContextNumber]=ScoreBlkGroup[i]->kbp[j];
        }
      }
    }
  }



  lbsbpBlastScoreBlk->kbp_gap=NULL;

  if (lintTotalNumberOfContexts>0) {
    if (ScoreBlkGroup[0]->kbp_gap) {
      lbsbpBlastScoreBlk->kbp_gap = (Blast_KarlinBlk**) calloc( lbsbpBlastScoreBlk->number_of_contexts, sizeof(Blast_KarlinBlk*) );
      int lintGlobalContextNumber=-1;

      for (int i=0; i<aintNumberOfThreadGroups; i++) {
        for (Int4 j=0; j < ScoreBlkGroup[i]->number_of_contexts; ++j) {
          lintGlobalContextNumber++;
          lbsbpBlastScoreBlk->kbp_gap[lintGlobalContextNumber]=ScoreBlkGroup[i]->kbp_gap[j];
        }
      }
    }
  }


  lbsbpBlastScoreBlk->gbp = NULL;

  for (int i=0; i<aintNumberOfThreadGroups; i++) {
    if (ScoreBlkGroup[i]->gbp) {
      lbsbpBlastScoreBlk->gbp=ScoreBlkGroup[i]->gbp;
    }
  }



  lbsbpBlastScoreBlk->kbp_std=NULL;

  if (lintTotalNumberOfContexts>0) {
    if (ScoreBlkGroup[0]->kbp_std) {
      lbsbpBlastScoreBlk->kbp_std = (Blast_KarlinBlk**) calloc( lbsbpBlastScoreBlk->number_of_contexts, sizeof(Blast_KarlinBlk*) );
      int lintGlobalContextNumber=-1;

      for (int i=0; i<aintNumberOfThreadGroups; i++) {
        for (Int4 j=0; j < ScoreBlkGroup[i]->number_of_contexts; ++j) {
          lintGlobalContextNumber++;
          lbsbpBlastScoreBlk->kbp_std[lintGlobalContextNumber]=ScoreBlkGroup[i]->kbp_std[j];
        }
      }
    }
  }



  lbsbpBlastScoreBlk->kbp_psi=NULL;

  if (lintTotalNumberOfContexts>0) {
    if (ScoreBlkGroup[0]->kbp_psi) {
      lbsbpBlastScoreBlk->kbp_psi = (Blast_KarlinBlk**) calloc( lbsbpBlastScoreBlk->number_of_contexts, sizeof(Blast_KarlinBlk*) );
      int lintGlobalContextNumber=-1;

      for (int i=0; i<aintNumberOfThreadGroups; i++) {
        for (Int4 j=0; j < ScoreBlkGroup[i]->number_of_contexts; ++j ) {
          lintGlobalContextNumber++;
          lbsbpBlastScoreBlk->kbp_psi[lintGlobalContextNumber] = ScoreBlkGroup[i]->kbp_psi[j];
        }
      }
    }
  }


  lbsbpBlastScoreBlk->kbp_gap_std = NULL;

  if (lintTotalNumberOfContexts>0) {
    if (ScoreBlkGroup[0]->kbp_gap_std) {
      lbsbpBlastScoreBlk->kbp_gap_std = (Blast_KarlinBlk**) calloc( lbsbpBlastScoreBlk->number_of_contexts, sizeof(Blast_KarlinBlk*) );
      int lintGlobalContextNumber = -1;

      for (int i=0; i<aintNumberOfThreadGroups; i++) 
      {
        for (Int4 j=0; j < ScoreBlkGroup[i]->number_of_contexts; ++j ) 
        {
          lintGlobalContextNumber++;
          lbsbpBlastScoreBlk->kbp_gap_std[lintGlobalContextNumber] = ScoreBlkGroup[i]->kbp_gap_std[j];
        }
      }
    }
  }


  lbsbpBlastScoreBlk->kbp_gap_psi=NULL;

  if (lintTotalNumberOfContexts>0) 
  {
    if (ScoreBlkGroup[0]->kbp_gap_psi) 
    {
      lbsbpBlastScoreBlk->kbp_gap_psi = (Blast_KarlinBlk**) calloc( lbsbpBlastScoreBlk->number_of_contexts, sizeof(Blast_KarlinBlk*) );
      int lintGlobalContextNumber = -1;

      for (int i=0; i < aintNumberOfThreadGroups; i++) 
      {
        for (Int4 j=0; j < ScoreBlkGroup[i]->number_of_contexts; ++j ) 
        {
          lintGlobalContextNumber++;
          lbsbpBlastScoreBlk->kbp_gap_psi[lintGlobalContextNumber] = ScoreBlkGroup[i]->kbp_gap_psi[j];
        }
      }
    }
  }


  lbsbpBlastScoreBlk->kbp_ideal = NULL;

  for (int i = 0; i < aintNumberOfThreadGroups; i++) 
  {
    if (ScoreBlkGroup[i]->kbp_ideal) 
    {
      lbsbpBlastScoreBlk->kbp_ideal = ScoreBlkGroup[i]->kbp_ideal;
    }
  }


  lbsbpBlastScoreBlk->ambig_size = ScoreBlkGroup[0]->ambig_size;
  lbsbpBlastScoreBlk->ambig_occupy = ScoreBlkGroup[0]->ambig_occupy;
  lbsbpBlastScoreBlk->ambiguous_res = ScoreBlkGroup[0]->ambiguous_res;
  lbsbpBlastScoreBlk->round_down = ScoreBlkGroup[0]->round_down;

  return lbsbpBlastScoreBlk;
}











BlastScoreBlk *gbsbpfncDuplicateBlastScoreBlk(BlastScoreBlk *ScoreBlkGroup)
{

#if 0
  for (index=0; index<sbp->number_of_contexts; index++) {
    if (sbp->sfp)
      sbp->sfp[index] = Blast_ScoreFreqFree(sbp->sfp[index]);
    if (sbp->kbp_std)
      sbp->kbp_std[index] = Blast_KarlinBlkFree(sbp->kbp_std[index]);
    if (sbp->kbp_gap_std)
      sbp->kbp_gap_std[index] = Blast_KarlinBlkFree(sbp->kbp_gap_std[index]);
    if (sbp->kbp_psi)
      sbp->kbp_psi[index] = Blast_KarlinBlkFree(sbp->kbp_psi[index]);
    if (sbp->kbp_gap_psi)
      sbp->kbp_gap_psi[index] = Blast_KarlinBlkFree(sbp->kbp_gap_psi[index]);
  }
  if (sbp->kbp_ideal)
    sbp->kbp_ideal = Blast_KarlinBlkFree(sbp->kbp_ideal);
  if (sbp->gbp) 
    sbp->gbp = s_BlastGumbelBlkFree(sbp->gbp);
  sfree(sbp->sfp);
  sfree(sbp->kbp_std);
  sfree(sbp->kbp_psi);
  sfree(sbp->kbp_gap_std);
  sfree(sbp->kbp_gap_psi);
  sbp->matrix = SBlastScoreMatrixFree(sbp->matrix);
  sbp->comments = ListNodeFreeData(sbp->comments);
  sfree(sbp->name);
  sbp->psi_matrix = SPsiBlastScoreMatrixFree(sbp->psi_matrix);
  sfree(sbp->ambiguous_res);
  sfree(sbp);
#endif


  if(ScoreBlkGroup==NULL){return NULL;}

  BlastScoreBlk* lbsbpBlastScoreBlk = (BlastScoreBlk*)calloc(1,sizeof(BlastScoreBlk));

  lbsbpBlastScoreBlk->protein_alphabet=ScoreBlkGroup->protein_alphabet;

  lbsbpBlastScoreBlk->alphabet_code=ScoreBlkGroup->alphabet_code;

  lbsbpBlastScoreBlk->alphabet_size=ScoreBlkGroup->alphabet_size;

  lbsbpBlastScoreBlk->alphabet_start=ScoreBlkGroup->alphabet_start;

  lbsbpBlastScoreBlk->name=ScoreBlkGroup->name;

  lbsbpBlastScoreBlk->comments=NULL;

  ListNode *llndpCurrentNode=NULL;

  ListNode *llndpPreviousNode=NULL;


 {

    if (ScoreBlkGroup->comments) {

      if (lbsbpBlastScoreBlk->comments) {

        llndpPreviousNode=lbsbpBlastScoreBlk->comments;

        llndpCurrentNode=lbsbpBlastScoreBlk->comments->next;

         while (llndpCurrentNode) {

           llndpPreviousNode=llndpCurrentNode;

           llndpCurrentNode=llndpCurrentNode->next;

         }

         llndpPreviousNode->next=ScoreBlkGroup->comments;

       } else {

         lbsbpBlastScoreBlk->comments=ScoreBlkGroup->comments;

      }

    }

  }


  lbsbpBlastScoreBlk->matrix = NULL;

 {

    if (ScoreBlkGroup->matrix) {

      lbsbpBlastScoreBlk->matrix=ScoreBlkGroup->matrix;
      //ceb make deep copy

    }

  }


  lbsbpBlastScoreBlk->psi_matrix = NULL;

  {
    if (ScoreBlkGroup->psi_matrix) {
      lbsbpBlastScoreBlk->psi_matrix=ScoreBlkGroup->psi_matrix;
    }
  }


  lbsbpBlastScoreBlk->matrix_only_scoring=ScoreBlkGroup->matrix_only_scoring;
  lbsbpBlastScoreBlk->complexity_adjusted_scoring=ScoreBlkGroup->complexity_adjusted_scoring;
  lbsbpBlastScoreBlk->loscore=ScoreBlkGroup->loscore;
  lbsbpBlastScoreBlk->hiscore=ScoreBlkGroup->hiscore;
  lbsbpBlastScoreBlk->penalty=ScoreBlkGroup->penalty;
  lbsbpBlastScoreBlk->reward=ScoreBlkGroup->reward;
  lbsbpBlastScoreBlk->scale_factor=ScoreBlkGroup->scale_factor;
  lbsbpBlastScoreBlk->read_in_matrix=ScoreBlkGroup->read_in_matrix;
  


  int lintTotalNumberOfContexts=0;

  {
    lintTotalNumberOfContexts+=ScoreBlkGroup->number_of_contexts;
    //lintTotalNumberOfContexts+=ScoreBlkGroup->number_of_contexts;
  }

  lbsbpBlastScoreBlk->number_of_contexts=lintTotalNumberOfContexts;

  lbsbpBlastScoreBlk->sfp=NULL;

  if (lbsbpBlastScoreBlk->number_of_contexts>0) { 

    if (ScoreBlkGroup->sfp) {

      lbsbpBlastScoreBlk->sfp = (Blast_ScoreFreq**) calloc( lbsbpBlastScoreBlk->number_of_contexts, sizeof(Blast_ScoreFreq*) );

      int lintGlobalContextNumber=-1;

      {

        for (Int4 j=0; j < ScoreBlkGroup->number_of_contexts; ++j) {

          lintGlobalContextNumber++;

          lbsbpBlastScoreBlk->sfp[lintGlobalContextNumber]=ScoreBlkGroup->sfp[j];

        }

      }

    }

  }



  lbsbpBlastScoreBlk->kbp=NULL;

  if (lintTotalNumberOfContexts>0) {
    if (ScoreBlkGroup->kbp) {
      lbsbpBlastScoreBlk->kbp = (Blast_KarlinBlk**) calloc( lbsbpBlastScoreBlk->number_of_contexts, sizeof(Blast_KarlinBlk*) );
      int lintGlobalContextNumber=-1;

       {
        for (Int4 j=0; j < ScoreBlkGroup->number_of_contexts; ++j) {
          lintGlobalContextNumber++;
          lbsbpBlastScoreBlk->kbp[lintGlobalContextNumber]=ScoreBlkGroup->kbp[j];
        }
      }
    }
  }



  lbsbpBlastScoreBlk->kbp_gap=NULL;

  if (lintTotalNumberOfContexts>0) {
    if (ScoreBlkGroup->kbp_gap) {
      lbsbpBlastScoreBlk->kbp_gap = (Blast_KarlinBlk**) calloc( lbsbpBlastScoreBlk->number_of_contexts, sizeof(Blast_KarlinBlk*) );
      int lintGlobalContextNumber=-1;

      {
        for (Int4 j=0; j < ScoreBlkGroup->number_of_contexts; ++j) {
          lintGlobalContextNumber++;
          lbsbpBlastScoreBlk->kbp_gap[lintGlobalContextNumber]=ScoreBlkGroup->kbp_gap[j];
        }
      }
    }
  }


  lbsbpBlastScoreBlk->gbp = NULL;

  {
    if (ScoreBlkGroup->gbp) {
      lbsbpBlastScoreBlk->gbp=ScoreBlkGroup->gbp;
    }
  }



  lbsbpBlastScoreBlk->kbp_std=NULL;

  if (lintTotalNumberOfContexts>0) {
    if (ScoreBlkGroup->kbp_std) {
      lbsbpBlastScoreBlk->kbp_std = (Blast_KarlinBlk**) calloc( lbsbpBlastScoreBlk->number_of_contexts, sizeof(Blast_KarlinBlk*) );
      int lintGlobalContextNumber=-1;

      {
        for (Int4 j=0; j < ScoreBlkGroup->number_of_contexts; ++j) {
          lintGlobalContextNumber++;
          lbsbpBlastScoreBlk->kbp_std[lintGlobalContextNumber]=ScoreBlkGroup->kbp_std[j];
        }
      }
    }
  }



  lbsbpBlastScoreBlk->kbp_psi=NULL;

  if (lintTotalNumberOfContexts>0) {
    if (ScoreBlkGroup->kbp_psi) {
      lbsbpBlastScoreBlk->kbp_psi = (Blast_KarlinBlk**) calloc( lbsbpBlastScoreBlk->number_of_contexts, sizeof(Blast_KarlinBlk*) );
      int lintGlobalContextNumber=-1;

       {
        for (Int4 j=0; j < ScoreBlkGroup->number_of_contexts; ++j ) {
          lintGlobalContextNumber++;
          lbsbpBlastScoreBlk->kbp_psi[lintGlobalContextNumber] = ScoreBlkGroup->kbp_psi[j];
        }
      }
    }
  }


  lbsbpBlastScoreBlk->kbp_gap_std = NULL;

  if (lintTotalNumberOfContexts>0) {
    if (ScoreBlkGroup->kbp_gap_std) {
      lbsbpBlastScoreBlk->kbp_gap_std = (Blast_KarlinBlk**) calloc( lbsbpBlastScoreBlk->number_of_contexts, sizeof(Blast_KarlinBlk*) );
      int lintGlobalContextNumber = -1;

 
      {
        for (Int4 j=0; j < ScoreBlkGroup->number_of_contexts; ++j ) 
        {
          lintGlobalContextNumber++;
          lbsbpBlastScoreBlk->kbp_gap_std[lintGlobalContextNumber] = ScoreBlkGroup->kbp_gap_std[j];
        }
      }
    }
  }


  lbsbpBlastScoreBlk->kbp_gap_psi=NULL;

  if (lintTotalNumberOfContexts>0) 
  {
    if (ScoreBlkGroup->kbp_gap_psi) 
    {
      lbsbpBlastScoreBlk->kbp_gap_psi = (Blast_KarlinBlk**) calloc( lbsbpBlastScoreBlk->number_of_contexts, sizeof(Blast_KarlinBlk*) );
      int lintGlobalContextNumber = -1;


      {
        for (Int4 j=0; j < ScoreBlkGroup->number_of_contexts; ++j ) 
        {
          lintGlobalContextNumber++;
          lbsbpBlastScoreBlk->kbp_gap_psi[lintGlobalContextNumber] = ScoreBlkGroup->kbp_gap_psi[j];
        }
      }
    }
  }


  lbsbpBlastScoreBlk->kbp_ideal = NULL;


  {
    if (ScoreBlkGroup->kbp_ideal) 
    {
      lbsbpBlastScoreBlk->kbp_ideal = ScoreBlkGroup->kbp_ideal;
    }
  }


  lbsbpBlastScoreBlk->ambig_size = ScoreBlkGroup->ambig_size;
  lbsbpBlastScoreBlk->ambig_occupy = ScoreBlkGroup->ambig_occupy;
  lbsbpBlastScoreBlk->ambiguous_res = ScoreBlkGroup->ambiguous_res;
  lbsbpBlastScoreBlk->round_down = ScoreBlkGroup->round_down;

  return lbsbpBlastScoreBlk;
}




















BlastHSPResults *gbhrpfncCombineThreadGroupBlastHSPResults(int aintNumberOfThreadGroups, BlastHSPResults **abhrppBlastHSPResults) 
{
  int lintTotalQueries = 0;
             
  for (int i=0; i < aintNumberOfThreadGroups; i++) {
    //std::cout<<"HSPResults["<<i<<"]->num_queries= "<<abhrppBlastHSPResults[i]->num_queries<<std::endl;  std::cout.flush();
    lintTotalQueries += abhrppBlastHSPResults[i]->num_queries;
  }

  int lintCount = 0;
  //create new BlastHSPResults object and initialize with total number of queries we will be merging 
  BlastHSPResults *lbhrpBlastHSPResults = Blast_HSPResultsNew(lintTotalQueries);
  //std::cout<<"lbhrpBlastHSPResults->num_queries=  " <<lbhrpBlastHSPResults->num_queries <<std::endl;  std::cout.flush();
 
  int lintQueryOffset = 0;
  //loop over thread groups
  for (int ithgrp=0; ithgrp < aintNumberOfThreadGroups; ithgrp++) 
  {
    //loop over number of hitlist arrays
    int lintBlastHSPResultsQueryIndex;
    for (lintBlastHSPResultsQueryIndex = 0; 
	 lintBlastHSPResultsQueryIndex < (abhrppBlastHSPResults[ithgrp]->num_queries); 
	 ++lintBlastHSPResultsQueryIndex) 
    {
#if USE_SHALLOW_COPY
      //shallow copy hitlist_array
      //loop over hsp lists (arrays) and add offset for each thread group
      BlastHitList* old_hlst = abhrppBlastHSPResults[ithgrp]->hitlist_array[lintBlastHSPResultsQueryIndex];
      lbhrpBlastHSPResults->hitlist_array[lintBlastHSPResultsQueryIndex + lintCount] = old_hlst;
      BlastHitList* new_hlst = lbhrpBlastHSPResults->hitlist_array[lintBlastHSPResultsQueryIndex + lintCount];
      for (int ihsplst=0; ihsplst < new_hlst->hsplist_count; ihsplst++) 
      {
        new_hlst->hsplist_array[ihsplst]->query_index += lintQueryOffset;
      }
#else
      //deep copy hitlist array
      BlastHitList* old_hlst = abhrppBlastHSPResults[ithgrp]->hitlist_array[lintBlastHSPResultsQueryIndex];
      //std::cout<<"gbhrpfncCombineThreadGroupBlastHSPResults checkpt 3\n"<<std::cout.flush();
      //Allocate memory for new hitlist
      BlastHitList* new_hlst = Blast_HitListNew(old_hlst->hsplist_count);
      new_hlst->hsplist_count   = old_hlst->hsplist_count;
      new_hlst->hsplist_max     = old_hlst->hsplist_max;
      new_hlst->worst_evalue    = old_hlst->worst_evalue;
      new_hlst->low_score       = old_hlst->low_score;
      new_hlst->heapified       = old_hlst->heapified;
      new_hlst->hsplist_current = old_hlst->hsplist_current;
      //#if 1
      //works
      //allocate memory for hsplist array
      new_hlst->hsplist_array = (BlastHSPList**)calloc( old_hlst->hsplist_current, sizeof(BlastHSPList*) );
      for (int ihsplst=0; ihsplst < new_hlst->hsplist_count; ihsplst++) 
      {
	new_hlst->hsplist_array[ihsplst]= BlastHSPListDuplicate(old_hlst->hsplist_array[ihsplst]);
	//add offset to list query indices
        new_hlst->hsplist_array[ihsplst]->query_index += lintQueryOffset;
      }
      //std::cout<<"gbhrpfncCombineThreadGroupBlastHSPResults checkpt 8\n"<<std::cout.flush();
      //#else
      ////works
      //new_hlst->hsplist_array = old_hlst->hsplist_array;
      //for (int ihsplst=0; ihsplst < new_hlst->hsplist_count; ihsplst++) 
      //{
      //new_hlst->hsplist_array[ihsplst]->query_index += lintQueryOffset;
      //}
      //#endif
      //std::cout<<"gbhrpfncCombineThreadGroupBlastHSPResults checkpt 9\n"<<std::cout.flush();
      lbhrpBlastHSPResults->hitlist_array[lintBlastHSPResultsQueryIndex + lintCount] = new_hlst;
      //std::cout<<"gbhrpfncCombineThreadGroupBlastHSPResults checkpt 10\n"<<std::cout.flush();
#endif
    }
    
    //increment starting query count for next thread group
    lintCount += lintBlastHSPResultsQueryIndex;
    lintQueryOffset += lbhrpBlastHSPResults->num_queries;
  }

  return lbhrpBlastHSPResults; 
}


BlastQueryInfo *gbqipfncCombineThreadGroupBlastQueryInfo(int aintNumberOfThreadGroups, BlastQueryInfo **abqippBlastQueryInfo) 
{
  int lintTotalNumberOfQueries = 0;
  for (int i=0; i < aintNumberOfThreadGroups; i++) {
    lintTotalNumberOfQueries += abqippBlastQueryInfo[i]->num_queries;
  }

  //if restart and we are skipping this batch return an empty struct
  if(lintTotalNumberOfQueries==0){ return  BlastQueryInfoNew(eBlastTypeBlastp,0);}

  int lintCount = 0;

  BlastQueryInfo *lbqipBlastQueryInfo = BlastQueryInfoNew(eBlastTypeBlastp,lintTotalNumberOfQueries);

  lbqipBlastQueryInfo->first_context = 0;
  lbqipBlastQueryInfo->last_context = lintTotalNumberOfQueries-1;
  lbqipBlastQueryInfo->num_queries = lintTotalNumberOfQueries;
  lbqipBlastQueryInfo->max_length = abqippBlastQueryInfo[0]->max_length;
  lbqipBlastQueryInfo->pattern_info = abqippBlastQueryInfo[0]->pattern_info;

  int lintQueryOffset = 0;
  int lintQueryIndex = 0;

  //ceb the following code does a deep copy of each QueryInfo object
  for (int fintThreadGroupIndex = 0; fintThreadGroupIndex<aintNumberOfThreadGroups; fintThreadGroupIndex++) 
  {
    int fintBlastQueryInfoQueryIndex;
    //For each QueryInfo object loop over num_queries or contexts (same in this case)
    for (fintBlastQueryInfoQueryIndex = 0; 
	 fintBlastQueryInfoQueryIndex < (abqippBlastQueryInfo[fintThreadGroupIndex]->num_queries); 
	 ++fintBlastQueryInfoQueryIndex) 
    {
      //ceb record query length for all contexts
      lbqipBlastQueryInfo->contexts[fintBlastQueryInfoQueryIndex + lintCount].query_length =
	abqippBlastQueryInfo[fintThreadGroupIndex]->contexts[fintBlastQueryInfoQueryIndex].query_length;//ceb

      //record new offset for combined context list
      lbqipBlastQueryInfo->contexts[fintBlastQueryInfoQueryIndex + lintCount].query_offset = lintQueryOffset;

      //increment offset
      lintQueryOffset += abqippBlastQueryInfo[fintThreadGroupIndex]->contexts[fintBlastQueryInfoQueryIndex].query_length;

      //compute max length of all combined contexts
      if ( (abqippBlastQueryInfo[fintThreadGroupIndex]->contexts[fintBlastQueryInfoQueryIndex].query_length) > lbqipBlastQueryInfo->max_length ) 
      {
        lbqipBlastQueryInfo->max_length = 
	  abqippBlastQueryInfo[fintThreadGroupIndex]->contexts[fintBlastQueryInfoQueryIndex].query_length;
      }

      lbqipBlastQueryInfo->contexts[fintBlastQueryInfoQueryIndex + lintCount].eff_searchsp = 
	abqippBlastQueryInfo[fintThreadGroupIndex]->contexts[fintBlastQueryInfoQueryIndex].eff_searchsp;

      lbqipBlastQueryInfo->contexts[fintBlastQueryInfoQueryIndex + lintCount].length_adjustment = 
	abqippBlastQueryInfo[fintThreadGroupIndex]->contexts[fintBlastQueryInfoQueryIndex].length_adjustment;

      lbqipBlastQueryInfo->contexts[fintBlastQueryInfoQueryIndex + lintCount].query_index = lintQueryIndex;

      //increment combined context index
      lintQueryIndex++;
      //copy frame (int)
      lbqipBlastQueryInfo->contexts[fintBlastQueryInfoQueryIndex + lintCount].frame = 
	abqippBlastQueryInfo[fintThreadGroupIndex]->contexts[fintBlastQueryInfoQueryIndex].frame;

      lbqipBlastQueryInfo->contexts[fintBlastQueryInfoQueryIndex + lintCount].is_valid = 
	abqippBlastQueryInfo[fintThreadGroupIndex]->contexts[fintBlastQueryInfoQueryIndex].is_valid;
    }
    lintCount += fintBlastQueryInfoQueryIndex;
  }
  return lbqipBlastQueryInfo; 
}




// Functions to combine Queries across replication groups
#if 0
//CRef<CBlastQueryVector> gbqvfncCombineReplicationGroupQueries(int aintNumberOfThreadGroups,
//							      CRef<CBlastQueryVector>* abqvpQueryBatchGroup)						    
//{
//  CRef<CBlastQueryVector> QueryBatchCombined(new CBlastQueryVector);           
//  for(int gid=0; gid< aintNumberOfThreadGroups; ++gid)
//  { 
//    CRef<CBlastQueryVector> QueryBatch= abqvpQueryBatchGroup[gid];
//    for(int i=0; i<QueryBatch->Size(); ++i)            
//    {
//      //std::cout<<"QueryId["<<i<<"]= "<<((*QueryBatch)[i]->GetQueryId())->AsFastaString()<<std::endl;
//      //std::cout<<"QuerySeqLoc["<<i<<"]= "<<((*QueryBatch)[i]->GetQuerySeqLoc())<<std::endl;
//      QueryBatchCombined->AddQuery((*QueryBatch)[i]);
//    }
//  }
//  return QueryBatchCombined;
//}
#endif

CRef<CBlastQueryVector> gbqvfncCombineReplicationGroupQueries(int aintNumberOfThreadGroups,
							      char** buffer,
							      int gbl_bsize,
							      int* lintpBatchStartPosition,
							      int* lintpBatchEndPosition,
							      CBlastInputSourceConfig* iconfig,
							      int QueryBatchSize,
							      int gbl_query_count,
							      CRef<CScope> gbl_scope)
{
  //std::cout<<"gbqvfncCombineReplicationGroupQueries checkpt1 QueryBatchSize="<< QueryBatchSize<<"\n";std::cout.flush();
  //std::cout<<"gbqvfncCombineReplicationGroupQueries checkpt2 gbl_query_count="<<gbl_query_count<<"\n";std::cout.flush();

  //For tid==0 set up global query input stream
  
  //copy thread local query buffers into single global query buffer
  char*gbl_buffer = (char*) calloc( gbl_bsize+1, sizeof(char));
  int sumSize=0;
  //tack each lcl_buffer onto the gbl_buffer
  for ( int i=0; i < num_thread_groups; ++i ){
    int segSize = (lintpBatchEndPosition[i] - lintpBatchStartPosition[i]);
    //std::cout<<"lintpBatchEndPosition["<<i<<"]="<<lintpBatchEndPosition[i] <<" lintpBatchStartPosition["<<i<<"]"<<lintpBatchStartPosition[i]<<std::endl<<std::cout.flush();
    //std::cout<<"segSize= "<<segSize<<std::endl;std::cout.flush();
    sumSize+=segSize;
    strncat(gbl_buffer,&(buffer[i][ lintpBatchStartPosition[i]]),segSize);
    //only do the next step if we are not at the end of the query set
    if (num_thread_groups>1)
    {
      //tack on carriage return between local buffer segments
      strncat(gbl_buffer,"\n",1);sumSize+=1;
      //if(i==num_thread_groups-1){strncat(gbl_buffer,">",1);sumSize+=1;}
    }
    //char*tmp_buffer = (char*) calloc( gbl_bsize+1, sizeof(char));
    //strncat(tmp_buffer,&(buffer[i][ lintpBatchStartPosition[i]]),segSize);
    //std::cout<<"\n\nlcl_buffer["<<i<<"]\n"<<tmp_buffer<<std::endl<<std::endl;std::cout.flush();
  }

  //std::cout<<"sumSize =\n"<<sumSize<<std::endl;std::cout.flush();
  //std::cout<<"\n\ngbl_buffer\n"<< gbl_buffer << std::endl;std::cout.flush();

//ceb
//gbl_buffer has all expected queries
//by the time GetNextSeqBatch is called it is one short of exepected number
//This is causing the bug downstream

  //buffer[i] contains full query input from each file[i]
  membuf gbl_query_membuf(&(gbl_buffer[0]),&(gbl_buffer[sumSize-1]));

  auto_ptr<CNcbiIstream> gbl_query_input_stream;
  gbl_query_input_stream.reset(new CNcbiIstream(&gbl_query_membuf));
  
  if(IsIStreamEmpty(*gbl_query_input_stream)){ERR_POST(Warning << "Query is Empty!");}
  
  iconfig->SetLocalIdCounterInitValue(gbl_query_count+1);//ceb
  
  CBlastFastaInputSource gbl_fasta(*gbl_query_input_stream, *iconfig);
  
  CBlastInput gbl_input(&gbl_fasta, QueryBatchSize*num_thread_groups);
  //test to make sure we get all batches in gbl_input
  //gbl_input.SetBatchSize(QueryBatchSize*num_thread_groups*10);
  //CRef<CBlastQueryVector> gbl_query_batch(gbl_input.GetNextSeqBatch(*gbl_scope));
  CRef<CBlastQueryVector> gbl_query_batch(gbl_input.GetAllSeqs(*gbl_scope));

  //std::cout<<"gbqvfncCombineReplicationGroupQueries checkpt3 gbl_query_batch->Size()="<<gbl_query_batch->Size()<<"\n";std::cout.flush();
  //ceb test 
  //for(int i=0;i<gbl_query_batch->Size();++i)
  //{
  //  std::cout<<"gbl_query_batch["<<i<<"].GetQueryId() ="<<(((gbl_query_batch.GetObject())[i].GetObject()).GetQueryId().GetObject()).AsFastaString()<<"\n"<<std::cout.flush();
  //  //std::cout<<"gbl_query_batch["<<i<<"].GetLength() ="<<(((gbl_query_batch.GetObject())[i].GetObject()).GetLength())<<"\n"<<std::cout.flush();
															//}

  return gbl_query_batch;
};




//Function to create a search resluts set from the  4 combined objects

CRef<CSearchResultSet> gsrsfncCollectResults(BlastHSPResults* hsp_results,
					     BlastQueryInfo* m_QueryInfo,
					     BlastScoreBlk* m_ScoreBlk,
					     CRef<IQueryFactory> m_QueryFactory,
					     const CBlastOptions* m_Options,
					     CRef<IBlastSeqInfoSrc> m_SeqInfoSrc)
{
  //std::cout<<"gsrsfncCollectResults checkpt 1 \n";std::cout.flush();
    int hitlist_size_backup = m_Options->GetHitlistSize();
    EResultType m_ResultType = eDatabaseSearch;
    TSearchMessages m_Messages;// = ;//TSearchMessages&

    // This is the data resulting from the traceback phase (before it is converted to ASN.1).
    // We wrap it this way so it is released even if an exception is thrown below.
    CRef< CStructWrapper<BlastHSPResults> > HspResults;
    HspResults.Reset(WrapStruct(hsp_results, Blast_HSPResultsFree));

    _ASSERT(m_SeqInfoSrc);
    _ASSERT(m_QueryFactory);

//ceb causing segfault downstream
    CRef<ILocalQueryData> qdata = m_QueryFactory->MakeLocalQueryData(m_Options);



    m_SeqInfoSrc->GarbageCollect();
    vector<TSeqLocInfoVector> subj_masks;
    //converts hsp object into align vector
    //std::cout<<"gsrsfncCollectResults checkpt 5 \n";std::cout.flush();

    TSeqAlignVector aligns =
        LocalBlastResults2SeqAlign(hsp_results,
                                   *qdata,
                                   *m_SeqInfoSrc,
                                   m_Options->GetProgramType(),
                                   m_Options->GetGappedMode(),
                                   m_Options->GetOutOfFrameMode(),
                                   subj_masks,
                                   m_ResultType);
    //std::cout<<"gsrsfncCollectResults checkpt 6 \n";std::cout.flush();

    //query factory
    vector< CConstRef<CSeq_id> > query_ids;
    query_ids.reserve(aligns.size());
    //std::cout<<"gsrsfncCollectResults checkpt 7 \n";std::cout.flush();

    for (size_t i = 0; i < qdata->GetNumQueries(); i++){
      //std::cout<<"gbl Q_id: "<<qdata->GetSeq_loc(i)->GetId()<<std::endl;
        query_ids.push_back(CConstRef<CSeq_id>(qdata->GetSeq_loc(i)->GetId()));
    }
    //std::cout<<"gsrsfncCollectResults checkpt 8 \n";std::cout.flush();

    return BlastBuildSearchResultSet(query_ids,
                                     m_ScoreBlk,//->GetPointer(),
                                     m_QueryInfo,
                                     m_Options->GetProgramType(),
                                     aligns,
                                     m_Messages,
                                     subj_masks,
                                     NULL,
                                     m_ResultType);
}


#endif /* SKIP_DOXYGEN_PROCESSING */
