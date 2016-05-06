#HPC-BLAST


HPC-BLAST was designed to span the performance space of a heterogenous, Intel Xeon Phi, cluster.
It has been optimized to run efficiently on a single Xeon Phi, a Single Xeon pair,
Multiple Xeon Phis on ther same node, both the Xeon pair and all the Xeon Phis on a single
node, any number of Xeon Phis on any number of nodes, any number of Xeon pairs
on any number of nodes, and on both any number of Xeon pairs and any number of
Xeon Phis on any number of nodes.

In this repository, you will find the source to build HPC-BLAST, the source to build the auxiliary tools associated with HPC-BLAST, a User Guide 
showing how to install and run HPC-BLAST and the auxiliary tools associated with HPC-BLAST, and a Best Practices Guide, showing how to tune HPC-BLAST. 



How to Get HPC-BLAST
--------------------

For the most up-to-date version of HPC-BLAST,
and associated auxiliary tools, we recommend that you
download from our Git repository. This can be accomplished via
cloning the repository from the command line, or by downloading a zip
from our GitHub page.

Git Repository Clone:

  Use the following command to clone HPC-BLAST to your machine:

  >$ git clone https://github.com/aace/HPC-BLAST.git

  Once cloned, you can update the code to the newest version using the following command (when in the HPC-BLAST directory):
   
  >$ git pull

Git Zip Download:
    
  Simply use the "Download Zip" option on our webpage at:
    
  https://github.com/aace/HPC-BLAST



Installation and Usage
----------------------

Installation and usage is covered in detail in: 

HPC-BLAST-User-Manual.pdf. in this repository.



Advanced Usage
--------------

Tuning HPC-BLAST is covered in detail in:

HPC-BLAST-Best-Practices.pdf, in this repository.



HPC-BLAST
---------

HPC-BLAST uses the message passing interface (MPI) to partition both the database and 
the query file roughly at the node level, and OpenMP to further 
partition both the database and the query file at the thead level. 

To the best of our knowledge, we are the first to implement a distributed, NCBI compliant,
BLAST+ (C++ toolkit) code, for Intel Xeon Phi clusters.

Best results to date indicate that HPC-BLAST runs approximately 40% faster for blastn 
and 45% faster for blastp than NCBI BLAST+ on Intel Xeon processors and around 773% 
faster for blastn and 225% for blastp than NCBI BLAST+ on Intel Xeon Phi coprocessors. 
HPC-BLAST also demonstrates near constant 90% parallel efficiency with respect to weak scaling 
across 128 Intel Xeon Phi coprocessors. Combined, these results indicate that HPC-BLAST 
offers substantial performance improvements over NCBI BLAST+ on highly parallel computing platforms.



Authors
-------

Shane E. Sawyer

Mitchel D. Horton

Chad Burdyshaw



Contact Information
-------------------

R. Glenn Brook 

glenn-brook@tennessee.edu

JICS, University of Tennessee


