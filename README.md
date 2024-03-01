# PipelineProject

This repo has Tessa Patton's COMP 438 Pipeline Project for Track 2: Genome Assembly. This script uses SPAdes to assemble transcriptome reads from four Human cytomegalovirus (HCMV) patients by filtering reads with Bowtie2 mapping, assembling with SPAdes, and blasting the lognest contig against the Betaherpesvirinae subfamily. Briefly, input files were obtained by retrieveing donor transcriptomes from SRA and converting to paired-end fastq files with fastq-dump. An index for HCMV was prepared for you with bowtie2-build from the ViralProj14559 HCMV fasta (NCBI accession NC_006273.2). 

Test data is included in the HCMV_input folder to test this script. This test data contains only the first 10000 reads from each fastq file.

# Dependencies
You should add all dependencies to your path prior to running this script 

- [Python](https://www.python.org/downloads/)
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) 
- [SPAdes](http://cab.spbu.ru/software/spades/) 
- [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) 

# Instructions
Download pipelineproject_28FEB24.py and the HCMV_input folder to your working directory. Make sure that your directory called HCMV_input contains all of the same files as mine (HCMV index files and sample .fastq files). To run the script use the command: `python PipelineProject_Tessa_Patton.py`

All outputs will be written to a directory called PipelineProject_Tessa_Patton including a log file called PipelineProject.log.
