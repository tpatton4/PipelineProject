# PipelineProject

This repo has Tessa Patton's COMP 438 Pipeline Project for Track 2. This script uses SPAdes to assemble transcriptome reads from four Human cytomegalovirus (HCMV) patients by filtering reads with Bowtie2 mapping, assembling with SPAdes, and blasting the lognest contig against the Betaherpesvirinae subfamily. Briefly, input files were obtained by retrieveing donor transcriptomes from SRA and converting to paired-end fastq files with fastq-dump. An index for HCMV was created with bowtie2-build from the ViralProj14559 HCMV fasta (NCBI accession NC_006273.2). 

# Dependencies
You should add all dependencies to your path prior to running this script  
Biopython  
BLAST+  
Bowtie  
SPAdes  

# Instructions
Download pipelineproject_28FEB24.py and the HCMV_input folder to your working directory. To run the script use the command: python PipelineProject_Tessa_Patton.py
All outputs will be written to a directory called PipelineProject_Tessa_Patton in the log file called PipelineProject.log.
