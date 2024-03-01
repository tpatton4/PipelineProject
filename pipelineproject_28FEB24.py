import os
import gzip
from Bio import SeqIO

def find_fastq_pairs(input_dir): #find pairs of fastq files downloaded input directory
    fastq_files = [f for f in os.listdir(input_dir) if f.endswith('.fastq')] #only .fastq files. there are index files here, too
    pairs = {} #make a list for the pairs
    for file in fastq_files:
        prefix = file.rsplit('_', 1)[0] #split by the _ and use the first part as the unique ID
        if prefix in pairs:
            pairs[prefix].append(file) #append file to list if there's another with same prefix
        else:
            pairs[prefix] = [file] #otherwise add it bc it's the first
    for prefix in pairs:
        pairs[prefix].sort() #sort so _1 comes before _2
    return pairs

def run_bowtie2(input_dir, index_prefix, log_file):
    output_files = [] #empty set for output files
    fastq_pairs = find_fastq_pairs(input_dir) #use prev function to look at the input directory from github
    for prefix, files in fastq_pairs.items():
        forward_read, reverse_read = [os.path.join(input_dir, f)  for f in files] #assign _1 and _2 to variables forward and reverse read
        aligned_output_prefix = os.path.join(input_dir, f"{prefix}_aligned_conc") #filename for aligned reads
        output_sam = f"{aligned_output_prefix}.sam" #filename for output sam file
        bowtie2_command = f"bowtie2 -x {index_prefix} -1 {forward_read} -2 {reverse_read} --al-conc-gz {aligned_output_prefix}%.fastq.gz -S {output_sam}" #make bowtie2 command with flag to save only mapped reads in new fastq.gz files
        os.system(bowtie2_command) #run the command!
        ####count reads####
        input_read_count = 0
        for file in [forward_read, reverse_read]:
            with gzip.open(file, "rt") if file.endswith('.gz') else open(file, "r") as handle: #should be .gz files, but maybe .fasta? 
                input_read_count += sum(1 for _ in SeqIO.parse(handle, "fastq")) / 2
        aligned_reads_1 = f"{aligned_output_prefix}1.fastq.gz"
        aligned_reads_2 = f"{aligned_output_prefix}2.fastq.gz"
        import gzip
        with gzip.open(aligned_reads_1, "rt") as handle:
            aligned_pairs_count = sum(1 for _ in SeqIO.parse(handle, "fastq")) #count the saved reads by parsing fastq files with SeqIO

        with open(log_file, "a") as log:
            log.write(f"Donor {prefix} had {int(input_read_count)} read pairs before Bowtie2 filtering and {aligned_pairs_count} read pairs after\n") #ugh i know the prefix is wrong - come back to this with a dict to link to metadata if time. 
        
        output_files.extend([aligned_reads_1, aligned_reads_2])
    return output_files

def run_spades(input_dir, output_dir, log_file):
    aligned_fastq_files = [f for f in os.listdir(input_dir) if f.endswith('_conc1.fastq.gz') or f.endswith('_conc2.fastq.gz')] #get pairs of fastq.gz files from bowtie2 output
    aligned_fastq_files.sort() #ensure pairs are adjacent
    pe_args = [] #to store pairs 
    for i in range(0, len(aligned_fastq_files), 2):
        pe_args.append(f"-1 {input_dir}/{aligned_fastq_files[i]} -2 {input_dir}/{aligned_fastq_files[i+1]}") #store filtered fastq files with the -1 and -2 flags for spades
    spades_command = f"spades.py -k 77,99,127 -t 2 -o {output_dir} --only-assembler --quiet " + " ".join(pe_args) #set up spades command - flags pulled from compbio ppt
    os.system(spades_command) #run the command
    with open(log_file, "a") as log:
        log.write(f"{spades_command}\n") #write this to the log file

def analyze_assembly(output_dir, log_file, db_name): #find contigs with length > 1000
    contigs_file = os.path.join(output_dir, "contigs.fasta") #get the contigs from output of SPADes
    contigs = list(SeqIO.parse(contigs_file, "fasta")) #parse so we can see length
    long_contigs = [contig for contig in contigs if len(contig) > 1000] #get contigs longer than 1000 bp
    longest_contig = max(contigs, key=lambda x: len(x)) #get the longest contig while we're here
    longest_contig_file = os.path.join(output_dir, "longest_contig.fasta")
    SeqIO.write(longest_contig, longest_contig_file, "fasta")
    total_length = sum(len(contig) for contig in long_contigs) #add the filtered contig lengths together to get total length
    makeblastdb_command = f"makeblastdb -in {longest_contig} -dbtype nucl -out {db_name}" #make the blast database command
    os.system(makeblastdb_command) #run it
    with open(log_file, "a") as log: #write those in the log file
        log.write(f"There are {len(long_contigs)} contigs > 1000 bp in the assembly\n")
        log.write(f"There are {total_length} bp in the assembly\n")
    return long_contigs #return list of long contigs to blast
#tomorrow - finish blasting longest contig to the database

def blast_longest_contig(long_contigs, db_name, log_file):
    longest_contig = max(long_contigs, key=lambda x: len(x)) #uhh not sure why I put this stuff in the other function but oh well lol
    with open("longest_contig.fasta", "w") as output_file:
        SeqIO.write(longest_contig, output_file, "fasta")
    blast_command = f"blastn -query longest_contig.fasta -db {db_name} -outfmt '6 sacc pident length qstart qend sstart send bitscore evalue stitle' -max_target_seqs 10 -max_hsps 1 > blast_output.txt" #use 6 for tab delimited output
    os.system(blast_command) #run blast
    with open("blast_output.txt", "r") as blast_output, open(log_file, "a") as log: #open the blast output to read and the log file to append that data
        log.write("sacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tbitscore\tevalue\tstitle\n")
        for line in blast_output:
            log.write(line) #put each line from blast output into the log file

directory = os.getcwd()
download_command = f"datasets download virus genome taxon betaherpesvirinae --refseq --include genome" #use NCBI datasets tool 
os.system(download_command) #run the command
os.system(f"unzip -o ncbi_dataset.zip") #unzip the file
os.system(f"makeblastdb -in ncbi_dataset/data/genomic.fna -out HCMV -title HCMV -dbtype nucl")

def main():
    input_dir = "HCMV_input"
    index_prefix = os.path.join(input_dir, "HCMV")
    output_dir = "PipelineProject_Tessa_Patton"
    log_file = "PipelineProject.log"
    db_name = "HCMV"
    
    #clear or create the log file
    open(log_file, 'w').close()
    #run Bowtie2
    bowtie2_outputs = run_bowtie2(input_dir, index_prefix, log_file)
    #run SPADes
    run_spades(input_dir, output_dir, log_file)
    #analyze the assembly
    long_contigs = analyze_assembly(output_dir, log_file, db_name)
    #blast the longest contig
    blast_longest_contig(long_contigs, db_name, log_file)

if __name__ == "__main__":
    main()

os.system(f"mv PipelineProject.log PipelineProject_Tessa_Patton/")
os.system(f"mv HCMV_input/*.fastq.gz PipelineProject_Tessa_Patton/")
