###
# DISCLAIMER:
# This was the workflow to find ONSEN new insertions in re-sequenced transpositionally-competent nrpd1-3 lines.
# I used standard open-source software in Linux command-line, and custom-made python scripts (I kept here the original names of the python script's files to match my records)
# It also involved some manual organization of files and manual assessment of results.
# Although it solved the problems faced at the time, note that this is a beginner's code, as I'm not a professional bioinformatician or data scientist.
# As such, it may probably look overly complicated; and other more advanced solutions are likely simpler or more effective.
# The workflow also display some historical contingencies, as I was learning while doing.
# None of this work would be possible without the collaboration with Dr. Herve Gaubert, and the kind bioinformatics support of Varodom Charoensawan, Hugo Tavares, Jeremy Gruel
#
# Dr. Diego H. Sanchez (diego.sanchez@slcu.cam.ac.uk)
###


#######################################################################
# STEP 1:
# Look for PE pairs, in which a mate hits the genome and the other mate
# hits ONSEN family of TEs.
#######################################################################

#########
# A) Mask genome with ONSEN sequences.
#########

# Generate a bed file with ONSEN coordinates and retrieve sequence as fasta file.
$bedtools getfasta -fi ~/PATH/TAIR10_chr_all.fas \
-bed ~/PATH/family_COPIA78_coordinates.txt \
-fo ~/PATH/fasta_ONSEN.fa -name &

# Blast ONSEN against Arabidopsis genome, then
# manually recover areas with at least 75 bp alignment lenght, and generate a .txt file with chr, start, end data.
$blastn -subject ~/PATH/TAIR10_chr_all.fas -query ~/PATH/fasta_ONSEN.fa \
-max_target_seqs 20000 -outfmt "6" -out ~/PATH/blast_table_ONSEN_vs_TAIR10_chr.blast &

# Generate a bed file of these areas; use python script to generate a .bed file
# with a coordinate as name from chr,start,end .txt file.
# FILE:	1_generates bed4 from chr_start_end list for LTRs.py

#Script:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import time
start_time   = time.time()  #to tell de time it take to process

in_filename  = "~/PATH/input.txt"
out_filename = "~/PATH/output_ONSEN.bed" # out_file in .bed4 format (chr,start,end,ID)
MIN_SIZE     = 75 # minimun size to consider

#####################################################################
## define the 'CountLines' function to count the lines a the file ##
#####################################################################
def CountLines (filename):    
    with open (filename, 'r') as myfile: # Open input gene exppression file (tab-separated)
        count = sum(1 for line in myfile)
    return count  

# open file and count
total_counts = CountLines (in_filename)
print total_counts, "lines in input file", "  --- %s minutes run---" % ((time.time() - start_time) / 60)

# open in and out files, and parse to generate the bed file
out_bed_file = open(out_filename,"w")     
in_list_file = open(in_filename, "r")     

for genes in in_list_file:
    
    chromosome, start, end = genes.strip("\n").split('\t')
    start_ = int(start)
    end_   = int(end)
    
    if start_ > end_: # Fix problems arising from start > end sometimes present in asseblies
        end   = start_
        start = end_
    
    start_name = str(start)
    end_name   = str(end)
    
    if int(end_name) - int(start_name) >= MIN_SIZE:
    
        name = chromosome + ':' + start_name + '-' + end_name # Generate ID
        print name
            
        out_line = [chromosome, '\t', start_name, '\t', end_name, '\t', name, '\n'] # Generates de components of the .bed4 file
        for item in out_line: 
            out_bed_file.write(item) # Writing each component to a new .tab delimited file, .bed6 format
            #out_bed_file.write('\n')  

in_list_file.close()
out_bed_file.close() 

total_counts = CountLines (out_filename)
print total_counts, "lines in output file", "  --- %s minutes run---" % ((time.time() - start_time) / 60)

print "Final: --- %s minutes run---" % ((time.time() - start_time) / 60)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Presort .bed by chromosome and then by start position, and then merge overlapping areas in the .bed file
$sort -k1,1 -k2,2n output_ONSEN.bed > output_ONSEN_sorted.bed &
$bedtools merge -d 1 -i output_ONSEN_sorted.bed > output_ONSEN_sorted_merged.bed &

# Use python script to add a name in the .bed file from chr,start,end.
# FILE: 1_generates bed4 from chr_start_end list for LTRs.py
Input file: 	'output_ONSEN_sorted_merged.bed'
Output file: 	'output_ONSEN_sorted_merged_name.bed'

# Mask genome with .bed file
$bedtools maskfasta -fi ~/PATH/TAIR10_chr_all.fas \
-bed ~/PATH/output_ONSEN_sorted_merged_name.bed \
-fo ~/PATH/TAIR10_chr_all_maskedONSEN.fa &

# Generate bowtie2 index of the masked genome
$bowtie2-build TAIR10_chr_all_maskedONSEN.fa TAIR10_chr_all_maskedONSEN &

#########
# B) Map PE sequencing to this masked genome using bowtie2
#########

#### Trim adaptors, using a fast file with them
$java -jar ~/PATH/trimmomatic-testing.jar PE -threads 4 \
~/PATH/SAMPLE_1.fq.gz \
~/PATH/SAMPLE_2.fq.gz \
~/PATH/SAMPLE_1_trimmo_paired.fq \
~/PATH/SAMPLE_1_trimmo_unpaired.fq \
~/PATH/SAMPLE_2_trimmo_paired.fq \
~/PATH/SAMPLE_2_trimmo_unpaired.fq \
ILLUMINACLIP:~/PATH/Trimmomatic-0.32/adapters/TruSeq2-PE.fa:2:10:5:1 &

# Check fastqc and gzip
$fastqc ~/PATH/SAMPLE_1_trimmo_paired.fq &
$fastqc ~/PATH/SAMPLE_1_trimmo_paired.fq &
# compress files
$gzip ~/PATH/SAMPLE_1_trimmo_paired.fq &
$gzip ~/PATH/SAMPLE_2_trimmo_paired.fq &

# Mapping with bowtie2 (very sensitive).
$bowtie2 -p 6 --quiet --very-sensitive -X 1000 --non-deterministic \
-x ~/PATH/TAIR10_chr_all_maskedONSEN \
-1 ~/PATH/Col_a_1_trimmo_paired.fq.gz \
-2 ~/PATH/Col_a_2_trimmo_paired.fq.gz \
-S ~/PATH/SAMPLE_trimmo_paired_bowtie2_verysensitive_TAIR10_chr_all_maskedONSEN.sam &

# Convert sam to bam, sort
$samtools view -bS ~/PATH/SAMPLE_trimmo_paired_bowtie2_verysensitive_TAIR10_chr_all_maskedONSEN.sam \
> ~/PATH/SAMPLE_trimmo_paired_bowtie2_verysensitive_TAIR10_chr_all_maskedONSEN.bam &

$samtools sort ~/PATH/SAMPLE_trimmo_paired_bowtie2_verysensitive_TAIR10_chr_all_maskedONSEN.bam \
~/PATH/SAMPLE_trimmo_paired_bowtie2_verysensitive_TAIR10_chr_all_maskedONSEN_sorted &

# Remove duplicates
$java -Xmx4g -jar ~/PATH/picard-tools-1.103/MarkDuplicates.jar \
INPUT=SAMPLE_trimmo_paired_bowtie2_verysensitive_TAIR10_chr_all_maskedONSEN_sorted.bam \
OUTPUT=SAMPLE_trimmo_paired_bowtie2_verysensitive_TAIR10_chr_all_maskedONSEN_sorted_rmdup_picard.bam \
METRICS_FILE=dup.txt VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true QUIET=true &

# Indexing
$samtools index ~/PATH/SAMPLE_trimmo_paired_bowtie2_verysensitive_TAIR10_chr_all_maskedONSEN_sorted.bam &

# Select for read with mapped pair and unmapped mate
#output .sam
$samtools view -f 8 -F 4 -h ~/PATH/SAMPLE_trimmo_paired_bowtie2_verysensitive_TAIR10_chr_all_maskedONSEN_sorted_rmdup_picard.bam \
> ~/PATH/Reads_SAMPLE_GenomethatoneUnmapped.sam &

#output .bam
$samtools view -f 8 -F 4 -b -h ~/PATH/SAMPLE_trimmo_paired_bowtie2_verysensitive_TAIR10_chr_all_maskedONSEN_sorted_rmdup_picard.bam \
> ~/PATH/Reads_SAMPLE_GenomethatoneUnmapped.bam &

# Transform bam to bed
$bamToBed -i ~/PATH/Reads_SAMPLE_GenomethatoneUnmapped.bam \
> ~/PATH/Reads_SAMPLE_GenomethatoneUnmapped.bed &

# Transform bam to fastq
$bedtools bamtofastq -i ~/PATH/Reads_SAMPLE_GenomethatoneUnmapped.bam \
-fq ~/PATH/Reads_SAMPLE_GenomethatoneUnmapped.fastq &

# Recover unmapped reads
&samtools view -b -f 4 ~/PATH/SAMPLE_trimmo_paired_bowtie2_verysensitive_TAIR10_chr_all_maskedONSEN_sorted_rmdup_picard.bam \
> ~/PATH/Reads_SAMPLE_GenomeUnmapped.bam &

# Generate .fastq and .fasta with the unmapped reads
$bedtools bamtofastq -i ~/PATH/Reads_SAMPLE_GenomeUnmapped.bam \
-fq ~/PATH/Reads_SAMPLE_GenomeUnmapped.fastq &
$awk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' Reads_SAMPLE_GenomeUnmapped.fastq > Reads_SAMPLE_GenomeUnmapped.fasta &

#########
# C) Map PE sequencing to ONSEN
#########

# Generate the reverse complement from ONSEN fasta file (fasta_ONSEN.fa) and
# include it in the fasta file (puts both foward and reverse strand for mapping).
# I learnt this is unneccesary, as bowtie2 seeds with both strands.
# However, just in case, I introduced this step in the original workflow.
# FILE:	6_input fasta sequences ONSEN generate reverse complement_include both in outcome_with SeqRecord.py
#
# Script:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import time
from Bio import SeqIO
from Bio.Seq import Seq

start_time = time.time() # To tell the time it take to process
count = 0 # To count the number of sequences

InFile  = "fasta_ONSEN.fa" # Input file
OutFile = "fasta_ONSEN_with_reverse.fasta"	# Output file

#  Opens the fasta
path          = "~/PATH/"
input_handle  = open("%s%s" % (path, InFile), "rU") # Opens .fasta with sequences
output_handle = open("%s%s" % (path, OutFile), "w") # Opens output file

for record in SeqIO.parse(input_handle, "fasta") :
    
    name = record.id[:]
    
    # write first the lead strand
    record.id          = name+'_foward'
    record.name        = name + '_foward'
    record.description = name + '_foward'
    SeqIO.write(record, output_handle, "fasta") # Writing the lead strand
    print record.id,
    
    # write the reverse strand
    sequence           = record.seq [:]
    reverse_sequence   = sequence.reverse_complement()
    record.seq         = reverse_sequence
    record.id          = name + '_reverse'
    record.name        = name + '_reverse'
    record.description = name + '_reverse'
    SeqIO.write(record, output_handle, "fasta") # Writing the reverse strand
    print record.id

    count += 1
    
input_handle.close()
output_handle.close()

print ' total number of sequences= ', count
print '.fasta file modified in         ', "--- %s minutes run ---" % ((time.time() - start_time) / 60)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Generate bowtie2 for the outfasta file
bowtie2-build fasta_ONSEN_with_reverse.fasta fasta_ONSEN &

#########
# C BIS) Map to ONSEN the unmapped reads from the genome
#########

# Mapping with bowtie2
$bowtie2 -p 6 --quiet --very-sensitive -X 1000 --non-deterministic \
-x ~/PATH/fasta_ONSEN \
-U ~/PATH/Reads_SAMPLE_GenomeUnmapped.fastq \
-S ~/PATH/Reads_SAMPLE_GenomeUnmapped_verysensitive_fasta_ONSEN.sam &

# Convert sam to bam
$samtools view -bS ~/PATH/Reads_SAMPLE_GenomeUnmapped_verysensitive_fasta_ONSEN.sam \
> ~/PATH/Reads_SAMPLE_GenomeUnmapped_verysensitive_fasta_ONSEN.bam &

# sort
$samtools sort ~/PATH/Reads_SAMPLE_GenomeUnmapped_verysensitive_fasta_ONSEN.bam \
~/PATH/Reads_SAMPLE_GenomeUnmapped_verysensitive_fasta_ONSEN_sorted &

# indexing
$samtools index ~/PATH/Reads_SAMPLE_GenomeUnmapped_verysensitive_fasta_ONSEN_sorted.bam &

# Select for reads that mapped to ONSEN (-b option of the outcome is in .bam)
# Output .sam (just in case) 
$samtools view -F 4 -h ~/PATH/Reads_SAMPLE_GenomeUnmapped_verysensitive_fasta_ONSEN_sorted.bam \
> ~/PATH/Reads_SAMPLE_GenomeUnmapped_verysensitive_fasta_ONSEN_sorted_MappedtoONSEN.sam &

# Output .bam
$samtools view -F 4 -b -h ~/PATH/Reads_SAMPLE_GenomeUnmapped_verysensitive_fasta_ONSEN_sorted.bam \
> ~/PATH/Reads_SAMPLE_GenomeUnmapped_verysensitive_fasta_ONSEN_sorted_MappedtoONSEN.bam &

# Output .bed 
$bedtools bamtobed -i ~/PATH/Reads_SAMPLE_GenomeUnmapped_verysensitive_fasta_ONSEN_sorted_MappedtoONSEN.bam \
> ~/PATH/Reads_SAMPLE_GenomeUnmapped_verysensitive_fasta_ONSEN_sorted_MappedtoONSEN.bed &

#########
# D) Compare the outcome between Reads_ONSEN_GenomethatoneUnmapped.sam and Reads_ONSEN_thatoneUnmapped.sam
# To make it fast and possible the use of my standar PC, mapping files were used in .bed format
#########

# Use custome made python script to find same read name in both data sets
# FILE:	12bis_withCounter_parse_2xbed6_files_output_list_of_shared_reads_bed6.py
# Script:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import time
start_time = time.time() # To tell de time it take to process
import numpy as np
from collections import Counter

# Files and path
FilenameBed = "Reads_SAMPLE_GenomethatoneUnmapped.bed" 						
Filename2   = "Reads_SAMPLE_GenomeUnmapped_verysensitive_fasta_ONSEN_sorted_MappedtoONSEN.bed"
OutFile     = 'report_SAMPLE_SAM_comparison_GenomeUnmapped_ONSENmapped.txt' # This is a text file with the results of the final parsing, showing best hits
OutBEDFile  = "Out_Reads_SAMPLE_GenomethatoneUnmapped_vs_GenomeUnmapped-but-ONSENmapped.bed"
path        = "~/PATH/"
d_type      = 'a48' # Max lenght of read names in number of characters

# define CountLines fx to count the lines a the file
def CountLines (filename):
    with open ("%s%s" % (path, filename), 'r') as myfile: # Open the input tab-file with gene exppression
        count = sum(1 for line in myfile)
    return count

# defines fx ParseBED6
def ParseBED6 (input_file, number_lines):

    BED_Matrix = np.array(range(number_lines), dtype = d_type) # Generates a matrix with strings of at most d_type characters, to store names of the reads
    i = 0 # To count position in matrix
    for record in input_file: 
        if record.startswith("@"): continue # Historical contingency again, disregards SAM header
        else:                      
            line          = record.strip().split("\t")
            name_bed      = line[3]
            name_bed      = name_bed.replace("/1", "")
            name_bed      = name_bed.replace("/2", "")
            BED_Matrix[i] = name_bed  
            i = i + 1
    BED_Matrix = np.unique(BED_Matrix) # Removes duplicated reads, sorts
    return BED_Matrix

#####################
# MAIN 1st part
#####################

# parse the 1st File
counts_1       = CountLines (FilenameBed)
print 'fist file count: ', counts_1
input_handle_1 = open("%s%s" % (path,FilenameBed), "r") # Opens .bed file
Data_File1     = ParseBED6 (input_handle_1, counts_1)  # Appplies ParseBED6 fx
input_handle_1.close()
print "fist file BED6 parsed ", "--- %s minutes run ---" % ((time.time() - start_time) / 60)

# parse the 2st File
counts_2       = CountLines (Filename2)
print 'second file count: ', counts_2
input_handle_2 = open("%s%s" %(path,Filename2), "r") # Opens .bed file
Data_File2     = ParseBED6 (input_handle_2, counts_2) # Appplies ParseBED6 fx
input_handle_2.close()
print "second file BED6 parsed ", "--- %s minutes run ---" % ((time.time() - start_time) / 60)
                   
# concatenate matrixes
Final_Matrix = np.concatenate((Data_File1, Data_File2), axis = 0)

# Match data
count = 0
Shared_array = np.array(([item for item,count in Counter(Final_Matrix).iteritems() if count > 1]), dtype = d_type)
print "Total reads match: ", len(Shared_array), " in time: --- %s minutes run ---" %((time.time()-start_time)/60)              

# write report
output_handle = open("%s%s" % (path,OutFile), "w") # Opens the file to be written
for out_line in Shared_array:       
            output_handle.write("%s\n" %out_line)# Writes line plus \n    
output_handle.close()

#####################
# MAIN 2nd part
#####################

# count lines in files
counts_text = CountLines (OutFile)
print 'Number of lines in text file: ', counts_text, "--- %s minutes run ---" %((time.time() - start_time) / 60)
counts_bed = CountLines (FilenameBed)
print 'Number of lines in BED file: ', counts_bed, "--- %s minutes run ---" %((time.time() - start_time) / 60)

# input txt data in matrix
Text_List = np.zeros( (counts_text), dtype = d_type) # Makes a matrix with at most counts_text items
input_text_handle = open("%s%s" %(path,OutFile), "r") # Opens the file to be written
count = 0
for record in input_text_handle:   
    record = record.strip()
    record = record.strip("\t")   
    record = record.strip("\n")
    Text_List[count] = record                              
    count += 1
input_text_handle.close()    
Text_List.sort()
print 'text data in matrix     ', "--- %s minutes run ---" %((time.time() - start_time) / 60)

# input data from .bed in matrix
Reads_List       = np.zeros( (counts_bed), dtype = d_type) # This matrix keeps the names, to find indices later
Reads_Data       = np.zeros( (counts_bed), dtype = np.ndarray) # This matrix keeps the data
input_bed_handle = open("%s%s" %(path, FilenameBed), "r") # Opens the file to be written

count = 0
for line in input_bed_handle:
    Chr,start,end,name,dot,strand = line.strip().split("\t") # unpack data from .bed6 columns to recover name of reads
    name              = name.replace("/1","")
    name              = name.replace("/2","")
    Reads_List[count] = str(name) # Keeps the read name in List matrix to use as index
    Reads_Data[count] = [line] # Writes the .bed6 line in Data matrix - list!
    count += 1
input_bed_handle.close()    
print 'bed data in matrix     ', "--- %s minutes run ---" %((time.time() - start_time) / 60)

# loop .txt data in the bed data
Final_List = np.zeros( (counts_bed), dtype=np.ndarray) # This matrix will store the shred reads to write, reshape later
index = 0
count = 0
for record in Text_List:
    if record in Reads_List: 
        index = np.where(Reads_List == record) # Find position of key in Reads_List
        for replicates in index: # Loops required to take into account possible replicated mapping
            for positions in replicates:
                new_index = int(positions) # Position of key in Reads_List
                Final_List[count] = Reads_Data[new_index]
                count += 1
                
Final_List.resize(count, 1)
Final_List.sort()     
print 'data matched betwen matrix     ', "--- %s minutes run ---" %((time.time() - start_time) / 60)
print 'number of out_put records: ', count
                        
# write final output file
output_handle = open("%s%s" %(path,OutBEDFile), "w") # Opens the file to be written
count = 0
for items in Final_List:
    for each in items:     
        output_handle.write("%s" %each[0])                            
        count +=1                                     
output_handle.close()

# end
print 'files parsed in         ', "--- %s minutes run ---" %((time.time() - start_time) / 60)
 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#######################################################################
# STEP 2:
# Search for junction reads (overlapping the edge of ONSEN and genome) in the unmapped reads.
# The strategy involves selecting unmapped reads having similarity to the 150 bp edges of ONSEN,
# then trimming them with the 30 bp edges of ONSEN, and finally re-map them to genome
#######################################################################

#########
# 1) Recognize those mapping to 150 bp edges, prepare the 5' and 3' edge with reverse compliment, and blast Unmapped reads to 150 bp edges
#########

# Manually prepare a tab-delimited file from the .bed file used to mask ('output_ONSEN_sorted_merged.bed'),
# with coordinates of 5' LTR start and add 150 bp as end; and then coordinates of end in the LTR and deleting 150 bp as start.
# Use python script to generate the final .bed file
# FILE: 1_generates bed4 from chr_start_end list for LTRs.py
# Input:	input_150edges_ONSEN.txt
# Output:	output_150edges_ONSEN.txt
# Input:	input_30edges_ONSEN.txt
# Output:	output_30edges_ONSEN.txt
# Perform the same for a file with 30 bp edges, for posterior trimming

# Get fasta
$bedtools getfasta -fi ~/PATH/TAIR10_chr_all.fas \
-bed ~/PATH/output_150edges_ONSEN.txt \
-fo ~/PATH/fasta_150edges_ONSEN.fa -name &

$bedtools getfasta -fi ~/PATH/TAIR10_chr_all.fas \
-bed ~/PATH/output_30edges_ONSEN.txt \
-fo ~/PATH/fasta_30edges_ONSEN.fa -name &


# Use python script to recover the reverse comliment and
# add it to the .fasta sequences.
# This step may be unnecessary but I didnt know the behaviour
# of Trimmomatic, so just in case I introduced this step in the original workflow.
# FILE:	6_input fasta sequences and generate reverse complement_include both in outcome_with SeqRecord.py
# Input:	fasta_150edges_ONSEN.fa
# Output:	fasta_150edges_ONSEN_with_reverse.fasta
# Input:	fasta_30edges_ONSEN.fa
# Output:	fasta_30edges_ONSEN_with_reverse.fasta

# Blast unmapped reads to ONSEN 150 bp edges
$blastn -subject ~/PATH/Reads_SAMPLE_GenomeUnmapped.fasta -query ~/PATH/fasta_150edges_ONSEN_with_reverse.fasta \
-max_target_seqs 10000 -outfmt "6" -out ~/PATH/SAMPLE_Blast_Unmapped_to_ONSEN150edges.blast &

# Manually check the blast outcome; see if it is reporting partial blasting and look
# for incomplete alignements (less than <read lenght> and more or equal than 30 bp)
# Manually recover reads after filter larger to smaller,
# 1st s.start, 2nd aligment lenght, and 3rd %identity
# As resequencing was 75 bp, use reads from 74bp aligment lenght with 100% identity
# (perfect blasting of the whole read lenght is of no use, becuase it does not represent junctions)
# Final file:	 SAMPLE_Blast_Unmapped_to_ONSEN150edges.txt
# Use custome-made python script to retrive the ID names from the blast report, and then recover the reads from a fastq file.
# FiLE:	7bis_parse 1xblast file_select reads from fastq file_output fastq.py
# Script:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import time
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq

start_time = time.time() # To tell de time it take to process

##Files and path
InputBlast = "SAMPLE_Blast_Unmapped_to_ONSEN150edges.txt" # This is the outcome from the blast analysis
InputFastQ = "Reads_SAMPLE_GenomeUnmapped.fastq" # This is fastq
OutFastQ   = 'SAMPLE_Blasted_to_150edges_ONSEN_fromUnmapped.fastq' # This is the text file with the results of the final parsing, showing best hits
path       = "~/PATH/"
d_type     = 'a48' # Max lenght of read names in number of characters

# define CountLines fx to count the lines a the file
def CountLines (filename):
    with open ("%s%s" %(path,filename), 'r') as myfile:	  
        count = sum(1 for line in myfile)
    return count

# defines fx ParseBlast
def ParseBlast (input_file, number_lines):

    Blast_Matrix = np.array(range(number_lines), dtype = d_type) # Generates a matrix with strings of at most d_type characters, to store names of reads
    i = 0 # Counts position in matrix
    for record in input_file: 
    # Historical contingency. As base for this function I used my previous functions opening SAM files;
    # this would disregard SAM header. Thus this is not necesarry here,
    # but it was used in the original workflow.
        if record.startswith("@"): continue
        else:                      
            line     = record.strip().split("\t")
            name_sam = line[1]               		
            evalue   = line [10]
            
            if evalue > 1:
                
                Blast_Matrix[i] = name_sam  
                i = i + 1
    Blast_Matrix.resize(i, 1)
    Blast_Matrix = np.unique(Blast_Matrix) # Removes duplicated reads and sorts
    return Blast_Matrix

################################################################################

############
### MAIN ###
############

# parse the 1st File
counts_1 = CountLines (InputBlast); print 'fist file count: ', counts_1 
input_handle_1 = open("%s%s" %(path,InputBlast), "r") # Opens the file
Data_File1 = ParseBlast (input_handle_1, counts_1) # Applies ParseBlast fx

input_handle_1.close()
print "fist blast file parsed ", "--- %s minutes run ---" %((time.time() - start_time) / 60)
                    
# Match data
input_handle = open("%s%s" %(path,InputFastQ), "rU")    # Opens .fastq file
output_handle = open("%s%s" %(path,OutFastQ), "w") 	# Opens the file to be written

count = 0
for record in SeqIO.parse(input_handle, "fastq") :
 
    if record.id in Data_File1:	# Search if the ID in .fastq file exist in blast report
        count=count+1; print count,
        SeqIO.write(record, output_handle, "fastq")

input_handle.close()
output_handle.close()  
      
print "Total reads match: ", count, " in time: --- %s minutes run ---" %((time.time() - start_time) / 60)

# end
print 'files parsed in         ', "--- %s minutes run ---" %((time.time() - start_time) / 60)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#########
# 2) Trim TEs sequences from the selected reads with the 30 bp edges of ONSEN file prepared before,
# and re-map "sof-clipped" to genome (acounts for TDS)
#########

# Run Trimmomatic with the selected reads in .fastq
$java -jar ~/PATH/trimmomatic-testing.jar SE -threads 4 \
~/PATH/SAMPLE_Blasted_to_150edges_ONSEN_fromUnmapped.fastq \
~/PATH/SAMPLE_Blasted_to_150edges_ONSEN_fromUnmapped_trimmoEDGES.fastq \
ILLUMINACLIP:~/PATH/fasta_30edges_ONSEN_with_reverse.fasta:2:10:5:1 &

# Recover .fasta from .fastq
$awk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' SAMPLE_Blasted_to_150edges_ONSEN_fromUnmapped_trimmoEDGES.fastq > SAMPLE_Blasted_to_150edges_ONSEN_fromUnmapped_trimmoEDGES.fasta &

# Remap with bowtie but in not in --end-to-end mode. Local mode allows for free tails (soft clipped at the end of the reads).
$bowtie2 -p 3 --quiet --local --very-sensitive-local --score-min L,5,0 --np 0 \
-x ~/PATH/TAIR10_chr_all_maskedONSEN \
-U ~/PATH/SAMPLE_Blasted_to_150edges_ONSEN_fromUnmapped_trimmoEDGES.fastq \
-S ~/PATH/SAMPLE_Blasted_to_150edges_ONSEN_fromUnmapped_trimmoEDGES_bowtie2_local_TAIR10_chr_all_maskedONSEN.sam &

# Convert sam to bam, sort
$samtools view -bS ~/PATH/SAMPLE_Blasted_to_150edges_ONSEN_fromUnmapped_trimmoEDGES_bowtie2_local_TAIR10_chr_all_maskedONSEN.sam \
> ~/PATH/SAMPLE_Blasted_to_150edges_ONSEN_fromUnmapped_trimmoEDGES_bowtie2_local_TAIR10_chr_all_maskedONSEN.bam &

$samtools sort ~/PATH/SAMPLE_Blasted_to_150edges_ONSEN_fromUnmapped_trimmoEDGES_bowtie2_local_TAIR10_chr_all_maskedONSEN.bam \
~/PATH/SAMPLE_Blasted_to_150edges_ONSEN_fromUnmapped_trimmoEDGES_bowtie2_local_TAIR10_chr_all_maskedONSEN_sorted &

# Recover mapped reads
$samtools view -F 4 -h ~/PATH/SAMPLE_Blasted_to_150edges_ONSEN_fromUnmapped_trimmoEDGES_bowtie2_local_TAIR10_chr_all_maskedONSEN_sorted.bam \
> ~/PATH/SAMPLE_MappedGenome_from_Unmapped-blasted_to_150edges_ONSEN.sam &

$samtools view -bS -h ~/PATH/SAMPLE_MappedGenome_from_Unmapped-blasted_to_150edges_ONSEN.sam \
> ~/PATH/SAMPLE_MappedGenome_from_Unmapped-blasted_to_150edges_ONSEN.bam &

#######################################################################
# Step 3:
# Althoght previous files can be visualized in a genome browser for manual assesment, I also did some
# intersection for ease of comparisons across many sequenced lines. Intersected 'discordant' PE reads hitting ONSEN,
# with 'junction' reads overlapping the edge of ONSEN-genome.
#######################################################################
# (i) Intersect outcomes
#
# Use python script to extend 450 bp the reads in .bed file
# 'Out_Reads_SAMPLE_GenomethatoneUnmapped_vs_GenomeUnmapped-but-ONSENmapped.bed'.
# But add to those in the + strand and decrease to those in the  - strand
# FILE:	9_input bed extend strands and generates new bed.py
#
# Script:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import time
import numpy as np

start_time = time.time() # To tell de time it take to process

# Files and path
InBEDFile  = "Out_Reads_SAMPLE_GenomethatoneUnmapped_vs_GenomeUnmapped-but-ONSENmapped.bed"
OutBEDFile = "Out_EXTENDED_Reads_SAMPLE_GenomethatoneUnmapped_vs_GenomeUnmapped-but-ONSENmapped.bed"
path       = "~/PATH/"

# define CountLines fx to count the lines a the file
def CountLines (filename):
    with open ("%s%s" %(path,filename), 'r') as myfile: # Open the input file
        count = sum(1 for line in myfile)
    return count   
################################################################################           

############
### MAIN ###
############

# count lines in files
counts_bed = CountLines (InBEDFile)
print 'Number of lines in BED file: ', counts_bed, "--- %s minutes run ---" %((time.time() - start_time) / 60)

# input data from .bed in matrix
input_bed_handle = open("%s%s" %(path,InBEDFile), "r") # Opens input file
output_handle = open("%s%s" %(path,OutBEDFile), "w") # Opens the file to be written

count = 0
for line in input_bed_handle:

    Chr,start,end,name,dot,strand = line.strip().split("\t") # Unpack data from .bed6
    count += 1
    print start,end,strand, "\t",
    start = int(start)
    end   = int(end)
    
    if strand == '+':
        end = end + 450
    if strand == '-':
        start = start - 450
    print start,end
    out_data = [Chr, '\t', str(start), '\t', str(end), '\t', name, '\t', dot, '\t', strand, '\n']
    for data in out_data:
        output_handle.write("%s" %data) # Write final output file
 
input_bed_handle.close()                                   
output_handle.close()
    
print 'number of out_put records: ', count
    
# end
print 'file parsed in         ', "--- %s minutes run ---" %((time.time()-start_time)/60)   

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Intersect 'SAMPLE_MappedGenome_from_Unmapped-blasted_to_150edges_ONSEN.bam' with
# the 'Out_EXTENDED_Reads_SAMPLE_GenomethatoneUnmapped_vs_GenomeUnmapped-but-ONSENmapped.bed
$bedtools intersect -wa -abam ~/PATH/SAMPLE_MappedGenome_from_Unmapped-blasted_to_150edges_ONSEN.bam \
-b ~/PATH/Out_EXTENDED_Reads_SAMPLE_GenomethatoneUnmapped_vs_GenomeUnmapped-but-ONSENmapped.bed \
> ~/PATH/Spanning_SAMPLE_GenomethatoneUnmapped_intersect_junctions.bam &










