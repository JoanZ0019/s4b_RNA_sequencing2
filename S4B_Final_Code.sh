#!/bin/bash

# Scripting for Biologists Final
# Collin Modelski, Emily Wilkins, and Qiong Zhang

#-------------------------------------------------#
#RNAseq Data Analysis
"""
The purpose of this script is to prime raw RNA sequences for differential gene expression analysis.
It is recommended to run this sequence through a high performance computer (like ASC).
We are using the Alabama Super Computer, so our loaded modules reflect the syntax of that system. Please adjust accordingly to your own system.

You will need the following in order to run the script:
-raw FASTA files
-Reference Genome (Genome sequence (GRCm39) (FASTA format)
    -- (we used the reference genome from https://www.gencodegenes.org/mouse/)
-Comprehensive gene annotation (GTF or GFF3 format)
    -- can also be found on https://www.gencodegenes.org/mouse/

We are using the data provided through the paper 'Granzyme B prevents aberrant IL-17 production and intestinal pathogenicity in CD4+ T cell'. When unzipping the data file 'RNASeq_Data.zip', there are two subfolders titled 'Case' and 'Control.' These subfolders contain the raw RNA sequences from the two conditions. For each sample there are two reads, R1 and R2 which are saved as .fastq files. The files are labeled with the specific sample ID, specification that it has been filtered, and what read direction it is:  SAMPLEID_filtered_R1.fastq SAMPLEID_filtered_R2.fastq

In this script we will first assess the quality of the reads through FASTQC to identify the adaptor sequence to remove. We will then use TrimGalore! to remove the adaptor sequence and re-run FASTQC to ensure that the data was trimmed properly. Next, we need to index the reference genome (provided by the link above) to prepare for alignment. After indexing, the script aligns the trimmed sequences (fastq files) to the reference genome.


FUTURE STEPS:
After finishing this script, the data is ready to be used to identify differentially expressed genes through programs like DEseq

"""

#-------------Step 1 - FASTQC on raw sequences --------------

"""
The purpose of running FASTQC is to check the quality of the reads as well as to identify the adaptor type used for sequencing. The output of this step will be a .zip and a .html file within the directory 'fastqc_results'
Packages needed for this step -FASTQC
If running on the Alabama Super Computer you need to make sure the 'module load fastqc/0.11.9' is in your script or else it will not run.
"""

function FASTQC_raw {

        #input files: *.fastq.gz in ~/RNASeq_Data/Case and ../Control - this is hardcoded specifically for our data, but this can be changed by the user
        #output files:  *.html to be viewed in a web browser and *.zip
                #This file will be located in ~/RNASeq_Data/fastqc_results as specified in the code. This can be modified by the user if they choose.
        #Packages: FastQC

        #----------------------------------------------------------------------------

        module load fastqc/0.11.9 #loading the fastqc module from ASC

        mkdir fastqc_results #make a new directory for output files
        cd Case #move into the directory with the fastq files
        fastqc *.fastq.gz -o ~/s4b-project/RNASeq_Data/fastqc_results #Use FastQC to perform a quality check of sequences (#Adjust the path of the output to your path on your system)
        cd ../Control #move into directory with the control fastq files
        fastqc *.fastq.gz -o ~/s4b-project/RNASeq_Data/fastqc_results #Use FastQC to perform a quality check of sequences (#Adjust the path of the output to your path on your system)

}

#------------------Step 2 - TrimGalore! to trim the sequences --------------------
"""
The purpose of running TrimGalore! is to trim the adaptor sequences used for sequencing. Using the data from FASTQC we are able to identify the quality of the sequence reads as well as the adaptor sequence used. TrimGalore! is able to automatically identify the adaptor type used and will trim accordingly. It may be necessary to trim additional sequences if the reads are lower.
Packages needed : trimgalore and anaconda
"""
function trim_reads {
        """
        #Input files needed: *.fastq.gz in ~/s4b-project/RNASeq_Data/Case and ../Control - hardcoded here, but this can be customized by the user.
        #Output files needed: *val_1.fq.gz (for R1) or *val_2.fq.gz (for R2) files
                #These files will be located in ~/s4b-project/RNASeq_Data/trimmed_reads
        #Packages: trimgalore
                #trimgalore will automatically detect and cut sequences at Illumina adapters
        """
        #---------------------------------------------------------------------------------------

        module load trimgalore/0.6.6 #loading trimgalore module from ASC

        #making new directories in ~/s4b-project/RNASeq_Data for the trimmed data to be stored
        mkdir trimmed_reads #create new dirextory
        mkdir trimmed_reads/Case
        mkdir trimmed_reads/Control

        cd Case #move into the directory with the fastq files
                #TrimGalore was struggling to pair files, so this had to be hard coded for each pair
        trim_galore --paired --output_dir ~/s4b-project/RNASeq_Data/trimmed_reads/Case 4040-KH-14.4040-KH-14_0_filtered_R1.fastq.gz 4040-KH-14.4040-KH-14_0_filtered_R2.fastq.gz
        trim_galore --paired --output_dir ~/s4b-project/RNASeq_Data/trimmed_reads/Case 4040-KH-16.4040-KH-16_0_filtered_R1.fastq.gz 4040-KH-16.4040-KH-16_0_filtered_R2.fastq.gz
        trim_galore --paired --output_dir ~/s4b-project/RNASeq_Data/trimmed_reads/Case 4040-KH-21.4040-KH-21_0_filtered_R1.fastq.gz 4040-KH-21.4040-KH-21_0_filtered_R2.fastq.gz
        trim_galore --paired --output_dir ~/s4b-project/RNASeq_Data/trimmed_reads/Case 4040-KH-22.4040-KH-22_0_filtered_R1.fastq.gz 4040-KH-22.4040-KH-22_0_filtered_R2.fastq.gz
        trim_galore --paired --output_dir ~/s4b-project/RNASeq_Data/trimmed_reads/Case 4040-KH-23.4040-KH-23_0_filtered_R1.fastq.gz 4040-KH-23.4040-KH-23_0_filtered_R2.fastq.gz
        trim_galore --paired --output_dir ~/s4b-project/RNASeq_Data/trimmed_reads/Case 4040-KH-24.4040-KH-24_0_filtered_R1.fastq.gz 4040-KH-24.4040-KH-24_0_filtered_R2.fastq.gz
        trim_galore --paired --output_dir ~/s4b-project/RNASeq_Data/trimmed_reads/Case 4040-KH-25.4040-KH-25_0_filtered_R1.fastq.gz 4040-KH-25.4040-KH-25_0_filtered_R2.fastq.gz

        cd ../Control #move into directory with the control fastq files
        trim_galore --paired --output_dir ~/s4b-project/RNASeq_Data/trimmed_reads/Control 4040-KH-1.4040-KH-1_0_filtered_R1.fastq.gz 4040-KH-1.4040-KH-1_0_filtered_R2.fastq.gz
        trim_galore --paired --output_dir ~/s4b-project/RNASeq_Data/trimmed_reads/Control 4040-KH-4.4040-KH-4_0_filtered_R1.fastq.gz 4040-KH-4.4040-KH-4_0_filtered_R2.fastq.gz
        trim_galore --paired --output_dir ~/s4b-project/RNASeq_Data/trimmed_reads/Control 4040-KH-5.4040-KH-5_0_filtered_R1.fastq.gz 4040-KH-5.4040-KH-5_0_filtered_R2.fastq.gz
        trim_galore --paired --output_dir ~/s4b-project/RNASeq_Data/trimmed_reads/Control 4040-KH-6.4040-KH-6_0_filtered_R1.fastq.gz 4040-KH-6.4040-KH-6_0_filtered_R2.fastq.gz
        trim_galore --paired --output_dir ~/s4b-project/RNASeq_Data/trimmed_reads/Control 4040-KH-17.4040-KH-17_0_filtered_R1.fastq.gz 4040-KH-17.4040-KH-17_0_filtered_R2.fastq.gz
        trim_galore --paired --output_dir ~/s4b-project/RNASeq_Data/trimmed_reads/Control 4040-KH-18.4040-KH-18_0_filtered_R1.fastq.gz 4040-KH-18.4040-KH-18_0_filtered_R2.fastq.gz
        
        }
 
 #-----------------------Step 3 - Re-run FASTQC on trimmed files -------------------------

"""
The purpose of running FASTQC is to check the quality of the reads as well as to identify the adaptor type used for sequencing. The output of this step will be a .zip and a .html file within the directory 'fastqc_results' which you need to make prior to running this step
Packages needed for this step -FASTQC. If running on the Alabama Super Computer you need to make sure the 'module load fastqc/0.11.9' is in your script or else it will not run.
"""
function FASTQC_trimmed {
        """
        input files: *.fastq.gz in ~/s4b-project/RNASeq_Data/trimmed_reads/Case and ../Control - hardcoded here, but this can be customized by the user.
        output files: *.html to be viewed in a web browser and *.zip
                #This file will be located in ~/s4b-project/RNASeq_Data/TrimmedReads/trimmed-FastQC as specified in the code. This can be modified by the user if they need.
        Packages: FastQC
        """
        #----------------------------------------------------------------------------------------

        module load fastqc/0.11.9 #loading the fastqc module from ASC

        mkdir trimmed_reads/fastqc_results_trimgalore #make a new directory for output files
        
        cd trimmed_reads/Case #move into the directory with the fastq files
        fastqc *.fq.gz -o ~/s4b-project/RNASeq_Data/trimmed_reads/fastqc_results_trimgalore #Use FastQC to perform a quality check of sequences
        cd ../Control #move into directory with the control fastq files
        fastqc *.fq.gz -o ~/s4b-project/RNASeq_Data/trimmed_reads/fastqc_results_trimgalore #Use FastQC to perform a quality check of sequences
}


#-----------------------Step 4 - Index the Reference Genome for Alignment ----------------------------------
"""
The purpose of indexing the reference genome is to allow the aligner to narrow down the potential origin of a sequence within the genome. This process helps to speed up the alignment process. For more information regarding the bowtie2 program, please check out the manual: http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#paired-inputs

"""

function INDEX_genome {

    """
    #Input File: GRCm39.genome.fa - reference genome
    #Output Files: 'mouse.1.bt2' 'mouse.2.bt2' 'mouse.3.bt2' 'mouse.4.bt2' 'mouse.rev.1.bt2' and 'mouse.rev.2.bt2'
    #Packages: bowtie2
    
    The set-up of the indexing script is:
        bowtie2-build [options]* <reference_in> <bt2_base>
    where bowtie2-build is the main command, -f means that the reference_in is a FASTA file, GRCm39.genome.fa is the reference genome. bt2_base is the base name for the output files. We used 'mouse' as the base name since this is the mouse genome. The output files will be used for the next alignment step.
    
    """

    #----------------------------------------------------------------------------------------------------------
    
    source /opt/asn/etc/asn-bash-profiles-special/modules.sh
    module load bowtie2/2.2.9
    
    bowtie2-build -f GRCm39.genome.fa mouse
    
}



#------------------------Step 5 - Align FASTQ Files to the Indexed Genome ------------------------------------------------
"""
The fifth step is to align trimmed FASTQ files to the new indexed genome. This allows for the ability to identify the specific region in the genome that each base pair in the sequence reads originate from. To align our FASTQ files to the indexed genome we generated in step 4, we will also use the bowtie2 program. Specific documentation for this program can be found in the link below. 
http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#paired-inputs 
"""

function Aligning_Reads {
    """
    Alignes FASTQ files to the indexed genome generated in the INDEX_genome function. In
    
    
    """
    source /opt/asn/etc/asn-bash-profiles-special/modules.sh
    module load tophat/2.1.2
    module load bowtie2/2.2.9

    cd trimmed_reads/Control
    #needed to hardcode the file names in order to get it to run
    bowtie2 -x mouse -1 4040-KH-14.4040-KH-14_0_filtered_R1_val_1.fq -2 4040-KH-14.4040-KH-14_0_filtered_R2_val_2.fq  -S aligned_genome_sequences14.sam
    bowtie2 -x mouse -1 4040-KH-1.4040-KH-1_0_filtered_R1_val_1.fq -2 4040-KH-1.4040-KH-1_0_filtered_R2_val_2.fq  -S aligned_genome_sequences1.sam
    bowtie2 -x mouse -1 4040-KH-18.4040-KH-18_0_filtered_R1_val_1.fq -2 4040-KH-18.4040-KH-18_0_filtered_R2_val_2.fq  -S aligned_genome_sequences18.sam
    bowtie2 -x mouse -1 4040-KH-5.4040-KH-5_0_filtered_R1_val_1.fq -2 4040-KH-5.4040-KH-5_0_filtered_R2_val_2.fq  -S aligned_genome_sequences5.sam
    bowtie2 -x mouse -1 4040-KH-6.4040-KH-6_0_filtered_R1_val_1.fq -2 4040-KH-6.4040-KH-6_0_filtered_R2_val_2.fq  -S aligned_genome_sequences6.sam
    
    cp mouse* ../Case
    cd ../Case
    bowtie2 -x mouse -1 4040-KH-16.4040-KH-16_0_filtered_R1_val_1.fq -2 4040-KH-16.4040-KH-16_0_filtered_R2_val_2.fq  -S aligned_genome_sequences16.sam
    bowtie2 -x mouse -1 4040-KH-21.4040-KH-21_0_filtered_R1_val_1.fq -2 4040-KH-21.4040-KH-21_0_filtered_R2_val_2.fq  -S aligned_genome_sequences21.sam
    bowtie2 -x mouse -1 4040-KH-22.4040-KH-22_0_filtered_R1_val_1.fq -2 4040-KH-22.4040-KH-22_0_filtered_R2_val_2.fq  -S aligned_genome_sequences22.sam
    bowtie2 -x mouse -1 4040-KH-23.4040-KH-23_0_filtered_R1_val_1.fq -2 4040-KH-23.4040-KH-23_0_filtered_R2_val_2.fq  -S aligned_genome_sequences23.sam
    bowtie2 -x mouse -1 4040-KH-24.4040-KH-24_0_filtered_R1_val_1.fq -2 4040-KH-24.4040-KH-24_0_filtered_R2_val_2.fq  -S aligned_genome_sequences24.sam
    bowtie2 -x mouse -1 4040-KH-25.4040-KH-25_0_filtered_R1_val_1.fq -2 4040-KH-25.4040-KH-25_0_filtered_R2_val_2.fq  -S aligned_genome_sequences25.sam
}


#---------------------------------Running All of the Functions -----------------------------

function main {
    FASTQC_raw home/aubemw001/s4b-project/RNASeq_Data
    trim_reads home/aubemw001/s4b-project/RNASeq_Data
    FASTQC_trimmed home/aubemw001/s4b-project/RNASeq_Data
    INDEX_genome home/aubemw001/s4b-project/RNASeq_Data
    Aligning_Reads home/aubemw001/s4b-project/RNASeq_Data
    
}

#--------------------------------------------Final Test ------------------------------------------------------

echo "Starting RNA Sequencing"

main

echo "Ready for DeSeq analysis"
