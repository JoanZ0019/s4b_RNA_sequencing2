# Analysis of Gene Expression Differentiation Following Environmental Stressors

Next generation sequencing has provided opportunities to understand the molecular and genomic consequences of climate change induced environmental stressors on organisms. Specifically, RNA-seq has enabled us to understand what genes are being up- or down-regulated in response to these stressors. By analyzing gene expression data, we have a much greater understanding of the physiological consequences that have been previously described in many model species. 


## Project Objectives

In this project, **bash** scripts were developed:

1. Check the quality of the reads
2. Clean reads
3. Index the genome
4. Map the reads to a reference genome
5. Get read counts per gene

These analyses will consist of RNA-seq data extracted from the common mouse species, *Mus musculus*, that was exposed to two separate treatment groups. 

## Table of Contents
[**Environment and Usage**](#environment-and-usage)  
[**Steps and functions**](#steps-and-functions)  
  1. Check the quality of the reads using function: [quality_check](#quality-check)
  2. Clean reads using function: [trim_read](#trim-read)  
  3. Check the quality of the trimmed reads using function: [qc_trimmed](#qc_trimmed)
  4. Mapping the genome using function: [mapping](#mapping)
   

# Environment and Usage
_**RNS_Sequencing requires Alabama Super Computer**. You can request an new or existing account from [ASC](https://www.asc.edu/hpc/ASA-HPC-Annual-Grant-Request-Form). Highly recommand to email them if you need to contact them._

- To run the RNA_Sequencing utility, **set up your ASC account** and run from ASC 

```  
$ ssh your_asc_account@dmc.asc.edu
```  
- **Upload** files from local to ASC using `scp`
```  
$ scp your_local_path_of_files your_asc_account@dmc.asc.edu:~/path_of_the_directory_you_want_to_paste_into
```
- Once you login, **submit** jobs to ASC to run your scripts

```  
$ chmod +x my_scripts
$ run_script my_scripts
``` 
- **Download** files from ASC to local (need to exit ASC)
```  
$ scp  your_asc_account@dmc.asc.edu:~/path_of_the_files_in_ASC path_of_dir_in_local
```
# Steps and functions
## 1. quality_check
| [Description](#description) | [Usage](#usage) | [Input Format](#input-format) | [Output Format](#output-format) |
|---|---|---|---|

### Description
`quality_check` load `fastqc` module
### Usage

Running the `quality_check` function requires  and an output path. 

```
$ quality_check /home/asc_account/path_to_your_folder  
``` 
### Input Format
**[Sample input file](sample_files/4040-KH-17.4040-KH-17_0_filtered_R1.fastq.gz)**

|   |Input should be formatted as a compressed fastq file showing as `.fastq.gz`. The fastq file format:|
|---|---|
|line1|Unique read ID.|
|line2|The sequence.|
|line3|A separator, which is simply a plus (+) sign.|
|line4|The base call quality scores.|
### Output Format
**[Sample output file](sample_files/4040-KH-17.4040-KH-14_0_filtered_R1_fastqc.html)**

|   |Output is  a `.html` file with the following columns. It can be downloaded from ASC to local and [viewed](file:///Users/qiongzhang/s4b_RNA_sequencing2/RNASeq_Data/FastQC/4040-KH-17.4040-KH-17_0_filtered_R1_fastqc.html)
in a web browser. |
|---|---|


## Box Folder

Here is the link to the Box Folder with our data: https://auburn.box.com/s/x17154cah63h3yp4wd14ak6p0fa3t5jz
