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

To run the RNA_Sequencing utility, set up your ASC account and run from ASC 

```  
$ ssh your_asc_account@dmc.asc.edu
```  
**Before** you loged in, copy and paste from local to ASC
```  
$ scp your_local_path_of_files your_asc_account@dmc.asc.edu:~/path_of_the_directory_you_want_to_paste_into
```
Once you login, **submit** jobs to ASC to run your scripts

```  
$ chmod +x my_scripts
$ run_script my_scripts
``` 

# Steps and functions
## 1. quality_check
| [Description](#description) | [Usage](#usage) | [Input Format](#input-format) | [Output Format](#output-format) |
|---|---|---|---|

### Description
`quality_check` 
### Usage

Running the `quality_check` function requires an input path, output path, and a minimum size argument. There are also three optional arguments which can be found in the table below.

```
$ snpbinner crosspoints --input PATH --output PATH (--min-length INT | --min-ratio FLOAT) [optional args]  
``` 
### Input Format
**[Sample input file](sample_files/crosspoints_in.fastq.gz)**

|   |Input should be formatted as a tab‑separated value (TSV) file with the following columns.|
|---|---|
|0|The SNP marker ID.|
|1|The position of the marker in base pairs from the start of the chromosome.|
|2+|RIL ID (header) and the called genotype of the RIL at each position.|

### Output Format
**[Sample output file](sample_files/crosspoints_out.tsv)**

|   |Output is formatted as a comma‑separated value (CSV) file with the following columns.|
|---|---|
|0|The RIL ID|
|Odd|Location of a crosspoint. (Empty after the chromosome ends.)|
|Even|Genotype in between the surrounding crosspoints. (Empty after the chromosome ends.)|

## Box Folder

Here is the link to the Box Folder with our data: https://auburn.box.com/s/x17154cah63h3yp4wd14ak6p0fa3t5jz
