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
[**Environment and Usage**](#installation-and-usage)  
[**Functions**](#commands)  
[**Quality check**](#Quality check)  
[**Reads clean**](#Reads clean)  
[**Genome index**](#Genome index) 
[**Reads mapping**](#Reads mapping)
[**Reads counts**](#Reads counts)

# Environment and Usage
_**RNS_Sequencing requires Alabama Super Computer**. You can request an new or existing account from [ASC](https://www.asc.edu/hpc/ASA-HPC-Annual-Grant-Request-Form). Highly recommand to email them if you need to contact them._

To run the RNA_Sequencing utility, set up your ASC account and run from ASC 

```  
$ ssh your_asc_account@@dmc.asc.edu
```  
**Before** you loged in, copy and paste from local to ASC
```  
$ scp your_local_path_of_files your_asc_account@dmc.asc.edu:~/path_of_the_directory_you_want_to_paste_into
```
Once you login, **submit** jobs to ASC to run your scripts
# Functions


## Box Folder

Here is the link to the Box Folder with our data: https://auburn.box.com/s/x17154cah63h3yp4wd14ak6p0fa3t5jz
