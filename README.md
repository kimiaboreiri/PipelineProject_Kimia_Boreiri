# PipelineProject_Kimia_Boreiri
##Projet Overview



This project is a Bioinformatics pipeline and designed to analyze HCMV (Human Cytomegalovirus) gene expression over time.
The main goal is to -> compare transcriptomes at 2 and 6 days post-infection.



##Project Workflow

The pipline folowing these steps:



1- Download sequencing data -> from SRA (HCMV transcriptomics at 2 and 6 days after infection)

2- Build a transcriptome index -> kallisto

3- Quantify transcript expression (TPM) -> for each CDS in the samples

4- Identify differentially expressed genes between two time points -> slueth R package

5- Filter sequencing reads mapped to HCMV  -> Bowtie2

6- Assemble transcriptomic data -> SPAes

7- Identify closest viral strans -> BLAST


##Required Tools

The following tools are used in this pipeline:


"Python" 


"Kallisto" for transcript quantification


"sleuth" for differential expression analysis


"Bowtie2" for read filtering


"SPAdes" fo sequence assembly


"BLAST" for starin identification

#Installation
'''bash
pip intall biooython pandas numpy
sudo apt intall kallisto bowtie2 spades blast+ sra-toolkit


 
##Download the SRA Datasets
using wget or prefetch(sra-toolkit) to download all raw sequencing file:
For data test I used head to extract a subset of reads



wget https://www.ncbi.nlm.nih.gov/sra/SRX2896360
wget https://www.ncbi.nlm.nih.gov/sra/SRX2896363
wget https://www.ncbi.nlm.nih.gov/sra/SRX2896374
wget https://www.ncbi.nlm.nih.gov/sra/SRX2896375
