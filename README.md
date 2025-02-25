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


#After downloding all the input fikes and the unzip them I got gbff file by accession number and write a code to make a .log file to store the result pr rerror in there
log_file = os.path.join(base_dir, "PipelineProject.log")

#after getting all the required files then everything is ready to start the extract the CDS part from sequnecs
#previously I define dynamic path for input files which is gbff file then I made the a fynction to open the file get protein id and count CDS 

then show a line about the number of CDS in HCMV genome in log file
####step3
##in step 3 we want to appply kallisto tCDSo quantify the abundances of transcripts from each CDS and make a table in log fike to show the minimum , median , mean , ans max for each results in .tsv table as kallisto output file
for doin kallisto first of all need to make kallisto index 
then kalisto code which is 
          kallisto index -i index.idx refrence_transcriptome.fa.gz file

##for easy following i make a dictionary and define the samples which are related to 2 days after infection and 2 other samples as 6 days after infection

after that based on the kallisto workflow we have to run kallisto quant for each sampe for quantification
code for kallisto quant:
(just befre running should be definet code to move through all the samples)
kallisto quant -i {index_file} -o {output_sample_dir} -b 100 -t 4 {sample_file1} {sample_file2}


##Nest step is doingsome statistics by R
##Sleuth  step4
sleuth is a R package that uses the kallisto outputto do statistical test to check the expression
for runnig R script in the python 
first need t define the path for input and output the need to define the path for R script and  ask the system to apply Rscript for all the input samples in R environment 
in this case i ask to make a table in the log file to show details about the transcript(FDR < 0.05)

##Bowtie2





###step5 
in this step I use SPades to for assembly ans analysis of sequencing data.
Spades is bioinfrmatics tools that revieve the fastq 1 and fastq 2 and do assembly







