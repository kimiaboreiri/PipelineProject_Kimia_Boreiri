import os
from Bio import Entrez, SeqIO
import subprocess
import pandas as pd

# Set Entrez email for NCBI requests
Entrez.email = "kboreiri@luc.edu"

#since we do not want hard code for this pipleline project I use these line 
#to make the path dynamically
#Define the base directory dynamically
base_dir = os.path.expanduser("~/PipelineProject_Kimia_Boreiri/output") 
os.makedirs(base_dir, exist_ok=True)  # make sure the directory exists
os.chdir(base_dir)  # Move into the directory

# Define genome accession number and paths dynamically
accession = "NC_006273.2"
gbff_file = os.path.join(base_dir, f"{accession}.gbff")
cds_fasta = os.path.join(base_dir, f"{accession}_CDS.fasta")
log_file = os.path.join(base_dir, "PipelineProject.log")  #this file same some information 
                                                          #to show a summary about the results

## Step 2: Fetch the HCMV genome
##in this step Fetch the HCMV genome from NCBI by accession number
def fetch_hcmv_genome(output_gbff, accession, log_file):
    #if the file aready exist, just log the info and skip downloading
    #this below line could be change based on the writer desire I preferd skip downloading if exist but show me the log info everytime
    if os.path.exists(output_gbff):
        with open(log_file, "a") as log:
            log.write(f"{accession} genome already exists.\n")  
        return  # Skip download
    #if it doesn't exist, download the file from NCBI by efetch
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
    with open(output_gbff, "w") as gbff_file:
        gbff_file.write(handle.read())
    handle.close()

    
## Then Extract CDS sequences
#to get CDS sequences define input file, output file and log file to save some info
def extract_cds(input_gbff, output_fasta, log_file):
    ##Extract CDS sequences from GenBank file and save to a FASTA file.
    cds_count = 0   #initialize with zero

    ##by these line parse the genebank file and extract the CDS seq
    with open(output_fasta, "w") as fasta_out:
        for record in SeqIO.parse(input_gbff, "genbank"):
            for feature in record.features:
                if feature.type == "CDS":
                    cds_count += 1
                    protein_id = feature.qualifiers.get("protein_id", [f"CDS{cds_count}"])[0]
                    sequence = feature.location.extract(record).seq
                    fasta_out.write(f">{protein_id}\n{sequence}\n")

    #  for this part I save the extracted CDS count
    with open(log_file, "a") as log:
        log.write(f"The HCMV genome {accession} has {cds_count} CDS.\n")

        
# Run the functions

fetch_hcmv_genome(gbff_file, accession, log_file)
extract_cds(gbff_file, cds_fasta, log_file)


##step3
#now i wanna apply kallisto for each sample to quantify the TMP of each CDS in each transcriptome

import os
import subprocess
##make kallisto index
# Define dynamic paths
index_file = os.path.join(base_dir, "NC_006273.2.idx")  
cds_fasta = os.path.join(base_dir, f"{accession}_CDS.fasta")  # Ensure the correct CDS FASTA file path



    # use kallisto index to align kalisto
cmd = f"kallisto index -i {index_file} {cds_fasta}"
result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    # Check if kallisto executed successfully (just in case t make sure kalisto successfuly create)
if result.returncode == 0:
    print("Kallisto index created successfully.")
    with open(log_file, "a") as log:
        log.write("Kallisto index successfully built.\n")
    


##Runnig kallisto quantification on each sample
# Define dynamic paths
fastq_dir = os.path.expanduser("~/project/fastq")  
kallisto_output_dir = os.path.expanduser("~/PipelineProject_Kimia_Boreiri/output") 
os.makedirs(kallisto_output_dir, exist_ok=True)  # make sure output directory exists

# Define input files dynamically
index_file = os.path.join(base_dir, "NC_006273.2.idx")  # Kallisto index file
log_file = os.path.join(base_dir, "PipelineProject.log")  # Log file

# Define samples as a dictionary with conditions
samples = {
    "SRR5660030": "2dpi",
    "SRR5660033": "6dpi",
    "SRR5660044": "2dpi",
    "SRR5660045": "6dpi",
}

# Run Kallisto quant for each sample
for sample, condition in samples.items():
    sample_file1 = os.path.join(fastq_dir, f"{sample}_1.fastq")
    sample_file2 = os.path.join(fastq_dir, f"{sample}_2.fastq")
    output_sample_dir = os.path.join(kallisto_output_dir, sample)

    print(f"Running Kallisto quant for {sample} (Condition: {condition})...")

    # kalisto code 
    cmd = f"kallisto quant -i {index_file} -o {output_sample_dir} -b 100 -t 4 {sample_file1} {sample_file2}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    
##then make tpm and other info ina log file
# Write TPM statistics to the log file
with open(log_file, "a") as log:
    log.write("\nsample\tcondition\tmin_tpm\tmed_tpm\tmean_tpm\tmax_tpm\n")

    for sample, condition in samples.items(): #make table and add info corresponding to each sample and condition
        abundance_file = os.path.join(kallisto_output_dir, sample, "abundance.tsv") #after runnig this line we get tsv file that help for sleuth

        if os.path.exists(abundance_file):
            df = pd.read_csv(abundance_file, sep="\t")

            min_tpm = df["tpm"].min()
            med_tpm = df["tpm"].median()
            mean_tpm = df["tpm"].mean()
            max_tpm = df["tpm"].max()

            log.write(f"{sample}\t{condition}\t{min_tpm:.2f}\t{med_tpm:.2f}\t{mean_tpm:.2f}\t{max_tpm:.2f}\n")

            print(f"Processed TPM statistics for {sample}")
        
##this step is optianl but helps to have better viion about the result
##lately I use it for Bowtie2 process
# Create sample_table.txt with sample and condition mapping
sample_table_file = os.path.join(kallisto_output_dir, "sample_table.txt")
with open(sample_table_file, "w") as sample_table:
    sample_table.write("sample\tcondition\n")
    for sample, condition in samples.items():
        sample_table.write(f"{sample}\t{condition}\n")

print(f"Sample table saved to {sample_table_file}")

import os
import subprocess

# Define base directory dynamically
base_dir = os.path.expanduser("~/PipelineProject_Kimia_Boreiri/output/kallisto_output")

# Define paths
##since sleuth is a R package I make a nano.R script and write code about the sleuth in R 
r_script_path = os.path.join(base_dir, "sleuth.R")  # Path to the R script
log_file = os.path.join(base_dir, "PipelineProject.log")  # Log file

# Run Sleuth R script (this line helps to make sure sleuth is working)
print("Running differential expression analysis using Sleuth...")

cmd = f"Rscript {r_script_path}"  # Command to run R script
result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

###to amke sure sleuth runnig
if result.returncode == 0:
    print("Sleuth analysis completed successfully.")
    with open(log_file, "a") as log:
        log.write("\ntarget_id\ttest_stat\tpval\tqval\n")
else:
    print(f"Error running Sleuth:\n{result.stderr}")  # Print full error message
    with open(log_file, "a") as log:
        log.write(f"\nError running Sleuth:\n{result.stderr}\n")



import os
import subprocess
##in this step i wanna apply bowtie 
# Define paths dynamically
base_dir = os.path.expanduser("~/PipelineProject_Kimia_Boreiri/output")
fastq_dir = os.path.expanduser("~/project/fastq")  # Path to raw fastq files
bowtie2_output_dir = os.path.join(base_dir, "bowtie2_output")
filtered_fastq_dir = os.path.join(base_dir, "filtered_fastq")  
os.makedirs(bowtie2_output_dir, exist_ok=True)
os.makedirs(filtered_fastq_dir, exist_ok=True)

# Define genome index paths
#fisrt of all amke index for bowtie2
accession = "NC_006273.2"
genome_fasta = os.path.join(base_dir, f"{accession}_CDS.fasta")  # The genome FASTA file
bowtie2_index = os.path.join(base_dir, "HCMV_bowtie2_index")  # Prefix for index files
log_file = os.path.join(base_dir, "PipelineProject.log")  # Log file

# Define samples and their corresponding donors
samples = {
    "SRR5660030": ("Donor 1", "2dpi"),
    "SRR5660033": ("Donor 1", "6dpi"),
    "SRR5660044": ("Donor 3", "2dpi"),
    "SRR5660045": ("Donor 3", "6dpi"),
}

## Build the Bowtie2 Index**
if not os.path.exists(f"{bowtie2_index}.1.bt2"):  # check to make sure index is exist for bowtie
    print("Building Bowtie2 index...")
    cmd = f"bowtie2-build {genome_fasta} {bowtie2_index}" #after making sure bowtie exist then run bowtie2
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    # Log index creation (this step is optioanl to make sure bowties runnig without problem)
    with open(log_file, 'a') as log:
        if result.returncode == 0:
            log.write("\n## Bowtie2 index successfully created ##\n")
        else:
            log.write("\n## Error creating Bowtie2 index ##\n")
            log.write(result.stderr + "\n")

## Map Reads Using Bowtie2 and Filter Mapped Reads
with open(log_file, 'a') as log:
    log.write("\n## Bowtie2 read mapping result ##\n")
#iterate over thesample files to make fastq file as bowtie2 output 
##then use this fastq file from bowtie output to run spades
#because spades need fastq file as input
for sample, (donor, timepoint) in samples.items():
    fastq1 = os.path.join(fastq_dir, f"{sample}_1.fastq")
    fastq2 = os.path.join(fastq_dir, f"{sample}_2.fastq")
    output_sam = os.path.join(bowtie2_output_dir, f"{sample}.sam")
    output_bam = os.path.join(bowtie2_output_dir, f"{sample}.bam")
    output_mapped_bam = os.path.join(bowtie2_output_dir, f"{sample}_mapped.bam")

    # Align reads using Bowtie2 for pairend
    print(f'Running Bowtie2 alignment for {sample}...')
    subprocess.run(["bowtie2", "-x", bowtie2_index, "-1", fastq1, "-2", fastq2, "-S", output_sam])

    # Convert SAM to BAM
    print(f"Converting SAM to BAM for {sample}...")
    subprocess.run(["samtools", "view", "-bS", output_sam, "-o", output_bam])

    # Filter mapped reads using Samtools
    print(f"Filtering mapped reads for {sample}...")
    subprocess.run(["samtools", "view", "-bS", "-F", "4", output_sam, "-o", output_mapped_bam])

    ## Convert Mapped BAM to FASTQ for Spades
    filtered_fastq1 = os.path.join(filtered_fastq_dir, f"{sample}_filtered_1.fastq")
    filtered_fastq2 = os.path.join(filtered_fastq_dir, f"{sample}_filtered_2.fastq")

    print(f"Extracting FASTQ reads for {sample}...")
    subprocess.run(["samtools", "sort", "-n", "-o", output_mapped_bam, output_mapped_bam])  # Sort BAM
    subprocess.run(["samtools", "fastq", "-1", filtered_fastq1, "-2", filtered_fastq2, output_mapped_bam])

    ## Count Reads Before and After Filtering
    total_reads = int(subprocess.check_output(["samtools", "view", "-c", output_sam]).strip())
    mapped_reads = int(subprocess.check_output(["samtools", "view", "-c", output_mapped_bam]).strip())

    #show result in log file
    with open(log_file, 'a') as log:
        log.write(f"{donor} ({timepoint}) had {total_reads} read pairs before Bowtie2 filtering and {mapped_reads} read pairs after.\n")

###step6
#start doing spades 
#make sure about the dynamic directory

base_dir = os.path.expanduser("~/PipelineProject_Kimia_Boreiri/output")  # Base output directory
spades_output_dir = os.path.join(base_dir, "spades_output")  # SPAdes output directory
filtered_fastq_dir = os.path.join(base_dir, "filtered_fastq")  # Directory where filtered FASTQ files are stored
log_file = os.path.join(base_dir, "PipelineProject.log")  # Log file path


os.makedirs(spades_output_dir, exist_ok=True)

# Samples grouped by donor
donors = {
    "Donor_1": ["SRR5660030", "SRR5660033"],
    "Donor_3": ["SRR5660044", "SRR5660045"],
}

# Run SPAdes for each donor
with open(log_file, "a") as log:
    log.write("\n## Running SPAdes Assemblies ##\n")

for donor, samples in donors.items():
    donor_fastq_files = []  # Stores FASTQ file paths for each donor

    for sample in samples:
        fq1 = os.path.join(filtered_fastq_dir, f"{sample}_filtered_1.fastq")
        fq2 = os.path.join(filtered_fastq_dir, f"{sample}_filtered_2.fastq")

        # Add FASTQ paths dynamically
        donor_fastq_files.extend(["-1", fq1, "-2", fq2])

    output_dir = os.path.join(spades_output_dir, donor)  # Define SPAdes output directory
    os.makedirs(output_dir, exist_ok=True)

    # Construct SPAdes command dynamically
    spades_cmd = ["spades.py", "--only-assembler", "-k", "77", "-o", output_dir] + donor_fastq_files

    # Run SPAdes and log the command used
    subprocess.run(spades_cmd)
    with open(log_file, "a") as log:
        log.write(f"\nSPAdes command for {donor}: {' '.join(spades_cmd)}\n")

print("SPAdes processing completed.")