# CellRanger code for generation of scRNAseq data using fastaq files 
# For better resuts and time manner, we used cellranger/7.2.0 pipeline for generating the raw gene expression data for each cells for all samples.
# https://www.10xgenomics.com/support
# need to download the cellranger in cluster, anaconda environment, then arrange data files for input and output locations. 
# fastaq sample generation and protocol can see in Main manuscript.
# repeat this code for each sample 
# Below is the real code to run simple pipeline, all depends on the capacity of cluster. It may take time to run pipelines 
# If everything is good then you can expect /Sample_1/outs/raw_feature_bc_matrix files in the format matrix, barcodes and cell IDs can use further downstream analysis. 
# or else you can also use raw_feature_bc_matrix.h5 file directly where all information at one place to use. 

#!/bin/bash
#SBATCH --job-name=cellranger-job
#SBATCH --mail-user=ID@institute.edu   # Email ID for job output 
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20  # Adjust as needed
#SBATCH --mem=125G           # Adjust as needed
#SBATCH --time=48:00:00
#SBATCH --chdir=/path/to/working/directory  # Set the working directory
#SBATCH --error=%x.%j.joberr  # Error file output
#SBATCH --output=%x.%j.jobout  # Job output log

# Load necessary modules
module load cellranger/7.2.0  

# Run Cell Ranger count command
cellranger count --id=Sample_1 \
                 --transcriptome=/path/to/refdata-gex-GRCh38-2020-A \  # Reference genome
                 --fastqs=/path/to/Multiome-GEX \  # FASTQ file directory
                 --sample=Sample_1 \
                 --chemistry=ARC-v1 \  # Verify chemistry version with Cell Ranger updates
                 --include-introns=true
