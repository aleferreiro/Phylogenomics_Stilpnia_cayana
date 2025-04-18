#########################################################################################################
######################################### Ipyrad Assembly Workflow ######################################

# Activate conda environment with ipyrad
conda env list                     # List all available environments
conda activate /home/gsc/miniconda2/envs/ipyradlast
ipyrad --version                  # Confirm that ipyrad is correctly loaded

# Step 0: Create parameter file for the new project
ipyrad -n Tangara_2023

# Edit the parameter file, especially [4] which defines the path to raw data
# Create a directory for the project
mkdir Demult_Tangara_2023

####################################################
###### Ingroup 83 samples dataset ##################

###### Step 1 ##################
# Project initialization and demultiplexing (if needed)
ipyrad -p params-Tangara_2023.txt -s 1

# Branch the project to exclude samples with low raw read counts
ipyrad -p params-Tangara_2023.txt -b Tangara_2023_83-II - UFG3804 UFG5634 UFG5660 UFG5667 UFG5690 UFG5692 UFG5764 UFG5767 UFG5770 UFG5871 UFG5930
###### Step 2 ##################
# Quality filtering of reads
ipyrad -p params-Tangara_2023_83-II.txt -s 2

###### Step 3 ##################
# Clustering and mapping within individuals
ipyrad -p params-Tangara_2023_83-II.txt -s 3

###### Step 4 ##################
# Estimate sequencing error rate and heterozygosity within clusters
# This step helps distinguish sequencing errors from true variation.
ipyrad -p params-Tangara_2023_83-II.txt -s 4

###### Step 5 ##################
# Generate consensus sequences per cluster
ipyrad -p params-Tangara_2023_83-II.txt -s 5

###### Step 6 ##################
# Inter-individual clustering and multiple sequence alignment
# Before step 7, update the parameter file:
# Set 'min_samples_locus' to 65 to retain loci present in ~80% of individuals.

# Run step 6 only after adjusting this parameter.
ipyrad -p params-Tangara_2023_83-II.txt -s 6


###### Step 7 ##################
# Output assembled files
ipyrad -p params-Tangara_2023_83-II.txt -s 7


####################################################
###### Ingroup 80 samples dataset ##################

# We perform a branching from the above project and remove 3 samples
ipyrad -p params-Tangara_2023_83-II.txt -b Tangara_2023_80-II.txt - UFG4028 UFG5668 UFG5753

# Check param file, correct 'min_samples_locus' to retain loci present in ~80% of individuals,
# and run from step 4 for this new branch
ipyrad -p params-Tangara_2023_80.txt -s 4567