#!/bin/bash
#SBATCH --job-name=py_get_perturb
#SBATCH --partition=engreitz
#SBATCH --time=24:00:00
#SBATCH --mail-user=ejagoda@broadinstitute.org
#SBATCH --mail-type=ALL
#SBATCH --mem=20gb



module load R/4.1.2

#Rscript /oak/stanford/groups/engreitz/Users/ejagoda/updated_crispri_scripts/all_upstream_of_mast_read_just_guide.R fresh_moi5A /oak/stanford/groups/engreitz/Users/ejagoda/230901_WTC11_ENCODE_DC_TAP_MOI5_FreshvsFrozen/Outputs_cellranger/Fresh-MOI5A/ 2 yes yes

#1 is the sample name
#2 is the 10x dir
#3 is the min treshold for guide calls
#4 yes/no use 10x calls
#5 yes/no testing flag

#add here the untar thing

if [ -d ${2}"/count" ]; then
	tar -C ${2} -xvf ${2}/count/crispr_analysis.tar.gz
	tar -C ${2} -xvf ${2}/count/filtered_feature_bc_matrix.tar.gz
else
	tar -C ${2} -xvf ${2}/crispr_analysis.tar.gz
        tar -C ${2} -xvf ${2}/filtered_feature_bc_matrix.tar.gz
fi	


Rscript /oak/stanford/groups/engreitz/Users/ejagoda/updated_crispri_scripts/all_upstream_of_mast_read_edits9.22.23_double_check.R $1 $2 $3 $4 $5

#sed -i 's/\./-/g' ${2}/outputs/${1}_perturb_status.txt

perturb_file=${2}/outputs/${1}calls_based_on_10x_${4}min_thresh_${3}_perturb_status.txt

sed -i 's/\./-/g' $perturb_file
