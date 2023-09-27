# evvie_crispr_scripts

guide_assignment_9.26.23_w_tar_and_reformat_perturb.sh is a bash script designed to take the output of a 10x cell ranger analysis of a CRISPRi dataset, compute guide assignmnet, produce guide assignment QCs, and create output files for downstream processes. 

The script can be run on a slurm-based system as follows:

sbatch guide_assignment_9.26.23_w_tar_and_reformat_perturb.sh SAMPLE_NAME PATH_to_10x_DIR_ENDING_IN/ MIN_GUIDE_UMI_THRESHOLD Y_N_USE_10x_GUIDE_CALLS Y_N_Testing

For example:
sbatch guide_assignment_9.26.23_w_tar_and_reformat_perturb.sh 230327_Encode_Tap_MOI5_sample7_novaseq_and_qc_seq /oak/stanford/groups/engreitz/Users/ejagoda/230327_Encode_K562_Tap_seq_full_seq/230327_Encode_Tap_MOI5_sample7_novaseq_and_qc_seq/ 2 yes no

The script will do as follows:

1. unzip the 10x cell ranger outputs crispr_analysis.tar.gz and filtered_feature_bc_matrix.tar.gz
2. call the script all_upstream_of_mast_read_edits9.22.23_double_check.R which computes guide assignments and creates the following output files described below
3. 1. **_protospacer_calls_per_cell.csv** - final guide calls per cell formatted in the cellranger format
2. **_guides_per_cell.txt** - CBC x #guides assigned
3. **_per_guide_summary_table.txt** - guide - number cells assigned - umi threshold used
4. **_perturb_status.txt** - binary matrix guide assignments in the form of CBC x guide
5. **_Cells_per_guide.png** - histogram of the number of cells assigned to each guide
6. **Guides_per_cell.png** - histogram of the number of guides assigned to each cell
7. **_umi_threshold_hist_log10.png** - hisogram of number of minimum guide umis needed to assign each guide to a cell
8. **_assinged_vs_not_assinged_umis_cdf.png** - per cell CDF plot of ratio of guide umis per cell that were attributed to an assigned guide vs a non-assigned guide (ambient rna)
9. **_assinged_vs_not_assinged_umis_boxplot.png** - boxplot of ratio of guide umis per cell that were attributed to an assigned guide vs a non-assigned guide (ambient rna)

4. modify the header of the "_perturb_status.txt" to make it useable for the MAST snakemake

Output files from step 2:
