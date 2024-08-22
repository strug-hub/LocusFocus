"""
Stages for running colocalization.

Expected order:
- create_session
- collect_user_input
- read_gwas_file
- read_secondary_datasets
- report_gtex_data
- get_ld_matrix
- ss_subset_gwas
- ss_colocalization X
- coloc2 X
- save_coloc_results X

"""