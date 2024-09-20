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
- coloc_simple_sum
- finalize_results
"""
