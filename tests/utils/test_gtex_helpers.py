"""
Test suite for helpers in api/utils/gtex.py
"""
from app.utils.apis.gtex import get_independent_eqtls
from app.utils.errors import InvalidUsage
from app.utils.gencode import get_genes_by_location


TISSUES = ['Adipose_Subcutaneous', 'Adipose_Visceral_Omentum', 'Adrenal_Gland', 'Artery_Aorta', 'Artery_Coronary', 'Artery_Tibial', 'Bladder', 'Brain_Amygdala', 'Brain_Anterior_cingulate_cortex_BA24', 'Brain_Caudate_basal_ganglia', 'Brain_Cerebellar_Hemisphere', 'Brain_Cerebellum', 'Brain_Cortex', 'Brain_Frontal_Cortex_BA9', 'Brain_Hippocampus', 'Brain_Hypothalamus', 'Brain_Nucleus_accumbens_basal_ganglia', 'Brain_Putamen_basal_ganglia', 'Brain_Spinal_cord_cervical_c-1', 'Brain_Substantia_nigra', 'Breast_Mammary_Tissue', 'Cells_Cultured_fibroblasts', 'Cells_EBV-transformed_lymphocytes', 'Cervix_Ectocervix', 'Cervix_Endocervix', 'Colon_Sigmoid', 'Colon_Transverse', 'Esophagus_Gastroesophageal_Junction', 'Esophagus_Mucosa', 'Esophagus_Muscularis', 'Fallopian_Tube', 'Heart_Atrial_Appendage', 'Heart_Left_Ventricle', 'Kidney_Cortex', 'Kidney_Medulla', 'Liver', 'Lung', 'Minor_Salivary_Gland', 'Muscle_Skeletal', 'Nerve_Tibial', 'Ovary', 'Pancreas', 'Pituitary', 'Prostate', 'Skin_Not_Sun_Exposed_Suprapubic', 'Skin_Sun_Exposed_Lower_leg', 'Small_Intestine_Terminal_Ileum', 'Spleen', 'Stomach', 'Testis', 'Thyroid', 'Uterus', 'Vagina', 'Whole_Blood']
GENES = ['NUCKS1', 'CDK18', 'RAB29', 'SLC41A1', 'ELK4', 'SLC45A3', 'PM20D1', 'SLC26A9', 'MFSD4A', 'SNORA72', 'RNU6-418P', 'RP4-681L3.2', 'RP11-6B6.3', 'RNU2-19P', 'Metazoa_SRP', 'RAB7B']
SNPS = ['1_205838376_C_T_b38', '1_205838625_T_A_b38', '1_205838885_C_A_b38', '1_205838990_C_G_b38', '1_205839087_C_T_b38', '1_205839580_T_C_b38', '1_205840093_A_G_b38', '1_205840514_G_C_b38', '1_205841502_G_A_b38', '1_205841819_C_T_b38', '1_205841927_C_T_b38', '1_205842700_A_T_b38', '1_205843175_C_T_b38', '1_205843372_G_A_b38', '1_205843486_C_T_b38', '1_205843784_A_G_b38', '1_205844544_T_C_b38', '1_205844549_T_C_b38', '1_205845060_G_A_b38', '1_205845357_G_A_b38', '1_205845567_A_G_b38', '1_205845607_A_G_b38', '1_205846259_A_G_b38', '1_205846692_A_C_b38', '1_205846802_A_T_b38', '1_205846864_C_T_b38', '1_205848190_T_G_b38', '1_205848446_T_C_b38', '1_205848498_T_A_b38', '1_205849051_T_G_b38', '1_205849244_T_C_b38', '1_205849911_C_T_b38', '1_205849976_G_A_b38', '1_205850526_G_C_b38', '1_205850710_G_A_b38', '1_205851774_G_A_b38', '1_205851804_A_T_b38', '1_205852026_G_A_b38', '1_205852680_T_A_b38', '1_205852947_T_C_b38', '1_205852970_C_T_b38', '1_205853070_T_C_b38', '1_205854903_A_G_b38', '1_205854937_T_A_b38', '1_205854985_C_G_b38', '1_205855020_A_G_b38', '1_205855131_A_G_b38', '1_205855507_A_C_b38', '1_205855537_C_T_b38', '1_205855758_A_G_b38']


def test_fetch_gtex_eqtls_v10(flask_app):
    """
    Test case for fetch_gtex_eqtls
    """
    with flask_app.app_context():
        from app.utils.gtex import fetch_gtex_eqtls
        gtex_eqtl_df = fetch_gtex_eqtls("gtex_v10", TISSUES, GENES, SNPS)
        assert len(gtex_eqtl_df) > 0


def test_get_working_query_for_eqtls(flask_app):
    """
    Test that's actually a function to generate a set of parameters that
    actually provide results from the ieqtl endpoint.
    """
    gtex_version = "gtex_v10"
    start, step, maximum = 0, 2_000_000, 160_000_000
    chrom_start, chrom_max = 1, 20
    with flask_app.app_context():
        from app.utils.gtex import fetch_gtex_eqtls

        working_query_found = False

        for chrom in range(chrom_start, chrom_max+1):
            for pos in range(start, maximum, step):
                print(f"chr{chrom}:{pos}-{pos+step}...", end=" ")
                # Get genes
                try:
                    genes = get_genes_by_location("hg38", chrom, pos, pos+step, gencode=True)
                except InvalidUsage:
                    continue
                if len(genes) == 0:
                    continue

                print(f"{len(genes)} genes...", end=" ")

                # Try and get eQTLs
                try:
                    eqtl_response = get_independent_eqtls(dataset_id=gtex_version, gencode_ids=genes, tissue_sites=TISSUES)
                    # gtex_eqtl_df = fetch_gtex_eqtls(gtex_version, TISSUES, genes, snps)
                    assert len(eqtl_response.data) > 0
                    working_query_found = True
                    print(f"Success! Region: chr{chrom}:{pos}-{pos+step}")
                    break
                except AssertionError:
                    print("nuthin")
                    continue
                except Exception as e:
                    print("error", e)
                    continue
            if working_query_found:
                break

        assert working_query_found

        print(f"!!!! Found working query for eQTLs in region !!!!")
        print(f"Genes: {genes}")
        print(f"Variants: {snps}")
