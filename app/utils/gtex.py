import pandas as pd
from flask import current_app

from app.utils.variants import get_variants_by_region
from app.utils.gencode import collapsed_genes_df_hg19, collapsed_genes_df_hg38
from app.utils.errors import InvalidUsage


def get_gtex(version, tissue, gene_id):
    """Fetch the merged eQTL + variant DataFrame for a tissue/gene pair.

    Returns a DataFrame with all eQTL columns plus rs_id, chr, pos, ref, alt.
    Returns a single-row error DataFrame (column "error") when no eQTL data
    exists for the gene/tissue combination.
    """
    if version.upper() == "V7":
        raise InvalidUsage(
            "Cannot standardize SNPs to hg19; GTEx V7 is no longer available."
        )

    gtex_db = current_app.extensions["gtex_db"]
    version = version.upper()
    tissue = tissue.replace(" ", "_")

    if tissue not in gtex_db.list_tissues(version):
        raise InvalidUsage(f"Tissue {tissue} not found", status_code=410)

    collapsed_genes_df = collapsed_genes_df_hg38

    if gene_id.startswith("ENSG"):
        if gene_id not in list(collapsed_genes_df["ENSG_name"]):
            raise InvalidUsage(f"Gene name {gene_id} not found", status_code=410)
        ensg_name = gene_id
    elif gene_id in list(collapsed_genes_df["name"]):
        i = list(collapsed_genes_df["name"]).index(gene_id)
        ensg_name = list(collapsed_genes_df["ENSG_name"])[i]
    else:
        raise InvalidUsage(f"Gene name {gene_id} not found", status_code=410)

    ensg_id_prefix = ensg_name.rsplit(".", 1)[0]

    result = gtex_db.get_eqtl_data(version, tissue, ensg_id_prefix)
    if result.empty:
        return pd.DataFrame([{"error": f"No eQTL data for {gene_id} in {tissue}"}])
    return result


def get_gtex_data(version, tissue, gene, snp_list, raiseErrors=False) -> pd.DataFrame:
    if version.upper() == "V7":
        raise InvalidUsage(
            "GTEx V7 is no longer available. Please use GTEx V8 or GTEx V10."
        )
    assert version.upper() in ["V8", "V10"]

    rsids = True
    rsid_snps = [x for x in snp_list if x.startswith("rs")]
    b37_snps = [x for x in snp_list if x.endswith("_b37")]
    b38_snps = [x for x in snp_list if x.endswith("_b38")]
    if len(rsid_snps) > 0 and (len(b37_snps) > 0 or len(b38_snps) > 0):
        raise InvalidUsage(
            "There is a mix of rsid and other variant id formats; please use a consistent format"
        )
    elif len(rsid_snps) > 0:
        rsids = True
    elif len(b37_snps) or len(b38_snps) > 0:
        rsids = False
    else:
        raise InvalidUsage(
            "Variant naming format not supported; ensure all are rs ID's are formatted as chrom_pos_ref_alt_b37 eg. 1_205720483_G_A_b37"
        )

    hugo_gene, ensg_gene = gene_names(gene, "hg38")
    response_df = get_gtex(version.upper(), tissue, gene)

    gtex_data = []
    if "error" not in response_df.columns:
        eqtl = response_df
        if rsids:
            snp_df = pd.DataFrame(snp_list, columns=["rs_id"])
            idx2 = pd.Index(list(eqtl["rs_id"]))
            eqtl = eqtl[~idx2.duplicated()]
            gtex_data = (
                snp_df.reset_index()
                .merge(eqtl, on="rs_id", how="left", sort=False)
                .sort_values("index")
            )
        else:
            snp_df = pd.DataFrame(snp_list, columns=["variant_id"])
            gtex_data = (
                snp_df.reset_index()
                .merge(eqtl, on="variant_id", how="left", sort=False)
                .sort_values("index")
            )
    else:
        try:
            gtex_data = pd.DataFrame({})
        except Exception:
            if raiseErrors:
                raise InvalidUsage(
                    "No response for tissue "
                    + tissue.replace("_", " ")
                    + " and gene "
                    + hugo_gene
                    + " ( "
                    + ensg_gene
                    + " )",
                    status_code=410,
                )
    return gtex_data  # type: ignore


def get_gtex_data_pvalues(eqtl_data, snp_list):
    rsids = True
    if snp_list[0].startswith("rs"):
        rsids = True
    elif snp_list[0].endswith("_b37"):
        rsids = False
    elif snp_list[0].endswith("_b38"):
        rsids = False
    else:
        raise InvalidUsage(
            "Variant naming format not supported; ensure all are rs ID's or formatted as chrom_pos_ref_alt_b37 eg. 1_205720483_G_A_b37"
        )
    if rsids:
        gtex_data = pd.merge(
            eqtl_data,
            pd.DataFrame(snp_list, columns=["rs_id"]),
            on="rs_id",
            how="right",
        )
    else:
        gtex_data = pd.merge(
            eqtl_data,
            pd.DataFrame(snp_list, columns=["variant_id"]),
            on="variant_id",
            how="right",
        )
    return list(gtex_data["pval"])


def get_gtex_snp_matches(stdsnplist, regiontxt, build, gtex_version="V10"):
    """Return the number of SNPs that can be found in the GTEx database for the given region."""
    assert gtex_version.upper() in ["V8", "V10"]
    from app.utils import parse_region_text

    chrom, startbp, endbp = parse_region_text(regiontxt, build)
    chrom = str(chrom).replace("23", "X")

    if build.lower() in ["hg19", "grch37"]:
        raise InvalidUsage(
            "Cannot use GTEx V7 variant table; GTEx V7 is no longer available."
        )
    variants_df = get_variants_by_region(
        int(startbp), int(endbp), str(chrom), gtex_version
    )
    gtex_std_snplist = list(variants_df["variant_id"])
    isInGTEx = [x for x in stdsnplist if x in gtex_std_snplist]
    return len(isInGTEx)


def gene_names(genename, build):
    """Given either an ENSG or HUGO gene name, return (HUGO, ENSG) names."""
    ensg_gene = genename
    if build.lower() in ["hg19", "grch37"]:
        collapsed_genes_df = collapsed_genes_df_hg19
    elif build.lower() in ["hg38", "grch38"]:
        collapsed_genes_df = collapsed_genes_df_hg38
    if genename in list(collapsed_genes_df["name"]):
        ensg_gene = collapsed_genes_df["ENSG_name"][
            list(collapsed_genes_df["name"]).index(genename)
        ]
    if genename in list(collapsed_genes_df["ENSG_name"]):
        genename = collapsed_genes_df["name"][
            list(collapsed_genes_df["ENSG_name"]).index(genename)
        ]
    return genename, ensg_gene
