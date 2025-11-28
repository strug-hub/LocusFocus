import pandas as pd
from flask import current_app as app
from pymongo import MongoClient


from app import mongo
from app.utils.variants import get_variants_by_region
from app.utils import parse_region_text
from app.utils.gencode import collapsed_genes_df_hg19, collapsed_genes_df_hg38
from app.utils.errors import InvalidUsage

client = None  # type: ignore

with app.app_context():
    client: MongoClient = mongo.cx  # type: ignore


# This is the main function to extract the data for a tissue and gene_id:
def get_gtex(version, tissue, gene_id):
    if version.upper() == "V8":
        db = client.GTEx_V8
    elif version.upper() == "V10":
        db = client.GTEx_V10
    elif version.upper() == "V7":
        raise InvalidUsage(
            "Cannot standardize SNPs to hg19; GTEx V7 is no longer available."
        )
    collapsed_genes_df = collapsed_genes_df_hg38

    tissue = tissue.replace(" ", "_")
    # gene_id = gene_id.upper()
    ensg_name = ""
    if tissue not in db.list_collection_names():
        raise InvalidUsage(f"Tissue {tissue} not found", status_code=410)
    collection = db[tissue]
    if gene_id.startswith("ENSG"):
        i = list(collapsed_genes_df["ENSG_name"]).index(gene_id)
        ensg_name = list(collapsed_genes_df["ENSG_name"])[i]
    elif gene_id in list(collapsed_genes_df["name"]):
        i = list(collapsed_genes_df["name"]).index(gene_id)
        ensg_name = list(collapsed_genes_df["ENSG_name"])[i]
    else:
        raise InvalidUsage(f"Gene name {gene_id} not found", status_code=410)
    ensg_name = ensg_name.rsplit(".", 1)[0]  # remove version
    results = list(collection.find({"gene_id": {"$regex": f"^{ensg_name}.*"}}))
    response = []
    try:
        response = results[0]["eqtl_variants"]
    except Exception:
        return pd.DataFrame([{"error": f"No eQTL data for {gene_id} in {tissue}"}])
    results_df = pd.DataFrame(response)
    chrom = int(list(results_df["variant_id"])[0].split("_")[0].replace("X", "23"))
    positions = [int(x.split("_")[1]) for x in list(results_df["variant_id"])]
    variants_df = get_variants_by_region(
        min(positions), max(positions), str(chrom), version.upper()
    )
    x = pd.merge(results_df, variants_df, on="variant_id")
    return x


# Function to merge the GTEx data with a particular snp_list
def get_gtex_data(version, tissue, gene, snp_list, raiseErrors=False) -> pd.DataFrame:
    if version.upper() == "V7":
        raise InvalidUsage(
            "GTEx V7 is no longer available. Please use GTEx V8 or GTEx V10."
        )
    assert version.upper() in ["V8", "V10"]
    build = "hg38"
    gtex_data = []
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
    hugo_gene, ensg_gene = gene_names(gene, build)
    response_df = get_gtex(version.upper(), tissue, gene)
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


# This function simply merges the eqtl_data extracted with the snp_list,
# then returns a list of the eQTL pvalues for snp_list (if available)
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
    """
    Return the number of SNPs that can be found in the GTEx database for the given region.
    """
    assert gtex_version.upper() in ["V8", "V10"]
    # Ensure valid region:
    chrom, startbp, endbp = parse_region_text(regiontxt, build)
    chrom = str(chrom).replace("23", "X")

    # Lookup variants in GTEx db
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
    # Given either ENSG gene name or HUGO gene name, returns both HUGO and ENSG names
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
