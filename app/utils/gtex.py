
import pandas as pd
from flask import current_app as app
from pymongo import MongoClient


from app import mongo
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
        collapsed_genes_df = collapsed_genes_df_hg38
    elif version.upper() == "V7":
        raise ValueError("Cannot standardize SNPs to hg19; GTEx V7 is no longer available.")

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
    results = list(collection.find({"gene_id": ensg_name}))
    response = []
    try:
        response = results[0]["eqtl_variants"]
    except:
        return pd.DataFrame([{"error": f"No eQTL data for {gene_id} in {tissue}"}])
    results_df = pd.DataFrame(response)
    chrom = int(list(results_df["variant_id"])[0].split("_")[0].replace("X", "23"))
    positions = [int(x.split("_")[1]) for x in list(results_df["variant_id"])]
    variants_query = db.variant_table.find(
        {
            "$and": [
                {"chr": chrom},
                {"variant_pos": {"$gte": min(positions), "$lte": max(positions)}},
            ]
        }
    )
    variants_list = list(variants_query)
    variants_df = pd.DataFrame(variants_list)
    variants_df = variants_df.drop(["_id"], axis=1)
    x = pd.merge(results_df, variants_df, on="variant_id")
    if version.upper() == "V7":
        raise ValueError("GTEx V7 is no longer available.")
    elif version.upper() == "V8":
        x.rename(columns={"rs_id_dbSNP151_GRCh38p7": "rs_id"}, inplace=True)
    return x


# Function to merge the GTEx data with a particular snp_list
def get_gtex_data(version, tissue, gene, snp_list, raiseErrors=False) -> pd.DataFrame:
    build = "hg19"
    if version.upper() == "V8":
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
    #    print(f'Gathering eQTL data for {hugo_gene} ({ensg_gene}) in {tissue}')
    response_df = pd.DataFrame({})
    if version.upper() == "V7":
        raise InvalidUsage("GTEx V7 is no longer available.")
    elif version.upper() == "V8":
        response_df = get_gtex("V8", tissue, gene)
    if "error" not in response_df.columns:
        eqtl = response_df
        if rsids:
            snp_df = pd.DataFrame(snp_list, columns=["rs_id"])
            # idx = pd.Index(list(snp_df['rs_id']))
            idx2 = pd.Index(list(eqtl["rs_id"]))
            # snp_df = snp_df[~idx.duplicated()]
            eqtl = eqtl[~idx2.duplicated()]
            # print('snp_df.shape' + str(snp_df.shape))
            gtex_data = (
                snp_df.reset_index()
                .merge(eqtl, on="rs_id", how="left", sort=False)
                .sort_values("index")
            )
            # print('gtex_data.shape' + str(gtex_data.shape))
            # print(gtex_data)
        else:
            snp_df = pd.DataFrame(snp_list, columns=["variant_id"])
            gtex_data = (
                snp_df.reset_index()
                .merge(eqtl, on="variant_id", how="left", sort=False)
                .sort_values("index")
            )
    else:
        try:
            error_message = list(response_df["error"])[0]
            gtex_data = pd.DataFrame({})
        except:
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


def get_gtex_snp_matches(stdsnplist, regiontxt, build):
    """
    Return the number of SNPs that can be found in the GTEx database for the given region.
    """
    # Ensure valid region:
    chrom, startbp, endbp = parse_region_text(regiontxt, build)
    chrom = str(chrom).replace("23", "X")

    # Load GTEx variant lookup table for region indicated
    db = client.GTEx_V8
    if build.lower() in ["hg19", "grch37"]:
        raise InvalidUsage("Cannot use GTEx V7 variant table; GTEx V7 is no longer available.")
    collection = db["variant_table"]
    variants_query = collection.find(
        {
            "$and": [
                {"chr": int(chrom.replace("X", "23"))},
                {"variant_pos": {"$gte": int(startbp), "$lte": int(endbp)}},
            ]
        }
    )
    variants_list = list(variants_query)
    variants_df = pd.DataFrame(variants_list)
    variants_df = variants_df.drop(["_id"], axis=1)
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
