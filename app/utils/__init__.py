"""
Common utility functions and classes shared by multiple routes in LocusFocus.
"""
import os
import re
from typing import List, Optional

import pysam
import pandas as pd
import numpy as np
from flask import current_app as app
from werkzeug.utils import secure_filename
from werkzeug.exceptions import RequestEntityTooLarge
from werkzeug.datastructures import ImmutableMultiDict, FileStorage

from app import mongo
from app.utils.errors import InvalidUsage

GENOMIC_WINDOW_LIMIT = 2e6


def get_session_filepath(filename: str) -> os.PathLike:
    """
    Given a desired filename, return the path to the file in the session folder.
    """
    filename = secure_filename(filename)
    with app.app_context():
        return os.path.join(app.config["SESSION_FOLDER"], filename)  # type: ignore


def get_upload_filepath(filename: str) -> os.PathLike:
    """
    Given a desired filename, return the path to the file in the upload folder.
    """
    filename = secure_filename(filename)
    with app.app_context():
        return os.path.join(app.config["UPLOAD_FOLDER"], filename)  # type: ignore


def download_file(file: FileStorage, check_only: bool = False) -> Optional[os.PathLike]:
    """
    Download the given file (from a request.files MultiDict) to a temporary file in the UPLOAD folder.

    If check_only is True, do not download the file.

    Return the path to the saved file.
    Raises an error if the file is too large.
    """
    if file.filename is None:
        # What causes this?
        return None
    filename = secure_filename(file.filename)
    filepath = get_upload_filepath(filename)
    if not check_only:
        file.save(filepath)
    if not os.path.isfile(filepath):
        raise RequestEntityTooLarge(f"File '{filename}' too large")

    return filepath


def get_file_with_ext(filepaths: List[os.PathLike], extensions: List[str]) -> Optional[os.PathLike]:
    """
    Grab the first file in the list that matches the given extensions. This assumes that the 
    files are already uploaded and stored in the upload folder.

    Return None if no such file exists.

    Extensions should not include the period. eg. `["html", "tsv", "txt"]`.
    """

    for filepath in filepaths:
        if os.path.isfile(filepath) and str(filepath).endswith(tuple(extensions)):
            return filepath
    return None


def decompose_variant_list(variant_list):
    """
    Parameters
    ----------
    variantid_list : list
        list of str standardized variants in chr_pos_ref_alt_build format

    Returns
    -------
    A pandas.dataframe with chromosome, pos, reference and alternate alleles columns
    """
    chromlist = [x.split("_")[0] if len(x.split("_")) == 5 else x for x in variant_list]
    chromlist = [int(x) if x not in ["X", "."] else x for x in chromlist]
    poslist = [
        int(x.split("_")[1]) if len(x.split("_")) == 5 else x for x in variant_list
    ]
    reflist = [x.split("_")[2] if len(x.split("_")) == 5 else x for x in variant_list]
    altlist = [x.split("_")[3] if len(x.split("_")) == 5 else x for x in variant_list]
    df = pd.DataFrame(
        {"CHROM": chromlist, "POS": poslist, "REF": reflist, "ALT": altlist,}
    )
    return df


def standardize_snps(variantlist, regiontxt, build):
    """
    Input: Variant names in any of these formats: rsid, chrom_pos_ref_alt, chrom:pos_ref_alt, chrom:pos_ref_alt_b37/b38
    Output: chrom_pos_ref_alt_b37/b38 variant ID format, but looks at GTEx variant lookup table first.
    In the case of multi-allelic variants (e.g. rs2211330(T/A,C)), formats such as 1_205001063_T_A,C_b37 are accepted
    If variant ID format is chr:pos, and the chr:pos has a unique biallelic SNV, then it will be assigned that variant
    """

    if all(x == "." for x in variantlist):
        raise InvalidUsage("No variants provided")

    if np.nan in variantlist:
        raise InvalidUsage(
            "Missing variant IDs detected in row(s): "
            + str([i + 1 for i, x in enumerate(variantlist) if str(x) == "nan"])
        )

    # Ensure valid region:
    chrom, startbp, endbp = parse_region_text(regiontxt, build)
    chrom = str(chrom).replace("23", "X")

    # Load GTEx variant lookup table for region indicated
    db = mongo.cx.GTEx_V7  # type: ignore
    rsid_colname = "rs_id_dbSNP147_GRCh37p13"
    if build.lower() in ["hg38", "grch38"]:
        db = mongo.cx.GTEx_V8  # type: ignore
        rsid_colname = "rs_id_dbSNP151_GRCh38p7"
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

    # Load dbSNP151 SNP names from region indicated
    dbsnp_filepath = ""
    suffix = "b37"
    if build.lower() in ["hg38", "grch38"]:
        suffix = "b38"
        dbsnp_filepath = os.path.join(
            app.config["LF_DATA_FOLDER"], "dbSNP151", "GRCh38p7", "All_20180418.vcf.gz"
        )
    else:
        suffix = "b37"
        dbsnp_filepath = os.path.join(
            app.config["LF_DATA_FOLDER"], "dbSNP151", "GRCh37p13", "All_20180423.vcf.gz"
        )

    # Load dbSNP file
    # delayeddf = delayed(pd.read_csv)(dbsnp_filepath,skiprows=getNumHeaderLines(dbsnp_filepath),sep='\t')
    # dbsnp = dd.from_delayed(delayeddf)
    tbx = pysam.TabixFile(dbsnp_filepath)  # type: ignore
    #    print('Compiling list of known variants in the region from dbSNP151')
    chromcol = []
    poscol = []
    idcol = []
    refcol = []
    altcol = []
    variantid = []  # in chr_pos_ref_alt_build format
    rsids = dict(
        {}
    )  # a multi-allelic variant rsid (key) can be represented in several variantid formats (values)
    for row in tbx.fetch(str(chrom), startbp, endbp):
        rowlist = str(row).split("\t")
        chromi = rowlist[0].replace("chr", "")
        posi = rowlist[1]
        idi = rowlist[2]
        refi = rowlist[3]
        alti = rowlist[4]
        varstr = "_".join([chromi, posi, refi, alti, suffix])
        chromcol.append(chromi)
        poscol.append(posi)
        idcol.append(idi)
        refcol.append(refi)
        altcol.append(alti)
        variantid.append(varstr)
        rsids[idi] = [varstr]
        altalleles = alti.split(
            ","
        )  # could have more than one alt allele (multi-allelic)
        if len(altalleles) > 1:
            varstr = "_".join([chromi, posi, refi, altalleles[0], suffix])
            rsids[idi].append(varstr)
            for i in np.arange(len(altalleles) - 1):
                varstr = "_".join([chromi, posi, refi, altalleles[i + 1], suffix])
                rsids[idi].append(varstr)

    #    print('Cleaning and mapping list of variants')
    variantlist = [
        asnp.split(";")[0].replace(":", "_").replace(".", "") for asnp in variantlist
    ]  # cleaning up the SNP names a bit
    stdvariantlist = []
    for variant in variantlist:
        if variant == "":
            stdvariantlist.append(".")
            continue
        variantstr = variant.replace("chr", "")
        if re.search("^23_", variantstr):
            variantstr = variantstr.replace("23_", "X_", 1)
        if variantstr.startswith("rs"):
            try:
                # Here's the difference from the first function version (we look at GTEx first)
                if variant in list(variants_df[rsid_colname]):
                    stdvar = (
                        variants_df["variant_id"]
                        .loc[variants_df[rsid_colname] == variant]
                        .to_list()[0]
                    )
                    stdvariantlist.append(stdvar)
                else:
                    stdvariantlist.append(rsids[variantstr][0])
            except:
                stdvariantlist.append(".")
        elif re.search(
            r"^\d+_\d+_[A,T,G,C]+_[A,T,C,G]+,*", variantstr.replace("X", "23")
        ):
            strlist = variantstr.split("_")
            strlist = list(filter(None, strlist))  # remove empty strings
            try:
                achr, astart, aend = parse_region_text(
                    strlist[0] + ":" + strlist[1] + "-" + str(int(strlist[1]) + 1),
                    build,
                )
                achr = str(achr).replace("23", "X")
                if achr == str(chrom) and astart >= startbp and astart <= endbp:
                    variantstr = (
                        variantstr.replace("_" + str(suffix), "") + "_" + str(suffix)
                    )
                    if len(variantstr.split("_")) == 5:
                        stdvariantlist.append(variantstr)
                    else:
                        raise InvalidUsage(
                            f"Variant format not recognizable: {variant}. Is it from another coordinate build system?",
                            status_code=410,
                        )
                else:
                    stdvariantlist.append(".")
            except:
                raise InvalidUsage(f"Problem with variant {variant}", status_code=410)
        elif re.search(r"^\d+_\d+_*[A,T,G,C]*", variantstr.replace("X", "23")):
            strlist = variantstr.split("_")
            strlist = list(filter(None, strlist))  # remove empty strings
            try:
                achr, astart, aend = parse_region_text(
                    strlist[0] + ":" + strlist[1] + "-" + str(int(strlist[1]) + 1),
                    build,
                )
                achr = str(achr).replace("23", "X")
                if achr == str(chrom) and astart >= startbp and astart <= endbp:
                    if len(strlist) == 3:
                        aref = strlist[2]
                    else:
                        aref = ""
                    stdvariantlist.append(fetch_snv(achr, astart, aref, build))
                else:
                    stdvariantlist.append(".")
            except:
                raise InvalidUsage(f"Problem with variant {variant}", status_code=410)
        else:
            raise InvalidUsage(
                f"Variant format not recognized: {variant}", status_code=410
            )
    return stdvariantlist


def parse_region_text(regiontext, build):
    if build not in ["hg19", "hg38"]:
        raise InvalidUsage(f"Unrecognized build: {build}", status_code=410)
    regiontext = regiontext.strip().replace(" ", "").replace(",", "").replace("chr", "")
    if not re.search(
        r"^\d+:\d+-\d+$", regiontext.replace("X", "23").replace("x", "23")
    ):
        raise InvalidUsage(
            f"Invalid coordinate format. '{regiontext}' e.g. 1:205,000,000-206,000,000",
            status_code=410,
        )
    chrom = regiontext.split(":")[0].lower().replace("chr", "").upper()
    pos = regiontext.split(":")[1]
    startbp = pos.split("-")[0].replace(",", "")
    endbp = pos.split("-")[1].replace(",", "")
    chromLengths = pd.read_csv(
        os.path.join(app.config["LF_DATA_FOLDER"], build + "_chrom_lengths.txt"),
        sep="\t",
        encoding="utf-8",
    )
    chromLengths.set_index("sequence", inplace=True)
    if chrom in ["X", "x"] or chrom == "23":
        chrom = 23
        maxChromLength = chromLengths.loc["chrX", "length"]
        try:
            startbp = int(startbp)
            endbp = int(endbp)
        except:
            raise InvalidUsage(
                f"Invalid coordinates input: '{regiontext}'", status_code=410
            )
    else:
        try:
            chrom = int(chrom)
            if chrom == 23:
                maxChromLength = chromLengths.loc["chrX", "length"]
            else:
                maxChromLength = chromLengths.loc["chr" + str(chrom), "length"]
            startbp = int(startbp)
            endbp = int(endbp)
        except:
            raise InvalidUsage(
                f"Invalid coordinates input '{regiontext}'", status_code=410
            )
    if chrom < 1 or chrom > 23:
        raise InvalidUsage("Chromosome input must be between 1 and 23", status_code=410)
    elif startbp > endbp:
        raise InvalidUsage(
            "Starting chromosome basepair position is greater than ending basepair position",
            status_code=410,
        )
    elif startbp > maxChromLength or endbp > maxChromLength:
        raise InvalidUsage("Start or end coordinates are out of range", status_code=410)
    elif (endbp - startbp) > GENOMIC_WINDOW_LIMIT:
        raise InvalidUsage(
            f"Entered region size is larger than {GENOMIC_WINDOW_LIMIT/1e6} Mbp",
            status_code=410,
        )
    else:
        return chrom, startbp, endbp


def fetch_snv(chrom, bp, ref, build):
    variantid = "."

    if ref is None or ref == ".":
        ref = ""

    # Ensure valid region:
    try:
        regiontxt = str(chrom) + ":" + str(bp) + "-" + str(int(bp) + 1)
    except:
        raise InvalidUsage(f"Invalid input for {str(chrom):str(bp)}")
    chrom, startbp, endbp = parse_region_text(regiontxt, build)
    chrom = str(chrom).replace("chr", "").replace("23", "X")

    # Load dbSNP151 SNP names from region indicated
    dbsnp_filepath = ""
    if build.lower() in ["hg38", "grch38"]:
        suffix = "b38"
        dbsnp_filepath = os.path.join(
            app.config["LF_DATA_FOLDER"], "dbSNP151", "GRCh38p7", "All_20180418.vcf.gz"
        )
    else:
        suffix = "b37"
        dbsnp_filepath = os.path.join(
            app.config["LF_DATA_FOLDER"], "dbSNP151", "GRCh37p13", "All_20180423.vcf.gz"
        )

    # Load variant info from dbSNP151
    tbx = pysam.TabixFile(dbsnp_filepath)  # type: ignore
    varlist = []
    for row in tbx.fetch(str(chrom), bp - 1, bp):
        rowlist = str(row).split("\t")
        chromi = rowlist[0].replace("chr", "")
        posi = rowlist[1]
        idi = rowlist[2]
        refi = rowlist[3]
        alti = rowlist[4]
        varstr = "_".join([chromi, posi, refi, alti, suffix])
        varlist.append(varstr)

    # Check if there is a match to an SNV with the provided info
    if len(varlist) == 1:
        variantid = varstr
    elif len(varlist) > 1 and ref != "":
        for v in varlist:
            if v.split("_")[2] == ref:
                variantid = v
                break
    return variantid


def x_to_23(l):
    """
    Given a list of chromosome strings,
    return list where all variations of string 'X' are converted to integer 23.
    Also checks that all values fall within integer range [1, 23], or is "."
    """
    newl = []
    validchroms = [str(i) for i in list(np.arange(1, 24))]
    validchroms.append(".")
    for x in l:
        if str(str(x).strip().lower().replace("chr", "").upper()) == "X":
            newl.append(23)
        elif str(str(x).strip().lower().replace("chr", "")) in validchroms:
            if x != ".":
                newl.append(int(str(x).strip().lower().replace("chr", "")))
            else:
                newl.append(".")
        else:
            raise InvalidUsage(f"Chromosome '{x}' unrecognized", status_code=410)
    return newl


def write_list(alist, filename):
    with open(filename, "w") as f:
        for item in alist:
            f.write("%s\n" % item)


def write_matrix(aMat, filename):
    aMat = np.matrix(aMat)
    with open(filename, "w") as f:
        for row in np.arange(aMat.shape[0]):
            for col in np.arange(aMat.shape[1] - 1):
                f.write("%s\t" % str(aMat[row, col]))
            f.write("%s\n" % str(aMat[row, -1]))


def getLeadSNPindex(leadsnpname, summaryStats, snpcol, pcol):
    lead_snp = leadsnpname
    snp_list = list(summaryStats.loc[:, snpcol])
    snp_list = [
        asnp.split(";")[0] for asnp in snp_list
    ]  # cleaning up the SNP names a bit
    if lead_snp == "":
        lead_snp = list(
            summaryStats.loc[
                summaryStats.loc[:, pcol] == min(summaryStats.loc[:, pcol])
            ].loc[:, snpcol]
        )[0].split(";")[0]
    if lead_snp not in snp_list:
        raise InvalidUsage(
            f"Lead SNP '{lead_snp}' not found in dataset", status_code=410
        )
    lead_snp_position_index = snp_list.index(lead_snp)
    return lead_snp_position_index


def clean_snps(variantlist, regiontext, build):
    """
    Parameters
    ----------
    variantlist : list
        list of variant IDs in rs id or chr_pos, chr_pos_ref_alt, chr_pos_ref_alt_build, etc formats
    regiontext : str
        the region of interest in chr:start-end format
    build : str
        build.lower() in ['hg19','hg38', 'grch37', 'grch38'] must be true

    Returns
    -------
    A cleaner set of SNP names
        rs id's are cleaned to contain only one,
        non-rs id formats are standardized to chr_pos_ref_alt_build format)
        any SNPs not in regiontext are returned as '.'
    """

    variantlist = [
        asnp.split(";")[0].replace(":", "_").replace(".", "") for asnp in variantlist
    ]  # cleaning up the SNP names a bit
    std_varlist = standardize_snps(variantlist, regiontext, build)
    final_varlist = [
        e if (e.startswith("rs") and std_varlist[i] != ".") else std_varlist[i]
        for i, e in enumerate(variantlist)
    ]

    return final_varlist
