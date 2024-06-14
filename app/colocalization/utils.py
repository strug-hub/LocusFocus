import os
import re
from typing import List, Optional

import pysam
import pandas as pd
import numpy as np
from flask import Request, current_app as app
from werkzeug.utils import secure_filename
from werkzeug.exceptions import RequestEntityTooLarge

from app import mongo
from app.routes import InvalidUsage

GENOMIC_WINDOW_LIMIT = 2e6


def download_file(request: Request, extensions: List[str]) -> Optional[str]:
    """
    Download the first file at 'files[]' that matches the given extensions,
    and return the filepath to the saved file. Return None if no such file exists.

    Extensions should not include the period. eg. `["html", "tsv", "txt"]`
    """
    if not ('files[]' in request.files):
        raise InvalidUsage(f"No files found in request")
    
    filenames = request.files.getlist('files[]')

    saved_filepath = None
    for file in filenames:
        if file.filename is None or file.filename.rsplit(".", 1)[-1].lower() not in extensions:
            continue

        filename = secure_filename(file.filename)
        filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
        file.save(filepath)
        if not os.path.isfile(filepath):
            raise RequestEntityTooLarge(f"File '{filename}' too large")
        
        saved_filepath = filepath
        break

    return saved_filepath


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
    chromlist = [x.split('_')[0] if len(x.split('_'))==5 else x for x in variant_list]
    chromlist = [int(x) if x not in ["X","."] else x for x in chromlist]
    poslist = [int(x.split('_')[1]) if len(x.split('_'))==5 else x for x in variant_list]
    reflist = [x.split('_')[2] if len(x.split('_'))==5 else x for x in variant_list]
    altlist = [x.split('_')[3] if len(x.split('_'))==5 else x for x in variant_list]
    df = pd.DataFrame({
        "CHROM": chromlist,
        "POS": poslist,
        "REF": reflist,
        "ALT": altlist,
    })
    return df


def standardize_snps(variantlist, regiontxt, build):
    """
    Input: Variant names in any of these formats: rsid, chrom_pos_ref_alt, chrom:pos_ref_alt, chrom:pos_ref_alt_b37/b38
    Output: chrom_pos_ref_alt_b37/b38 variant ID format, but looks at GTEx variant lookup table first.
    In the case of multi-allelic variants (e.g. rs2211330(T/A,C)), formats such as 1_205001063_T_A,C_b37 are accepted
    If variant ID format is chr:pos, and the chr:pos has a unique biallelic SNV, then it will be assigned that variant
    """

    if all(x=='.' for x in variantlist):
        raise InvalidUsage('No variants provided')

    if np.nan in variantlist:
        raise InvalidUsage('Missing variant IDs detected in row(s): ' + str([ i+1 for i,x in enumerate(variantlist) if str(x) == 'nan' ]))

    # Ensure valid region:
    chrom, startbp, endbp = parse_region_text(regiontxt, build)
    chrom = str(chrom).replace('23',"X")

    # Load GTEx variant lookup table for region indicated
    db = mongo.cx.GTEx_V7 # type: ignore
    rsid_colname = 'rs_id_dbSNP147_GRCh37p13'
    if build.lower() in ["hg38", "grch38"]:
        db = mongo.cx.GTEx_V8 # type: ignore
        rsid_colname = 'rs_id_dbSNP151_GRCh38p7'
    collection = db['variant_table']
    variants_query = collection.find(
        { '$and': [
            { 'chr': int(chrom.replace('X','23')) },
            { 'variant_pos': { '$gte': int(startbp), '$lte': int(endbp) } }
            ]}
        )
    variants_list = list(variants_query)
    variants_df = pd.DataFrame(variants_list)
    variants_df = variants_df.drop(['_id'], axis=1)


    # Load dbSNP151 SNP names from region indicated
    dbsnp_filepath = ''
    suffix = 'b37'
    if build.lower() in ["hg38", "grch38"]:
        suffix = 'b38'
        dbsnp_filepath = os.path.join(app.config["LF_DATA_FOLDER"], 'dbSNP151', 'GRCh38p7', 'All_20180418.vcf.gz')
    else:
        suffix = 'b37'
        dbsnp_filepath = os.path.join(app.config["LF_DATA_FOLDER"], 'dbSNP151', 'GRCh37p13', 'All_20180423.vcf.gz')


    # Load dbSNP file
    #delayeddf = delayed(pd.read_csv)(dbsnp_filepath,skiprows=getNumHeaderLines(dbsnp_filepath),sep='\t')
    #dbsnp = dd.from_delayed(delayeddf)
    tbx = pysam.TabixFile(dbsnp_filepath) # type: ignore
#    print('Compiling list of known variants in the region from dbSNP151')
    chromcol = []
    poscol = []
    idcol = []
    refcol = []
    altcol = []
    variantid = [] # in chr_pos_ref_alt_build format
    rsids = dict({}) # a multi-allelic variant rsid (key) can be represented in several variantid formats (values)
    for row in tbx.fetch(str(chrom), startbp, endbp):
        rowlist = str(row).split('\t')
        chromi = rowlist[0].replace('chr','')
        posi = rowlist[1]
        idi = rowlist[2]
        refi = rowlist[3]
        alti = rowlist[4]
        varstr = '_'.join([chromi, posi, refi, alti, suffix])
        chromcol.append(chromi)
        poscol.append(posi)
        idcol.append(idi)
        refcol.append(refi)
        altcol.append(alti)
        variantid.append(varstr)
        rsids[idi] = [varstr]
        altalleles = alti.split(',') # could have more than one alt allele (multi-allelic)
        if len(altalleles)>1:
            varstr = '_'.join([chromi, posi, refi, altalleles[0], suffix])
            rsids[idi].append(varstr)
            for i in np.arange(len(altalleles)-1):
                varstr = '_'.join([chromi, posi, refi, altalleles[i+1], suffix])
                rsids[idi].append(varstr)

#    print('Cleaning and mapping list of variants')
    variantlist = [asnp.split(';')[0].replace(':','_').replace('.','') for asnp in variantlist] # cleaning up the SNP names a bit
    stdvariantlist = []
    for variant in variantlist:
        if variant == '':
            stdvariantlist.append('.')
            continue
        variantstr = variant.replace('chr','')
        if re.search("^23_",variantstr): variantstr = variantstr.replace('23_','X_',1)
        if variantstr.startswith('rs'):
            try:
                # Here's the difference from the first function version (we look at GTEx first)
                if variant in list(variants_df[rsid_colname]):
                    stdvar = variants_df['variant_id'].loc[ variants_df[rsid_colname] == variant].to_list()[0]
                    stdvariantlist.append(stdvar)
                else:
                    stdvariantlist.append(rsids[variantstr][0])
            except:
                stdvariantlist.append('.')
        elif re.search(r"^\d+_\d+_[A,T,G,C]+_[A,T,C,G]+,*", variantstr.replace('X','23')):
            strlist = variantstr.split('_')
            strlist = list(filter(None, strlist)) # remove empty strings
            try:
                achr, astart, aend = parse_region_text(strlist[0]+":"+strlist[1]+"-"+str(int(strlist[1])+1), build)
                achr = str(achr).replace('23','X')
                if achr == str(chrom) and astart >= startbp and astart <= endbp:
                    variantstr = variantstr.replace("_"+str(suffix),"") + "_"+str(suffix)
                    if len(variantstr.split('_')) == 5:
                        stdvariantlist.append(variantstr)
                    else:
                        raise InvalidUsage(f'Variant format not recognizable: {variant}. Is it from another coordinate build system?', status_code=410)
                else:
                    stdvariantlist.append('.')
            except:
                raise InvalidUsage(f'Problem with variant {variant}', status_code=410)
        elif re.search(r"^\d+_\d+_*[A,T,G,C]*", variantstr.replace('X','23')):
            strlist = variantstr.split('_')
            strlist = list(filter(None, strlist)) # remove empty strings
            try:
                achr, astart, aend = parse_region_text(strlist[0]+":"+strlist[1]+"-"+str(int(strlist[1])+1), build)
                achr = str(achr).replace('23','X')
                if achr == str(chrom) and astart >= startbp and astart <= endbp:
                    if len(strlist)==3:
                        aref=strlist[2]
                    else:
                        aref=''
                    stdvariantlist.append(fetch_snv(achr, astart, aref, build))
                else:
                    stdvariantlist.append('.')
            except:
                raise InvalidUsage(f'Problem with variant {variant}', status_code=410)
        else:
            raise InvalidUsage(f'Variant format not recognized: {variant}', status_code=410)
    return stdvariantlist


def parse_region_text(regiontext, build):
    if build not in ['hg19', 'hg38']:
        raise InvalidUsage(f'Unrecognized build: {build}', status_code=410)
    regiontext = regiontext.strip().replace(' ','').replace(',','').replace('chr','')
    if not re.search(r"^\d+:\d+-\d+$", regiontext.replace('X','23').replace('x','23')):
       raise InvalidUsage(f'Invalid coordinate format. {regiontext} e.g. 1:205,000,000-206,000,000', status_code=410)
    chrom = regiontext.split(':')[0].lower().replace('chr','').upper()
    pos = regiontext.split(':')[1]
    startbp = pos.split('-')[0].replace(',','')
    endbp = pos.split('-')[1].replace(',','')
    chromLengths = pd.read_csv(os.path.join(app.config["LF_DATA_FOLDER"], build + '_chrom_lengths.txt'), sep="\t", encoding='utf-8')
    chromLengths.set_index('sequence',inplace=True)
    if chrom in ['X','x'] or chrom == '23':
        chrom = 23
        maxChromLength = chromLengths.loc['chrX', 'length']
        try:
            startbp = int(startbp)
            endbp = int(endbp)
        except:
            raise InvalidUsage(f"Invalid coordinates input: {regiontext}", status_code=410)
    else:
        try:
            chrom = int(chrom)
            if chrom == 23:
                maxChromLength = chromLengths.loc['chrX', 'length']
            else:
                maxChromLength = chromLengths.loc['chr'+str(chrom), 'length']
            startbp = int(startbp)
            endbp = int(endbp)
        except:
            raise InvalidUsage(f"Invalid coordinates input {regiontext}", status_code=410)
    if chrom < 1 or chrom > 23:
        raise InvalidUsage('Chromosome input must be between 1 and 23', status_code=410)
    elif startbp > endbp:
        raise InvalidUsage('Starting chromosome basepair position is greater than ending basepair position', status_code=410)
    elif startbp > maxChromLength or endbp > maxChromLength:
        raise InvalidUsage('Start or end coordinates are out of range', status_code=410)
    elif (endbp - startbp) > GENOMIC_WINDOW_LIMIT:
        raise InvalidUsage(f'Entered region size is larger than {GENOMIC_WINDOW_LIMIT/1e6} Mbp', status_code=410)
    else:
        return chrom, startbp, endbp
    

def fetch_snv(chrom, bp, ref, build):
    variantid = '.'

    if ref is None or ref=='.':
        ref=''

    # Ensure valid region:
    try:
        regiontxt = str(chrom) + ":" + str(bp) + "-" + str(int(bp)+1)
    except:
        raise InvalidUsage(f'Invalid input for {str(chrom):str(bp)}')
    chrom, startbp, endbp = parse_region_text(regiontxt, build)
    chrom = str(chrom).replace('chr','').replace('23',"X")

    # Load dbSNP151 SNP names from region indicated
    dbsnp_filepath = ''
    if build.lower() in ["hg38", "grch38"]:
        suffix = 'b38'
        dbsnp_filepath = os.path.join(app.config["LF_DATA_FOLDER"], 'dbSNP151', 'GRCh38p7', 'All_20180418.vcf.gz')
    else:
        suffix = 'b37'
        dbsnp_filepath = os.path.join(app.config["LF_DATA_FOLDER"], 'dbSNP151', 'GRCh37p13', 'All_20180423.vcf.gz')

    # Load variant info from dbSNP151
    tbx = pysam.TabixFile(dbsnp_filepath) #type: ignore
    varlist = []
    for row in tbx.fetch(str(chrom), bp-1, bp):
        rowlist = str(row).split('\t')
        chromi = rowlist[0].replace('chr','')
        posi = rowlist[1]
        idi = rowlist[2]
        refi = rowlist[3]
        alti = rowlist[4]
        varstr = '_'.join([chromi, posi, refi, alti, suffix])
        varlist.append(varstr)

    # Check if there is a match to an SNV with the provided info
    if len(varlist) == 1:
        variantid = varstr
    elif len(varlist) > 1 and ref != '':
        for v in varlist:
            if v.split('_')[2] == ref:
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
    validchroms = [str(i) for i in list(np.arange(1,24))]
    validchroms.append('.')
    for x in l:
        if str(str(x).strip().lower().replace('chr','').upper()) == "X":
            newl.append(23)
        elif str(str(x).strip().lower().replace('chr','')) in validchroms:
            if x!='.':
                newl.append(int(str(x).strip().lower().replace('chr','')))
            else:
                newl.append('.')
        else:
            raise InvalidUsage('Chromosome unrecognized', status_code=410)
    return newl