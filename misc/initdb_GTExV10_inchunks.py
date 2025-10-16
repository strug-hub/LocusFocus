#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan  2 21:20:01 2021

@author: naim
Adapted from: https://maxhalford.github.io/blog/pandas-streaming-groupby/
Works well in Ubuntu
"""

import cProfile
import itertools
# import multiprocessing as mp # Mongo does not support several write requests
import pandas as pd
import numpy as np
import os
from pathlib import Path
import pymongo
from pymongo import MongoClient
import subprocess
from datetime import datetime


def stream_groupby_parquet(path, key, agg, chunk_size=1e6, pool=None, **kwargs):

    # Make sure path is a list
    if not isinstance(path, list):
        path = [path]

    # Chain the chunks
    kwargs['chunksize'] = chunk_size
    chunks = itertools.chain(*[
        pd.read_parquet(p, **kwargs)
        for p in path
    ])

    results = []
    orphans = pd.DataFrame()
        
    for chunk in chunks:

        # Add the previous orphans to the chunk
        chunk = pd.concat((orphans, chunk))

        # Determine which rows are orphans
        last_val = chunk[key].iloc[-1]
        is_orphan = chunk[key] == last_val

        # Put the new orphans aside
        chunk, orphans = chunk[~is_orphan], chunk[is_orphan]

        # If a pool is provided then we use apply_async
        if pool:
            results.append(pool.apply_async(agg, args=(chunk,)))
        else:
            results.append(agg(chunk))

    # If a pool is used then we have to wait for the results
    if pool:
        results = [r.get() for r in results]
    
    results.append(agg(orphans)) # ensure last chunk (gene) is pushed as well!

    return pd.concat(results)


def agg(chunk):
    """lambdas can't be serialized so we need to use a function"""
    chunk.set_index('gene_id', inplace=True)
    return chunk.groupby('gene_id').apply(push_variant_dict)


def push_variant_dict(collection):
    def wrapped(gene_df):
        
        variant_id = [str(x).replace('chr','').encode('utf-8') for x in list(gene_df['variant_id'])]
        pval = list(gene_df['pval_nominal'])
        beta = list(gene_df['slope'])
        se = list(gene_df['slope_se'])
        ma_samples = list(gene_df['ma_samples'])
        ma_count = list(gene_df['ma_count'])
        sample_maf = list(gene_df['af'])
        geneid = gene_df.reset_index()['gene_id'][0]
        variants_list = []
        for row in np.arange(len(variant_id)):
            variants_list.append({
                'variant_id': variant_id[row].decode('utf-8')
                , 'pval': float(pval[row])
                , 'beta': float(beta[row])
                , 'se': float(se[row])
                , 'ma_samples': float(ma_samples[row])
                , 'ma_count': float(ma_count[row])
                , 'sample_maf': float(sample_maf[row])                     
                })
        gene_dict = {'gene_id': geneid, 'eqtl_variants': variants_list }
        collection.insert_one(gene_dict)
    return wrapped


##########################################
# MAIN
##########################################

conn = "mongodb://localhost:27017"
client = MongoClient(conn)
db = client.GTEx_V10

# Adipose_Subcutaneous.v10.allpairs.chr1.parquet
# <tissue>.v10.allpairs.<chrom>.parquet

GTEX_ROOT = os.path.join('data', 'GTEx_v10_eQTL')
files_list = list(map(lambda name: os.path.join(GTEX_ROOT, name), filter(lambda x: x.endswith('.parquet'), os.listdir(GTEX_ROOT))))
tissues = list(map(lambda x: x.split('.')[0].replace(' ','_'), files_list))
completed_files_file = os.path.join(GTEX_ROOT, 'completed_files.txt')
if not os.path.isfile(completed_files_file):
    # touch file
    Path(completed_files_file).touch()

completed_filenames = set()
with open(completed_files_file, 'r', encoding='utf-8') as f:
    for line in f:
        completed_filenames.add(line.strip())

files_df = pd.DataFrame(
    {
        "filename": files_list,
        "tissue": [x.split("/")[-1].split('.')[0].replace(' ','_') for x in files_list],
        "chrom": [x.split('.')[3] for x in files_list]
    }
)

test_tile = None
test_file = files_df.iloc[202:203]

def write_to_db(df):
    """filename df groupby helper"""
    tissue = df.iloc[0]['tissue']
    if tissue not in db.list_collection_names():
        collection = db.create_collection(tissue)
        print(tissue + ' collection created')
    collection = db[tissue]
    chroms = df["chrom"].unique()
    for chrom in chroms:
        file = df[df['chrom'] == chrom]['filename'].iloc[0]
        print('Parsing file ' + file)
        if file in completed_filenames:
            print('Skipping file ' + file)
            continue
        
        tissue_eqtls = pd.read_parquet(
            file,
            columns=['gene_id', 'variant_id', 'af', 'ma_samples', 'ma_count', 'pval_nominal', 'slope', 'slope_se']
        )
        tissue_eqtls.groupby('gene_id').apply(push_variant_dict(collection))  # type: ignore
        with open("completed_files.txt", "a", encoding='utf-8') as f:
            f.write(file + '\n')

        del tissue_eqtls  # try and free up memory
    print(datetime.now().strftime('%c'))
    print('Now indexing by gene_id')
    print(datetime.now().strftime('%c'))
    collection.create_index('gene_id')
    print('Indexing done')
    print(datetime.now().strftime('%c'))
    print('Done with tissue ' + tissue)
    print(datetime.now().strftime('%c'))


#test_file.apply(write_to_db)
cProfile.run("test_file.groupby('tissue').apply(write_to_db)", filename="test.prof")
#files_df.groupby("tissue").apply(write_to_db) # type: ignore

# no variant lookup table
