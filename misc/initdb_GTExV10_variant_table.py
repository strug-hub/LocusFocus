#!/usr/bin/env python3

import pandas as pd
import pymongo
from pymongo import MongoClient
from datetime import datetime
import os

# Initialize MongoDB database
conn = "mongodb://localhost:27017"
client = MongoClient(conn)
db = client.GTEx_V10

# create the variant lookup table
print('Reading variant lookup file GTEx_Analysis_2021-02-11_v10_WholeGenomeSeq_953Indiv.lookup_table.txt.gz')
collection = db['variant_table']
tbl_chunk = pd.read_csv(os.path.join('data','GTEx_Analysis_2021-02-11_v10_WholeGenomeSeq_953Indiv.lookup_table.txt.gz'),
                        sep="\t", chunksize=100, encoding='utf-8')
print('Pushing variant information into variant_table collection by chunks')
for tbl in tbl_chunk:
    tbl['chr'] = [int(str(x).replace('chr','').replace('X','23')) for x in list(tbl['chr'])]
    tbl['variant_id'] = [x.replace('chr','') for x in list(tbl['variant_id'])]
    tbl_dict = tbl.to_dict(orient='records')
    collection.insert_many(tbl_dict)
print('Variant collection created')
print(datetime.now().strftime('%c'))

print('Now indexing by variant_id')
collection.create_index('variant_id')
print('Indexing by variant_id done')
print(datetime.now().strftime('%c'))
print('Now indexing by chr and pos (ascending) order')
collection.create_index([("chr", pymongo.ASCENDING),
                         ("variant_pos", pymongo.ASCENDING)])
print('Indexing by chr and pos done')
print(datetime.now().strftime('%c'))
