#!/bin/env python3
import os
import sys
import datetime
import time
import numpy as np
import pandas as pd

new_dl = True
data_url = "https://data.cdc.gov/api/views/jr58-6ysp/rows.tsv"
data_file = "SARS-CoV-2_Variant_Proportions.tsv"
start_date = (datetime.datetime.today() - datetime.timedelta(days=70)).strftime('%Y-%m-%d ')
if new_dl:
    os.system(f"mv -f {data_file} {data_file}.old")
    os.system(f"curl {data_url} -v -o {data_file}")


exit()

linage_rank_dict = {}

data_df = pd.read_csv(data_file, sep='\t',)

data_df = data_df[data_df["usa_or_hhsregion"] == "USA"]
data_df = data_df[data_df["modeltype"] == "smoothed"]

data_df = data_df[['week_ending', 'variant', 'share']]

data_df = data_df[data_df["week_ending"] > start_date]
data_df = data_df[data_df["share"] > 0]
data_df["repeat"] = data_df.groupby(['week_ending', 'variant']).cumcount()
data_df['tag'] = data_df["week_ending"] + data_df["repeat"].astype(str)
data_df = data_df.pivot(index='variant', columns='tag', values='share')

data_df["avg"] = data_df.mean(axis=1)
data_df = data_df.reset_index()

data_df = data_df[['variant', 'avg']]

data_df = data_df.sort_values(by='avg', ascending=False)

# data_df['variant'] = data_df['variant'] + "*"

ranked_vars = data_df['variant'].tolist()

print(ranked_vars)

