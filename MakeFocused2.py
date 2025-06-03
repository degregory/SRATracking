#!/bin/env python3
import os
import json
import datetime
import time
import numpy as np
import pandas as pd

data_file = "SARS-CoV-2_Variant_Proportions.tsv"
start_date = (datetime.datetime.today() - datetime.timedelta(days=70)).strftime('%Y-%m-%d ')

data_df = pd.read_csv(data_file, sep='\t',)

data_df = data_df[data_df["usa_or_hhsregion"] == "USA"]
data_df = data_df[data_df["modeltype"] == "smoothed"]

data_df = data_df[['week_ending', 'variant', 'share']]

data_df = data_df[data_df["week_ending"] > start_date]
data_df = data_df[data_df["share"] > 0]
lins = set(data_df['variant'].tolist())

del_dict = {}

del_in = "VariantPMsLong.tsv"
data_df = pd.read_csv(del_in, sep='\t',)
data_df = data_df[["NT Change", "ORFs", "AA Change"]] #
data_df.drop_duplicates(inplace=True)

del_df = data_df[data_df['NT Change'].str.contains('del')]
del_df['NT Change'] = del_df['NT Change'].str.strip("ATCG")
del_df['NT Change'] = 'Del:' + del_df['NT Change'].str.strip("del")
del_df['AA Change'] = del_df['ORFs'] + ":" + del_df['AA Change']

del_dict = pd.Series(del_df["AA Change"].tolist() ,index=del_df['NT Change']).to_dict()

del_dict['Del:21653-21655'] = 'S:NS30-31Ndel'

data_df['AA_changes'] = data_df['ORFs'] + ":" + data_df['AA Change']
AA_changes = set(data_df['AA_changes'].tolist())

match_dict = {}
for change in AA_changes:
    change = change.split(":")
    if not change[0] in match_dict:
        match_dict[change[0]] = {}
    match_dict[change[0]][change[1][1:]] = change[1]

PM_dict = {}

for var in lins:
    if var == "Other":
        continue
    path = ""
    for c in var:
        if c == ".":
            path += "/"
        elif not c.isnumeric():
            path += c + "/"
        else:
            path += c
    path = f"/mnt/d/Tree/sars-cov-2-lineage-dominant-mutations/{path}/{var}-muts.txt"
    if not os.path.isfile(path):
        print(f"{var} not found")
        print(path)
    else:
        with open(path, "r") as in_fh:
            PM_dict[var] = in_fh.read().strip().split("\n")
            PM_dict[var] = [entry for entry in PM_dict[var] if not entry.startswith("Nuc")]

if not "JN.1" in PM_dict:
    with open("/mnt/d/Tree/sars-cov-2-lineage-dominant-mutations/J/N/1/JN.1-muts.txt", "r") as in_fh:
        PM_dict['JN.1'] = in_fh.read().strip().split("\n")

all_dels = []
matched_pms = []

lin_dict = {}
for lin in PM_dict:
    # print(lin)
    if not lin == 'JN.1':
        PM_dict[lin] = [entry for entry in PM_dict[lin] if not entry in PM_dict['JN.1']]
        all_dels += [entry for entry in PM_dict[lin] if entry.startswith("Del:")]

    for deletion in del_dict:
        if deletion in PM_dict[lin]:
            PM_dict[lin].append(del_dict[deletion])

    lin_dict[lin] = {}
    for ent in PM_dict[lin]:
        ent = ent.split(":")
        if ent[0] in ("Del", "Ins"):
            continue
        if "ORF1" in ent[0]:
            ent[0] = "ORF1"
        if len(ent) < 2:
            print(f"{lin} {ent}")
            continue
        try:
            lin_dict[lin][ent[0]]
        except:
            lin_dict[lin][ent[0]] = []
        if "del" in ent[1]:
            lin_dict[lin][ent[0]].append(ent[1])
        else:
            try:
                lin_dict[lin][ent[0]].append(match_dict[ent[0]][ent[1]])
                matched_pms.append(f"{ent[0]} {ent[1]}")
            except:
                if not lin == "JN.1":
                    print(f"{ent[0]} {ent[1]} not matched")

for d in set(all_dels):
    try:
        del_dict[d]
    except:
        print(f"{d} not found")
        print(" ")

with open("major_lineages.json", "w") as j_out:
    j_out.write("[\n")
    lin_outs = []
    for lin in lin_dict:
        if lin_dict[lin]:
            lin_str = "\t{\n"
            lin_str += f'\t"lineage": "{lin}",\n'
            lin_str += '\t"mutations": "'
            mut_list = []
            for orf in lin_dict[lin]:
                for mut in lin_dict[lin][orf]:
                    mut_list.append(f"{orf}:{mut}")
            lin_str += ','.join(mut_list)
            lin_str += '"\n'
            lin_str += "\t}"
            lin_outs.append(lin_str)
    j_out.write((",\n").join(lin_outs))
    j_out.write("\n]\n")