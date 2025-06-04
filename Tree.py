#!/bin/env python3

import os
import sys
from Bio import Phylo

dl_file = 'public-latest.all.masked.pb.gz'

refresh = True

if refresh:
    os.system(f"wget -O {dl_file} http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/public-latest.all.masked.pb.gz")
    os.system(f"matUtils extract -i {dl_file} -C LineageDefinitions.tsv")
    os.system(f"matUtils extract -i {dl_file} -t sars.nwk")



    ## "matUtils summary --translate translate_output.tsv -i {dl_file} -g SARS2.gtf -f SARS2.fasta"
    ## os.system("wget -O sars.nwk.gz https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/public-latest.all.nwk.gz")



    if os.path.isfile("sars.nwk.gz"):
        os.system("gzip -f -d sars.nwk.gz")

    SARS_Tree = Phylo.read("sars.nwk", "newick")

    with open("SARS2Tree.txt", "w") as out_fh:
        out_fh.write(str(SARS_Tree))




node = ""
with open("LineageDefinitions.tsv", "r") as in_fh:
    node = in_fh.read().split("JN.1\t")[1].split("\t")[0]

with open("SARS2Tree.txt", "r") as in_fh, open("JN.1.Tree.txt", "w") as out_fh:
    depth = 0
    rel_depth = 0
    in_node = False
    for line in in_fh:
        if node in line:
            in_node = True
            depth = len(line.split("Clade")[0])

        elif in_node and "node" in line:
            rel_depth = len(line.split("Clade")[0]) - depth
            if rel_depth <= 0:
                in_node = False
        if in_node and "node" in line:
            out_fh.write(" "*rel_depth)
            out_fh.write(line.strip())
            out_fh.write("\n")

nodes = []
with open("JN.1.Tree.txt", "r") as in_fh:
    for line in in_fh:
        nodes.append(line.split("name=")[-1].strip().strip("')"))

node_dict = {}
with open("LineageDefinitions.tsv", "r") as in_fh:
    for line in in_fh:
        line = line.split("\t")
        if line[1] in nodes:
            node_dict[line[1]] = line[0]

with open("JN.1.Tree.txt", "r") as in_fh, open("JN.1.Tree.ann.txt", "w") as out_fh:
    for line in in_fh:
        out_fh.write(line.strip("\n"))
        try:
            out_fh.write(node_dict[line.split("name=")[-1].strip().strip("')")])
        except:
            pass
        out_fh.write("\n")


var_nodes_dict = {}
with open("JN.1.Tree.ann.txt", "r") as in_fh:
    depth = 0
    cur_nodes = []
    for line in in_fh:
        depth = int(len(line.split("Clade")[0]) / 4)
        line = line.strip()
        if cur_nodes and depth < len(cur_nodes):
            cur_nodes = cur_nodes[:depth]
        node = line.split("name=")[-1].split(")")[0].strip("'")
        if not line.endswith(")"):
            var = line.split(")")[-1]
            var_nodes = []
            for entry in cur_nodes:
                if " " in entry:
                    var_nodes.append(entry.split(" ")[1])
            var_nodes_dict[var] = var_nodes
            node += " " + var
        cur_nodes.append(node)

PM_dict = {}

for var in var_nodes_dict:
    path = ""
    for c in var:
        if c == ".":
            path += "/"
        elif not c.isnumeric():
            path += c + "/"
        else:
            path += c
    path = f"sars-cov-2-lineage-dominant-mutations/{path}/{var}-muts.txt"
    if not os.path.isfile(path):
        print(f"{var} not found")
        print(path)
    else:
        with open(path, "r") as in_fh:
            PM_dict[var] = in_fh.read().strip().split("\n")


with open("JN.1.AAPM.tsv", "w") as out_fh:
    for var in var_nodes_dict:
        if not var in PM_dict:
            continue
        prev_PMs = []
        if var_nodes_dict[var]:
            for lin in var_nodes_dict[var]:
                prev_PMs += PM_dict[lin]
        PMs = []
        for PM in PM_dict[var]:
            if not PM in prev_PMs:
                PMs.append(PM)
        PMs = ";".join(PMs)
        out_fh.write(f"{var_nodes_dict[var]}\t{var}\t{PMs}\n")

# with open("Highlight_PMs.tsv", "w") as out_fh:
    # for var in var_nodes_dict:
        # if not var in PM_dict:
            # continue
        # prev_PMs = []
        # if var_nodes_dict[var]:
            # for lin in var_nodes_dict[var]:
                # prev_PMs += PM_dict[lin]
                # break
        # PMs = []
        # for PM in PM_dict[var]:
            # if not PM in prev_PMs:
                # PMs.append(PM)
        # PMs = ";".join(PMs)
        # out_fh.write(f"{var_nodes_dict[var]}\t{var}\t{PMs}\n")
