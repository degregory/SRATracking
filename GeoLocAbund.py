#!/bin/env python3
import os
import sys
import datetime
import time
import pysam
import multiprocessing as mp
import numpy as np
import pandas as pd
import gzip

pysam.set_verbosity(0)


regions_dict = {
    1 : "Connecticut, Maine, Massachusetts, New Hampshire, Rhode Island, Vermont",
    2 : "New Jersey, New York, Puerto Rico, the Virgin Islands",
    3 : "Delaware, District of Columbia, Maryland, Pennsylvania, Virginia, West Virginia",
    4 : "Alabama, Florida, Georgia, Kentucky, Mississippi, North Carolina, South Carolina, Tennessee",
    5 : "Illinois, Indiana, Michigan, Minnesota, Ohio, Wisconsin",
    6 : "Arkansas, Louisiana, New Mexico, Oklahoma, Texas",
    7 : "Iowa, Kansas, Missouri, Nebraska",
    8 : "Colorado, Montana, North Dakota, South Dakota, Utah, Wyoming",
    9 : "Arizona, California, Hawaii, Nevada, American Samoa, Commonwealth of the Northern Mariana Islands, Federated States of Micronesia, Guam, Marshall Islands, Republic of Palau",
    10 : "Alaska, Idaho, Oregon, Washington",
}

s_r_dict = {}
for region in regions_dict:
    for st in regions_dict[region].split(", "):
        s_r_dict[st] = region

time_var_dict = {}

for file in os.listdir():
    if file.endswith(".passed.pms.txt"):
        # print(file)
        lweek = datetime.datetime.strptime(file.split(".")[0], '%Y-%m-%d')
        weeks = []
        for i in range(0, 6):
            weeks.append((lweek - datetime.timedelta(days=(7*i))).strftime('%Y-%m-%d'))

        for week in weeks:
            if not week in time_var_dict:
                time_var_dict[week] = []
        with open(file, "r") as pms_in:
            for line in pms_in:
                pm = line.strip()
                for week in weeks:
                    time_var_dict[week].append(pm)

weeks = sorted(time_var_dict.keys(), reverse=True)
dtweeks = [datetime.datetime.strptime(week, '%Y-%m-%d') for week in weeks]



samp_dict = {}

for date in time_var_dict:
    time_var_dict[date] = set(time_var_dict[date])
    samp_dict[date] = {}

accs = {}

if os.path.isfile("found_meta.tsv"):
    with open("found_meta.tsv", "r") as met_in:
        for line in met_in:
            sline = line.split("\t")
            loc = sline[5]
            if loc.startswith("USA"):
                st = loc.split(":")[1].split(",")[0].strip()
                if st == "Calif":
                    st = "California"
                if not st in s_r_dict:
                    print(f"{st} not recognized in regions")
                    print(line)
                    exit()
                date = datetime.datetime.strptime(sline[4], '%Y-%m-%d')
                for i in range(0, len(dtweeks)):
                    if date >= dtweeks[i]:
                        accs[sline[0]] = weeks[i]
                        samp_dict[weeks[i]][sline[0]] = {"region" : s_r_dict[st]}
                        break

elif os.path.isfile("found_meta.tsv.gz"):
    with gzip.open("found_meta.tsv.gz", "rb") as met_in:
        met_in.readline()
        for line in met_in:
            line = line.decode()
            sline = line.split("\t")
            loc = sline[5]
            if loc.startswith("USA"):
                st = loc.split(":")[1].split(",")[0].strip()
                if st == "Calif":
                    st = "California"
                if not st in s_r_dict:
                    print(f"{st} not recognized in regions")
                    print(line)
                    exit()
                date = datetime.datetime.strptime(sline[4], '%Y-%m-%d')
                for i in range(0, len(dtweeks)):
                    if date >= dtweeks[i]:
                        accs[sline[0]] = weeks[i]
                        samp_dict[weeks[i]][sline[0]] = {"region" : s_r_dict[st]}
                        break

else:
    print("meta not found")
    exit()


for subdir, dirs, files in os.walk("/mnt/f/SRA/SARS2/Wastewater/"):
    for file in files:
        acc = file.split(".")[0]
        if acc in accs:
            dt = ""
            if file.endswith("covars.tsv.gz"):
                dt = "covar"
            if file.endswith("nt_calls.tsv.gz"):
                dt = "nts"
            if dt:
                cut = "uncut"
                if "cut" in file:
                    cut = "cut"
                try:
                    samp_dict[accs[acc]][acc][dt]
                except KeyError:
                    samp_dict[accs[acc]][acc][dt] = {}
                try:
                    samp_dict[accs[acc]][acc][dt][cut]
                except KeyError:
                    samp_dict[accs[acc]][acc][dt][cut] = []
                samp_dict[accs[acc]][acc][dt][cut].append(os.path.join(subdir, file))


def get_pos(PM):
    POS = PM.split("\t")[0].split(" ")[0].split(",")[0].split("|")[0].split("(")[0].split("-")[0].strip("ATGCNDel")
    return POS


def proc_week(week):
    covar_dict = {}
    cov_dict = {}
    for covar in time_var_dict[week]:
        cov_dict[get_pos(covar)] = 0

    for acc in samp_dict[week]:
        try:
            samp_dict[week][acc]["covar"]["cut"][0]
        except KeyError as err:
            pass
        else:
            with gzip.open(samp_dict[week][acc]["covar"]["cut"][0], "rb") as in_fh:
                for line in in_fh:
                    line = line.decode()
                    if "Reference" in line:
                        continue
                    line = line.split("\t")
                    if line[0] in time_var_dict[week]:
                        try:
                            covar_dict[line[0]]
                        except:
                            covar_dict[line[0]] = {}
                        try:
                            covar_dict[line[0]][samp_dict[week][acc]["region"]] +=int(line[1])
                        except:
                            covar_dict[line[0]][samp_dict[week][acc]["region"]] = int(line[1])
            try:
                samp_dict[week][acc]["nts"]["cut"][0]
            except KeyError as err:
                print(f"no nt file for {acc}")
            else:
                with gzip.open(samp_dict[week][acc]["nts"]["cut"][0], "rb") as nt_in:
                    nt_in.readline()
                    nt_in.readline()
                    for line in nt_in:
                        line = line.decode()
                        line = line.split("\t")
                        if line[0] in cov_dict:
                            cov_dict[line[0]] += int(line[9])

    with open(f"HHS/{week}.covars.tsv", "w") as out_fh:
        for covar in covar_dict:
            for region in covar_dict[covar]:
                out_fh.write(f"{covar}\t{region}\t{covar_dict[covar][region]}\t{cov_dict[get_pos(covar)]}\t{covar_dict[covar][region]/cov_dict[get_pos(covar)]}\n")


with mp.Manager() as manager:
    with mp.Pool(6) as pool:

        procings = []

        for week in weeks:
            proc = pool.apply_async(proc_week, (week, ))
            procings.append(proc)

        for proc in procings:
            proc.get()

        pool.close()
        pool.join()

with open("HHS.long.tsv", "w") as out_fh:
    out_fh.write("Position\tNT Change\tORFs\tAA Change\tWeek\tRegion\tCount\tTotal Coverage\tAbundance\n")
    for file in os.listdir("HHS"):
        with open(os.path.join("HHS", file), "r") as in_fh:
            date = file.split(".")[0]
            changes_dict = {}
            for line in in_fh:
                line = line.strip().split("\t")
                change = line[0].split("|")
                try:
                    changes_dict[change[1]]["nts"].append(change[0])
                except KeyError:
                    changes_dict[change[1]] = {
                        "nts" : [change[0]],
                    }
                except:
                    print(line)
                    exit
                changes_dict[change[1]]["counts"] = line[1:]
            for change in changes_dict:
                if len(changes_dict[change]["nts"]) == 1:
                    changes_dict[change]["nts"] = changes_dict[change]["nts"][0]
                else:
                    changes_dict[change]["nts"] = ', '.join(sorted(set(changes_dict[change]["nts"]), key=lambda x: int(get_pos(x))))
                changes_dict[change]["pos"] = get_pos(changes_dict[change]["nts"])

            for change in changes_dict:
                pm_change = change.strip("()")
                AA_change = pm_change.split("(")[-1]
                orf = pm_change.split(":")[0]
                if "ORF1a" in orf:
                    orf = "ORF1"
                elif "ORF" in orf:
                    orf = orf.split("_")[0]
                else:
                    orf = orf[0].upper()
                out_fh.write(f"{changes_dict[change]["pos"]}\t{changes_dict[change]["nts"]}\t")
                out_fh.write(f"{orf}\t{AA_change}\t{date}")
                for item in changes_dict[change]["counts"]:
                    out_fh.write(f"\t{item}")
                out_fh.write("\n")

exit()
