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

delta_thresh = .02
count_thresh = 100
file_search = True
reproc_rs = False
reproc_fin = True
update_cir = "CDC"

cur_date = datetime.datetime.today()
weekly = datetime.datetime.strptime('2024-10-06', '%Y-%m-%d') - datetime.timedelta(days=(7*3))
weektime = datetime.timedelta(days=7)
weeks = []
acc_dict = {}
found_acc = {}
found_cram = {}

while cur_date > weekly:
    date_str = weekly.strftime('%Y-%m-%d')
    acc_dict[date_str] = []
    found_acc[date_str] = []
    found_cram[date_str] = []
    weeks.append(weekly)
    weekly += weektime

weeks = weeks[::-1]

meta_dict = {}

with open("/mnt/f/SRA/SARS2/Wastewater/Meta/SRA_meta.tsv", "r") as meta_in, open("weeks_meta.tsv", "w") as out_fh:
    for line in meta_in:
        sline = line.split("\t")
        col_date = sline[4].split("/")[0]
        if 'PRJNA748354' in line:
            continue
        if 'PRJEB44932' in line:
            continue
        if "Maryland" in line and '78365' in line:
            continue
        try:
            col_dt = datetime.datetime.strptime(col_date, '%Y-%m-%d')
        except:
            continue
        if  col_dt >= weeks[-1]:
            for week in weeks:
                if col_dt >= week:
                    week_date = week.strftime('%Y-%m-%d')
                    acc_dict[week_date].append(sline[0])
                    out_fh.write(f"{week_date}\t{line}")
                    meta_dict[sline[0]] = line
                    break

data_int = False

# check data integrity
for week in acc_dict:
    data_int = False
    if os.path.isfile(f"{week}.nts.tsv") and os.path.isfile(f"{week}.covars.tsv"):
        with open(f"{week}.nts.tsv", "r") as in1, open(f"{week}.covars.tsv", "r") as in2:
            if "\nall end\n" in in1.read() and "\nall end\n" in in2.read():
                data_int = True
    if data_int == False:
        print(f"error or no data for {week}")

last_week = ''
count = 0

weeks = [week.strftime('%Y-%m-%d') for week in weeks]

total_samp = 0
for week in weeks:
    print(f"{week} : {len(acc_dict[week])}")
    total_samp += len(acc_dict[week])
print(f"total : {total_samp}")

read_acc = []
if os.path.isfile("procced_acc.txt"):
    with open("procced_acc.txt", "r") as read_in:
        for line in read_in:
            read_acc.append(line.strip())

file_dict = {}
cram_dict = {}
new = 0
for week in weeks:
    file_dict[week] = []
    cram_dict[week] = []
    for acc in acc_dict[week]:
        if not acc in read_acc:
            new += 1

if new == 0:
    print("no new samples to proc")
else:
    print(f"{new} samples to process")

if file_search:
    paths = ("/mnt/f/SRA/SARS2/Wastewater/2025/20250522/", "/mnt/f/SRA/SARS2/Wastewater/FromWI/", "/mnt/f/SRA/SARS2/Wastewater/2024/", "/mnt/f/SRA/SARS2/Wastewater/2025/")
    for path in paths:
        WI_count = 0
        for subdir, dirs, files in os.walk(path):
            for file in files:
                if file.endswith(".cut_covars.tsv.gz"):
                    acc = file.split(".")[0]
                    for week in weeks:
                        if acc in acc_dict[week] and not acc in found_acc[week]:
                            found_acc[week].append(acc)
                            if not acc in read_acc:
                                fp = os.path.join(subdir, file)
                                file_dict[week].append(fp)
                                read_acc.append(acc)
                                WI_count += 1
                                break
                elif file.endswith(".cut.cram"):
                    acc = file.split(".")[0]
                    for week in weeks:
                        if acc in acc_dict[week] and not acc in found_cram[week]:
                            found_cram[week].append(acc)

        print(WI_count)

    unfound = 0
    for week in found_acc:
        missed = len(acc_dict[week]) - len(found_acc[week])
        print(f"{week} : {len(found_acc[week])} with {missed} missing")
        unfound += missed
    print(f"{unfound} samples not found")

    with open("found_meta.tsv", "w") as out_fh, open("breaks.txt", "w") as out2_fh:
        out_fh.write("Sample\tBioProject\tBioSample\tSubmitter\tCollection Date\tLocation\tPopulation\tReads\tReleaseDate\tLoadDate\n")
        for week in found_acc:
            for acc in set(found_acc[week] + found_cram[week]):
                out_fh.write(meta_dict[acc])
            for acc in set(found_cram[week]):
                if not acc in found_acc[week]:
                    out2_fh.write(f"{acc}\n")

def proc_week(week):
    if file_dict[week]:
        print(f"Starting SR file processing for {week}")
        for file in file_dict[week]:
            sfile = file.split("/")
            acc = sfile[-1].split(".")[0]
            path = "/".join(sfile[:-1])
            if os.path.isfile(os.path.join(path, f"{acc}.cut_covars.tsv.gz")) and os.path.isfile(os.path.join(path, f"{acc}.cut_nt_calls.tsv.gz")):
                with gzip.open(os.path.join(path, f"{acc}.cut_covars.tsv.gz"), "rb") as seqs_in, open(f"{week}.covars.tsv", "a") as seq_out:
                    seqs_in.readline()
                    seqs_in.readline()
                    seq_out.write(f"{acc} start\n")
                    for line in seqs_in:
                        line = line.decode()
                        if "Reference" in line:
                            continue

                        line = line.split("\t")
                        try:
                            int(line[1])
                            seq_out.write(f"{line[0]}\t{line[1]}\t\n")
                        except:
                            continue
                    seq_out.write(f"{acc} end\n")

                time.sleep(1)
                samp_dict = {}
                with gzip.open(os.path.join(path, f"{acc}.cut_nt_calls.tsv.gz"), "rb") as nt_in, open(f"{week}.nts.tsv", "a") as nt_out:
                    nt_in.readline()
                    nt_in.readline()
                    nt_out.write(f"{acc} start\n")
                    for line in nt_in:
                        line = line.decode()
                        line = line.split("\t")
                        if line[2]:
                            nt_out.write(f"{line[0]}\t{line[9]}\n")
                    nt_out.write(f"{acc} end\n")
                time.sleep(1)
                with open("procced_acc.txt", "a") as procced:
                    procced.write(f"{acc}\n")

                continue

            else:
                print(f"No NT call file for {file} (size : {os.path.getsize(file)})")
                continue

    else:
        print(f"no new files for {week}")

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

def rs(week):
    if len(file_dict[week]) > 0 or (reproc_rs and os.path.isfile(f"{week}.covars.tsv")):
        print(f"starting new samples incorporation for {week}")
        covars_dict = {}
        with open(f"{week}.covars.tsv", "r") as in_fh:
            acc_dict = {}
            accs = []
            cur_acc = ''

            for line in in_fh:
                line = line.strip()
                if line.endswith("start"):
                    cur_acc = line.split(' ')[0]
                    acc_dict = {}
                elif line.endswith("end"):
                    if line.split(' ')[0] == cur_acc:
                        if acc_dict and not cur_acc in accs:
                            if cur_acc == "all":
                                for entry in acc_dict:
                                    covars_dict[entry] = []
                                    for i in (0, 1, 2, 3): # read count, samp count, 100 read samp count, 1000 read samp count
                                        covars_dict[entry].append(int(acc_dict[entry][i]))
                                accs.append(cur_acc)
                            else:
                                for entry in acc_dict:
                                    cur_stats = [int(acc_dict[entry][0]), 1, 0, 0]
                                    if cur_stats[0] >= 100:
                                        cur_stats[2] = 1
                                        if cur_stats[0] >= 1000:
                                            cur_stats[3] = 1
                                    try:
                                        for i in (0, 1, 2, 3):
                                            covars_dict[entry][i] += cur_stats[i]
                                    except:
                                        covars_dict[entry] = cur_stats
                                accs.append(cur_acc)

                elif "|" in line:
                    line = line.split("\t")
                    acc_dict[line[0]] = line[1:]

        with open(f"{week}.covars.tsv", "w") as out_fh:
            out_fh.write("all start\n")
            for entry in covars_dict:
                out_fh.write(entry)
                for data in covars_dict[entry]:
                    out_fh.write(f"\t{data}")
                out_fh.write("\n")
            out_fh.write("all end\n")
        nts_dict = {}

        with open(f"{week}.nts.tsv", "r") as in_fh:
            acc_dict = {}
            finished = True
            accs = []
            cur_acc = ''
            for line in in_fh:
                line = line.strip()
                if line.endswith("start"):
                    if finished and acc_dict and not cur_acc in accs:
                        for entry in acc_dict:
                            try:
                                nts_dict[entry] += acc_dict[entry]
                            except:
                                nts_dict[entry] = acc_dict[entry]
                        accs.append(cur_acc)
                    cur_acc = line.split(' ')[0]
                    finished = False
                    acc_dict = {}
                elif line.endswith("end") and line.split(' ')[0] == cur_acc:
                        finished = True
                else:
                    line = line.split("\t")
                    try:
                        acc_dict[line[0]] += int(line[1])
                    except:
                        acc_dict[line[0]] = int(line[1])

            if finished and acc_dict and not cur_acc in accs:
                for entry in acc_dict:
                    try:
                        nts_dict[entry] += acc_dict[entry]
                    except:
                        nts_dict[entry] = acc_dict[entry]

        with open(f"{week}.nts.tsv", "w") as out_fh:
            out_fh.write("all start\n")
            for entry in nts_dict:
                out_fh.write(f"{entry}\t{nts_dict[entry]}\n")
            out_fh.write("all end\n")
    else:
        print(f"no new samples for {week}")

with mp.Manager() as manager:
    with mp.Pool(6) as pool:

        refinings = []

        for week in weeks:
            refined = pool.apply_async(rs, (week, ))
            refinings.append(refined)

        for refined in refinings:
            refined.get()

        pool.close()
        pool.join()

JN_changes = []
JN_AAchanges = {}
with open("JN.1.covars.tsv", "r") as JN_in:
    JN_in.readline()
    JN_in.readline()
    for line in JN_in:
        line = line.split("\t")
        if "|" in line[0]:
            JN_changes.append(line[0])
            for pm_change in line[0].split("|")[1:]:
                pm_change = pm_change.strip("()")
                orf = pm_change.split(":")[0]
                try:
                    JN_AAchanges[orf].append(pm_change.split("(")[-1])
                except:
                    JN_AAchanges[orf] = [pm_change.split("(")[-1].strip(")")]


def get_pos(PM):
    POS = int(PM.split("\t")[0].split(" ")[0].split(",")[0].split("|")[0].split("(")[0].split("-")[0].strip("ATGCNDel"))
    return POS


JN_clade_PM_dict = {}
with open("/mnt/d/Tree/JN.1.AAPM.tsv", "r") as in_fh:
    in_fh.readline()
    for line in in_fh:
        line = line.strip().split("\t")
        if not line[1].endswith("*"):
            line[1] += "*"
        if len(line) < 3:
            continue
        for pm in line[2].split(";"):
            pm = pm.split(":")
            if "ORF1" in pm[0]:
                pm[0] = "ORF1"
            try:
                JN_clade_PM_dict[pm[0]]
            except:
                JN_clade_PM_dict[pm[0]] = {}
            try:
                JN_clade_PM_dict[pm[0]][pm[1]].append(line[1])
            except:
                JN_clade_PM_dict[pm[0]][pm[1]] = [line[1]]

ranked_vars = []

if update_cir == "WHO":
    data_url = "https://data.who.int/dashboards/covid19/variants"
    data_file = "WHOVar.html"

    os.system(f"curl {data_url} -s -v -o {data_file}")

    with open(data_file, "r") as in_fh:
        who_data = in_fh.read().split("Pango lineage")
        for entry in who_data[1:]:
            ranked_vars.append(entry.split("</strong>")[0].split("<strong>")[-1] + "*")

if update_cir == "CDC":
    data_file = "SARS-CoV-2_Variant_Proportions.tsv"
    data_df = pd.read_csv(data_file, sep='\t',)

    data_df = data_df[data_df["usa_or_hhsregion"] == "USA"]
    data_df = data_df[data_df["modeltype"] == "smoothed"]

    data_df = data_df[['week_ending', 'variant', 'share']]

    start_date = (datetime.datetime.today() - datetime.timedelta(days=70)).strftime('%Y-%m-%d ')
    data_df = data_df[data_df["week_ending"] > start_date]
    data_df = data_df[data_df["share"] > 0]
    data_df["repeat"] = data_df.groupby(['week_ending', 'variant']).cumcount()
    data_df['tag'] = data_df["week_ending"] + data_df["repeat"].astype(str)
    data_df = data_df.pivot(index='variant', columns='tag', values='share')

    data_df["avg"] = data_df.mean(axis=1)
    data_df = data_df.reset_index()

    data_df = data_df[['variant', 'avg']]
    data_df['variant'] = data_df['variant'] + "*"

    data_df = data_df.sort_values(by='avg', ascending=False)

    ranked_vars = data_df['variant'].tolist()

header = "Position\tNT Change\tORFs\tAA Change\tAssociated Variants"
for i in (0, 1):
    header += "\tTime Span\tCount\tAbundance\tTotal Coverage\tPos. Samples\t1k+ read samples\tSamples Processed"
header += "\tdiff\tProcess end date\n"

def weekly_out(j):
    new = False

    for i in range(0, 6):
        if file_dict[weeks[i+j]]:
            new = True
    if (not new) and not reproc_fin:
        print(f"nothing new to process for {weeks[j]} windows")
        return

    triweeks = [
        f"{(datetime.datetime.strptime(weeks[j], '%Y-%m-%d') - datetime.timedelta(days=14)).strftime('%Y-%m-%d')}--{(datetime.datetime.strptime(weeks[j], '%Y-%m-%d') + datetime.timedelta(days=6)).strftime('%Y-%m-%d')}",
        f"{(datetime.datetime.strptime(weeks[j+3], '%Y-%m-%d') - datetime.timedelta(days=14)).strftime('%Y-%m-%d')}--{(datetime.datetime.strptime(weeks[j+3], '%Y-%m-%d') + datetime.timedelta(days=6)).strftime('%Y-%m-%d')}",
    ]
    pm_dict = {}
    coverage_dict = {}
    for triweek in triweeks:
        pm_dict[triweek] = {}
        coverage_dict[triweek] = {}

    total_proc_samps = {}

    for i in range(0, 6):
        if i < 3:
            triweek = triweeks[0]
        else:
            triweek = triweeks[1]
        try:
            total_proc_samps[triweek]
        except:
            total_proc_samps[triweek] = 0

        i = i + j
        total_proc_samps[triweek] += len(found_acc[weeks[i]])

        with open(f"{weeks[i]}.nts.tsv", "r") as in_fh:
            for line in in_fh:
                if line.startswith("all"):
                    continue
                line = line.strip().split("\t")
                try:
                    coverage_dict[triweek][int(line[0])] += int(line[1])
                except:
                    coverage_dict[triweek][int(line[0])] = int(line[1])
        with open(f"{weeks[i]}.covars.tsv", "r") as in_fh:
            for line in in_fh:
                line = line.strip().split("\t")
                if "|" in line[0] and not line[0] in JN_changes:
                    if line[0].split("|")[0][-1] == "N":
                        continue
                    # try:
                        # pm_dict[triweek]
                    # except:
                        # pm_dict[triweek] = {}
                    try:
                        pm_dict[triweek][line[0]]
                    except:
                        pm_dict[triweek][line[0]] = [0, 0, 0, 0]
                    for k in (0, 1, 2, 3):
                        pm_dict[triweek][line[0]][k] += int(line[k+1])


    passed_pms = []
    for triweek in pm_dict:
        for change in pm_dict[triweek]:

            try:
                if pm_dict[triweek][change][0] / coverage_dict[triweek][int(get_pos(change))] > delta_thresh:
                    passed_pms.append(change)
            except:
                pass

    passed_pms = list(set(passed_pms))

    with open(f"{weeks[j]}.passed.pms.txt", "w") as pms_out:
        pms_out.write("\n".join(passed_pms))


    passed_pm_dict = {}
    for change in passed_pms:
        schange = change.split("|")
        try:
            passed_pm_dict[schange[1]]["nts"].append(schange[0])
        except:
            passed_pm_dict[schange[1]] = {
            "nts" : [schange[0]],
            }
        passed_pm_dict[schange[1]]["full"] = change

    for change in passed_pm_dict:
        if len(passed_pm_dict[change]["nts"]) == 1:
            passed_pm_dict[change]["nts"] = passed_pm_dict[change]["nts"][0]
        else:
            passed_pm_dict[change]["nts"] = ', '.join(sorted(set(passed_pm_dict[change]["nts"]), key=lambda x: get_pos(x)))
        passed_pm_dict[change]["pos"] = get_pos(passed_pm_dict[change]["nts"])

    for change in passed_pm_dict:
        for triweek in triweeks:
            try:
                passed_pm_dict[change][triweek] = pm_dict[triweek][passed_pm_dict[change]['full']]
            except:
                pass

    pm_dict = {}
    for change in passed_pm_dict:

        try:
            pm_dict[passed_pm_dict[change]["pos"]]
        except:
            pm_dict[passed_pm_dict[change]["pos"]] = {}
        for triweek in triweeks:
            try:
                pm_dict[passed_pm_dict[change]["pos"]][triweek]
            except:
                pm_dict[passed_pm_dict[change]["pos"]][triweek] ={}
        for triweek in triweeks:
            try:
                pm_dict[passed_pm_dict[change]["pos"]][triweek][change] = passed_pm_dict[change][triweek]
            except Exception as err:
                pass

    plot_dict = {}
    header = ""
    to_write = []
    matched_lin_dict = {}
    with open(f"{weeks[j]}.VariantPMs.tsv", "w") as var_out:
        var_out.write(header)

        for pos in sorted(pm_dict.keys()):
            changes = []
            for triweek in triweeks:
                try:
                    changes += pm_dict[pos][triweek].keys()
                except:
                    pass

            for change in set(changes):
                orfs = []
                AAchanges = []
                pm_change = change.strip("()")
                AA_change = pm_change.split("(")[-1]
                orf = pm_change.split(":")[0]
                if orf in JN_AAchanges and AA_change in JN_AAchanges[orf]:
                    continue
                AAchanges.append(AA_change)

                if "ORF1a" in orf:
                    orf = "ORF1"
                elif "ORF" in orf:
                    orf = orf.split("_")[0]

                else:
                    orf = orf[0].upper()

                orfs.append(orf)
                if not orfs:
                    continue
                orfs = ", ".join(list(set(orfs)))
                orf_entry = orfs

                silent = 0
                nonsilent = []
                for AAchange in list(set(AAchanges)):
                    if not AAchange[0] == AAchange[-1]:
                        nonsilent.append(AAchange)
                if len(nonsilent) == 0:
                    continue

                nts_change = passed_pm_dict[change]["nts"]
                AAchanges = ", ".join(nonsilent)
                to_write.append(f"{pos}\t{nts_change}\t{orfs}\t{AAchanges}\t")
                try:
                    okey = orf_entry
                    pmkey = AAchanges[1:]
                    if "del" in nts_change:
                        okey = "Del"
                        pmkey = nts_change.strip("ATCGdel")

                    matched_vars = list(set(JN_clade_PM_dict[okey][pmkey]))
                    ranked_matched_vars = []
                    for lin in ranked_vars:
                        if lin in matched_vars:
                            ranked_matched_vars.append(lin)
                    for lin in matched_vars:
                        if not lin in ranked_matched_vars:
                            ranked_matched_vars.append(lin)
                    matched_lin_dict[AAchanges] = ranked_matched_vars
                    to_write.append(",".join(ranked_matched_vars))
                except:
                    pass

                abunds = []
                pos_1k = []
                for triweek in triweeks:
                    to_write.append(f"\t{triweek}")
                    try:
                        pm_dict[pos][triweek][change]
                    except:
                        abund = 0
                        to_write.append(f"\t\t\t\t\t\t{total_proc_samps[triweek]}")
                    else:
                        abund = pm_dict[pos][triweek][change][0] / coverage_dict[triweek][pos]
                        to_write.append(f"\t{pm_dict[pos][triweek][change][0]}\t{abund:05f}\t{coverage_dict[triweek][pos]}")
                        to_write.append(f"\t{pm_dict[pos][triweek][change][1]}\t{pm_dict[pos][triweek][change][3]}")
                        to_write.append(f"\t{total_proc_samps[triweek]}")
                        pos_1k.append(pm_dict[pos][triweek][change][3])
                    abunds.append(abund)
                to_write.append(f"\t{abs(abunds[0] - abunds[1]):05f}\t{cur_date}\n")

                if sum(pos_1k) > 1:
                    try:
                        plot_dict[orf_entry]
                    except:
                        plot_dict[orf_entry] = {}
                    try:
                        plot_dict[orf_entry][AAchanges].append(abunds)
                    except:
                        plot_dict[orf_entry][AAchanges] = [abunds]
        var_out.writelines(to_write)

    with open("VariantPMsComp.tsv", "a") as var_out2:
        var_out2.writelines(to_write)



if not os.path.isfile("VariantPMsComp.tsv"):
    with open("VariantPMsComp.tsv", "w") as var_out2:
        var_out2.write(header)

with mp.Manager() as manager:
    with mp.Pool(8) as pool:

        outtings = []

        for j in range(0, len(weeks)-5):
            print(weeks[j])
            if len(found_acc[weeks[j]]) < 20:
                continue

            outted = pool.apply_async(weekly_out, (j, ))
            outtings.append(outted)

        for outted in outtings:
            outted.get()

        pool.close()
        pool.join()


file = "VariantPMsComp.tsv"
df = pd.read_csv(file, sep='\t', encoding='utf-8')

df = df.drop("diff", axis=1)

keep_static = ['Position', 'NT Change', 'ORFs', 'AA Change', 'Associated Variants', 'Process end date']
to_melt = [
    'Time Span',
    'Count',
    'Abundance',
    'Total Coverage',
    'Pos. Samples',
    '1k+ read samples',
    'Samples Processed',
]
to_melt2 = [
    'Time Span.1',
    'Count.1',
    'Abundance.1',
    'Total Coverage.1',
    'Pos. Samples.1',
    '1k+ read samples.1',
    'Samples Processed.1'
]

df1 = df[keep_static + to_melt ]
df2 = df[keep_static + to_melt2 ]
df2.columns = keep_static + to_melt

df = (pd.concat([df1, df2])
         .sort_index(level=0)
         .reset_index(level=0, drop=True)
         .reset_index())

df_runs = df[['Time Span', 'Process end date']]
df_runs.drop_duplicates(subset=['Time Span'], keep='last', inplace=True)
df_runs.drop_duplicates(subset=['Process end date'], keep='last', inplace=True)

good_procs = df_runs['Process end date'].to_list()

df = df[df['Process end date'].isin(good_procs)]
df.drop_duplicates(subset=['NT Change', 'Time Span'], keep='last', inplace=True)

df = df.drop("index", axis=1)


df[['Time Span Start', 'Time Span End']] = df['Time Span'].astype("string").str.split('--',expand=True) # .str.replace("[()]", "", regex=True)

df = df.drop('Time Span', axis=1)
df = df[['Position', 'NT Change', 'ORFs', 'AA Change', 'Associated Variants', 'Time Span Start', 'Time Span End', 'Count', 'Abundance', 'Total Coverage', 'Pos. Samples', '1k+ read samples', 'Samples Processed', 'Process end date']]
df['Process end date'] = df['Process end date'].astype("string").str.split(' ',expand=True)[0]

df.to_csv("VariantPMsLong.tsv", index=False, sep='\t')

exit()
