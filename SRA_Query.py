#!/bin/env python3

import os
import sys
import xml.parsers.expat
import datetime

# date search range
# start
s_y = '2025'
s_m = '06'
s_d = '03'
# end
e_y = '2030'
e_m = '12'
e_d = '31'

searchstr = f'(sars-cov-2%20wastewater)%20AND%20(%22{s_y}%2F{s_m}%2F{s_d}%22%5BPublication%20Date%5D%20%3A%20%22{e_y}%2F{e_m}%2F{e_d}%22%5BPublication%20Date%5D)'

os.system(f"curl -A 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:50.0) Gecko/20100101 Firefox/50.0' 'https://www.ncbi.nlm.nih.gov/sra/?term={searchstr}' -o SRA_Search_Results.html")

MCID = ""
Key = ""

ts = datetime.datetime.now().strftime("%y.%m.%d.%H.%M")
print(ts)


with open("SRA_Search_Results.html", "r") as search_res:
    MCID = search_res.read().split('value="MCID_')[1].split('"')[0]
    search_res.seek(0)
    Key = search_res.read().split("query_key:&quot;")[1].split("&quot")[0]
print(MCID)
print(Key)
os.system(f"curl 'https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/sra-db-be.cgi?rettype=exp&WebEnv=MCID_{MCID}&query_key={Key}' -L -o SRAmetadata.{ts}.xml")
os.system(f"curl 'https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/sra-db-be.cgi?rettype=acclist&WebEnv=MCID_{MCID}&query_key={Key}' -L -o FullSRRList.{ts}.txt")
os.system(f"curl 'https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/sra-db-be.cgi?rettype=runinfo&WebEnv=MCID_{MCID}&query_key={Key}' -L -o RunInfo.{ts}.csv")

release_dict = {}

with open(f"RunInfo.{ts}.csv", "r") as info_in:
    for line in info_in:
        line = line.split(",")
        release_dict[line[0]] = (line[1], line[2])

with open(f"SRAmetadata.{ts}.xml", "r") as xml_fh:
    with open(f"SRA_meta.{ts}.txt", "w") as txt_fh:
        with open(f"SRA_meta.{ts}.tsv", "w") as tsv_fh:

            elements = []
            val_dict = {
                "SRR_acc" : "",
                "BioProject" : "",
                "BioSamp" : "",
                "Submitter" : "",
                "Col_Date" : "",
                "GeoLoc" : "",
                "Population" : "",
                "Reads" : "",
            }
            flag_dict = {
                "SRR_acc" : 0,
                "BioProject" : 0,
                "BioSamp" : 0,
                "Submitter" : 0,
                "Col_Date" : 0,
                "GeoLoc" : 0,
                "Population" : 0,
            }

            parse_xml = xml.parsers.expat.ParserCreate()
            def start_element(name, attrs):
                elements.append(name)
                if len(elements) > 1:
                    if len(elements) > 2:
                        for ele_ in elements[2:]:
                            txt_fh.write("    ")
                    txt_fh.write(name)
                    txt_fh.write(" :")
                    if attrs:
                        txt_fh.write(" ")
                        txt_fh.write(str(attrs))
                    txt_fh.write("\n")
                    if "BioProject" in str(attrs):
                        flag_dict["BioProject"] = 1
                    if "BioSamp" in str(attrs):
                        flag_dict["BioSamp"] = 1
                    if name == "RUN" and attrs:
                        val_dict["SRR_acc"] = attrs["accession"]
                        try:
                            val_dict["Reads"] = attrs["total_spots"]
                        except:
                            pass
                    if name == "Member" and attrs:
                        try:
                            val_dict["Reads"] = attrs["spots"]
                        except:
                            pass
                    if name == "SUBMISSION" and attrs:
                        try:
                            val_dict["Submitter"] = attrs["center_name"]
                        except:
                            pass
                    if name == "SUBMITTER_ID" and attrs and not val_dict["Submitter"]:
                        try:
                            val_dict["Submitter"] = attrs["namespace"]
                        except:
                            pass

            def end_element(name):
                elements.pop()
                if elements and len(elements) < 2:
                    tsv_fh.write("\t".join(val_dict.values()))
                    try:
                        date_str = "\t".join(release_dict[val_dict["SRR_acc"]])
                        tsv_fh.write(f"\t{date_str}")
                    except KeyError:
                        tsv_fh.write("\t0\t0")
                    tsv_fh.write("\n")
                    for entry in val_dict:
                        val_dict[entry] = ""
                    txt_fh.write("\n")

            def char_data(data):
                if data and data.strip():
                    if len(elements) > 1:
                        for ele in elements[1:]:
                            txt_fh.write("    ")
                    txt_fh.write(data)
                    txt_fh.write("\n")
                    if flag_dict["BioProject"] == 1 and data.startswith("PR"):
                        flag_dict["BioProject"] = 0
                        val_dict["BioProject"] = data
                    if flag_dict["BioSamp"] == 1 and data.startswith("SAM"):
                        flag_dict["BioSamp"] = 0
                        val_dict["BioSamp"] = data
                    if flag_dict["Col_Date"] == 1:
                        val_dict["Col_Date"] = data
                        flag_dict["Col_Date"] = 0
                    if data in ("collection_date", "collection date"):
                        flag_dict["Col_Date"] = 1
                    if flag_dict["GeoLoc"] == 1:
                        val_dict["GeoLoc"] += data + ", "
                        flag_dict["GeoLoc"] = 0
                    if data in ("geo_loc_name", "geo loc name") or "geographic location" in data:
                        flag_dict["GeoLoc"] = 1
                    if flag_dict["Population"] == 1:
                        if not val_dict["Population"]:
                            val_dict["Population"] = data
                        elif val_dict["Population"] in (">", "<") and data.isnumeric():
                            val_dict["Population"] += data
                        else:
                            flag_dict["Population"] = 0
                    if data in ("ww_population", ):
                        flag_dict["Population"] = 1

            parse_xml.StartElementHandler = start_element
            parse_xml.EndElementHandler = end_element
            parse_xml.CharacterDataHandler = char_data
            xml_data = xml_fh.read().split("</EXPERIMENT_PACKAGE_SET>\n<EXPERIMENT_PACKAGE_SET>")
            parse_xml.Parse("\n".join(xml_data))
            


os.system(f"cat SRA_meta.{ts}.tsv >> SRA_meta.tsv")
os.system("python3 clean_meta.py") # script that removes duplicate entries
