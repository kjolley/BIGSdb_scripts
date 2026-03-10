# Neisseria gonorrhoeae azithromycin resistance determination.
# Inputs are JSON files from the PubMLST API for the following
# loci: 23S_rRNA, pro_NEIS1635, and NEIS1633.
# You also need a JSON input containing a list of objects, each
# consisting of id, 23S_rRNA (allele number), NEIS1633 (allele number),
# and pro_NEIS1635 (allele number)

# Written by Made Krisna.
# Modified by Keith Jolley.

import os
import json
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument(
    "--rRNA_23S",
    required=True,
    help="Path of JSON input file for 23S_rRNA alleles.",
)
parser.add_argument(
    "--pro_NEIS1635",
    required=True,
    help="Path of JSON input file for pro_NEIS1635 alleles.",
)
parser.add_argument(
    "--NEIS1633",
    required=True,
    help="Path of JSON input file for NEIS1633 alleles.",
)
parser.add_argument(
    "--isolates",
    required=True,
    help="Path of isolate data JSON input file. This should contain "
    "values for id, 23S_rRNA, NEIS1633, and pro_NEIS1635.",
)
parser.add_argument(
    "--output",
    required=False,
    help="Output file path. Send to STDOUT if not specified.",
)
args = parser.parse_args()

loci = ["23S_rRNA", "pro_NEIS1635", "NEIS1633"]
variants = [["SNP2599", "SNP2047"], ["SNP12"], ["SAV823", "SAV854", "SAV59"]]

file_paths = [args.rRNA_23S, args.pro_NEIS1635, args.NEIS1633]

input_a = []
for path in file_paths:
    if not os.path.exists(path):
        exit(f"{path} does not exist.")
    with open(path, "r") as file:
        data = json.load(file)
        input_a.append(data)

if os.path.exists(args.isolates):
    with open(args.isolates, "r") as file:
        input_b = json.load(file)
else:
    exit(f"{args.isolates} does not exist.")

for idx, locus_info in enumerate(input_a):
    new_locus_info = []
    for allele_info in locus_info:
        complete = list(allele_info.keys())
        retain = variants[idx]
        retain.append("locus")
        retain.append("allele_id")
        remove = list(set(complete) - set(retain))
        for item in remove:
            del allele_info[item]

isolate_ids = []
var23_SNP2047s = []
var23_SNP2599s = []
varpro_SNP12s = []
varmtrd_SAV59s = []
varmtrd_SAV823s = []
varmtrd_SAV854s = []

for index, isolate_dict in enumerate(input_b):
    isolate_id = isolate_dict["id"]

    allele_23s = isolate_dict[loci[0]]
    allele_pro = isolate_dict[loci[1]]
    allele_mtrd = isolate_dict[loci[2]]

    if allele_23s != "":
        for idx, dict_allele in enumerate(input_a[0]):
            if dict_allele["allele_id"] == allele_23s:
                var23_SNP2047 = dict_allele["SNP2047"]
                var23_SNP2599 = dict_allele["SNP2599"]
    else:
        var23_SNP2047 = "no value"
        var23_SNP2599 = "no value"

    if allele_pro != "":
        for idx, dict_allele in enumerate(input_a[1]):
            if dict_allele["allele_id"] == allele_pro:
                varpro_SNP12 = dict_allele["SNP12"]
    else:
        varpro_SNP12 = "no value"

    if allele_mtrd != "":
        for idx, dict_allele in enumerate(input_a[2]):
            if dict_allele["allele_id"] == allele_mtrd:
                varmtrd_SAV59 = dict_allele["SAV59"]
                varmtrd_SAV823 = dict_allele["SAV823"]
                varmtrd_SAV854 = dict_allele["SAV854"]
    else:
        varmtrd_SAV59 = "no value"
        varmtrd_SAV823 = "no value"
        varmtrd_SAV854 = "no value"

    isolate_ids.append(isolate_id)
    var23_SNP2047s.append(var23_SNP2047)
    var23_SNP2599s.append(var23_SNP2599)
    varpro_SNP12s.append(varpro_SNP12)
    varmtrd_SAV59s.append(varmtrd_SAV59)
    varmtrd_SAV823s.append(varmtrd_SAV823)
    varmtrd_SAV854s.append(varmtrd_SAV854)

df_main = pd.DataFrame(
    {
        "id": isolate_ids,
        "23s_A2059G": var23_SNP2047s,
        "23s_C2599T": var23_SNP2599s,
        "12AC_subs": varpro_SNP12s,
        "NEIS1633_G59D": varmtrd_SAV59s,
        "NEIS1633_K823D/E": varmtrd_SAV823s,
        "NEIS1633_F854L": varmtrd_SAV854s,
    }
)

results = []
results_idtrack = []

for index, row in df_main.iterrows():
    if row["12AC_subs"] == "WT (A)":
        if row["23s_C2599T"] == "C2599T":  # SNP2599
            results.append("R")
            results_idtrack.append(row["id"])
        elif row["23s_C2599T"] != "WT (C)" and row["23s_C2599T"] != "C2599T":
            if row["23s_C2599T"] == "no value":
                results.append("insufficient data")
                results_idtrack.append(row["id"])
            else:
                results.append("S")
                results_idtrack.append(row["id"])
        else:
            if row["NEIS1633_K823D/E"] != "WT (K)":
                if row["NEIS1633_K823D/E"] == "no value":
                    results.append("insufficient data")
                    results_idtrack.append(row["id"])
                else:
                    results.append("S")
                    results_idtrack.append(row["id"])
            else:
                if row["23s_A2059G"] == "WT (A)":  # SNP2047 on PubMLST
                    results.append("S")
                    results_idtrack.append(row["id"])
                elif row["23s_A2059G"] == "no value":
                    results.append("insufficient data")
                    results_idtrack.append(row["id"])
                else:
                    results.append("R")
                    results_idtrack.append(row["id"])
    else:
        if row["NEIS1633_F854L"] == "WT (F)":
            results.append("S")
            results_idtrack.append(row["id"])
        else:
            if row["23s_C2599T"] == "WT (C)":
                if row["NEIS1633_G59D"] == "WT (G)":
                    results.append("S")
                    results_idtrack.append(row["id"])
                elif row["NEIS1633_G59D"] == "no value":
                    results.append("insufficient data")
                    results_idtrack.append(row["id"])
                else:
                    results.append("R")
                    results_idtrack.append(row["id"])
            elif row["23s_C2599T"] == "no value":
                results.append("insufficient data")
                results_idtrack.append(row["id"])
            else:
                results.append("R")
                results_idtrack.append(row["id"])

df_result = pd.DataFrame({"id": results_idtrack, "prediction": results})

if args.output:
    df_result.to_json(args.output, orient="records", compression="infer")
else:
    print(df_result.to_json(orient="records", compression="infer"))
