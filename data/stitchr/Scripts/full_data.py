import pandas as pd
import Scripts.stitchr as stitchr


def drop_hle_head(d: str) -> str:
    d = d.replace("HLA-", "")
    return d


def drop_vj_head(d: str) -> str:
    d = str(d)
    if d.startswith("TCRBV"):
        return d.replace("TCRBV", "TRBV")
    if d.startswith("TCRBJ"):
        return d.replace("TCRBJ", "TRBJ")
    if d.startswith("TCRAV"):
        return d.replace("TCRAV", "TRAV")
    if d.startswith("TCRAJ"):
        return d.replace("TCRAJ", "TRAJ")
    return d


def handle_iedb(dir):
    iedb = pd.read_csv(dir, sep=',', dtype=str)
    iedb = iedb[iedb['Description'].str.isalpha()]

    harmonized_iedb = iedb[['Chain 1 CDR3 Curated',
                            'Chain 1 CDR3 Calculated',
                            'Chain 2 CDR3 Curated',
                            'Chain 2 CDR3 Calculated',
                            'Description',
                            'MHC Allele Names',
                            'Curated Chain 1 V Gene',
                            'Calculated Chain 1 V Gene',
                            'Curated Chain 1 J Gene',
                            'Calculated Chain 1 J Gene',
                            'Curated Chain 2 V Gene',
                            'Calculated Chain 2 V Gene',
                            'Curated Chain 2 J Gene',
                            'Calculated Chain 2 J Gene']]
    # Coalesce teh 'curated' and 'calculated' features into single columns, using 'calculated' to fill in missing values
    harmonized_iedb['Chain 1 CDR3 Curated'].fillna(harmonized_iedb['Chain 1 CDR3 Calculated'], inplace=True)
    harmonized_iedb['Chain 2 CDR3 Curated'].fillna(harmonized_iedb['Chain 2 CDR3 Calculated'], inplace=True)
    harmonized_iedb.drop(['Chain 1 CDR3 Calculated',
                          'Chain 2 CDR3 Calculated',
                          'Calculated Chain 1 V Gene',
                          'Calculated Chain 1 J Gene',
                          'Calculated Chain 2 V Gene',
                          'Calculated Chain 2 J Gene'], axis=1, inplace=True)
    harmonized_columns = ['CDR3A', 'CDR3B', 'Antigen', 'HLA', 'VA', 'JA', 'VB', 'JB']
    harmonized_iedb.columns = harmonized_columns
    # harmonized_iedb.drop(["cdr3a", "v_a", "j_a"], axis=1, inplace=True)
    harmonized_iedb.drop_duplicates(inplace=True)
    # re = harmonized_iedb.dropna()
    re = harmonized_iedb.reset_index(drop=True)
    drop_index = []
    # for i, HLA in enumerate(re["HLA"]):
    #     if "H2" in str(HLA):
    #         drop_index.append(i)
    for i, Antigen in enumerate(re["Antigen"]):
        if len(Antigen) > 15:
            drop_index.append(i)
    re = re.drop(drop_index)
    re = re.reset_index(drop=True)

    # re["VB"] = re.apply(lambda x: drop_vj_head(x["VB"]), axis=1)
    # re["JB"] = re.apply(lambda x: drop_vj_head(x["JB"]), axis=1)
    # re["VA"] = re.apply(lambda x: drop_vj_head(x["VA"]), axis=1)
    # re["JA"] = re.apply(lambda x: drop_vj_head(x["JA"]), axis=1)

    # CDR3 = re["CDR3"].tolist()
    # cdr3a = re["cdr3a"].tolist()
    # Vb = re["v_b"].tolist()
    # Jb = re["j_b"].tolist()
    # Va = re["v_a"].tolist()
    # Ja = re["j_a"].tolist()
    # b_seq = []
    # cdr3b_re = []
    # cdr3a_re = []
    # a_seq = []
    # for i in range(len(CDR3)):
    #     cdr3 = str(CDR3[i])
    #     if cdr3.startswith("A"):
    #         cdr3 = "C" + cdr3
    #     if cdr3.endswith("Y") or cdr3.endswith("T") or cdr3.endswith("TH"):
    #         cdr3 = cdr3 + "F"
    #     cdr3b_re.append(cdr3)
    #     # try:
    #     #     b_seq.append(stitchr.pre_full_seq(cdr3, Vb[i], Jb[i]))
    #     # except Exception as e:
    #     #     print(e)
    #     #     print(CDR3[i], Vb[i], Jb[i])
    #     #     b_seq.append("nan")
    #     cdr3 = str(cdr3a[i])
    #     if cdr3.startswith("A"):
    #         cdr3 = "C" + cdr3
    #     if not cdr3.endswith("F"):
    #         cdr3 = cdr3 + "F"
    #     cdr3a_re.append(cdr3)
    #     # try:
    #     #     a_seq.append(stitchr.pre_full_seq(cdr3, Va[i], Ja[i]))
    #     # except Exception as e:
    #     #     print(e)
    #     #     print(cdr3a[i], Va[i], Ja[i])
    #     #     a_seq.append("nan")
    # CDR3_PD = pd.DataFrame(
    #     {"CDR3b": cdr3b_re, "CDR3a": cdr3a_re, "Antigen": re['Antigen'],
    #      "HLA": re['HLA'], "a_seq": a_seq, "b_seq": b_seq},
    # )
    return re


def handle_vdj(dir):
    data = pd.read_csv(dir, sep="\t", index_col=False)
    # data = data[data["Gene"] == "TRB"]
    CDR3_PD = pd.DataFrame(
        {"CDR3": data['CDR3'], "Antigen": data['Epitope'], "HLA": data['MHC A'], "V": data["V"], "J": data["J"],
         "Gene": data["Gene"]},
    )
    CDR3_B = CDR3_PD[CDR3_PD["Gene"] == "TRB"]
    CDR3_A = CDR3_PD[CDR3_PD["Gene"] == "TRA"]
    CDR3_B = CDR3_B.rename(columns={'V': 'VB', 'J': 'JB', "CDR3": "CDR3B"})
    CDR3_A = CDR3_A.rename(columns={'V': 'VA', 'J': 'JA', "CDR3": "CDR3A"})
    re = CDR3_B.append(CDR3_A)
    return re


import tidytcells as it


def full(HLA_data):
    CDR3B = HLA_data["CDR3B"].tolist()
    cdr3A = HLA_data["CDR3A"].tolist()
    Vb = HLA_data["VB"].tolist()
    Jb = HLA_data["JB"].tolist()
    Va = HLA_data["VA"].tolist()
    Ja = HLA_data["JA"].tolist()
    HLA = HLA_data["HLA"].tolist()
    b_seq = []
    cdr3b_re = []
    cdr3a_re = []
    a_seq = []
    for i in range(len(CDR3B)):
        print(i)
        cdr3 = str(CDR3B[i])
        if not cdr3.endswith("F"):
            cdr3 = cdr3 + "F"
        try:
            if HLA[i].startswith("H2"):
                b_seq.append(stitchr.pre_full_seq(cdr3, it.tr.standardise(Vb[i]), it.tr.standardise(Jb[i]), "MOUSE"))
            else:
                b_seq.append(stitchr.pre_full_seq(cdr3, it.tr.standardise(Vb[i]), it.tr.standardise(Jb[i])))
            # b_seq.append(stitchr.pre_full_seq("CAASRPPNFGNEKLT", "TRAV13-1", "TRAJ48"))
        except Exception as e:
            print(e)
            print(CDR3B[i], Vb[i], Jb[i])
            b_seq.append("nan")

        cdr3 = str(cdr3A[i])
        if not cdr3.endswith("F"):
            cdr3 = cdr3 + "F"
        try:
            if HLA[i].startswith("H2"):
                a_seq.append(stitchr.pre_full_seq(cdr3, it.tr.standardise(Va[i]), it.tr.standardise(Ja[i]), "MOUSE"))
            else:
                a_seq.append(stitchr.pre_full_seq(cdr3, it.tr.standardise(Va[i]), it.tr.standardise(Ja[i])))
        except Exception as e:
            print(e)
            print(cdr3A[i], Va[i], Ja[i])
            a_seq.append("nan")

    HLA_data["TCRB"] = b_seq
    HLA_data["TCRA"] = a_seq
    # HLA_data["CDR3B"] = cdr3b_re
    # HLA_data["CDR3A"] = cdr3a_re
    return HLA_data


def read_data():
    HLA_data = handle_vdj("/Users/ljx/Code/zkTcrPmhc/tcr_pmhc/data3/vdjdb/SearchTable-2022-08-22 02_07_34.028.tsv")
    data_PMT = handle_iedb("/Users/ljx/Code/zkTcrPmhc/tcr_pmhc/data3/iedb/2022-01-25.csv")
    HLA_data = HLA_data.append(data_PMT, ignore_index=True)
    HLA_data = HLA_data.dropna(subset=["CDR3B", "CDR3A"], how="all")
    HLA_data = HLA_data.dropna(subset=["Antigen"], how="all")
    HLA_data = HLA_data.dropna(subset=["VA", "JA", "VB", "JB"], how="all")
    HLA_data = HLA_data.drop_duplicates(keep='first')
    HLA_data = full(HLA_data)
    HLA_data.to_csv("tcr_save", index=False)


if __name__ == '__main__':
    read_data()
