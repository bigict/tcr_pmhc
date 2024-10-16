import pandas as pd
import numpy as np
import math

# TCRB,TCRA
re = pd.read_csv("tcr_save")
mhc_str = []
re = re.dropna(subset=["TCRB", "TCRA"], how='all')

HLA = pd.read_csv("MHC_key")
name = HLA["name"].tolist()
seq = HLA["seq"].tolist()
HLA_seq_lib = {}
for i in range(len(name)):
    n = name[i].replace("*", "").replace(":", "")
    HLA_seq_lib[n] = seq[i]
MHC_list = re["HLA"].tolist()

for i in range(len(MHC_list)):
    MHC = ""
    if isinstance(MHC_list[i], float):
        mhc_str.append(MHC)
        continue
    MHC_n = MHC_list[i].split(",")[0]
    H = MHC_n.replace("*", "").replace(":", "").replace("HLA-", "")
    if H not in HLA_seq_lib.keys():
        for k in HLA_seq_lib.keys():
            if H in k:
                MHC = HLA_seq_lib[k]
                break
    else:
        MHC = HLA_seq_lib[H]
    mhc_str.append(MHC)
re["MHC_str"] = mhc_str
re.to_csv("tcr_mhc", index=False)


def clear(cdr3b, cdr3):
    re = []
    for i in range(len(cdr3b)):
        if (isinstance(cdr3b[i], float) and np.isnan(cdr3b[i])) or (isinstance(cdr3[i], float) and np.isnan(cdr3[i])):
            re.append(np.nan)
            continue
        index = cdr3b[i].find(cdr3[i])
        if index > 90:
            re.append(cdr3b[i][index - 90:index + len(cdr3[i]) + 10])
        else:
            re.append(cdr3b[i][0:index + len(cdr3[i]) + 10])
    return re


re = pd.read_csv("tcr_mhc")
re = re.replace('nan', np.nan)
# CDR3b,CDR3a
re = re.dropna(subset=["TCRB", "TCRA"], how='all')
re = re.dropna(subset=["Antigen"], how='all')
re["TCRB"] = clear(re["TCRB"].tolist(), re["CDR3B"].tolist())
re["TCRA"] = clear(re["TCRA"].tolist(), re["CDR3A"].tolist())
re = re.dropna(subset=["TCRB", "TCRA"], how='all')

pep_clusters = pd.read_csv("pep_cluster.csv")
# pep_clusters = pd.read_csv("pep_test")
p_data = pd.merge(re, pep_clusters, how='left', on="Antigen")
p_data = p_data.drop_duplicates()
p_data = p_data.dropna(subset=["Cluster"], how='all')

Cluster = list(set(pep_clusters["Cluster"].tolist()))
import random


def random_split_list_into_groups(input_list, num_groups):
    shuffled_list = random.sample(input_list, len(input_list))
    group_size = len(input_list) // num_groups
    groups = [shuffled_list[i:i + group_size] for i in range(0, len(shuffled_list), group_size)]

    return groups


groups = random_split_list_into_groups(Cluster, 5)


#
# CDR3b,CDR3a,Antigen,HLA,a_seq,b_seq,MHC_str TCRB,TCRA

def n_data(data, num=5):
    CDR3b = data["CDR3B"].tolist()
    CDR3a = data["CDR3A"].tolist()
    c = data["Cluster"].tolist()
    a_seq = data["TCRA"].tolist()
    b_seq = data["TCRB"].tolist()

    Antigen_data = data.drop_duplicates(subset=["Antigen"])
    Antigen = Antigen_data["Antigen"].tolist()
    Antigen_c = Antigen_data["Cluster"].tolist()
    HLA = Antigen_data["HLA"].tolist()
    MHC_str = Antigen_data["MHC_str"].tolist()

    CDR3b_s = []
    CDR3a_s = []
    Antigen_s = []
    HLA_s = []
    a_seq_s = []
    b_seq_s = []
    MHC_str_s = []
    for i, t in enumerate(CDR3b):
        for j, a in enumerate(Antigen):
            if c[i] != Antigen_c[j] and random.randint(1, 100) % 40 == 0:
                CDR3b_s.append(CDR3b[i])
                CDR3a_s.append(CDR3a[i])

                a_seq_s.append(a_seq[i])
                b_seq_s.append(b_seq[i])
                Antigen_s.append(a)
                HLA_s.append(HLA[j])
                MHC_str_s.append(MHC_str[j])
    n_data = pd.DataFrame(
        {"CDR3B": CDR3b_s, "CDR3A": CDR3a_s, "Antigen": Antigen_s, "HLA": HLA_s, "TCRA": a_seq_s, "TCRB": b_seq_s,
         "MHC_str": MHC_str_s})
    n_data = n_data.drop_duplicates()
    n_data = n_data.sample(data.shape[0] * num)
    return n_data


# CDR3b,CDR3a,Antigen,HLA,a_seq,b_seq,MHC_str

# data = n_data(p_data, 1)
# data.to_csv(f"n_data", index=False)
train_list = [[1, 2, 3, 4], [0, 2, 3, 4], [1, 0, 3, 4], [1, 2, 0, 4], [1, 2, 3, 0]]
for i in range(5):
    test = groups[i]
    train = []
    for j in train_list[i]:
        train.extend(groups[j])
    test_data = p_data[p_data["Cluster"].isin(test)]
    test_data_n = n_data(test_data)
    test_data["y"] = 1
    test_data_n["y"] = 0
    # test_data = test_data.append(test_data_n)
    test_data = pd.concat([test_data, test_data_n], ignore_index=True)

    train_data = p_data[p_data["Cluster"].isin(train)]
    train_data_n = n_data(train_data)
    train_data["y"] = 1
    train_data_n["y"] = 0
    # train_data = train_data.append(train_data_n)
    train_data = pd.concat([train_data, train_data_n], ignore_index=True)
    test_data.to_csv(f"train_data/test{i}", index=False)
    train_data.to_csv(f"train_data/train{i}", index=False)
