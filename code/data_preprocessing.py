import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter

data = pd.read_csv("../data_ext2/Extremophiles_GTDB_Radio.tsv", sep = "\t")
new_data = pd.read_excel("../data_ext2/radiophiles.xlsx")

h_index = new_data[new_data["Temp Label"] == 'Hyperthermophiles'].index
n_index = new_data[new_data["Temp Label"].isna()].index

new_labels = []
for n, v in enumerate(new_data["Temp Label"].values):
    if n in h_index:
        new_labels.append("Thermophiles")
    elif n in n_index:
        new_labels.append("none")
    else:
        new_labels.append(v)




new_new_data = pd.DataFrame()
for c in data.columns:
    new_new_data[c] = data[c]

new_new_data["Temp_label"] = new_labels

new_new_data.to_csv("../data_ext2/Extremophiles_GTDB_Radio_n.csv")

d = pd.read_csv("../data_ext2/new_data.csv")
print(d.columns)
#
# h_index = data[data["Temp_label"].isna()].index
# data["Temp_label"].iloc[h_index] = "none"
#
# data.to_csv("../data_ext2/Extremophiles_GTDB_Radio.tsv")
#
# print(Counter(data["Temp_label"]))
#

