import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter

data = pd.read_csv("../data_ext2/Extremophiles_GTDB_Radio.tsv")
new_data = pd.read_excel("../data_ext2/radiophiles.xlsx")


h_index = new_data[new_data["Temp Label"] == 'Hyperthermophiles'].index
new_data["Temp Label"].iloc[h_index] = "Thermophiles"

n_index = new_data[new_data["Temp Label"].isna()].index
new_data["Temp Label"].iloc[n_index] = "none"


data["Temp_label"] = list(new_data["Temp Label"])
data.to_csv("../data_ext2/Extremophiles_GTDB_Radio_new.csv")

#
# h_index = data[data["Temp_label"].isna()].index
# data["Temp_label"].iloc[h_index] = "none"
#
# data.to_csv("../data_ext2/Extremophiles_GTDB_Radio.tsv")
#
# print(Counter(data["Temp_label"]))
#

