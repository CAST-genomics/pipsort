import pandas as pd
import numpy as np

import sys

processedfiles = sys.argv[1]
outfile = sys.argv[2]


f = open(processedfiles, 'r')

lines = f.readlines()
num_files = len(lines)

count = 0
assert(len(lines) >= 2)
df_all_pos = pd.read_csv(lines[0].strip(), sep=" ")[["rsid"]]
print(df_all_pos.shape)
print(df_all_pos)
for idx in range(1, len(lines)):
    count += 1
    df_next = pd.read_csv(lines[idx].strip(), sep = " ")[["rsid"]]
    print(df_next.shape)
    df_all_pos = df_all_pos.merge(df_next, how="outer", on=["rsid"])
    print(df_all_pos)

print("df all pos done")
print(df_all_pos.shape)
print(df_all_pos)
#df_all_pos.sort_values(by="pos", inplace=True)

f.close()
f = open(processedfiles, 'r')
lines = f.readlines()
print(len(lines))

num_studies = 0
for idx in range(len(lines)):
    num_studies += 1
    df = pd.read_csv(lines[idx].strip(), sep=" ")[["rsid","pos"]]
    print(df.shape)
    print(df)
    df_all_pos = df_all_pos.merge(df, on=["rsid"], how="left")
    print("df all pos after merge")
    print(df_all_pos.shape)
    col_label = "study" + str(idx)
    #df_all_pos[col_label] = df['pos']
    df_all_pos.rename(columns={"pos":col_label}, inplace=True)

print("all merging done")
print(df_all_pos[:10])
#df_all_pos.sort_values(by=["pos","rsid"], inplace=True)

cols_to_min_over = []
for i in range(num_studies):
    cols_to_min_over.append("study"+str(i))
df_all_pos['min'] = df_all_pos[cols_to_min_over].min(axis=1)
df_all_pos.sort_values("min", inplace=True, ignore_index=True)
#df_all_pos.sort_values("pos", inplace=True, ignore_index=True)
idx = np.arange(df_all_pos.shape[0], dtype=int)
curr_pos = df_all_pos["min"][0]
curr_idx = 0
for i in range(1, df_all_pos.shape[0]):
    next_pos = df_all_pos["min"][i]
    if (next_pos != curr_pos):
        curr_idx += 1
    idx[i] = curr_idx
    curr_pos = next_pos
df_all_pos["idx"] = idx
df_all_pos.sort_values(["idx","rsid"], inplace=True, ignore_index=True)
print("after column-wise sorting")
print(df_all_pos[:10])
df_all_pos.drop(columns=["min","idx"],inplace=True)
print(df_all_pos[:10])

def mk_idx(col):
    counter = 0
    newcol = col.copy()
    print(newcol)
    for i in range(len(col)):
        print("this")
        print(type(newcol[i]))
        if newcol[i] == pd.NA:
            newcol[i] = -1
        else:
            newcol[i] = counter
            counter += 1
    return newcol

print("num files", num_files)
for i in range(num_files):
    col_label = "study" + str(i)
    print("fixing", col_label)
    na_idxs = df_all_pos[col_label].isna()
    not_na = df_all_pos[col_label].notna()
    idx = np.arange(np.sum(not_na))
    df_all_pos[col_label][not_na] = idx
    df_all_pos[col_label][na_idxs] = -1
    df_all_pos[col_label] = df_all_pos[col_label].astype('int')

print(df_all_pos[:10])
#df_all_pos.drop(columns=['pos'], inplace=True)
df_all_pos.to_csv(outfile, sep=',', header=False, index=False)
