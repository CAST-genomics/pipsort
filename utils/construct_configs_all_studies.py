import sys
import pandas as pd
import numpy as np
import itertools
import os

imp_snps_file = sys.argv[1]
num_snps_file = sys.argv[2]
outfile = sys.argv[3]
dim_file = sys.argv[4]

f1 = open(imp_snps_file, "r")
lines = f1.readlines()
num_studies = 0
dfs = []
all_empty = True
for line in lines:
    num_studies += 1
    filename = line.strip()
    if ((not os.path.exists(filename)) or (os.path.getsize(filename) == 0)):
        dfs.append(None)
        print("empty file name")
    else:
        all_empty = False
        df = pd.read_csv(filename, sep="\t", header=None)
        dfs.append(df)
f1.close()

if all_empty:
    np.array(-1, dtype=np.int16).tofile(outfile)
    with open(dim_file, 'w') as f:
        f.write("1\n")
        f.write("1")
    exit(0)


num_snps_per_study = np.zeros(num_studies)
f2 = open(num_snps_file, "r")
lines = f2.readlines()
cntr = 0
for line in lines:
    num_snps_per_study[cntr] = int(np.loadtxt(line.strip()))
    cntr += 1
f2.close()
assert(cntr == num_studies)


def get_arrs_for_study(imp_snps, global_offset):
    sets, set_sizes = np.unique(imp_snps[3], return_counts=True)
    num_sets = len(sets)
    arrs = []
    study_offset = 0
    for i in range(num_sets):
        a = np.zeros(set_sizes[i]+1, dtype=np.int16) #+1 for the "no causal in this group" aka "all zeros for this group" config
        a[0] = -1 #undef index will be -1, could change to 0 for memory optimization purposes but then remaining snp pos need to be 1-indexed instead of the current 0-indexing
        for j in range(set_sizes[i]):
            a[j+1] = imp_snps.iloc[study_offset+j][2] + global_offset
        study_offset += set_sizes[i]
        arrs.append(a)
        print(a)
    return arrs

all_arrs = []
global_offset = 0
for i in range(num_studies):
    imp_snps = dfs[i]
    if i > 0: #is important that this steps happens even if imp_snps is empty
        global_offset += num_snps_per_study[i-1] #set global offset to position at end of prev study
    if imp_snps is None:
        #all_arrs += [np.array([-1])]
        continue
    else:
        all_arrs += get_arrs_for_study(imp_snps, global_offset)


'''

vec = np.zeros(num_snps)

mats = []
offset = 0
for i in range(num_sets):
    mat = np.zeros((set_sizes[i]+1,num_snps),dtype=np.int8) #+1 for an all zeros config
    for j in range(set_sizes[i]):
        #print(j+1)
        #print(imp_snps.iloc[j][2])
        mat[j+1,imp_snps.iloc[offset+j][2]] = 1 #+1 because first config is all zeros config
    #print(mat)
    offset += set_sizes[i]
    #print("offset")
    #print(offset)
    mats.append(mat)

for mat in mats:
    print(mat.shape)

num_all_configs = np.prod(set_sizes+1)
print("num_all_configs: ", num_all_configs)
configs = np.zeros((num_all_configs, num_snps))
'''


#https://stackoverflow.com/questions/11144513/cartesian-product-of-x-and-y-array-points-into-single-array-of-2d-points/49445693#49445693
def cartesian_product_transpose_pp(arrays):
   la = len(arrays)
   #dtype = np.result_type(*arrays)
   arr = np.empty((la, *map(len, arrays)), dtype=np.int16)
   idx = slice(None), *itertools.repeat(None, la)
   for i, a in enumerate(arrays):
      arr[i, ...] = a[idx[:la-i]]
   return arr.reshape(la, -1).T

print(all_arrs)
all_configs = cartesian_product_transpose_pp(all_arrs)
print(all_configs[:10])
print(all_configs.shape)
prod = np.prod([len(a) for a in all_arrs])
assert(prod == all_configs.shape[0])

'''
arrs = []
for i in range(num_sets):
    arrs.append(np.arange(set_sizes[i]+1))
all_combos = cartesian_product_transpose_pp(arrs)
print("all combos")
print(all_combos)
print(len(all_combos))
'''

def sum_all(idxs):
    r = np.zeros(num_snps)
    for i in range(num_sets):
        r += mats[i][idxs[i]]
    return r

def special_sort(idxs):
    non_neg = idxs[np.where(idxs >= 0)]
    non_neg = np.sort(non_neg)
    j = 0
    for i in range(idxs.shape[0]):
        if idxs[i] >= 0:
            idxs[i] = non_neg[j]
            j += 1
    return idxs



sorted_all_configs = np.apply_along_axis(special_sort, 1, all_configs)
print(sorted_all_configs[:10])
#configs = np.apply_along_axis(sum_all, 1, all_combos)
#print(configs)
#print(configs.shape)
#np.savetxt("temp_all_configs_not_binary.txt", sorted_all_configs, delimiter=',', fmt='%5d')
sorted_all_configs.astype(np.int16,copy=False).tofile(outfile)
print(sorted_all_configs.dtype)
with open(dim_file, 'w') as f:
  f.write('%d\n' % sorted_all_configs.shape[0])
  f.write('%d' % sorted_all_configs.shape[1])

#all_configs.tofile(outfile)
