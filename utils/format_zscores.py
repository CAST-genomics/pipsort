import sys
import pandas as pd

infile = sys.argv[1]
outfile = sys.argv[2]

df = pd.read_csv(infile, sep=" ")
df = df[["rsid","Zscore"]]
df.to_csv(outfile, header=False, index=False, sep="\t")


