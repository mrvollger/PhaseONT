#!/usr/bin/env python
import argparse
import pysam
import pandas as pd

def read_in_lenghts(file, tag):
	reads = pysam.FastxFile(file)
	out = {"name":[], "length":[], "tag":[]}
	for rec in reads:
		out["length"].append(len(rec.sequence))
		out["name"].append(rec.name)
		out["tag"].append(tag)
	return(out)

def phased_by_length(df, length):
	df = df[df["length"]>=length]
	if df.shape[0] == 0:
		return
	c = df["tag"].append(pd.Series(["mat","pat","unk"])).value_counts()-1
	f= c/sum(c)*100
	print(f'{length}\t{f.mat:.1f}\t{c.mat}\t{f.pat:.1f}\t{c.pat}\t{f.unk:.1f}\t{c.unk}')

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("mat", help="positional input")
	parser.add_argument("pat", help="positional input")
	parser.add_argument("unk", help="positional input")
	args = parser.parse_args()

	df = pd.concat([
		pd.DataFrame(read_in_lenghts(args.mat, "mat")),
		pd.DataFrame(read_in_lenghts(args.pat, "pat")),
		pd.DataFrame(read_in_lenghts(args.unk, "unk"))],
		ignore_index=True)
	print("min_read_length\tmat\tmat_count\tpat\tpat_count\tunknown\tunknown_count")
	for l in range(0, 200000, 10000):
		phased_by_length(df, l)
