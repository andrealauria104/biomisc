#!/usr/bin/env python3
import argparse, sys

# functions ---
def read_nreads(fnreads):
	d = {}
	with open(fnreads,"r") as f:
		for lines in f:
			key, value = lines.strip("\n").split("\t")
			d[key] = int(value)
	return d

def create_parser():
	parser = argparse.ArgumentParser(
		description = "Merge libraries from Illumina runs."
		, epilog = "Author: Andrea Lauria")
	
	parser.add_argument("-r1","--run_1"
		, dest = "run_1"
		, action = "store"
		, default = "nreads_1.txt"
		, help = "Path to run 1.")
	
	parser.add_argument("-r2","--run_2"
		, dest = "run_2"
		, action = "store"
		, default = "nreads_2.txt"
		, help = "Path to run 2.")
	
	parser.add_argument("-m","--merged"
		, dest = "merged"
		, action = "store"
		, default = "nreads_merged.txt"
		, help = "Path to merged run.")
	
	return parser
	
# read data ---
parser = create_parser()

if len(sys.argv) == 1:
	parser.print_help()
	sys.exit()

args = parser.parse_args()

r1 = read_nreads(args.run_1)
r2 = read_nreads(args.run_2)
m  = read_nreads(args.merged)

# check merging
print("\n -- Check merging ...")
print("\n")
check = {}
for k in m.keys():
	res = r1[k]+r2[k]
	check[k] = m[k] == r1[k]+r2[k]
	if check[k] is True:
		toprint = str(r1[k]) + " + " + str(r2[k]) + " = " + str(res) + " == "+ str(m[k])
		print(toprint)

checkres = all(check[k] == True for k in check.keys())
if checkres == True:
	print("\nn = " + str(len(m)))
	print("\n [+] Merging is correct!\n")
else:
	print("\n [!] Something is wrong, uncorrect merging. \n")
