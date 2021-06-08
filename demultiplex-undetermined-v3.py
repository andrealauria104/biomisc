#!/usr/bin/env python3
############################
# demultiplex-undetermined #
############################
"""
Demultiplexing Undetermined_S0_R1_001.fastq.gz files produced by bcl2fastq.
Reports reads with up to one mismatch in barcode sequence.
"""
import argparse, gzip, os, re, regex, sys, time
import pandas as pd
from multiprocessing import Pool

# Functions ----
def sub(text, label):
    SUBSTITUTIONS = {
    "A": "TCGN",
    "T": "ACGN",
    "C": "ATGN",
    "G": "ATCN"
    }
    possibilities = [c + SUBSTITUTIONS[c] for c in text]
    words = []
    for i in range(0,len(possibilities)):
        for j in range(0,len(possibilities[i])):
            words.append(text[:i] + possibilities[i][j] + text[i + 1:])
    return dict(zip(words,[label]*len(words)))

def create_parser():
	parser = argparse.ArgumentParser(
		description = "Demultiplexing Undetermined_S0_R1_001.fastq.gz files produced by bcl2fastq."
		, epilog = "Author: Andrea Lauria")
	
	parser.add_argument("-s","--sample-sheet"
		, dest = "sample_sheet"
		, action = "store"
		, default = None
		, required=True
		, help = "Path to SampleSheet.")
	
	parser.add_argument("-u","--undetermined"
		, dest = "undetermined"
		, action = "store"
		, default = None
		, required=True
		, help = "Path to Undetermined_S0_R1_001.fastq.gz.")
	
	parser.add_argument("-m","--mismatch"
		, dest = "mismatch"
		, action = "store"
		, default = 1
		, type=int
		, help = "Max allowed mismatches in barcode [Default = 1].")

	parser.add_argument("-p","--n-processes"
		, dest = "nprocess"
		, action = "store"
		, default = 8
		, type=int
		, help = "Number of processes [Default = 8].")

	parser.add_argument("-l","--log-report"
		, dest = "log_report"
		, action = "store"
		, default = "demultiplex-undetermined_report.txt"
		, help = "Demultiplexing report file. [Default = demultiplex-undetermined_report.txt]")

	return parser

# 0. Parse arguments ----
parser = create_parser()

if len(sys.argv) == 1:
	parser.print_help()
	sys.exit()

args = parser.parse_args()

path_samplesheet = args.sample_sheet
path_undetermined = args.undetermined
mismatch = args.mismatch
nprocess = args.nprocess
path_log_report = args.log_report

# 1. Read sample sheet ----
samplesheet = pd.read_csv(path_samplesheet, sep="\t")
if "index2" in samplesheet.columns:
    samplesheet["barcode"] = samplesheet["index"] + samplesheet["index2"]
else:
    samplesheet["barcode"] = samplesheet["index"]
     
samplesheet = samplesheet[["Sample_Name","barcode"]]

# 2. Define dict/hash for barcodes ----
samplesheet_dict = dict(zip(samplesheet["Sample_Name"],samplesheet["barcode"]))
countreads_dict = dict.fromkeys(samplesheet["Sample_Name"],0)
barcode_dict = {}
for k,v in samplesheet_dict.items():
    barcode_dict.update(sub(v,k))

# 3. Process unaligned ----
start_time = time.time()
undetermined = gzip.open(path_undetermined,'rt')
countline = 0
countreads = 0
while undetermined:
    r_id  = undetermined.readline()
    # break loop at EOF
    if not r_id:
        break
    r_seq  = undetermined.readline()
    r_plus = undetermined.readline()
    r_qual = undetermined.readline()

    barcode = re.sub(".*\\:|\n|\\+","",r_id)
    try:
        outfile = barcode_dict[barcode] + "_R1.fastq.gz"
        with gzip.open(outfile, 'at') as out:
            out.write(r_id)
            out.write(r_seq)
            out.write(r_plus)
            out.write(r_qual)
            out.close()
            countreads_dict[barcode_dict[barcode]] = countreads_dict[barcode_dict[barcode]] + 1
    except:
        pass
       
    countline = countline+4
    if countline % 4000 == 0:
        print(str(int(countline/4)) + " reads processed")

undetermined.close()
print("--- %s seconds ---" % (time.time() - start_time))

print(" -- writing report: " + path_log_report)
with open(path_log_report, 'w') as out_report:
	out_report.write("Sample\tBarcode\tNreads\n")
	for k in samplesheet_dict.keys():
		out_report.write(k + "\t" + samplesheet_dict[k] + "\t" + str(countreads_dict[k]) + "\n")
out_report.close()
