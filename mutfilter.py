#!/usr/bin/env python3
'''
##
# BAM mutation filter tool.
##
# Author: Andrea Lauria, 2019
##
# Born as: Filter for PARCLIP aligments: allow only T => C mutations
##
'''
import pysam, re, subprocess, argparse, sys

# Global ---
alphabet = {"A":"T","T":"A","C":"G","G":"C"}

# Function definitions ---
def parse_mdtag(bampath):
	'''
	# Read BAM file and parse MD tag
	'''
	print(" --- reading ", bampath)
	bamfile = pysam.AlignmentFile(bampath, "rb")
	mdtags  = dict()

	for read in bamfile.fetch():
		if read.is_unmapped:
			continue
		md = read.get_tag("MD")
		mdtags.update({read.query_name: {'md':md, 'seq':read.query_sequence, 'reverse':read.is_reverse}})

	bamfile.close()

	return mdtags


def get_mismatch(tagsv):
	'''
	# Retrieve mismatches from MD tag
	# Returns a tuple containing:
	- list of mismatches in the form ref=>alt
	- MD tag
	- Sequence in read
	- Is reverse-complemented
	'''

	dd = list()
	nt = list()

	for match in re.finditer(r"\d+",tagsv['md']):
		dd.append((int(match.group()), match.start()))

	for match in re.finditer(r"[A,T,C,G]",tagsv['md']):
	    nt.append((match.group(), match.start()))

	mismatch = list()

	pos = 0
	for i in range(0,len(nt)):
		pos += dd[i][0]
		mismatch.append(nt[i][0] + '=>' + tagsv['seq'][pos])
		pos += 1
			
	return (mismatch, tagsv)


def filter_mismatch(mismatch, mutation, report_exact_matches, directional):
	'''
	# Filter mismatch for specific mutations
	'''
	reverse_mutation = alphabet[mutation.split('=>')[0]] + '=>' + alphabet[mutation.split('=>')[1]]

	if report_exact_matches == True:
		if directional == True:
			if mismatch[1]['reverse'] == True:
				filtmismatch = all(x == reverse_mutation for x in mismatch[0]) or len(mismatch[0])==0
			else:
				filtmismatch = all(x == mutation for x in mismatch[0]) or len(mismatch[0])==0
		else:
			filtmismatch = all(x == mutation for x in mismatch[0]) or len(mismatch[0])==0
	else:
		if directional == True:
			if mismatch[1]['reverse'] == True:
				filtmismatch = all(x == reverse_mutation for x in mismatch[0]) and len(mismatch[0])!=0
			else:
				filtmismatch = all(x == mutation for x in mismatch[0]) and len(mismatch[0])!=0
		else:
			filtmismatch = all(x == mutation for x in mismatch[0]) and len(mismatch[0])!=0

	return filtmismatch


def filter_bam(bampath, mutation, report_exact_matches, directional):
	'''
	# Filter and write to BAM file
	# Sort and index BAM
	'''
	filtbampath = re.sub(r'\.bam$','_tmp.filt.bam', bampath)
	sortbampath = re.sub(r'_tmp','', filtbampath)
	
	print(" --- reading ", bampath)
	print(" --- writing ", sortbampath)

	bamfile     = pysam.AlignmentFile(bampath, "rb")
	filtbamfile = pysam.AlignmentFile(filtbampath, "wb", template = bamfile)

	for read in bamfile.fetch():
		if read.is_unmapped:
			continue
		# create input for mismatch retrieval: MD tag, sequence, strandness
		tagsv = {"md": read.get_tag("MD"), "seq": read.query_sequence, "reverse": read.is_reverse}
		
		mismatch = get_mismatch(tagsv)
		
		if filter_mismatch(mismatch, mutation, report_exact_matches, directional):
			filtbamfile.write(read)

	filtbamfile.close()
	bamfile.close()

	pysam.sort("-o",sortbampath, filtbampath)
	pysam.index(sortbampath)
	subprocess.run(["rm",filtbampath])


def create_parser():
	'''
	Command line argument parser
	'''
	parser = argparse.ArgumentParser(
		description='BAM mutation filter tool.'
	   , epilog = 'Author: Andrea Lauria'
		)
	
	parser.add_argument('bampath'
		, metavar = '</path/to/bam>'
		, action = 'store'
		, help = 'Path to BAM file')
	
	parser.add_argument('-m','--mutation'
		, dest='mutation'
		, action='store'
		, default='T=>C'
		, help='Mutation to filter BAM with. Default: T=>C')

	parser.add_argument('-s','--stdout'
		, dest='stdout'
		, action='store_true'
		, default=None
		, help='If set, report mutations to stdout instead of writing to BAM.')

	parser.add_argument('-e','--exact'
		, dest='report_exact_matches'
		, action='store_true'
		, default=None
		, help='Include reads with no mutations.')

	parser.add_argument('-d','--directional'
		, dest='directional'
		, action='store'
		, default=True
		, help='Assume directional (stranded) library (Default: True).')

	return parser


def main(bampath, mutation, stdout, report_exact_matches, directional):
	'''
	Main program function
	'''
	print("\n[*] Mutation filter [*]\n")

	if stdout:
		# Print mutations to stdout
		mdtags   = parse_mdtag(bampath)
		mismatch = dict((k, get_mismatch(v)) for k,v in mdtags.items())
		print(report_exact_matches)
		print(mutation)
		for k,v in mismatch.items():
			if filter_mismatch(v, mutation, report_exact_matches, directional):
				print(k, " :: ", v)
			else:
				continue
	else:
		# Write filtered BAM
		filter_bam(bampath, mutation, report_exact_matches, directional)

	sys.exit("\n[*] done.\n")


if __name__ == '__main__':

	parser = create_parser()

	if len(sys.argv) == 1:
		parser.print_help()
		sys.exit()

	# Parse arguments	
	args = parser.parse_args()

	# Execution
	main(bampath = args.bampath
		, mutation = args.mutation
		, stdout = args.stdout
		, report_exact_matches = args.report_exact_matches
		, directional = args.directional)
