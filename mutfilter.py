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
		mdtags.update({read.query_name: (read.query_sequence, md)})

	bamfile.close()

	return mdtags


def get_mismatch(tagsv):
	'''
	# Retrieve mismatches from MD tag
	# Returns a tuple containing:
	- list of mismatches in the form ref=>alt
	- MD tag
	- Sequence in read
	'''
	md  = tagsv[1]
	seq = tagsv[0]

	dd = list()
	nt = list()

	for match in re.finditer(r"\d+",md):
		dd.append((int(match.group()), match.start()))

	for match in re.finditer(r"[A,T,C,G]",md):
	    nt.append((match.group(), match.start()))

	mismatch = list()

	pos = 0
	for i in range(0,len(nt)):
		pos += dd[i][0]
		mismatch.append(nt[i][0] + '=>' + seq[pos])
		pos += 1
			
	return (mismatch, md, seq)


def filter_mismatch(mismatch, mutation, report_exact_matches):
	'''
	# Filter mismatch for specific mutations
	'''
	
	if report_exact_matches:
		filtmismatch = all(x == mutation for x in mismatch[0]) or len(mismatch[0])==0
	else:
		filtmismatch = all(x == mutation for x in mismatch[0]) and len(mismatch[0])!=0

	return filtmismatch


def filter_bam(bampath, mutation, report_exact_matches):
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
		md  = read.get_tag("MD")
		seq = read.query_sequence

		mismatch = get_mismatch((seq, md))
		
		if filter_mismatch(mismatch, mutation, report_exact_matches):
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
		, help='Mutation to filter BAM with. Default: None')

	parser.add_argument('-s','--stdout'
		, dest='stdout'
		, action='store_true'
		, default=None
		, help='If set, report mutations to stdout instead of writing to BAM.')

	parser.add_argument('-e','--exact'
		, dest='report_exact_matches'
		, action='store_true'
		, default=None
		, help='Report reads with no mismatches.')

	return parser


def main(bampath, mutation, stdout, report_exact_matches):
	'''
	Main program function
	'''
	print("\n[*] Mutation filter [*]\n")

	if stdout:
		# Print mutations to stdout
		mdtags   = parse_mdtag(bampath)
		mismatch = dict((k, get_mismatch(v)) for k,v in mdtags.items())

		for k,v in mismatch.items():
			print(k, " :: ", v)

	else:
		# Write filtered BAM
		filter_bam(bampath, mutation, report_exact_matches)

	sys.exit("\n[*] done.\n")


if __name__ == '__main__':

	parser = create_parser()

	if len(sys.argv) == 1:
		parser.print_help()
		sys.exit()

	# Parse arguments	
	args = parser.parse_args()

	# Execution
	main(bampath 				= args.bampath
		, mutation              = args.mutation
		, stdout 				= args.stdout
		, report_exact_matches  = args.report_exact_matches)
