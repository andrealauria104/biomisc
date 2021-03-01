#!/usr/bin/env python3
import argparse, sys, os,re

# functions ---
def create_parser():
	parser = argparse.ArgumentParser(
		description = "Merge libraries from Illumina runs."
		, epilog = "Author: Andrea Lauria")
	
	parser.add_argument("-r1","--run_1"
		, dest = "run1_dir"
		, action = "store"
		, default = None
		, required=True
		, help = "Path to run 1 directory.")
	
	parser.add_argument("-r2","--run_2"
		, dest = "run2_dir"
		, action = "store"
		, default = None
		, required=True
		, help = "Path to run 2 directory.")
	
	parser.add_argument("-m","--merged"
		, dest = "mergedir"
		, action = "store"
		, default = None
		, required=True
		, help = "Path to merged run directory.")

	parser.add_argument("-l","--link"
		, dest = "link"
		, action='store_true'
		, help='Enable symlink if not merging (N.B.: be careful on relative paths! If you are not sure, use full directory paths).')

	parser.add_argument("-n","--dry"
		, dest = "dry"
		, action='store_true'
		, help='Enable dry run.')

	return parser

# read data ---
if __name__ == '__main__':
	
	parser = create_parser()

	if len(sys.argv) == 1:
		parser.print_help()
		sys.exit()

	args = parser.parse_args()

	run1_dir = args.run1_dir
	run2_dir = args.run2_dir
	mergedir = args.mergedir
	link     = args.link
	dry      = args.dry

	if not os.path.exists(mergedir):
		os.makedirs(mergedir)

	run1_files = [f for f in os.listdir(run1_dir) if re.match(r'.*\.(fastq|fq)(\.gz|)$', f)]
	run2_files = [f for f in os.listdir(run2_dir) if re.match(r'.*\.(fastq|fq)(\.gz|)$', f)]

	files = set(run1_files + run2_files)

	for f in files:
		merged_file = mergedir + "/" + f
		f1 = run1_dir + "/" + f
		if os.path.isfile(f1):
			cmd1 = "cat " + f1
		else:
			cmd1 = ""
		f2 = run2_dir + "/" + f
		if os.path.isfile(f2):
			if cmd1 == "":
				cmd2 = "cat " + f2
			else:
				cmd2 = " " + f2
		else:
			cmd2 = ""

		if link == True and (cmd1 == "" or cmd2 == ""):
			cmd = "ln -s " + cmd1 + cmd2 + " " + merged_file
			cmd = re.sub("cat ","",cmd)
		else:
			cmd = cmd1 + cmd2 + " > " + merged_file
		print(cmd)
		if dry == False:
			os.system(cmd)
