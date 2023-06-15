#!/usr/bin/env python3
import argparse
import os
import re
import shlex
import subprocess
import sys
import warnings
import pandas as pd

special_char = '[@!#$%^&*()<>?|}{~:];:`"\' '

# functions ---
def create_parser():

	parser = argparse.ArgumentParser(
		description = "Merge libraries from Illumina sequencing runs.",
		epilog = "Author: Andrea Lauria",
	)

	parser.add_argument(
		"-i",
		"--config",
		dest = "config",
		action = "store",
		default = None,
		help = "Path to src/dest configuration file (tab-delimited). N.B.: It must contain full paths.",
	)

	parser.add_argument(
		"-r1",
		"--run_1",
		dest = "run_1",
		action = "store",
		default = None,
		help = "Path to run 1 (Works if config is missing, merges by filenames).",
	)
	
	parser.add_argument(
		"-r2",
		"--run_2",
		dest = "run_2",
		action = "store",
		default = None,
		help = "Path to run 2 (Works if config is missing, merges by filenames).",
	)
	
	parser.add_argument(
		"-m",
		"--merged",
		dest = "merged",
		action = "store",
		default = None,
		help = "Path to merged run (Works if config is missing, merge by filenames).",
	)

	parser.add_argument(
		"--dry",
		dest="dry_run",
		action="store_true",
		help="Enable dry run.",
	)

	parser.add_argument(
		"--soft-link",
		dest="soft_link",
		action="store_true",
		help="Create soft links for samples present only once (Works if config is provided).",
	)

	return parser

def config_merge(path_config, dry_run = False, soft_link = False):
	
	config_df = pd.read_csv(path_config, sep='\t')

	# validate config
	validate_config = []
	for index, row in config_df.iterrows():
		validate_config.append(any([i in special_char for i in row['SRC']]) or any([i in special_char for i in row['DEST']]))

	if any(validate_config):
		print('\n[!] Special characters not allowed.')
		#print('\nAny of {} was found.'.format(shlex.quote(special_char)))
		print('\nCheck offending rows in configuration file:')
		print(config_df[validate_config])
	else:
		for index, row in config_df.iterrows():
			if not os.path.isfile(row['SRC']):
				print('[!] SRC file {} does not exist.'.format(row['SRC']))
				pass
			else:
				if index == 0 or config_df.loc[index]['DEST'] != config_df.loc[index-1]['DEST']:
					if soft_link == True and config_df['DEST'].value_counts()[row['DEST']]==1:
						cmd = 'ln -sf {} {}'.format(shlex.quote(row['SRC']), shlex.quote(row['DEST']))
					else:
						cmd = 'cat {} > {}'.format(shlex.quote(row['SRC']), shlex.quote(row['DEST']))
				else:
					cmd = 'cat {} >> {}'.format(shlex.quote(row['SRC']), shlex.quote(row['DEST']))
				print(cmd)
				if dry_run == False:
					subprocess.run(cmd, shell=True, check=True)

def run_merge(run1_dir, run2_dir, mergedir, dry_run = False, soft_link = False):
	
	if run1_dir == None or run2_dir == None or mergedir == None:
		print("\n[!] Missing required arguments.\n")
		parser.print_help()
		sys.exit()
	else:
		if not os.path.exists(mergedir):
			os.makedirs(mergedir)

		run1_files = [f for f in os.listdir(run1_dir) if re.match(r'.*\.(fastq|fq)(\.gz|)$', f)]
		run2_files = [f for f in os.listdir(run2_dir) if re.match(r'.*\.(fastq|fq)(\.gz|)$', f)]

		files = set(run1_files + run2_files)

		for f in files:
			merged_file = '{}/{}'.format(mergedir,f)
			f1 = '{}/{}'.format(run1_dir,f)
			if os.path.isfile(f1):
				cmd1 = "cat " + shlex.quote(f1)
			else:
				cmd1 = ""
			f2 = '{}/{}'.format(run2_dir,f)
			if os.path.isfile(f2):
				if cmd1 == "":
					cmd2 = "cat " + shlex.quote(f2)
				else:
					cmd2 = " " + shlex.quote(f2)
			else:
				cmd2 = ""

			if soft_link == True and (cmd1 == "" or cmd2 == ""):
				warnings.warn("\nSoft links might break with relative paths. Use config file.\n")
				cmd = cmd1 + cmd2 + " " + shlex.quote(merged_file)
				cmd = re.sub("cat","ln -sf",cmd)
			else:
				cmd = cmd1 + cmd2 + " > " + shlex.quote(merged_file)
			print(cmd)
			if dry_run == False:
				subprocess.run(cmd, shell=True, check=True)


def main(path_config, run1_dir, run2_dir, mergedir, dry_run, soft_link):
	print("\n[*] Merge Illumina Runs [*]\n")
	if path_config != None:
		config_merge(path_config, dry_run, soft_link)
	else:
		run_merge(run1_dir, run2_dir, mergedir, dry_run, soft_link)

# execution ---
if __name__ == '__main__':
	parser = create_parser()

	if len(sys.argv) == 1:
		parser.print_help()
		sys.exit()

	args = parser.parse_args()

	main(
		path_config = args.config,
		run1_dir = args.run_1,
		run2_dir = args.run_2,
		mergedir = args.merged,
		dry_run = args.dry_run,
		soft_link = args.soft_link,
	)
