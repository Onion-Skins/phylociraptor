#!/usr/bin/env python
# written by Philipp Resl
# big revision by Rodger Wang in July 2023 to
# - Adapt to new Busco version(works with BUSCO 5.3.2)
# - Add a new way to eliminate the need of single_copy_sequence.txt
# - Add a new option to specify desired sequence file type: faa or fna
# - Add a new option to specify the location of squence files location under busco result directory

import os
import sys
import argparse


if sys.version_info[0] < 3:
    raise Exception("Must be using Python 3")

pars = argparse.ArgumentParser(prog="extract_busco_table.py", description = """This script will get all busco3 hmms and look in the busco results all specified genomes for the presence of the files""", epilog = """written by Philipp Resl""")
pars.add_argument('--hmm', dest="hmm_dir", required=True, help="Directory of the BUSCO hmms.")
pars.add_argument('--busco_results', dest="busco_results", required=True, help="Results directory containing all BUSCO runs.")
pars.add_argument('-o', dest="out", required=True, help="BUSCO table output file.")
pars.add_argument('--type' , dest="type", required=True, help="Type of sequences (aa or nu).")
pars.add_argument('--genome_busco_sequence_subdir', dest="genome_busco_sequence_subdir", default="", type=str, help="""sub-directory path under BUSCO result dirrectory which contains all genome BUSCO sequence files, defaults to '', but could be orthology/busco/busco_runs.""")
pars.add_argument('--method', dest="method", default="no-tar-file", type=str, help="""method used to create sequence files. Could be 'use-tar-file' or 'no-tar-file' (default).""")
args=pars.parse_args()

extension=""
if args.type == "nu":
	extension = ".fna"
else:
	extension = ".faa" 

hmms = os.listdir(args.hmm_dir)
hmms = [filename.strip(".hmm") for filename in hmms]

genomes = os.listdir(args.busco_results)
genomes = [ f for f in genomes if 'BUSCO_' in f ]      # filter out any directory/file not genome

with open(args.out, "w") as outfile:
	header = "species\t"
	header += "\t".join(hmms)
	header += "\tpercent_complete" 
	print(header, file= outfile)
	for species in genomes:
		ones = 0
		zeros = 0
		outstring = species
		print("Extracting HMMs for", species, file=sys.stderr)
		if args.method == "use-tar-file":
			print('Use precompiled single copy sequence list file:')
			file_list = []
			file_path = os.path.join(args.busco_results, args.genome_busco_sequence_subdir, species, "run_busco/single_copy_busco_sequences.txt")
			#print(file_path)
			with open(file_path, "r") as file_list_file:
				 for line in file_list_file:
				 	line = line.strip()
				 	filename = line.split(' ')[-1]
				 	file_list.append(filename)
		else:  # default method - no tar file
			dir_path = os.path.join(args.busco_results, species, args.genome_busco_sequence_subdir, "run_eutheria_odb10/busco_sequences/single_copy_busco_sequences")
			file_list = os.listdir(dir_path)
		
		# use specified extension (faa or fna) file only. This should be enough. Strip file suffix at the same time
		buscos = [filename.strip(extension) for filename in file_list if extension in filename]
		#print('buscos', species + ':', buscos)
		try:  
			for hmm in hmms:
				if hmm in buscos:
					outstring += "\t"
					outstring += "1"
					ones +=1
				else:
					outstring += "\t"
					outstring += "0"
					zeros +=1
			percent = ones / (ones+zeros)
			outstring += "\t"
			outstring += str(percent)
			print(outstring, file=outfile)
		except:
			out = species + " not found. Skipped.\n"
			print(out, file=sys.stderr)
