#!/usr/bin/env python
# First written by Philipp Resl
# big revision by Rodger Wang in July 2023 to  
# - adapt to new Busco version (works with BUSCO 5.3.2);
# - fixed a bug by in Philipp's code which may cause duplicated sequences when feeding data to iq-tree later
# - vastly improve speed in the order of hundred times by untarring tar file in advance
# - Add a new way to eliminate the need of tar file
# - Add a new option to specify the location of sub-directory containing squence files under busco result directory

import os
import sys
import pandas as pd
from Bio import SeqIO
import argparse
import tarfile
from io import StringIO
from io import TextIOWrapper
import shutil

if sys.version_info[0] < 3:
	raise Exception("Must be using Python 3")

pars = argparse.ArgumentParser(prog="create_sequence_files.py", description = """This script will create fasta files for all the buscos from all species with >% of buscos present""", epilog = """First written by Philipp Resl, then revised by Rodger Wang""")
pars.add_argument('--busco_table', dest="busco_table", required=True, help="Path to BUSCO table.")
pars.add_argument('--busco_results', dest="busco_results", required=True, help="Results directory containing all BUSCO genomes.")
pars.add_argument('--cutoff', dest="cutoff", default=0, required=True, help="Percent cutoff %% for BUSCO presence. Species below this threshold will be excluded.")
pars.add_argument('--outdir', dest="outdir", required=True, help="Path to output directory.")
pars.add_argument('--minsp', dest="minsp", required=True, help="Minimum number of species which have to be present to keep the sequences.")
pars.add_argument('--type' , dest="type", required=True, help="Type of sequences (aa or nu).")
pars.add_argument('--genome_statistics' , dest="genome_stats", required=True, help="Path to genome statistics output file.")
pars.add_argument('--gene_statistics' , dest="gene_stats", required=True, help="Path to gene statistics output file.")
pars.add_argument('--batchID' , dest="batchid", default=0, type=int, help="Batch ID (start for subsampling)")
pars.add_argument('--nbatches', dest="nbatches", default=1, type=int, help="Total number of batches (step size for subsampling)")
pars.add_argument('--genome_busco_sequence_subdir', dest="genome_busco_sequence_subdir", default="", type=str, help="""sub-directory path under BUSCO result dirrectory which contains all genome BUSCO sequence files, defaults to '', but could be something like 'orthology/busco/busco_runs'.""")
pars.add_argument('--method', dest="method", default="no-tar-file", type=str, help="""method used to create sequence files. Could be 'use-tar-file' or 'no-tar-file' (default).""")


args=pars.parse_args()

extension=""
if args.type == "nu":
	extension = ".fna"
else:
	extension = ".faa" 
	
temp_dir = "temp"

busco_overview = pd.read_csv(args.busco_table, sep="\t")
#print(busco_overview)
print("Settings:")
print("cutoff: ", args.cutoff)
print("minsp: ", args.minsp)
print("type: ", args.type)
print("outdir: ", args.outdir)
print("batchID: %i / %i" %(args.batchid, args.nbatches))

species_list = busco_overview.species.tolist()
print("Original number of species:", len(species_list))
#print(species_list)
#first remove species with too low busco coverage
busco_overview = busco_overview.set_index("species")

with open(args.genome_stats, "w") as genome_file:
	for sp in species_list:
		if busco_overview.loc[sp, "percent_complete"] < float(args.cutoff):
			out = sp + "\tFAILED" + "\t" + str(busco_overview.loc[sp, "percent_complete"]) + "\t" + str(float(args.cutoff))
			print(out, file=genome_file)
			busco_overview = busco_overview.drop([sp])
		else:
			out = sp + "\tOK" + "\t" + str(busco_overview.loc[sp, "percent_complete"]) + "\t" + str(float(args.cutoff)) 
			print(out, file=genome_file)
	species_list =  list(busco_overview.index)
	print("Species remaining after applying cutoff:", len(species_list))

#now loop through each busco and extract sequence for each species
buscos = list(busco_overview.columns.values)
buscos.remove("percent_complete")

gene_file = open(args.gene_stats, "w").close()   # clear the gene_stats file if already exists

with open(args.gene_stats, "a") as gene_file:
	# first get gene file lists of all genomes
	genome_gene_file_lists = {}
	if args.method == "use-tar-file":
			# remove old temp dir and recreate a empty one
			if os.path.isdir(temp_dir):
				print("Temp dir exists. Empty it.")
				shutil.rmtree(temp_dir)
			os.mkdir(temp_dir)

	for species in species_list:
		print('create file list for', species+'.')
		if args.method == "use-tar-file":
			file_list = []
			file_path = os.path.join(args.busco_results, args.genome_busco_sequence_subdir, species, "run_busco/single_copy_busco_sequences.txt")
			#print(file_path)
			with open(file_path, "r") as file_list_file:
				 for line in file_list_file:
				 	line = line.strip()
				 	filename = line.split(' ')[-1]
				 	file_list.append(filename)
			
			# untar the tar file to temp dir at the same time
			print('Untarring file for', species+'.')
			species_dir = os.path.join(temp_dir, species)
			os.mkdir(species_dir)
			tf = tarfile.open(os.path.join(args.busco_results, args.genome_busco_sequence_subdir, species, "run_busco/single_copy_busco_sequences.tar"))
			tf.extractall(species_dir)     # extract tar file to temp folder
			tf.close()
			print('Done untarring.')
		else:   # default method - no tar file
			dir_path = os.path.join(args.busco_results, args.genome_busco_sequence_subdir, species, "run_eutheria_odb10/busco_sequences/single_copy_busco_sequences")
			file_list = os.listdir(dir_path)
		
		filtered_list = [f for f in file_list if extension in f]
		genome_gene_file_lists[species] = filtered_list
		print('Done file list creation.')
		
	# now we begin to loop through genes to filter based on minsp
	# and create fas files along the way
	i = args.batchid     # the start if the subsampling. defaults to 0
	while i < len(buscos):
		busco = buscos[i]
		print("Processing: " + busco)
		outstring = ""
		for species in species_list:
			filename = busco + extension
			if filename in genome_gene_file_lists[species]:
				if args.method == "use-tar-file":
					file_path = os.path.join(temp_dir, species, filename)
				else:  #defualt method - no tar file
					file_path = os.path.join(args.busco_results , args.genome_busco_sequence_subdir, species, "run_eutheria_odb10/busco_sequences/single_copy_busco_sequences", filename)
				if os.path.exists(file_path):
					for seq_record in SeqIO.parse(file_path, "fasta"):
						name = ">" +species+"\n"
						outstring += name
						outstring += (str(seq_record.seq) + "\n")
				else:
					print(busco, "not found for", species)
							
		if outstring.count(">") >= int(args.minsp):	# only keep sequences if total number is larger than specified cutoff above.		
			print(busco + "\t" + "OK" + "\t" + str(outstring.count(">")) +"\t" + str(int(args.minsp)), file=gene_file)
			with open(os.path.join(args.outdir, busco+"_all.fas"), "w") as outfile:
				outfile.write(outstring)
		else:
			print(busco + "\t" + "FAILED" + "\t" + str(outstring.count(">")) +"\t" + str(int(args.minsp)), file=gene_file)
		i += args.nbatches  # for subsampling. nbatches defaults to 1 which is no subsmpling
