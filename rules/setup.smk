# needs to be run first before other rules can be run.
rule download_genomes:
        output:
                checkpoint = "results/checkpoints/download_genomes.done",
		download_overview = "results/downloaded_genomes/download_overview.txt",
		success = "results/downloaded_genomes/successfully_downloaded.txt"
        singularity:
                "docker://reslp/biomartr:0.9.2"
        params:
                species = get_species_names,
                wd = os.getcwd()
        shell:
                """
                Rscript bin/genome_download.R {params.species} {params.wd}
                touch {output}
                """

rule rename_assemblies:
	input:
		rules.download_genomes.output.success
	output:
		checkpoint = "results/checkpoints/rename_assemblies.done"
	params:
		#downloaded_species = get_species_names_rename,
		local_species = get_local_species_names_rename,
		wd = os.getcwd()
	shell:
		"""
		#have to first remove this folder
		#rm -rf results/assemblies
		mkdir -p results/assemblies
		for spe in $(cat {input}); do
			if [[ -f {params.wd}/results/assemblies/"$spe".fna ]]; then
				continue
			else
				if [[ ! -f {params.wd}/results/downloaded_genomes/"$spe"_genomic_genbank.fna ]]; then
					echo "Species not found: $spe. Maybe it was not downloaded."
					continue
				else
					ln -s {params.wd}/results/downloaded_genomes/"$spe"_genomic_genbank.fna {params.wd}/results/assemblies/"$spe".fna
				fi
			fi
		done	
		for spe in  {params.local_species}; do
			sparr=(${{spe//,/ }})
			if [[ -f {params.wd}/results/assemblies/"${{sparr[0]}}".fna ]]; then
				continue
			else
				ln -s {params.wd}/"${{sparr[1]}}" {params.wd}/results/assemblies/"${{sparr[0]}}".fna
			fi
		done
		touch {output.checkpoint}
		"""



rule download_busco_set:
        output:
                busco_set = directory("results/busco_set"),
                checkpoint = "results/checkpoints/download_busco_set.done"
        params:
                set = config["busco"]["set"]
        shell:
                """
                echo {params.set}
                wget http://busco.ezlab.org/v2/datasets/{params.set}.tar.gz
                tar xfz {params.set}.tar.gz
                rm {params.set}.tar.gz
                if [ -d {output.busco_set} ]; then rm -rf {output.busco_set}; fi
                mv {params.set} {output.busco_set}
                touch {output.checkpoint}
                """

rule prepare_augustus:
        output:
                augustus_config_path = directory("results/augustus_config_path"),
                checkpoint = "results/checkpoints/prepare_augustus.done"
        singularity:
                "docker://reslp/busco:3.0.2"
        shell:
                """
                cp -r /opt/conda/config {output.augustus_config_path}
                touch {output.checkpoint}
                """
