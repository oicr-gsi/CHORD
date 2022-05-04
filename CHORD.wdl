version 1.0

workflow CHORD {
	input {
		File smallsVcfFile
		File smallsVcfFileIndex
		File structuralVcfFile
		String basename = basename("~{structuralVcfFile}", ".vcf.gz")
		String rScript
	}

	call filterSmalls {
		input:
		smallsVcfFile = smallsVcfFile,
		smallsVcfFileIndex = smallsVcfFileIndex
	}

	call filterStructural {
		input: 
		structuralVcfFile = structuralVcfFile
	}

	call hrdResults {
		input:	
		smallsVcffiltered = filterSmalls.smallsVcffiltered,
		structuralbedpe = filterStructural.structuralbedpe,
		rScript = rScript,
		basename = basename
	}

	parameter_meta {
		smallsVcfFile: "Input VCF file of small mutations, eg from mutect2"
		smallsVcfFileIndex: "Index for VCF file of small mutations"
		structuralVcfFile: "Input VCF file"
		rScript: "Temporary variable to call the .R script containing CHORD, will be modulated. default: ~/CHORD/CHORD_results.R"
	}

	meta {
		author: "Felix Beaudry, Alex Fortuna"
		email: "fbeaudry@oicr.on.ca"
		description: "Homolog Recombination Deficiency Prediction Workflow"
		dependencies: 
		[
			{
				name: "gatk/4.2.0.0",
				url: "https://github.com/broadinstitute/gatk/releases"
			},
			{
				name: "bcftools/1.9",
				url: "https://samtools.github.io/bcftools/"
			},
			{
				name: "bedtools/2.27.1",
				url: "https://bedtools.readthedocs.io/en/latest/content/installation.html"
			},
			{
				name: "R/4.1.2",
				url: "https://cran.r-project.org/mirrors.html"
			}
		]
		output_meta: {
			hrdOutput: "HRD score and confidence as .txt file as estimated by CHORD"
		}
	}

	output {
		File hrdOutput = "~{basename}.CHORD.hrd.txt"
	}
}

task filterStructural {
	input {
		File structuralVcfFile 
		String basename = basename("~{structuralVcfFile}", ".vcf.gz")
		String modules = "bcftools/1.9"
		Int jobMemory = 5
		Int threads = 1
		Int timeout = 1
	}

	parameter_meta {
		structuralVcfFile: "Vcf input file"
		basename: "Base name"
		modules: "Required environment modules"
		jobMemory: "Memory allocated for this job (GB)"
		threads: "Requested CPU threads"
		timeout: "Hours before task timeout"
	}

	command <<<
		set -euo pipefail

		$BCFTOOLS_ROOT/bin/bcftools view -f 'PASS' ~{structuralVcfFile} |\
		$BCFTOOLS_ROOT/bin/bcftools filter -e 'INFO/SVTYPE = "BND"' |\
		$BCFTOOLS_ROOT/bin/bcftools query -f "%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\n" |\
		awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$3-$2} ' >~{basename}.lengths.bed 

	>>>

	runtime {
		modules: "~{modules}"
		memory:  "~{jobMemory} GB"
		cpu:     "~{threads}"
		timeout: "~{timeout}"
	}

	output {
		File structuralbedpe = "~{basename}.lengths.bed"
	}

	meta {
		output_meta: {
			structuralbedpe: "filtered structural .bedpe with lengths",
		}
	}
}

task filterSmalls {
	input {
		File smallsVcfFile
		File smallsVcfFileIndex
		String smallsFilter = "'PASS,clustered_events,slippage'"
		String basename = basename("~{smallsVcfFile}", ".vcf.gz")
		String modules = "gatk/4.2.0.0 tabix/1.9 bcftools/1.9 hg38/p12 grch38-alldifficultregions/3.0"
		String genome = "$HG38_ROOT/hg38_random.fa"
		String? difficultRegions
		String VAF
		Int jobMemory = 10
		Int threads = 1
		Int timeout = 2
	}

	parameter_meta {
		smallsVcfFile: "Vcf input file"
		basename: "Base name"
		modules: "Required environment modules"
		genome: "Path to loaded genome"
		difficultRegions: "Path to .bed of difficult regions to align to, string must include the --exclude-intervals flag, eg: --exclude-intervals $GRCH38_ALLDIFFICULTREGIONS_ROOT/GRCh38_alldifficultregions.bed"
		VAF: "VAF for indels"
		jobMemory: "Memory allocated for this job (GB)"
		threads: "Requested CPU threads"
		timeout: "Hours before task timeout"
	}

	command <<<
		set -euo pipefail

		gatk SelectVariants \
		-V ~{smallsVcfFile} \
		-R ~{genome} ~{difficultRegions} \
		-O ~{basename}.vcf  

		$BCFTOOLS_ROOT/bin/bcftools view -f ~{smallsFilter} ~{basename}.vcf  |  $BCFTOOLS_ROOT/bin/bcftools filter -i "(FORMAT/AD[0:1])/(FORMAT/AD[0:0]+FORMAT/AD[0:1]) >= 0.~{VAF}" >~{basename}.VAF.vcf

		bgzip ~{basename}.VAF.vcf
		tabix -p vcf ~{basename}.VAF.vcf.gz

		zcat ~{smallsVcfFile} | awk '$1 !~ "#" {print}'  | wc -l >~{basename}.filteringReport.txt
		awk '$1 !~ "#" {print}' ~{basename}.vcf | wc -l >>~{basename}.filteringReport.txt
		zcat ~{basename}.VAF.vcf.gz | awk '$1 !~ "#" {print}'  | wc -l >>~{basename}.filteringReport.txt

	>>> 

	runtime {
		modules: "~{modules}"
		memory: "~{jobMemory} GB"
		cpu: "~{threads}"
		timeout: "~{timeout}"
	}

	output {
		File smallsVcffiltered = "~{basename}.VAF.vcf.gz"
		File smallsVcfIndexfiltered = "~{basename}.VAF.vcf.gz.tbi"
		File smallsFilteringReport = "~{basename}.filteringReport.txt"
	}

	meta {
		output_meta: {
			indelVcfOutput: "filtered INDEL .vcf",
			indelVcfIndexOutput: "filtered INDEL .vcf.tbi indexed",
			indelFilteringReport: "counts of variants pre and post filtering"
		}
	}
}

task hrdResults {
	input {
		File smallsVcffiltered
		File structuralbedpe
		String basename 
		String rScript
		String modules = "chord/2.0"
		Int jobMemory = 32
		Int threads = 4
		Int timeout = 16
	}

	parameter_meta {
		smallsVcffiltered: "filtered SNV Vcf input file"
		structuralbedpe: "filtered structural variant bedpe input file"
		modules: "Required environment modules"
		jobMemory: "Memory allocated for this job (GB)"
		threads: "Requested CPU threads"
		timeout: "Hours before task timeout"
	}

	command <<<
		set -euo pipefail

		Rscript --vanilla ~{rScript} ~{smallsVcffiltered} ~{structuralbedpe} ~{basename}

	>>> 

	runtime {
		modules: "~{modules}"
		memory:  "~{jobMemory} GB"
		cpu:     "~{threads}"
		timeout: "~{timeout}"
	}

	output {
		File hrdOutput = "~{basename}.CHORD.hrd.txt"
	}

	meta {
		output_meta: {
			hrdOutput: "HRD score and confidence as .txt file as estimated by CHORD"
		}
	}
}