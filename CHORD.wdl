workflow CHORD {

	input {
    	File snvVcfFile
    	File snvVcfIndex
    	File structuralVcfFile
    	File structuralVcfIndex
	}

	call filterSNVs {
		input: vcfFile = snvVcfFile
	}

	call filterStructural {
		input: vcfFile = structuralVcfFile
	}

	call hrdResults {
		input:	snvVcfFiltered = filterSNVs.snvVcfOutput
				structuralVcfFiltered = filterStructural.structuralVcfOutput
	}

	parameter_meta {
    	vcfFile: "Input VCF file"
    	vcfIndex: "Input VCF index file"
    	targetBed: "Target bed file"
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
   			
    	}
	}

	output {
		

  }
}

task filterSNVs {
 	input {
	    File vcfFile 
	    String basename = basename("~{vcfFile}", ".vcf.gz")
	    String ncbiBuild
	    String referenceFasta
	    String modules = "gatk/4.2.0.0 bcftools/1.9"
	    Int jobMemory = 32
	    Int threads = 4
	    Int timeout = 16
	 }

	parameter_meta {
	    vcfFile: "Vcf input file"
	    basename: "Base name"
	    ncbiBuild: "The assembly version"
	    referenceFasta: "Reference fasta file"
	    modules: "Required environment modules"
	    jobMemory: "Memory allocated for this job (GB)"
	    threads: "Requested CPU threads"
	    timeout: "Hours before task timeout"
	}

	command <<<
	    set -euo pipefail

		gatk SelectVariants -R ~{referenceFasta} \
			--exclude-intervals GRCh38_alldifficultregions.bed \
			-V ~{SNVvcfFile} \
			-sn ~{tumor}  -O ~{basename}.wellMapped.vcf

		bcftools filter -i '(FORMAT/AD[0:1]*100)/(FORMAT/AD[0:0]+FORMAT/AD[0:1]) >= 10'\
			~{basename}.wellMapped.vcf >~{basename}.wellMapped.MAF.vcf

	>>> 

	runtime {
	    modules: "~{modules}"
	    memory:  "~{jobMemory} GB"
	    cpu:     "~{threads}"
	    timeout: "~{timeout}"
	}

	output {
    	File snvVcfOutput = "~{basename}.wellMapped.MAF.vcf"
	}
	
	meta {
	    	output_meta: {
      		snvVcfOutput: "filtered SNV VCF output"
    	}
  	}
}


task filterStructural {
	input {
		File vcfFile 
		String basename = basename("~{vcfFile}", ".vcf.gz")
		String modules = "bedtools/2.27.1"
		Int jobMemory = 32
		Int threads = 4
		Int timeout = 16
	}

	parameter_meta {
		vcfFile: "Vcf input file"
		basename: "Base name"
		modules: "Required environment modules"
		jobMemory: "Memory allocated for this job (GB)"
		threads: "Requested CPU threads"
		timeout: "Hours before task timeout"
	}

	command <<<
		set -euo pipefail

		zcat ${vcfFile} | awk '$7 ~ "PASS" {print}'  | \
			awk '  split($8,a,";") split(a[4],b,"=") {print $1"\t"$2"\t"b[2]"\t"$5"\t"b[2]-$2}' | \
			awk '$4 !~ ":" {print}' | awk '$4 !~ "INS" {print}' | \
			sed 's/<//g; s/>//g' >~{basename}.lengths.bed

		bedtools intersect -a ~{basename}.lengths.bed \
			-b GRCh38_notinalldifficultregions.bed -u  >~{basename}.lengths.wellMapped.bed 

	>>> 

	runtime {
	    modules: "~{modules}"
	    memory:  "~{jobMemory} GB"
	    cpu:     "~{threads}"
	    timeout: "~{timeout}"
	}

	output {
    	File structuralVcfOutput = "~{basename}.lengths.wellMapped.bed"
	}

	meta {
    	output_meta: {
      		structuralVcfOutput: "filtered structural VCF output"
    	}
  	}
}

task hrdResults {
	input {
		file snvVcfFiltered
		file structuralVcfFiltered
		String modules = "R"
		Int jobMemory = 32
		Int threads = 4
		Int timeout = 16
	}

	parameter_meta {
		snvVcfFiltered: "filtered SNV Vcf input file"
		structuralVcfFiltered: "filtered structural variant Vcf input file"
		modules: "Required environment modules"
		jobMemory: "Memory allocated for this job (GB)"
		threads: "Requested CPU threads"
		timeout: "Hours before task timeout"
	}

	command <<<
		set -euo pipefail

		Rscript --vanilla CHORD_results.R ~{snvVcfFiltered} ~{structuralVcfFiltered}

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

