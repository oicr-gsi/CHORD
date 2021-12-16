version 1.0

workflow CHORD {
    meta {
        author: "Vivek Alamuri"
        email: "valamuri@oicr.on.ca"
        description: "Workflow that runs CHORD to predict HRD"
        dependencies: [
            {
                name: "rstats/4.0",
                url: "https://www.r-project.org/"
            }
        ]
    }

    input {
        File snvFile
        File indelFile
        File svFile
        File script
    }

    call runCHORD {
        input:
            snvFile = snvFile,
            indelFile = indelFile,
            svFile = svFile,
            script = script
    }
}

task runCHORD {
    input {
        File snvFile
        File indelFile
        File svFile
        File script
        String modules = "chord/2.0 hg38-refgene/p12 hg19-refgene/p13"
        Int threads = 8
        Int jobMemory = 64
        Int timeout = 72
    }

    parameter_meta {
        snvFile: "SNV VCF filepath"
        indelFile: "indel VCF filepath"
        svFile: "SV VCF filepath"
        script: "Rscript filepath"
        modules: "Names and versions of modules to load"
        threads: "Requested CPU threads"
        jobMemory: "Memory allocated for this job"
        timeout: "Hours before task timeout"
    }

    command <<<
        Rscript ~{script} -n ~{snvFile} -i ~{indelFile} -v ~{svFile}
    >>>

    runtime {
        memory:  "~{jobMemory} GB"
        modules: "~{modules}"
        cpu:     "~{threads}"
        timeout: "~{timeout}"
    }

    meta {
        output_meta {
            contexts = "contexts text file"
            predictions = "predictions text file"
        }
    }

    output {
        File contexts = "merged_contexts.txt"
        File predictions = "chord_pred.txt"
    }
}
