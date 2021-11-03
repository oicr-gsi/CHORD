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
    }

    call runCHORD {
        input:
            snvFile = snvFile,
            indelFile = indelFile,
            svFile = svFile
    }
}

task runCHORD {
    input {
        File snvFile
        File indelFile
        File svFile
    }

    parameter_meta {
        snvFile: "SNV VCF filepath"
        indelFile: "indel VCF filepath"
        svFile: "SV VCF filepath"
    }

    command {
        Rscript runCHORDScript.R -n ~{snvFile} -i ~{indelFile} -v ~{svFile}
    }
}

# to validate
# java -jar womtool-69.jar validate wdls/CHORD-v1_wkflow_inputs.wdl
