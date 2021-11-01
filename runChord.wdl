version 1.0

workflow CHORD {
    # TODO: add Rscript dependencies under meta 
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
    call runCHORD

    # TODO: does not allow to pass inputs from worflow level 
    # input {
    #     String wkdir
    #     String snvPattern
    #     String indelPattern
    #     String svPattern
    # }

    # call runCHORD {
    #     input:
    #         wkdir = wkdir
    #         snvPattern = snvPattern
    #         indelPattern = indelPattern
    #         svPattern = svPattern
    # }
}

task runCHORD {
    input {
        String wkdir
        String snvPattern
        String indelPattern
        String svPattern
    }

    parameter_meta {
        wkdir: "required; working directory where vcfs folder is located"
        snvPattern: "SNV VCF files pattern; default: *snp.vcf"
        indelPattern: "indel VCF files pattern; default: *indel*"
        svPattern: "SV VCF files pattern; default: *delly.merged.vcf.gz"
    }

    # echo variables for now, will be passed to R script later
    command {
        echo 'path: -w ~{wkdir}, snvPattern: -n ~{snvPattern}, indelPattern: -i ~{indelPattern},
        svPattern: -v ~{svPattern}'
    }

}

