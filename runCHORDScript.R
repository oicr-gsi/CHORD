library(optparse)

option_list = list(
  make_option(c("-w", "--wkdir"), type="character", default=NULL,
              help="**required** working directory where vcfs folder is located", metavar="character"),
  
  make_option(c("-n", "--snvPattern"), type="character", default="*snp.vcf",
              help="SNV files pattern; default: *snp.vcf", metavar="character"),
  
  make_option(c("-v", "--svPattern"), type="character", default="*delly.merged.vcf.gz",
              help="SV files pattern; default: *delly.merged.vcf.gz", metavar="character"),
  
  make_option(c("-i", "--indelPattern"), type="character", default="*indel*",
              help="indel files pattern; default: *indel*", metavar="character")
  # make_option(c("-i", "--snv"), type="character", default="chord_pred.txt",
  #             help="predictions output file name; default: chord_pred.txt", metavar="character")
  # make_option(c("-v", "--sv"), type="character", default="chord_pred.txt",
  #             help="predictions output file name; default: chord_pred.txt", metavar="character")
);

# get command line arguments
opt_parser = OptionParser(option_list=option_list);
args = parse_args(opt_parser)

workingDir <- args$wkdir
snvPattern <- args$snvPattern
svPattern <- args$svPattern
indelPattern <- args$indelPattern

mergedContextsFileName <- "merged_contexts.txt"
predictionOutputFileName <- "chord_predictions.txt"

# libraries imported after command line arguments so that when --help is called, these are not loaded
library(CHORD)
library(BSgenome.Hsapiens.UCSC.hg19)
options(stringsAsFactors=F)

setwd(workingDir)


vcf_files <- data.frame(
  snv=list.files('vcf', pattern=snvPattern, full.names=T),
  sv=list.files('vcf', pattern=svPattern, full.names=T),
  indel=list.files('vcf', pattern=indelPattern, full.names = T)
)

vcf_files$sample <- sapply(strsplit(basename(vcf_files$snv),'_'),`[`,1)

dir.create("output/", recursive = T)
contexts_dir <- 'output/contexts/'
dir.create(contexts_dir, recursive=T)

## Extract signatures for all samples
for(i in 1:nrow(vcf_files)){
  params <- as.list(vcf_files[i,])
  out_path <- paste0(contexts_dir,'/',params$sample,'_contexts.txt')
  
  if(!file.exists(out_path)){
    extractSigsChord(
      vcf.snv=params$snv,
      vcf.indel = params$indel,
      vcf.sv=params$sv, 
      sv.caller='manta',
      sample.name=params$sample,
      output.path=out_path, verbose=F, ref.genome = BSgenome.Hsapiens.UCSC.hg19,
      
      vcf.filters=list(snv="PASS",indel="PASS",sv="PASS") 
    )
  }
}

context_files <- list.files(contexts_dir, full.names=T)

l_contexts <- lapply(context_files, function(i){
  read.delim(i, check.names=F)
})

## Merge the contexts into a matrix
merged_contexts <- do.call(rbind, l_contexts)

## Write to output directory
write.table(merged_contexts, paste0("output", "/", mergedContextsFileName), sep='\t', quote=F)

chord_output <- chordPredict(merged_contexts, do.bootstrap=T, verbose=F)
write.table(chord_output, paste0("output", "/", predictionOutputFileName), sep='\t', quote=F)

if (file.exists(paste0("output", "/", predictionOutputFileName))) {
  paste("Prediction outputs written to", predictionOutputFileName)
}
