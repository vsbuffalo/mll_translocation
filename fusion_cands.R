## fusion_cands.R - find fusion candidates from data.
##
## This begins by looking at the
## BARCODE-output/stats/reliable-rearrangement-counts.txt file. It
## produces:
##
## - "intensities" for rearrangements
##
## - candidate translocated regions from the split-mates data
##
## - FASTA files of the unmapped reads for each rearrangement
##   candidate.

suppressMessages(require(ggplot2))

if (!interactive()) {
  args <- commandArgs()
  arg.delim <- which(args == '--args') + 1
  args <- args[arg.delim:length(args)]
  outdir <- args[1]
} else {
  outdir <- "CAGTACT-output"
}

# This is specified in find_split_mates.py as an option; TODO gather from there.
min.mqual <- 30

## Initial analysis: produce plots of rearrangement intensities
ra.file <- paste(outdir, "stats", "reliable-rearrangement-counts.txt", sep='/')
ra.d <- read.csv(ra.file, header=FALSE, sep='\t', col.names=c('rearrangement', 'count'))

# get alternate chromosome name
tmp <- lapply(strsplit(as.character(ra.d$rearrangement), '-'), function(x) x[x != "chr11"])
ra.d$alt.chr <- unlist(tmp)

# temporarily unbin data for histogram - a bit hacky.
local({
  tmp <- data.frame(chromosome=unlist(apply(ra.d, 1, function(x) rep(x[3], x[2]))))
  # order chromosomes
  tmp$chromosome <- factor(tmp$chromosome, levels=paste('chr', c(1:22, 'X', 'Y'), sep=''))
  p <- ggplot(tmp, aes(x=chromosome))
  p <- p + geom_histogram()
  p <- p + opts(title="Rearrangement Candidates\n(counts could include fusion at restriction enzyme site)")
  p <- p + ylab(sprintf("count of paired-end read mates mapped,\n(other mate in chromosome 11, mapping quality > %s)", min.mqual))
  print(p)
})
