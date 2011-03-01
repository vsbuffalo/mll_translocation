## split_read_stats.R - statistics (and FASTA file output) of
## paired-end reads with one mate mapped and the other mate unmapped
## (and thus a candidate for containing a hybrid read).

suppressMessages(require(ggplot2))

if (!interactive()) {
  args <- commandArgs()
  arg.delim <- which(args == '--args') + 1
  args <- args[arg.delim:length(args)]
  outdir <- args[1]
}

if (!exists('outdir'))
  outdir <- "CAGTACT-output"

writeFasta =
# Write a fastafile, given headers and sequences.
function(headers, sequences, filename) {
  if (length(headers) != length(sequences))
    error("Arguments for headers and sequences must be same length.")

  con <- file(filename, open='w')
  for (j in 1:length(headers)) {
    cat(sprintf(">%s\n%s\n", headers[j], sequences[j]), file=con)
  }
  close(con)
}

processSinglesFile =
# Given a file in which one mate mapped to chr11, and another is
# unmapped, process, gather statistics, and write FASTA file per a
# mapped chromosome's unmapped mate.
function(filename, mqual.thresh=30, outputdir=file.path(outdir, 'fusion-reads')) {

  if (!file.exists(outputdir)) {
    message(sprintf("Creating output directory '%s'", outputdir))
    system(sprintf("mkdir %s", outputdir))
  }

  d <- read.csv(filename, header=FALSE, sep='\t',
                colClasses=c('character', 'character', 'integer', 'integer',
                  'factor', 'character', 'integer', 'integer'),
                col.names=c('name', 'useq', 'upos', 'umqual', 'mchr', 'mseq', 'mpos', 'mmqual'))

  by.chr <- split(d, d$mchr)
  for (chr in names(by.chr)) {
    tmp <- by.chr[[chr]]
    tmp <- tmp[tmp$mmqual > mqual.thresh, ]
    if (nrow(tmp) > 0) {
      writeFasta(tmp$name, tmp$useq,
                 file.path(outputdir, sprintf("%s-fusion-candidates.fasta", chr)))
    } else {
      message(sprintf("No output for %s", chr))
    }
  }

}

processSinglesFile(file.path(outdir, "split-mates", "singles.txt"))
