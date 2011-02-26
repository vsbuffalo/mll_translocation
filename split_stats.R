## split_stats.R - Gather statistics about split-mate files; output
## candidate rearrangement subset.

if (!interactive()) {
  args <- commandArgs()
  arg.delim <- which(args == '--args') + 1
  args <- args[arg.delim:length(args)]
  outdir <- args[1]
}

if (!exists('outdir'))
  outdir <- "CAGTACT-output"

split.mates.dir <- 'split-mates-reliable'

pathJoin =
# Make a path out of directory and file parts.
function(...) {
  args <- list(...)
  args$sep = '/'
  do.call(paste, args)
}

getBasename =
# Given a path to a file like test/dir/filename.ext, extract
# "filename".
function(filename) {
  parts <- unlist(strsplit(filename, '/'))
  filename <- parts[length(parts)]
  basename <- unlist(strsplit(filename, '.', fixed=TRUE))[1]
  return(basename)
}

writeFasta =
# Write a fastafile, given headers and sequences.
function(headers, seqeunces, filename) {
  if (length(header) != length(seqeunces))
    error("Arguments for headers and sequences must be same length.")

  con <- file(filename, open='w')
  for (j in 1:length(headers)) {
    cat(sprintf(">%s\n%s\n", headers[j], sequences[j]), file=con)
  }
  close(con)
}

processSplitMateFile =
# Given a split-mate file, gather statistics and output a subset of
# the original file of entries that have mates that map to the forward
# strand of chromosome 11.
function(filename, outputdir=pathJoin(outdir,'fusion-reads')) {
  if (!file.exists(outputdir)) {
    message(sprintf("Creating output directory '%s'", outputdir))
    system(sprintf("mkdir %s", outputdir))
  }
    
  d <- read.csv(filename, header=FALSE, sep='\t')
  basename <- getBasename(filename)

  chromosomes <- local({
    tmp <- unlist(strsplit(filename, '[\\.-]'))
    return(tmp[grep('chr', tmp)])
  })

  colnames(d) <- c('name', paste(rep(c('seq', 'pos', 'strand', 'mqual'), each=2),
                                 chromosomes, sep='.'))

  chralt <- chromosomes[chromosomes != 'chr11']
  stats <- 
  stats$chralt <- chralt

  ## Get strand counts
  chralt.strands <- table(d[, paste('strand', chralt, sep='.')])
  stats$strands.chralt.forward <- chralt.strands['forward']
  stats$strands.chralt.reverse <- chralt.strands['reverse']


  chr11.strands <- table(d$strand.chr11)
  stats$strands.chr11.forward <- chr11.strands['forward']
  stats$strands.chr11.reverse <- chr11.strands['reverse']
  ## stats$pos.chr11 <- d$pos.chr11
  ## stats$pos.chralt <- d[, paste('pos', chralt, sep='.')]

  candidates <- subset(d, strand.chr11 == 'forward')
  stats$num.candidates <- nrow(candidates)
  
  ## write.table(candidates,
  ##             file=pathJoin(outputdir, sprintf("%s-subset.txt", basename)),
  ##             quote=FALSE, row.names=FALSE, sep='\t')

  return(stats)
}

chr11.split.mates.files <- dir(pathJoin(outdir, split.mates.dir), pattern="chr11")
a = lapply(chr11.split.mates.files, function(fn) processSplitMateFile(pathJoin(outdir, split.mates.dir, fn)))

do.call(rbind, a)
