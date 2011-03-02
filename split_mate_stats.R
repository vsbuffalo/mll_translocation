## split_stats.R - Gather statistics about split-mate files; output
## candidate rearrangement subset.


suppressMessages(require(ggplot2))

if (!interactive()) {
  args <- commandArgs()
  arg.delim <- which(args == '--args') + 1
  args <- args[arg.delim:length(args)]
  outdir <- args[1]
}

if (!exists('outdir'))
  outdir <- "CAGTACT-output"

split.mates.dir <- 'split-mates-reliable'

writeFasta =
# Write a fastafile, given headers and sequences.
function(headers, sequences, filename) {
  if (length(headers) != length(sequences))
    stop("Arguments for headers and sequences must be same length.")
  con <- file(filename, open='w')
  for (j in 1:length(headers)) {
    cat(sprintf(">%s\n%s\n", headers[j], sequences[j]), file=con)
  }
  close(con)
}

checkDir =
# if a directory doesn't exist, make it
function(dir) {
  if (!file.exists(dir)) {
    message(sprintf("Making directory '%s'", dir))
    system(sprintf("mkdir -p %s", dir))
  }
  dir
}

clusterByPosition =
# Cluster and output FASTA file sequences by position, only pairs with
# forward-mapped mate on chr11.
function(pos.stats, count.thresh=30, num.groups=10000) {
  x <- table(pos.stats$pos.chralt)

  ## check for lack of peaks
  if (sum(x) <= count.thresh || 0.8*length(x) <= length(which(x >= quantile(x, probs=0.9)))) {
    message(sprintf("No peaks in %s", pos.stats$chralt))
    return()
  }
  
  y <- as.integer(names(x))
  names(y) <- y
  h <- hclust(dist(y))

  groups <- cutree(h, h=num.groups)
  groups.table <- table(groups)

  top.decile <- quantile(groups.table, probs=0.9)
  keep <- as.integer(names(groups.table))[groups.table > top.decile]
  if (length(keep) != 1)
    stop("More than one cluster passes threshold.")
  group.range <- range(as.integer(names(groups)[groups == keep]))
  lower <- group.range[1]
  upper <- group.range[2]

  d <- pos.stats

  headers = d$pos.chralt[d$strand.chr11 == 'forward' & d$pos.chralt > lower & d$pos.chralt < upper]
  if (length(headers) > 0)
    sequences = d$seq.chralt[d$strand.chr11 == 'forward' & d$pos.chralt > lower & d$pos.chralt < upper]
  else
    return()
  
  rootname <- sprintf("%s-grouped-seqs.fasta", pos.stats$chralt)
  dir <- checkDir(file.path(outdir, "assembled-mates"))
  filename <- file.path(dir, rootname)
  writeFasta(headers, sequences, filename)
}


processSplitMateFile =
# Given a split-mate file, gather statistics and output a subset of
# the original file of entries that have mates that map to the forward
# strand of chromosome 11.
function(filename) {
  d <- read.csv(filename, header=FALSE, sep='\t')
  basename <- basename(filename)

  chromosomes <- local({
    tmp <- unlist(strsplit(basename(filename), '[\\.-]'))
    return(tmp[grep('chr', tmp)])
  })

  colnames(d) <- c('name', paste(rep(c('seq', 'pos', 'strand', 'mqual'), each=2),
                                 chromosomes, sep='.'))

  chralt <- chromosomes[chromosomes != 'chr11']
  strand.stats <- list()
  pos.stats <- list()
  pos.stats$chralt <- chralt
  strand.stats$chralt <- chralt

  ## Get strand counts
  chralt.strands <- table(d[, paste('strand', chralt, sep='.')])
  strand.stats$chralt.forward <- chralt.strands['forward']
  strand.stats$chralt.reverse <- chralt.strands['reverse']


  chr11.strands <- table(d$strand.chr11)
  strand.stats$chr11.forward <- chr11.strands['forward']
  strand.stats$chr11.reverse <- chr11.strands['reverse']

  pos.stats$strand.chralt <- d[, paste('strand', chralt, sep='.')]
  pos.stats$strand.chr11 <- d$strand.chr11
  
  pos.stats$name <- d$name
  pos.stats$pos.chr11 <- d$pos.chr11
  pos.stats$pos.chralt <- d[, paste('pos', chralt, sep='.')]
  pos.stats$mqual.chr11 <- d$mqual.chr11
  pos.stats$mqual.chralt <- d[, paste('mqual', chralt, sep='.')]
  pos.stats$seq.chralt <- d[, paste('seq', chralt, sep='.')]

  candidates <- subset(d, strand.chr11 == 'forward')
  strand.stats$num.candidates <- nrow(candidates)
  
  ## write.table(candidates,
  ##             file=file.path(outputdir, sprintf("%s-subset.txt", basename)),
  ##             quote=FALSE, row.names=FALSE, sep='\t')

  ## Cluster by position, output grouped sequences.
  clusterByPosition(pos.stats)
  
  return(list(strand.stats=strand.stats, pos.stats=pos.stats))
}

chr11.split.mates.files <- dir(file.path(outdir, split.mates.dir), pattern="chr11")
output = lapply(chr11.split.mates.files, function(fn) processSplitMateFile(file.path(outdir, split.mates.dir, fn)))

out.name <- file.path(outdir, "stats",
                      sprintf("%s-candidate-summary.txt", unlist(strsplit(outdir, '-'))[1]))

# Output lists contain statistics for position and strand statistics -
# unaggregate these.
strand.stats <- lapply(output, function(x) x$strand.stats)
pos.stats <- lapply(output, function(x) x$pos.stats)

write.table(do.call(rbind, strand.stats), file=out.name, quote=FALSE, row.names=FALSE, sep='\t')

