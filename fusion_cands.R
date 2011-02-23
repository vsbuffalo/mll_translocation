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
}

if (!exists('outdir'))
  outdir <- "CAGTACT-output"


# Number of rearrangements to purse
N <- 4

# Basename
basename <- unlist(strsplit(outdir, "-"))[1]

# split-mates directory
split.mates.dir <- "split-mates"

# This is specified in find_split_mates.py as an option; TODO gather from there.
min.mqual <- 30

# Statistics output directory
stats.dir <- paste(outdir, "stats", sep='/')

## Initial analysis: produce plots of rearrangement intensities
ra.file <- paste(outdir, "stats", "reliable-rearrangement-counts.txt", sep='/')
ra.d <- read.csv(ra.file, header=FALSE, sep='\t', col.names=c('rearrangement', 'count'))

# get alternate chromosome name
tmp <- lapply(strsplit(as.character(ra.d$rearrangement), '-'), function(x) x[x != "chr11"])
ra.d$alt.chr <- unlist(tmp)




########## Functions ##########
output = 
# output graphic to file
function(e, file) {
  png(file, width=1200, height=800, res=100)
  eval(e)
  dev.off()
}

getRange =
# from ?cut.
function(labs) {
  lower <- as.numeric(sub("\\((.+),.*", "\\1", labs))
  upper <- as.numeric(sub("[^,]*,([^]]*)\\]", "\\1", labs))
  c(lower=lower, upper=upper)
}


combinePeaks = 
# Check if peaks found are consecutive; this is a bit hacky.
function(x) {
  all <- numeric(length(x)*2)
  j <- 1
  for (i in 1:length(x)) {
    tmp <- getRange(names(x)[i])
    all[j] <- tmp['lower']
    all[j+1] <- tmp['upper']
    j <- j + 2
  }
  all.l <- length(all)
  all.uniq <- unique(all)
  if (length(all.uniq) < all.l) {
    out <- range(all.uniq)
    names(out) <- c('lower', 'upper')
    return(out)
  }
  FALSE
}


findBinModes =
# Find modes of a table of bins by first keeping only counts above q%
# quantile (after removing zeros), and then combining neighboring bins
# with high counts.
function(xtbl, q=0.95) {
  nonzero <- xtbl[xtbl > 0]
  tmp <- nonzero[nonzero > quantile(nonzero, probs=q)]

  # Check for neighboring bins
  cat(sprintf("initial number of non-zero bins: %s\n", length(nonzero)))
  cat(sprintf("number of bins above %s quantile: %s\n", q, length(tmp)))
  combined.tbl <- combinePeaks(tmp)
  if (length(combined.tbl) > 1 && !combined.tbl) {
    message("Could not collapse bins into single peak.")
    browser()
  }  
  return(combined.tbl)
}

findRegions =
# Given a split mates file's data, produce histograms and search for
# translocation candidate regions.
#
# Side effects: produces plots.
function(d, chr.name, alt=NULL) {
  pos.colname <- sprintf("pos.%s", chr.name)
  chr.num <- gsub('chr(\\d+)', '\\1', chr.name)

  # Chr11 will be overwritten; so we pass in alt for special naming
  if (!is.null(alt))
    special <- sprintf("%s-%s", chr.num, alt)
  else
    special <- chr.num  
  
  ## Largest possible range histogram
  p <- ggplot(d, aes_string(x=pos.colname))
  p <- p + geom_histogram(binwidth=diff(range(d[, which(colnames(d) == pos.colname)]))/30)
  p <- p + xlab("position")
  p <- p + ylab(sprintf("count of mates mapped\n(mapping quality > %s)", min.mqual))
  p <- p + opts(title=sprintf("Distribution of reads mapped in chromosome %s", chr.num))
  output({
    print(p)}, file=paste(stats.dir, paste("positions-chr", special, ".png", sep=''), sep='/'))

  ## Do binning
  width.par <- 100
  bins <- cut(d[, pos.colname], breaks=diff(range(d[, pos.colname]))/100)
  peak.range <- findBinModes(table(bins))
  d.peak <- d[d[, pos.colname] >= peak.range['lower'] & d[, pos.colname] <= peak.range['upper'], ]
  
  ## Plot narrowed histogram
  p <- ggplot(d.peak, aes_string(x=pos.colname))
  p <- p + geom_histogram()
  p <- p + xlab("position")
  p <- p + ylab(sprintf("count of mates mapped\n(mapping quality > %s)", min.mqual))
  p <- p + opts(title=sprintf("Distribution of reads mapped in chromosome %s, after selection", chr.num))
  output({
    print(p)}, file=paste(stats.dir, paste("selected-positions-chr", special, ".png", sep=''), sep='/'))
  
  return(list(data=d, data.peak=d.peak, peak.range=peak.range))
}

# temporarily unbin data for histogram - a bit hacky.
local({
  tmp <- data.frame(chromosome=unlist(apply(ra.d, 1, function(x) rep(x[3], x[2]))))
  # order chromosomes
  tmp$chromosome <- factor(tmp$chromosome, levels=paste('chr', c(1:22, 'X', 'Y'), sep=''))
  p <- ggplot(tmp, aes(x=chromosome))
  p <- p + geom_histogram()
  p <- p + opts(title="Rearrangement Candidates\n(counts could include fusion at restriction enzyme site)")
  p <- p + ylab(sprintf("count of paired-end read mates mapped,\n(other mate in chromosome 11, mapping quality > %s)", min.mqual))
  output({
    print(p)}, file=paste(stats.dir, "rearrangement-counts.png", sep='/'))
})


## Investigate each chromosome candidate
split.mate.classes <- c('character', 'character', 'character', 'integer',
                        'integer', 'factor', 'factor', 'integer', 'integer')

## Load Singles File
fusion.read.dir <- paste(outdir, "fusion-reads", sep='/')
if (!file.exists(fusion.read.dir))
  system(sprintf("mkdir %s", fusion.read.dir))

singles.file <- paste(outdir, "split-mates", sprintf("%s-singles.txt", basename), sep='/')

ds <- read.csv(singles.file, header=FALSE, sep='\t',
               colClasses=c('character', 'character', 'integer', 'integer',
                 'factor', 'character', 'integer', 'integer'))
colnames(ds) <- c('name', 'u.seq', 'u.pos', 'u.mqual', 'chr', 'm.seq', 'm.pos', 'm.mqual')

ds <- ds[, -9]

# For each of the top candidate rearrangements (now just top N), get
# and output regions and histograms.
for (i in 4:4) {
  chr.ra.name <- ra.d$rearrangement[i]
  cat(sprintf("Processing rearrangement: %s\n", chr.ra.name))
  split.file.dir <- paste(outdir, split.mates.dir, sprintf("%s-%s.txt", basename, chr.ra.name), sep='/')

  # Load file
  d <- read.csv(split.file.dir, header=FALSE, sep='\t', colClasses=split.mate.classes)
  d <- d[, -10]
  colnames(d) <- c('name', paste(rep(c('seq', 'pos', 'strand', 'mqual'), each=2),
                                 c('chr11', ra.d$alt.chr[i]), sep='.'))

  chr11.regions <- findRegions(d, 'chr11', alt=ra.d$alt.chr[i])
  chralt.regions <- findRegions(d, ra.d$alt.chr[i])
  chralt.range <- chralt.regions$peak.range
  
  tmp <- ds[ds$chr == ra.d$alt.chr[i], ]
  chralt.singles <- tmp[tmp$m.pos >= chralt.range['lower'] & tmp$m.pos <= chralt.range['upper'], ]
  con <- file(sprintf("%s/%s-fusion-reads.fasta", fusion.read.dir, ra.d$alt.chr[i]), open='w')
  for (j in 1:nrow(chralt.singles)) {
    cat(sprintf(">%s\n%s\n", chralt.singles$name[j], chralt.singles$u.seq[j]),
        file=con)
  }
  close(con)
}
