# split_stats.R - gather statistics on spit files.
suppressMessages(require(ggplot2))

if (!interactive()) {
  args <- commandArgs()
  arg.delim <- which(args == '--args') + 1
  args <- args[arg.delim:length(args)]
  file <- args[1]
  singles.file <- args[2]
  outdir <- args[3]
} else {
  file <- "CAGTACT-output/split-mates/CAGTACT-chr11-chr9.txt"
  singles.file <- "CAGTACT-output/split-mates/CAGTACT-singles.txt"
  outdir <- "CAGTACT-output/stats/"
}


# grab which chromosomes are are in this file
parts <- unlist(strsplit(file, '[-\\.]', perl=TRUE))
chrs <- parts[grep('chr', parts)]

# load file
d <- read.csv(file, header=FALSE, sep='\t',
              colClasses=c('character', 'character', 'character',
                'integer', 'integer', 'factor', 'factor', 'integer', 'integer'))
colnames(d) <- c('name', paste(rep(c('seq', 'pos', 'strand', 'mqual'), each=2), chrs, sep='.'))

# get positions of position columns
pos.chr11.col <- grep('pos.chr11', colnames(d))
tmp <- grep('pos.chr', colnames(d))
pos.altchr.col <- tmp[which(tmp != pos.chr11.col)]
chralt.name <- chrs[chrs != 'chr11']
chralt.pos <- paste("pos", chralt.name, sep='.')

# find chr11 peak
chr11.t <- table(cut(d[, pos.chr11.col], 10000))
chr11.peaks <- chr11.t[chr11.t > mean(unique(chr11.t))]

# find chralt peak
chralt.t <- table(cut(d[, pos.altchr.col], 10000))
chralt.peaks <- chralt.t[chralt.t > mean(unique(chralt.t))]

#cat("number of chr11 peaks about mean unique bin counts", length(chr11.peaks), "\n", sep='\t')

if (length(chr11.peaks) > 1)
  stop("More than one chr11 peak found")

if (length(chralt.peaks) > 1)
  stop("More than one chralt peak found")

#cat("peak bin", names(chr11.peaks), "\n", sep='\t')

getRange <- function(labs) {
  lower <- as.numeric(sub("\\((.+),.*", "\\1", labs))
  upper <- as.numeric(sub("[^,]*,([^]]*)\\]", "\\1", labs))
  c(lower=lower, upper=upper)
}

## chr11.cords <- getRange(names(chr11.peaks))
## candidate.pairs <- d[chr11.cords['lower'] >= d[, pos.chr11.col] &
##                      d[, pos.chr11.col] <= chr11.cords['upper'], ]

chralt.cords <- getRange(names(chralt.peaks))
candidate.pairs <-  d[chralt.cords['lower'] <= d[, pos.altchr.col] &
                      d[, pos.altchr.col] <= chralt.cords['upper'], ]

chralt.range <- range(candidate.pairs[, chralt.pos])
cat(sprintf("selected %s lower", chralt.name), chralt.range[1], "\n", sep='\t')
cat(sprintf("selected %s upper", chralt.name), chralt.range[2], "\n", sep='\t')

chr11.range <- range(candidate.pairs$pos.chr11)
cat("selected chr11 lower", chr11.range[1], "\n", sep='\t')
cat("selected chr11 upper", chr11.range[2], "\n", sep='\t')

# chr11 hist
png(paste(outdir, "chr11-histgram.png", sep='/'),
    width=900, height=900, res=100)
p <- ggplot(data=candidate.pairs)
p <- p + geom_histogram(aes(x=pos.chr11),
                        binwidth=diff(range(candidate.pairs$pos.chr11))/60)
p + xlab("position") + opts(title='chr11')
garbage <- dev.off()

# chralt hist
png(paste(outdir, sprintf("%s-histgram.png", chralt.name), sep='/'),
    width=900, height=900, res=100)
p <- ggplot(data=candidate.pairs)
p <- p + geom_histogram(aes(x=candidate.pairs[, chralt.pos]),
                        binwidth=diff(range(candidate.pairs[, chralt.pos]))/60)
p + xlab("position") + opts(title=chralt.name)
garbage <- dev.off()


## Load singles file
ds <- read.csv(singles.file, header=FALSE, sep='\t',
              colClasses=c('character', 'character', 'integer', 'integer',
                'factor', 'character', 'integer', 'integer'))
colnames(ds) <- c('name', 'u.seq', 'u.pos', 'u.mqual', 'chr', 'm.seq', 'm.pos', 'm.mqual')


## Get relevant sequences from chralt, and within range
tmp <- ds[ds$chr == chralt.name, ]

# in range
chralt.singles <- tmp[tmp$m.pos >= chralt.range[1] &  tmp$m.pos <= chralt.range[2],]

# only keep *perfect* mapping quality of 37
chralt.singles <- chralt.singles[chralt.singles$m.mqual == 37, ]

for (i in 1:nrow(chralt.singles)) {
  cat(sprintf(">%s\n%s\n", chralt.singles$name[i], chralt.singles$u.seq[i]),
      file="out.fasta", append=TRUE)
}
