# split_stats.R - gather statistics on spit files.
require(ggplot2)

args = commandArgs()

# grab which chromosomes are are in this file
file <- args[[1]]
parts <- unlist(strsplit(file, '[-\\.]', perl=TRUE))
chrs <- parts[grep('chr', parts)]

# load file
d <- read.csv(file, header=FALSE, sep='\t',
              colClasses=c('character', 'character', 'character',
                'integer', 'integer', 'factor', 'factor'))
colnames(d) <- c('name', paste(rep(c('seq', 'pos', 'strand'), each=2), chrs, sep='.'))

# get positions of position columns
pos.chr11.col <- grep('pos.chr11', colnames(d))
tmp <- grep('pos.chr', colnames(d))
pos.altchr.col <- tmp[which(tmp != pos.chr11.col)]

# bin by alt chromosome's position
cat("altchr range", range(d[, pos.altchr.col]), "\n", sep='\t')

# find chr11 peak
chr11.t <- table(cut(d[, pos.chr11.col], 10000))
chr11.peaks <- chr11.t[chr11.t > mean(unique(t))]

# find chralt peak
chralt.t <- table(cut(d[, pos.altchr.col], 10000))
chralt.peaks <- chralt.t[chralt.t > mean(unique(t))]

cat("number of chr11 peaks about mean unique bin counts", length(chr11.peaks), "\n", sep='\t')

if (length(chr11.peaks))
  stop("More than one chr11 peak found")

if (length(chralt.peaks))
  stop("More than one chralt peak found")

cat("peak bin", names(chr11.peaks), "\n", sep='\t')

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

tmp <- range(candidate.pairs$pos.chr9)
cat("selected chralt lower", tmp[1], "\n", sep='\t')
cat("selected chralt upper", tmp[2], "\n", sep='\t')

tmp <- range(candidate.pairs$pos.chr11)
cat("selected chr11 lower", tmp[1], "\n", sep='\t')
cat("selected chr11 upper", tmp[2], "\n", sep='\t')

# chr11 hist
p <- ggplot(data=candidate.pairs)
p <- p + geom_histogram(aes(x=pos.chr11))
p + xlab("position") + opts(title='chr11')

# chralt hist
chralt.name <- colnames(d)[pos.altchr.col]
p <- ggplot(data=candidate.pairs)
p <- p + geom_histogram(aes(x=candidate.pairs[,chralt.name]))
p + xlab("position") + opts(title=chrs[chrs != 'chr11'])

