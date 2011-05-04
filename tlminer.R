## tlminer.R - mine translocations

library(Rsamtools)
library(RSQLite)


### Find BAM file
TEST.MODE <- system('uname -s', intern=TRUE) == 'Darwin'
# For testing and org-mode usesage.
if (interactive() && TEST.MODE) {
  bamfile <- "CAGTACT.sorted.bam"
} else if(!interactive()) {
  args <- commandArgs()
  arg.delim <- which(args == '--args') + 1
  args <- args[arg.delim:length(args)]
  bamfile <- args[1]
}

### Configure database
dbfile <- sprintf("%s.sqlite", strsplit(basename(bamfile), '.', fixed=TRUE)[1])
message(sprintf("Establishing database connection to %s", dbfile))

# establish database
drv <- dbDriver("SQLite")
con <- dbConnect(drv, dbname=dbfile)

### Split-Mates
# First, we find all BAM entries in which there is a mate that maps to
# the MLL gene. Then, we join these IDs on all others in the database.
mll.region <- RangesList(chr11=IRanges(118307205, 118395934))
splitmate.param <- ScanBamParam(which=mll.region,
                                what=c("qname", "pos", "mrnm", "mpos", "seq", "flag"),
                                flag=scanBamFlag(isUnmappedQuery=FALSE, isMinusStrand=FALSE))
splitmates <- scanBam(bamfile, param=splitmate.param)

