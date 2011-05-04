## tlminer.R - mine translocations
suppressMessages({
  library(Rsamtools)
  library(RSQLite)
  library(org.Hs.eg.db)
})

source("plotting-methods.R")

### Config
mll.region <- RangesList(chr11=IRanges(118307205, 118395934))

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
## dbfile <- sprintf("%s.sqlite", strsplit(basename(bamfile), '.', fixed=TRUE)[1])
## message(sprintf("Establishing database connection to %s", dbfile))

## # establish database
## drv <- dbDriver("SQLite")
## con <- dbConnect(drv, dbname=dbfile)

### Alignment Statistics
unmapped <- countBam(bamfile, param=ScanBamParam(flag=scanBamFlag(isUnmappedQuery=TRUE)))
mapped <- countBam(bamfile, param=ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE)))
mapped.mll.forward <- countBam(bamfile, param=ScanBamParam(which=mll.region,
                                          flag=scanBamFlag(isMinusStrand=FALSE,
                                            isUnmappedQuery=FALSE)))
mapped.mll.reverse <- countBam(bamfile, param=ScanBamParam(which=mll.region,
                                          flag=scanBamFlag(isMinusStrand=TRUE,
                                            isUnmappedQuery=FALSE)))
### MLL Region coverage plots



### Split-Mates
# First, we find all BAM entries in which there is a mate that maps to
# the MLL gene. Then, we join these IDs on all others in the database.
mll.forward.param <- ScanBamParam(which=mll.region,
                                what=c("rname", "pos", "mrnm", "mpos", "seq", "flag"),
                                flag=scanBamFlag(isUnmappedQuery=FALSE, isMinusStrand=FALSE))
mll.forward <- scanBam(bamfile, param=mll.forward.param)

# Find all entries with mrnm that is chr11 and remove these.
keep <- with(mll.forward[[1]], mrnm != 'chr11')
splitmates <- lapply(mll.forward[[1]], function(x) x[keep])


d <- with(splitmates, data.frame(rname, pos, mrnm, mpos))
aggregate(d$mpos, list(d$mrnm, d$mpos), length)

### Summarize the split-mates



#keep <- with(splitmates[[1]], cbind(rname, pos, mrnm, mpos, seq, flag))



### Mapped-Unmapped pairs
mappedunmapped.param <- ScanBamParam(what=c("qname", "rname", "pos", "mrnm", "mpos", "seq", "flag"),
                                     flag=scanBamFlag(isUnmappedQuery=TRUE, hasUnmappedMate=FALSE))
mappedunmapped <- scanBam(bamfile, param=mappedunmapped.param)
