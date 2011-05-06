## tlminer.R - mine translocations
suppressMessages({
  library(Rsamtools)
  library(RSQLite)
  library(org.Hs.eg.db)
})

source("plotting-methods.R")
source("utils.R")

### Config
mll.region <- RangesList(chr11=IRanges(118307205, 118395934))
ref <- "mll_template/mll.fasta"
checkBWA(dirname(ref))

### Find BAM file
TEST.MODE <- system('uname -s', intern=TRUE) == 'Darwin'
# For testing and org-mode usesage.
if (interactive() || TEST.MODE) {
  bamfile <- "CAGTACT.sorted.bam"
} else if(!interactive()) {
  args <- commandArgs()
  arg.delim <- which(args == '--args') + 1
  args <- args[arg.delim:length(args)]
  bamfile <- args[1]
}

rootname <- getRootname(bamfile)

### Set up directories and check references
dirs <- makeResultsDir(rootname)
bwacmd <- inPath(c("bwa", "/usr/local/bin/bwa"))
cdhitcmd <- inPath(c("cd-hit", "/usr/local/bin/cd-hit"))
samtoolscmd <- inPath('samtools')

if (!TEST.MODE) {
  message("Checking human genome reference is properly indexed for BWA.")
  hg.refdir <- "/classico/jfass/projects/Vaughan/ref_hg19"
  hg.ref <- file.path(hg.refdir, "hg19_24chrom.fasta")
  checkBWA(hg.refdir)
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
                                what=c("rname", "pos", "qwidth", "mrnm", "mpos", "seq", "flag"),
                                flag=scanBamFlag(isUnmappedQuery=FALSE, isMinusStrand=FALSE))
mll.forward <- scanBam(bamfile, param=mll.forward.param)

# Find all entries with mrnm that is chr11 and remove these.
keep <- with(mll.forward[[1]], mrnm != 'chr11')
splitmates <- lapply(mll.forward[[1]], function(x) x[keep])


### Summarize the split-mates
d <- with(splitmates, data.frame(rname, pos, mrnm, mpos, qwidth))
aggregate(d$mpos, list(d$mrnm, d$mpos), length)

### Look for islands
altchr.mappings = local({
  tmp <- split(d, list(d$mrnm))
  lapply(tmp, function(x) {
    IRanges(start=x$mpos, width=x$qwidth)    
  })
})

findIslands <- function(x, minCoverage=60) {
  lapply(altchr.mappings, function(x) {
    slice(coverage(x), minCoverage)
  })
}

# Find islands of coverage, and their max coverage
splitmate.islands <- findIslands(altchr.mappings)
splitmate.max.cov <- lapply(splitmate.islands, viewMaxs)


### Mapped-Unmapped pairs
mappedunmapped.param <- ScanBamParam(what=c("qname", "rname", "pos", "mrnm", "mpos", "seq", "flag"),
                                     flag=scanBamFlag(isUnmappedQuery=TRUE, hasUnmappedMate=FALSE))
mappedunmapped <- scanBam(bamfile, param=mappedunmapped.param)

## Group by chromosome and write to file
with(mappedunmapped[[1]], {
  seqs <- split(seq, mrnm)
  headers <- split(qname, mrnm)
  mpos <- split(mpos, mrnm)
  mstrand <- sapply(split(flag, mrnm), function(x) c("+", "-")[(bitAnd(x, 32L) > 0L) + 1 ])
  lapply(names(seqs), function(chr) {
    tmp <- seqs[[chr]]
    names(tmp) <- sprintf("%s;;mpos=%s;;mstrand=%s", headers[[chr]], mpos[[chr]], mstrand[[chr]])
    fn <- file.path(dirs$unmapped.mates, sprintf("%s.fasta", chr))
    write.XStringSet(tmp, file=fn, format="fasta")
    TRUE
  })
})

## Map mapped-unmapped pairs
# run long read aligner on all samples
message("Running BWA bwasw on unmapped mates.")
bwarun <- "%s bwasw -T 10 -c 5 -t 3 %s %s > %s 2> /dev/null"
for (fasta.file in dir(dirs$unmapped.mates, pattern="\\.fasta")) {
  chr <- getRootname(fasta.file)
  samfile <- file.path(dirs$aln, sprintf("%s.sam", chr))
  message(sprintf(" - mapping: %s", fasta.file))
  system(sprintf(bwarun, bwacmd, ref, file.path(dirs$unmapped.mates, fasta.file), samfile))

  message(sprintf(" - converting to BAM: %s", samfile))
  bamfile <- file.path(dirs$aln, sprintf("%s.bam", chr))
  system(sprintf("%s view -S -b -o %s %s 2> /dev/null", samtoolscmd, bamfile, samfile))
}

## Read in alignments for each chromosome and search for fusion sites
message("Writing tail sequences to file.")
for (bamfile in dir(dirs$aln, pattern="\\.bam")) {
  chr <- getRootname(bamfile)
  param <- ScanBamParam(simpleCigar=FALSE, flag=scanBamFlag(isUnmappedQuery=FALSE),
                        what=c("qname", "rname", "pos", "strand", "seq", "cigar"))
  aln <- scanBam(file.path(dirs$aln, bamfile), param=param)[[1]]
  fd <- data.frame(as.character(aln$seq), aln$cigar, aln$pos, stringsAsFactors=FALSE)
  tmp <- apply(fd, 1, function(x) {
    f <- findFusion(x[1], x[2], as.numeric(x[3]))
    unlist(f)
  })
  no.fusion <- sapply(tmp, is.null)
  fusions <- tmp[!no.fusion]
  fusions <- data.frame(qname=aln$qname[!no.fusion], aln$rname[!no.fusion], 
                        do.call(rbind, fusions), stringsAsFactors=FALSE)
  tailseqs <- as(fusions$unmapped, "XStringSet")
  names(tailseqs) <- sprintf("%s;;split=%s", fusions$qname, fusions$break.pos)
  fn <- file.path(dirs$tailseqs, sprintf("%s.fasta", chr))
  write.XStringSet(tailseqs, file=fn, format="fasta")
}

message("Running cd-hit on fusion candidate reads.")
for (fasta.file in dir(dirs$tailseqs, pattern="\\.fasta$")) {
  chr <- getRootname(fasta.file)
  fn <- file.path(dirs$tailseqs, fasta.file)
  cfn <- file.path(dirs$cluster, sprintf("%s-clusters.fasta", chr))
  ok <- system(sprintf("%s -i %s -o %s -g 1 -d 0 > /dev/null", cdhitcmd, fn, cfn))
  stopifnot(ok == 0)
}

## Process clustering results
all.clusters <- list()
message("Processing cd-hit clustering results.")
# Grab representative sequences and corresponding clstr files and process.
for (fasta.file in dir(dirs$cluster, pattern="\\-clusters.fasta$")) {
  chr <- strsplit(getRootname(fasta.file), '-', fixed=TRUE)[[1]][1]
  fn <- file.path(dirs$cluster, fasta.file)
  clstr.file <- file.path(dirs$cluster, sprintf("%s.clstr", fasta.file))
  
  d <- processClstrFile(clstr.file)
    if (length(d) == 0)
      next()
  
  # process .clstr files
  rep.seqs <- as.data.frame(cbind(d))
  
  # read FASTA file
  clusters <- local({
    tmp <- readFASTA(fn, strip.descs=TRUE)
    headers <- lapply(tmp, function(x) x[[1]])
    seqs <- lapply(tmp, function(x) x[[2]])
    
    tmp <- as.data.frame(cbind(seqs))
    rownames(tmp) <- headers
    tmp
  })
  
  # match counts and sequence
  clusters <- merge(clusters, rep.seqs, by.x=0, by.y=0)
  clusters <- cbind(chr, clusters)
  colnames(clusters) <- c('chr', 'name', 'seq', 'count')

  # extract out mate mapping position and break point from clusters
  clusters$mate.pos <- as.numeric(sub(".*;;mpos=(\\d+).*", "\\1", clusters$name))
  clusters$split <- as.numeric(sub(".*;;split=(\\d+).*", "\\1", clusters$name))
  clusters$mate.strand <- sub(".*;;mstrand=([\\+\\-]).*", "\\1", clusters$name)
  
  all.clusters[[chr]] <- clusters
}
  
# bit of name mangling here...
clusters <- do.call(rbind, all.clusters)
rownames(clusters) <- NULL
cluster.cands <- local({
  d <- clusters[clusters$count > 10, c('chr', 'seq', 'count', 'split', 'mate.pos', 'mate.strand')]
  fn <- file.path(dirs$results, "clustered-tailseqs.txt")
  write.table(as.matrix(d), fn, row.names=FALSE, quote=FALSE, sep='\t')
  d
})

### Align clustered tail sequences to human genome
if (!TEST.MODE) {
  fn <- file.path(dirs$cluster.aln, "cluster-seqs.fasta")
  local({
    tmp <- as(unlist(cluster.cands$seq), "XStringSet")
    names(tmp) <- with(cluster.cands, sprintf("chr=%s;;split=%s;;mpos=%s;;mstrand=%s", chr, split, mate.pos, mate.strand))
    write.XStringSet(tmp, file=fn, format="fasta")
  })
  edit.dist <- 8
  
  message("Running BWA aln on cluster sequences (to human genome).")
  bwarun <- "%s aln -n %s %s %s > %s 2> /dev/null"
  aln.file <- file.path(dirs$cluster.aln, "cluster-seqs.sai")
  system(sprintf(bwarun, bwacmd, edit.dist, hg.ref, fn, aln.file))
    
  #bwa samse database.fasta aln_sa.sai short_read.fastq > aln.sam
  cluster.sam <- file.path(dirs$cluster.aln, "cluster-seqs.sam")
  system(sprintf("%s samse %s %s %s > %s 2> /dev/null",
                 bwacmd, hg.ref, aln.file, fn, cluster.sam))

  # convert to bam
  bamfile <- file.path(dirs$cluster.aln, "cluster-seqs.bam")
  system(sprintf("%s view -S -b -o %s %s 2> /dev/null", samtoolscmd, bamfile, cluster.sam))
}


### Query out overlaps between aligned tail sequences and split-mate regions
tailseqs.aln.param <- ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE), 
                                   what=c("rname", "qname", "pos", "cigar", "qwidth", "mapq", "strand", "seq"))
if (!TEST.MODE) {
  aln <- scanBam(bamfile, param=tailseqs.aln.param)[[1]]
  
  tailseqs <- local({

    # Fit the BAM lists into a nice dataframe, extracting the split
    # from the qname header info. Note that this contains the columns
    # "chr", which is the chromosome the tail sequence mapped to, and
    # "mate.chr", which is the chromosome the mapping mate mapped to.
    d <- with(aln, {
      mate.chr <- sub(".*chr=(chr[0-9XYM]+).*", "\\1", qname)
      mate.strand <- sub(".*;;mstrand=([\\+\\-]).*", "\\1", qname)
      split <- as.numeric(sub(".*split=([0-9]+).*", "\\1", qname))
      pos <- as.numeric(sub(".*pos=([0-9]+).*", "\\1", qname))
      data.frame(chr=as.character(rname), start=pos, width=qwidth, strand=strand, seq=as.character(seq), 
                 split=split, mate.chr=mate.chr, mate.pos=mate.pos, mapq=mapq, cigar=cigar, stringsAsFactors=FALSE)
    })

    tmp <- apply(d, 1, function(x)
                 GRanges(x[1], IRanges(as.numeric(x[2]), width=as.numeric(x[3]))))
    gr <- GRangesList(tmp)
  })

  no.chr11 <- names(altchr.mappings)[names(altchr.mappings) != 'chr11']
  splitmates.gr <- lapply(no.chr11, function(chr) {
    GRanges(chr, altchr.mappings[[chr]])
  })

  overlaps = findOverlaps(tailseqs, GRangesList(splitmates.gr))
  as.character(seqnames(tailseqs))[as.matrix(oo)[,1]]
  as.character(seqnames(GRangesList(splitmates.gr)))


  ## alt
  has.islands <- names(splitmate.islands)[sapply(splitmate.islands, function(x) length(x) != 0)]
  islands.gr <- GRangesList(lapply(has.islands, function(chr) {
    GRanges(chr, splitmate.islands[[chr]])
  }))
  
}


### add mapped mate position to d
