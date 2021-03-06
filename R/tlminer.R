## tlminer.R - mine translocations
suppressMessages({
  library(Rsamtools)
  library(ShortRead)
  library(bitops)
  library(multicore)
})

source("R/plotting-methods.R")
source("R/utils.R")

### Config
mll.region <- RangesList(chr11=IRanges(118307205, 118395934))
ref <- "data/mll_template/mll.fasta"
checkBWA(dirname(ref))

### Find BAM file
groupdir <- NULL
TEST.MODE <- system('uname -s', intern=TRUE) == 'Darwin'
# For testing and org-mode usesage.  The exists() check allows us to
# do bamfile = "..." and then source() this.
if (!exists("mapfile")) {
  if (interactive()) {
    message("Running in interactive/test-mode.")
    mapfile <- "data/barcoded_run_01/CAGTACT.sorted.bam" ## test file
    groupdir <- "test"
  } else {
    message("Running in command line mode.")
    cargs <- commandArgs()
    arg.delim <- which(cargs == '--args') + 1
    cargs <- cargs[arg.delim:length(cargs)]
    mapfile <- cargs[1]
    groupdir <- cargs[2]
  }
}

rootname <- getRootname(mapfile)

### Set up directories and check references
dirs <- makeResultsDir(rootname, groupdir=groupdir)
cat(sprintf("Using output directory: %s\n", dirs$base))
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
unmapped <- countBam(mapfile, param=ScanBamParam(flag=scanBamFlag(isUnmappedQuery=TRUE)))
mapped <- countBam(mapfile, param=ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE)))
mapped.mll.forward <- countBam(mapfile, param=ScanBamParam(which=mll.region,
                                          flag=scanBamFlag(isMinusStrand=FALSE,
                                            isUnmappedQuery=FALSE)))
mapped.mll.reverse <- countBam(mapfile, param=ScanBamParam(which=mll.region,
                                          flag=scanBamFlag(isMinusStrand=TRUE,
                                            isUnmappedQuery=FALSE)))
### MLL Region coverage plots



### Split-Mates
# First, we find all BAM entries in which there is a mate that maps to
# the MLL gene. Then, we join these IDs on all others in the database.
mll.forward.param <- ScanBamParam(which=mll.region,
                                  what=c("rname", "pos", "qwidth", "qname", "mrnm", "mpos", "seq", "flag"),
                                  flag=scanBamFlag(isUnmappedQuery=FALSE, isMinusStrand=FALSE))
mll.forward <- scanBam(mapfile, param=mll.forward.param)

# Find all entries with mrnm that is chr11 and remove these.
keep <- which(with(mll.forward[[1]], mrnm != 'chr11'))
splitmates <- lapply(mll.forward[[1]], function(x) x[keep])


### Summarize the split-mates
splitmates <- d <- with(splitmates, data.frame(rname, pos, mrnm, mpos, qwidth)) ## d was used originally as a temp var
splitmates$mpos <- as.numeric(splitmates$mpos)
write.table(splitmates, quote=FALSE, sep="\t", 
            file=file.path(dirs$base, "splitread-candidates.txt"))

#aggregate(d$mpos, list(d$mrnm, d$mpos), length)

### Look for islands
altchr.mappings = local({
  tmp <- split(d, list(d$mrnm))
  lapply(tmp, function(x) {
    IRanges(start=x$mpos, width=x$qwidth)    
  })
})

findIslands <- function(x, minCoverage=40) {
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
mappedunmapped <- scanBam(mapfile, param=mappedunmapped.param)

## Group by chromosome and write to file
ok <- with(mappedunmapped[[1]], {
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

stopifnot(all(ok))

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
    if (FALSE %in% x)
      return(NULL)
    f <- findFusion(x[1], x[2], as.numeric(x[3]))
    unlist(f)
  })
  no.fusion <- sapply(tmp, is.null)
  fusions <- tmp[!no.fusion]
  if (!is.list(fusions) | is.null(fusions))
    next()
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
  if (!system(sprintf('test $(wc -l %s | cut -f1 -d" ") -gt 2', fn))) { ## ! because call returns 1 on error
    ok <- system(sprintf('%s -i %s -o %s -g 1 -d 0 > /dev/null', cdhitcmd, fn, cfn))
    if (ok != 0)
      warning(sprintf("Skipping FASTA file '%s' - error with cd-hit", fasta.file))
  } else {
    warning(sprintf("Skipping FASTA file '%s' - too short (wc <=2)", fasta.file))
  }
  
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
    headers <- sapply(tmp, function(x) x[[1]])
    seqs <- sapply(tmp, function(x) x[[2]])

    ## In rare cases, we could have the same header twice. Rownames
    ## need to be unique, so some are dicarded.
    tbl.headers <- table(headers)
    if (any(tbl.headers > 1)) {
      warning("Cluster sequence headers are not unique.")
      dups <- names(tbl.headers)[tbl.headers > 1]
      remove <- which(dups == rownames(headers))
      seqs <- seqs[-remove]
      headers <- headers[-remove]
    }
    
    tmp <- data.frame(seqs=seqs, stringsAsFactors=FALSE)
    rownames(tmp) <- headers
    tmp
  })

  if (!nrow(clusters))
    next()
  
  # match counts and sequence
  clusters <- merge(clusters, rep.seqs, by.x=0, by.y=0)
  clusters <- data.frame(chr, clusters, stringsAsFactors=FALSE)
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
aln <- scanBam(bamfile, param=tailseqs.aln.param)[[1]]

# Fit the BAM lists into a nice dataframe, extracting the split
# from the qname header info. Note that this contains the columns
# "chr", which is the chromosome the tail sequence mapped to, and
# "mate.chr", which is the chromosome the mapping mate mapped to.
d <- with(aln, {
  mate.chr <- sub(".*chr=(chr[0-9XYM]+).*", "\\1", qname)
  mate.strand <- sub(".*;;mstrand=([\\+\\-]).*", "\\1", qname)
  split <- as.numeric(sub(".*split=([0-9]+).*", "\\1", qname))
  mate.pos <- as.numeric(sub(".*mpos=([0-9]+).*", "\\1", qname))
  data.frame(chr=as.character(rname), start=pos, width=qwidth, strand=strand, seq=as.character(seq), 
             split=split, mate.chr=mate.chr, mate.pos=mate.pos, mate.strand=mate.strand,
             mapq=mapq, cigar=cigar, stringsAsFactors=FALSE)
})

## From the above dataframe, make GRanges of tailsequences
tailseqs <- unlist(local({
  tmp <- apply(d, 1, function(x) {
    GRanges(x[1], IRanges(as.numeric(x[2]), width=as.numeric(x[3])), strand=x[4])
  })
  gr <- GRangesList(tmp)
}))
elementMetadata(tailseqs) = d[, c("split", "mate.chr", "mate.pos", "mate.strand", "mapq", "cigar")]


## Find all split mate islands, covert to GRanges
has.islands <- names(splitmate.islands)[sapply(splitmate.islands, function(x) length(x) != 0)]
if (!length(has.islands)) {
  warning("No islands found with 40 coverage, dropping to 20.")
  splitmate.islands <- findIslands(altchr.mappings, 10)
  splitmate.max.cov <- lapply(splitmate.islands, viewMaxs)
  has.islands <- names(splitmate.islands)[sapply(splitmate.islands, function(x) length(x) != 0)]
  if (!length(has.islands))
    write(sprintf("No islands with 20 coverage; manually run sample '%s'\n", mapfile), stderr())
}

islands <- unlist(GRangesList(lapply(has.islands, function(chr) {
  tmp <- GRanges(chr, splitmate.islands[[chr]])
})))
elementMetadata(islands) <- data.frame(meanCoverage=unlist(sapply(splitmate.islands[has.islands], viewMeans)))

## Find all overlaps between the tail sequences and the split mate
## islands. First, we need to expand the split mate islands (so
## matching is a bit fuzzier).
islands.fuzzy <- islands
buffer <- 400L
start(islands.fuzzy) <- start(islands.fuzzy) - buffer
end(islands.fuzzy) <- end(islands.fuzzy) + buffer

top.candidates <- subsetByOverlaps(tailseqs, islands.fuzzy)
tc.metadata <- local({
  i <- as.numeric(names(top.candidates))
  
  # Use row names to find corersponding metadata in islands 
  tmp <- as.matrix(findOverlaps(tailseqs, islands.fuzzy))
  mc <- elementMetadata(islands)[tmp[match(i, tmp[, 1]), 2],]
  
  mdf <- elementMetadata(top.candidates)
  mdf$splitmate.island.mean.coverage <- mc
  mdf
})
elementMetadata(top.candidates) <- tc.metadata

write.table(as.data.frame(top.candidates), quote=FALSE, sep="\t", 
            file=file.path(dirs$base, "top-candidates.txt"))


write(sprintf("File '%s' complete.", mapfile), stderr())


#### Secondary analysis - splitread technique
# read in all unmapped reads, convert to workable format, save to FASTA for mapping
unmapped <- scanBam(mapfile, param=ScanBamParam(flag=scanBamFlag(isUnmappedQuery=TRUE)))
unmapped.seqs <- as(unmapped[[1]]$seq, "XStringSet")
names(unmapped.seqs) <- unmapped[[1]]$qname
fn <- file.path(dirs$splitread, "unmapped-reads.fasta")
write.XStringSet(unmapped.seqs, file=fn, format="fasta")

message("Running BWA bwasw on unmapped reads (possible splitreads).")
bwarun <- "%s bwasw -T 10 -c 5 -t 3 %s %s > %s 2> /dev/null"

samfile <- file.path(dirs$splitread, "unmapped-reads.sam")
message(sprintf(" - mapping: %s", fasta.file))
system(sprintf(bwarun, bwacmd, ref, fn, samfile))

message(sprintf(" - converting to BAM: %s", samfile))
bamfile <- file.path(dirs$splitread, "unmapped-reads.bam")
system(sprintf("%s view -S -b -o %s %s 2> /dev/null", samtoolscmd, bamfile, samfile))

## Now look for split reads using findFusion
param <- ScanBamParam(simpleCigar=FALSE, flag=scanBamFlag(isUnmappedQuery=FALSE),
                      what=c("qname", "rname", "pos", "strand", "seq", "cigar"))
unmapped.mll.aln <- scanBam(file.path(dirs$splitread, "unmapped-reads.bam"), param=param)[[1]]

fd <- with(unmapped.mll.aln, data.frame(qname=as.character(qname), seq=as.character(seq), cigar, pos, stringsAsFactors=FALSE))
## tmp <- mcapply(fd, 1, function(x) {
##   Cif (FALSE %in% x)
##     return(NULL)
##   f <- findFusion(x[1], x[2], as.numeric(x[3]))
##   unlist(f)
## })

cigar.split = mclapply(fd$cigar, extractCigar, mc.cores=6)

is.splitread <- mclapply(cigar.split, function(x) {
  if (nrow(x) != 2)
    return(FALSE)
  if (x$op[1] == "M" && x$op[2] == "S")
    return(TRUE)
  return(FALSE)
}, mc.cores=6)

keep <- unlist(is.splitread)

fd$mapped <- unlist(mapply(function(is.splitread, seq, cigar) {
  if (!is.splitread)
    return(FALSE)
  substr(seq, 1, cigar$length[1])
}, is.splitread, fd$seq, cigar.split))

fd$unmapped <- unlist(mapply(function(is.splitread, seq, cigar) {
  if (!is.splitread)
    return(FALSE)
  substr(seq, cigar$length[1]+1, cigar$length[1]+cigar$length[2])
}, is.splitread, fd$seq, cigar.split))

fd$break.pos <- unlist(mapply(function(is.splitread, seq, pos, cigar) {
  if (!is.splitread)
    return(-1)
  pos+cigar$length[1]-1
}, is.splitread, fd$seq, fd$pos, cigar.split))

## Now cluster unmapped tail sequences
tail.seqs <- as(fd$unmapped[keep], "XStringSet")
names(tail.seqs) <- fd$qname[keep]
fn <- file.path(dirs$splitread, "tailseqs.fasta")
write.XStringSet(tail.seqs, file=fn, format="fasta")

## Run CD-hit
cfn <- file.path(dirs$splitread, "tailseqs-clusters.fasta")
if (!system(sprintf('test $(wc -l %s | cut -f1 -d" ") -gt 2', fn))) { ## ! because call returns 1 on error
  ok <- system(sprintf('%s -i %s -o %s -g 1 -d 0 > /dev/null', cdhitcmd, fn, cfn))
  if (ok != 0)
    warning(sprintf("Skipping FASTA file '%s' - error with cd-hit", fn))
} else {
  warning(sprintf("Skipping FASTA file '%s' - too short (wc <=2)", fn))
}


## Process clustering results
# Grab representative sequences and corresponding clstr files and process.
clstr.file <- file.path(dirs$splitread, "tailseqs-clusters.fasta.clstr")
  
d <- processClstrFile(clstr.file)

# process .clstr files
rep.seqs <- as.data.frame(cbind(d))

## match qname of clusters with sequence
tmp <-fd[keep, 'unmapped'][match(rownames(rep.seqs), fd$qname[keep])]

# match counts and sequence
clusters <- data.frame(rownames(rep.seqs), rep.seqs, tmp, stringsAsFactors=FALSE)
colnames(clusters) <- c('qname', 'count', 'seq')
rownames(clusters) <- NULL

### Align clustered tail sequences to human genome
fn <- file.path(dirs$splitread, "tailseqs-clusters-for-aln.fasta")
local({
  tmp <- as(clusters$seq, "XStringSet")
  names(tmp) <- clusters$qname
  write.XStringSet(tmp, file=fn, format="fasta")
})
edit.dist <- 8
indel.end <- 0

message("Running BWA aln on cluster sequences (to human genome).")
bwarun <- "%s aln -n %s -i %s %s %s > %s 2> /dev/null"
aln.file <- file.path(dirs$splitread, "cluster-seqs.sai")
system(sprintf(bwarun, bwacmd, edit.dist, indel.end, hg.ref, fn, aln.file))

#bwa samse database.fasta aln_sa.sai short_read.fastq > aln.sam
cluster.sam <- file.path(dirs$splitread, "cluster-seqs.sam")
system(sprintf("%s samse %s %s %s > %s 2> /dev/null",
               bwacmd, hg.ref, aln.file, fn, cluster.sam))

# convert to bam
bamfile <- file.path(dirs$cluster.aln, "cluster-seqs.bam")
system(sprintf("%s view -S -b -o %s %s 2> /dev/null", samtoolscmd, bamfile, cluster.sam))

### Query out overlaps between aligned tail sequences and split-mate regions
tailseqs.aln.param <- ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE), 
                                   what=c("rname", "qname", "pos", "cigar", "qwidth", "mapq", "strand", "seq"))
aln <- scanBam(bamfile, param=tailseqs.aln.param)[[1]]

d <- with(aln, {
  tmp <- match(qname, fd$qname)
  data.frame(qname=qname, start=pos, chr=rname, width=qwidth, strand=strand, seq=as.character(seq), 
             mapq=mapq, cigar=cigar, count=clusters$count[match(qname, clusters$qname)],
             mapped=fd$mapped[tmp], unmapped=fd$unmapped[tmp], break.pos=fd$break.pos[tmp], stringsAsFactors=FALSE)
})

write.table(d, quote=FALSE, sep="\t", file=file.path(dirs$base, "splitread-top-candidates.txt"))
## output column information:
# qname: query name (FASTQ header)
# start: genomic coordinate of the mapped clustered tail sequence
# chr: which chromosome the clustered tail sequence mapped to
# width: mapping width
# strand: mapping strand
# seq: the sequence that mapped (could be RC of sequence of mapped on reverse strand)
# mapq: the mapping quality of the mapped clustered tail sequence
# count: the number of sequences that clustered into this sequence
# mapped: the mapped part of the splitread (FALSE if CIGAR string didn't match xMyS pattern)
# unmapped: the unmapped part (tail sequence) of the splitread (FALSE if CIGAR string didn't match xMyS pattern)
# break.pos: The break position in the MLL template
