## splitread-miner.R - mine translocations from splitreads
suppressMessages({
  library(Rsamtools)
  library(bitops)
})

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
message(sprintf(" - mapping: %s", fn))
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
