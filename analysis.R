
## Prerequistes

## The input data is a SQLite database, produced by =sampe2db.py=. This
## file should have a meaningful name, such as the barcode of a
## particular sample.

suppressMessages({
  require(ggplot2)
  require(Biostrings)
  require(RSQLite)
})
TEST.MODE <- system('uname -s', intern=TRUE) == 'Darwin'
# For testing and org-mode usesage.
if (interactive() && TEST.MODE) {
  dbfile <- "CAGTACT.db"
} else if(!interactive()) {
  args <- commandArgs()
  arg.delim <- which(args == '--args') + 1
  args <- args[arg.delim:length(args)]
  dbfile <- args[1]
}
  
message(sprintf("Establishing database connection to %s", dbfile))
# establish database
drv <- dbDriver("SQLite")
con <- dbConnect(drv, dbname=dbfile)

## ** Utility functions

## These functions are used frequently in the code.

check.dir =
# if a directory doesn't exist, make it
function(dir) {
  if (!file.exists(dir)) {
    message(sprintf("Making directory '%s'", dir))
    system(sprintf("mkdir -p %s", dir))
  }
  dir
}
   
in.path =
 # Use a system call to which to check if a program is in the
 # path. cmd can be a vector of places to look. This is because Emacs
 # on OS X doesn't load $PATH, so it will check all of these places
 # and return a string of the first place that returns a result.
 function(cmd) {
   if (version$minor < 12) {
     paths <- sapply(cmd, function(c) system(sprintf("which %s", c)) == 0)
   } else {
     paths <- sapply(cmd, function(c) system(sprintf("which %s", c), 
                                           ignore.stdout=TRUE) == 0)
   }
   
   if (!any(paths))
     NA
   cmd[min(which(paths))]
 }

check.bwa = 
# check that the reference is properly indexed
function(refdir) {
  index.ext <- unlist(strsplit("amb;;ann;;bwt;;pac;;rbwt;;rpac;;rsa;;sa;;fasta", ';;'))
  contents <- dir(refdir)
  m <- sapply(contents, function(x) {
    tmp <- strsplit(basename(x), '\\.')[[1]]
    tmp[length(tmp)] %in% index.ext
  })
  
  if (!any(m))
    stop(sprintf("Reference in '%s' does not appear to be indexed.", refdir))
}
  
writeFASTA =
# Write a FASTA file, given headers and sequences. Mimics the version
# from Biostrings, but faster
function(x, file="", desc=NULL) {
  if (length(desc) != length(x))
    stop("Arguments for headers and sequences must be same length.")
  con <- file(file, open='w')
  for (j in 1:length(desc)) {
    cat(sprintf(">%s\n%s\n", desc[j], x[j]), file=con)
  }
  close(con)
}

## ** Directory structure

## Everything in this analysis is contained in a directory named after
## the basename of the SQLite file and "-output", i.e. =CAGTACT-output=.

sample.name <- strsplit(basename(dbfile), '\\.')[[1]][1]
outdir <- check.dir(sprintf("%s-output", sample.name))

# make a results directory
results.dir <- check.dir(file.path(outdir, "results"))

## ** Check required programs are in path

message("Checking required programs are in path.")
# check explicit path - allows this to work with org-mode with Emacs
bwa.cmd <- in.path(c("bwa", "/usr/local/bin/bwa"))
py.cmd <- in.path(c("python", "/usr/bin/python"))
cdhit.cmd <- in.path(c("cd-hit", "/usr/local/bin/cd-hit"))

required.cmd <- c(bwa=bwa.cmd, python=py.cmd, cdhit=cdhit.cmd)
if (any(is.na(required.cmd))) {
  failed <- paste(names(required.cmd)[which(is.na(required.cmd))], collapse=',')
  stop(sprintf("The following required programs could not be found in path: %s", failed))
}

## ** Check BWA reference is indexed

message("Checking BWA reference is properly indexed.")
# check the reference has been indexed - this should be in the main directory
refdir <- "mll_template"
ref <- file.path(refdir, "mll.fasta")
check.bwa(refdir)

## ** Raw counts of reads mapped with one forward mate to chr11 and another mate mapped elsewhere.

query <- "
SELECT chr_1, chr_2, strand_1, strand_2, count(*) as count
FROM split_mates
WHERE (chr_1 = 'chr11' OR chr_2='chr11') AND (mqual_1 > 30 AND mqual_2 > 30) 
AND (strand_1 = 'forward')
GROUP BY chr_1, chr_2, strand_1, strand_2;
"
message("Querying split_mate for rearrangement candidate counts.")
all.counts <- dbGetQuery(con, query)

## ** Instate basic count threshold: candidates with more than 10 counts

count.thresh <- 10
counts <- all.counts[all.counts$count > count.thresh, ]
rownames(counts) <- NULL

fn <- file.path(results.dir, "split-mate-candidate-counts.txt")
write.table(counts, fn, quote=FALSE, row.names=FALSE)
print(counts)

## ** Positions of rearrangement candidate reads

## Are there consistent positions of mapped reads in each rearrangement
## candidate? Hierarchical clustering is used to group by distance.

extractCandidates = 
# Given rows from the split_mates table subset for a candidate
# rearrangement (same chr_2, other requirements met), cluster the
# mapped alternate chromosome positions to form clusters of mapped
# reads. Take a subset of these with a mapping count above the
# threshold, extract their position range and total counts.
function(reads.df, clust.member.thresh=2) {
  if (nrow(reads.df) == 0)
    return(NULL)
  pos <- reads.df$pos_2
  names(pos) <- pos
  message("    running hclust() and dist() - this can take a while.")

  pos <- unique(pos)
  names(pos) <- pos
  h = hclust(dist(pos))
  groups <- cutree(h, h=10000)
  groups.counts <- table(groups)
  keep <- groups.counts[groups.counts > clust.member.thresh]
  
  if (length(keep) == 0)
    return(NULL)
  
  
  candidate.pos <- lapply(as.integer(names(keep)), function(x) {
    y <- groups == x
    r <- range(as.integer(names(groups))[y])
    if (any(is.na(r)))
      stop("Range contains NA!")
    return(list(range=r, count=sum(y)))
  })

  return(candidate.pos)
}

# template query for grabbing split mate rows
read.query <- "
SELECT *
FROM split_mates
WHERE chr_1 = 'chr11' AND mqual_1 > 30 AND mqual_2 > 30
AND strand_1 ='forward'
AND chr_2 = '%s' AND strand_2 = '%s';
"

message("Extracting candidate rearrangement positions for:")
# process all candidates from count thresholding step
cands = apply(counts, 1, function(row) {
  message(sprintf("  %s", row[2]))
  d <- dbGetQuery(con, sprintf(read.query, row[2], row[4]))
  return(extractCandidates(d))
})

cands <- cands[which(!is.null(cands))]

names(cands) <- counts$chr_2

## ** Build a mapped mate assembled consensus sequence from mapping positions

## The mates mapped to the translocated sequence from the rearrangement
## chromosome have mapping positions, which can be used to reassemble
## this sequence from the reads. The program =assemble.py= does this
## quickly (as it doesn't need to map to the entire human genome).

## First, we must query all mates mapped in the regions of interest and
## output them to FASTA files.

query <- "
SELECT chr_2, name, seq_2, pos_2 FROM split_mates
WHERE pos_2 >= %s AND pos_2 <= %s AND chr_1 = 'chr11' AND chr_2 = '%s'
AND strand_1 = 'forward' AND mqual_1 > 30 AND mqual_2 > 30;
"

# Make FASTA files for each region of interst's reads. This can take a
# while, so we'll output a message.
message("Writing sequences for consensus assemblies.")
assembly.dir <- check.dir(file.path(outdir, "assembly"))
for (chr in names(cands)) {
  regions <- cands[[chr]]
  for (region in regions) {
    region.name <- paste(region$range, collapse="-")
    fn <- file.path(assembly.dir, sprintf("%s-%s.fasta", chr, region.name))
    results <- dbGetQuery(con, sprintf(query, min(region$range), max(region$range), chr))
    writeFASTA(results$seq_2, file=fn, desc=as.character(results$pos_2))
  }
}

message("Running assemble.py to build consensus assemblies.")
cons.dir <- check.dir(file.path(assembly.dir, "consensuses"))
for (fn in dir(assembly.dir, pattern="\\.fasta$")) {
  chr <- strsplit(fn, '-')[[1]][1]
  out.fn <- file.path(cons.dir, sprintf("%s.fasta", chr))
  system(sprintf("%s assemble.py %s > %s", py.cmd, file.path(assembly.dir, fn), out.fn))
}

## * Split Reads: extracting possible fusion sites and confirming rearrangement partners

## The other information in the paired end reads mapped to the entire
## human genome are those that have one mate mapped and another
## unmapped. BWA's short read aligner (unlike its long read aligner) will
## not align only part of a read. Thus a read containing the fusion site
## somewhere in the middle of its sequence will likely not map, since the
## sequence will contain a large section of translocation chromosome.

## The =unmapped_mates= table contains all reads in which one mate is
## unmapped. Ordering by count, we see evidence of the same rearragement
## partners as with the split-mates data:

query <- "
SELECT mapped_chr, count(*) AS count FROM unmapped_mates 
WHERE mapped_mqual > 30 GROUP BY mapped_chr ORDER BY count DESC;"

message("Querying composition of mapped mates with unmapped partner mate.")
fusion.counts <- dbGetQuery(con, query)

write.table(fusion.counts, file.path(results.dir, "fusion-counts.txt"),
            quote=FALSE, row.names=FALSE)
print(fusion.counts)

## #+results:
## #+begin_example
##    mapped_chr count
## 1       chr11 28218
## 2        chr2  6862
## 3        chr9  3590
## 4       chr21  3535
## 5        chr8  1223
## 6        chr6   767
## 7        chr1   563
## 8        chr5   422
## 9       chr12   414
## 10      chr22   357
## 11       chr4   350
## 12      chr17   295
## 13       chr3   294
## 14      chr16   265
## 15       chr7   238
## 16      chr20   237
## 17      chr14   199
## 18      chr15   160
## 19      chr19   155
## 20      chr18   103
## 21      chr10    94
## 22       chrX    93
## 23      chr13    89
## 24       chrY     3
## #+end_example

## The presumption here is that the unmapped mate will contain some
## chromosome 11 (specifically MLL) sequence. We extract and map the
## unmapped mates, keeping them grouped by the chromosome of their mapped
## mate (which, if this were a true rearrangement, would be the
## rearrangement partner).

query <- "
SELECT mapped_chr as chr, name as header, unmapped_seq as seq 
FROM unmapped_mates WHERE mapped_mqual > 30 ORDER BY mapped_chr;"

# Get all Unmapped mates
message("Querying all unmapped mates.")
unmapped.df <- dbGetQuery(con, query)
unmapped.by.chr <- split(unmapped.df, unmapped.df$chr)

fusion.read.dir <- check.dir(file.path(outdir, "fusion-reads"))

message("Writing unmapped sequences to FASTA files.")
for (chr in names(unmapped.by.chr)) {
  fn <- file.path(fusion.read.dir, sprintf("%s-fusion-candidates.fasta", chr))
  d <- unmapped.by.chr[[chr]]
  writeFASTA(d$seq, fn, desc=d$header)
}

## ** BWA BWASW alignment of unmapped sequences

# run long read aligner on all samples
message("Running BWA bwasw on unmapped mates.")
bwarun <- "%s bwasw -T 10 -c 5 -t 3 %s %s > %s 2> /dev/null"
aln.dir <- check.dir(file.path(fusion.read.dir, "alignments"))
for (fasta.file in dir(fusion.read.dir, pattern="\\.fasta")) {
  chr <- unlist(strsplit(fasta.file, '-'))[1]
  aln.file <- file.path(aln.dir, sprintf("%s.sam", chr))
  system(sprintf(bwarun, bwa.cmd, ref, file.path(fusion.read.dir, fasta.file), aln.file))
}

## ** Processing alignment results with =find_fusion.py=

## Now, we must parse the SAM results and find the fusion sites from
## mapped reads with a CIGAR string of the format *x*M*y*S where *x* and
## *y* are integers and M and S indicate mapped and soft-clipped bases.

## This is done with =find_fusion.py= which uses pysam.

message("Running find_fusion.py on BWA bwasw results.")
system(sprintf("ls %s/*sam | xargs -n1 %s find_fusion.py", aln.dir, py.cmd))

## ** Statistical analysis of fusion sites

## We load each of these alignment files into the =hybid_candidates= table.

tbl.name <- "hybrid_candidates"

message("Loading find_fusion results into database table.")
# Remove any existing tables; otherwise we could load duplicates.
if (dbExistsTable(con, tbl.name))
  dbRemoveTable(con, tbl.name)

# Build a table
cols <- c(chr='text', name='text', split='integer', mapped='text',
          softclipped='text', strand='text', mqual='integer')
tbl.query <- dbBuildTableDefinition(drv, tbl.name, NULL, field.types=cols)
dbGetQuery(con, tbl.query)

# Load each dataframe into table
for (f in dir(aln.dir, pattern="fusion-candidates\\.txt")) {
  chr <- strsplit(f, '-')[[1]][1]

  # wrap in try because some files will be completely empty.
  d <- try({read.csv(file.path(aln.dir, f), header=FALSE, sep='\t')}, TRUE)
  if (is(d, 'try-error') || !nrow(d))
    next()
  d <- cbind(chr, d)
  ok <- dbWriteTable(con, tbl.name, d, append=TRUE, row.names=FALSE)
  stopifnot(ok)
}

## What does the distribution of split points look like? Chromosome 11 is
## presently excluded because its level of counts are much higher than
## those of the other chromosomes.

query <- "SELECT chr AS chromosome, split, count(*) AS count
FROM hybrid_candidates
WHERE strand = 'forward'
GROUP BY chr, split
HAVING count > 20 ORDER BY count DESC;"

split.df <- dbGetQuery(con, query)
p <- ggplot(subset(split.df, chromosome!='chr11'), aes(x=split))
p <- p + geom_histogram(aes(y=count, fill=chromosome), size=3, stat="identity", position='dodge')
p <- p + scale_fill_brewer() #+ scale_x_continuous(limit=c(200, 250))
p <- p + xlab("position") + ylab("count") 
if (interactive()) {
  print(p)
} else {
  png(file.path(results.dir, "break-points.png"), width=900, height=900)
  print(p)
  dev.off()
}

## #+results:

## The chromosomes above are all rearrangement candidates. Now, write
## FASTA files for each of these (again, for the moment excluding chr11) to
## cluster.

split.cands <- split.df$chromosome[split.df$chromosome!='chr11']

query <- "
SELECT chr, name, softclipped FROM hybrid_candidates
WHERE chr IN (%s) GROUP BY name;"

message("Querying fusion candidate hybrid reads from database, writing FASTA files.")
tmp <- sapply(split.cands, function(x) sprintf("'%s'", x))
seqs <- dbGetQuery(con, sprintf(query, paste(tmp, collapse=', ')))

# divide by chromosome, write FASTA files
seqs.by.chr <- split(seqs, seqs$chr)

cluster.dir <- check.dir(file.path(outdir, "clusters"))
for (chr in names(seqs.by.chr)) {
  fn <- file.path(cluster.dir, sprintf("%s-clipped.fasta", chr))
  s <- seqs.by.chr[[chr]]
  if (nrow(s) > 0)
    writeFASTA(s$softclipped, fn, desc=s$name)
}

## ** Soft-clipped sequence clustering

## Use =cd-hit= to cluster soft-clipped sequences (which in a
## rearrangement will be the alternate chromosome). =cd-hit= produces
## FASTA files of representative sequences, as well as .clstr files that
## indicate cluster membership. To see how many sequences are clustered
## into a single representative sequence, we extract the information from
## the .clstr file.

message("Running cd-hit on fusion candidate reads.")
for (fasta.file in dir(cluster.dir, pattern="\\-clipped.fasta$")) {
  chr <- strsplit(fasta.file, '-')[[1]][1]
  fn <- file.path(cluster.dir, fasta.file)
  cfn <- file.path(cluster.dir, sprintf("%s-clusters.fasta", chr))
  ok <- system(sprintf("%s -i %s -o %s -g 1 -d 200 > /dev/null", cdhit.cmd, fn, cfn))
  stopifnot(ok == 0)
}

processClstrFile =
# cd-hit produces .clstr files with information on the cluster
# sequences. This is the only way to (1) get the number of sequences
# clustered in a representative sequence, and (2) get the header of
# the representative sequence.
function(filename) {
  contents <- readLines(filename)
  clusters <- list()
  first <- TRUE
  for (line in contents) {    
    if (length(grep('>Cluster', line))) {
      if (!first) {
        # push result to clusters list
        clusters[[header]] <- count
      } else {
        first <- FALSE
      }
      count <- 0
    } else {
      count <- count + 1
      parts <- unlist(strsplit(line, '\\s+'))
      if ('*' %in% parts) {
        header <- gsub('[^>]+>([^\\.]+).*', '\\1', line)
      }
    }
  }

  return(unlist(clusters))
}

all.clusters <- list()
message("Processing cd-hit clustering results.")
# Grab representative sequences and corresponding clstr files and process.
for (fasta.file in dir(cluster.dir, pattern="\\-clusters.fasta$")) {
  chr <- strsplit(fasta.file, '-')[[1]][1]
  fn <- file.path(cluster.dir, fasta.file)
  clstr.file <- file.path(cluster.dir, sprintf("%s.clstr", fasta.file))

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
  all.clusters[[chr]] <- clusters
  print(clusters[clusters$count > 10, c('chr', 'seq', 'count')])
}

# bit of name mangling here...
clusters <- do.call(rbind, all.clusters)
rownames(clusters) <- NULL

## * Final Output

## Now, output all results. The =basename-output= directory contains:

##  - =assembly=
##    - =*.fasta*= files of the mates in which the other mate maps to
##      chromosome 11 (forward strand), with headers of the mapping
##      position. We may want to just use the mapped mate in split-read
##      mates.
##    - =consensuses= assemblies of the mates done by position.
##  - =fusion-reads=
##    - =*.fasta= files containing the unmapped mates (mapped mate
##      chromosome is in the filename).
##    - =alignments=
##      - =*.sam= =BWA bwasw= alignment of the unmapped mates to the
##        reference (MLL) template. These should produce
##        partially-aligned sequences, with soft-clipped sequence being
##        that of the rearrangement partner in the case that the read
##        spans the fusion site.
##      - =*-fusion-candidates.txt= tab delimited output from
##        =find_fusion.py=, which parses the SAM mapping file to find
##        mappings with the cigar string in a *x*M*y*S format, which
##        corresponds to a soft-clipped 3'-end sequence. The strand in
##        the translocated case should be 'forward', and the position
##        around 215.
##      - =*-hybrids.fasta= The soft-clipped portions of alignment, used
##        for clustering.
##  - =clusters=
##    - =*-clusters.fasta= clustered (with =cd-hit= tail
##      sequences). Longer sequences from the rearrangement partner
##      should BLAT to the Human Genome, with a match in the same
##      location as the mate's chromosome.
##    - =*-clusters.fasta.clstr= Cluster summary file produced by
##      =cd-hit=. This contains representative sequences, and those
##      sequences clustered with them. This file is used to get a count
##      of the cluster density. Translocation partners should have high
##      cluster density.
##    - =*-clipped.fasta= duplicates =*-hybrids.fasta=.

sm.cands <- local({
  # Make ragged list cands more output-friendly
  tmp <- sapply(cands, function(x) do.call(rbind, x))
  d <- do.call(rbind, lapply(names(tmp), function(x) {
    d <- tmp[[x]]
    ranges <- do.call(rbind, d[, 1])
    counts <- do.call(rbind, d[, 2])
    cbind(x, ranges, counts)
  }))
  
  d <- as.data.frame(d, row.names=FALSE)
  colnames(d) <- c('chromosome', 'lower.pos', 'upper.pos', 'count')
  
  fn <- file.path(results.dir, 'split-mates-candiatates.txt')
  write.table(as.matrix(d), fn, row.names=FALSE, quote=FALSE, sep='\t')
  d
})
print(sm.cands)

## Output clustered sequences with mapped mate and counts.

cluster.cands <- local({
  d <- clusters[clusters$count > 10, c('chr', 'seq', 'count')]
  fn <- file.path(results.dir, "clustered-seqs.txt")
  write.table(as.matrix(d), fn, row.names=FALSE, quote=FALSE, sep='\t')
  d
})
print(cluster.cands)
