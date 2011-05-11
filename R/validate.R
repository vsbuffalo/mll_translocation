## validate.R - Summarize and visualize the MLL grep'ing results.
# After top-candidates.txt are produced by tlminer.R, these results
# are fed into Joe's validation script which forms the translocated
# reads, forms different substrings spanning the hypothesis fusion
# site, and greps for these in the original data. Offset versions are
# created as a way of testing for false positives.
#

library(lattice)
source("R/utils.R")

if (interactive()) {
  validdir <- "data/barcoded_run_01/validation"
  groupdir <- "barcoded_run_01"
} else {
  cargs <- commandArgs()
  arg.delim <- which(cargs == '--args') + 1
  cargs <- cargs[arg.delim:length(cargs)]
  validdir <- cargs[1]
  groupdir <- cargs[2]
}


## Data munging to load all validation results for each barcode

extractMetaInfo <-
# Extracts meta information about the validation from the filename.
function(filename) {
  tmp <- gsub("IPCRtemplate:1-(\\d+)\\.(chr[0-9XY]+):(\\d+)([\\+\\-])\\.output", '\\1;;;\\2;;;\\3;;;\\4', filename)
  parts <- unlist(strsplit(tmp, ";;;"))
  stopifnot(length(parts) == 4)
  names(parts) <- c("mll.breakpoint", "chr.alt", "chr.alt.pos", "chr.alt.strand")
  parts
}

## Build a list of dataframes with attributes that contain information about the validation
barcodes <- list.files(validdir, pattern="[ATCG]+")
d <- lapply(barcodes, function(f) {
  vnames <- list.files(file.path(validdir, f), pattern=".output$")
  vparts <- lapply(vnames, extractMetaInfo)
  names(vparts) <- vnames
  
  tmp <- lapply(vnames, function(x) {
    df <- read.csv(file.path(validdir, f, x), header=TRUE)
    mi <- vparts[[x]] ## get metainfo
    attributes(df) <- c(attributes(df), chr.alt=as.character(mi['chr.alt']), chr.alt.pos=as.integer(mi['chr.alt.pos']),
                        mll.breakpoint=as.integer(mi['mll.breakpoint']), chr.alt.strand=as.character(mi['chr.alt.strand']))
    df
  })

  names(tmp) <- vnames
  tmp
})
names(d) <- barcodes


## Reformat list of lists into single datafram with two keys
dd <- list()
for (b in names(d)) {
  # for each dataframe, add the barcode info
  tmp <- lapply(d[[b]], function(x) {
    atr <- attributes(x)
    x$barcode <- b
    x$chr.alt = atr['chr.alt']
    x$chr.alt.strand = atr['chr.alt.strand']
    x$chr.alt.pos = atr['chr.alt.pos']
    x$mll.breakpoint = atr['mll.breakpoint']
    x
  })
  collapsed <- do.call(rbind, tmp)
  dd[[b]] <- collapsed
}
dd <- do.call(rbind, dd)

rownames(dd) <- NULL

dd <- transform(dd, chr.alt=factor(unlist(chr.alt)),
                chr.alt.strand=factor(unlist(chr.alt.strand)),
                chr.alt.pos=as.integer(unlist(chr.alt.pos)),
                mll.breakpoint=as.integer(unlist(mll.breakpoint)),
                barcode=factor(unlist(barcode)))

## Data visualization
## lapply(names(d), function(b) {
##   lapply(names(d[[b]]), function(v) {
##     odir <- checkDir(file.path("results", groupdir, "validation", b))

##     # we want to remove .output and change to .png
##     fn <- gsub(".output", ".png", v, fixed=TRUE)
   
##     fp <- file.path(odir, fn)
##     png(fp, width=800, height=800, res=100)
##     message(sprintf("creating plot for '%s', saving to '%s'", v, fp))
##     p <- levelplot(count ~ MLL.offset + ChrAlt.offset, data=d[[b]][[v]])
##     print(p)
##     dev.off()
##   })
## })

