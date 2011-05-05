## utils.R - utilities for tlminer.R

check.dir =
# if a directory doesn't exist, make it
function(dir) {
  if (!file.exists(dir)) {
    message(sprintf("Making directory '%s'", dir))
    system(sprintf("mkdir -p %s", dir))
  }
  dir
}

getRootname =
# Given a file, gets the part of the basename before the extension
function(file)
  sapply(strsplit(basename(file), ".", fixed=TRUE), function(x) x[1])

makeResultsDir =
# make the directory structure
function(rootname) {
  base <- check.dir(sprintf("%s-run", rootname))
  unmapped.mates <- check.dir(file.path(base, "unmapped-mate"))
  aln <- check.dir(file.path(base, "unmapped-mate-alns"))
  tailseqs <- check.dir(file.path(base, "tailseqs"))
  cluster <- check.dir(file.path(base, "tailseqs-clustered"))
  cluster.aln <- check.dir(file.path(base, "clustered-tailseqs-aln"))
  results <- check.dir(file.path(base, "results"))
  
  list(base=base, unmapped.mates=unmapped.mates, aln=aln,
       tailseqs=tailseqs, cluster=cluster, cluster.aln=cluster.aln,
       results=results)
}

 
checkBWA = 
# check that the reference is properly indexed
function(refdir) {
  index.ext <- unlist(strsplit("amb;;ann;;bwt;;pac;;rbwt;;rpac;;rsa;;sa;;fasta", ';;'))
  contents <- dir(refdir)
  m <- sapply(contents, function(x) {
    tmp <- strsplit(basename(x), '\\.')[[1]]
    tmp[length(tmp)]
  })
  
  if (any(is.na(match(index.ext, m))))
    stop(sprintf("Reference in '%s' does not appear to be indexed.", refdir))
}   

inPath =
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

extractCigar =
function(cigar) {
  nums <- as.numeric(unlist(strsplit(cigar, '[MIDNSHPX=]')))
  type <- unlist(strsplit(cigar, '[0-9]+'))
  data.frame(op=type[-1], length=nums, stringsAsFactors=FALSE)
}

findFusion =
# Split a read into the mapped and unmapped part, returning chunks of
# sequence. Only reads with cigar strings in the form xMyS or ySxM are
# handled.
function(seq, cigar, pos) {
  cigar = extractCigar(cigar)
  if (nrow(cigar) != 2)
    return(NULL)
  if (cigar$op[1] == "M" && cigar$op[2] == "S")
    return(list(mapped=substr(seq, 1, cigar$length[1]),
                unmapped=substr(seq, cigar$length[1]+1, cigar$length[1]+cigar$length[2]),
                break.pos=pos+cigar$length[1]-1))
  return(NULL)
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
