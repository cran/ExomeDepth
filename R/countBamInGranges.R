

countBamInGRanges.exomeDepth <- function (bam.file, granges, min.mapq = 1, read.width = 1)  {

  rds.counts <- numeric(length(granges))
  seq.names <- seqlevels(granges)
  seq.names.in.bam <- names(scanBamHeader(bam.file)[[1]]$targets)
        
  for (seq.name in seq.names) {
    if (seq.name %in% seq.names.in.bam) {
      granges.subset <- granges[seqnames(granges) == seq.name]
      strand(granges.subset) <- "*"

      empty <- TRUE
      ############################################################################# read paired end 
      rds <- scanBam(bam.file,
                     param = ScanBamParam(flag = scanBamFlag(isDuplicate = FALSE, isPaired = TRUE, isProperPair = TRUE), what = c("pos", "mpos", "mapq"), which = range(granges.subset)))
      mapq.test <- (rds[[1]]$mapq >= min.mapq) & !is.na(rds[[1]]$pos) & (abs(rds[[1]]$mpos - rds[[1]]$pos) < 1000)
      
      if (sum(mapq.test) > 0 && !is.na(sum(mapq.test) )) {
        empty <- FALSE
        rds.ranges <- GRanges(seq.name, IRanges(start = 0.5*(rds[[1]]$pos[mapq.test] + rds[[1]]$mpos[mapq.test]) - read.width, width = read.width))
        rds.counts.seq.name <- countOverlaps(granges.subset, rds.ranges)
        rds.counts[as.logical(seqnames(granges) == seq.name)] <- rds.counts.seq.name
      }

      ############################################################################# read single end 
      rds <- scanBam(bam.file,
                     param = ScanBamParam(flag = scanBamFlag(isPaired = FALSE), what = c("pos", "mapq"), which = range(granges.subset)))
      mapq.test <- (rds[[1]]$mapq >= min.mapq) & !is.na(rds[[1]]$pos) 
      
      if (sum(mapq.test) > 0 && !is.na(sum(mapq.test))) {
        empty <- FALSE
        rds.ranges <- GRanges(seq.name, IRanges(start = rds[[1]]$pos[mapq.test] - read.width, width = read.width))
        rds.counts.seq.name <- countOverlaps(granges.subset, rds.ranges)
        rds.counts[as.logical(seqnames(granges) == seq.name)] <- rds.counts.seq.name
      }
######################
      if (empty) rds.counts[as.logical(seqnames(granges) == seq.name)] <- 0  ## do I need that?
      #message('Sequence ', seq.name, ' ',  rds.counts[as.logical(seqnames(granges) == seq.name)], '\n')
    }
  }
  rds.counts
}



getBamCounts <- function(bed.frame = NULL, bed.file = NULL, bam.files, min.mapq = 20, read.width = 300, include.chr = FALSE, referenceFasta = NULL) {
  require(GenomicRanges)
  require(Rsamtools)
  
  if (is.null(bed.frame)) {
    if (is.null(bed.file)) {
      stop("If no bed data frame is provided there must be a link to a bed file")
    }
    bed.frame <- read.delim(file = bed.file, header =  FALSE, stringsAsFactors = FALSE)
  }

  names(bed.frame)[1] <- 'seqnames'
  names(bed.frame)[2] <- 'start'
  names(bed.frame)[3] <- 'end'

  if (include.chr) bed.frame$seqnames <- paste('chr', bed.frame$seqnames, sep = '')
  chr.names.used <- unique(as.character(bed.frame$seqnames))
  chr.levels <- c(as.character(seq(1, 22)), subset( chr.names.used, ! chr.names.used %in% as.character(seq(1, 22))))

  bed.frame$seqnames <- factor(bed.frame$seqnames, levels = chr.levels)  ####specifying the levels is important here to not mess up the order
  bed.frame <- bed.frame[ order(bed.frame$seqnames, bed.frame$start + bed.frame$end), ]  ##order the data frame by position

  target <- GRanges(seqnames = bed.frame$seqnames,  
                    IRanges(start=bed.frame$start+1,end=bed.frame$end))
  
                                        #if (exomeCopy) {target <- subdivideGRanges(target)}  ##exomeCopy recommends some splitting
  
  rdata <- RangedData(space=seqnames(target),
                      ranges=ranges(target))

  if  ((ncol(bed.frame) >= 4) && (class(bed.frame[,4]) %in% c('character', 'factor'))) {    
    row.names(rdata) <- make.unique(bed.frame[,4])  ##add exon names if available
  }
  
############################################################################# add GC content
  if (!is.null(referenceFasta)) {
    message('Reference fasta file provided so exomeDepth will compute the GC content in each window')
    target.dnastringset <- scanFa(referenceFasta, target)
  
    getGCcontent <- function(x) {
      GC.count <- letterFrequency(x,"GC")
      all.count <- letterFrequency(x,"ATGC")
      as.vector(ifelse(all.count==0,NA,GC.count/all.count))
    }
    rdata[["GC"]] <- getGCcontent(target.dnastringset)
  }

############################################################################# Parse BAM files
  message('Parse BAM files')
  for (bam in bam.files) {
    message('Parsing ', bam)
    rdata[[ basename(bam) ]] <- countBamInGRanges.exomeDepth(bam,target, min.mapq = min.mapq, read.width = read.width)
  }

  #if (!exomeCopy) row.names(rdata) <- make.unique(bed.frame$names)
  return(rdata)
}
  
