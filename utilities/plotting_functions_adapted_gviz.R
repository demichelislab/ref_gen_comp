extract_index_mism_reads <- function(sequence_track, alignment_track, region, pos, var_base){
  ref <- as.character(as(subseq(sequence_track, start = start(region), end = end(region)), "Rle"))
  cm <- consensusMatrix(alignment_track@sequences, as.prob = FALSE, baseOnly = TRUE)[-5, ]
  cmm <- colMaxs(cm)
  css <- colSums(cm)
  cmp <- rbind(t(t(cm) / css), 0)
  rownames(cmp)[5] <- "N"
  sel <- is.na(cmp["A", ])
  cmp[, sel] <- 0
  cmp["N", sel] <- 1
  consStr <- strsplit(consensusString(cmp), "")[[1]]
  varRegs <- which(cmm != css | (consStr != "N" & consStr != ref))
  rvg <- ref[varRegs]
  sel <- rvg != "-" & rvg != "N"
  varRegs <- varRegs[sel]
  rvg <- rvg[sel]
  mmTab <- do.call(rbind, lapply(varRegs, function(x) as.character(subseq(alignment_track@sequences, x, width = 1))))
  isMm <- t(rvg != "-" & mmTab != "+" & mmTab != "-" & mmTab != rvg)
  mmRelPos <- col(isMm)[which(isMm)]
  mmPos <- varRegs[mmRelPos] + start(region) - 1
  mmSampInd <- row(isMm)[which(isMm)]
  mmSamp <- rownames(isMm)[mmSampInd]
  mmSeq <- mmTab[ncol(isMm) * (mmSampInd - 1) + mmRelPos]
  total_bases <- sum(cm[,end(region) - pos +1])
  print(paste0("Coverage: ", total_bases))
  id_reads <- mmSampInd[mmPos == pos & mmSeq == var_base]
  print(paste0("VAF:", round((length(id_reads) / total_bases) * 100, 2), "%"))
  
  return(id_reads)
}

plot_coverage_track <- function(path_bam, bs_genome, chrom, pos, width, lwd_line = 0.7){
  sTrack <- SequenceTrack(bs_genome)
  alTrack <- AlignmentsTrack(path_bam)
  
  gTrack <- GenomeAxisTrack()
  
  displayPars(gTrack) <- list( ticksAt = c(pos-width+10, pos, pos + width-10),
                               #labelPos = "below",
                               littleTicks = F,
                               lwd = 0.3,
                               cex = 0.8)
  displayPars(sTrack) <- list(  lwd = 0.05,
                                min.width = 0)
  
  displayPars(alTrack) <- list(background.title = "white", 
                               lwd.coverage = 0.8,
                               col.axis = "black",
                               cex.axis = 0.8,
                               lwd.axis = 0.1,
                               showTitle = F,
                               col = NULL,
                               cex = 0.8,
                               min.height = 0,
                               fill.coverage = "white",
                               showIndels = T,
                               type = "coverage"
  )
  
  htrack <- HighlightTrack(trackList = list(gTrack, alTrack), start = pos, end = pos, chromosome = chrom,
                           inBackground = F, fill = NA, col = "#a12d69", lwd =  lwd_line)
  
  plotTracks(list(htrack), from = pos - width , to = pos + width, 
             chromosome = chrom,
             #margin = c(-0.5, 3, 0, 0),
             margin = 20,  # Left margin for y-axis (adjust as needed)
             innerMargin = 0,
             add = T,
             sizes=c(2, 4))
  
}

plot_variant_reads <- function(path_bam, bs_genome, chrom, pos, width, defined_mar = c(0, 40, 0, 0) ,
                               lwd_line = 0.7, random_id_ns = 0){
  bf <- BamFile(path_bam, asMates = TRUE)
  coords <- paste0(chrom, ":", pos - width, "-", pos + width)
  selection <-  GRanges(coords)
  param <- ScanBamParam(which = selection, 
                        what = scanBamWhat(), 
                        flag = scanBamFlag(isUnmappedQuery = FALSE))  #isSecondaryAlignment =F), mapqFilter=60)
  zz <- scanBam(bf, param = param)[[1]]
  
  layed_seq <- sequenceLayer(zz$seq, zz$cigar)
  region <- unlist(bamWhich(param), use.names = FALSE)
  reads <- if (as.character(seqnames(selection)[1]) %in% names(scanBamHeader(bf)$targets)) scanBam(bf, param = param)[[1]] else list()
  ans <- stackStrings(layed_seq, start(region), end(region), shift = reads$pos - 1L, Lpadding.letter = "+", Rpadding.letter = "+")
  names(ans) <- seq_along(reads$qname)
  
  gr <- GRanges(zz$rname, IRanges(zz$pos, width = zz$qwidth), zz$strand, 
                id = zz$qname, cigar = zz$cigar, mapq = zz$mapq, flag = zz$flag, isize = zz$isize, 
                groupid = zz$groupid, status = zz$mate_status)
  
  sTrack <- SequenceTrack( bs_genome, chromosome = chrom, start = start(region), end = end(region))
  full_alTrack <- AlignmentsTrack(gr, referenceSequence = sTrack, seqs = ans,  showIndels = T)
  ids_mism <- extract_index_mism_reads(sTrack, full_alTrack, region, pos)
  # to_pick <- 1:length(gr)[-ids_mism]
  # set.seed(123)
  # random_ids <- sample(to_pick, random_id_ns)
  # subset_reads <- gr$id[c(ids_mism, random_ids)]
  subset_reads <- gr$id[c(ids_mism)]
  
  z <- lapply(zz, function(y) y[zz$qname %in% subset_reads])
  layed_seq <- sequenceLayer(z$seq, z$cigar)
  region <- unlist(bamWhich(param), use.names = FALSE)
  reads <- if (as.character(seqnames(selection)[1]) %in% names(scanBamHeader(bf)$targets)) scanBam(bf, param = param)[[1]] else list()
  reads <- lapply(zz, function(y) y[zz$qname %in% subset_reads])
  ans <- stackStrings(layed_seq, start(region), end(region), shift = reads$pos - 1L, Lpadding.letter = "+", Rpadding.letter = "+")
  names(ans) <- seq_along(reads$qname)
  
  ## Now turn that into a GRanges object
  gr <- GRanges(z$rname, IRanges(z$pos, width = z$qwidth), z$strand, 
                id = z$qname, cigar = z$cigar, mapq = z$mapq, flag = z$flag, isize = z$isize, 
                groupid = z$groupid, status = z$mate_status)
  
  sTrack <- SequenceTrack( bs_genome, chromosome = chrom, start = start(region), end = end(region))
  alTrack <- AlignmentsTrack(gr, referenceSequence = sTrack, seqs = ans)
  
  gTrack <- GenomeAxisTrack()
  
  displayPars(gTrack) <- list(labelPos = "below",
                              littleTicks = F,
                              lwd = 0.3,
                              cex = 0.3)
  
  displayPars(sTrack) <- list(lwd = 0.05,
                              min.width = 0)
  
  displayPars(alTrack) <- list(background.title = "white", 
                               showTitle = F,
                               col = NULL,
                               cex = 0.8,
                               min.height = 20,
                               showIndels = T,
                               max.height = 20,
                               type = "pileup")
  
  
  htrack <- HighlightTrack(trackList = list(alTrack, sTrack), start = pos, end = pos, chromosome = chrom,
                           inBackground = F, fill = NA, col = "#a12d69", lwd = lwd_line)
  
  plotTracks(list(htrack), from = start(region) , to = end(region),
             chromosome = chrom,
             
             #margin = 40,  # Left margin for y-axis (adjust as needed)
             #innerMargin = 0,
             margin = defined_mar,
             sizes=c(5,2), add = T)
}

plot_variant_reads <- function(path_bam, bs_genome, chrom, pos, width, var_base, defined_mar = c(0, 45.5, -17, 16.3),
                               lwd_line = 0.7, random_id_ns = 0){
  bf <- BamFile(path_bam, asMates = TRUE)
  coords <- paste0(chrom, ":", pos - width, "-", pos + width)
  selection <-  GRanges(coords)
  param <- ScanBamParam(which = selection, 
                        what = scanBamWhat(), 
                        flag = scanBamFlag(isUnmappedQuery = FALSE,  isSecondaryAlignment = F,
                                           isNotPassingQualityControls = F), mapqFilter=60L)
  zz <- scanBam(bf, param = param)[[1]]
  
  layed_seq <- sequenceLayer(zz$seq, zz$cigar)
  region <- unlist(bamWhich(param), use.names = FALSE)
  reads <- if (as.character(seqnames(selection)[1]) %in% names(scanBamHeader(bf)$targets)) scanBam(bf, param = param)[[1]] else list()
  ans <- stackStrings(layed_seq, start(region), end(region), shift = reads$pos - 1L, Lpadding.letter = "+", Rpadding.letter = "+")
  names(ans) <- seq_along(reads$qname)
  
  gr <- GRanges(zz$rname, IRanges(zz$pos, width = zz$qwidth), zz$strand, 
                id = zz$qname, cigar = zz$cigar, mapq = zz$mapq, flag = zz$flag, isize = zz$isize, 
                groupid = zz$groupid, status = zz$mate_status)
  
  sTrack <- SequenceTrack( bs_genome, chromosome = chrom, start = start(region), end = end(region))
  full_alTrack <- AlignmentsTrack(gr, referenceSequence = sTrack, seqs = ans,  showIndels = T)
  ids_mism <- extract_index_mism_reads(sTrack, full_alTrack, region, pos, var_base)
  #to_pick <- 1:length(gr)[-ids_mism]
  #set.seed(123)
  #random_ids <- sample(to_pick, random_id_ns)
  subset_reads <- gr$id[c(ids_mism)]
  print(subset_reads)
  
  z <- lapply(zz, function(y) y[zz$qname %in% subset_reads])
  layed_seq <- sequenceLayer(z$seq, z$cigar)
  region <- unlist(bamWhich(param), use.names = FALSE)
  reads <- if (as.character(seqnames(selection)[1]) %in% names(scanBamHeader(bf)$targets)) scanBam(bf, param = param)[[1]] else list()
  reads <- lapply(zz, function(y) y[zz$qname %in% subset_reads])
  ans <- stackStrings(layed_seq, start(region), end(region), shift = reads$pos - 1L, Lpadding.letter = "+", Rpadding.letter = "+")
  names(ans) <- seq_along(reads$qname)
  
  ## Now turn that into a GRanges object
  gr <- GRanges(z$rname, IRanges(z$pos, width = z$qwidth), z$strand, 
                id = z$qname, cigar = z$cigar, mapq = z$mapq, flag = z$flag, isize = z$isize, 
                groupid = z$groupid, status = z$mate_status)
  
  sTrack <- SequenceTrack( bs_genome, chromosome = chrom, start = start(region), end = end(region))
  alTrack <- AlignmentsTrack(gr, referenceSequence = sTrack, seqs = ans)
  
  z <- lapply(zz, function(y) y[!zz$qname %in% subset_reads])
  layed_seq <- sequenceLayer(z$seq, z$cigar)
  region <- unlist(bamWhich(param), use.names = FALSE)
  reads <- if (as.character(seqnames(selection)[1]) %in% names(scanBamHeader(bf)$targets)) scanBam(bf, param = param)[[1]] else list()
  reads <- lapply(zz, function(y) y[!zz$qname %in% subset_reads])
  ans <- stackStrings(layed_seq, start(region), end(region), shift = reads$pos - 1L, Lpadding.letter = "+", Rpadding.letter = "+")
  names(ans) <- seq_along(reads$qname)
  
  ## Now turn that into a GRanges object
  gr <- GRanges(z$rname, IRanges(z$pos, width = z$qwidth), z$strand, 
                id = z$qname, cigar = z$cigar, mapq = z$mapq, flag = z$flag, isize = z$isize, 
                groupid = z$groupid, status = z$mate_status)
  
  collapse_al <- AlignmentsTrack(gr, referenceSequence = sTrack, seqs = ans)
  
  gTrack <- GenomeAxisTrack()
  
  displayPars(gTrack) <- list(labelPos = "below",
                              littleTicks = F,
                              lwd = 0.3,
                              cex = 0.3)
  
  displayPars(sTrack) <- list(lwd = 0.05,
                              min.width = 0)
  
  displayPars(alTrack) <- list(background.title = "white", 
                               showTitle = F,
                               col = NULL,
                               cex = 0.8,
                               min.height = 10,
                               showIndels = T,
                               max.height = 10,
                               type = "pileup"
  )
  
  
  displayPars(collapse_al) <- list(background.title = "white", 
                                   lwd.coverage = 0.8,
                                   col.axis = "black",
                                   cex.axis = 0.7,
                                   lwd.axis = 0.1,
                                   showTitle = F,
                                   col = NULL,
                                   cex = 0.8,
                                   min.height = 0,
                                   max.height = 1,
                                   fill.coverage = "white",
                                   showIndels = T,
                                   type = "pileup"
  )
  displayPars(full_alTrack) <- list(background.title = "white", 
                                    lwd.coverage = 0.8,
                                    col.axis = "black",
                                    cex.axis = 0.8,
                                    lwd.axis = 0.1,
                                    showTitle = F,
                                    col = NULL,
                                    cex = 0.8,
                                    min.height = 0,
                                    fill.coverage = "white",
                                    showIndels = T,
                                    type = "coverage"
  )
  # htrack <- HighlightTrack(trackList = list(alTrack, sTrack), start = pos, end = pos, chromosome = chrom,
  #                        inBackground = F, fill = NA, col = "#a12d69", lwd = lwd_line)
  # 
  # plotTracks(list(htrack), from = start(region) , to = end(region), 
  #          chromosome = chrom,
  #          margin = c(0, 0, 0, 0) ,
  #          sizes=c(5,2))
  
  htrack <- HighlightTrack(trackList = list(alTrack, collapse_al, sTrack), start = pos, end = pos, chromosome = chrom,
                           inBackground = F, fill = NA, col = "#a12d69", lwd = lwd_line ,
                           stackHeight=1, shape="box")
  
  plotTracks(list(htrack), from = start(region) , to = end(region), 
             chromosome = chrom,
             margin = defined_mar,
             sizes=c(5,10,2),
             add = T
           #  margin = 20,  # Left margin for y-axis (adjust as needed)
          #   innerMargin = 0
  )
}

plot_all_reads <- function(path_bam, bs_genome, chrom, pos, width){
  sTrack <- SequenceTrack(bs_genome)
  alTrack <- AlignmentsTrack(path_bam)
  
  gTrack <- GenomeAxisTrack()
  
  displayPars(gTrack) <- list( labelPos = "below",
                               littleTicks = F,
                               lwd = 0.3,
                               cex = 0.3)
  displayPars(sTrack) <- list(  lwd = 0.05,
                                min.width = 0)
  
  displayPars(alTrack) <- list(background.title = "white", 
                               lwd.coverage = 0.8,
                               col.axis = "black",
                               cex.axis = 0.7,
                               lwd.axis = 0.1,
                               showTitle = F,
                               col = NULL,
                               cex = 0.8,
                               min.height = 0,
                               fill.coverage = "white",
                               showIndels = T)
  
  htrack <- HighlightTrack(trackList = list(alTrack, sTrack), start = pos, end = pos, chromosome = chrom,
                           inBackground = F, fill = NA, col = "#a12d69", lwd = 0.7)
  
  plotTracks(list(gTrack, htrack), from = pos - width , to = pos + width, 
             chromosome = chrom,
             #margin = c(0, 0, -15, 0) ,
             margin = 20,  # Left margin for y-axis (adjust as needed)
             innerMargin = 0,
             sizes=c(1,12,1))
}

plot_compact_all_reads <- function(path_bam, bs_genome, chrom, pos, width){
  sTrack <- SequenceTrack(bs_genome)
  alTrack <- AlignmentsTrack(path_bam)
  
  gTrack <- GenomeAxisTrack()
  
  displayPars(gTrack) <- list( labelPos = "below",
                               littleTicks = F,
                               lwd = 0.3,
                               cex = 0.3)
  displayPars(sTrack) <- list(  lwd = 0.05,
                                min.width = 0)
  
  displayPars(alTrack) <- list(background.title = "white", 
                               lwd.coverage = 0.8,
                               col.axis = "black",
                               cex.axis = 0.7,
                               lwd.axis = 0.1,
                               showTitle = F,
                               col = NULL,
                               cex = 0.8,
                               min.height = 0,
                               max.height = 1,
                               fill.coverage = "white",
                               showIndels = T,
                               type = "pileup")
  
  htrack <- HighlightTrack(trackList = list(alTrack, sTrack), start = pos, end = pos, chromosome = chrom,
                           inBackground = F, fill = NA, col = "#a12d69", lwd = 0.7)
  
  plotTracks(list(htrack), from = pos - width , to = pos + width, 
             chromosome = chrom,
             #margin = c(0, 0, -15, 0) ,
             add = T,
             sizes=c(1,1),
             margin = 20,  # Left margin for y-axis (adjust as needed)
             innerMargin = 0)
}