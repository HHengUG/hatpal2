#' Map APAc with prebed to count region
#'
#' @param RE_prebed RE file of one strand generated from ExtractAPA()
#' @param SA_prebed SA file of one strand generated from ExtractAPA()
#' @param APAc_file APA clusters file from FilterAPAc()
#' @param out_dir The directory where the output files are stored
#' @param blocksize Threshold for isolated nodes and clusters (default 250)
#' @param APAdist_cut Threshold for clustering the APA (default 25)
#' @param SA_weight The weight increase of SA file (default 3)
#'
#' @return NULL, writes the output "APA_map.out" in the out_dir.
#' @export
#'
#' @import data.table Ckmeans.1d.dp
#' @rawNamespace import(GenomicAlignments, except = c(first, second, last))
#'
#' @examples MapAPAc("rc_end.prebed", "sa.prebed", "APA_clusters.filtered.out", "output/")
MapAPAc <- function(RE_prebed, SA_prebed, APAc_file, out_dir, blocksize = 250, APAdist_cut = 25, SA_weight = 3){


  ### load and bind prebeds
  print("MapAPAc: load prebed and APAc")
  cat("MapAPAc: load prebed and APAc", file = paste0(out_dir, "runstat.o"), append = TRUE)
  RE <- fread(RE_prebed, header = F, drop = 4)
  SA <- fread(SA_prebed, header = F, drop = 4)
  SA <- SA[, V3 := V3 * SA_weight]
  APAc <- fread(APAc_file, header = FALSE)
  input_strand <- fread(SA_prebed, header = F, nrows = 1)$V4
  # input_strand <- "-" ### test


  ### chr loop
  alist <- list()
  for (chrno in Reduce(intersect,list(SA$V1, RE$V1, APAc$V1))) {

    ## by chr
    print(paste0("map the prebed and APAc in ", chrno, " to define APA count region"))
    cat(paste0("map the prebed and APAc in ", chrno, " to define APA count region"),
        file = paste0(out_dir, "runstat.o"), append = TRUE)
    SA_chr <- SA[V1 == chrno]
    RE_chr <- RE[V1 == chrno]
    APAchr <- APAc[V1 == chrno]
    setnames(APAchr, c("V1", "V2", "V3"), c("chr", "start", "end"))
    setorder(RE_chr, V2)
    # nodes with high observations

    ## remove isolated nodes
    # find gap
    if (nrow(RE_chr) >= 3){
      gap <- RE_chr[2:.N, V2] - RE_chr[1:(.N-1), V2]
      gaploc <- which(gap > blocksize)
      # filter isolated nodes
      gaploc_on <- gaploc[which((gaploc %in% (gaploc + 1)))]
      if (length(gaploc_on != 0)) {
        RE_chr <- RE_chr[-gaploc_on]
      }
      # supplementary adjustment
      setorder(RE_chr, V2)
      # check head/tail gap
      if (nrow(RE_chr) >= 3){
        if (RE_chr[2, V2] - RE_chr[1, V2] > (blocksize)) {
          RE_chr <- RE_chr[-1]
        }
        if (RE_chr[.N, V2] - RE_chr[.N-1, V2] > (blocksize)) {
          RE_chr <- RE_chr[-(.N)]
        }
      }

      ## clean environment
      rm(gap)
      rm(gaploc)
      rm(gaploc_on)
      gc()

      ## bind SA and RC
      RE_chr <- RE_chr[, 1:3]
      CLH_chr <- rbindlist(list(RE_chr, SA_chr))
      setkey(CLH_chr, V2)
      setorder(CLH_chr, V2)

    }else if (nrow(RE_chr) == 2) { #  == 2
      if (RE_chr[2, V2] - RE_chr[1, V2] < blocksize) {
        ## bind SA and RC
        RE_chr <- RE_chr[, 1:3]
        CLH_chr <- rbindlist(list(RE_chr, SA_chr))
        setkey(CLH_chr, V2)
        setorder(CLH_chr, V2)
      }else{
        ## SA only
        CLH_chr <- SA_chr[1:.N]
        setkey(CLH_chr, V2)
        setorder(CLH_chr, V2)
      }
    }else{ # == 1
      CLH_chr <- SA_chr[1:.N]
      setkey(CLH_chr, V2)
      setorder(CLH_chr, V2)
    }

    ## clean
    rm(SA_chr)
    rm(RE_chr)

    ## check strand
    if (input_strand == "-") {
      maxposiaddone <- CLH_chr[, max(V2)] + 1
      CLH_chr[, V2 := maxposiaddone - V2]
      setorder(CLH_chr, V2)
      APAchr[, start := maxposiaddone - start]
      APAchr[, end := maxposiaddone - end]
      setnames(APAchr, c("start", "end"), c("end", "start"))
      setcolorder(APAchr, c("start", "end"))
      setorder(APAchr, start)
    }

    ## calculate count
    CLH_chr <- CLH_chr[, sum(V3), by = V2]
    setnames(CLH_chr, c("V2", "V1"), c("loc", "count"))

    if (nrow(CLH_chr) >= 2){
      ## merge the nodes by gap distance
      gap <- CLH_chr[2:.N, loc] - CLH_chr[1:(.N-1), loc]
      gaploc <- which(gap > blocksize)
      # add head and tail
      endloc <- c(gaploc, nrow(CLH_chr))
      startloc <- c(1, gaploc + 1)

      ## define rough APAc
      REc <- data.table(start_index = startloc, end_index = endloc)

      ## clean
      rm(startloc)
      rm(endloc)

      ## separate the rough clusters
      print(paste0("cluster the nodes of ",chrno," to get APA clusters"))
      c_start <- c()
      c_end <- c()
      for (i in 1 : nrow(REc)) {
        # load the nodes of the rough clusters
        clu_dt <- CLH_chr[REc[i, start_index]: REc[i, end_index]]
        # detect dump region
        gap <- clu_dt[2:.N, loc] - clu_dt[1:(.N-1), loc]
        gaploc <- which(gap > APAdist_cut)
        blockloc <- which(gap > blocksize)
        RE_dist <- clu_dt[.N, loc] - clu_dt[1, loc] - sum(gap[gaploc])
        # too small to separate
        if (length(gaploc)  == 0 | (RE_dist < APAdist_cut)) {
          c_start <- c(c_start, clu_dt[1, loc])
          c_end <- c(c_end, clu_dt[.N, loc])
          # c_start <- c(c_start, 0)
          # end <- c(end, clu_dt[, sum(count)])
        } else {
          # eval the range of k, to make long block separated and long dump ignored
          k_l <- max(length(blockloc) + 1, 1)
          k_r <- max(k_l+1, length(blockloc) +  max(ceiling(RE_dist/blocksize), 1))
          k_range <- c(k_l, k_r)
          # cluster
          ckmean_1 <- Ckmeans.1d.dp::Ckmeans.1d.dp(as.numeric(na.omit(clu_dt$loc)),
                                                   k = k_range, as.numeric(na.omit(clu_dt$count)))

          # get the start & end of APA clusters
          clu_dt[, group := ckmean_1$cluster]
          c_start <- c(c_start, clu_dt[, min(loc), by = group]$V1)
          c_end <- c(c_end, clu_dt[, max(loc), by = group]$V1)
          # start <- c(start, length(unique(ckmean_1$cluster)))
          # end <- c(end, clu_dt[, sum(count), by = group]$V1)

          # clean
          rm(k_l)
          rm(k_r)
          rm(k_range)
          rm(ckmean_1)
        }

        ## clean
        rm(gap)
        rm(gaploc)
        rm(RE_dist)
        rm(blockloc)
        rm(clu_dt)
      }
      # clean
      rm(REc)

    }else{
      c_start <- c(CLH_chr$loc)
      c_end <- c(CLH_chr$loc)
    }

    ## clean
    rm(CLH_chr)


    # fwrite(REc[, .(chr = chrno, startloc = CLH_chr[start_index, loc], endloc = CLH_chr[end_index, loc])],
    #        "bed/pbmcroughtest.bed", quote = FALSE, sep = "\t", col.names = FALSE, scipen=9999)
    # REc[, .(chr = chrno, startloc = CLH_chr[start_index, loc], endloc = CLH_chr[end_index, loc], .I)][startloc == 14038335]


    ## get the data.table of APA clusters
    # chr_clu_bed <- data.table(Chr = chrno, clu_startposi = start, clu_endposi = end)
    # fwrite(chr_clu_bed, "bed/pbmctest.bed", quote = FALSE, sep = "\t", col.names = FALSE, scipen=9999)

    ## map clusters to APAc
    clu_todo <- data.table(clu_startposi = c_start, clu_endposi = c_end)
    apa_todo <- APAchr[, .(start, end)]
    setkey(clu_todo, clu_startposi,clu_endposi)
    setkey(apa_todo, start, end)
    # findoverlaps & map
    mapindex <- foverlaps(apa_todo, clu_todo, type = "any", which=TRUE, nomatch = NA)
    chr_map <- data.table(apa_todo[mapindex$xid], clu_todo[mapindex$yid])

    ## resize
    # hold
    chr_map[start < clu_startposi, clu_startposi := start]
    chr_map[end > clu_endposi, clu_endposi := end]
    chr_map[is.na(clu_startposi), clu_startposi := start]
    chr_map[is.na(clu_endposi), clu_endposi := end]
    # cbind
    chr_map <- chr_map[, .(clu_l = min(clu_startposi), clu_r = max(clu_endposi)),
                       by = .(start, end)]
    # cut end
    chr_map[, clu_r := end]

    if (nrow(chr_map) >= 2) {
      ## overlap
      gap_keep <- which(chr_map[2:.N, clu_l] - chr_map[1:(.N-1), clu_r] <= 0)
      # gap_ol <- gap_keep + 1
      # chr_map[gap_keep + 1, clu_l] - chr_map[gap_keep, clu_r]
      chr_map[gap_keep + 1, clu_l := chr_map[gap_keep, clu_r] + 1]
      # clean
      rm(gap_keep)
    } # == 1, pass

    ## check strand
    if (input_strand == "-") {
      chr_map[, start := maxposiaddone - start]
      chr_map[, end := maxposiaddone - end]
      chr_map[, clu_l := maxposiaddone - clu_l]
      chr_map[, clu_r := maxposiaddone - clu_r]
      setnames(chr_map, c("start", "end", "clu_l", "clu_r"),
               c("end", "start", "clu_r", "clu_l"))
      setcolorder(chr_map, c("start", "end", "clu_l", "clu_r"))
      setorder(chr_map, start)
    }

    ## write list
    alist[[chrno]] <- chr_map[, .(chr = chrno, start, end, clu_l, clu_r)]

    ## clean
    rm(APAchr)
    rm(chr_map)
    rm(mapindex)
    rm(apa_todo)
    rm(clu_todo)
    rm(c_start)
    rm(c_end)

  }


  ### bind
  APAc <- rbindlist(alist)[complete.cases(start, clu_l)]
  fwrite(APAc, paste(out_dir, "APA_map.out", sep = ""),
         quote = FALSE, sep = "\t", col.names = FALSE)
  print("MapAPAc: finish")
  cat("MapAPAc: finish", file = paste0(out_dir, "runstat.o"), append = TRUE)

}
