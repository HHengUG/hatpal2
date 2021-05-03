#' Cluster the two-dimensional coordinate files (SA and RC/SA) to define the APA clusters
#'
#' @param SA_prebed SA file of one strand generated from ExtractAPA()
#' @param uni_prebed rc_uni or SA file of one strand generated from ExtractAPA()
#' @param out_dir The directory where the output files are stored
#' @param gapdist_cut Threshold for filtering out isolated nodes (default 5)
#' @param APAdist_cut Threshold for the gap between APA (default 25)
#' @param RC_cutoff High observation nodes (default 3)
#' @param SA_weight The weight of SA file (default 3)
#'
#' @return NULL, writes the output "APA_clusters.out" in the out_dir.
#' @export
#'
#' @import data.table
#' @rawNamespace import(GenomicAlignments, except = c(first, second, last))
#'
#' @examples ClusterAPA("sa.prebed", "rc_uni.prebed", "output/")
ClusterAPA <- function(SA_prebed, uni_prebed, out_dir, gapdist_cut = 5, APAdist_cut = 25, RC_cutoff = 3, SA_weight = 3){


  ### load prebeds
  print("ClusterAPA: load prebeds")
  cat("ClusterAPA: load prebeds", file = paste0(out_dir, "runstat.o"), append = TRUE)
  uni <- fread(uni_prebed, header = FALSE)
  SA <- fread(SA_prebed, header = F, drop = 4)
  # SA weight
  SA <- SA[, V3 := V3 * SA_weight]
  # define strand
  input_strand <- fread(SA_prebed, header = F, nrows = 1)$V4


  ### chr loop
  alist <- list()
  for (chrno in union(uni$V1, SA$V1)) {

    ## load
    # chrno <- "chr1"
    print(paste0("load and process ", chrno, " to define APA clusters"))
    cat(paste0("load and process ", chrno, " to define APA clusters"), file = paste0(out_dir, "runstat.o"), append = TRUE)
    SA_chr <- SA[V1 == chrno, .(V2)]
    uni_chr <- uni[V1 == chrno, .(V2, V3)]
    setorder(uni_chr, V2)
    # nodes with high observations
    uni_chr_rc <- uni_chr[V3 >= RC_cutoff] # RC >= 3

    ## remove isolated nodes
    if (nrow(uni_chr) >= 3) {
      # find gap
      gap <- uni_chr[2:.N, V2] - uni_chr[1:(.N-1), V2]
      gaploc <- which(gap > gapdist_cut)
      # filter isolated nodes
      gaploc_on <- gaploc[which((gaploc %in% (gaploc + 1)))]
      if (length(gaploc_on != 0)) {
        uni_chr <- uni_chr[-gaploc_on]
      }

      if (nrow(uni_chr) >= 3) {
        # check head/tail gap
        if ((uni_chr[2, V2] - uni_chr[1, V2] > gapdist_cut) & (uni_chr[1, V3] < RC_cutoff)) {
          uni_chr <- uni_chr[-1]
        }
        if ((uni_chr[.N, V2] - uni_chr[.N-1, V2] > gapdist_cut) & (uni_chr[.N, V3] < RC_cutoff)) {
          uni_chr <- uni_chr[-(.N)]
        }
      }

      ## supplementary adjustment
      uni_chr <- unique(rbind(uni_chr, uni_chr_rc))
      setorder(uni_chr, V2)

      ## clean
      rm(gap)
      rm(gaploc)
      rm(gaploc_on)

      ## nrow
    }else if (nrow(uni_chr) == 2) { #  == 2
      if (uni_chr[2, V2] - uni_chr[1, V2] > gapdist_cut) {
        uni_chr <- uni_chr[V3 >= RC_cutoff]
      }
    }else{ # == 1
      uni_chr <- uni_chr[V3 >= RC_cutoff]
    }

    ## bind SA and RC
    # RC_chr <- uni_chr[, 1:3]
    CLH_chr <- rbindlist(list(uni_chr[, .(V2)], SA_chr))
    setkey(CLH_chr, V2)
    setorder(CLH_chr, V2)
    # calculate count
    # CLH_chr <- CLH_chr[, sum(V3), by = V2]
    # setnames(CLH_chr, c("V2", "V1"), c("loc", "count"))
    CLH_chr <- CLH_chr[!duplicated(V2)]
    setnames(CLH_chr, c("V2"), c("loc"))

    if (nrow(CLH_chr) >= 2) {
      ## merge the nodes by gap distance
      gap <- CLH_chr[2:.N, loc] - CLH_chr[1:(.N-1), loc]
      gaploc <- which(gap > APAdist_cut)
      # add head and tail
      endloc <- c(gaploc, nrow(CLH_chr))
      startloc <- c(1, gaploc + 1)

      ## define rough APAchr
      APAchr <- data.table(chr = chrno, start = CLH_chr[startloc, loc], end = CLH_chr[endloc, loc])

      ## clean
      rm(gap)
      rm(gaploc)
      rm(startloc)
      rm(endloc)

    }else{
      APAchr <- data.table(chr = chrno, start = CLH_chr[1, loc], end = CLH_chr[1, loc])
    }

    ## list
    alist[[chrno]] <- APAchr

    ## clean environment
    rm(uni_chr_rc)
    rm(SA_chr)
    rm(uni_chr)
    rm(CLH_chr)
    rm(APAchr)
  }

  ### bind
  APAc <- rbindlist(alist)[complete.cases(start, end)]
  fwrite(APAc, paste(out_dir, "APA_clusters.out", sep = ""),
         quote = FALSE, sep = "\t", col.names = FALSE)
  print("ClusterAPA: finish")
  cat("ClusterAPA: finish", file = paste0(out_dir, "runstat.o"), append = TRUE)
}
