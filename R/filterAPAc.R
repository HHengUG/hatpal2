#' Filter internal APAc
#'
#' @param APAc_file APA clusters file from ClusterAPA()
#' @param fasta_file fasta file with an index
#' @param input_strand strand
#' @param out_dir The directory where the output files are stored
#' @param dist Scanning range (default c(20, 50))
#' @param AT_rich A rich motif (default "AAAAAAA")
#'
#' @return NULL, writes the output "APA_clusters.filtered.out" in the out_dir.
#' @export
#'
#' @import data.table Biostrings
#' @rawNamespace import(GenomicAlignments, except = c(first, second, last))
#'
#' @examples FilterAPAc("APA_clusters.out", "hg19.fa", "-", "output/")
FilterAPAc <- function(APAc_file, fasta_file, input_strand, out_dir, dist = c(20, 50), AT_rich = "AAAAAAA"){


  ### load fasta and APAc
  print("FilterAPAc: load fasta and APAc")
  cat("FilterAPAc: load fasta and APAc", file = paste0(out_dir, "runstat.o"), append = TRUE)
  APAc <- fread(APAc_file, header = FALSE)
  fa_seq <- Biostrings::readDNAStringSet(fasta_file)
  if (input_strand == "-") {
    switch <- strsplit(chartr("ATGC","TACG",AT_rich), split = "")[[1]]
    AT_rich <- paste(rev(switch), collapse = "")
    rm(switch)
  }

  ### filter APAc by A/T rich region
  ## chr loop
  print(paste0("There are ", length(intersect(APAc$V1, names(fa_seq))), " chrs found both in APAc and Fasta"))
  cat(paste0("There are ", length(intersect(APAc$V1, names(fa_seq))), " chrs found both in APAc and Fasta"),
      file = paste0(out_dir, "runstat.o"), append = TRUE)
  if (length(intersect(APAc$V1, names(fa_seq))) == 0) {
    warning("Fasta file does not match the APAc file, check the chr names")
    return(message("ERROR: check the chr names of APAc and Fasta"))
  }
  alist <- list()
  for (chrno in intersect(APAc$V1, names(fa_seq))) {

    ## APAchr by chr
    print(paste0("Filter APAc in ", chrno))
    cat(paste0("Filter APAc in ", chrno), file = paste0(out_dir, "runstat.o"), append = TRUE)
    APAchr <- APAc[V1 == chrno]

    ## reduce negative
    # APAV2 <- c(APAchr$V2)
    # APAV3 <- c(APAchr$V3)
    # APAchr[V2 <= dist[2], V2 := dist[2] + 0]
    # APAchr[V3 <= dist[1], V3 := dist[1] + 0]
    fa_len <- length(fa_seq[[chrno]])

    ## extract sequences by APAchr
    if (input_strand == "+"){
      APAc_seq <- lapply(APAchr$V3, function(x)
        as.character(fa_seq[[chrno]][max((x-dist[1]), 0):min((x+dist[2]), fa_len)]))
    }else if (input_strand == "-") {
      APAc_seq <- lapply(APAchr$V2, function(x)
        as.character(fa_seq[[chrno]][max((x-dist[2]), 0):min((x+dist[1]), fa_len)]))
    }


    ## identify and filter A/T rich region
    richindex <- grepl(AT_rich, APAc_seq)
    richloc <- which(richindex == TRUE)

    ## write list
    if (length(richloc) != 0) {
      alist[[chrno]] <- APAchr[-richloc]
    }else{
      alist[[chrno]] <- APAchr[1:.N]
    }

    ## clean
    rm(APAchr)
    rm(APAc_seq)
    rm(richindex)
    rm(richloc)
  }


  ### bind
  APAc <- rbindlist(alist)[complete.cases(V2, V3)]
  fwrite(APAc, paste(out_dir, "APA_clusters.filtered.out", sep = ""),
         quote = FALSE, sep = "\t", col.names = FALSE)
  print("FilterAPAc: finish")
  cat("FilterAPAc: finish", file = paste0(out_dir, "runstat.o"), append = TRUE)

}
