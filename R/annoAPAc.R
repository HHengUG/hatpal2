#' Annotate the APA clusters from a GTF file
#'
#' @param APAc_file  APA clusters file from MapAPAc()
#' @param gtf_file A reference GTF file
#' @param w_strand which strand, "+" or "-"
#' @param add_chr whether to add character "chr" on the chr_name of the GTF (default FALSE)
#'
#' @return Null, An annotated file in the same dir.
#' @export
#'
#' @import data.table rtracklayer
#' @rawNamespace import(GenomicAlignments, except = c(first, second, last))
#'
#'
#' @examples AnnoAPAc("APA_map.out", "hg19.gtf", "-")
AnnoAPAc <- function(APAc_file, gtf_file, w_strand, add_chr = FALSE){


  ### load files
  print("AnnoAPAc: load gtf and APAc.map")
  cat("AnnoAPAc: load gtf and APAc.map", file = paste0(APAc_file, ".runstat.o"), append = TRUE)
  gtf <- rtracklayer::import(gtf_file)
  gtf <- as.data.table(gtf)
  APAc <- fread(APAc_file, header = F)
  anno_order <- c("three_prime_utr", "stop_codon", "UTR", "exon", "five_prime_utr",
                  "gene", "transcript", "protein", "CDS", "start_codon")

  ### filter gtf
  ## rename gtf and APAc
  gtf <- gtf[, .(seqnames, type, start, end, strand, gene_name)]
  setnames(gtf, 1:6, c("Chr","GType","Start","End","Strand","Gene"))
  setnames(APAc, c("Chr", "clu_startposi", "clu_endposi", "map_start", "map_end"))
  setkey(gtf, Chr, Start, End)
  setkey(APAc, Chr, clu_startposi, clu_endposi)

  ## strand and chr
  gtf <- gtf[Strand == w_strand]
  # add string "chr" to gtf or not
  if (add_chr == TRUE){gtf[, Chr := paste0("chr", Chr)]}
  # standard ref
  if (w_strand == "+") {
    gtf$GType <- factor(gtf$GType, levels = anno_order)
    setorder(gtf, Chr, Gene, GType, -End)
    gtf <- gtf[, .SD[1], by = .(Chr, Gene, GType)]
  }else if(w_strand == "-"){
    gtf$GType <- factor(gtf$GType, levels = anno_order)
    setorder(gtf, Chr, Gene, GType, Start)
    gtf <- gtf[, .SD[1], by = .(Chr, Gene, GType)]
  }
  setorder(gtf, Chr, Start, End)


  ### another gtf to anno the dist clusters
  naimpu <- gtf[GType %in% c("three_prime_utr", "UTR", "exon"), .(Chr, Start, End, Gene, GType)]
  setkey(naimpu, Chr, Start, End)


  ### chr loop for anno
  all_clu_anno <- list()
  for (chrno in intersect(gtf$Chr, APAc$Chr)) {

    ## define two extra data.table to match locations
    print(paste0("Annotate the APAc.map in ", chrno))
    cat(paste0("Annotate the APAc.map in ", chrno), file = paste0(APAc_file, ".runstat.o"), append = TRUE)
    gtf_todo <- gtf[Chr == chrno, .(Start,End)]
    bed_todo <- APAc[Chr == chrno, .(clu_startposi, clu_endposi)]
    setkey(gtf_todo, Start, End)
    setkey(bed_todo, clu_startposi, clu_endposi)

    ## findoverlaps & map
    mapindex <- foverlaps(bed_todo, gtf_todo, which=TRUE, nomatch = NA)
    chr_anno <- data.table(Chr=chrno,
                           APAc[Chr == chrno, ][mapindex[["xid"]],
                                                .(clu_startposi, clu_endposi, map_start, map_end)],
                           gtf[Chr == chrno][mapindex[["yid"]], .(Gene, GType, Start, End)])
    # clean
    rm(gtf_todo)
    rm(bed_todo)
    rm(mapindex)

    ## set the order of annotation
    # relative distance
    if (w_strand == "+") {
      chr_anno[, loc := (clu_endposi - Start) / (End - Start)]
    }else if(w_strand == "-"){
      chr_anno[, loc := (clu_startposi - End) / (Start - End)]
    }
    chr_anno[ ,c("Start", "End") := NULL]
    # ordered
    chr_anno$GType <- factor(chr_anno$GType, levels = anno_order)
    setorder(chr_anno, clu_startposi, GType, -loc)

    ## get the first anno of all annotations
    chr_anno <- chr_anno[, .SD[1], by = clu_startposi]
    setkey(chr_anno, Chr, clu_startposi, clu_endposi)
    setcolorder(chr_anno, c("Chr", "clu_startposi", "clu_endposi", "Gene", "GType"))

    ## anno the na cluster with the closest anno ahead
    nagtf <- naimpu[Chr == chrno]
    if (nrow(nagtf) != 0){
      nagtf$GType <- factor(nagtf$GType,levels = c("three_prime_utr", "UTR", "exon"))
      if (w_strand == "+") {
        # positive strand
        nagtf <- nagtf[, .SD[1], by = End]
        setkey(nagtf, End, GType)
        setorder(nagtf, End, GType)
        # find the closest anno before the cluster
        remap <- chr_anno[is.na(Gene),
                          .(sapply(clu_startposi,
                                   function(x) max(which((x - nagtf$End) >= 0), 1)))]
        dist <- (chr_anno[is.na(Gene), clu_startposi] - nagtf[remap$V1, End]) /
          nagtf[remap$V1, End - Start]
      } else if (w_strand == "-") {
        # positive strand
        nagtf <- nagtf[, .SD[1],by = Start]
        setkey(nagtf, Start, GType)
        setorder(nagtf, Start, GType)
        # find the closest anno after the cluster
        remap <- chr_anno[is.na(Gene),
                          .(sapply(clu_endposi,
                                   function(x) min(which((x - nagtf$Start) <= 0), nrow(nagtf))))]
        dist <- (nagtf[remap$V1, Start] - chr_anno[is.na(Gene), clu_endposi]) /
          nagtf[remap$V1, End - Start]
      }

      ## impute the clusters without anno with a closest anno in another col
      chr_anno[is.na(Gene),
               c("NRgene", "NRGtype", "Dist") :=
                 .(nagtf[remap$V1, Gene], nagtf[remap$V1, GType], dist)]

      ## clean
      rm(nagtf)
      rm(remap)
      rm(dist)
    }

    ## get the data.table of annotated APA clusters
    all_clu_anno[[chrno]] <- chr_anno

    ## clean
    rm(chr_anno)
  }


  ### write down the merged file
  all_clu_anno_DT <- rbindlist(all_clu_anno)
  all_clu_anno_DT[, Dist := round(Dist, 2)]
  all_clu_anno_DT[, loc := round(loc, 2)]
  fwrite(all_clu_anno_DT, paste0(APAc_file,".anno"),
         quote = FALSE, sep = "\t", na = "NA", col.names = FALSE)
  print("AnnoAPAc: finish")
  cat("AnnoAPAc: finish", file = paste0(APAc_file, ".runstat.o"), append = TRUE)
}
