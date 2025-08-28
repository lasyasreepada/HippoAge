#Functions to perform CpG - gene linkage analysis 

getGenes <- function(platform) {
  genome <- sesameData_check_genome(NULL,platform=platform)
  sesameData_getTxnGRanges(genome, merge2gene = TRUE)
}

getProbeGR <- function(qry,platform,gr=NULL,genome="hg38") {
  if(is.null(gr)) {
    gr <- sesameData_getManifestGRanges(
      platform=platform, 
      genome = genome
    )
  }
  qry2 <- qry[qry %in% names(gr)]
  if(!identical(qry2,qry)) {
    missing <- length(qry) - length(qry2)
    message(paste("Lost ",as.character(missing), "probes from query"))
  }
  probe_gr <- gr[qry2,]
}

overlapGenes <- function(gene_gr,probe_gr,distance=10000) {
  hits <- findOverlaps(
    gene_gr,
    probe_gr + distance,
    ignore.strand = TRUE
  )
  id_df <- data.frame(
    Probe_ID=names(probe_gr)[subjectHits(hits)],
    ENST=names(gene_gr)[queryHits(hits)],
    Gene=gene_gr$gene_name[queryHits(hits)]
  )
  probe_coord_df <-  data.frame(
    chr_probe = seqnames(probe_gr)[subjectHits(hits)],
    start_probe = start(probe_gr)[subjectHits(hits)],
    end_probe = end(probe_gr)[subjectHits(hits)]
  )
  gene_coord_df <-  data.frame(
    chr_gene = seqnames(gene_gr)[queryHits(hits)],
    start_gene = start(gene_gr)[queryHits(hits)],
    end_gene = end(gene_gr)[queryHits(hits)]
  )
  cbind(id_df,probe_coord_df,gene_coord_df)
}

linkGenes <- function(qry,platform,gene_gr=NULL,distance) {
  if(is.character(qry)) {
    probe_gr <- getProbeGR(qry=qry,platform=platform)
  } else {
    probe_gr <- qry
  }
  if(is.null(gene_gr)) {
    gene_gr <- getGenes(platform=platform)
  }
  overlapGenes(
    gene_gr=gene_gr,
    probe_gr = probe_gr,
    distance = distance
  )
}