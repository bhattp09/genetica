setwd("/data/Priya/genetica/")
options(stringsAsFactors = F)

## Dependent Libraries
library(plyr)
library(stringr)
library(tidyr)
library(xlsx)
library(biomaRt)
library(limma)
library(goseq)
library(GO.db)
library(openxlsx)


getMQLresults_short <- function(myInterpretList, results_fname_xlsx){
  
  #extract name data
  extract.name <- lapply(myInterpretList, function(y) y[["name"]])
  extract.name <- do.call(rbind, extract.name)
  
  #extract location data
  extract.location <- lapply(myInterpretList, function(y) y[["location"]])
  extract.location <- do.call(rbind, extract.location)
  names(extract.location)[names(extract.location) == "MAPINFO"] <- "BP"
  
  final_name_location_df <- merge(extract.location, extract.name[c(2:6)], by.x = "row.names", by.y = "row.names")
  rownames(final_name_location_df) <- final_name_location_df$Row.names
  final_name_location_df$Row.names <- NULL
  colnames(final_name_location_df)[1:3] <- c("probe_id", "probe_chr", "probe_pos")
  
  #############
  # mQTL Membership
  #############
  
  extract.mQTL <- lapply(myInterpretList, function(y) y[["mqtls"]])
  extract.mQTL <- do.call(rbind, extract.mQTL)
  names(extract.mQTL)[names(extract.mQTL) == "File"] <- "mQTL_dataset"
  
  if(dim(extract.mQTL)[1]!= 0){
    extract.mQTL <- extract.mQTL[c(1,6:8,2:5,9:10)]
    extract.mQTL <- extract.mQTL[order(extract.mQTL$mQTL_dataset, extract.mQTL$DNA_M_PROBEID, extract.mQTL$SNP_CHR, extract.mQTL$SNP_BP),]
    rownames(extract.mQTL) <- NULL
  }
  master <- data.frame(final_name_location_df)
  master$has_mQTLs <- master$probe_id %in% extract.mQTL$DNA_M_PROBEID
  
  results <- list(master, extract.mQTL)
  names(results) <- c("Summary", "mQTLs")
  
  for (i in 1:length(results)) {
    
    if (nrow(results[[i]]) == 0){
      
      results[[i]] <- paste("There were no", names(results)[i], "results for your query", sep = " ")
      
    }
    
  }
  
  ########
  # Write in one excel file
  ########
  
  openxlsx::write.xlsx(results, file = results_fname_xlsx, keepNA = TRUE)
  
  
  return(results)
}

geneticaList <- function(queryList, input_genetic_data_type = c("gene", "snp", "probe")) {
  
  ## input genetic type = gene
  if (input_genetic_data_type == "gene") {
    
    interpret_headers <- c("name", "location", "goCategories", "gtex_cisEqtls", "gwasCatalog", "pgc_top_hits")
    toInterpret <- rep(list(NA),length(interpret_headers))
    names(toInterpret) = interpret_headers
    
    queryList <- sapply(queryList, function(x) {alias2Symbol(x, species = "Hs")}, USE.NAMES = FALSE)
    queryList <- queryList[lapply(queryList,length)>0]
    num_queryList <- length(queryList)
    
    myInterpretList = rep(list(toInterpret), num_queryList)
    names(myInterpretList) = queryList
    
    
    genemart = biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                                host = "grch37.ensembl.org", path = "/biomart/martservice", 
                                dataset = "hsapiens_gene_ensembl")
    
    genemart38 <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
    
    # bm_genemart <- biomaRt::getBM(attributes = c("external_gene_name", 
    #                                              "ensembl_gene_id", "description"), filters = "external_gene_name", 
    #                               values = names(myInterpretList), mart = genemart)
    # 
    # bm_locmart <- biomaRt::getBM(attributes = c("external_gene_name","chromosome_name", "start_position", "end_position", "strand"), filters = "external_gene_name", 
    #                              values = names(myInterpretList), mart = genemart)
    
    bm_genemart <- biomaRt::getBM(attributes = c("external_gene_name", "ensembl_gene_id", "description","chromosome_name", "start_position", "end_position", "strand"), 
                                  filters = "external_gene_name", 
                                  values = names(myInterpretList), 
                                  mart = genemart)
    
    bm_genemart$chromosome_name <- as.character(bm_genemart$chromosome_name)
    chr_check <- c(1:22,"X","Y")
    bm_genemart <- bm_genemart[bm_genemart$chromosome_name %in% chr_check,]
    
    
    bm_genemart38 <- biomaRt::getBM(attributes = c("external_gene_name", "chromosome_name", "start", "end"), 
                                    filters = "external_gene_name", 
                                    values = names(myInterpretList), 
                                    mart = genemart38)
    
    bm_genemart38$chromosome_name <- as.character(bm_genemart38$chromosome_name)
    bm_genemart38 <- bm_genemart38[bm_genemart38$chromosome_name %in% chr_check,]
    colnames(bm_genemart38) = c("external_gene_name", "chromosome_name_hg38", "start_hg38", "end_hg38")
    
    bm_genemart <- merge(bm_genemart, bm_genemart38, by = "external_gene_name", all = TRUE)
    no_interpretation_available <- list()
    for (g in queryList) {
      
      if (g %in% bm_genemart[, 1]) {
        
        bm_g <- bm_genemart[which(bm_genemart[, 1] == g), ]
        
        myInterpretList[[g]]$name <- bm_g[c("external_gene_name", "ensembl_gene_id", "description")]
        myInterpretList[[g]]$location <- bm_g[c("external_gene_name","chromosome_name", "start_position", "end_position", "strand", "chromosome_name_hg38", "start_hg38", "end_hg38")]
      }
      
      else {
        
        no_interpretation_available[[g]] <- g
        myInterpretList[[g]]$name <- data.frame(external_gene_name = g, ensembl_gene_id = NA, description  = NA)
        myInterpretList[[g]]$location <- data.frame(external_gene_name = g, chromosome_name = NA, start_position  = NA, end_position = NA, strand  = NA, chromosome_name_hg38 = NA, start_hg38 = NA, end_hg38 = NA)
      }
      
    }
    
    
    if (length(no_interpretation_available) != 0) {
      
      warning("Gene names missing: ", paste(no_interpretation_available, 
                                            collapse = ", "))
    }
  }
  
  ## input genetic type = snp
  else if (input_genetic_data_type == "snp") {
    
    interpret_headers <- c("name", "location", "goCategories", "gtex_cisEqtls", "gwasCatalog", "pgc_top_hits", "mqtls")
    toInterpret <- rep(list(NA),length(interpret_headers))
    names(toInterpret) = interpret_headers
    
    
    num_queryList <- length(queryList)
    
    myInterpretList = rep(list(toInterpret), num_queryList)
    names(myInterpretList) = queryList
    
    
    snpmart = biomaRt::useMart(biomart = "ENSEMBL_MART_SNP", 
                               host = "grch37.ensembl.org", path = "/biomart/martservice", 
                               dataset = "hsapiens_snp")
    
    genemart = biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                                host = "grch37.ensembl.org", path = "/biomart/martservice", 
                                dataset = "hsapiens_gene_ensembl")
    
    genemart38 <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
    
    bm_snpmart <- getBM(attributes = c("refsnp_id", "ensembl_gene_stable_id","allele"), 
                        filters = "snp_filter", 
                        values = names(myInterpretList), 
                        mart = snpmart)
    
    bm_genemart <- getBM(attributes = c("ensembl_gene_id", "external_gene_name","chromosome_name", "start_position", "end_position"), 
                         filters = "ensembl_gene_id", 
                         values = bm_snpmart$ensembl_gene_stable_id, 
                         mart = genemart)
    
    bm_genemart38 <- biomaRt::getBM(attributes = c("external_gene_name", "chromosome_name", "start", "end"), 
                                    filters = "external_gene_name", 
                                    values = names(myInterpretList), 
                                    mart = genemart38)
    
    bm_genemart38$chromosome_name <- as.character(bm_genemart38$chromosome_name)
    bm_genemart38 <- bm_genemart38[bm_genemart38$chromosome_name %in% chr_check,]
    colnames(bm_genemart38) = c("external_gene_name", "chromosome_name_hg38", "start_hg38", "end_hg38")
    
    bm_genemart <- merge(bm_genemart, bm_genemart38, by = "external_gene_name", all = TRUE)
    
    
    bm_snpmart <- merge(bm_snpmart, bm_genemart, by.x = 'ensembl_gene_stable_id', by.y = 'ensembl_gene_id', all.x = T)
    
    
    bm_snploc <- biomaRt::getBM(attributes = c("refsnp_id", "chr_name", "chrom_start", "chrom_end", "chrom_strand"), filters = "snp_filter", 
                                values = names(myInterpretList), mart = snpmart)
    names(bm_snploc)[names(bm_snploc) == "refsnp_id"] <- "snp_id"
    
    bm <- ddply(bm_snpmart, .(refsnp_id), summarize,
                ensembl_gene_stable_id = paste(unique(ensembl_gene_stable_id),collapse=";"),
                snp_id = paste(unique(refsnp_id),collapse=";"),
                allele = paste(unique(allele),collapse=";"),
                gene_name = paste(external_gene_name,collapse=";"),
                gene_chr37 = paste(unique(chromosome_name),collapse=";"),
                gene_start_position37 = paste(unique(start_position),collapse=";"),
                gene_end_position37 = paste(unique(end_position),collapse=";"),
                gene_chr38 = paste(unique(chromosome_name_hg38),collapse=";"),
                gene_start_position38 = paste(unique(start_hg38),collapse=";"),
                gene_end_position38 = paste(unique(end_hg38),collapse=";"))
    
    bm <- bm[c(2:4,1,5:10)]
    
    no_interpretation_available <- list()
    for (g in queryList) {
      
      if (g %in% bm[, c("snp_id")]) {
        
        myInterpretList[[g]]$name <- bm[which(bm[, c("snp_id")] == g), ]
        myInterpretList[[g]]$location <- bm_snploc[which(bm_snploc[, 1] == g), ]
      }
      
      else {
        
        no_interpretation_available[[g]] <- g
        myInterpretList[[g]]$name <- data.frame(snp_id = g, allele = NA, gene_name  = NA, ensembl_gene_stable_id = NA, gene_chr37  = NA, gene_start_position37 = NA, gene_end_position37  = NA, gene_chr38  = NA, gene_start_position38 = NA, gene_end_position38  = NA)
        myInterpretList[[g]]$location <- data.frame(snp_id = g, chr_name = NA, chrom_start  = NA, chrom_end = NA, chrom_strand  = NA)
      }
      
    }
    
    
    if (length(no_interpretation_available) != 0) {
      
      warning("SNP names missing: ", paste(no_interpretation_available, 
                                           collapse = ", "))
    }
    
    
  }
  
  ## input genetic type = probe
  else if (input_genetic_data_type == "probe") {
    interpret_headers <- c("name", "location", "illumina.manifest", "mqtls", "goCategories", "gtex_cisEqtls", "gwasCatalog", "pgc_top_hits")
    toInterpret <- rep(list(NA),length(interpret_headers))
    names(toInterpret) = interpret_headers
    
    toInterpret$gtex_cisEqtls <- list(gtex_mqtlSNPs = NA, gtex_mappedGenes = NA)
    toInterpret$gwasCatalog <- list(gwasCatalog_mqtlSNPs = NA, gwasCatalog_mappedGenes = NA)
    toInterpret$pgc_top_hits <- list(pgc_top_hits_mqtlSNPs = NA, pgc_top_hits_mappedGenes = NA) 
    
    num_queryList <- length(queryList)
    myInterpretList = rep(list(toInterpret), num_queryList)
    names(myInterpretList) = queryList
    
    illum_manifest <- read.delim("mqtl_files/MethylationEPIC_v-1-0_B2_v3.txt", header = T, stringsAsFactors = F)
    illum_manifest <- illum_manifest[which(illum_manifest$IlmnID %in% names(myInterpretList)),]
    
    illum_manifest_colnames <- names(illum_manifest)
    
    empty_illumina_manifest <- data.frame(matrix(ncol = length(illum_manifest_colnames), nrow = 1))
    colnames(empty_illumina_manifest) <- illum_manifest_colnames
    
    for (i in 1:length(illum_manifest$UCSC_RefGene_Name)){
      if (grepl(";", illum_manifest$UCSC_RefGene_Name[[i]])) {
        probe_genes <- sort(unique(strsplit(illum_manifest$UCSC_RefGene_Name,";| ]")[[i]]))
        
        illum_manifest$UCSC_RefGene_Name[[i]] <- paste(probe_genes, collapse = ";")
      }
    }
    
    genemart = biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                                host = "grch37.ensembl.org", path = "/biomart/martservice", 
                                dataset = "hsapiens_gene_ensembl")
    
    
    bm_genemart <- biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name"), 
                                  filters = "external_gene_name", 
                                  values = unlist(unique(strsplit(illum_manifest$UCSC_RefGene_Name,";| ]"))), 
                                  mart = genemart)
    
    genemart38 <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
    
    
    bm_genemart38 <- biomaRt::getBM(attributes = c("external_gene_name", "chromosome_name", "start", "end"),
                                  filters = "external_gene_name", 
                                  values = unlist(unique(strsplit(illum_manifest$UCSC_RefGene_Name,";| ]"))), 
                                  mart = genemart38)
    
    colnames(bm_genemart38) = c("external_gene_name", "chromosome_name_hg38", "start_hg38", "end_hg38")
    bm_genemart <- merge(bm_genemart, bm_genemart38, by = "external_gene_name", all = TRUE)
    
    bm_genemart_updated <- data.frame()
    for (i in 1:length(illum_manifest$UCSC_RefGene_Name)) {
      val = illum_manifest$UCSC_RefGene_Name[i]
      if (grepl(";", val)) {
        
        split_i <- strsplit(val,";| ]")[[1]]
        split_i <- sort(split_i)
        bm_genemart_split_i <- bm_genemart[bm_genemart$external_gene_name %in% split_i,]
        bm_genemart_split_i <- bm_genemart_split_i[order(bm_genemart_split_i$external_gene_name),]
        bm_genemart_split_i <- t(data.frame(apply(rbind( bm_genemart_split_i[c('ensembl_gene_id','external_gene_name')]), 2, paste, collapse=";")))
        bm_genemart_updated <- rbind(bm_genemart_updated,bm_genemart_split_i,make.row.names = FALSE)
      }
      
      else {
        bm_genemart_split_i = bm_genemart[bm_genemart$external_gene_name == val,]
        bm_genemart_updated <- rbind(bm_genemart_updated,bm_genemart_split_i,make.row.names = FALSE)
      }
      
    }
    
    illum_manifest <- merge(illum_manifest, bm_genemart_updated, by.x = "UCSC_RefGene_Name", by.y = "external_gene_name", all.x = T)
    
    no_interpretation_available <- list()
    for (g in queryList) {
      
      if (g %in% illum_manifest[, c("IlmnID")]) {
        
        name_df <- illum_manifest[which(illum_manifest[, "IlmnID"] == g), 
                       c("IlmnID", "Genome_Build", "UCSC_RefGene_Name", "UCSC_RefGene_Accession", "UCSC_RefGene_Group", "ensembl_gene_id")]
        
        name_df_collapsed <- plyr::ddply(name_df, .(IlmnID), summarize,
                    Genome_Build = paste(unique(Genome_Build),collapse=";"),
                    UCSC_RefGene_Name = paste(unique(UCSC_RefGene_Name),collapse=";"),
                    UCSC_RefGene_Accession = paste(unique(UCSC_RefGene_Accession),collapse=";"),
                    UCSC_RefGene_Group = paste(unique(UCSC_RefGene_Group),collapse=";"),
                    ensembl_gene_id = paste(unique(ensembl_gene_id),collapse=";"))
        
        myInterpretList[[g]]$name <- name_df_collapsed
        
        location_df <- illum_manifest[which(illum_manifest[, "IlmnID"] == g), 
                                  c("IlmnID", "CHR", "MAPINFO")]
        
        location_df_collapsed <- ddply(location_df, .(IlmnID), summarize,
                             CHR = paste(unique(CHR),collapse=";"),
                         POS = paste(unique(MAPINFO),collapse=";"))
        
        
        myInterpretList[[g]]$location <- location_df_collapsed
        
        myInterpretList[[g]]$illumina.manifest <- illum_manifest[which(illum_manifest[, "IlmnID"] == g),]
      }
      
      else {
        
        no_interpretation_available[[g]] <- g
        myInterpretList[[g]]$name <- data.frame(IlmnID = g, Genome_Build = NA, UCSC_RefGene_Name = NA, UCSC_RefGene_Accession  = NA, UCSC_RefGene_Group  = NA, ensembl_gene_id = NA)
        myInterpretList[[g]]$location <- data.frame(IlmnID = g, CHR  = NA, MAPINFO = NA)
        myInterpretList[[g]]$illumina.manifest <- empty_illumina_manifest
        myInterpretList[[g]]$illumina.manifest$IlmnID <- g
      }
      
    }
    
    
    
    if (length(no_interpretation_available) != 0) {
      
      warning("Probe names missing: ", paste(no_interpretation_available, 
                                             collapse = ", "))
    }
    
  }
  
  ## returns warning if input genetic type is != c("gene", "snp", "probe")
  else {
    
    warning(paste0("input genetic data type not supported:  ", input_genetic_data_type))
  }
  
  return(myInterpretList)
}
getGtexCisEqtls <- function(myInterpretList, data_dir = "./", input_genetic_data_type = c("gene", "snp", "probe")) {
  
  files <- list.files(paste(data_dir, "/GTEx_Analysis_v7_eQTL", 
                            sep = ""), pattern = "signif_variant_gene_pairs")
  
  gtexNames <- stringr::word(files, 1, sep = stringr::fixed("_Analysis"))
  gtexNames <- gsub(".v7.signif_variant_gene_pairs.txt.gz", "", gtexNames)
  
  varIDmap <- read.table(paste0(data_dir, "GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt"), header = T, stringsAsFactors = F)
  varIDmap <- varIDmap[c("chr", "variant_pos", "variant_id", "ref", "alt", "rs_id_dbSNP147_GRCh37p13")]
  names(varIDmap)[names(varIDmap) == "rs_id_dbSNP147_GRCh37p13"] <- "snp_id"
  
  
  if (input_genetic_data_type == "gene"){
    ensgs <- unlist(lapply(myInterpretList, function(x) {x$name$ensembl_gene_id} ))
    ensgs_clean <- ensgs[!is.na(ensgs)]
    ensgs_clean <- ensgs_clean[!duplicated(ensgs_clean)]
    ensgs_concat <- paste(ensgs_clean, collapse = "|")
    ensgs_concat <- paste("\'",ensgs_concat,"\'", sep="")
    
    zgrep_code <- unlist(lapply(files, function(x) paste0("zgrep -E ", ensgs_concat, 
                                                          " ", data_dir, "GTEx_Analysis_v7_eQTL/", x)))
    
    gtex_eqtls <- list()
    for (i in 1:length(zgrep_code)) {
      tryCatch({
        gtex_eqtls[[i]] <- system(zgrep_code[[i]], intern = TRUE)
        gtex_eqtls[[i]] <- as.data.frame(gtex_eqtls[[i]])
        gtex_eqtls[[i]] <- tidyr::separate(data = gtex_eqtls[[i]], col = 1, into = c("variant_id","gene_id",	"tss_distance",	
                                                                                     "ma_samples",	"ma_count",	"maf",	"pval_nominal",	"slope",	
                                                                                     "slope_se",	"pval_nominal_threshold",	"min_pval_nominal",	"pval_beta"), sep = "\t")
        gtex_eqtls[[i]] <- cbind(File = gtexNames[i], gtex_eqtls[[i]])
      }, error = function(i) {})
    }
    names(gtex_eqtls) <- gtexNames
    gtex_eqtls_df <- do.call("rbind", gtex_eqtls)
    gtex_eqtls_df <- merge(gtex_eqtls_df, varIDmap, by.x = "variant_id", by.y = "variant_id", all.x = T)
    gtex_eqtls_df <- gtex_eqtls_df[c(2,1,3:18)]
    
    
    gtex_eqtls_split <- list()
    for (i in ensgs_clean) {
      gtex_eqtls_split[[i]] <- gtex_eqtls_df[grep(i, gtex_eqtls_df$gene_id), ]
      gtex_eqtls_split[[i]] <- gtex_eqtls_split[[i]][order(gtex_eqtls_split[[i]]$File),]
    }
    
    names(gtex_eqtls_split) <- names(ensgs_clean)
    
    
  } ##  if (input_genetic_data_type == "gene") 
  
  else if (input_genetic_data_type == "snp") {
    ensgs <- names(myInterpretList)
    
    subset_varIDMap <- varIDmap[which(varIDmap$snp_id %in% ensgs),]
    subset_varIDMap <- subset_varIDMap[match(ensgs, subset_varIDMap$snp_id),]
    
    variantIDs_query <- subset_varIDMap$variant_id
    
    ensgs_clean <- variantIDs_query[!is.na(variantIDs_query)]
    ensgs_concat <- paste(ensgs_clean, collapse = "|")
    ensgs_concat <- paste("\'",ensgs_concat,"\'", sep="")
    
    zgrep_code <- unlist(lapply(files, function(x) paste0("zgrep -wE ", ensgs_concat, 
                                                          " ", data_dir, "GTEx_Analysis_v7_eQTL/", x)))
    
    
    gtex_eqtls <- list()
    for (i in 1:length(zgrep_code)) {
      tryCatch({
        gtex_eqtls[[i]] <- system(zgrep_code[[i]], intern = TRUE)
        gtex_eqtls[[i]] <- as.data.frame(gtex_eqtls[[i]])
        gtex_eqtls[[i]] <- tidyr::separate(data = gtex_eqtls[[i]], col = 1, into = c("variant_id","gene_id",	"tss_distance",	
                                                                                     "ma_samples",	"ma_count",	"maf",	"pval_nominal",	"slope",	
                                                                                     "slope_se",	"pval_nominal_threshold",	"min_pval_nominal",	"pval_beta"), sep = "\t")
        gtex_eqtls[[i]] <- cbind(File = gtexNames[i], gtex_eqtls[[i]])
      }, error = function(i) {})
    }
    names(gtex_eqtls) <- gtexNames
    gtex_eqtls_df <- do.call("rbind", gtex_eqtls)
    gtex_eqtls_df <- merge(gtex_eqtls_df, subset_varIDMap, by.x = "variant_id", by.y = "variant_id", all.x = T)
    gtex_eqtls_df <- gtex_eqtls_df[c(2,1,3:18)]
    
    genemart = biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                                host = "grch37.ensembl.org", path = "/biomart/martservice", 
                                dataset = "hsapiens_gene_ensembl")
    
    bm_genemart <- getBM(attributes = c("ensembl_gene_id_version", "external_gene_name","chromosome_name", "start_position", "end_position"), 
                         filters = "ensembl_gene_id_version", 
                         values = gtex_eqtls_df$gene_id, 
                         mart = genemart)
    
    bm_genemart  <- bm_genemart[match(gtex_eqtls_df$gene_id, bm_genemart$ensembl_gene_id_version),]
    names(bm_genemart) <- c("ensembl_gene_id_version", "expressed_gene_name", "expressed_gene_chr", "expressed_gene_start_pos", "expressed_gene_end_pos")
    
    gtex_eqtls_df <- merge(gtex_eqtls_df, bm_genemart, by.x = "gene_id", by.y = "ensembl_gene_id_version", all.x = T)
    gtex_eqtls_df <- gtex_eqtls_df[!duplicated(gtex_eqtls_df), ]
    
    gtex_eqtls_df <- gtex_eqtls_df[c(2,18,3,1,4:17,19:22)]
    
    
    gtex_eqtls_split <- list()
    for (i in ensgs) {
      gtex_eqtls_split[[i]] <- gtex_eqtls_df[grep(i, gtex_eqtls_df$snp_id), ]
      gtex_eqtls_split[[i]] <- gtex_eqtls_split[[i]][order(gtex_eqtls_split[[i]]$File),]
    }
    
    
  } ## else if (input_genetic_data_type == "snp") 
  
  else if (input_genetic_data_type == "probe") {
    
    genemart = biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                                host = "grch37.ensembl.org", path = "/biomart/martservice", 
                                dataset = "hsapiens_gene_ensembl")
    
    
    ### STEP 1: Mapped GENES
    ensgs <- unlist(lapply(myInterpretList, function(x) {x$name$ensembl_gene_id} ))
    
    ensgs_clean <- ensgs[!is.na(ensgs)]
    ensgs_clean <- ensgs_clean[ensgs_clean != ""]
    ensgs_clean <- ensgs_clean[ensgs_clean != "NA"]
    ensgs_concat <- paste(ensgs_clean, collapse = "|")
    ensgs_concat <- gsub(";", "|", ensgs_concat)
    ensgs_concat <- paste("\'",ensgs_concat,"\'", sep="")
    
    zgrep_code <- unlist(lapply(files, function(x) paste0("zgrep -E ", ensgs_concat, 
                                                          " ", data_dir, "GTEx_Analysis_v7_eQTL/", x)))
    
    gtex_eqtls_mappedGENES <- list()
    for (i in 1:length(zgrep_code)) {
      tryCatch({
        gtex_eqtls_mappedGENES[[i]] <- system(zgrep_code[[i]], intern = TRUE)
        gtex_eqtls_mappedGENES[[i]] <- as.data.frame(gtex_eqtls_mappedGENES[[i]])
        gtex_eqtls_mappedGENES[[i]] <- tidyr::separate(data = gtex_eqtls_mappedGENES[[i]], col = 1, into = c("variant_id", "gene_id",	"tss_distance",	
                                                                                                             "ma_samples",	"ma_count",	"maf",	"pval_nominal",	"slope",	
                                                                                                             "slope_se",	"pval_nominal_threshold",	"min_pval_nominal",	"pval_beta"), sep = "\t")
        gtex_eqtls_mappedGENES[[i]] <- cbind(tissue = gtexNames[i], gtex_eqtls_mappedGENES[[i]])
      }, error = function(i) {})
    }
    names(gtex_eqtls_mappedGENES) <- gtexNames
    gtex_eqtls_mappedGENES_df <- do.call("rbind", gtex_eqtls_mappedGENES)
    gtex_eqtls_mappedGENES_df <- merge(gtex_eqtls_mappedGENES_df, varIDmap, by.x = "variant_id", by.y = "variant_id", all.x = T)
    
    
    bm_genemart <- getBM(attributes = c("ensembl_gene_id_version", "external_gene_name"), 
                         filters = "ensembl_gene_id_version", 
                         values = gtex_eqtls_mappedGENES_df$gene_id, 
                         mart = genemart)
    bm_genemart <- bm_genemart[!duplicated(bm_genemart), ]
    
    #bm_genemart  <- bm_genemart[match(gtex_eqtls_mappedGENES_df$gene_id, bm_genemart$ensembl_gene_id_version),]
    names(bm_genemart) <- c("ensembl_gene_id_version", "mapped_gene_name")
    
    gtex_eqtls_mappedGENES_df <- merge(gtex_eqtls_mappedGENES_df, bm_genemart, by.x = "gene_id", by.y = "ensembl_gene_id_version", all.x = T)
    gtex_eqtls_mappedGENES_df <- gtex_eqtls_mappedGENES_df[!duplicated(gtex_eqtls_mappedGENES_df), ]
    
    gtex_eqtls_mappedGENES_df <- gtex_eqtls_mappedGENES_df[c(3,1:2,19,4:18)]
    gtex_eqtls_mappedGENES_df <- gtex_eqtls_mappedGENES_df[order(gtex_eqtls_mappedGENES_df$mapped_gene_name, gtex_eqtls_mappedGENES_df$tissue),]
    gtex_eqtls_mappedGENES_df$gene_id2 <- gsub("\\..*","",gtex_eqtls_mappedGENES_df$gene_id)
    
    
    df <- data.frame(ensgs = unlist(strsplit(ensgs_clean, split = ";")))
    df$probeid <- unlist(lapply(df$ensgs, function(i) {names(grep(i, ensgs_clean, value = T))}))
    row.names(df) <- NULL
    
    gtex_eqtls_mappedGENES_df <- merge(df,gtex_eqtls_mappedGENES_df,  by.x = "ensgs", by.y = "gene_id2", all.x = TRUE)
    gtex_eqtls_mappedGENES_df$ensgs <- NULL
    
    gtex_eqtls_mappedGENES_df <- gtex_eqtls_mappedGENES_df[order(gtex_eqtls_mappedGENES_df$probeid, gtex_eqtls_mappedGENES_df$pval_nominal),]
    
    for (j in 1:length(names(myInterpretList))) {
      
      item <- names(myInterpretList)[j]
      
      if (item %in% unique(gtex_eqtls_mappedGENES_df$probeid)) {
        
        myInterpretList[[item]]$gtex_cisEqtls$gtex_mappedGenes <- gtex_eqtls_mappedGENES_df[gtex_eqtls_mappedGENES_df$probeid == item,]
        
      }   #if (item %in% unique(gtex_eqtls_mappedGENES_df$probeid))
      
      else{
        
        myInterpretList[[item]]$gtex_cisEqtls$gtex_mappedGenes <- data.frame(matrix(ncol = ncol(gtex_eqtls_mappedGENES_df), nrow = 1))
        colnames(myInterpretList[[item]]$gtex_cisEqtls$gtex_mappedGenes) <-  names(gtex_eqtls_mappedGENES_df)
        
        
      } # end else
      
    } # for (j in 1:length(names(myInterpretList)))
    
    
    
    
    ### STEP 2: mqtl SNPs
    
    gtex_eqtls_mqtlSNPs_colnames <- c("variant_id", "snp_id", "tissue", "gene_id", "tss_distance", "ma_samples", "ma_count", "maf", "pval_nominal", 
                                      "slope", "slope_se", "pval_nominal_threshold", "min_pval_nominal", "pval_beta", "chr", "variant_pos", 
                                      "ref", "alt", "ensembl_gene_id_version", "expressed_gene_name", "expressed_gene_chr", "expressed_gene_start_pos", "expressed_gene_end_pos")
    
    
    for (i in 1:length(myInterpretList)) {
      
      probe_id <- names(myInterpretList)[i]
      
      if (length(myInterpretList[[probe_id]]$mqtls$SNP) == 0) {
        
        gtex_eqtls_mqtlSNPs_df <- data.frame(matrix(ncol = length(gtex_eqtls_mqtlSNPs_colnames), nrow = 0))
        colnames(gtex_eqtls_mqtlSNPs_df) <- gtex_eqtls_mqtlSNPs_colnames
        myInterpretList[[probe_id]]$gtex_cisEqtls$gtex_mqtlSNPs <- gtex_eqtls_mqtlSNPs_df
        
      }  # end of if (length(myInterpretList[[probe_id]]$mqtls$SNP) == 0)
      
      else {
        
        mqtlSNPs <- myInterpretList[[probe_id]]$mqtls$SNP
        mqtlSNPs <- mqtlSNPs[mqtlSNPs != ""]
        mqtlSNPs <- mqtlSNPs[mqtlSNPs != "NA"]
        
        subset_varIDMap <- varIDmap[which(varIDmap$snp_id %in% mqtlSNPs),]
        subset_varIDMap <- subset_varIDMap[match(mqtlSNPs, subset_varIDMap$snp_id),]
        
        variantIDs_query <- subset_varIDMap$variant_id
        variantIDs_query <- variantIDs_query[!is.na(variantIDs_query)]
        
        mqtlSNPs_concat <- paste(variantIDs_query, collapse = "|")
        mqtlSNPs_concat <- paste("\'",mqtlSNPs_concat,"\'", sep="")
        
        
        zgrep_code <- unlist(lapply(files, function(x) paste0("zgrep -wE ", mqtlSNPs_concat, 
                                                              " ", data_dir, "GTEx_Analysis_v7_eQTL/", x)))
        
        gtex_eqtls_mqtlSNPs <- list()
        
        for (i in 1:length(zgrep_code)) {
          tryCatch({
            gtex_eqtls_mqtlSNPs[[i]] <- system(zgrep_code[[i]], intern = TRUE)
            gtex_eqtls_mqtlSNPs[[i]] <- as.data.frame(gtex_eqtls_mqtlSNPs[[i]])
            gtex_eqtls_mqtlSNPs[[i]] <- tidyr::separate(data = gtex_eqtls_mqtlSNPs[[i]], col = 1, into = c("variant_id","gene_id",	"tss_distance",	
                                                                                                           "ma_samples",	"ma_count",	"maf",	"pval_nominal",	"slope",	
                                                                                                           "slope_se",	"pval_nominal_threshold",	"min_pval_nominal",	"pval_beta"), sep = "\t")
            gtex_eqtls_mqtlSNPs[[i]] <- cbind(File = gtexNames[i], gtex_eqtls_mqtlSNPs[[i]])
          }, error = function(i) {})
        } # for (i in 1:length(zgrep_code))
        names(gtex_eqtls_mqtlSNPs) <- gtexNames
        gtex_eqtls_mqtlSNPs_df <- do.call("rbind", gtex_eqtls_mqtlSNPs)
        gtex_eqtls_mqtlSNPs_df <- merge(gtex_eqtls_mqtlSNPs_df, subset_varIDMap, by.x = "variant_id", by.y = "variant_id", all.x = T)
        gtex_eqtls_mqtlSNPs_df <- gtex_eqtls_mqtlSNPs_df[!duplicated(gtex_eqtls_mqtlSNPs_df), ]
        
        if (length(gtex_eqtls_mqtlSNPs_df$variant_id) != 0) {
          snpmart = biomaRt::useMart(biomart = "ENSEMBL_MART_SNP", 
                                     host = "grch37.ensembl.org", path = "/biomart/martservice", 
                                     dataset = "hsapiens_snp")
          
          
          bm_snpmart <- biomaRt::getBM(attributes = c("refsnp_id", "ensembl_gene_stable_id"), 
                                       filters = "snp_filter", 
                                       values = gtex_eqtls_mqtlSNPs_df$snp_id, 
                                       mart = snpmart)
          
          bm_genemart <- biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name","chromosome_name", "start_position", "end_position"), 
                                        filters = "ensembl_gene_id", 
                                        values = bm_snpmart$ensembl_gene_stable_id, 
                                        mart = genemart)
          
          bm_snpmart <- merge(bm_snpmart, bm_genemart, by.x = 'ensembl_gene_stable_id', by.y = 'ensembl_gene_id', all.x = T)
          
          names(bm_snpmart) <- c("ensembl_gene_id_version", "snp_id","expressed_gene_name", "expressed_gene_chr", "expressed_gene_start_pos", "expressed_gene_end_pos")
          
          gtex_eqtls_mqtlSNPs_df <- merge(gtex_eqtls_mqtlSNPs_df, bm_snpmart, by.x = "snp_id", by.y = "snp_id", all.x = T)
          gtex_eqtls_mqtlSNPs_df <- gtex_eqtls_mqtlSNPs_df[!duplicated(gtex_eqtls_mqtlSNPs_df), ]
          
          gtex_eqtls_mqtlSNPs_df <- gtex_eqtls_mqtlSNPs_df[c(2,1,3:length(gtex_eqtls_mqtlSNPs_df))]
          
          gtex_eqtls_mqtlSNPs_df <- gtex_eqtls_mqtlSNPs_df[order(gtex_eqtls_mqtlSNPs_df$snp_id, gtex_eqtls_mqtlSNPs_df$File),]
          names(gtex_eqtls_mqtlSNPs_df)[names(gtex_eqtls_mqtlSNPs_df) == "File"] <- "tissue"
        }
        
        else{
          
          gtex_eqtls_mqtlSNPs_df <- data.frame(matrix(ncol = length(gtex_eqtls_mqtlSNPs_colnames), nrow = 0))
          colnames(gtex_eqtls_mqtlSNPs_df) <- gtex_eqtls_mqtlSNPs_colnames
          
        }
        
        myInterpretList[[probe_id]]$gtex_cisEqtls$gtex_mqtlSNPs <- gtex_eqtls_mqtlSNPs_df
      } # end of else
      
    } # end of for (i in 1:length(myInterpretList)) {
    
  } ## else if (input_genetic_data_type == "probe") 
  
  else {
    warning(paste0("input genetic data type not supported:  ", input_genetic_data_type))
  }
  
  if (input_genetic_data_type %in% c("gene", "snp")) {
    for (j in 1:length(names(myInterpretList))) {
      
      item <- names(myInterpretList)[j]
      
      if (item %in% names(gtex_eqtls_split)) {
        myInterpretList[[item]]$gtex_cisEqtls <- gtex_eqtls_split[[item]]
        
      }
      
      else {
        
        myInterpretList[[item]]$gtex_cisEqtls <- data.frame(matrix(ncol = 23, nrow = 1))
        colnames(myInterpretList[[item]]$gtex_cisEqtls) <-  c("variant_id", "snp_id", "tissue","gene_id",	"tss_distance",	
                                                              "ma_samples",	"ma_count",	"maf",	"pval_nominal",	
                                                              "slope",	"slope_se",	"pval_nominal_threshold",	
                                                              "min_pval_nominal",	"pval_beta", "chr","variant_pos","ref", "alt",
                                                              "ensembl_gene_id_version","expressed_gene_name","expressed_gene_chr",
                                                              "expressed_gene_start_pos", "expressed_gene_end_pos")
      }
    } # if (input_genetic_data_type %in% c("gene", "snp"))
  }
  rm(varIDmap)
  
  return(myInterpretList)
} # DONE
getGwasCatalog <- function(myInterpretList, rangeval = 0, data_dir = "./", input_genetic_data_type = c("gene", "snp", "probe")) {
  
  
  #gwas <- read.delim(url("http://www.ebi.ac.uk/gwas/api/search/downloads/full"), sep="\t", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
  r_val = rangeval
  #gwas$CHR_POS <- as.numeric(gwas$CHR_POS)
  if (input_genetic_data_type == "gene") {
    
    for (i in names(myInterpretList)){
      i_bound <- paste0("\\b", i, "\\b")
      gwasCatalogresults <- gwas[with(gwas, grepl(i_bound, REPORTED.GENE.S.)|grepl(i_bound, MAPPED_GENE)),]
      #gwasCatalogresults <- gwas[with(gwas, grepl(i_bound, REPORTED.GENE.S.)|grepl(i_bound, MAPPED_GENE)|grepl(i_bound, UPSTREAM_GENE_ID)|grepl(i_bound, DOWNSTREAM_GENE_ID)),]
      #gwasCatalogresults_range <- gwas[which((gwas$CHR_ID == as.numeric(myInterpretList[[i]]$location$chromosome_name)) & (gwas$CHR_POS >= myInterpretList[[i]]$location$start_position - rangeval) & (as.numeric(gwas$CHR_POS) <= myInterpretList[[i]]$location$end_position + rangeval)),]
      #gwasCatalogresults <- rbind(gwasCatalogresults, gwasCatalogresults_range)
      gwasCatalogresults <- gwasCatalogresults[!duplicated(gwasCatalogresults),]
      rownames(gwasCatalogresults) <- NULL
      
      myInterpretList[[i]]$gwasCatalog <- gwasCatalogresults
    }
  }
  
  else if (input_genetic_data_type == "snp") {
    for (i in names(myInterpretList)){
      gwasCatalogresults <- gwas[with(gwas, grepl(paste0("^", i,"-.$"), STRONGEST.SNP.RISK.ALLELE)|grepl(paste0("^", i,"$"), SNPS)),]
      #gwasCatalogresults_range <- gwas[which(gwas$CHR_ID == myInterpretList[[i]]$location$chr_name & (gwas$CHR_POS >= myInterpretList[[i]]$location$chrom_start - rangeval & gwas$CHR_POS <= myInterpretList[[i]]$location$chrom_end + rangeval)),]
      #gwasCatalogresults <- rbind(gwasCatalogresults, gwasCatalogresults_range)
      gwasCatalogresults <- gwasCatalogresults[!duplicated(gwasCatalogresults),]
      rownames(gwasCatalogresults) <- NULL
      
      myInterpretList[[i]]$gwasCatalog <- gwasCatalogresults
    }
  }
  
  else if (input_genetic_data_type == "probe") {
    
    ensgs <- unlist(lapply(myInterpretList, function(x) {x$name$ensembl_gene_id} ))
    ensgs_clean <- ensgs[!is.na(ensgs)]
    ensgs_clean <- ensgs_clean[ensgs_clean != ""]
    ensgs_clean <- ensgs_clean[ensgs_clean != "NA"]
    
    df <- data.frame(ensgs = unlist(strsplit(ensgs_clean, split = ";")))
    df$probeid <- unlist(lapply(df$ensgs, function(i) {names(grep(i, ensgs_clean, value = T))}))
    row.names(df) <- NULL
    
    
    genemart = biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                                host = "grch37.ensembl.org", path = "/biomart/martservice", 
                                dataset = "hsapiens_gene_ensembl")
    
    bm_genemart <- biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name","chromosome_name", "start_position", "end_position"), 
                                  filters = "ensembl_gene_id", 
                                  values = df$ensgs, 
                                  mart = genemart)
    
    df <- merge(df,bm_genemart, by.x = 'ensgs', by.y = 'ensembl_gene_id',all.x = T)
    names(df)[names(df) == "ensgs"] <- 'ensembl_gene_id'
    
    gwasCatalogresults_mappedGenes_cols <- c(names(df), names(gwas))
    
    
    for (i in names(myInterpretList)){
      
      
      df_subset <- df[df$probeid == i,]
      
      if(length(df_subset$probeid) == 1) {
        
        geneName = df_subset$external_gene_name
        geneName_bound <- paste0("\\b", geneName, "\\b")
        gwasCatalogresults <- gwas[with(gwas, grepl(geneName_bound, REPORTED.GENE.S.)|grepl(geneName_bound, MAPPED_GENE)),]
        #gwasCatalogresults <- gwas[with(gwas, grepl(geneName_bound, REPORTED.GENE.S.)|grepl(geneName_bound, MAPPED_GENE)|grepl(geneName_bound, UPSTREAM_GENE_ID)|grepl(geneName_bound, DOWNSTREAM_GENE_ID)),]
        #gwasCatalogresults_mappedGenes_range <- gwas[which(gwas$CHR_ID == df_subset$chromosome_name & (as.numeric(gwas$CHR_POS) >= df_subset$start_position) - rangeval & (as.numeric(gwas$CHR_POS) <= as.numeric(df_subset$end_position) + rangeval)),]
        #gwasCatalogresults_mappedGenes <- rbind(gwasCatalogresults_mappedGenes, gwasCatalogresults_mappedGenes_range)
        gwasCatalogresults_mappedGenes <- gwasCatalogresults_mappedGenes[!duplicated(gwasCatalogresults_mappedGenes),]
        rownames(gwasCatalogresults_mappedGenes) <- NULL
        
        df_subset_w_gwas <- cbind(df_subset, gwasCatalogresults_mappedGenes)
        
        
      } else if (length(df_subset$probeid) > 1) {
        
        df_subset_w_gwas = data.frame()
        
        for (j in df_subset$external_gene_name){
          
          df_subset2 <- df_subset[df_subset$external_gene_name == j,]
          j_bound <- paste0("\\b", j, "\\b")
          gwasCatalogresults <- gwas[with(gwas, grepl(j_bound, REPORTED.GENE.S.)|grepl(j_bound, MAPPED_GENE)),]
          #gwasCatalogresults_mappedGenes <- gwas[with(gwas, grepl(j, REPORTED.GENE.S.)|grepl(j, MAPPED_GENE)|grepl(j, UPSTREAM_GENE_ID)|grepl(j, DOWNSTREAM_GENE_ID)),]
          #gwasCatalogresults_mappedGenes_range <- gwas[which(gwas$CHR_ID == df_subset2$chromosome_name & (as.numeric(gwas$CHR_POS) >= df_subset2$start_position) - rangeval & (as.numeric(gwas$CHR_POS) <= as.numeric(df_subset2$end_position) + rangeval)),]
          #gwasCatalogresults_mappedGenes <- rbind(gwasCatalogresults_mappedGenes, gwasCatalogresults_mappedGenes_range)
          gwasCatalogresults_mappedGenes <- gwasCatalogresults_mappedGenes[!duplicated(gwasCatalogresults_mappedGenes),]
          rownames(gwasCatalogresults_mappedGenes) <- NULL
          
          df_subset_temp <- cbind(df_subset2, gwasCatalogresults_mappedGenes)
          df_subset_w_gwas<- rbind(df_subset_w_gwas, df_subset_temp)
          
        } # for (i in df_subset$external_gene_name)
        
        
        
      } else {
        
        df_subset_w_gwas <- data.frame(matrix(ncol = length(gwasCatalogresults_mappedGenes_cols), nrow = 1))
        colnames(df_subset_w_gwas) <-  gwasCatalogresults_mappedGenes_cols
        
      } # else
      
      
      myInterpretList[[i]]$gwasCatalog$gwasCatalog_mappedGenes <- df_subset_w_gwas
      
      
      #### GWAS CATALOG RESULTS FOR MQTL SNPS
      
      mqtls = myInterpretList[[i]]$mqtls
      mqtls = mqtls[mqtls$SNP != "NA",]
      if (length(mqtls$SNP) != 0) {
        
        gwasCatalogresults_mqtlSNPs_ALL <- data.frame()
        for (j in 1:length(mqtls$SNP)) {
          gwasCatalogresults_mqtlSNPs <- gwas[with(gwas, grepl(paste0("^", mqtls$SNP[j],"-.$"), STRONGEST.SNP.RISK.ALLELE)|grepl(paste0("^", mqtls$SNP[j],"$"), SNPS)),]
          #gwasCatalogresults_mqtlSNPs_range <- gwas[which(gwas$CHR_ID == mqtls$SNP_CHR[[j]] & (as.numeric(gwas$CHR_POS) >= as.numeric(mqtls$SNP_BP[[j]]) - rangeval & as.numeric(gwas$CHR_POS) <= as.numeric(mqtls$SNP_BP[[j]]) + rangeval)),]
          #gwasCatalogresults_mqtlSNPs <- rbind(gwasCatalogresults_mqtlSNPs, gwasCatalogresults_mqtlSNPs_range)
          gwasCatalogresults_mqtlSNPs <- gwasCatalogresults_mqtlSNPs[!duplicated(gwasCatalogresults_mqtlSNPs),]
          rownames(gwasCatalogresults_mqtlSNPs) <- NULL
          gwasCatalogresults_mqtlSNPs_ALL <- rbind(gwasCatalogresults_mqtlSNPs_ALL,gwasCatalogresults_mqtlSNPs)
        }
        
        myInterpretList[[i]]$gwasCatalog$gwasCatalog_mqtlSNPs <- gwasCatalogresults_mqtlSNPs_ALL
        
        
      }
      
      else {
        
        myInterpretList[[i]]$gwasCatalog$gwasCatalog_mqtlSNPs <- data.frame(matrix(ncol = ncol(gwas), nrow = 1))
        colnames(myInterpretList[[i]]$gwasCatalog$gwasCatalog_mqtlSNPs) <-  names(gwas)
      }
      
    }# for (p in names(myInterpretList)){
  } ## elseif
  
  
  else {
    
    warning(paste0("input genetic data type not supported:  ", input_genetic_data_type))
    
  }
  
  return(myInterpretList)
} # DONE
getGoCategories <- function(myInterpretList, input_genetic_data_type  = c("gene", "snp", "probe")) {
  
  if(input_genetic_data_type == "gene") {
    
    ensgs <- unlist(lapply(myInterpretList, function(x) {x$name$external_gene_name} ))
    
    goIDlist <- getgo(ensgs, genome = 'hg19', id = "geneSymbol", fetch.cats=c("GO:CC","GO:BP","GO:MF"))
    
    for (i in 1:length(names(myInterpretList))) {
      
      gName = toupper(names(myInterpretList)[[i]])
      
      if (gName %in% names(goIDlist)) {
        
        myInterpretList[[i]]$goCategories <- select(GO.db, keys = goIDlist[[gName]], columns = c("TERM","ONTOLOGY"), keytype="GOID")
      }
      
      else {
        
        myInterpretList[[i]]$goCategories <- data.frame(GOID = NA, TERM = NA, ONTOLOGY = NA)
      }
    }
  }
  
  else if (input_genetic_data_type == "snp") {
    
    ensgs <- unlist(lapply(myInterpretList, function(x) {x$name$gene_name} ))
    ensgs_clean <- ensgs[which(!(ensgs == "NA"))]
    ensgs_clean <- strsplit(ensgs_clean, split = ";")

    for (i in 1:length(names(ensgs_clean)) ){
      snpname <- names(ensgs_clean)[i]
      
      goIDlist <- getgo(ensgs_clean[[snpname]], genome = 'hg19', id = "geneSymbol", fetch.cats=c("GO:CC","GO:BP","GO:MF"))
      goIDlist <- goIDlist[which(!(is.null(goIDlist)))]
      
      if (snpname == names(myInterpretList)[[i]]){
        
        
        goCatdf <- data.frame()
        for (j in 1:length(goIDlist)) {
          
          goCattemp <- select(GO.db, keys = goIDlist[[j]], columns = c("TERM","ONTOLOGY"), keytype="GOID")
          goCattemp$GENE <- names(goIDlist)[j]
          goCattemp <- goCattemp[c(4,1:3)] 
          
          goCatdf <- rbind(goCatdf, goCattemp)
        } # end of j for loop
        
        myInterpretList[[i]]$goCategories <- goCatdf
      } #end of if snpname ==
      
      else {
        
        myInterpretList[[i]]$goCategories <- data.frame(GENE = NA, GOID = NA, TERM = NA, ONTOLOGY = NA)
        
      } #end of else
    } # end of i for loop
    
  } # end of else if
  
  else if (input_genetic_data_type == "probe") {
    
    ensgs <- unlist(lapply(myInterpretList, function(x) {x$name$UCSC_RefGene_Name}) )
    df <- data.frame(ensgs = unlist(strsplit(ensgs, split = ";")))
    df$probeid <- unlist(lapply(df$ensgs, function(i) {names(grep(i, ensgs, value = T))}))
    row.names(df) <- NULL
    
    goIDlist <- getgo(df$ensgs, genome = 'hg19', id = "geneSymbol", fetch.cats=c("GO:CC","GO:BP","GO:MF"))
    goIDdf <- data.frame(mappedGene = rep(names(goIDlist), sapply(goIDlist, length)),
                         goIDs = unlist(goIDlist))
    rownames(goIDdf) <- NULL
    goIDdf <- merge(df, goIDdf, by.x = "ensgs", by.y = "mappedGene", all.y = TRUE)
    goCategoriesDF <- select(GO.db, keys = goIDdf$goIDs, columns = c("TERM","ONTOLOGY"), keytype="GOID")
    goCategoriesDF <- goCategoriesDF[!duplicated(goCategoriesDF),]
    goCategoriesDF <- merge(goIDdf, goCategoriesDF, by.x = "goIDs", by.y = "GOID", all.x = T)
    names(goCategoriesDF)[names(goCategoriesDF) == "ensgs"] <- "mappedGene"
    
    
    for (i in names(myInterpretList)) {
      
      if (i %in% unique(goCategoriesDF$probeid)) {
        
        probe_go_subset <- goCategoriesDF[goCategoriesDF$probeid == i,]
        probe_go_subset <- probe_go_subset[order(probe_go_subset$mappedGene, probe_go_subset$ONTOLOGY, probe_go_subset$goIDs),]
        myInterpretList[[i]]$goCategories <- probe_go_subset[c("probeid", "mappedGene", "goIDs", "TERM", "ONTOLOGY")]
        
        
      }
      
      else {
        
        myInterpretList[[i]]$goCategories <- data.frame(probeid = NA, mappedGene = NA, goIDs = NA, TERM = NA, ONTOLOGY = NA)
        
      }
      
    }
    
  } ## else if probe
  
  
  else {
    
    warning(paste0("input genetic data type not supported:  ", input_genetic_data_type))
    
  }
  
  return(myInterpretList)
}   # DONE
getmQTLs <- function(myInterpretList, data_dir = "/data/Priya/genetica/mqtl_files/", input_genetic_data_type = c("snp", "probe")){
  
  if (!(input_genetic_data_type %in% c('snp', 'probe'))) {
    
    stop("input_genetic_data_type not supported: ", input_genetic_data_type)
  }
  
  mqtl_files <- list.files(paste(data_dir, "/final_mqtl_files/", 
                                 sep = ""), pattern = "_mQTLs_edited_cols")
  
  mqtlNames <- gsub("_mQTLs_edited_cols.txt", "", mqtl_files)
  
  ensgs <- names(myInterpretList)
  ensgs_concat <- paste(ensgs, collapse = "|")
  ensgs_concat <- paste("\'",ensgs_concat,"\'", sep="")
  
  
  grep_code <- unlist(lapply(mqtl_files, function(x) paste0("grep -wE ", ensgs_concat, 
                                                            " ", data_dir, "final_mqtl_files/", x)))
  
  mqtls <- list()
  for (i in 1:length(grep_code)) {
    tryCatch({
      mqtls[[i]] <- system(grep_code[[i]], intern = TRUE)
      mqtls[[i]] <- as.data.frame(mqtls[[i]])
      mqtls[[i]] <- tidyr::separate(data = mqtls[[i]], col = 1, into = c("SNP", "SNP_CHR", "SNP_BP", "SNP_CHR_POS", "DNA_M_PROBEID", "DNA_M_PROBE_CHR", "DNA_M_PROBE_BP", "MQTL_BETA", "MQTL_P_VALUE"), sep = "\t")
      mqtls[[i]] <- cbind(File = mqtlNames[i], mqtls[[i]])
    }, error = function(i) {})
  }
  names(mqtls) <- mqtlNames
  mqtls_df <- do.call("rbind", mqtls)
  row.names(mqtls_df) <- c(1:length(mqtls_df$File))
  
  
  
  if (input_genetic_data_type == "snp") {
    
    mqtls_split <- list()
    for (i in ensgs) {
      mqtls_split[[i]] <- mqtls_df[grep(i, mqtls_df$SNP), ]
      mqtls_split[[i]] <- mqtls_split[[i]][order(mqtls_split[[i]]$File),]
    }
    
  } ##  if (input_genetic_data_type == "snp")
  
  if (input_genetic_data_type == "probe") {
    
    mqtls_split <- list()
    for (i in ensgs) {
      mqtls_split[[i]] <- mqtls_df[grep(i, mqtls_df$DNA_M_PROBEID), ]
      mqtls_split[[i]] <- mqtls_split[[i]][order(mqtls_split[[i]]$File),]
    } 
  } ## if (input_genetic_data_type == "probe")
  
  
  names(mqtls_split) <- ensgs
  
  for (j in 1:length(names(myInterpretList))) {
    
    item <- names(myInterpretList)[j]
    
    if (item %in% names(mqtls_split)) {
      myInterpretList[[item]]$mqtls <- mqtls_split[[item]]
      
    }
  }
  
  return(myInterpretList)
} ## DONE
getPGCTopHits <- function(myInterpretList, snp_disorders = c("adhd", "asd", "bip", "mdd", "scz"), data_dir = "/data/Priya/genetica/top_pgc_results_files/", input_genetic_data_type = c("gene", "snp", "probe")) {
  
  
  if (input_genetic_data_type  == "gene") {
    
    gene_pgc_top_hits <- read.table(paste0(data_dir,"Top_gene_hits_adhd_and_cross_disorders.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    
    for (i in 1:length(myInterpretList)) {
      gName = unique(myInterpretList[[i]]$name$external_gene_name)
      myInterpretList[[i]]$pgc_top_hits <- gene_pgc_top_hits[gene_pgc_top_hits$GENE == gName,]
      
    }
  } ## end if (input_genetic_data_type  = "gene")
  
  
  else if (input_genetic_data_type  == "snp") {
    
    sName = names(myInterpretList)
    
    pgc_files <- list.files(paste(data_dir, "meta_analysis_results", 
                                  sep = ""), pattern = "_eur_PB_edited_fname.txt")
    
    pgcNames <- gsub("_eur_PB_edited_fname.txt", "", pgc_files)
    file_cols <- c("SNP", "CHR", "BP", "A1", "A2", "INFO", "OR", "P", "SE")
    
    pgcToQuery <- sapply(snp_disorders, function(i) {grep(i, pgc_files, value = TRUE)})
    
    pgcgwas <- lapply(pgcToQuery, function(x) { read.delim(paste0(data_dir,"meta_analysis_results/",x), header = TRUE, stringsAsFactors = FALSE)})
    pgcgwas <- lapply(pgcgwas, function(x) { x[file_cols]})
    
    
    for (i in 1:length(names(myInterpretList))){
      sName = names(myInterpretList)[i]
      pgc_hit <- lapply(pgcgwas, function(x) { x[x$SNP == sName,]})
      
      for (j in 1:length(pgc_hit)) {
        
        d <- toupper(names(pgc_hit[j]))
        colnames(pgc_hit[[j]])[c(6:9)] <- paste(d, file_cols[c(6:9)], sep = "_")
      }
      
      pgc_hits_results <- Reduce(function(x, y)  merge(x, y, by= c("SNP", "CHR", "BP", "A1", "A2"),  all=TRUE), pgc_hit)
      
      tophits <- c()
      for (p in grep("_P",names(pgc_hits_results), value = T)){
        tophits[[p]] <- ifelse(pgc_hits_results[[p]] < 1e-05 , d, NA)
        
      }
      
      if  (all(is.na(tophits))){
        
        pgc_hits_results$TopHitDisorders = NA
      }
      
      else {
        
        tophits <- tophits[is.na(tophits)]
        pgc_hits_results$TopHitDisorders  <- paste(names(tophits), collapse = ";")
        
      }
      
      pgc_hits_results <- pgc_hits_results[c(1:5,26,6:25)]
      
      myInterpretList[[i]]$pgc_top_hits <- pgc_hits_results
    } ## end for (i in 1:length(names(myInterpretList)))
    
    
  } ## end else if (input_genetic_data_type  = "snp")
  
  else if (input_genetic_data_type  == "probe") {
    
    
    ###STEP 1 GENES:
    
    gene_pgc_top_hits <- read.table(paste0(data_dir,"Top_gene_hits_adhd_and_cross_disorders.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    
    gNames <- unlist(lapply(myInterpretList, function(x) {x$name$UCSC_RefGene_Name} ))
    df <- data.frame(mappedGenes = unlist(strsplit(gNames, split = ";")))
    df$probeid <- unlist(lapply(df$mappedGenes, function(i) {names(grep(i, gNames, value = T))}))
    row.names(df) <- NULL
    
    tophits_df <- merge(df,gene_pgc_top_hits,  by.x = "mappedGenes", by.y = "GENE")
    
    for (i in names(myInterpretList)) {
      
      if (i %in% unique(tophits_df$probeid)) {
        
        probe_pgc_subset <- tophits_df[tophits_df$probeid == i,]
        #probe_pgc_subset <- probe_pgc_subset[order(probe_go_subset$mappedGene, probe_go_subset$ONTOLOGY, probe_go_subset$goIDs),]
        #myInterpretList[[i]]$pgc_top_hits$pgc_top_hits_mappedGenes <- probe_pgc_subset
        
        
      }
      
      else {
        
        myInterpretList[[i]]$pgc_top_hits$pgc_top_hits_mappedGenes <- data.frame(matrix(ncol = ncol(tophits_df), nrow = 1))
        colnames(myInterpretList[[i]]$pgc_top_hits$pgc_top_hits_mappedGenes) <-  names(tophits_df)
        
      }
      
    }
    
    
    ##STEP 2 SNPS
    pName = names(myInterpretList)
    
    pgc_files <- list.files(paste(data_dir, "meta_analysis_results", 
                                  sep = ""), pattern = "_eur_PB_edited_fname.txt")
    
    pgcNames <- gsub("_eur_PB_edited_fname.txt", "", pgc_files)
    file_cols <- c("SNP", "CHR", "BP", "A1", "A2", "INFO", "OR", "P", "SE")
    
    pgcToQuery <- sapply(snp_disorders, function(i) {grep(i, pgc_files, value = TRUE)})
    
    pgcgwas <- lapply(pgcToQuery, function(x) { read.delim(paste0(data_dir,"meta_analysis_results/",x), header = TRUE, stringsAsFactors = FALSE)})
    pgcgwas <- lapply(pgcgwas, function(x) { x[file_cols]})
    
    for (i in 1:length(names(myInterpretList))){
      
      
      sName <- myInterpretList[[i]]$mqtls$SNP
      pgc_hit <- lapply(pgcgwas, function(x) { x[x$SNP %in% sName,]})
      
      for (j in 1:length(pgc_hit)) {
        
        d <- toupper(names(pgc_hit[j]))
        colnames(pgc_hit[[j]])[c(6:9)] <- paste(d, file_cols[c(6:9)], sep = "_")
      }
      
      pgc_hits_results <- Reduce(function(x, y)  merge(x, y, by= c("SNP", "CHR", "BP", "A1", "A2"),  all=TRUE), pgc_hit)
      
      
      pgc_hits_results$TopHitDisorders <- apply(pgc_hits_results[grep("_P",names(pgc_hits_results))], 1, function(x) paste((names(which(x < 1e-05 ))), collapse = ";"))
      pgc_hits_results$TopHitDisorders <- ifelse(pgc_hits_results$TopHitDisorders == "", NA, pgc_hits_results$TopHitDisorders)
      
      
      pgc_hits_results <- pgc_hits_results[c(1:5,26,6:25)]
      
      myInterpretList[[i]]$pgc_top_hits$pgc_top_hits_mqtlSNPs <- pgc_hits_results
    } ## for (i in 1:length(names(myInterpretList)))
    
    
    
    
    
  } ## else if (input_genetic_data_type  = "probe")
  
  else {
    
    stop("input_genetic_data_type not supported: ", input_genetic_data_type)
  } 
  
  return(myInterpretList)
  
}

combine_genetica_genes <- function(myInterpretList, results_fname_xlsx = "genetica_results_genes.xlsx"){
  
  #############
  # NAME LOC
  #############
  
  #extract name data
  extract.name <- lapply(myInterpretList, function(y) y[["name"]])
  extract.name <- do.call(rbind, extract.name)
  
  #extract 
  extract.location <- lapply(myInterpretList, function(y) y[["location"]])
  extract.location <- do.call(rbind, extract.location)
  
  final_name_location_df <- merge(extract.name, extract.location[c(2:8)], by.x = "row.names", by.y = "row.names")
  rownames(final_name_location_df) <- final_name_location_df$Row.names
  final_name_location_df$Row.names <- NULL
  colnames(final_name_location_df) <- c("gene", "ensembl_gene_id", "description","gene_chr_hg19", "gene_start_hg19", "gene_end_hg19", "gene_strand", "gene_chr_hg38", "gene_start_hg38", "gene_end_hg38")
  
  #############
  # gTex results
  #############
  extract.gtex <- lapply(myInterpretList, function(y) y[["gtex_cisEqtls"]])
  
  for (i in 1:length(extract.gtex)) {
    if(dim(extract.gtex[[i]])[1] == 0) {
      
    colnames(extract.gtex[[i]]) =  c("File","variant_id" ,"gene_id", "tss_distance","ma_samples","ma_count","maf",
                                     "pval_nominal", "slope","slope_se","pval_nominal_threshold", "min_pval_nominal","pval_beta",  "chr", "variant_pos",
                                     "ref","alt","snp_id")
    }
    
    if(is.na(extract.gtex[[i]]$variant_id[1])) {
      
      extract.gtex[[i]] = data.frame(File = NA,variant_id = NA ,gene_id = NA, tss_distance = NA, ma_samples = NA, ma_count = NA,maf = NA,
                                       pval_nominal = NA, slope = NA,slope_se = NA,pval_nominal_threshold = NA, min_pval_nominal = NA, pval_beta = NA, 
                                      chr = NA, variant_pos = NA, ref = NA,alt = NA, snp_id = NA)
    }
    
    
  }
  extract.gtex <- do.call(rbind, extract.gtex)
  extract.gtex$gene <- gsub("\\..*", "", row.names(extract.gtex))
  
  #In the case of SNPs with multiple gtex expression results, order the pvalues by most significant to least for each SNP
  extract.gtex <- extract.gtex[order(extract.gtex$gene,extract.gtex$snp_id, extract.gtex$pval_nominal), ]
  
  extract.gtex <- ddply(extract.gtex, .(snp_id), summarize,
                        gene = paste(unique(gene),collapse=";"),
                        variant_id = paste(unique(variant_id),collapse=";"), 
                        chr = paste(unique(chr),collapse=";"), 
                        variant_pos = paste(unique(variant_pos),collapse=";"),
                        ref = paste(unique(ref),collapse=";"),
                        alt = paste(unique(alt),collapse=";"),
                        tss_distance = paste(unique(tss_distance),collapse=";"),
                        pval_nominal = paste(pval_nominal,collapse=";"),
                        slope = paste(slope,collapse=";"),
                        tissue = paste(File,collapse=";")
  )
  extract.gtex <- merge(final_name_location_df, extract.gtex, by.x = "gene", by.y = "gene", all.y = T)
  extract.gtex <- extract.gtex[c(11,1:10,12:20)]
  #rename columns
  colnames(extract.gtex) <- c("snp_rs_id", "gene_expressed", "ensembl_gene_id", "description","gene_chr_hg19", "gene_start_hg19", "gene_end_hg19", "gene_strand", 
                              "gene_chr_hg38", "gene_start_hg38", "gene_end_hg38",
                              "variant_id", "variant_chr", "variant_pos", "variant_ref", "variant_alt", 
                              "gtex_tss_distance", "gtex_pvalue", "gtex_NES", "gtex_tissue")
  
  #Order the results based on SNP position
  extract.gtex <- extract.gtex[order(extract.gtex$gene_expressed, extract.gtex$variant_chr, extract.gtex$variant_pos), ]
  
  #Flag indicating if the variant is within the gene 
  extract.gtex$is_SNP_in_gene <- ifelse(as.numeric(extract.gtex$variant_pos) > as.numeric(extract.gtex$gene_start_hg19) & as.numeric(extract.gtex$variant_pos) < as.numeric(extract.gtex$gene_end_hg19), TRUE, FALSE) 
  
  #############
  # GWAS Catalog
  #############
  extract.gwascatalog <- lapply(myInterpretList, function(y) y[["gwasCatalog"]])
  extract.gwascatalog <- do.call(rbind, extract.gwascatalog)
  extract.gwascatalog$gene <- gsub("\\..*", "", row.names(extract.gwascatalog))
  extract.gwascatalog <- extract.gwascatalog[c(length(extract.gwascatalog), 1: (length(extract.gwascatalog) - 1))]
  
  extract.gwascatalog <- merge(final_name_location_df, extract.gwascatalog, by.x = "gene", by.y = "gene", all.y = T)
  
  #############
  # GO Categories
  #############
  extract.gocategories <- lapply(myInterpretList, function(y) y[["goCategories"]])
  extract.gocategories <- do.call(rbind, extract.gocategories)
  extract.gocategories$gene <- gsub("\\..*", "", row.names(extract.gocategories))
  extract.gocategories <- merge(final_name_location_df, extract.gocategories, by.x = "gene", by.y = "gene", all.y = T)
  
  
  #############
  # PGC Top Hits
  #############  
  extract.pgchits <- lapply(myInterpretList, function(y) y[["pgc_top_hits"]])
  extract.pgchits <- do.call(rbind, extract.pgchits)
  extract.pgchits <- merge(final_name_location_df, extract.pgchits, by.x = "gene", by.y = "GENE", all.y = T)
  
  
  
  ##########
  # SUMMARY PAGE
  ##########
  # A true or false table indicating if results are present in Gtex, gwas catalog or pgc top hits for each queried gene
  
  master <- data.frame(final_name_location_df)
  
  master$has_GTex <- master$gene %in% extract.gtex$gene_expressed
  master$has_GwasCatalog <- master$gene %in% extract.gwascatalog$gene
  master$has_goCategories <- master$gene %in% extract.gocategories$gene
  master$has_PGC_hit <- master$gene %in% extract.pgchits$GENE
  
  ########
  # Notification of no results
  ########
  
  results <- list(master,extract.gtex,extract.gwascatalog,extract.gocategories,extract.pgchits)
  names(results) <- c("Summary", "GTex_hg19", "GWASCatalog_hg38", "GOCategories_hg19", "PGCHits_hg19")
  
  for (i in 1:length(results)) {
    
    if (nrow(results[[i]]) == 0){
      x <- paste("There were no", names(results)[i], "results for your query", sep = " ")
      results[[i]] <- data.frame(x)
      
    }
    
  }
  
  ########
  # Write in one excel file
  ########
  
  openxlsx::write.xlsx(results, file = results_fname_xlsx)
  
  
  return(results)
  
} ## END OF FUNCTION
combine_genetica_snps <- function(myInterpretList, results_fname_xlsx = "genetica_results_snps.xlsx"){
  
  #############
  # NAME LOC
  #############
  
  #extract name data
  extract.name <- lapply(myInterpretList, function(y) y[["name"]])
  extract.name <- do.call(rbind, extract.name)
  
  #extract 
  extract.location <- lapply(myInterpretList, function(y) y[["location"]])
  extract.location <- do.call(rbind, extract.location)
  extract.location$chrom_end_hg19 <- NULL 
  extract.location$chrom_end_hg38 <- NULL 
  
  final_name_location_df <- merge(extract.location, extract.name[c(2:10)], by.x = "row.names", by.y = "row.names")
  rownames(final_name_location_df) <- final_name_location_df$Row.names
  final_name_location_df$Row.names <- NULL
  colnames(final_name_location_df) <- c("snp_id", "snp_chr_hg19", "snp_pos_hg19", "snp_strand", "snp_chr_hg38", "snp_pos_hg38","allele", 
                                        "ensembl_gene_stable_id", "gene_name", "gene_chr_hg19", "gene_start_hg19", "gene_end_hg19", "gene_chr_hg38", "gene_start_hg38", "gene_end_hg38")
  
  #############
  # mQTL Membership
  #############
  
  extract.mQTL <- lapply(myInterpretList, function(y) y[["mqtls"]])
  extract.mQTL <- do.call(rbind, extract.mQTL)
  names(extract.mQTL)[names(extract.mQTL) == "File"] <- "mQTL_dataset"
  
  #extract.mQTL$snp_id <- gsub("\\..*", "", row.names(extract.mQTL))
  
  
  #############
  # gTex results
  #############
  extract.gtex <- lapply(myInterpretList, function(y) y[["gtex_cisEqtls"]])
  extract.gtex <- do.call(rbind, extract.gtex)
  #extract.gtex$snp_ <- gsub("\\..*", "", row.names(extract.gtex))
  
  extract.gtex <- extract.gtex[order(extract.gtex$snp_id, extract.gtex$pval_nominal), ]
  extract.gtex <- merge(final_name_location_df, extract.gtex, by.x = "snp_id", by.y = "snp_id")
  #extract.gtex <- extract.gtex[c(1:23,28:36)]
  
  #############
  # GWAS Catalog
  #############
  extract.gwascatalog <- lapply(myInterpretList, function(y) y[["gwasCatalog"]])
  extract.gwascatalog <- do.call(rbind, extract.gwascatalog)
  extract.gwascatalog$snp_id <- gsub("\\..*", "", row.names(extract.gwascatalog))
  extract.gwascatalog <- extract.gwascatalog[c(length(extract.gwascatalog), 1: (length(extract.gwascatalog) - 1))]
  
  extract.gwascatalog <- merge(final_name_location_df, extract.gwascatalog, by.x = "snp_id", by.y = "snp_id", all.y = T)
  
  #############
  # GO Categories
  #############
  extract.gocategories <- lapply(myInterpretList, function(y) y[["goCategories"]])
  extract.gocategories <- do.call(rbind, extract.gocategories)
  extract.gocategories$snp_id <- gsub("\\..*", "", row.names(extract.gocategories))
  extract.gocategories <- merge(final_name_location_df, extract.gocategories[c(2:5)], by.x = "snp_id", by.y = "snp_id", all.y = T)
  
  
  #############
  # PGC Top Hits
  #############  
  extract.pgchits <- lapply(myInterpretList, function(y) y[["pgc_top_hits"]])
  extract.pgchits <- do.call(rbind, extract.pgchits)
  
  
  
  ##########
  # SUMMARY PAGE
  ##########
  # A true or false table indicating if results are present in Gtex, gwas catalog or pgc top hits for each queried gene
  
  master <- data.frame(final_name_location_df)
  
  master$has_GTex <- master$snp_id %in% extract.gtex$snp_id
  master$has_mQTLs <- master$snp_id %in% extract.mQTL$SNP
  master$has_GwasCategories <- master$snp_id %in% extract.gwascatalog$snp_id
  master$has_goCategories <- master$snp_id %in% extract.gocategories$snp_id
  
  master <- merge(master, extract.pgchits[c("SNP", "TopHitDisorders")], by.x = "snp_id", by.y = "SNP", all.x = TRUE)
  names(master)[names(master) == "TopHitDisorders"] <- "has_PGC_hit"
  master$has_PGC_hit <- ifelse(is.na(master$has_PGC_hit), FALSE, TRUE)
  
  
  ########
  # Notification of no results
  ########
  
  results <- list(master,extract.mQTL,extract.gtex,extract.gwascatalog,extract.gocategories,extract.pgchits)
  names(results) <- c("Summary", "mQTLs", "GTex", "GWASCatalog", "GOCategories", "PGCHits")
  
  for (i in 1:length(results)) {
    
    if (nrow(results[[i]]) == 0){
      
      results[[i]] <- paste("There were no", names(results)[i], "results for your query", sep = " ")
      
    }
    
  }
  
  ########
  # Write in one excel file
  ########
  
  write.xlsx(results, file = results_fname_xlsx)
  
  
  return(results)
  
} ## END OF FUNCTION
combine_genetica_probes <- function(myInterpretList, results_fname_xlsx = "genetica_results_probes.xlsx"){
  
  #############
  # NAME LOC
  #############
  
  #extract name data
  extract.name <- lapply(myInterpretList, function(y) y[["name"]])
  extract.name <- do.call(rbind, extract.name)
  
  #extract location data
  extract.location <- lapply(myInterpretList, function(y) y[["location"]])
  extract.location <- do.call(rbind, extract.location)
  names(extract.location)[names(extract.location) == "MAPINFO"] <- "BP"
  
  final_name_location_df <- merge(extract.location, extract.name[c(2:6)], by.x = "row.names", by.y = "row.names")
  rownames(final_name_location_df) <- final_name_location_df$Row.names
  final_name_location_df$Row.names <- NULL
  colnames(final_name_location_df)[1:3] <- c("probe_id", "probe_chr", "probe_pos")
  
  #############
  # mQTL Membership
  #############
  
  extract.mQTL <- lapply(myInterpretList, function(y) y[["mqtls"]])
  extract.mQTL <- do.call(rbind, extract.mQTL)
  names(extract.mQTL)[names(extract.mQTL) == "File"] <- "mQTL_dataset"
  extract.mQTL <- extract.mQTL[c(1,6:8,2:5,9:15)]
  extract.mQTL <- extract.mQTL[order(extract.mQTL$mQTL_dataset, extract.mQTL$DNA_M_PROBEID, extract.mQTL$SNP_CHR, extract.mQTL$SNP_BP),]
  rownames(extract.mQTL) <- NULL
  #############
  # gTex results
  #############
  
  ### MQTL SNPS ###
  extract.gtex.snps <- lapply(myInterpretList, function(y) y$gtex_cisEqtls[["gtex_mqtlSNPs"]])
  extract.gtex.snps <- do.call(rbind, extract.gtex.snps)
  extract.gtex.snps <- extract.gtex.snps[complete.cases(extract.gtex.snps$variant_id),]
  extract.gtex.snps$probe_id <- gsub("\\..*", "", row.names(extract.gtex.snps))
  
  extract.gtex.snps$pval_nominal <- as.numeric(extract.gtex.snps$pval_nominal)
  extract.gtex.snps <- extract.gtex.snps[order(extract.gtex.snps$snp_id, extract.gtex.snps$pval_nominal), ]
  
  extract.gtex.snps <- ddply(extract.gtex.snps, .(snp_id), summarize,
                             probe_id = paste(unique(probe_id),collapse=";"),
                             mapped_gene = paste(unique(expressed_gene_name),collapse=";"),
                             variant_id = paste(unique(variant_id),collapse=";"), 
                             variant_chr = paste(unique(chr),collapse=";"), 
                             variant_pos = paste(unique(variant_pos),collapse=";"),
                             ref = paste(unique(ref),collapse=";"),
                             alt = paste(unique(alt),collapse=";"),
                             tss_distance = paste(unique(tss_distance),collapse=";"),
                             pval_nominal = paste(pval_nominal,collapse=";"),
                             slope = paste(slope,collapse=";"),
                             tissue = paste(tissue,collapse=";"),
                             ensembl_mapped_gene_id_version = paste(unique(ensembl_gene_id_version),collapse=";"),
                             mapped_gene_chr = paste(unique(expressed_gene_chr),collapse=";"),
                             mapped_gene_start_pos = paste(unique(expressed_gene_start_pos),collapse=";"),
                             mapped_gene_end_pos = paste(unique(expressed_gene_end_pos),collapse=";"))
  if (length(extract.gtex.snps$variant_id)!= 0){
    extract.gtex.snps <- extract.gtex.snps[order(extract.gtex.snps$probe_id, extract.gtex.snps$variant_chr, extract.gtex.snps$variant_pos), ]
  }
  
  ### MAPPED GENES ###
  extract.gtex.mappedgenes <- lapply(myInterpretList, function(y) y$gtex_cisEqtls[["gtex_mappedGenes"]])
  extract.gtex.mappedgenes <- do.call(rbind, extract.gtex.mappedgenes)
  extract.gtex.mappedgenes <- extract.gtex.mappedgenes[complete.cases(extract.gtex.mappedgenes$tissue),]
  
  extract.gtex.mappedgenes$probe_id <- gsub("\\..*", "", row.names(extract.gtex.mappedgenes))
  
  #In the case of SNPs with multiple gtex expression results, order the pvalues by most significant to least for each SNP
  extract.gtex.mappedgenes <- extract.gtex.mappedgenes[order(extract.gtex.mappedgenes$probe_id,extract.gtex.mappedgenes$snp_id, extract.gtex.mappedgenes$pval_nominal), ]
  
  extract.gtex.mappedgenes <- ddply(extract.gtex.mappedgenes, .(snp_id), summarize,
                                    probe_id = paste(unique(probe_id),collapse=";"),
                                    mapped_gene_name = paste(unique(mapped_gene_name),collapse=";"),
                                    gene_id = paste(unique(gene_id),collapse=";"),
                                    variant_id = paste(unique(variant_id),collapse=";"), 
                                    variant_chr = paste(unique(chr),collapse=";"), 
                                    variant_pos = paste(unique(variant_pos),collapse=";"),
                                    ref = paste(unique(ref),collapse=";"),
                                    alt = paste(unique(alt),collapse=";"),
                                    tss_distance = paste(unique(tss_distance),collapse=";"),
                                    pval_nominal = paste(pval_nominal,collapse=";"),
                                    slope = paste(slope,collapse=";"),
                                    tissue = paste(tissue,collapse=";")
  )
  
  #Order the results based on probe id and SNP position
  if (length(extract.gtex.mappedgenes$variant_id)!= 0){
    extract.gtex.mappedgenes <- extract.gtex.mappedgenes[order(extract.gtex.mappedgenes$probe_id, extract.gtex.mappedgenes$variant_chr, extract.gtex.mappedgenes$variant_pos), ]
  }
  
  #############
  # GWAS Catalog
  #############
  
  ### MQTL SNPS ###
  extract.gwascatalog.mqtlSNPs <- lapply(myInterpretList, function(y) y$gwasCatalog[["gwasCatalog_mqtlSNPs"]])
  extract.gwascatalog.mqtlSNPs <- do.call(rbind, extract.gwascatalog.mqtlSNPs)
  extract.gwascatalog.mqtlSNPs <- extract.gwascatalog.mqtlSNPs[complete.cases(extract.gwascatalog.mqtlSNPs$DATE.ADDED.TO.CATALOG),]
  extract.gwascatalog.mqtlSNPs$probe_id <- gsub("\\..*", "", row.names(extract.gwascatalog.mqtlSNPs))
  extract.gwascatalog.mqtlSNPs <- extract.gwascatalog.mqtlSNPs[c(length(extract.gwascatalog.mqtlSNPs), 1: (length(extract.gwascatalog.mqtlSNPs) - 1))]
  
  extract.gwascatalog.mqtlSNPs <- merge(final_name_location_df, extract.gwascatalog.mqtlSNPs, by.x = "probe_id", by.y = "probe_id", all.y = T)
  
  ### MAPPED GENES ###
  extract.gwascatalog.mappedgenes <- lapply(myInterpretList, function(y) y$gwasCatalog[["gwasCatalog_mappedGenes"]])
  extract.gwascatalog.mappedgenes <- do.call(rbind, extract.gwascatalog.mappedgenes)
  extract.gwascatalog.mappedgenes <- extract.gwascatalog.mappedgenes[!(is.na(extract.gwascatalog.mappedgenes$ensembl_gene_id)),]
  extract.gwascatalog.mappedgenes$probe_id <- gsub("\\..*", "", row.names(extract.gwascatalog.mappedgenes))
  extract.gwascatalog.mappedgenes <- extract.gwascatalog.mappedgenes[c(length(extract.gwascatalog.mappedgenes), 1: (length(extract.gwascatalog.mappedgenes) - 1))]
  extract.gwascatalog.mappedgenes$probeid <- NULL
  extract.gwascatalog.mappedgenes <- merge(final_name_location_df, extract.gwascatalog.mappedgenes, by.x = "probe_id", by.y = "probe_id", all.y = T)
  
  
  #############
  # GO Categories
  #############
  extract.gocategories <- lapply(myInterpretList, function(y) y[["goCategories"]])
  extract.gocategories <- do.call(rbind, extract.gocategories)
  extract.gocategories$probe_id <- gsub("\\..*", "", row.names(extract.gocategories))
  extract.gocategories <- merge(final_name_location_df, extract.gocategories, by.x = "probe_id", by.y = "probe_id", all.y = T)
  extract.gocategories <- extract.gocategories[complete.cases(extract.gocategories$goIDs),]
  
  
  #############
  # PGC Top Hits
  #############  
  
  ### MQTL SNPS ###
  extract.pgchits.mqtlsnps <- lapply(myInterpretList, function(y) y$pgc_top_hits[["pgc_top_hits_mqtlSNPs"]])
  extract.pgchits.mqtlsnps <- do.call(rbind, extract.pgchits.mqtlsnps)
  extract.pgchits.mqtlsnps$probe_id <- gsub("\\..*", "", row.names(extract.pgchits.mqtlsnps))
  extract.pgchits.mqtlsnps <- merge(final_name_location_df, extract.pgchits.mqtlsnps, by.x = "probe_id", by.y = "probe_id", all.y = T)
  
  ### MAPPED GENES ###
  extract.pgchits.mappedgenes <- lapply(myInterpretList, function(y) y$pgc_top_hits[["pgc_top_hits_mappedGenes"]])
  extract.pgchits.mappedgenes <- do.call(rbind, extract.pgchits.mappedgenes)
  extract.pgchits.mappedgenes <- extract.pgchits.mappedgenes[complete.cases(extract.pgchits.mappedgenes$mappedGenes),]
  extract.pgchits.mappedgenes$probe_id <- gsub("\\..*", "", row.names(extract.pgchits.mappedgenes))
  extract.pgchits.mappedgenes <- merge(final_name_location_df, extract.pgchits.mappedgenes, by.x = "probe_id", by.y = "probe_id", all.y = T)
  
  
  
  ##########
  # SUMMARY PAGE
  ##########
  # A true or false table indicating if results are present in Gtex, gwas catalog or pgc top hits for each queried gene
  
  master <- data.frame(final_name_location_df)
  master$has_mQTLs <- master$probe_id %in% extract.mQTL$DNA_M_PROBEID
  
  master$has_goCategories_mappedgenes <- master$probe_id %in% extract.gocategories$probe_id
  
  master$has_GTex_mqtlSNPs <- master$probe_id %in% extract.gtex.snps$probe_id
  master$has_GTex_mappedgenes <- master$probe_id %in% extract.gtex.mappedgenes$probe_id
  
  master$has_GwasCategories_mqtlSNPs <- master$probe_id %in% extract.gwascatalog.mqtlSNPs$probe_id
  master$has_GwasCategories_mappedgenes <- master$probe_id %in% extract.gwascatalog.mappedgenes$probe_id
  
  PGChit_mqtlSNPs <- ddply(extract.pgchits.mqtlsnps, .(probe_id), summarize, 
                           has_PGC_hit_mqtlSNPs = paste(unique(TopHitDisorders),collapse=";"))
  
  PGChit_mqtlSNPs$has_PGC_hit_mqtlSNPs <- ifelse(PGChit_mqtlSNPs$has_PGC_hit_mqtlSNPs == "NA", NA, PGChit_mqtlSNPs$has_PGC_hit_mqtlSNPs)
  PGChit_mqtlSNPs$has_PGC_hit_mqtlSNPs <- ifelse(is.na(PGChit_mqtlSNPs$has_PGC_hit_mqtlSNPs), FALSE, TRUE)
  
  master$has_PGC_hit_mappedgenes <- master$probe_id %in% extract.pgchits.mappedgenes$probe_id
  
  master <- merge(master, PGChit_mqtlSNPs, by = 'probe_id',all.x=T)
  master$has_PGC_hit_mqtlSNPs <- ifelse(is.na(master$has_PGC_hit_mqtlSNPs), FALSE, master$has_PGC_hit_mqtlSNPs)
  ########
  # Notification of no results
  ########
  results <- list(master, extract.mQTL, extract.gocategories, extract.gtex.snps, extract.gtex.mappedgenes,
                  extract.gwascatalog.mqtlSNPs, extract.gwascatalog.mappedgenes, extract.pgchits.mqtlsnps, extract.pgchits.mappedgenes)
  names(results) <- c("Summary", "mQTLs", "GOCategories_mappedgenes", "GTex_mqtlSNPs","GTex_mappedgenes", 
                      "GWASCatalog_mqtlSNPs","GWASCatalog_mappedgenes", "PGCHits_mqtlSNPs","PGCHits_mappedgenes")
  
  
  for (i in 1:length(results)) {
    
    if (nrow(results[[i]]) == 0){
      
      results[[i]] <- paste("There were no", names(results)[i], "results for your query", sep = " ")
      
    }
    
  }
  
  ########
  # Write in one excel file
  ########
  
  openxlsx::write.xlsx(results, file = results_fname_xlsx, keepNA = TRUE)
  
  
  return(results)
  
} ## END OF FUNCTION
geneticaList_temp <- function(queryList, input_genetic_data_type = c("gene", "snp", "probe")) {
  
  ## input genetic type = gene
  if (input_genetic_data_type == "gene") {
    
    interpret_headers <- c("name", "location", "goCategories", "gtex_cisEqtls", "gwasCatalog", "pgc_top_hits")
    toInterpret <- rep(list(NA),length(interpret_headers))
    names(toInterpret) = interpret_headers
    
    queryList <- sapply(queryList, function(x) {alias2Symbol(x, species = "Hs")}, USE.NAMES = FALSE)
    queryList <- queryList[lapply(queryList,length)>0]
    num_queryList <- length(queryList)
    
    myInterpretList = rep(list(toInterpret), num_queryList)
    names(myInterpretList) = queryList
    
    
    genemart = biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                                host = "grch37.ensembl.org", path = "/biomart/martservice", 
                                dataset = "hsapiens_gene_ensembl")
    
    genemart38 <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
    
    
    bm_genemart <- biomaRt::getBM(attributes = c("external_gene_name", "ensembl_gene_id", "description","chromosome_name", "start_position", "end_position", "strand"), 
                                  filters = "external_gene_name", 
                                  values = names(myInterpretList), 
                                  mart = genemart)
    
    bm_genemart$chromosome_name <- as.character(bm_genemart$chromosome_name)
    chr_check <- c(1:22,"X","Y")
    bm_genemart <- bm_genemart[bm_genemart$chromosome_name %in% chr_check,]
    
    
    bm_genemart38 <- biomaRt::getBM(attributes = c("external_gene_name", "chromosome_name", "start_position", "end_position"), 
                                    filters = "external_gene_name", 
                                    values = names(myInterpretList), 
                                    mart = genemart38)
    
    bm_genemart38$chromosome_name <- as.character(bm_genemart38$chromosome_name)
    bm_genemart38 <- bm_genemart38[bm_genemart38$chromosome_name %in% chr_check,]
    colnames(bm_genemart38) = c("external_gene_name", "chromosome_name_hg38", "start_hg38", "end_hg38")
    
    bm_genemart <- merge(bm_genemart, bm_genemart38, by = "external_gene_name", all = TRUE)
    no_interpretation_available <- list()
    for (g in queryList) {
      
      if (g %in% bm_genemart[, 1]) {
        
        bm_g <- bm_genemart[which(bm_genemart[, 1] == g), ]
        
        myInterpretList[[g]]$name <- bm_g[c("external_gene_name", "ensembl_gene_id", "description")]
        myInterpretList[[g]]$location <- bm_g[c("external_gene_name","chromosome_name", "start_position", "end_position", "strand", "chromosome_name_hg38", "start_hg38", "end_hg38")]
      }
      
      else {
        
        no_interpretation_available[[g]] <- g
        myInterpretList[[g]]$name <- data.frame(external_gene_name = g, ensembl_gene_id = NA, description  = NA)
        myInterpretList[[g]]$location <- data.frame(external_gene_name = g, chromosome_name = NA, start_position  = NA, end_position = NA, strand  = NA, chromosome_name_hg38 = NA, start_hg38 = NA, end_hg38 = NA)
      }
      
    }
    
    
    if (length(no_interpretation_available) != 0) {
      
      warning("Gene names missing: ", paste(no_interpretation_available, 
                                            collapse = ", "))
    }
  }
  
  ## input genetic type = snp
  else if (input_genetic_data_type == "snp") {
    
    interpret_headers <- c("name", "location", "goCategories", "gtex_cisEqtls", "gwasCatalog", "pgc_top_hits", "mqtls")
    toInterpret <- rep(list(NA),length(interpret_headers))
    names(toInterpret) = interpret_headers
    
    
    num_queryList <- length(queryList)
    
    myInterpretList = rep(list(toInterpret), num_queryList)
    names(myInterpretList) = queryList
    
    chr_check <- c(1:22,"X","Y")
    snpmart = biomaRt::useMart(biomart = "ENSEMBL_MART_SNP", 
                               host = "grch37.ensembl.org", path = "/biomart/martservice", 
                               dataset = "hsapiens_snp")
    
    genemart = biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                                host = "grch37.ensembl.org", path = "/biomart/martservice", 
                                dataset = "hsapiens_gene_ensembl")
    
    snpmart38 = biomaRt::useMart(biomart = "ENSEMBL_MART_SNP", dataset="hsapiens_snp")
    
    genemart38 <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl")
    
    bm_snpmart <- getBM(attributes = c("refsnp_id", "ensembl_gene_stable_id","allele"), 
                        filters = "snp_filter", 
                        values = names(myInterpretList), 
                        mart = snpmart)
    
    bm_snpmart38 <- getBM(attributes = c("refsnp_id", "ensembl_gene_stable_id","allele"), 
                          filters = "snp_filter", 
                          values = names(myInterpretList), 
                          mart = snpmart38)
    
    bm_genemart <- getBM(attributes = c("ensembl_gene_id", "external_gene_name","chromosome_name", "start_position", "end_position"), 
                         filters = "ensembl_gene_id", 
                         values = bm_snpmart$ensembl_gene_stable_id, 
                         mart = genemart)
    
    bm_genemart38 <- biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name", "start_position", "end_position"), 
                                    filters = "ensembl_gene_id", 
                                    values = bm_snpmart38$ensembl_gene_stable_id, 
                                    mart = genemart38)
    
    bm_genemart38$chromosome_name <- as.character(bm_genemart38$chromosome_name)
    bm_genemart38 <- bm_genemart38[bm_genemart38$chromosome_name %in% chr_check,]
    colnames(bm_genemart38) = c("ensembl_gene_id", "external_gene_name", "chromosome_name_hg38", "start_hg38", "end_hg38")
    
    colnames(bm_genemart) = c("ensembl_gene_id", "external_gene_name", "chromosome_name_hg19", "start_hg19", "end_hg19")
    
    bm_genemart <- merge(bm_genemart, bm_genemart38, by = c("external_gene_name","ensembl_gene_id") , all = TRUE)
    
    
    
    bm_snpmart <- merge(bm_snpmart, bm_genemart, by.x = 'ensembl_gene_stable_id', by.y = 'ensembl_gene_id', all.x = T)
    
    
    bm_snploc <- biomaRt::getBM(attributes = c("refsnp_id", "chr_name", "chrom_start", "chrom_end", "chrom_strand"), filters = "snp_filter", 
                                values = names(myInterpretList), mart = snpmart)
    
    bm_snploc38 <- biomaRt::getBM(attributes = c("refsnp_id", "chr_name", "chrom_start", "chrom_end"), filters = "snp_filter", 
                                  values = names(myInterpretList), mart = snpmart38)
    
    colnames(bm_snploc) <- c("refsnp_id", "chr_name_hg19", "chrom_start_hg19", "chrom_end_hg19", "chrom_strand")
    colnames(bm_snploc38) <- c("refsnp_id", "chr_name_hg38", "chrom_start_hg38", "chrom_end_hg38")
    
    bm_snploc <- merge(bm_snploc,bm_snploc38, by= "refsnp_id", all = TRUE)
    names(bm_snploc)[names(bm_snploc) == "refsnp_id"] <- "snp_id"
    
    bm <- ddply(bm_snpmart, .(refsnp_id), summarize,
                ensembl_gene_stable_id = paste(unique(ensembl_gene_stable_id),collapse=";"),
                snp_id = paste(unique(refsnp_id),collapse=";"),
                allele = paste(unique(allele),collapse=";"),
                gene_name = paste(external_gene_name,collapse=";"),
                gene_chr_hg19 = paste(unique(chromosome_name_hg19),collapse=";"),
                gene_start_position_hg19 = paste(unique(start_hg19),collapse=";"),
                gene_end_position_hg19 = paste(unique(end_hg19),collapse=";"),
                gene_chr_hg38 = paste(unique(chromosome_name_hg38),collapse=";"),
                gene_start_position_hg38 = paste(unique(start_hg38),collapse=";"),
                gene_end_position_hg38 = paste(unique(end_hg38),collapse=";"))
    bm$refsnp_id <- NULL
    bm <- bm[c(2:3,1,4:10)]
    
    no_interpretation_available <- list()
    for (g in queryList) {
      
      if (g %in% bm[, c("snp_id")]) {
        
        myInterpretList[[g]]$name <- bm[which(bm[, c("snp_id")] == g), ]
        myInterpretList[[g]]$location <- bm_snploc[which(bm_snploc[, 1] == g), ]
      }
      
      else {
        
        no_interpretation_available[[g]] <- g
        myInterpretList[[g]]$name <- data.frame(snp_id = g, allele = NA, gene_name  = NA, ensembl_gene_stable_id = NA, gene_chr37  = NA, gene_start_position37 = NA, gene_end_position37  = NA, gene_chr38  = NA, gene_start_position38 = NA, gene_end_position38  = NA)
        myInterpretList[[g]]$location <- data.frame(snp_id = g, chr_name = NA, chrom_start  = NA, chrom_end = NA, chrom_strand  = NA)
      }
      
    }
    
    
    if (length(no_interpretation_available) != 0) {
      
      warning("SNP names missing: ", paste(no_interpretation_available, 
                                           collapse = ", "))
    }
    
    
  }
  
  ## input genetic type = probe
  else if (input_genetic_data_type == "probe") {
    interpret_headers <- c("name", "location", "illumina.manifest", "mqtls", "goCategories", "gtex_cisEqtls", "gwasCatalog", "pgc_top_hits")
    toInterpret <- rep(list(NA),length(interpret_headers))
    names(toInterpret) = interpret_headers
    
    toInterpret$gtex_cisEqtls <- list(gtex_mqtlSNPs = NA, gtex_mappedGenes = NA)
    toInterpret$gwasCatalog <- list(gwasCatalog_mqtlSNPs = NA, gwasCatalog_mappedGenes = NA)
    toInterpret$pgc_top_hits <- list(pgc_top_hits_mqtlSNPs = NA, pgc_top_hits_mappedGenes = NA) 
    
    num_queryList <- length(queryList)
    myInterpretList = rep(list(toInterpret), num_queryList)
    names(myInterpretList) = queryList
    
    illum_manifest <- read.delim("mqtl_files/MethylationEPIC_v-1-0_B2_v3.txt", header = T, stringsAsFactors = F)
    illum_manifest <- illum_manifest[which(illum_manifest$IlmnID %in% names(myInterpretList)),]
    
    illum_manifest_colnames <- names(illum_manifest)
    
    empty_illumina_manifest <- data.frame(matrix(ncol = length(illum_manifest_colnames), nrow = 1))
    colnames(empty_illumina_manifest) <- illum_manifest_colnames
    
    for (i in 1:length(illum_manifest$UCSC_RefGene_Name)){
      if (grepl(";", illum_manifest$UCSC_RefGene_Name[[i]])) {
        probe_genes <- sort(unique(strsplit(illum_manifest$UCSC_RefGene_Name,";| ]")[[i]]))
        
        illum_manifest$UCSC_RefGene_Name[[i]] <- paste(probe_genes, collapse = ";")
      }
    }
    
    genemart = biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                                host = "grch37.ensembl.org", path = "/biomart/martservice", 
                                dataset = "hsapiens_gene_ensembl")
    
    
    bm_genemart <- biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name"), 
                                  filters = "external_gene_name", 
                                  values = unlist(unique(strsplit(illum_manifest$UCSC_RefGene_Name,";| ]"))), 
                                  mart = genemart)
    
    genemart38 <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
    
    
    bm_genemart38 <- biomaRt::getBM(attributes = c("external_gene_name", "chromosome_name", "start_position", "end_position"),
                                    filters = "external_gene_name", 
                                    values = unlist(unique(strsplit(illum_manifest$UCSC_RefGene_Name,";| ]"))), 
                                    mart = genemart38)
    
    colnames(bm_genemart38) = c("external_gene_name", "chromosome_name_hg38", "start_hg38", "end_hg38")
    bm_genemart <- merge(bm_genemart, bm_genemart38, by = "external_gene_name", all = TRUE)
    
    bm_genemart_updated <- data.frame()
    for (i in 1:length(illum_manifest$UCSC_RefGene_Name)) {
      val = illum_manifest$UCSC_RefGene_Name[i]
      if (grepl(";", val)) {
        
        split_i <- strsplit(val,";| ]")[[1]]
        split_i <- sort(split_i)
        bm_genemart_split_i <- bm_genemart[bm_genemart$external_gene_name %in% split_i,]
        bm_genemart_split_i <- bm_genemart_split_i[order(bm_genemart_split_i$external_gene_name),]
        bm_genemart_split_i <- t(data.frame(apply(rbind( bm_genemart_split_i), 2, paste, collapse=";")))
        bm_genemart_updated <- rbind(bm_genemart_updated,bm_genemart_split_i,make.row.names = FALSE)
      }
      
      else {
        bm_genemart_split_i = bm_genemart[bm_genemart$external_gene_name == val,]
        bm_genemart_updated <- rbind(bm_genemart_updated,bm_genemart_split_i,make.row.names = FALSE)
      }
      
    }
    
    illum_manifest <- merge(illum_manifest, bm_genemart_updated, by.x = "UCSC_RefGene_Name", by.y = "external_gene_name", all.x = T)
    
    no_interpretation_available <- list()
    for (g in queryList) {
      
      if (g %in% illum_manifest[, c("IlmnID")]) {
        
        name_df <- illum_manifest[which(illum_manifest[, "IlmnID"] == g), 
                                  c("IlmnID", "Genome_Build", "UCSC_RefGene_Name", "UCSC_RefGene_Accession", "UCSC_RefGene_Group", "ensembl_gene_id")]
        
        name_df_collapsed <- plyr::ddply(name_df, .(IlmnID), summarize,
                                         Genome_Build = paste(unique(Genome_Build),collapse=";"),
                                         UCSC_RefGene_Name = paste(unique(UCSC_RefGene_Name),collapse=";"),
                                         UCSC_RefGene_Accession = paste(unique(UCSC_RefGene_Accession),collapse=";"),
                                         UCSC_RefGene_Group = paste(unique(UCSC_RefGene_Group),collapse=";"),
                                         ensembl_gene_id = paste(unique(ensembl_gene_id),collapse=";"))
        
        myInterpretList[[g]]$name <- name_df_collapsed
        
        location_df <- illum_manifest[which(illum_manifest[, "IlmnID"] == g), 
                                      c("IlmnID", "CHR", "MAPINFO")]
        
        location_df_collapsed <- ddply(location_df, .(IlmnID), summarize,
                                       CHR = paste(unique(CHR),collapse=";"),
                                       POS = paste(unique(MAPINFO),collapse=";"))
        
        
        myInterpretList[[g]]$location <- location_df_collapsed
        
        myInterpretList[[g]]$illumina.manifest <- illum_manifest[which(illum_manifest[, "IlmnID"] == g),]
      }
      
      else {
        
        no_interpretation_available[[g]] <- g
        myInterpretList[[g]]$name <- data.frame(IlmnID = g, Genome_Build = NA, UCSC_RefGene_Name = NA, UCSC_RefGene_Accession  = NA, UCSC_RefGene_Group  = NA, ensembl_gene_id = NA)
        myInterpretList[[g]]$location <- data.frame(IlmnID = g, CHR  = NA, MAPINFO = NA)
        myInterpretList[[g]]$illumina.manifest <- empty_illumina_manifest
        myInterpretList[[g]]$illumina.manifest$IlmnID <- g
      }
      
    }
    
    
    
    if (length(no_interpretation_available) != 0) {
      
      warning("Probe names missing: ", paste(no_interpretation_available, 
                                             collapse = ", "))
    }
    
  }
  
  ## returns warning if input genetic type is != c("gene", "snp", "probe")
  else {
    
    warning(paste0("input genetic data type not supported:  ", input_genetic_data_type))
  }
  
  return(myInterpretList)
}
getGwasCatalog_temp <- function(myInterpretList, rangeval = 0, data_dir = "./", input_genetic_data_type = c("gene", "snp", "probe")) {
  
  
  gwas <- read.delim(url("http://www.ebi.ac.uk/gwas/api/search/downloads/full"), sep="\t", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
  r_val = rangeval
  gwas$UNIQID <- c(1:length(gwas$PUBMEDID))
  
  x <- colnames(gwas)
  x <- gsub(" ",".", x)
  x <- gsub("\\(",".", x)
  x <- gsub(")",".", x)
  x <- gsub("\\[",".", x)
  x <- gsub("]",".", x)
  x <- gsub("/",".", x)
  x <- gsub("-",".", x)
  x <- gsub("\\.\\.",".", x)
  colnames(gwas) <- x
  
  gwas_loc_table <- gwas %>% separate_rows(CHR_ID, CHR_POS)
  gwas_loc_table <- gwas_loc_table[gwas_loc_table$CHR_ID != 'x',]
  gwas_loc_table$CHR_POS <- ifelse(gwas_loc_table$CHR_POS == '', NA, gwas_loc_table$CHR_POS)
  gwas_loc_table$CHR_POS <- as.numeric(gwas_loc_table$CHR_POS)
  
  if (input_genetic_data_type == "gene") {
    
    for (i in names(myInterpretList)){
      i_bound <- paste0("\\b", i, "\\b")
      gwasCatalogresults <- gwas_loc_table[with(gwas_loc_table, grepl(i_bound, REPORTED.GENE.S.)|grepl(i_bound, MAPPED_GENE)),]
      gwasCatalogresults_range <- gwas_loc_table[which((gwas_loc_table$CHR_ID == as.character(myInterpretList[[i]]$location$chromosome_name_hg38)) & (gwas_loc_table$CHR_POS >= myInterpretList[[i]]$location$start_hg38 - rangeval) & (as.numeric(gwas_loc_table$CHR_POS) <= myInterpretList[[i]]$location$end_hg38 + rangeval)),]
      gwasCatalogresults_UNIQIDs <- unique(c(gwasCatalogresults$UNIQID, gwasCatalogresults_range$UNIQID))
      gwasCatalogresults_FINAL <- gwas[gwas$UNIQID %in% gwasCatalogresults_UNIQIDs,]
      rownames(gwasCatalogresults_FINAL) <- NULL
      
      names(gwasCatalogresults_FINAL)[names(gwasCatalogresults_FINAL) == 'CHR_ID'] <- "CHR_ID_hg38"
      names(gwasCatalogresults_FINAL)[names(gwasCatalogresults_FINAL) == 'CHR_POS'] <- "CHR_POS_hg38"
      myInterpretList[[i]]$gwasCatalog <- gwasCatalogresults_FINAL
    }
  }
  
  else if (input_genetic_data_type == "snp") {
    for (i in names(myInterpretList)){
      gwasCatalogresults <- gwas[with(gwas, grepl(paste0("^", i,"-.$"), STRONGEST.SNP.RISK.ALLELE)|grepl(paste0("^", i,"$"), SNPS)),]
      gwasCatalogresults_range <- gwas_loc_table[which(gwas_loc_table$CHR_ID ==  as.character(myInterpretList[[i]]$location$chr_name38) & (gwas_loc_table$CHR_POS >= myInterpretList[[i]]$location$chrom_start38 - rangeval & gwas_loc_table$CHR_POS <= myInterpretList[[i]]$location$chrom_end38 + rangeval)),]
      gwasCatalogresults_UNIQIDs <- unique(c(gwasCatalogresults$UNIQID, gwasCatalogresults_range$UNIQID))
      gwasCatalogresults_FINAL <- gwas[gwas$UNIQID %in% gwasCatalogresults_UNIQIDs,]
      rownames(gwasCatalogresults_FINAL) <- NULL
      myInterpretList[[i]]$gwasCatalog <- gwasCatalogresults_FINAL
    }
  }
  
  else if (input_genetic_data_type == "probe") {
    
    ensgs <- unlist(lapply(myInterpretList, function(x) {x$name$ensembl_gene_id} ))
    ensgs_clean <- ensgs[!is.na(ensgs)]
    ensgs_clean <- ensgs_clean[ensgs_clean != ""]
    ensgs_clean <- ensgs_clean[ensgs_clean != "NA"]
    
    df <- data.frame(ensgs = unlist(strsplit(ensgs_clean, split = ";")))
    df$probeid <- unlist(lapply(df$ensgs, function(i) {names(grep(i, ensgs_clean, value = T))}))
    row.names(df) <- NULL
    
    
    genemart = biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                                host = "grch37.ensembl.org", path = "/biomart/martservice", 
                                dataset = "hsapiens_gene_ensembl")
    
    genemart38 <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
    
    
    bm_genemart <- biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name","chromosome_name", "start_position", "end_position"), 
                                  filters = "ensembl_gene_id", 
                                  values = df$ensgs, 
                                  mart = genemart)
    
    bm_genemart38 <- biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name","chromosome_name", "start_position", "end_position"), 
                                    filters = "ensembl_gene_id", 
                                    values = df$ensgs, 
                                    mart = genemart38)
    
    df <- merge(df,bm_genemart, by.x = 'ensgs', by.y = 'ensembl_gene_id',all.x = T)
    names(df)[names(df) == "ensgs"] <- 'ensembl_gene_id'
    
    gwasCatalogresults_mappedGenes_cols <- c(names(df), names(gwas))
    
    
    for (i in names(myInterpretList)){
      
      
      df_subset <- df[df$probeid == i,]
      
      if(length(df_subset$probeid) == 1) {
        
        geneName = df_subset$external_gene_name
        geneName_bound <- paste0("\\b", geneName, "\\b")
        gwasCatalogresults <- gwas[with(gwas, grepl(geneName_bound, REPORTED.GENE.S.)|grepl(geneName_bound, MAPPED_GENE)),]
        #gwasCatalogresults_mappedGenes_range <- gwas[which(gwas$CHR_ID == df_subset$chromosome_name & (as.numeric(gwas$CHR_POS) >= df_subset$start_position) - rangeval & (as.numeric(gwas$CHR_POS) <= as.numeric(df_subset$end_position) + rangeval)),]
        #gwasCatalogresults_mappedGenes <- rbind(gwasCatalogresults_mappedGenes, gwasCatalogresults_mappedGenes_range)
        gwasCatalogresults_mappedGenes <- gwasCatalogresults_mappedGenes[!duplicated(gwasCatalogresults_mappedGenes),]
        rownames(gwasCatalogresults_mappedGenes) <- NULL
        
        df_subset_w_gwas <- cbind(df_subset, gwasCatalogresults_mappedGenes)
        
        
      } else if (length(df_subset$probeid) > 1) {
        
        df_subset_w_gwas = data.frame()
        
        for (j in df_subset$external_gene_name){
          
          df_subset2 <- df_subset[df_subset$external_gene_name == j,]
          j_bound <- paste0("\\b", j, "\\b")
          gwasCatalogresults <- gwas[with(gwas, grepl(j_bound, REPORTED.GENE.S.)|grepl(j_bound, MAPPED_GENE)),]
          #gwasCatalogresults_mappedGenes_range <- gwas[which(gwas$CHR_ID == df_subset2$chromosome_name & (as.numeric(gwas$CHR_POS) >= df_subset2$start_position) - rangeval & (as.numeric(gwas$CHR_POS) <= as.numeric(df_subset2$end_position) + rangeval)),]
          #gwasCatalogresults_mappedGenes <- rbind(gwasCatalogresults_mappedGenes, gwasCatalogresults_mappedGenes_range)
          gwasCatalogresults_mappedGenes <- gwasCatalogresults_mappedGenes[!duplicated(gwasCatalogresults_mappedGenes),]
          rownames(gwasCatalogresults_mappedGenes) <- NULL
          
          df_subset_temp <- cbind(df_subset2, gwasCatalogresults_mappedGenes)
          df_subset_w_gwas<- rbind(df_subset_w_gwas, df_subset_temp)
          
        } # for (i in df_subset$external_gene_name)
        
        
        
      } else {
        
        df_subset_w_gwas <- data.frame(matrix(ncol = length(gwasCatalogresults_mappedGenes_cols), nrow = 1))
        colnames(df_subset_w_gwas) <-  gwasCatalogresults_mappedGenes_cols
        
      } # else
      
      
      myInterpretList[[i]]$gwasCatalog$gwasCatalog_mappedGenes <- df_subset_w_gwas
      
      
      #### GWAS CATALOG RESULTS FOR MQTL SNPS
      
      mqtls = myInterpretList[[i]]$mqtls
      mqtls = mqtls[mqtls$SNP != "NA",]
      if (length(mqtls$SNP) != 0) {
        
        gwasCatalogresults_mqtlSNPs_ALL <- data.frame()
        for (j in 1:length(mqtls$SNP)) {
          gwasCatalogresults_mqtlSNPs <- gwas[with(gwas, grepl(paste0("^", mqtls$SNP[j],"-.$"), STRONGEST.SNP.RISK.ALLELE)|grepl(paste0("^", mqtls$SNP[j],"$"), SNPS)),]
          #gwasCatalogresults_mqtlSNPs_range <- gwas[which(gwas$CHR_ID == mqtls$SNP_CHR[[j]] & (as.numeric(gwas$CHR_POS) >= as.numeric(mqtls$SNP_BP[[j]]) - rangeval & as.numeric(gwas$CHR_POS) <= as.numeric(mqtls$SNP_BP[[j]]) + rangeval)),]
          #gwasCatalogresults_mqtlSNPs <- rbind(gwasCatalogresults_mqtlSNPs, gwasCatalogresults_mqtlSNPs_range)
          gwasCatalogresults_mqtlSNPs <- gwasCatalogresults_mqtlSNPs[!duplicated(gwasCatalogresults_mqtlSNPs),]
          rownames(gwasCatalogresults_mqtlSNPs) <- NULL
          gwasCatalogresults_mqtlSNPs_ALL <- rbind(gwasCatalogresults_mqtlSNPs_ALL,gwasCatalogresults_mqtlSNPs)
        }
        
        myInterpretList[[i]]$gwasCatalog$gwasCatalog_mqtlSNPs <- gwasCatalogresults_mqtlSNPs_ALL
        
        
      }
      
      else {
        
        myInterpretList[[i]]$gwasCatalog$gwasCatalog_mqtlSNPs <- data.frame(matrix(ncol = ncol(gwas), nrow = 1))
        colnames(myInterpretList[[i]]$gwasCatalog$gwasCatalog_mqtlSNPs) <-  names(gwas)
      }
      
    }# for (p in names(myInterpretList)){
  } ## elseif
  
  
  else {
    
    warning(paste0("input genetic data type not supported:  ", input_genetic_data_type))
    
  }
  
  return(myInterpretList)
} # DONE


dec14 <- c("SLC7A8",
           "MARK2",
           "PDLIM5",
           "VPS28",
           "ZNF706",
           "RBPMS",
           "FAM59A",
           "RBPMS-AS1")
myGeneList <- geneticaList_temp(queryList = dec14, input_genetic_data_type = "gene")
myGeneList <- getGtexCisEqtls(myInterpretList = myGeneList, data_dir = "./gtex_db_files/", input_genetic_data_type = "gene") 
myGeneList <- getGwasCatalog_temp(myInterpretList = myGeneList, data_dir = "./gwas_catalog_db_files/", input_genetic_data_type = "gene")
myGeneList <- getGoCategories(myInterpretList = myGeneList, input_genetic_data_type = "gene") 
myGeneList <- getPGCTopHits(myInterpretList  = myGeneList, data_dir = "/data/Priya/genetica/top_pgc_results_files/", input_genetic_data_type = "gene")
gene_results <- combine_genetica_genes(myInterpretList = myGeneList, results_fname_xlsx = " Dec15_SLC7A8_fullmodel_genetica_results_genes.xlsx")

geneticaList_temp <- function(queryList, input_genetic_data_type = c("gene", "snp", "probe")) {
  
  ## input genetic type = gene
  if (input_genetic_data_type == "gene") {
    
    interpret_headers <- c("name", "location", "goCategories", "gtex_cisEqtls", "gwasCatalog", "pgc_top_hits")
    toInterpret <- rep(list(NA),length(interpret_headers))
    names(toInterpret) = interpret_headers
    
    queryList <- sapply(queryList, function(x) {alias2Symbol(x, species = "Hs")}, USE.NAMES = FALSE)
    queryList <- queryList[lapply(queryList,length)>0]
    num_queryList <- length(queryList)
    
    myInterpretList = rep(list(toInterpret), num_queryList)
    names(myInterpretList) = queryList
    
    
    genemart = biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                                host = "grch37.ensembl.org", path = "/biomart/martservice", 
                                dataset = "hsapiens_gene_ensembl")
    
    genemart38 <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
    
    # bm_genemart <- biomaRt::getBM(attributes = c("external_gene_name", 
    #                                              "ensembl_gene_id", "description"), filters = "external_gene_name", 
    #                               values = names(myInterpretList), mart = genemart)
    # 
    # bm_locmart <- biomaRt::getBM(attributes = c("external_gene_name","chromosome_name", "start_position", "end_position", "strand"), filters = "external_gene_name", 
    #                              values = names(myInterpretList), mart = genemart)
    
    bm_genemart <- biomaRt::getBM(attributes = c("external_gene_name", "ensembl_gene_id", "description","chromosome_name", "start_position", "end_position", "strand"), 
                                  filters = "external_gene_name", 
                                  values = names(myInterpretList), 
                                  mart = genemart)
    
    bm_genemart$chromosome_name <- as.character(bm_genemart$chromosome_name)
    chr_check <- c(1:22,"X","Y")
    bm_genemart <- bm_genemart[bm_genemart$chromosome_name %in% chr_check,]
    
    
    bm_genemart38 <- biomaRt::getBM(attributes = c("external_gene_name", "chromosome_name", "start_position", "end_position"), 
                                    filters = "external_gene_name", 
                                    values = names(myInterpretList), 
                                    mart = genemart38)
    
    bm_genemart38$chromosome_name <- as.character(bm_genemart38$chromosome_name)
    bm_genemart38 <- bm_genemart38[bm_genemart38$chromosome_name %in% chr_check,]
    colnames(bm_genemart38) = c("external_gene_name", "chromosome_name_hg38", "start_hg38", "end_hg38")
    
    bm_genemart <- merge(bm_genemart, bm_genemart38, by = "external_gene_name", all = TRUE)
    no_interpretation_available <- list()
    for (g in queryList) {
      
      if (g %in% bm_genemart[, 1]) {
        
        bm_g <- bm_genemart[which(bm_genemart[, 1] == g), ]
        
        myInterpretList[[g]]$name <- bm_g[c("external_gene_name", "ensembl_gene_id", "description")]
        myInterpretList[[g]]$location <- bm_g[c("external_gene_name","chromosome_name", "start_position", "end_position", "strand", "chromosome_name_hg38", "start_hg38", "end_hg38")]
      }
      
      else {
        
        no_interpretation_available[[g]] <- g
        myInterpretList[[g]]$name <- data.frame(external_gene_name = g, ensembl_gene_id = NA, description  = NA)
        myInterpretList[[g]]$location <- data.frame(external_gene_name = g, chromosome_name = NA, start_position  = NA, end_position = NA, strand  = NA, chromosome_name_hg38 = NA, start_hg38 = NA, end_hg38 = NA)
      }
      
    }
    
    
    if (length(no_interpretation_available) != 0) {
      
      warning("Gene names missing: ", paste(no_interpretation_available, 
                                            collapse = ", "))
    }
  }
  
  ## input genetic type = snp
  else if (input_genetic_data_type == "snp") {
    
    interpret_headers <- c("name", "location", "goCategories", "gtex_cisEqtls", "gwasCatalog", "pgc_top_hits", "mqtls")
    toInterpret <- rep(list(NA),length(interpret_headers))
    names(toInterpret) = interpret_headers
    
    
    num_queryList <- length(queryList)
    
    myInterpretList = rep(list(toInterpret), num_queryList)
    names(myInterpretList) = queryList
    
    
    snpmart = biomaRt::useMart(biomart = "ENSEMBL_MART_SNP", 
                               host = "grch37.ensembl.org", path = "/biomart/martservice", 
                               dataset = "hsapiens_snp")
    
    genemart = biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                                host = "grch37.ensembl.org", path = "/biomart/martservice", 
                                dataset = "hsapiens_gene_ensembl")
    
    snpmart38 = biomaRt::useMart(biomart = "ENSEMBL_MART_SNP", dataset="hsapiens_snp")
    
    genemart38 <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl")
    
    bm_snpmart <- getBM(attributes = c("refsnp_id", "ensembl_gene_stable_id","allele"), 
                        filters = "snp_filter", 
                        values = names(myInterpretList), 
                        mart = snpmart)
    
    bm_snpmart38 <- getBM(attributes = c("refsnp_id", "ensembl_gene_stable_id","allele"), 
                        filters = "snp_filter", 
                        values = names(myInterpretList), 
                        mart = snpmart38)
    
    bm_genemart <- getBM(attributes = c("ensembl_gene_id", "external_gene_name","chromosome_name", "start_position", "end_position"), 
                         filters = "ensembl_gene_id", 
                         values = bm_snpmart$ensembl_gene_stable_id, 
                         mart = genemart)
    
    bm_genemart38 <- biomaRt::getBM(attributes = c("external_gene_name", "chromosome_name", "start_position", "end_position"), 
                                    filters = "external_gene_name", 
                                    values = names(myInterpretList), 
                                    mart = genemart38)
    
    bm_genemart38$chromosome_name <- as.character(bm_genemart38$chromosome_name)
    bm_genemart38 <- bm_genemart38[bm_genemart38$chromosome_name %in% chr_check,]
    colnames(bm_genemart38) = c("external_gene_name", "chromosome_name_hg38", "start_hg38", "end_hg38")
    
    bm_genemart <- merge(bm_genemart, bm_genemart38, by = "external_gene_name", all = TRUE)
    
    
    
    bm_snpmart <- merge(bm_snpmart, bm_genemart, by.x = 'ensembl_gene_stable_id', by.y = 'ensembl_gene_id', all.x = T)
    
    
    bm_snploc <- biomaRt::getBM(attributes = c("refsnp_id", "chr_name", "chrom_start", "chrom_end", "chrom_strand"), filters = "snp_filter", 
                                values = names(myInterpretList), mart = snpmart)
    
    bm_snploc38 <- biomaRt::getBM(attributes = c("refsnp_id", "chr_name", "chrom_start", "chrom_end"), filters = "snp_filter", 
                                values = names(myInterpretList), mart = snpmart38)
    
    colnames(bm_snploc38) <- c("refsnp_id", "chr_name38", "chrom_start38", "chrom_end38")
    
    bm_snploc <- merge(bm_snploc,bm_snploc38, by= "refsnp_id", all = TRUE)
    names(bm_snploc)[names(bm_snploc) == "refsnp_id"] <- "snp_id"
    
    bm <- ddply(bm_snpmart, .(refsnp_id), summarize,
                ensembl_gene_stable_id = paste(unique(ensembl_gene_stable_id),collapse=";"),
                snp_id = paste(unique(refsnp_id),collapse=";"),
                allele = paste(unique(allele),collapse=";"),
                gene_name = paste(external_gene_name,collapse=";"),
                gene_chr37 = paste(unique(chromosome_name),collapse=";"),
                gene_start_position37 = paste(unique(start_position),collapse=";"),
                gene_end_position37 = paste(unique(end_position),collapse=";"),
                gene_chr38 = paste(unique(chromosome_name_hg38),collapse=";"),
                gene_start_position38 = paste(unique(start_hg38),collapse=";"),
                gene_end_position38 = paste(unique(end_hg38),collapse=";"))
    
    bm <- bm[c(2:4,1,5:10)]
    
    no_interpretation_available <- list()
    for (g in queryList) {
      
      if (g %in% bm[, c("snp_id")]) {
        
        myInterpretList[[g]]$name <- bm[which(bm[, c("snp_id")] == g), ]
        myInterpretList[[g]]$location <- bm_snploc[which(bm_snploc[, 1] == g), ]
      }
      
      else {
        
        no_interpretation_available[[g]] <- g
        myInterpretList[[g]]$name <- data.frame(snp_id = g, allele = NA, gene_name  = NA, ensembl_gene_stable_id = NA, gene_chr37  = NA, gene_start_position37 = NA, gene_end_position37  = NA, gene_chr38  = NA, gene_start_position38 = NA, gene_end_position38  = NA)
        myInterpretList[[g]]$location <- data.frame(snp_id = g, chr_name = NA, chrom_start  = NA, chrom_end = NA, chrom_strand  = NA)
      }
      
    }
    
    
    if (length(no_interpretation_available) != 0) {
      
      warning("SNP names missing: ", paste(no_interpretation_available, 
                                           collapse = ", "))
    }
    
    
  }
  
  ## input genetic type = probe
  else if (input_genetic_data_type == "probe") {
    interpret_headers <- c("name", "location", "illumina.manifest", "mqtls", "goCategories", "gtex_cisEqtls", "gwasCatalog", "pgc_top_hits")
    toInterpret <- rep(list(NA),length(interpret_headers))
    names(toInterpret) = interpret_headers
    
    toInterpret$gtex_cisEqtls <- list(gtex_mqtlSNPs = NA, gtex_mappedGenes = NA)
    toInterpret$gwasCatalog <- list(gwasCatalog_mqtlSNPs = NA, gwasCatalog_mappedGenes = NA)
    toInterpret$pgc_top_hits <- list(pgc_top_hits_mqtlSNPs = NA, pgc_top_hits_mappedGenes = NA) 
    
    num_queryList <- length(queryList)
    myInterpretList = rep(list(toInterpret), num_queryList)
    names(myInterpretList) = queryList
    
    illum_manifest <- read.delim("mqtl_files/MethylationEPIC_v-1-0_B2_v3.txt", header = T, stringsAsFactors = F)
    illum_manifest <- illum_manifest[which(illum_manifest$IlmnID %in% names(myInterpretList)),]
    
    illum_manifest_colnames <- names(illum_manifest)
    
    empty_illumina_manifest <- data.frame(matrix(ncol = length(illum_manifest_colnames), nrow = 1))
    colnames(empty_illumina_manifest) <- illum_manifest_colnames
    
    for (i in 1:length(illum_manifest$UCSC_RefGene_Name)){
      if (grepl(";", illum_manifest$UCSC_RefGene_Name[[i]])) {
        probe_genes <- sort(unique(strsplit(illum_manifest$UCSC_RefGene_Name,";| ]")[[i]]))
        
        illum_manifest$UCSC_RefGene_Name[[i]] <- paste(probe_genes, collapse = ";")
      }
    }
    
    genemart = biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                                host = "grch37.ensembl.org", path = "/biomart/martservice", 
                                dataset = "hsapiens_gene_ensembl")
    
    
    bm_genemart <- biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name"), 
                                  filters = "external_gene_name", 
                                  values = unlist(unique(strsplit(illum_manifest$UCSC_RefGene_Name,";| ]"))), 
                                  mart = genemart)
    
    genemart38 <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
    
    
    bm_genemart38 <- biomaRt::getBM(attributes = c("external_gene_name", "chromosome_name", "start_position", "end_position"),
                                    filters = "external_gene_name", 
                                    values = unlist(unique(strsplit(illum_manifest$UCSC_RefGene_Name,";| ]"))), 
                                    mart = genemart38)
    
    colnames(bm_genemart38) = c("external_gene_name", "chromosome_name_hg38", "start_hg38", "end_hg38")
    bm_genemart <- merge(bm_genemart, bm_genemart38, by = "external_gene_name", all = TRUE)
    
    bm_genemart_updated <- data.frame()
    for (i in 1:length(illum_manifest$UCSC_RefGene_Name)) {
      val = illum_manifest$UCSC_RefGene_Name[i]
      if (grepl(";", val)) {
        
        split_i <- strsplit(val,";| ]")[[1]]
        split_i <- sort(split_i)
        bm_genemart_split_i <- bm_genemart[bm_genemart$external_gene_name %in% split_i,]
        bm_genemart_split_i <- bm_genemart_split_i[order(bm_genemart_split_i$external_gene_name),]
        bm_genemart_split_i <- t(data.frame(apply(rbind( bm_genemart_split_i[c('ensembl_gene_id','external_gene_name')]), 2, paste, collapse=";")))
        bm_genemart_updated <- rbind(bm_genemart_updated,bm_genemart_split_i,make.row.names = FALSE)
      }
      
      else {
        bm_genemart_split_i = bm_genemart[bm_genemart$external_gene_name == val,]
        bm_genemart_updated <- rbind(bm_genemart_updated,bm_genemart_split_i,make.row.names = FALSE)
      }
      
    }
    
    illum_manifest <- merge(illum_manifest, bm_genemart_updated, by.x = "UCSC_RefGene_Name", by.y = "external_gene_name", all.x = T)
    
    no_interpretation_available <- list()
    for (g in queryList) {
      
      if (g %in% illum_manifest[, c("IlmnID")]) {
        
        name_df <- illum_manifest[which(illum_manifest[, "IlmnID"] == g), 
                                  c("IlmnID", "Genome_Build", "UCSC_RefGene_Name", "UCSC_RefGene_Accession", "UCSC_RefGene_Group", "ensembl_gene_id")]
        
        name_df_collapsed <- plyr::ddply(name_df, .(IlmnID), summarize,
                                         Genome_Build = paste(unique(Genome_Build),collapse=";"),
                                         UCSC_RefGene_Name = paste(unique(UCSC_RefGene_Name),collapse=";"),
                                         UCSC_RefGene_Accession = paste(unique(UCSC_RefGene_Accession),collapse=";"),
                                         UCSC_RefGene_Group = paste(unique(UCSC_RefGene_Group),collapse=";"),
                                         ensembl_gene_id = paste(unique(ensembl_gene_id),collapse=";"))
        
        myInterpretList[[g]]$name <- name_df_collapsed
        
        location_df <- illum_manifest[which(illum_manifest[, "IlmnID"] == g), 
                                      c("IlmnID", "CHR", "MAPINFO")]
        
        location_df_collapsed <- ddply(location_df, .(IlmnID), summarize,
                                       CHR = paste(unique(CHR),collapse=";"),
                                       POS = paste(unique(MAPINFO),collapse=";"))
        
        
        myInterpretList[[g]]$location <- location_df_collapsed
        
        myInterpretList[[g]]$illumina.manifest <- illum_manifest[which(illum_manifest[, "IlmnID"] == g),]
      }
      
      else {
        
        no_interpretation_available[[g]] <- g
        myInterpretList[[g]]$name <- data.frame(IlmnID = g, Genome_Build = NA, UCSC_RefGene_Name = NA, UCSC_RefGene_Accession  = NA, UCSC_RefGene_Group  = NA, ensembl_gene_id = NA)
        myInterpretList[[g]]$location <- data.frame(IlmnID = g, CHR  = NA, MAPINFO = NA)
        myInterpretList[[g]]$illumina.manifest <- empty_illumina_manifest
        myInterpretList[[g]]$illumina.manifest$IlmnID <- g
      }
      
    }
    
    
    
    if (length(no_interpretation_available) != 0) {
      
      warning("Probe names missing: ", paste(no_interpretation_available, 
                                             collapse = ", "))
    }
    
  }
  
  ## returns warning if input genetic type is != c("gene", "snp", "probe")
  else {
    
    warning(paste0("input genetic data type not supported:  ", input_genetic_data_type))
  }
  
  return(myInterpretList)
}
getGwasCatalog_temp <- function(myInterpretList, rangeval = 0, data_dir = "./", input_genetic_data_type = c("gene", "snp", "probe")) {
  
  
  gwas <- read.delim(url("http://www.ebi.ac.uk/gwas/api/search/downloads/full"), sep="\t", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
  r_val = rangeval
  gwas$UNIQID <- c(1:length(gwas$PUBMEDID))
  
  x <- colnames(gwas)
  x <- gsub(" ",".", x)
  x <- gsub("\\(",".", x)
  x <- gsub(")",".", x)
  x <- gsub("\\[",".", x)
  x <- gsub("]",".", x)
  x <- gsub("/",".", x)
  x <- gsub("-",".", x)
  x <- gsub("\\.\\.",".", x)
  colnames(gwas) <- x
  
  gwas_loc_table <- gwas %>% separate_rows(CHR_ID, CHR_POS)
  gwas_loc_table <- gwas_loc_table[gwas_loc_table$CHR_ID != 'x',]
  gwas_loc_table$CHR_POS <- ifelse(gwas_loc_table$CHR_POS == '', NA, gwas_loc_table$CHR_POS)
  gwas_loc_table$CHR_POS <- as.numeric(gwas_loc_table$CHR_POS)
  
  if (input_genetic_data_type == "gene") {
    
    for (i in names(myInterpretList)){
      i_bound <- paste0("\\b", i, "\\b")
      gwasCatalogresults <- gwas_loc_table[with(gwas_loc_table, grepl(i_bound, REPORTED.GENE.S.)|grepl(i_bound, MAPPED_GENE)),]
      gwasCatalogresults_range <- gwas_loc_table[which((gwas_loc_table$CHR_ID == as.character(myInterpretList[[i]]$location$chromosome_name_hg38)) & (gwas_loc_table$CHR_POS >= myInterpretList[[i]]$location$start_hg38 - rangeval) & (as.numeric(gwas_loc_table$CHR_POS) <= myInterpretList[[i]]$location$end_hg38 + rangeval)),]
      gwasCatalogresults_UNIQIDs <- unique(c(gwasCatalogresults$UNIQID, gwasCatalogresults_range$UNIQID))
      gwasCatalogresults_FINAL <- gwas[gwas$UNIQID %in% gwasCatalogresults_UNIQIDs,]
      rownames(gwasCatalogresults_FINAL) <- NULL
      
      names(gwasCatalogresults_FINAL)[names(gwasCatalogresults_FINAL) == 'CHR_ID'] <- "CHR_ID_hg38"
      names(gwasCatalogresults_FINAL)[names(gwasCatalogresults_FINAL) == 'CHR_POS'] <- "CHR_POS_hg38"
      myInterpretList[[i]]$gwasCatalog <- gwasCatalogresults_FINAL
    }
  }
  
  else if (input_genetic_data_type == "snp") {
    for (i in names(myInterpretList)){
      gwasCatalogresults <- gwas[with(gwas, grepl(paste0("^", i,"-.$"), STRONGEST.SNP.RISK.ALLELE)|grepl(paste0("^", i,"$"), SNPS)),]
      gwasCatalogresults_range <- gwas_loc_table[which(gwas_loc_table$CHR_ID ==  as.character(myInterpretList[[i]]$location$chr_name38) & (gwas_loc_table$CHR_POS >= myInterpretList[[i]]$location$chrom_start38 - rangeval & gwas_loc_table$CHR_POS <= myInterpretList[[i]]$location$chrom_end38 + rangeval)),]
      gwasCatalogresults_UNIQIDs <- unique(c(gwasCatalogresults$UNIQID, gwasCatalogresults_range$UNIQID))
      gwasCatalogresults_FINAL <- gwas[gwas$UNIQID %in% gwasCatalogresults_UNIQIDs,]
      rownames(gwasCatalogresults_FINAL) <- NULL
      myInterpretList[[i]]$gwasCatalog <- gwasCatalogresults_FINAL
    }
  }
  
  else if (input_genetic_data_type == "probe") {
    
    ensgs <- unlist(lapply(myInterpretList, function(x) {x$name$ensembl_gene_id} ))
    ensgs_clean <- ensgs[!is.na(ensgs)]
    ensgs_clean <- ensgs_clean[ensgs_clean != ""]
    ensgs_clean <- ensgs_clean[ensgs_clean != "NA"]
    
    df <- data.frame(ensgs = unlist(strsplit(ensgs_clean, split = ";")))
    df$probeid <- unlist(lapply(df$ensgs, function(i) {names(grep(i, ensgs_clean, value = T))}))
    row.names(df) <- NULL
    
    
    genemart = biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                                host = "grch37.ensembl.org", path = "/biomart/martservice", 
                                dataset = "hsapiens_gene_ensembl")
    
    genemart38 <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
    
    
    bm_genemart <- biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name","chromosome_name", "start_position", "end_position"), 
                                  filters = "ensembl_gene_id", 
                                  values = df$ensgs, 
                                  mart = genemart)
    
    bm_genemart38 <- biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name","chromosome_name", "start_position", "end_position"), 
                                  filters = "ensembl_gene_id", 
                                  values = df$ensgs, 
                                  mart = genemart38)
    
    df <- merge(df,bm_genemart, by.x = 'ensgs', by.y = 'ensembl_gene_id',all.x = T)
    names(df)[names(df) == "ensgs"] <- 'ensembl_gene_id'
    
    gwasCatalogresults_mappedGenes_cols <- c(names(df), names(gwas))
    
    
    for (i in names(myInterpretList)){
      
      
      df_subset <- df[df$probeid == i,]
      
      if(length(df_subset$probeid) == 1) {
        
        geneName = df_subset$external_gene_name
        geneName_bound <- paste0("\\b", geneName, "\\b")
        gwasCatalogresults <- gwas[with(gwas, grepl(geneName_bound, REPORTED.GENE.S.)|grepl(geneName_bound, MAPPED_GENE)),]
        #gwasCatalogresults_mappedGenes_range <- gwas[which(gwas$CHR_ID == df_subset$chromosome_name & (as.numeric(gwas$CHR_POS) >= df_subset$start_position) - rangeval & (as.numeric(gwas$CHR_POS) <= as.numeric(df_subset$end_position) + rangeval)),]
        #gwasCatalogresults_mappedGenes <- rbind(gwasCatalogresults_mappedGenes, gwasCatalogresults_mappedGenes_range)
        gwasCatalogresults_mappedGenes <- gwasCatalogresults_mappedGenes[!duplicated(gwasCatalogresults_mappedGenes),]
        rownames(gwasCatalogresults_mappedGenes) <- NULL
        
        df_subset_w_gwas <- cbind(df_subset, gwasCatalogresults_mappedGenes)
        
        
      } else if (length(df_subset$probeid) > 1) {
        
        df_subset_w_gwas = data.frame()
        
        for (j in df_subset$external_gene_name){
          
          df_subset2 <- df_subset[df_subset$external_gene_name == j,]
          j_bound <- paste0("\\b", j, "\\b")
          gwasCatalogresults <- gwas[with(gwas, grepl(j_bound, REPORTED.GENE.S.)|grepl(j_bound, MAPPED_GENE)),]
          #gwasCatalogresults_mappedGenes_range <- gwas[which(gwas$CHR_ID == df_subset2$chromosome_name & (as.numeric(gwas$CHR_POS) >= df_subset2$start_position) - rangeval & (as.numeric(gwas$CHR_POS) <= as.numeric(df_subset2$end_position) + rangeval)),]
          #gwasCatalogresults_mappedGenes <- rbind(gwasCatalogresults_mappedGenes, gwasCatalogresults_mappedGenes_range)
          gwasCatalogresults_mappedGenes <- gwasCatalogresults_mappedGenes[!duplicated(gwasCatalogresults_mappedGenes),]
          rownames(gwasCatalogresults_mappedGenes) <- NULL
          
          df_subset_temp <- cbind(df_subset2, gwasCatalogresults_mappedGenes)
          df_subset_w_gwas<- rbind(df_subset_w_gwas, df_subset_temp)
          
        } # for (i in df_subset$external_gene_name)
        
        
        
      } else {
        
        df_subset_w_gwas <- data.frame(matrix(ncol = length(gwasCatalogresults_mappedGenes_cols), nrow = 1))
        colnames(df_subset_w_gwas) <-  gwasCatalogresults_mappedGenes_cols
        
      } # else
      
      
      myInterpretList[[i]]$gwasCatalog$gwasCatalog_mappedGenes <- df_subset_w_gwas
      
      
      #### GWAS CATALOG RESULTS FOR MQTL SNPS
      
      mqtls = myInterpretList[[i]]$mqtls
      mqtls = mqtls[mqtls$SNP != "NA",]
      if (length(mqtls$SNP) != 0) {
        
        gwasCatalogresults_mqtlSNPs_ALL <- data.frame()
        for (j in 1:length(mqtls$SNP)) {
          gwasCatalogresults_mqtlSNPs <- gwas[with(gwas, grepl(paste0("^", mqtls$SNP[j],"-.$"), STRONGEST.SNP.RISK.ALLELE)|grepl(paste0("^", mqtls$SNP[j],"$"), SNPS)),]
          #gwasCatalogresults_mqtlSNPs_range <- gwas[which(gwas$CHR_ID == mqtls$SNP_CHR[[j]] & (as.numeric(gwas$CHR_POS) >= as.numeric(mqtls$SNP_BP[[j]]) - rangeval & as.numeric(gwas$CHR_POS) <= as.numeric(mqtls$SNP_BP[[j]]) + rangeval)),]
          #gwasCatalogresults_mqtlSNPs <- rbind(gwasCatalogresults_mqtlSNPs, gwasCatalogresults_mqtlSNPs_range)
          gwasCatalogresults_mqtlSNPs <- gwasCatalogresults_mqtlSNPs[!duplicated(gwasCatalogresults_mqtlSNPs),]
          rownames(gwasCatalogresults_mqtlSNPs) <- NULL
          gwasCatalogresults_mqtlSNPs_ALL <- rbind(gwasCatalogresults_mqtlSNPs_ALL,gwasCatalogresults_mqtlSNPs)
        }
        
        myInterpretList[[i]]$gwasCatalog$gwasCatalog_mqtlSNPs <- gwasCatalogresults_mqtlSNPs_ALL
        
        
      }
      
      else {
        
        myInterpretList[[i]]$gwasCatalog$gwasCatalog_mqtlSNPs <- data.frame(matrix(ncol = ncol(gwas), nrow = 1))
        colnames(myInterpretList[[i]]$gwasCatalog$gwasCatalog_mqtlSNPs) <-  names(gwas)
      }
      
    }# for (p in names(myInterpretList)){
  } ## elseif
  
  
  else {
    
    warning(paste0("input genetic data type not supported:  ", input_genetic_data_type))
    
  }
  
  return(myInterpretList)
} # DONE



