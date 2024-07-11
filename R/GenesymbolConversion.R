#' @title Gene Naming Conversion Functions
#' @description These functions facilitate the conversion between human and mouse gene symbols and Ensembl IDs, and vice versa. They leverage both local and remote databases to provide fast and reliable gene identifier conversions, supporting a wide range of genetic studies.
#' @param human_genes A vector or matrix of human gene symbols to be converted to mouse gene symbols or Ensembl IDs. If a matrix is provided, the conversion is applied to the row names representing genes.
#' @param mouse_genes A vector or matrix of mouse gene symbols to be converted to human gene symbols or Ensembl IDs. If a matrix is provided, the conversion is applied to the row names representing genes.
#' @param Ensembl A vector or matrix of Ensembl IDs to be converted to gene symbols. If a matrix is provided, the conversion is applied to the row names representing genes.
#' @param Genesymbol A vector or matrix of gene symbols to be converted to Ensembl IDs. If a matrix is provided, the conversion is applied to the row names representing genes.
#' @param spe Specifies the species for which the conversion is being performed. Possible values are 'human' or 'mouse', depending on the function used.
#' @param mirror Optional parameter to specify an alternative BioMart mirror if the main site is inaccessible, when `local.mode` is set to FALSE. Possible mirrors include 'www', 'useast', and 'asia'. Default: NULL.
#' @param local.mode Indicates whether to use a local database for conversions. Using local databases is recommended for faster response times and increased reliability. Default: TRUE.
#' @param keep.seq Boolean flag to determine whether to maintain a one-to-one mapping in the output. This is useful when precise correspondence between input and output identifiers is required. Default: FALSE.
#' @param match Boolean flag to specify whether to return a matching table or a simple vector of results. When TRUE, the function returns a data frame showing how inputs match to outputs; when FALSE, it returns only the output identifiers. Default: TRUE.
#' @param keep.orig.id Only applicable for the EnsemblToGenesymbol function. This parameter is useful when converting Ensembl IDs to gene symbols, especially in cases where some Ensembl IDs do not have corresponding gene symbols. Setting this to TRUE allows the function to retain rows with original Ensembl IDs in the gene expression matrix, which is useful for downstream analyses such as PCA or clustering. Default: FALSE.
#' @return Depending on the function used and parameters set, the output can vary:
#'
#' - A vector of converted gene identifiers, if `match` is set to FALSE.
#'
#' - A named vector where names are the input identifiers and values are the converted identifiers, if `keep.seq` is set to TRUE.
#'
#' - A data.frame showing detailed matches between input and output identifiers, if `match` is set to TRUE.
#'
#' - A matrix with converted gene identifiers as row names, if the input was a matrix and `keep.seq` is TRUE. This output type is particularly useful for converting entire gene expression matrices for subsequent analyses.
#' @details These functions are essential tools for researchers working across genomic databases or conducting comparative studies between human and mouse. By utilizing both local and remote databases, these functions ensure that gene naming conversions are both accurate and adaptable to a variety of research needs.
#' @examples
#' # Load SeuratExtend and prepare a vector of human gene symbols
#' library(SeuratExtend)
#' human_genes <- VariableFeatures(pbmc)[1:6]
#' print(human_genes)
#'
#' # Default usage: Convert human gene symbols to mouse gene symbols and return a detailed match table
#' HumanToMouseGenesymbol(human_genes)
#'
#' # Simplified output without match details
#' MouseToHumanGenesymbol(human_genes, match = FALSE)
#'
#' # Ensure one-to-one correspondence with named vector output
#' GenesymbolToEnsembl(human_genes, keep.seq = TRUE)
#'
#' # Convert a gene expression matrix from human to mouse gene symbols
#' human_matr <- GetAssayData(pbmc)[human_genes, 1:4]
#' print(human_matr)
#' mouse_matr <- HumanToMouseGenesymbol(human_matr)
#' print(mouse_matr)
#'
#' # Convert mouse gene symbols to human gene symbols
#' mouse_genes <- c("Cd14", "Cd3d", "Cd79a")
#' MouseToHumanGenesymbol(mouse_genes, match = FALSE)
#'
#' # Convert human gene symbols to Ensembl IDs and back, ensuring one-to-one correspondence
#' ens_ids <- GenesymbolToEnsembl(human_genes, spe = "human", keep.seq = TRUE)
#' print(ens_ids)
#' back_to_symbols <- EnsemblToGenesymbol(ens_ids, spe = "human", keep.seq = TRUE)
#' print(back_to_symbols)
#'
#' # Fetch Ensembl IDs using an online BioMart database when local conversion is not sufficient
#' online_ens_ids <- GenesymbolToEnsembl(human_genes, spe = "human", local.mode = FALSE, keep.seq = TRUE)
#' print(online_ens_ids)
#' @rdname gene-naming-conversions
#' @export

HumanToMouseGenesymbol <- function(
    human_genes,
    mirror = NULL,
    local.mode = T,
    keep.seq = F,
    match = T) {
  result <- biomartConvertGeneral(
    item = human_genes,
    match.table = mouse_human_genesymbols,
    col.source = "HGNC.symbol",
    col.target = "MGI.symbol",
    biomart.par = c(
      dataset1 = "hsapiens_gene_ensembl",
      dataset2 = "mmusculus_gene_ensembl",
      attr1 = "hgnc_symbol",
      attr2 = "mgi_symbol"
    ),
    mirror = mirror,
    local.mode = local.mode,
    keep.seq = keep.seq,
    match = match
  )
  return(result)
}

#' @rdname gene-naming-conversions
#' @export

MouseToHumanGenesymbol <- function(
    mouse_genes,
    mirror = NULL,
    local.mode = T,
    keep.seq = F,
    match = T) {
  result <- biomartConvertGeneral(
    item = mouse_genes,
    match.table = mouse_human_genesymbols,
    col.source = "MGI.symbol",
    col.target = "HGNC.symbol",
    biomart.par = c(
      dataset1 = "mmusculus_gene_ensembl",
      dataset2 = "hsapiens_gene_ensembl",
      attr1 = "mgi_symbol",
      attr2 = "hgnc_symbol"
    ),
    mirror = mirror,
    local.mode = local.mode,
    keep.seq = keep.seq,
    match = match
  )
  return(result)
}

#' @rdname gene-naming-conversions
#' @export

EnsemblToGenesymbol <- function(
    Ensembl,
    spe = getOption("spe"),
    mirror = NULL,
    local.mode = T,
    keep.seq = F,
    match = T,
    keep.orig.id = F) {
  check_spe(spe)
  mappings <- list(
    human = list(
      match.table = human_ensembl_genesymbol,
      col.target = "hgnc_symbol",
      Dataset = "hsapiens_gene_ensembl",
      attr = "hgnc_symbol"
    ),
    mouse = list(
      match.table = mouse_ensembl_genesymbol,
      col.target = "mgi_symbol",
      Dataset = "mmusculus_gene_ensembl",
      attr = "mgi_symbol"
    )
  )
  selected_mapping <- mappings[[spe]]
  result <- biomartConvertGeneral(
    item = Ensembl,
    match.table = selected_mapping$match.table,
    col.source = "ensembl_gene_id",
    col.target = selected_mapping$col.target,
    biomart.par = c(
      dataset1 = selected_mapping$Dataset,
      attr1 = selected_mapping$attr
    ),
    mirror = mirror,
    local.mode = local.mode,
    keep.seq = keep.seq,
    match = match,
    keep.orig.id = keep.orig.id
  )
  return(result)
}

#' @rdname gene-naming-conversions
#' @export

GenesymbolToEnsembl <- function(
    Genesymbol,
    spe = getOption("spe"),
    mirror = NULL,
    local.mode = T,
    keep.seq = F,
    match = T) {
  check_spe(spe)
  mappings <- list(
    human = list(
      match.table = human_ensembl_genesymbol,
      col.source = "hgnc_symbol",
      Dataset = "hsapiens_gene_ensembl",
      attr = "hgnc_symbol"
    ),
    mouse = list(
      match.table = mouse_ensembl_genesymbol,
      col.source = "mgi_symbol",
      Dataset = "mmusculus_gene_ensembl",
      attr = "mgi_symbol"
    )
  )
  selected_mapping <- mappings[[spe]]
  result <- biomartConvertGeneral(
    item = Genesymbol,
    match.table = selected_mapping$match.table,
    col.source = selected_mapping$col.source,
    col.target = "ensembl_gene_id",
    biomart.par = c(
      dataset1 = selected_mapping$Dataset,
      attr1 = selected_mapping$attr
    ),
    mirror = mirror,
    local.mode = local.mode,
    keep.seq = keep.seq,
    match = match
  )
  return(result)
}

#' @title Convert UniProt IDs to Gene Symbols
#' @description Converts UniProt IDs to corresponding gene symbols for specified species using local databases. This function supports integration between proteomic and genomic data, enhancing the utility of proteomic identifiers in genomic analyses.
#' @param Uniprot A vector of UniProt IDs that are to be converted to gene symbols.
#' @param spe The species for which the conversion should be performed. The default is set to the user's global option (`getOption("spe")`), but can be explicitly specified as "human" or "mouse".
#' @return Returns a vector of gene symbols corresponding to the input UniProt IDs.
#' @details `UniprotToGenesymbol` utilizes a locally stored version of the UniProt database to convert UniProt IDs to gene symbols. This approach ensures fast and reliable conversions without the need for online database access. The local database includes comprehensive mappings from UniProt IDs to gene symbols for both human and mouse, facilitating seamless integration of proteomic and genomic data.
#' @examples
#' # Converting UniProt IDs to human gene symbols
#' UniprotToGenesymbol(c("Q8NF67", "Q9NPB9"), spe = "human")
#'
#' # Converting UniProt IDs to mouse gene symbols
#' UniprotToGenesymbol(c("Q9R1C8", "Q9QY84"), spe = "mouse")
#' @rdname UniprotToGenesymbol
#' @export

UniprotToGenesymbol <- function(Uniprot, spe = getOption("spe")) {
  check_spe(spe)
  if(any(!Uniprot %in% uniprot_db[[spe]]$Entry)) {
    Unknown <- setdiff(Uniprot, uniprot_db[[spe]]$Entry)
    message(paste0(length(Unknown), " item(s) not found in Uniprot ID list: \n",
                   paste0(head(Unknown, 10), collapse = "\n"),
                   ifelse(length(Unknown)>10, "\n...", "")))
  }
  return(uniprot_db[[spe]][Uniprot,"Gene names  (primary )"])
}

# internal ----------------------------------------------------------------

biomartConvertGeneral <- function(
    item,
    match.table,
    col.source,
    col.target,
    biomart.par, # c("dataset1", "dataset2", "attr1", "attr2")
    mirror = NULL,
    local.mode = T,
    keep.seq = F,
    match = T,
    keep.orig.id = F){
  if(is.vector(item)){
    if(local.mode){
      genelist <- match.table[match.table[[col.source]] %in% item, ]
    } else {
      item <- unique(item)
      import("biomaRt")
      message("Posible mirrors: 'www', 'useast', 'asia'.")
      if(!is.na(biomart.par['dataset2'])) {
        dataset1 <- useDataset(biomart.par[["dataset1"]], useEnsembl(biomart = "ensembl", mirror = mirror))
        dataset2 <- useDataset(biomart.par[["dataset2"]], useEnsembl(biomart = "ensembl", mirror = mirror))
        genelist <- getLDS(
          attributes = biomart.par[["attr1"]],
          filters = biomart.par[["attr1"]],
          values = item, mart = dataset1,
          attributesL = biomart.par[["attr2"]],
          martL = dataset2, uniqueRows=T)
      } else {
        mart <- useDataset(biomart.par[["dataset1"]], useEnsembl(biomart = "ensembl", mirror = mirror))
        genelist <- getBM(
          attributes = c(biomart.par[["attr1"]], "ensembl_gene_id"),
          filters = col.source,
          values = item,
          mart = mart)
      }
    }
    if(nrow(genelist) == 0) message("No homologous genes found")
    if(keep.seq) {
      genelist2 <- rep(NA, length(item))
      names(genelist2) <- item
      gl <- genelist
      gl <- gl[gl[[col.source]] != "" & gl[[col.target]] != "",]
      while(!is.null(gl)) {
        tmplist <- split(
          gl,
          f = duplicated(gl[[col.source]]) | duplicated(gl[[col.target]])
        )
        glf <- tmplist[["FALSE"]]
        glt <- tmplist[["TRUE"]]
        genelist2[glf[[col.source]]] <- glf[[col.target]]
        gl <- glt[!glt[[col.source]] %in% glf[[col.source]] & !glt[[col.target]] %in% glf[[col.target]], ]
      }
      if(keep.orig.id) genelist2[is.na(genelist2)] <- names(genelist2)[is.na(genelist2)]
      return(genelist2)
    }
    if(!match) genelist <- unique(genelist[[col.target]])
    return(genelist)
  }else{
    if(is.data.frame(item)) item <- as.matrix(item)
    rownames(item) <- biomartConvertGeneral(
      item = rownames(item),
      match.table = match.table,
      col.source = col.source,
      col.target = col.target,
      biomart.par = biomart.par, # c("dataset1", "dataset2", "attr1", "attr2")
      mirror = mirror,
      local.mode = local.mode,
      keep.seq = T,
      keep.orig.id = keep.orig.id)
    item <- item[!is.na(rownames(item)),]
    return(item)
  }
}
