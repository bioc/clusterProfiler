##' ORA analysis for Pathway Commons
##'
##' This function performs over-representation analysis using  Pathway Commons
##' @title enrichPC
##' @param gene a vector of genes (either hgnc symbols or uniprot IDs)
##' @param ... additional parameters, see also the parameters supported by the enricher() function
##' @return A \code{enrichResult} instance
##' @export
enrichPC <- function(gene, ...) {
    pcdata <- get_pc_data(output = 'gson')
    res <- enricher(gene, gson = pcdata, ...)

    if (is.null(res)) return(res)

    res@ontology <- pcdata@gsname
    res@organism <- pcdata@species
    res@keytype <-  keyType

    return(res)
}

##' GSEA analysis for  Pathway Commons
##'
##' This function performs GSEA using  Pathway Commons
##' @title gsePC
##' @param geneList a ranked gene list
##' @param ... additional parameters, see also the parameters supported by the GSEA() function
##' @importFrom rlang check_installed
##' @return A \code{gseaResult} instance
##' @export
gsePC <- function(geneList, source, keyType, ...) {
    pcdata <- get_pc_data(output = 'gson')
    res <- GSEA(geneList, gson = pcdata, ...)

    if (is.null(res)) return(res)

    res@ontology <- pcdata@gsname
    res@organism <- pcdata@species
    res@keytype <-  keyType

    return(res)
}

prepare_pc_data <- function(source, keyType) {
    pc2gene <- get_pc_data(source, keyType, output = 'data.frame')
    ##TERM2GENE
    pcid2gene <- pc2gene[, c("id", "gene")]
    ##TERM2NAME
    pcid2name <- unique(pc2gene[, c("id", "name")])

    list(PCID2GENE = pcid2gene,
        PCID2NAME = pcid2name)
}

get_pc_gmtfile <- function() {
    pc2 <- 'https://download.baderlab.org/PathwayCommons/PC2/'
    con <- readLines(pc2)
    pattern <- ".*>v(\\d+)/</a>.*"
    latest_version <- con[grep(pattern, con)] |>
        sub(pattern, "\\1", x = _) |>
        as.numeric() |>
        max()

    pcurl <- sprintf('%sv%s/', pc2, latest_version)

    x <- readLines(pcurl)
    y <- x[grep('\\.gmt.gz',x)]
    file <- sub(".*>(.*\\.gmt\\.gz)</a>.*", "\\1", y)
    sprintf("%s%s", pcurl, file)
}


read.gmt.pc_internal <- function(gmtfile) {
    # x <- readLines(gmtfile)
    
    check_installed('readr', 'for `read.gmt.pc_internal()`, which is an internal function.')
    
    x <- yread(gmtfile, readr::read_lines)
    
    y <- strsplit(x, "\t")
    id <- vapply(y, `[`, 1, FUN.VALUE = character(1))
    pcid <- sub(".*/", "", id)

    url <- sub(pcid[1], "", id[1]) # can be used to restored the url for web browse.

    nn <- vapply(y, `[`, 2, FUN.VALUE = character(1))
    names(y) <- sprintf("id: %s; %s", pcid, nn)

    y <- lapply(y, "[", -c(1:2))
  
    ont2gene <- stack(y)
    ont2gene <- ont2gene[, c("ind", "values")]
    colnames(ont2gene) <- c("term", "gene")
    return(ont2gene)
    # res <- list(ont2gene = ont2gene, pcid = pcid, url = url)
    # return(res)
}

##' Parse gmt file from Pathway Common
##'
##' This function parse gmt file downloaded from Pathway common
##' @title read.gmt.pc
##' @param gmtfile A gmt file   
##' @param output one of 'data.frame' or 'GSON'
##' @return A data.frame or A GSON object depends on the value of 'output'
##' @importFrom rlang .data
##' @importFrom tidyr separate
##' @export
read.gmt.pc <- function(gmtfile, output = "data.frame") {
    output <- match.arg(output, c("data.frame", "gson", "GSON"))

    pcdata <- read.gmt.pc_internal(gmtfile)
    x <- tidyr::separate(pcdata, .data$term, c("id", "name","datasource","organism","idtype"), "; ")
    x <- lapply(x, function(col) sub("\\w+:\\s*", "", col)) |> as.data.frame()
    if (output == "data.frame") {
        return(x)
    }
    

    gsid2gene <- data.frame(gsid=x$id, gene=x$gene)
    gsid2name <- unique(data.frame(gsid=x$id, name=x$name))
    organism <- taxID2name(x$organism[1])
    gson(gsid2gene = gsid2gene, 
        gsid2name = gsid2name, 
        gsname = "Pathway Commons", 
        species = organism)
}


get_pc_data <- function(output = "data.frame") {
    url <- get_pc_gmtfile()
    read.gmt.pc(url, output = output)
}
