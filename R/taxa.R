
kegg_taxa <- function(append = TRUE) {
    kegg_species_url <- "https://rest.kegg.jp/list/organism"
    # Read the files with caching
    kegg_species <- yread_tsv(kegg_species_url, params = list(header = FALSE))

    sn <- sub("\\s\\(.*$", "", kegg_species[,3])
    kegg.code = kegg_species[,2]

    if (append) {
        kegg_taxa <- readRDS(system.file("extdata/kegg_taxa.rds", package="clusterProfiler"))
        i <- which(!sn %in% kegg_taxa$kegg.name)
        sn <- sn[i]   
        kegg.code <- kegg.code[i]
    }

    tid <- character(length(sn))

    # pb <- progress::progress_bar$new(
    #         format = "[:bar] :current/:total (:percent)", 
    #         total = length(sn)
    #     )

    for (i in seq_along(sn)) {
        # pb$tick()
        if (tid[i] != "") next
        
        tid[i] <- tryCatch(
            getTaxInfo(sn[i], source = 'ensembl')$id,
            error = function(e) ""
        )
    }
    res <- data.frame(
        kegg.code = kegg.code,
        kegg.name = sn,
        kegg.taxa = tid)
    res <- res[res$kegg.taxa != "",]
    if (append) {
        kegg_taxa <- rbind(kegg_taxa, res)
    } else {
        kegg_taxa <- res
    }

    saveRDS(kegg_taxa, file = 'inst/extdata/kegg_taxa.rds')    
    invisible(kegg_taxa)
}

stringdb_version <- function(current = TRUE) {
    address <- "https://string-db.org"
    # Version info
    current_version <- read.table(url(paste(address, "/api/tsv/version", sep = "")), header = TRUE)
    if (current) return(current_version)
    available_version <- read.table(url(paste(address, "/api/tsv/available_api_versions", sep = "")), header = TRUE)
    return(available_version)
}


#' getTaxInfo
#'
#' Query taxonomy information from `stringdb` or `ensembl` web services
#' @param species scientific name of a species
#' @param source one of `stringdb` or `ensembl`
#' @importFrom rlang check_installed
#' @return a `data.frame` of query information
#' @author Guangchuang Yu
#' @export
getTaxInfo <- function(species, source = "stringdb") {
    source <- match.arg(source, c("stringdb", "ensembl"))
 
    check_installed('jsonlite', 'for `getTaxInfo()`.')
    
    if (source == "ensembl") {
        url <- paste0(
            "https://rest.ensembl.org/taxonomy/id/",
            gsub(" ", "%20", species),
            "?application/json"
        )

        res <- jsonlite::fromJSON(url)
        return(res)
    }

    stringdb_species <- getTaxInfo_stringdb()
    stringdb_species[stringdb_species$official_name_NCBI %in% species,]
}

getTaxInfo_stringdb <- function() {
    ver <- stringdb_version()
    stringdb_species_url <- paste0(
        "https://stringdb-static.org/download/species.v",
        ver$string_version, ".txt")
    stringdb_species <- yread_tsv(stringdb_species_url, params = list(header = TRUE))
    names(stringdb_species)[1] <- 'taxon_id'
    return(stringdb_species)
}

# getTaxInfo('Homo sapiens', source='ensembl')


#' getTaxID
#'
#' Convert species scientific name to taxonomic ID
#' @param species scientific name of a species
#' @return taxonomic ID
#' @author Guangchuang Yu
#' @export
getTaxID <- function(species) {
    if (inherits(species, "OrgDb")) {
        m <- AnnotationDbi::metadata(species)
        res <- m$value[m$name == "TAXID"]
        return(res)
    }

    ## load cached taxonomy ID information
    ## kegg.code kegg.name kegg.taxa
    ## e.g. hsa Homo Sapiens 9606
    ##
    ## use this function (will call `getTaxaInfo`) to prepare the cached data
    ## species list from https://rest.kegg.jp/list/organism
    ## extract species from the third column
    ## use second column as kegg.code 
    # tryCatch(utils::data(list='kegg_taxa', package='clusterProfiler'))
    # kegg_taxa <- get("kegg_taxa")
    kegg_taxa <- readRDS(system.file("extdata/kegg_taxa.rds", package="clusterProfiler"))
    
    i <- which(kegg_taxa[,1] == species | kegg_taxa[,2] == species)
    if (length(i) != 0)
    return(kegg_taxa[i,3])

    res <- getTaxInfo(species)
    return(res$id)
}



taxID2name <- function(taxID) {
    kegg_taxa <- readRDS(system.file("extdata/kegg_taxa.rds",
        package = "clusterProfiler"))
    kegg_taxa$kegg.name[kegg_taxa$kegg.taxa == taxID]
}

