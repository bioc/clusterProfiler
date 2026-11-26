#' @importFrom utils packageDescription
.onAttach <- function(libname, pkgname) {
    # pkgVersion <- packageDescription(pkgname, fields="Version")
    # msg <- paste0(pkgname, " v", pkgVersion, "  ",
    #               "For help: https://yulab-smu.top/biomedical-knowledge-mining-book/", "\n\n")
    
    # The cache directory is deliberately NOT resolved here. user_dir() creates
    # what it returns, and attaching the package must not write to disk: the
    # CRAN packages depending on clusterProfiler load it during their own
    # check, and a directory created at attach time is precisely what the
    # out-of-tempdir write complaint was about. The option is read lazily by
    # mechanism_cache_root() and get_ppi_network() when a download first needs
    # a directory, and a value set in .Rprofile still wins.
    packageStartupMessage(yulab.utils::yulab_msg(pkgname))
}


