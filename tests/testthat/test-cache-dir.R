## Attaching the package must not touch the disk.
##
## The CRAN packages depending on clusterProfiler load it during their own
## R CMD check. A cache directory created at attach time -- whatever its path
## resolves to -- is a write outside tempdir() from a check, which is what the
## CRAN complaint against clusterProfiler was about. So the default directory
## is resolved on first download (mechanism_cache_root(), get_ppi_network()),
## never in .onAttach.

test_that("attaching clusterProfiler writes nothing", {
    ud <- withr::local_tempdir()
    withr::local_options(list(clusterProfiler_cache_dir = NULL))
    withr::local_envvar(c(R_USER_DATA_DIR = ud))

    ## the attach hook is what runs when a CRAN dependent loads the package;
    ## call it directly rather than detaching the namespace the tests run in
    getNamespace("clusterProfiler")$.onAttach(ud, "clusterProfiler")

    expect_length(list.files(ud, recursive = TRUE, all.files = TRUE, no.. = TRUE), 0)
    expect_null(getOption("clusterProfiler_cache_dir"))
})

test_that("the default cache location is tools::R_user_dir()'s", {
    ud <- withr::local_tempdir()
    withr::local_envvar(c(R_USER_DATA_DIR = ud))
    withr::local_options(list(clusterProfiler_cache_dir = NULL))

    expect_equal(clusterProfiler:::mechanism_cache_root(),
                 tools::R_user_dir("clusterProfiler", "data"))
})

test_that("the cache directory is only created when a download needs it", {
    ## mechanism_cache_root() is called from code paths that read the cache too;
    ## the write belongs to .mechanism_ensure_dir(), next to the download.
    ud <- withr::local_tempdir()
    withr::local_envvar(c(R_USER_DATA_DIR = ud))
    dir <- file.path(ud, "new-subdir")

    expect_false(dir.exists(dir))
    clusterProfiler:::.mechanism_ensure_dir(dir)
    expect_true(dir.exists(dir))
    ## an existing directory is left alone, and a path that cannot be made
    ## (a file stands where the directory should be) is reported, not ignored
    expect_silent(clusterProfiler:::.mechanism_ensure_dir(dir))
    blocker <- withr::local_tempfile()
    writeLines("not a directory", blocker)
    expect_error(clusterProfiler:::.mechanism_ensure_dir(file.path(blocker, "cache")),
                 "cache directory")
})