.First.lib <- function(lib, pkg) {
    if (isTRUE(options('stringsAsFactors'))) {
        options(stringsAsFactors=FALSE);
        warning("The option stringsAsFactors has been set to FALSE.",
            call.=FALSE);
    }
}
