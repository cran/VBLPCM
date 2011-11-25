.First.lib <- function(lib,pkg) {
    library.dynam("VBLPCM", pkg, lib)

    ehelp <- library(help="VBLPCM",lib.loc=lib,character.only=TRUE)[["info"]][[1]]
    cat("\n\n")
    cat(paste(substring(ehelp[4],first=16),"\n",
              "Version ",substring(ehelp[2],first=16),
              " created on ",
              substring(ehelp[3],first=16),".\n", sep=""))

    cat('For citation information type \'citation("VBLPCM")\'\n')
    cat('Type \'help(VBLPCM)\' to get started.\n')
    cat('Some worked examples are given by \'example(VBLPCM)\' \n')
} 


.onLoad <- function(lib, pkg) {library.dynam("VBLPCM", pkg, lib)}
