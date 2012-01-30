.First.lib <- function(lib,pkg) {
    library.dynam("VBLPCM", pkg, lib)

    ehelp <- library(help="VBLPCM",lib.loc=lib,character.only=TRUE)[["info"]][[1]]
    packageStartupMessage("\n\n")
    packageStartupMessage(paste(substring(ehelp[4],first=16),"\n",
              "Version ",substring(ehelp[2],first=16),
              " created on ",
              substring(ehelp[3],first=16),".\n", sep=""))

    packageStartupMessage('For citation information type \'citation("VBLPCM")\'\n')
    packageStartupMessage('Type \'help(VBLPCM)\' to get started.\n')
    packageStartupMessage('Some worked examples are given by \'example(VBLPCM)\' \n')
} 


.onLoad <- function(lib, pkg) {library.dynam("VBLPCM", pkg, lib)}
