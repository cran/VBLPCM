.onAttach <- function(lib,pkg) {
    library.dynam("VBLPCM", pkg, lib)
    ehelp <- packageDescription("VBLPCM")
    packageStartupMessage("\n\n")
    packageStartupMessage(paste("\n",ehelp[3],"\n"),paste(ehelp[1], "Version", ehelp[4], 
        "Created on", ehelp[5]), paste("\n Created and maintained by ", ehelp[6]), "\n")
    packageStartupMessage('For citation information type \'citation("VBLPCM")\'\n')
    packageStartupMessage('Type \'help(VBLPCM)\' to get started.\n')
    packageStartupMessage('Some worked examples are given by \'example(VBLPCM)\' \n')
} 

