.onLoad <- function(lib,pkg) {
  library.dynam("VBLPCM", pkg, lib)
  invisible()
}

.onAttach <- function(lib,pkg) {
  packageStartupMessage("\n\n")
  packageStartupMessage(paste("\n", "Variational Bayes Latent Position Cluster Model for networks.","\n"),
			paste("VBLPCM", "Version", "2.4.9", 
			      "Created on", "2023-03-22"), paste("\n Created and maintained by ", "Michael Salter-Townshend"), "\n")
  packageStartupMessage('For citation information type \'citation("VBLPCM")\'\n')
  packageStartupMessage('Type \'help(VBLPCM)\' to get started.\n')
  packageStartupMessage('Some worked examples are given by \'example(VBLPCM)\' \n')
} 

#.onUnload <- function(libpath){
#  library.dynam.unload("VBLPCM", libpath)
#  invisible()
#} 
