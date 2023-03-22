# VBLPCM 2.4.9

* Added a `NEWS.md` file to track changes to the package.

* fixed issue of "DLL requires the use of native symbols" by changing strings to symbols and using  .fixes="C_" in useDynLib _2023/03_    

# VBLPCM 2.4.8

* mirroring on github _2022/05_    

* updated configure script using autoupdate _2021/10_    

* removed use of externs in src .c files, now passing pointers to structs _2020/12_    

* using useDynLib now for compatibility with new R versions _2020/12_    

* had to change name of log_like_forces.c to ll_layout.c as was not compiling. Odd error resolved but not understood.  _2020/10_    

* changed use of deprecated plot.gofobject function from ergm to plot.gof and added cleanup script. _2018/10_    

* updated GSL requirements; added imports of 'sna' functions; added network to depends; shortened BIC example.  _2015/12_    

* removed inv_sigma02 term from omega2 update (thanks to Robin Gong).  _2014/01_    

* fixed bug in use of start function when using largest connected component only and edge covariates. _2014/01_    

* checked for leaks using AddressSanitizer, moved some Depends to Imports. _2013/10_    

* now if mclust fails Fruchterman-Reingold is (re)called until mclust can find clusters.  _2013/05_    

* changed setting of seed. A single set.seed is now used.  _2013/05_    

* plugged memory leak due do diagonal of sociomatrix and corrected case-control usage. _2013/04_    

* fixed bug in convergence check for V_xi_* and V_psi2_* (thanks for Triona Ryan). _2013/04_    

* fixed bug in random effects updates and tweaked plot function. _2013/02_    

* fixed bug in predict function for random effects models. _2013/02_    

* changed onAttach to onLoad in zzz.R _2013/01_    

* fixed bugs in handling of nodal covariates and cast same as sender, receiver or social.  _2012/08_    

* working startup message.  _2012/08_    

* fixed V_xi_e and V_psi2_e update bug (thanks to Sebastian Leist).  _2012/08_    

* changed to mclustBIC in vblpcmstart function for compatibility with mclust 4.  _2012/07_    
