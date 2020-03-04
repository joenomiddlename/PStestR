# PStestR
PStestR is an R package for testing for preferential sampling in spatio-temporal data. This includes discrete (i.e. areal) and continuous (i.e. point-referenced) data. The package is compatible with many popular R packages for handling spatial data including the spatstat, sp, and sf packages.

PStestR uses two functions: PSTestInit and PSTestRun. The first is a helper function for creating an object in the neccessary format for the test. The second function implements the test. 

To install PStestR, make sure the devtools package is loaded in R and run install_github('joenomiddlename/PStestR', dependencies=T, build_vignettes=T)

To use PStestR, read the tutorials found in the comprehensive vignette provided by typing vignette('PStestR') into R.
