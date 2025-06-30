
# Setup Action
setwd( R'(C:\Users\James.Thorson\Desktop\Git\ecostate)' )
usethis::use_coverage()

# Run locally
setwd( R'(C:\Users\James.Thorson\Desktop\Git\ecostate)' )
library(covr)
#dir = tempdir()
#dir.create(dir)
covr::codecov() # install_path = dir ) #path = getwd(),
               #install_path = file.path(getwd(),"scratch", "tmp") )
package_coverage( path = getwd(),
                  type = "none",
                  install_path = tempdir() ) #covr:::temp_file("R_LIBS") )
