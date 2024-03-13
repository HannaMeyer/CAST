# This script creates the cookfarm dataset
#


cookfarm = readRDS("inst/extdata/Cookfarm.RDS")
save(cookfarm, file = "data/cookfarm.rda", compress = "xz")
