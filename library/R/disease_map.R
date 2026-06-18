## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

DetermineGroup <- function(primary.disease) {

  # data(diseaseMap)
  load(file.path(pkg_dir, "data", "diseaseMap.RData"))

  if(is.na(primary.disease)) { return(NA) }

  group <- try(get(primary.disease, disease_map), silent=TRUE)
  if (inherits(group, "try-error")) {
    ## It doesn't exist, just return the primary disease, it'll
    ## fall through
    print(paste("Disease type", primary.disease, "not in diseaseMap; falling through."))
    return(primary.disease)
  } else {
    print(paste("Disease type mapped:", primary.disease, "->", group))
    return(group)
  }
}
