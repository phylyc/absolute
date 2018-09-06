## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

CreateReviewObject = function(obj.name, absolute.files, copy_num_type, plot.modes=TRUE, plot.mode.review=TRUE, num_solutions_plotted=NA, verbose=FALSE)
{
  
  if (copy_num_type == "total") {
    set_total_funcs()
  } else if (copy_num_type == "allelic") {
    set_allelic_funcs()
  } else {
    stop("Unsupported copy number type: ", copy_num_type)
  }
  

  ## read in processed SEG / MODES results and assemble
  segobj.list = vector(mode = "list", length = length(absolute.files))
  failed.list = vector(mode="list", length=length(absolute.files))
  so.ix = 1
  fa.ix = 1
  
  for (i in seq_along(absolute.files))  
  #  for (i in 1:10)
  {
    ## read in absolute rda
    seg.out.fn = absolute.files[i]
    if (!file.exists(seg.out.fn)) {
      if (verbose) {
        cat("\n")
        print(paste("sample #", i, " result not found", sep = ""))
      }
      next
    }

    ## provides seg.dat
    load(absolute.files[i])
    if (is.null(seg.dat$array.name)) {
      seg.dat$array.name = seg.dat$sample.name
    }
    SID = seg.dat$array.name
   
    if (is.na(seg.dat[["mode.res"]][["mode.flag"]])) {
      segobj.list[[so.ix]] = seg.dat
      names(segobj.list)[so.ix] = SID
      so.ix = so.ix + 1
      
      if (verbose) {
        cat(".")
      }
    } else {
      failed.list[[fa.ix]] = seg.dat
      names(failed.list)[fa.ix] = SID
      failed.list[[fa.ix]][["sample.name"]] = seg.dat[["sample.name"]]
      fa.ix = fa.ix + 1
      if (verbose) {
        cat("-")
      }
    }
  }
  
  segobj.list = segobj.list[c(1:(so.ix - 1))]
  failed.list = failed.list[c(1:(fa.ix - 1))]
  

  ## sort samples by diff in score between best and 2nd best mode
  mode.ent = rep(NA, length(segobj.list))
  names(mode.ent) = names(segobj.list)
  for (i in seq_along(segobj.list)) 
  {
    mtab = segobj.list[[i]][["mode.res"]][["mode.tab"]]
#
    if( any( is.na(mtab[, "combined_LL"]) ) ) 
    {
       cat("Bad sample found: ")
       cat(segobj.list[[i]][["sample.name"]] )
       cmd = paste( "echo ", segobj.list[[i]][["sample.name"]], " >> bad_samples", sep="")
       system(cmd)
       cat("\n")
       next
    }
#
#    ix = which.max(mtab[, "combined_LL"])
#    mode.ent[i] = mtab[ix, "entropy"]
#
     print(str(mtab))
     write.table(mtab, paste0("/cga/scarter/ncamarda/sds/ABSOLUTE_results/mtab", i, ".tsv"), quote = F, sep = "\t")
     ix = order(mtab[,"combined_LL"], decreasing=TRUE)[1:2]
     mode.ent[i] = mtab[ix[1], "combined_LL"] - mtab[ix[2], "combined_LL"]
  }

  mode.ent[is.na(mode.ent)] = Inf   ## handle case of only 1 solution

  samples = names(sort(mode.ent))
  segobj.list = segobj.list[samples]

  return( list(segobj.list=segobj.list, failed.list=failed.list) )
}
