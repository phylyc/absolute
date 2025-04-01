## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

write_called_seg_maf = function(called_segobj_list, pp_calls, out_dir, verbose=FALSE) {
  dir.create(out_dir, showWarnings = FALSE)
  
  for (i in seq_along(called_segobj_list)) {
    called_segobj = called_segobj_list[[i]] 
    s_name = called_segobj$array.name
    
    abs_seg = GetAbsSegDat(called_segobj)
    WriteAbsSegtab(list(abs_seg), s_name, file.path(out_dir, paste(s_name, "segtab.txt", sep=".")))
    WriteIGVSegtab(called_segobj, abs_seg, s_name, out_dir)

    ## MAF
    if (!is.null(called_segobj$mode.res$modeled.muts)) {
      maf_out_fn = file.path(out_dir, paste(s_name, "ABS_MAF.txt", sep="."))
      WriteMAF(called_segobj, abs_seg, maf_out_fn)
    }
    if (verbose) {   
      cat(".")
    }
  }
}

WriteAbsSegtab <- function(seg, s_name, out_fn) {
  for (s in seq_along(seg)) {
    sample <- rep(s_name, nrow(seg[[s]]))
    s_tab <- cbind(sample, seg[[s]])
    
    ## colames only for 1st sample
    app <- s > 1   
    col <- s == 1
    
    write.table(s_tab, file=out_fn, col.names=col, append=app, row.names=FALSE,
                quote=FALSE, sep="\t")
  }   
}

WriteIGVSegtab <- function(segobj, seg, s_name, out_dir) {
  # Write IGV segtab:
  ploidy = segobj[["mode.res"]][["mode.tab"]][1, "genome mass"]
  seg[, "sample"] = s_name
  seg[, "Segment_Mean"] = log2(seg[, "rescaled_total_cn"] + 1e-4) - log2(ploidy)

  X.ix = (seg[,"Chromosome"] == "X")
  Y.ix = (seg[,"Chromosome"] == "Y")

  gender = segobj$gender
  if (!is.na(gender) && gender == "Male")
  {
    seg[X.ix, "Segment_Mean"] = seg[X.ix, "Segment_Mean"] + 1
    seg[Y.ix, "Segment_Mean"] = seg[Y.ix, "Segment_Mean"] + 1
  }
  # if( !is.na(gender) && gender == "Female" )
  # {
  #   seg[Y.ix, "Segment_Mean"] = seg[Y.ix, "Segment_Mean"] + 1
  # }
  write.table( seg[, c("sample", "Chromosome", "Start.bp", "End.bp", "Segment_Mean", "rescaled_total_cn")], file=file.path( out_dir, paste( s_name, "IGV.seg.txt", sep=".")), row.names=FALSE, sep="\t", quote=FALSE )
}

WriteMAF <- function(called_segobj, seg, out_fn) {
  modeled = called_segobj$mode.res$modeled.muts[[1]]
  mut_dat = cbind(called_segobj$mut.cn.dat, modeled)

  SSNV_ccf_dens = called_segobj[["mode.res"]][["SSNV.ccf.dens"]][1,,]
  maf = cbind(mut_dat, SSNV_ccf_dens)
  # revert to standard MAF
  cols <- colnames(maf)
  if ("ref" %in% cols) { colnames(maf)[which(cols %in% c("ref"))] <- c("t_ref_count") }
  if ("alt" %in% cols) { colnames(maf)[which(cols %in% c("alt"))] <- c("t_alt_count") }
  if ("dbSNP" %in% cols) { colnames(maf)[which(cols %in% c("dbSNP"))] <- c("dbSNP_Val_Status") }
  if ("sample" %in% cols) { colnames(maf)[which(cols %in% c("sample"))] <- c("Tumor_Sample_Barcode") }

  # Add local allelic CN to variants:
  # Convert to data.table and ensure maf has a start and end position (required by foverlaps)
  maf_dt <- as.data.table(maf)
  seg_dt <- as.data.table(seg[, c("Chromosome", "Start.bp", "End.bp", "rescaled.cn.a1", "rescaled.cn.a2")])

  # Ensure Chromosome is character in both tables
  maf_dt[, Chromosome := as.character(Chromosome)]
  seg_dt[, Chromosome := as.character(Chromosome)]

  # Set keys for fast binary search join
  setkey(maf_dt, Chromosome, Start_position, End_position)
  setkey(seg_dt, Chromosome, Start.bp, End.bp)

  # Perform the overlap join
  new_maf <- data.table::foverlaps(maf_dt, seg_dt, by.x = c("Chromosome", "Start_position", "End_position"), type = "within", nomatch = 0L)
  data.table::setnames(new_maf, old = "rescaled.cn.a1", new = "local_cn_a1")
  data.table::setnames(new_maf, old = "rescaled.cn.a2", new = "local_cn_a2")
  new_maf[, c("Start.bp", "End.bp") := NULL]
  setcolorder(new_maf, c(names(maf), c("local_cn_a1", "local_cn_a2")))

  write.table(file=out_fn, new_maf, row.names=FALSE, sep="\t", quote=FALSE)
}
