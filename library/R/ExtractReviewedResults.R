## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

ExtractReviewedResults = function( called.segobj.list, analyst.id, out.dir.base, obj.name, verbose=FALSE) 
{
 ## agg MAF
  cat("Outputting aggregate MAF...")
  MAF_list_fn = file.path(out.dir.base, "reviewed", paste(obj.name, ".MAF_list.Rds", sep=""))
  if( !file.exists(MAF_list_fn) )
  {
     MAF_list = get_MAF_list_from_called_seglist_obj( called.segobj.list )
     saveRDS( MAF_list, file=MAF_list_fn )
  } else{ MAF_list = readRDS(MAF_list_fn) }

  fn = file.path(out.dir.base, "reviewed", paste(obj.name, ".aggregate_MAF.Rds", sep=""))
  if( !file.exists(fn) )
  {
     AGG_MAF = aggregate_sample_MAF_list( MAF_list )
     saveRDS( AGG_MAF, file=fn)
  } else{ AGG_MAF = readRDS(fn) }
  cat("done\n")
 ##

 ## significance of SSNVs
#  cat("Calculating significance of SSNVs...")
#  out.base = paste( file.path(out.dir.base, "reviewed", obj.name), "_", sep="")
#  ABS_class_mutsig( AGG_MAF, called.segobj.list, out.base )
#  cat("done\n")
 ##


 ## SEG_MAFs
  cat("Extracting SEG_MAF files...")
  seg.maf.dir = file.path(out.dir.base, "reviewed", "SEG_MAF")
  dir.create(seg.maf.dir, recursive=TRUE)
  write_called_seg_maf(called.segobj.list, pp.calls, seg.maf.dir)
  cat("done\n")
 ##

 ## Called summary plot
  pdf.fn = file.path(out.dir.base, "reviewed", 
                     paste(obj.name, ".called.ABSOLUTE.plots.pdf", sep=""))
  
#  PlotModes(called.segobj.list, chr.arms.dat, pdf.fn, n.print=1)

## This is useless because it does not include sample names in the plot!
if( FALSE )
{
  cat("Plotting called mode for matched samples")
  pdf( pdf.fn, 17.5, 18.5 )
  PlotModes_layout()
  for( i in 1:length(called.segobj.list)) 
  {        
     PlotModes(called.segobj.list[[i]], chr.arms.dat, n.print=1)
     cat(".")
  }
  dev.off()
  cat("done\n") 
}

  ## Called indv. RData files
   cat("Extracting RData called mode files for matched samples")
   indv.called.dir = file.path(out.dir.base, "reviewed", "samples")
   dir.create(indv.called.dir, recursive=TRUE)

   file.base = file.path(paste(names(called.segobj.list), ".ABSOLUTE.", analyst.id, 
                               ".called", sep = ""))
   called.files= file.path(indv.called.dir, paste(file.base, "RData", sep = "."))
   foreach (i=seq_along(called.files)) %dopar% {
      seg.obj = called.segobj.list[[i]]
      save(seg.obj, file=called.files[i])
      cat(".")
   }
   cat("done\n")
 
 ## print detailed SSNV plots for called solutions (1 plot per sample)
## These are intended to help debug SSNV on subclonal SCNA analysis
#   called_detailed_SSNV_plots( called.segobj.list, out.dir.base )


## Note this is inefficient as we are about to read in the files we just wrote out.
## Procuce matrix of corrected gene-level copy-numbers
   
   SCNA_thresholds = get_SCNA_thresholds( amp.CN.threshold = 7, H.amp.CN.threshold = 10 )

   out.dir.base = file.path( "ABSOLUTE_results", obj.name )
   transcript_GRs = get_GENCODE_transcript_GRs()
   gene_SCNA_calls = genotype_transcript_SCNAs_in_called_ABS_files( indv.called.dir, transcript_GRs, SCNA_thresholds, analyst_id = analyst.id )
   saveRDS( gene_SCNA_calls, file.path(out.dir.base, "reviewed", paste(obj.name, "_gene_SCNA_data.Rds", sep="")) )

   CN_dat = gene_SCNA_calls[["SCNA_event_dat"]][["amp.gene.data"]][,,"rescaled_total_cn"]
   write.table( round(CN_dat,2), file=file.path(out.dir.base, "reviewed", paste(obj.name, "_gene_corrected_CN.txt", sep="")), quote=FALSE, sep="\t" )

# write a merged-sample IGV file with log2 corrected copy-ratios
   write_IGV_segfile( seg.maf.dir, file.path( out.dir.base, "reviewed", paste( obj.name, "_rescaled_total_cn.IGV.seg.txt", sep="")) )

}

apply_review_and_extract = function( pp.review.fn=NA, pp.solution.num=NA, obj.name, analyst.id, pp.calls_ploidy_colname = "ploidy", ploidy_colname="genome mass",
 copy_num_type = "allelic", genome_build = "hg19", verbose=TRUE )
{
   if( copy_num_type == "allelic" )  {  set_allelic_funcs() }
   if( copy_num_type == "total" )  {  set_total_funcs() }
   if( !(copy_num_type %in% c("allelic", "total") ) ) { stop( "apply_review_and_extract: copy_num_type must be either 'allelic' or 'total'!") }

  if ( genome_build == "hg18") { chr.arms.dat.file = file.path(pkg_dir, "data", "hg18_ChrArmsDat.RData") }
  else if ( genome_build == "hg19") { chr.arms.dat.file = file.path(pkg_dir, "data", "hg19_ChrArmsDat.RData") }
  else if ( genome_build == "hg38") { chr.arms.dat.file = file.path(pkg_dir, "data", "hg38_ChrArmsDat.RData") }
  else if ( genome_build == "mm9" ) { chr.arms.dat.file = file.path(pkg_dir, "data", "mm9_ChrArmsDat.RData") }
  else {}
  print(paste("Loading", chr.arms.dat.file))
  load(chr.arms.dat.file)

  if (copy_num_type == "total") {
    set_total_funcs()
  } else if (copy_num_type == "allelic") {
    set_allelic_funcs()
  } else {
    stop("Unsupported copy number type: ", copy_num_type)
  }

   out.dir.base = file.path( "ABSOLUTE_results", obj.name )

   modesegs.fn = file.path(out.dir.base, paste0(obj.name, ".PP-modes.data.RData"))

   if ( file.exists(pp.review.fn) ) {
      called.segobj.list = run_PP_calls_liftover(pp.review.fn, analyst.id, modesegs.fn, out.dir.base, obj.name, chr.arms.dat, pp.calls_ploidy_colname, ploidy_colname, verbose=verbose )
   } else if (!isna(pp.solution.num)) {
     called.segobj.list = run_PP_calls_liftover_from_num(pp.solution.num, analyst.id, modesegs.fn, out.dir.base, obj.name, chr.arms.dat, pp.calls_ploidy_colname, ploidy_colname, verbose=verbose )
   } else {
     stop("pp.review.fn or pp.solution.num does not exist!")
   }

   if( length(called.segobj.list) > 0 )
   {
      ExtractReviewedResults( called.segobj.list, analyst.id , out.dir.base, obj.name, verbose=TRUE )
   }
   else{ stop("called.segobj.list has length 0!") }
}


called_detailed_SSNV_plots = function( called.segobj.list, out.dir.base )
{
  plot_dir = file.path(out.dir.base, "reviewed", "SSNV_detail")
  dir.create( plot_dir, recursive=TRUE)

  for ( i in 1:length(called.segobj.list) )
  {
     SID = names(called.segobj.list)[i]
     pdf.fn = file.path( plot_dir, paste(SID, ".SSNV.detail.plot.pdf", sep=""))

     seg.dat = called.segobj.list[[i]]
     mut.cn.dat <- seg.dat[["mut.cn.dat"]]
     modeled <- seg.dat[["mode.res"]][["modeled.muts"]][[1]]
     modeled.mut.dat <- cbind(mut.cn.dat, modeled)
#     SSNV_skew = modeled.mut.dat[1,"SSNV_skew"]
     SSNV_model = called.segobj.list[[i]][["mode.res"]][["mode_SSNV_models"]][[1]]

     pdf( pdf.fn, 10, 10 )

     detailed_SSNV_on_subclonal_HSCN_plot( SSNV_model, seg.dat, modeled.mut.dat, mode.ix=1, verbose=TRUE )
     
     dev.off()
  }
}


PP_liftover_plot = function( match.dat, PDF_OUT_FN )
{
   pdf( PDF_OUT_FN, 8, 4 )
   par(mfrow=c(1,2))
   par(bty='n')
   par(las=1)

   N = nrow(match.dat)
   PCH= rep(NA, N)
   col= rep(NA, N)

   nix = is.na(match.dat[,"mode.ix"])
   PCH=16
   col[nix] ="red"
   col[!nix]="black"

   CEX=0.8
   

   plot( match.dat[,"called_purity"], match.dat[,"best_purity"], xlim=c(0,1), ylim=c(0,1), main="", xlab="Called purity", ylab="Best match", pch=PCH, col=col )
   abline( a=0, b=1, col="grey", cex=CEX)


   MAXP=max(match.dat[,c("called_ploidy","best_ploidy")], na.rm=TRUE)
   plot( match.dat[,"called_ploidy"], match.dat[,"best_ploidy"], xlim=c(1,MAXP), ylim=c(1,MAXP), main="", xlab="Called ploidy", ylab="Best match",  pch=PCH, col=col, cex=CEX )
   abline( a=0, b=1, col="grey")

   dev.off()
}

