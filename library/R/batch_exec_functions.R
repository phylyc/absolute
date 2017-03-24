#source("~scarter/CGA/R/broad_utils/scatter_jobs.R")
#source("~scarter/CGA/R/broad_utils/qq_pval.R")

#batch_exec_ABSOLUTE = function( ABSOLUTE_argv, obj.name, MAF_DIR, SIF_FN, sample_list_FN, MAF_SIF_FN, queue="hour", overwrite=FALSE, dry_run=FALSE )
batch_exec_ABSOLUTE = function( ABSOLUTE_argv, obj.name, var_bsub_argv, R_STUB_FN, queue="hour", overwrite=FALSE, dry_run=FALSE, num_solutions_plotted=NA, plot.mode.review=FALSE, wait=TRUE, summarize=TRUE, EXE_ENGINE="LSF", groupname="" )
{
   ABSOLUTE_argv = insert_default_args( ABSOLUTE_argv )

  ## setup job scatter control args
   control_argv = list()
	control_argv$OVERWRITE= overwrite
	control_argv$DRY_RUN = dry_run
	control_argv$wait = wait
	control_argv$QUEUE = queue
#	control_argv$R_STUB_FN = file.path(CGA_DIR, "broad_utils/ABSOLUTE_stub.R" )
	control_argv$R_STUB_FN = R_STUB_FN
        control_argv$BJOB = obj.name
        control_argv$EXE_ENGINE=EXE_ENGINE
        control_argv$groupname=groupname

        OUT_DIR_base = file.path( "ABSOLUTE_results", obj.name )
        control_argv$R_DIR = file.path( OUT_DIR_base, "R/" )

        RESULTS_DIR =  file.path( OUT_DIR_base  )
        ABSOLUTE_argv$results.dir = file.path( OUT_DIR_base, "results" )

   scatter_jobs( control_argv, ABSOLUTE_argv, var_bsub_argv )


 #### Step 2: collect results for summarization:
   if( summarize )
   {
      file.base = paste( var_bsub_argv[, "output.fn.base"], ".ABSOLUTE", sep = "")
      absolute.files= file.path( ABSOLUTE_argv$results.dir, paste(file.base, "RData", sep = "."))

      if( !file.exists(RESULTS_DIR)) {
         dir.create(RESULTS_DIR, recursive=TRUE, showWarnings=FALSE)
      }
      nm = file.path(RESULTS_DIR, paste(obj.name, ".PP-modes", sep = ""))
      modesegs.fn = paste(nm, ".data.RData", sep = "")
      failed.pdf.fn = paste(nm, "FAILED_plots.pdf", sep = ".")
      failed.tab.fn = paste(nm, "FAILED_tab.txt", sep = ".")
      call.tab.fn = file.path(RESULTS_DIR, paste(obj.name, "PP-calls_tab.txt", sep = "."))

      agg_res = CreateReviewObject( obj.name, absolute.files, ABSOLUTE_argv[["copy_num_type"]], plot.modes=TRUE, plot.mode.review=plot.mode.review, num_solutions_plotted=num_solutions_plotted, verbose=TRUE) 
#      save(agg_res, file = modesegs.fn)
      segobj.list = agg_res[["segobj.list"]]
      save(segobj.list, file = modesegs.fn)




#modesegs.fn = "ABSOLUTE_results/uncalled_test_df10___PR_TCGA_LUAD_PAIR_Capture_All_Pairs_QCPASS_v3__Illumina_WES/uncalled_test_df10___PR_TCGA_LUAD_PAIR_Capture_All_Pairs_QCPASS_v3__Illumina_WES.PP-modes.data.RData"
#load(modesegs.fn)
#pdf.fn = "test_dens_mode_summary_review_plot.pdf"
#num_solutions_plotted = NA

      PrintPpCallTable(segobj.list, call.tab.fn)
  
      if( plot.mode.review ) 
      {
         pdf.fn = paste(nm, ".dens.mode-review.plots.pdf", sep = "")
         pdf(pdf.fn, 17.5, 18.5 )
#         pdf(pdf.fn, 8.5, 11 )
         cat( paste("Plotting mode-review summary for ", length(segobj.list), " samples", sep="") )
         for( i in 1:length(segobj.list) )
         {
#            PlotModes_review_layout()  ## redo layout to make sure each sample starts at top of new page
#            dens_PlotMode_review_summary( segobj.list[[i]], n.print=4 )
#
            PlotModes_layout()
            PlotModes( segobj.list[[i]], n.print=4 )
            cat(".")
         }
         dev.off()
         cat("done.\n")
      }

      if (!is.null( agg_res[["failed.list"]][[1]])) 
      {
         try( PlotFailedSamples(agg_res[["failed.list"]], failed.pdf.fn) )
         PrintFailedTable(agg_res[["failed.list"]], failed.tab.fn)
      }

      rm(segobj.list);  gc() 
   }
#####
}
 

insert_default_args = function( ABSOLUTE_argv )
{
   default_abs_args = list( min.ploidy=1.1, max.ploidy=6,
                        max.as.seg.count=5000, max.non.clonal=0.99, max.neg.genome=0.05,
                        maf.fn=NA, indel.maf.fn=NA, min.mut.af = 0,
                        output.fn.base=NA, min_probes=1, max_sd=100, sigma.h=0.01, 
			SSNV_skew=1, filter_segs=TRUE, force.alpha=NA, force.tau=NA, allelic_capseg_rds=NA,
                        N_threads=1,
		        verbose=TRUE )

   for( i in 1:length(default_abs_args) )
   {
      if( !( names(default_abs_args)[i] %in% names(ABSOLUTE_argv) ) ) 
      {
         ABSOLUTE_argv[[  names(default_abs_args)[i]  ]] = default_abs_args[[i]] 
      }
   }

   return( ABSOLUTE_argv )
}

