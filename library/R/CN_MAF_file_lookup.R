## Assumes SIF columns derived from firehose annotations
firehose_CAPSEG_SIF = function( SIF, PP_CALLS_FN=NA, FORCE_CALL=FALSE, EXCLUDE_CALLED=FALSE, EXCLUDE_PASSED=FALSE, CONFIRM_RESULT_FILES=FALSE, MAF_COL=NA )
{
   if( (is.na(PP_CALLS_FN) || !file.exists(PP_CALLS_FN)) & any(c(FORCE_CALL, EXCLUDE_CALLED, EXCLUDE_PASSED)) ) 
   {
      msg = paste("Invalid PP_CALLS_FN supplied: ", PP_CALLS_FN, sep="")
      stop(msg)
   }

   if(EXCLUDE_CALLED | EXCLUDE_PASSED ) 
   {
      PP_calls = read.delim( PP_CALLS_FN, header=TRUE, check.names=FALSE, stringsAsFactors=FALSE, row.names=1 )
   }

#   if( !("indel.maf.fn" %in% colnames(SIF)) ) {
#      SIF = cbind(SIF, "indel.maf.fn"=NA)
#   }

   if( !("gender" %in% colnames(SIF)) ) {
      SIF = cbind(SIF, "gender"=NA)
   }

 ## needed by ABSOLUTE
   ACS_Rds_col = "allelic_capseg_rds"
   if( !is.null(SIF[[ACS_Rds_col]]) )
   {
      cols = c( "sample.name", "seg.dat.fn", "output.fn.base", "SSNV_skew", "gender", ACS_Rds_col )
   }
   else { cols =  c( "sample.name", "seg.dat.fn", "output.fn.base", "SSNV_skew", "gender" ) }

   print(cols)

   FH_maf_annots = c( "maf_file_capture_master_filter_removed", "union_maf_file_forcecalled", "maf_union_forcecalls", "maf_file_ffpeBias_capture", "maf_file_oxoG3_capture", "maf_file_capture", "maf_file", "maf_file_SSNV"  )
   FH_indel_annots = c( "strelka_passed_somatic_indel_maf_file_capture_pair", "maf_file_capture_strelka_forcecalled_newversion", "maf_file_capture_indel_forcecalled", "indel_maf_file_capture", "maf_file_indel" )

   FH_seg_annots = c("alleliccapseg_tsv", "allelic_scottseg_tsv", "HAPSEG_result", "seg.dat.fn")
   FH_skew_annots = c("AllelicCapseg_skew", "allelic_scottseg_skew")

   if( is.na(MAF_COL) )
   {
      maf_col = intersect( FH_maf_annots, colnames(SIF) ) [1]
   }
   else{ maf_col = MAF_COL }

   indel_col = intersect( FH_indel_annots, colnames(SIF) ) [1]
   seg_col = intersect( FH_seg_annots, colnames(SIF) ) [1]
   skew_col = intersect( FH_skew_annots, colnames(SIF) ) [1] 

   if( is.na(skew_col) )   ## ignored by ABS
   {
#      print("No allelic skew column found.  Filling in default value of 0.95")
      skew_col = FH_skew_annots[1]
      SIF[[skew_col]] = 0.95
   }
   if( is.na(seg_col) )  ##  ignored by ABS
   {
      seg_col = "alleliccapseg_tsv"
      SIF[[seg_col]] = NA 
   }


## indel and maf files are optional
   if( !is.null(SIF[[ACS_Rds_col]]) )
   {
      mat = data.frame( rownames(SIF), SIF[,seg_col], rownames(SIF), SIF[,skew_col], SIF[,"gender"], SIF[,ACS_Rds_col], stringsAsFactors=FALSE )
   }
   else
   {
      mat = data.frame( rownames(SIF), SIF[,seg_col], rownames(SIF), SIF[,skew_col], SIF[,"gender"], stringsAsFactors=FALSE )
   }

   colnames(mat) = cols


   if( !is.na(indel_col) )
   {
      mat[["indel.maf.fn"]] = SIF[,indel_col]
   }

   if( !is.na(maf_col) )
   {
      mat[["maf.fn"]] = SIF[,maf_col]
   }


   if( EXCLUDE_CALLED ) 
   {
      noncalled= mat[,"sample.name"][ !(mat[,"sample.name"] %in% rownames(PP_calls) ) ]
      N = nrow(mat)
      N_rem = N - length(noncalled)
      print( paste( "Removing ", N_rem, " of ", N, " samples with PP calls.", sep="" ))

      mat = mat[ mat[,"sample.name"] %in% noncalled, ]
   }

   if( EXCLUDE_PASSED ) 
   {
      int = intersect( mat[,"sample.name"], rownames(PP_calls) ) 
      keep = int[ !(PP_calls[int, "call status"] %in% c("called") )  ]

      mat = mat[ mat[,"sample.name"] %in% keep, ]
   }

   if( FORCE_CALL ) ## Lookup old PP modes and add t0 mat
   {
      pp = read.delim(PP_CALLS_FN, row.names=1, check.names=FALSE, stringsAsFactors=FALSE )

      SIDs = intersect( mat[,1], rownames(pp) )

      pp.dat =  data.frame(matrix( NA, ncol=2, nrow=nrow(mat) ))
      rownames(pp.dat) = mat[,1]
      pp.dat[SIDs,] = pp[SIDs, c("purity","tau")] 
      colnames(pp.dat) = c("force.alpha", "force.tau") 
      mat = cbind(mat, pp.dat)
   }


## check for missing samples... (essential cols)
#   nix = apply( is.na( mat[, c("seg.dat.fn", "maf.fn", "SSNV_skew")]), 1, any )

#   nix = is.na(mat[, ACS_Rds_col])
#   mat = mat[ !nix, ]
#   print( paste( "Removing ", sum(nix), " samples with missing seg.dat.fn", sep="") )
 
   if( CONFIRM_RESULT_FILES )
   {
      print("Confirming existence of result files needed for ABSOLUTE..." )
      nixmat = cbind( file.exists( mat[,"seg.dat.fn"] ), 
                      file.exists( mat[,"maf.fn"] ))
      nix = apply( nixmat, 1, sum ) < 2 
      print( paste( "Removing ", sum(nix), " samples with missing seg/MAF result files.", sep="") )
      cat("\n")
      if( sum(nix)>0 ) { print( paste( "Missing: ", mat[nix,1], sep="")) }
      mat = mat[ !nix, ]
   }


   return(mat)
}




get_batch_args_from_SIF_with_HAPSEG_map = function( SIF, snp6_map, PP_CALLS_FN=NA, FORCE_CALL=FALSE, EXCLUDE_CALLED=FALSE, EXCLUDE_PASSED=FALSE )
{
   sif  = firehose_CAPSEG_SIF( SIF, PP_CALLS_FN, FORCE_CALL=FORCE_CALL, EXCLUDE_CALLED=EXCLUDE_CALLED, EXCLUDE_PASSED=EXCLUDE_PASSED, CONFIRM_RESULT_FILES=FALSE )

 ## replace seg.dat.fn col with SNP6 hapseg result, instead of capseg result
#   snp6_map = read.delim(SNP6_HAPSEG_SIF_FN , stringsAsFactors=FALSE, check.names=FALSE )
#  "~scarter/Projects/pancan2/TCGA_SNP_ABS_tables/pancan2_snp6_mapped_uniq_patients.txt"

   m.ix = match( sif[,"sample.name"], snp6_map[,"pair_id"])

   sif[ !is.na(m.ix), "seg.dat.fn"] = snp6_map[m.ix[!is.na(m.ix)], "HAPSEG_result"]
   sif[  is.na(m.ix), "seg.dat.fn"] = NA

   return(sif)
}





## This is a temporary solution for CAPSEG results - will eventually be replaced by SIF_sample_query once CAPSEG segfile annotations are propagated into master SIF via Firehose.
simple_CAPSEG_SIF = function( SIF, PP_CALLS_FN=NA, FORCE_CALL=FALSE, EXCLUDE_CALLED=FALSE, EXCLUDE_PASSED=FALSE, default_SSNV_skew=1.0 )
{
   if( (is.na(PP_CALLS_FN) || !file.exists(PP_CALLS_FN)) & any(c(FORCE_CALL, EXCLUDE_CALLED, EXCLUDE_PASSED)) ) 
   {
      msg = paste("Invalid PP_CALLS_FN supplied: ", PP_CALLS_FN, sep="")
      stop(msg)
   }

   if(EXCLUDE_CALLED | EXCLUDE_PASSED ) 
   {
      PP_calls = read.delim( PP_CALLS_FN, header=TRUE, check.names=FALSE, stringsAsFactors=FALSE, row.names=1 )
   }

   if( !("indel.maf.fn" %in% colnames(SIF)) ) {
      SIF = cbind(SIF, "indel.maf.fn"=NA)
   }

   if( !("sample_SSNV_skew" %in% colnames(SIF)) ) {
      SIF = cbind(SIF, "sample_SSNV_skew"=default_SSNV_skew)
   }

 ## needed by ABSOLUTE
   cols = c( "sample.name", "seg.dat.fn", "output.fn.base", "maf.fn", "indel.maf.fn", "SSNV_skew", "gender" )

   mat = data.frame( rownames(SIF), SIF[,"capseg_segment_file"], rownames(SIF), SIF[,"oxogFiltered_oncotated_maf"], SIF[,"indel.maf.fn"], SIF[,"sample_SSNV_skew"], SIF[,"gender"], stringsAsFactors=FALSE )

   colnames(mat) = cols

   if( EXCLUDE_CALLED ) 
   {
      noncalled= mat[,"sample.name"][ !(mat[,"sample.name"] %in% rownames(PP_calls) ) ]
      N = nrow(mat)
      N_rem = N - length(noncalled)
      print( paste( "Removing ", N_rem, " of ", N, " samples with PP calls.", sep="" ))

      mat = mat[ mat[,"sample.name"] %in% noncalled, ]
   }

   if( EXCLUDE_PASSED ) 
   {
      int = intersect( mat[,"sample.name"], rownames(PP_calls) ) 
      keep = int[ !(PP_calls[int, "call status"] %in% c("called") )  ]

      mat = mat[ mat[,"sample.name"] %in% keep, ]
   }

   if( FORCE_CALL ) ## Lookup old PP modes and add t0 mat
   {
      pp = read.delim(PP_CALLS_FN, row.names=1, check.names=FALSE, stringsAsFactors=FALSE )

      SIDs = intersect( mat[,1], rownames(pp) )

      pp.dat =  data.frame(matrix( NA, ncol=2, nrow=nrow(mat) ))
      rownames(pp.dat) = mat[,1]
      pp.dat[SIDs,] = pp[SIDs, c("purity","tau")] 
      colnames(pp.dat) = c("force.alpha", "force.tau") 
      mat = cbind(mat, pp.dat)
   }

   return(mat)
}

HapSeg_res_lookup = function( array_info )
{
   print( paste( nrow(array_info), " arrays requested", sep="" ) )

   SIF = array_info

#   if( "maf_file_oxoG_capture" %in% colnames(SIF) ) { MAF_col = "maf_file_oxoG_capture" }
#   else { MAF_col = "maf_file_capture" }
   MAF_col = "maf_file_capture" 

   cols = c( "sample.name", "seg.dat.fn", "output.fn.base", "maf.fn", "gender" )
   mat = data.frame( SIF[,"sample_id"], SIF[,"snp6_hapseg_segdat"], SIF[,"snp6_array_name"], SIF[,MAF_col ],  SIF[,"BIRDSEED_GENDER"], stringsAsFactors=FALSE )
   colnames(mat) = cols

#oxogFiltered_oncotated_maf

 ## normalize gender calls.  Note that %in% evaluates to FALSE on comparison with NA, unlike '==', which returns NA.
   mat[mat[,"gender"] %in% "NoCall", "gender"] = NA
   mat[mat[,"gender"] %in% "M", "gender"] = "Male"
   mat[mat[,"gender"] %in% "F", "gender"] = "Female"
   mat[mat[,"gender"] %in% FALSE, "gender"] = "Female"  ## this can happen if the SIF has only females entered as "F"
 ## 

   array_info = mat
   missing_MAF = rep(FALSE, nrow(array_info))

   for( i in 1:nrow(array_info) )
   {
      array_dat_fn = array_info[i, "seg.dat.fn"]
      if( file.exists( array_dat_fn ) )
      {
         array_info[i, "seg.dat.fn"] = array_dat_fn
      }
      else { array_info[i, "seg.dat.fn"]=NA }

      if( is.na(array_info[i,"maf.fn"]) || !file.exists( array_info[i,"maf.fn"] ) )
      {
#         array_info[i, "maf.fn"]=NA 
         missing_MAF[i] = TRUE
      }
   }


   missing_ix = is.na(array_info[ , "seg.dat.fn"])  
   M = sum(!missing_ix)

   print( paste( M, " array results found", sep="" ) )
   if( M < nrow(array_info) )
   {
      print("NOT FOUND: ")
      print( array_info[ missing_ix, "output.fn.base" ] )
   }


   M = sum(!missing_MAF)
   print( paste( M, " MAF files found", sep=""))
   if( M < nrow(array_info) )
   { 
      print("Missing MAFs: ")
      print( array_info[ missing_MAF, "maf.fn" ] )
   }

   ###

   array_info = array_info[ !missing_ix, ]
   return(array_info)
}





SIF_sample_query = function( SIF_FN, PP_CALLS_FN, EXCLUDE_CALLED=FALSE, EXCLUDE_PASSED=FALSE, sample_type="", tumor_type="", GROUP="", sample_list_FN=NA )
{
   print("Reading in SIF" )
   sif = read.delim( SIF_FN, header=TRUE, check.names=FALSE, stringsAsFactors=FALSE, row.names=1 )

   if(EXCLUDE_CALLED | EXCLUDE_PASSED ) 
   {
      PP_calls = read.delim( PP_CALLS_FN, header=TRUE, check.names=FALSE, stringsAsFactors=FALSE, row.names=1 )
   }

   if( is.na(sample_list_FN) )
   {
      if( sample_type != "" )
      {
         sample_list = rownames( sif ) [ sif[,"sample_type"] %in% sample_type ] 
         sif = sif[ sample_list, ]
      }

      if( tumor_type != "" )
      {
         sample_list = rownames( sif ) [ sif[,"tumor_type"] %in% tumor_type ] 
         sif = sif[ sample_list, ]
      }

      if( GROUP != "" )
      {
         sample_list = rownames( sif ) [ sif[,"GROUP"] %in% GROUP ] 
         sif = sif[ sample_list, ]
      }
   }
   else
   {
      sample_list= rownames( read.delim(sample_list_FN, row.names=1) )

      missing = sample_list[ !(sample_list %in% rownames(sif)) ]
      if( length(missing) > 0 )
      {
         print( paste( length(missing), " requested arrays not in SIF: ", sep="") )
         print(missing) 
      }

      sn = intersect(sample_list, rownames(sif) )
      sif = sif[ sn,]

      norm.ix = sif[, "sample_type"] == "Normal"
      if( sum(norm.ix) > 0 )
      {
         print( paste( sum(norm.ix), " requested arrays marked Normal: ", sep="") )
         print(rownames(sif)[norm.ix]) 

         sif = sif[!norm.ix,] 
      }
      
   }

   sample_info = cbind( rep(NA,nrow(sif)), sif )
   sample_info[,1] = rownames(sif) 

   colnames(sample_info)[1] = "sample_id"
   rownames(sample_info) = sample_info[,1]

   if( EXCLUDE_CALLED ) 
   {
      noncalled= rownames(sample_info)[ !(sample_info[,"snp6_array_name"] %in% rownames(PP_calls) ) ]
      N = nrow(sample_info)
      N_rem = N - length(noncalled)

      print( paste( "Removing ", N_rem, " of ", N, " arrays with PP calls.", sep="" ))

      sample_info = sample_info[ noncalled, ]
   }

   if( EXCLUDE_PASSED ) 
   {
      int = intersect( rownames(sample_info), rownames(PP_calls) ) 

      keep = int[ !(PP_calls[int, "call status"] %in% c("called") )  ]
      sample_info = sample_info[ keep, ]
   }


   if( nrow(sample_info) == 0 ) { stop("no samples found!") }

   return( sample_info )
}



