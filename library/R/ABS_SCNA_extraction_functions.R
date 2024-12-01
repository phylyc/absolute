get_old_SCNA_thresholds = function()
{
 ## Used in Brastianos et al. Can Disc 2015  (brain met paper 1)
  ## Default threshold parameters for calling SCNAs
   SCNA_thresholds = list()
   SCNA_thresholds[["amp_thresh_slope"]] = -1/5
   SCNA_thresholds[["amp_delta_k"]] = 0.1
   SCNA_thresholds[["amp_thresh_foc"]] = 0.98
   SCNA_thresholds[["amp_thresh_log2_rCN"]] = log( 5, 2 )

   SCNA_thresholds[["del_thresh_rCN"]] = 0.25  
   SCNA_thresholds[["del_thresh_foc"]] = 0.995

   SCNA_thresholds[["H.amp_thresh_slope"]] = -1/5
   SCNA_thresholds[["H.amp_delta_k"]] = 0.1
   SCNA_thresholds[["H.amp_thresh_foc"]] = 0.98
   SCNA_thresholds[["H.amp_thresh_log2_rCN"]] = log(7, 2)

   return( SCNA_thresholds )
}

get_SCNA_thresholds = function( amp.CN.threshold = 7, H.amp.CN.threshold = 10 )
{
## Threshold parameters for calling SCNAs
   SCNA_thresholds = list()
## This works OK - hom dels in CDKN2A can be up to 11MB! (0.5% genome)
   SCNA_thresholds[["del_thresh_rCN"]] = 0.25  
   SCNA_thresholds[["del_thresh_foc"]] = 0.995

   SCNA_thresholds[["amp_thresh_slope"]] = -1/5
   SCNA_thresholds[["amp_delta_k"]] = 0.1
   SCNA_thresholds[["amp_thresh_foc"]] = 0.98
   SCNA_thresholds[["amp_thresh_log2_rCN"]] = log(amp.CN.threshold, 2)

   SCNA_thresholds[["H.amp_thresh_slope"]] = -1/5
   SCNA_thresholds[["H.amp_delta_k"]] = 0.1
   SCNA_thresholds[["H.amp_thresh_foc"]] = 0.98
   SCNA_thresholds[["H.amp_thresh_log2_rCN"]] = log( H.amp.CN.threshold, 2 )

   return(SCNA_thresholds)
}



build_gene_GR_data = function()
{
   #genes_FN = "refGene.hg19.20100825.sorted.txt"
   #genes_FN = "/xchip/cga/reference/annotation/db/ucsc/hg19/gene_table.txt"
#   genes_FN = "/xchip/gistic/variables/hg19/rg_20120227_dump.txt"
#  gene_dat = read.delim(genes_FN, check.names=FALSE, stringsAsFactors=FALSE )

   data(refgene)
   gene_dat = refgene 

   gene_dat[,"chr"] =  gsub("chr", "", gene_dat[,"chr"] )
   gene_footprints = GRanges( gene_dat[,"chr"], IRanges(gene_dat[,"start"], gene_dat[,"end"]) )
#   colnames(gene_dat)[3] = "name"   ## symb -> name
   names(gene_footprints) = gene_dat[,3]
   
   #name	longname	chr	strand	tx_start	tx_end
   #gene_footprints = GRanges( gene_dat[,"chr"], IRanges(gene_dat[,"tx_start"], gene_dat[,"tx_end"]) )
   #gene_footprints = GRanges( gene_dat[,"chr"], IRanges(gene_dat[,"gene_start"], gene_dat[,"gene_end"]) )

  return( gene_footprints )
}

load_GISTIC_peak_GR_data = function( GISTIC_peaks_fn, reg_tab=NA )
{ 
   GISTIC_peaks = read.delim( GISTIC_peaks_fn, check.names=FALSE, stringsAsFactors=FALSE, header=0, row.names=1, nrows=5 )
   GISTIC_peaks = GISTIC_peaks[,-ncol(GISTIC_peaks)]   ## extra col
   peak_strs = as.character( GISTIC_peaks["wide peak boundaries",] )
   
   peak_mat =  matrix( unlist(strsplit( peak_strs, ":|-")), ncol=3, byrow=TRUE )
   peak_mat[,1] = gsub( "chr", "", peak_mat[,1] )
   
   peak_GR = GRanges( peak_mat[,1], IRanges(as.integer(peak_mat[,2]), as.integer(peak_mat[,3])) )
   names(peak_GR) =  GISTIC_peaks["wide peak boundaries",]

   peak_GR$cytoband <- as.character(GISTIC_peaks["cytoband",] )

# try to rename regs
   if( !is.na(reg_tab))
   {
      r.ix = match( names(peak_GR), reg_tab[,"Peak region"] )
      if( any(is.na(r.ix)) ) { stop("unmatched reg name!") }

      reg_names = reg_tab[r.ix, "Peak Name"]
      names(peak_GR) = reg_names
   }

   return(peak_GR)
} 


compute_focality_score = function( segtab, mode )
{
   if( !(mode %in% c("amp", "del") ) ) { stop() }

   cn = segtab[,"rescaled_total_cn"]
   WD = 0.1
   br = seq( min(cn), max(cn)+WD, by=WD)

   bin = rep(NA, length(cn)) 
   for (i in 1:length(cn) )
   {
      bin[i] <- max(which(br <= cn[i]))
   }

   levels = sort(unique(bin))

   if( mode == "del" ) { levels = rev(levels) }

   genome_frac = rep(NA, length(levels)) 
   for( i in 1:length(levels) )
   {
      genome_frac[i] = sum(segtab[bin==levels[i],"W"])
   }  

   cum_genome_frac = cumsum(c(0,genome_frac) )
   focality = rep(NA, nrow(segtab))
   for( i in 1:length(levels) )
   {
     focality[bin==levels[i]] = cum_genome_frac[i] 
   }  

   return(focality)
}


## Calls all segments
## used in Brain met 2015 paper
original_call_genome_wide_ABSOLUTE_SCNAs = function( ABS.dat, SCNA_thresholds )
{
   segtab = AllelicGetAbsSegDat(ABS.dat)

   seg_amp_focality = compute_focality_score(segtab, "amp" )
   seg_del_focality = compute_focality_score(segtab, "del" )

   total.only.ix = is.na( segtab[,"HZ"] )
   hzdel.ix = rep( FALSE, nrow(segtab) )
   hzdel.ix = segtab[,"rescaled_total_cn"] < SCNA_thresholds[["del_thresh_rCN"]] & 
              seg_del_focality > SCNA_thresholds[["del_thresh_foc"]]

   del = as.logical(hzdel.ix)

## AMPS
   ploidy = ABS.dat[["mode.res"]][["mode.tab"]][1,"genome mass"]
   purity = ABS.dat[["mode.res"]][["mode.tab"]][1,"alpha"]
   WGD = as.integer(ABS.dat[["mode.res"]][["mode.tab"]][1,"WGD"])
#   amp = seg_amp_focality > 0.98 & segtab[,"rescaled_total_cn"] >= 5.0

   delta_k = SCNA_thresholds[["amp_delta_k"]]
   t_f = SCNA_thresholds[["amp_thresh_foc"]]
   t_c = SCNA_thresholds[["amp_thresh_log2_rCN"]]
   thresh_slope = SCNA_thresholds[["amp_thresh_slope"]]
## y-intercepts
   b1 = t_f - (thresh_slope * t_c)
   b2 = b1 - delta_k

#      foc1 = amp.gene.data[gene_list[i], pair_SIDs[,1], "amp_foc"]
#      CN1 = amp.gene.data[gene_list[i], pair_SIDs[,1], "rescaled_total_cn"]

   foc = seg_amp_focality
   CN = segtab[,"rescaled_total_cn"]

   ## classify gene amp/foc scores:
   Y_N_1 = thresh_slope * log(CN, 2) + b2
   Y_Q_1 = thresh_slope * log(CN, 2) + b1

## very likely Y_N_1 > Y_Q_1 - written this way for consistency with pairwise / multi-sample callers
   amp = (foc >=  Y_N_1 & foc >= Y_Q_1)

   called_segtab = cbind( segtab, "amp.call"=amp, "del.call"=del, "amp_foc"=seg_amp_focality, "del_foc"=seg_del_focality, "ploidy"=ploidy, "purity"=purity, "WGD"=WGD )

   return( called_segtab )
}


## Convert calls on segs to calls on defined regions e.g. peaks, genes
call_regions_ABSOLUTE_SCNAs = function( segtab, regions )
{
   seg_GRs = GRanges( segtab[,"Chromosome"], IRanges(segtab[,"Start.bp"], segtab[,"End.bp"]) )
#   reg_to_seg = gr.match( regions, seg_GRs )  ## only returns 1st match!!  not cool

# gr.findoverlaps(gr1, gr2) is a faster and less clumsy reimplementation of GRanges findOverlaps
#   reg_to_seg = gr.findoverlaps( regions, seg_GRs )

   res = findOverlaps( regions, seg_GRs )
   reg_to_seg = list( "query.id"=res@queryHits, "subject.id"=res@subjectHits )
   

## Need to resolve multiple matches from regions to genomic segs.   Break ties using CN (highest match for amps, lowest for dels)
   amp.seg.ix = rep(NA, length(regions) )
   del.seg.ix = rep(NA, length(regions) )
   for( i in 1:length(regions) )
   {
      ix = which( reg_to_seg$query.id == i )
      if( length(ix) == 0 ) { next }
      CN_vals = segtab[ reg_to_seg$subject.id[ ix ], "corrected_total_cn"]


## Modify these to require a certain fraction of exons included in the event...
      amp.seg.ix[i] = reg_to_seg$subject.id[ ix[ which.max(CN_vals) ] ]
      del.seg.ix[i] = reg_to_seg$subject.id[ ix[ which.min(CN_vals) ] ]
   }
   
# Some regions may not hit any segs because the reg is between two segs
   missing.ix = which( is.na(amp.seg.ix) )
   for( i in missing.ix )
   {
      dd = distance( regions[i], seg_GRs )
# indices of 2 closest seg matches
      ix2 = order(dd, decreasing=FALSE, na.last=TRUE)[c(1,2)]
 
      CN_vals = segtab[ ix2, "corrected_total_cn"]
      amp.seg.ix[i] = ix2[ which.max(CN_vals) ] 
      del.seg.ix[i] = ix2[ which.min(CN_vals) ] 
   }

# create segtabs containing unique 'best' segment matching each input region
   amp.reg.segtab = segtab[ amp.seg.ix, ]
   del.reg.segtab = segtab[ del.seg.ix, ]

   return( list("amp.reg.segtab"=amp.reg.segtab, "del.reg.segtab"=del.reg.segtab) )
}





## used for Brain met can disc 2015 paper
original_genotype_transcript_SCNAs_in_called_ABS_files = function( ABS_BASE_DIR, regs, SCNA_thresholds, samples=NA, analyst_id="SLC" )
{
   fn_exts = paste(".ABSOLUTE.", analyst_id, ".called.RData", sep="")

   if( is.na(samples) )
   {
      files = grep(  fn_exts, dir(ABS_BASE_DIR, full.names=FALSE), value=TRUE)
      snames = gsub( fn_exts, "", files )
   } 
   else
   {
      files = file.path( ABS_BASE_DIR, paste( samples, fn_exts ) )

      ix = file.exists( files )
      if( any(!ix))
      {
         print( "Missing sample RData files: ")
         print( files[!ix] )
      }

      files = files[ix]
      snames = samples[ix]
   }

   amp_ev_mat = matrix( NA, nrow=length(regs), ncol=length(snames) )
   del_ev_mat = matrix( NA, nrow=length(regs), ncol=length(snames) )

   rownames(amp_ev_mat) = names(regs)
   rownames(del_ev_mat) = names(regs)
   colnames(amp_ev_mat) = colnames(del_ev_mat) = snames


   if( length(files) == 0 ) { stop( paste("No ABSOLUTE result files found in ABS_BASE_DIR ", ABS_BASE_DIR, sep="")) }
### Debugging:
#   i = 24
#   i = grep( "PB0274-TM", files )
#   load( file.path(ABS_BASE_DIR, files[i]) ) 
#   ABS.dat = seg.obj
#   called.segtab = call_genome_wide_ABSOLUTE_SCNAs( ABS.dat )
#   del.reg.segtab = call_regions_ABSOLUTE_SCNAs( called.segtab, regs ) [["del.reg.segtab"]]

#   called.segtabs = list()
#   amp.segtabs = list()
#   del.segtabs = list()
#   for( i in 1:length(files) )
   res = foreach( i = 1:length(files)) %dopar%
   {
      load( file.path(ABS_BASE_DIR, files[i]) ) 
      ABS.dat = seg.obj

      called.segtab = original_call_genome_wide_ABSOLUTE_SCNAs( ABS.dat, SCNA_thresholds )
      res = call_regions_ABSOLUTE_SCNAs( called.segtab, regs ) 

      del.reg.segtab = res[["del.reg.segtab"]]
      amp.reg.segtab = res[["amp.reg.segtab"]]

      cat(".")
      return( list("called.segtab"=called.segtab, "amp.segtab"=amp.reg.segtab, "del.segtab"=del.reg.segtab) )

#      called.segtabs[[i]] = called.segtab
#      amp.segtabs[[i]] = amp.reg.segtab
#      del.segtabs[[i]] = del.reg.segtab
   }

   called.segtabs = lapply( res, "[[", "called.segtab" )
   amp.segtabs = lapply( res, "[[", "amp.segtab" )
   del.segtabs = lapply( res, "[[", "del.segtab" )

   for( i in 1:length(files) )
   {
      amp_ev_mat[,i] = amp.segtabs[[i]][,"amp.call"]
      del_ev_mat[,i] = del.segtabs[[i]][,"del.call"]
   }
   cat("done\n")

   names(called.segtabs) = names(amp.segtabs) = names(del.segtabs) = snames
   SCNA_event_dat = list("amp_ev_mat"=amp_ev_mat, "del_ev_mat"=del_ev_mat, "amp_regs"=regs, "del_regs"=regs )  

# build 3D array: genes X samples X annotations for amp SCNAs  X  refgene transcripts (genes)
   genes = names(regs)  ## exactly the same as del_regs
   N_samps = length(called.segtabs)
   amp.gene.data = array( NA, dim=c(length(genes), N_samps, 5) )
   dimnames(amp.gene.data)[[1]] = genes
   dimnames(amp.gene.data)[[2]] = names(called.segtabs)
   cols = c("amp_foc", "corrected_total_cn", "rescaled_total_cn", "amp.call", "length")
   dimnames(amp.gene.data)[[3]] = cols
   
   for( i in 1:length(cols) )
   {
      amp.gene.data[,,i] = matrix( unlist( lapply( amp.segtabs, "[", cols[i])), nrow=length(genes), ncol=N_samps, byrow=FALSE )
   }


# 3D array for deletions X refgene
   N_samps = length(called.segtabs)
   del.gene.data = array( NA, dim=c(length(genes), N_samps, 5) )
   dimnames(del.gene.data)[[1]] = genes
   dimnames(del.gene.data)[[2]] = names(called.segtabs)
   cols = c("del_foc", "corrected_total_cn", "rescaled_total_cn", "del.call", "length")
   dimnames(del.gene.data)[[3]] = cols
   for( i in 1:length(cols) )
   {
      del.gene.data[,,i] = matrix( unlist( lapply( del.segtabs, "[", cols[i])), nrow=length(genes), ncol=N_samps, byrow=FALSE )
   }

   SCNA_event_dat[["amp.gene.data"]] = amp.gene.data
   SCNA_event_dat[["del.gene.data"]] = del.gene.data
   
   return( list( "called.segtabs"=called.segtabs, "amp.segtabs"=amp.segtabs, "del.segtabs"=del.segtabs, "SCNA_event_dat"=SCNA_event_dat ))
}



select_samples_from_gene_SCNA_calls = function( gene_SCNA_calls, samples )
{

#stop("fix me H.amp")
#[1] "called.segtabs" "amp.segtabs"    "del.segtabs"    "SCNA_event_dat"
   gene_SCNA_calls[["called.segtabs"]] = gene_SCNA_calls[["called.segtabs"]][samples]
   gene_SCNA_calls[["amp.segtabs"]] = gene_SCNA_calls[["amp.segtabs"]][samples]
   gene_SCNA_calls[["del.segtabs"]] = gene_SCNA_calls[["del.segtabs"]][samples]

   gene_SCNA_calls[["amp_ev_mat"]] = gene_SCNA_calls[["amp_ev_mat"]][,samples, drop=FALSE ]
   gene_SCNA_calls[["del_ev_mat"]] = gene_SCNA_calls[["del_ev_mat"]][,samples, drop=FALSE ]

   gene_SCNA_calls[["amp.gene.data"]] = gene_SCNA_calls[["amp.gene.data"]][ , samples, , drop=FALSE ]
   gene_SCNA_calls[["del.gene.data"]] = gene_SCNA_calls[["del.gene.data"]][ , samples, , drop=FALSE ]

   return(gene_SCNA_calls)
}





get_TARGET_SCNA_genes = function()
{
   ## use other cancer DBs to auto-name peaks
   # data("VanAllen2014_TARGET", package="ABSOLUTE")  ## provides TARGET
  load(file.path(pkg_dir, "data", "VanAllen2014_TARGET.RData"))
   target_genes = TARGET

   crit_col = "Types_of_recurrent_alterations"
   del_genes.ix =  c( grep( "Biallelic Inactivation", target_genes[, crit_col], ignore.case=TRUE),
                      grep( "Deletions", target_genes[, crit_col], ignore.case=TRUE )  )  
   amp_genes.ix = grep( "Amplification",  target_genes[, crit_col], ignore.case=TRUE ) 

   del_TARGET = unique( target_genes[ del_genes.ix, "Gene"] )
   amp_TARGET = unique( target_genes[ amp_genes.ix, "Gene"] )

   return( list("amp_genes"=amp_TARGET, "del_genes"=del_TARGET) )
}

get_GISTIC_SCNA_genes = function()
{
#   data( "refgene.hg19.genes", package="ABSOLUTE" )   # provides refgene
#     data( "gencode.hg19.genes", package="ABSOLUTE" )   # provides GENCODE
  load(file.path(pkg_dir, "data", "gencode.hg19.genes.RData"))
    txdb = gencode
    genelist = unique(txdb[,"HGNC"]) 


## provide regs
   # data("Zack2013_GISTIC_regions", package="ABSOLUTE")
  load(file.path(pkg_dir, "data", "Zack2013_GISTIC_regions.RData"))
   amp_GISTIC = names(regs[["pancan"]][["amps"]])
   amp_GISTIC = gsub( "\\[", "", amp_GISTIC )
   amp_GISTIC = gsub( "\\]", "", amp_GISTIC )

   del_GISTIC = names(regs[["pancan"]][["dels"]])
   del_GISTIC = gsub( "\\[", "", del_GISTIC )
   del_GISTIC = gsub( "\\]", "", del_GISTIC )

   amp_GISTIC = intersect( amp_GISTIC, genelist )
   del_GISTIC = intersect( del_GISTIC, genelist )

# Annotated GISTIC peaks in various cancer pubs - curated by amaro
#    data("CN_genes", package="ABSOLUTE")
  load(file.path(pkg_dir, "data", "CN_genes.RData"))

## this is a hack to assign genes curated by aramo to amp or del status
# pare them down to only those genes not in curated GISTIC / TARGET sets
   res = get_TARGET_SCNA_genes()
   amp_TARGET = res$amp_genes
   del_TARGET = res$del_genes
   univ = unique(  c( amp_GISTIC, del_GISTIC, amp_TARGET, del_TARGET) )
   others = c()
   for( i in 1:length(CN_genes) )
   {
      others = c( others, setdiff( CN_genes[[i]], univ ) )
   }
   others = unique(others)
   others = intersect(others, genelist)
   
# fix typos / aliasis in list
# A2BP1 = RBFOX1
   others = c( others, "MIR483", "RBFOX1", "AGBL4")
      
## ASSUME THEY ARE ALL AMPS ??
   amp_genes = unique( c(amp_GISTIC, others ) )
   del_genes = unique( del_GISTIC )

   return( list("amp_genes"=amp_genes, "del_genes"=del_genes) )
}

get_refgene_transcript_GRs = function( genelist=NA )
{
   # data( "refgene.hg19.genes", package="ABSOLUTE" )   # provides refgene
  load(file.path(pkg_dir, "data", "refgene.hg19.genes.RData"))

   if( is.na(genelist) ) { genelist = unique(refgene[,"symb"]) }

   ix = match( genelist, refgene[,"symb"] )
   chr = refgene[ix,"chr"]
   Y.ix = chr %in% c("chrY")
   if( any(Y.ix) ) { print( "Dropping Y chromosome genes:"); print( genelist[Y.ix] ) }
   ix = ix[ !(chr %in% c("chrY"))]
   gene_regs = GRanges( gsub("chr", "", refgene[ix,"chr"]), IRanges(refgene[ix,"start"], refgene[ix,"end"]) )
   names(gene_regs) = refgene[ix,"symb"]

   return(gene_regs)
}


get_GENCODE_transcript_GRs = function( genelist=NA, dropY=TRUE )
{
   # data( "gencode.hg19.genes", package="ABSOLUTE" )   # provides GENCODE
  load(file.path(pkg_dir, "data", "gencode.hg19.genes.RData"))
   txdb = gencode

   if( is.na(genelist) ) { genelist = unique(txdb[,"HGNC"]) }

   if( dropY  )
   {
      ix = match( genelist, txdb[,"HGNC"] )
      chr = txdb[ix,"Chr"]
      Y.ix = chr %in% c("Y")
      if( any(Y.ix) ) { print( "Dropping Y chromosome genes:"); print( genelist[Y.ix] ) }
      ix = ix[ !(chr %in% c("Y"))]

      gene_regs = GRanges( txdb[ix,"Chr"], IRanges(txdb[ix,"Start"], txdb[ix,"End"]) )
      names(gene_regs) = txdb[ix,"HGNC"]
   }
   else
   {
      gene_regs = GRanges( gsub("Chr", "", txdb[,"Chr"]), IRanges(txdb[,"Start"], txdb[,"End"]) )
      names(gene_regs) = txdb[,"HGNC"]
   }
  
   return(gene_regs)
}


get_all_annotated_reg_GRs = function( )
{
x1 = "/xchip/cga/reference/hg19/gencode.v12.gc.txt"

x2 = "/xchip/cga/reference/hg19/gencode.v12.annotation.patched_contigs.gtf"
}







seg_focality_plots = function( PP_CALLS, SCNA_calls, drivers )
{
  ## has the genomic scnas 
   called.segtabs = SCNA_calls[["called.segtabs"]]

   wgd = PP_CALLS[ names(called.segtabs), "Genome doublings" ]
   ploidy = PP_CALLS[ names(called.segtabs), "ploidy"]
   names(ploidy) = rownames(PP_CALLS)

   amp_foc_0 =  unlist(lapply( called.segtabs[wgd==0], "[", "amp_foc" ))
   amp_foc_1 =  unlist(lapply( called.segtabs[wgd>0], "[", "amp_foc" ))
   cn_0 = unlist(lapply( called.segtabs[wgd==0], "[", "corrected_total_cn" ))
   cn_1 = unlist(lapply( called.segtabs[wgd>0], "[", "corrected_total_cn" ))

   rcn_0 = unlist(lapply( called.segtabs[wgd==0], "[", "rescaled_total_cn" ))
   rcn_1 = unlist(lapply( called.segtabs[wgd>0], "[", "rescaled_total_cn" ))

   pdf("seg_cn_vs_foc.pdf", 12, 12 )
   par(mfrow=c(2,2))
   par(  bty="n", las=1  )

# amps
   called_segs_0 = unlist(lapply(called.segtabs[wgd==0], "[", "amp.call"))
   called_segs_1 = unlist(lapply(called.segtabs[wgd>0], "[", "amp.call"))

   seg.colors_0 = rep("black", length(called_segs_0) )
   seg.colors_0[called_segs_0] = "red" 

   seg.colors_1 = rep("black", length(called_segs_1) )
   seg.colors_1[called_segs_1] = "red" 

# corrected_cn (has weighted comp mixture prior)
   plot( log(cn_0,2), amp_foc_0, main="", xlim=c(0, max(log(cn_0,2))), xlab="Log2 corrected CN", ylab="Focality", pch=".", col=seg.colors_0  )
   abline( v=log(7,2), lty=3, col=2 )

   plot( log(cn_1,2), amp_foc_1, main="", xlim=c(0, max(log(cn_1,2))), xlab="Log2 corrected CN", ylab="Focality", pch=".", col=seg.colors_1  )
   abline( v=log(7,2), lty=3, col=2 )


## rescaled_cn (linear rescale of copy-ratio using purity / ploidy)
   plot( log(rcn_0,2), amp_foc_0, main="", xlim=c(0, max(log(rcn_0,2),na.rm=TRUE)), xlab="Log2 rescaled CN", ylab="Focality", pch=".", col=seg.colors_0  )
   abline( v=log(7,2), lty=3, col=2 )

   plot( log(rcn_1,2), amp_foc_1, main="", xlim=c(0, max(log(rcn_1,2),na.rm=TRUE)), xlab="Log2 rescaled CN", ylab="Focality", pch=".", col=seg.colors_1  )
   abline( v=log(7,2), lty=3, col=2 )


## deletions
   called_segs_0 = unlist(lapply(called.segtabs[wgd==0], "[", "del.call"))
   called_segs_1 = unlist(lapply(called.segtabs[wgd>0], "[", "del.call"))

   seg.colors_0 = rep("black", length(called_segs_0) )
   seg.colors_0[called_segs_0] = "dodgerblue" 

   seg.colors_1 = rep("black", length(called_segs_1) )
   seg.colors_1[called_segs_1] = "dodgerblue" 

#   del.segtabs = SCNA_calls[["del.segtabs"]]
   del_foc_0 =  unlist(lapply( called.segtabs[wgd==0], "[", "del_foc" ))
   del_foc_1 =  unlist(lapply( called.segtabs[wgd>0], "[", "del_foc" ))

   plot( cn_0, del_foc_0, main="", xlim=c(0, 4), xlab="corrected CN", ylab="Focality", pch=".", col=seg.colors_0 )
   abline( v=0, lty=3, col=2 )

   plot( cn_1, del_foc_1, main="", xlim=c(0, 4), xlab="corrected CN", ylab="Focality", pch=".", col=seg.colors_1 )
   abline( v=0, lty=3, col=2 )

   plot( rcn_0, del_foc_0, main="", xlim=c(0, 4), xlab="rescaled CN", ylab="Focality", pch=".", col=seg.colors_0 )
   abline( v=0, lty=3, col=2 )

   plot( rcn_1, del_foc_1, main="", xlim=c(0, 4), xlab="rescaled CN", ylab="Focality", pch=".", col=seg.colors_1 )
   abline( v=0, lty=3, col=2 )

   dev.off()







## gene-level plots

   gene_SCNA_dat = SCNA_calls[["SCNA_event_dat"]]
   amp.gene.data = gene_SCNA_dat[["amp.gene.data"]]
   del.gene.data = gene_SCNA_dat[["del.gene.data"]]



   pdf("drivergene_cn_vs_foc.pdf", 12, 12 )
   par(mfrow=c(2,2))
   par(  bty="n", las=1  )

## cut down to driver genes only
#   print("Driver genes missing:")
#   print( setdiff(drivers, genes) )
#   drivers = intersect(drivers, genes)

## subset to specified samples - make 2 plots
   plot_genes_SCNA_dat = function( gene.data, samples, main )
   {
      gene.data = gene.data[, samples, ]

      amp.call = as.vector(  (gene.data[,,"amp.call"]) )
      amp.cn = as.vector( (gene.data[,,"corrected_total_cn"]) )

      
      amp.relative.cn = as.vector( gene.data[,,"corrected_total_cn"] / matrix(ploidy[samples], nrow=dim(gene.data)[1], ncol=dim(gene.data)[2], byrow=TRUE)  )


      amp.foc =  as.vector( (gene.data[,,"amp_foc"]) )
      amp.len.kb =  as.vector( (gene.data[,,"length"]) ) / 1e3

      gene_color = rep( "black", length(amp.call) )
      gene_color[as.logical(amp.call)] = "red"

#   plot( log(amp.cn, 2), amp.foc, main="", xlim=c(0, max(log(amp.cn,2))),  xlab="Log2 corrected CN", ylab="Focality", pch=".", col=col)
      plot( 0, type="n", main=main, xlim=c(0, max(log(amp.cn,2))), ylim=c(0.975,1), xlab="Log2 corrected CN", ylab="Focality")
      text( x=log(amp.cn, 2), y=amp.foc, labels=gene_list, font=3, cex=0.5, col=gene_color)

      plot( 0, type="n", main=main, xlim=c(0, max(log(amp.relative.cn,2))), ylim=c(0.975,1), xlab="Log2 (corrected CN / ploidy)", ylab="Focality")
      text( x=log(amp.relative.cn, 2), y=amp.foc, labels=gene_list, font=3, cex=0.5, col=gene_color)

#      plot( 0, type="n", main=main, xlim=c(0, max(log(amp.cn,2))), ylim=c(0,5000), xlab="Log2 corrected CN", ylab="Length (kb)")
#      text( x=log(amp.cn, 2), y=amp.len.kb, labels=gene_list, font=3, cex=0.5, col=gene_color)
   }  


# cut down to freq amps
   gene_amp_freq = rowSums( amp.gene.data[,,"amp.call"] )
   amp_genes = names(gene_amp_freq)[gene_amp_freq >= 5]
   gene_list = intersect( amp_genes, drivers) 

   driver.amp.gene.data = amp.gene.data[gene_list,,]

   wgd_1_snames = dimnames(amp.gene.data)[[2]][wgd==1]
   wgd_0_snames = dimnames(amp.gene.data)[[2]][wgd==0]

   plot_genes_SCNA_dat( driver.amp.gene.data, wgd_0_snames, "Amps in non-doubled genomes"  )
   plot_genes_SCNA_dat( driver.amp.gene.data,  wgd_1_snames, "Amps in doubled (1) genomes"  )


## Output 1 table for each del gene:   del.dat X sample
   out.dir = ("CN/genetabs")
   dir.create(out.dir, recursive = TRUE)
   for( i in 1:length(gene_list) )
   {
      fn = file.path( out.dir, paste( gene_list[i], ".amp.table.txt", sep="") )
      gd = driver.amp.gene.data[gene_list[i],,]
      ix= order(gd[,"amp.call"], gd[,"corrected_total_cn"], gd[,"amp_foc"])
      write.table( gd[ix,], file=fn, sep="\t", quote=FALSE )
   }


#   dev.off()


## Deletions
## cut down to driver genes only
#   print("Driver genes missing:")
#   print( setdiff(drivers, genes) )
#   drivers = intersect(drivers, genes)

# cut down to freq dels
   gene_del_freq = rowSums( del.gene.data[,,"del.call"] )
   del_genes = names(gene_del_freq)[gene_del_freq >= 5]
   gene_list = intersect( del_genes, drivers) 
#   gene_list = grep("TTT", del_genes, invert=TRUE, value=TRUE)
   



   driver.del.gene.data = del.gene.data[gene_list,,]
   del.call = as.vector(  (driver.del.gene.data[,,"del.call"]) )
   del.cn = as.vector( (driver.del.gene.data[,,"corrected_total_cn"]) )
   del.rcn = as.vector( (driver.del.gene.data[,,"rescaled_total_cn"]) )
   del.foc =  as.vector( (driver.del.gene.data[,,"del_foc"]) )
   del.len.kb =  as.vector( (driver.del.gene.data[,,"length"]) ) / 1e3

   gene_color = rep( "black", length(del.call) )
   gene_color[as.logical(del.call)] = "dodgerblue"

   plot( 0, type="n", main="", xlim=c(-0.25, 1.25), ylim=c(0.975, 1), xlab="Corrected CN", ylab="Focality")
   text( x=del.cn, y=del.foc, labels=gene_list, font=3, cex=0.5, col=gene_color)

   plot( 0, type="n", main="", xlim=c(min(del.rcn), 1.25), ylim=c(0.975, 1), xlab="Rescaled CN", ylab="Focality")
   text( x=del.rcn, y=del.foc, labels=gene_list, font=3, cex=0.5, col=gene_color)


   plot( 0, type="n", main="", xlim=c(-0.25, 1.25), ylim=c(min(del.len.kb), 500), xlab="Corrected CN", ylab="Segment length (kp)")
   text( x=del.cn, y=del.len.kb, labels=gene_list, font=3, cex=0.5, col=gene_color)

   plot( 0, type="n", main="", xlim=c(min(del.rcn), 1.25), ylim=c(min(del.len.kb), 500), xlab="Rescaled CN", ylab="Segment length (kp)")
   text( x=del.rcn, y=del.len.kb, labels=gene_list, font=3, cex=0.5, col=gene_color)
  
   
  # plot del len vs. focality
   plot( 0, type="n", main="", xlim=c(0.975, 1), ylim=c(min(del.len.kb), 500), xlab="Focality", ylab="Segment length (kp)")
   text( x=del.foc, y=del.len.kb, labels=gene_list, font=3, cex=0.5, col=gene_color)
 


   dev.off()


## Output 1 table for each del gene:   del.dat X sample
   out.dir = ("CN/genetabs")
   dir.create(out.dir)
   for( i in 1:length(gene_list) )
   {
      fn = file.path( out.dir, paste( gene_list[i], ".del.table.txt", sep="") )
      gd = driver.del.gene.data[gene_list[i],,]
      ix= order(gd[,"del.call"], gd[,"corrected_total_cn"], gd[,"del_foc"])
      write.table( gd[ix,], file=fn, sep="\t", quote=FALSE )
   }
}



# also calls ABS segs: temp solution until this is done in ABS at extraction
# produce matrix of sample X region  true/false values + metadata
# Takes direction into account - e.g. only call amps in amp regs

# NOT USED currently - qualify events at the gene level instead
genotype_amp_and_del_SCNAs_in_called_ABS_files = function( ABS_BASE_DIR, amp_regs, del_regs )
{
stop("deprecated")
   files = grep(  ".ABSOLUTE.SLC.called.RData", dir(ABS_BASE_DIR, full.names=FALSE), value=TRUE)
   snames = gsub( ".ABSOLUTE.SLC.called.RData", "", files )

   amp_ev_mat = matrix( NA, nrow=length(amp_regs), ncol=length(snames) )
   del_ev_mat = matrix( NA, nrow=length(del_regs), ncol=length(snames) )

   rownames(amp_ev_mat) = names(amp_regs)
   rownames(del_ev_mat) = names(del_regs)
   colnames(amp_ev_mat) = colnames(del_ev_mat) = snames

### Debugging:
#   i = 24
   i = grep( "PB0274-TM", files )
   load( file.path(ABS_BASE_DIR, files[i]) ) 
   ABS.dat = seg.obj
   called.segtab = call_genome_wide_ABSOLUTE_SCNAs( ABS.dat )
   del.reg.segtab = call_regions_ABSOLUTE_SCNAs( called.segtab, del_regs ) [["del.reg.segtab"]]
###

   res = foreach( i = 1:length(files)) %dopar%
   {
      load( file.path(ABS_BASE_DIR, files[i]) ) 
      ABS.dat = seg.obj

      called.segtab = call_genome_wide_ABSOLUTE_SCNAs( ABS.dat )
      del.reg.segtab = call_regions_ABSOLUTE_SCNAs( called.segtab, del_regs ) [["del.reg.segtab"]]
      amp.reg.segtab = call_regions_ABSOLUTE_SCNAs( called.segtab, amp_regs ) [["amp.reg.segtab"]]

      cat(".")
      return( list("called.segtab"=called.segtab, "amp.segtab"=amp.reg.segtab, "del.segtab"=del.reg.segtab) )
   }

   called.segtabs = lapply( res, "[[", "called.segtab" )
   amp.segtabs = lapply( res, "[[", "amp.segtab" )
   del.segtabs = lapply( res, "[[", "del.segtab" )

   for( i in 1:length(files) )
   {
      amp_ev_mat[,i] = amp.segtabs[[i]][,"amp.call"]
      del_ev_mat[,i] = del.segtabs[[i]][,"del.call"]
   }
   cat("done\n")

   names(called.segtabs) = names(amp.segtabs) = names(del.segtabs) = snames
   SCNA_event_dat = list("amp_ev_mat"=amp_ev_mat, "del_ev_mat"=del_ev_mat, "amp_regs"=amp_regs, "del_regs"=del_regs )  
   
   return( list( "called.segtabs"=called.segtabs, "amp.segtabs"=amp.segtabs, "del.segtabs"=del.segtabs, "SCNA_event_dat"=SCNA_event_dat) )
}



## updated to support H.amp calling in single samples
genotype_transcript_SCNAs_in_called_ABS_files = function( ABS_BASE_DIR, regs, SCNA_thresholds, sample_ids=NA, sample_names=NA, analyst_id="SLC" )
{
   fn_exts = paste(".ABSOLUTE.", analyst_id, ".called.RData", sep="")

   if( is.na(sample_ids) )
   {
      files = grep(  fn_exts, dir(ABS_BASE_DIR, full.names=FALSE), value=TRUE)
      SIDs = gsub( fn_exts, "", files )
      snames=SIDs
      files = file.path(ABS_BASE_DIR, files )
#print(SIDs)
   } 
   else
   {
      files = file.path( ABS_BASE_DIR, paste( sample_ids, fn_exts, sep="" ) )

      ix = file.exists( files )
      if( any( !ix ) )
      {
         print( "Missing sample RData files: ")
         print( files[!ix] )
      }

      files = files[ix] 
      SIDs = sample_ids[ix]
      snames = sample_names[ix]
   }

   N_samps = length(SIDs)

   amp_ev_mat = matrix( NA, nrow=length(regs), ncol=N_samps )
   H.amp_ev_mat = matrix( NA, nrow=length(regs), ncol=N_samps )
   del_ev_mat = matrix( NA, nrow=length(regs), ncol=N_samps )

   rownames(amp_ev_mat) = names(regs)
   rownames(H.amp_ev_mat) = names(regs)
   rownames(del_ev_mat) = names(regs)
   colnames(amp_ev_mat) = colnames(H.amp_ev_mat) = colnames(del_ev_mat) = snames

   called.segtabs = list()
   amp.segtabs = list()
   del.segtabs = list()

   if( length(files) == 0 ) { stop( paste("No ABSOLUTE result files found in ABS_BASE_DIR ", ABS_BASE_DIR, sep="")) }

   res = foreach( i = 1:length(files)) %dopar%
#   for( i in 1:length(files) )
   {
      load( files[i] )  
      ABS.dat = seg.obj

      called.segtab = call_genome_wide_ABSOLUTE_SCNAs( ABS.dat, SCNA_thresholds )
      r = call_regions_ABSOLUTE_SCNAs( called.segtab, regs ) 

      del.reg.segtab = r[["del.reg.segtab"]]
      amp.reg.segtab = r[["amp.reg.segtab"]]

      cat(".")
      return( list("called.segtab"=called.segtab, "amp.segtab"=amp.reg.segtab, "del.segtab"=del.reg.segtab) )

#      called.segtabs[[i]] = called.segtab
#      amp.segtabs[[i]] = amp.reg.segtab
#      del.segtabs[[i]] = del.reg.segtab
   }

   called.segtabs = lapply( res, "[[", "called.segtab" )
   amp.segtabs = lapply( res, "[[", "amp.segtab" )
   del.segtabs = lapply( res, "[[", "del.segtab" )

   for( i in 1:length(files) )
   {
      amp_ev_mat[,i] = amp.segtabs[[i]][,"amp.call"]
      H.amp_ev_mat[,i] = amp.segtabs[[i]][,"H.amp.call"]
      del_ev_mat[,i] = del.segtabs[[i]][,"del.call"]
   }
   cat("done\n")

   names(called.segtabs) = names(amp.segtabs) = names(del.segtabs) = snames

# build 3D array: genes X samples X annotations for amp SCNAs  X  refgene transcripts (genes)
   genes = names(regs)  ## exactly the same as del_regs
   N_samps = length(called.segtabs)
   amp.gene.data = array( NA, dim=c(length(genes), N_samps, 5) )
   dimnames(amp.gene.data)[[1]] = genes
   dimnames(amp.gene.data)[[2]] = names(called.segtabs)
   cols = c("amp_foc", "corrected_total_cn", "rescaled_total_cn", "amp.call", "length")
   dimnames(amp.gene.data)[[3]] = cols
   
   for( i in 1:length(cols) )
   {
      amp.gene.data[,,i] = matrix( unlist( lapply( amp.segtabs, "[", cols[i])), nrow=length(genes), ncol=N_samps, byrow=FALSE )
   }

#H.AMPS
   H.amp.gene.data = array( NA, dim=c(length(genes), N_samps, 5) )
   dimnames(H.amp.gene.data)[[1]] = genes
   dimnames(H.amp.gene.data)[[2]] = names(called.segtabs)
   cols = c("amp_foc", "corrected_total_cn", "rescaled_total_cn", "amp.call", "length")
   dimnames(H.amp.gene.data)[[3]] = cols
   
   for( i in 1:length(cols) )
   {
      H.amp.gene.data[,,i] = matrix( unlist( lapply( amp.segtabs, "[", cols[i])), nrow=length(genes), ncol=N_samps, byrow=FALSE )
   }


# DELS
   del.gene.data = array( NA, dim=c(length(genes), N_samps, 5) )
   dimnames(del.gene.data)[[1]] = genes
   dimnames(del.gene.data)[[2]] = names(called.segtabs)
   cols = c("del_foc", "corrected_total_cn", "rescaled_total_cn", "del.call", "length")
   dimnames(del.gene.data)[[3]] = cols
   for( i in 1:length(cols) )
   {
      del.gene.data[,,i] = matrix( unlist( lapply( del.segtabs, "[", cols[i])), nrow=length(genes), ncol=N_samps, byrow=FALSE )
   }

## make H.amp and amp calls exclusive
   mask = H.amp.gene.data[,,"amp.call"]
   mask[is.na(mask)] = FALSE
   amp.gene.data[,,"amp.call"][ mask ] = FALSE

   SCNA_event_dat = list("amp_ev_mat"=amp_ev_mat, "H.amp_ev_mat"=H.amp_ev_mat, "del_ev_mat"=del_ev_mat, "amp_regs"=regs, "del_regs"=regs )  

   SCNA_event_dat[["amp.gene.data"]] = amp.gene.data
   SCNA_event_dat[["H.amp.gene.data"]] = H.amp.gene.data
   SCNA_event_dat[["del.gene.data"]] = del.gene.data

   
   return( list( "called.segtabs"=called.segtabs, "amp.segtabs"=amp.segtabs, "del.segtabs"=del.segtabs, "SCNA_event_dat"=SCNA_event_dat ))
}

call_genome_wide_ABSOLUTE_SCNAs = function( ABS.dat, SCNA_thresholds )
{
   segtab = AllelicGetAbsSegDat(ABS.dat)

   seg_amp_focality = compute_focality_score(segtab, "amp" )
   seg_del_focality = compute_focality_score(segtab, "del" )

   total.only.ix = is.na( segtab[,"HZ"] )
   hzdel.ix = rep( FALSE, nrow(segtab) )
   hzdel.ix = segtab[,"rescaled_total_cn"] < SCNA_thresholds[["del_thresh_rCN"]] & 
              seg_del_focality > SCNA_thresholds[["del_thresh_foc"]]

   del = as.logical(hzdel.ix)

## AMPS
   ploidy = ABS.dat[["mode.res"]][["mode.tab"]][1,"genome mass"]
   purity = ABS.dat[["mode.res"]][["mode.tab"]][1,"alpha"]
   WGD = as.integer(ABS.dat[["mode.res"]][["mode.tab"]][1,"WGD"])

   delta_k = SCNA_thresholds[["amp_delta_k"]]
   t_f = SCNA_thresholds[["amp_thresh_foc"]]
   t_c = SCNA_thresholds[["amp_thresh_log2_rCN"]]
   thresh_slope = SCNA_thresholds[["amp_thresh_slope"]]
## y-intercepts
   b1 = t_f - (thresh_slope * t_c)
   b2 = b1 - delta_k

   foc = seg_amp_focality
   CN = segtab[,"rescaled_total_cn"]

   ## classify gene amp/foc scores:
   Y_N_1 = thresh_slope * log(CN, 2) + b2
   Y_Q_1 = thresh_slope * log(CN, 2) + b1

## very likely Y_N_1 > Y_Q_1 - written this way for consistency with pairwise / multi-sample callers
   amp = (foc >=  Y_N_1 & foc >= Y_Q_1)



##  H.amps
   delta_k = SCNA_thresholds[["H.amp_delta_k"]]
   t_f = SCNA_thresholds[["H.amp_thresh_foc"]]
   t_c = SCNA_thresholds[["H.amp_thresh_log2_rCN"]]
   thresh_slope = SCNA_thresholds[["H.amp_thresh_slope"]]
## y-intercepts
   b1 = t_f - (thresh_slope * t_c)
   b2 = b1 - delta_k

   Y_N_1 = thresh_slope * log(CN, 2) + b2
   Y_Q_1 = thresh_slope * log(CN, 2) + b1
   H.amp = (foc >=  Y_N_1 & foc >= Y_Q_1)


## make amp and H.amp exclusive
   amp[ H.amp ] = FALSE


   called_segtab = cbind( segtab, "amp.call"=amp, "H.amp.call"=H.amp, "del.call"=del, "amp_foc"=seg_amp_focality, "del_foc"=seg_del_focality, "ploidy"=ploidy, "purity"=purity, "WGD"=WGD )

   return( called_segtab )
}

