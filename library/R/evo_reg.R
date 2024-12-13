get_post_DP_chrarm_states = function( seg.obj, SCNA_model, subclonal_scna_tab, chr.arms.dat )
{
   get_segs_modal_CN = function( CN, w )
   {
      na.ix = is.na(CN)
      CN = CN[!na.ix]
      w = w[!na.ix]

      if(length(CN)==0) { return(NA) }

      uCN = unique(CN)
      CN_w = rep(0, length(uCN))
      for( i in 1:length(uCN)) { CN_w[i] = sum(w[CN==uCN[i]]) }
 
      res = uCN[which.max(CN_w)]
      if( is.na(res)) { stop("NA")}

      return( res )

#      o.ix = order(CN)
#      CN = CN[o.ix]
#      w = w[o.ix]

#      rr = rle(CN)
#      CN_w = rep(NA, length(rr$values) )
#      for( j in 1:length(rr$values) ) { CN_w[j] = sum( w[CN == rr$values[j]] ) }
#      res =  rr$values[ which.max(CN_w) ] 
#
#      if( is.na(res)) { stop("NA")}
#      return( rr$values[ which.max(CN_w) ] )
   }

   get_segs_modal_clust_assign = function( assign, w )
   {
      if(all(is.na(assign))) { return(NA) }

      seg_arm_assign = cbind( "assign"=assign, "w"=w )
      wsum_assign = aggregate( w~assign, data=seg_arm_assign, sum)
      ix = which.max(wsum_assign[,"w"]) 
      return( wsum_assign[ix,"assign"])
   }

# TODO: X-chr
#   chr.arms.dat = chr.arms.dat[ setdiff(rownames(chr.arms.dat), c("Xp", "Xq")), ]
   arm.seg.ix.dat = get_seg_chrarm_ix( seg.obj, chr.arms.dat )

   tree_clust = resort_tree_clusters( SCNA_model, SCNA_model[["seg_CCF_DP"]][["tree_clust"]] )
   SCNA_model[["seg_CCF_DP"]][["tree_clust"]] = tree_clust
#   SCNA_model[["seg_CCF_DP"]][["seg_clust_tab"]] = get_seg_clust_tab( SCNA_model )

   n.seg = nrow(SCNA_model[["seg.ix.tab"]])
   n.arm = length(arm.seg.ix.dat) 
 
## DP decon
#  obs = AllelicExtractSampleObs( seg.obj)
#  allele.segs = get_hom_pairs_segtab_internal( obs, seg.obj[["as.seg.dat"]] )
#  segs_d0 = deconstruct_SCNAs( SCNA_model, obs, allele.segs, b=0, delta=0.5 )

   DP_CN_states = get_DP_anc_der_CN_states( SCNA_model )
   clonal_segs = !DP_CN_states[,"subclonal"] 
   if(any(is.na(clonal_segs))) { stop("NA in clonal_segs") }

   clonal_clust_num = SCNA_model[["seg_CCF_DP"]][["tree_clust"]][["CCF_order"]][1]
   seg_DP_assign = SCNA_model[["seg_CCF_DP"]][["tree_clust"]][["assign"]]
   if( length(seg_DP_assign) != n.seg ) { stop("segDP len mismatch!") }

## treat assign=NA as a distinct state
#   seg_remove.ix = apply( SCNA_model[["seg.ix.tab"]][,c("high.sem.ix", "bi.W.ix")], 1, any)
   seg_amp.ix = SCNA_model[["seg.ix.tab"]][,"amp.ix"]
   seg_neg.ix = SCNA_model[["seg.ix.tab"]][,"neg.ix"]

## below does not use SCNA_model

   amp.assign.num = max(seg_DP_assign, na.rm=TRUE) + 1
   seg_DP_assign[seg_amp.ix] = amp.assign.num 
   seg_DP_assign[seg_neg.ix] = -1
## allow NA assign for uninformative segs
#   seg_DP_assign[seg_remove.ix] = amp.assign.num+1  ## for now
#   if( any( is.na(seg_DP_assign))) { stop("Still NA in seg_DP_assign") }


 ## SC seg deconstruction
## DP version
   CN_states = DP_CN_states

#   CN_states[seg_amp.ix, "qs"] = amp.assign.num

   major_arm_clonal = rep(NA, n.arm)
   minor_arm_clonal = rep(NA, n.arm)

   major_arm_CN = rep(NA, n.arm)
   minor_arm_CN = rep(NA, n.arm)

   major_clust_assign = rep(NA, n.arm)
   minor_clust_assign = rep(NA, n.arm)

   minor_arm_der_anc = matrix( NA, ncol=2, nrow=n.arm )
   major_arm_der_anc = matrix( NA, ncol=2, nrow=n.arm )
   colnames(minor_arm_der_anc) = c("minor_ancestral_CN", "minor_derived_CN")
   colnames(major_arm_der_anc) = c("major_ancestral_CN", "major_derived_CN")
 
   for( i in 1:n.arm)
   {
      maj.ix = arm.seg.ix.dat[[i]][["high.ix"]]
      min.ix = arm.seg.ix.dat[[i]][["low.ix"]]
      w = arm.seg.ix.dat[[i]][["int.w"]]

## decide which CCF clust chrarm is assigned to: choose clust with highest summed w over segs
      major_clust_assign[i] = get_segs_modal_clust_assign( seg_DP_assign[ maj.ix ], w )
      minor_clust_assign[i] = get_segs_modal_clust_assign( seg_DP_assign[ min.ix ], w )

      if( is.na(major_clust_assign[i]) | is.na(minor_clust_assign[i])  )  {   next   }

##  decide on CCF==1 or subclonal for each homolog
      major_clonal_w_sum = sum( w * clonal_segs[maj.ix] )
      minor_clonal_w_sum = sum( w * clonal_segs[min.ix] )

      major_arm_clonal[i] = major_clonal_w_sum > 0.5 | major_clust_assign[i] == clonal_clust_num
      minor_arm_clonal[i] = minor_clonal_w_sum > 0.5 | minor_clust_assign[i] == clonal_clust_num

      if( !minor_arm_clonal[i] & minor_clust_assign[i] != -1 )
      {
         ix = !clonal_segs[min.ix] 
         minor_arm_der_anc[i, "minor_ancestral_CN"] = get_segs_modal_CN(  CN_states[min.ix[ix], "qc"], w[ix] )
         minor_arm_der_anc[i, "minor_derived_CN"] = get_segs_modal_CN(  CN_states[min.ix[ix], "qs"], w[ix] )
      }
      if( minor_arm_clonal[i] )
      {
         ix = clonal_segs[min.ix] 
         minor_arm_CN[i] = get_segs_modal_CN( CN_states[ min.ix[ix]], w[ix]   )
      }


      if( !major_arm_clonal[i] & major_clust_assign[i] != -1  )
      {
         ix = !clonal_segs[maj.ix] 

         if( any( !is.na( CN_states[maj.ix[ix]] ) ) )
         {
            major_arm_der_anc[i, "major_ancestral_CN"] = get_segs_modal_CN(  CN_states[maj.ix[ix], "qc"], w[ix] )
            major_arm_der_anc[i, "major_derived_CN"] = get_segs_modal_CN(  CN_states[maj.ix[ix], "qs"], w[ix] )
         }
      }
      if( major_arm_clonal[i] )
      {
         ix = clonal_segs[maj.ix] 
         major_arm_CN[i] = get_segs_modal_CN( CN_states[ maj.ix[ix]], w[ix]   )
#         if(major_arm_CN[i] == 0) { stop("Clonal assignment of major arm to CN 0!") }
      }
   }

   CN_chrarm_states = cbind( minor_arm_CN, major_arm_CN, minor_arm_clonal, major_arm_clonal, minor_arm_der_anc, major_arm_der_anc, major_clust_assign, minor_clust_assign )
   rownames(CN_chrarm_states) = rownames(chr.arms.dat)


## reassign amps to CN = Q+1
   maj_amp = CN_chrarm_states[,"major_clust_assign"] == amp.assign.num
   CN_chrarm_states[maj_amp, "major_arm_CN"] = SCNA_model[["kQ"]]
   CN_chrarm_states[maj_amp, "major_arm_clonal"] = 1

   min_amp = CN_chrarm_states[,"minor_clust_assign"] == amp.assign.num
   CN_chrarm_states[min_amp, "minor_arm_CN"] = SCNA_model[["kQ"]]
   CN_chrarm_states[min_amp, "minor_arm_clonal"] = 1

   c.ix = as.logical(CN_chrarm_states[,"major_arm_clonal"])
   z.der.ix = CN_chrarm_states[,"major_derived_CN"] == 0
   n.ix = CN_chrarm_states[, "major_clust_assign"] == -1

   remove.ix = apply( is.na(CN_chrarm_states), 1, all )
   if( any(remove.ix) ) { print( paste("Removing ", sum(remove.ix), " of ", length(remove.ix), " chrarms with all NA CN_chrarm_states.", sep="")) }

   CN_chrarm_states = CN_chrarm_states[!remove.ix,]

#   if( any( is.na(!c.ix & !z.der.ix & !n.ix) ) ) { stop("!")}

   if( any(abs(CN_chrarm_states[,"major_ancestral_CN"] - CN_chrarm_states[,"major_derived_CN"]) > 1, na.rm=TRUE)) { warning("Invalid major anc/der state!", immediate.=TRUE)}
   if( any(abs(CN_chrarm_states[,"minor_ancestral_CN"] - CN_chrarm_states[,"minor_derived_CN"]) > 1, na.rm=TRUE)) { warning("Invalid minor anc/der state!", immediate.=TRUE)}

   return(CN_chrarm_states)
}


compute_subclone_CCF_volume = function( SCNA_model, CN_chrarm_states )
{
   tree_clust = resort_tree_clusters( SCNA_model, SCNA_model[["seg_CCF_DP"]][["tree_clust"]] )

   nc = ncol(SCNA_model[["collapsed_CCF_dens"]])
#   GRID = cumsum( c(0, rep(1/(nc-1),(nc-1))))
   log_CCF_prior = log(rep(1/nc, nc))
   clonal_clust = tree_clust[["CCF_order"]][1]
   Q = SCNA_model[["kQ"]]

   arm_clusts = c(CN_chrarm_states[,"major_clust_assign"], CN_chrarm_states[,"minor_clust_assign"] )
   arm_amp = c(CN_chrarm_states[,"major_arm_CN"] == Q, CN_chrarm_states[,"minor_arm_CN"] == Q )
   arm_amp[is.na(arm_amp)] = FALSE
   arm_clusts = arm_clusts[!arm_amp]

   u_arm_clusts = unique( arm_clusts )
# remove clonal cluster and neg arms  ?? AMP?
   u_arm_clusts = setdiff( u_arm_clusts, c(clonal_clust,-1) )
#   if( length(u_arm_clusts)==0 ) { return(0) }   ## all arms assigned to clonal clust

#   arm_SC_clust_log_CCF = log(tree_clust[["CCF_dens"]][u_arm_clusts,,drop=FALSE])
# count all subclone clusters, even if no arms assigned
   arm_SC_clust_log_CCF = log(tree_clust[["CCF_dens"]][,,drop=FALSE])

#   arm_SC_clust_log_Z = tree_clust[["cluster_log_Z"]][u_arm_clusts] 
#   unnorm_CCF_LL = arm_SC_clust_log_CCF + arm_SC_clust_log_Z

  ## compute integral of ratio of prior to posterior 
   lograt = altrat = rep(NA, nrow(arm_SC_clust_log_CCF))
   for( i in 1:length(lograt))
   {
#      lograt[i] = sum( exp(arm_SC_clust_log_CCF[i,] - log_CCF_prior )) * (1/nc)
#      lograt[i] = exp( sum( unnorm_CCF_LL[i,] - log_CCF_prior )) * (1/nc)
      lograt[i] = log( sum(exp(arm_SC_clust_log_CCF[i,]) >= 1/nc) / nc )
      altrat[i] = LogAdd( arm_SC_clust_log_CCF[i,] + log(1/nc) )
   }

#   E_clust_CCF = apply( GRID * SCNA_model[["seg_CCF_DP"]][["tree_clust"]][["CCF_dens"]], 1, sum )
   if( any(is.na(lograt))) { stop("NA lograt!") }

   return(sum(lograt))
}



compute_chrarm_ev_score = function( CN_chrarm_states, SCNA_model, WGD0_Prs, WGD1_Prs, verbose=TRUE )
{
   ev_result = list()
#   WGD0_Prs = get_WGD0_CN_prob_vector(SCNA_model)
#   WGD1_Prs = get_WGD1_CN_prob_vector(SCNA_model)
#   WGD2_Prs = get_WGD2_CN_prob_vector(SCNA_model)

## to score the ancestral copy_state
   N_reg = nrow(CN_chrarm_states)
   minor_ancestral_CN = matrix( NA, nrow=N_reg, ncol=2)
   rownames(minor_ancestral_CN)=rownames(CN_chrarm_states)
   major_ancestral_CN = minor_ancestral_CN
   minor_derived_CN = matrix( NA, nrow=N_reg, ncol=2)
   rownames(minor_derived_CN)=rownames(CN_chrarm_states)
   major_derived_CN = minor_derived_CN
#
   n.ix = CN_chrarm_states[, "minor_clust_assign"] == -1; n.ix = ifelse(is.na(n.ix),T,n.ix)
   c.ix = as.logical(CN_chrarm_states[,"minor_arm_clonal"]); c.ix = ifelse(is.na(c.ix), F, c.ix)
   minor_ancestral_CN[c.ix,1] =  CN_chrarm_states[c.ix,"minor_arm_CN"]
   minor_ancestral_CN[!c.ix,1] = CN_chrarm_states[!c.ix,"minor_ancestral_CN"]  
   z.der.ix = CN_chrarm_states[,"minor_derived_CN"] == 0; z.der.ix = ifelse(is.na(z.der.ix), F, z.der.ix)
   minor_ancestral_CN[!c.ix & !n.ix & !z.der.ix, 2] = CN_chrarm_states[!c.ix & !n.ix & !z.der.ix, "minor_derived_CN"] 
#
   minor_derived_CN[!c.ix,1] =  CN_chrarm_states[!c.ix,"minor_derived_CN"]  
   minor_derived_CN[!c.ix & !n.ix & !z.der.ix, 2] = CN_chrarm_states[!c.ix & !n.ix & !z.der.ix, "minor_ancestral_CN"] 

#
   n.ix = CN_chrarm_states[, "major_clust_assign"] == -1; n.ix = ifelse(is.na(n.ix), T, n.ix)
   c.ix = as.logical(CN_chrarm_states[,"major_arm_clonal"]); c.ix = ifelse(is.na(c.ix), F, c.ix)
 
   major_ancestral_CN[c.ix,1] =  CN_chrarm_states[c.ix,"major_arm_CN"]
   major_ancestral_CN[!c.ix,1] = CN_chrarm_states[!c.ix,"major_ancestral_CN"]  
   z.der.ix = CN_chrarm_states[,"major_derived_CN"] == 0; z.der.ix = ifelse(is.na(z.der.ix), F, z.der.ix)
   major_ancestral_CN[!c.ix & !n.ix & !z.der.ix ,2] = CN_chrarm_states[!c.ix & !n.ix & !z.der.ix,"major_derived_CN"] 
#
   major_derived_CN[!c.ix,1] =  CN_chrarm_states[!c.ix,"major_derived_CN"]  
   major_derived_CN[!c.ix & !n.ix & !z.der.ix, 2] = CN_chrarm_states[!c.ix & !n.ix & !z.der.ix, "major_ancestral_CN"] 

#if( any(major_ancestral_CN==0, na.rm=TRUE)) { stop("??")}

## USE ANCESTRAL CN
   WGD_ev_result = rep(NA,3)
#   WGD_ev_result[1] = anc_ambig_compute_integrated_WGD_LL( minor_ancestral_CN, major_ancestral_CN, WGD0_Prs )
#   WGD_ev_result[2] = anc_ambig_compute_integrated_WGD_LL( minor_ancestral_CN, major_ancestral_CN, WGD1_Prs )
#   WGD_ev_result[3] = compute_WGD_LL( minor_ancestral_CN, major_ancestral_CN, WGD2_Prs )
   WGD_ev_result[3] = -Inf

## now work on derived genome
   major_SC_events = major_derived_CN - major_ancestral_CN 
   major_SC_events[ is.na(major_SC_events)] = 0
   minor_SC_events = minor_derived_CN - minor_ancestral_CN 
   minor_SC_events[ is.na(minor_SC_events)] = 0
# needs some hacking due to the way chrarm segs are constructed
   minor_SC_events[ minor_SC_events < -1 ] = -1
   minor_SC_events[ minor_SC_events > 1 ] = 1
   major_SC_events[ major_SC_events < -1 ] = -1
   major_SC_events[ major_SC_events > 1 ] = 1

## USE DERIVED CN
   WGD_ev_result[1] = anc_ambig_compute_integrated_WGD_LL( minor_ancestral_CN + minor_SC_events, major_ancestral_CN + major_SC_events, WGD0_Prs )
   WGD_ev_result[2] = anc_ambig_compute_integrated_WGD_LL( minor_ancestral_CN + minor_SC_events, major_ancestral_CN + major_SC_events, WGD1_Prs )


## kind of arbitrary for now:
#   major_SC_LL = -sum(abs(major_SC_events[,1]))
#   minor_SC_LL = -sum(abs(minor_SC_events[,1]))

   sc_Pr = get_subclonal_CN_prob_vector(SCNA_model)
   N1 = dim(sc_Pr)[1] 
   N2 = dim(sc_Pr)[2]
   SC_LL = matrix( NA, nrow=N1, ncol=N2)
   for( i in 1:N1 )
   {
      for( j in 1:N2 ) 
      {
         SC_LL[i,j] = sum( c(log(sc_Pr[i,j,minor_SC_events[,1] + 2]), log(sc_Pr[i,j,major_SC_events[,1] + 2])))
      }
   }
#   iSC_LL = LogAdd(as.vector(SC_LL) + log(1/(N1*N2)))    ## Prior over grid is uniform
   iSC_LL = 0

## compute information-theoretic complexity penalty for CCF clusters:
   arm_clust_assign = c( CN_chrarm_states[,"major_clust_assign"],
		         CN_chrarm_states[,"minor_clust_assign"] )

   N = length(arm_clust_assign) 
   xx = rle(sort(arm_clust_assign))
   assign_lens = xx$lengths

   ev_result[["subclone_CCF_volume"]] = compute_subclone_CCF_volume( SCNA_model, CN_chrarm_states )
   ev_result[["SC_ev_steps_LL"]] = iSC_LL #major_SC_LL + minor_SC_LL
 #+ subclone_CCF_volume

   alpha = rep(1, length(assign_lens) )
#   logev = dPolya( assign_lens, alpha, log.p=TRUE )
   A = sum(alpha)
#   ev_result[["log_multinom_coef"]]  = log(N) + lbeta(A,N) - sum( log(assign_lens) + lbeta(alpha, assign_lens) )
   ev_result[["log_multinom_coef"]] =  -(lgamma(N+1) - sum( lgamma(assign_lens+1) ))


   if( verbose ) {
      print( paste("Calculated log_multinom_coef =", round(ev_result[["log_multinom_coef"]],5), " for arm-assignment to ", length(assign_lens), " clusters.", sep=""))
   }

   num_neg_chrarms = sum(CN_chrarm_states[, "major_clust_assign"] == -1, na.rm=T) + 
                     sum(CN_chrarm_states[, "minor_clust_assign"] == -1, na.rm=T) 

   if( verbose ) {
      print( paste("Found ", num_neg_chrarms, " chr-arms assigned < 0.", sep=""))
   }
   ev_result[["neg_arm_LL"]] = (SCNA_model[["NA_chrarm_penalty"]] * num_neg_chrarms) 

   WGD_score =  WGD_ev_result 
 ## normalize WGD-scores to Pr and compute model average
   ev_result[["WGD_Pr"]] = exp(WGD_score - LogAdd(WGD_score))
   WGD_score[ev_result[["WGD_Pr"]]==0] = 0
   ev_result[["average_WGD_score"]] = sum( ev_result[["WGD_Pr"]] * WGD_score )

   best_ev_WGD = which.max(ev_result[["WGD_Pr"]]) - 1
   if(verbose & best_ev_WGD != SCNA_model[["WGD"]] ) { print(paste("Best WGD-model ", best_ev_WGD, " score differs from ABS call: ", SCNA_model[["WGD"]], sep="")) }

   ev_result[["score"]] = ev_result[["average_WGD_score"]] + 
		ev_result[["log_multinom_coef"]] + 
		ev_result[["SC_ev_steps_LL"]] + 
		ev_result[["subclone_CCF_volume"]] +
                ev_result[["neg_arm_LL"]]
	        
#   ev_result = list("score"=mode_score, "WGD_Pr"=WGD_Pr, "ancestral_LL"=average_WGD_score, "log_multinom_coef"=log_multinom_coef, "SC_ev_steps_LL"=SC_ev_steps_LL, "SC_ev_steps_LL"=SC_ev_steps_LL, "subclone_CCF_volume"=subclone_CCF_volume, "num_NA"=num_neg_chrarms) 

   if( is.na(ev_result[["score"]]) ) { stop("NA mode ev score!") }

   return(ev_result) 
}


anc_ambig_compute_integrated_WGD_LL = function( minor_CN, major_CN, WGD_Prs )
{
   get_log_CN_Pr = function( CN, CN_pr )
   {
         gLL = rep(NA, nrow(CN))
         ambig.ix = !is.na(CN[,2])
         gLL[!ambig.ix] = log(CN_pr[CN[!ambig.ix,1]+1])
#         if( any(CN[ambig.ix,2] == 0) ) { stop("modeled 0 anc CN!") }

         if( any(ambig.ix))
         {
            res = cbind( log(CN_pr[CN[ambig.ix,1]+1]),
                         log(CN_pr[CN[ambig.ix,2]+1]) ) + log(1/2)  # prior over amig CN config is uniform
            gLL[ambig.ix] = LogAdd(res)
         }
      return(gLL)
   }

#   if(any(major_CN==0, na.rm=TRUE)) { stop("major_CN==0") }
#sum(log(pr[ dat[,1]+1])*dat[,2] )  ==  sum( log(pr[CN+1]))     ## TRUE
 #
#
   ## Pr major has 1e-5 pr of CN=0 - renormalize
   major_WGD_Prs = WGD_Prs

if( FALSE )
{
   zero_pr =  major_WGD_Prs[,,1]
   nonzero_pr = 1-zero_pr
   nonzero_pr[nonzero_pr<0] = 0
   major_WGD_Prs[,,1] = 1e-5   #  too sensitive to small artifacts in allelic capseg f!=0
   for(i in 1:dim(major_WGD_Prs)[3] )
   {
      major_WGD_Prs[,,i] = major_WGD_Prs[,,i] / nonzero_pr
   }
}

   N1 = dim(WGD_Prs)[1]
   N2 = dim(WGD_Prs)[2]
   LL = array(NA, dim=c(N1,N2))

   for( i in 1:N1 )
   {
      for( j in 1:N2 )
      {
         CN_pr = WGD_Prs[i,j,]
         minor_LL = get_log_CN_Pr( minor_CN, CN_pr )
         major_LL = get_log_CN_Pr( major_CN, major_WGD_Prs[i,j,] )
         LL[i,j] = sum(minor_LL, na.rm=T) + sum(major_LL, na.rm=T)
      }
   }

   LL = LL + log(1/(N1*N2))    ## Prior over grid is uniform
   iLL = LogAdd(as.vector(LL))
   if( is.na(iLL) ) { stop("NA iLL!")} 

   return(iLL)
}




compute_integrated_WGD_LL = function( CN, WGD_Prs )
{
 #
   xx = rle(sort(CN))
   y = cbind(xx[["values"]], xx[["lengths"]])
   dat = cbind( c(0:7), 0 )
   dat[ y[,1]+1, 2] = y[,2]
#
   N1 = dim(WGD_Prs)[1]
   N2 = dim(WGD_Prs)[2]
   LL = array(NA, dim=c(N1,N2))

   for( i in 1:N1 )
   {
      for( j in 1:N2 )
      {
         LL[i,j] = dmultinom( dat[,2], prob=WGD_Prs[i,j,], log=TRUE )
      }
   }

   LL = LL + log(1/(N1*N2))    ## Prior over grid is uniform
   iLL = LogAdd(as.vector(LL))
   return(iLL)
}



compute_WGD_LL = function( CN, WGD_Prs )
{
 #
   xx = rle(sort(CN))
   y = cbind(xx[["values"]], xx[["lengths"]])
   dat = cbind( WGD_Prs[,"allelic_CN"], 0 )
   dat[ y[,1]+1, 2] = y[,2]
#
   LL = dmultinom( dat[,2], prob=WGD_Prs[,"Pr"], log=TRUE )

   return(LL)
}



get_subclonal_CN_prob_vector = function( SCNA_model )
{
   grid = seq(0.1, 2.0, by=0.1 )
   N_grid = length(grid)
#
#   CN_events = c(-1:(SCNA_model[["kQ"]]-1))
   CN_events = c(-1:1)
   score = rep(NA, length(CN_events))
   CN_Prs = array(NA, dim=c( N_grid, N_grid, length(CN_events)) )
#
   for( i in 1:N_grid )
   {
      score[CN_events <= 0] = grid[i]
      for( j in 1:N_grid )
      {
         score[CN_events > 0] = grid[j] 
         LL = -abs(CN_events) * score
         CN_Prs[i,j,] = exp(LL - LogAdd(LL))
      }
   }
   return(CN_Prs)
}





get_WGD0_CN_prob_vector = function( SCNA_model )
{
#   Q = SCNA_model[["kQ"]]

   grid = seq(0.1, 2.0, by=0.1 )
   N_grid = length(grid)
#
   CN_events = c(-1:(SCNA_model[["kQ"]]-1))
   score = rep(NA, length(CN_events))
   CN_Prs = array(NA, dim=c( N_grid, N_grid, length(CN_events)) )
#
   for( i in 1:N_grid )
   {
      score[CN_events <= 0] = grid[i]
      for( j in 1:N_grid )
      {
         score[CN_events > 0] = grid[j] 
         LL = -abs(CN_events) * score
         CN_Prs[i,j,] = exp(LL - LogAdd(LL))
      }
   }
   return(CN_Prs)
#
#   LL = -abs(CN_events) * score
#   obs_Prs = cbind( "allelic_CN"=CN_events+1, "Pr"=exp(LL - LogAdd(LL)) )
#
#   return(obs_Prs)
}



get_WGD1_CN_prob_vector = function( SCNA_model )
{
# integrate over gain/loss rates
   grid = seq(0.1, 2.0, by=0.1 )
   N_grid = length(grid)
#
   Q = SCNA_model[["kQ"]]
   ev_steps = c(-14:7)
   N_steps = length(ev_steps)
#
   CN_after_WGD1 = 2*(ev_steps+1)
   CN_after_WGD1 [CN_after_WGD1 < 0] = Inf
#
   final_cn = outer( CN_after_WGD1, ev_steps, "+" )
   final_cn[final_cn < 0] = Inf
#
   evtab = melt( final_cn )
   ev_events = cbind( "obs"=evtab[,3], "g0"=ev_steps[evtab[,1]], "g1"=ev_steps[evtab[,2]])
#
   z0.ix = 2*(1+ev_events[,"g0"]) == 0
   n.ix =  z0.ix & ev_events[,"g1"] != 0 
   ev_events = ev_events[ !n.ix,]
#
   n.ix = ev_events[,1] < 0 | !is.finite(ev_events[,1])
   ev_events = ev_events[ !n.ix,]
#
#
   obs_Prs = array(NA, dim=c( N_grid, N_grid, length(unique(ev_events[,"obs"])) ) )
   ev_events = cbind( ev_events, "g0_score"=NA, "g1_score"=NA, "Pr"=NA)
   neg_g0 = ev_events[,"g0"] <=0
   neg_g1 = ev_events[,"g1"] <=0
#
   for( i in 1:N_grid )
   {
      ev_events[neg_g0,"g0_score"] = grid[i] * -abs(ev_events[neg_g0,"g0"])
      ev_events[neg_g1,"g1_score"] = grid[i] * -abs(ev_events[neg_g1,"g1"])
#
      for( j in 1:N_grid )
      {
         ev_events[!neg_g0,"g0_score"] = grid[j] * -abs(ev_events[!neg_g0,"g0"])
         ev_events[!neg_g1,"g1_score"] = grid[j] * -abs(ev_events[!neg_g1,"g1"])
#
         LL = ev_events[,"g0_score"] + ev_events[,"g1_score"]
         ev_events[,"Pr"] = exp(LL-LogAdd(LL))
#         CN_Prs[i,j,] = exp(LL - LogAdd(LL))
#
         obs_Prs[i,j,] = dcast( data=as.data.frame(ev_events), formula=obs~., value.var="Pr", fun.aggregate=sum )[,2]
      }
   }

## accumulate PR > Q (amp) in to the Q+1th Pr state
   amp_pr = 1 - apply( obs_Prs[,,1:Q], c(1,2), sum )
   amp_pr[amp_pr<0] = 0  ## roundoff error

   obs_Prs[,,Q+1] = amp_pr
   obs_Prs = obs_Prs[,,1:(Q+1)]

   return(obs_Prs)


#   num_events = rowSums( abs(ev_events[, c("g0", "g1")]) )
#   LL = -num_events
#   Pr = exp(LL-LogAdd(LL))
#   ev_events = cbind(ev_events, Pr)
#
#   obs_Prs = dcast( data=as.data.frame(ev_events), formula=obs~., value.var="Pr", fun.aggregate=sum )
#   rownames(obs_Prs) = NULL
#   colnames(obs_Prs) = c("allelic_CN", "Pr")

## not normalized
#   obs_Prs = obs_Prs[1:Q,]

#   return(obs_Prs)
}


get_WGD2_CN_prob_vector = function( SCNA_model )
{
   Q = SCNA_model[["kQ"]]
   ev_steps = c(-14:7)
   N_steps = length(ev_steps)

   CN_after_WGD1 = 2*(ev_steps+1)
   CN_after_WGD1 [CN_after_WGD1 < 0] = Inf

   CN_after_WGD2 = 2 * outer( CN_after_WGD1, ev_steps, "+" )
   CN_after_WGD2[CN_after_WGD2 < 0] = Inf
   final_cn = array( NA, dim=c(N_steps,N_steps,N_steps) )

   for( i in 1:N_steps )
   {
      final_cn[,,i] = CN_after_WGD2 + ev_steps[i]
   }
#
#   z.ix=which(CN_after_WGD1==0)
## don't come back from 0 after 2nd doubling
#   final_pos = final_cn[ z.ix,,] > 0
#   final_finite = is.finite(final_cn[ z.ix,,])
#   final_cn[z.ix,,][ final_pos & final_finite  ] = 0

   evtab = melt( final_cn )
   ev_events = cbind( "obs"=evtab[,4], "g0"=ev_steps[evtab[,1]], "g1"=ev_steps[evtab[,2]], "g2"=ev_steps[evtab[,3]] ) 

#only moves of 0 are allowed after reaching 0
   if(TRUE)
   {
      z0.ix = 2*(1+ev_events[,"g0"]) == 0
      z1.ix = !z0.ix & 2*(1+ev_events[,"g0"]) + ev_events[,"g1"] == 0
      ev_events[z0.ix,c("obs","g1","g2")]
      ev_events[z1.ix,]
#
      n.ix = z1.ix & ev_events[,"g2"] != 0      |     z0.ix & (ev_events[,"g1"] != 0 | ev_events[,"g2"] != 0)
      ev_events = ev_events[ !n.ix,]
   }
#
   n.ix = ev_events[,"obs"] < 0 | !is.finite(ev_events[,"obs"])
   ev_events = ev_events[ !n.ix,]
#
   num_events = rowSums( abs(ev_events[, c("g0", "g1", "g2")]) )
   LL = -num_events
   Pr = exp(LL-LogAdd(LL))
   ev_events = cbind(ev_events, Pr)
#
   obs_Prs = dcast( data=as.data.frame(ev_events), formula=obs~., value.var="Pr", fun.aggregate=sum )
   rownames(obs_Prs) = NULL
   colnames(obs_Prs) = c("allelic_CN", "Pr")
## not normalized
   obs_Prs = obs_Prs[1:Q,]
#
   return(obs_Prs)
# multinomial draw from this dist:
#   draw = sample( size=78, x=nrow(obs_Prs), prob=obs_Prs[,"Pr"], replace=TRUE ) - 1
}





get_seg_chrarm_ix = function( seg.obj, chr.arms.dat )
{
  n.arm <- nrow(chr.arms.dat)
  
  arm.seg.ix.dat = list()

  for (i in 1:n.arm) 
  {
    chr.dat <- GetAllelicChrArmSegs(seg.obj, chr.arms.dat[i, ])   

    # low.ix, high.ix, int.w
     arm.seg.ix.dat[[i]] = chr.dat
  }   
    
  return(arm.seg.ix.dat)
}


get_DP_anc_der_CN_states = function( SCNA_model )
{
   CN_states = SCNA_model[["CN_states"]]

   tree_clust = SCNA_model[["seg_CCF_DP"]][["tree_clust"]]
   seg_clust_tab = SCNA_model[["seg_CCF_DP"]][["seg_clust_tab"]]
   n.ix = apply( SCNA_model[["seg.ix.tab"]][ , c("amp.ix", "neg.ix" ) ], 1, any )

   seg_modal_clust = tree_clust[["assign"]]
#   seg_modal_clust = rep( NA, nrow(seg_clust_tab) )
#   if( any( !n.ix ) )
#   {
#      seg_modal_clust[!n.ix] = apply( seg_clust_tab[!n.ix,, drop=FALSE], 1, which.max )
#   }

   clonal_clust_num = tree_clust[["CCF_order"]][1]

   clonal.ix = seg_modal_clust == clonal_clust_num 
   clonal.ix[  SCNA_model[["seg.ix.tab"]][,"clonal.ix"] ] = TRUE 

   clonal.ix[is.na(clonal.ix)] = FALSE

   ng = ncol(SCNA_model[["seg_CCF_dens"]])
   clonal.bit = rep(FALSE, length(clonal.ix) )
   clonal.bit[clonal.ix] = SCNA_model[["seg_CCF_dens"]][clonal.ix,1] > SCNA_model[["seg_CCF_dens"]][clonal.ix,ng]

   subclonal.ix = !n.ix & !clonal.ix
   subclonal.ix[is.na(subclonal.ix)] = FALSE

## select ancestral or derived allele setting for clonal segs
   segs_d0 = rep(NA, length(clonal.ix) )
   segs_d0[clonal.ix & clonal.bit] = CN_states[ clonal.ix & clonal.bit , "qc"]
   segs_d0[clonal.ix & !clonal.bit] = CN_states[ clonal.ix & !clonal.bit , "qs"]


   new_qc = rep(NA, length(clonal.ix) )
   new_qs = rep(NA, length(clonal.ix) )
   new_qc[clonal.ix] = segs_d0[clonal.ix]

   new_qc[subclonal.ix] = CN_states[subclonal.ix, "qc"]  
   new_qs[subclonal.ix] = CN_states[subclonal.ix, "qs"]  

   return( cbind( "qc"=new_qc, "qs"=new_qs, "subclonal"=!clonal.ix ) )
}
