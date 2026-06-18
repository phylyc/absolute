## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

total_get_copy_ratio_comb = function(Q, delta, b, error_model) {
  xx = (delta * (c(1:Q) - 1) + b)
  return(xx)
}


## Comb for male sex chromosomes (X and Y in males), where the germline copy number is 1.
## The correct correction depends on how the input copy ratio was normalized:
##
##  - TOTAL CR (GATK denoised CR, sex-matched / per-contig normalization): a normal male X/Y
##    reads at the SAME baseline as a normal autosome (germline = 1.0 on the diploid-relative
##    scale). One X copy therefore spans a full CR unit, so the comb keeps the autosome
##    intercept b but DOUBLES the per-copy slope (2*delta). With that, germline X (CN 1) maps
##    to the same CR as autosome diploid (CN 2): 2*delta*1 + b == delta*2 + b.
##
##  - allelic / autosome-normalized data (original ABSOLUTE / HAPSEG): a normal male X reads at
##    HALF the diploid baseline, so the comb keeps slope delta but halves the intercept (b/2).
##
## We branch on obs$data.type ("TOTAL" only for the total CR path); everything else keeps the
## historical b/2 behavior so the allelic path is unchanged.
get_male_sex_chr_comb = function( Q, delta, b, obs )
{
  if( identical(obs[["data.type"]], "TOTAL") ) {
    return( GetCopyRatioComb(Q, 2 * delta, b, obs[["error.model"]]) )
  }
  return( GetCopyRatioComb(Q, delta, b / 2, obs[["error.model"]]) )
}


total_get_cr_grid_from_ccf_grid = function(qc, qs, delta, comb, ccf_grid) {
  d = qs - qc
  d_d = delta * d
  cr_grid = (d_d * ccf_grid) + comb[qc + 1]
  
  if (d < 0) { 
    cr_grid = rev(cr_grid)
  }
  
  return(cr_grid)
}


get_tcr_chr_arm_segs = function(seg.obj, chr_arm_dat) {
  seg.dat = seg.obj$segtab
  
  int_w = c()
  ix = c()
  
  arm_len_bp = chr_arm_dat["End.bp"] - chr_arm_dat["Start.bp"]
  
  for (i in 1:nrow(seg.dat)) {
    if( seg.dat[i,"Chromosome"] != chr_arm_dat["chr"]) { 
      next 
    }
    int_start = max(as.numeric(c(seg.dat[i, "Start.bp"], chr_arm_dat["Start.bp"])))
    int_end = min(as.numeric(c(seg.dat[i, "End.bp"], chr_arm_dat["End.bp"])))
    
    ## seg does not overlap region
    if (int_start > int_end) { 
      next 
    }
    int_len_bp = int_end - int_start
    int_w = c(int_w, int_len_bp / arm_len_bp)
    
    ix = c(ix, i)
  }
  
  int_w = as.numeric(int_w)
  
  return(list("int_W"=int_w, "ix"=ix ))
}


total_calc_chr_arm_distr = function(seg.obj, seg_q, chr_arms_dat) {
  n_arm = nrow(chr_arms_dat)
  chr_arm_tab = array(NA, dim=c(1, n_arm, ncol(seg_q)))
  
  for (i in seq_len(n_arm)) {
    chr_dat = get_tcr_chr_arm_segs(seg.obj, chr_arms_dat[i, ])      
    
    if (length(chr_dat$int_W) == 0 ) { 
      next 
    }
    
    chr_arm = array(NA, dim=c(1, ncol(seg_q)))
    chr_arm[1, ] = colSums(seg_q[chr_dat$ix, , drop=FALSE] * chr_dat$int_W)
    
    chr_arm_tab[, i, ] = 0
    chr_arm_tab[1, i, which.max(chr_arm[1,])] = 1
  }
  
  return(chr_arm_tab)
}


## ---- Total-CR chr-arm minimum-event (parsimony) score ----
## ABSOLUTE breaks the purity/ploidy (WGD) degeneracy by preferring, among modes that fit
## the data comparably, the one implying the fewest copy-number events. The allelic path does
## this with compute_chrarm_ev_score, which needs per-homolog ancestral/derived CN. Total
## copy-ratio data has no homolog split, so this is a total-CN analog:
##
##   events = sum over chr-arms of |modal_total_CN - ancestral_CN|, minimised over a
##            non-doubled genome (ancestral CN 2), one whole-genome doubling (ancestral CN 4),
##            and two doublings (ancestral CN 8) -- each doubling charged one extra "doubling"
##            event -- plus a penalty (lambda_sc) on the genome fraction forced non-clonal
##            (frac.het).
##
## NOTE on WGD2: with the default kQ=8 (CN states 0..7) the WGD2 ancestral (CN 8) sits at/above
## the cap, so it only wins for genomes whose arms sit near the top of the CN range; it rarely
## triggers but is included for completeness so genuinely 2x-doubled (~octoploid) samples can
## be expressed rather than forced into a WGD1 call.
##
## Higher (less negative) score = more parsimonious. WeighSampleModes ranks total-CR modes on
## this score (mode.tab[,"SCNA_min_chrarm_events"]) exactly as the allelic path does.
##
## Without it, total-CR ranking falls back to SCNA_LL alone, which monotonically prefers
## higher ploidy and collapses onto a whole-genome-doubled solution.
##
## CALIBRATION NOTE: lambda_sc=2.5 was calibrated on one sample, where the selected mode is
## stable for lambda_sc in [1.5, 4] and insensitive to wgd_cost. It still needs validation
## across additional total-CR samples before being treated as production-tuned.

total_get_chrarm_modal_CN = function( chr_arm_tab )
{
   n_arm = dim(chr_arm_tab)[2]
   arm_CN = rep(NA, n_arm)
   for( a in 1:n_arm )
   {
      v = chr_arm_tab[1, a, ]
      if( all(is.na(v)) ) { next }       ## arm with no overlapping segments
      arm_CN[a] = which.max(v) - 1        ## modal total CN (0-based)
   }
   return( arm_CN[ !is.na(arm_CN) ] )
}

total_compute_chrarm_ev_score = function( arm_CN, frac_het, SCNA_model, lambda_sc=2.5, wgd_cost=1 )
{
   ev_result = list()
   kQ = SCNA_model[["kQ"]]

   if( length(arm_CN) == 0 || !is.finite(frac_het) ) {
      ev_result[["score"]] = NA
      ev_result[["WGD"]] = 0
      ev_result[["events"]] = NA
      ev_result[["e_wgd"]] = rep(NA_real_, 3)
      return( ev_result )
   }

   arm_CN[ arm_CN > (kQ - 1) ] = kQ - 1
   arm_CN[ arm_CN < 0 ] = 0

   ## events relative to a non-doubled (anc 2), once-doubled (anc 4) or twice-doubled (anc 8)
   ## genome, each doubling charged wgd_cost extra events
   e_wgd <- c( sum(abs(arm_CN - 2)),
               sum(abs(arm_CN - 4)) + wgd_cost,
               sum(abs(arm_CN - 8)) + 2 * wgd_cost )
   events = min( e_wgd )

   ev_result[["WGD"]] = as.integer( which.min(e_wgd) - 1L )   ## 0, 1, or 2 doublings
   ev_result[["events"]] = events
   ev_result[["e_wgd"]] = e_wgd                               ## per-hypothesis event counts (for plotting)
   ## non-clonal genome mass penalty (total-CN analog of the allelic negative-arm penalty)
   ev_result[["score"]] = -(events) - lambda_sc * frac_het * length(arm_CN)

   if( is.na(ev_result[["score"]]) ) { stop("NA total mode ev score!") }
   return( ev_result )
}


## This version is for data from raw CAPSEG - no seg-level std errs!
total_get_abs_seg_dat = function(segobj) {
  seg_q_tab <- segobj[["mode.res"]][["mode_SCNA_models"]][[1]] [["seg.q.tab"]]
  seg_qz_tab <- segobj[["mode.res"]][["mode_SCNA_models"]][[1]] [["seg.qz.tab"]]

  ## Derive the number of clonal CN states from the model rather than hardcoding it.
  ## seg.qz.tab = cbind(clonal CN states, subclonal/z column), so its last column is z.
  Q = ncol(seg_q_tab)
  qq = Q

  # Get column number of the max of each row and the expected
  max_mat <- apply(seg_qz_tab, MARGIN=1, function(x) which.max(x))
  subclonal_ix = (max_mat == ncol(seg_qz_tab))   ## last (z) column dominates -> subclonal

  max_mat = apply(seg_q_tab, MARGIN=1, function(x) which.max(x))

  exp_mat = apply(seg_q_tab, MARGIN=1,
                  function(x) {
                    x <- x[1:qq] / sum(x[1:qq])
                    return(sum(x * c(1:qq)))
                  })

  # seg_list is relevant seg table
  seg_list <- segobj$segtab
  Chromosome = seg_list[,"Chromosome"]
  seg_list = as.matrix( seg_list[,-1] )  ## HACK - take out non-numeric Chromosome field to fix below

  # make vectors of 0s for columns
  modal_cn = vector(mode="numeric", length=nrow(seg_list))
  expected_cn = vector(mode="numeric", length=nrow(seg_list))
  hz = vector(mode="numeric", length=nrow(seg_list))
  subclonal = vector(mode="numeric", length=nrow(seg_list))
  copy_ratio = vector(mode="numeric", length=nrow(seg_list))

  cancer_cell_frac = rep(NA_real_, nrow(seg_list))
  ccf_ci95_low = rep(NA_real_, nrow(seg_list))
  ccf_ci95_high = rep(NA_real_, nrow(seg_list))

  ## Absolute total copy number obtained by inverting the comb (CR = delta*CN + b => CN =
  ## (CR - b)/delta). Needed by the gene-level SCNA genotyper (rescaled_total_cn / HZ) and IGV
  ## output. Male X/Y use the doubled-slope sex-chromosome comb (2*delta), matching the model.
  rescaled_total_cn = rep(NA_real_, nrow(seg_list))
  corrected_total_cn = rep(NA_real_, nrow(seg_list))
  HZ = vector(mode="numeric", length=nrow(seg_list))

  alpha = segobj[["mode.res"]][["mode.tab"]][1, "alpha"]
  tau   = segobj[["mode.res"]][["mode.tab"]][1, "tau"]
  bd = get_b_and_delta(alpha, tau); b = bd$b; delta = bd$delta
  gender = segobj[["gender"]]
  is_male_sex = (!is.na(gender) && gender == "Male") & (Chromosome %in% c("X", "Y", "chrX", "chrY"))

  ## CCF summary for the top mode. In total CR mode the primary subclonal_SCNA_tab slot
  ## holds the total-CN results (see fit_modes_SCNA_models); fall back to NA if unset.
  sc_tab = segobj[["mode.res"]][["subclonal_SCNA_res"]][["subclonal_SCNA_tab"]][1, , ]
  have_sc = !is.null(sc_tab) && is.matrix(sc_tab) && !all(is.na(sc_tab))

  # print modal and expected values and check for hz and LOH
  for (i in seq_len(nrow(seg_list))) {
    cn = seg_list[i, "copy_num" ]
    copy_ratio[i] = round(cn, 5)   ## cn is already tau/2 (copy-ratio scale); do not halve again
    modal_cn[i] = max_mat[i] - 1
    expected_cn[i] = round(exp_mat[i] - 1, 5)
    subclonal[i] = as.numeric(subclonal_ix[i])
    if (have_sc) {
      cancer_cell_frac[i] = round(sc_tab[i, "CCF_hat"], 5)
      ccf_ci95_low[i] = round(sc_tab[i, "CI95_low"], 5)
      ccf_ci95_high[i] = round(sc_tab[i, "CI95_high"], 5)
    }

    use_delta = if (is_male_sex[i]) 2 * delta else delta
    rescaled_total_cn[i] = round(max(0, (copy_ratio[i] - b) / use_delta), 5)
    corrected_total_cn[i] = expected_cn[i]

    if (modal_cn[i] == 0) {
      hz[i] = 1
      HZ[i] = 1
    }
  }

  # round and delete appropriate fields from existing seg table
  ix = which(colnames(seg_list) %in% c("copy_num"))

  tab = round(seg_list[, c(-ix)], 5)
## Add back text Chromosome col
  tab = data.frame( Chromosome, tab, stringsAsFactors=FALSE )

  return(cbind(tab, copy_ratio, modal_cn, expected_cn, subclonal, cancer_cell_frac,
               ccf_ci95_low, ccf_ci95_high, hz, rescaled_total_cn, corrected_total_cn, HZ))
}

