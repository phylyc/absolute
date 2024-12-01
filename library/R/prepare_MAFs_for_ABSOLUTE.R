classic_CreateMutCnDat <- function(maf, indel.maf, seg.dat, min.mut.af=0, verbose=FALSE) 
{
	maf <- rename_maf_colnames(maf);
	indel.maf <- rename_maf_colnames(indel.maf);

  indel_filters = function(maf)
  {
     class = maf[, "Variant_Classification"]
     n.ix = class %in% c("IGR", "Intron", "5'UTR", "3'UTR", "RNA", "5'Flank")

     msg = paste("Removing ", sum(n.ix), " of ", length(n.ix), " indels in IGR / Intron / UTR / Flank / RNA", sep="")
     print(msg)
     maf = maf[!n.ix,]

## Turn off "don't ask don't tell" for indels; seems hard to hallucinate these reads due to seq errors
     if( FALSE & nrow(maf) > 0)
     {
#        ix = is.na(maf[,"i_judgement"])  ## only here due to forced-calling
#        ix = maf[,"i_judgement"] == "REJECT"  ## only here due to forced-calling
        ix = maf[,"force_called_site"] == "YES"
        ix[is.na(ix)] = FALSE 
        maf[ix,"t_alt_count"] = 0
     }

#     type = maf[, "Variant_Type"] ## "INS" or "DEL"
     return(maf)
  }


  if( !is.null(indel.maf) )
  {
     indel_missing_cols = setdiff(colnames(maf), colnames(indel.maf))
     imc = matrix( NA, nrow=nrow(indel.maf), ncol=length(indel_missing_cols))
     colnames(imc)= indel_missing_cols
     indel.maf = cbind(indel.maf, imc)
     nc = intersect(colnames(indel.maf), colnames(maf)) 
     filtered.indel.maf = indel_filters(indel.maf)
     maf = rbind(maf, filtered.indel.maf[,nc] )
  }

  mut.cn.dat <- maf
  
  if ("total_normals_called" %in% colnames(mut.cn.dat)) {
    ix <- mut.cn.dat[, "total_normals_called"] > 1
    if (verbose) {
      print(paste("Removing ", sum(ix), " of ", length(ix),
                  " mutations due to seen in > 1 normals", 
                  sep = ""))
    }
    mut.cn.dat <- mut.cn.dat[!ix, ]
  }
  


  class = mut.cn.dat[, "Variant_Classification"]
  n.ix = class %in% c("IGR", "Intron", "5'UTR", "3'UTR", "RNA", "5'Flank")
  
  if( sum(n.ix) > 1000 & sum(!n.ix) > 1000 )
  {
     msg = paste("Removing ", sum(n.ix), " of ", length(n.ix), " SSNVs in IGR / Intron / UTR / Flank / RNA", sep="")
     print(msg)
     mut.cn.dat = mut.cn.dat[!n.ix,]
  }


  if ("dbSNP_Val_Status" %in% colnames(mut.cn.dat)) {
    mut.cn.dat[["dbSNP_Val_Status"]][is.na(mut.cn.dat[["dbSNP_Val_Status"]])] <- ""
  }
  
  cols <- colnames(mut.cn.dat)
  
  cix <- which(cols %in% c("i_t_ref_count", "t_ref_count"))
  colnames(mut.cn.dat)[cix] <- c("ref")
  
  cix <- which(cols %in% c("i_t_alt_count", "t_alt_count"))
  colnames(mut.cn.dat)[cix] <- c("alt")
  
  cix <- which(cols %in% c("dbSNP_Val_Status"))
  colnames(mut.cn.dat)[cix] <- "dbSNP"
  
  cix <- which(cols %in% c("Tumor_Sample_Barcode"))
  if (length(cix) > 0) {
    colnames(mut.cn.dat)[cix] <- "sample"
  }

  if( "failure_reasons" %in% cols ) {
     colnames(mut.cn.dat)[ which(colnames(mut.cn.dat)=="failure_reasons") ] = "i_failure_reasons" 
     cols[ which(cols=="failure_reasons") ] = "i_failure_reasons" 
  }

## save actual observed count of alt reads
  mut.cn.dat[["observed_alt"]] = mut.cn.dat[["alt"]]

## Turn on/off forced-calling of SSNVs
  if( FALSE & "i_failure_reasons" %in% cols )
  {
    ix1 = grep( "fstar_tumor_lod", mut.cn.dat[, "i_failure_reasons"] )
    ix2 = which(mut.cn.dat[,"alt"] > 0)
    ix = intersect( ix1, ix2)
    
    if( length(ix) > 0 ) {
       mut.cn.dat[ix, "alt"] = 0
    }

    if( verbose ) {
       msg = paste( length(ix), " mutations rejected for fstar_tumor_lod, alt set to 0", sep="" )
       print(msg)
    }

    is.important = identify_potential_clinically_actionable_mutations( mut.cn.dat )
    save.ix = is.important & mut.cn.dat[,"observed_alt"] > mut.cn.dat[,"alt"] 
    save.ix[is.na(save.ix)] = FALSE
 
    if( any(save.ix) )
    {
       mut.cn.dat[save.ix,"alt"] = mut.cn.dat[save.ix, "observed_alt"]
    }
    if(verbose) 
    {
       print( paste( "Reverting ", sum(save.ix), " important mutations to force-called alt reads:", sep=""))
       if( any(save.ix) ) { print(mut.cn.dat[save.ix, c("Hugo_Symbol", "Variant_Classification")] )  }
    }
  }

  na.contig.ix = is.na(mut.cn.dat[,"Chromosome"])
  if( any(na.contig.ix) ) 
  {
     print( paste( "Removing ", sum(na.contig.ix), " of ", length(na.contig.ix), " mutations with NA Chromosome", sep=""))
     mut.cn.dat = mut.cn.dat[!na.contig.ix,,drop=FALSE]
  }
  
#  na.ix <- apply(is.na(mut.cn.dat[, c("ref", "alt")]), 1, sum) > 0

# rklein-debugging; it appears that mut.cn.dat does not contain columns ref and alt as purported below
# instead it contains columns: ref_count and alt_count
# Altering those column names appropriately
  check_colnames_src = colnames(mut.cn.dat)
  check_colnames_suspected = c("ref_count","alt_count")
  check_colnames_suspected2 = c("ref", "alt")
  if(all(check_colnames_suspected %in% check_colnames_src) && all(check_colnames_suspected2 %in% check_colnames_src)) {
    stop("Assertion Error: Both ref_count, alt_count and alt, ref columns exist; unsure which to use as alt and ref for purposes of count computation")
  }

  if(all(check_colnames_suspected %in% check_colnames_src)) {
    colnames(mut.cn.dat)[which(names(mut.cn.dat) == "ref_count")] <- "ref"
    colnames(mut.cn.dat)[which(names(mut.cn.dat) == "alt_count")] <- "alt"  
  }

# it is also necessary to convert these columns into integers as they are represented as character by default
# which breaks the below arithmetic
  mut.cn.dat[, "alt"] = as.integer(mut.cn.dat[, "alt"])
  mut.cn.dat[, "ref"] = as.integer(mut.cn.dat[, "ref"])

  na.ix <- is.na(mut.cn.dat[, "alt"] + mut.cn.dat[, "ref"])
  if (verbose) {
    print(paste("Removing ", sum(na.ix), " of ", length(na.ix),
                " mutations with NA coverage", 
                sep = ""))
  }
  mut.cn.dat <- mut.cn.dat[!na.ix, ]

## check for negative read-counts (yes, this happens sometimes for indels)
  neg.ix = mut.cn.dat[, "alt"] < 0 
  if (verbose) { print(paste("Setting ", sum(neg.ix), " of ", length(neg.ix), " alt read counts < 0 to 0", sep = "")) }
  mut.cn.dat[neg.ix,"alt"] = 0
  neg.ix = mut.cn.dat[, "ref"] < 0
  if (verbose) { print(paste("Setting ", sum(neg.ix), " of ", length(neg.ix), " ref read counts < 0 to 0", sep = "")) }
  mut.cn.dat[neg.ix,"ref"] = 0
  
  af <- mut.cn.dat[, "alt"] / (mut.cn.dat[, "alt"] + mut.cn.dat[, "ref"])
  ix <- af < min.mut.af
  ix[  (mut.cn.dat[, "alt"] + mut.cn.dat[, "ref"]) == 0 ] = FALSE 

  if (verbose) {
    print(paste("Removing ", sum(ix), " of ", length(ix),
                " mutations due to allelic fraction < ", 
                min.mut.af, sep = ""))
  }

  if (sum(!ix) == 0) {  
    stop("no mutations left!") 
  }
  mut.cn.dat = mut.cn.dat[!ix, , drop=FALSE]
    
  mut.seg.ix <- GetMutSegIx(mut.cn.dat, seg.dat[["segtab"]])  
  T.seg.ix = total_get_mut_seg_ix(mut.cn.dat, seg.dat[["total.seg.dat"]] ) 
  ix <- apply(is.na(mut.seg.ix), 1, sum) == 0 | !is.na(T.seg.ix)

#  if( any(!ix)) { stop("unmapped mutation?") }



   normal_allele_count = rep(2, nrow(mut.cn.dat))
   X.ix = (mut.cn.dat[,"Chromosome"]=="X")
   Y.ix = (mut.cn.dat[,"Chromosome"]=="Y")

   gender = seg.dat$gender
   if( !is.na(gender) && gender == "Male" ) 
   {
      normal_allele_count[X.ix] = 1 
      normal_allele_count[Y.ix] = 1
   }
   if( !is.na(gender) && gender == "Female" ) 
   {
      normal_allele_count[X.ix] = 2
      normal_allele_count[Y.ix] = 0
   }

  mut.cn.dat = data.frame( mut.cn.dat, "normal_allele_count"=normal_allele_count )

  if( any(is.na(normal_allele_count))){ stop("NA normal_allele_count") }

  if (verbose) {
#    print( paste( "Mapped ", sum(male_X), " mutations with male_X status", sep=""))
    print(paste("Removing ", sum(!ix), " unmapped mutations on Chrs: ", sep = ""))
    print(mut.cn.dat[!ix, "Chromosome"])
  }
  if (sum(ix) == 0) {
    stop("No mutations left")
  }

  mut.cn.dat <- cbind(mut.cn.dat, mut.seg.ix, "T.seg.ix"=T.seg.ix )
  mut.cn.dat <- mut.cn.dat[ix, , drop=FALSE]

  mut.cn.dat = select_protein_change_annot_using_COSMIC( mut.cn.dat, verbose=verbose )
  
  return(mut.cn.dat)
}









minimal_CreateMutCnDat = function(maf, indel.maf, seg.dat, verbose=FALSE) 
{
	maf <- rename_maf_colnames(maf);
	indel.maf <- rename_maf_colnames(indel.maf);

  mut.cn.dat = maf

  mut.seg.ix <- GetMutSegIx(mut.cn.dat, seg.dat[["segtab"]])  
  ix <- apply(is.na(mut.seg.ix), 1, sum) == 0

  mut.cn.dat <- mut.cn.dat[ix, ]
  mut.seg.ix <- mut.seg.ix[ix, , drop = FALSE]
  mut.cn.dat <- cbind(mut.cn.dat, mut.seg.ix)


  normal_allele_count = seg.dat[["obs.scna"]][["normal_allele_count"]][mut.seg.ix[,1]] 
  normal_allele_count[ !ix ] = FALSE

  mut.cn.dat = data.frame( mut.cn.dat, "normal_allele_count"=normal_allele_count )

  if( any(is.na(normal_allele_count))){ stop("Found NA normal_allele_count") }

  if (verbose) {
#    print( paste( "Mapped ", sum(male_X), " mutations with male_X status", sep=""))
    print(paste("Removing ", sum(!ix), " unmapped mutations on Chrs: ", sep = ""))
    print(mut.cn.dat[!ix, "Chromosome"])
  }

  if (sum(ix) == 0) {
    stop("No mutations left")
  }

## support synonyms for ref and alt read counts
  cols <- colnames(mut.cn.dat)
  
  cix <- which(cols %in% c("i_t_ref_count", "t_ref_count"))
  colnames(mut.cn.dat)[cix] <- c("ref")
  
  cix <- which(cols %in% c("i_t_alt_count", "t_alt_count"))
  colnames(mut.cn.dat)[cix] <- c("alt")



  return(mut.cn.dat)
  
}










identify_potential_clinically_actionable_mutations = function( maf )
{
   # data("VanAllen2014_TARGET", package="ABSOLUTE")  ## provides TARGET
  load(file.path(pkg_dir, "data", "VanAllen2014_TARGET.RData"))

   target_genes = TARGET

   silent_classes = c("Silent", "3'UTR", "3'Flank", "5'UTR", "5'Flank", "IGR", "Intron", "RNA", "Targeted_Region", "De_novo_Start_InFrame") 

#   LOF_classes = c("Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation") 
   crit_col = "Types_of_recurrent_alterations"

   mut_genes.ix = grep( "Mutation", target_genes[, crit_col] )
   biallelic_genes.ix = grep( "Biallelic Inactivation", target_genes[,crit_col] )
   gene.ix = union(mut_genes.ix, biallelic_genes.ix)
   genes= target_genes[gene.ix, "Gene"]
   
   is.important = maf[,"Hugo_Symbol"] %in% genes & !(maf[,"Variant_Classification"] %in% silent_classes)

   return(is.important)
}


## Selects between Gencode and UNIPROT AA change using COSMIC
select_protein_change_annot_using_COSMIC = function( MAF, verbose=FALSE )
{
   get_COSMIC_count = function( Hugo_Symbol, Protein_Change, pankey_counts )
   {
# count # of identical codon changes in COSMIC
      keys = paste(Hugo_Symbol, "__", Protein_Change, sep="")
      var_count = rep(0, length(keys) )
      ix = which( keys %in% names(pankey_counts) )
      var_count[ix] = pankey_counts[ keys[ix] ]

      return(var_count)
   }

   # data("COSMIC_protein_change_counts_v67_241013", package="ABSOLUTE")
  load(file.path(pkg_dir, "data", "COSMIC_protein_change_counts_v67_241013.RData"))
   pankey_counts = COSMIC_protein_change_counts

# count # of identical codon changes in COSMIC
   keys = paste( MAF[, "Hugo_Symbol"], "__", MAF[, "Protein_Change"], sep="")
   gencode_count = get_COSMIC_count( MAF[, "Hugo_Symbol"], MAF[, "Protein_Change"], pankey_counts )

   if( !all(is.na(MAF[,"Protein_Change"])) )
   {
      res = strsplit( MAF[,"Protein_Change"], "[0-9]+" )
      A1 = unlist( lapply( res, "[", 1 ))
      A2 = unlist( lapply( res, "[", 2 ))
      UNIPROT_Protein_Change = paste( A1, MAF[,"UniProt_AApos"], A2, sep="" )
   }
   else {UNIPROT_Protein_Change = rep(NA, nrow(MAF)) }

   uniprot_count = get_COSMIC_count( MAF[, "Hugo_Symbol"], UNIPROT_Protein_Change, pankey_counts )

   switch_to_uniprot = uniprot_count > gencode_count & MAF[,"Protein_Change"] != "" & !is.na(MAF[,"UniProt_AApos"])

   if( verbose ) 
   {
      print( paste("Switching from GENCODE to UNIPROT for ", sum(switch_to_uniprot), " mutations based on COSMIC counts ", sep=""))

      cat("GENCODE:\n")
      print( cbind( MAF[ switch_to_uniprot, "Hugo_Symbol"], 
                    MAF[ switch_to_uniprot, "Protein_Change"],
                    gencode_count[ switch_to_uniprot ] ) )

      cat("UNIPROT:\n")
      print( cbind( MAF[ switch_to_uniprot, "Hugo_Symbol"], 
                    UNIPROT_Protein_Change[switch_to_uniprot],
                    uniprot_count[ switch_to_uniprot ] ))
   }
   MAF[ switch_to_uniprot, "Protein_Change"] = UNIPROT_Protein_Change[switch_to_uniprot]
   COSMIC_count = gencode_count
   COSMIC_count[ switch_to_uniprot ] = uniprot_count

   MAF = cbind( MAF, "Number of times codon change is in COSMIC"=COSMIC_count )

   return(MAF)
}

# Start_Position and End_Position are fields defined in the MAF
# specificiation; however, ABSOLUTE internally uses different field names.
# Therefore, rename the fields of the maf data.frame
rename_maf_colnames <- function(maf) {
	if (!is.null(maf)) {
		colnames(maf)[colnames(maf) == "Start_Position"] <- "Start_position";
		colnames(maf)[colnames(maf) == "End_Position"] <- "End_position";
	}
	maf
}

