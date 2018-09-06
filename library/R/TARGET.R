
get_LOF_var_classes = function()
{
   LOF_classes = c("Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", "Splice_Site", "In_frame_Indel")
   return(LOF_classes)
}


prioritize_clinically_actionable_variants = function( maf )
{
   #data("VanAllen2014_TARGET", package="ABSOLUTE")  ## provides TARGET
   #target_genes = TARGET
	target_genes = read.delim("/xchip/scarter/ncamarda/projects/paad_ccpm/PDAC_TARGET_v11_83117.txt",sep="\t", header=T, stringsAsFactors=F)
   LOF_classes = get_LOF_var_classes()
   silent_classes = reject_mutation_classes()

   crit_col = "Types_of_recurrent_alterations"
   del_genes.ix =  c( grep( "Biallelic Inactivation", target_genes[, crit_col], ignore.case=TRUE),
                      grep( "Deletions", target_genes[, crit_col], ignore.case=TRUE )  )  
   amp_genes.ix = grep( "Amplification",  target_genes[, crit_col], ignore.case=TRUE ) 

   del_TARGET = unique( target_genes[ del_genes.ix, "Gene"] )
   amp_TARGET = unique( target_genes[ amp_genes.ix, "Gene"] )

   mut_genes.ix = grep( "Mutation", target_genes[, crit_col] )
   biallelic_genes.ix = grep( "Biallelic Inactivation", target_genes[,crit_col] )

# If a gene needs to be mutated but not biallically, then require a hotspot
   hotspot_TARGET = target_genes[ setdiff(mut_genes.ix, biallelic_genes.ix), "Gene"]
   biallelic_TARGET = target_genes[ biallelic_genes.ix, "Gene"]


##  Now add priority column to MAF
   maf = cbind( maf, "clinical actionability priority"=0 )

#   consequence_mat =  get_variant_consequence_anntoations(maf)
   del.ix = maf[,"Variant_Classification"] == "homozygous deletion" & maf[,"Hugo_Symbol"] %in% c(del_TARGET, biallelic_TARGET)
   amp.ix = maf[,"Variant_Classification"] == "amplification" & maf[,"Hugo_Symbol"] %in% c(amp_TARGET)
   H.amp.ix = maf[,"Variant_Classification"] == "high-level amplification" & maf[,"Hugo_Symbol"] %in% c(amp_TARGET)

   hotspot.ix = maf[,"Number of times codon change is in COSMIC"] > 0 & maf[,"Hugo_Symbol"] %in% c(hotspot_TARGET)
   hotspot.ix[ is.na(maf[,"Number of times codon change is in COSMIC"]) ] = 0
  
# three kinds of biallelic inactivation (in addition to hom del) :
# 1) Het del &  mutation
   maf[ is.na(maf[,"homozygous.ix"]), "homozygous.ix"] = FALSE
   biallelic.ix = maf[,"homozygous.ix"]  & maf[,"Hugo_Symbol"] %in% c(biallelic_TARGET)
print("executing clinical actionability analysis")
   # 2) two het muts in same gene / patient with at least one of them either LOF or COSMIC
# NDC
# 3) Germline mutation and somatic LOH

   if( any( biallelic_TARGET %in% maf[,"Hugo_Symbol"] ) )
   {
      gmaf = maf
      gmaf$Hugo_Symbol = factor(maf$Hugo_Symbol, levels=biallelic_TARGET)
      gmaf$count = 1
      grouped = aggregate(count ~ Hugo_Symbol + pair_id, sum, na.rm=TRUE, data=gmaf)

      if( any(grouped[,"count"]>1) )
      {
         bi.ix = which( grouped[,"count"]>1 )
         grouped[,"Hugo_Symbol"] = as.character( grouped[,"Hugo_Symbol"] )

         for( i in 1:length(bi.ix))
         {
            maf.ix = maf[,"Hugo_Symbol"] == grouped[bi.ix[i], "Hugo_Symbol"] & maf[,"pair_id"] == grouped[bi.ix[i], "pair_id"] 
            if( sum( maf[maf.ix, "Variant_Classification"] %in% LOF_classes ) > 0 |
                max( maf[maf.ix, "Number of times codon change is in COSMIC"], na.rm=TRUE ) > 0 )
            {
               biallelic.ix[maf.ix] = TRUE 
            }
         }
      }
   }
   else{ biallelic.ix = rep(FALSE, nrow(maf)) }

   maf[ del.ix | amp.ix | H.amp.ix | hotspot.ix | biallelic.ix, "clinical actionability priority"] = 1    
   maf[ maf[,"Variant_Classification"] %in% silent_classes,  "clinical actionability priority"] = 0

   return(maf)
}


