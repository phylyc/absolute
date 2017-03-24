write_IGV_segfile = function(BASE_DIR, OUT_FN)
{   
   #cols = c("sample", "Chromosome", "Start.bp", "End.bp", "total_copy_ratio" )
   cols = c("sample", "Chromosome", "Start.bp", "End.bp", "rescaled_total_cn" )
   
   fn_list = grep( "segtab.txt", dir(BASE_DIR, full.names=TRUE), value=TRUE )
   
   segtab = data.frame()
   for( i in 1:length(fn_list))
   {
      dat = read.delim( fn_list[i], check.names=FALSE, stringsAsFactors=FALSE )
   
      segtab = rbind( segtab, dat[,cols] )
   
      cat(".")
   }
   cat("\n")
   
   write.table( segtab, file=OUT_FN, row.names=FALSE, sep="\t", quote=FALSE )
}

write_paired_sample_IGV_annot_file = function( SIF, OUT_FN, c1="Primary SID", c2 = "Met SID")
{
   new_SIF = rbind( SIF, SIF )
   TRACK_ID = c( SIF[,c1], SIF[,c2] )
   primet = c( rep("Primary", nrow(SIF)), rep("Met", nrow(SIF)) )

   new_SIF = cbind( "TRACK_ID"=TRACK_ID, "primet"=primet, new_SIF )

   write.table( new_SIF, file=OUT_FN, row.names=FALSE,  sep="\t", quote=FALSE )
}


write_single_sample_IGV_annot_file = function( SIF, OUT_FN, id_col="metIDs" )
{
   write.table( SIF, file=OUT_FN, row.names=FALSE,  sep="\t", quote=FALSE )
}
