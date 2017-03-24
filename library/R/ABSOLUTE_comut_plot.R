## Just the genes x samples matrix
ABSOLUTE_mutation_VCF_matrix_plot = function( ABSOLUTE_VCF_data, power_threshold=0.99, print_sample_names=FALSE )
{
   plot_mut_symbol = function( xleft, xwd, ytop, ywd, ypos, Variant_Classification_Num, annot_classes, annot_cex, annot_pch, hom_ssnv )
   {
      if( ypos == "top" ) { ycrd = ytop - ywd*0.2 }
      if( ypos == "bottom" ) { ycrd = ytop - ywd*0.8 }
      if( ypos == "center" ) { ycrd = ytop - ywd*0.5 }

      var_levels = get_variant_class_levels()
      var_class = var_levels[Variant_Classification_Num]

      if( !var_class %in% annot_classes ) { stop( paste("Undefined var_class: ", var_class, sep="")) }
      k = which( annot_classes == var_class )

      color = ifelse( var_class == "homozygous deletion" | hom_ssnv, "black", "white" )

      if( var_class == "homozygous deletion" )
      {
         points( xleft + xwd*0.5, ycrd, bg=color, col=NA, pch=annot_pch[[k]], cex=annot_cex[[k]] )
      }
      else
      {
         points( xleft + xwd*0.5, ycrd, col=color, pch=annot_pch[[k]], cex=annot_cex[[k]] )
      }

   }

   VCF = ABSOLUTE_VCF_data
   var_mat = VCF[["var_mat"]]
   secondary_var_mat = VCF[["sec_var_mat"]]
#   underpowered_VCF = VCF[["underpowered_matrix"]]

   N_genes = dim(var_mat)[1]
   N_samps = dim(var_mat)[2]
#   N_branches = dim(var_mat)[3]

## look for collisions -  same gene/samp mutated on >1 ev branch:
   mut_count = matrix(NA, N_genes, N_samps)
   for( i in 1:N_samps)
   {
      mut_count[,i] = apply( var_mat[,i,drop=FALSE]>1, 1, sum )
   }

   annot_classes = c("homozygous deletion", "high-level amplification", "amplification", "Frame_Shift", "Nonsense", "Missense_COSMIC_site", "In_frame_Indel", "Splice_Site", "Missense", "Other_non_syn.", "lincRNA", "Silent" )
   annot_pch = list( 25, "^", 17, "!", "*", 18, "@", "$", 18, "o", 'l', '' )
   annot_cex = list( 0.5, 0.65, 0.5, 0.7, 0.8, 0.65, 0.5, 0.5, 0.25, 0.6, 0.6, 0 )

## Just draw a blank matrix 1st
   image( t(matrix(1,nrow=N_genes, ncol=N_samps) ), col=c("grey94"), xlab="", ylab="", axes=FALSE)

   yc = par("usr")[c(3,4)]
   yrng = yc[2] - yc[1]
   xc = par("usr")[c(1,2)]
   xrng = xc[2] - xc[1]

   ## Now draw rectangles for each mutation..
   for( i in 1:N_genes )
   {
      for( j in 1:N_samps )
      {
         xleft  = xc[1] + xrng/N_samps * (j-1)
         xright = xc[1] + xrng/N_samps * j
         xwd = xrng/N_samps
         ybottom = yc[1] + yrng/N_genes * (i-1)
         ytop    = yc[1] + yrng/N_genes * i
         ywd = yrng/N_genes

         if( mut_count[i,j] > 0 )   ## No collision
         {
#            br.ix = which(var_mat[i,j,] > 1 ) 
#            if( underpowered_VCF[i,j] == 0 )
            if( TRUE )
            {
#               color = ev_branch_colors[ br.ix ]
               color =  "grey50"
               rect( xleft = xleft, xright = xright, ybottom = ybottom, ytop = ytop, col = color, border=NA )

               if( secondary_var_mat[i,j] <= 1 )  ## only 1 mutation
               {
                  plot_mut_symbol( xleft, xwd, ytop, ywd, "center", var_mat[i,j], annot_classes, annot_cex, annot_pch, VCF[["hom_ssnv_matrix"]][i,j] )
               }
               else # >1 mutation in this gene,sample,branch; plot two
               {
                  plot_mut_symbol( xleft, xwd, ytop, ywd, "top", var_mat[i,j], annot_classes, annot_cex, annot_pch, VCF[["hom_ssnv_matrix"]][i,j] )
                  plot_mut_symbol( xleft, xwd, ytop, ywd, "bottom", secondary_var_mat[i,j], annot_classes, annot_cex, annot_pch, VCF[["sec_hom_ssnv_matrix"]][i,j])
               }
            }  
         }  # else { stop("Branch collision???") }
      } 
   }

   axis( side=2, at= yc[1] + yrng/N_genes * c(1:N_genes-0.5), labels=dimnames(var_mat)[[1]], font=3, las=1, cex.axis=1.0, tck=-0.005, mgp=c(3, 0.3, 0), lwd=0.5 )

   if( print_sample_names ) 
   {
      axis( side=1, at= xc[1] + xrng/N_samps * c(1:N_samps-0.5), labels=dimnames(var_mat)[[2]], font=2, las=3, cex.axis=0.25, tck=-0.005, mgp=c(3, 0.3, 0), lwd=0.25 )
   }
}




ABSOLUTE_comut_plot = function( ABS_VCF_data, gene_list, clinical_data, clin_tracks, track_colors, sort_column, PDF_OUT_FN )
{
   samples = colnames(ABS_VCF_data[["var_mat"]])
   genes = rownames(ABS_VCF_data[["var_mat"]])

   clinical_tracks = get_clinical_tracks( samples, clinical_data, clin_tracks )
   sample_color_matrix = get_clinical_color_matrix( clinical_tracks, clin_tracks, track_colors )
   rownames(sample_color_matrix) = clinical_data[,"ID"]

   m.ix = match( samples, clinical_data[,"ID"] )
   keys = clinical_data[m.ix, sort_column ]

   gene_list = rev(gene_list)   ## natural display order - 1st gene on top

   VCF = reduce_Phylogic_VCF( ABS_VCF_data, gene_list, samples )
   VCF = remove_empty_genes_Phylogic_VCF( VCF ) 

      N_genes = dim(VCF[["var_mat"]])[1]
      N_samps = dim(VCF[["var_mat"]])[2]
      genes = rownames(VCF[["var_mat"]])

      o.ix = sort_VCF_by_genemuts( VCF, keys )
      o.samples = samples[o.ix]
      VCF = reduce_Phylogic_VCF( VCF, genes, o.samples )

   ## points at which histology changes
#      rr = rle( sample_color_matrix[o.samples,"Histology group"] )
#      br.pts = cumsum(rr[["lengths"]])
#      br.pts = br.pts[ -length(br.pts)]
#      if( length(br.pts)==0 ) { br.pts = NA }



#      simple_Phylogic_comut_plot( VCF, sample_color_matrix[o.ix,, drop=FALSE], ev_branch_colors, PDF_OUT_FN, sample_SSNV_rates_by_branch[, o.samples],  breaks=br.pts  )




   pdf(PDF_OUT_FN, 8.5, 11 )

   ABSOLUTE_mutation_VCF_matrix_plot( VCF, power_threshold=0.99, print_sample_names=TRUE )


   dev.off()

}



reduce_Phylogic_VCF = function( Phylogic_VCF_data, genes, samples, verbose=TRUE )
{
   if( any( !(genes %in% dimnames(Phylogic_VCF_data[["var_mat"]])[[1]])  ) ) 
   {
      if(verbose)
      {
         print("Dropping genes:") 
         print( setdiff( genes, dimnames(Phylogic_VCF_data[["var_mat"]])[[1]]) ) 
      }  
      genes = intersect( genes,  dimnames(Phylogic_VCF_data[["var_mat"]])[[1]] ) 
#     stop("Missing genes")
   }
   if( any( !(samples %in% dimnames(Phylogic_VCF_data[["var_mat"]])[[2]])  ) ) { stop("Missing samples") }

   if( !is.na(dim(Phylogic_VCF_data[[1]])[3]) )
   {
      for( i in 1:length(Phylogic_VCF_data) )
      {
         Phylogic_VCF_data[[i]] = Phylogic_VCF_data[[i]][ genes, samples, , drop=FALSE ]
      }
   }
   else
   {
      for( i in 1:length(Phylogic_VCF_data) )
      {
         Phylogic_VCF_data[[i]] = Phylogic_VCF_data[[i]][ genes, samples, drop=FALSE ]
      }
   }


   return( Phylogic_VCF_data )
}



assemble_ABSOLUTE_VCF_data = function( combined_ABS_VCF, sample_names )
{
   pri_VCF = combined_ABS_VCF[["primary_multi_VCF"]]
   sec_VCF = combined_ABS_VCF[["secondary_multi_VCF"]]

   N_genes = nrow(pri_VCF[,,1])
   N_samps = ncol(pri_VCF[,,1])

   
## create matrix of variant types in  genes X samples X phy brs 
   var_mat = array( NA, dim=c(N_genes, N_samps ) )
   hom_ssnv_matrix = array(0, dim=c(N_genes, N_samps ) )

   sec_var_mat = array( NA, dim=c(N_genes, N_samps ) )
   sec_hom_ssnv_matrix = array(0, dim=c(N_genes, N_samps ) )

   clinically_actionable_matrix = array(0, dim=c(N_genes, N_samps ) )
   genes = rownames(pri_VCF)

   var_mat = pri_VCF[,,"var_class"]
   hom_ssnv_matrix =  pri_VCF [,,"homozygous.ix"] == 2 

   sec_var_mat = sec_VCF[,,"var_class"]
   sec_hom_ssnv_matrix =  sec_VCF[,,"homozygous.ix"] == 2 

   clin_act = array( NA, dim=c( N_genes, N_samps, 2) )
   clin_act[,,1] = pri_VCF[,,"clinical actionability priority"]
   clin_act[,,2] = sec_VCF[,,"clinical actionability priority"]
   clinically_actionable_matrix = apply( clin_act, c(1,2), max )


   VCF_data = list("var_mat"=var_mat, "hom_ssnv_matrix"=hom_ssnv_matrix, "sec_var_mat"=sec_var_mat, "sec_hom_ssnv_matrix"=sec_hom_ssnv_matrix, "clinically_actionable_matrix"=clinically_actionable_matrix )

## insert dimnames into all VCF dims
   for( j in 1:length(VCF_data) )
   {
      dimnames( VCF_data[[j]])[[1]] = genes
      dimnames( VCF_data[[j]])[[2]] = sample_names
   }

   VCF = remove_empty_genes_Phylogic_VCF( VCF_data ) 

   return(VCF) 
}


## returns sample sort order by keys and then mut#
sort_VCF_by_genemuts = function( VCF, pri.keys )
{
   gene_samp_mut_counts = apply( VCF[["var_mat"]] > 1, c(1,2), sum )
   gene_mut_counts = rowSums(gene_samp_mut_counts>0)  ## gene totals
   g.o.ix = order(gene_mut_counts, names(gene_mut_counts), decreasing=TRUE )
   sorted_gene_samp_counts = gene_samp_mut_counts[g.o.ix, ,drop=FALSE]

   key_df = data.frame( pri.keys, t(-sorted_gene_samp_counts) )
   o.ix = do.call( order, key_df )

   return(o.ix)
}

