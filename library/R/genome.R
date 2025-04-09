
GetChrLens <- function(chr.arms.dat, x=FALSE) {
  agg = aggregate(End.bp ~ chr, data = chr.arms.dat, FUN = max)
  chr_order <- unique(chr.arms.dat$chr)
  agg$chr <- factor(agg$chr, levels = chr_order)
  lens <- agg[order(agg$chr), "End.bp"]

  if (x == FALSE) {
    lens <- lens[c(1:22)]
  } else {
    lens <- lens[c(1:23)]
  }

  return(lens)
}

GetCentromerePos <- function(chr.arms.dat, x=FALSE) {
  chrarm_names <- paste(c(1:22), "q", sep = "")
  
  if (x) {
    chrarm_names <- c(chrarm_names, "Xq")
  }
  
  cent_pos <- chr.arms.dat[chrarm_names, "Start.bp"]

  if (length(cent_pos) == 0) {
    chrarm_names <- paste0("chr", chrarm_names)
    cent_pos <- chr.arms.dat[chrarm_names, "Start.bp"]
  }
  
  return(cent_pos)
}


chromosome_labels = function(x=FALSE) 
{
   labs = as.character( c(1:22) )
   if(x) { labs = c(labs, "X") }
   return(labs)
}

chr2int = function(chr)
{
  labs = as.character( c(1:22) )
  labs = c(labs, "X", "Y")
  ints = c(1:24)
  names(ints) = labs
  name = sub("^chr", "", chr)
  res = ints[name]

  names(res)=NULL
  return(res)
}
