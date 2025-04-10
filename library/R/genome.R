
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
  agg = aggregate(Start.bp ~ chr, data = chr.arms.dat, FUN = max)
  chr_order <- unique(chr.arms.dat$chr)
  agg$chr <- factor(agg$chr, levels = chr_order)
  pos <- agg[order(agg$chr), "Start.bp"]

  if (x == FALSE) {
    pos <- pos[c(1:22)]
  } else {
    pos <- pos[c(1:23)]
  }

  return(pos)
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
