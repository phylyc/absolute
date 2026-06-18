
## Number of autosomes for a given chr-arm table. We assume the arm table lists the
## autosomes first, then X, then Y (true for the hg18/hg19/hg38 and mm9/mm10/mm39 reference
## tables shipped in library/data). So n_autosomes = (#distinct chromosomes) - 2.
n_autosomes_from_chr_arms = function(chr.arms.dat)
{
  n_chr <- length(unique(chr.arms.dat$chr))
  return(n_chr - 2L)
}

GetChrLens <- function(chr.arms.dat, x=FALSE, y=FALSE) {
  agg = aggregate(End.bp ~ chr, data = chr.arms.dat, FUN = max)
  chr_order <- unique(chr.arms.dat$chr)
  agg$chr <- factor(agg$chr, levels = chr_order)
  lens <- agg[order(agg$chr), "End.bp"]

  n_auto <- length(chr_order) - 2L     ## autosomes; X and Y assumed last
  if (x & y) { lens <- lens[c(1:(n_auto + 2L))] }
  else if (x) { lens <- lens[c(1:(n_auto + 1L))] }
  else { lens <- lens[c(1:n_auto)] }

  return(lens)
}

GetCentromerePos <- function(chr.arms.dat, x=FALSE, y=FALSE) {
  agg = aggregate(Start.bp ~ chr, data = chr.arms.dat, FUN = max)
  chr_order <- unique(chr.arms.dat$chr)
  agg$chr <- factor(agg$chr, levels = chr_order)
  pos <- agg[order(agg$chr), "Start.bp"]

  n_auto <- length(chr_order) - 2L     ## autosomes; X and Y assumed last
  if (x & y) { pos <- pos[c(1:(n_auto + 2L))] }
  else if (x) { pos <- pos[c(1:(n_auto + 1L))] }
  else { pos <- pos[c(1:n_auto)] }

  return(pos)
}


chromosome_labels = function(x=FALSE, y=FALSE, n_auto=22)
{
   labs = as.character( c(1:n_auto) )
   if(x) { labs = c(labs, "X") }
   if(y) { labs = c(labs, "Y") }
   return(labs)
}

chr2int = function(chr, n_auto=22)
{
  labs = as.character( c(1:n_auto) )
  labs = c(labs, "X", "Y")
  ints = c(1:(n_auto + 2L))
  names(ints) = labs
  name = sub("^chr", "", chr)
  res = ints[name]

  names(res)=NULL
  return(res)
}
