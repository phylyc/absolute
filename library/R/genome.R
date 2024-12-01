
GetChrLens <- function(x=FALSE) {
  lens <- c(247249719, 242951149, 199501827, 191273063,
            180857866, 170899992, 158821424, 
            146274826, 140273252, 135374737, 134452384,
            132349534, 114142980, 106368585, 
            100338915, 88827254, 78774742, 76117153, 63811651,
            62435964, 46944323, 49691432, 154913754)
    
  if (x == FALSE) {
    lens <- lens[c(1:22)]
  }
    
  return(lens)
}

GetCentromerePos <- function(chr.arms.dat, x=FALSE) {
  chrarm_names <- paste(c(1:22), "q", sep = "")
  
  if (x) {
    chrarm_names <- c(chrarm_names, "Xq")
  }
  
  cent_pos <- chr.arms.dat[chrarm_names, "Start.bp"]
  
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
   labs = c(labs, "X") 
   ints = c(1:23)
   names(ints)=labs  
   res = ints[chr]

   names(res)=NULL
   return(res)
}
