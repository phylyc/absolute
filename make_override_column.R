#!/usr/bin/Rscript
# the only part of this script that changes is the obj.name -> this is the input into the script

library(dplyr)

make_override_column = function(obj.name) {
init.dir = paste(getwd(),"/",sep="")
setwd(init.dir)
results.dir = "ABSOLUTE_results/"
obj.dir = paste(obj.name, "/",sep="")

copy_dir = paste(init.dir,results.dir,obj.dir,sep="")
#print(copy_dir)
all.files = list.files(copy_dir)


system(sprintf("cp %s%s.PP-calls_tab.txt %stemp_review.txt", copy_dir,obj.name,copy_dir))
cat("\nCopying PP.calls file and creating temp_review.txt...")

# read the man_review file and make a new column
man_review = read.delim(file=paste(copy_dir, "temp_review.txt",sep=""), header=TRUE, sep ='\t')
override_man_review = cbind(override = "",man_review)
override_man_review$override = as.character(override_man_review$override)
override_man_review = arrange(override_man_review, sample)

write.table(override_man_review, quote=FALSE,file = paste(init.dir,results.dir,obj.dir,"man_review.txt",sep=""), sep="\t", row.names = FALSE)
cat(sprintf("\nWrote man_review.txt in %s \n",copy_dir))
}
