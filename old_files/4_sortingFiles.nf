#!/usr/bin/env nextflow


process sort_files {
	
	input:

	output:

	script:
	"""

	#!/usr/bin/Rscript
	setwd("/home/jack/Downloads/temp/data")     
	x <- list.files()
    for(num in 1:length(x)){
    	x[num] <- paste(getwd(), "/", x[num], "/quant.sf" ,sep = "")
    }
    set = 0  #set this to 0 so can count the files in R so have an idea of how far through it is
    for(file in x){
    	set = set + 1
    	import <- read.table(file, header = TRUE, sep= "\t")
    	import[,1] <- substr(import[,1], 1, 15)
    	write.table(import, file = file, quote = F, col.names =  T, row.names = F, sep = "\t")
    	print(set)
    }
	"""
}
