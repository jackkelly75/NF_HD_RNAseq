#!/usr/bin/env nextflow

params.transcriptome = "$baseDir/hsapien.fa.gz"
params.reads = "$baseDir/1_FastqPuri/*_{1,2}_good.fq.gz"
params.outdir = "results"
params.pdata = "$baseDir/data/pData.csv"

log.info """\
 N F - H D - R N A S E Q  P I P E L I N E
 ===================================
 reads        : ${params.reads}
 outdir       : ${params.outdir}
 """
transcriptome_file = file(params.transcriptome)
pdata_file = file(params.pdata)

Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .into { read_pairs_ch; read_pairs2_ch }

process buildIndex {
    tag "$transcriptome.simpleName"

    input:
    file transcriptome from transcriptome_file

    output:
    file 'index' into transcriptome_index

    script:
    """
    salmon index -t $transcriptome -i index -k 31
    """
}  

process quant {
    cpus = 2
    tag "$pair_id"
    publishDir '2_quant'

    input:    
    file index from transcriptome_index
    tuple val(pair_id), path(reads) from read_pairs_ch

    output:
    file(pair_id) into quant_ch

    script:
    """
    salmon quant -l A -i $index -1 ${reads[0]} -2 ${reads[1]} -o $pair_id --validateMappings --seqBias --gcBias
    """
}


process sort_files {
    
    input:
    file("results_list") from quant_ch.collect()
    file pdata from pdata_file
    
    output:

    script:
    """
    #!/usr/bin/Rscript
    
    #import libraries
    library(DESeq2)
    library(biomaRt)
    library(EnsDb.Hsapiens.v86)
    library(IHW)
    library(DMwR)
    library(tximport)
    library(stringr)
    

    pData <- read.csv("$pdata", header=T, row.names = 1)
    sampleTable <- data.frame(pData[,c(1,2,9,10)]) 
    #replace missing pmi data point using KNN
    temp <- sampleTable
    temp[,1] <- gsub("Huntingtons", "1", temp[,1])
    temp[,1] <- gsub("Neurologically_normal", "0", temp[,1])
    for(i in 1:ncol(temp)){
        temp[,i] <- as.numeric(temp[,i])
    }
    
    setwd("$baseDir")    
    set.seed(420)
    png("elbow_plot.png", height = 800, width = 1000)
    temp_na <- na.omit(temp)
    wss <- (nrow(temp_na)-1)*sum(apply(temp_na,2,var))
    for (i in 2:15) wss[i] <- sum(kmeans(temp_na, centers=i)[[4]])
    plot(1:15, wss, type = "b", xlab="Number of Clusters", ylab="Within groups sum of squares")
    dev.off()
    knnOutput <- knnImputation(temp, k = 5)
    sampleTable[,3] <- round(knnOutput[,3], digits = 2)
    sampleTable[,1] <- factor(sampleTable[,1])
    #bin the ages into 3 age bins
    sampleTable[,5] <- cut(sampleTable[,2], breaks=c(0,55,71,200), right = FALSE)
    colnames(sampleTable)[5] <- "binned_age"
    
    setwd("$baseDir/2_quant")  
    #get list of file names
    x <- list.files()
    for(num in 1:length(x)){
        x[num] <- paste(getwd(), "/", x[num], "/quant.sf" ,sep = "")
    }
    #import the quant.sf files and edit gene name to be ENSG number. set prints how many files have been fixed
    set = 0
    for(file in x){
        set = set + 1
        import <- read.table(file, header = TRUE, sep= "\t")
        import[,1] <- substr(import[,1], 1, 15)
        write.table(import, file = file, quote = F, col.names =  T, row.names = F, sep = "\t")
        print(set)
    }
    edbx <- EnsDb.Hsapiens.v86
    tx2gene <- transcripts(edbx, columns = c("tx_id", "gene_id"), return.type = "DataFrame")
    x <- list.files()
    for(num in 1:length(x)){
        x[num] <- paste(getwd(), "/", x[num], "/quant.sf", sep = "")
    }
    txi <- tximport(files = x, type = "salmon", tx2gene = tx2gene)
    for(num in 1:length(x)){
        x[num] <- str_match(x[num], "SRR(.*?)/quant.sf")[,1]
    }
    colnames(txi[[2]]) <- x
    
    MAD <- vector(mode="numeric", length=0)
    for( i in 1:nrow(txi[[2]])){                
      MAD[i] <- mad(txi[[2]][i,1:69])
    }
    ExprsMAD <- cbind(MAD, txi[[2]])
    ExprsMAD <- as.data.frame(ExprsMAD)
    hgnc <- select(EnsDb.Hsapiens.v86, key=rownames(ExprsMAD), columns=c("SYMBOL"), keytype="GENEID")
    rownames(hgnc) <- hgnc[,1]
    hgnc <- hgnc[rownames(ExprsMAD),] #put it into the same order as the results object
    ExprsMAD <- cbind(hgnc[,2], ExprsMAD)
    ExprsMAD = ExprsMAD[order(ExprsMAD[,1], abs(ExprsMAD[,2]), decreasing = TRUE), ]
    entrezID = unique(ExprsMAD[,1])
    id = match(entrezID, ExprsMAD[,1])
    ExprsMAD = ExprsMAD[id[!is.na(id)], ]
    sum(duplicated(ExprsMAD[,1]), na.rm = TRUE) #0
    txi[[2]] <- txi[[2]][rownames(ExprsMAD),]
    txi[[1]] = txi[[1]][rownames(ExprsMAD),]
    txi[[3]] = txi[[3]][rownames(ExprsMAD),]
    
    setwd("$baseDir")
    save(counts, file = "counts.Rdata")
    save(txi, file = "txi.Rdata")
    colnames(txi[[2]]) <- rownames(sampleTable)
    sampleTable[,6] <- relevel(sampleTable[,1], "Neurologically_normal")
    colnames(sampleTable)[6] <- 'condition'
    dds <- DESeqDataSetFromTximport(txi, sampleTable, ~rin + pmi + binned_age + condition)
    dds <- DESeq(dds)
    png("expression_boxplots.png", width = 1000, height = 500)
    par(mar = c(8,5,2,2))
    boxplot(log10(assays(dds)[["cooks"]]), range = 0 , las =2)
    dev.off()
    res <- results(dds)
    resOrdered <- res[order(res[,5]),]
    save(dds, file = "dds.Rdata")
    hgnc <- select(EnsDb.Hsapiens.v86, key=rownames(resOrdered), columns=c("SYMBOL"), keytype="GENEID")
    rownames(hgnc) <- hgnc[,1]
    hgnc <- hgnc[rownames(resOrdered),] #put it into the same order as the results object
    resOrdered[,7] <- hgnc[,2]
    colnames(resOrdered)[7] <- "genes"
    write.csv(as.data.frame(resOrdered), file="DESeq2_HD_results_2.csv")
    resSig <- subset(resOrdered, padj < 0.05)
    write.csv(as.data.frame(resSig), file="DESeq2_HD_results_adjPvalue_0.05_2.csv")
    resSig <- subset(resOrdered, pvalue < 0.01)
    write.csv(as.data.frame(resSig), file="DESeq2_HD_results_Pvalue_0.01_2.csv")
    resSig <- subset(resOrdered, pvalue < 0.05)
    write.csv(as.data.frame(resSig), file="DESeq2_HD_results_Pvalue_0.05_2.csv")
    resIHW <- results(dds, filterFun=ihw)
    summary(resIHW)
    sum(resIHW[,6] < 0.05, na.rm=TRUE)
    ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version=86)
    hgnc <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol'),filters = 'ensembl_gene_id', values = rownames(resIHW), mart = ensembl)
    rownames(hgnc) <- hgnc[,1]
    hgnc <- hgnc[rownames(resIHW),] #put it into the same order as the results object
    resIHW[,8] <- hgnc[,2]
    colnames(resIHW)[8] <- "genes"
    resIHW <- resIHW[order(resIHW[,6]),]
    save(resIHW, file = "resIHW.Rdata")
    write.csv(as.data.frame(resIHW), file="DESeq2_HD_results_IHW.csv")
    resSig <- subset(resIHW, padj < 0.05)
    write.csv(as.data.frame(resSig), file="DESeq2_HD_results_IHW_0.05.csv")
    resSig[,9] <- 2 ^ resSig[,2]
    colnames(resSig)[9] <- "Fold_change"
    x <- cbind(rownames(resSig), resSig)
    colnames(x)[1] <- "Ensembl ID" 
    colnames(x)[9] <- "Symbol"
    colnames(x)[7] <- "IHW pvalue"
    write.csv(x, file = "DEGs.csv", row.names = F)
    """
}
