# NF_HD_RNAseq
RNA sequencing analysis pipeline for HD data reproducibility (work in progress) 

### Introduction
The RNA-seq data (GEO identifier: GSE64810)

Still to add to pipeline
* Add data to github

<pre>
[=======================  ] 90%
</pre>

#### Install Nexflow
Install `nextflow` using [their tutorial](https://www.nextflow.io/docs/latest/getstarted.html)

#### Run NF_HD_RNAseq pipeline
Check if java is installed using `java --version`. \
Can install java using `sudo apt install default-jre`. \
`wget -qO- https://get.nextflow.io | bash` \
The easiest way is to install NF_HD_RNAseq is using nextflow's `pull` command:

```
nextflow run jackkelly75/HF_HD_RNAseq
```

### Processes

*buildindex*\
Uses salmon to build index
(Need to have the human transcriptome downloaded already from this bash line "curl ftp://ftp.ensembl.org/pub/release-96/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz -o hsapien.fa.gz")
Index is created using the human transcriptome downloaded and k-mer length of 31 (largest possible size, optimed for reads over ~75bp).


*trimfilter*\
FastqPuri (https://github.com/jengelmann/FastqPuri) used to filter data (Pérez-Rubio *et al*., 2019). trimFilterPE is called to filter the paired end data. --trim ENDSFRAC is called to remove low quality base callings at the end and beggining of reads until the read is above a quality threshold. Accept the trimmed read if the number of low quality nucleotides does not exceed 5%. --trimN ENDS is called to trim N's if found at ends of reads. If the trimmed read length is smaller than the -m set, it is discarded. Trimming and minimum read length are used together as recommended by Williams *et al*., 2016. -l is set to 101 as is read length.

If desired Qreports can be used to create quality reports in HTML. The RNA-seq data has 48 tiles, however 58 samples are run on two lanes and so with the Qreports call these 58 samples must be run with -t (number of tiles) of 96. The results of this were run previously and supplied in case needed rather than run through NextFlow.


*quant*\
Uses salmon to quantify the filtered RNA-seq results. Index was created previously  in *buildindex*/
Salmon is run using:
<pre>
salmon quant -i $index -l A \
            -1 ${reads[0]} -2 ${reads[1]} \
            --seqBias  \  #learn and correct for sequence-specific biases in the input data
            --gcBias \  #learn and correct for fragment-level GC biases in the input data. Does not impact on results if GC bias is not present, only marginally increases run time
            --maxMMPExtension 7 \  # limits the length that a mappable prefix of a fragment may be extended before another search along the fragment is started. Smaller values improve the sensitivity but increase run time.
            --validateMappings\  #selective alignment that is more sensitive
            --threads $task.cpus  -o $pair_id
</pre>



*sortfiles*\
Uses R to sort the salmon quant files and get DEGs.

replace missing pmi data point using KNN
to remove the .*N* at end of the annotation so can be annotated to gene symbol
change ensembl annotation to symbol (any duplicate symbol genes, keep only the one with the highest MAD)
convert transcript to gene controlling for pmi and rin
differential expression using DESeq2 (IHW correction)





### References

Pérez-Rubio P, Lottaz C & Engelmann JC (2019). FastqPuri: high-performance preprocessing of RNA-seq data. *BMC Bioinformatics* 20: Article number 226
Williams CR, Baccarella A, Parrish JZ & Kim CC (2016). Trimming of sequence reads alters RNA-Seq gene expression estimates. *BMC Bioinformatics* 17: Article number 103
