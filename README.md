# NF_HD_RNAseq
RNA sequencing analysis pipeline for HD data reproducibility

## Introduction
The RNA-seq data (GEO identifier: GSE64810)

3 steps
#filtering of samples using FastqPuri
#quantification using salmon
#identifying DEGs using DESeq2


## Analysing data
\

### Download data

link to download data

### Data pre-processing
*FastqPuri*\
\
FastqPuri (https://github.com/jengelmann/FastqPuri) used to filter data (Pérez-Rubio *et al*., 2019). trimFilterPE is called to filter the paired end data. --trim ENDSFRAC is called to remove low quality base callings at the end and beggining of reads until the read is above a quality threshold. Accept the trimmed read if the number of low quality nucleotides does not exceed 5%. --trimN ENDS is called to trim N's if found at ends of reads. If the trimmed read length is smaller than the -m set, it is discarded. Trimming and minimum read length are used together as recommended by Williams *et al*., 2016. -l is set to 101 as is read length.

If desired Qreports can be used to create quality reports in HTML. The RNA-seq data has 48 tiles, however 58 samples are run on two lanes and so with the Qreports call these 58 samples must be run with -t (number of tiles) of 96. The results of this were run previously and supplied in case needed rather than run through NextFlow.

FastqPuri does not work well through NextFlow so is run through docker instead.

### Install docker
This link describes it well:
https://www.digitalocean.com/community/tutorials/how-to-install-and-use-docker-on-ubuntu-18-04



```
mkdir nf_hd_rnaseq
cd nf_hd_rnaseq
mkdir data
cd data
#download the fastq files into a folder named data.
cd ..
cd ..
#should to be in dir that contains nf_hd_rnaseq
#install docker and run with nf_hd_rnaseq folder as volume
sudo docker run -v $(pwd)/nf_hd_rnaseq:/home/nf_hd_rnaseq -it jackkelly75/nf_hd_rnaseq
#time taken - <4 mins ; size -  3.25GB
cd nf_hd_rnaseq

mkdir 1_FastqPuri
FILES=/home/nf_hd_rnaseq/data/*_1.fastq.gz
for fn in $FILES;
do
	echo "Processing sample $fn"
	a=$(echo ${fn} | sed -e 's/_1.f/_2.f/')
	ln=${fn##*/}
	v2=${ln::-10}
	trimFilterPE -f $fn:$a -l 101 --trimQ ENDSFRAC --trimN ENDS -m 31 -o ./1_FastqPuri/${v2}
	#-m is minumum read length to keep. 
done

```


#### Install Nexflow
Install `nextflow` using [their tutorial](https://www.nextflow.io/docs/latest/getstarted.html)

#### Run NF_HD_RNAseq pipeline
Check if java is installed using `java --version`. \
Can install java using `sudo apt install default-jre`. \
`wget -qO- https://get.nextflow.io | bash` \
The easiest way is to install NF_HD_RNAseq is using nextflow's `pull` command:

```
sudo ./nextflow pull jackkelly75/NF_HD_RNAseq
sudo ./nextflow run jackkelly75/NF_HD_RNAseq --reads '1_FastqPuri/*_{1,2}_good.fq.gz' --transcriptome 'hsapien.fa.gz'
```

with --reads set to where reads are stored on your system, and command run in working directory that contains nextflow (if not in PATH)


### Processes

*buildindex*\
Uses salmon to build index
(Need to have the human transcriptome downloaded already from this bash line 
```
curl ftp://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz -o hsapien.fa.gz
```
Index is created using the human transcriptome downloaded and k-mer length of 31 (largest possible size, optimed for reads over ~75bp).


*quant*\
Uses salmon to quantify the filtered RNA-seq results. Index was created previously  in *buildindex*/
Salmon is run using:
<pre>
salmon quant -i $index -l A \
            -1 ${reads[0]} -2 ${reads[1]} \
            --seqBias  \  #learn and correct for sequence-specific biases in the input data
            --gcBias \  #learn and correct for fragment-level GC biases in the input data. Does not impact on results if GC bias is not present, only marginally increases run time
            --posBias \
            --validateMappings\  #selective alignment that is more sensitive
            -o $pair_id
</pre>



*sortfiles*\
Uses R to sort the salmon quant files and get DEGs.

replace missing pmi data point using KNN
to remove the .*N* at end of the annotation so can be annotated to gene symbol
change ensembl annotation to symbol (any duplicate symbol genes, keep only the one with the highest MAD)
convert transcript to gene controlling for pmi and rin
differential expression using DESeq2 (IHW correction)


### Download phenodata and rscript files
```
curl -o pData.csv https://raw.githubusercontent.com/jackkelly75/NF_HD_RNAseq/master/data/pData.csv
curl -o DESeq.R https://raw.githubusercontent.com/jackkelly75/NF_HD_RNAseq/master/DESeq.R
cd ..
sudo docker run -v $(pwd)/nf_hd_rnaseq:/home/nf_hd_rnaseq \
                -it jackkelly75/nf_hd_rnaseq
cd nf_hd_rnaseq
mkdir 3_DESeq
chmod +x DESeq.R
./DESeq.R
```




### References

Pérez-Rubio P, Lottaz C & Engelmann JC (2019). FastqPuri: high-performance preprocessing of RNA-seq data. *BMC Bioinformatics* 20: Article number 226
Williams CR, Baccarella A, Parrish JZ & Kim CC (2016). Trimming of sequence reads alters RNA-Seq gene expression estimates. *BMC Bioinformatics* 17: Article number 103
