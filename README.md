# NF_HD_RNAseq
RNA sequencing analysis pipeline for HD data reproducibility

## Introduction
The RNA-seq data (GEO identifier: GSE64810)

3 steps
* filtering of samples using FastqPuri
* quantification using salmon
* identifying DEGs using DESeq2


## Analysing data

### Download data

Create directory named *nf_hd_rnaseq*
```
mkdir nf_hd_rnaseq
```
Download the fastq files into this folder.
Paired data was downloaded from https://www.ebi.ac.uk/ena/data/view/PRJNA271929 for this study.


### Data pre-processing

#### Install docker
This link describes it well:
https://www.digitalocean.com/community/tutorials/how-to-install-and-use-docker-on-ubuntu-18-04

#### Use FastqPuri to preprocess reads
*FastqPuri*\
\
FastqPuri (https://github.com/jengelmann/FastqPuri) used to filter data (Pérez-Rubio *et al*., 2019). trimFilterPE is called to filter the paired end data. --trim ENDSFRAC is called to remove low quality base callings at the end and beggining of reads until the read is above a quality threshold. Accept the trimmed read if the number of low quality nucleotides does not exceed 5%. --trimN ENDS is called to trim N's if found at ends of reads. If the trimmed read length is smaller than the -m set, it is discarded. Trimming and minimum read length are used together as recommended by Williams *et al*., 2016. -l is set to 101 as is read length.

If desired Qreports can be used to create quality reports in HTML. The RNA-seq data has 48 tiles, however 58 samples are run on two lanes and so with the Qreports call these 58 samples must be run with -t (number of tiles) of 96. The results of this were run previously and supplied in case needed rather than run through NextFlow.

FastqPuri does not work well through NextFlow so is run through docker instead.

```
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

### Quantifying reads using Salmon
This nextflow pipeline uses the pre-processed data from above and quantifies it using salmon. It is composed of two steps.

*buildindex*\
Uses salmon to build index
Index is created using the human transcriptome downloaded and k-mer length of 31 (largest possible size, optimed for reads over ~75bp).

*quant*\
Uses salmon to quantify the filtered RNA-seq results. Index was created previously  in *buildindex*/
Salmon is run using:
<pre>
salmon quant -i $index -l A \
            -1 ${reads[0]} -2 ${reads[1]} \
            --seqBias  \  #learn and correct for sequence-specific biases in the input data
            --gcBias \  #learn and correct for fragment-level GC biases in the input data. Does not impact on results if GC bias is not present, only marginally increases run time
            --validateMappings\  #selective alignment that is more sensitive
            -o $pair_id
</pre>


#### Install Nexflow
Check if java is installed using `java --version`. \
If not, can install java using `sudo apt install default-jre`. \
\
Ensure you are within your nf_hd_rnaseq folder when installing nextflow (however if not the executable can be moved to this folder when downloaded).
Install `nextflow` using [their tutorial](https://www.nextflow.io/docs/latest/getstarted.html)  \


#### Run NF_HD_RNAseq pipeline
Download the human reference transcriptome into the nf_hd_rnaseq folder \
```
curl ftp://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz -o hsapien.fa.gz
```
\
The easiest way is to install and run NF_HD_RNAseq is using nextflow's `pull` command:
```
sudo ./nextflow pull jackkelly75/NF_HD_RNAseq
sudo ./nextflow run jackkelly75/NF_HD_RNAseq --reads '1_FastqPuri/*_{1,2}_good.fq.gz' --transcriptome 'hsapien.fa.gz'
```
`--reads` is set to where reads are stored on your system (which should be in the 1_FastqPuri folder)


### Identifying DEGs using DESeq2

Uses R to sort the salmon quant files and get DEGs. \
The main steps are:
* Replace missing pmi data point using KNN
* bin the ages so can be controlled for in model
* remove the .*N* at end of the annotation so can be annotated to gene symbol
* change ensembl annotation to symbol (any duplicate symbol genes, keep only the one with the highest MAD)
* convert transcript to gene controlling for pmi and rin
* differential expression using DESeq2 (IHW correction)

#### Download phenodata and Rscript files and run
Download the phenodata and rscript into the nf_hd_rnaseq folder. Return to the parent directory of nf_hd_rnaseq so it can be included as a volume in the docker.
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
