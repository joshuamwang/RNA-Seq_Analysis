## Introduction
Among the most groundbreaking advancements developed in the 21st century, genome sequencing as a [disruptive technology](https://www.mckinsey.com/business-functions/digital-mckinsey/our-insights/disruptive-technologies) has expanded new fields of scientific thought to probe at the genetic complexity underlying human disease from different perspectives. In this tutorial, we walk through a specific approach to analyze [RNA-Sequencing datasets](https://en.wikipedia.org/wiki/RNA-Seq), with examples and parameters specific to [Candida albicans](https://en.wikipedia.org/wiki/Candida_albicans).  


## Workflow
![alt text](https://github.com/joshuamwang/RNA-Seq_Analysis/blob/master/figure-animated.gif "workflow")

We will largely be focused on identifying differentially expressed (DE) genes, which follows the steps highlighted in red.


## Installation
These installation steps have been designed for and tested on a [Ubuntu 16.04.3 LTS](http://releases.ubuntu.com/16.04/) workstation.

**1\) [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to check read quality:**
```bash
$ wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.7.zip
$ unzip fastqc_v0.11.7.zip 
$ cd FastQC
$ sudo chmod +x fastqc # enter admin password when prompted
$ sudo ln -s /path/to/FastQC/fastqc /usr/bin/fastqc 

# the path to FastQC can be determined as follows:
$ pwd
/home/wang.10983/FastQC # your output will be different
# in this example, my path to FastQC would be /home/wang.10983/FastQC/fastqc
```

**2\) [Spliced Transcripts Alignment to a Reference (STAR)](https://github.com/alexdobin/STAR) to align the FASTQ reads:**
```bash
# get STAR source using git
$ git clone https://github.com/alexdobin/STAR.git
$ cd STAR/source

# Build STAR
$ make STAR
```

**Optional:** After aligning, [Integrated Genome Viewer (IGV)](http://software.broadinstitute.org/software/igv/) can be used to visualize where reads aligned across the genome.

a) Download [IGV](http://software.broadinstitute.org/software/igv/download) to local machine, in addition to [Java 8](https://java.com/en/download/) (or later) if not already locally pre-installed.

b) Install [Samtools](http://www.htslib.org/) to server.
```bash
$ wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 -O samtools.tar.bz2
$ tar -xjvf samtools.tar.bz2
$ cd samtools-{version}
$ make
$ sudo make prefix=/usr/local/bin install
$ sudo ln -s /usr/local/bin/bin/samtoos /usr/bin/samtools
```

**3\) [HTSeq](https://htseq.readthedocs.io/en/master/install.html) to create a counts table:**
```bash
$ sudo apt-get update
$ sudo apt-get install build-essential python2.7-dev python-numpy python-matplotlib python-pysam python-htseq
```

**4\) [Empirical Analysis of Digital Gene Expression Data in R (edgeR)](https://bioconductor.org/packages/release/bioc/html/edgeR.html) to identify differentially expressed genes:**
```bash
$ sudo apt-get update
$ sudo apt-get install r-base
$ sudo su -c "R -e \"source('https://bioconductor.org/biocLite.R'); biocLite('edgeR')\""
```


## Steps to Run

**0\) Materials**

In this tutorial, we will run through a sample Candida albicans RNA Sequencing file which can be downloded from [here](https://drive.google.com/file/d/1AP_xnwmqXobrquf7kHYEOdo0FtI9h7fo/view?usp=sharing).

**1\) Validate Read Quality**

"FastQC aims to provide a simple way to do some quality control checks on raw sequence data coming from high throughput sequencing pipelines. It provides a modular set of analyses which you can use to give a quick impression of whether your data has any problems of which you should be aware before doing any further analysis." [reference](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
```bash
$ fastqc CA_fastq.gz
```
This should output an html file along named "CA_fastqc.html". Download this html file using [Filezilla](https://filezilla-project.org/download.php?type=client) or [Cyberduck](https://cyberduck.io/?l=en) and open the html file locally on your computer using [Chrome](https://www.google.com/chrome/browser/desktop/index.html) or [Firefox](https://www.mozilla.org/en-US/firefox/new/). 
In particular, the per base sequence quality graph reports the quality score across all bases:
<p align="center">
  <img src="https://github.com/joshuamwang/RNA-Seq_Analysis/blob/master/fastqc.png">
</p>

We want to aim for quality scores >= 28. **Low quality read scores require a revisit to wet lab validation. Only read samples that pass FastQC should procede forward.**

**2\) Align Reads to Reference Genome**

We choose the STAR aligner for its balance of computational efficiency and mapping accuracy. "STAR outperforms other aligners by a factor of >50 in mapping speed, aligning to the human genome 550 million 2 × 76 bp paired-end reads per hour on a modest 12-core server, while at the same time improving alignment sensitivity and precision." [reference](https://academic.oup.com/bioinformatics/article/29/1/15/272537)

The reference genome will need to be downloaded separately. Here, we demonstrate the steps to align a Candida Albicans RNA-Sequencing sample to the [SC5314 v21 genome reference](http://candidagenome.org/download/sequence/C_albicans_SC5314/Assembly21/current/).
```bash
# download & extract reference genome
$ mkdir genome
$ cd genome
$ wget http://candidagenome.org/download/sequence/C_albicans_SC5314/Assembly21/current/C_albicans_SC5314_A21_current_chromosomes.fasta.gz
$ gunzip C_albicans_SC5314_A21_current_chromosomes.fasta.gz
$ cd ../
```

The first time STAR is executed, we need to create a genome index as follows:
```bash
$ STAR --runThreadN 12 --runMode genomeGenerate --genomeDir genome --genomeFastaFiles $ genome/C_albicans_SC5314_A21_current_chromosomes.fasta --sjdbGTFtagExonParentTranscript ID
```
a) ---genomeDir = name of the folder containing the genome sequence (in our case, it is just genome)

b) ---genomeFastaFiles = location of the reference fasta file we previously downloded

Next, we can align the fastq file against the reference:
```bash
# aligning to STAR with recommended parameters
$ STAR --runThreadN 12 --outSAMstrandField intronMotif --genomeDir genome --readFilesIn CA.fastq.gz --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --alignIntronMin 30 --alignIntronMax 1000 
```
**Note:** Remove the "--readFilesCommand zcat" flag if the fastq file is not gzipped.

After a few minutes, the aligner should output the following files:
```bash
$ ls
Aligned.sortedByCoord.out.bam
Log.final.out
Log.out
Log.progress.out
SJ.out.tab
```
Log.final.out provides summary statistics of the alignment process, namely the percentage of uniquely mapped reads (93.97% in our example). The Aligned.sortedByCoord.out.bam is the output alignment file. Let's rename it for convenience:
```bash
mv Aligned.sortedByCoord.out.bam CA.bam
```

**Optional:** To visualize the newly created alignment file in [IGV](http://software.broadinstitute.org/software/igv/), we first need to create its index.
```bash
$ samtools index CA.bam
# the terminal will probably hang for a few seconds.
```
a) Transfer the CA.bam, CA.bam.bai, and C_albicans_SC5314_A21_current_chromosomes.fasta files to your local computer. Both files should be stored in the same folder. 

b) Launch [IGV](http://software.broadinstitute.org/software/igv/)

c) Click File from the upper left corner -> Load from file -> select the CA.bam file. 

d) Select Genomes -> Load Genomes from File -> select C_albicans_SC5314_A21_current_chromosomes.fasta

![alt text](https://github.com/joshuamwang/RNA-Seq_Analysis/blob/master/IGV.png "IGV")

**3\) Count Gene Features**
Next, we will use [HTSeq](https://htseq.readthedocs.io/en/master/install.html) to quantify the number of reads that aligned to each annotated feature. Supplying an annotated [.gff file](https://en.wikipedia.org/wiki/General_feature_format) to [HTSeq](https://htseq.readthedocs.io/en/master/install.html) is required.
```bash
# Obtain C. albicans gff file
wget http://candidagenome.org/download/gff/C_albicans_SC5314/Assembly21/C_albicans_SC5314_A21_current_features.gff

# run HTSeq
$ htseq-count -t gene -i ID -f bam CA.bam C_albicans_SC5314_A21_current_features.gff > CA_count.txt
```

**4\) Identify DE Genes**

**Optional:** This section involves executing ```R``` code -- Users are *highly recommended* to use [RStudio](https://www.rstudio.com/products/rstudio/download/#download) to facilitate with the debugging process. 

A gene differential analysis would ideally involve the comparison between at least two different biological conditions with multiple replicates for each. To simulate, let's pretend we have 2 biological conditions labeled A and B, with 2 replicates of each as follows:
  
| Condition/Replicate | Alignment File | Counts File |
| :-------------: |:-------------:|:-------------:|
| A/1 | A_1.bam | A_1_counts.txt |
| A/2 | A_2.bam | A_2_counts.txt |
| B/1 | B_1.bam | B_1_counts.txt |
| B/2 | B_2.bam | B_2_counts.txt |


```R
# The abridged tutorial below is based on a more complete version located [here](https://web.stanford.edu/class/bios221/labs/rnaseq/lab_4_rnaseq.html).
library('edgeR')

compList <- c("A_1_counts.txt","A_2_counts.txt","B_1_counts.txt","B_2_counts.txt")

for(comp in compList){
  comp1 <- unlist(strsplit(comp,split="_"))[1]
  comp2 <- unlist(strsplit(comp,split="_"))[2]
  
  a <- read.table(as.character(compList[1]),row.names=1,col.names=c("Symbol","Counts"),stringsAsFactors=F)
  b <- read.table(as.character(compList[1]),row.names=1,col.names=c("Symbol","Counts"),stringsAsFactors=F)
  d <- read.table(as.character(compList[1]),row.names=1,col.names=c("Symbol","Counts"),stringsAsFactors=F)
  e <- read.table(as.character(compList[1]),row.names=1,col.names=c("Symbol","Counts"),stringsAsFactors=F)

  x <- cbind(a,b[rownames(a),],d[rownames(a),],e[rownames(a),])
  group <- factor(c(1,1,2,2))
  
  y <- DGEList(counts=x,group=group)
  
  keep <- rowSums(cpm(y)>1) >= 2
  y <- y[keep, , keep.lib.sizes=FALSE]
  
  y <- calcNormFactors(y)
  design <- model.matrix(~group)
  y <- estimateCommonDisp(y)
  y <- estimateTagwiseDisp(y)
  
  all <- exactTest(y,pair=c(1,2))
  
  write.csv(data.frame(topTags(all,n=10000)),paste0(comp1,"vs",comp2,"_ALL.csv"))
  temp <- data.frame(topTags(all,sort.by="PValue",p.value=0.05,n=10000))
  temp <- temp[abs(temp$logFC)>=2,]
  
  write.csv(temp,paste0(comp1,"vs",comp2,"_FILTERED.csv"))
}
 ```


