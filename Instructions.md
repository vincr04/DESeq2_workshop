# DESeq2_workshop
December 2022

In this workshop, we will use DESeq2 to analyse differential gene expression in an RNA-seq dataset. 
***

## 1. Preparation

### Install the software

DESeq2 runs in R and is available as a package via [Bioconductor](bioconductor.org), which is a large-scale project to develop, support, and disseminate open source software for bioinformatic data analysis. Many tools used by computational biologists are available there. 

You should all have [RStudio](https://posit.co/downloads/) installed on your computer. Google 'DESeq2 Bioconductor' to find out how to install DESeq2. Once DESeq2 is installed, you need R to load the DESeq2 package by typing `library("DESeq2")`. *Remember, code is ALWAYS case-sensitive*. Once this is done, move to the next step.

### Get the data

Bioconductor also contains some 'pre-packaged' datasets that can be easily downloaded and used as examples. Today, we will work on the 'pasilla' dataset ([Brooks et al., Genome Research 2011](https://pubmed.ncbi.nlm.nih.gov/20921232/)) which explores the effects of RNAi knockdown of *pasilla*, a nuclear RNA binding protein implicated in splicing and ortholog of mammalian NOVA1 and NOVA2, on the transcriptome of cultured *Drosophila melanogaster* cells.

As above, use Google to find out how to install the 'pasilla' dataset, and once this is done, don't forget to type `library("pasilla")` to load the package and the data. 

### Additional packages required

We will use [tidyverse](https://www.tidyverse.org/) to manipulate datasets. Tidyverse is available via CRAN, which means that you can install it from the 'Packages' tab on the top right of your screen in RStudio. Then, load the package. 

[ggplot2](https://ggplot2.tidyverse.org/reference/ggplot.html) might be an old friend of yours. If it's already installed on your computer you can just load the package. Otherwise, you can also obtain ggplot2 *via* CRAN. 

## 2. Generate a *DESeqDataSet*

### Extract the data
A *DESeqDataSet* is an object class in R that stores the read counts and other metadata, as well as the results of intermediate calculations in a single R object. We need to create a new *DESeqDataSet* from the 'pasilla' data. We first create an object containing the **counts data as a matrix** from the 'pasilla' package that we have previously installed and loaded:
```
pasCts <- system.file("extdata",
                      "pasilla_gene_counts.tsv",
                      package="pasilla", mustWork=TRUE)
cts <- as.matrix(read.csv(pasCts,sep="\t",row.names="gene_id"))
```
And we also create a separate object with the **sample information**:
```
pasAnno <- system.file("extdata",
                       "pasilla_sample_annotation.csv",
                       package="pasilla", mustWork=TRUE)
coldata <- read.csv(pasAnno, row.names=1)
coldata <- coldata[,c("condition","type")]
```
We can now visualise the **count matrix**, to make sure it actually contains our data.

>The count matrix likely contains several thousands of rows (one per gene). What command would you use to display only a few rows? Can you tell how many rows and columns this table has? Use Google if you need. 

We can also visualise the **sample information** (it is a small table so we can just display it entirely). 

It is essential that the column headers of the **counts matrix** correspond to the sample names in the **sample information** table, and are arranged in the same order. We can easily see that this is not the case: the sample names in the **sample information** have an extra 'fb' that needs to be removed, and the columns in the **counts matrix** needs to be rearranged. So we need to edit the sample names
```
rownames(coldata) <- sub("fb", "", rownames(coldata))
```
and arrange the order of the columns in the **counts matrix**
```
cts <- cts[, rownames(coldata)]
```

> What do these commands do? Try and visualise the counts matrix and sample info again. What has changed?

### Build the *DESeqDataSet* 
We can now load the **count matrix** and **sample info** into a new *DESeqDataSet*, which we will call `dds`. For this, we use the *DESeqDataSetFromMatrix* function, which is part of the DESeq2 package:
```
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
```

>Here, we created a new object that contains count data and sample info. Type `dds` to display some info about the new object you created. What kind of information do you get? 
>DESeq2 allows easy access to the key elements of a *DESeqDataSet*. For example, you can type `counts(dds)` to quickly view the **count matrix**.

## 3. Pre-filtering
The genes with the lowest count numbers are likely to be almost absent from the cell culture that we are analysing. Removing these genes is not 100% essential, however having less gene might speed up the analysis. Here, we only do minimal filtering by removing genes for which we have less than 10 reads across all samples.
```
keep <- rowSums(counts(dds)) >= 10 #gives a TRUE or FALSE value for each gene (row)
dds <- dds[keep,]
```

>How many genes were removed?

## 4. Running the differential expression analysis






