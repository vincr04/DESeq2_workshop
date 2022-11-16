# DESeq2_workshop
December 2022

In this workshop, we will use DESeq2 to analyse differential gene expression in an RNA-seq dataset. This simplified exercise is inspired by the [original DESeq2 vignette by Love, Anders and Huber](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
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

>The count matrix likely contains several thousands of rows (one per gene). What command would you use to display only a few rows?  
>Can you tell how many rows and columns this table has? Use Google if you need. 

We can also visualise the **sample information** (it is a small table so we can just display it entirely). 

It is essential that the column headers of the **counts matrix** correspond to the sample names in the **sample information** table, and are arranged in the same order. We can easily see that this is not the case: the sample names in the **sample information** have an extra 'fb' that needs to be removed, and the columns in the **counts matrix** needs to be rearranged. So we need to edit the sample names
```
rownames(coldata) <- sub("fb", "", rownames(coldata))
```
and arrange the order of the columns in the **counts matrix**
```
cts <- cts[, rownames(coldata)]
```

> What do these commands do?  
> Try and visualise the counts matrix and sample info again. What has changed?

### Build the *DESeqDataSet* 
We can now load the **count matrix** and **sample info** into a new *DESeqDataSet*, which we will call `dds`. For this, we use the *DESeqDataSetFromMatrix* function, which is part of the DESeq2 package:
```
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
```

>Here, we created a new object that contains count data and sample info. Type `dds` to display some info about the new object you created. What kind of information do you get?  
>DESeq2 allows easy access to the key elements of a *DESeqDataSet*. For example, you can type `counts(dds)` to quickly view the **count matrix**.

## 3. Pre-filtering and changing reference levels
The genes with the lowest count numbers are likely to be almost absent from the cell culture that we are analysing. Removing these genes is not 100% essential, however having less gene might speed up the analysis. Here, we only do minimal filtering by removing genes for which we have less than 10 reads across all samples.
```
keep <- rowSums(counts(dds)) >= 10 #gives a TRUE or FALSE value for each gene (row): TRUE if there are more than 10 reads
dds <- dds[keep,] #keep only genes with a TRUE value
```

>How many genes were removed?  

****Note**** by default, conditions are considered in alphabetical order, with the first one being assigned as the 'reference'. So, if for example your conditions are labelled 'untreated' and 'treated', by default the software will assign the 'treated' group as the control, which would be wrong. Conditions are saved as *factors* in R and if you type `dds$condition`, you will notice that indeed 'treated' comes first in the list of levels. Use the following command to make sure that 'untreated' is the first level
```
dds$condition <- relevel(dds$condition, ref = "untreated")
```
>Type `dds$condition` again. What has happened?

## 4. Running the differential expression analysis

The actual command to run DESeq2 is pretty simple:
```
dds <- DESeq(dds)
```
This can take up to a minute, depending on your computer. The results are saved in the `dds` object and you can visualise a preview of the results by typing 
```
res <- results(dds)
res
```

>Take a moment to look at and/or play around with the table in front of you. What information is in the header? Can you make sense of what each column represents?  
>If you are struggling, the command `mcols(res)$description` conveniently provides info about each column.  
>Can you find the gene with the highest/lowest *log2FoldChange*? And the gene with the lowest adjusted p-value? Hint: you can sort the table by the values of selected columns.  
>Why are some *log2FoldChange* values negative?  
>OPTIONAL - to avoid false positives (type I errors), the p-value is corrected with the Benjamini-Hochberg (BH) method. Use Google to find more about the BH method.  

We can add options to the `results` command. For example, `contrast` can be used to use a different reference sample. 

>Try the two following commands. How do the results change?  
>`results(dds, contrast=c("condition","treated","untreated"))`  
>`results(dds, contrast=c("condition","untreated","treated"))`

`summary(res, alpha=0.05)` provides some interesting statistics about the data

## 5. Visualising the data

### MA-plot
DESeq2 offers simple commands to visualise the data. For example, the `plotMA` function displays *log2FoldChange* over the mean of normalised counts, while highlighting genes with a significant adjusted p-value.
```
plotMA(res, ylim=c(-4,4), alpha=0.05)
```
>What does this plot tell you? What happens if you use `alpha=0.01`?  
>OPTIONAL - Can you recreate this plot using ggplots? 

### Volcano plot
A volcano plot is a scatter plot representing the negative log of the *padj* over the *log2FoldChange*. As for the MA-plot, each gene is represented by a point, and a visualisation of the distribution for these two values can be obtained. DESeq2 does not include a function to make volcano plots, but we can easily make one using `ggplot`. 
```
#We can start with a very basic representation
p <- ggplot(as.data.frame(res), aes(x=log2FoldChange, y=-log(padj)))
p + geom_point()

# And make it look a bit nicer
p + geom_point(color=ifelse(res$padj<0.05, "red", "grey")) # adding some colour
  + lims(y=c(0,40), x=c(-2.5,2.5)) # cropping X and Y axes

# OPTIONAL: we can also add gene labels, using the `ggrepel` package
install.packages("ggrepel")
library("ggrepel")
p + geom_point(color=ifelse(res$padj<0.05, "red", "grey")) # adding some colour
  + lims(y=c(0,40), x=c(-2.5,2.5)) # cropping X and Y axes
  + geom_text_repel(aes(x=log2FoldChange, y=-log(padj), label=ifelse(rownames(res)=="FBgn0038198", rownames(res), "")))
```

### P-value histogram
Now let's plot a histogram of the adjusted p-values
```
hist(res$padj,breaks = 100); abline(v=0.05,col="red")
```
>What does this plot tell us?  

### Plot normalised counts for specific genes
Using the following command, you can plot the normalised counts between conditions for any gene you like (just replace *XXX* with the name of your gene of choice). 
```
plotCounts(dds, gene=XXX, intgroup="condition")
```
To plot the gene with the lowest *padj*, you can replace *XXX* with `which.min(res$padj)`. 

>What command would you use to plot the gene with highest *log2FoldChange*?  
>OPTIONAL - how would you plot this with ggplots?

## 6. Log fold change shrinkage
Genes with low counts are more likely to have high *log2FoldChange* values because the natural variation between samples may create artificial differences between samples. For example, in the table below, the *mean* count value for the treated samples is more than 4x lower than for the untreated sample for **Gene A**, but almost equal for **Gene B**, even though the *difference* in read counts between both samples is the same. This leads to a bias towards low-count genes having a high *log2FoldChange*, with **Gene A** having a higher *log2FoldChange* than **Gene B**. 

|Gene  |untreated 1|untreated 2|untreated 3|untreated 4|treated 1|treated 2|treated 3|
|------|-----------|-----------|-----------|-----------|---------|---------|---------|
|Gene A|10         |8          |4          |5          |1        |3        |1        |
|Gene B|1010       |1008       |1004       |1005       |1001     |1003     |1001     |

To correct for this, we can apply a DESeq2 function called `lfcShrink`. This function requires the user to specify a `coef`, which is an argument specifying the comparison to extract. We can use the `resultsNames(dds)` command to display available `coef` options. We will select `"condition_treated_vs_untreated"`. 

We can now run the `lfcShrink` function. The 'type' argument specifies the LFC shrinkage method - in general, the more recent `apeglm` method is recommended.
```
resLFC <- lfcShrink(dds, coef="condition_treated_vs_untreated", type="apeglm")
plotMA(resLFC, ylim=c(-4,4), alpha=0.05)
```
>What has changed between this plot and the previous one? 






