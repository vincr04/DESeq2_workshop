# DESeq2 workshop
December 2022

In this workshop, we will use DESeq2 to analyse differential gene expression in an example RNA-seq dataset. This simplified exercise is inspired by the [original DESeq2 vignette by Love, Anders and Huber](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
***

## 1. Preparation

### Install the software

DESeq2 runs in R and is available as a package via [Bioconductor](bioconductor.org), which is a large-scale project to develop, support, and disseminate open source software for bioinformatic data analysis. Many tools used by computational biologists are available there. 

You should all have [RStudio](https://posit.co/downloads/) installed on your computer. Google 'DESeq2 Bioconductor' to find out how to install DESeq2. Once DESeq2 is installed, you need R to load the DESeq2 package by typing `library("DESeq2")`. *Remember, code is ALWAYS case-sensitive*. Once this is done, move to the next step.

### Get the data

Bioconductor also contains some 'pre-packaged' datasets that can be easily downloaded and used as examples. Today, we will work on the 'pasilla' dataset ([Brooks et al., Genome Research 2011](https://pubmed.ncbi.nlm.nih.gov/20921232/)) which explores the effects of RNAi knockdown of *pasilla*, a nuclear RNA binding protein implicated in splicing and ortholog of mammalian NOVA1 and NOVA2, on the transcriptome of cultured *Drosophila melanogaster* cells.

The data from this experiment is available as an R package. As above, use Google to find out how to install the `pasilla` package from Bioconductor, and once this is done, don't forget to type `library("pasilla")` to load the package and the data. 

### Additional packages to install

We will use [tidyverse](https://www.tidyverse.org/) to manipulate data. Tidyverse is available via CRAN, which means that you can install it from the 'Packages' tab on the top right of your screen in RStudio. Then, load the package. 

[ggplot2](https://ggplot2.tidyverse.org/reference/ggplot.html) might be an old friend of yours. If it's already installed on your computer you can just load the package. Otherwise, you can also obtain `ggplot2` *via* CRAN. Then, load the package. 

<br/>

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

<br/>

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

<br/>

## 4. Running the differential expression analysis

The actual command to run DESeq2 is pretty simple:
```
dds <- DESeq(dds)
```
This can run for up to a minute, depending on your computer. The results are saved in the `dds` object and you can visualise a preview of the results by typing 
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

<br/>

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
p + geom_point(color=ifelse(res$padj<0.05, "red", "grey")) + # adding some colour
    lims(y=c(0,40), x=c(-2.5,2.5)) # cropping X and Y axes

# OPTIONAL: we can also add gene labels, using the `ggrepel` package
install.packages("ggrepel")
library("ggrepel")
p + geom_point(color=ifelse(res$padj<0.05, "red", "grey")) + # adding some colour
    lims(y=c(0,40), x=c(-2.5,2.5)) + # cropping X and Y axes
    geom_text_repel(aes(x=log2FoldChange, y=-log(padj), label=ifelse(rownames(res)=="FBgn0038198", rownames(res), "")))
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
>Can you confirm that the knock down of this gene was indeed successful? Note: you need to search [FlyBase](flybase.org) to find the FlyBase ID (FBgn) for the *pasilla* gene.  
>OPTIONAL - how would you plot this with ggplots?

### Heatmap
Often a heatmap is an ideal way of representing variation in expression across several genes. Here is an example (using `ggplot`) with the top 10 genes with the highest *log2FoldChange* and significant adjusted p-value. 
```
res.sig <- res[!is.na(res$padj) & res$padj<0.05,] # Filtering only significant results
top10genes <- rownames(res.sig[rev(order(res.sig$log2FoldChange)),][1:10,]) # Selecting the genes with highest log2FoldChange
top10counts <- counts(dds,normalized=TRUE)[top10genes,] # extract normalised counts for each gene
top10counts <- as.data.frame(t(scale(t(top10counts)))) # calculate z-scores
top10counts <- cbind(top10counts, Gene=rownames(top10counts)) # add extra column with gene name
top10counts <- gather(top10counts, Treatment, Value, -Gene) # reshape data to use with ggplot

ggplot(top10counts, aes(x=Treatment, y=Gene, fill=Value)) + # create a plot
    geom_tile() + # add tiles (for heatmap)
    scale_fill_gradient2(low="navy", mid="linen", high="darkred", na.value="transparent") # add colour scheme
```
> Can you recreate this heatmap for the 10 most downregulated genes? 
> OPTIONAL - can you understand every part of the code above? Can you improve the process?

### PCA
Principal component analysis can be used to establish how different samples are from each other. Conveniently, the `plotPCA` function is included in DESeq2, but first, we need to transform the raw count data using *variance stabilising transformations* (*VST*), which produces normalised, log2 scale values. 
```
vsd <- vst(dds, blind=FALSE) # perform the variance stabilising transformations. We use blind=FALSE to calculate the 'within group' variabliity.
head(assay(vsd), 3) # print the calculated values for the first three genes
plotPCA(vsd, intgroup=c("condition", "type"))
```
> What do you conclude from this plot?

<br/>

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
>Repeat some of the plots above with the LFC-shrinked data. Can you see any differences? 

<br/>

## 7. Gene-set enrichment analysis (GSEA)

We can now tell which genes are differently regulated in *pasilla* knock-down cells, however this in itself does not provide much biological information. What would be interesting at this stage is to switch the focus of our analysis from individual genes onto biological pathways. For this, we can use gene-set enrichment analysis (GSEA), also called pathway enrichment analysis.  

GSEA uses existing databases with genome-wide information about the characteristics of each gene, which enables grouping these genes into specific categories. The best known example of such database is the Gene Ontology (GO), which categorises (almost) all known genes according to their cellular component, molecular function, and biological process.  

Here, we will use the `clusterProfiler` package to perform GSEA, and the `msigdbr` package, which contains gene sets from the [Molecuar Signatures Database (MSigDB)](https://www.gsea-msigdb.org/gsea/msigdb/index.jsp). We will also need the `org.Dm.eg.db` package to convert *Drosophila* gene identifiers.  

First, install these three packages (from Bioconductor) and load the packages.  

### Generate gene set

The `msigdbr_species()` command tells you what species are available in MSigDB. *Drosophila melanogaster* is one of them. However, within MSigDB there are several gene set collections, all generated differently. The 'hallmark' collection is a good starting point, as it was made from aggregating other gene sets to obtain a comprehensive and coherent gene set. We will therefore generate an R object for this gene set. 
```
Dm_hallmark_sets <- msigdbr(species = "Drosophila melanogaster", category = "H")
head(Dm_hallmark_sets)
```
> Take a moment to look at what information is available in this table.  

### Prepare the gene list

Before performing GSEA, we need to create a vector with sorted *log2FoldChange* values, and the corresponding gene names. 
```
lfc_vector <- res$log2FoldChange
names(lfc_vector) <- rownames(res)
lfc_vector <- sort(lfc_vector, decreasing = TRUE)
```

### Run the GSEA

We use the `GSEA` function on the list of *log2FoldChange* values to 
```
gsea_results <- GSEA(geneList = lfc_vector, pvalueCutoff = 1, TERM2GENE = dplyr::select(Dm_hallmark_sets, gs_name, ensembl_gene))
gsea_result_df <- data.frame(gsea_results@result) # make a data.frame with the results
gsea_result_df %>% arrange(pvalue) %>% select(1:10) %>% head() # visualise the gene sets with the lowest p-values
```
We can see that only one gene set has a significant p-value, however after p-value adjustment this is no longer significant. 
>Does that really mean that there is no significant gene set specifically enriched in our data? Could you get significant enrichment using another MSigDB collection? Or another gene set database?

### Plot the results

Even if there is no significantly enriched gene set, we can always visualise GSEA results for the most enriched ones. 
```
topNES <- gsea_result_df %>% arrange(pvalue) %>% select(ID) %>% head(1) %>% pull() # shows the gene set with the highest NES and lowest p-value
enrichplot::gseaplot(gsea_results, geneSetID = topNES, title = topNES)
```
> How do you interpret this plot?

<br/>

## 8. OPTIONAL - Differential isoform expression

The *pasilla* gene encodes a nuclear RNA binding protein implicated in mRNA splicing. Therefore, we expect to record differential splicing events between the treated and untreated samples. To investigate this, download DEXSeq from Bioconductor, and follow the steps in the [DEXSeq vignette](https://bioconductor.org/packages/release/bioc/vignettes/DEXSeq/inst/doc/DEXSeq.html) (you can start from step 3.1.). 
> How many genes can you find with significantly different exon usage?  
> Read the [original *pasilla* paper](https://genome.cshlp.org/content/21/2/193.long). Does your analysis help you understand the findings from this paper?

<br/>
<br/>

<h1 align="center">The end</h1>
<p align="center">(you are now an expert!)</p>

