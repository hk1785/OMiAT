# R package: OMiAT

Version: 6.0

Date: 2020-9-26

Title: Optimal Microbiome-based Association Test (OMiAT)

Author: Hyunwook Koh

Maintainer: Hyunwook Koh <hyunwook.koh@stonybrook.edu>

Description: This software package provides facilities for 1) optimal microbiome-based association test (OMiAT) which tests the association between the entire community (e.g., kingdom) or indivisual upper-level taxa (e.g., phylum, class, order, family, genus) and a host phenotype of interest (or disease status). 2) microbiome comprehensive association mapping (MiCAM) which tests the association for all microbial taxa through a breadth of taxonomic levels. Continuous/binary responses (e.g., BMI, disease status) can be surveyed with/without covariate adjustments (e.g., age, gender).

NeedsCompilation: No

Depends: R(>= 3.5.3)

Imports: ape, BiasedUrn, cluster, CompQuadForm, dirmult, ecodist, GUniFrac, phangorn, phyloseq, robustbase, robCompositions

License: GPL-2

## Reference

* Koh H, Blaser MJ, Li H. (2017) A powerful microbiome-based association test and a microbial taxa discovery framework for comprehensive association mapping. _Microbiome_ 5:45.

* DOI: https://doi.org/10.1186/s40168-017-0262-x

## Troubleshooting Tips

If you have any problems for using this R package, please report in Issues (https://github.com/hk1785/OMiAT/issues) or email Hyunwook Koh (hkoh@jhu.edu).

* Tip 1. Depending on your pre-installed R libraries, this R package can require you to install additional R packages such as "gh", "usethis", "cli", etc using the command: install.packages("package_name").
* Tip 2. Please make sure if you have the most recent package version.

## Prerequites

ape
```
install.packages("ape")
```
BiasedUrn
```
install.packages("BiasedUrn")
```
cluster
```
install.packages("cluster")
```
CompQuadForm
```
install.packages("CompQuadForm")
```
dirmult
```
install.packages("dirmult")
```
ecodist
```
install.packages("ecodist")
```
GUniFrac
```
install.packages("GUniFrac")
```
phangorn
```
install.packages("phangorn")
```
robustbase
```
install.packages("robustbase")
```
robCompositions
```
install.packages("robCompositions")
```
phyloseq
```
source("https://bioconductor.org/biocLite.R")
biocLite("phyloseq")
```

## Installation
```
library(devtools)
install_github("hk1785/OMiAT", force=T)
```

## Data format
```
library(phyloseq)
```
* URL: https://joey711.github.io/phyloseq/

---------------------------------------------------------------------------------------------------------------------------------------

# Manual
This R package includes two core functions, OMiAT and MiCAM. Please find the details below.

## :mag: OMiAT

### Description
OMiAT tests the association between the microbial composition and a host phenotype of interest (or disease status) with/without covariate adjustments (e.g., age, gender). For the microbial composition, the entire community (e.g., kingdom) or indivisual upper-level taxa (e.g., phylum, class, order, family, genus) can be surveyed.

### Usage
```
OMiAT(Y, otu.tab, cov = NULL, tree, total.reads = NULL, model = c("gaussian", "binomial"), pow = c(1:4, Inf), g.unif.alpha = c(0.5), n.perm = 3000)
```

### Arguments
* _Y_ - A numeric vector for continuous or binary responses (e.g., BMI, disease status). 
* _otu.tab_ - A matrix of the OTU table. Notice 1: rows are subjects and columns are OTUs. Notice 2: Monotone/singletone OTUs need to be removed. 
* _cov_ - A data frame for covariate adjustment(s) (e.g., age, gender). Notice: rows are subjects and columns are covariate variables. Default is NULL. 
* _tree_ - A rooted phylogenetic tree. 
* _total.reads_ - A numeric vector for total reads per sample in the entire community. If you survey the entire community, you do not need to specify this numeric vector. If you test an upper-level taxon, you need to specify this arguments additionally. See the example below. Default is NULL. 
* _model_ - "gaussian" is for the linear regression model and "binomial" is for the logistic regression model.
* _pow_ - A set of candidate gamma values. Default is c(1:4, Inf).
* _g.unif.alpha_ - A set of alpha parameter values for generalized UniFrac distance (e.g., c(0.25, 0.5)). Default is c(0.5). 
* _n.perm_ - A number of permutations. Default is 3000.

### Values
_$SPU.pvs_ - A vector of the p-values for individual SPU and Adaptive SPU (aSPU) tests.

_$MiRKAT.pvs_ - A vector of the p-values for individual MiRKAT and Optimal MiRKAT (Opt.MiRKAT) tests.

_$OMiAT.pvalue_ - The p-value for OMiAT.

### References
* Koh et al. (2017) A powerful microbiome-based association test and a microbial taxa discovery framework for comprehensive association mapping. Microbiome. 5:45.

* Pan et al. (2014) A powerful and adaptive association test for rare variants. Genetics. 197(4);1081-95.

* Zhao et al. (2015) Testing in microbiome-profiling studies with MiRKAT, the microbiome regression-based kernel association test. American Journal of Human Genetics. 96(5);797-807.

### Example
Import requisite R packages
```
library(ape)
library(BiasedUrn)
library(cluster)
library(CompQuadForm)
library(dirmult)
library(ecodist)
library(GUniFrac)
library(phangorn)
library(phyloseq)
library(robustbase)
library(robCompositions)
library(OMiAT)
```
Import example microbiome data (phylseq format)
```
data(MiData)
otu.tab <- otu_table(MiData)
tax.tab <- tax_table(MiData)
tree <- phy_tree(MiData)
y.con <- sample_data(MiData)[[1]]
y.bin <- sample_data(MiData)[[2]]
x1 <- sample_data(MiData)[[3]]
x2 <- sample_data(MiData)[[4]]
cov <- as.data.frame(cbind(x1, x2))
cov[,1] <- as.factor(cov[,1])
```

To test the entire community
```
set.seed(123)
OMiAT(Y=y.con, otu.tab=otu.tab, tree=tree, cov=cov, model="gaussian")

set.seed(123)
OMiAT(Y=y.bin, otu.tab=otu.tab, tree=tree, cov=cov, model="binomial")
```

To test the upper-level phylum, Firmicutes
```
# Notice: Create total read counts in the entire community.
total.reads <- rowSums(otu.tab)

ind.Firmicutes <- which(tax.tab[,2] == "p__Firmicutes")
otu.tab.Firmicutes <- otu.tab[,ind.Firmicutes]
tree.Firmicutes <- prune_taxa(colnames(otu.tab.Firmicutes), tree)

set.seed(123)
OMiAT(Y=y.con, otu.tab=otu.tab.Firmicutes, total.reads=total.reads, 
tree=tree.Firmicutes, cov=cov, model="gaussian")

set.seed(123)
OMiAT(Y=y.bin, otu.tab=otu.tab.Firmicutes, total.reads=total.reads, 
tree=tree.Firmicutes, cov=cov, model="binomial")
```

## :mag: MiCAM

### Description
MiCAM tests the association between the microbial composition and a host phenotype of interest (or disease status) with/without covariate adjustments (e.g., age, gender). For the microbial composition, all the microbial taxa through a breadth of taxonomic levels (e.g., kingdom, phylum, class, order, family, genus, species) are surveyed.

### Usage
```
MiCAM(Y, otu.tab, cov=NULL, tax.tab, tree, model = c("gaussian", "binomial"), tax.levels = c("phylum to genus","kingdom to species"), multi.corr = c("BY", "BH"), pow = c(1:4, Inf), g.unif.alpha = c(0.5), n.perm = 3000, ...)
```

### Arguments
* _Y_ - A numeric vector for continuous or binary responses (e.g., BMI, disease status). 
* _otu.tab_ - A matrix of the OTU table. Notice 1: rows are subjects and columns are OTUs. Notice 2: Monotone/singletone OTUs need to be removed. 
* _cov_ - A data frame for covariate adjustment(s) (e.g., age, gender). Notice: rows are subjects and columns are covariate variables. Default is NULL. 
* _tree_ - A rooted phylogenetic tree. 
* _total.reads_ - A numeric vector for total reads per sample in the entire community. If you survey the entire community, you do not need to specify this numeric vector. If you test an upper-level taxon, you need to specify this arguments additionally. See the example below. Default is NULL. 
* _model_ - "gaussian" is for the linear regression model and "binomial" is for the logistic regression model.
* _pow_ - A set of candidate gamma values. Default is c(1:4, Inf).
* _g.unif.alpha_ - A set of alpha parameter values for generalized UniFrac distance (e.g., c(0.25, 0.5)). Default is c(0.5). 
* _n.perm_ - A number of permutations. Default is 3000.
* _tax.levels_ - "phylum to genus" tests from phylum to genus and "kingdom to species" tests from kingdom to species. "phylum to genus" is recommended for 16S data.
* _multi.corr_ - A procedure for multiple testing correction. "BY" is for the Benjamini-Yekutieli procedure, Benjamini & Yekutieli (2001), and "BH" is for the Benjamini-Hochberg procedure, Benjamini & Hochberg (1995).  
* _filename.rel.abundance_ - A file name for the table of relative abundances.
* _filename.unadj.pvs_ - A file name for the table of unadjusted p-values.
* _filename.adj.pvs_ - A file name for the table of adjusted p-values.
* _filename.sig.assoc.taxa_ - A file name for the table of significantly associated taxonomic names.
* _filename.graph_ - A file name for the hierarchical visualization.

### Values
Five output files will be exported to your current working directory (i.e., getwd()). Refer to the arguments "filename.rel.abundance", "filename.unadj.pvs", "filename.adj.pvs", "filename.sig.assoc.taxa", and "filename.graph").

### References
* Koh et al. (2017) A powerful microbiome-based association test and a microbial taxa discovery framework for comprehensive association mapping. Microbiome. 5:45.

* Pan et al. (2014) A powerful and adaptive association test for rare variants. Genetics. 197(4);1081-95.

* Zhao et al. (2015) Testing in microbiome-profiling studies with MiRKAT, the microbiome regression-based kernel association test. American Journal of Human Genetics. 96(5);797-807.

<Multiple testing correction>

* Benjamini & Yekutieli (2001). The control of the false discovery rate in multiple testing under dependency. Annals of Statistics. 29:1165–1188. 

* Benjamini & Hochberg (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society Series B. 57:289–300. 

### Example
Import requisite R packages
```
library(ape)
library(BiasedUrn)
library(cluster)
library(CompQuadForm)
library(dirmult)
library(ecodist)
library(GUniFrac)
library(phangorn)
library(phyloseq)
library(robustbase)
library(robCompositions)
library(OMiAT)
```
Import example microbiome data (phylseq format)
```
data(MiData)
otu.tab <- otu_table(MiData)
tax.tab <- tax_table(MiData)
tree <- phy_tree(MiData)
y.con <- sample_data(MiData)[[1]]
y.bin <- sample_data(MiData)[[2]]
x1 <- sample_data(MiData)[[3]]
x2 <- sample_data(MiData)[[4]]
cov <- as.data.frame(cbind(x1, x2))
cov[,1] <- as.factor(cov[,1])
```

To test from phylum to genus
```
set.seed(123)
MiCAM(Y=y.con, otu.tab=otu.tab, cov=cov, tax.tab=tax.tab, 
tree=tree, model="gaussian", tax.levels="phylum to genus", 
multi.corr="BH", 
filename.rel.abundance="pg.con.rel.abundance.txt", 
filename.unadj.pvs="pg.con.unadj.pvalues.txt", 
filename.adj.pvs="pg.con.adj.pvalues.txt", 
filename.sig.assoc.taxa="pg.con.sig.assoc.taxa.txt", 
filename.graph="pg.con.graph.pdf")

set.seed(123)
MiCAM(Y=y.bin, otu.tab=otu.tab, cov=cov, tax.tab=tax.tab, 
tree=tree, model="binomial", tax.levels="phylum to genus", 
multi.corr="BH", 
filename.rel.abundance="pg.bin.rel.abundance.txt", 
filename.unadj.pvs="pg.bin.unadj.pvalues.txt", 
filename.adj.pvs="pg.bin.adj.pvalues.txt", 
filename.sig.assoc.taxa="pg.bin.sig.assoc.taxa.txt", 
filename.graph="pg.bin.graph.pdf")
```
