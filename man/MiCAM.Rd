\name{MiCAM}
\alias{MiCAM}
\title{
Microbiome Comprehensive Association Mapping (MiCAM)
}

\description{
MiCAM tests the association between the microbial composition and a host phenotype of interest (or disease status) with/without covariate adjustments (e.g., age, gender). For the microbial composition, all the microbial taxa through a breadth of taxonomic levels (e.g., kingdom, phylum, class, order, family, genus, species) are surveyed.
}

\usage{
MiCAM(Y, otu.tab, cov=NULL, tax.tab, tree, 
model = c("gaussian", "binomial"), 
tax.levels = c("phylum to genus","kingdom to species"), 
multi.corr = c("BY", "BH"), pow = c(1:4, Inf), 
g.unif.alpha = c(0.5), n.perm = 3000, ...)
}

\arguments{
  \item{Y}{
A numeric vector for continuous or binary responses (e.g., BMI, disease status). 
}
  \item{otu.tab}{
A matrix of the OTU table. Notice 1: rows are subjects and columns are OTUs. Notice 2: Monotone/singletone OTUs need to be removed. 
}
  \item{cov}{
A data frame for covariate adjustment(s) (e.g., age, gender). Notice: rows are subjects and columns are covariate variables. 
}
  \item{tax.tab}{
A matrix of the taxonomic table. Notice 1: rows are subjects. Notice 2: columns are 7 taxonomic levels (i.e., Kingdom, Phylum, Class, Order, Family, Genus, Species). 7 columns are required. 
}
  \item{tree}{
A rooted phylogenetic tree. 
}
  \item{model}{
"gaussian" is for the linear regression model and "binomial" is for the logistic regression model.
}
  \item{tax.levels}{
"phylum to genus" tests from phylum to genus and "kingdom to species" tests from kingdom to species. "phylum to genus" is recommended for 16S data.
}
  \item{multi.corr}{
A procedure for multiple testing correction. "BY" is for the Benjamini-Yekutieli procedure, Benjamini & Yekutieli (2001), and "BH" is for the Benjamini-Hochberg procedure, Benjamini & Hochberg (1995).  
}
  \item{pow}{
A set of candidate gamma values. Default is c(1:4, Inf).
}
  \item{g.unif.alpha}{
A set of alpha values for the generalized UniFrac distances to be used. Default is c(0.5).
}
  \item{n.perm}{
The number of permutations. Default is 3000.
}
  \item{filename.rel.abundance}{
A file name for the table of relative abundances.
}
  \item{filename.unadj.pvs}{
A file name for the table of unadjusted p-values.
}
  \item{filename.adj.pvs}{
A file name for the table of adjusted p-values.
}
  \item{filename.sig.assoc.taxa}{
A file name for the table of significantly associated taxonomic names.
}
  \item{filename.graph}{
A file name for the hierarchical visualization.
}
}

\value{
The 5 outcome files will be exported to your current working directory (i.e., getwd()). Refer to the arguments "filename.rel.abundance", 
"filename.unadj.pvs", "filename.adj.pvs", "filename.sig.assoc.taxa", and "filename.graph").
}

\references{
Koh et al. (2017) A powerful microbiome-based association test and a microbial taxa discovery framework for comprehensive association mapping. Microbiome. 5:45.

Pan et al. (2014) A powerful and adaptive association test for rare variants. Genetics. 197(4);1081-95.

Zhao et al. (2015) Testing in microbiome-profiling studies with MiRKAT, the microbiome regression-based kernel association test. American Journal of Human Genetics. 96(5);797-807.


Multiple testing correction:


Benjamini & Yekutieli (2001). The control of the false discovery rate in multiple testing under dependency. Annals of Statistics. 29:1165–1188. 

Benjamini & Hochberg (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society Series B. 57:289–300. 
}

\author{
Hyunwook Koh
}

\examples{
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

### An example data 
#(phyloseq format: URL: https://joey711.github.io/phyloseq/)

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

### To test from phylum to genus

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

}
