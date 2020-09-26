\name{OMiAT}
\alias{OMiAT}
\title{
Optimal Microbiome-based Association Test (OMiAT)
}
\description{
OMiAT tests the association between the microbial composition and a host phenotype of interest (or disease status) with/without covariate adjustments (e.g., age, gender). For the microbial composition, the entire community (e.g., kingdom) or indivisual upper-level taxa (e.g., phylum, class, order, family, genus) can be surveyed.
}
\usage{
OMiAT(Y, otu.tab, cov = NULL, tree, total.reads = NULL, 
model = c("gaussian", "binomial"), pow = c(1:4, Inf), 
g.unif.alpha = c(0.5), n.perm = 3000)
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
  \item{tree}{
A rooted phylogenetic tree. 
}
  \item{total.reads}{
A numeric vector for total reads per sample in the entire community. If you survey the entire community, you do not need to specify this numeric vector. If you test an upper-level taxon, you need to specify this arguments additionally. See the example below. Default is NULL. 
}
  \item{model}{
"gaussian" is for the linear regression model and "binomial" is for the logistic regression model.
}
  \item{pow}{
A set of candidate gamma values. Default is c(1:4, Inf).
}
  \item{g.unif.alpha}{
A set of alpha parameter values for generalized UniFrac distance (e.g., c(0.25, 0.5)). Default is c(0.5). 
}
  \item{n.perm}{
The number of permutations. Default is 3000.
}
}

\value{
SPU.pvs: A vector of the p-values for individual SPU and Adaptive SPU (aSPU) tests.

MiRKAT.pvs: A vector of the p-values for individual MiRKAT and Optimal MiRKAT (Opt.MiRKAT) tests.

OMiAT.pvalue: The p-value for OMiAT.
}

\references{
Koh et al. (2017) A powerful microbiome-based association test and a microbial taxa discovery framework for comprehensive association mapping. Microbiome. 5:45.

Pan et al. (2014) A powerful and adaptive association test for rare variants. Genetics. 197(4);1081-95.

Zhao et al. (2015) Testing in microbiome-profiling studies with MiRKAT, the microbiome regression-based kernel association test. American Journal of Human Genetics. 96(5);797-807.
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

### To test the entire community

set.seed(123)
OMiAT(Y=y.con, otu.tab=otu.tab, tree=tree, cov=cov, model="gaussian")

set.seed(123)
OMiAT(Y=y.bin, otu.tab=otu.tab, tree=tree, cov=cov, model="binomial")

### To test the upper-level phylum, Firmicutes

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

}
