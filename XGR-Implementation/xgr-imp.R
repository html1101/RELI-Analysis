source("http://bioconductor.org/biocLite.R")

BiocManager::install(c("hexbin","ape","supraHex","graph","Rgraphviz","igraph","Biobase","limma","survival","foreach","doMC","devtools"))

install.packages("dnet",repos="http://cran.r-project.org",type="source")

BiocManager::install("hfang-bristol/XGR", dependencies=T)

## reload the installed package
detach(package:XGR, unload=T)

library(XGR)

install.packages("igraph")

anno <- xRDataLoader(RData='GWAS2EF', RData.location=RData.location)
allSNPs <- rownames(anno)
data <- sample(allSNPs,8)
data
RData.location <- "http://galahad.well.ox.ac.uk/bigdata/"


sim <- xSocialiserSNPs(data=data, include.LD="EUR", LD.r2=0.8, RData.location=RData.location)

# c) save similarity results to the file called 'EF_similarity.txt'
output <- igraph::get.data.frame(sim, what="edges")
utils::write.table(output, file="EF_similarity.txt", sep="\t",
                   row.names=FALSE)

# d) visualise the SNP network
## extract edge weight (with 2-digit precision)
x <- signif(as.numeric(E(sim)$weight), digits=2)
## rescale into an interval [1,4] as edge width
edge.width <- 1 + (x-min(x))/(max(x)-min(x))*3
## do visualisation
xVisNet(g=sim, vertex.shape="sphere", edge.width=edge.width,
        edge.label=x, edge.label.cex=0.7)

