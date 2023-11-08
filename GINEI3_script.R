library(GENIE3)
set.seed(123) # For reproducibility of results

exprMatr <- matrix(sample(1:10, 100, replace=TRUE), nrow=20)
rownames(exprMatr) <- paste("Gene", 1:20, sep="")
colnames(exprMatr) <- paste("Sample", 1:5, sep="")
head(exprMatr)

weightMat <- GENIE3(exprMatr)
dim(weightMat)
weightMat[1:5,1:5]

# Genes that are used as candidate regulators
regulators <- c(2, 4, 7)
# Or alternatively:
regulators <- c("Gene2", "Gene4", "Gene7")
weightMat <- GENIE3(exprMatr, regulators=regulators)

regulatorsList <- list("Gene1"=rownames(exprMatr)[1:10],
                       "Gene2"=rownames(exprMatr)[10:20],
                       "Gene20"=rownames(exprMatr)[15:20])
weightList <- GENIE3(exprMatr, nCores=1, targets=names(regulatorsList), regulators=regulatorsList, returnMatrix=FALSE)
# Use Extra-Trees (ET) method
# 7 randomly chosen candidate regulators at each node of a tree
# 5 trees per ensemble
weightMat <- GENIE3(exprMatr, treeMethod="ET", K=7, nTrees=50)

# Get the list of the regulatory links

linkList <- getLinkList(weightMat)
dim(linkList)
head(linkList)

# Get only the top-ranked links
linkList <- getLinkList(weightMat, reportMax=5)

# Get only the links with a weight higher than some threshold
linkList <- getLinkList(weightMat, threshold=0.1)