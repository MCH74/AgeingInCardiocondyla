##obtain datExpr in main script "WGCNA_script.R"


##remove genes with expression less than 10 in > 12 samples (90%)
datExpr_red <- datExpr[ ,colSums(datExpr >= 10) >= 2 ]

sft = pickSoftThreshold(datExpr_red, powerVector = 1:20, verbose = 5, networkType = "signed", corFnc = "bicor", corOptions = list(maxPOutliers = 0.05) )
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

##pick 14 as power
softPower = 14
adjacency = adjacency(datExpr_red, power = softPower, type = "signed", corFnc = "bicor", corOptions = list(maxPOutliers = 0.05));
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM
# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
#sizeGrWindow(12,9)
#plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04);

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
deepSplit = 2, pamRespectsDendro = FALSE,
minClusterSize = minModuleSize);
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
#table(dynamicColors)
# Plot the dendrogram and colors underneath
#sizeGrWindow(8,6)
#plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05,main = "Gene dendrogram and module colors")

##merging of modules that are similar in expression profiles

# Calculate eigengenes
MEList = moduleEigengenes(datExpr_red, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
xlab = "", sub = "")

##set cut-off
MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr_red, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;


sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
c("Dynamic Tree Cut", "Merged dynamic"),
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)
#dev.off()

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;


###relate modules to traits

# Define numbers of genes and samples
nGenes = ncol(datExpr_red);
nSamples = nrow(datExpr_red);
# Recalculate MEs with color labels
moduleTraitCor = cor(MEs, cobs.design[,4], use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

pdf("cobs_signed_module_traits.pdf")
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3), cex=0.5);
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
xLabels = c("age"),
yLabels = names(MEs),
ySymbols = names(MEs),
colorLabels = FALSE,
colors = greenWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.5,
zlim = c(-1,1),
main = paste("Module-trait relationships Cobs"))

dev.off()


###summarise gene info in a table

geneTraitSignificance = as.data.frame(cor(datExpr_red, cobs.design[,4], use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", "Age", sep="");
names(GSPvalue) = paste("p.GS.", "Age", sep="");


geneInfo_signed = data.frame(genes = colnames(datExpr_red),
                    moduleColor = moduleColors,
                    geneTraitSignificance,
                    GSPvalue)


modOrder = order(-abs(cor(MEs, cobs.design[,4], use = "p")));
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(signedKME(datExpr_red,  MEs, corFnc = "bicor", corOptions = list(maxPOutliers = 0.05)));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

for (mod in 1:ncol(geneModuleMembership))
{
oldNames = names(geneInfo_signed)
geneInfo_signed = data.frame(geneInfo_signed, geneModuleMembership[, modOrder[mod]],
MMPvalue[, modOrder[mod]]);
names(geneInfo_signed) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
paste("p.MM.", modNames[modOrder[mod]], sep=""))
}

geneOrder = order(geneInfo_signed$moduleColor, -abs(geneInfo_signed$GS.Age));
geneInfo_signed = geneInfo_signed[geneOrder, ]

write.table(geneInfo_signed, "geneInfo_cobs_signed.tab", sep="\t",quote=F, row.names=F)


##calculate connectivity
cobs.connectivity <- softConnectivity( datExpr_red, power = softPower, type = "signed", corFnc = "bicor", corOptions = list(maxPOutliers = 0.05)  )
geneInfo_signed <- cbind(geneInfo_signed, connectivity = cobs.connectivity)


##intramodular connectivity
kIM <- intramodularConnectivity(adjacency, moduleColors, scaleByMax = TRUE) 
geneInfo_signed <- cbind(geneInfo_signed, IntraConnectivity = kIM[,2])



##plot modules against changein connectivity
pdf("module_conn_change.pdf",width=14,height=7)
plot(module_traits[,3],module_traits[,1],cex=module_traits[,4]/80,bg=rownames(module_traits),lwd=2,pch=21, xlab="log2 fold-change in connectivity (old versus young)",ylab="Module correlation with age")
points(module_traits[ module_traits[,2] < 0.05,3],module_traits[ module_traits[,2] < 0.05,1],cex=module_traits[ module_traits[,2] < 0.05,4]/80,col="green",bg = rownames(module_traits[ module_traits[,2] < 0.05,]),pch=21,lwd=3)
abline(h=0)
abline(v=0)
dev.off()


##for reduced data set - only genes with at least one TOM of 0.3, all modules except blue
TOM_red <- TOM[rowSums(TOM >= 0.1) > 1, colSums(TOM >= 0.1) > 1]
modules = c("paleturquoise", "red", "greenyellow", "grey60", "darkred", "cyan", "magenta", "violet")#, "yellowgreen")
probes = colnames(TOM_red)
moduleColor_red <- as.character(geneInfo_signed[ geneInfo_signed$genes %in% colnames(TOM_red),"moduleColor"])
inModule = is.finite(match(moduleColor_red, modules));
modProbes = probes[inModule];
modTOM = TOM_red[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,edgeFile = "Cytoscape_edges.txt", nodeFile = "Cytoscape_nodes.txt",weighted = TRUE, nodeNames = modProbes, nodeAttr = moduleColor_red[inModule], threshold=0.01)

##Show network only for top hubs - IntraConnectivity > 0.5 - doesn't work, destroys the network!
modules = c("paleturquoise", "red", "greenyellow", "grey60", "darkred", "cyan", "magenta", "violet", "yellowgreen")
probes = geneInfo_signed[ geneInfo_signed$moduleColor %in% modules & geneInfo_signed$IntraConnectivity > 0.4,1]
moduleColor_red <- as.character(geneInfo_signed[ geneInfo_signed$genes %in% probes,"moduleColor"])
modTOM = TOM[probes, probes];
cyt = exportNetworkToCytoscape(modTOM,edgeFile = "Cytoscape_edges_hubs.txt", nodeFile = "Cytoscape_nodes_hubs.txt",weighted = TRUE, nodeNames = probes, nodeAttr = moduleColor_red, threshold=0.01)

##for reduced data set - all modules but reduce further with higher threshold
TOM_red <- TOM[rowSums(TOM >= 0.1) > 1, colSums(TOM >= 0.1) > 1]
probes = colnames(TOM_red)
node_attr <- cbind(geneInfo_signed[ geneInfo_signed$genes %in% colnames(TOM_red),c("moduleColor","DE")],young_old.CV_connectivity_expression[ young_old.CV_connectivity_expression$Gene %in% colnames(TOM_red),8:9])
cyt = exportNetworkToCytoscape(TOM_red,edgeFile = "Cytoscape_edges_all.txt", nodeFile = "Cytoscape_nodes_all.txt",weighted = TRUE, nodeNames = probes, nodeAttr = node_attr, threshold=0.1)

##Show network for all genes in modules of interest
#get all genes from modules of interest plus limited genes from rest
modules = c("paleturquoise", "red", "greenyellow", "grey60", "darkred", "cyan", "magenta", "violet", "yellowgreen")
module_genes =  row.names(geneInfo_signed[ geneInfo_signed$moduleColor %in% modules,])
limited_genes = row.names(TOM[rowSums(TOM >= 0.15) > 1, colSums(TOM >= 0.15) > 1])
probes = unique(c(module_genes,limited_genes))


probes <- row.names(geneInfo_signed[ (geneInfo_signed$moduleColor %in% modules) | (geneInfo_signed$IntraConnectivity > 0.3),])
TOM_red = TOM[probes,probes]
moduleColor_red <- as.character(geneInfo_signed[ geneInfo_signed$genes %in% probes,"moduleColor"])
cyt = exportNetworkToCytoscape(TOM_red,edgeFile = "Cytoscape_edges_hubs.txt", nodeFile = "Cytoscape_nodes_hubs.txt",weighted = TRUE, nodeNames = probes, nodeAttr = moduleColor_red, threshold=0.05)


