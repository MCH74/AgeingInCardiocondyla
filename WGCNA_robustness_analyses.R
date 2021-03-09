##increased data set
design <- read.table("../cobs.design")
shams <- c("SRR2177530","SRR2177529","SRR2177528","SRR2177527","SRR2177542","SRR2177543","SRR2177544")
design_sham <- design
#create design table
for(i in 1:length(shams)){
sample15 = shams[i]
levels(design_sham$Sample) <- c(levels(design_sham$Sample), sample15)
levels(design_sham$mating) <- c(levels(design_sham$mating), "sham")
design_sham <- rbind(design_sham, c(sample15,"sham","old",1))
}
rownames(design_sham) <- design_sham$Sample


#get normalised counts
counts_sham <- counts[,rownames(design_sham)]

dds <- DESeqDataSetFromMatrix(countData = counts_sham, colData = design_sham, design= ~age)
dds <- DESeq(dds)                                                   
datExpr <- t(counts(dds, normalized=TRUE))   

geneInfo_14 <- read.table("../geneInfo_cobs_signed.tab",header=T)  
datExpr_red_full <- datExpr[ ,colnames(datExpr) %in% geneInfo_14$genes ]


## now make nework for each sham

for(i in 15:21){
sample15 = rownames(design_sham)[i]
trait_file <- paste0("moduleTraitCor_",sample15)
geneinfo_file <- paste0("geneInfo_", sample15)

datExpr_red <- datExpr_red_full[ c(1:14,i) ,]


##pick 14 as power
softPower = 14
adjacency = adjacency(datExpr_red, power = softPower, type = "signed", corFnc = "bicor", corOptions = list(maxPOutliers = 0.05));
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM
geneTree = hclust(as.dist(dissTOM), method = "average");

minModuleSize = 30;
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
deepSplit = 2, pamRespectsDendro = FALSE,
minClusterSize = minModuleSize);
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)

# Calculate eigengenes
MEList = moduleEigengenes(datExpr_red, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");

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

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;


###relate modules to traits

# Define numbers of genes and samples
nGenes = ncol(datExpr_red);
nSamples = 14 #nrow(datExpr_red);  we're excluding the virgin sample
# Recalculate MEs with color labels
moduleTraitCor = cor(MEs[1:14,], as.numeric(design_15[1:14,4]), use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

moduletrait <- cbind(moduleTraitCor ,moduleTraitPvalue)
colnames(moduletrait) <- c("Cor","p")
write.table(moduletrait, trait_file, quote=F, sep="\t")

###summarise gene info in a table

geneTraitSignificance = as.data.frame(cor(datExpr_red, as.numeric(design_15[,4]), use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", "Age", sep="");
names(GSPvalue) = paste("p.GS.", "Age", sep="");


geneInfo_signed = data.frame(genes = colnames(datExpr_red),
                    moduleColor = moduleColors,
                    geneTraitSignificance,
                    GSPvalue)


modOrder = order(-abs(cor(MEs, as.numeric(design_15[,4]), use = "p")));
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

cobs.connectivity <- softConnectivity( datExpr_red, power = softPower, type = "signed", corFnc = "bicor", corOptions = list(maxPOutliers = 0.05)  )
geneInfo_signed <- cbind(geneInfo_signed, connectivity = cobs.connectivity)

kIM <- intramodularConnectivity(adjacency, moduleColors, scaleByMax = TRUE) 
geneInfo_signed <- cbind(geneInfo_signed, IntraConnectivity = kIM[,2])

write.table(geneInfo_signed, geneinfo_file, sep="\t",quote=F, row.names=F)


#######
#now run preservation stats

setLabels = c("small", "large");
multiExpr = list(small= list(data = datExpr_red[1:14,]), large = list(data = datExpr_red)); # make sure genes are in same order
multiColor = list(small = geneInfo_14$moduleColor, large = geneInfo_signed$moduleColor);

system.time( {
mp_14v15 = modulePreservation(multiExpr, multiColor,
referenceNetworks = 1, #small
nPermutations = 200,
randomSeed = 1,
quickCor = 0,
verbose = 3)
} );


ref = 1
test = 2
statsObs = cbind(mp_14v15$quality$observed[[ref]][[test]][, -1], mp_14v15$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp_14v15$quality$Z[[ref]][[test]][, -1], mp_14v15$preservation$Z[[ref]][[test]][, -1]);
stats.tab <- cbind(statsObs[, c("medianRank.pres", "medianRank.qual")], signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2))

stats_file <- paste0(sample15, "_pres.stats")
write.table(stats.tab, stats_file, sep="\t",quote=F, row.names=T)



}

