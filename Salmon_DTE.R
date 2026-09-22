setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(readr)
library(dplyr)
#BiocManager::install("tximeta")
library(tximeta)
#BiocManager::install("fishpond")
library(fishpond)
suppressPackageStartupMessages(library(SummarizedExperiment))
library(ggplot2)

#define objects for building column data table
sampleNames=paste0("T",rep(c("2","3"),each=4),rep(c("16","Y"),each=2),rep(c("_DMSO","_GW"),each=1))
fileNames=paste0("/Users/Christina/Documents/RNA_seq_25Apr2024_Salmon/Salmon_quant/","T",rep(c("2","3"),each=4),rep(c("16","Y"),each=2),rep(c("_DMSO","_GW"),each=1),"_salmon/","quant.sf")
hpvstatus=paste0(rep(c("pos","neg"),each=2))
donor=paste0(rep(c("2","3"),each=4))
cells=paste0(rep(c("T216","T2Y","T316","T3Y"),each=2))
treatment=paste0(rep(c("DMSO","GW"),each=1))

#build column data table
coldata <- data.frame(files=fileNames, names=sampleNames, Treatment=factor(treatment), HPVstatus=factor(hpvstatus), Donor=factor(donor), Cells=factor(cells), stringsAsFactors=FALSE)


all(file.exists(coldata$files))



y <- tximeta(coldata, skipMeta = T) # reads in counts and inf reps
y
df <- assay(y) #counts by default, raw counts from Salmon
data <- df[78:100,] #HPV16 genes only
write.csv(data, "HPVcounts.csv")

y
df2 <- assays(y)$infRep1 #now look at infrep1
df2[78:100,]

plotInfReps(y, idx="HPV16_W12_REF_Species_a", x="Treatment", legend = TRUE) #unscaled, raw (see y-axis label)

y <- scaleInfReps(y) # scales counts
plotInfReps(y, idx="HPV16_W12_REF_Species_a", "Treatment", legend = TRUE) 
# "highly variable infRep" probably related to low abundance and close to detection limit
plotInfReps(y, idx="HPV16_W12_REF_Species_d", "Cells", legend = TRUE) #similar range of error around the mean, just more of a scale dominator at low range
plotInfReps(y, idx="HPV16_W12_REF_Species_d", "Cells", legend = TRUE) 
plotInfReps(y, idx="HPV16_W12_REF_Species_e", "Cells", legend = TRUE)
plotInfReps(y, idx="HPV16_W12_REF_Species_i", "Cells", legend = TRUE) 
plotInfReps(y, idx="HPV16_W12_REF_Species_k", "Cells", legend = TRUE) 
plotInfReps(y, idx="HPV16_W12_REF_Species_s", "Cells", legend = TRUE)  
plotInfReps(y, idx="HPV16_W12_REF_Species_u", "Cells", legend = TRUE) 
plotInfReps(y, idx="HPV16_W12_REF_Species_w", "Cells", legend = TRUE) 

#late promoted, there are few counts
plotInfReps(y, idx="HPV16_W12_REF_Species_r", "Cells", legend = TRUE) 
plotInfReps(y, idx="HPV16_W12_REF_Species_s", "Cells", legend = TRUE) 
plotInfReps(y, idx="HPV16_W12_REF_Species_t", "Cells", legend = TRUE) 
plotInfReps(y, idx="HPV16_W12_REF_Species_p", "Cells", legend = TRUE) 
plotInfReps(y, idx="HPV16_W12_REF_Species_q", "Cells", legend = TRUE) 


#intermediate promoter...?
plotInfReps(y, idx="HPV16_W12_REF_Species_u", "Cells", legend = TRUE) 
plotInfReps(y, idx="HPV16_W12_REF_Species_v", "Cells", legend = TRUE) 
plotInfReps(y, idx="HPV16_W12_REF_Species_w", "Cells", legend = TRUE) 

# The following lines of code will perform a basic transcript-level `swish` two group analysis of bulk RNA-seq.
#`swish` has three steps: scaling the inferential replicates, labeling the rows with sufficient counts for running differential
# expression, and then calculating the statistics
#y <- scaleInfReps(y) #I did this above
y <- labelKeep(y) # labels features to keep (ones with sufficient counts)
set.seed(1) #`swish` uses pseudo-random number generation in breaking ties and in calculating permutations, so to obtain identical results, one needs to set a random seed before running `swish()`
# simplest Swish ("SAMseq With Inferential Samples Helps") case (a nonparametric differential analysis)
y <- swish(y, "Treatment") # x must be factor!!!!
# head(assay(y))


#ppara transcripts
plotInfReps(y, idx="ENST00000440343.5|ENSG00000186951.17|OTTHUMG00000150443.7|OTTHUMT00000318747.4|PPARA-205|PPARA|368|protein_coding|", "Cells", legend = TRUE)
#this is expressed highly
plotInfReps(y, idx="ENST00000407236.6|ENSG00000186951.17|OTTHUMG00000150443.7|OTTHUMT00000318129.3|PPARA-202|PPARA|10118|protein_coding|", "Cells", legend = TRUE)
plotInfReps(y, idx="ENST00000407236.6|ENSG00000186951.17|OTTHUMG00000150443.7|OTTHUMT00000318129.3|PPARA-202|PPARA|10118|protein_coding|", "Treatment", legend = TRUE)

plotInfReps(y, idx="ENST00000415785.5|ENSG00000186951.17|OTTHUMG00000150443.7|OTTHUMT00000318764.4|PPARA-203|PPARA|531|protein_coding|", "Cells", legend = TRUE)
plotInfReps(y, idx="ENST00000420804.5|ENSG00000186951.17|OTTHUMG00000150443.7|OTTHUMT00000318748.3|PPARA-204|PPARA|572|protein_coding|", "Cells", legend = TRUE)
plotInfReps(y, idx="ENST00000460086.1|ENSG00000186951.17|OTTHUMG00000150443.7|OTTHUMT00000318749.2|PPARA-206|PPARA|433|protein_coding_CDS_not_defined|", "Cells", legend = TRUE)
plotInfReps(y, idx="ENST00000496865.1|ENSG00000186951.17|OTTHUMG00000150443.7|OTTHUMT00000318765.1|PPARA-210|PPARA|449|protein_coding_CDS_not_defined|", "Cells", legend = TRUE)

plotInfReps(y, idx="ENST00000624793.1|ENSG00000186951.17|OTTHUMG00000150443.7|OTTHUMT00000479013.1|PPARA-211|PPARA|1556|retained_intron|", "Cells", legend = TRUE)
#none
plotInfReps(y, idx="ENST00000484619.1|ENSG00000186951.17|OTTHUMG00000150443.7|OTTHUMT00000318752.1|PPARA-208|PPARA|278|protein_coding_CDS_not_defined|", "Cells", legend = TRUE)
#none
plotInfReps(y, idx="ENST00000481567.5|ENSG00000186951.17|OTTHUMG00000150443.7|OTTHUMT00000318751.1|PPARA-207|PPARA|567|protein_coding_CDS_not_defined|", "Cells", legend = TRUE)
#none
plotInfReps(y, idx="ENST00000402126.1|ENSG00000186951.17|OTTHUMG00000150443.7|OTTHUMT00000318753.2|PPARA-201|PPARA|1558|protein_coding|", "Cells", legend = TRUE)
#low
plotInfReps(y, idx="ENST00000493286.1|ENSG00000186951.17|OTTHUMG00000150443.7|OTTHUMT00000318750.1|PPARA-209|PPARA|1564|retained_intron|", "Cells", legend = TRUE)
plotInfReps(y, idx="", "Cells", legend = TRUE)
plotInfReps(y, idx="", "Cells", legend = TRUE)


#The results can be found in `mcols(y)`. 
#calculate the number of genes passing a 5% FDR threshold:
table(mcols(y)$qvalue < .05)

# #rank features by their effect size, to break ties in the `qvalue`
# most.sig <- with(mcols(y),
#                  order(qvalue, -abs(log2FC)))
# View(mcols(y)[head(most.sig),c("log2FC","qvalue")])
# plotMASwish(y, alpha=0.5)
# 
# #construct two vectors that give the significant genes with the lowest (most negative) and highest (most positive) log2 fold changes.
# with(mcols(y),
#      table(sig=qvalue < .05, sign.lfc=sign(log2FC))
# )
# sig <- mcols(y)$qvalue < .05
# head(sig)
# lo <- order(mcols(y)$log2FC * sig)
# hi <- order(-mcols(y)$log2FC * sig)

#Here we print a small table with just the calculated statistics for the large positive log fold change transcripts (up-regulation):
# top_up <- mcols(y)[head(hi),]
# names(top_up)
# cols <- c("counts","log10mean","log2FC","pvalue","qvalue")
# print(as.data.frame(top_up)[,cols], digits=3)
# 
# #Likewise for the largest negative log fold change transcripts (down-regulation): 
# top_down <- mcols(y)[head(lo),]
# print(as.data.frame(top_down)[,cols], digits=3)



# #plot scaled counts grouping by covariates
# plotInfReps(y, idx=hi[1] x="Treatment", legend = TRUE)
# 
# 
# plotInfReps(y, idx=hi[1], x="Donor", cov="Treatment", legend = TRUE)




#make a table containing the counts of each gene (rows) for each sample (columns)
counttable <- assay(y)
head(counttable)


#list of E1^E4 containing transcripts
transcripts=c("HPV16_W12_REF_Species_b", "HPV16_W12_REF_Species_c", "HPV16_W12_REF_Species_d", "HPV16_W12_REF_Species_e", "HPV16_W12_REF_Species_q", "HPV16_W12_REF_Species_r", "HPV16_W12_REF_Species_s")

#count table for just E1^E4 containing transcripts
E4counttable <- assay(y)[rownames(y) %in% transcripts, ]
View(E4counttable)
allE4counts <- colSums (E4counttable, dims = 1)
View(allE4counts)

#plot scaled counts grouping by covariates
ggplot(allE4counts)#fix this....

#preserve the data frame
### i dont know how to add the values together within the data frame...
E4countmatrix <- y[rownames(y) %in% transcripts, ]
View(E4countmatrix)

#c and d are the dominant E1^E4 containing species
plotInfReps(y, idx="HPV16_W12_REF_Species_c", x="Donor", cov="Treatment", legend = TRUE, xaxis=FALSE, xlab="Treatment")
plotInfReps(y, idx="HPV16_W12_REF_Species_d", x="Treatment", cov="Donor", legend = TRUE)





#Redo with just HPV+ samples


#define objects for building column data table
sampleNames_HPV=paste0("T",rep(c("2","3"),each=2),"16",rep(c("_DMSO","_GW"),each=1))
fileNames_HPV=paste0("/Users/Christina/Documents/RNA_seq_25Apr2024_Salmon/Salmon_quant/","T",rep(c("2","3"),each=2),"16",rep(c("_DMSO","_GW"),each=1),"_salmon/","quant.sf")
donor_HPV=paste0(rep(c("2","3"),each=2))
cells_HPV=paste0(rep(c("T216","T316"),each=1))
treatment_HPV=paste0(rep(c("DMSO","GW"),each=1))

#build column data table
coldata_HPV <- data.frame(files=fileNames_HPV, names=sampleNames_HPV, Treatment=factor(treatment_HPV), Donor=factor(donor_HPV), Cells=factor(cells_HPV), stringsAsFactors=FALSE)

all(file.exists(coldata_HPV$files))

y_H <- tximeta(coldata_HPV, skipMeta=T) # reads in counts and inf reps
y_H
head(assay(y_H))


# The following lines of code will perform a basic transcript-level `swish` two group analysis of bulk RNA-seq.
#`swish` has three steps: scaling the inferential replicates, labeling the rows with sufficient counts for running differential
# expression, and then calculating the statistics
y_H <- scaleInfReps(y_H) # scales counts
y_H <- labelKeep(y_H) # labels features to keep (ones with sufficient counts)
set.seed(1) #`swish` uses pseudo-random number generation in breaking ties and in calculating permutations, so to obtain identical results, one needs to set a random seed before running `swish()`
# simplest Swish ("SAMseq With Inferential Samples Helps") case (a nonparametric differential analysis)
y_H <- swish(y_H, "Treatment") # x must be factor!!!!
head(assay(y_H))



#c and d are the dominant E1^E4 containing species
plotInfReps(y_H, idx="HPV16_W12_REF_Species_c", x="Treatment", cov="Donor", legend = TRUE)
plotInfReps(y_H, idx="HPV16_W12_REF_Species_d", x="Treatment", cov="Donor", legend = TRUE)





# thought I needed a local index, but actually I don't need the index because that is just for mapping reads to the genome
#and I only want the transcript counts not the gene counts
# makeLinkedTxome(
#   indexDir="/xdisk/vandoorslaer/indices/SALMON/Human44.HPV18.transcripts.fas_index",
#   source="myGENCODE",
#   organism="Homo sapiens",
#   release="29",
#   genome="GRCh38",
#   fasta="/xdisk/vandoorslaer/genome_transcripts_GTF/Human44.HPV18.genome.transcripts.fas",
#   gtf="",
#   write=FALSE
# )
