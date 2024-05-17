library(fishpond)
library(tximport)
library(GenomicFeatures)
library(readr)
library(ggplot2)


library(devtools)
library(tximeta)
library(tidytable)
library(ggplot2)
library(DESeq2)
library(dplyr)


setwd("~/kvdlab/data/RNA-seq")
makeLinkedTxome(
  indexDir = "Human44.HPV18.transcripts.fas_index",
  source = "de-novo",
  organism = "Homo sapiens",
  release = "110",
  genome = "GRCh38",
  fasta = "Human44.HPV18.genome.transcripts.fas",
  gtf = "Human44.HPV18.transcripts10.gtf",
  write = FALSE)

coldata = read.csv("samples.csv")
coldata$condition = factor(coldata$condition,
                            levels = c("inhibitor","vehicle"))
coldata$siRNA = factor(coldata$siRNA, levels = c("mock","KD"))

y = tximeta(coldata, useHub=FALSE) # reads in counts and inf reps
#y = tximeta(coldata, useHub=FALSE, skipMeta=TRUE) # reads in counts and inf reps
assayNames(y)
(rownames(y))


library(fishpond)

y = scaleInfReps(y)
y = labelKeep(y)
y = y[mcols(y)$keep,]
set.seed(1)
y = swish(y, x="condition", pair="replicate")

data = (mcols(y))
data = data[order(data$log2FC), ]
data@rownames
data.HPV18 = as.data.frame(data[data$tx_name %like% "HPV18",])
data.HPV18$tx_name


library(tidyverse)

col_number = max(str_count(data.HPV18$tx_name, "_") + 1)
data.HPV18$tx_name_dup = data.HPV18$tx_name
data.HPV18= data.HPV18 %>% separate(tx_name_dup,
                                    into = paste0("tx_name", seq_len(col_number)),
                                    sep = "_")
data.HPV18$transcript = paste(data.HPV18$tx_name3, data.HPV18$tx_name4,sep = "_")
data.HPV18.2 = subset(data.HPV18, transcript %in% agg_tbl$transcript)

intron = read.csv("intron_count.csv") #update

#for each experiment count the total number of mapped reads as a normalizer
normalizer = dat %>% group_by(replicate, treatment) %>%
  summarise(normalizer=sum(count))


#add intron data and a column that allows for propoer graphing
df_list = list(data.HPV18.2,intron) #combine the intron count data

data.HPV18.2 = df_list %>% reduce(full_join, by='transcript')
data.HPV18.2$graph = paste(data.HPV18.2$bar,data.HPV18.2$tx_name)
data.HPV18.2  = subset(data.HPV18.2, intron.count != 'NA')

colourCount = length(unique(data.HPV18.2$transcript))
getPalette = colorRampPalette(brewer.pal(12, "Dark2"))

ggplot(data.HPV18.2) +
  geom_bar( aes(y=graph, x=log2FC, fill = transcript), stat="identity") +
  geom_text(data = data.HPV18.2, aes(y=graph, x=log2FC,label = round(as.numeric(qvalue), digits=2) ), colour = "white", hjust = -0.2) +
  geom_vline(xintercept=0, color = "black", linewidth = 0.75) +
  theme_classic() +
  theme(legend.position = "none")

ggsave("log2FC.pdf")

#write.csv(data, "all_TRANSCRIPTS.csv")

write.csv(apply(data,2,as.character), "all_TRANSCRIPTS.csv")
significant = data[data[, "qvalue"] < 0.05, ]
significant = significant[order(significant$log2FC), ]
write.csv(apply(significant,2,as.character), "significant_TRANSCRIPTS.csv")

up = significant[significant[, "log2FC"] > 0, ]
up

write.csv(apply(up,2,as.character), "upregulated_TRANSCRIPTS.csv")
down = significant[significant[, "log2FC"] < 0, ]
write.csv(apply(down,2,as.character), "downregulated_TRANSCRIPTS.csv")
down

data = as_tidytable(mcols(y))
data

threshold = as.factor(ifelse(data$pvalue <= 0.05 & abs(data$log2FC) >= log2(1.5) , ifelse(data$log2FC >= log2(1.5) ,'Up','Down'),'Not'))
pdf("plot1.pdf")

ggplot(data,aes(x=log2FC,y=-log10(pvalue),colour=threshold)) +
  xlab("log2(Fold Change)")+ylab("qvalue") +
  geom_point(size = 2,alpha=1) +
  ylim(0,7) + xlim(-5,5) +
  scale_color_manual(values=c("blue","grey", "red"))+
  geom_vline(xintercept = c(-log2(1.5), log2(1.5)), lty = 2,colour="#000000")+
  geom_hline(yintercept = c(-log10(0.05)), lty = 2,colour="#000000") +
  ggtitle("B") +
  guides(fill = guide_legend(reverse = F))+
  theme(plot.title = element_text(hjust = -0.06,size = 28, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 15, face = "bold"),
        legend.position = 'right',
        legend.key.size=unit(0.4,'cm'))+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text=element_text(size=12,face = "bold"),axis.title.x=element_text(size=15),axis.title.y=element_text(size=15))

gse = summarizeToGene(y)
gy = gse
gy = scaleInfReps(gy)
gy = labelKeep(gy)
gy = gy[mcols(gy)$keep,]
set.seed(1)
gy = swish(gy, x="condition", pair="replicate")

data = (mcols(gy))
data
data = data[order(data$log2FC), ]

significant = data[data[, "qvalue"] < 0.05, ]
significant = significant[order(significant$log2FC), ]

up = significant[significant[, "log2FC"] > 0, ]
up
write.csv(apply(up,2,as.character), "upregulated_GENES.csv")
down = significant[significant[, "log2FC"] < 0, ]
down

write.csv(apply(down,2,as.character), "downregulated_GENES.csv")

plotInfReps(gy, "ENSG00000111206.13", x="condition", cov="replicate") #FoxM1
plotInfReps(gy, "early", x="condition", cov="replicate") #AURKB

#ENSG00000073111

colData(gy)

iso = isoformProportions(y)
iso = swish(iso, x="condition", pair="replicate")

data = (mcols(iso))
data
data = data[order(data$log2FC), ]

significant = data[data[, "qvalue"] < 0.05, ]
significant = significant[order(significant$log2FC), ]

up = significant[significant[, "log2FC"] > 0, ]
up
write.csv(apply(up,2,as.character), "upregulated_DTU.csv")

down = significant[significant[, "log2FC"] < 0, ]
down
write.csv(apply(down,2,as.character), "downregulated_DTU.csv")