#This script is supposed to extract all the CAZymes data from the dbCAN website:
#The cazymes were obtained from the funannotate pipleline using the following command:
#for i in *; do seqtk subseq ../combined.proteins/$i/predict_results/*.fa $i/annotate_misc/*.txt | sed "s/>/>$i\_/" >> /Users/user/Downloads/all_CAZymes.fa;done
library("tidyr")
library("ggplot2")
library("ape")
library("gridExtra")
library("ggfortify")
library("phytools")
library("ggtree")
library("treeio")
library("ggnewscale")
library("ggrepel")
library("gridExtra")
library("ggimage")



setwd("/Users/ben/Desktop/lennart_from_ben")
#genome summary stats----
#Now we want to load the data
stats <- read.csv("genome.stats.csv")

summary(lm(stats$Million.Reads ~ stats$CAZymes))
summary(lm(stats$Million.Reads ~ stats$Buscos))
summary(lm(stats$Million.Reads ~ stats$total.proteins))
summary(lm(stats$CAZymes ~ stats$total.proteins))

p1 <- ggplot() + theme_bw() + xlim(0,35) + ylim(0,500) +
  geom_smooth(aes(x=stats$Million.Reads,y=stats$CAZymes),method="lm",level=0,lty=1)+
  geom_point(aes(x=stats$Million.Reads,y=stats$CAZymes))+
  annotate("text",size=6,x=20,y=100,label=expression(paste(italic(R)^"2","=0.001")))+
  geom_segment(aes(x=11.0,y=175,xend=13.45,yend=200),arrow=arrow(length=unit(0.1,"inches")))+
  labs(x="Reads per sample (millions)",y="CAZyme count",title="A) CAZyme vs. Read Depth")
p2 <- ggplot() + theme_bw() + xlim(0,35) + ylim(0,1200) +
  geom_smooth(aes(x=stats$Million.Reads,y=stats$Buscos),method="lm",level=0,lty=1)+
  geom_point(aes(x=stats$Million.Reads,y=stats$Buscos)) +
  annotate("text",size=6,x=20,y=250,label=expression(paste(italic(R)^"2","=0.033")))+
  geom_segment(aes(x=11.8,y=1200,xend=13.5,yend=1150),arrow=arrow(length=unit(0.1,"inches")))+
  labs(x="Reads per sample (millions)",y="BUSCO count",title="B) BUSCO vs. Read Depth")
p3 <- ggplot() + theme_bw() + xlim(0,35) + ylim(0,20000) +
  geom_smooth(aes(x=stats$Million.Reads,y=stats$total.proteins),method="lm",level=0,lty=1)+
  geom_point(aes(x=stats$Million.Reads,y=stats$total.proteins)) +
  annotate("text",size=6,x=20,y=5000,label=expression(paste(italic(R)^"2","=0.008")))+
  labs(x="Reads per sample (millions)",y="Total # predicted proteins",title="C) Predicted Proteins vs. Read Depth")
p4 <- ggplot() + theme_bw() + xlim(0,20000) + ylim(0,500) +
  geom_point(aes(x=stats$total.proteins,y=stats$CAZymes)) + 
  geom_smooth(aes(x=stats$total.proteins,y=stats$CAZymes),method="lm",level=0,lty=1) +
  annotate("text",size=6,x=10000,y=100,label=expression(paste(italic(R)^"2","=0.225")))+
  labs(x="Total # predicted proteins",y="CAZyme count",title="D) CAZyme vs. Predicted Proteins")
grid.arrange(p1,p2,p3,p4,nrow=2)
pdf("/Users/ben/Desktop/lennart_from_ben/lennart.supplement.genome.stats.pdf")
grid.arrange(p1,p2,p3,p4,nrow=2)
dev.off()

#cazyme data prep----
##now lets parse the dbCAN results from the webserver
dbCAN.download <- read.delim("dbCAN.results.Nov2.txt",header=T,sep="\t")
tree <- read.tree("rooted.concord.cf.tree")
data <- read.csv("Table_morphologicaltraitsV2_onlyincludedtaxa.csv")
tree$tip.label <- gsub("\\.1","",tree$tip.label)
data$D. <- gsub("\\.1","",data$D.)
#the # of tools gets a funny name, so we fix it.
colnames(dbCAN.download)[6] <- "num_of_tools"

tree$tip.label
tree$tip.label %in% data$D.

#replace tree tips with data values
tree$tip.label <- data$work_name[match(tree$tip.label,data$D.)]
#lets remove the bootstrap values, since they are all 100
tree$node.label <- gsub("^100/","",tree$node.label)
tree$node.label <- gsub("\\.[0-9]","",tree$node.label,perl=T)
tree$node.label

#lets get a summary of how many identified iwth each tools
sum(dbCAN.download$num_of_tools == 0)
sum(dbCAN.download$num_of_tools == 1)
sum(dbCAN.download$num_of_tools == 2)
sum(dbCAN.download$num_of_tools == 3)

#so less than half are identified by all three, and only 3/4 are identified by two, but lets summarize by that criteria
dbCAN.twotools <- dbCAN.download[dbCAN.download$num_of_tools > 1,]

#some (<20) have no name from HMMER, so we will transfer from DIAMOND
dbCAN.twotools$HMMER[dbCAN.twotools$HMMER == "N"]<- dbCAN.twotools$DIAMOND[dbCAN.twotools$HMMER=="N"]

#now we want to split enzymes that have more than one domian from HMMER
dbCAN.twotools <- separate_rows(dbCAN.twotools,HMMER,sep = "\\+")

#now we want to simplify the names from HMMER and split the proteins into just the organims
dbCAN.twotools$HMMER <- gsub("\\s*\\([^\\)]+\\)","",dbCAN.twotools$HMMER)
#we need to strip everything after the underscore, but need to switch NCBI names to keep underscore, then switch back
dbCAN.twotools$Gene.ID <- sub("GCA_","GCAx",dbCAN.twotools$Gene.ID)
dbCAN.twotools$Gene.ID <- gsub("\\_.*","",dbCAN.twotools$Gene.ID)
dbCAN.twotools$Gene.ID <- sub("GCAx","GCA_",dbCAN.twotools$Gene.ID)
dbCAN.twotools$Gene.ID <- gsub("\\.1","",dbCAN.twotools$Gene.ID)

#assign genus names
dbCAN.twotools$genus <- "Other"
dbCAN.twotools$genus[dbCAN.twotools$Gene.ID %in% c("D8", "D9", "D10","D11","D12","D13",
                                                   "D14","D15","D16","D17","D19",
                                                   "D20","D21","D22","D25","D26",
                                                   "D27","D28","D36","D41","D42",
                                                   "D43","D44","D45","D46","D52",
                                                   "D53","D54","D55","D58","D59")] <- "Termitomyces"
dbCAN.twotools$genus[dbCAN.twotools$Gene.ID %in% c("GCA_001263195.1","GCA_001972325.1","GCA_003313055.1",
                                                   "GCA_003313075.1","GCA_003313675.1","GCA_003313785.1",
                                                   "GCA_003316525.1")] <- "Termitomyces"

dbCAN.twotools$genus[dbCAN.twotools$Gene.ID %in% c("D32","D50")] <- "Arthromyces"
dbCAN.twotools$genus[dbCAN.twotools$Gene.ID %in% c("D30","D24")] <- "Blastosporella"
dbCAN.twotools$genus[dbCAN.twotools$Gene.ID %in% c("D31","D33","D64","D39",
                                                   "D62","D18","D29","D63")] <- "Termitomycetoid"


dbCAN.twotools$Gene.ID <- data$work_name[match(dbCAN.twotools$Gene.ID,data$D.)]

#now make short names
dbCAN.twotools$HMMER_longname <- dbCAN.twotools$HMMER
#alternate name cleanings
dbCAN.twotools$HMMER <- gsub("AA.*","AA",dbCAN.twotools$HMMER)
dbCAN.twotools$HMMER <- gsub("GH.*","GH",dbCAN.twotools$HMMER)
dbCAN.twotools$HMMER <- gsub("PL.*","PL",dbCAN.twotools$HMMER)
dbCAN.twotools$HMMER <- gsub("CBM.*","CBM",dbCAN.twotools$HMMER)
dbCAN.twotools$HMMER <- gsub("CE.*","CE",dbCAN.twotools$HMMER)
dbCAN.twotools$HMMER <- gsub("GT.*","GT",dbCAN.twotools$HMMER)

#now we want to keep the subgroups
dbCAN.twotools$HMMER_longname <- gsub("_Chitin_synth*","",dbCAN.twotools$HMMER_longname)
dbCAN.twotools$HMMER_longname <- gsub("_Glyco_trans.*","",dbCAN.twotools$HMMER_longname)
dbCAN.twotools$HMMER_longname <- gsub("_Glycos_trans.*","",dbCAN.twotools$HMMER_longname)
dbCAN.twotools$HMMER_longname <- gsub("_.*","",dbCAN.twotools$HMMER_longname)

dbCAN.twotools$HMMER <- as.factor(dbCAN.twotools$HMMER)
dbCAN.twotools$HMMER_longname <- as.factor(dbCAN.twotools$HMMER_longname)

dbCAN.twotools.wide <- as.data.frame(table(genus=dbCAN.twotools$genus,organism=dbCAN.twotools$Gene.ID,dbCAN.twotools$HMMER))
dbCAN.twotools.wide <- pivot_wider(dbCAN.twotools.wide,names_from = "Var3",values_from = "Freq")
#for some reason it is making lots of rows with zero values, because the of the interaction between the genus column, I don't know where the error is so,
#lets just delete them
dbCAN.twotools.wide <- dbCAN.twotools.wide[rowSums(dbCAN.twotools.wide[,3:8]) > 0,]

dbCAN.twotools.wide_longname <- as.data.frame(table(genus=dbCAN.twotools$genus,organism=dbCAN.twotools$Gene.ID,dbCAN.twotools$HMMER_longname))
dbCAN.twotools.wide_longname <- pivot_wider(dbCAN.twotools.wide_longname,names_from = "Var3",values_from = "Freq")
dbCAN.twotools.wide_longname <- dbCAN.twotools.wide_longname[rowSums(dbCAN.twotools.wide_longname[,c(3:100)])!=0,]
dbCAN.twotools.wide_longname$total <- rowSums(dbCAN.twotools.wide_longname[,c(3:119)])

dbCAN.twotools.trans_longname <- as.data.frame(t(dbCAN.twotools.wide_longname))
dbCAN.twotools.trans_longname <- cbind(rownames(dbCAN.twotools.trans_longname),dbCAN.twotools.trans_longname)
dbCAN.twotools.trans_longname <- dbCAN.twotools.trans_longname[2:nrow(dbCAN.twotools.trans_longname),]
dbCAN.twotools.trans_longname <- cbind(rep(as.character("(null)"),nrow(dbCAN.twotools.trans_longname)),dbCAN.twotools.trans_longname)
dbCAN.twotools.trans_longname[[1]][[1]] <- "desc"
dbCAN.twotools.trans_longname[[2]][[1]] <- "ID"
write.table(dbCAN.twotools.trans_longname[,2:ncol(dbCAN.twotools.trans_longname)],"/Users/ben/Desktop/lennart_from_ben/dbCAN.totals.csv",col.names=F,sep=",",row.names = F)
#Now scale data by column
dbCAN.twotools.scaled <- as.data.frame(cbind(as.character(dbCAN.twotools.wide$organism), scale(dbCAN.twotools.wide[,3:8])))
colnames(dbCAN.twotools.scaled)[1] <- "organism"
#row.names(dbCAN.twotools.scaled) <- dbCAN.twotools.wide$organism
#colMeans(as.numeric(dbCAN.twotools.scaled[,2:7]))
#apply(as.numeric(dbCAN.twotools.scaled[,2:7]), 2, sd)

#first we need it in long format
dbCAN.twotools.scaled.longer <- pivot_longer(dbCAN.twotools.scaled,cols=c(2:7))
dbCAN.twotools.scaled.longer$organism <- factor(dbCAN.twotools.scaled.longer$organism,levels=tree$tip.label)
dbCAN.twotools.scaled.longer$value <- as.numeric(dbCAN.twotools.scaled.longer$value)
dbCAN.twotools.scaled.longer <- dbCAN.twotools.scaled.longer[dbCAN.twotools.scaled.longer$name != "N",]
dbCAN.twotools.wide$organism <- as.character(dbCAN.twotools.wide$organism)
dbCAN.twotools.raw.longer    <- pivot_longer(dbCAN.twotools.wide,cols=c(2:7))

#pdf("split_groups..pdf")
#split_groups <- ggplot(dbCAN.twotools.scaled.longer,
#       aes(y=organism,
#           x=name,
#           fill=value)) + theme_bw()+
#         geom_tile() + scale_fill_distiller(type="div",aesthetics="fill",palette=7) #+
  #geom_label(aes(label=dbCAN.twotools.raw.longer$value),
  #           label.padding=unit(0.05,"lines"),
  #           label.r = unit(0.00, "lines"),
  #           label.size = 0,
  #           size=3,
  #           fill="white")
split_groups
dev.off()
#pdf("term_phylo.Feb12.pdf")
#plot.phylo(tree,show.node.label=T,edge.width=1.5,cex=0.5,align.tip.label=T)
#dev.off()

pPCA <- phyl.pca(tree,dbCAN.twotools.wide_longname[,c(3:119)])
#need to calculate variance for PCs
diag <- diag(pPCA$Eval)
diag <- diag/sum(diag)
plot_pPCA <- ggplot() + geom_point(aes(x=pPCA$S[,1],y=pPCA$S[,2],shape=dbCAN.twotools.wide_longname$genus,color=dbCAN.twotools.wide_longname$genus),size=2.5) + 
  theme_bw() + labs(x = paste0("Phylogenetic PC1 = ",round(diag[1]*100,1),"%"),
                    y = paste0("Phylogenetic PC2 = ",round(diag[2]*100,1),"%"),
                    title="Phylogenetic Corrected PCA of CAZyme Content") +
  scale_color_discrete(name="Genus")+
  scale_shape_discrete(name="Genus")+
  theme(axis.title = element_text(size=15),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15))
pdf("phylogenetic_pca.pdf",width=7,height=5)
plot_pPCA
dev.off()


pca_res <- prcomp(dbCAN.twotools.wide_longname[,c(3:119)])
pca_CAZy <- autoplot(pca_res,data = dbCAN.twotools.wide_longname,shape="genus",colour = "genus",size=2.5,label.label="organism",label=F,label.size=3,label.hjust=0,label.vjust=0) + theme_bw() + 
    scale_shape_manual(values=c(15,16,1,17,18,19))+
    scale_color_manual(values=c("steelblue2","steelblue2","grey30","tomato1","steelblue2"))+
    scale_alpha_manual(values=0.3)+
    theme(axis.title = element_text(size=8),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8))
#pdf("pca.Feb12.pdf",width=7,height=5)
pca_CAZy
#dev.off()

morphology <- data[,c(1,9,10,11,12)]
morphology[morphology == "0"] <- "101"
morphology[morphology == "?"] <- "100"
morphology[morphology == "1"] <- "99"
morphology_long <- pivot_longer(morphology,cols=c(2:5))
morphology_long$work_name <- factor(morphology_long$work_name,levels=tree$tip.label)
morphology_long$value <- as.numeric(morphology_long$value)
rownames(morphology) <- morphology$work_name

tree <- read.tree("rooted.concord.cf.tree")
tree$node.label <- gsub("^100/","",tree$node.label)
data <- read.csv("Table_morphologicaltraitsV2_onlyincludedtaxa.csv")
data$D. <- gsub("\\.1","",data$D.)

tree$tip.label <- data$work_name[match(tree$tip.label,data$D.)]
rownames(dbCAN.twotools.scaled) <- dbCAN.twotools.wide$organism
dbCAN.twotools.scaled[,2:7] <- sapply(dbCAN.twotools.scaled[,2:7],as.numeric)

#lets make some pie charts
split_support <- data.frame(
  node = seq(1+46,tree$Nnode+46),
  start = tree$node.label)
split_support <- separate(split_support,start,into=c("first","second"),sep="/")
split_support$first <- as.numeric(split_support$first)
split_support$second <- as.numeric(split_support$second)
split_support$not_first <- 100-split_support$first
split_support$not_second <- 100-split_support$second
head(split_support)
first_pie <- nodepie(split_support,cols=c(2,4),color=c("black","white"))
second_pie <- nodepie(split_support,cols=c(3,5),color=c("white","black"))

#old offsets
#geom_tiplab 0.133
#morpho 0.13
#cazy 0.20

phylo <- ggtree(tree) %<+% data +
  geom_hilight(node=49,fill="steelblue2")+ 
  geom_hilight(node=54,fill="tomato1",alpha=1) +
  #geom_rect(xmin=0.4265,xmax=0.581,ymin=0,ymax=47,fill="white",color="white") +
  geom_tree()+
  #old tiplabel
  geom_tiplab(align=T,size=2,offset=0.1545,linesize=0.15,hjust=0,geom="label",label.size=0,label.padding=unit(0,"inches")) + 
  #new tiplabel from https://groups.google.com/g/bioc-ggtree/c/BA7g-npY1BM
  #geom_tiplab(aes(label=paste0('italic(', Genus, ')~italic(', Species, ')~bold(', Strain, ')')),parse=T)+
  geom_rootedge() + geom_nodelab(size=1.7,nudge_x=0.0,hjust=0) + geom_nodepoint(size=1.2,aes(x=x-0.003))
phylo
#phylo_one <- ggtree::inset(phylo,first_pie,width=0.08,height=0.1,hjust=0.005,vjust=-0.45)
#phylo_two <- ggtree::inset(phylo_one,second_pie,width=0.1,height=0.1,hjust=0.005,vjust=0.565)
#phylo_one
#phylo_two
#morpho <- gheatmap(phylo,morphology,offset=1.4,width=1.75,low="red",high="blue",color="white",colnames=T,legend_title="Trait",colnames_angle = 90) + scale_x_ggtree()
morpho <- gheatmap(phylo,morphology[2:5],color="grey50",offset=-0.0063,width=0.15,colnames=F,legend_title="Trait",colnames_angle = 90) + 
  scale_fill_manual(breaks=c(99,100,101),values=c("grey30","grey80","white"),name="Trait",labels=c("Present","Unknown","Absent")) + scale_x_ggtree()
#morpho
cazy <- morpho + new_scale_fill()
final <- gheatmap(cazy,dbCAN.twotools.scaled[2:7],offset=0.062,width=0.2,colnames_angle = 90) + scale_fill_viridis_c(option="D", name="CAZyme\nRelative\nCount") +
  theme(legend.margin = margin(c(0,0,0,0))) + xlim(0,0.69)
final
pdf("combined_phylo_Mar3.pdf",width=8,height=7)
grid.arrange(
  final,pca_CAZy,
  heights=c(1.9,1),
  widths = c(1.1, 1),
  layout_matrix = rbind(c(1, 1),
                        c(NA, 2)))
dev.off()


#output for the CAFE analysis
#ultra_tree <- chronos(tree,0.5)
#plot.phylo(ultra_tree)
#ape::write.tree(ultra_tree,"ultra_tree.Dec15.tre")

#now we read the ASTRAL tree in
astral_tree <- read.tree("ASTRAL.BUSCO_25plus_Gblocks_cat.tre")
astral_tree <- root(astral_tree,"D51")
data <- read.csv("Table_morphologicaltraitsV2_onlyincludedtaxa.csv")
data$D. <- gsub("\\.1","",data$D.)

astral_tree$tip.label <- data$work_name[match(astral_tree$tip.label,data$D.)]
#Termitomyces MRCA
getMRCA(astral_tree,c("Termitomyces striatus","Termitomyces fragilis"))
#Termitomycetoid MRCA
getMRCA(astral_tree,c("Termitomyces striatus","Myochromella boudieri"))
astral <- ggtree(astral_tree) + xlim(0,25)+
  geom_hilight(node=81,fill="orange",alpha=0.3)+ 
  geom_hilight(node=68,fill="green",alpha=0.3) +
  geom_tiplab(align=F,size=3,offset=0.01) + geom_rootedge() + geom_nodelab(size=4,nudge_x=0.25) + geom_nodepoint(size=1.2,aes(x=x-0.003))
#astral
pdf("supplemental_astral_tree.pdf")
astral
dev.off()


#read CAFE results
cafe_trees <- read.nexus("/Users/ben/Desktop/lennart_from_ben/Base_asr.tre")
sample_tree <- cafe_trees$AA1
sample_tree$tip.label <- gsub("<.*","",sample_tree$tip.label)
worked_trees <- cafe_trees

pdf("/Users/ben/Desktop/lennart_from_ben/CAFE_results.pdf")
par(mfrow=c(2,3),mai=c(0,0,00.5,0))
plot.phylo(sample_tree,cex=0.75)
title("Species Phylogeny")

worked_trees$AA1$node.label[grepl("\\*",cafe_trees$AA1$node.label) == F] <- NA
worked_trees$AA1$node.label[grepl("\\*",cafe_trees$AA1$node.label) == T] <- "X"
worked_trees$AA1$tip.label[grepl("\\*",cafe_trees$AA1$tip.label) == F] <- NA
worked_trees$AA1$tip.label[grepl("\\*",cafe_trees$AA1$tip.label) == T] <- "X"
plot.phylo(worked_trees$AA1,show.node.label = T,cex=0.7,edge.width=0.15)
title("AA1")

worked_trees$AA3$node.label[grepl("\\*",cafe_trees$AA3$node.label) == F] <- NA
worked_trees$AA3$node.label[grepl("\\*",cafe_trees$AA3$node.label) == T] <- "X"
worked_trees$AA3$tip.label[grepl("\\*",cafe_trees$AA3$tip.label) == F] <- NA
worked_trees$AA3$tip.label[grepl("\\*",cafe_trees$AA3$tip.label) == T] <- "X"
plot.phylo(worked_trees$AA3,show.node.label = T,cex=0.7,edge.width=0.15)
title("AA3")

worked_trees$AA9$node.label[grepl("\\*",cafe_trees$AA9$node.label) == F] <- NA
worked_trees$AA9$node.label[grepl("\\*",cafe_trees$AA9$node.label) == T] <- "X"
worked_trees$AA9$tip.label[grepl("\\*",cafe_trees$AA9$tip.label) == F] <- NA
worked_trees$AA9$tip.label[grepl("\\*",cafe_trees$AA9$tip.label) == T] <- "X"
plot.phylo(worked_trees$AA9,show.node.label = T,cex=0.7,edge.width=0.15)
title("AA9")

worked_trees$GH16$node.label[grepl("\\*",cafe_trees$GH16$node.label) == F] <- NA
worked_trees$GH16$node.label[grepl("\\*",cafe_trees$GH16$node.label) == T] <- "X"
worked_trees$GH16$tip.label[grepl("\\*",cafe_trees$GH16$tip.label) == F] <- NA
worked_trees$GH16$tip.label[grepl("\\*",cafe_trees$GH16$tip.label) == T] <- "X"
plot.phylo(worked_trees$GH16,show.node.label = T,cex=0.7,edge.width=0.15)
title("GH16")

worked_trees$GH5$node.label[grepl("\\*",cafe_trees$GH5$node.label) == F] <- NA
worked_trees$GH5$node.label[grepl("\\*",cafe_trees$GH5$node.label) == T] <- "X"
worked_trees$GH5$tip.label[grepl("\\*",cafe_trees$GH5$tip.label) == F] <- NA
worked_trees$GH5$tip.label[grepl("\\*",cafe_trees$GH5$tip.label) == T] <- "X"
plot.phylo(worked_trees$GH5,show.node.label = T,cex=0.7,edge.width=0.15)
title("GH5")
dev.off()


