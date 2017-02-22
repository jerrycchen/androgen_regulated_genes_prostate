
library(affycoretools);
library(rgl);
library(limma);
library(illuminaMousev2.db);
library(annaffy);
library(GOstats);
library(WGCNA);
library(Heatplus);
library(statmod);
library(KEGGREST);
library(doBy);

##============STEP 0: Prepare the dataset============
Data_Input <- read.csv("results_PrEC1.2_PrSC1.3_of_model3.csv");
Data_Input <- Data_Input[ ,-1];

# table(Data_Input$Str_Code, Data_Input$Epi_Code);
#         -1     0     1
#   -1   145   543     6
#   0    202 26006   307
#   1      4   514   159

##============STEP 1: Add gene EntrezID and gene symbol============
# Based on Jenny's original code
# Read in annotation from Illumina (based on MouseWG-6_V2_0_R3_11278593_A.bgx dated 6/30/2010, which was read into GenomeStudio and output as .txt file)

annot <- read.delim("Illumina_MouseWG6v2_annot_02apr12.txt");
names(annot);
names(annot)[9:11] <- c("GO_CC","GO_BP","GO_MF");
annot[is.na(annot)] <- "";
sum(!Data_Input$ProbeID %in% annot$ProbeID); # 0, all probes are in annotation file
rownames(annot) <- annot$ProbeID;

Data_Input$ProbeID <- as.character(Data_Input$ProbeID);

Output <- cbind(annot[Data_Input$ProbeID, ],Data_Input);

# Save the probe file
write.csv(x = Output, file = "results_PrEC1.2_PrSC1.3_of_model3_with_entrez.csv", row.names = FALSE);
rm(annot);
rm(Data_Input);
rm(Output);


##============STEP 2: Collapsed into genewise report============
# Load the datafile
Data_Full <- read.csv("results_PrEC1.2_PrSC1.3_of_model3_with_entrez.csv", header=TRUE);

# Delete some of the variables
Data_Reduced <- Data_Full[ ,c(1,2,3,6,14,17,19,20,23,25)];

# Remove unknown probes which have no EntrezID
Data_Reduced <- Data_Reduced[!is.na(Data_Reduced$ENTREZ_GENE_ID), ]; # There are 20835 probes left;
Data_Reduced <- Data_Reduced[order(Data_Reduced$ENTREZ_GENE_ID), ];

# Group the probes into different folders
attach(Data_Reduced);
Data_Reduced$Group[Str_Code==0 & Epi_Code==0] <- "1_None";
Data_Reduced$Group[Str_Code==1 & Epi_Code==0] <- "2_Str_Up";
Data_Reduced$Group[Str_Code==0 & Epi_Code==1] <- "3_Epi_Up";
Data_Reduced$Group[Str_Code==1 & Epi_Code==1] <- "4_Both_Up";
Data_Reduced$Group[Str_Code==-1 & Epi_Code==0] <- "5_Str_Down";
Data_Reduced$Group[Str_Code==0 & Epi_Code==-1] <- "6_Epi_Down";
Data_Reduced$Group[Str_Code==-1 & Epi_Code==-1] <- "7_Both_Down";
Data_Reduced$Group[Str_Code * Epi_Code == -1 ] <- "8_Differential";

table(Data_Reduced$Group);


detach(Data_Reduced);

# Grouped into entrez:group
Summary_Gene <- summaryBy( Str_LogFC + Epi_LogFC ~ Group + ENTREZ_GENE_ID, FUN = c(mean, min, max), data = Data_Reduced);
write.csv(Summary_Gene, "summary_genelist_by_group_and_entrez.csv", row.names = FALSE);
table(Summary_Gene$Group);
# Group    # Gene
# 1        13676
# 2        259
# 3        213
# 4        117
# 5        402
# 6        134
# 7        101
# 8        8

sum(duplicated(Summary_Gene$ENTREZ_GENE_ID));
# 431 genes show up in different categories

# Focus on the conflict genes that show at least some androgen response
Summary_ARG <- Summary_Gene[which(Summary_Gene$Group!="1_None"), ];
sum(duplicated(Summary_ARG$ENTREZ_GENE_ID));
# 37 genes show up in different categories
Summary_ARG_conflict <- Summary_ARG[which(duplicated(Summary_ARG$ENTREZ_GENE_ID)), ];
List_conflict <- Summary_ARG_conflict$ENTREZ_GENE_ID;
Data_conflict <- Data_Full[which(Data_Full$ENTREZ_GENE_ID %in% List_conflict), c(6,3,1,2,14:25)];
Data_conflict <- Data_conflict[ ,c(1,2,10,16,3:9,11:15)]; # Reorder the output;
write.csv(x = Data_conflict, file = "summary_genelist_by_group_and_entrez_conflict_to_be_reviewed.csv", row.names = FALSE);

length(unique(Summary_Gene$ENTREZ_GENE_ID))  # 14479



