library(GO.db);
library(topGO);
library(illuminaMousev2.db);

##============ STEP 0: Prepare the dataset ============##
Data_raw <- read.csv("results_PrEC1.2_PrSC1.3_of_model3_with_entrez.csv");
Data_raw <- Data_raw[ ,-c(4,5,7,8,12)];
Data_raw <- Data_raw[ ,-8];

attach(Data_raw);     # Add the group variable
Group <- c();
Group[Str_Code==0 & Epi_Code==0] <- "1_None";
Group[Str_Code==1 & Epi_Code==0] <- "2_Str_Up";
Group[Str_Code==0 & Epi_Code==1] <- "3_Epi_Up";
Group[Str_Code==1 & Epi_Code==1] <- "4_Both_Up";
Group[Str_Code==-1 & Epi_Code==0] <- "5_Str_Down";
Group[Str_Code==0 & Epi_Code==-1] <- "6_Epi_Down";
Group[Str_Code==-1 & Epi_Code==-1] <- "7_Both_Down";
Group[Str_Code * Epi_Code == -1 ] <- "8_Differential";
detach(Data_raw);

Data1 <- data.frame(Group, Data_raw);  # This is the file used in the future;
rownames(Data1) <- as.character(Data1$ProbeID);
colnames(Data1)[5:8] <- c("ENTREZ_ID", "CC", "BP", "MF");  # The GO terms here is useless;

Data2 <- Data1[!is.na(Data1$ENTREZ_ID), ];  # Remove the probes without Entrez ID;


##============ STEP 1: Edit the analysis parameters ============##

# Form the gene universe which must be a name vector;
Entrez_unique <- unique(Data2$ENTREZ_ID);  # This vector will be used as the name of gene universe;

# Modify the "which" function to select the list of Entrez ID in which the analysis will be performed;
Entrez_interest <- unique(Data2[which(Data2$Group 
                                      %in% c("3_Epi_Up", "4_Both_Up", "6_Epi_Down", "7_Both_Down", "8_Differential")), "ENTREZ_ID"]);

# Make the gene universe named vector;
Entrez_select <- as.integer(Entrez_unique %in% Entrez_interest);
Gene_universe <- Entrez_select;
names(Gene_universe) <- Entrez_unique;

# Define gene selection function;
Gene_selfun <- function(x) {x==1};  # loaded with the Entrez ID values;

# Set up the following parameters;
GO_type <- "BP";  # select in "BP", "MF" or "CC";
Gene_group <- "PrEC All ARGs";  # Describe the group of interests;
Node_size <- 10;  # minimum number of genes in each gene ontology terms;
Description <- paste(GO_type, "Analysis on", Gene_group);  # provide some information of the analysis;

Output_number <- 100;
Output_file <- paste("output/", Description, ".csv", sep = "");

##============ STEP 2: Form the proper topGO object and perform corresponding analysis ============##
# Form the object used for enrichment analysis;
Data_GO <- new("topGOdata",
               description = Description,
               ontology = GO_type,
               allGenes = Gene_universe,
               geneSelectionFun = Gene_selfun,
               nodeSize = Node_size,
               annotationFun = annFUN.org, ID = "entrez", mapping = "org.Mm.eg");

# Do the corresponding GO analysis;
Result_GO <- runTest(Data_GO, algorithm = "classic", statistic = "fisher");

# Output the analysis result;
Output_GO <- GenTable(Data_GO, classicFisher = Result_GO, topNodes = Output_number);
write.csv(Output_GO, file = Output_file, row.names = FALSE);
# write.csv(Output_GO, file = other output options, row.names = FALSE);
