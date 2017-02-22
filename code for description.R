##============Step 0: Preparation============#

library(ggplot2);

Input <- read.csv("results_no_filter_of_model3.csv", header=TRUE);
Input <- Input[ ,-1];


##============Step 1: Histograms for response============#
# 1-1 LogFC in PrSC
x11();
h1 <- ggplot(Input, aes(Str_LogFC));
h1 + geom_histogram(bins = 1000, fill = "red", alpha = 0.2) +
     geom_vline(xintercept = log2(1), linetype = "longdash", colour = "grey" ) +
     geom_vline(xintercept = c(log2(1.2), -log2(1.2)), alpha = 0.5, linetype = "longdash", colour = "blue") +
     geom_vline(xintercept = c(log2(1.3), -log2(1.3)), alpha = 0.5, linetype = "longdash", colour = "purple") +
     xlim(-2, 2) +
     theme_grey();
ggsave("histogram_PrSC_LogFC_line1.2-1.3.png", plot = last_plot(), width = 6, height = 4);

# 1-2 LogFC in PrEC
x11();
h2 <- ggplot(Input, aes(Epi_LogFC));
h2 + geom_histogram(bins = 1000, fill = "cyan", alpha = 0.2) +
     geom_vline(xintercept = log2(1), linetype = "longdash", colour = "grey" ) +
     geom_vline(xintercept = c(log2(1.2), -log2(1.2)), alpha = 0.5, linetype = "longdash", colour = "blue") +
     xlim(-1, 1) +
     theme_grey();
ggsave("histogram_PrEC_LogFC_line1.2.png", plot = last_plot(), width = 6, height = 4);

# 1-3 LogFC of both in one graph

LogFC_both <- c(Input$Str_LogFC, Input$Epi_LogFC);
Cell_Type <- as.factor(c(rep("PrSC", 27886), rep("PrEC", 27886)));
Cell_Type <- relevel(Cell_Type, ref = "PrSC");
Temp_both <- data.frame(LogFC_both, Cell_Type);

x11();
h3 <- ggplot(Temp_both, aes(LogFC_both, fill = Cell_Type));
h3 + geom_histogram(aes(y=..density..), bins = 1000, alpha = 0.2, position = "identity") +
     xlim(-2, 2)
ggsave("histogram_LogFC_Both.png", plot = last_plot(), width = 6, height = 4);

# 1-4 LogFC of both stacked
x11();
h4 <- ggplot(Temp_both, aes(LogFC_both, fill = Cell_Type));
h4 + geom_histogram(aes(y=..density..), bins = 1000, alpha = 0.2) + 
     facet_grid(Cell_Type~.) + 
     geom_vline(xintercept = c(-log2(1.3), log2(1.3)), alpha = 0.5, colour = "red", linetype = "longdash") + 
     geom_vline(xintercept = c(-log2(1.2), log2(1.2)), alpha = 0.5, colour = "green", linetype = "longdash") + 
     xlim(-2, 2);
ggsave("histogram_LogFC_Both_line_1.2_1.3.png", plot = last_plot(), width = 6, height = 4);

##============Step 2: Scatter plot for response============#
# 2-1 No filter plot
Data_sp <- Input[ , c("Str_LogFC", "Epi_LogFC", "Str_Code", "Epi_Code")];  #Filter is not applied yet

attach(Data_sp);
Data_sp$Group[Str_Code==0 & Epi_Code==0] <- "1_None (24391)"
Data_sp$Group[Str_Code==1 & Epi_Code==0] <- "2_Str_Up (1055)"
Data_sp$Group[Str_Code==0 & Epi_Code==1] <- "3_Epi_Up (430)"
Data_sp$Group[Str_Code==1 & Epi_Code==1] <- "4_Both_Up (272)"
Data_sp$Group[Str_Code==-1 & Epi_Code==0] <- "5_Str_Down (1152)"
Data_sp$Group[Str_Code==0 & Epi_Code==-1] <- "6_Epi_Down (309)"
Data_sp$Group[Str_Code==-1 & Epi_Code==-1] <- "7_Both_Down (237)"
Data_sp$Group[Str_Code * Epi_Code == -1 ] <- "8_Differential (40)"

table(Str_Code, Epi_Code);
table(Group);

#         Epi_Code
#Str_Code    -1     0     1
#      -1   237  1152    19
#      0    309 24391   430
#      1     21  1055   272
# Group
# Both_Down      Both_Up Differential     Epi_Down       Epi_Up         None     Str_Down       Str_Up 
#       237          272           40          309          430        24391         1152         1055 

detach(Data_sp);

x11();
sp1 <- ggplot(Data_sp, aes(Epi_LogFC, Str_LogFC, colour = factor(Group)));
sp1 + geom_point() +
      xlim(-4, 4) +
      ylim(-4, 4) + 
      scale_color_manual(values = c("#CCCCCC", "#FF9900", "#FF6666", "#FF0000", "#66FF00","#33FFFF", "#3333FF", "#660099"))
ggsave("scatter_LogFC_filter1.0.png", plot = last_plot(), width = 8, height = 7);

# 2-2 Apply the filters
Data_sp2 <- Input[ , c("Str_LogFC", "Epi_LogFC", "Str_Code", "Epi_Code")];  #Filter is not applied yet

attach(Data_sp2);
Data_sp2$Str_Code[abs(Str_LogFC)<log2(1.3)] <- 0;
Data_sp2$Epi_Code[abs(Epi_LogFC)<log2(1.2)] <- 0;
detach(Data_sp2);


attach(Data_sp2);
Data_sp2$Group[Str_Code==0 & Epi_Code==0] <- "1_None (26006)";
Data_sp2$Group[Str_Code==1 & Epi_Code==0] <- "2_Str_Up (514)";
Data_sp2$Group[Str_Code==0 & Epi_Code==1] <- "3_Epi_Up (307)";
Data_sp2$Group[Str_Code==1 & Epi_Code==1] <- "4_Both_Up (159)";
Data_sp2$Group[Str_Code==-1 & Epi_Code==0] <- "5_Str_Down (543)";
Data_sp2$Group[Str_Code==0 & Epi_Code==-1] <- "6_Epi_Down (202)";
Data_sp2$Group[Str_Code==-1 & Epi_Code==-1] <- "7_Both_Down (145)";
Data_sp2$Group[Str_Code * Epi_Code == -1 ] <- "8_Differential (10)";

table(Data_sp2$Group);

# 1_None      2_Str_Up        3_Epi_Up      4_Both_Up     5_Str_Down     6_Epi_Down    7_Both_Down 
# 26006            514            307            159            543            202            145 
# 8_Differential 
# 10 

detach(Data_sp2);

x11();
sp2 <- ggplot(Data_sp2, aes(Epi_LogFC, Str_LogFC, colour = factor(Group)));
sp2 + geom_point() +
      xlim(-4, 4) +
      ylim(-4, 4) + 
      scale_color_manual(values = c("#CCCCCC", "#FF9900", "#FF6666", "#FF0000", "#66FF00","#33FFFF", "#3333FF", "#660099"))
ggsave("scatter_LogFC_Sti_filter1.3_Epi_filter_1.2.png", plot = last_plot(), width = 8, height = 7);


##============Step 3: Generate heatmap using raw data============#

library(Heatplus);

ARG_No_Filter <- subset(Input, Str_Code!=0 | Epi_Code!=0, select = "ProbeID");
List_No_Filter <- ARG_No_Filter[ ,1];

ARG_Filter_Used <- subset(Input, (abs(Str_LogFC)>log2(1.3) & Str_Code!=0) | 
                                 (abs(Epi_LogFC)>log2(1.2) & Epi_Code!=0),
                          select = "ProbeID");
List_Filter_Used <- ARG_Filter_Used[ ,1];

Expr_No_Filter <- read.csv("Filtered_Data.csv");
rownames(Expr_No_Filter) <- Expr_No_Filter[ ,1];
Expr_No_Filter <- Expr_No_Filter[ ,c(-1, -2)];
Expr_Heat_No_Filter <- Expr_No_Filter[as.character(List_No_Filter), ];
                                  

Expr_Filtered <- read.csv("Filtered_Data.csv");
rownames(Expr_Filtered) <- Expr_Filtered[ ,1];
Expr_Filtered <- Expr_Filtered[ ,c(-1, -2)];
Expr_Heat_Filtered <- Expr_Filtered[as.character(List_Filter_Used), ];

rm(ARG_Filter_Used);
rm(ARG_No_Filter);

x11();
png("heatmap_ARG_no_filter.png", width = 350, height = 350);
hmp1 <- regHeatmap(as.matrix(Expr_Heat_No_Filter));
hmp1;
plot(hmp1);
dev.off();


x11();
png("heatmap_ARG_PrSC1.3_PrEC1.2.png", width = 350, height = 350);
hmp2 <- regHeatmap(as.matrix(Expr_Heat_Filtered));
hmp2;
plot(hmp2);
dev.off();

##============END OF STEPS============



