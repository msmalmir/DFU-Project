##############GSEA
library(devtools)
library(presto)


library(msigdbr)
library(fgsea)
library(dplyr)
library(ggplot2)
library(tidyverse)

dge<-read.csv('C:/Users/Wally/Documents/Research/DFU/GSEA/Non_Healing_gsea.csv')

msigdbr_show_species()
m_df<- msigdbr(species = "Homo sapiens", category = "C2",subcategory = 'CP:WIKIPATHWAYS')#subcategory = 'GO:BP'

head(m_df)

#write.csv(m_df,'C:/Users/Wally/Documents/LiverCancer/New-analysis/cc/Selected-G1/Cell_adhesion/BP.csv')

fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

LISTA<-fgsea_sets$GOBP_CELL_CELL_ADHESION

#dge$<-row.names(dge)
#row.names(dge)<-dg
big.genes<- dge %>%
  dplyr::select(X_index,AS)
ranks<- deframe(big.genes)

head(ranks)

fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 5000)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))


fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()




#fgseaResTidy1<- fgseaResTidy %>% filter(padj < 0.05) %>% head(n= 100)

fgseaResTidy$pathway <- factor(fgseaResTidy$pathway, levels = unique(fgseaResTidy$pathway[order(fgseaResTidy$NES)]))

ggplot(fgseaResTidy , aes(x = NES, y = pathway, color = padj, size = size)) + 
  geom_point(stat = 'identity') + 
  xlab("NES") + ylab("path") + ggtitle("GO BP: Control vs KO") + 
  theme_bw()+ theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

#plotEnrichment(fgsea_sets[["GOBP_CYTOPLASMIC_TRANSLATION"]],
#              ranks) + labs(title="GOBP_CYTOPLASMIC_TRANSLATION")

fwrite(fgseaResTidy, file ='C:/Users/Wally/Documents/LiverCancer/RNA-SEQ-AC-S/S-AC-GSEA-BP.csv')
#fwrite(fgseaResTidy, file ='C:/Users/Wally/Documents/LiverCancer/New-analysis/GSEA/Results_vs_control/GSEA-GO-BP-HIGH-c-s10.csv')





#################### ENRICHR
dge<-c('COL6A2', 'GSN', 'C1S', 'SFRP2', 'IGFBP4', 'COL6A3', 'SOCS3', 'EPS8',
       'MXRA8', 'CTSZ', 'COL4A2', 'MAF', 'CYBA', 'NOTCH3', 'PCOLCE', 'MMP2',
       'THBS2', 'FBLN2', 'GNG11', 'PRRX1', 'NBL1', 'COL5A2', 'VCAN', 'CAV1',
       'NFKBIA', 'LGI4', 'ASRGL1', 'MYL9')

library(devtools)
library(enrichR)
websiteLive <- getOption("enrichR.live")
if (websiteLive) {
  listEnrichrSites()
  setEnrichrSite("Enrichr") # Human genes   
}
if (websiteLive) dbs <- listEnrichrDbs()
dbs <- c( "MSigDB_Hallmark_2020")
if (websiteLive) {
  enriched <- enrichr(c(dge), dbs)
}
##

if (websiteLive) enriched[["MSigDB_Hallmark_2020"]]
if (websiteLive) {
  plotEnrich(enriched[[1]], showTerms = 60, numChar = 40, y = "Count", orderBy = "P.value")
}
bp <- enriched[["MSigDB_Hallmark_2020"]]
bp <- bp[bp$Adjusted.P.value <= 0.05, ]
bp <- bp[order(-bp$Odds.Ratio), ]

bp$Term <- factor(bp$Term, levels = unique(bp$Term[order(bp$Odds.Ratio)]))

ggplot(bp , aes(x = Odds.Ratio, y = Term, color = Adjusted.P.value, size = Combined.Score)) + 
  geom_point(stat = 'identity') + 
  xlab("Ratio") + ylab("path") + ggtitle("Hallmark: Healing") + 
  theme_bw()+ theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
write.csv(bp,'C:/Users/Wally/Documents/Research/DFU/GSEA/Healing/Hallmark3.csv')
#################################################################################################
websiteLive <- getOption("enrichR.live")
if (websiteLive) {
  listEnrichrSites()
  setEnrichrSite("Enrichr") # Human genes   
}
if (websiteLive) dbs <- listEnrichrDbs()
dbs <- c( "GO_Molecular_Function_2023")
if (websiteLive) {
  enriched <- enrichr(c(dge), dbs)
}
##

if (websiteLive) enriched[["GO_Molecular_Function_2023"]]
if (websiteLive) {
  plotEnrich(enriched[[1]], showTerms = 30, numChar = 40, y = "Count", orderBy = "P.value")
}
bp <- enriched[["GO_Molecular_Function_2023"]]
bp <- bp[bp$Adjusted.P.value <= 0.05, ]
bp <- bp[order(-bp$Odds.Ratio), ]

bp$Term <- factor(bp$Term, levels = unique(bp$Term[order(bp$Odds.Ratio)]))


# Plot using the new truncated Term column
ggplot(bp, aes(x = Odds.Ratio, y = Term, color = Adjusted.P.value, size = Combined.Score)) + 
  geom_point(stat = 'identity') + 
  xlab("Ratio") + 
  ylab("path") + 
  ggtitle("GO MF: Non-Healing") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
write.csv(bp,'C:/Users/Wally/Documents/Research/DFU/GSEA/Non_healing/MF3.csv')

#################################################################################################
websiteLive <- getOption("enrichR.live")
if (websiteLive) {
  listEnrichrSites()
  setEnrichrSite("Enrichr") # Human genes   
}
if (websiteLive) dbs <- listEnrichrDbs()
dbs <- c( "GO_Biological_Process_2023")
if (websiteLive) {
  enriched <- enrichr(c(dge), dbs)
}
##

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
if (websiteLive) {
  plotEnrich(enriched[[1]], showTerms = 30, numChar = 40, y = "Count", orderBy = "P.value")
}
bp <- enriched[["GO_Biological_Process_2023"]]
bp <- bp[bp$Adjusted.P.value <= 0.05, ]
bp <- bp[order(-bp$Odds.Ratio), ]

bp$Term <- factor(bp$Term, levels = unique(bp$Term[order(bp$Odds.Ratio)]))

ggplot(bp , aes(x = Odds.Ratio, y = Term, color = Adjusted.P.value, size = Combined.Score)) + 
  geom_point(stat = 'identity') + 
  xlab("Ratio") + ylab("path") + ggtitle("GO BP: Non-Healing") + 
  theme_bw()+ theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
write.csv(bp,'C:/Users/Wally/Documents/Research/DFU/GSEA/Non_healing/BP3.csv')

#################################################################################################
websiteLive <- getOption("enrichR.live")
if (websiteLive) {
  listEnrichrSites()
  setEnrichrSite("Enrichr") # Human genes   
}
if (websiteLive) dbs <- listEnrichrDbs()
dbs <- c( "WikiPathways_2013")
if (websiteLive) {
  enriched <- enrichr(c(dge), dbs)
}
##

if (websiteLive) enriched[["WikiPathways_2013"]]
if (websiteLive) {
  plotEnrich(enriched[[1]], showTerms = 30, numChar = 40, y = "Count", orderBy = "P.value")
}
bp <- enriched[["WikiPathways_2013"]]
bp <- bp[bp$Adjusted.P.value <= 0.05, ]
bp <- bp[order(-bp$Odds.Ratio), ]

bp$Term <- factor(bp$Term, levels = unique(bp$Term[order(bp$Odds.Ratio)]))

ggplot(bp , aes(x = Odds.Ratio, y = Term, color = Adjusted.P.value, size = Combined.Score)) + 
  geom_point(stat = 'identity') + 
  xlab("Ratio") + ylab("path") + ggtitle("GO WIKI: Healing") + 
  theme_bw()+ theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
write.csv(bp,'C:/Users/Wally/Documents/Research/DFU/GSEA/Healing/WIKI.csv')

#################################################################################################
websiteLive <- getOption("enrichR.live")
if (websiteLive) {
  listEnrichrSites()
  setEnrichrSite("Enrichr") # Human genes   
}
if (websiteLive) dbs <- listEnrichrDbs()
dbs <- c( "KEGG_2021_Human")
if (websiteLive) {
  enriched <- enrichr(c(dge), dbs)
}
##

if (websiteLive) enriched[["KEGG_2021_Human"]]
if (websiteLive) {
  plotEnrich(enriched[[1]], showTerms = 30, numChar = 40, y = "Count", orderBy = "P.value")
}
bp <- enriched[["KEGG_2021_Human"]]
bp <- bp[bp$Adjusted.P.value <= 0.05, ]
bp <- bp[order(-bp$Odds.Ratio), ]

bp$Term <- factor(bp$Term, levels = unique(bp$Term[order(bp$Odds.Ratio)]))

ggplot(bp , aes(x = Odds.Ratio, y = Term, color = Adjusted.P.value, size = Combined.Score)) + 
  geom_point(stat = 'identity') + 
  xlab("Ratio") + ylab("path") + ggtitle("KEGG: Non-Healing") + 
  theme_bw()+ theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
write.csv(bp,'C:/Users/Wally/Documents/Research/DFU/GSEA/Healing/KEGG.csv')


#################################################################################################
websiteLive <- getOption("enrichR.live")
if (websiteLive) {
  listEnrichrSites()
  setEnrichrSite("Enrichr") # Human genes   
}
if (websiteLive) dbs <- listEnrichrDbs()
dbs <- c( "Diabetes_Perturbations_GEO_2022")
if (websiteLive) {
  enriched <- enrichr(c(dge), dbs)
}
##

if (websiteLive) enriched[["Diabetes_Perturbations_GEO_2022"]]
if (websiteLive) {
  plotEnrich(enriched[[1]], showTerms = 30, numChar = 40, y = "Count", orderBy = "P.value")
}
bp <- enriched[["Diabetes_Perturbations_GEO_2022"]]
bp <- bp[bp$Adjusted.P.value <= 0.05, ]
bp <- bp[order(-bp$Odds.Ratio), ]

bp$Term <- factor(bp$Term, levels = unique(bp$Term[order(bp$Odds.Ratio)]))

ggplot(bp , aes(x = Odds.Ratio, y = Term, color = Adjusted.P.value, size = Combined.Score)) + 
  geom_point(stat = 'identity') + 
  xlab("Ratio") + ylab("path") + ggtitle("Diabetes Perturbations: Non-Healing") + 
  theme_bw()+ theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
write.csv(bp,'C:/Users/Wally/Documents/Research/DFU/GSEA/Non_healing/Diabetes3.csv')


#################################################################################################
websiteLive <- getOption("enrichR.live")
if (websiteLive) {
  listEnrichrSites()
  setEnrichrSite("Enrichr") # Human genes   
}
if (websiteLive) dbs <- listEnrichrDbs()
dbs <- c( "Reactome_2022")
if (websiteLive) {
  enriched <- enrichr(c(dge$X_index), dbs)
}
##

if (websiteLive) enriched[["Reactome_2022"]]
if (websiteLive) {
  plotEnrich(enriched[[1]], showTerms = 30, numChar = 40, y = "Count", orderBy = "P.value")
}
bp <- enriched[["Reactome_2022"]]
bp <- bp[bp$Adjusted.P.value <= 0.05, ]
bp <- bp[order(-bp$Odds.Ratio), ]

bp$Term <- factor(bp$Term, levels = unique(bp$Term[order(bp$Odds.Ratio)]))

ggplot(bp , aes(x = Odds.Ratio, y = Term, color = Adjusted.P.value, size = Combined.Score)) + 
  geom_point(stat = 'identity') + 
  xlab("Ratio") + ylab("path") + ggtitle("GO BP: Non-Healing") + 
  theme_bw()+ theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
write.csv(bp,'C:/Users/Wally/Documents/Research/DFU/GSEA/Non_healing/Reactome.csv')
