#20240123 pathway orphans re-analysis

library(tidyverse)
library(clusterProfiler)
library(rtracklayer)
library(wesanderson)

####Pathway databases download and analysis####

#hgnc updated January 22, 2024 with 19,258 protein-coding genes
prot.cod <- read.delim("gene_with_protein_product.txt")

#read and prepare GO datasets from msigdb human
#October 2023 release
#genes are in entrez id form
hall <- read.gmt("h.all.v2023.2.Hs.entrez.gmt")
c2 <- read.gmt("c2.all.v2023.2.Hs.entrez.gmt")
msigdb.subset <- rbind(hall, c2)
colnames(msigdb.subset) <- c("Pathway", "Entrez_ID") #561242 entries

#reactome v87 updated December 4, 2023
reactome <- read.delim("NCBI2Reactome.txt", header = FALSE)
reactome %>% dplyr::filter(V6 == "Homo sapiens") -> reactome
reactome %>% dplyr::select(V4, V1) -> reactome
colnames(reactome) <- c("Pathway", "Entrez_ID")

#wikipathways updated January 11, 2024
wiki <- read.gmt("wikipathways-20240101-gmt-Homo_sapiens.gmt")
colnames(wiki) <- c("Pathway", "Entrez_ID")

#kegg release 109 January 1, 2024
kegg <- read_csv("hsa.csv")

#join pathway database entries
pathways <- rbind(msigdb.subset, reactome, wiki, kegg)

#count pathway hits per gene

pathways %>% group_by(Entrez_ID) %>% summarize(Freq=n()) -> pathway.count
pathway.count[!is.na(as.numeric(pathway.count$Entrez_ID)), ] -> pathway.count
pathway.count$Entrez_ID <- as.numeric(pathway.count$Entrez_ID)

left_join(prot.cod, pathway.count, by = c("entrez_id" = "Entrez_ID")) -> prot.cod.pathway.count
sum(is.na(prot.cod.pathway.count$Freq)) #633 pathway orphans based on these 4 databases


#gene ontology
go <- read_csv("GO_BPMF.csv")

length(intersect(prot.cod.pathway.count$symbol, go$DB_Object_Symbol))

prot.cod.pathway.count %>% mutate(GO = case_when(symbol %in% go$DB_Object_Symbol ~ "Y",
                                                 !(symbol %in% go$DB_Object_Symbol) ~ "N")) -> prot.cod.pathway.count

prot.cod.pathway.count$Freq <- as.numeric(prot.cod.pathway.count$Freq)
sum(is.na(prot.cod.pathway.count$Freq))
sum(is.na(prot.cod.pathway.count$Freq) & prot.cod.pathway.count$GO == "Y")


#CPDB
cpdb.sym = read.delim("0331_CPDB_pathways_genes.tab")

#count pathway hits per gene
pathway_occurences = c()
for (i in 1:nrow(prot.cod)) {
  gene = prot.cod[i,"symbol"]
  count = sum(grepl(gene, cpdb.sym$hgnc_symbol_ids))
  pathway_occurences = c(pathway_occurences, count)
}

prot.cod.pathway.count$cpdb_occurences <- pathway_occurences

#check how many protein-coding genes have no pathway hits
sum(is.na(prot.cod.pathway.count$Freq) & prot.cod.pathway.count$GO == "N" & prot.cod.pathway.count$cpdb_occurences == 0)

#publication count
gene2pubmed <- read_csv("gene2pubmed_hs.csv")

left_join(prot.cod.pathway.count, gene2pubmed, by = c("entrez_id" = "Entrez_ID")) -> prot.cod.pathway.count

length(setdiff(prot.cod.pathway.count$entrez_id, gene2pubmed$Entrez_ID))

#pathway orphan identification

prot.cod.pathway.count %>% dplyr::filter(is.na(prot.cod.pathway.count$Freq) & prot.cod.pathway.count$GO == "N" & prot.cod.pathway.count$cpdb_occurences == 0) -> orphans

#save info on pathway hits
write.csv(prot.cod.pathway.count, "prot_cod_pathway_count.csv")
write.csv(orphans, "orphans.csv")


#housekeeping genes
housekeeping <- read_delim("housekeeping_info.txt")


#classify protein-coding genes
prot.cod.pathway.count %>% mutate(Type = case_when(entrez_id %in% orphans$entrez_id ~ "Pathway orphan",
                                                   entrez_id %in% housekeeping$`NCBI gene ID` ~ "Housekeeping",
                                                   TRUE ~ "Non-pathway orphan")) -> prot.cod.pathway.count



prot.cod.pathway.count %>% dplyr::filter(is.na(prot.cod.pathway.count$Freq) & prot.cod.pathway.count$GO == "N" & prot.cod.pathway.count$cpdb_occurences == 0) -> orphans.info

####Characterize pathway orphans####

#Process gencode data for gene length

gencode <- readGFF("gencode.v45.basic.annotation.gtf")

gencode$gene_id <- str_extract(gencode$gene_id, ".*(?=\\.)")

gencode %>% dplyr::filter(type == "gene") -> gencode

gencode %>% mutate(gene_length = end - start) -> gencode

write_csv(gencode, "gencode.csv")

gencode %>% dplyr::select(gene_id, gene_length) -> gencode

left_join(prot.cod.pathway.count, gencode, by = c("ensembl_gene_id" = "gene_id")) -> prot.cod.pathway.count


#Process UniProt data for protein length and mass
uniprot <- read_csv("uniprot_nodup.csv")

left_join(prot.cod.pathway.count, uniprot, by = c("entrez_id" = "From")) -> prot.cod.pathway.count


#Chromosomal location data
chr <- read_delim("chromosome.txt")

left_join(prot.cod.pathway.count, chr, by = c("hgnc_id" = "HGNC ID")) -> prot.cod.pathway.count

table(prot.cod.pathway.count$Chromosome) %>% write.csv("chr.csv")

#Kidney expression data

kidney.exp <- read_csv("tpm_consensus_kidney.csv")

left_join(prot.cod.pathway.count, kidney.exp, by = c("ensembl_gene_id" = "Gene")) -> prot.cod.pathway.count


#Kidney enrichment data

kidney.z <- read_csv("gtex_zscore_dianne_kidney.csv")
distinct(kidney.z, Description, .keep_all = TRUE)
kidney.z <- kidney.z[!duplicated(kidney.z$Description),]

left_join(prot.cod.pathway.count, kidney.z, by = c("ensembl_gene_id" = "Description")) -> prot.cod.pathway.count

####Pathway orphans in the kidney####

#Identify which pathway orphans are expressed in the kidney

gtex.tpm <- read_delim(file = "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct", skip = 2)
gtex.tpm %>% dplyr::select("Name", "Description", "Kidney - Cortex", "Kidney - Medulla") -> gtex.tpm.kidney
gtex.tpm.kidney$Name <- str_extract(gtex.tpm.kidney$Name, ".*(?=\\.)")

gtex.tpm.kidney <- gtex.tpm.kidney[!duplicated(gtex.tpm.kidney$Description),]
left_join(prot.cod.pathway.count, gtex.tpm.kidney, by = c("ensembl_gene_id" = "Name")) -> prot.cod.pathway.count

hpa <- read_tsv("rna_tissue_hpa.tsv")
hpa %>% dplyr::filter(Tissue == "kidney") -> hpa
left_join(prot.cod.pathway.count, hpa, by = c("ensembl_gene_id" = "Gene")) -> prot.cod.pathway.count

write.csv(prot.cod.pathway.count, "prot_cod_info.csv")

prot.cod.pathway.count <- read_csv("prot_cod_info.csv")
write.csv(orphans.info, "orphans_info.csv")

#get prot cod in kidney

gtex.cortex.all <- read_delim(file = "gene_tpm_kidney_cortex.gct", skip = 2)
gtex.medulla.all <- read_delim(file = "gene_tpm_kidney_medulla.gct", skip = 2)


gtex.cortex <- gtex.cortex.all[rowSums(gtex.cortex.all[,c(4:88)] > 0) >= 63,]

gtex.medulla <- gtex.medulla.all[rowSums(gtex.medulla.all[,c(4:7)] > 0) >= 3,]

gtex.cortex$Name <- str_extract(gtex.cortex$Name, ".*(?=\\.)")
gtex.medulla$Name <- str_extract(gtex.medulla$Name, ".*(?=\\.)")

prot.cod.info %>% dplyr::filter(!is.na(ensembl_gene_id)) -> protcod.w
prot.cod.info %>% dplyr::filter(is.na(ensembl_gene_id)) -> protcod.wo

gtex.cortex %>% dplyr::filter(Name %in% protcod.w$ensembl_gene_id) -> protcod.w.cortex
gtex.medulla %>% dplyr::filter(Name %in% protcod.w$ensembl_gene_id) -> protcod.w.medulla

gtex.cortex %>% dplyr::filter(Description %in% protcod.wo$symbol) -> protcod.wo.cortex
gtex.medulla %>% dplyr::filter(Description %in% protcod.wo$symbol) -> protcod.wo.medulla

protcod.w.cortex %>% dplyr::select(2,3) -> cortex
protcod.w.medulla %>% dplyr::select(2,3) -> medulla

rbind(cortex, medulla) -> kidney.gtex

kidney.gtex <- kidney.gtex[!duplicated(kidney.gtex$Name),]

length(intersect(orphans.info$ensembl_gene_id, kidney.gtex$Name))

prot.cod.info %>% dplyr::filter(ensembl_gene_id %in% kidney.gtex$Name) -> protcod.kidney

write.csv(protcod.kidney, "protcod_kidney.csv")

length(intersect(orphans.info$ensembl_gene_id, cortex$Name))
length(intersect(orphans.info$ensembl_gene_id, medulla$Name))

orphans.info %>% dplyr::filter(ensembl_gene_id %in% kidney.gtex$Name) -> orphans.kidney
orphans.info %>% dplyr::filter(ensembl_gene_id %in% cortex$Name) -> orphans.cortex
orphans.info %>% dplyr::filter(ensembl_gene_id %in% medulla$Name) -> orphans.medulla

protcod.kidney %>% dplyr::filter(Type == "Non-pathway orphan") -> nonorp.kidney
protcod.kidney %>% dplyr::filter(Type == "Housekeeping") -> hk.kidney

summary(hk.kidney$gene_length)
summary(nonorp.kidney$gene_length)
summary(orphans.kidney$gene_length)

summary(hk.kidney$Length)
summary(nonorp.kidney$Length)
summary(orphans.kidney$Length)

summary(hk.kidney$Mass)
summary(nonorp.kidney$Mass)
summary(orphans.kidney$Mass)

protcod.kidney$nTPM <- as.numeric(protcod.kidney$nTPM)


#############



prot.cod.pathway.count %>% dplyr::filter(Type == "Housekeeping") -> hk
prot.cod.pathway.count %>% dplyr::filter(Type == "Non-pathway orphan") -> nonpo

summary(orphans$gene_length)
summary(hk$gene_length)
summary(nonpo$gene_length)

summary(orphans$Length)
summary(hk$Length)
summary(nonpo$Length)

#Uniprot data

uniprotall <- read_tsv("uniprotredo.tsv")
distinct(uniprotall, From, .keep_all = TRUE) -> uniprotclean

left_join(prot.cod.pathway.count, uniprotclean, by = c("hgnc_id" = "From")) -> prot.cod.info
summary(prot.cod.info$Length.y)

write_csv(prot.cod.info, "prot_cod_info.csv")
prot.cod.info <- read_csv("prot_cod_info.csv")

prot.cod.info %>% dplyr::filter(Type == "Housekeeping") -> hk
prot.cod.info %>% dplyr::filter(Type == "Pathway orphan") -> orp
prot.cod.info %>% dplyr::filter(Type == "Non-pathway orphan") -> nonorp

summary(orp$gene_length)
summary(nonorp$gene_length)
summary(hk$gene_length)
summary(orp$Length.y)
summary(nonorp$Length.y)
summary(hk$Length.y)
summary(orp$Mass.y)
summary(nonorp$Mass.y)
summary(hk$Mass.y)

####Plotting####

chr.all <- read_csv("chr_work.csv")

ggplot(data = chr.all, aes(fill = name, x = Chromosome, y = `Gene Count` )) + 
  geom_bar(position = "stack", stat = "identity") +
  theme_classic() + scale_fill_manual(values = wes_palette("Cavalcanti1")) + 
  theme(legend.position = "top", axis.text = element_text(size = 20), axis.title = element_text(size=20)) +
  theme(legend.text = element_text(size=20)) + 
  labs(fill="")

chr.orp <- read_csv("chr_orphans.csv")
chr.orp$Chromosome <- factor(chr.orp$Chromosome, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12",
                                                            "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"))


ggplot(chr.orp) + geom_bar(aes(x=Chromosome, y=`Gene count`), stat = "identity", fill = "#3B9AB2") +
  theme_classic() + 
  theme(legend.position = "top", axis.text = element_text(size = 30), axis.title = element_text(size=30)) +
  theme(legend.text = element_text(size=30)) + 
  labs(fill="")

dev.copy(png,"chr_orp.png", width=1200, height=600)
dev.off()

ggplot(data = prot.cod.info) + geom_boxplot(aes(x=Type, y=log10(Length), fill=Type)) + scale_fill_manual(values = c("#3B9AB2", "#EBCC2A", "#F21A00")) +
  theme_classic() + 
  theme(legend.position = "none", axis.text = element_text(size = 20), axis.title = element_text(size=30)) +
  labs(fill="") + scale_y_continuous("log10(Protein length in AA)") +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
dev.copy(png,"prot_length2.png", width=600, height=600)
dev.off()

kruskal.test(prot.cod.info$Length.y, prot.cod.info$Type)

ggplot(data = prot.cod.info) + geom_boxplot(aes(x=Type, y=log10(Mass), fill=Type)) + scale_fill_manual(values = c("#3B9AB2", "#EBCC2A", "#F21A00")) +
  theme_classic() + 
  theme(legend.position = "none", axis.text = element_text(size = 20), axis.title = element_text(size=30)) +
  labs(fill="") + scale_y_continuous("log10(Protein mass in Da)") +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
dev.copy(png,"prot_mass2.png", width=600, height=600)
dev.off()
kruskal.test(prot.cod.info$Mass.y, prot.cod.info$Type)

ggplot(data = prot.cod.info) + geom_boxplot(aes(x=Type, y=log10(gene_length), fill=Type)) + scale_fill_manual(values = c("#3B9AB2", "#EBCC2A", "#F21A00")) +
  theme_classic() + 
  theme(legend.position = "none", axis.text = element_text(size = 20), axis.title = element_text(size=30)) +
  labs(fill="") + scale_y_continuous("log10(Gene length in bp)") +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
dev.copy(png,"gene_length2.png", width=600, height=600)
dev.off()
kruskal.test(prot.cod.info$gene_length, prot.cod.info$Type)

ggplot(data = prot.cod.info) + geom_boxplot(aes(x=Type, y=log10(Publications), fill=Type)) + scale_fill_manual(values = c("#3B9AB2", "#EBCC2A", "#F21A00")) +
  theme_classic() + 
  theme(legend.position = "none", axis.text = element_text(size = 20), axis.title = element_text(size=30)) +
  labs(fill="") + scale_y_continuous("log10(Publications)") +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
dev.copy(png,"publications.png", width=600, height=600)
dev.off()


####DeepLoc####

deeploc <- read_csv("deeploc.csv")

deeploc %>% group_by(Localizations) %>% summarize(Freq=n())  %>% View()
table(deeploc$Localizations) %>% View()
deeploc %>% group_by(Localizations) %>% summarize(Freq=n()) %>% write_csv("deeploc_unproc.csv")
deeploc.summary <- read_csv("deeploc_summary.csv")

ggplot(deeploc.summary) + geom_bar(aes(x=Localizations, y=Freq), stat = "identity", fill = "#3B9AB2") +
  theme_classic() + scale_y_continuous("Pathway orphan count") +
  theme(legend.position = "top", axis.text = element_text(size = 20), axis.title = element_text(size=20)) +
  theme(legend.text = element_text(size=20)) + 
  labs(fill="") + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))

dev.copy(png,"sub.png", width=600, height=600)
dev.off()

####InterPro####


a <- read.delim("orphans1.tsv", header = FALSE)
b <- read.delim("orphans2.tsv", header = FALSE)
c <- read.delim("orphans3.tsv", header = FALSE)
rbind(a,b,c) -> interpro

interpro <- subset(interpro, V4!= "MobiDBLite")
interpro %>% dplyr::filter(V4 == "Pfam") -> d
d <- d[!duplicated(d$V1),]
d %>% group_by(V6) %>% summarize(Freq=n()) -> e
write_csv(e, "interpro_summary.csv")
interpro.summary <- read_csv("interpro_summary.csv")

####Plotting####

prot.cod.info %>% dplyr::filter(entrez_id %in% kidney$entrez_id) -> prot.cod.info.kidney

ggplot(data = prot.cod.info.kidney) + geom_boxplot(aes(x=Type, y=log10(nTPM.y), fill=Type)) + scale_fill_manual(values = c("#3B9AB2", "#EBCC2A", "#F21A00")) +
  theme_classic() + 
  theme(legend.position = "none", axis.text = element_text(size = 20), axis.title = element_text(size=30)) +
  labs(fill="") + scale_y_continuous("log10(Expression in TPM)") +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
dev.copy(png,"exp_hpa.png", width=600, height=600)
dev.off()

ggplot(data = prot.cod.info.kidney) + geom_boxplot(aes(x=Type, y=log10(`Kidney - Cortex.y`), fill=Type)) + scale_fill_manual(values = c("#3B9AB2", "#EBCC2A", "#F21A00")) +
  theme_classic() + 
  theme(legend.position = "none", axis.text = element_text(size = 20), axis.title = element_text(size=30)) +
  labs(fill="") + scale_y_continuous("log10(Expression in TPM)") +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
dev.copy(png,"exp_cortex.png", width=600, height=600)
dev.off()

ggplot(data = prot.cod.info.kidney) + geom_boxplot(aes(x=Type, y=log10(`Kidney - Medulla.y`), fill=Type)) + scale_fill_manual(values = c("#3B9AB2", "#EBCC2A", "#F21A00")) +
  theme_classic() + 
  theme(legend.position = "none", axis.text = element_text(size = 20), axis.title = element_text(size=30)) +
  labs(fill="") + scale_y_continuous("log10(Expression in TPM)") +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
dev.copy(png,"exp_medulla.png", width=600, height=600)
dev.off()

ggplot(data = prot.cod.info.kidney) + geom_boxplot(aes(x=Type, y=log10(pTPM), fill=Type)) + scale_fill_manual(values = c("#3B9AB2", "#EBCC2A", "#F21A00")) +
  theme_classic() + 
  theme(legend.position = "none", axis.text = element_text(size = 20), axis.title = element_text(size=30)) +
  labs(fill="") + scale_y_continuous("log10(Expression in pTPM)") +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
dev.copy(png,"exp_prot.png", width=600, height=600)
dev.off()

ggplot(data = prot.cod.info.kidney) + geom_density(aes(x=`Kidney - Cortex.x`, fill = Type), alpha = 0.75) +
  theme_classic() + scale_fill_manual(values = c("#3B9AB2", "#EBCC2A", "#F21A00")) +
  xlab("Cortex enrichment (z-score)") + theme(legend.position = c(.8, .9)) + 
  ylab("Density") + labs(fill="") + 
  theme(axis.text = element_text(size = 20), axis.title = element_text(size=20), legend.text = element_text(size = 20))

dev.copy(png,"zscore_cortex.png", width=600, height=600)
dev.off()

ggplot(data = prot.cod.info.kidney) + geom_density(aes(x=`Kidney - Medulla.x`, fill = Type), alpha = 0.75) +
  theme_classic() + scale_fill_manual(values = c("#3B9AB2", "#EBCC2A", "#F21A00")) +
  xlab("Medulla enrichment (z-score)") + theme(legend.position = c(.8, .9)) + 
  ylab("Density") + labs(fill="") + 
  theme(axis.text = element_text(size = 20), axis.title = element_text(size=20), legend.text = element_text(size = 20))

dev.copy(png,"zscore_medulla.png", width=600, height=600)
dev.off()

prot.cod.info.kidney %>% dplyr::filter(Type == "Pathway orphan") -> orp
prot.cod.info.kidney %>% dplyr::filter(Type == "Housekeeping") -> hk
prot.cod.info.kidney %>% dplyr::filter(Type == "Non-pathway orphan") -> nonorp

summary(orp$nTPM.y)
summary(hk$nTPM.y)
summary(nonorp$nTPM.y)

summary(orp$`Kidney - Cortex.y`)
summary(hk$`Kidney - Cortex.y`)
summary(nonorp$`Kidney - Cortex.y`)

summary(orp$`Kidney - Medulla.y`)
summary(hk$`Kidney - Medulla.y`)
summary(nonorp$`Kidney - Medulla.y`)
