# ### Phylogenetic analyses
library(phytools)
library(castor)
library(phangorn)
library(seqinr)

### Visualization
library(ggplot2)
library(ggtree)
library(ggtreeExtra)
library(gggenomes)
library(patchwork)
library(RColorBrewer)
library(ggnewscale)
library(ComplexHeatmap)
library(circlize)

### Miscellaneous
library(data.table)
library(tidyverse)
library(igraph)

# ### Statistics
library(rstatix)
library(ggpubr)

setwd("/Users/Chainorato/Desktop/Poxvirus_evol_project")

############################## Genome statistics of poxviruses ##############################
########## Gather metadata of poxviruses from NCBI
poxv_ncbi_ds.f <- "241008_Poxviridae_and_representative_MPXV_NCBI_datasets_taxonomy.tsv"
poxv_ncbi_ds.df <- fread(poxv_ncbi_ds.f, header=T, sep="\t", quote="", check.names=T, data.table=FALSE)
poxv_ncbi_ds.df <- poxv_ncbi_ds.df %>% mutate(tip_label = paste(Assembly.Accession, gsub(" ", "_", Organism.Name), sep="_"))
poxv_ncbi_ds.df$tip_label <- ifelse(poxv_ncbi_ds.df$tip_label=="GCF_000857045.1_Monkeypox_virus", "GCF_000857045.1_Monkeypox_virus_Zaire-96-I-16_Ia",
                             ifelse(poxv_ncbi_ds.df$tip_label=="GCA_039269415.1_Monkeypox_virus", "GCA_039269415.1_Monkeypox_virus_RDC-NKV-GOM-MPOX-010_Ib",
                             ifelse(poxv_ncbi_ds.df$tip_label=="GCA_006465585.1_Monkeypox_virus", "GCA_006465585.1_Monkeypox_virus_Sierra_Leone_IIa",
                             ifelse(poxv_ncbi_ds.df$tip_label=="GCF_014621545.1_Monkeypox_virus", "GCF_014621545.1_Monkeypox_virus_MPXV-M5312_HM12_Rivers_IIb", poxv_ncbi_ds.df$tip_label))))
poxv_ncbi_ds.df <- poxv_ncbi_ds.df %>% mutate(tip_label_dummy=str_split(tip_label, ".1_", simplify = TRUE)[,2], clade=Genus)
poxv_ncbi_ds.df <- poxv_ncbi_ds.df %>% arrange(Subfamily, clade)

poxv_color.clade.v <- c(brewer.pal(9, "Set1")[1:5], brewer.pal(9, "Set1")[7:9], brewer.pal(8, "Accent")[1:3], brewer.pal(8, "Accent")[5:7], brewer.pal(8, "Dark2"))
names(poxv_color.clade.v) <- unique(poxv_ncbi_ds.df$clade)

########## Gather genome statistics of poxviruses
poxv_genome_stats.f <- "251117_genes_per_genome_length.tsv"
poxv_genome_stats.df <- fread(poxv_genome_stats.f, header=T, sep="\t", quote="", check.names=T, data.table=FALSE)
poxv_genome_stats.df$species_name <- ifelse(poxv_genome_stats.df$species_name=="GCF_000857045.1_Monkeypox_virus", "GCF_000857045.1_Monkeypox_virus_Zaire-96-I-16_Ia",
                                     ifelse(poxv_genome_stats.df$species_name=="GCA_039269415.1_Monkeypox_virus", "GCA_039269415.1_Monkeypox_virus_RDC-NKV-GOM-MPOX-010_Ib",
                                     ifelse(poxv_genome_stats.df$species_name=="GCA_006465585.1_Monkeypox_virus", "GCA_006465585.1_Monkeypox_virus_Sierra_Leone_IIa",
                                     ifelse(poxv_genome_stats.df$species_name=="GCF_014621545.1_Monkeypox_virus", "GCF_014621545.1_Monkeypox_virus_MPXV-M5312_HM12_Rivers_IIb", poxv_genome_stats.df$species_name))))
poxv_genome_stats.df <- poxv_genome_stats.df %>% mutate(accession_no=paste(str_split(assembly_name, ".1_", simplify=TRUE)[,1], ".1", sep=""))
poxv_genome_stats.df <- poxv_genome_stats.df %>% mutate(Subfamily=poxv_ncbi_ds.df[match(species_name, poxv_ncbi_ds.df$tip_label), "Subfamily"], 
                                                        clade=poxv_ncbi_ds.df[match(species_name, poxv_ncbi_ds.df$tip_label), "clade"])
poxv_genome_stats.df <- poxv_genome_stats.df %>% mutate(species_name=str_split(species_name, paste(str_split(assembly_name, "_", simplify=TRUE)[,1], "_", str_split(assembly_name, "_", simplify=TRUE)[,2], "_", sep=""), simplify=TRUE)[,2])
poxv_genome_stats.df <- poxv_genome_stats.df %>% arrange(Subfamily, clade, species_name)
poxv_genome_stats.df$species_name <- factor(poxv_genome_stats.df$species_name, levels=poxv_genome_stats.df$species_name)
poxv_genome_stats.df$clade <- factor(poxv_genome_stats.df$clade, levels=unique(poxv_genome_stats.df$clade))


########## ***** Figure 1B: Correlation between gene number and genome size *****
poxv_gene_no_vs_genome_size_cor.res <- cor.test(poxv_genome_stats.df$gene_number, poxv_genome_stats.df$genome_length, method = "spearman")
poxv_gene_no_vs_genome_size_cor.res

poxv_gene_no_vs_genome_size.p <- ggplot(poxv_genome_stats.df, aes(x=genome_length, y=gene_number, color=clade))
poxv_gene_no_vs_genome_size.p <- poxv_gene_no_vs_genome_size.p + geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "black")
poxv_gene_no_vs_genome_size.p <- poxv_gene_no_vs_genome_size.p + geom_point()
poxv_gene_no_vs_genome_size.p <- poxv_gene_no_vs_genome_size.p + geom_text(aes(label=substr(clade, 1, 3)), hjust=-0.13, vjust=0, size=7)
poxv_gene_no_vs_genome_size.p <- poxv_gene_no_vs_genome_size.p + scale_color_manual("clade", name="Clade", values=poxv_color.clade.v, breaks=names(poxv_color.clade.v), na.value = NA)
poxv_gene_no_vs_genome_size.p <- poxv_gene_no_vs_genome_size.p + labs(title="Relationship between gene number and genome size", 
                                                                            subtitle=paste("r = ", round(poxv_gene_no_vs_genome_size_cor.res$estimate, 2), " ; P-value = ", 
                                                                                           sprintf("%.2e", poxv_gene_no_vs_genome_size_cor.res$p.value), sep=""), x="Virus genome size", y="Virus gene number")
poxv_gene_no_vs_genome_size.p <- poxv_gene_no_vs_genome_size.p + guides(color="none")
poxv_gene_no_vs_genome_size.p <- poxv_gene_no_vs_genome_size.p + theme_classic()
poxv_gene_no_vs_genome_size.p <- poxv_gene_no_vs_genome_size.p + theme(text=element_text(size = 24), plot.title=element_text(size = 32), axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
poxv_gene_no_vs_genome_size.p <- poxv_gene_no_vs_genome_size.p + scale_x_continuous(breaks=seq(120000, 360000, by=30000), limits=c(120000,360000))
poxv_gene_no_vs_genome_size.p <- poxv_gene_no_vs_genome_size.p + scale_y_continuous(limits=c(100,350))
poxv_gene_no_vs_genome_size.p

# pdf("Figures/gene_number_vs_genome_size.pdf", width=16, height=12)
# print(poxv_gene_no_vs_genome_size.p)
# dev.off()


########## Gather information of virus-host relationship
virus_host_db.f <- "250327_VirusHostDB/virushostdb.tsv"
virus_host_db.df <- fread(virus_host_db.f, header=T, sep="\t", quote="", check.names=T, data.table=FALSE)

virus_host_db.filtered.df <- virus_host_db.df %>% mutate(virus.name=ifelse(virus.name=="Parapoxvirus red deer/HL953", "Parapoxvirus red deer HL953",
                                                                    ifelse(virus.name=="Choristoneura rosaceana entomopoxvirus 'L'", "Choristoneura rosaceana entomopoxvirus L",
                                                                    ifelse(virus.name=="Mythimna separata entomopoxvirus 'L'", "Mythimna separata entomopoxvirus L",
                                                                    ifelse(virus.name=="Adoxophyes honmai entomopoxvirus 'L'", "Adoxophyes honmai entomopoxvirus L", virus.name)))))
virus_host_db.filtered.df <- virus_host_db.filtered.df %>% filter(virus.name %in% poxv_ncbi_ds.df$Organism.Name, evidence!="UniProt") %>% select(virus.tax.id, virus.name, refseq.id, host.tax.id, host.name, pmid, evidence)

virus_host_db.filtered.f <- "250327_VirusHostDB/virushostdb_filtered.tsv"
# write.table(virus_host_db.filtered.df, virus_host_db.filtered.f, col.names=T, row.names=F, sep="\t", quote=F)

# virus_host_db.filtered.edited.f <- "250327_VirusHostDB/virushostdb_filtered_edited.tsv"
# virus_host_db.filtered.edited.df <- read.delim(virus_host_db.filtered.edited.f, check.names = FALSE)

##### Manually edit host information as follows
### 1) Remove hosts that were inferred from experimental infection or using sentinel animals
### 2) Replace hosts whose genome statistics are unavailable with their most closely related animals whose genome statistics are available

##### Gather the manually curated information of virus-host relationship
virus_host_db.filtered.edited.final.f <- "250327_VirusHostDB/virushostdb_filtered_edited_final.tsv"
virus_host_db.filtered.edited.final.df <- read.delim(virus_host_db.filtered.edited.final.f, check.names = FALSE)

host_tax_id.ncbi.v <- sort(unique(c(virus_host_db.filtered.edited.final.df %>% filter(!is.na(host.tax.id.edited.genome_size)) %>% pull(host.tax.id.edited.genome_size),
                                    virus_host_db.filtered.edited.final.df %>% filter(!is.na(host.tax.id.edited.inferred)) %>% pull(host.tax.id.edited.inferred))))
host_tax_id.ncbi.cmm_sep <- paste(host_tax_id.ncbi.v, collapse=",")

##### Create an indexing file for retrieving genome statistics of hosts collected from DNAzoo database
host_tax_id.dna_zoo.f <- "250327_VirusHostDB/260219_virushostdb_dna_zoo.tsv"
host_tax_id.dna_zoo.df <- data.frame(species_name=gsub(" ", "_", sort(unique(virus_host_db.filtered.edited.final.df %>% filter(DNA_Zoo==TRUE) %>% pull(host.name.edited.genome_size)))))
host_tax_id.dna_zoo.df <- host_tax_id.dna_zoo.df %>% mutate(assembly_name=species_name)
# write.table(host_tax_id.dna_zoo.df, host_tax_id.dna_zoo.f, col.names=T, row.names=F, sep="\t", quote=F) # Manually edit assembly name

##### Gather genome statistics of hosts of poxviruses
host_genome_info.f <- "260220_host_genes_per_genome_length.tsv"
host_genome_info.df <- read.delim(host_genome_info.f, check.names = FALSE)
host_genome_info.df <- host_genome_info.df %>% mutate(Assembly.Accession = ifelse(str_detect(assembly_name, "GCA_|GCF_"), paste(str_split(assembly_name, "_", simplify = TRUE)[,1], "_", str_split(assembly_name, "_", simplify = TRUE)[,2], sep=""), assembly_name))

host_genome_ncbi_info.df <- host_genome_info.df %>% filter(str_detect(assembly_name, "GCF_|GCA_"))
host_genome_non_ncbi_info.df <- host_genome_info.df %>% filter(!str_detect(assembly_name, "GCF_|GCA_"))

##### Gather metadata of hosts of poxviruses from NCBI
host_ncbi_ds.f <- "260219_Poxviridae_and_representative_MPXV_host_NCBI_datasets.tsv"
host_ncbi_ds.df <- fread(host_ncbi_ds.f, header=T, sep="\t", quote="", check.names=T, data.table=FALSE)

poxv_ncbi_ds.virus_host_db.f <- "260219_Poxviridae_and_representative_MPXV_host_NCBI_dataset_wt_host_info.tsv"
poxv_ncbi_ds.virus_host_db.df <- poxv_ncbi_ds.df %>% left_join(virus_host_db.filtered.edited.final.df, by=c("Organism.Name"="virus.name"), relationship="many-to-many")
poxv_ncbi_ds.virus_host_db.df <- poxv_ncbi_ds.virus_host_db.df %>% mutate(host_assembly_accession_gs=host_genome_non_ncbi_info.df$assembly_name[match(gsub(" ", "_", host.name.edited.genome_size), host_genome_non_ncbi_info.df$species_name)])
poxv_ncbi_ds.virus_host_db.df <- poxv_ncbi_ds.virus_host_db.df %>% mutate(host_genome_length=host_genome_non_ncbi_info.df$genome_length[match(host_assembly_accession_gs, host_genome_non_ncbi_info.df$Assembly.Accession)])
poxv_ncbi_ds.virus_host_db.df <- poxv_ncbi_ds.virus_host_db.df %>% mutate(host_assembly_accession_gs=ifelse(is.na(host_assembly_accession_gs), host_ncbi_ds.df$Assembly.Accession[match(host.name.edited.genome_size, host_ncbi_ds.df$Organism.Name)], host_assembly_accession_gs))
poxv_ncbi_ds.virus_host_db.df <- poxv_ncbi_ds.virus_host_db.df %>% mutate(host_genome_length=ifelse(is.na(host_genome_length), host_genome_ncbi_info.df[match(host_assembly_accession_gs, host_genome_ncbi_info.df$Assembly.Accession), "genome_length"], host_genome_length))

poxv_ncbi_ds.virus_host_db.df <- poxv_ncbi_ds.virus_host_db.df %>% mutate(host_assembly_accession_if=host_ncbi_ds.df$Assembly.Accession[match(host.name.edited.inferred, host_ncbi_ds.df$Organism.Name)])
poxv_ncbi_ds.virus_host_db.df <- poxv_ncbi_ds.virus_host_db.df %>% mutate(host_gene_num=host_genome_ncbi_info.df[match(host_assembly_accession_if, host_genome_ncbi_info.df$Assembly.Accession), "protein_coding_gene_number"],
                                                                          host_genome_length_if=host_genome_ncbi_info.df[match(host_assembly_accession_if, host_genome_ncbi_info.df$Assembly.Accession), "genome_length"],
                                                                          host_total_gene_length=host_genome_ncbi_info.df[match(host_assembly_accession_if, host_genome_ncbi_info.df$Assembly.Accession), "total_gene_length"],
                                                                          host_total_cds_length=host_genome_ncbi_info.df[match(host_assembly_accession_if, host_genome_ncbi_info.df$Assembly.Accession), "total_cds_length"])
poxv_ncbi_ds.virus_host_db.df <- poxv_ncbi_ds.virus_host_db.df %>% filter(!is.na(host_assembly_accession_gs))
# write.table(poxv_ncbi_ds.virus_host_db.df, poxv_ncbi_ds.virus_host_db.f, col.names=T, row.names=F, sep="\t", quote=F) # Manually edit assembly name

##### Create an output table pooling genome statistics of each virus-host pair
virus_host_db.filtered.edited.final_wt_info.f <- "250327_VirusHostDB/virushostdb_filtered_edited_final_wt_host_info.tsv"
virus_host_db.filtered.edited.final_wt_info.df <- virus_host_db.filtered.edited.final.df %>% mutate(host_assembly_accession_gs=host_genome_non_ncbi_info.df$assembly_name[match(gsub(" ", "_", host.name.edited.genome_size), host_genome_non_ncbi_info.df$species_name)])
virus_host_db.filtered.edited.final_wt_info.df <- virus_host_db.filtered.edited.final_wt_info.df %>% mutate(host_genome_length=host_genome_non_ncbi_info.df$genome_length[match(host_assembly_accession_gs, host_genome_non_ncbi_info.df$Assembly.Accession)])
virus_host_db.filtered.edited.final_wt_info.df <- virus_host_db.filtered.edited.final_wt_info.df %>% mutate(host_assembly_accession_gs=ifelse(is.na(host_assembly_accession_gs), host_ncbi_ds.df$Assembly.Accession[match(host.name.edited.genome_size, host_ncbi_ds.df$Organism.Name)], host_assembly_accession_gs))
virus_host_db.filtered.edited.final_wt_info.df <- virus_host_db.filtered.edited.final_wt_info.df %>% mutate(host_genome_length=ifelse(is.na(host_genome_length), host_genome_ncbi_info.df[match(host_assembly_accession_gs, host_genome_ncbi_info.df$Assembly.Accession), "genome_length"], host_genome_length))

virus_host_db.filtered.edited.final_wt_info.df <- virus_host_db.filtered.edited.final_wt_info.df %>% mutate(host_assembly_accession_if=host_ncbi_ds.df$Assembly.Accession[match(host.name.edited.inferred, host_ncbi_ds.df$Organism.Name)])
virus_host_db.filtered.edited.final_wt_info.df <- virus_host_db.filtered.edited.final_wt_info.df %>% mutate(host_gene_num=host_genome_ncbi_info.df[match(host_assembly_accession_if, host_genome_ncbi_info.df$Assembly.Accession), "protein_coding_gene_number"],
                                                                          host_genome_length_if=host_genome_ncbi_info.df[match(host_assembly_accession_if, host_genome_ncbi_info.df$Assembly.Accession), "genome_length"],
                                                                          host_total_gene_length=host_genome_ncbi_info.df[match(host_assembly_accession_if, host_genome_ncbi_info.df$Assembly.Accession), "total_gene_length"],
                                                                          host_total_cds_length=host_genome_ncbi_info.df[match(host_assembly_accession_if, host_genome_ncbi_info.df$Assembly.Accession), "total_cds_length"])
# write.table(virus_host_db.filtered.edited.final_wt_info.df, virus_host_db.filtered.edited.final_wt_info.f, col.names=T, row.names=F, sep="\t", quote=F) # Manually edit assembly name

virus_host_db.filtered.edited.final_wt_info.df <- virus_host_db.filtered.edited.final_wt_info.df %>% mutate(clade=poxv_ncbi_ds.df$clade[match(gsub(" ", "_", virus.name), poxv_ncbi_ds.df$tip_label_dummy)])

##### Compare genome size of hosts and animals closely related to hosts
poxv_host_genome_size_cp.df <- virus_host_db.filtered.edited.final_wt_info.df %>% select(host.name.edited, host.name.edited.genome_size, host.name.edited.inferred, host_genome_length, host_genome_length_if, clade)
poxv_host_genome_size_cp.df <- poxv_host_genome_size_cp.df %>% filter(!is.na(host.name.edited) & !is.na(host.name.edited.genome_size) & !is.na(host.name.edited.inferred) & host.name.edited == host.name.edited.genome_size)
poxv_host_genome_size_cp.df <- poxv_host_genome_size_cp.df %>% mutate(g1 = str_split(host.name.edited.genome_size, " ", simplify = TRUE)[,1],
                                                                      g2 = str_split(host.name.edited.inferred,    " ", simplify = TRUE)[,1],
                                                                      cond1 = host.name.edited.genome_size != host.name.edited.inferred & g1 == g2,
                                                                      cond2 = host.name.edited.genome_size != host.name.edited.inferred & g1 != g2,
                                                                      condition = case_when(cond1 ~ "same_genus", cond2 ~ "different_genus", TRUE  ~ NA_character_)) %>% filter(!is.na(condition))

poxv_host_genome_size_cp.p <- ggplot(poxv_host_genome_size_cp.df, aes(x=host_genome_length, y=host_genome_length_if, color=condition))
poxv_host_genome_size_cp.p <- poxv_host_genome_size_cp.p + geom_point()
poxv_host_genome_size_cp.p <- poxv_host_genome_size_cp.p + geom_text(aes(label=substr(clade, 1, 3)), hjust=-0.13, vjust=0, size=7)
poxv_host_genome_size_cp.p <- poxv_host_genome_size_cp.p + geom_abline(slope = 1, intercept = 0, color = "black")
poxv_host_genome_size_cp.p <- poxv_host_genome_size_cp.p + labs(title="Genome size of host and animals closely related to the hostsd", x="Genome size of hosts", y="Genome size of animals closely related to the hosts")
poxv_host_genome_size_cp.p <- poxv_host_genome_size_cp.p + theme_classic()
poxv_host_genome_size_cp.p <- poxv_host_genome_size_cp.p + theme(text=element_text(size = 24), plot.title=element_text(size = 32), axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
poxv_host_genome_size_cp.p <- poxv_host_genome_size_cp.p + coord_fixed()
poxv_host_genome_size_cp.p

##### Extract host information from Virus-Host DB database
virushostdb_poxv.f <- "250917_virushostdb_poxviridae.tsv"
virushostdb_poxv.df <- read.delim(virushostdb_poxv.f, check.names = FALSE)
colnames(virushostdb_poxv.df) <- gsub(" ", "_", colnames(virushostdb_poxv.df))

virushostdb_poxv_long.df <- virushostdb_poxv.df %>% select(virus_tax_id, virus_name, host_tax_id, host_name, host_lineage) %>% separate_rows(host_lineage, sep = "; ") %>% mutate(host_lineage = trimws(host_lineage))
virushostdb_poxv_long.df <- virushostdb_poxv_long.df %>% mutate(virus_name=gsub(" ", "_", virus_name)) %>% mutate(virus_name=gsub("'", "", virus_name))
virushostdb_poxv_long.df <- virushostdb_poxv_long.df %>% mutate(virus_name=ifelse(virus_name=="Parapoxvirus_red_deer/HL953", "Parapoxvirus_red_deer_HL953", 
                                                                           ifelse(virus_name=="Monkeypox_virus_Zaire-96-I-16", "Monkeypox_virus_Zaire-96-I-16_Ia",
                                                                           ifelse(virus_name=="Monkeypox_virus", "Monkeypox_virus_MPXV-M5312_HM12_Rivers_IIb",
                                                                           ifelse(str_detect(virus_name, "Vaccinia_virus"), "Vaccinia_virus", virus_name)))))
virushostdb_poxv_long.df <- rbind(virushostdb_poxv_long.df, virushostdb_poxv_long.df %>% filter(virus_name=="Monkeypox_virus_Zaire-96-I-16_Ia") %>% mutate(virus_name="Monkeypox_virus_RDC-NKV-GOM-MPOX-010_Ib"))
virushostdb_poxv_long.df <- rbind(virushostdb_poxv_long.df, virushostdb_poxv_long.df %>% filter(virus_name=="Monkeypox_virus_Zaire-96-I-16_Ia") %>% mutate(virus_name="Monkeypox_virus_Sierra_Leone_IIa"))



############################## Phylogenetic analysis ##############################
########## Retrieve tree information
poxv_t.f <- "241008_Poxviridae_and_representative_MPXV_tree_reconstruction/supermatrix_partition_protein.nex.treefile"
poxv_t <- read.tree(poxv_t.f) # The tree is already rooted.
poxv_t.rt <- phytools::reroot(poxv_t, node=75, position=0.5*poxv_t$edge.length[which(poxv_t$edge[,2]==75)]) # Root at midpoint of specified and earlier nodes

poxv_t.rt.nolabs <- poxv_t.rt
poxv_t.rt.nolabs$node.label <- NULL
poxv_t.rt.nolabs.f <- "241008_Poxviridae_and_representative_MPXV_tree_reconstruction/supermatrix_partition_protein_rooted.nex.treefile"
# write.tree(poxv_t.rt.nolabs, file=poxv_t.rt.nolabs.f, digits = 10)

poxv_t.rt.withlabs <- poxv_t.rt
poxv_t.rt.withlabs$node.label <- c("Root", ((ggtree(poxv_t.rt.withlabs)$data %>% filter(isTip!=TRUE))$node)[2:length((ggtree(poxv_t.rt.withlabs)$data %>% filter(isTip!=TRUE))$label)])

poxv_t.rt.withlabs.f <- "241008_Poxviridae_and_representative_MPXV_tree_reconstruction/supermatrix_partition_protein_rooted_pastml.nex.treefile"
# write.tree(poxv_t.rt.withlabs, file=poxv_t.rt.withlabs.f, digits = 10)

########## Visualize a species tree and generate tree files from analyses by OrthoFinder and PastML
poxv_t.rt.rn <- poxv_t.rt
poxv_t.rt.rn$tip.label <- ifelse(poxv_t.rt.rn$tip.label=="GCF_000857045.1_Monkeypox_virus", "GCF_000857045.1_Monkeypox_virus_Zaire-96-I-16_Ia",
                          ifelse(poxv_t.rt.rn$tip.label=="GCA_039269415.1_Monkeypox_virus", "GCA_039269415.1_Monkeypox_virus_RDC-NKV-GOM-MPOX-010_Ib",
                          ifelse(poxv_t.rt.rn$tip.label=="GCA_006465585.1_Monkeypox_virus", "GCA_006465585.1_Monkeypox_virus_Sierra_Leone_IIa",
                          ifelse(poxv_t.rt.rn$tip.label=="GCF_014621545.1_Monkeypox_virus", "GCF_014621545.1_Monkeypox_virus_MPXV-M5312_HM12_Rivers_IIb", poxv_t.rt.rn$tip.label))))
poxv_t.rt.rn$tip.label <- unname(sapply(poxv_t.rt.rn$tip.label, function(x) {
  parts <- str_split(x, "_", simplify = TRUE)
  paste(parts[3:length(parts)], collapse = "_")  # Join everything after the first two parts
}))

poxv_t.rt.info.df <- ggtree(poxv_t.rt.rn)$data
poxv_t.rt.info.df <- poxv_t.rt.info.df %>% mutate(clade=poxv_ncbi_ds.df[match(label, poxv_ncbi_ds.df$tip_label_dummy), "clade"])

##### Calculate proportion of nodes whose SH-aLRT > 80 and UFBoot2 > 95
nrow(poxv_t.rt.info.df %>% filter(!isTip, parent!=node) %>% mutate(label=as.numeric(str_split(label, "/", simplify=T)[,1])) %>% filter(label >= 80))/nrow(poxv_t.rt.info.df %>% filter(!isTip, parent!=node) %>% mutate(label=as.numeric(str_split(label, "/", simplify=T)[,1])))
nrow(poxv_t.rt.info.df %>% filter(!isTip, parent!=node) %>% mutate(label=as.numeric(str_split(label, "/", simplify=T)[,2])) %>% filter(label >= 95))/nrow(poxv_t.rt.info.df %>% filter(!isTip, parent!=node) %>% mutate(label=as.numeric(str_split(label, "/", simplify=T)[,2])))

########## ***** Figure 1A: Poxvirus tree with gene number and host information (Tree and Gene number) *****
poxv_t.p <- ggtree(poxv_t.rt.rn, linewidth=1)
poxv_t.p <- poxv_t.p %<+% poxv_t.rt.info.df
poxv_t.p <- poxv_t.p + geom_tiplab(size=5)
# poxv_t.p <- poxv_t.p + geom_text2(aes(subset=!isTip, label=label), size=5, hjust=1.3)
poxv_t.p <- poxv_t.p + geom_text2(aes(subset=!isTip & as.numeric(str_split(label, "/", simplify=T)[,2]) < 95, label=label), size=5, hjust=1.3)
# poxv_t.p <- poxv_t.p + geom_text2(aes(subset=!isTip & as.numeric(str_split(label, "/", simplify=T)[,2]) < 95, label=str_split(label, "/", simplify=T)[,2]), size=5, hjust=1.3)
# poxv_t.p <- poxv_t.p + geom_text2(aes(subset=!isTip, label=node), size=5, hjust=1.3)
poxv_t.p <- poxv_t.p + geom_point2(shape=16, size=4, aes(subset=isTip, color=clade))
poxv_t.p <- poxv_t.p + scale_color_manual("clade", name="Clade", values=poxv_color.clade.v, breaks=names(poxv_color.clade.v), na.value = NA)
poxv_t.p <- poxv_t.p + theme_tree()
poxv_t.p <- poxv_t.p + theme(text = element_text(size = 20))
poxv_t.p <- poxv_t.p + geom_treescale(width=0.1)
poxv_t.p <- poxv_t.p + xlim(0,6)
poxv_t.p <- poxv_t.p + geom_cladelab(node=60, label="Chordopoxvirus\n(Vertebrate poxvirus)", offset=.45, offset.text=.001, fontsize=5)
poxv_t.p <- poxv_t.p + geom_cladelab(node=68, label="Orthopoxvirus", offset=.025, offset.text=.001, fontsize=5)
poxv_t.p <- poxv_t.p + geom_cladelab(node=83, label="Centapoxvirus", offset=.025, offset.text=.001, fontsize=5)
poxv_t.p <- poxv_t.p + geom_cladelab(node=90, label="Capripoxvirus", offset=.025, offset.text=.001, fontsize=5)
poxv_t.p <- poxv_t.p + geom_cladelab(node=93, label="Leporipoxvirus", offset=.025, offset.text=.001, fontsize=5)
poxv_t.p <- poxv_t.p + geom_cladelab(node=94, label="Oryzopoxvirus", offset=.025, offset.text=.001, fontsize=5)
poxv_t.p <- poxv_t.p + geom_cladelab(node=95, label="Yatapoxvirus", offset=.025, offset.text=.001, fontsize=5)
poxv_t.p <- poxv_t.p + geom_cladelab(node=98, label="Parapoxvirus", offset=.025, offset.text=.001, fontsize=5)
poxv_t.p <- poxv_t.p + geom_cladelab(node=103, label="Avipoxvirus", offset=.025, offset.text=.001, fontsize=5)
poxv_t.p <- poxv_t.p + geom_cladelab(node=108, label="Macropopoxvirus", offset=.025, offset.text=.001, fontsize=5)
poxv_t.p <- poxv_t.p + geom_cladelab(node=109, label="Entomopoxvirus\n(Insect poxvirus)", offset=.3, offset.text=.001, fontsize=5)
poxv_t.p <- poxv_t.p + geom_cladelab(node=112, label="Betaentomopoxvirus", offset=.025, offset.text=.001, fontsize=5)
poxv_t.p <- poxv_t.p + scale_fill_manual("clade", name="Clade", values=poxv_color.clade.v, breaks=names(poxv_color.clade.v), na.value = NA)
poxv_t.p <- poxv_t.p + geom_fruit(data = poxv_genome_stats.df %>% rename(label = species_name) %>% rename(clade_bar = clade),
                                                                         geom = geom_col, mapping = aes(x = gene_number, y = label, fill = clade_bar),
                                                                         orientation = "y", offset = 0.3, width = 0.6, pwidth = 0.2, grid.params=list(),
                                                                         axis.params = list(axis = "x", title = "Gene number", title.size = 6, text.size = 4, text.angle = 90, hjust = 1, limits = c(0, 400)))
poxv_t.p <- poxv_t.p + guides(fill=guide_legend(ncol=1), color="none")
poxv_t.p <- poxv_t.p + labs(title="Poxvirus gene number", fill = "Clade", x = "Gene number")
poxv_t.p <- poxv_t.p + coord_cartesian(clip="off")
poxv_t.p

# pdf("Figures/poxvirus_phylogeny_wt_gene_number.pdf", width=16, height=16)
# print(poxv_t.p)
# dev.off()

########## ***** Figure 1A: Poxvirus tree with gene number and host information (Host information) *****
virushostdb_poxv_long_tree.df <- virushostdb_poxv_long.df %>% filter(host_lineage %in% c("Mammalia","Insecta","Actinopteri","Aves","Crocodylia"))

virus_host_mat <- virushostdb_poxv_long_tree.df %>% distinct(virus_name, host_lineage) %>% mutate(value = "1") %>% pivot_wider(names_from = host_lineage, values_from = value, values_fill = list(value = "0")) %>% column_to_rownames("virus_name")
virus_host_mat <- virus_host_mat %>% select(Insecta, Actinopteri, Aves, Crocodylia, Mammalia)
virus_host_mat <- virus_host_mat[poxv_t.p$data %>% filter(isTip) %>% arrange(desc(y)) %>% pull(label),]

virus_host_mat_tree <- virus_host_mat %>% filter(row.names(.) %in% poxv_t.rt.rn$tip.label)
virus_host_mat_tree_ht <- Heatmap(virus_host_mat_tree, cluster_rows = FALSE, cluster_columns = FALSE, row_names_side = "left", column_names_side = "top",
              row_gap = unit(5, "mm"), column_gap = unit(5, "mm"), col = c("white", "red"), heatmap_legend_param = list(title = "Host taxon"))   
virus_host_mat_tree_ht

# pdf("Figures/poxvirus_phylogeny_host_heatmap.pdf", width=6.5415, height=24)
# draw(virus_host_mat_tree_ht)
# dev.off()


########## PIC analysis for poxvirus gene number and genome size
poxv_genome_stats.df.host_ord <- tibble(species_name = poxv_t.rt.rn$tip.label) %>% left_join(poxv_genome_stats.df, by = "species_name")

gene_no_pic <- pic(poxv_genome_stats.df.host_ord$gene_number, poxv_t.rt.rn)
genome_l_pic <- pic(poxv_genome_stats.df.host_ord$genome_length, poxv_t.rt.rn)

poxv_gene_no_vs_genome_size_cont_picModel <- lm(gene_no_pic ~ genome_l_pic - 1)
summary(poxv_gene_no_vs_genome_size_cont_picModel)
anova(poxv_gene_no_vs_genome_size_cont_picModel)

poxv_gene_no_vs_genome_size_cont_picModel.df <- data.frame(gene_no_contrast = gene_no_pic, genome_l_contrast = genome_l_pic)

########## ***** Supplementary Figure 1B: Correlation between PICs of poxvirus gene number and genome size *****
poxv_gene_no_vs_genome_size_cont.p <- ggplot(poxv_gene_no_vs_genome_size_cont_picModel.df, aes(x=genome_l_contrast, y=gene_no_contrast))
poxv_gene_no_vs_genome_size_cont.p <- poxv_gene_no_vs_genome_size_cont.p + geom_abline(slope = coef(poxv_gene_no_vs_genome_size_cont_picModel), intercept = 0, color = "black")
poxv_gene_no_vs_genome_size_cont.p <- poxv_gene_no_vs_genome_size_cont.p + geom_point()
poxv_gene_no_vs_genome_size_cont.p <- poxv_gene_no_vs_genome_size_cont.p + labs(title="Relationship between gene number and genome size", 
                                                                                subtitle=paste("r = ", round(sign(coef(poxv_gene_no_vs_genome_size_cont_picModel)) * sqrt(summary(poxv_gene_no_vs_genome_size_cont_picModel)$r.squared), 2), " ; P-value = ",
                                                                                               sprintf("%.2e", summary(poxv_gene_no_vs_genome_size_cont_picModel)$coefficients[1, 4]), sep=""), x="Contrast in virus genome size", y="Contrast in virus gene number")
poxv_gene_no_vs_genome_size_cont.p <- poxv_gene_no_vs_genome_size_cont.p + guides(color="none")
poxv_gene_no_vs_genome_size_cont.p <- poxv_gene_no_vs_genome_size_cont.p + theme_classic()
poxv_gene_no_vs_genome_size_cont.p <- poxv_gene_no_vs_genome_size_cont.p + theme(text=element_text(size = 24), plot.title=element_text(size = 32), axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
poxv_gene_no_vs_genome_size_cont.p <- poxv_gene_no_vs_genome_size_cont.p + scale_x_continuous(breaks=seq(-300000, 3000000, by=100000), limits=c(-300000,300000))
poxv_gene_no_vs_genome_size_cont.p <- poxv_gene_no_vs_genome_size_cont.p + scale_y_continuous(breaks=seq(-400, 500, by=100), limits=c(-400,500))
poxv_gene_no_vs_genome_size_cont.p

# pdf("Figures/gene_number_vs_genome_size_contrast.pdf", width=16, height=12)
# print(poxv_gene_no_vs_genome_size_cont.p)
# dev.off()



############################## Ortholog group gain and loss analysis ##############################
poxv_noc.v <- c(73,75,3,2,74,1,50)
poxv_noc.df <- data.frame(node=poxv_noc.v, node_name=c("MPXV","I","Ia","Ib","II","IIa","IIb"))
# 73: MPXV
# 75: Clade I
# 3: Clade Ia
# 2: Clade Ib
# 74: Clade II
# 1: Clade IIa
# 50: Clade IIb
# Node 67: Orthopoxvirus and Centapoxvirus
# Node 73: Monkeypoxvirus
# Node 94: Oryzopoxvirus
# Node 95: Yatapoxvirus
# Node 98: Parapoxvirus
# Node 103: Avipoxvirus (no gene name)
# Node 112: Betaentomopoxvirus

########## Identify unassigned proteins (single-gene ortholog group)
hog_res.f <- "251208_gain_and_loss_analysis/N0_pooled_additional_assigned.tsv"
hog_res.df <- read.delim(hog_res.f, check.names = FALSE)

hog_res_ct.df <- hog_res.df
hog_res_ct.df[c("OG","Gene Tree Parent Clade")] <- NULL

hog_res_ct.df <- hog_res_ct.df %>% select(-1) %>% mutate(across(everything(), ~ if_else(is.na(.x) | str_trim(.x) == "", 0L, lengths(str_split(.x, fixed(", "))))))
hog_res_ct.df <- bind_cols(hog_res.df[1], hog_res_ct.df)

hog_res_ct_binary.df <- hog_res_ct.df %>% mutate(across(-HOG, ~ as.integer(.x > 0))) %>% mutate(total = rowSums(across(-1), na.rm = TRUE))

poxv_single_sp_hog.v <- gsub("\\.", "", (hog_res_ct_binary.df %>% filter(total==1))$HOG)

########## Plot results and count ortholog group gain and loss events in the poxvirus lineage
poxv_vacv_gene_name.f <- "251208_gain_and_loss_analysis/N0_pooled_additional_assigned_named_summarized.tsv"
poxv_vacv_gene_name.df <- read.delim(poxv_vacv_gene_name.f, check.names = FALSE, header = T)
poxv_vacv_gene_name.df <- poxv_vacv_gene_name.df %>% mutate(og_id_dummy = gsub("\\.", "", og_id))
poxv_vacv_gene_name.v <- poxv_vacv_gene_name.df$Vaccinia_virus

poxv_state_tr_pooled.df <- data.frame()
poxv_state_tr.p_list <- list()

##### Plot ortholog group gain and loss events in the poxvirus lineage identified by MPPA, MAP, and DOWNPASS methods
for (a in c("MPPA","MAP","DOWNPASS")) {
  poxv_pastml_rs.f <- paste("251208_gain_and_loss_analysis/251208_cs_matrix_pooled_additional_pastml", a, "output.tsv", sep="_")
  poxv_pastml_rs <- fread(poxv_pastml_rs.f, header=T, sep="\t", quote="", check.names=T, data.table=FALSE)
  
  poxv_pastml_rs_state.df <- poxv_pastml_rs
  poxv_pastml_rs_state.df[which(poxv_pastml_rs_state.df$node=="Root"),"node"] <- find_root(poxv_t.rt.rn)
  
  poxv_pastml_rs_state.df <- poxv_pastml_rs_state.df %>% mutate(node=ifelse(node=="GCF_000857045.1_Monkeypox_virus", "GCF_000857045.1_Monkeypox_virus_Zaire-96-I-16_Ia",
                                                                     ifelse(node=="GCA_039269415.1_Monkeypox_virus", "GCA_039269415.1_Monkeypox_virus_RDC-NKV-GOM-MPOX-010_Ib",
                                                                     ifelse(node=="GCA_006465585.1_Monkeypox_virus", "GCA_006465585.1_Monkeypox_virus_Sierra_Leone_IIa",
                                                                     ifelse(node=="GCF_014621545.1_Monkeypox_virus", "GCF_014621545.1_Monkeypox_virus_MPXV-M5312_HM12_Rivers_IIb",node)))))
  poxv_pastml_rs_state.df <- poxv_pastml_rs_state.df %>% mutate(node=ifelse(str_detect(node, "GCF_|GCA_"), unname(sapply(node, function(x) {
    parts <- str_split(x, "_", simplify = TRUE)
    paste(parts[3:length(parts)], collapse = "_")  # Join everything after the first two parts
    })), node))
  
  poxv_state_tr_dir <- paste("251208_gain_and_loss_analysis/state_transitions_", a, "/", sep="")
  if (!file.exists(poxv_state_tr_dir)) {dir.create(poxv_state_tr_dir, recursive = TRUE)}
  
  poxv_state_tr.df <- as.data.frame(ggtree(poxv_t.rt.rn)$data)
  poxv_state_tr.df <- poxv_state_tr.df %>% mutate(node_name = ifelse(isTip==FALSE, node, label))
  
  poxv_og_name.v <- colnames(poxv_pastml_rs_state.df[,2:ncol(poxv_pastml_rs_state.df)])
  
  for (i in 1:length(poxv_og_name.v)) {
    poxv_state_tr_og_dir <- paste("251208_gain_and_loss_analysis/state_transitions_", a, "/", poxv_og_name.v[i], "/", sep="")
    if (!file.exists(poxv_state_tr_og_dir)) {dir.create(poxv_state_tr_og_dir, recursive = TRUE)}
    
    poxv_pastml_rs_state_og.df <- poxv_pastml_rs_state.df[,c("node",poxv_og_name.v[i])]
    
    poxv_state_tr_og.df <- poxv_state_tr.df
    poxv_state_tr_og.df$state_parent <- as.character(poxv_pastml_rs_state_og.df[match(poxv_state_tr_og.df$parent, poxv_pastml_rs_state_og.df$node), poxv_og_name.v[i]])
    
    if (poxv_og_name.v[i] %in% poxv_single_sp_hog.v) {poxv_state_tr_og.df$state_parent <- "0"}
    
    poxv_state_tr_og.df$state_desc <- as.character(poxv_pastml_rs_state_og.df[match(poxv_state_tr_og.df$node_name, poxv_pastml_rs_state_og.df$node), poxv_og_name.v[i]])
    
    poxv_state_tr_og.df$og_name <- poxv_og_name.v[i]
    poxv_state_tr_og.df$transition <- paste(poxv_state_tr_og.df$state_parent, poxv_state_tr_og.df$state_desc, sep="->")
    
    ### Save transition branch info
    poxv_state_tr_og.df.f <- paste(poxv_state_tr_og_dir, poxv_og_name.v[i], "_transition_info.tsv", sep="")
    write.table(poxv_state_tr_og.df, poxv_state_tr_og.df.f, col.names=T, row.names=F, sep="\t", quote=F)
    
    poxv_state_tr_og.df <- poxv_state_tr_og.df %>% mutate(noc_name=ifelse(node %in% poxv_noc.v, poxv_noc.df[match(node, poxv_noc.df$node), "node_name"], ""),
                                                                state_desc=ifelse(state_desc=="0", "Absent",
                                                                           ifelse(state_desc=="0|1", "Undeterminable",
                                                                           ifelse(state_desc=="1", "Present", NA))))
    
    poxv_state_tr_og.p <- ggtree(poxv_t.rt.rn, linewidth=0.25)
    poxv_state_tr_og.p <- poxv_state_tr_og.p %<+% poxv_state_tr_og.df
    poxv_state_tr_og.p <- poxv_state_tr_og.p + theme_tree()
    poxv_state_tr_og.p <- poxv_state_tr_og.p + geom_tiplab(size=2, align=FALSE, linesize=.5)
    poxv_state_tr_og.p <- poxv_state_tr_og.p + geom_point2(shape=16, size=3, aes(color=state_desc))
    poxv_state_tr_og.p <- poxv_state_tr_og.p + scale_color_manual("state_desc", name="State", values=c("lightsteelblue1","grey","red"), breaks=c("Absent","Undeterminable","Present"), na.value = NA)
    poxv_state_tr_og.p <- poxv_state_tr_og.p + labs(title=paste(poxv_og_name.v[i], ": ", poxv_vacv_gene_name.v[i], sep=""), subtitle=paste(a, ": Transition no. = ", nrow(poxv_state_tr_og.df %>% filter(transition %in% c("0->1", "1->0"))), sep=""))
    poxv_state_tr_og.p <- poxv_state_tr_og.p + theme(text = element_text(size = 20), plot.title = element_text(size = 20), plot.subtitle = element_text(size = 16))
    poxv_state_tr_og.p <- poxv_state_tr_og.p + xlim(0,4.5)
    poxv_state_tr_og.p <- poxv_state_tr_og.p + geom_treescale(width=0.1)
    
    # ### Save the tree plot
    # pdf(paste(poxv_state_tr_og_dir, poxv_og_name.v[i], "_transition_plot.pdf", sep=""), width=12, height=12)
    # poxv_state_tr_og_to_save.p <- poxv_state_tr_og.p + geom_text2(aes(subset=(transition %in% c("0->1", "1->0")), label=transition), size=3, hjust=1.3)
    # print(poxv_state_tr_og.p)
    # dev.off()
    
    ### Pool data from summarizing
    poxv_state_tr_og.df <- poxv_state_tr_og.df %>% mutate(method = a)
    poxv_state_tr_pooled.df <- rbind(poxv_state_tr_pooled.df, poxv_state_tr_og.df)
    
    ### Record all plots
    poxv_state_tr_og.p <- poxv_state_tr_og.p + geom_text2(aes(subset=((transition %in% c("0->1", "1->0")) & (node %in% poxv_noc.v)), 
                                                              label=ifelse(!is.na(noc_name) & !is.na(transition), paste(noc_name, transition, sep=":"), "")), size=3, hjust=1.3)
    poxv_state_tr.p_list[[poxv_og_name.v[i]]][[paste(poxv_og_name.v[i] ,"_" , a, sep="")]] <- poxv_state_tr_og.p
  }
}


##### Plot counts of ortholog group gain and loss events identified by MPPA, MAP, and DOWNPASS methods
poxv_no_gain_loss.p_list <- list()

for (a in c("MPPA","MAP","DOWNPASS")) {
  poxv_state_tr_pooled_ct.df <- poxv_state_tr_pooled.df %>% filter(parent!=node, transition %in% c("0->1", "1->0"), method==a) %>% group_by(parent, node, transition) %>% summarise(transition_count = n()) %>% ungroup()

  ########## ***** Supplementary Figure 2A: Number of ortholog group gain events identified by MPPA, MAP, and DOWNPASS methods *****
  poxv_state_tr_pooled_gain_ct.info.df <- poxv_t.rt.info.df %>% left_join(poxv_state_tr_pooled_ct.df %>% filter(transition == "0->1") %>% select(parent, node, transition_count), by = c("parent", "node"))
  poxv_state_tr_pooled_gain_ct.info.df <- poxv_state_tr_pooled_gain_ct.info.df %>% mutate(transition_count=ifelse(is.na(transition_count), 0, transition_count))
  
  poxv_state_tr_pooled_gain_ct.p <- ggtree(poxv_t.rt.rn, linewidth=2, aes(color=transition_count))
  poxv_state_tr_pooled_gain_ct.p <- poxv_state_tr_pooled_gain_ct.p %<+% poxv_state_tr_pooled_gain_ct.info.df
  poxv_state_tr_pooled_gain_ct.p <- poxv_state_tr_pooled_gain_ct.p + theme_tree()
  # poxv_state_tr_pooled_gain_ct.p <- poxv_state_tr_pooled_gain_ct.p + geom_text2(aes(label=transition_count), color="black", size=5, hjust=1.3)
  poxv_state_tr_pooled_gain_ct.p <- poxv_state_tr_pooled_gain_ct.p + scale_color_gradient("transition_count", trans="pseudo_log", name="Number of ortholog group gain events", low="black", high="yellow", breaks=c(0,3,5,10,30,50,100,max(poxv_state_tr_pooled_gain_ct.info.df$transition_count)))
  poxv_state_tr_pooled_gain_ct.p <- poxv_state_tr_pooled_gain_ct.p + labs(title=a)
  poxv_state_tr_pooled_gain_ct.p <- poxv_state_tr_pooled_gain_ct.p + theme(text = element_text(size = 24), plot.title = element_text(size = 30), plot.subtitle = element_text(size = 16))
  poxv_state_tr_pooled_gain_ct.p <- poxv_state_tr_pooled_gain_ct.p + xlim(0,4.5)

  ########## ***** Supplementary Figure 2B: Number of ortholog group loss events identified by MPPA, MAP, and DOWNPASS methods *****
  poxv_state_tr_pooled_loss_ct.info.df <- poxv_t.rt.info.df %>% left_join(poxv_state_tr_pooled_ct.df %>% filter(transition == "1->0") %>% select(parent, node, transition_count), by = c("parent", "node"))
  poxv_state_tr_pooled_loss_ct.info.df <- poxv_state_tr_pooled_loss_ct.info.df %>% mutate(transition_count=ifelse(is.na(transition_count), 0, transition_count))
  
  poxv_state_tr_pooled_loss_ct.p <- ggtree(poxv_t.rt.rn, linewidth=2, aes(color=transition_count))
  poxv_state_tr_pooled_loss_ct.p <- poxv_state_tr_pooled_loss_ct.p %<+% poxv_state_tr_pooled_loss_ct.info.df
  poxv_state_tr_pooled_loss_ct.p <- poxv_state_tr_pooled_loss_ct.p + theme_tree()
  # poxv_state_tr_pooled_loss_ct.p <- poxv_state_tr_pooled_loss_ct.p + geom_text2(aes(label=transition_count), color="black", size=5, hjust=1.3)
  poxv_state_tr_pooled_loss_ct.p <- poxv_state_tr_pooled_loss_ct.p + scale_color_gradient("transition_count", trans="pseudo_log", name="Number of ortholog group loss events", low="black", high="deepskyblue", breaks=c(0,3,5,10,30,50,100,max(poxv_state_tr_pooled_loss_ct.info.df$transition_count)))
  poxv_state_tr_pooled_loss_ct.p <- poxv_state_tr_pooled_loss_ct.p + labs(title=a)
  poxv_state_tr_pooled_loss_ct.p <- poxv_state_tr_pooled_loss_ct.p + theme(text = element_text(size = 24), plot.title = element_text(size = 30), plot.subtitle = element_text(size = 16))
  poxv_state_tr_pooled_loss_ct.p <- poxv_state_tr_pooled_loss_ct.p + xlim(0,4.5)

  poxv_no_gain_loss.p_list[[a]] <- poxv_state_tr_pooled_gain_ct.p + poxv_state_tr_pooled_loss_ct.p + plot_layout(ncol=1)
}

# pdf("Figures/gain_and_loss_all_methods.pdf", width=24, height=18)
# wrap_plots(poxv_no_gain_loss.p_list) + plot_layout(nrow=1)
# dev.off()


##### Plot counts of ortholog group gain and loss events correspondingly identified by all three methods
poxv_state_tr_all_m.df <- poxv_state_tr_pooled.df %>% filter(parent!=node) %>% group_by(parent, node, transition, og_name) %>% summarise(transition_count = n()) %>% ungroup()
poxv_state_tr_all_m.df <- poxv_state_tr_all_m.df %>% mutate(transition_count = transition_count/3) %>% filter(transition_count == 1) # Pick only transition events estimated correspondingly by all methods

table((poxv_state_tr_all_m.df %>% filter(transition %in% c("0->1", "1->0")))$transition) # Numbers of ortholog group gain and loss events

poxv_state_tr_all_m_ct.df <- poxv_state_tr_all_m.df %>% group_by(parent, node, transition) %>% summarise(transition_count = n()) %>% ungroup()

poxv_state_tr_all_m_total_ct.df <- poxv_state_tr_all_m_ct.df %>% filter(transition %in% c("0->1", "1->0")) %>% mutate(transition_count=ifelse(transition=="1->0", transition_count*-1, transition_count))
poxv_state_tr_all_m_total_ct.df <- poxv_state_tr_all_m_total_ct.df %>% group_by(parent, node) %>% summarise(transition_total = sum(transition_count)) %>% ungroup()

########## ***** Figure 1E, left: Number of ortholog group gain events correspondingly identified by all methods *****
poxv_state_tr_all_m_total_gain_ct.info.df <- poxv_t.rt.info.df %>% left_join(poxv_state_tr_all_m_ct.df %>% filter(transition == "0->1") %>% select(parent, node, transition_count), by = c("parent", "node"))
poxv_state_tr_all_m_total_gain_ct.info.df <- poxv_state_tr_all_m_total_gain_ct.info.df %>% mutate(transition_count=ifelse(is.na(transition_count), 0, transition_count))

poxv_state_tr_all_m_gain_ct.p <- ggtree(poxv_t.rt.rn, linewidth=2, aes(color=transition_count))
poxv_state_tr_all_m_gain_ct.p <- poxv_state_tr_all_m_gain_ct.p %<+% poxv_state_tr_all_m_total_gain_ct.info.df
poxv_state_tr_all_m_gain_ct.p <- poxv_state_tr_all_m_gain_ct.p + theme_tree()
# poxv_state_tr_all_m_gain_ct.p <- poxv_state_tr_all_m_gain_ct.p + geom_text2(aes(label=round((transition_count), 1)), color="black", size=5, hjust=1.3)
poxv_state_tr_all_m_gain_ct.p <- poxv_state_tr_all_m_gain_ct.p + scale_color_gradient("transition_count", trans="pseudo_log", name="Number of ortholog group gain events", low="black", high="yellow", breaks=c(0,3,5,10,30,50,100,max(poxv_state_tr_all_m_total_gain_ct.info.df$transition_count)))
poxv_state_tr_all_m_gain_ct.p <- poxv_state_tr_all_m_gain_ct.p + labs(title="All methods", subtitle="Ortholog group gain")
poxv_state_tr_all_m_gain_ct.p <- poxv_state_tr_all_m_gain_ct.p + theme(text = element_text(size = 20), plot.title = element_text(size = 30), plot.subtitle = element_text(size = 24))
poxv_state_tr_all_m_gain_ct.p <- poxv_state_tr_all_m_gain_ct.p + xlim(0,4.5)
poxv_state_tr_all_m_gain_ct.p

########## ***** Figure 1E, right: Number of ortholog group loss events correspondingly identified by all methods *****
poxv_state_tr_all_m_total_loss_ct.info.df <- poxv_t.rt.info.df %>% left_join(poxv_state_tr_all_m_ct.df %>% filter(transition == "1->0") %>% select(parent, node, transition_count), by = c("parent", "node"))
poxv_state_tr_all_m_total_loss_ct.info.df <- poxv_state_tr_all_m_total_loss_ct.info.df %>% mutate(transition_count=ifelse(is.na(transition_count), 0, transition_count))

poxv_state_tr_all_m_loss_ct.p <- ggtree(poxv_t.rt.rn, linewidth=2, aes(color=transition_count+1))
poxv_state_tr_all_m_loss_ct.p <- poxv_state_tr_all_m_loss_ct.p %<+% poxv_state_tr_all_m_total_loss_ct.info.df
poxv_state_tr_all_m_loss_ct.p <- poxv_state_tr_all_m_loss_ct.p + theme_tree()
# poxv_state_tr_all_m_loss_ct.p <- poxv_state_tr_all_m_loss_ct.p + geom_text2(aes(label=transition_count), color="black", size=5, hjust=1.3)
poxv_state_tr_all_m_loss_ct.p <- poxv_state_tr_all_m_loss_ct.p + scale_color_gradient("transition_count", trans="pseudo_log", name="Number of ortholog group loss events", low="black", high="deepskyblue", breaks=c(0,3,5,10,30,50,100,max(poxv_state_tr_all_m_total_loss_ct.info.df$transition_count)))
poxv_state_tr_all_m_loss_ct.p <- poxv_state_tr_all_m_loss_ct.p + labs(title="All methods", subtitle="Ortholog group loss")
poxv_state_tr_all_m_loss_ct.p <- poxv_state_tr_all_m_loss_ct.p + theme(text = element_text(size = 20), plot.title = element_text(size = 30), plot.subtitle = element_text(size = 24))
poxv_state_tr_all_m_loss_ct.p <- poxv_state_tr_all_m_loss_ct.p + xlim(0,4.5)
poxv_state_tr_all_m_loss_ct.p

# pdf("Figures/gain_and_loss_summarized.pdf", width=24, height=10)
# print(poxv_state_tr_all_m_gain_ct.p + poxv_state_tr_all_m_loss_ct.p + plot_layout(nrow=1))
# dev.off()


########## Analyze rate of net change in ortholog group numbers identified by MPPA, MAP, and DOWNPASS
poxv_no_gain_loss_rt.p_list <- list()
poxv_root_to_tip_pooled_m_ct.p_list <- list()

for (a in c("MPPA","MAP","DOWNPASS")) {
  poxv_state_tr_pooled_ct.df <- poxv_state_tr_pooled.df %>% filter(parent!=node, transition %in% c("0->1", "1->0"), method==a) %>% group_by(parent, node, transition) %>% summarise(transition_count = n()) %>% ungroup()

  ##### Rates of ortholog group gain identified by MPPA, MAP, and DOWNPASS methods
  poxv_state_tr_pooled_gain_ct.info.df <- poxv_t.rt.info.df %>% left_join(poxv_state_tr_pooled_ct.df %>% filter(transition == "0->1") %>% select(parent, node, transition_count), by = c("parent", "node"))
  poxv_state_tr_pooled_gain_ct.info.df <- poxv_state_tr_pooled_gain_ct.info.df %>% mutate(transition_count=ifelse(is.na(transition_count), 0, transition_count))
  poxv_state_tr_pooled_gain_ct.info.df <- poxv_state_tr_pooled_gain_ct.info.df %>% mutate(rate_of_change=transition_count/branch.length)
  
  if (nrow(poxv_state_tr_pooled_gain_ct.info.df %>% filter(rate_of_change < 1 & rate_of_change > -1 & rate_of_change!=0)) > 0) {
    print("Error: Rate of net change falls in the range of (-1,1)")
    break
    }
  
  poxv_state_tr_pooled_gain_rt.p <- ggtree(poxv_t.rt.rn, linewidth=2, aes(color=transition_count/branch.length))
  poxv_state_tr_pooled_gain_rt.p <- poxv_state_tr_pooled_gain_rt.p %<+% poxv_state_tr_pooled_gain_ct.info.df
  poxv_state_tr_pooled_gain_rt.p <- poxv_state_tr_pooled_gain_rt.p + theme_tree()
  # poxv_state_tr_pooled_gain_rt.p <- poxv_state_tr_pooled_gain_rt.p + geom_text2(aes(label=round((transition_count/branch.length), 1)), color="black", size=5, hjust=1.3)
  poxv_state_tr_pooled_gain_rt.p <- poxv_state_tr_pooled_gain_rt.p + scale_color_gradient("transition_count", name="Rate of ortholog group gain", low="black", high="yellow", breaks=seq(0, 10000, by=500))
  poxv_state_tr_pooled_gain_rt.p <- poxv_state_tr_pooled_gain_rt.p + labs(title=a, subtitle="Rate of ortholog group gain")
  poxv_state_tr_pooled_gain_rt.p <- poxv_state_tr_pooled_gain_rt.p + theme(text = element_text(size = 24), plot.title = element_text(size = 30), plot.subtitle = element_text(size = 16))
  poxv_state_tr_pooled_gain_rt.p <- poxv_state_tr_pooled_gain_rt.p + xlim(0,4.5)

  ##### Rates of ortholog group loss identified by MPPA, MAP, and DOWNPASS methods
  poxv_state_tr_pooled_loss_ct.info.df <- poxv_t.rt.info.df %>% left_join(poxv_state_tr_pooled_ct.df %>% filter(transition == "1->0") %>% select(parent, node, transition_count), by = c("parent", "node"))
  poxv_state_tr_pooled_loss_ct.info.df <- poxv_state_tr_pooled_loss_ct.info.df %>% mutate(transition_count=ifelse(is.na(transition_count), 0, transition_count))
  
  poxv_state_tr_pooled_loss_rt.p <- ggtree(poxv_t.rt.rn, linewidth=2, aes(color=transition_count/branch.length))
  poxv_state_tr_pooled_loss_rt.p <- poxv_state_tr_pooled_loss_rt.p %<+% poxv_state_tr_pooled_loss_ct.info.df
  poxv_state_tr_pooled_loss_rt.p <- poxv_state_tr_pooled_loss_rt.p + theme_tree()
  # poxv_state_tr_pooled_loss_rt.p <- poxv_state_tr_pooled_loss_rt.p + geom_text2(aes(label=round((transition_count/branch.length), 1)), color="black", size=5, hjust=1.3)
  poxv_state_tr_pooled_loss_rt.p <- poxv_state_tr_pooled_loss_rt.p + scale_color_gradient("transition_count", name="Rate of ortholog group loss", low="black", high="deepskyblue", breaks=seq(0, 10000, by=1000))
  poxv_state_tr_pooled_loss_rt.p <- poxv_state_tr_pooled_loss_rt.p + labs(title=a, subtitle="Rate of ortholog group loss")
  poxv_state_tr_pooled_loss_rt.p <- poxv_state_tr_pooled_loss_rt.p + theme(text = element_text(size = 24), plot.title = element_text(size = 30), plot.subtitle = element_text(size = 16))
  poxv_state_tr_pooled_loss_rt.p <- poxv_state_tr_pooled_loss_rt.p + xlim(0,4.5)

  ##### Rates of net change in ortholog group number identified by MPPA, MAP, and DOWNPASS methods
  poxv_state_tr_pooled_total_ct.df <- poxv_state_tr_pooled_ct.df %>% mutate(transition_count=ifelse(transition=="1->0", transition_count*-1, transition_count))
  poxv_state_tr_pooled_total_ct.df <- poxv_state_tr_pooled_total_ct.df %>% group_by(parent, node) %>% summarise(transition_total = sum(transition_count)) %>% ungroup()
  poxv_state_tr_pooled_total_ct.info.df <- poxv_t.rt.info.df %>% left_join(poxv_state_tr_pooled_total_ct.df, by = c("parent", "node"))
  poxv_state_tr_pooled_total_ct.info.df <- poxv_state_tr_pooled_total_ct.info.df %>% mutate(transition_total=ifelse(is.na(transition_total), 0, transition_total))

  poxv_state_tr_pooled_m_total_rt.p <- ggtree(poxv_t.rt.rn, linewidth=2, aes(color=(log10(abs(transition_total/branch.length))*sign(transition_total))))
  poxv_state_tr_pooled_m_total_rt.p <- poxv_state_tr_pooled_m_total_rt.p %<+% poxv_state_tr_pooled_total_ct.info.df
  poxv_state_tr_pooled_m_total_rt.p <- poxv_state_tr_pooled_m_total_rt.p + theme_tree()
  # poxv_state_tr_pooled_m_total_rt.p <- poxv_state_tr_pooled_m_total_rt.p + geom_text2(aes(label=round((transition_count/branch.length), 1)), color="black", size=5, hjust=1.3)
  poxv_state_tr_pooled_m_total_rt.p <- poxv_state_tr_pooled_m_total_rt.p + scale_color_gradient2("transition_total", name="log10 rate of net change in ortholog group number", low="deepskyblue", mid="black", midpoint=0, high="yellow", na.value="grey", breaks=seq(-3, 3, by=1))
  poxv_state_tr_pooled_m_total_rt.p <- poxv_state_tr_pooled_m_total_rt.p + labs(title=a, subtitle="Rate of net change in ortholog group number")
  poxv_state_tr_pooled_m_total_rt.p <- poxv_state_tr_pooled_m_total_rt.p + theme(text = element_text(size = 24), plot.title = element_text(size = 30), plot.subtitle = element_text(size = 16))
  poxv_state_tr_pooled_m_total_rt.p <- poxv_state_tr_pooled_m_total_rt.p + xlim(0,4.5)

  poxv_no_gain_loss_rt.p_list[[a]] <- poxv_state_tr_pooled_gain_rt.p + poxv_state_tr_pooled_loss_rt.p + poxv_state_tr_pooled_m_total_rt.p + plot_layout(ncol=1)
  
  ### Measure changes in gene number from root to tip estimated by MPPA, MAP, and DOWNPASS methods
  tip_nodes <- setdiff(poxv_state_tr_pooled_total_ct.info.df$node, poxv_state_tr_pooled_total_ct.info.df$parent)
  poxv_t.rt.rn_root_node <- find_root(poxv_t.rt.rn)
  
  for (i in 1:length(tip_nodes)) {
    node_path <- nodepath(poxv_t.rt.rn, from = poxv_t.rt.rn_root_node, to = tip_nodes[i])
    transition_total_count <- 0
    for (j in 1:(length(node_path)-1)) {
      poxv_state_tr_pooled_total_ct.info.df.filtered <- poxv_state_tr_pooled_total_ct.info.df %>% filter(parent==node_path[j], node==node_path[j+1])
      transition_total_count <- transition_total_count + poxv_state_tr_pooled_total_ct.info.df.filtered$transition_total
    }
    if (i==1) {
      poxv_root_to_tip_pooled_m_ct.df <- data.frame(node=tip_nodes[i], transition_count_to_tip=transition_total_count)
    } else {
      poxv_root_to_tip_pooled_m_ct.df <- bind_rows(poxv_root_to_tip_pooled_m_ct.df, data.frame(node=tip_nodes[i], transition_count_to_tip=transition_total_count))
    }
  }
  
  ########## ***** Supplementary Figure 2C: Change in ortholog group number identified by each method *****
  poxv_root_to_tip_pooled_m_ct.df <- poxv_root_to_tip_pooled_m_ct.df %>% mutate(species_name=(as.data.frame(poxv_state_tr_pooled_total_ct.info.df))[match(node, (as.data.frame(poxv_state_tr_pooled_total_ct.info.df))$node), "label"])
  poxv_root_to_tip_pooled_m_ct.df <- poxv_root_to_tip_pooled_m_ct.df %>% mutate(Subfamily=poxv_ncbi_ds.df[match(species_name, poxv_ncbi_ds.df$tip_label_dummy), "Subfamily"], clade=poxv_ncbi_ds.df[match(species_name, poxv_ncbi_ds.df$tip_label_dummy), "clade"])
  poxv_root_to_tip_pooled_m_ct.df <- poxv_root_to_tip_pooled_m_ct.df %>% arrange(Subfamily, clade, species_name)
  poxv_root_to_tip_pooled_m_ct.df$species_name <- factor(poxv_root_to_tip_pooled_m_ct.df$species_name, levels=ggtree(poxv_t.rt.rn)$data %>% filter(isTip) %>% arrange(desc(y)) %>% pull(label))
  
  poxv_root_to_tip_pooled_m_ct.p <- ggplot(poxv_root_to_tip_pooled_m_ct.df, aes(x=species_name, y=transition_count_to_tip, fill=clade))
  poxv_root_to_tip_pooled_m_ct.p <- poxv_root_to_tip_pooled_m_ct.p + geom_bar(stat="identity", width=1)
  poxv_root_to_tip_pooled_m_ct.p <- poxv_root_to_tip_pooled_m_ct.p + geom_text(aes(label=substr(clade, 1, 3), col=clade), vjust=-0.8, size=7)
  poxv_root_to_tip_pooled_m_ct.p <- poxv_root_to_tip_pooled_m_ct.p + scale_fill_manual("clade", name="Clade", values=poxv_color.clade.v, breaks=names(poxv_color.clade.v), na.value = NA)
  poxv_root_to_tip_pooled_m_ct.p <- poxv_root_to_tip_pooled_m_ct.p + scale_color_manual("clade", name="Clade", values=poxv_color.clade.v, breaks=names(poxv_color.clade.v), na.value = NA)
  poxv_root_to_tip_pooled_m_ct.p <- poxv_root_to_tip_pooled_m_ct.p + labs(title="Cumulative net change in ortholog group number from root to each tip", 
                                                                          subtitle=a, x="Poxvirus", y="Cumulative net change in ortholog group number from root to each tip")
  poxv_root_to_tip_pooled_m_ct.p <- poxv_root_to_tip_pooled_m_ct.p + guides(fill=guide_legend(ncol=1), col="none")
  poxv_root_to_tip_pooled_m_ct.p <- poxv_root_to_tip_pooled_m_ct.p + theme_classic()
  poxv_root_to_tip_pooled_m_ct.p <- poxv_root_to_tip_pooled_m_ct.p + theme(text=element_text(size = 18), plot.title=element_text(size = 32), axis.text.x=element_blank())
  poxv_root_to_tip_pooled_m_ct.p <- poxv_root_to_tip_pooled_m_ct.p + scale_y_continuous(breaks=seq(-50, 250, by=50), limits=c(-50,250))
  
  poxv_root_to_tip_pooled_m_ct.p_list[[a]] <- poxv_root_to_tip_pooled_m_ct.p
}

# pdf("Figures/gain_and_loss_rate_all_methods.pdf", width=24, height=20)
# wrap_plots(poxv_no_gain_loss_rt.p_list) + plot_layout(nrow=1)
# dev.off()
# 
# pdf("Figures/total_change_in_gene_number_all_methods.pdf", width=24, height=20)
# wrap_plots(poxv_root_to_tip_pooled_m_ct.p_list) + plot_layout(ncol=1, guides="collect")
# dev.off()


########## Plot branchwise rates of net change in gene number events correspondingly determined by all methods
########## ***** Figure 1F: Branchwise rates of net change in ortholog group number correspondingly determined by all methods *****
poxv_state_tr_all_m_total_ct.info.df <- poxv_t.rt.info.df %>% left_join(poxv_state_tr_all_m_total_ct.df, by = c("parent", "node"))
poxv_state_tr_all_m_total_ct.info.df <- poxv_state_tr_all_m_total_ct.info.df %>% mutate(transition_total=ifelse(is.na(transition_total), 0, transition_total)) # All rates are < -1 and > 1
poxv_state_tr_all_m_total_ct.info.df <- poxv_state_tr_all_m_total_ct.info.df %>% mutate(rate_of_change=transition_total/branch.length)

if (nrow(poxv_state_tr_all_m_total_ct.info.df %>% filter(rate_of_change < 1 & rate_of_change > -1 & rate_of_change!=0)) > 0) {print("Error: Rate of net change falls in the range of (-1,1)")}

poxv_state_tr_all_m_total_rt.p <- ggtree(poxv_t.rt.rn, linewidth=2, aes(color=(log10(abs(transition_total/branch.length)+1)*sign(transition_total))))
poxv_state_tr_all_m_total_rt.p <- poxv_state_tr_all_m_total_rt.p %<+% poxv_state_tr_all_m_total_ct.info.df
poxv_state_tr_all_m_total_rt.p <- poxv_state_tr_all_m_total_rt.p + theme_tree()
poxv_state_tr_all_m_total_rt.p <- poxv_state_tr_all_m_total_rt.p + geom_text2(aes(label=round((transition_total/branch.length), 0)), color="black", size=5, hjust=1.3)
poxv_state_tr_all_m_total_rt.p <- poxv_state_tr_all_m_total_rt.p + scale_fill_manual("clade", name="Clade", values=poxv_color.clade.v, breaks=names(poxv_color.clade.v), na.value = NA)
poxv_state_tr_all_m_total_rt.p <- poxv_state_tr_all_m_total_rt.p + scale_color_gradient2("transition_total", name="Rate of net change\nin ortholog group number", low="deepskyblue", mid="black", midpoint=0, high="yellow", na.value="grey", breaks=seq(-3, 3, by=1))
poxv_state_tr_all_m_total_rt.p <- poxv_state_tr_all_m_total_rt.p + new_scale_color()

##### Measure net changes in ortholog group number from root to each tip estimated correspondingly by all methods
tip_nodes <- setdiff(poxv_state_tr_all_m_total_ct.info.df$node, poxv_state_tr_all_m_total_ct.info.df$parent)
poxv_t.rt.rn_root_node <- find_root(poxv_t.rt.rn)

for (i in 1:length(tip_nodes)) {
  node_path <- nodepath(poxv_t.rt.rn, from = poxv_t.rt.rn_root_node, to = tip_nodes[i])
  transition_total_count <- 0
  for (j in 1:(length(node_path)-1)) {
    poxv_state_tr_all_m_total_ct.info.df.filtered <- poxv_state_tr_all_m_total_ct.info.df %>% filter(parent==node_path[j], node==node_path[j+1])
    transition_total_count <- transition_total_count + poxv_state_tr_all_m_total_ct.info.df.filtered$transition_total
  }
  if (i==1) {
    poxv_root_to_tip_all_m_ct.df <- data.frame(node=tip_nodes[i], transition_count_to_tip=transition_total_count)
  } else {
    poxv_root_to_tip_all_m_ct.df <- bind_rows(poxv_root_to_tip_all_m_ct.df, data.frame(node=tip_nodes[i], transition_count_to_tip=transition_total_count))
  }
}

poxv_root_to_tip_all_m_ct.df <- poxv_root_to_tip_all_m_ct.df %>% mutate(species_name=(as.data.frame(poxv_state_tr_all_m_total_ct.info.df))[match(node, (as.data.frame(poxv_state_tr_all_m_total_ct.info.df))$node), "label"])
poxv_root_to_tip_all_m_ct.df <- poxv_root_to_tip_all_m_ct.df %>% mutate(Subfamily=poxv_ncbi_ds.df[match(species_name, poxv_ncbi_ds.df$tip_label_dummy), "Subfamily"], clade=poxv_ncbi_ds.df[match(species_name, poxv_ncbi_ds.df$tip_label_dummy), "clade"])
poxv_root_to_tip_all_m_ct.df <- poxv_root_to_tip_all_m_ct.df %>% arrange(Subfamily, clade, species_name)
poxv_root_to_tip_all_m_ct.df$species_name <- factor(poxv_root_to_tip_all_m_ct.df$species_name, levels=poxv_root_to_tip_all_m_ct.df$species_name)

poxv_state_tr_all_m_total_rt.p <- poxv_state_tr_all_m_total_rt.p + geom_fruit(data = poxv_root_to_tip_all_m_ct.df %>% rename(label = species_name) %>% rename(clade_bar = clade), 
                                                                              geom = geom_col, mapping = aes(x = transition_count_to_tip, y = label, fill = clade_bar),
                                                                              orientation = "y", offset = 0.15, width = 0.6, pwidth = 0.2, grid.params=list(),
                                                                              axis.params = list(axis = "x", title = "Cumulative net change in ortholog group number\nfrom root to each tip", title.size = 6, text.size = 4, text.angle = 90, hjust = 1, limits = c(0, 150)))
poxv_state_tr_all_m_total_rt.p <- poxv_state_tr_all_m_total_rt.p + theme(text = element_text(size = 20), plot.title = element_text(size = 30), plot.subtitle = element_text(size = 24))
poxv_state_tr_all_m_total_rt.p <- poxv_state_tr_all_m_total_rt.p + xlim(0,6)
poxv_state_tr_all_m_total_rt.p <- poxv_state_tr_all_m_total_rt.p + guides(fill=guide_legend(ncol=1), color="none")
poxv_state_tr_all_m_total_rt.p <- poxv_state_tr_all_m_total_rt.p + labs(title="All methods", subtitle="Rate of net change in ortholog group number", fill = "Clade", x = "Ortholog group number")
poxv_state_tr_all_m_total_rt.p <- poxv_state_tr_all_m_total_rt.p + coord_cartesian(clip="off")
poxv_state_tr_all_m_total_rt.p

# pdf("Figures/gain_and_loss_rate_summarized.pdf", width=16, height=12)
# poxv_state_tr_all_m_total_rt.p
# dev.off()


########## Assign status of ortholog group gain and loss to all poxvirus genes
##### Gained genes
####### The status of root is 0
##### Genes with gain events (for analysis of domain gain and loss associated with ortholog group gain and loss)
####### 1) There is a 0->1 state change in the tree
####### 2) The status of root is 0
##### Genes with loss events (for analysis of domain gain and loss associated with ortholog group gain and loss)
####### 1) There is a 1->0 state change in the tree
####### 2) The status of MRCA is 1 and there is a missing descendant

poxv_state_root_status.df <- as.data.frame(poxv_state_tr_pooled.df %>% filter(parent==node) %>% group_by(parent, node, state_parent, og_name) %>% summarise(state_count = n()) %>% ungroup() %>% mutate(state_count = state_count/3) %>% filter(state_count == 1))

poxv_state_absent_tip.df <- poxv_state_tr_pooled.df %>% group_by(parent, node, transition, og_name) %>% summarise(transition_count = n()) %>% ungroup()
poxv_state_absent_tip.df <- poxv_state_absent_tip.df %>% mutate(transition_count = transition_count/3) %>% filter(transition_count == 1) # Pick only transition events estimated correspondingly by all methods

poxv_state_absent_tip_ct.df <- poxv_state_tr_pooled.df %>% filter(isTip==TRUE, method=="MPPA", state_desc == "Present")
poxv_state_absent_tip_ct.df <- poxv_state_absent_tip_ct.df %>% group_by(og_name) %>% summarise(tip_set = list(unique(node)), num_tip = length(tip_set[[1]]), .groups = "drop")
poxv_state_absent_tip_ct.df <- poxv_state_absent_tip_ct.df %>% mutate(ca_node = map_int(tip_set, ~ (if (length(.x) < 2) {.x[[1]]} else {getMRCA(poxv_t.rt.rn, .x)})),
                                                                      num_desc = map_int(ca_node, ~ length(Descendants(poxv_t.rt.rn, .x, type = "tips")[[1]])))
poxv_state_absent_tip_ct.df <- poxv_state_absent_tip_ct.df %>% left_join(poxv_state_absent_tip.df %>% rename(ca_node = node) %>% select(og_name, ca_node, transition), by = c("og_name", "ca_node"))
poxv_state_absent_tip_ct.df <- poxv_state_absent_tip_ct.df %>% mutate(ca_node_status = ifelse(!is.na(transition), str_split(transition, "->", simplify=TRUE)[,2], transition))
poxv_state_absent_tip_ct.df <- poxv_state_absent_tip_ct.df %>% mutate(loss_status = ifelse(!is.na(ca_node_status) & ca_node_status==1, ifelse(num_tip < num_desc, "some_absent", "conserved_after_gain"), "undeterminable")) # Consider only the case that CA is 1
poxv_state_absent_tip_ct.df <- poxv_state_absent_tip_ct.df %>% mutate(loss_status = ifelse(ca_node==find_root(poxv_t.rt.rn) & !is.na(ca_node_status) & ca_node_status==1 & num_desc==length(poxv_t.rt.rn$tip.label), "conserved", loss_status))

poxv_og_gain_loss_status.df <- data.frame(poxv_og_name=gsub("\\.", "", poxv_vacv_gene_name.df$og_id))
poxv_og_gain_loss_status.df <- poxv_og_gain_loss_status.df %>% mutate(root_status=poxv_state_root_status.df$state_parent[match(poxv_og_name, poxv_state_root_status.df$og_name)],
                                                                      tip_status=poxv_state_absent_tip_ct.df$loss_status[match(poxv_og_name, poxv_state_absent_tip_ct.df$og_name)])
poxv_og_gain_loss_status.df <- poxv_og_gain_loss_status.df %>% mutate(detected_gain_event=ifelse(poxv_og_name %in% (poxv_state_tr_all_m.df %>% filter(transition=="0->1"))$og_name, "yes", "no"),
                                                                      detected_loss_event=ifelse(poxv_og_name %in% (poxv_state_tr_all_m.df %>% filter(transition=="1->0"))$og_name, "yes", "no"))
poxv_og_gain_loss_status.df <- poxv_og_gain_loss_status.df %>% mutate(gain_loss_status=ifelse(tip_status=="conserved", "conserved", "undeterminable"))
poxv_og_gain_loss_status.df <- poxv_og_gain_loss_status.df %>% mutate(gain_loss_status=ifelse(detected_gain_event=="yes" | (!is.na(root_status) & root_status==0),
                                                                                       ifelse(detected_loss_event=="yes" | tip_status=="some_absent", "both", "gain"),
                                                                                       ifelse(detected_loss_event=="yes" | tip_status=="some_absent", "loss", gain_loss_status)))
poxv_og_gain_loss_status.df <- poxv_og_gain_loss_status.df %>% mutate(gained_status=ifelse(!is.na(root_status) & root_status==0, "gained",
                                                                                        ifelse(!is.na(root_status) & root_status==1, "pre-existing","undeterminable")))
poxv_og_gain_loss_status.df$gain_loss_status <- factor(poxv_og_gain_loss_status.df$gain_loss_status, levels=c("gain","loss","both","conserved","undeterminable"))
poxv_og_gain_loss_status.df$gained_status <- factor(poxv_og_gain_loss_status.df$gained_status, levels=c("gained","pre-existing","undeterminable"))

table(poxv_og_gain_loss_status.df$gained_status) # Number of gained, pre-existing, and undeterminable ortholog groups


########## Calculate number of ortholog group loss events per branch length in gained and pre-existing ortholog groups
##### Calculate sum of branch length of gained and pre-existing ortholog groups where no loss events occurred
poxv_og_branch_l.df <- poxv_state_tr_all_m.df %>% mutate(ca_node=poxv_state_absent_tip_ct.df$ca_node[match(og_name, poxv_state_absent_tip_ct.df$og_name)])
poxv_og_branch_l.df <- poxv_og_branch_l.df %>% filter(transition=="1->1")

poxv_tree_branch_l.df <- data.frame(parent=poxv_t.rt.rn$edge[,1], node=poxv_t.rt.rn$edge[,2], branch.length=poxv_t.rt.rn$edge.length)
poxv_tree_branch_l.df <- poxv_tree_branch_l.df %>% mutate(branch_id=paste(parent, node, sep="_"))

poxv_og_branch_l.df <- poxv_og_branch_l.df %>% mutate(branch_id=paste(parent, node, sep="_"))
poxv_og_branch_l.df <- poxv_og_branch_l.df %>% mutate(branch.length=poxv_tree_branch_l.df$branch.length[match(branch_id, poxv_tree_branch_l.df$branch_id)])
poxv_og_branch_l.df <- poxv_og_branch_l.df %>% rowwise() %>% mutate(is_under_ca=if(is.na(ca_node)) FALSE else parent %in% c(ca_node, Descendants(poxv_t.rt.rn, ca_node, type="all"))) %>% ungroup()
poxv_og_branch_l.df <- poxv_og_branch_l.df %>% mutate(og_name=gsub("N0HOG", "N0.HOG", og_name, fixed=TRUE))
poxv_og_branch_l.df <- poxv_og_branch_l.df %>% mutate(gained_status=poxv_og_gain_loss_status.df$gained_status[match(og_name, gsub("N0HOG", "N0.HOG", poxv_og_gain_loss_status.df$poxv_og_name, fixed=TRUE))])
poxv_og_branch_l.df <- poxv_og_branch_l.df %>% filter(is_under_ca==TRUE, gained_status %in% c("gained","pre-existing"))

poxv_og_branch_l_sum.df <- poxv_og_branch_l.df %>% group_by(og_name, gained_status) %>% summarise(weight=sum(branch.length, na.rm=TRUE)) %>% ungroup()

poxv_gained_status_branch_l.df <- poxv_og_branch_l_sum.df %>% group_by(gained_status) %>% summarise(weight=sum(weight, na.rm=TRUE)) %>% ungroup()

##### Calculate number of ortholog group loss events per branch length
poxv_state_tr_all_m_loss_norm.df <- poxv_state_tr_all_m.df %>% mutate(og_name=gsub("N0HOG", "N0.HOG", og_name, fixed=TRUE))
poxv_state_tr_all_m_loss_norm.df <- poxv_state_tr_all_m_loss_norm.df %>% filter(transition=="1->0")

poxv_state_tr_all_m_loss_norm_ct.df <- poxv_state_tr_all_m_loss_norm.df %>% group_by(og_name) %>% summarise(count=n()) %>% ungroup()
poxv_state_tr_all_m_loss_norm_ct.df <- poxv_og_branch_l_sum.df %>% select(og_name, gained_status) %>% left_join(poxv_state_tr_all_m_loss_norm_ct.df, by="og_name")
poxv_state_tr_all_m_loss_norm_ct.df <- poxv_state_tr_all_m_loss_norm_ct.df %>% mutate(count=ifelse(is.na(count), 0, count))

poxv_loss_events_norm.df <- poxv_state_tr_all_m_loss_norm_ct.df %>% group_by(gained_status) %>% summarise(count=sum(count)) %>% ungroup()
poxv_loss_events_norm.df <- poxv_gained_status_branch_l.df %>% left_join(poxv_loss_events_norm.df, by="gained_status")
poxv_loss_events_norm.df <- poxv_loss_events_norm.df %>% mutate(count=ifelse(is.na(count), 0, count), norm_count=count/weight)

########## ***** Figure 2J: Proportion of ortholog group loss events in gained and pre-existing ortholog groups *****
poxv_loss_events_bar.p <- ggplot(poxv_loss_events_norm.df, aes(x=gained_status, y=norm_count, fill=gained_status))
poxv_loss_events_bar.p <- poxv_loss_events_bar.p + geom_col(width=0.5)
poxv_loss_events_bar.p <- poxv_loss_events_bar.p + scale_fill_manual("gained_status", values = c("gained"="orange","pre-existing"="navyblue","undeterminable"="lightgrey"))
poxv_loss_events_bar.p <- poxv_loss_events_bar.p + labs(title="Number of ortholog group loss events per branch length", x="Category", y="Loss events per branch length")
poxv_loss_events_bar.p <- poxv_loss_events_bar.p + theme_classic()
poxv_loss_events_bar.p <- poxv_loss_events_bar.p + theme(text = element_text(size = 18), plot.title = element_text(size = 24), plot.subtitle = element_text(size = 22))
poxv_loss_events_bar.p <- poxv_loss_events_bar.p + scale_y_continuous(breaks=seq(0, 0.4, by=0.1), limits=c(0,0.4))
poxv_loss_events_bar.p

# pdf("Figures/norm_loss_ct_sep_by_gained_status.pdf", width=16, height=12)
# poxv_loss_events_bar.p
# dev.off()

  

############################## Gene duplication analysis ##############################
gene_dup_events.f <- "251208_Gene_Duplication_Events/Duplications_HOGs.tsv"
gene_dup_events.df <- read.delim(gene_dup_events.f, check.names = FALSE, header = T)

colnames(gene_dup_events.df)[2] <- "Species_Tree_Node"

gene_dup_events_ct.df <- gene_dup_events.df %>% group_by(Species_Tree_Node) %>% summarise(count = n()) %>% ungroup()

poxv_t.rt.rn.wt_of2_node.f <- "251208_Gene_Duplication_Events/SpeciesTree_rooted_node_labels.txt"
poxv_t.rt.rn.wt_of2_node <- read.tree(poxv_t.rt.rn.wt_of2_node.f) # The tree is already rooted.

poxv_t.rt.rn.wt_of2_node.info.df <- ggtree(poxv_t.rt.rn.wt_of2_node)$data

poxv_t.rt.rn.wt_of2_node$tip.label <- ifelse(poxv_t.rt.rn.wt_of2_node$tip.label=="GCF_000857045.1_Monkeypox_virus", "GCF_000857045.1_Monkeypox_virus_Zaire-96-I-16_Ia",
                                      ifelse(poxv_t.rt.rn.wt_of2_node$tip.label=="GCA_039269415.1_Monkeypox_virus", "GCA_039269415.1_Monkeypox_virus_RDC-NKV-GOM-MPOX-010_Ib",
                                      ifelse(poxv_t.rt.rn.wt_of2_node$tip.label=="GCA_006465585.1_Monkeypox_virus", "GCA_006465585.1_Monkeypox_virus_Sierra_Leone_IIa",
                                      ifelse(poxv_t.rt.rn.wt_of2_node$tip.label=="GCF_014621545.1_Monkeypox_virus", "GCF_014621545.1_Monkeypox_virus_MPXV-M5312_HM12_Rivers_IIb", poxv_t.rt.rn.wt_of2_node$tip.label))))
poxv_t.rt.rn.wt_of2_node$tip.label <- unname(sapply(poxv_t.rt.rn.wt_of2_node$tip.label, function(x) {
  parts <- str_split(x, "_", simplify = TRUE)
  paste(parts[3:length(parts)], collapse = "_")  # Join everything after the first two parts
}))

poxv_t.rt.rn.wt_of2_node.info.df <- poxv_t.rt.rn.wt_of2_node.info.df %>% mutate(dup_count=gene_dup_events_ct.df$count[match(label, gene_dup_events_ct.df$Species_Tree_Node)])
poxv_t.rt.rn.wt_of2_node.info.df <- poxv_t.rt.rn.wt_of2_node.info.df %>% mutate(dup_count=ifelse(is.na(dup_count), 0, dup_count))

poxv_t.rt.rn.wt_of2_node.info.df$label <- ifelse(poxv_t.rt.rn.wt_of2_node.info.df$label=="GCF_000857045.1_Monkeypox_virus", "GCF_000857045.1_Monkeypox_virus_Zaire-96-I-16_Ia",
                                          ifelse(poxv_t.rt.rn.wt_of2_node.info.df$label=="GCA_039269415.1_Monkeypox_virus", "GCA_039269415.1_Monkeypox_virus_RDC-NKV-GOM-MPOX-010_Ib",
                                          ifelse(poxv_t.rt.rn.wt_of2_node.info.df$label=="GCA_006465585.1_Monkeypox_virus", "GCA_006465585.1_Monkeypox_virus_Sierra_Leone_IIa",
                                          ifelse(poxv_t.rt.rn.wt_of2_node.info.df$label=="GCF_014621545.1_Monkeypox_virus", "GCF_014621545.1_Monkeypox_virus_MPXV-M5312_HM12_Rivers_IIb", poxv_t.rt.rn.wt_of2_node.info.df$label))))
poxv_t.rt.rn.wt_of2_node.info.df$label <- unname(sapply(poxv_t.rt.rn.wt_of2_node.info.df$label, function(x) {
  parts <- str_split(x, "_", simplify = TRUE)
  if (dim(parts)[2]==1) {x} else {paste(parts[3:length(parts)], collapse = "_")}  # Join everything after the first two parts
}))

poxv_t.rt.rn.wt_of2_node.info.df <- poxv_t.rt.rn.wt_of2_node.info.df %>% mutate(clade=poxv_ncbi_ds.df[match(label, poxv_ncbi_ds.df$tip_label_dummy), "clade"])

########## ***** Figure 2G, left: Number of gene duplication events in each branch *****
gene_dup_events.p <- ggtree(poxv_t.rt.rn.wt_of2_node, linewidth=1, aes(color=dup_count))
gene_dup_events.p <- gene_dup_events.p %<+% poxv_t.rt.rn.wt_of2_node.info.df
gene_dup_events.p <- gene_dup_events.p + geom_point2(aes(subset=(node==find_root(poxv_t.rt.rn)), fill=dup_count, size=5))
# gene_dup_events.p <- gene_dup_events.p + geom_text2(aes(subset=!isTip, label=paste(label, dup_count, sep=":")), size=5, hjust=1.3)
gene_dup_events.p <- gene_dup_events.p + scale_color_gradient("dup_count", trans="pseudo_log", name="Number of\nduplication\nevents", low="black", high="red", na.value="black", breaks=c(0,3,5,10,30,50,100,max(poxv_t.rt.rn.wt_of2_node.info.df$dup_count)))
gene_dup_events.p <- gene_dup_events.p + theme_tree()
gene_dup_events.p <- gene_dup_events.p + theme(text = element_text(size = 20))
gene_dup_events.p <- gene_dup_events.p + geom_treescale(width=0.1)
gene_dup_events.p <- gene_dup_events.p + labs(title="Number of gene duplication events")
gene_dup_events.p <- gene_dup_events.p + coord_cartesian(clip="off")
gene_dup_events.p

# pdf("Figures/duplication_events.pdf", width=16, height=12)
# gene_dup_events.p
# dev.off()


########## Calculate number of gene duplication events per branch length in gained and pre-existing ortholog groups
##### Note: This calculation does not consider duplication of ortholog groups
gene_dup_events_sep.df <- gene_dup_events.df %>% filter(HOG_genes_1==HOG_genes_2)
gene_dup_events_sep.df <- gene_dup_events_sep.df %>% mutate(gained_status=poxv_og_gain_loss_status.df$gained_status[match(HOG_genes_1, gsub("N0HOG", "N0.HOG", poxv_og_gain_loss_status.df$poxv_og_name, fixed=TRUE))])
gene_dup_events_sep.df <- gene_dup_events_sep.df %>% filter(gained_status %in% c("gained","pre-existing"))

gene_dup_events_norm.df <- gene_dup_events_sep.df %>% group_by(gained_status) %>% summarise(count=n()) %>% ungroup()

gene_dup_events_norm.df <- poxv_gained_status_branch_l.df %>% left_join(gene_dup_events_norm.df, by="gained_status")
gene_dup_events_norm.df <- gene_dup_events_norm.df %>% mutate(count=ifelse(is.na(count), 0, count), norm_count=count/weight)

########## ***** Figure 2I: Number of gene duplication events per branch length in gained and pre-exising ortholog groups *****
gene_dup_events_bar.p <- ggplot(gene_dup_events_norm.df, aes(x=gained_status, y=norm_count, fill=gained_status))
gene_dup_events_bar.p <- gene_dup_events_bar.p + geom_col(width=0.5)
gene_dup_events_bar.p <- gene_dup_events_bar.p + scale_fill_manual("gained_status", values = c("gained"="orange","pre-existing"="navyblue","undeterminable"="lightgrey"))
gene_dup_events_bar.p <- gene_dup_events_bar.p + labs(title="Normalized count of gene duplication events per branch length", x="Category", y="Duplication events per branch length")
gene_dup_events_bar.p <- gene_dup_events_bar.p + theme_classic()
gene_dup_events_bar.p <- gene_dup_events_bar.p + theme(text = element_text(size = 18), plot.title = element_text(size = 24), plot.subtitle = element_text(size = 22))
gene_dup_events_bar.p <- gene_dup_events_bar.p + scale_y_continuous(breaks=seq(0, 3, by=1), limits=c(0,3))
gene_dup_events_bar.p

# pdf("Figures/norm_duplication_ct_sep_by_gained_status.pdf", width=16, height=12)
# gene_dup_events_bar.p
# dev.off()



############################## Analysis of domain gain and loss associated with ortholog group gain and loss in the poxvirus lineage ##############################
########## Plot results and count domain gain and loss events associated with ortholog group gain and loss events in the poxvirus lineage
poxv_sign_state_tr_pooled.df <- data.frame()
poxv_sign_state_tr.p_list <- list()

##### Plot domain gain and loss events associated with ortholog group gain and loss events in the poxvirus lineage identified by MPPA, MAP, and DOWNPASS methods
for (a in c("MPPA","MAP","DOWNPASS")) {
  poxv_pastml_rs_sign.f <- paste("251208_gain_and_loss_domain_analysis/251208_cs_matrix_pooled_additional_domain_pastml", a, "output.tsv", sep="_")
  poxv_pastml_rs_sign <- fread(poxv_pastml_rs_sign.f, header=T, sep="\t", quote="", check.names=F, data.table=FALSE)
  
  poxv_pastml_rs_sign_state.df <- poxv_pastml_rs_sign
  poxv_pastml_rs_sign_state.df[which(poxv_pastml_rs_sign_state.df$node=="Root"),"node"] <- find_root(poxv_t.rt.rn)
  
  poxv_pastml_rs_sign_state.df <- poxv_pastml_rs_sign_state.df %>% mutate(node=ifelse(node=="GCF_000857045.1_Monkeypox_virus", "GCF_000857045.1_Monkeypox_virus_Zaire-96-I-16_Ia",
                                                                               ifelse(node=="GCA_039269415.1_Monkeypox_virus", "GCA_039269415.1_Monkeypox_virus_RDC-NKV-GOM-MPOX-010_Ib",
                                                                               ifelse(node=="GCA_006465585.1_Monkeypox_virus", "GCA_006465585.1_Monkeypox_virus_Sierra_Leone_IIa",
                                                                               ifelse(node=="GCF_014621545.1_Monkeypox_virus", "GCF_014621545.1_Monkeypox_virus_MPXV-M5312_HM12_Rivers_IIb",node)))))
  poxv_pastml_rs_sign_state.df <- poxv_pastml_rs_sign_state.df %>% mutate(node=ifelse(str_detect(node, "GCF_|GCA_"), unname(sapply(node, function(x) {
    parts <- str_split(x, "_", simplify = TRUE)
    paste(parts[3:length(parts)], collapse = "_")  # Join everything after the first two parts
    })), node))
  
  poxv_sign_state_tr_dir <- paste("251208_gain_and_loss_domain_analysis/state_transitions_", a, "/", sep="")
  if (!file.exists(poxv_sign_state_tr_dir)) {dir.create(poxv_sign_state_tr_dir, recursive = TRUE)}
  
  poxv_sign_state_tr.df <- as.data.frame(ggtree(poxv_t.rt.rn)$data)
  poxv_sign_state_tr.df <- poxv_sign_state_tr.df %>% mutate(node_name = ifelse(isTip==FALSE, node, label))
  
  poxv_og_sign_name.v <- colnames(poxv_pastml_rs_sign_state.df[,2:ncol(poxv_pastml_rs_sign_state.df)])
  
  for (i in 1:length(poxv_og_sign_name.v)) {
    poxv_sign_state_tr_og_dir <- paste("251208_gain_and_loss_domain_analysis/state_transitions_", a, "/", poxv_og_sign_name.v[i], "/", sep="")
    if (!file.exists(poxv_sign_state_tr_og_dir)) {dir.create(poxv_sign_state_tr_og_dir, recursive = TRUE)}
    
    poxv_pastml_rs_sign_state_og.df <- poxv_pastml_rs_sign_state.df[,c("node",poxv_og_sign_name.v[i])]
    poxv_sign_state_tr_og.df <- poxv_sign_state_tr.df
    poxv_sign_state_tr_og.df$state_parent <- as.character(poxv_pastml_rs_sign_state_og.df[match(poxv_sign_state_tr_og.df$parent, poxv_pastml_rs_sign_state_og.df$node), poxv_og_sign_name.v[i]])
    poxv_sign_state_tr_og.df$state_desc <- as.character(poxv_pastml_rs_sign_state_og.df[match(poxv_sign_state_tr_og.df$node_name, poxv_pastml_rs_sign_state_og.df$node), poxv_og_sign_name.v[i]])
    
    if (nrow(poxv_sign_state_tr_og.df %>% filter(isTip==TRUE, state_desc=="1"))==1) {poxv_sign_state_tr_og.df$state_parent <- "0"}
    
    poxv_sign_state_tr_og.df$og_name <- poxv_og_sign_name.v[i]
    poxv_sign_state_tr_og.df$transition <- paste(poxv_sign_state_tr_og.df$state_parent, poxv_sign_state_tr_og.df$state_desc, sep="->")
    
    ### Save transition branch info
    poxv_sign_state_tr_og.df.f <- paste(poxv_sign_state_tr_og_dir, poxv_og_sign_name.v[i], "_transition_info.tsv", sep="")
    # write.table(poxv_sign_state_tr_og.df, poxv_sign_state_tr_og.df.f, col.names=T, row.names=F, sep="\t", quote=F)
    
    poxv_sign_state_tr_og.df <- poxv_sign_state_tr_og.df %>% mutate(noc_name=ifelse(node %in% poxv_noc.v, poxv_noc.df[match(node, poxv_noc.df$node), "node_name"], ""),
                                                                    state_desc=ifelse(state_desc=="0", "Absent", ifelse(state_desc=="0|1", "Undeterminable", ifelse(state_desc=="1", "Present", NA))))
    
    poxv_sign_state_tr_og.p <- ggtree(poxv_t.rt.rn, size=0.25)
    poxv_sign_state_tr_og.p <- poxv_sign_state_tr_og.p %<+% poxv_sign_state_tr_og.df
    poxv_sign_state_tr_og.p <- poxv_sign_state_tr_og.p + theme_tree()
    poxv_sign_state_tr_og.p <- poxv_sign_state_tr_og.p + geom_tiplab(size=2, align=FALSE, linesize=.5)
    poxv_sign_state_tr_og.p <- poxv_sign_state_tr_og.p + geom_point2(shape=16, size=3, aes(color=state_desc))
    poxv_sign_state_tr_og.p <- poxv_sign_state_tr_og.p + scale_color_manual("state_desc", name="State", values=c("lightsteelblue1","grey","red"), breaks=c("Absent","Undeterminable","Present"), na.value = NA)
    poxv_sign_state_tr_og.p <- poxv_sign_state_tr_og.p + labs(title=paste(poxv_og_sign_name.v[i], ": ", (poxv_vacv_gene_name.df %>% filter(og_id_dummy==str_split(poxv_og_sign_name.v[i], '\\.', simplify = TRUE)[,1]))$Vaccinia_virus, sep=""), 
                                                              subtitle=paste(a, ": Transition no. = ", nrow(poxv_sign_state_tr_og.df %>% filter(transition %in% c("0->1", "1->0"))), sep=""))
    poxv_sign_state_tr_og.p <- poxv_sign_state_tr_og.p + theme(text = element_text(size = 20), plot.title = element_text(size = 20), plot.subtitle = element_text(size = 16))
    poxv_sign_state_tr_og.p <- poxv_sign_state_tr_og.p + xlim(0,4.5)
    poxv_sign_state_tr_og.p <- poxv_sign_state_tr_og.p + geom_treescale(width=0.1)
    
    # ### Save the tree plot
    # pdf(paste(poxv_sign_state_tr_og_dir, poxv_og_sign_name.v[i], "_transition_plot.pdf", sep=""), width=12, height=12)
    # poxv_sign_state_tr_og_to_save.p <- poxv_sign_state_tr_og.p + geom_text2(aes(subset=(transition %in% c("0->1", "1->0")), label=transition), size=3, hjust=1.3)
    # print(poxv_sign_state_tr_og.p)
    # dev.off()
    
    ### Pool data from summarizing
    poxv_sign_state_tr_og.df <- poxv_sign_state_tr_og.df %>% mutate(method = a)
    poxv_sign_state_tr_pooled.df <- rbind(poxv_sign_state_tr_pooled.df, poxv_sign_state_tr_og.df)
    
    ### Record all plots
    poxv_sign_state_tr_og.p <- poxv_sign_state_tr_og.p + geom_text2(aes(subset=((transition %in% c("0->1", "1->0")) & (node %in% poxv_noc.v)), 
                                                                        label=ifelse(!is.na(noc_name) & !is.na(transition), paste(noc_name, transition, sep=":"), "")), size=3, hjust=1.3)
    poxv_sign_state_tr.p_list[[poxv_og_sign_name.v[i]]][[paste(poxv_og_sign_name.v[i] ,"_" , a, sep="")]] <- poxv_sign_state_tr_og.p
  }
}


##### Plot counts of domain gain and loss events associated with ortholog group gain and loss events identified by MPPA, MAP, and DOWNPASS methods
poxv_sign_state_tr_pooled.df <- poxv_sign_state_tr_pooled.df %>% mutate(sign_name = str_remove(og_name, "^N0HOG\\d{7}\\^"), og_name=str_split(og_name, "\\^", simplify=TRUE)[,1])
poxv_sign_state_tr_pooled.df <- poxv_sign_state_tr_pooled.df %>% mutate(sign_database = str_split(sign_name, "\\*", simplify = TRUE)[,1], sign_id = str_split(sign_name, "\\*", simplify = TRUE)[,2])

poxv_sign_no_gain_loss.p_list <- list()

for (a in c("MPPA","MAP","DOWNPASS")) {
  poxv_sign_state_tr_pooled_ct.df <- poxv_sign_state_tr_pooled.df %>% filter(parent!=node, transition %in% c("0->1", "1->0"), method==a, sign_database=="Pfam")
  poxv_sign_state_tr_pooled_ct.df <- poxv_sign_state_tr_pooled_ct.df %>% left_join(poxv_state_tr_all_m.df %>% select(parent, node, transition, og_name) %>% mutate(by_gene_transition = TRUE), by = c("parent", "node", "transition", "og_name")) %>% 
    mutate(by_gene_transition = if_else(is.na(by_gene_transition), FALSE, TRUE)) %>% filter(by_gene_transition==TRUE)
  poxv_sign_state_tr_pooled_ct.df <- poxv_sign_state_tr_pooled_ct.df %>% group_by(parent, node, transition) %>% summarise(transition_count = n()) %>% ungroup()

  ##### Number of domain gain events associated with ortholog group gain events identified by MPPA, MAP, and DOWNPASS methods
  poxv_sign_state_tr_pooled_gain_ct.info.df <- poxv_t.rt.info.df %>% left_join(poxv_sign_state_tr_pooled_ct.df %>% filter(transition == "0->1") %>% select(parent, node, transition_count), by = c("parent", "node"))
  poxv_sign_state_tr_pooled_gain_ct.info.df <- poxv_sign_state_tr_pooled_gain_ct.info.df %>% mutate(transition_count=ifelse(is.na(transition_count), 0, transition_count))
  
  poxv_sign_state_tr_pooled_gain_ct.p <- ggtree(poxv_t.rt.rn, linewidth=2, aes(color=transition_count))
  poxv_sign_state_tr_pooled_gain_ct.p <- poxv_sign_state_tr_pooled_gain_ct.p %<+% poxv_sign_state_tr_pooled_gain_ct.info.df
  poxv_sign_state_tr_pooled_gain_ct.p <- poxv_sign_state_tr_pooled_gain_ct.p + theme_tree()
  # poxv_sign_state_tr_pooled_gain_ct.p <- poxv_sign_state_tr_pooled_gain_ct.p + geom_text2(aes(label=transition_count), color="black", size=5, hjust=1.3)
  poxv_sign_state_tr_pooled_gain_ct.p <- poxv_sign_state_tr_pooled_gain_ct.p + scale_color_gradient("transition_count", trans="pseudo_log", name="Number of protein domain gain events", low="black", high="yellow", breaks=c(0,3,5,10,30,50,100,max(poxv_sign_state_tr_pooled_gain_ct.info.df$transition_count)))
  poxv_sign_state_tr_pooled_gain_ct.p <- poxv_sign_state_tr_pooled_gain_ct.p + labs(title=a)
  poxv_sign_state_tr_pooled_gain_ct.p <- poxv_sign_state_tr_pooled_gain_ct.p + theme(text = element_text(size = 24), plot.title = element_text(size = 30), plot.subtitle = element_text(size = 16))
  poxv_sign_state_tr_pooled_gain_ct.p <- poxv_sign_state_tr_pooled_gain_ct.p + xlim(0,4.5)
  
  ##### Number of domain loss events associated with ortholog group loss events identified by MPPA, MAP, and DOWNPASS methods
  poxv_sign_state_tr_pooled_loss_ct.info.df <- poxv_t.rt.info.df %>% left_join(poxv_sign_state_tr_pooled_ct.df %>% filter(transition == "1->0") %>% select(parent, node, transition_count), by = c("parent", "node"))
  poxv_sign_state_tr_pooled_loss_ct.info.df <- poxv_sign_state_tr_pooled_loss_ct.info.df %>% mutate(transition_count=ifelse(is.na(transition_count), 0, transition_count))
  
  poxv_sign_state_tr_pooled_loss_ct.p <- ggtree(poxv_t.rt.rn, linewidth=2, aes(color=transition_count))
  poxv_sign_state_tr_pooled_loss_ct.p <- poxv_sign_state_tr_pooled_loss_ct.p %<+% poxv_sign_state_tr_pooled_loss_ct.info.df
  poxv_sign_state_tr_pooled_loss_ct.p <- poxv_sign_state_tr_pooled_loss_ct.p + theme_tree()
  # poxv_sign_state_tr_pooled_loss_ct.p <- poxv_sign_state_tr_pooled_loss_ct.p + geom_text2(aes(label=transition_count), color="black", size=5, hjust=1.3)
  poxv_sign_state_tr_pooled_loss_ct.p <- poxv_sign_state_tr_pooled_loss_ct.p + scale_color_gradient("transition_count", trans="pseudo_log", name="Number of protein domain loss events", low="black", high="deepskyblue", breaks=c(0,3,5,10,30,50,100,max(poxv_sign_state_tr_pooled_loss_ct.info.df$transition_count)))
  poxv_sign_state_tr_pooled_loss_ct.p <- poxv_sign_state_tr_pooled_loss_ct.p + labs(title=a)
  poxv_sign_state_tr_pooled_loss_ct.p <- poxv_sign_state_tr_pooled_loss_ct.p + theme(text = element_text(size = 24), plot.title = element_text(size = 30), plot.subtitle = element_text(size = 16))
  poxv_sign_state_tr_pooled_loss_ct.p <- poxv_sign_state_tr_pooled_loss_ct.p + xlim(0,4.5)
  
  poxv_sign_no_gain_loss.p_list[[a]] <- poxv_sign_state_tr_pooled_gain_ct.p + poxv_sign_state_tr_pooled_loss_ct.p + plot_layout(ncol=1)
}

# pdf("Figures/domain_gain_and_loss_all_methods.pdf", width=24, height=18)
# wrap_plots(poxv_sign_no_gain_loss.p_list) + plot_layout(nrow=1)
# dev.off()


##### Plot counts of domain gain and loss events associated with ortholog group gain and loss events correspondingly identified by all three methods
poxv_sign_state_tr_all_m.df <- poxv_sign_state_tr_pooled.df %>% filter(parent!=node, transition %in% c("0->1", "1->0")) %>% group_by(parent, node, transition, og_name, sign_database, sign_id) %>% summarise(transition_count = n()) %>% ungroup()
poxv_sign_state_tr_all_m.df <- poxv_sign_state_tr_all_m.df %>% mutate(transition_count = transition_count/3) %>% filter(transition_count == 1) # Pick only transition events estimated correspondingly by all methods
poxv_sign_state_tr_all_m.df <- poxv_sign_state_tr_all_m.df %>% mutate(poxv_sign_name=paste(og_name, '^', sign_database, '*', sign_id, sep=""))
poxv_sign_state_tr_all_m.df <- poxv_sign_state_tr_all_m.df %>% left_join(poxv_state_tr_all_m.df %>% select(parent, node, transition, og_name) %>% mutate(by_gene_transition = TRUE), by = c("parent", "node", "transition", "og_name")) %>% 
  mutate(by_gene_transition = if_else(is.na(by_gene_transition), FALSE, TRUE)) %>% filter(by_gene_transition==TRUE)

poxv_sign_state_tr_all_m_ct_by_db.df <- poxv_sign_state_tr_all_m.df %>% group_by(parent, node, transition, sign_database) %>% summarise(transition_count = n()) %>% ungroup()

poxv_sign_state_gain_and_loss_by_db.p_list <- list()
poxv_sign_state_db.v <- sort(unique(poxv_sign_state_tr_all_m_ct_by_db.df$sign_database))

for (c in 1:length(poxv_sign_state_db.v)) {
  ########## ***** Figure 5A: Number of protein domain gain events associated with ortholog group gain events correspondingly identified by all methods *****
  poxv_sign_state_tr_all_m_gain_ct_by_db.info.df <- poxv_t.rt.info.df %>% left_join(poxv_sign_state_tr_all_m_ct_by_db.df %>% filter(transition == "0->1", sign_database == poxv_sign_state_db.v[c]) %>% select(parent, node, transition_count), by = c("parent", "node"))
  poxv_sign_state_tr_all_m_gain_ct_by_db.info.df <- poxv_sign_state_tr_all_m_gain_ct_by_db.info.df %>% mutate(transition_count=ifelse(is.na(transition_count), 0, transition_count))
  
  poxv_sign_state_tr_all_m_gain_ct_by_db.p <- ggtree(poxv_t.rt.rn, size=2, aes(color=transition_count))
  poxv_sign_state_tr_all_m_gain_ct_by_db.p <- poxv_sign_state_tr_all_m_gain_ct_by_db.p %<+% poxv_sign_state_tr_all_m_gain_ct_by_db.info.df
  poxv_sign_state_tr_all_m_gain_ct_by_db.p <- poxv_sign_state_tr_all_m_gain_ct_by_db.p + theme_tree()
  # poxv_sign_state_tr_all_m_gain_ct_by_db.p <- poxv_sign_state_tr_all_m_gain_ct_by_db.p + geom_text2(aes(label=transition_count), color="black", size=5, hjust=1.3)
  poxv_sign_state_tr_all_m_gain_ct_by_db.p <- poxv_sign_state_tr_all_m_gain_ct_by_db.p + scale_color_gradient("transition_count", trans="pseudo_log", name="Number of ortholog group gain", low="black", high="yellow", breaks=c(0,3,5,10,30,50,100,max(poxv_sign_state_tr_all_m_gain_ct_by_db.info.df$transition_count)))
  poxv_sign_state_tr_all_m_gain_ct_by_db.p <- poxv_sign_state_tr_all_m_gain_ct_by_db.p + labs(title="All methods: Domain gain", subtitle=poxv_sign_state_db.v[c])
  poxv_sign_state_tr_all_m_gain_ct_by_db.p <- poxv_sign_state_tr_all_m_gain_ct_by_db.p + theme(text = element_text(size = 20), plot.title = element_text(size = 30), plot.subtitle = element_text(size = 24))
  poxv_sign_state_tr_all_m_gain_ct_by_db.p <- poxv_sign_state_tr_all_m_gain_ct_by_db.p + xlim(0,4.5)

  ########## ***** Supplementary Figure 5: Number of protein domain loss events associated with ortholog group loss events correspondingly identified by all methods
  poxv_sign_state_tr_all_m_loss_ct_by_db.info.df <- poxv_t.rt.info.df %>% left_join(poxv_sign_state_tr_all_m_ct_by_db.df %>% filter(transition == "1->0", sign_database == poxv_sign_state_db.v[c]) %>% select(parent, node, transition_count), by = c("parent", "node"))
  poxv_sign_state_tr_all_m_loss_ct_by_db.info.df <- poxv_sign_state_tr_all_m_loss_ct_by_db.info.df %>% mutate(transition_count=ifelse(is.na(transition_count), 0, transition_count))
  
  poxv_sign_state_tr_all_m_loss_ct_by_db.p <- ggtree(poxv_t.rt.rn, size=2, aes(color=transition_count))
  poxv_sign_state_tr_all_m_loss_ct_by_db.p <- poxv_sign_state_tr_all_m_loss_ct_by_db.p %<+% poxv_sign_state_tr_all_m_loss_ct_by_db.info.df
  poxv_sign_state_tr_all_m_loss_ct_by_db.p <- poxv_sign_state_tr_all_m_loss_ct_by_db.p + theme_tree()
  # poxv_sign_state_tr_all_m_loss_ct_by_db.p <- poxv_sign_state_tr_all_m_loss_ct_by_db.p + geom_text2(aes(label=transition_count), color="black", size=5, hjust=1.3)
  poxv_sign_state_tr_all_m_loss_ct_by_db.p <- poxv_sign_state_tr_all_m_loss_ct_by_db.p + scale_color_gradient("transition_count", trans="pseudo_log", name="Number of ortholog group gain", low="black", high="deepskyblue", breaks=c(0,3,5,10,30,50,100,max(poxv_sign_state_tr_all_m_gain_ct_by_db.info.df$transition_count)))
  poxv_sign_state_tr_all_m_loss_ct_by_db.p <- poxv_sign_state_tr_all_m_loss_ct_by_db.p + labs(title="All methods: Domain loss", subtitle=poxv_sign_state_db.v[c])
  poxv_sign_state_tr_all_m_loss_ct_by_db.p <- poxv_sign_state_tr_all_m_loss_ct_by_db.p + theme(text = element_text(size = 20), plot.title = element_text(size = 30), plot.subtitle = element_text(size = 24))
  poxv_sign_state_tr_all_m_loss_ct_by_db.p <- poxv_sign_state_tr_all_m_loss_ct_by_db.p + xlim(0,4.5)

  poxv_sign_state_gain_and_loss_by_db.p_list[[poxv_sign_state_db.v[c]]] <- poxv_sign_state_tr_all_m_gain_ct_by_db.p + poxv_sign_state_tr_all_m_loss_ct_by_db.p + plot_layout(nrow = 1)
}

# pdf("Figures/domain_gain_and_loss_summarized.pdf", width=24, height=10)
# wrap_plots(poxv_sign_state_gain_and_loss_by_db.p_list[c("Pfam")])
# dev.off()



########## Analyze enrichment of protein domains in gained and lost ortholog groups
##### Count number of gain and loss events of each protein domain occurred in the poxvirus lineage
### Create Pfam/Interpro ID mapping table
interpro_res.f <- "251208_pooled_translated_cds_renamed_interproscan_result.tsv"
interpro_res.df <- read.delim(interpro_res.f, check.names = FALSE, header = F)

interpro_res_dict.df <- interpro_res.df %>% select(V4, V5, V6, V12, V13)
colnames(interpro_res_dict.df) <- c("analysis","sign_acc_id","sign_desc","interpro_annot","interpro_desc")

interpro_res_dict.df <- interpro_res_dict.df %>% filter(analysis=="Pfam") %>% arrange(sign_acc_id) %>% distinct()

### Note: We didn't analyze the exact numbers of gained and lost protein domains (only states of presence/absence of protein domains in gained and lost ortholog groups)
poxv_sign_state_tr_all_m_sign_ct.df <- poxv_sign_state_tr_all_m.df %>% filter(sign_database == "Pfam") %>% group_by(transition, sign_id) %>% summarise(transition_count = n(), .groups = "drop")

poxv_sign_state_tr_all_m_ct.df <- poxv_sign_state_tr_all_m_sign_ct.df %>% group_by(sign_id) %>% summarise(sum_transition_count = sum(transition_count)) %>% ungroup() %>% arrange(desc(sum_transition_count))

poxv_sign_state_tr_all_m_sign_ct.df <- poxv_sign_state_tr_all_m_sign_ct.df %>% select(sign_id, transition, transition_count) %>% pivot_wider(names_from  = transition, values_from = transition_count, values_fill = 0)
colnames(poxv_sign_state_tr_all_m_sign_ct.df) <- c("sign_id","gain_num","loss_num")

poxv_sign_state_tr_all_m_sign_ct.df <- poxv_sign_state_tr_all_m_sign_ct.df %>% mutate(sign_desc=interpro_res_dict.df$sign_desc[match(sign_id, interpro_res_dict.df$sign_acc_id)])

top_gain <- poxv_sign_state_tr_all_m_sign_ct.df %>% arrange(desc(gain_num)) %>% top_n(3, gain_num)
top_loss <- poxv_sign_state_tr_all_m_sign_ct.df %>% arrange(desc(loss_num)) %>% top_n(3, loss_num)
top_union <- distinct(bind_rows(top_gain, top_loss))

########## ***** Figure 5B: Number of gain and loss events of each protein domain occurred in the poxvirus lineage *****
poxv_sign_state_tr_all_m_sign_ct.p <- ggplot(poxv_sign_state_tr_all_m_sign_ct.df, aes(x=gain_num, y=loss_num))
poxv_sign_state_tr_all_m_sign_ct.p <- poxv_sign_state_tr_all_m_sign_ct.p + geom_point()
# poxv_sign_state_tr_all_m_sign_ct.p <- poxv_sign_state_tr_all_m_sign_ct.p + geom_jitter(alpha=0.1, shape=16, size=3, position=position_jitter(width=0.2, height=0.2))
poxv_sign_state_tr_all_m_sign_ct.p <- poxv_sign_state_tr_all_m_sign_ct.p + geom_text(data=top_union, aes(label=sign_desc), hjust=-0.13, vjust=0, size=7)
poxv_sign_state_tr_all_m_sign_ct.p <- poxv_sign_state_tr_all_m_sign_ct.p + labs(title="Domain gain and loss from poxvirus genes", x="Number of protein domain gain events", y="Number of protein domain loss events")
poxv_sign_state_tr_all_m_sign_ct.p <- poxv_sign_state_tr_all_m_sign_ct.p + theme_classic()
poxv_sign_state_tr_all_m_sign_ct.p <- poxv_sign_state_tr_all_m_sign_ct.p + theme(rect=element_rect(fill="transparent"), text=element_text(size = 24), plot.title=element_text(size = 32))
poxv_sign_state_tr_all_m_sign_ct.p <- poxv_sign_state_tr_all_m_sign_ct.p + scale_x_continuous(breaks=seq(0, 20, by=5), limits=c(0,20))
poxv_sign_state_tr_all_m_sign_ct.p <- poxv_sign_state_tr_all_m_sign_ct.p + scale_y_continuous(breaks=seq(0, 20, by=5), limits=c(0,20))
poxv_sign_state_tr_all_m_sign_ct.p <- poxv_sign_state_tr_all_m_sign_ct.p + coord_fixed()
poxv_sign_state_tr_all_m_sign_ct.p

# pdf("Figures/domain_gain_and_loss_num_events.pdf", width=10, height=10)
# print(poxv_sign_state_tr_all_m_sign_ct.p)
# dev.off()


##### Count number of gain or loss events of each protein domain occurred in the poxvirus lineage
poxv_sign_state_tr_all_m_sign_sep_ct.df <- poxv_sign_state_tr_all_m.df %>% filter(sign_database == "Pfam") %>% group_by(transition, sign_id) %>% summarise(transition_count = n()) %>% ungroup()
poxv_sign_state_tr_all_m_sign_sep_ct.df <- poxv_sign_state_tr_all_m_sign_sep_ct.df %>% arrange(transition, desc(transition_count))
poxv_sign_state_tr_all_m_sign_sep_ct.df <- poxv_sign_state_tr_all_m_sign_sep_ct.df %>% mutate(sign_desc=interpro_res_dict.df$sign_desc[match(sign_id, interpro_res_dict.df$sign_acc_id)])

########## ***** Supplementary Table 4: Numbers of gain and loss events of Pfam protein domains accompanied by ortholog group gain and loss *****
poxv_sign_state_tr_all_m_sign_sep_ct.f <- "Tables/251208_domain_gain_and_loss_num_separated_by_id.tsv"
# write.table(poxv_sign_state_tr_all_m_sign_sep_ct.df, file = poxv_sign_state_tr_all_m_sign_sep_ct.f, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)



########## Enrichment analysis of protein domain gained or lost in the poxvirus lineage
poxv_sign_state_tr_all_m_enrich.df <- poxv_sign_state_tr_pooled.df %>% filter(parent!=node) %>% group_by(parent, node, transition, og_name, sign_database, sign_id) %>% summarise(transition_count = n()) %>% ungroup()
poxv_sign_state_tr_all_m_enrich.df <- poxv_sign_state_tr_all_m_enrich.df %>% mutate(transition_count = transition_count/3) %>% filter(transition_count == 1) # Pick only transition events estimated correspondingly by all methods
poxv_sign_state_tr_all_m_enrich.df <- poxv_sign_state_tr_all_m_enrich.df %>% left_join(poxv_state_tr_all_m.df %>% select(parent, node, transition, og_name) %>% mutate(by_gene_transition = TRUE), by = c("parent", "node", "transition", "og_name")) %>% 
  mutate(by_gene_transition = if_else(is.na(by_gene_transition), FALSE, TRUE)) %>% filter(by_gene_transition==TRUE)

poxv_sign_state_tr_all_m_pfam_enrich.df <- poxv_sign_state_tr_all_m_enrich.df %>% filter(sign_database=="Pfam")
poxv_sign_state_tr_all_m_pfam_enrich.tab <- poxv_sign_state_tr_all_m_pfam_enrich.df %>% group_by(sign_id, transition) %>% summarise(count = n(), .groups = "drop") %>% pivot_wider(names_from = transition, values_from = count, values_fill = list(count = 0))

##### Perform enrichment analysis for domain gained in or lost from the poxvirus lineage
poxv_sign_state_tr_all_m_pfam_enrich_gain_or_loss.tab <- poxv_sign_state_tr_all_m_pfam_enrich.tab
poxv_sign_state_tr_all_m_pfam_enrich_gain_or_loss.tab$gained_or_lost <- poxv_sign_state_tr_all_m_pfam_enrich_gain_or_loss.tab$"0->1" + poxv_sign_state_tr_all_m_pfam_enrich_gain_or_loss.tab$"1->0"
poxv_sign_state_tr_all_m_pfam_enrich_gain_or_loss.tab$unchanged <- poxv_sign_state_tr_all_m_pfam_enrich_gain_or_loss.tab$"0->0" + poxv_sign_state_tr_all_m_pfam_enrich_gain_or_loss.tab$"1->1"
poxv_sign_state_tr_all_m_pfam_enrich_gain_or_loss.tab <- poxv_sign_state_tr_all_m_pfam_enrich_gain_or_loss.tab[,c("sign_id","gained_or_lost","unchanged")] 

poxv_sign_state_tr_all_m_pfam_enrich_gain_or_loss_p_values <- sapply(1:nrow(poxv_sign_state_tr_all_m_pfam_enrich_gain_or_loss.tab), function(i) { # Perform Fisher's Exact Test for each sign_id
  poxv_sign_state_tr_all_m_pfam_enrich_gain_or_loss_sign.tab <- matrix(c(poxv_sign_state_tr_all_m_pfam_enrich_gain_or_loss.tab$gained_or_lost[i], poxv_sign_state_tr_all_m_pfam_enrich_gain_or_loss.tab$unchanged[i], 
                                                                 sum(poxv_sign_state_tr_all_m_pfam_enrich_gain_or_loss.tab$gained_or_lost) - poxv_sign_state_tr_all_m_pfam_enrich_gain_or_loss.tab$gained_or_lost[i], 
                                                                 sum(poxv_sign_state_tr_all_m_pfam_enrich_gain_or_loss.tab$unchanged) - poxv_sign_state_tr_all_m_pfam_enrich_gain_or_loss.tab$unchanged[i]), nrow = 2)
  fisher.test(poxv_sign_state_tr_all_m_pfam_enrich_gain_or_loss_sign.tab)$p.value
})

poxv_sign_state_tr_all_m_pfam_enrich_gain_or_loss_adj_p_values <- p.adjust(poxv_sign_state_tr_all_m_pfam_enrich_gain_or_loss_p_values, method = "fdr")

poxv_sign_state_tr_all_m_pfam_enrich_gain_or_loss_odd_ratio <- sapply(1:nrow(poxv_sign_state_tr_all_m_pfam_enrich_gain_or_loss.tab), function(i) { # Perform Fisher's Exact Test for each sign_id
  poxv_sign_state_tr_all_m_pfam_enrich_gain_or_loss_sign.tab <- matrix(c(poxv_sign_state_tr_all_m_pfam_enrich_gain_or_loss.tab$gained_or_lost[i], poxv_sign_state_tr_all_m_pfam_enrich_gain_or_loss.tab$unchanged[i], 
                                                                         sum(poxv_sign_state_tr_all_m_pfam_enrich_gain_or_loss.tab$gained_or_lost) - poxv_sign_state_tr_all_m_pfam_enrich_gain_or_loss.tab$gained_or_lost[i], 
                                                                         sum(poxv_sign_state_tr_all_m_pfam_enrich_gain_or_loss.tab$unchanged) - poxv_sign_state_tr_all_m_pfam_enrich_gain_or_loss.tab$unchanged[i]), nrow = 2)
  fisher.test(poxv_sign_state_tr_all_m_pfam_enrich_gain_or_loss_sign.tab)$estimate
})

poxv_sign_state_tr_all_m_pfam_enrich_gain_or_loss.tab <- poxv_sign_state_tr_all_m_pfam_enrich_gain_or_loss.tab %>%
  mutate(p_value = poxv_sign_state_tr_all_m_pfam_enrich_gain_or_loss_p_values, adj_p_value = poxv_sign_state_tr_all_m_pfam_enrich_gain_or_loss_adj_p_values, odd_ratio = poxv_sign_state_tr_all_m_pfam_enrich_gain_or_loss_odd_ratio) %>% arrange(adj_p_value, p_value)
poxv_sign_state_tr_all_m_pfam_enrich_gain_or_loss.tab <- poxv_sign_state_tr_all_m_pfam_enrich_gain_or_loss.tab %>% mutate(sign_desc=interpro_res_dict.df$sign_desc[match(sign_id, interpro_res_dict.df$sign_acc_id)])

head(poxv_sign_state_tr_all_m_pfam_enrich_gain_or_loss.tab)

########## ***** Supplementary Table 5. Enrichment analysis of protein domains from ortholog groups gained or lost during poxvirus evolution *****
poxv_sign_state_tr_all_m_pfam_enrich_gain_or_loss.f <- "Tables/251208_domain_gain_or_loss_enrichment_analysis.tsv"
# write.table(poxv_sign_state_tr_all_m_pfam_enrich_gain_or_loss.tab, file = poxv_sign_state_tr_all_m_pfam_enrich_gain_or_loss.f, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


##### Perform enrichment analysis for domain gained in the poxvirus lineage
poxv_sign_state_tr_all_m_pfam_enrich_gain.tab <- poxv_sign_state_tr_all_m_pfam_enrich.tab[,c("sign_id","0->1","0->0")] 
colnames(poxv_sign_state_tr_all_m_pfam_enrich_gain.tab)[2:3] <- c("gained", "non-gained")

poxv_sign_state_tr_all_m_pfam_enrich_gain_p_values <- sapply(1:nrow(poxv_sign_state_tr_all_m_pfam_enrich_gain.tab), function(i) { # Perform Fisher's Exact Test for each sign_id
  poxv_sign_state_tr_all_m_pfam_enrich_gain_sign.tab <- matrix(c(poxv_sign_state_tr_all_m_pfam_enrich_gain.tab$gained[i], poxv_sign_state_tr_all_m_pfam_enrich_gain.tab$"non-gained"[i], 
                                                                       sum(poxv_sign_state_tr_all_m_pfam_enrich_gain.tab$gained) - poxv_sign_state_tr_all_m_pfam_enrich_gain.tab$gained[i], 
                                                                       sum(poxv_sign_state_tr_all_m_pfam_enrich_gain.tab$"non-gained") - poxv_sign_state_tr_all_m_pfam_enrich_gain.tab$"non-gained"[i]), nrow = 2)
  fisher.test(poxv_sign_state_tr_all_m_pfam_enrich_gain_sign.tab)$p.value
})

poxv_sign_state_tr_all_m_pfam_enrich_gain_adj_p_values <- p.adjust(poxv_sign_state_tr_all_m_pfam_enrich_gain_p_values, method = "fdr")

poxv_sign_state_tr_all_m_pfam_enrich_gain_odd_ratio <- sapply(1:nrow(poxv_sign_state_tr_all_m_pfam_enrich_gain.tab), function(i) { # Perform Fisher's Exact Test for each sign_id
  poxv_sign_state_tr_all_m_pfam_enrich_gain_sign.tab <- matrix(c(poxv_sign_state_tr_all_m_pfam_enrich_gain.tab$gained[i], poxv_sign_state_tr_all_m_pfam_enrich_gain.tab$"non-gained"[i], 
                                                                 sum(poxv_sign_state_tr_all_m_pfam_enrich_gain.tab$gained) - poxv_sign_state_tr_all_m_pfam_enrich_gain.tab$gained[i], 
                                                                 sum(poxv_sign_state_tr_all_m_pfam_enrich_gain.tab$"non-gained") - poxv_sign_state_tr_all_m_pfam_enrich_gain.tab$"non-gained"[i]), nrow = 2)
  fisher.test(poxv_sign_state_tr_all_m_pfam_enrich_gain_sign.tab)$estimate
})

poxv_sign_state_tr_all_m_pfam_enrich_gain.tab <- poxv_sign_state_tr_all_m_pfam_enrich_gain.tab %>%
  mutate(p_value = poxv_sign_state_tr_all_m_pfam_enrich_gain_p_values, adj_p_value = poxv_sign_state_tr_all_m_pfam_enrich_gain_adj_p_values, odd_ratio = poxv_sign_state_tr_all_m_pfam_enrich_gain_odd_ratio) %>% arrange(adj_p_value, p_value)
poxv_sign_state_tr_all_m_pfam_enrich_gain.tab <- poxv_sign_state_tr_all_m_pfam_enrich_gain.tab %>% mutate(sign_desc=interpro_res_dict.df$sign_desc[match(sign_id, interpro_res_dict.df$sign_acc_id)])

head(poxv_sign_state_tr_all_m_pfam_enrich_gain.tab)


##### Perform enrichment analysis for domain lost from the poxvirus lineage
poxv_sign_state_tr_all_m_pfam_enrich_loss.tab <- poxv_sign_state_tr_all_m_pfam_enrich.tab[,c("sign_id","1->0","1->1")] 
colnames(poxv_sign_state_tr_all_m_pfam_enrich_loss.tab)[2:3] <- c("lost", "non-lost")

poxv_sign_state_tr_all_m_pfam_enrich_loss_p_values <- sapply(1:nrow(poxv_sign_state_tr_all_m_pfam_enrich_loss.tab), function(i) { # Perform Fisher's Exact Test for each sign_id
  poxv_sign_state_tr_all_m_pfam_enrich_loss_sign.tab <- matrix(c(poxv_sign_state_tr_all_m_pfam_enrich_loss.tab$lost[i], poxv_sign_state_tr_all_m_pfam_enrich_loss.tab$"non-lost"[i], 
                                                                      sum(poxv_sign_state_tr_all_m_pfam_enrich_loss.tab$lost) - poxv_sign_state_tr_all_m_pfam_enrich_loss.tab$lost[i], 
                                                                      sum(poxv_sign_state_tr_all_m_pfam_enrich_loss.tab$"non-lost") - poxv_sign_state_tr_all_m_pfam_enrich_loss.tab$"non-lost"[i]), nrow = 2)
  fisher.test(poxv_sign_state_tr_all_m_pfam_enrich_loss_sign.tab)$p.value
})

poxv_sign_state_tr_all_m_pfam_enrich_loss_adj_p_values <- p.adjust(poxv_sign_state_tr_all_m_pfam_enrich_loss_p_values, method = "fdr")

poxv_sign_state_tr_all_m_pfam_enrich_loss_odd_ratio <- sapply(1:nrow(poxv_sign_state_tr_all_m_pfam_enrich_loss.tab), function(i) { # Perform Fisher's Exact Test for each sign_id
  poxv_sign_state_tr_all_m_pfam_enrich_loss_sign.tab <- matrix(c(poxv_sign_state_tr_all_m_pfam_enrich_loss.tab$lost[i], poxv_sign_state_tr_all_m_pfam_enrich_loss.tab$"non-lost"[i], 
                                                                 sum(poxv_sign_state_tr_all_m_pfam_enrich_loss.tab$lost) - poxv_sign_state_tr_all_m_pfam_enrich_loss.tab$lost[i], 
                                                                 sum(poxv_sign_state_tr_all_m_pfam_enrich_loss.tab$"non-lost") - poxv_sign_state_tr_all_m_pfam_enrich_loss.tab$"non-lost"[i]), nrow = 2)
  fisher.test(poxv_sign_state_tr_all_m_pfam_enrich_loss_sign.tab)$estimate
})

poxv_sign_state_tr_all_m_pfam_enrich_loss.tab <- poxv_sign_state_tr_all_m_pfam_enrich_loss.tab %>%
  mutate(p_value = poxv_sign_state_tr_all_m_pfam_enrich_loss_p_values, adj_p_value = poxv_sign_state_tr_all_m_pfam_enrich_loss_adj_p_values, odd_ratio = poxv_sign_state_tr_all_m_pfam_enrich_loss_odd_ratio) %>% arrange(adj_p_value, p_value)
poxv_sign_state_tr_all_m_pfam_enrich_loss.tab <- poxv_sign_state_tr_all_m_pfam_enrich_loss.tab %>% mutate(sign_desc=interpro_res_dict.df$sign_desc[match(sign_id, interpro_res_dict.df$sign_acc_id)])

head(poxv_sign_state_tr_all_m_pfam_enrich_loss.tab)


# ##### Assign status of gain and loss to each protein domain
# poxv_sign_state_root_status.df <- as.data.frame(poxv_sign_state_tr_pooled.df %>% filter(parent==node) %>% group_by(parent, node, state_parent, og_name, sign_database, sign_id) %>% summarise(state_count = n()) %>% ungroup() %>% mutate(state_count = state_count/3) %>% filter(state_count == 1))
# 
# poxv_sign_state_absent_tip.df <- poxv_sign_state_tr_pooled.df %>% group_by(parent, node, transition, og_name, sign_database, sign_id) %>% summarise(transition_count = n()) %>% ungroup()
# poxv_sign_state_absent_tip.df <- poxv_sign_state_absent_tip.df %>% mutate(transition_count = transition_count/3) %>% filter(transition_count == 1) # Pick only transition events estimated correspondingly by all methods
# 
# poxv_sign_state_absent_tip_ct.df <- poxv_sign_state_tr_pooled.df %>% filter(isTip==TRUE, method=="MPPA", state_desc == "Present")
# poxv_sign_state_absent_tip_ct.df <- poxv_sign_state_absent_tip_ct.df %>% group_by(og_name, sign_database, sign_id) %>% summarise(tip_set = list(unique(node)), num_tip = length(tip_set[[1]]), .groups = "drop")
# poxv_sign_state_absent_tip_ct.df <- poxv_sign_state_absent_tip_ct.df %>% mutate(ca_node = map_int(tip_set, ~ (if (length(.x) < 2) {.x[[1]]} else {getMRCA(poxv_t.rt.rn, .x)})),
#                                                                       num_desc = map_int(ca_node, ~ length(Descendants(poxv_t.rt.rn, .x, type = "tips")[[1]])))
# poxv_sign_state_absent_tip_ct.df <- poxv_sign_state_absent_tip_ct.df %>% left_join(poxv_sign_state_absent_tip.df %>% rename(ca_node = node) %>% select(og_name, sign_database, sign_id, ca_node, transition), by = c("og_name", "sign_database", "sign_id", "ca_node"))
# poxv_sign_state_absent_tip_ct.df <- poxv_sign_state_absent_tip_ct.df %>% mutate(ca_node_status = ifelse(!is.na(transition), str_split(transition, "->", simplify=TRUE)[,2], transition))
# poxv_sign_state_absent_tip_ct.df <- poxv_sign_state_absent_tip_ct.df %>% mutate(loss_status = ifelse(!is.na(ca_node_status) & ca_node_status==1, ifelse(num_tip < num_desc, "some_absent", "conserved_after_gain"), "undeterminable")) # Consider only the case that CA is 1
# poxv_sign_state_absent_tip_ct.df <- poxv_sign_state_absent_tip_ct.df %>% mutate(loss_status = ifelse(ca_node==find_root(poxv_t.rt.rn) & !is.na(ca_node_status) & ca_node_status==1 & num_desc==length(poxv_t.rt.rn$tip.label), "conserved", loss_status))
# 
# poxv_sign_gain_loss_status.df <- data.frame(poxv_sign_name=names(poxv_sign_state_tr.p_list))
# poxv_sign_gain_loss_status.df <- poxv_sign_gain_loss_status.df %>% mutate(og_name = str_split(poxv_sign_name, "\\^", simplify = TRUE)[,1], 
#                                                                           sign_database = str_split(str_split(poxv_sign_name, "\\^", simplify = TRUE)[,2], "\\*", simplify = TRUE)[,1], 
#                                                                           sign_id = str_split(str_split(poxv_sign_name, "\\^", simplify = TRUE)[,2], "\\*", simplify = TRUE)[,2])
# poxv_sign_gain_loss_status.df <- poxv_sign_gain_loss_status.df %>% left_join(poxv_sign_state_root_status.df %>% select(og_name, sign_database, sign_id, state_parent), by = c("og_name","sign_database","sign_id")) %>% rename(root_status=state_parent)
# poxv_sign_gain_loss_status.df <- poxv_sign_gain_loss_status.df %>% left_join(poxv_sign_state_absent_tip_ct.df %>% select(og_name, sign_database, sign_id, loss_status), by = c("og_name","sign_database","sign_id")) %>% rename(tip_status=loss_status)
# poxv_sign_gain_loss_status.df <- poxv_sign_gain_loss_status.df %>% mutate(detected_gain_event=ifelse(poxv_sign_name %in% (poxv_sign_state_tr_all_m.df %>% filter(transition=="0->1"))$poxv_sign_name, "yes", "no"),
#                                                                           detected_loss_event=ifelse(poxv_sign_name %in% (poxv_sign_state_tr_all_m.df %>% filter(transition=="1->0"))$poxv_sign_name, "yes", "no"))
# poxv_sign_gain_loss_status.df <- poxv_sign_gain_loss_status.df %>% mutate(gain_loss_status=ifelse(tip_status=="conserved", "conserved", "undeterminable"))
# poxv_sign_gain_loss_status.df <- poxv_sign_gain_loss_status.df %>% mutate(gain_loss_status=ifelse(detected_gain_event=="yes" | (!is.na(root_status) & root_status==0),
#                                                                                                   ifelse(detected_loss_event=="yes" | tip_status=="some_absent", "both", "gain"),
#                                                                                                   ifelse(detected_loss_event=="yes" | tip_status=="some_absent", "loss", gain_loss_status)))
# poxv_sign_gain_loss_status.df$gain_loss_status <- factor(poxv_sign_gain_loss_status.df$gain_loss_status, levels=c("gain","loss","both","conserved","undeterminable"))
# 
# table(poxv_sign_gain_loss_status.df$gain_loss_status)



############################## Selective pressure analysis ##############################
poxv_cdn_alm_ds <- "241008_Poxviridae_and_representative_MPXV_codon_alignment"
poxv_og_id.v <- str_split(list.files(paste(poxv_cdn_alm_ds, "/og_cds", sep=""), pattern="*.mafft.treefile$"), "_cds", simplify=TRUE)[,1]

# ############################## Comparing selective pressures across different type of genes ##############################
# ### Distribution of levels of selective pressure in poxviruses genes
# sel_st_rs_gene <- data.frame(og_name=poxv_og_id.v)
# # sel_st_rs_gene$prop_positive <- (as.data.frame(poxv_meme_st_ct.df))[match(sel_st_rs_gene$og_name, (as.data.frame(poxv_meme_st_ct.df))$og_name), "prop_positive"]
# sel_st_rs_gene$prop_positive <- (as.data.frame(poxv_fel_st_ct.df %>% filter(selection_class=="Diversifying")))[match(sel_st_rs_gene$og_name, (as.data.frame(poxv_fel_st_ct.df %>% filter(selection_class=="Diversifying")))$og_name), "prop_site"]
# sel_st_rs_gene$prop_negative <- (as.data.frame(poxv_fel_st_ct.df %>% filter(selection_class=="Purifying")))[match(sel_st_rs_gene$og_name, (as.data.frame(poxv_fel_st_ct.df %>% filter(selection_class=="Purifying")))$og_name), "prop_site"]
# sel_st_rs_gene$prop_neutral <- (as.data.frame(poxv_fel_st_ct.df %>% filter(selection_class=="Neutral")))[match(sel_st_rs_gene$og_name, (as.data.frame(poxv_fel_st_ct.df %>% filter(selection_class=="Neutral")))$og_name), "prop_site"]
# sel_st_rs_gene$prop_invariable <- (as.data.frame(poxv_fel_st_ct.df %>% filter(selection_class=="Invariable")))[match(sel_st_rs_gene$og_name, (as.data.frame(poxv_fel_st_ct.df %>% filter(selection_class=="Invariable")))$og_name), "prop_site"]
# sel_st_rs_gene[is.na(sel_st_rs_gene)] <- 0
# sel_st_rs_gene <- sel_st_rs_gene %>% mutate(vacv_gene_name=poxv_vacv_gene_name.df[match(og_name, poxv_vacv_gene_name.df$og_id), "Vaccinia_virus"])
# sel_st_rs_gene <- sel_st_rs_gene %>% mutate(og_name=gsub("\\.", "", og_name))
# sel_st_rs_gene <- sel_st_rs_gene %>% mutate(gain_loss_status=poxv_og_gain_loss_status.df[match(og_name, poxv_og_gain_loss_status.df$poxv_og_name), "gain_loss_status"])
# sel_st_rs_gene <- sel_st_rs_gene %>% mutate(gained_status=poxv_og_gain_loss_status.df[match(og_name, poxv_og_gain_loss_status.df$poxv_og_name), "gained_status"])
# 
# poxv_state_tr_all_m_ct_by_og.df <- as.data.frame(poxv_state_tr_all_m.df %>% filter(transition %in% c("0->1","1->0")) %>% group_by(og_name) %>% summarise(transition_count = n()) %>% ungroup())
# 
# sel_st_rs_gene <- sel_st_rs_gene %>% mutate(transition_count=ifelse(!is.na(poxv_state_tr_all_m_ct_by_og.df[match(og_name, poxv_state_tr_all_m_ct_by_og.df$og_name), "transition_count"]),
#                                             poxv_state_tr_all_m_ct_by_og.df[match(og_name, poxv_state_tr_all_m_ct_by_og.df$og_name), "transition_count"], 0), og_name=gsub("N0HOG", "N0.HOG", og_name))
# 

########## Calculate dN/dS using ape
poxv_og_dnds.df <- data.frame(HOG=character(0), dNdS=numeric(0), stringsAsFactors=FALSE)

for (a in 1:length(poxv_og_id.v)) {
  poxv_gene_pal2nal.f <- paste(poxv_cdn_alm_ds, "/og_translated_cds/", poxv_og_id.v[a], "_translated_cds.pal2nal", sep="")
  if (file.exists(poxv_gene_pal2nal.f) == FALSE) {next}
  
  # keep DNAbin for dnds()
  poxv_gene_pal2nal.aln <- read.alignment(poxv_gene_pal2nal.f, format = "fasta")

  res <- kaks(poxv_gene_pal2nal.aln, forceUpperCase=TRUE, rmgap=TRUE)

  ka <- as.vector(res$ka)
  ks <- as.vector(res$ks)
  
  valid <- is.finite(ka) & is.finite(ks) & ks > 0 & ks < 9.999 & ka >= 0
  # Note: when the alignment does not contain enough information (i.e. close to saturation), the Ka and Ks values are forced to 10 (more exactly to 9.999999).

  omega <- sum(ka[valid]) / sum(ks[valid])

  poxv_og_dnds.df <- rbind(poxv_og_dnds.df, data.frame(HOG=poxv_og_id.v[a], dNdS=omega, stringsAsFactors=FALSE))
}

poxv_og_dnds.df <- poxv_og_dnds.df %>% mutate(vacv_gene_name=poxv_vacv_gene_name.df[match(HOG, poxv_vacv_gene_name.df$og_id), "Vaccinia_virus"])
poxv_og_dnds.df <- poxv_og_dnds.df %>% mutate(HOG=gsub("\\.", "", HOG))
poxv_og_dnds.df <- poxv_og_dnds.df %>% mutate(gain_loss_status=poxv_og_gain_loss_status.df[match(HOG, poxv_og_gain_loss_status.df$poxv_og_name), "gain_loss_status"])
poxv_og_dnds.df <- poxv_og_dnds.df %>% mutate(gained_status=poxv_og_gain_loss_status.df[match(HOG, poxv_og_gain_loss_status.df$poxv_og_name), "gained_status"])

### Test difference in dN/dS
poxv_og_dnds.df %>% filter(!is.nan(dNdS), !gained_status=="undeterminable") %>% group_by(gained_status) %>% shapiro_test(dNdS)
poxv_og_dnds.df %>% filter(!is.nan(dNdS), !gained_status=="undeterminable") %>% levene_test(dNdS~gained_status)

### Use Kruskal-Wallis test since the normality assumption is violated
poxv_og_dnds_res.kruskal <- poxv_og_dnds.df %>% filter(!is.nan(dNdS), !gained_status=="undeterminable") %>% kruskal_test(dNdS~gained_status)
poxv_og_dnds_res.kruskal

poxv_og_dnds_stats <- poxv_og_dnds.df %>% filter(!is.nan(dNdS), !gained_status=="undeterminable") %>% dunn_test(dNdS~gained_status, p.adjust.method="fdr")
poxv_og_dnds_stats <- poxv_og_dnds_stats %>% add_xy_position(x="gained_status", fun="max")
poxv_og_dnds_stats <- poxv_og_dnds_stats %>% mutate(y.position=log10(y.position))
poxv_og_dnds_stats

########## ***** Figure 2D: dN/dS ratios of gained and pre-existing ortholog groups *****
poxv_og_dnds_box.p <- ggplot(poxv_og_dnds.df %>% filter(!is.nan(dNdS), !gained_status=="undeterminable"), aes(x=gained_status, y=log10(dNdS), col=gained_status), alpha=0.5)
poxv_og_dnds_box.p <- poxv_og_dnds_box.p + geom_hline(yintercept=0, linetype="dashed", color="black", size=0.5)
poxv_og_dnds_box.p <- poxv_og_dnds_box.p + geom_jitter(height=0, width=0.25, alpha=0.1, shape=16, size=3)
poxv_og_dnds_box.p <- poxv_og_dnds_box.p + geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.75)
poxv_og_dnds_box.p <- poxv_og_dnds_box.p + scale_color_manual("gained_status", values = c("gained"="orange","pre-existing"="navyblue","undeterminable"="lightgrey"))
poxv_og_dnds_box.p <- poxv_og_dnds_box.p + labs(title="dN/dS", x="Category", y="log10 dN/dS")
poxv_og_dnds_box.p <- poxv_og_dnds_box.p + theme_classic()
poxv_og_dnds_box.p <- poxv_og_dnds_box.p + theme(text = element_text(size = 18), plot.title = element_text(size = 24), plot.subtitle = element_text(size = 22))
poxv_og_dnds_box.p <- poxv_og_dnds_box.p + stat_pvalue_manual(poxv_og_dnds_stats, label="p.adj", hide.ns=TRUE, inherit.aes=FALSE, size=3)
poxv_og_dnds_box.p <- poxv_og_dnds_box.p + scale_y_continuous(breaks=c(-2,-1,0,1,2), limits=c(-2,2))
poxv_og_dnds_box.p

# pdf("Figures/dNdS_ratio_genes.pdf", width=5, height=8)
# poxv_og_dnds_box.p
# dev.off()



############################## Plot gain and loss events of protein domains of interests in each ortholog group harboring them ##############################
pfof_pastml_rs_sign.f <- paste("251208_gain_and_loss_domain_analysis/251208_cs_matrix_pooled_additional_domain_pastml_MPPA_output.tsv", sep="_")
pfof_pastml_rs_sign <- fread(pfof_pastml_rs_sign.f, header=T, sep="\t", quote="", check.names=F, data.table=FALSE)

pfof_pastml_rs_sign_state.df <- pfof_pastml_rs_sign
pfof_pastml_rs_sign_state.df[which(pfof_pastml_rs_sign_state.df$node=="Root"),"node"] <- find_root(poxv_t.rt.rn)
pfof_pastml_rs_sign_state.df <- pfof_pastml_rs_sign_state.df %>% mutate(node=ifelse(node=="GCF_000857045.1_Monkeypox_virus", "GCF_000857045.1_Monkeypox_virus_Zaire-96-I-16_Ia",
                                                                             ifelse(node=="GCA_039269415.1_Monkeypox_virus", "GCA_039269415.1_Monkeypox_virus_RDC-NKV-GOM-MPOX-010_Ib",
                                                                             ifelse(node=="GCA_006465585.1_Monkeypox_virus", "GCA_006465585.1_Monkeypox_virus_Sierra_Leone_IIa",
                                                                             ifelse(node=="GCF_014621545.1_Monkeypox_virus", "GCF_014621545.1_Monkeypox_virus_MPXV-M5312_HM12_Rivers_IIb",node)))))
pfof_pastml_rs_sign_state.df <- pfof_pastml_rs_sign_state.df %>% mutate(node=ifelse(str_detect(node, "GCF_|GCA_"), unname(sapply(node, function(x) {
  parts <- str_split(x, "_", simplify = TRUE)
  paste(parts[3:length(parts)], collapse = "_")  # Join everything after the first two parts
})), node))

pfof_sign_state_tr.df <- as.data.frame(ggtree(poxv_t.rt.rn)$data)
pfof_sign_state_tr.df <- pfof_sign_state_tr.df %>% mutate(node_name = ifelse(isTip==FALSE, node, label))

pfof.v <- c("PF06227")
for (a in 1:length(pfof.v)) {
  pfof_sign_state_tr.p_list <- list()

  pfof_og_sign_name.v <- colnames(pfof_pastml_rs_sign_state.df[,2:ncol(pfof_pastml_rs_sign_state.df)])
  pfof_og_sign_name.v <- pfof_og_sign_name.v[str_detect(pfof_og_sign_name.v, paste("Pfam", "\\*", pfof.v[a], sep=""))]
  
  if (pfof.v[a]=="PF06227") {pf06227_og_name.v <- gsub("N0HOG", "N0.HOG", str_split(pfof_og_sign_name.v, '\\^', simplify = TRUE)[,1])}
  
  pfof_sign_state_tr_all_m.df <- poxv_sign_state_tr_all_m.df %>% filter(og_name %in% str_split(pfof_og_sign_name.v, "\\^", simplify=TRUE)[,1], sign_id==pfof.v[a])
  
  for (i in 1:length(pfof_og_sign_name.v)) {
    pfof_sign_state_tr_og_dir <- paste("251208_gain_and_loss_domain_analysis/state_transitions_MPPA/", pfof_og_sign_name.v[i], "/", sep="")
    if (!file.exists(pfof_sign_state_tr_og_dir)) {dir.create(pfof_sign_state_tr_og_dir, recursive = TRUE)}
    
    pfof_pastml_rs_sign_state_og.df <- pfof_pastml_rs_sign_state.df[,c("node",pfof_og_sign_name.v[i])]
    pfof_sign_state_tr_og.df <- pfof_sign_state_tr.df
    pfof_sign_state_tr_og.df$state_parent <- as.character(pfof_pastml_rs_sign_state_og.df[match(pfof_sign_state_tr_og.df$parent, pfof_pastml_rs_sign_state_og.df$node), pfof_og_sign_name.v[i]])
    pfof_sign_state_tr_og.df$state_desc <- as.character(pfof_pastml_rs_sign_state_og.df[match(pfof_sign_state_tr_og.df$node_name, pfof_pastml_rs_sign_state_og.df$node), pfof_og_sign_name.v[i]])
    
    if (nrow(pfof_sign_state_tr_og.df %>% filter(isTip==TRUE, state_desc=="1"))==1) {pfof_sign_state_tr_og.df$state_parent <- "0"}
    
    pfof_sign_state_tr_og.df$og_name <- pfof_og_sign_name.v[i]
    pfof_sign_state_tr_og.df$transition <- paste(pfof_sign_state_tr_og.df$state_parent, pfof_sign_state_tr_og.df$state_desc, sep="->")
    
    pfof_sign_state_tr_og.df <- pfof_sign_state_tr_og.df %>% mutate(noc_name=ifelse(node %in% poxv_noc.v, poxv_noc.df[match(node, poxv_noc.df$node), "node_name"], ""),
                                                                    state_desc=ifelse(state_desc=="0", "Absent", ifelse(state_desc=="0|1", "Undeterminable", ifelse(state_desc=="1", "Present", NA))))
    
    pfof_sign_state_tr_og_all_m.df <- pfof_sign_state_tr_all_m.df %>% filter(og_name %in% str_split(pfof_og_sign_name.v[i], "\\^", simplify=TRUE)[,1])
    # if (nrow(pfof_sign_state_tr_og_all_m.df %>% filter(transition=="0->1"|transition=="1->0"))==0) {next}
    
    pfof_sign_state_tr_og.df <- pfof_sign_state_tr_og.df %>% mutate(concensus=ifelse((parent %in% pfof_sign_state_tr_og_all_m.df$parent) & (node %in% pfof_sign_state_tr_og_all_m.df$node), "yes", "no"))
    pfof_sign_state_tr_og.df <- pfof_sign_state_tr_og.df %>% mutate(gain_loss_status=ifelse(transition=="0->1", "G",
                                                                                     ifelse(transition=="1->0", "L", "")))

    # pfof_gene_dup_events.df <- gene_dup_events.df %>% filter(HOG_genes_1==HOG_genes_2, HOG_genes_1==gsub("N0HOG", "N0.HOG", str_split(pfof_og_sign_name.v[i], '\\^', simplify = TRUE)[,1]))
    pfof_gene_dup_events.df <- gene_dup_events.df %>% filter(HOG_genes_1==gsub("N0HOG", "N0.HOG", str_split(pfof_og_sign_name.v[i], '\\^', simplify = TRUE)[,1]))
    pfof_gene_dup_events.df <- pfof_gene_dup_events.df %>% mutate(Species_Tree_Node=ifelse(Species_Tree_Node %in% poxv_ncbi_ds.df$tip_label, poxv_ncbi_ds.df$tip_label_dummy[match(Species_Tree_Node, poxv_ncbi_ds.df$tip_label)], Species_Tree_Node))
    
    ### Count nodes of duplication within HOG
    if (nrow(pfof_gene_dup_events.df %>% filter(HOG_genes_1==HOG_genes_2)) > 0) {
      pfof_gene_dup_events_ct.df <- pfof_gene_dup_events.df %>% filter(HOG_genes_1==HOG_genes_2) %>% group_by(Species_Tree_Node) %>% summarise(count = n()) %>% ungroup()
      
      get_clade_keys <- function(tree) {
        tree <- as.phylo(tree)
        n_tip <- length(tree$tip.label)
        children <- split(tree$edge[,2], tree$edge[,1])
        get_descendant_tips <- function(node) {
          if (node <= n_tip) {return(tree$tip.label[node])}
          unlist(lapply(children[[as.character(node)]], get_descendant_tips),use.names=FALSE)
        }
        tibble(node=seq_len(n_tip + tree$Nnode), clade_key=sapply(seq_len(n_tip + tree$Nnode), function(node) {paste(sort(get_descendant_tips(node)), collapse="|")}))
      }
      
      of2_node_key.df <- get_clade_keys(poxv_t.rt.rn.wt_of2_node)
      sign_node_key.df <- get_clade_keys(poxv_t.rt.rn)
      
      pfof_poxv_t.rt.rn.wt_of2_node.info.df <- ggtree(poxv_t.rt.rn.wt_of2_node)$data
      pfof_poxv_t.rt.rn.wt_of2_node.info.df <- pfof_poxv_t.rt.rn.wt_of2_node.info.df %>% mutate(dup_count=pfof_gene_dup_events_ct.df$count[match(label, pfof_gene_dup_events_ct.df$Species_Tree_Node)])
      pfof_poxv_t.rt.rn.wt_of2_node.info.df <- pfof_poxv_t.rt.rn.wt_of2_node.info.df %>% mutate(dup_count=ifelse(is.na(dup_count), 0, dup_count))
      
      dup_count_by_clade.df <- pfof_poxv_t.rt.rn.wt_of2_node.info.df %>% select(node, dup_count) %>% left_join(of2_node_key.df, by="node") %>% select(clade_key, dup_count)
      
      pfof_sign_state_tr_og.df <- pfof_sign_state_tr_og.df %>% left_join(sign_node_key.df, by="node") %>% left_join(dup_count_by_clade.df, by="clade_key") %>% select(-clade_key)
      pfof_sign_state_tr_og.df %>% filter(is.na(dup_count))
    } else {
      pfof_sign_state_tr_og.df$dup_count <- 0
    }
    
    ### Count nodes of duplication across HOG
    if (nrow(pfof_gene_dup_events.df %>% filter(HOG_genes_1!=HOG_genes_2)) > 0) {
      pfof_gene_dup_events_ct.df <- pfof_gene_dup_events.df %>% filter(HOG_genes_1!=HOG_genes_2) %>% group_by(Species_Tree_Node) %>% summarise(count = n()) %>% ungroup()
      
      get_clade_keys <- function(tree) {
        tree <- as.phylo(tree)
        n_tip <- length(tree$tip.label)
        children <- split(tree$edge[,2], tree$edge[,1])
        get_descendant_tips <- function(node) {
          if (node <= n_tip) {return(tree$tip.label[node])}
          unlist(lapply(children[[as.character(node)]], get_descendant_tips),use.names=FALSE)
        }
        tibble(node=seq_len(n_tip + tree$Nnode), clade_key=sapply(seq_len(n_tip + tree$Nnode), function(node) {paste(sort(get_descendant_tips(node)), collapse="|")}))
      }
      
      of2_node_key.df <- get_clade_keys(poxv_t.rt.rn.wt_of2_node)
      sign_node_key.df <- get_clade_keys(poxv_t.rt.rn)
      
      pfof_poxv_t.rt.rn.wt_of2_node.info.df <- ggtree(poxv_t.rt.rn.wt_of2_node)$data
      pfof_poxv_t.rt.rn.wt_of2_node.info.df <- pfof_poxv_t.rt.rn.wt_of2_node.info.df %>% mutate(dup_count_acr_og=pfof_gene_dup_events_ct.df$count[match(label, pfof_gene_dup_events_ct.df$Species_Tree_Node)])
      pfof_poxv_t.rt.rn.wt_of2_node.info.df <- pfof_poxv_t.rt.rn.wt_of2_node.info.df %>% mutate(dup_count_acr_og=ifelse(is.na(dup_count_acr_og), 0, dup_count_acr_og))
      
      dup_count_acr_og_by_clade.df <- pfof_poxv_t.rt.rn.wt_of2_node.info.df %>% select(node, dup_count_acr_og) %>% left_join(of2_node_key.df, by="node") %>% select(clade_key, dup_count_acr_og)
      
      pfof_sign_state_tr_og.df <- pfof_sign_state_tr_og.df %>% left_join(sign_node_key.df, by="node") %>% left_join(dup_count_acr_og_by_clade.df, by="clade_key") %>% select(-clade_key)
      pfof_sign_state_tr_og.df %>% filter(is.na(dup_count_acr_og))
    } else {
      pfof_sign_state_tr_og.df$dup_count_acr_og <- 0
    }
    
    pfof_sign_state_tr_og.df <- pfof_sign_state_tr_og.df %>% mutate(dup_status=ifelse(dup_count > 0, "yes", "no"))
    pfof_sign_state_tr_og.df <- pfof_sign_state_tr_og.df %>% mutate(dup_acr_og_status=ifelse(dup_count_acr_og > 0, "yes", "no"))
    
    pfof_sign_state_tr_og.p <- ggtree(poxv_t.rt.rn, size=0.25)
    pfof_sign_state_tr_og.p <- pfof_sign_state_tr_og.p %<+% pfof_sign_state_tr_og.df
    pfof_sign_state_tr_og.p <- pfof_sign_state_tr_og.p + theme_tree()
    # pfof_sign_state_tr_og.p <- pfof_sign_state_tr_og.p + geom_tiplab(size=4, align=FALSE, linesize=.5)
    pfof_sign_state_tr_og.p <- pfof_sign_state_tr_og.p + geom_point2(shape=21, size=3, aes(color=dup_status, fill=state_desc))
    pfof_sign_state_tr_og.p <- pfof_sign_state_tr_og.p + scale_color_manual("dup_status", name="Duplication status", values=c("white","yellow"), breaks=c("no","yes"))
    pfof_sign_state_tr_og.p <- pfof_sign_state_tr_og.p + guides(col="none")
    pfof_sign_state_tr_og.p <- pfof_sign_state_tr_og.p + new_scale_color()
    pfof_sign_state_tr_og.p <- pfof_sign_state_tr_og.p + scale_color_manual("dup_acr_og_status", name="Duplication of HOG", values=c("white","blue"), breaks=c("no","yes"))
    pfof_sign_state_tr_og.p <- pfof_sign_state_tr_og.p + scale_fill_manual("state_desc", name="State", values=c("lightsteelblue1","grey","red"), breaks=c("Absent","Undeterminable","Present"), na.value = NA)
    pfof_sign_state_tr_og.p <- pfof_sign_state_tr_og.p + labs(title=str_split(pfof_og_sign_name.v[i], "\\^", simplify=TRUE)[,1], 
                                                                    subtitle=(paste((poxv_vacv_gene_name.df %>% filter(og_id_dummy==str_split(pfof_og_sign_name.v[i], '\\^', simplify = TRUE)[,1]))$Vaccinia_virus, 
                                                                              "\ndNdS = ", round((poxv_og_dnds.df %>% filter(HOG==gsub("N0HOG", "N0.HOG", str_split(pfof_og_sign_name.v[i], '\\^', simplify = TRUE)[,1])))$dNdS, 2))))
    pfof_sign_state_tr_og.p <- pfof_sign_state_tr_og.p + guides(fill="none", col="none")
    pfof_sign_state_tr_og.p <- pfof_sign_state_tr_og.p + theme(text = element_text(size = 20), plot.title = element_text(size = 20), plot.subtitle = element_text(size = 16))
    pfof_sign_state_tr_og.p <- pfof_sign_state_tr_og.p + xlim(0,4.5)
    pfof_sign_state_tr_og.p <- pfof_sign_state_tr_og.p + coord_cartesian(clip="off")
    
    pfof_sign_state_tr_og.p <- pfof_sign_state_tr_og.p + geom_text2(aes(subset=((transition %in% c("0->1", "1->0")) & (concensus=="yes")), label=ifelse(!is.na(transition), gain_loss_status, "")), size=5, hjust=1.3)
    pfof_sign_state_tr_og.p <- pfof_sign_state_tr_og.p + geom_text2(aes(subset=(dup_count > 0), label=dup_count), size=5, hjust=1.3)
    pfof_sign_state_tr_og.p <- pfof_sign_state_tr_og.p + geom_text2(aes(subset=(dup_count_acr_og > 0), label=paste("*", dup_count_acr_og, sep="")), size=5, hjust=1.3)
    pfof_sign_state_tr.p_list[[pfof_og_sign_name.v[i]]] <- pfof_sign_state_tr_og.p
  }
  
  # ########## ***** Figure 5C: Gain and loss events of protein domains of interest (only poxvirus Bcl-2 like protein domain) in each associated ortholog group *****
  # pdf(paste("Figures/", pfof.v[a], "_gain_and_loss_events.pdf", sep=""), width=24, height=10)
  # print(wrap_plots(pfof_sign_state_tr.p_list) + plot_layout(nrow=1))
  # dev.off()
}

##### Test different between dN/dS ratio of PF06227 and pre-existing ortholog groups
pf06227_poxv_og_dnds.df <- poxv_og_dnds.df %>% filter(gsub("N0HOG", "N0.HOG", HOG) %in% pf06227_og_name.v | gained_status=="pre-existing")
pf06227_poxv_og_dnds.df <- pf06227_poxv_og_dnds.df %>% mutate(group=ifelse(gained_status=="pre-existing", "pre-existing", "PF06227"))
pf06227_poxv_og_dnds.df <- pf06227_poxv_og_dnds.df %>% filter(!is.nan(dNdS))

pf06227_poxv_og_dnds.df %>% group_by(group) %>% shapiro_test(dNdS)
pf06227_poxv_og_dnds.df %>% levene_test(dNdS~group)

### Use Kruskal-Wallis test since normality assumption is violated
pf06227_poxv_og_dnds_res.kruskal <- pf06227_poxv_og_dnds.df %>% kruskal_test(dNdS~group)
pf06227_poxv_og_dnds_res.kruskal

pf06227_poxv_og_dnds_stats <- pf06227_poxv_og_dnds.df %>% dunn_test(dNdS~group, p.adjust.method="fdr")
pf06227_poxv_og_dnds_stats <- pf06227_poxv_og_dnds_stats %>% add_xy_position(x="group", fun="max")
pf06227_poxv_og_dnds_stats

########## ***** Figure 5E: dN/dS ratios of ortholog groups with poxvirus Bcl-2-like protein domains and pre-existing ortholog groups *****
pf06227_poxv_og_dnds_box.p <- ggplot(pf06227_poxv_og_dnds.df, aes(x=group, y=dNdS, col=group), alpha=0.5) 
pf06227_poxv_og_dnds_box.p <- pf06227_poxv_og_dnds_box.p + geom_jitter(alpha=0.1, shape=16, size=3, position=position_jitter(0.25))
pf06227_poxv_og_dnds_box.p <- pf06227_poxv_og_dnds_box.p + geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.75)
pf06227_poxv_og_dnds_box.p <- pf06227_poxv_og_dnds_box.p + scale_color_manual("group", values=c("PF06227"="green","pre-existing"="navyblue"))
pf06227_poxv_og_dnds_box.p <- pf06227_poxv_og_dnds_box.p + labs(title="dN/dS comparison of PF06227 ortholog groups and pre-existing ortholog groups", x="Category", y="dN/dS")
pf06227_poxv_og_dnds_box.p <- pf06227_poxv_og_dnds_box.p + theme_classic()
pf06227_poxv_og_dnds_box.p <- pf06227_poxv_og_dnds_box.p + theme(text = element_text(size = 18), plot.title = element_text(size = 24), plot.subtitle = element_text(size = 22))
pf06227_poxv_og_dnds_box.p <- pf06227_poxv_og_dnds_box.p + scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8), limits=c(0,0.8))
pf06227_poxv_og_dnds_box.p <- pf06227_poxv_og_dnds_box.p + stat_pvalue_manual(pf06227_poxv_og_dnds_stats, label="p.adj", hide.ns=TRUE, inherit.aes=FALSE, size=3)
pf06227_poxv_og_dnds_box.p

# pdf("Figures/pf06227_dNdS.pdf", width=6, height=8)
# print(pf06227_poxv_og_dnds_box.p)
# dev.off()



############################## Horizontal gene transfer analysis ##############################
########## Gather information of proteins assigned to ortholog groups
##### Read the original wide-format HOG result file
hog_res.f <- "251208_HGT_analysis/N0_pooled_additional_assigned.tsv"
hog_res.df <- read.delim(hog_res.f, check.names = FALSE)

hog_res_long.df <- hog_res.df %>% pivot_longer(cols=-c(HOG, OG, `Gene Tree Parent Clade`), names_to = "species_name", values_to = "prot_accessions")
poxv_og_name.v <- unique(hog_res_long.df$HOG)

hog_res_long.df <- hog_res_long.df %>% mutate(prot_accessions = strsplit(as.character(prot_accessions), ",\\s*")) %>% unnest(prot_accessions) %>% filter(prot_accessions != "" & !is.na(prot_accessions)) %>% select(HOG, species_name, prot_accessions)
hog_res_long.df$species_name <- ifelse(hog_res_long.df$species_name=="GCF_000857045.1_Monkeypox_virus", "GCF_000857045.1_Monkeypox_virus_Zaire-96-I-16_Ia",
                                ifelse(hog_res_long.df$species_name=="GCA_039269415.1_Monkeypox_virus", "GCA_039269415.1_Monkeypox_virus_RDC-NKV-GOM-MPOX-010_Ib",
                                ifelse(hog_res_long.df$species_name=="GCA_006465585.1_Monkeypox_virus", "GCA_006465585.1_Monkeypox_virus_Sierra_Leone_IIa",
                                ifelse(hog_res_long.df$species_name=="GCF_014621545.1_Monkeypox_virus", "GCF_014621545.1_Monkeypox_virus_MPXV-M5312_HM12_Rivers_IIb", hog_res_long.df$species_name))))
hog_res_long.df <- hog_res_long.df %>% mutate(species_name=poxv_ncbi_ds.df[match(species_name, poxv_ncbi_ds.df$tip_label), "tip_label_dummy"])
hog_res_long.df <- hog_res_long.df %>% mutate(Subfamily=poxv_ncbi_ds.df[match(species_name, poxv_ncbi_ds.df$tip_label_dummy), "Subfamily"])
hog_res_long.df <- hog_res_long.df %>% mutate(Genus=poxv_ncbi_ds.df[match(species_name, poxv_ncbi_ds.df$tip_label_dummy), "Genus"])

##### Note: Gained ortholog groups may contain ortholog groups arose from gene duplication, which may not come directly from horizontal gene transfer
### Gather information of gene duplication
gene_dup_events.new.hog.df <- gene_dup_events.df %>% filter(HOG_genes_1!=HOG_genes_2) %>% select(HOG_genes_1, HOG_genes_2)

gene_dup_rs.edges <- gene_dup_events.new.hog.df %>% mutate(HOG_genes_1 = strsplit(HOG_genes_1, "\\s*,\\s*"), HOG_genes_2 = strsplit(HOG_genes_2, "\\s*,\\s*")) %>% rowwise() %>% mutate(pairs = list(expand.grid(HOG_genes_1, HOG_genes_2, stringsAsFactors = FALSE))) %>% ungroup() %>%
  unnest(pairs) %>% transmute(from = Var1, to = Var2) %>% mutate(from = trimws(from), to = trimws(to)) %>% filter(from != "", to != "")

gene_dup_rs.g <- graph_from_data_frame(gene_dup_rs.edges, directed = FALSE)

components.df <- components(gene_dup_rs.g)$membership %>% tibble(HOG = names(.), component = as.integer(.))
component_lists <- components.df %>% group_by(component) %>% summarise(hog_set = paste(sort(HOG), collapse = ","), .groups = "drop")

gene_dup_events_rs.new.hog.df <- components.df %>% left_join(component_lists, by = "component") %>% select(HOG, hog_set) %>% arrange(HOG)
gene_dup_events_rs.new.hog.df <- gene_dup_events_rs.new.hog.df %>% mutate(root_status = as.numeric(poxv_state_root_status.df$state_parent[match(HOG, gsub("N0HOG", "N0.HOG", poxv_state_root_status.df$og_name, fixed=TRUE))]))
gene_dup_events_rs.new.hog.df <- gene_dup_events_rs.new.hog.df %>% group_by(hog_set) %>% mutate(pooled_root_status = as.character(max(root_status))) %>% ungroup()

hog_set_res_long.df <- hog_res_long.df %>% mutate(hog_set=ifelse(HOG %in% gene_dup_events_rs.new.hog.df$HOG, gene_dup_events_rs.new.hog.df$hog_set[match(HOG, gene_dup_events_rs.new.hog.df$HOG)], HOG))

##### Count number of proteins in each species in each ortholog group set during poxvirus evolution
hog_sp_prot_ct.df <- hog_set_res_long.df %>% group_by(species_name, hog_set) %>% summarise(count = n()) %>% ungroup()

##### Count number of species shared in each poxvirus subfamily or genus
hog_subfam_sp_ct.df <- hog_set_res_long.df %>% group_by(Subfamily, hog_set) %>% summarise(count = length(unique(species_name))) %>% ungroup()
hog_genus_sp_ct.df <- hog_set_res_long.df %>% group_by(Genus, hog_set) %>% summarise(count = length(unique(species_name))) %>% ungroup()

##### Count number of ortholog groups, number of ortholog group sets, number of ortholog group sets that were gained
hog_ct_per_sp.df <- hog_res_long.df %>% mutate(hog_set = ifelse(HOG %in% gene_dup_events_rs.new.hog.df$HOG, gene_dup_events_rs.new.hog.df$hog_set[match(HOG, gene_dup_events_rs.new.hog.df$HOG)], HOG)) %>% 
  group_by(species_name, HOG) %>% summarise(count = n()) %>% ungroup() %>% group_by(species_name) %>% summarise(og_count = n()) %>% ungroup()

hog_set_ct_per_sp.df <- hog_set_res_long.df %>% group_by(species_name, hog_set) %>% summarise(count = n()) %>% ungroup() %>% group_by(species_name) %>% summarise(hog_set_count = n()) %>% ungroup()

hog_ct_per_sp.df <- hog_ct_per_sp.df %>% mutate(og_count_pooled=hog_set_ct_per_sp.df$hog_set_count[match(species_name, hog_set_ct_per_sp.df$species_name)])

hog_ct_per_sp.df <- hog_ct_per_sp.df %>% mutate(clade=poxv_ncbi_ds.df[match(species_name, poxv_ncbi_ds.df$tip_label_dummy), "clade"])
hog_ct_per_sp.df$species_name <- factor(hog_ct_per_sp.df$species_name, levels=poxv_genome_stats.df$species_name)
hog_ct_per_sp.df$clade <- factor(hog_ct_per_sp.df$clade, levels=levels(poxv_genome_stats.df$clade))


########## Assign states of gain and loss to each ortholog group set
poxv_og_set_state_absent_tip.df <- poxv_state_tr_pooled.df %>% group_by(parent, node, transition, og_name) %>% summarise(transition_count = n()) %>% ungroup()
poxv_og_set_state_absent_tip.df <- poxv_og_set_state_absent_tip.df %>% mutate(transition_count = transition_count/3) %>% filter(transition_count == 1) # Pick only transition events estimated correspondingly by all methods
poxv_og_set_state_absent_tip.df <- poxv_og_set_state_absent_tip.df %>% mutate(HOG=gsub("N0HOG", "N0.HOG", og_name))
poxv_og_set_state_absent_tip.df <- poxv_og_set_state_absent_tip.df %>% mutate(hog_set=ifelse(HOG %in% gene_dup_events_rs.new.hog.df$HOG, gene_dup_events_rs.new.hog.df$hog_set[match(HOG, gene_dup_events_rs.new.hog.df$HOG)], HOG))

poxv_og_set_transition_list.df <- poxv_og_set_state_absent_tip.df %>% group_by(hog_set, parent, node) %>% summarise(transition_set = list(unique(transition)), .groups = "drop")

poxv_og_set_state_absent_tip_ct.df <- poxv_state_tr_pooled.df %>% filter(isTip==TRUE, method=="MPPA", state_desc == "Present") %>% mutate(HOG=gsub("N0HOG", "N0.HOG", og_name))
poxv_og_set_state_absent_tip_ct.df <- poxv_og_set_state_absent_tip_ct.df %>% mutate(hog_set=ifelse(HOG %in% gene_dup_events_rs.new.hog.df$HOG, gene_dup_events_rs.new.hog.df$hog_set[match(HOG, gene_dup_events_rs.new.hog.df$HOG)], HOG))
poxv_og_set_state_absent_tip_ct.df <- poxv_og_set_state_absent_tip_ct.df %>% group_by(hog_set) %>% summarise(tip_set = list(unique(node)), num_tip = length(tip_set[[1]]), .groups = "drop")
poxv_og_set_state_absent_tip_ct.df <- poxv_og_set_state_absent_tip_ct.df %>% mutate(ca_node = map_int(tip_set, ~ (if (length(.x) < 2) {.x[[1]]} else {getMRCA(poxv_t.rt.rn, .x)})),
                                                                                    num_desc = map_int(ca_node, ~ length(Descendants(poxv_t.rt.rn, .x, type = "tips")[[1]])))
poxv_og_set_state_absent_tip_ct.df <- poxv_og_set_state_absent_tip_ct.df %>% left_join(poxv_og_set_transition_list.df %>% rename(ca_node = node) %>% select(hog_set, ca_node, transition_set), by = c("hog_set", "ca_node"))
poxv_og_set_state_absent_tip_ct.df <- poxv_og_set_state_absent_tip_ct.df %>% mutate(transition = map_chr(transition_set, ~ if (length(.x) == 1) {.x} else {paste(.x, collapse = ",")}))
poxv_og_set_state_absent_tip_ct.df <- poxv_og_set_state_absent_tip_ct.df %>% mutate(ca_node_status = map_int(transition, ~ {if (is.na(.x) || .x == "") {return(NA_integer_)}
      parts <- str_split(.x, ",")[[1]]
      rhs <- str_split(parts, "->", simplify = TRUE)[, 2]
      max(as.integer(rhs))
    })
  )
poxv_og_set_state_absent_tip_ct.df <- poxv_og_set_state_absent_tip_ct.df %>% mutate(loss_status = ifelse(!is.na(ca_node_status) & ca_node_status==1, ifelse(num_tip < num_desc, "some_absent", "conserved_after_gain"), "undeterminable")) # Consider only the case that CA is 1
poxv_og_set_state_absent_tip_ct.df <- poxv_og_set_state_absent_tip_ct.df %>% mutate(loss_status = ifelse(ca_node==find_root(poxv_t.rt.rn) & !is.na(ca_node_status) & ca_node_status==1 & num_desc==length(poxv_t.rt.rn$tip.label), "conserved", loss_status))

poxv_og_set_gain_loss_status.df <- poxv_og_gain_loss_status.df %>% mutate(poxv_og_name=gsub("N0HOG", "N0.HOG", poxv_og_name)) %>% rename(HOG=poxv_og_name)
poxv_og_set_gain_loss_status.df <- poxv_og_set_gain_loss_status.df %>% mutate(hog_set=ifelse(HOG %in% gene_dup_events_rs.new.hog.df$HOG, gene_dup_events_rs.new.hog.df$hog_set[match(HOG, gene_dup_events_rs.new.hog.df$HOG)], HOG))
poxv_og_set_gain_loss_status.df <- poxv_og_set_gain_loss_status.df %>% left_join(gene_dup_events_rs.new.hog.df %>% select(HOG, hog_set, pooled_root_status), by=c("HOG","hog_set"))
poxv_og_set_gain_loss_status.df <- poxv_og_set_gain_loss_status.df %>% mutate(pooled_root_status=ifelse(HOG==hog_set, root_status, pooled_root_status))
poxv_og_set_gain_loss_status.df <- poxv_og_set_gain_loss_status.df %>% group_by(hog_set) %>% mutate(pooled_tip_status=poxv_og_set_state_absent_tip_ct.df$loss_status[match(hog_set, poxv_og_set_state_absent_tip_ct.df$hog_set)],
                                                                                                    pooled_detected_gain_event=ifelse("yes" %in% detected_gain_event, "yes", "no"),
                                                                                                    pooled_detected_loss_event=ifelse("yes" %in% detected_loss_event, "yes", "no")) %>% ungroup()
poxv_og_set_gain_loss_status.df <- poxv_og_set_gain_loss_status.df %>% mutate(pooled_gain_loss_status=ifelse(pooled_tip_status=="conserved", "conserved", "undeterminable"))
poxv_og_set_gain_loss_status.df <- poxv_og_set_gain_loss_status.df %>% mutate(pooled_gain_loss_status=ifelse(pooled_detected_gain_event=="yes" | (!is.na(pooled_root_status) & pooled_root_status==0),
                                                                                               ifelse(pooled_detected_loss_event=="yes" | pooled_tip_status=="some_absent", "both", "gain"),
                                                                                               ifelse(pooled_detected_loss_event=="yes" | pooled_tip_status=="some_absent", "loss", pooled_gain_loss_status)))
poxv_og_set_gain_loss_status.df <- poxv_og_set_gain_loss_status.df %>% select(hog_set, pooled_root_status, pooled_tip_status, pooled_detected_gain_event, pooled_detected_loss_event, pooled_gain_loss_status)
poxv_og_set_gain_loss_status.df <- poxv_og_set_gain_loss_status.df %>% distinct()
poxv_og_set_gain_loss_status.df <- poxv_og_set_gain_loss_status.df %>% rename(root_status=pooled_root_status, tip_status=pooled_tip_status, detected_gain_event=pooled_detected_gain_event, detected_loss_event=pooled_detected_loss_event, gain_loss_status=pooled_gain_loss_status)
poxv_og_set_gain_loss_status.df <- poxv_og_set_gain_loss_status.df %>% mutate(gained_status=ifelse(!is.na(root_status) & root_status==0, "gained",
                                                                                            ifelse(!is.na(root_status) & root_status==1, "pre-existing","undeterminable")))
poxv_og_set_gain_loss_status.df$gain_loss_status <- factor(poxv_og_set_gain_loss_status.df$gain_loss_status, levels=c("gain","loss","both","conserved","undeterminable"))
poxv_og_set_gain_loss_status.df$gained_status <- factor(poxv_og_set_gain_loss_status.df$gained_status, levels=c("gained","pre-existing","undeterminable"))


########## Assign names to each ortholog group set
poxv_vacv_og_set_name.df <- poxv_vacv_gene_name.df %>% mutate(hog_set=ifelse(og_id %in% gene_dup_events_rs.new.hog.df$HOG, gene_dup_events_rs.new.hog.df$hog_set[match(og_id, gene_dup_events_rs.new.hog.df$HOG)], og_id))
poxv_vacv_og_set_name.df <- poxv_vacv_og_set_name.df %>% mutate(vaccinia_tok = str_split(Vaccinia_virus, "\\|")) %>% unnest(vaccinia_tok) %>% mutate(vaccinia_tok = str_trim(vaccinia_tok))
poxv_vacv_og_set_name.df <- poxv_vacv_og_set_name.df %>% distinct(hog_set, vaccinia_tok) %>% group_by(hog_set) %>% summarise(pooled_vv_name = {
  toks <- vaccinia_tok[vaccinia_tok != ""]
  if (length(toks) == 0) {""} else {paste(toks, collapse = "|")}},
  .groups = "drop"
  )

poxv_vacv_og_set_abbr_name.df <- poxv_vacv_gene_name.df %>% mutate(hog_set=ifelse(og_id %in% gene_dup_events_rs.new.hog.df$HOG, gene_dup_events_rs.new.hog.df$hog_set[match(og_id, gene_dup_events_rs.new.hog.df$HOG)], og_id))
poxv_vacv_og_set_abbr_name.df <- poxv_vacv_og_set_abbr_name.df %>% mutate(abbr_tok = str_split(abbreviated_name, "\\|")) %>% unnest(abbr_tok) %>% mutate(abbr_tok = str_trim(abbr_tok))
poxv_vacv_og_set_abbr_name.df <- poxv_vacv_og_set_abbr_name.df %>% distinct(hog_set, abbr_tok) %>% group_by(hog_set) %>% summarise(hog_set_name = {
  toks <- abbr_tok[abbr_tok != ""]
  if (length(toks) == 0) {""} else {paste(toks, collapse = "|")}},
  .groups = "drop"
)

poxv_vacv_og_set_name.df <- poxv_vacv_og_set_name.df %>% mutate(pooled_abbr_name=poxv_vacv_og_set_abbr_name.df$hog_set_name[match(hog_set, poxv_vacv_og_set_abbr_name.df$hog_set)])


########## Count number of gained ortholog groups, gained ortholog group sets, and gained ortholog group sets with duplication events
gained_hog_ct_per_sp.df <- hog_res_long.df %>% mutate(gained_status = poxv_og_gain_loss_status.df$gained_status[match(HOG, gsub("N0HOG", "N0.HOG", poxv_og_gain_loss_status.df$poxv_og_name, fixed=TRUE))]) %>% filter(gained_status=="gained") %>% 
  group_by(species_name, HOG) %>% summarise(count = n()) %>% ungroup() %>% group_by(species_name) %>% summarise(gained_hog_count = n()) %>% ungroup()
# gained_hog_set_ct_per_sp.df <- hog_set_res_long.df %>% mutate(gained_status = poxv_og_set_gain_loss_status.df$gained_status[match(hog_set, poxv_og_set_gain_loss_status.df$hog_set)]) %>% filter(gained_status=="gained") %>% 
#   group_by(species_name, hog_set) %>% summarise(count = n()) %>% ungroup() %>% group_by(species_name) %>% summarise(gained_hog_set_count = n()) %>% ungroup()

hog_ct_per_sp.df <- hog_ct_per_sp.df %>% mutate(gained_hog_num=gained_hog_ct_per_sp.df$gained_hog_count[match(species_name, gained_hog_ct_per_sp.df$species_name)])
hog_ct_per_sp.df <- hog_ct_per_sp.df %>% mutate(prop_gained_hog_num=gained_hog_num/og_count)

# hog_ct_per_sp.df <- hog_ct_per_sp.df %>% mutate(gained_hog_set_num=gained_hog_set_ct_per_sp.df$gained_hog_set_count[match(species_name, gained_hog_set_ct_per_sp.df$species_name)])
# hog_ct_per_sp.df <- hog_ct_per_sp.df %>% mutate(prop_gained_hog_set_num=gained_hog_set_num/og_count)


########## ***** Figure 2B: Proportion of gained, pre-existing, and undeterminable genes *****
gene_gained_status_hog_ct_per_sp.df <- hog_res_long.df %>% mutate(gained_status = poxv_og_gain_loss_status.df$gained_status[match(HOG, gsub("N0HOG", "N0.HOG", poxv_og_gain_loss_status.df$poxv_og_name, fixed=TRUE))]) %>% 
  group_by(species_name, HOG, gained_status) %>% summarise(count = n()) %>% ungroup() %>% group_by(species_name, gained_status) %>% summarise(hog_count = n()) %>% ungroup()

gene_gained_status_hog_ct_per_sp.df$gained_status <- factor(gene_gained_status_hog_ct_per_sp.df$gained_status, levels=c("undeterminable","pre-existing","gained"))
gene_gained_status_hog_ct_per_sp.df$species_name <- factor(gene_gained_status_hog_ct_per_sp.df$species_name, levels=rev(poxv_t.p$data %>% filter(isTip) %>% arrange(desc(y)) %>% pull(label)))

gene_gained_status_hog_ct_per_sp.df <- gene_gained_status_hog_ct_per_sp.df %>% mutate(og_count=hog_ct_per_sp.df$og_count[match(species_name, hog_ct_per_sp.df$species_name)])

summary(gene_gained_status_hog_ct_per_sp.df %>% mutate(prop_gained_hog_num=hog_count/og_count) %>% filter(gained_status=="gained") %>% pull(prop_gained_hog_num))*100 ### Use this to count proportions of gained genes

prop_poxv_gene_by_gained_status.p <- ggplot(gene_gained_status_hog_ct_per_sp.df, aes(x=hog_count, y=species_name, fill=gained_status))
prop_poxv_gene_by_gained_status.p <- prop_poxv_gene_by_gained_status.p + geom_bar(stat="identity", position="fill", width=0.5)
prop_poxv_gene_by_gained_status.p <- prop_poxv_gene_by_gained_status.p + scale_fill_manual("gained_status", values = c("gained"="orange","pre-existing"="navyblue","undeterminable"="lightgrey"))
prop_poxv_gene_by_gained_status.p <- prop_poxv_gene_by_gained_status.p + labs(title="Proportion of gained, pre-existing, and undeterminable genes", x="Length", y="Virus")
prop_poxv_gene_by_gained_status.p <- prop_poxv_gene_by_gained_status.p + theme_classic()
prop_poxv_gene_by_gained_status.p <- prop_poxv_gene_by_gained_status.p + theme(text = element_text(size = 18), plot.title = element_text(size = 24), plot.subtitle = element_text(size = 22))
prop_poxv_gene_by_gained_status.p

# pdf("Figures/hog_prop_by_gained_status.pdf", width=16, height=16)
# print(prop_poxv_gene_by_gained_status.p)
# dev.off()

##### Count number of gained ortholog groups that further underwent gene duplication
gained_dupl.df <- hog_set_res_long.df %>% mutate(gained_status = poxv_og_set_gain_loss_status.df$gained_status[match(hog_set, poxv_og_set_gain_loss_status.df$hog_set)]) %>%
  filter(gained_status=="gained", hog_set %in% unique(c(component_lists$hog_set, gene_dup_events.df$HOG_genes_1, gene_dup_events.df$HOG_genes_2)))

gained_dupl_ct.df <- gained_dupl.df %>% group_by(species_name, HOG) %>%
    summarise(count = n()) %>% ungroup() %>% group_by(species_name) %>% summarise(gained_dupl_hog_set_num = n()) %>% ungroup()

hog_ct_per_sp.df <- hog_ct_per_sp.df %>% mutate(gained_dupl_hog_set_num=gained_dupl_ct.df$gained_dupl_hog_set_num[match(species_name, gained_dupl_ct.df$species_name)])


########## Check source of poxvirus genes
##### Visualize number of genes generated from HGT
hgt_res.f <- "251208_HGT_analysis/251208_HGT_result.tsv"
hgt_res.df <- fread(hgt_res.f, header=T, sep="\t", quote="", check.names=T, data.table=FALSE)
hgt_res.df <- hgt_res.df %>% mutate(species_name=ifelse(species_name=="N/A",
                                                        as.character(poxv_genome_stats.df$species_name)[match(paste(str_split(gsub("_[^_]+$", "", gsub(".*vs_", "", qseqid)), ".1_", simplify=TRUE)[,1], ".1", sep=""), poxv_genome_stats.df$accession_no)],
                                                        as.character(poxv_genome_stats.df$species_name)[match(species_name, poxv_genome_stats.df$assembly_name)]))
hgt_res.df$species_name <- factor(hgt_res.df$species_name, levels=as.vector(poxv_genome_stats.df$species_name))

##### Check whether the best hits fall into non-poxvirus virus or non-virus category
hgt_res_cprs.df <- hgt_res.df %>% select(qseqid, best_type, bitscore) %>% group_by(qseqid, best_type) %>% slice_max(order_by=bitscore, n=1, with_ties=FALSE) %>% ungroup()
hgt_res_cprs_mat <- as.data.frame.matrix(as.matrix(xtabs(bitscore ~ qseqid + best_type, data = hgt_res_cprs.df)))
hgt_res_cprs_mat <- hgt_res_cprs_mat %>% mutate(hgt_status=ifelse(best_non_virus_non_poxvirus > best_virus_non_poxvirus, "hgt_non_virus",
                                                           ifelse(best_non_virus_non_poxvirus < best_virus_non_poxvirus, "hgt_virus", "N/A")))

##### Filter for hits which are likely derived from HGT
##### Assigned likely HGT state to each protein
hgt_res.np.df <- hgt_res.df %>% mutate(hgt_status=hgt_res_cprs_mat[match(qseqid, rownames(hgt_res_cprs_mat)),"hgt_status"])
hgt_res.np.df <- hgt_res.np.df %>% filter((best_type == "best_non_virus_non_poxvirus" & hgt_status == "hgt_non_virus") | (best_type == "best_virus_non_poxvirus" & hgt_status == "hgt_virus"))
hgt_res.np.df <- hgt_res.np.df %>% mutate(HOG_ID=hog_res_long.df$HOG[match(qseqid, hog_res_long.df$prot_accessions)])

# ##### Filter for only hits to genes with gain events (0->1 or root = 0)
# ### Exclude hits from duplicated ortholog groups whose root is not present (root !=0 or root=NA)
# gained_prot_long.df <- hog_set_res_long.df %>% filter(hog_set %in% (poxv_og_set_gain_loss_status.df %>% filter(gained_status=="gained") %>% pull(hog_set)))

##### Identify ortholog group set likely to originate from HGT
hog_prot_ct.df <- hog_set_res_long.df %>% group_by(hog_set) %>% summarise(prot_count = n()) %>% ungroup()

hgt_og_set_status.df <- hgt_res.np.df %>% mutate(hog_set = ifelse(HOG_ID %in% gene_dup_events_rs.new.hog.df$HOG, gene_dup_events_rs.new.hog.df$hog_set[match(HOG_ID, gene_dup_events_rs.new.hog.df$HOG)], HOG_ID))
hgt_og_set_status.df <- hgt_og_set_status.df %>% group_by(qseqid, sscinames) %>% slice_head(n=1) %>% ungroup()                                                                                                    # Pick only best hit from a subject species
hgt_og_set_status.df <- hgt_og_set_status.df %>% filter(qcovhsp > 50) %>% group_by(hog_set, qseqid) %>% summarise(hit_count=n()) %>% ungroup()                                                                    # Filter by % coverage and count hit numbers
hgt_og_set_status.df <- hgt_og_set_status.df %>% mutate(hgt_stat=ifelse(hit_count >= 1, 1, 0)) %>% group_by(hog_set) %>% summarise(hgt_prot=sum(hgt_stat))                                                        # Count proteins with hits in each HOG set
hgt_og_set_status.df <- hgt_og_set_status.df %>% left_join(hog_prot_ct.df, by=c("hog_set")) %>% mutate(prop_hgt = hgt_prot/prot_count)                                                                            # Calculate proportion of proteins with HGT hits
hgt_og_set_status.df <- hgt_og_set_status.df %>% mutate(gained_status=poxv_og_set_gain_loss_status.df$gained_status[match(hog_set, poxv_og_set_gain_loss_status.df$hog_set)])
hgt_og_set_status.df <- hgt_og_set_status.df %>% mutate(hgt_status=ifelse(prop_hgt > 0.5, "yes", "no"))
hgt_og_set_status.df <- hgt_og_set_status.df %>% mutate(gained_hgt_status=ifelse(gained_status=="gained" & hgt_status=="yes", "yes", "no"))

hist(hgt_og_set_status.df$prop_hgt, breaks = 10) # Distribution of percentage coverage
hist(hgt_og_set_status.df$prop_hgt[hgt_og_set_status.df$prop_hgt > 0.5], breaks = 10) # Distribution of percentage coverage

##### Count number of ortholog group sets with HGT events
gained_hgt_ct_per_sp.df <- hog_set_res_long.df %>% mutate(gained_status = poxv_og_set_gain_loss_status.df$gained_status[match(hog_set, poxv_og_set_gain_loss_status.df$hog_set)], hgt_status = hgt_og_set_status.df$hgt_status[match(hog_set, hgt_og_set_status.df$hog_set)]) %>%
  filter(gained_status=="gained", hgt_status=="yes") %>% group_by(species_name, hog_set) %>% summarise(count = n()) %>% ungroup() %>% group_by(species_name) %>% summarise(gained_hgt_hog_set_num = n()) %>% ungroup()

hog_ct_per_sp.df <- hog_ct_per_sp.df %>% mutate(gained_hgt_hog_set_num=gained_hgt_ct_per_sp.df$gained_hgt_hog_set_num[match(species_name, gained_hgt_ct_per_sp.df$species_name)])

# gained_dupl_hgt_ct_per_sp.df <- hog_set_res_long.df %>% mutate(gained_status = poxv_og_set_gain_loss_status.df$gained_status[match(hog_set, poxv_og_set_gain_loss_status.df$hog_set)], hgt_status = hgt_og_set_status.df$hgt_status[match(hog_set, hgt_og_set_status.df$hog_set)]) %>% 
#   filter(gained_status=="gained", hog_set %in% unique(c(component_lists$hog_set, gene_dup_events.df$HOG_genes_1, gene_dup_events.df$HOG_genes_2)), hgt_status=="yes") %>% group_by(species_name, hog_set) %>% summarise(count = n()) %>% ungroup() %>% group_by(species_name) %>% summarise(gained_dupl_hgt_hog_set_num = n()) %>% ungroup()
# 
# hog_ct_per_sp.df <- hog_ct_per_sp.df %>% mutate(gained_dupl_hgt_hog_set_num=gained_dupl_hgt_ct_per_sp.df$gained_dupl_hgt_hog_set_num[match(species_name, gained_dupl_hgt_ct_per_sp.df$species_name)])
# hog_ct_per_sp.df <- hog_ct_per_sp.df %>% mutate(gained_dupl_hgt_hog_set_num=ifelse(is.na(gained_dupl_hgt_hog_set_num), 0, gained_dupl_hgt_hog_set_num))
# 
hog_ct_per_sp.df <- hog_ct_per_sp.df %>% mutate(prop_gained_dupl_hog_set_num=gained_dupl_hog_set_num/og_count)
hog_ct_per_sp.df <- hog_ct_per_sp.df %>% mutate(prop_gained_hgt_hog_set_num=gained_hgt_hog_set_num/og_count)
# hog_ct_per_sp.df <- hog_ct_per_sp.df %>% mutate(prop_gained_dupl_hgt_hog_set_num=gained_dupl_hgt_hog_set_num/og_count)

poxv_genome_stats.df <- poxv_genome_stats.df %>% mutate(hgt_likely_number=hog_ct_per_sp.df$gained_hgt_hog_set_num[match(species_name, hog_ct_per_sp.df$species_name)])
poxv_genome_stats.df <- poxv_genome_stats.df %>% mutate(hog_number=hog_ct_per_sp.df$og_count[match(species_name, hog_ct_per_sp.df$species_name)])


##### Poxvirus tree with proportion of gained and pre-existing genes contributing to genome expansion
gene_dup_prot_list.v <- gene_dup_events.df %>% select('Genes 1', 'Genes 2') %>% unlist(use.names = FALSE) %>% strsplit(",\\s*") %>% unlist(use.names = FALSE) %>% unique()
gene_dup_prot_list.v <- unname(sapply(gene_dup_prot_list.v, function(x) {
  parts <- str_split(x, "_", simplify = TRUE)
  num_parts <- length(parts)
  
  # Check if the second-to-last part is "NP", "XP", or "YP"
  if (parts[num_parts - 1] %in% c("NP", "XP", "YP")) {
    paste(parts[num_parts - 1], "_", parts[num_parts], sep="")
  } else if (TRUE %in% str_detect(x, "N0.HOG")) {
    paste(parts[(which(str_detect(parts, "N0.HOG")==TRUE)):num_parts], collapse = "_")
  } else {
    parts[num_parts]
  }
}))


##### Proportion of ortholog groups expanding poxvirus genomes in each species
gene_sum_status_hog_ct_per_sp.df <- hog_set_res_long.df %>% mutate(gained_status=poxv_og_set_gain_loss_status.df$gained_status[match(hog_set, gsub("N0HOG", "N0.HOG", poxv_og_set_gain_loss_status.df$hog_set, fixed=TRUE))])
gene_sum_status_hog_ct_per_sp.df <- gene_sum_status_hog_ct_per_sp.df %>% mutate(hgt_status=ifelse(hog_set %in% (hgt_og_set_status.df %>% filter(gained_hgt_status=="yes") %>% pull(hog_set)), "hgt", "non-hgt"))
gene_sum_status_hog_ct_per_sp.df <- gene_sum_status_hog_ct_per_sp.df %>% mutate(dupl_status=ifelse(prot_accessions %in% gene_dup_prot_list.v, "dupl", "not-dupl"))
gene_sum_status_hog_ct_per_sp.df <- gene_sum_status_hog_ct_per_sp.df %>% mutate(sum_status=paste(gained_status, hgt_status, dupl_status, sep="_"))
gene_sum_status_hog_ct_per_sp.df <- gene_sum_status_hog_ct_per_sp.df %>% group_by(species_name, HOG, sum_status) %>% summarise(count = n()) %>% ungroup() %>% group_by(species_name, sum_status) %>% summarise(hog_count = n()) %>% ungroup()
gene_sum_status_hog_ct_per_sp.df <- gene_sum_status_hog_ct_per_sp.df %>% filter(!sum_status %in% c("undeterminable_non-hgt_dupl","undeterminable_non-hgt_not-dupl","pre-existing_non-hgt_not-dupl"))

gene_sum_status_hog_ct_per_sp.df$sum_status <- factor(gene_sum_status_hog_ct_per_sp.df$sum_status, levels=c("pre-existing_non-hgt_dupl","gained_non-hgt_dupl","gained_hgt_dupl","gained_hgt_not-dupl","gained_non-hgt_not-dupl"))
gene_sum_status_hog_ct_per_sp.df$species_name <- factor(gene_sum_status_hog_ct_per_sp.df$species_name, levels=rev(poxv_t.p$data %>% filter(isTip) %>% arrange(desc(y)) %>% pull(label)))

gene_sum_status_hog_ct_per_sp_gained.df <- gene_sum_status_hog_ct_per_sp.df %>% filter(sum_status %in% c("gained_non-hgt_dupl","gained_hgt_dupl","gained_hgt_not-dupl","gained_non-hgt_not-dupl"))
gene_sum_status_hog_ct_per_sp_gained.df <- gene_sum_status_hog_ct_per_sp_gained.df %>% group_by(species_name) %>% mutate(sum_gained = sum(hog_count)) %>% ungroup()

##### Calculate average percentage of gained ortholog groups that further underwent gene duplication
summary(gene_sum_status_hog_ct_per_sp_gained.df %>% mutate(prop_hog_num=hog_count/sum_gained) %>% filter(sum_status %in% c("gained_hgt_dupl","gained_non-hgt_dupl")) %>% group_by(species_name) %>% summarise(sum_prop_hog_num = sum(prop_hog_num)) %>% ungroup() %>% pull(sum_prop_hog_num))*100

##### Calculate average percentage of gained ortholog groups derived from HGT
summary(gene_sum_status_hog_ct_per_sp_gained.df %>% mutate(prop_hog_num=hog_count/sum_gained) %>% filter(sum_status %in% c("gained_hgt_dupl","gained_hgt_not-dupl")) %>% group_by(species_name) %>% summarise(sum_prop_hog_num = sum(prop_hog_num)) %>% ungroup() %>% pull(sum_prop_hog_num))*100

##### Proportion of ortholog groups expanding poxvirus genomes (average)
gene_sum_status_hog_ct_per_sp_gained_av.df <- gene_sum_status_hog_ct_per_sp_gained.df %>% mutate(sum_status=factor(sum_status, levels=c("gained_non-hgt_not-dupl","gained_non-hgt_dupl","gained_hgt_dupl","gained_hgt_not-dupl")))
gene_sum_status_hog_ct_per_sp_gained_av.df <- gene_sum_status_hog_ct_per_sp_gained_av.df  %>% group_by(species_name, sum_status) %>% summarise(hog_count=sum(hog_count), .groups = "drop") %>% tidyr::complete(species_name, sum_status, fill = list(hog_count = 0)) %>% group_by(species_name) %>%
  mutate(prop = hog_count / sum(hog_count)) %>% ungroup() %>% group_by(sum_status) %>% summarise(av_prop = mean(prop), .groups = "drop")

gene_sum_status_hog_ct_per_sp_gained_av.df$sum_status <- factor(gene_sum_status_hog_ct_per_sp_gained_av.df$sum_status, levels=c("gained_non-hgt_not-dupl","gained_non-hgt_dupl","gained_hgt_dupl","gained_hgt_not-dupl"))

sum_status_color.v <- brewer.pal(5, "Paired")
names(sum_status_color.v) <- c("pre-existing_non-hgt_dupl","gained_non-hgt_dupl","gained_hgt_dupl","gained_hgt_not-dupl","gained_non-hgt_not-dupl")

########## ***** Figure 2H: Proportion of gained ortholog groups classified by duplication and HGT status *****
av_prop_poxv_gene_by_sum_status.p <- ggplot(gene_sum_status_hog_ct_per_sp_gained_av.df, aes(x="Average across species", y=av_prop, fill=sum_status))
av_prop_poxv_gene_by_sum_status.p <- av_prop_poxv_gene_by_sum_status.p + geom_bar(stat="identity", width=0.5)
av_prop_poxv_gene_by_sum_status.p <- av_prop_poxv_gene_by_sum_status.p + scale_fill_manual("sum_status", values=sum_status_color.v)
av_prop_poxv_gene_by_sum_status.p <- av_prop_poxv_gene_by_sum_status.p + labs(title="Average proportion of gained ortholog groups", x="Proportion", y="")
av_prop_poxv_gene_by_sum_status.p <- av_prop_poxv_gene_by_sum_status.p + theme_classic()
av_prop_poxv_gene_by_sum_status.p <- av_prop_poxv_gene_by_sum_status.p + theme(text = element_text(size = 18), plot.title = element_text(size = 24), plot.subtitle = element_text(size = 22))
av_prop_poxv_gene_by_sum_status.p

# pdf("Figures/hog_av_prop_by_sum_status.pdf", width=16, height=16)
# print(av_prop_poxv_gene_by_sum_status.p)
# dev.off()

##### Filter for only hits of HGT-derived ortholog groups
hgt_gained_res.np.df <- hgt_res.np.df %>% mutate(hog_set = ifelse(HOG_ID %in% gene_dup_events_rs.new.hog.df$HOG, gene_dup_events_rs.new.hog.df$hog_set[match(HOG_ID, gene_dup_events_rs.new.hog.df$HOG)], HOG_ID))
hgt_gained_res.np.df <- hgt_gained_res.np.df %>% filter(hog_set %in% (hgt_og_set_status.df %>% filter(gained_hgt_status=="yes") %>% pull(hog_set)))

########## ***** Supplementary Figure 6: Number of HGT-derived proteins separated by percentage query coverage per hsp *****
hgt_gained_by_perc_cov.df.p <- ggplot(hgt_gained_res.np.df %>% filter(qcovhsp > 50), aes(x=qcovhsp))
hgt_gained_by_perc_cov.df.p <- hgt_gained_by_perc_cov.df.p + geom_histogram(binwidth = 5, position = "identity")
hgt_gained_by_perc_cov.df.p <- hgt_gained_by_perc_cov.df.p + theme_classic()
hgt_gained_by_perc_cov.df.p <- hgt_gained_by_perc_cov.df.p + theme(text=element_text(size = 18), plot.title=element_text(size = 32))
hgt_gained_by_perc_cov.df.p <- hgt_gained_by_perc_cov.df.p + labs(x="Counts of HGT-derived proteins separated by qcovhsp")
hgt_gained_by_perc_cov.df.p <- hgt_gained_by_perc_cov.df.p + scale_x_continuous(breaks=seq(50, 100, by=10))
hgt_gained_by_perc_cov.df.p <- hgt_gained_by_perc_cov.df.p + scale_y_continuous(breaks=seq(0, 3000, by=500), limits=c(0,3000))
hgt_gained_by_perc_cov.df.p <- hgt_gained_by_perc_cov.df.p + geom_vline(xintercept = mean((hgt_gained_res.np.df %>% filter(qcovhsp > 50))$qcovhsp), linetype = "dashed", color = "black")

# pdf("Figures/qcovhsp_protein_encoded_from_putative_hgt_genes.pdf", width=16, height=12)
# print(hgt_gained_by_perc_cov.df.p)
# dev.off()

gained_hgt_og.v.f <- "260107_HGT_analysis/260107_gained_hog_list.tsv"
# write.table(unique(sort(hgt_gained_res.np.df$HOG_ID)), gained_hgt_og.v.f, col.names=F, row.names=F, sep="\n", quote=F)

gained_hgt_og_set.v.f <- "260107_HGT_analysis/260107_gained_hog_set_list.tsv"
# write.table(hgt_og_set_status.df %>% filter(gained_hgt_status=="yes") %>% pull(hog_set), gained_hgt_og_set.v.f, col.names=F, row.names=F, sep="\n", quote=F)


########## ***** Figure 2G, right: Number of HGT-derived proteins separated by percentage query coverage per hsp *****
poxv_state_tr_all_m_hgt.df <- poxv_state_tr_all_m.df %>% mutate(og_name=gsub("N0HOG", "N0.HOG", og_name, fixed=TRUE))
poxv_state_tr_all_m_hgt.df <- poxv_state_tr_all_m_hgt.df %>% mutate(hog_set=ifelse(og_name %in% gene_dup_events_rs.new.hog.df$HOG, gene_dup_events_rs.new.hog.df$hog_set[match(og_name, gene_dup_events_rs.new.hog.df$HOG)], og_name))
poxv_state_tr_all_m_hgt.df <- poxv_state_tr_all_m_hgt.df %>% filter(hog_set %in% (hgt_og_set_status.df %>% filter(gained_hgt_status=="yes") %>% pull(hog_set)))
poxv_state_tr_all_m_hgt.df <- poxv_state_tr_all_m_hgt.df %>% filter(transition == "0->1")

poxv_state_tr_all_m_hgt_ct.df <- poxv_state_tr_all_m_hgt.df %>% mutate(og_num=lengths(str_split(hog_set, "\\s*,\\s*")))
poxv_state_tr_all_m_hgt_ct.df <- poxv_state_tr_all_m_hgt_ct.df %>% group_by(hog_set) %>% mutate(ca_node=ifelse(unique(og_num)>1 & length(node)>1, getMRCA(poxv_t.rt.rn, node), NA)) %>% ungroup()
poxv_state_tr_all_m_hgt_ct.df <- poxv_state_tr_all_m_hgt_ct.df %>% group_by(hog_set) %>% filter(!(n() > 1 & any(!is.na(ca_node))) | parent==ca_node) %>% ungroup()
poxv_state_tr_all_m_hgt_ct.df <- poxv_state_tr_all_m_hgt_ct.df %>% group_by(hog_set) %>% mutate(other_parent_is_anc = sapply(seq_len(n()), function(i){any(parent[-i] %in% unlist(Ancestors(poxv_t.rt.rn, parent[i], type = "all")))})) %>% ungroup()
poxv_state_tr_all_m_hgt_ct.df <- poxv_state_tr_all_m_hgt_ct.df %>% filter(other_parent_is_anc==FALSE) %>% group_by(hog_set, parent, node) %>% summarise(num_transition=1) %>% ungroup()
poxv_state_tr_all_m_hgt_ct.df <- poxv_state_tr_all_m_hgt_ct.df %>% group_by(parent, node) %>% summarise(transition_count = n()) %>% ungroup()

poxv_state_tr_all_m_total_hgt_ct.info.df <- poxv_t.rt.info.df %>% left_join(poxv_state_tr_all_m_hgt_ct.df, by = c("parent", "node"))
poxv_state_tr_all_m_total_hgt_ct.info.df <- poxv_state_tr_all_m_total_hgt_ct.info.df %>% mutate(transition_count=ifelse(is.na(transition_count), 0, transition_count))

########## ***** Figure 2G, right: Number of HGT events in each branch *****
poxv_state_tr_all_m_hgt_ct.p <- ggtree(poxv_t.rt.rn, linewidth=2, aes(color=transition_count))
poxv_state_tr_all_m_hgt_ct.p <- poxv_state_tr_all_m_hgt_ct.p %<+% poxv_state_tr_all_m_total_hgt_ct.info.df
poxv_state_tr_all_m_hgt_ct.p <- poxv_state_tr_all_m_hgt_ct.p + theme_tree()
poxv_state_tr_all_m_hgt_ct.p <- poxv_state_tr_all_m_hgt_ct.p + geom_text2(aes(label=round((transition_count), 1)), color="black", size=5, hjust=1.3)
poxv_state_tr_all_m_hgt_ct.p <- poxv_state_tr_all_m_hgt_ct.p + scale_color_gradient("transition_count", trans="pseudo_log", name="Number of HGT events", low="black", high="violet", breaks=c(0,3,5,10,30,50,100,max(poxv_state_tr_all_m_total_hgt_ct.info.df$transition_count)))
poxv_state_tr_all_m_hgt_ct.p <- poxv_state_tr_all_m_hgt_ct.p + labs(title="All methods", subtitle="HGT")
poxv_state_tr_all_m_hgt_ct.p <- poxv_state_tr_all_m_hgt_ct.p + theme(text = element_text(size = 20), plot.title = element_text(size = 30), plot.subtitle = element_text(size = 24))
poxv_state_tr_all_m_hgt_ct.p <- poxv_state_tr_all_m_hgt_ct.p + xlim(0,4.5)
poxv_state_tr_all_m_hgt_ct.p

# pdf("Figures/hgt_summarized.pdf", width=12, height=10)
# print(poxv_state_tr_all_m_hgt_ct.p)
# dev.off()


########## Visualize percentage identity of all genes separated by Subfamily
hgt_gained_p_iden.filtered.df <- hgt_gained_res.np.df %>% group_by(qseqid, sscinames) %>% slice_head(n=1) %>% filter(qcovhsp > 50) %>% ungroup()                                                                       # Pick only best hit from a subject species, filter by % coverage
hgt_gained_p_iden.filtered.df <- hgt_gained_p_iden.filtered.df %>% mutate(clade=poxv_genome_stats.df$clade[match(species_name, poxv_genome_stats.df$species_name)])
hgt_gained_p_iden.filtered.df <- hgt_gained_p_iden.filtered.df %>% mutate(Subfamily=poxv_genome_stats.df$Subfamily[match(species_name, poxv_genome_stats.df$species_name)])

hgt_gained_p_iden.filtered.w.df <- hgt_gained_p_iden.filtered.df %>% group_by(qseqid) %>% mutate(w_qseqid = 1/n()) %>% ungroup()                                                                                       # Calculate weight for each query protein to normalize counts by number of hits for each query protein
hgt_gained_p_iden.filtered.w.df <- hgt_gained_p_iden.filtered.w.df %>% left_join(hog_sp_prot_ct.df, by=c("species_name","hog_set")) %>% rename(hog_sp_prot_ct = count) %>% mutate(w_hog_set_ct = 1/hog_sp_prot_ct)     # Calculate weight for each HOG set to normalize counts by number of proteins in each HOG set

hgt_gained_p_iden.filtered.w.css.class.df <- hgt_gained_p_iden.filtered.w.df %>% mutate(class=ifelse(order=="Crocodylia", "Crocodylia", class), pooled_w = w_qseqid*w_hog_set_ct) %>% group_by(hog_set, class) %>% 
  summarise(total_w = sum(pooled_w), av_class_pident = mean(pident), av_class_qcovhsp = mean(qcovhsp)) %>% ungroup()
hgt_gained_p_iden.filtered.w.css.order.df <- hgt_gained_p_iden.filtered.w.df %>% mutate(class=ifelse(order=="Crocodylia", "Crocodylia", class), pooled_w = w_qseqid*w_hog_set_ct) %>% group_by(hog_set, class, order) %>% summarise(total_w = sum(pooled_w)) %>% ungroup()
hgt_gained_p_iden.filtered.w.css.family.df <- hgt_gained_p_iden.filtered.w.df %>% mutate(class=ifelse(order=="Crocodylia", "Crocodylia", class), pooled_w = w_qseqid*w_hog_set_ct) %>% group_by(hog_set, class, order, family) %>% summarise(total_w = sum(pooled_w)) %>% ungroup()

hgt_gained_p_iden.filtered.w.css.mj.df <- hgt_gained_p_iden.filtered.w.css.class.df %>% group_by(hog_set) %>% mutate(total_sum = sum(total_w), max_w = max(total_w), n_max = sum(total_w == max_w)) %>%
  # summarise(class = class[which.max(total_w)], av_class_pident = av_class_pident[which.max(total_w)], av_class_qcovhsp = av_class_qcovhsp[which.max(total_w)], .groups = "drop")
  summarise(class = ifelse(n_max[1] > 1 | max_w[1] / total_sum[1] <= 0.5, NA, class[which.max(total_w)]), av_class_pident = av_class_pident[which.max(total_w)], av_class_qcovhsp = av_class_qcovhsp[which.max(total_w)], .groups = "drop") %>%
  mutate(class = ifelse(class == "N/A", NA, class))

order_assign.df <- hgt_gained_p_iden.filtered.w.css.order.df %>% inner_join(hgt_gained_p_iden.filtered.w.css.mj.df %>% select(hog_set, class), by = c("hog_set", "class")) %>% group_by(hog_set, class) %>%
  mutate(total_sum = sum(total_w), max_w = max(total_w), n_max = sum(total_w == max_w)) %>% summarise(order = ifelse(is.na(class[1]) | n_max[1] > 1 | max_w[1] / total_sum[1] <= 0.5, NA, order[which.max(total_w)]), .groups = "drop")
hgt_gained_p_iden.filtered.w.css.mj.df <- hgt_gained_p_iden.filtered.w.css.mj.df %>% left_join(order_assign.df, by = c("hog_set", "class"))

family_assign.df <- hgt_gained_p_iden.filtered.w.css.family.df %>% inner_join(hgt_gained_p_iden.filtered.w.css.mj.df %>% select(hog_set, class, order), by = c("hog_set", "class", "order")) %>% group_by(hog_set, class, order) %>%
  mutate(total_sum = sum(total_w), max_w = max(total_w), n_max = sum(total_w == max_w)) %>% summarise(family = ifelse(is.na(class[1]) | is.na(order[1]) | n_max[1] > 1 | max_w[1] / total_sum[1] <= 0.5, NA, family[which.max(total_w)]), .groups = "drop")
hgt_gained_p_iden.filtered.w.css.mj.df <- hgt_gained_p_iden.filtered.w.css.mj.df %>% left_join(family_assign.df, by = c("hog_set", "class", "order"))

########## ***** Supplementary Table 3: Putative donors of HGT-derived ortholog groups *****
hgt_gained_p_iden.filtered.w.css.mj.f <- "260107_HGT_analysis/260107_HGT_donor_taxa_majority_rule.tsv"
# write.table(hgt_gained_p_iden.filtered.w.css.mj.df, file = hgt_gained_p_iden.filtered.w.css.mj.f, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

gained_hgt_og_set_majority.v.f <- "260107_HGT_analysis/260107_gained_hog_set_list_majority_rule.tsv"
# write.table(hgt_gained_p_iden.filtered.w.css.mj.df %>% filter(!is.na(class)) %>% pull(hog_set), gained_hgt_og_set_majority.v.f, col.names=F, row.names=F, sep="\n", quote=F)


########## Plot number of HGT events from putative donors of each taxonomic class in each branch
poxv_state_tr_all_m_hgt_split.df <- poxv_state_tr_all_m_hgt.df %>% mutate(og_num=lengths(str_split(hog_set, "\\s*,\\s*")))
poxv_state_tr_all_m_hgt_split.df <- poxv_state_tr_all_m_hgt_split.df %>% group_by(hog_set) %>% mutate(ca_node=ifelse(unique(og_num)>1 & length(node)>1, getMRCA(poxv_t.rt.rn, node), NA)) %>% ungroup()
poxv_state_tr_all_m_hgt_split.df <- poxv_state_tr_all_m_hgt_split.df %>% group_by(hog_set) %>% filter(!(n() > 1 & any(!is.na(ca_node))) | parent==ca_node) %>% ungroup()
poxv_state_tr_all_m_hgt_split.df <- poxv_state_tr_all_m_hgt_split.df %>% group_by(hog_set) %>% mutate(other_parent_is_anc = sapply(seq_len(n()), function(i) {any(parent[-i] %in% unlist(Ancestors(poxv_t.rt.rn, parent[i], type = "all")))})) %>% ungroup()
poxv_state_tr_all_m_hgt_split.df <- poxv_state_tr_all_m_hgt_split.df %>% filter(other_parent_is_anc==FALSE) %>% group_by(hog_set, parent, node) %>% summarise(num_transition=1) %>% ungroup()
poxv_state_tr_all_m_hgt_split.df <- poxv_state_tr_all_m_hgt_split.df %>% left_join(hgt_gained_p_iden.filtered.w.css.mj.df %>% select(hog_set, host_class = class, host_av_class_pident = av_class_pident, host_av_class_qcovhsp = av_class_qcovhsp, host_order = order, host_family = family), by = "hog_set")

poxv_state_tr_all_m_hgt_split.df %>% group_by(hog_set) %>% summarise(sum_transition=sum(num_transition), .groups = "drop") %>% filter(sum_transition > 1, hog_set %in% (hgt_gained_p_iden.filtered.w.css.mj.df %>% filter(!is.na(class)) %>% pull(hog_set)))

poxv_state_tr_all_m_hgt_split_pooled_name.df <- poxv_state_tr_all_m_hgt_split.df %>% mutate(host_name=paste(hog_set, host_class, host_order, host_family, sep="%")) %>% group_by(parent, node) %>% summarise(annot_name=paste(host_name, collapse=";"), .groups = "drop")

poxv_state_tr_all_m_hgt_split_ct.df <- poxv_state_tr_all_m_hgt_split.df %>% group_by(parent, node, host_class) %>% summarise(transition_count = n()) %>% ungroup()

poxv_state_tr_all_m_hgt_split_ct_list <- list()

for (a in (poxv_state_tr_all_m_hgt_split_ct.df %>% group_by(host_class) %>% summarise(sum_transition_count=sum(transition_count), .groups = "drop") %>% arrange(desc(sum_transition_count)))$host_class) {
  poxv_state_tr_all_m_total_hgt_split_ct.info.df <- poxv_t.rt.info.df %>% left_join(poxv_state_tr_all_m_hgt_split_ct.df %>% filter(host_class==a), by = c("parent", "node"))
  poxv_state_tr_all_m_total_hgt_split_ct.info.df <- poxv_state_tr_all_m_total_hgt_split_ct.info.df %>% mutate(transition_count=ifelse(is.na(transition_count), 0, transition_count))
  poxv_state_tr_all_m_total_hgt_split_ct.info.df <- poxv_state_tr_all_m_total_hgt_split_ct.info.df %>% left_join(poxv_state_tr_all_m_hgt_split_pooled_name.df, by = c("parent", "node"))
  poxv_state_tr_all_m_total_hgt_split_ct.info.df <- poxv_state_tr_all_m_total_hgt_split_ct.info.df %>% mutate(annot_name = map_chr(str_split(annot_name, ";"), function(x) {
    x <- x[str_split_i(x, "%", 2) == a]
    paste(x, collapse = ";")}))
  
  poxv_state_tr_all_m_hgt_split_ct.p <- ggtree(poxv_t.rt.rn, linewidth=2, aes(color=transition_count))
  poxv_state_tr_all_m_hgt_split_ct.p <- poxv_state_tr_all_m_hgt_split_ct.p %<+% poxv_state_tr_all_m_total_hgt_split_ct.info.df
  poxv_state_tr_all_m_hgt_split_ct.p <- poxv_state_tr_all_m_hgt_split_ct.p + theme_tree()
  poxv_state_tr_all_m_hgt_split_ct.p <- poxv_state_tr_all_m_hgt_split_ct.p + geom_text2(aes(subset=transition_count>0, label=paste(transition_count, "|", annot_name, sep="")), color="darkviolet", size=5, hjust=1.3)
  poxv_state_tr_all_m_hgt_split_ct.p <- poxv_state_tr_all_m_hgt_split_ct.p + scale_color_gradient("transition_count", name="Number of HGT genes", low="black", high="violet", breaks = seq(0, max(poxv_state_tr_all_m_hgt_split_ct.df$transition_count, na.rm = TRUE), by=1))
  poxv_state_tr_all_m_hgt_split_ct.p <- poxv_state_tr_all_m_hgt_split_ct.p + labs(title="HGT", subtitle=a)
  poxv_state_tr_all_m_hgt_split_ct.p <- poxv_state_tr_all_m_hgt_split_ct.p + theme(text = element_text(size = 20), plot.title = element_text(size = 30), plot.subtitle = element_text(size = 24))
  poxv_state_tr_all_m_hgt_split_ct.p <- poxv_state_tr_all_m_hgt_split_ct.p + xlim(0,4.5)
  
  poxv_state_tr_all_m_hgt_split_ct_list[[a]] <- poxv_state_tr_all_m_hgt_split_ct.p
}

########## ***** Figure 4C: Number of HGT events from putative donors of each taxonomic class in each branch *****
wrap_plots(poxv_state_tr_all_m_hgt_split_ct_list[c(1,3:7)])

# pdf("Figures/gene_origin_hgt_gene_by_class.pdf", width=60, height=20)
# wrap_plots(poxv_state_tr_all_m_hgt_split_ct_list[c(1,3:6,8)]) + plot_layout(nrow=1)
# dev.off()


########## Visualize proportions of HGT-derived ortholog groups
tax_rank.v <- c("phylum", "class", "order", "family", "genus", "species")
subfamily.v <- c("all","Chordopoxvirinae","Entomopoxvirinae")

tax_rank.ht_list <- list()
hgt_res_tax_test_list <- list()
hgt_res_tax_test_sep_host_list <- list()
hgt_res_tax_pooled_mat_ht_list <- list()

num_top_clade <- 3
norm_count_cutoff_ht <- 1

for (a in tax_rank.v) {
  for (b in subfamily.v) {
      hgt_gained_res_vnp.v <- sort(unique((hgt_gained_res.np.df %>% filter(sskingdoms=="Viruses"))[,a]))
      hgt_gained_res_nvnp.v <- sort(unique((hgt_gained_res.np.df %>% filter(sskingdoms!="Viruses"))[,a]))
      
      # ### Summarise proportion of HGT genes (qcovhsp > 80 or 90 still provide similar results) using majority rule
      # hgt_res_tax_ct.df <- hgt_gained_res.np.df %>% group_by(qseqid, sscinames) %>% slice_head(n=1) %>% ungroup()                                                                                                                                # Pick only best hit from a subject species
      # hgt_res_tax_ct.df <- hgt_res_tax_ct.df %>% filter(qcovhsp > 50) %>% group_by(species_name, hog_set, qseqid, !!sym(a)) %>% summarise(count = n(), .groups = "keep")                                                                         # Filter by % coverage, count hits for each query protein
      # hgt_res_tax_ct.df <- hgt_res_tax_ct.df %>% group_by(species_name, hog_set, qseqid) %>% mutate(total=sum(count), prot_norm_count = count/total)                                                                                             # Normalize counts for each query protein
      # hgt_res_tax_ct.df <- hgt_res_tax_ct.df %>% group_by(species_name, hog_set, !!sym(a)) %>% summarise(count = sum(prot_norm_count), .groups = "keep")                                                                                         # Sum normalized counts for each HOG set
      # hgt_res_tax_ct.df <- hgt_res_tax_ct.df %>% group_by(species_name, hog_set, !!sym(a)) %>% mutate(total_count = sum(count), max_count = max(count), n_max = sum(count == max_count)) %>% ungroup()
      # hgt_res_tax_ct.df <- hgt_res_tax_ct.df %>% group_by(species_name, hog_set) %>% summarise(!!sym(a) := if(first(n_max) > 1 || first(max_count)/first(total_count) <= 0.5) {NA} else {.data[[a]][which.max(total_count)]}, .groups = "drop")  # Apply the majority rule for each ortholog group in each species
      # hgt_res_tax_ct.df <- hgt_res_tax_ct.df %>% mutate(!!sym(a) := ifelse(!!sym(a) == "N/A", NA, !!sym(a)), count = 1)
      # hgt_res_tax_ct.df <- hgt_res_tax_ct.df %>% group_by(species_name, !!sym(a)) %>% summarise(tax_count = sum(count), .groups = "keep")                                                                                                        # Sum counts for each virus
      # hgt_res_tax_ct.df <- hgt_res_tax_ct.df %>% mutate(og_num=hog_ct_per_sp.df$og_count[match(species_name, hog_ct_per_sp.df$species_name)], tax_prop = tax_count/og_num) %>% ungroup()                                                         # Calculate proportions of genes likely originated from non-poxvirus (divided by total number of HOG)
      # hgt_res_tax_ct.df$species_name <- factor(hgt_res_tax_ct.df$species_name, levels = levels(poxv_genome_stats.df$species_name))
  
      ### Summarise proportion of HGT genes (qcovhsp > 80 or 90 still provide similar results)
      hgt_res_tax_ct.df <- hgt_gained_res.np.df %>% group_by(qseqid, sscinames) %>% slice_head(n=1) %>% ungroup()                                                                                                 # Pick only best hit from a subject species
      hgt_res_tax_ct.df <- hgt_res_tax_ct.df %>% filter(qcovhsp > 50) %>% group_by(species_name, hog_set, qseqid, !!sym(a)) %>% summarise(count = n(), .groups = "keep")                                          # Filter by % coverage, count hits for each query protein
      hgt_res_tax_ct.df <- hgt_res_tax_ct.df %>% group_by(species_name, hog_set, qseqid) %>% mutate(total=sum(count), prot_norm_count = count/total)                                                              # Normalize counts for each query protein
      hgt_res_tax_ct.df <- hgt_res_tax_ct.df %>% group_by(species_name, hog_set, !!sym(a)) %>% summarise(count = sum(prot_norm_count), .groups = "keep")                                                          # Sum normalized counts for each HOG set
      hgt_res_tax_ct.df <- hgt_res_tax_ct.df %>% left_join(hog_sp_prot_ct.df, by=c("species_name","hog_set"), suffix = c("", "_total")) %>% rename(total = count_total) %>% mutate(og_norm_count = count/total)   # Normalize counts for each HOG
      hgt_res_tax_ct.df <- hgt_res_tax_ct.df %>% group_by(species_name, !!sym(a)) %>% summarise(tax_count = sum(og_norm_count), .groups = "keep")                                                                 # Sum normalized counts for each virus
      hgt_res_tax_ct.df <- hgt_res_tax_ct.df %>% mutate(og_num=hog_ct_per_sp.df$og_count[match(species_name, hog_ct_per_sp.df$species_name)], tax_prop = tax_count/og_num) %>% ungroup()                          # Calculate proportions of genes likely originated from non-poxvirus (divided by total number of HOG)
      hgt_res_tax_ct.df$species_name <- factor(hgt_res_tax_ct.df$species_name, levels = levels(poxv_genome_stats.df$species_name))

      if (b=="all") {
        hgt_res_tax_ct_m_box.df <- hgt_res_tax_ct.df %>% filter(!!sym(a) != "N/A") %>% rowwise() %>% 
          mutate(host_status = ifelse(any(virushostdb_poxv_long.df$virus_name == species_name & virushostdb_poxv_long.df$host_lineage == !!sym(a)), !!sym(a), "Mismatch")) %>%
          group_by(species_name, host_status) %>% summarise(av_tax_prop = mean(tax_prop)) %>% ungroup()
        
        ### Test difference in proportions of putative HGT genes by host match status separated by donor taxa 
        print(paste("Test at", a, "level", sep=" "))
        
        hgt_res_tax_ct_m_box.p <- ggplot(hgt_res_tax_ct_m_box.df, aes(x=host_status, y=av_tax_prop, fill=host_status))
        hgt_res_tax_ct_m_box.p <- hgt_res_tax_ct_m_box.p + geom_jitter(alpha=0.25, shape=16, size=3, position=position_jitter(0))
        hgt_res_tax_ct_m_box.p <- hgt_res_tax_ct_m_box.p + geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.75)
        hgt_res_tax_ct_m_box.p <- hgt_res_tax_ct_m_box.p + guides(fill="none")
        hgt_res_tax_ct_m_box.p <- hgt_res_tax_ct_m_box.p + theme_classic()
        hgt_res_tax_ct_m_box.p <- hgt_res_tax_ct_m_box.p + theme(text = element_text(size = 24), plot.title = element_text(size = 30))
        hgt_res_tax_ct_m_box.p <- hgt_res_tax_ct_m_box.p + labs(title=a, x="Host status", y="Proportion of genes likely\noriginating from HGT")
        # hgt_res_tax_ct_m_box.p <- hgt_res_tax_ct_m_box.p + scale_y_continuous(breaks=seq(-0.02, 0.18, by=0.02), limits=c(-0.02,0.18))
        
        if (length(table((hgt_res_tax_ct_m_box.df %>% group_by(host_status))$host_status))!=1) {
          # print(hgt_res_tax_ct_m_box.df %>% group_by(host_status) %>% shapiro_test(av_tax_prop))
          # print(hgt_res_tax_ct_m_box.df %>% levene_test(av_tax_prop~host_status))
          
          ### Use Kruskal-Wallis test since normality assumption is violated
          hgt_res_tax_ct_m_box_res.kruskal <- hgt_res_tax_ct_m_box.df %>% kruskal_test(av_tax_prop~host_status)
          print(hgt_res_tax_ct_m_box_res.kruskal)
          
          hgt_res_tax_ct_m_box_stats <- hgt_res_tax_ct_m_box.df %>% wilcox_test(av_tax_prop~host_status, p.adjust.method="fdr")
          hgt_res_tax_ct_m_box_stats <- hgt_res_tax_ct_m_box_stats %>% filter(group1=="Mismatch" | group2=="Mismatch")
          hgt_res_tax_ct_m_box_stats <- hgt_res_tax_ct_m_box_stats %>% add_xy_position(x="host_status", fun="max")

          hgt_res_tax_ct_m_box.p <- hgt_res_tax_ct_m_box.p + labs(subtitle=paste("p = ", hgt_res_tax_ct_m_box_res.kruskal$p))
          hgt_res_tax_ct_m_box.p <- hgt_res_tax_ct_m_box.p + geom_line(aes(group = species_name), color = "gray40", alpha = 0.35)
          hgt_res_tax_ct_m_box.p <- hgt_res_tax_ct_m_box.p + stat_pvalue_manual(hgt_res_tax_ct_m_box_stats, label="p", hide.ns=TRUE, inherit.aes=FALSE, size=8)
        }
        
        hgt_res_tax_test_sep_host_list[[a]] <- hgt_res_tax_ct_m_box.p
      }
      
      if (b=="all") {
        hgt_res_tax_ct_m_box.df <- hgt_res_tax_ct.df %>% filter(!!sym(a) != "N/A") %>% rowwise() %>% 
          mutate(host_status = ifelse(any(virushostdb_poxv_long.df$virus_name == species_name & virushostdb_poxv_long.df$host_lineage == !!sym(a)), "Match","Mismatch")) %>%
          group_by(species_name, host_status) %>% summarise(av_tax_prop = mean(tax_prop)) %>% ungroup()
        
        ### Test difference in proportions of putative HGT genes separated by host match status
        print(paste("Test at", a, "level", sep=" "))
        
        hgt_res_tax_ct_m_box.p <- ggplot(hgt_res_tax_ct_m_box.df, aes(x=host_status, y=av_tax_prop, fill=host_status))
        hgt_res_tax_ct_m_box.p <- hgt_res_tax_ct_m_box.p + geom_jitter(alpha=0.25, shape=16, size=3, position=position_jitter(0))
        hgt_res_tax_ct_m_box.p <- hgt_res_tax_ct_m_box.p + geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.75)
        hgt_res_tax_ct_m_box.p <- hgt_res_tax_ct_m_box.p + guides(fill="none")
        hgt_res_tax_ct_m_box.p <- hgt_res_tax_ct_m_box.p + theme_classic()
        hgt_res_tax_ct_m_box.p <- hgt_res_tax_ct_m_box.p + theme(text = element_text(size = 24), plot.title = element_text(size = 30))
        hgt_res_tax_ct_m_box.p <- hgt_res_tax_ct_m_box.p + labs(title=a, x="Host status", y="Proportion of genes likely\noriginating from HGT")
        hgt_res_tax_ct_m_box.p <- hgt_res_tax_ct_m_box.p + scale_y_continuous(breaks=seq(-0.02, 0.18, by=0.02), limits=c(-0.02,0.18))
        
        if (length(table((hgt_res_tax_ct_m_box.df %>% group_by(host_status))$host_status))!=1) {
          print(hgt_res_tax_ct_m_box.df %>% group_by(host_status) %>% shapiro_test(av_tax_prop))
          print(hgt_res_tax_ct_m_box.df %>% levene_test(av_tax_prop~host_status))
          
          ### Use Kruskal-Wallis test since normality assumption is violated
          hgt_res_tax_ct_m_box_res.kruskal <- hgt_res_tax_ct_m_box.df %>% kruskal_test(av_tax_prop~host_status)
          print(hgt_res_tax_ct_m_box_res.kruskal)
          
          # hgt_res_tax_ct_m_box_stats <- hgt_res_tax_ct_m_box.df %>% wilcox_test(av_tax_prop~host_status, p.adjust.method="fdr")
          # hgt_res_tax_ct_m_box_stats <- hgt_res_tax_ct_m_box_stats %>% add_xy_position(x="host_status", fun="max")
          # print(hgt_res_tax_ct_m_box_stats)
          
          hgt_res_tax_ct_m_box.p <- hgt_res_tax_ct_m_box.p + labs(subtitle=paste("p = ", hgt_res_tax_ct_m_box_res.kruskal$p))
          hgt_res_tax_ct_m_box.p <- hgt_res_tax_ct_m_box.p + geom_line(aes(group = species_name), color = "gray40", alpha = 0.35)
          # hgt_res_tax_ct_m_box.p <- hgt_res_tax_ct_m_box.p + stat_pvalue_manual(hgt_res_tax_ct_m_box_stats, label="p", hide.ns=FALSE, inherit.aes=FALSE, size=8)
        }
        
        hgt_res_tax_test_list[[a]] <- hgt_res_tax_ct_m_box.p
      }
      
      hgt_res_tax_mat <- as.matrix(xtabs(reformulate(c("species_name", a), response = "tax_prop"), data = hgt_res_tax_ct.df))
      if (b!="all") {hgt_res_tax_mat <- hgt_res_tax_mat[rownames(hgt_res_tax_mat) %in% (poxv_ncbi_ds.df %>% filter(Subfamily==b))$tip_label_dummy, , drop = FALSE]}

      ## Vector of clades for each row of hgt_res_tax_mat
      row_clade <- poxv_genome_stats.df$clade[match(rownames(hgt_res_tax_mat), poxv_genome_stats.df$species_name)]
      stopifnot(!any(is.na(row_clade)))
      
      ## Pool by clade (sum across viruses in the same clade)
      hgt_res_tax_mat_clade <- rowsum(as.matrix(hgt_res_tax_mat), group = row_clade, reorder = FALSE)
      
      # Get top n columns per row for both non-poxvirus viral and non-viral sources (filtered by num_top_clade)
      top_cols <- unique(unlist(apply(hgt_res_tax_mat_clade, 1, function(x) {
        valid_idx_vnp <- which(colnames(hgt_res_tax_mat_clade) %in% hgt_gained_res_vnp.v & colnames(hgt_res_tax_mat_clade) != "N/A")
        x_valid_vnp <- x[valid_idx_vnp]
        ord_vnp <- order(x_valid_vnp, decreasing = TRUE, na.last = NA)
        ord_vnp <- ord_vnp[x_valid_vnp[ord_vnp] > 0]
        top_vnp <- head(ord_vnp, num_top_clade)

        valid_idx_nvnp <- which(colnames(hgt_res_tax_mat_clade) %in% hgt_gained_res_nvnp.v & colnames(hgt_res_tax_mat_clade) != "N/A")
        x_valid_nvnp <- x[valid_idx_nvnp]
        ord_nvnp <- order(x_valid_nvnp, decreasing = TRUE, na.last = NA)
        ord_nvnp <- ord_nvnp[x_valid_nvnp[ord_nvnp] > 0]
        top_nvnp <- head(ord_nvnp, num_top_clade)

        c(colnames(hgt_res_tax_mat_clade)[valid_idx_vnp[top_vnp]], colnames(hgt_res_tax_mat_clade)[valid_idx_nvnp[top_nvnp]])
      })))

      # # Get top columns for showing in the heatmap (filtered by norm_count_cutoff_ht)
      # hgt_gained_ct4ht.df <- hgt_gained_res.np.df %>% group_by(qseqid, sscinames) %>% slice_head(n=1) %>% ungroup() %>% filter(qcovhsp > 50)                                                             # Pick only best hit from a subject species, filter by % coverage
      # hgt_gained_ct4ht.df <- hgt_gained_ct4ht.df %>% mutate(clade=poxv_genome_stats.df$clade[match(species_name, poxv_genome_stats.df$species_name)])
      # hgt_gained_ct4ht.df <- hgt_gained_ct4ht.df %>% mutate(Subfamily=poxv_genome_stats.df$Subfamily[match(species_name, poxv_genome_stats.df$species_name)])
      # 
      # hgt_gained_ct4ht.w.df <- hgt_gained_ct4ht.df %>% group_by(qseqid) %>% mutate(w_qseqid = 1/n()) %>% ungroup()                                                                                       # Calculate weight for each query protein to normalize counts by number of hits for each query protein
      # hgt_gained_ct4ht.w.df <- hgt_gained_ct4ht.w.df %>% left_join(hog_sp_prot_ct.df, by=c("species_name","hog_set")) %>% rename(hog_sp_prot_ct = count) %>% mutate(w_hog_set_ct = 1/hog_sp_prot_ct)     # Calculate weight for each HOG set to normalize counts by number of proteins in each HOG set
      # 
      # hgt_gained_ct4ht.w.subfam.df <- hgt_gained_ct4ht.w.df %>% left_join(hog_subfam_sp_ct.df, by=c("Subfamily","hog_set")) %>% rename(subfam_sp_ct = count) %>% mutate(w_subfam_sp_ct = 1/subfam_sp_ct) # Calculate weight for each subfamily to normalize counts by number of species shared in each subfamily
      # 
      # hgt_gained_ct4ht.w.subfam_cdpx_ct.df <- hgt_gained_ct4ht.w.subfam.df %>% filter(Subfamily=="Chordopoxvirinae") %>% select(qseqid, sym(a), hog_set, w_qseqid, w_hog_set_ct, w_subfam_sp_ct)
      # hgt_gained_ct4ht.w.subfam_cdpx_ct.df <- hgt_gained_ct4ht.w.subfam_cdpx_ct.df %>% group_by(!!sym(a)) %>% summarise(count = sum(w_qseqid*w_hog_set_ct*w_subfam_sp_ct)) %>% ungroup() %>% filter(!(!!sym(a))=="N/A") %>% arrange(desc(count))
      # 
      # hgt_gained_ct4ht.w.subfam_etpx_ct.df <- hgt_gained_ct4ht.w.subfam.df %>% filter(Subfamily=="Entomopoxvirinae") %>% select(qseqid, sym(a), hog_set, w_qseqid, w_hog_set_ct, w_subfam_sp_ct)
      # hgt_gained_ct4ht.w.subfam_etpx_ct.df <- hgt_gained_ct4ht.w.subfam_etpx_ct.df %>% group_by(!!sym(a)) %>% summarise(count = sum(w_qseqid*w_hog_set_ct*w_subfam_sp_ct)) %>% ungroup() %>% filter(!(!!sym(a))=="N/A") %>% arrange(desc(count))
      # 
      # top_cols <- union((hgt_gained_ct4ht.w.subfam_cdpx_ct.df %>% filter(count > norm_count_cutoff_ht) %>% pull(!!sym(a))), (hgt_gained_ct4ht.w.subfam_etpx_ct.df %>% filter(count > norm_count_cutoff_ht) %>% pull(!!sym(a)))) # For heatmap

      # Subset to valid and unique columns only
      host_lineages_subset <- virushostdb_poxv_long.df$host_lineage[virushostdb_poxv_long.df$virus_name %in% rownames(hgt_res_tax_mat)]
      top_cols_wt_host <- intersect(colnames(hgt_res_tax_mat), union(top_cols, sort(unique(host_lineages_subset))))
      
      # Subset matrices
      hgt_res_tax_mat <- hgt_res_tax_mat[, top_cols_wt_host, drop = FALSE]

      hgt_res_tax_vnp.v <- colnames(hgt_res_tax_mat)[colnames(hgt_res_tax_mat) %in% hgt_gained_res_vnp.v]
      hgt_res_tax_vnp.v <- hgt_res_tax_vnp.v[hgt_res_tax_vnp.v!="N/A"]
      hgt_res_tax_nvnp.v <- colnames(hgt_res_tax_mat)[colnames(hgt_res_tax_mat) %in% hgt_gained_res_nvnp.v]
      hgt_res_tax_nvnp.v <- hgt_res_tax_nvnp.v[hgt_res_tax_nvnp.v!="N/A"]
      
      hgt_res_tax_mat <- hgt_res_tax_mat[,c(hgt_res_tax_nvnp.v, hgt_res_tax_vnp.v)]
      
      # Define cell annotation function for top n
      cell_label_fun <- function(j, i, x, y, width, height, fill) {
        poxv_name  <- rownames(hgt_res_tax_mat)[i]
        row_values <- hgt_res_tax_mat[i, ]
        
        n_nvnp  <- length(hgt_res_tax_nvnp.v)
        n_total <- ncol(hgt_res_tax_mat)
        
        ## ---- Partition 1 ----
        idx1 <- order(row_values[1:n_nvnp], decreasing = TRUE)
        top1 <- head(idx1[row_values[1:n_nvnp][idx1] > 0], num_top_clade)
        labs1 <- as.character(seq_along(top1))
        
        ## ---- Partition 2 ----
        idx2 <- order(row_values[(n_nvnp+1):n_total], decreasing = TRUE)
        top2 <- head(idx2[row_values[(n_nvnp+1):n_total][idx2] > 0], num_top_clade) + n_nvnp
        labs2 <- as.character(seq_along(top2))
        
        ## ---- Combine ----
        top_indices <- c(top1, top2)
        rank_labels  <- c(labs1, labs2)
        
        ## ---- Default: no label ----
        label <- ""
        
        ## ---- Add rank if top n ----
        # if (j %in% top_indices) {label <- rank_labels[which(top_indices == j)]}
        
        ## ---- Add "*" if lineage matches ----
        col_name <- colnames(hgt_res_tax_mat)[j]
        lineages <- virushostdb_poxv_long.df$host_lineage[virushostdb_poxv_long.df$virus_name == poxv_name]
        if (col_name %in% lineages) {label <- "*"}
        # if (col_name %in% lineages) {label <- paste0(label, "*")}
        if (label != "") {grid.text(label, x = x, y = y, vjust = 0.75, gp = gpar(col = "black", fontsize = 24))}
      }
      
      # Tip order as plotted
      tip_order <- poxv_t.p$data %>% filter(isTip) %>% arrange(desc(y)) %>% pull(label)
      
      if (b != "all") {
        keep <- (poxv_ncbi_ds.df %>% filter(Subfamily == b))$tip_label_dummy
        tip_order <- tip_order[tip_order %in% keep]
      }
      
      # Z-score by each column (i.e., across viruses/species)
      hgt_res_tax_mat_z <- apply(hgt_res_tax_mat, 2, function(x) {
        if (sd(x) == 0) {
          rep(0, length(x))  # Avoid NaN if all values are the same
        } else {
          (x - mean(x)) / sd(x)
        }
      })
      rownames(hgt_res_tax_mat_z) <- rownames(hgt_res_tax_mat)
      
      hgt_res_tax_mat <- hgt_res_tax_mat[tip_order, , drop = FALSE]
      
      # Extract clade information in same order and proportion of genes likely derived by HGT
      row_clade <- poxv_genome_stats.df[match(rownames(hgt_res_tax_mat), poxv_genome_stats.df$species_name), "clade"]
      names(row_clade) <- rownames(hgt_res_tax_mat)
      
      row_prop_hgt <- poxv_genome_stats.df[match(rownames(hgt_res_tax_mat), poxv_genome_stats.df$species_name), "hgt_likely_number"]/poxv_genome_stats.df[match(rownames(hgt_res_tax_mat), poxv_genome_stats.df$species_name), "hog_number"]
      names(row_prop_hgt) <- rownames(hgt_res_tax_mat)
      
      max_row_prot_hgt <- max(row_prop_hgt)
      
      if (b=="all") {
        clade_colors <- poxv_color.clade.v
      } else {
        row_clade <- row_clade[row_clade %in% (poxv_ncbi_ds.df %>% filter(Subfamily==b))$clade]
        clade_colors <- poxv_color.clade.v[names(poxv_color.clade.v) %in% (poxv_ncbi_ds.df %>% filter(Subfamily==b))$clade]
        
        row_prop_hgt <- row_prop_hgt[names(row_prop_hgt) %in% (poxv_ncbi_ds.df %>% filter(Subfamily==b))$tip_label_dummy]
      }
      
      # Create annotation
      row_ha <- rowAnnotation(clade = row_clade, prop_hgt = row_prop_hgt, col = list(clade = clade_colors, prop_hgt = colorRamp2(c(0, max_row_prot_hgt), c("white", "navyblue"))), show_annotation_name = TRUE)
      
      # Create heatmap
      if (b=="all") {
        if (a=="class") {class_max <- max(hgt_res_tax_mat)}
        if (a=="order") {
          ht <- Heatmap(hgt_res_tax_mat, name = paste("PropNormCt\n(", a, ")", sep = ""), cluster_rows = FALSE, cluster_columns = FALSE, row_names_side = "left", column_names_side = "top",
                        column_split=c(rep("Non-virus", length(hgt_res_tax_nvnp.v)),rep("Virus", length(hgt_res_tax_vnp.v))), row_gap = unit(5, "mm"), column_gap = unit(5, "mm"),
                        col = colorRamp2(c(0, class_max), c("white", "red")), heatmap_legend_param = list(title = "Proportion of genes likely\nhaving a non-poxviral origin"), cell_fun = cell_label_fun, right_annotation = row_ha)      
        } else {
          # Create order-rank plot using maximum value of the class-rank plot for extracting Crocodylia taxon
          ht <- Heatmap(hgt_res_tax_mat, name = paste("PropNormCt\n(", a, ")", sep = ""), cluster_rows = FALSE, cluster_columns = FALSE, row_names_side = "left", column_names_side = "top",
                        column_split=c(rep("Non-virus", length(hgt_res_tax_nvnp.v)),rep("Virus", length(hgt_res_tax_vnp.v))), row_gap = unit(5, "mm"), column_gap = unit(5, "mm"),
                        col = colorRamp2(c(0, max(hgt_res_tax_mat)), c("white", "red")), heatmap_legend_param = list(title = "Proportion of genes likely\nhaving a non-poxviral origin"), cell_fun = cell_label_fun, right_annotation = row_ha)      
        }
      } else {
        ht <- Heatmap(hgt_res_tax_mat, name = paste("PropNormCt\n(", a, ")", sep = ""), cluster_rows = FALSE, cluster_columns = FALSE, row_names_side = "left", column_names_side = "top",
                      column_split=c(rep("Non-virus", length(hgt_res_tax_nvnp.v)),rep("Virus", length(hgt_res_tax_vnp.v))), row_gap = unit(5, "mm"), column_gap = unit(5, "mm"),
                      col = colorRamp2(c(0, max(hgt_res_tax_mat)), c("white", "red")), heatmap_legend_param = list(title = "Proportion of genes likely\nhaving a non-poxviral origin"), cell_fun = cell_label_fun, right_annotation = row_ha)        
      }
      
      # Create a pooled heatmap at class level
      if (b=="all") {
        if (a %in% c("class","order")) {
          ##### Extract host information for adding in the phylogenetic tree
          virushostdb_poxv_long_summarized.df <- virus_host_mat %>% as.data.frame() %>% rownames_to_column("virus") %>% pivot_longer(cols = -virus, names_to = "host_taxa", values_to = "presence") %>% filter(presence=="1")
          virushostdb_poxv_long_summarized.df <- virushostdb_poxv_long_summarized.df %>% mutate(host_taxa_gene_num=poxv_genome_stats.df$hog_number[match(virus, poxv_genome_stats.df$species_name)])
          
          virushostdb_poxv_long_summarized_gene_ct.df <- virushostdb_poxv_long_summarized.df %>% group_by(host_taxa) %>% summarise(sum_host_taxa_gene_num=sum(host_taxa_gene_num)) %>% ungroup()
          
          hgt_res_tax_ct_pooled.df <- hgt_res_tax_ct.df %>% mutate(host_taxa=virushostdb_poxv_long_summarized.df$host_taxa[match(species_name, virushostdb_poxv_long_summarized.df$virus)])
          # hgt_res_tax_ct_pooled.df <- hgt_res_tax_ct_pooled.df %>% group_by(host_taxa, !!sym(a)) %>% summarise(sum_tax_count=sum(tax_count)) %>% ungroup() 
          # hgt_res_tax_ct_pooled.df <- hgt_res_tax_ct_pooled.df %>% mutate(sum_og_num=virushostdb_poxv_long_summarized_gene_ct.df$sum_host_taxa_gene_num[match(host_taxa, virushostdb_poxv_long_summarized_gene_ct.df$host_taxa)], sum_tax_prop=sum_tax_count/sum_og_num)
          hgt_res_tax_ct_pooled.df <- hgt_res_tax_ct_pooled.df %>% group_by(host_taxa) %>% complete(species_name, !!sym(a), fill=list(tax_count=0, tax_prop=0)) %>% ungroup()
          hgt_res_tax_ct_pooled.df <- hgt_res_tax_ct_pooled.df %>% group_by(host_taxa, !!sym(a)) %>% summarise(mean_tax_prop=mean(tax_prop, na.rm=TRUE)) %>% ungroup() 
          
          hgt_res_tax_pooled_mat <- as.matrix(xtabs(reformulate(c("host_taxa", a), response = "mean_tax_prop"), data = hgt_res_tax_ct_pooled.df))
          hgt_res_tax_pooled_mat <- hgt_res_tax_pooled_mat[, sort(colnames(virus_host_mat)[colnames(virus_host_mat) %in% colnames(hgt_res_tax_pooled_mat)]), drop=FALSE]
          
          if (a=="class") {
            hgt_res_tax_pooled_class_max <- max(hgt_res_tax_pooled_mat)
            hgt_res_tax_ct_class.df <- hgt_res_tax_ct.df

            hgt_res_tax_ct_class_pooled.df <- hgt_res_tax_ct_class.df %>% mutate(host_taxa=virushostdb_poxv_long_summarized.df$host_taxa[match(species_name, virushostdb_poxv_long_summarized.df$virus)])
            # hgt_res_tax_ct_class_pooled.df <- hgt_res_tax_ct_class_pooled.df %>% group_by(host_taxa, !!sym(a)) %>% summarise(sum_tax_count=sum(tax_count)) %>% ungroup()
            # hgt_res_tax_ct_class_pooled.df <- hgt_res_tax_ct_class_pooled.df %>% mutate(sum_og_num=virushostdb_poxv_long_summarized_gene_ct.df$sum_host_taxa_gene_num[match(host_taxa, virushostdb_poxv_long_summarized_gene_ct.df$host_taxa)], sum_tax_prop=sum_tax_count/sum_og_num)
            hgt_res_tax_ct_class_pooled.df <- hgt_res_tax_ct_class_pooled.df %>% group_by(host_taxa) %>% complete(species_name, !!sym(a), fill=list(tax_count=0, tax_prop=0)) %>% ungroup()
            hgt_res_tax_ct_class_pooled.df <- hgt_res_tax_ct_class_pooled.df %>% group_by(host_taxa, !!sym(a)) %>% summarise(mean_tax_prop=mean(tax_prop, na.rm=TRUE)) %>% ungroup() 

            hgt_res_tax_ct_class_pooled.df %>% filter(class %in% colnames(virus_host_mat))
          }
          
          hgt_res_tax_pooled_mat_ht <- Heatmap(hgt_res_tax_pooled_mat, cluster_rows = FALSE, cluster_columns = FALSE, row_names_side = "left", column_names_side = "top", row_gap = unit(5, "mm"), column_gap = unit(5, "mm"),
                                               row_title="Host taxa", column_title="Donor taxa", col = colorRamp2(c(0, hgt_res_tax_pooled_class_max), c("white", "red")), heatmap_legend_param = list(title = "Proportion of HGT-derived genes"))
          hgt_res_tax_pooled_mat_ht_list[[paste(a, "_", b, sep="")]] <- hgt_res_tax_pooled_mat_ht
        }
      }
      
      # Ensure row order is consistent between matrix and annotation
      hgt_res_tax_mat_z <- hgt_res_tax_mat_z[tip_order, , drop = FALSE]
      common_species <- intersect(rownames(hgt_res_tax_mat_z), poxv_genome_stats.df$species_name)
      hgt_res_tax_mat_z <- hgt_res_tax_mat_z[common_species, , drop=FALSE]
      if (b!="all") {hgt_res_tax_mat_z <- hgt_res_tax_mat_z[rownames(hgt_res_tax_mat_z) %in% (poxv_ncbi_ds.df %>% filter(Subfamily==b))$tip_label_dummy, , drop = FALSE]}
      
      # Extract clade information in same order
      row_clade <- poxv_genome_stats.df[match(rownames(hgt_res_tax_mat_z), poxv_genome_stats.df$species_name), "clade"]
      names(row_clade) <- rownames(hgt_res_tax_mat_z)
      if (b=="all") {
        clade_colors <- poxv_color.clade.v
        } else {
        row_clade <- row_clade[row_clade %in% (poxv_ncbi_ds.df %>% filter(Subfamily==b))$clade]
        clade_colors <- poxv_color.clade.v[names(poxv_color.clade.v) %in% (poxv_ncbi_ds.df %>% filter(Subfamily==b))$clade]
      }
      
      # Create annotation
      row_ha_z <- rowAnnotation(clade = row_clade, col = list(clade = clade_colors), show_annotation_name = TRUE)
      
      # Heatmap with top n labels in cells
      if (b=="all") {
        ht_z <- Heatmap(hgt_res_tax_mat_z, name = "Z-score", cluster_rows = FALSE, cluster_columns = FALSE, row_names_side = "left", column_names_side = "top", col = colorRamp2(c(-3, 0, 3), c("blue", "white", "red")),
                        column_split=c(rep("Non-virus", length(hgt_res_tax_nvnp.v)),rep("Virus", length(hgt_res_tax_vnp.v))), row_gap = unit(5, "mm"), column_gap = unit(5, "mm"),
                        heatmap_legend_param = list(title = "Z-score"), right_annotation = row_ha_z)        
      } else {
        ht_z <- Heatmap(hgt_res_tax_mat_z, name = "Z-score", cluster_rows = FALSE, cluster_columns = FALSE, row_names_side = "left", column_names_side = "top", col = colorRamp2(c(-3, 0, 3), c("blue", "white", "red")),
                        column_split=c(rep("Non-virus", length(hgt_res_tax_nvnp.v)),rep("Virus", length(hgt_res_tax_vnp.v))), column_gap = unit(5, "mm"),
                        heatmap_legend_param = list(title = "Z-score"), right_annotation = row_ha_z)
      }
            
      tax_rank.ht_list[[paste(a, "_", b, sep="")]] <- ht
      tax_rank.ht_list[[paste(a, "_", b, "_z", sep="")]] <- ht_z
  }
}

##### Combine desired heatmaps vertically
# draw(tax_rank.ht_list[["phylum_all"]] + tax_rank.ht_list[["phylum_all_z"]], padding = unit(c(5, 40, 5, 5), "mm"))
# draw(tax_rank.ht_list[["class_all"]] + tax_rank.ht_list[["class_all_z"]], padding = unit(c(5, 40, 5, 5), "mm"))
# draw(tax_rank.ht_list[["order_all"]] + tax_rank.ht_list[["order_all_z"]], padding = unit(c(5, 40, 5, 5), "mm"))

########## ***** Supplementary Figure 4: Proportions of HGT-derived orthogroups from each putative donor in each poxvirus genome *****
# pdf(paste("Figures/gene_origin_heatmap_top", num_top_clade, ".pdf", sep=""), width=24, height=12)
# pdf(paste("Figures/gene_origin_heatmap_ct_cutoff", norm_count_cutoff_ht, ".pdf", sep=""), width=24, height=12) # Use this for the publication
# draw(tax_rank.ht_list[["class_all"]] + tax_rank.ht_list[["order_all"]], padding = unit(c(5, 40, 5, 5), "mm"))
# dev.off()

########## ***** Figure 4B: Proportions of HGT-derived orthogroups stratified by whether the putative donor taxon matched the documented host taxon at different taxonomic ranks *****
# pdf("Figures/gene_origin_host_match_comparison.pdf", width=12, height=20)
# print(wrap_plots(hgt_res_tax_test_list))
# dev.off()

########## Proportions of HGT-derived orthogroups from each putative donor, showing whether the putative donor taxon matched the documented host taxon at different taxonomic ranks *****
# pdf("Figures/gene_origin_host_match_comparison_host_sep.pdf", width=16, height=20)
# print(wrap_plots(hgt_res_tax_test_sep_host_list))
# dev.off()

########## ***** Figure 4A: Host–donor taxonomic associations of HGT-derived ortholog groups *****
# pdf("Figures/gene_origin_pooled_heatmap.pdf", width=12, height=10)
# draw(hgt_res_tax_pooled_mat_ht_list[["class_all"]] + hgt_res_tax_pooled_mat_ht_list[["order_all"]])
# dev.off()



############################## Analyze correlation between virus and host genomic parameters and that between their contrasts ##############################
mean_zero_rule <- function(x){
  x <- x[!is.na(x)]        # remove NA
  if(length(x) == 0) return(0)   # only NA case
  mean(x)
}

poxv_virus_host.df <- poxv_ncbi_ds.virus_host_db.df %>%
  group_by(Organism.Name) %>%
  summarise(
    av_host_gene_num = mean_zero_rule(host_gene_num),
    av_host_genome_length = mean_zero_rule(as.numeric(host_genome_length)),
    av_host_genome_length_if = mean_zero_rule(as.numeric(host_genome_length_if)),
    av_host_total_gene_length = mean_zero_rule(host_total_gene_length),
    av_host_total_cds_length = mean_zero_rule(host_total_cds_length)
  ) %>%
  ungroup()

poxv_virus_host.df <- poxv_virus_host.df %>% mutate(av_host_gene_density=av_host_gene_num/av_host_genome_length_if,
                                                    av_host_prop_gene_length=av_host_total_gene_length/av_host_genome_length_if,
                                                    av_host_prop_cds_length=av_host_total_cds_length/av_host_genome_length_if)
poxv_virus_host.df <- poxv_virus_host.df %>% mutate(clade=poxv_ncbi_ds.df[match(Organism.Name, poxv_ncbi_ds.df$Organism.Name), "clade"])
poxv_virus_host.df <- poxv_virus_host.df %>% mutate(Organism.Name=gsub(" ", "_", Organism.Name))
poxv_virus_host.df <- poxv_virus_host.df %>% mutate(virus_gene_num=poxv_genome_stats.df[match(Organism.Name, gsub(" ", "_", poxv_genome_stats.df$species_name)), "gene_number"],
                                                    virus_genome_length=poxv_genome_stats.df[match(Organism.Name, gsub(" ", "_", poxv_genome_stats.df$species_name)), "genome_length"],
                                                    virus_total_gene_length=poxv_genome_stats.df[match(Organism.Name, gsub(" ", "_", poxv_genome_stats.df$species_name)), "total_gene_length"],
                                                    virus_total_cds_length=poxv_genome_stats.df[match(Organism.Name, gsub(" ", "_", poxv_genome_stats.df$species_name)), "total_cds_length"])
poxv_virus_host.df <- poxv_virus_host.df %>% mutate(virus_gene_num=as.numeric(virus_gene_num), virus_genome_length=as.numeric(virus_genome_length), virus_total_gene_length=as.numeric(virus_total_gene_length), virus_total_cds_length=as.numeric(virus_total_cds_length),
                                                    virus_og_count_pooled=as.numeric(hog_ct_per_sp.df$og_count_pooled[match(Organism.Name, hog_ct_per_sp.df$species_name)]),
                                                    virus_gained_gene_num=as.numeric(hog_ct_per_sp.df$gained_hog_num[match(Organism.Name, hog_ct_per_sp.df$species_name)]),
                                                    virus_gained_dupl_gene_num=as.numeric(hog_ct_per_sp.df$gained_dupl_hog_set_num[match(Organism.Name, hog_ct_per_sp.df$species_name)]),
                                                    virus_gained_hgt_gene_num=as.numeric(hog_ct_per_sp.df$gained_hgt_hog_set_num[match(Organism.Name, hog_ct_per_sp.df$species_name)]),
                                                    virus_gained_gene_prop=as.numeric(hog_ct_per_sp.df$prop_gained_hog_num[match(Organism.Name, hog_ct_per_sp.df$species_name)]),
                                                    virus_gained_dupl_gene_prop=as.numeric(hog_ct_per_sp.df$prop_gained_dupl_hog_set_num[match(Organism.Name, hog_ct_per_sp.df$species_name)]),
                                                    virus_gained_hgt_gene_prop=as.numeric(hog_ct_per_sp.df$prop_gained_hgt_hog_set_num[match(Organism.Name, hog_ct_per_sp.df$species_name)]))

poxv_virus_host.df[which(poxv_virus_host.df$Organism.Name=="Monkeypox_virus"),"virus_gene_num"] <- poxv_genome_stats.df %>% filter(str_detect(species_name, "Monkeypox_virus")) %>% pull(gene_number) %>% mean()
poxv_virus_host.df[which(poxv_virus_host.df$Organism.Name=="Monkeypox_virus"),"virus_genome_length"] <- poxv_genome_stats.df %>% filter(str_detect(species_name, "Monkeypox_virus")) %>% pull(genome_length) %>% mean()
poxv_virus_host.df[which(poxv_virus_host.df$Organism.Name=="Monkeypox_virus"),"virus_total_gene_length"] <- poxv_genome_stats.df %>% filter(str_detect(species_name, "Monkeypox_virus")) %>% pull(total_gene_length) %>% mean()
poxv_virus_host.df[which(poxv_virus_host.df$Organism.Name=="Monkeypox_virus"),"virus_total_cds_length"] <- poxv_genome_stats.df %>% filter(str_detect(species_name, "Monkeypox_virus")) %>% pull(total_cds_length) %>% mean()
poxv_virus_host.df[which(poxv_virus_host.df$Organism.Name=="Monkeypox_virus"),"virus_og_count_pooled"] <- hog_ct_per_sp.df %>% filter(str_detect(species_name, "Monkeypox_virus")) %>% pull(og_count_pooled) %>% mean()
poxv_virus_host.df[which(poxv_virus_host.df$Organism.Name=="Monkeypox_virus"),"virus_gained_gene_num"] <- hog_ct_per_sp.df %>% filter(str_detect(species_name, "Monkeypox_virus")) %>% pull(gained_hog_num) %>% mean()
poxv_virus_host.df[which(poxv_virus_host.df$Organism.Name=="Monkeypox_virus"),"virus_gained_dupl_gene_num"] <- hog_ct_per_sp.df %>% filter(str_detect(species_name, "Monkeypox_virus")) %>% pull(gained_dupl_hog_set_num) %>% mean()
poxv_virus_host.df[which(poxv_virus_host.df$Organism.Name=="Monkeypox_virus"),"virus_gained_hgt_gene_num"] <- hog_ct_per_sp.df %>% filter(str_detect(species_name, "Monkeypox_virus")) %>% pull(gained_hgt_hog_set_num) %>% mean()
poxv_virus_host.df[which(poxv_virus_host.df$Organism.Name=="Monkeypox_virus"),"virus_gained_gene_prop"] <- hog_ct_per_sp.df %>% filter(str_detect(species_name, "Monkeypox_virus")) %>% pull(prop_gained_hog_num) %>% mean()
poxv_virus_host.df[which(poxv_virus_host.df$Organism.Name=="Monkeypox_virus"),"virus_gained_dupl_gene_prop"] <- hog_ct_per_sp.df %>% filter(str_detect(species_name, "Monkeypox_virus")) %>% pull(prop_gained_dupl_hog_set_num) %>% mean()
poxv_virus_host.df[which(poxv_virus_host.df$Organism.Name=="Monkeypox_virus"),"virus_gained_hgt_gene_prop"] <- hog_ct_per_sp.df %>% filter(str_detect(species_name, "Monkeypox_virus")) %>% pull(prop_gained_hgt_hog_set_num) %>% mean()

poxv_virus_host.df <- poxv_virus_host.df %>% mutate(virus_prop_gene_length=as.numeric(virus_total_gene_length/virus_genome_length),
                                                    virus_prop_cds_length=as.numeric(virus_total_cds_length/virus_genome_length))

poxv_virus_host.df <- as.data.frame(poxv_virus_host.df)

poxv_t.rt.rn.host <- drop.tip(poxv_t.rt.rn, setdiff(poxv_t.rt.rn$tip.label[which(!str_detect(poxv_t.rt.rn$tip.label, "Monkeypox"))], poxv_virus_host.df$Organism.Name))
poxv_t.rt.rn.host <- drop.tip(poxv_t.rt.rn.host, poxv_t.rt.rn.host$tip.label[which(str_detect(poxv_t.rt.rn.host$tip.label, "Monkeypox"))], trim.internal = FALSE)
poxv_t.rt.rn.host <- drop.tip(poxv_t.rt.rn.host, "98/100")
poxv_t.rt.rn.host$tip.label[which(poxv_t.rt.rn.host$tip.label=="")] <- "Monkeypox_virus"

poxv_virus_host.df.host_ord <- as.data.frame(tibble(Organism.Name = poxv_t.rt.rn.host$tip.label) %>% left_join(poxv_virus_host.df, by = "Organism.Name"))

virus_host_cor.df <- data.frame(factor_1 = character(0), factor_2 = character(0), rho = character(0), p_value = character(0), pic_rho = character(0), pic_p_value = character(0), stringsAsFactors = FALSE)

host_param_order.v <- c("av_host_gene_num","av_host_genome_length","av_host_gene_density","av_host_prop_gene_length","av_host_prop_cds_length")
virus_param_order.v <- c("virus_gene_num","virus_og_count_pooled","virus_gained_gene_prop","virus_gained_dupl_gene_prop","virus_gained_hgt_gene_prop","virus_genome_length","virus_prop_gene_length","virus_prop_cds_length")

for (a in host_param_order.v) {
  for (b in virus_param_order.v) {
    poxv_t.rt.rn.host.test <- drop.tip(poxv_t.rt.rn.host, poxv_virus_host.df.host_ord %>% filter(!!sym(a)==0 | is.nan(!!sym(a))) %>% pull(Organism.Name))
    poxv_virus_host.df.host_ord.test <- poxv_virus_host.df.host_ord %>% filter(!!sym(a)!=0, !is.nan(!!sym(a)))
    
    virus_host_cor.res <- cor.test(poxv_virus_host.df.host_ord.test[,a], poxv_virus_host.df.host_ord.test[,b], method = "spearman", exact = FALSE)
    virus_host_cor_sep.df <- data.frame(factor_1 = character(0), factor_2 = character(0), rho = character(0), p_value = character(0), pic_rho = character(0), pic_p_value = character(0), stringsAsFactors = FALSE)
    virus_host_cor_sep.df[1,1] <- a
    virus_host_cor_sep.df[1,2] <- b
    virus_host_cor_sep.df[1,3] <- virus_host_cor.res$estimate
    virus_host_cor_sep.df[1,4] <- virus_host_cor.res$p.value
    
    test_1_pic <- pic(poxv_virus_host.df.host_ord.test[,a], poxv_t.rt.rn.host.test)
    test_2_pic <- pic(poxv_virus_host.df.host_ord.test[,b], poxv_t.rt.rn.host.test)
    
    test_1_vs_test_2_cont_picModel <- lm(test_1_pic ~ test_2_pic - 1)
    
    virus_host_cor_sep.df[1,5] <- sign(coef(test_1_vs_test_2_cont_picModel)) * sqrt(summary(test_1_vs_test_2_cont_picModel)$r.squared)
    virus_host_cor_sep.df[1,6] <- summary(test_1_vs_test_2_cont_picModel)$coefficients[1, 4]
    
    virus_host_cor.df <- bind_rows(virus_host_cor.df, virus_host_cor_sep.df)
  }
}

virus_host_cor.df <- virus_host_cor.df %>% mutate(rho=as.numeric(rho), p_value=as.numeric(p_value)) %>% arrange(pic_p_value)
virus_host_cor.df <- virus_host_cor.df %>% mutate(rho = as.numeric(rho), p_value = as.numeric(p_value), pic_rho = as.numeric(pic_rho), pic_p_value = as.numeric(pic_p_value))

virus_host_cor_rho_mat <- virus_host_cor.df %>% select(factor_1, factor_2, rho) %>% pivot_wider(names_from = factor_2, values_from = rho) %>% column_to_rownames("factor_1") %>% as.matrix()
virus_host_cor_rho_mat <- t(virus_host_cor_rho_mat[host_param_order.v,virus_param_order.v])

virus_host_cor_p_mat <- virus_host_cor.df %>% select(factor_1, factor_2, p_value) %>% pivot_wider(names_from = factor_2, values_from = p_value) %>% column_to_rownames("factor_1") %>% as.matrix()
virus_host_cor_p_mat <- t(virus_host_cor_p_mat[host_param_order.v,virus_param_order.v])

virus_host_cor_rho_mat_ht <- Heatmap(virus_host_cor_rho_mat, name = "rho", cluster_rows = FALSE, cluster_columns = FALSE, row_names_side = "left", column_names_side = "top",
                                     row_gap = unit(5, "mm"), column_gap = unit(5, "mm"), col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")), na_col = "black", heatmap_legend_param = list(title = "rho"),
                                     # cell_fun = function(j, i, x, y, w, h, fill) {if (virus_host_cor_p_mat[i, j] < 0.05) {grid.text("*", x, y, gp = gpar(fontsize = 24))}})
                                     cell_fun = function(j, i, x, y, w, h, fill) {if (virus_host_cor_p_mat[i, j] < 0.05) {grid.text(round(virus_host_cor_rho_mat[i, j],2), x, y, gp = gpar(fontsize = 24))}})

virus_host_cor_pic_rho_mat <- virus_host_cor.df %>% select(factor_1, factor_2, pic_rho) %>% pivot_wider(names_from = factor_2, values_from = pic_rho) %>% column_to_rownames("factor_1") %>% as.matrix()
virus_host_cor_pic_rho_mat <- t(virus_host_cor_pic_rho_mat[host_param_order.v,virus_param_order.v])

virus_host_cor_pic_p_mat <- virus_host_cor.df %>% select(factor_1, factor_2, pic_p_value) %>% pivot_wider(names_from = factor_2, values_from = pic_p_value) %>% column_to_rownames("factor_1") %>% as.matrix()
virus_host_cor_pic_p_mat <- t(virus_host_cor_pic_p_mat[host_param_order.v,virus_param_order.v])

virus_host_cor_pic_rho_mat_ht <- Heatmap(virus_host_cor_pic_rho_mat, name = "pic_rho", cluster_rows = FALSE, cluster_columns = FALSE, row_names_side = "left", column_names_side = "top",
                                         row_gap = unit(5, "mm"), column_gap = unit(5, "mm"), col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")), na_col = "black", heatmap_legend_param = list(title = "pic_rho"),
                                         # cell_fun = function(j, i, x, y, w, h, fill) {if (!is.nan(virus_host_cor_pic_rho_mat[i, j]) && virus_host_cor_pic_p_mat[i, j] < 0.05) {grid.text("*", x, y, gp = gpar(fontsize = 24))}})
                                         cell_fun = function(j, i, x, y, w, h, fill) {if (!is.nan(virus_host_cor_pic_rho_mat[i, j]) && virus_host_cor_pic_p_mat[i, j] < 0.05) {grid.text(round(virus_host_cor_pic_rho_mat[i, j],2), x, y, gp = gpar(fontsize = 24))}})

virus_host_cor_rho_mat_ht + virus_host_cor_pic_rho_mat_ht

########## ***** Supplementary Figure 1C: Correlation between virus and host genomic parameters and that between their contrasts *****
# pdf("Figures/spearman_corr_and_pic_corr_heatmap.pdf", width=24, height=16)
# print(virus_host_cor_rho_mat_ht + virus_host_cor_pic_rho_mat_ht)
# dev.off()


########## Check correlation between virus gene number and host protein-coding gene number
poxv_gene_no_vs_gene_no_host_cor.res <- cor.test((poxv_virus_host.df %>% filter(!clade %in% c("Salmonpoxvirus"), av_host_gene_num!=0))$av_host_gene_num, 
                                                 (poxv_virus_host.df %>% filter(!clade %in% c("Salmonpoxvirus"), av_host_gene_num!=0))$virus_gene_num, 
                                                 method = "spearman", exact = FALSE)
poxv_gene_no_vs_gene_no_host_cor.res

poxv_gene_no_vs_gene_no_host_lm.res <- lm(virus_gene_num ~ av_host_gene_num, data = poxv_virus_host.df %>% filter(!clade %in% c("Salmonpoxvirus"), av_host_gene_num!=0))

########## ***** Figure 1C: Correlation between virus gene number and host protein-coding gene number *****
poxv_gene_no_vs_gene_no_host.p <- ggplot(poxv_virus_host.df %>% filter(av_host_gene_num!=0), aes(x=av_host_gene_num, y=virus_gene_num, color=clade))
poxv_gene_no_vs_gene_no_host.p <- poxv_gene_no_vs_gene_no_host.p + geom_point()
poxv_gene_no_vs_gene_no_host.p <- poxv_gene_no_vs_gene_no_host.p + geom_text(aes(label=substr(clade, 1, 3)), hjust=-0.13, vjust=0, size=7)
poxv_gene_no_vs_gene_no_host.p <- poxv_gene_no_vs_gene_no_host.p + geom_abline(slope = coef(poxv_gene_no_vs_gene_no_host_lm.res)[2], intercept = coef(poxv_gene_no_vs_gene_no_host_lm.res)[1], color = "black")
poxv_gene_no_vs_gene_no_host.p <- poxv_gene_no_vs_gene_no_host.p + scale_color_manual("clade", name="Clade", values=poxv_color.clade.v, breaks=names(poxv_color.clade.v), na.value = NA)
poxv_gene_no_vs_gene_no_host.p <- poxv_gene_no_vs_gene_no_host.p + labs(title="Relationship between host and virus gene number", 
                                                                        subtitle=paste("r = ", round(poxv_gene_no_vs_gene_no_host_cor.res$estimate, 2), " ; P-value = ", poxv_gene_no_vs_gene_no_host_cor.res$p.value), 
                                                                        x="Host protein-coding gene number", y="Virus gene number")
poxv_gene_no_vs_gene_no_host.p <- poxv_gene_no_vs_gene_no_host.p + guides(color="none")
poxv_gene_no_vs_gene_no_host.p <- poxv_gene_no_vs_gene_no_host.p + theme_classic()
poxv_gene_no_vs_gene_no_host.p <- poxv_gene_no_vs_gene_no_host.p + theme(text=element_text(size = 24), plot.title=element_text(size = 32), axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
poxv_gene_no_vs_gene_no_host.p <- poxv_gene_no_vs_gene_no_host.p + scale_x_continuous(breaks=seq(10000, 50000, by=10000), limits=c(10000,50000))
poxv_gene_no_vs_gene_no_host.p <- poxv_gene_no_vs_gene_no_host.p + scale_y_continuous(limits=c(100,350))
poxv_gene_no_vs_gene_no_host.p

# pdf("Figures/virus_gene_number_vs_host_protein_coding_gene_number.pdf", width=16, height=12)
# print(poxv_gene_no_vs_gene_no_host.p)
# dev.off()


########## Check correlation between virus gene number and host genome size
poxv_gene_no_vs_genome_l_host_cor.res <- cor.test((poxv_virus_host.df %>% filter(!clade %in% c("Deltaentomopoxvirus")))$av_host_genome_length, 
                                                  (poxv_virus_host.df %>% filter(!clade %in% c("Deltaentomopoxvirus")))$virus_gene_num, 
                                                  method = "spearman", exact = FALSE)
poxv_gene_no_vs_genome_l_host_cor.res

poxv_gene_no_vs_genome_l_host_lm.res <- lm(virus_gene_num ~ av_host_genome_length, data = poxv_virus_host.df %>% filter(!clade %in% c("Deltaentomopoxvirus")))

########## ***** Supplementary Figure 1D: Correlation between virus gene number and host genome size *****
poxv_gene_no_vs_genome_l_host.p <- ggplot(poxv_virus_host.df, aes(x=av_host_genome_length, y=virus_gene_num, color=clade))
poxv_gene_no_vs_genome_l_host.p <- poxv_gene_no_vs_genome_l_host.p + geom_point()
poxv_gene_no_vs_genome_l_host.p <- poxv_gene_no_vs_genome_l_host.p + geom_text(aes(label=substr(clade, 1, 3)), hjust=-0.13, vjust=0, size=7)
poxv_gene_no_vs_genome_l_host.p <- poxv_gene_no_vs_genome_l_host.p + geom_abline(slope = coef(poxv_gene_no_vs_genome_l_host_lm.res)[2], intercept = coef(poxv_gene_no_vs_genome_l_host_lm.res)[1], color = "black")
poxv_gene_no_vs_genome_l_host.p <- poxv_gene_no_vs_genome_l_host.p + scale_color_manual("clade", name="Clade", values=poxv_color.clade.v, breaks=names(poxv_color.clade.v), na.value = NA)
poxv_gene_no_vs_genome_l_host.p <- poxv_gene_no_vs_genome_l_host.p + labs(title="Relationship between host genome size and virus gene number", 
                                                                          subtitle=paste("r = ", round(poxv_gene_no_vs_genome_l_host_cor.res$estimate, 2), " ; P-value = ", poxv_gene_no_vs_genome_l_host_cor.res$p.value), 
                                                                          x="Host genome size", y="Virus gene number")
poxv_gene_no_vs_genome_l_host.p <- poxv_gene_no_vs_genome_l_host.p + guides(color="none")
poxv_gene_no_vs_genome_l_host.p <- poxv_gene_no_vs_genome_l_host.p + theme_classic()
poxv_gene_no_vs_genome_l_host.p <- poxv_gene_no_vs_genome_l_host.p + theme(text=element_text(size = 24), plot.title=element_text(size = 32), axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
poxv_gene_no_vs_genome_l_host.p <- poxv_gene_no_vs_genome_l_host.p + scale_x_continuous(breaks=seq(0, 1e+10, by=2e+9), limits=c(0,1e+10))
poxv_gene_no_vs_genome_l_host.p <- poxv_gene_no_vs_genome_l_host.p + scale_y_continuous(limits=c(100,350))
poxv_gene_no_vs_genome_l_host.p

# pdf("Figures/virus_gene_number_vs_host_genome_size.pdf", width=16, height=12)
# print(poxv_gene_no_vs_genome_l_host.p)
# dev.off()



############################## Synteny analysis ##############################
########## Create a genome length data frame
s0 <- tibble(poxv_genome_stats.df %>% select(species_name, genome_length)) %>% rename(seq_id = species_name, length = genome_length)
s0 <- s0 %>% mutate(seq_id = factor(seq_id, levels = poxv_t.p$data %>% filter(isTip) %>% arrange(-y) %>% pull(label))) %>% arrange(seq_id)

########## Create a gene position data frame
for (a in 1:nrow(poxv_genome_stats.df)) {
  poxv_region.f <- paste("251208_synteny_analysis/", poxv_genome_stats.df$assembly_name[a], "_prot_gene_map.tsv", sep="")
  poxv_region <- read.delim(poxv_region.f, check.names = FALSE)
  
  if (ncol(str_split(poxv_region$coded_by, "\\..", simplify=TRUE)) > 3) {
    print(paste("ERROR: There is an intron annotation in the file", poxv_region.f, sep="")) # Check presences of any intron annotations
    break
  }
  
  poxv_region <- poxv_region %>% mutate(strand=ifelse(str_detect(coded_by, "complement"), "-", "+"), coded_by=gsub("<|>|)", "", gsub("complement(", "", coded_by, fixed=TRUE)),
                                        start=as.numeric(str_split(str_split(coded_by, ":", simplify=TRUE)[,2], "\\..", simplify=TRUE)[,1]),
                                        stop=as.numeric(str_split(str_split(coded_by, ":", simplify=TRUE)[,2], "\\..", simplify=TRUE)[,2])) %>% select(start, stop, version)
  colnames(poxv_region) <- c("start", "end", "gene_ID")

  poxv_region <- poxv_region %>% mutate(species = poxv_genome_stats.df$accession_no[a])

  attribute <- read.delim(poxv_region.f, check.names = FALSE)
  attribute <- attribute %>% mutate(strand=ifelse(str_detect(coded_by, "complement"), "-", "+"), group=hog_res_long.df$HOG[match(version, hog_res_long.df$prot_accessions)]) %>% select(version, gene_name, group, strand)
  colnames(attribute) <- c("gene_ID", "symbol", "group", "strand")
  
  if (a==1) {poxv_regions <- poxv_region} else {poxv_regions <- rbind(poxv_regions, poxv_region)}
  if (a==1) {attributes <- attribute} else {attributes <- rbind(attributes, attribute)}
}  

attributes <- attributes %>% filter(!is.na(group))

poxv_regions <- poxv_regions %>% left_join(attributes, by = "gene_ID")

tblastn_bed.f <- "251208_synteny_analysis/260207_tblastn_result_3ntplus.bed"
tblastn_bed <- read.delim(tblastn_bed.f, check.names = FALSE, header = FALSE)

colnames(tblastn_bed) <- c("accession_id","start","end","gene_ID","score","strand")

tblastn_regions <- tblastn_bed %>% select(start, end, gene_ID, "strand")
tblastn_regions <- tblastn_regions %>% mutate(start=start+1, species=paste(str_split(str_split(gene_ID, "_vs_", simplify = TRUE)[,2], "\\.", simplify = TRUE)[,1], ".1", sep=""))
tblastn_regions <- tblastn_regions %>% mutate(symbol=gene_ID, group=str_split(gene_ID, "_vs_", simplify = TRUE)[,1]) %>% select(start, end, gene_ID, species, symbol, group, strand)

poxv_regions <- rbind(poxv_regions, tblastn_regions)
poxv_regions <- poxv_regions %>% mutate(gene_mid = start + (end-start)/2)
poxv_regions <- poxv_regions %>% mutate(group = case_when(is.na(group) ~ "Other", T ~ group)) 
poxv_regions <- poxv_regions %>% mutate(clade=poxv_ncbi_ds.df$clade[match(species, poxv_ncbi_ds.df$Assembly.Accession)]) 
poxv_regions <- poxv_regions %>% mutate(Subfamily=poxv_ncbi_ds.df$Subfamily[match(species, poxv_ncbi_ds.df$Assembly.Accession)]) 
poxv_regions <- poxv_regions %>% mutate(species=poxv_ncbi_ds.df$tip_label_dummy[match(species, poxv_ncbi_ds.df$Assembly.Accession)])
poxv_regions <- poxv_regions %>% mutate(min = 1, max = s0$length[match(species, s0$seq_id)], seq_center=(max+1)/2, max_dist_center=(max-1)/2)
poxv_regions <- poxv_regions %>% mutate(rel_loc=abs(gene_mid-seq_center)/max_dist_center)
poxv_regions <- poxv_regions %>% mutate(species.f = factor(species, levels = poxv_t.p$data %>% filter(isTip) %>% arrange(y) %>% pull(label))) %>% mutate(species.n = as.numeric(species.f))

og_occupancy.df <- poxv_regions %>% filter(!is.na(group)) %>% distinct(group, species) %>% group_by(group) %>% summarize(og_sp_total=n()) %>% mutate(occupancy=og_sp_total/n_distinct(poxv_regions$species))

poxv_regions <- poxv_regions %>% mutate(og_sp_total = og_occupancy.df$og_sp_total[match(group, og_occupancy.df$group)])
poxv_regions <- poxv_regions %>% mutate(gene_status = poxv_og_gain_loss_status.df$gained_status[match(group, gsub("N0HOG", "N0.HOG", poxv_og_gain_loss_status.df$poxv_og_name, fixed=TRUE))])
poxv_regions <- poxv_regions %>% filter(group != "Other")
poxv_regions <- poxv_regions %>% arrange(species, start, end)


########## Define boundaries of central and terminal genomic regions using a sliding-window approach
window_width <- 0.10
window_step <- 0.01
minimum_genes <- 10L
minimum_significant_run <- 5L

##### Prepare gene-level data
poxv_genes.df <- poxv_regions %>% left_join(og_occupancy.df, by="group") %>% mutate(x=(gene_mid-min)/(max-min), is_preexisting=case_when(gene_status=="pre-existing" ~ TRUE, gene_status=="gained" ~ FALSE, TRUE ~ NA))

##### Calculate sliding-window statistics
make_windows <- function(d, width=window_width, step=window_step, min_genes=minimum_genes) {
  window_mid.v <- seq(width/2, 1-width/2, by=step)
  map_dfr(window_mid.v, function(window_mid) {
    in_window <- !is.na(d$x) & d$x>=window_mid-width/2 & d$x<window_mid+width/2
    classifiable <- in_window & !is.na(d$is_preexisting)
    mean_occupancy <- if(any(in_window & !is.na(d$occupancy))) {mean(d$occupancy[in_window], na.rm=TRUE)} else {NA_real_}
    tibble(species=first(d$species), clade=first(d$clade), Subfamily=first(d$Subfamily), window_mid=window_mid, n_classifiable=sum(classifiable), n_preexisting=sum(d$is_preexisting[classifiable]), prop_preexisting=ifelse(sum(classifiable)>=min_genes, mean(d$is_preexisting[classifiable]), NA_real_), mean_occupancy=mean_occupancy)
  })
}

poxv_windows.df <- poxv_genes.df %>% group_by(species) %>% group_split() %>% map_dfr(~make_windows(.x, width=window_width, step=window_step, min_genes=minimum_genes))

##### Define the middle 20% as the reference central region
central_reference.df <- poxv_genes.df %>% filter(x>=0.4, x<=0.6, !is.na(is_preexisting)) %>% group_by(species) %>% summarise(core_pre=sum(is_preexisting), core_nonpre=sum(!is_preexisting), .groups="drop") %>% mutate(core_prop=core_pre/(core_pre+core_nonpre))

##### Test whether terminal windows have fewer pre-existing genes
window_tests.df <- poxv_windows.df %>% inner_join(central_reference.df, by="species") %>% filter(!is.na(prop_preexisting), window_mid<=0.35 | window_mid>=0.65) %>% rowwise() %>%
  mutate(p_value=fisher.test(matrix(c(n_preexisting, n_classifiable-n_preexisting, core_pre, core_nonpre), nrow=2, byrow=TRUE), alternative="less")$p.value) %>% ungroup() %>% group_by(species) %>% mutate(adj_p=p.adjust(p_value, method="BH"), significantly_lower=adj_p<0.05 & prop_preexisting<core_prop) %>% ungroup()

### Add test results to all windows
poxv_windows_tested.df <- poxv_windows.df %>% left_join(window_tests.df %>% select(species, window_mid, adj_p, significantly_lower), by=c("species", "window_mid")) %>% left_join(central_reference.df %>% select(species, core_prop), by="species")

##### Plot significant terminal windows by clade
poxv_windows_tested.p <- ggplot(poxv_windows_tested.df %>% filter(!is.na(prop_preexisting)), aes(x=window_mid, y=prop_preexisting, group=species))
poxv_windows_tested.p <- poxv_windows_tested.p + geom_line(alpha=0.25)
poxv_windows_tested.p <- poxv_windows_tested.p + geom_point(data=poxv_windows_tested.df %>% filter(significantly_lower %in% TRUE, !is.na(prop_preexisting)), color="red", size=1)
poxv_windows_tested.p <- poxv_windows_tested.p + facet_wrap(~clade)
poxv_windows_tested.p <- poxv_windows_tested.p + scale_x_continuous(breaks=seq(0, 1, 0.2))
poxv_windows_tested.p <- poxv_windows_tested.p + scale_y_continuous(limits=c(0, 1))
poxv_windows_tested.p <- poxv_windows_tested.p + labs(x="Normalized genomic position", y="Proportion of pre-existing genes")
poxv_windows_tested.p <- poxv_windows_tested.p + theme_classic()
poxv_windows_tested.p

##### Identify runs of consecutive significant windows
mark_long_runs <- function(x, min_run=3L) {
  x[is.na(x)] <- FALSE
  runs <- rle(x)
  rep(runs$values & runs$lengths>=min_run, runs$lengths)
}

##### Estimate left and right central-region boundaries
estimate_boundaries <- function(d, min_run=minimum_significant_run, step=0.01) {
  d <- d %>% arrange(window_mid)
  left.df <- d %>% filter(window_mid<0.5)
  right.df <- d %>% filter(window_mid>0.5)
  left_keep <- mark_long_runs(left.df$significantly_lower, min_run=min_run)
  right_keep <- mark_long_runs(right.df$significantly_lower, min_run=min_run)
  left_boundary <- if(any(left_keep)) {min(0.5, max(left.df$window_mid[left_keep])+step/2)} else {NA_real_}
  right_boundary <- if(any(right_keep)) {max(0.5, min(right.df$window_mid[right_keep])-step/2)} else {NA_real_}
  tibble(species=first(d$species), clade=first(d$clade), Subfamily=first(d$Subfamily), left_boundary=left_boundary, right_boundary=right_boundary)
}

##### Estimate boundaries for every genome
boundary.df <- poxv_windows_tested.df %>% group_by(species) %>% group_split() %>% map_dfr(~estimate_boundaries(.x, min_run=minimum_significant_run, step=window_step))

##### Convert normalized boundaries to nucleotide coordinates
genome_limits.df <- poxv_genes.df %>% group_by(species) %>% summarise(genome_min=first(min), genome_max=first(max), .groups="drop")

boundary.df <- boundary.df %>% left_join(genome_limits.df, by="species") %>% mutate(left_boundary_bp=(genome_min+left_boundary*(genome_max-genome_min)), right_boundary_bp=(genome_min+right_boundary*(genome_max-genome_min)), species_label=str_wrap(str_replace_all(species, "_", " "), width=35))
boundary.df %>% select(species, left_boundary, right_boundary, left_boundary_bp, right_boundary_bp) %>% arrange(species)

##### Prepare data for plotting individual genomes
plot.df <- poxv_windows_tested.df %>% mutate(species_label=str_wrap(str_replace_all(species, "_", " "), width=35))
core_plot.df <- plot.df %>% filter(!is.na(core_prop)) %>% distinct(species_label, core_prop)
central_rect.df <- boundary.df %>% filter(!is.na(left_boundary), !is.na(right_boundary))
boundary_long.df <- boundary.df %>% select(species, species_label, clade, Subfamily, left_boundary, right_boundary) %>% pivot_longer(cols=c(left_boundary, right_boundary), names_to="boundary_side", values_to="boundary") %>% filter(!is.na(boundary))

##### Plot every genome in separate panels
genome_ct_regions.p <- ggplot(plot.df %>% filter(!is.na(prop_preexisting)), aes(x=window_mid, y=prop_preexisting))
genome_ct_regions.p <- genome_ct_regions.p + geom_rect(data=central_rect.df, aes(xmin=left_boundary, xmax=right_boundary, ymin=-Inf, ymax=Inf), inherit.aes=FALSE, fill="steelblue", alpha=0.12)
genome_ct_regions.p <- genome_ct_regions.p + geom_line(linewidth=0.7)
genome_ct_regions.p <- genome_ct_regions.p + geom_hline(data=core_plot.df, aes(yintercept=core_prop), inherit.aes=FALSE, linetype="dashed", color="blue", linewidth=0.5)
genome_ct_regions.p <- genome_ct_regions.p + geom_point(data=plot.df %>% filter(significantly_lower %in% TRUE, !is.na(prop_preexisting)), color="red", size=1.5)
genome_ct_regions.p <- genome_ct_regions.p + geom_vline(data=boundary_long.df, aes(xintercept=boundary), inherit.aes=FALSE, linetype="dashed", color="black", linewidth=0.6)
genome_ct_regions.p <- genome_ct_regions.p + facet_wrap(~species_label, ncol=4)
genome_ct_regions.p <- genome_ct_regions.p + scale_x_continuous(breaks=seq(0, 1, 0.2))
genome_ct_regions.p <- genome_ct_regions.p + scale_y_continuous(limits=c(0, 1))
genome_ct_regions.p <- genome_ct_regions.p + labs(x="Normalized genomic position", y="Proportion of pre-existing genes", caption="Red points: significantly depleted windows; black lines: inferred boundaries; blue shading: inferred central region")
genome_ct_regions.p <- genome_ct_regions.p + theme_classic()
genome_ct_regions.p <- genome_ct_regions.p + theme(strip.text=element_text(size=7))
genome_ct_regions.p


########## Plot syntenic relationship across orthologous genes in ortholog groups
g0 <- tibble(poxv_regions %>% select(species, start, end, strand, group, gene_ID)) %>% rename(seq_id = species, start = start, end = end, feat_id = group)
g0 <- g0 %>% mutate(gene_status = poxv_og_gain_loss_status.df$gained_status[match(feat_id, gsub("N0HOG", "N0.HOG", poxv_og_gain_loss_status.df$poxv_og_name, fixed=TRUE))])
g0 <- g0 %>% mutate(dupl_status=ifelse(gene_ID %in% gene_dup_prot_list.v, "dupl", "non-dupl"))
g0 <- g0 %>% mutate(hgt_status=ifelse(gene_ID %in% hgt_gained_res.np.df$qseqid, "hgt", "non-hgt"))

g0 <- g0 %>% mutate(border=ifelse(feat_id=="N0.HOG0000000", "N0.HOG0000000", 
                           ifelse(feat_id=="N0.HOG0000007", "N0.HOG0000007", "no")))

attributes <- rbind(attributes, tblastn_regions %>% select(gene_ID, symbol, group, strand))

l0 <- attributes %>% group_by(group) %>% filter(n() >= 2) %>% summarise(pairs = list(combn(gene_ID, 2, simplify = FALSE)), .groups = "drop") %>%
  unnest(pairs) %>% mutate(seq_id = map_chr(pairs, 1), seq_id2 = map_chr(pairs, 2)) %>% select(group, seq_id, seq_id2)
l0 <- l0 %>% mutate(og_sp_total = og_occupancy.df$og_sp_total[match(group, og_occupancy.df$group)])
l0 <- l0 %>% mutate(strand = poxv_regions$strand[match(seq_id, poxv_regions$gene_ID)],
                    start = ifelse(strand == "+", poxv_regions$start[match(seq_id, poxv_regions$gene_ID)], poxv_regions$end[match(seq_id, poxv_regions$gene_ID)]), 
                    end = ifelse(strand == "+", poxv_regions$end[match(seq_id, poxv_regions$gene_ID)], poxv_regions$start[match(seq_id, poxv_regions$gene_ID)]),
                    strand2 = poxv_regions$strand[match(seq_id2, poxv_regions$gene_ID)],
                    start2 = ifelse(strand == "+", poxv_regions$start[match(seq_id2, poxv_regions$gene_ID)], poxv_regions$end[match(seq_id2, poxv_regions$gene_ID)]),
                    end2 = ifelse(strand == "+", poxv_regions$end[match(seq_id2, poxv_regions$gene_ID)], poxv_regions$start[match(seq_id2, poxv_regions$gene_ID)]),
                    seq_id = poxv_regions$species[match(seq_id, poxv_regions$gene_ID)],
                    seq_id2 = poxv_regions$species[match(seq_id2, poxv_regions$gene_ID)])
l0 <- l0 %>% select(-group, -strand, -strand2, seq_id, start, end, seq_id2, start2, end2, og_sp_total)
l0 <- tibble(l0)

########## ***** Figure 3A: Gene synteny across poxvirus genomes *****
poxv_synteny_gggenomes.p <- gggenomes(genes=g0, seqs=s0, links=l0)
poxv_synteny_gggenomes.p <- poxv_synteny_gggenomes.p + geom_seq() 
poxv_synteny_gggenomes.p <- poxv_synteny_gggenomes.p + geom_link(aes(fill = og_sp_total, color = og_sp_total), alpha = 1, offset = 0)
poxv_synteny_gggenomes.p <- poxv_synteny_gggenomes.p + scale_fill_gradient("og_sp_total", name="Number of shared species", low="white", high="blue", na.value="grey", breaks=c(0,min(poxv_regions$og_sp_total),10,20,30,40,50,max(l0$og_sp_total)), limits = c(0, max(l0$og_sp_total)))
poxv_synteny_gggenomes.p <- poxv_synteny_gggenomes.p + scale_color_gradient("og_sp_total", name="Number of shared species", low=NA, high=NA, na.value=NA)
poxv_synteny_gggenomes.p <- poxv_synteny_gggenomes.p + new_scale_fill()
poxv_synteny_gggenomes.p <- poxv_synteny_gggenomes.p + new_scale_color()
poxv_synteny_gggenomes.p <- poxv_synteny_gggenomes.p + geom_gene(aes(fill = gene_status, color = border), size = 8)
poxv_synteny_gggenomes.p <- poxv_synteny_gggenomes.p + scale_fill_manual("gene_status", values = c("gained"="orange","pre-existing"="navyblue","undeterminable"="lightgrey"))
poxv_synteny_gggenomes.p <- poxv_synteny_gggenomes.p + scale_color_manual("border", values = c("N0.HOG0000000"="hotpink", "N0.HOG0000007"="brown", "no"=NA))
poxv_synteny_gggenomes.p <- poxv_synteny_gggenomes.p |> align(.justify = "center")

# pdf("Figures/poxvirus_synteny.pdf", width=20, height=30)
# print(poxv_synteny_gggenomes.p)
# dev.off()


########## Add information for central/terminal genomic regions
boundary_region.df <- bind_rows(
  boundary.df %>% transmute(seq_id=species, region="left_t", start=genome_min, end=left_boundary_bp),
  boundary.df %>% transmute(seq_id=species, region="central", start=left_boundary_bp, end=right_boundary_bp),
  boundary.df %>% transmute(seq_id=species, region="right_t", start=right_boundary_bp, end=genome_max)
) %>% tidyr::drop_na(start, end)

##### Generate every pair of species
species.v <- unique(boundary_region.df$seq_id)
pair.m <- t(combn(species.v, 2))

b0 <- as_tibble(pair.m, .name_repair="minimal")
colnames(b0) <- c("seq_id", "seq_id2")

##### Assign matching region coordinates to both species
b0 <- tidyr::crossing(b0, region=c("left_t", "central", "right_t")) %>% left_join(boundary_region.df, by=c("seq_id", "region"))
b0 <- b0 %>% left_join(boundary_region.df %>% transmute(seq_id2=seq_id, region, start2=start, end2=end), by=c("seq_id2", "region"))
b0 <- b0 %>% select(seq_id, seq_id2, region, start, end, start2, end2)
b0 <- b0 %>% filter(!is.na(start) & !is.na(end) & !is.na(start2) & !is.na(end2))
b0 <- tibble(b0)


########## ***** Supplementary Figure 3A: Central/Terminal genomic regions *****
poxv_synteny_boundary.p <- gggenomes(genes=g0, seqs=s0, feats=boundary_region.df)
poxv_synteny_boundary.p <- poxv_synteny_boundary.p + geom_rect(data=feats(), aes(xmin=pmin(x, xend), xmax=pmax(x, xend), ymin=y-0.35, ymax=y+0.35, fill=region), color=NA, alpha=1, inherit.aes=FALSE)
poxv_synteny_boundary.p <- poxv_synteny_boundary.p + scale_fill_manual("region", values=c("left_t"="hotpink", "central"="#0080FF", "right_t"="hotpink"), breaks=c("left_t", "central", "right_t"), labels=c("Left terminal", "Central", "Right terminal"))
poxv_synteny_boundary.p <- poxv_synteny_boundary.p + new_scale_fill()
poxv_synteny_boundary.p <- poxv_synteny_boundary.p + geom_seq()
poxv_synteny_boundary.p <- poxv_synteny_boundary.p + geom_gene(aes(fill = gene_status, color = gene_status), size = 8)
poxv_synteny_boundary.p <- poxv_synteny_boundary.p + scale_fill_manual("gene_status", values = c("gained"="orange","pre-existing"="navyblue","undeterminable"="lightgrey"))
poxv_synteny_boundary.p <- poxv_synteny_boundary.p + scale_color_manual("gene_status", values = c("gained"=NA,"pre-existing"=NA,"undeterminable"=NA))
poxv_synteny_boundary.p <- poxv_synteny_boundary.p |> align(.justify = "center")

# pdf("Figures/poxvirus_synteny_boundary.pdf", width=20, height=30)
# print(poxv_synteny_boundary.p)
# dev.off()


########## ***** Figure 5C: Genomic location of Bcl-2-like gene family across poxvirus genomes *****
g0_pf06227 <- g0 %>% mutate(gene_status=ifelse(feat_id %in% pf06227_og_name.v, "PF06227", "other"))
g0_pf06227 <- g0_pf06227 %>% mutate(border=ifelse(feat_id %in% pf06227_og_name.v, "PF06227", "other"))

poxv_synteny_pf06227.p <- gggenomes(genes=g0_pf06227, seqs=s0, feats=boundary_region.df)
poxv_synteny_pf06227.p <- poxv_synteny_pf06227.p + geom_rect(data=feats(), aes(xmin=pmin(x, xend), xmax=pmax(x, xend), ymin=y-0.35, ymax=y+0.35, fill=region), color=NA, alpha=1, inherit.aes=FALSE)
poxv_synteny_pf06227.p <- poxv_synteny_pf06227.p + scale_fill_manual("region", values=c("left_t"="hotpink", "central"="#0080FF", "right_t"="hotpink"), breaks=c("left_t", "central", "right_t"), labels=c("Left terminal", "Central", "Right terminal"))
poxv_synteny_pf06227.p <- poxv_synteny_pf06227.p + new_scale_fill()
poxv_synteny_pf06227.p <- poxv_synteny_pf06227.p + geom_seq() 
poxv_synteny_pf06227.p <- poxv_synteny_pf06227.p + geom_gene(aes(fill = gene_status, color = border), size = 8)
poxv_synteny_pf06227.p <- poxv_synteny_pf06227.p + scale_fill_manual("gene_status", values = c("PF06227"="black","other"=NA))
poxv_synteny_pf06227.p <- poxv_synteny_pf06227.p + scale_color_manual("border", values = c("PF06227"="black", "other"=NA))
poxv_synteny_pf06227.p <- poxv_synteny_pf06227.p |> align(.justify = "center")

# pdf("Figures/poxvirus_synteny_pf06227.pdf", width=20, height=30)
# print(poxv_synteny_pf06227.p)
# dev.off()


########## Analyze ITRs
mummer_alm_dir <- "260823_mummer_alignment_res"
mummer_alm_file_list <- list.files(mummer_alm_dir, pattern="_nucmer_res\\.coords$", full.names=TRUE)

mummer_coord_col_names.v <- c("ref_start", "ref_end", "qry_start", "qry_end", "ref_aln_length", "qry_aln_length", "identity", "ref_length", "qry_length", "ref_coverage", "qry_coverage", "ref_id", "qry_id")

mummer_alm.df <- map_dfr(mummer_alm_file_list, function(mummer_alm_file) {
    read_tsv(mummer_alm_file, col_names=mummer_coord_col_names.v, show_col_types=FALSE, progress=FALSE) %>% mutate(source_file=basename(mummer_alm_file))
})

mummer_alm.df <- mummer_alm.df %>% mutate(genome=ref_id, orientation=if_else(qry_start > qry_end, "inverted", "direct"), alignment_length=pmin(ref_aln_length, qry_aln_length), 
                                          left_terminal_gap=ref_start-1, right_terminal_gap=qry_length-qry_start, total_terminal_gap=left_terminal_gap+right_terminal_gap,
                                          assembly_name=paste(str_split(source_file, ".1_", simplify=TRUE)[,1], '.1', sep=""))

##### Obtain genome information
mummer_genome_info.df <- mummer_alm.df %>% group_by(genome) %>% summarise(genome_length=max(ref_length), .groups="drop")
mummer_genome_info.df <- mummer_genome_info.df %>% mutate(assembly_name=mummer_alm.df$assembly_name[match(genome, mummer_alm.df$genome)],
                                                          species_name=poxv_ncbi_ds.df$tip_label_dummy[match(assembly_name, poxv_ncbi_ds.df$Assembly.Accession)],
                                                          clade=poxv_ncbi_ds.df$clade[match(assembly_name, poxv_ncbi_ds.df$Assembly.Accession)])
mummer_genome_info.df$species_name <- factor(mummer_genome_info.df$species_name, levels=(poxv_t.p$data %>% filter(isTip) %>% arrange(desc(y)) %>% pull(label)))

##### Parameters for identifying ITRs
terminal_tolerance <- 100
minimum_itr_length <- 1000
minimum_identity <- 90

##### Identify left-terminal-to-right-terminal inverted alignments
itr_candidate.df <- mummer_alm.df %>% filter(orientation=="inverted", left_terminal_gap <= terminal_tolerance, right_terminal_gap <= terminal_tolerance, alignment_length >= minimum_itr_length, identity >= minimum_identity)

##### Select the longest candidate ITR alignment from each genome
itr_all.df <- itr_candidate.df %>% arrange(genome, desc(alignment_length), desc(identity), total_terminal_gap) %>% group_by(genome) %>% slice_head(n=1) %>% ungroup()
itr_all.df <- itr_all.df %>% transmute(genome, genome_length=ref_length, left_itr_start=ref_start, left_itr_end=ref_end, right_itr_start=pmin(qry_start, qry_end),
                                       right_itr_end=pmax(qry_start, qry_end), itr_length=alignment_length, identity, left_terminal_gap, right_terminal_gap, source_file)
itr_all.df <- itr_all.df %>% mutate(assembly_name=mummer_alm.df$assembly_name[match(genome, mummer_alm.df$genome)],
                                    species_name=poxv_ncbi_ds.df$tip_label_dummy[match(assembly_name, poxv_ncbi_ds.df$Assembly.Accession)],
                                    clade=poxv_ncbi_ds.df$clade[match(assembly_name, poxv_ncbi_ds.df$Assembly.Accession)])

if (nrow(itr_all.df)==0) {stop("No candidate ITRs passed the specified thresholds.")}

##### Identify genomes without a detected ITR
no_itr.df <- mummer_genome_info.df %>% anti_join(itr_all.df, by="genome")

if (nrow(no_itr.df)>0) {
    message(nrow(no_itr.df), " genome(s) had no ITR passing the specified thresholds.")
    print(no_itr.df)
}

genome_order.v <- mummer_genome_info.df %>% arrange(species_name) %>% pull(genome)
genome_order.v <- genome_order.v[genome_order.v %in% mummer_genome_info.df$genome]
genome_order.v <- c(genome_order.v, setdiff(mummer_genome_info.df$genome, genome_order.v))

##### Order ITRs
itr_all.df <- itr_all.df %>% mutate(genome=factor(genome, levels=genome_order.v)) %>% arrange(genome) %>% mutate(genome=as.character(genome))

itr_diagnostic.df <- mummer_alm.df %>% group_by(genome) %>% summarise(n_inverted=sum(orientation=="inverted"), n_terminal=sum(orientation=="inverted" & left_terminal_gap <= terminal_tolerance & right_terminal_gap <= terminal_tolerance), 
                                                                      n_terminal_long=sum(orientation=="inverted" & left_terminal_gap <= terminal_tolerance & right_terminal_gap <= terminal_tolerance & alignment_length >= minimum_itr_length), 
                                                                      n_passing=sum(orientation=="inverted" & left_terminal_gap <= terminal_tolerance & right_terminal_gap <= terminal_tolerance & alignment_length >= minimum_itr_length & identity >= minimum_identity), .groups="drop") %>%
  mutate(status=case_when(n_inverted==0 ~ "No inverted alignment found", n_terminal==0 ~ "Inverted alignment does not reach both termini", n_terminal_long==0 ~ "Terminal inverted alignment is shorter than minimum_itr_length", n_passing==0 ~ "Terminal inverted alignment is below minimum_identity", TRUE ~ "ITR detected"))

table(itr_diagnostic.df$status)
itr_diagnostic.df %>% filter(status!="ITR detected")


itr_links.df <- itr_all.df %>% transmute(seq_id=species_name, seq_id2=species_name, start=left_itr_start, end=left_itr_end, start2=right_itr_start, end2=right_itr_end, strand="-", identity, itr_length) %>% filter(seq_id %in% s0$seq_id)
itr_region.df <- bind_rows(itr_links.df %>% transmute(seq_id, itr_side="left", itr_start=pmin(start, end), itr_end=pmax(start, end)), itr_links.df %>% transmute(seq_id, itr_side="right", itr_start=pmin(start2, end2), itr_end=pmax(start2, end2)))
g0_itr <- g0 %>% mutate(gene_mid=(start+end)/2) %>% semi_join(itr_region.df, by=join_by(seq_id, gene_mid>=itr_start, gene_mid<=itr_end))

terminal_region.df <- boundary_region.df %>% filter(region %in% c("left_t", "right_t")) %>% transmute(seq_id, region, terminal_start=pmin(start, end), terminal_end=pmax(start, end))
g0_terminal <- g0 %>% mutate(gene_mid=(start+end)/2) %>% semi_join(terminal_region.df, by=join_by(seq_id, gene_mid>=terminal_start, gene_mid<=terminal_end))

make_itr_ribbon <- function(d, depth=1, n=100) {
  t <- seq(0, 1, length.out=n)
  
  control_y1 <- d$y+4*depth/3
  control_y2 <- d$yend+4*depth/3
  
  curved_y <- (1-t)^3*d$y+
    3*(1-t)^2*t*control_y1+
    3*(1-t)*t^2*control_y2+
    t^3*d$yend
  
  outer_x <- (1-t)^3*d$x+
    3*(1-t)^2*t*d$x+
    3*(1-t)*t^2*d$xmin+
    t^3*d$xmin
  
  inner_x <- (1-t)^3*d$xend+
    3*(1-t)^2*t*d$xend+
    3*(1-t)*t^2*d$xmax+
    t^3*d$xmax
  
  tibble(ribbon_id=d$ribbon_id,
         ribbon_x=c(outer_x, rev(inner_x)),
         ribbon_y=c(curved_y, rev(curved_y)))
}

poxv_synteny_itr.p <- gggenomes(genes=g0_terminal, seqs=s0, links=itr_links.df, feats=boundary_region.df, adjacent_only=FALSE)
poxv_synteny_itr.p <- poxv_synteny_itr.p |> align(.track_id="seqs", .justify="center")

itr_link_layout.df <- pull_links(poxv_synteny_itr.p) %>% mutate(ribbon_id=row_number())

itr_ribbon.df <- map_dfr(seq_len(nrow(itr_link_layout.df)), function(i) make_itr_ribbon(itr_link_layout.df[i, ], depth=1))



########## Plot paralogous relationships among all duplicated genes
gene_dup_events_seq.df <- gene_dup_events.df
gene_dup_events_seq.df <- gene_dup_events_seq.df %>% mutate(hog_set=ifelse(HOG_genes_1!=HOG_genes_2,
                                                                    ifelse(HOG_genes_1 %in% gene_dup_events_rs.new.hog.df$HOG, gene_dup_events_rs.new.hog.df$hog_set[match(HOG_genes_1, gene_dup_events_rs.new.hog.df$HOG)],
                                                                    ifelse(HOG_genes_2 %in% gene_dup_events_rs.new.hog.df$HOG, gene_dup_events_rs.new.hog.df$hog_set[match(HOG_genes_2, gene_dup_events_rs.new.hog.df$HOG)], NA)),
                                                                    ifelse(HOG_genes_1 %in% gene_dup_events_rs.new.hog.df$HOG, gene_dup_events_rs.new.hog.df$hog_set[match(HOG_genes_1, gene_dup_events_rs.new.hog.df$HOG)], HOG_genes_1)))

clean_gene_dup_prot <- function(gene_dup_prot_list.v) {
  gene_dup_prot_list.v <- unlist(str_split(gene_dup_prot_list.v, ", "))
  gene_dup_prot_list.v <- unname(sapply(gene_dup_prot_list.v, function(x) {
    parts <- str_split(x, "_", simplify = TRUE)
    num_parts <- length(parts)
    if (parts[num_parts - 1] %in% c("NP", "XP", "YP")) {
      paste(parts[num_parts - 1], "_", parts[num_parts], sep="")
    } else if (TRUE %in% str_detect(x, "N0.HOG")) {
      paste(parts[(which(str_detect(parts, "N0.HOG")==TRUE)):num_parts], collapse = "_")
    } else {
      parts[num_parts]
    }
  }))
  return(gene_dup_prot_list.v)
}

##### Define left and right ITR intervals
itr_region_class.df <- bind_rows(
  itr_links.df %>% transmute(seq_id, region="left_itr", itr_start=pmin(start, end), itr_end=pmax(start, end)),
  itr_links.df %>% transmute(seq_id=seq_id2, region="right_itr", itr_start=pmin(start2, end2), itr_end=pmax(start2, end2))
)

### Classify genes by their midpoint
gene_region.df <- g0 %>% transmute(seq_id, gene_ID, gene_mid=(start+end)/2, gene_key=paste(seq_id, gene_ID, sep="___"))

left_itr_gene_key.v <- gene_region.df %>% inner_join(itr_region_class.df %>% filter(region=="left_itr"), by=join_by(seq_id, gene_mid>=itr_start, gene_mid<=itr_end)) %>% pull(gene_key) %>% unique()
right_itr_gene_key.v <- gene_region.df %>% inner_join(itr_region_class.df %>% filter(region=="right_itr"), by=join_by(seq_id, gene_mid>=itr_start, gene_mid<=itr_end)) %>% pull(gene_key) %>% unique()

left_terminal_gene_key.v <- gene_region.df %>% inner_join(boundary_region.df %>% filter(region=="left_t"), by=join_by(seq_id, gene_mid>=start, gene_mid<=end)) %>% pull(gene_key) %>% unique()
right_terminal_gene_key.v <- gene_region.df %>% inner_join(boundary_region.df %>% filter(region=="right_t"), by=join_by(seq_id, gene_mid>=start, gene_mid<=end)) %>% pull(gene_key) %>% unique()

central_gene_key.v <- gene_region.df %>% inner_join(boundary_region.df %>% filter(region=="central"), by=join_by(seq_id, gene_mid>=start, gene_mid<=end)) %>% pull(gene_key) %>% unique()

gene_region.df <- gene_region.df %>% mutate(gene_region=case_when(
  gene_key %in% left_itr_gene_key.v ~ "Left ITR",
  gene_key %in% right_itr_gene_key.v ~ "Right ITR",
  gene_key %in% left_terminal_gene_key.v ~ "Left non-ITR terminal",
  gene_key %in% right_terminal_gene_key.v ~ "Right non-ITR terminal",
  gene_key %in% central_gene_key.v ~ "Central",
  TRUE ~ "Unclassified"
))


l0_dupl <- gene_dup_events_seq.df %>% transmute(hog_set, genes1=str_split(`Genes 1`, ",\\s*"), genes2=str_split(`Genes 2`, ",\\s*"))
l0_dupl <- l0_dupl %>% unnest_longer(genes1, values_to="gene_ID1") %>% unnest_longer(genes2, values_to="gene_ID2")
l0_dupl <- l0_dupl %>% mutate(gene_ID1=clean_gene_dup_prot(gene_ID1), gene_ID2=clean_gene_dup_prot(gene_ID2))
l0_dupl <- l0_dupl %>% mutate(gene_ID1.temp=pmin(gene_ID1, gene_ID2), gene_ID2.temp=pmax(gene_ID1, gene_ID2))
l0_dupl <- l0_dupl %>% select(-gene_ID1, -gene_ID2) %>% rename(gene_ID1=gene_ID1.temp, gene_ID2=gene_ID2.temp)
l0_dupl <- l0_dupl %>% distinct(gene_ID1, gene_ID2, .keep_all=TRUE)
l0_dupl <- l0_dupl %>% mutate(seq_id=g0$seq_id[match(gene_ID1, g0$gene_ID)], seq_id2=g0$seq_id[match(gene_ID2, g0$gene_ID)])
l0_dupl <- l0_dupl %>% filter(seq_id==seq_id2)
l0_dupl <- l0_dupl %>% mutate(group1=attributes$group[match(gene_ID1, attributes$gene_ID)], group2=attributes$group[match(gene_ID2, attributes$gene_ID)])
l0_dupl <- l0_dupl %>% mutate(dupl1=g0$dupl_status[match(gene_ID1, g0$gene_ID)], dupl2=g0$dupl_status[match(gene_ID2, g0$gene_ID)])
l0_dupl <- l0_dupl %>% mutate(strand = poxv_regions$strand[match(gene_ID1, poxv_regions$gene_ID)],
                              start = ifelse(strand == "+", poxv_regions$start[match(gene_ID1, poxv_regions$gene_ID)], poxv_regions$end[match(gene_ID1, poxv_regions$gene_ID)]),
                              end = ifelse(strand == "+", poxv_regions$end[match(gene_ID1, poxv_regions$gene_ID)], poxv_regions$start[match(gene_ID1, poxv_regions$gene_ID)]),
                              strand2 = poxv_regions$strand[match(gene_ID2, poxv_regions$gene_ID)],
                              start2 = ifelse(strand2 == "+", poxv_regions$start[match(gene_ID2, poxv_regions$gene_ID)], poxv_regions$end[match(gene_ID2, poxv_regions$gene_ID)]),
                              end2 = ifelse(strand2 == "+", poxv_regions$end[match(gene_ID2, poxv_regions$gene_ID)], poxv_regions$start[match(gene_ID2, poxv_regions$gene_ID)]))
l0_dupl <- l0_dupl %>% filter(dupl1=="dupl", dupl2=="dupl")
l0_dupl <- l0_dupl %>% mutate(length = s0$length[match(seq_id, s0$seq_id)])

l0_dupl_ct.temp <- l0_dupl

linked_gene_key.v <- unique(c(paste(l0_dupl_ct.temp$seq_id, l0_dupl_ct.temp$gene_ID1, sep="___"), paste(l0_dupl_ct.temp$seq_id2, l0_dupl_ct.temp$gene_ID2, sep="___")))

g0_dupl <- g0 %>% mutate(gene_key=paste(seq_id, gene_ID, sep="___"))
g0_dupl <- g0_dupl %>% mutate(has_link=gene_key %in% linked_gene_key.v)
g0_dupl <- g0_dupl %>% mutate(gene_status = as.character(gene_status), 
                              gene_status=ifelse(dupl_status=="non-dupl" | has_link==FALSE, "not_show", gene_status), 
                              border=ifelse(dupl_status=="non-dupl" | has_link==FALSE, "no", 
                                     ifelse(feat_id=="N0.HOG0000000", "N0.HOG0000000", 
                                     ifelse(feat_id=="N0.HOG0000007", "N0.HOG0000007", "yes"))))

l0_dupl <- l0_dupl %>% mutate(gene_key1=paste(seq_id, gene_ID1, sep="___"), gene_key2=paste(seq_id2, gene_ID2, sep="___"))
l0_dupl <- l0_dupl %>% mutate(gene_region1=gene_region.df$gene_region[match(gene_key1, gene_region.df$gene_key)], gene_region2=gene_region.df$gene_region[match(gene_key2, gene_region.df$gene_key)])

gene_region_levels.v <- c("Left ITR", "Left non-ITR terminal", "Central", "Right non-ITR terminal", "Right ITR", "Unclassified")
l0_dupl <- l0_dupl %>% mutate(region_order1=match(gene_region1, gene_region_levels.v), region_order2=match(gene_region2, gene_region_levels.v))
l0_dupl <- l0_dupl %>% mutate(link_region_status=ifelse(region_order1<=region_order2, paste(gene_region1, gene_region2, sep="-"), paste(gene_region2, gene_region1, sep="-")))

l0_dupl <- l0_dupl %>% mutate(link_region_status=ifelse(link_region_status %in% c("Left ITR-Left ITR","Right ITR-Right ITR","Left non-ITR terminal-Left non-ITR terminal","Right non-ITR terminal-Right non-ITR terminal","Left ITR-Left non-ITR terminal","Right non-ITR terminal-Right ITR"), "Within individual terminal",
                                                 ifelse(link_region_status %in% c("Left ITR-Central","Central-Right ITR","Left non-ITR terminal-Central", "Central-Right non-ITR terminal"), "Terminal-Central", 
                                                 ifelse(link_region_status %in% c("Left ITR-Right ITR","Left non-ITR terminal-Right non-ITR terminal","Left ITR-Right non-ITR terminal","Left non-ITR terminal-Right ITR"), "Across opposite terminal",
                                                 ifelse(link_region_status=="Central-Central", "Within central", link_region_status)))))

l0_dupl <- l0_dupl %>% select(-hog_set, -dupl1, -dupl2, -strand, -strand2, -gene_ID1, -gene_ID2, -group1, -group2, -gene_key1, -gene_key2, -region_order1, -region_order2, seq_id, start, end, seq_id2, start2, end2, length, gene_region1, gene_region2, link_region_status)
l0_dupl <- tibble(l0_dupl)

########## ***** Supplemetary Figure 3B and Figure 3D: Paralogous relationship among all duplicated genes *****
poxv_synteny_dupl_gggenomes.p <- gggenomes(genes=g0_dupl, seqs=s0, links=l0_dupl, feats=boundary_region.df, adjacent_only = FALSE)
poxv_synteny_dupl_gggenomes.p <- poxv_synteny_dupl_gggenomes.p + geom_rect(data=feats(), aes(xmin=pmin(x, xend), xmax=pmax(x, xend), ymin=y-0.35, ymax=y+0.35, fill=region), color=NA, alpha=1, inherit.aes=FALSE)
poxv_synteny_dupl_gggenomes.p <- poxv_synteny_dupl_gggenomes.p + scale_fill_manual("region", values=c("left_t"="hotpink", "central"="#0080FF", "right_t"="hotpink"), breaks=c("left_t", "central", "right_t"), labels=c("Left terminal", "Central", "Right terminal"))
poxv_synteny_dupl_gggenomes.p <- poxv_synteny_dupl_gggenomes.p + new_scale_fill()
poxv_synteny_dupl_gggenomes.p <- poxv_synteny_dupl_gggenomes.p + geom_seq()
poxv_synteny_dupl_gggenomes.p <- poxv_synteny_dupl_gggenomes.p + geom_curve(data = links(), aes(x = pmin((x+xend)/2, (xmin+xmax)/2), xend = pmax((x+xend)/2, (xmin+xmax)/2), y = y, yend = yend, color = link_region_status), curvature = -0.05, alpha = 0.05)
poxv_synteny_dupl_gggenomes.p <- poxv_synteny_dupl_gggenomes.p + scale_color_manual("link_region_status", values = c("Within individual terminal"="red","Terminal-Central"="green","Across opposite terminal"="blue","Within central"="violet"), na.value = "grey")
poxv_synteny_dupl_gggenomes.p <- poxv_synteny_dupl_gggenomes.p + new_scale_color()
poxv_synteny_dupl_gggenomes.p <- poxv_synteny_dupl_gggenomes.p + geom_polygon(data=itr_ribbon.df, aes(x=ribbon_x, y=ribbon_y, group=ribbon_id), fill="purple", alpha=1, linewidth=1, inherit.aes=FALSE)
poxv_synteny_dupl_gggenomes.p <- poxv_synteny_dupl_gggenomes.p + new_scale_fill()
poxv_synteny_dupl_gggenomes.p <- poxv_synteny_dupl_gggenomes.p + geom_gene(aes(fill = gene_status, color = border), size = 8)
poxv_synteny_dupl_gggenomes.p <- poxv_synteny_dupl_gggenomes.p + scale_fill_manual("gene_status", values = c("gained"="orange","pre-existing"="navyblue","undeterminable"="lightgrey","not_show"=NA))
poxv_synteny_dupl_gggenomes.p <- poxv_synteny_dupl_gggenomes.p + scale_color_manual("border", values = c("N0.HOG0000000"="magenta", "N0.HOG0000007"="brown", "yes"=NA, "no"=NA))
poxv_synteny_dupl_gggenomes.p <- poxv_synteny_dupl_gggenomes.p + scale_y_continuous(expand = expansion(add = c(0.5, 10)))
poxv_synteny_dupl_gggenomes.p <- poxv_synteny_dupl_gggenomes.p |> align(.justify = "center")

# pdf("Figures/poxvirus_synteny_dupl_status_boundary.pdf", width=60, height=30)
# print(poxv_synteny_dupl_gggenomes.p)
# dev.off()


########## Compare relative location of gained and pre-existing ortholog groups
##### Test difference in positions of genes 
poxv_regions %>% filter(gene_status!="gained") %>% filter(gene_status!="undeterminable") %>% group_by(gene_status) %>% shapiro_test(rel_loc) # Gain group has more than 5,000 data points
# ggplot(poxv_regions, aes(x=rel_loc, fill=gene_status)) + geom_density() + facet_grid(.~gene_status) # Unlikely normally distributed data
poxv_regions %>% filter(gene_status!="undeterminable") %>% levene_test(rel_loc~gene_status)

##### Use Kruskal-Wallis test since normality assumption is violated
poxv_rel_loc_pos_res.kruskal <- poxv_regions %>% filter(gene_status!="undeterminable") %>% kruskal_test(rel_loc~gene_status)
poxv_rel_loc_pos_res.kruskal

poxv_rel_loc_pos_stats <- poxv_regions %>% filter(gene_status!="undeterminable") %>% dunn_test(rel_loc~gene_status, p.adjust.method="fdr")
poxv_rel_loc_pos_stats <- poxv_rel_loc_pos_stats %>% add_xy_position(x="gene_status", fun="max")
poxv_rel_loc_pos_stats

########## ***** Figure 3C: Relative location of gained and pre-existing ortholog groups *****
poxv_rel_loc.p <- ggplot(poxv_regions %>% filter(gene_status!="undeterminable"), aes(x=gene_status, y=rel_loc, fill=gene_status), alpha=0.5) 
poxv_rel_loc.p <- poxv_rel_loc.p + geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.75)
poxv_rel_loc.p <- poxv_rel_loc.p + labs(title="Relative location of genes from the center of the genome", x="Gene group", y="Relative location from the center of genome")
poxv_rel_loc.p <- poxv_rel_loc.p + scale_fill_manual("gene_status", values = c("gained"="orange","pre-existing"="navyblue","undeterminable"="lightgrey"))
poxv_rel_loc.p <- poxv_rel_loc.p + theme_classic()
poxv_rel_loc.p <- poxv_rel_loc.p + theme(text = element_text(size = 18), plot.title = element_text(size = 24), plot.subtitle = element_text(size = 22))
poxv_rel_loc.p <- poxv_rel_loc.p + stat_pvalue_manual(poxv_rel_loc_pos_stats, label="p.adj", hide.ns=TRUE, inherit.aes=FALSE, size=3)
poxv_rel_loc.p <- poxv_rel_loc.p + scale_y_continuous(breaks=seq(0, 1.25, by=0.25), limits=c(0,1.25))
poxv_rel_loc.p

# pdf("Figures/relative_location_of_poxvirus_genes.pdf", width=6, height=8)
# print(poxv_rel_loc.p)
# dev.off()


########## Compare proportion of homologs of gained ortholog groups in central and terminal regions
species_order <- poxv_t.p$data %>% filter(isTip) %>% arrange(y) %>% pull(label)

poxv_regions_bdr <- poxv_regions %>% left_join(boundary.df %>% select(species, genome_max, left_boundary_bp, right_boundary_bp), by="species")

poxv_gene_test <- poxv_regions_bdr %>% mutate(region = ifelse(!is.na(left_boundary_bp) & !is.na(right_boundary_bp), ifelse(gene_mid > left_boundary_bp & gene_mid < right_boundary_bp, "central", "terminal"), NA))
poxv_gene_test <- poxv_gene_test %>% mutate(region = factor(region, levels = c("terminal", "central")), gene_status = factor(gene_status, levels = c("pre-existing", "gained", "undeterminable")), species = factor(species, levels = species_order))

prop_gained_table <- poxv_gene_test %>% group_by(species, region, gene_status) %>% summarise(n_gene_status=n(), .groups="drop")
prop_gained_table <- prop_gained_table %>% complete(species, region, gene_status, fill = list(n_gene_status = 0))
prop_gained_table <- prop_gained_table %>% group_by(species, region) %>% mutate(n_total = sum(n_gene_status), prop = n_gene_status / n_total) %>% ungroup()
prop_gained_table <- prop_gained_table %>% mutate(clade = poxv_ncbi_ds.df$clade[match(as.character(species), poxv_ncbi_ds.df$tip_label_dummy)])

prop_gained_table$gene_status <- factor(prop_gained_table$gene_status, levels=c("undeterminable","pre-existing","gained"))

########## ***** Figure 3B: Proportion of homologs of gained ortholog groups in central and terminal regions *****
gained_prop_term_cent.p <- ggplot(prop_gained_table %>% filter(gene_status=="gained"), aes(x=n_gene_status, y=species, fill=region))
gained_prop_term_cent.p <- gained_prop_term_cent.p + geom_bar(stat="identity", position="fill", width=0.5)
gained_prop_term_cent.p <- gained_prop_term_cent.p + scale_fill_manual("region", values = c("terminal"="hotpink","central"="#0080FF"))
gained_prop_term_cent.p <- gained_prop_term_cent.p + labs(title="Proportion of gained genes", x="Proportion", y="Virus")
gained_prop_term_cent.p <- gained_prop_term_cent.p + theme_classic()
gained_prop_term_cent.p <- gained_prop_term_cent.p + theme(text = element_text(size = 18), plot.title = element_text(size = 24), plot.subtitle = element_text(size = 22))
gained_prop_term_cent.p

# pdf("Figures/prop_gained_genes_in_central_terminal_regions.pdf", width=16, height=16)
# print(gained_prop_term_cent.p)
# dev.off()



############################## Estimating genome length increased due to gene gain ##############################
##### Find proportion of genome length contributed from gene gain and gene duplication
merge_intervals <- function(df) {
  df %>% arrange(species, start, end) %>% group_by(species) %>% mutate(running_end = cummax(end), new_block = start > lag(running_end, default = -Inf), block_id = cumsum(new_block)) %>%
    group_by(species, block_id) %>% summarise(start = min(start), end = max(end), .groups = "drop")}

overlapped_ranges <- poxv_regions %>% select(species, gene_ID, start, end) %>% inner_join(poxv_regions %>% select(species, gene_ID, start, end), by = "species", suffix = c("_1", "_2"), relationship = "many-to-many") %>%
  filter(gene_ID_1 < gene_ID_2) %>% filter(start_1 <= end_2, start_2 <= end_1) %>% transmute(species, start = pmax(start_1, start_2), end = pmin(end_1, end_2)) %>% merge_intervals()

per_interval_unique <- poxv_regions %>% mutate(gene_row = row_number(), interval_len = end - start + 1) %>% left_join(overlapped_ranges %>% rename(ovl_start = start, ovl_end = end), by = "species", relationship = "many-to-many") %>%
  filter(is.na(ovl_start) | (start <= ovl_end & ovl_start <= end)) %>% mutate(overlap_part = ifelse(is.na(ovl_start), 0, pmin(end, ovl_end) - pmax(start, ovl_start) + 1)) %>% group_by(gene_row) %>%
  summarise(species = first(species), gene_ID = first(gene_ID), start = first(start), end = first(end), interval_len = first(interval_len), 
            overlapped_len = sum(overlap_part), unique_len = interval_len - overlapped_len,.groups = "drop") %>% select(-gene_row)

poxv_regions_l.df <- poxv_regions %>% left_join(per_interval_unique, by=c("species","gene_ID","start","end"))
poxv_regions_l.df <- poxv_regions_l.df %>% mutate(length=end-start+1)
poxv_regions_l.df <- poxv_regions_l.df %>% mutate(non_ov_length=ifelse(is.na(unique_len), end-start+1, unique_len))
poxv_regions_l.df <- poxv_regions_l.df %>% mutate(gained_status=poxv_og_gain_loss_status.df$gained_status[match(group, gsub("N0HOG", "N0.HOG", poxv_og_gain_loss_status.df$poxv_og_name, fixed=TRUE))])
poxv_regions_l.df <- poxv_regions_l.df %>% mutate(dupl_status=ifelse(gene_ID %in% gene_dup_prot_list.v, "dupl", "not-dupl"))
poxv_regions_l.df <- poxv_regions_l.df %>% mutate(hgt_status=ifelse(gene_ID %in% hgt_gained_res.np.df$qseqid, "hgt", "non-hgt"))
poxv_regions_l.df <- poxv_regions_l.df %>% mutate(sum_status=paste(gained_status, hgt_status, dupl_status, sep="_"))

poxv_regions_l.df %>% group_by(group) %>% summarise(av_size=mean(length)) %>% arrange(desc(av_size)) # Average size of genes 

poxv_regions_l_gain_status_ct.df <- poxv_regions_l.df %>% filter(!sum_status=="pre-existing_non-hgt_not-dupl") %>% group_by(species, group, gained_status) %>% summarise(count=n(), total_length=sum(length), total_non_ov_length=sum(non_ov_length), .groups = "drop")
poxv_regions_l_gain_status_ct.df <- poxv_regions_l_gain_status_ct.df %>% mutate(added_non_ov_length=ifelse(gained_status=="pre-existing" & count > 1, total_non_ov_length-(total_length/count), total_non_ov_length)) # Remain the total non-overlapped length if the count is one
poxv_regions_l_gain_status_ct.df <- poxv_regions_l_gain_status_ct.df %>% group_by(species, gained_status) %>% summarise(total_non_ov_length=sum(added_non_ov_length))
  
summary(poxv_regions_l_gain_status_ct.df %>% group_by(species) %>% mutate(exp_length = sum(total_non_ov_length)) %>% ungroup() %>% mutate(prop_exp_length=total_non_ov_length/exp_length) %>% filter(gained_status=="gained") %>% pull(prop_exp_length))*100 # % of genomic length expanded by gained genes

poxv_regions_l_gain_status_ct.df$gained_status <- factor(poxv_regions_l_gain_status_ct.df$gained_status, levels=c("undeterminable","pre-existing","gained"))
poxv_regions_l_gain_status_ct.df$species <- factor(poxv_regions_l_gain_status_ct.df$species, levels=rev(poxv_t.p$data %>% filter(isTip) %>% arrange(desc(y)) %>% pull(label)))

########## ***** Figure 2C: Genome length contributed from gained, pre-existing, and undeterminable ortholog groups *****
poxv_added_non_ov_length.p <- ggplot(poxv_regions_l_gain_status_ct.df, aes(x=total_non_ov_length, y=species, fill=gained_status))
poxv_added_non_ov_length.p <- poxv_added_non_ov_length.p + geom_bar(stat="identity", width=0.5)
poxv_added_non_ov_length.p <- poxv_added_non_ov_length.p + scale_fill_manual("gained_status", values = c("gained"="orange","pre-existing"="navyblue","undeterminable"="lightgrey"))
poxv_added_non_ov_length.p <- poxv_added_non_ov_length.p + labs(title="Genome length contributed from gained, pre-existing, and undeterminable ortholog groups", x="Length", y="Virus")
poxv_added_non_ov_length.p <- poxv_added_non_ov_length.p + theme_classic()
poxv_added_non_ov_length.p <- poxv_added_non_ov_length.p + theme(text = element_text(size = 18), plot.title = element_text(size = 24), plot.subtitle = element_text(size = 22))
poxv_added_non_ov_length.p <- poxv_added_non_ov_length.p + scale_x_continuous(breaks=seq(0, 300000, by=75000), limits=c(0,300000))
poxv_added_non_ov_length.p

# pdf("Figures/genome_added_gained_status.pdf", width=12, height=8)
# print(poxv_added_non_ov_length.p)
# dev.off()

########## Proportion of genome length contributed from gained, pre-existing, and undeterminable ortholog groups
poxv_prop_added_non_ov_length.p <- ggplot(poxv_regions_l_gain_status_ct.df, aes(x=total_non_ov_length, y=species, fill=gained_status))
poxv_prop_added_non_ov_length.p <- poxv_prop_added_non_ov_length.p + geom_bar(stat="identity", position="fill", width=0.5)
poxv_prop_added_non_ov_length.p <- poxv_prop_added_non_ov_length.p + geom_vline(xintercept = 0.5, linetype = "dashed", color = "black")
poxv_prop_added_non_ov_length.p <- poxv_prop_added_non_ov_length.p + scale_fill_manual("gained_status", values = c("gained"="orange","pre-existing"="navyblue","undeterminable"="lightgrey"))
poxv_prop_added_non_ov_length.p <- poxv_prop_added_non_ov_length.p + labs(title="Proportion of genome length contributed from gained, pre-existing, and undeterminable ortholog groups", x="Proportion of length", y="Virus")
poxv_prop_added_non_ov_length.p <- poxv_prop_added_non_ov_length.p + theme_classic()
poxv_prop_added_non_ov_length.p <- poxv_prop_added_non_ov_length.p + theme(text = element_text(size = 18), plot.title = element_text(size = 24), plot.subtitle = element_text(size = 22))
poxv_prop_added_non_ov_length.p

# pdf("Figures/genome_added_prop_gained_status.pdf", width=12, height=8)
# print(poxv_prop_added_non_ov_length.p)
# dev.off()



############################## Analyze HGT events from gene trees ##############################
poxv_hgt_dir <- "241008_Poxviridae_and_representative_MPXV_hgt_gene_tree_top_1_mincov_90_family_trimmed"
poxv_hgt_gene_t.rt.p_list <- list()

prot_tax.f <- paste(poxv_hgt_dir, "all_init_prot_tax.tsv", sep="/")
prot_tax.df <- fread(prot_tax.f, header=T, sep="\t", quote="", check.names=T, data.table=FALSE)

poxv_all_hgt_gene_t_info.df <- data.frame(og_name = character(0), num_tip = numeric(0), num_virus_sp = numeric(0), stringsAsFactors = FALSE)

poxv_hgt_og_set_id.v <- hgt_og_set_status.df %>% filter(gained_hgt_status=="yes") %>% pull(hog_set)
poxv_hgt_og_set_id.v <- gsub(",", "_", poxv_hgt_og_set_id.v)

##### Check trees of mammalian ortholog groups acquired in the Avipoxvirus lineage
poxv_hgt_og_set_id_on_fc.v <- c(poxv_state_tr_all_m_hgt_split.df %>% filter(node %in% c(103, getDescendants(poxv_t.rt.rn, 103))) %>% filter(!is.na(host_class), host_class=="Mammalia") %>% pull(unique(hog_set)))

# > poxv_hgt_og_set_id_on_fc.v
# [1] "N0.HOG0000123" "N0.HOG0000137" "N0.HOG0000139" "N0.HOG0000149" "N0.HOG0000183" "N0.HOG0000382" "N0.HOG0000455" "N0.HOG0000972"
# 
# poxv_hgt_og_set_id_on_fc_done.v <- c("not_done_rc","not_done_rc","not_done_rc","not_done_rc","not_done_rc","not_done_rc","not_done_rc","not_done_rc")
# poxv_hgt_og_set_id_on_fc_root.v <- c(NA,NA,NA,NA,NA,NA,NA,NA)

poxv_hgt_og_set_id_on_fc_done.v <- c("done","not_done_rc","not_done_rc","done","done","done","done","done")
poxv_hgt_og_set_id_on_fc_root.v <- c(2337,NA,NA,581,1058,225,968,445)

poxv_hgt_og_set_id_on_fc.df <- data.frame(hog_set=poxv_hgt_og_set_id_on_fc.v, done=poxv_hgt_og_set_id_on_fc_done.v, root_node=as.numeric(poxv_hgt_og_set_id_on_fc_root.v))

for (a in 1:length(poxv_hgt_og_set_id.v)) {
  poxv_hgt_gene_t.f <- paste(poxv_hgt_dir, "/", poxv_hgt_og_set_id.v[a], "/", poxv_hgt_og_set_id.v[a], "_init_prot.trimal.treefile", sep="")
  if (file.exists(poxv_hgt_gene_t.f) == FALSE) {next}
  
  poxv_hgt_gene_t_info.df <- data.frame(og_name = character(0), num_tip = numeric(0), num_virus_sp = numeric(0), stringsAsFactors = FALSE)
  
  poxv_hgt_gene_t <- read.tree(poxv_hgt_gene_t.f)

  poxv_hgt_gene_t$tip.label <- unname(sapply(poxv_hgt_gene_t$tip.label, function(x) {
    parts <- str_split(x, "_", simplify = TRUE)
    num_parts <- length(parts)
    
    # Check if the second-to-last part is "NP", "XP", or "YP"
    if (parts[num_parts - 1] %in% c("NP", "XP", "YP")) {
      # Recombine all but the last three parts with underscores
      prefix <- paste(parts[1:(num_parts - 2)], collapse = "_")
      # Append an asterisk before "NP" or "YP" part and add the remaining parts
      paste(prefix, "*", paste(parts[num_parts - 1], "_", parts[num_parts], sep=""), sep="")
    } else if (TRUE %in% str_detect(x, "N0.HOG")) {
      prefix <- paste(parts[1:(which(str_detect(parts, "N0.HOG")==TRUE)-1)], collapse = "_")
      paste(prefix, "*", paste(parts[(which(str_detect(parts, "N0.HOG")==TRUE)):num_parts], collapse = "_"), sep="")
    } else {
      # For all other cases, recombine all but the last part with underscores
      prefix <- paste(parts[1:(num_parts - 1)], collapse = "_")
      # Append an asterisk before the last part
      paste(prefix, "*", parts[num_parts], sep="")
    }
  }))
  
  poxv_hgt_gene_t$tip.label <- ifelse(str_detect(poxv_hgt_gene_t$tip.label, "GCF_000857045.1_Monkeypox_virus"), gsub("GCF_000857045.1_Monkeypox_virus*", "GCF_000857045.1_Monkeypox_virus_Zaire-96-I-16_Ia*", poxv_hgt_gene_t$tip.label, fixed=TRUE),
                               ifelse(str_detect(poxv_hgt_gene_t$tip.label, "GCA_039269415.1_Monkeypox_virus"), gsub("GCA_039269415.1_Monkeypox_virus*", "GCA_039269415.1_Monkeypox_virus_RDC-NKV-GOM-MPOX-010_Ib*", poxv_hgt_gene_t$tip.label, fixed=TRUE),
                               ifelse(str_detect(poxv_hgt_gene_t$tip.label, "GCA_006465585.1_Monkeypox_virus"), gsub("GCA_006465585.1_Monkeypox_virus*", "GCA_006465585.1_Monkeypox_virus_Sierra_Leone_IIa*", poxv_hgt_gene_t$tip.label, fixed=TRUE),
                               ifelse(str_detect(poxv_hgt_gene_t$tip.label, "GCF_014621545.1_Monkeypox_virus"), gsub("GCF_014621545.1_Monkeypox_virus*", "GCF_014621545.1_Monkeypox_virus_MPXV-M5312_HM12_Rivers_IIb*", poxv_hgt_gene_t$tip.label, fixed=TRUE), poxv_hgt_gene_t$tip.label))))
  
  # poxv_hgt_gene_t.tip_to_drop <- (poxv_hgt_gene_t$tip.label)[!(str_split(poxv_hgt_gene_t$tip.label, "\\*", simplify=TRUE)[,1]) %in% host_t.rt.nolabs.combined$tip.label]
  # poxv_hgt_gene_t <- drop.tip(poxv_hgt_gene_t, poxv_hgt_gene_t.tip_to_drop)

  poxv_hgt_gene_t.rn.f <- paste(poxv_hgt_dir, "/", poxv_hgt_og_set_id.v[a], "/", poxv_hgt_og_set_id.v[a], "_init_prot.mafft.renamed.treefile", sep="")
  if (file.exists(poxv_hgt_gene_t.rn.f) == FALSE) {write.tree(poxv_hgt_gene_t, file=poxv_hgt_gene_t.rn.f, digits = 10)}
  
  poxv_hgt_gene_t.rn <- read.tree(poxv_hgt_gene_t.rn.f)
  
  poxv_hgt_gene_t.rt.dir <- paste(poxv_hgt_dir, "/", poxv_hgt_og_set_id.v[a], sep="")

  poxv_hgt_gene_t.rt.f <- poxv_hgt_gene_t.rn.f
  if (file.exists(poxv_hgt_gene_t.rt.f) == FALSE) {next}
  
  poxv_hgt_gene_t.rt <- read.tree(poxv_hgt_gene_t.rt.f)

  if (poxv_hgt_og_set_id.v[a] %in% poxv_hgt_og_set_id_on_fc.df$hog_set) {
    set_index <- which(poxv_hgt_og_set_id_on_fc.df$hog_set==poxv_hgt_og_set_id.v[a])
    if (poxv_hgt_og_set_id_on_fc.df[set_index,"done"]=="not_done") {break}
    else {
      if (!is.na(poxv_hgt_og_set_id_on_fc.df$root_node[set_index])) {
        poxv_hgt_gene_t.rt <- phytools::reroot(poxv_hgt_gene_t.rt, node=poxv_hgt_og_set_id_on_fc.df$root_node[set_index],
                                               position=0.5*poxv_hgt_gene_t.rt$edge.length[which(poxv_hgt_gene_t.rt$edge[,2]==poxv_hgt_og_set_id_on_fc.df$root_node[set_index])])
      }
    }
  }

  poxv_hgt_gene_t.rt.df <- ggtree(poxv_hgt_gene_t.rt)$data
  poxv_hgt_gene_t.rt.df <- poxv_hgt_gene_t.rt.df %>% mutate(type=ifelse(isTip, ifelse(str_detect(label, "NVNP_"), "NVNP", ifelse(str_detect(label, "VNP_"), "VNP", "POXV")), NA))
  poxv_hgt_gene_t.rt.df <- poxv_hgt_gene_t.rt.df %>% mutate(prot_id=str_split(label, "\\*", simplify=TRUE)[,2])
  poxv_hgt_gene_t.rt.df <- poxv_hgt_gene_t.rt.df %>% mutate(phylum=prot_tax.df$phylum[match(prot_id, prot_tax.df$prot_id)])
  poxv_hgt_gene_t.rt.df <- poxv_hgt_gene_t.rt.df %>% mutate(class=prot_tax.df$class[match(prot_id, prot_tax.df$prot_id)])
  poxv_hgt_gene_t.rt.df <- poxv_hgt_gene_t.rt.df %>% mutate(class=ifelse(isTip, ifelse(str_detect(label, "NVNP_"), class, ifelse(str_detect(label, "VNP_"), class, "Pokkesviricetes")), NA))
  poxv_hgt_gene_t.rt.df <- poxv_hgt_gene_t.rt.df %>% mutate(virus_name=ifelse(type=="POXV", str_split(label, "\\*", simplify=TRUE)[,1], NA))
  poxv_hgt_gene_t.rt.df <- poxv_hgt_gene_t.rt.df %>% mutate(clade=poxv_ncbi_ds.df[match(virus_name, poxv_ncbi_ds.df$tip_label), "clade"])
  
  poxv_hgt_gene_t.rt_ht.df <- poxv_hgt_gene_t.rt.df %>% filter(isTip) %>% select(label, type, clade, phylum, class) %>% column_to_rownames("label")

  poxv_hgt_gene_t_info.df[1,1] <- poxv_hgt_og_set_id.v[a]
  poxv_hgt_gene_t_info.df[1,2] <- length(poxv_hgt_gene_t.rt$tip.label)
  poxv_hgt_gene_t_info.df[1,3] <- length(unique(poxv_hgt_gene_t.rt.df %>% filter(!is.na(virus_name)) %>% pull(virus_name)))
  
  poxv_all_hgt_gene_t_info.df <- rbind(poxv_all_hgt_gene_t_info.df, poxv_hgt_gene_t_info.df)
  
  ##### Plot the results
  poxv_hgt_gene_t.rt.p <- ggtree(as.phylo(poxv_hgt_gene_t.rt), linewidth=0.5)
  poxv_hgt_gene_t.rt.p <- poxv_hgt_gene_t.rt.p + theme_tree()
  poxv_hgt_gene_t.rt.p <- poxv_hgt_gene_t.rt.p %<+% poxv_hgt_gene_t.rt.df
  poxv_hgt_gene_t.rt.p <- poxv_hgt_gene_t.rt.p + geom_tiplab(aes(subset=isTip, color=type), size=2, show.legend = FALSE) # Poxvirus
  poxv_hgt_gene_t.rt.p <- poxv_hgt_gene_t.rt.p + geom_text2(aes(subset=!isTip, label=node), size=5, hjust=1.3)
  poxv_hgt_gene_t.rt.p <- poxv_hgt_gene_t.rt.p + geom_point2(aes(subset=isTip, color=type), size=2)
  poxv_hgt_gene_t.rt.p <- poxv_hgt_gene_t.rt.p + guides(colour=guide_legend(override.aes = list(label="", shape=16), ncol=1))
  poxv_hgt_gene_t.rt.p <- poxv_hgt_gene_t.rt.p + scale_color_manual("type", values=c("POXV"="red", "NVNP"="blue", "VNP"="green"), name="Sources")
  poxv_hgt_gene_t.rt.p <- poxv_hgt_gene_t.rt.p + theme(text = element_text(size = 20), plot.title = element_text(size = 20), plot.subtitle = element_text(size = 16))
  poxv_hgt_gene_t.rt.p <- poxv_hgt_gene_t.rt.p + labs(title=poxv_hgt_og_set_id.v[a], subtitle=paste(poxv_vacv_og_set_name.df %>% filter(gsub(",", "_", hog_set)==poxv_hgt_og_set_id.v[a]) %>% pull(pooled_vv_name),
                                                                                                "\ntip_num = ", poxv_hgt_gene_t_info.df[1,2], ", sp_num = ", poxv_hgt_gene_t_info.df[1,3], sep=""))
  poxv_hgt_gene_t.rt.p <- poxv_hgt_gene_t.rt.p + geom_treescale(width=0.5)
  
  base_offset <- 0
  strip_w <- 0.05
  gap <- 0.05
  
  poxv_hgt_gene_t.rt_wt_ht.p <- gheatmap(poxv_hgt_gene_t.rt.p, poxv_hgt_gene_t.rt_ht.df[, "type", drop=FALSE], offset=base_offset, width=0.125, color=NA, colnames=TRUE, colnames_angle=90, colnames_position="top", colnames_offset_y=20, hjust=0)
  poxv_hgt_gene_t.rt_wt_ht.p <- poxv_hgt_gene_t.rt_wt_ht.p + scale_fill_manual("type", values=c("POXV"="red", "NVNP"="blue", "VNP"="green"))
  poxv_hgt_gene_t.rt_wt_ht.p <- poxv_hgt_gene_t.rt_wt_ht.p + new_scale_fill()
  poxv_hgt_gene_t.rt_wt_ht.p <- gheatmap(poxv_hgt_gene_t.rt_wt_ht.p, poxv_hgt_gene_t.rt_ht.df[, "clade", drop=FALSE], offset=base_offset+(strip_w+gap)*10, width=0.125, color=NA, colnames=TRUE, colnames_angle=90, colnames_position="top", colnames_offset_y=20, hjust=0)
  poxv_hgt_gene_t.rt_wt_ht.p <- poxv_hgt_gene_t.rt_wt_ht.p + scale_fill_manual("clade", labels=names(poxv_color.clade.v)[which(names(poxv_color.clade.v) %in% poxv_hgt_gene_t.rt_ht.df$clade)], 
                                                                               values=poxv_color.clade.v[which(names(poxv_color.clade.v) %in% poxv_hgt_gene_t.rt_ht.df$clade)], na.value = NA)
  poxv_hgt_gene_t.rt_wt_ht.p <- poxv_hgt_gene_t.rt_wt_ht.p + new_scale_fill()
  poxv_hgt_gene_t.rt_wt_ht.p <- gheatmap(poxv_hgt_gene_t.rt_wt_ht.p, poxv_hgt_gene_t.rt_ht.df[, "class", drop=FALSE], offset=base_offset+(strip_w+gap)*20, width=0.125, color=NA, colnames=TRUE, colnames_angle=90, colnames_position="top", colnames_offset_y=20, hjust=0)
  poxv_hgt_gene_t.rt_wt_ht.p <- poxv_hgt_gene_t.rt_wt_ht.p + scale_fill_discrete("class")
  poxv_hgt_gene_t.rt_wt_ht.p <- poxv_hgt_gene_t.rt_wt_ht.p + new_scale_fill()
  poxv_hgt_gene_t.rt_wt_ht.p <- gheatmap(poxv_hgt_gene_t.rt_wt_ht.p, poxv_hgt_gene_t.rt_ht.df[, "phylum", drop=FALSE], offset=base_offset+(strip_w+gap)*30, width=0.125, color=NA, colnames=TRUE, colnames_angle=90, colnames_position="top", colnames_offset_y=20, hjust=0)
  poxv_hgt_gene_t.rt_wt_ht.p <- poxv_hgt_gene_t.rt_wt_ht.p + scale_fill_discrete("phylum", na.value=NA)
  poxv_hgt_gene_t.rt_wt_ht.p <- poxv_hgt_gene_t.rt_wt_ht.p + guides(fill=guide_legend(override.aes = list(label="", shape=16), ncol=1))
  poxv_hgt_gene_t.rt_wt_ht.p <- poxv_hgt_gene_t.rt_wt_ht.p + vexpand(.01, -1)
  poxv_hgt_gene_t.rt_wt_ht.p <- poxv_hgt_gene_t.rt_wt_ht.p + coord_cartesian(clip="off")

  poxv_hgt_gene_t.rt.p_list[[poxv_hgt_og_set_id.v[a]]] <- poxv_hgt_gene_t.rt_wt_ht.p
}

########## Plot all gene trees, including
########## ***** Figure 4D: Gene trees of mammalian-derived ortholog groups (N0.HOG0000382 and N0.HOG0000455) *****
# pdf("Figures/HGT_gene_tree_plots_family_trimmed_unrooted.pdf", width=30, height=30)
# pdf("Figures/HGT_gene_tree_plots_family_trimmed_rooted.pdf", width=30, height=30) # Only ortholog groups of interest were manually rooted
# for (a in 1:length(poxv_hgt_og_set_id.v)) {
#   # if (poxv_hgt_og_set_id.v[a] %in% (hgt_og_set_status.df %>% filter(gained_hgt_status=="yes") %>% pull(hog_set))) {print(poxv_hgt_gene_t.rt.p_list[[poxv_hgt_og_set_id.v[a]]])}
#   print(poxv_hgt_gene_t.rt.p_list[[poxv_hgt_og_set_id.v[a]]])
# }
# dev.off()
