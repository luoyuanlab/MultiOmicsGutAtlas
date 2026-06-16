# CITE-seq: AUCell Annotation and Treg Pseudotime Bin Transfer
#
# Sections:
# 1. Build Seurat object from RNA + ADT matrices
# 2. AUCell myeloid subtype annotation (second refined gene sets)
# 3. AUCell T cell subtype annotation (CD atlas reference gene sets)
# 4. Treg pseudotime bin transfer via Seurat label transfer
#
# Inputs:
#   RNA_raw_counts.csv        – cell × gene raw count matrix
#   protein_raw_counts.csv    – cell × protein raw count matrix
#   annotation_table.csv      – per-cell RNA_annotation / Protein_annotation
#   myeloid/norm_counts.csv   – normalized myeloid RNA counts (from 02_cite_seq_mac.ipynb)
#   myeloid/myeloid_subtypes_third_round_AUC.csv – first-round AUCell myeloid labels
#   tcell/good_anno_t_wtcr_norm_counts.csv       – normalized T cell counts (from 03_cite_seq_treg_tcr.ipynb)
#   treg_pt_bin.RDS           – CD scRNA-seq Treg pseudotime reference object
#   tcell/notstrict_treg_raw_counts.csv          – CITE-seq Treg raw counts for PT transfer

library(Seurat)
library(ggplot2)
library(dplyr)
library(data.table)
library(AUCell)

# ── Paths ──────────────────────────────────────────────────────────────────────
DATA_DIR   <- "/path/to/cite-seq/output"
SCRNA_DIR  <- "/path/to/cd/scrna/output/merged_cd"
MYE_DIR    <- file.path(DATA_DIR, "myeloid")
TCELL_DIR  <- file.path(DATA_DIR, "tcell")

# ── 1. Build Seurat object ─────────────────────────────────────────────────────
rna        <- read.csv(file.path(DATA_DIR, "RNA_raw_counts.csv"),     row.names = 1)
adt        <- read.csv(file.path(DATA_DIR, "protein_raw_counts.csv"), row.names = 1)
annotation <- read.csv(file.path(DATA_DIR, "annotation_table.csv"),   row.names = 1)

cite <- CreateSeuratObject(counts = t(rna))
cite[["ADT"]] <- CreateAssay5Object(counts = t(adt))
cite <- AddMetaData(cite, metadata = annotation)

# Subset and save myeloid cells
myeloid <- subset(cite, subset = RNA_annotation == "myeloid")
saveRDS(myeloid, file.path(MYE_DIR, "myeloid_cite_anno.RDS"))
saveRDS(cite,    file.path(DATA_DIR, "cite_anno.RDS"))

# ── 2. AUCell myeloid subtype annotation (refined gene sets) ──────────────────
# Load normalized myeloid counts (output of 02_cite_seq_mac.ipynb → norm_counts.csv)
myeloid_norm <- as.data.frame(fread(file.path(MYE_DIR, "norm_counts.csv")))
rownames(myeloid_norm) <- myeloid_norm[[1]]
myeloid_norm <- myeloid_norm[, -1]

# Gene sets (refined, based on normalized counts)
Mac_S_SG      <- c('SIGLEC1','F13A1','LILRB5','CD209','CD163L1','MSR1','STAB1','TREM2','ATP6V0D2','ADORA3','MS4A4A','FOLR2','C1QB','C1QA','FPR3','C1QC','SDS','OTOA','CD163','MS4A7','MRC1','VSIG4','SLCO2B1','GPNMB','CMKLR1','CSF1R','MERTK','CCL18','MS4A6A','GFRA2','ADAP2','SLC40A1','SELENOP','SLC7A8','APOC1','DNASE1L3','RAB42','MPEG1','LACC1','PLA2G7','FUCA1','LIPA','DAPK1','LGMN','ITGA9','DAB2','VMO1','MMP12','IGF1','TTYH3')
Trans_cDC2_Mac <- c('CLEC10A','CD1E','VSIG4','FCER1A','CPVL','CSF1R','LYZ','SDS','MS4A6A','C1QC','C1QA','DNASE1L3','C1QB','P2RY6','FPR3','MPEG1','PLD4','CD33','HLA-DQB2','CSF2RA','SLAMF8','CFP','CST3','MRC1','AIF1','CD209','VMO1','HLA-DRB5','HLA-DRB1','CSTA','HLA-DPB1','HLA-DPA1','HLA-DRA','FGL2','MS4A7','SEMA6B','TYROBP','HLA-DQA1','STAB1','PID1','DAPK1','LST1','SLC8A1','PILRA','LRRC25','ADAP2','CLEC7A','MAFB','PLA2G7','CD1D')
CCR7_DC        <- c('CCL22','LAMP3','HMSD','CD1E','TVP23A','FSCN1','EBI3','IDO1','TIFAB','ARHGAP22','GPR157','CCL19','CSF2RA','SLCO5A1','CD80','RASSF4','CFP','CSTA','CXCL9','BCL2L14','RAMP1','SLC7A11','TXN','LAD1','CLIC2','KYNU','IL4I1','VMO1','NRP2','CD86','CD274','CCR7','LGALS2','MARCKSL1','HLA-DQB2','PTGIR','TNFAIP2','TUBB6','TSPAN33','ARHGAP10','CD200','HLA-DPB1','DAPP1','TBC1D4','FAM49A','CCDC28B','LY75','GPR137B','BIRC3','POGLUT1')
cDC2           <- c('CD207','CD1E','CACNA2D3','CLEC10A','FCER1A','CFP','CD1C','CPVL','CSF2RA','HLA-DQB2','CD33','PLD4','LGALS2','MYCL','LYZ','MS4A6A','HLA-DPB1','CD1D','LST1','CD86','CSTA','HLA-DQA1','CST3','HLA-DPA1','HLA-DRA','HLA-DRB5','HLA-DRB1','HLA-DQB1','CLEC4A','IFI30','ZNF385A','SLC8A1','P2RY6','FAM110A','SPI1','VSIG4','DNASE1L3','CSF1R','PID1','HLA-DQA2','SLAMF8','P2RY13','AIF1','FPR3','EREG','HCK','FGR','HLA-DMB','IDO1','TNFAIP8L2')
Inf_Mo_Mac     <- c('MARCO','CD300E','IL27','SUCNR1','CXCL11','LILRB2','CLEC5A','IDO1','LILRA1','LILRA5','INHBA','SLC16A10','CXCL10','RETN','APOBEC3A','CXCL9','PTX3','EREG','TNIP3','AC005224.2','RNF144B','SLC1A3','LDLRAD3','LILRB1','CD274','FCN1','CALHM6','CLEC4D','DOCK4','LILRA6','KYNU','GPR84','ADGRE2','CLEC4E','EPB41L3','TLR8','IL6','ANKRD22','CD163','CCL3L1','MIR3945HG','LILRB4','MMP19','IL1A','ATP13A3','OLR1','IL1RN','RIN2','ADGRE1','DSE')
Mo_Mac         <- c('RETN','FCN1','CD300E','EREG','LILRA1','CSTA','LILRA5','AC020656.1','OLR1','MCEMP1','APOBEC3A','ASGR1','VCAN','LILRB2','NLRP3','CFP','FCGR1A','SLC11A1','CLEC12A','LYZ','IFI30','S100A9','LILRA2','C5AR1','CLEC4E','HRH2','CD163','HBEGF','S100A12','CYBB','CSF1R','LRRC25','IL1B','EMILIN2','LILRB1','AIF1','THBS1','ADGRE2','CALHM6','MPEG1','CPVL','EPB41L3','FGR','DMXL2','CD93','TLR2','PILRA','KYNU','SEMA6B','CXCL2')
Mac_S_XS       <- c('C1QC','MMP12','PLA2G7','C1QB','C1QA','MS4A4A','DNASE1L3','APOC1','MS4A6A','FUCA1','AIF1','SELENOP','CD68','GPNMB','MS4A7','SLC40A1','STAB1','LYZ','CD14','LIPA','C1orf54','TYROBP','CPVL','IGSF6','LGMN','FCER1G','RNASE1','ACP5','MPEG1','FTL','DAB2','IFI30','CTSL','CTSB','FGL2','SPI1','RNASE6','CTSD','CTSZ','LST1','HLA-DPA1','HLA-DPB1','CREG1','HLA-DRA','HLA-DQA1','HLA-DRB1','HLA-DRB5','HLA-DQA2','PSAP','HLA-DMB')
Cycl_Myeloid   <- c('CDK1','SPC25','UBE2C','RRM2','CEP55','AURKB','GTSE1','CENPA','TYMS','CCNA2','BIRC5','TOP2A','SHCBP1','UHRF1','TK1','CDCA3','RAD51AP1','CCNB2','CDCA5','HMMR','CDC20','C1QC','ZWINT','TPX2','C1QB','ASF1B','MKI67','CDKN3','VSIG4','C1QA','MS4A4A','FAM111B','NUSAP1','MS4A6A','PRC1','ADORA3','CLSPN','CCNB1','SDS','LYZ','DNASE1L3','CDT1','P2RY6','STMN1','AIF1','CENPF','UBE2T','CSF1R','ADAP2','MMP12')
Mac_S_M_S      <- c('MMP12','CCL18','MMP9','CXCL5','CLEC5A','SLC7A11','SPP1','SLC1A3','MMP10','CXCL9','PLA2G7','MRC1','TREM2','CD163','MSR1','FPR3','LILRB4','SLAMF8','CTSL','CD68','TNFSF13','BCAT1','SLC16A10','CTSB','OLFML2B','AC020656.1','NR1H3','VMO1','INHBA','CXCL10','APOC1','MMP19','CSF1R','TFEC','CD80','CD209','CALHM6','CMKLR1','GPNMB','C1QB','FCGR1A','SUCNR1','HMOX1','DOCK4','LILRB5','LHFPL2','MS4A4A','FCGR3A','SH3PXD2B','EPB41L3')
cDC1           <- c('CLEC9A','XCR1','SERPINF2','CPVL','FLT3','IDO1','TACSTD2','C1orf54','LGALS2','DNASE1L3','BATF3','SHTN1','S100B','CST3','MYCL','HLA-DPB1','ENPP1','WDFY4','HLA-DQB2','TNNI2','VMO1','CLNK','P2RY6','SLAMF8','HLA-DPA1','LYZ','FKBP1B','MPEG1','CCND1','IRF8','CSF2RA','RAB32','HLA-DRB5','HLA-DQA1','HLA-DRB1','HLA-DQB1','HLA-DRA','SPI1','HLA-DQA2','VAC14','HCK','PLCD1','GCSAM','CPNE3','SNX3','ASB2','CLIC2','CD74','PPT1','BID')
Mac_S_M_P      <- c('PLA2G2D','IL22RA2','TMEM163','MMP9','LILRB4','TIFAB','DNASE1L3','ENPP2','PDE6G','CSTA','CCL18','RAB42','AC020656.1','FUCA1','APOC1','IL18','MMP12','PLA2G7','PTGDS','CXCL9','GPNMB','NPL','NR1H3','FOLR2','CD68','LYZ','KLHDC8B','LGALS2','SLAMF8','ATOX1','C1QC','SELENOP','CD14','APOE','CSF1R','MS4A6A','MS4A4A','RASSF4','CYP27A1','CSF2RA','C1QA','C1orf54','IL4I1','LINC00996','AIF1','CTSL','C1QB','GPR34','CUL9','SLCO2B1')
Neutrophil     <- c('FCGR3B','CMTM2','CXCR2','HCAR3','CXCR1','PROK2','FFAR2','ADGRG3','CSF3R','AL034397.3','FPR2','HCAR2','KCNJ15','CXCL8','SMIM25','AZIN1-AS1','LUCAT1','AQP9','S100A8','FPR1','TREM1','S100A12','CYP4F3','IL1RN','MIR3945HG','G0S2','S100A9','KCNJ2','LILRB3','AC015912.3','MNDA','DGAT2','AC023157.3','TLR4','PLEK','P2RY13','MMP25','ACSL1','MME','TLR2','IL1B','ALPL','LRRK2','NAMPT','STEAP4','NCF2','TNFAIP6','OSM','C5AR1','CLEC4E')

gut_set_myeloid <- GeneSetCollection(list(
  Mac_S_SG       = GeneSet(Mac_S_SG,       setName = "Mac S+SG+"),
  Trans_cDC2_Mac = GeneSet(Trans_cDC2_Mac, setName = "Trans cDC2/Mac"),
  CCR7_DC        = GeneSet(CCR7_DC,        setName = "CCR7+ DC"),
  cDC2           = GeneSet(cDC2,           setName = "cDC2"),
  Inf_Mo_Mac     = GeneSet(Inf_Mo_Mac,     setName = "Inf Mo-Mac"),
  Mo_Mac         = GeneSet(Mo_Mac,         setName = "Mo-Mac"),
  Mac_S_XS       = GeneSet(Mac_S_XS,       setName = "Mac S+XS-"),
  Cycl_Myeloid   = GeneSet(Cycl_Myeloid,   setName = "Cycl Myeloid"),
  Mac_S_M_S      = GeneSet(Mac_S_M_S,      setName = "Mac S+M+S+"),
  cDC1           = GeneSet(cDC1,           setName = "cDC1"),
  Mac_S_M_P      = GeneSet(Mac_S_M_P,      setName = "Mac S+M+P+"),
  Neutrophil     = GeneSet(Neutrophil,     setName = "Neutrophil")
))

myeloid_so  <- CreateSeuratObject(counts = t(myeloid_norm))
counts      <- GetAssayData(myeloid_so, slot = "counts")
gut_sub     <- subsetGeneSets(gut_set_myeloid, rownames(counts))

set.seed(123)
cells_AUC     <- AUCell_run(counts, gut_sub)
save(cells_AUC, file = file.path(MYE_DIR, "pr_rna_myeloid_aucell.RData"))

cells_AUC_mat <- t(getAUC(cells_AUC))
cell_labels   <- colnames(cells_AUC_mat)[max.col(cells_AUC_mat, "first")]
result        <- data.frame(cells = rownames(cells_AUC_mat), label = cell_labels)
print(table(result$label))
write.csv(result, file.path(MYE_DIR, "pr_rna_myeloid_aucell.csv"))

# ── 3. AUCell T cell subtype annotation (CD atlas reference) ──────────────────
tcell_norm <- as.data.frame(fread(file.path(TCELL_DIR, "good_anno_t_wtcr_norm_counts.csv")))
rownames(tcell_norm) <- tcell_norm[[1]]
tcell_norm <- tcell_norm[, -1]

CD4_Effector   <- c('CCL20','IL4I1','CCR6','CXCR6','LAG3','FURIN','AC058791.1','RORA','CD40LG','KLRB1','TNFRSF25','PRR5','SPOCK2','AC092580.4','CTLA4','TMIGD2','IL12RB1','PBX4','TRBC1','PDE4D','IL7R','CD2','MAF','SH2D2A','GPR171','CD3D','JAML','TNFRSF4','SLAMF1','CD6','CD3G','ADAM19','BHLHE40','MGAT4A','SYTL2','TNFRSF18','PHTF2','TTC39C','CD96','CD3E','STAT4','PBXIP1','SIRPG','AQP3','CD4','BATF','CD28','ATP2B4','TNFSF13B','ICOS')
CD4_Trm        <- c('CD40LG','IL7R','KLRB1','ANXA1','PTGER2','GATA3','GPR171','TRAT1','CD3E','PTGER4','RORA','TRAC','TTC39C','CD2','CD3D','FLT3LG','CD69','TNFAIP3','STAT4','CD96','CD6','ID2','TC2N','GPR183','ICOS','SYTL3','OXNAD1','TRBC2','TRBC1','RCAN3','KLF6','PHLDA1','SH2D2A','RGCC','CAMK4','ODF2L','PARP8','IL32','LYAR','MYADM','LCK','CD3G','ADGRE5','CLEC2D','DOK2','SIT1','LEPROTL1','TSC22D3','S100A4','GIMAP7')
GZMK_Neg_CD8_Trm <- c('XCL1','XCL2','GPR15','CD160','CD8B','CD8A','CCL5','KLRD1','HOPX','LDLRAD4','SLA2','GZMB','IFNG','TRGC2','THEMIS','GZMA','PTPN22','AOAH','AC092580.4','ITGA1','SYTL3','CD96','TMIGD2','CCL4','MT-ATP8','CXCR3','SH2D2A','CD7','CRY1','GABARAPL1','SCML4','CTSW','PRMT9','CD3E','PDE4A','PDE4D','GNG2','CBLB','RUNX3','PTGER4','ARHGAP9','TNFAIP3','CXCR6','STAT4','NKG7','CAMK4','GPR171','SPRY1','LINC-PINT','GATA3')
Treg           <- c('FOXP3','RTKN2','IL2RA','LAIR2','CTLA4','TNFRSF4','BATF','TBC1D4','ICA1','TNFRSF18','MIR4435-2HG','TIGIT','ZNRF1','CXCR6','LAG3','TRBC1','SIRPG','TRIB2','ICOS','TNFRSF1B','MAF','IL32','CD28','SPOCK2','PHTF2','MAGEH1','PBXIP1','RASGRP1','IL2RB','GBP5','TNFRSF25','CD3D','TRBC2','CD2','GBP2','DUSP4','CD247','SKAP1','SLA','CD3E','SLAMF1','TRAC','CD4','PYHIN1','LAT','UCP2','FAS','CD7','LTB','SIT1')
Gd_T_Vd2g9    <- c('FGFBP2','CX3CR1','S1PR5','GZMH','FCRL6','TRGC1','TRDC','NKG7','FCGR3A','GNLY','PRF1','KLRF1','KLRD1','SAMD3','KLRG1','GZMB','CTSW','TRGC2','PATL2','CST7','MATK','LINC00861','CLIC3','GZMM','MYO1F','CCL5','CD8B','C12orf75','ITGB2','CMC1','HOPX','CD8A','APOBEC3G','MYO1G','LYAR','CDC25B','GZMA','P2RY8','HCST','SPN','AOAH','BIN2','APMAP','FLNA','CCL4','PYHIN1','SELPLG','SPON2','ZAP70','DGKZ')
GZMK_CD8_Trm  <- c('GZMK','CRTAM','GZMH','DTHD1','CST7','CD8B','CCL4L2','CD8A','NKG7','SAMD3','CCL4','OASL','SH2D1A','KLRG1','IFNG','CCL5','XCL2','CMC1','AOAH','GZMM','APOBEC3G','GZMA','CXCR3','GZMB','PIK3R1','TIGIT','PRF1','LAG3','LYST','DUSP4','DUSP2','CTSW','TUBA4A','CD3E','PYHIN1','ADGRE5','NFATC2','RNF19A','ITGAL','CBLB','FYN','LYAR','PSTPIP1','RNF125','ARAP2','PRKCH','HNRNPLL','MT-ATP8','SRRT','HCST')
CD8_Effector  <- c('IL26','TNFRSF9','IL23R','LAYN','ATP8B4','CXCR6','SRGAP3','NCR3','AC092580.4','FASLG','CD8A','GZMB','GZMA','DUSP4','GNLY','LRRN3','HAVCR2','TMIGD2','PRF1','CCL20','MATR3.1','GEM','SIRPG','CBLB','UQCRHL','SH2D2A','LDLRAD4','ASB2','PDE4D','PTPN22','GTF3C1','LINC-PINT','SLFN12L','LAG3','IL18RAP','JAML','IL21R','CD226','MAST4','CASS4','CTLA4','MT-ATP8','KLRD1','LAX1','SLA2','CHRM3-AS2','CCR6','CD96','DUSP16','IKZF3')
MAIT          <- c('SLC4A10','NCR3','KLRG1','ZBTB16','KLRB1','PRR5','LYAR','TRGC2','CTSW','PRF1','GZMK','SAMD3','TMIGD2','AC092580.4','NKG7','MATK','MYO1F','GPR171','IL7R','HOPX','AQP3','TC2N','GIMAP1','HCST','TNFRSF25','BIN2','LTB','CD3E','GZMA','PCED1B-AS1','S100A4','CST7','SPOCK2','RORA','GZMM','ANXA2R','ZAP70','ARL4C','GIMAP4','S1PR4','CCL5','SH2D1A','FLT3LG','PTPRCAP','C12orf75','GYG1','PARP8','PTPN4','APOBEC3G','TGFB1')
Resting_T     <- c('LEF1','CCR7','MAL','TCF7','LINC00861','LDLRAP1','SELL','PIK3IP1','CD3E','TRAT1','GIMAP1','BCL11B','NOSIP','FLT3LG','GIMAP7','ITK','GIMAP4','CAMK4','GIMAP2','RCAN3','SATB1','LDHB','TRBC1','AAK1','TRAC','CD6','OXNAD1','TRAF3IP3','CD3D','IFITM1','IKZF1','RASGRP2','C1orf162','ICOS','GIMAP5','LAT','CCND3','IL7R','SARAF','DENND2D','LEPROTL1','CD3G','TRBC2','LTB','LCK','SLC2A3','ZAP70','RPS29','ATM','FXYD5')
Tfh           <- c('TOX2','ICA1','PDCD1','PASK','FBLN7','PGM2L1','TBC1D4','TIGIT','CTLA4','TOX','MAGEH1','TRBC1','MAF','TSHZ2','SESN3','ICOS','ITM2A','LAT','CD40LG','SH2D1A','IL6R','FKBP5','RNF19A','MAL','CORO1B','CD4','BATF','NFATC1','TCF7','SIRPG','TNFRSF4','CD28','IKZF3','DGKA','FYB1','SPOCK2','INPP4B','SMCO4','BCL11B','NR3C1','FCMR','GPRIN3','TRBC2','CD3D','PYHIN1','PHACTR2','CD3E','CD5','SFXN1','CD247')
Gd_T          <- c('KIR2DL4','KLRC3','KLRC2','TRDC','ATP8B4','GZMA','CD160','TRGC1','GNLY','TRGC2','XCL1','TMIGD2','LAYN','IKZF2','HOPX','FASLG','SLA2','XCL2','KLRD1','LDLRAD4','CCL5','GPR15','CD7','IL2RB','ITGA1','CD96','AOAH','TIGIT','CTSW','TNFRSF9','GZMB','GEM','PTPN22','AC092580.4','SIRPG','MATK','PRF1','ASB2','CLIC3','RIN3','ABCB1','SH2D2A','PIK3R1','CD8A','TXK','CBLB','NKG7','DUSP4','IL21R','PDE3B')
Cycl_T        <- c('DLGAP5','HIST1H1B','ASPM','MKI67','RRM2','PCLAF','CENPA','MCM10','KIF15','TOP2A','CDCA2','KIF2C','HJURP','CKAP2L','CDCA5','PKMYT1','CEP55','SPC25','CCNA2','NCAPG','GTSE1','TPX2','TYMS','UBE2C','BIRC5','CDC45','KIF11','CENPF','NUF2','NUSAP1','AURKB','ESCO2','CCNB2','HMMR','UHRF1','CDCA8','DTL','CLSPN','DEPDC1B','CENPE','SPC24','MND1','CDK1','PLK1','PBK','CDCA3','CDT1','ASF1B','TROAP','CENPU')
DN_T          <- c('MCUB','CRYBG1','NOP53','FYB1','LINC01871','PAXX','LIME1','RACK1','REX1BD','NSD3','METTL26','TRIR','ATP5PB','IL7R','SELENOT','ATP5MC2','SELENOF','FLT3LG','IFITM1','ATP5MG','ATP5F1B','ATP2B1-AS1','SELENOK','MT-ATP8','VSIR','CRIP1','TENT5C','NDUFAF8','RTRAF','MESD','LAT','NSMCE3','ELOC','SMIM26','STMP1','ATP5F1D','CDK17','ATP5F1C','ATP5F1A','ATP5PD','MT-ND4L','JPT1','ATP5MD','SELENOH','SARAF','INPP4B','COPS9','ADGRE5','RTF2','CD3E')
NK            <- c('NKG7','KLRB1','CD247','GZMM','TRBC2','TRAC','CTSW','CD7','GIMAP7','LCK','CCL5','CST7','CD2','IL7R','FYN','HCST','GZMA','RORA','SPOCK2','GIMAP4','EVL','CD96','CD3G','IL32','ACAP1','LEPROTL1','ARL4C','DUSP2','PTGER4','CCL4','ITM2A','TRAF3IP3','CREM','SOCS1','PTPN7','RGCC','PPP2R5C','TBC1D10C','ID2','CCND3','PCED1B-AS1','PTPRCAP','LTB','TNFAIP3','CD3D','EML4','ITGB2','CXCR4','CORO1A','ARHGEF1')

tgsc <- GeneSetCollection(list(
  CD4_Effector     = GeneSet(CD4_Effector,     setName = "CD4_Effector"),
  CD4_Trm          = GeneSet(CD4_Trm,          setName = "CD4_Trm"),
  GZMK_Neg_CD8_Trm = GeneSet(GZMK_Neg_CD8_Trm, setName = "GZMK_Neg_CD8_Trm"),
  Treg             = GeneSet(Treg,             setName = "Treg"),
  Gd_T_Vd2g9      = GeneSet(Gd_T_Vd2g9,      setName = "Gd_T_Vd2g9"),
  GZMK_CD8_Trm     = GeneSet(GZMK_CD8_Trm,    setName = "GZMK_CD8_Trm"),
  CD8_Effector     = GeneSet(CD8_Effector,     setName = "CD8_Effector"),
  MAIT             = GeneSet(MAIT,             setName = "MAIT"),
  Resting_T        = GeneSet(Resting_T,        setName = "Resting_T"),
  Tfh              = GeneSet(Tfh,              setName = "Tfh"),
  Gd_T             = GeneSet(Gd_T,            setName = "Gd_T"),
  Cycl_T           = GeneSet(Cycl_T,          setName = "Cycl_T"),
  DN_T             = GeneSet(DN_T,            setName = "DN_T"),
  NK               = GeneSet(NK,              setName = "NK")
))

tcell_so  <- CreateSeuratObject(counts = t(tcell_norm))
counts    <- GetAssayData(tcell_so, slot = "counts")
tgut_sub  <- subsetGeneSets(tgsc, rownames(counts))

set.seed(123)
cells_AUC     <- AUCell_run(counts, tgut_sub)
save(cells_AUC, file = file.path(TCELL_DIR, "pr_rna_tcr_tcell_aucell.RData"))

cells_AUC_mat <- t(getAUC(cells_AUC))
cell_labels   <- colnames(cells_AUC_mat)[max.col(cells_AUC_mat, "first")]
result        <- data.frame(cells = rownames(cells_AUC_mat), label = cell_labels)
print(table(result$label))
write.csv(result, file.path(TCELL_DIR, "pr_rna_tcr_tcell_aucell.csv"))

# ── 4. Treg pseudotime bin transfer (Seurat label transfer) ───────────────────
# Reference: CD scRNA-seq Treg object with pt_bin (B0–B4) from cd/12_cd_treg_pseudotime.ipynb
scrna <- readRDS(file.path(SCRNA_DIR, "treg_pt_bin.RDS"))

# Query: CITE-seq Treg raw counts (exported from 03_cite_seq_treg_tcr.ipynb)
tcr_counts <- read.csv(file.path(TCELL_DIR, "notstrict_treg_raw_counts.csv"), row.names = 1)
tcr_counts <- t(tcr_counts)

tcr <- CreateSeuratObject(counts = as.matrix(tcr_counts), assay = "RNA")
tcr <- NormalizeData(tcr)
tcr <- FindVariableFeatures(tcr, selection.method = "vst", nfeatures = 2000)
tcr <- ScaleData(tcr)
tcr <- RunPCA(tcr, npcs = 30, verbose = FALSE)
tcr <- RunUMAP(tcr, reduction = "pca", dims = 1:30)

# Gene overlap QC
common_genes <- intersect(rownames(scrna), rownames(tcr))
cat("Common genes:", length(common_genes), "\n")

# Find anchors and transfer pt_bin labels
transfer_anchors <- FindTransferAnchors(
  reference           = scrna,
  query               = tcr,
  dims                = 1:30,
  reference.reduction = "pca",
  normalization.method = "LogNormalize"
)
cat("Anchors found:", nrow(transfer_anchors@anchors), "\n")

scrna$pt_bin <- as.character(scrna$pt_bin)
predictions  <- TransferData(
  anchorset        = transfer_anchors,
  refdata          = scrna$pt_bin,
  dims             = 1:30,
  weight.reduction = tcr[["pca"]]
)
tcr <- AddMetaData(tcr, metadata = predictions)
saveRDS(tcr, file.path(TCELL_DIR, "notstrict_treg_ptbin_anno.RDS"))

cat("Prediction score summary:\n")
print(summary(tcr$prediction.score.max))
print(table(tcr$predicted.id))

# QC plots
ggplot(tcr@meta.data, aes(x = prediction.score.max)) +
  geom_histogram(bins = 50, fill = "steelblue") +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "red") +
  labs(title = "Label Transfer Prediction Scores",
       x = "Max Prediction Score", y = "Number of Cells") +
  theme_classic()
ggsave(file.path(TCELL_DIR, "tregs_prediction_scores_histogram.pdf"), width = 8, height = 6)

ggplot(tcr@meta.data, aes(x = predicted.id, y = prediction.score.max)) +
  geom_boxplot(fill = "steelblue") +
  labs(title = "Prediction Confidence by Pseudotime Bin",
       x = "Predicted Bin", y = "Prediction Score") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(TCELL_DIR, "tregs_prediction_scores_by_bin.pdf"), width = 8, height = 6)

# Export
write.csv(tcr@meta.data, file.path(TCELL_DIR, "notstrict_treg_pt_transfer_labels.csv"))
