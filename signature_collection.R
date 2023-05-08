## signature collection


# 1 Melanoma --------------------------------------------------------------


# 1.1 CRMA.sig ------------------------------------------------------------

# ref: Cancer-Germline Antigen Expression Discriminates Clinical Outcome to CTLA-4 Blockade
# algorithms: Geometric mean of gene expression
CRMA.sig = c('MAGEA2', 'MAGEA3', 'MAGEA6', 'MAGEA12', 'CSAG1', 'CSAG2', 'CSAG3', 'MAGEA2B')


# 1.2 IMPRES.Sig ----------------------------------------------------------
# ref: Robust prediction of response to immune checkpoint blockade therapy in metastatic melanoma
# algorithms: Comparison of 15 gene pairs
IMPRES.sig = c('PDL-1/VISTA','CD28/CD276','CD86/OX40L','CD86/CD200','CTLA-4/OX40L','PD-1/OX40L',
                'CD80/CD137L','CD86/TIM3','CD28/CD86','CD27/PD-1','CD40/PDL-1','CD40/CD80',
                'CD40/CD28','CD40/PD-1','HVEM/CD86')


# 1.3 ImmmunCells.sig -----------------------------------------------------

# ref: A gene expression signature of TREM2hi macrophages and γδT cells predicts immunotherapy response
# algorithms: Gene expression signature of TREM2hi macrophages and γδ T cells

ImmunCells.sig = c(  "JMJD7",  "TRAF3IP2","UBE2C", "CDCA5",    "TM4SF19",   "CLNK",      "TMEM171",   "CLDN7",     "CR2",      
                     "SMEK3P", "SPC24",   "CILP2", "SYT6",     "ENTHD1",    "PRUNE2",    "ALDH1L2",   "STOML3",    "NUDT10",   
                     "KLHDC8B","FBLN1",   "FBLN2", "C6orf223", "FOXI1",     "FMOD",      "FOLR2",     "LYSMD2",    "ASPM",     
                     "NUPR1",  "PPA2",    "GPR31", "GRIA1",    "GRM7",      "APBB2",     "BIRC5",     "APOC2",     "ITGA3",    
                     "KRT4",   "LALBA",   "MATN3", "MFAP2",    "SCGB2A2",   "MKI67",     "MMP12",     "MT1G",      "MUSK",     
                     "MYL1",   "CEACAM6", "ROR1",  "LEF1",     "DUSP13",    "ZNF219",    "RASL12",    "TREM2",     "CYTL1",    
                     "MXRA8",  "MAP2K5",  "PRPH",  "CD177",    "SHZ3",     "RNASE1",    "RRM2",      "CCL18",     "SEPP1",    
                     "SPP1",   "STC1",    "TK1",   "TRPC4",    "TYMS",      "CACNG1",    "KIRREL2",   "TEAD2",     "MAEL",     
                     "STC2",   "ADAM21",  "DLK1",  "SLC16A3",  "PKDCC",     "KIAA0101",  "CDH1",      "ARSF",      "CD244",    
                     "CRTAM",  "GBP1P1",  "GIMAP4","KIR2DL4",  "LINC00243", "MYO1G",     "OTOF",      "SH2D2A",    "SHC3",     
                     "SPATA13","TDRD15",  "TUBA8", "RPL36AP41","MMP9",      "UNC80",     "NACA2",     "ZNF462",    "CORO7",    
                     "NACA3P", "DHRS9",   "ASAH2", "GDF1",     "ZNF610",    "PLA2G2D",   "EIF4A2",    "RIMS2",     "ZNF880" )


# 1.4 IMS.sig -------------------------------------------------------------

# ref: Ratio of the interferon-γsignature to the immunosuppression signature predicts anti-PD-1 therapy response in melanoma
# algorithms: Ratio of IMS scores / INF-γ scores arithmetic mean of the log2-transformed, housekeeping gene-normalized expression levels 
# note: different cancer type --> different prognostic effect
IMS.score = c('FAP', 'PDGFRB','CD163', 'SIGLEC1', 'IL10', 'CCL2', 'CCL8', 
              'CCL13', 'INHBA', 'VCAN', 'AXL', 'TWIST2','ADAM12', 'COL6A3',
              'STC1', 'ISG15', 'BCAT1', 'OLFML2B')
INF.score = c('IFNG','STAT1','CCR5','CXCL9','CXCL10','CXCL11','IDO1','PRF1','GZMA','HLA-DRA')
IMS.sig = INF.score/IMS.score 


# 1.5 IPRES.sig -----------------------------------------------------------

# ref: Genomic and Transcriptomic Features of Response to Anti-PD-1 Therapy in Metastatic Melanoma
# algorithms: Mean z-scores of 26 GSVA scores

IPRES.sig = c('JAEGER_METASTASIS_UP','MAPKi_INDUCED_EMT','LEF1_UP.V1_UP',
              'POOLA_INVASIVE_BREAST_CANCER_UP','YE_METASTATIC_LIVER_CANCER','ANASTASSIOU_MULTICANCER_INVASIVENESS_SIGNATURE',
              'CHARAFE_BREAST_CANCER_BASAL_VS_MESENCHYMAL_DN','MAHADEVAN_GIST_MORPHOLOGICAL_SWITCH','VECCHI_GASTRIC_CANCER_ADVANCED_VS_EARLY_UP',
              'LIEN_BREAST_CARCINOMA_METAPLASTIC','LU_TUMOR_VASCULATURE_UP','LU_TUMOR_ANGIOGENESIS_UP','LU_TUMOR_ENDOTHELIAL_MARKERS_UP',
              'ROY_WOUND_BLOOD_VESSEL_UP','MAPKi_INDUCED_ANGIOGENESIS','EP_BLOOD_VESS_DEVEL_DN_IN_R','WESTON_VEGFA_TARGETS_6HR',
              'WESTON_VEGFA_TARGETS_12HR','MAINA_VHL_TARGETS_UP','MS_RESP_TO_HYPOXIA_UP_IN_MAPKi_aPDL1_NR',
              'HARRIS_HYPOXIA','KARAKAS_TGFB1_SIGNALING','JEON_SMAD6_TARGETS_DN','POST_OP_WOUNDHEALING',
              'MISHRA_CARCINOMA_ASSOCIATED_FIBROBLAST_UP','MS_RESP_TO_WOUNDING_UP_IN_MAPKi_aPDL1_NR')


# 1.6 TcellExc.sig --------------------------------------------------------

# ref: A Cancer Cell Program Promotes T Cell Exclusion and Resistance to Checkpoint Blockade
# algorithms: Overall Expression

TcellExc.sig = load('TcellExc.Sig.RData')
# FINAL_Tcell_Exclusion.sig.RData

# 1.7 TRS.sig -------------------------------------------------------------

# ref: Single-Cell Transcriptomic Analysis Reveals a Tumor-Reactive T Cell Signature Associated With Clinical Outcome and Immunotherapy Response In Melanoma
# algorithms: GSVA

TRS.Sig <- c('CTLA4', 'CXCR6', 'LYST', 'CD38', 'GBP2','HLA-DRB5')


# 2 Pan-cancer ------------------------------------------------------------



# 2.1 IFNG.sig -------------------------------------------------------------

# ref: IFN-γ–related mRNA profile predicts clinical response to PD-1 blockade
# algorithms: Average gene expression

IFNG.sig = c('IFNG', 'STAT1' , 'CCR5', 'CXCL9', 'CXCL10', 'CXCL11',  'IDO1',  'PRF1',  'GZMA', 'HLA-DRA')


# 2.2 T.cell.inflamed.Sig -------------------------------------------------

# ref: IFN-γ–related mRNA profile predicts clinical response to PD-1 blockade
# algorithms: Average gene expression

T.cell.inflamed.sig = c('CD3D', 'IL2RG', 'IDO1', 'NKG7', 'CIITA', 'HLA-E', 'CD3E', 'CXCR6',
                        'CCL5', 'LAG3', 'GZMK', 'TAGAP',  'CD2', 'CXCL10','HLA-DRA', 'STAT1',
                        'CXCL13', 'GZMB')



# 2.3 PDL1.sig ------------------------------------------------------------

# ref: Safety, Activity, and Immune Correlates of Anti–PD-1 Antibody in Cancer
# algorithms: Expression of PD-L1 

PDL1.Sig <- c('PDL1','PDCD1')


# 2.4 LRRC15.CAF.sig ----------------------------------------------------------

# ref:  Single-Cell RNA  Sequencing Reveals Stromal Evolution into LRRC15 +  Myofi  broblasts as a Determinant of Patient Response to Cancer Immunotherapy 
# algorithms: eigenWeightedMean method provided by R package multiGSEA

LRRC15.CAF.Sig <- c('MMP11','COL11A1','C1QTNF3','CTHRC1','COL12A1','COL10A1','COL5A2','GJB2','THBS2','AEBP1','MFAP2','LRRC15','PLAU','ITGA11') # Alias for 'PRLHR' are:'GR3','GPR10','PrRPR'


# 2.5 Cytotoxic.Sig -------------------------------------------------------
# ref: Molecular and Genetic Properties of Tumors Associated with Local Immune Cytolytic Activity
# algorithms: geometric mean of GZMA and PRF1

CYT.sig = c('GZMA','PRF1')



# 2.6 NLRP3.sig -----------------------------------------------------------

# ref: Pan-cancer analysis of NLRP3 inflammasome with potential implications in prognosis and immunotherapy in human cancer
# algorithms: ssGSEA

NLRP3.Sig <- c('ARRDC1-AS1','CARD8','GSDMD','ATAT1','CD36','CPTP','DHX33','EIF2AK2','GBP5','NLRC3','PYDC2','SIRT2','TLR4','TLR6','USP50','APP','CASP1','HSP90AB1','MEFV','NFKB1','NFKB2','NLRP3','P2RX7','PANX1','PSTPIP1','PYCARD','RELA','SUGT1','TXN','TXNIP')


# 2.7 Stem.sig ------------------------------------------------------------

# ref: Integrated analysis of single-cell and bulk RNA sequencing data reveals a pan-cancer stemness signature predicting immunotherapy response
# algorithms: Naive-bayes
# not: can't repetition


# 2.8 IRG.Sig ----------------------------------------------------

# ref: Identification of a prognostic immune signature for cervical cancer to predict survival and response to immune checkpoint inhibitors
# algorithms: glmnet '(0.32196 * LEPR) + (−0.64921* PRLHR) + (−0.32677 *NR2F2) + (0.23573*PRL) + (0.39005*NRP1) + (0.02975*TNFRSF10B) + (0.39830*TNFRSF10A) + (0.14607*PLAU) + (−0.68625 *IFI30) + (0.38166 *ANGPTL5) + (−0.03522*IGF1)'

IRG.sig = c('LEPR', 'PRLHR', 'NR2F2', 'PRL', 'NRP1', 'TNFRSF10B', 'TNFRSF10A', 'PLAU', 'IFI30', 'ANGPTL5', 'IGF1')


# 2.9 Inflammatory.Sig ----------------------------------------------------

# ref: Gene signatures of tumor inflammation and epithelial-to-mesenchymal transition (EMT) predict responses to immune checkpoint blockade in lung cancer with high accuracy 
# algorithms: summing (unweighted) log2 Z scores 

Inflammatory.sig = c('CCL5','CCR5','CD274','CD3D','CD3E','CD8A','CIITA','CTLA4','CXCL10','CXCL11','CXCL13','CXCL9','GZMA','GZMB','HLA-DRA','HLA-DRB1',
                     'HLA-E','IDO1','IL2RG','ITGAL','LAG3','NKG7','PDCD1','PRF1','PTPRC','STAT1','TAGAP')


# 2.10 EMT.sig ------------------------------------------------------------

# ref: Gene signatures of tumor inflammation and epithelial-to-mesenchymal transition (EMT) predict responses to immune checkpoint blockade in lung cancer with high accuracy 
# algorithms: sum of the log2 Z scores of 6 established mesenchymal genes (AGER, FN1, MMP2, SNAI2, VIM, ZEB2) and subtracting the sum of the log2 Z scores of 6 established epithelial genes (CDH1, CDH3, CLDN4, EPCAM, MAL2, and ST1

E.sig = c('CDH1','CDH3','CLDN4','EPCAM','ST14','MAL2')
T.sig = c('VIM','SNAI2','SEB2','FN1','MMP2','AGER')
EMT.sig = E.sig - T.sig


# 2.11 Blood.Sig ----------------------------------------------------------

# ref: Whole-blood RNA transcript-based models can predict clinical response in two large independent clinical studies of patients with advanced melanoma treated with the checkpoint inhibitor, tremelimumab
# algorithms: 24.009-(0.8697xADAM17) + (0.7486xCDK2)-(0.5885xCDKN2A) + (0.3462xDPP4)-(0.2401xERBB2) + (1.7427xHLA-DRA) + (0.2481xICOS)-(1.1975xITGA4)-(1.0184xLARGE) + (1.1721xMYC)-(0.6531xNAB2)-(1.1491xNRAS) + (0.7377xRHOC)-(1.0585xTGFB1) + (0.8328xTIMP1)

Blood.sig = c('ADAM17', 'CDK2', 'CDKN2A', 'DPP4', 'ERBB2', 'HLA-DRA', 'ICOS', 'ITGA4', 'LARGE', 'MYC', 'NAB2', 'NRAS', 'RHOC', 'TGFB1', 'TIMP1')


# 3 Other sig -------------------------------------------------------------

# ref: library(hacksig)

IFNG.6.sig = c('CXCL9', 'CXCL10', 'IDO1', 'IFNG', 'HLA-DRA', 'STAT1')