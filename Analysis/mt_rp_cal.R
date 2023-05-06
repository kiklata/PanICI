seu = raw_cohort1_tumor
seu$percent.mt = PercentageFeatureSet(seu,pattern = '^MT-')
seu$percent.rp = PercentageFeatureSet(seu,pattern = '^RP[SL]')+
  PercentageFeatureSet(seu,pattern = '^MRP[SL]')

seu = subset(seu,expansion !='n/a')
p1 = VlnPlot(seu,'percent.rp',split.by = 'expansion',group.by = 'timepoint',pt.size = 0)
p2 = VlnPlot(seu,'percent.mt',split.by = 'expansion',group.by = 'timepoint',pt.size = 0)

seu = raw_cohort2_tumor

p3 = VlnPlot(seu,'percent.rp',split.by = 'expansion',group.by = 'timepoint',pt.size = 0)
p4 = VlnPlot(seu,'percent.mt',split.by = 'expansion',group.by = 'timepoint',pt.size = 0)

seu = raw_BCC_tumor
raw_BCC_tumor$response = if_else(raw_BCC_tumor$patient %in% c('su001','su002','su003','su004','su009','su010'),
                                 'Yes','No')

p5 = VlnPlot(seu,'percent.rp',split.by = 'response',group.by = 'treatment',pt.size = 0)
p6 = VlnPlot(seu,'percent.mt',split.by = 'response',group.by = 'treatment',pt.size = 0)


icb = subset(raw_tumor, ICB_Exposed == 'ICB')
icb$response = if_else(icb$ICB_Response %in% c('ICB_PR','ICB_SD'),'Yes',
                       if_else(icb$ICB_Response == 'ICB_PD','No','NE'))
icb = subset(icb,response !='NE')
seu = icb
p7 = VlnPlot(seu,'percent.rp',group.by = 'response',pt.size = 0)
p8 = VlnPlot(seu,'percent.mt',group.by = 'response',pt.size = 0)


