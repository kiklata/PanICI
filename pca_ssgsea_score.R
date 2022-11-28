raw_BCC_tumor$response = if_else(raw_BCC_tumor$patient %in% c('su001','su002','su003','su004','su009','su010'),
                                 'Yes','No')

library(factoextra)
score.ccRCC = t(score.ccRCC)
score.ccRCC = as.data.frame(score.ccRCC)
raw_tumor = subset(raw_tumor,ICB_Response %in% c('ICB_PD','ICB_PR','ICB_SD'))
score.ccRCC.norm = score.ccRCC[colnames(raw_tumor),]
raw_tumor$response = if_else(raw_tumor$ICB_Response %in% c('ICB_PD'),'No','Yes')
pca.res <- prcomp(score.ccRCC.norm,center = T,scale. = T)

fviz_eig(pca.res, addlabels = TRUE, ylim = c(0, 50))

fviz_pca_ind(pca.res,geom.ind = "point",col.ind = raw_tumor$response,
             addEllipses = F,repel = T)


score.ccRCC.norm = as.data.frame(score.ccRCC.norm)
score.ccRCC.norm$pt = substring(rownames(score.ccRCC.norm),18)

pt.names = names(table(score.ccRCC.norm$pt))

score.mean.pt = list()

for(i in 1:length(pt.names)){
  score = filter(score.ccRCC.norm,pt == pt.names[i])
  pt.score = apply(score[,-ncol(score)], 2, mean)
  pt.score = as.data.frame(pt.score)
  score.mean.pt[[pt.names[i]]] = as.data.frame(t(pt.score))
}

score.pt = data.table::rbindlist(score.mean.pt)
rownames(score.pt) = pt.names
score.pt$response = if_else(rownames(score.pt) == 'p906','No','Yes')
pca.res <- prcomp(score.pt[,-86],center = T,scale. = F)

fviz_eig(pca.res, addlabels = TRUE, ylim = c(0, 50))

fviz_pca_ind(pca.res,col.ind = score.pt$response,
             addEllipses = F,repel = T)

