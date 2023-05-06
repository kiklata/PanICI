

# SCP--------------------------

p <- ggplot() +
  geom_point(data = df[!df[["DE"]] & !df[["border"]], ], aes(x = x, y = y, color = .data[[color_by]]), size = pt.size, alpha = pt.alpha) +
  geom_point(data = df[!df[["DE"]] & df[["border"]], ], aes(x = x, y = y, color = .data[[color_by]]), size = pt.size, alpha = pt.alpha, position = jitter) +
  geom_point(data = df[df[["DE"]] & !df[["border"]], ], aes(x = x, y = y), color = cols.highlight, size = sizes.highlight + stroke.highlight, alpha = alpha.highlight) +
  geom_point(data = df[df[["DE"]] & df[["border"]], ], aes(x = x, y = y), color = cols.highlight, size = sizes.highlight + stroke.highlight, alpha = alpha.highlight, position = jitter) +
  geom_point(data = df[df[["DE"]] & !df[["border"]], ], aes(x = x, y = y, color = .data[[color_by]]), size = pt.size, alpha = pt.alpha) +
  geom_point(data = df[df[["DE"]] & df[["border"]], ], aes(x = x, y = y, color = .data[[color_by]]), size = pt.size, alpha = pt.alpha, position = jitter) +
  geom_hline(yintercept = 0, color = "black", linetype = 1) +
  geom_vline(xintercept = 0, color = "grey", linetype = 2) +
  geom_text_repel(
    data = df[df[["label"]], ], aes(x = x, y = y, label = gene),
    min.segment.length = 0, max.overlaps = 100, segment.colour = "grey40",
    color = label.fg, bg.color = label.bg, bg.r = label.bg.r, size = label.size, force = 20,
    nudge_x = ifelse(df[df[["label"]], "y"] >= 0, -x_nudge, x_nudge)
  ) +
  labs(x = xlab, y = ylab) +
  scale_color_gradientn(
    name = ifelse(x_metric == "diff_pct", "log2FC", "diff_pct"), colors = palette_scp(palette = palette, palcolor = palcolor),
    values = rescale(unique(c(min(df[, color_by], 0), 0, max(df[, color_by])))),
    guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", title.hjust = 0, order = 1)
  ) +
  scale_y_continuous(labels = abs) +
  facet_wrap(~group1) +
  theme_scp(aspect.ratio = aspect.ratio)


# grey

marker1$gene = rownames(marker1)

library(dplyr) 
marker1 <- marker1 %>% mutate(Difference = pct.1 - pct.2) 

library(ggplot) 
library(ggrepel)

ggplot(marker1, aes(x=Difference, y=avg_log2FC)) + 
  geom_point(size=0.5, color="#999999") + 
  geom_label_repel(
    data=subset(marker1, avg_log2FC >= 1 & Difference >= 0.2 & pct.2 <= 0.05), 
    aes(label=names), label.padding = 0.1, fill="tomato2", 
    segment.size = 0.25, size=2.5)+ theme_classic() 

marker1$group=0 

for (i in 1:nrow(marker1)){ 
  if (
    marker1$avg_log2FC[i] >= 1 & marker1$Difference[i] >= 0.2 & marker1$pct.2[i] <= 0.05)
  { marker1$group[i]='up' } 
  else if(
    marker1$avg_log2FC[i] <= -1 & marker1$Difference[i] <= -0.2 & marker1$pct.1[i] <= 0.05)
  { marker1$group[i]='down' } else { marker1$group[i]='no' } 
  }

ggplot(marker1, aes(x=Difference, y=avg_log2FC)) + 
  geom_point(size=0.5,aes(color=group)) + 
  scale_color_manual(values=c('blue','grey','red'))+ 
  geom_label_repel(
    data=subset(marker1, group !='no'), aes(label=names), 
    segment.size = 0.25, size=2.5)+ geom_vline(xintercept = 0.0,linetype=2)+ 
  geom_hline(yintercept = 0,linetype=2)+ theme_classic()