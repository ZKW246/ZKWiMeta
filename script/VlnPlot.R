
# input file
seu_data <- readRDS(seurat_file)

#  Gene expression differences between groups using violin plots
VlnPlot(seu_data, features = gene,pt.size = 0,
        group.by = 'Type', raster = FALSE, cols = color) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 12),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", size = 0.8, fill = NA, linetype = "solid"),  
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        axis.title.x = element_blank(),
        # axis.title.y = element_blank(),
        axis.title.y = element_text(size = 14,face = 'bold'),
        axis.line = element_line(colour = "black", size = 0.8),  
        axis.ticks = element_line(colour = "black", size = 0.5),  
        axis.title = element_text(size = 12, face = "bold"),  
        axis.text = element_text(size = 12, face = "bold")  
  ) +
  stat_compare_means(aes(group = seu_data$Type),
                     method = "wilcox.test",
                     label = "p.format",label.x=1.3) + 
  NoLegend() +
  scale_fill_manual(values = color)