
data = read.table()
cluster_rows = T
cluster_cols = F
scale = "row"
treeheight_row = 40
treeheight_col = 40
color = colorRampPalette(c("navy","#366dbf", "#f1f7f3", "#e6412c","firebrick3"))(1000)
cellheight = 100
cellwidth = 100
fontsize_row = 100
fontsize_col = 100
annotation_row = NA
annotation_col = NA
border_color = F
show_rownames = T
show_colnames = T

#### heatmap作图函数使用####
pheatmap::pheatmap(data,
treeheight_row = treeheight_row, # 聚类树的列长度
treeheight_col = treeheight_col, # 聚类树行长度
scale = scale, # 矩阵有没有进行标准化
cluster_cols = cluster_cols, # 按列聚类
cluster_rows = cluster_rows, # 按行聚类
fontsize_row = fontsize_row, # 行字体大小
fontsize_col = fontsize_col, # ；列字体大小
cellwidth = cellwidth, cellheight = cellheight, # 设置图片大小
color = color, # 设置渐变色
border_color = border_color, # 每个小块间是否要用颜色分开
annotation_row = annotation_row,
annotation_col = annotation_col,
height = height,
width = width,
show_rownames = show_rownames,
show_colnames = show_colnames,
filename = NA
)
