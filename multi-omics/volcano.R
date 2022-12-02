data = read.table()
FC = 1.5
PValue = 0.05
pointcolor = c("red","grey","blue")

data$sig[(-1*log10(data$PValue) < -1*log10(PValue)|data$PValue=="NA")|(log2(data$FC) < log2(FC))& log2(data$FC) > -log2(FC)] <- "Normal"
data$sig[-1*log10(data$PValue) >= -1*log10(PValue) & log2(data$FC) >= log2(FC)] <- "up-regulated"
data$sig[-1*log10(data$PValue) >= -1*log10(PValue) & log2(data$FC) <= -log2(FC)] <- "down-regulated"
  

# 火山图标记 #
library(ggrepel)
pic_repel <- geom_text_repel(aes(x = log2(FC),                        # geom_text_repel 标记函数
							   y = -1*log10(PValue),          
							   label = label),                       
						   max.overlaps = 10000,                    # 最大覆盖率，当点很多时，有些标记会被覆盖，调大该值则不被覆盖，反之
						   size = 3,                                # 字体大小
						   box.padding=unit(0.5,'lines'),           # 文本框周边填充
						   point.padding=unit(0.1, 'lines'),        # 点周围填充
						   segment.color='black',                   # 标记线条的颜色
						   show.legend=FALSE)

x = min(abs(log2(sort(data$FC[which(data$FC != Inf)],decreasing = T)[1])),
	  abs(log2(sort(data$FC[which(data$FC != Inf)],decreasing = F)[1])))


# 火山图绘制 #
library(ggplot2)
plotpic <- ggplot(data,aes(log2(FC),-1*log10(PValue))) +                                 # 加载数据，定义横纵坐标
geom_point(aes(color = sig)) +                                                         # 绘制散点图，分组依据是数据框的sig列
labs(title="Volcano Plot",                                                             # 定义标题，x轴，y轴名称
	 x="log[2](FC)", 
	 y="-log[10](PValue)") + 
theme(plot.title = element_text(hjust = 0.5)) +                                        #设置标题居中
scale_color_manual(values = pointcolor) +                                              # 自定义颜色，将values更改成你想要的三个颜色
geom_hline(yintercept=-log10(PValue),linetype=2) +                                     # 在图上添加虚线
geom_vline(xintercept=c(-log2(FC),log2(FC)),linetype=2) +                              # 在图上添加虚线
theme(plot.margin=unit(rep(1,4),'lines'), panel.spacing=unit(c(4,4,4,4), "cm")) + 
theme(panel.border = element_rect(fill=NA,color="black", size=0.5,linetype="solid")) + 
theme(panel.background = element_blank()) + theme(legend.key=element_blank()) +        # 移除图像背景和图例键 同theme_bw()
theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +            # 移除灰色背景和网格
ylim(c(0,-log10(sort(data$PValue[which(data$PValue != 0)], decreasing = F)[1])*1.2)) + # 设置y轴的取值范围
xlim(c(-floor(x*1.2),floor(x*1.2))) +                                                  # 设置x轴的取值范围
theme(legend.justification = c(.98,.98), legend.position = c(.98,.98)) +               # 设置图例嵌入图片的位置
theme(legend.background = element_rect(fill = 'white', colour = 'black')) +            # 设置图例边框
theme(legend.title = element_blank())
