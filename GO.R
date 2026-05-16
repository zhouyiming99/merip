
a <- fread("/data1/zhouyiming/qianqian/merip2/exomePeak_output/diff_peak.txt") 
b <- as.data.table(a)
c <- b[diff == "Up", unique(name)]
write.table(c, file = "name_up.txt", sep = "\t", quote = FALSE, row.names = FALSE)

/data1/liuqi/miniconda3/envs/R4.3.0/bin/R
library("data.table")
library(ggplot2)
# 读取数据
go_enrich_1 <- fread("/data1/zhouyiming/DT4/crtl-1-1/mod_diff/xpore/results/filtered.txt")


# 筛选出 term_type 为 "P" 的条目
p_terms <- go_enrich[go_enrich$Category == "GOTERM_BP_DIRECT", ]
f_terms <- go_enrich[go_enrich$Category == "GOTERM_MF_DIRECT", ]
c_terms <- go_enrich[go_enrich$Category == "GOTERM_CC_DIRECT", ]
# 根据 FDR 列进行升序排序
p_terms <- p_terms[order(p_terms$FDR, decreasing = FALSE), ]


# 提取排序后的前 10 个条目
p_terms <- head(p_terms, 10)




# 根据 fdr 值降序排列并提取前10个 F 类别的条目
f_terms <- head(go_enrich[go_enrich$term_type == "F", ], 10)
f_terms <- f_terms[order(f_terms$FDR, decreasing = FALSE), ]




# 根据 fdr 值降序排列并提取前10个 C 类别的条目
c_terms <- head(go_enrich[go_enrich$term_type == "C", ], 10)
c_terms <- c_terms[order(c_terms$FDR, decreasing = FALSE), ]




# 合并 P、F 和 C 类别的条目
all_terms <- rbind(p_terms, f_terms, c_terms)
write.table(all_terms, '/data1/zhouyiming/qianqian/merip2/exomePeak_output/GO/GO_FDR10_up.txt', sep = '\t', quote = FALSE, row.names = FALSE)、、


go_enrich_1 <- fread("/data1/zhouyiming/DT4/crtl-1-1/mod_diff/xpore/results/filtered_bp.txt")

# 使用reorder函数调整term的顺序
go_enrich_1$term <- reorder(go_enrich_1$Term, go_enrich_1$Count)


p1 <- ggplot(go_enrich_1,
       aes(x=Term,y=Count, fill=Category)) +
  geom_bar(stat="identity", width=0.8) +  
  scale_fill_manual(values = c("#6666FF", "#33CC33", "#FF6666") ) +  
  coord_flip() +  
  xlab("GO term") +  
  ylab("Gene Number") +  
  labs(title = "DEG AS GO Terms Enrich") +  
  theme_bw()
# 根据ONTOLOGY分类信息添加分组框，并按照Count列的大小排列y轴
p1 <- p1 + facet_grid(Category ~ ., scale = 'free_y', space = 'free_y')
ggsave("GO.png",plot=p1, width = 10, height = 2, dpi = 300)

# 根据ONTOLOGY分类信息添加分组框，并按照Count列的大小排列y轴
p1 <- p1 + facet_grid(Category ~ ., scale = 'free_y', space = 'free_y')
library(gridExtra)

# 合并 p3 和 p4 到一张图上
combined_plot_2 <- grid.arrange(p3, p4, ncol=2)

# 保存合并后的图形为 PNG 文件
ggsave("GO_fdr_combined.png", combined_plot_2, width = 20, height = 6, dpi = 300)

#纵向柱状图——-根据pvalue值绘制#

# 根据 fdr 值降序排列并提取前10个 P 类别的条目
p_terms <- head(go_enrich[go_enrich$term_type == "P", ], 10)
p_terms <- p_terms[order(p_terms$FDR, decreasing = TRUE), ]


# 根据 fdr 值降序排列并提取前10个 F 类别的条目
f_terms <- head(go_enrich[go_enrich$term_type == "F", ], 10)
f_terms <- f_terms[order(f_terms$FDR, decreasing = TRUE), ]


# 根据 fdr 值降序排列并提取前10个 C 类别的条目
c_terms <- head(go_enrich[go_enrich$term_type == "C", ], 10)
c_terms <- c_terms[order(c_terms$FDR, decreasing = TRUE), ]


# 合并 P、F 和 C 类别的条目
all_terms <- rbind(p_terms, f_terms, c_terms)

# 将横坐标的向量转换为因子 all_terms$queryitem <- as.factor(all_terms$queryitem)

#纵向柱状图——-根据pvalue值绘制#
p4 <- ggplot(go_enrich_1, aes(y = Term, x = Count, fill = PValue))+  #x、y轴定义；根据pvalue填充颜色
  geom_bar(stat = "identity",width=0.8)+ #柱状图宽度设置
  scale_fill_gradient(low = "red",high ="blue" )+
  labs(title = "Up GO Terms Enrich",  #设置标题、x轴和Y轴名称
       x = "Gene number",
       y = "GO Terms")+
  theme(axis.title.x = element_text(face = "bold",size = 16),
        axis.title.y = element_text(face = "bold",size = 16),
        legend.title = element_text(face = "bold",size = 16))+
  theme_bw()
p3
#根据ONTOLOGY分类信息添加分组框#
p4 <- p4 + facet_grid(term_type ~ ., scale = 'free_y', space = 'free_y')



ggsave("GO_fdr_Enrich_down.png", plot=p3, dpi=300)





2. df<-fread("./data/results/B5oHnvXQUHUuCMMU/pathway.down.txt")：使用fread函数读取指定路径下的文件，并将其存储为数据框df。
3. df <- df[!GO_acc %in% c("R-MMU-163200"),]：去除ID列包含特定数值的行。
4. df <- df[order(-term_type,-Count)]：按照Type和Count列的降序排列数据框df。
5. df$Description <- factor(df$Description, levels =rev(df$Description))：将Description列转换为因子型变量，并按照逆序重新排列。
6. ggplot(df, aes_string(x = "Description", y = "Count", fill = "p_adjust"))：创建一个基础图形对象，并设置x轴为Description列，y轴为Count列，填充颜色根据p_adjust列变化。
7. theme(text = element_text(size=12))：设置主题格式，包括文本大小等。
8. scale_fill_continuous(low="red", high="blue", name = 'p.adjust',guide=guide_colorbar(reverse=TRUE))：设置填充色调，从红色到蓝色，标注名称为p.adjust，颜色栏反转。
9. geom_bar(stat = "identity")：添加柱状图层，并使用Count列的值作为柱子高度。
10. coord_flip()：交换x轴和y轴的位置。
11. ggtitle("") + xlab(NULL) + ylab(NULL) + theme_bw()：设置标题为空白，x轴和y轴标签为空白，使用黑白主题。
12. ggsave("./data/results/B5oHnvXQUHUuCMMU/pathway.barplot.all.pdf",width=10,height=8)：将绘制好的图表保存为PDF格式，并指定保存路径和大小。
#气泡图#
P5 <- ggplot(go_enrich,
       aes(y=Description,x=Count))+
  geom_point(aes(size=Count,color=padj))+
  scale_color_gradient(low = "red",high ="blue")+
  labs(color=expression(padj,size="Count"),
       x="Gene Ratio",y="Description",title="GO Enrichment")+
  theme_bw()
P5
#根据ONTOLOGY分类信息添加分组框#
P5 + facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y')




