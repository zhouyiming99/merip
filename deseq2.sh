library('DESeq2')


#注意是count数据，只能有基因ID和野生型突变体重复的一个文件，要提出来
base <- read.table("/data2/users/zhouyiming/pacbio_liu/correct/rna_seq/merged_bams/C12_NC.txt",row.names = 1 ,header=T,sep = '\t')  # 若遇到报错“多字节字符串错”，加上参数 fileEncoding="UTF-16LE"即可。

condition <- factor(c(rep("C",2), rep("a",2)))
#第一个rep：control数据，2个重复；第二个rep：处理（突变）数据，两个重复  
coldata <- data.frame(row.names = colnames(base), condition)
coldata    #显示coldata值,看看分组信息与真实数据是否一致，很关键！！！！

dds <- DESeqDataSetFromMatrix(countData=base, colData=coldata, design=~condition)
keep <- rowSums(counts(dds) > 0) >= 1  #过滤低表达基因，至少有1个样品都满足0个以上的reads数  
dds <- dds[keep, ]

dds1 <- DESeq(dds)    # 将数据标准化，必要步骤！！！
resultsNames(dds1)    # 查看结果的名称。
dds1$condition        #默认后者的处理组比前面的对照组。
res <- results(dds1)  # 必要步骤！！！
summary(res)          #看一下结果的概要信息，p值默认小于0.1。
table(res$padj < 0.05)        #padj 即矫正后的P值。看看有多少差异基因满足所设的P值要求。TRUE的数值为满足要求的基因个数。

res <- res[order(res$padj),]  #按照padj 进行升序排列
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds1, normalized=TRUE)),by="row.names",sort=FALSE)
write.table(resdata, file = "C4_diff.txt", sep = "\t", quote = FALSE)
