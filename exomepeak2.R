library("exomePeak2")
d1_IP <-  "/data1/zhouyiming/liu_merip/BAM/C8MS1_Ip.uni.bam"
d1_input <-  "/data1/zhouyiming/liu_merip/BAM/C8MS1_Input.uni.bam"
d2_IP <-  "/data1/zhouyiming/liu_merip/BAM/C8MS2_Ip.uni.bam"
d2_input <-  "/data1/zhouyiming/liu_merip/BAM/C8MS2_Input.uni.bam"
# 野生型c样本
c1_IP <-  "/data1/zhouyiming/liu_merip/BAM/NCMS1_Ip.uni.bam"
c1_input <-  "/data1/zhouyiming/liu_merip/BAM/NCMS1_Input.uni.bam"
c2_IP <-  "/data1/zhouyiming/liu_merip/BAM/NCMS2_Ip.uni.bam"
c2_input <-  "/data1/zhouyiming/liu_merip/BAM/NCMS2_Input.uni.bam"
# 假设IP_BAM和INPUT_BAM变量是用于存储所有样本的IP和input BAM文件路径的向量。
IP_BAM <- c(c1_IP, c2_IP)
INPUT_BAM <- c(c1_input, c2_input)
# 假设TREATED_IP_BAM和TREATED_INPUT_BAM变量是用于存储处理过的IP和input  BAM文件路径的向量。
TREATED_IP_BAM <- c(d1_IP, d2_IP)
TREATED_INPUT_BAM <- c(d1_input, d2_input)
# 指定基因注释文件的路径和基因组版本
GENE_ANNO_GTF <-  "/data1/zhouyiming/liu_merip/Oryza_rufipogon.OR_W1943.60.gtf"
output_dir <- "/data1/zhouyiming/liu_merip"
# 调用exomePeak2函数
result <- exomePeak2(bam_ip = IP_BAM,
                     bam_input = INPUT_BAM,
                     bam_treated_input = TREATED_INPUT_BAM,
                     bam_treated_ip = TREATED_IP_BAM,
                     gff_dir = GENE_ANNO_GTF, 
                     save_dir = output_dir)
