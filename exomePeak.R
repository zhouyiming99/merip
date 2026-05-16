library("exomePeak")
gtf <- "/data1/zhouyiming/liu_merip/Oryza_rufipogon.OR_W1943.60.gtf"
f1 <- "/data1/zhouyiming/liu_merip/BAM/C8MS1_Ip.uni.bam"
f2 <- "/data1/zhouyiming/liu_merip/BAM/C8MS1_Input.uni.bam"
f3 <- "/data1/zhouyiming/liu_merip/BAM/C8MS2_Ip.uni.bam"
f4 <- "/data1/zhouyiming/liu_merip/BAM/C8MS2_Input.uni.bam"
f5 <- "/data1/zhouyiming/liu_merip/BAM/NCMS1_Ip.uni.bam"
f6 <- "/data1/zhouyiming/liu_merip/BAM/NCMS1_Input.uni.bam"
f7 <- "/data1/zhouyiming/liu_merip/BAM/NCMS2_Ip.uni.bam"
f8 <- "/data1/zhouyiming/liu_merip/BAM/NCMS2_Input.uni.bam"



result <- exomepeak(GENE_ANNO_GTF = gtf,
                    IP_BAM = c(f5, f7),
                    INPUT_BAM = c(f6, f8),
                    TREATED_IP_BAM = c(f1,f3),
                    TREATED_INPUT_BAM = c(f2,f4))
