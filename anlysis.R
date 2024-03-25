library(ggplot2)
library(openxlsx)
library(patchwork)
library(ggsignif)
library(dplyr)
library(ggbreak)
library(stringr)
library(UpSetR)
library(reshape2)

########################################################################
#SNPs_per
########################################################################

a = read.table("vcf/merge/merged.vcf", sep = "\t", header = F)
all = a
strain_name=list("strain" = c(rep("AB", 1), 
                              rep("TU", 1),
                              rep("ABS", 3),
                              rep("India", 3),
                              rep("M-AB", 3),
                              rep("IM", 3)))

hetero = NULL
homo = NULL
for (i in 1:14) {
  #hetero
  b = a[grep("0/1",a[,i+9]),]
  hetero_count = c(hetero_count,nrow(b))
  c = (nrow(b)/1345101831)*100
  d = c(as.data.frame(strain_name$strain[i]),c)
  e = as.data.frame(d)
  colnames(e)[1] = "name"
  colnames(e)[2] = "hetero_genome"
  hetero = rbind(hetero,e)
  #homo
  b = a[grep("1/1",a[,i+9]),]
  homo_count = c(homo_count,nrow(b))
  c = (nrow(b)/1345101831)*100
  d = c(as.data.frame(strain_name$strain[i]),c)
  e = as.data.frame(d)
  colnames(e)[1] = "name"
  colnames(e)[2] = "homo_genome"
  homo = rbind(homo,e)
  rm(b)
  rm(c)
  rm(d)
  rm(e)
  cat(i)
}
strain = unique(strain_name$strain)

#hetero
hetero$name = factor(hetero$name,levels = strain)

hetero_mean  = group_by(hetero, name) %>% 
  summarise_all(funs(mean))

sd  = group_by(hetero, name) %>% 
  summarise_all(funs(sd))

se = sd[,2]/sqrt(3)

hetero_mean = cbind (hetero_mean, sd$hetero_genome, se)
hetero_mean = as.data.frame(hetero_mean)
colnames(hetero_mean) = c("name", "mean", "sd", "se")

hetero_mean$name = factor(hetero_mean$name,levels = strain)

hetero_g = ggplot(hetero_mean, aes(x = name, y = mean))+
  geom_bar(stat = "identity", fill ="white", color = "black")+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),width = .1)+
  ggtitle("Hetero%SNPs")+
  xlab("Labstoks")+
  ylab("Hetero%SNPs")+
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  theme_bw() +
  theme(title = element_text(size = 13)) +
  theme(axis.text.x = element_text(size = 15)) +
  theme(axis.text.y = element_text(size = 15)) +
  theme(axis.title.x = element_text(size = 13)) +
  theme(axis.title.y = element_text(size = 13))

hetero_g = hetero_g + 
  geom_jitter(data = hetero, aes(x = name, y = hetero_genome), size = 1) +
  theme_bw()+
  theme(legend.position = "none")

pdf("hetero_per.pdf", height = 5, width = 8)
hetero_g
dev.off()

#homo
homo$name = factor(homo$name,levels = strain)

homo_mean  = group_by(homo, name) %>% 
  summarise_all(funs(mean))

sd  = group_by(homo, name) %>% 
  summarise_all(funs(sd))

se = sd[,2]/sqrt(3)

homo_mean = cbind (homo_mean, sd$homo_genome, se)
homo_mean = as.data.frame(homo_mean)
colnames(homo_mean) = c("name", "mean", "sd", "se")

homo_mean$name = factor(homo_mean$name,levels = strain)

homo_g = ggplot(homo_mean, aes(x = name, y = mean))+
  geom_bar(stat = "identity", fill ="white", color = "black")+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),width = .1)+
  ggtitle("Homo%SNPs")+
  xlab("Labstoks")+
  ylab("homo%SNPs")+
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  theme_bw() +
  theme(title = element_text(size = 13)) +
  theme(axis.text.x = element_text(size = 15)) +
  theme(axis.text.y = element_text(size = 15)) +
  theme(axis.title.x = element_text(size = 13)) +
  theme(axis.title.y = element_text(size = 13))

homo_g = homo_g + 
  geom_jitter(data = homo, aes(x = name, y = homo_genome), size = 1) +
  theme_bw()+
  theme(legend.position = "none")

pdf("homo_per.pdf", height = 5, width = 8)
homo_g
dev.off()

########################################################################
#Denisity of Heterozygous SNPs
########################################################################

#all = read.table("vcf/merge/merged.vcf", sep = "\t", header = F)
base = read.table("zebrafish_genome_size.txt")
all = all[,c(1:2,seq(10,23, by=1))]
name = c("chr", "pos","AB","TU","AB*_1", "*AB_2", "*AB_3",
         "India_1", "India_2", "India_3",
         "M-AB_1", "M-AB_2", "M-AB_3",
         "IM_1", "IM_2", "IM_3")
colnames(all) = name
all = a

dt_list = list(c("ABM_all_1","ABM_all_2","ABM_all_3","IM_all_1","IM_all_2","IM_all_3"))
i=2
dt =NULL
for (i in 1:6) {
  b = 10+i
  c = all[, c(1:2, b)]
  dt[[1]][i] = list(c)
}
names(dt[[1]]) = c("ABM_all_1","ABM_all_2","ABM_all_3","IM_all_1","IM_all_2","IM_all_3")

dt_res = NULL
for (i in 1:6) {
  p = as.data.frame(dt[[1]][i])
  q = sapply(unique(p[,3]), function(x){sum(p[,3]==x)})
  max = dim(p)[1]
  a = t(apply(p[1:max,], 1, FUN = function(x){
    x = c(as.character(x[3]), c(names(q)))
    x = table(x)
    x = x-1
    return(x)
  }))
  a = as_data_frame(a)
  new = cbind(all[1:2],a$`0/1`)
  dt_res[[1]][i] = list(new)
}

head(dt_res[[1]][1])

t=1
chr=1
window = as.integer(10^5)
step = as.integer(10^5)
hetro_density = NULL
for (t in 1:6) {
  density = NULL
  pos.dt = as.data.frame(dt_res[[1]][t])
  for(chr in 1:25){
    max = chr.len[chr]
    pos.tmp = pos.dt[pos.dt$chr==chr,1:3]
    for(pos.end in seq(window, max,step)){
      cat(t,chr,":",pos.end/10^5,"\n")
      pos.start = pos.end - window +1
      a = pos.tmp$a..0.1.[pos.start<=pos.tmp$pos & pos.tmp$pos<pos.end]
      a = as.data.frame(a)
      a = sum(a)
      hetero = a
      tmp2 = data.frame("chromosome" = chr, "pos.start" = pos.start, "pos.end" = pos.end,"hetero" = hetero)
      density = rbind(density,tmp2)
    }
  }
  hetro_density[[1]][t] = list(density)
}
t=1
hetro_density_new = NULL
for (t in 1:6) {
  density = cbind(as.data.frame(hetro_density[[1]][t]),c(1:nrow(as.data.frame(hetro_density[[1]][t]))))
  name = c("chromosome", "pos.start", "pos.end","hetero","x")
  colnames(density) = name
  hetro_density_new[[1]][t] = list(density)
}

gl = list("gl_1", "gl_2", "gl_3", "gl_4", "gl_5", "gl_6")
original <- colorRampPalette(c("#23b16a", 
                               "#2bd5d5", 
                               "#2b80d5", 
                               "#2b2bd5", 
                               "#802bd5", 
                               "#d52bd5", 
                               "#d52b80", 
                               "#d52b2b",
                               "#d5802b", 
                               "#d5d52b", 
                               "#80d52b", 
                               "#2bd52b"))

for (t in 1:6) {
  hetro_density_in = as.data.frame(hetro_density_new[[1]][t])
  border = c(0,cumsum(table(hetro_density_in$chromosome))+0.5)
  gl[[t]]= ggplot(hetro_density_in, aes(x = x, y = hetero, color = as.factor(chromosome)))+
    scale_color_manual(palette = original)+
    geom_vline(xintercept = border, colour="grey", size =0.2)+
    ylab("Hetero_SNPs")+
    ylim(c(0,1300))+
    xlab("Chr")+
    geom_point(size = 0.01)+
    guides(colour = "none")+
    theme_bw()
}

p = gl[[1]] / gl[[2]] / gl[[3]] / gl[[4]] / gl[[5]] / gl[[6]]

pdf("hetero_SNPs_all.pdf",width = 10, height = 8)
p
dev.off()


sum(hetro_density_in$hetero)
sum(density$hetero)
sum(a$a..0.1..2)
a = as.data.frame(dt_res)

hetro_density_new
b = as.data.frame(hetro_density[[1]][1])

t=1
i=1
all = sum(chr.len)
new = NULL
func <- function(x) { rep(x, 1) }
new_list = NULL
for (t in 2:6) {
  hetro_density_in = as.data.frame(hetro_density_new[[1]][t])
  hetero_all = sum(hetro_density_in$hetero)
  pvalue = NULL
  pvalue <- matrix(0,length(hetro_density_in$hetero))
  for (i in 1:length(hetro_density_in$hetero)) {
    a = hetro_density_in[i,4]
    b = 10^5
    c = hetero_all
    d = all
    Sig_hetero = matrix(c(c, a, d, b), nrow = 2, dimnames = list(SNPs = c("all", "window"), Type = c("heteo", "others")))
    testResult = fisher.test(Sig_hetero)
    pvalue[i,] = func(testResult$p.value)
    cat(t,i,"\n")
  }
  t(pvalue)
  new = cbind(hetro_density_in,pvalue)
  new_list[[1]][t] = list(new) 
}

t(pvalue)
new = cbind(hetro_density_in,pvalue)


pvalue <- matrix(0,length(hetro_density_in$hetero))

cat(sprintf("%s/%s","end_",t))

class(new)

t=1
end_list_1 = NULL
end_list = NULL
min_list = NULL
new_new_pvalue = NULL
for (t in 1:6) {
  a = as.data.frame(new_list[[1]][t])
  new_new_pvalue = p.adjust(a$pvalue, "bonferroni")
  new_new_pvalue = cbind(a[,1:5],new_new_pvalue)
  end_list_1[[t]] = list(new_new_pvalue)
  new_new_pvalue = new_new_pvalue[new_new_pvalue$new_new_pvalue < 0.05,]
  new_new_pvalue = new_new_pvalue[new_new_pvalue$hetero != 0,]
  end_list[[t]] = list(new_new_pvalue)
  min_list[[t]] = list(min(new_new_pvalue$hetero))
} 

t=1
gll = NULL
for (t in 1:6) {
  gll[[t]] = gl[[t]] + geom_hline(aes(yintercept = min_list[[t]][[1]]),color = "red")
}

a = gll[[1]] / gll[[2]] / gll[[3]] / gll[[4]] / gll[[5]] / gll[[6]] 

pdf("hetero_SNPs_density.pdf",width = 10)
a
dev.off()

a = as.data.frame(end_list_1[[1]])
a = a[,c(5,c(1:4),6)]
for (t in 2:6) {
  b = as.data.frame(end_list_1[[t]])
  a = cbind(a,b[,c(4,6)])
}

write.csv(a,"Hetero_density.csv")

##chr

gl_chr = NULL
gl_chr[[t]][[i]] = list(sprintf("%s%s", "gl_chr_",i)) 

a_all = NULL
for (i in 1:25) {
  a = sprintf("%s_%s", "Chr", i)
  a_all = c(a_all,a)
}
a_all = list(a_all)  
name_list = list("M-AB_1","M-AB_2","M-AB_3", "IM_1", "IM_2", "IM_3")
t=2
i=1
gl_chr = NULL
for (t in 1:6) {
  hetro_density_in = as.data.frame(hetro_density_new[[1]][t])
  for(i in 1:25) {
    hetro_density_in_chr = as.data.frame(hetro_density_in[hetro_density_in$chromosome==i,])
    a = seq(1,length(hetro_density_in_chr$chromosome))
    hetro_density_in_chr = cbind(hetro_density_in_chr,a)
    gl_chr[[1]][[(t-1)*25+i]]= ggplot(hetro_density_in_chr, aes(x = a, y = hetero))+
      scale_color_manual(palette = original)+
      ylab(name_list[[t]])+
      ylim(c(0,1300))+
      xlab(sprintf(""))+
      geom_point(size = 0.01)+
      guides(colour = "none")+
      theme(axis.title.x=element_blank())+
      theme_bw()
    gl_chr[[1]][[(t-1)*25+i]] = gl_chr[[1]][[(t-1)*25+i]] + geom_hline(aes(yintercept = min_list[[t]][[1]]),color = "red")
    if (t==1) {
      gl_chr[[1]][[(t-1)*25+i]] = gl_chr[[1]][[(t-1)*25+i]] + ggtitle(a_all[[1]][i])
    }
  }
}

gl_chr_cmp = NULL
for(i in 1:25) {
  gl_chr_cmp[[i]] = gl_chr[[1]][[i]]/gl_chr[[1]][[i+25*1]]/gl_chr[[1]][[i+25*2]]/gl_chr[[1]][[i+25*3]]/gl_chr[[1]][[i+25*4]]/gl_chr[[1]][[i+25*5]]
}

pdf("chr_all_denity.pdf", width = 8, height = 10)
for (i in 1:25) {
  print(gl_chr_cmp[[i]])
}
dev.off()

########################################################################
# Anotation
########################################################################

ano_high_path = "vcf/"
ano_high_list = list(list.files(ano_high_path, recursive=, pattern="ano_all_HIGH.vcf.gz"))

vcf.ano2genename = function(snpeff){
  tmp = strsplit(snpeff$V8,"\\|")
  tmp = unlist(lapply(tmp, FUN = function(x){
    return(x[4])
  }))
  return(unique(tmp))
}
vcf.ano2geneID = function(snpeff){
  tmp = strsplit(snpeff$V8,"\\|")
  tmp = unlist(lapply(tmp, FUN = function(x){
    return(x[5])
  }))
  return(unique(tmp))
}
vcf.transcriptID = function(snpeff){
  tmp = strsplit(snpeff$V8,"\\|")
  tmp = unlist(lapply(tmp, FUN = function(x){
    return(x[7])
  }))
  return(unique(tmp))
}
gtf = read.table("../snpEff_v4_3t_core/snpEff/data/GRCz11.101/genes.gtf",sep = "\t", header = F)
gene_all = NULL
gene = gtf[gtf$V3=="exon",]
for (i in 1:length(gene$V1)) {
  gene_4 = NULL
  gene_1 = t(as.data.frame(unlist(strsplit(gene$V9[i], "; "))))
  gene_2 = t(as.data.frame(unlist(strsplit(gene_1[,6]," "))))
  gene_2 = gene_2[,2]
  gene_3 = t(as.data.frame(unlist(strsplit(gene_1[,3]," "))))
  gene_3 = gene_3[,2]
  gene_4 = cbind(gene_2, gene_3) 
  gene_all = rbind(gene_all,gene_4)
  cat(i,"\n")
}

gene_all = data.frame(gene_all)

gene_all_new = gene_all %>% distinct(gene_2,gene_3, .keep_all = T)

a = table(unlist(gene_all_new[,1]))
a = as.data.frame(a)
name = c("Gene_name", "Splicing_variant")
colnames(a) = name
Splicing_variant_name = a

#ABM
ABM_hetero_high = NULL
x = c(9,10,11)
for (i in x) {
  a = read.table(sprintf("vcf/%s", ano_high_list[[1]][i]),sep = "\t", header = F)
  b = a[grep("0/1",a$V10),]
  b = list(b)
  names(b) = sprintf("ABM_hetero_high_%s", i)
  ABM_hetero_high = c(ABM_hetero_high, b)
}

ABM_homo_high = NULL
x = c(9,10,11)
for (i in x) {
  a = read.table(sprintf("vcf/%s", ano_high_list[[1]][i]),sep = "\t", header = F)
  b = a[grep("1/1",a$V10),]
  b = list(b)
  names(b) = sprintf("ABM_homo_high_%s", i)
  ABM_homo_high = c(ABM_homo_high, b)
}



ABM_upset = list("ABM._high_homo_1" =  vcf.ano2genename(ABM_homo_high[[1]]), 
                 "ABM._high_homo_2" =  vcf.ano2genename(ABM_homo_high[[2]]),
                 "ABM._high_homo_3" =  vcf.ano2genename(ABM_homo_high[[3]]),
                 "ABM._high_hetero_1" =  vcf.ano2genename(ABM_hetero_high[[1]]), 
                 "ABM._high_hetero_2" =  vcf.ano2genename(ABM_hetero_high[[2]]),
                 "ABM._high_hetero_3" =  vcf.ano2genename(ABM_hetero_high[[3]]))

upset(fromList(ABM_upset),nsets = 8, nintersects = 8,order.by = "freq")

ABM_homo_upset = list("ABM._high_homo_1" =  vcf.ano2genename(ABM_homo_high[[1]]), 
                      "ABM._high_homo_2" =  vcf.ano2genename(ABM_homo_high[[2]]),
                      "ABM._high_homo_3" =  vcf.ano2genename(ABM_homo_high[[3]]))

i = upset(fromList(ABM_homo_upset),nsets = 8, 
          nintersects = 8,
          order.by = "freq",
          mainbar.y.max = 800, 
          sets = rev(c("ABM._high_homo_1","ABM._high_homo_2","ABM._high_homo_3")),
          keep.order = TRUE)

pdf("ABM_high_homoupset.pdf", height = 4 ,width = 6)
i
dev.off()

ABM_hetero_upset = list("ABM._high_hetero_1" =  vcf.ano2genename(ABM_hetero_high[[1]]), 
                        "ABM._high_hetero_2" =  vcf.ano2genename(ABM_hetero_high[[2]]),
                        "ABM._high_hetero_3" =  vcf.ano2genename(ABM_hetero_high[[3]]))

i = upset(fromList(ABM_hetero_upset),nsets = 100, 
          nintersects = 100,order.by = "freq",mainbar.y.max = 100,
          sets = rev(c("ABM._high_hetero_1","ABM._high_hetero_2","ABM._high_hetero_3")), 
          keep.order = TRUE, 
          boxplot.summary = NULL)

pdf("ABM_high_heteroupset.pdf", height = 3  ,width = 6)
i
dev.off()

#IM
IM_hetero_high = NULL
x = c(12,13,14)
for (i in x) {
  a = read.table(sprintf("vcf/%s", ano_high_list[[1]][i]),sep = "\t", header = F)
  b = a[grep("0/1",a$V10),]
  b = list(b)
  names(b) = sprintf("IM_hetero_high_%s", i)
  IM_hetero_high = c(IM_hetero_high, b)
}

IM_homo_high = NULL
x = c(12,13,14)
for (i in x) {
  a = read.table(sprintf("vcf/%s", ano_high_list[[1]][i]),sep = "\t", header = F)
  b = a[grep("1/1",a$V10),]
  b = list(b)
  names(b) = sprintf("IM_homo_high_%s", i)
  IM_homo_high = c(IM_homo_high, b)
}

IM_upset = list("IM._high_homo_1" =  vcf.ano2genename(IM_homo_high[[1]]), 
                "IM._high_homo_2" =  vcf.ano2genename(IM_homo_high[[2]]),
                "IM._high_homo_3" =  vcf.ano2genename(IM_homo_high[[3]]),
                "IM._high_hetero_1" =  vcf.ano2genename(IM_hetero_high[[1]]), 
                "IM._high_hetero_2" =  vcf.ano2genename(IM_hetero_high[[2]]),
                "IM._high_hetero_3" =  vcf.ano2genename(IM_hetero_high[[3]]))

upset(fromList(IM_upset),nsets = 8, nintersects = 100,order.by = "freq")

ABM_IM_upset = list("IM._high_homo_1" =  vcf.ano2genename(IM_homo_high[[1]]), 
                    "IM._high_homo_2" =  vcf.ano2genename(IM_homo_high[[2]]),
                    "IM._high_homo_3" =  vcf.ano2genename(IM_homo_high[[3]]),
                    "IM._high_hetero_1" =  vcf.ano2genename(IM_hetero_high[[1]]), 
                    "IM._high_hetero_2" =  vcf.ano2genename(IM_hetero_high[[2]]),
                    "IM._high_hetero_3" =  vcf.ano2genename(IM_hetero_high[[3]]),
                    "ABM._high_homo_1" =  vcf.ano2genename(ABM_homo_high[[1]]), 
                    "ABM._high_homo_2" =  vcf.ano2genename(ABM_homo_high[[2]]),
                    "ABM._high_homo_3" =  vcf.ano2genename(ABM_homo_high[[3]]),
                    "ABM._high_hetero_1" =  vcf.ano2genename(ABM_hetero_high[[1]]), 
                    "ABM._high_hetero_2" =  vcf.ano2genename(ABM_hetero_high[[2]]),
                    "ABM._high_hetero_3" =  vcf.ano2genename(ABM_hetero_high[[3]]))

upset(fromList(ABM_IM_upset),nsets = 12, nintersects = 28,order.by = "freq")

ABM_IM_hetero_upset = list("ABM._high_hetero_1" =  vcf.ano2genename(ABM_hetero_high[[1]]), 
                           "ABM._high_hetero_2" =  vcf.ano2genename(ABM_hetero_high[[2]]),
                           "ABM._high_hetero_3" =  vcf.ano2genename(ABM_hetero_high[[3]]),
                           "IM._high_hetero_1" =  vcf.ano2genename(IM_hetero_high[[1]]), 
                           "IM._high_hetero_2" =  vcf.ano2genename(IM_hetero_high[[2]]),
                           "IM._high_hetero_3" =  vcf.ano2genename(IM_hetero_high[[3]]))

a = upset(fromList(ABM_IM_hetero_upset),nsets = 12, nintersects = 10,order.by = "freq", sets = rev(c("ABM._high_hetero_1","ABM._high_hetero_2","ABM._high_hetero_3","IM._high_hetero_1","IM._high_hetero_2","IM._high_hetero_3")), keep.order = TRUE)

ABM_IM_homo_upset = list("ABM._high_homo_1" =  vcf.ano2genename(ABM_homo_high[[1]]), 
                         "ABM._high_homo_2" =  vcf.ano2genename(ABM_homo_high[[2]]),
                         "ABM._high_homo_3" =  vcf.ano2genename(ABM_homo_high[[3]]),
                         "IM._high_homo_1" =  vcf.ano2genename(IM_homo_high[[1]]), 
                         "IM._high_homo_2" =  vcf.ano2genename(IM_homo_high[[2]]),
                         "IM._high_homo_3" =  vcf.ano2genename(IM_homo_high[[3]]))

i = upset(fromList(ABM_IM_homo_upset),nsets = 12, nintersects = 10,order.by = "freq", sets = rev(c("ABM._high_homo_1","ABM._high_homo_2","ABM._high_homo_3","IM._high_homo_1","IM._high_homo_2","IM._high_homo_3")), keep.order = TRUE)

pdf("ABM_IM_homo.pdf", height = 5 ,width = 6)
i
dev.off()

i = set_intersection(vcf.ano2genename(ABM_hetero_high[[1]]),
                     vcf.ano2genename(ABM_hetero_high[[2]]), 
                     vcf.ano2genename(ABM_hetero_high[[3]]))


ABM_set_gene = as.matrix(i)

i = set_intersection(vcf.ano2genename(IM_hetero_high[[1]]),
                     vcf.ano2genename(IM_hetero_high[[2]]), 
                     vcf.ano2genename(IM_hetero_high[[3]]))
IM_set_gene = as.matrix(i)


write.csv(ABM_set_gene,file = "ABM_set_gene.csv")
write.csv(IM_set_gene,file = "IM_set_gene.csv")

i = set_intersection(vcf.ano2genename(ABM_hetero_high[[1]]), 
                     vcf.ano2genename(ABM_hetero_high[[2]]),
                     vcf.ano2genename(ABM_hetero_high[[3]]),
                     vcf.ano2genename(IM_hetero_high[[1]]), 
                     vcf.ano2genename(IM_hetero_high[[2]]),
                     vcf.ano2genename(IM_hetero_high[[3]]))

pdf("ABM_IM_hetero.pdf", height = 5 ,width = 6)
a
dev.off()

t = 1
i = 1
d = NULL
ABM_homo_high_list = NULL
ABM_homo_high_list_name = NULL
for (t in 1:3) {
  d = NULL
  for (i in 1:length(ABM_homo_high[[t]]$V1)) {
    a = cbind(vcf.ano2genename(ABM_homo_high[[t]][i,]),vcf.transcriptID(ABM_homo_high[[t]][i,]))
    b = Splicing_variant_name[Splicing_variant_name[,1]==a[,1],2]
    c = length(b)
    if (c > 0) {
      if (b==1) {
        a = cbind(a,b)
        d = rbind(d,a)
      }
    }
  }
  ABM_homo_high_list[t] = list(unlist(d)) 
  ABM_homo_high_list_name[t] = list(unique(unlist(d[,1])))
}

IM_homo_high_list = NULL
IM_homo_high_list_name = NULL
for (t in 1:3) {
  d = NULL
  for (i in 1:length(IM_homo_high[[t]]$V1)) {
    a = cbind(vcf.ano2genename(IM_homo_high[[t]][i,]),vcf.transcriptID(IM_homo_high[[t]][i,]))
    b = Splicing_variant_name[Splicing_variant_name[,1]==a[,1],2]
    c = length(b)
    if (c > 0) {
      if (b==1) {
        a = cbind(a,b)
        d = rbind(d,a)
      }
    }
  }
  IM_homo_high_list[t] = list(unlist(d)) 
  IM_homo_high_list_name[t] = list(unique(unlist(d[,1])))
}

names(ABM_homo_high_list_name) = c("ABM_1","ABM_2","ABM_3")
class(b)

names(IM_homo_high_list_name) = c("IM_1","IM_2","IM_3")
class(b)

ABM_IM_list = c(ABM_homo_high_list_name,IM_homo_high_list_name)
upset(fromList(ABM_IM_list),nsets = 6 ,order.by = "freq", keep.order = TRUE)

a = as.character(ABM_IM_list$IM_1)
b = as.character(ABM_IM_list$IM_2)
c = as.character(ABM_IM_list$IM_3)
a[is.element(a,b)]
d = as.character(ABM_IM_list$ABM_1)
e = as.character(ABM_IM_list$ABM_2)
f = as.character(ABM_IM_list$ABM_3)
g = unlist(set_intersection(ABM_IM_list$IM_1,
                            ABM_IM_list$IM_2,
                            ABM_IM_list$IM_3))
g = g[!is.element(g,d)]
g = g[!is.element(g,e)]
g = g[!is.element(g,f)]

g = g[grep("si",g,invert = T)]
g = g[grep("zgc",g,invert = T)]
gnl = g[grep("[A-Z]",g,invert = T)]

a = IM_homo_high$IM_homo_high_37
out = NULL
for (gn in gnl) {
  b = a[grep(gn,a$V8),]
  b = b[1,]
  c = cbind(b$V1,b$V2)
  b = b$V8
  b = unlist(lapply(strsplit(b,"\\|"),FUN = function(x){return(x[c(2,7,4,10,13,14)])}))
  b = matrix(b,nrow = 1)
  b = cbind(c,b)
  out = rbind(out,b)
}


out = NULL
ABM_homo_high_list = NULL
for (z in 1:3) {
  a = ABM_homo_high[[z]]
  for (gn in 1:length(a$V1)) {
    b = a[gn,]
    p = cbind(b$V1,b$V2)
    b = b$V8
    q = unlist(strsplit(b,"\\|"))
    i = unlist(grep("HIGH",q))
    g = NULL
    for (n in c(i)) {
      c = b
      e = unlist(lapply(strsplit(c,"\\|"),FUN = function(x){return(x[c(n-1, n+1, n+4, n+7, n+8, n+9, n+10)])}))
      e = matrix(e,nrow = 1)
      d = Splicing_variant_name[Splicing_variant_name[,1]==e[2],2]
      if (length(d) == 0) {
        d =0
      }
      f = NULL
      f = cbind(p,e,d)
      g = rbind(g,f)
    }
    u = as.data.frame(table(g[,4]))
    for (t in 1:length(g[,1])) {
      for (i in 1:length(u[,1])) {
        x = g[g[t,4]==u[i,1],10]
        if (length(x) == 0) {
          x =0
        }
        if (x == u[,2]) {
          out = rbind(out,g[t,])
        }
      } 
    }
  }
  out = as.data.frame(out)
  ABM_homo_high_list[[z]] = list(out)
}

ABM_homo_high_unique_list = NULL
for (i in 1:3) {
  g = unique(ABM_homo_high_list[[i]][[1]][[4]])
  g = g[grep("si",g,invert = T)]
  g = g[grep("zgc",g,invert = T)]
  gnl = g[grep("[A-Z]",g,invert = T)]
  ABM_homo_high_unique_list[[i]] = list(gnl) 
}

IM_homo_high_list = NULL
for (z in 1:3) {
  out = NULL
  a = IM_homo_high[[z]]
  for (gn in 1:length(a$V1)) {
    b = a[gn,]
    p = cbind(b$V1,b$V2)
    b = b$V8
    q = unlist(strsplit(b,"\\|"))
    i = unlist(grep("HIGH",q))
    g = NULL
    for (n in c(i)) {
      c = b
      e = unlist(lapply(strsplit(c,"\\|"),FUN = function(x){return(x[c(n-1, n+1, n+4, n+7, n+8, n+9, n+10)])}))
      e = matrix(e,nrow = 1)
      d = Splicing_variant_name[Splicing_variant_name[,1]==e[2],2]
      if (length(d) == 0) {
        d =0
      }
      f = NULL
      f = cbind(p,e,d)
      g = rbind(g,f)
    }
    u = as.data.frame(table(g[,4]))
    for (t in 1:length(g[,1])) {
      for (i in 1:length(u[,1])) {
        x = g[g[t,4]==u[i,1],10]
        if (length(x) == 0) {
          x =0
        }
        if (x == u[,2]) {
          out = rbind(out,g[t,])
        }
      } 
    }
  }
  out = as.data.frame(out)
  IM_homo_high_list[[z]] = list(out)
}

IM_homo_high_unique_list = NULL
for (i in 1:3) {
  g = unique(IM_homo_high_list[[i]][[1]][[4]])
  g = g[grep("si",g,invert = T)]
  g = g[grep("zgc",g,invert = T)]
  gnl = g[grep("[A-Z]",g,invert = T)]
  IM_homo_high_unique_list[[i]] = list(gnl) 
}

a = as.matrix(set_intersection(IM_homo_high_unique_list[[1]][[1]], IM_homo_high_unique_list[[2]][[1]], IM_homo_high_unique_list[[3]][[1]]))

a[grep("nbas",a)]

ABM_IM_homo_upset = NULL
ABM_IM_homo_upset = list("ABM._high_homo_1" = unlist(ABM_homo_high_unique_list[[1]]),
                         "ABM._high_homo_2" = unlist(ABM_homo_high_unique_list[[2]]),
                         "ABM._high_homo_3" = unlist(ABM_homo_high_unique_list[[3]]),
                         "IM._high_homo_1" =  unlist(IM_homo_high_unique_list[[1]]),
                         "IM._high_homo_2" =  unlist(IM_homo_high_unique_list[[2]]),
                         "IM._high_homo_3" = unlist(IM_homo_high_unique_list[[3]]))

a = upset(fromList(ABM_IM_homo_upset),nsets = 12,order.by = "freq", keep.order = TRUE)

pdf("Homo_high.pdf",width = 10)
a
dev.off()


##IM_Homo

d = as.character(ABM_IM_homo_upset$ABM._high_homo_1)
e = as.character(ABM_IM_homo_upset$ABM._high_homo_2)
f = as.character(ABM_IM_homo_upset$ABM._high_homo_3)
g = unlist(set_intersection(ABM_IM_homo_upset$IM._high_homo_1,
                            ABM_IM_homo_upset$IM._high_homo_2,
                            ABM_IM_homo_upset$IM._high_homo_3))
g = g[!is.element(g,d)]
g = g[!is.element(g,e)]
gnl = g[!is.element(g,f)]
z=1

c = NULL
for (i in 1:length(gnl)) {
  a = as.data.frame(IM_homo_high_list[[z]])
  b = a[grep(gnl[i],a$V4),]
  c = rbind(c,b)
}

c = c %>% 
  distinct(V1, V2, V3, V4, V5, V6, V8, V9, d ,.keep_all=TRUE)
length(c$V1)
unique(c$V4)

c = c[,c(1:7, 9, 10)]
length()

write.csv(c,"Homo_IM.csv")

##ABM_Homo

d = as.character(ABM_IM_homo_upset$IM._high_homo_1)
e = as.character(ABM_IM_homo_upset$IM._high_homo_2)
f = as.character(ABM_IM_homo_upset$IM._high_homo_3)
g = unlist(set_intersection(ABM_IM_homo_upset$ABM._high_homo_1,
                            ABM_IM_homo_upset$ABM._high_homo_2,
                            ABM_IM_homo_upset$ABM._high_homo_3))
g = g[!is.element(g,d)]
g = g[!is.element(g,e)]
gnl = g[!is.element(g,f)]
z=1

c = NULL
for (i in 1:length(gnl)) {
  a = as.data.frame(ABM_homo_high_list[[z]])
  b = a[grep(gnl[i],a$V4),]
  c = rbind(c,b)
}

c = c %>% 
  distinct(V1, V2, V3, V4, V5, V6, V8, V9, d ,.keep_all=TRUE)
length(c$V1)
unique(c$V4)

c = c[,c(1:7, 9, 10)]
length()

write.csv(c,"Homo_ABM.csv")


#ABM_IM_Homo

g = unlist(set_intersection(ABM_IM_homo_upset$IM._high_homo_1,
                            ABM_IM_homo_upset$IM._high_homo_2,
                            ABM_IM_homo_upset$IM._high_homo_3,
                            ABM_IM_homo_upset$ABM._high_homo_1,
                            ABM_IM_homo_upset$ABM._high_homo_2,
                            ABM_IM_homo_upset$ABM._high_homo_3))

gnl = g
z=1

c = NULL
for (i in 1:length(gnl)) {
  a = as.data.frame(ABM_homo_high_list[[z]])
  b = a[grep(gnl[i],a$V4),]
  c = rbind(c,b)
}

c = c %>% 
  distinct(V1, V2, V3, V4, V5, V6, V8, V9, d ,.keep_all=TRUE)
length(c$V1)
unique(c$V4)

c = c[,c(1:7, 9, 10)]
length()

write.csv(c,"Homo_ABM_IM.csv")
