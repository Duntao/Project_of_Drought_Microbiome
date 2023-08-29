#The following core code are used for drought microbiome analysis.
#Pls don not hestitate to email me if your want to conduct deeply communication with us
#The loading R packages upon yourself system status and dadaset


#Figure1################################################################################
library(ggplot2)
library(ggpubr)
library(agricolae)
library(vegan)
library(ggsci)
library(gridExtra)
library(cowplot)
library(tidyverse)
####bacteria####
otu_bac <- read.delim("otu_taxon_bac.xls",row.names = 1,header=T,comment.char='',sep='\t',check.names = F)
otu_bac=otu_bac[,-ncol(otu_bac)]
otu_bac=t(otu_bac)
shannon_bac=diversity(otu_bac,"shannon")
simpson_bac=diversity(otu_bac,"simpson")
Chao1_bac <- estimateR(otu_bac)[2,]
Ace_bac <- estimateR(otu_bac)[4,]
####fungi####
otu_fun <- read.delim("otu_taxon_fun.xls",row.names = 1,header=T,comment.char='',sep='\t',check.names = F)
otu_fun=otu_fun[,-ncol(otu_fun)]
otu_fun=t(otu_fun)
shannon_fun=diversity(otu_fun,"shannon")
simpson_fun=diversity(otu_fun,"simpson")
Chao1_fun <- estimateR(otu_fun)[2,]
Ace_fun <- estimateR(otu_fun)[4,]
alpha_bac=data.frame(shannon_bac,simpson_bac,Chao1_bac,Ace_bac)
colnames(alpha_bac) <-c("shannon","simpson","Chao","Ace")
alpha_fun <- data.frame(shannon_fun,simpson_fun,Chao1_fun,Ace_fun)
colnames(alpha_fun) <-c("shannon","simpson","Chao","Ace")
write.table(alpha_bac,'alpha-summary_bac.tsv',sep = '\t',quote=F)
write.table(alpha_fun,'alpha-summary_fun.tsv',sep = '\t',quote=F)
group<-read.table('sample.csv',row.names = 1,
                header = T,sep=',',comment.char='',check.names=F)
alpha_bac$Original_name <- rownames(alpha_bac)
alpha_group_bac <- left_join(group,alpha_bac,by="Original_name")
alpha_fun$Original_name <- rownames(alpha_fun)
alpha_group_fun <- left_join(group,alpha_fun,by="Original_name")
alpha_group <- rbind(alpha_group_bac,alpha_group_fun)
alpha_group$Type <- c(rep("bacteria",40),rep("fungi",40))

df <- alpha_group
p_shannon_fungi <- ggplot(data=df[df$Type=="fungi",],aes(x=Group,y=shannon,fill=Group))+
stat_boxplot(geom = "errorbar",
             width=0.3,
             position = position_dodge(0.9))+
geom_boxplot(position = position_dodge(0.9),
outlier.colour = "red")+
  geom_jitter(position = position_jitter(0.17),size=2,alpha=1)+
  stat_compare_means(method = "kruskal.test",label.y = 4,label.x = 0.8,label = "p.format")+
  geom_text(data=df_LSD,mapping=aes(x=Group,y=shannon_fun+std_fun+0.15,label=label_2))+ 
  theme(text = element_text(size = 10,face = "bold"),
        axis.text.x = element_text(angle = 0),
        panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        legend.position = "none",
        plot.margin=unit(rep(0.5,4),'cm'))+
  labs(x="")+
  theme_bw()+
  scale_fill_npg(palette = c("nrc"), alpha = 1)
plot_grid(p_shannon_bacteria,p_shannon_fungi,p_Chao_bacteria,p_Chao_fungi,
          labels = c("A","B","C","D"),ncol = 2,nrow = 2)
ggsave("diversity.pdf",width = 10,height = 9)

library(tidyverse)
library(scales)
library(ggh4x)
library(patchwork)
library(magrittr)
library(cowplot)
computed_persent <- function(path) {
  data <- path %>%
    read.delim(check.names = F,sep=",",row.names = 1) %>% 
    as.data.frame()
  data2 <- data %>%
    dplyr::mutate(sum = rowSums(.), persent = sum / sum(sum) * 100, 
           sum = NULL,) %>%
    rbind(filter(., persent < 0.1) %>% colSums()) %>% 
    dplyr::mutate(Taxa = c(data %>% rownames(), "others"))
  filter(data2[1:(nrow(data2) - 1),], persent > 0.1) %>%
    rbind(data2[nrow(data2),]) %>%
    dplyr::select(ncol(.), 1:(ncol(.) - 2)) %>%
    set_rownames(seq_len(nrow(.))) %>%
    return()
}

aaa <- computed_persent("SQ-Bulk-rhi-consensus-phylum.csv")[-12,]
aaa_1 <- aaa[,-1]
colsum.otu<-apply(aaa_1,2,sum)
t_otu_rel<-t(aaa_1)/colsum.otu
otu_rel<-t(t_otu_rel)
colSums(otu_rel)
b <- as.data.frame(aaa$Taxa)
colnames(b) <- "Taxa"
aaa_2 <-data.frame(b,as.data.frame(otu_rel))
otu_taxa <- aaa_2 %>% 
  pivot_longer(cols = !Taxa,names_to = "Samples",
               values_to = "number") %>% arrange(desc(number))

sample_group <- read.delim("sample_group.csv",check.names = F,sep=",")
sample_group_SQ<-sample_group[which(sample_group$Location=="SQ"),]
meta_taxa <- sample_group_SQ %>% 
  inner_join(.,otu_taxa,by="Samples")
write.csv(meta_taxa,file = "SQ_relative.csv")
meta_taxa$Taxa <- factor(meta_taxa$Taxa,levels = unique(meta_taxa$Taxa))
palette <-c("#00545b","#ff856d","#640025","#3ddda5","#cdffaa","#150e00","#bae278",
            "#007a98","#ffe093","#00533f","#90f0ff","#6d3c00","#004f17")

p2 <- ggplot(meta_taxa,aes(Samples,number,fill=Taxa))+
  geom_col(position="stack") +
  facet_nested(.~Location+Treatment+Species,drop=T,
               scale="free",space="free",switch="x")+
  scale_fill_manual(values=palette)+
  labs(x=NULL, y="Percent Phyla Abundance")+
  scale_y_continuous(expand = c(0,0),labels=scales::percent)+
  theme(strip.background = element_rect(fill="white",color="black"),
        panel.spacing = unit(0,"lines"),
        strip.text.x = element_text(size=12,color="black"),
        axis.text.y=element_text(size=12),
        axis.title.y = element_text(size=12,color="black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  labs(fill="Phylum");p2

g <- ggplot_gtable(ggplot_build(p2))
strips <- which(grepl('strip-', g$layout$name))
pal <- c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF",
         "#FF0000","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF",
         "#F8AFA8","#4DBBD5FF","#B09C85FF","#3C5488FF","#F39B7FFF","#B09C85FF","#91D1C2FF",
         "#D3DDDC","#00A087FF","#E6A0C4","#3C5488FF")
for (i in seq_along(strips)) {
  k <- which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  l <- which(grepl('titleGrob', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- pal[i]
}
p3 <- plot(g)
ggsave("SQ_bac_species_consist.pdf",width = 10,height = 6)



#Figure2#######################################################################################
library(dplyr)
core_filter_16s<-function(otu){
  ban.otu<-ceiling(otu/max(otu))
  occur<-rowSums(ban.otu)
  otu_filter<-otu[which(occur>=38),]
  df_rel<-apply(otu,1,sum)/sum(apply(otu,1,sum))
  df_rel<-as.data.frame(df_rel)
  rownames(df_rel)<-rownames(otu)
  df_rel$otu<-rownames(otu)
  colnames(df_rel)<-c("rel","otu")
  df_rel2<-df_rel[order(df_rel[,1],decreasing=T),]
  df_rel_core<-df_rel2[1:floor(0.1*nrow(df_rel2)),]
  df_core2<-otu[match(rownames(df_rel_core),rownames(otu)),]
  core_df<-intersect(df_core2,otu_filter)
  result<-core_df
  return(result)
}
core_16s<-core_filter_16s(df5)

core_filter_its<-function(otu){
  ban.otu<-ceiling(otu/max(otu))
  occur<-rowSums(ban.otu)
  otu_filter<-otu[which(occur>=5),]
  df_rel<-apply(otu,1,sum)/sum(apply(otu,1,sum))
  df_rel<-as.data.frame(df_rel)
  rownames(df_rel)<-rownames(otu)
  df_rel$otu<-rownames(otu)
  colnames(df_rel)<-c("rel","otu")
  df_rel2<-df_rel[order(df_rel[,1],decreasing=T),]
  df_rel_core<-df_rel2[1:floor(0.4*nrow(df_rel2)),]
  df_core2<-otu[match(rownames(df_rel_core),rownames(otu)),]
  core_df<-intersect(df_core2,otu_filter)
  result<-core_df
  return(result)
}
core_its<-core_filter_its(df6)
df_core_cs_bac<-dplyr::select(core_16s,contains("CST"))
cs_yl_bac <- dplyr::select(core_16s,contains("YL_CST"))
cs_sq_bac <- dplyr::select(core_16s,contains("SQ_CST"))
df_core_cs_its <- dplyr::select(core_its,contains("CST"))
df_core_cs <- rbind(df_core_cs_bac,df_core_cs_its)
cs_yl_bac_its <- dplyr::select(df_core_cs,contains("YL_CST"))
cs_sq_bac_its <- dplyr::select(df_core_cs,contains("SQ_CST"))

df_core_yh_bac<-dplyr::select(core_16s,contains("YHT"))
yh_yl_bac <- dplyr::select(core_16s,contains("YL_YHT"))
yh_sq_bac <- dplyr::select(core_16s,contains("SQ_YHT"))
df_core_yh_its <- dplyr::select(core_its,contains("YHT"))
df_core_yh <- rbind(df_core_yh_bac,df_core_yh_its)
yh_yl_bac_its <- dplyr::select(df_core_yh,contains("YL_YHT"))
yh_sq_bac_its <- dplyr::select(df_core_yh,contains("SQ_YHT"))

library(Hmisc)
library(igraph)
library(WGCNA)
co_network<-function(otu_all,df_x){
  b<-rcorr(t(otu_all),type="spearman")
  str(b)
  r<-b$r
  p<-b$P
  diag(p) <- 0
  p[is.na(p)] <- 0
  r[abs(r)<0.8]<-0
  p_qvalue<-qvalue(p)
  ###str(p_qvalue)
  q<-p_qvalue$qvalues
  q[q<=0.01]<-1
  q[q>0.01&q<1]<-0
  g<-r*q
  #p<-p.adjust(p,method="BH")
  diag(g)<-0
  g[is.na(g)] <- 0
  g2<-graph.adjacency(g, weighted=TRUE, mode="undirected")
  g1<-igraph::simplify(g2)
  is.simple(g1)
  g<-delete.vertices(g1,names(degree(g1)[degree(g1)==0]))
  V(g);E(g)
  
  rownames_type<-c(rownames(df_x))
  type<-c(rep("Bacteria",nrow(df_x)))
  otu_type<-data.frame(rownames_type,type)
  nodes_type<-otu_type[match(V(g)$name,otu_type[,1]),]
  #dim(nodes_type)
  V(g)$nodes_type<-as.character(nodes_type[,2])
  
  ##define the relative abundance of Vertice
  dfz<-df56_rel[match(V(g)$name,rownames(df56_rel)),]
  V(g)$rel<-as.numeric(10000*dfz$rel)
  typeof(V(g)$rel)
  
  ##node topological
  g_whole.weight<-E(g)$weight
  sum(g_whole.weight>0)#number of postive correlation
  sum(g_whole.weight<0)#number of negative correlation
  E(g)$weight<-as.numeric(100*E(g)$weight)
  g_whole.weight<-E(g)$weight#边属性
  g_whole.weight[g_whole.weight>0]<-"positive"
  g_whole.weight[g_whole.weight<0]<-"negative"
  g_whole.weight<-as.character(g_whole.weight)
  edge_attr(g, "Correlation") <- g_whole.weight
  V(g)$label<-V(g)$name
  #write.graph(g,"g.gml", format="gml")
  return(g)
}
g_cs<-co_network(df_core_cs_bac,df_core_cs_bac)
write.graph(g_cs,"cs.gml", format="gml")
g_yh<-co_network(df_core_yh_bac,df_core_yh_bac)
write.graph(g_yh,"yh.gml", format="gml")
g_yh_yl<-co_network(yh_yl_bac,yh_yl_bac)
write.graph(g_yh_yl,"yh_yl.gml", format="gml")
g_yh_sq<-co_network(yh_sq_bac,yh_sq_bac)
write.graph(g_yh_sq,"yh_sq.gml", format="gml")

#Figure3######################################################################################
df<-read.csv("df_z.csv",row.names=1)
rownames(df)[2:length(rownames(df))]
rownames(df_root) <- rownames(df)[2:length(rownames(df))]
df_root$sample <- rownames(df_root)

df_all <- merge(df_root,df_otu_filt_re_t,by = "sample")
rownames(df_all) <- df_all$sample

df_all_28<-df_all[,c(2,8)]
colnames(df_all_28) <- c("RL","Alternaria")
colnames(df_all_28)
model_1<-lm(RL~g__Alternaria,df_all_28)
summary(model_1)
p_28<-ggplot2::ggplot(df_all_28,aes(x=RL,y=Alternaria,color=Alternaria))+
  geom_point(size=1.8,alpha=0.5)+
  stat_smooth(method=lm,se=TRUE)+
  annotate("text",label="R2=0.6454 P=0.0001",x=50,y=0.2)+
  xlab("Root length (cm)")+ylab("Relative abundance of Alternaria (%)")+
  theme_bw()



#Figure4#############################################################################
df  %<>% as_tibble() %>% mutate(Group=as.factor(group))
anova <- aov(fresh_weight~Group,data=df)
summary(anova)
tukey <- TukeyHSD(anova)
library(multcompView)
cld <- multcompLetters4(anova,tukey)

dt <- df %>% group_by(Group) %>% 
  summarise(mean = mean(fresh_weight), sd = sd(fresh_weight),max = max(fresh_weight)) %>%
  arrange(desc(mean)) %>%
  ungroup() %>%
  left_join(.,as.data.frame.list(cld$Group) %>% select(1) %>% 
              rownames_to_column("Group"))
unique(df$Group)
p_jg2 <- df %>% ggplot(aes(x = factor(Group,levels = c(unique(Group)[c(8,1:7)])),y=fresh_weight)) + 
  stat_boxplot(geom="errorbar",position = position_dodge(0.5),width=0.15)+
  geom_boxplot(aes(fill = Group),width = 0.4,show.legend = T) +
  geom_point(aes(color = Group), position=position_jitterdodge(0.12),
             alpha = 0.6) +
  geom_text(data = dt, aes(x = Group, y = max, label = Letters), nudge_y = 0.01)+
  theme_bw()+
  scale_fill_manual(values = mycolor)+
  scale_color_manual(values = mycolor)+
  theme(axis.title.x = element_blank(),
        # legend.position = "none",
        panel.grid.major = element_line(linewidth  = 1),
        panel.grid.minor = element_line(linewidth = 0.35),
        panel.border = element_rect(fill = NA,colour = "black",linewidth  = 1,linetype =NULL ),
        axis.text.x = element_text(angle = 45))+
  annotate('text', label = sprintf('P == %.10f', 1.03e-07), 2.0, y = 0.6, size = 4, parse = TRUE)+
  stat_summary(fun=mean, geom="line", aes(group=1))  + 
  stat_summary(fun=mean, geom="point")+
  labs( y = "fresh_weight");p_jg2


#Figure5############################################################################
library(pheatmap)
library(readxl)
df1 <- read_xlsx("drought gene2.xlsx",sheet = "Sheet1",col_names = T)
df1$K <- paste(df1$KO,df1$`KEGG Name`,sep = " ")
df1 <- df1[,-c(1,2,15)]
df1 <- as.data.frame(df1[,sort(colnames(df1))])
rownames(df1) <- df1$K
df1 <- df1[,-1]
df1 <- log2(df1+0.0005)
df1 <-df1[apply(df1,1,var)!=0,]  
data_scale<-as.data.frame(t(apply(df1,1,scale)))
colnames(data_scale)<- colnames(df1)
col = colorRampPalette(colors = c("blue","white","red"))(80)
gaps_col = c(3,6,9)
pheatmap(data_scale,
         cluster_cols = F,
         cutree_rows = 4,
         cutree_cols = 3,
         show_colnames = F,
         color = col,
         cellwidth = 12,
         cellheight = 10,
         annotation_col = df2,
         border = F,
         treeheight_row = 30,
         gaps_col = c(3,6,9),
         filename = "pheatmap.pdf")
anno_colors<-list(Group_1=c(A="red",B="Blue",C="green"),
                  Group_2=c(A_gene="orange",B_gene="skyblue"))
pheatmap(df1,
         annotation_col = df2,
         annotation_row = df3,
         cutree_rows = 2,
         cutree_cols = 3,
         annotation_colors = anno_colors)
