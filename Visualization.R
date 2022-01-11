require(ggplot2)
require(ggtree)
require(phylobase)
require(ggnewscale)
require(cowplot)
require(ggpubr)
require(stringr)
require(here)
require(reshape2)

#### Figure 5 ####

# load and pre-process known target for monotherapy
here()
target_mat<-matrix(NA,nrow=0,ncol=2)
cancer.vec<-c("BRCA","CRC","CM","PDAC","NSCLC")
for(cancerIdx in 1:length(cancer.vec)){
  ele<-readRDS(paste0(here("MH_posteriorTree",cancer.vec[cancerIdx],paste0(cancer.vec[cancerIdx],"_trtTargetGene.RDS"))))
  ele<-ele[ele$Treatment!="untreated",]
  ele<-ele[order(as.character(ele$Treatment)),]
  target_mat<-rbind(target_mat,ele)
}
target_mat$Treatment<-as.character(target_mat$Treatment)
target_mat$Treatment.target<-as.character(target_mat$Treatment.target)
target_mat<-target_mat[!duplicated(target_mat$Treatment),]
target_mat[grepl("\\+",target_mat$Treatment),"Target.Family"]<-"Combination Therapy"
selIdx<-(!grepl("\\+",target_mat$Treatment)) & (!grepl("\\,",target_mat$Treatment.target))
target_mat[selIdx,"Target.Family"]<-target_mat[selIdx,"Treatment.target"]
target_mat[is.na(target_mat$Target.Family),"Target.Family"]<-c("MAPK","PI3K","PI3K","CDK","NTRK")

Receptor_gr<-c("ALK","EGFR","ERBB2","ERBB3","ESR1","FGFR","FGFR2/4","DR5","IGF1R","MET","NTRK","PORCN","PRLR","SMO")
PI3K_MAPK_CDK_gr<-c("PI3K","PIK3CA","CDK","MAPK")
JAK_gr<-"JAK"
BRAF_gr<-"BRAF"
Other_gr<-c("chemotherapy","HSP90","IAP","PIM","TNKS","Tubulin")

target_mat$Target.Gr<-ifelse(target_mat$Target.Family %in% Receptor_gr, "Receptor",
                             ifelse(target_mat$Target.Family %in% PI3K_MAPK_CDK_gr, "PI3K-MAPK-CDK",
                                    ifelse(target_mat$Target.Family %in% JAK_gr, "JAK",
                                           ifelse(target_mat$Target.Family %in% BRAF_gr, "BRAF",
                                                  ifelse(target_mat$Target.Family %in% Other_gr,"Other",target_mat$Target.Family)))))
target_mat$Target.Gr<-factor(target_mat$Target.Gr,levels=c("PI3K-MAPK-CDK","Receptor","MDM2","JAK","BRAF","Combination Therapy","Other"))
saveRDS(target_mat,here("MH_posteriorTree","All_target.RDS"))


# generate figure 5
target_mat<-readRDS(here("MH_posteriorTree","All_target.RDS"))

cancer.vec<-c("BRCA","CRC","CM","NSCLC","PDAC")
g_plot<-list()
for(cancerIdx in 1:length(cancer.vec)){
  cancer.type<-cancer.vec[cancerIdx]
  MAP_phylo4d<-readRDS(here("MH_posteriorTree",cancer.type,paste0(cancer.type,"_MAP.RDS")))
  
  ### extract the order of MAP tree to re-order the pairwise iPCP and correlation
  apeMAP<-as(extractTree(MAP_phylo4d),"phylo") 
  is_tip <- apeMAP$edge[,2] <= length(apeMAP$tip.label)
  ordered_tips <- apeMAP$edge[is_tip, 2]
  
  iPCP_mat<-readRDS(here("MH_posteriorTree",cancer.type,paste0(cancer.type,"_pairIPCP.RDS")))
  iPCP_mat_sym<-iPCP_mat[iPCP_mat$trt1 != iPCP_mat$trt2, c("trt2","trt1","iPCP")]
  colnames(iPCP_mat_sym)<-c("trt1","trt2","iPCP")
  iPCP_mat<-rbind(iPCP_mat,iPCP_mat_sym)
  iPCP_mat$trt.x<-factor(iPCP_mat$trt2,levels=apeMAP$tip.label[ordered_tips])
  iPCP_mat$trt.x_num<-as.numeric(iPCP_mat$trt.x)
  
  raw_df<-readRDS(here("PDX_Data",paste0(cancer.type,"_BAR.RDS")))
  smpCorr_sq<-matrix(NA,nrow=nrow(raw_df),ncol=nrow(raw_df))
  rownames(smpCorr_sq)<-colnames(smpCorr_sq)<-rownames(raw_df)
  for(idx1 in 1:(nrow(raw_df)-1)){
    for(idx2 in (idx1+1):(nrow(raw_df))){
      trt1<-rownames(smpCorr_sq)[idx1]
      trt2<-rownames(smpCorr_sq)[idx2]
      tmpCor<-cor(as.numeric(as.matrix(raw_df[idx1,-1])),as.numeric(as.matrix(raw_df[idx2,-1])))
      smpCorr_sq[idx1,idx2]<-smpCorr_sq[idx2,idx1]<-tmpCor
    }
  }  
  diag(smpCorr_sq)<-1
  smpCorr_sq_reSc<-(smpCorr_sq+1)/2
  smpCorr_reSc<-reshape2::melt(smpCorr_sq_reSc)
  smpCorr_reSc$trt.x<-factor(smpCorr_reSc$Var2,levels=apeMAP$tip.label[ordered_tips])
  smpCorr_reSc$trt.x_num<-as.numeric(smpCorr_reSc$trt.x)
  
  p <- ggtree(extractTree(MAP_phylo4d),ladderize = F) 
  if(cancer.type=="BRCA" | cancer.type=="PDAC"){
    p <- p %<+% target_mat + geom_tippoint(aes(shape=Target.Gr),size=6) + 
      geom_tiplab(offset = .25, hjust = .3,size=7.5) + 
      scale_shape_manual("Target Pathways",guide = guide_legend(order = 1),values=c(15,1,17,18,12,3))
  }else if(cancer.type=="CRC"| cancer.type=="CM"){
    p <- p %<+% target_mat + geom_tippoint(aes(shape=Target.Gr),size=6) + 
      geom_tiplab(offset = .25, hjust = .3,size=7.5) + 
      scale_shape_manual("Target Pathways",guide = guide_legend(order = 1),values=c(15,1,17,16,12,3))
  }else if(cancer.type=="NSCLC"){
    p <- p %<+% target_mat + geom_tippoint(aes(shape=Target.Gr),size=6) + 
      geom_tiplab(offset = .25, hjust = .3,size=7.5) + 
      scale_shape_manual("Target Pathways",guide = guide_legend(order = 1),values=c(15,1,17,12,3))
  }
  p1<-facet_plot(p=p,panel = "Pairwise iPCP", data = iPCP_mat, geom = geom_tile,
                 mapping=aes(x = trt.x_num, fill = iPCP)) + labs(fill="iPCP") +
    scale_fill_gradientn(colors = RColorBrewer::brewer.pal(3, "YlOrRd"))
  p2 <- p1 + ggnewscale::new_scale_fill()
  p3<-facet_plot(p=p2 + xlim_tree(1.7) ,panel = "Rescaled Pearson Correlation", data = smpCorr_reSc, geom = geom_tile,
                 mapping=aes(x = trt.x_num, fill = value)) +labs(fill="Correlation")+
    scale_fill_gradientn(colors = RColorBrewer::brewer.pal(3, "PuBu"),guide = guide_legend(order = 3))
  p_final<-p3+theme(legend.position=c(.05, .6),text=element_text(size=27))
  p_final<-facet_labeller(p_final, c(Tree = "MAP Rx Tree"))
  g_plot[[cancerIdx]]<-p_final
  names(g_plot)[cancerIdx]<-cancer.type
}

g_final<-plot_grid(g_plot$BRCA, g_plot$CRC,g_plot$CM, nrow = 3,labels = c("BRCA", "CRC","CM"),label_size = 25)
ggsave(plot=g_final,width=34,height=34,file=here("Figure5.png"),dpi=320)

# g_final_s<-plot_grid(g_plot$NSCLC, g_plot$PDAC, nrow = 2,labels = c("NSCLC", "PDAC"),label_size = 25)
# ggsave(plot=g_final_s,width=36,height=24,file="~/Documents/Research_Veera/RxTree/ggtree2_supp.png")


#### Figure 6 ####
cancer.vec<-c("BRCA","CRC","CM","PDAC","NSCLC")
g_lt<-list()
iPCP_lt<-list()
for(cancerIdx in 1:length(cancer.vec)){
  print(cancerIdx)
  cancer.type<-cancer.vec[cancerIdx]
  iPCP_mat<-readRDS(here("MH_posteriorTree",cancer.type,paste0(cancer.type,"_pairIPCP.RDS")))
  iPCP_mat$name=paste0(iPCP_mat$trt1,", ",iPCP_mat$trt2)
  iPCP_threshold<-0.7
  iPCP_mat<-iPCP_mat[iPCP_mat$iPCP>=iPCP_threshold & iPCP_mat$trt1 != iPCP_mat$trt2,]
  iPCP_mat$cmb<-ifelse(str_detect(iPCP_mat$trt1,"\\+") | str_detect(iPCP_mat$trt2,"\\+"), "Combination Therapy","Monotherapy")
  
  iPCP_mono<-iPCP_mat[!str_detect(iPCP_mat$trt1,"\\+") & !str_detect(iPCP_mat$trt2,"\\+"),]
  iPCP_mono$name<- sapply(iPCP_mono$name,function(x) paste(sort(unlist(str_split(x,", "))),collapse = ", "))
  iPCP_cmb<-iPCP_mat[str_detect(iPCP_mat$trt1,"\\+") & str_detect(iPCP_mat$trt2,"\\+"),]
  iPCP_cmb$name<- sapply(iPCP_cmb$name,function(x) paste(sort(unlist(str_split(x,", "))),collapse = ", "))
  iPCP_final<-rbind(iPCP_mono,iPCP_cmb)
  iPCP_lt[[cancerIdx]]<-iPCP_final
  names(iPCP_lt)[cancerIdx]<-cancer.type
  g<-ggplot(iPCP_final,aes(y=reorder(name,iPCP),x=iPCP,fill=cmb))+geom_col()+  
    labs(y="Treatments Pairs",fill="Treatment Type",x="iPCP") + coord_cartesian(xlim=c(0,1)) +
    theme(text = element_text(size=30,face="bold"))
  g_lt[[cancerIdx]]<-g
  names(g_lt)[cancerIdx]<-cancer.type
}

col_tab<-data.frame(Colour=c("Blue","Yellow","Green","Orange","Grey"),
                    Target=c("PI3K", "MAPK", "(PI3K, MAPK) or (CDK, MAPK) or (PI3K, CDK)", "MDM2", 
                             "(Tubulin, PI3K-MAPK-CDK)\n or (ERBB3, PI3K-MAPK-CDK)"))
colnames(col_tab)<-c("Block Color","Treatment Target Pathways")
stable<-ggtexttable(col_tab, rows = NULL, theme = ttheme("mBlue",base_size=30))

g_final<-ggarrange(g_lt$BRCA,g_lt$CM,g_lt$CRC,stable,common.legend = T,
                   ncol=2,nrow=2,heights=c(max(nrow(iPCP_lt$BRCA),nrow(iPCP_lt$CM))/nrow(iPCP_lt$CRC),1),
                   labels=c("(A) BRCA","(B) CM","(C) CRC", " "),vjust=c(10,10,-0.2,-0.2),font.label = list(size = 30))
ggsave(file=here("Figure6.png"),plot=g_final,width=36,height=28,dpi=600)

