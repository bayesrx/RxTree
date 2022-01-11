require(ape)
require(phylobase)
require(abc)
require(ggplot2)
require(here)

here()
source("ABC_SyntheticData.R")
source("Inference.R")
source("Tree_Summary.R")
#source("Visualization.R")

#### data pre-processing ####
rawData<-load(here("PDX_Data","rawData.rda"))
cancer.type="BRCA"
outcome = "BAR"

dat = split.cm.data$BRCA; clinical = dat[ , 1:17]
clinical$RespScaled = -clinical$RespScaled  

new.resp.mat = matrix(NA, nrow = length(unique(clinical$Model)), 
                      ncol = length(unique(clinical$Treatment)))
rownames(new.resp.mat) = unique(clinical$Model)
colnames(new.resp.mat) = unique(clinical$Treatment)
for (dim1 in 1:nrow(new.resp.mat)) {
  for (dim2 in 1:ncol(new.resp.mat)) {
    row = which(clinical$Model == rownames(new.resp.mat)[dim1] & 
                clinical$Treatment == colnames(new.resp.mat)[dim2])
    if (length(row) != 0) {
      new.resp.mat[dim1, dim2] = clinical$RespScaled[row]
    }
  }
} 
clinical = new.resp.mat 
df_out<-t(clinical)
df_out_<-df_out[,!apply(df_out,2,function(x) sum(is.na(x))/nrow(df_out)>0.4 )]
df_out_impute<-bnstruct::knn.impute(df_out_)
obsDf<-df_out_impute[rownames(df_out_impute)!="untreated",] -
  matrix(rep(df_out_impute[rownames(df_out_impute)=="untreated",],nrow(df_out_impute)-1),
         ncol=ncol(df_out_impute),byrow=T)
#saveRDS(obsDf,file=here("PDX_Data",paste0(cancer.type,"_",outcome,".RDS")))

#### ABC synthetic data ####
sim_sumstat<-ABC_SyntheticData(nLeaf=nrow(obsDf), nPDX=ncol(obsDf), nSyn= 600000)
#sim_sumstat<-readRDS(here("ABC_Synthetic_Data",paste0(cancer.type,"_sim_sumstat.RDS")))

#### ABC inference ####
#obsDf<-readRDS(here("PDX_Data",paste0(cancer.type,"_",outcome,".RDS")))
post_s2<-abc_s2( simulation_sumstat =sim_sumstat, df = obsDf, d=0.005)
post_c<-abc_c( simulation_sumstat = sim_sumstat, df = obsDf, d=0.005)
postMed_c<-summary(post_c, print=F)["Weighted Median:",1]
postMed_s2<-summary(post_s2, print=F)["Weighted Median:",1]

#### MH Stage ####
options(warn=-1)
method_vec=c("ward.D", "ward.D2", "single", "complete", "mcquitty")
MH_res_lt<-list()
MCMC_Num<-10000
for(medIdx in 1:length(method_vec)){ # each cancer runs five chains and each chain was initialized by different tree initialization
  MH_res<-Tr_twoStage(post_c=postMed_c,post_s2=postMed_s2,iteNum = MCMC_Num, obsDf = obsDf, hclust_method = method_vec[medIdx])
  MH_res_lt[[medIdx]]<-MH_res
}
#saveRDS(MH_res_lt,here("MH_posteriorTree",cancer.type,paste0(cancer.type,"_MH_res.RDS")))

#### Posterior Tree Summary: MAP ####
burnIn<-9000
MH_res_lt<-readRDS(here("MH_posteriorTree",cancer.type,paste0(cancer.type,"_MH_res.RDS")))
MAP_lt<-list() # find the MAP tree of each chain
postTr<-list() # discard the first 9,000 trees and stroe the rest trees
maxLlh_vec<-rep(NA,length(MH_res_lt))
for(medIdx in 1:length(method_vec)){
  chain_llh<-MH_res_lt[[medIdx]]$llh[(burnIn+1):MCMC_Num,]
  chain_tr<-MH_res_lt[[medIdx]]$tree_lt[(burnIn+1):MCMC_Num]
  postTr<-c(postTr,chain_tr)
  MAP_info<-getMAP(chain_llh,chain_tr)
  MAP_lt[[medIdx]]<-MAP_info$MAP_Tree
  maxLlh_vec[medIdx]<-MAP_info$Max_llh
}
MAP_Tr<-MAP_lt[[which.max(maxLlh_vec)]]
#saveRDS(MAP_Tr,here("MH_posteriorTree",cancer.type,paste0(cancer.type,"_MAP.RDS")))
#saveRDS(postTr,here("MH_posteriorTree",cancer.type,paste0(cancer.type,"_postTr.RDS")))

#### Posterior Tree Summary: Pairwise iPCP ####
postTr<-readRDS(here("MH_posteriorTree",cancer.type,paste0(cancer.type,"_postTr.RDS")))
pair_iPCP<-matrix(NA,nrow=0,ncol=3)
colnames(pair_iPCP)<-c("trt1","trt2","iPCP")
pair_iPCP<-data.frame(pair_iPCP)
obsDf<-readRDS(file=here("PDX_Data",paste0(cancer.type,"_",outcome,".RDS")))

trt_names<-rownames(obsDf)
counter<-1
for(trt1_idx in 1:(nrow(obsDf)-1)){
  for(trt2_idx in (trt1_idx+1):nrow(obsDf)){
    trt1<-trt_names[trt1_idx]
    trt2<-trt_names[trt2_idx]
    iPCP_tmp<-iPCP(treatments = c(trt1, trt2), postTr_lt = postTr)
    pair_iPCP[counter,]<-c(trt1,trt2,iPCP_tmp)
    counter<-counter+1
  }
}

iPCP_one<-data.frame(trt1=trt_names,trt2=trt_names,iPCP=1)
pair_iPCP<-rbind(pair_iPCP,iPCP_one)
saveRDS(pair_iPCP,here("MH_posteriorTree",cancer.type,paste0(cancer.type,"_pairIPCP.RDS")))

# the pairwise iPCP and the posterior tree samples are stored and visualized by the Visualization.R file.


