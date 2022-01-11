library(ape)
library(phylobase)

getMAP<-function(llh_mat, postTr_lt){
  Accept_llh<-llh_mat[llh_mat$Accept==1,]
  Rej_llh<-llh_mat[llh_mat$Accept==0,]
  if(max(Accept_llh$llh_Prop) >= max(Rej_llh$llh_old)){
    MAP_llhIdx<-which(max(Accept_llh$llh_Prop)==llh_mat$llh_Prop)
    maxllh<-max(Accept_llh$llh_Prop)
  }else if(max(Accept_llh$llh_Prop) < max(Rej_llh$llh_old)){
    MAP_llhIdx<-which.max(Rej_llh$llh_old)
    maxllh<-max(Rej_llh$llh_old)
  }else if(nrow(Accept_llh)==0){
    MAP_llhIdx<-which.max(llh_mat$llh_old)
    maxllh<-max(llh_mat$llh_old)
  }else{
    stop("CHECK!")
  }
  MAPtr<-postTr_lt[[MAP_llhIdx]]
  return(list(Max_llh=maxllh, MAP_Tree=MAPtr))
}

MRCA_Time<-function(leaves, phylo4d_lt){
  MRCATime_vec<-sort(sapply(phylo4d_lt,function(x) 
    nodeHeight(x,node=MRCA(x,leaves),"root")))
  tot<-length(MRCATime_vec)
  tmp_tab<-table(MRCATime_vec)
  y_val<-c(1,(tot-cumsum(tmp_tab))/tot)
  x_val<-as.numeric(names(tmp_tab))
  sfun<-stepfun(x_val,y_val,f=0,right=F)
  return(list(df=data.frame(x=c(0,x_val),y=y_val),MRCATime_vec=MRCATime_vec,
              sfun=sfun))
}

getAUC<-function(step_df){
  auc_df<-data.frame(dx=step_df[-1,"x"]-step_df[-nrow(step_df),"x"],dy=step_df[-nrow(step_df),"y"])
  auc<-sum(auc_df$dx * auc_df$dy)
  return(auc)
}

iPCP<-function(treatments, postTr_lt){
  tmp_<-MRCA_Time(leaves = treatments, postTr_lt)
  tmp_<-rbind(tmp_$df,c(1,0))
  rownames(tmp_)<-NULL
  iPCP<-getAUC(tmp_)
  return(iPCP)
}
