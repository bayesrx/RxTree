#### libraries loaded ####
library(ape)
library(phylobase)
library(mvtnorm)
library(MASS)
library(abc)

#### ABC Inference  ####
sumstat_gen<-function(df,nLeaf,nPDX){
  sumstat_o<-data.frame(matrix(NA,nrow=1, ncol=(1+5+5)))
  colnames(sumstat_o)<-c("Summary_sigma",
                         paste0("dist_rel",c(10,25,50,75,90)),
                         paste0("hclu_tip_BrLen",c(10,25,50,75,90)))
  tip_dist<-dist(df)
  S_sigma<-sum(diag(as.matrix(df) %*% as.matrix(t(df))))/(nLeaf*nPDX)
  dist_rel_pct<-quantile(as.vector(tip_dist),c(0.1,0.25,0.5,0.75,0.9))
  hclu<-as.phylo(hclust(tip_dist))
  hclu_tip_BrLen<-quantile(edgeLength(phylo4(hclu))[getEdge(phylo4(hclu),1:nLeaf,"descendant")],c(0.1,0.25,0.5,0.75,0.9))
  sumstat_o[1,c("Summary_sigma",paste0("dist_rel",c(10,25,50,75,90)),paste0("hclu_tip_BrLen",c(10,25,50,75,90)))]<-
    c(S_sigma,dist_rel_pct,hclu_tip_BrLen)
  return(sumstat_o)
}

abc_s2<-function(simulation_sumstat,df,d=0.005){
  obs_sumstat<-sumstat_gen(df,nLeaf=nrow(df),nPDX=ncol(df))
  res<-abc(target=as.matrix(obs_sumstat["Summary_sigma"]),
           sumstat=as.matrix(simulation_sumstat[,"Summary_sigma"]),
           param=as.matrix(simulation_sumstat[,"sigma_sq"]),
           tol=d,transf = rep("none",1),method ="loclinear")
  return(res)
}

abc_c<-function(simulation_sumstat,df,d=0.005){
  obs_sumstat<-sumstat_gen(df,nLeaf=nrow(df),nPDX=ncol(df))
  res<-abc(target=as.matrix(obs_sumstat[1,c(paste0("hclu_tip_BrLen",c(10,25,50,75,90)),paste0("dist_rel",c(10,25,50,75,90)))]),
           sumstat=simulation_sumstat[,c(paste0("hclu_tip_BrLen",c(10,25,50,75,90)),paste0("dist_rel",c(10,25,50,75,90)))],
           param=as.matrix(simulation_sumstat[,"c_hazard"]),
           tol=d,transf = rep("none",1),method ="loclinear")
  return(res)
}


#### MH Inference  ####
is.invertible <- function(m){
  class(try(solve(m),silent=T))=="matrix"
} 

Marginal_loc_log_llh<-function(tr_Phylo4d,cmn_Sigma2,tre_sigma_old,nPDX){
  nLeaf<-nTips(tr_Phylo4d)
  loc_<-tr_Phylo4d@data[1:nLeaf,1:nPDX]
  rownames(loc_)<-tipLabels(tr_Phylo4d)
  loc_<-loc_[order(as.numeric(substring(rownames(loc_),2))),]
  if(missing(tre_sigma_old)){
    tre_sigma<-matrix(NA,nrow=nLeaf,ncol=nLeaf)
    N_grid<-expand.grid(1:nLeaf,1:nLeaf)
    N_grid<-N_grid[N_grid$Var1>=N_grid$Var2,]
    tre_sigma[lower.tri(tre_sigma,diag = T)]<-apply(N_grid,1,function(x) nodeHeight(tr_Phylo4d,MRCA(tr_Phylo4d,paste0("N",x[1]),paste0("N",x[2])),"root"))
    tre_sigma[upper.tri(tre_sigma)]<-t(tre_sigma)[upper.tri(tre_sigma)]
  }else{
    tre_sigma=tre_sigma_old
  }
  if(!is.invertible(tre_sigma)){
    return(list(singularMat=TRUE))
  }else{
    return(list(llh=sum(mvtnorm::dmvnorm(t(loc_),mean=rep(0,nLeaf),sigma=cmn_Sigma2*tre_sigma,log=T)),MRCA_mat=tre_sigma))
  }
}

H_n<-function(n_){
  if(n_==0){
    return(0)
  }else{
    return(sum(1/1:n_))
  }
}

J_n<-function(l_v,r_v){
  return(H_n(n_=(r_v+l_v-1))-H_n(n_=(l_v-1))-H_n(n_=(r_v-1)))
}

divTime_DDT_log_llh<-function(c_,l_v,r_v,t_v){
  J=J_n(l_v=l_v,r_v=r_v)
  return(log(c_)+(c_*J-1)*log(1-t_v))
}

treeStruct_DDT_log_llh<-function(l_v,r_v){
  return(lfactorial(l_v-1)+ lfactorial(r_v-1) - lfactorial(r_v+l_v-1))
}

DDT_log_llh<-function(c_,cmn_sigma2,tr_phylo4d,trStruct_old,tre_sigma_old,nPDX_){
  nLeaf<-nTips(tr_phylo4d)
  if(missing(trStruct_old)){
    tr_phylo4<-extractTree(tr_phylo4d)
    all_divT<-c(nodeHeight(tr_phylo4,nodeLabels(tr_phylo4),from="root"),rep(1,nLeaf))
    names(all_divT)[(nLeaf+1):(2*nLeaf)]<-1:nLeaf
    from_NumPt<-c(sapply(descendants(tr_phylo4,nodeLabels(tr_phylo4),"tip"),length),rep(1,nLeaf))
    names(from_NumPt)[(nLeaf+1):(2*nLeaf)]<-1:nLeaf
    trStruct_df<-data.frame(node=names(all_divT),
                            divT=all_divT,
                            m_v=from_NumPt)
    all_edge<-getEdge(tr_phylo4)
    allEdge_dashIdx<-unlist(gregexpr(pattern="-",all_edge))
    trStruct<-data.frame(from_=as.numeric(substring(all_edge,1,allEdge_dashIdx-1)),
                         to_=as.numeric(substring(all_edge,allEdge_dashIdx+1)))
    trStruct<-trStruct[trStruct$from_!=0,]
    trStruct<-trStruct[trStruct$to_ >nLeaf,]
    trStruct<-merge(x=trStruct,y=trStruct_df[,c("node","divT")],by.x="to_",by.y="node",all.x=T)
    trStruct<-merge(x=trStruct,y=trStruct_df[,c("node","m_v")],by.x="to_",by.y="node",all.x=T)
    trStruct_dup<-merge(x=trStruct,y=trStruct[,c("from_","m_v")],by.x="to_",by.y="from_",all.x=T,all.y=F)
    trStruct_dup<-trStruct_dup[!duplicated(trStruct_dup[,1:3]),]
    colnames(trStruct_dup)[4:5]<-c("m_v","l_v")
    trStruct_dup$l_v[trStruct_dup$m_v==2]=1
    trStruct_dup$r_v<-trStruct_dup$m_v-trStruct_dup$l_v
    trStruct_dup$trStruct_llh<-apply(trStruct_dup,1,function(x) treeStruct_DDT_log_llh(l_v=x[5],r_v=x[6]))
    trStruct_dup$divT_llh<-apply(trStruct_dup,1,function(x) divTime_DDT_log_llh(c_=c_,l_v=x[5],r_v=x[6],t_v=x[3]))
  }else{
    trStruct_dup=trStruct_old
  }
  location_llh<-Marginal_loc_log_llh(tr_Phylo4d =  tr_phylo4d,cmn_Sigma2 =  cmn_sigma2, tre_sigma_old,nPDX = nPDX_)
  if(length(location_llh)==1){
    return(location_llh)
  }else{
    res_log_llh<-sum(trStruct_dup$trStruct_llh,trStruct_dup$divT_llh) +location_llh$llh
    return(list(res_log_llh=res_log_llh, trStruct=trStruct_dup, tre_sigma=location_llh$MRCA_mat))
  }
}

a_t<-function(t,c_){
  return(c_/(1-t))
}

A_t<-function(t,c_){
  res<-c_*(-1)*log(1-t)
  return(res)
}

A_inv<-function(A,c_){
  return(1-exp(-A/c_))
}

div_time<-function(m_v,t_u,c_,theta_,alpha_){
  u=runif(1)
  input=A_t(t_u, c_=c_)-exp(lgamma(m_v+1+theta_)-lgamma(m_v-alpha_))*log(1-u)
  return(A_inv(input, c_=c_))
}

mean_divTime<-function(a, t_u){
  return((1-pbeta(t_u,2,a))/((a+1)*(1-t_u)^a))
}

add_tip<-function(tr_old,div_t,to_time,to_node,label){
  tip<-list(edge=matrix(c(2,1),1,2),
            tip.label=label,
            edge.length=1-div_t,
            Nnode=1)
  class(tip)<-"phylo"
  new_tr<-bind.tree(tr_old,tip,where=to_node,position=to_time-div_t)
  return(new_tr)
}

add_root<-function(tr,root_edge,rootLabel,tipLabel){
  tip<-list(edge=matrix(c(2,1),1,2),
            node.label=rootLabel,
            tip.label=tipLabel,
            edge.length=root_edge,
            Nnode=1)
  class(tip)<-"phylo"
  rooted_Ape<-bind.tree(tip,tr,where=1)
  return(rooted_Ape)
}

Random_Detach<-function(tr_phylo4){
  nLeaf<-nTips(tr_phylo4)
  root_name<-names(rootNode(tr_phylo4))
  detachNode<-sample( setdiff(paste0("N",(1:(nNodes(tr_phylo4)+nTips(tr_phylo4)))),
                              c(root_name,names(descendants(tr_phylo4,root_name,"children")))),
                      1,F)
  paNode<-ancestor(tr_phylo4,detachNode)
  paEdge<-edgeLength(tr_phylo4)[getEdge(tr_phylo4,detachNode,"descendant")]
  paDivTime<-nodeHeight(tr_phylo4,paNode,"root")
  divTime<-nodeHeight(tr_phylo4,detachNode,"root")
  subTree_tip<-names(descendants(tr_phylo4,detachNode,"tips"))
  #if(names(paNode)==paste0("N",N+2)){}
  if(length(subTree_tip)==1){
    subTree<-detachNode
  }else{
    subTree<-as(subset(tr_phylo4,node.subtree=detachNode),"phylo")  
  }
  # remaining tree contains the singleton from thr original root 
  rmnTree_tip<-setdiff(paste0("N",1:nLeaf),subTree_tip)
  if(length(rmnTree_tip)==1){
    rmnTree<-list(edge=matrix(c(2,1),1,2),
                  node.label=paste0("N",nLeaf+1),
                  tip.label=rmnTree_tip,
                  edge.length=1,
                  Nnode=1)
    class(rmnTree)<-"phylo"
  }else{
    rmnTree<-add_root(tr=as(subset(tr_phylo4,tips.include=rmnTree_tip),"phylo"),
                      root_edge = nodeHeight(tr_phylo4,names(MRCA(tr_phylo4,rmnTree_tip)),from="root"),
                      rootLabel = paste0("N",nLeaf+1),
                      tipLabel = names(rootNode(tr_phylo4)))
  }
  return(list(subTr=subTree,rmnTr=rmnTree,paNodeLab=names(paNode),paDivT=paDivTime,divT=divTime, divLab=detachNode))
}

ReAttach_Pt<-function(rmnTree,divTime,c_,theta_,alpha_){
  rmnTree_tip<-Ntip(rmnTree)
  new_tip<-list()
  if(rmnTree_tip==1){
    t<-div_time(1,0,c_=c_ ,theta_ = theta, alpha_= alpha )
    return(list(divT=t,divRoot=rmnTree$node.label,divTo=rmnTree$tip.label,distToDiv=1-t))
  }else{
    rmnTree_phylo4<-phylo4(rmnTree,check.node.labels="keep")
    Cur_n<-rmnTree_tip
    root_time=0
    root_<-names(rootNode(rmnTree_phylo4))
    root_to<-names(descendants(rmnTree_phylo4, root_, "child"))
    edge_len<-edgeLength(rmnTree_phylo4)[getEdge(rmnTree_phylo4, root_to)]
    dist_root_to<-root_time+edge_len
    while(length(new_tip)==0){
      t<-div_time(Cur_n,root_time, c_ ,theta_, alpha_)
      if(t<dist_root_to){
        return(list(divT=t,divRoot=root_,divTo=root_to,distToDiv=dist_root_to-t))
      }
      root_div<-names(descendants(rmnTree_phylo4, root_to, type="child"))
      K=length(root_div)
      n_k<-unlist(lapply(descendants(rmnTree_phylo4, root_div, type="tip"),length))
      names(n_k)<-root_div
      prob_vec<-c(n_k-alpha_,theta_+alpha_*K)/(sum(n_k)+theta_)
      path_idx_sel<-which.max(rmultinom(1,1,prob_vec))
      if(path_idx_sel==length(prob_vec)){
        return(list(divT=dist_root_to,divRoot=root_,divTo=root_to,distToDiv=dist_root_to-t))
      }else{
        selected_node<-root_div[path_idx_sel]
        m_v<-n_k[path_idx_sel]
        if(m_v==1){
          t<-div_time(1,dist_root_to, c_,theta_, alpha_)
          return(list(divT=t,divRoot=root_to,divTo=selected_node,distToDiv=1-t))
        }else{
          root_time=dist_root_to
          Cur_n<-m_v
          root_<-root_to
          root_to<-selected_node
          edge_len<-edgeLength(rmnTree_phylo4)[getEdge(rmnTree_phylo4, root_to)]
          dist_root_to<-root_time+edge_len
        }
      }
    }
  }
}

attach_subTree<-function(subTree, rmnTree, old_divTime, old_subTr_paLab, c_, theta_, alpha_){
  if(is.character(subTree)){
    newDiv_pt<-ReAttach_Pt(rmnTree=rmnTree, c_=c_, theta_=theta_, alpha_=alpha_)
    rmnTree_phylo4<-phylo4(rmnTree,check.node.labels="keep")
    while(newDiv_pt$divT>=old_divTime){
      newDiv_pt<-ReAttach_Pt(rmnTree=rmnTree, c_=c_, theta_=theta_, alpha_=alpha_)
    }
    new_phylo4<-phylo4(add_tip(rmnTree,
                               div_t = newDiv_pt$divT, 
                               to_time = nodeHeight(rmnTree_phylo4,newDiv_pt$divTo,"root"), 
                               to_node=getNode(rmnTree_phylo4,newDiv_pt$divTo),
                               label = subTree))
    nodeLabels(new_phylo4)[is.na(nodeLabels(new_phylo4))]<-old_subTr_paLab
    return(list(new_phylo4=new_phylo4, attach_root=newDiv_pt$divRoot,attach_to=newDiv_pt$divTo,new_divT=newDiv_pt$divT))
  }else if(Ntip(rmnTree)==1){
    t<-div_time(1,0,c_=c_ ,theta_ = theta_, alpha_= alpha_ )
    while(t>=old_divTime){
      t<-div_time(1,0,c_=c_ ,theta_ = theta_, alpha_= alpha_ )
    }
    new_subTree<-add_root(tr=subTree,
                          root_edge = old_divTime - t ,
                          rootLabel = old_subTr_paLab,
                          tipLabel = names(rootNode(phylo4(subTree))))
    
    new_phylo4<-phylo4(bind.tree(rmnTree,
                                 new_subTree,
                                 where=1, 
                                 position = 1-t))
    nodeLabels(new_phylo4)[is.na(nodeLabels(new_phylo4))]<-old_subTr_paLab
    return(list(new_phylo4=new_phylo4, attach_root=rmnTree$node.label,attach_to=rmnTree$tip.label,new_divT=t))
  }else{
    newDiv_pt<-ReAttach_Pt(rmnTree=rmnTree, c_=c_, theta_=theta_, alpha_=alpha_)
    rmnTree_phylo4<-phylo4(rmnTree,check.node.labels="keep")
    while(newDiv_pt$divT>=old_divTime){
      newDiv_pt<-ReAttach_Pt(rmnTree=rmnTree, c_=c_, theta_=theta_, alpha_=alpha_)
    }
    new_subTree<-add_root(tr=subTree,
                          root_edge = old_divTime - newDiv_pt$divT ,
                          rootLabel = old_subTr_paLab,
                          tipLabel = names(rootNode(phylo4(subTree))))
    new_phylo4<-phylo4(bind.tree(rmnTree,
                                 new_subTree,
                                 where=getNode(rmnTree_phylo4,newDiv_pt$divTo), 
                                 position = newDiv_pt$distToDiv))
    nodeLabels(new_phylo4)[is.na(nodeLabels(new_phylo4))]<-old_subTr_paLab
    return(list(new_phylo4=new_phylo4, attach_root=newDiv_pt$divRoot,attach_to=newDiv_pt$divTo,new_divT=newDiv_pt$divT))
  }
}

prop_log_prob<-function(orgTree_phylo4 ,rmnTree, old_detach_PaTime, old_detach_PaLab, old_detach_divLab, new_divTime, new_AttachRoot, new_AttachTo,c_){
  if(nTips(rmnTree)==1){
    q_RmnToNew=-A_t(new_divTime,c_=c_)+log(a_t(new_divTime,c_=c_))
    q_RmnToOld=-A_t(old_detach_PaTime,c_=c_)+log(a_t(old_detach_PaTime,c_=c_))
    return(c(q_RmnToNew=q_RmnToNew, q_RmnToOld=q_RmnToOld))
  }else{
    rmnTree_phylo4<-phylo4(rmnTree)
    rmnTree_root_<-names(rootNode(rmnTree_phylo4))
    
    old_DeTachRoot<-names(ancestor(orgTree_phylo4,old_detach_PaLab))
    old_DeTachTo<-setdiff(names(descendants(orgTree_phylo4,old_detach_PaLab,"children")),old_detach_divLab)
    old_pathRoot<-unique(c(rmnTree_root_,names(shortestPath(rmnTree_phylo4,old_DeTachRoot,rmnTree_root_)),old_DeTachRoot,old_DeTachTo))
    if(length(old_pathRoot)==2){
      m_v<-length(descendants(rmnTree_phylo4,old_pathRoot[-1],"tips"))
      names(m_v)<-old_pathRoot[-1]
      divT<-old_detach_PaTime
      names(divT)<-"old_detachPaTime"
    }else{
      m_v<-sapply(descendants(rmnTree_phylo4,old_pathRoot[-1],"tips"),length)
      names(m_v)<-old_pathRoot[-1]
      divT<-c(nodeHeight(rmnTree_phylo4,names(m_v)[-length(m_v)],"root"),old_detach_PaTime)
      names(divT)<-c(names(m_v)[-length(m_v)],"old_detachPaTime")
    }
    A_t_divT<-A_t(divT,c_=c_)
    A_t_tu<-c(0,A_t_divT[-length(A_t_divT)])
    A_t_tv<-A_t_divT
    final_A_t<-sum((A_t_tu-A_t_tv)/m_v)
    frac<-sum(log(m_v[-1]/m_v[-length(m_v)]))
    final_a_t<-log(a_t(old_detach_PaTime,c_=c_)/m_v[length(m_v)])
    q_RmnToOld=sum(final_A_t,frac,final_a_t)
    
    new_pathRoot<-unique(c(rmnTree_root_,names(shortestPath(rmnTree_phylo4,new_AttachRoot,rmnTree_root_)),new_AttachRoot,new_AttachTo))
    if(length(new_pathRoot)==2){
      m_v<-length(descendants(rmnTree_phylo4,new_pathRoot[-1],"tips"))
      names(m_v)<-new_pathRoot[-1]
      divT<-new_divTime
      names(divT)<-"new_divTime"
    }else{
      m_v<-sapply(descendants(rmnTree_phylo4,new_pathRoot[-1],"tips"),length)
      names(m_v)<-new_pathRoot[-1]
      divT<-c(nodeHeight(rmnTree_phylo4,names(m_v)[-length(m_v)],"root"),new_divTime)
      names(divT)<-c(names(m_v)[-length(m_v)],"new_divTime")
    }
    A_t_divT<-A_t(divT,c_=c_)
    A_t_tu<-c(0,A_t_divT[-length(A_t_divT)])
    A_t_tv<-A_t_divT
    final_A_t<-sum((A_t_tu-A_t_tv)/m_v)
    frac<-sum(log(m_v[-1]/m_v[-length(m_v)]))
    final_a_t<-log(a_t(new_divTime,c_=c_)/m_v[length(m_v)])
    q_RmnToNew=sum(final_A_t,frac,final_a_t)
    
    return(c(q_RmnToNew=q_RmnToNew, q_RmnToOld=q_RmnToOld))
  }
}


Tr_twoStage<-function(post_c,post_s2,obsDf,iteNum,hclust_method="average"){
  llh_chk<-data.frame(idx=numeric(0),
                      Accept=numeric(0),
                      llh_Prop=numeric(0),
                      q_ToOld=numeric(0),
                      llh_old=numeric(0),
                      q_ToNew=numeric(0),
                      c=numeric(0),
                      sigma2=numeric(0))
  tree_lt<-list()
  MRCA_mat_lt<-list()
  trStruct_lt<-list()
  nLeaf<-nrow(obsDf)
  nPDX<-ncol(obsDf)
  cov_est<-as.matrix(obsDf) %*% t(as.matrix(obsDf))/nPDX
  sing_brLen<-min(cov_est)/max(diag(cov_est))
  if(sing_brLen<0){
    sing_brLen<-0.01
  }
  
  CorrHclu_tmp<-phylo4(as.phylo(hclust(as.dist(1-cor(t(obsDf),use="pa")),method = hclust_method)))
  CorrHclu_totLen<-nodeHeight(CorrHclu_tmp,1,"root")
  edgeLength(CorrHclu_tmp)<-(1-sing_brLen)*edgeLength(CorrHclu_tmp)/CorrHclu_totLen
  CorrHclu<-phylo4(add_root(tr=as(CorrHclu_tmp,"phylo"), root_edge = sing_brLen, rootLabel = paste0("N",nLeaf+1), tipLabel = paste0("N","N+2")))
  obs_hclu_phylo4d<-addData(CorrHclu,tip.data=obsDf)
  nodeLabels(obs_hclu_phylo4d)<-paste0("N",(nLeaf+1):(2*nLeaf))
  nameTbl<-data.frame(ID=paste0("N",1:nrow(obsDf)),Name=tipLabels(obs_hclu_phylo4d), stringsAsFactors = F)
  tipLabels(obs_hclu_phylo4d)<-paste0("N",1:nrow(obsDf))
  tmp_obsDf<-obsDf
  rownames(tmp_obsDf)<-nameTbl$ID[order(match(nameTbl$Name,rownames(tmp_obsDf)))]
  
  for(i in 1:iteNum){
    #print(i)
    if(i==1){
      old_phylo4d<-obs_hclu_phylo4d
    }
    old_phylo4<-extractTree(old_phylo4d)
    rnd_Detach<-Random_Detach(tr_phylo4 = old_phylo4)
    prop_newTree<-attach_subTree(subTree = rnd_Detach$subTr, rmnTree = rnd_Detach$rmnTr,
                                 old_divTime = rnd_Detach$divT, old_subTr_paLab = rnd_Detach$paNodeLab, 
                                 c_ = post_c, theta_=0, alpha_=0)
    while(nTips(prop_newTree$new_phylo4)!=nLeaf){
      prop_newTree<-attach_subTree(subTree = rnd_Detach$subTr, rmnTree = rnd_Detach$rmnTr,
                                   old_divTime = rnd_Detach$divT, old_subTr_paLab = rnd_Detach$paNodeLab, 
                                   c_ = post_c, theta_=0, alpha_=0)
    }
    q_prob<-prop_log_prob(orgTree_phylo4 = old_phylo4, rmnTree = rnd_Detach$rmnTr,
                          old_detach_PaTime = rnd_Detach$paDivT, old_detach_PaLab = rnd_Detach$paNodeLab, old_detach_divLab = rnd_Detach$divLab,
                          new_divTime = prop_newTree$new_divT, new_AttachRoot = prop_newTree$attach_root, new_AttachTo = prop_newTree$attach_to,
                          c_=post_c)
    
    new_phylo4d<-phylo4d(prop_newTree$new_phylo4,tip.data=tmp_obsDf)
    
    if(i==1){
      oldTree_info<-DDT_log_llh(c_=post_c, cmn_sigma2 = post_s2, tr_phylo4d = old_phylo4d, nPDX_ = nPDX)
      newTree_info<-DDT_log_llh(c_=post_c, cmn_sigma2 = post_s2, tr_phylo4d = new_phylo4d, nPDX_ = nPDX)
    }else{
      oldTree_info<-DDT_log_llh(c_=post_c, cmn_sigma2 = post_s2, tr_phylo4d = old_phylo4d,
                                trStruct_old = trStruct_lt[[(i-1)]] , tre_sigma_old = MRCA_mat_lt[[(i-1)]], nPDX_ = nPDX)
      newTree_info<-DDT_log_llh(c_=post_c, cmn_sigma2 = post_s2, tr_phylo4d = new_phylo4d, nPDX_ = nPDX)
    }
    
    while(length(newTree_info)==1 | nTips(prop_newTree$new_phylo4)!=nLeaf){
      prop_newTree<-attach_subTree(subTree = rnd_Detach$subTr, rmnTree = rnd_Detach$rmnTr,
                                   old_divTime = rnd_Detach$divT, old_subTr_paLab = rnd_Detach$paNodeLab, 
                                   c_ = post_c, theta_=0, alpha_=0)
      q_prob<-prop_log_prob(orgTree_phylo4 = old_phylo4, rmnTree = rnd_Detach$rmnTr,
                            old_detach_PaTime = rnd_Detach$paDivT, old_detach_PaLab = rnd_Detach$paNodeLab, old_detach_divLab = rnd_Detach$divLab,
                            new_divTime = prop_newTree$new_divT, new_AttachRoot = prop_newTree$attach_root, new_AttachTo = prop_newTree$attach_to,
                            c_=post_c)
      
      new_phylo4d<-phylo4d(prop_newTree$new_phylo4,tip.data=tmp_obsDf)
      
      if(i==1){
        oldTree_info<-DDT_log_llh(c_=post_c, cmn_sigma2 = post_s2, tr_phylo4d = old_phylo4d, nPDX_ = nPDX)
        newTree_info<-DDT_log_llh(c_=post_c, cmn_sigma2 = post_s2, tr_phylo4d = new_phylo4d, nPDX_ = nPDX)
      }else{
        oldTree_info<-DDT_log_llh(c_=post_c, cmn_sigma2 = post_s2, tr_phylo4d = old_phylo4d,
                                  trStruct_old = trStruct_lt[[(i-1)]] , tre_sigma_old = MRCA_mat_lt[[(i-1)]], nPDX_ = nPDX)
        newTree_info<-DDT_log_llh(c_=post_c, cmn_sigma2 = post_s2, tr_phylo4d = new_phylo4d, nPDX_ = nPDX)
      }
    }
    
    llh_chk[i,-2]<-c(i,newTree_info$res_log_llh,q_prob["q_RmnToOld"],
                     oldTree_info$res_log_llh,q_prob["q_RmnToNew"],
                     post_c, post_s2)
    u<-runif(1)
    alpha_Accp<-llh_chk[i,"llh_Prop"]+llh_chk[i,"q_ToOld"]-llh_chk[i,"llh_old"]-llh_chk[i,"q_ToNew"]
    if(u<exp(alpha_Accp)){
      llh_chk[i,2]<-1
      tmp_<-new_phylo4d
      tipLabels(tmp_)<-nameTbl$Name[order(match(nameTbl$ID,tipLabels(new_phylo4d)))]
      tree_lt[[i]]<-tmp_
      MRCA_mat_lt[[i]]<-newTree_info$tre_sigma
      trStruct_lt[[i]]<-newTree_info$trStruct
      old_phylo4d<-new_phylo4d
    }else{
      llh_chk[i,2]<-0
      tmp_<-old_phylo4d
      tipLabels(tmp_)<-nameTbl$Name[order(match(nameTbl$ID,tipLabels(old_phylo4d)))]
      tree_lt[[i]]<-tmp_
      MRCA_mat_lt[[i]]<-oldTree_info$tre_sigma
      trStruct_lt[[i]]<-oldTree_info$trStruct
    }
  }
  return(list(llh=llh_chk,tree_lt=tree_lt))
}

