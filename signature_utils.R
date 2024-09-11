library(dplyr)
library(Seurat)
library(knitr)
library(ggplot2)
library(plotly)
library(stringr)
library(tidyverse)  
library(reshape2)
library(pheatmap)
library(dendsort)
library(RColorBrewer)
library(grid)
library(igraph)
library(MAST)
library(circlize)
library(viridis)
library(ComplexHeatmap)
library(gridBase)
library(formattable)
library(kableExtra)
library(reticulate)
#setwd(r"(C:\Users\Ron\Desktop\MHCII)")
use_python(use_python("/Users/kerenreshef/opt/anaconda3/bin/python", required = T))
source_python("/Users/kerenreshef/Desktop/TAU/PhD/Madi_lab/Carmit_sc/Signature_propo/signature_utils.py")
#source_python("/home/kerenre/UV_mouse/Signature_propo/signature_utils.py")
#source_python(r"(C:\Users\Ron\Desktop\Ligand-Receptor-Pipeline\files\signature_utils.py)")
Sys.setenv('R_MAX_VSIZE'=32000000000)
jet.colors = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
th = theme_bw() +
  theme(
    plot.title = element_text(size = 25),
    axis.text.x = element_text(size = 25, color='black'),
    axis.text.y = element_text(size = 25, color='black'),
    strip.text.x = element_text(size = 25, color='black'),
    strip.text.y = element_text(size = 25, color='black',),
    axis.title.x = element_text(size = 25, color='black'),
    axis.title.y = element_text(size = 25,angle=90),
    legend.text=element_text(size=25))

unit_sig  = function(obj,lfeature){
  rData = obj@meta.data
  rData$t1 = unlist(rData[lfeature[1]])
  rData$t2 = unlist(rData[lfeature[2]])
  
  if (mean(abs(rData$t1)) > mean(abs(rData$t2))){
    scaleFactor =  mean(abs(rData$t1)) / mean(abs(rData$t2))
    rData$t2 =  rData$t2 * scaleFactor
    rRank =   rData$t1  + rData$t2
    return(rRank)
    
  }else{
    scaleFactor =  mean(abs(rData$t2))/  mean(abs(rData$t1))
    rData$t1 = rData$t1 * scaleFactor
    rRank =   rData$t1  + rData$t2
    return(rRank)
    
  }
}
unite_sig_graph = function(obj, feature,redction = "tSNE",title = NULL){
  reduct_small_cup = tolower(redction)
  if (is.null(title)){
    rplot = FeaturePlot(obj, features = feature,pt.size=1, label = FALSE,label.size = 10,reduction = reduct_small_cup ) + scale_colour_gradientn(colours = jet.colors(20)) + th
  }
  rplot = FeaturePlot(obj, features = feature,pt.size=1, label = FALSE,label.size = 10,reduction = reduct_small_cup ) + scale_colour_gradientn(colours = jet.colors(20)) + th + ggtitle(title) 
  rtable = rplot$data
  rtable = rtable[rtable[feature] >=  mean(unlist(rtable[feature])) + 2*sd(unlist(rtable[feature])),]
  rplot = rplot + stat_density2d(data = rtable,aes(x = unlist(rtable[paste(sep = "",redction,"_1")]), y = unlist(rtable[paste(sep = "",redction,"_2")])),geom='polygon',colour='red',size=0.9, bins = 5 ,alpha =0.1)
  return(rplot)
}


#' Make signature
#' 
#' creates signature plot and p-values
#' 
#' @param obj - Seurat object
#' @param up_sig - signature of up regulated genes
#' @param idents - a list of idents if there are multiple groups else - Null
#'
#' @return list with: the objects with the signature df and the graph of the dimensional reduction of the signature
#' 
make_signature = function(obj, up_sig, idents ,is_active=FALSE, title = "signature score", down_sig = NULL, format_flag=TRUE, ingerated_flag=FALSE){
  exp = obj@assays$RNA$data  %>% as.data.frame()
  sigs_scores = signature_values(exp, up_sig, format_flag, down_sig) # wilcoxon score
  if (ingerated_flag){
    graph = obj@graphs$integrated_snn %>% as.data.frame() # knn graph
  }else{
    graph = obj@graphs$RNA_snn %>% as.data.frame()
  }
  sigs_scores = propagation(sigs_scores,graph)
  obj[["SigUint"]] = sigs_scores
  usg = unite_sig_graph(obj,"SigUint")

  # generate box plot graph
  if (!is.null(idents)){
    pvalg = pValue_grpah(obj,"SigUint",title,idents,is_active = is_active)
    if (length(idents) == 2)
    {
      return(list(obj,usg,pvalg))
    }else{
      return(list(obj,usg,pvalg[[1]],pvalg[[2]]))
    }
  }
  return(list(obj,usg))
}

make_random_signature = function(obj,idents ,is_active=FALSE,title = "signature score"){
  exp = obj@assays$RNA$data   %>% as.data.frame()
  graph = obj@graphs$RNA_snn %>% as.data.frame()
  sigs_scores = random_siganture(exp,graph)
  obj[["SigUint"]] = sigs_scores[[1]]
  usg = unite_sig_graph(obj,"SigUint")
  if (!is.null(idents)){
    pvalg = pValue_grpah(obj,"SigUint",title,idents,is_active = is_active)
    return(list(obj,usg,pvalg,sigs_scores[[2]]))
  }
  return(list(obj,usg,sigs_scores[[2]]))
}
organize_gene = function (gene){
  return(substr(gene,gregexpr(pattern ="[A-Z]",gene)[[1]][1],nchar(gene)))
}

violine_t.test = function(obj,feture,idents,slotName = "data",is_ActiceIdent = FALSE,dots = TRUE){
  if (!is_ActiceIdent){
    orig = ifelse(grepl(idents[1],names(obj@active.ident)),idents[1],idents[2])
    names(orig) = names(obj@active.ident)
    obj[["Orig"]] = orig
    obj[["Clusters"]] = obj@active.ident
    obj@active.ident = as.factor(obj$Orig)
  }
  obj[["is_ident1"]] = obj@active.ident == idents[1]
  sub_t1 = subset(obj,is_ident1)
  sub_t2 = subset(obj,is_ident1 == FALSE)
  rna1 = GetAssayData(object = sub_t1, slot = slotName) %>% as.data.frame()
  rna2 = GetAssayData(object = sub_t2, slot = slotName) %>% as.data.frame()
  t1 = as.vector(rna1[row.names(rna1)==feture,])
  t2 = as.vector(rna2[row.names(rna2)==feture,])
  p_val = t.test(t1,t2)$p.value
  lim1 = mean(as.numeric(t1)) + 2 * sd(as.numeric(t1))
  lim2 = mean(as.numeric(t2)) + 2 * sd(as.numeric(t2))
  if(dots){
    rplot1 = VlnPlot(obj, features = feture,cols = c("#007ACC","#FF6600"),slot = slotName)+
      geom_boxplot(width=0.1,fill="white",outlier.size=-1)+
      scale_y_continuous(limits = (c(0, max(lim1,lim2))))+
      labs(title =  feture, subtitle = paste("P.value = " ,round(p_val,20) )) +
      th
  }else{
    rplot1 = VlnPlot(obj, features = feture,cols = c("#007ACC","#FF6600"),slot = slotName, pt.size = 0)+
      geom_boxplot(width=0.1,fill="white",outlier.size=-1)+
      scale_y_continuous(limits = (c(0, max(lim1,lim2))))+
      labs(title =  feture, subtitle = paste("P.value = " ,round(p_val,20) )) +
      th
  }
  if (!is_ActiceIdent){
    obj@active.ident = as.factor(obj[["Clusters"]])
  }
  return(rplot1)
}
violine_f.test = function(obj,feture,idents,slotName = "data",is_ActiceIdent = FALSE,dots = TRUE){
  if (!is_ActiceIdent){
    orig = ifelse(grepl(idents[1],names(obj@active.ident)),idents[1],idents[2])
    names(orig) = names(obj@active.ident)
    obj[["Orig"]] = orig
    obj[["Clusters"]] = obj@active.ident
    obj@active.ident = as.factor(obj$Orig)
  }
  vecs = list()
  lim = c()
  for(ident in idents){
    sub_obj = subset(obj,idents = ident)
    rna = GetAssayData(object = sub_obj, slot = slotName) %>% as.data.frame()
    vec = as.vector(rna[row.names(rna)==feture,])%>% as.numeric()
    vecs[[ident]] = vec
    lim = c(lim, mean(as.numeric(vec)) + 4 * sd(as.numeric(vec)))
  }
  p_val = run_f_test(vecs)
  posthoc = post_hock_anaylsis(vecs)
  if(dots){
    rplot1 = VlnPlot(obj, features = feture,slot = slotName)+
      geom_boxplot(width=0.1,fill="white",outlier.size=-1)+
      scale_y_continuous(limits = (c(0, max(lim))))+
      labs(title =  feture, subtitle = paste("P.value = " ,round(p_val,20))) +
      th
  }else{
    rplot1 = VlnPlot(obj, features = feture,cols = c("#007ACC","#FF6600"),slot = slotName, pt.size = 0)+
      geom_boxplot(width=0.1,fill="white",outlier.size=-1)+
      scale_y_continuous(limits = (c(0, max(lim))))+
      labs(title =  feture, subtitle = paste("P.value = " ,round(p_val,20) )) +
      th
  }
  if (!is_ActiceIdent){
    obj@active.ident = as.factor(obj[["Clusters"]])
  }
  return(list(rplot1,posthoc))
}
combine_violine_grpah = function(obj,genes,idents,slotName = "counts",is_ActiceIdent = FALSE,dots = TRUE){
  plotlist = list()
  phl = list()
  flag = length(idents) == 2

  for(gene in genes){
    tryCatch(expr = {
      if (flag){
        pl = violine_t.test(obj,gene,idents,slotName,is_ActiceIdent,dots = dots)
        plotlist[[gene]] = pl
      }
      else{
        plst = violine_f.test(obj,gene,idents,slotName,is_ActiceIdent,dots = dots)
        plotlist[[gene]] = plst[[1]]
        phl[[gene]] = plst[[2]]
      }
    },error = function(e){}
    )
  }

  allplot = CombinePlots(plotlist)
  if(flag){
    return(allplot)
  }
  return(list(allplot,phl))
}
Gene_format = function(x){
  x = as.character(x)
  fl = substr(x,1,1)
  sl = substr(x,2,nchar(x))
  sl = tolower(sl)
  return(paste(fl,sl,sep = ""))
}

pValue_grpah = function(obj,feature, ytitle, idents,is_active){
  col_name = "orig.ident"
  if (is_active){
    obj[["active"]] = as.vector(obj@active.ident)
    col_name = "active"
  }
  flag = length(idents) == 2

  if (flag == TRUE){
    tdata =  obj@meta.data
    tdata =  obj@meta.data
    t1 = tdata[grepl(idents[1],tdata[[col_name]]),][feature]
    t2 = tdata[grepl(idents[2],tdata[[col_name]]),][feature]
    p_val = t.test(t1,t2)$p.value
  }else{
    vecs = list()
    for(ident in idents){
      tdata =  obj@meta.data
      vec = tdata[grepl(ident,tdata$orig.ident),][feature]
      vecs[[ident]] = vec
    }
    p_val = run_f_test(vecs)
    posthoc = post_hock_anaylsis(vecs)
  }


  rplot1 = ggplot(tdata, aes(x = tdata[[col_name]] , y = unlist(tdata[feature])))+
    #geom_boxplot(alpha =.8,outlier.shape = NA, fill = c("#EDBE2A","#00B6EB",jet.colors(length(unique(tdata$orig.ident))-2)))+
    geom_boxplot(alpha =.8,outlier.shape = NA, fill = c("#7fbf7b", "#af8dc3",jet.colors(length(unique(tdata[[col_name]]))-2)))+
    xlab("Cell source")+
    ylab("Signature score")+
    labs(title =  "Signature score", subtitle = paste("P.value = " ,round(p_val,20) )) +
    ylim(low = mean(unlist(tdata[feature])) - 2* sd(unlist(tdata[feature])), high = mean(unlist(tdata[feature])) + 2*sd(unlist(tdata[feature]))) +
    th

  if(flag == FALSE){
    return(list(rplot1,posthoc))
  }
  return(rplot1)
}

box_plot_grphs = function(obj,lfeature, ytitle, Lorig = c("NC1,shSELP2")){
  return(CombinePlots(plots = list(pValue_grpah(obj,lfeature[1],ytitle,Lorig),pValue_grpah(obj,lfeature[2],ytitle,Lorig))))
}
subClusters = function(obj,identList,dm = NULL,res = 0.5){
  idents = obj@active.ident %in% identList
  names(idents) = names(obj@active.ident)
  obj[["IsidentList"]] = idents

  sub_obj = subset(obj,subset = IsidentList)
  sub_obj = FindVariableFeatures(sub_obj, selection.method = "vst", nfeatures = 2000)
  sub_obj = RunPCA(sub_obj, features = VariableFeatures(object = sub_obj))
  plot(ElbowPlot(obj, ndims = 40))

  if(is.null(dm)){
    dm = as.numeric(readline("Enter the dims number you want"))
  }

  sub_obj = FindNeighbors(sub_obj, dims = 1:dm)
  sub_obj = FindClusters(sub_obj, resolution = res)
  sub_obj = RunTSNE(sub_obj, dims = 1:dm)

  sub_obj = RunUMAP(sub_obj,dims = 1:dm)
  return(sub_obj)

}