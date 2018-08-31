celltype enrichment


celltr<-function(clusters_list,runname="",inpath='~/Dropbox/SHARED/tools/Data_to_load_CellFET/',outpath=getwd(),selection=F){ #inpath="~/Dropbox/SHARED/tools/Data_to_load_CellFET",
cat('\tNOTE:\tclusters_list only ENSG ids currently supported')
cat("\tCell background is the genes with one2one human orthologous of mice genes used to build the list of cell class enriched genes by Zeisel et al 2015(Science)\n")
  library(MetaDE)
  library('parallel')

  ###load data
 Load(paste(inpath,"/HUMQb.Rdata",sep=""))            #"HUMQb" human ENSid orthologous genes of mice background genes
 Load(paste(inpath,"/hmscDF_neuron.Rdata",sep=""))    #"hmscDF" human ENSid orthologous of mice single cell enriched by class dataframe
  
  ### create a matrix for results
  cFET=matrix(nrow=length(clusters_list), ncol=14)
  row.names(cFET)=names(clusters_list)
  colnames(cFET)=c("cell class","FET p.value","FET FDR","OR","[95% OR CI] inf","OR [95% OR CI] sup",
                      "module cell enriched genes","module out cell enriched genes",
                      "non module cell enriched genes","non module out cell enriched genes",
                      "gene names of in modules cell enriched genes","module size","cell enriched genes size",
                      "% of genes in module in cell background"
  )

  
  resMsc=list()
  for(ccl in 1:length(hmscDF)){ # ccl: cell class
    cat('\t=====================',names(hmscDF)[ccl],'=====================',ccl,' of ',length(hmscDF),'\n')
    cclENS=hmscDF[[ccl]]
    ### function to fill the matrix of results
    #for(i in 1:length(clusters_list)){
    FUNC=function(i){
      Ms=length(clusters_list[[i]]) #Ms: module size
      CB=HUMQb[,'hsapiens_homolog_ensembl_gene'] #CB: cell background
      Cs=length(cclENS) #Cs: cell enriched genes size
      MCBp=length(intersect(CB,clusters_list[[i]]))/Ms #MCBp: % of genes in module in cell background

      #cFET
      cat('\t\t',names(clusters_list)[i],'\n')
      #calculate the number Mc of module i cell enriched genes(Mc: in module AND in cell class)
      Mc=length(intersect(cclENS,clusters_list[[i]]))
      McID=paste(unlist(HUMQb[which(CB %in% intersect(cclENS,clusters_list[[i]])),'external_gene_name']),collapse=", ")
      #calculate the number NMc of remaining genes not in module but in cell class
      NMc=length(cclENS)-Mc
      #calculate the number Mnc of genes in module but not in cell class
      Mnc=length(intersect(CB,clusters_list[[i]]))-Mc
      #calculate the number NMnc of genes out of module AND not in cell class
      NMnc=length(CB)-(Mc+NMc+Mnc)
      # contingency matrice for Fisher Exact Test FET all DNMs and ns DNMs
      matr=matrix(c(Mc,NMc,Mnc,NMnc), nrow=2)
      #FET
      #FisherM=fisher.test(matr,alternative="greater")
      FisherM=fisher.test(matr)
      Fisher.p=FisherM$p.value
      Fisher.or=FisherM$estimate
      Fisher.cinf=FisherM$conf.int[1]
      Fisher.cis=FisherM$conf.int[2]
      cFET[i,]=c(names(hmscDF)[ccl],Fisher.p,NA,Fisher.or,Fisher.cinf,Fisher.cis,Mc,Mnc,NMc,NMnc,McID,Ms,Cs,MCBp)
    }
#    cfet=mclapply(1:length(clusters_list),FUNC,mc.cores=detectCores())
    cfet=lapply(1:length(clusters_list),FUNC)
    #The cfet output object of the mclapply function is a list of n vectors cFET[i,] in the good order
    for(i in 1:length(clusters_list)){
      cFET[i,]=cfet[[i]]
    }
    cFET[,"FET FDR"]=p.adjust(cFET[,"FET p.value"],method="fdr")
#    write.table(cFET, sep='\t', file=paste(outpath,"/",names(hmscDF)[ccl],"_cFET_",runname,".txt",sep=""), row.names=TRUE, quote=FALSE, col.names=NA)
    resMsc[[ccl]]=cFET
  }  
  names(resMsc)=names(hmscDF)
  if(selection==T){
    select=which(as.numeric(resMsc[[1]][,"FET FDR"]) < 0.2)
    cat("number of selected modules for ", names(resMsc)[1]," :",length(select),'\n')
    if(length(select)==1){
      SignifT=resMsc[[1]][c(select,NA),]
      SignifT=SignifT[-which(is.na(rownames(SignifT))==T),]
    }else{
      SignifT=resMsc[[1]][select,]
    }
    for(ccl in 2:length(resMsc)){
      select=which(as.numeric(resMsc[[ccl]][,"FET FDR"]) < 0.2)
      cat("number of selected modules for ", names(resMsc)[ccl]," :",length(select),'\n')
      if(length(select)==1){
        SignifT=rbind(SignifT,resMsc[[ccl]][c(select,NA),])
        SignifT=SignifT[-which(is.na(rownames(SignifT))==T),]
      }else{
        SignifT=rbind(SignifT,resMsc[[ccl]][select,])
      }          
    }
#    write.table(SignifT, sep='\t', file=paste(outpath,'/significant_cFET_',runname,'.txt',sep=''), row.names=TRUE, quote=FALSE, col.names=NA) 
  }else{
    allT=resMsc[[1]]
    for(ccl in 2:length(resMsc)){
      allT=rbind(allT,resMsc[[ccl]])
    }
#    write.table(allT, sep='\t', file=paste(outpath,'/ALL_cFET_',runname,'.txt',sep=''), row.names=TRUE, quote=FALSE, col.names=NA) 
  }
  return(as.data.frame(lapply(resMsc,function(x){x[,'FET p.value']})))     
}

