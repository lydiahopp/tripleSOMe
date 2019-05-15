#' Preparation of data for multiSOMe analysis
#'
#' @param indata.exp  gene expression matrix, Should have the same colnames as indata.meth and indata.cnv
#' @param indata.meth  DNA methylation matrix, Should have the same colnames as indata.exp and indata.cnv
#' @param indata.cnv  copy number variation matrix, Should have the same colnames as indata.exp and indata.meth
#' @param modpar  see manual of modpar function
#'
#' @return modified matrix for input
#' @export

ModSOM <- function(indata.exp,indata.meth,indata.cnv,group.labels,modpar)
  {
  indata.meth=indata.meth[,colnames(indata.exp)]
  indata.cnv=indata.cnv[,colnames(indata.exp)]


  meth.type=modpar$meth.type
  mean.norm.exp=modpar$mean.norm.exp
  mean.norm.meth=modpar$mean.norm.meth
  mean.norm.cnv=modpar$mean.norm.cnv
  quant.cnv=modpar$quant.cnv
  quant.exp=modpar$quant.exp
  quant.meth=modpar$quant.meth
  M=modpar$M
  row.ids.exp=modpar$row.ids.exp
  row.ids.meth=modpar$row.ids.meth
  row.ids.cnv=modpar$row.ids.cnv
  return.ids=modpar$return.ids
  database.dataset=modpar$database.dataset
  database.biomart=modpar$database.biomart
  weight_exp=modpar$weight_exp
  weight_meth=modpar$weight_meth
  weight_cnv=modpar$weight_cnv


  if(meth.type=="27K")
  {
    data(gene_annotation_27K)
    annotation=annotation_27K

    indata.meth=indata.meth[na.omit(match(annotation[,2],rownames(indata.meth))),]

    if(any(is.na(match(annotation[,2],rownames(indata.meth)))) )
      annotation=annotation[-which(is.na(match(annotation[,2],rownames(indata.meth)))),]

    annotation_vec=annotation[,1]
    names(annotation_vec)=annotation[,2]
    indata.meth=do.call(rbind, by(indata.meth,annotation_vec, colMeans))[unique(annotation_vec),]
  }
  if(meth.type=="450K")
    {

    library(ChAMP)


    myNorm <- champ.norm(beta=indata.meth,arraytype="450K",cores=5)

    RSobject <- RatioSet(myNorm, annotation = c(array = "IlluminaHumanMethylation450k",
                                                annotation = "ilmn12.hg19"))
  }
  if(meth.type=="EPIC")
  {

    library(ChAMP)

    myNorm <- champ.norm(beta=indata.meth,arraytype="EPIC",cores=5)

    RSobject <- RatioSet(beta, annotation = c(array = "IlluminaHumanMethylationEPIC",annotation = "ilm10b4.hg19"))

  }
  if(meth.type %in% c("EPIC","450K"))
  {
    probe.features <- getAnnotation(RSobject)

    anno.table=cbind(rep(rownames(probe.features),sapply(probe.features$UCSC_RefGene_Accession,function(x){length(strsplit(x,";")[[1]])})),
                     unlist(lapply(probe.features$UCSC_RefGene_Accession,function(x){as.character(strsplit(x,";")[[1]])})),
                     unlist(lapply(probe.features$UCSC_RefGene_Group,function(x){as.character(strsplit(x,";")[[1]])})))

    anno.table=anno.table[grep("TSS",anno.table[,3]),]



    indata.meth=myNorm[anno.table[,1],]
    nm.vec=anno.table[,2]
    names(nm.vec)=anno.table[,1]

    indata.meth=do.call(rbind,by(indata.meth,nm.vec,colMeans))
    row.ids.meth="refseq_mrna"
  }

  if(any(indata.meth==0))
    indata.meth[which(indata.meth==0)]=0.00001


  if(row.ids.meth!=row.ids.exp && row.ids.cnv!=row.ids.exp)
  {
    library("biomaRt" )

    mart<-useMart(biomart = database.biomart, host ="www.ensembl.org") #   "jul2015.archive.ensembl.org"
    mart<-useDataset(database.dataset ,mart=mart)

    biomart.table = getBM( c( row.ids.meth, row.ids.exp, row.ids.cnv, return.ids ) ,row.ids.meth ,rownames(indata.meth), mart, checkFilters=F )

    if(length(which(apply(biomart.table,1,function(x){any(x=="")})==TRUE))>0)
      biomart.table= biomart.table[-which(apply(biomart.table,1,function(x){any(x=="")})==TRUE),]

    if(any(is.na(match(biomart.table[,2],rownames(indata.exp)))))
      biomart.table=biomart.table[-which(is.na(match(biomart.table[,2],rownames(indata.exp)))),]

    if(any(is.na(match(biomart.table[,3],rownames(indata.cnv)))))
      biomart.table=biomart.table[-which(is.na(match(biomart.table[,3],rownames(indata.cnv)))),]

    indata.meth=  indata.meth[biomart.table[,1],]
    indata.exp=  indata.exp[biomart.table[,2],]
    indata.cnv=  indata.cnv[biomart.table[,3],]

    rn=biomart.table[,4]
    names(rn)=biomart.table[,2]
    indata.exp=do.call(rbind,by(indata.exp,rn,colMeans))[unique(rn),]

    names(rn)=biomart.table[,1]
    indata.meth=do.call(rbind,by(indata.meth,rn,colMeans))[unique(rn),]

    names(rn)=biomart.table[,3]
    indata.cnv=do.call(rbind,by(indata.cnv,rn,colMeans))[unique(rn),]

  }else
  {
    rn=intersect(intersect(rownames(indata.exp),rownames(indata.meth)),rownames(indata.cnv))

    indata.meth= indata.meth[rn,]
    indata.exp=indata.exp[rn,]
    indata.cnv=indata.cnv[rn,]
  }


    library("biomaRt" )

    mart<-useMart(biomart = database.biomart, host ="www.ensembl.org" ) #   "jul2015.archive.ensembl.org"
    mart<-useDataset(database.dataset ,mart=mart)

    biomart.table = getBM( c( row.ids.meth,  "chromosome_name","strand" ) ,row.ids.meth ,rn, mart, checkFilters=F )

    if((length(which(biomart.table$chromosome_name=="X"))+length(which(biomart.table$chromosome_name=="Y")))>0)
    {
      indata.meth=indata.meth[- match(biomart.table[c(which(biomart.table$chromosome_name=="X"),which(biomart.table$chromosome_name=="Y")),1],rn),]

      rn=intersect(rownames(indata.exp),rownames(indata.meth))

      indata.meth= indata.meth[rn,]
      indata.exp=indata.exp[rn,]
      indata.cnv=indata.cnv[rn,]
    }



  if(M==T)
  {
    indata.meth[which(indata.meth==0)]=0.00001
    M=log10(indata.meth/(1-indata.meth))
    indata.meth=M

    cat( "\nTransform Beta into M-values\n" ); flush.console()
  }

  if(quant.exp)
  {
    indata.exp = oposSOM:::Quantile.Normalization( indata.exp )

    cat( "\nQuantile normalization indata.exp\n" ); flush.console()
  }


  if(quant.meth)
  {
    indata.meth = oposSOM:::Quantile.Normalization( indata.meth )
    cat( "\nQuantile normalization indata.meth\n" ); flush.console()

  }

  if(quant.cnv)
  {
    indata.cnv = oposSOM:::Quantile.Normalization( indata.cnv )

    cat( "\nQuantile normalization indata.cnv\n" ); flush.console()
  }



  if(mean.norm.exp)
  {
    indata.exp = indata.exp - rowMeans( indata.exp )
    cat( "\nZentralization indata.exp\n" ); flush.console()

  }


  if(mean.norm.meth)
  {
    indata.meth = indata.meth - rowMeans( indata.meth )
    cat( "\nZentralization indata.meth\n" ); flush.console()
  }

  if(mean.norm.cnv)
  {
    indata.cnv = indata.cnv - rowMeans( indata.cnv )
    cat( "\nZentralization indata.cnv\n" ); flush.console()

  }



  w_meth=mean(unlist(abs(indata.meth)))
  indata.meth= indata.meth*weight_meth/w_meth
  w_exp=mean(unlist(abs(indata.exp)))
  indata.exp=indata.exp*weight_exp/w_exp
  w_cnv=mean(unlist(abs(indata.cnv)))
  indata.cnv=indata.cnv*weight_cnv/w_cnv


  cat( "\nHarmonization\n" ); flush.console()


  colnames(indata.meth)=paste(colnames(indata.meth),"meth",sep="_")
  colnames(indata.exp)=paste(colnames(indata.exp),"exp",sep="_")
  colnames(indata.cnv)=paste(colnames(indata.cnv),"cnv",sep="_")

  indata = cbind(indata.exp,indata.meth,indata.cnv)


  group.labels = c(group.labels,group.labels,group.labels)
  names(group.labels)=colnames(indata)

  indata=as.matrix(indata[,order(group.labels)])
  group.labels=sort(group.labels)


  mod_par=list(indata,w_meth,w_exp,w_cnv,indata.meth,indata.exp,indata.cnv,group.labels,weight_exp,weight_meth,weight_cnv)
  names(mod_par)=c("indata","w_meth","w_exp","w_cnv","indata.meth","indata.exp","indata.cnv","group.labels","weight_exp","weight_meth","weight_cnv")

  return(mod_par)
}

