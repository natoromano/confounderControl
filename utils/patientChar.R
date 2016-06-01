# Utility script wrote by Yen Low
# 
# Functions:
# 	1. matchStats: compares patient features and outputs to Excel file
#				   (numerical variables compared by Welch t-test;
#					categorical variables compared by Chisq or Fisher's test)
# 
# 03-Feb-14 Yen Low
#
# To use:
# source("/mnt/hgfs/Dropbox/scripts/R/patientChar.R")
###############################################################################
#TODO: createTable, make sure missing levels are in
#RScriptPath=Sys.getenv("RScriptPath")
#RScriptPath="~/scripts"
# source(paste(RScriptPath,"/R/utils.R",sep="")) #installnewpackage

# installnewpackage(c("xlsx"))
require(xlsx)

#compare patient characteristics
#categorical variables are compared by chisq (chisq.test)
#numerical variables are compared by Welch (t.test)
matchStats<-function(numvar,catvar,treatment,data,
                     labelsNumvar=NULL,ordervar=NULL,groups=NULL,outXlsx=NULL,verbose=TRUE){
  
  ######numerical variables
  numvardataframe=matrix(NA,ncol=8,nrow=length(numvar))
  colnames(numvardataframe)=c("label","","meanX","sdX","meanY","sdY","p.value","smd")
  for(j in 1:length(numvar)){
    if(verbose==TRUE) cat("\n","\n",numvar[j],"\n")
    x=data[data[,treatment]==0,numvar[j]]
    y=data[data[,treatment]==1,numvar[j]]
    ttest=try(t.test(x,y,alternative="two.sided"))
    smdval=smd(data,exposed=treatment,variable=numvar[j],categorical=F,verbose=verbose)
    if (class(ttest) != "try-error"){
      if(verbose==TRUE) print(ttest)
      #put in data.frame
      numvardataframe[j,3:8]=c(ttest$estimate[1],sd(x,na.rm=T),ttest$estimate[2],sd(x,na.rm=T),ttest$p.value,smdval)
    }else{  #if ttest is invalid
      numvardataframe[j,3:8]=c(rep(NA,5),smdval)
    }
  }
  #manipulate data.frame
  rownames(numvardataframe)=numvar
  if(is.null(labelsNumvar)) labelsNumvar=numvar
  numvardataframe[,1]=labelsNumvar
  if(verbose==TRUE) print(numvardataframe)
  
  ######categorical variables
  catvarlist=list()
  for(i in 1:length(catvar)){
    if(verbose==TRUE) cat("\n","\n",catvar[i],"\n")
    tab=table(data[,catvar[i]],data[,treatment],useNA="ifany")
    pctTab=prop.table(tab,2)
    mat=cbind(tab[,1],pctTab[,1],tab[,2],pctTab[,2])
    colnames(mat)=c("nX","%X","nY","%Y")
    if(sum(is.na(data[,catvar[i]]))>0) rownames(mat)[nrow(mat)]="Unknown"
    #calc overall independence pvalue from chisq or fisher's test
    if(min(tab)>=5) pvalue=chisq.test(tab)$p.value else pvalue=fisher.test(tab, workspace=2e8)$p.value
    #calc pvalue for each level/row proportions (binomial test but approach Normal)
    ngroups=colSums(tab)
    smdval=smd(data,exposed=treatment,variable=catvar[i],categorical=T,verbose=verbose)
    pval.indiv=c()
    for(j in 1:nrow(tab)){
      pval.indiv[j]=prop.test(x=tab[j,],n=ngroups)$p.value
    }
    catvarlist[[catvar[i]]]=list(p.value=c(pvalue,pval.indiv),mat=as.data.frame(mat),smd=smdval)
    if(verbose==TRUE) print(catvarlist[[i]])
  }
  
  if(!is.null(outXlsx)){
    if(is.null(groups)) groups=names(table(data[,treatment]))
    outwb <- createWorkbook()
    sheet <- createSheet(outwb, sheetName = "patientChar")
    cb=CellBlock(sheet,1,1,1000,10)
    titleStyle=CellStyle(outwb,font=Font(outwb,isBold=TRUE),alignment=Alignment(h="ALIGN_CENTER",wrapText=T))
    CB.setRowData(cb,c("","",groups[1],"",groups[2]),1,rowStyle=titleStyle)	
    CB.setRowData(cb,c("","","N or mean","% or SD","N or mean","% or SD","P value","SMD"),2,rowStyle=titleStyle)	
    row=3 #min 1
    if(is.null(ordervar)) ordervar=c(numvar,catvar)
    for(k in 1:length(ordervar)){
      if(ordervar[k] %in% numvar){  #for numvar
        CB.setRowData(cb,numvardataframe[ordervar[k],],row, showNA=F)
        row=row+1
      }else if(ordervar[k] %in% catvar){   #for catvar
        if(nrow(catvarlist[[ordervar[k]]]$mat)>2){  #category has 3 levels or more
          #fill in overall category info
          newrow=createRow(sheet,rowIndex=row)
          cellsHeader=createCell(newrow,colIndex=1:8)
          setCellValue(cellsHeader[[1,1]],ordervar[k]) #fill in category name
          setCellValue(cellsHeader[[1,7]],catvarlist[[ordervar[k]]]$p.value[1]) #fill in pvalue
          #fill in indiv level info
          #CB.setMatrixData(cb,catvarlist[[ordervar[k]]]$mat,row,1)
          addDataFrame(cbind(catvarlist[[ordervar[k]]]$mat,
                             catvarlist[[ordervar[k]]]$p.value[-1],
                             catvarlist[[ordervar[k]]]$smd),
                       sheet,col.names=F,row.names=T,startRow=row+1,startColumn=2)
          row=row+nrow(catvarlist[[ordervar[k]]]$mat)+1
        }else if(nrow(catvarlist[[ordervar[k]]]$mat)==2){
          CB.setRowData(cb,c(ordervar[k],"",catvarlist[[ordervar[k]]]$mat[2,],
                             catvarlist[[ordervar[k]]]$p.value[1],catvarlist[[ordervar[k]]]$smd[2]),row)
          row=row+1
        }
      }
    }  #end of for loop 
    saveWorkbook(outwb,file=outXlsx)
  }
  list(numvar=numvardataframe,catvar=catvarlist)
}
