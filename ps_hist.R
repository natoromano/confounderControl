# To simulated data 2000x10, do expertPS, hdPS, lasso logit, RF,
# Euclidean dist, Jaccard, Dice, Cosine similarities, Pearson and Spearman correl
# Nsim=1000
# Collect pval, smd of baseline variables
#
# 20-Oct-14 use rfPS regression instead (Malley 2012)
# 11-Sep-14 adjust rfPS by OOB error by sourcing rfFunctions.R
# 26-Aug-14 Yen Low
##############################################################################
#
# Read in the original patients file 
require(rJava)
require(Matching)
require(Epi) #clogistic
require(glmnet)
require(randomForest)
require(vegan)
require(FNN)
require(Matrix)
#require(doParallel)
#require(foreach)

#############SET PARAMETERS######################
#RScriptPath=Sys.getenv("RScriptPath") #Rscriptpath is an environment variable, path="/mnt/hgfs/Dropbox"
#ProjectPath=Sys.getenv("ProjectPath") #Project is an environment variable, path="/mnt/hgfs/projects"
RScriptPath="~/scripts"
ProjectPath="~/projects"
setwd(paste(ProjectPath,"/psm/models/sim2000x10",sep=""))
getwd()

baseDir = getwd();
tempDir = paste(baseDir,"tmp",sep="/") #output directory
patientfile=paste(baseDir,"patients.txt",sep="")
dir.create(tempDir)

# initialize the Java subsystem.  include a jdbc object
.jinit(classpath="/home/yenlow/software/pharmacoepi.jar", force.init = TRUE,parameters="-Xmx2g");
.jclassPath()

#required functions
source(paste(RScriptPath,"/R/logit.R",sep=""))
source(paste(RScriptPath,"/R/patientChar.R",sep=""))
source("../../hdpsFunctions.R")
source("../../rfFunctions.R")

id="id"
exposed="exposure" #name of exposure variable (must be string)
outcome="outcome"
expertvariables=c("x1","x2","x3","x4","x5","x6")  #exclude x7, IV
outcomemodelvariables_exclmatchingvar=c("x8","x9","x10") 
#empirical variables may be matched depending on HD-PS
empvariables=paste("x",1:10,sep="")
empvariables_cat=c("x1","x3","x5","x6","x8","x9")
empvariables_num=setdiff(empvariables,empvariables_cat)

#formula for expertPS model: exposure ~ var1 + var 2 + ...
fPS = as.formula(paste(exposed," ~ ",paste(expertvariables,collapse=" + "),sep=""))

#formula for outcome after expertPS
fmod=paste("outcome ~ exposure +", paste(outcomemodelvariables_exclmatchingvar,collapse=" + ",sep=""))

fPS
fmod

desiredOrder=c("n1","ORmatched","ORlow_matched","ORupp_matched","coeff_matched","se_matched","bias_matched","inCI_matched",
               "n0","ORadj","ORlow_adj","ORupp_adj",            "coeff_adj","se_adj",       "bias_adj","inCI_adj",
               "SMD","KS","KL","Cstat")

#to store consolidated results
Nsim=2
Ndraw=100
Nmethods=14 #14 methods (incl. 3 rdm sampling),
resultsArray=array(dim=c(Nmethods,20,Nsim))  #12methods (incl rdm downsampling), 20 metrics
#drawArray=array(dim=c(Nmethods,16,Ndraw))  #3methods (incl rdm downsampling), 20-4 metrics
pmatArray=array(dim=c(Nmethods+1,10,Nsim))  #12methods (incl rdm downsampling)+before, 10 variables
smdmatArray=pmatArray  #12methods (incl rdm downsampling), 10 variables
Noutcomesmat=matrix(nrow=Nsim,ncol=2)
time=array(dim=c(Nmethods,3,Nsim))
dimnames(time)[[1]]=c("expertPS","hdPS","lassoPS","rfPS","euclidean",
                      "jaccard","dice","cosine","pearson","spearman","lassoMV","rdmSamp1x","rdmBag","rdmJK")
dimnames(time)[[2]]=c("user","system","elapsed")
dimnames(resultsArray)[[1]]=dimnames(time)[[1]]
#dimnames(drawArray)[[1]]=dimnames(time)[[1]]
#dimnames(drawArray)[[2]]=c("n0","n1","KS","KL","Cstat","SMD",
#                           "ORadj","ORlow_adj","ORupp_adj","coeff_adj","se_adj",
#                           "ORmatched","ORlow_matched","ORupp_matched","coeff_matched","se_matched")

rfImpt=matrix(nrow=10,ncol=Nsim)  #10 rows for x1-x10
rownames(rfImpt)=paste("x",1:10,sep="")
lassoPSBeta=rfImpt
lassoMVvar=wantedVarList=matchedID=var_corOutcomeList=rdmBag=rdmJK=list()
Nexposed=trueORvec=trueSMDvec=trueCoeffvec=bestlambda=c()


#foreach(i=1:Nsim) %doPar%{
for(i in 1:Nsim){
	
cat("\n\n\n------------------------",i,"---------------------\n\n")

#generate simulated data
source("../../data/2000by10.R")
#load("../../data/2000x10.RData",verbose=T)
#rm(list=ls()[!(ls() %in% "ds10")]) #remove all objects except ds10
ds10$id=1:nrow(ds10)
nminor=min(table(ds10[,exposed]))

trueORvec[i]=trueOR
trueSMDvec[i]=trueSmd[2]
trueCoeffvec[i]=trueCoeff
Nexposed[i]=sum(ds10[,exposed])
Noutcomesmat[i,]=Noutcomes

###### expertPSM ######
#scale continuous variables (already scaled to standard normal)

#calculate estimated PS
print("############### expertPS ############")
expertPsMod=glm(fPS, data=ds10[,c(id,exposed,expertvariables)], family="binomial")
expertResults=extractResults(ps=expertPsMod$fitted.values,exposurevec=expertPsMod$y,fmod=fmod,
							data=ds10,id=id,exposed=exposed, outcome=outcome,logitFlag=TRUE,
							outfile="PShist_expertPS.pdf", verbose=TRUE,printvalues=FALSE,ylim=2.4)


####### HDPS ######################
print("############### hdPS ############")
#prepare data into format required by hdps
#read data file as single string (to use addPatientsFromBuffer)
patientheader=paste(id,exposed,outcome,sep="\t")
patientstring=paste(paste(ds10[,id],ds10[,exposed],ds10[,outcome],sep="\t"),collapse="\n")
#for reading in whole file as single string
#patientstring=readChar(patientfile, file.info(patientfile)$size)
#patientstring=gsub("\"","",patientstring)
#patients=dataset[,c("pid_org",exposed,outcome[1])] #get patient id, exposure, outcome
#write.table(patients,file=patientfile,sep="\t",col.names=T,row.names=F,na="") #output data to patient.txt;
datainstring=paste(patientheader, patientstring, sep="\n")
variables=empvariables_cat
dimdata=ds10[,c(id,variables)]

### INVOKE pharmacoepi.jar ### (handles cat variables only)
hdpsobj=hdps(datainstring,dimdata,outDir=tempDir,Nmostfreq=10,k=10,stratifyDim=FALSE,
		outfile="output_cohort.txt",FullOutput=TRUE,verbose=T,ZeroCellCorrection=F)
hdpsobj$selectedvariables #(see Fig2, Schneeweiss 2008 for recoded definitions)
wantedvar=hdpsobj$selectedvariables[grep("1Once$",hdpsobj$selectedvariables)]

var_corExposed=empvariables_num[abs(cor(ds10[,exposed],ds10[,empvariables_num]))>0.05]
var_corOutcome=empvariables_num[abs(cor(ds10[,outcome],ds10[,empvariables_num]))>0.05] #include in PS
IV=setdiff(var_corExposed,var_corOutcome) #exclude from PS

# Estimate the PS (force numerical variables into PS model)
dataPS=cbind(ds10[,c(id,exposed,var_corOutcome)],hdpsobj$hdpsdata[,wantedvar])
hdPsMod=glm(paste(exposed,"~ . -",id), data=dataPS, family="binomial")
names(hdPsMod$fitted.values)=as.character(dataPS[,id])
summary(hdPsMod$fitted.values)

hdResults=extractResults(ps=hdPsMod$fitted,exposurevec=hdPsMod$y,fmod=NULL,
                      		data=ds10,id=id,exposed=exposed, outcome=outcome,logitFlag=TRUE,
                         outfile="PShist_hdPS.pdf", verbose=TRUE,printvalues=FALSE,ylim=4.2)


####### lasso ##################
print("############### lassoPS ############")
mat=as.matrix(ds10[,paste("x",1:10,sep="")])
#penalty.factor=c(0,rep(1,ncol(xmat)-1)) #set 0 to force variable into model
#tune lambda for glmnet using 5-fold CV
lassoPsmod=cv.glmnet(xmat,ds10[,exposed],alpha=1,family="binomial",standardize=F,nfold=5)
bestlambda[i]=lassoPsmod$lambda.1se
lassoPSBeta[,i]=coeffAtlambda(lassoPsmod)[-1]  #exclude intercept

#get estimated ps
psl=unlogit(as.numeric(predict(lassoPsmod,xmat,s=lassoPsmod$lambda.1se))) #in logit form
names(psl)=ds10[,id]

lassoResults=extractResults(ps=psl,exposurevec=ds10[,exposed],fmod=NULL,
                        		data=ds10,id=id,exposed=exposed, outcome=outcome,logitFlag=TRUE,
                        		outfile="PShist_lassoPS.pdf", verbose=TRUE,printvalues=FALSE,ylim=2.5)


####### RF ##################
print("############### rfPS ############")
rfPsMod=randomForest(xmat,ds10[,exposed],ntree=100,importance=T,nodesize=100)  #reg mode for rel class freq
#ps=predict(rfPsMod, type="prob")[,2]
#ps=scales::rescale(ps,to=c(0.001,0.999))   #Method 1: rescale to avoid zeros and 1 which become Inf after logit
#ps=padjfromRF(rfPsMod)   #Method 2: oob-adjust prob from RF to avoid 100% prob
ps=rfPsMod$predicted
names(ps)=ds10[,id]

rfImpt[,i]=rfPsMod$importance[,1]
tryobj=try(extractResults(ps=ps,exposurevec=ds10[,exposed],
                         data=ds10,fmod=NULL,id=id,exposed=exposed, outcome=outcome,logitFlag=TRUE,
                         outfile="PShist_rfPS.pdf",verbose=TRUE,printvalues=FALSE,ylim=3))
if(class(tryobj)!="try-error") rfResults=tryobj else rfResults=NA
}
