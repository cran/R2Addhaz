#' @title R2Addhaz
#'
#' @description This package estimates an R2 measure of explained variation under the Additive Hazards Model
#'
#' @param data=data frame of survival data. First column needs to be the censored failure time, second column needs to be the failure indicator (1=observed failure, 0=censored), the other columns are covariates.
#'
#'
#' @return R=estimated R2
#'
#' @examples
#'
#' @export R2addhaz

#R2 for one-dimensional covariate
Rsquaredone<-function(data,beta,Lambda0,time,Zbar)
{prec=0.01
n<-nrow(data)
eps=0.01
function1<-function(j)
{s1=exp(-Lambda0-c(beta%*%data[j,3])*time)
diff=c(diff(time))
idx=diff>eps
time2<-c()
idx2<-c()
for (i in which(idx))
{ seq2<-seq(time[i],time[i+1],by=eps)
time2<-c(time2,seq2[-c(1)])
idx2<-c(idx2,rep(i,length(seq2)-1))
}
time3<-time[idx2]
Zbar3<-Zbar[idx2+1]
Lambda03<-Lambda0[idx2]
s2=exp(-Lambda03-c(beta*data[j,3])*time2+c(beta*Zbar3)*(time2-time3))
s=c(s1,s2)
timef=c(time,time2)
s=s[order(timef)]
timef=timef[order(timef)]
sadj=cummin(s)
etz=caTools::trapz(timef,sadj)
norm=sadj[length(sadj)]
etz=etz/(1-norm)-timef[length(timef)]*norm/(1-norm)
etz2=caTools::trapz(timef,timef*sadj)
norm=sadj[length(sadj)]
etz2=2*etz2/(1-norm)-timef[length(timef)]^2*norm/(1-norm)
return(c(etz,etz2))}
res<-as.vector(unlist(lapply(1:n,function1)))
ETZ<-res[c(TRUE, FALSE)]
ET2Z<-res[c(FALSE, TRUE)]
varT<-mean(ET2Z)-(mean(ETZ))^2
varTZ=mean(ET2Z-(ETZ)^2)
R=1-(varTZ/varT)
differenza=10
eps=eps/2
while (differenza>0.01)
{function1<-function(j)
{s1=exp(-Lambda0-c(beta%*%data[j,3])*time)
diff=c(diff(time))
idx=diff>eps
time2<-c()
idx2<-c()
for (i in which(idx))
{ seq2<-seq(time[i],time[i+1],by=eps)
time2<-c(time2,seq2[-c(1)])
idx2<-c(idx2,rep(i,length(seq2)-1))
}
time3<-time[idx2]
Zbar3<-Zbar[idx2+1]
Lambda03<-Lambda0[idx2]
s2=exp(-Lambda03-c(beta*data[j,3])*time2+c(beta*Zbar3)*(time2-time3))
s=c(s1,s2)
timef=c(time,time2)
s=s[order(timef)]
timef=timef[order(timef)]
sadj=cummin(s)
etz=caTools::trapz(timef,sadj)
norm=sadj[length(sadj)]
etz=etz/(1-norm)-timef[length(timef)]*norm/(1-norm)
etz2=caTools::trapz(timef,timef*sadj)
norm=sadj[length(sadj)]
etz2=2*etz2/(1-norm)-timef[length(timef)]^2*norm/(1-norm)
return(c(etz,etz2))}
res<-as.vector(unlist(lapply(1:n,function1)))
ETZ<-res[c(TRUE, FALSE)]
ET2Z<-res[c(FALSE, TRUE)]
varT<-mean(ET2Z)-(mean(ETZ))^2
varTZ=mean(ET2Z-(ETZ)^2)
Rnew=1-(varTZ/varT)
differenza=abs(Rnew-R)
R=Rnew
if (differenza>0.01) {eps=eps/2}}
return(R)}

#R2 for multi-dimensional covariates
Rsquaredmult<-function(data,beta,Lambda0,time,Zbar)
{n<-nrow(data)
datacov<-data[,which(names(data)!='time' & names(data)!='status')]

prec=0.01
eps=0.01
function1<-function(j)
{s1=exp(-Lambda0-c(t(beta)%*%as.numeric(datacov[j,]))*time)
diff=c(diff(time))
idx=diff>eps
time2<-c()
idx2<-c()
for (i in which(idx))
{ seq2<-seq(time[i],time[i+1],by=eps)
time2<-c(time2,seq2[-c(1)])
idx2<-c(idx2,rep(i,length(seq2)-1))
}
time3<-time[idx2]
Zbar3<-Zbar[idx2+1,]
Lambda03<-Lambda0[idx2]
s2=exp(-Lambda03-c(t(beta)%*%as.numeric(datacov[j,]))*time2+c(t(beta)%*%t(Zbar3))*(time2-time3))
s=c(s1,s2)
timef=c(time,time2)
s=s[order(timef)]
timef=timef[order(timef)]
sadj=cummin(s)
etz=caTools::trapz(timef,sadj)
norm=sadj[length(sadj)]
if (norm==1){norm=0.99999999999}
etz=etz/(1-norm)-timef[length(timef)]*norm/(1-norm)
etz2=caTools::trapz(timef,timef*sadj)
norm=sadj[length(sadj)]
if (norm==1){norm=0.99999999999}
etz2=2*etz2/(1-norm)-timef[length(timef)]^2*norm/(1-norm)
return(c(etz,etz2))}
res<-as.vector(unlist(lapply(1:n,function1)))
ETZ<-res[c(TRUE, FALSE)]
ET2Z<-res[c(FALSE, TRUE)]
varT<-mean(ET2Z)-(mean(ETZ))^2
varTZ=mean(ET2Z-(ETZ)^2)
R=1-(varTZ/varT)

differenza=10
eps=eps/2
while (differenza>0.01)
{function1<-function(j)
{s1=exp(-Lambda0-c(t(beta)%*%as.numeric(datacov[j,]))*time)
diff=c(diff(time))
idx=diff>eps
time2<-c()
idx2<-c()
for (i in which(idx))
{ seq2<-seq(time[i],time[i+1],by=eps)
time2<-c(time2,seq2[-c(1)])
idx2<-c(idx2,rep(i,length(seq2)-1))
}
time3<-time[idx2]
Zbar3<-Zbar[idx2+1,]
Lambda03<-Lambda0[idx2]
s2=exp(-Lambda03-c(t(beta)%*%as.numeric(datacov[j,]))*time2+c(t(beta)%*%t(Zbar3))*(time2-time3))
s=c(s1,s2)
timef=c(time,time2)
s=s[order(timef)]
timef=timef[order(timef)]
sadj=cummin(s)
etz=caTools::trapz(timef,sadj)
norm=sadj[length(sadj)]
if (norm==1){norm=0.99999999999}
etz=etz/(1-norm)-timef[length(timef)]*norm/(1-norm)
etz2=caTools::trapz(timef,timef*sadj)
norm=sadj[length(sadj)]
if (norm==1){norm=0.99999999999}
etz2=2*etz2/(1-norm)-timef[length(timef)]^2*norm/(1-norm)
return(c(etz,etz2))}
res<-as.vector(unlist(lapply(1:n,function1)))
ETZ<-res[c(TRUE, FALSE)]
ET2Z<-res[c(FALSE, TRUE)]
varT<-mean(ET2Z)-(mean(ETZ))^2
varTZ=mean(ET2Z-(ETZ)^2)
Rnew=1-(varTZ/varT)
differenza=abs(Rnew-R)
R=Rnew
if (differenza>0.01) {eps=eps/2}}

return(R)}

R2addhaz<-function(data)
{data=data[order(data[,1]),]
fit<-ahaz::ahaz(survival::Surv(data[,1],data[,2]), as.matrix(data[,-c(1,2)])) #fit additive hazards model
cumahaz <- predict(fit, type="cumhaz") #estimated cumulative baseline
beta<-as.numeric(coef(fit)) #estimated coefficient
n=nrow(data)
Lambda0<-cumahaz$cumhaz #estimated cumulative baseline

#creating Zbar (mean of the covariates in the risk set)
r=n:1
datacov<-data[,-c(1,2)]
m<-ncol(datacov)
if (is.null(m)) {s_1=rev(cumsum(rev(datacov)))
Zbar=s_1/r #Zbar
Zbar<-rbind(Zbar[1],Zbar)}
if (!is.null(m)) {s_1<-datacov
for (j in 1:m){s_1[,j]=rev(cumsum(rev(datacov[,j])))}
Zbar=s_1/r #Zbar
Zbar<-rbind(Zbar[1,],Zbar)}

if(length(beta)==1){R<-Rsquaredone(data,beta,Lambda0,c(0,data[,1]),Zbar)} #computing Rsquared for model with only one covariates
if(length(beta)>1){R<-Rsquaredmult(data,beta,Lambda0,c(0,data[,1]),Zbar)} #computing Rsquared for model with more than one covariates

return(R)}

