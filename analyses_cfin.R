#R script to analyze data on C. finmarchicus growth and development in response to food level and presence of predator cues

#Kristina Kvile, July 2020

#Load libraries and data file:
library(mgcv)
library(plyr)
load("img.rda")

#Experimental treatments:
treatments<-levels(img$treatment)
#Define stages for analyses (a few observations of C3 are removed):
stages<-levels(img$stage)[-1]

#Set color scales for plots:
col_stages<-data.frame(stages=levels(img$stage),
                       col=c("#ffffcc","#a1dab4","#41b6c4","#2c7fb8","#253494")) 
col_treatments<-data.frame(treatments=treatments,
                           col=c("#4daf4a","#984ea3","#377eb8","#e41a1c"))
img$col<-col_treatments$col[match(img$treatment,col_treatments$treatments)]


####Analyses of development####

#Calculate mean developmental stage per day and tank:
img_mean<-ddply(img,.(day,tank),summarise,stage_num=mean(stage_num),
                food=unique(food),predation=unique(predation))

#Set predictor variables:
yvar<-list(
  "food",
  "predation",
  "s(tank,bs='re')",
  "s(day,k=4,by=food)",
  "s(day,k=4,by=predation)"
)

#Set response variable:
Xvar<-"stage_num~"

#Fit GAM:
mod_dev<-gam(formula(paste(Xvar,paste(yvar,collapse="+"))),data=img_mean)

#Create table for R2, factor coefficients and P-values:
modtab_dev<-rep(NA,5+length(yvar))
names(modtab_dev)<-c("R2","Food","Predation",yvar[1:3],
                     "DayxLoF","DayxHiF","DayxLoP","DayxHiP")
modtab_dev[1]<-summary(mod_dev)$r.sq #R2
modtab_dev[2:3]<-summary(mod_dev)$p.table[2:3,1] #Factor terms coefficients
modtab_dev[unlist(yvar[1:2])]<-summary(mod_dev)$p.table[2:3,4] #P-values factor terms 
modtab_dev[unlist(yvar[3])]<-summary(mod_dev)$s.table[1,4] #P-values simple smooth term
modtab_dev[c("DayxLoF","DayxHiF","DayxLoP","DayxHiP")]<-summary(mod_dev)$s.table[2:5,4] #P-values varying coefficient terms

#Make predictions of stage number for new data frame 
#'newdat' contains all combinations of day and treatment (tank is random):
newdat<-data.frame(day=rep(1:24,4),food=rep(levels(img$food),each=24*2),
                   predation=rep(rep(levels(img$predation),each=24),2),tank="E")
pred<-predict.gam(mod_dev,newdata=newdat,type="response",se.fit=T)
newdat$dev<-pred$fit #Predicted value
newdat$dev_hi<-pred$fit+1.96*pred$se.fit #Upper confidence interval
newdat$dev_lo<-pred$fit-1.96*pred$se.fit #Lower confidence interval

#Plot proportions per stage and predictions together:
layout(matrix(c(1,2,5,3,4,6),nrow=2,ncol=3,byrow=TRUE),widths=c(0.3,0.3,0.4))
par(mar=c(0.5,0.15,1.5,0.7),oma=c(3,3,0,0),lwd=1.5)
for(x in 1:length(treatments)){
  subdat<-img[img$treatment==treatments[x],] #Data for treatment x
  freq<-table(subdat$stage,by=subdat$day) #Observations of stages per day
  perc<-as.data.frame(apply(freq,2,function(x) x/sum(x))) #Percentages of stages per day
  missing<-seq(1,24)[!seq(1,24)%in%subdat$day] #Days without observations
  perc[,as.character(missing)]<-NA #Add columns for without observations with NAs
  perc<-as.matrix(perc[,order(as.numeric(colnames(perc)))]) #Order by day
  #Plot per treatment, with/without x and y axes:
  if(x==1){barplot(perc,legend=FALSE,col=adjustcolor(col_stages$col),xaxt="n",border=NA,axis.lty=1)}
  if(x==2){barplot(perc,legend=FALSE,col=adjustcolor(col_stages$col),xaxt="n",yaxt="n",border=NA,axis.lty=1)}
  if(x==3){barplot(perc,legend=FALSE,col=adjustcolor(col_stages$col),border=NA,axis.lty=1)}
  if(x==4){barplot(perc,legend=FALSE,col=adjustcolor(col_stages$col),yaxt="n",border=NA,axis.lty=1)}
  mtext(side=3,paste0(LETTERS[x],"."),outer=F,line=0.25,adj=0.15)
  mtext(side=3,treatments[x],col=adjustcolor(col_treatments$col[col_treatments$treatments==treatments[x]]),
        outer=F,line=0.25,adj=0.35)
  if(x%in% c(3,4)){ mtext(side=1,"Day",outer=FALSE,line=2)}
}
mtext(side=2,"Fraction of stages",outer=TRUE,line=1.75)

#Add legend:
plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='')
legend("left",legend=levels(img$stage),horiz=F,bty="n",
       pch=22,pt.bg=adjustcolor(col_stages$col),cex=1.25)

#Plot predictions per treatment:
par(mar=c(0.5,2,0,0.25))
plot(1,type="n",xlim=range(img$day),ylim=c(4,6),ylab="",xlab="")
for(x in 1:length(treatments)){
  subdat<-newdat[newdat$food==unique(img$food[img$treatment==treatments[x]]) &
                   newdat$predation==unique(img$predation[img$treatment==treatments[x]]),]  
  polygon(x=c(subdat$day,rev(subdat$day)),y=c(subdat$dev_lo,rev(subdat$dev_hi)),
          col=adjustcolor(col_treatments$col[col_treatments$treatments==treatments[x]],alpha.f=0.2),border=NA)
  points(subdat$day,subdat$dev,type="l",lty=x,lwd=1.5,col=adjustcolor(col_treatments$col[col_treatments$treatments==treatments[x]]))
}
legend("topleft",legend=treatments,horiz=F,bty="n",cex=1.5,
       col=adjustcolor(col_treatments$col),pch=20,lty=1:4,x.intersp = 0.2)
mtext(side=1,"Day",outer=FALSE,line=2)
mtext(side=2,"Mean stage",outer=FALSE,line=2)
mtext(side=3,"E. Model predictions",outer=FALSE,line=0.5)


####Analyses of growth#####

#Select variables for analyses:
vars<-colnames(img)[c(9,11,15,21)]
varnames<-c("Prosome area","Lipid fullness","C:N","ln(RNA:DNA)")

#Define predictor variables:
yvar<-list(
  "food",
  "predation",
  "s(tank,bs='re')",
  "s(day,k=4,by=food)",
  "s(day,k=4,by=predation)"
)

#Create table for R2, factor coefficients and P-values:
modtab_growth<-array(NA,dim=c(length(stages),5+length(yvar),length(vars)))
dimnames(modtab_growth)<-list(stages,c("R2","Food","Predation",yvar[1:3],
                                       "DayxLoF","DayxHiF","DayxLoP","DayxHiP"),varnames)

#Fit models and fill table for all response variables and stages:
for(x in 1:length(vars)){
  for(stage in stages){
    subdat<-img[img$stage==stage & !is.na(img[,vars[x]]),] #Subset data for stage
    #For table, set response variable to standardized within stage:
    subdat[,vars[x]]<-(subdat[,vars[x]]-mean(subdat[,vars[x]],na.rm=T))/
      sd(subdat[,vars[x]],na.rm=T)
    Xvar<-paste0(vars[x],"~") #Define response variable
    #Fit GAM if there are data available for 3 unique sampling days:
    if(length(unique(subdat$day))>3){
      gam_xvar<-gam(formula(paste(Xvar,paste(yvar,collapse="+"))),data=subdat)
      modtab_growth[stage,1,varnames[x]]<-summary(gam_xvar)$r.sq #R2
      modtab_growth[stage,2:3,varnames[x]]<-summary(gam_xvar)$p.table[2:3,1] #Factor term coefficients
      modtab_growth[stage,unlist(yvar[1:2]),varnames[x]]<-summary(gam_xvar)$p.table[2:3,4] #P-values factor terms 
      modtab_growth[stage,unlist(yvar[3]),varnames[x]]<-summary(gam_xvar)$s.table[1,4] #P-values simple smooth terms 
      modtab_growth[stage,7:10,varnames[x]]<-summary(gam_xvar)$s.table[2:5,4] #P-values varying coefficient terms
    }
  }
}

#Create empty dataframe to store model predictions:
#'newdat' contains all combinations of day and treatment (tank is random):
newdat<-data.frame(day=rep(1:24,4),food=rep(levels(img$food),each=24*2),
                   predation=rep(rep(levels(img$predation),each=24),2),tank="G")

#Plot variation in observations + model predictions 
panels<-array(c(LETTERS[5:6],c(NA,NA),LETTERS[7:16]), #Panel numbering
              dim=c(length(vars),length(stages)),
              dimnames = list(vars,stages))

xvals<-c(0.5,5.5,11,16.5) #X-axis location per stage for boxplot

layout(matrix(c(1:20),nrow=4,ncol=5,byrow=TRUE),widths=c(0.42,rep(0.145,4)))
par(mar=rep(0,4),oma=c(3.5,3.5,3,1))
for(x in 1:length(vars)){
  ylims<-range(img[,vars[x]],na.rm=T) #Limits of y-axis depends on variable
  #1: Make boxplot+stripchart per response variable:
  plot(1,type="n",xlim=c(0,20),ylim=ylims,axes=F,frame=T,ylab="",xlab="")
  for(i in 1:length(stages)){
    xvals_i<-seq(xvals[i],xvals[i]+3,length.out=4) #Set positions along x-axis
    subdat<-img[img$stage==stages[i],] #Subset data for stage
    stripchart(subdat[,vars[x]]~subdat$treatment,add=T,vertical=T,pch=20,cex=0.4,col=adjustcolor(col_treatments$col),
               method="jitter",jitter=0.1,ylim=ylims,axes=F,frame=T,at=xvals_i,lwd=0)
    boxplot(subdat[,vars[x]]~subdat$treatment,outline=FALSE,add=T,col=adjustcolor(col_treatments$col,alpha.f=0.2),
            lwd=0.5,axes=F,frame=T,at=xvals_i)
  }
  mtext(side=2,varnames[x],outer=F,line=2)
  mtext(side=3,paste0(LETTERS[x],"."),outer=F,line=-2,adj=0.05,col="grey40")
  Axis(side=2)
  if(vars[x]==vars[4]){Axis(side=1,at=c(2,7,12.5,18),labels=stages);
    mtext(side=1,"Stages",line=2.25)}
  if(vars[x]==vars[1]){mtext(side=3,"A. Overall variation",outer=F,line=0)}
  #2: Plot model predictions per stage:
  Xvar<-paste0(vars[x],"~") #Define response variable
  for(stage in stages){
    subdat<-img[img$stage==stage & !is.na(img[,vars[x]]),] #Subset data per stage
    #Draw empty plot if less than 3 days with data:
    if(length(unique(subdat$day))<=3){ 
      plot(1,type="n",xlim=c(0.5,4.5),ylim=ylims,axes=F,frame=F,ylab="",xlab="")
      if(vars[x]==vars[4]){
        legend("top",legend=treatments,horiz=F,bty="n",cex=1.5,
               col=adjustcolor(col_treatments$col),pch=20,lty=1:4,x.intersp = 0.2)}
    }
    #Otherwise, draw predictions:
    if(length(unique(subdat$day))>3){
      xlims<-range(img$day[img$stage==stage]) #Limits of x-axis are stage-specific
      plot(1,type="n",xlim=xlims,ylim=ylims,axes=F,frame=T,ylab="",xlab="") #Plot background
      gam_xvar<-gam(formula(paste(Xvar,paste(yvar,collapse="+"))),data=subdat) #GAM fit
      pred<-predict.gam(gam_xvar,newdata=newdat,type="response",se.fit=T) #Model predictions
      newdat$pred<-pred$fit
      newdat$pred_hi<-pred$fit+1.96*pred$se.fit
      newdat$pred_lo<-pred$fit-1.96*pred$se.fit
      #Plot predictions for different treatments:
      for(i in 1:length(treatments)){
        daylim<-range(subdat$day[subdat$treatment==treatments[i]]) #Only plot predictions for days with observations in treatment
        sub_newdat<-newdat[newdat$day>=min(daylim) & newdat$day<=max(daylim) &
                             newdat$food==unique(img$food[img$treatment==treatments[i]]) &
                             newdat$predation==unique(img$predation[img$treatment==treatments[i]]),]  
        polygon(x=c(sub_newdat$day,rev(sub_newdat$day)),y=c(sub_newdat$pred_lo,rev(sub_newdat$pred_hi)),
                col=adjustcolor(col_treatments$col[col_treatments$treatments==treatments[i]],alpha.f=0.2),border=NA)
        points(sub_newdat$day,sub_newdat$pred,type="l",lty=i,lwd=1.5,col=adjustcolor(col_treatments$col[col_treatments$treatments==treatments[i]]))
      }
      mtext(side=3,paste0(panels[vars[x],stage],"."),outer=F,line=-2,adj=0.05,col="grey40")
      if(vars[x]==vars[1]){mtext(side=3,stage,outer=F,line=0)}
      if(vars[x]==vars[4] & stage!="C5" | vars[x]==vars[2] & stage=="C4"){
        Axis(side=1,at=seq(min(xlims),max(xlims-2),length.out=3))}
      if(vars[x]==vars[4] & stage=="C5"){
        Axis(side=1,at=seq(min(xlims),max(xlims-3),length.out=3))}
    }
  }
}
mtext(side=1,"Day",outer=T,line=2.25,adj=0.75)
mtext(side=3,"B. Temporal variation",outer=T,line=1.75,adj=0.8)
#Draw line on top of plot
par(mfrow=c(1,1),xpd=NA,new=TRUE) 
plot(1:10,1:10,type="n",axes=FALSE,xlab="",ylab="",main="")
lines(c(4.55,10.45),c(11.35,11.35))

####Supplementary analyses: Model diagnostics#####

#Model diagnostics for response variables, including RNA:DNA before/after log-transformation

#Define predictor variables:
yvar<-list(
  "food",
  "predation",
  "s(tank,bs='re')",
  "s(day,k=4,by=food)",
  "s(day,k=4,by=predation)"
)

#For mean developmental stage (Eq. 1):
par(mfrow=c(1,5),mar=rep(1.5,4),oma=c(1,2,4,0),xpd=FALSE)
# Distribution of data:
hist(img_mean$stage_num,main="")
mtext(side=3,"A. Data")
#Plot model diagnostics:
resid <- residuals(mod_dev)
linpred <- napredict(mod_dev$na.action, mod_dev$linear.predictors)
observed.y <- napredict(mod_dev$na.action, mod_dev$y)
qq.gam(mod_dev)
mtext(side=3,"B. Q-Q plot")
hist(resid, xlab = "", main = "")
mtext(side=3,"C. Residuals")
plot(linpred, resid)
mtext(side=3,"D. Residuals vs.",line=1)
mtext(side=3,"linear pred.")
plot(fitted(mod_dev), observed.y)
mtext(side=3,"E. Response vs.",line=1)
mtext(side=3,"fitted")
mtext(side=3,"Model diagnostics: Mean developmental stage",outer=TRUE,line=2,font=2)

#For growth (Eq. 2):
vars<-colnames(img)[c(9,11,15,20,21)]
varnames<-c("Prosome area","Lipid fullness","C:N","RNA:DNA","log(RNA:DNA)")

#Plots for body area and lipid fraction (all stages):
for(x in 1:2){
  par(mfrow=c(length(stages),5),mar=rep(1.5,4),oma=c(1,2,4,0))
  for(stage in stages){
    #Plot distribution of data:
    subdat<-img[img$stage==stage & !is.na(img[,vars[x]]),]
    hist(subdat[,vars[x]],main="")
    if(stage==stages[1]){mtext(side=3,"A. Data")}
    mtext(side=2,outer=F,stage,line=2)
    #Fit GAM:
    Xvar<-paste0(vars[x],"~")
    gam_xvar<-gam(formula(paste(Xvar,paste(yvar,collapse="+"))),data=subdat)
    #Plot model diagnostics:
    resid <- residuals(gam_xvar)
    linpred <- napredict(gam_xvar$na.action,gam_xvar$linear.predictors)
    observed.y <- napredict(gam_xvar$na.action,gam_xvar$y)
    qq.gam(gam_xvar)
    if(stage==stages[1]){mtext(side=3,"B. Q-Q plot")}
    hist(resid,xlab="",main="")
    if(stage==stages[1]){mtext(side=3,"C. Residuals")}
    plot(linpred,resid)
    if(stage==stages[1]){mtext(side=3,"D. Residuals vs.",line=1);
      mtext(side=3,"linear pred.")}
    plot(fitted(gam_xvar),observed.y)
    if(stage==stages[1]){mtext(side=3,"E. Response vs.",line=1);
      mtext(side=3,"fitted")}
  }
  mtext(side=3,paste0("Model diagnostics: ",varnames[x]),outer=TRUE,line=2,font=2)
}

# C:N and RNA:DNA / RNA:DNA-log  (excluding C4)
for(x in 3:5){
  par(mfrow=c(length(stages[-1]),5),mar=rep(1.5,4),oma=c(1,2,4,0))
  for(stage in stages[-1]){
    #Plot distribution of data:
    subdat<-img[img$stage==stage & !is.na(img[,vars[x]]),]
    hist(subdat[,vars[x]],main="")
    if(stage==stages[2]){mtext(side=3,"A. Data")}
    mtext(side=2,outer=F,stage,line=2)
    #Fit GAM:
    Xvar<-paste0(vars[x],"~")
    gam_xvar<-gam(formula(paste(Xvar,paste(yvar,collapse="+"))),data=subdat)
    #Plot model diagnostics:
    resid <- residuals(gam_xvar)
    linpred <- napredict(gam_xvar$na.action,gam_xvar$linear.predictors)
    observed.y <- napredict(gam_xvar$na.action,gam_xvar$y)
    qq.gam(gam_xvar)
    if(stage==stages[2]){mtext(side=3,"B. Q-Q plot")}
    hist(resid,xlab="",main="")
    if(stage==stages[2]){mtext(side=3,"C. Residuals")}
    plot(linpred,resid)
    if(stage==stages[2]){mtext(side=3,"D. Residuals vs.",line=1);
      mtext(side=3,"linear pred.")}
    plot(fitted(gam_xvar),observed.y)
    if(stage==stages[2]){mtext(side=3,"E. Response vs.",line=1);
      mtext(side=3,"fitted")}
  }
  mtext(side=3,paste0("Model diagnostics: ",varnames[x]),outer=TRUE,line=2,font=2)
}


####Supplementary analyses: C, N, DNA and RNA as μg or % of body weight#####
vars<-colnames(img)[c(13,14,18,19,16,17,22,23)]
varnames<-c("C (μg)","N (μg)","DNA (μg)","RNA (μg)","%C","%N","%DNA","%RNA")

#Define predictor variables:
yvar<-list(
  "food",
  "predation",
  "s(tank,bs='re')",
  "s(day,k=4,by=food)",
  "s(day,k=4,by=predation)"
)

#Create table for R2, factor coefficients and P-values:
modtab_growth_supp<-array(NA,dim=c(length(stages)-1,5+length(yvar),length(vars)))
dimnames(modtab_growth_supp)<-list(stages[-1],c("R2","Food","Predation",yvar[1:3],
                                       "DayxLoF","DayxHiF","DayxLoP","DayxHiP"),varnames)

#Fit models and fill table for all response variables and stages:
for(x in 1:length(vars)){
  for(stage in stages){
    subdat<-img[img$stage==stage & !is.na(img[,vars[x]]),] #Subset data for stage
    #For table, set response variable to standardized within stage:
    subdat[,vars[x]]<-(subdat[,vars[x]]-mean(subdat[,vars[x]],na.rm=T))/
      sd(subdat[,vars[x]],na.rm=T)
    Xvar<-paste0(vars[x],"~") #Define response variable
    #Fit GAM if there are data available for 3 unique sampling days:
    if(length(unique(subdat$day))>3){
      gam_xvar<-gam(formula(paste(Xvar,paste(yvar,collapse="+"))),data=subdat)
      modtab_growth_supp[stage,1,varnames[x]]<-summary(gam_xvar)$r.sq #R2
      modtab_growth_supp[stage,2:3,varnames[x]]<-summary(gam_xvar)$p.table[2:3,1] #Factor term coefficients
      modtab_growth_supp[stage,unlist(yvar[1:2]),varnames[x]]<-summary(gam_xvar)$p.table[2:3,4] #P-values factor terms 
      modtab_growth_supp[stage,unlist(yvar[3]),varnames[x]]<-summary(gam_xvar)$s.table[1,4] #P-values simple smooth terms 
      modtab_growth_supp[stage,7:10,varnames[x]]<-summary(gam_xvar)$s.table[2:5,4] #P-values varying coefficient terms
    }
  }
}

#Create empty dataframe to store model predictions:
#'newdat' contains all combinations of day and treatment (tank is random):
newdat<-data.frame(day=rep(1:24,4),food=rep(levels(img$food),each=24*2),
                   predation=rep(rep(levels(img$predation),each=24),2),tank="G")

#Plot variation in observations + model predictions 
xvals<-c(0.5,5.5,11,16.5) #X-axis location per stage for boxplot

layout(matrix(c(1:32),nrow=8,ncol=4,byrow=TRUE),widths=c(0.55,rep(0.15,3)))
par(mar=c(0.25,0,0.1,0),oma=c(3,3.5,3,1))
for(x in 1:length(vars)){
  ylims<-range(img[,vars[x]],na.rm=T) #Limits of y-axis depends on variable
  #1: Make boxplot+stripchart per response variable:
  plot(1, type="n",xlim=c(0,20),ylim=ylims,axes=F,frame=T,ylab="",xlab="")
  #box(lwd=2)
  for(i in 1:length(stages)){
    xvals_i<-seq(xvals[i],xvals[i]+3,length.out=4) #Set positions along x-axis
    subdat<-img[img$stage==stages[i],] #Subset data for stage
    stripchart(subdat[,vars[x]]~subdat$treatment,add=T,vertical=T,pch=20,cex=0.4,col=adjustcolor(col_treatments$col),
               method="jitter",jitter=0.1,ylim=ylims,axes=F,frame=T,at=xvals_i,lwd=0)
    boxplot(subdat[,vars[x]]~subdat$treatment,outline=FALSE,add=T,col=adjustcolor(col_treatments$col,alpha.f=0.2),
            lwd=0.5,axes=F,frame=T,at=xvals_i)
  }
  mtext(side=2,varnames[x],outer=F,line=2)
  Axis(side=2)
  if(vars[x]==vars[length(vars)]){Axis(side=1,at=c(2,7,12.5,18),labels=stages);
    mtext(side=1,"Stages",line=2)}
  if(vars[x]==vars[1]){mtext(side=3,"A. Overall variation",outer=F,line=0);
    legend("topleft",legend=treatments,horiz=F,bty="n",cex=1.5,
           col=adjustcolor(col_treatments$col),pch=20,lty=1:4,x.intersp = 0.2,y.intersp = 0.75)}
  #2: Plot model predictions per stage (excluding C4):
  Xvar<-paste0(vars[x],"~") #Define response variable
  for(stage in stages[-1]){
    subdat<-img[img$stage==stage & !is.na(img[,vars[x]]),] #Subset data per stage
    xlims<-range(img$day[img$stage==stage]) #Limits of x-axis are stage-specific
    plot(1, type="n",xlim=xlims,ylim=ylims,axes=F,frame=T,ylab="",xlab="") #Plot background
    gam_xvar<-gam(formula(paste(Xvar, paste(yvar, collapse="+"))),data=subdat) #GAM fit
    pred<-predict.gam(gam_xvar,newdata=newdat,type="response",se.fit = T) #Model predictions
    newdat$pred<-pred$fit
    newdat$pred_hi<-pred$fit+1.96*pred$se.fit
    newdat$pred_lo<-pred$fit-1.96*pred$se.fit
    #Plot predictions for different treatments:
    for(i in 1:length(treatments)){
      daylim<-range(subdat$day[subdat$treatment==treatments[i]]) #Only plot predictions for days with observations in treatment
      sub_newdat<-newdat[newdat$day>=min(daylim) & newdat$day<=max(daylim) &
                           newdat$food==unique(img$food[img$treatment==treatments[i]]) &
                           newdat$predation==unique(img$predation[img$treatment==treatments[i]]),]  
      polygon(x=c(sub_newdat$day,rev(sub_newdat$day)),y=c(sub_newdat$pred_lo,rev(sub_newdat$pred_hi)),
              col=adjustcolor(col_treatments$col[col_treatments$treatments==treatments[i]],alpha.f=0.2),border=NA)
      points(sub_newdat$day,sub_newdat$pred,type="l",lty=i,lwd=1.5,col=adjustcolor(col_treatments$col[col_treatments$treatments==treatments[i]]))
    }
    if(vars[x]==vars[1]){mtext(side=3,stage,outer=F,line=0)}
    if(vars[x]==vars[length(vars)]){Axis(side=1,at=seq(1,24,by=6))}
  }
}
mtext(side=1,"Day",outer=T,line=2,adj=0.8)
mtext(side=3,"B. Temporal variation",outer=T,line=1.75,adj=0.85)
#Draw line on top of plot:
par(mfrow=c(1,1),xpd = NA,new=TRUE) 
plot(1:10,1:10, type = "n", axes = FALSE, xlab = "", ylab = "", main = "")
lines(c(5.9,10.45),c(11,11))

