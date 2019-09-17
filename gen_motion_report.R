infilepath<-file.choose()
dvar_thresh<-20
fd_thresh<-0.9


getMotion_report<-function(infilepath=file.choose(),dvar_thresh=20,fd_thresh=0.9){
masterdemo<-bsrc::bsrc.checkdatabase2(protocol = ptcs$masterdemo,batch_size=1000L)
dfa<-read.csv(infilepath,stringsAsFactors = F)
dfa$X<-NULL
sp_a<-split(dfa,paste(dfa$ID,dfa$session,sep = "_"))
outdf<-do.call(rbind,lapply(sp_a,function(dfb){
  data.frame(per_fd=length(which(as.numeric(dfb$fd) >= fd_thresh))/nrow(dfb),
             per_dvar=length(which(as.numeric(dfb$dvars) <= dvar_thresh))/nrow(dfb),
             max_fd = max(as.numeric(dfb$fd)),
             max_dvar = max(as.numeric(dfb$dvar)),
             modelname=unique(dfb$modelname),ID=unique(dfb$ID),session=unique(dfb$session),stringsAsFactors = F)

}))
outdf<-outdf[outdf$max_fd <= 10,]
outdf<-bsrc::bsrc.findid(df = outdf,idmap = masterdemo$data[c("registration_redcapid","registration_wpicid","registration_soloffid",
                                                              "registration_group","registration_lethality")],id.var = "ID")
outdf$Group<-outdf$registration_group
outdf$Lethality<-outdf$registration_lethality
sp_b<-split(outdf,outdf$ID)
outdf2<-do.call(rbind,lapply(sp_b,function(dfc){
  ga<-cbind(as.data.frame(as.list(apply(dfc[1:2],2,mean))),dfc[1,3:ncol(dfc)])
  ga$session<-"mean"
  gat<-reshape2::melt(data = rbind(dfc,ga),id.vars=c("ID","Group","Lethality","modelname","session"))
  gat$vari<-paste(gat$variable,gat$session,sep = "_")
  reshape2::dcast(data = gat,formula = ID+modelname+Group+Lethality~vari,value.var = "value")
}))
outdf2$Group[which(outdf2$Group=="")]<-NA
return(outdf2)
}


