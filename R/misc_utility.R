###miscellaneous utility functions:
###This script will now house the misc. functions that some major functions depends on
##Mostly helper functions
cleanuplist<-function(listx){
  if (any(sapply(listx, is.null))){
    listx[sapply(listx, is.null)] <- NULL}
  return(listx)
}

fsl_2_sys_env<-function(bashprofilepath=NULL,force=T){
  if (is.null(bashprofilepath)){bashprofilepath<-file.path(Sys.getenv("HOME"),".bash_profile")}
  if (length(system("env | grep 'FSL' ",intern = T))<1 | force) {
    if(file.exists(bashprofilepath)) {
      print("Using user .bashprofile")
      fslinfo<-cfg_info(bashprofilepath)
      pathx<-unique(strsplit(fslinfo$PATH,":")[[1]])
      afni_path<-pathx[which(grepl("/abin",pathx) | grepl("/afni",pathx))]
      if(length(afni_path)==1) {fslinfo$AFNIDIR <- afni_path}
      info_to_sysenv(fslinfo)
    }
  }
}

info_to_sysenv<-function(info=NULL) {
  if (is.null(info)) {stop("HEY! NO INFO!")}
  for (kz in names(info)) {
    print(kz)
    tryCatch({
      eval(parse(text=paste0("Sys.setenv(",kz,"","='",info[[kz]],"')")))},error=function(e){message(e)}
    )
  }  
}

recondrop<-function(datax){
  if (is.list(datax)){
    lapply(datax, function(x) {
      recondrop(x)
    })
  } else if (is.null(dim(datax))) {return(datax)
  } else if (length(dim(datax))>2) {
    drop(datax)->datay
    recondrop(datay)
  } 
}

recon.list<-function(mat.raw=NULL){
  mat.better<-sapply(mat.raw, function(x){
    emp<-list()
    drop(x)->x
    for (i in 1:length(x)) {
      rownames(x)[i]->whichone
      emp[[whichone]]<-drop(x[i])
    }
    return(emp)
  })
  names(mat.better)<-names(mat.raw)
  return(mat.better)
}

recon.array<-function(x=NULL) {
  emp<-list()
  for (i in 1:length(x)) {
    rownames(x)[i]->whichone
    if (length(emp[[whichone]][[1]])==1) {
      emp[[whichone]]<-drop(unlist((x[[i]])))
    } else if (ifelse(is.null(dim(x[[i]])[1]),-1,dim(x[[i]])[1]) == 1) {
      emp[[whichone]]<-drop(unlist(x[[i]]))
    } else {emp[[whichone]]<-drop((x[[i]]))}
  }
  #names(emp)<-names(x)
  return(emp)
}

recon.mat<-function(data.raw=NULL) {
  #do recurisive if class of x is list, if it's array then use 
  if (class(data.raw)=="list") {
    lapply(data.raw, function(x) {
      
    }) 
  } else if (class(data.raw)=="array") {
    data.raw<-recon.array(x=data.raw)
  } else if (!lapply(huy$muX,class) %in% c("list","array")) {}
  
}

make_eprime_struct<-function(rawtext=file.choose(),rawdf=NULL,breakpointname="BreakProc") {
  if (is.null(rawdf)){
    rawdf<-read.csv(rawtext,stringsAsFactors = F,sep = "\t")
  } 
  wherearethebreaks<-grep(breakpointname,rawdf$Procedure)
  rawdf$blocknum<-NA
  
  for (i in 1:(length(wherearethebreaks)+1)) {
    print(i)
    if (i==1) {
      rawdf$blocknum[i:wherearethebreaks[i]]<-i
      
    } else if (i==(length(wherearethebreaks)+1)) {
      rawdf$blocknum[wherearethebreaks[i-1]:length(rawdf$Procedure)]<-i
      
    } else {
      rawdf$blocknum[wherearethebreaks[i-1]:wherearethebreaks[i]]<-i}
  }
  df<-rawdf[-wherearethebreaks,]
  df$trailnum<-seq_along(df$Procedure)
  df$trialnum.block<-unsplit(lapply(split(df$blocknum,df$blocknum), seq_along),df$blocknum)
  return(df)
}

sigmatransform<-function(x) {
  (1./(1+exp(x)))->y
  return(y)
}

depthoflist <- function(list,thisdepth=0){
  if(!is.list(list)){
    return(thisdepth)
  }else{
    return(max(unlist(lapply(list,depthoflist,thisdepth=thisdepth+1))))    
  }
}

addcenterscaletolist<-function(list) {
  test<-lapply(list, scale,center=T)
  names(test)<-paste(names(list),"centerscaled",sep = "_")
  newlist<-c(list,test) 
  return(newlist)
}

cfg_info<-function(cfgpath=NULL,noproc=F) {
  if (is.null(cfgpath)) {stop("No cfg file supplied!")}
  pre.sym<-system("env",intern = T)
  sysm.temp<-system(paste("source",cfgpath,"\n","env"),intern = T)
  sysm<-sysm.temp[which(!sysm.temp %in% pre.sym)]
  sysm<-sysm[!grepl("()",sysm,fixed = T)]
  larg<-regmatches(sysm, regexpr("=", sysm), invert = TRUE)
  xout<-as.environment(list())
  NULL.x<-lapply(larg,function(y) {
    if (length(y)>1) {
      if (length(grep("-* ",y[2]))>0 & !noproc) {
        tc<-strsplit(y[2],split = " -")[[1]]
        if (regexpr("-",tc[1])[[1]] > 0) {
          tc[1]<-substr(tc[1],start = 2,stop = nchar(tc[1]))}
        tx<-strsplit(tc,split = " ")
        etx<-as.environment(list())
        NUx<-lapply(tx, function(x) {
          assign(x[1],x[2],envir = etx)
        })
        x.2x<-as.list(etx)
      } else {x.2x<-y[2]} 
      assign(y[1],x.2x,envir = xout)
    } else {return(NULL)}
  })
  return(as.list(xout))
}

findbox<-function(usebek=F) {
  if(usebek){boxdir<-"/Volumes/bek/Box Sync"} else if (dir.exists("~/Box")) {
    boxdir <- "~/Box"
  } else {
    boxdir<-system("find ~ -iname 'Box*' -maxdepth 2 -type d",intern = T)}
  return(boxdir)
}

getMotion_report<-function(infilepath=file.choose(),dvar_thresh=20,fd_thresh=0.9,filter_name=NULL,protocol=ptcs$masterdemo){
  masterdemo<-bsrc::bsrc.checkdatabase2(protocol = protocol,batch_size=1000L,forceskip = T,online = T)
  dfa<-read.csv(infilepath,stringsAsFactors = F)
  dfa$X<-NULL
  
  if(is.null(dfa$session)){
    dfa$session <- dfa$run
    dfa$modelname <- dfa$func_name
  } 
  if(!is.null(filter_name)){
    dfa <- dfa[which(dfa$modelname %in% filter_name),]
  }
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
                                                                "registration_group","registration_lethality")],id.var = "ID",onlyoutput = T)
  outdf$ogid<-NULL;outdf$ifexist<-NULL;outdf$registration_wpicid<-NULL;outdf$registration_soloffid<-NULL;
  names(outdf)[which(names(outdf)=="registration_redcapid")]<-"DNPL_ID"
  names(outdf)[which(names(outdf)=="registration_group")]<-"Group"
  names(outdf)[which(names(outdf)=="registration_lethality")]<-"Lethality"
  
  
  sp_b<-split(outdf,outdf$ID)
  outdf2<-do.call(rbind,lapply(sp_b,function(dfc){
    ga<-cbind(as.data.frame(as.list(apply(dfc[1:2],2,mean))),dfc[1,3:ncol(dfc)])
    ga$session<-"mean"
    gat<-reshape2::melt(data = rbind(dfc,ga),id.vars=c("ID","DNPL_ID","Group","Lethality","modelname","session"))
    gat$vari<-paste(gat$variable,gat$session,sep = "_")
    reshape2::dcast(data = gat,formula = ID+modelname+DNPL_ID+Group+Lethality~vari,value.var = "value")
  }))
  outdf2$Group[which(outdf2$Group=="")]<-NA
  return(outdf2)
}

get_matrix<-function(raw_text,heading="/Matrix",ending=NULL,colnames=NULL,split=" "){
  if(is.null(ending)){end_pos<-length(raw_text)}else{end_pos<-(grep(ending,raw_text)-1)}
  if(is.null(heading)){end_pos<-0}else{start_pos<-(grep(heading,raw_text)+1)}
  raw_mx<-raw_text[start_pos:end_pos]
  matrix_df<-as.data.frame(do.call(rbind,strsplit(raw_mx,split = split)))
  if(!is.null(colnames) && length(colnames)==ncol(matrix_df)){names(matrix_df)<-colnames} 
  return(matrix_df)
}

check_argu<-function(argu=NULL,init=F){
  message("Incomplete function!")
  requiresargus<-c("path.cfg","funcimg.namepatt","path.outroot","path.lvl1grid","ssub.func","ssub.datalist","modelname")
  if(is.null(argu)) {argu<-list()}
  if(!init & any(!requiresargus %in% names(argu)) ) {
    stop("These arguments are required and not present: ",paste(requiresargus[!requiresargus %in% names(argu)],collapse = ", "))
  }
  defaultargus<-list(ncpus=4,runstep=NULL,forcererun=F,proc_id_subs="",
                     #Level 1 arguements:
                     lvl1.CenterScaleAll=F,lvl1.MotionMethod="",lvl1.ZThresh=1.96,lvl.PThresh=0.05,
                     
                     #Level 2 arguments:
                     
                     #Level 3 argumetns:
                     lvl3.Test = "FLAME1+2",lvl3.AutoContrast=T,lvl3.SepGroupVari=F
  )
  
  for(xi in names(defaultargus)){
    if(is.null(argu[[xi]])){
      argu[[xi]]<-defaultargus[[xi]]
    }
  }
  
  
  if(!init){
    #Path config:
    argu$path.regroot<-file.path(argu$path.outroot,"regressors",argu$modelname)
    argu$path.ssubroot<-file.path(argu$path.outroot,"single_subject",argu$modelname)
    argu$path.grproot<-file.path(argu$path.outroot,"group_analysis",argu$modelname)
    #######lvl 1 motion sensoring##########
    if(argu$lvl1.MotionMethod=="spike_regression") {
      message("First level motion method is selected to be spike_regression.")
      if(is.null(argu$lvl1.SpikeRegType)){
        message("Motion spike regressor's type is not provided, using default 'fd'")
        argu$lvl1.SpikeRegType<-"fd"
      }
      if(argu$lvl1.SpikeRegType=="fd" & is.null(argu$SpikeRegThres)){
        message("Motion spike regressor's (type FD) threshold is not provided, using default 0.9")
        argu$SpikeRegThres<-0.9
      }
      
      if(argu$lvl1.SpikeRegType=="dvar" & is.null(argu$SpikeRegThres)){
        message("Motion spike regressor's (type DVAR) threshold is not provided, using default 20")
        argu$SpikeRegThres<-20
      }
    }
    
    argu$cfg<-cfg_info(cfgpath = argu$path.cfg)
    argu$lvl1grid<-read.table(argu$path.lvl1grid,header = T,sep = c(","),stringsAsFactors = F,strip.white = T,skipNul = T)
    ##############END OF LVL1################
  }
  #lvl 2
  
  #lvl 3
  
  
  
  
  
  
  # 
  # #DO NOT PROVIDE THOSE TWO AND IT WILL BE FINE;
  # argu$randomize_p_threshold<-0.001
  # argu$randomize_thresholdingways<-c("tfce","voxel-based","cluster-based-extent","cluster-based-mass")
  # argu$ss_zthreshold<-3.2  #This controls the single subject z threshold (if enabled in template)
  # argu$ss_pthreshold<-0.05 #This controls the single subject p threshold (if enabled in template)
  
  
  return(argu)
}







