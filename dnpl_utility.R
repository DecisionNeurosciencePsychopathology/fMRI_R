#######
##
##
devtools::source_url("https://raw.githubusercontent.com/DecisionNeurosciencePsychopathology/fMRI_R/master/fslpipe.R")
##Here's all the functions that helps with the fsl pipe function;

fsl_2_sys_env<-function(bashprofilepath=NULL){
  if (is.null(bashprofilepath)){bashprofilepath<-file.path(Sys.getenv("HOME"),".bash_profile")}
  if (length(system("env | grep 'fsl' ",intern = T))<1) {
  if(file.exists(bashprofilepath)) {
    print("Using user .bashprofile")
    fslinfo<-cfg_info(bashprofilepath)
    info_to_sysenv(fslinfo)
  }
  }
}

info_to_sysenv<-function(info=NULL) {
  if (is.null(info)) {stop("HEY! NO INFO!")}
  for (kz in names(info)) {
    print(kz)
    eval(parse(text=paste0("Sys.setenv(",kz,"='",info[[kz]],"')")))
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

#TO DO: Try to make recursive
recon.mat<-function(data.raw=NULL) {
  #do recurisive if class of x is list, if it's array then use 
  if (class(data.raw)=="list") {
    lapply(data.raw, function(x) {
      
    }) 
  } else if (class(data.raw)=="array") {
    data.raw<-recon.array(x=data.raw)
  } else if (!lapply(huy$muX,class) %in% c("list","array")) {}
    
  }


make_eprime_struct<-function(rawtext=file.choose(),
                           rawdf=NULL,
                           breakpointname="BreakProc") {
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

#End

#bandit specific functions:


#Make single signal as a list:
makesignal.single<-function(output,ename,norm="none",normargu=c("durmax_1","evtmax_1"),valuefrom,modifier=NA,style="subset",nonah=T) {
  jrz<-list()
  jrz[["event"]]<-ename
  jrz[["normalization"]]<-norm
  jrz[normargu]<-TRUE
  
  subelist<-output$event.list$allconcat[which(output$event.list$allconcat$event==ename),]
  value<-output$value[[valuefrom]]
  
  value.df<-data.frame(
    run=subelist$run,
    trial=subelist$trial,
    value=value
  )
  if (!is.na(modifier)) {
    #INVERTING THE LOGICAL STATEMENT SO 1 for valid trials and 0 for non-valid!
    if (nonah) {
    modx<-!output$output.df[modifier] & !is.na(subelist$onset)
    } else { modx<-!output$output.df[modifier]}
    switch (style,
            "multiply" = {value.df$value<-value.df$value * as.numeric(modx)},
            "subset"   = {value.df<-value.df[modx ,]}
    )
    
  }
  jrz[["value"]]<-value.df
  return(jrz)
}

#Make all of them with a grid:
make_signal_with_grid<-function(dsgrid=NULL,dsgridpath="grid.csv", outputdata=NULL,nona=F, expectedtn=300,add_taskness=F) {
  if (is.null(dsgrid)) {dsgrid<-read.csv(dsgridpath, stringsAsFactors = FALSE)}
  if (is.null(outputdata)) {stop("NO DATA")}
  dsgrid.og<-dsgrid
  allofthem<-new.env(parent = emptyenv())
  dsgrid[dsgrid=="NA"]<-NA
  if (length(grep("evt",dsgrid.og$valuefrom))>0){
  dsgrid<-dsgrid.og[-grep("evt",dsgrid.og$valuefrom),]} else {dsgrid.og->dsgrid}
  for (i in 1:length(dsgrid$ename)) {
    print(paste("making...",dsgrid$valuefrom[i],sep=""))
    #could have totally do a do.call and make assignment within single function 
    #but that ONE datainput variable is hard to get by with do.call so......loop is fine
    #Forcing all arguments so it's kinda bad....
    result<-makesignal.single(output = outputdata,
                              ename = dsgrid$ename[i],
                              norm = dsgrid$norm[i],
                              normargu = dsgrid$normargu[i],
                              valuefrom = dsgrid$valuefrom[i],
                              modifier = dsgrid$modifier[i],
                              style = dsgrid$style[i],
                              nonah = nona)
    assign(dsgrid$name[i],result,envir = allofthem)
  }
  
  #change it to list:
  allofthemlist<-as.list(allofthem)
  dsgrid.og->dsgrid
  #Taskness varibale:
  if (add_taskness) {
  for (taskname in unique(dsgrid$ename)) {
    
    if (any(is.na(outputdata$event.list[[taskname]]))) {
      tempdf<-outputdata$event.list[[taskname]][which(!is.na(outputdata$event.list[[taskname]]$onset)),]
      tempdf.a<-subset.data.frame(tempdf,select = c("run","trial"))
      tempdf.a$value<-1
      allofthemlist[[paste(taskname,"_evt",sep = "")]] <- list(value=tempdf.a,event=taskname,normalization="none")
    }else {allofthemlist[[paste(taskname,"_evt",sep = "")]] <- list(value=1,event=taskname,normalization="none")}
  }
}
  return(allofthemlist)
  
}

cfg_info<-function(cfgpath=NULL) {
  if (is.null(cfgpath)) {stop("No cfg file supplied!")}
  pre.sym<-system("env",intern = T)
  sysm.temp<-system(paste("source",cfgpath,"\n","env"),intern = T)
  sysm<-sysm.temp[which(!sysm.temp %in% pre.sym)]
  larg<-regmatches(sysm, regexpr("=", sysm), invert = TRUE)
  xout<-as.environment(list())
  NULL.x<-lapply(larg,function(y) {
    if (length(y)>1) {
      if (length(grep("-* ",y[2]))>0) {
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
    }
    assign(y[1],x.2x,envir = xout)
  })
  return(as.list(xout))
}

get_nuisance_preproc<-function(id=NULL,
                              cfgfilepath="/Volumes/bek/autopreprocessing_pipeline/Learn/bandit_oldPreCMMR.cfg",
                              returnas=c("path","data.frame")
) {
  cfg<-cfg_info(cfgpath = cfgfilepath)
  lpath<-lapply(1:cfg$n_expected_funcruns, function(i) {
    list(
    nuisance=
    file.path(cfg$loc_mrproc_root,id,cfg$preprocessed_dirname,paste(cfg$paradigm_name,i,sep = ""),cfg$preproc_call$nuisance_file),
    motion=
    file.path(cfg$loc_mrproc_root,id,cfg$preprocessed_dirname,paste(cfg$paradigm_name,i,sep = ""),"motion.par"))
    })
  names(lpath)<-paste("run",1:cfg$n_expected_funcruns,sep = "")
  if (returnas=="path") {
  return(lpath)} else if (returnas=="data.frame") {
    ldf<-lapply(lpath,function(x) {
      nui<-read.table(x$nuisance)
      names(nui)<-unlist(strsplit(cfg$preproc_call$nuisance_compute,split = ","))
      mo<-read.table(x$motion)
      names(mo)<-paste0("motion_V",1:length(mo))
      combo<-cbind(nui,mo)
      return(combo)
    })
    
  }
}

####################
get_volume_run<-function(id=NULL,
                         cfgfilepath=NULL,
                         reg.nii.name="swudktm*[0-9].nii.gz",
                         returnas=c("path","numbers")){
  cfg<-cfg_info(cfgpath = cfgfilepath)
  if (returnas=="path"){
  lpath<-lapply(1:cfg$n_expected_funcruns, function(i) {
    file.path(cfg$loc_mrproc_root,id,cfg$preprocessed_dirname,paste(cfg$paradigm_name,i,sep = ""))->procpath
    system(paste0("find ",procpath," -iname ",reg.nii.name," -maxdepth 2 -mindepth 1"),intern = T)
   # file.path(,nii.name)
  })
  return(unlist(lpath))
  }
  if (returnas=="numbers"){
    lnum<-lapply(1:cfg$n_expected_funcruns, function(i) {
      length(readLines(file.path(cfg$loc_mrproc_root,id,
                cfg$preprocessed_dirname,
                paste(cfg$paradigm_name,i,sep = ""),
                "motion_info","fd.txt")
      ))
      # file.path(,nii.name)
    })
    return(unlist(lnum))
  }
}

make_heatmap_with_design<-function(design=NULL) {
  return(dependlab::cor_heatmap(as.data.frame(dependlab::concat_design_runs(design))))
}

findbox<-function() {
  if (Sys.getenv("USER")=="jiazhouchen") {boxdir <- "/Users/jiazhouchen/Box Sync"
  } else if (Sys.getenv("USER")=="jiazhou") {boxdir <- "/Volumes/bek/Box Sync"} else {
    boxdir<-system("find ~ -iname 'Box*' -maxdepth 2 -type d",intern = T)}
  return(boxdir)
}

######General function for Single subject loop: (can be ready for lapply or do call)
do.all.subjs<-function(
  tid=NULL,
  do.prep.call="prep.son1",
  do.prep.arg=list(son1_single=son1_single),
  cfgpath=NULL,
  regpath=NULL,
  gridpath="grid.csv",
  func.nii.name="swudktm*[0-9].nii.gz",
  proc_id_subs=NULL,    #Put "" for nothing.
  wrt.timing=c("convolved", "FSL","AFNI"),
  model.name=NULL,
  model.varinames=NULL,
  add.nuisa=TRUE,
  assigntoenvir=NULL) {
  
  #Read config file:
  cfg<-cfg_info(cfgpath)
  
  #Prep the data into generally acceptable output object;
  output<-do.call(what = do.prep.call,args = do.prep.arg)
  #assign("output",do.call(what = do.prep.call,args = do.prep.arg),envir=globalenv())
  
  dsgrid<-read.csv(gridpath,stringsAsFactors = F)
  #if (length(grep("evt",dsgrid.og$valuefrom))>0){
  #  dsgrid<-dsgrid.og[-grep("evt",dsgrid.og$valuefrom),]} else {dsgrid.og->dsgrid}
  #Generate signal with make signal with grid function (grid.csv need to be in working directory or specified otherwise)
  signals<-make_signal_with_grid(outputdata = output,add_taskness = T,dsgrid = dsgrid)
  
  if (length(grep("evt",dsgrid$valuefrom))>0){
    dxgrid<-dsgrid[grep("evt",dsgrid$valuefrom),]
    for (u in 1:length(dxgrid$name)) {
      signals[dxgrid$name[u]]<-signals[dxgrid$valuefrom[u]]
    }
  }  
  #Get nuissance regressor: 
  #Still concat nuisa regressor together
  nuisa<-get_nuisance_preproc(id=paste0(tid,proc_id_subs),cfgfilepath = cfgpath,returnas = "data.frame") 
  if (add.nuisa) {
    nuisa.x<-nuisa
  } else {
    nuisa.x<-NULL
  }
  
  #Get the actual volume by run:
  run_volum<-get_volume_run(id=paste0(tid,proc_id_subs),
                            cfgfilepath = cfgpath,
                            returnas = "numbers",
                            reg.nii.name = func.nii.name)
  
  #Create  models:
  model<-signals[model.varinames]
  
  #Use Michael's package to generate design matrix and correlation graph;
  design<-dependlab::build_design_matrix(
    events = output$event.list$allconcat, #Load the task info
    signals = model,     #Load the Model
    write_timing_files = wrt.timing, #Output timing files to FSL style
    tr=as.numeric(cfg$preproc_call$tr), #Grab the tr from cfg instead of hard coding it...
    plot = F,
    run_volumes = run_volum,
    #tr=1 second, maybe need to double check, I'm kinda sure....
    output_directory = file.path(regpath,model.name,tid), #Where to output the timing files, default is the working directory
    nuisance_regressors = nuisa.x #Maybe could add in nuisance_regressors from pre-proc
  )
  
  design$heatmap<-make_heatmap_with_design(design)
  design$volume<-run_volum
  design$nuisan<-nuisa
  design$ID<-tid
  design$preprocID<-paste0(tid,proc_id_subs)
  design$regpath<-file.path(regpath,model.name,tid)
  
  if (!is.null(nuisa)){
    for (k in 1:length(nuisa)) {
      write.table(as.matrix(nuisa[[k]]),file.path(regpath,model.name,tid,
                                                  paste0("run",k,"_nuisance_regressor_with_motion.txt")),
                  row.names = F,col.names = FALSE)
    }}
  
  if (is.environment(assigntoenvir)) {assign(as.character(tid),design,envir = assigntoenvir)
  } else {return(design)}
  
}
######Modify fsl template with variable switch
change_fsl_template<-function(fsltemplate=NULL,begin="ARG_",end="_END",searchenvir=xarg) {
  for (tofind in grep(paste0(begin,"*"),fsltemplate) ){
    tryCatch(
      {
      varixma<-substr(fsltemplate[tofind],regexpr(paste0(begin,"*"),fsltemplate[tofind])+nchar(begin),
                       regexpr(paste0("*",end),fsltemplate[tofind])-1)
      fsltemplate[tofind]<-gsub(paste0(begin,varixma,end),searchenvir[[varixma]],fsltemplate[tofind])
      },
      error=function(x){print("something went wrong...")})
  }
  return(fsltemplate)
}

#####Generate reg path from model name:

gen_reg<-function(vmodel=NULL,regpath=NULL,idx=NULL,runnum=NULL,env=NULL,regtype=NULL) {
  NUP<-lapply(vmodel<-argu$model.varinames, function(x) {
    assign(paste0(x,"reg"),file.path(regpath,idx,paste0("run",runnum,"_",x,regtype)),envir = env)
  })
}               
#Use data.frame"
 #Name, Begining, Ending, alter = name, beg, dura

#directory to nuisances
# 
# con <- dbConnect(odbc(),
#                  Driver = "SQLServer",
#                  Server = "OACDATA4",
#                  Database = "crc",
#                  UID = "chenj20",
#                  PWD = rstudioapi::askForPassword("Database password"),
#                  Port = 1433)

feat_w_template<-function(templatepath=NULL,
                          fsltemplate=NULL,
                          beg="ARG_",
                          end="_END",
                          fsf.path=NULL,
                          envir=NULL) {
  if (is.null(fsltemplate)) {fsltemplate<-readLines(templatepath)}
  subbyrunfeat<-change_fsl_template(fsltemplate = fsltemplate,begin = beg,end=end,searchenvir = envir)
  fsfpath<-fsf.path
  writeLines(subbyrunfeat,fsfpath)
  message("starting to do feat...")
  system(paste0("feat ",fsfpath),intern = T)
  message("feat completed")
}

plot_image_all<-function(rootpath=NULL,
                         templatedir=NULL,
                         model.name=NULL,
                         patt=NULL,
                         threshold=0.99,
                         outputdir=getwd(),
                         colour="red") {
  dirs<-system(paste0("find ",file.path(rootpath,model.name)," -iname '",patt,"' -maxdepth 3 -mindepth 1 -type f"),intern = T)
  for (sdir in dirs) {
    spdir<-strsplit(sdir,.Platform$file.sep) 
    spdir[[1]][sapply(spdir, function(x) {grep(patt,x)})-1]->filex
    tep<-readNIfTI(templatedir)
    img<-readNIfTI(sdir)
    mask<-tep
    in_mask<- img > threshold
    mask[in_mask] <- 1
    mask[!in_mask] <- NA
    jpeg(filename = file.path(outputdir,paste0(filex,".jpeg")),width = 2560, height = 1440,quality = 90,pointsize = 20)
    orthographic(x = tep, y = mask, col.y = colour)
    dev.off()
  }
}






