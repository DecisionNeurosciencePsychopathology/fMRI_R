#######
##
##
#devtools::source_url("https://raw.githubusercontent.com/DecisionNeurosciencePsychopathology/fMRI_R/master/fslpipe.R")
#devtools::source_url("https://raw.githubusercontent.com/DecisionNeurosciencePsychopathology/fMRI_R/master/prep_for_second_lvl.R")
##Here's all the functions that helps with the fsl pipe function;

cleanuplist<-function(listx){
  if (any(sapply(listx, is.null))){
    listx[sapply(listx, is.null)] <- NULL}
  return(listx)
}

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

#Make single signal as a list:
makesignal.single<-function(output,ename,norm=c("none","durmax_1","evtmax_1"),normargu="convmax_1",valuefrom,modifier=NA,style="subset",nonah=T) {
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
  #0->value.df$value[is.na(value.df$value)]
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
  if (length(grep("_evt",dsgrid.og$valuefrom))>0){
    dsgrid<-dsgrid.og[-grep("_evt",dsgrid.og$valuefrom),]} else {dsgrid.og->dsgrid}
  
  if (dim(dsgrid)[1]>0) {
    for (i in 1:length(dsgrid$ename)) {
      message(paste("making...",dsgrid$name[i],sep=""))
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
      # output = outputdata
      # ename = dsgrid$ename[i]
      # norm = dsgrid$norm[i]
      # normargu = dsgrid$normargu[i]
      # valuefrom = dsgrid$valuefrom[i]
      # modifier = dsgrid$modifier[i]
      # style = dsgrid$style[i]
      # nonah = nona
      assign(dsgrid$name[i],result,envir = allofthem)
    }
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

cfg_info<-function(cfgpath=NULL,noproc=F) {
  if (is.null(cfgpath)) {stop("No cfg file supplied!")}
  pre.sym<-system("env",intern = T)
  sysm.temp<-system(paste("source",cfgpath,"\n","env"),intern = T)
  sysm<-sysm.temp[which(!sysm.temp %in% pre.sym)]
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
    }
    assign(y[1],x.2x,envir = xout)
  })
  return(as.list(xout))
}

get_nuisance_preproc<-function(id=NULL,cfgfilepath="/Volumes/bek/autopreprocessing_pipeline/Learn/bandit_oldPreCMMR.cfg",
                               returnas=c("path","data.frame"),dothese=c("nuisance","motion_par","motion_outlier"),
                               type="fd",threshold="default") {
  if (type == "dvar") {
    moutn<-"dvars.txt"
    if (threshold == "default") {threshold<-20}
  }
  if (type == "fd") {
    moutn<-"fd.txt"
    if (threshold == "default") {threshold<-0.9}
  }
  
  cfg<-cfg_info(cfgpath = cfgfilepath)
  lpath<-lapply(1:cfg$n_expected_funcruns, function(i) {
    list(
      nuisance=
        file.path(cfg$loc_mrproc_root,id,cfg$preprocessed_dirname,paste(cfg$paradigm_name,i,sep = ""),cfg$preproc_call$nuisance_file),
      motionpar=
        file.path(cfg$loc_mrproc_root,id,cfg$preprocessed_dirname,paste(cfg$paradigm_name,i,sep = ""),"motion.par"),
      motionoutlier=
        file.path(cfg$loc_mrproc_root,id,cfg$preprocessed_dirname,paste(cfg$paradigm_name,i,sep = ""),"motion_info",moutn))
  })
  names(lpath)<-paste("run",1:cfg$n_expected_funcruns,sep = "")
  if (returnas=="path") {
    return(lpath)} else if (returnas=="data.frame") {
      ldf<-lapply(lpath,function(x) {
        combo<-list()
        if ("nuisance" %in% dothese) {
          nui<-read.table(x$nuisance)
          names(nui)<-unlist(strsplit(cfg$preproc_call$nuisance_compute,split = ","))
        } else {nui<-data.frame()}
        if ("motion_par" %in% dothese) {
          mopar<-read.table(x$motionpar)
          names(mopar)<-paste0("motion_V",1:length(mopar))
        } else {mopar<-data.frame()}
        if ("motion_outlier" %in% dothese) {
          moout<-read.table(x$motionoutlier)
          moout$logi<-moout$V1>threshold
          if (any(moout$logi)) {
            xj<-matrix(data = 0,nrow = length(moout$logi),ncol = length(which(moout$logi)))
            for (xn in 1:length(which(moout$logi))) {
              xj[which(moout$logi)[xn],xn]<-1
            }
            moout<-as.data.frame(xj)
            names(moout)<-paste0("motionoutlier_",1:length(moout))
          } else {moout<-data.frame()}
        } else {moout<-data.frame()}
        
        combo<-c(nui,mopar,moout)
        combo<-as.data.frame(combo)
        return(combo)
      })
      
    }
}

get_volume_run<-function(id=NULL,cfgfilepath=NULL,reg.nii.name="swudktm*[0-9].nii.gz",returnas=c("path","numbers")){
  cfg<-cfg_info(cfgpath = cfgfilepath)
  if (returnas=="path"){
    lpath<-lapply(1:cfg$n_expected_funcruns, function(i) {
      file.path(cfg$loc_mrproc_root,id,cfg$preprocessed_dirname,paste(cfg$paradigm_name,i,sep = ""))->procpath
      system(paste0("find -L ",procpath," -iname ",reg.nii.name," -maxdepth 2 -mindepth 1"),intern = T)
      # file.path(,nii.name)
    })
    return(unlist(lpath))
  }
  if (returnas=="numbers"){
    lnum<-lapply(1:cfg$n_expected_funcruns, function(i) {
      fdpath<-file.path(cfg$loc_mrproc_root,id,
                        cfg$preprocessed_dirname,
                        paste(cfg$paradigm_name,i,sep = ""),
                        "motion_info","fd.txt")
      if (file.exists(fdpath)) {
        length(readLines(fdpath))
      } else return(NA)
      
    })
    return(unlist(lnum))
  }
}

make_heatmap_with_design<-function(design=NULL) {
  return(dependlab::cor_heatmap(as.data.frame(dependlab::concat_design_runs(design))))
}

findbox<-function(usebek=F) {
  if(usebek){boxdir<-"/Volumes/bek/Box Sync"} else if (dir.exists("~/Box")) {
    boxdir <- "~/Box"
  } else {
    boxdir<-system("find ~ -iname 'Box*' -maxdepth 2 -type d",intern = T)}
  return(boxdir)
}
prepare4secondlvl<-function(ssana.path=NULL,preproc.path=NULL,
                            standardbarin.path="/Volumes/bek/Newtemplate_may18/fsl_mni152/MNI152_T1_2mm_brain.nii",
                            proc.name=NULL,taskname=NULL,dir.filter=NULL,overwrite=TRUE,outputmap=FALSE,paralleln=NULL) {
  if (is.null(ssana.path) | is.null(preproc.path) | is.null(standardbarin.path) | is.null(proc.name) | is.null(taskname)){stop("not enough info to run")}
  #Step 1: Step up parameters, but it's a function so do it outside please
  #Step 2: Go to find all the feat. directory:
  featlist<-system(paste0("find ",ssana.path," -iname '*output.feat' -maxdepth 4 -mindepth 1 -type d"),intern = T)
  #Break them down&take out all old_template_results
  strsplit(featlist,split = "/")->s.featlist
  s.featlist.p<-lapply(s.featlist,function(x) {
    if (any(x %in% dir.filter)) {
      x<-NULL
    } else {x}
  })
  s.featlist.p[sapply(s.featlist.p,is.null)] <- NULL
  maxlength<-length(s.featlist.p[[1]])
  #Step 3: Get all the necessary info from the breakdown list:
  linkmap<-data.frame(id=sapply(s.featlist.p, "[[",maxlength-1), runword=sapply(s.featlist.p, "[[",maxlength))
  linkmap$num<-lapply(strsplit(as.character(linkmap$runword),split =''), 
                      function(x) {suppressWarnings(x[which(!is.na(as.numeric(x)))])})
  linkmap$num<-as.numeric(linkmap$num)
  na.omit(linkmap)->linkmap 
  linkmap$origin<-paste(taskname,linkmap$num,sep = "")
  
  #Step 4, make original directory and destination directory:
  linkmap$destination<-file.path(ssana.path,linkmap$id,linkmap$runword,"reg")
  linkmap$destination_standard<-file.path(ssana.path,linkmap$id,linkmap$runword,"reg_standard")
  linkmap$originplace<-file.path(ssana.path,linkmap$id,linkmap$runword,"masktostandtransforms.mat")
  
  if (is.null(paralleln)) {
    for (i in 1:length(linkmap$id)) {
      if (overwrite) {
        if (file.exists(linkmap$destination[i])) {file.remove(linkmap$destination[i])}
        if (file.exists(linkmap$originplace[i])) {file.remove(linkmap$originplace[i])}
      }
      st2<-dir.create(showWarnings = F,path = linkmap$destination[i])
      if (!file.exists(linkmap$originplace[i]))  {
        fsl_2_sys_env()
        system(paste("${FSLDIR}/bin/flirt",
                     "-in",file.path(ssana.path,linkmap$id[i],linkmap$runword[i],"mask.nii.gz"),
                     "-ref",standardbarin.path,
                     "-omat",linkmap$originplace[i],
                     "-usesqform",sep = " "),intern = F)
      }
      if (!file.exists(file.path(linkmap$destination,"example_func2standard.mat")[i])){
        file.symlink(from = file.path(linkmap$originplace)[i],to = file.path(linkmap$destination,"example_func2standard.mat")[i])
        file.symlink(from = file.path(linkmap$originplace)[i],to = file.path(linkmap$destination,"standard2example_func.mat")[i])
        file.symlink(from = standardbarin.path,to = file.path(linkmap$destination,"standard.nii.gz")[i])
      } else {message("meh,already there, if you want to overwirite, do overwrite...")}
      
    }} else {
      cluster_prep2ndlvl<-makeCluster(paralleln,outfile=file.path(argu$ssub_outputroot,argu$model.name,"prep42ndlvl_log.txt"),type = "FORK")
      NU<-parSapply(cluster_prep2ndlvl,1:length(linkmap$id),function(i) {
        fsl_2_sys_env()
        if (overwrite) {
          if (file.exists(linkmap$destination[i])) {file.remove(linkmap$destination[i])}
          if (file.exists(linkmap$originplace[i])) {file.remove(linkmap$originplace[i])}
        }
        st2<-dir.create(showWarnings = F,path = linkmap$destination[i])
        if (!file.exists(linkmap$originplace[i]))  {
          
          system(paste("${FSLDIR}/bin/flirt",
                       "-in",file.path(ssana.path,linkmap$id[i],linkmap$runword[i],"mask.nii.gz"),
                       "-ref",standardbarin.path,
                       "-omat",linkmap$originplace[i],
                       "-usesqform",sep = " "),intern = F)
        }
        if (!file.exists(file.path(linkmap$destination,"example_func2standard.mat")[i])){
          file.symlink(from = file.path(linkmap$originplace)[i],to = file.path(linkmap$destination,"example_func2standard.mat")[i])
          file.symlink(from = file.path(linkmap$originplace)[i],to = file.path(linkmap$destination,"standard2example_func.mat")[i])
          file.symlink(from = standardbarin.path,to = file.path(linkmap$destination,"standard.nii.gz")[i])
        } else {message("meh,already there, if you want to overwirite, do overwrite...")}
        
      })
      stopCluster(cluster_prep2ndlvl)
    }
  
  if(outputmap) {return(linkmap)}
  print("DONE")
}

######General function for Single subject loop: (can be ready for lapply or do call)
do.all.subjs<-function(tid=NULL,do.prep.call="prep.son1",do.prep.arg=list(son1_single=son1_single),cfgpath=NULL,
                       regpath=NULL,gridpath="grid.csv",func.nii.name="swudktm*[0-9].nii.gz",proc_id_subs=NULL, 
                       wrt.timing=c("convolved", "FSL","AFNI"),model.name=NULL,model.varinames=NULL,
                       nuisa_motion=c("nuisance","motion_par","motion_outlier"),motion_type="fd",centerscaleall=FALSE,
                       motion_threshold="default",convlv_nuisa=F,argu=NULL) {
  
  #Read config file:
  cfg<-cfg_info(cfgpath)
  argu$cfg<-cfg
  #Prep the data into generally acceptable output object;
  output<-do.call(what = do.prep.call,args = do.prep.arg,envir = argu)
  #assign("output",do.call(what = do.prep.call,args = do.prep.arg),envir=globalenv())
  dsgrid<-read.table(gridpath,header = T,sep = c(","),stringsAsFactors = F,strip.white = T,skipNul = T)
  #if (length(grep("evt",dsgrid.og$valuefrom))>0){
  #  dsgrid<-dsgrid.og[-grep("evt",dsgrid.og$valuefrom),]} else {dsgrid.og->dsgrid}
  #Generate signal with make signal with grid function (grid.csv need to be in working directory or specified otherwise)
  signals<-make_signal_with_grid(outputdata = output,add_taskness = T,dsgrid = dsgrid,nona = T)
  
  if (length(grep("evt",dsgrid$valuefrom))>0){
    dxgrid<-dsgrid[grep("evt",dsgrid$valuefrom),]
    for (u in 1:length(dxgrid$name)) {
      signals[dxgrid$name[u]]<-signals[dxgrid$valuefrom[u]]
    }
  }  
  #Get nuissance regressor: 
  #Still concat nuisa regressor together
  nuisa<-get_nuisance_preproc(id=paste0(tid,proc_id_subs),
                              cfgfilepath = cfgpath,
                              returnas = "data.frame",
                              dothese=nuisa_motion,
                              type=motion_type,
                              threshold=motion_threshold) 
  if(any(is.na(nuisa))) {stop("Can't get nuisance files from proc dirs")}
  if (convlv_nuisa) {
    nuisa.x<-nuisa
  } else {
    nuisa.x<-NULL
  }
  
  #Get the actual volume by run:
  run_volum<-get_volume_run(id=paste0(tid,proc_id_subs),
                            cfgfilepath = cfgpath,
                            returnas = "numbers",
                            reg.nii.name = func.nii.name)
  if(any(is.na(run_volum))) {stop("Can't get volume number from proc dirs")}
  
  
  if(is.null(model.varinames)){
    model.varinames<-dsgrid$name
    argu$model.varinames<-model.varinames
    }
  #Create  models:
  model<-signals[model.varinames]
  
  #Use Michael's package to generate design matrix and correlation graph;
  design<-dependlab::build_design_matrix(
    center_values=centerscaleall,
    events = output$event.list$allconcat, #Load the task info
    signals = model,     #Load the Model
    write_timing_files = wrt.timing, #Output timing files to FSL style
    tr=as.numeric(cfg$preproc_call$tr), #Grab the tr from cfg instead of hard coding it...
    plot = F,
    run_volumes = run_volum,
    output_directory = file.path(regpath,model.name,tid), #Where to output the timing files, default is the working directory
    nuisance_regressors = nuisa.x #Maybe could add in nuisance_regressors from pre-proc
  )
  
  design$heatmap<-make_heatmap_with_design(design)
  design$volume<-run_volum
  design$nuisan<-nuisa
  design$ID<-tid
  design$preprocID<-paste0(tid,proc_id_subs)
  design$regpath<-file.path(regpath,model.name,tid)
  # design$grid<-dsgrid
  
  if (!is.null(nuisa)){
    for (k in 1:length(nuisa)) {
      write.table(as.matrix(nuisa[[k]]),file.path(regpath,model.name,tid,
                                                  paste0("run",k,"_nuisance_regressor_with_motion.txt")),
                  row.names = F,col.names = FALSE)
    }}
  
  return(design)
  
}
######Modify fsl template with variable switch
change_fsl_template<-function(fsltemplate=NULL,begin="ARG_",end="_END",searchenvir=xarg,focus=T) {
  # for (tofind in grep(paste0(begin,"*"),fsltemplate) ){
  #   tryCatch(
  #     {
  #     varixma<-substr(fsltemplate[tofind],regexpr(paste0(begin,"*"),fsltemplate[tofind])+nchar(begin),
  #                      regexpr(paste0("*",end),fsltemplate[tofind])-1)
  #     if (!focus | varixma %in% objects(searchenvir)) {
  #       fsltemplate[tofind]<-gsub(paste0(begin,varixma,end),searchenvir[[varixma]],fsltemplate[tofind])
  #     } 
  #     },
  #     error=function(e){print(paste0("something went wrong in ",varixma," :",e))})
  # }
  # return(fsltemplate)
  x1<-paste0(begin,".*",end)
  x2<-paste0(".*",begin,"*(.*?) *",end,".*")
  if(focus){}
  fsltemplate[which(grepl(x1,fsltemplate))]<-unlist(lapply(fsltemplate[which(grepl(x1,fsltemplate))],function(stx){
    otx<-sub(x=stx,pattern = x2,replacement="\\1")
    ##print(otx)
    #data.frame(org=otx,replacement=searchenvir[[otx]])
    if(!exists(x = otx,envir = searchenvir) || is.null(searchenvir[[otx]]) || is.na(searchenvir[[otx]])) {stop(paste0("Variable ",otx," is NULL or NA."))}
    gsub(x1,replacement = searchenvir[[otx]],x = stx)
  }))
  
  return(fsltemplate)
}






#####Generate reg path from model name:

gen_reg<-function(vmodel=NULL,regpath=NULL,idx=NULL,runnum=NULL,env=NULL,regtype=NULL) {
  #message("Depreciated! Is probably useless; Will keep for 3 versions before deleting")
  NUP<-lapply(vmodel, function(x) {
    assign(paste0(x,"reg"),file.path(regpath,idx,paste0("run",runnum,"_",x,regtype)),envir = env)
  })
} 

populate_within<-function(chunk_within=NULL,xvlist=NULL,variname=NULL,firstorder=NULL){
  message("Depreciated! Is probably useless; Will keep for 3 versions before deleting")
  seq_lines<-lapply(xvlist,function(xjk){
    temp<-as.environment(list())
    assign(variname,xjk,envir = temp)
    result<-change_fsl_template(fsltemplate = chunk_within,begin = "XG_",end = "_EX",searchenvir = temp,focus=firstorder)
    if (firstorder) {
      result<-change_fsl_template(fsltemplate = result,begin = "BxG_",end = "_ExN",searchenvir = temp,focus=firstorder)
    }
    rm(temp)
    return(result)
  })
  re_within<-do.call(c,seq_lines)  
  return(re_within)
}

rep_within<-function(adptemplate=NULL,searchenvir=NULL){
  if(any(grepl("^SEC_.*_START$",adptemplate))){
    startnum<-grep("^SEC_.*_START$",adptemplate)[1]
    vari_x<-gsub("_START","",gsub("SEC_","",adptemplate[startnum]))
    endnum<-grep(paste0("SEC_",vari_x,"_END"),adptemplate)[1]
    chunk_pre<-adptemplate[1:startnum-1]
    chunk_post<-adptemplate[(endnum+1):length(adptemplate)]
    chunk_within<-adptemplate[(startnum+1):(endnum-1)]
    lsx<-tryCatch({
      lsx<-get(vari_x,envir = searchenvir)
    },error = function(e) {
      vari_y<-gsub("XG_","",gsub("_EX","",gsub("BxG_","",gsub("ExN","",vari_x))))
      return(get(vari_y,envir = searchenvir))
    })
    new_chunk<-do.call(c,lapply(lsx,function(repl){
      gsub(pattern = vari_x,replacement = repl,x = chunk_within,fixed = T)
    }))
    adptemplate<-c(chunk_pre,new_chunk,chunk_post)
    writeLines(adptemplate,"~/Desktop/test.fsf")
    rep_within(adptemplate=adptemplate,searchenvir = searchenvir)
  } else {return(adptemplate)}
}


adopt_feat<-function(adptemplate_path=NULL,searenvir=NULL,firstorder=F) {
  readLines(adptemplate_path)->adptemplate
  #find the area that needs to be replaced:
  #For real tho, the order should be fine but just to be 1000000 % sure, let's match them
  while(length(grep("SEC_.*_START",adptemplate))>0) {
    startnum<-grep("SEC_.*_START",adptemplate)[1]
    name<-gsub("_START","",gsub("SEC_","",adptemplate[startnum]))
    endnum<-grep(paste0("SEC_",name,"_END"),adptemplate)
    #Let's not use lapply cuz combining them would be a headache..
    #print(as.character(name))
    chunk_pre<-adptemplate[1:startnum-1]
    chunk_post<-adptemplate[(endnum+1):length(adptemplate)]
    #Okay now we constructure what's within.
    chunk_within<-adptemplate[(startnum+1):(endnum-1)]
    
    variname<-substr(chunk_within,regexpr(paste0("XG_","*"),chunk_within)+nchar("XG_"),
                     regexpr(paste0("*","_EX"),chunk_within)-1)
    unique(variname[variname!=""])->finavariname
    if(firstorder) {
      if(any(grepl("evnum",finavariname))) {finavariname<-"evnum"} else if (any(grepl("copenum",finavariname))) {
        finavariname<-"copenum"}  else {
          finavariname<-finavariname[1]}
    }
    get(finavariname,envir = searenvir)->xvlist
    if (length(finavariname)==1) {
      new_within<-populate_within(chunk_within = chunk_within,xvlist = xvlist,variname = finavariname,firstorder=T)
    } else {stop("Something Went Wrong! This is adopt_feat, first step")}
    
    if (length(new_within)>0) {
      adptemplate<-c(chunk_pre,new_within,chunk_post)
    } else {stop(paste0("Something Went Wrong! This is adopt_feat, 2nd step"))}
  } #End while 
  return(adptemplate)
}

feat_w_template<-function(templatepath=NULL,fsltemplate=NULL,beg="ARG_",end="_END",fsfpath=NULL,envir=NULL,outcommand=F) {
  if (is.null(fsltemplate)) {fsltemplate<-readLines(templatepath)}
  subbyrunfeat<-change_fsl_template(fsltemplate = fsltemplate,begin = beg,end=end,searchenvir = envir)
  #fsfpath<-fsf.path
  writeLines(subbyrunfeat,fsfpath)
  if(!outcommand){
    message("starting to do feat...")
    system(paste0("feat ",fsfpath),intern = T)
    message("feat completed")
  } else {return(paste0("feat ",fsfpath))}
}

plot_image_all<-function(rootpath=NULL,templatedir=NULL,model.name=NULL,patt=NULL,threshold=0.99,colour="red") {
  dirs<-system(paste0("find ",file.path(rootpath,model.name)," -iname '",patt,"' -maxdepth 4 -mindepth 1 -type f"),intern = T)
  for (sdir in dirs) {
    spdir<-strsplit(sdir,.Platform$file.sep) 
    spdir[[1]][sapply(spdir, function(x) {grep(patt,x)})-1]->filex
    paste(spdir[[1]][-c(length(spdir[[1]]),length(spdir[[1]])-1)],collapse = .Platform$file.sep)->outputdir
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

#####Group Lvl Func
glvl_all_cope<-function(rootdir="/Volumes/bek/neurofeedback/sonrisa1/nfb/ssanalysis/fsl",
                        outputdir="/Volumes/bek/neurofeedback/sonrisa1/nfb/grpanal/fsl",
                        modelname="PE_8C_old",grp_sep=argu$group_id_sep,copestorun=1:8,
                        paralleln=NULL,onesamplet_pergroup=T,pairedtest=T,unpairedtest=T,
                        thresh_cluster_mass=NULL,thresh_cluster_extent=NULL,pvalue=0.001, ifDeMean=T,
                        usethesetest=c("tfce","voxel-based","cluster-based-mass","cluster-based-extent")) {
  if ( is.null(modelname) ) {stop("Must specify a model name other wise it will be hard to find all copes")}
  if(!exists("supplyidmap",envir = argu)){argu$supplyidmap<-sortid(dix=file.path(rootdir,modelname),idgrep=grp_sep,dosymlink=F)}
  #Creating directory in case they are not there;
  dir.create(file.path(outputdir,modelname),showWarnings = FALSE,recursive = T)
  #Ensure fsl are in path:
  fsl_2_sys_env()
  
  raw<-system(paste0("find ",file.path(rootdir,modelname,"*/average.gfeat")," -iname '*.feat' -maxdepth 2 -mindepth 1 -type d"),intern = T)
  strsplit(raw,split = "/") ->raw.split
  df.ex<-data.frame(ID=unlist(lapply(raw.split,function(x) {
    x[grep("average.gfeat",x)-1]
  })),
  COPENUM=unlist(lapply(raw.split,function(x) {
    x[grep("average.gfeat",x)+1]
  })),
  PATH=file.path(raw,"stats","cope1.nii.gz")
  )
  df.ex$COPENUM<-substr(df.ex$COPENUM,start=regexpr("[0-9]",df.ex$COPENUM),stop = regexpr(".feat",df.ex$COPENUM)-1)
  if (max(aggregate(COPENUM~ID,data = df.ex,max)$COPENUM)<max(copestorun)) {stop("HEY! There's not that many copes to run! Change argument!")}
  noIDpos<-which(aggregate(COPENUM~ID,data = df.ex,max)$COPENUM!=max(aggregate(COPENUM~ID,data = df.ex,max)$COPENUM) & aggregate(COPENUM~ID,data = df.ex,max)$COPENUM<max(copestorun))
  if (length(noIDpos)>0){
    noID<-aggregate(COPENUM~ID,data = df.ex,max)$ID[noIDpos]
    message(paste("This ID:",noID,"does not have enough COPES, will be removed from running...."))
    df.ex[which(!df.ex$ID %in% noID),]->df.ex
  } else {message("All Good!")}
  message("Now will run fslmerge and randomise, function will fail if FSL ENVIR is not set up. (Should not happen since it's guarded by func)")
  
  #we now add in the group level stuff here, in hope that it will have more flexiblity;
  if (length(grp_sep)>1) {
    df.jx<-merge(df.ex,do.call(rbind,lapply(argu$supplyidmap,function(x) {data.frame(ID=x$ID,GROUP=x$name)})),by = "ID",all.x = T)
  } else {df.ex$GROUP<-"ONE"
  df.jx<-df.ex}
  
  allcopecomx<-as.environment(list())
  
  #DO the options:
  option_grp<-""
  if(ifDeMean) {option_grp<-paste0(option_grp," -D")}
  if("tfce" %in% usethesetest){option_grp<-paste0(option_grp," -T")}
  if("voxel-based" %in% usethesetest){option_grp<-paste0(option_grp," -x")}
  if("cluster-based-extent" %in% usethesetest){
    if(is.null(thresh_cluster_extent)) {thresh_cluster_extent<-round(abs(qt(p = pvalue,df = (length(unique(df.ex$ID))-1) )),4)}
    option_grp<-paste0(option_grp," -c ",thresh_cluster_extent)
  }
  if("cluster-based-mass" %in% usethesetest){
    if(is.null(thresh_cluster_mass)) {thresh_cluster_mass<-round(abs(qt(p = pvalue,df = (length(unique(df.ex$ID))-1) )),4)}
    option_grp<-paste0(option_grp," -C ",thresh_cluster_mass)
  }
  
  #Now we face this problem of runing bunch of them...
  #One sample T test;
  if (length(unique(df.jx$GROUP))==1){ 
    #Make Commands;
    cope.fslmerge<-lapply(copestorun,function(x) {
      outputroot<-file.path(outputdir,modelname,paste0("cope",x,"_randomize_onesample_ttest"))
      dir.create(outputroot, showWarnings = FALSE,recursive = T)
      if (length(list.files(pattern = "*_tstat1",path = outputroot,no.. = T))<1) {
        copefileconcat<-paste(df.jx$PATH[which(df.jx$COPENUM==x)],collapse = " ")
        return(paste0("fslmerge -t ",outputroot,"/OneSamp4D"," ",
                      copefileconcat
                      ," \n ",
                      "randomise -i ",outputroot,"/OneSamp4D -o ",outputroot,"/OneSampT -1",option_grp
        ))
      }else {return(NULL)}
    })
    cleanuplist(cope.fslmerge)->cope.fslmerge
    assign(x = "onesamplet_onegroup",value = cope.fslmerge,envir = allcopecomx)
  } else {
    if (onesamplet_pergroup) {
      #Make symoblic link first
      NX<-sortid(dix=file.path(rootdir,modelname),idgrep=grp_sep,dosymlink=F)
      copexgroup<-do.call(c,lapply(copestorun,function(zt) {unlist(paste(zt,grp_sep,sep = "_"))}))
      cope.fslmerge<-lapply(copexgroup,function(x) {
        unlist(strsplit(x,split = "_"))->cope_group
        outputroot<-file.path(outputdir,modelname,cope_group[2],paste0("cope",x,"_randomize_onesample_ttest"))
        dir.create(outputroot, showWarnings = FALSE,recursive = T)
        if (length(list.files(pattern = "*_tstat1",path = outputroot,no.. = T))<1) {
          copefileconcat<-paste(as.character(df.jx$PATH[which(df.jx$COPENUM==cope_group[1] & df.jx$GROUP==cope_group[2])]),collapse = " ")
          return(paste0("fslmerge -t ",outputroot,"/OneSamp4D"," ",
                        copefileconcat
                        ," \n ",
                        "randomise -i ",outputroot,"/OneSamp4D -o ",outputroot,"/OneSampT -1",option_grp
          ))
        }else {return(NULL)}
      })
      cleanuplist(cope.fslmerge)->cope.fslmerge
      assign(x = "onesamplet_pergroup",value = cope.fslmerge,envir = allcopecomx)
    } 
    if (pairedtest) {
      dir.create(file.path(outputdir,modelname),showWarnings = F,recursive = T)
      if(!exists("supplyidmap",envir = argu)){argu$supplyidmap<-sortid(dix=file.path(rootdir,modelname),idgrep=grp_sep,dosymlink=F)}
      
      commonid<-prep_paired_t(idsep = argu$supplyidmap, outpath = file.path(outputdir,modelname))
      cidindex<-data.frame(ID=do.call(c,lapply(commonid,function(xm) {paste(xm,grp_sep,sep = "_")})),notag=do.call(c,lapply(commonid,function(xm) {paste(xm,c("",""),sep = "")})))
      df.kh<-merge(df.jx,cidindex,all = T)
      df.jk<-df.kh[which(!is.na(df.kh$notag)),]
      cope.fslmerge<-lapply(copestorun,function(x) {
        outputroot<-file.path(outputdir,modelname,paste0("cope",x,"_randomize_paired_ttest"))
        dir.create(outputroot, showWarnings = FALSE,recursive = T)
        if (length(list.files(pattern = "*tfce_corrp_tstat1",path = outputroot,no.. = T))<1) {
          onecope<-df.jk[which(df.jk$COPENUM==x),]
          onecope_re<-onecope[with(onecope,order(GROUP,notag)),]
          copefileconcat<-paste(onecope_re$PATH,collapse = " ")
          return(paste0("fslmerge -t ",outputroot,"/PairedT4D"," ",
                        copefileconcat
                        ," \n ",
                        "randomise -i ",outputroot,"/PairedT4D -o ",outputroot,"/PairedT -d ",
                        file.path(outputdir,modelname),"/design.mat -t ",
                        file.path(outputdir,modelname),"/design.con -e ",
                        file.path(outputdir,modelname),"/design.grp",option_grp
          ))
        }else {return(NULL)}
      })
      cleanuplist(cope.fslmerge)->cope.fslmerge
      assign(x = "pairedtests",value = cope.fslmerge,envir = allcopecomx)
    }
    if (unpairedtest) {
      dir.create(file.path(outputdir,modelname),showWarnings = F,recursive = T)
      if(!exists("supplyidmap",envir = argu)){argu$supplyidmap<-sortid(dix=file.path(rootdir,modelname),idgrep=grp_sep,dosymlink=F)}
      IDorder<-prep_unpaired_t(idsep = argu$supplyidmap, outpath = file.path(outputdir,modelname))

      df.jr<-df.jx[which(df.jx$ID %in% IDorder),]
      if(length(IDorder[which(is.na(match(IDorder,df.jx$ID)))])>0){stop("Some of the ID supplied does not have data, they are: ",IDorder[which(is.na(match(IDorder,df.jx$ID)))])}
      cope.fslmerge<-lapply(copestorun,function(x) {
        outputroot<-file.path(outputdir,modelname,paste0("cope",x,"_randomize_unpaired_ttests"))
        dir.create(outputroot, showWarnings = FALSE,recursive = T)
        if (length(list.files(pattern = "*tfce_corrp_tstat1",path = outputroot,no.. = T))<1) {
          onecope<-df.jr[which(df.jr$COPENUM==x),]
          onecope_re<-onecope[match(IDorder,onecope$ID),]
          copefileconcat<-paste(onecope_re$PATH,collapse = " ")
          return(paste0("fslmerge -t ",outputroot,"/UnpairedT4D"," ",
                        copefileconcat
                        ," \n ",
                        "randomise -i ",outputroot,"/UnpairedT4D -o ",outputroot,"/UnpairedT -d ",
                        file.path(outputdir,modelname),"/design_unpaired.mat -t ",
                        file.path(outputdir,modelname),"/design_unpaired.con ",option_grp
          ))
        }else {return(NULL)}
      })
      cleanuplist(cope.fslmerge)->cope.fslmerge
      assign(x = "pairedtests",value = cope.fslmerge,envir = allcopecomx)
    }
  }
  XNN<-eapply(env = allcopecomx, FUN = function(cope.fslmerge) {
    sink(file=file.path(outputdir,modelname,"glvl_log.txt"),split=TRUE)
    #Do Parallel or nah
    if (!is.null(paralleln)){
      cj1<-makeCluster(paralleln,outfile="",type = "FORK")
      NU<-parSapply(cj1,cope.fslmerge,function(x) {
        message(paste0("Now running ",x))
        pb<-txtProgressBar(min = 0,max = 100,char = "|",width = 50,style = 3)
        numdx<-which(x==cope.fslmerge)
        indx<-suppressWarnings(split(1:length(cope.fslmerge),1:paralleln))
        pindx<-grep(paste0("\\b",numdx,"\\b"),indx)
        setTxtProgressBar(pb,(which(numdx==indx[[pindx]]) / length(indx[[pindx]]))*100)
        system(command = x,intern = T,ignore.stdout = F,ignore.stderr = F)
        message("completed")
      })
      stopCluster(cj1)
    } else {
      lapply(cope.fslmerge,function(x) {
        message(paste0("Now running ",x))
        pb<-txtProgressBar(min = 0,max = 100,char = "|",width = 50,style = 3)
        numdx<-which(x==cope.fslmerge)
        indx<-suppressWarnings(split(1:length(cope.fslmerge),1))
        pindx<-grep(paste0("\\b",numdx,"\\b"),indx)
        setTxtProgressBar(pb,(which(numdx==indx[[pindx]]) / length(indx[[pindx]]))*100)
        
        system(command = x,intern = T,ignore.stdout = F,ignore.stderr = F)
        
        message("completed")
      })
    }
  })
  message("DONE")
}


sortid<-function(dix=file.path(argu$ssub_outputroot,argu$model.name),idgrep=argu$group_id_sep,dosymlink=ifdosymlink){
  system(paste0("find -L ",dix," -iname 'average.gfeat' -maxdepth 2"),intern = T)->dirxs
  alldirs<-sapply(strsplit(dirxs,split = .Platform$file.sep),function(x) {x[(length(x)-1)]})
  split_dirx<-lapply(idgrep,function(x){
    j<-alldirs[grep(x,alldirs)]
    if(dosymlink) {
      dir.create(file.path(argu$ssub_outputroot,argu$model.name,x),showWarnings = F)
      file.symlink(from = file.path(argu$ssub_outputroot,argu$model.name,alldirs[grep(x,alldirs)]),
                   to = file.path(argu$ssub_outputroot,argu$model.name,x,alldirs[grep(x,alldirs)]))
    }
    return(list(ID=j,name=x))
  })
  names(split_dirx)<-idgrep
  return(split_dirx)    
}

prep_unpaired_t<-function(idsep=NULL,outpath=NULL){
  ngrp<-length(idsep)
  nsub<-sapply(idsep,function(xr){length(xr$ID)})
  allnames<-sapply(idsep,function(xr){xr$name})
  dmat<-do.call(rbind,lapply(1:ngrp,function(rz){
    dmatx<-matrix(data = 0,ncol = ngrp,nrow = nsub[rz])
    dmatx[,rz]<-1 
    return(dmatx)
  }))
  
  pairs<-levels(interaction(as.factor(sapply(idsep,function(s){s$name})),as.factor(sapply(idsep,function(s){s$name})),sep = ">"))
  cmat_ls<-lapply(cleanuplist(lapply(strsplit(pairs,split = ">"),function(x){
    ifelse(any(duplicated(x)),return(NULL),return(x))})),function(sr){
      cmat<-matrix(nrow = 1,ncol = ngrp,data = 0)
      cmat[,match(sr[1],allnames)]<-1
      cmat[,match(sr[2],allnames)]<- -1
      return(list(name=paste(sr,collapse = ">"),cmat=cmat))
      })
  cmat<-do.call(rbind,lapply(cmat_ls,function(r){r$cmat}))
  gmat<-do.call(rbind,lapply(1:ngrp,function(io){matrix(ncol = 1,nrow = nsub[io],data = io)}))
  
  write.table(dmat,file = file.path(outpath,"design.mat.txt"),row.names = F,col.names = F)
  system(paste0("${FSLDIR}/bin/Text2Vest ",file.path(outpath,"design.mat.txt")," ",file.path(outpath,"design_unpaired.mat")))
  file.remove(file.path(outpath,"design.mat.txt"))
  
  write.table(cmat,file = file.path(outpath,"design.con.txt"),row.names = F,col.names = F)
  system(paste0("${FSLDIR}/bin/Text2Vest ",file.path(outpath,"design.con.txt")," ",file.path(outpath,"design_unpaired.con")))
  file.remove(file.path(outpath,"design.con.txt"))
  
  write.table(gmat,file = file.path(outpath,"design.grp.txt"),row.names = F,col.names = F)
  system(paste0("${FSLDIR}/bin/Text2Vest ",file.path(outpath,"design.grp.txt")," ",file.path(outpath,"design_unpaired.grp")))
  file.remove(file.path(outpath,"design.grp.txt"))
  
  write.table(data.frame(contrastnum=1:nrow(cmat),constrastname=unlist(lapply(cmat_ls, function(x){x$name}))),
              file = file.path(outpath,"unpaired_t_contrastnames.csv"))
  return(as.character(unlist(sapply(idsep,function(x){x$ID}))))    
  #return(list(dmat,cmat,gmat))
}

prep_paired_t<-function(idsep=NULL,outpath=NULL){
  
  commonid<-Reduce(intersect, lapply(names(idsep),function(xj){
    idsep[[xj]]$ID->xk
    gsub(pattern = paste0("_",xj),replacement = "",x = xk)
  }))
  
  
  write.table(
    do.call(rbind, lapply(c(-1,1,2:100)[1:length(idsep)], function(x){
      cbind(matrix(nrow = length(commonid),ncol = 1,data = x),diag(x = 1,nrow = length(commonid),ncol = length(commonid)))
    })),file = file.path(outpath,"design.mat.txt"),row.names = F,col.names = F)
  system(paste0("${FSLDIR}/bin/Text2Vest ",file.path(outpath,"design.mat.txt")," ",file.path(outpath,"design.mat")))
  file.remove(file.path(outpath,"design.mat.txt"))
  
  write.table(
    rbind(
      diag(x = 1,nrow = 1, ncol = length(commonid)+1),
      diag(x = -1,nrow = 1, ncol = length(commonid)+1)
    )
    ,file = file.path(outpath,"design.con.txt"),row.names = F,col.names = F)
  system(paste0("${FSLDIR}/bin/Text2Vest ",file.path(outpath,"design.con.txt")," ",file.path(outpath,"design.con")))
  file.remove(file.path(outpath,"design.con.txt"))
  
  write.table(do.call(rbind, lapply(1:length(idsep), function(x){
    matrix(nrow = length(commonid),ncol = 1,data = seq_along(commonid))
  })),file = file.path(outpath,"design.grp.txt"),row.names = F,col.names = F)
  system(paste0("${FSLDIR}/bin/Text2Vest ",file.path(outpath,"design.grp.txt")," ",file.path(outpath,"design.grp")))
  file.remove(file.path(outpath,"design.grp.txt"))
  
  return(commonid)
  
}


get_motion_info<-function(configpath=NULL,type="fd",threshold="default"){
  if (!file.exists(configpath)){stop("Config File Does NOT Exist")}
  cfg<-cfg_info(configpath)
  idlist<-list.dirs(cfg$loc_mrproc_root,recursive = F,full.names = F)
  NX<-lapply(idlist,function(X) {
    tryCatch(
      {return(get_nuisance_preproc(X,cfgfilepath=configpath,
                                   returnas=c("data.frame"),
                                   dothese=c("motion_outlier"),
                                   type=type,
                                   threshold=threshold))
      },error=function(e) {return(NULL)})
  })
  names(NX)<-idlist
  NX<-cleanuplist(NX)
  NU<-lapply(1:length(NX),function(i){
    yx<-suppressMessages(reshape2::melt(as.data.frame(lapply(NX[[i]],length))))
    names(yx)<-c("run","outlier")
    yx$run<-as.numeric(gsub("[a-z]*[A-Z]*","",yx$run))
    IDx<-names(NX)[i]
    yx$totalvol<-get_volume_run(IDx,cfgfilepath = configpath,returnas = "numbers")
    yx$out_per<-yx$outlier / yx$totalvol
    yx$ID<-IDx
    return(yx)
  })
  return(do.call(rbind,NU))
}

amputate_run<-function(small.sub=NULL,cfgpath=NULL,type="fd",threshold="default") {
  
  stop("UNFINISHED FUNC")
  
}

####QC


gen_model_arg<-function(cfgpath=NULL,func.nii.name="nfswudktm*[0-9]_[0-9].nii.gz",mni_template=NULL,QC_auxdir="/Volumes/bek/QC_fsl",parallen=4,fullmodel=F,...){
  cfg<-cfg_info(cfgpath = cfgpath)
  npaths<-lapply(c("ssanalysis/fsl","regs","grpanal/fsl"),function(xx) {file.path(cfg$loc_root,cfg$paradigm_name,xx)})
  NX<-lapply(npaths,dir.create,showWarnings = F,recursive = T)
  if(fullmodel) {addarg<-list(adaptive_gfeat=TRUE,gsub_fsl_templatepath=file.path(QC_auxdir,"fsl_gfeat_general_adaptive_template.fsf"),
                              glvl_output=npaths[[3]],hig_lvl_path_filter=NULL,graphic.threshold=0.95,convlv_nuisa=F,
                              nuisa_motion=c("nuisance","motion_par"),motion_type="fd",motion_threshold="default")} else {addarg=list()} 
  argu<-as.environment(c(addarg,list(nprocess=parallen,onlyrun=1:3,proc_id_subs=NULL,
                                     regtype=".1D",ifnuisa=FALSE,ifoverwrite_secondlvl=F,cfgpath=cfgpath,
                                     forcereg=F,regpath=npaths[[2]],func.nii.name=func.nii.name,ssub_outputroot=npaths[[1]],
                                     templatedir=mni_template,
                                     nuisa_motion=c("nuisance","motion_par"),convlv_nuisa=F,
                                     QC_auxdir=QC_auxdir,motion_type="fd",motion_threshold="default",
                                     gridpath=file.path(QC_auxdir,"qc_grid.csv"),
                                     model.name="QC",
                                     model.varinames=c("QC","QC_none"),
                                     ssub_fsl_templatepath=file.path(QC_auxdir,"qc_fsl_template.fsf"),
                                     ...
                                     
  )))
  return(argu)
}

check_incomplete_preproc<-function(cfgpath=NULL,enforce=F,verbose=T) {
  cfg<-cfg_info(cfgpath)
  idstocheck<-list.files(cfg$loc_mrproc_root)
  runnums<-as.numeric(cfg$n_expected_funcruns)
  outerrors<-lapply(idstocheck,function(cid) {
    #print(cid)
    proc_num<-get_volume_run(id = cid,cfgfilepath = cfgpath,returnas = "numbers")
    if (any(is.na(proc_num))) {proc_num[which(is.na(proc_num))]<-0}
    func_dir_raw<-system(paste0("find ",file.path(cfg$loc_mrraw_root,cid)," -iname ",cfg$functional_dirpattern," -maxdepth 2 -mindepth 1"),intern = T)
    TJ<-lapply(func_dir_raw,list.files,recursive=F)
    raw_num<-sapply(TJ, length)
    if(length(raw_num)>0) {
      if(any(proc_num!=raw_num) & verbose) {
        message("#################################################")
        message(cid)
        message(paste("run",which(proc_num!=raw_num)))
        message(paste("Raw:",paste(raw_num,collapse = " ")))
        message(paste("Proc:",paste(proc_num,collapse = " ")))
        message("#################################################")
      }
      outerror<-data.frame(proc_num,raw_num)
      outerror$Run<-1:length(raw_num)
      outerror$ID<-cid
      return(outerror)
    } else {return(NULL)}
  })
  allerrors<-do.call(rbind,outerrors)
  suballerrors<-allerrors[which(allerrors$proc_num != allerrors$raw_num),]
  return(suballerrors)
  if (enforce) {
    extd<-suballerrors[suballerrors$proc_num!=0,]
    print(extd)
    ifrun<-as.logical(readline(prompt = "Review above info, please type T/TRUE to continue or F/FALSE to stop: "))
    if (ifrun) {
      for (o in 1:length(extd$ID)) {
        cxtd<-extd[o,]
        unlink(file.path(cfg$loc_mrproc_root,cxtd$ID,paste0(cfg$paradigm_name,"_proc"),paste0(cfg$paradigm_name,cxtd$Run)),recursive = T,force = T)
      }
    }
  }
}
###############

# create_roimask_atlas<-function(atlas_name=NULL,atlas_xml=NULL,target=NULL,outdir=tempdir(),atlas_root=NULL,
#                                fsl_dir=Sys.getenv("FSLDIR"),volxsize="2mm",type="",singlemask=T)



create_roimask_atlas<-function(atlas_name=NULL,atlas_xml=NULL,target=NULL,outdir=tempdir(),atlas_root=NULL,
                               fsl_dir=Sys.getenv("FSLDIR"),volxsize="2mm",type="",singlemask=T,atlas_readtype=c("fsl","spm"),output_main=F) {
  if (is.null(atlas_root)){
    if (is.null(fsl_dir)) {
      fsl_2_sys_env(); fsl_dir=Sys.getenv("FSLDIR")}
    atlas_dir<-file.path(fsl_dir,"data","atlases")    
    maybefsl<-TRUE
  } else {atlas_dir<-atlas_root}
  
  if(is.null(atlas_readtype)){
    if(maybefsl){source_type<-"fsl"}else{source_type<-"spm"}
  }else{
    if(length(atlas_readtype)==1 & c("fsl" %in% atlas_readtype)){source_type<-"fsl"
    } else if (length(atlas_readtype)==1 & c("spm" %in% atlas_readtype)){source_type<-"spm"
    } else {stop("unknown source type")}}
  
  if (is.null(atlas_xml)) {
    if (is.null(atlas_name)) {
      message("Below are the available atlases in atlas folder:")
      print(latlas<-gsub(".xml","",list.files(path = atlas_dir,pattern = "*.xml"))) 
      wic<-readline(prompt = "Which atlas to use? Type in exact: ")
      if (!wic %in% latlas) {stop(paste0("No atlas named ",wic," found in default fsl atlas directory. Try again"))
      } else {atlas_name<-wic} 
    }
    atlas_xml<-file.path(atlas_dir,paste0(atlas_name,".xml"))
  } else if (!grepl("/",atlas_xml)) {atlas_xml<-file.path(atlas_dir,atlas_xml)}
  atlas_info<-XML::xmlToList(XML::xmlParse(file = atlas_xml))
  images<-atlas_info$header[grep("images",names(atlas_info$header))]
  imagex<-cleanuplist(lapply(images,function(img) {if(length(grep(volxsize,img$imagefile))>0) {return(img)} else {return(NULL)}}))
  #fsl
  if(source_type=="fsl"){
    atrx<-do.call(rbind,lapply(atlas_info$data,function(dt) {return(cbind(data.frame(target=dt$text),data.frame(as.list(dt$.attrs),stringsAsFactors = F)))}))
    atrx$index<-as.numeric(atrx$index)+1
    target_imag<-file.path(atlas_dir,imagex$images$summaryimagefile)
    atrx$thres<-1:nrow(atrx)
  }
  #spm
  if(source_type=="spm"){
    atrx<-do.call(rbind,lapply(atlas_info$data,function(dt) {
      return(data.frame(index=dt$index,target=dt$name,stringsAsFactors = F))
    }))
    target_imag<-file.path(atlas_dir,imagex$images$imagefile)
    atrx$thres<-atrx$index
    volxsize<-""
  }
  
  atrx$maskdir<-NA
  atrx$total_voxel<-NA
  
  if (is.null(target)){target<-as.numeric(atrx$index)}
  if (any(is.character(target))) {
    tarindxs<-as.character(atrx$index[grep(pattern = paste(target,collapse = "|"),x = atrx$target)])
  } else {tarindxs<-as.character(target)}
  
  dir.create(path = outdir,recursive = T,showWarnings = F)
  for (sindsk in tarindxs) {
    tarwhich<-which(atrx$index==sindsk)
    tarthres <- atrx$thres[tarwhich]
    tartext<-gsub(" ","_",atrx$target[tarwhich])
    outfile<-file.path(outdir,paste(atlas_name,volxsize,tartext,"bin.nii.gz",sep="_"))
    if(!file.exists(outfile)){
      opt<-paste0("-thr ",tarthres," -uthr ",tarthres," -bin \"",outfile,"\"")
      cmd<-paste("fslmaths",target_imag,opt)
      system(cmd,intern = F)
    }
    atrx$maskdir[tarwhich]<-outfile
    atrx$total_voxel[tarwhich]<-voxel_count(cfile = outfile)[1]
  }
  
  if (singlemask && length(tarindxs)>1) {
    singlemask<-file.path(outdir,paste(atlas_name,volxsize, gsub(" ","_",paste(atrx$target[tarindxs],collapse = "_")),"bin.nii",sep="_"))
    cmd<-paste("${FSLDIR}/bin/fslmaths",paste(na.omit(atrx$maskdir),collapse = " -add "),singlemask)
    system(cmd,intern = F)
  } else {singlemask<-NULL}
  if(output_main){return(list(indexdf=atrx,singlemask=singlemask,main_img=target_imag))
  }else{return(list(indexdf=atrx,singlemask=singlemask))}
}

voxel_count<-function(cfile=NULL,mask=NULL,nonzero=T) {
  if (!is.null(mask)) {maskargu=paste0("-k ",mask)}else{maskargu<-""}
  if (nonzero) {ty<-"-V"} else {ty<-"-v"}
  cmd<-paste("${FSLDIR}/bin/fslstats",cfile,maskargu,ty,sep = " ")
  resultx<-system(cmd,intern = T)
  lres<-as.numeric(strsplit(resultx,split = " ")[[1]])
  names(lres)<-c("voxels","volumes")
  return(lres)
}

get_voxel_count<-function(cfile=NULL,stdsfile=NULL,intmat=NULL,mask=NULL,outdir=tempdir()) {
  if (outdir==tempdir() | !file.exists(file.path(outdir,"temp_1.nii"))) {
    cmd<-paste("${FSLDIR}/bin/flirt",
               "-ref",stdsfile,
               "-in",cfile,
               "-out",file.path(outdir,"temp_1.nii"),
               "-applyxfm","-init",intmat,"-interp","trilinear","-datatype","float")
    system(cmd,intern = F)}
  resultx<-voxel_count(cfile = file.path(outdir,"temp_1.nii"),mask = mask)
  return(resultx)
}

########################roi extraction;
whichfile<-function(textx=NULL,dirx=NULL,fullnameq=T){
  if("tstat" %in% textx){fxj<-list.files(dirx,pattern = "*T_tstat.*nii.gz",full.names = fullnameq)}
  if("voxel-based" %in% textx){fxj<-list.files(dirx,pattern = "*T_vox_corrp_tstat.*nii.gz",full.names = fullnameq)}
  if("tfce" %in% textx){fxj<-list.files(dirx,pattern = "*T_tfce_corrp_tstat.*nii.gz",full.names = fullnameq)}
  if("cluster-based-extent" %in% textx){fxj<-list.files(dirx,pattern = "*T_clustere_corrp_tstat.*nii.gz",full.names = fullnameq)}
  if("cluster-based-mass" %in% textx){fxj<-list.files(dirx,pattern = "*T_clusterm_corrp_tstat.*nii.gz",full.names = fullnameq)}
  return(fxj)
}  

# rootdir="/Volumes/bek/neurofeedback/sonrisa1/nfb/ssanalysis/fsl"
# grproot="/Volumes/bek/neurofeedback/sonrisa1/nfb/grpanal/fsl"
# modelname="Value1n"
# basemask="tstat"
# corrp_mask="tfce"
# saveclustermap=TRUE
# VersionCode="tfce0.95"
# corrmaskthreshold=0.98
# roimaskthreshold=0.0001
# voxelnumthres=10
# clustertoget=NULL
# copetoget=NULL
# maxcore=4

# featdir=featdir
# base=basemask
# corrp_mask=corrp_mask
# outdir = cmoutdir
# VersionCode = Version
# maskthresholdvalue=corrmaskthreshold
# roimaskthreshold=roimaskthreshold
# overwrite=!saveclustermap

gen_cluster_mask<-function(featdir="/Volumes/bek/neurofeedback/sonrisa1/nfb/grpanal/fsl/alignment6/cope12_randomize_onesample_ttest",useMMcor=T,
                           outdir=NULL,VersionCode=NULL,base="tstat",corrp_mask="tstat",maskthresholdvalue=3.0,roimaskthreshold=0.0001,
                           overwrite=T,writecsv=T,savetempdir=NULL) {
  if(is.null(VersionCode)){versionCode<-""}else{versionCode<-paste0("_",VersionCode)}
  if(is.null(outdir)){outdir<-file.path(featdir,paste0("ROI_masks",versionCode))}
  if(is.null(savetempdir)){tempdirx<-file.path(featdir,paste0("ROI_masks",versionCode))}else{tempdirx<-savetempdir}
  #if(is.null(templatebrain)) {templatebrain<-file.path(Sys.getenv("FSLDIR"),"data","standard","MNI152_T1_2mm_brain_mask.nii.gz")}
  dir.create(tempdirx,showWarnings = F,recursive = F)
  step1path<-file.path(tempdirx,"signifi_thres.nii.gz")
  cmd1<-paste("fslmaths",whichfile (corrp_mask,featdir,T),"-thr",maskthresholdvalue,
              "-bin","-mul",whichfile(base,featdir,T),step1path,sep = " ")
  system(cmd1)
  
  if(useMMcor){corOpt <- " --mm"}else{corOpt<-""}
  if(voxel_count(step1path)[1]<2) {
    nullmask<-TRUE
    message("This set up result in less than 2 voxels.")
  } else {nullmask<-FALSE}
  #Step 1; get only the 
  if (!nullmask){
    step2path<-file.path(tempdirx,"cluster_index.nii.gz")
    cmd2<-paste(sep = " ","cluster","-i",step1path,"-t ",roimaskthreshold,"-o",step2path,corOpt)
    s2raw<-system(cmd2,intern = T)
    s2df<-read.table(text = s2raw,sep = "\t",header = T,check.names = F)
    s2df$name<-s2df$`Cluster Index` #Here's where you can do something about getting their indexs
    if(overwrite){unlink(outdir,recursive = T)}
    dir.create(path = outdir,recursive = T,showWarnings = F)
    outpaths<-lapply(s2df$`Cluster Index`,function(cix){
      outfile<-file.path(outdir,paste0("cluster_mask_",cix,".nii.gz"))
      if(!file.exists(outfile)){
        opt<-paste0("-thr ",cix," -uthr ",cix," -bin \"",outfile,"\"")
        cmd<-paste("fslmaths",step2path,opt)
        system(cmd,intern = F)}
      return(outfile)
    })
    names(outpaths)<-s2df$name
    s2df$maskpath<-unlist(outpaths)
    if(writecsv) {write.csv(s2df,file.path(outdir,"index.csv"),row.names = F)}
    return(s2df)
  } else {return(NULL)}
  if(is.null(savetempdir)){unlink(tempdirx)}
}

#Right
#paste0("fslmaths ",templatebrain," -roi 0 45 0 -1 -1 -1 -1 -1 righthem_standard.nii")
#Left
#paste0("fslmaths ",templatebrain," -roi 45 -1 0 -1 0 -1 0 -1 lefthem_standard.nii")

# rootdir="/Volumes/bek/neurofeedback/sonrisa1/nfb/ssanalysis/fsl"
# grproot="/Volumes/bek/neurofeedback/sonrisa1/nfb/grpanal/fsl"
# modelname="Value1n"
# basemask="tstat"
# corrp_mask="tfce"
# saveclustermap=TRUE
# VersionCode="tfce0.95"
# corrmaskthreshold=0.98
# roimaskthreshold=0.0001
# voxelnumthres=10
# clustertoget=NULL
# copetoget=NULL
# maxcore=4

roi_getvalue<-function(rootdir=argu$ssub_outputroot,grproot=argu$glvl_outputroot,modelname=argu$model.name,grp_identif=NA,forceconcat=F,
                       basemask="tstat",corrp_mask="tstat",saveclustermap=TRUE,Version="t_t",corrmaskthreshold=3.0,useMMcor=F,
                       roimaskthreshold=0.0001, voxelnumthres=5, clustertoget=NULL,copetoget=NULL,maxcore=4,saverdata=T,...){
  #clustertoget=list(`12`=c(43,44),`13`=c(26,25)),copetoget=12:13){ #This is completely optional
  raw_avfeat<-system(paste0("find ",file.path(rootdir,modelname,"*/average.gfeat")," -iname '*.feat' -maxdepth 2 -mindepth 1 -type d"),intern = T)
  fsl_2_sys_env()
  if(is.null(Version)){Version<-paste0(corrp_mask,corrmaskthreshold)}
  strsplit(raw_avfeat,split = "/") ->raw.split
  df.ex<-data.frame(ID=unlist(lapply(raw.split,function(x) {
    x[grep("average.gfeat",x)-1]
  })),
  COPENUM=unlist(lapply(raw.split,function(x) {
    x[grep("average.gfeat",x)+1]
  })),
  PATH=file.path(raw_avfeat,"stats","cope1.nii.gz")
  )
  df.ex$COPENUM<-substr(df.ex$COPENUM,start=regexpr("[0-9]",df.ex$COPENUM),stop = regexpr(".feat",df.ex$COPENUM)-1)
  if(!is.na(grp_identif)){
    if(length(grp_identif)>1){stop("Function do not want to handle more than one grp_identifier at a time, do a lapply or loop.")}
    truerootdir<-file.path(grproot,modelname,grp_identif)
  } else {truerootdir<-file.path(grproot,modelname)}
  
  rnddirs<-list.dirs(path = truerootdir,recursive = F)
  if(saveclustermap){cmoutdir<-NULL}else{cmoutdir<-base::tempdir()}
  if(is.null(copetoget)){copetoget<-unique(as.character(df.ex$COPENUM))}
  if(length(copetoget)<maxcore & length(copetoget)>1){coresx<-length(copetoget)}else{coresx<-4}
  sharedproc<-parallel::makeCluster(coresx,outfile="",type = "FORK")
  copes_roivalues<-parallel::parLapply(cl=sharedproc,X = copetoget,function(copenum){
    #message(copenum)
    df.idx<-df.ex[df.ex$COPENUM==copenum,]
    featdir<-list.files(path = truerootdir,pattern = paste0("cope",paste0(copenum,"_"),".*randomize"),full.names = T)
    featdir<-featdir[-grep(".jpeg",featdir)]
    cmindx<-gen_cluster_mask(featdir=featdir,base=basemask,corrp_mask=corrp_mask,outdir = cmoutdir,VersionCode = Version,
                             maskthresholdvalue=corrmaskthreshold,roimaskthreshold=roimaskthreshold,useMMcor=useMMcor,
                             overwrite=!saveclustermap)
    cmindx<-cmindx[cmindx$Voxels>voxelnumthres,]
    
    clx<-clustertoget[[as.character(copenum)]]
    if(is.null(clx)){clx<-cmindx$`Cluster Index`}
    clx<-clx[clx %in% cmindx$`Cluster Index`]
    if (length(clx)>0) {
      concatimg<-list.files(pattern = ".*4D.nii.gz$",path = featdir,full.names = T)
      if(any(grepl("concat4D.nii.gz",concatimg)) & !forceconcat) {concatimg<-concatimg[grepl("concat4D.nii.gz",concatimg)]
      }else if(length(concatimg)<1 | grepl("PairedT",concatimg) | !is.na(grp_identif) | forceconcat){
        concatimg<-file.path(featdir,"concat4D.nii.gz")
        concatcmd<-paste(sep=" ","fslmerge -t",concatimg,paste(df.idx$PATH,collapse = " "))
        system(concatcmd,intern = F)
      } else {stop("unexpected error at imaging concating...")}
      roivalues<-as.data.frame(do.call(cbind,lapply(clx, function(clz){
        #print(clz)
        # roivalue<-sapply(1:length(df.idx$ID), function(iz){
        #   print(iz)
        # cmdx<-paste(sep=" ","fslstats",as.character(df.idx$PATH[iz]),
        #            "-k",cmindx$maskpath[cmindx$`Cluster Index`==clz],"-M")
        # system(cmdx,intern = T)})
        #Use fslstat timeserires to calculate;
        cmdx<-paste(sep=" ","fslstats -t",concatimg,
                    "-k",cmindx$maskpath[cmindx$`Cluster Index`==clz],"-M")
        as.numeric(system(cmdx,intern = T))
      })),stringsAsFactors = F)
      roivalues<-as.data.frame(lapply(roivalues,as.numeric),stringsAsFactors = F)
      names(roivalues)<-paste(paste0("cope",copenum),"cluster",clx,sep = "_")
      roivalues$ID<-df.idx$ID
      return(list(roivalues=roivalues,index=cmindx[cmindx$`Cluster Index`%in% clx,-grep("maskpath",names(cmindx))],corrthreshold=corrmaskthreshold))
    } else return(NULL)
  })
  parallel::stopCluster(sharedproc)
  names(copes_roivalues)<-paste("cope",copetoget,sep = "_")
  if(saverdata) {
    tempenvir<-as.environment(list())
    if(file.exists(file.path(truerootdir,"extracted_roi.rdata"))) {
      load(file = file.path(truerootdir,"extracted_roi.rdata"),envir = tempenvir)}
    assign(gsub("-","_",paste0(modelname,"_",corrp_mask,corrmaskthreshold,"_","roi")),copes_roivalues,envir = tempenvir)
    save(list=objects(tempenvir),
         file = file.path(truerootdir,"extracted_roi.rdata"),envir = tempenvir)
  }
  return(copes_roivalues)
}

#####Time series data extraction:
get_dim_single<-function(imagepath){
  fsl_2_sys_env()
  fslinfoout<-system(paste("${FSLDIR}/bin/fslinfo",imagepath,sep = " "),intern = T)
  kx<-lapply(strsplit(fslinfoout," "),function(x){x[x!=""]})
  rx<-as.data.frame(t(sapply(kx,`[[`,2)),stringsAsFactors = F)
  names(rx)<-sapply(kx, `[[`,1)
  return(rx)
}

###Get time serires from imaging data
get_timeserires_single<-function(imagepath=NULL,maskpath=NULL,outputpn=NULL,templatepath=NULL,subjectmask=NULL,notrans=F,output=T,forcererun=F){
  fsl_2_sys_env()
  rtp<-strsplit(imagepath,split = .Platform$file.sep)[[1]]
  imagepath_fd<-paste(rtp[1:(length(rtp)-1)],collapse = .Platform$file.sep)
  imagepath_fn<-rtp[end(rtp)[1]]
  if(system(paste("${FSLDIR}/bin/fslstats",maskpath,"-r",sep = " "),intern = T)!="0.000000 1.000000 "){stop("Mask is not binary. Can only use binary mask")}
  voxz_img<-unique(unlist(get_dim_single(imagepath)[c("pixdim1","pixdim2","pixdim3")]))
  dim_img<-unlist(get_dim_single(imagepath)[c("dim1","dim2","dim3")])
  
  voxz_mask<-unique(unlist(get_dim_single(maskpath)[c("pixdim1","pixdim2","pixdim3")]))
  dim_mask<-unlist(get_dim_single(maskpath)[c("dim1","dim2","dim3")])
  
  if(any(!dim_img == dim_mask)){
    if(notrans){stop("function terminated because mask and image has different dimension, and specified not to transform any imagery.")}
    if(is.null(templatepath)){stop("template not supplied. Will not be able to transform")}
    
    if(!file.exists(file.path(imagepath_fd,paste0("trans_",imagepath_fn))) | forcererun){
      
      voxz_tp<-unique(unlist(get_dim_single(templatepath)[c("pixdim1","pixdim2","pixdim3")]))
      dim_tp<-unlist(get_dim_single(templatepath)[c("dim1","dim2","dim3")])
      if(!any(!dim_tp == dim_mask)){
        message("mask and image has different dimension, trying to wrap image to standard space")
        if(is.null(subjectmask)){stop("Requires a binary subject mask.")}
        NRJ<-system(paste("${FSLDIR}/bin/flirt","-in",subjectmask,"-ref",templatepath,"-omat",file.path(imagepath_fd,"ts_transout.mat"),"-usesqform"),intern = T)
        NRK<-system(paste("${FSLDIR}/bin/flirt","-in",imagepath,"-ref",templatepath,"-out",file.path(imagepath_fd,paste0("trans_",imagepath_fn)),
                          "-init",file.path(imagepath_fd,"ts_transout.mat"),"-applyxfm"),intern = T)
        imagepath<-file.path(imagepath_fd,paste0("trans_",imagepath_fn))
      } else {
        stop("mask and image has different dimension and mask and standard space has different dimension, sort it out.")
      }
    }else {
      imagepath<-file.path(imagepath_fd,paste0("trans_",imagepath_fn))
    }
  }
  TNX<-system(paste("fslmeants","-i",imagepath,"-m",maskpath),intern = T)
  TNX<-as.numeric(TNX[TNX!=""])
  if(!is.null(outputpn)){writeLines(text = TNX,con=outputpn)}
  if(output){
    return(TNX)
  }
}

get_timeserires<-function(ssub_root=NULL,maskpath=NULL,templatepath=NULL,tarname=NULL,submaskname=NULL,parallen=4,type="preproc",depthcontr=3,forcerunmatch=F) {
  lx<-lapply(list.dirs(ssub_root,recursive = F,full.names = F),function(x){
    rx<-system(paste("find -L",file.path(ssub_root,x),"-name",tarname,"-maxdepth",depthcontr,"-type f"),intern = T)
    if(length(rx)>0){
      rxp<-strsplit(rx,split = .Platform$file.sep)
      lapply(1:length(rxp),function(rx){
        rtp<-rxp[[rx]]
        rx_fd<-paste(rtp[1:(length(rtp)-1)],collapse = .Platform$file.sep)
        return(list(root=rx_fd,imagepath=paste(rtp,collapse = .Platform$file.sep),ID=x,num=rx,
                    subjectmask=list.files(pattern = submaskname,path = rx_fd,full.names = T,no.. = T)))
      })
    
    }else{return(list(list(ID=x,num=0)))}
  })
  glx<-unlist(lx,recursive = F)
  names(glx) <- sapply(glx,function(gx){paste(gx$ID,gx$num,sep="_")})
  FUNX<-function(xrz){
    if(!is.null(xrz$imagepath)){
      return(get_timeserires_single(imagepath = xrz$imagepath,maskpath = maskpath,output = T,templatepath = templatepath,subjectmask = xrz$subjectmask,forcererun = forcerunmatch))
    }else{return(NULL)}
  }
  if(is.null(parallen)){
    RESULTX<-lapply(cl = clusterx,X = glx,fun = FUNX)
  } else {
    clusterx<-parallel::makeCluster(parallen,type = "FORK")
    RESULTX<-parallel::parLapply(cl = clusterx,X = glx,fun = FUNX)
    parallel::stopCluster(clusterx)
  }
  
  return(RESULTX)
}



deconv_timeseries<-function(datalist=NULL,tslist=NULL,func.proc=NULL,evtname=NULL,tr=NULL,num.calibrate=0.5,variname="ts_beta",func.deconv=mean){
  xz<-lapply(datalist,do.call,what=func.proc)
  if(!is.null(num.calibrate)){message("Calibration number is set to be ",num.calibrate,", this will be added to before and after each epoch to capture more volumes.")
    }else{num.calibrate=0}
  tr<-as.numeric(tr)
  rxz<-lapply(1:length(xz),function(ki){
    kx<-xz[[ki]]
    ky<-kx$event.list[[evtname]]
    ky$ID<-names(xz)[[ki]]
    kys<-split(ky,ky$run)
    })
  txz<-unlist(rxz,recursive = F,use.names = F)
  names(txz)<-sapply(txz,function(xi){paste0(unique(xi$ID),"_",unique(xi$run))})
  
  frx<-lapply(names(txz),function(rtx){
    #print(rtx)
    TTbyRun<-txz[[rtx]]
    TSbyRun<-tslist[[rtx]]
    if(is.null(TSbyRun)){return(NULL)}
    ts_df<-data.frame(tsbeta=TSbyRun,t_start=seq(from=0,to=(length(TSbyRun)-1)*tr,by = tr),stringsAsFactors = F)
    ts_df$t_end<-ts_df$t_start + tr
    ts_df$whichtrial<-NA
    
    for(ui in 1:nrow(TTbyRun)){
      uidx<-which(ts_df$t_start > TTbyRun$onset[ui]-num.calibrate & ts_df$t_end < TTbyRun$onset[ui] + TTbyRun$duration[ui] +num.calibrate)
      if(length(uidx)==0){message("Failed to caputre any volumes for subject ",rtx," and trial ",ui,". Consider increase the calibration number. Putting in NA for now")}
      if(any(!is.na(ts_df$whichtrial[uidx]))){message("subject ",rtx," and trial ",ui,": Trial event and clibration num might be algined in a way that overlaps two epoch and result in uninteneded mapping. Try reduce calibration number.")} 
      ts_df$whichtrial[uidx]<-ui
    }  
    ts_df_cl<-ts_df[!is.na(ts_df$whichtrial),]
    if(any(!TTbyRun$trial %in% ts_df_cl$whichtrial)){
      ts_df_cl<-rbind(ts_df_cl,
      data.frame(tsbeta=NA,t_start=NA,t_end=NA,whichtrial=TTbyRun$trial[!TTbyRun$trial %in% ts_df_cl$whichtrial],stringsAsFactors = F)
      )
      }
    TTbyRun[[variname]]<-sapply(lapply(split(ts_df_cl,ts_df_cl$whichtrial),function(r){r$tsbeta}),func.deconv)
    return(TTbyRun[c("ID","run","trial",variname)])
  })
  grx<-do.call(rbind,frx)
  return(split(grx,grx$ID))
}










