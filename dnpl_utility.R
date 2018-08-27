#######
##
##
devtools::source_url("https://raw.githubusercontent.com/DecisionNeurosciencePsychopathology/fMRI_R/master/fslpipe.R")
#devtools::source_url("https://raw.githubusercontent.com/DecisionNeurosciencePsychopathology/fMRI_R/master/prep_for_second_lvl.R")
##Here's all the functions that helps with the fsl pipe function;

cleanuplist<-function(listx){
  listx[sapply(listx, is.null)] <- NULL
  return(listx)}

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
                              returnas=c("path","data.frame"),
                              dothese=c("nuisance","motion")
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
      combo<-list()
      if ("nuisance" %in% dothese) {
      nui<-read.table(x$nuisance)
      names(nui)<-unlist(strsplit(cfg$preproc_call$nuisance_compute,split = ","))
      } else {nui<-data.frame()}
      if ("motion" %in% dothese) {
      mo<-read.table(x$motion)
      names(mo)<-paste0("motion_V",1:length(mo))
      } else {mo<-data.frame()}
      combo<-c(nui,mo)
      combo<-as.data.frame(combo)
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

#############Make reg for second level:
prepare4secondlvl<-function(ssana.path=NULL,    
                                #Single Sub Analysis Folder, should contain ID/*.feat
                                preproc.path=NULL,
                                #Preproc folder, should contain ID/${taskname}/${proc.name}[0-9]
                                standardbarin.path="/Volumes/bek/Newtemplate_may18/fsl_mni152/MNI152_T1_2mm_brain.nii",
                                #Path to standard template, take either .nii or .nii.gz
                                proc.name=NULL,
                                #See preproc path comment, this should always followed by a number 
                                taskname=NULL,
                                #See preproc path comment, this should be the dir that contains multiple runs
                                dir.filter=NULL,
                                #Filter out some of the unwanted folders in ssana.path, important especailly for those contain .feat structure
                                overwrite=TRUE,
                                #Logical Statement, if to overwirte existing link
                                outputmap=FALSE,
                                #Logical Statement, if true to return linkmap
                                paralleln=NULL
                                #Number of parallel jobs on going
) {
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
  #Step 5, make symb.link for each subject and each directory pair:
  for (i in 1:length(linkmap$id)) {
    if (overwrite) {
      st1<-suppressWarnings(system(paste0("rm -r ", linkmap$destination[i]),intern = T))
      #st1_5<-suppressWarnings(system(paste0("rm -r ", linkmap$destination_standard[i]),intern = T))
      st1_9<-suppressWarnings(system(paste0("rm -r ", linkmap$originplace[i]),intern = T))
    }
    st2<-dir.create(showWarnings = F,path = linkmap$destination[i])
    if (!file.exists(linkmap$originplace[i]))  {
      fsl_2_sys_env()
      system(paste0("${FSLDIR}/bin/flirt -in ",file.path(ssana.path,linkmap$id[i],linkmap$runword[i],"mask.nii.gz")," -ref ",
                    standardbarin.path," -out /Volumes/bek/neurofeedback/.temp.nii -omat ",linkmap$originplace[i],
                    " -bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12  -interp trilinear"))
    }
    
    #Make actual links now....
    #invisible(
    if (!file.exists(file.path(linkmap$destination,"example_func2standard.mat")[i])){
      file.symlink(from = file.path(linkmap$originplace)[i],to = file.path(linkmap$destination,"example_func2standard.mat")[i])
      file.symlink(from = file.path(linkmap$originplace)[i],to = file.path(linkmap$destination,"standard2example_func.mat")[i])
      file.symlink(from = standardbarin.path,to = file.path(linkmap$destination,"standard.nii.gz")[i])
    } else {message("meh,already there, if you want to overwirite, do overwrite...")}
  
    #)
  }} else {
    cluster_prep2ndlvl<-makeCluster(paralleln,outfile=file.path(argu$ssub_outputroot,argu$model.name,"prep42ndlvl_log.txt"),type = "FORK")
    NU<-parSapply(cluster_prep2ndlvl,1:length(linkmap$id),function(i) {
      fsl_2_sys_env()
      if (overwrite) {
        st1<-suppressWarnings(system(paste0("rm -r ", linkmap$destination[i]),intern = T))
        #st1_5<-suppressWarnings(system(paste0("rm -r ", linkmap$destination_standard[i]),intern = T))
        st1_9<-suppressWarnings(system(paste0("rm -r ", linkmap$originplace[i]),intern = T))
      }
      st2<-dir.create(showWarnings = F,path = linkmap$destination[i])
      if (!file.exists(linkmap$originplace[i]))  {
        system(paste0("${FSLDIR}/bin/flirt -in ",file.path(ssana.path,linkmap$id[i],linkmap$runword[i],"mask.nii.gz")," -ref ",
                      standardbarin.path," -out /Volumes/bek/neurofeedback/.temp.nii -omat ",linkmap$originplace[i],
                      " -bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12  -interp trilinear"))
      }
      
      #Make actual links now...
      if (!file.exists(file.path(linkmap$destination,"example_func2standard.mat")[i])){
        file.symlink(from = file.path(linkmap$originplace)[i],to = file.path(linkmap$destination,"example_func2standard.mat")[i])
        file.symlink(from = file.path(linkmap$originplace)[i],to = file.path(linkmap$destination,"standard2example_func.mat")[i])
        file.symlink(from = standardbarin.path,to = file.path(linkmap$destination,"standard.nii.gz")[i])
      } else {message("meh,already there, if you want to overwirite, do overwrite...")}
  
      })
  }
    
  if(outputmap) {return(linkmap)}
  print("DONE")
}

#######################################

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
  add.nuisa=TRUE) {
  
  #Read config file:
  cfg<-cfg_info(cfgpath)
  
  #Prep the data into generally acceptable output object;
  output<-do.call(what = do.prep.call,args = do.prep.arg)
  #assign("output",do.call(what = do.prep.call,args = do.prep.arg),envir=globalenv())
  
  dsgrid<-read.csv(gridpath,stringsAsFactors = F)
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
  
  return(design)
  
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


populate_within<-function(chunk_within=NULL,xvnum=NULL,variname=NULL){
  seq_lines<-lapply(1:xvnum,function(xjk){
      temp<-as.environment(list())
      assign(variname,xjk,envir = temp)
      result<-change_fsl_template(fsltemplate = chunk_within,begin = "XG_",end = "_EX",searchenvir = temp)
      rm(temp)
      return(result)
  })
  re_within<-do.call(c,seq_lines)  
  return(re_within)
}


adopt_gfeat<-function(adptemplate_path=NULL,searenvir=NULL) {
  readLines(adptemplate_path)->adptemplate
  #find the area that needs to be replaced:
  #For real tho, the order should be fine but just to be 1000000 % sure, let's match them
  while(length(grep("SEC_.*_START",adptemplate))>0) {
  startnum<-grep("SEC_.*_START",adptemplate)[1]
  name<-sapply(strsplit(adptemplate[startnum],split = "_"),"[[",2)
  endnum<-grep(paste0("SEC_",name,"_END"),adptemplate)
  #Let's not use lapply cuz combining them would be a headache..
    print(as.character(name))
    chunk_pre<-adptemplate[1:startnum-1]
    chunk_post<-adptemplate[(endnum+1):length(adptemplate)]
    #Okay now we constructure what's within.
    chunk_within<-adptemplate[(startnum+1):(endnum-1)]
    variname<-substr(chunk_within,regexpr(paste0("XG_","*"),chunk_within)+nchar("XG_"),
                      regexpr(paste0("*","_EX"),chunk_within)-1)
    unique(variname[variname!=""])->finavariname
    get(finavariname,envir = searenvir)->xvnum
    if (length(finavariname)==1) {
      new_within<-populate_within(chunk_within = chunk_within,xvnum = xvnum,variname = finavariname)
    } else {stop("Something Went Wrong! This is adopt_gfeat, first step")}
    if (length(new_within)>0) {
    adptemplate<-c(chunk_pre,new_within,chunk_post)
    } else {stop(paste0("Something Went Wrong! This is adopt_gfeat, 2nd step"))}
  } #End while 
  return(adptemplate)
}

###############################
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
                          fsfpath=NULL,
                          envir=NULL) {
  if (is.null(fsltemplate)) {fsltemplate<-readLines(templatepath)}
  subbyrunfeat<-change_fsl_template(fsltemplate = fsltemplate,begin = beg,end=end,searchenvir = envir)
  #fsfpath<-fsf.path
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
                         colour="red") {
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
                        modelname="PE_8C_old",
                        grp_sep=argu$group_id_sep,
                        copestorun=1:8,
                        paralleln=NULL,
                        onesamplet_pergroup=T,
                        pairedtest=T,
                        thresh_cluster_siz=3
) {
  if ( is.null(modelname) ) {stop("Must specify a model name other wise it will be hard to find all copes")}
  
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
    df.jx<-merge(df.ex,do.call(rbind,lapply(sortid(dix=file.path(rootdir,modelname),idgrep=grp_sep,dosymlink=F),function(x) {data.frame(ID=x$ID,GROUP=x$name)})),by = "ID",all.x = T)
  } else {df.ex$GROUP<-"ONE"
  df.jx<-df.ex}
  
  allcopecomx<-as.environment(list())
  #Now we face this problem of runing bunch of them...
  #One sample T test;
  if (length(unique(df.jx$GROUP))==1){ 
    #Make Commands;
    cope.fslmerge<-lapply(copestorun,function(x) {
      outputroot<-file.path(outputdir,modelname,paste0("cope",x,"_randomize_onesample_ttest"))
      dir.create(outputroot, showWarnings = FALSE,recursive = T)
      if (length(list.files(pattern = "*tfce_corrp_tstat1",path = outputroot,no.. = T))<1) {
      copefileconcat<-paste(df.jx$PATH[which(df.jx$COPENUM==x)],collapse = " ")
      return(paste0("fslmerge -t ",outputroot,"/OneSamp4D"," ",
             copefileconcat
             ," \n ",
             "randomise -i ",outputroot,"/OneSamp4D -o ",outputroot,"/OneSampT -1 -T -x -c ",thresh_cluster_siz
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
        if (length(list.files(pattern = "*tfce_corrp_tstat1",path = outputroot,no.. = T))<1) {
          copefileconcat<-paste(as.character(df.jx$PATH[which(df.jx$COPENUM==cope_group[1] & df.jx$GROUP==cope_group[2])]),collapse = " ")
          return(paste0("fslmerge -t ",outputroot,"/OneSamp4D"," ",
                        copefileconcat
                        ," \n ",
                        "randomise -i ",outputroot,"/OneSamp4D -o ",outputroot,"/OneSampT -1 -T -x -c ",thresh_cluster_siz
          ))
        }else {return(NULL)}
      })
      cleanuplist(cope.fslmerge)->cope.fslmerge
      assign(x = "onesamplet_pergroup",value = cope.fslmerge,envir = allcopecomx)
    } 
    if (pairedtest) {
      dir.create(file.path(outputdir,modelname),showWarnings = F,recursive = T)
      commonid<-prep_paired_t(idsep = sortid(dix=file.path(rootdir,modelname),idgrep=grp_sep,dosymlink=F), outpath = file.path(outputdir,modelname))
      cidindex<-data.frame(ID=do.call(c,lapply(commonid,function(xm) {paste(xm,grp_sep,sep = "_")})),notag=do.call(c,lapply(commonid,function(xm) {paste(xm,c("",""),sep = "")})))
      df.kh<-merge(df.jx,cidindex,all = T)
      df.jk<-df.kh[which(!is.na(df.kh$notag)),]
      cope.fslmerge<-lapply(copestorun,function(x) {
        outputroot<-file.path(outputdir,modelname,paste0("cope",x,"_randomize_paired_ttest"))
        list.files(pattern = "*tfce_corrp_tstat1",path = outputroot)
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
                        file.path(outputdir,modelname),"/design.grp -1 -T -x -c ",thresh_cluster_siz
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
    cj1<-makeCluster(paralleln,outfile=file.path(outputdir,modelname,"glvl_log.txt"),type = "FORK")
    NU<-parSapply(cj1,cope.fslmerge,function(x) {
      message(paste0("Now running ",x))
      system(command = x,intern = T,ignore.stdout = F,ignore.stderr = F)
      message("completed")
    })
    stopCluster(cj1)
  } else {
    lapply(cope.fslmerge,function(x) {
      message(paste0("Now running ",x))
      system(command = x,intern = T,ignore.stdout = F,ignore.stderr = F)
      message("completed")
    })
  }
  })
  message("DONE")
}


sortid<-function(dix=file.path(argu$ssub_outputroot,argu$model.name),idgrep=argu$group_id_sep,dosymlink=ifdosymlink){
  system(paste0("find ",dix," -iname 'average.gfeat' -maxdepth 2"),intern = T)->dirxs
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
    diag(x = 1,nrow = 1, ncol = length(commonid)+1),file = file.path(outpath,"design.con.txt"),row.names = F,col.names = F)
  system(paste0("${FSLDIR}/bin/Text2Vest ",file.path(outpath,"design.con.txt")," ",file.path(outpath,"design.con")))
  file.remove(file.path(outpath,"design.con.txt"))
  
  write.table(do.call(rbind, lapply(1:length(idsep), function(x){
    matrix(nrow = length(commonid),ncol = 1,data = seq_along(commonid))
  })),file = file.path(outpath,"design.grp.txt"),row.names = F,col.names = F)
  system(paste0("${FSLDIR}/bin/Text2Vest ",file.path(outpath,"design.grp.txt")," ",file.path(outpath,"design.grp")))
  file.remove(file.path(outpath,"design.grp.txt"))
  
  return(commonid)
  
}




