######FSL PIPE:
#####Upcoming: Level 1 re-make for better efficiency


fsl_pipe<-function(argu=NULL, #This is the arguments environment, each model should have a different one;
                   prep.call.func="", #This should be a character string that's the name of the prep proc function
                   prep.call.allsub=list(ID=list(data=NULL)),#List of ID list of arguments for prep.call.
                   initonly=F
) {
  
  require("devtools")
  if("dependlab" %in% installed.packages()){
    system(paste("echo","GREAT, DEPENDLAB PACK IS INSTALLED"))
  }else{
    devtools::install_github("PennStateDEPENdLab/dependlab")
  }
  require("parallel")
  
  fsl_2_sys_env(force = T)
  
  #############STEP 0: Initializing argu #########
  #argu$cfg<-cfg_info(cfgpath = argu$cfgpath)
  argu$dsgrid<-read.table(argu$gridpath,header = T,sep = c(","),stringsAsFactors = F,strip.white = T,skipNul = T)
  if(is.null(argu$dsgrid$AddNeg)){argu$dsgrid$AddNeg<-FALSE}
  if(is.null(argu$dsgrid$RunGrpLvl)){argu$dsgrid$RunGrpLvl<-TRUE}
  argu$dsgrid$AddNeg<-as.logical(argu$dsgrid$AddNeg)
  argu$dsgrid$RunGrpLvl <- as.logical(argu$dsgrid$RunGrpLvl)
  # if(is.null(argu$model.varinames)) {argu$model.varinames<-argu$dsgrid$name}
  argu$home_dir <- system("echo ~",intern = T)
  ###Getting Paths:
  if (!is.null(argu$rootpath_output)) {
    argu$rootpath_output <- gsub("~",argu$home_dir,argu$rootpath_output)
    argu$lvl1path_output <- file.path(argu$rootpath_output,"single_subject")
    argu$lvl1path_reg    <- file.path(argu$rootpath_output,"reg")
    argu$lvl3path_output <- file.path(argu$rootpath_output,"group_level")
  }
  
  ###Version upgrade safe keeping
  
  #Replacing old variables names from previous versions
  replaceLS<-list(ifnuisa="convlv_nuisa",onlyrun="run_steps",centerscaleall="lvl1_centervalues",lvl1_run_on_pbs="run_on_pbs",
                  regpath = "lvl1path_reg",ssub_outputroot="lvl1path_output",subj_outputroot="lvl1path_output",
                  glvl_outputroot = "lvl3path_output",glvl_output = "lvl3path_output", func.nii.name = "name_func_nii",
                  model.name="model_name",templatedir="templatebrain_path")
  
  for(og in names(replaceLS)){
    if (exists(og,envir = argu) & !exists(replaceLS[[og]],envir = argu)) {
      message(og," variable is now depreciated, please use ",replaceLS[[og]],".")
      argu[[replaceLS[[og]]]]<-argu[[og]]}
  }
  
  if (!exists("lvl1_cmat",envir = argu)) {
    message("Single subject contrast matrix is not specified, will use automatic generated one by using grid.")
    ogLength<-length(argu$dsgrid$name)
    negNum<-(which(argu$dsgrid$AddNeg))
    if(length(negNum)>0){
      negMat<-diag(x=-1,ogLength)
      negMat<-negMat[negNum,]
      argu$lvl1_cmat<-rbind(diag(x=1,ogLength),negMat)
      colnames(argu$lvl1_cmat)<-argu$dsgrid$name
      rownames(argu$lvl1_cmat)<-c(argu$dsgrid$name,paste(argu$dsgrid$name[negNum],"neg",sep = "_"))
    } else {
      argu$lvl1_cmat<-diag(x=1,ogLength)
      rownames(argu$lvl1_cmat)<-argu$dsgrid$name
      colnames(argu$lvl1_cmat)<-argu$dsgrid$name
    }
  } 
  
  default_ls<-list(lvl2_prep_refit=FALSE,lvl1_centervalues=FALSE,run_on_pbs=FALSE,lvl1_centervalues=TRUE,lvl1_forcegenreg=FALSE,
                   qsub_limits=20,lvl2_force_prep=FALSE,lvl1_retry=FALSE,lvl1_afnify=F,lvl2_afnify=F,lvl3_afnify=T,
                   nuisa_motion=c("nuisance","motion_par"),lvl3_lowlvlfeatreg="average.gfeat",motion_type="fd",
                   motion_threshold="default",job_per_qsub=as.numeric(argu$cfg$n_expected_funcruns),run_pipeline=TRUE,
                   lvl3_type="flame",adaptive_ssfeat=TRUE)
  default_ls<-default_ls[!names(default_ls) %in% names(argu)]
  if (length(default_ls)>0){
    for(lx in 1:length(default_ls)){
      message("Variable: '",names(default_ls)[lx],"' is not set, will use default value: ",default_ls[[lx]])
    }
    argu<-list2env(default_ls,envir = argu)
  }
  
  ####SYSTEM ENV LS
  ###Get from the system environment, placed after default so the default value can be taken care of;
  sysenvLS<-c(run_pipeline=TRUE,run_steps=NULL)
  for (sysenvX in sysenvLS) {
    if(!exists("run_pipeline",envir = argu)){
      if(Sys.getenv("run_pipeline")==""){argu$run_pipeline<-TRUE}else{argu$run_pipeline<-Sys.getenv("run_pipeline")}
    }
    if(!exists("run_steps",envir = argu)){
      if(Sys.getenv("run_steps")==""){argu$run_steps<-NULL}else{argu$run_steps<-Sys.getenv("run_steps")}
    }
  } 
  
  ##Parallel processes number
  if (is.null(argu$nprocess)){
    if (detectCores()>12){
      num_cores<-12 
    } else {num_cores<-detectCores()-2} 
  } else {argu$nprocess->num_cores}
  
  #Renaming;
  if(argu$adaptive_ssfeat){argu$ssub_fsl_templatepath<-system.file("extdata", "fsl_ssfeat_general_adaptive_template_R.fsf", package="fslpipe")}
  if(is.null(argu$tr)) {argu$tr <- argu$cfg$preproc_call$tr}
  if(!argu$run_pipeline){return(NULL)}
  if(initonly) {return(argu)}
  
  ###Construct project level information from cfg:
  dir.create(file.path(argu$lvl1path_output,argu$model_name,"misc_info"),showWarnings = F,recursive = T)
  
  
  if(is.null(argu$lvl1_volinfo)) {
    message("####!!!!No project configuration object found in argu environment. Will use default one for DependLab Option!!!####")
    argu$lvl1_infodf <- gen_project_config_wCFG(cfg = argu$cfg,
                                                bID_array = names(prep.call.allsub),
                                                input_nii_pattern = argu$name_func_nii,
                                                add_nuisance = argu$lvl1_proc_nuisance)
    write.csv(argu$lvl1_infodf,file = file.path(argu$lvl1path_output,argu$model_name,"misc_info","step_0_info.csv"),row.names = F)
    argu$lvl1_volinfo <- argu$lvl1_infodf[names(argu$lvl1_infodf) %in% c("ID","behavioral_data") | grepl("path",x = names(argu$lvl1_infodf))]
    argu$lvl1_volinfo <- reshape2::melt(argu$lvl1_volinfo,id.vars=c("ID","behavioral_data"))
    argu$lvl1_volinfo$variable <-gsub("path_","",argu$lvl1_volinfo$variable)
    names(argu$lvl1_volinfo) <- c("ID","behavioral_data","run","path")
  }
  
  

  
  #Get preproc info:
  #############STEP 1: Regressor generation####################
  #GENERATE REGRESSOR USING DEPENDLAB PIPELINE:
  stepnow<-1
  if (is.null(argu$run_steps) | stepnow %in% argu$run_steps) {
    #Create the directory if not existed
    #print(argu$lvl1path_output, argu$model_name)
    argu$lvl1_datalist<-prep.call.allsub
    argu$lvl1_procfunc<-get(prep.call.func)
    argu$lvl1_volinfo$vol <- get_volume_run2(paths = argu$lvl1_volinfo$path)
    #argu$lvl1_datalist<-argu$lvl1_datalist[which(names(argu$lvl1_datalist) %in% IDTORUN)]
    dir.create(file.path(argu$lvl1path_output,argu$model_name),showWarnings = FALSE,recursive = T)
    # #load the design rdata file if exists;
    # lvl1_datalist=argu$lvl1_datalist
    # lvl1_proc_func = argu$lvl1_procfunc
    # lvl1_volinfo = argu$lvl1_volinfo
    # forcererun = argu$lvl1_forcegenreg
    # retry=argu$lvl1_retry
    # model_name=argu$model_name
    # dsgrid=argu$dsgrid
    # center_values=argu$lvl1_centervalues
    # nprocess=argu$nprocess
    # reg_rootpath=argu$lvl1path_reg
    # tr = argu$tr
    
    step1_cmd<-substitute({
      allsub_design<-do_all_first_level(lvl1_datalist=argu$lvl1_datalist,lvl1_proc_func = argu$lvl1_procfunc,lvl1_volinfo = argu$lvl1_volinfo,
                                        forcererun = argu$lvl1_forcegenreg,retry=argu$lvl1_retry,
                                        model_name=argu$model_name,dsgrid=argu$dsgrid,center_values=argu$lvl1_centervalues,
                                        reg_rootpath=argu$lvl1path_reg,nprocess=argu$nprocess) 
      save(allsub_design,file = file.path(argu$lvl1path_output,argu$model_name,"design.rdata"))
    })

    if(argu$run_on_pbs){
      workingdir<-file.path(argu$lvl1path_output,argu$model_name,"lvl1_misc")
      dir.create(workingdir,showWarnings = F,recursive = F)
      setwd(workingdir)
      save(list = ls(),file = "curwd.rdata")
      writeLines("library(fslpipe);load(\"curwd.rdata\");eval(step1_cmd)",con = "temp.r")
      pbs_torun<-get_pbs_default();pbs_torun$cmd<-"Rscript temp.r";pbs_torun$ppn<-argu$nprocess
      writeLines(do.call(pbs_cmd,pbs_torun),"pbs_temp.sh")
      dependlab::wait_for_job(dependlab::qsub_file("pbs_temp.sh"))
      load(file.path(argu$lvl1path_output,argu$model_name,"design.rdata"))
    } else {
      eval(step1_cmd)
    }
    
    #dfc <- read.csv(file.path(argu$lvl1path_output,argu$model_name,"misc_info","step_0_info.csv"),stringsAsFactors = F)
    #dfe <- merge(dfc,data.frame(ID=names(allsub_design)[!sapply(allsub_design,is.null)],gen_reg = T,stringsAsFactors = F),by = "ID",all = T)
    #write.csv(dfe,file = file.path(argu$lvl1path_output,argu$model_name,"misc_info","step_1_info.csv"),row.names = F)
    
    message("Step ", stepnow ," Ended")
  }
  
  #############Step 2: LVL1: Single Subject PARALLEL##########
  #Now we do the single sub processing using FSL and the regressor that was generated
  stepnow<-2
  if (is.null(argu$run_steps) | stepnow %in% argu$run_steps) {
    
    
    if (file.exists(file.path(argu$lvl1path_output,argu$model_name,"design.rdata"))) {
      load(file.path(argu$lvl1path_output,argu$model_name,"design.rdata"))
    } else {stop("No design rdata file found, must re-run step 1")}
    
    
    #let's subset this 
    small.sub<-lapply(as.list(allsub_design), function(x) {
      list(
        design=x$design,
        ID=x$ID,
        run_volumes=x$design$run_volumes,
        lvl1path_reg=x$lvl1path_reg,
        preprocID=x$preprocID)
    })
    
    #small.sub <- small.sub[which(sapply(small.sub,`[[`,'ID') %in% IDTORUN)]
    #IT's the same for all participants!!!! WHY DO YOU RE RUN IT FOR EVERYONE!!!
    xarg<-as.environment(list())
    xarg$templatebrain<-argu$templatedir
    xarg$tr<-argu$tr
    
    if(argu$adaptive_ssfeat){
      argu$model_varinames<-argu$dsgrid$name
      argu$copenames<-rownames(argu$lvl1_cmat)
      xarg$evnum<-1:ncol(argu$lvl1_cmat)
      xarg$copenum<-1:nrow(argu$lvl1_cmat)
      xarg$maxevnum<-ncol(argu$lvl1_cmat)
      xarg$maxcopenum<-nrow(argu$lvl1_cmat)
      argu$maxcopenum<-nrow(argu$lvl1_cmat)
      
      for (xy in 1:nrow(argu$lvl1_cmat)){
        assign(paste0("copemat",xy),value = 1:ncol(argu$lvl1_cmat),envir = xarg)
        assign(paste0("copetitle",xy),argu$copenames[xy],envir = xarg)
        assign(paste0("cope_lessnum",xy),(1:length(argu$copenames))[-xy],envir = xarg)
        for(xx in 1:ncol(argu$lvl1_cmat)){
          assign(paste0("copevalue",xy,"_",xx),value = argu$lvl1_cmat[xy,xx],envir = xarg)
        }
      } 
      for (mv in 1:ncol(argu$lvl1_cmat)) {
        assign(paste0("evtitle",mv),argu$model_varinames[mv],envir = xarg)
        assign(paste0("orthocombo",mv),paste(mv,(0:length(argu$model_varinames)),sep = "."),envir = xarg)
      }
      #PUT NEW FUNCTION HERE
      ssfsltemp<-rep_within(adptemplate = readLines(argu$ssub_fsl_templatepath),searchenvir = xarg)
    } else {ssfsltemp<-readLines(argu$ssub_fsl_templatepath)}
    
    step2commands<-unlist(lapply(small.sub,function(x) {
      idx<-x$ID
      aarg<-xarg
      cmmd<-unlist(lapply(1:length(x$run_volumes), function(runnum) {
        vol_info_rlvl <- argu$lvl1_volinfo[which(argu$lvl1_volinfo$ID==idx & argu$lvl1_volinfo$run == paste0("run",runnum)),]
        if(nrow(vol_info_rlvl)!=1) {
          message("Unable to initialize participant: ",idx,", and run: ",runnum," because lvl1 volume information data frame has incorrect number of rows: ",nrow(vol_info_rlvl))
          return(NULL)  
        }
        if(is.na(vol_info_rlvl$path)) {
          message(paste0("No input data file for participant: ",idx,", and run: ",runnum))
          return(NULL)
        }
        aarg$outputpath<-file.path(argu$lvl1path_output,argu$model_name,idx,paste0("run",runnum,"_output"))
        if (!file.exists(file.path(paste0(aarg$outputpath,".feat"),"stats","zstat1.nii.gz")) ) {
          message(paste0("Initializing feat for participant: ",idx,", and run: ",runnum))
          if(dir.exists(file.path(paste0(aarg$outputpath,".feat"))) ){
            message("Found an incomplete folder...Removing...")
            unlink(file.path(paste0(aarg$outputpath,".feat")),recursive = T,force = T)}
          if(is.null(argu$ss_zthreshold)) {aarg$zthreshold<-3.2} else {aarg$zthreshold<-argu$ss_zthreshold}
          if(is.null(argu$ss_pthreshold)) {aarg$pthreshold<-0.05} else {aarg$pthreshold<-argu$ss_pthreshold}
          
          aarg$runnum <- runnum   
          aarg$volumes <- vol_info_rlvl$vol
          aarg$funcfile <- vol_info_rlvl$path
          aarg$ifnuisa <- as.numeric(!is.null(vol_info_rlvl$nuisance) && !is.na(vol_info_rlvl$nuisance))
          if(as.logical(aarg$ifnuisa)) {
            aarg$nuisa<-vol_info_rlvl$nuisance
          } else {
            aarg$nuisa<-""
          }
          
          
          if(argu$adaptive_ssfeat){
            for (mv in 1:ncol(argu$lvl1_cmat)) {
              assign(paste0("evreg",mv),file.path(file.path(argu$lvl1path_reg,argu$model_name),idx,paste0("run",runnum,"_",argu$model_varinames[mv],".1D")),envir = aarg)
            }
          } else {
            gen_reg(vmodel=argu$model_varinames,lvl1path_reg=file.path(argu$lvl1path_reg,argu$model_name),idx=idx,runnum=runnum,env=aarg,regtype=argu$regtype)
          }
          if (any(unlist(eapply(aarg,is.na)))) {stop("NA exists in one of the arguments; please double check!")}
          #gen_reg(vmodel=argu$model.varinames,lvl1path_reg=file.path(argu$lvl1path_reg,argu$model_name),idx=idx,runnum=runnum,env=xarg,regtype = argu$regtype)
          
          cmmd<-feat_w_template(fsltemplate = ssfsltemp,beg = "ARG_",end = "_END",
                                fsfpath = file.path(argu$lvl1path_reg,argu$model_name,idx,paste0("run",runnum,"_",argu$model_name,".fsf")),
                                envir = aarg,outcommand = T)
          rm(aarg)
          return(cmmd)
        } else {
          message(paste("ID:",idx,"RUN:",runnum,",already exists,","to re-run, remove the directory."))
          return(NULL)}
      }))
      return(cmmd)
    }))
    
    if(length(step2commands)>0){
      if(argu$run_on_pbs){
        #PBS
        lvl1_workingdir<-file.path(argu$lvl1path_output,argu$model_name,"lvl1_misc",paste0(gsub(":","",gsub("-","_",gsub(pattern = " ","_",Sys.time()))),"log"))
        qsub_commands(cmds = step2commands,jobperqsub = argu$job_per_qsub,workingdir = lvl1_workingdir,
                      tagname = "lvl1",ppn = 4,qsublimit = argu$qsub_limits)
        
      }else{
        #run localy
        cluster_step2<-makeCluster(num_cores,outfile="",type = "FORK")
        NX<-parSapply(cluster_step2,step2commands,function(yx) {
          fsl_2_sys_env()
          message(paste0("starting to run /n ",yx))
          tryCatch(
            {system(command = yx,intern = T)
              pb<-txtProgressBar(min = 0,max = 100,char = "|",width = 50,style = 3)
              numdx<-which(yx==step2commands)
              indx<-suppressWarnings(split(1:length(step2commands),1:num_cores))
              pindx<-grep(paste0("\\b",numdx,"\\b"),indx)
              setTxtProgressBar(pb,(which(numdx==indx[[pindx]]) / length(indx[[pindx]]))*100)
              message("DONE")
            }, error=function(e){stop(paste0("feat unsuccessful...error: ", e))}
          )
          
        })
        stopCluster(cluster_step2)
        
      }
    } else {message("Nothing to run on lvl 1.")}
    
    ###Generate Table Output for Processing Status:
    #dfe <- read.csv(file.path(argu$lvl1path_output,argu$model_name,"misc_info","step_1_info.csv"),stringsAsFactors = F)
    #lvl2_featlist<-system(paste0("find ",file.path(argu$lvl1path_output,argu$model_name)," -iname ","*output.feat"," -maxdepth 2 -mindepth 1 -type d"),intern = T)
    #ncount<-sapply(split(basename(dirname(lvl2_featlist)),basename(dirname(lvl2_featlist))),length)
    
    #dff <- merge(dfe,data.frame(ID=names(ncount),nruns=ncount,stringsAsFactors = F),by = "ID",all = T)
    #write.csv(dfe,file = file.path(argu$lvl1path_output,argu$model_name,"misc_info","step__info.csv"),row.names = F)
    
    message("Step ", stepnow ," Ended")
  }
  
  #############Step 3: Prep for Higher Level #######################
  #Now we make the symbolic link for template matching...so they are not misaligned anymore...
  stepnow<-3
  if (is.null(argu$run_steps) | stepnow %in% argu$run_steps) {
    if(argu$lvl2_prep_refit){
      #cfg<-cfg_info(cfgpath = argu$cfgpath)
      gtax<-prepare4secondlvl(
        ssana.path=file.path(argu$lvl1path_output,argu$model_name),            
        standardbarin.path=argu$templatedir, featfoldername = "*output.feat",
        dir.filter=argu$hig_lvl_path_filter,                                                
        overwrite=argu$lvl2_force_prep,
        outputmap=TRUE,
        paralleln = num_cores)           
    } else {
      lvl2_featlist<-prep_session_lvl(subj_rootpath = file.path(argu$lvl1path_output,argu$model_name),subj_folderreg = "*output.feat",overwrite = argu$lvl2_force_prep,
                                      template_brainpath = argu$templatebrain_path)
    }
    message("Step ", stepnow ," Ended")
  }
  
  #############Step 4: LVL2: Fixed Effect for Single Subject PARALLEL ###############
  #This starts averaging for each subject:
  stepnow<-4
  if (is.null(argu$run_steps) | stepnow %in% argu$run_steps) {
    
    lowlvl_featreg = "*output.feat"
    if(!exists("lvl2_featlist")) {
      lvl2_featlist<-system(paste0("find ",file.path(argu$lvl1path_output,argu$model_name)," -iname ",lowlvl_featreg," -maxdepth 2 -mindepth 1 -type d"),intern = T)
    }
    raw.split <- strsplit(lvl2_featlist,split = .Platform$file.sep) 
    lvl2_rawdf<-do.call(rbind,lapply(raw.split,function(x){
      data.frame(
        uID=x[grep(lowlvl_featreg,x)-1],
        RUNNUM = gsub("run","",gsub("_output.feat","",x[grep(lowlvl_featreg,x)],fixed = T)),
        ID=paste(x[grep(lowlvl_featreg,x)-1],gsub("run","",gsub("_output.feat","",x[grep(lowlvl_featreg,x)],fixed = T)),sep = "_"),
        PATH = paste(x,collapse = .Platform$file.sep),NAME="average",
        OUTPUTPATH = dirname(paste(x,collapse = .Platform$file.sep)),
        SUBJMEAN=1,
        stringsAsFactors = F)
    })
    )
    
    
    if(!is.null(argu$lvl2_limitrun)){
      lvl2_rawdf<-lvl2_rawdf[which(as.character(lvl2_rawdf$RUNNUM) %in% as.character(argu$lvl2_limitrun)),]
    }
    
    lvl2_raw_sp<-split(lvl2_rawdf,lvl2_rawdf$uID)
    lvl2_default<-list(flame_type = 3, #FLAME 1 for sess_lvl and FLAME 1+2 for grp_lvl
                       thresh_type = 3, #0 : None \n # 1 : Uncorrected \n# 2 : Voxel \n # 3 : Cluster \n"
                       z_thresh = 2.3, #1.96 for subj lvl, 2.3 for sess lvl and 3.09 for grp lvl
                       p_thresh = 0.05, #0.05 default for both;
                       overwrite = F,Pairred_Group=FALSE,
                       covariate_names=c("SUBJMEAN")
    )
    lvl2_arg<-lapply(names(lvl2_default),function(xa){
      if(exists(paste0("lvl2_",xa),envir = argu)){get(paste0("lvl2_",xa),envir = argu)}else{lvl2_default[[xa]]}
    })
    names(lvl2_arg) <- names(lvl2_default)
    lvl2_arg$proc_ls_fsf <- lvl2_raw_sp 
    lvl2_arg$template_brain <- argu$templatebrain_path
    lvl2_arg$fsltemplate <- readLines(system.file("extdata", "fsl_flame_general_adaptive_template.fsf", package="fslpipe"))
    lvl2_alldf <- do.call(gen_fsf_highlvl,lvl2_arg)
    
    #lvl2_alldf[which(lvl2_alldf$uID %in% IDTORUN),]
    
    if(nrow(lvl2_alldf)>0){
      if(argu$run_on_pbs){
        #PBS
        lvl2_workingdir<-file.path(argu$lvl1path_output,argu$model_name,"lvl2_misc",paste0(gsub(":","",gsub("-","_",gsub(pattern = " ","_",Sys.time()))),"log"))
        pbs_args <- get_pbs_default(); pbs_args$ppn<-4; pbs_args$walltime="40:00:00";
        qsub_commands(cmds = paste("feat",unique(lvl2_alldf$FSF_PATH)),jobperqsub = argu$job_per_qsub,pbs_args=pbs_args,
                      workingdir = lvl2_workingdir,tagname = "lvl2",qsublimit = argu$qsub_limits)
        
      } else {
        lvl2_cluster<-parallel::makeCluster(argu$nprocess,outfile="",type = "FORK")
        NU<-parallel::parSapply(lvl2_cluster,unique(lvl2_alldf$FSF_PATH), function(y) {
          fsl_2_sys_env()
          message("starting feat /n ",y)
          tryCatch(
            {system(command = paste("feat",y,sep = " "),intern = T)
              message("DONE")
            }, error=function(e){stop(paste0("feat unsuccessful...error: ", e))}
          )
          
        })
        parallel::stopCluster(lvl2_cluster)
      }
    }
    
    message("Step ", stepnow ," Ended")
  }
  
  #############Step 5: LVL3: Higher Level (Randomize/FLAME) ##PARALLEL by function#########
  
  stepnow<-5
  if (is.null(argu$run_steps) | stepnow %in% argu$run_steps) {
    ssfsltemp<-readLines(argu$ssub_fsl_templatepath)
    
    #Randomize here:
    if(argu$lvl3_type=="randomize"){
      onesamplet_pergroup<-F
      pairedtest<-F
      unpairedtest<-F
      if (is.null(argu$group_id_sep) | !exists('group_id_sep',envir = argu)) {argu$group_id_sep<-""} 
      if (is.null(argu$cluster_thresh) | !exists('cluster_thresh',envir = argu)) {argu$cluster_thresh<-3} 
      if (is.null(argu$whichttest) | !exists('whichttest',envir = argu)) {argu$whichttest<-"onesample"}
      if (exists("supplyidmap",envir = argu)) {unpairedtest<-T}
      if (exists('group_id_sep',envir = argu) & "unpaired" %in% argu$whichttest) {unpairedtest<-T} 
      if ("onesample" %in% argu$whichttest) {onesamplet_pergroup<-T}
      if ("paired" %in% argu$whichttest) {pairedtest<-T}
      #To adopt the new chnages made in adaptive ss 
      if(argu$adaptive_ssfeat) {maxcopenum<-1:nrow(argu$lvl1_cmat)} else {
        maxcopenum<-1:max(as.numeric(gsub(".*?([0-9]+).*", "\\1", ssfsltemp[grep("# Title for contrast",ssfsltemp)])))
      }
      #Start Group Level Analysis:
      glvl_all_cope(rootdir=argu$lvl1path_output,
                    outputdir=argu$lvl3path_output,
                    modelname=argu$model_name,
                    grp_sep=argu$group_id_sep,
                    onesamplet_pergroup=onesamplet_pergroup,
                    pairedtest=pairedtest,
                    copestorun=maxcopenum,
                    thresh_cluster_mass=argu$thresh_cluster_mass,
                    thresh_cluster_extent=argu$thresh_cluster_extent,
                    pvalue=argu$randomize_p_threshold,unpairedtest=unpairedtest,
                    usethesetest=argu$randomize_thresholdingways,
                    ifDeMean=argu$randomize_demean,
                    paralleln = num_cores)
      # Use for debugging:
      # rootdir=argu$lvl1path_output
      # outputdir=argu$lvl3path_output
      # modelname=argu$model_name
      # grp_sep=argu$group_id_sep
      # onesamplet_pergroup=onesamplet_pergroup
      # pairedtest=pairedtest
      # thresh_cluster_siz=argu$cluster_thresh
      # copestorun = maxcopenum
      # thresh_cluster_mass=argu$thresh_cluster_mass
      # thresh_cluster_extent=argu$thresh_cluster_extent
      # pvalue=argu$randomize_p_threshold
      # usethesetest=argu$randomize_thresholdingways
      # ifDeMean=argu$randomize_demean
      # paralleln = num_cores
    } else if(tolower(argu$lvl3_type)=="flame") {
      #Run flame here:
      
      lowlvl_featreg<-argu$lvl3_lowlvlfeatreg
      
      raw<-system(paste0("find ",
                         file.path(argu$lvl1path_output,argu$model_name,"*/",lowlvl_featreg),
                         " -iname '*.feat' -maxdepth 2 -mindepth 1 -type d"),intern = T)
      raw.split <- strsplit(raw,split = .Platform$file.sep)  
      lvl3_rawdf<-do.call(rbind,lapply(raw.split,function(x){
        data.frame(
          ID=x[grep(lowlvl_featreg,x)-1],
          COPENUM=gsub("cope","",gsub(".feat","",dplyr::last(x))),
          PATH = paste(x,collapse = .Platform$file.sep),
          Intercept=1,
          stringsAsFactors = F)
      })
      )
      # NAME = ?
      # OUTPUTPATH = ?
      
      if(!is.null(argu$run_these_ID)) {
        lvl3_rawdf <- lvl3_rawdf[which(lvl3_rawdf$ID %in% argu$run_these_ID),]
      }
      
      if(!is.null(argu$exclude_these_ID)) {
        message("The following IDs are excluded from the final analysis: ",paste(argu$exclude_these_ID,collapse = ", "))
        lvl3_rawdf <- lvl3_rawdf[which(!lvl3_rawdf$ID %in% argu$exclude_these_ID),]
      }
      
      
      lvl3_raw_sp<-split(lvl3_rawdf,lvl3_rawdf$COPENUM)
      #lvl3_raw_sp<-lvl3_raw_sp[which(argu$dsgrid$RunGrpLvl)]
      
      lvl3_raw_sp<-lapply(lvl3_raw_sp,function(x){
        if(!is.null(argu$lvl3_ref_df)){x<-merge(x,argu$lvl3_ref_df,all.x=T,by="ID")}
        x$NAME = unique(readLines(file.path(x$PATH[1],"design.lev")))
        x$OUTPUTPATH = file.path(argu$lvl3path_output,paste(argu$model_name,paste(argu$lvl3_covariate_names,collapse = "_"),sep = "_"))
        return(x)
      })
      
      lvl3_default<-list(flame_type = 1, #FLAME 1 for sess_lvl and FLAME 1+2 for grp_lvl
                         thresh_type = 3, #0 : None \n # 1 : Uncorrected \n# 2 : Voxel \n # 3 : Cluster \n"
                         z_thresh = 3.09, #1.96 for subj lvl, 2.3 for sess lvl and 3.09 for grp lvl
                         p_thresh = 0.05, #0.05 default for both;
                         overwrite = F,Pairred_Group=FALSE,custom_evmat=NULL,custom_ctmat=NULL,
                         covariate_names=c("Intercept")
      )
      
      # lvl3_reg_default<-list(reg2main=0,reg2initial=0,reg2standard=1)
      
      lvl3_arg<-lapply(names(lvl3_default),function(xa){
        if(exists(paste0("lvl3_",xa),envir = argu)){get(paste0("lvl3_",xa),envir = argu)}else{lvl3_default[[xa]]}
      })
      names(lvl3_arg) <- names(lvl3_default)
      lvl3_arg$proc_ls_fsf <- lvl3_raw_sp 
      lvl3_arg$template_brain <- argu$templatebrain_path
      lvl3_arg$fsltemplate <- readLines(system.file("extdata", "fsl_flame_general_adaptive_template.fsf", package="fslpipe"))
      #lvl3_arg$covariate_names<-argu$lvl3_covarnames
      
      lvl3_alldf <- do.call(gen_fsf_highlvl,lvl3_arg)
      
      #Now there isn't a good way to deal with this problem, we have to have two arguments;
      #Sometimes two sample has to be combined but ran the same group level, but needs exclusion
      #Sometimes it's easy to include only certarin people (HC only)
      
      
      save(lvl3_alldf,file = file.path(unique(lvl3_alldf$OUTPUTPATH),"lvl3_alldf.rdata"))
      #lvl3_alldf <- lvl3_alldf[!grepl("_evt",lvl3_alldf$NAME),]
      # xaj<-ls()
      # save(xaj,file = "~/debug_lvl3.rdata")
      if(argu$run_on_pbs){
        #PBS
        #stop()
        message("Running LEVEL 3 analysis.")
        lvl3_workingdir<-file.path(argu$lvl1path_output,argu$model_name,"lvl3_misc",paste0(gsub(":","",gsub("-","_",gsub(pattern = " ","_",Sys.time()))),"log"))
        pbs_args <- get_pbs_default(); pbs_args$ppn<-4; pbs_args$walltime="40:00:00";
        qsub_commands(cmds = paste("feat",unique(lvl3_alldf$FSF_PATH)),jobperqsub = argu$job_per_qsub,
                      workingdir = lvl3_workingdir,tagname = "lvl3",qsublimit = argu$qsub_limits)
      } else {
        lvl3_cluster<-parallel::makeCluster(argu$nprocess,outfile="",type = "FORK")
        NU<-parallel::parSapply(lvl3_cluster,unique(lvl3_alldf$FSF_PATH), function(y) {
          fsl_2_sys_env()
          message("starting lvl3 feat: /n ",y)
          tryCatch(
            {system(command = paste("feat",y,sep = " "),intern = T)
              message("DONE")
            }, error=function(e){stop(paste0("feat unsuccessful...error: ", e))}
          )
          
        })
        parallel::stopCluster(lvl3_cluster)
      }
      
      
      
      
      
      
      
    } else {stop("Higher level type ",argu$lvl3_type," is not supported, only 'randomize' or 'flame' is currently supported")}
    
    message("Step ", stepnow ," Ended")
  } 
  
  #############Step 6: AFNIfy and Simple Extraction of Informaiton ###################
  stepnow<-6
  fsl_2_sys_env()
  if (is.null(argu$run_steps) | stepnow %in% argu$run_steps) {
    
    ###AFNIFY HERE
    ss_rootdir <- file.path(argu$lvl1path_output,argu$model_name,"misc_info")
    
    if(argu$lvl1_afnify){
      if(!exists("lvl2_featlist")) {
        lvl2_featlist<-system(paste0("find ",file.path(argu$lvl1path_output,argu$model_name)," -iname ","*output.feat"," -maxdepth 2 -mindepth 1 -type d"),intern = T)
      }
      
      ss_dirs<-data.frame(ss_path=lvl2_featlist,stringsAsFactors = F)
      ss_dirs$ID<-basename(dirname(ss_dirs$ss_path))
      ss_dirs$run<-gsub("_output.feat","",basename(ss_dirs$ss_path))
      
      if(nrow(ss_dirs)<1) {message("No single subject folder found. Skip")} else {
        dir.create(file.path(ss_rootdir,"ss_afni_view"),recursive = T,showWarnings = F)
        file.copy(from = file.path(Sys.getenv("FSLDIR"),"data/standard/MNI152_T1_1mm_brain.nii.gz"),
                  to = file.path(ss_rootdir,"ss_afni_view","template_brain.nii.gz"),overwrite = T)
        
        for (ssx in 1:nrow(ss_dirs)) {
          message("Start to process ID:",ss_dirs$ID[ssx],", Run: ",ss_dirs$run[ssx])
          feat2afni_single(feat_dir = ss_dirs$ss_path[ssx],include_copestat = T,include_varcope = T,include_auxstats = F,AFNIdir = argu$AFNIPATH,
                           outputdir = file.path(ss_rootdir,"ss_afni_view"),prefix = paste(ss_dirs$ID[ssx],ss_dirs$run[ssx],sep = "_"))
          message("\n")
        }
      }
    }
    
    if(argu$lvl2_afnify){
      avg_dirs<-system(paste0("find ",file.path(argu$lvl1path_output,argu$model_name)," -iname ","average.gfeat"," -maxdepth 2 -mindepth 1 -type d"),intern = T)
      if(length(avg_dirs)<1) {message("No subject average gfeat folder found. Skip")} else {
        dir.create(file.path(ss_rootdir,"avg_afni_view"),recursive = T,showWarnings = F)
        file.copy(from = file.path(Sys.getenv("FSLDIR"),"data/standard/MNI152_T1_1mm_brain.nii.gz"),to = file.path(ss_rootdir,"avg_afni_view","template_brain.nii.gz"),overwrite = T)
        for (avgx in avg_dirs) {
          message("Start to process: ",avgx)
          gfeat2afni(gfeat_dir = avgx,include_varcope = F,copy_subj_cope = F,outputdir = file.path(ss_rootdir,"avg_afni_view"),
                     prefix = paste0("avg_",basename(dirname(avgx))),AFNIdir = argu$AFNIPATH,verbos = F)
          message("\n")
        }
      }
    }
    
    if(argu$lvl3_afnify){
      fsf_ls<-list.files(path = file.path(file.path(argu$lvl3path_output,argu$model_name),"fsf_files"),pattern = ".*.fsf",full.names = T,recursive = F)
      if(length(fsf_ls)>1){
        gxroot<-unique(file.path(dirname(dirname(fsf_ls)),"grp_afni_view"))
        dir.create(gxroot,recursive = T,showWarnings = F)
        ds_path<-file.path(gxroot,"design_files")
        dir.create(ds_path,showWarnings = F,recursive = T)
        file.copy(from = file.path(Sys.getenv("FSLDIR"),"data/standard/MNI152_T1_1mm_brain.nii.gz"),to = file.path(unique(file.path(dirname(dirname(fsf_ls)),"grp_afni_view")),"template_brain.nii.gz"),overwrite = T)
        file.copy(from = file.path(unique(lvl3_alldf$OUTPUTPATH),"lvl3_alldf.rdata"),
                  to = file.path(ds_path,"lvl3_design.rdata"),overwrite = T)
        #system.file("extdata", "my_raw_data.csv", package="my_package")
        for (fsf in fsf_ls){
          
          fsf_name<-gsub(".fsf","",basename(fsf))
          message("Start to process: ",fsf_name)
          retuxa<-gfeat2afni(gfeat_dir = file.path(dirname(dirname(fsf)),gsub(".fsf",".gfeat",basename(fsf)))
                             ,include_varcope = F,copy_subj_cope = T,outputdir = gxroot,AFNIdir = argu$AFNIPATH,
                             prefix = gsub(".fsf","_grpstat",basename(fsf)),verbos = F)
          if(!is.null(retuxa)){
            file.copy(from = file.path(dirname(dirname(fsf)),gsub(".fsf",".gfeat",basename(fsf)),"design.png"),
                      to = file.path(ds_path,paste0(fsf_name,"design.png")),overwrite = T)
            dscon_raw<-readLines(file.path(dirname(dirname(fsf)),gsub(".fsf",".gfeat",basename(fsf)),"design.mat"))
            dscon_matrix<-get_matrix(raw_text = dscon_raw,heading = "/Matrix",ending = NULL,split = "\t",
                                     colnames = sub("/ContrastName\\d+\\s+([\\w_.]+).*", "\\1", grep("/ContrastName", dscon_raw, value=TRUE), perl=TRUE))
            dscon_matrix<-as.data.frame(apply(dscon_matrix,2,as.numeric))
            dscon_mk<-rbind(dscon_matrix,apply(dscon_matrix,2,mean),apply(dscon_matrix,2,sum))
            rownames(dscon_mk)[(nrow(dscon_matrix)+1):nrow(dscon_mk)]<-c("mean","sum")
            write.csv(x = dscon_mk,file = file.path(ds_path,paste0(fsf_name,"design.csv")))
          }
          message("\n")
        }
      }
    }
    
    
    
    message("Step ", stepnow ," Ended")
  }
  
  
  
  
  # stepnow<-6
  # if (is.null(argu$run_steps) | stepnow %in% argu$run_steps) {
  #   library(oro.nifti)
  #   ssfsltemp<-readLines(argu$ssub_fsl_templatepath)
  #   
  #   plot_image_all(rootpath=argu$lvl3path_output,
  #                  templatedir=argu$templatedir,
  #                  model.name=argu$model_name,
  #                  patt="*_tfce_corrp_tstat1.nii.gz",
  #                  threshold=argu$graphic.threshold,
  #                  colour="red")
  #   
  #   #Create cope index; regardless of the paths and stuff, it should be fine...
  #   if(argu$adaptive_ssfeat){
  #     xout<-rbind(
  #       data.frame(copenum=seq(argu$dsgrid$name),copename=(argu$dsgrid$name)),
  #       data.frame(copenum=seq(from=length(argu$dsgrid$name)+1,along.with = which(argu$dsgrid$AddNeg)),
  #                  copename=paste0(argu$dsgrid$name[which(argu$dsgrid$AddNeg)],"_neg"))
  #     )
  #     write.table(xout,file = file.path(argu$lvl3path_output,argu$model_name,"cope_title_index.txt"),row.names = F)
  #   }else{
  #     write.table(data.frame(copenum=paste0("cope ",as.numeric(gsub(".*?([0-9]+).*", "\\1", ssfsltemp[grep("# Title for contrast_orig",ssfsltemp)]))),
  #                            title=gsub("\"","",gsub(pattern = "[0-9]*) \"",replacement = "",
  #                                                    x = gsub(pattern = "set fmri(conname_orig.",replacement = "",
  #                                                             x = gsub(pattern = "set fmri(conname_orig.",replacement = "",
  #                                                                      x = ssfsltemp[grep("# Title for contrast_orig",ssfsltemp)+1],fixed = T),fixed = T),fixed = F))
  #     ),file = file.path(argu$lvl3path_output,argu$model_name,"cope_title_index.txt"),row.names = F)
  #   }
  #   #End of Step 6
  # }
  
  
  
  
  #############End of function fsl_pipe#####################
}
