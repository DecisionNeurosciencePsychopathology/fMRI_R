######FSL PIPE:
#New features:

#You can now supply xmat to argu for design matrix 


fsl_pipe<-function(argu=NULL, #This is the arguments environment, each model should have a different one;
                   prep.call.func="", #This should be a character string that's the name of the prep proc function
                   prep.call.allsub=list(ID=list(data=NULL)),#List of ID list of arguments for prep.call.
                   initonly=F
) {
  
  ###STEP 0: Source necessary scripts:
  #devtools::source_url("https://raw.githubusercontent.com/DecisionNeurosciencePsychopathology/fMRI_R/master/dnpl_utility.R")
  
  require("devtools")
  if("dependlab" %in% installed.packages()){"GREAT, DEPENDLAB PACK IS INSTALLED"}else{devtools::install_github("PennStateDEPENdLab/dependlab")}
  fsl_2_sys_env(force = T)
  require("parallel")
  
  if(!exists("run_pipeline",envir = argu)){
    if(Sys.getenv("run_pipeline")==""){argu$run_pipeline<-TRUE}else{argu$run_pipeline<-Sys.getenv("run_pipeline")}
  }
  if(!exists("run_steps",envir = argu)){
    if(Sys.getenv("run_steps")==""){argu$run_steps<-NULL}else{argu$run_steps<-Sys.getenv("run_steps")}
  }
  
  if (is.null(argu$nprocess)){
    if (detectCores()>12){
      num_cores<-12 
    } else {num_cores<-detectCores()-2} 
  } else {argu$nprocess->num_cores}
  
  ###Initializing argu;
  argu$cfg<-cfg_info(cfgpath = argu$cfgpath)
  argu$dsgrid<-read.table(argu$gridpath,header = T,sep = c(","),stringsAsFactors = F,strip.white = T,skipNul = T)
  if(is.null(argu$dsgrid$AddNeg)){argu$dsgrid$AddNeg<-FALSE}
  argu$dsgrid$AddNeg<-as.logical(argu$dsgrid$AddNeg)
  # if(is.null(argu$model.varinames)) {argu$model.varinames<-argu$dsgrid$name}
  
  ###Version upgrade safe keeping
  
  #Replacing old variables names from previous versions
  replaceLS<-list(ifnuisa="convlv_nuisa",onlyrun="run_steps",centerscaleall="lvl1_centervalues",lvl1_run_on_pbs="run_on_pbs",
                  model.name="model_name",ssub_outputroot="subj_outputroot",templatedir="templatebrain_path")
  
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
  
  default_ls<-list(lvl2_prep_refit=FALSE,lvl1_centervalues=FALSE,run_on_pbs=FALSE,lvl1_centervalues=TRUE,lvl1_forcegenreg=FALSE,qsub_limits=20,
                   nuisa_motion=c("nuisance","motion_par"),motion_type="fd",motion_threshold="default",job_per_qsub=as.numeric(argu$cfg$n_expected_funcruns),
                   lvl3_type="flame",adaptive_ssfeat=TRUE)
  default_ls<-default_ls[!names(default_ls) %in% names(argu)]
  if (length(default_ls)>0){
    for(lx in 1:length(default_ls)){
      message("Variable: '",names(default_ls)[lx],"' is not set, will use default value: ",default_ls[[lx]])
    }
    argu<-list2env(default_ls,envir = argu)
  }
  
  #Renaming;
  if(argu$adaptive_ssfeat){argu$ssub_fsl_templatepath<-system.file("extdata", "fsl_ssfeat_general_adaptive_template_R.fsf", package="fslpipe")}
  if(!argu$run_pipeline){return(NULL)}
  if(initonly) {return(argu)}
  
  
  
  
  #############STEP 1: Regressor generation####################
  #GENERATE REGRESSOR USING DEPENDLAB PIPELINE:
  stepnow<-1
  if (is.null(argu$run_steps) | stepnow %in% argu$run_steps) {
    #Create the directory if not existed
    #print(argu$subj_outputroot, argu$model_name)
    argu$lvl1_datalist<-prep.call.allsub;argu$lvl1_procfunc<-get(prep.call.func)
    dir.create(file.path(argu$subj_outputroot,argu$model_name),showWarnings = FALSE,recursive = T)
    #load the design rdata file if exists;
    step1_cmd<-substitute({
      allsub_design<-do_all_first_level(lvl1_datalist=argu$lvl1_datalist,lvl1_proc_func = argu$lvl1_procfunc,forcererun = argu$lvl1_forcegenreg,
                                        dsgrid=argu$dsgrid,func_nii_name=argu$func.nii.name,nprocess=argu$nprocess,
                                        cfg=argu$cfg,proc_id_subs=argu$proc_id_subs,model_name=argu$model_name,
                                        reg_rootpath=argu$regpath,center_values=argu$lvl1_centervalues,nuisance_types=argu$nuisa_motion) 
      save(allsub_design,file = file.path(argu$subj_outputroot,argu$model_name,"design.rdata"))
    })
    
    if(argu$run_on_pbs){
      workingdir<-file.path(argu$subj_outputroot,argu$model_name,"lvl1_misc")
      dir.create(workingdir,showWarnings = F,recursive = F)
      setwd(workingdir)
      save(list = ls(),file = "curwd.rdata")
      writeLines("library(fslpipe);load(\"curwd.rdata\");eval(step1_cmd)",con = "temp.r")
      pbs_torun<-get_pbs_default();pbs_torun$cmd<-"Rscript temp.r";pbs_torun$ppn<-argu$nprocess
      writeLines(do.call(pbs_cmd,pbs_torun),"pbs_temp.sh")
      dependlab::wait_for_job(dependlab::qsub_file("pbs_temp.sh"))
    } else {
      eval(step1_cmd)
    }
    #End of Step 1
  }
  
  
  #############Step 2: LVL1: Single Subject PARALLEL##########
  #Now we do the single sub processing using FSL and the regressor that was generated
  stepnow<-2
  if (is.null(argu$run_steps) | stepnow %in% argu$run_steps) {
    
    
    if (file.exists(file.path(argu$subj_outputroot,argu$model_name,"design.rdata"))) {
      load(file.path(argu$subj_outputroot,argu$model_name,"design.rdata"))
    } else {stop("No design rdata file found, must re-run step 1")}
    
    
    #let's subset this 
    small.sub<-lapply(as.list(allsub_design), function(x) {
      list(
        design=x$design,
        ID=x$ID,
        run_volumes=x$design$run_volumes,
        regpath=x$regpath,
        preprocID=x$preprocID)
    })
    
    #IT's the same for all participants!!!! WHY DO YOU RE RUN IT FOR EVERYONE!!!
    xarg<-as.environment(list())
    xarg$templatebrain<-argu$templatedir
    xarg$tr<-argu$cfg$preproc_call$tr
    
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
        aarg$outputpath<-file.path(argu$subj_outputroot,argu$model_name,idx,paste0("run",runnum,"_output"))
        if (!file.exists(paste0(aarg$outputpath,".feat")) ) {
          if(is.null(argu$ss_zthreshold)) {aarg$zthreshold<-3.2} else {aarg$zthreshold<-argu$ss_zthreshold}
          if(is.null(argu$ss_pthreshold)) {aarg$pthreshold<-0.05} else {aarg$pthreshold<-argu$ss_pthreshold}
          message(paste0("Initializing feat for participant: ",idx,", and run: ",runnum))
          aarg$runnum<-runnum   
          aarg$volumes<-x$run_volumes[runnum]
          aarg$funcfile<-get_volume_run(id=paste0(idx,argu$proc_id_subs),cfg = argu$cfg,reg.nii.name = argu$func.nii.name,returnas = "path")[runnum]
          if(!length(aarg$funcfile)>0 || any(is.na(aarg$funcfile))){message("ID: ",idx," RUN: ",runnum,", failed to find a functional image. Terminate.");return(NULL)}
          aarg$nuisa<-file.path(argu$regpath,argu$model_name,idx,paste0("run",runnum,"_nuisance_regressor_with_motion.txt"))
          
          if(argu$adaptive_ssfeat){
            for (mv in 1:ncol(argu$lvl1_cmat)) {
              assign(paste0("evreg",mv),file.path(file.path(argu$regpath,argu$model_name),idx,paste0("run",runnum,"_",argu$model_varinames[mv],argu$regtype)),envir = aarg)
            }
          } else {
            gen_reg(vmodel=argu$model_varinames,regpath=file.path(argu$regpath,argu$model_name),idx=idx,runnum=runnum,env=aarg,regtype=argu$regtype)
          }
          if (any(unlist(eapply(aarg,is.na)))) {stop("NA exists in one of the arguments; please double check!")}
          #gen_reg(vmodel=argu$model.varinames,regpath=file.path(argu$regpath,argu$model_name),idx=idx,runnum=runnum,env=xarg,regtype = argu$regtype)
          
          cmmd<-feat_w_template(fsltemplate = ssfsltemp,beg = "ARG_",end = "_END",
                                fsfpath = file.path(argu$regpath,argu$model_name,idx,paste0("run",runnum,"_",argu$model_name,".fsf")),
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
        lvl1_workingdir<-file.path(argu$subj_outputroot,argu$model_name,"lvl1_misc",paste0(gsub(":","",gsub("-","_",gsub(pattern = " ","_",Sys.time()))),"log"))
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
    #End of Step 2
  }
  
  #############Step 3: Prep for Higher Level #######################
  #Now we make the symbolic link for template matching...so they are not misaligned anymore...
  stepnow<-3
  if (is.null(argu$run_steps) | stepnow %in% argu$run_steps) {
    if(argu$lvl2_prep_refit){
      #cfg<-cfg_info(cfgpath = argu$cfgpath)
      gtax<-prepare4secondlvl(
        ssana.path=file.path(argu$subj_outputroot,argu$model_name),            
        standardbarin.path=argu$templatedir, featfoldername = "*output.feat",
        dir.filter=argu$hig_lvl_path_filter,                                                
        overwrite=argu$lvl2_force_prep,
        outputmap=TRUE,
        paralleln = num_cores)           
    } else {
      lvl2_featlist<-prep_session_lvl(subj_rootpath = file.path(argu$subj_outputroot,argu$model_name),subj_folderreg = "*output.feat",overwrite = argu$lvl2_force_prep,
                                      template_brainpath = argu$templatebrain_path)
    }
    ##End of Step 3
  }
  
  #############Step 4: LVL2: Fixed Effect for Single Subject PARALLEL ###############
  #This starts averaging for each subject:
  stepnow<-4
  if (is.null(argu$run_steps) | stepnow %in% argu$run_steps) {
    
    lowlvl_featreg = "*output.feat"
    if(!exists("lvl2_featlist")) {
      lvl2_featlist<-system(paste0("find ",file.path(argu$subj_outputroot,argu$model_name)," -iname ",lowlvl_featreg," -maxdepth 4 -mindepth 1 -type d"),intern = T)
    }
    raw.split <- strsplit(lvl2_featlist,split = .Platform$file.sep) 
    lvl2_rawdf<-do.call(rbind,lapply(raw.split,function(x){
      data.frame(
        uID=x[grep(lowlvl_featreg,x)-1],
        ID=paste(x[grep(lowlvl_featreg,x)-1],gsub("run","",gsub("_output.feat","",x[grep(lowlvl_featreg,x)],fixed = T)),sep = "_"),
        PATH = paste(x,collapse = .Platform$file.sep),NAME="average",
        OUTPUTPATH = dirname(paste(x,collapse = .Platform$file.sep)),
        SUBJMEAN=1,
        stringsAsFactors = F)
    })
    )
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
    
    if(nrow(lvl2_alldf)>0){
      if(argu$run_on_pbs){
        #PBS
        lvl2_workingdir<-file.path(argu$subj_outputroot,argu$model_name,"lvl2_misc",paste0(gsub(":","",gsub("-","_",gsub(pattern = " ","_",Sys.time()))),"log"))
        qsub_commands(cmds = paste("feat",unique(lvl2_alldf$FSF_PATH)),jobperqsub = argu$job_per_qsub,
                      workingdir = lvl2_workingdir,tagname = "lvl2",ppn = 4,qsublimit = argu$qsub_limits)
        
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

    ################END of step 4
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
      if(argu$adaptive_ssfeat) {maxcopenum<-1:nrow(argu$xmat)} else {
        maxcopenum<-1:max(as.numeric(gsub(".*?([0-9]+).*", "\\1", ssfsltemp[grep("# Title for contrast",ssfsltemp)])))
      }
      #Start Group Level Analysis:
      glvl_all_cope(rootdir=argu$subj_outputroot,
                    outputdir=argu$glvl_outputroot,
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
      # rootdir=argu$subj_outputroot
      # outputdir=argu$glvl_outputroot
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
      
      lowlvl_featreg<-"average.gfeat"
      raw<-system(paste0("find ",
                         file.path(argu$subj_outputroot,argu$model_name,"*/",lowlvl_featreg),
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
      
      lvl3_raw_sp<-split(lvl3_rawdf,lvl3_rawdf$COPENUM)
      
      lvl3_raw_sp<-lapply(lvl3_raw_sp,function(x){
        if(!is.null(argu$lvl3_ref_df)){x<-merge(x,argu$lvl3_ref_df,all.x=T,by="ID")}
        x$NAME = unique(readLines(file.path(x$PATH[1],"design.lev")))
        x$OUTPUTPATH = file.path(argu$glvl_output,argu$model_name)
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
      lvl3_alldf <- lvl3_alldf[!grepl("_evt",lvl3_alldf$NAME),]
      # xaj<-ls()
      # save(xaj,file = "~/debug_lvl3.rdata")
      if(argu$run_on_pbs){
        #PBS
        #stop()
        message("Running LEVEL 3 analysis.")
        lvl3_workingdir<-file.path(argu$subj_outputroot,argu$model_name,"lvl3_misc",paste0(gsub(":","",gsub("-","_",gsub(pattern = " ","_",Sys.time()))),"log"))
        qsub_commands(cmds = paste("feat",unique(lvl3_alldf$FSF_PATH)),jobperqsub = argu$job_per_qsub,
                      workingdir = lvl3_workingdir,tagname = "lvl3",ppn = 4,qsublimit = argu$qsub_limits)
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
    
    #End Step 5
  } 
  
  #############Step 6: Simple Graph and Info Extraction ###################
  # Will not be compatible with new pipeline at all.
  
  # stepnow<-6
  # if (is.null(argu$run_steps) | stepnow %in% argu$run_steps) {
  #   library(oro.nifti)
  #   ssfsltemp<-readLines(argu$ssub_fsl_templatepath)
  #   
  #   plot_image_all(rootpath=argu$glvl_outputroot,
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
  #     write.table(xout,file = file.path(argu$glvl_outputroot,argu$model_name,"cope_title_index.txt"),row.names = F)
  #   }else{
  #     write.table(data.frame(copenum=paste0("cope ",as.numeric(gsub(".*?([0-9]+).*", "\\1", ssfsltemp[grep("# Title for contrast_orig",ssfsltemp)]))),
  #                            title=gsub("\"","",gsub(pattern = "[0-9]*) \"",replacement = "",
  #                                                    x = gsub(pattern = "set fmri(conname_orig.",replacement = "",
  #                                                             x = gsub(pattern = "set fmri(conname_orig.",replacement = "",
  #                                                                      x = ssfsltemp[grep("# Title for contrast_orig",ssfsltemp)+1],fixed = T),fixed = T),fixed = F))
  #     ),file = file.path(argu$glvl_outputroot,argu$model_name,"cope_title_index.txt"),row.names = F)
  #   }
  #   #End of Step 6
  # }
  
  
  #############End of function fsl_pipe#####################
}


#In development:
if (FALSE) {
  #Flame
  
}
