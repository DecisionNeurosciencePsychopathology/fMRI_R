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
  
  
  #IDTORUN <- names(prep.call.allsub)
  
  ###Initializing argu;
  argu$cfg<-cfg_info(cfgpath = argu$cfgpath)
  argu$dsgrid<-read.table(argu$gridpath,header = T,sep = c(","),stringsAsFactors = F,strip.white = T,skipNul = T)
  if(is.null(argu$dsgrid$AddNeg)){argu$dsgrid$AddNeg<-FALSE}
  if(is.null(argu$dsgrid$RunGrpLvl)){argu$dsgrid$RunGrpLvl<-TRUE}
  argu$dsgrid$AddNeg<-as.logical(argu$dsgrid$AddNeg)
  argu$dsgrid$RunGrpLvl <- as.logical(argu$dsgrid$RunGrpLvl)
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
  
  default_ls<-list(lvl2_prep_refit=FALSE,lvl1_centervalues=FALSE,run_on_pbs=FALSE,lvl1_centervalues=TRUE,lvl1_forcegenreg=FALSE,
                   qsub_limits=20,lvl2_force_prep=FALSE,lvl1_retry=FALSE,lvl1_afnify=F,lvl2_afnify=F,lvl3_afnify=T,
                   nuisa_motion=c("nuisance","motion_par"),lvl3_lowlvlfeatreg="average.gfeat",motion_type="fd",
                   motion_threshold="default",job_per_qsub=as.numeric(argu$cfg$n_expected_funcruns),
                   lvl3_type="flame",adaptive_ssfeat=TRUE)
  default_ls<-default_ls[!names(default_ls) %in% names(argu)]
  if (length(default_ls)>0){
    for(lx in 1:length(default_ls)){
      message("Variable: '",names(default_ls)[lx],"' is not set, will use default value: ",default_ls[[lx]])
    }
    argu<-list2env(default_ls,envir = argu)
  }
  argu$lvl3_lowlvlfeatreg
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
    
    #argu$lvl1_datalist<-argu$lvl1_datalist[which(names(argu$lvl1_datalist) %in% IDTORUN)]
    dir.create(file.path(argu$subj_outputroot,argu$model_name),showWarnings = FALSE,recursive = T)
    #load the design rdata file if exists;
    step1_cmd<-substitute({
      allsub_design<-fslpipe::do_all_first_level(lvl1_datalist=argu$lvl1_datalist,lvl1_proc_func = argu$lvl1_procfunc,forcererun = argu$lvl1_forcegenreg,
                                                 dsgrid=argu$dsgrid,func_nii_name=argu$func.nii.name,nprocess=argu$nprocess,
                                                 cfg=argu$cfg,proc_id_subs=argu$proc_id_subs,model_name=argu$model_name,retry=argu$lvl1_retry,
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
    
    #small.sub <- small.sub[which(sapply(small.sub,`[[`,'ID') %in% IDTORUN)]
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
        if (!file.exists(file.path(paste0(aarg$outputpath,".feat"),"stats","zstat1.nii.gz")) ) {
          message(paste0("Initializing feat for participant: ",idx,", and run: ",runnum))
          if(dir.exists(file.path(paste0(aarg$outputpath,".feat"))) ){
            message("Found an incomplete folder...Removing...")
            unlink(file.path(paste0(aarg$outputpath,".feat")),recursive = T,force = T)}
          if(is.null(argu$ss_zthreshold)) {aarg$zthreshold<-3.2} else {aarg$zthreshold<-argu$ss_zthreshold}
          if(is.null(argu$ss_pthreshold)) {aarg$pthreshold<-0.05} else {aarg$pthreshold<-argu$ss_pthreshold}
          
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
      lvl2_featlist<-system(paste0("find ",file.path(argu$subj_outputroot,argu$model_name)," -iname ",lowlvl_featreg," -maxdepth 2 -mindepth 1 -type d"),intern = T)
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
        lvl2_workingdir<-file.path(argu$subj_outputroot,argu$model_name,"lvl2_misc",paste0(gsub(":","",gsub("-","_",gsub(pattern = " ","_",Sys.time()))),"log"))
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
                    outputdir=argu$glvl_output,
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
      # outputdir=argu$glvl_output
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
      #lvl3_raw_sp<-lvl3_raw_sp[which(argu$dsgrid$RunGrpLvl)]
      
      lvl3_raw_sp<-lapply(lvl3_raw_sp,function(x){
        if(!is.null(argu$lvl3_ref_df)){x<-merge(x,argu$lvl3_ref_df,all.x=T,by="ID")}
        x$NAME = unique(readLines(file.path(x$PATH[1],"design.lev")))
        x$OUTPUTPATH = file.path(argu$glvl_output,paste(argu$model_name,paste(argu$lvl3_covariate_names,collapse = "_"),sep = "_"))
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
      
      if(!is.null(argu$run_these_ID)) {
        lvl3_alldf <- lvl3_alldf[which(lvl3_alldf$ID %in% argu$run_these_ID),]
      }
      
      if(!is.null(argu$exclude_these_ID)) {
        message("The following IDs are excluded from the final analysis: ",paste(argu$exclude_these_ID,collapse = ", "))
        lvl3_alldf <- lvl3_alldf[which(!lvl3_alldf$ID %in% argu$exclude_these_ID),]
        
      }
      save(lvl3_alldf,file = file.path(unique(lvl3_alldf$OUTPUTPATH),"lvl3_alldf.rdata"))
      #lvl3_alldf <- lvl3_alldf[!grepl("_evt",lvl3_alldf$NAME),]
      # xaj<-ls()
      # save(xaj,file = "~/debug_lvl3.rdata")
      if(argu$run_on_pbs){
        #PBS
        #stop()
        message("Running LEVEL 3 analysis.")
        lvl3_workingdir<-file.path(argu$subj_outputroot,argu$model_name,"lvl3_misc",paste0(gsub(":","",gsub("-","_",gsub(pattern = " ","_",Sys.time()))),"log"))
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
    
    #End Step 5
  } 
  
  #############Step 6: AFNIfy and Simple Extraction of Informaiton ###################
  stepnow<-6
  if (is.null(argu$run_steps) | stepnow %in% argu$run_steps) {
    fsl_2_sys_env()
    if(argu$lvl1_afnify){
      if(!exists("lvl2_featlist")) {
        lvl2_featlist<-system(paste0("find ",file.path(argu$subj_outputroot,argu$model_name)," -iname ","*output.feat"," -maxdepth 2 -mindepth 1 -type d"),intern = T)
      }
      
      ss_dirs<-data.frame(ss_path=lvl2_featlist,stringsAsFactors = F)
      ss_dirs$ID<-basename(dirname(ss_dirs$ss_path))
      ss_dirs$run<-gsub("_output.feat","",basename(ss_dirs$ss_path))
      if(nrow(ss_dirs)<1) {message("No single subject folder found. Skip")} else {
        dir.create(file.path(ss_rootdir,"ss_afni_view"),recursive = T,showWarnings = F)
        file.copy(from = file.path(Sys.getenv("FSLDIR"),"data/standard/MNI152_T1_1mm_brain.nii.gz"),to = file.path(ss_rootdir,"ss_afni_view","template_brain.nii.gz"),overwrite = T)
        
        for (ssx in 1:nrow(ss_dirs)) {
          message("Start to process ID:",ss_dirs$ID[ssx],", Run: ",ss_dirs$run[ssx])
          feat2afni_single(feat_dir = ss_dirs$ss_path[ssx],include_copestat = T,include_varcope = T,include_auxstats = F,AFNIdir = argu$AFNIPATH,
                           outputdir = file.path(ss_rootdir,"ss_afni_view"),prefix = paste(ss_dirs$ID[ssx],ss_dirs$run[ssx],sep = "_"))
          message("\n")
        }
      }
    }
    
    if(argu$lvl2_afnify){
      avg_dirs<-system(paste0("find ",file.path(argu$subj_outputroot,argu$model_name)," -iname ","average.gfeat"," -maxdepth 2 -mindepth 1 -type d"),intern = T)
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
      fsf_ls<-list.files(path = file.path(file.path(argu$glvl_output,argu$model_name),"fsf_files"),pattern = ".*.fsf",full.names = T,recursive = F)
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
    
  }
  
  get_matrix<-function(raw_text,heading="/Matrix",ending=NULL,colnames=NULL,split=" "){
    if(is.null(ending)){end_pos<-length(raw_text)}else{end_pos<-(grep(ending,raw_text)-1)}
    if(is.null(heading)){end_pos<-0}else{start_pos<-(grep(heading,raw_text)+1)}
    raw_mx<-raw_text[start_pos:end_pos]
    matrix_df<-as.data.frame(do.call(rbind,strsplit(raw_mx,split = split)))
    if(!is.null(colnames) && length(colnames)==ncol(matrix_df)){names(matrix_df)<-colnames} 
    return(matrix_df)
  }
  
  
  # stepnow<-6
  # if (is.null(argu$run_steps) | stepnow %in% argu$run_steps) {
  #   library(oro.nifti)
  #   ssfsltemp<-readLines(argu$ssub_fsl_templatepath)
  #   
  #   plot_image_all(rootpath=argu$glvl_output,
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
  #     write.table(xout,file = file.path(argu$glvl_output,argu$model_name,"cope_title_index.txt"),row.names = F)
  #   }else{
  #     write.table(data.frame(copenum=paste0("cope ",as.numeric(gsub(".*?([0-9]+).*", "\\1", ssfsltemp[grep("# Title for contrast_orig",ssfsltemp)]))),
  #                            title=gsub("\"","",gsub(pattern = "[0-9]*) \"",replacement = "",
  #                                                    x = gsub(pattern = "set fmri(conname_orig.",replacement = "",
  #                                                             x = gsub(pattern = "set fmri(conname_orig.",replacement = "",
  #                                                                      x = ssfsltemp[grep("# Title for contrast_orig",ssfsltemp)+1],fixed = T),fixed = T),fixed = F))
  #     ),file = file.path(argu$glvl_output,argu$model_name,"cope_title_index.txt"),row.names = F)
  #   }
  #   #End of Step 6
  # }
  
  
  #############End of function fsl_pipe#####################
}
