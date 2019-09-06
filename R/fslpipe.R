######FSL PIPE:
#New features:

#You can now supply xmat to argu for design matrix 


fsl_pipe<-function(argu=NULL, #This is the arguments environment, each model should have a different one;
                   prep.call.func="", #This should be a character string that's the name of the prep proc function
                   prep.call.allsub=list(ID=list(data=NULL)) #List of ID list of arguments for prep.call.
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
  
  if (is.null(argu$nprocess)){
    if (detectCores()>12){
      num_cores<-12 #Use 8 cores to minimize burden; if on throndike 
      #Or if you are running this on laptop; whatever cores minus 2; I guess if it's a dual core...let's just don't do that (zero core will not paralle anything)
    } else {num_cores<-detectCores()-2} 
  } else {argu$nprocess->num_cores}
  
  ###Initializing argu;
  argu$cfg<-cfg_info(cfgpath = argu$cfgpath)
  argu$dsgrid<-read.table(argu$gridpath,header = T,sep = c(","),stringsAsFactors = F,strip.white = T,skipNul = T)
  if(is.null(argu$dsgrid$AddNeg)){argu$dsgrid$AddNeg<-FALSE}
  argu$dsgrid$AddNeg<-as.logical(argu$dsgrid$AddNeg)
  if(is.null(argu$model.varinames)) {argu$model.varinames<-argu$dsgrid$name}
  
  ###Version upgrade safe keeping
  default_ls<-list(lvl2_prep_refit=FALSE,centerscaleall=FALSE,
                     nuisa_motion=c("nuisance","motion_par"),motion_type="fd",motion_threshold="default",
                     lvl3_type="flame",adaptive_ssfeat=TRUE)
  default_ls<-default_ls[!names(default_ls) %in% names(argu)]
  for(lx in 1:length(default_ls)){
    message("Variable: '",names(default_ls)[lx],"' is not set, will use default value: ",default_ls[[lx]])
  }
  argu<-list2env(default_ls,envir = argu)
  
  #RE-config
  if (exists("ifnuisa",envir = argu) & !exists("convlv_nuisa",envir = argu)) {
    message("ifnuisa variable is now depreciated, please use convlv_nuisa to control if the pipeline should convolve nuissance regressor")
    argu$convlv_nuisa<-argu$ifnuisa}
  
  if (exists("onlyrun",envir = argu) & !exists("run_steps",envir = argu)) {
    message("onlyrun variable is now depreciated, please use run_steps")
    argu$run_steps<-argu$onlyrun}
  
  if (!exists("xmat",envir = argu)) {
    message("Single subject design matrix is not specified, will use automatic generated one by using grid.")
    ogLength<-length(argu$dsgrid$name)
    negNum<-(which(argu$dsgrid$AddNeg))
    if(length(negNum)>0){
      negMat<-diag(x=-1,ogLength)
      negMat<-negMat[negNum,]
      argu$xmat<-rbind(diag(x=1,ogLength),negMat)
      
    } else {
      argu$xmat<-diag(x=1,ogLength)
    }
  } 
  
  #Renaming;
  argu$model_name <-   argu$model.name
  argu$subj_outputroot <-  argu$ssub_outputroot
  argu$templatebrain_path <- argu$templatedir

  if(!argu$run_pipeline){return(NULL)}
  #############STEP 1: Regressor generation####################
  #GENERATE REGRESSOR USING DEPENDLAB PIPELINE:
  stepnow<-1
  if (is.null(argu$run_steps) | stepnow %in% argu$run_steps) {
    #Create the directory if not existed
    dir.create(file.path(argu$ssub_outputroot,argu$model.name),showWarnings = FALSE,recursive = T)
    #load the design rdata file if exists;
    if (file.exists(file.path(argu$ssub_outputroot,argu$model.name,"design.rdata"))) {
      tryCatch({
        load(file.path(argu$ssub_outputroot,argu$model.name,"design.rdata"))
      }, error=function(e) {
        message(paste0("load not successful, have to re-run step 1...message: ",e))
        assign('allsub.design',as.environment(list()),envir = globalenv())
      })
      
    } else {allsub.design<-as.environment(list())}
    #Take out people who have already been processed;
    if (length(names(allsub.design))>0 & !argu$forcereg) {
      idtodo<-as.character(names(prep.call.allsub)[which(! names(prep.call.allsub) %in% names(allsub.design))])
      message(paste("These IDs already has regressors: ",paste(names(allsub.design),collapse=" "),sep=" ",collapse = " "))
    } else {idtodo<-names(prep.call.allsub)}
    
    assign(prep.call.func,get(prep.call.func),envir = argu)
    if (length(idtodo)>0) {
      for (xid in idtodo) {
        prep.call.list<-prep.call.allsub[[xid]]
        tryCatch(
          {
            message(xid)
            assign(as.character(xid),
                   do.all.subjs(
                     tid=xid,
                     do.prep.call=prep.call.func,
                     do.prep.arg=prep.call.list,
                     cfgpath=argu$cfgpath,
                     regpath=argu$regpath,
                     gridpath=argu$gridpath,
                     func.nii.name=argu$func.nii.name,
                     proc_id_subs=argu$proc_id_subs,    #Put "" for nothing.
                     wrt.timing=c("convolved", "FSL"),
                     model.name=argu$model.name,
                     model.varinames=argu$model.varinames,
                     nuisa_motion=argu$nuisa_motion,
                     motion_type=argu$motion_type,
                     motion_threshold=argu$motion_threshold,
                     convlv_nuisa=argu$convlv_nuisa,
                     centerscaleall=argu$centerscaleall,
                     argu=argu
                   ),envir = allsub.design
            )
            # tid=xid
            # do.prep.call=prep.call.func
            # do.prep.arg=prep.call.list
            # cfgpath=argu$cfgpath
            # regpath=argu$regpath
            # gridpath=argu$gridpath
            # func.nii.name=argu$func.nii.name
            # proc_id_subs=argu$proc_id_subs   #Put "" for nothing.
            # wrt.timing=c("convolved", "FSL")
            # model.name=argu$model.name
            # func.nii.name=argu$func.nii.name
            # proc_id_subs=argu$proc_id_subs
            # model.name=argu$model.name
            # model.varinames=argu$model.varinames
            # add.nuisa=argu$ifnuisa
            # nuisa_motion=argu$nuisa_motion
            # motion_type=argu$motion_type
            # motion_threshold=argu$motion_threshold
            # convlv_nuisa=argu$convlv_nuisa
            # centerscaleall=argu$centerscaleall
            
          },error=function(e) {message("failed regressor generation...go investigate: ",e)}
        )
      }
      
      save("allsub.design",file = file.path(argu$ssub_outputroot,argu$model.name,"design.rdata"))
    } else {message("NO NEW DATA NEEDED TO BE PROCESSED")}
    
    #End of Step 1
  }
  

  #############Step 2: LVL1: Single Subject PARALLEL##########
  #Now we do the single sub processing using FSL and the regressor that was generated
  stepnow<-2
  if (is.null(argu$run_steps) | stepnow %in% argu$run_steps) {
    
    if (!is.null(argu$run_steps) & !1 %in% argu$run_steps) {
      if (file.exists(file.path(argu$ssub_outputroot,argu$model.name,"design.rdata"))) {
        load(file.path(argu$ssub_outputroot,argu$model.name,"design.rdata"))
      } else {stop("No design rdata file found, must re-run step 1")}
    }  
    
    #let's subset this 
    small.sub<-eapply(allsub.design, function(x) {
      list(
        design=x$design,
        ID=x$ID,
        run_volumes=x$run_volumes,
        regpath=x$regpath,
        preprocID=x$preprocID)
    })
    
    #IT's the same for all participants!!!! WHY DO YOU RE RUN IT FOR EVERYONE!!!
    xarg<-as.environment(list())
    xarg$templatebrain<-argu$templatedir
    xarg$tr<-argu$cfg$preproc_call$tr
    
    if(argu$adaptive_ssfeat){
      argu$model.varinames<-argu$dsgrid$name
      argu$copenames<-c(argu$model.varinames,paste0(argu$dsgrid$name[which(argu$dsgrid$AddNeg)],"_neg"))
      xarg$evnum<-1:ncol(argu$xmat)
      xarg$copenum<-1:nrow(argu$xmat)
      xarg$maxevnum<-ncol(argu$xmat)
      xarg$maxcopenum<-nrow(argu$xmat)
      argu$maxcopenum<-nrow(argu$xmat)
      
      for (xy in 1:nrow(argu$xmat)){
        assign(paste0("copemat",xy),value = 1:ncol(argu$xmat),envir = xarg)
        assign(paste0("copetitle",xy),argu$copenames[xy],envir = xarg)
        assign(paste0("cope_lessnum",xy),(1:length(argu$copenames))[-xy],envir = xarg)
        for(xx in 1:ncol(argu$xmat)){
          assign(paste0("copevalue",xy,"_",xx),value = argu$xmat[xy,xx],envir = xarg)
        }
      } 
      for (mv in 1:ncol(argu$xmat)) {
        assign(paste0("evtitle",mv),argu$model.varinames[mv],envir = xarg)
        assign(paste0("orthocombo",mv),paste(mv,(0:length(argu$model.varinames)),sep = "."),envir = xarg)
      }
      #PUT NEW FUNCTION HERE
      ssfsltemp<-rep_within(adptemplate = readLines(argu$ssub_fsl_templatepath),searchenvir = xarg)
    } else {ssfsltemp<-readLines(argu$ssub_fsl_templatepath)}
    
    step2commands<-unlist(lapply(small.sub,function(x) {
      idx<-x$ID
      aarg<-xarg
      cmmd<-unlist(lapply(1:length(x$run_volumes), function(runnum) {
        aarg$outputpath<-file.path(argu$ssub_outputroot,argu$model.name,idx,paste0("run",runnum,"_output"))
        if (!file.exists(paste0(aarg$outputpath,".feat")) ) {
          if(is.null(argu$ss_zthreshold)) {aarg$zthreshold<-3.2} else {aarg$zthreshold<-argu$ss_zthreshold}
          if(is.null(argu$ss_pthreshold)) {aarg$pthreshold<-0.05} else {aarg$pthreshold<-argu$ss_pthreshold}
          
          message(paste0("Initializing feat for participant: ",idx,", and run: ",runnum))
          aarg$runnum<-runnum   
          aarg$volumes<-x$run_volumes[runnum]
          aarg$funcfile<-get_volume_run(id=paste0(idx,argu$proc_id_subs),cfgfilepath = argu$cfgpath,reg.nii.name = argu$func.nii.name,returnas = "path")[runnum]
          aarg$nuisa<-file.path(argu$regpath,argu$model.name,idx,paste0("run",runnum,"_nuisance_regressor_with_motion.txt"))
          
          if(argu$adaptive_ssfeat){
            for (mv in 1:ncol(argu$xmat)) {
              assign(paste0("evreg",mv),file.path(file.path(argu$regpath,argu$model.name),idx,paste0("run",runnum,"_",argu$model.varinames[mv],argu$regtype)),envir = aarg)
            }
          } else {
            gen_reg(vmodel=argu$model.varinames,regpath=file.path(argu$regpath,argu$model.name),idx=idx,runnum=runnum,env=aarg,regtype=argu$regtype)
          }
          if (any(unlist(eapply(aarg,is.na)))) {stop("NA exists in one of the arguments; please double check!")}
          #gen_reg(vmodel=argu$model.varinames,regpath=file.path(argu$regpath,argu$model.name),idx=idx,runnum=runnum,env=xarg,regtype = argu$regtype)
          
          cmmd<-feat_w_template(fsltemplate = ssfsltemp,beg = "ARG_",end = "_END",
                                fsfpath = file.path(argu$regpath,argu$model.name,idx,paste0("run",runnum,"_",argu$model.name,".fsf")),
                                envir = aarg,outcommand = T)
          rm(aarg)
          return(cmmd)
        } else {
          message(paste("ID:",idx,"RUN:",runnum,",already exists,","to re-run, remove the directory."))
          return(NULL)}
      }))
      return(cmmd)
    }))
    
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
    
    #End of Step 2
  }
  
  #############Step 3: Prep for Higher Level #######################
  #Now we make the symbolic link for template matching...so they are not misaligned anymore...
  stepnow<-3
  if (is.null(argu$run_steps) | stepnow %in% argu$run_steps) {
    if(argu$lvl2_prep_refit){
      #cfg<-cfg_info(cfgpath = argu$cfgpath)
      gtax<-prepare4secondlvl(
        ssana.path=file.path(argu$ssub_outputroot,argu$model.name),            
        standardbarin.path=argu$templatedir, featfoldername = "*output.feat",
        dir.filter=argu$hig_lvl_path_filter,                                                
        overwrite=argu$ifoverwrite_secondlvl,
        outputmap=TRUE,
        paralleln = num_cores)           
    } else {
      lvl2_featlist<-prep_session_lvl(subj_rootpath = file.path(argu$subj_outputroot,argu$model_name),subj_folderreg = "*output.feat",
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
                       overwrite = F,
                       vari_to_run=c("SUBJMEAN")
    )
    lvl2_arg<-lapply(names(lvl2_default),function(xa){
      if(exists(paste0("lvl2_",xa),envir = argu)){get(paste0("lvl2_",xa),envir = argu)}else{lvl2_default[[xa]]}
    })
    names(lvl2_arg) <- names(lvl2_default)
    lvl2_arg$proc_ls_fsf <- lvl2_raw_sp 
    lvl2_arg$template_brain <- argu$templatebrain_path
    lvl2_arg$fsltemplate <- readLines(system.file("extdata", "fsl_flame_general_adaptive_template.fsf", package="fslpipe"))
    lvl2_alldf <- do.call(gen_fsf_highlvl,lvl2_arg)
    
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
      glvl_all_cope(rootdir=argu$ssub_outputroot,
                    outputdir=argu$glvl_outputroot,
                    modelname=argu$model.name,
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
      # rootdir=argu$ssub_outputroot
      # outputdir=argu$glvl_outputroot
      # modelname=argu$model.name
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
    } else if (argu$lvl3_type=="flame") {
      #Run flame here:

      lowlvl_featreg<-"average.gfeat"
      raw<-system(paste0("find ",
                         file.path(argu$subj_outputroot,argu$model_name,"*/",lowlvl_featreg),
                         " -iname '*.feat' -maxdepth 2 -mindepth 1 -type d"),intern = T)
      raw.split <- strsplit(raw,split = .Platform$file.sep)  
      lvl3_rawdf<-do.call(rbind,lapply(raw.split,function(x){
        data.frame(
          ID=x[grep(lowlvl_featreg,x)-1],
          COPENUM=gsub("cope","",gsub(".feat","",last(x))),
          PATH = paste(x,collapse = .Platform$file.sep),
          INTERCEPT=1,
          stringsAsFactors = F)
      })
      )
      # NAME = ?
      # OUTPUTPATH = ?
        
      lvl3_raw_sp<-split(lvl3_rawdf,lvl3_rawdf$COPENUM)
      
      lvl3_raw_sp<-lapply(lvl3_raw_sp,function(x){
        x$NAME = unique(readLines(file.path(x$PATH[1],"design.lev")))
        x$OUTPUTPATH = file.path(argu$glvl_output,argu$model_name)
        return(x)
      })
      
      lvl3_default<-list(flame_type = 1, #FLAME 1 for sess_lvl and FLAME 1+2 for grp_lvl
                         thresh_type = 3, #0 : None \n # 1 : Uncorrected \n# 2 : Voxel \n # 3 : Cluster \n"
                         z_thresh = 3.09, #1.96 for subj lvl, 2.3 for sess lvl and 3.09 for grp lvl
                         p_thresh = 0.05, #0.05 default for both;
                         overwrite = F,
                         vari_to_run=c("INTERCEPT")
      )
      lvl3_arg<-lapply(names(lvl3_default),function(xa){
        if(exists(paste0("lvl3_",xa),envir = argu)){get(paste0("lvl3_",xa),envir = argu)}else{lvl3_default[[xa]]}
      })
      names(lvl3_arg) <- names(lvl3_default)
      lvl3_arg$proc_ls_fsf <- lvl3_raw_sp 
      lvl3_arg$template_brain <- argu$templatebrain_path
      lvl3_arg$fsltemplate <- readLines(system.file("extdata", "fsl_flame_general_adaptive_template.fsf", package="fslpipe"))
      
      lvl3_alldf <- do.call(gen_fsf_highlvl,lvl3_arg)
      lvl3_alldf <- lvl3_alldf[!grepl("_evt",lvl3_alldf$NAME),]
      
      lvl3_cluster<-parallel::makeCluster(argu$nprocess,outfile="",type = "FORK")
      NU<-parallel::parSapply(lvl3_cluster,unique(lvl3_alldf$FSF_PATH), function(y) {
        fsl_2_sys_env()
        message("starting feat /n ",y)
        tryCatch(
          {system(command = paste("feat",y,sep = " "),intern = T)
            message("DONE")
          }, error=function(e){stop(paste0("feat unsuccessful...error: ", e))}
        )
        
      })
      parallel::stopCluster(lvl3_cluster)

      
      #Variables to get:
      #outputpath:
      #
      
      
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
  #                  model.name=argu$model.name,
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
  #     write.table(xout,file = file.path(argu$glvl_outputroot,argu$model.name,"cope_title_index.txt"),row.names = F)
  #   }else{
  #     write.table(data.frame(copenum=paste0("cope ",as.numeric(gsub(".*?([0-9]+).*", "\\1", ssfsltemp[grep("# Title for contrast_orig",ssfsltemp)]))),
  #                            title=gsub("\"","",gsub(pattern = "[0-9]*) \"",replacement = "",
  #                                                    x = gsub(pattern = "set fmri(conname_orig.",replacement = "",
  #                                                             x = gsub(pattern = "set fmri(conname_orig.",replacement = "",
  #                                                                      x = ssfsltemp[grep("# Title for contrast_orig",ssfsltemp)+1],fixed = T),fixed = T),fixed = F))
  #     ),file = file.path(argu$glvl_outputroot,argu$model.name,"cope_title_index.txt"),row.names = F)
  #   }
  #   #End of Step 6
  # }

  
  #############End of function fsl_pipe#####################
}


#In development:
if (FALSE) {
  #Flame
  
}