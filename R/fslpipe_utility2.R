####This script supports version 2 of the fslpipe package. 
####3 levels:
#subj
#sess
#grp


check_argu<-function(argu=NULL,init=F){
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

make_signal_with_grid<-function(outputdata=NULL,dsgrid=NULL,...) {
  #This new version is cleaner and more efficient than the old one!
  if (is.null(outputdata)) {stop("NO DATA SUPPLIED TO MAKE SIGNAL WITH GRID")}
  sp_grid<-split(dsgrid,grepl("_evt",dsgrid$valuefrom))
  message("###Making signal...###")
  if(!is.null(sp_grid$`TRUE`)){
    #Do the evts regressor signals:
    evtGrid<-sp_grid$`TRUE`
    evtOutlist<-lapply(1:nrow(evtGrid),function(ixb){
      message("making...",evtGrid$name[ixb])
      if (any(is.na(outputdata$event.list[[evtGrid$ename[ixb]]]))) {
        tempdf<-outputdata$event.list[[evtGrid$ename[ixb]]][which(!is.na(outputdata$event.list[[evtGrid$ename[ixb]]]$onset)),]
        tempdf.a<-subset.data.frame(tempdf,select = c("run","trial"))
        tempdf.a$value<-1
        list(value=tempdf.a,event=evtGrid$ename[ixb],normalization="none")
      }else {list(value=1,event=evtGrid$ename[ixb],normalization="none")}
    })
    names(evtOutlist)<-evtGrid$name 
  } else {evtOutlist<-list()}
  
  if(!is.null(sp_grid$`FALSE`)){
    #Do the parametric modulated regressors signals
    paraGrid<-sp_grid$`FALSE`
    if(any(!paraGrid$valuefrom %in% names(outputdata$value))) {stop("Make Signal with Grid error: data supplied does not contain all VALUEFROM values")}
    paraOutlist<-lapply(1:nrow(paraGrid),function(ixa){
      curGrid<-paraGrid[ixa,]
      message("making...",curGrid$name)
      jrz<-list()
      jrz[["event"]]<-curGrid$ename
      jrz[["normalization"]]<-curGrid$norm
      jrz[curGrid$normargu]<-TRUE
      value.df<-data.frame(
        run=outputdata$event.list[[jrz$event]]$run,
        trial=outputdata$event.list[[jrz$event]]$trial,
        value=outputdata$value[[curGrid$valuefrom]]
      )
      if (!is.na(curGrid$modifier)) {
        if(any(!curGrid$modifier %in% names(outputdata$value))) {stop("Make Signal with Grid error: data supplied does not contain all MODIFIER values")}
        switch (curGrid$style,
                "multiply" = {value.df$value<-as.numeric(outputdata$value[[curGrid$valuefrom]]) * as.numeric(outputdata$value[[curGrid$modifier]])},
                "subset"   = {value.df<-value.df[which(as.logical(outputdata$value[[curGrid$modifier]])) ,]})
      }
      jrz[["value"]]<-value.df
      return(jrz)
    })
    names(paraOutlist)<-paraGrid$name
  } else {paraOutlist<-list()}
  
  allOutlist<-c(paraOutlist,evtOutlist)
  #DONE!
  return(allOutlist)
  
}

get_preproc_info<-function(id=NULL,cfg=argu$cfg,reg.nii.name=argu$funcimg.namepatt,reg.nii.runnum=argu$funcimg.runnumpatt){
  
  expectdir<-file.path(cfg$loc_mrproc_root,id,cfg$preprocessed_dirname)
  nopreproc<-!dir.exists(expectdir)
  if(!nopreproc){
    mrfiles<-system(paste0("find -L ",expectdir," -iname ",reg.nii.name," -maxdepth 3 -mindepth 1"),intern = T)
    yesmrfiles<-length(mrfiles)>0
    if(yesmrfiles){
      asexpectedrun<-length(mrfiles) == cfg$n_expected_funcruns
      runlength <- unname(sapply(mrfiles, function(x) { 
        fslinfoout<-system(paste("fslinfo",x),intern = T)
        fslinfoout<-do.call(rbind,strsplit(gsub("\\s+","REMOVETHIS",fslinfoout),"REMOVETHIS"))
        as.numeric(fslinfoout[which(fslinfoout[,1]=="dim4"),2])
      }))
      runnum<-gsub(reg.nii.runnum,"\\1",basename(mrfiles))
      
      nuisance<-rep(NA,length(mrfiles))
      nuisance[which(file.exists(file.path(dirname(mrfiles),"nuisance_regressors.txt")))]<-file.path(dirname(mrfiles),"nuisance_regressors.txt")[which(file.exists(file.path(dirname(mrfiles),"nuisance_regressors.txt")))]
      
    } else {
      message("No fMRI file found for ID: ",id)
      asexpectedrun<-NA
      mrfiles<-NA
      nuisance<-NA
    }
  } else {
    message("No preprocessing folder found for ID: ",id)
    mrfiles<-NA
    asexpectedrun<-NA
    nuisance<-NA
  }
  
  return(data.frame(ID=id,mrfiles=mrfiles,runnum=runnum,runlength=runlength,nuisance,asexpectedrunnum=asexpectedrun,stringsAsFactors = F))
}

prep_session_lvl<-function(subj_rootpath=NULL,subj_folderreg=NULL,template_brainpath=NULL,overwrite=T) {
  if (is.null(subj_rootpath) | is.null(template_brainpath) ){stop("not enough info to run")}
  featlist<-system(paste0("find ",subj_rootpath," -iname ",subj_folderreg," -maxdepth 4 -mindepth 1 -type d"),intern = T)
  if(overwrite){unlink(file.path(featlist,"reg"),recursive = T,force = T)}
  NXU<-lapply(file.path(featlist,"reg"),dir.create,showWarnings = F,recursive = T)
  NXU<-file.symlink(from = file.path(featlist,"example_func.nii.gz"),to = file.path(featlist,"reg","example_func.nii.gz"))
  NXU<-file.symlink(from = file.path(featlist,"example_func.nii.gz"),to = file.path(featlist,"reg","example_func2standard.nii.gz"))
  fslpipe::fsl_2_sys_env()
  fsldir<-Sys.getenv("FSLDIR")
  NXU<-file.symlink(from = file.path(fsldir,"/etc/flirtsch/ident.mat"),to = file.path(featlist,"reg","example_func2standard.mat"))
  NXU<-file.symlink(from = file.path(fsldir,"/etc/flirtsch/ident.mat"),to = file.path(featlist,"reg","example_standard2example_func.mat"))
  NXU<-file.symlink(from = template_brainpath,to = file.path(featlist,"reg","standard.nii"))
  return(featlist)
}

pasteFSF<-function(fsfvari="",value="",addComment=NULL,quotevalue=F,featfile=F){
  if(quotevalue) {value<-paste0("\"",value,"\"")}
  if(featfile) {syx<-"feat_files("} else {syx<-"fmri("}
  c(addComment,paste0("set ",syx,fsfvari,")"," ",value))
}

emat_cmat_FSF<-function(ev_mat=NULL,ct_mat=NULL){
  
  Ev_text_re<-unlist(lapply(1:ncol(ev_mat),function(evnum) {
    base_text<-c(paste0("# EV ",evnum),
                 pasteFSF(fsfvari = paste0("evtitle",evnum),value = colnames(ev_mat)[evnum],addComment =NULL,quotevalue = T),
                 pasteFSF(fsfvari = paste0("shape",evnum),value = 2,addComment =NULL,quotevalue = F),
                 pasteFSF(fsfvari = c(paste0("convolve",evnum),paste0("convolve_phase",evnum),paste0("tempfilt_yn",evnum),paste0("deriv_yn",evnum)
                 ),value = 0,addComment =NULL,quotevalue = F),
                 pasteFSF(fsfvari = paste0("custom",evnum),value = "dummy",addComment =NULL,quotevalue = T)
    )
    ortho_text<-unlist(lapply(0:ncol(ev_mat),function(nc){pasteFSF(fsfvari = paste0("ortho",evnum,".",nc),value = 0,addComment =NULL,quotevalue = F)}))
    input_text<-unlist(lapply(1:nrow(ev_mat),function(ns){
      pasteFSF(fsfvari = paste0("evg",ns,".",evnum),value = ev_mat[ns,evnum],addComment =NULL,quotevalue = F)
    }))
    return(c(base_text,ortho_text,input_text))
  })
  )
  Ev_text<-c(Ev_text_re,pasteFSF(fsfvari = c("evs_orig","evs_real"),value = ncol(ev_mat),addComment =NULL,quotevalue = F),
             pasteFSF(fsfvari = "evs_vox",value = 0,addComment =NULL,quotevalue = F)
  )
  #Contrast
  Contrast_text<-c(unlist(lapply(1:nrow(ct_mat),function(ctnum){
    c(paste0("# Contrast ",ctnum),
      pasteFSF(fsfvari = paste("conpic_real",ctnum,sep = "."),value = 1,addComment =NULL,quotevalue = F),
      pasteFSF(fsfvari = paste("conname_real",ctnum,sep = "."),value = rownames(ct_mat)[ctnum],addComment =NULL,quotevalue = T),
      unlist(lapply(1:ncol(ct_mat),function(ny){pasteFSF(fsfvari = paste0("con_real",ctnum,".",ny),value = ct_mat[ctnum,ny],addComment =NULL,quotevalue = F)})),
      pasteFSF(fsfvari = paste0("conmask",ctnum,"_",which(!1:nrow(ct_mat) %in% ctnum)),value = 0,addComment ="##F-Test Variables",quotevalue = F)
    )
  })),
  pasteFSF(fsfvari = c("con_mode_old","con_mode"),value = "real",addComment = "######### Display images for contrast_real",quotevalue = F),
  pasteFSF(fsfvari = c("ncon_orig","ncon_real"),value = nrow(ct_mat),addComment = "######### number of contrasts",quotevalue = F),
  pasteFSF(fsfvari = c("nftests_orig","nftests_real"),value = 0,addComment = "######### number of F tests",quotevalue = F)
  )
  return(c(Ev_text,Contrast_text))
}


gen_fsf_highlvl<-function(proc_ls_fsf=NULL,flame_type = 3, thresh_type = 3,z_thresh = 2.3, p_thresh = 0.05,covariate_names=c("SUBJMEAN"),
                          Pairred_Group=FALSE, custom_evmat=NULL,custom_ctmat=NULL,
                          template_brain = "/Volumes/bek/Newtemplate_may18/fsl_mni152/MNI152_T1_2mm_brain.nii",lowlvlcopenum=NULL,overwrite=F,
                          fsltemplate=readLines("/Volumes/bek/helper_scripts/fsl_pipe/templates/fsl_flame_general_adaptive_template.fsf")){
  if(length(fsltemplate)<1) {stop("No template provided.")}
  Head_text<-c(
    pasteFSF(fsfvari = "thresh",value = thresh_type,
             addComment = "# Thresholding \n # 0 : None \n # 1 : Uncorrected \n# 2 : Voxel \n # 3 : Cluster \n",quotevalue = F),
    pasteFSF(fsfvari = "prob_thresh",value = p_thresh,addComment = "# P threshold",quotevalue = F),
    pasteFSF(fsfvari = "mixed_yn",value = flame_type,
             addComment = "# Higher-level modelling # 3 : Fixed effects # 0 : Mixed Effects: Simple OLS # 2 : Mixed Effects: FLAME 1 # 1 : Mixed Effects: FLAME 1+2",
             quotevalue = F),
    pasteFSF(fsfvari = "z_thresh",value = z_thresh,addComment = "# Z threshold",quotevalue = F),
    pasteFSF(fsfvari = "overwrite_yn",value = 0,addComment = "# Z threshold",quotevalue = F), #We can't really overwrite because it causes problem in higher levels.
    pasteFSF(fsfvari = "conmask1_1",value = 0,addComment = "# Do contrast masking at all?",quotevalue = F),
    pasteFSF(fsfvari = "regstandard",value = template_brain,addComment = "# Standard image",quotevalue = T)
    # #registration tri
    # pasteFSF(fsfvari = "reginitial_highres_yn",value = reg2initial,addComment = "# Registration to initial structural",quotevalue = F),
    # pasteFSF(fsfvari = "reghighres_yn",value = reg2main,addComment = "# Registration to main structural",quotevalue = F),
    # pasteFSF(fsfvari = "regstandard_yn",value = reg2standard,addComment = "  # Registration to standard image?",quotevalue = F)
    
  )
  # xaj<-ls()
  # save(list = xaj,file = "~/debug_fsl_3.rdata")
  #split info into single fsf
  alldf<-do.call(rbind,lapply(proc_ls_fsf,function(gvar_cope_df){
    if(any(!covariate_names %in% names(gvar_cope_df))){stop("One or more vaariables to run is not included in the input data frame")}
    if(flame_type %in% c(1,2)){f_text<-"FLAME"}else{f_text<-"HIGHER LEVEL"}
    if(is.null(gvar_cope_df$Group_Membership)) {gvar_cope_df$Group_Membership<-1}
    
    
    num_lowerlvl<-unique(as.numeric(sapply(lapply(sapply(gvar_cope_df$PATH,list.files,pattern = "design.con",full.name=T),readLines),function(x){
      as.numeric(gsub("\\D", "", x[grepl("/NumContrasts",x)]))
    })))
    if(is.null(lowlvlcopenum)){
      if(length(num_lowerlvl)!=1 | is.na(num_lowerlvl)){stop("unable to detect lower level contrast number. specify please.")}else{lowlvlcopenum<-num_lowerlvl}
    }
    #gvar_cope_df <- merge(single_fsf,input_df,by = "ID",all.x = T)
    if(length(unique(gvar_cope_df$OUTPUTPATH))>1){stop("Incorrect output path length")}
    if(file.exists( file.path(unique(gvar_cope_df$OUTPUTPATH),paste0(unique(gvar_cope_df$NAME),".gfeat"),"cope1.feat","stats","zstat1.nii.gz") )){
      if(overwrite){
        message("For IDs: ",paste(unique(gvar_cope_df$ID),collapse = ", "),
                "\n","Found ",f_text," for '",unique(gvar_cope_df$NAME),"' and overwrite is set to TRUE Will REMOVE & RE-RUN",
                "\n")
        unlink(file.path(unique(gvar_cope_df$OUTPUTPATH),paste0(unique(gvar_cope_df$NAME),".gfeat")),recursive = T,force = T)
      } else {
        message("Found ",f_text," for '",unique(gvar_cope_df$NAME),"' and overwrite is set to FALSE. Will skip.",
                "\n","IDs: ",paste(unique(gvar_cope_df$ID),collapse = ", "),"\n")
        return(NULL)
      }
    }
    if(any(is.na(gvar_cope_df))) {
      message("Found NA in the entry, the whole data point will be removed. If wish to include, change the NA in the input data frame to 0")
      gvar_cope_df<-na.omit(gvar_cope_df)
    }
    
    
    
    if(!is.null(custom_ctmat) && !is.null(custom_evmat)) {
      message("Skipping EV and Contrast making. Custom versions supplied.")
      ev_mat<-custom_ctmat;ct_mat<-custom_ctmat
    } else {
      #Intercept only
      if(length(covariate_names)==1) {
        #Do one sample here:
        ev_mat<-as.matrix(gvar_cope_df[covariate_names])
        ct_mat<-  diag(x = 1,nrow = ncol(gvar_cope_df[covariate_names]),ncol = ncol(gvar_cope_df[covariate_names]))
        colnames(ct_mat)<-colnames(ev_mat)
        rownames(ct_mat)<-covariate_names
      } else if(Pairred_Group){
        message("RUNNING MOD: Pairred Group")
        if(is.null(gvar_cope_df$uID)){stop("'uID' variable is required in the group variable data frame. Please make sure it is included.")}
        if(is.null(gvar_cope_df$Group)){stop("'Group' variable is required in the group variable data frame. Please make sure it is included.")}
        
        paired_m<-diag(x = 1,nrow = length(gvar_cope_df$uID),ncol = length(gvar_cope_df$uID))
        colnames(paired_m)<-rownames(paired_m)<-gvar_cope_df$uID
        paired_mx<-lapply(unique(colnames(paired_m)),function(x){
          if(length(which(colnames(paired_m)==x))>1){
            apply(paired_m[,which(colnames(paired_m)==x)],1,sum,na.rm=T)
          }else{return(NULL)}
        })
        names(paired_mx)<-unique(colnames(paired_m))
        paired_my<-do.call(cbind,paired_mx)
        ev_mat<-cbind(ev_mat,paired_my)
        ev_mat<-ev_mat[-which(!apply(paired_my,1,function(x){any(x!=0)})),]
        
        ct_mat<-cbind(ct_mat,matrix(ncol = ncol(paired_my),nrow = nrow(ct_mat),data = 0))
        
        gvar_cope_df<-gvar_cope_df[-which(gvar_cope_df$uID %in% names(which(!apply(paired_my,1,function(x){any(x!=0)})))),]
      } else {
        gvar_cope_df$dummy<-rnorm(nrow(gvar_cope_df))
        formula_dum<-as.formula(paste("dummy~",paste(covariate_names,collapse = "+"),sep = "+"))
        dummy_model <- lm(formula_dum, gvar_cope_df)
        ev_mat<-as.data.frame(model.matrix(dummy_model))
        if("(Intercept)" %in% colnames(ev_mat) && "Intercept" %in% colnames(ev_mat)){
          ev_mat<-ev_mat[-which(colnames(ev_mat) %in% "(Intercept)")]
        }
        
        
        if(!any(attr(dummy_model$terms,"dataClasses") %in% c("character","factor"))){
          ct_mat<-  diag(x = 1,nrow = ncol(gvar_cope_df[covariate_names]),ncol = ncol(gvar_cope_df[covariate_names]))
          colnames(ct_mat)<-colnames(ev_mat)
          rownames(ct_mat)<-covariate_names
        } else {
          factorvar = attr(dummy_model$terms,"term.labels")[which(attr(dummy_model$terms,"dataClasses") %in% c("character","factor"))-1]
          v <- emmeans::emmeans(dummy_model, list(as.formula(paste0("pairwise ~ ",paste(factorvar,collapse = "+")))))
          condmeans <- v[[1]]@linfct #condition means
          #(condmeans)
          rownames(condmeans) <- summary(v[[1]])$group
          contrasts <- v[[2]]@linfct #emo diffs (pairwise)
          rownames(contrasts) <- sub(" - ", "_gt_", summary(v[[2]])$contrast, fixed=TRUE)
          ct_mat<-as.data.frame(contrasts,stringsAsFactors = F)
          if("(Intercept)" %in% colnames(ct_mat) && "Intercept" %in% colnames(ct_mat)){
            ct_mat<-ct_mat[-which(colnames(ct_mat) %in% "(Intercept)")]
            ct_mat[nrow(ct_mat)+1,]<-as.numeric(colnames(ct_mat) %in% "Intercept")
            rownames(ct_mat)[nrow(ct_mat)]<-"Intercept"
          }
          for (xname in names(attr(dummy_model$terms,"dataClasses"))[!names(attr(dummy_model$terms,"dataClasses")) %in% c(factorvar,"dummy","Intercept")]) {
            ct_mat[nrow(ct_mat)+1,]<-as.numeric(colnames(ct_mat) %in% xname)
            rownames(ct_mat)[nrow(ct_mat)]<-xname
          }
        }
        
        
      }
    }
    
    message("Setting up ",f_text," for '",unique(gvar_cope_df$NAME),"'. Number of data point: ",nrow(gvar_cope_df),
            ".\n","For IDs: ",paste(unique(gvar_cope_df$ID),collapse = ", "),"\n",
            "With covariate: ",paste(covariate_names,collapse = ", "))
    
    EvContrast_text<-emat_cmat_FSF(ev_mat = ev_mat, ct_mat = ct_mat)
    
    #Group input
    Groupinput_text<-c(
      pasteFSF(fsfvari = c(paste("copeinput",1:lowlvlcopenum,sep = ".")),value = 1,
               addComment = "# Number of lower-level copes feeding into higher-level analysis;
         # Change lv2_copenum to do more",quotevalue = F),
      pasteFSF(fsfvari = "ncopeinputs", value = lowlvlcopenum, 
               addComment = "######### # ncopeinputs ",quotevalue = F),
      pasteFSF(fsfvari = paste("groupmem",1:nrow(gvar_cope_df),sep = "."),
               value = gvar_cope_df$Group_Membership,addComment = "######### # Group membership for input ",quotevalue = F),
      pasteFSF(fsfvari = c("npts","multiple"),
               value = nrow(gvar_cope_df),addComment = "######### Number of first-level analyses",quotevalue = F),
      pasteFSF(fsfvari = 1:nrow(gvar_cope_df),
               value = gvar_cope_df$PATH,addComment = "######### # 4D AVW data or FEAT directory ",quotevalue = T,featfile = T)
    )
    
    fsf_final<-c(fsltemplate,
                 pasteFSF(fsfvari = "outputdir",value = file.path(unique(gvar_cope_df$OUTPUTPATH),unique(gvar_cope_df$NAME)),addComment = "# Output directory",quotevalue = T),
                 Head_text,Groupinput_text,EvContrast_text)
    gvar_cope_df$FSF_PATH<-file.path(unique(gvar_cope_df$OUTPUTPATH),"fsf_files",paste0(unique(gvar_cope_df$NAME),".fsf"))
    dir.create(dirname(unique(gvar_cope_df$FSF_PATH)),recursive = T,showWarnings = F)
    writeLines(fsf_final,con = unique(gvar_cope_df$FSF_PATH))
    return(gvar_cope_df)
  })
  )
  return(alldf)
}

do_all_first_level<-function(lvl1_datalist=NULL,lvl1_proc_func=NULL,dsgrid=NULL,func_nii_name=NULL,cfg=NULL,proc_id_subs=NULL,model_name=NULL,nprocess=4,forcererun=F,
                             reg_rootpath=NULL,center_values=TRUE,nuisance_types=c("nuisance","motion_par")) {
  ls_out<-lapply(lvl1_datalist,do.call,what=lvl1_proc_func)
  message("The lvl1 proc did not finish for the following participant(s): ",names(ls_out)[sapply(ls_out,is.null)])
  ls_out<-ls_out[!sapply(ls_out,is.null)]
  ls_signals<-lapply(ls_out,make_signal_with_grid,add_taskness=T,dsgrid=dsgrid)
  ls_signals<-lapply(names(ls_signals),function(ID){lsa<-ls_signals[[ID]];lsa$ID<-ID;return(lsa)})
  names(ls_signals)<-names(ls_out)
  #Parallel
  cluster_step1<- parallel::makeCluster(nprocess,outfile="",type = "FORK")
  ls_ds_matrix<-parallel::parSapply(cluster_step1,ls_signals,function(signalx){
    output<-list(ID=signalx$ID)
    ID = output$ID
    signalx$ID<-NULL
    output$regpath<-file.path(reg_rootpath,model_name,ID)
    dir.create(output$regpath,recursive = T,showWarnings = F)
    if(file.exists(file.path(output$regpath,"gendesign_failed")) && !forcererun) {
      system(paste0("echo This person: ",ID," has failed regressor generation previously. Will Skip."))
      return(NULL)}
    if(file.exists(file.path(output$regpath,"design_output.rdata")) && !forcererun) {
      system(paste0("echo Loading: ",ID))
      outx<-new.env()
      load(file.path(output$regpath,"design_output.rdata"),envir = outx)
      if(!is.null(outx$output)){return(outx$output)}
    }
    system(paste0("echo Deconvolving: ",ID))
    tryCatch({
      run_volum<-get_volume_run(id=paste0(ID,proc_id_subs),cfg = cfg,returnas = "numbers",reg.nii.name = func_nii_name)
      output$nuisan<-get_nuisance_preproc(id=paste0(ID,proc_id_subs),
                                          cfg = cfg,
                                          returnas = "data.frame",
                                          dothese=nuisance_types) 
      output$design<-dependlab::build_design_matrix(center_values=center_values,signals = signalx,
                                                    events = ls_out[[ID]]$event.list$allconcat,write_timing_files = c("convolved", "FSL","AFNI"),
                                                    tr=as.numeric(argu$cfg$preproc_call$tr),plot = F,run_volumes = run_volum,
                                                    output_directory = file.path(reg_rootpath,model_name,ID))
    },error=function(e){print(e);writeLines("FAILED",con = file.path(output$regpath,"gendesign_failed"));return(NULL)})
    if(is.null(output$design)){writeLines("FAILED",con = file.path(output$regpath,"gendesign_failed"));return(NULL)}
    if (!is.null(output$nuisan)){
      for (k in 1:length(output$nuisan)) {
        write.table(as.matrix(output$nuisan[[k]]),file.path(reg_rootpath,model_name,ID,
                                                            paste0("run",k,"_nuisance_regressor_with_motion.txt")),
                    row.names = F,col.names = FALSE)
      }}
    
    output$heatmap<-make_heatmap_with_design(output$design)
    output$volume<-run_volum
    output$preprocID<-paste0(ID,proc_id_subs)
    save(output,file = file.path(output$regpath,"design_output.rdata"))
    return(output)
  })
  parallel::stopCluster(cluster_step1)                 
  ls_ds_matrix<-ls_ds_matrix[!sapply(ls_ds_matrix,is.null)]  
  message("Total of '",length(ls_ds_matrix),"' participants finished regressor generation: /n",paste(names(ls_ds_matrix),collapse = ", "))
  return(ls_ds_matrix)
}


get_pbs_default<-function(){
  pbs_default<-list(account="mnh5174_a_g_hc_default",nodes=1,ppn=4,memory=8,walltime="40:00:00",titlecmd="cd $PBS_O_WORKDIR

export G=/gpfs/group/mnh5174/default

module use $G/sw/modules
module load r/3.6.0
module load fsl/6.0.1
module load afni/19.0.26
module load gsl/2.5

ni_tools=\"$G/lab_resources\"

#location of MRI template directory for preprocessing
MRI_STDDIR=\"${ni_tools}/standard\"

#add preprocessing scripts that may be called in this pipeline
PATH=\"${ni_tools}/c3d-1.1.0-Linux-x86_64/bin:${ni_tools}/fmri_processing_scripts:${ni_tools}/fmri_processing_scripts/autopreproc:${ni_tools}/bin:${PATH}\"

export PATH MRI_STDDIR
",morecmd="")
  return(pbs_default)
}
pbs_cmd<-function(account,nodes,ppn,memory,walltime,titlecmd,morecmd,cmd,wait_for=""){
  heading<-c("#!/usr/bin/env sh",
             "",
             ifelse(wait_for != "", paste0("#PBS -W depend=afterok:", wait_for), ""),
             paste0("#PBS -A ",account),
             paste0("#PBS -l nodes=",nodes,":ppn=",ppn),
             paste0("#PBS -l pmem=",memory,"gb"),
             paste0("#PBS -l walltime=",walltime),
             "#PBS -W group_list=mnh5174_collab",
             "#PBS -j oe",
             "#PBS -m n",
             "",titlecmd,morecmd,cmd)
  return(heading)
}



qsub_commands<-function(cmds=NULL,jobperqsub=NULL,workingdir=NULL,tagname="lvl1",qsublimit=30,ppn=4,pbs_args=NULL) {
  dir.create(workingdir,showWarnings = F,recursive = T)
  setwd(workingdir)
  if((length(cmds)/jobperqsub)>qsublimit) {
    jobperqsub = round(length(cmds)/qsublimit,digits = 0)
    message("Maximum qsub job can submit is set to ",qsublimit,". Rejecting the job_per_qsub argument and will use the recalculated value: ",jobperqsub)
  }
  df <- data.frame(cmd=cmds, job=rep(1:(length(cmds)/jobperqsub),each=jobperqsub,length.out=length(cmds)), stringsAsFactors=FALSE)
  sp_df <- split(df,df$job)
  joblist<-c()
  for (ix in 1:length(sp_df)) {
    cmdx<-sp_df[[ix]]
    message("Setting up job#: ",unique(cmdx$job))
    outfile <- paste0(workingdir, "/qsub_",tagname,"_featsep_", basename(tempfile()), ".pbs")
    if(is.null(pbs_args)){pbs_args<-get_pbs_default();}
    pbs_torun$cmd<-cmdx$cmd
    writeLines(do.call(pbs_cmd,pbs_torun),outfile)
    joblist[ix]<-dependlab::qsub_file(outfile)
    cmdx<-NULL
  }
  dependlab::wait_for_job(joblist)
  return(NULL)
}




