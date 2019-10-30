####This script supports version 2 of the fslpipe package. 
####3 levels:
#subj
#sess
#grp


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
    if(overwrite){unlink(file.path(featlist,"reg"),recursive = T,force = T);unlink(file.path(featlist,"reg_standard"),recursive = T,force = T)}
  NXU<-lapply(file.path(featlist,"reg"),dir.create,showWarnings = F,recursive = T)
  NXU<-file.copy(from = file.path(featlist,"example_func.nii.gz"),to = file.path(featlist,"reg","example_func.nii.gz"),overwrite = overwrite)
  NXU<-file.copy(from = file.path(featlist,"example_func.nii.gz"),to = file.path(featlist,"reg","example_func2standard.nii.gz"),overwrite = overwrite)
  fslpipe::fsl_2_sys_env()
  fsldir<-Sys.getenv("FSLDIR")
  NXU<-file.copy(from = file.path(fsldir,"/etc/flirtsch/ident.mat"),to = file.path(featlist,"reg","example_func2standard.mat"),overwrite = overwrite)
  NXU<-file.copy(from = file.path(fsldir,"/etc/flirtsch/ident.mat"),to = file.path(featlist,"reg","example_standard2example_func.mat"),overwrite = overwrite)
  if(tools::file_ext(template_brainpath)=="gz"){ex_t=".nii.gz"}else{ex_t=".nii"}
  
  NXU<-file.copy(from = template_brainpath,to = file.path(featlist,"reg",paste0("standard",ex_t)),overwrite = overwrite)
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
    
    if (file.exists( file.path(unique(gvar_cope_df$OUTPUTPATH),paste0(unique(gvar_cope_df$NAME),".gfeat") ) ) ) {
      message("For IDs: ",paste(unique(gvar_cope_df$ID),collapse = ", "),
                "\n","Found ",f_text," folder but not completed, for '",unique(gvar_cope_df$NAME),"' Will REMOVE & RE-RUN",
                "\n")
      unlink(file.path(unique(gvar_cope_df$OUTPUTPATH),paste0(unique(gvar_cope_df$NAME),".gfeat")),recursive = T,force = T)
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
                             reg_rootpath=NULL,center_values=TRUE,nuisance_types=c("nuisance","motion_par"),retry=F) {
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
    if(file.exists(file.path(output$regpath,"gendesign_failed")) && !retry) {
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
    
    output$heatmap<-NA
    output$volume<-run_volum
    output$preprocID<-paste0(ID,proc_id_subs)
    save(output,file = file.path(output$regpath,"design_output.rdata"))
    return(output)
  })
  parallel::stopCluster(cluster_step1)                 
  ls_ds_matrix<-ls_ds_matrix[!sapply(ls_ds_matrix,is.null)]  
  message("Total of [",length(ls_ds_matrix),"] participants has completed regressor generation: /n",paste(names(ls_ds_matrix),collapse = ", "))
  return(ls_ds_matrix)
}


get_pbs_default<-function(){
  pbs_default<-list(account="mnh5174_c_g_sc_default",nodes=1,ppn=4,memory=8,walltime="40:00:00",titlecmd="cd $PBS_O_WORKDIR

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
    pbs_args$cmd<-cmdx$cmd
    writeLines(do.call(pbs_cmd,pbs_args),outfile)
    joblist[ix]<-dependlab::qsub_file(outfile)
    cmdx<-NULL
  }
  dependlab::wait_for_job(joblist)
  return(NULL)
}

feat2afni_single<-function(feat_dir=NULL,include_copestat=T,include_varcope=F,include_auxstats=F,outputdir=NULL,prefix=NULL,AFNIdir=NULL,template_path=NULL,verbos=F){
  #This is a functionalized version of Dr. Michael Hallquists' Rscript. 
  #Original script can be found here:
  #https://github.com/PennStateDEPENdLab/fmri_processing_scripts
  
  #Safe keeping;
  if(is.null(AFNIdir)){AFNIdir<-Sys.getenv("AFNIDIR")}
  if(is.null(AFNIdir) | AFNIdir==""){AFNIdir<-"~/abin"}
  if(!dir.exists(AFNIdir)) {stop("AFNI directory can not be found, please specify.")}
  if(is.null(feat_dir) || !dir.exists(feat_dir) || length(feat_dir)>1){stop("Feat directory is either not supplied, not found or has a length greater than 1.")}
  if(is.null(prefix)){message("Prefix is set to default");prefix="ss"}
  
  if (include_copestat){
    prefix_statsfiles = paste(prefix,"stats",sep = "_")
    prefix_auxfiles = paste(prefix,"auxstats",sep = "_")
    
    zfiles <- list.files(file.path(feat_dir,"stats"), pattern="zstat[0-9]+\\.nii.*", full.names=TRUE)
    
    if(length(zfiles)<1){message("No zstat file identified. Terminating.");return(NULL)}
    
    statnums <- as.numeric(sub(".*zstat(\\d+)\\.nii.*", "\\1", zfiles, perl=TRUE))
    
    design_contrasts <- readLines(file.path(feat_dir,"design.con"))
    contrast_names <- sub("/ContrastName\\d+\\s+([\\w_.]+).*", "\\1", grep("/ContrastName", design_contrasts, value=TRUE), perl=TRUE)
    
    df_zstats<-data.frame(contrast_names=contrast_names,statnums=statnums,z=zfiles,stringsAsFactors = F)
    df_zstats<-df_zstats[order(df_zstats$statnums),]
    df_zstats$coef<-gsub("/zstat","/cope",df_zstats$z)
    df_zstats$var<-gsub("/zstat","/varcope",df_zstats$z)
    
    
    df_stats_melt<-reshape2::melt(df_zstats,id.vars=c("contrast_names","statnums"))
    df_stats_melt$variable<-factor(df_stats_melt$variable,levels = c("coef","z","var"))
    df_stats_melt<-df_stats_melt[order(df_stats_melt$statnums,df_stats_melt$variable),]
    
    nstats <- nrow(df_zstats)
    message("Found ",nstats," stats files.")
    
    #outputdir<-file.path(feat_dir,"afni_stats")
    dir.create(outputdir,showWarnings = F,recursive = T)
    
    if(include_varcope){toget_string<-c("coef","z","var")} else {toget_string<-c("coef","z")}
    
    df_stats_melt <- df_stats_melt[which(as.character(df_stats_melt$variable) %in% toget_string),]
    df_stats_melt$bricknum<-1:nrow(df_stats_melt)
    
    
    tcatcall <- paste("3dTcat -overwrite -prefix", file.path(outputdir,prefix_statsfiles),
                      paste(df_stats_melt$value,collapse = " ")
                      ,sep=" ")
    
    
    refitcall <- paste("3drefit -fbuc", 
                       paste("-substatpar", df_stats_melt$bricknum[which(as.character(df_stats_melt$variable) == "z")], "fizt", collapse=" "), 
                       "-relabel_all_str", paste0("'",paste(df_stats_melt$contrast_names,df_stats_melt$variable,sep = ":",collapse = " "),"'"), 
                       file.path(outputdir,paste0(prefix_statsfiles,"+tlrc"))) 
    
    system(command = paste(AFNIdir,tcatcall,sep = "/"),intern = F,ignore.stdout = !verbos,ignore.stderr = !verbos)
    system(command = paste(AFNIdir,refitcall,sep = "/"),intern = F,ignore.stdout = !verbos,ignore.stderr = !verbos)
    o_statsfile<-file.path(outputdir,paste0(prefix_statsfiles,"+tlrc"))
  } else {
    o_statsfile<-NA
    } 
  
  if (include_auxstats) {
    
    ##read auxiliary files (PEs + error)
    
    auxfiles<-list(
      pe = list.files(file.path(feat_dir,"stats"), pattern="^pe.*\\.nii.*", full.names=TRUE),
      threshz = list.files(path = feat_dir,pattern="^thresh_zstat.*\\.nii.*", full.names=TRUE),
      zfstat = list.files(path = file.path(feat_dir,"stats"), pattern="zfstat.*\\.nii.*", full.names=TRUE)
    )
    if ( file.exists(file.path(feat_dir,"stats","sigmasquareds.nii.gz")) ) {
      auxfiles$sigmasquared <- file.path(feat_dir,"stats","sigmasquareds.nii.gz")
    }
    auxfiles = auxfiles[sapply(auxfiles,length)>0]
    
    aux_df<-do.call(rbind,lapply(names(auxfiles),function(x){
      xa<-data.frame(type=x,file=auxfiles[[x]],stringsAsFactors = F)
      xa$ifZ <- ifelse(x %in% c("threshz","zfstat"),T,F)
      return(xa)
    }))
    aux_df$bricknum<-1:nrow(aux_df)
    aux_df$name<-gsub(".nii.gz","",basename(aux_df$file))
    
    tcat_auxcall <- paste("3dTcat -overwrite -prefix", file.path(outputdir,prefix_auxfiles),paste(aux_df$file,collapse = " "),sep = " ")
    
    refit_auxcall <- paste("3drefit -fbuc", 
                           paste("-substatpar", aux_df$bricknum[which(aux_df$ifZ)], "fizt", collapse=" "), 
                           paste0("-relabel_all_str '", paste(aux_df$name,collapse = " "), "' "),
                           file.path(outputdir,paste0(prefix_auxfiles,"+tlrc")),sep = " "
    )
    
    
    system(paste(AFNIdir,tcat_auxcall,sep = "/"),intern = F,ignore.stdout = !verbos,ignore.stderr = !verbos)
    system(paste(AFNIdir,refit_auxcall,sep = "/"),intern = F,ignore.stdout = !verbos,ignore.stderr = !verbos)
    o_auxstatsfile <- file.path(outputdir,paste0(prefix_auxfiles,"+tlrc"))
  } else {
    o_auxstatsfile <- NA
  }
  
  
  if(is.null(template_path)){
    file.copy(from = template_path,to = file.path(outputdir,"template_brain.nii.gz"),overwrite = T)
  }
  #message("Completed")
  return(list(statsfile = o_statsfile,auxstatsfile=o_auxstatsfile))
}

gfeat2afni <- function(gfeat_dir=NULL,include_varcope=F,copy_subj_cope=F,outputdir=NULL,prefix=NULL,AFNIdir=NULL,template_path=NULL,verbos=F){
  #Safe keeping;
  if(is.null(AFNIdir)){AFNIdir<-Sys.getenv("AFNIDIR")}
  if(is.null(AFNIdir) | AFNIdir==""){AFNIdir<-"~/abin"}
  if(!dir.exists(AFNIdir)) {stop("AFNI directory can not be found, please specify.")}
  if(is.null(gfeat_dir) || !dir.exists(gfeat_dir) || length(gfeat_dir)>1){stop("Feat directory is either not supplied, not found or has a length greater than 1.")}
  if(is.null(prefix)){message("Prefix is set to default");prefix="grp"}
  
  dir.create(outputdir,showWarnings = F,recursive = T)
  
  copedirs <- grep("/cope[0-9]+\\.feat", list.dirs(path=gfeat_dir, full.names=TRUE, recursive=FALSE), value=TRUE, perl=TRUE)
  
  if(length(copedirs)<1) {
    message("Failed to locate any completed feat directory.")
    return(NULL)
  } else {
    message("Located ",length(copedirs)," cope directories")
  }
  allafniout<-sapply(copedirs,function(dxa){
    message("Processing: ",dxa)
    afniout<-suppressMessages(feat2afni_single(feat_dir = dxa,include_copestat = T,include_varcope = include_varcope,include_auxstats = F,outputdir = dxa,
                              prefix = "sfeat",verbos = verbos))$statsfile
    if(is.null(afniout)){return(NULL)}
    copename <- gsub(" ","_",readLines(file.path(dxa, "design.lev"))) #contains the L2 effect name (e.g., clock_onset)
    #afniout <- file.path(dxa, "feat_stats+tlrc")
    briklabels<-system(command = paste(AFNIdir,paste("3dinfo -label", afniout),sep = "/"),intern = T)
    briklabels <- paste(copename, strsplit(briklabels, "|", fixed=TRUE)[[1]], sep="_", collapse=" ")
    
    ##need to add prefix for each cope to make the stats unique
    system(command = paste(AFNIdir,
                           paste0("3drefit -relabel_all_str '", briklabels, "' ", afniout)
                           ,sep = "/"),intern = F,ignore.stdout = !verbos,ignore.stderr = !verbos)
    
    ##for now, eliminate the aux file (now handled by --no_auxstats above)
    #system(paste("rm", file.path(copedirs[d], "feat_aux+tlrc*")))
    
    if (copy_subj_cope) {
      ##filtered_func_data contains the cope from the lower level. Copy to output directory and rename
      dir.create(path = file.path(outputdir,"subj_coef"),showWarnings = F,recursive = T)
      system(command = paste(AFNIdir,
                             paste("3dcopy -overwrite", file.path(dxa, "filtered_func_data.nii.gz"), file.path(outputdir,"subj_coef",paste0(prefix, "_", copename, "_cope.nii.gz")),sep = " ")
                             ,sep = "/"),intern = F,ignore.stdout = !verbos,ignore.stderr = F)
      if(include_varcope){
        system(command = paste(AFNIdir,
                               paste("3dcopy -overwrite", file.path(dxa, "var_filtered_func_data.nii.gz"), file.path(outputdir,"subj_coef",paste0(prefix, "_", copename, "_varcope.nii.gz")),sep = " ")
                               ,sep = "/"),intern = F,ignore.stdout = !verbos,ignore.stderr = F)
      }
      
    }
    
    return(afniout)
  })
  allafniout <- allafniout[which(!sapply(allafniout,is.null))]
  
  #glue together the stats files
  system(command = paste(AFNIdir,
                         paste("3dTcat -overwrite -prefix", file.path(outputdir,prefix), paste(allafniout, collapse=" "))
                         ,sep = "/"),intern = F,ignore.stdout = !verbos,ignore.stderr = !verbos)
  
  if (file.exists(file.path(outputdir,paste0(prefix, "+tlrc.BRIK")))) {
    system(paste0("gzip ",file.path(outputdir,paste0(prefix, "+tlrc.BRIK"))))
  }
  
  #cleanup ingredients of individual cope aggregation
  system(paste("rm", paste0(allafniout, "*", collapse=" ")))
  
  if(is.null(template_path)){
    file.copy(from = template_path,to = file.path(outputdir,"template_brain.nii.gz"),overwrite = T)
  }
  return(file.path(outputdir,paste0(prefix, "+tlrc.BRIK")))
}
