#QC Clock
#2*pnorm(-abs(as.numeric(scale(voxinfo$voxel_count,center = T,scale = T)))) < 0.05

# file_pattern<-"fMRIEmoClock_.*.txt"
# QC_func<-"clock_qc"
# cfgpath="/Volumes/bek/autopreprocessing_pipeline/BPD_7341/clockBPD7341.cfg"
# preproc.nii.patt="nfswudktm*[0-9]_[0-9].nii.gz"
# hdtemplate="/Volumes/bek/Newtemplate_may18/fsl_mni152/MNI152_T1_2mm_brain.nii"


# lapply(,function(dir){
#   strp<-unlist(strsplit(dir,split = .Platform$file.sep))
#   ID<-strp[length(strp)]
#   resultx<-do.call(proc_func,as.list(list.files(path = dir,pattern = file_pattern,recursive = F,full.names = T)))
#   
# })
qc_getinfo<-function(cfgpath=NULL,ssub_dir=file.path(argu$ssub_outputroot,argu$model.name),
                     qc_var="PxH",ssub_template=argu$ssub_fsl_templatepath,stdspace=argu$templatedir,
                     mask=roi_indx_df$singlemask,tempdir=T,...){
  cfg<-cfg_info(cfgpath = cfgpath)
  tp<-readLines(ssub_template)
  qc_evnum<-unique(as.numeric(gsub("^.*([0-9]+).*$", "\\1", tp[grep(qc_var,tp)])))
  dirs<-list.dirs(path = ssub_dir,recursive = F)
  voxlist<-lapply(dirs, function(dir){
    strp<-unlist(strsplit(dir,split = .Platform$file.sep))
    ID<-strp[length(strp)]
    cfl<-lapply(1:cfg$n_expected_funcruns,function(runnum) {
      if(file.exists(file.path(dir,paste0("run",runnum,"_output.feat")))){
        blkdir<-file.path(dir,paste0("run",runnum,"_output.feat"))
        if (!tempdir) {
          outdirx<-file.path(blkdir,"QC")
          dir.create(path = outdirx,showWarnings = F,recursive = F)
        } else {outdirx<-tempdir()}
        voxvol<-get_voxel_count(cfile = file.path(blkdir,paste0("thresh_zstat",qc_evnum,".nii")),stdsfile = stdspace,
                                intmat = file.path(blkdir,"masktostandtransforms.mat"),mask = mask,outdir = outdirx)
        data.frame(ID=ID,run=runnum,voxel_count=voxvol[1],volume_count=voxvol[2])
      }else {data.frame(ID=ID,run=runnum,voxel_count=NA,volume_count=NA)}
    })
    cfvox<-do.call(rbind,cleanuplist(cfl))
    rownames(cfvox)<-NULL
    totalmaskvox<-voxel_count(cfile = mask)
    cfvox$total_mask_voxel<-totalmaskvox[1]
    cfvox$total_mask_volumes<-totalmaskvox[2]
    return(cfvox)
  })
  
  voxdf<-do.call(rbind,voxlist)
  voxdf$voxsuv_per<-voxdf$voxel_count / voxdf$total_mask_voxel
  motiondf<-get_motion_info(configpath = cfgpath,type = argu$motion_type,argu$motion_threshold,...)
  
  voxinfo<-merge(voxdf,motiondf,by = c("ID","run"),all=T)
  return(voxinfo)
  
}


QC_pipe<-function(cfgpath=NULL,QC_func=NULL,bhav_datapath=NULL,bhav_file_patt=NULL,hdtemplate=NULL,
                  QC_auxdir="/Volumes/bek/QC_fsl",nparalle=NULL,supplylist=NULL,preproc.nii.patt="nfswudktm*[0-9]_[0-9].nii.gz",
                  atlas_name="MNI",atlas_index=c(5,6),...){
  if(is.null(nparalle)) {nparalle<-4}
  cfg<-cfg_info(cfgpath = cfgpath)
  argu<<-gen_model_arg(cfgpath=cfgpath,func.nii.name=preproc.nii.patt,mni_template=hdtemplate,QC_auxdir=QC_auxdir,fullmodel=F,
                     QC_func=QC_func,cfg=cfg)
  if(is.null(supplylist)) {
  if(is.null(bhav_datapath)){stop("No way of generating the list for fsl_pipe, you can custom generate one outside QC_pipe if desires so")}
  dirs<-as.list(list.files(path = bhav_datapath,recursive = F,include.dirs = T,full.names = T))
  names(dirs)<-sapply(dirs,function(n) {
    strsplit(n,.Platform$file.sep)[[1]]->nx
    nx[[length(nx)]]
  })
  supplylist<-lapply(dirs,as.list)
  }
  #run it through fsl pipe step 1:3;
  fsl_pipe(argu = argu,prep.call.func = QC_func,prep.call.allsub = supplylist)
  roi_indx_df<-create_roimask_atlas(atlas_name="MNI",target=c(5,6),outdir=argu$QC_auxdir,
                                    fsl_dir=Sys.getenv("FSLDIR"),volxsize="2mm",type="",
                                    singlemask=T,atlas_readtype = "fsl",...) 
  qc_info<-qc_getinfo(cfgpath=argu$cfgpath,ssub_dir=file.path(argu$ssub_outputroot,argu$model.name),qc_var="QC",
                      ssub_template=argu$ssub_fsl_templatepath,stdspace=argu$templatedir,
                      mask=roi_indx_df$singlemask,tempdir = F)
  qc_info$paradigm<-cfg$paradigm_name
  #qc_info$studyname<-cfg$study
  return(qc_info)
}

