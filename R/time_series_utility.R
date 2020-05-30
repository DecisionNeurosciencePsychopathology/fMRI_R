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

