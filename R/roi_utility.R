##Functions to facilitate ROI EXTRACTION:

gen_4D_from_subj_level <- function(subj_rootdir=NULL,gfeat_name=NULL,step2_cope=1,output_dir = NULL){
  raw<-system(paste0("find ",file.path(subj_rootdir,gfeat_name)," -iname '*.feat' -maxdepth 2 -mindepth 1 -type d"),intern = T)
  #attr(raw,"status")
  strsplit(raw,split = "/") ->raw.split
  df.ex<-data.frame(ID=unlist(lapply(raw.split,function(x) {
    x[grep(gfeat_name,x)-1]
  })),
  COPENUM=unlist(lapply(raw.split,function(x) {
    x[grep(gfeat_name,x)+1]
  })),
  PATH=file.path(raw,"stats",paste0("cope",step2_cope,".nii.gz"))
  ,stringsAsFactors = F)
  
  df.ex$COPENUM<-substr(df.ex$COPENUM,start=regexpr("[0-9]",df.ex$COPENUM),stop = regexpr(".feat",df.ex$COPENUM)-1)
  df.ex$exists <- sapply(df.ex$PATH,file.exists)
  df_sp <- split(df.ex,df.ex$COPENUM)
  
  if(is.null(output_dir)) {output_dir = file.path(gsub("*","",subj_rootdir,fixed = T),"4D_nii")}
  dir.create(path=output_dir,showWarnings = F,recursive = T)
  
  ind_info<-lapply(df_sp,function(dea){
    dea <- dea[dea$exists,]
    if(nrow(dea)<1) {return(NULL)}
    concatcmd<-paste(sep=" ","fslmerge -t",file.path(output_dir,paste0("4DConcat_lvl1_cope",unique(dea$COPENUM),".nii.gz")),paste(dea$PATH,collapse = " "))
    system(concatcmd,intern = F)
    write.csv(dea,file = file.path(output_dir,paste0("4DConcat_lvl1_cope",unique(dea$COPENUM),".csv")),row.names = F)
    return(list(copenum = unique(dea$COPENUM),index_df=dea,
                concat_img = file.path(output_dir,paste0("4DConcat_lvl1_cope",unique(dea$COPENUM),".nii.gz")),
                index_csv = file.path(output_dir,paste0("4DConcat_lvl1_cope",unique(dea$COPENUM),".csv"))
    ))
  })
  ind_info <- cleanuplist(ind_info)
  return(ind_info)
}


extract_roi_masks <- function(concat_img = NULL,ID_seq=NULL,mask_dir = NULL,search_pattern="*.nii.gz",
                              max_ncpu = 8,masks_to_exclude=c("cluster_combined.nii")) {
  all_masks<-list.files(pattern = search_pattern,path = mask_dir,full.names = T,recursive = F,all.files = F,include.dirs = F,ignore.case = T)
  
  ncpu_to_use <- ifelse(length(all_masks)>max_ncpu,max_ncpu,length(all_masks))
  roi_para_fork <- parallel::makeForkCluster(nnodes = ncpu_to_use)
  
  ls_roi <- parallel::parLapply(cl = roi_para_fork,X = all_masks, function(mask) {
    if(basename(mask) %in% masks_to_exclude) {return(NULL)}
    system(paste0("echo running ",basename(mask)))
    cmdx<-paste(sep=" ",
                "fslstats",
                "-t",concat_img,
                "-k",mask,
                "-M")
    
    dfa<-data.frame(as.numeric(system(cmdx,intern = T)))
    dfa[[1]][dfa[[1]]==0] <- NA
    names(dfa)<-gsub(".nii.gz","",basename(mask))
    if(nrow(dfa)<1){
      return(NULL)
    }
    return(dfa)
  })
  
  parallel::stopCluster(roi_para_fork)
  ls_roi <- cleanuplist(ls_roi)
  roivalues<-do.call(cbind,ls_roi)
  roivalues$ID <- ID_seq
  return(roivalues)
}