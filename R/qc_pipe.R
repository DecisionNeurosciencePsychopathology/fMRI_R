#QC Clock
2*pnorm(-abs(as.numeric(scale(voxinfo$voxel_count,center = T,scale = T)))) < 0.05

datapath<-"/Users/jiazhouchen/Box Sync/skinner/data/matlab task data/bpd_clock"
file_pattern<-"fMRIEmoClock_.*.txt"
proc_func<-"clock_qc"
cfgpath="/Volumes/bek/autopreprocessing_pipeline/BPD_7341/clockBPD7341.cfg"
func.nii.name="nfswudktm*[0-9]_[0-9].nii.gz"
mni_template="/Volumes/bek/Newtemplate_may18/fsl_mni152/MNI152_T1_2mm_brain.nii"
# lapply(,function(dir){
#   strp<-unlist(strsplit(dir,split = .Platform$file.sep))
#   ID<-strp[length(strp)]
#   resultx<-do.call(proc_func,as.list(list.files(path = dir,pattern = file_pattern,recursive = F,full.names = T)))
#   
# })





clock_qc<-function(x){
  strp<-unlist(strsplit(x,split = .Platform$file.sep))
  ID<-strp[length(strp)-1]
  xf<-list.files(path = x,pattern = argu$file_pattern,full.names = T)
  dt<-read.table(xf,header = T)
  dt$Trial<-as.numeric(unlist(lapply(split(dt$Trial,dt$Run),seq_along)))
  dt$duration<-as.numeric(dt$RT/1000)
  dt$QC_OUT<-as.numeric(dt$Score)
  finalist<-list(QC=data.frame(event="QC",
                                     onset=dt[[argu$QC_onset]],
                                     duration=dt[["duration"]],
                                     run=dt[["Run"]],
                                     trial=dt[["Trial"]]))
  for (i in 1:length(finalist)) {
    if (i==1) {ktz<-finalist[[i]]} else {
      ktz<-rbind(ktz,finalist[[i]])}
  }
  
  finalist[["allconcat"]]<-ktz
  output<-list(event.list=finalist,output.df=dt,value=dt)
}

argu<-gen_qc_model(cfgpath=cfgpath,func.nii.name="nfswudktm*[0-9]_[0-9].nii.gz",mni_template=mni_template,QC_auxdir="/Volumes/bek/QC_fsl",
                   QC_onset="Clock_Onset",file_pattern="fMRIEmoClock_.*.txt")
dirs<-as.list(list.files(path = datapath,recursive = F,include.dirs = T,full.names = T))
names(dirs)<-sapply(dirs,function(n) {
  strsplit(n,.Platform$file.sep)[[1]]->nx
  nx[[length(nx)]]
  })
dirs<-lapply(dirs,as.list)
#run it through fsl pipe step 1:3;

fsl_pipe(argu = argu,prep.call.func = "clock_qc",prep.call.allsub = dirs)

