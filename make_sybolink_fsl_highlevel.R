####################################################
####Make softlink for fsl higher level analysis#####
####Author: Jiazhou Chen ; Last Change: May 2018####
####################################################

#This function can't be run in R Studio anymore because of Mac's inheritence policy
#Alternatively, use command line to open Rstuido, 'open -a Rstudio' so profile will be inherented
#Function also use fsl's FLIRT to create first level flirt matrix if needed; 
#if made somewhere else, soft link to each .feat folder as "masktostandtransforms.mat" matrix must be in FLIRT style

########Start Function:
son.prepare4secondlvl<-function(ssana.path=NULL,preproc.path=NULL,standardbarin.path="/Volumes/bek/Newtemplate_may18/fsl_mni152/MNI152_T1_2mm_brain.nii",
  proc.name=NULL,taskname=NULL,overwrite=TRUE) {
  if (is.null(ssana.path) | is.null(preproc.path) | is.null(standardbarin.path) | is.null(proc.name) | is.null(taskname)){stop("not enough info to run")}
#Step 1: Step up parameters, but it's a function so do it outside please
  #Step 2: Go to find all the feat. directory:
  featlist<-system(paste0("find ",ssana.path," -iname '*.feat' -maxdepth 4 -mindepth 1 -type d"),intern = T)
  #Break them down&take out all old_template_results
  strsplit(featlist,split = "/")->s.featlist
  s.featlist.p<-lapply(s.featlist,function(x) {
    if (any(x %in% "old_template_results")) {
      x<-NULL
    } else {x}
  })
  s.featlist.p[sapply(s.featlist.p,is.null)] <- NULL
  maxlength<-length(s.featlist.p[[1]])
#Step 3: Get all the necessary info from the breakdown list:
  linkmap<-data.frame(id=sapply(s.featlist.p, "[[",maxlength-1), runword=sapply(s.featlist.p, "[[",maxlength))
  linkmap$num<-lapply(strsplit(as.character(linkmap$runword),split =''), 
                             function(x) {suppressWarnings(x[which(!is.na(as.numeric(x)))])})
  linkmap$num<-as.numeric(linkmap$num)
  na.omit(linkmap)->linkmap 
  linkmap$origin<-paste(taskname,linkmap$num,sep = "")

#Step 4, make original directory and destination directory:
  linkmap$destination<-file.path(ssana.path,linkmap$id,linkmap$runword,"reg")
  linkmap$destination_standard<-file.path(ssana.path,linkmap$id,linkmap$runword,"reg_standard")
  linkmap$originplace<-file.path(ssana.path,linkmap$id,linkmap$runword,"masktostandtransforms.mat")

#Step 5, make symb.link for each subject and each directory pair:
  for (i in 1:length(linkmap$id)) {
    if (overwrite) {
      st1<-suppressWarnings(system(paste0("rm -r ", linkmap$destination[i]),intern = T))
      st1_5<-suppressWarnings(system(paste0("rm -r ", linkmap$destination_standard[i]),intern = T))
    }
      st2<-system(paste0("mkdir -m 775 ", linkmap$destination[i]),intern = T)
    if (!file.exists(linkmap$originplace[i]))  {
      system("source ~/.bash_profile")
      system(paste0("/usr/local/ni_tools/fsl/bin/flirt -in ",file.path(ssana.path,linkmap$id[i],linkmap$runword[i],"mask.nii.gz")," -ref ",
                    standardbarin.path," -out /Volumes/bek/neurofeedback/.temp.nii -omat ",linkmap$originplace[i],
                    " -bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12  -interp trilinear"))
    }
    if (length(st2)==0) {
      #Make actual links now....
      #invisible(
      tryCatch(
      system(paste0("ln -sv ",file.path(linkmap$originplace)[i]," ",
             file.path(linkmap$destination,"example_func2standard.mat")[i]),intern = T),
      system(paste0("ln -sv ",file.path(linkmap$originplace)[i]," ",
                    file.path(linkmap$destination,"standard2example_func.mat")[i])),
      system(paste0("ln -sv ",standardbarin.path," ",
                    file.path(linkmap$destination,"standard.nii.gz")[i])), error = function(x) {
                      print("something went wrong")
                    }
      )
      #)
      
    } else if (!attr(st2,"status")==1) {print("already exists, skip")}
  }
  return(linkmap)
  print("DONE")
}

########End Function

#####Actually run it:
test<-son.prepare4secondlvl(
                      ssana.path="/Volumes/bek/neurofeedback/sonrisa1/nfb/ssanalysis/fsl",
                      preproc.path="/Volumes/bek/neurofeedback/sonrisa1/proc",
                      standardbarin.path="/Volumes/bek/Newtemplate_may18/fsl_mni152/MNI152_T1_2mm_brain.nii",
                      proc.name="nfb_proc",
                      taskname<-"nfb",
                      overwrite<-TRUE
                      )
