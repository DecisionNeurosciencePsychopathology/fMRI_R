####################################################
####Make softlink for fsl higher level analysis#####
####Author: Jiazhou Chen ; Last Change: May 2018####
####################################################

#This function can't be run in R Studio anymore because of Mac's inheritence policy
#NONSENSEEEEEE YOU CAN JUST USE THE fsl_2_sys_env() command!

#Set up fsl enviroment just in case we are using RStudio GUI
#fsl_2_sys_env()
########Start Function:
son.prepare4secondlvl<-function(ssana.path=NULL,    
                                #Single Sub Analysis Folder, should contain ID/*.feat
                                preproc.path=NULL,
                                #Preproc folder, should contain ID/${taskname}/${proc.name}[0-9]
                                standardbarin.path="/Volumes/bek/Newtemplate_may18/fsl_mni152/MNI152_T1_2mm_brain.nii",
                                #Path to standard template, take either .nii or .nii.gz
                                proc.name=NULL,
                                #See preproc path comment, this should always followed by a number 
                                taskname=NULL,
                                #See preproc path comment, this should be the dir that contains multiple runs
                                dir.filter=NULL,
                                #Filter out some of the unwanted folders in ssana.path, important especailly for those contain .feat structure
                                overwrite=TRUE,
                                #Logical Statement, if to overwirte existing link
                                outputmap=FALSE
                                #Logical Statement, if true to return linkmap
                                ) {
  if (is.null(ssana.path) | is.null(preproc.path) | is.null(standardbarin.path) | is.null(proc.name) | is.null(taskname)){stop("not enough info to run")}
#Step 1: Step up parameters, but it's a function so do it outside please
  #Step 2: Go to find all the feat. directory:
  featlist<-system(paste0("find ",ssana.path," -iname '*output.feat' -maxdepth 4 -mindepth 1 -type d"),intern = T)
  #Break them down&take out all old_template_results
  strsplit(featlist,split = "/")->s.featlist
  s.featlist.p<-lapply(s.featlist,function(x) {
    if (any(x %in% dir.filter)) {
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
      #st1_5<-suppressWarnings(system(paste0("rm -r ", linkmap$destination_standard[i]),intern = T))
      st1_9<-suppressWarnings(system(paste0("rm -r ", linkmap$originplace[i]),intern = T))
    }
      st2<-dir.create(showWarnings = F,path = linkmap$destination[i])
    if (!file.exists(linkmap$originplace[i]))  {
      fsl_2_sys_env()
      system(paste0("${FSLDIR}/bin/flirt -in ",file.path(ssana.path,linkmap$id[i],linkmap$runword[i],"mask.nii.gz")," -ref ",
                    standardbarin.path," -out /Volumes/bek/neurofeedback/.temp.nii -omat ",linkmap$originplace[i],
                    " -bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12  -interp trilinear"))
    }
    
      #Make actual links now....
      #invisible(
      if (!file.exists(file.path(linkmap$destination,"example_func2standard.mat")[i])){
      file.symlink(from = file.path(linkmap$originplace)[i],to = file.path(linkmap$destination,"example_func2standard.mat")[i])
      file.symlink(from = file.path(linkmap$originplace)[i],to = file.path(linkmap$destination,"standard2example_func.mat")[i])
      file.symlink(from = standardbarin.path,to = file.path(linkmap$destination,"standard.nii.gz")[i])
      } else {message("meh,already there, if you want to overwirite, do overwrite...")}
     
      #)
  }
  if(outputmap) {return(linkmap)}
  print("DONE")
}

########End Function

#####Actually run it:
#####You can get all the warnings by assign the result of the function to an object, and exam the object: test<-son.prepare4secondlvl(...)
if (FALSE) {

son.prepare4secondlvl(
                      ssana.path="/Volumes/bek/neurofeedback/sonrisa2/con_framing/ssanalysis/fsl",            
                      preproc.path="/Volumes/bek/neurofeedback/sonrisa2/proc",                                
                      standardbarin.path="/Volumes/bek/Newtemplate_may18/fsl_mni152/MNI152_T1_2mm_brain.nii", 
                      dir.filter="SC_FirstLevels_oldTemplate",                                                
                      proc.name="cf",                                                                         
                      taskname<-"cf_proc",                                                                   
                      overwrite<-TRUE                                                                         
                      )

}

# Below are arguments for nfb:

# ssana.path="/Volumes/bek/neurofeedback/sonrisa1/nfb/ssanalysis/fsl",
# preproc.path="/Volumes/bek/neurofeedback/sonrisa1/proc",
# standardbarin.path="/Volumes/bek/Newtemplate_may18/fsl_mni152/MNI152_T1_2mm_brain.nii",
# proc.name="nfb",
# taskname<-"nfb_proc",
# dir.filter=c("old_template_results","old"),
# overwrite<-TRUE
