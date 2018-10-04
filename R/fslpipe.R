######FSL PIPE:



fsl_pipe<-function(argu=NULL, #This is the arguments environment, each model should have a different one;
                   prep.call.func="", #This should be a character string that's the name of the prep proc function
                   prep.call.allsub=list(ID=list(data=NULL)) #List of ID list of arguments for prep.call.
                   ) {

###STEP 0: Source necessary scripts:
#devtools::source_url("https://raw.githubusercontent.com/DecisionNeurosciencePsychopathology/fMRI_R/master/dnpl_utility.R")

require("devtools")
if("dependlab" %in% installed.packages()){"GREAT, DEPENDLAB PACK IS INSTALLED"}else{devtools::install_github("PennStateDEPENdLab/dependlab")}
fsl_2_sys_env()

require("parallel")
if (is.null(argu$nprocess)){
  if (detectCores()>12){
    num_cores<-12 #Use 8 cores to minimize burden; if on throndike 
    #Or if you are running this on laptop; whatever cores minus 2; I guess if it's a dual core...let's just don't do that (zero core will not paralle anything)
  } else {num_cores<-detectCores()-2} 
} else {argu$nprocess->num_cores}


#######STEP 1:
#GENERATE REGRESSOR USING DEPENDLAB PIPELINE:
stepnow<-1
if (is.null(argu$onlyrun) | stepnow %in% argu$onlyrun) {
#Create the directory if not existed
dir.create(file.path(argu$ssub_outputroot,argu$model.name),showWarnings = FALSE)
#load the design rdata file if exists;
if (file.exists(file.path(argu$ssub_outputroot,argu$model.name,"design.rdata"))) {
  tryCatch({load(file.path(argu$ssub_outputroot,argu$model.name,"design.rdata"))},
           error=function(e) {
             message(paste0("load not successful, have to re-run step 1...message: ",e))
             assign('allsub.design',as.environment(list()),envir = globalenv())
    })
         
} else {allsub.design<-as.environment(list())}
#Take out people who have already been processed;
if (length(names(allsub.design))>0 & !argu$forcereg) {
  idtodo<-as.character(names(prep.call.allsub)[which(! names(prep.call.allsub) %in% names(allsub.design))])
} else {idtodo<-names(prep.call.allsub)}
  #Version upgrade safe keeping
  if (exists("ifnuisa",envir = argu) & !exists("convlv_nuisa",envir = argu)) {
    message("ifnuisa variable is now depreciated, please use convlv_nuisa to control if the pipeline should convolve nuissance regressor")
    argu$convlv_nuisa<-argu$ifnuisa}
  if (!exists("nuisa_motion",envir = argu)) {
    message("argument nuisa_motion doesn't exist, using default options: nuisance and motion parameters")
    argu$nuisa_motion<-c("nuisance","motion_par")
  }
  if (!exists("motion_type",envir = argu)) {
    message("argument motion_type doesn't exist, using default options: fd")
    argu$motion_type<-"fd"}
  if (!exists("motion_threshold",envir = argu)) {
    message("argument motion_threshold doesn't exist, using default options, (fd 0.9; dvar 20)")
    argu$motion_threshold<-"default"}
  
if (length(idtodo)>0) {
  for (xid in idtodo) {
    prep.call.list<-prep.call.allsub[[xid]]
    tryCatch(
      {
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
          convlv_nuisa=argu$convlv_nuisa
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

      },error=function(e) {message("failed regressor generation...go investigate: ",e)}
    )
  }
  
  save("allsub.design",file = file.path(argu$ssub_outputroot,argu$model.name,"design.rdata"))
} else {message("NO NEW DATA NEEDED TO BE PROCESSED")}

#End of Step 1
}

#Step 2: #PARALLEL
#Now we do the single sub processing using FSL and the regressor that was generated
stepnow<-stepnow+1
if (is.null(argu$onlyrun) | stepnow %in% argu$onlyrun) {

  if (!is.null(argu$onlyrun) & !1 %in% argu$onlyrun) {
    if (file.exists(file.path(argu$ssub_outputroot,argu$model.name,"design.rdata"))) {
      load(file.path(argu$ssub_outputroot,argu$model.name,"design.rdata"))
    } else {stop("No design rdata file found, must re-run step 1")}
  }  
  
#let's subset this 
small.sub<-eapply(allsub.design, function(x) {
  list(
    ID=x$ID,
    run_volumes=x$run_volumes,
    regpath=x$regpath,
    preprocID=x$preprocID)
})

step2commands<-unlist(lapply(small.sub,function(x) {
  idx<-x$ID
  cmmd<-unlist(lapply(1:length(x$run_volumes), function(runnum) {
    xarg<-as.environment(list())
    xarg$runnum<-runnum    
    xarg$outputpath<-file.path(argu$ssub_outputroot,argu$model.name,idx,paste0("run",runnum,"_output"))
    xarg$templatebrain<-argu$templatedir
    if (!file.exists(paste0(xarg$outputpath,".feat")) ) {
      message(paste0("Initializing feat for participant: ",idx,", and run: ",runnum))
      xarg$volumes<-x$run_volumes[runnum]
      xarg$funcfile<-get_volume_run(id=paste0(idx,argu$proc_id_subs),cfgfilepath = argu$cfgpath,reg.nii.name = argu$func.nii.name,returnas = "path")[runnum]
      xarg$nuisa<-file.path(argu$regpath,argu$model.name,idx,paste0("run",runnum,"_nuisance_regressor_with_motion.txt"))
      if (any(unlist(eapply(xarg,is.na)))) {stop("NA exists in one of the arguments; please double check!")}
      gen_reg(vmodel=argu$model.varinames,regpath=file.path(argu$regpath,argu$model.name),idx=idx,runnum=runnum,env=xarg,regtype = argu$regtype)
      cmmd<-feat_w_template(templatepath = argu$ssub_fsl_templatepath,beg = "ARG_",end = "_END",
                      fsfpath = file.path(argu$regpath,argu$model.name,idx,paste0("run",runnum,"_",argu$model.name,".fsf")),
                      envir = xarg,outcommand = T)
      return(cmmd)
    } else {
      message(paste("ID:",idx,"RUN:",runnum,",already exists,","to re-run, remove the directory."))
      return(NULL)}
  }))
  return(cmmd)
}))

cluster_step2<-makeCluster(num_cores,outfile=file.path(argu$ssub_outputroot,argu$model.name,"step2_log.txt"),type = "FORK")
NX<-parSapply(cluster_step2,step2commands,function(yx) {
          fsl_2_sys_env()
          message(paste0("starting to run /n ",yx))
          tryCatch(
            {system(command = yx,intern = T)
            message("done")
              }, error=function(e){stop(paste0("feat unsuccessful...error: ", e))}
          )
          
  })
stopCluster(cluster_step2)

#End of Step 2
}



#Step 3: 
#Now we make the symbolic link for template matching...so they are not misaligned anymore...
stepnow<-stepnow+1
if (is.null(argu$onlyrun) | stepnow %in% argu$onlyrun) {

#This one runs fast enough that it should be fine to not parallel it
cfg<-cfg_info(cfgpath = argu$cfgpath)
prepmap<-prepare4secondlvl(
  ssana.path=file.path(argu$ssub_outputroot,argu$model.name),            
  preproc.path=cfg$loc_mrproc_root,                                
  standardbarin.path=argu$templatedir, 
  dir.filter=argu$hig_lvl_path_filter,                                                
  proc.name=cfg$paradigm_name,                                                                         
  taskname=cfg$preprocessed_dirname,                                                                   
  overwrite=argu$ifoverwrite_secondlvl,
  outputmap=TRUE,
  paralleln = num_cores)           
##End of Step 3
}


#Step 4: #PARALLEL
#This starts averaging for each subject:
stepnow<-stepnow+1
if (is.null(argu$onlyrun) | stepnow %in% argu$onlyrun) {
  
  if (!is.null(argu$onlyrun) & !2 %in% argu$onlyrun) {
    if (file.exists(file.path(argu$ssub_outputroot,argu$model.name,"design.rdata"))) {
      load(file.path(argu$ssub_outputroot,argu$model.name,"design.rdata"))
    } else {stop("No design rdata file found")}
    small.sub<-lapply(allsub.design, function(x) {
      list(
        ID=x$ID,
        run_volumes=x$run_volumes,
        regpath=x$regpath,
        preprocID=x$preprocID)
    })
  }
cfg<-cfg_info(cfgpath = argu$cfgpath)
cfg$n_expected_funcruns->runnum
featlist<-lapply(small.sub,function(x) {
  x$ID->idz
  emp<-list()
  for (runnum in 1:length(x$run_volumes)) {
    emp[[paste0("feat",runnum)]]<-file.path(argu$ssub_outputroot,argu$model.name,idz,paste0("run",runnum,"_output.feat"))
  }
  small.sub[[idz]]$featlist<-emp
  #assign("small.sub",small.sub,envir = globalenv())
  small.sub<<-small.sub
  return(emp)
})

clusterjobs1<-makeCluster(num_cores,outfile=file.path(argu$ssub_outputroot,argu$model.name,"step4_log.txt"),type = "FORK")
NU<-parSapply(clusterjobs1,small.sub, function(y) {
  larg<-as.environment(list())
  y$ID->larg$idx
  larg$outputpath<-file.path(argu$ssub_outputroot,argu$model.name,larg$idx,"average")
  larg<-list2env(y$featlist,envir = larg)
  larg$templatedir<-argu$templatedir
  if (argu$adaptive_gfeat) {
    larg$maxrunnum<-length(y$featlist)
    ssfsltemp<-readLines(argu$ssub_fsl_templatepath)
    larg$maxcopenum<-max(as.numeric(gsub(".*?([0-9]+).*", "\\1", ssfsltemp[grep("# Title for contrast",ssfsltemp)])))
    #PUT NEW FUNCTION HERE
    studyfsltemp<-adopt_gfeat(adptemplate_path = argu$gsub_fsl_templatepath,searenvir=larg)
  } else {studyfsltemp<-readLines(argu$gsub_fsl_templatepath)}
  if (!file.exists(paste0(larg$outputpath,".gfeat"))) {
    message(paste0("Initializing gfeat for participant: ",larg$idx))
    feat_w_template(fsltemplate = studyfsltemp,
                    beg = "ARG_",
                    end = "_END",
                    fsfpath = file.path(argu$regpath,argu$model.name,larg$idx,paste0("gfeat_temp",".fsf")),
                    envir = larg)
  } else {message("This person already got average done!")}
  
})
stopCluster(clusterjobs1)
################END of step 4
}


###############################
#Step 5: #PARALLEL by function#
#######Do Group Level########## 
###############################

stepnow<-stepnow+1
if (is.null(argu$onlyrun) | stepnow %in% argu$onlyrun) {
  
  ssfsltemp<-readLines(argu$ssub_fsl_templatepath)
  
  onesamplet_pergroup<-F
  pairedtest<-F
  if (is.null(argu$group_id_sep) | !exists('group_id_sep',envir = argu)) {argu$group_id_sep<-""} 
  if (is.null(argu$cluster_thresh) | !exists('cluster_thresh',envir = argu)) {argu$cluster_thresh<-3} 
  if (is.null(argu$whichttest) | !exists('whichttest',envir = argu)) {argu$whichttest<-"onesample"}
  if ("onesample" %in% argu$whichttest) {onesamplet_pergroup<-T}
  if ("paired" %in% argu$whichttest) {pairedtest<-T}
  #Start Group Level Analysis:
  glvl_all_cope(rootdir=argu$ssub_outputroot,
                outputdir=argu$glvl_outputroot,
                modelname=argu$model.name,
                grp_sep=argu$group_id_sep,
                onesamplet_pergroup=onesamplet_pergroup,
                pairedtest=pairedtest,
                thresh_cluster_siz=argu$cluster_thresh,
                copestorun=1:max(as.numeric(gsub(".*?([0-9]+).*", "\\1", ssfsltemp[grep("# Title for contrast",ssfsltemp)]))),
                paralleln = num_cores)
  
  
  #Use for debugging:
  # rootdir=argu$ssub_outputroot
  # outputdir=argu$glvl_outputroot
  # modelname=argu$model.name
  # grp_sep=argu$group_id_sep
  # onesamplet_pergroup=onesamplet_pergroup
  # pairedtest=pairedtest
  # thresh_cluster_siz=argu$cluster_thresh
  # copestorun=1:max(as.numeric(gsub(".*?([0-9]+).*", "\\1", ssfsltemp[grep("# Title for contrast",ssfsltemp)])))
  # paralleln = num_cores
  
  #End Step 5
}

########################################
#############Step 6: ###################
######simple graph & extra Info ########
########################################

stepnow<-stepnow+1
if (is.null(argu$onlyrun) | stepnow %in% argu$onlyrun) {
library(oro.nifti)
ssfsltemp<-readLines(argu$ssub_fsl_templatepath)

plot_image_all(rootpath=argu$glvl_outputroot,
               templatedir=argu$templatedir,
               model.name=argu$model.name,
               patt="*_tfce_corrp_tstat1.nii.gz",
               threshold=argu$graphic.threshold,
               colour="red")

#Create cope index; regardless of the paths and stuff, it should be fine...
write.table(data.frame(copenum=paste0("cope ",as.numeric(gsub(".*?([0-9]+).*", "\\1", ssfsltemp[grep("# Title for contrast_orig",ssfsltemp)]))),
           title=gsub("\"","",gsub(pattern = "[0-9]*) \"",replacement = "",
                  x = gsub(pattern = "set fmri(conname_orig.",replacement = "",
                           x = gsub(pattern = "set fmri(conname_orig.",replacement = "",
                                    x = ssfsltemp[grep("# Title for contrast_orig",ssfsltemp)+1],fixed = T),fixed = T),fixed = F))
),file = file.path(argu$glvl_outputroot,argu$model.name,"cope_title_index.txt"),row.names = F)

#End of Step 6
}

######################################################################
#End of function 
}


#In development:
if (FALSE) {
  #Adaptive fsl cope;
  #To see if make sense to just create a new 
  splitAt <- function(x, pos) unname(split(x, cumsum(seq_along(x) %in% pos)))
  chuckofev<-fsltemplate[min(grep("EV [0-9]* title",fsltemplate)):grep("# Contrast & F-tests mode",fsltemplate,fixed = T)]
  byev<-splitAt(chuckofev,grep("EV [0-9]* title",chuckofev))
  length(byev)->nev
  
}
