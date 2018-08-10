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
  load(file.path(argu$ssub_outputroot,argu$model.name,"design.rdata"))
} else {allsub.design<-as.environment(list())}
#Take out people who have already been processed;
if (length(names(allsub.design))>0 & !argu$forcereg) {
  idtodo<-as.character(names(prep.call.allsub)[which(! names(prep.call.allsub) %in% names(allsub.design))])
} else {idtodo<-names(prep.call.allsub)}

if (length(idtodo)>0) {
  for (xid in idtodo) {
    prep.call.list<-prep.call.allsub[[xid]]
    tryCatch(
      {
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
          add.nuisa=argu$ifnuisa,
          assigntoenvir=allsub.design)
      },error=function(x) {paste0(xid,": This person failed regressor generation...go investigate")}
    )
  }
  
  save(list = "allsub.design",file = file.path(argu$ssub_outputroot,argu$model.name,"design.rdata"))
} else {message("NO NEW DATA NEEDED TO BE PROCESSED")}

#End of Step 1
}

#Step 2: #PARALLEL
#Now we do the single sub processing using FSL and the regressor that was generated
stepnow<-stepnow+1
if (is.null(argu$onlyrun) | stepnow %in% argu$onlyrun) {

#let's subset this 
small.sub<-eapply(allsub.design, function(x) {
  list(
    ID=x$ID,
    run_volumes=x$run_volumes,
    regpath=x$regpath,
    preprocID=x$preprocID)
})

#This part takes a long time...Let's paralle it:

clusterjobs<-makeCluster(num_cores,outfile="step2_log.txt",type = "FORK")
#clusterExport(clusterjobs,c("argu","gen_reg","small.sub","get_volume_run",
#                            "cfg_info","change_fsl_template","fsl_2_sys_env",
#                            "feat_w_template","info_to_sysenv"),envir = environment())

NU<-parSapply(clusterjobs,small.sub,function(x) {
  fsl_2_sys_env()
  idx<-x$ID
  for (runnum in 1:length(x$run_volumes)) {
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
      feat_w_template(templatepath = argu$ssub_fsl_templatepath,
                      beg = "ARG_",
                      end = "_END",
                      fsf.path = file.path(argu$regpath,argu$model.name,idx,paste0("run",runnum,"_",argu$model.name,".fsf")),
                      envir = xarg)

    } else {message(paste("ID:",idx,"RUN:",runnum,",already exists,","to re-run, remove the directory."))}
  }
  
})

stopCluster(clusterjobs)

#End of Step 2
}

#Step 3: 
#Now we make the symbolic link for template matching...so they are not misaligned anymore...
stepnow<-stepnow+1
if (is.null(argu$onlyrun) | stepnow %in% argu$onlyrun) {

#This one runs fast enough that it should be fine to not parallel it
cfg<-cfg_info(cfgpath = argu$cfgpath)
prepmap<-son.prepare4secondlvl(
  ssana.path=file.path(argu$ssub_outputroot,argu$model.name),            
  preproc.path=cfg$loc_mrproc_root,                                
  standardbarin.path=argu$templatedir, 
  dir.filter=argu$hig_lvl_path_filter,                                                
  proc.name=cfg$paradigm_name,                                                                         
  taskname=cfg$preprocessed_dirname,                                                                   
  overwrite=argu$ifoverwrite_secondlvl,
  outputmap=TRUE)           
##End of Step 3
}


#Step 4: #PARALLEL
#This starts averaging for each subject:
stepnow<-stepnow+1
if (is.null(argu$onlyrun) | stepnow %in% argu$onlyrun) {
  
  if (!is.null(argu$onlyrun) & !argu$onlyrun %in% 2) {
    if (file.exists(file.path(argu$ssub_outputroot,argu$model.name,"design.rdata"))) {
      load(file.path(argu$ssub_outputroot,argu$model.name,"design.rdata"))
    } else {stop("No design rdata file found")}
    small.sub<-eapply(allsub.design, function(x) {
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

clusterjobs1<-makeCluster(num_cores,outfile="step4_log.txt",type = "FORK")
#clusterExport(clusterjobs1,c("cfg",
#                             "argu",
#                             "small.sub",
#                             "get_volume_run",
#                             "cfg_info",
#                             "change_fsl_template",
#                             "fsl_2_sys_env",
#                             "feat_w_template"),envir = environment())

NU<-parSapply(clusterjobs1,small.sub, function(y) {
  larg<-as.environment(list())
  y$ID->larg$idx
  larg$outputpath<-file.path(argu$ssub_outputroot,argu$model.name,larg$idx,"average")
  larg<-list2env(y$featlist,envir = larg)
  if (!file.exists(paste0(larg$outputpath,".gfeat"))) {
    message(paste0("Initializing gfeat for participant: ",larg$idx))
    feat_w_template(templatepath = argu$gsub_fsl_templatepath,
                    beg = "ARG_",
                    end = "_END",
                    fsf.path = file.path(argu$regpath,argu$model.name,larg$idx,paste0("gfeat_temp",".fsf")),
                    envir = larg)
  } else {message("This person already got average done!")}
  
})
stopCluster(clusterjobs1)
################END of step 4
}

#Step 5: #PARALLEL by function
#Do Group Level 
stepnow<-stepnow+1
if (is.null(argu$onlyrun) | stepnow %in% argu$onlyrun) {
#########Start step group lvl analysis
fsltemplate.GL<-readLines(argu$gsub_fsl_templatepath)
#Start Group Level Analysis:
glvl_all_cope(rootdir=argu$ssub_outputroot,
              outputdir=argu$glvl_outputroot,
              modelname=argu$model.name,
              copestorun=1:as.numeric(gsub(".*?([0-9]+).*", "\\1", fsltemplate.GL[grep("ncopeinputs",fsltemplate.GL)])),
              paralleln = num_cores
)
#End Step 5
}

#Step 6: 
#simple graph
stepnow<-stepnow+1
if (is.null(argu$onlyrun) | stepnow %in% argu$onlyrun) {
library(oro.nifti)
plot_image_all(rootpath=argu$glvl_outputroot,
               templatedir=argu$templatedir,
               model.name=argu$model.name,
               patt="OneSampT_tfce_corrp_tstat1.nii.gz",
               threshold=argu$graphic.threshold,
               outputdir=file.path(argu$glvl_outputroot,argu$model.name),
               colour="red")

#End of Step 6
}


#End of function 
}





#In development:
if (FALSE) {
  #To see if make sense to just create a new 
  splitAt <- function(x, pos) unname(split(x, cumsum(seq_along(x) %in% pos)))
  chuckofev<-fsltemplate[min(grep("EV [0-9]* title",fsltemplate)):grep("# Contrast & F-tests mode",fsltemplate,fixed = T)]
  byev<-splitAt(chuckofev,grep("EV [0-9]* title",chuckofev))
  length(byev)->nev
  
  
  
  
  
}
