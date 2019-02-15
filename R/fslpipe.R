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

###Initializing argu;
argu$cfg<-cfg_info(cfgpath = argu$cfgpath)
argu$dsgrid<-read.table(argu$gridpath,header = T,sep = c(","),stringsAsFactors = F,strip.white = T,skipNul = T)
if(is.null(argu$dsgrid$AddNeg)){argu$dsgrid$AddNeg<-FALSE}
argu$dsgrid$AddNeg<-as.logical(argu$dsgrid$AddNeg)
if(is.null(argu$model.varinames)) {argu$model.varinames<-argu$dsgrid$name}
#Version upgrade safe keeping
if (!exists("centerscaleall",envir = argu)) {
  message("centerscaleall does not exist, using default options: FALSE")
  argu$centerscaleall<-FALSE}
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
#Version upgrade safe keeping: RE: FALME Support
if (!exists("higherleveltype",envir = argu)) {
  message("Higher level modeling type is not specified, will use randomize by default.")
  argu$higherleveltype<-"randomize"
} 
if (!exists("adaptive_ssfeat",envir = argu)) {
  message("Single Subject Level type is not specified, will use non-adaptive version by default.")
  argu$adaptive_ssfeat<-FALSE
} 
if (!exists("xmat",envir = argu)) {
  message("Single Subject Level type is not specified, will use non-adaptive version by default.")
  ogLength<-length(argu$dsgrid$name)
  argu$xmat<-rbind(diag(x=1,ogLength),diag(x=-1,ogLength))
} 


#############STEP 1: Regressor generation #####################
#GENERATE REGRESSOR USING DEPENDLAB PIPELINE:
stepnow<-1
if (is.null(argu$onlyrun) | stepnow %in% argu$onlyrun) {
  #Create the directory if not existed
  dir.create(file.path(argu$ssub_outputroot,argu$model.name),showWarnings = FALSE,recursive = T)
  #load the design rdata file if exists;
  if (file.exists(file.path(argu$ssub_outputroot,argu$model.name,"design.rdata"))) {
    tryCatch({
      load(file.path(argu$ssub_outputroot,argu$model.name,"design.rdata"))
    }, error=function(e) {
      message(paste0("load not successful, have to re-run step 1...message: ",e))
      assign('allsub.design',as.environment(list()),envir = globalenv())
    })
    
  } else {allsub.design<-as.environment(list())}
  #Take out people who have already been processed;
  if (length(names(allsub.design))>0 & !argu$forcereg) {
    idtodo<-as.character(names(prep.call.allsub)[which(! names(prep.call.allsub) %in% names(allsub.design))])
    message(paste("These IDs already has regressors: ",paste(names(allsub.design),collapse=" "),sep=" ",collapse = " "))
  } else {idtodo<-names(prep.call.allsub)}
  
  assign(prep.call.func,get(prep.call.func),envir = argu)
  if (length(idtodo)>0) {
    for (xid in idtodo) {
      prep.call.list<-prep.call.allsub[[xid]]
      tryCatch(
        {
          message(xid)
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
                   convlv_nuisa=argu$convlv_nuisa,
                   centerscaleall=argu$centerscaleall,
                   argu=argu
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
          # centerscaleall=argu$centerscaleall

        },error=function(e) {message("failed regressor generation...go investigate: ",e)}
      )
    }
    
    save("allsub.design",file = file.path(argu$ssub_outputroot,argu$model.name,"design.rdata"))
  } else {message("NO NEW DATA NEEDED TO BE PROCESSED")}
  
#End of Step 1
}

#############Step 2: Single Subject PARALLEL#######################
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
    design=x$design,
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
    xarg$tr<-argu$cfg$preproc_call$tr
    
    if(is.null(argu$ss_zthreshold)) {xarg$zthreshold<-3.2} else {xarg$zthreshold<-argu$ss_zthreshold}
    if(is.null(argu$ss_pthreshold)) {xarg$pthreshold<-0.05} else {xarg$pthreshold<-argu$ss_pthreshold}
    if (!file.exists(paste0(xarg$outputpath,".feat")) ) {
      message(paste0("Initializing feat for participant: ",idx,", and run: ",runnum))
      xarg$volumes<-x$run_volumes[runnum]
      xarg$funcfile<-get_volume_run(id=paste0(idx,argu$proc_id_subs),cfgfilepath = argu$cfgpath,reg.nii.name = argu$func.nii.name,returnas = "path")[runnum]
      xarg$nuisa<-file.path(argu$regpath,argu$model.name,idx,paste0("run",runnum,"_nuisance_regressor_with_motion.txt"))
      if (any(unlist(eapply(xarg,is.na)))) {stop("NA exists in one of the arguments; please double check!")}
      #gen_reg(vmodel=argu$model.varinames,regpath=file.path(argu$regpath,argu$model.name),idx=idx,runnum=runnum,env=xarg,regtype = argu$regtype)
    
      if(argu$adaptive_ssfeat){
        argu$model.varinames<-argu$dsgrid$name
        argu$copenames<-c(argu$model.varinames,paste0(argu$dsgrid$name[which(argu$dsgrid$AddNeg)],"_neg"))
        xarg$evnum<-1:ncol(argu$xmat)
        xarg$copenum<-1:nrow(argu$xmat)
        xarg$maxevnum<-ncol(argu$xmat)
        xarg$maxcopenum<-nrow(argu$xmat)
        argu$maxcopenum<-nrow(argu$xmat)

        for (xy in 1:nrow(argu$xmat)){
          assign(paste0("copemat",xy),value = 1:ncol(argu$xmat),envir = xarg)
          assign(paste0("copetitle",xy),argu$copenames[xy],envir = xarg)
          assign(paste0("cope_lessnum",xy),(1:length(argu$copenames))[-xy],envir = xarg)
          for(xx in 1:ncol(argu$xmat)){
            assign(paste0("copevalue",xy,"_",xx),value = argu$xmat[xy,xx],envir = xarg)
          }
        } 
        for (mv in 1:ncol(argu$xmat)) {
          assign(paste0("evtitle",mv),argu$model.varinames[mv],envir = xarg)
          assign(paste0("orthocombo",mv),paste(mv,(0:length(argu$model.varinames)),sep = "."),envir = xarg)
          assign(paste0("evreg",mv),file.path(file.path(argu$regpath,argu$model.name),idx,paste0("run",runnum,"_",argu$model.varinames[mv],argu$regtype)),envir = xarg)
        }
        #PUT NEW FUNCTION HERE
        ssfsltemp<-rep_within(adptemplate = readLines(argu$ssub_fsl_templatepath),searchenvir = xarg)
      } else {ssfsltemp<-readLines(argu$ssub_fsl_templatepath)}
      cmmd<-feat_w_template(fsltemplate = ssfsltemp,beg = "ARG_",end = "_END",
                      fsfpath = file.path(argu$regpath,argu$model.name,idx,paste0("run",runnum,"_",argu$model.name,".fsf")),
                      envir = xarg,outcommand = T)
      return(cmmd)
    } else {
      message(paste("ID:",idx,"RUN:",runnum,",already exists,","to re-run, remove the directory."))
      return(NULL)}
  }))
  return(cmmd)
}))

cluster_step2<-makeCluster(num_cores,outfile="",type = "FORK")
NX<-parSapply(cluster_step2,step2commands,function(yx) {
          fsl_2_sys_env()
          message(paste0("starting to run /n ",yx))
          tryCatch(
            {system(command = yx,intern = T)
            pb<-txtProgressBar(min = 0,max = 100,char = "|",width = 50,style = 3)
            numdx<-which(yx==step2commands)
            indx<-suppressWarnings(split(1:length(step2commands),1:num_cores))
            pindx<-grep(paste0("\\b",numdx,"\\b"),indx)
            setTxtProgressBar(pb,(which(numdx==indx[[pindx]]) / length(indx[[pindx]]))*100)
            message("DONE")
              }, error=function(e){stop(paste0("feat unsuccessful...error: ", e))}
          )
          
  })
stopCluster(cluster_step2)

#End of Step 2
}



#############Step 3: Prep for Higher Level #######################
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


#############Step 4: Fixed Effect for Single Subject PARALLEL ###############
#This starts averaging for each subject:
stepnow<-stepnow+1
if (is.null(argu$onlyrun) | stepnow %in% argu$onlyrun) {
  
  if (!is.null(argu$onlyrun) & !2 %in% argu$onlyrun) {
    if (file.exists(file.path(argu$ssub_outputroot,argu$model.name,"design.rdata"))) {
      load(file.path(argu$ssub_outputroot,argu$model.name,"design.rdata"))
    } else {stop("No design rdata file found")}
    small.sub<-lapply(allsub.design, function(x) {
      list(
        design=x$design,
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
IDs<-sapply(small.sub,function(j){j$ID})
clusterjobs1<-makeCluster(num_cores,outfile="",type = "FORK")
NU<-parSapply(clusterjobs1,small.sub, function(y) {
  larg<-as.environment(list())
  y$ID->larg$idx
  larg$outputpath<-file.path(argu$ssub_outputroot,argu$model.name,larg$idx,"average")
  larg<-list2env(y$featlist,envir = larg)
  larg$templatedir<-argu$templatedir
  if (argu$adaptive_gfeat) {
    larg$maxrunnum<-1:length(y$featlist)
    ssfsltemp<-readLines(argu$ssub_fsl_templatepath)
    larg$maxcopenum<-1:nrow(argu$xmat)
    
    #PUT NEW FUNCTION HERE
    studyfsltemp<-adopt_feat(adptemplate_path = argu$gsub_fsl_templatepath,searenvir=larg,firstorder=F)
    larg$maxrunnum<-max(larg$maxrunnum)
    larg$maxcopenum<-max(larg$maxcopenum)
    
  } else {
  studyfsltemp<-readLines(argu$gsub_fsl_templatepath)
  larg$maxcopenum<-1:max(as.numeric(gsub(".*?([0-9]+).*", "\\1", ssfsltemp[grep("# Title for contrast",ssfsltemp)])))
  }
  if (!file.exists(paste0(larg$outputpath,".gfeat"))) {
    message(paste0("Initializing gfeat for participant: ",larg$idx))
    feat_w_template(fsltemplate = studyfsltemp,
                    beg = "ARG_",
                    end = "_END",
                    fsfpath = file.path(argu$regpath,argu$model.name,larg$idx,paste0("gfeat_temp",".fsf")),
                    envir = larg)
    pb<-txtProgressBar(min = 0,max = 100,char = "|",width = 50,style = 3)
    numdx<-which(larg$idx==IDs)
    indx<-suppressWarnings(split(1:length(IDs),1:num_cores))
    pindx<-grep(paste0("\\b",numdx,"\\b"),indx)
    setTxtProgressBar(pb,(which(numdx==indx[[pindx]]) / length(indx[[pindx]]))*100)
  } else {message("This person already got average done!")}
  
})
stopCluster(clusterjobs1)
################END of step 4
}



#############Step 5a: Higher Level (Randomize) ##PARALLEL by function#########

stepnow<-stepnow+1
if (is.null(argu$onlyrun) | stepnow %in% argu$onlyrun) {
  ssfsltemp<-readLines(argu$ssub_fsl_templatepath)
  
  #Randomize here:
  if(argu$higherleveltype=="randomize"){
  onesamplet_pergroup<-F
  pairedtest<-F
  if (is.null(argu$group_id_sep) | !exists('group_id_sep',envir = argu)) {argu$group_id_sep<-""} 
  if (is.null(argu$cluster_thresh) | !exists('cluster_thresh',envir = argu)) {argu$cluster_thresh<-3} 
  if (is.null(argu$whichttest) | !exists('whichttest',envir = argu)) {argu$whichttest<-"onesample"}
  if ("onesample" %in% argu$whichttest) {onesamplet_pergroup<-T}
  if ("paired" %in% argu$whichttest) {pairedtest<-T}
  #To adopt the new chnages made in adaptive ss 
  if(argu$adaptive_ssfeat) {maxcopenum<-1:length(argu$model.varinames)} else {
    maxcopenum<-1:max(as.numeric(gsub(".*?([0-9]+).*", "\\1", ssfsltemp[grep("# Title for contrast",ssfsltemp)])))
    }
  #Start Group Level Analysis:
  glvl_all_cope(rootdir=argu$ssub_outputroot,
                outputdir=argu$glvl_outputroot,
                modelname=argu$model.name,
                grp_sep=argu$group_id_sep,
                onesamplet_pergroup=onesamplet_pergroup,
                pairedtest=pairedtest,
                copestorun=maxcopenum,
                thresh_cluster_mass=argu$thresh_cluster_mass,
                thresh_cluster_extent=argu$thresh_cluster_extent,
                pvalue=argu$randomize_p_threshold,
                usethesetest=argu$randomize_thresholdingways,
                ifDeMean=argu$randomize_demean,
                paralleln = num_cores)
  #Use for debugging:
  # rootdir=argu$ssub_outputroot
  # outputdir=argu$glvl_outputroot
  # modelname=argu$model.name
  # grp_sep=argu$group_id_sep
  # onesamplet_pergroup=onesamplet_pergroup
  # pairedtest=pairedtest
  # thresh_cluster_siz=argu$cluster_thresh
  # copestorun = maxcopenum
  # thresh_cluster_mass=argu$thresh_cluster_mass
  # thresh_cluster_extent=argu$thresh_cluster_extent
  # pvalue=argu$randomize_p_threshold
  # usethesetest=argu$randomize_thresholdingways
  # ifDeMean=argu$randomize_demean
  # paralleln = num_cores
  } else if (argu$higherleveltype=="flame") {
    #Run flame here:
    if(!exists("whichflame",envir = argu)) {
      message("whichflame variable not specified, will use 1+2 by default.")
      argu$whichflame<-"1+2"
    }
    if(argu$whichflame=="1+2"){flametype<-1} else if(argu$whichflame=="1"){flametype<-2} else {stop("Feat's flame can only support '1+2' or '1'")}
    
    #Variables to get:
      #outputpath:
      #
    
    
  } else {stop("Higher level type ",argu$higherleveltype," is supported, only 'randomize' or 'flame' is currently supported")}
  
  #End Step 5
} 


#############Step 6: Simple Graph and Info Extraction ###################

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
if(argu$adaptive_ssfeat){
  xout<-rbind(
  data.frame(copenum=seq(argu$dsgrid$name),copename=(argu$dsgrid$name)),
  data.frame(copenum=seq(from=length(argu$dsgrid$name)+1,along.with = which(argu$dsgrid$AddNeg)),
             copename=paste0(argu$dsgrid$name[which(argu$dsgrid$AddNeg)],"_neg"))
  )
  write.table(xout,file = file.path(argu$glvl_outputroot,argu$model.name,"cope_title_index.txt"),row.names = F)
}else{
write.table(data.frame(copenum=paste0("cope ",as.numeric(gsub(".*?([0-9]+).*", "\\1", ssfsltemp[grep("# Title for contrast_orig",ssfsltemp)]))),
           title=gsub("\"","",gsub(pattern = "[0-9]*) \"",replacement = "",
                  x = gsub(pattern = "set fmri(conname_orig.",replacement = "",
                           x = gsub(pattern = "set fmri(conname_orig.",replacement = "",
                                    x = ssfsltemp[grep("# Title for contrast_orig",ssfsltemp)+1],fixed = T),fixed = T),fixed = F))
),file = file.path(argu$glvl_outputroot,argu$model.name,"cope_title_index.txt"),row.names = F)
}
#End of Step 6
}

#############End of function fsl_pipe#####################
}


#In development:
if (FALSE) {
#Flame
  
}
