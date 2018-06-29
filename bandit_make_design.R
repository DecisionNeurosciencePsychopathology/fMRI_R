####Bandit script
source("./gen_eventsignal.R")

if (F) {
  start_time <- Sys.time()
  R.matlab::readMat(filena)->testmat
  #doesn't even work
  rmatio::read.mat(filena)->testmat.2 #Use #C which is way more efficient, although it seems like it comes with error message 
  #Problem with C complier in Macs 
  
  end_time <- Sys.time()
}

ds.bandit<-list(nblock=3, ntrial=300, names.event=c("decision","feedback"),
                duration_feedback.mm=1000,ifkind=c("myst","comp"))

bandit.proc<-function(data.mat=data.mat,design=ds.bandit){
  if (is.null(data.mat)) {stop("NO INPUT")}
  ###########do df structure
  
  #NOTE: ALL THE LOGICAL VARIABLES ARE FALSE FOR ABSENT / TRUE FOR EXIST; 
  
  #C0: Create empty environment to put stuff in so that they don't just float around function space but rather saved
  #Also it'd be easier to grab them in later function or loop; since environments are on search path
  #Also, with the envir2df function, it can be converted to df with identical length! #bsrc::dnpl.envir2df(wking)
  wking<-new.env(parent = emptyenv())
  
  #C1: get myst and comp trials
  checkzx<-lapply(design$ifkind, function(x) {
    logiccheck<-!as.logical(seq_along(data.mat$b$protocol.type))
    logiccheck[grep(x,data.mat$b$protocol.type)]<-TRUE
    logiccheck
  })
  names(checkzx)<-ifkind
  #End 
  
  #C2: Make Choice and feedback sensor
  assign("motor_sign",plyr::mapvalues(as.numeric(data.mat$b$stim.RESP),from = c(7,2,3),to = c("leftindex","rightinfex","rightmiddle"),warn_missing = F),envir = wking)
  assign("hand_reg",plyr::mapvalues(as.numeric(data.mat$b$stim.RESP),from = c(7,2,3),to = c("left","right","right"),warn_missing = F),envir = wking)
  assign("zerort_logical", {as.logical(unlist(data.mat$b$stim.RT)=="0" | unlist(data.mat$b$chosen.stim)=="999")},
         envir = wking)
  assign("choice_logical", {as.logical(zerort.logical | checkzx$comp)}, envir = wking)
  assign("feedback_logical", {as.logical(zerort.logical | checkzx$myst)}, envir = wking)
  #END 
  
  #C3: Make trial number & block number timing variables:
  blocknum<-rep(1:design$nblock,each=(design$ntrial / design$nblock))
  assign("run",blocknum,envir = wking)
  assign("trialbyblock",unsplit(lapply(split(blocknum,blocknum), seq_along),blocknum),envir = wking)
  #Get NA in feedback to nothing;
  drop(data.mat$b$feedback.OnsetTime)->data.mat$b$feedback.OnsetTime
  which(is.na(data.mat$b$feedback.OnsetTime))->napos
  data.mat$b$feedback.OnsetTime[napos]<-(drop(data.mat$b$stim.OnsetTime) + drop(data.mat$b$stim.jitter1))[napos]
  
  glmm<-lapply(split(seq(1:design$ntrial),blocknum), function(x) {
    xrange<-min(x):max(x)
    #time<-data.mat$b$stim.OnsetTime
    return(list(decision_onset=(data.mat$b$stim.OnsetTime[xrange]-data.mat$b$stim.OnsetTime[min(x)]) / 1000,
                decision_end=(data.mat$b$stim.OnsetTime[xrange]-data.mat$b$stim.OnsetTime[min(x)]+data.mat$b$stim.RT[xrange]) /1000,
                feedback_onset=(data.mat$b$feedback.OnsetTime[xrange]-data.mat$b$stim.OnsetTime[min(x)]) / 1000,
                feedback_end=(data.mat$b$feedback.OnsetTime[xrange]-data.mat$b$stim.OnsetTime[min(x)]+design$duration_feedback.mm) /1000#,
                #trial_onset=data.mat$b$stim.OnsetTime[xrange]-data.mat$b$stim.OnsetTime[min(x)],
                #trial_end=data.mat$b$stim.OffsetTime[xrange]-data.mat$b$stim.OffsetTime[min(x)]
    ))
  })
  
  
  tes.en<-new.env(parent = emptyenv())
  assign("lzy",list(),envir = tes.en)
  tes<-lapply(glmm, function(x) {
    # assign("tempx",x,envir = tes.en)
    # lapply(names(x),function(y){
    #   xz<-get("tempx",envir = tes.en)
    #   
    # })
    get("lzy",envir = tes.en)->lzy
    for (y in names(x)) {
      lzy[[y]]<-c(lzy[[y]],x[[y]])
    }
    assign("lzy",lzy,envir = tes.en)
  })
  get("lzy",envir = tes.en)->lzy
  for (i in names(lzy)) {
    assign(i,lzy[[i]],envir = wking)
  }
  
  #Last step, get all the ones from 
  finalist<-list()
  for (nev in design$names.event) { #Loop around 
    objects(envir = wking)[grep(nev,objects(envir = wking))]->whichones
    finalist[[nev]]<-data.frame(
      event=nev,
      onset=get(whichones[grep("onset",whichones)],envir = wking),
      duration=get(whichones[grep("end",whichones)],envir = wking) - get(whichones[grep("onset",whichones)],envir = wking),
      run=wking$run,
      trial=wking$trialbyblock
    )
  }
  
  for (i in 1:length(finalist)) {
    if (i==1) {ktz<-finalist[[i]]} else {
      ktz<-rbind(ktz,finalist[[i]])}
  }
  ktz[ktz$allconcat=="NaN"]<-NA
  #ktz<-na.omit(ktz)
  finalist[["allconcat"]]<-ktz
  
  vba<-recon.array(data.mat$out$suffStat)
  
  #Additional Variables, create here, used in gird:
  
  vba$value.chosen.diff.sigmoid<-(1./(1+exp(as.numeric(unlist(vba["value.chosen.diff"])))))
  vba$nullres<-as.numeric(wking$choice_logical)
  vba$choice.num<-as.numeric(!wking$choice_logical)
  vba$feedback.num<-as.numeric(!wking$feedback_logical)
  vba$motor.left<-as.numeric(wking$hand_reg=="left")
  vba$motor.right<-as.numeric(wking$hand_reg=="right")
  vba$comp.trials<-as.numeric(data.mat$b$comp.index)
  vba$myst.trials<-as.numeric(data.mat$b$myst.index)
  
  
  output<-list(event.list=finalist,output.df=bsrc::dnpl.envir2df(wking),value=vba)
  return(output)
}
#End

if (F) {
  output<-bandit.proc(data.mat = data.mat,design = ds.bandit)
  final<-makesignalwithgrid(outputdata = output,nona = F)
  final.multi->final
  #final.subset->final
  
  
  
  
  test.x<-list(feedback=final$feedback.num,
               decision=final$choice.num,
               unsingedPE=final$signedPEs,
               
               computer=final$computer_trials,
               myst=final$myst_trials,
               rewardMagnitudeFeedbackAligned_MC=final$rewardMagnitudeFeedbackAligned_MC,
               #valueDecisionAligned_diff=final$valueDecisionAligned_diff,
               #valueDecisionAligned_chosen=final$valueDecisionAligned_chosen,
               valueDecisionAligned_chosen_diff_sigmoid=final$valueDecisionAligned_chosen_diff_sigmoid,
               valueFeedbackAligned_chosen=final$valueFeedbackAligned_chosen,
               #valueDecisionAligned=final$valueDecisionAligned,
               stakeFeedbackAligned=final$stakeFeedbackAligned,
               stakeDecisionAligned_MC=final$stakeDecisionAligned_MC,
               chosenPosPEs=final$chosenPosPEs,
               chosenNegPEs=final$chosenNegPEs
  )
  #test.x$pe$add_deriv <- TRUE
  #test.x$ev$add_deriv <- TRUE
  
  
  design<-dependlab::build_design_matrix(events = output$event.list$allconcat, 
                                         signals = test.x,
                                         write_timing_files = c("convolved", "AFNI", "FSL"),
                                         tr=1.0)
  design$collin_convolve$run3$vif
  
  
  design.full.2<-dependlab::build_design_matrix(events = output$event.list$allconcat, 
                                                signals = final,
                                                #write_timing_files = c("convolved", "AFNI", "FSL"),
                                                tr=1.0)
}
