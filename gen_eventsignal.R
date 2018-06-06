#######
##Event Creation and Signal Object Generation for Each subject:

#Overall strcture: 
  #Get to dir, loop over each subject use ID
  #Single sub data proc: use post VBA results and b
  #Assign them to environment as each ID as an list object that contains Event and Signal object which are compatible with pipeline

#Objectives:
  #Read .mat files and see if we need to restructure that

#Compare performance between readMat and read.mat
R.matlab::readMat(filena)->testmat #doesn't even work
rmatio::read.mat(filena)->testmat.2 #Use #C which is way more efficient, although it seems like it comes with error message 

#make a configuration file so organizer knows what to grab

#Let's make sense of the process:
  #1, time bin size; calculated by frequency scaling 
  #2, bandit design matrix? Change and modify it
  
dnpl.makeeprimestruct<-function(rawtext=file.choose(),
                           rawdf=NULL,
                           breakpointname="BreakProc") {
  if (is.null(rawdf)){
    rawdf<-read.csv(rawtext,stringsAsFactors = F,sep = "\t")
  } 
  wherearethebreaks<-grep(breakpointname,rawdf$Procedure)
  rawdf$blocknum<-NA

  for (i in 1:(length(wherearethebreaks)+1)) {
    print(i)
    if (i==1) {
      rawdf$blocknum[i:wherearethebreaks[i]]<-i
      
    } else if (i==(length(wherearethebreaks)+1)) {
      rawdf$blocknum[wherearethebreaks[i-1]:length(rawdf$Procedure)]<-i
     
    } else {
      rawdf$blocknum[wherearethebreaks[i-1]:wherearethebreaks[i]]<-i}
  }
  df<-rawdf[-wherearethebreaks,]
  df$trialnum.block<-unsplit(lapply(split(df$blocknum,df$blocknum), seq_along),df$blocknum)
  return(df)
}


reg.makeconfig<-function(layout=NULL,
                         onset="b$",
                         offset="b$",
                         design=NULL,
                         value.variable=NULL,
                         parametric.variable=NULL,
                         additional=NULL) {
  
}

mat.organize<-function(mat.raw=NULL,configuration=NULL) {
  #get required item
  #behaviroal onset time:
  
  
}


