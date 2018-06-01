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
reg.makeconfig<-function(layout=NULL,
                         onset=NULL,
                         offset=NULL,
                         value.variable=NULL,
                         parametric.variable=NULL,
                         additional=NULL)

mat.organize<-function(mat.raw=NULL,configuration=NULL) {
  
  #get required item
  #behaviroal onset time:
  
}


