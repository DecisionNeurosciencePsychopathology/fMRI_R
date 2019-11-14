####Preliminary Reprots:

aggregate_motion_info<-function(procpath=NULL,type_name=NULL,output_name=NULL){
  
  rawlist<-system(paste0("find ",procpath," -maxdepth 4 -name motion_info"),intern = T)
  print(rawlist)
  xa<-do.call(rbind,lapply(rawlist,function(dirxa){
    if(file.exists(file.path(dirxa,"fd.txt"))){
      fd<-readLines(file.path(dirxa,"fd.txt"))
    }  else {fd<-NA}
    if(file.exists(file.path(dirxa,"dvars.txt"))){
      dvars<-readLines(file.path(dirxa,"dvars.txt"))
    }  else {dvars<-NA}
    modelname<-basename(dirname(dirname(dirxa)))
    ID<-basename(dirname(dirname(dirname(dirxa))))
    session<-basename(dirname(dirxa))
    volume_num <- 1:length(fd)
    data.frame(fd,dvars,modelname,ID,session,volume_num,stringsAsFactors = F)
  }))
  if(!is.null(type_name)){xa$type = type_name}else {xa$type = NA}
  if(!is.null(output_name)){
    write.csv(xa,file = output_name,row.names = F)
  }
  return(xa)
}