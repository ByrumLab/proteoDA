make_contrasts <- function(file=NULL, design){

  param<-stats<-list()

  if(is.null(file)){file=file.choose()}
  ## check that file exists in specified location. if file exits and is csv/tsv/txt
  file <-gsub("\\./","",file)
  file.dir <- dirname(file);file.dir
  file <- file.path(file.dir, match.arg(basename(file), list.files(file.dir),several.ok=F))
  filext <- tools::file_ext(file);filext
  filext <- match.arg(filext, c("csv","txt","tsv"),several.ok=F)
  if(filext=="txt" | filext=="tsv"){ sep="\t" }
  if(filext=="csv"){ sep=","}

  ## imports contrast file
  contrast.vec <- utils::read.csv(file=file, sep=sep, stringsAsFactors=F, header=F)
  contrast.vec <- contrast.vec[,1]
  contrast.vec <-gsub(" ", "", contrast.vec)

  ## extract groups included in each contrast
  contrast.grps <-extract_contrast_groups(contrast.vec=contrast.vec)

  ## if the groups defined in the contrast file do no match the design
  ## then return the imported contrast.vec.
  if(all(unique(unlist(contrast.grps))%in%colnames(design))==FALSE){
    message("the contrast file imported successfully, however, some of the groups ",
            "included in the contrasts do not match the design matrix.",
            "The contrast matrix will need to be defined manually using the makeContrasts() function.")
    message("contrasts <- makeContrasts(contrasts=contrast.vec, levels=design)\n",
            "colnames(contrasts) <- gsub('=.*','',colnames(contrasts))")
    message("Returning the imported contrast.vec for additional processing ...")
    return(contrast.vec)
  }
  ## define contrasts
  contrasts <- makeContrasts(contrasts=contrast.vec, levels=design)
  colnames(contrasts) <- gsub("=.*","",colnames(contrasts))

  param[["file"]]<-file
  stats[["no_contrasts"]] <-ncol(contrasts)
  title="CONTRASTS"
  logs<-make_log(param=param, stats=stats, title=title, save=TRUE)

  data2=list(contrasts=contrasts, contrast.vec=contrast.vec,param=logs$param, stats=logs$stats)
  return(data2)

}


extract_contrast_groups <-function(contrast.vec){

  contrast.vec <- gsub(" ", "", as.character(contrast.vec)) ## remove blank spaces
  mat   <- matrix(unlist(strsplit(contrast.vec,"=")), ncol=2, byrow=TRUE)
  names <- mat[,1]

  vec <- mat[,2]
  vec <- gsub(" ","",as.character(vec)) ## removes blank spaces
  vec <- gsub("\\/+[[:digit:]]","", as.character(vec))  ## removes division by one digit number (/2)
  vec <- gsub("\\-+[[:digit:]]", "", as.character(vec)) ## removes subtraction by one digit number (-2)
  vec <- gsub("\\++[[:digit:]]", "", as.character(vec)) ## removes addition of one digit number (+2)
  vec <- gsub("\\*+[[:digit:]]", "", as.character(vec)) ## removes multiply by one digit number (*2)
  vec <- gsub("\\+", "-", as.character(vec)) ## replaces + with - (+ -> -)
  vec <- gsub("\\(|\\)", "", as.character(vec)) ## removes parentheses ()
  vec <- gsub("\\/", "-", as.character(vec)) ## replace division sign / with -
  vec <- gsub(" ", "", as.character(vec)) ## removes blank spaces
  vec2 <- strsplit(vec, "-")
  names(vec2) <- names

  vec3 <- NULL
  for(i in 1:length(vec2)){ vec3[[i]] <- unique(vec2[[i]]) }
  names(vec3) <- names

  return(vec3)

}
