#' Load BigWig files into R
#' 
#' @description
#' Load a list of BigWig files into R to prepare data for generating coverage track plot.
#' 
#' @param file_path A list or vector containing a series of BigWig file paths.
#' @param file_name A list or vector containing the sample label of each BigWig file.
#' 
#' @return A data frame contains the position, coverage score and sample label of each genome tile.
load_BigWig <- function(file_path,
                        file_name = NULL){
  
  #unlist parameter
  file_path <- base::unlist(file_path)
  if(!base::is.null(file_name)){
    file_name <- base::unlist(file_name)
  }
  
  #check parameter
  if(!base::is.null(file_name)){
    if(base::length(file_path) != base::length(file_name)){
      base::stop('file_name must match file_path!')
    }
  }
  
  #load BigWig
  bw_table <- base::do.call(what = base::rbind,args = base::lapply(X = 1:(base::length(file_path)),FUN = function(x){
    
    #load single file
    bw_file <- rtracklayer::import.bw(con = x)
    bw_file <- rtracklayer::as.data.frame(bw_file)
    
    #assign sample name
    if(base::is.null(file_name)){
      temp_char <- base::unlist(base::strsplit(x = x,split = '/|.bw|.bigwig|.BigWig',fixed = FALSE))
      idx <- base::which(temp_char != '')
      temp_char <- temp_char[idx]
      bw_file$Sample <- temp_char[base::length(temp_char)]
    }else{
      bw_file$Sample <- file_name[x]
    }
    
    #return
    return(bw_file)
  }))
  
  #return
  return(bw_table)
}