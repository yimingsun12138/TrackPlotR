#' Load BigWig files into R
#' 
#' @description
#' Load a list of BigWig files into R to prepare data for generating coverage track plot.
#' 
#' @param file_path A list or vector containing a series of BigWig file paths.
#' @param file_name A list or vector containing the sample label of each BigWig file.
#' 
#' @return A data frame contains the position, coverage score and sample label of each genome tile.
#' 
#' @export
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
    bw_file <- rtracklayer::import.bw(con = file_path[x])
    bw_file <- rtracklayer::as.data.frame(bw_file)
    
    #assign sample name
    if(base::is.null(file_name)){
      temp_char <- base::unlist(base::strsplit(x = file_path[x],split = '/|.bw|.bigwig|.BigWig',fixed = FALSE))
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

#' Basic function for coverage track visualization
#' 
#' @description
#' Generate coverage track plot using GRanges-like data frame, 
#' score column stores the normalized coverage value, 
#' Sample column stores the sample label of each genome tile.
#' 
#' @param coverage_table A GRanges-like data frame that contains the position, coverage score and sample label of each genome tile.
#' @param region Genome region to generate the coverage track plot, must be a GRanges object.
#' @param y_lim The numerical y-axis limit for coverage track plot.
#' @param sample_order The sample order when generating the coverage track.
#' @param col_pal A custom palette used to override coloring for samples.
#' 
#' @return A ggplot object.
#' 
#' @export
coverage_vis_basic <- function(coverage_table,
                               region,
                               y_lim = NULL,
                               sample_order = NULL,
                               col_pal = NULL){
  
  #check parameter
  if(base::sum(!(c('seqnames','start','end','score','Sample') %in% base::colnames(coverage_table))) > 0){
    base::stop('missing columns in coverage_table!')
  }
  
  if(base::class(region) != 'GRanges'){
    base::stop('region must be a GRanges object!')
  }else{
    if(base::length(region) != 1){
      base::stop('only 1 region required!')
    }else{
      chr <- base::as.character(region@seqnames)
      start_site <- region@ranges@start
      end_site <- region@ranges@start + region@ranges@width - 1
    }
  }
  
  #subset coverage table by overlap
  idx <- base::which((coverage_table$seqnames == chr) & (coverage_table$end >= start_site) & (coverage_table$start <= end_site))
  coverage_table <- coverage_table[idx,]
  coverage_table <- tidyr::gather(data = coverage_table,key = 'side',value = 'position',start,end)
  base::gc()
  
  #order sample
  if(!base::is.null(sample_order)){
    coverage_table$Sample <- base::factor(x = coverage_table$Sample,levels = sample_order)
  }
  
  #set y limit
  if(base::is.null(y_lim)){
    y_lim <- base::max(coverage_table$score)
  }else{
    if(base::class(y_lim) != 'numeric'){
      base::stop('y_lim must be numerical!')
    }
  }
  
  #basic plot
  gg_object <- ggplot2::ggplot(data = coverage_table,mapping = ggplot2::aes(x = position,y = score,color = Sample,fill = Sample)) + 
    ggplot2::geom_area(stat = 'identity') + 
    ggplot2::facet_wrap(facets = ~ Sample,strip.position = 'right',ncol = 1) + 
    ggplot2::scale_x_continuous(limits = c(start_site,end_site),expand = c(0,0)) + 
    ggplot2::scale_y_continuous(limits = c(0,y_lim),expand = c(0,0)) + 
    ggplot2::ylab(base::paste0('Grouped Coverage\nRange : 0 - ',base::round(x = y_lim,digits = 3))) + 
    ggplot2::xlab(base::paste0(chr,' : ',start_site,' - ',end_site)) + 
    ggplot2::theme_bw() + 
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   panel.spacing = grid::unit(0,'lines'),
                   axis.ticks.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   strip.text.y = ggplot2::element_text(angle = 0))
  
  #color
  if(!base::is.null(col_pal)){
    gg_object <- gg_object + 
      ggplot2::scale_color_manual(values = col_pal) + 
      ggplot2::scale_fill_manual(values = col_pal)
  }
  
  #return
  return(gg_object)
}

#' Basic function for genome range visualization
#' 
#' @description
#' Generate range track plot to visualize cis-elements on genome, like peak region, human accelerated region (HAR) and so on.
#' 
#' @export
range_vis_basic <- function(Ranges,
                            region,
                            collapse_range = FALSE,
                            col_pal = NULL,
                            overlap_col = 'darkgrey',
                            segment_size = 2){
  
  #check parameter
  if(base::class(region) != 'GRanges'){
    base::stop('region must be a GRanges object!')
  }else{
    if(base::length(region) != 1){
      base::stop('only 1 region required!')
    }else{
      chr <- base::as.character(region@seqnames)
      start_site <- region@ranges@start
      end_site <- region@ranges@start + region@ranges@width - 1
    }
  }
  
  #convert to GRanges list
  if(!base::grepl(pattern = 'GRanges',x = base::class(Ranges),fixed = TRUE)){
    base::stop('Ranges must be a GRanges or GRangesList object!')
  }else{
    if(base::grepl(pattern = 'List$',x = base::class(Ranges),fixed = FALSE)){
      Ranges <- Ranges
    }else{
      Ranges <- S4Vectors::SimpleList(Feature = Ranges)
    }
  }
  
  #assign range list name
  if(base::is.null(base::names(Ranges))){
    base::names(Ranges) <- base::paste('Feature',base::as.character(c(1:(base::length(Ranges)))),sep = '_')
  }
  
  #get overlapped range
  if((collapse_range) & (base::length(Ranges) > 1)){
    overlapped_range <- base::Reduce(f = IRanges::intersect,x = Ranges)
    if(base::length(overlapped_range) == 0){
      overlapped_range <- methods::as(object = base::paste0(chr,':','0-0'),Class = 'GRanges')
    }else{
      overlapped_range <- IRanges::subsetByOverlaps(x = overlapped_range,ranges = region)
      if(base::length(overlapped_range) == 0){
        overlapped_range <- methods::as(object = base::paste0(chr,':','0-0'),Class = 'GRanges')
      }
    }
    S4Vectors::mcols(overlapped_range) <- NULL
    overlapped_range <- rtracklayer::as.data.frame(overlapped_range)
    overlapped_range$track <- 'Feature'
  }
  
  #convert to data frame
  Ranges_table <- base::do.call(what = base::rbind,args = base::lapply(X = base::names(Ranges),FUN = function(x){
    temp_Range <- Ranges[[x]]
    S4Vectors::mcols(temp_Range) <- NULL
    
    #subset Range by overlap
    temp_Range <- IRanges::subsetByOverlaps(x = temp_Range,ranges = region)
    if(base::length(temp_Range) == 0){
      temp_Range <- methods::as(object = base::paste0(chr,':','0-0'),Class = 'GRanges')
    }
    
    temp_Range <- rtracklayer::as.data.frame(temp_Range)
    temp_Range$Feature <- x
    
    #return
    return(temp_Range)
  }))
  
  #order track
  if((collapse_range) & (base::length(Ranges) > 1)){
    Ranges_table$track <- 'Feature'
  }else{
    Ranges_table$track <- base::factor(x = Ranges_table$Feature,levels = base::rev(base::names(Ranges)))
  }
  
  #basic plot
  gg_object <- ggplot2::ggplot() + 
    ggplot2::geom_segment(data = Ranges_table,mapping = ggplot2::aes(x = start,xend = end,y = track,yend = track,color = Feature),linewidth = segment_size)
  
  if((collapse_range) & (base::length(Ranges) > 1)){
    gg_object <- gg_object + 
      ggplot2::geom_segment(data = overlapped_range,mapping = ggplot2::aes(x = start,xend = end,y = track,yend = track),color = overlap_col,linewidth = segment_size)
  }
  
  gg_object <- gg_object + 
    ggplot2::ylab('Feature') + 
    ggplot2::xlab(base::paste0(chr,' : ',start_site,' - ',end_site)) + 
    ggplot2::scale_x_continuous(limits = c(start_site,end_site),expand = c(0,0)) + 
    ggplot2::theme_bw() + 
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank())
  
  #color
  if(!is.null(col_pal)){
    gg_object <- gg_object + 
      ggplot2::scale_color_manual(values = col_pal)
  }
  
  #return
  return(gg_object)
}
