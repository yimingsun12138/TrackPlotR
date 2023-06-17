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
    ggplot2::coord_cartesian(xlim = c(start_site,end_site),ylim = c(0,y_lim),expand = FALSE) + 
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
#' Generate range track plot to visualize cis-elements (features) on genome, like peak region, human accelerated region (HAR) and so on.
#' 
#' @param Ranges A GRanges or GRangesList object that stores the coordinates (ranges) of cis-elements (features) on the genome.
#' @param region Genome region to generate the range track plot, must be a GRanges object.
#' @param collapse_range Whether collapse the range track, default is FALSE.
#' @param col_pal A custom palette used to override coloring for features.
#' @param overlap_col Color for overlapped region of different features if collapse_range is set to TRUE.
#' @param segment_size The thickness (linewidth) of the feature segment.
#' 
#' @return A ggplot object.
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
    idx <- utils::combn(x = base::names(Ranges),m = 2,FUN = NULL,simplify = FALSE)
    overlapped_range <- base::do.call(what = base::c,args = base::lapply(X = idx,FUN = function(x){
      Range_1 <- Ranges[[x[1]]]
      Range_2 <- Ranges[[x[2]]]
      return(IRanges::intersect(x = Range_1,y = Range_2))
    }))
    if(base::length(overlapped_range) == 0){
      overlapped_range <- methods::as(object = base::paste0(chr,':','0-0'),Class = 'GRanges')
    }else{
      overlapped_range <- IRanges::subsetByOverlaps(x = overlapped_range,ranges = region)
      if(base::length(overlapped_range) == 0){
        overlapped_range <- methods::as(object = base::paste0(chr,':','0-0'),Class = 'GRanges')
      }
    }
    overlapped_range <- IRanges::reduce(x = overlapped_range,drop.empty.ranges = FALSE)
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
    ggplot2::coord_cartesian(xlim = c(start_site,end_site),expand = FALSE) + 
    ggplot2::theme_bw() + 
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank())
  
  #color
  if(!base::is.null(col_pal)){
    gg_object <- gg_object + 
      ggplot2::scale_color_manual(values = col_pal)
  }
  
  #return
  return(gg_object)
}

#' Basic function for transcript visualization
#' 
#' @description
#' Generate transcript track plot to visualize GTF-like genome annotation.
#' 
#' @export
transcript_vis_basic <- function(anno,
                                 region,
                                 transcript_width = 0.3,
                                 transcript_color = 'black',
                                 exon_width = 3,
                                 exon_color = '#333399',
                                 CDS_width = 5,
                                 CDS_color = '#333399',
                                 arrow_break = 0.04,
                                 arrow_length = 0.08,
                                 arrow_color = 'darkgrey',
                                 show_name = c('gene','none','transcript'),
                                 name_size = 3){
  
  #check parameter
  if(base::class(anno) != 'GRanges'){
    base::stop('anno must be a GRanges object!')
  }else{
    anno <- rtracklayer::as.data.frame(anno)
  }
  
  if(base::sum(!(c('type','gene_id','gene_name','transcript_id','cluster') %in% base::colnames(anno))) > 0){
    base::stop('anno must contain 5 columns: type, gene_id, gene_name, transcript_id, cluster.')
  }else{
    anno$gene_id <- base::as.character(anno$gene_id)
    anno$gene_name <- base::as.character(anno$gene_name)
  }
  
  if(base::sum(base::is.na(anno$gene_id)) > 0){
    base::stop('NA exist in gene_id!')
  }else{
    idx <- base::which(base::is.na(anno$gene_name))
    anno[idx,"gene_name"] <- anno[idx,"gene_id"]
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
  
  if((arrow_break >= 1) | (arrow_break <= 0)){
    base::stop('The value of arrow_break must be between 0 and 1!')
  }
  
  #subset anno by overlap
  idx <- base::which((anno$seqnames == chr) & (anno$end >= start_site) & (anno$start <= end_site))
  
  if(base::length(idx) == 0){
    transcript_table <- base::data.frame(seqnames = chr,start = 0,end = 0,strand = '+',type = 'transcript',gene_id = 'none',gene_name = 'none',transcript_id = 'none',cluster = 1)
    transcript_table$mid_point <- base::round(x = (transcript_table$start + transcript_table$end)/2,digits = 0)
    exon_table <- base::data.frame(seqnames = chr,start = 0,end = 0,strand = '+',type = 'exon',gene_id = 'none',gene_name = 'none',transcript_id = 'none',cluster = 1)
    CDS_table <- base::data.frame(seqnames = chr,start = 0,end = 0,strand = '+',type = 'CDS',gene_id = 'none',gene_name = 'none',transcript_id = 'none',cluster = 1)
  }else{
    anno <- anno[idx,]
    
    #transcript table
    idx <- base::which(anno$type == 'transcript')
    if(base::length(idx) == 0){
      transcript_table <- base::data.frame(seqnames = chr,start = 0,end = 0,strand = '+',type = 'transcript',gene_id = 'none',gene_name = 'none',transcript_id = 'none',cluster = 1)
    }else{
      transcript_table <- anno[idx,]
    }
    transcript_table$mid_point <- base::round(x = (transcript_table$start + transcript_table$end)/2,digits = 0)
    
    #exon table
    idx <- base::which(anno$type == 'exon')
    if(base::length(idx) == 0){
      exon_table <- base::data.frame(seqnames = chr,start = 0,end = 0,strand = '+',type = 'exon',gene_id = 'none',gene_name = 'none',transcript_id = 'none',cluster = 1)
    }else{
      exon_table <- anno[idx,]
    }
    
    #CDS table
    idx <- base::which(anno$type == 'CDS')
    if(base::length(idx) == 0){
      CDS_table <- base::data.frame(seqnames = chr,start = 0,end = 0,strand = '+',type = 'CDS',gene_id = 'none',gene_name = 'none',transcript_id = 'none',cluster = 1)
    }else{
      CDS_table <- anno[idx,]
    }
  }
  
  #create arrow table
  arrow_table <- base::do.call(what = base::rbind,args = base::lapply(X = 1:(base::nrow(transcript_table)),FUN = function(x){
    single_transcript <- transcript_table[x,,drop = FALSE]
    break_length <- base::round(x = (region@ranges@width * arrow_break),digits = 0)
    arrow_number <- base::floor((single_transcript$end - single_transcript$start)/break_length)
    
    if(arrow_number > 0){
      breaked_transcript <- base::do.call(what = base::rbind,args = base::lapply(X = 1:arrow_number,FUN = function(y){
        single_arrow <- single_transcript
        
        if(single_arrow$strand == '-'){
          single_arrow$start <- single_arrow$end
          single_arrow$end <- single_arrow$start - (break_length * y)
        }else{
          single_arrow$end <- single_arrow$start + (break_length * y)
        }
        
        return(single_arrow)
      }))
      
      return(breaked_transcript)
    }else{
      breaked_transcript <- single_transcript
      
      breaked_transcript$start <- 0
      breaked_transcript$end <- 0
      
      return(breaked_transcript)
    }
  }))
  
  #plot arrow
  gg_object <- ggplot2::ggplot() + 
    ggplot2::geom_segment(data = arrow_table,mapping = ggplot2::aes(x = start,xend = end,y = cluster,yend = cluster),linewidth = transcript_width,color = arrow_color,arrow = grid::arrow(length = grid::unit(x = arrow_length,units = 'inches')))
  
  #plot transcript
  gg_object <- gg_object + 
    ggplot2::geom_segment(data = transcript_table,mapping = ggplot2::aes(x = start,xend = end,y = cluster,yend = cluster),linewidth = transcript_width,color = transcript_color)
  
  #plot exon
  gg_object <- gg_object + 
    ggplot2::geom_segment(data = exon_table,mapping = ggplot2::aes(x = start,xend = end,y = cluster,yend = cluster),linewidth = exon_width,color = exon_color)
  
  #plot CDS
  gg_object <- gg_object + 
    ggplot2::geom_segment(data = CDS_table,mapping = ggplot2::aes(x = start,xend = end,y = cluster,yend = cluster),linewidth = CDS_width,color = CDS_color)
  
  #plot name
  show_name <- show_name[1]
  if(show_name != 'none'){
    if(show_name == 'gene'){
      if(base::sum(base::is.na(transcript_table$gene_name)) > 0){
        base::stop('NA exist in gene_name!')
      }else{
        gg_object <- gg_object + 
          ggplot2::geom_text(data = transcript_table,mapping = ggplot2::aes(x = mid_point,y = cluster,label = gene_name),size = name_size,vjust = 1.5)
      }
    }
    if(show_name == 'transcript'){
      if(base::sum(base::is.na(transcript_table$transcript_id)) > 0){
        base::stop('NA exist in transcript_id!')
      }else{
        gg_object <- gg_object + 
          ggplot2::geom_text(data = transcript_table,mapping = ggplot2::aes(x = mid_point,y = cluster,label = transcript_id),size = name_size,vjust = 1.5)
      }
    }
  }
  
  #modify ggplot object
  gg_object <- gg_object + 
    ggplot2::ylab('Transcript') + 
    ggplot2::xlab(base::paste0(chr,' : ',start_site,' - ',end_site)) + 
    ggplot2::coord_cartesian(xlim = c(start_site,end_site),expand = FALSE) + 
    ggplot2::theme_bw() + 
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank())
  
  #return
  return(gg_object)
}
