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
    bw_file <- rtracklayer::as.data.frame(bw_file,row.names = NULL)
    
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
#' @param region Genome region used to generate the coverage track plot, must be a GRanges object.
#' @param y_lim The numerical y-axis limit for coverage track plot.
#' @param sample_order The sample order when generating the coverage track.
#' @param col_pal A custom palette used to override coloring for samples.
#' 
#' @return A ggplot object.
coverage_vis_basic <- function(coverage_table,
                               region,
                               y_lim = NULL,
                               sample_order = NULL,
                               col_pal = NULL){
  
  #check parameter
  if(!base::all(c('seqnames','start','end','score','Sample') %in% base::colnames(coverage_table))){
    base::stop('missing columns in coverage_table!')
  }
  
  if(base::class(region) != 'GRanges'){
    base::stop('region must be a GRanges object!')
  }else{
    if(base::length(region) != 1){
      base::stop('only 1 region required!')
    }else{
      chr <- base::as.character(region@seqnames)
      start_site <- GenomicRanges::start(region)
      end_site <- GenomicRanges::end(region)
    }
  }
  
  #subset coverage table by overlap
  idx <- base::which((coverage_table$seqnames == chr) & (coverage_table$end >= start_site) & (coverage_table$start <= end_site))
  if(base::length(idx) == 0){
    base::stop('no coverage signal overlaps with provided region!')
  }
  coverage_table <- coverage_table[idx,,drop = FALSE]
  coverage_table <- tidyr::gather(data = coverage_table,key = 'side',value = 'position',start,end)
  base::gc()
  
  #order sample
  if(!base::is.null(sample_order)){
    coverage_table$Sample <- base::factor(x = coverage_table$Sample,levels = sample_order)
  }else{
    sample_order <- base::sort(x = base::unique(coverage_table$Sample),decreasing = FALSE)
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
    ggplot2::ylab(base::paste0('Grouped Coverage\nRange : 0 - ',base::format(x = y_lim,scientific = TRUE,digits = 2))) + 
    ggplot2::xlab(base::paste0(chr,' : ',base::format(x = start_site,scientific = FALSE),' - ',base::format(x = end_site,scientific = FALSE))) + 
    ggplot2::theme_bw() + 
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   panel.spacing = grid::unit(0,'lines'),
                   axis.ticks.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   strip.text.y = ggplot2::element_text(angle = 0))
  
  #color
  if(!base::is.null(col_pal)){
    gg_object <- gg_object + 
      ggplot2::scale_color_manual(values = col_pal,breaks = sample_order) + 
      ggplot2::scale_fill_manual(values = col_pal,breaks = sample_order)
  }else{
    gg_object <- gg_object + 
      ggplot2::scale_color_discrete(breaks = sample_order) + 
      ggplot2::scale_fill_discrete(breaks = sample_order)
  }
  
  #return
  return(gg_object)
}

#' Basic function for genome feature visualization
#' 
#' @description
#' Generate feature track plot to visualize cis-elements on genome, like peak region, human accelerated region (HAR) and so on.
#' 
#' @param Ranges A GRanges or GRangesList object that stores the coordinates of cis-elements (features) on the genome.
#' @param region Genome region used to generate the feature track plot, must be a GRanges object.
#' @param collapse_range Whether collapse the feature track, default is FALSE.
#' @param col_pal A custom palette used to override coloring for features.
#' @param overlap_col Color for overlapped region of different features if collapse_range is set to TRUE.
#' @param segment_size The thickness (linewidth) of the feature segment.
#' 
#' @return A ggplot object.
feature_vis_basic <- function(Ranges,
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
      start_site <- GenomicRanges::start(region)
      end_site <- GenomicRanges::end(region)
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
  
  #assign GRanges list name
  if(base::is.null(base::names(Ranges))){
    base::names(Ranges) <- base::paste('Feature',base::as.character(c(1:(base::length(Ranges)))),sep = '_')
  }
  
  #get overlapped range
  if((collapse_range) & (base::length(Ranges) > 1)){
    idx <- utils::combn(x = base::names(Ranges),m = 2,FUN = NULL,simplify = FALSE)
    overlapped_range <- base::do.call(what = base::c,args = base::lapply(X = idx,FUN = function(x){
      Range_1 <- Ranges[[x[1]]]
      Range_2 <- Ranges[[x[2]]]
      return(GenomicRanges::intersect(x = Range_1,y = Range_2,ignore.strand = TRUE))
    }))
    
    if(base::length(overlapped_range) == 0){
      overlapped_range <- methods::as(object = base::paste0(chr,':','0-0'),Class = 'GRanges')
    }else{
      overlapped_range <- IRanges::subsetByOverlaps(x = overlapped_range,ranges = region,ignore.strand = TRUE)
      if(base::length(overlapped_range) == 0){
        overlapped_range <- methods::as(object = base::paste0(chr,':','0-0'),Class = 'GRanges')
      }else{
        #truncate range
        idx <- base::which(GenomicRanges::start(overlapped_range) < start_site)
        if(base::length(idx) > 0){
          GenomicRanges::start(overlapped_range)[idx] <- start_site
        }
        idx <- base::which(GenomicRanges::end(overlapped_range) > end_site)
        if(base::length(idx) > 0){
          GenomicRanges::end(overlapped_range)[idx] <- end_site
        }
      }
    }
    
    overlapped_range <- GenomicRanges::reduce(x = overlapped_range,drop.empty.ranges = FALSE,ignore.strand = TRUE)
    S4Vectors::mcols(overlapped_range) <- NULL
    overlapped_range <- rtracklayer::as.data.frame(overlapped_range,row.names = NULL)
    overlapped_range$track <- 'Feature'
  }
  
  #convert to data frame
  Ranges_table <- base::do.call(what = base::rbind,args = base::lapply(X = base::names(Ranges),FUN = function(x){
    temp_Range <- Ranges[[x]]
    S4Vectors::mcols(temp_Range) <- NULL
    
    #subset Range by overlap
    temp_Range <- IRanges::subsetByOverlaps(x = temp_Range,ranges = region,ignore.strand = TRUE)
    if(base::length(temp_Range) == 0){
      temp_Range <- methods::as(object = base::paste0(chr,':','0-0'),Class = 'GRanges')
    }else{
      #truncate range
      idx <- base::which(GenomicRanges::start(temp_Range) < start_site)
      if(base::length(idx) > 0){
        GenomicRanges::start(temp_Range)[idx] <- start_site
      }
      idx <- base::which(GenomicRanges::end(temp_Range) > end_site)
      if(base::length(idx) > 0){
        GenomicRanges::end(temp_Range)[idx] <- end_site
      }
    }
    
    temp_Range <- rtracklayer::as.data.frame(temp_Range,row.names = NULL)
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
    ggplot2::xlab(base::paste0(chr,' : ',base::format(x = start_site,scientific = FALSE),' - ',base::format(x = end_site,scientific = FALSE))) + 
    ggplot2::scale_x_continuous(limits = c(start_site,end_site),expand = c(0,0)) + 
    ggplot2::theme_bw() + 
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank())
  
  #color
  if(!base::is.null(col_pal)){
    gg_object <- gg_object + 
      ggplot2::scale_color_manual(values = col_pal,breaks = base::names(Ranges))
  }else{
    gg_object <- gg_object + 
      ggplot2::scale_color_discrete(breaks = base::names(Ranges))
  }
  
  #return
  return(gg_object)
}

#' Basic function for transcript visualization
#' 
#' @description
#' Generate transcript track plot to visualize GTF-like gene annotation.
#' 
#' @param anno A GTF-like gene annotation GRanges object, which contains at least two columns: type and cluster.
#' Column type stores the type of element on the genome (gene, transcript, exon, CDS...), 
#' column cluster indicates the position on y axis.
#' @param region Genome region used to generate the transcript track plot, must be a GRanges object.
#' @param transcript_width Line width used to draw the transcript.
#' @param transcript_color Color used to draw the transcript.
#' @param exon_width Line width used to draw the exon.
#' @param exon_color Color used to draw the exon.
#' @param CDS_width Line width used to draw the CDS.
#' @param CDS_color Color used to draw the CDS.
#' @param arrow_break The gap between neighbor arrows equals to region width times arrow_break.
#' @param arrow_length Line length used to draw the arrow.
#' @param arrow_color Color used to draw the arrow.
#' @param show_name Which column in anno stores the transcript names to be plotted? Set to NULL and no transcript name will be shown.
#' @param name_size Text size used to draw the transcript name.
#' 
#' @return A ggplot object.
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
                                 show_name = NULL,
                                 name_size = 3){
  
  #check parameter
  if(base::class(anno) != 'GRanges'){
    base::stop('anno must be a GRanges object!')
  }else{
    anno <- rtracklayer::as.data.frame(anno,row.names = NULL)
  }
  
  if(!base::all(c('type','cluster') %in% base::colnames(anno))){
    base::stop('anno must contain 2 columns: type and cluster!')
  }else{
    anno$type <- base::as.character(anno$type)
    anno$cluster <- base::as.character(anno$cluster)
  }
  
  if(!base::is.null(show_name)){
    if(!(show_name %in% base::colnames(anno))){
      base::stop(base::paste(show_name,'column not found in anno!',sep = ' '))
    }else{
      anno$anno_name <- base::as.character(anno[,show_name])
    }
  }
  
  if(base::class(region) != 'GRanges'){
    base::stop('region must be a GRanges object!')
  }else{
    if(base::length(region) != 1){
      base::stop('only 1 region required!')
    }else{
      chr <- base::as.character(region@seqnames)
      start_site <- GenomicRanges::start(region)
      end_site <- GenomicRanges::end(region)
    }
  }
  
  if((arrow_break >= 1) | (arrow_break <= 0)){
    base::stop('The value of arrow_break must be between 0 and 1!')
  }
  
  #subset anno by overlap
  idx <- base::which((anno$seqnames == chr) & (anno$end >= start_site) & (anno$start <= end_site))
  
  if(base::length(idx) == 0){
    transcript_table <- base::data.frame(seqnames = chr,start = 0,end = 0,strand = '+',type = 'transcript',cluster = 1,anno_name = NA,row.names = NULL)
    transcript_table$mid_point <- base::round(x = (transcript_table$start + transcript_table$end)/2,digits = 0)
    exon_table <- base::data.frame(seqnames = chr,start = 0,end = 0,strand = '+',type = 'exon',cluster = 1,anno_name = NA,row.names = NULL)
    CDS_table <- base::data.frame(seqnames = chr,start = 0,end = 0,strand = '+',type = 'CDS',cluster = 1,anno_name = NA,row.names = NULL)
  }else{
    anno <- anno[idx,,drop = FALSE]
    
    #truncate range
    idx <- base::which(anno$start < start_site)
    if(base::length(idx) > 0){
      anno[idx,'start'] <- start_site
    }
    idx <- base::which(anno$end > end_site)
    if(base::length(idx) > 0){
      anno[idx,'end'] <- end_site
    }
    
    #transcript table
    idx <- base::which(anno$type == 'transcript')
    if(base::length(idx) == 0){
      transcript_table <- base::data.frame(seqnames = chr,start = 0,end = 0,strand = '+',type = 'transcript',cluster = 1,anno_name = NA,row.names = NULL)
    }else{
      transcript_table <- anno[idx,,drop = FALSE]
    }
    transcript_table$mid_point <- base::round(x = (transcript_table$start + transcript_table$end)/2,digits = 0)
    
    #exon table
    idx <- base::which(anno$type == 'exon')
    if(base::length(idx) == 0){
      exon_table <- base::data.frame(seqnames = chr,start = 0,end = 0,strand = '+',type = 'exon',cluster = 1,anno_name = NA,row.names = NULL)
    }else{
      exon_table <- anno[idx,,drop = FALSE]
    }
    
    #CDS table
    idx <- base::which(anno$type == 'CDS')
    if(base::length(idx) == 0){
      CDS_table <- base::data.frame(seqnames = chr,start = 0,end = 0,strand = '+',type = 'CDS',cluster = 1,anno_name = NA,row.names = NULL)
    }else{
      CDS_table <- anno[idx,,drop = FALSE]
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
        
        if(base::as.character(single_arrow$strand) == '-'){
          arrow_start_site <- single_arrow$end
          single_arrow$start <- arrow_start_site - (break_length * (y - 1))
          single_arrow$end <- arrow_start_site - (break_length * y)
        }else if(base::as.character(single_arrow$strand) == '+'){
          arrow_start_site <- single_arrow$start
          single_arrow$start <- arrow_start_site + (break_length * (y - 1))
          single_arrow$end <- arrow_start_site + (break_length * y)
        }else{
          single_arrow$start <- 0
          single_arrow$end <- 0
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
  
  #reduce arrow table
  arrow_table <- base::do.call(what = base::rbind,args = base::lapply(X = base::unique(arrow_table$cluster),FUN = function(x){
    break_length <- base::round(x = (region@ranges@width * arrow_break),digits = 0)
    idx <- base::which(arrow_table$cluster == x)
    
    if(base::length(idx) < 2){
      arrow_cluster <- arrow_table[idx,,drop = FALSE]
      
      return(arrow_cluster)
    }else{
      arrow_cluster <- arrow_table[idx,,drop = FALSE]
      
      keep_idx <- c(1)
      for(i in 2:(base::nrow(arrow_cluster))){
        diff_list <- base::abs(arrow_cluster[keep_idx,"end"] - arrow_cluster[i,"end"])
        if(base::sum(diff_list < break_length) == 0){
          keep_idx <- c(keep_idx,i)
        }
      }
      arrow_cluster <- arrow_cluster[keep_idx,,drop = FALSE]
      
      return(arrow_cluster)
    }
  }))
  
  #plot arrow
  gg_object <- ggplot2::ggplot() + 
    ggplot2::geom_segment(data = arrow_table,mapping = ggplot2::aes(x = start,xend = end,y = cluster,yend = cluster),linewidth = transcript_width,color = arrow_color,arrow = grid::arrow(length = grid::unit(x = arrow_length,units = 'inches'))) + 
    ggplot2::geom_segment(data = arrow_table,mapping = ggplot2::aes(x = start,xend = end,y = cluster,yend = cluster),linewidth = transcript_width,color = transcript_color)
  
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
  if(!base::is.null(show_name)){
    gg_object <- gg_object + 
      ggplot2::geom_text(data = transcript_table,mapping = ggplot2::aes(x = mid_point,y = cluster,label = anno_name),size = name_size,vjust = 1.5)
  }
  
  #modify ggplot object
  gg_object <- gg_object + 
    ggplot2::ylab('Transcript') + 
    ggplot2::xlab(base::paste0(chr,' : ',base::format(x = start_site,scientific = FALSE),' - ',base::format(x = end_site,scientific = FALSE))) + 
    ggplot2::scale_x_continuous(limits = c(start_site,end_site),expand = c(0,0)) + 
    ggplot2::theme_bw() + 
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank())
  
  #return
  return(gg_object)
}

#' Basic function for grouping transcripts into clusters
#' 
#' @description
#' Group transcripts into clusters for better visualization.
#' 
#' @param gene_anno A GRanges object which stores the gene annotation information.
#' @param column_name Which column in gene_anno stores the unique id for each transcript/gene?
#' 
#' @return Attaching cluster ids to gene_anno.
group_transcripts <- function(gene_anno,
                              column_name){
  
  #check parameter
  if(base::class(gene_anno) != 'GRanges'){
    base::stop('gene_anno must be a GRanges object!')
  }
  
  if(base::length(gene_anno) == 0){
    base::stop('no data in gene_anno!')
  }
  
  if(!base::all(c('type',column_name) %in% base::colnames(gene_anno@elementMetadata))){
    base::stop(base::paste0('gene_anno must contain two columns: type and ',column_name,'!'))
  }
  
  if(!base::all(!base::is.na(gene_anno@elementMetadata[,column_name]))){
    base::stop(base::paste0('NA is not allowed in gene_anno column: ',column_name,'!'))
  }
  
  #get transcript region
  transcript_list <- base::unique(base::as.character(gene_anno@elementMetadata[,column_name]))
  transcript_list <- base::do.call(what = base::c,args = base::lapply(X = transcript_list,FUN = function(x){
    idx <- base::which(gene_anno@elementMetadata[,column_name] == x)
    temp_anno <- gene_anno[idx]
    
    chr <- base::unique(base::as.character(temp_anno@seqnames))
    if(base::length(chr) != 1){
      base::stop('check the seqname for each transcript/gene!')
    }
    start_site <- base::as.character(base::min(GenomicRanges::start(temp_anno)))
    end_site <- base::as.character(base::max(GenomicRanges::end(temp_anno)))
    
    temp_anno <- methods::as(object = base::paste0(chr,':',start_site,'-',end_site),Class = 'GRanges')
    temp_anno$unique_name <- x
    
    return(temp_anno)
  }))
  
  #group transcript
  idx <- base::order(GenomicRanges::start(transcript_list),decreasing = FALSE)
  transcript_list <- transcript_list[idx]
  
  ordered_transcript_list <- transcript_list[1]
  ordered_transcript_list$cluster <- 1
  transcript_list <- transcript_list[-1]
  
  cluster_list <- c(1)
  cluster_end <- c(GenomicRanges::end(ordered_transcript_list))
  
  while(base::length(transcript_list) > 0){
    temp_transcript <- transcript_list[1]
    transcript_list <- transcript_list[-1]
    
    idx <- base::which(GenomicRanges::start(temp_transcript) > cluster_end)
    if(base::length(idx) == 0){
      temp_transcript$cluster <- base::max(cluster_list) + 1
      cluster_list <- c(cluster_list,temp_transcript$cluster)
      cluster_end <- c(cluster_end,GenomicRanges::end(temp_transcript))
    }else{
      idx <- base::min(idx)
      temp_transcript$cluster <- cluster_list[idx]
      cluster_end[idx] <- GenomicRanges::end(temp_transcript)
    }
    
    ordered_transcript_list <- c(ordered_transcript_list,temp_transcript)
  }
  
  #group gene_anno
  cluster_list <- ordered_transcript_list$cluster
  base::names(cluster_list) <- ordered_transcript_list$unique_name
  gene_anno$cluster <- cluster_list[base::as.character(gene_anno@elementMetadata[,column_name])]
  
  #return
  return(gene_anno)
}
