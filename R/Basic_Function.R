#' Truncate GRanges object by specified genome region
#' 
#' @description
#' All segments beyond the specified genome region will be truncated.
#' 
#' @param single_Range A single GRanges object to be truncated.
#' @param chr Chromosome name of the specified genome region.
#' @param start_site Start position of the specified genome region.
#' @param end_site End position of the specified genome region.
#' @param must_return Regulate the returned object in the scenario where the GRanges object doesn't intersect with the specified genome region. If set to TRUE, a GRanges object with a 0-0 range would be returned. If set to FALSE, this function will directly return NULL.
#' 
#' @return A truncated GRanges object.
truncate_GRanges <- function(single_Range,
                             chr,
                             start_site,
                             end_site,
                             must_return = TRUE){
  
  #check parameter
  if(base::class(single_Range) != 'GRanges'){
    base::stop('single_Range must be a single GRanges object!')
  }
  if(base::length(single_Range) == 0){
    base::stop('single_Range is empty!')
  }
  
  #subset by overlap
  region <- methods::as(object = base::paste0(chr,':',start_site,'-',end_site),Class = 'GRanges')
  subset_Range <- IRanges::subsetByOverlaps(x = single_Range,ranges = region,type = 'any',invert = FALSE,ignore.strand = TRUE)
  
  #return
  if(base::length(subset_Range) == 0){
    if(must_return == TRUE){
      temp <- S4Vectors::mcols(single_Range)[1,,drop = FALSE]
      subset_Range <- methods::as(object = base::paste0(chr,':','0-0'),'GRanges')
      S4Vectors::mcols(subset_Range) <- temp
      return(subset_Range)
    }else{
      return(NULL)
    }
  }else{
    #truncate start site
    each_start_site <- GenomicRanges::start(subset_Range)
    idx <- base::which(each_start_site < start_site)
    if(base::length(idx) > 0){
      each_start_site[idx] <- start_site
    }
    
    #truncate end site
    each_end_site <- GenomicRanges::end(subset_Range)
    idx <- base::which(each_end_site > end_site)
    if(base::length(idx) > 0){
      each_end_site[idx] <- end_site
    }
    
    GenomicRanges::start(subset_Range) <- each_start_site
    GenomicRanges::end(subset_Range) <- each_end_site
    return(subset_Range)
  }
}

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
      chr <- base::as.character(GenomicRanges::seqnames(region))
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
      chr <- base::as.character(GenomicRanges::seqnames(region))
      start_site <- GenomicRanges::start(region)
      end_site <- GenomicRanges::end(region)
    }
  }
  
  #convert to truncated GRanges list
  if(!base::grepl(pattern = 'GRanges',x = base::class(Ranges),fixed = TRUE)){
    base::stop('Ranges must be a GRanges or GRangesList object!')
  }else{
    if(base::grepl(pattern = 'List$',x = base::class(Ranges),fixed = FALSE)){
      Ranges_name <- base::names(Ranges)
      Ranges <- base::do.call(what = S4Vectors::SimpleList,args = base::lapply(X = c(1:(base::length(Ranges))),FUN = function(x){
        temp_Range <- Ranges[[x]]
        temp_Range <- truncate_GRanges(single_Range = temp_Range,chr = chr,start_site = start_site,end_site = end_site,must_return = TRUE)
      }))
      names(Ranges) <- Ranges_name
    }else{
      Ranges <- truncate_GRanges(single_Range = Ranges,chr = chr,start_site = start_site,end_site = end_site,must_return = TRUE)
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
  }
  
  if(!base::all(c('type','cluster') %in% base::colnames(S4Vectors::mcols(anno)))){
    base::stop('anno must contain 2 columns: type and cluster!')
  }else{
    anno$type <- base::as.character(anno$type)
    anno$cluster <- base::as.character(anno$cluster)
  }
  
  if(!base::is.null(show_name)){
    if(!(show_name %in% base::colnames(S4Vectors::mcols(anno)))){
      base::stop(base::paste(show_name,'column not found in anno!',sep = ' '))
    }else{
      anno$anno_name <- base::as.character(S4Vectors::mcols(anno)[,show_name])
    }
  }
  
  if(base::class(region) != 'GRanges'){
    base::stop('region must be a GRanges object!')
  }else{
    if(base::length(region) != 1){
      base::stop('only 1 region required!')
    }else{
      chr <- base::as.character(GenomicRanges::seqnames(region))
      start_site <- GenomicRanges::start(region)
      end_site <- GenomicRanges::end(region)
    }
  }
  
  if((arrow_break >= 1) | (arrow_break <= 0)){
    base::stop('The value of arrow_break must be between 0 and 1!')
  }
  
  #subset anno by overlap
  anno <- truncate_GRanges(single_Range = anno,chr = chr,start_site = start_site,end_site = end_site,must_return = FALSE)
  
  if(base::length(anno) == 0){
    transcript_table <- base::data.frame(seqnames = chr,start = 0,end = 0,strand = '+',type = 'transcript',cluster = 1,anno_name = NA,row.names = NULL)
    transcript_table$mid_point <- base::round(x = (transcript_table$start + transcript_table$end)/2,digits = 0)
    exon_table <- base::data.frame(seqnames = chr,start = 0,end = 0,strand = '+',type = 'exon',cluster = 1,anno_name = NA,row.names = NULL)
    CDS_table <- base::data.frame(seqnames = chr,start = 0,end = 0,strand = '+',type = 'CDS',cluster = 1,anno_name = NA,row.names = NULL)
  }else{
    anno <- rtracklayer::as.data.frame(anno,row.names = NULL)
    
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
          transcript_start_site <- single_transcript$end
          single_arrow$end <- transcript_start_site - (break_length * y)
          single_arrow$start <- single_arrow$end + 1
        }else if(base::as.character(single_arrow$strand) == '+'){
          transcript_start_site <- single_transcript$start
          single_arrow$end <- transcript_start_site + (break_length * y)
          single_arrow$start <- single_arrow$end - 1
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
      for (i in 2:(base::nrow(arrow_cluster))) {
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
  
  if(!base::all(c('type',column_name) %in% base::colnames(S4Vectors::mcols(gene_anno)))){
    base::stop(base::paste0('gene_anno must contain two columns: type and ',column_name,'!'))
  }
  
  if(!base::all(!base::is.na(S4Vectors::mcols(gene_anno)[,column_name]))){
    base::stop(base::paste0('NA is not allowed in gene_anno column: ',column_name,'!'))
  }
  if(!base::all(!base::is.null(S4Vectors::mcols(gene_anno)[,column_name]))){
    base::stop(base::paste0('NULL is not allowed in gene_anno column: ',column_name,'!'))
  }
  
  #get transcript region
  transcript_list <- base::unique(base::as.character(S4Vectors::mcols(gene_anno)[,column_name]))
  transcript_list <- base::do.call(what = base::c,args = base::lapply(X = transcript_list,FUN = function(x){
    idx <- base::which(S4Vectors::mcols(gene_anno)[,column_name] == x)
    temp_anno <- gene_anno[idx]
    
    chr <- base::unique(base::as.character(GenomicRanges::seqnames(temp_anno)))
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
  cluster_end <- c(GenomicRanges::end(ordered_transcript_list))
  
  while(base::length(transcript_list) > 0){
    temp_transcript <- transcript_list[1]
    transcript_list <- transcript_list[-1]
    
    idx <- base::which(GenomicRanges::start(temp_transcript) > cluster_end)
    if(base::length(idx) == 0){
      temp_transcript$cluster <- base::length(cluster_end) + 1
      cluster_end <- c(cluster_end,GenomicRanges::end(temp_transcript))
    }else{
      idx <- base::min(idx)
      temp_transcript$cluster <- idx
      cluster_end[idx] <- GenomicRanges::end(temp_transcript)
    }
    
    ordered_transcript_list <- c(ordered_transcript_list,temp_transcript)
  }
  
  #group gene_anno
  cluster_list <- ordered_transcript_list$cluster
  base::names(cluster_list) <- ordered_transcript_list$unique_name
  gene_anno$cluster <- cluster_list[base::as.character(S4Vectors::mcols(gene_anno)[,column_name])]
  
  #return
  return(gene_anno)
}

#' Basic function for genome position linkage visualization
#' 
#' @description
#' Generate directional curves to characterize the linkage between cis-elements on genome, like the activating effect of enhancers to promoters, and so on.
#' 
#' @param linkage A GRanges object with each range's start and end points representing two loci on the genome where a linkage exists, and the strand represents the direction of the linkage.
#' @param region Genome region used to generate the linkage track plot, must be a GRanges object.
#' @param color_by Which column in the linkage stores the intensity of each linkage that is used to control the color? Set to NULL and each linkage will be the same color.
#' @param col_pal A custom palette used to override coloring for each linkage.
#' @param allow_truncated Whether to display linkages that are truncated due to exceeding the region.
#' @param curve_width Line width used to draw the linkage.
#' @param max_arrow_length Max line length used to draw the arrow.
#' 
#' @return A ggplot object.
linkage_vis_basic <- function(linkage,
                              region,
                              color_by = NULL,
                              col_pal = c('#E6E7E8','#3A97FF','#8816A7','#000000'),
                              allow_truncated = FALSE,
                              curve_width = 0.5,
                              max_arrow_length = 0.08){
  
  #check parameter
  if(base::class(linkage) != 'GRanges'){
    base::stop('linkage must be a GRanges object!')
  }
  
  if(!base::is.null(color_by)){
    if(!(color_by %in% base::colnames(S4Vectors::mcols(linkage)))){
      base::stop(base::paste0('linkage must contain the column: ',color_by,'!'))
    }else{
      linkage$value <- base::as.numeric(S4Vectors::mcols(linkage)[,color_by])
    }
  }else{
    linkage$value <- 1
  }
  
  if(base::class(region) != 'GRanges'){
    base::stop('region must be a GRanges object!')
  }else{
    if(base::length(region) != 1){
      base::stop('only 1 region required!')
    }else{
      chr <- base::as.character(GenomicRanges::seqnames(region))
      start_site <- GenomicRanges::start(region)
      end_site <- GenomicRanges::end(region)
    }
  }
  
  #subset linkage by region
  if(allow_truncated){
    linkage <- IRanges::subsetByOverlaps(x = linkage,ranges = region,type = 'any',invert = FALSE,ignore.strand = TRUE)
  }else{
    linkage <- IRanges::subsetByOverlaps(x = linkage,ranges = region,type = 'within',invert = FALSE,ignore.strand = TRUE)
  }
  if(base::length(linkage) == 0){
    linkage <- methods::as(object = base::paste0(chr,':','0-0'),Class = 'GRanges')
    linkage$value <- 1
  }
  
  #break the linkage curve into small segments
  angles <- base::seq(from = base::pi,to = 2*(base::pi),length.out = 100)
  radius_list <- base::unlist(base::lapply(X = 1:(base::length(linkage)),FUN = function(x){
    temp_linkage <- linkage[x]
    if(base::as.character(GenomicRanges::strand(temp_linkage)) == '-'){
      temp_radius <- (GenomicRanges::start(temp_linkage) - GenomicRanges::end(temp_linkage))/2
    }else{
      temp_radius <- (GenomicRanges::end(temp_linkage) - GenomicRanges::start(temp_linkage))/2
    }
    return(temp_radius)
  }))
  if(base::max(base::abs(radius_list)) == 0){
    max_radius <- 1
  }else{
    max_radius <- base::max(base::abs(radius_list))
  }
  circle_center_list <- c(GenomicRanges::start(linkage) + base::abs(radius_list))
  
  breaked_segments <- base::do.call(what = base::rbind,args = base::lapply(X = 1:(base::length(linkage)),FUN = function(x){
    segment_x_list <- (radius_list[x])*base::cos(angles) + circle_center_list[x]
    segment_y_list <- (base::abs(radius_list[x])/max_radius)*base::sin(angles) + 0
    segments_table <- base::data.frame(x = segment_x_list,
                                       y = segment_y_list,
                                       group = base::paste0('l',x),
                                       radius_ratio = base::abs(radius_list[x])/max_radius,
                                       value = linkage$value[x],
                                       row.names = NULL)
    return(segments_table)
  }))
  
  #truncate linkage segments
  idx_to_remove <- base::which((breaked_segments$x < start_site) | (breaked_segments$x > end_site))
  if(base::length(idx_to_remove) > 0){
    breaked_segments <- breaked_segments[-c(idx_to_remove),,drop = FALSE]
  }
  
  #basic plot
  gg_object <- ggplot2::ggplot(data = breaked_segments,mapping = ggplot2::aes(x = x,y = y,group = group,color = value))
  for (i in base::unique(breaked_segments$group)) {
    temp_breaked_segments <- breaked_segments[base::which(breaked_segments$group == i),,drop = FALSE]
    gg_object <- gg_object + 
      ggplot2::geom_path(data = temp_breaked_segments,
                         mapping = ggplot2::aes(x = x,y = y,group = group,color = value),
                         arrow = grid::arrow(angle = 15,length = grid::unit(x = max_arrow_length*base::unique(base::as.numeric(temp_breaked_segments$radius_ratio)),units = 'inches'),ends = 'last'),
                         linewidth = curve_width)
  }
  
  gg_object <- gg_object + 
    ggplot2::ylab('Linkage') + 
    ggplot2::xlab(base::paste0(chr,' : ',base::format(x = start_site,scientific = FALSE),' - ',base::format(x = end_site,scientific = FALSE))) + 
    ggplot2::scale_x_continuous(limits = c(start_site,end_site),expand = c(0,0)) + 
    ggplot2::theme_bw() + 
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank())
  
  #color
  if(!base::is.null(col_pal)){
    gg_object <- gg_object + 
      ggplot2::scale_color_gradientn(colours = base::as.character(col_pal))
  }
  
  #return
  return(gg_object)
}

#' Convert a linkage data frame to GRanges object that can be directly used for linkage visualization
#' 
#' @description
#' The GRanges object requires the start position to be less than the end position, which is a bit cumbersome for characterizing the direction of linkage.
#' Here provides a function that helps to convert a data frame that stores information about pairwise interactions of different loci on the genome into a GRanges object, which is convenient for linkage visualization.
#' 
#' @param linkage A data frame that must contain the following three columns: 
#' 1. The "chrom" column indicates the chromosome where the linkage occurs.
#' 2. The "start" column represents the starting point of the linkage.
#' 3. The "end" column denotes the ending point of the linkage.
#' 
#' @return A GRanges object.
#' 
#' @export
linkage_table_to_GRanges <- function(linkage){
  
  #check parameter
  linkage <- base::as.data.frame(x = linkage,row.names = NULL)
  if(!base::all(c('chrom','start','end') %in% base::colnames(linkage))){
    base::stop('linkage must contain three columns: chrom, start and end!')
  }
  if(base::nrow(linkage) < 1){
    base::stop('No linkage detected!')
  }
  
  #create GRanges object
  linkage_GRanges <- base::do.call(what = base::c,args = base::lapply(X = 1:(base::nrow(linkage)),FUN = function(x){
    chr <- base::as.character(linkage[x,'chrom'])
    if(base::as.numeric(linkage[x,'end']) >= base::as.numeric(linkage[x,'start'])){
      start_site <- base::as.numeric(linkage[x,'start'])
      end_site <- base::as.numeric(linkage[x,'end'])
      strand_sign <- '+'
    }else{
      start_site <- base::as.numeric(linkage[x,'end'])
      end_site <- base::as.numeric(linkage[x,'start'])
      strand_sign <- '-'
    }
    temp <- methods::as(object = base::paste0(chr,':',start_site,'-',end_site,':',strand_sign),Class = 'GRanges')
    return(temp)
  }))
  
  #add elementMetadata
  idx_to_remove <- base::which(base::colnames(linkage) %in% c('chrom','start','end'))
  meta_data <- linkage[,-c(idx_to_remove),drop = FALSE]
  S4Vectors::mcols(linkage_GRanges) <- meta_data
  
  #return
  return(linkage_GRanges)
}

#' Highlight specified genome region on the generated track plots
#' 
#' @description
#' Highlight specified genome region on the generated track plots (for example, coverage track, feature track and transcript track).
#' 
#' @param gg_object Track plot as a ggplot object.
#' @param start Start position for the genome region to be highlighted.
#' @param end End position for the genome region to be highlighted.
#' @param color Highlight color for the specified genome region.
#' 
#' @return A ggplot object.
#' 
#' @export
highlight_region <- function(gg_object,
                             start,
                             end,
                             color = '#FFEDAC'){
  
  #check parameter
  if(!base::all(base::class(gg_object) %in% c('gg','ggplot'))){
    base::stop('gg_object must be a ggplot object!')
  }
  
  if(!((base::length(start) == 1) & (base::length(end) == 1))){
    base::stop('only 1 region required!')
  }
  if(base::length(color) != 1){
    base::stop('only 1 color required!')
  }
  
  start <- base::as.numeric(start)
  end <- base::as.numeric(end)
  color <- base::as.character(color)
  
  #add highlight region layer
  highlight_layer <- ggplot2::geom_rect(xmin = start,xmax = end,ymin = -Inf,ymax = Inf,color = NA,fill = color,stat = 'identity',position = 'identity',show.legend = FALSE)
  gg_object$layers <- base::append(x = gg_object$layers,values = highlight_layer,after = 0)
  
  #return
  return(gg_object)
}
