## for handling R CMD check on non standard evaluation calls
utils::globalVariables(c("reference", "size","strand",
                         "start5p","end3p","sizeGenome",
                         "firstBase","acf","lag", "n"))

##########################################################################################################################aux function
#' multiplot
#'
#' plot multiple something
#' function created by Winston Chang,
#' in cookbook for R
#' http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
#'
#' @param ... xxx
#'
#' @param plotlist xxx
#'
#' @param file xxx
#'
#' @param cols xxx
#'
#' @param layout xxx
#'
#' @return The sum of \code{x} and \code{y}.
#'
#' @export
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  ##library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = grid::viewport(layout.pos.row = matchidx$row,layout.pos.col = matchidx$col))
    }
  }
}

#######################################################################################################################loading sam file
#' Read and load SAM file
#'
#' read and load a sam file if the input is a valid path
#'
#' @param f SAM file path.
#'
#' @param sizeStart Minimum value of size to evaluate
#'
#' @param sizeEnd Maximum value of size to evaluate
#'
#' @export
loadSam <- function(f, sizeStart=11, sizeEnd=1000){

  if(!file.exists(f)){
    print(paste("File ", f, " does not exists!", sep=""))
    return(FALSE)
  }else{
    print(paste("Loading ", f, " file...", sep=""))
    aux <- utils::read.delim(f, sep="\t", header=FALSE, comment.char = "@")

    samfile = aux[, c(1, 2, 3, 4, 10)] # getting important fields [read name, strand, reference name, start mapping, sequence]
    remove(aux) #remove temp object

    colnames(samfile) <- c("read","strand","reference","start5p","sequence")
    samfile$sequence <- as.character(samfile$sequence)
    samfile$size <- nchar(samfile$sequence)
    samfile$end3p <- ifelse(samfile$strand == 0, samfile$start5p + samfile$size - 1, samfile$start5p)
    samfile$start5p <- ifelse(samfile$strand == 16, samfile$end3p + samfile$size - 1, samfile$start5p)

    samfile <- subset(samfile, size >= sizeStart & size <= sizeEnd) #getting reads in interval sizeStart>= reads <= sizeEnd
    samfile$firstBase <- substr(samfile$sequence, 1, 1)


    return(samfile)
  }
}

#######################################################################################################################loading sam file
#' Read a SAM file and return it in BED format
#'
#' read a SAM file and convert the structure to a BED file
#'
#' @param f SAM file path.
#'
#' @param sizeStart Minimum value of size to evaluate
#'
#' @param sizeEnd Maximum value of size to evaluate
#'
#' @export
loadSamLikeBed <- function(f,sizeStart=11,sizeEnd=1000){

  if(!file.exists(f)){
    print(paste("File ",f," does not exists!",sep=""))
    return(FALSE)
  }else{
    print(paste("Loading ",f," file...",sep=""))
    aux = utils::read.delim(f,sep="\t",header=FALSE,comment.char = "@")

    samfile = aux[,c(1,2,3,4,10)] # getting important fields [read name, strand, reference name, start mapping, sequence]
    remove(aux) #remove temp file

    colnames(samfile)<- c("read","strand","reference","start5p","sequence")
    samfile$sequence=as.character(samfile$sequence)
    samfile$size = nchar(samfile$sequence)
    samfile$end3p = ifelse(samfile$strand == 0, samfile$start5p + samfile$size -1, samfile$start5p)
    samfile$start5p = ifelse(samfile$strand == 16, samfile$end3p + samfile$size -1, samfile$start5p -1)

    samfile = subset(samfile, "size" >= sizeStart & "size" <= sizeEnd ) #getting reads in interval sizeStart>= reads <= sizeEnd
    samfile$firstBase = substr(samfile$sequence,1,1)


    return(samfile)
  }
}

################################################################################################################plot distribution per base
#' Size distribution plots
#'
#' This function returns a ggplot object with ....
#'
#' @param samObject XXXX
#
#' @param ref xxx
#'
#' @param stranded xxx

#' @param norm xxx
#'
#' @param sizeStart xxx
#'
#' @param sizeEnd xxx
#'
#' @param ymax xxx
#'
#' @param ymin xxx
#'
#' @export
plotSizeDistribution <- function(samObject, ref=NULL, stranded=TRUE, norm=NULL,
                                 sizeStart=NULL, sizeEnd=NULL, ymax=NULL, ymin=NULL){

  ylegend = "Number of small RNAs"
  cbPalette <- c( "#D55E00","#F0E442","#0072B2","#009E73","#999999", "#E69F00", "#56B4E9", "#CC79A7")

  if(!is.null(ref)){
    samObject = subset(samObject, reference == ref)
  }

  g = dplyr::group_by(samObject, reference, strand, firstBase, size)
  m = dplyr::count(g, reference, strand, size, firstBase)

  #defining size interval to plot
  if(!is.null(sizeStart) && !is.null(sizeEnd)){
    xmin = sizeStart -1
    xmax=sizeEnd + 1
  }else{
    xmin = min(m$size)
    xmax = max(m$size)
  }

  #normalizing
  if(!is.null(norm) ){
    m$n = (m$n / norm ) * 1000000
    ylegend = "Number of small RNAs (RPM)"
  }

  ##defining limits
  if(!is.null(ymax) && !is.null(ymin)){

    limit1_min = ymin
    limit1_max = ymax

  }else{

    m2 = dplyr::group_by(m, size, strand)
    s1 = dplyr::summarise(m2, n = sum(n))
    limit1_min = - max(s1$n)
    limit1_max = max(s1$n)

    m3=dplyr::group_by(m, size)
    s2 = dplyr::summarise(m3, n = sum(n))
    limit2 = max(s2$n)
  }

  ##multiple value by -1 for negative strand
  m[m$strand==16,"n"] = m[m$strand==16,"n"] * -1
  utils::head(m)


  if(stranded==TRUE){
  ###both strands
  plot <- ggplot2::ggplot() +
    ggplot2::theme_classic(base_size = 16) +
    ggplot2::geom_bar(data = m, ggplot2::aes(size, n, fill = firstBase), stat = "identity", width=.8) +
    ggplot2::scale_fill_manual("5' base\npreference", values = cbPalette) +
    ggplot2::ylim(limit1_min,limit1_max) +
    ggplot2::xlab("small RNA length (nt)") +
    ggplot2::ylab(ylegend) +
    ggplot2::theme(panel.border = ggplot2::element_blank()) +
    ggplot2::xlim(xmin,xmax) +
    ggplot2::theme(axis.line = ggplot2::element_line(color = 'black'))  +
    ggplot2::facet_grid(~reference)
  }else{
  #sum strands
  plot <- ggplot2::ggplot() +
    ggplot2::theme_classic( base_size = 16) +
    ggplot2::geom_bar( data = m, ggplot2::aes(size,abs(n),fill=firstBase),stat="identity",width=.8) +
    ggplot2::scale_fill_manual("5' base\npreference",values = cbPalette) +
    ggplot2::ylim(0,limit2) +
    ggplot2::xlab("small RNA length (nt)") +
    ggplot2::ylab(ylegend) +
    ggplot2::theme(panel.border = ggplot2::element_blank()) +
    ggplot2::xlim(xmin,xmax) +
    ggplot2::theme(axis.line = ggplot2::element_line(color = 'black'))   +
    ggplot2::facet_grid(~reference)
  }

  return(plot)
}

###############################################################################################################calculating Density Per Base
#' Calculate
#'
#' calculating Density Per Base
#'
#' @param samObject XXXX
#
#' @param ref xxx
#'
#' @param readSize xxx

#' @param bin xxx
#'
#' @param refSize xxx
#'
#' @param norm xxx
#'
#' @param xmin xxx
#'
#' @param xmax xxx
#'
#' @param ymax xxx
#'
#' @param ymin xxx
#'
#' @export
calcDensityPerBase <- function(samObject, ref, readSize=0, bin=1,
                               refSize=NULL, norm=NULL, xmin=NULL,
                               xmax=NULL, ymax=NULL, ymin=NULL){


  ##defining reads to use
  if(readSize==0){
    samObject = subset( samObject, reference == ref )
  }else{
    samObject = subset(samObject, reference == ref & size %in% readSize)
  }


  if(is.null(refSize)){
    sizeGenome = max(samObject$start5p)
  }else{
    sizeGenome = refSize
  }


  #calc density per base
  cov_positive = integer(sizeGenome)
  cov_negative = integer(sizeGenome)

  for (i in 1:nrow(samObject) ){
    start = samObject[i,"start5p"]
    end = samObject[i,"end3p"]
    strand = samObject[i,"strand"]

    center = (end - ((end - start)/2))
    round = as.integer(center)
    if (round > sizeGenome) {
      print("Error: Read mapped outside frontiers of genome...")
    }

    if (strand == 0 ) {
      for (e in start:end) {
        cov_positive[e]=cov_positive[e]+1
      }
    }
    if (strand == 16 ) {
      for (e in start:end) {
        cov_negative[e]=cov_negative[e]+1
      }
    }
  }

  #calc density per base [BINS]
  bin=100
  reference="VSV"
  df_cov=data.frame()
  for(i in seq(1,sizeGenome,by=bin)){

    count=0
    totalPos=0
    totalNeg=0
    sum=0
    start=i

    while(count<=bin){
      totalPos = totalPos + cov_positive[(i + count)]
      totalNeg = totalNeg + cov_negative[(i + count)]
      sum = sum + totalPos + totalNeg;
      count=count+1
    }
    end = (i + count - 1)

    df_cov=rbind(df_cov,data.frame(reference,start,end,totalPos,totalNeg,sum,stringsAsFactors = FALSE))

  }

  df_cov[is.na(df_cov)] = 0 # replace NA -> 0

  #normalizing
  if(!is.null(norm) ){
    df_cov$totalPos = (df_cov$totalPos/ norm ) * 1000000
    df_cov$totalNeg = (df_cov$totalNeg/ norm ) * 1000000
  }

  ##defining limits
  if(!is.null(ymax) && !is.null(ymin)){

    limit_ymin = ymin
    limit_ymax = ymax

  }else{

    limit_ymin = -max(df_cov$totalNeg,df_cov$totalPos)
    limit_ymax = max(df_cov$totalNeg,df_cov$totalPos)
  }

  #defining X and Y limits
  if(!is.null(xmin) && !is.null(xmax)){
    limit_xmin = xmin
    limit_xmax=  xmax
  }else{
    limit_xmin = 1
    limit_xmax= sizeGenome
  }


  ggplot2::ggplot(data=df_cov, ggplot2::aes(x=start)) +
    ggplot2::ylab("Small RNA density") +
    ggplot2::xlab("Reference sequence (nt)") +
    ggplot2::ylim(limit_ymin,limit_ymax) +
    ggplot2::theme_classic(base_size=16) +
    ggplot2::xlim(limit_xmin,limit_xmax) +
    ggplot2::geom_density(ggplot2::aes(y=totalPos),fill="#56B4E9",colour="#56B4E9",size=0.1,stat="identity") +
    ggplot2::geom_density(ggplot2::aes(y=totalNeg*-1),colour="#D55E00",fill="#D55E00",size=0.1,stat="identity")  +
    ggplot2::theme(panel.border = ggplot2::element_blank()) +
    ggplot2::theme(axis.line = ggplot2::element_line(color = 'black'))

}

#######################################################################################################################
#' weblogo function
#'
#' plot web logo
#'
#' @param samObject XXXX
#
#' @param ref xxx
#'
#' @param sizeStart xxx
#'
#' @param sizeEnd xxx
#'
#' @export
createWebLogo <- function(samObject, ref, sizeStart=18, sizeEnd=30){

  #selecting sequence target
  target = subset(samObject,reference==ref & size >=sizeStart & sizeEnd <=sizeEnd)
  #trimming sequence to the same size (minimum size = sizeStart)
  target$sequence = substr(target$sequence, 1, sizeStart)

  #counting number of reads in each strand
  nreadsSense = nrow(target[target$strand ==0,])
  nreadsAntisense = nrow(target[target$strand ==16,])

  #make sure sequences are strings, upper case and replace T to U
  target$sequence = as.character(toupper(target$sequence))
  target$sequence <- gsub("T", "U", target$sequence)

  #weblogo (+) reads
  p1 <- ggplot2::ggplot(data = target) +
    ggseqlogo::geom_logo(target[target$strand==0,"sequence"], method = 'bits',seq_type = 'rna') +
    ggseqlogo::theme_logo() +
    ggplot2::theme_classic() +
    ggplot2::ylim(0,1)+
    ggplot2::xlab(paste("nucleotide position - #sequences (+) : ",nreadsSense,sep=""))

  #weblogo (-) reads
  p2 <- ggplot2::ggplot(data = target) +
    ggseqlogo::geom_logo(target[target$strand==16,"sequence"], method = 'bits',seq_type = 'rna') +
    ggseqlogo::theme_logo() +
    ggplot2::theme_classic() +
    ggplot2::ylim(0,1) +
    ggplot2::xlab(paste("nucleotide position - #sequences (-) : ",nreadsAntisense,sep=""))

  #multiplot (+) and (-) weblogos
  return (multiplot(p1, p2, cols=1))

}

#####################################################################################################################
#' ACF
#'
#' calculate ACF
#'
#' @param samObject XXXX
#
#' @param ref xxx
#'
#' @param sizeStart xxx
#'
#' @param sizeEnd xxx
#'
#' @param regIni xxx
#'
#' @param regEnd xxx
#'
#' @export
calcACF <- function(samObject, ref, sizeStart=15, sizeEnd=30, regIni=1, regEnd=100){

  #selecting sequence target
  target_pos = dplyr::count(subset(samObject,reference==ref & size >=sizeStart & sizeEnd <=sizeEnd & strand==0),reference,start5p)
  target_neg = dplyr::count(subset(samObject,reference==ref & size >=sizeStart & sizeEnd <=sizeEnd & strand==16),reference,start5p)

  pos=integer()
  neg=integer()
  total=integer()

  for (i in regIni:regEnd){

    p=ifelse(is.element(i,target_pos$start5p),as.numeric(as.character(target_pos[is.element(target_pos$start5p,i), "n"])),0)
    n=ifelse(is.element(i,target_neg$start5p),as.numeric(as.character(target_neg[is.element(target_neg$start5p,i), "n"])),0)

    pos = append(pos,p)
    neg = append(n,neg)
    total = append(total,p+n)
  }


  posacf <- stats::acf(pos, plot = FALSE)
  posacf2 <- with(posacf, data.frame(lag, acf))

  acfPOS <- ggplot2::ggplot(data = posacf2, mapping = ggplot2::aes(x = lag, y = acf)) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = 0)) +
    ggplot2::geom_segment(mapping = ggplot2::aes(xend = lag, yend = 0)) +
    ggplot2::theme_bw(base_size=12) +
    ggplot2::xlab(paste("lag (reads sense: ",sum(target_pos$n),")",sep=""))

  negacf <- stats::acf(neg, plot = FALSE)
  negacf2 <- with(negacf, data.frame(lag, acf))

  acfNEG <- ggplot2::ggplot(data = negacf2, mapping = ggplot2::aes(x = lag, y = acf)) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = 0)) +
    ggplot2::geom_segment(mapping = ggplot2::aes(xend = lag, yend = 0)) +
    ggplot2::theme_bw(base_size=12) +
    ggplot2::xlab(paste("lag (reads antisense: ",sum(target_neg$n),")",sep=""))

  totacf <- stats::acf(total, plot = FALSE)
  totacf2 <- with(totacf, data.frame(lag, acf))

  acfTOTAL <- ggplot2::ggplot(data = totacf2, mapping = ggplot2::aes(x = lag, y = acf)) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = 0)) +
    ggplot2::geom_segment(mapping = ggplot2::aes(xend = lag, yend = 0)) +
    ggplot2::theme_bw(base_size=12) +
    ggplot2::xlab(paste("lag (all reads : ",sum(target_pos$n,target_neg$n),")",sep=""))

  multiplot(acfPOS,acfNEG,acfTOTAL,cols = 1)
}


#####################################################################################################################
#' Calculate phasing value
#'
#' calculate phasing
#'
#' @param samObject XXXX
#
#' @param ref xxx
#'
#' @param extend_upstream xxx
#'
#' @param extend_downstream xxx
#'
#' @param readSize xxx
#'
#' @param bin xxx
#'
#' @param refSize xxx
#'
#' @param norm xxx
#'
#' @param region_start xxx
#'
#' @param region_end xxx
#'
#' @param ymax xxx
#'
#' @param ymin xxx
#'
#' @export
calcPhasing <- function(samObject, ref, extend_upstream=20,
                        extend_downstream=20, readSize=0, bin=1,
                        refSize=NULL, norm=NULL, region_start=NULL,
                        region_end=NULL, ymax=NULL, ymin=NULL){

  #1 : Parse command line
  #2: set_filter_region
  #3: load_external_intervals

  ##defining reads to use
  if(!is.null(min) && !is.null(max)){
    region_start = min
    region_end=  max
  }else{
    region_start = 1
    region_end= sizeGenome
  }

  #set_filter_region
  if(readSize==0){
    samObject = subset(samObject,reference==ref & (start5p %in% region_start:region_end))
  }else{
    samObject = subset(samObject,reference==ref & (size %in% readSize) & (start5p %in% region_start:region_end) )
  }

  region_size = region_end - region_start + 1

  #initialize_empty_windows
  window_positive_sequences = integer(region_size)
  window_negative_sequences = integer(region_size)
  window_positive_reads     = integer(region_size)
  window_negative_reads     = integer(region_size)


  #calculate_region_coverage
  collapsed = dplyr::count(samObject,reference,strand,size,start5p,end3p)
  collapsed = as.data.frame(collapsed)
  #utils::head(collapsed)

  for( i in 1:nrow(collapsed)){
    chrom  = collapsed[i,"reference"]
    start  = collapsed[i,"start5p"]
    end    = collapsed[i,"end3p"] -1
#    name   = collapsed[i,"start5p"]
    strand = collapsed[i,"strand"]
    multiplicity=collapsed[i,"n"]

    if ( strand == 0 ) {
      iindex = as.numeric(start - region_start) + 1 # add 1 to index because 0 possibility

      window_positive_sequences[iindex]= window_positive_sequences[iindex] + 1
      window_positive_reads[iindex]    = window_positive_reads[iindex] + multiplicity
      } else {
        iindex = as.numeric(end - region_start) + 1 # add 1 to index because 0 possibility

        window_negative_sequences[iindex] = window_negative_sequences[iindex] + 1
        window_negative_reads[iindex]     = window_negative_reads[iindex] + multiplicity
    }
  }

  utils::head(window_positive_sequences)
  utils::head(window_positive_reads)

  ####sort_intervals_to_process
  collapsed = collapsed[order(collapsed$strand,collapsed$start5p,collapsed$end3p),]

  result = data.frame()

  #scan_region_windows()
  for( i in 1:nrow(collapsed)){
    chrom  = collapsed[i,"reference"]
    start  = collapsed[i,"start5p"]
    end    = collapsed[i,"end3p"]
    #    name   = collapsed[i,"start5p"]
    strand = collapsed[i,"strand"]
    multiplicity=collapsed[i,"n"]

     if(strand==0){
       nuc_position = start -region_start
     } else{
       nuc_position = end -region_start
     }

    genomic_nuc_position = region_start + nuc_position
      ##
      ## Calculate the window around this read (up/down stream is relative to the interval's strand)
      ##

      if ( strand == 0 ) {
        extended_start = nuc_position - extend_upstream
        extended_end   = nuc_position + extend_downstream
      } else {
        extended_start = nuc_position - extend_downstream
        extended_end   = nuc_position + extend_upstream
      }

      extended_start = ifelse( extended_start <= 0 , 1, extended_start)
      extended_end   = ifelse(extended_end > region_size-1,  region_size-1 ,extended_end+1)


      #print "extended region (start= $extended_start end= $extended_end nuc_pos=$nuc_position ) \n";
      ## Scan the window around (before/after) the nucleotide position
      for (window_index in extended_start:extended_end ) {
        #my ($sequences_count_same, $sequences_count_opposixte,$reads_count_same, $reads_count_opposite, $offset);
        sequences_count_same=0
        sequences_count_opposite=0
        reads_count_same=0
        reads_count_opposite=0
        offset=0


        ##
        ## Find the counts of reads/intervals in this window position
        ##
        if ( strand == 0 ) {
          sequences_count_same      = window_positive_sequences[window_index]
          reads_count_same          = window_positive_reads[window_index]
          sequences_count_opposite  = window_negative_sequences[window_index]
          reads_count_opposite      = window_negative_reads[window_index]

          offset = window_index - nuc_position
        } else {
          sequences_count_same =      window_negative_sequences[window_index]
          reads_count_same =          window_negative_reads[window_index]
          sequences_count_opposite =  window_positive_sequences[window_index]
          reads_count_opposite =      window_positive_reads[window_index]

          offset = nuc_position - window_index
        }
        #print "indice = $window_index  ---off =  $offset\n";
        # At offset 0, subtract the current interval from the counts
        # (only if NOT using external intervals )
        if ( offset == 0 )  {
          sequences_count_same =sequences_count_same - 1
          reads_count_same =reads_count_same - multiplicity
        }

        print(paste(start,end,extended_start,extended_end,nuc_position,genomic_nuc_position ,sequences_count_same, reads_count_same,offset,sep="#" ))

        if ( sequences_count_same != 0 || reads_count_same != 0 || sequences_count_opposite != 0 || reads_count_opposite != 0) {
         line = data.frame(ref,genomic_nuc_position,"+",  nuc_position,offset,   ## This is the window offset
                  sequences_count_same ,reads_count_same,sequences_count_opposite,reads_count_opposite);
          utils::write.table(line,"output.tab",append = TRUE,col.names = FALSE,row.names = FALSE)
          # result = rbind (result,c(ref,genomic_nuc_position,"+",  nuc_position,offset,   ## This is the window offset
        #  sequences_count_same ,reads_count_same,sequences_count_opposite,reads_count_opposite));

          #result = rbind (result,c(ref, "\t",genomic_nuc_position, "\t","+", "\t", nuc_position, "\t", offset, "\t",  ## This is the window offset
          #sequences_count_same, "\t",reads_count_same, "\t",sequences_count_opposite, "\t",reads_count_opposite));
        }
      }
    }
  }
