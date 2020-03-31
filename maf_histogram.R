#' @description Runs PCA and makes biplot based on genotype values 
#'
#' @param genos           Genotype matrix, where rownames are lines and columns are markers (n x m)
#' @param linegroups      Dtaframe including lines as first column and their desired grouping factor in the 2nd col
#' @param pattern         A character vector specifying patterns to find and use for assigning a groupings.. 
#'                        i.e. pattern=c("pat1","pat2","pat3"). Last item is what to call the leftovers w/o a match
#'                        assumes line names are named with some sore of meaningful prefix, suffix, etc.. 
#' @param overlaid        T/F, whether to overlay hist or plot separately
#' @param commonmarkers   marker set should be subsetted for common marker prior to putting into genos arg... 
#'                        unless groups correspond to genotyping runs 
#'                        Turning common markers to TRUE is only helpful if groups are different genotyping runs. 
#'                        i.e. if groups (via linegroups or pttern) combine multiple genotyping runs (different sets of markers)
#'                        then all will be used and MAF is calculated only from the non NA entries, not all indiiduals in the 
#'                        desired groups... resulting MAF is NOT correct.
#' @param do.plotly       Whether to return the plot as a plotly version. takes longer esp for large data
#' @param MAFcutoff       MAF to use in calculating the num to be removed and printed when stats=T
#' @param includestats    Whether to print the mean MAF and number of markers that would be removed with a given MAF
#' @param filter          Whether to filter out the markers with maf around 0.5 (if not, difficult to see differences in bins <0.5)
#' @param do.plotly       Whether to return the plot as a plotly version. takes longer esp for large data
#' @param colors          User specified colors                                                                                                      
#'
#' 
#' @example 
#' \dontrun{
#' 
#
#' }
#'   
#'@import tidyverse
#'@import adegenet
#'@import viridis
#'@import plotly
#'@import gridExtra
#'@import purrr
#'@import ggbiplot
#'
#' @return                a (plotly) ggplot2 plot (overlaid=T) or a 
#'                        grid.arrange of multiple plots (overladid=F)
 maf_histogram<-function(genos=genos,
                        linegroups=NULL,
                        pattern=c(),
                        overlaid=FALSE,
                        commonmarkers=FALSE,
                        MAFcutoff=0.05,
                        includestats=TRUE,
                        filter=TRUE,
                        do.plotly=FALSE,
                        colors=NULL,
                        alpha=0.3,
                        binwidth=0.02){
library(ggbiplot)
library(tidyverse)
library(plotly)
library(viridis)
library(purrr)
library(gridExtra)
library(adegenet)
   
   
   # convert to format for adegenet  
   mydata1<-df2genind(genos, sep="")
   
   # function I found on stack overflow for conditional piping
   pif <- function(x, p, true, false = identity){
     if(!requireNamespace("purrr")) 
       stop("Package 'purrr' needs to be installed to use function 'pif'")
     
     if(inherits(p,     "formula"))
       p     <- purrr::as_mapper(
         if(!is.list(x)) p else update(p,~with(...,.)))
     if(inherits(true,  "formula"))
       true  <- purrr::as_mapper(
         if(!is.list(x)) true else update(true,~with(...,.)))
     if(inherits(false, "formula"))
       false <- purrr::as_mapper(
         if(!is.list(x)) false else update(false,~with(...,.)))
     
     if ( (is.function(p) && p(x)) || (!is.function(p) && p)){
       if(is.function(true)) true(x) else true
     }  else {
       if(is.function(false)) false(x) else false
     }
   } 
   
   # If no linegroups are provided, determine groups based on pattern
   if(is.null(linegroups)){
     ## split my data into different group based on some pattern or given list of names
     rownames(mydata1$tab)%>% 
       as.data.frame() %>%
       mutate(group=rep(NA))->linegroups
     names(linegroups)[1]<-c("line_name")
     
     for(n in 1:length(linegroups$line_name)){
       for (t in 1:length(pattern)){
         if (str_detect(linegroups$line_name[n],pattern[t])){
           linegroups$group[n]<-pattern[t]
           break
         } # if detect the pattern, apply name in pattern column and exit loop
         if (t == length(pattern)){
           linegroups$group[n]<-pattern[t]
         } # if on the last pattern in list and it didn't get renamed above, then call it the last item in list
       } # end of p loop
     }# end n loop
   } #end if linegroups is null
   
   # If linegroups is provided
   if(!is.null(linegroups)){
     linegroups %>%
       rename(line_name=1,group=2)->linegroups
   } ## if linegroups is provided, just change the names so the rest of this works
   
   # set the number of / names of groups
   levels(as.factor(linegroups$group))->groups
   
   names(mydata1$loc.n.all)%>%
     as.data.frame()%>%
     rename("line"=1)->allmafs
   
   # calculate MAF for each of the specified groups
   for (g in 1:length(groups)){
     selected.lines<- linegroups$line_name[linegroups$group==groups[g]]
     mf<-minorAllele(mydata1[row.names(mydata1$tab)%in% selected.lines])
     
     mf %>%
       as.data.frame() %>%
       rownames_to_column(var="line")->mf
     names(mf)[2]<-groups[g]
     mf[,2]<-as.numeric(mf[,2])
     
     allmafs<-merge(allmafs,mf)
     # note.. not all markers will have a value, some markers are missing from some groups. 
     # make an option to only use common markers.
     
   }#end of g loop
   
   cat("There are", dim(allmafs)[1], "total markers in the data, of which",
       dim(allmafs[complete.cases(allmafs),])[1], "are common among all groups.")
   
   # Common markers T/F
   if (commonmarkers==T){
     allmafs<-allmafs[complete.cases(allmafs),]
     cat("Trimming down to", dim(allmafs[complete.cases(allmafs),])[1],"markers.\n")
     warning("If groups do not correspond to separate genotyping runs, this may not be correct. \n
          Subset marker data in genos command to include common markers between genotyping runs.")
   }
   
   if (commonmarkers==F){
     cat("Using ", dim(allmafs)[1], "markers.\n")
   }
   
   # Calculate mean maf
   allmafs %>%
      pivot_longer(2:(length(groups)+1), names_to=c("group"))%>%
      pif(filter==TRUE, filter(.,value<0.5),) %>%
     pivot_wider(names_from=group) %>%
      select(-1) %>%
     summarise_all(mean, na.rm=TRUE) %>%
     as.data.frame() %>%
     pivot_longer(1:length(groups)) %>%
     rename(mean=value,group=name)->maf_means
   
   # Calculate how many markers would be removed at a given MAF cutoff
   allmafs %>%
     pivot_longer(2:(length(groups)+1), names_to=c("group"))%>%
      pif(filter==TRUE, filter(.,value<0.5),) %>%
     filter(value<MAFcutoff)%>%
     group_by(group)%>%
     tally()-> maf_markers
   
   # merge these ^ into "stats"
   merge(maf_markers,maf_means)->stats
   stats$mean<-round(stats$mean,3)
   names(stats)[2]<-paste0("n<",MAFcutoff)
   
   print(stats)
   
   
   
   # if overlaid=T
   if(overlaid==TRUE){
     
     # make histograms, with each group as a diff color, overlaid
     allmafs %>%
       pivot_longer(2:(length(groups)+1), names_to=c("group"))%>%
       pif(filter==TRUE, filter(.,value<0.5),) %>%
       ggplot(aes(x=value,fill=group))+
       geom_histogram(alpha=alpha,position="identity",binwidth = binwidth,color="black")+
       # scale_color_viridis(option="D",discrete = T)+
       theme_classic()+
       theme(legend.title=element_blank(), legend.position="top",
             legend.spacing.x=unit(1,"mm"),legend.direction="horizontal")+
       guides(fill=guide_legend(label.position="top",nrow=1))->p
     
     if(is.null(colors)){
       p+  scale_fill_viridis(option="D",discrete = T) ->p}
     if(!is.null(colors)){
       p+ scale_fill_manual(values=colors)->p
     }
     
     
     if(includestats==T){
       p+annotation_custom(tableGrob(t(stats),
                                     theme=ttheme_minimal(base_size = 9,
                                                          padding=unit(c(1,1),"mm"))),
                           0.25,0.25,max(layer_data(p)$y))->p
     }
     
     # if plotly version is desired:
     # note: pltoly does not support horizontal legend or custom annotation
     # includestates=T will NOT work if done in plotly
     if(do.plotly==TRUE){
       return(ggplotly(p))
     }else{
       return(p)
     }
   } #end if overlaid==T
   
   
   # if overlaid=F
   if (overlaid==FALSE){
     # make histograms, with each group as a separate graph
     plots<-list()
     
     if(is.null(colors)){
       plotcolors<-viridis(length(groups),alpha=alpha,begin=0,end=1, option="D")}
     if(!is.null(colors)){
       plotcolors<-colors
     }
     
     for (h in 1:length(groups)){
       
       allmafs %>%
         select(1,h+1) %>%
         rename(value=2) %>%
         pif(filter==TRUE, filter(.,value<0.5),) %>%
         ggplot(aes(x=value, fill=groups[h]))+
         geom_histogram(alpha=alpha,position="identity",
                        binwidth = binwidth,color="black", fill=plotcolors[h])+
         theme_classic()+
         theme(legend.position="none") + 
         labs(title=(paste0((groups[h]))),
              x="MAF", y="Number of markers")->plots[[h]]
     } # end of h loop (making the hists for each group)
     
     # if table with means and number of markers excldued b/c <MAF is desired
     # adds table to list of plots and puts in last spot in grid arrange
     if(includestats==T){
       plots[[h+1]]<-(tableGrob((stats),theme=ttheme_minimal(base_size=9,padding = unit(c(1,1),"mm"))))
     } #end if include stats=T
     
     ncol<-ceiling(sqrt(length(groups)))
     return(do.call("grid.arrange",c(plots, ncol=ncol)))
     #do.call("grid.arrange",c(plots, ncol=ncol))
   } # end if overlaid=FALSE
   
   
  } # end of function


