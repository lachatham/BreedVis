#' @description Runs PCA and makes biplot based on genotype values 
#'
#' @param genos           The $imputed output of A.mat from rrBLUP
#' @param color_by        How to color the PCA: either "user_defined" specified or based on some "pattern" in the line names
#' @param linegroups      A data frame, first olumn should be line names and second column specifies the grouping factor to use for coloring
#' @param pattern         A character vector specifying patterns to find and use for assigning a groupings.. i.e. pattern=c("pat1","pat2","pat3"). Last item is what to call the leftovers w/o a match
#' @param colors          Optional list of colors to plot with
#' @param scree           Whether to make a scree plot (T/F)
#' @param do.plotly       Whether to return the plot as a plotly version. takes longer esp for large data
#' 
#' @example 
#' \dontrun{
#' 
#' to install ggplot: 
#' devtools::install_github("vqv/ggbiplot")
#' 
#'   genos_pca(genos=genos$imputed,
#'   color_by="pattern",
#'   # color_by="user_defined",
#'   pattern=c("MN17","MN18","MN191","MN192","Founder"),
#'   # pattern=c("MN","SD","ND","WI","IL","misc"),
#'   linegroups=linegroups)
#'  # colors=c("pink","blue","green","black","yellow","red"))
#' 
#' }
#'  
#' 
#'@import tidyverse
#'@import ggbiplot
#'@import viridis
#'@import plotly
#'
#' @return                a plotly ggplot2 plot
#' 
#' 
genos_pca<-function(genos=genos,
        color_by=c("pattern"),
        linegroups=linegroups,
        pattern=c(),
        colors=NULL,
        scree=F,
        do.plotly=TRUE,
        size=2,
        pch=16,
        alpha=0.5){
  

## Color by genotyping file
    # require a dataframe with 2 columns (line, factor indicating original geno file (i.e. USDA, UMGC, F3))
if (color_by=="user_defined"){
  
  linegroups %>%
    dplyr::rename(line_name=1, linegroups=2)%>%
    unique()->linegroups
  
  genos %>% as.data.frame() %>%
    rownames_to_column(var="line_name") %>%
    unique() ->lines.genos
  
  if(!is.na(table(lines.genos$line_name %in% linegroups$line_name)[2])){
  warning("Not all genotype lines are in the supplied linegroups data frame")
  }
  
  if(!is.na(table(linegroups$line_name %in% lines.genos$line_name)[2])){
  warning("Not all lines in the supplied linegroups data frame are in the genotypes lines")
  }
  
genos %>% as.data.frame() %>%
  rownames_to_column(var="line_name")%>%
  left_join(linegroups,., by="line_name")%>%
  unique()->geno_pca
  
    # Show how many are in each group
    geno_pca %>%
      select(line_name,linegroups)%>%
      group_by(linegroups)%>%
      tally() %>%
      print()
    
       # PCA   
    geno_pca %>% 
      select(-line_name,-linegroups)%>%
      na.omit() %>%
      prcomp()->prin1

    # Set colors
    if(is.null(colors)){
      vircolorC<-viridis(10,alpha=0.7,begin=0,end=1, option="C")
      vircolorD<-viridis(10,alpha=0.7,begin=0,end=1, option="D")
      colors<-c(vircolorC[1],vircolorD[7], vircolorC[7],vircolorD[9],
                vircolorC[5],vircolorD[4],vircolorC[9],vircolorC[3],
                vircolorD[5],vircolorD[1],vircolorD[10])
      
    }
    
    ## Biplot
    p1<-ggbiplot(prin1,choices=c(1,2)) + 
      geom_text(aes(label=geno_pca$line_name,group=2),color=NA)+
      geom_point(aes(color=geno_pca$linegroups,text=geno_pca$line_name,group=1),alpha=alpha,pch=pch,size=size)+
      theme_classic() + scale_colour_manual(values=colors)+
      theme(legend.title=element_blank())
    p1$layers<-p1$layers[c(5,4)]
    
  
}
  
## Color by Generation
    # requires a character vector indicating the list of patterns you want to use, detected from line_name
    # (e.g. pattern=c("pattern1","pattern2","pattern3"))
    # the last item should be what to call any leftovers that dont match any of the previous entries
if (color_by =="pattern"){ #color by pattern is default


  ## Default is color by pattern, if none given, will plot everything as 1 color:  
  if(is.null(pattern)){
    warning("No pattern provided to match with line names. All lines will be plotted in the same color")
    
    genos %>% as.data.frame() %>%
      rownames_to_column(var="line_name")%>%
      unique()->geno_pca
    
        # PCA   
    
    geno_pca %>% 
      select(-line_name)%>%
      na.omit()%>%
      prcomp()->prin1
    
    # Set colors
    if(is.null(colors)){
      vircolorC<-viridis(10,alpha=0.7,begin=0,end=1, option="C")
      vircolorD<-viridis(10,alpha=0.7,begin=0,end=1, option="D")
      colors<-c(vircolorC[1],vircolorD[7], vircolorC[7],vircolorD[9],
                vircolorC[5],vircolorD[4],vircolorC[9],vircolorC[3],
                vircolorD[5],vircolorD[1],vircolorD[10])
      
    }
    
     p1<-ggbiplot(prin1,choices=c(1,2)) + 
      geom_text(aes(label=geno_pca$line_name,group=2),color=NA)+
      geom_point(aes(text=geno_pca$line_name,group=1),alpha=alpha,pch=pch,size=size)+
      theme_classic() + scale_colour_manual(values=colors)+ 
      theme(legend.position='none')
    p1$layers<-p1$layers[c(5,4)] 
    
    
    } # end of if no pattern given

# If pattern IS given.. this will create a new column in genos_pca that specifies the group to color by
## set a new column in geno_pca that specifies groups to color by based on some pattern, suffix, whatever in the line names
if(!is.null(pattern)){
  
  genos %>% as.data.frame() %>%
    rownames_to_column(var="line_name")%>%
    mutate(pattern=rep(NA))%>%
    unique()->geno_pca
  
  for(n in 1:length(geno_pca$line_name)){
for (p in 1:length(pattern)){
  if (str_detect(geno_pca$line_name[n],pattern[p])){
    geno_pca$pattern[n]<-pattern[p]
    break
  } # if you detect the pattern, apply name in pattern column and exit loop
  if (p == length(pattern)){
    geno_pca$pattern[n]<-pattern[p]
  } # if on the last pattern in list and it didn't get renamed above, then call it the last item in list
}
  }

    # Show how many are in each group
      geno_pca %>%
        select(line_name,pattern)%>%
        group_by(pattern)%>%
        tally() %>%
        print()
      
      # PCA   
      geno_pca %>% 
        select(-line_name,-pattern)%>%
        na.omit() %>%
        prcomp()->prin1
      
      # Set colors
       if(is.null(colors)){
        vircolorC<-viridis(10,alpha=0.7,begin=0,end=1, option="C")
        vircolorD<-viridis(10,alpha=0.7,begin=0,end=1, option="D")
        colors<-c(vircolorC[1],vircolorD[7], vircolorC[7],vircolorD[9],
                  vircolorC[5],vircolorD[4],vircolorC[9],vircolorC[3],
                  vircolorD[5],vircolorD[1],vircolorD[10])
     
       }
      
       ## Biplot
        p1<-ggbiplot(prin1,choices=c(1,2)) + 
        geom_text(aes(label=geno_pca$line_name,group=2),color=NA)+
        geom_point(aes(color=geno_pca$pattern,text=geno_pca$line_name,group=1),alpha=alpha,pch=pch,size=size)+
        theme_classic() + scale_colour_manual(values=colors)+ 
        theme(legend.title=element_blank())
      p1$layers<-p1$layers[c(5,4)]
      
}
     
}
 
  if(do.plotly==TRUE){
  p1<-ggplotly(p1)
  }
  
  return(p1)
  
  
  # Scree plot, if desired  
  if(scree==TRUE){
        windows()
        return(screeplot(prin1,type="lines",main="genos_PCA_scree"))
        # will open in a new window, so doesnt replace pca
  }
  rm(colors)
  rm(vircolorC)
  rm(vircolorD)
  rm(geno_pca)
}

### To Do:
####### Test with barley data
####### Get rid of warning messages:
############ Scale for 'colour' is already present. Adding another scale for 'colour', which will replace the existing scale.
############ Warning message:  Ignoring unknown aesthetics: text 
############ Change colors to colorblind friendly?
