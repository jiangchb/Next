suppressPackageStartupMessages(library('optparse'))
option_list = list(
  make_option(c("--input","-i"),type = "character", 
              help= "Gene expression matrix, fpkm.xls; or GSVA result, GO_GSVA_enrichment.xls; or differential result, diffexp_result.xls"),
  make_option(c('--gmt', "-g"), type = "character",
              help = "gmt file"),
  make_option(c("--min_sz", "-s"), type = "integer",  default= 2,
              help = "Minimum size of the resulting gene sets."),
  make_option(c("--max_sz", "-S"), type = "integer", default= 99999,
              help = "Maximum size of the resulting gene sets."),
  make_option(c("--parallel", "-t"), type = "integer", default = 6,
              help = "Number of processors to use when doing the"),
  make_option(c("--pvalue","-p"), type = "character", default = 0.05,
              help = "The threshold of pvalue"),
  make_option(c("--FDR","-F"), type = "character", default = 1,
              help = "The threshold of Banjamimi-Hochberg adjusted pvalue"),
  make_option(c("--FoldChange","-f"), type = "character", default = 0,
              help = "The threshold of FoldChange"),
  make_option(c("--method", "-M"), type = "character", default = 'gsva',
              help = "Method to employ in the estimation of gene-set
                    enrichment scores per sample. By default this is set
                    to gsva (Hanzelmann et al, 2013) and other options
                    are ssgsea (Barbie et al, 2009), zscore (Lee et
                    al, 2008) or plage (Tomfohr et al, 2005)."),
  make_option(c("--kcdf","-k"), type = "character", default= "Gaussian",
              help = "Character string denoting the kernel to use during the
                    non-parametric estimation of the cumulative distribution
                    function of expression levels across samples when method='gsva'.
                    The option can be Guassian, Possion and none.
                    By default, kcdf='Gaussian' which is suitable when input
                    expression values are continuous, such as microarray
                    fluorescent units in logarithmic scale, RNA-seq log-
                    CPMs, log-RPKMs or log-TPMs. When input expression
                    values are integer counts, such as those derived from
                    RNA-seq experiments, then this argument should be set
                    to kcdf='Poisson'. This argument supersedes arguments
                    rnaseq and kernel, which are deprecated and will be
                    removed in the next release."),  
  make_option( c( "--mx_diff", "-x" ), type = "logical", default = T,
              help = "Offers two approaches to calculate the enrichment (ES)
                    from the KS random walk statistic.mx_diff=FALSE: ES is
                    calculated as the maximum distance of the random walk
                    from 0. mx_diff=TRUE (default): ES is calculated as the
                    magnitude difference between the largest positive and
                    negative random walk deviations."),
  make_option(c("--OUTDIR", "-o"), type = "character",
              help = "Ouput directory"),
  make_option(c("--gsva_fn", "-r"), type = "character", default = "GSVA_enrichement_scores.xls",
              help = "filename of gsva output"),
  #For differential analysis
  make_option(c("--group", "-c"), type = "character",
              help = "Table of the group which each smaple belongs to, for example A,4,B,4" ),
  make_option(c("--limma_all", "-A"), type = "character", default = "allexp_genesets_GSVA_score.xls",
              help = "All differential output filename"),
  make_option(c("--limma_diff", "-D"), type = "character", default = "diffexp_genesets_GSVA_score.xls",
              help = "Significantly differential output filename"),
  #For plot
  make_option(c("--topn", "-N"), type = "integer", default = 100, 
              help = "the maximium number of ranked pathway terms on the top for each cluster"),
  make_option(c("--barplot", "-B"), default = "diffexp_genesets_GSVA_score",
              help = "Bar plot filname "),
  make_option(c("--height", "-H"), type = "double", default = 40,
              help = "the height of bar plot"),
  make_option(c("--width", "-W"), type = "double", default = 20,
              help = "the width of bar plot"),
  make_option(c("--sig_reg", "-G"), type = "character", default = "reg",
              help = "Color bars by significance or up_down, sig or reg")
)
opt_parser = OptionParser(option_list=option_list)
opts = parse_args(opt_parser)

######################### Output set ###################################
if ( is.null(opts$OUTDIR) ){
  print("NO output directory specified,the current directory will be used!")
  output_dir = getwd()
}else{
  if ( file.exists(opts$OUTDIR) ){
    output_dir = opts$OUTDIR
  }else{
    output_dir = opts$OUTDIR
    dir.create(output_dir, recursive = T)
  }
}
output_dir = normalizePath(output_dir )

if ( is.null(opts$limma_diff) ){
  limma_diff <- "diff_expression_result.xls"
}else{
  limma_diff= opts$limma_diff
}

if ( is.null(opts$limma_all) ){
  limma_all = "all_expression_result.xls"
}else{
  limma_all=opts$limma_all
}

if ( is.null(opts$gsva_fn) ){
  gsva_filename = "GSVA_enrichment_scores.xls"
}else{
  gsva_filename=opts$gsva_fn
}

############################## GSVA parameter ######################################
if ( is.null(opts$min_sz) ){
  min.sz = 2
}else{
  min.sz=as.numeric(opts$min_sz)
}

if( is.null(opts$max_sz)) {
  max.sz = 99999
} else {
  max.sz = as.numeric(opts$max_sz)
}

if ( is.null(opts$parallel_sz) ){
  parallel.sz= 6
}else{
  parallel.sz=as.numeric(opts$parallel_sz)
}

if ( !is.null(opts$method )){
  method = opts$method
}else{
  print("NO enrichment algrithms is AVAILABLE!GSVA will be used as default!")
  method = "gsva"
}

if ( is.null(opts$mx_dff) ){
  if(method == "gsva") {
    mx.diff= T
  } else {
    mx.diff= F
  }
}

if ( is.null( opts$kcdf) ){
  kcdf = "Gaussian"
}else{
  kcdf=opts$kcdf
}



################################ Limma parameter ######################################
if( is.null(opts$pval)) {
  pval = 0.05
} else {
  pval = as.numeric(opts$pval)
}

if( is.null(opts$FDR)) {
  adjp = 1
} else {
  adjp = as.numeric(opts$FDR)
}

if( is.null(opts$FoldChange)){
  fc = 0
  lfc = 0
  if(fc != 0){lfc = log2(0)}
} else {
  fc = as.numeric(opts$FoldChange)
  lfc = log2(fc)
}

####################### Visualize #############################
if( is.null(opts$topn)) {
  topn = 30
} else {
  topn = as.numeric(opts$topn)
}

if( is.null(opts$barplot)) {
  bar_plot = "diffexp_genesets_GSVA_score"
} else {
  bar_plot = as.character(opts$barplot)
}

if( is.null(opts$height)) {
  height = 40
} else {
  height = as.numeric(opts$height)
}

if( is.null(opts$width)) {
  width = 20
} else {
  width = as.numeric(opts$width)
}

if (is.null(opts$sig_reg)){
  sig_reg = 'sig'
} else {
  sig_reg = as.character(opts$sig_reg)
}


save_file <- function(df, output_dir, filename){
  result <- cbind(rownames(df),df)
  colnames(result)[1]<- 'geneset'
  output <- file.path(output_dir,filename)
  output <- file(output,'wb')
  write.table(result,output, quote = FALSE, col.names = TRUE, row.names = FALSE,sep = "\t")
  return(result)
}

# Convert Inf to min or max
InfConvert <- function(x){
  tmp <- x[x != "-Inf" & x != "Inf"]
  max = max(tmp)
  min = min(tmp)
  x <- as.numeric(sub("-Inf", min,x))
  x <- as.numeric(sub("Inf", max, x))
  return(x)
}

########################## GSVA #############################################
if(!is.null(opts$input) & !is.null(opts$gmt)){
  gmt <- opts$gmt
  expset <- opts$input
  suppressPackageStartupMessages(library('GSVA'))
  suppressPackageStartupMessages(library('GSEABase'))
  geneSet <- getGmt(gmt)
  #Get expset table and 
  expset <- read.table(expset, sep = '\t', header = TRUE, row.names = 1)
  if (kcdf == "Guassian"){
    log(expset)
    expset <- apply(expset, 2,InfConvert)
  }
  matrix_expset <- data.matrix(expset)
  #Get gsva enrichment matrix
  gsva_scores <- gsva(matrix_expset, geneSet, method = method, min.sz = min.sz , max.sz = max.sz, kcdf = kcdf, parallel.sz = parallel.sz)
  save_file(gsva_scores, output_dir, gsva_filename)

  #gsva_enrich_output <- file.path(output_dir,gsva_filename)
  #gsva_enrich_output <- file(gsva_enrich_output,'wb')
  #write.table(gsva_result,gsva_enrich_output, quote = FALSE, col.names = TRUE, row.names = FALSE,sep = "\t")
  print(paste0("GSVA enrichment scores result has been saved in ", file.path(output_dir,gsva_filename)))
}

############################# limma ######################
if (is.null(opts$group)){
    print("No group input")
} else {
  if (!exists("gsva_scores")){
    gsva_scores <-  read.table(opts$input, sep = '\t', quote = "", header = TRUE, row.names = 1)
  }
  group <- opts$group
  group <- strsplit(group,",")
  groups <- factor(c(rep(group[[1]][1],group[[1]][2]),rep(group[[1]][3],group[[1]][4])))
  suppressPackageStartupMessages(library("limma"))
  print("Differential expression analysis start")
  design <- model.matrix(~ groups + 0)
  print(colnames(gsva_scores))
  rownames(design) <- colnames(gsva_scores)
  case <- paste0("groups",group[[1]][1])
  control <- paste0("groups", group[[1]][3])
  comparE <- makeContrasts(paste0(case," - ", control), levels=design)
  # Limma linear fitness
  fit <- lmFit(gsva_scores, design)
  fit <- contrasts.fit(fit, comparE)
  fit <- eBayes(fit)
  # All results
  fit_all <- topTable(fit, coef=1,adjust = "BH",n = Inf)
  # Filter by thresholds
  fit_cutoff <- topTable(fit, coef=1,adjust = "BH",n = Inf, p.value = adjp)
  if(nrow(fit_cutoff) != 0){
    fit_cutoff <- cbind(rownames(fit_cutoff),fit_cutoff)
    fit_cutoff <-fit_cutoff[fit_cutoff[,"P.Value"] <= pval,]
    fit_cutoff <- fit_cutoff[fit_cutoff[,"logFC"] >= lfc | fit_cutoff[,"logFC"] <= lfc,]
    # Save files
    cutoff_output <- file.path(output_dir,limma_diff)
    cutoff_output <- file(cutoff_output,'wb')
    colnames(fit_cutoff) = c("geneset", "logFC", "avgExp", "t", "pval", "FDR", "B")
    write.table(fit_cutoff,cutoff_output, quote = FALSE, col.names = TRUE, row.names = FALSE,sep = "\t")
  }else{
    print("No pathway passes the threshold")
  }
  print(colnames(fit_all))
  colnames(fit_all) = c("logFC", "avgExp", "t", "pval", "FDR", "B")
  save_file(fit_all, output_dir,limma_all)
  print("Differential anaylysis has done!")
}

  
######################### Bar plot #####################################################
if(exists('fit_cutoff')){
    if(nrow(fit_cutoff)>0){
      diff_result <- fit_cutoff
    } else if (exists('fit_all')) {
      diff_result <- fit_all
      print(paste0("Plotting all genesets due to no significant geneset"))
    }
}  else if (!is.null(opts$input) & is.null(opts$gmt)){
  diff_result <- read.table(opts$input, sep = '\t', quote = "", header = TRUE)
} 
if (exists('diff_result')){
  print("Plotting...")
  suppressPackageStartupMessages(library(ggplot2))
  if (adjp != 1){
    print(adjp)
    diff_result$is.sig <- diff_result$FDR < adjp
    diff_result <- head(diff_result[order(diff_result[,"FDR"]),], topn)
  } else {
    print(pval)
    diff_result$is.sig <- diff_result$pval < pval
    diff_result <- head(diff_result[order(diff_result[,"pval"]),], topn)
  }
  diff_result$just = ifelse(diff_result$t<0,0,1)
  diff_result$up_down = ifelse(diff_result$t<0,'Down','Up')
  pp = ggplot(diff_result, aes(reorder(geneset, t), t)) 
  if(sig_reg == 'sig'){
    pp = pp + geom_col(aes(fill=is.sig),width=0.9) + scale_fill_manual(values=c("#6CC570","#2A5078"))
  } else {
    pp = pp + geom_col(aes(fill=up_down),width=0.9) + scale_fill_manual(values=c("#6CC570","#2A5078"))
  }
  pp = pp +  coord_flip() +
    labs(x="Gene set", y="t value of GSVA Score" ) +
    theme_minimal() +
    geom_text( aes(x= geneset, y=0, label = geneset), hjust = diff_result$just, size = 3.5) +
    theme(axis.text.y=element_blank()) +
    theme(panel.grid =element_blank())
  if(sig_reg == 'sig'){
    if (adjp == 1){
      pp = pp + labs( fill = paste0("p value < ", pval) )
    }else{
      pp = pp + labs( fill = paste0("FDR < ", adjp))
    }
  } else { 
    pp = pp + labs( fill = "Up_Down") + theme(axis.line.x = element_line(colour = "black"))+
    theme(axis.text.x=element_text(size=8))+theme(legend.text=element_text(size=8))
  }
  ggsave(file.path(output_dir,
                   paste0(bar_plot, ".pdf")),
         height = height,width = width,plot = pp, limitsize = FALSE)
  print(paste0(bar_plot,".pdf has done"))
  ggsave(file.path(output_dir,
                   paste0(bar_plot,".png")),
         height = height, width = width, plot = pp, limitsize = FALSE)
  print(paste0(bar_plot,".png has done"))
}


