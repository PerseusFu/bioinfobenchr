# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'
#' @import Seurat SeuratObject tidyverse Matrix ggsci ggpubr harmony SeuratData
#'








hellobioinfobenr <- function() {
  print("Hello, world! by Y.P.Fu 2021.11.04")
}




#scRNAseq functions

#' Pre-progress for seurat object.
#' @description This function include 3steps: percentfeature, normalize and reducion.
#'
#' @param s.obj seurat object
#' @param reduce_method choose from "pca","umap","tsne"
#' @param batch_label if there may be batch effects, use meta data label
#' @param label if output with label
#'
#' @return DimPlot, type depend on para "reduce_method"
#' @export
#'
#' @examples
#' data("pbmc3k")
#' scrna_pre_process(pbmc3k)
#'
scrna_pre_process <- function(sobj,reduce_method="pca",batch_label=NULL,label=F){
  print("For any further help, mail to fan_pku@163.com")
  print("Step 1 Percentage Features Detected")
  sobj <- Seurat::PercentageFeatureSet(sobj, pattern = "^MT-", col.name = "percent.mt")
  sobj <- Seurat::PercentageFeatureSet(sobj, pattern = "^RPS", col.name = "percent.rps")
  sobj <- Seurat::PercentageFeatureSet(sobj, pattern = "^RPL", col.name = "percent.rpl")
  sobj <- Seurat::PercentageFeatureSet(sobj, pattern = "^RNA\\d8S5", col.name = "percent.rrna")
  print("Step 2 Normalization and Scale")
  sobj <- Seurat::NormalizeData(sobj,verbose=F)
  sobj <- Seurat::FindVariableFeatures(sobj,selection.method = "vst", nfeatures = 2000)
  sobj<- Seurat::ScaleData(sobj, features = rownames(sobj))
  sobj <- Seurat::ScaleData(sobj, vars.to.regress = c("percent.rps", "percent.rpl", "nCount_RNA","percent.rrna", "nFeature_RNA"))
  print("Step 3 Run Reduction")
  sobj <- Seurat::RunPCA(sobj, verbose = F)
print(reduce_method)

  if(reduce_method == "pca"){
    print(Seurat::DimPlot(sobj,label=label,reduction = "pca"))
    return("Steps Done")
    }
  else{
    if( !is.null(batch_label)){
     sobj<-harmony::RunHarmony(s.obj,batch_label)
    }
  sobj<- FindNeighbors(sobj, reduction = "pca",verbose = FALSE)
  sobj<- FindClusters(sobj,verbose = FALSE)
  if(reduce_method=="umap"){
  sobj <- RunUMAP(sobj,dims = 1:20,reduction = "pca",verbose = F)
  print(Seurat::DimPlot(sobj,reduction = "umap",label=label))}
  if(reduce_method=="tsne"){
  sobj <- RunTSNE(sobj,dims = 1:20,reduction = "pca",verbose = F)
  print(Seurat::DimPlot(sobj,reduction = "tsne",label=label))}
  }


}



#' Sample Seurat Object
#' @description  can provide a subset which are randomly select from a given Seurat Object with a very rate.
#' @param sobj seurat object
#' @param sample_rate sample rate. default is 1
#' @param rep whether want repeated, default is False
#'
#' @return
#' @export
#'
#' @examples
#' data("pbmc3k")
#' sample_seurat(pbmc3k,sample_rate=0.5)
#'
sample_seurat <- function(sobj,sample_rate=1,rep=F){
  cell.ident <- colnames(sobj)
  sample_cell.ident <- sample(cell.ident,sample_rate,replace = rep)
  return(sobj[sample_cell.ident,])
}



#other functions
#' Sparse Matrix change to Normal Matrix
#' @description  this function turns sparse matrix into norm matrix, especially for those matrixs which can't be changed by as.matrix()
#' @param mat
#'
#' @return
#' @export
#'
#' @examples
#' i <- c(1,3:8);
#' j <- c(2,9,6:10);
#' x <- 7 * (1:7)
#' (A <- sparseMatrix(i, j, x = x))
#' A.norm <- as_matrix(A)
#'
as_matrix <- function(mat){

  tmp <- matrix(data=0L, nrow = mat@Dim[1], ncol = mat@Dim[2])

  row_pos <- mat@i+1
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])+1
  val <- mat@x

  for (i in seq_along(val)){
    tmp[row_pos[i],col_pos[i]] <- val[i]
  }

  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}
