run.som <- function(modpar=modpar,modsom=modsom)
{
  library(oposSOM)
  library(igraph)
  library(ape)
  library(tsne)
  library(fastICA)

  #source("multiSOMe.r")


  env <- opossom.new(list(dataset.name = modpar$dataset.name,

                          dim.1stLvlSom = modpar$dim.1stLvlSom,
                          dim.2ndLvlSom = 20,

                          training.extension = 1,
                          rotate.SOM.portraits = 0,
                          flip.SOM.portraits = F,

                          database.biomart = "ENSEMBL_MART_ENSEMBL",
                          database.host = "aug2017.archive.ensembl.org",#"www.ensembl.org",
                          database.dataset = "hsapiens_gene_ensembl",
                          database.id.type =modpar$return.ids,

                          activated.modules = list( "reporting" = TRUE,
                                                    "primary.analysis" = TRUE,
                                                    "sample.similarity.analysis" = TRUE,
                                                    "geneset.analysis" = modpar$geneset.analysis,
                                                    "geneset.analysis.exact" = FALSE,
                                                    "group.analysis" = TRUE,
                                                    "difference.analysis" = TRUE ),

                          standard.spot.modules = "dmap",

                          spot.coresize.modules = 3,
                          spot.threshold.modules = 0.95,
                          spot.coresize.groupmap = 3,
                          spot.threshold.groupmap = 0.75,

                          feature.centralization = F,
                          sample.quantile.normalization = F,

                          pairwise.comparison.list = list() ) )


  # Load input data
  env$indata <-modsom$indata

  # Define sample groups
  env$group.labels <- modsom$group.labels


  tripleSOMe.run(env)
  eachOME.run(env)
}
