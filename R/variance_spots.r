variance.spots=function()
{
  environment(modules.CSV.sheets) <- environment()
  environment(modules.report.sheets) <- environment()
  environment(modules.profiles) <- environment()
  environment(modules.chromosomes) <- environment()
  environment(modules.relations) <- environment()

  uh <- matrix(log10(apply(metadata, 1, var)), preferences$dim.1stLvlSom, preferences$dim.1stLvlSom)

  metadata.scaled <- (uh-min(uh)) / (max(uh)-min(uh))

  ##### variance Spots ######

  # extract sample modules
  sample.spot.list <- list()
  sample.spot.core.list <- list()
  n.sample.modules <- 0
  preferences$spot.threshold.modules=0.62
  preferences$spot.coresize.modules=2
  #for (m in 1:ncol(indata))
  #{
  # define bigger core regions
  core <- matrix(NA, preferences$dim.1stLvlSom, preferences$dim.1stLvlSom)
  core[which(metadata.scaled> preferences$spot.threshold.modules)] <- -1

  spot.i = 0

  while (nrow(which(core == -1, arr.ind=TRUE)) > 0)
  {
    start.pix <- which(core == -1, arr.ind=TRUE)[1,]
    spot.i <- spot.i + 1
    core <- col.pix(core, start.pix[1], start.pix[2], spot.i, preferences$dim.1stLvlSom)
  }

  # shrink each separate region to core size
  for (s.i in 1:max(core, na.rm=TRUE))
  {
    if (sum(core == s.i, na.rm=TRUE) > preferences$spot.coresize.modules)
    {
      core.metagenes <- which(core == s.i)
      o <- order(uh[core.metagenes], decreasing=TRUE)[1:preferences$spot.coresize.modules]
      core[setdiff(core.metagenes, core.metagenes[o])] <- NA
    }
  }

  core[which(!is.na(core))] <- -1
  spot.i <- 0

  while (nrow(which(core == -1, arr.ind=TRUE)) > 0)
  {
    start.pix <- which(core == -1, arr.ind=TRUE)[1,]
    spot.i <- spot.i + 1
    core <- col.pix(core, start.pix[1], start.pix[2], spot.i, preferences$dim.1stLvlSom)
  }

  # define spot area around cores
  for (s.i in 1:max(core,na.rm=TRUE))
  {
    n.sample.modules <- n.sample.modules + 1

    sample.spot.core.list[[n.sample.modules]] <-
      matrix(NA, preferences$dim.1stLvlSom, preferences$dim.1stLvlSom)

    sample.spot.core.list[[n.sample.modules]][which(core == s.i)] <-
      uh[which(core == s.i)]

    spot <- matrix(NA, preferences$dim.1stLvlSom, preferences$dim.1stLvlSom)
    spot[which(metadata.scaled > preferences$spot.threshold.modules)] <- -1

    start.pix <- which(!is.na(sample.spot.core.list[[n.sample.modules]]), arr.ind=TRUE)
    start.pix <- start.pix[which.max(sample.spot.core.list[[n.sample.modules]][start.pix]),]

    spot <- col.pix(spot, start.pix[1], start.pix[2], 1, preferences$dim.1stLvlSom)

    sample.spot.list[[n.sample.modules]] <- matrix(NA, preferences$dim.1stLvlSom, preferences$dim.1stLvlSom)
    sample.spot.list[[n.sample.modules]][which(spot == 1)] <- metadata.scaled[which(spot == 1)]
  }


  # filter
  remove <- c()

  for (i in 1:n.sample.modules)
  {
    if (sum(!is.na(sample.spot.list[[i]])) <= 1)
    {
      # empty, i.e. does not exceed threshold -> remove
      remove <- c(remove, i)
    } else if (sum(!is.na(sample.spot.list[[i]])) < sum(!is.na(sample.spot.core.list[[i]])))
    {
      # core larger than spot -> shrink core
      sample.spot.core.list[[i]] <- sample.spot.list[[i]]
    }
  }

  if (length(remove) > 0)
  {
    sample.spot.list <- sample.spot.list[-remove]
    sample.spot.core.list <- sample.spot.core.list[-remove]
    n.sample.modules <- length(sample.spot.list)
  }

  o <- order(sapply(sample.spot.core.list, function(x) { mean(x, na.rm=TRUE) }), decreasing=TRUE)
  sample.spot.list <- sample.spot.list[o]
  sample.spot.core.list <- sample.spot.core.list[o]


  ## define variance spots ##
  spot.list.variance <<- list()

  spot.list.variance$overview.map <<-uh

  spot.list.variance$overview.mask <<- rep(NA, preferences$dim.1stLvlSom ^ 2)
  spot.list.variance$spots <<- list()

  for (i in seq_along(sample.spot.list))
  {
    spot.metagenes <- which(!is.na(sample.spot.list[[i]]))
    spot.genes <- rownames(indata)[which(som.result$feature.BMU %in% spot.metagenes)]

    if (length(spot.genes) > 0)
    {
      spot.list.variance$overview.mask[which(!is.na(sample.spot.list[[i]]))] <<- i
      spot.list.variance$spots[[LETTERS[i]]] <<- list()
      spot.list.variance$spots[[LETTERS[i]]]$metagenes <<- spot.metagenes
      spot.list.variance$spots[[LETTERS[i]]]$genes <<- spot.genes
      spot.list.variance$spots[[LETTERS[i]]]$mask <<- rep(NA, preferences$dim.1stLvlSom * preferences$dim.1stLvlSom)
      spot.list.variance$spots[[LETTERS[i]]]$mask[spot.metagenes] <<- 1

      spot.list.variance$spots[[LETTERS[i]]]$position <<-
        colMeans(apply(som.result$node.summary[spot.metagenes, 1:2] + 1, 2, range))

      spot.list.variance$spots[[LETTERS[i]]]$beta.statistic <<-
        get.beta.statistic(set.data=metadata[spot.list.variance$spots[[LETTERS[i]]]$metagenes,,drop=FALSE],
                           weights=som.result$node.summary[spot.list.variance$spots[[LETTERS[i]]]$metagenes,]$n.features)
    }
  }

  o <- order(sapply(spot.list.variance$spots, function(x)
  {
    which.max(apply(metadata[x$metagenes,,drop=FALSE], 2, mean))
  }))

  spot.list.variance$spots <<- spot.list.variance$spots[o]
  names(spot.list.variance$spots) <<- LETTERS[seq_along(spot.list.variance$spots)]

  spot.list.variance$overview.mask[!is.na(spot.list.variance$overview.mask)] <<-
    match(spot.list.variance$overview.mask[!is.na(spot.list.variance$overview.mask)], o)

  spot.list.variance$spotdata <<-
    t(sapply(spot.list.variance$spots, function(x)
    {
      if (length(x$genes > 0))
      {
        colMeans(indata[x$genes,,drop=FALSE],na.rm=TRUE)
      } else
      {
        rep(0, ncol(indata))
      }
    }))

  colnames(spot.list.variance$spotdata) <<- colnames(indata)


  spot.fisher.p <- function(spot)
  {
    spot$Fisher.p <- GeneSet.Fisher(unique(na.omit(gene.info$ids[spot$genes])),
                                    unique.protein.ids, gs.def.list, sort=TRUE)

    return(spot)
  }

  spot.list.variance$spots =lapply( spot.list.variance$spots, spot.fisher.p)

  dirname <- file.path( paste(files.name, "- Results"),"Summary Sheets - Modules","Variance Spots" )
  util.info("Writing:", file.path(dirname, "*.pdf"))
  dir.create(dirname, showWarnings=FALSE, recursive=TRUE)

  modules.report.sheets(spot.list=spot.list.variance, main="Variance Spots", path=file.path(dirname,"Report.pdf") )
  modules.profiles(spot.list=spot.list.variance,main="Variance Spots", path=file.path(dirname,"Profiles.pdf") )
  modules.chromosomes(spot.list=spot.list.variance, main="Variance Spots", path=file.path(dirname,"Chromosomes.pdf") )

}
