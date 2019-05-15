variance.spots=function()
{
  environment(modules.CSV.sheets) <- environment()
  environment(modules.report.sheets) <- environment()
  environment(modules.profiles) <- environment()
  environment(modules.chromosomes) <- environment()
  environment(modules.relations) <- environment()

  
  
  ##### variance Spots ######
  uh <- matrix(log10(apply(metadata, 1, var)), preferences$dim.1stLvlSom, preferences$dim.1stLvlSom)

  peak.matrix <- matrix( FALSE, preferences$dim.1stLvlSom, preferences$dim.1stLvlSom )

  for( pos.x in 1:preferences$dim.1stLvlSom )
    for( pos.y in 1:preferences$dim.1stLvlSom )
    {
      neighbor.dists <- sapply( get.neighbors(pos.x, pos.y, preferences$dim.1stLvlSom), function(x) uh[x[1],x[2]] )
      peak.matrix[pos.x, pos.y] <- all( uh[ pos.x, pos.y] <= neighbor.dists )
    }
  
  spot.matrix <- matrix(0, preferences$dim.1stLvlSom, preferences$dim.1stLvlSom )
  for( sel.minimum in which(peak.matrix)[order( uh[ which(peak.matrix) ], decreasing=TRUE )] )
  {
    spot.members <- c( sel.minimum )
    spot.d <- c( mean( uh[spot.members] ) ) 
    
    while(TRUE)
    {
      spot.neighbors <- unique( unlist(   sapply( spot.members, get.neighbors, dim=preferences$dim.1stLvlSom )  ) )
      spot.neighbors <- setdiff( spot.neighbors, spot.members )
      
      expansion.d <- sapply( spot.neighbors, function(x) dist( metadata[c(x,sel.minimum),] ) )
      expansion <- spot.neighbors[ which.min( expansion.d ) ]
      
      if( ( mean( uh[c(spot.members,expansion)] ) < spot.d[length(spot.d)] && 
            spot.d[length(spot.d)] < spot.d[length(spot.d)-1] && 
            spot.d[length(spot.d)-1] < spot.d[length(spot.d)-2] ) || 
          length(spot.members) > preferences$dim.1stLvlSom^2*0.1 ) 
      {
        break
      } else
      {
        spot.members <- c( spot.members, expansion )
        spot.d <- c( spot.d, mean( uh[spot.members] ) ) 
      }
      
    }
#    if( mean(uh[spot.members]) > mean(uh) ) 
    spot.matrix[spot.members] <- max(spot.matrix,na.rm=TRUE)+1
  }
  spot.matrix[which(spot.matrix==0)] <- NA
  
  
  
  spot.list.variance <<- list()
  spot.list.variance$overview.map <<- uh
  spot.list.variance$overview.mask <<- rep(NA, preferences$dim.1stLvlSom ^ 2)
  spot.list.variance$spots <<- list()
  
  count.cluster <- 1
  for (i in seq_along(sort(unique(na.omit(as.vector(spot.matrix))))) )
  {
    spot.metagenes <- which(spot.matrix==i)
    spot.genes <- rownames(indata)[which(som.result$feature.BMU %in% spot.metagenes)]
    
    if (length(spot.genes) > 0)
    {
      spot.list.variance$overview.mask[spot.metagenes] <<- count.cluster
      spot.list.variance$spots[[LETTERS[count.cluster]]] <<- list()
      spot.list.variance$spots[[LETTERS[count.cluster]]]$metagenes <<- spot.metagenes
      spot.list.variance$spots[[LETTERS[count.cluster]]]$genes <<- spot.genes
      spot.list.variance$spots[[LETTERS[count.cluster]]]$mask <<- rep(NA, preferences$dim.1stLvlSom * preferences$dim.1stLvlSom)
      spot.list.variance$spots[[LETTERS[count.cluster]]]$mask[spot.metagenes] <<- 1
      
      spot.list.variance$spots[[LETTERS[count.cluster]]]$position <<-
        colMeans(apply(som.result$node.summary[spot.metagenes, 1:2]+1, 2, range))
      
      spot.list.variance$spots[[LETTERS[count.cluster]]]$beta.statistic <<-
        get.beta.statistic(set.data=metadata[spot.list.variance$spots[[LETTERS[count.cluster]]]$metagenes,,drop=FALSE],
                           weights=som.result$node.summary[spot.list.variance$spots[[LETTERS[count.cluster]]]$metagenes,]$n.features)
      
      count.cluster <- count.cluster + 1
    }
  }
  
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
  
  sig.spots <- which( apply( spot.list.variance$spotdata, 1, function(x) sd(x) > sd(spot.list.variance$spotdata,na.rm=T) ) )
  if( length(sig.spots) > 0 )
  {
    spot.list.variance$spots <<- spot.list.variance$spots[sig.spots]
    spot.list.variance$overview.mask[which(!spot.list.variance$overview.mask%in%sig.spots)] <<- NA
    spot.list.variance$overview.mask[!is.na(spot.list.variance$overview.mask  )] <<-
      match(spot.list.variance$overview.mask[!is.na(spot.list.variance$overview.mask)], sort(unique(na.omit(as.vector(spot.list.variance$overview.mask)))))
  }
  
  start.spot <- which.min( apply( sapply(spot.list.variance$spots, function(x) x$position ), 2, min ) )
  
  spot.arcs <- sapply(spot.list.variance$spots, function(x)
  {
    -atan2(x$position['y'] - preferences$dim.1stLvlSom / 2, x$position['x'] - preferences$dim.1stLvlSom / 2)
  })
  
  spot.arcs <- spot.arcs - spot.arcs[start.spot]
  
  if (any(spot.arcs<0))
  {
    spot.arcs[which(spot.arcs<0)] <- spot.arcs[which(spot.arcs<0)] + (2 * pi)
  }
  
  
  
  o <- order(spot.arcs)
  
  spot.list.variance$spots <<- spot.list.variance$spots[o]
  names(spot.list.variance$spots) <<- LETTERS[seq_along(spot.list.variance$spots)]
  
  spot.list.variance$overview.mask[!is.na(spot.list.variance$overview.mask  )] <<-
    match(spot.list.variance$overview.mask[!is.na(spot.list.variance$overview.mask)], sort(unique(na.omit(as.vector(spot.list.variance$overview.mask))))[o])
  
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

  #spot.list.group.underexpression2=spot.list.group.underexpression

  #rm(list=ls()[-which(ls()=="spot.list.group.underexpression2")])
  #detach(env)
  # load("ScoV.RData")
  #env$spot.list.group.underexpression=NULL
  #env$spot.list.group.underexpression=spot.list.group.underexpression2
  # save(env, file="ScoV.RData")

  # rm(list=ls())
  # load("ScoV.RData")
  #attach(env)
  #library(oposSOM)

  dirname <- file.path( paste(files.name, "- Results"),"Summary Sheets - Modules","Variance Spots" )
  util.info("Writing:", file.path(dirname, "*.pdf"))
  dir.create(dirname, showWarnings=FALSE, recursive=TRUE)

  modules.report.sheets(spot.list=spot.list.variance, main="Group neg Correlation Spots", path=file.path(dirname,"Report.pdf") )
  modules.profiles(spot.list=spot.list.variance,main="Group neg Correlation Spots", path=file.path(dirname,"Profiles.pdf") )
  modules.chromosomes(spot.list=spot.list.variance, main="Group neg Correlation Spots", path=file.path(dirname,"Chromosomes.pdf") )

}