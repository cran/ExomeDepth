


setClass("ExomeDepth",
         representation(test = "numeric",
                        reference = "numeric",
                        formula = "character",
                        expected = "numeric",
                        phi = "numeric",
                        likelihood = "matrix",
                        annotations = "data.frame",
                        CNV.calls = "data.frame"))




#############################################################################
setMethod("initialize", "ExomeDepth", function(.Object,
                                               data = NULL,
                                               test,
                                               reference,
                                               formula = 'cbind(test, reference) ~ 1',
                                               subset.for.speed = NULL) {
  if (length(test) != length(reference)) stop("Length of test and numeric must match")

  

  if (is.null(data)) data <- data.frame(intercept = rep(1, length(test)))

  if (!is.null(subset.for.speed)) {
    test <- test[ subset.for.speed ]
    reference <- reference[ subset.for.speed ]
    data <- subset(data, 1:nrow(data) %in% subset.for.speed)
  }

  require(aod)
  message('Now fitting the beta-binomial model: this step can take a few minutes.')
  mod <- betabin( data = data, formula = as.formula(formula), random = ~ 1, link = 'logit')
  
  .Object@formula <- formula
  .Object@test <- test
  .Object@reference <- reference
  .Object@expected <- aod::fitted(mod)
  .Object@phi <- mod@param[[ 'phi.(Intercept)']]
  .Object@annotations <- data.frame()

  message('Now computing the likelihood for the different copy number states')
  .Object@likelihood <- .Call("get_loglike_matrix",
                              phi = .Object@phi,
                              expected = .Object@expected,
                              total = as.integer(.Object@reference + .Object@test),
                              observed = as.integer(.Object@test))
  .Object
})


#############################################################################
if (!isGeneric("show")) {
  if (is.function("show"))
    fun <- show
  else fun <- function(object) standardGeneric("show")
  setGeneric("show", fun)
}

setMethod("show", "ExomeDepth", function(object) {
  cat('Number of data points: ', length(object@test), '\n')
  cat('Formula: ', object@formula, '\n')
  cat('Phi parameter: ', object@phi, '\n')
  if (ncol(object@likelihood) == 3) cat("Likelihood computed\n") else cat("Likelihood not computed\n")
})

#############################################################################


#############################################################################
#if (!isGeneric("ComputeLikelihood")) {
#  if (is.function("ComputeLikelihood"))
#    fun <- ComputeLikelihood
#  else fun <- function(object) standardGeneric("ComputeLikelihood")
#  setGeneric("ComputeLikelihood", fun)
#}

#setMethod("ComputeLikelihood", "ExomeDepth", function(object) {  
#  object@likelihood <- .Call("get_loglike_matrix",
#                        phi = object@phi,
#                        expected = object@expected,
#                        total = as.integer(object@reference + object@test),
#                        observed = as.integer(object@test),
#                        PACKAGE = 'ExoDepth')
#  dimnames(object@likelihood)[[2]] <- c('loglike.del', 'loglike.normal', 'loglike.dup')
#  object
#})

#############################################################################

#############################################################################
setGeneric("TestCNV", def = function(x, chromosome, start, end, type) standardGeneric('TestCNV'))

setMethod("TestCNV", "ExomeDepth", function(x, chromosome, start, end, type) {
  if (! type %in% c('deletion', 'duplication')) stop("type must be either duplication or deletion\n")
  if (length(chromosome) != 1 || length(start) != 1 || length(end) != 1 || length(type) != 1) stop("The arguments chromosome, start, end and type must all be of length 1")
  if (class(chromosome) == 'factor') chromosome <- as.character(chromosome)
  if (class(chromosome) != 'character') stop('The input chromosome must be a character or a factor')


  
  which.exons <- which((x@annotations$chromosome == chromosome) & (x@annotations$start >= start) & (x@annotations$end <= end))

  if (type == 'deletion') log.ratio <- sum(x@likelihood[ which.exons, 1] - x@likelihood[ which.exons, 2])
  if (type == 'duplication') log.ratio <- sum(x@likelihood[ which.exons, 3] - x@likelihood[ which.exons, 2])

  return  (log.ratio)
})




setGeneric("CallCNVs", def = function(x, chromosome, start, end, name, transition.probability = 0.0001) standardGeneric('CallCNVs'))


setMethod("CallCNVs", "ExomeDepth", function( x, chromosome, start, end, name, transition.probability) {

  if ( length(start) != length(chromosome) || length(end) != length(chromosome) || length(name) != length(chromosome) ) stop('Chromosome, start and end vector must have the same lengths.\n')
  if (nrow(x@likelihood) != length(chromosome) ) stop('The annotation vectors must have the same length as the data in the ExomeDepth x')

  ### Try to get the chromosome order right
  chr.names.used <- unique(as.character(chromosome))
  chr.levels <- c(as.character(seq(1, 22)), subset( chr.names.used, ! chr.names.used %in% as.character(seq(1, 22))))
  chr.levels <- subset(chr.levels, chr.levels %in% chr.names.used)
  
  x@annotations <- data.frame(name = name, chromosome = factor(chromosome, levels = chr.levels), start = start, end = end)
  my.new.order <-  order(x@annotations$chromosome, 0.5*(x@annotations$start + x@annotations$end) )

  if (sum( my.new.order != 1:nrow(x@annotations) ) > 0) {
    message('Positions of exons seem non ordered, so ExomeDepth will reorder the data according to chromosome and position')
    x@test <- x@test[ my.new.order ]
    x@reference <- x@reference[ my.new.order ]
    x@annotations <- x@annotations[ my.new.order, ]
    x@likelihood <- x@likelihood[ my.new.order, ]
  }

  
  total <- x@test + x@reference
  transitions <- matrix(nrow = 3, ncol = 3,
                        c( 1. - transition.probability, transition.probability/2., transition.probability/2.,
                          0.5, 0.5, 0.,
                          0.5, 0, 0.5),
                        byrow = TRUE)
  
  my.breaks <- which(diff(as.numeric(x@annotations$chromosome)) != 0) + 1
  x@likelihood[ my.breaks,1 ] <- - Inf
  x@likelihood[ my.breaks,3  ] <- - Inf
  
  my.calls <- viterbi.hmm (transitions, loglikelihood = x@likelihood[, c(2, 1, 3)], positions = NA)

  ################################ Now make it look better, add relevant info
  if (nrow(my.calls$calls) > 0) {

    my.calls$calls$start <- x@annotations$start[ my.calls$calls$start.p ]
    my.calls$calls$end <- x@annotations$end[ my.calls$calls$end.p ]
    my.calls$calls$chromosome <- as.character(x@annotations$chromosome[ my.calls$calls$start.p ])
  
    my.calls$calls$id <- paste('chr', my.calls$calls$chromosome, ':',  my.calls$calls$start, '-',  my.calls$calls$end, sep = '')
    my.calls$calls$type <- c('deletion', 'duplication')[ my.calls$calls$type ]
  
########## make things pretty
    my.calls$calls$BF <- NA
    my.calls$calls$reads.expected <- NA
    my.calls$calls$reads.observed <- NA


    for (ir in 1:nrow(my.calls$calls)) {
      
      if (my.calls$calls$type[ir] == 'duplication') my.calls$calls$BF[ir] <-  sum(x@likelihood [ my.calls$calls$start.p[ir] : my.calls$calls$end.p[ir],3 ] - x@likelihood [ my.calls$calls$start.p[ir] : my.calls$calls$end.p[ir],2 ])
      
      if (my.calls$calls$type[ir] == 'deletion') my.calls$calls$BF[ir] <-  sum(x@likelihood [ my.calls$calls$start.p[ir] : my.calls$calls$end.p[ir], 1 ] - x@likelihood [ my.calls$calls$start.p[ir] : my.calls$calls$end.p[ir],2  ])
      
      my.calls$calls$reads.expected[ ir ] <-  sum( total [my.calls$calls$start.p[ir] : my.calls$calls$end.p[ir] ] * x@expected [my.calls$calls$start.p[ir] : my.calls$calls$end.p[ ir ] ])
      my.calls$calls$reads.observed[ ir ] <-  sum( x@test [my.calls$calls$start.p[ir] : my.calls$calls$end.p[ir] ] )      
    }
    
    my.calls$calls$reads.expected <- as.integer( my.calls$calls$reads.expected)
    my.calls$calls$reads.ratio <-  signif(my.calls$calls$reads.observed / my.calls$calls$reads.expected, 3)
    my.calls$calls$BF <- signif( log10(exp(1))*my.calls$calls$BF, 3)
  }
  
  x@CNV.calls <- my.calls$calls
  return (x)

})

