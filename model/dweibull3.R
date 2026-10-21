

dmyweib3 <- nimbleFunction(
  run = function(x     = double(0),
                 shape = double(0),
                 scale = double(0),
                 shift = double(0),
                 log   = integer(0, default = 0)) {
    returnType(double(0))
    
    if(shape <= 0) { if(log == 1) return(-Inf) else return(0) }
    if(scale <= 0) { if(log == 1) return(-Inf) else return(0) }
    if(x <= shift)     { if(log == 1) return(-Inf) else return(0) }
    
    u <- x - shift
    z <- u / scale
    logdens <- log(shape) - log(scale) +
      (shape - 1.0) * log(z) -
      pow(z, shape)
    
    if(log == 1) return(logdens) else return(exp(logdens))
  }
)

rmyweib3 <- nimbleFunction(
  run = function(n     = integer(0),
                 shape = double(0),
                 scale = double(0),
                 shift = double(0)) {
    returnType(double(0))  # <-- scalar draw
    
    u <- runif(1, 0, 1)
    x <- shift + scale * pow(-log(1.0 - u), 1.0 / shape)
    return(x)
  }
)