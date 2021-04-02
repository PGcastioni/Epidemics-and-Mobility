
# Samples n integers numbers with variance p*sd.max and sum X_start. Be careful: the precision with which you get close to 
# you target sd depends on how large is n ( ex. if n = 100 you got a resolution of p = 0.01)
sample_sd <- function(n, X_start, p) {
    # Maximum value of standard deviation given a certain X_start and n
    sd.max <- sd( c( X_start, rep(0, n-1) ) ) 
    # Vector with minimum standard deviation
    x <- rep( floor( X_start / n ), n) + c( rep( 1, X_start %% n) , rep( 0, n - X_start %% n) )
    try <- x
    
    current_sd <- sd(try)
    target_sd <- p * sd.max
    while ( current_sd < target_sd ) {
        from <- sample.int( n , size = 1, prob = x > 0)
        to <- sample.int( n , size = 1, prob = x )
        try[from] <- x[from] - 1
        try[to] <- try[to] + 1
        current_sd <- sd(try)
        if ( current_sd > sd(x) ) {
            x <- try
        } else {
            try <- x
        }
    }
    return(x)
}