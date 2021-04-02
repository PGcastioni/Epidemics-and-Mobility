sample_kind <- function(x) {
    if (sum(x) == 0) {
        a<-0
    } else {
        a <- sample( c('s', 'x', 'r'), size = 1, replace = F, prob = x )
    }
    return(a)
}

library(deSolve)

# MASTER EQUATION =======================================================================
SIR_ode.iter <- function(times, yini, params){
    beta <- params[[1]] 
    mu <- params[[2]]
    eps <- params[[3]]
    R_ij <- params[[4]]
    dt <- params[[5]]
    n <- dim(R_ij)[1]
    
    s_i <- yini[1:n]
    x_i <- yini[(n+1):(2*n)]
    r_i <- yini[(2*n+1):(3*n)]
    n_i <- s_i + x_i + r_i
    
    ds <- s_i - dt*( beta * x_i / n_i * s_i ) + dt*( s_i %*% R_ij - s_i )/eps
    dx <- x_i + dt*( beta * x_i / n_i * s_i - mu * x_i ) +  dt*( x_i %*% R_ij - x_i )/eps
    dr <- r_i + dt*(mu * x_i )  +  dt*( r_i %*% R_ij - r_i  )/eps
    
    ds <- as.vector(ds)
    dx <- as.vector(dx)
    dr <- as.vector(dr)
    
    return(list( c(ds, dx, dr) ))
}

SIR_ode2.iter <- function(times, yini, params){
    beta <- params[[1]] 
    mu <- params[[2]]
    eps <- params[[3]]
    R_ij <- params[[4]]
    dt <- params[[5]]
    n <- dim(R_ij)[1]
    
    # To avoid problems with this differential equation one should always chose the time step "dt" so that the product of it with
    # any rate ("beta", "mu" or "R_ij") is lower than one
    
    s_i <- yini[1:n]
    x_i <- yini[(n+1):(2*n)]
    r_i <- yini[(2*n+1):(3*n)]
    n_i <- s_i + x_i + r_i
    
    ds <- s_i - dt*( beta * x_i / n_i * s_i )*min(1,eps) + dt*( s_i %*% R_ij - s_i )*min(1,1/eps)
    dx <- x_i + dt*( beta * x_i / n_i * s_i - mu * x_i )*min(1,eps)  +  dt*( x_i %*% R_ij - x_i )*min(1,1/eps)
    dr <- r_i + dt*(mu * x_i )*min(1,eps)  +  dt*( r_i %*% R_ij - r_i  )*min(1,1/eps)
    
    ds <- as.vector(ds)
    dx <- as.vector(dx)
    dr <- as.vector(dr)
    
    return(list( c(ds, dx, dr) ))
}

SIR_ode2 <- function(net, beta, mu, eps, dt) {
    
    # First of all we extract from the variable net the number of agents in each compartment...
    x_i <- net$x
    r_i <- net$r
    n_i <- net$n
    s_i <- n_i - r_i - x_i
    N <- sum(n_i)
    
    # Then we define the time
    t_max <- 200
    times <- seq(0, t_max, dt)
    
    # The transition matrix 
    W_ij <- as_adjacency_matrix(net, attr = 'weight', sparse = F)
    R_ij <- W_ij / rowSums(W_ij)
    R_ij[is.nan(R_ij)] <- 0
    
    
    yini <- c(s_i, x_i, r_i)
    params <- list(beta, mu, eps, R_ij, dt)
    stationarity_indicator <- 10^100
    R <- X <- numeric(0)
    while ( stationarity_indicator > 1e-4 ) {
        out <- ode(yini, times, func = SIR_ode.iter, params, method = 'iteration')
        out <- as.matrix(out)
        # rowSums(out[,seq(2,3*n+1)])
        
        R <- c(R, rowSums(out[,seq(2*n+2,3*n+1)]) )
        X <- c(X, rowSums(out[,seq(n+2,2*n+1)]) )
        
        stationarity_indicator <- abs(R[length(R)] - R[length(R)-199]) / N
        yini <- as.vector( out[dim(out)[1], seq(2, 3*n+1)] )
        # print(stationarity_indicator)
    }
    
    return( list( X, R) )
}


#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------


