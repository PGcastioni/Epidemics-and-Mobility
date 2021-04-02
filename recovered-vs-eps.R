suppressPackageStartupMessages({
    library(ggplot2)
    library(gganimate)
    library(gifski)
    library(dplyr)
    library(Matrix)
    library(igraph)
    library(RColorBrewer)
    library(reshape2)
    library(viridis)
    library(mapproj)
    library(latex2exp)
    library(tidyr)
    library(ggpubr)
    library(DescTools)
})

setwd("C:/Users/pierg/Documents/R/Epidemics-and-Mobility")
# setwd("C:/Users/pierg/Documents/Fisica Magistrale/Tesi Magistrale/REAL WORLD DATA")

source("reaction-diffusion.R")
source("sample_sd.R")
source("OD_matrices.R")


het_parameter <- function(x) {
    a <- ( quantile(strength(x))[4] - quantile(strength(x))[2] ) / max(strength(x))
    return( a )
}

prefix <- "FRA"
top_percentage <- 0
out <- mobility_matrix(prefix, top_percentage)
vcount(out[[1]])
net <- out[[1]]
reciprocity(net)
nodes <- out[[2]]
sum(nodes$pop)
edges <- out[[3]]
# hist(strength(net))
# hist(nodes$pop[-which(max(nodes$pop)==nodes$pop)], xlim = c(0, 4e5), breaks = 100)
# net <- as.undirected(net, mode = "each")

# Initialize the network and its initial conditions
{
    # net <- simplify(net)
    n <- vcount(net)
    if (is.null(nodes$pop)) {
        net$n <- rep(1e4, n)
    } else {
        net$n <- round(nodes$pop)
    }
    # net$n <- sample(seq(100, 1e4,20), n, replace = T)
    
    W_ij <- as_adjacency_matrix(net, attr = 'weight', sparse = F)
    R_ij <- W_ij / rowSums(W_ij)
    
    # net$x <- rep( 0, n )
    # net$x[which( betas == max(betas))] <- 20
    # net$x <- sample_sd( n , X_start, target_sd)
    # net$x <- sample(c(0,10), n, replace = T)
    # net$x <- rep(1000, n)
    net$x <- rep(10, n) #????
    
    net$r <- rep(0, n)
    N <- sum(net$n)
    
    mu <- 0.1
    # betas <- runif(n) * mu
    # betas[betas == max(betas)] <- 3 * mu
    # betas <- 1.5*mu
    if (is.null(nodes$pop)){
        betas <- rlnorm(n, meanlog = 2, sdlog = 2)
    } else {
        betas <- net$n
    }
    
    target <- 0.5
    min <- 10^100
    while (min > target*1.05 | min < target*0.95){
        betas <- betas / min * target
        betas[betas > 2 * mu] <- 2* mu
        betas[max(betas)== betas] <- 2*mu
        min <- as.numeric(betas %*% strength(net, weights = E(net)$weight) / mu / 2 / sum(E(net)$weight) ) 
    }    
    
    
    
    # as.numeric(betas %*% degree(net) / mu / 2 / ecount(net) )    

}

# Attenzione, rinormalizziamo la popolazione
# net$n <- round(net$n / 100)
# N <- sum(net$n)


#### ORA SI CALCOLA L'ANDAMENTO DEI RECOVERED AL VARIARE DI EPSILON

epsilons <- 10^seq(0 , 4 , 0.1)
# epsilons <- 10000
R_max <- numeric (length(epsilons))
df <- numeric(0)
for (eps in epsilons) {
    out <- SIR_ode2(net, betas, mu, eps, dt = min(1,eps)/10) 
    R_max[eps == epsilons] <- max(out[[2]])/N
    print(which(eps == epsilons))
    
    X_ode <- out[[1]]
    R_ode <- out[[2]]
    t_max <- length(X_ode)
    # plot( c(0, t_max) , c(1,1), 'l', lty='dashed',
    #       xlab = 'time', ylab = 'x(t)', xlim = c(0,t_max), ylim = c(-0, 1.1) )
    # lines(seq(1, t_max), 1 - X_ode/N - R_ode/N, col ='black', lty =2, lwd = 3)
    # lines(seq(1, t_max), X_ode/N, 'l', col='red', lty =2, lwd = 3)
    # lines(seq(1, t_max), R_ode/N, col = 'blue', lty =2, lwd = 3)
    
    idx <- max(X_ode)==X_ode
    idx2 <- min( 2*which(idx) , t_max )
    idx5 <- min( 5*which(idx) , t_max )
    idx10 <- min( 10*which(idx) , t_max )
    
    
    df <- rbind(df, c(eps, max(R_ode)/N, 
                      R_ode[idx10]/N,
                      R_ode[idx5]/N,
                      R_ode[idx2]/N,
                      R_ode[idx]/N ) )
                      # R_ode[round(t_max/2)]/N, 
                      # R_ode[round(t_max/4)]/N, 
                      # R_ode[round(t_max/10)]/N,
                      # R_ode[round(t_max/20)]/N ) )
                      
                      
    
    
}

df <- as.data.frame(df)
colnames(df) <- c("eps", "1", "2", "3", "4", "5")
df6 <- gather(df, "stage", "recovered", -eps)
pal <- brewer.pal(7 , "YlOrRd")[ 7:3 ]

stages <- c("long-time limit", "10x infectious peak", "5x infectious peak", "2x infectious peak", "infectious peak")
f6 <- ggplot(df6, aes(x=eps, y = recovered, color = as.factor(stage))) + scale_x_log10(  ) +
    labs(x = TeX("$\\epsilon$"), y = "Recovered", color = "Epidemic stage", title = "Spain") +
    geom_line(size = 2) + 
    scale_color_manual(values = pal, labels = stages) +
    geom_segment(aes(x = min(eps), y = 0, xend = max(eps), yend = 0), lty = "dashed", col = "black") +
    theme(text = element_text(size = 15),
          panel.grid.major = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.key = element_rect(fill = "white", colour = "black"),
          legend.title = element_text(color = 'black', size = 15), 
          legend.text = element_text(color = "black", size = 13),
          legend.position = c(0.2, 0.85),
          legend.key.size = unit(1, "cm"),
          legend.key.width = unit(0.8,"cm"),
          plot.title = element_text(hjust = 0.5))
f6

ggarrange(f1, f2, f3, f4, f5, f6, ncol = 3, nrow = 2, legend = 'bottom', common.legend = TRUE)

sum(net$n[betas>mu])/N

plot(epsilons, rep(0, length(epsilons)), 'l', lty = 'dashed', ylim=c(0, max(df[5,])/N), log='x')
lines(epsilons, df[5,]/N)
lines(epsilons, df[4,]/N)
lines(epsilons, df[3,]/N)
lines(epsilons, df[2,]/N)
lines(epsilons, df[6,]/N)



plot(epsilons, R_max, 'l',log='x')
points(epsilons, R_max, log='x', pch = 1)
df <- as.data.frame(cbind(aa, bb, col = cc))
aa <- c(aa, R_max)
bb <- c(bb, epsilons)
cc <- c(cc, rep(5, length(R_max)))
d <- ggplot(data = df) + scale_x_log10() + xlab(TeX("$\\epsilon$")) + ylab(TeX("$R_{max}$")) +
    geom_point( aes(x = bb, y = aa,  shape = factor(col)), size = 3) +
    labs(shape = "Dataset") +
    theme(  text = element_text(size = 20),
            plot.background=element_rect(fill="transparent",colour=NA),
            legend.background = element_rect(fill = "transparent"),
            legend.key = element_rect(fill = "white", colour = "black"),
            legend.title = element_text("dasdsa<", color = 'black', size = 20), 
            legend.text = element_text(color = "black", size = 15) )



######################################################################################################################################
######################################################################################################################################
######################################################################################################################################


# FIRST WE TURN ALL THE DATA INTO MOBILITY MATRICES

ISOs <- unique( substring( list.files("C:/Users/pierg/Documents/Fisica Magistrale/Tesi Magistrale/REAL WORLD DATA/Mobility_Flows/data"), 1, 3 ))
country_names <- read.csv("https://raw.githubusercontent.com/lukes/ISO-3166-Countries-with-Regional-Codes/master/all/all.csv") %>% filter(alpha.3 %in% ISOs) %>% select(name)
continent <- read.csv("https://raw.githubusercontent.com/lukes/ISO-3166-Countries-with-Regional-Codes/master/all/all.csv") %>% filter(alpha.3 %in% ISOs) %>% select(region)
    
ISOs <- c("FRA", "ESP", "PRT", "ITA", "SEN", "CIV", "NGA", "ZMB", "CHN", "IND", "AFG", "AZE", "TUR", "MEX", "BRA", "PHL")
top_percentage <- 0
G <- IQR <- numeric()
epsilons <- 10^seq(0 , 4 , 0.1)
nets <- nodes <- edges <- list()
    # plots <- recovered <- list()
for ( iso in ISOs ) {
    
    out <- mobility_matrix(iso, top_percentage)
    idx <- which(iso == ISOs)
    nets[[idx]] <- out[[1]]
    nodes[[idx]] <- out[[2]]
    edges[[idx]] <- out[[3]]
    # net <- as.undirected(net, mode = "each")
    IQR[idx] <- het_parameter(nets[[idx]])
    G[idx] <- Gini(strength(nets[[idx]]))
    
}




######################################################################################################################################
######################################################################################################################################
######################################################################################################################################



# FIRST WE TURN ALL THE DATA INTO MOBILITY MATRICES

# ISOs <- unique( substring( list.files("C:/Users/pierg/Documents/Fisica Magistrale/Tesi Magistrale/REAL WORLD DATA/Mobility_Flows/data"), 1, 3 ))
# ISOs <- c("FRA", "ESP", "PRT", "ITA", "SEN", "CIV", "NGA", "ZMB", "CHN", "IND", "AZE", "MMR", "TUR", "MEX", "BRA", "PHL")
ISOs <- c("FRA", "ESP", "PRT", "ITA", "SEN", "CIV", "NGA", "ZMB", "CHN", "IND", "AZE", "MMR")
country_names <- read.csv("https://raw.githubusercontent.com/lukes/ISO-3166-Countries-with-Regional-Codes/master/all/all.csv") %>% filter(alpha.3 %in% ISOs) %>% select(name)
continent <- read.csv("https://raw.githubusercontent.com/lukes/ISO-3166-Countries-with-Regional-Codes/master/all/all.csv") %>% filter(alpha.3 %in% ISOs) %>% select(region)


top_percentage <- 0
epsilons <- 10^seq(0 , 4 , 0.1)
nets <- nodes <- edges <- list()
for ( iso in ISOs ) {
    
    out <- mobility_matrix(iso, top_percentage)
    idx <- which(iso == ISOs)
    nets[[idx]] <- out[[1]]
    nodes[[idx]] <- out[[2]]
    edges[[idx]] <- out[[3]]
    # net <- as.undirected(net, mode = "each")
    
    
}

# SOME CUSTOMARY MODIFICATION

# sapply(nodes, function(x){length(x[,1])} )
# populations <- sapply(nodes, function(x) sum(x$pop) )
# idx <- sapply(nets, function(x) vcount(x) ) %in% sort( sapply(nets, function(x) vcount(x) ), decreasing = T )[1:12]
# idx <- populations %in% sort(populations, decreasing = T)[1:16]
# idx <- idx & sapply(nets, function(x) reciprocity(x) ) > 0.9
# idx <- idx & sapply(nets, function(x) vcount(x) ) > 20
# nodes <- nodes[ idx ]
# edges <- edges[ idx ]
# nets <- nets[ idx ]
# ISOs <- ISOs[ idx ]
# country_names <- country_names[ idx, ]
# continent <- continent[ idx, ]
# names <- as.character(seq_len(sum(!idx)))


# recovered <- data <- plots <- list()

for ( iso in ISOs ) {
    
    net <- nets[[which(iso == ISOs)]]
    node <- nodes[[which(iso == ISOs)]]
    n <- vcount(net)
    net$n <- round(node$pop)
    net$x <- rep(10, n)
    net$r <- rep(0, n)
    mu <- 0.1
    N <- sum(net$n)
    
    
    betas <- net$n
    target <- 0.5
    min <- 10^100
    while (min > target | min < target*0.9){
        betas <- betas / min * target
        betas[betas > 2 * mu] <- 2* mu
        betas[max(betas)== betas] <- 2*mu
        min <- as.numeric(betas %*% strength(net, weights = E(net)$w) / mu / 2 / sum(E(net)$w) ) 
    }  
    
    # while (min > target | min < target*0.95){
    #     min <- as.numeric(betas %*% strength(net, weights = E(net)$w) / mu / 2 / sum(E(net)$w) )    
    #     betas <- betas / min * target
    #     betas[betas > 2 * mu] <- 2* mu
    # }    
    # betas[max(betas)== betas] <- 2*mu
    
    R_max <- numeric (length(epsilons))
    df <- numeric(0) 
    for (eps in epsilons) {
        out <- SIR_ode2(net, betas, mu, eps, dt = min(1,eps)/10) 
        R_max[eps == epsilons] <- max(out[[2]])/N
        cat(iso, "\t", which(eps == epsilons), "\n")
        
        X_ode <- out[[1]]
        R_ode <- out[[2]]
        t_max <- length(X_ode)
        # plot( c(0, t_max) , c(1,1), 'l', lty='dashed',
        #       xlab = 'time', ylab = 'x(t)', xlim = c(0,t_max), ylim = c(-0, 1.1) )
        # lines(seq(1, t_max), 1 - X_ode/N - R_ode/N, col ='black', lty =2, lwd = 3)
        # lines(seq(1, t_max), X_ode/N, 'l', col='red', lty =2, lwd = 3)
        # lines(seq(1, t_max), R_ode/N, col = 'blue', lty =2, lwd = 3)
        
        idx <- max(X_ode)==X_ode
        idx2 <- min( 2*which(idx) , t_max )
        idx5 <- min( 5*which(idx) , t_max )
        idx10 <- min( 10*which(idx) , t_max )
        
        
        df <- rbind(df, c(eps, max(R_ode)/N, 
                          R_ode[idx10]/N,
                          R_ode[idx5]/N,
                          R_ode[idx2]/N,
                          R_ode[idx]/N ) )
        
    }
    
    df <- as.data.frame(df)
    recovered[[which(iso== ISOs)]] <- df
    colnames(df) <- c("eps", "1", "2", "3", "4", "5")
    df6 <- gather(df, "stage", "recovered", -eps)
    pal <- brewer.pal(7 , "YlOrRd")[ 7:3 ]
    data[[which(iso==ISOs)]] <- df6
    names(data)[which(iso==ISOs)] <- iso
    
    stages <- c("long-time limit", "10x infectious peak", "5x infectious peak", "2x infectious peak", "infectious peak")
    f <- ggplot(df6, aes(x=eps, y = recovered, color = as.factor(stage))) + scale_x_log10(  ) +
        labs(x = TeX("$\\epsilon$"), y = "Attack rate", color = "Epidemic stage", title = iso) +
        geom_line(size = 1.5) + 
        scale_color_manual(values = pal, labels = stages) +
        geom_segment(aes(x = min(eps), y = 0, xend = max(eps), yend = 0), lty = "dashed", col = "black") +
        theme(text = element_text(size = 12),
              panel.grid.major = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(colour = "black"),
              legend.key = element_rect(fill = "white", colour = "black"),
              legend.title = element_text(color = 'black', size = 15), 
              legend.text = element_text(color = "black", size = 13),
              legend.position = c(0.2, 0.85),
              legend.key.size = unit(1, "cm"),
              legend.key.width = unit(0.8,"cm"),
              plot.title = element_text(hjust = 0.5))
    
    plots[[which(iso == ISOs)]] <- f
}


ggarrange(plotlist = plots, nrow=4, ncol=4, legend = "bottom", common.legend = T)

sapply(recovered, function(x) x$V1[x$V2 > 0.02][1] )
drop <- lapply( data, function(x) {x <- x$recovered[x$stage == 1]; epsilons[ which( diff(x) == min(diff(x) ) ) + 1] } )
ind <- sapply(nodes, function(x) max(x$pop))


## USA VERSION


setwd("C:/Users/pierg/Documents/Fisica Magistrale/Tesi Magistrale/REAL WORLD DATA/USA_ods")
out <- mobility_matrix("USA", top_percentage)
nodes <- out[[2]]
edges <- out[[3]]

# Each US state is
states <- c(12, 48, 6, 8, 17, 42)
states_name <- c("USA (Florida)", "USA (Texas)", "USA (California)", "USA (Colorado)", "USA (Illinois)", "USA (Pennsylvania)")

for (state in states) {
    node <- nodes[nodes$id > state*1000 & nodes$id < (state+1)*1000, ]
    edge <- edges[edges$from %in% node$id, ]
    edge <- edge[edge$to %in% node$id, ]
    net <- graph_from_data_frame(d = dplyr::select(edge, from, to, flow), vertices = node$id, directed = T)
    E(net)$weight <- edge$flow
    net <- as.undirected(net, mode = "each")
    
    n <- vcount(net)
    net$n <- round(node$population)
    net$x <- rep(10, n)
    net$r <- rep(0, n)
    mu <- 0.1
    N <- sum(net$n)
    
    
    betas <- net$n
    target <- 0.5
    min <- 10^100
    while (min > target | min < target*0.9){
        betas <- betas / min * target
        betas[betas > 2 * mu] <- 2* mu
        betas[max(betas)== betas] <- 2*mu
        min <- as.numeric(betas %*% strength(net, weights = E(net)$w) / mu / 2 / sum(E(net)$w) ) 
    }  
    
    # while (min > target | min < target*0.95){
    #     min <- as.numeric(betas %*% strength(net, weights = E(net)$w) / mu / 2 / sum(E(net)$w) )    
    #     betas <- betas / min * target
    #     betas[betas > 2 * mu] <- 2* mu
    # }    
    # betas[max(betas)== betas] <- 2*mu
    
    R_max <- numeric (length(epsilons))
    df <- numeric(0) 
    for (eps in epsilons) {
        out <- SIR_ode2(net, betas, mu, eps, dt = min(1,eps)/10) 
        R_max[eps == epsilons] <- max(out[[2]])/N
        cat(states_name[state==states], "\t", which(eps == epsilons), "\n")
        
        X_ode <- out[[1]]
        R_ode <- out[[2]]
        t_max <- length(X_ode)
        # plot( c(0, t_max) , c(1,1), 'l', lty='dashed',
        #       xlab = 'time', ylab = 'x(t)', xlim = c(0,t_max), ylim = c(-0, 1.1) )
        # lines(seq(1, t_max), 1 - X_ode/N - R_ode/N, col ='black', lty =2, lwd = 3)
        # lines(seq(1, t_max), X_ode/N, 'l', col='red', lty =2, lwd = 3)
        # lines(seq(1, t_max), R_ode/N, col = 'blue', lty =2, lwd = 3)
        
        idx <- max(X_ode)==X_ode
        idx2 <- min( 2*which(idx) , t_max )
        idx5 <- min( 5*which(idx) , t_max )
        idx10 <- min( 10*which(idx) , t_max )
        
        
        df <- rbind(df, c(eps, max(R_ode)/N, 
                          R_ode[idx10]/N,
                          R_ode[idx5]/N,
                          R_ode[idx2]/N,
                          R_ode[idx]/N ) )
        
    }
    
    df <- as.data.frame(df)
    # recovered[[which(iso== ISOs)]] <- df
    colnames(df) <- c("eps", "1", "2", "3", "4", "5")
    df6 <- gather(df, "stage", "recovered", -eps)
    pal <- brewer.pal(7 , "YlOrRd")[ 7:3 ]
    data[[which(state==states)]] <- df6
    names(data)[which(state==states)] <- states_name[state==states]
    
    stages <- c("long-time limit", "10x infectious peak", "5x infectious peak", "2x infectious peak", "infectious peak")
    f <- ggplot(df6, aes(x=eps, y = recovered, color = as.factor(stage))) + scale_x_log10(  ) +
        labs(x = TeX("$\\epsilon$"), y = "Attack rate", color = "Epidemic stage", title = "USA (Colorado)") +
        geom_line(size = 1.5) + 
        scale_color_manual(values = pal, labels = stages) +
        geom_segment(aes(x = min(eps), y = 0, xend = max(eps), yend = 0), lty = "dashed", col = "black") +
        theme(text = element_text(size = 12),
              panel.grid.major = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(colour = "black"),
              legend.key = element_rect(fill = "white", colour = "black"),
              legend.title = element_text(color = 'black', size = 15), 
              legend.text = element_text(color = "black", size = 13),
              legend.position = c(0.2, 0.85),
              legend.key.size = unit(1, "cm"),
              legend.key.width = unit(0.8,"cm"),
              plot.title = element_text(hjust = 0.5))
    
    plots[[which(state==states)]] <- f
}


ggarrange(plotlist = plots, nrow=4, ncol=4, legend = "bottom", common.legend = T)

sapply(recovered, function(x) x$V1[x$V2 > 0.02][1] )
drop <- lapply( data, function(x) {x <- x$recovered[x$stage == 1]; epsilons[ which( diff(x) == min(diff(x) ) ) + 1] } )
ind <- sapply(nodes, function(x) max(x$pop))

for(ii in n$id) {
    e$flow[e$from == ii] <- e$flow[e$from == ii]*n$pop[n$id == ii]
}

sapply(n$id, function(x) {
    e$flow[e$from == x] <- e$flow[e$from == x]*n$pop[n$id == x]
})


###########################################################################################################################
## IMAGES FROM IMPORTED DATA
###########################################################################################################################


# string <- TeX('$\\dfrac{\\langle s^2 \\rangle}{\\langle s \\rangle} =$')
plots = list()
setwd("C:/Users/pierg/Documents/Fisica Magistrale/Tesi Magistrale/REAL WORLD DATA/Mobility_Flows")
a <- read.csv("countries3.csv", sep = ",")

ISOs <- c("FRA", "ESP", "PRT", "ITA", "SEN", "CIV", "NGA", "ZMB", "CHN", "IND", "USA (Florida)", "USA (Colorado)", "TUR", "MEX", "BRA", "PHL")
for (iso in ISOs) {
    
    df <- a[,(1:6)+6*(which(iso==ISOs)-1)]
    colnames(df) <- c("eps", "1", "2", "3", "4", "5")
    df6 <- gather(df, "stage", "recovered", -eps)
    pal <- brewer.pal(7 , "YlOrRd")[ 7:3 ]
    string <- TeX( paste0("G = ", toString(round( G[which(iso==ISOs)] , digits = 2)) ) )
    
    stages <- c("long-time limit", "10x infectious peak", "5x infectious peak", "2x infectious peak", "infectious peak")
    f <- ggplot(df6, aes(x=eps, y = recovered, color = as.factor(stage))) + scale_x_log10(  ) +
        labs(x = TeX("$\\epsilon$"), y = "Attack rate", color = "Epidemic stage", title = iso) +
        geom_line(size = 6) +
        scale_color_manual(values = pal, labels = stages) +
        geom_vline(xintercept = 10, lty = "dashed", col = "black", size = 3.5) +
        # annotate("text", x = 2, y = max(df6$recovered)*4/5, label = expression(beta = round( het[which(iso==ISOs)]) ) , size = 20, col = "darkgrey"  ) +
        annotate("text", x = 2.5, y = max(df6$recovered)*4/5, label =  string, size = 20, col = "black" , parse = T ) +
        # scale_x_continuous(breaks= c(1, 10, 100, 1000)) + 
        # geom_segment(aes(x = min(eps), y = 0, xend = max(eps), yend = 0), lty = "dashed", col = "black") +
        theme(text = element_text(size = 12),
              panel.grid.major = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(colour = "black", size = 2),
              axis.text.x = element_text(size=70, margin = margin(t = 1.5, unit = "cm"), hjust = .65),
              axis.text.y = element_text(size=70, margin = margin(r = 1.5, unit = "cm")),
              axis.title.x = element_text(size=85),
              axis.title.y = element_text(size=80),
              axis.ticks.length.x = unit(-0.5, "cm"),
              axis.ticks.length.y = unit(-0.5, "cm"),
              axis.ticks.x = element_line(size=2),
              axis.ticks.y = element_line(size=2),
              legend.key = element_rect(fill = "white", colour = "black"),
              legend.title = element_text(color = 'black', size = 95),
              legend.text = element_text(color = "black", size = 80),
              legend.position = c(0.2, 0.85),
              legend.key.size = unit(4, "cm"),
              legend.key.width = unit(4,"cm"),
              plot.title = element_text(hjust = 0.5, size = 100) )
    
    plots[[which(iso == ISOs)]] <- f
}
png("prova.png", width = 1600*4, height = 775*4)
ggarrange(plotlist = plots, nrow=4, ncol=4, legend = "bottom", common.legend = T)
dev.off()

# het <- c(3.046942e+04, 3.658262e+04, 1.543556e+04, 282672.1, 1.505035e+05, 1.772538e+04, 2.178923e+05, 5.180110e+04, 1.446811e+06, 1.207787e+06, 1.791333e+06, 1.056991e+06, 3.984609e+05, 5.428241e+05, 9.463640e+05, 4.411050e+05)
# G <- c(0.4772961, 0.3942572, 0.6404427, 0.3841042, 0.4466481, 0.6020811, 0.2752073, 0.3671779, 0.2175603, 0.5140182, 0.8411913, 0.3095361, 0.9044770, 0.4092287, 0.5372721, 0.5455268 )

# write.csv(aa, "countries3.csv", col.names = F)
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################

