suppressPackageStartupMessages({
    library(igraph)
    library(compiler)
    library(matrixStats)
    library(ggplot2)
    library(dqrng)
    library(Matrix)
    library(latex2exp)
    library(viridis)
    library(tidyr)
    library(ggpubr)
    library(dplyr)
})

setwd("C:/Users/pierg/Documents/R/Epidemics-and-Mobility")
source("reaction-diffusion.R")
source("sample_sd.R")

# HERE STARTS THE REAL PROGRAM ##########################################################
{
    eps <- 10^0
    n <- 100L
    N <- n * 2000L
    dt <- 0.1
    t_max <- 50
    times <- seq(0, t_max, dt)
    
    # Weighted and directed network. Stop if not connected
    k <- 20
    # net <- make_star(n, mode = 'undirected', center = 1)
    net <- make_full_graph(n, directed = FALSE, loops = FALSE)
    # net <- sample_gnm(n, ceiling(2 * log(n)/n * choose(n,2)), directed = F )
    # net <- sample_pa(n, directed = F)
    if ( !is.connected(net) ) {print('DISCONNESSO'); stop()}
    # E(net)$w <- sample( 50 , size = gsize(net), replace = T)
    E(net)$w <- sample( 1 , size = gsize(net), replace = T)
    
    W_ij <- as_adjacency_matrix(net, attr = 'w', sparse = T)
    R_ij <- W_ij / rowSums(W_ij)
    
    # Number of infected, recovered and agents in general per node in the initial condition
    X_start <- n 
    R_start <- 0
    net$x <- sample_sd( n, X_start, 0 )
    net$r <- rep( R_start / n , n )
    net$n <- rep( N / n , n)
    
    
}

table <- read.table("results")
results <- list()
results[[1]] <- as.matrix( table[1:5] )
results[[2]] <- as.matrix( table[6:10] )


# ----------------------------- HEATMAP FOR METACRITICAL POINT -------------------------------------
n <- 500
N <- n * 2000L

mu <- 2
R0 <- numeric(0)

net <- sample_gnm(n, ceiling(3 * log(n)/n * choose(n,2)), directed = F )
if ( !is.connected(net) ) {
    net <- sample_gnm(n, ceiling(3 * log(n)/n * choose(n,2)), directed = F )
}
E(net)$w <- sample( 50 , size = gsize(net), replace = T)

W_ij <- as_adjacency_matrix(net, attr = 'w', sparse = T)
R_ij <- W_ij / rowSums(W_ij)

# Number of infected, recovered and agents in general per node in the initial condition
X_start <- n 
R_start <- 0
net$x <- sample_sd( n, X_start, 0 )
net$r <- rep( R_start / n , n )
net$n <- rep( N / n , n)

b <- seq(0.6, 1.4, length.out = 80)
b_max <- 6*mu
beta <-  c(b_max,  rnorm(n-1, mean = 3*mu , 0.5) )
while ( any(beta < 0 | beta > b_max) ) {
    idx <- beta < 0 | beta > b_max
    beta[idx] <- rnorm(sum(idx) , mean = 3*mu , 0.5)
}
R0_min <- as.numeric( beta %*% degree(net) / 2 / ecount(net) / mu ) 
beta <- beta / R0_min
# epsilons <- 10^seq(-0.5, -2, length.out = 30)
epsilons <- 10^seq(-0.5, -2, length.out = 80)

df <- numeric(0)

for (bb in b) {
    betas <- beta * bb
    betas[betas == max(betas)] <- b_max
        
    for (eps in epsilons) {
    
        out <- SIR_ode2(net, betas, mu, eps, dt = min(1,eps)) 
        # R_ode <- out[[2]]
        df <- rbind(df, c(max(out[[2]]), eps, bb) )  
        
    cat(which(bb == b), "\t\t", which(eps==epsilons), "\n")
    }
}

df <- as.data.frame(df)
colnames(df) <- c("recovered", "eps", "R0_min")
df$eps <- log10(df$eps)

## Troviamo al volo il valore critico di epsilon:
# Prima troviamo R0(eps)

R0s <- numeric(0)
min <- 0.6
epsilons <- 10^seq(-2, 2, length.out =40)
for(eps in epsilons){
    # betas <- beta * min
    # betas[max(beta) == beta] <- b_max
    
    id <- diag(1, n)
    
    V <- eps * mu * id + id - R_ij
    V_inv <- solve(V)
    K <- eps * diag(betas) %*% V_inv
    
    R0s <- c(R0s , max(Re(eigen( K )$values)) )
    print(which(eps==epsilons))
}
plot(epsilons, R0s, log = 'x')

# Poi per trovare eps* fittiamo con la funzione che sappiamo noi

fit.p <- seq(0.1, 5, 0.05)
fit.q <- seq_along(epsilons)
S <- 10^100
model_best <- p_best <- q_best <- 0
for (p in fit.p) {
    for( q in fit.q){
        model <- c( rep(min, q) , 
                    ( b_max*epsilons[-(1:q)] / (p+mu*epsilons[-(1:q)]) ) * (b_max/mu - min) / (b_max/mu) + min )
        S_new <- sum ((R0s - model)^2/model)
        if (S_new < S){
            S <- S_new
            model_best <- model
            p_best <- p
            q_best <- q
        }
    }
}
print(p_best)
print(q_best)
log10(epsilons[q_best])


df_heatmap <- read.csv("C:/Users/pierg/Documents/R/Epidemics-and-Mobility/figure3/heatmap_metapunto3.csv")
png("prova.png", width = 1147*4, height = 900*4)
ggplot(df_heatmap, aes(x = eps, y = R0_min, fill = recovered/N )) + geom_tile() + 
    xlab(TeX("$\\log_{10}(\\epsilon)$")) + ylab(TeX("$R_0^{min}$")) +
    scale_fill_gradient2(low=alpha("black", 1), mid=alpha("red", 0.8), high = alpha("yellow",0.8) , midpoint = 0.3) +
    guides(fill = guide_colorbar( title.position = "top",
                               title = "Attack rate") ) +
    theme(legend.key.size = unit(5, "cm"),
          legend.text = element_text(size = 120)) +
    geom_point(aes(x = -1, y = 1), col = "blue", size =30, pch = 17, show.legend = F) +
    # annotate("text", x = 1.8, y = -1.8, label = "sd") +
    scale_y_continuous(breaks= c(0.6, 1, 1.4)) + 
    scale_x_continuous(breaks= c(-2, -1, -0.5), labels = c("-2", TeX("$\\epsilon{*}$"), "-0.5") ) + 
    geom_segment(aes(y = 1, x = -2, yend = 1, xend = log10(1/(b_max-mu))), col = "blue", size = 8, lty = "dotted") +
    geom_segment(aes(y = min(b), x = log10(1/(b_max-mu)), yend = 1, xend = log10(1/(b_max-mu))), col = "darkviolet", size = 8, lty = "dotted") +
    theme(text = element_text(size = 20),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_blank(),
          panel.background = element_blank(),
          axis.text.y = element_text(colour=c("black", "blue", "black"), size=130, margin = margin(r=-3, unit = "cm")),
          axis.text.x = element_text(colour=c("black", "darkviolet", "black"), size=130, margin = margin(t=-3, unit = "cm")),
          axis.title.x = element_text(size=130),
          axis.title.y = element_text(size=130),
          # axis.ticks.length.x = unit(-3, "cm"),
          # axis.ticks.length.y = unit(-3, "cm"),
          # axis.ticks.x = element_line(size=6, colour=c("black", "darkviolet", "black")),
          # axis.ticks.y = element_line(size=6, colour = c("black", "blue", "black")),
          legend.title = element_text( angle = 0, size=130),
          legend.position = "right") 
     # scale_fill_viridis_c(option = "B")
dev.off()
d


## Adesso occupiamoci della seconda figura

df1 <- epsilons
epsilons <- 10^seq(-3, -0.5, length.out =50)
b <- seq(0.2)

for (bb in b) {
    betas <- beta * bb
    betas[betas == max(betas)] <- b_max
    
    R0s <- numeric(0)
    for (eps in epsilons) {
        
        id <- diag(1, n)
        
        V <- eps * mu * id + id - R_ij
        V_inv <- solve(V)
        K <- eps * diag(betas) %*% V_inv
        
        R0s <- c( R0s, max(Re(eigen( K )$values)))
        cat(which(bb == b), "\t\t", which(eps==epsilons), "\n")
    }
    df1 <- cbind(df1, R0s)
}

colnames(df1) <- c("R0", "eps", "R0_min")

ggplot(df1, aes(x = eps ))  + scale_x_log10() + geom_line( aes(y = R0s), color = "red", size = 2 ) + ylim(c(0,2.5)) +
    geom_line( aes(y = R0s.1), color = "blue", size = 2 ) + 
    geom_line( aes(y = R0s.2), color = "black", size = 2 ) +
    geom_segment(aes(x = 0.001, y = 1, xend = 10^-0.5, yend = 1), col = "black", size = 1, lty = "dashed") + 
    geom_point(aes(x = 0.1, y = 1), col = "red", size =4, pch = 17) +
    ylab(TeX("$R_0$")) + xlab(TeX("$\\epsilon$")) +
    theme(text = element_text(size = 20),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_blank(),
          panel.background = element_blank())
        
    
    scale_fill_gradient2(low="black", mid="red", high = "yellow", midpoint = 0.3) +
    theme(text = element_text(size = 20),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(colour=c("black", "blue", "black")), 
          axis.text.y = element_text(colour=c("black", "red", "black")),
          legend.title = element_text( angle = -90)) + 
    guides(fill = guide_legend(title.position = "right",
                               title = "Recovered")) +
    scale_x_continuous(breaks= c(0.6, 1, 1.4)) + 
    scale_y_continuous(breaks= round(c(-2, log10(epsilons[q_best]), -0.5), digits = 1)) + 
    geom_segment(aes(x = 1, y = min(log10(epsilons)), xend = 1, yend = log10(epsilons[q_best])), col = "blue", size = 2, lty = "dotted") +
    geom_segment(aes(x = min(b), y = log10(epsilons[q_best]), xend = 1, yend = log10(epsilons[q_best])), col = "red", size = 2, lty = "dotted")


    
    
    
    
    
    
    
    
    
# ----------------------------- ESTIMATE R0 USING MAX EIGENVALUE ----------------------------------------------

    # Controlla che gli input siano gli stessi della heatmap
    
n <- 500
net <- sample_gnm(n, ceiling(3 * log(n)/n * choose(n,2)), directed = F )
if ( !is.connected(net) ) {
    net <- sample_gnm(n, ceiling(3 * log(n)/n * choose(n,2)), directed = F )
}
# net <- sample_pa(n, m = 1, directed = F)
    
E(net)$w <- sample( 50 , size = gsize(net), replace = T)


W_ij <- as_adjacency_matrix(net, attr = 'w', sparse = T)
R_ij <- W_ij / rowSums(W_ij)

mu = 2
b_max <- 2*mu
beta <-  c(b_max,  rnorm(n-1, mean = mu , 0.2) )
# while ( any(beta < 0 | beta > b_max) ) {
#     idx <- beta < 0 | beta > b_max
#     beta[idx] <- rnorm(sum(idx) , mean = 3*mu , 0.5)
# }
# R0_min <- as.numeric(beta %*% strength(net, weights = E(net)$w) / mu / 2 / sum(E(net)$w) )
#     
# beta <- beta / R0_min
    
results <- list ( matrix(nrow = length(epsilons), ncol = 5) , matrix(nrow = length(epsilons), ncol = 5) )

epsilons <- 10^seq(-3, 3, length.out = 100)
N_iter <- 1
fractions <- 1

for (fraction in fractions){
R0s <- numeric(0)
a <- numeric(0)
for (eps in epsilons) {
    
    R0  <- summation <- 0
    for (iter in seq_len(N_iter)) {
        
        min <- 10^100
        betas <- beta
        while (min > fraction*1.05 | min < fraction*0.95){
            min <- as.numeric(betas %*% strength(net, weights = E(net)$w) / mu / 2 / sum(E(net)$w) )    
            betas <- betas / min * fraction
            betas[betas > b_max] <- b_max
        }    
        betas[max(betas)== betas] <- b_max
        
        id <- diag(1, n)
        
        V <- eps * mu * id + id - R_ij
        V_inv <- solve(V)
        K <- eps * diag(betas) %*% V_inv
        A_ij <- R_ij %*% solve( id*(1 + mu * eps) + R_ij )
        
        # R0 <- R0 + max(Re(eigen( K )$values))
        # R0 <- R0 + max(Re(eigen( K )$values)) 
        idx <- which(max(betas)==betas)[1]
        summation <- summation + sum( betas[-idx] *A_ij[idx,-idx] * A_ij[-idx, idx] / (b_max-betas[-idx]) )
    }
    
    # R0s <- append(R0s, R0 / N_iter)
    a <- append( a, summation / N_iter)
    print(which(eps==epsilons))
}

}


#### SINGOLA LINEA -------

first <- b_max * epsilons / ((1-R_ij[1,1]) + mu * epsilons ) 
second <- b_max *mu * epsilons^2 / ( mu*epsilons + 1 )^2 * a
model <- first + second
df <- as.data.frame( cbind( epsilons, model, R0s) ) 
colnames(df) <- c("eps", "model", "R0")

# Stimiamo al volo l'epsilon critico
mu<-2
b_max<-2*mu
fit.p <- seq(0.01, 10, 0.01)
fit.q <- seq_along(epsilons)
fit.r <- seq(0.01, min*mu, 0.005)
S <- 10^100
model_best <- p_best <- q_best <- r_best <- 0
for( q in fit.q){
    
    for (r in fit.r) {
        eps <-epsilons[q+1]
        p <- eps*(b_max - min*mu ) / (1-r/mu)
        model <- c( rep(min, q) , 
                    ( (b_max-r)*epsilons[-(1:q)] / (p+mu*epsilons[-(1:q)]) + r/mu ) ) # * (b_max/mu - min) / (b_max/mu) + min )
        S_new <- sum ((R0s - model)^2/R0s)
        if (S_new < S){
            S <- S_new
            model_best <- model
            p_best <- p
            q_best <- q
            r_best <- r
        
        }
    }
}
plot(epsilons, R0s, log="x")
lines(epsilons, model_best)
eps_crit1 <- epsilons[q_best+1]
print(r_best)
print(q_best)
# df1 <- gather(df, "type", "R0s", -eps)
df2 <- df

png("BA.png", height=346*4, width=712*4 )
# par(mar=c(0,0,0,0))
ggplot(data = df2, aes(x = eps, y =model) ) + #scale_x_log10() +# ylim(1, 2) +
    geom_line(linetype = "dashed", col = alpha('red', 1), size = 5.5) +
    geom_line(aes(x=eps, y=R0 ), size = 6, col = alpha('black', 0.5), )  + 
    geom_vline(xintercept =  eps_crit2, linetype = "dashed", lwd = 5, col = "darkviolet") +
    geom_hline(yintercept = 1, linetype = "dotted", lwd = 5) +
    geom_hline(yintercept = 2, linetype = "dotted", lwd = 5) +
    labs( x = TeX("$\\epsilon$"),
          y = TeX("$R_0$") ) +
          # title = "Barabasi-Albert Network" ) + 
    theme_classic() + 
    scale_y_continuous( breaks= c(1, 2), labels = c(TeX("$R_0^{min}$"), TeX("$\\frac{\\beta_{max}}{\\mu}$")) ) +
    coord_cartesian(ylim=c(1,2.1))+
    scale_x_continuous(breaks= c(0.01, eps_crit2, 1, 100), labels = c(0.01, TeX("$\\epsilon{*}$"), 1, 100), trans = "log10" ) +
    theme(  text = element_text(size = 35),
            legend.background = element_rect(fill = "transparent"),
            legend.key = element_rect(fill = "white", colour = "black"),
            legend.title = element_text(color = 'black', size = 20), 
            legend.text = element_text(color = "black", size = 35),
            legend.position = c(5/6, 0.5),
            legend.key.size = unit(1.5, "cm"),
            legend.key.width = unit(1.5,"cm"),
            axis.text.x = element_text(colour=c("black", "darkviolet", "black", "black"), size =95, margin = margin(t = 1, unit = "cm")),
            axis.text.y = element_text(size = 95, margin = margin(r = 1, unit = "cm"), colour = "black"),
            axis.title.x = element_text(size=95, colour = "black"),
            axis.title.y = element_text(size=95, colour = "black"),
            axis.ticks.x = element_line(size = 4, colour = c("black", "darkviolet", "black", "black") ),
            axis.ticks.y = element_line(size = 4, colour = "black"),
            axis.ticks.length = unit(0.7,"cm"),
            axis.line = element_line(size=1)) 
    # annotate("text", x = 1, y = 1.1, label = "P.T.", size = 6, col = 'red') 
    # annotate("text", x = 0.01, y = 1.1, label = TeX("spectral radius of \\textbf{K}"), size = 6, col = 'black') +
    # annotate("text", x = 100, y = 1.5, label = TeX("(a)"), size = 10, col = 'black')
    # scale_x_continuous(breaks= c(0.6, 1, 1.4))
dev.off()
 
d2

ggarrange(d1,d2, nrow = 2, ncol = 1)

#### MULTIPLE LINEE ---------


df <- as.data.frame( cbind( rep(epsilons, length(fractions)) , as.vector(results[[1]]), rep(b_max*fractions/mu, each = length(epsilons))))
colnames(df) <- c("eps", "R0", "sh")
df$model <- first

plot(epsilons, model, 'l', lty = 2, col='red', log='x', lwd = 2, xlab = TeX("$\\epsilon$"))
points(epsilons, R0s)

color = matrix('black', nrow = dim(results[[1]])[1], ncol = dim(results[[1]])[2])
color[abs(results[[1]] - first ) < 0.05] <- 'red'
col <- as.vector(color)

d <- ggplot(data = df, aes(x=eps, y=R0, shape = factor(sh) ) ) + scale_x_log10() + ylim(0, 2) + 
    geom_point( size = 3, color = col ) +
    geom_line(aes(x = eps, y =model), linetype = "dotted", col = 'red', size = 1) + 
    xlab(TeX("$\\epsilon$")) + ylab(TeX("$R_0$")) +
    labs( shape = TeX(" $\\bar{ R }_0 $") )
d + theme(  text = element_text(size = 20),
            plot.background=element_rect(fill="transparent",colour=NA),
            legend.background = element_rect(fill = "transparent"),
            legend.key = element_rect(fill = "white", colour = "black"),
            legend.title = element_text(color = 'black', size = 20), 
            legend.text = element_text(color = "black", size = 15),
            legend.position = c(5/6, 0.5),
            legend.key.size = unit(1.5, "cm"),
            legend.key.width = unit(1.5,"cm") )+

plot(epsilons, R0s, ylim = c(0,b_max/mu), lty = 1, log='x' ,
     col = color  )

# lines(epsilons, first, col = 'red', lwd = 0.5 )
lines(epsilons, first+second, col = 'red', lty = 4 )
abline(h = min )
abline(v=1, lt = 2)

title(TeX('Barabasi-Albert weighted, n = 1000, $\\beta\\in (\\bar{\\beta}-r,\\bar{\\beta}+r)$'))
text(3e-3, 2,TeX('$\\beta_{medio} = \\dfrac{3}{10}\\beta_{max}$'))
legend <- c('r = 0.05', 'r = 0.1', 'r = 0.2', 'r = 0.3', 'r = 0.4', 'r=0.5', 'r=1', 'r=2')
legend(10, 1, legend, pch = c(5:1,6,7,8))

# ------------------ FIT USING THE SIGMOID-------------------------------

fit.p <- seq(0.1, 5, 0.05)
fit.q <- seq_along(R0s)
S <- 10^100
model_best <- p_best <- q_best <- 0
for (p in fit.p) {
    for( q in fit.q){
        model <- c( rep(min, q) , 
                    ( b_max*epsilons[-(1:q)] / (p+mu*epsilons[-(1:q)]) ) * (b_max/mu - min) / (b_max/mu) + min )
        S_new <- sum ((R0s - model)^2/model)
        if (S_new < S){
            S <- S_new
            model_best <- model
            p_best <- p
            q_best <- q
        }
    }
}
print(p_best)
print(q_best)
color = rep_len('black', length(epsilons))
color[abs(R0s - model_best ) < 0.05] <- 'red'
plot(epsilons, R0s, ylim = c(0,b_max/mu), log ='x', col = color )
title('Barabasi-Albert')
text(3e-3, 2,TeX('$\\beta_{medio} = \\beta_{max}/2$'))
lines(epsilons, model_best)
abline(v = 0.05*p_best / (mu - 0.05*mu) )


n <- 1000
n_i <- degree( sample_pa(n, m = 10, directed = F) )
b<-0
for (i in seq_along(n_i)) {
    
    b[i] <- sum(degree(sample_pa(n_i[i], m = 10, directed = F)))
}
plot(n_i,b,log='xy')
S <- 10^100
for (j in seq(1, 2, 0.0001) ){
    model <- n_i^j
    model <- model + min(b) - min(model)
    S_new <- sum((model - b)^2 )
    if(S_new < S){
        S<-S_new
        model_best <- model
        exp <- j
    }
}
lines(n_i, model_best)
print(exp)



# ----------------------------- FIND THE CRITICAL EPSILONS ----------------------------------------------

n <- 1000
epsilons <- 10^seq(-2, 3, length.out = 20)
# fraction <- seq(0, 1, 0.1)
fraction <- 0.1
betas <- result <- numeric()
R0s <- mean_betas <- numeric()
mu <- 2
b_max <- 4
N_iter <- 1

for (eps in epsilons) {
    
    R0 <- 0
    a <- min <- eps_critical <- 0
    for (iter in seq_len(N_iter)) {
        
        # Weighted and directed network. Choose the one you want
        # Tree
        # net <-
        #     make_tree(n, children = 25, mode = 'undirected') %>%
        #     add.edges(., c(2,3, 3,4, 2,4))
        # net <-
        #     make_tree(n, children = 5, mode = 'undirected')
        # Star-like
        # net <- make_star(n, mode = 'undirected', center = 1)
        # Completo
        # net <- make_full_graph(n, directed = FALSE, loops = FALSE)
        # Erdos-Renyi
        p <- 0.5
        # net <- sample_gnp(n, ceiling(3 * log(n)/n * choose(n,2)), directed = F )
        net <- sample_gnp(n, p, directed = T )
        if ( !is.connected(net) ) {
            net <- sample_gnp(n, p, directed = T )
        }
        # Barabasi-Albert
        # net <- sample_pa(n, m = 2, directed = F)
        # Small-world
        # net <- sample_smallworld(1, n, 4, 0.1)
        # Stochastic Block Model
        # pref.matrix <- matrix(runif(n),10,10)
        # net <- sample_sbm(n, (pref.matrix + t(pref.matrix))/2 , block.sizes = rep(n/10, 10))
        # Double star
        # net <-
        #     make_star(n-5, mode = 'undirected', center = 1) %du% make_star(5, mode = 'undirected', center = 1) %>%
        #     add_edges(., c(1,n-5 + 1))
        
        # # E(net)$w <- sample( 50 , size = gsize(net), replace = T)
        E(net)$w <- sample( 50 , size = gsize(net), replace = T)
        # 
        W_ij <- as_adjacency_matrix(net, attr = 'w', sparse = T)
        # # diag(W_ij) <- 100
        R_ij <- W_ij / rowSums(W_ij)
        
        beta <- sample_beta(n, b_max, fraction)
        # beta <- c( b_max, rep(b_max*fraction, n-1) )
        # beta <- c(b_max, runif(n-1, b_max*fraction-2, b_max*fraction+2) )
        
        id <- diag(1, n)
        
        V <- eps * mu * id + id - R_ij
        V_inv <- solve(V)
        K <- eps * diag(beta) %*% V_inv
        
        R0 <- R0 + max(Re(eigen( K )$values))
        a <- a + sum(beta[-1]*R_ij[1,-1]*R_ij[-1,1]/(b_max-beta[-1]))
        # a <- a + sum(b*R_ij[1,-1]*R_ij[-1,1]/(b_max-b))
        min <- min + sum(degree(net)*beta)/mu/2/gsize(net)
    }
    
    min <- min / N_iter
    print(which(eps==epsilons))
    
    if (R0 / N_iter - min < 0.1){
        eps_critical <- eps
    }
}

par(mfrow=c(1,2))
a <- sample_gnp(20, 0.15)
b <- sample_pa(20, m=1, directed = F)
plot(a, vertex.size = 15, vertex.color = 'red',
     vertex.label.color='white',
     vertex.label.font = 2,
     # layout=layout_on_sphere(a),
     main = "a) Erdos-Renyi network")
plot(b, vertex.size = 15, vertex.color = 'red',
     vertex.label.color='white',
     vertex.label.font = 2,
     main = "b) Barabasi-Albert network")

