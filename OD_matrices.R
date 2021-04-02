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
})

alpha_distribution <- function(a, flex){
    if (flex>1 |flex<0) stop("flex must be a fraction")
    flex <- quantile(a, flex)
    return(tanh(0.1*(a-flex))*0.3 + 0.7)
}

# Possible prefixes: "ITA", "USA", "FRA", "ESP", "PRT", "SEN", "CIV"
mobility_matrix <- function(prefix, top_percentage){
    temporal_network <- F
    if ( prefix == "USA" | prefix == "ITA") temporal_network <- T # Solo IT e US sono network temporali
    
    if (prefix == "USA"){
        file <- paste0("C:/Users/pierg/Documents/Fisica Magistrale/Tesi Magistrale/REAL WORLD DATA/", prefix, "_ods")
    } else {
        file <- "C:/Users/pierg/Documents/Fisica Magistrale/Tesi Magistrale/REAL WORLD DATA/Mobility_Flows/data"
    }
    setwd(file)
    
    # Import nodes and edges and remove rows where flow is NA or zero
    nodes <- read.csv(paste0(prefix,"_nodes.csv"), header = T)
    edges <- read.csv(paste0(prefix,"_edges.csv"), header = T, fill = T )
    edges <- edges %>% filter(!is.na(flow) & flow != 0)
    edges <- edges %>% filter(from != to)
    # If you don't have a temporal network is like the time is constant
    if (!temporal_network) time <- 1
    
    if (prefix == "TUR" ){
        nodes <- nodes[nodes$pop > 5000,]
        edges <- edges[which(edges$from %in% nodes$id),]
        edges <- edges[which(edges$to %in% nodes$id),]
    }
    
    
    # Take only the a fraction of the top edges to increase speed
    
    edges_top <- edges %>% filter(flow > quantile(flow, top_percentage), time == 1)
    nodes_top <- nodes %>% filter(id %in% unique( c(edges_top$from, edges_top$to) ))
    net_top <- graph_from_data_frame(d = dplyr::select(edges_top, from, to, flow), vertices = nodes_top$id, directed = T)
    E(net_top)$weight <- edges_top$flow
    
    return(list(net_top, nodes_top, edges_top))
}


# prefix <- "ES"
# top_percentage <- 0.3
# out <- mobility_matrix(prefix, top_percentage)
# net <- out[[1]]
# nodes <- out[[2]]
# edges <- out[[3]]




