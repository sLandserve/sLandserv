#'Create network summary statistics
#'
#'@description `network_analyse` takes as input a network object (output of `network_create`) and calculates the number and density of each type of nodes
#'
#'@param net network object (output of `network_create`)
#'
#'@return Updated parameters with number and density of nodes
#'
#'@keywords ecosystem services, spatial, social ecological system, neutral landscape model
#'
#'@export
network_analyse <- function(net) {
  #number of supply nodes
  x$params$num_supply <- length(which(network$node_type == "supply"))
  #number of demand nodes
  x$params$num_demand <- length(which(network$node_type == "demand"))
  #density of ecological-ecological (supply-supply) network
  if (x$params$num_supply > 1) {
    EE_network <-  network::network(as.matrix(network$net_links[which(network$node_type == "supply"), which(network$node_type == "supply")]), directed=FALSE, loops=FALSE)
    x$params$ee_density <- network::network.density(EE_network)
  } else {
    x$params$ee_density <- as.matrix(network$net_links[which(network$node_type == "supply"), which(network$node_type == "supply")])[1,1]
  }

  #density of social-ecological (demand-supply) bipartite network
  SE_matrix <- matrix(0,nrow=x$params$num_demand + x$params$num_supply, ncol=x$params$num_demand + x$params$num_supply)
  SE_matrix[1:x$params$num_demand, (x$params$num_demand + 1):ncol(SE_matrix)] <- network$net_links[which(network$node_type == "demand"), which(network$node_type == "supply")]
  SE_matrix[(x$params$num_demand + 1):nrow(SE_matrix), 1:x$params$num_demand] <- network$net_links[which(network$node_type == "supply"), which(network$node_type == "demand")]
  SE_network <-  network::network(SE_matrix, bipartite=x$params$num_demand, directed = FALSE)
  x$params$se_density <- network::network.density(SE_network, discount.bipartite=TRUE)

  return(x$params)
}
