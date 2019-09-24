#'Create social-ecological networks
#'
#'@description `create_es_network` takes as input a spatial data layer (either real data or result of
#'  `ls_create`) and derives the underlying supply (ecological-ecological) network
#'
#'@param ls_supply polygon containing ecosystem service supply areas
#'
#'@param ls_demand polygon containing ecosystem service demand areas
#'
#'@param es_thresh distance threshold for the social-ecological links
#'
#'@param supply_area name of the column containing the supply area measure
#'
#'@param demand_area name of the column containing the demand area measure
#'
#'#'@param e2e logical. If `TRUE` edge-to-edge distances between patches are calculated. If `FALSE` centroid-to-centroid distances are calculated (the latter is much quicker)
#'
#'@param params vector containing the parameters used to generate the landscape (if landscape is simulated, default = NULL)
#'
#'@return A list containing the network (and its attributes) and the parameters used to create the network
#'
#'@keywords ecosystem services, spatial, ecological system, neutral landscape model
#'
#'@export
#'
#'@import magrittr

create_es_network <- function(ls_supply,
                              ls_demand,
                              es_thresh,
                              supply_area,
                              demand_area,
                              e2e = TRUE,
                              params = NULL) {

  # if no parameters are input, start the table here
  if(is.null(params)) {
    params <- data.frame(es_thresh = es_thresh)
  } else {
    params$es_thresh <- es_thresh
  }

  # turn into sf object
  if(is(ls_supply, "Spatial")) ls_supply <- sf::st_as_sf(ls_supply)
  if(is(ls_demand, "Spatial")) ls_demand <- sf::st_as_sf(ls_demand)

  # add on an ID column
  ls_supply <- dplyr::mutate(ls_supply, ID = 1:n()) %>%
    dplyr::rename(area = !!supply_area)
  ls_demand <- dplyr::mutate(ls_demand, ID = 1:n()) %>%
    dplyr::rename(area = !!demand_area)

  # calculate all pairwise distances

  # if e2e is false, convert to centroid points before calculating distances
  if(!e2e) {
    sf::st_agr(ls_supply) = "constant" # this removes the warning message
    sf::st_agr(ls_demand) = "constant" # this removes the warning message
    ls_supply <- sf::st_centroid(ls_supply)
    ls_demand <- sf::st_centroid(ls_demand)

  }

  net_links <- sf::st_distance(ls_supply, ls_demand)
  net_links <- ifelse(net_links <= es_thresh, 1, 0)

  #number of demand nodes
  params$num_demand <- nrow(ls_demand)

  # escape from function and return NA if no social-ecological links
  if(sum(net_links) == 0) {
    params$es_density <- NA
    params$es_density_supply <- NA
    params$es_density_demand <- NA
    params$es_centr_degree <- NA
    params$es_centr_degree_supply <- NA
    params$es_centr_degree_demand <- NA
    params$es_edge_per_node_mean <- NA
    params$es_edge_per_node_mean_supply <- NA
    params$es_edge_per_node_mean_demand <- NA
    params$es_edge_per_node_sd <- NA
    params$es_edge_per_node_sd_supply <- NA
    params$es_edge_per_node_sd_demand <- NA
    return(list(network = NA, params = params))
  }

  # calculate some network metrics for the entire network, the supply nodes, and the demand nodes
  es_network <-  igraph::graph_from_incidence_matrix(net_links, directed = FALSE, multiple = FALSE)
  # assign types to the network (supply = FALSE, demand = TRUE)
  es_network$type <- c(logical(length = nrow(net_links)), !logical(length = ncol(net_links)))
  # calculate density as the number of edges / total number of possible edges in the bipartite network
  params$es_density <- igraph::gsize(es_network) / (sum(es_network$type)*sum(!es_network$type))
  params$es_density_supply <- igraph::gsize(es_network) / (sum(!es_network$type))
  params$es_density_demand <- igraph::gsize(es_network) / (sum(es_network$type))

  # OLD LOGIC - NEED TO CHECK
  # on assumption that each node of type a can have maximum degree = num nodes of type b
  # for type a, there are N_a nodes each with degree = N_b ??
  # for type b, there are N_b nodes each with degree = N_a ??
  # so theoretical maximum degree centralisation for the whole graph is 2*N_a*N_b ??
  # degree_tmax <- 2*sum(es_network$type)*sum(!es_network$type)
  # degree_vals <- igraph::degree(es_network)
  # params$es_centr_degree <- igraph::centralize(degree_vals, theoretical.max = degree_tmax, normalized = TRUE)
  # edges per node based on degree centrality normalised by the number of possible edges
  # params$es_edge_per_node_mean <- mean(c(igraph::degree(es_network, loops = FALSE, normalized = FALSE)[es_network$type] /
  #                                                    sum(!es_network$type),
  #                                        igraph::degree(es_network, loops = FALSE, normalized = FALSE)[!es_network$type] /
  #                                                    sum(es_network$type)))
  # params$es_edge_per_node_sd <- sd(c(igraph::degree(es_network, loops = FALSE, normalized = FALSE)[es_network$type] /
  #                                                    sum(!es_network$type),
  #                                        igraph::degree(es_network, loops = FALSE, normalized = FALSE)[!es_network$type] /
  #                                                    sum(es_network$type)))

  # NEW LOGIC - NEED TO CHECK
  # for the whole network, maximum degree centralisation is when one node of the type with the lowest number of nodes
  # is connected to all nodes of the other type and the other nodes of the type with the lowest number
  # of nodes are connected to nothing; in this case the maximum degree centralisation (sum(d_max - d)) is
  # max(N_a,N_b) * [max(N_a,N_b) - 1] + [min(N_a,N_b) - 1] * max(N_a,N_b)
  degree_tmax <- max(sum(es_network$type),sum(!es_network$type))) * (1 - max(sum(es_network$type),sum(!es_network$type)))) +
                 (min(sum(es_network$type),sum(!es_network$type))) - 1) * max(sum(es_network$type),sum(!es_network$type)))
  # for only one node type, the maximum degree centralisation is when one node of focal type is connected to all nodes of
  # the other type and all other nodes are connected to nothing; in this case the maximum degree centralisation (sum(d_max - d)) is
  # (N_a - 1) * N_b for node type a and (N_b - 1) * N_a for node type b
  degree_tmax_supply <- (sum(!es_network$type) - 1) * sum(es_network$type)
  degree_tmax_demand <- (sum(es_network$type) - 1) * sum(!es_network$type)
  degree_vals <- igraph::degree(es_network, loops = FALSE, normalized = FALSE)
  degree_vals_supply <- igraph::degree(es_network, loops = FALSE, normalized = FALSE)[!es_network$type]
  degree_vals_demand <- igraph::degree(es_network, loops = FALSE, normalized = FALSE)[es_network$type]
  params$es_centr_degree <- igraph::centralize(degree_vals, theoretical.max = degree_tmax, normalized = TRUE)
  params$es_centr_degree_supply <- igraph::centralize(degree_vals_supply, theoretical.max = degree_tmax_supply, normalized = TRUE)
  params$es_centr_degree_demand <- igraph::centralize(degree_vals_demand, theoretical.max = degree_tmax_demand, normalized = TRUE)
  # edges per node based on degree centrality normalised by the number of possible edges
  params$es_edge_per_node_mean <- mean(c(igraph::degree(es_network, loops = FALSE, normalized = FALSE)[es_network$type] /
                                                     sum(!es_network$type),
                                         igraph::degree(es_network, loops = FALSE, normalized = FALSE)[!es_network$type] /
                                                     sum(es_network$type)))
  params$es_edge_per_node_sd <- sd(c(igraph::degree(es_network, loops = FALSE, normalized = FALSE)[es_network$type] /
                                                     sum(!es_network$type),
                                         igraph::degree(es_network, loops = FALSE, normalized = FALSE)[!es_network$type] /
                                                     sum(es_network$type)))
  params$es_edge_per_node_mean_supply <- mean(igraph::degree(es_network, loops = FALSE, normalized = FALSE)[!es_network$type] /
                                                     sum(es_network$type))
  params$es_edge_per_node_sd_supply <- sd(igraph::degree(es_network, loops = FALSE, normalized = FALSE)[!es_network$type] /
                                                     sum(es_network$type))
  params$es_edge_per_node_mean_demand <- mean(igraph::degree(es_network, loops = FALSE, normalized = FALSE)[es_network$type] /
                                                     sum(!es_network$type))
  params$es_edge_per_node_sd_demand <- sd(igraph::degree(es_network, loops = FALSE, normalized = FALSE)[es_network$type] /
                                                     sum(!es_network$type)))

  # get network in correct format
  network <- net_links %>% tibble::as_tibble() %>%
    tibble::rownames_to_column("node_supply") %>%
    tidyr::gather(node_demand, link, -node_supply) %>%
    dplyr::mutate(node_demand = stringr::str_replace(node_demand, "V", ""),
                  node_supply = as.integer(node_supply),
                  node_demand = as.integer(node_demand)) %>%
    dplyr::inner_join(ls_supply %>% sf::st_set_geometry(NULL),
                      by = c("node_supply" = "ID")) %>%
    dplyr::inner_join(ls_demand %>% sf::st_set_geometry(NULL),
                      by = c("node_demand" = "ID"),
                      suffix = c("_supply", "_demand")) %>%
    dplyr::filter(link == 1) %>%
    dplyr::select(-link)

  return(list(network = network, params = params))
}
