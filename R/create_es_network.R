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
    params$es_centr_close <- NA
    params$es_centr_close_supply <- NA
    params$es_centr_close_demand <- NA
    params$es_centr_betw <- NA
    params$es_centr_betw_supply <- NA
    params$es_centr_betw_demand <- NA
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
  # the denominators for these indices for bipartite networks are taken from
  # Borgatti & Everett (1997). Social Networks 19:243-269.
  es_network <-  igraph::graph_from_incidence_matrix(net_links, directed = FALSE, multiple = FALSE)
  # assign types to the network (supply = FALSE, demand = TRUE)
  es_network$type <- c(logical(length = nrow(net_links)), !logical(length = ncol(net_links)))
  # calculate number of each node type
  n_supply <- sum(!es_network$type)
  n_demand <- sum(es_network$type)

  # edge density
  # calculate density as the number of edges / total number of possible edges in the bipartite network
  params$es_density <- igraph::gsize(es_network) / (sum(es_network$type) * sum(!es_network$type))

  # closeness centralisation
  if (n_supply > n_demand){
    close_tmax <- 2 * (n_demand - 1) * ((n_demand + n_supply - 2) / (3 * n_demand + 4 * n_supply - 8)) +
                          2 * (n_supply - n_demand) * ((2 * n_demand - 1) / (5 * n_demand + 2 * n_supply - 6)) +
                          2 * (n_demand - 1) * ((n_supply - 2) / (2 * n_demand + 3 * n_supply - 6)) +
                          2 * ((n_demand - 1) / (n_supply + 4 * n_demand - 4))
  } else {
    close_tmax <- 2 * (n_supply - 1) * ((n_demand + n_supply - 4) / (3 * n_demand + 4 * n_supply - 8)) +
                          2 * (n_supply - 1) * ((n_supply - 2) / (2 * n_demand + 3 * n_supply - 6)) +
                          2 * (n_supply - 1) * ((n_demand - n_supply + 1) / (2 * n_demand + 3 * n_supply - 4))
  }
  if (n_supply > n_demand){
    close_tmax_supply <- ((n_demand - 1) * (n_supply - 2)) / (2 * n_supply - 3) +
                          ((n_demand - 1) * (n_supply - n_demand)) / (n_supply + n_demand - 2)
  } else {
    close_tmax_supply <- ((n_supply - 2) * (n_supply - 1)) / (2 * n_supply - 3)
  }
  if (n_demand > n_supply){
    close_tmax_demand <- ((n_supply - 1) * (n_demand - 2)) / (2 * n_demand - 3) +
                          ((n_supply - 1) * (n_demand - n_supply)) / (n_demand + n_supply - 2)
  } else {
    close_tmax_demand <- ((n_demand - 2) * (n_demand - 1)) / (2 * n_demand - 3)
  }
  close_vals <- igraph::closeness(es_network, normalized = FALSE)
  close_vals_supply <- igraph::closeness(es_network, normalized = FALSE)[!es_network$type]
  close_vals_demand <- igraph::closeness(es_network, normalized = FALSE)[es_network$type]
  params$es_centr_close <- igraph::centralize(close_vals, theoretical.max = close_tmax, normalized = TRUE)
  params$es_centr_close_supply <- igraph::centralize(close_vals_supply, theoretical.max = close_tmax_supply, normalized = TRUE)
  params$es_centr_close_demand <- igraph::centralize(close_vals_demand, theoretical.max = close_tmax_demand, normalized = TRUE)

  # betweenness centralisation
  if (n_supply > n_demand){
    betw_tmax <- 2 * (n_supply - 1) * (n_demand - 1) * (n_supply + n_demand - 1) -
                          (n_demand - 1) * (n_supply + n_demand - 2) -
                          0.5 * (n_supply - n_demand) * (n_supply + 3 * n_demand - 3)
  } else {
    betw_tmax <- (0.5 * n_demand * (n_demand - 1) + 0.5 * (n_supply - 1) * (n_supply - 2) + (n_supply - 1) * (n_demand - 2)) *
                          (n_supply + n_demand - 1) +
                          (n_supply - 1)
  }
  if (n_supply > n_demand){
    betw_tmax_supply <- 2 * (n_supply - 1) ^ 2 * (n_demand - 1)
  } else {
    betw_tmax_supply <- (n_supply - 1) *
                          (0.5 * n_demand * (n_demand - 1) + 0.5 * (n_supply - 1) * (n_supply - 2) + (n_supply - 1) * (n_demand - 1))
  }
  if (n_demand > n_supply){
    betw_tmax_demand <- 2 * (n_demand - 1) ^ 2 * (n_supply - 1)
  } else {
    betw_tmax_demand <- (n_demand - 1) *
                          (0.5 * n_supply * (n_supply - 1) + 0.5 * (n_demand - 1) * (n_demand - 2) + (n_demand - 1) * (n_supply - 1))
  }
  betw_vals <- igraph::betweenness(es_network, directed = FALSE, normalized = FALSE)
  betw_vals_supply <- igraph::betweenness(es_network, directed = FALSE, normalized = FALSE)[!es_network$type]
  betw_vals_demand <- igraph::betweenness(es_network, directed = FALSE, normalized = FALSE)[es_network$type]
  params$es_centr_betw <- igraph::centralize(betw_vals, theoretical.max = betw_tmax, normalized = TRUE)
  params$es_centr_betw_supply <- igraph::centralize(betw_vals_supply, theoretical.max = betw_tmax_supply, normalized = TRUE)
  params$es_centr_betw_demand <- igraph::centralize(betw_vals_demand, theoretical.max = betw_tmax_demand, normalized = TRUE)

  # degree centralisation
  degree_tmax <- (max(n_supply, n_demand) + min(n_supply, n_demand)) * max(n_supply, n_demand) -
                    2 * (max(n_supply, n_demand) + min(n_supply, n_demand) - 1)
  degree_tmax_supply <- (n_supply - 1) * (n_demand - 1)
  degree_tmax_demand <- (n_supply - 1) * (n_demand - 1)
  degree_vals <- igraph::degree(es_network, loops = FALSE, normalized = FALSE)
  degree_vals_supply <- igraph::degree(es_network, loops = FALSE, normalized = FALSE)[!es_network$type]
  degree_vals_demand <- igraph::degree(es_network, loops = FALSE, normalized = FALSE)[es_network$type]
  params$es_centr_degree <- igraph::centralize(degree_vals, theoretical.max = degree_tmax, normalized = TRUE)
  params$es_centr_degree_supply <- igraph::centralize(degree_vals_supply, theoretical.max = degree_tmax_supply, normalized = TRUE)
  params$es_centr_degree_demand <- igraph::centralize(degree_vals_demand, theoretical.max = degree_tmax_demand, normalized = TRUE)

  # mean edges per node
  params$es_edge_per_node_mean <- mean(c(igraph::degree(es_network, loops = FALSE, normalized = FALSE)[es_network$type] /
                                                     sum(!es_network$type),
                                         igraph::degree(es_network, loops = FALSE, normalized = FALSE)[!es_network$type] /
                                                     sum(es_network$type)))
  params$es_edge_per_node_mean_supply <- mean(igraph::degree(es_network, loops = FALSE, normalized = FALSE)[!es_network$type] /
                                                     sum(es_network$type))
  params$es_edge_per_node_mean_demand <- mean(igraph::degree(es_network, loops = FALSE, normalized = FALSE)[es_network$type] /
                                                     sum(!es_network$type))

  # sd edges per node
  params$es_edge_per_node_sd <- sd(c(igraph::degree(es_network, loops = FALSE, normalized = FALSE)[es_network$type] /
                                                     sum(!es_network$type),
                                         igraph::degree(es_network, loops = FALSE, normalized = FALSE)[!es_network$type] /
                                                     sum(es_network$type)))
  params$es_edge_per_node_sd_supply <- sd(igraph::degree(es_network, loops = FALSE, normalized = FALSE)[!es_network$type] /
                                                     sum(es_network$type))

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
