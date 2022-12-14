#'Create social-ecological networks
#'
#'@description `create_es_network` takes as input a spatial data layer (either real data or result of
#' `ls_create`) and derives the underlying supply (ecological-ecological) network
#'
#'@param ls_supply Polygons containing ecosystem service supply areas
#'
#'@param ls_demand Polygons containing ecosystem service demand areas
#'
#'@param es_thresh Distance threshold for the social-ecological links
#'
#'@param excludable Whether the ecosystem service is excludable or not. If `TRUE` then only one demand node
#' (chosen at random) is connected to each supply node. If `FALSE' then any number of demand nodes
#' can be connected to each supply node
#'
#'@param supply_area Name of the column containing the supply area measure
#'
#'@param demand_area Name of the column containing the demand area measure
#'
#'@param e2e Logical. If `TRUE` edge-to-edge distances between patches are calculated. If `FALSE`
#' centroid-to-centroid distances are calculated (the latter is much quicker)
#'
#'@param params Vector containing the parameters used to generate the landscape
#' (if simulated landscapes used, default = NULL)
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
                              excludable,
                              supply_area,
                              demand_area,
                              e2e = TRUE,
                              params = NULL) {

  # if no parameters are input, start the table here
  if(is.null(params)) {
    params <- data.frame(es_thresh = es_thresh, excludable = excludable)
  } else {
    params$es_thresh <- es_thresh
    params$excludable <- excludable
  }

  # turn into sf object
  if(is(ls_supply, "Spatial")) ls_supply <- sf::st_as_sf(ls_supply)
  if(is(ls_demand, "Spatial")) ls_demand <- sf::st_as_sf(ls_demand)

  # add on an ID column
  ls_supply <- dplyr::mutate(ls_supply, ID = 1:dplyr::n()) %>%
    dplyr::rename(area_node = !!supply_area)
  ls_demand <- dplyr::mutate(ls_demand, ID = 1:dplyr::n()) %>%
    dplyr::rename(area_node = !!demand_area)

  # calculate all pairwise distances

  # if e2e is false, convert to centroid points before calculating distances
  if(!e2e) {
    sf::st_agr(ls_supply) <- "constant" # this removes the warning message
    sf::st_agr(ls_demand) <- "constant" # this removes the warning message
    ls_supply <- sf::st_centroid(ls_supply)
    ls_demand <- sf::st_centroid(ls_demand)
  }

  # get distance between supply and demand nodes
  net_links <- sf::st_distance(ls_supply, ls_demand)
  if (excludable == FALSE) {
      # not excludable - retain all supply-demand links
      net_links <- ifelse(net_links <= es_thresh, 1, 0)
  }
  else {
      # excludable - retain only one supply-demand link (chosen at random) per supply node
      net_links <- ifelse(net_links <= es_thresh, 1, 0) %>%
        apply(1, FUN = function(X){X[-sample(which(X == 1), 1)] <- 0; return(X)}) %>% t()
  }

  # number of demand nodes
  params$num_demand <- nrow(ls_demand)

  # mean of area of demand nodes
  params$mean_demand_area <- mean(ls_demand$area_node)

  # standard deviation of demand nodes
  params$sd_demand_area <- sd(ls_demand$area_node)

  # escape from function and return NA if no social-ecological links
  if(sum(net_links) == 0) {
    params$es_density <- NA
    params$es_edge_per_node_mean <- NA
    params$es_edge_per_node_sd <- NA
    params$es_centr_close <- NA
    params$es_centr_close_supply <- NA
    params$es_centr_close_demand <- NA
    params$es_centr_betw <- NA
    params$es_centr_betw_supply <- NA
    params$es_centr_betw_demand <- NA
    params$es_centr_degree <- NA
    params$es_centr_degree_supply <- NA
    params$es_centr_degree_demand <- NA

    return(list(network = NA, params = params))
  }

  # calculate some network metrics for the entire network, just the supply nodes, and just the demand nodes
  # with the normalisation denominators for these indices for bipartite networks taken from
  # Borgatti & Everett (1997). Social Networks 19:243-269.
  # also assign types to the network (supply = FALSE, demand = TRUE)
  es_network <-  igraph::graph_from_incidence_matrix(net_links, directed = FALSE, multiple = FALSE) %>%
                  igraph::set_vertex_attr("type", value = c(logical(length = nrow(net_links)), !logical(length = ncol(net_links))))
  # calculate number of each node type
  n_supply <- sum(!igraph::vertex_attr(es_network)$type)
  n_demand <- sum(igraph::vertex_attr(es_network)$type)

  # edge density
  # calculate edge density as the number of edges / total number of possible edges in the bipartite network
  params$es_density <- igraph::gsize(es_network) / (n_supply * n_demand)

  # mean edges per node
  params$es_edge_per_node_mean <- mean(c(igraph::degree(es_network, loops = FALSE, normalized = FALSE)[!igraph::vertex_attr(es_network)$type] /
                                                     n_demand,
                                         igraph::degree(es_network, loops = FALSE, normalized = FALSE)[igraph::vertex_attr(es_network)$type] /
                                                     n_supply))

  # sd edges per node
  params$es_edge_per_node_sd <- sd(c(igraph::degree(es_network, loops = FALSE, normalized = FALSE)[!igraph::vertex_attr(es_network)$type] /
                                                     n_demand,
                                         igraph::degree(es_network, loops = FALSE, normalized = FALSE)[igraph::vertex_attr(es_network)$type] /
                                                     n_supply))

  # for centrality measures we only calculate for nodes that are connected to at least one other node because we are
  # interested in centrality for the nodes that link supply and demand (e.g., to identify cases when one or a few supply
  # nodes are connected to many demand nodes, or when one or a few demand nodes are connected to many supply nodes)
  es_network_con <- igraph::delete_vertices(es_network, which(igraph::degree(es_network) == 0))
  n_supply_con <- sum(!igraph::vertex_attr(es_network_con)$type)
  n_demand_con <- sum(igraph::vertex_attr(es_network_con)$type)

  # closeness centralisation - note this is only valid for a fully connected network, so use with caution
  if (n_supply_con > n_demand_con){
    close_tmax <- 2 * (n_demand_con - 1) * ((n_demand_con + n_supply_con - 2) / (3 * n_demand_con + 4 * n_supply_con - 8)) +
                          2 * (n_supply_con - n_demand_con) * ((2 * n_demand_con - 1) / (5 * n_demand_con + 2 * n_supply_con - 6)) +
                          2 * (n_demand_con - 1) * ((n_supply_con - 2) / (2 * n_demand_con + 3 * n_supply_con - 6)) +
                          2 * ((n_demand_con - 1) / (n_supply_con + 4 * n_demand_con - 4))
  } else {
    close_tmax <- 2 * (n_supply_con - 1) * ((n_demand_con + n_supply_con - 4) / (3 * n_demand_con + 4 * n_supply_con - 8)) +
                          2 * (n_supply_con - 1) * ((n_supply_con - 2) / (2 * n_demand_con + 3 * n_supply_con - 6)) +
                          2 * (n_supply_con - 1) * ((n_demand_con - n_supply_con + 1) / (2 * n_demand_con + 3 * n_supply_con - 4))
  }
  if (n_supply_con > n_demand_con){
    close_tmax_supply <- ((n_demand_con - 1) * (n_supply_con - 2)) / (2 * n_supply_con - 3) +
                          ((n_demand_con - 1) * (n_supply_con - n_demand_con)) / (n_supply_con + n_demand_con - 2)
  } else {
    close_tmax_supply <- ((n_supply_con - 2) * (n_supply_con - 1)) / (2 * n_supply_con - 3)
  }
  if (n_demand_con > n_supply_con){
    close_tmax_demand <- ((n_supply_con - 1) * (n_demand_con - 2)) / (2 * n_demand_con - 3) +
                          ((n_supply_con - 1) * (n_demand_con - n_supply_con)) / (n_demand_con + n_supply_con - 2)
  } else {
    close_tmax_demand <- ((n_demand_con - 2) * (n_demand_con - 1)) / (2 * n_demand_con - 3)
  }
  close_vals <- igraph::closeness(es_network_con, normalized = FALSE)
  close_vals_supply <- igraph::closeness(es_network_con, normalized = FALSE)[!igraph::vertex_attr(es_network_con)$type]
  close_vals_demand <- igraph::closeness(es_network_con, normalized = FALSE)[igraph::vertex_attr(es_network_con)$type]
  params$es_centr_close <- igraph::centralize(close_vals, theoretical.max = close_tmax, normalized = TRUE)
  params$es_centr_close_supply <- igraph::centralize(close_vals_supply, theoretical.max = close_tmax_supply, normalized = TRUE)
  params$es_centr_close_demand <- igraph::centralize(close_vals_demand, theoretical.max = close_tmax_demand, normalized = TRUE)

  # betweenness centralisation - note this is only valid for a fully connected network, so use with caution
  if (n_supply_con > n_demand_con){
    betw_tmax <- 2 * (n_supply_con - 1) * (n_demand_con - 1) * (n_supply_con + n_demand_con - 1) -
                          (n_demand_con - 1) * (n_supply_con + n_demand_con - 2) -
                          0.5 * (n_supply_con - n_demand_con) * (n_supply_con + 3 * n_demand_con - 3)
  } else {
    betw_tmax <- (0.5 * n_demand_con * (n_demand_con - 1) + 0.5 * (n_supply_con - 1) * (n_supply_con - 2) + (n_supply_con - 1) * (n_demand_con - 2)) *
                          (n_supply_con + n_demand_con - 1) +
                          (n_supply_con - 1)
  }
  if (n_supply_con > n_demand_con){
    betw_tmax_supply <- 2 * (n_supply_con - 1) ^ 2 * (n_demand_con - 1)
  } else {
    betw_tmax_supply <- (n_supply_con - 1) *
                          (0.5 * n_demand_con * (n_demand_con - 1) + 0.5 * (n_supply_con - 1) * (n_supply_con - 2) + (n_supply_con - 1) * (n_demand_con - 1))
  }
  if (n_demand_con > n_supply_con){
    betw_tmax_demand <- 2 * (n_demand_con - 1) ^ 2 * (n_supply_con - 1)
  } else {
    betw_tmax_demand <- (n_demand_con - 1) *
                          (0.5 * n_supply_con * (n_supply_con - 1) + 0.5 * (n_demand_con - 1) * (n_demand_con - 2) + (n_demand_con - 1) * (n_supply_con - 1))
  }
  betw_vals <- igraph::betweenness(es_network_con, directed = FALSE, normalized = FALSE)
  betw_vals_supply <- igraph::betweenness(es_network_con, directed = FALSE, normalized = FALSE)[!igraph::vertex_attr(es_network_con)$type]
  betw_vals_demand <- igraph::betweenness(es_network_con, directed = FALSE, normalized = FALSE)[igraph::vertex_attr(es_network_con)$type]
  params$es_centr_betw <- igraph::centralize(betw_vals, theoretical.max = betw_tmax, normalized = TRUE)
  params$es_centr_betw_supply <- igraph::centralize(betw_vals_supply, theoretical.max = betw_tmax_supply, normalized = TRUE)
  params$es_centr_betw_demand <- igraph::centralize(betw_vals_demand, theoretical.max = betw_tmax_demand, normalized = TRUE)

  # degree centralisation
  degree_tmax <- (min(n_supply_con, n_demand_con) + min(n_supply_con, n_demand_con)) * max(n_supply_con, n_demand_con) -
                    2 * (max(n_supply_con, n_demand_con) + min(n_supply_con, n_demand_con) - 1) # note: corrected from Borgatti & Everett (1997) p.259
                                                                                                # equation should be (no + no) * ni - 2 * (ni + no - 1)
  degree_tmax_supply <- (n_supply_con - 1) * (n_demand_con - 1)
  degree_tmax_demand <- (n_supply_con - 1) * (n_demand_con - 1)
  degree_vals <- igraph::degree(es_network_con, loops = FALSE, normalized = FALSE)
  degree_vals_supply <- igraph::degree(es_network_con, loops = FALSE, normalized = FALSE)[!igraph::vertex_attr(es_network_con)$type]
  degree_vals_demand <- igraph::degree(es_network_con, loops = FALSE, normalized = FALSE)[igraph::vertex_attr(es_network_con)$type]
  params$es_centr_degree <- igraph::centralize(degree_vals, theoretical.max = degree_tmax, normalized = TRUE)
  params$es_centr_degree_supply <- igraph::centralize(degree_vals_supply, theoretical.max = degree_tmax_supply, normalized = TRUE)
  params$es_centr_degree_demand <- igraph::centralize(degree_vals_demand, theoretical.max = degree_tmax_demand, normalized = TRUE)

  # get network in correct format
  network <- net_links %>% tibble::as_tibble(.name_repair = "unique") %>%
    tibble::rownames_to_column("node_supply") %>%
    tidyr::gather(node_demand, link, -node_supply) %>%
    dplyr::mutate(node_demand = stringr::str_replace(node_demand, "...", ""),
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
