#'Create ecological-ecological networks
#'
#'@description `create_ee_network` takes as input a spatial data layer (either real data or result of
#'  `ls_create`) and derives the underlying supply (ecological-ecological) network
#'
#'@param ls_supply Polygons containing ecosystem service supply areas
#'
#'@param ee_thresh Distance threshold for the ecological-ecological links
#'
#'@param supply_area Name of the column containing the supply area measure
#'
#'@param e2e Logical. If `TRUE` edge-to-edge distances between patches are calculated. If `FALSE` centroid-to-centroid distances are calculated (the latter is much quicker)
#'
#'@param params Vector containing the parameters used to generate the landscape if landscape is simulated (default = NULL)
#'
#'@return A list containing the network (and its attributes) and the parameters used to create the network
#'
#'@keywords ecosystem services, spatial, ecological system, neutral landscape model
#'
#'@export
#'
#'@import magrittr

create_ee_network <- function(ls_supply,
                              ee_thresh,
                              supply_area,
                              e2e = TRUE,
                              params = NULL) {

  # if no parameters are input, start the table here
  if(is.null(params)) {
    params <- data.frame(ee_thresh = ee_thresh)
  } else {
    params$ee_thresh <- ee_thresh
  }

  # turn into sf object
  if(is(ls_supply, "Spatial")) ls_supply <- sf::st_as_sf(ls_supply)

  # add on an ID column and rename the area column
  ls_supply <- dplyr::mutate(ls_supply, ID = 1:dplyr::n()) %>%
    dplyr::rename(area = !!supply_area) %>%
    sf::st_as_sf()

  # calculate all pairwise distances
  if(!e2e) {
    sf::st_agr(ls_supply) = "constant" # this removes the warning message
    ls_supply <- sf::st_centroid(ls_supply)
  }

  net_links <- sf::st_distance(ls_supply)
  ee_thresh <- ifelse(is.na(ee_thresh), -1, ee_thresh)
  net_links <- ifelse(net_links <= ee_thresh, 1, 0)

  # number of supply nodes
  params$num_supply <- nrow(ls_supply)

  # mean of area of supply nodes
  params$mean_supply_area <- mean(ls_supply$area)

  # standard deviation of area of supply nodes
  params$sd_supply_area <- sd(ls_supply$area)

  # calculate some network metrics
  ee_network <-  igraph::graph_from_adjacency_matrix(net_links, diag = TRUE, mode = "undirected")

  # mean edges per node
  params$ee_edge_per_node_mean <- mean(igraph::degree(ee_network, loops = TRUE, normalized = TRUE))

  # sd edges per node
  params$ee_edge_per_node_sd <- sd(igraph::degree(ee_network, loops = TRUE, normalized = TRUE))

  # edge density
  params$ee_density <- igraph::edge_density(ee_network, loops = TRUE)

  # for centrality measures we only calculate for nodes that are connected to at least one other node because we are
  # interested in the structure of the network for connected habitat patches to detect whether there is a single or a
  # only a few habitat patches responsible for connectivity (i.e., highly centralised), or not (i.e., not centralised)
  ee_network_con <- igraph::delete.vertices(ee_network, which(degree(ee_network) == 2))

  # closeness centralisation - note this is only valid for a fully connected network, so use with caution
  params$ee_centr_close <- igraph::centr_clo(ee_network_con, normalized = TRUE)$centralization

  # betweenness centralisation - note this is only valid for a fully connected network, so use with caution
  params$ee_centr_betw <- igraph::centr_betw(ee_network_con, directed = FALSE, normalized = TRUE)$centralization

  # degree centralisation
  params$ee_centr_degree <- igraph::centr_degree(ee_network_con, loops = TRUE, normalized = TRUE)$centralization

  # get network in correct format
  network <- net_links %>% tibble::as_tibble() %>%
    tibble::rownames_to_column("node1") %>%
    tidyr::gather(node2, link, -node1) %>%
    dplyr::mutate(node2 = stringr::str_replace(node2, "V", ""),
                  node1 = as.integer(node1),
                  node2 = as.integer(node2)) %>%
    dplyr::inner_join(ls_supply %>% sf::st_set_geometry(NULL), by = c("node1" = "ID")) %>%
    dplyr::inner_join(ls_supply %>% sf::st_set_geometry(NULL), by = c("node2" = "ID"),
                      suffix = c("_node1", "_node2"))

  return(list(network = network, params = params))
}
