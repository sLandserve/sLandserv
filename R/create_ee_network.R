#'Create ecological-ecological networks
#'
#'@description `create_ee_network` takes as input a spatial data layer (either real data or result of
#'  `ls_create`) and derives the underlying supply (ecological-ecological) network
#'
#'@param ls_supply polygon containing ecosystem service supply areas
#'
#'@param ee_thresh distance threshold for the ecological-ecological links
#'
#'@param supply_area name of the column containing the supply area measure
#'
#'@param e2e logical. If `TRUE` edge-to-edge distances between patches are calculated. If `FALSE` centroid-to-centroid distances are calculated (the latter is much quicker)
#'
#'@param params vector containing the parameters used to generate the landscape if landscape is simulated (default = NULL)
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
  ls_supply <- dplyr::mutate(ls_supply, ID = 1:n()) %>%
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

  #number of supply nodes
  params$num_supply <- nrow(ls_supply)

  # calculate network density
  ee_network <-  igraph::graph_from_adjacency_matrix(net_links, diag = TRUE, mode = "undirected")
  params$ee_density <- igraph::edge_density(ee_network, loops = TRUE)
  params$ee_centr_betw <- igraph::centr_betw(ee_network, directed = FALSE)$centralization
  params$ee_centr_degree <- igraph::centr_degree(ee_network, loops = TRUE)$centralization
  params$ee_edge_per_node_mean <- mean(igraph::degree(ee_network, loops = TRUE))
  params$ee_edge_per_node_sd <- sd(igraph::degree(ee_network, loops = TRUE))

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
