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
#'@param params vector containing the parameters used to generate the landscape if landscape is simulated (default = NULL)
#'
#'@return A list containing the network (and its attributes) and the parameters used to create the network
#'
#'@keywords ecosystem services, spatial, ecological system, neutral landscape model
#'
#'@export
create_ee_network <- function(ls_supply, ee_thresh, supply_area, params = NULL) {

  # turn into sf object
  if(is(ls_supply, "Spatial")) ls_supply <- sf::st_as_sf(ls_supply)

  # add on an ID column and rename the area column
  ls_supply <- dplyr::mutate(ls_supply, ID = 1:n()) %>%
    rename(area = !!supply_area)

  # calculate all pairwise distances
  sf::st_agr(ls_supply) = "constant" # this removes the warning message
  pts <- sf::st_centroid(ls_supply)
  net_links <- sf::st_distance(pts)

  net_links <- ifelse(net_links <= ee_thresh, 1, 0)

  #number of supply nodes
  num_nodes <- nrow(ls_supply)

  ee_network <-  network::network(as.matrix(net_links, directed=FALSE, loops=TRUE))
  ee_density <- network::network.density(ee_network)

  if(is.null(params)) {
    params <- data.frame(ee_thresh = ee_thresh, num_supply = num_nodes, ee_density = ee_density)
  } else {
    params <- data.frame(params, ee_thresh = ee_thresh, num_supply = num_nodes, ee_density = ee_density)
  }

  # get network in correct format
  network <- net_links %>% tibble::as_tibble() %>%
    tibble::rownames_to_column("node1") %>%
    tidyr::gather(node2, link, -node1) %>%
    dplyr::mutate(node2 = stringr::str_replace(node2, "V", ""),
                  node1 = as.integer(node1),
                  node2 = as.integer(node2)) %>%
    dplyr::inner_join(ls_supply %>% sf::st_set_geometry(NULL), by = c("node1" = "ID")) %>%
    dplyr::inner_join(ls_supply %>% sf::st_set_geometry(NULL), by = c("node2" = "ID"),
                      suffix = c("_node1", "_node2")) %>%
    dplyr::filter(link == 1) %>%
    dplyr::select(-link)

  return(list(network = network, params = params))
}
