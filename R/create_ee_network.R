#'Create ecological-ecological networks
#'
#'@description `create_ee_network` takes as input a spatial data layer (either real data or result of
#'  `ls_create`) and derives the underlying supply (ecological-ecological) network
#'
#'@param ls_supply polygon containing ecosystem service supply areas
#'
#'@param ee_thresh distance threshold for the ecological-ecological links
#'
#'@param area_col column which contains the area information (NB this can be area, or modified area, such as by quality)
#'
#'@param params vector containing the parameters used to generate the landscape (if landscape is simulated, default = NULL)
#'
#'@return A list containing the network (and its attributes) and the parameters used to create the network
#'
#'@keywords ecosystem services, spatial, ecological system, neutral landscape model
#'
#'@export
create_ee_network <- function(ls_supply, ee_thresh, area_col = NULL, params = NULL) {

  # turn into sf object
  if(class(ls_supply) == "sp") ls_supply <- sf::st_as_sf(ls_supply)

  # rename the area column for comparison
  if(!is.null(area_col)) ls_supply <- dplyr::mutate(ls_supply, patch_area = get(area_col))

  # calculate all pairwise distances
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
    dplyr::inner_join(ls_supply, by = c("node1" = "ID")) %>%
    dplyr::inner_join(ls_supply, by = c("node2" = "ID")) %>%
    dplyr::select(node1, node2, link, node1_area = patch_area.x, node2_area = patch_area.y)

  return(list(network = network, params = params))
}
