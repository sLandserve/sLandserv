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
#'@param area_col column which contains the area information (NB this can be area, or modified area, such as by quality)
#'
#'@param params vector containing the parameters used to generate the landscape (if landscape is simulated, default = NULL)
#'
#'@return A list containing the network (and its attributes) and the parameters used to create the network
#'
#'@keywords ecosystem services, spatial, ecological system, neutral landscape model
#'
#'@export
create_es_network <- function(ls_supply, ls_demand, es_thresh, area_col = NULL, params = NULL) {

  # calculate all pairwise distances
  pts_supply <- sf::st_centroid(ls_supply)
  pts_demand <- sf::st_centroid(ls_demand)
  net_links <- sf::st_distance(pts_supply, pts_demand)

  net_links <- ifelse(net_links <= es_thresh, 1, 0)

  #number of supply nodes
  num_demand <- nrow(ls_demand)

  # network density
  es_network <-  network::network(as.matrix(net_links, directed=FALSE, loops=FALSE))
  es_density <- network::network.density(es_network)

  if(is.null(params)) {
    params <- data.frame(es_thresh = es_thresh, num_demand = num_demand, es_density = es_density)
  } else {
    params <- data.frame(params, es_thresh = es_thresh, num_demand = num_demand, es_density = es_density)
  }

  # get network in correct format
  network <- net_links %>% tibble::as_tibble() %>%
    tibble::rownames_to_column("node_supply") %>%
    tidyr::gather(node_demand, link, -node_supply) %>%
    dplyr::mutate(node_demand = stringr::str_replace(node_demand, "V", ""),
                  node_supply = as.integer(node_supply),
                  node_demand = as.integer(node_demand)) %>%
    dplyr::inner_join(ls_supply, by = c("node_supply" = "ID")) %>%
    dplyr::inner_join(ls_demand, by = c("node_demand" = "ID")) %>%
    dplyr::select(node_supply, node_demand, link, supply_area = patch_area.x, demand_area = patch_area.y)

  return(list(network = network, params = params))
}
