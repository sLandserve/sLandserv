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
#'@param params vector containing the parameters used to generate the landscape (if landscape is simulated, default = NULL)
#'
#'@return A list containing the network (and its attributes) and the parameters used to create the network
#'
#'@keywords ecosystem services, spatial, ecological system, neutral landscape model
#'
#'@export
create_es_network <- function(ls_supply,
                              ls_demand,
                              es_thresh,
                              supply_area,
                              demand_area,
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
    rename(area = !!supply_area)
  ls_demand <- dplyr::mutate(ls_demand, ID = 1:n()) %>%
    rename(area = !!demand_area)

  # calculate all pairwise distances
  sf::st_agr(ls_supply) = "constant" # this removes the warning message
  sf::st_agr(ls_demand) = "constant" # this removes the warning message
  pts_supply <- sf::st_centroid(ls_supply)
  pts_demand <- sf::st_centroid(ls_demand)
  net_links <- sf::st_distance(pts_supply, pts_demand)

  net_links <- ifelse(net_links <= es_thresh, 1, 0)

  #number of supply nodes
  params$num_demand <- nrow(ls_demand)

  # escape from function and return NA if no social-ecological links
  if(sum(net_links) == 0) {
    params$es_density <- NA
    return(list(network = NA, params = params))
  }

  es_network <-  network::network(as.matrix(net_links, directed=FALSE, loops=FALSE))
  params$es_density <- network::network.density(es_network)

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
