#'Create social-ecological networks from raster landscapes
#'
#'@description `network_create` takes as input a raster (either real data or result of
#'  `ls_create` and derives the underlying social ecological network
#'
#'@param ls raster containing cells classified as supply, demand and neutral
#'
#'@param params vector containing the parameters used to generate the landscape (if landscape is simulated)
#'
#'@param es_thresh distance threshold for the ecological-social links
#'
#'@param ee_thresh distance threshold for the ecological-ecological links (default value of `NULL` means no ecological-ecological links)
#'
#'@param ss_thresh distance threshold for the social-social links (default value of `NULL` means no social-social links)
#'
#'@return A list containing the network (and its attributes) and the parameters used to create the network
#'
#'@keywords ecosystem services, spatial, social ecological system, neutral landscape model
#'
#'@export
network_create <- function(ls, params = NULL, es_thresh, ee_thresh = NA, ss_thresh = NA) {
  # create attribute table
  all_nodes <- raster::rasterToPolygons(ls, dissolve=TRUE) %>%
    raster::disaggregate() %>%
    sf::st_as_sf() %>%
    dplyr::mutate(ID = as.character(1:n()),
           patch_type = dplyr::case_when(layer == 0 ~ "supply",
                                  layer == 1 ~ "neutral",
                                  TRUE ~ "demand"),
           patch_area = st_area(.)) %>%
    dplyr::select(-layer) %>%
    dplyr::filter(patch_type != "neutral") %>%
    dplyr::mutate(patch_code = dplyr::case_when(patch_type == "supply" ~ 0,
                                  patch_type == "demand" ~ 1,
                                  TRUE ~ NaN))

  # create the social-ecological network based on ee_thresh, es_thresh, and ss_thresh (ss_thresh currently not implemented)
  # this generates a list of the node attributes (node colour) and the link presence/absence
  # between each node
  supply_nodes <- which(all_nodes$patch_type == "supply")
  demand_nodes <- which(all_nodes$patch_type == "demand")

  net_links <- sf::st_distance(all_nodes)

  net_links[supply_nodes, demand_nodes] <- ifelse(net_links[supply_nodes, demand_nodes] <= es_thresh, 1, 0)
  net_links[demand_nodes, supply_nodes] <- ifelse(net_links[demand_nodes, supply_nodes] <= es_thresh, 1, 0)

  if(!is.na(ee_thresh)) {
    net_links[supply_nodes, supply_nodes] <- ifelse(net_links[supply_nodes, supply_nodes] <= ee_thresh, 1, 0)
  } else {
    net_links[supply_nodes, supply_nodes] <- 0
  }

  if(!is.na(ss_thresh)) {
    net_links[demand_nodes, demand_nodes] <- ifelse(net_links[demand_nodes, demand_nodes] <= ee_thresh, 1, 0)
  } else {
    net_links[demand_nodes, demand_nodes] <- 0
  }

  network <- list(node_code = all_nodes$patch_code, node_type = all_nodes$patch_type, node_areas = all_nodes$patch_area, net_links = net_links)

  if(is.null(params)) {
    params <- data.frame(es_thresh = es_thresh, ee_thres = ee_thresh, ss_thresh = ss_thresh)
  } else {
    params <- data.frame(params, es_thresh = es_thresh, ee_thres = ee_thresh, ss_thresh = ss_thresh)
  }

  return(list(network = network, params = params))
}
