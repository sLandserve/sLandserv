#'Create nodes from a raster
#'
#'@description `raster_to_node` takes as input a raster layer and converts to nodes of one type
#'
#'@param ls raster containing the supply-demand landscape
#'
#'@param type_id Identifier for the cells to be converted to nodes
#'
#'@return An sf object containing nodes of the specified type. Area of each node is also calculated and added to the attributes
#'
#'@keywords ecosystem services, spatial, ecological system, neutral landscape model
#'
#'@export
raster_to_node <- function(ls,
                           type_id) {
  nodes <- raster::rasterToPolygons(ls$ls, dissolve=TRUE) %>%
    raster::disaggregate() %>%
    sf::st_as_sf() %>%
    dplyr::filter(layer == type_id) %>%
    dplyr::mutate(patch_area = sf::st_area(.)) %>%
    dplyr::select(-layer)

  return(nodes)
}
