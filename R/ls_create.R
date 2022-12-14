#'Create landscapes of ecosystem service supply (ecological) and demand (social) patches
#'
#'@description `ls_create` simulates ecosystem service supply and demand landscapes based
#'  on a number of parameters: proportion of supply/demand; fragmentation of
#'  supply/demand; steepness of supply gradient; interspersion of supply/demand
#'
#'@param nrow Number of rows in the raster
#'
#'@param ncol Number of columns in the raster
#'
#'@param p_supply Proportion of the landscape which is ES supply
#'
#'@param p_demand Proportion of the landscape which is ES demand
#'
#'@param f_supply Fragmentation level of the ES supply (in the range [0, 1] with 1 being
#'  the most fragmented)
#'
#'@param f_demand Fragmentation level of the ES demand (in the range [0, 1] with 1 being
#'  the most fragmented)
#'
#'@param grad Steepness of spatial gradient for supply and demand (in the range [0, 1] with 1 being
#'  the steepest gradient)
#'
#'@param inter Interpersion between ES supply and demand along the spatial gradient (in
#'  the range [0, 1] with 1 being most interspersed)
#'
#'@return A list containing a polygon containing the supply patches, a
#'  polygon containing the demand patches, and the parameters used to generate them
#'
#'@keywords ecosystem services, spatial, social ecological system, neutral landscape model
#'
#'@export
#'
#'@import magrittr

ls_create <- function(nrow,
                      ncol,
                      p_supply,
                      p_demand,
                      f_supply,
                      f_demand,
                      grad,
                      inter) {

  max_dim <- max(nrow, ncol)
  N <- as.integer(ceiling(base::log(max_dim - 1, 2)))
  size <- 2 ^ N + 1
  nrow <- size
  ncol <- size

  params <- data.frame(nrow = nrow, ncol = ncol, p_supply = p_supply, p_demand = p_demand, f_supply = f_supply, f_demand = f_demand, inter = inter,
                grad = grad)

  # create a gradient surface for supply and then a gradient surface for demand
  #relative to the supply gradient based on the level of interspersion
  g_supply <- NLMR::nlm_planargradient(ncol,
                          nrow)
  g_demand <- 0.5 + ((-1 + inter * 2) * (g_supply - 0.5))

  # create supply and demand surfaces
  # here we control the level of fragmentation
  suppressWarnings(supply <- NLMR::nlm_mpd(ncol + 2, nrow + 2,
                    roughness = f_supply)) # add 2 to ncol and nrow since nlm_mdp removes outer edge
  suppressWarnings(demand <- NLMR::nlm_mpd(ncol + 2, nrow + 2,
                    roughness = f_demand)) # add 2 to ncol and nrow since nlm_mdp removes outer edge

  # check supply and demand rasters are the same size as g_supply and g_demand rasters
  if (raster::extent(supply) != raster::extent(g_supply)) {
    warning("supply gradient and factal landscapes not the same extent")
  }
  if (raster::extent(demand) != raster::extent(g_demand)) {
    warning("demand gradient and factal landscapes not the same extent")
  }

  # create the analysis landscape: this takes 3 steps:
  # 1. merge with gradients and control gradient steepness
  ls_supply <- landscapetools::util_merge(supply,
                   g_supply,
                   scalingfactor = grad)
  ls_demand <- landscapetools::util_merge(demand,
                   g_demand,
                   scalingfactor = grad)
  # 2. classify
  ls_supply <- landscapetools::util_classify(ls_supply,
                      weighting = c(1 - p_supply, p_supply),
                      level_names = c("neutral", "supply"))
  ls_demand <- landscapetools::util_classify(ls_demand,
                      weighting = c(1 - p_demand, p_demand),
                      level_names = c("neutral", "demand"))

  # 3. polygonise and split out supply and demand
  ls_supply <- raster::rasterToPolygons(ls_supply, dissolve=TRUE) %>%
    raster::disaggregate() %>%
    sf::st_as_sf() %>%
    dplyr::mutate(patch_area = sf::st_area(.)) %>%
    dplyr::filter(layer == 2) %>%
    dplyr::select(-layer)

  ls_demand <- raster::rasterToPolygons(ls_demand, dissolve=TRUE) %>%
    raster::disaggregate() %>%
    sf::st_as_sf() %>%
    dplyr::mutate(patch_area = sf::st_area(.)) %>%
    dplyr::filter(layer == 2) %>%
    dplyr::select(-layer)

  return(list(ls_supply = ls_supply, ls_demand = ls_demand, params = data.frame(params)))
}
