#'Create landscapes of ecosystem service supply (ecological) and demand (social) patches
#'
#'@description `ls_create` simulates ecosystem service supply and demand landscapes based
#'  on a number of parameters: proportion of supply/demand; fragmentation of
#'  supply/demand; interspersion of supply/demand
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
#'@param f_demand Fragmentation level of the ES supply (in the range [0, 1] with 1 being
#'  the most fragmented)
#'
#'@param inter Interpersion between ES supply and demand (in the range [0, 1] with 1 being
#'  completely interspersed)
#'
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
                      inter) {

  max_dim <- max(nrow, ncol)
  N <- as.integer(ceiling(base::log(max_dim - 1, 2)))
  size <- 2 ** N + 1
  nrow = size
  ncol = size

  params <- data.frame(nrow = nrow, ncol = ncol, p_supply = p_supply, p_demand = p_demand, f_supply = f_supply, f_demand = f_demand, inter = inter)

  # create a gradient surface for supply and a mirror image gradient surface for demand
  g_supply <- NLMR::nlm_planargradient(ncol,
                          nrow)
  g_demand <- 1 - g_supply

  # create supply and demand surfaces
  # here we control the fragmentation and the amount
  supply <- NLMR::nlm_mpd(ncol,
                    nrow,
                    roughness = f_supply,
                    verbose = FALSE)
  demand <- NLMR::nlm_mpd(ncol,
                    nrow,
                    roughness = f_demand,
                    verbose = FALSE)

  # create the analysis landscape: this takes 3 steps:
  # 1. merge gradients
  ls_supply <- landscapetools::util_merge(supply,
                   g_supply,
                   scalingfactor = 1 - inter)
  ls_demand <- landscapetools::util_merge(demand,
                   g_demand,
                   scalingfactor = 1 - inter)
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
    dplyr::filter(layer == 1) %>%
    dplyr::select(-layer)

  ls_demand <- raster::rasterToPolygons(ls_demand, dissolve=TRUE) %>%
    raster::disaggregate() %>%
    sf::st_as_sf() %>%
    dplyr::mutate(patch_area = sf::st_area(.)) %>%
    dplyr::filter(layer == 1) %>%
    dplyr::select(-layer)

  return(list(ls_supply = ls_supply, ls_demand = ls_demand, params = data.frame(params)))
}
