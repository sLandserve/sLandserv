#'Calculate ecosystem service benefit from a network
#'
#'@description `calculate_benefit` calculates the ecosystem services benefits from a given
#'  network and a given set of rules.
#'
#'@param ee_network ecological-ecological network from which to calculate supply (created using `create_ee_network`)
#'
#'@param es_network social-ecological network from which to calculate benefit (created using `create_es_network`)
#'
#'@param rival whether the services are rival (TRUE) or non-rival (FALSE)
#'
#'@param alpha the rate of production of the potential ecosystem service per unit area at supply nodes (value > 0)
#'
#'@param beta the influence of connected supply nodes (i.e., the ecological-ecological links) on the rate of production of the potential ecosystem service at supply nodes (value > 0)
#'
#'@param gamma marginal utility of the service at zero service used
#'
#'@return A vector containing the ecosystem service benefit and the parameters used to generate the network
#'
#'@keywords ecosystem services, spatial, social ecological system, neutral landscape model
#'
#'@export

calculate_benefit <- function(ee_network, es_network, rival, alpha, beta, gamma, params = NULL) {

  # if no parameters are input, start the table here
  if(is.null(params)) {
    params <- data.frame(rival = rival, alpha = alpha, beta = beta, gamma = gamma)
  } else {
    params <- data.frame(params, rival = rival, alpha = alpha, beta = beta, gamma = gamma)
  }

  # if there is no social-ecological network, escape function and return 0
  if(!is(es_network, "data.frame")) {
    params$benefit <- 0
    return(params)
  }

  # 1. calculate per node supply ----
  supply <- ee_network %>%
    dplyr::group_by(node1) %>%
    dplyr::summarise(area_node1 = mean(area_node1))

  connected_supply <- ee_network %>%
    dplyr::filter(link == 1) %>%
    dplyr::group_by(node1) %>%
    dplyr::summarise(connected_supply = sum(area_node2))

  supply <- left_join(supply, connected_supply) %>%
    dplyr::mutate(connected_supply = case_when(is.na(connected_supply) ~ 0,
                                               TRUE ~ connected_supply)) %>%
    dplyr::mutate(supply = (area_node1^alpha)*exp(beta*connected_supply))

  # 2. calculate per node benefit ----
  demand <- es_network %>%
    dplyr::inner_join(supply, by = c("node_supply" = "node1")) %>%
    dplyr::select(-connected_supply)

  if(rival) {
    # rival
    competing <- demand %>%
      dplyr::group_by(node_supply) %>%
      dplyr::summarise(connected_demand_area = sum(area_demand),
                       connected_demand_no = n())

    benefit <- demand %>%
      dplyr::inner_join(competing, by = "node_supply") %>%
      dplyr::mutate(prop_benefit = area_demand / connected_demand_area) %>%
      dplyr::group_by(node_demand, area_demand, node_supply, supply) %>%
      dplyr::summarise(prop_benefit = sum(prop_benefit)) %>%
      dplyr::mutate(supply = supply*prop_benefit) %>%
      dplyr::group_by(node_demand, area_demand) %>%
      dplyr::summarise(supply = sum(supply)) %>%
      dplyr::mutate(benefit = (area_demand / gamma) * (1 - exp(-gamma * supply)))

  } else {
    # non-rival
    benefit <- demand %>%
      dplyr::group_by(node_demand, area_demand) %>%
      dplyr::summarise(supply = sum(supply)) %>%
      dplyr::mutate(benefit = (area_demand / gamma) * (1 - exp(-gamma * supply)))
  }

  # 3. calculate total benefit ----
  params$benefit <- sum(benefit$benefit)

  # 4. output parameters ----
  return(params)
}
