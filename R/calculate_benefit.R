#'Calculate ecosystem service benefit from a network
#'
#'@description `calculate_benefit` calculates the ecosystem services benefits from a given
#'  network and a given set of rules.
#'
#'@param net social-ecological network from which to calculate benefit (created using `network_create`)
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

calculate_benefit <- function(ee_network, es_network, alpha, beta, gamma, params = NULL) {
  # 1. calculate per node supply ----
  supply <- ee_network %>%
    dplyr::group_by(node1) %>%
    dplyr::summarise(node1_area = mean(node1_area),
                     connected_area = sum(node2_area)) %>%
    dplyr::mutate(supply = (node1_area^alpha)*exp(beta*connected_area))

  # 2. calculate per node benefit ----
  demand <- es_network %>%
    dplyr::inner_join(supply, by = c("node_supply" = "node1")) %>%
    dplyr::select(-connected_area)

  if(rival) {
    # rival
    competing <- demand %>% group_by(node_supply) %>%
      summarise(connected_demand = sum(demand_area))

    benefit <- demand %>%
      dplyr::inner_join(competing) %>%
      dplyr::mutate(prop_benefit = 1 / connected_demand) %>%
      group_by(node_demand, demand_area, node_supply, supply_area) %>%
      summarise(prop_benefit = sum(prop_benefit)) %>%
      mutate(supply = supply_area*prop_benefit) %>%
      group_by(node_demand, demand_area) %>%
      summarise(supply = sum(supply)) %>%
      mutate(benefit = (demand_area / gamma) * (1 - exp(-gamma * supply)))

  } else {
    # non-rival
    benefit <- demand %>%
      group_by(node_demand, demand_area) %>%
      summarise(supply = sum(supply_area)) %>%
      mutate(benefit = (demand_area / gamma) * (1 - exp(-gamma * supply)))
  }

  # 3. calculate total benefit ----
  benefit <- sum(benefit$benefit)

  # 4. output parameters ----
  params <- data.frame(params, benefit = benefit)
  return(params)
}
