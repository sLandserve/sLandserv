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

calculate_benefit <- function(net, rival, alpha, beta, gamma) {
  # generates ES benefit values for the networks
  # net = output from the function network_create() for one landscape
  # rival = whether the service is rival (TRUE) or non-rival (FALSE)
  # alpha = rate of production of the potential ecosystem service per unit area at small isolated supply nodes
  # beta = influence of connected supply nodes on the rate of production of the potential ecosystem service
  # gamma = marginal utility of the service at zero service used
  # eta = rate of decline of the marginal utility with increasing use of the service

  # get the ES production at each of the supply nodes
  # get supply and demand ids - note here that supply/demand nodes are not considered, only supply or demand
  supply_ids <- which(net$network$node_type == "supply")
  demand_ids <- which(net$network$node_type == "demand")
  #calculate ES production for each supply node
  ss_links <- net$network$net_links[supply_ids, supply_ids]
  s_areas <- net$network$node_areas[supply_ids]
  supply_ES <- apply(as.matrix(1:length(supply_ids)),
                     MARGIN = 1,
                     FUN = function(x, ss_links, s_areas, alpha, beta)
                       {if(!is.null(dim(ss_links))) {s_areas[x]^alpha * exp(beta * sum(ss_links[x,] * s_areas))} else {s_areas[x]^alpha + exp(beta * sum(ss_links * s_areas))}}, ss_links = ss_links, s_areas = s_areas, alpha = alpha, beta = beta)
  #calculate ES utility for each demand node
  ds_links <- net$network$net_links[demand_ids, supply_ids]
  d_areas <- net$network$node_areas[demand_ids]
  if (rival == FALSE) {
    # in this case the services are non-rival
    benefit_ES <- apply(as.matrix(1:length(demand_ids)),1,function(x, supply_ES, ds_links, s_areas, gamma) {
      if(!is.null(dim(ds_links))) {
        (d_areas[x] / gamma) * (1 - exp(-gamma * sum(ds_links[x,] * supply_ES)))
      } else {
        (d_areas[x] / gamma) * (1 - exp(-gamma * sum(ds_links * supply_ES)))
      }
    }, supply_ES = supply_ES, ds_links = ds_links, s_areas = s_areas, gamma = gamma)
  } else if (rival == TRUE) {
    # in this case the services are rival
    benefit_ES <- apply(as.matrix(1:length(demand_ids)), 1, function(x, supply_ES, ds_links, s_areas, gamma) {
      if(!is.null(dim(ds_links))) {
        (d_areas[x] / gamma) * (1 - exp(-gamma * sum(d_areas[x] * apply(ds_links, MARGIN = 2, FUN = function(y, d_areas) {
          if (sum(y) > 0) {1 / sum(y * d_areas)
        } else {0}}, d_areas = d_areas) * ds_links[x,] * supply_ES)))
      } else {
        (d_areas[x] / gamma) * (1 - exp(-gamma * sum(ds_links * supply_ES)))
      }
    }, supply_ES = supply_ES, ds_links = ds_links, s_areas = s_areas, gamma = gamma)
  }
  benefit <- sum(benefit_ES)
  output <- data.frame(net$params, rival = rival, alpha = alpha, beta = beta, gamma = gamma, benefit = benefit)

  return(output)
}
