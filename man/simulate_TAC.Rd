% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_TAC.R
\name{simulate_TAC}
\alias{simulate_TAC}
\title{Simulate Total Allowable Catch Projections
#'
Determine the total allowable catch (TAC) in the next simulation
year based on the fishing mortality rate identified by a harvest
control rule (HCR). TAC is computed by projecting the population
forward one year under the specific level of F and computing the
the amount of catch that would be taken from that projected 
population structure based on the joint fishery selectivity curve.}
\usage{
simulate_TAC(hcr_F, naa, recruitment, joint_sel, dem_params)
}
\arguments{
\item{hcr_F}{fishing mortality rate identified by HCR}

\item{naa}{numbers-at-age vector (dimensions [1, nages, nsexes, nregions])}

\item{recruitment}{projected recruitment in next year}

\item{dem_params}{demographic parameter matrices subsetted to 1 year}

\item{join_sel}{joint fishery selectivity}
}
\description{
Simulate Total Allowable Catch Projections
#'
Determine the total allowable catch (TAC) in the next simulation
year based on the fishing mortality rate identified by a harvest
control rule (HCR). TAC is computed by projecting the population
forward one year under the specific level of F and computing the
the amount of catch that would be taken from that projected 
population structure based on the joint fishery selectivity curve.
}
