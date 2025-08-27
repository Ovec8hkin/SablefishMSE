#' Default Operating Model for Alaska Sablefish
#' 
#' Dataset containing demographic matrices, observation
#' process parameters, and initial numbers-at-age for 
#' the alaska sablefish population. This data is used as
#' the base for all operating models used in this 
#' package.
#' 
#' @format list containing OM components
#' 
"sable_om"

#' Sablefish assessment
#'
#' A dataset containing sablefish assessment results.
#' This data is used for defining the historical recruitment
#' timeseries.
#'
#' @format A list containing assessment data
#' 
"sablefish_assessment_data"

#' Random Recruitment OM Object
#' 
#' An operating model object that will generate future
#' recruitment by randomly resampling from teh historical
#' recruitment timeseries
#' 
#' @format An operating model object
#' 
"om_rand_recruit"

#' Cylic Recruitment OM Object
#' 
#' An operating model object that will generate future
#' recruitment via a cyclic pattern of Beverton-Holt SRRs.
#' The first cycle is a low recruitment phase lasting 20 
#' years, the second cycle is a high recruitment phase
#' lasting 5 years.
#' 
#' @format An operating model object
#' 
"om_bhcyclic_recruit"

#' Crash Recruitment OM Object
#' 
#' An operating model object that will generate future
#' recruitment via a recruitment crash for the first 20
#' years, and then resort to resampling from the historical
#' recruitment timeseries.
#' 
#' @format An operating model object
#' 
"om_crash_recruit"

# F40 Managament Procedure
#'
#' An example management procedure that replicates
#' the NPFMC F40 harvest control rule policy that
#' is operationally used to manage sablefish in
#' Alaska.
#' 
#' @format A management procedure object.
#' 
"mp_f40"

# Harvest Cap Managament Procedure
#'
#' An example management procedure that applies
#' a 25,000 mt cap on annual harvest to the 
#' operational F40% harvest control rule policy.
#' 
#' @format A management procedure object.
#' 
"mp_25cap"

# Stability Constraint Managament Procedure
#'
#' An example management procedure that applied
#' a symmetric 5% stability constraint to the
#' operational F40% harvest control rule policy.
#' 
#' @format A management procedure object.
#' 
"mp_5perc"