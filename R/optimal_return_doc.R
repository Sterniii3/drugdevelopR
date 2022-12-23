#' @name optimal_return_doc
#' 
#' Generates the documentation of the return value
#' of an optimal function including some custom text.
#' 
#' @param type string deciding whether this is return text for 
#'  normal, binary or time-to-event endpoints
#' @param setting string containing the setting, i.e. "basic", "bias",
#' "multitrial"
#'
#' @return string containing the documentation of the return value.
optimal_return_doc <- function(type,
                               setting = "basic"){
  custom_threshold = ""
  custom_further = ""
  steps1_str = ""
  stepm1_str = ""
  stepl1_str = ""
  if(type == "normal"){
    custom_threshold = "\\item{Kappa}{optimal threshold value for the decision rule to go to phase III}"
  }
  if(type == "binary"){
    custom_threshold = "\\item{RRgo}{optimal threshold value for the decision rule to go to phase III}"
  }
  if(type == "tte"){
    custom_threshold = "\\item{HRgo}{optimal threshold value for the decision rule to go to phase III}
                        \\item{d2}{optimal total number of events for phase II}
                        \\item{d3}{total expected number of events for phase III; rounded to next natural number}
                        \\item{d}{total expected number of events in the program; d = d2 + d3}"
  }
  if(setting == "bias"){
    custom_further = "\\item{Method}{Type of adjustment: \"multipl.\" (multiplicative adjustment of effect size), \"add.\" (additive adjustment of effect size), \"multipl2.\" (multiplicative adjustment of effect size and threshold), \"add2.\" (additive adjustment of effect size and threshold)}\n \\item{Adj}{optimal adjustment parameter (lambda or alphaCI according to Method)}"
  }
  if(setting == "multitrial"){
    custom_further = "\\item{Case}{Case: \"number of significant trials needed\"}\\item{Strategy}{Strategy: \"number of trials to be conducted in order to achieve the goal of the case\"}"
    if(type == "tte" | type == "binary"){
      steps1_str = " (lower boundary in HR scale is set to 1, as proposed by IQWiG (2016))"
      stepm1_str = " (lower boundary in HR scale is set to 0.95, as proposed by IQWiG (2016))"
      stepl1_str = " (lower boundary in HR scale is set to 0.85, as proposed by IQWiG (2016))"
    }
    if(type == "normal"){
      steps1_str = " (lower boundary in HR scale is set to 0, as proposed by Cohen (1988))"
      stepm1_str = " (lower boundary in HR scale is set to 0.5, as proposed Cohen (1988))"
      stepl1_str = " (lower boundary in HR scale is set to 0.8, as proposed Cohen (1988))"
    }
  }
  if(setting == "multiarm"){
    custom_further = "\\item{Strategy}{Strategy, 1: \"only best promising\" or 2: \"all promising\"}"
  }
  if(setting == "multiple"){
    if(type == "tte"){
      custom_further = paste0(
        "\\item{Strategy}{Strategy, 1: \"only best promising\" or 2: \"all promising\"}",
        "\\item{OP}{probability that one endpoint is significant}")
    }
  }
  return(paste0("The output of the function is a `data.frame` object containing the optimization results:
                 \\describe{",
                custom_further,
                "\\item{u}{maximal expected utility under the optimization constraints, i.e. the expected utility of the optimal sample size and threshold value}",  
                custom_threshold,
                "\\item{n2}{total sample size for phase II; rounded to the next even natural number}
                 \\item{n3}{total sample size for phase III; rounded to the next even natural number}
                 \\item{n}{total sample size in the program; n = n2 + n3}
                 \\item{K}{maximal costs of the program}
                 \\item{pgo}{probability to go to phase III}
                 \\item{sProg}{probability of a successful program}
                 \\item{sProg1}{probability of a successful program with \"small\" treatment effect in phase III", steps1_str, "}
                 \\item{sProg2}{probability of a successful program with \"medium\" treatment effect in phase III", stepm1_str,"}
                 \\item{sProg3}{probability of a successful program with \"large\" treatment effect in phase III", stepl1_str,"}
                 \\item{K2}{expected costs for phase II}
                 \\item{K3}{expected costs for phase III}}
                 and further input parameters. Taking `cat(comment())` of the
                 data frame lists the used optimization sequences, start and 
                 finish date of the optimization procedure."))
}