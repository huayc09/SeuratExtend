% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SearchDatabase.R
\name{SearchDatabase}
\alias{SearchDatabase}
\title{FUNCTION_TITLE}
\usage{
SearchDatabase(
  item,
  type = c("gene", "SetID", "SetName"),
  database = c("GO", "Reactome"),
  spe = getOption("spe"),
  export.to.data.frame = F,
  n.min = 1,
  n.max = Inf,
  only.end.terms = F
)
}
\arguments{
\item{item}{PARAM_DESCRIPTION}

\item{type}{PARAM_DESCRIPTION, Default: c("gene", "SetID", "SetName")}

\item{database}{PARAM_DESCRIPTION, Default: c("GO", "Reactome")}

\item{spe}{PARAM_DESCRIPTION, Default: getOption("spe")}

\item{export.to.data.frame}{PARAM_DESCRIPTION, Default: F}

\item{n.min}{PARAM_DESCRIPTION, Default: 1}

\item{n.max}{PARAM_DESCRIPTION, Default: Inf}

\item{only.end.terms}{PARAM_DESCRIPTION, Default: F}
}
\value{
OUTPUT_DESCRIPTION
}
\description{
FUNCTION_DESCRIPTION
}
\details{
DETAILS
}
\examples{
\dontrun{
if(interactive()){
 #EXAMPLE1
 }
}
}
