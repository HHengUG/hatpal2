% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clusterAPA.R
\name{ClusterAPA}
\alias{ClusterAPA}
\title{Cluster the two-dimensional coordinate files (SA and RC/SA) to define the APA clusters}
\usage{
ClusterAPA(
  SA_prebed,
  uni_prebed,
  out_dir,
  gapdist_cut = 5,
  APAdist_cut = 25,
  RC_cutoff = 3,
  SA_weight = 3
)
}
\arguments{
\item{SA_prebed}{SA file of one strand generated from ExtractAPA()}

\item{uni_prebed}{rc_uni or SA file of one strand generated from ExtractAPA()}

\item{out_dir}{The directory where the output files are stored}

\item{gapdist_cut}{Threshold for filtering out isolated nodes (default 5)}

\item{APAdist_cut}{Threshold for the gap between APA (default 25)}

\item{RC_cutoff}{High observation nodes (default 3)}

\item{SA_weight}{The weight of SA file (default 3)}
}
\value{
NULL, writes the output "APA_clusters.out" in the out_dir.
}
\description{
Cluster the two-dimensional coordinate files (SA and RC/SA) to define the APA clusters
}
\examples{
ClusterAPA("sa.prebed", "rc_uni.prebed", "output/")
}
