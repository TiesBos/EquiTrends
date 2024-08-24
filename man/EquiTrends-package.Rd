\name{EquiTrends-package}
\alias{EquiTrends-package}
\alias{EquiTrends}
\docType{package}
\title{
  Equivalence Testing for Pre-Trends in Difference-in-Differences Designs
}
\description{
  Testing for parallel trends is crucial in the Difference-in-Difference framework. EquiTrends is an R package for equivalence testing in the context of Difference-in-Differences estimation. It allows users to test if pre-treatment trends in the treated group are “equivalent” to those in the control group. Here, “equivalence” means that rejection of the null hypothesis implies that a function of the pre-treatment placebo effects (maximum absolute, average or root mean squared value) does not exceed a pre-specified threshold below which trend differences are considered negligible. The package is based on the theory
developed in Dette & Schumann (2024) <doi: 10.1080/07350015.2024.2308121>. 
}
\details{
  The package contains the functions maxEquivTest to perform the testing procedure surrounding the maximum placebo coefficient (see equation (3.1) of Dette & Schumann (2024)), meanEquivTest to perform the testing procedure surrounding the mean placebo coefficient (see equation (3.2) of Dette & Schumann (2024)) and rmsEquivTest to perform the testing procedure surrounding the root mean squared placebo coefficient (see equation (3.3) and (3.4) of Dette & Schumann (2024)). Furthermore, the package contains the function sim_paneldata to simulate a paneldataset for such testing purposes.
}
\author{
Maintainer: Ties Bos <tc.bos@student.maastrichtuniversity.nl>
}
\references{
  Dette H., & Schumann M. (2024). “Testing for Equivalence of Pre-Trends
in Difference-in-Differences Estimation.” *Journal of Business &
Economic Statistics*, 1–13. DOI:
[10.1080/07350015.2024.2308121](https://doi.org/10.1080/07350015.2024.2308121)
}
\keyword{ package }


