\name{PhenoPro}
\alias{PhenoPro}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	This is PhenoPro
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
PhenoPro(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
RawData <- data(sds)
PhenoPro(RawData, x = "Fs", y = "LIGHT", label = "R4DX.F",
	defaultLabel = "D", block = "RowRun", orderby = "Time.n",
	method = "SVM",
	step = 1, width = 5,
	cvNumber = 100,
	testBlockProp = 0.2,
	visualization = FALSE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
