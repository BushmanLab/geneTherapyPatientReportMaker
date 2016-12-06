#    This source code file is a component of the larger INSPIIRED genomic analysis software package.
#    Copyright (C) 2016 Frederic Bushman
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

codeDir <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly=FALSE), value=T)))
if( length(codeDir)!=1 ) codeDir <- list.files(path="~", pattern="geneTherapyPatientReportMaker$", recursive=TRUE, include.dirs=TRUE, full.names=TRUE)
stopifnot(file.exists(file.path(codeDir, "GTSPreport.css")))
stopifnot(file.exists(file.path(codeDir, "GTSPreport.Rmd")))

#### load up require packages + objects #### 
library("RMySQL") #also loads DBI
library("plyr")
library("dplyr")
library("stringr")
library("markdown")
library("knitr")
library("PubMedWordcloud")
library("hiAnnotator")
library("ggplot2")
library("reldist")
library("sonicLength")
library("reshape2")
library("scales")
library("intSiteRetriever") 
## intSiteRetriever package was installed from github as follows
## git clone https://github.com/BushmanLab/intSiteRetriever.git
## cd intSiteRetriever
## git checkout remove_multihitClusterID
## R
## devtools::document()
## devtools::install()

source(file.path(codeDir, "utilities.R"))
source(file.path(codeDir, "specimen_management.R"))
source(file.path(codeDir, "estimatedAbundance.R"))
source(file.path(codeDir, "dereplicateSites.R"))
source(file.path(codeDir, "standardizeSites.R"))
source(file.path(codeDir, "read_site_totals.R"))
source(file.path(codeDir, "populationInfo.R"))
source(file.path(codeDir, "abundanceFilteringUtils.R"))

load("debug.RData")

knit(file.path(codeDir, "GTSPreport.Rmd"), output=mdfile)
markdownToHTML(mdfile, htmlfile, extensions=c('tables'),
               options=c(markdownHTMLOptions(defaults=T), "toc"),
               stylesheet=file.path(codeDir, "GTSPreport.css") )
cmd <- paste("Rscript ~/geneTherapyPatientReportMaker/printReportToPdf.R", htmlfile)
system(cmd)
pdffile <- sub("html$", "pdf", htmlfile)

system(paste("scp", htmlfile, "microb215:Sites/share/GTSPReports/bushmanlab/"))
system(paste("scp", pdffile, "microb215:Sites/share/GTSPReports/bushmanlab/"))
              
