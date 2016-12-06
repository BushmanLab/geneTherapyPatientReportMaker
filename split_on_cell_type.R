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

library(argparse)
library(assertthat)
library(dplyr)

parser <- ArgumentParser(description=
  paste0('Split sample-gtsp file on cell type ',
  'and write separate files for each cell type. File Format: celltype_SAMPLE-GTSP-FILENAME'))
parser$add_argument("sample_gtsp", nargs=1, default="sampleName_GTSP.csv")

arguments <- parser$parse_args()

sample_gtsp_filename <- arguments$sample_gtsp
assert_that(file.exists(sample_gtsp_filename))

sampleName_GTSP <- read.csv(sample_gtsp_filename)

get_cell_types <- function(sampleName_GTSP) {
  return(unique(sampleName_GTSP$celltype))
}

cell_types <- get_cell_types(sampleName_GTSP)
print('Cell types present:')
print(cell_types)

out <- sapply(cell_types, function(cell_type) {
    message(cell_type)
    write.csv(filter(sampleName_GTSP, celltype==cell_type), 
        file=paste0(cell_type, '_', sample_gtsp_filename),
        quote=FALSE, row.names = FALSE
    )          
})


