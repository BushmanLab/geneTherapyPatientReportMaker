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


