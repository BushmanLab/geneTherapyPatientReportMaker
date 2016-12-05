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

# access patient metadata for GTSP

#' @param GTSP vector of GTSP 
#' @return df with columns: Trial, GTSP, Patient, Timepoint, CellType, FragMethod, VCN.
#' @note all GTSP in th vector should be unique
get_metadata_for_GTSP <- function(GTSP, db_group) {
    stopifnot(length(GTSP) == length(unique(GTSP)))
    GTSP = paste0("^", GTSP, "$", collapse="|")
    
    GTSPDBconn <- dbConnect(MySQL(), group=db_group)
    
    query = paste0("SELECT Trial, SpecimenAccNum, Patient, Timepoint, CellType,
                           SamplePrepMethod, VCN
                   FROM gtsp
                   WHERE SpecimenAccNum
                   REGEXP ", dbQuoteString(GTSPDBconn, GTSP), ";")
    
    sets <- dbGetQuery(GTSPDBconn, query)
    dbDisconnect(GTSPDBconn)
    
    names(sets) <- c("Trial", "GTSP", "Patient", "Timepoint", "CellType", "FragMethod", "VCN")
    sets    
}
