#### geneTherapyPatientReportMaker
For a specific patient, report integration near oncogenes and potentially expanded clones across cell types and multiple time points. Abundance by default is based on unique fragment length.
	
#### Input
A csv file such as `sampleName_GTSP.csv` to describe the replicates and the samples.
```
head sampleName_GTSP.csv

sampleName,GTSP
GTSP0308-1,GTSP0308
GTSP0308-2,GTSP0308
GTSP0308-3,GTSP0308
GTSP0308-4,GTSP0308
GTSP0309-1,GTSP0309
GTSP0309-2,GTSP0309

#or 
Rscript path/to/check_patient_GTSP.R pFR03

sampleName,GTSP,patient
GTSP0308-1,GTSP0308,pFR03
GTSP0308-2,GTSP0308,pFR03
GTSP0308-3,GTSP0308,pFR03
```

* only `sampleName`, `GTSP` columns are necessary and the rest are ignored,
* the `GTSPxxxx` names must correspond to the same patient,
* the `GTSPxxxx` names must be in the `specimen_management.gtsp` database,
* all sites should be computed based on one reference genome.

An option is available to generate reports without having the integration sites loaded into the database. The structure of the data through needs to be in the terminal format from `intSiteCaller`. For this method, the allSites.RData and multihitData.RData files need to be present for each sample to be processed. An example is below:
```
Data/
    GTSP0308-1/
    	allSites.RData
    	multihitData.RData
    	...
    GTSP0308-2/
    	allSites.RData
    	multihitData.RData
    	...
    ...
```

#### Output
`$trial.$patient.$today.html`

#### Code example
- 1. `check_patient_gtsp.R`: get available datasets, generate a csv file for a patient 
```
 Rscript path/to/check_patient_GTSP.R                    #get all processed samples
 Rscript path/to/check_patient_GTSP.R pFR03              #get data sets for patient pFR03 and output to csv format
 Rscript path/to/check_patient_GTSP.R pFR03 > pFR03.csv  #get data sets for patient pFR03 and output to tmp.csv
```

- 2. `makeGeneTherapyPatientReport.R`: generate report for a patient from the csv file 
```
Rscript makeGeneTherapyPatientReport.R                     #read in sampleName_GTSP.csv by default
Rscript path/to/makeGeneTherapyPatientReport.R pFR03.csv   #generated above
Rscript path/to/makeGeneTherapyPatientReport.R pFR03.csv -s #determine abundance by sonicLength package (Berry, C. 2012)
Rscript path/to/makeGeneTherapyPatientReport.R pFR02.csv -d path/to/Data/ #Generate report from processed data set
```

- 3 `check_gtsp_patient.R`: get trial and patient information for the GTSPxxxx folders
```
Rscript path/to/makeGeneTherapyPatientReport.R                         #check current folder
Rscript path/to/makeGeneTherapyPatientReport.R  ~/Frances/run20150505  #check a run folder
```

Reference genome is specified by `--ref_genome` or `-r` option with default of `hg18`:
```
Rscript path/to/makeGeneTherapyPatientReport.R  ~/Frances/run20150505  --ref_genome hg19
```

To connect to databases .my.cnf should be present in ~ and group can be changed with `--sites_group` option
with default of `intsites_miseq` for integration sites DB:
```
Rscript path/to/makeGeneTherapyPatientReport.R  ~/Frances/run20150505  --sites_group test_db
```

Metadata for GTSP is held in specimen_management DB and can be changed with `--gtsp_group` option
with default "specimen_management":
```
Rscript path/to/makeGeneTherapyPatientReport.R  ~/Frances/run20150505  --gtsp_group gtsp_group_in_my_cnf
```



#### Database config file location

config file should be in home directory and called .my.cnf,
e.g. ~/.my.cnf

The .my.cnf format is as follows:

```
[GROUP_NAME]
user=YYYYYYY
password=XXXXXX
host=microbYYYY.med.upenn.edu
port=3309
database=intsites_miseq
```

#### Dependencies

...

#### Testing

Run in the R console:

```bash
library(testthat)
test_dir(".")
```
