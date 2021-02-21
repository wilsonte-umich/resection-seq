
This repository contains all script and reference files
required to analyze data from the DSB resection sequencing
approach reported by Bazzano et al.

<link is pending>

The main file to execute the pipeline is 'process_sample.q',
which must be called by the q pipeline manager:

https://github.com/wilsonte-umich/q-pipeline-manager

Alternatively, individual scripts can be run directly
as long as all appropriate environment variables have been set.

The pipeline must be provided with appropriate DSB and control
allele reference sequences and files that describe them.
Appropriate files for all sequences in the original manuscript
are found in folder Bazzano_2021.

The pipeline must also be provided with file
ResectionMasterStrainTable.txt; once again, an appropriate file
is provided for the original paper in folder Bazzano_2021.

After pipeline steps 'align', 'collapse' and 'crosstab' have
been successfully run, data can be visualized in the R Shiny web tools
found in the '_server' directory. These web tools peform all
final data normalization steps and generate figures and
provide access to download data tables.

See the Methods section of the source paper for additional details. 

