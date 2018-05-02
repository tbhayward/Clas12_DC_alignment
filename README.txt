Currently using 3 groovy scripts:

elastic_peak.groovy
Arguments: 
0 - directory with hipo files to analyze, 
1 - number of files to analyze (if not given assumes all), 
2 - minimum bin in W^2 (GeV^2) (if not given assumes 0.65 GeV^2), 
3 - maximum bin W^2 (GeV^2) (if not given assumes 1.2 GeV^2), 
4 - number of bins (if not given assumes 250)

residuals.groovy
Arguments: 
0 - directory with hipo files to analyze, 
1 - number of files to analyze (if not given assumes all), 
2 - minimum theta of particle track to allow (if not given, assumes 0 degrees)
3 - maximum theta of particle track to allow (if not given, assumes 90 degrees)
4 - int 1 (sector vs superlayer plots) or int 2 (sector vs layer plots)

update_alignment_table.groovy 
Arguments: