# PACAS

PACAS consists of three files: Pacas.py, Pacas_FB.pl, and Pacas_CP.pl. In the current version, these all need to present in the same directory as the data files. The input files consist of a fasta file that contains the accession numbers for protein sequences, and the sequences themselves. The sequences need to be properly aligned, on a single line, and of equal length. It also requires an *.out file identically named to the *.fasta file, which indicates the low complexity regions of interest. These can be created as the output of programs such as seg, or manually.

To run, place the scripts and input files in the same directory, and type “python Pacas.py <name of file pair> <length of sliding windows>” . For example, using the AMA1.fasta and AMA1.out files given as an example and desiring a length of ten, enter “python Pacas.py AMA1 10”. 

The output consists of an htm file <file name>_comparison.htm which details each of the LCR regions for each sequence in a pair-wise comparison to every other sequence, listing mismatches and gaps for each pair and the adjoining regions. There is also a *.csv file with the same information. There is also an Excel file that shows all of sequences with the consensus sequence and the total number of mismatches at every site position.

In terms of graphic output, there are two sets of *.png files. In one, the Shannon entropy of the entire sequence is computed, and then a heat map is produced from the student-t comparison of each sequence. This allows for a quick visual comparison of which sequence is showing the most variance. There are also heatmaps for each of the LCRs individually. Note: when comparing a gap region to another the heat map shows a blank square, since comparing the entropy of null strings is not meaningful.

The other output is a pair of *.png files showing the entropy value of each of the sliding window sections along the length of the sequence, plotted as a stacked line graph with the LCRs highlighted, and with and without the flanks highlighted. This too allows for visual inspections of regions with a greater than average degree of variability.

Additionally, the contents of the dataframes that underly the graphs are written out as Excel or *.csv files, for those whose wish to perform their own numeric analysis.

