#Python code to visualize some stats for the pacas files
#Kenneth DeMonn
#20210201
#v07

# at the top of the file, before other imports
import warnings

warnings.filterwarnings('ignore')

# no warnings will be printed from now on

#imports matplotlib
import sys, subprocess
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd
from collections import Counter
import seaborn as sns
import math
import scipy.stats as stats
from scipy.stats import ttest_ind



filepath1 = str(sys.argv[1]) + '.fasta'
print('\nFilepath1: ', filepath1)


filepath2 = 'config.prm'
print('\nFilepath2: ', filepath2)

filepath3 = str(sys.argv[1])
print('\nFilepath3: ', filepath3)

#Window size
width = str(sys.argv[2])
print("\nWindow Size: ", width)

#Hardcoded to be 'out' at present
extension = 'out'


#Read in file to be processed, assign path to variables


#
#
# Section zero: call perl programs
#
#

#Call perl program that reads in *.fasta and *.out files and generates the *.csv and *.htm files

print("\n\nCalling perl program Pacas_FB to create *.csv and *.htm files")

subprocess.call(['perl', 'Pacas_FB.pl', '-d',  width, '-e', extension, *sys.argv[1:]])

#Call perl progfram that reads *.htm file and generates the config.prm file

print("\n\nCalling perl program Pacas_CB to create config.prm file")

subprocess.call(['perl', 'Pacas_CP.pl'])

print("\nDone!")



#These two lists hold the accession numbers from the fasta file to use as headers, and the actual sequences.

accessions = []
sequences = []

#Initialize
flanking = False
number_of_LCRs = 0



#
#
#   Section One: Read in data.
#
#
print("\n\nReading files.")


#Read in parameter file, assign values to variables

with open(filepath2) as fp2:
    print('\nFasta file: ', filepath1)
    #filepath2 = x_prm.rstrip()
    print('\nParameter file: ', filepath2)
    #filepath3 = x1.rstrip()
    #print('\nFilepath3: ', filepath3)

    for line in fp2:
        x = list(line.split())
        #print("\nX: ", x)

        #Python case statement
        if x[0] == 'Window_size:':
            windows_size = int(x[1])
            print("\nWindow_size: ", windows_size)

        if x[0] == 'Flanking:':
            if x[1].upper() == 'Y' or 'YES':
                flanking = True
            #print("\nFlanking: ", flanking)

        if x[0] == 'Number_of_LCRs:':
            number_of_LCRs = int(x[1])
            #print("\nNumber_of_LCRs: ", number_of_LCRs)

        if x[0] == "LCR_range1:":
            LCR_range1_lower = int(x[1])
            LCR_range1_upper = int(x[2])
            #print("\nLCR_range1: ", LCR_range1_lower, " - ", LCR_range1_upper)
        if x[0] == "LCR_range2:":
            LCR_range2_lower = int(x[1])
            LCR_range2_upper = int(x[2])
            #print("\nLCR_range2: ", LCR_range2_lower, " - ", LCR_range2_upper)
        if x[0] == "LCR_range3:":
            LCR_range3_lower = int(x[1])
            LCR_range3_upper = int(x[2])
            #print("\nLCR_range3: ", LCR_range3_lower, " - ", LCR_range3_upper)
        if x[0] == "LCR_range4:":
            LCR_range4_lower = int(x[1])
            LCR_range4_upper = int(x[2])
            #print("\nLCR_range4: ", LCR_range4_lower, " - ", LCR_range4_upper)
            
            
        if x[0] == "LCR_range5:":
            LCR_range5_lower = int(x[1])
            LCR_range5_upper = int(x[2])
            #print("\nLCR_range5: ", LCR_range5_lower, " - ", LCR_range5_upper)
        if x[0] == "LCR_range6:":
            LCR_range6_lower = int(x[1])
            LCR_range6_upper = int(x[2])
            #print("\nLCR_range6: ", LCR_range6_lower, " - ", LCR_range6_upper)
        if x[0] == "LCR_range7:":
            LCR_range7_lower = int(x[1])
            LCR_range7_upper = int(x[2])
            #print("\nLCR_range7: ", LCR_range7_lower, " - ", LCR_range7_upper)
        if x[0] == "LCR_range8:":
            LCR_range8_lower = int(x[1])
            LCR_range8_upper = int(x[2])
            #print("\nLCR_range8: ", LCR_range8_lower, " - ", LCR_range8_upper)



#Read in fasta file and place the accession numbers and sequences
#in their respective lists

with open(filepath1) as fp1:
    #print("\nFilepath1", fp1)
#   for cnt, line in enumerate(fp):
    for line in fp1:
       if line[0] == '>' :
           #append to accessions
           accessions.append(line.lstrip('>').rstrip('\n'))
       else:
           #append to sequences
           sequences.append(line.rstrip('\n'))

#print("\n\nAccessions: ", accessions)
#print("\n\nSequences: ", sequences)





#
#
#       Section Two: Calculate Consensus Sequence and Mismatches.
#
#
print("\n\nCalculating consensus and mismatches.")





#Create blank Pandas dataframe
df_Mismatch = pd.DataFrame(data=None, index=None, columns=['Read_In', 'Cumulative', 'Consensus', 'Mismatch', 'Temp'], dtype=None, copy=True)


#Loop thru sequences, reading in each line. Convert to list, and place in
#first column. Then concatenate to second column. To wit, each line is
#read and placed vertically in a column, it is appended to the second column
#to create a structure whose rows each correspond to a position
#on the sequence.

for i in range(len(sequences)):
    sequence = sequences[i]
    df_Mismatch['Read_In'] = list(sequence)
    df_Mismatch['Cumulative'] = df_Mismatch['Read_In'] + df_Mismatch['Cumulative'].astype(str)

#print("\nDataFrame after file read: \n\n", df_Mismatch)

#Remove 'nan's from end of second column.
df_Mismatch['Cumulative'] = df_Mismatch['Cumulative'].str.rstrip('nan')
#print("\nDataFrame after nans removed: \n\n", df_Mismatch)


#Mismatch function

def calculate_mismatches(Cumulative):
    mismatches = 0

    for i in Cumulative:
        for j in Cumulative:
            #print("\nI: ", i, "J: ", j)
            if i != j:
                mismatches +=1

    return mismatches /2


#Consensus function

def find_most_common_char(input):

    q = Counter(input).most_common(1)[0]
    return(q)


#Calculate the mismatches, place in dataframe
df_Mismatch['Mismatch'] = df_Mismatch.Cumulative.apply(calculate_mismatches)

#Calculate the consensus sequence, place in dataframe
df_Mismatch['Temp'] = df_Mismatch.Cumulative.apply(find_most_common_char)
#Returns the tuple for some reason. Temporarily stored to subsequently strip
df_Mismatch['Consensus'] = list(map(lambda x: x[0], df_Mismatch['Temp']))
#print("\nDataFrame after most common character: \n\n ", df_Mismatch)






#
#
#       Section Three: Calculate Entropies and Z-Scores
#
#       
#
print("\n\nCalculating Entropies and Z-scores.")



#Create empty dataframes
df_Entropies = pd.DataFrame()
df_AbsZScores = pd.DataFrame()



#Create necessary functions

def window(fseq,fn):
    alpha=[fseq[i:i+fn] for i in range(len(fseq)-(fn-1))]
    return alpha


def entropy(s):
     p, lns = Counter(s), float(len(s))
     return -sum( count/lns * math.log(count/lns, 2) for count in p.values()) 
     

def sliding_window(s, w):
    entropy_sequence = entropy(s)
    #print("\n\nEntropy of sequence: ", entropy_sequence)
    for i in s:
        substring = s[i:i + w]
        #print("\n\nSubstring: ", substring)




#Go thru the stored sequences, and for each one in turn, create a list of
#subranges of size windows_size

for i in range(len(sequences)):
    #print("\n\nLoop iteration outermost: ", i)
    sequence = sequences[i]
    #print("\n\nLoop iteration outermost: ", i, "Sequence: ", sequence)
    y = entropy(sequence)
    #print ("\n\nEntropy of entire sequence ", i,": ", y)

    #Creates a list x of sliding subsegments of length windows_size
    x = window(sequence, windows_size)
    #print("\n\nReturn from window function call: ", x)


    #Holds substrings and their entropies
    shannon_tuple = []
    #print("\nShannon_tuple: ", shannon_tuple)

    for k in range(len(x)):
        z = x[k], entropy(x[k])
        #Both values are kept mainly for possible diagnostics
        shannon_tuple.append(z)
        #print("\n\nLoop iteration: ", k)
    #print("\n\nShannon_tuple: ", shannon_tuple)

    #Import the tuple of substrings with entropies into dataframe
    df_Temp = pd.DataFrame(shannon_tuple, columns=['Substring', 'Entropies'])
    #Create a column of the absolute value of the zscore of the column of entropies of substrings
    df_Temp['AbsZScores'] = abs(stats.zscore(df_Temp['Entropies']))


    #Print the dataframe
    #print("\n\ndf_Temp: \n", df_Temp)

    #Assign the assession numbers as column headers, then import the columns
    #from the tempory dataframe as new columns in the respective permanent
    #dataframes
    column_header = accessions[i]
    #print("\n\nColumn_header: ", column_header, type(column_header))
    df_Entropies[column_header] = df_Temp['Entropies']





    df_AbsZScores[column_header] = df_Temp['AbsZScores']


#print("\n\ndf_Entropies: \n",df_Entropies)

df_Entropies.to_excel("Entropies.xlsx")

#print("\n\ndf_AbsZScores: \n",df_AbsZScores)

df_AbsZScores.to_excel("AbsZScores.xlsx")




#
#
#       Section Four: Student T Entropy Heatmaps
#
#
print("\n\nCreating Heatmaps")



df_Student_T = pd.DataFrame(index=accessions)
df_Student_T_LCR1 = pd.DataFrame(index=accessions)
df_Student_T_LCR2 = pd.DataFrame(index=accessions)
df_Student_T_LCR3 = pd.DataFrame(index=accessions)
df_Student_T_LCR4 = pd.DataFrame(index=accessions)
df_Student_T_LCR5 = pd.DataFrame(index=accessions)
df_Student_T_LCR6 = pd.DataFrame(index=accessions)
df_Student_T_LCR7 = pd.DataFrame(index=accessions)
df_Student_T_LCR8 = pd.DataFrame(index=accessions)


#Perform a pairwise two value student t test for the mean of each column
#of entropies

for j in range(len(accessions)):
    temp_Student_T =[]
    for k in range(len(accessions)):
        l = accessions[j]
        #print("\nL: ", l)
        m = accessions[k]
        #print("\nM: ", m)
        z1 = ttest_ind(df_Entropies[l], df_Entropies[m], nan_policy='omit')
        #if z1 == nan:
        #    z1 = 1
        #print("\n\nStudent-T of entropy column", l, " versus ", m, ":", z1)
        temp_Student_T.append(z1[1])
        #print("\n\ntemp_Student_T: ", temp_Student_T)
    df_Student_T[accessions[j]] = temp_Student_T
    #print("\n\ndf_Student_T: ",df_Student_T)

#df_Student.set_index(accession)
#print("\n\ndf_Student_T: ",df_Student_T)
df_Student_T.to_excel("Student_T_all.xlsx")

#Do each of the subranges in turn

if number_of_LCRs >= 1:
	df_Entropies_LCR1 = df_Entropies[LCR_range1_lower:LCR_range1_upper+1].copy()
	#print("\ndf_LCR1: ", df_Entropies_LCR1)
	df_Entropies_LCR1.to_csv("Entropies_LCR1.csv")

	#Perform a pairwise two value student t test for the mean each column of
	#entropies
	for j in range(len(accessions)):
		temp_Student_T =[]
		for k in range(len(accessions)):
			l = accessions[j]
			#print("\nL: ", l)
			m = accessions[k]
			#print("\nM: ", m)
			z1 = ttest_ind(df_Entropies_LCR1[l], df_Entropies_LCR1[m])
			#print("\n\nStudent-T of entropy column", l, " versus ", m, ":", z1)
			temp_Student_T.append(z1[1])
			#print("\n\ntemp_Student_T: ", temp_Student_T)
		df_Student_T_LCR1[accessions[j]] = temp_Student_T
		#print("\n\ndf_Student_T_LCR1: ",df_Student_T_LCR1)

	df_Student_T_LCR1.to_excel("Student_T_LCR1.xlsx")  
	#print("\n\ndf_Student_T_LCR1: ",df_Student_T_LCR1)
	#This step was for the case of gaps returning a NaN student-t.
	#Left it as is because the blank square give useful information
	#df_Student_T_LCR1 = df_Student_T_LCR1.replace(np.nan,1)
	#print("\nReplace nans \n")
	#print("\n\ndf_Student_T_LCR1: ",df_Student_T_LCR1)
     

if number_of_LCRs >= 2:
	df_Entropies_LCR2 = df_Entropies[LCR_range2_lower:LCR_range2_upper+1].copy()
	#print("\ndf_LCR2: ", df_Entropies_LCR2)
	df_Entropies_LCR2.to_csv("Entropies_LCR2.csv")
    
	#Perform a pairwise two value student t test for the mean each column of
	#entropies
	for j in range(len(accessions)):
		temp_Student_T =[]
		for k in range(len(accessions)):
			l = accessions[j]
			#print("\nL: ", l)
			m = accessions[k]
			#print("\nM: ", m)
			z1 = ttest_ind(df_Entropies_LCR2[l], df_Entropies_LCR2[m])
			#print("\n\nStudent-T of entropy column", l, " versus ", m, ":", z1)
			temp_Student_T.append(z1[1])
			#print("\n\ntemp_Student_T: ", temp_Student_T)
		df_Student_T_LCR2[accessions[j]] = temp_Student_T
		#print("\n\ndf_Student_T: ",df_Student_T)

	df_Student_T_LCR2.to_excel("Student_T_LCR2.xlsx")
	#print("\n\ndf_Student_T_LCR2: ",df_Student_T_LCR2)

if number_of_LCRs >= 3:
	df_Entropies_LCR3 = df_Entropies[LCR_range3_lower:LCR_range3_upper+1].copy()
	#print("\ndf_LCR3: ", df_Entropies_LCR3)
	df_Entropies_LCR3.to_csv("Entropies_LCR3.csv")
    
	#Perform a pairwise two value student t test for the mean each column of
	#entropies
	for j in range(len(accessions)):
		temp_Student_T =[]
		for k in range(len(accessions)):
			l = accessions[j]
			#print("\nL: ", l)
			m = accessions[k]
			#print("\nM: ", m)
			z1 = ttest_ind(df_Entropies_LCR3[l], df_Entropies_LCR3[m])
			#print("\n\nStudent-T of entropy column", l, " versus ", m, ":", z1)
			temp_Student_T.append(z1[1])
			#print("\n\ntemp_Student_T: ", temp_Student_T)
		df_Student_T_LCR3[accessions[j]] = temp_Student_T
		#print("\n\ndf_Student_T: ",df_Student_T)

	df_Student_T_LCR3.to_excel("Student_T_LCR3.xlsx")
	#print("\n\ndf_Student_T_LCR3: ",df_Student_T_LCR3)

if number_of_LCRs >= 4:
	df_Entropies_LCR4 = df_Entropies[LCR_range4_lower:LCR_range4_upper+1].copy()
	#print("\ndf_LCR4: ", df_Entropies_LCR4)
	df_Entropies_LCR4.to_csv("Entropies_LCR4.csv")
    
	#Perform a pairwise two value student t test for the mean each column of
	#entropies
	for j in range(len(accessions)):
		temp_Student_T =[]
		for k in range(len(accessions)):
			l = accessions[j]
			#print("\nL: ", l)
			m = accessions[k]
			#print("\nM: ", m)
			z1 = ttest_ind(df_Entropies_LCR4[l], df_Entropies_LCR4[m])
			#print("\n\nStudent-T of entropy column", l, " versus ", m, ":", z1)
			temp_Student_T.append(z1[1])
			#print("\n\ntemp_Student_T: ", temp_Student_T)
		df_Student_T_LCR4[accessions[j]] = temp_Student_T
		#print("\n\ndf_Student_T: ",df_Student_T)

	df_Student_T_LCR4.to_excel("Student_T_LCR4.xlsx")
	#print("\n\ndf_Student_T_LCR4: ",df_Student_T_LCR4)

    

if number_of_LCRs >= 5:
	df_Entropies_LCR5 = df_Entropies[LCR_range5_lower:LCR_range5_upper+1].copy()
	#print("\ndf_LCR5: ", df_Entropies_LCR5)
	df_Entropies_LCR5.to_csv("Entropies_LCR5.csv")
    
	#Perform a pairwise two value student t test for the mean each column of
	#entropies
	for j in range(len(accessions)):
		temp_Student_T =[]
		for k in range(len(accessions)):
			l = accessions[j]
			#print("\nL: ", l)
			m = accessions[k]
			#print("\nM: ", m)
			z1 = ttest_ind(df_Entropies_LCR5[l], df_Entropies_LCR5[m])
			#print("\n\nStudent-T of entropy column", l, " versus ", m, ":", z1)
			temp_Student_T.append(z1[1])
			#print("\n\ntemp_Student_T: ", temp_Student_T)
		df_Student_T_LCR5[accessions[j]] = temp_Student_T
		#print("\n\ndf_Student_T: ",df_Student_T)

	df_Student_T_LCR5.to_excel("Student_T_LCR5.xlsx")
	#print("\n\ndf_Student_T_LCR5: ",df_Student_T_LCR5)
    
    
if number_of_LCRs >= 6:
	df_Entropies_LCR6 = df_Entropies[LCR_range6_lower:LCR_range6_upper+1].copy()
	#print("\ndf_LCR6: ", df_Entropies_LCR6)
	df_Entropies_LCR6.to_csv("Entropies_LCR6.csv")
    
	#Perform a pairwise two value student t test for the mean each column of
	#entropies
	for j in range(len(accessions)):
		temp_Student_T =[]
		for k in range(len(accessions)):
			l = accessions[j]
			#print("\nL: ", l)
			m = accessions[k]
			#print("\nM: ", m)
			z1 = ttest_ind(df_Entropies_LCR6[l], df_Entropies_LCR6[m])
			#print("\n\nStudent-T of entropy column", l, " versus ", m, ":", z1)
			temp_Student_T.append(z1[1])
			#print("\n\ntemp_Student_T: ", temp_Student_T)
		df_Student_T_LCR6[accessions[j]] = temp_Student_T
		#print("\n\ndf_Student_T: ",df_Student_T)

	df_Student_T_LCR6.to_excel("Student_T_LCR6.xlsx")
	#print("\n\ndf_Student_T_LCR6: ",df_Student_T_LCR6)


if number_of_LCRs >= 7:
	df_Entropies_LCR7 = df_Entropies[LCR_range7_lower:LCR_range7_upper+1].copy()
	#print("\ndf_LCR7 : ", df_Entropies_LCR7)
	df_Entropies_LCR7.to_csv("Entropies_LCR7.csv")

	#Perform a pairwise two value student t test for the mean each column of
	#entropies
	for j in range(len(accessions)):
		temp_Student_T =[]
		for k in range(len(accessions)):
			l = accessions[j]
			#print("\nL: ", l)
			m = accessions[k]
			#print("\nM: ", m)
			z1 = ttest_ind(df_Entropies_LCR7[l], df_Entropies_LCR7[m])
			#print("\n\nStudent-T of entropy column", l, " versus ", m, ":", z1)
			temp_Student_T.append(z1[1])
			#print("\n\ntemp_Student_T: ", temp_Student_T)
		df_Student_T_LCR7[accessions[j]] = temp_Student_T
		#print("\n\ndf_Student_T: ",df_Student_T)

	df_Student_T_LCR7.to_excel("Student_T_LCR7.xlsx")
	#print("\n\ndf_Student_T_LCR7: ",df_Student_T_LCR7)


if number_of_LCRs >= 8:
	df_Entropies_LCR8 = df_Entropies[LCR_range8_lower:LCR_range8_upper+1].copy()
	#print("\ndf_LCR8: ", df_Entropies_LCR8)
	df_Entropies_LCR8.to_csv("Entropies_LCR8.csv")
	#entropies    

	#Perform a pairwise two value student t test for the mean each column of
	#entropies
	for j in range(len(accessions)):
		temp_Student_T =[]
		for k in range(len(accessions)):
			l = accessions[j]
			#print("\nL: ", l)
			m = accessions[k]
			#print("\nM: ", m)
			z1 = ttest_ind(df_Entropies_LCR8[l], df_Entropies_LCR8[m])
			#print("\n\nStudent-T of entropy column", l, " versus ", m, ":", z1)
			temp_Student_T.append(z1[1])
			#print("\n\ntemp_Student_T: ", temp_Student_T)
		df_Student_T_LCR8[accessions[j]] = temp_Student_T
		#print("\n\ndf_Student_T: ",df_Student_T)

	df_Student_T_LCR8.to_excel("Student_T_LCR8.xlsx")
	#print("\n\ndf_Student_T_LCR8: ",df_Student_T_LCR8)





#
#
#       Print PNGs
#
#
print("\n\nWriting out PNG files.")



#
#
#       Section Four: Printing Mismatches
#
#

# Originally plotted out as png, but didn't yield much
# visually useful info, so written out as excel file

df_Mismatch.to_excel("Mismatches.xlsx")


#consensus = list(df_Mismatch['Consensus'])
#z = df_Mismatch.info
#print("\nnMismatch structure: ", z)
#print("\n\nConsensus = ", consensus)
#Without flanks
#values = list(df_Mismatch['Mismatch'])
#print("\nValues: ", values)
#position = range(1, len(values) +1)
#print("\nPosition: ", position)
#plt.figure(figsize=(30,10))
#plt.title("Mismatches at each site position")
#sns.color_palette('GnBu')
#ax = sns.scatterplot(x = position, y = values, marker = "+")
#ax = df_Mismatch['Mismatch'].plot.bar(cmap = 'BrBG')
#ax.set_xlabel("Site Position")
#ax.set_ylabel("Number of Mismatches")
#plt.set_xlabel("Site Position")
#plt.set_ylabel("Number of Mismatches")


#Color in the LCR ranges

#if number_of_LCRs >= 1:
#	ax.axvspan(LCR_range1_lower, LCR_range1_upper, facecolor='green', alpha=0.4)
#if number_of_LCRs >= 2:
#	ax.axvspan(LCR_range2_lower, LCR_range2_upper, facecolor='green', alpha=0.4)
#if number_of_LCRs >= 3:
#	ax.axvspan(LCR_range3_lower, LCR_range3_upper, facecolor='green', alpha=0.4)
#if number_of_LCRs >= 4:
#	ax.axvspan(LCR_range4_lower, LCR_range4_upper, facecolor='green', alpha=0.4)
#if number_of_LCRs >= 5:
#	ax.axvspan(LCR_range5_lower, LCR_range5_upper, facecolor='green', alpha=0.4)
#if number_of_LCRs >= 6:
#	ax.axvspan(LCR_range6_lower, LCR_range6_upper, facecolor='green', alpha=0.4)
#if number_of_LCRs >= 7:
#	ax.axvspan(LCR_range7_lower, LCR_range7_upper, facecolor='green', alpha=0.4)
#if number_of_LCRs >= 8:
#	ax.axvspan(LCR_range8_lower, LCR_range8_upper, facecolor='green', alpha=0.4)

#outputfilename = str(filepath3) + "_Mismatch.png"
#print("\nFilename: ", outputfilename)
#plt.savefig(outputfilename)

#plt.show()
#plt.clf()






#       Plot Entropy Heatmaps

plt.figure(figsize=(10,10))
az = sns.heatmap(df_Student_T,  annot = True, cmap = 'GnBu', vmin = 0, vmax = 1)
#cmap = 'GnBu',
plt.title("Heat map of p-value for Student-T cross-comparison of entropies")

outputfilename = str(filepath3) + "_Student_T_All.png"
print("\nFilename: ", outputfilename)
az.figure.tight_layout()
plt.savefig(outputfilename)
plt.clf()



if number_of_LCRs >= 1:
    plt.figure(figsize=(10,10))
    az1 = sns.heatmap(df_Student_T_LCR1, annot = True, cmap = 'GnBu', vmin = 0, vmax = 1)
    plt.title("Heat map of p-value for Student-T LCR1")
    outputfilename = str(filepath3) + "_Student_T_LCR1.png"
    print("\nFilename: ", outputfilename)
    az1.figure.tight_layout()
    plt.savefig(outputfilename)
    plt.clf()

if number_of_LCRs >= 2:
    plt.figure(figsize=(10,10))
    az2 = sns.heatmap(df_Student_T_LCR2, annot = True, cmap = 'GnBu', vmin = 0, vmax = 1)
    plt.title("Heat map of p-value for Student-T LCR2")
    outputfilename = str(filepath3) + "_Student_T_LCR2.png"
    print("\nFilename: ", outputfilename)
    az2.figure.tight_layout()
    plt.savefig(outputfilename)
    plt.clf()
    
    
    
if number_of_LCRs >= 3:
    plt.figure(figsize=(10,10))
    az3 = sns.heatmap(df_Student_T_LCR3, annot = True, cmap = 'GnBu', vmin = 0, vmax = 1)
    plt.title("Heat map of p-value for Student-T LCR3")
    outputfilename = str(filepath3) + "_Student_T_LCR3.png"
    print("\nFilename: ", outputfilename)
    az3.figure.tight_layout()
    plt.savefig(outputfilename)
    plt.clf()

if number_of_LCRs >= 4:
    plt.figure(figsize=(10,10))
    az4 = sns.heatmap(df_Student_T_LCR4, annot = True, cmap = 'GnBu', vmin = 0, vmax = 1)
    plt.title("Heat map of p-value for Student-T LCR4")
    outputfilename = str(filepath3) + "_Student_T_LCR4.png"
    print("\nFilename: ", outputfilename)
    az4.figure.tight_layout()
    plt.savefig(outputfilename)
    plt.clf()


if number_of_LCRs >= 5:
    plt.figure(figsize=(10,10))
    az5 = sns.heatmap(df_Student_T_LCR5, annot = True, cmap = 'GnBu', vmin = 0, vmax = 1)
    plt.title("Heat map of p-value for Student-T LCR5")
    outputfilename = str(filepath3) + "_Student_T_LCR5.png"
    print("\nFilename: ", outputfilename)
    az5.figure.tight_layout()
    plt.savefig(outputfilename)
    plt.clf()

if number_of_LCRs >= 6:
    plt.figure(figsize=(10,10))
    az6= sns.heatmap(df_Student_T_LCR6, annot = True, cmap = 'GnBu', vmin = 0, vmax = 1)
    plt.title("Heat map of p-value for Student-T LCR6")
    outputfilename = str(filepath3) + "_Student_T_LCR6.png"
    print("\nFilename: ", outputfilename)
    az6.figure.tight_layout()
    plt.savefig(outputfilename)
    plt.clf()

if number_of_LCRs >= 7:
    plt.figure(figsize=(10,10))
    az7= sns.heatmap(df_Student_T_LCR7, annot = True, cmap = 'GnBu', vmin = 0, vmax = 1)
    plt.title("Heat map of p-value for Student-T LCR7")
    outputfilename = str(filepath3) + "_Student_T_LCR7.png"
    print("\nFilename: ", outputfilename)
    az7.figure.tight_layout()
    plt.savefig(outputfilename)
    plt.clf()

if number_of_LCRs >= 8:
    plt.figure(figsize=(10,10))
    az8= sns.heatmap(df_Student_T_LCR8, annot = True, cmap = 'GnBu', vmin = 0, vmax = 1)
    plt.title("Heat map of p-value for Student-T LCR8")
    outputfilename = str(filepath3) + "_Student_T_LCR8.png"
    print("\nFilename: ", outputfilename)
    az8.figure.tight_layout()
    plt.savefig(outputfilename)
    plt.clf()






#
# Plot Z-Scores
#





# Without flanks

plt.clf()
plt.figure(figsize=(40,40))

ay1 = df_AbsZScores.plot.area(cmap = 'GnBu')
ay1.set_xlabel("Site Position")
ay1.set_ylabel("Absolute Z Score For Sliding Window Entropy")

if number_of_LCRs >= 1:
	ay1.axvspan(LCR_range1_lower, LCR_range1_upper, facecolor='green', alpha=0.4)
if number_of_LCRs >= 2:
	ay1.axvspan(LCR_range2_lower, LCR_range2_upper, facecolor='green', alpha=0.4)
if number_of_LCRs >= 3:
	ay1.axvspan(LCR_range3_lower, LCR_range3_upper, facecolor='green', alpha=0.4)
if number_of_LCRs >= 4:
	ay1.axvspan(LCR_range4_lower, LCR_range4_upper, facecolor='green', alpha=0.4)
if number_of_LCRs >= 5:
	ay1.axvspan(LCR_range5_lower, LCR_range5_upper, facecolor='green', alpha=0.4)
if number_of_LCRs >= 6:
	ay1.axvspan(LCR_range6_lower, LCR_range6_upper, facecolor='green', alpha=0.4)
if number_of_LCRs >= 7:
	ay1.axvspan(LCR_range7_lower, LCR_range7_upper, facecolor='green', alpha=0.4)
if number_of_LCRs >= 8:
	ay1.axvspan(LCR_range8_lower, LCR_range8_upper, facecolor='green', alpha=0.4)    
    
    

outputfilename = str(filepath3) + "_Z_score_entropy.png"
print("\nFilename: ", outputfilename)
plt.title("Input file: " + str(filepath1))
ay1.figure.tight_layout()
plt.savefig(outputfilename)




# With flanks


plt.clf()
plt.figure(figsize=(40,40))
plt.title(str(filepath1) )
#outputfilename = str(filepath3) + "_Z_score_entropy_flanks.png"


ay2 = df_AbsZScores.plot.area(cmap = 'GnBu')

ay2.set_xlabel("Site Position")
ay2.set_ylabel("Absolute Z Score For Sliding Window Entropy")


if number_of_LCRs >= 1:
	ay2.axvspan(LCR_range1_lower - windows_size, LCR_range1_lower, facecolor='green', alpha=0.2)
	ay2.axvspan(LCR_range1_lower, LCR_range1_upper, facecolor='green', alpha=0.4)
	ay2.axvspan(LCR_range1_upper, LCR_range1_upper + windows_size, facecolor='green', alpha=0.2)

if number_of_LCRs >= 2:
	ay2.axvspan(LCR_range2_lower - windows_size, LCR_range2_lower, facecolor='green', alpha=0.2)
	ay2.axvspan(LCR_range2_lower, LCR_range2_upper, facecolor='green', alpha=0.4)
	ay2.axvspan(LCR_range2_upper, LCR_range2_upper + windows_size, facecolor='green', alpha=0.2)

if number_of_LCRs >= 3:
	ay2.axvspan(LCR_range3_lower - windows_size, LCR_range3_lower, facecolor='green', alpha=0.2)
	ay2.axvspan(LCR_range3_lower, LCR_range3_upper, facecolor='green', alpha=0.4)
	ay2.axvspan(LCR_range3_upper, LCR_range3_upper + windows_size, facecolor='green', alpha=0.2)

if number_of_LCRs >= 4:
	ay2.axvspan(LCR_range4_lower - windows_size, LCR_range4_lower, facecolor='green', alpha=0.2)
	ay2.axvspan(LCR_range4_lower, LCR_range4_upper, facecolor='green', alpha=0.4)
	ay2.axvspan(LCR_range4_upper, LCR_range4_upper + windows_size, facecolor='green', alpha=0.2)


if number_of_LCRs >= 5:
	ay2.axvspan(LCR_range5_lower - windows_size, LCR_range5_lower, facecolor='green', alpha=0.2)
	ay2.axvspan(LCR_range5_lower, LCR_range5_upper, facecolor='green', alpha=0.4)
	ay2.axvspan(LCR_range5_upper, LCR_range5_upper + windows_size, facecolor='green', alpha=0.2)

if number_of_LCRs >= 6:
	ay2.axvspan(LCR_range6_lower - windows_size, LCR_range6_lower, facecolor='green', alpha=0.2)
	ay2.axvspan(LCR_range6_lower, LCR_range6_upper, facecolor='green', alpha=0.4)
	ay2.axvspan(LCR_range6_upper, LCR_range6_upper + windows_size, facecolor='green', alpha=0.2)

if number_of_LCRs >= 7:
	ay2.axvspan(LCR_range7_lower - windows_size, LCR_range7_lower, facecolor='green', alpha=0.2)
	ay2.axvspan(LCR_range7_lower, LCR_range7_upper, facecolor='green', alpha=0.4)
	ay2.axvspan(LCR_range7_upper, LCR_range7_upper + windows_size, facecolor='green', alpha=0.2)

if number_of_LCRs >= 8:
	ay2.axvspan(LCR_range8_lower - windows_size, LCR_range8_lower, facecolor='green', alpha=0.2)
	ay2.axvspan(LCR_range8_lower, LCR_range8_upper, facecolor='green', alpha=0.4)
	ay2.axvspan(LCR_range8_upper, LCR_range8_upper + windows_size, facecolor='green', alpha=0.2)



outputfilename = str(filepath3) + "_Z_score_entropy_flanks.png"
print("\nFilename: ", outputfilename)
plt.title("Input file: " + str(filepath1))
ay2.figure.tight_layout()
plt.savefig(outputfilename)






















