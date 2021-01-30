#!/usr/bin/perl #-w
use strict;
use Getopt::Long;
use List::Util qw[min max];

print "\n";
print "PACAS- Pairwise Comparison of Aligned Subsequences   Version 7.0 - Oakland University\n";
print "Usage: PACAS.pl -d DesiredFlankingRegionLength -e FileExtensionToUse\n";
print "------------------------------------------------------------------------------\n";
print "\n";

#Declare the perl command line flags/options we want to allow ...
my ($FlankingRegionLength, $FileExtensionToUse, $CombinedFileName, $Help, $Brief, $dataDirectory, $tabOutput);

GetOptions(	"d=i", \$FlankingRegionLength,
			"e=s", \$FileExtensionToUse,
#			"c:s", \$CombinedFileName,
			"h", \$Help,
#			"b", \$Brief,
			"dir:s", \$dataDirectory,
			"t", \$tabOutput )
or die("Error in command line arguments!\nUse -h to display help ...\n");

if ($Help)
{
	print "-h Command Line Option detected ...\n";
	QuickStart();
	exit;
}
if (not ($FlankingRegionLength))
{
	print "-d DesiredFlankingRegionLength Command Line Parameter missing!\n-h To see Help.\n";
	exit;
}
if (! ($FileExtensionToUse))
{
	print "-e FileExtensionToUse Command Line Parameter missing!\n-h To see Help.\n";
	exit;
}

#if ($Brief)
#{
#	print "-b Command Line Option detected ...\n";
#}
#if ($CombinedFileName)
#{
#	print "-c $CombinedFileName Command Line Option detected ...\n";
#}
if ($tabOutput)
{
	print "-t Command Line Option detected ...\n";
}

if ($dataDirectory)
{
	chdir "$dataDirectory";
	opendir (DIR, $dataDirectory) or die ("Cannot open $dataDirectory");
	print "Data Directory -=>$dataDirectory\n"
}

#Looks like you cannot choose '0' as a flanking input
my $isFlanking;

if ($FlankingRegionLength eq "0") {
	$isFlanking = 'No';
} else {
	$isFlanking = 'Yes';
}

open my $fhConfig, '>', 'config.prm' or die $!;
print $fhConfig "Window_size: $FlankingRegionLength";
print $fhConfig "\nFlanking: $isFlanking";
close $fhConfig; 



my @DIR =<*$dataDirectory>;

print "\nDesired Flanking Region Length: $FlankingRegionLength\n";
print "         File Extension to Use: $FileExtensionToUse\n\n";

#Determine what type of LCR/SEGMENT input file to use: HP or .txt ...
my $FileType = $FileExtensionToUse;

print "Looking for $FileExtensionToUse files in the folder to process ...\n";
my @files = <*$FileType>;
if (scalar(@files) > 0)
{
	print "Found! Processing $FileType files ...\n\n";
}
else
{
	print "Not Found! Aborting ...\n";
	exit;
}

#my $CombinedFileOpenFlag = 0; #Flag to open the combined file once ...

foreach my $file (@files) {

	print $file . "\n";

	my ($offset, $InputFilePrefix, $InputFile);
	$offset = index($file,"."); #Look for the first . in file name ...
	$InputFilePrefix = substr($file, 0, $offset);
	$InputFile = $InputFilePrefix . ".fasta";

	my @species_names = ();
	my @species_protein = ();

	#Read associated FASTA file into arrays ...
	#Sample identifier line ...
	#>pberghei.fasta|gi|PBANKA_060540|ref|PBANKA_060540|_protein__|_length_677
	#Followed with Sequence Data line(s) ...
	#
	Read_Into_Arrays($InputFile, \@species_names, \@species_protein); #Pass by value and reference ...

	my ($OutputFile);
	my ($OutputFile0);
	my ($OutputFile1);
	my ($OutputFile2);

	#if ($tabOutput)
	#{
	#	$OutputFile0 = $InputFile."_Comparison.tab";
	#}
	#else
	#{
	#	$OutputFile1 = $InputFile."_Comparison.htm";
	#

	#}




	if ($tabOutput)
	{
		$OutputFile0 = $InputFile."_Comparison.htm";
		$OutputFile1 = $InputFile."_Comparison.tab";
		open FILE00, ">$OutputFile0" or die "Can't open output file: $OutputFile0\n";
		open FILE01, ">$OutputFile1" or die "Can't open output file: $OutputFile1\n";

	}
	else
	{
		$OutputFile0 = $InputFile."_Comparison.htm";
		$OutputFile2 = $InputFile."_Comparison.csv";
		open FILE00, ">$OutputFile0" or die "Can't open output file: $OutputFile0\n";
		open FILE02, ">$OutputFile2" or die "Can't open output file: $OutputFile2\n";
	}





	if ($tabOutput)
	{
		print FILE00 "<h1>PACAS - Reporting Matches/Mismatches with a desired flanking region length of $FlankingRegionLength</h1>\n";
		print FILE01 "PACAS - Reporting Matches/Mismatches with a desired flanking region length of $FlankingRegionLength\n";
#		print FILE02 "PACAS - Reporting Matches/Mismatches with a desired flanking region length of $FlankingRegionLength\n";
	}
	else
	{
		print FILE00 "<h1>PACAS - Reporting Matches/Mismatches with a desired flanking region length of $FlankingRegionLength</h1>\n";
#		print FILE01 "PACAS - Reporting Matches/Mismatches with a desired flanking region length of $FlankingRegionLength\n";
		print FILE02 "PACAS - Reporting Matches/Mismatches with a desired flanking region length of $FlankingRegionLength\n";
	}


	print "Processing $InputFile ...\n";

	if ($tabOutput)
	{
		print FILE00 "<h3>Segment File: $file</h3>\n";
		print FILE00 "<h3>FASTA File: $InputFile</h3>\n";
		print FILE00 "<font face='courier'><table border=3>";
		print FILE01 "Segment File: $file\n";
		print FILE01 "FASTA File: $InputFile\n";
#		print FILE02 "Segment File: $file\n";
#		print FILE02 "FASTA File: $InputFile\n";
	}
	else
	{
		print FILE00 "<h3>Segment File: $file</h3>\n";
		print FILE00 "<h3>FASTA File: $InputFile</h3>\n";
		print FILE00 "<font face='courier'><table border=3>";
#		print FILE01 "Segment File: $file\n";
#		print FILE01 "FASTA File: $InputFile\n";
		print FILE02 "Segment File: $file\n";
		print FILE02 "FASTA File: $InputFile\n";
	}

#	if (not $Brief) #Verbose, original output ...

	if (1) #Verbose, original output ...
	{
		if ($tabOutput)
		{
			print FILE00 "<tr><th>File Handle</th><th>Accession #1</th><th>ID #1</th><th>Accession #2</th><th>Segment Start</th><th>Segment End</th><th>Window Start</th><th>Window End</th><th>Window #1 with Flanks</th><th>Window #2 with Flanks</th><th>Non-Gap Length #1</th><th>Window Length</th><th>Matches</th><th>Mismatches</th><th>Gaps</th><th>Gap Matches</th><th>Left Flank Matches</th><th>Left Flank Mismatches</th><th>Left Flank Gaps</th><th>Left Flank Gap Matches</th><th>Right Flank Matches</th><th>Right Flank Mismatches</th><th>Right Flank Gaps</th><th>Right Flank Gap Matches</th><th>Gapped Flanks</th><th>Rest of sequence comparison Matches</th><th>Rest of sequence comparison Mismatches</th><th>Rest of sequence comparison Gaps</th><th>Rest of sequence comparison Gap Matches</th></tr>\n";

			print FILE01 "File Handle\tAccession #1\tID #1\tAccession #2\tSegment Start\tSegment End\tWindow Start\tWindow End\tWindow #1 with Flanks\tWindow #2 with Flanks\tNon-Gap Length #1\tWindow Length\tMatches\tMismatches\tGaps\tGap Matches\tLeft Flank Matches\tLeft Flank Mismatches\tLeft Flank Gaps\tLeft Flank Gap Matches\tRight Flank Matches\tRight Flank Mismatches\tRight Flank Gaps\tRight Flank Gap Matches\tGapped Flanks\tRest of sequence comparison Matches\tRest of sequence comparison Mismatches\tRest of sequence comparison Gaps\tRest of sequence comparison Gap Matches\n";


#			print FILE02 "File Handle,Accession #1,ID #1,Accession #2,Segment Start,Segment End,Aligned Start,Aligned End,Window Start,Window End,Window #1 with Flanks,Window #2 with Flanks,Non-Gap Length #1,Window Length,Matches,Mismatches,Gaps,Gap Matches,Left Flank Matches,Left Flank Mismatches,Left Flank Gaps,Left Flank Gap Matches,Right Flank Matches,Right Flank Mismatches,Right Flank Gaps,Right Flank Gap Matches,Alignment Corrupted,Comparison A Matches,Comparison A Mismatches,Comparison A #Gaps,Comparison A Gap Matches,Comparison B Matches,Comparison B Mismatches,Comparison B Gaps,Comparison B Gap Matches,Comparison C Matches,Comparison C Mismatches,Comparison C Gaps,Comparison C Gap Matches,Comparison D Matches,Comparison D Mismatches,Comparison D Gaps,Comparison A Gap Matches\n";


		}
		else
		{
			print FILE00 "<tr><th>File Handle</th><th>Accession #1</th><th>ID #1</th><th>Accession #2</th><th>Segment Start</th><th>Segment End</th><th>Window Start</th><th>Window End</th><th>Window #1 with Flanks</th><th>Window #2 with Flanks</th><th>Non-Gap Length #1</th><th>Window Length</th><th>Matches</th><th>Mismatches</th><th>Gaps</th><th>Gap Matches</th><th>Left Flank Matches</th><th>Left Flank Mismatches</th><th>Left Flank Gaps</th><th>Left Flank Gap Matches</th><th>Right Flank Matches</th><th>Right Flank Mismatches</th><th>Right Flank Gaps</th><th>Right Flank Gap Matches</th><th>Gapped Flanks</th><th>Rest of sequence comparison Matches</th><th>Rest of sequence comparison Mismatches</th><th>Rest of sequence comparison Gaps</th><th>Rest of sequence comparison Gap Matches</th></tr>\n";


#			print FILE01 "File Handle\tAccession #1\tID #1\tAccession #2\tSegment Start\tSegment End\tAligned Start\tAligned End\tWindow Start\tWindow End\tWindow #1 with Flanks\tWindow #2 with Flanks\tNon-Gap Length #1\tWindow Length\tMatches\tMismatches\tGaps\tGap Matches\tLeft Flank Matches\tLeft Flank Mismatches\tLeft Flank Gaps\tLeft Flank Gap Matches\tRight Flank Matches\tRight Flank Mismatches\tRight Flank Gaps\tRight Flank Gap Matches\tAlignment Corrupted\tComparison A Matches\tComparison A #Mismatches\tComparison A Gaps\tComparison A Gap Matches\tComparison B Matches\tComparison B Mismatches\tComparison B Gaps\tComparison B Gap Matches\tComparison C Matches\tComparison C Mismatches\tComparison C Gaps\tComparison C Gap Matches\tComparison D Matches\tComparison D Mismatches\tComparison D Gaps\tComparison A Gap Matches\n";



			print FILE02 "File Handle,Accession #1,ID #1,Accession #2,Segment Start,Segment End,Window Start,Window End,Window #1 with Flanks,Window #2 with Flanks,Non-Gap Length #1,Window Length,Matches,Mismatches,Gaps,Gap Matches,Left Flank Matches,Left Flank Mismatches,Left Flank Gaps,Left Flank Gap Matches,Right Flank Matches,Right Flank Mismatches,Right Flank Gaps,Right Flank Gap Matches,Gapped Flanks,Rest of sequence comparison Matches,Rest of sequence comparison Mismatches,Rest of sequence comparison Gaps,Rest of sequence comparison Gap Matches,\n";


#
#			}
		}
	}


	#Read LCR/SEGMENT file into parallel arrays Species, Genes, Location, and MatchString ...
	my (@Species, @Genes, @StartLocation, @EndLocation, @MatchString);
	my (@AlignedStart, @AlignedEnd, @LCRStart, @LCREnd, @SkipFlag);
	@Species = ();
	@Genes = ();
	@StartLocation = ();
	@EndLocation = ();
	@MatchString = (); #Will ignore it since it does not show gaps (if any) ...
	@AlignedStart = ();
	@AlignedEnd = ();
	@LCRStart = ();
	@LCREnd = ();
	@SkipFlag = (); #Will use to flag entries printed as a part of merged intervals ...

	open FILE, $file or die "Can't read source file: $file\n";

	my $FileFormatFirstLineFlag = 0; #Flag to determine the format of LCR/SEGMENT file (HP or txt) from the first line ...
	$FileType = "TBD";

	my ($Specie2Compare, $StartLocation2Compare, $EndLocation2Compare);
	my ($LCR2Compare, $LCR2CompareLong, $LCR2CompareStart, $LCR2CompareEnd); #AlignedStart and AlignedEnd ...

	L1: while (my $line = <FILE>)
	{
		chomp $line;

		next L1 if ($line eq ""); #Discard blank lines ...

		if ($FileFormatFirstLineFlag == 0)
		{
			if ($line =~ /^>/) #Can look for comparison as well ...
			{
				$FileType = ".out"; #HP
			}
			else
			{
				$FileType = ".txt"; #.txt
			}
			$FileFormatFirstLineFlag = 1;
		}

		if ($FileType eq ".out")
		{
			if ($line =~ /^>/)
			{
				$offset = index($line,".fasta"); #Look for .fasta ...
				push @Species, substr($line, 1, $offset-1);

				$Specie2Compare = substr($line, 1, $offset-1);

				$offset = index($line,"gi");
				my $offset2 = index($line,"ref");
				push @Genes, substr($line, $offset+3, $offset2 - $offset -4);

				$offset = index($line,"(");
				$offset2 = index($line,"-");
				my $offset3 = index($line,")");
				push @StartLocation, substr($line, $offset+1, $offset2 - $offset -1);
				push @EndLocation, substr($line, $offset2+1, $offset3 - $offset2 -1);

				$StartLocation2Compare = substr($line, $offset+1, $offset2 - $offset -1);
				$EndLocation2Compare = substr($line, $offset2+1, $offset3 - $offset2 -1);

				$line = <FILE>;
				chomp $line;
				$line =~ tr/A-Za-z/a-zA-Z/; #Convert to uppercase ...
				push @MatchString, $line;
			}
		}

		if ($FileType eq ".txt")
		{
			my @varSplit = split(",", $line);

			push @Species, trim($varSplit[0]);
			push @Genes, trim($varSplit[1]);
			push @StartLocation, trim($varSplit[2]);
			push @EndLocation, trim($varSplit[3]);
			push @MatchString, $line; #Ignored ...

			$Specie2Compare = trim($varSplit[0]);
			$StartLocation2Compare = trim($varSplit[2]);
			$EndLocation2Compare = trim($varSplit[3]);
		}
		#print ("EndLocation2Compare: $EndLocation2Compare\n");

		#Find AlignedStart and AlignedEnd ...
		#Find its protein sequence in @species_protein read from fasta file ...
		my $Protein2Compare = ();

		L2: for (my $j=0; $j < scalar(@species_names); $j++)
		{
			my ($varOffset, $varSpecie);

			$varOffset = index(@species_names[$j],".fasta"); #Look for .fasta ...
			$varSpecie = substr(@species_names[$j], 1, $varOffset-1);

			if (lc $Specie2Compare eq lc $varSpecie) #Compare in lowercase to avoid issues ...
			{
				$Protein2Compare = @species_protein[$j];
				last L2;
			}
		}

		#Find its LCR/SEGMENT aligned start and end ...
		$LCR2Compare = ();
		my $ProteinCounter = 0;
		print ("\nProtein2Compare334: $Protein2Compare\n");
		L3: for (my $j=0; $j < length($Protein2Compare); $j++)
		{
			if (substr($Protein2Compare, $j, 1) ne "&")
			{
				$ProteinCounter++;
				if (($ProteinCounter >= $StartLocation2Compare) && ($ProteinCounter <= $EndLocation2Compare))
				{
					$LCR2Compare = $LCR2Compare . substr($Protein2Compare, $j, 1);
				}

				if ($ProteinCounter == $StartLocation2Compare)
				{
					$LCR2CompareStart = $j;
				}

				if ($ProteinCounter == $EndLocation2Compare)
				{
					$LCR2CompareEnd = $j;
					#print ("LCR2CompareEnd: $LCR2CompareEnd\n");
					last L3;
				}
			}
		}
		push @AlignedStart, $LCR2CompareStart;
		push @AlignedEnd, $LCR2CompareEnd;
	} #L1 ...

	close FILE;

	#Try to determine the largest window of comparison when overlap exists in AlignedStart's and AlignedEnd's ...
	my @IntervalStart = ();
	my @IntervalEnd = ();
	my $Temp;

	for (my $i=0; $i < scalar(@AlignedStart); $i++)
	{
		push @IntervalStart, $AlignedStart[$i];
		push @IntervalEnd, $AlignedEnd[$i];
	}
	#print ("IntervalEnd arr: @IntervalEnd\n");

	#Sort Intervals ...
	for (my $i=0; $i < scalar(@IntervalStart)-1; $i++)
	{
		for (my $j=$i+1; $j < scalar(@IntervalStart); $j++)
		{
			if ($IntervalStart[$i] > $IntervalStart[$j])
			{
				$Temp = $IntervalStart[$i];
				$IntervalStart[$i] = $IntervalStart[$j];
				$IntervalStart[$j] = $Temp;

				$Temp = $IntervalEnd[$i];
				$IntervalEnd[$i] = $IntervalEnd[$j];
				$IntervalEnd[$j] = $Temp;
			}
			else
			{
				if ($IntervalStart[$i] == $IntervalStart[$j])
				{
					if ($IntervalEnd[$i] > $IntervalEnd[$j])
					{
						$Temp = $IntervalEnd[$i];
						$IntervalEnd[$i] = $IntervalEnd[$j];
						$IntervalEnd[$j] = $Temp;
					}
				}
			}
		}
	}

	#Merge Overlapping Intervals ...
	my @MergedStart = ();
	my @MergedEnd = ();
	my $MergedIndex = 0;

	for (my $i=0; $i < scalar(@IntervalStart); $i++)
	{
		if ( ($MergedIndex != 0) && ($MergedEnd[$MergedIndex - 1] >= $IntervalStart[$i]) && ($IntervalEnd[$i] >= $MergedStart[$MergedIndex - 1]) )
		{
			#The two intervals overlap, merge them ...
	        $MergedEnd[$MergedIndex - 1] = max($MergedEnd[$MergedIndex - 1], $IntervalEnd[$i]);
			$MergedStart[$MergedIndex - 1] = min($MergedStart[$MergedIndex - 1], $IntervalStart[$i]);
		}
		else
		{
			#Add new interval ...
			$MergedStart[$MergedIndex] = $IntervalStart[$i];
			$MergedEnd[$MergedIndex] = $IntervalEnd[$i];
			$MergedIndex++;
		}
	}

	#Now, assign the overlapping window to use ...
	for (my $i=0; $i < scalar(@AlignedStart); $i++)
	{
		my ($X, $Y);
		my ($WindowStart, $WindowEnd);

		$X = $AlignedStart[$i];
		$Y = $AlignedEnd[$i];
		#print ("MergedIndex: $MergedIndex\n");
		#Find merged interval $X, $Y belong to ...
		for (my $j=0; $j < $MergedIndex; $j++)
		{
			if ( ($X >= $MergedStart [$j]) && ($Y <= $MergedEnd[$j]) )
			{
				$WindowStart = $MergedStart [$j];
				$WindowEnd = $MergedEnd[$j];

				last;
			}
		}
		print ("\nWindowEnd: $WindowEnd\n");
		$LCRStart[$i] = $WindowStart;
		$LCREnd[$i] = $WindowEnd;
	}

	#Display the windows ...
	# for (my $i=0; $i < scalar(@LCRStart); $i++)
	# {
		# print "$AlignedStart[$i]-$AlignedEnd[$i] and $LCRStart[$i]--$LCREnd[$i]\n";
	# }
	# exit;

	#Ok. Prepare for printing one row for multiple intervals that are merged ...
	#

	#Sort our parallel arrays ... Getting more inefficient ...
	#First by Species, then by StartLocation ...

	for (my $i=0; $i < scalar(@Species)-1; $i++)
	{
		for (my $j=$i+1; $j < scalar(@Species); $j++)
		{
			if ($Species[$i] gt $Species[$j])
			{
				$Temp = $Species[$i];
				$Species[$i] = $Species[$j];
				$Species[$j] = $Temp;

				$Temp = $Genes[$i];
				$Genes[$i] = $Genes[$j];
				$Genes[$j] = $Temp;

				$Temp = $StartLocation[$i];
				$StartLocation[$i] = $StartLocation[$j];
				$StartLocation[$j] = $Temp;

				$Temp = $EndLocation[$i];
				$EndLocation[$i] = $EndLocation[$j];
				$EndLocation[$j] = $Temp;

				$Temp = $AlignedStart[$i];
				$AlignedStart[$i] = $AlignedStart[$j];
				$AlignedStart[$j] = $Temp;

				$Temp = $AlignedEnd[$i];
				$AlignedEnd[$i] = $AlignedEnd[$j];
				$AlignedEnd[$j] = $Temp;

				$Temp = $LCRStart[$i];
				$LCRStart[$i] = $LCRStart[$j];
				$LCRStart[$j] = $Temp;

				$Temp = $LCREnd[$i];
				$LCREnd[$i] = $LCREnd[$j];
				$LCREnd[$j] = $Temp;

				$Temp = $MatchString[$i];
				$MatchString[$i] = $MatchString[$j];
				$MatchString[$j] = $Temp;
			}
			else
			{
				if ( ($Species[$i] eq $Species[$j]) && ($StartLocation[$i] > $StartLocation[$j]) )
				{
					$Temp = $Species[$i];
					$Species[$i] = $Species[$j];
					$Species[$j] = $Temp;

					$Temp = $Genes[$i];
					$Genes[$i] = $Genes[$j];
					$Genes[$j] = $Temp;

					$Temp = $StartLocation[$i];
					$StartLocation[$i] = $StartLocation[$j];
					$StartLocation[$j] = $Temp;

					$Temp = $EndLocation[$i];
					$EndLocation[$i] = $EndLocation[$j];
					$EndLocation[$j] = $Temp;

					$Temp = $AlignedStart[$i];
					$AlignedStart[$i] = $AlignedStart[$j];
					$AlignedStart[$j] = $Temp;

					$Temp = $AlignedEnd[$i];
					$AlignedEnd[$i] = $AlignedEnd[$j];
					$AlignedEnd[$j] = $Temp;

					$Temp = $LCRStart[$i];
					$LCRStart[$i] = $LCRStart[$j];
					$LCRStart[$j] = $Temp;

					$Temp = $LCREnd[$i];
					$LCREnd[$i] = $LCREnd[$j];
					$LCREnd[$j] = $Temp;

					$Temp = $MatchString[$i];
					$MatchString[$i] = $MatchString[$j];
					$MatchString[$j] = $Temp;
				}
			}
		}
	}

	#Now, loop through to determine which duplicate rows to skip ...
	my ($strSpecie, $intWStart, $intWEnd, $strStart, $strEnd, $strAStart, $strAEnd, $j);
	my ($Temp1);

	for (my $i=0; $i < scalar(@Species); $i++)
	{
		if ($SkipFlag[$i] eq "")
		{
			$strSpecie = $Species[$i];
			$intWStart = $LCRStart[$i];
			$intWEnd = $LCREnd[$i];
			$strStart = $StartLocation[$i];
			$strEnd = $EndLocation[$i];
			$strAStart = $AlignedStart[$i]+1; #To accomodate indexing by 1 instead of 0 for print purposes ...
			$strAEnd = $AlignedEnd[$i]+1;

			$j = $i + 1;

			while ($j < scalar(@Species))

				{
					if ( ($Species[$j] eq $strSpecie) && ($LCRStart[$j] == $intWStart ) && ($LCREnd[$j] == $intWEnd ) )
					{
						$strStart = $strStart . "," . $StartLocation[$j];
						$strEnd = $strEnd . "," . $EndLocation[$j];

						$Temp1 = $AlignedStart[$j] + 1;
						$strAStart = $strAStart . "," . $Temp1; #To accomodate indexing by 1 instead of 0 for print purposes ...
						$Temp1 = $AlignedEnd[$j] + 1;
						$strAEnd = $strAEnd . "," . $Temp1; #To accomodate indexing by 1 instead of 0 for print purposes ...

						$SkipFlag[$j] = "Skip!";
					}

					if ($Species[$j] gt $strSpecie)
					{
						last;
					}

					$j++;
				}

			$SkipFlag[$i] = $strStart . " " . $strEnd . " " . $strAStart . " " . $strAEnd;
		}
	}
	# for (my $i=0; $i < scalar(@Species); $i++)
	# {
		# print "$SkipFlag[$i] \n";
	# }

	# exit;

	#Sort our parallel arrays AGAIN .. Getting more inefficient ...
	#First by LCRStart, then by Speces ...

	for (my $i=0; $i < scalar(@Species)-1; $i++)
	{
		for (my $j=$i+1; $j < scalar(@Species); $j++)
		{
			if ($LCRStart[$i] > $LCRStart[$j])
			{
				$Temp = $Species[$i];
				$Species[$i] = $Species[$j];
				$Species[$j] = $Temp;

				$Temp = $Genes[$i];
				$Genes[$i] = $Genes[$j];
				$Genes[$j] = $Temp;

				$Temp = $StartLocation[$i];
				$StartLocation[$i] = $StartLocation[$j];
				$StartLocation[$j] = $Temp;

				$Temp = $EndLocation[$i];
				$EndLocation[$i] = $EndLocation[$j];
				$EndLocation[$j] = $Temp;

				$Temp = $AlignedStart[$i];
				$AlignedStart[$i] = $AlignedStart[$j];
				$AlignedStart[$j] = $Temp;

				$Temp = $AlignedEnd[$i];
				$AlignedEnd[$i] = $AlignedEnd[$j];
				$AlignedEnd[$j] = $Temp;

				$Temp = $LCRStart[$i];
				$LCRStart[$i] = $LCRStart[$j];
				$LCRStart[$j] = $Temp;

				$Temp = $LCREnd[$i];
				$LCREnd[$i] = $LCREnd[$j];
				$LCREnd[$j] = $Temp;

				$Temp = $MatchString[$i];
				$MatchString[$i] = $MatchString[$j];
				$MatchString[$j] = $Temp;

				$Temp = $SkipFlag[$i];
				$SkipFlag[$i] = $SkipFlag[$j];
				$SkipFlag[$j] = $Temp;
			}
			else
			{
				if ( ($LCRStart[$i] == $LCRStart[$j]) && ($Species[$i] gt $Species[$j]) )
				{
					$Temp = $Species[$i];
					$Species[$i] = $Species[$j];
					$Species[$j] = $Temp;

					$Temp = $Genes[$i];
					$Genes[$i] = $Genes[$j];
					$Genes[$j] = $Temp;

					$Temp = $StartLocation[$i];
					$StartLocation[$i] = $StartLocation[$j];
					$StartLocation[$j] = $Temp;

					$Temp = $EndLocation[$i];
					$EndLocation[$i] = $EndLocation[$j];
					$EndLocation[$j] = $Temp;

					$Temp = $AlignedStart[$i];
					$AlignedStart[$i] = $AlignedStart[$j];
					$AlignedStart[$j] = $Temp;

					$Temp = $AlignedEnd[$i];
					$AlignedEnd[$i] = $AlignedEnd[$j];
					$AlignedEnd[$j] = $Temp;

					$Temp = $LCRStart[$i];
					$LCRStart[$i] = $LCRStart[$j];
					$LCRStart[$j] = $Temp;

					$Temp = $LCREnd[$i];
					$LCREnd[$i] = $LCREnd[$j];
					$LCREnd[$j] = $Temp;

					$Temp = $MatchString[$i];
					$MatchString[$i] = $MatchString[$j];
					$MatchString[$j] = $Temp;

					$Temp = $SkipFlag[$i];
					$SkipFlag[$i] = $SkipFlag[$j];
					$SkipFlag[$j] = $Temp;
				}
			}
		}
	}

	#Now, start LCR comparisons ... Skip duplicates ...
	L4: for (my $i=0; $i < scalar(@Species); $i++)
	{
		next L4 if ($SkipFlag[$i] eq "Skip!");

		my $Specie2Compare = $Species[$i];
		my $Gene2Compare = $Genes[$i];

		#Find its protein sequence in @species_protein read from fasta file ...
		my $Protein2Compare = ();

		for (my $j=0; $j < scalar(@species_names); $j++)
		{
			my ($varOffset, $varSpecie);

			$varOffset = index(@species_names[$j],".fasta"); #Look for .fasta ...
			$varSpecie = substr(@species_names[$j], 1, $varOffset-1);

			if (lc $Specie2Compare eq lc $varSpecie) #Compare in lowercase to avoid issues ...
			{
				$Protein2Compare = @species_protein[$j];
				print ("Protein2Compare722: $Protein2Compare\n");
				last;
			}
		}

		#Find its LCR ...
		my ($LCR2Compare, $LCR2CompareLong, $LCR2CompareStart, $LCR2CompareEnd, $LCR2WFlanks, $LCR2QueryLength);
		my $ProteinCounter = 0;

		$LCR2CompareStart = $LCRStart[$i];
		$LCR2CompareEnd = $LCREnd[$i];
		print ("LCR2CompareStart: $LCR2CompareStart\n");
		print ("LCR2CompareEnd: $LCR2CompareEnd\n");
		print ("\nProtein2Compare735: $Protein2Compare\n");
		$LCR2CompareLong = substr($Protein2Compare, $LCR2CompareStart, ($LCR2CompareEnd - $LCR2CompareStart) +1);
		print ("lcr2comparelong: $LCR2CompareLong\n");
		my $LCRLength2Compare = length($LCR2CompareLong);

		my $ProtenLength2Compare = 0;
		for (my $j=0; $j < length($Protein2Compare); $j++)
		{
			if (substr($Protein2Compare, $j, 1) ne "-")
			{
				$ProtenLength2Compare++;
			}
		}

		my ($LeftFlank2Compare, $RightFlank2Compare, $LeftStart, $RightEnd);

		$LeftStart = $LCR2CompareStart - $FlankingRegionLength;
		if ($LeftStart < 0)
		{
			$LeftStart = 0;
		}
		$LeftFlank2Compare = substr($Protein2Compare, $LeftStart, ($LCR2CompareStart - $LeftStart));

		$RightEnd = $LCR2CompareEnd + $FlankingRegionLength;
		if ($RightEnd >= length($Protein2Compare))
		{
			$RightEnd  = length($Protein2Compare) - 1;
		}
		$RightFlank2Compare = substr($Protein2Compare, $LCR2CompareEnd+1 , ($RightEnd - $LCR2CompareEnd));

		#With flanks ...
		if ($tabOutput)
		{
			$LCR2WFlanks = $LeftFlank2Compare . "(" . $LCR2CompareLong . ")" . $RightFlank2Compare;
		}
		else
		{
			$LCR2WFlanks = $LeftFlank2Compare . "<b>" . $LCR2CompareLong . "</b>" . $RightFlank2Compare;
		}

		#Query.Length ...
		$LCR2QueryLength = TrimLeadingAndTrailingLength($LCR2CompareLong);

		####################
		#Loop through .fasta file for actual LCR comparisons ...
		####################
		for (my $j=0; $j < scalar(@species_names); $j++)
		{
			my ($MySpecie2Compare, $MyGene2Compare, $MyProtein2Compare, $MyLCR2Compare, $MyLCR2CompareLong, $MyLCR2WFlanks, $MyLCR2SubjectLength);
			my ($MyLeftFlank2Compare, $MyRightFlank2Compare);

			$offset = index($species_names[$j],".fasta"); #Look for .fasta ...
			$MySpecie2Compare = substr($species_names[$j], 1, $offset-1);

			$offset = index($species_names[$j],"gi");
			my $offset2 = index($species_names[$j],"ref");
			$MyGene2Compare = substr($species_names[$j], $offset+3, $offset2 - $offset -4);

			#Ignore if same Specie ...
			if (lc $MySpecie2Compare ne lc $Specie2Compare)
			{

				#if ($Flag == 1)
				if (1 == 1)
				{
					#Find its protein sequence in @species_protein read from fasta file ...

					for (my $k=0; $k < scalar(@species_names); $k++)
					{
						my ($varOffset, $varSpecie);
						$varOffset = index(@species_names[$k],".fasta"); #Look for .fasta ...
						$varSpecie = substr(@species_names[$k], 1, $varOffset-1);

						if (lc $MySpecie2Compare eq lc $varSpecie)
						{
							$MyProtein2Compare = @species_protein[$k];
							last;
						}
					}

					#Find its LCR ... Using AlignedStart and AlignedEnd we found ...
					#Potential error if substring is looking outside the range ...

					$MyLCR2CompareLong = substr($MyProtein2Compare, $LCR2CompareStart, ($LCR2CompareEnd - $LCR2CompareStart) +1);

					my $FindMatches= Find_Matches( $LCR2CompareLong, $MyLCR2CompareLong );
					my $FindMismatches= Find_MisMatches( $LCR2CompareLong, $MyLCR2CompareLong );
					my $FindGaps= Find_Gaps ( $LCR2CompareLong, $MyLCR2CompareLong);
					my $FindGapMatches= Find_GapMatches ( $LCR2CompareLong, $MyLCR2CompareLong );

					#Determine Flanks ... same potential error ...
					$MyLeftFlank2Compare = substr($MyProtein2Compare, $LeftStart, ($LCR2CompareStart - $LeftStart));

					$MyRightFlank2Compare = substr($MyProtein2Compare, $LCR2CompareEnd+1 , ($RightEnd - $LCR2CompareEnd));

					my $FlankLeftMatches= Find_Matches( $LeftFlank2Compare, $MyLeftFlank2Compare );
					my $FlankLeftMismatches= Find_MisMatches( $LeftFlank2Compare, $MyLeftFlank2Compare );
					my $FlankLeftGaps= Find_Gaps ( $LeftFlank2Compare, $MyLeftFlank2Compare );
					my $FlankLeftGapMatches= Find_GapMatches ( $LeftFlank2Compare, $MyLeftFlank2Compare );

					my $FlankRightMatches= Find_Matches( $RightFlank2Compare, $MyRightFlank2Compare );
					my $FlankRightMismatches= Find_MisMatches( $RightFlank2Compare, $MyRightFlank2Compare );
					my $FlankRightGaps= Find_Gaps ( $RightFlank2Compare, $MyRightFlank2Compare );
					my $FlankRightGapMatches= Find_GapMatches ( $RightFlank2Compare, $MyRightFlank2Compare );
					#Check for Alignment Corrupted ...
					my $AlignmentCorrupted;
					if (($FlankLeftGaps > 0) or ($FlankRightGaps > 0))
					{
						$AlignmentCorrupted = "YES"
					}
					else
					{
						$AlignmentCorrupted = "NO"
					}

					#With flanks ...
					if ($tabOutput)
					{
						$MyLCR2WFlanks = $MyLeftFlank2Compare . "(" . $MyLCR2CompareLong . ")" . $MyRightFlank2Compare;
					}
					else
					{
						$MyLCR2WFlanks = $MyLeftFlank2Compare . "<b>" . $MyLCR2CompareLong . "</b>" . $MyRightFlank2Compare;
					}

					#Subjects.Length ...
					$MyLCR2SubjectLength = TrimLeadingAndTrailingLength($MyLCR2CompareLong);
					print ("lcrlength2compare: $LCRLength2Compare\n");
					my $LenDifLenRatio = abs( $LCR2QueryLength - $MyLCR2SubjectLength ) / $LCRLength2Compare;

					#Additional calculations for Brief output ...
					my $MisMatchesLengthRatio = $FindMismatches / $LCRLength2Compare;
					my $FlankMisMatches = $FlankLeftMismatches + $FlankRightMismatches;
					my $FlankGaps = $FlankLeftGaps + $FlankRightGaps;
					my $FlankGapMatches = $FlankLeftGapMatches +$FlankRightGapMatches;

					########################################################################################################################
					#A-Comparisons of entire sequence minus the current LCR AlignedStart and End ...
					my ($A_Matches, $A_MisMatches, $A_Gaps, $A_GapMatches);
					$A_Matches = Find_Matches_Exclude($Protein2Compare, $MyProtein2Compare, $AlignedStart[$i], $AlignedEnd[$i]);
					$A_MisMatches = Find_MisMatches_Exclude($Protein2Compare, $MyProtein2Compare, $AlignedStart[$i], $AlignedEnd[$i]);
					$A_Gaps = Find_Gaps_Exclude($Protein2Compare, $MyProtein2Compare, $AlignedStart[$i], $AlignedEnd[$i]);
					$A_GapMatches = Find_Gap_Matches_Exclude ($Protein2Compare, $MyProtein2Compare, $AlignedStart[$i], $AlignedEnd[$i]);
					my ($B_Matches, $B_MisMatches, $B_Gaps, $B_GapMatches);
					$B_Matches = Find_Matches_Exclude($Protein2Compare, $MyProtein2Compare, $LCRStart[$i], $LCREnd[$i]);
					$B_MisMatches = Find_MisMatches_Exclude($Protein2Compare, $MyProtein2Compare, $LCRStart[$i], $LCREnd[$i]);
					$B_Gaps = Find_Gaps_Exclude($Protein2Compare, $MyProtein2Compare, $LCRStart[$i], $LCREnd[$i]);
					$B_GapMatches= Find_Gap_Matches_Exclude	($Protein2Compare, $MyProtein2Compare, $LCRStart[$i], $LCREnd[$i]);
=pod					
					#Now each LCR should be excluded ...
					#This is NOT the most desired way to do this ... But, perhaps most easily understood/done ...
					#
					my ($C_Matches, $C_MisMatches, $C_Gaps, $C_GapMatches);
					my ($D_Matches, $D_MisMatches, $D_Gaps, $D_GapMatches);
					my $ExcludeFlag;

					$C_Matches = 0;
					for (my $k=0; $k < length($Protein2Compare); $k++)
					{
						if ( (substr($Protein2Compare, $k, 1) ne "-") && (substr($Protein2Compare, $k, 1) eq substr($MyProtein2Compare, $k, 1)) )
						{
							$ExcludeFlag = 0;

							for (my $L=0; $L < scalar(@Species); $L++)
							{
								if ( ($Species[$L] eq $Species[$i]) && ($Genes[$L] eq $Genes[$i]) )
								{
									if ( ($k >= $AlignedStart[$L]) && ($k <= $AlignedEnd[$L]) ) #Exclude range ...
									{
										$ExcludeFlag = 1;
										last;
									}
								}
							}

							if ($ExcludeFlag != 1) #Exclude flag value is 1 ...
							{
								$C_Matches++;
							}
						}
					}

					$C_MisMatches = 0;
					for (my $k=0; $k < length($Protein2Compare); $k++)
					{
						if ( (substr($Protein2Compare, $k, 1) ne "-") && (substr($MyProtein2Compare, $k, 1) ne "-") && (substr($Protein2Compare, $k, 1) ne substr($MyProtein2Compare, $k, 1)) )
						{
							$ExcludeFlag = 0;

							for (my $L=0; $L < scalar(@Species); $L++)
							{
								if ( ($Species[$L] eq $Species[$i]) && ($Genes[$L] eq $Genes[$i]) )
								{
									if ( ($k >= $AlignedStart[$L]) && ($k <= $AlignedEnd[$L]) ) #Exclude range ...
									{
										$ExcludeFlag = 1;
										last;
									}
								}
							}

							if ($ExcludeFlag != 1) #Exclude flag value is 1 ...
							{
								$C_MisMatches++;
							}
						}
					}

					$C_Gaps = 0;
					for (my $k=0; $k < length($Protein2Compare); $k++)
					{
						if ( (substr($Protein2Compare, $k, 1) eq "-") && (substr($MyProtein2Compare, $k, 1) ne "-") )
						{
							$ExcludeFlag = 0;

							for (my $L=0; $L < scalar(@Species); $L++)
							{
								if ( ($Species[$L] eq $Species[$i]) && ($Genes[$L] eq $Genes[$i]) )
								{
									if ( ($k >= $AlignedStart[$L]) && ($k <= $AlignedEnd[$L]) ) #Exclude range ...
									{
										$ExcludeFlag = 1;
										last;
									}
								}
							}

							if ($ExcludeFlag != 1) #Exclude flag value is 1 ...
							{
								$C_Gaps++;
							}
						}
						if ( (substr($MyProtein2Compare, $k, 1) eq "-") && (substr($Protein2Compare, $k, 1) ne "-") )
						{
							$ExcludeFlag = 0;

							for (my $L=0; $L < scalar(@Species); $L++)
							{
								if ( ($Species[$L] eq $Species[$i]) && ($Genes[$L] eq $Genes[$i]) )
								{
									if ( ($k >= $AlignedStart[$L]) && ($k <= $AlignedEnd[$L]) ) #Exclude range ...
									{
										$ExcludeFlag = 1;
										last;
									}
								}
							}

							if ($ExcludeFlag != 1) #Exclude flag value is 1 ...
							{
								$C_Gaps++;
							}
						}
					}
				$C_GapMatches= 0;
					for (my $k=0; $k < length($Protein2Compare); $k++)
					{
						if ( (substr($Protein2Compare, $k, 1) eq "-") && (substr($MyProtein2Compare, $k, 1) eq "-") )
						{
							$ExcludeFlag = 0;

							for (my $L=0; $L < scalar(@Species); $L++)
							{
								if ( ($Species[$L] eq $Species[$i]) && ($Genes[$L] eq $Genes[$i]) )
								{
									if ( ($k >= $AlignedStart[$L]) && ($k <= $AlignedEnd[$L]) ) #Exclude range ...
									{
										$ExcludeFlag = 1;
										last;
									}
								}
							}

							if ($ExcludeFlag != 1) #Exclude flag value is 1 ...
							{
								$C_GapMatches++;
							}
						}
					}
					$D_Matches = 0;
					for (my $k=0; $k < length($Protein2Compare); $k++)
					{
						if ( (substr($Protein2Compare, $k, 1) ne "-") && (substr($Protein2Compare, $k, 1) eq substr($MyProtein2Compare, $k, 1)) )
						{
							$ExcludeFlag = 0;

							for (my $L=0; $L < scalar(@Species); $L++)
							{
								if ( ($Species[$L] eq $Species[$i]) && ($Genes[$L] eq $Genes[$i]) )
								{
									if ( ($k >= $LCRStart[$L]) && ($k <= $LCREnd[$L]) ) #Exclude range ...
									{
										$ExcludeFlag = 1;
										last;
									}
								}
							}

							if ($ExcludeFlag != 1) #Exclude flag value is 1 ...
							{
								$D_Matches++;
							}
						}
					}

					$D_MisMatches = 0;
					for (my $k=0; $k < length($Protein2Compare); $k++)
					{
						if ( (substr($Protein2Compare, $k, 1) ne "-") && (substr($MyProtein2Compare, $k, 1) ne "-") && (substr($Protein2Compare, $k, 1) ne substr($MyProtein2Compare, $k, 1)) )
						{
							$ExcludeFlag = 0;

							for (my $L=0; $L < scalar(@Species); $L++)
							{
								if ( ($Species[$L] eq $Species[$i]) && ($Genes[$L] eq $Genes[$i]) )
								{
									if ( ($k >= $LCRStart[$L]) && ($k <= $LCREnd[$L]) ) #Exclude range ...
									{
										$ExcludeFlag = 1;
										last;
									}
								}
							}

							if ($ExcludeFlag != 1) #Exclude flag value is 1 ...
							{
								$D_MisMatches++;
							}
						}
					}

					$D_Gaps = 0;
					for (my $k=0; $k < length($Protein2Compare); $k++)
					{
						if ( (substr($Protein2Compare, $k, 1) eq "-") && (substr($MyProtein2Compare, $k, 1) ne "-") )
						{
							$ExcludeFlag = 0;

							for (my $L=0; $L < scalar(@Species); $L++)
							{
								if ( ($Species[$L] eq $Species[$i]) && ($Genes[$L] eq $Genes[$i]) )
								{
									if ( ($k >= $LCRStart[$L]) && ($k <= $LCREnd[$L]) ) #Exclude range ...
									{
										$ExcludeFlag = 1;
										last;
									}
								}
							}

							if ($ExcludeFlag != 1) #Exclude flag value is 1 ...
							{
								$D_Gaps++;
							}
						}
						if ( (substr($MyProtein2Compare, $k, 1) eq "-") && (substr($Protein2Compare, $k, 1) ne "-") )
						{
							$ExcludeFlag = 0;

							for (my $L=0; $L < scalar(@Species); $L++)
							{
								if ( ($Species[$L] eq $Species[$i]) && ($Genes[$L] eq $Genes[$i]) )
								{
									if ( ($k >= $LCRStart[$L]) && ($k <= $LCREnd[$L]) ) #Exclude range ...
									{
										$ExcludeFlag = 1;
										last;
									}
								}
							}

							if ($ExcludeFlag != 1) #Exclude flag value is 1 ...
							{
								$D_Gaps++;
							}
						}
					}
				$D_GapMatches= 0;
					for (my $k=0; $k < length($Protein2Compare); $k++)
					{
						if ( (substr($Protein2Compare, $k, 1) eq "-") && (substr($MyProtein2Compare, $k, 1) eq "-") )
						{
							$ExcludeFlag = 0;

							for (my $L=0; $L < scalar(@Species); $L++)
							{
								if ( ($Species[$L] eq $Species[$i]) && ($Genes[$L] eq $Genes[$i]) )
								{
									if ( ($k >= $AlignedStart[$L]) && ($k <= $AlignedEnd[$L]) ) #Exclude range ...
									{
										$ExcludeFlag = 1;
										last;
									}
								}
							}

							if ($ExcludeFlag != 1) #Exclude flag value is 1 ...
							{
								$D_GapMatches++;
							}
						}
					}
=cut
					##### Debug Prints #####
					#print "$Species[$i] and $Genes[$i]\n";
					#print "$Protein2Compare\n";
					#print "$MyProtein2Compare\n";

					#for (my $k=0; $k < scalar(@LCRStart); $k++)
					#{
					#	print "$Species[$k], $Genes[$k], $AlignedStart[$k]-$AlignedEnd[$k] and $LCRStart[$k]--$LCREnd[$k]\n";
					#}

					#print "$A_Matches, $A_MisMatches, $A_Gaps\n";
					#print "$B_Matches, $B_MisMatches, $B_Gaps\n";
					#print "$C_Matches, $C_MisMatches, $C_Gaps\n";
					#print "$D_Matches, $D_MisMatches, $D_Gaps\n";

					#exit;
					########################################################################################################################

					#Handlig printing after introduction of $SkipFlag to avoid essentially duplicate rows ...
					my @FieldsSplit = split(" ", $SkipFlag[$i]);
					#These are all strings ...
					#$FieldsSplit[0] is $StartLocation ...
					#$FieldsSplit[1] is $EndLocation ...
					#$FieldsSplit[2] is $AlignedStart ...
					#$FieldsSplit[3] is $AlignedEnd ...

#					if (not $Brief)
					if (1)
					{
						if ($tabOutput)
						{
							#print FILE00 "<tr><td>$InputFilePrefix</td><td>$Specie2Compare</td><td>$Gene2Compare</td><td>$MySpecie2Compare</td><td>$FieldsSplit[0]</td><td>$FieldsSplit[1]</td><td><center>$FieldsSplit[2]</center></td><td><center>$FieldsSplit[3]</center></td><td>",$LCR2CompareStart+1,"</td><td>",$LCR2CompareEnd+1,"</td><td>$LCR2WFlanks</td><td>$MyLCR2WFlanks</td><td>$ProtenLength2Compare</td><td>$LCRLength2Compare</td><td>$FindMatches</td><td>$FindMismatches</td><td>$FindGaps</td><td>$FindGapMatches</td><td>$FlankLeftMatches</td><td>$FlankLeftMismatches</td><td>$FlankLeftGaps</td><td>$FlankLeftGapMatches</td><td>$FlankRightMatches</td><td>$FlankRightMismatches</td><td>$FlankRightGaps</td><td>$FlankRightGapMatches</td><td>$AlignmentCorrupted</td>";
							print FILE00 "<tr><td>$InputFilePrefix</td><td>$Specie2Compare</td><td>$Gene2Compare</td><td>$MySpecie2Compare</td><td>$FieldsSplit[0]</td><td>$FieldsSplit[1]</td><td>",$LCR2CompareStart+1,"</td><td>",$LCR2CompareEnd+1,"</td><td>$LCR2WFlanks</td><td>$MyLCR2WFlanks</td><td>$ProtenLength2Compare</td><td>$LCRLength2Compare</td><td>$FindMatches</td><td>$FindMismatches</td><td>$FindGaps</td><td>$FindGapMatches</td><td>$FlankLeftMatches</td><td>$FlankLeftMismatches</td><td>$FlankLeftGaps</td><td>$FlankLeftGapMatches</td><td>$FlankRightMatches</td><td>$FlankRightMismatches</td><td>$FlankRightGaps</td><td>$FlankRightGapMatches</td><td>$AlignmentCorrupted</td>";
							print FILE00 "<td>$A_Matches</td><td>$A_MisMatches</td><td>$A_Gaps</td><td>$A_GapMatches</td>\n";#<td>$B_Matches</td><td>$B_MisMatches</td><td>$B_Gaps</td><td>$B_GapMatches</td><td>$C_Matches</td><td>$C_MisMatches</td><td>$C_Gaps</td><td>$C_GapMatches</td><td>$D_Matches</td><td>$D_MisMatches</td><td>$D_Gaps</td><td>$D_GapMatches</td></tr>\n";


							#print FILE01 "$InputFilePrefix\t$Specie2Compare\t$Gene2Compare\t$MySpecie2Compare\t$FieldsSplit[0]\t$FieldsSplit[1]\t$FieldsSplit[2]\t$FieldsSplit[3]\t",$LCR2CompareStart+1,"\t",$LCR2CompareEnd+1,"\t$LCR2WFlanks\t$MyLCR2WFlanks\t$ProtenLength2Compare\t$LCRLength2Compare\t$FindMatches\t$FindMismatches\t$FindGaps\t$FindGapMatches\t$FlankLeftMatches\t$FlankLeftMismatches\t$FlankLeftGaps\t$FlankLeftGapMatches\t$FlankRightMatches\t$FlankRightMismatches\t$FlankRightGaps\t$FlankRightGapMatches\t$AlignmentCorrupted\t";
							print FILE01 "$InputFilePrefix\t$Specie2Compare\t$Gene2Compare\t$MySpecie2Compare\t$FieldsSplit[0]\t$FieldsSplit[1]\t",$LCR2CompareStart+1,"\t",$LCR2CompareEnd+1,"\t$LCR2WFlanks\t$MyLCR2WFlanks\t$ProtenLength2Compare\t$LCRLength2Compare\t$FindMatches\t$FindMismatches\t$FindGaps\t$FindGapMatches\t$FlankLeftMatches\t$FlankLeftMismatches\t$FlankLeftGaps\t$FlankLeftGapMatches\t$FlankRightMatches\t$FlankRightMismatches\t$FlankRightGaps\t$FlankRightGapMatches\t$AlignmentCorrupted\t";
							print FILE01 "$A_Matches\t$A_MisMatches\t$A_Gaps\t$A_GapMatches\n";#\t$B_Matches\t$B_MisMatches\t$B_Gaps\t$B_GapMatches\t$C_Matches\t$C_MisMatches\t$C_Gaps\t$C_GapMatches\t$D_Matches\t$D_MisMatches\t$D_Gaps\t$D_GapMatches\n";



#						Comment out csv when tab option is used

#							print FILE02 #"$InputFilePrefix,$Specie2Compare,$Gene2Compare,$MySpecie2Compare,$FieldsSplit[0],$FieldsSplit[1],$FieldsSplit[2],$FieldsSplit[3],",$LCR2CompareStart+1,",",$LCR2CompareEnd+1,",$LCR2WFlanks,$MyLCR2WFlanks,$ProtenLength2Compare,$LCRLength2Compare,$FindMatches,$FindMismatches,$FindGaps,$FindGapMatches,$FlankLeftMatches,$FlankLeftMismatches,$FlankLeftGaps,$FlankLeftGapMatches,$FlankRightMatches,$FlankRightMismatches,$FlankRightGaps,$FlankRightGapMatches,$AlignmentCorrupted,";
#							print FILE02 "$A_Matches,$A_MisMatches,$A_Gaps,$A_GapMatches,$B_Matches,$B_MisMatches,$B_Gaps,$B_GapMatches,$C_Matches,$C_MisMatches,$C_Gaps,$C_GapMatches,$D_Matches,$D_MisMatches,$D_Gaps,$D_GapMatches\n";






						}
						else
						{


							#print FILE00 "<tr><td>$InputFilePrefix</td><td>$Specie2Compare</td><td>$Gene2Compare</td><td>$MySpecie2Compare</td><td>$FieldsSplit[0]</td><td>$FieldsSplit[1]</td><td><center>$FieldsSplit[2]</center></td><td><center>$FieldsSplit[3]</center></td><td>",$LCR2CompareStart+1,"</td><td>",$LCR2CompareEnd+1,"</td><td>$LCR2WFlanks</td><td>$MyLCR2WFlanks</td><td>$ProtenLength2Compare</td><td>$LCRLength2Compare</td><td>$FindMatches</td><td>$FindMismatches</td><td>$FindGaps</td><td>$FindGapMatches</td><td>$FlankLeftMatches</td><td>$FlankLeftMismatches</td><td>$FlankLeftGaps</td><td>$FlankLeftGapMatches</td><td>$FlankRightMatches</td><td>$FlankRightMismatches</td><td>$FlankRightGaps</td><td>$FlankRightGapMatches</td><td>$AlignmentCorrupted</td>";
							print FILE00 "<tr><td>$InputFilePrefix</td><td>$Specie2Compare</td><td>$Gene2Compare</td><td>$MySpecie2Compare</td><td>$FieldsSplit[0]</td><td>$FieldsSplit[1]</td><td>",$LCR2CompareStart+1,"</td><td>",$LCR2CompareEnd+1,"</td><td>$LCR2WFlanks</td><td>$MyLCR2WFlanks</td><td>$ProtenLength2Compare</td><td>$LCRLength2Compare</td><td>$FindMatches</td><td>$FindMismatches</td><td>$FindGaps</td><td>$FindGapMatches</td><td>$FlankLeftMatches</td><td>$FlankLeftMismatches</td><td>$FlankLeftGaps</td><td>$FlankLeftGapMatches</td><td>$FlankRightMatches</td><td>$FlankRightMismatches</td><td>$FlankRightGaps</td><td>$FlankRightGapMatches</td><td>$AlignmentCorrupted</td>";
							print FILE00 "<td>$A_Matches</td><td>$A_MisMatches</td><td>$A_Gaps</td><td>$A_GapMatches</td>\n";#<td>$B_Matches</td><td>$B_MisMatches</td><td>$B_Gaps</td><td>$B_GapMatches</td><td>$C_Matches</td><td>$C_MisMatches</td><td>$C_Gaps</td><td>$C_GapMatches</td><td>$D_Matches</td><td>$D_MisMatches</td><td>$D_Gaps</td><td>$D_GapMatches</td></tr>\n";



#						Comment out tab when option is not used, use csv by default

#							print FILE01 #"$InputFilePrefix\t$Specie2Compare\t$Gene2Compare\t$MySpecie2Compare\t$FieldsSplit[0]\t$FieldsSplit[1]\t$FieldsSplit[2]\t$FieldsSplit[3]\t",$LCR2CompareStart+1,"\t",$LCR2CompareEnd+1,"\t$LCR2WFlanks\t$MyLCR2WFlanks\t$ProtenLength2Compare\t$LCRLength2Compare\t$FindMatches\t$FindMismatches\t$FindGaps\t$FindGapMatches\t$FlankLeftMatches\t$FlankLeftMismatches\t$FlankLeftGaps\t$FlankLeftGapMatches\t$FlankRightMatches\t$FlankRightMismatches\t$FlankRightGaps\t$FlankRightGapMatches\t$AlignmentCorrupted\t";
#							print FILE01 "$A_Matches\t$A_MisMatches\t$A_Gaps\t$A_GapMatches\t$B_Matches\t$B_MisMatches\t$B_Gaps\t$B_GapMatches\t$C_Matches\t$C_MisMatches\t$C_Gaps\t$C_GapMatches\t$D_Matches\t$D_MisMatches\t$D_Gaps\t$D_GapMatches\n";

							$Gene2Compare =~ s/\r//;
							print FILE02 "$InputFilePrefix,$Specie2Compare,$Gene2Compare,$MySpecie2Compare,$FieldsSplit[0],$FieldsSplit[1],",$LCR2CompareStart+1,",",$LCR2CompareEnd+1,",$LCR2WFlanks,$MyLCR2WFlanks,$ProtenLength2Compare,$LCRLength2Compare,$FindMatches,$FindMismatches,$FindGaps,$FindGapMatches,$FlankLeftMatches,$FlankLeftMismatches,$FlankLeftGaps,$FlankLeftGapMatches,$FlankRightMatches,$FlankRightMismatches,$FlankRightGaps,$FlankRightGapMatches,$AlignmentCorrupted,";
							#print FILE02 "$InputFilePrefix,$Specie2Compare,$Gene2Compare,$MySpecie2Compare,$FieldsSplit[0],$FieldsSplit[1],",$LCR2CompareStart+1,",",$LCR2CompareEnd+1,",$LCR2WFlanks,$MyLCR2WFlanks,$ProtenLength2Compare,$LCRLength2Compare,$FindMatches,$FindMismatches,$FindGaps,$FindGapMatches,$FlankLeftMatches,$FlankLeftMismatches,$FlankLeftGaps,$FlankLeftGapMatches,$FlankRightMatches,$FlankRightMismatches,$FlankRightGaps,$FlankRightGapMatches,$AlignmentCorrupted,";
							print FILE02 "$A_Matches,$A_MisMatches,$A_Gaps,$A_GapMatches\n";#,$B_Matches,$B_MisMatches,$B_Gaps,$B_GapMatches,$C_Matches,$C_MisMatches,$C_Gaps,$C_GapMatches,$D_Matches,$D_MisMatches,$D_Gaps,$D_GapMatches\n";







						}
					}


				}
			}
		}
	}

	#if (not $tabOutput)
	#{
	print FILE00 "</font></table>";

#}
	close FILE00;
	#Next file in the folder ...
}


print "\nDone!\n";

exit;

################################################################################
# Subroutines
################################################################################

sub Read_Into_Arrays {

    # Collect array arguments passed by reference ...
    my($InputFile, $A, $B) = @_;

	open FILE, $InputFile or die "Can't read source file: $InputFile\n";
	my @species = <FILE>;

	#Parse the file data into 2 parallel arrays ...

	my $i = 0;

	while ($i < scalar (@species))
		{	if ($species[$i] =~ /^>/)

				{
					chomp($species[$i]);
					$species[$i] =~ s/\r//; #LinuxFix
										push @$A, $species[$i];

					my $protein = "";

					$i = $i+1;
					do
						{	my $line = chomp ($species[$i]);
							$species[$i] =~ s/\r//;

							$protein = $protein . $species[$i];
							$i = $i + 1;
						} while (($i < scalar @species) && !($species[$i] =~ /^>/));

					push @$B, $protein;
				}
		}
	close FILE;
}

sub Find_Matches {

	my($A, $B) = @_;
	my $Matches = 0;

	for (my $i=0; $i < length($A); $i++)
	{
		if ( (substr($A, $i, 1) ne "-") && (substr($A, $i, 1) eq substr($B, $i, 1)) )
		{
			$Matches++;
		}
	}
	return $Matches;
}


sub Find_Matches_Exclude {

	my($A, $B, $x, $y) = @_;
	my $Matches = 0;

	for (my $i=0; $i < length($A); $i++)
	{
		if ( (substr($A, $i, 1) ne "-") && (substr($A, $i, 1) eq substr($B, $i, 1)) )
		{
			if (($i < $x) || ($i > $y)) #Exclude range ...
			{
				$Matches++;
			}
		}
	}
	return $Matches;
}

sub Find_MisMatches {

	my ($A, $B)=@_;
	my $MisMatches = 0;
	for (my $i=0; $i < length($A); $i++)
	{
		if ( (substr($A, $i, 1) ne "-") && (substr($B, $i, 1) ne "-") && (substr($A, $i, 1) ne substr ($B, $i, 1)) )
			{
				$MisMatches++;
			}
	}
	return $MisMatches;
}


sub Find_MisMatches_Exclude {

	my($A, $B, $x, $y) = @_;
	my $MisMatches = 0;

	for (my $i=0; $i < length($A); $i++)
	{
		if ( (substr($A, $i, 1) ne "-") && (substr($B, $i, 1) ne "-") && (substr($A, $i, 1) ne substr ($B, $i, 1)) )
		{
			if (($i < $x) || ($i > $y)) #Exclude range ...
			{
				$MisMatches++;
			}
		}
	}
	return $MisMatches;
}

sub Find_Gaps {
	my ($A, $B)=@_;
	my $Gaps = 0;
	for (my $i=0; $i< length($A); $i++)
	{
		if ( (substr($A, $i, 1) eq "-") && (substr($B, $i, 1) ne "-") )
			{
				$Gaps++;
			}
		if ( (substr($B, $i, 1) eq "-") && (substr($A, $i, 1) ne "-") )
			{
				$Gaps++;
			}
	}
	return $Gaps;
}
sub Find_GapMatches {
	my ($A, $B)=@_;
	my $GapMatches = 0;
	for (my $i=0; $i< length($A); $i++)
	{
		if ( (substr($A, $i, 1) eq "-") && (substr($B, $i, 1) eq "-") )
			{
				$GapMatches++;
			}
	}
	return $GapMatches;
}
sub Find_Gaps_Exclude {
	my ($A, $B, $x, $y)=@_;
	my $Gaps = 0;
	for (my $i=0; $i< length($A); $i++)
	{
		if ( (substr($A, $i, 1) eq "-") && (substr($B, $i, 1) ne "-") )
			{
				if (($i < $x) || ($i > $y)) #Exclude range ...
				{
					$Gaps++;
				}
			}
		if ( (substr($B, $i, 1) eq "-") && (substr($A, $i, 1) ne "-") )
			{
				if (($i < $x) || ($i > $y)) #Exclude range ...
				{
					$Gaps++;
				}
			}
	}
	return $Gaps;
}
sub Find_Gap_Matches_Exclude {
	 my ($A, $B, $x, $y)=@_;
	 my $GapMatches = 0;
	for (my $i=0; $i< length($A); $i++)
	 {
		if ( (substr($A, $i, 1) eq "-") && (substr($B, $i, 1) eq "-") )
			 {
				if (($i < $x) || ($i > $y)) #Exclude range ...
				 {
					$GapMatches++;
				 }
			 }
	 }
	 return $GapMatches;
 }

sub TrimLeadingAndTrailingLength {
	my ($A)=@_;
	my $TrimmedStart = 0;
	my $TrimmedEnd = length($A) -1;

	while ( ($TrimmedStart < length($A)) && (substr($A, $TrimmedStart, 1) eq "-") )
	{
		$TrimmedStart++;
	}
	if ( $TrimmedStart == length($A) )
	{
		#All - ...
		return 0;
	}
	else
	{
		while ( ($TrimmedEnd >= 0) && (substr($A, $TrimmedEnd, 1) eq "-") )
		{
			$TrimmedEnd--;
		}

		return $TrimmedEnd - $TrimmedStart +1;
	}
}

# Perl trim function to remove whitespace from the start and end of the string ...
sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

#---------------------------#
# Quick Start Documentation #
#---------------------------#
sub QuickStart {
print "#---------------------------#\n";
print "# Quick Start Documentation #\n";
print "#---------------------------#\n";
print "\n";
printf "%14s %-60s\n", "Purpose:", "To compare segments identified in a .fasta file against others in the file";
printf "%14s %-60s\n", "Input Files:", "Pairs of files (a .fasta file and a file identifying segments)";
printf "%14s %-60s\n", "Output Files:", "An HTML file displaying the comparisons";
printf "%14s %-60s\n", "Notes:", "Two types of files identifying segments are supported";
printf "%14s %-60s\n", ":", "All files with the specified extension in the folder will be";
printf "%14s %-60s\n", ":", "processed. Please see sample expected file names and data format.";
printf "%14s %-60s\n", "Parameters:", "-b Brief output format. Optional.";
printf "%14s %-60s\n", ":", "-c filename.htm Produces Combined output file. Optional.";
printf "%14s %-60s\n", ":", "-d DesiredFlankingRegionLength.";
printf "%14s %-60s\n", ":", "-dir DataDirectoryFolder. Optional.";
printf "%14s %-60s\n", ":", "-e FileExtensionToUse Identifies extension (e.g., .txt) of segment files.";
printf "%14s %-60s\n", ":", "-t produces tab-delimited output. Optional.";
printf "%14s %-60s\n", ":", "-h displays help. Optional.";


print "\n";
print "Sample .fasta file:\n\n";
print "RodentOG490Mod.fasta\n";
print ">pyoelii_17x.fasta\n";
print "M-SIGVPVTKNNNNNNNNNNNFNKLKKYPYNNNNNNINKFDRRETNEY----PNNNNNMK\n";
print "... etc.\n\n";

print "Sample segment file:\n\n";
print "RodentOG490MOD.HPrs6-0-0.out\n";
print ">pyoelii_17x.fasta|gi|PY17X_0102400|ref||(10-20) complexity=-0.00 (6/0.00/0.00)\n";
print "... etc.\n\n";

print "Sample .txt (comma delimited) segment file:\n\n";
print "RodentOG490MOD.txt\n";
print "pyoelii_17x, PY17X_0102400,  10, 20\n";
print "... etc.\n\n";

print "\n";
print "#---------------------------#\n";
print "\n";

return 0;
}
