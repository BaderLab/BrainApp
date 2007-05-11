#!/usr/bin/perl
# **********************************************************************************
# File: ConvertBrainPeptideFile.pl 
# Converts Brain peptide profiles from version 1.0 to version 1.1 
# Created: 2007 April 24 (Moyez Dharsee) 
# Updated: 2007 April 25 (Moyez Dharsee)
#               
# **********************************************************************************

# Globals
my $ScriptName = "ConvertBrainPeptideFile.pl";
my $Debug = 0; 		# set to 1 for verbose debugging output 
my $Comment = "";	# this string will be inserted in the Comment field
my $FilesConverted = 0;	# counter - # of files actually converted

##
# MAIN
##

##
# Read command-line args
# 1st argument is the file or directory name to process
# 2nd argument is the output directory

#if ($ARGC < 2) { 
#	print "\nERROR - Invalid arguments: $ARGV[0] $ARGV[1] $ARGC\n";
#	PrintUsage(); exit(0); 
#}
my $pathin = $ARGV[0]; 
my $pathout = $ARGV[1]; 
if (!defined($pathin) || (!defined($pathout)) || (length($pathin)<1) || (length($pathout)<1)) {
	print "You must provide an input file or directory path, and an output directory path.\n"; 
	exit(1);
}

# process the specified file or directory
if (-d $pathin) {
	ConvertDirectory($pathin, $pathout);
}
elsif (-f $pathin) {
	ConvertFile($pathin, $pathout);
}
else {
	print "Unable to process the path provided: $pathin\n";
	exit(1);
}
print "$FilesConverted files converted.\n";

##
# End MAIN
##

##
# Prints usage info to STDOUT 
##
sub PrintUsage() {
	print "\n";
	print "This script converts Brain peptide files from version 1.0 to version 1.1.\n";
	print "Original files are not modified. New files are written to the specified directory\n";
	print "\n";
	print "Usage: $ScriptName input_path output_path\n";
	print "\tinput_path\tPath to peptide file or directory containing peptide file(s)\n";
	print "\toutput_path\tPath to output directory (will be created if it doesn't exist)\n";
	print "\n";
}

##
# Converts peptide files in the given directory path
##
sub ConvertDirectory() {
	
	my $inPath = $_[0];
	if (!defined($inPath) || (length($inPath) == 0)) {
		die "Inavalid input directory argument";
	} 
	my $outPath = $_[1];
	if (!defined($outPath) || (length($outPath) == 0)) {
		die "Inavalid output directory argument";
	} 
	
	# iterate over files in input path
	opendir(DIR, $inPath) || die "Cannot open directory $inPath: $!";
	foreach $name (sort readdir(DIR)) {
		if (($name eq ".") || ($name eq "..")) { next; }
		ConvertFile($name, $outPath, $inPath);
	}
	closedir(DIR);

} # end sub


##
# Converts given input peptide file 
##
sub ConvertFile() {

	# input file
	my $filein = $_[0];
	if (!defined($filein) || (length($filein) == 0)) {
		die "Inavalid input file argument";
	} 
	# output directory
	my $dirout = $_[1];
	if (!defined($dirout) || (length($dirout) == 0)) {
		die "Inavalid output directory argument";
	} 
	# input directory - optional; set to empty string if not provided
	my $dirin = $_[2];
	my $haveDirin = 1;
	if (!defined($dirin) || (length($dirin) == 0)) {
		$dirin = "";
		$haveDirin = 0;
	} 

	print "Converting $filein\n";

	if ($Debug) {
		print "[debug] Input file: $filein\n[debug] Input Dir: $dirin\n[debug] Output Dir: $dirout\n";
	}
	
	##
	# glob file data
	my $filepath = $filein;
	if ($haveDirin) {
		$filepath = $dirin . "/" . $filepath;	# include input directory if provided
	}
	open(FILEIN, "<", $filepath) || die "Unable to open peptide file $filepath\n";
	my @datain = <FILEIN>;
	close(FILEIN);

	##
	# Convert the file
	my @dataout = undef;
	my $lineCount = 0;
	my $headerChanged = 0;
	my $isValid = 0;	# set to 1 if it's a valid peptide file
	for (@datain) {
		@lineSplit = split("\t", $_);
		$cols = (scalar @lineSplit);

		# ignore comments and blank lines
		if ((index($_,"#") == 0) || (length($_) < 1)) {
			next;
		}
		# if 2 columns, insert new fields at end of header 
		if ($cols == 2) {
			# a "valid" peptide file will contain the field "Gene Name" in first of 2 columns
			# this is a very loose validity check, but will do for now
			if ($isValid == 0) {
				if ($lineSplit[0] eq "Gene Name") {
					$isValid = 1;
				}
			}
			push (@dataout, $_);
			if ($lineSplit[0] eq "Domain Sequence") {
				push (@dataout, "Domain Range\t\n");
				push (@dataout, "Comment\t$Comment\n");
			}
		}
		# if 3 columns, replace peptide section column names 
		elsif ($cols == 3) {
			if ($headerChanged == 0) {
				$header = $_;
				chomp $header;
				$header = $header . "\tQuantData\tExternalIdentifier\n";
				push (@dataout, $header);
				$headerChanged = 1;
			}
			else {
				push (@dataout, $_);
			}
		}
	}

	# do nothing if file is not valid
	if ($isValid == 0) {
		print "NOTE: Skipping file $filein - not a valid peptide file.\n";
		return;
	}

	##
	# Write new peptide file

	# create output directory if necessary
	my $createdDir = 0;
	if (-d $dirout) {
		# nothing to do
		if ($Debug) {
			print "[debug] Directory '$dirout' exists.\n"; #debug
		}
	}
	else {
		mkdir $dirout || die "Unable to create output directory $dirout: $!";
		$createdDir = 1;
	}

	# open output file
	$fileout = $dirout . "/" . $filein;
	unless (open(FILEOUT, ">", $fileout)) {
		print STDERR "Can't open output file $fileout : $!\n";
		if ($createdDir) {
			rmdir $dirout || die "Unable to remove output directory $dirout: $!";
		}
		return;
	}

	# write it out
	for (@dataout) {
		print FILEOUT $_;
	}
	close(FILEOUT);

	# finish
	$FilesConverted++;
	print "=> $fileout\n";

} # end sub
