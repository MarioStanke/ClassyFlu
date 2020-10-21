#!/usr/bin/perl
# Sandra Van der Auwera
# March 07, 2012
# training-algorithm for the hmm-files to get a trained HMM database
# all cluster files cl.fa



use strict;
use warnings;

use Getopt::Long;

my ($help, $lic);				# variables for the arguments
my $args="cat ";
my $i;
my @seq;
my @seqn;
my $line='';
my $c;
my $j;
my $f;
my @seqtemp;
my $ctemp;
my @var;
my @eval;
my @score;


GetOptions('help!'=>\$help,
	   'lic!'=>\$lic);


my $gpl = <<'ENDGPL';

hmmClassify  Copyright (C) 2012  author: Sandra Van der Auwera
    This program comes with ABSOLUTELY NO WARRANTY; for details type --lic.
    This is free software, and you are welcome to redistribute it
    under certain conditions; see GNU_GPL for details.

ENDGPL

print $gpl;


my $usage= <<'ENDUSAGE';

##############################
#                            #
# README FOR hmmTraining.pl  #
#                            #
##############################

DESCRIPTION: hmmTraining.pl is a tool to create a trained HMM database

REQUIRED SOFTWARE:
	- PERL
	- MUSCLE
	- HMMER

REQUIRED FILES:
	- the cluster files (unaligned FASTA format with extension .fa or .fasta) must be in the folder (names: cl1.fa, cl2.fa, ...) 

SYNOPSIS:
	hmmTraineval.pl cl1.fa cl2.fa ... 


OPTIONS:
	--help outputs this help message
	--lic outputs the WARRANTY

OUTPUT:
	- trained HMM database (trained_hmmdb)

ENDUSAGE


my $warranty = <<'ENDWARRANTY';

  15. Disclaimer of Warranty.

  THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY
APPLICABLE LAW.  EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT
HOLDERS AND/OR OTHER PARTIES PROVIDE THE PROGRAM "AS IS" WITHOUT WARRANTY
OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE.  THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM
IS WITH YOU.  SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF
ALL NECESSARY SERVICING, REPAIR OR CORRECTION.

  16. Limitation of Liability.

  IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MODIFIES AND/OR CONVEYS
THE PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES, INCLUDING ANY
GENERAL, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE
USE OR INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED TO LOSS OF
DATA OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD
PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER PROGRAMS),
EVEN IF SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF
SUCH DAMAGES.

  17. Interpretation of Sections 15 and 16.

  If the disclaimer of warranty and limitation of liability provided
above cannot be given local legal effect according to their terms,
reviewing courts shall apply local law that most closely approximates
an absolute waiver of all civil liability in connection with the
Program, unless a warranty or assumption of liability accompanies a
copy of the Program in return for a fee.

ENDWARRANTY


if ($lic){
	print $warranty;
	exit(1);
}

if($help){
	print $usage;
	exit(1);
}

if (!defined(@ARGV)){
    	print "Missing input cluster files.\n$usage";
    	exit(1);
}

my $nbcl=@ARGV;			# number of clusters in array @ARGV

# building multiple sequence alignments of the clusters
for($i=0; $i<$nbcl; $i++){
	
	system("muscle -in $ARGV[$i] -out $ARGV[$i].msa.fa -maxiters 2");
	&fasta2stockholm("$ARGV[$i].msa.fa", "$ARGV[$i].sto");
	system("hmmbuild --dna hmm.$ARGV[$i] $ARGV[$i].sto");
}


# a loop through all clusters
for(my $b=0; $b<$nbcl; $b++){
	$c=0;
	# sequences of the cluster are stored in @seq and their names in @seqn
	open (CLUSTER, "<$ARGV[$b]");
	while(<CLUSTER>){
		chomp($_);
		if($_=~m/^>/){
			$seqn[$c]=$_;
			$seq[$c]=$line;
			$c++;
			$line='';
		}
		else{
			$line=$line.$_;
		}
	}
	$seq[$c]=$line;
	close CLUSTER;
	$line="cat ";
	for($j=0; $j<$nbcl; $j++){
		$line=$line."hmm.$ARGV[$j] ";
	}
	$line=$line."> hmmdb";

	# a loop through the sequences
	for($i=0; $i<$c; $i++){
		# loop to adjust the hmm-parameters
		open (TESTSEQ, ">testseq.fa");
		print TESTSEQ ("$seqn[$i]\n");
		print TESTSEQ ("$seq[$i+1]\n");
		close TESTSEQ;

		@var=&hmmer($nbcl, $line, $ARGV[$b]);
	
		# correcting the parameters if the classification was wrong (for the real and the wrong hmm)
		if($var[0]!=0){
			$var[0]=$var[0];
			&hmmcorrect("hmm.$ARGV[$b]", $var[2], $var[3], $var[0], $var[1]);
			$var[0]=(-1)*$var[0];
			&hmmcorrect("hmm.$var[4]", $var[6], $var[7], $var[0], $var[5]);	
		}
	}
}
system("rm testseq.fa");
for($i=0; $i<@ARGV; $i++){
	$args= $args."hmm.$ARGV[$i] ";
	system("rm $ARGV[$i].msa.fa $ARGV[$i].sto");
} 
system("$args > trained_hmmdb");
for($i=0; $i<@ARGV; $i++){
	system("rm hmm.$ARGV[$i]");
} 		

##################################################

# convert the MSA from fasta to stockholm format

sub fasta2stockholm {
# loop through FASTA files

	# read FASTA file
	my %seq;
	my @name;
	my $name;
	open FASTA, "<$_[0]";
	while (<FASTA>) {
		if (/^\s*>\s*(\S+)/) {
			$name = $1;
			die "Duplicate name: $name" if defined $seq{$name};
			push(@name, $name);
		} 
		else {
			if (/\S/ && !defined $name) {
				warn "Ignoring: $_";
			} 
			else {
				s/\s//g;
				$seq{$name} .= $_;
			}
		}
	}
	close FASTA;

	# check all seqs are same length
	my $length;
	my $lname;
	foreach my $name (@name) {
		my $l = length $seq{$name};
		if (defined $length) {
			die "Sequences not all same length ($lname is $length, $name is $l)" unless $length == $l;
		} 
		else {
			$length = length $seq{$name};
			$lname = $name;
		}
	}
	open (OUT, ">$_[1]");

	# print Stockholm output
	print OUT ("# STOCKHOLM 1.0\n");
	foreach my $name (@name) {
		print OUT ($name, " ", $seq{$name}, "\n");
	}
	print OUT ("//\n");
	close OUT;
}

##################################################

sub hmmer{
		
	my @ret;
	my $k=15;
	my @temp1;
	my $temp2;
	my $cluster=$_[2];
	# execute the HMMER commands
	my @args1=$_[1];
	system(@args1);
	@args1=("hmmpress", "hmmdb");
	system(@args1);
	@args1=("hmmscan", "-o", "out", "hmmdb", "testseq.fa");
	system(@args1);
	open (SCAN, "<out");
	my @scan=<SCAN>;
	close SCAN;
	@args1=("rm", "out", "hmmdb", "hmmdb.h3f", "hmmdb.h3i", "hmmdb.h3m", "hmmdb.h3p");
	system(@args1);
	$ret[0]=0;
	$ret[1]=0;

	# finding the maximum score
	@args1=split(/ +/, $scan[$k]);
	while(defined $args1[2]){	
		if($args1[2]>$ret[0]){
			$ret[0]=$args1[2];
			$temp1[0]=$args1[9];				# hmm wrong
		}
		if($args1[9]=~m/$cluster/){
			$ret[1]=$args1[2];
		}
		$k++;
		@args1=split(/ +/, $scan[$k]);
	}

	$ret[0]=$ret[0]-$ret[1];
	$k=$k+1;
	$temp2=$k;
	if($ret[0]!=0){							# if the classification was wrong
		until($scan[$k]=~m/^(>> $cluster)/){
			$k++;
		}
		@args1=split(/ +/, $scan[$k+3]);
		$ret[1]=$args1[7];					# hmm start
		$ret[0]=$ret[0]/(($args1[8]-$args1[7]+1)*2);		# eta
		$k=$k+7;
		@args1=split(/ +/, $scan[$k]);
		$ret[2]='';
		$ret[3]='';

		# saving the MSA of the hmm and the sequence
		while(exists $args1[3]){
			$ret[2]=$ret[2].$args1[3];			# consensus hmm
			@args1=split(/ +/, $scan[$k+2]);
			$ret[3]=$ret[3].$args1[3];			# sequence
			$k=$k+5;
			@args1=split(/ +/, $scan[$k]);
		}
		$ret[2]=~tr/AaCcGgTt\./223344556/;
		$ret[2]=~tr/[A-Z]/7/;
		$ret[3]=~tr/AaCcGgTt-/223344556/;
		$ret[3]=~tr/[A-Z]/7/;
		$ret[4]=$temp1[0]; 

		until($scan[$temp2]=~m/^(>> $temp1[0])/){
			$temp2++;
		}
		@args1=split(/ +/, $scan[$temp2+3]);
		$ret[5]=$args1[7];					# hmm wrong start
		$temp2=$temp2+7;
		@args1=split(/ +/, $scan[$temp2]);
		$ret[6]='';
		$ret[7]='';

		# saving the MSA of the hmm and the sequence
		while(exists $args1[3]){
			$ret[6]=$ret[6].$args1[3];			# consensus hmm wrong
			@args1=split(/ +/, $scan[$temp2+2]);
			$ret[7]=$ret[7].$args1[3];			# sequence
			$temp2=$temp2+5;
			@args1=split(/ +/, $scan[$temp2]);
		}
		$ret[6]=~tr/AaCcGgTt\./223344556/;
		$ret[6]=~tr/[A-Z]/7/;
		$ret[7]=~tr/AaCcGgTt-/223344556/;
		$ret[7]=~tr/[A-Z]/7/;
	}
	# @ret=(eta, hmm start, hmm wrong, consensus rigt, sequence right, hmm wrong start, consensus wrong, sequence wrong)
	return(@ret);
}

##################################################

sub hmmcorrect{
	open (HMM, "<$_[0]");
	my @hmm=<HMM>;
	close HMM;
	my @args2=("rm", $_[0]);
	system(@args2);
	my @tempseq1=split(//, $_[1]);
	my @tempseq2=split(//, $_[2]);	
	
	# correcting each aligned column (tempseq1=hmm, tempseq2=targetseq)		
	for(my $f=0; $f<length($_[1]); $f++){
		@args2=split(/ +/, $hmm[16+(3*$_[4])]);		
		if($tempseq1[$f]!=6 && $tempseq2[$f]!=6){
			if($tempseq2[$f] !=7){
				$args2[$tempseq2[$f]]=$args2[$tempseq2[$f]]-$_[3];
				$args2[$tempseq2[$f]]=substr($args2[$tempseq2[$f]], 0, 7);
				my $l="      $_[4]";
				$l=substr($l, -7);
				my $templine=$l."   ".$args2[2]."  ".$args2[3]."  ".$args2[4]."  ".$args2[5]."      ".$args2[6]." - -\n";			
				$hmm[16+3*$_[4]]=$templine;
			}	
			$_[4]++;
		}
		if($tempseq2[$f]==6){
			$_[4]++;
		}					
	}
	# saving the new hmm-file
	open (HMM, ">$_[0]");
	print HMM ("@hmm");
	close HMM;
}
