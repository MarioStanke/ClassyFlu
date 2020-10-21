#!/usr/bin/perl
# Sandra Van der Auwera
# Mai 08, 2012
# classification using HMM profile database
# usage: hmmClassify.pl input.fa outputname 


use strict;
use warnings;

use Getopt::Long;

my ($seqfile, $output, $help, $db, $lic, $profileDir);	# variables for argument
my $i;
my $c;
my @seqn=();
my @seq=();
my $line='';
my @arg;
my @eval;
my @score;
my $m;
my $l;
my $temp;
my @arg1=();
my $check=1;
my $le;

GetOptions('in=s'=>\$seqfile,
	   'out=s'=>\$output,
	   'help!'=>\$help,
	   'db=s'=>\$db,
	   'lic!'=>\$lic,
	   'profileDir=s'=>\$profileDir);


my $gpl = <<'ENDGPL';

hmmClassify  Copyright (C) 2012  author: Sandra Van der Auwera
    This program comes with ABSOLUTELY NO WARRANTY; for details type --lic.
    This is free software, and you are welcome to redistribute it
    under certain conditions; see GNU_GPL for details.

ENDGPL

print $gpl;

my $usage = <<'ENDUSAGE';

#########################################
#                                       #
# README FOR hmmClassify.pl             #
#                                       #
#########################################

DESCRIPTION: hmmClassify.pl is a tool to classify influenza A HA-gene sequences based on hidden markov profiles (HMM profile database)
		
REQUIRED SOFTWARE:
	- PERL
	- HMMER


SYNOPSIS:

hmmClassify.pl --in input.fa --out output

	input.fa is an unaligned or aligned sequence file in FASTA format
	output is the name one can choose for the output file (result)


OPTIONS:

    	--help          output this help message
    	--db		specifying the hmm database to use (by default the program uses the file hmmClassifydb_new)
	--lic 		outputs the WARRANTY
	--profileDir	the directory of the hmm database and the cutoff file (if it's not the same as for the program)
	                       

EXAMPLE:

	hmmClassify.pl --in sequences.fa --out results --cut 0 --profileDir usr/me/myfolder

                 
OUTPUT:	

	text file with two rows per input sequence: 
	- sequence name; length of the sequence, highest score, predicted model (tab-delimited) if a match was found
	- sequence name; "No Hit against one of the models." if no match was found
	- sequence name; "Sequence si to short for classification.", if the input sequence is to short. 

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

if ($help){
    	print $usage;
    	exit(1);
}
if (!defined($seqfile)){
    	print "Missing input sequence file.\n$usage";
    	exit(1);
}

if (!defined($output)){
    	print "No output file name specified.\n$usage";
    	exit(1);
}

if (!defined($db)){
	$db="hmmdb_full";
}

# the current directory $path and the directory of the program $profileDir

if (!defined($profileDir)){
	$temp=`which hmmClassify.pl`;
	if($temp eq ''){
		$profileDir=`pwd`;
	}
	else{
		chomp($temp);
		@arg=split(/\//, $temp);
		$l=@arg;
		for($i=0; $i<($l-1); $i++){
			push(@arg1, $arg[$i]);
		}
		$profileDir=join("\/", @arg1);
	}
}
my $path=`pwd`;

# sequences for classification (sequencenames are stored in @seqn and the sequences in @seq)
open (SEQ, "<$seqfile");
$c=0;
while(<SEQ>){
	chomp($_);
	$_=~tr/-//d;				# eliminate gap characters (in aligned FASTA files)
	$_ =~ s/\r|\n//g;	
	if($check==1){
		if($_=~m/^>/){
			$seqn[$c]=$_;
			$seq[$c]=$line;
			$c++;
			$line='';
		}
		elsif($_=~m/\W/){
			$check=0;
		}
		else{
			if($_=~m/[EFIKLOQ]/i){
				$check=0;
			}
			else{
				$line=$line.$_;
			}			
		}
	}
}
$seq[$c]=$line;
close SEQ;
$temp=`ls`;
$temp=~s/.~/X/;
$temp=~s/.~/X/;

if($temp =~m/hmmClassify.pl/){
	if($path eq $profileDir){
		if($temp !~m/$db/){	
			print ("The hmm database file can not be found in the current working directory $path \n");
    			exit(1);
		}
	}
}
else{
	chdir($profileDir) or die "Can't chdir to $profileDir";
	$temp=`ls`;
	$temp=~s/.~/X/;
	if($temp !~m/$db/){
		print ("The hmm database file can not be found in the current working directory $profileDir \n");
   		exit(1);
	}
}


	

if($check==1 && $c>0){

	# prepare the hmm-database with hmmpress
	system("hmmpress $db");
	open (RES, ">$output");
	print RES ("sequence name\t length of sequence \t score\t H-type/subfamily\n");

	# a loop through the test sequences
	for($i=0; $i<$c; $i++){
	
		$le=length($seq[$i+1]);
		if($le>5){

			open (TESTSEQ, ">testseq.fa");
			print TESTSEQ ("$seqn[$i]\n");
			print TESTSEQ ("$seq[$i+1]\n");
			close TESTSEQ;
		
	
			# execute hmmscan
			system("hmmscan --noali -o eval $db testseq.fa");

			# read hmmscan output
			open (EVAL, "<eval");
			@eval=<EVAL>;
			close EVAL;
			system("rm eval");

			# if one got hits, find the best score (and compare it to the cutoffs!)
			$score[0]=0;
			$score[1]='';
			$m=16;
			@arg=split(/ +/, $eval[$m]);
			if($arg[2] =~m/-/){
				$m=$m+1
			}
			@arg=split(/ +/, $eval[$m]);
			while(defined $arg[2]){	
				if($arg[2]>$score[0]){
					$score[0]=$arg[2];
					$score[1]=$arg[9];	
				}
				$m++;
				@arg=split(/ +/, $eval[$m]);
			}
	
			# write output: match or no match
			chomp($seqn[$i]);
		
			if($score[0]!=0){
				$seqn[$i]=~ s/>//g;
				print RES ("$seqn[$i]\t $le \t $score[0]\t H-type: $score[1]\n");
			}
			elsif($score[0]==0){
				$seqn[$i]=~ s/>//g;
				print RES ("$seqn[$i]\t $le\t 0\t No hit against one of the models.\n");
			}	
		}
		else{
			$seqn[$i]=~ s/>//g;
			print RES ("$seqn[$i]\t $le \t - \t Sequence is too short for classification.\n");
		}
	}

	system("rm testseq.fa $db.h3f $db.h3i $db.h3m $db.h3p");
	close RES;
	
	if($path !~m/$profileDir/){
		system("mv -f $output $path");
	}
}
else{
	print ("The sequence file is in the wrong format.");
}
