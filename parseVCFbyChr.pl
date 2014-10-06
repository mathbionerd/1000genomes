#!/usr/bin/perl -w

use strict;

my $line;
my $chr;
my $chrStart = $ARGV[0];
my $chrStop = $ARGV[1];
my $NUM_NONINDIV_COL = 9;
my $ID_COL = 0;
my $PAT_COL = 1;
my $MAT_COL = 2;
my $POP_COL = 3;

my %nt = ();
$nt{"A"} = 1;
$nt{"C"} = 1;
$nt{"G"} = 1;
$nt{"T"} = 1;

my @temp = ();
my @tempRefs = ();
my $nRefs;
my @tempAlts = ();
my $nAlts;
my $i;
my $numCols;
my @indivIDs = ();
my $numIndivs = 0;

my $dChr;
my $dPos;
my $dID;
my $dRef;
my $dAlt;
my $dQual;
my $dFilter;
my $dInfo;
my $dFormat;
my @dAlleles;
my $flagBad;
my @alleleInfo = ();

my @indivInfo = ();
my %indivInfoH = ();

open(IN, "phase1_samples_integrated_20101123.ped") or die("ERROR: Could not open phase1_samples_integrated_20101123.ped");
$line = <IN>; # Header
while(<IN>) {
	$line = $_;
	
	if(defined $line) {
		($indivInfo[$numIndivs][$ID_COL]) = ($line =~ m/\S+\s+(\S+)/);
		($indivInfo[$numIndivs][$PAT_COL]) = ($line =~ m/\S+\s+\S+\s+(\S+)/);
		($indivInfo[$numIndivs][$MAT_COL]) = ($line =~ m/\S+\s+\S+\s+\S+\s+(\S+)/);
		($indivInfo[$numIndivs][$POP_COL]) = ($line =~ m/\S+\s+\S+\s+\S+\s+\S+\s+\d+\s+\d+\s+(\S+)/);

		$indivInfo[$numIndivs][$POP_COL] = "1000G\_" . $indivInfo[$numIndivs][$POP_COL];
		$indivInfoH{$indivInfo[$numIndivs][$ID_COL]} = $numIndivs;
		$numIndivs++;
	}
}
close(IN);

for($chr = $chrStart; $chr <= $chrStop; $chr++) {
	$numIndivs = 0;
	my $infile = sprintf("ALL.chr$chr.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf", $chr);
	my $outfile = sprintf("1000Genomes.genotypes.$chr.txt", $chr);

	open(IN, "$infile") or die("ERROR: Could not open $infile");
	open(OUT, ">$outfile") or die("ERROR: Could not open $outfile");

	while(<IN>) {
		$line = $_;

		if(defined $line) {
			chomp($line);

			if($line =~ m/^##/) {  # Header line
				next;
			}
			elsif($line =~ m/^#CHROM/) { # Column information
				@temp = split("\t", $line);		
				$numCols = @temp;	
	
				for($i = $NUM_NONINDIV_COL; $i < $numCols; $i++) {
					$indivIDs[$numIndivs] = $temp[$i];
					$numIndivs++;	

					print OUT $temp[$i], " "; 
				} 				
				print OUT "\n";
				for($i = $NUM_NONINDIV_COL; $i < $numCols; $i++) {
					print OUT $indivInfo[$indivInfoH{$temp[$i]}][$POP_COL], " ";
				}
				print OUT "\n";
			}
			else {
				
			
				($dChr, $dPos, $dID, $dRef, $dAlt, $dQual, $dFilter, $dInfo, $dFormat, @dAlleles) = split("\t", $line);

				if($dFilter ne "PASS") {
					next;
				}

				if(!($dInfo =~ m/VT=SNP/)) {
					next;
				}
			
				if(!(exists $nt{$dRef}) or !(exists $nt{$dAlt})) { # If either the reference or the alternate is not a nt {	
					next;
	#				@tempRefs = split(",", $dRef);
	#				@tempAlts = split(",", $dAlt);
	#				$nRefs = @tempRefs;
	#				$nAlts = @tempAlts;			
				
	#				if($nRefs == 1 and $nAlts == 1) { # Site is not biallelic for nucleotides
	#					next;
	#				}
	#				elsif((exists $nt{$dRef}) and ($nAlts > 1)) { # Site is not biallelic
	#					next;
	#				}
	#				elsif(($nRefs > 1) and (exists $nt{$dAlt})) { # Site is not biallelic
	#					next;
	#				}
	#				elsif(($nRefs > 1) and ($nAlts > 1)) { # Site is not biallelic
	#					next;
	#				}	
	#				elsif($dRef eq "N") {
	#					if($nAlts != 2) {
	#						next;
	#					}
						
	#					$flagBad = 0;
	#					for($i = 0; $i < $nAlts; $i++) {
	#						if(!(exists $nt{$tempAlts[$i]})) { # Not a biallelic site of nucleotides
	#							$flagBad = 1;
#							} 
#						}

#						if($flagBad) {
#							next;
#						}
#					}
#					elsif($dAlt eq "N") {
#						if($nRefs != 2) {
#							next;
#						}

#						$flagBad = 0;
     #                                           for($i = 0; $i < $nRefs; $i++) {
    #                                                    if(!(exists $nt{$tempRefs[$i]})) { # Not a biallelic site of nucleotides
   #                                                             $flagBad = 1;
  #                                                      }
 #                                               }

#						if($flagBad) {
#							next;
#						}
#					}
				}
				
				print OUT "$dChr $dID $dPos $dRef $dAlt";

				for($i = 0; $i < $numIndivs; $i++) {
					@alleleInfo = split(":", $dAlleles[$i]);

					if($alleleInfo[0] eq "0|0") {
						print OUT " 1";
					}
					elsif($alleleInfo[0] eq "0|1") {
						print OUT " 2";
					}
					elsif($alleleInfo[0] eq "1|0") {
                                                print OUT " 2";
                                        }
					elsif($alleleInfo[0] eq "1|1") {
                                                print OUT " 3";
                                        }
					elsif($alleleInfo[0] eq ".|.") {
                                                print OUT " 0";
                                        }
					else {
						print "ERROR: strange genotype $dAlleles[$i]\n";
						die;
					}
				}
				print OUT "\n";

			#	die;			
			}
		}
	}

	close(IN);
	close(OUT);
}
