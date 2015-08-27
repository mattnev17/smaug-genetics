#!/usr/local/bin/perl

##############################################################################
# Old version of getNonMut.pl, using sequence conservation scores as an 
# additional covariate for logistic regression.  Drops too many sites to be
# useful
##############################################################################

use strict;
use warnings;
use POSIX;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use File::Path qw(make_path);
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use List::MoreUtils 'pairwise';
use Cwd;
use Benchmark;
use Tie::File;

# Set options and inputs
my $wdir=getcwd;
my $parentdir="/net/bipolar/jedidiah/mutation";

my $baseopt;
my $chr;
my $categ;
my $bw = 100;
my $f_covs = "$parentdir/output/logmod_data/${bw}kb_mut_cov2.txt";

GetOptions ('b=s'=> \$baseopt,
'chr=s'=> \$chr,
'categ=s' => \$categ,
'bw=i' => \$bw,
'covs=s' => \$f_covs) or pod2usage(1);

my $b1;
my $b2;
if($baseopt eq "AT"){
	$b1="A";
	$b2="T";
} else {
	$b1="C";
	$b2="G";
}

my $mask_flag=0;
my $adj=2;
my $binwidth=$bw*1000;

my $nextchr;
if ($chr<22) {
	$nextchr=$chr+1;
} elsif ($chr==22) {
	$nextchr="X";
}

my $subseq=2;
if ($adj!=0) {
	$subseq = $adj*2+1;
}

# initialize output file
my $outfile = "$parentdir/output/logmod_data/chr${chr}_${categ}_sites.txt";
open(OUT, '>', $outfile) or die "can't write to $outfile: $!\n";

# initialize covariate data
open my $covs, '<', $f_covs or die "can't open $f_covs: $!";

# initialize singleton file
my $f_positions = "$parentdir/output/logmod_data/chr${chr}_${categ}_pos_examples.txt"; #main line for full processing
open my $positions, '<', $f_positions or die "can't open $f_positions: $!";

# initialize phastCons data
my $f_cons = "$parentdir/reference_data/chr$chr.phastCons46way.primates.wigFix";
open my $cons, '<', $f_cons or die "can't open $f_cons: $!";

# Get reference sequence
my $seq=&getRef();
my $altseq=$seq;
$altseq =~ tr/ACGT/TGCA/;

# my $seqlength=length($seq);
# print "seqlength of chr$chr: $max\n"; #<-used to validate that getRef() returns correct seq length

my $printheader=0;
if($printheader==1){
	print OUT "CHR \t POS \t BIN \t Sequence \t mut \n"; #<-add header to output, if needed
}

# Create hash keyed by Chr/Bin pairs, with row of PCs as value
print "Indexing chr${chr} covariate data...\n";
our %hash=();
while (<$covs>){
	chomp;
	my @line=split(/\t/, $_);
	my $key=join("\t", @line[0 .. 1]);
	my $pcs=join("\t", @line[2 .. $#line]);
	
	$hash{$key}=$pcs;
}

my $key=join("\t", 20, 100);
print "$hash{$key}\n";

# Create hash keyed by singleton positions, with input line as value
print "Indexing chr${chr}: ${categ} singleton positions...\n";
my @POS;
my %poshash=();
while (<$positions>) {
	chomp;
	my @line=split(/\t/, $_);
	my $key=$line[2];
	push (@POS, $key);
	
	$poshash{$key}=$_;
}

# Reads conservation score data, writes line of output if site meets criteria
print "Writing chr${chr}: ${categ} data file...\n";
my $pos;
while (<$cons>){
	chomp;
	
	# Reads descriptor line of wigFix file and extracts starting position
	if($_ =~ /fixedStep/){
		my @head = split(/[=,\s]+/, $_);
		$pos=$head[4];
	}elsif($pos>1){
		my $base = substr($seq, $pos-1, 1);

		# if statement evaluated for sites not in singleton list; elsif statement evaluated for singletons
		if(($base =~ /$b1|$b2/) & (!exists $poshash{$pos})){
			# push (@POS, $pos); # add position to exclusion list
			my $localseq = substr($seq, $pos-$adj-1, $subseq);
			my $altlocalseq = reverse substr($altseq, $pos-$adj-1, $subseq);
			my $bin = ceil($pos/$binwidth);
			
			# Coerce local sequence info to format used in R
			my $sequence;
			if(substr($localseq,$adj,1) lt substr($altlocalseq,$adj,1)){
				$sequence = $localseq . '(' . $altlocalseq . ')';
			} else {
				$sequence = $altlocalseq . '(' . $localseq . ')';
			}
			
			my $key2=join("\t", $chr, $bin);
			
			# write line if site has non-N context and exists in covariate data
			if (($sequence !~ /N/) && (exists $hash{$key2})) {
				# my $key2=join("\t", $chr, $bin);
				# my $hv=$hash{$key2};
				# print "$key2\t$hv\n";
				my $covs=&updateCovs($chr, $bin, $pos);
				# my $covs=$hv;
				print OUT "$chr\t$bin\t$pos\t$sequence\t 0 \t$covs\t$_\n";	 
			} 
		} elsif((exists $poshash{$pos})){
			my $bin = ceil($pos/$binwidth);
			# print "$bin\t$pos\n";
			# my $covs=&updateCovs($chr, $bin, $pos);
			my $key2=join("\t", $chr, $bin);
			# my $covs=$hash{$key2};
			if(exists $hash{$key2}){
				print "$key2 exists \n";
				my $covs=&updateCovs($chr, $bin, $pos);
				print OUT "$poshash{$pos}\t$covs\t$_\n";
			}		
		}
		
		$pos++;
	}
}
print "Done\n";

sub getRef{
	my $f_fasta;
	if($mask_flag==1){
		$f_fasta = "$parentdir/reference_data/human_g1k_v37.mask.fasta";
	} else {
		$f_fasta = "$parentdir/reference_data/human_g1k_v37.fasta";
	}

	if (-e $f_fasta) {
		print "Using reference genome: $f_fasta\n";
	} else {
		print "Reference genome not found in parent directory. Would you like to download one? (y/n): ";
		my $choice = <>;
		chomp $choice;
		if ($choice eq "y") {
			my $dlcmd="wget -P $parentdir/ ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz";
			&forkExecWait($dlcmd);
			my $unzipcmd="gunzip $parentdir/human_g1k_v37.fasta";
			&forkExecWait($unzipcmd);
		} else {
			die "Please upload an appropriate reference genome to the parent directory\n";
		}
	}

	open my $fasta, '<', $f_fasta or die "can't open $f_fasta: $!";

	##############################################################################
	# Retrieve reference sequence for selected chromosome
	# -also returns symmetric sequence to be used in local sequence analysis
	##############################################################################

	print "Getting reference sequence for chromosome $chr...\n";

	my $seq;
	if($mask_flag==1){
		while (<$fasta>) {
			chomp;
			if (/^>$chr$/../^>$nextchr$/) {
				next if /^>$chr$/ || /^>$nextchr$/;
				$seq .=$_;
			}
		}
	} else {
		while (<$fasta>) {
			chomp;
			if (/>$chr /../>$nextchr /) {
				next if />$chr / || />$nextchr /;
				$seq .=$_;
			}
		}
	}
	
	return $seq;
}


sub updateCovs{
	
	my $CHR=shift;
	my $BIN=shift;
	my $pos=shift;
	
	# Get keys for current and neighboring bins
	my $c_linekey=join("\t", $CHR, $BIN);
	my $p_linekey=join("\t", $CHR, $BIN-1);
	my $n_linekey=join("\t", $CHR, $BIN+1);
	
	# my $pos=$line[2];
	my $posmin=$BIN*$binwidth-$binwidth;
	
	# Get relative position
	my $relpos = ($pos-$posmin)/$binwidth;

	# Calculate proportion 
	my $prop_c_bin = -abs($relpos-0.5)+1;
	
	# Get covariates of containing bin
	my $o_line = $hash{$c_linekey};
	
	# print "$c_linekey\n";
	# print "$o_line\n";
	
	my @c_feats = split(/\t/, $o_line);
	
	# Calculate covariates in current bin proportional to position
	foreach my $x (@c_feats) { $x = $x * $prop_c_bin; }

	my @sum;
	# Repeat for adjacent window
	if(($relpos-0.5<=0) && exists($hash{$p_linekey})){
		my $prop_p_bin = -$relpos+0.5;
		my @p_feats = split(/\t/, $hash{$p_linekey});
		foreach my $x (@p_feats) { $x = $x * $prop_p_bin; }
		@sum = pairwise { $a + $b } @c_feats, @p_feats;
		
	} elsif(($relpos-0.5>0) && exists($hash{$n_linekey})){
		my $prop_n_bin = $relpos-0.5;
		my @n_feats = split(/\t/, $hash{$n_linekey});
		foreach my $x (@n_feats) { $x = $x * $prop_n_bin; }
		@sum = pairwise { $a + $b } @c_feats, @n_feats;
	}
	
	my $covs=join("\t", @sum[0 .. $#sum]);
	print "$covs\n";
	return $covs;
}

