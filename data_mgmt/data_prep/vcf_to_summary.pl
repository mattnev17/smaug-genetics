#!/usr/bin/perl
################################################################################
# Get mutation type/context information from input vcfs
# To avoid the possibility of modifying raw vcf data, includes option to copy
# raw vcfs to a working project directory
################################################################################
use strict;
use warnings;
use POSIX;
use File::Basename;
use File::Path qw(make_path);
use File::Find::Rule;
use FindBin;
use YAML::XS 'LoadFile';
use feature 'say';

my $relpath = $FindBin::Bin;
my $configpath = dirname(dirname($relpath));
my $config = LoadFile("$configpath/_config.yaml");

my $mac = $config->{mac};
my $parentdir = $config->{parentdir};
my $vcftoolsdir = $config->{vcftoolsdir};
my $rawvcfdir = $config->{rawvcfdir};
my $rawvcfext = $config->{rawvcfext};

use lib "$FindBin::Bin/../lib";
use SmaugFunctions qw(forkExecWait getRef);

# Default summary directory
my $outdir="$parentdir/summaries";
make_path($outdir);

# Running script with argument 'copy'
# copies VCFs from another directory to avoid overwriting raw data
my $makecopy='';
if(defined($ARGV[0])){
  $makecopy = $ARGV[0];
}

# Overwrite the MAC value in YAML config with command line arg
if(defined($ARGV[1])){
  $mac = $ARGV[1];
}

# Specify project folder for VCFs
my $vcfdir="$parentdir/vcfs";
make_path($vcfdir);

################################################################################
# Reads raw vcfs from original location, annotates with necessary info,
# and copies to input directory. Original vcfs are preserved and not modified.
################################################################################
if ($makecopy eq "copy") {
	my @rawvcfs = File::Find::Rule->file()
                            ->name("*.$rawvcfext")
                            ->maxdepth(3)
                            ->in($rawvcfdir);

	foreach my $rawvcf (@rawvcfs) {
		my $filename=fileparse($rawvcf);

    if($filename !~ /chrX/){
      my @parts = split(/\.vcf.gz/, $filename);
      my $basename = $parts[0];
      my @nameparts = split(/\./, $basename);
      my $i = $nameparts[0];
      $i =~ s/chr//g;

      my $newvcf = "$parentdir/vcfs/$basename.ma.aa.$mac.vcf.gz";
      my $fasta = "$parentdir/reference_data/Homo_sapiens.GRCh38_15/chr$i.fasta.gz";

      # last command adds Motif and Category info fields
      my $infoparse = "perl $relpath/fill_motif.pl -a $fasta";

      # pipe commands and execute
      my $pipe;
      if($mac eq "singletons"){
        my $filter = "bcftools view -i 'AC==1' -f PASS";
        $pipe = "$filter $rawvcf | $infoparse | bgzip -c > $newvcf";
      } elsif($mac eq "common"){
        my $filter = "bcftools view -i 'AC>=10' -f PASS";
        $pipe = "$filter $rawvcf | $infoparse | bgzip -c > $newvcf";
      } elsif($mac eq "rare"){
        my $filter = "bcftools view --max-af 0.0001404 -f PASS";
        $pipe = "$filter $rawvcf | $infoparse | bgzip -c > $newvcf";
      } elsif($mac eq "all"){
        my $filter = "bcftools view -f PASS";
        $pipe = "$filter $rawvcf | $infoparse | bgzip -c > $newvcf";
      }

      print STDERR "Input file: $rawvcf\n";
      print STDERR "Writing to: $newvcf...\n";
      forkExecWait($pipe);
      print STDERR "Done\n";
    }
	}

  print STDERR "Operation complete\n";
}

################################################################################
# Scans input directory for vcfs and outputs summary file to get full genome summary file
################################################################################
my $script = 1;

if ($script==1){
  my @vcfs = File::Find::Rule->file()
                            ->name("*$mac.vcf.gz")
                            ->maxdepth(1)
                            ->in($vcfdir);
  my $header = "\"CHR\tPOS\tREF\tALT\tAC\tAN\tMotif\tCategory\"";
  my $summout = "$outdir/$mac.full.summary";
  my $headercmd = "echo $header > $summout";
  forkExecWait($headercmd);

	foreach my $file (@vcfs) {
    print STDERR "Getting summary for $file...\n";
    unless(-e "$file.tbi"){
      print STDERR "$file not indexed--indexing now...";
      my $tabixcmd = "tabix -p vcf $file";
      forkExecWait($tabixcmd);
      print STDERR "Done\n";
    }

    my $filename=fileparse($file);
    my @parts = split(/\./, $filename);
    my $chr = $parts[0];
    $chr =~ s/chr//;
    my $bcfquery;
    my $outputcols;

    $bcfquery = "bcftools query -r chr$chr";
    $outputcols = "'%CHROM\t%POS\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\t%INFO/Motif\t%INFO/Category\n'";

    my $cmd = "$bcfquery -f $outputcols $file >> $summout";
		forkExecWait($cmd);
	}
  print STDERR "Operation complete\n";
}
