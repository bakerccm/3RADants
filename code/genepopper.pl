#!/usr/bin/perl
use strict; use warnings; use Getopt::Long;

### This code adapted from Jack Boyle's genepopper.pl by Chris Baker.

### This program looks at a genepop-formatted output file from the populations
### step of stacks.  It outputs the information you need to do another
### populations step.

### If you give it a cut-off value, it will also write out a populations file
### including only those individuals with genotypes for that proportion of SNPs

###-------------------------------------------------------------------------###
###                                Initialize                               ###
###-------------------------------------------------------------------------###

my ($input_file, $output_log_file, $cutoff, $output_popmap_file, $help); my %data_hash;
GetOptions ('input|i:s' => \$input_file,
        'output|o:s' => \$output_log_file,
        'cutoff|c:f' => \$cutoff,
        'popmap|p:s' => \$output_popmap_file,
        'help|h' => \$help);

if ($help) {die "This script looks at the genepop-formatted output file from
the populations step of stacks and outputs some important information.
  -i (required) filepath for genepop file to parse.
  -o (required) filepath for output log file.
  -c (optional) writes a stacks population map including only those
     individuals with genotypes for >= the provided proportion of SNPs.
  -p (optional; required if -c supplied) filepath for output population map file.
  -h display this message.\n"}

open (INPUT, '<', $input_file) or die "
Unable to open the genepop file $input_file.\n";

open (OUTPUT_LOG, '>', "$output_log_file") or die "
Unable to open the file $output_log_file.\n";

if (defined $cutoff) {
  unless (defined $output_popmap_file) {die "
If you use the cutoff option, you must specify an output file path with -o.\n"}
  unless (($cutoff <= 1) and ($cutoff >= 0)) {die "
The cutoff value must be a proportion between 0 and 1.\n"}
  open (OUTPUT_POPMAP, '>', "$output_popmap_file") or die "
Cannot open the file $output_popmap_file.\n";
}

###-------------------------------------------------------------------------###
###                        Read in the genepop file                         ###
###-------------------------------------------------------------------------###
<INPUT>; <INPUT>; # discards the first two lines
my $length;
while (<INPUT>) {
  chomp $_;
  if ($_ =~ m/pop$/) {next}
  my @rray = split (/\t/, $_);
  my $sample = shift (@rray);
  $sample =~ m/(.*),/;
  $sample = $1;
  $length = @rray;
  my $counter = 0;
  foreach (@rray) {
    if ($_ eq '0000') {
      $counter += 1
    }
  }
  $data_hash{$sample} = ($length - $counter) / $length
}

###-------------------------------------------------------------------------###
###                      Report a few basic statistics                      ###
###-------------------------------------------------------------------------###
my $individuals = keys %data_hash;
my $missing_data;
foreach (keys %data_hash) {
  $missing_data += $data_hash{$_}
}
my $print_this = 1 - $missing_data / $individuals;
print OUTPUT_LOG "Input file: $input_file
The total number of SNPs is $length.
The total number of individuals is $individuals.
The total proportion of missing data is $print_this\n";
foreach my $proportion (.1, .2, .5, .8, .9, 1) {
  my $counter = 0;
  my $percentage = $proportion * 100;
  foreach (keys %data_hash) {
    if ($data_hash{$_} >= $proportion) {
      $counter += 1
    }
  }
  print OUTPUT_LOG "$counter individuals have a genotype for $percentage% of these SNPs.\n"
}

###-------------------------------------------------------------------------###
###                     Write out the populations file                      ###
###-------------------------------------------------------------------------###
if (defined $cutoff) {
  foreach (keys %data_hash) {
    if ($data_hash{$_} >= $cutoff) {
      print OUTPUT_POPMAP "$_\t1\n"
    }
  }
}

###-------------------------------------------------------------------------###
###                               Denitialize                               ###
###-------------------------------------------------------------------------###
close OUTPUT_LOG; close OUTPUT_POPMAP; close INPUT
