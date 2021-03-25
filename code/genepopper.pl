#!/usr/bin/perl
use strict; use warnings; use Getopt::Long;

### This program looks at a genepop-formatted output file from the populations
### step of stacks.  It outputs the information you need to do another
### populations step.

### If you give it a cut-off value, it will also write out a populations file
### including only those individuals with genotypes for that proportion of SNPs
### If you give it an -r parameter, it will write a new sbatch file to run
### populations with that -r parameter

###-------------------------------------------------------------------------###
###                                Initialize                               ###
###-------------------------------------------------------------------------###

my ($input_file, $help, $cutoff, $id, $run_id, $r_parameter); my %data_hash;
GetOptions ('input|i:s' => \$input_file,
            'run_id|u:s' => \$run_id,
	    'help|h' => \$help,
	    'cutoff|c:f' => \$cutoff,
        'output|o:s' => \$output_file,
	    'r_parameter|r:f' => \$r_parameter);

if ($help) {die "This script looks at the genepop-formatted output file from
the populations step of stacks and outputs some important information.
  -i provides the genepop file to look at.
  -c (optional) writes a stacks populations map including only those
     individuals with genotypes for >= the provided proportion of SNPs.
  -o (optional; required if -c supplied) filepath for output population map file.
  -r (optional) writes a qsub file to run stacks' population program with
     the -r parameter given.
  -R (optional) the run ID for out files from the qsub file.
  -h displays this message.\n"}

open (INPUT, '<', $input_file) or die "
Unable to open the genepop file $input_file.\n";

if (defined $cutoff) {
  unless (defined $output_file) {die "
If you use the cutoff option, you must specify an output file path with -o.\n"}
  unless (($cutoff <= 1) and ($cutoff >= 0)) {die "
The cutoff value must be a proportion between 0 and 1.\n"}
  if (-e "$output_file") {die "
Cannot write to $output_file: this file already exists.\n"}
  open (POPULATIONS, '>', "$output_file") or die "
Cannot open the file $output_file.\n";
}

if (defined $r_parameter) {
  unless (defined $id) {die "
If you use the r_parameter option, you must specify an id with -n.\n"}
  unless (($r_parameter <= 1) and ($r_parameter >=0)) {die "
The -r parameter must be a proportion between 0 and 1.\n"}
  if (-e "run_populations_$id.sbatch") {die "
Cannot write a file named run_populations_$id.sbatch:
this file already exists.\n"}
  open (SBATCH, '>', "run_$run_id" . "_pop$id.qsub") or die "
Cannot open the file run_populations_$id.sbatch.\n"
}

###-------------------------------------------------------------------------###
###                        Read in the genepop file                         ###
###-------------------------------------------------------------------------###
<INPUT>; <INPUT>; #Discards the first two lines
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
print STDOUT "The total number of SNPs is $length.
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
  print STDOUT "$counter individuals have a genotype for $percentage% of these SNPs.\n"
}
### END TESTED PART ###

###-------------------------------------------------------------------------###
###                     Write out the populations file                      ###
###-------------------------------------------------------------------------###
if (defined $cutoff) {
  foreach (keys %data_hash) {
    if ($data_hash{$_} >= $cutoff) {
      print POPULATIONS "$_\t1\n"
    }
  }
}

###-------------------------------------------------------------------------###
###                        Write out the sbatch file                        ###
###-------------------------------------------------------------------------###
if (defined $r_parameter) {
  print SBATCH "#!/bin/tcsh
#PBS -l nodes=1:hima:ppn=1
#PBS -l walltime=72:00:00
#PBS -N S$run_id
#PBS -j oe
#PBS -m ae
#PBS -M john.h.boyle\@gmail.com

# This is step $run_id in the stacks pipeline

cd \$PBS_O_WORKDIR
module load gcc/7.2.0

/sciclone/home2/jhboyle/programs/stacks/bin/populations -P ref_map_output \\
  -b 1 \\
  -O . \\
  -M population_map_$id \\
  -r $r_parameter \\
  --genepop\n\n"
}

###-------------------------------------------------------------------------###
###                               Denitialize                               ###
###-------------------------------------------------------------------------###
close POPULATIONS; close SBATCH; close INPUT
