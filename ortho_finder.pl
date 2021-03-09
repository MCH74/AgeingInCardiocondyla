#!/usr/bin/perl
use strict; use warnings;use IO::Handle;use diagnostics;

##find reciprocal best blast hits

my ($blast1, $blast2) = @ARGV; #two blast output files

open OUTPUT, '>', "output.txt" or die $!;
open ERROR,  '>', "error.txt"  or die $!;

STDOUT->fdopen( \*OUTPUT, 'w' ) or die $!;
STDERR->fdopen( \*ERROR,  'w' ) or die $!;

open ONE, "<$blast1";
chomp(my @one = <ONE>);

#my $line = join '', @one;
#my @lines1 = split "Query=", $line;

my $dmel_gene;
my %dmel_match;
my %dmel_score;
my %alt_name;

foreach my $one (@one){
   if($one =~ /^Query=\s+(\S+)\s+/){
	   $dmel_gene = $1;
   }
   elsif($one =~ /(FBgn\d+)/){
	   $alt_name{$dmel_gene} = $1;
   }
   elsif($one =~ /^\s+gnl\|Cobs_1\.4\|(Cobs_\d+)-mRNA-\d+\s+\S+\s+(\S+)\s*$/){
      unless(exists $dmel_match{$dmel_gene}){
         $dmel_match{$dmel_gene} = $1;
         $dmel_score{$dmel_gene} = $2;
      }
   }
}



open TWO, "<$blast2";
chomp(my @two = <TWO>);

#$line = join '', @two;
#my @lines2 = split "Query=", $line;

my $cobs_gene;
my %cobs_match;

foreach my $two (@two){
   if($two =~ /^Query=\s+gnl\|Cobs_1\.4\|(Cobs_\d+)\-mRNA/){
      $cobs_gene = $1;
   }
   elsif($two =~ /^\s+(FBpp\S+)\s+/){
      unless(exists $cobs_match{$cobs_gene}){
         $cobs_match{$cobs_gene} = $1;
      }
  }
}

open OUT, ">ortholog_pairs.txt";

print OUT "Dmel\tCobs\tDmelScore\tAltName\n";

foreach my $key (sort keys %dmel_match){
   if($key eq $cobs_match{ $dmel_match{$key}  } && $dmel_score{$key} < 0.001){
      print OUT "$key\t$dmel_match{$key}\t$dmel_score{$key}\t$alt_name{$key}\n";
   }
}

