#!/usr/bin/perl -w
#use strict;
use Bio::SeqIO;

my $infile = $ARGV[0];
my $outfile = $ARGV[1];
my $out = Bio::SeqIO->new(-file => ">$outfile" ,'-format' => 'Fasta');

print $infile . "\n";
$in  = Bio::SeqIO->new(-file => "$infile" , '-format' => 'GenBank');
while(my $seq = $in->next_seq() ){
		print $seq->display_id();
		$out->write_seq($seq)
}






