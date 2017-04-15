#!/usr/bin/perl

###############################################################################
##	converts genbank file to table
#	23.08.13 	A. Manolov

use Bio::SeqIO;
use File::Basename;

print "name\tid\tprod\tcontig\tstart\tstop\tprok\n";


my $fin = shift or die "[genbank file]";

$name = fileparse($fin, qr/\.[^.]*/);

$in = Bio::SeqIO->new(-file => $fin, -format => 'Genbank');
$id = 1;
while (my $seq = $in->next_seq()) {
	$contig = $seq->display_id();
	my @cds = grep { $_->primary_tag eq 'CDS' } $seq->get_SeqFeatures();
	foreach my $feature (@cds) {
		if($feature->has_tag('translation') ){
			my @fseq = $feature->get_tag_values('translation');
			my $prod = ($feature->get_tag_values('product'))[0];
			my $ltag = ($feature->get_tag_values('locus_tag'))[0];
			print "$name\t$id\t$prod\t$contig\t".$feature->start."\t".$feature->end."\t$ltag\n";
			$id++;
		}
	}
}
