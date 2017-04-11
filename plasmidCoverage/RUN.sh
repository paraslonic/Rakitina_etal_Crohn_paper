mkdir blastout
for f in ref/*.fasta
do 
	name=$(basename $f .fasta)
	blastn -query plf82.fasta -subject $f -outfmt '6 std qlen slen qcovs' -evalue 1e-5 > blastout/$name
done

Rscript calcCoverage.r
