mkdir -p blastout_lf
mkdir -p blastout_jj
for f in ref/*.fasta
do 
	name=$(basename $f .fasta)
	blastn -query plf82.fasta -subject $f -outfmt '6 std qlen slen qcovs' -evalue 1e-5 > blastout_lf/$name
	blastn -query pJJ1886_4.fasta -subject $f -outfmt '6 std qlen slen qcovs' -evalue 1e-5 > blastout_jj/$name
done

Rscript calcCoverage.r
