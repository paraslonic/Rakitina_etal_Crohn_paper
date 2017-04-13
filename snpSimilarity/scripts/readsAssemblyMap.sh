# requirements: bowtie2 samtools varscan

outdir=../readsAssembliesMap
rdir=../reads
adir=../assemblies

mkdir -p $outdir
mkdir -p $outdir/bwt

rm covered alength snps

for a in `ls -1 $adir/*.fna`; do
	aname=$( basename $a .fna )
	echo $aname
	bowtie2-build $adir/$aname.fna $outdir/bwt/$aname
done
 
for r in `ls -1 $rdir/*.fastq`; do
	rname=$( basename $r .fastq )
	for a in $adir/*.fna; do
		aname=$( basename $a .fna )
		out=$outdir/${rname}Reads2${aname}
		bowtie2 -x $outdir/bwt/$aname -U $r -p 22 2> $out.log | samtools view -bS - > ${out}_us.bam 
		samtools sort ${out}_us.bam $out
		samtools index $out.bam
		rm ${out}_us.bam
		samtools mpileup -f $a $out.bam > $out.mpileup
		varscan mpileup2snp $out.mpileup -p-value 1e-5 --min-var-freq 0.9 --min-coverage 4 > $out.snp
		wc -l $out.snp | cut -d' ' -f1 | xargs -I {} printf "$rname\t$aname\t{}\n" >> SNP_COUNT
		awk '/^>/ {if (seqlen){print seqlen};seqlen=0;next; }{ seqlen = seqlen +length($0)}END{print seqlen}' $a | paste -sd+ | bc | xargs -I {} printf "$aname\t{}\n" >> LENGTH 
		samtools depth $out.bam | awk '{ if ($3 > 4) print }' | wc -l | cut -d' ' -f1 | xargs -I {} printf "$rname\t$aname\t{}\n" >> COVERED
	done
done
