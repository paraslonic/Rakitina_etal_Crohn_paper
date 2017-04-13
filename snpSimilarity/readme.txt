requirements: bowtie2, samtools, varscan, R with libraries: gplots, reshape2, RColorBrewer
input data: read files in fastq format (not included) and assemblies in fasta format (included)
output: Figure1 (heatmap of snps calculated by reads to assembly mapping. SNPs on heatmap are normalized to covered genome length)

comments: RUN.sh ignores readsAssemblyMap, because it requires not included data

