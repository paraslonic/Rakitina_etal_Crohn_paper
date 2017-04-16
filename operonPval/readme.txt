Calculates p-value of operon overrepresentation by random reshufflings and by using Poisson distibution

requirements: perl with bioperl, R

input data: annotation of lf82 strain Escherichia_coli_LF82_uid161965.gbf and the result of operon predicion by DOOR door_output.txt, ongenome.tab and table.txt are the output files of ogEnrichment scripts.

output: statsignificant_operons.pdf, statsignificant_operon_info.txt, comparison between observed operon enrichment and the expected by chance;
statsignificant_operon_info_Poisson.txt - operons with adjusted p-value < 0.05 (test assumin Poisson distribution)
