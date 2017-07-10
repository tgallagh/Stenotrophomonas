This folder contains output from the edgeR package to find statistically significant genes (.genes. ext) and bins (.bins. ext)

To add gene names to edgeR output table:

 join -1 1 -2 1 -o 2.1,2.2,2.3,2.4,2.5,2.6,1.2 <(sort -k 1b,1 Id_name_convert.txt) <(sort -k 1b,1 de7_9_filter.genes.txt) | column -t > de7_9filter.genes.names.txt 
echo -e "gene\tlogFC\tlogCPM\tLR\tPValue\tFDR\tKEGG\tName" | cat - de5_7filter.genes.names.txt > de5_7filter.genes.names.final.txt 

### Summary of edgeR genes output
## pH 5 versus pH 7:
using FDR < 0.05
Numer of  genes upregulated: 98
	Number hypothetical : 40
	Number of operons : 7 
Number of genes  downregulated: 86
	Number of hypothetical: 18 
	Number of operons: 3

## pH 9 versus pH 7:
using FDR < 0.05 
Number of genes upregulated: 5
	Number hypothetical: 0
	Number of operons: 1
Number of genes downregulated: 2
	Number of hypoethetical: 0
	Number of operons: 0

### Overlaps between pH 5 and pH 9 comparisons
fig|6666666.230262.peg.2616 His Repressor is downregulated in pH 5 but upregulated in pH 9
fig|6666666.230262.peg.3167 Succinate-semialdehyde_dehydrogenase_ is downregulated in both
fig|6666666.230262.peg.3404 Galactosamine-6-phosphate_isomerase_(EC_5.3.1.-) is upreg in both


	

