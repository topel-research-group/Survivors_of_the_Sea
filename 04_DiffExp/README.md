# Readme file for final results

Note - Tables of DE results for [R05](R05_FINAL/LFCs.R05.lfc0.585.reformatted.tsv) and 
[GF04](GF04_FINAL/LFCs.GF04.lfc0.585.reformatted.tsv) use the old, work-in-progress gene naming conventions.  
For naming conventions matching the publication and the NCBI upload, see [NameConversion.tsv](NameConversion.tsv)

## Upregulation

* Genes significantly upregulated in at least one time point and strain: 6,975
`cut -f1 *_FINAL/lfc_0.585/0_*d/0_*d_SigResults.*d_Up.tsv | grep -v "base" | sort | uniq | wc -l`

* Genes significantly upregulated in all twelve comparisons: 1,329
`cut -f1 *_FINAL/lfc_0.585/0_*d/0_*d_SigResults.*d_Up.tsv | grep -v "base" | sort | uniq -c | grep -c " 12 "`

## Downregulation

* Genes significantly downregulated in at least one time point and strain: 7,644
`cut -f1 *_FINAL/lfc_0.585/0_*d/0_*d_SigResults.*d_Down.tsv | grep -v "base" | sort | uniq | wc -l`

* Genes significantly downregulated in all twelve comparisons: 1,911
`cut -f1 *_FINAL/lfc_0.585/0_*d/0_*d_SigResults.*d_Down.tsv | grep -v "base" | sort | uniq -c | grep -c " 12 "`

###

## Significantly differentially expressed genes in resting stages vs. vegetative

| Strain |     49d    |     56d    |     72d    |     91d    |    126d    |    189d    |
|--------|------------|------------|------------|------------|------------|------------|
| R05    | 3,982 up   | 3,701 up   | 4,004 up   | 4,081 up   | 4,449 up   | 4,151 up   |
|        | 4,693 down | 4,456 down | 4,836 down | 4,917 down | 5,187 down | 4,313 down |
|--------|------------|------------|------------|------------|------------|------------|
| GF04   | 3,655 up   | 3,004 up   | 3,254 up   | 3,129 up   | 3,340 up   | 3,167 up   |
|        | 4,501 down | 3,782 down | 4,264 down | 3,934 down | 4,228 down | 3,870 down |
|--------|------------|------------|------------|------------|------------|------------|
| Shared | 2,525 up   | 2,001 up   | 2,266 up   | 2,204 up   | 2,445 up   | 2,146 up   |
|        | 3,280 down | 2,766 down | 3,174 down | 3,054 down | 3,257 down | 2,638 down |

```
for i in 49 56 72 91 126 189; do
	echo "R05: 0 vs. ${i}"
	UP=$(grep -c "Sm_" R05_FINAL/lfc_0.585/0_${i}d/0_${i}d_SigResults.${i}d_Up.tsv)
	echo $UP
	DOWN=$(grep -c "Sm_" R05_FINAL/lfc_0.585/0_${i}d/0_${i}d_SigResults.${i}d_Down.tsv)
	echo $DOWN
	echo ""

	echo "GF04: 0 vs. ${i}"
	UP=$(grep -c "Sm_" GF04_FINAL/lfc_0.585/0_${i}d/0_${i}d_SigResults.${i}d_Up.tsv)
	echo $UP
	DOWN=$(grep -c "Sm_" GF04_FINAL/lfc_0.585/0_${i}d/0_${i}d_SigResults.${i}d_Down.tsv)
	echo $DOWN
	echo ""

	echo "Shared: 0 vs. ${i}"
	UP=$(grep "Sm_" *_FINAL/lfc_0.585/0_${i}d/0_${i}d_SigResults.${i}d_Up.tsv | cut -f1 | cut -f2 -d':' | sort | uniq -c | grep " 2 " | wc -l)
	echo $UP
	DOWN=$(grep "Sm_" *_FINAL/lfc_0.585/0_${i}d/0_${i}d_SigResults.${i}d_Down.tsv | cut -f1 | cut -f2 -d':' | sort | uniq -c | grep " 2 " | wc -l)
	echo $DOWN
	echo ""
done
```
