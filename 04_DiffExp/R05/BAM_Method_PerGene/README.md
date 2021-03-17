# BAM-based method, per-gene analysis

Samples seem to cluster together
* 189d samples a little separated on PCA?

```
for i in *tsv; do UP=${i%.*}; tail -n +2 $i | awk '$3 > 0' > ${UP}.DormantDown.tsv; tail -n +2 $i | awk '$3 < 0' > ${UP}.DormantUp.tsv; done
```

## Significant results, adjusted p < 0.05

			Control
	0d	49d	56d	72d	91d	126d	189d
0d	.	.	.	.	.	.	.
49d	14507	.	.	.	.	.	.
56d	14316	6770	.	.	.	.	.
72d	14893	6682	3264	.	.	.	.
91d	15254	4191	3412	334	.	.	.
126d	15797	5927	7760	5123	3068	.	.
189d	13859	9936	9828	8711	7158	4773	.

Comparing 0d to all other samples, 19,054 unique genes were found
* `cut -f1 0_*/*_SigResults.tsv | grep "TRINITY" | sort | uniq > veg_vs_rest.lst`

4,030 of these are not found in any of the resting-vs-resting comparisons
* `cut -f1 {1,4,5,7,9}*/*SigResults.tsv | grep "TRINITY" | sort | uniq > rest_vs_rest.lst`
* `comm -2 -3 veg_vs_rest.lst rest_vs_rest.lst > veg_vs_rest_only.lst`

## Upregulated

			Control
	0d	49d	56d	72d	91d	126d	189d
0d	.	.	.	.	.	.	.
49d	6997	.	.	.	.	.	.
56d	6760	3390	.	.	.	.	.
72d	7143	3527	1724	.	.	.	.
91d	7251	2053	1570	84	.	.	.
126d	7605	2999	3927	2511	1656	.	.
189d	7112	5158	5145	4581	3814	2590	.

## Downregulated
			Control
	0d	49d	56d	72d	91d	126d	189d
0d	.	.	.	.	.	.	.
49d	7510	.	.	.	.	.	.
56d	7556	3380	.	.	.	.	.
72d	7750	3155	1540	.	.	.	.
91d	8003	2138	1842	250	.	.	.
126d	8192	2928	3833	2612	1412	.	.
189d	6747	4778	4683	4130	3344	2183	.
