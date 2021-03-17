# BAM-based method, per-gene analysis

Heatmap shows P18257_136 clustering with 189 day samples, not 126 day
* Doesn't seem to be reflected on the PCA?
* Regenerating the figures without the controls gives the shows the same pattern in both PCA and heatmap

```
for i in *tsv; do UP=${i%.*}; tail -n +2 $i | awk '$3 > 0' > ${UP}.DormantDown.tsv; tail -n +2 $i | awk '$3 < 0' > ${UP}.DormantUp.tsv; done
```

## Significant results, adjusted p < 0.05

			Control
	0d	49d	56d	72d	91d	126d	189d
0d	.	.	.	.	.	.	.
49d	14225	.	.	.	.	.	.
56d	12780	9270	.	.	.	.	.
72d	13437	8244	4625	.	.	.	.
91d	13226	7748	6486	4598	.	.	.
126d	13615	8063	7897	4678	999	.	.
189d	13012	7552	8385	6845	3923	2579	.

Comparing 0d to all other samples, 17,469 unique genes were found
* `cut -f1 0_*/*_SigResults.tsv | grep "TRINITY" | sort | uniq > veg_vs_rest.lst`

2,930 of these are not found in any of the resting-vs-resting comparisons
* `cut -f1 {1,4,5,7,9}*/*SigResults.tsv | grep "TRINITY" | sort | uniq > rest_vs_rest.lst`
* `comm -2 -3 veg_vs_rest.lst rest_vs_rest.lst > veg_vs_rest_only.lst`


## Upregulated

			Control
	0d	49d	56d	72d	91d	126d	189d
0d	.	.	.	.	.	.	.
49d	6875	.	.	.	.	.	.
56d	5967	4833	.	.	.	.	.
72d	6318	4205	2159	.	.	.	.
91d	6328	4004	3180	2362	.	.	.
126d	6571	4088	3849	2332	399	.	.
189d	6386	3699	3874	3193	1596	1005	.

## Downregulated

			Control
	0d	49d	56d	72d	91d	126d	189d
0d	.	.	.	.	.	.	.
49d	7350	.	.	.	.	.	.
56d	6813	4437	.	.	.	.	.
72d	7119	4039	2466	.	.	.	.
91d	6898	3744	3306	2236	.	.	.
126d	7044	3975	4048	2346	600	.	.
189d	6626	3853	4511	3652	2327	1574	.
