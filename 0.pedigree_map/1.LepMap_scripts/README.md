Care for LepMap path and path to files as input to scripts


Command to change OrderMarkers output into more easily usable table 
``` 
for file in *gp.txt; do awk -vFS="\t" -vOFS="\t" '(NR==FNR){s[NR-1]=$0}(NR!=FNR){if ($1 in s) $1=s[$1];print}' Call2_tf1_maf05_miss05_1kb_pruned_f.call.snps.txt ${file} | tail -n +4 | cut -f 1,2,3,4 > ${file}.perf; done 
```
