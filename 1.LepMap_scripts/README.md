Care for LepMap path and path to files as input to scripts


Command to change OrderMarkers output into more easily usable table 
`awk -vFS="\t" -vOFS="\t" '(NR==FNR){s[NR-1]=$0}(NR!=FNR){if ($1 in s) $1=s[$1];print}' LM_new_f.call.snps.txt order_50.txt  | tail -n +4 | cut -f 1,2,3,4 > order_50.perf`
