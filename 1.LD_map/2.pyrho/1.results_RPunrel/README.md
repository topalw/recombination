The folders contain results (true for 1. and 2. where results are small) 
or example results (true for 3. where results for all SS would be too big)


For make windows the first step is making the total pyrho file for all scaffolds like this 
```
for FILE in *Super-Scaffold*
do
SS=`echo ${FILE} | cut -d '_' -f 3,4`
ALLPY=`echo ${FILE} | cut -d '_' -f 1,2`
awk '{print ss "\t" $1 "\t" $2 "\t" $3 }' ss="$SS" ${FILE} >> "${ALLPY}_total.bed"
done
```
