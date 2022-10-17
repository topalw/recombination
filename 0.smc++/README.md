# smc++ implementation for the barn owls

I have been working a bit on these scripts that were pretty ancient and highly non-replicable. 
TO DO: 
- curate vcf2smc, estimate and slurm so that they can be applied to any input. *care what I assume about filenames*

The scipts work in sequence. 
First we use the 0.aux_scripts/0.1_prepare_smcpp_input.slurm ../1.vcfs/CH_miss10.vcf.gz more_than_100_snps.table which uses 
 vcftools/0.1.14; bcftools/1.15.1 and plink-ng/1.9.20200712 to prepare the file / scaffold and create the distinguished individual list and the csv of samples to feed to smc++. 
 Then we use smc v1.15.4.dev18+gca077da.d20210316 to make the vcf into smc format for each distinguished (5/pop). 
 Then we estimate using the 2.estimate.slurm and finally we plot to pdf but also use the -c parameter to output the csv that pyrho needs
