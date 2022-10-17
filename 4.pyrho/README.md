# pyrho 

I am also working on restructuring the pyrho pipeline and making efficient scripts that handle the outputs. 
While the pipeline is still TO DO the make_windows process has been sped up immensely and can be applied to most inputs (especially with different window sizes). 
I need to explain it in more detail and make the script agnostic to input name (ex. I kinda need them named ch_m146_S-S-# for it to work).  

	old_scripts/ is a directory with ... well... old scripts 
	0.aux_scripts/ has the script needed to prepare the files before the pyrho pipeline (not an actual pipeline)
	x.pyrho_scripts/ has the scripts for the pyrho pipeline 
	1.results_RPunrel has results (and example results) for the dataset used after identifying unrelated individuals for the refpanel
