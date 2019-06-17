#!/usr/bin/env Rscript
# run with:
# nohup ./analysis.r &

## set up R environment
	library(ape)
	library(parallel)
	library(data.table)
	source("../../R_functions/assembly_model_functions.r")

## global variables
	global_ncores <- 16
	global_n_disp <- 500
	global_n_neut <- 500

## read in and format inputs
	inputs <- prep_check_inputs(
		otutable=read.table("../model_inputs/esv_table.txt", sep='\t', header=T, row.names=1, comment.char=""),
		metadata=read.table("../model_inputs/movingpix_metadata.txt", stringsAsFactors=F, sep='\t', header=T, row.names=1, comment.char=""),
		md_order_col=11,
		tree=read.tree("../model_inputs/esv_seqs_filt_aln_pfiltered.fasta.treefile"),
		distmat=NULL, rarefy=TRUE, rare_depth=5000, 
		merge_tree_0dists=TRUE, drop_unobserved_otus=TRUE, name="movingpix"
	)

	# subset so all 8 subject/site combos are sampled over the same days
	inputs$metadata$subjsite <- paste(inputs$metadata$individual, "_", inputs$metadata$body_site, sep="")
	shared_days <- unique(inputs$metadata$days_since_exp_start)
	for(ss in unique(inputs$metadata$subjsite)){
		shared_days <- shared_days[shared_days %in% inputs$metadata$days_since_exp_start[inputs$metadata$subjsite==ss]]
	}
	# subset_inputs does NOT rarefy by default, which is good since data were already rarefied above.
	inputs <- subset_inputs(
		inputs_list=inputs,
		sub_tf=inputs$metadata$days_since_exp_start %in% shared_days,
		merge_tree_0dists=TRUE,
		drop_unobserved_otus=TRUE,
		name="movingpix_shareddays"
	)

	

## male right hand contiguous (all days 13-31 inclusive)
	inputs_M3_R_palm <- subset_inputs(
		inputs_list=inputs,
		sub_tf=inputs$metadata$days_since_exp_start %in% 13:31 & inputs$metadata$subjsite == "M3_R_palm",
		merge_tree_0dists=TRUE,
		drop_unobserved_otus=TRUE,
		name="M3_R_palm"
	)
	results_M3_R_palm <- pd_assemb_windowed_workflow(
		inputs=inputs_M3_R_palm,
		window_size=5,
		out_dir="auto",
		ncores=12
	)

## male left hand contiguous (all days 13-31 inclusive)
	inputs_M3_L_palm <- subset_inputs(
		inputs_list=inputs,
		sub_tf=inputs$metadata$days_since_exp_start %in% 13:31 & inputs$metadata$subjsite == "M3_L_palm",
		merge_tree_0dists=TRUE,
		drop_unobserved_otus=TRUE,
		name="M3_L_palm"
	)
	results_M3_L_palm <- pd_assemb_windowed_workflow(
		inputs=inputs_M3_L_palm,
		window_size=5,
		out_dir="auto",
		ncores=12
	)

## female right hand contiguous (all days 13-31 inclusive)
	inputs_F4_R_palm <- subset_inputs(
		inputs_list=inputs,
		sub_tf=inputs$metadata$days_since_exp_start %in% 13:31 & inputs$metadata$subjsite == "F4_R_palm",
		merge_tree_0dists=TRUE,
		drop_unobserved_otus=TRUE,
		name="F4_R_palm"
	)
	results_F4_R_palm <- pd_assemb_windowed_workflow(
		inputs=inputs_F4_R_palm,
		window_size=5,
		out_dir="auto",
		ncores=12
	)


## female left hand contiguous (all days 13-31 inclusive)
	inputs_F4_L_palm <- subset_inputs(
		inputs_list=inputs,
		sub_tf=inputs$metadata$days_since_exp_start %in% 13:31 & inputs$metadata$subjsite == "F4_L_palm",
		merge_tree_0dists=TRUE,
		drop_unobserved_otus=TRUE,
		name="F4_L_palm"
	)
	results_F4_L_palm <- pd_assemb_windowed_workflow(
		inputs=inputs_F4_L_palm,
		window_size=5,
		out_dir="auto",
		ncores=12
	)


## female feces contiguous (all days 13-31 inclusive)
	inputs_F4_feces <- subset_inputs(
		inputs_list=inputs,
		sub_tf=inputs$metadata$days_since_exp_start %in% 13:31 & inputs$metadata$subjsite == "F4_feces",
		merge_tree_0dists=TRUE,
		drop_unobserved_otus=TRUE,
		name="F4_feces"
	)
	results_F4_feces <- pd_assemb_windowed_workflow(
		inputs=inputs_F4_feces,
		window_size=5,
		out_dir="auto",
		ncores=12
	)
