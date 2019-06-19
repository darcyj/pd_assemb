#!/usr/bin/env Rscript
# run with:
# nohup ./analysis.r &

## set up R environment
	library(ape)
	library(parallel)
	library(data.table)
	source("../../R_functions/assembly_model_functions.r")

## global variables
	global_ncores <- 12
	global_n_disp <- 500
	global_n_neut <- 500

## read in and format inputs
	inputs <- prep_check_inputs(
		otutable=read.table("../model_inputs/esv_table.txt", sep='\t', header=T, row.names=1, comment.char=""),
		metadata=read.table("../model_inputs/movingpix_metadata.txt", stringsAsFactors=F, sep='\t', header=T, row.names=1, comment.char=""),
		md_order_col=11,
		tree=read.tree("../model_inputs/esv_seqs_filt_aln_pfiltered.fasta.treefile"),
		distmat=NULL, rarefy=TRUE, rare_depth=5000, 
		merge_tree_0dists=TRUE, drop_unobserved_otus=TRUE,
		name="all_movingpix_inputs"
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
		name="filtered_movingpix_inputs"
	)



## for each subject/site combo (ss), do the whole workflow
	
	ss_combos <- unique(inputs$metadata$subjsite)
	results_list <- vector("list", length=length(ss_combos))

	for(ssi in 1:length(ss_combos)){
		ss <- ss_combos[ssi]
		banner(paste("Analysis of", ss))

		# subset inputs to only contain samples from ss
		inputs_ss <- subset_inputs(
			inputs_list=inputs,
			sub_tf=inputs$metadata$subjsite == ss,
			merge_tree_0dists=TRUE,
			drop_unobserved_otus=TRUE,
			name=ss
		)

		# run workflow, saving in directory ss
		results_list[[ssi]] <- c(results_list, pd_assemb_workflow(
			inputs_ss, out_dir=ss, 
			ncores=global_ncores, 
			n_dispersions=global_n_disp, 
			n_neutral_sim=global_n_neut
		))

		# garbage collect
		temp <- gc()
	}


## save results as rdata
	save(list="results_list", file=paste(rdata_fp, "4x2_results.rdata", sep=""))
	# file was large so it was deleted so that analysis could be put on GitHub.
	# individual rdata files for each subjectXsite can be found in its respective folder

## re-order results_list
	results_list <- results_list[order(ss_combos)]

## look at P-values
	results_summary <- do.call("rbind", lapply(X=results_list, FUN=summarize_workflow_results))
	results_summary$sex <- ifelse(startsWith(results_summary$name, "F"), "Female", "Male")
	results_summary$site <- sapply(X=results_summary$name, FUN=function(x){tail(unlist(strsplit(x, split="_")), n=1)})
	results_summary$site <- tolower(results_summary$site)

## make plot
	pdf("results_violin.pdf", useDingbats=F)
	plot_dispersion_ests_violin(results_list, fill_var=results_summary$site, stroke_var=results_summary$sex,
		fill_cols=c("peru", "bisque", "lightcoral"), stroke_cols=c("darkred", "blue"))

	dev.off()