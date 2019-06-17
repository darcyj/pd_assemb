#!/usr/bin/env Rscript
# run with:
# nohup ./analysis.r &

## set up R environment
	library(ape)
	library(parallel)
	library(data.table)
	source("../dev/assembly_model_functions_v9.r")
	global_ncores <- 8
	global_n_disp <- 1000
	global_n_neut <- 1000

## create folder for model outputs
	dir.create("model_outputs")

## read in and format inputs
	inputs <- prep_check_inputs(
		otutable=read.table("model_inputs/zOTUtable_filt_1000.txt", skip=1, sep='\t', header=T, row.names=1, comment.char=""),
		metadata=read.table("model_inputs/metadata.txt", stringsAsFactors=F, sep='\t', header=T, row.names=1, comment.char=""),
		md_order_col=5,
		tree=read.tree("model_inputs/zotu_seeds_tree.tre"),
		distmat=NULL,
		merge_tree_0dists=TRUE,
		rarefy=FALSE,
		name="Infant"
	)

## remove mother sample
	inputs <- subset_inputs(inputs, sub_tf=inputs$metadata$sex == "male", name="Infant")

## generate empirical PD curve for whole data set
	empirical_PDs <- parallel_pd(inputs$cum_otutable, inputs$tree, ncores=8)
	empirical_PDs_scaled <- scale_PDs(empirical_PDs)
	pdf("model_outputs/PD_accumulation_whole_dataset.pdf", useDingbats=FALSE)
		plot(empirical_PDs_scaled ~ inputs$metadata$age, type="l", xlab="Age (days)", ylab="Proportion of new phylodiversity observed")
		abline(v=161, col="red", lty=2)
		segments(x0=inputs$metadata$age, x1=inputs$metadata$age, y0=0, y1=-0.03, col="blue")
		abline(h=0)
	dev.off()

## run model for days 0-161 (before formula)
	inputs_0_161 <- subset_inputs(inputs, sub_tf=(inputs$metadata$age < 161), name="Before formula")
	# run workflow
	results_0_161 <- pd_assemb_workflow(
		inputs_0_161, out_dir="model_outputs/results_0_161", 
		ncores=2, 
		n_dispersions=global_n_disp, 
		n_neutral_sim=global_n_neut
	)

## run model for days 161-300
	inputs_161_300 <- subset_inputs(inputs, sub_tf=(inputs$metadata$age >= 161 & inputs$metadata$age <= 300), name="After formula")
	# run workflow
	results_161_300 <- pd_assemb_workflow(
		inputs_161_300, out_dir="model_outputs/results_161_300", 
		ncores=2, 
		n_dispersions=global_n_disp, 
		n_neutral_sim=global_n_neut
	)


## make plot
	est_list <- list(results_0_161, results_161_300)

	pdf("model_outputs/dispersion_comparison_violinplot.pdf", useDingbats=FALSE)
	plot_dispersion_ests_violin(est_list)
	dev.off()

## 5. Save workspace
	save(file="koenig_infant_analysis.rdata", list=ls())

