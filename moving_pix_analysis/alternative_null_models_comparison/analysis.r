#!/usr/bin/env Rscript

## set up R environment
	library(ape)
	library(parallel)
	library(data.table)
	source("../../R_functions/assembly_model_functions.r")
	# functions for alternative null models (takes a second due to compilation of c++ code)
	source("alt_null_functions.r")

# set a seed
	set.seed(12345)

## read in and format inputs for all subjectxsites
	inputs <- prep_check_inputs(
		otutable=read.table("../model_inputs/esv_table.txt", sep='\t', header=T, row.names=1, comment.char=""),
		metadata=read.table("../model_inputs/movingpix_metadata.txt", stringsAsFactors=F, sep='\t', header=T, row.names=1, comment.char=""),
		md_order_col=11,
		tree=read.tree("../model_inputs/esv_seqs_filt_aln_pfiltered.fasta.treefile"),
		distmat=NULL, rarefy=TRUE, rare_depth=5000, 
		merge_tree_0dists=TRUE, drop_unobserved_otus=TRUE,
		name="all_movingpix_inputs"
	)

## wrapper function for null comparison (m = sample midpoint, in paper)
	null_wrapper <- function(bodysite, individual, m=54, nperm=500){
		# subset inputs	
		inputs_sub <- subset_inputs(inputs, inputs$metadata$body_site==bodysite & inputs$metadata$individual==individual)
		# run indiv null
		ind_nuls <- replicate(nperm, indiv_null(inputs_sub$otutable, inputs_sub$tree, nsamp=m))
		# get empirical pd
		emp <- one_sample_pd(inputs_sub$cum_otutable[,m], ids=rownames(inputs_sub$cum_otutable), tree=inputs_sub$tree)
		# calculate p-value
		pval <- (sum(ind_nuls < emp) + 1) / length(ind_nuls)
		# make plot
		plottitle <- paste0(individual, ", ", bodysite, ", P=", pval)
		xrange <- range(c(emp, ind_nuls)) + c(-1 * sd(ind_nuls), sd(ind_nuls))
		hist(ind_nuls, xlim=xrange, main=plottitle, col="red", xlab="PDm")
		abline(v=emp, col="blue", lwd=3)
	}

# for each subxsite combo, run null_wrapper
	pdf("ind_null_results.pdf")
	par(mfrow=c(4,2))
	for(bs in unique(inputs$metadata$body_site)){
		for(ind in unique(inputs$metadata$individual)){
			null_wrapper(bs, ind)
		}
	}
	dev.off()
