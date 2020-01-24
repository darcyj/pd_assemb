#!/usr/bin/env Rscript

## set up R environment
	library(ape)
	library(parallel)
	library(data.table)
	source("../../R_functions/assembly_model_functions.r")
	# functions for alternative null models (takes a second due to compilation of c++ code)
	source("alt_null_functions.r")

## set a seed
	set.seed(12345)

## create folder for outputs
	dir.create("outputs")

## read in and format inputs for all subjectxsites
	# this was copied from the regular analysis
	inputs <- prep_check_inputs(
		otutable=read.table("../model_inputs/zOTUtable_filt_1000.txt", skip=1, sep='\t', header=T, row.names=1, comment.char=""),
		metadata=read.table("../model_inputs/metadata.txt", stringsAsFactors=F, sep='\t', header=T, row.names=1, comment.char=""),
		md_order_col=5,
		tree=read.tree("../model_inputs/zotu_seeds_tree.tre"),
		distmat=NULL,
		merge_tree_0dists=TRUE,
		rarefy=FALSE,
		name="Infant"
	)


## wrapper function for null comparison (m = sample midpoint, in paper)
	null_wrapper <- function(daymin, daymax, nperm=1000){
		subject <- paste0(daymin, "-", daymax)
		# subset inputs	
		inputs_sub <- subset_inputs(inputs, inputs$metadata$age >= daymin & inputs$metadata$age <= daymax )
		# calculate m
		m <- floor(ncol(inputs_sub$otutable)/2)
		# run indiv null
		ind_nuls <- replicate(nperm, indiv_null(inputs_sub$otutable, inputs_sub$tree, nsamp=m))
		# get empirical pd
		emp <- one_sample_pd(inputs_sub$cum_otutable[,m], ids=rownames(inputs_sub$cum_otutable), tree=inputs_sub$tree)
		# calculate p-value
		p_perm <- (sum(ind_nuls < emp) + 1) / length(ind_nuls)
		p_fit <- pnorm(q=emp, mean=mean(ind_nuls), sd=sd(ind_nuls))

		# make plot
		pdf(paste0("outputs/", subject, ".pdf"))
		plottitle <- paste0("days ", daymin, "-", daymax, ", P=", p_perm)
		xrange <- range(c(emp, ind_nuls)) + c(-1 * sd(ind_nuls), sd(ind_nuls))
		hist(ind_nuls, xlim=xrange, main=plottitle, col="red", xlab="PDm")
		abline(v=emp, col="blue", lwd=3)
		dev.off()
		return(data.frame(subject=subject, p_perm=p_perm, p_fit=p_fit))
	}

## for both day ranges, run null_wrapper
	
	null_wrapper(daymin=0, daymax=146)

	null_wrapper(daymin=161, daymax=297)



	# both are significantly underdispersed