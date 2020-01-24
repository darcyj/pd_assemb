#!/usr/bin/env Rscript
# run with:
# nohup ./analysis.r &

## set up R environment
	library(ape)
	library(parallel)
	library(data.table)
	library(ggplot2)
	source("../../R_functions/assembly_model_functions.r")

## global variables
	global_ncores <- 16
	global_n_disp <- 2000
	global_n_neut <- 2000

## read in and format inputs
	inputs <- prep_check_inputs(
		otutable=read.table("../model_inputs/esv_table.txt", sep='\t', header=T, row.names=1, comment.char=""),
		metadata=read.table("../model_inputs/metadata.txt", stringsAsFactors=F, sep='\t', header=T, row.names=1, comment.char=""),
		md_order_col=4, md_names_col=1,
		tree=read.tree("../model_inputs/esv_seqs_filt_aln_pfiltered.fasta.treefile"),
		distmat=NULL, rarefy=TRUE, rare_depth=5000, 
		merge_tree_0dists=TRUE, drop_unobserved_otus=TRUE
	)



## run workflow for each subject
	n_subj <- length(unique(inputs$metadata$subject))
	outputs_list <- vector("list", length = n_subj)
	dir.create("model_outputs")
	for(i in 1:n_subj){
		set.seed(123456)

		subj <- unique(inputs$metadata$subject)[i]
		message("#####")
		message(paste("# Processing subject", i, ":", subj))
		message("#####")
		inputs_subj <- subset_inputs(
			inputs_list=inputs,
			sub_tf=inputs$metadata$subject == subj,
			merge_tree_0dists=TRUE,
			drop_unobserved_otus=TRUE
		)

		out_fp <- paste("model_outputs/", subj, sep="")
		dir.create(out_fp)
		outputs_list[[i]] <- pd_assemb_workflow(
			inputs_subj, out_dir=out_fp, 
			ncores=global_ncores, 
			n_dispersions=global_n_disp, 
			n_neutral_sim=global_n_neut
		)

		# garbage collect
		temp <- gc()
	}

## make data frame of estimates, p-values, and treatment info

	results <- data.frame(
		subj = sapply(X=outputs_list, FUN=function(x){ sub(pattern=".*/", replacement="", x=x$info$out_dir) }),
		d = sapply(X=outputs_list, FUN=function(x){ x$model_fitting$estimate$est[1] }),
		d_lci = sapply(X=outputs_list, FUN=function(x){ x$model_fitting$estimate$est[2] }),
		d_hci = sapply(X=outputs_list, FUN=function(x){ x$model_fitting$estimate$est[3] }),
		n_times = sapply(X=outputs_list, FUN=function(x){ x$info$n_samples }),
		end_time = sapply(X=outputs_list, FUN=function(x){ x$neutral_testing$end_step }),
		n_otus = sapply(X=outputs_list, FUN=function(x){ x$info$n_otus }),
		total_pd = sapply(X=outputs_list, FUN=function(x){ x$info$total_tree_pd }),
		pval = sapply(X=outputs_list, FUN=function(x){ x$neutral_testing$prob }),
		stringsAsFactors=F
	)

## add subject metadata to results
	which_order_a2b <- function(unsorted, template){ as.integer(sapply(X=template, FUN=function(x){which(x==unsorted)})) }
	subject_metadata <- read.table("../model_inputs/subject_data.txt", sep='\t', header=T, stringsAsFactors=F)
	subject_metadata <- subject_metadata[which_order_a2b(subject_metadata$subject, results$subj), ]
	# double check:
	all(subject_metadata$subject == subject_metadata$subj)
	# combine (ggplot needs everything together....)
	subject_metadata <- data.frame(subject_metadata, results)
	subject_metadata <- subject_metadata[subject_metadata$treatment_group %in% c("Antibiotics", "Control"),]


## make plots
	# antibiotics
	pdf("antibiotics.pdf")
	ggplot(subject_metadata, aes(x=treatment_group, y=d, fill=treatment_group)) + 
	geom_dotplot(binaxis="y", stackdir="center") + 
	labs(title="",x="", y = "D estimate")+ 
	theme(legend.position="none")
	dev.off()

	# which are not significant?
	subject_metadata[which(subject_metadata$pval > 0.05), ]

	# test if mean Ds are different between antibiotics vs control
	# actually a Mann-Whitney U test
	wilcox.test(
		x=subject_metadata$d[subject_metadata$treatment_group == "Antibiotics"],
		y=subject_metadata$d[subject_metadata$treatment_group == "Control"],
	)
	# P = 0.3615
	



	# data crap
	pdf("d_vs_other_vars.pdf")
	par(mfrow=c(2,2))
	plot(d ~ n_otus, ylab="D estimate", xlab="Number of zOTUs", data=subject_metadata, pch=20)
	plot(d ~ total_pd, ylab="D estimate", xlab="Total phylodiversity", data=subject_metadata, pch=20)
	plot(d ~ n_times, ylab="D estimate", xlab="Number of time points", data=subject_metadata, pch=20)

	dev.off()


