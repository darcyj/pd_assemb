# does initial richness or PD predict D?


## set up R environment
	library(ape)
	library(parallel)
	library(data.table)
	library(ggplot2)
	source("../../R_functions/assembly_model_functions.r")

## read in and format inputs
	inputs <- prep_check_inputs(
		otutable=read.table("../model_inputs/esv_table.txt", sep='\t', header=T, row.names=1, comment.char=""),
		metadata=read.table("../model_inputs/metadata.txt", stringsAsFactors=F, sep='\t', header=T, row.names=1, comment.char=""),
		md_order_col=4, md_names_col=1,
		tree=read.tree("../model_inputs/esv_seqs_filt_aln_pfiltered.fasta.treefile"),
		distmat=NULL, rarefy=TRUE, rare_depth=5000, 
		merge_tree_0dists=TRUE, drop_unobserved_otus=TRUE
	)

## for each subject, get D and initial alpha diversity. store in data.frame.
	subjects <- system("ls model_outputs", intern=T)
	dval_hci <- dval_lci <- dval_est <- rep(-11, length(subjects))
	init_sobs <- init_pd <- rep(-11, length(subjects))
	for(i in 1:length(subjects)){
		s <- subjects[i]
		# get D
		s_fp <- paste("model_outputs", s, "model_objects.rdata", sep="/")
		load(s_fp)
		dval_est[i] <- model_objects$model_fitting$estimate$est[1]
		dval_lci[i] <- model_objects$model_fitting$estimate$est[2]
		dval_hci[i] <- model_objects$model_fitting$estimate$est[3]

		rm(model_objects)

		# figure out what initial richness was
		inputs_i <- subset_inputs(
			inputs_list=inputs,
			sub_tf=inputs$metadata$subject == s,
			merge_tree_0dists=TRUE,
			drop_unobserved_otus=TRUE, name=s
		)
		init_pd[i] <- one_sample_pd(abunds=inputs_i$otutable[,1], ids=rownames(inputs_i$otutable), tree=inputs_i$tree)
		init_sobs[i] <- sum(inputs_i$otutable[,1] > 0)
	}

## make output plots
	pdf("init_pd_vs_D.pdf")
	par(mfrow=c(2,2))
	plot(1,1)
	plot(1,1)
	plot(dval_est ~ init_pd, ylab="D estimate", xlab="Initial phylodiversity", pch=20)
	plot(dval_est ~ init_pd, ylab="D estimate", xlab="Initial number of zOTUs", pch=20)
	dev.off()

	# above plot and d_vs_other_vars.pdf combined into one figure using illustrator