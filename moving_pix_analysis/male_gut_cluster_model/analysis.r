#!/usr/bin/env Rscript

## NOTE - THIS ANALYSIS IS USING AN UN-RAREFIED OTU TABLE

## set up R environment
	library(ape)
	library(parallel)
	library(data.table)
	source("../../R_functions/assembly_model_functions.r")
	# functions for reading, writing and re-clustering sequences
	source("recluster_functions.r")

	# filepath for storing rdata files
	rdata_fp <- "../../moving_pix/rdata_backups/"

## create folder for model outputs
	dir.create("model_outputs")

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
	# keep only male gut
	inputs <- subset_inputs(inputs, inputs$metadata$body_site=="feces" & inputs$metadata$individual=="M3")

	# read in unaligned ESV rep seqs and get rid of unused seqs
	esv_seqs <- read.fastx(file="../model_inputs/esv_seqs_deint.fasta", type="fasta")
	esv_seqs <- esv_seqs[names(esv_seqs) %in% rownames(inputs$otutable)]

# vector of identity values to iterate over
	ids <- seq(from=1.000, to=0.94, by=-0.002)

# list to store results
	results <- list()

# run workflow for each identity value
	for(i in 1:length(ids)){
		
		# directory for i:
		dir_i <- paste("model_outputs/identity_", sprintf("%.3f", ids[i]), sep="")

		# re-cluster otutable and remake inputs, don't bother if ids[i] == 1
		if(ids[i] == 1){
			inputs_i <- inputs
		}else{
			uctable_i <- vsearch_r(esv_seqs, id_cutoff=ids[i])
			inputs_i <- prep_check_inputs(
				otutable=collapse_otutable_uctable(inputs$otutable, uctable_i),
				metadata=inputs$metadata,
				md_order_col=11,
				tree=inputs$tree,
				distmat=NULL
			)
		}

		# run workflow
		results[[i]] <- c(results, pd_assemb_workflow(
			inputs_i, out_dir=dir_i, 
			ncores=10, 
			n_dispersions=1000, 
			n_neutral_sim=100
		))

		# garbage collect
		temp <- gc()

		print(paste(ids[i], ";", i, "DONE"))
		
	}	

# extract Ds and their confidence intervals
	D <- sapply(X=results, FUN=function(x){x$model_fitting$estimate$est[1]})
	D_lci <- sapply(X=results, FUN=function(x){x$model_fitting$estimate$est[2]})
	D_hci <- sapply(X=results, FUN=function(x){x$model_fitting$estimate$est[3]})

	n_otus <- sapply(X=results, FUN=function(x){x$info$n_otus})

# make plot
	pdf("cluster_model_plots.pdf", useDingbats=F)

	plot(D ~ ids, ylim=c(-1, 0), pch=20, xlab="Re-clustering identity (vsearch)", ylab="D estimate")
	segments(x0=ids, y0=D, y1=D_hci)
	segments(x0=ids, y0=D, y1=D_lci)

	plot(n_otus ~ ids, ylab="Number of OTUs", xlab="Re-clustering identity (vsearch)", pch=20)

	dev.off()