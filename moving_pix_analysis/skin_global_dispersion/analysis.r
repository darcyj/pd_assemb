#!/usr/bin/env Rscript

## NOTE - THIS ANALYSIS IS USING AN UN-RAREFIED OTU TABLE

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
		merge_tree_0dists=TRUE, drop_unobserved_otus=TRUE
	)

# subset to only include skin environment (palms)
	inputs_palms <- subset_inputs(
		inputs_list=inputs, 
		sub_tf=inputs$metadata$body_site=="R_palm" | inputs$metadata$body_site=="L_palm", 
		merge_tree_0dists=TRUE,
		drop_unobserved_otus=TRUE # only want to sample from species observed in the palms
	)

## male right hand both with and without "metacommunity"
	sub_tf <- inputs_palms$metadata$body_site=="R_palm" & inputs_palms$metadata$individual=="M3"
	MRH_meta <- subset_inputs( inputs_list=inputs_palms, sub_tf=sub_tf, drop_unobserved_otus=FALSE )
	MRH_self <- subset_inputs( inputs_list=inputs_palms, sub_tf=sub_tf, drop_unobserved_otus=TRUE )

	MRH_self_results <- pd_assemb_workflow(MRH_self, "MRH_self", ncores=global_ncores, n_dispersions=global_n_disp, n_neutral_sim=global_n_neut )
	MRH_meta_results <- pd_assemb_workflow(MRH_meta, "MRH_meta", ncores=global_ncores, n_dispersions=global_n_disp, n_neutral_sim=global_n_neut )

	rm(MRH_meta, MRH_self)
	gc()

## male left hand
	sub_tf <- inputs_palms$metadata$body_site=="L_palm" & inputs_palms$metadata$individual=="M3"
	MLH_meta <- subset_inputs( inputs_list=inputs_palms, sub_tf=sub_tf, drop_unobserved_otus=FALSE )
	MLH_self <- subset_inputs( inputs_list=inputs_palms, sub_tf=sub_tf, drop_unobserved_otus=TRUE )

	MLH_self_results <- pd_assemb_workflow(MLH_self, "MLH_self", ncores=global_ncores, n_dispersions=global_n_disp, n_neutral_sim=global_n_neut)
	MLH_meta_results <- pd_assemb_workflow(MLH_meta, "MLH_meta", ncores=global_ncores, n_dispersions=global_n_disp, n_neutral_sim=global_n_neut)

	rm(MLH_meta, MLH_self)
	gc()

## female right hand both with and without "metacommunity"
	sub_tf <- inputs_palms$metadata$body_site=="R_palm" & inputs_palms$metadata$individual=="F4"
	FRH_meta <- subset_inputs( inputs_list=inputs_palms, sub_tf=sub_tf, drop_unobserved_otus=FALSE )
	FRH_self <- subset_inputs( inputs_list=inputs_palms, sub_tf=sub_tf, drop_unobserved_otus=TRUE )

	FRH_self_results <- pd_assemb_workflow(FRH_self, "FRH_self", 
		ncores=global_ncores, n_dispersions=global_n_disp, 
		n_neutral_sim=global_n_neut)
	FRH_meta_results <- pd_assemb_workflow(FRH_meta, "FRH_meta", ncores=global_ncores, n_dispersions=global_n_disp, n_neutral_sim=global_n_neut)

	rm(FRH_meta, FRH_self)
	gc()

## female left hand
	sub_tf <- inputs_palms$metadata$body_site=="L_palm" & inputs_palms$metadata$individual=="F4"
	FLH_meta <- subset_inputs( inputs_list=inputs_palms, sub_tf=sub_tf, drop_unobserved_otus=FALSE )
	FLH_self <- subset_inputs( inputs_list=inputs_palms, sub_tf=sub_tf, drop_unobserved_otus=TRUE )

	FLH_self_results <- pd_assemb_workflow(FLH_self, "FLH_self", ncores=global_ncores, n_dispersions=global_n_disp, n_neutral_sim=global_n_neut)
	FLH_meta_results <- pd_assemb_workflow(FLH_meta, "FLH_meta", ncores=global_ncores, n_dispersions=global_n_disp, n_neutral_sim=global_n_neut)

	rm(FLH_meta, FLH_self)
	gc()

## make plot
	get_plot_data <- function(toget){
		attach(paste(toget, "/model_objects.rdata", sep=""), name=toget)
		if(startsWith(toget, "F")){sex <- "Female"}else{sex <- "Male"}
		if(endsWith(toget, "meta")){type <- "meta"}else{type <- "self"}
		output <- data.frame( id=toget, sex=sex, type=type,
			est=model_objects$model_fitting$estimate$est[1],
			est_lci=model_objects$model_fitting$estimate$est[2],
			est_hci=model_objects$model_fitting$estimate$est[3],
			pval=model_objects$neutral_testing$prob
		)
		detach()
		return(output)
	}

	folder_names <- c( "FLH_self", "FRH_self", "MLH_self", "MRH_self", "FLH_meta", "FRH_meta", "MLH_meta", "MRH_meta")
	folder_names <- folder_names[order(folder_names)]
	
	plotdata <- do.call("rbind", lapply(X=folder_names, FUN=get_plot_data))

	library(ggplot2)

	pdf("summary_plot.pdf")
	ggplot(plotdata, aes(x=id, y=est, color=sex, fill=type)) +
	scale_fill_manual(values=c("gray", "black")) +
	scale_color_manual(values=c("red", "blue")) +
	geom_bar(stat="identity", size=0.5) + 
	geom_errorbar(aes(ymin=est_lci, ymax=est_hci), width=0.2, position=position_dodge(.9))
	dev.off()
