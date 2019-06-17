## rarefy_otu_table
	# x: an data frame which contains otu counts. Columns are samples, rows are species.
	# will DROP columns where col sum is less than depth.
	rarefy_otu_table <- function(otu_table, depth, seed=12345){
		set.seed(seed)

		# function to rarefy one sample
		rarefy_sample <- function(abunds, depth){
			obs_2_samp <- rep.int(1:length(abunds), abunds)
			counts_rare <- table(sample(x=obs_2_samp, size=depth, replace=FALSE))
			output <- rep(0, length(abunds))
			output[as.integer(names(counts_rare))] <- counts_rare
			return(output)
		} 

		# remove columns of OTU table with insufficient sample depth
		otu_table <- otu_table[ , colSums(otu_table) >= depth]

		# rarefy remaining columns
		output_table <- apply(X=otu_table, MAR=2, FUN=rarefy_sample, depth=depth)
		dimnames(output_table) <- dimnames(otu_table)
		return(output_table)
	}



## prep_check_inputs
	# takes raw inputs for model and checks for completeness, sorts them, 
	# and formats them properly for downstream functions. 
	# otutable - matrix or df. rows=otus, cols=samples. REQUIRED.
		# colnames must be within metadata.
	# metadata - matrix or df. must contain rownames or a column containing all otutable rownames.
	# md_order_col - integer. which column of metadata contains the order? col be useable with 
		# R's order() function. 0=no re-ordering. REQUIRED.
	# tree - phylo object. tip names must be in otutable rownames. optional.
	# distmat - square matrix with row+colnames in otutable rownames. optional.
	# md_names_col - integer. which column of metadata contains names? use if no rownames. optional.
	# rev_md_order - logical. should order be reversed? default=FALSE.
	# merge_0dists - logical. should identical otus be merged? default=TRUE, but ignored if no tree/distmat provided.
	# rarefy - logical. Should otutable be rarefied? 
	# rare_depth - integer. To what depth should otutable be rarefied? 0 means automatic, which will rarefy to the lowest col sum. default=0
	# rare_seed - integer. RNG seed for rarefaction.
	# name - string. The name for this inputs object. If auto, object name of otutable will be used.
	prep_check_inputs <- function(otutable, metadata, md_order_col, tree=NULL, distmat=NULL, 
		md_names_col=NULL, rev_md_order=FALSE, merge_tree_0dists=TRUE, merge_distmat_0dists=TRUE, 
		drop_unobserved_otus=TRUE, rarefy=TRUE, rare_depth=0, rare_seed=12345, name="auto"){
		require(ape)

		# handle auto name
		if(name == "auto"){
			name <- deparse(substitute(otutable))
		}

		# handle metadata sample names if they are supplied as a column instead of as rownames
		if(! is.null(md_names_col)){ 
			if(any(substring(metadata[,md_names_col], 1, 1) %in% as.character(0:9))){
				stop("ERROR: Some sampleIDs in metadata names column start with numbers. That's not allowed.")
			}
			rownames(metadata) <- metadata[,md_names_col] 
		}

		# function for simultaneously subsetting and reordering
		which_order_a2b <- function(unsorted, template){ as.integer(sapply(X=template, FUN=function(x){which(x==unsorted)})) }

		# figure out which samples appear in both otutable and metadata, then subset both to only contain union
		message("Removing samples not present in both metadata and otutable")
		samps <- unique(c(rownames(metadata), colnames(otutable)))
		samps_union <- samps[ (samps %in% colnames(otutable)) & (samps %in% rownames(metadata))]
		metadata <- metadata[ which_order_a2b(rownames(metadata), samps_union), ]
		otutable <- otutable[ , which_order_a2b(colnames(otutable), samps_union) ]

		# order metadata and otutable by md_order_col
		otutable <- otutable[ , order(metadata[,md_order_col])]
		metadata <- metadata[ order(metadata[,md_order_col]), ]

		# check that everything worked above:
		if( ! all(rownames(metadata) == colnames(otutable))){
			stop("ERROR: otutable and metadata sorting failure.")
		}

		# figure out which otus appear in both otutable and tree/distmat, then subset each to only contain union
		message("Removing OTUs not present in both tree/distmat and otutable")
		otus <- unique(c(rownames(otutable), rownames(distmat), tree$tip.label))
		if(!is.null(distmat) && !is.null(tree)){
			# case where both tree and distmat are provided
			otus_union <- otus[(otus %in% rownames(otutable)) & (otus %in% rownames(distmat)) & (otus %in% tree$tip.label)]
		}else if(!is.null(tree)){
			# case where only tree is provided
			otus_union <- otus[(otus %in% rownames(otutable)) & (otus %in% tree$tip.label)]
		}else if(!is.null(distmat)){
			# case where only distmat is provided
			otus_union <- otus[(otus %in% rownames(otutable)) & (otus %in% rownames(distmat))]
		}else{
			# case where neither tree or distmat is provided
			otus_union <- otus
		}
		otutable <- otutable[which_order_a2b(rownames(otutable), otus_union), ]
		if(!is.null(distmat)){ distmat <- distmat[ which_order_a2b(rownames(distmat), otus_union), which_order_a2b(colnames(distmat), otus_union)]} 
		if(!is.null(tree)){ tree <- drop_extra(tree, keep=otus_union)}


		# do rarefaction, which may drop samps, so fix otutable and metadata again.
		if(rarefy == TRUE){
			if(rare_depth==0){
				rare_depth <- min(colSums(otutable))
			}
			message(paste("Rarefying otutable to a depth of", rare_depth, "observations"))
			nsamps_bef <- ncol(otutable)
			otutable <- rarefy_otu_table(otutable, depth=rare_depth, seed=rare_seed)
			nsamps_aft <- ncol(otutable)
			# check if any samples were dropped and give warning
			if(nsamps_aft < nsamps_bef){
				warning(paste("WARNING:", nsamps_bef - nsamps_aft, "samples were dropped due to rarefaction."))
			}
			metadata <- metadata[ which_order_a2b(rownames(metadata), colnames(otutable)), ]
		}


		# remove empty OTUs
		if(drop_unobserved_otus == TRUE){
			message("Removing empty OTUs")
			observed <- rowSums(otutable) > 0
			otutable <- otutable[observed, ]
			if(!is.null(distmat)){ distmat <- distmat[observed, observed] }
			if(!is.null(tree)){ tree <- drop_extra(tree, keep=rownames(otutable)) }
		}

		# merge tree 0-dists 
		if(merge_tree_0dists == TRUE && (!is.null(tree)) ){
			message("Merging 0-distance OTUs")
			otutable <- combine_0dist_otus(otutable=otutable, tree=tree, distmat=NULL)
			tree <- drop_extra(tree, keep=rownames(otutable))
			if(!is.null(distmat)){ distmat <- distmat[ 
				which_order_a2b(rownames(distmat), rownames(otutable)), 
				which_order_a2b(colnames(distmat), rownames(otutable))
			]}
		}

		# merge distmat 0-dists 
		if(merge_distmat_0dists == TRUE && (!is.null(distmat)) ){
			message("Reconciling distmat with merged OTUs")
			otutable <- combine_0dist_otus(otutable=otutable, tree=NULL, distmat=distmat)
			distmat <- distmat[ 
				which_order_a2b(rownames(distmat), rownames(otutable)), 
				which_order_a2b(colnames(distmat), rownames(otutable))
			]
			if(!is.null(tree)){ tree <- drop_extra(tree, keep=rownames(otutable)) }
		}

		# sort metadata and otutable by md_order_col
		message("Sorting objects by metadata order column")
		if(md_order_col > 0){ metadata <- metadata[ order(metadata[,md_order_col], decreasing=rev_md_order), ] }
		if(! all(colnames(otutable) == rownames(metadata))){
			otutable <- otutable[, which_order_a2b(colnames(otutable), rownames(metadata)  )]
		}

		# accumulate otu table
		message("Accumulating OTU table")
		cum_otutable <- accumulate_otutable(otutable, binary=TRUE)

		# double-check all sorting is propper
		message("Double-checking all objects match up")
		all_good <- TRUE
		if(! all( rownames(metadata) == colnames(otutable) )){all_good <- FALSE}
		if(! all( rownames(metadata) == colnames(cum_otutable) )){all_good <- FALSE}
		if( (!is.null(tree)) && (! all(rownames(otutable) %in% tree$tip.label)) ){all_good <- FALSE}
		if( !is.null(distmat) && (! all(rownames(otutable) == colnames(distmat))) ){all_good <- FALSE}
		if( !is.null(distmat) && (! all(rownames(otutable) == rownames(distmat))) ){all_good <- FALSE}
		if(all_good == FALSE){
			stop("ERROR: unknown catastrophic error.")
		}else{
			message("It's all good.")
		}

		return(list(metadata=metadata, cum_otutable=cum_otutable, otutable=otutable, tree=tree, distmat=distmat, md_order_col=md_order_col, name=name))
	}

## subset_inputs
	# subsets list of inputs 
	# inputs_list - list, from prep_check_inputs()
	# sub_tf - logical vector, same length as nrow(list$metadata)
	# name - if "auto", inputs_list$name will be used.
	# ... - extra arguments passed to prep_check_inputs
	# note: cum_otutable input will be ignored and re-calculated from otutable
	subset_inputs <- function(inputs_list, sub_tf, rarefy=FALSE, name="auto", ...){
		# check that data are sorted first
		if( ! all(colnames(inputs_list$otutable) == rownames(inputs_list$metadata)) ){
			stop("ERROR: otutable colnames and metadata rownames do not match exactly.")
		}else if( ! length(sub_tf) == ncol(inputs_list$otutable)){
			stop("ERROR: sub_tf length does not match otutable columns.")
		}

		if(name == "auto"){
			name <- inputs_list$name
		}

		# subset otutable, so everything else will get subset when prep_check_inputs is run below
		inputs_list$otutable <- inputs_list$otutable[,sub_tf]

		return(prep_check_inputs(
			otutable=inputs_list$otutable,
			metadata=inputs_list$metadata,
			tree=inputs_list$tree,
			distmat=inputs_list$distmat, 
			rarefy=rarefy, 
			md_order_col=inputs$md_order_col,
			name=name,
			...
		))
	}


## workflow script running pd accumulation model for one time-series data set.
	# uses inputs object generated with prep_check_inputs and subset for just one
	# time-series (e.g. one location) using subset_inputs.
	# very important that temporal sorting is done right in prep_check_inputs.
	# out_dir is a directory where figures and data will be stored. if "auto", 
	# the name parameter of inputs_i will be used (in pwd).
	# other arguments should be self-explanatory
	pd_assemb_workflow <- function(inputs_i, out_dir="auto", n_dispersions=2000, n_neutral_sim=2000, ncores=2, end_step=0){
		# remove terminal / on out_dir if there is one
		if(out_dir == "auto"){ out_dir <- inputs_i$name}
		while(endsWith(out_dir, "/")){ out_dir <- substr(out_dir, 1, nchar(out_dir) - 1) }
		dir.create(out_dir)

		# Generate empirical PD accumulation
		message("Calculating empirical PDs...")
		empirical_PDs <- parallel_pd(inputs_i$cum_otutable, inputs_i$tree, ncores=ncores)
		empirical_PDs_scaled <- scale_PDs(empirical_PDs)
		message("...done.")

		# Figure out end step (0 means auto, means half way)
		if(end_step == 0){
			end_step <- round(nrow(inputs_i$metadata)/2)
		}

		# Generate PDs for parameter estimation
		dispersions <- seq(from=-5, to=5, length.out=n_dispersions)
		dispersion_PDs <- multiple_phylo_assembly(
			cum_otutable=inputs_i$cum_otutable, 
			tree=inputs_i$tree, 
			dispersions=dispersions,
			ncores=ncores, noisy=T, make_distmat=T,
			endpoint_only=F, use_real_t1=T,
			end_step=end_step # half way point
		)
		dispersion_PDs_scaled <- scale_PDs(dispersion_PDs, ref=empirical_PDs)

		# Make bananagram plot and do model fitting
		pdf(paste(out_dir, "/", "model_fitting.pdf", sep=""), useDingbats=FALSE)
		fit_i <- fit_phylo_assemb(
			pd_sim=dispersion_PDs_scaled,
			pd_emp=empirical_PDs_scaled,
			dispersions=dispersions,
			times=inputs_i$metadata[,inputs_i$md_order_col],
			ncores=ncores
		)
		dev.off()

		# Generate PDs for significance testing
		neutral_dispersions <- rep(0, n_neutral_sim)
		neutral_PDs <- multiple_phylo_assembly(
			cum_otutable=inputs_i$cum_otutable, 
			tree=inputs_i$tree, 
			dispersions=neutral_dispersions,
			ncores=ncores, noisy=T, make_distmat=T,
			endpoint_only=F, use_real_t1=T,
			end_step=end_step # half way point
		)
		neutral_PDs_scaled <- scale_PDs(neutral_PDs, ref=empirical_PDs)

		# 2.5 test against neutral model (D=0)
		pdf(paste(out_dir, "/", "neutral_testing.pdf", sep=""), useDingbats=FALSE)
		p_i <- test_against_neutral(
			pd_sim=neutral_PDs_scaled,
			pd_emp=empirical_PDs_scaled,
			times=inputs_i$metadata[,inputs_i$md_order_col],
			dispersions=neutral_dispersions
		)
		dev.off()

		model_objects <- list(
			info=list(
				name=inputs_i$name,
				out_dir=out_dir,
				n_dispersions=n_dispersions,
				n_neutral_sim=n_neutral_sim,
				n_otus=nrow(inputs_i$cum_otutable),
				n_samples=ncol(inputs_i$cum_otutable),
				ncores=ncores,
				total_tree_pd=sum(inputs_i$tree$edge.length)
			), empirical=list(
				empirical_PDs=empirical_PDs, 
				empirical_PDs_scaled=empirical_PDs_scaled
			), model_fitting=list(
				dispersions=dispersions,
				end_step=end_step,
				dispersion_PDs=dispersion_PDs,
				dispersion_PDs_scaled=dispersion_PDs_scaled,
				estimate=fit_i$results,
				estimate_boot=fit_i$boot_ests
			), neutral_testing=list(
				dispersions=neutral_dispersions,
				end_step=end_step,
				neutral_PDs=neutral_PDs,
				neutral_PDs_scaled=neutral_PDs_scaled,
				prob=p_i
			)
		)

		# store stuff as rdata in out_dir
		save(list="model_objects", file=paste(out_dir, "/", "model_objects.rdata", sep=""))

		# return stuff
		return(model_objects)
	}	

# little function for banners
	banner <- function(x, char="█"){
		a <- paste(rep(char, nchar(x) + 6), collapse="")
		message(a)
		message(paste(char, char, " ", x, " ", char, char, sep=""))
		message(a)
	}

# pd_assemb_windowed_workflow
	# Does workflow using a sliding window approach. Each window is treated as its
	# own analysis, with evaluation occurring in the middle of the window. Thus, the
	# number of samples per window must be odd. 
	# inputs: an inputs list, as generated by prep_check_inputs() or subset_inputs()
		# must NOT have duplicate time points.
	# window_size: size of sliding window, in samples (i.e. 5 samples per window)
	# out_dir: place to store outputs. "auto" will use inputs$name.
	# ncores: number of CPU cores for parallel processing.
	# ndisp: number of dispersions for parameter estimation.
	# nneut: number of neutral replicates for testing against neutral model.
	pd_assemb_windowed_workflow <- function(inputs, window_size=5, out_dir="auto", ncores=2, ndisp=500, nneut=500){
		# auto out_dir name:
		if(out_dir == "auto"){
			out_dir <- paste(inputs$name, "/WS", window_size, sep="")
		}
		dir.create(out_dir, recursive=TRUE)

		# check that all time points are unique:
		if(max(table(inputs$metadata[, inputs$md_order_col])) > 1){
			stop("ERROR: Input metadata has duplicate times.")
		}

		# check that window size is odd:
		if(window_size %% 2 == 0){
			stop("ERROR: window_size must be odd.")
		}

		# sliding window stuff
		time_inds <- 1:nrow(inputs$metadata)
		start_inds <- 1:(length(time_inds) - window_size + 1)

		# list to hold results stuff
		results_list <-  vector("list", length=length(start_inds))

		# run model for each window
		for( w in start_inds){

			banner(paste("Time point ", w, "/", length(start_inds), sep=""))

			inds_w <- w:(w + window_size - 1)
			# subset inputs to only window
			inputs_w <- subset_inputs(
				inputs_list=inputs,
				sub_tf=(1:nrow(inputs$metadata)) %in% inds_w,
				merge_tree_0dists=TRUE,
				drop_unobserved_otus=TRUE,
				name=paste(out_dir, "/W", w, sep="")
			)
			# run workflow
			results_list[[w]] <- pd_assemb_workflow(
				inputs_w, out_dir="auto", n_dispersions=ndisp, n_neutral_sim=nneut, 
				ncores=ncores, end_step=ceiling(window_size / 2)
			)

		}

		# get stuff to plot
		offset <- (window_size - 1) / 2
		mid_inds <- start_inds + offset
		mid_times <- inputs$metadata[mid_inds, inputs$md_order_col]
		min_times <- inputs$metadata[mid_inds - offset, inputs$md_order_col]
		max_times <- inputs$metadata[mid_inds + offset, inputs$md_order_col]

		results_df <- data.frame(
			Dest=sapply(X=results_list, FUN=function(x){x$model_fitting$estimate$est[1]}),
			Dlci=sapply(X=results_list, FUN=function(x){x$model_fitting$estimate$est[2]}),
			Dhci=sapply(X=results_list, FUN=function(x){x$model_fitting$estimate$est[3]}),
			Pval=sapply(X=results_list, FUN=function(x){x$neutral_testing$prob}),
			TotalPD=sapply(X=results_list, FUN=function(x){x$info$total_tree_pd}),
			n_otus=sapply(X=results_list, FUN=function(x){x$info$n_otus}),
			mid_inds=mid_inds,
			mid_times=mid_times,
			min_times=min_times,
			max_times=max_times
		)
		results_df$significant=results_df$Pval <= 0.05

		# save plot data
		save(file=paste(out_dir, "/results_list.rdata", sep=""), list=c("results_df", "results_list"))

		# make plot
		pdf(paste(out_dir, "/dplot.pdf", sep=""), useDingbats=F)

		point_colors <- rep("black", nrow(results_df))
		point_colors[results_df$significant == FALSE] <- "white"


		plot(0, xlim=range(c(results_df$min_times, results_df$max_times)), ylim=range(c(results_df$Dlci, results_df$Dhci)), 
			xlab=colnames(inputs$metadata)[inputs$md_order_col], ylab="Estimated D")
		
		# error bars
		segments(x0=results_df$mid_times, y0=results_df$Dlci, y1=results_df$Dhci, x1=results_df$mid_times, col="gray")
		segments(x0=results_df$min_times, y0=results_df$Dest, y1=results_df$Dest, x1=results_df$max_times, col="gray")

		# points
		points(results_df$Dest ~ results_df$mid_times, col="black", bg=point_colors, pch=23, cex=2.5)

		dev.off()

		# return stuff
		return(results_df)

	}


## plot_dispersion_ests
	# plots dispersion estimates generated by fit_phylo_assemb() as vioilins
	# returns a ggplot object
	# est_list: list, where each item in the list is an output of pd_assemb_workflow.
	# fill_var: string vector, same length as est_list. variable to determine violin fill.
	# stroke_var: same as above but for stroke.
	# fill_cols: colors to use for fill. length of unique elements in fill_var.
	# stroke_cols: same as above but for stroke and with stroke_var.
	# For vioilin labels, name fields from each item in list will be used.
	plot_dispersion_ests_violin <- function(est_list, fill_var=NULL, stroke_var=NULL,
		fill_cols=NULL, stroke_cols=NULL){

		require(ggplot2)

		n <- length(est_list)

		# handle fill and stroke
		fill_leg <- stroke_leg <- TRUE
		if(is.null(fill_var)){ fill_leg <- FALSE; fill_var <- as.character(1:n); fill_cols <- rep("gray", n) }
		if(is.null(stroke_var)){ stroke_leg <- FALSE; stroke_var <- as.character(1:n); fill_cols <- rep("black", n) }

		if(length(fill_var) != n || length(stroke_var) != n){
			stop("Length of fill_var or stroke_var is bad. Must be 0 or length(est_list).")
		}


		# make df for each item in list.
		df <- vector("list", length=n)
		for(i in 1:n){
			nb <- nrow(est_list[[i]]$model_fitting$estimate_boot)
			df[[i]] <- data.frame(
				est_list[[i]]$model_fitting$estimate_boot,
				name=rep(est_list[[i]]$info$name, nb),
				fill_var=rep(fill_var[i], nb),
				stroke_var=rep(stroke_var[i], nb)
			)
		}
		df <- do.call("rbind", df)


		# get 95% CIs and means
		df2 <- data.frame(
			name=sapply(X=est_list, FUN=function(x){x$info$name}),
			mean=sapply(X=est_list, FUN=function(x){x$model_fitting$estimate$est[1]}),
			lci=sapply(X=est_list, FUN=function(x){x$model_fitting$estimate$est[2]}),
			hci=sapply(X=est_list, FUN=function(x){x$model_fitting$estimate$est[3]}),
			col=ifelse(stroke_leg, stroke_cols, rep("black", n))
		)


		if(fill_leg == TRUE && stroke_leg == TRUE){
			p <- ggplot(df, aes(x=name, y=est, fill=fill_var, color=stroke_var))
			p <- p + geom_violin(trim=F, scale="width")
			p <- p + suppressWarnings(geom_errorbar(aes(y=mean, ymin=lci, ymax=hci), data=df2, width = 0.2))
			p <- p + geom_point(aes(y=mean), data=df2, size=2)
			if(!is.null(stroke_cols)){ p <- p + scale_color_manual(values=stroke_cols)}
			if(!is.null(fill_cols)){ p <- p + scale_fill_manual(values=fill_cols)}
		}else if (fill_leg == TRUE && stroke_leg==FALSE){
			p <- ggplot(df, aes(x=name, y=est, fill=fill_var))
			p <- p + geom_violin(trim=F, scale="width", color="black")
			p <- p + geom_point(aes(y=mean), data=df2, size=2, color="black")
			p <- p + suppressWarnings(geom_errorbar(aes(y=mean, ymin=lci, ymax=hci), data=df2, color="black", width = 0.2))
			if(!is.null(fill_cols)){ p <- p + scale_fill_manual(values=fill_cols)}
		}else if (fill_leg == FALSE && stroke_leg==TRUE){
			p <- ggplot(df, aes(x=name, y=est, color=stroke_var))
			p <- p + geom_violin(trim=F, scale="width", fill="gray")
			p <- p + suppressWarnings(geom_errorbar(aes(y=mean, ymin=lci, ymax=hci), data=df2, width = 0.2))
			p <- p + geom_point(aes(y=mean), data=df2, size=2)
			if(!is.null(stroke_cols)){ p <- p + scale_color_manual(values=stroke_cols)}
		}else{
			p <- ggplot(df, aes(x=name, y=est))
			p <- p + geom_violin(trim=F, scale="width", fill="lightgray", color="black")
			p <- p + suppressWarnings(geom_errorbar(aes(y=mean, ymin=lci, ymax=hci), data=df2, color="black", width = 0.2))
			p <- p + geom_point(aes(y=mean), data=df2, size=2, color="black")
		}

		# classic theme best theme
		p <- p + theme_classic()
		# legend titles are stupid
		p <- p + theme(legend.title = element_blank()) 

		p

		return(p)

	}

## summarize summarize_workflow_results
	# summarizes output of pd_assemb_workflow()
	summarize_workflow_results <- function(x){
		data.frame(
			name=x$info$name,
			n_otus=x$info$n_otus,
			n_samp=x$info$n_samples,
			total_pd=x$info$total_tree_pd,
			estimate=x$model_fitting$estimate$est[1],
			est_lci=x$model_fitting$estimate$est[2],
			est_hci=x$model_fitting$estimate$est[3],
			p_value=x$neutral_testing$prob,
			stringsAsFactors=F
		)
	}

