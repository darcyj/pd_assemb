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
	prep_check_inputs <- function(otutable, metadata, md_order_col, tree=NULL, distmat=NULL, 
		md_names_col=NULL, 	rev_md_order=FALSE, merge_tree_0dists=TRUE, merge_distmat_0dists=TRUE, 
		drop_unobserved_otus=TRUE, rarefy=TRUE, rare_depth=0, rare_seed=12345){
		require(ape)

		# handle metadata sample names if they are supplied as a column instead of as rownames
		if(! is.null(md_names_col)){ 
			if(substring(metadata[,md_names_col], 1, 1) %in% as.character(0:9)){
				warning("WARNING: Some sampleIDs in metadata names column start with numbers. They will be discarded.")
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

		# do rarefaction, which may drop samps, so fix otutable and metadata again.
		if(rarefy == TRUE){
			if(rare_depth==0){
				rare_depth <- min(colSums(otutable))
			}
			message(paste("Rarefying otutable to a depth of", rare_depth, "observations".))
			nsamps_bef <- ncol(otutable)
			otutable <- rarefy_otu_table(otutable, depth=rare_depth, seed=rare_seed)
			nsamps_aft <- ncol(otutable)
			# check if any samples were dropped and give warning
			if(nsamps_aft < nsamps_bef){
				warning(paste("WARNING:", nsamps_bef - nsamps_aft, "samples dropped due to rarefaction."))
			}
			metadata <- metadata[ which_order_a2b(rownames(metadata), colnames(otutable)), ]
		}

		# order metadata and otutable by md_order_col
		otutable <- otutable[ , order(metadata[,md_order_col])]
		metadata <- metadata[ order(metadata[,md_order_col]), ]

		# check that everything worked above:
		if( ! all(rownames(metadata) == colnames(otutable))){
			stop("ERROR: otutable and metadata sorting failure.")
		}

		# figure out which otus appear in both otutable and tree/distmat, then subset each to only contain union
		message("Removing samples not present in both tree/distmat and otutable")
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

		return(list(metadata=metadata, cum_otutable=cum_otutable, otutable=otutable, tree=tree, distmat=distmat))
	}

## subset_inputs
	# subsets list of inputs 
	# inputs_list - list, from prep_check_inputs()
	# sub_tf - logical vector, same length as nrow(list$metadata)
	# ... - extra arguments passed to prep_check_inputs
	# note: cum_otutable input will be ignored and re-calculated from otutable
	subset_inputs <- function(inputs_list, sub_tf, rarefy=FALSE, ...){
		# check that data are sorted first
		if( ! all(colnames(inputs_list$otutable) == rownames(inputs_list$metadata)) ){
			stop("ERROR: otutable colnames and metadata rownames do not match exactly.")
		}else if( ! length(sub_tf) == ncol(inputs_list$otutable)){
			stop("ERROR: sub_tf length does not match otutable columns.")
		}

		# subset otutable, so everything else will get subset when prep_check_inputs is run below
		inputs_list$otutable <- inputs_list$otutable[,sub_tf]

		return(prep_check_inputs(
			otutable=inputs_list$otutable,
			metadata=inputs_list$metadata,
			tree=inputs_list$tree,
			distmat=inputs_list$distmat, 
			rarefy=rarefy, ...
		))
	}


## workflow script running pd accumulation model for one time-series data set.
	# uses inputs object generated with prep_check_inputs and subset for just one
	# time-series (e.g. one location) using subset_inputs.
	# out_dir is a directory where figures and data will be stored.
	# other arguments should be self-explanatory
	pd_assemb_workflow <- function(inputs_i, out_dir, n_dispersions=2000, n_neutral_sim=2000, ncores=2){
		# remove terminal / on out_dir if there is one
		while(endsWith(out_dir, "/")){ out_dir <- substr(out_dir, 1, nchar(out_dir) - 1) }
		dir.create(out_dir)

		# Generate empirical PD accumulation
		empirical_PDs <- parallel_pd(inputs_i$cum_otutable, inputs_i$tree, ncores=ncores)
		empirical_PDs_scaled <- scale_PDs(empirical_PDs)
		
		# Generate PDs for parameter estimation
		dispersions <- seq(from=-5, to=5, length.out=n_dispersions)
		end_step <- ceiling(nrow(inputs_i$metadata)/2) + 1
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
		est_i <- fit_phylo_assemb(
			pd_sim=dispersion_PDs_scaled,
			pd_emp=empirical_PDs_scaled,
			dispersions=dispersions,
			times=inputs_i$metadata$days_since_exp_start,
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
			times=inputs_i$metadata$days_since_exp_start,
			dispersions=neutral_dispersions
		)
		dev.off()

		model_objects <- list(
			empirical=list(
				empirical_PDs=empirical_PDs, 
				empirical_PDs_scaled=empirical_PDs_scaled
			), model_fitting=list(
				dispersions=dispersions,
				end_step=end_step,
				dispersion_PDs=dispersion_PDs,
				dispersion_PDs_scaled=dispersion_PDs_scaled,
				estimate=est_i
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


