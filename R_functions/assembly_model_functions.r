# assembly_model_functions.r, Version 9.0
	# John L. Darcy, June 2019
	require(parallel)
	require(ape)
	require(ggplot2)

## drop_extra
	# drops extra tips from a tree
	# tree: phylo object, to be pruned
	# keep: character vector, tip names to retain
	# test code:
	# tree1 <- read.tree(text="(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);")
	# drop_extra(tree1, keep=c("A", "B", "D"))
	drop_extra <- function(tree, keep){
		treeabsent <- tree$tip.label[! tree$tip.label %in% keep]
		return(drop.tip(tree, treeabsent))
	}


## accumulate_otutable
	# takes a regular time-series OTU table and accumulates it 
	# otutable: matrix of observations, rows=otus, cols=samples
	# binary: T/F, should output be transformed to presence/absence? default is TRUE
	accumulate_otutable <- function(otutable, binary=TRUE){
		output <- t(apply(X=otutable, MAR=1, FUN=cumsum))
		rownames(output) <- rownames(otutable)
		colnames(output) <- colnames(otutable)
		if(binary){
			output <- (output > 0) * 1
		}
		return(output)
	}


## combine_0dist_otus
	# uses phylogenetic tree or distance matrix to combine OTUs that are phylogenetically identical
	# 0-distances don't work with the model because 0^-D is undefined. 
	# otutable: matrix or data.frame object with rownames. Does NOT need to be sorted to distmat.
	# tree: a phylo object, containing tips corresponding to rownames of otutable
	# distmat: a distance matrix, with rownames and colnames corresponding to rownames of otutable
	combine_0dist_otus <- function(otutable, tree=NULL, distmat=NULL, diagnostic=FALSE){
		require(ape)
		if( ! xor(is.null(tree), is.null(distmat)) ){stop("ERROR: must provide tree xor distmat.")}

		if( !is.null(tree) ){
			# make distmat if tree was provided
			if( ! all(rownames(otutable) %in% tree$tip.label) ){stop("ERROR: not all OTUs are in tree.")}
			tree <- drop_extra(tree, keep=rownames(otutable))
			distmat <- cophenetic.phylo(tree)
		}else{
			# check and simplify distmat if it was provided instead of tree
			if( ! all( rownames(distmat) == colnames(distmat))){ stop("ERROR: distmat isn't square and/or isn't sorted.") }
			if( ! all(rownames(otutable) %in% rownames(distmat)) ){stop("ERROR: not all OTUs are in distmat.")}
			distmat <- distmat[rownames(distmat) %in% rownames(otutable), colnames(distmat) %in% rownames(otutable)]
		}

		# make distmat boolean. TRUE means otus are identical. also make diag and upper tri all false
		distmat <- (distmat <= 0)
		diag(distmat) <- FALSE
		distmat[upper.tri(distmat)] <- FALSE

		# if there are no zeroes in distmat, just return otutable and be done.
		if( sum(distmat) == 0 ){
			return(otutable)
		}else{
			# assign categories to each otu
			categories <- 1:nrow(distmat)
			# for each cell in distmat i,j where distmat is TRUE
			for(i in 1:nrow(distmat)){
				# make sure otu i and otu j have the same category
				categories[distmat[i,]==TRUE] <- categories[i]
			}
			# aggregate
			which_order_a2b <- function(unsorted, template){ as.integer(sapply(X=template, FUN=function(x){which(x==unsorted)})) }
			otu_df <- data.frame(otu=rownames(distmat), cat=categories, newotu=rep("", length(categories)), stringsAsFactors=FALSE)
			otu_df <- otu_df[which_order_a2b(otu_df$otu, rownames(otutable)), ]
			if(!  all(rownames(otutable) == otu_df$otu)){stop("ERROR: critical sorting error.")}
			for(ci in unique(otu_df$cat)){
				otu_df$newotu[otu_df$cat==ci] <- otu_df$otu[otu_df$cat==ci][1]
			}

			if(diagnostic){
				for(o in unique(otu_df$newotu)){
					otus_in_o <- otu_df$otu[otu_df$newotu==o]
					if(length(otus_in_o) > 1){
						dists_o <- distmat[rownames(distmat) %in% otus_in_o, colnames(distmat) %in% otus_in_o]
						dists_o <- dists_o[lower.tri(dists_o)]
						print(paste(o, ";   ", paste(otus_in_o, collapse=","), sep=""))
						print(dists_o)
						print("  ")
					}
				}
			}

			otutable2 <- t(sapply(X=unique(otu_df$newotu), FUN=function(x){colSums(otutable[otu_df$newotu==x,])} ) )

			# some final checking, just in case
			if( nrow(otutable2) != length(unique(otu_df$newotu))){
				stop("ERROR: catastrophic aggregation error.")
			}

			return(otutable2)
		}
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
	prep_check_inputs <- function(otutable, metadata, md_order_col, tree=NULL, distmat=NULL, md_names_col=NULL, rev_md_order=FALSE, merge_0dists=TRUE){
		require(ape)
		# metadata rownames
		if(! md_order_col %in% 0:ncol(metadata)){ stop("ERROR: invalid md_order_col.") }
		if(! is.null(md_names_col)){ rownames(metadata) <- metadata[,md_names_col] }

		# check completeness
		if(! all(colnames(otutable) %in% rownames(metadata))){
			stop("ERROR: not all otutable colnames are in metadata rownames (or md_names_col).")
		}
		if( (!is.null(tree)) && (! all(rownames(otutable) %in% tree$tip.label)) ){
			stop("ERROR: not all otutable rownames are in tree tip labels.")
		}
		if( (!is.null(distmat)) && (! all(rownames(otutable) %in% rownames(distmat))) ){
			stop("ERROR: not all otutable rownames are in distmat rownames.")
		}
		if( (!is.null(distmat)) && (! all(rownames(otutable) %in% colnames(distmat))) ){
			stop("ERROR: not all otutable rownames are in distmat colnames")
		}
		if(!is.null(md_names_col) && ! md_names_col %in% 1:ncol(metadata)){
			stop("ERROR: invalid md_names_col") 
		}

		# drop empty rows from otutable
		otutable <- otutable[rowSums(otutable) > 0, ]

		# merge 0-dists from otutable maybe
		if(merge_0dists == TRUE){
			if( is.null(distmat) && is.null(tree) ){
				warning("WARNING: 0-dist OTUs NOT combined because no tree or distmat was provided.")
			}else if(!is.null(distmat)){
				otutable <- combine_0dist_otus(otutable=otutable, tree=NULL, distmat=distmat)
			}else{
				otutable <- combine_0dist_otus(otutable=otutable, tree=tree, distmat=NULL)
			}
		}

		# order objects to match
		if(md_order_col > 0){ metadata <- metadata[ order(metadata[,md_order_col], decreasing=rev_md_order), ] }
		which_order_a2b <- function(unsorted, template){ as.integer(sapply(X=template, FUN=function(x){which(x==unsorted)})) }
		metadata <- metadata[which_order_a2b(rownames(metadata), colnames(otutable)), ]
		if( !is.null(tree) ){ tree <- drop_extra(tree=tree, keep=rownames(otutable)) }
		if( !is.null(distmat)){ distmat <- distmat[which_order_a2b(rownames(distmat), rownames(otutable)), which_order_a2b(colnames(distmat), rownames(otutable))] }

		# accumulate otu table
		cum_otutable <- accumulate_otutable(otutable, binary=TRUE)

		# double-check all sorting is propper
		all_good <- TRUE
		if(! all( rownames(metadata) == colnames(otutable) )){all_good <- FALSE}
		if(! all( rownames(metadata) == colnames(cum_otutable) )){all_good <- FALSE}
		if( (!is.null(tree)) && (! all(rownames(otutable) %in% tree$tip.label)) ){all_good <- FALSE}
		if( !is.null(distmat) && (! all(rownames(otutable) == colnames(distmat))) ){all_good <- FALSE}
		if( !is.null(distmat) && (! all(rownames(otutable) == rownames(distmat))) ){all_good <- FALSE}
		if(all_good == FALSE){stop("ERROR: unknown catastrophic error.")}

		return(list(metadata=metadata, cum_otutable=cum_otutable, otutable=otutable, tree=tree, distmat=distmat))
	}


## subset_inputs
	# subsets list of inputs 
	# inputs_list - list, from prep_check_inputs()
	# sub_tf - logical vector, same length as nrow(list$metadata)
	# note: cum_otutable input will be ignored and re-calculated from otutable
	subset_inputs <- function(inputs_list, sub_tf){
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
			md_order_col=0,
			tree=inputs_list$tree,
			distmat=inputs_list$distmat
		))
	}


## phylo_assemb_model_step
	# basically, this function takes the overall model (next function) from i to i+1
	# OBS AND PDISTMAT MUST BE SORTED TO EACH OTHER!!!!
	# obs: boolean vector, which species have been observed so far
	# delta_sr: positive integer, how many new OTUs to sample for time i+1
	# pdistmat: square matrix, pair-wise phylogenetic distances between species. 
		# MUST HAVE SAME NROW AND NCOL AS LENGTH(OBS)
	# disp: real, the dispersion parameter (D) to use. D=0 means neutral model.
	phylo_assemb_model_step <- function(obs, delta_sr, pdistmat, disp){
		if(delta_sr <= 0){
			# if delta_sr is 0, just return obs. nothing changes.
			return(obs)
		}else{
			# get nyo or Not Yet Observed - boolean vector just like obs
			nyo <- !obs
			nyo_inds <- which(nyo)
			obs_inds <- which(obs)

			if(disp == 0){
				# disp is 0, so it's the neutral model. special case where 
				# we get to skip some computation.
				sampled_nyo_inds <- sample(nyo_inds, size=delta_sr)
			}else{
				# disp is NOT 0, so model must be used to calculate probabilities
				# calculate minimum distances from each nyo to obs
				min_ds <- apply(X=pdistmat[obs_inds, nyo_inds], MAR=2, FUN=min)
				# sample nyo_inds after raising min_ds to disp
				sampled_nyo_inds <- sample(nyo_inds, size=delta_sr, prob=min_ds^disp)
			}
			# add sampled nyo back to obs.

			obs[sampled_nyo_inds] <- TRUE
		}
		return(obs)
	}

## phylo_assemb_model
	# Simulates a cumulative OTU table using the phylo aseemb model
	# cum_otutable: matrix, a cumulative OTU table. see accumulate_otutable().
	# pdistmat: square matrix, pair-wise phylogenetic distances between species
	# disp: real, the dispersion parameter (D) to use. D=0 means neutral model.
	# end_step: integer, the last step to simulate. default (0) is all steps.
	# use_real_t1: bool, true means that t1 is real. false simulates t1.
	phylo_assemb_model <- function(cum_otutable, pdistmat, disp, end_step=0, use_real_t1=FALSE){

		# resolve end_step==0, which means to use all time points
		if(end_step == 0){ end_step <- ncol(cum_otutable) }

		# create empty OTU table to hold simulated data
		new_cum_otutable <- matrix(nrow=nrow(cum_otutable), ncol=end_step, data=0)
		rownames(new_cum_otutable) <- rownames(cum_otutable)
		colnames(new_cum_otutable) <- colnames(cum_otutable)[1:ncol(new_cum_otutable)]

		# calculate cumulative species richness (sr) of cum_otutable at each time point
		cum_srs <- colSums(cum_otutable > 0)

		# either create t1 or use real t1
		if(use_real_t1 == FALSE){
			new_cum_otutable[,1] <- sample(c(rep(1, cum_srs[1]), rep(0, nrow(cum_otutable) - cum_srs[1])), size=nrow(cum_otutable))
		}else{
			new_cum_otutable[,1] <- cum_otutable[,1]
		}

		for(j in 2:end_step){
			# calculate change in SR between empirical j-1 and j
			delta_sr_j <- as.vector(cum_srs[j] - cum_srs[j-1])

			# creating new column at j using j-1, keeping empirical delta sr.
			# using empirical delta sr is important because otherwise we can't tell the difference between sr and pd 
			new_cum_otutable[,j] <- phylo_assemb_model_step(
				obs=new_cum_otutable[,(j-1)] > 0,
				delta_sr=delta_sr_j,
				pdistmat=pdistmat,
				disp=disp
			)
		}

		return(new_cum_otutable)		
	}

# mclapply2 function
	# this is a version of mclapply that has a PROGRESS INDICATOR
	# THANKS SO MUCH to wannymahoots on stackoverflow, what a HERO.
	# this is over twice as fast as pblapply from the pbapply package.
	# https://stackoverflow.com/questions/10984556/is-there-way-to-track-progress-on-a-mclapply
	mclapply2 <- function(X, FUN, ..., mc.preschedule = TRUE, mc.set.seed = TRUE, mc.silent = FALSE, mc.cores = getOption("mc.cores", 2L), mc.cleanup = TRUE, mc.allow.recursive = TRUE, mc.progress=TRUE, mc.style=3) {
		if(version$major==3 && version$minor==5.1){
			message("Ignore \"cannot wait for child\" warnings. Only in R 3.5.1, they are not a problem.")
		}
		if (!is.vector(X) || is.object(X)) X <- as.list(X)

		if (mc.progress) {
			f <- fifo(tempfile(), open="w+b", blocking=T)
			p <- parallel:::mcfork()
			pb <- txtProgressBar(0, length(X), style=mc.style)
			setTxtProgressBar(pb, 0) 
			progress <- 0
			if (inherits(p, "masterProcess")) {
				while (progress < length(X)) {
					readBin(f, "double")
					progress <- progress + 1
					setTxtProgressBar(pb, progress) 
				}
				cat("\n")
				parallel:::mcexit()
			}
		}
		tryCatch({
			result <- mclapply(X, ..., function(...) {
					res <- FUN(...)
					if (mc.progress) writeBin(1, f)
					res
				}, 
				mc.preschedule = mc.preschedule, mc.set.seed = mc.set.seed,
				mc.silent = mc.silent, mc.cores = mc.cores,
				mc.cleanup = mc.cleanup, mc.allow.recursive = mc.allow.recursive
			)

		}, finally = {
			if (mc.progress) close(f)
		})
		result
	}

# one_sample_pd
	# this is a simplified version of the pd function from picante.
	# it can be parallelized across multiple samples (see parallel_pd), but
	# is less feature-rich than picante's implementation.
	# Because this is intended to go FAST, it assumes (and does not check):
		# tree tips and abund names are THE SAME. don't have to be in
		# the same order, but they have to have the same items within.
		# THAT'S VERY IMPORTANT!!!!!!!
		# This IS checked by the wrapper function parallel_pd.
	# returns a single PD value
	# abunds: vector of abundances (can be 0/1s)
	# ids: names for abunds, all must be in tree
	# tree: phylo object
	one_sample_pd <- function(abunds, ids, tree){
		treeabsent <- ids[abunds <= 0]
		SR <- length(ids) - length(treeabsent)
		if(SR == length(ids)){
			# whole tree pd, since nothing's missing
			PD <- sum(tree$edge.length)
		}else{
			subtree <- drop.tip(tree, treeabsent)
			PD <- sum(subtree$edge.length)
		}
		return(PD)
	}

# This is a parallelized wrapper of one_sample_pd()
	# it's also different because it takes an OTU table as input,
	# which is transposed versus the sample matrix picante uses. 
	# otutable should be a matrix object with both row and col labels.
		# rows = otus
		# cols = samples
	# returns a vector of PD values
	# otutable: matrix, an accumulated otu table (see accumulate_otutable)
	# tree: phylo object, must contain all rownames of otutable as tips
	# ncores: integer, number of cores to use for parallelization
	parallel_pd <- function(otutable, tree, ncores=2, quiet=F){
		require(parallel)
		# make sure tree and otutable have SAME OTUs inside
		otuids <- rownames(otutable)
		if(all(otuids %in% tree$tip.label) == FALSE){ # if there are OTUs in the otutable that AREN'T in the tree, error.
			stop("Some otus are not in the tree.")
		} else if(length(tree$tip.label) > length(otuids)){ # if tree has extra OTUs, get rid of them
			treeextras <- tree$tip.label[!tree$tip.label %in% otuids]
			tree <- drop.tip(tree, treeextras)
		}
		
		# make a list out of the OTU table columns. gotta have a list for mclcapply.
		otulist <- as.list(as.data.frame(otutable))
		
		# make sure we're not using more cores than we need
		if(ncores > length(otulist)){
			ncores <- length(otulist)
		}
		
		# parallel apply one_sample_pd
		if(ncores > 1 && quiet==FALSE){
			pd_list <- mclapply2(X=otulist, FUN=one_sample_pd, ids=rownames(otutable), tree=tree, mc.cores=ncores)
		}else if(ncores > 1 && quiet ==TRUE){
			pd_list <- mclapply(X=otulist, FUN=one_sample_pd, ids=rownames(otutable), tree=tree, mc.cores=ncores)
		}else{
			pd_list <- lapply(X=otulist, FUN=one_sample_pd, ids=rownames(otutable), tree=tree)
		}
		pds <- simplify2array(pd_list)

		return(pds)	
	}


## scale_PDs
	# scales PD values to 0:1 scale.
	# x: a vector or matrix of PDs (REQUIRED)
	# ref: a vector of PD values to use as a reference, i.e. to calculate minmax from.
	# using ref overrides minmax if both are present.
	scale_PDs <- function(x, ref=NULL){
		if(is.null(ref)){
			lo_pd <- min(x)
			delta <- max(x) - lo_pd
		}else{
			lo_pd <- min(ref)
			delta <- max(ref) - lo_pd
		}
		return( (x - lo_pd) / delta )	
	}


## multiple_phylo_assembly
	# runs phylo_assemb_model once for each value in dispersions
	# can be run in two modes: endpoint_only or all.
		# endpoint_only: only returns PD at endpoint (vector output)
		# all: returns PD for each timepoint (matrix output)
	# if all dispersions are 0 (neutral replicates), distmat can be NULL
	# otutable: matrix, cumulative otutable (see accumulate_otutable)
	# tree: phylo object, must have tips for each rowname of otutable
	# dispersions: vector, dispersions parameters for model runs
	# ncores: integer, number of cores to use for parallelization
	# noisy: logical, do you want a progress indicator? (default=T)
	# make_distmat: logical, if T distmat will be made from tree,
		# if F, distmat must be provided but tree can be NULL
		# default = TRUE
	# distmat: matrix, square distance matrix (default=NULL)
	# end_step: see phylo_assemb_model
	# use_real_t1: see phylo_assemb_model
	multiple_phylo_assembly <- function(cum_otutable, tree, dispersions, ncores=2, noisy=T, make_distmat=T, distmat=NULL, endpoint_only=T, end_step=0, use_real_t1=T){
		# make distance matrix if required
		if(make_distmat){
			if(noisy){message("Converting tree into phylogenetic distance matrix...")}
			distmat <- cophenetic.phylo(tree)
			if(noisy){message("   ...done.")}
		}

		# chceck and sort distmat to match cum_otutable
		if(all(rownames(distmat) == rownames(cum_otutable) && all(colnames(distmat) == rownames(cum_otutable)))){
			# everything is sorted nicely
			if(noisy){message("Otutable and distmat are properly sorted.")}
		}else{
			# everything is NOT sorted nicely (expected)
			if(noisy){message("Sorting cum_otutable and distmat to match...")}
			which_order_a2b <- function(unsorted, template){ as.integer(sapply(X=template, FUN=function(x){which(x==unsorted)})) }
			cum_otutable <- cum_otutable[which_order_a2b(unsorted=rownames(cum_otutable), template=rownames(distmat)), ]
			if(noisy){message("   ...done.")}
		}

		# run model and make otutables. this outputs a LIST of OTUTABLES. 
		if(noisy){message("Running model in parallel across dispersion values...")}
		pa_otutables_list_cum <- mclapply2(dispersions, FUN=function(x){phylo_assemb_model(cum_otutable=cum_otutable, pdistmat=distmat, disp=x, end_step=end_step, use_real_t1=use_real_t1)}, mc.cores=ncores)
		if(noisy){message("   ...done.")}

		# calculate PD (output)
		if(noisy){message("Calculating phylogenetic diversities...")}
		if(endpoint_only==TRUE){
			lastcols <- simplify2array(lapply(X=pa_otutables_list_cum, FUN=function(x){x[,end_step]}))

			# use parallel_pd on lastcols
			pds_out <- parallel_pd(lastcols, tree=tree, ncores=ncores)
			names(pds_out) <- dispersions
			if(noisy){message("   ...done.")}
			return(pds_out)

		}else{
			# for each otu table in pa_otutables_list_cum, get PDs
			pds_mat <- t(simplify2array(mclapply2(X=pa_otutables_list_cum, FUN=parallel_pd, tree=tree, ncores=1, mc.cores=ncores)))
			if(noisy){message("   ...done.")}
			# make that into a nice matrix and return it
			# rownames(pds_mat) <- dispersions
			return(pds_mat)
		}
	}

## fit_phylo_assemb
	# calculates error between simulated and empirical data, fits error model,
	# makes plots, returns list containing:
		# results: 
	# pd_sim: matrix of simulated pds, from multiple_phylo_assembly.
		# rownames can be dispersions, if not dispersions must be supplied separately.
		# can also be a column vector of pds from the comparative time point (required)
	# pd_emp: vector of empirical pds corresponding to rows of pd_sim (may be longer)
		# or a single pd value for the comparative time point (required)
	# comp_samp: integer, which sample index is the comparative time point?
		# default is to use final column/value of pd_sim.
	# times: vector, metadata for time, corresponding to  pd_emp (optional)
	# dispersions: vector of dispersions corresponding to pd_sim (optional)
	# make_plot: logical. make plots or not? (default=TRUE)
	# disp_col_range: range of dispersion parameters to color in plot (default=-1 to 1)
	# boot_prop: proportion of points to sample when bootstrapping (default=0.10)
	# boot_n: number of bootstraps to do (default=1000)
	# ncores: number of cpu cores to use for bootstrapping (default=2)
	fit_phylo_assemb <- function(pd_sim, pd_emp, comp_samp=0, times=NULL, dispersions=NULL, make_plot=TRUE, disp_col_range=c(-1,1), boot_prop=0.10, boot_n=1000, ncores=2){
		# determine what kind of data pd_sim are, and make universal inputs for fitting
		if(is.matrix(pd_sim)){
			# determine what comp_samp and dispersions should be if they were default
			if(comp_samp == 0){comp_samp <- ncol(pd_sim)}
			if(is.null(dispersions)){dispersions <- as.numeric(rownames(pd_sim))}
			# make inputs for fitting
			pd_sim_t <- as.vector(pd_sim[,comp_samp])
			if(length(pd_emp) != 1){
				pd_emp_t <- pd_emp[comp_samp]
			}else{
				pd_emp_t <- pd_emp
			}
		}else if(is.vector(pd_sim)){
			# determine what comp_samp and dispersions should be if they were default
			if(is.null(dispersions)){dispersions <- as.numeric(names(pd_sim))}
			if(length(pd_emp) > 1){
				if(comp_samp == 0){stop("ERROR: comp_samp must be specified for this input type.")}
				pd_emp_t <- pd_emp[comp_samp]
				pd_sim_t <- pd_sim
			}else{
				pd_emp_t <- pd_emp
				pd_sim_t <- pd_sim
			}
		}else{
			stop("ERROR: pd_sim must be a matrix or vector.")
		}

		# come up with colors for dispersions
		if(is.null(dispersions)){stop("ERROR: no dispersions provided.")}
		# 200 bins for dispersion colors (should be enough)
		disp_inc <- (max(disp_col_range) - min(disp_col_range)) / 200
		disp_bin_df <- data.frame(
			min=seq(from=min(disp_col_range), to=max(disp_col_range)+disp_inc, length=201)[1:200],
			max=seq(from=min(disp_col_range), to=max(disp_col_range)+disp_inc, length=201)[2:201],
			col=as.character(rev(rainbow(200, start=2/6, end=4/6))),

			stringsAsFactors=FALSE
		)

		# for each dispersion, figure out which bin its in and assign a color.
		disp_colors <- rep("gray", length(dispersions))
		for(i in 1:length(disp_colors)){
			x <- dispersions[i]
			if(x >= min(disp_col_range) && x < max(disp_col_range)){
				disp_colors[i] <- disp_bin_df$col[which(x >= disp_bin_df$min & x < disp_bin_df$max)]
			}
		}

		# make bananagram if pd_sim matrix was provided
		if(is.matrix(pd_sim) && make_plot==TRUE){
			if(length(pd_emp) != length(times)){
				warning("WARNING: can't make bananagram because pd_emp and times are not same length.")
			}else{
				# set R to do 2 plots
				par(mfrow=c(2,1))
				# make times for pd sim matrix (since it may be smaller)
				times_sim <- times[1:ncol(pd_sim)]
				# plot everything
				plot(pd_emp ~ times, ylim=range(c(pd_emp, pd_sim)), type="l", col="black", lwd=1, xlab="Time", ylab="PD accumulation")
				# plot lines in random order so that they don't "stack" as much
				set.seed(12345)
				# for each line, only plot if color isn't gray
				for(i in sample(1:nrow(pd_sim))){if(disp_colors[i] != "gray"){
					points(pd_sim[i,] ~ times_sim, col=disp_colors[i], type="l", lwd=1)
				}}
				# plot empirical again, thicker 
				points(pd_emp ~ times, type="l", col="black", lwd=4 )
				points(pd_emp ~ times, type="l", col="white", lwd=2, lty=2 )
				
				# calculate legend limits
				# far right of x-axis, 10% width, starts at 10% height, ends at 60% height
				leg_x0 <- min(times) + 0.70 * (max(times) - min(times)) # 90%
				leg_x1 <- leg_x0 + (max(times) - min(times)) * 0.10 # 80%
				leg_y1 <- 0.6 * max(pd_emp) - min(pd_emp) # 60%
				leg_y0 <- leg_y1 / 6 # 10%
				recths <- seq(from=0, to=(leg_y1 - leg_y0), length.out=nrow(disp_bin_df))

				# actually draw legend for real
				for(i in 1:nrow(disp_bin_df)){
					rect( xleft=leg_x0, xright=leg_x1, ybottom=leg_y0 + recths[i], ytop=leg_y1, col=disp_bin_df$col[i],	border=NA )
				}
				# border around legend
				rect(xleft=leg_x0, xright=leg_x1, ybottom=leg_y0, ytop=leg_y1, col=NULL, border="black", lwd=2)
				# legend labels
				labelys <- seq(from=leg_y0, to=leg_y1, length=5)
				tickwidth <- (max(times) - min(times)) * 0.02
				# ticks
				segments( x0=leg_x1, x1=leg_x1+tickwidth, y0=labelys, y1=labelys, col="black", lwd=2 )
				# labels
				text(x=leg_x1 + 1.5 * tickwidth, y=labelys, labels=format(round(as.vector(quantile(disp_bin_df$min)), 2), nsmall = 2), adj=0 )
				# x-axis tics for sampling points
				segments( x0=times, x1=times, y0=0, y1=-0.04, col="blue", lwd=2 )
			}
		}

		# deltas are how different simulated and empirical data are
		deltas <- (pd_emp_t - pd_sim_t)

		# fit logistic model. par <- c(r, i, ymin, K)
		logistic_fitter <- function(par, x=dispersions, y=deltas, ret="ssq"){
			r<-par[1] ; i<-par[2]; ymin<-par[3]; K<-par[4]
			pred <- ((K-ymin) / (1+exp( r*(x-i) ) ) ) +ymin
			# Y = (K / 1+exp(-r * x-i))+b
			if(ret=="pred"){ return(pred) }
			ssq <- sum((pred - y)^2)
			if(ret=="ssq"){ return(ssq)
			}else if(ret=="r2"){ return( 1 - (ssq / sum((data$y - mean(data$y))^2)))
			}else if(ret=="resid"){ return(data$y - pred)
			}
		}

		# this function used to solve logistic model for y=0 or y=whatever. Just use par=y, and par_known as parameters.
		logistic_reverse <- function(par, par_known){
			x <- par; r <- par_known[1]; i<-par_known[2]; ymin<-par_known[3]; K<-par_known[4]
			pred <- ((K-ymin) / (1+exp( r*(x-i) ) ) ) +ymin
			return(abs(0 - pred))
		}

		# this function wraps the two above, so bootstraps can be parallelized.
		logistic_boot <- function(seed, dispersions, deltas, boot_prop){
			set.seed(seed)
			par_init <- c(r=-3, i=0, ymin=min(deltas), K=max(deltas))
			boot_pos <- sample(1:length(dispersions), round(boot_prop * length(dispersions)), replace=FALSE)
			mod_boot <- optim(par_init, x=dispersions[boot_pos], y=deltas[boot_pos], logistic_fitter)
			# solve for x
			est <- optimize(f=logistic_reverse, interval=c(-2, 2), maximum=FALSE, par_known=mod_boot$par)$minimum
			return(c(mod_boot$par, est=est))
		}

		# do bootstrapping
		set.seed(123456)
		seeds <- replicate(boot_n, as.numeric(paste(sample(0:9, 6, replace=TRUE), collapse="")))
		message("Generating error model via logistic model bootstrapping...")
		boot_ests <- t(simplify2array(mclapply2(X=seeds, FUN=logistic_boot, dispersions=dispersions, deltas=deltas, boot_prop=boot_prop, mc.cores=ncores)))
		message("   ...done")
		# get results as mean, lci, and hci of each output
		get_summ <- function(x){
			x <- x[order(x)]
			lci025 <- x[ round(length(x) * 0.025) ]
			hci975 <- x[ round(length(x) * 0.975) ]
			return(c(mean=mean(x), lci025=lci025, hci975=hci975))
		}
		results <- as.data.frame(apply(X=boot_ests, MAR=2, FUN=get_summ))

		# generate mean logistic curve overlay
		fake_xs <- seq(from=min(dispersions), to=max(dispersions), length.out=1000)
		fake_par<- c( r=results$r[1], i=results$i[1], ymin=results$ymin[1], K=results$K[1] )
		fake_ys <- logistic_fitter(par=fake_par, x=fake_xs, ret="pred")

		# make model fitting plot
		plot(deltas ~ dispersions, col=disp_colors, pch=20, cex=1.5, xaxs="i", yaxs="i")
		points(fake_ys ~ fake_xs, type="l", lwd=3, col="black")
		arrowlen <- 0.015 * (max(dispersions) - min(dispersions))
		arrows(x0=min(dispersions), x1=results$est[1], y0=0, y1=0, col="red", lwd=3, length=arrowlen)
		arrows(x0=results$est[1], x1=results$est[1], y0=0, y1=min(deltas), col="red", lwd=3, length=arrowlen)

		# return(results)
		return(list(results=results, boot_ests=boot_ests))
	}

## test_against_neutral
	# Tests against neutral model (D=0). 
	# pd_sim: matrix of simulated pds where D=0, from multiple_phylo_assembly.
		# can also be a column vector of pds from the comparative time point (required)
	# pd_emp: vector of empirical pds corresponding to rows of pd_sim (may be longer)
		# or a single pd value for the comparative time point (required)
	# comp_samp: integer, which sample index is the comparative time point?
		# default is to use final column/value of pd_sim.
	# times: vector, metadata for time, corresponding to  pd_emp (optional)
	# make_plot: logical. make plots or not? (default=TRUE)
	test_against_neutral <- function(pd_sim, pd_emp, comp_samp=0, times=NULL, make_plot=TRUE, dispersions=NULL){
		# determine what kind of data pd_sim are, and make universal inputs for fitting
		if(is.matrix(pd_sim)){
			# determine what comp_samp and dispersions should be if they were default
			if(comp_samp == 0){comp_samp <- ncol(pd_sim)}
			if(is.null(dispersions)){dispersions <- as.numeric(rownames(pd_sim))}
			# make inputs for fitting
			pd_sim_t <- as.vector(pd_sim[,comp_samp])
			if(length(pd_emp) != 1){
				pd_emp_t <- pd_emp[comp_samp]
			}else{
				pd_emp_t <- pd_emp
			}
		}else if(is.vector(pd_sim)){
			# determine what comp_samp and dispersions should be if they were default
			if(is.null(dispersions)){dispersions <- as.numeric(names(pd_sim))}
			if(length(pd_emp) > 1){
				if(comp_samp == 0){stop("ERROR: comp_samp must be specified for this input type.")}
				pd_emp_t <- pd_emp[comp_samp]
				pd_sim_t <- pd_sim
			}else{
				pd_emp_t <- pd_emp
				pd_sim_t <- pd_sim
			}
		}else{
			stop("ERROR: pd_sim must be a matrix or vector.")
		}

		# make bananagram if pd_sim matrix was provided
		if(is.matrix(pd_sim) && length(pd_emp) > 1 && make_plot==TRUE){
			if(length(pd_emp) != length(times)){
				warning("WARNING: can't make bananagram because pd_emp and times are not same length.")
			}else{
				# set R to do 2 plots
				par(mfrow=c(2,1))
				# make times for pd sim matrix (since it may be smaller)
				times_sim <- times[1:ncol(pd_sim)]
				# plot everything
				plot(pd_emp ~ times, ylim=range(c(pd_emp, pd_sim)), type="l", col="black", lwd=1, xlab="Time", ylab="PD accumulation")
				set.seed(12345)
				# function to make a random red tone
				randomred <- function(){
					light <- sample(c(TRUE, FALSE), 1)
					change <- sample(0:255, 1)
					if(light){ return(rgb(255, change, change, maxColorValue=255))
					}else{return(rgb(change, 0, 0, maxColorValue=255))}
				}
				# plot each neutral replicate, use random shade of red
				for(i in sample(1:nrow(pd_sim))){
					points(pd_sim[i,] ~ times_sim, col=randomred(), type="l", lwd=1)
				}
				# plot empirical again, thicker 
				points(pd_emp ~ times, type="l", col="black", lwd=4 )
				points(pd_emp ~ times, type="l", col="white", lwd=2, lty=2 )
			}
		}

		# calculate P-value
		pval <- 1
		n_sim <- length(pd_sim_t)
		if(pd_emp_t >= mean(pd_sim_t)){
			pval <- sum(pd_emp_t <= pd_sim_t) / n_sim
		}else if(pd_emp_t < mean(pd_sim_t)){
			pval <- sum(pd_emp_t >= pd_sim_t) / n_sim
		}
		if(pval == 0){
			pval <- 1 / n_sim
		}

		# make histogram
		n_breaks <- round(n_sim * 0.05) # sturges' is too conservative
		ptitle <- paste("P =", round(pval, 3))
		hist(pd_sim_t, col="red", xlim=range(c(pd_emp_t, pd_sim_t)), main=ptitle, xlab="PD", breaks=n_breaks, lty="blank")
		abline(v=pd_emp_t, col="blue", lwd=3)

		return(pval)

	}

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

