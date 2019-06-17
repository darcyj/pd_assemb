# assembly_model_functions.r
	require(parallel)
	require(ape)

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

