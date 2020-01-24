

## function for rarefying a sample
	# written for a different project, kind of slow.
	rarefy_sample <- function(abunds, names=NULL, depth){
		# this does rarefaction by indexing observations to species
		names_alias <- 1:length(abunds)
		obs_alias <- factor(unlist(apply(X=data.frame(names_alias, abunds), MAR=1, FUN=function(x){rep(x[1], x[2])} )), levels=names_alias)
		table_rare <- table( sample(obs_alias, size=depth, replace = FALSE) )
		# because obs_alias is a factor, table_rare will ALWAYS have colnames == names_alias
		return(as.vector(table_rare))
	}

	# re-written in c++ to go faster because why wait?
	# x: vector of species abunds
	# d: depth to which you'd like to rarefy
	# note that this will just take d down to sum(x) if d > x.
	library(Rcpp)
	cppFunction('IntegerVector rarefy_sample_cpp(const IntegerVector x, const int d) {
		// x is the "sample", i.e. column from otutable
		// need to expand x to be an array of size sum(x)
		int xsum = sum(x);
		IntegerVector xexp(xsum);
		int nspec = x.size();
		int pos = 0;
		int ind = 0;
		// for each position in x
		for(int i = 0; i < nspec; i++){
			int si = x[i];
			// for each observation at x[i]
			for(int j = 0; j < si; j++){
				xexp[pos] = ind;
				pos++;
			}
			ind++;
		}
		// shuffle xexp
		std::random_shuffle(xexp.begin(), xexp.end());
		// build output
		IntegerVector output(x.size());
		for(int i = 0; i < d; i++){
			// prevent errors by making sure d is reasonable!
			if(i < xsum){
				output[xexp[i]]++;
			}else{
				break;
			}
		}
		return output;
	}')


## indiv_null function, which creates an otutable that resamples
	# individual observations instead of species. 
	# x: an otu table (prep_check_inputs()$otutable)
	# tree: a tree (prep_check_inputs()$tree)
	# nsamp: sample at which to calculate PD
	indiv_null <- function(x, tree, nsamp){

		# this simply rarefies the totals in x to the total combined depth expected
		# at nsamp, and calculates pd. This is just a fast way of getting cumulative
		# pd at sample n. 
		x1 <- x[,1]
		totals <- rowSums(x) - x1
		depth2n <- sum(colSums(x)[2:nsamp])
		# totalnull <- rarefy_sample(totals, depth=depth2n) + x1
		totalnull <- rarefy_sample_cpp(x=totals, d=depth2n) + x1

		return(one_sample_pd(totalnull, rownames(x), tree))
	}

## This is just a much slower version of the function above, used to show
	# that they get the same answer (sanity check). Don't use it unless Rcpp doesn't
	# work for you, or you're REALLY curious.
	indiv_null_slow <- function(x, tree, nsamp){

		# this simply rarefies the totals in x to the total combined depth expected
		# at nsamp, and calculates pd. This is just a fast way of getting cumulative
		# pd at sample n. 
		output <- matrix
		x1 <- x[,1]
		totals <- rowSums(x) - x1
		output <- matrix(0, nrow=nrow(x), ncol=ncol(x))
		output[,1] <- x1
		for(i in 2:nsamp){
			di <- sum(x[,i])
			# oi <- rarefy_sample(totals, depth=di)
			oi <- rarefy_sample_cpp(x=totals, d=di)
			totals <- totals - oi
			output[,i] <- oi
		}
		output <- accumulate_otutable(output)

		return(one_sample_pd(output[,nsamp], rownames(x), tree))
	}

