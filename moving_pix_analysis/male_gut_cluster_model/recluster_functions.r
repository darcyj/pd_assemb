# function to read in fasta or fastq file as list
# NOTE: ONLY WORKS ON SEQUENTIAL FASTX FILES. 
# Hope you aren't still using interleaved files after like 2005
read.fastx <- function(file, type=c("fastq", "fasta")){
	# read in raw data
	lines <- scan(file, what="character", sep='\n')
	
	# check to make sure data are OK-ish
	nlines <- length(lines)
	if(type=="fastq" & nlines %% 4 > 0){
		stop("CRITICAL ERROR: number of lines in fastq file not divisible by 4.")
	}else if(type == "fasta" & nlines %% 2 > 0){
		stop("CRITICAL ERROR: number of lines in fasta file not divisible by 2")
	}

	# make indices for item types (1=header, 2=seq, 3="+", 4=qual)
	if(type=="fastq"){
		line_inds <- rep(1:4, nlines/4)
	}else if(type=="fasta"){
		line_inds <- rep(1:2, nlines/2)
	}

	# make names generic
	headers <- as.vector(sapply(X=lines[line_inds == 1], FUN=function(x) substring(x, 2) ) )

	# make output list object
	if(type=="fastq"){
		output <- mapply(FUN=c, lines[line_inds == 2], lines[line_inds == 4], SIMPLIFY=FALSE)
	}else if(type=="fasta"){
		output <- mapply(FUN=c, lines[line_inds == 2], SIMPLIFY=FALSE)
	}

	# add names to output list
	names(output) <- headers

	# all done
	return(output)

}


# turns the above list into a simple string array of headers and seqs
# that can be written to text as a fasta file
fastalist2char <- function(fastalist){
	outvec <- rep("error", length(fastalist) * 2)
	for(i in 1:length(fastalist)){
		writepos <- (i * 2) - 1
		# write header
		outvec[writepos] <- paste(">", names(fastalist)[i], sep="")
		# write seq
		outvec[writepos+1] <- fastalist[[i]][1]
	}
	return(outvec)
}


# runs vsearch on sequences
vsearch_r <- function(seqs, id_cutoff, tmpfile_prefix="vsearch_tmp_seqs"){
	tmpfile <- paste(tmpfile_prefix, "_", Sys.getpid(), ".fasta", sep="")
	require(data.table)
	# write seqs to tmpfile
	# the "input" argument for system() is NOT working properly...
	fwrite(list(fastalist2char(seqs)), file=tmpfile)

	# make vsearch string
	v_string <- paste("vsearch --cluster_fast ", tmpfile, " --id", id_cutoff, "--uc -")

	# run vsearch, capture uc output!
	uctable <- system(command=v_string, intern=T, ignore.stderr=T)
	uctable <- read.table(text=uctable, sep='\t', header=F, stringsAsFactors=F)

	# remove temp file
	rmstring <- paste("rm", tmpfile)
	system(command=rmstring)

	return(uctable)
}


collapse_otutable_uctable <- function(otutable, uctable){
	# prune extraneous info from uctable, make otus self-referential
	uctable <- uctable[uctable[,1] != "S",9:10]
	uctable$V10[uctable$V10=="*"] <- uctable$V9[uctable$V10=="*"]
	# uctable and otutable should be SAME number of rows right now. check:
	if(nrow(uctable) != nrow(otutable)){stop("uctable and otutable not compatible!")}
	# sort uctable to match otutable
	which_order_a2b <- function(unsorted, template){ as.integer(sapply(X=template, FUN=function(x){which(x==unsorted)})) }
	uctable <- uctable[which_order_a2b(unsorted=uctable$V9, template=rownames(otutable)), ]
	# aggregate each column of otutable by uctable$V10
	agg_fun <- function(x){	agg_tab <- aggregate(x, by=list(uctable$V10), FUN=sum); out <- agg_tab[,2]; names(out) <- agg_tab[,1]; return(out)}
	return(apply(X=otutable, MAR=2, FUN=agg_fun))
}
