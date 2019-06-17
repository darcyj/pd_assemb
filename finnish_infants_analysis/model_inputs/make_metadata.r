## create metadata using sample names
	otutable_temp <- read.table("esv_table.txt", sep='\t', header=T, row.names=1)
	samples <- colnames(otutable_temp)
	subjects <- sub(pattern="_.*", replacement="", x=samples)
	months <- sub(pattern=".*?_", replacement="", x=samples)
	months <- as.numeric(sub(pattern="_", replacement=".", x=months))
	
	metadata <- data.frame(
		SampleID=samples,
		subject=subjects,
		age_months=months,
		age_days=months * 30.5, 
		stringsAsFactors=F
	)

## create metadata for subjects
	# data that are per-subject, not per-sample
	subject_data <- read.table("../download/subject_metadata_guittar.txt", sep='\t', header=T, stringsAsFactors=F)

	# reconcile metadata and subject data
	subject_data <- subject_data[subject_data$subject %in% metadata$subject,]
	metadata <- metadata[metadata$subject %in% subject_data$subject, ]

## write out metadata
	write.table(metadata, file="metadata.txt", quote=F, sep='\t')

## write out subject data
	write.table(subject_data, file="subject_data.txt", quote=F, sep='\t')