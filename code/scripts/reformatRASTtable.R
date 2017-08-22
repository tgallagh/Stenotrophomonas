### Rob and Joann's script to reformat RAST categories table

setwd("~/GoogleDrive/Stenotrophomonas/output/tables/")
file <- read.delim('seedannotation.tsv')
file <- as.data.frame(file)
data <- file

data$Features <- as.character(data$Features)
indices <- lapply(1:nrow(data), function(x) {
	splitter <- strsplit(data[x,'Features'], 'fig')[[1]]
	splitter <- splitter[-1]
	unsplit <- paste('fig', splitter, sep = '')
	category <- rep(data$Category[x], length(splitter))
	subcategory <- rep(data$Subcategory[x], length(splitter))
	subsystem <- rep(data$Subsystem[x], length(splitter))
	role <- rep(data$Role[x], length(splitter))
	features <- unsplit
	df1 <- data.frame(Category = category, Subcategory = subcategory, Subsystem = subsystem, Role = role, Features = unsplit, stringsAsFactors = FALSE) 
	return(df1)} )

data1 <- do.call(rbind.data.frame, indices)

data1$Features <- unlist(lapply(data1$Features, function(x) {
	splitter <- strsplit(x, ',')[[1]][1]
	return(splitter) } ) )

write.table(data1, 'reformmatedseedannotations.txt', sep = '\t', quote = F, col.names = T, row.names = F)
