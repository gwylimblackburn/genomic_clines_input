library(dplyr)

### Access data.

# Set path.
name_path <- "path_leading_to_your_parental_input_files/"

# Read sample data files.
ptl0_sample_data <- read.table(paste0(name_path, "ptl0_sample_data.txt", sep=""), header=TRUE)
ptl0_sample_data$sample_id <- as.character(ptl0_sample_data$sample_id)

ptl1_sample_data <- read.table(paste0(name_path, "ptl1_sample_data.txt", sep=""), header=TRUE)
ptl1_sample_data$sample_id <- as.character(ptl1_sample_data$sample_id)

### Make a separate genotype data frame for each parental type.
ptl0_genotypes_pairwise <- ptl0_sample_data[,-1]
ptl1_genotypes_pairwise <- ptl1_sample_data[,-1]

### Turn individual rows of genotypes into two rows of alleles.

# Change genotype calls to allele counts for two separate allele objects. Genotypes
# not matching "0", "1", or "2" are assigned a value of "0" so that they don't alter
# the tally.
ptl0_allele_count1 <- data.frame(ifelse(ptl0_genotypes_pairwise == 0,2,
                                        ifelse(ptl0_genotypes_pairwise == 1,1,
                                               ifelse(ptl0_genotypes_pairwise == 2,0,0))))
ptl0_allele_count1$allele <- 1
ptl0_allele_count1$individual <- 1:nrow(ptl0_allele_count1)
ptl0_allele_count1$sample_id <- ptl1_sample_data$sample_id
ptl0_allele_count1 <- ptl0_allele_count1[,c((length(ptl0_allele_count1)-2):length(ptl0_allele_count1), 1:length(ptl0_genotypes_pairwise))]

ptl0_allele_count2 <- data.frame(ifelse(ptl0_genotypes_pairwise == 0,0,
                                        ifelse(ptl0_genotypes_pairwise == 1,1,
                                               ifelse(ptl0_genotypes_pairwise == 2,2,0))))
ptl0_allele_count2$allele <- 2
ptl0_allele_count2$individual <- 1:nrow(ptl0_allele_count2)
ptl0_allele_count2$sample_id <- ptl1_sample_data$sample_id
ptl0_allele_count2 <- ptl0_allele_count2[,c((length(ptl0_allele_count2)-2):length(ptl0_allele_count2), 1:length(ptl0_genotypes_pairwise))]

ptl1_allele_count1 <- data.frame(ifelse(ptl1_genotypes_pairwise == 0,2,
                                            ifelse(ptl1_genotypes_pairwise == 1,1,
                                                   ifelse(ptl1_genotypes_pairwise == 2,0,0))))
ptl1_allele_count1$allele <- 1
ptl1_allele_count1$individual <- 1:nrow(ptl1_allele_count1)
ptl1_allele_count1$sample_id <- ptl1_sample_data$sample_id
ptl1_allele_count1 <- ptl1_allele_count1[,c((length(ptl1_allele_count1)-2):length(ptl1_allele_count1), 1:length(ptl1_genotypes_pairwise))]


ptl1_allele_count2 <- data.frame(ifelse(ptl1_genotypes_pairwise == 0,0,
                                            ifelse(ptl1_genotypes_pairwise == 1,1,
                                                   ifelse(ptl1_genotypes_pairwise == 2,2,0))))
ptl1_allele_count2$allele <- 2
ptl1_allele_count2$individual <- 1:nrow(ptl1_allele_count2)
ptl1_allele_count2$sample_id <- ptl1_sample_data$sample_id
ptl1_allele_count2 <- ptl1_allele_count2[,c((length(ptl1_allele_count2)-2):length(ptl1_allele_count2), 1:length(ptl1_genotypes_pairwise))]

# Reduce the data frames above to vectors of alleles, each with a "marker" column.
ptl0_data_allele_only_count1 <- subset(ptl0_allele_count1, select=c(4:length(ptl0_allele_count1)))
ptl0_data_allele_only_count1 <- data.frame(t(ptl0_data_allele_only_count1))
rownames(ptl0_data_allele_only_count1) <- 1:nrow(ptl0_data_allele_only_count1)
colnames(ptl0_data_allele_only_count1) <- as.character(1:ncol(ptl0_data_allele_only_count1))
ptl0_data_allele_only_count1$marker <- 1:nrow(ptl0_data_allele_only_count1)
ptl0_data_allele_only_count1$marker <- as.numeric(ptl0_data_allele_only_count1$marker)
ptl0_data_allele_only_count1 <- ptl0_data_allele_only_count1[,c(length(ptl0_data_allele_only_count1),1:nrow(ptl0_sample_data))]

ptl0_data_allele_only_count2 <- subset(ptl0_allele_count2, select=c(4:length(ptl0_allele_count2)))
ptl0_data_allele_only_count2 <- data.frame(t(ptl0_data_allele_only_count2))
rownames(ptl0_data_allele_only_count2) <- 1:nrow(ptl0_data_allele_only_count2)
colnames(ptl0_data_allele_only_count2) <- as.character(1:ncol(ptl0_data_allele_only_count2))
ptl0_data_allele_only_count2$marker <- 1:nrow(ptl0_data_allele_only_count2)
ptl0_data_allele_only_count2$marker <- as.numeric(ptl0_data_allele_only_count2$marker)
ptl0_data_allele_only_count2 <- ptl0_data_allele_only_count2[,c(length(ptl0_data_allele_only_count2),1:nrow(ptl0_sample_data))]

ptl1_data_allele_only_count1 <- subset(ptl1_allele_count1, select=c(4:length(ptl1_allele_count1)))
ptl1_data_allele_only_count1 <- data.frame(t(ptl1_data_allele_only_count1))
rownames(ptl1_data_allele_only_count1) <- 1:nrow(ptl1_data_allele_only_count1)
colnames(ptl1_data_allele_only_count1) <- as.character(1:ncol(ptl1_data_allele_only_count1))
ptl1_data_allele_only_count1$marker <- 1:nrow(ptl1_data_allele_only_count1)
ptl1_data_allele_only_count1$marker <- as.numeric(ptl1_data_allele_only_count1$marker)
ptl1_data_allele_only_count1 <- ptl1_data_allele_only_count1[,c(length(ptl1_data_allele_only_count1),1:nrow(ptl1_sample_data))]

ptl1_data_allele_only_count2 <- subset(ptl1_allele_count2, select=c(4:length(ptl1_allele_count2)))
ptl1_data_allele_only_count2 <- data.frame(t(ptl1_data_allele_only_count2))
rownames(ptl1_data_allele_only_count2) <- 1:nrow(ptl1_data_allele_only_count2)
colnames(ptl1_data_allele_only_count2) <- as.character(1:ncol(ptl1_data_allele_only_count2))
ptl1_data_allele_only_count2$marker <- 1:nrow(ptl1_data_allele_only_count2)
ptl1_data_allele_only_count2$marker <- as.numeric(ptl1_data_allele_only_count2$marker)
ptl1_data_allele_only_count2 <- ptl1_data_allele_only_count2[,c(length(ptl1_data_allele_only_count2),1:nrow(ptl1_sample_data))]

### Create allele counts in BGC format.

# Isolate each individual's allele status, assign it a name "individual_i", and
# combine the data for all individuals.
for(i in 1:nrow(ptl0_sample_data)) {
        ptl0_allele_1 <- subset(ptl0_data_allele_only_count1, select=c("marker", i))
        ptl0_allele_2 <- subset(ptl0_data_allele_only_count2, select=c("marker", i))
        assign(paste0("ptl0_individual_", i), left_join(ptl0_allele_1, ptl0_allele_2, "marker"))
}

ptl0_parental_list_raw <- ptl0_individual_1

for(i in 2:(length(ptl0_data_allele_only_count1)-1)){
        ptl0_parental_list_raw <- rbind(ptl0_parental_list_raw, setNames(get(paste0("ptl0_individual_",i)), names(ptl0_parental_list_raw)))
}

names(ptl0_parental_list_raw)[names(ptl0_parental_list_raw) == "1.x"] <- "x"
names(ptl0_parental_list_raw)[names(ptl0_parental_list_raw) == "1.y"] <- "y"

for(i in 1:nrow(ptl1_sample_data)) {
        ptl1_allele_1 <- subset(ptl1_data_allele_only_count1, select=c("marker", i))
        ptl1_allele_2 <- subset(ptl1_data_allele_only_count2, select=c("marker", i))
        assign(paste0("ptl1_individual_", i), left_join(ptl1_allele_1, ptl1_allele_2, "marker"))
}

ptl1_parental_list_raw <- ptl1_individual_1

for(i in 2:(length(ptl1_data_allele_only_count1)-1)){
        ptl1_parental_list_raw <- rbind(ptl1_parental_list_raw, setNames(get(paste0("ptl1_individual_",i)), names(ptl1_parental_list_raw)))
}

names(ptl1_parental_list_raw)[names(ptl1_parental_list_raw) == "1.x"] <- "x"
names(ptl1_parental_list_raw)[names(ptl1_parental_list_raw) == "1.y"] <- "y"

# Sort the output.
ptl0_parental_list_sorted <- ptl0_parental_list_raw[order(ptl0_parental_list_raw$marker),]
rownames(ptl0_parental_list_sorted) <- 1:nrow(ptl0_parental_list_sorted)
ptl0_parental_list_sorted[,1] <- as.character(ptl0_parental_list_sorted[,1])
ptl0_parental_list_sorted[,2] <- as.character(ptl0_parental_list_sorted[,2])
ptl0_parental_list_sorted[,3] <- as.character(ptl0_parental_list_sorted[,3])

ptl1_parental_list_sorted <- ptl1_parental_list_raw[order(ptl1_parental_list_raw$marker),]
rownames(ptl1_parental_list_sorted) <- 1:nrow(ptl1_parental_list_sorted)
ptl1_parental_list_sorted[,1] <- as.character(ptl1_parental_list_sorted[,1])
ptl1_parental_list_sorted[,2] <- as.character(ptl1_parental_list_sorted[,2])
ptl1_parental_list_sorted[,3] <- as.character(ptl1_parental_list_sorted[,3])

# Sum allele counts per locus.
ptl0_parental_list_alter <- ptl0_parental_list_sorted
ptl0_parental_list_alter$marker <- as.numeric(ptl0_parental_list_alter$marker)
ptl0_parental_list_alter$x <- as.numeric(ptl0_parental_list_alter$x)
ptl0_parental_list_alter$y <- as.numeric(ptl0_parental_list_alter$y)
ptl0_parental_list_allelesums <- as.data.frame(ptl0_parental_list_alter %>% group_by(marker) %>% summarise_each(funs(sum)))

ptl1_parental_list_alter <- ptl1_parental_list_sorted
ptl1_parental_list_alter$marker <- as.numeric(ptl1_parental_list_alter$marker)
ptl1_parental_list_alter$x <- as.numeric(ptl1_parental_list_alter$x)
ptl1_parental_list_alter$y <- as.numeric(ptl1_parental_list_alter$y)
ptl1_parental_list_allelesums <- as.data.frame(ptl1_parental_list_alter %>% group_by(marker) %>% summarise_each(funs(sum)))

### Add locus names.
newlocusnames <- matrix(c("locus","locus","locus"),nrow=1,ncol=3) # new locus names to insert

ptl0_locustally <- length(ptl0_genotypes_pairwise)
ptl0_locusnames_spot <- as.integer(seq(1,ptl0_locustally,1)) # indices of the new rows

ptl0_insertrow_loc <- function(ptl0_parental_list_allelesums, newlocusnames) {
        ptl0_new_locusnames_spot <- sort(ptl0_locusnames_spot) + seq(0, length(ptl0_locusnames_spot) - 1) # adjust for rows shifting by one for each previous insertion
        ptl0_old_rows_parental_list_allelesums <- seq(nrow(ptl0_parental_list_allelesums) + length(ptl0_new_locusnames_spot))[-ptl0_new_locusnames_spot] # get indices for the old rows
        ptl0_parental_list_allelesums[ptl0_old_rows_parental_list_allelesums,] <- ptl0_parental_list_allelesums # assign old rows
        ptl0_parental_list_allelesums[ptl0_new_locusnames_spot,] <- newlocusnames # assign new rows
        ptl0_parental_list_allelesums
}

ptl0_parental_list_allelesums_locusnames <- ptl0_insertrow_loc(ptl0_parental_list_allelesums, newlocusnames)

ptl1_locustally <- length(ptl1_genotypes_pairwise)
ptl1_locusnames_spot <- as.integer(seq(1,ptl1_locustally,1)) # indices of the new rows

ptl1_insertrow_loc <- function(ptl1_parental_list_allelesums, newlocusnames) {
        ptl1_new_locusnames_spot <- sort(ptl1_locusnames_spot) + seq(0, length(ptl1_locusnames_spot) - 1) # adjust for rows shifting by one for each previous insertion
        ptl1_old_rows_parental_list_allelesums <- seq(nrow(ptl1_parental_list_allelesums) + length(ptl1_new_locusnames_spot))[-ptl1_new_locusnames_spot] # get indices for the old rows
        ptl1_parental_list_allelesums[ptl1_old_rows_parental_list_allelesums,] <- ptl1_parental_list_allelesums # assign old rows
        ptl1_parental_list_allelesums[ptl1_new_locusnames_spot,] <- newlocusnames # assign new rows
        ptl1_parental_list_allelesums
}

ptl1_parental_list_allelesums_locusnames <- ptl1_insertrow_loc(ptl1_parental_list_allelesums, newlocusnames)

### Clean up "ptl1_parental_list_allelesums_locusnames".

# Remove "locus" and "pop 0" entries where they are not needed.
ptl0_parental_list_allelesums_list <- subset(ptl0_parental_list_allelesums_locusnames, select=c(x,y))
ptl0_parental_list_allelesums_list$y[ptl0_parental_list_allelesums_list$y == "locus"] <- ""
ptl0_parental_list_allelesums_list$y[ptl0_parental_list_allelesums_list$y == "pop 0"] <- ""

ptl1_parental_list_allelesums_list <- subset(ptl1_parental_list_allelesums_locusnames, select=c(x,y))
ptl1_parental_list_allelesums_list$y[ptl1_parental_list_allelesums_list$y == "locus"] <- ""
ptl1_parental_list_allelesums_list$y[ptl1_parental_list_allelesums_list$y == "pop 0"] <- ""

# Change "locus" cell values to sequentially numbered values.
ptl0_parental_list_allelesums_list <- as.data.frame(ptl0_parental_list_allelesums_list %>% group_by(x) %>% mutate(count = sequence(n())))
ptl0_parental_list_allelesums_list$count <- paste0("locus ", ptl0_parental_list_allelesums_list$count, sep="")
ptl0_parental_list_allelesums_list$x[ptl0_parental_list_allelesums_list$x == "locus"] <- ptl0_parental_list_allelesums_list$count[ptl0_parental_list_allelesums_list$x == "locus"]
ptl0_parental_list_allelesums_list$count <- NULL

ptl1_parental_list_allelesums_list <- as.data.frame(ptl1_parental_list_allelesums_list %>% group_by(x) %>% mutate(count = sequence(n())))
ptl1_parental_list_allelesums_list$count <- paste0("locus ", ptl1_parental_list_allelesums_list$count, sep="")
ptl1_parental_list_allelesums_list$x[ptl1_parental_list_allelesums_list$x == "locus"] <- ptl1_parental_list_allelesums_list$count[ptl1_parental_list_allelesums_list$x == "locus"]
ptl1_parental_list_allelesums_list$count <- NULL

### Write to file.
write.table(ptl0_parental_list_allelesums_list, paste(name_path, "ptl0_in.txt", sep=""), col.names=FALSE, row.names=FALSE, quote = FALSE, sep = " ")
write.table(ptl1_parental_list_allelesums_list, paste(name_path, "ptl1_in.txt", sep=""), col.names=FALSE, row.names=FALSE, quote = FALSE, sep = " ")
