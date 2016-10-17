library(dplyr)

### Access data.

# Set path.
name_path <- "path_leading_to_your_input_files/"

# Read sample data files.
axd_sample_data <- read.table(paste0(as.character(name_path), "axd_sample_data.txt", sep=""), header=TRUE)
axd_genotypes <- subset(axd_sample_data, select = (c(2:length(axd_sample_data))))
axd_locus <- subset(axd_sample_data, select = "sample_id")

### Creat allele objects, based on the genotypes in "axd_genotypes".
axd_samples_allele_count1 <- data.frame(ifelse(axd_genotypes == 0, "2",
                                                 ifelse(axd_genotypes == 1,"1",
                                                        ifelse(axd_genotypes == 2,"0","-9"))))
axd_samples_allele_count1$allele <- as.character(1)
axd_samples_allele_count1$individual <- as.character(1:nrow(axd_samples_allele_count1))
axd_samples_allele_count1$sample_id <- as.character(axd_locus$sample_id)
axd_samples_allele_count1$location <- axd_locus$location
axd_samples_allele_count1 <- axd_samples_allele_count1[,c((length(axd_samples_allele_count1)-3):length(axd_samples_allele_count1),1:length(axd_genotypes))]

axd_samples_allele_count2 <- data.frame(ifelse(axd_genotypes == 0, "0",
                                                 ifelse(axd_genotypes == 1,"1",
                                                        ifelse(axd_genotypes == 2,"2","-9"))))
axd_samples_allele_count2$allele <- as.character(2)
axd_samples_allele_count2$individual <- as.character(1:nrow(axd_samples_allele_count2))
axd_samples_allele_count2$sample_id <- as.character(axd_locus$sample_id)
axd_samples_allele_count2$location <- axd_locus$location
axd_samples_allele_count2 <- axd_samples_allele_count2[,c((length(axd_samples_allele_count2)-3):length(axd_samples_allele_count2),1:length(axd_genotypes))]

# Reduce the data frames above to vectors of alleles, with a "marker" column.
axd_locustally <- ncol(axd_samples_allele_count1)-4

axd_samples_allele_only_count1 <- subset(axd_samples_allele_count1, select=c(5:ncol(axd_samples_allele_count1)))
axd_samples_allele_only_count1 <- data.frame(t(axd_samples_allele_only_count1))
rownames(axd_samples_allele_only_count1) <- 1:nrow(axd_samples_allele_only_count1)
colnames(axd_samples_allele_only_count1) <- as.character(1:ncol(axd_samples_allele_only_count1))
axd_samples_allele_only_count1$marker <- 1:axd_locustally
axd_samples_allele_only_count1$marker <- as.numeric(axd_samples_allele_only_count1$marker)
axd_samples_allele_only_count1 <- axd_samples_allele_only_count1[,c(ncol(axd_samples_allele_only_count1),1:(ncol(axd_samples_allele_only_count1)-1))]

axd_samples_allele_only_count2 <- subset(axd_samples_allele_count2, select=c(5:ncol(axd_samples_allele_count2)))
axd_samples_allele_only_count2 <- data.frame(t(axd_samples_allele_only_count2))
rownames(axd_samples_allele_only_count2) <- 1:nrow(axd_samples_allele_only_count2)
colnames(axd_samples_allele_only_count2) <- as.character(1:ncol(axd_samples_allele_only_count2))
axd_samples_allele_only_count2$marker <- 1:axd_locustally
axd_samples_allele_only_count2$marker <- as.numeric(axd_samples_allele_only_count2$marker)
axd_samples_allele_only_count2 <- axd_samples_allele_only_count2[,c(ncol(axd_samples_allele_only_count2),1:(ncol(axd_samples_allele_only_count2)-1))]

### Create genotypes in BGC format.

# Isolate each individual's genotype, assign it a name "individual_i", and
# combine the data for all individuals.
indivtally <- nrow(axd_genotypes)

for(i in 1:indivtally) {
        allele_1 <- subset(axd_samples_allele_only_count1, select=c("marker", i))
        allele_2 <- subset(axd_samples_allele_only_count2, select=c("marker", i))
        assign(paste0("individual_axd_", i), left_join(allele_1, allele_2, "marker"))
}

admixlist_raw <- individual_axd_1

for(i in 2:indivtally){
        admixlist_raw <- rbind(admixlist_raw, setNames(get(paste0("individual_axd_",i)), names(admixlist_raw)))
}
names(admixlist_raw)[names(admixlist_raw) == "1.x"] <- "x"
names(admixlist_raw)[names(admixlist_raw) == "1.y"] <- "y"

# Sort the output.
admixlist_sorted <- admixlist_raw[order(admixlist_raw$marker),]
rownames(admixlist_sorted) <- 1:nrow(admixlist_sorted)
admixlist_sorted[,1] <- as.character(admixlist_sorted[,1])
admixlist_sorted[,2] <- as.character(admixlist_sorted[,2])
admixlist_sorted[,3] <- as.character(admixlist_sorted[,3])

### Add population and locus names.

# Population names.
newpopnames <- matrix(c("pop 0","pop 0","pop 0"),nrow=1,ncol=3) # new pop names to insert

axd_locustally <- length(colnames(axd_genotypes))
popnames_spot <- as.integer(seq(1,axd_locustally*indivtally-indivtally+1,indivtally)) # indices of the new rows

insertrow_pop <- function(admixlist_sorted, newpopnames) {
        new_popnames_spot <- sort(popnames_spot) + seq(0, length(popnames_spot) - 1) # adjust for rows shifting by one for each previous insertion
        old_rows_admixlist_sorted <- seq(nrow(admixlist_sorted) + length(new_popnames_spot))[-new_popnames_spot] # get indices for the old rows
        admixlist_sorted[old_rows_admixlist_sorted,] <- admixlist_sorted # assign old rows
        admixlist_sorted[new_popnames_spot,] <- newpopnames # assign new rows
        admixlist_sorted
}

admixlist_sorted_popnames <- insertrow_pop(admixlist_sorted, newpopnames)

# Locus names.
newlocusnames <- matrix(c("locus","locus","locus"),nrow=1,ncol=3) # new locus names to insert
locusnames_spot <- as.integer(seq(1,(axd_locustally)*(indivtally+1)-indivtally,indivtally+1)) # indices of the new rows

insertrow_loc <- function(admixlist_sorted_popnames, newlocusnames) {
        new_locusnames_spot <- sort(locusnames_spot) + seq(0, length(locusnames_spot) - 1) # adjust for rows shifting by one for each previous insertion
        old_rows_admixlist_sorted_popnames <- seq(nrow(admixlist_sorted_popnames) + length(new_locusnames_spot))[-new_locusnames_spot] # get indices for old rows
        admixlist_sorted_popnames[old_rows_admixlist_sorted_popnames,] <- admixlist_sorted_popnames # assign old rows
        admixlist_sorted_popnames[new_locusnames_spot,] <- newlocusnames # assign new rows
        admixlist_sorted_popnames
}

admixlist_sorted_popnames_locusnames <- insertrow_loc(admixlist_sorted_popnames, newlocusnames)

### Clean up "admixlist_sorted_popnames_locusnames".

# Remove "locus" and "pop 0" entries where they are not needed.
admixlist <- admixlist_sorted_popnames_locusnames
admixlist$marker[admixlist$y == "locus"] <- ""
admixlist$y[admixlist$y == "locus"] <- ""
admixlist$y[admixlist$y == "pop 0"] <- ""

# Change "locus" cell values to sequentially numbered values.
admixlist <- as.data.frame(admixlist %>% group_by(x) %>% mutate(count = sequence(n())))
admixlist$count <- paste0("locus ", admixlist$count, sep="")
admixlist$x[admixlist$x == "locus"] <- admixlist$count[admixlist$x == "locus"]
admixlist$count <- NULL
admixlist$marker <- NULL

# Write to file.
write.table(admixlist, paste(name_path, "axd_in.txt", sep=""), col.names=FALSE, row.names=FALSE, quote = FALSE, sep = " ")

