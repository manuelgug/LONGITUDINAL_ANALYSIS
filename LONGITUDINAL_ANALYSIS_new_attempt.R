library(dplyr)
library(ggplot2)
library(purrr)

#####################################################################################3

####### IMPORTS AND FORMATTING #######

#import metadata
metadata <- readxl::read_xlsx("epi_sero_traj_data_filled_20240502.xlsx")

# subset useful columns
metadata <- metadata[c("Numero de estudo", "NIDA", "Visita", "traj_c3_58", "traj_c3_1w")]
metadata <- metadata[metadata$NIDA != "NA", ] #remove the one NA in NIDA
metadata <- metadata[metadata$Visita != "NPG", ] #remove NPGs
metadata <- metadata[metadata$traj_c3_58 != "Ex",] #remove excluded samples


# filtered allele data import and formatting
allele_data1 <- read.csv("BOH22_NextSeq04_RESULTS_v0.1.8_070324_FILTERED/allele_data_global_max_0_filtered.csv")
allele_data2 <- read.csv("ASINT_NextSeq01_RESULTS_FILTERED/allele_data_global_max_0_filtered.csv")
allele_data3 <- read.csv("ASINT_NextSeq02_RESULTS_v0.1.8_FILTERED_exclude_file/allele_data_global_max_0_filtered.csv")

# merge runs (WHEN ALL RUNS ARE READY)
allele_data <- rbind(allele_data1, allele_data2, allele_data3)


# reformat NIDA column
allele_data$sampleID <- sub("_S.*$", "", allele_data$sampleID) #remove _S* from sampleID
allele_data$sampleID <- sub("N", "", allele_data$sampleID)
allele_data$sampleID <- sub("_", ".", allele_data$sampleID)
allele_data$sampleID <- sub("\\.0", "", allele_data$sampleID)
colnames(allele_data)[1] <- "NIDA"

#keep pools 1A and 5 only.
allele_data <- allele_data %>%
  filter(grepl("-1A$|-1B$", locus)) # confirm which pools!!!!

# allele name formatting
allele_data$allele <- paste0(allele_data$locus, "___", allele_data$pseudo_cigar)

#merge paired_samples and allele_data by NIDA
merged_dfs <- merge(metadata, allele_data, by="NIDA", all.x = TRUE)


# NIDAS without sequencing
missing <- merged_dfs[is.na(merged_dfs$locus),]
missing_nidas <- unique(missing$NIDA)
print(missing_nidas)

#remove missing NIDAs from data
merged_dfs <- merged_dfs[!is.na(merged_dfs$locus),]

#####################################################################################3

####### POOL V1 AND V5 VISITS INTO V1V5 #######
merged_dfs_locifil_allefil <- merged_dfs %>%
  mutate(Visita = case_when(
    Visita %in% c("V1", "V5") ~ "V1V5",
    TRUE ~ Visita
  ))

merged_dfs_locifil_allefil <- merged_dfs_locifil_allefil %>%
  group_by(`Numero de estudo`, Visita, locus, Category, allele, pseudo_cigar) %>%
  summarize(reads = sum(reads))

# recalculate in-sample allele freqs
k<- merged_dfs_locifil_allefil %>% 
  group_by(`Numero de estudo`,Visita,locus,Category) %>%
  summarize(norm.reads.locus = reads / sum(reads))

merged_dfs_locifil_allefil <- cbind(k, merged_dfs_locifil_allefil[c("allele", "pseudo_cigar", "reads")])

#pseudo nida
merged_dfs_locifil_allefil$NIDA <- paste0(merged_dfs_locifil_allefil$`Numero de estudo`, "_", merged_dfs_locifil_allefil$Visita)


#####################################################################################
# #check for presence of ama-1 amplicon
# ama1<- merged_dfs_locifil_allefil[merged_dfs_locifil_allefil$locus == "Pf3D7_11_v3-1294284-1294580-1B",] #está en pocas muestras
# length(unique(ama1$NIDA))
# 
# median(ama1$reads)
# hist(ama1$reads)
#####################################################################################3


#####################################################################################3

####### 1) ALLELE FILTERING #######
MAF <- 0.01
min_reads <- 100

merged_dfs_locifil_allefil <- merged_dfs_locifil_allefil[merged_dfs_locifil_allefil$norm.reads.locus >= MAF,] #MAF
merged_dfs_locifil_allefil <- merged_dfs_locifil_allefil[merged_dfs_locifil_allefil$reads >= min_reads,] # 100 reads minimum
merged_dfs_locifil_allefil <- merged_dfs_locifil_allefil[!grepl("I=", merged_dfs_locifil_allefil$allele),] #remove alleles with I (insertion)
merged_dfs_locifil_allefil <- merged_dfs_locifil_allefil[!grepl("D=", merged_dfs_locifil_allefil$allele),] #remove alleles with D (deletion)


#####################################################################################
# #check for presence of ama-1 amplicon
# ama1<- merged_dfs_locifil_allefil[merged_dfs_locifil_allefil$locus == "Pf3D7_11_v3-1294284-1294580-1B",] #está en pocas muestras
# length(unique(ama1$NIDA))
# 
# median(ama1$reads)
# hist(ama1$reads)
# 
# #GONE
#####################################################################################3



###### 2) SAMPLE FILTERING #######
unique_loci_list <- merged_dfs_locifil_allefil %>%
  group_by(`Numero de estudo`, NIDA, Visita) %>%
  summarize(unique_loci = unique(locus)) %>%
  summarize(unique_loci_list = list(unique_loci))

print(unique_loci_list, n = 9999)

#remove samples with less than 50 loci sequenced
unique_loci_list$unique_loci_count <- sapply(unique_loci_list$unique_loci_list, length)

samples_to_remove <- unique_loci_list[unique_loci_list$unique_loci_count < 50,]$NIDA

unique_loci_list <- unique_loci_list[!unique_loci_list$NIDA %in% samples_to_remove,]
hist(unique_loci_list$unique_loci_count)

merged_dfs_locifil_allefil <- merged_dfs_locifil_allefil[!merged_dfs_locifil_allefil$NIDA %in% samples_to_remove,]
length(unique(merged_dfs_locifil_allefil$NIDA))


####### 3) PATIENT FILTERING #######
#remove patients that don't have all 4 visits
merged_dfs_locifil_allefil <- merged_dfs_locifil_allefil %>%
  group_by(`Numero de estudo`) %>%
  filter(n_distinct(Visita) == 4) %>%
  ungroup()

#check
l<- merged_dfs_locifil_allefil %>%
  group_by(`Numero de estudo`)%>%
  summarize(length(unique(Visita)))

print(paste0(length(unique(l$`Numero de estudo`))  ," patients with ", unique(l$`length(unique(Visita))`), " visits"))
#####################################################################################3



###### LOCI SELECTION #######
unique_loci_list <- merged_dfs_locifil_allefil %>%
  group_by(`Numero de estudo`, NIDA, Visita) %>%
  summarize(unique_loci = unique(locus)) %>%
  summarize(unique_loci_list = list(unique_loci))

print(unique_loci_list, n = 9999)

unique_loci_list$unique_loci_count <- sapply(unique_loci_list$unique_loci_list, length)

#loci present in all visits of all patients
all_elements <- unlist(unique_loci_list$unique_loci_list)
element_counts <- as.data.frame(table(all_elements))
element_counts <- element_counts[order(-element_counts$Freq), ]

common_loci <- Reduce(intersect, unique_loci_list$unique_loci_list) #no common loci across all visits/patients....

# keep loci present in at least 95% of samples (visits)
threshold <- round(0.995 * length(unique(merged_dfs_locifil_allefil$NIDA))) 
good_loci <- as.character(element_counts[element_counts$Freq >= threshold,]$all_elements)

#ama-1 as only good locus
# good_loci <- "Pf3D7_11_v3-1294284-1294580-1B"

# PENDING: FILTER OUT LOCI THAT ARE VERY HIGH IN SOME OF THE NEG CONTROLS (see bar plots of FILTERED data used here) =  unreliable loci.
unreliable_loci <- c("PmUG01_12_v1-1397996-1398245-1AB", "Pf3D7_03_v3-653961-654206-1A")

common_loci <- common_loci[!common_loci %in% unreliable_loci]
good_loci <- good_loci[!good_loci %in% unreliable_loci]


# select loci to use
if (length(common_loci) == 0){
  
  print("WARNING: no common loci across all visits of all patients. Choosing good loci instead:")
  cat("\n")
  
  #keep loci present in 95% of samples
  merged_dfs_locifil_allefil <- merged_dfs_locifil_allefil[merged_dfs_locifil_allefil$locus %in% good_loci,]
  
  print(good_loci)
  
} else {
  
  print(paste0(length(common_loci), " common loci found across all samples."))
  
  #keep common loci
  merged_dfs_locifil_allefil <- merged_dfs_locifil_allefil[merged_dfs_locifil_allefil$locus %in% common_loci,]
  
}


#####################################################################################3

#area plots:

unique(merged_dfs_locifil_allefil$allele)

# Convert `Visita` to a factor to ensure proper ordering in the plot
merged_dfs_locifil_allefil <- merged_dfs_locifil_allefil %>%
  mutate(Visita = factor(Visita, levels = unique(Visita)))

merged_dfs_locifil_allefil$Visita

merged_dfs_locifil_allefil <- merged_dfs_locifil_allefil %>%
  group_by(Visita) %>%
  arrange(norm.reads.locus, .by_group = TRUE) %>%
  ungroup()

#coloring
set.seed(69)
num_colors <- length(unique(merged_dfs_locifil_allefil$allele))
color_jump <- 1
rand_indices <- sample(1:num_colors, num_colors, replace = FALSE)
selected_colors <- rainbow(num_colors * color_jump)[rand_indices]


stack_bar <- ggplot(merged_dfs_locifil_allefil, aes(x = Visita, y = norm.reads.locus, fill = allele)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Visit", y = "Normalized Reads per Locus", fill = "Allele") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")+
  facet_wrap(~`Numero de estudo`, ncol = 8)+
  scale_fill_manual(values = selected_colors)

stack_bar

ggsave(paste0("stacked_barplot_MAF_", MAF, "_minreads_", min_reads, ".png"), stack_bar, dpi = 300, height = 12, width = 17, bg = "white")



###### ANALYSIS #######

#calculate the number of alleles for each locus per visit
unique_alleles <- merged_dfs_locifil_allefil %>%
  group_by(`Numero de estudo`, Visita) %>%
  summarize(alleles = list(allele),
            allele_count = length(allele))

print(unique_alleles, n =9999)


### CUMNULATIVE FINDING OF NEW ALLELES THROUGHOUT VISITS ###                                   ESTÁ MAL!!!!!!!!!!1

# Function to calculate the different alleles between two visits
get_different_alleles <- function(alleles1, alleles2) {
  setdiff(alleles2, alleles1)
}

participants <- unique(unique_alleles$`Numero de estudo`)

#participant <- "ASINT2-0044" #troubleshooting

allele_Accumulation <- list()

# Loop over each participant
for (participant in participants) {
  
  subset <- unique_alleles[unique_alleles$`Numero de estudo` == participant, ]
  
  # Check if there's only one visit for this participant
  if (n_distinct(subset$Visita) == 1) {
    subset$diff_from_previous <- NA # If there's only one visit, insert NA and move to the next participant
  } else {
    subset <- subset %>%
      arrange(Visita) %>%
      mutate(diff_from_previous = map(1:n(), function(i) {
        if (i == 1) {
          return(NA)
        } else {
          all_previous_alleles <- unique(unlist(subset$alleles[1:(i-1)]))
          get_different_alleles(all_previous_alleles, subset$alleles[[i]])
        }
      }))
  }
  
  allele_Accumulation[[participant]] <- subset
}

allele_Accumulation <- do.call(rbind, allele_Accumulation)

#count alleles
# Function to check the class of each element in the list column
check_class <- function(cell) {
  if (is.logical(cell)) {
    return(0)
  } else {
    return(length(cell))
  }
}

allele_Accumulation$diff_from_previous_counts <- unlist(map(allele_Accumulation$diff_from_previous, check_class))
print(allele_Accumulation)

## cummulative sum of alleles throughout time
cs <- allele_Accumulation %>% 
  group_by(`Numero de estudo`) %>%
  summarize(cumsum_diff_from_previous_counts = cumsum(diff_from_previous_counts))

allele_Accumulation <- cbind(allele_Accumulation, cumsum_diff_from_previous_counts = cs$cumsum_diff_from_previous_counts)


# #cumulative check by hand for "ASINT2-0002". is it equal to the loop? YES
# subset <- unique_alleles[unique_alleles$`Numero de estudo` == participants[1], ]
# 
# length(setdiff(subset$alleles[[2]], subset$alleles[[1]]))
# length(setdiff(subset$alleles[[3]], c(subset$alleles[[2]], subset$alleles[[1]])))
# length(setdiff(subset$alleles[[4]], c(subset$alleles[[3]], subset$alleles[[2]], subset$alleles[[1]])))
# length(setdiff(subset$alleles[[5]], c(subset$alleles[[4]], subset$alleles[[3]], subset$alleles[[2]], subset$alleles[[1]])))


#categorize infections
infection_results <- allele_Accumulation %>%
  group_by(`Numero de estudo`) %>%
  summarize(infection_status = ifelse(length(diff_from_previous_counts) ==  1, NA, # there was only 1 visit
                                      ifelse(sum(diff_from_previous_counts) > 0, "NEW_INFECTION", "same_infection")))


write.csv(infection_results, paste0("infection_results_MAF_", MAF, "_minreads_", min_reads,".csv"), row.names = F)
  

#visualization
allele_Accumulation_status <- merge(allele_Accumulation, infection_results, by=c("Numero de estudo"))

csplot <- ggplot(allele_Accumulation_status, aes(x = Visita, y = cumsum_diff_from_previous_counts, group = `Numero de estudo`, color = infection_status)) +
  geom_line(linewidth = 1.5, alpha = 0.5) +
  labs(x = "Visit", y = "Allele Accumulation Throughout Visits", title = paste0("MAF = ", MAF, "; Min reads: ", min_reads)) +
  theme_minimal() +
  facet_wrap(~`Numero de estudo`, ncol = 8)

csplot

ggsave(paste0("allele_Accumulation_MAF_", MAF, "_minreads_", min_reads, ".png"), csplot, dpi = 300, height = 12, width = 17, bg = "white")

#####################################################################################3

# # PCoA plot of cmmulative alleles through visits.
# 
# 
# pcoa_data <- allele_Accumulation_status[c(1, 2, 7,length(colnames(allele_Accumulation_status)))]
# 
# library(vegan)
# library(ape)
# library(tidyr)
# library(dplyr)
# library(ggplot2)
# 
# 
# # Reshape the data to wide format
# data_wide <- pcoa_data %>%
#   select(-infection_status) %>%
#   pivot_wider(names_from = Visita, values_from = cumsum_diff_from_previous_counts, values_fill = 0) %>%
#   as.data.frame()
# 
# #data_wide <- data_wide[,-2]
# 
# # Extract the infection status for coloring
# infection_status <- pcoa_data %>%
#   distinct(`Numero de estudo`, infection_status)
# 
# # Calculate the distance matrix using vegan's vegdist function
# dist_matrix <- vegdist(data_wide %>% select(-`Numero de estudo`), method = "bray")
# 
# # Perform PCoA using ape's pcoa function
# pcoa_result <- pcoa(dist_matrix)
# 
# # Create a data frame for the PCoA results
# pcoa_data <- data.frame(
#   Numero_de_estudo = data_wide$`Numero de estudo`,
#   PC1 = pcoa_result$points[, 1],
#   PC2 = pcoa_result$points[, 2]
# )
# 
# # Merge the PCoA results with the infection status
# colnames(infection_status)[1] <- "Numero_de_estudo"
# pcoa_data <- merge(pcoa_data, infection_status, by = "Numero_de_estudo")
# 
# # Plot the PCoA results
# ggplot(pcoa_data, aes(x = PC1, y = PC2, label = Numero_de_estudo, color = infection_status)) +
#   geom_point(size = 3, alpha =0.5) +
#   geom_text(vjust = 1.5, hjust = 1.5) +
#   labs(title = "PCoA of Allele Accumulation",
#        x = "Principal Coordinate 1",
#        y = "Principal Coordinate 2") +
#   theme_minimal() +
#   scale_color_discrete(name = "Infection Status")
# 









# #####################################################################################3
# 
# ######## FORMAT ######## 
# input_df <- merged_dfs %>%
#   group_by(PairsID, NIDA, locus, Site, time_point, Drug, days) %>%
#   summarize(total_n_alleles = length(locus)) %>%
#   ungroup()
# 
# #for ggplot 
# input_df <- input_df %>%
#   group_by(PairsID) %>%
#   mutate(max_days = max(days)) %>%
#   ungroup()
# 
# ## CALCULATE COI ##
# 
# allele_data$sample_id <- allele_data$NIDA
# 
# #saveRDS(allele_data, "allele_data_global_max_0_filtered_for_moire.RDS")
# #allele_data <- readRDS("allele_data_global_max_0_filtered_for_moire.RDS")
# 
# #library(moire)
# 
# #dat_filter <- moire::load_long_form_data(allele_data)
# 
# # set MOIRE parameterss
# #burnin <- 1e4
# #num_samples <- 1e4
# #pt_chains <- seq(1, .5, length.out = 20)
# 
# #mcmc_results <- moire::run_mcmc(
# #  dat_filter, is_missing = dat_filter$is_missing,
# #  verbose = T, burnin = burnin, samples_per_chain = num_samples,
# #  pt_chains = pt_chains, pt_num_threads = length(pt_chains),
# #  thin = 10
# #)
# 
# #saveRDS(mcmc_results, "allele_data_global_max_0_filtered_MOIRE-RESULTS.RDS")
# 
# mcmc_results <- readRDS("allele_data_global_max_0_filtered_MOIRE-RESULTS.RDS")
# eff_coi <- moire::summarize_effective_coi(mcmc_results)
# naive_coi <- moire::summarize_coi(mcmc_results)
# relatedness <- moire::summarize_relatedness(mcmc_results)
# 
# #change NA and NaN to 1 (means the strain was or became monoclonal)
# relatedness <- relatedness %>%
#   mutate_all(~ ifelse(is.na(.x) | is.nan(.x), 1, .x))
# 
# colnames(eff_coi)[1]<-"NIDA"
# colnames(naive_coi)[1]<-"NIDA"
# colnames(relatedness)[1]<-"NIDA"
# 
# input_df <- merge(input_df[, c("PairsID","NIDA", "time_point")], eff_coi[, c("NIDA","post_effective_coi_med")], by="NIDA")
# input_df <- merge(input_df, naive_coi[, c("NIDA","naive_coi", "post_coi_med")], by="NIDA")
# input_df <- merge(input_df, relatedness[, c("NIDA","post_relatedness_med")], by="NIDA")
# 
# # Plot the post_effective_coi_mean with dodged bars
# coicoi <- ggplot(input_df, aes(x = factor(PairsID), y = post_effective_coi_med, fill = time_point)) +
#   geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6) +
#   labs(x = "PairsID", y = "post_effective_coi_med") +
#   #scale_fill_manual(values = c("D0_alleles" = "darkgreen", "Dx_alleles" = "darkorange", "shared_alleles" = "purple")) +
#   theme_minimal() +
#   theme(legend.position = "top")
# 
# coicoi
# 
# ggsave("eCOI.png", coicoi, width = 16, height = 9, bg = "white")
# 
# #############################################################################
# 
# 
# ### HACER TABLA CON PAIRID, DX, PROVINCE, DRUG, %SHARED_ALLELES (D0/shared_alleles), %CHANGE_ALLELES, %CHANGE_COI, %CHANGE_RELATEDNESS
# 
# #change in post_effective_coi_med and post_relatedness_med
# ok_all<-merged_dfs %>%
#   group_by(PairsID, NIDA, Site, time_point, Drug, days) %>%
#   summarize(total_n_alleles = length(locus)) %>%
#   ungroup()%>%
#   arrange(PairsID, time_point)
# 
# ok <- merge(ok_all, eff_coi[, c("NIDA","post_effective_coi_med")], by="NIDA")
# #ok <- merge(ok, relatedness[, c("NIDA","post_relatedness_med")], by="NIDA")
# 
# ok <- ok %>%
#   group_by(PairsID) %>%
#   summarise(eCOI_change = (post_effective_coi_med[time_point == 'Dx'] - post_effective_coi_med[time_point == 'D0'])/ post_effective_coi_med[time_point == 'D0']) #,
#             #relatedness_change = (post_relatedness_med[time_point == 'Dx'] - post_relatedness_med[time_point == 'D0'])/ post_relatedness_med[time_point == 'D0'])
# 
# #change in allele counts
# #ok$alleles_change <- (changes_per_pairs$Dx_alleles - changes_per_pairs$D0_alleles) /  changes_per_pairs$D0_alleles
# 
# #percentage of shared alleles from D0 (INITIAL ALLELES)
# ok$percentage_shared_alleles_from_D0 <- changes_per_pairs$shared_alleles/changes_per_pairs$D0_alleles
# 
# #percentage of shared alleles from Dx (NEW ALLELES)
# ok$percentage_shared_alleles_from_Dx <- changes_per_pairs$shared_alleles/changes_per_pairs$Dx_alleles
# 
# #lost alleles
# ok$lost_alleles <- changes_per_pairs$D0_alleles-changes_per_pairs$shared_alleles
# 
# #gained alleles
# ok$gained_alleles <- changes_per_pairs$Dx_alleles-changes_per_pairs$shared_alleles
# 
# #lost_gained_alleles_ratio
# #ok$lost_gained_alleles_ratio <- (ok$lost_alleles-ok$gained_alleles)/ok$gained_alleles
# 
# #ecoi before
# ecoias<- input_df %>%
#   select(PairsID, time_point, post_effective_coi_med)%>%
#   distinct()%>%
#   arrange(PairsID, time_point)
# 
# #ok$eCOI_D0<- ecoias$post_effective_coi_med[ecoias$time_point == "D0"]
# 
# #ecoi after
# #ok$eCOI_Dx <- ecoias$post_effective_coi_med[ecoias$time_point == "Dx"]
# 
# #percentages_ratio
# ok$percentage_alleles_ratio <- ok$percentage_shared_alleles_from_Dx/ok$percentage_shared_alleles_from_D0
# 
# #net_allele_change
# #ok$net_allele_change <- ok$gained_alleles-ok$lost_alleles
# 
# 
# # Create a data frame from the correlation matrix
# cor_matrix <- cor(ok[, c(-1)])
# cor_matrix
# cor_data <- as.data.frame(cor_matrix)
# cor_data <- melt(cor_matrix, id.vars = "PairsID")
# 
# # Create a correlation heatmap
# ggplot(cor_data, aes(x = Var1, y = Var2, fill = value)) +
#   geom_tile() +
#   theme_minimal() +
#   labs(title = "Correlation Heatmap")
# 
# hist(ok$eCOI_change); median(ok$eCOI_change); quantile(ok$eCOI_change, c(0.95, 0.05))
# hist(ok$percentage_alleles_ratio); median(ok$percentage_alleles_ratio); quantile(ok$percentage_alleles_ratio, c(0.95, 0.05))
# 
# plot(ok$eCOI_change, ok$percentage_alleles_ratio) 
# 
# hist(ok$lost_alleles); median(ok$lost_alleles); quantile(ok$lost_alleles, c(0.95, 0.05))
# hist(ok$gained_alleles); median(ok$gained_alleles); quantile(ok$gained_alleles, c(0.95, 0.05))
# 
# 
# # Loop through each column and calculate statistics for the thresholds
# summary_stats <- data.frame()
# 
# for (col_name in names(ok)) {
#   if (is.numeric(ok[[col_name]])) {
#     stats <- ok %>%
#       summarise(
#         Variable = col_name,
#         Mean = mean(.data[[col_name]], na.rm = TRUE),
#         SD = sd(.data[[col_name]], na.rm = TRUE),
#         lower_sigma = Mean - SD,
#         upper_sigma = Mean + SD,
#         Median = median(.data[[col_name]], na.rm = TRUE),
#         q25 = Median - quantile(.data[[col_name]], 0.25),
#         q75 = Median + quantile(.data[[col_name]], 0.25)
#       )
#     summary_stats <- bind_rows(summary_stats, stats)
#   }
# }
# 
# summary_stats<-summary_stats[-1,]
# summary_stats
# 
# #all distributions mostly normal (use mean, sd), except eCOI_change (use median, q25-q75). 
# #sapply(ok, hist)
# 
# ## BUILD CATEGORY TABLE
# ecoi_cats <- ifelse(ok$eCOI_change > summary_stats$q25[1], "increase", 
#                     ifelse(ok$eCOI_change < summary_stats$q75[1], "decrease", "mostly_neutral"))
#   
# lost_alleles_cats <- ifelse(ok$lost_alleles > summary_stats$upper_sigma[4], "sharp",
#                             ifelse(ok$lost_alleles < summary_stats$lower_sigma[4], "weak", "mid"))
# 
# gained_alleles_cats <- ifelse(ok$gained_alleles > summary_stats$upper_sigma[5], "sharp",
#                             ifelse(ok$gained_alleles < summary_stats$lower_sigma[5], "weak", "mid"))
# 
# percentage_allele_ratio_cats <- ifelse(ok$percentage_alleles_ratio > summary_stats$upper_sigma[6], "high",
#                               ifelse(ok$percentage_alleles_ratio < summary_stats$lower_sigma[6], "low", "mostly_neutral"))
# 
# CATEGORY_TABLE <- as.data.frame(cbind(PairsID = ok$PairsID, ecoi_cats, lost_alleles_cats, gained_alleles_cats, percentage_allele_ratio_cats))
# CATEGORY_TABLE
# 
# ### ADD INTERPRETATIONS TO CATEGORY TABLE
# 
# #unique combinations in the TES22 run
# un<-c(NULL)
# for (r in 1: length(rownames(CATEGORY_TABLE))){
#   on <- paste(CATEGORY_TABLE$ecoi_cats[r], 
#               CATEGORY_TABLE$lost_alleles_cats[r], 
#               CATEGORY_TABLE$gained_alleles_cats[r], 
#               CATEGORY_TABLE$percentage_allele_ratio_cats[r], sep = "_")
#   un <-c(on, un)
# }
# 
# un <- unique(un)
# un<-as.data.frame(sort(un))
# 
# #without filtering
# interpretations <-c("recrudescence", "recrudescence", "recrudescence_STRONG", "recrudescence", 
#                     "recrudescence_STRONG", "rec/reinf(additive)", "reinfection(additive)", "rec/reinf(additive)", 
#                     "reinfection(additive)_STRONG", "rec/reinf(additive)","rec/reinf(additive)", "recrudescence_TOTAL")
#  
# #with filtering: CHECK FOR SPECIFIC ONE OR HARD CODE THE EQUIVALENCES
# #interpretations <-c("recrudescence", "recrudescence", "recrudescence_STRONG", "recrudescence", 
# #                    "recrudescence_STRONG", "rec/reinf(additive)", "reinfection(additive)", "rec/reinf(additive)", "rec/reinf(additive)",
# #                    "reinfection(additive)_STRONG", "rec/reinf(additive)","rec/reinf(additive)", "recrudescence_TOTAL")
# 
# interpretations <- as.data.frame(cbind(changes=un, interpretations))
# colnames(interpretations)[1]<-"changes"
# 
# CATEGORY_TABLE$full <- paste(CATEGORY_TABLE$ecoi_cats, 
#       CATEGORY_TABLE$lost_alleles_cats, 
#       CATEGORY_TABLE$gained_alleles_cats, 
#       CATEGORY_TABLE$percentage_allele_ratio_cats, sep = "_")
# 
# 
# 
# matches <- match(CATEGORY_TABLE$full, interpretations$changes)
# matching_rows <- !is.na(matches)
# CATEGORY_TABLE$infection[matching_rows] <- interpretations$interpretations[matches[matching_rows]]
# 
# 
# #############################################################################
# # ADDING PHASED RESISTANCE HAPLOTYPES TO THE CATEGORIES
# #############################################################################
# phased_haplos <- read.csv("TESS22_phased_haplos_correct_only.csv")
# phased_haplos <- phased_haplos[,c(-2:-8, -12:-17)]
# colnames(phased_haplos)[1]<- "NIDA"
# 
# phased_haplos$NIDA <- sub("_S.*", "", phased_haplos$NIDA)
# 
# phased_haplos <- merge(phased_haplos, merged_dfs[, c("NIDA", "PairsID", "time_point")], by="NIDA")
# 
# phased_haplos <- distinct(phased_haplos)
# 
# phased_haplos <- phased_haplos %>% arrange(PairsID, time_point)
# 
# haplos_plot <- ggplot(phased_haplos, aes(x = time_point, y = HAPLO_FREQ_RECALC, fill = haplotype)) +
#   geom_bar(stat = "identity", position = "stack") +
#   labs(x = "Time Point", y = "HAPLO_FREQ_RECALC") +
#   theme_minimal() +
#   facet_wrap(~PairsID)
# 
# haplos_plot
# 
# ggsave("haplos.png", haplos_plot, width = 14, height = 12, bg = "white")
# 
# # Create an empty dataframe to store the results
# phased_haplos_cats <- data.frame(PairsID = character(0), overlap = character(0), change = character(0))
# 
# # Iterate over unique PairsID values
# for (pair in unique(phased_haplos$PairsID)) {
#   subsettt <- phased_haplos[phased_haplos$PairsID == pair, ]
#   D0 <- subsettt[subsettt$time_point == "D0", "haplotype"]
#   Dx <- subsettt[subsettt$time_point == "Dx", "haplotype"]
#   
#   # Determine the values of i and j
#   i <- ifelse(length(intersect(unique(D0), unique(Dx)) > 0), "overlapping", "non-overlapping")
#   j <- ifelse(length(unique(D0)) < length(unique(Dx)), "diversity_increase", 
#               ifelse(length(unique(D0)) == length(unique(Dx)), "diversity_same", "diversity_decrease"))
#   
#   # Add the values to the results dataframe
#   phased_haplos_cats <- rbind(phased_haplos_cats, data.frame(PairsID = pair, overlap = i, change = j))
# }
# 
# #haplotype categories by change: overlapping between times? same number of haplos between times?
# phased_haplos_cats$final_cats <- paste(phased_haplos_cats$overlap, phased_haplos_cats$change, sep ="_")
# 
# #update CATEGORY_TABLE: if no haplos shared between infections, category must be changed to reinfection_TOTAL
# reinfection_TOTAL <- CATEGORY_TABLE[phased_haplos_cats$overlap == "non-overlapping",]["infection"]
# 
# idx_reinfection_TOTAL <- as.numeric(rownames(reinfection_TOTAL))
# CATEGORY_TABLE[idx_reinfection_TOTAL,]["infection"] <- "reinfection_TOTAL"
# 
# #some other categorical variables:
# ok_all_subset<-ok_all[ok_all$time_point == "Dx",]
# ok_all_subset$Drug_days <- paste(ok_all_subset$Drug, ok_all_subset$days, sep ="_") 
# 
# 
# #export variables table
# FINAL_VARIABLES_TABLE<- cbind(ok, CATEGORY_TABLE[,-1], ok_all_subset[,-1:-2], phased_haplos_cats[,-1])
# write.csv(FINAL_VARIABLES_TABLE, "FINAL_VARIABLES_TABLE.csv", row.names=F)
# 
# 
# # ##############################
# # #reinfection_TOTAL pairs may be taken as a separate category from the rest because there are no shared haplos between the samples and thus it's CLEAR they are reinfections
# # #calculate the pca without them (if want to include them in pca, comment this chunk of code)
# # CATEGORY_TABLE <- CATEGORY_TABLE[-idx_reinfection_TOTAL, ]
# # length(row.names(CATEGORY_TABLE))
# # ok <- ok[-idx_reinfection_TOTAL, ]
# # length(row.names(ok))
# # ############
# # #keep overlapping haplos only. if needed to include reinfection_TOTAL category (non-overlapping haplos), remove...
# # # ...this code and change the ggbiplot command on the shape variable for: shape = factor(phased_haplos_cats$final_cats)
# # overlapping_haplos<-phased_haplos_cats$final_cats[-idx_reinfection_TOTAL]
# # ############
# # ###############################
# # 
# # #pca
# # pca_data <- ok[, c(-1, -3, -4)] ## EXCLUDE ALLELES_CHANGE SINCE IS 97% CORRELATED WITH ECOI_CHANGE
# # 
# # # Standardize the data (centering and scaling)
# # pca_data <- scale(pca_data)
# # 
# # # Perform PCA
# # pca_result <- prcomp(pca_data)
# # 
# # # Access the principal component scores
# # pca_scores <- pca_result$x
# # 
# # # View the summary of the PCA
# # summary(pca_result)
# # 
# # # Variance explained by each principal component
# # variance_explained <- (pca_result$sdev^2) / sum(pca_result$sdev^2)
# # variance_explained
# # 
# # pca_df <- as.data.frame(pca_scores)
# # 
# # # Add PairsID if you want to distinguish data points by PairsID
# # pca_df$eCOI_change <- ok$eCOI_change
# # pca_df$PairsID <- ok$PairsID
# # 
# # #plot
# # biplot_ggplot_2 <- ggplot(pca_df, aes(PC1, PC2, color =CATEGORY_TABLE$infection, shape = factor(overlapping_haplos))) + #shape = factor(phased_haplos_cats$final_cats)
# #   geom_point(size=6, alpha =0.65) +
# #   geom_hline(yintercept = 0, linetype = "dashed", alpha =0.25) +
# #   geom_vline(xintercept = 0, linetype = "dashed",  alpha =0.25) +
# #   geom_text(aes(label = PairsID), hjust = -0.4, vjust = -0.3) +
# #   labs(x = paste("PC1 (", round(variance_explained[1] * 100, 2), "%)"),
# #        y = paste("PC2 (", round(variance_explained[2] * 100, 2), "%)")) +
# #   theme_minimal()
# # 
# # biplot_ggplot_2
# # 
# # ggsave("PCA_res_haplo-cat.png", biplot_ggplot_2, width = 16, height = 9, bg = "white")
# 
# #################################################################################################3
# 
# 
# #######################################################################
# #################   AMPLICON BY AMPLICON METHOD #######################
# #######################################################################
# 
# #GRUNENBERG et al., 2019-LIKE (see sources folder)
# 
# 
# # 0.- filter out amplicons with MAF < 0.05 and with indels (see pseudocigar)
# # 1.- calculate He for each amplicon using all samples from D0 (could be by province or use all runs to see variance) (CURRENTLY DOING ONLY ALLELE COUNTS 'CAUSE CALCULATING HE REQUIRES A LOT OF FORMATTING...)
# # 2.- rank amplicons by He descendingly
# # 3.- pick most variable amplicons (USING >= Q95 atm)
# # 4.- plot shit up
# # 5.- categorize infections (use grunenberg 2019's definitions and come up with own ones: recrudescence (additive[aka reinfection], neutral, substractive), reinfection) 
# 
# 
# # 0.- filter out amplicons with MAF < 0.05 and with indels (see pseudocigar)
# merged_dfs_filtered <- merged_dfs[merged_dfs$norm.reads.locus >= 0.05,]
# merged_dfs_filtered <- merged_dfs_filtered[!grepl("I=", merged_dfs_filtered$pseudo_cigar),] #remove alleles with I (insertion)
# merged_dfs_filtered <- merged_dfs_filtered[!grepl("D=", merged_dfs_filtered$pseudo_cigar),] #remove alleles with D (deletion)
# 
# merged_dfs_D0_only <- merged_dfs_filtered[merged_dfs_filtered$time_point=="D0",]
# 
# amp_count <- merged_dfs_D0_only %>% #amplicons below length(unique(merged_dfs_D0_only$SampleID)) are not suited to use as markers because they can be anset sometimes.
#   group_by(locus)%>%
#   summarize(n_amplicons = length(unique(SampleID)))
# 
# amp_count_OK <- amp_count$locus[amp_count$n_amplicons == max(amp_count$n_amplicons)]
# 
# #clean amplicons
# merged_dfs_D0_only <- merged_dfs_D0_only[merged_dfs_D0_only$locus %in% amp_count_OK,]
# 
# #########################################3
# #OPTION HETEROZYGOSITY
# #############################################
# 
# # 1.- calculate Heterozygosity for each amplicon using all samples from D0 (could be by province or use all runs to see variance) (CURRENTLY DOING ONLY ALLELE COUNTS 'CAUSE CALCULATING HE REQUIRES A LOT OF FORMATTING...)
# heterozygosity <- merged_dfs_D0_only %>%
#   group_by(SampleID, locus) %>%
#   summarize(
#     heterozygosity = 1 - sum(norm.reads.locus^2)
#   )
# 
# mean_HE<- tapply(heterozygosity$heterozygosity, heterozygosity$locus,  function(x) mean(x))
# 
# mean_HE <- as.data.frame(mean_HE)
# mean_HE$locus<- rownames(mean_HE)
# 
# #remove amplicons w/ meanHe = 0
# mean_HE <- mean_HE[mean_HE$mean_HE != 0,] 
# 
# # 2.- rank amplicons by H descendingly
# mean_HE<- mean_HE[order(-mean_HE$mean_HE), ]
# 
# # 3.- pick most variable amplicons #benchmarking BELOW!!
# 
# # Plot the histogram
# png("histogram.png", width = 800, height = 600)
# hist(mean_HE$mean_HE, main = "Mean He per amplicon", xlab = "mean He", ylab = "# Amplicons")
# #abline(v = thresh_HE, col = "red", lty = 1)
# #text(thresh_HE + 0.02, 60, "0.9 quantile", col = "red")
# dev.off()
# 
# #this is according to benchmarking, which is at the end of the script
# variable_alleles_D0 <- mean_HE$locus[1:184] #mean_HE$locus[mean_HE$mean_HE >= thresh_HE]
# 
# merged_dfs_filtered_variable_loci <- merged_dfs_filtered[merged_dfs_filtered$locus %in% variable_alleles_D0, ]
# 
# amp_Categories<-unique(merged_dfs_filtered_variable_loci[merged_dfs_filtered_variable_loci$locus %in% variable_alleles_D0,][c("locus", "Category")]) #only diversity and immune alleles (ofc also amplicons) left. good.
# amp_Categories<- merge(amp_Categories,mean_HE, by="locus" )
# amp_Categories<- amp_Categories[order(-amp_Categories$mean_HE), ]
# 
# write.csv(amp_Categories, "most_variable_amplicons_benchmarked.csv", row.names = F)
# 
# #############################################
# #OPTION AMOUNT OF ALLELES
# #############################################
# 
# # # 1.- calculate amount of alleles for each amplicon using all samples from D0 (could be by province or use all runs to see variance) (CURRENTLY DOING ONLY ALLELE COUNTS 'CAUSE CALCULATING HE REQUIRES A LOT OF FORMATTING...)
# # 
# # alleles_per_locus <- merged_dfs_D0_only %>%
# #   group_by(locus) %>%
# #   summarise(length(unique(pseudo_cigar)))
# # 
# # # 2.- rank amplicons by amount of alleles descendingly
# # alleles_per_locus<- alleles_per_locus[order(-alleles_per_locus$`length(unique(pseudo_cigar))`), ]
# # 
# # 
# # # 3.- pick most variable amplicons (USING >= Q90 atm)
# # thresshhh <- quantile(alleles_per_locus$`length(unique(pseudo_cigar))`, 0.90) #above this, alleles are considered highly variable
# # 
# # png("histogram.png", width = 800, height = 600)
# # # Plot the histogram
# # hist(alleles_per_locus$`length(unique(pseudo_cigar))`, main = "Unique alleles per each amplicon", xlab = "# Alleles", ylab = "# Amplicons")
# # abline(v = thresshhh, col = "red", lty = 1)
# # text(thresshhh + 2.2, 120, "0.9 quantile", col = "red")
# # dev.off()
# # 
# # variable_alleles_D0 <- alleles_per_locus$locus[alleles_per_locus$`length(unique(pseudo_cigar))`>= thresshhh]
# # 
# # merged_dfs_filtered_variable_loci <- merged_dfs_filtered[merged_dfs_filtered$locus %in% variable_alleles_D0, ]
# # 
# # amp_Categories<-unique(merged_dfs_filtered_variable_loci[merged_dfs_filtered_variable_loci$locus == variable_alleles_D0,][c("locus", "Category")]) #only diversity and immune alleles (ofc also amplicons) left. good.
# # amp_Categories<- merge(amp_Categories,alleles_per_locus, by="locus" )
# # amp_Categories<- amp_Categories[order(-amp_Categories$`length(unique(pseudo_cigar))`), ]
# # colnames(amp_Categories)[3]<- "unique_alleles"
# # 
# # write.csv(amp_Categories, "most_variable_amplicons_q90.csv", row.names = F)
# 
# #hist(merged_dfs_filtered_variable_loci$norm.reads.locus)
# #hist(merged_dfs_filtered_variable_loci$n.alleles)
# #########################################3 #########################################3 #########################################3 #########################################3
# 
# # 4.- plotting shit up
# #### BY AMPLICON ###
# #plot_list_amplicons <- list()
# 
# #for (amp in unique(merged_dfs_filtered_variable_loci$locus)) {
# #  subset_merged_dfs_filtered_variable_loci <- merged_dfs_filtered_variable_loci[merged_dfs_filtered_variable_loci$locus == amp,]
# 
# # Create individual ggplot
# #  single_plot <- ggplot(subset_merged_dfs_filtered_variable_loci, aes(x = time_point, y = norm.reads.locus, fill = pseudo_cigar)) +
# #    geom_bar(stat = "identity", position = "stack") +
# #    labs(x = "Time Point", y = "norm.read.locus") +
# #    theme_minimal() +
# #    facet_wrap(~PairsID) +
# #    ggtitle(paste("Amplicon:", amp))+
# #    guides(fill = guide_legend(ncol = 2))
# 
# # Append the plot to the list
# #  plot_list_amplicons[[length(plot_list_amplicons) + 1]] <- single_plot
# #}
# #plot_list_amplicons[[6]]
# 
# ### BY PAIR OF SAMPLES ###
# plot_list_PairsID <- list()
# 
# for (pair in unique(merged_dfs_filtered_variable_loci$PairsID)) {
#   
#   subset_merged_dfs_filtered_variable_loci <- merged_dfs_filtered_variable_loci[merged_dfs_filtered_variable_loci$PairsID == pair,]
#   
#   #coloring
#   set.seed(69)
#   num_colors <- length(unique(subset_merged_dfs_filtered_variable_loci$pseudo_cigar))
#   color_jump <- 1
#   rand_indices <- sample(1:num_colors, num_colors, replace = FALSE)
#   selected_colors <- rainbow(num_colors * color_jump)[rand_indices]
#   
#   #plotting
#   single_plot <- ggplot(subset_merged_dfs_filtered_variable_loci, aes(x = time_point, y = norm.reads.locus, fill = pseudo_cigar)) +
#     geom_bar(stat = "identity", position = "stack") +
#     labs(x = "Time Point", y = "norm.read.locus") +
#     theme_minimal() +
#     facet_wrap(~locus) +
#     ggtitle(paste("PairID:", pair))+
#     guides(fill = guide_legend(ncol = 2))+
#     scale_fill_manual(values = selected_colors)+
#     guides(fill = "none") #remove to see alleles in legend
#   
#   #ggsave(paste("PairID_", pair,".png"), single_plot, width = 30, height = 30, bg = "white")
#   
#   # Append the plot to the list
#   plot_list_PairsID[[length(plot_list_PairsID) + 1]] <- single_plot
# }
# 
# #plot_list_PairsID[[17]] #46 is plot_list_PairsID[[3]], 47 es plot_list_PairsID[[9]]
# 
# 
# # 5.- categorize infections (use grunenberg 2019's definitions and come up with own ones: recrudescence (additive[aka reinfection], neutral, substractive), reinfection) 
# 
# common_alleles_table <-merged_dfs_filtered_variable_loci %>%
#   filter(time_point %in% c("D0", "Dx")) %>%
#   group_by(PairsID, locus) %>%
#   summarise(alleles_D0 = length(unique(pseudo_cigar[time_point == "Dx"])),
#             alleles_Dx = length(unique(pseudo_cigar[time_point == "Dx"])),
#             common_alleles_D0_Dx = n_distinct(pseudo_cigar[pseudo_cigar %in% intersect(unique(pseudo_cigar[time_point == "D0"]), unique(pseudo_cigar[time_point == "Dx"]))]))
# 
# common_alleles_table
# 
# #criteria:
# #what does each marker say? (using NI for "New Infection" == "Reinfection")
# #Classification of recrudescence (R) required that at least one of the haplotypes occurred in both, the pre- and post-treatment sample, with a minimum haplotype frequency of 1% in at least two of the three independent replicates performed. 
# #A new infection (NI) was defined by the occurrence of only new haplotypes in the post-treatment sample with a minimum haplotype frequency of 1%.
# common_alleles_table$infection <- ifelse(common_alleles_table$common_alleles_D0_Dx > 0, "R", "NI")
# 
# #what does each marker say?
# infection_results <- common_alleles_table %>%
#   group_by(PairsID) %>%
#   summarise(Recrudescence_markers = sum(infection == "R"),
#            Reinfection_markers = sum(infection == "NI"),
#            perc_Recrudescence_markers = Recrudescence_markers/(Recrudescence_markers + Reinfection_markers),
#            perc_Reinfection_markers = Reinfection_markers/(Recrudescence_markers + Reinfection_markers))
# 
# # add phasing data
# infection_results$phasing_dhfr_dhps_haplo <- phased_haplos_cats$final_cats
# 
# # in case markers disagree...
# # (i) The WHO/MMV recommended approach1, where a recurrent parasitemia is classified in the overall outcome as NI, if at least one of the three markers had given a NI result
# infection_results$WHO_MMV_cireteria <- ifelse(infection_results$Reinfection_markers > 0, "NI", "R")
# 
# # (ii) The alternative, newer approach termed “2/3 algorithm”8. This algorithm uses the consensus of two markers as final PCR-correction outcome.
# infection_results$two_of_three_algo <- ifelse(infection_results$perc_Reinfection_markers >= (2/3), "NI", 
#                                               ifelse(infection_results$perc_Recrudescence_markers >= (2/3), "R", "?"))
# 
# #two_three_phasing algo:
# # 1) if two_of_three_algo == ? AND phasing_dhfr_dhps_haplo == "non-overlapping_diversity_same", "NI"
# # 2) if two_of_three_algo == ? AND phasing_dhfr_dhps_haplo == "overlapping_diversity_same", "R"
# # 3) if two_of_three_algo == ? AND phasing_dhfr_dhps_haplo == "overlapping_diversity_increase", "NI" 
# # 4) if two_of_three_algo == ? AND phasing_dhfr_dhps_haplo == "overlapping_diversity_decrease", "R"
# 
# infection_results$two_three_phasing <- infection_results$two_of_three_algo 
# 
# infection_results$two_three_phasing <- ifelse(infection_results$two_of_three_algo == "?" & infection_results$phasing_dhfr_dhps_haplo == "non-overlapping_diversity_same", "NI", infection_results$two_three_phasing) 
# #infection_results$two_three_phasing <- ifelse(infection_results$two_of_three_algo == "?" & infection_results$phasing_dhfr_dhps_haplo == "overlapping_diversity_same", "?", infection_results$two_three_phasing) 
# #infection_results$two_three_phasing <- ifelse(infection_results$two_of_three_algo == "?" & infection_results$phasing_dhfr_dhps_haplo == "overlapping_diversity_increase", "NI", infection_results$two_three_phasing) 
# #infection_results$two_three_phasing <- ifelse(infection_results$two_of_three_algo == "?" & infection_results$phasing_dhfr_dhps_haplo == "overlapping_diversity_decrease", "R", infection_results$two_three_phasing) 
# 
# #two_three_majority_rule_phasing algo:
# infection_results$two_three_majority_rule_phasing <- infection_results$two_three_phasing 
# 
# infection_results$two_three_majority_rule_phasing <- ifelse(infection_results$two_three_majority_rule_phasing == "?",
#   ifelse(infection_results$perc_Recrudescence_markers > infection_results$perc_Reinfection_markers, "R","NI"),
#   infection_results$two_three_majority_rule_phasing)
# 
# #two_three_majority_rule_phasing algo:
# infection_results$two_three_phasing_WHO_MMV <- infection_results$two_three_phasing 
# 
# infection_results$two_three_phasing_WHO_MMV <- ifelse(infection_results$two_three_phasing == "?", infection_results$WHO_MMV_cireteria, infection_results$two_three_majority_rule_phasing )
# 
# infection_results
# 
# write.csv(infection_results, "infection_results.csv", row.names =F)
# 
# ####################################################################################################################3
# ###PCA CON NUEVO INPUT ###
# ####################################################################################################################3
# 
# ### HACER TABLA filtrada CON PAIRID, DX, PROVINCE, DRUG, %SHARED_ALLELES (D0/shared_alleles), %CHANGE_ALLELES, %CHANGE_COI, %CHANGE_RELATEDNESS
# 
# ok_2<-merged_dfs_filtered_variable_loci %>%
#   group_by(PairsID, NIDA, Site, time_point, Drug, days) %>%
#   summarize(total_n_alleles = length(locus)) %>%
#   ungroup()%>%
#   arrange(PairsID, time_point)
# 
# ok_2 <- merge(ok_2, eff_coi[, c("NIDA","post_effective_coi_med")], by="NIDA")
# #ok <- merge(ok, relatedness[, c("NIDA","post_relatedness_med")], by="NIDA")
# 
# ok_2 <- ok_2 %>%
#   group_by(PairsID) %>%
#   summarise(eCOI_change = (post_effective_coi_med[time_point == 'Dx'] - post_effective_coi_med[time_point == 'D0'])/ post_effective_coi_med[time_point == 'D0']) #,
# 
# #percentage of shared alleles from D0 (INITIAL ALLELES)
# ok_2$percentage_shared_alleles_from_D0 <- changes_per_pairs$shared_alleles/changes_per_pairs$D0_alleles
# 
# #percentage of shared alleles from Dx (NEW ALLELES)
# ok_2$percentage_shared_alleles_from_Dx <- changes_per_pairs$shared_alleles/changes_per_pairs$Dx_alleles
# 
# #lost alleles
# ok_2$lost_alleles <- changes_per_pairs$D0_alleles-changes_per_pairs$shared_alleles
# 
# #gained alleles
# ok_2$gained_alleles <- changes_per_pairs$Dx_alleles-changes_per_pairs$shared_alleles
# 
# #ecoi before
# ecoias_2<- input_df %>%
#   select(PairsID, time_point, post_effective_coi_med)%>%
#   distinct()%>%
#   arrange(PairsID, time_point)
# 
# #percentages_ratio
# ok_2$percentage_alleles_ratio <- ok_2$percentage_shared_alleles_from_Dx/ok_2$percentage_shared_alleles_from_D0
# 
# 
# #pca
# pca_data <- ok_2[, c(-1, -5:-7)] ## EXCLUDE ALLELES_CHANGE SINCE IS 97% CORRELATED WITH ECOI_CHANGE
# 
# # Standardize the data (centering and scaling)
# pca_data <- scale(pca_data)
# 
# # Perform PCA
# pca_result <- prcomp(pca_data)
# 
# # Access the principal component scores
# pca_scores <- pca_result$x
# 
# # View the summary of the PCA
# summary(pca_result)
# 
# # Variance explained by each principal component
# variance_explained <- (pca_result$sdev^2) / sum(pca_result$sdev^2)
# variance_explained
# 
# pca_df <- as.data.frame(pca_scores)
# 
# # Add PairsID if you want to distinguish data points by PairsID
# pca_df$eCOI_change <- ok_2$eCOI_change
# pca_df$PairsID <- ok_2$PairsID
# 
# #plot
# biplot_ggplot_3 <- ggplot(pca_df, aes(PC1, PC2, color = infection_results$WHO_MMV_cireteria)) + #shape = factor(phased_haplos_cats$final_cats)
#   geom_point(size=6, alpha =0.65) +
#   geom_hline(yintercept = 0, linetype = "dashed", alpha =0.25) +
#   geom_vline(xintercept = 0, linetype = "dashed",  alpha =0.25) +
#   geom_text(aes(label = PairsID), hjust = -0.4, vjust = -0.3) +
#   labs(x = paste("PC1 (", round(variance_explained[1] * 100, 2), "%)"),
#        y = paste("PC2 (", round(variance_explained[2] * 100, 2), "%)")) +
#   ggtitle("WHO_MMV_NGS")+
#   theme_minimal() + 
#   guides(color = FALSE, shape = FALSE)
# 
# biplot_ggplot_3
# 
# biplot_ggplot_4 <- ggplot(pca_df, aes(PC1, PC2, color = infection_results$two_three_phasing_WHO_MMV)) + #shape = factor(phased_haplos_cats$final_cats)
#   geom_point(size=6, alpha =0.65) +
#   geom_hline(yintercept = 0, linetype = "dashed", alpha =0.25) +
#   geom_vline(xintercept = 0, linetype = "dashed",  alpha =0.25) +
#   geom_text(aes(label = PairsID), hjust = -0.4, vjust = -0.3) +
#   labs(x = paste("PC1 (", round(variance_explained[1] * 100, 2), "%)"),
#        y = paste("PC2 (", round(variance_explained[2] * 100, 2), "%)")) +
#   ggtitle("two_three_phasing_WHO_MMV_NGS")+
#   theme_minimal() +
#   guides(color = FALSE, shape = FALSE)
# 
# biplot_ggplot_4
# 
# #ggsave("PCA_res_haplo_NI_R_WHO_MMV_cireteria.png", biplot_ggplot_3, width = 16, height = 9, bg = "white")
# ggsave("PCA_res_haplo_NI_R_WHO_MMV.png", biplot_ggplot_3, width = 16, height = 9, bg = "white")
# 
# ##check Arlindo's CE categories
# CE_cats <- read.csv("CE_categories_arlindo.csv", sep ="\t")
# CE_cats <- CE_cats[c(-20,-30),]
# 
# biplot_ggplot_5 <- ggplot(pca_df, aes(PC1, PC2, color = CE_cats$WHO_2_3_CE)) + #shape = factor(phased_haplos_cats$final_cats)
#   geom_point(size=6, alpha =0.65) +
#   geom_hline(yintercept = 0, linetype = "dashed", alpha =0.25) +
#   geom_vline(xintercept = 0, linetype = "dashed",  alpha =0.25) +
#   geom_text(aes(label = PairsID), hjust = -0.4, vjust = -0.3) +
#   labs(x = paste("PC1 (", round(variance_explained[1] * 100, 2), "%)"),
#        y = paste("PC2 (", round(variance_explained[2] * 100, 2), "%)")) +
#   ggtitle("WHO_2_3_CE")+
#   theme_minimal() + 
#   guides(color = FALSE, shape = FALSE)
# 
# biplot_ggplot_5
# 
# biplot_ggplot_6 <- ggplot(pca_df, aes(PC1, PC2, color = CE_cats$WHO_3_3_CE)) + #shape = factor(phased_haplos_cats$final_cats)
#   geom_point(size=6, alpha =0.65) +
#   geom_hline(yintercept = 0, linetype = "dashed", alpha =0.25) +
#   geom_vline(xintercept = 0, linetype = "dashed",  alpha =0.25) +
#   geom_text(aes(label = PairsID), hjust = -0.4, vjust = -0.3) +
#   labs(x = paste("PC1 (", round(variance_explained[1] * 100, 2), "%)"),
#        y = paste("PC2 (", round(variance_explained[2] * 100, 2), "%)")) +
#   ggtitle("WHO_3_3_CE")+
#   theme_minimal() + 
#   guides(color = FALSE, shape = FALSE)
# 
# biplot_ggplot_6
# 
# panels_plots <- grid.arrange(biplot_ggplot_3, biplot_ggplot_4, biplot_ggplot_6, biplot_ggplot_5, ncol = 2)
# ggsave("PCA_NGS-CE_final.png", panels_plots, width = 16, height = 9, bg = "white")
# 
# ####################################################################################################################3
# 
# #how many markers to use? 
# #test changes in NI and R by removing markers
# mean_HE<- mean_HE[order(-mean_HE$mean_HE), ]
# 
# #list of dfs removing the last row iteratively
# df_list <- list()
# df_list[[1]] <- mean_HE
# n_rows <- nrow(mean_HE)
# 
# for (i in 1:(n_rows - 1)) {
#   df_list[[i + 1]] <- mean_HE[1:(n_rows - i), ]
# }
# 
# df_list <- rev(df_list)
# 
# #benchmark
# var_amps <-c()
# infection_results_list <- list()
# 
# for (amps in df_list){
#   var_amps <- merged_dfs_filtered[merged_dfs_filtered$locus %in% amps$locus, ]
#   
#   common_alleles_table <-var_amps %>%
#     filter(time_point %in% c("D0", "Dx")) %>%
#     group_by(PairsID, locus) %>%
#     summarise(alleles_D0 = length(unique(pseudo_cigar[time_point == "Dx"])),
#               alleles_Dx = length(unique(pseudo_cigar[time_point == "Dx"])),
#               common_alleles_D0_Dx = n_distinct(pseudo_cigar[pseudo_cigar %in% intersect(unique(pseudo_cigar[time_point == "D0"]), unique(pseudo_cigar[time_point == "Dx"]))]))
#   
#   #criteria:
#   #what does each marker say? (using NI for "New Infection" == "Reinfection")
#   #Classification of recrudescence (R) required that at least one of the haplotypes occurred in both, the pre- and post-treatment sample, with a minimum haplotype frequency of 1% in at least two of the three independent replicates performed. 
#   #A new infection (NI) was defined by the occurrence of only new haplotypes in the post-treatment sample with a minimum haplotype frequency of 1%.
#   common_alleles_table$infection <- ifelse(common_alleles_table$common_alleles_D0_Dx > 0, "R", "NI")
#   
#   #what does each marker say?
#   infection_results <- common_alleles_table %>%
#     group_by(PairsID) %>%
#     summarise(Recrudescence_markers = sum(infection == "R"),
#               Reinfection_markers = sum(infection == "NI"),
#               perc_Recrudescence_markers = Recrudescence_markers/(Recrudescence_markers + Reinfection_markers),
#               perc_Reinfection_markers = Reinfection_markers/(Recrudescence_markers + Reinfection_markers))
#   
#   # add phasing data
#   infection_results$phasing_dhfr_dhps_haplo <- phased_haplos_cats$final_cats
#   
#   # in case markers disagree...
#   # (i) The WHO/MMV recommended approach1, where a recurrent parasitemia is classified in the overall outcome as NI, if at least one of the three markers had given a NI result
#   infection_results$WHO_MMV_cireteria <- ifelse(infection_results$Reinfection_markers > 0, "NI", "R")
#   
#   # (ii) The alternative, newer approach termed “2/3 algorithm”8. This algorithm uses the consensus of two markers as final PCR-correction outcome.
#   infection_results$two_of_three_algo <- ifelse(infection_results$perc_Reinfection_markers >= (2/3), "NI", 
#                                                 ifelse(infection_results$perc_Recrudescence_markers >= (2/3), "R", "?"))
#   
#   infection_results$two_three_phasing <- infection_results$two_of_three_algo 
#   
#   infection_results$two_three_phasing <- ifelse(infection_results$two_of_three_algo == "?" & infection_results$phasing_dhfr_dhps_haplo == "non-overlapping_diversity_same", "NI", infection_results$two_three_phasing) 
# 
#   #two_three_majority_rule_phasing algo:
#   infection_results$two_three_majority_rule_phasing <- infection_results$two_three_phasing 
#   
#   infection_results$two_three_majority_rule_phasing <- ifelse(infection_results$two_three_majority_rule_phasing == "?",
#                                                               ifelse(infection_results$perc_Recrudescence_markers > infection_results$perc_Reinfection_markers, "R","NI"),
#                                                               infection_results$two_three_majority_rule_phasing)
#   
#   #two_three_majority_rule_phasing algo:
#   infection_results$two_three_phasing_WHO_MMV <- infection_results$two_three_phasing 
#   
#   infection_results$two_three_phasing_WHO_MMV <- ifelse(infection_results$two_three_phasing == "?", infection_results$WHO_MMV_cireteria, infection_results$two_three_majority_rule_phasing )
#   
#   infection_results_list[[length(infection_results_list) + 1]] <- infection_results
# }
# 
# 
# #create criteria lists and plots
# df_x <- data.frame(row.names = c(1, 2))
# infections_results_ALL_iterations_WHO_MMV <-data.frame() 
# 
# for (i in seq_along(infection_results_list)) {
#   
#   current_results <- infection_results_list[[i]]
#   df_x$infection <- names(table(current_results$WHO_MMV_cireteria))
#   df_x$samples <- as.numeric(table(current_results$WHO_MMV_cireteria))
#   df_x$num_amplicons <- i
#   
#   infections_results_ALL_iterations_WHO_MMV <- rbind(infections_results_ALL_iterations_WHO_MMV, df_x)
#   
# }
# 
# a <- ggplot(infections_results_ALL_iterations_WHO_MMV, aes(x = as.factor(num_amplicons), y = samples, fill = infection)) +
#   geom_bar(stat = "identity") +
#   labs(x = "Amplicons", y = "Samples", fill = "Infection") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#         panel.grid.major.x = element_blank(),
#         panel.grid.minor.x = element_blank())+
#   ggtitle("WHO_MMV")
# 
# 
# df_x <- data.frame(row.names = c(1, 2))
# infections_results_ALL_iterations_two_three_phasing_WHO_MMV <-data.frame() 
# 
# for (i in seq_along(infection_results_list)) {
#   
#   current_results <- infection_results_list[[i]]
#   df_x$infection <- names(table(current_results$two_three_phasing_WHO_MMV))
#   df_x$samples <- as.numeric(table(current_results$two_three_phasing_WHO_MMV))
#   df_x$num_amplicons <- i
#   
#   infections_results_ALL_iterations_two_three_phasing_WHO_MMV <- rbind(infections_results_ALL_iterations_two_three_phasing_WHO_MMV, df_x)
#   
# }
# 
# b <- ggplot(infections_results_ALL_iterations_two_three_phasing_WHO_MMV, aes(x = as.factor(num_amplicons), y = samples, fill = infection)) +
#   geom_bar(stat = "identity") +
#   labs(x = "Amplicons", y = "Samples", fill = "Infection") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#         panel.grid.major.x = element_blank(),
#         panel.grid.minor.x = element_blank())+
#   ggtitle("two_three_phasing_WHO_MMV")
# 
# ## el amplicon solo es el que tiene mayor He
# amplicon_plots <- grid.arrange(a, b, nrow = 2)
# ggsave("amplicon_stacked_bar.png", amplicon_plots, width = 24, height = 12, bg = "white")
# ggsave("amplicon_stacked_bar_WHO_MMV.png", a, width = 20, height = 8, bg = "white")
