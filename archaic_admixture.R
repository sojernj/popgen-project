## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, fig.width = 14, fig.height = 5)
library(tidyverse)
library(cowplot)
library(knitr)
theme_set(theme_bw(base_size = 12))
#setwd("~/Documents/git/popgen-project/")


## ----data-initialization------------------------------------------------------

archaic_segments <- read_tsv("ArchaicAdmixture/ArchaicSegments.txt") %>% 
  mutate(length = end-start+1000) %>% 
  filter(length>0)

test_segments <- archaic_segments %>% slice_sample(n=10000) %>% mutate(chrom=factor(chrom, levels= c(1:22, 'X')))

# snps <- read_tsv("ArchaicAdmixture/SNPs.txt") %>%  
#     mutate(origin = factor(
#       case_when(
#         archaic %in% c("AltaiNea_nofil","Vi33.19_nofil") ~ "Neanderthal",
#         archaic == "Denisova_nofil" ~ "Denisova",
#         .default = 'human'),
#       levels = c('human', 'Neanderthal', 'Denisova')))

# test_snps <- snps %>% slice_sample(n=10000)
# test_snps %>% head() %>% kable
# length(unique(test_snps$start))

# archaic_snps <- snps %>% filter(archaic != "human") %>% select(-c("archaic"))
# nrow(archaic_snps)
# nrow(archaic_segments)


## ---- eval=FALSE-----------------------------------------------------------
## test_segments %>% filter(name %in% c(individual_1,individual_2)) %>% group_by(chrom, name)



## ----function-find-shared-segments--------------------------------------------

# Function to find shared segments between two individuals
find_shared_segments <- function(individual_1, individual_2, segments_data) {
  
  df <- segments_data %>% filter(name %in% c(individual_1,individual_2)) %>% 
    mutate(origin = case_when(
      Shared_with_Altai > Shared_with_Denisova | Shared_with_Vindija > Shared_with_Denisova ~ "Neanderthal",
      Shared_with_Altai < Shared_with_Denisova & Shared_with_Vindija < Shared_with_Denisova ~ "Denisovan",
      Shared_with_Altai + Shared_with_Vindija + Shared_with_Denisova == 0 ~ "Unclassified"
    )) %>% select(name, chrom, start, end, origin)
  df1 <- df %>% filter(name==individual_1)
  df2 <- df %>% filter(name==individual_2)
  
  #cat("Comparing individuals: ", individual_1, ":", individual_2)
  
  shared_segments <- data.frame()
  
  for (chrom in c(1:22, 'X')) {
    #cat("Chromosome:", chrom, "\t")
    
    df1_chrom <- df1 %>% filter(chrom == !!chrom)
    df2_chrom <- df2 %>% filter(chrom == !!chrom)
    
    overlaps <- fuzzyjoin::interval_inner_join(df1_chrom, df2_chrom, c("start", "end")) %>% 
      filter((origin.x == origin.y) & (origin.x %in% c("Neanderthal", "Denisova")))
    
    overlaps <- overlaps %>% mutate(
      overlap_length = case_when(nrow(overlaps)>0 ~ pmin(end.x, end.y) - pmax(start.x, start.y),
                                 .default = 0))
    
    
    #cat("--> length:", sum(overlaps$overlap_length), "\n")
    
    shared_segments <- rbind(shared_segments, overlaps)
  }
  
  shared_segments <- shared_segments %>% mutate(chrom=chrom.x, origin=origin.x) %>% select(name.x, name.y, chrom, start.x, start.y, end.x, end.y, origin, overlap_length)
  
  #cat("\t total overlap: ", sum(shared_segments$overlap_length), "\n")
  
  return(shared_segments )
}

# Test function 
#  find_shared_segments("S_Hungarian-1", "ERS700536", test_segments)
 
# # Test validity 
# find_shared_segments("ERS700536", "ERS700536", test_segments) %>% 
#     summarise(total=sum(overlap_length))


# test_segments %>% filter(name=="ERS700536") %>% summarise(total=sum(length))



## ----use-function, eval=FALSE-------------------------------------------------
## 
## # Example: Calculate shared SNPs between two individuals
## individual_1 <- sample(unique(test_segments$name), size = 1)
## individual_2 <- sample(unique(test_segments$name), size = 1)
## shared_segments <- find_shared_segments(individual_1, individual_2, test_segments)
## shared_segments
## 


## ----use-the-cluster----------------------------------------------------------
cat("Starting the multiprocessing\n")
# Get all unique individuals
individuals <- unique(test_segments$name)

# Define a function to be applied to each pair
find_shared_segments_for_pair <- function(pair) {
  shared_segments <- find_shared_segments(pair[1], pair[2], archaic_segments)
  
  # Calculate the total length of overlapping segments
  total_length <- sum(shared_segments$overlap_length)
  
  # Return a data frame with the results
  return(data.frame(individual_1 = pair[1], 
                    individual_2 = pair[2], 
                    total_length = total_length))
}

# Generate all possible pairs of unique individuals and apply the function
library(future)
library(furrr)

pairs <- combn(individuals, 2, simplify = TRUE)

plan(cluster)
t0 <- Sys.time()
results <- 1:ncol(pairs) %>%
  future_map_dfr(\(i) find_shared_segments_for_pair(pairs[,i]), .progress = TRUE, .options = furrr_options(stdout = TRUE))

t1 <- Sys.time()
# cat("took", t1-t0, "sec")

# results <- combn(individuals, 2, FUN = find_shared_segments_for_pair, simplify = FALSE)

# # Combine the results into a single data frame
# results <- do.call(rbind, results)

# # Write the results to a CSV file
 write.csv(results, file = "results_diff.csv", row.names = FALSE)



## ----plot-the-counted-snps-derived-from-each-archaic-ancestry, fig.width=8 , eval=FALSE----
## 
## plot_admix_prop <- snps %>%
##   group_by(name, region, origin) %>%
##   summarize(shared_snps = n()) %>%
##   ungroup() %>%
##   group_by(name) %>%
##   reframe(region,
##           origin,
##           shared_snps,
##           shared_prop = shared_snps#/sum(shared_snps)
##             ) %>%
##   arrange(origin) %>%
##   ggplot(aes(x=reorder(name, shared_snps), y=shared_prop, fill=origin)) +
##   geom_bar(stat='identity', position="fill", width=1) +
##   scale_x_discrete(expand=c(0,0), name = "Individuals") +
##   scale_y_continuous(expand=c(0,0), name = "Shared ancestry") +
##   theme(
##     axis.text.x = element_blank(),
##     axis.ticks.x = element_blank()
##   )
## plot_admix_prop
## 


## ----save-whole-plot, eval=FALSE----------------------------------------------
## ggsave(plot = plot_admix_prop, device = 'svg',
##        path = "graphics/",
##        filename = "admixture_proportion_whole.svg", width = 6.5, height = 2,
##        dpi=2000)


## ----table-mean-admixture-proportions, eval=FALSE-----------------------------
## mean_admix <- snps %>%
##   group_by(name, region, origin) %>%
##   summarize(shared_snps = n()) %>%
##   ungroup() %>%
##   group_by(name) %>%
##   reframe(region,
##           origin,
##           shared_snps,
##           shared_prop = shared_snps/sum(shared_snps)
##           ) %>%
##   group_by(origin) %>%
##   summarise(full = mean(shared_prop)) %>%
##   pivot_longer(cols = 2, names_to="region",values_to = "mean_prop")
## 
## mean_admix
## 
## mean_admix_by_region <- snps %>%
##   group_by(name, region, origin) %>%
##   summarize(shared_snps = n()) %>%
##   ungroup() %>%
##   group_by(name) %>%
##   reframe(region,
##           origin,
##           shared_snps,
##           shared_prop = shared_snps/sum(shared_snps)
##           ) %>%
##   group_by(origin,region) %>%
##   summarise(mean_prop = mean(shared_prop))
## 
## mean_admix <- bind_rows(mean_admix, mean_admix_by_region)
## 
## mean_admix %>% pivot_wider(names_from=region, values_from = mean_prop) %>% kable(format = "latex", digits=3)
## 
## # I don't know of this is relevant, but here goes ANOVA one way
## aov(mean_prop~region, data=mean_admix_by_region) %>% summary(.)
## aov(mean_prop~region+origin, data=mean_admix_by_region) %>% summary(.)
## 


## ----plot the counted snps derived from each archaic ancestry per region, eval=FALSE----
## 
## # Count the number of individuals in each facet group
## facet_counts <- snps %>%
##   group_by(region) %>%
##   summarize(num_individuals = n_distinct(name))
## 
## # Make a labelling to be called in facet_wrap
## include_n <- function(name) {
##   number <- facet_counts$num_individuals[facet_counts$region == name]
##   paste0(name, " (", number, ")")
## }
## 
## plot_admix_prop_facet <- plot_admix_prop +
##   facet_wrap(~region, scales='free_x', nrow = 1,
##              labeller = as_labeller(include_n)) +
##   theme(panel.spacing.x = unit(1, "mm"))
## plot_admix_prop_facet
## 


## ----save-facetted-plot, eval=FALSE-------------------------------------------
## ggsave(plot = plot_admix_prop_facet, device = 'svg',
##        path = "graphics/",
##        filename = "admixture_proportion_byorigin.svg", width = 6.5, height = 2,
##        dpi=2000)
## 


## ----save the combined plot, fig.height=10, eval=FALSE------------------------
## 
## admix_proportion_combined <- plot_grid(plot_admix_prop,
##                                        plot_admix_prop_facet +  theme(legend.position = 'none'),
##                                        labels = c('A','B'), nrow = 2)
## admix_proportion_combined
## 


## ----save-combined-plot, eval=FALSE-------------------------------------------
## ggsave(plot = plot_admix_prop_facet, device = 'svg',
##        path = "graphics/",
##        filename = "admixture_proportion_byorigin.svg", width = 6.5, height = 2,
##        dpi=2000)


## ----count-how-many-times-each-snp-is-shared-in-an-individual, eval=FALSE-----
## 
## shared <- snps %>%
##   filter(origin != "human") %>%
##   group_by(chrom, start, origin) %>%
##   summarise(n_shared = n()) %>%
##   mutate(chrom=factor(chrom, levels=c(1:22,'X')))
## 
## head(shared)
## summary(shared$n_shared) %>% as.array() %>% as.data.frame() %>% pivot_wider(names_from=Var1, values_from = Freq) %>% kable(format='latex')
## 
## # Make a density plot
## plot_snp_counts <- shared %>%
##   ggplot(aes(x=n_shared, fill=origin)) +
##   geom_bar(stat='density', position='dodge', width=1, color=NA) +
##   #geom_density(aes(alpha=origin, color=origin)) +
##   #ggtitle("Distribution of shared SNP positions based on inferred origin") +
##   xlab("Number of times a SNP position is repeated") +
##   ylab("Density") +
##   coord_cartesian(xlim=c(1,150)) +
##   theme(legend.position = c(0.7,0.7),
##         legend.margin = margin(.1,.1,.1,.1),
##         legend.box.margin = margin(.1,.1,.1,.1),
##         legend.background = element_blank()) +
##   NULL
## plot_snp_counts


## ----save-plot-snp-counts, eval=FALSE-----------------------------------------
## ggsave(plot_snp_counts,
##        path = "graphics/",
##        filename = "snp-counts-density-1_150.svg", width = 6.5, height = 2)


## ----boxplot-snp-counts, eval=FALSE-------------------------------------------
## 
## # Prepare for plotting jittered outliers on a boxplot
## shared <- shared %>%
##   group_by(origin) %>%
##   mutate(outlier.high = n_shared > quantile(n_shared, .75) + 1.50*IQR(n_shared),
##          outlier.low = n_shared < quantile(n_shared, .25) - 1.50*IQR(n_shared)) %>%
##   mutate(outlier.color = case_when(outlier.high ~ "firebrick4",
##                                    outlier.low ~"steelblue"))
## 
## # Make a boxplot (also not really worth showing)
## shared %>%
##   ggplot(aes(x=origin, y=n_shared, fill=origin)) +
##   #geom_boxplot(outlier.shape=NA) +
##   #geom_jitter(color=shared$outlier.color, width=0.35, alpha=0.1) +
##   geom_jitter(aes(color=origin), alpha=0.1)
## 
## # Calculate the mean shared number of snps across origin
## shared %>%
##   group_by(origin) %>%
##   summarise(min = min(n_shared),
##             mean = mean(n_shared),
##             median = median(n_shared),
##             max = max(n_shared),
##             n = n(),
##             sum(n_shared <= 5)/n,
##             sum(n_shared <= 10)/n,
##             sum(n_shared <= 50)/n,
##             prop_above_50   = sum((n_shared > 50))/n) %>%
##   kable(format = 'latex', digits=2)


## ---- eval=FALSE--------------------------------------------------------------
## rnorm()


## ----total-archaic-segments, eval=FALSE---------------------------------------
## head(archaic_segments)
## mean(archaic_segments$length)
## summary(archaic_segments$length)
## 
## total_seg_ind <- archaic_segments %>%
##   group_by(name, pop, region) %>%
##   summarise(total_length = sum(length)) %>%
##   arrange(region)
## 
## summary(total_seg_ind$total_length)
## 
## plot_total_seg_ind <- total_seg_ind %>%
##   ggplot(aes(x=total_length, fill=region)) +
##   geom_histogram(bins=50, color='black', linewidth=0.2) +
##   geom_vline(aes(xintercept = median(total_seg_ind$total_length)), linewidth=.8, color='black') +
##   annotate('text', x=median(total_seg_ind$total_length)*1.2, y=58, label="Median") +
##   xlab("Total length") +
##   ylab("Count") +
##   scale_x_continuous(expand=c(0,0)) +
##   scale_y_continuous(expand=c(0,0)) +
##   coord_cartesian(ylim = c(0,60))
## plot_total_seg_ind


## ----save-total-archaic-segments, eval=FALSE----------------------------------
## ggsave(plot_total_seg_ind, path="graphics/", filename="total_segment_distribution.svg",
##        width = 6.5, height = 2)
## 


## ----skip-the-rest, eval=FALSE------------------------------------------------
## knit_exit()


## ---- eval=FALSE--------------------------------------------------------------
## 
## population_comparison <- archaic_segments %>%
##   # Perform analysis here
##   # ...
## 
## # Visualize population comparison results
## population_plot <- ggplot(data = population_comparison, aes(x = ..., y = ..., color = ...)) +
##   geom_point() +
##   theme_minimal() +
##   labs(title = "Population Comparison")
## 
## # Print population comparison plot
## print(population_plot)
## 


## ---- eval=FALSE--------------------------------------------------------------
## region_comparison <- archaic_segments %>%
##   # Perform analysis here
## 
## # Visualize region comparison results
## region_plot <- ggplot(data = region_comparison, aes(x = ..., y = ..., color = ...)) +
##   geom_point() +
##   theme_minimal() +
##   labs(title = "Geographical Region Comparison")
## 
## # Print region comparison plot
## print(region_plot)
## 


## ---- eval=FALSE--------------------------------------------------------------
## epas1_analysis <- archaic_segments %>%
##   # Perform analysis here
## 
## # Visualize EPAS1 analysis results
## epas1_plot <- ggplot(data = epas1_analysis, aes(x = ..., y = ..., color = ...)) +
##   geom_point() +
##   theme_minimal() +
##   labs(title = "EPAS1 Gene Analysis")
## 
## # Print EPAS1 analysis plot
## print(epas1_plot)
## 


## ---- eval=FALSE--------------------------------------------------------------
## admixture_total <- archaic_segments %>%
##   # Perform analysis here
## 
## # Visualize total admixture results
## admixture_plot <- ggplot(data = admixture_total, aes(x = ..., y = ..., color = ...)) +
##   geom_point() +
##   theme_minimal() +
##   labs(title = "Total Admixture Calculation")
## 
## # Print total admixture plot
## print(admixture_plot)
## 


## ---- eval=FALSE--------------------------------------------------------------
## 
## correlation_admixture_analysis <- archaic_segments %>%
##   # Perform analysis here
## 
## # Visualize correlation and admixture total results
## correlation_admixture_plot <- ggplot(data = correlation_admixture_analysis, aes(x = ..., y = ..., color = ...)) +
##   geom_point() +
##   theme_minimal() +
##   labs(title = "Correlation and Admixture Total Analysis")
## 
## # Print correlation and admixture total plot
## print(correlation_admixture_plot)
## 


## ---- eval=FALSE--------------------------------------------------------------
## additional_analyses <- archaic_segments %>%
##   # Perform analysis here
## 
## # Visualize additional analysis results
## additional_plot <- ggplot(data = additional_analyses, aes(x = ..., y = ..., color = ...)) +
##   geom_point() +
##   theme_minimal() +
##   labs(title = "Additional Analyses")
## 
## # Print additional analysis plot
## print(additional_plot)
## 


## ---- eval=FALSE--------------------------------------------------------------
## 

