# Load necessary libraries
library(vegan)  # For Morisita-Horn dissimilarity
library(ggplot2)  # For plotting
library(dplyr)  # For data manipulation
library(tidyr)

# Load dataset
data <- read.csv("bee_sex_plant_morisita_horn.csv")


# Function to calculate Morisita-Horn index
calculate_morisita_horn <- function(data) {
  
  # Create plant species visit matrix
  visit_matrix <- data %>%
    group_by(bee_species, bee_sex, plant_species) %>%
    summarise(visits = n(), .groups = "drop") %>%
    pivot_wider(names_from = plant_species, values_from = visits, values_fill = 0)
  
  # Set result variables 
  species <- unique(visit_matrix$bee_species)
  results <- data.frame(Bee_Species = species, Observed_MH = NA)
  
  # Calculate the index and save it to results
  for (bee in species) {
    subset_data <- filter(visit_matrix, bee_species == bee)
    if (n_distinct(subset_data$bee_sex) == 2) {
      male_vector <- subset_data %>% filter(bee_sex == "M") %>% select(-bee_species, -bee_sex) %>% as.matrix()
      female_vector <- subset_data %>% filter(bee_sex == "F") %>% select(-bee_species, -bee_sex) %>% as.matrix()
      results$Observed_MH[results$Bee_Species == bee] <- as.numeric(vegdist(rbind(male_vector, female_vector), method = "horn"))
    }
  }
  return(results)
}


# Function to calculate Morisita-Horn index for the Null model with shuffling
calculate_null_morisita_horn <- function(data, n_perm = 999, seed = 123) {
  
  # Set seed ONCE for reproducibility across iterations
  set.seed(seed)
  
  species <- unique(data$bee_species)
  null_results_list <- vector("list", n_perm)  # Store all iterations
  
  # Iterate the calculation of the index (shuffled data every time)
  for (i in 1:n_perm) {
  

    # data_shuffled <- data %>%
    #   mutate(bee_sex = sample(bee_sex))
    
    # Shuffle the bee_sex column while maintaining M/F counts per plant_species
    data_shuffled <- data %>%
      group_by(bee_species) %>%
      mutate(bee_sex = sample(bee_sex)) %>%  # Shuffle within plant species
      ungroup()
    
    # Create plant species visit matrix
    visit_matrix <- data_shuffled %>%
      group_by(bee_species, bee_sex, plant_species) %>%
      summarise(visits = n(), .groups = "drop") %>%
      pivot_wider(names_from = plant_species, values_from = visits, values_fill = 0)
    
    # Create dataframe to store results per species
    iteration_results <- data.frame(Bee_Species = species, Null_MH = NA)
    
    # Calculate the index
    for (bee in species) {
      subset_data <- filter(visit_matrix, bee_species == bee)
      if (n_distinct(subset_data$bee_sex) == 2) {
        male_vector <- subset_data %>% filter(bee_sex == "M") %>% select(-bee_species, -bee_sex) %>% as.matrix()
        female_vector <- subset_data %>% filter(bee_sex == "F") %>% select(-bee_species, -bee_sex) %>% as.matrix()
        iteration_results$Null_MH[iteration_results$Bee_Species == bee] <- as.numeric(vegdist(rbind(male_vector, female_vector), method = "horn"))
      }
    }
    
    null_results_list[[i]] <- iteration_results  # Store each iteration
  }
  
  # Combine all iterations into a single dataframe
  null_results_df <- bind_rows(null_results_list)
  
  # Calculate Mean, SD, 95% CI, and Normality Tests per species
  final_results <- null_results_df %>%
    group_by(Bee_Species) %>%
    summarise(
      Null_MH_Mean = mean(Null_MH, na.rm = TRUE),
      simSD = sd(Null_MH, na.rm = TRUE),  # Added Standard Deviation Calculation
      CI_95_Low = quantile(Null_MH, 0.025, na.rm = TRUE),
      CI_95_Up = quantile(Null_MH, 0.975, na.rm = TRUE),
      
      # Kolmogorov-Smirnov Test for Normality
      ks_p = ks.test(Null_MH, "pnorm", mean = Null_MH_Mean, sd = simSD)$p.value,
      
      # Shapiro-Wilk Test for Normality (Direct Computation)
      shapiro_p = shapiro.test(Null_MH)$p.value,
      
      # Empirical P-Value Calculation
      empirical_p = (sum(Null_MH >= Null_MH_Mean, na.rm = TRUE) + 1) / (n() + 1),
      
      # Assign significance stars based on empirical p-values
      significance_empirical = case_when(
        empirical_p < 0.001 ~ "***",  # Highly significant
        empirical_p < 0.01  ~ "**",   # Very significant
        empirical_p < 0.05  ~ "*",    # Significant
        TRUE               ~ "ns"     # Not significant
      ),
      
      .groups = "drop"
    )
  
  return(final_results)
}

# Compute observed Morisita-Horn index
observed_results <- calculate_morisita_horn(data)

# Compute null model
null_results <- calculate_null_morisita_horn(data, n_perm = 999, seed = 42)

# Merge observed and null results
summary_model <- merge(observed_results, null_results, by = "Bee_Species")

# Save results
write.csv(summary_model, "Morisita_Horn_Null_Model_Results.csv", row.names = FALSE)

# Compute Z-score and p-value
summary_model <- summary_model %>%
  mutate(
    z_score = (Observed_MH - Null_MH_Mean) / simSD,  # Compute Z-score
    p_value = 2 * (1 - pnorm(abs(z_score))),  # Two-tailed p-value from Z-score
    significance = case_when(
      p_value < 0.001 ~ "***",   # Highly significant
      p_value < 0.01  ~ "**",    # Very significant
      p_value < 0.05  ~ "*",     # Significant
      TRUE            ~ "ns"     # Not significant
    )
  )

# Print the updated summary_model with p-values
print(summary_model)

#Plot version 3

ggplot(summary_model, aes(x = reorder(Bee_Species, -Observed_MH))) +  
  geom_errorbar(aes(ymin = CI_95_Low, ymax = CI_95_Up), width = 0.2, color = "black") +  
  geom_point(aes(y = Null_MH_Mean), color = "blue", size = 2) +  
  geom_point(aes(y = Observed_MH), color = "red", size = 3) + 
  geom_text(aes(y = Observed_MH + 0.05, label = significance), size = 5, fontface = "bold") +  # Add significance labels
  labs(title = "Observed vs Null Morisita-Horn Index",  
       x = "Bee Species",  
       y = "Morisita-Horn Dissimilarity Index") +  
  theme_minimal() +  
  theme(
    text = element_text(size = 16),
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1, face = "italic"),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
  )


#Normality test for null model distributions-----

#Shapiro-Wilk normality test
shapiro_test <- shapiro.test(null_results$Null_MH_Mean)
print(shapiro_test)

# Compute normality tests per bee species
normality_results <- null_results %>%
  group_by(Bee_Species) %>%
  summarise(
    shapiro_p = ifelse(n() >= 3 & n() <= 5000, shapiro.test(Null_MH)$p.value, NA),  # Shapiro-Wilk Test
    ks_p = ks.test(Null_MH, "pnorm", mean(Null_MH, na.rm = TRUE), sd(Null_MH, na.rm = TRUE))$p.value  # KS Test
  ) %>%
  ungroup()

# Print results
print(normality_results)


#Kolmogorov-Smirnov test
ks_test <- ks.test(null_results$Null_MH_Mean, "pnorm", mean(null_results$Null_MH_Mean, na.rm = TRUE), sd(null_results$Null_MH_Mean, na.rm = TRUE))
print(ks_test)

# Histogram with density curve
ggplot(null_results, aes(x = Null_MH_Mean)) + 
  geom_histogram(binwidth = 0.05, color = "black", fill = "lightblue", alpha = 0.7) + 
  geom_density(aes(y = ..density.. * length(null_results) * 0.05), color = "red", size = 1) +  
  labs(title = "Distribution of Null Morisita-Horn Indices",
       x = "Null Morisita-Horn Index",
       y = "Frequency") +  
  theme_minimal()

# QQ plot
qqnorm(null_results$Null_MH_Mean)
qqline(null_results$Null_MH_Mean, col = "red")

# Compute empirical p-value for non-parametric data
summary_model <- summary_model %>%
  rowwise() %>%
  mutate(
    p_value_empirical = (sum(null_results$Null_MH >= Observed_MH, na.rm = TRUE) + 1) / (nrow(null_results) + 1),
    significance = case_when(
      p_value_empirical < 0.001 ~ "***",  # Highly significant
      p_value_empirical < 0.01  ~ "**",   # Very significant
      p_value_empirical < 0.05  ~ "*",    # Significant
      TRUE                      ~ "ns"    # Not significant
    )
  ) %>%
  ungroup()

#add the empirical p-value significance to the plot:

ggplot(summary_model, aes(x = reorder(Bee_Species, -Observed_MH))) +  
  geom_errorbar(aes(ymin = CI_95_Low, ymax = CI_95_Up), width = 0.2, color = "black") +  
  geom_point(aes(y = Null_MH_Mean), color = "blue", size = 2) +  
  geom_point(aes(y = Observed_MH), color = "red", size = 3) + 
  geom_text(aes(y = Observed_MH + 0.05, label = significance), size = 5, fontface = "bold") +  # Add significance labels
  labs(title = "Observed vs Null Morisita-Horn Index",  
       x = "Bee Species",  
       y = "Morisita-Horn Dissimilarity Index") +  
  theme_minimal() +  
  theme(
    text = element_text(size = 16),
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1, face = "italic"),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
  )

#Compute a Rank-Based Effect Size (Mann-Whitney U Alternative)

summary_model <- summary_model %>%
  rowwise() %>%
  mutate(
    rank_based_z = (rank(c(Observed_MH, null_results$Null_MH))[1] - median(rank(null_results$Null_MH))) /
      sd(rank(null_results$Null_MH), na.rm = TRUE)
  ) %>%
  ungroup()

#plot the rank-based Z-score

ggplot(summary_model, aes(x = reorder(Bee_Species, -rank_based_z), y = rank_based_z)) +  
  geom_bar(stat = "identity", fill = "steelblue", color = "black", alpha = 0.7) +  
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1) +  # Reference line at 0
  labs(title = "Rank-Based Z-Scores of Morisita-Horn Index",  
       x = "Bee Species",  
       y = "Rank-Based Z-Score") +  
  theme_minimal() +  
  theme(
    text = element_text(size = 16),
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1, face = "italic"),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
  )


# Ensure necessary libraries are loaded
library(dplyr)

# Compute Wilcoxon Signed-Rank Test per bee species
summary_model <- summary_model %>%
  rowwise() %>%
  mutate(
    wilcox_p = wilcox.test(
      x = null_results %>% filter(Bee_Species == Bee_Species) %>% pull(Null_MH_Mean),  # Extracting the null model mean values
      mu = Observed_MH,  # Comparing against the observed Morisita-Horn index
      alternative = "two.sided"  # Two-tailed test
    )$p.value,
    
    # Assign significance stars based on Wilcoxon p-values
    significance_wilcox = case_when(
      wilcox_p < 0.001 ~ "***",  # Highly significant
      wilcox_p < 0.01  ~ "**",   # Very significant
      wilcox_p < 0.05  ~ "*",    # Significant
      TRUE             ~ "ns"     # Not significant
    )
  ) %>%
  ungroup()

# Print results with Wilcoxon p-values
print(summary_model)


ggplot(summary_model, aes(x = reorder(Bee_Species, -Observed_MH))) +  
  geom_errorbar(aes(ymin = CI_95_Low, ymax = CI_95_Up), width = 0.2, color = "black") +  
  geom_point(aes(y = Null_MH_Mean), color = "blue", size = 2) +  
  geom_point(aes(y = Observed_MH), color = "red", size = 3) + 
  geom_text(aes(y = Observed_MH + 0.05, label = significance_wilcox), size = 5, fontface = "bold") +  # Wilcoxon significance labels
  labs(title = "Wilcoxon Test: Observed vs Null Morisita-Horn Index",  
       x = "Bee Species",  
       y = "Morisita-Horn Dissimilarity Index") +  
  theme_minimal() +  
  theme(
    text = element_text(size = 16),
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1, face = "italic"),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
  )


# Compute empirical p-value per bee species
summary_model <- summary_model %>%
  rowwise() %>%
  mutate(
    empirical_p = (sum(null_results %>% filter(Bee_Species == Bee_Species) %>% pull(Null_MH_Mean) >= Observed_MH, na.rm = TRUE) + 1) / 
      (nrow(null_results %>% filter(Bee_Species == Bee_Species)) + 1),
    
    # Assign significance stars based on empirical p-values
    significance_empirical = case_when(
      empirical_p < 0.001 ~ "***",  # Highly significant
      empirical_p < 0.01  ~ "**",   # Very significant
      empirical_p < 0.05  ~ "*",    # Significant
      TRUE               ~ "ns"     # Not significant
    )
  ) %>%
  ungroup()

# Print results with empirical p-values
print(summary_model)

ggplot(summary_model, aes(x = reorder(Bee_Species, -Observed_MH))) +  
  geom_errorbar(aes(ymin = CI_95_Low, ymax = CI_95_Up), width = 0.2, color = "black") +  
  geom_point(aes(y = Null_MH_Mean), color = "blue", size = 2) +  
  geom_point(aes(y = Observed_MH), color = "red", size = 3) + 
  geom_text(aes(y = Observed_MH + 0.05, label = significance_empirical), size = 5, fontface = "bold") +  # Empirical significance labels
  labs(title = "Empirical P-Value Test: Observed vs Null Morisita-Horn Index",  
       x = "Bee Species",  
       y = "Morisita-Horn Dissimilarity Index") +  
  theme_minimal() +  
  theme(
    text = element_text(size = 16),
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1, face = "italic"),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
  )
