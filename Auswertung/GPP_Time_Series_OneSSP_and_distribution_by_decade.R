# Load required libraries
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr"))   install.packages("dplyr")
if (!require("zoo"))     install.packages("zoo")
if (!require("ggpubr"))  install.packages("ggpubr")

library(ggplot2)
library(dplyr)
library(zoo)
library(ggpubr)

# Read CSV data
setwd("C:/Users/philipp/Documents/TropForLPJmL_outputs/")
globalflux_old <- read.csv("globalflux_oldLU126_neu.csv", sep = ",", skip = 2)
globalflux_new <- read.csv("globalflux_newLU126_neu.csv", sep = ",", skip = 2)

# Define column names
columns <- c(
  "Year", "NEP", "GPP", "NPP", "RH", "estabc", "negc_fluxes", "firec", "area",
  "harvestc", "prod_turnoverc", "NBP", "transp", "evap", "interc", "wd",
  "discharge", "prec", "SoilC", "SoilC_slow", "LitC", "VegC", "ProductC",
  "nuptake", "ndemand", "nlosses", "ninflux", "estabn", "negn_fluxes",
  "firen", "harvestN", "prod_turnovern", "SoilN", "SoilN_slow", "LitN",
  "VegN", "ProductN", "estab_storageC", "estab_storageN"
)
colnames(globalflux_old) <- columns
colnames(globalflux_new) <- columns

# Add scenario & combine data
globalflux_old$Scenario <- "Old Land Use"
globalflux_new$Scenario <- "New Land Use from 2006"
globalflux_combined <- bind_rows(globalflux_old, globalflux_new)

# Filter and process data
filtered_data <- globalflux_combined %>%
  filter(Year >= 1901 & Year <= 2100) %>%
  filter(!is.na(GPP) & !is.na(NBP) & !is.na(RH)) %>%
  group_by(Scenario) %>%
  arrange(Year, .by_group = TRUE) %>%
  mutate(
    GPP_MA20 = rollmean(GPP, 20, fill = "extend", align = "center"),
    GPP_NBP_Ratio = GPP / NBP
  ) %>%
  ungroup() %>%
  mutate(
    Decade = cut(
      Year,
      breaks = seq(1900, 2110, by = 10),
      labels = paste(seq(1900, 2100, by = 10), "s", sep = ""),
      right = FALSE
    )
  )

# Statistical analysis
t_test_result <- t.test(GPP ~ Scenario, data = filtered_data %>% filter(Year >= 2006))
p_value <- t_test_result$p.value

mean_gpp <- filtered_data %>% group_by(Scenario) %>% summarize(Mean_GPP = mean(GPP, na.rm = TRUE))

# Define color palette and theme
line_colors <- c("Old Land Use" = "#377EB8", "New Land Use from 2006" = "#E41A1C")
fill_colors <- c("Old Land Use" = "#A6CEE3", "New Land Use from 2006" = "#FB9A99")

custom_theme <- theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    text = element_text(color = "black", family = "sans"),  # Standard Schriftart ohne Kapitälchen
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold", size = 14),
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    legend.position = "bottom"
  )

# GPP trend plot with black curve before 2006
gpp_trend_plot <- ggplot() +
  # Black line before 2006
  geom_line(data = filtered_data %>% filter(Year < 2006), 
            aes(x = Year, y = GPP), color = "black", size = 0.8) +
  
  # Colored lines after 2006
  geom_line(data = filtered_data %>% filter(Year >= 2006), 
            aes(x = Year, y = GPP, color = Scenario), size = 0.8) +
  
  # Moving Average vor 2006 schwarz (für beide Szenarien gemeinsam)
  geom_line(data = filtered_data %>% filter(Year < 2006), 
            aes(x = Year, y = GPP_MA20), 
            linetype = "dashed", size = 0.6, color = "black", show.legend = FALSE) +
  
  # Moving Average ab 2006 für Old Land Use (grau)
  geom_line(data = filtered_data %>% filter(Year >= 2006 & Scenario == "Old Land Use"), 
            aes(x = Year, y = GPP_MA20), 
            linetype = "dashed", size = 0.6, color = "#808080", show.legend = FALSE) +
  
  # Moving Average ab 2006 für New Land Use (grau)
  geom_line(data = filtered_data %>% filter(Year >= 2006 & Scenario == "New Land Use from 2006"), 
            aes(x = Year, y = GPP_MA20), 
            linetype = "dashed", size = 0.6, color = "#B0B0B0", show.legend = FALSE) +
  
  
  # Display p-value
  #annotate("text", x = 2020, y = max(filtered_data$GPP, na.rm = TRUE) * 0.95, 
           #label = paste("p-value (2006+):", round(p_value, 5)), size = 5, hjust = 0)+
  
  # Define colors for scenarios
  scale_color_manual(values = line_colors) +
  
  # Axis formatting
  scale_x_continuous(breaks = seq(1900, 2100, by = 10)) +
  
  # Labels and theme
  labs(title = "GPP Time Series (SSP126, 1901-2100)", x = "", y = "GPP (10^15 g C / Year)") +
  custom_theme

# GPP boxplot without jitter points in legend
gpp_boxplot <- ggplot(filtered_data, aes(x = Decade, y = GPP, fill = Scenario)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, size = 0.5) +
  geom_jitter(aes(color = Scenario), width = 0.2, alpha = 0.3, size = 1, show.legend = FALSE) +  # Hide points in legend
  scale_fill_manual(values = line_colors, name = "Land Use Scenario") +
  labs(title = "GPP Distribution per Decade (SSP126, 1900-2100)", x = "", y = "GPP (10^15 g C / Year)") +
  custom_theme 
  #stat_compare_means(aes(label = ..p.format..), label.x = 1.5, label.y = max(filtered_data$GPP))

# Display plots
print(gpp_trend_plot)
print(gpp_boxplot)

# Save plots
ggsave("GPP_Trend_Plot_HighQuality.png", plot = gpp_trend_plot, dpi = 600, width = 14, height = 8, bg = "white")
ggsave("GPP_Boxplot_HighQuality.png", plot = gpp_boxplot, dpi = 600, width = 14, height = 8, bg = "white")
