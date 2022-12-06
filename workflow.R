# Workflow for batch processing Fluorocam output
# Including summarizing, standardizing, and Tcrit analysis
# Intended as a workflow script with 2 inputs
# 1: Fluorocam output
# 2: Identifying variables for each well position

# Load libraries
library(readr)
library(dplyr)
library(tidyr)
library(segmented)
library(ggplot2)

# Create directories if needed

if(!dir.exists("data_raw/")) {
  dir.create("data_raw/")
}
if(!dir.exists("data_processed/")) {
  dir.create("data_processed/")
}
if(!dir.exists("data_labels/")) {
  dir.create("data_labels/")
}


#### Add and detect new raw files ####

# NOTE: Manually add Fluorcam .txt file to 'data_raw/'
# NOTE: Manually add well names .csv file to 'data_labels/
# Both must have identical file names


new <- setdiff(list.files("data_raw/"), list.files("data_processed/"))

#### Load in data from data_raw/ ####
fn <- new[1] # Assuming only one new file

in_dat <- read_table(paste0("data_raw/", fn),
                     skip = 2) %>%
  rename(time = 1,
         temp = 2)

in_dat2 <- in_dat %>%
  mutate(temp_bin = rep(1:(nrow(in_dat)/5), each = 5)) %>% # Group every 5 measurements
  relocate(temp_bin, .after = temp)

#### Summarize mean fluorescence by temp and well ####
# Note: actual temp has been replaced by the average of the 5 temps

sum_dat <- in_dat2 %>%
  pivot_longer(-1:-3, 
               names_to = "well",
               values_to = "fluor") %>%
  group_by(well, temp_bin) %>%
  summarize(temp = mean(temp),
            fluor_mean = mean(fluor)) %>%
  ungroup() %>%
  dplyr::select(-temp_bin)

#### Detect outliers in fluorescence ####

outlier_dat <- sum_dat %>% 
  group_by(well) %>%
  mutate(fluor_thresh = median(fluor_mean) + 3 * IQR(fluor_mean), # threshold set at 3 IQRs above the median
         outlier = ifelse(fluor_mean > fluor_thresh, TRUE, FALSE)) %>%
  ungroup()

#### Scale fluorescence between 0 and 1 by well ####

scale_dat <- outlier_dat %>%
  filter(outlier == FALSE) %>% # exclude outliers
  group_by(well) %>%
  mutate(fluor_max = max(fluor_mean),
         Tmax = temp[which.max(fluor_mean)],
         fluor_min = min(fluor_mean[which(temp <= Tmax)]), # takes minimum of fluor_mean where temp was < Tmax
         Tmin = temp[which(fluor_mean == fluor_min)],
         check = ifelse(Tmax > Tmin, TRUE, FALSE),
         fluor_scale = (fluor_mean - fluor_min) / (fluor_max - fluor_min)) %>%
  ungroup()

#### Evaluate T50 and Tmax for scaled fluorescence by well ####

param_dat <- scale_dat %>% 
  dplyr::select(well, temp, fluor_scale) %>%
  group_by(well) %>%
  mutate(Tmax = temp[which(fluor_scale == 1)], # temp at scaled fluorescence = 1
         T50 = temp[which.min(abs(fluor_scale[1:which(fluor_scale == 1)] - 0.5))]) # temp when scaled fluorescence is closest to 0.5 and occurs prior to Tmax

#### Evaluate Tcrit for scaled fluorescence by well ####
# still needs conditional if segmented model does not converge

# Setting up data according to Arnold et al. 2021
xlow <- 30 # degrees C
xhigh <- 60 # degrees C
x_prebuffer <- 12 # degrees C; not sure why 12 was picked
x_postbuffer <- 1 # degrees C

# Unique wells to loop over
wells <- unique(param_dat$well)

# Dataframe to fill
Tcrits <- data.frame(well = wells,
                     Tcrit = NA,
                     Tcrit_se = NA)

# Create diretory with same name as input file in data_processed/
if(!dir.exists(paste0("data_processed/", fn))) {
  dir.create(paste0("data_processed/", fn))
}

for(w in wells) {
  # Select data
  sub_well <- param_dat %>%
    filter(well == w) # for plotting only
  
  sub <- param_dat %>%
    filter(well == w, # one well at a time
           temp > T50 - x_prebuffer, # include temps within the prebuffer of T50
           temp < T50 + x_postbuffer) # and within the postbuffer of T50
  
  # Fit breakpoint model
  m1 <- lm(fluor_scale ~ temp, data = sub)
  m2 <- segmented(m1, 
                  seg.Z = ~ temp,
                  npsi = 1)
  
  # Store breakpoint as Tcrit
  Tcrits$Tcrit[which(Tcrits$well == w)] <- round(m2$psi[2], 2)
  Tcrits$Tcrit_se[which(Tcrits$well == w)] <- round(m2$psi[3], 2)
  
  # Create fitted dataframe
  fit_dat <- data.frame(temp = sub$temp,
                        fitted = fitted(m2))
  
  # Create plot
  fig <- ggplot() +
    geom_line(data = sub_well,
               aes(x = temp,  y = fluor_scale, 
                   color = "observed")) +
    geom_line(data = fit_dat,
               aes(x = temp, y = fitted, 
                   color = "fitted")) +
    # geom_rect(aes(xmin = m2$psi[2] - m2$psi[3],
    #               xmax = m2$psi[2] + m2$psi[3],
    #               ymin = -Inf,
    #               ymax = Inf),
    #           alpha = 0.2) +
    geom_vline(aes(xintercept = m2$psi[2], 
                   color = "Tcrit"),
               lty = 2) +
    geom_vline(aes(xintercept = unique(sub_well$T50), 
                   color = "T50"),
               lty = 2) +
    geom_vline(aes(xintercept = unique(sub_well$Tmax), 
                   color = "Tmax"),
               lty = 2) +
    scale_x_continuous("Temperature (°C)") +
    scale_y_continuous("Scaled fluorescence",
                       breaks = seq(0, 1, 0.25)) +
    scale_color_manual(values = c("blue", "black",
                                  "orange", "pink", "red")) +
    theme_classic(base_size = 10)
  
  ggsave(filename = paste0("data_processed/", fn, "/", w, ".png"),
         plot = fig,
         width = 4,
         height = 3)
  
 print(paste("Well", w, "complete")) 
}

#### Load well labels ####

labs <- read_csv(paste0("data_labels/", strsplit(fn, ".TXT"), ".csv"))

#### Create final product ####

out_dat <- scale_dat %>% 
  dplyr::select(well, temp, fluor_scale) %>%
  group_by(well) %>%
  summarize(Tmax = temp[which(fluor_scale == 1)], 
            T50 = temp[which.min(abs(fluor_scale[1:which(fluor_scale == 1)] - 0.5))]) %>%
  left_join(Tcrits, by = "well") %>%
  left_join(labs, by = "well") %>%
  relocate(Tmax, T50, Tcrit, Tcrit_se, .after = last_col())


write_csv(out_dat,
          file = paste0("data_processed/", strsplit(fn, ".TXT"), ".csv"))
