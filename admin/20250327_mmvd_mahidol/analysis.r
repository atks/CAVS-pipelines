library(ggplot2)
library(dplyr)
library(tidyr)

#setwd("/home/atks/analysis/admin/20250327_mmvd_mahidol")
setwd("/Users/atks/analysis/20250327_mmvd_mahidol")

# preparation of mmvd.txt
#data is translated into long format with 8 variables
#
#investigator e.g. Sirada
#group e.g. Self-study vet student
#dog_id e.g. MU-1113-01
#video_id e.g. v15. (the video ID that appends to dog_id for.)
#frame_id e.g. f1 (the different frames in a video)
#method e.g. rishniw 
#la e.g. measurement of lateral atrium
#ao e.g. measurement of the aorta based on the different methods

data <- read.table("mmvd.txt",
    sep = "\t", header = TRUE, quote = "",
    fill = TRUE, comment.char = ""
)

#compute the ratios of LA/AO
data$r <- data$la / data$ao

#figure 1
#have a quick plot of the ratios for each dog.
#each dog is expected to have values somewhat that are somewhat close together.
ggplot(data, aes(x = dog_id, y = r)) +
    geom_boxplot()

#pick up outliers
subset(data, r>7)

############################################################
#try on your own, remove the outliers and replot the boxplot
#
#try plotting
#only AO measurements
#only rishniw measurements
#only residents
#compare the variance of the measurements for each group, this is a measure of consistency
############################################################

head(data)

subset(data, investigator=="Alisa" & dog_id=="MU-1113-01")

#aggregate matched measurements into a wide table
data_wide <- as.data.frame(data %>%
  select(-r) %>%
  group_by(investigator, dog_id, video_id, frame_id) %>%
  mutate() %>%  
  pivot_wider(
    names_from = method,
    values_from = ao,
  ))

#populate with the corresponding LA/AO
data_wide$rishniw_r = data_wide$la / data_wide$rishniw
data_wide$hansson_r = data_wide$la / data_wide$hansson
head(data_wide, 20)

subset(data_wide, is.na(hansson))

#figure 2
#view the correlation between both methods
#each dog is expected to have values somewhat that are somewhat close together.
ggplot(data_wide, aes(x=rishniw_r, y=hansson_r)) + geom_point()
  
data_wide <- data %>%
  group_by(dog_id, video_id, frame_id) %>%
  mutate(obs_id = row_number()) %>%         # Add a row number within each group
  pivot_wider(
    names_from = obs_id,
    values_from = value,
    names_prefix = "value"
  )

#to examine consistency of either methods, 
data.rishniw <- subset(data, method == "rishniw")
ggplot(data.rishniw, aes(x = dog_id, y = r)) +
    geom_boxplot()

