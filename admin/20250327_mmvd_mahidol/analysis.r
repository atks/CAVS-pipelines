library(ggplot2)

setwd("/home/atks/analysis/admin/20250327_mmvd_mahidol")

# preparation of mmvd2.txt
# remove "'" in the header line
# replace spaces with _

data <- read.table("mmvd.txt",
    sep = "\t", header = TRUE, quote = "",
    fill = TRUE, comment.char = ""
)

data$r <- data$la / data$ao

subset(data, r > 10)

ggplot(data, aes(x = dog_id, y = r)) +
    geom_boxplot()

data.rishniw <- subset(data, method == "rishniw")



nrow(data)

nrow(data.rishniw)

ggplot(data.rishniw, aes(x = dog_id, y = r)) +
    geom_boxplot()

data.hansson <- subset(data, method == "hansson")

ggplot(data.hansson, aes(x = dog_id, y = r)) +
    geom_boxplot()


# group by dog_id, video_id and frame_id





head(data)
colnames(data)

# examine data type
unique(data$Investigator)

unique(data$Group)

unique(data$Echocardiogram_code_.HN.LOOP._)

data$Echocardiogram_code_.HN.LOOP._

strsplit(data$Echocardiogram_code_.HN.LOOP._, split = "-")[1:3, ]

# extract the dog
data$dog <- sapply(
    strsplit(data$Echocardiogram_code_.HN.LOOP._, "-"),
    function(x) paste(x[1:3], collapse = "-")
)


# plot dog

length(unique(data$Echocardiogram_code_.HN.LOOP._))
length(unique(data$dog))



# difference between Rishiniw and Hansson



unique(data$dog)

which(count.fields("mmvd2.txt", sep = "\t") != 30)


count.fields("mmvd2.txt", sep = "\t")



# variance due to frames (within video)

# variance due to video taken (between videos)
