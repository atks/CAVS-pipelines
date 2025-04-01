library(ggplot2)


setwd("/home/atks/analysis/admin/20250327_mmvd_mahidol")

# preparation of mmvd2.txt
# remove "'" in the header line
# replace spaces with _

data <- read.table("mmvd2.txt",
    sep = "\t", header = TRUE, quote = "",
    fill = TRUE, comment.char = ""
)

# possible dog ID errors
#     25 MU-653-01
#     11 MU-653-10
#
#     11 records to be renamed as follows
#     MU-653-10-08  => MU-653-01-08
#

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

head(data)

# plot dog

length(unique(data$Echocardiogram_code_.HN.LOOP._))
length(unique(data$dog))



# difference between Rishiniw and Hansson



unique(data$dog)

which(count.fields("mmvd2.txt", sep = "\t") != 30)


count.fields("mmvd2.txt", sep = "\t")



# variance due to frames (within video)

# variance due to video taken (between videos)
