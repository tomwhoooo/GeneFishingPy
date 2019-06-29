library(stringr)

load('~/data/normalized.GTex.with.ncRNA.median.larger.than.0.1.RData')
names_data_ori <- names(normalized.GTex.data)
names_data <- str_replace_all(names(normalized.GTex.data), " ", "_") # avoid problems in linux

dir.create(file.path("~/data", "csv_data"), showWarnings = FALSE)

for(i in 1:length(names_data)){
    write.csv(as.data.frame(normalized.GTex.data[[names_data_ori[i]]]), paste0("~/data/csv_data/", names_data[i], ".csv"))
} 
#use paste0 to wipe out the blank space