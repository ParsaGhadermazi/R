# Load Packages
library(readr)
library(dplyr)
library(magrittr)
library(tximport)
Sample_Table=read_csv("SraRunTable.txt") %>%
  select('Sample Name',source_name,treatment,Cell_Line,Cell_type,time_point) %>%
  slice(seq(1,48,by=4))
tximport(files = ,type = "salmon",tx2gene = ,
           ignoreTxVersion = TRUE)