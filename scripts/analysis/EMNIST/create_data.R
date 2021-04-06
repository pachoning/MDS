library(tidyverse)
library(gtools)

# Import all data ----
train <- read_csv(file = "data/EMNIST/emnist-byclass-train.csv", col_names = FALSE)
test <- read_csv(file = "data/EMNIST/emnist-byclass-test.csv", col_names = FALSE)

mapping <- read_delim(file = "data/EMNIST/emnist-byclass-mapping.txt", col_names = FALSE, delim = " ")
colnames(mapping) <- c("id", "ascii")

# Join data ----
train_label <- train %>% 
  left_join(mapping, by = c("X1" = "id")) 

test_label <- test %>% 
  left_join(mapping, by = c("X1" = "id")) 

all_data <- rbind(train_label, test_label) %>% 
  mutate(value = chr(ascii))

all_data_pixels <- all_data[, 2:785]
colnames(all_data_pixels) <- paste0("X", 1:(28*28))
target <- all_data$value

# Save data ----
save(all_data_pixels, target, file = "data/EMNIST/all_data.RData")
