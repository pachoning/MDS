source("tools/load_libraries.R")
df_2000 = fread(
  file.path(
    "data/us-population-by-zip-code",
    "population_by_zip_2000.csv"
  ),
  sep = ',',
  stringsAsFactors = FALSE
)


dim(df_2000)
tail(df_2000)
df_2000$year = 2000



df_2010 = read_delim(
  file.path(
    "data/us-population-by-zip-code",
    "population_by_zip_2010.csv"
  ),
  delim = ','
)

dim(df_2010)
head(df_2010)
df_2010$year = 2010


df_pop = rbind(
  df_2000,
  df_2010
)
dim(df_pop)
save(df_pop, file = "data/us-population-by-zip-code/df_pop.RData")
