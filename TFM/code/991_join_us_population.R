df_2000 = fread(
  file.path(
    "data/us-population-by-zip-code",
    "population_by_zip_2000.csv"
  ),
  sep = ',',
  stringsAsFactors = FALSE
)


nrow(df_2000)
tail(df_2000)
df_2000$year = 2000
dim(df_2000)


df_2010 = read_delim(
  file.path(
    "data/us-population-by-zip-code",
    "population_by_zip_2010.csv"
  ),
  delim = ','
)

head(df_2010)
df_2010$year = 2010
dim(df_2010)
