load("/Volumes/LACIE SHARE/TFM/ws_2018_12_10_18_09_12_1000000_part4.Rproj")


df_input_all = rbind(
  df_input_all,
  df
  
)

df_summary_all = rbind(
  df_summary_all,
  df_summary
)




to_delete = setdiff(
  ls(),
  c("df_input_all", "df_summary_all")
)

rm(list = to_delete)



char_time = gsub(
  pattern = "-|:| ",
  replacement = '_',
  x = Sys.time()
)

save.image(
  file = paste0("ws_gr_",char_time, ".Rproj")
)
