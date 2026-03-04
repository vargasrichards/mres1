devtools::load_all()
res <- MixSwitchEpi::example_comp(
  res = 4, target_r0 = 2, pathogen_string = "sars-cov-2", calibration_mode = "switch"
)
df <- res$long_df
df_r0 <- df[df$variable == "r0" & df$panel == "metric", ]
mat <- tapply(df_r0$value, list(epsilon = df_r0$epsilon, switch = df_r0$switch_rate), mean)
print("R0 values (rows=epsilon, cols=switch_rate):")
print(mat)
