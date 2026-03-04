devtools::load_all()
res <- MixSwitchEpi::example_comp(
  res = 4, target_r0 = 2, pathogen_string = "sars-cov-2", calibration_mode = "switch"
)
print("res columns:")
print(colnames(res))
if ("r0" %in% names(res)) {
  print(head(as.data.frame(res)[, c("epsilon", "switch_rate", "r0")]))
  print(table(res$epsilon, res$r0))
}
