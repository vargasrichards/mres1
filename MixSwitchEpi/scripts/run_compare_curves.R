# Run the main comparison figure: faceted curves (left) + grouped bars (right)

devtools::load_all()

# Parameters
N <- 1e6
r0 <- 2
num_activity <- 5
haslemere_activity_scores <- c(0.936825, 2.182864, 3.492136, 5.271117, 9.467058)

# Run and save outputs
out <- MixSwitchEpi::compare_curves_main(
  N = N,
  r0 = r0,
  num_activity = num_activity,
  haslemere_activity_scores = haslemere_activity_scores,
  save = TRUE,
  output_dir = "output/curve_compare"
)

print(names(out))
