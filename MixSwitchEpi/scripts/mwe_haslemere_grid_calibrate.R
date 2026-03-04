# mwe_haslemere_grid_calibrate.R
#
# Minimal working example: Haslemere empirical 1-D scan with and without
# per-residence-time recalibration. Save each row (R0) and calibration setup
# to separate PDF files.

pkgload::load_all(quiet = TRUE)

# ── Parameters ────────────────────────────────────────────────────────────────
R0_LIST <- c(1.5, 2, 2.5, 3) 

# ── Run the grid ───────────────────────────────────
message("Running Haslemere grid scan (calibrate_each = TRUE) …")
grid_df <- haslemere_grid(r0_list = R0_LIST, calibrate_each = TRUE)

# ── Save separate file for each row and calibration condition ───
out_dir <- "output/haslemere_1d_separated_panels"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
message(sprintf("\nSaving separated panel rows to %s/ ...", out_dir))

r0s <- sort(unique(grid_df$calib_at_r0))
calibs <- unique(grid_df$calibration)

# ── Save full grid plots (all rows together) ───
for (cal in calibs) {
  if (is.na(cal)) next
  sub_df_full <- grid_df[grid_df$calibration == cal, ]
  if (nrow(sub_df_full) > 0) {
    cal_str <- tolower(gsub(" ", "_", cal))
    fname_full <- file.path(out_dir, sprintf("haslemere_all_R0_%s.pdf", cal_str))
    plot_haslemere_grid(sub_df_full, save_path = fname_full)
  }
}
fname_both_full <- file.path(out_dir, "haslemere_all_R0_both.pdf")
plot_haslemere_grid(grid_df, save_path = fname_both_full)

for (r in r0s) {
  # 1) Save specific calibrations separately (Fixed p / Recalibrated)
  for (cal in calibs) {
    if (is.na(cal)) next
    
    sub_df <- grid_df[grid_df$calib_at_r0 == r & grid_df$calibration == cal, ]
    if (nrow(sub_df) == 0) next
    
    # format safe filename
    cal_str <- tolower(gsub(" ", "_", cal))
    fname <- file.path(out_dir, sprintf("haslemere_R0_%.1f_%s.pdf", r, cal_str))
    
    plot_haslemere_grid(sub_df, save_path = fname)
  }
  
  # 2) Also optionally save the combined row with both lines plotted together
  sub_df_both <- grid_df[grid_df$calib_at_r0 == r, ]
  if (nrow(sub_df_both) > 0) {
    fname_both <- file.path(out_dir, sprintf("haslemere_R0_%.1f_both.pdf", r))
    plot_haslemere_grid(sub_df_both, save_path = fname_both)
  }
}

message("Done! Check ", out_dir)
