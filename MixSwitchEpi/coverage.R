library(covr)
c <- package_coverage()
df <- as.data.frame(c)
agg <- aggregate(df$value > 0, list(df$filename), sum)
colnames(agg) <- c("filename", "covered")

total <- aggregate(df$value >= 0, list(df$filename), sum)
colnames(total) <- c("filename", "total")

res <- merge(agg, total)
res$percent <- res$covered / res$total * 100
res <- res[order(res$percent), ]
print("Coverage by file:")
print(res)
write.csv(res, "coverage_summary.csv", row.names = FALSE)
