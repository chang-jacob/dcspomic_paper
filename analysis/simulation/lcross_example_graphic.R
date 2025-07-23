library(spatstat)

# Set output directory
output_dir <- "output/simulation/lcross_example"
dir.create(output_dir, showWarnings = FALSE)

set.seed(123)

# Homogeneous Poisson
pp <- rpoispp(lambda = 100, win = square(1))

# Save point pattern
png(file.path(output_dir, "poisson_pattern.png"), width = 600, height = 600)
svgg(file.path(output_dir, "poisson_pattern.png"), width = 600, height = 600)

plot(pp, main = "Homogeneous Poisson Point Pattern")
dev.off()

# Save L-function
png(file.path(output_dir, "poisson_l_function.png"), width = 600, height = 600)
plot(Lest(pp), main = "L-function for Homogeneous Process")
dev.off()

# Thomas Process
pp_clustered <- rThomas(kappa = 10, scale = 0.05, mu = 10, win = square(1))

# Save clustered point pattern
png(file.path(output_dir, "thomas_pattern.png"), width = 600, height = 600)
plot(pp_clustered, main = "Thomas Process Point Pattern")
dev.off()

# Save L-function for clustered process
png(file.path(output_dir, "thomas_l_function.png"), width = 600, height = 600)
plot(Lest(pp_clustered), main = "L-function for Thomas Process")
dev.off()
