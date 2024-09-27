setwd("/Users/cmdb/qbb2024-answers/9_13")

# step 1.5

# create a png file
png("ex1_30x_cov.png", width = 800, height = 600)
# fetch the 30x data
coverage_data <- scan("/Users/cmdb/qbb2024-answers/9_13/coverage_30x.txt")

hist(coverage_data,
     breaks = seq(-0.5, max(coverage_data) + 0.5, by = 1),
     xlab = "Coverage",
     ylab = "Frequency",
     main = "Coverage Depth",
     col = "white",
     border = "black")

# set the parameters
coverage_values <- 0:max(coverage_data)
genome_size <- length(coverage_data)
# Poisson
lambda <- 30
# Normal
normal_mean <- 30
normal_stdv <- sqrt(30)

# the probes & counts for poisson and normal distribution
poisson_probs <- dpois(coverage_values, lambda)
poisson_counts <- poisson_probs * genome_size
normal_probs <- dnorm(coverage_values, normal_mean, normal_stdv)
normal_counts <- normal_probs * genome_size

# Poisson distribution
lines(coverage_values, poisson_counts, col = "red", lwd = 2)
# Normal distribution
lines(coverage_values, normal_counts, col = "blue", lwd = 2)

# legend
legend("topright", 
       legend = c("Poisson Distribution", "Normal Distribution"), 
       col = c("red", "blue"),
       lwd = 2)

dev.off()

