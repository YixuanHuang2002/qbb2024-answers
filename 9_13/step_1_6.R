# step 1.6
# Define the file name for the PNG output
png(filename = "~/qbb2024-answers/9_13/ex1_30x_cov.png", width = 800, height = 600)
# Load the coverage data into R
coverage_data <- read.table("~/qbb2024-answers/9_13/coverage_30x.txt", header = TRUE)

# Total number of genome positions
N <- 1e6


# Plot histogram of coverage
hist(coverage_data$X0, breaks = 70, probability = TRUE, main = "Coverage Histogram", xlab = "Coverage Depth", ylab = "Frequency", col = "lightgreen", border = "black")

# Generate a sequence of coverage depths
coverage_depths <- 0:max(coverage_data$X0)

# Poisson distribution with lambda = 3
lambda <- 30
poisson_probs <- dpois(coverage_depths, lambda)
lines(coverage_depths, poisson_probs, col = "red", lwd = 2, type = "b")

# Normal distribution with mean = 3 and std dev = sqrt(3)
mean_normal <- 30
sd_normal <- sqrt(30)
normal_probs <- dnorm(coverage_depths, mean = mean_normal, sd = sd_normal)
lines(coverage_depths, normal_probs, col = "blue", lwd = 2, type = "b")

legend("topright", legend = c("Poisson Distribution", "Normal Distribution"), 
       col = c("red", "blue"), lwd = 3)

# Count the number of positions with 0x coverage
zero_coverage <- sum(coverage_data$X0 == 0)

# Calculate the percentage of the genome with 0x coverage
percentage_zero_coverage <- (zero_coverage / N) * 100

poisson_zero_coverage_prob <- dpois(0, lambda)


# Save the plot and close the device
dev.off()
