library(CVXR)
library(ggplot2)
n = 50
d = 100
N = 500

xi_0 = c(rep(0.3, 5), rep(0, d-5)) * sqrt(d)

norm_vec <- function(x) sqrt(sum(x^2))
proxy_Y <- function(n, d, Y, Z){
  Yhat <- Variable(n)
  objective <- Maximize(sum(Yhat * Y))
  constraint_1 <- sum(Yhat^2) <= sum(Y^2)
  vector = rep(0,d)
  for (i in 1:d){
    vector[i] = norm_vec(Z[,i])
  }
  constraint_2 <- norm_inf(t(Z) %*% Yhat) <=   1.5 * max(vector) * sqrt(log (d))
  prob <- Problem(objective, constraints = list(constraint_1, constraint_2))
  solution <- solve(prob)
  solution$getValue(Yhat)
}

ratio_1 = ratio_2 = rep(0, N)
starttime = Sys.time()
count_1 = count_2 = count_3 = 0
for (i in 1:N){
  if (i%%10 == 0){
    print(i)
    print(Sys.time() - starttime)
    starttime = Sys.time()
  }
  Z = matrix(rnorm(n * d), n, d) / sqrt(d)
  Y = c(Z %*% xi_0) + rnorm(n)
  Yhat = proxy_Y(n, d, Y, Z)
  r_1 = norm_vec(Y- c(Z %*% xi_0))
  r_2 = sum(proxy_Y(n, d, Y, Z)*(Y- c(Z %*% xi_0)))/norm_vec(proxy_Y(n, d, Y, Z))
  r_3 = sum(Y*(Y- c(Z %*% xi_0)))/norm_vec(Y)
  ratio_1[i] = r_1/r_2
  ratio_2[i] = r_1/r_3
}

data <- data.frame(
  ratio = c(ratio_1, ratio_2), 
  ratio = rep(c("proxy_Y", "original_Y"), each = N)
)
data$method = factor(data$method, levels = c("proxy_Y", "original_Y"))


summary_stats <- data %>%
  group_by(group) %>%
  summarize(
    mean = mean(value),
    sd = sd(value),
    .groups = 'drop'
  )

# Plot with density curves and annotation
ggplot(data, aes(x = value, fill = group)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Vector 1" = "blue", "Vector 2" = "red")) +
  labs(
    title = "Density Plot of Two Vectors",
    x = "Value",
    y = "Density",
    fill = "Vector"
  ) +
  # Add annotations with summary statistics
  annotate(
    "text", x = 6, y = 0.3, hjust = 1, 
    label = paste("Vector 1: Mean =", round(summary_stats$mean[1], 2), 
                  ", SD =", round(summary_stats$sd[1], 2), "\n",
                  "Vector 2: Mean =", round(summary_stats$mean[2], 2), 
                  ", SD =", round(summary_stats$sd[2], 2))
  ) +
  theme_minimal()

vector = rep(0,d)
for (i in 1:d){
  vector[i] = norm_vec(Z[,i])
}
max(vector)
