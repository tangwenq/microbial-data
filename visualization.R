# Load the data for m=50, m=100, m=200, and m=400
pro.errors200 <- read.table("//fileservices.ad.jyu.fi/homes/wentang/My Documents/results/gllvm results/sim1 m=200 k=500a.txt")
pro.errors400 <- read.table("//fileservices.ad.jyu.fi/homes/wentang/My Documents/results/gllvm results/sim1 m=400 k=500a.txt")
pro.errors100 <- read.table("//fileservices.ad.jyu.fi/homes/wentang/My Documents/results/gllvm results/sim1 m=100 k=500a.txt")
pro.errors50 <- read.table("//fileservices.ad.jyu.fi/homes/wentang/My Documents/results/gllvm results/sim1 m=50 k=500a.txt")

# Create data frames for each simulation result
simres50 <- data.frame(err = as.numeric(unlist(pro.errors50)))
simres100 <- data.frame(err = as.numeric(unlist(pro.errors100)))
simres200 <- data.frame(err = as.numeric(unlist(pro.errors200)))
simres400 <- data.frame(err = as.numeric(unlist(pro.errors400)))

# Add dimension and method columns for each dataset
simres50$dimension <- "m=50"
simres100$dimension <- "m=100"
simres200$dimension <- "m=200"
simres400$dimension <- "m=400"

# Combine all data into a single data frame
simres <- rbind(simres50, simres100, simres200, simres400)

# Create method labels
methods <- c(rep("C-ZINB", B), rep("C-NB", B), rep("LVM-ZINB", B),
             rep("LVM-NB", B), rep("clr + PCA", B), rep("nMDS", B))

simres$method <- rep(methods, times = 4)  # Repeat methods for each dimension

# Rename columns
names(simres) <- c("error", "dimension", "method")

# Reorder the levels of 'dimension' factor to ensure m=50 comes first
simres$dimension <- factor(simres$dimension, levels = c("m=50", "m=100", "m=200", "m=400"))

# Reorder 'method' factor to ensure "clr + PCA" is last
simres$method <- factor(simres$method, levels = c("C-ZINB", "C-NB", "LVM-ZINB", "LVM-NB", "nMDS", "clr + PCA"))

# Plot with log scale on Y-axis, and facets for different dimensions (m=50, 100, 200, 400)
ggplot(simres, aes(x = method, y = error, fill = method)) +
  geom_boxplot(alpha = 0.3) +
  ylab("Mean Procrustes error(log)") + 
  xlab("") +  # Use log scale for Y-axis
  theme(legend.position = "none") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate X-axis labels for readability
  facet_wrap(~ dimension, scales = "fixed") +  # Create 2x2 grid with fixed Y-axis scale
  labs(fill = "Method") +
  coord_trans(y = "log") + 
  scale_y_continuous(limits = c(0.7,40)) + 
  ggtitle("gllvm") # Apply log scale to the y-axis