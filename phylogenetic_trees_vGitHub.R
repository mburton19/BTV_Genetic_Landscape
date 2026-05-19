# Load libraries
library(phytools)   # For plotting tanglegrams
library(ape)        # For reading trees
library(readxl)     # To rename trees

# Step 1: Read in new trees
# Make sure these are the newest tree files
Segment_1 <- read.tree("Segment_1_Alignment_01022026.fasta.treefile")
Segment_2 <- read.tree("Segment_2_Alignment_01022026.fasta.treefile")
Segment_3 <- read.tree("Segment_3_Alignment_01022026.fasta.treefile")
Segment_4 <- read.tree("Segment_4_Alignment_01022026.fasta.treefile")
Segment_5 <- read.tree("Segment_5_Alignment_01102026.fasta.treefile")
Segment_6 <- read.tree("Segment_6_Alignment_01022026.fasta.treefile")
Segment_7 <- read.tree("Segment_7_Alignment_01022026.fasta.treefile")
Segment_8 <- read.tree("Segment_8_Alignment_01022026.fasta.treefile")
Segment_9 <- read.tree("Segment_9_Alignment_01022026.fasta.treefile")
Segment_10 <- read.tree("Segment_10_Alignment_01082026.fasta.treefile")

# Step 2: Rename trees
lookup <- read_excel("new_tree_ID.xlsx")
colnames(lookup) <- c("tree_label", "new_ID")
new_names <- setNames(lookup$new_ID, lookup$tree_label)

rename_tips <- function(tree, map) {
  tree$tip.label <- ifelse(
    tree$tip.label %in% names(map),
    map[tree$tip.label],
    tree$tip.label
  )
  return(tree)
}

segment_names <- paste0("Segment_", 1:10)
segments <- mget(segment_names)
segments <- lapply(segments, rename_tips, map = new_names)
list2env(segments, envir = .GlobalEnv)

# Step 3: Midpoint root trees
Segment_1_rooted <- midpoint_root(Segment_1)
Segment_2_rooted <- midpoint_root(Segment_2)
Segment_3_rooted <- midpoint_root(Segment_3)
Segment_4_rooted <- midpoint_root(Segment_4)
Segment_5_rooted <- midpoint_root(Segment_5)
Segment_6_rooted <- midpoint_root(Segment_6)
Segment_7_rooted <- midpoint_root(Segment_7)
Segment_8_rooted <- midpoint_root(Segment_8)
Segment_9_rooted <- midpoint_root(Segment_9)
Segment_10_rooted <- midpoint_root(Segment_10)

# Step 4: Check to make sure tips match
# Return value of TRUE, good to move forward
setequal(Segment_2_rooted$tip.label, Segment_1_rooted$tip.label)
setequal(Segment_2_rooted$tip.label, Segment_3_rooted$tip.label)
setequal(Segment_2_rooted$tip.label, Segment_4_rooted$tip.label)
setequal(Segment_2_rooted$tip.label, Segment_5_rooted$tip.label)
setequal(Segment_2_rooted$tip.label, Segment_6_rooted$tip.label)
setequal(Segment_2_rooted$tip.label, Segment_7_rooted$tip.label)
setequal(Segment_2_rooted$tip.label, Segment_8_rooted$tip.label)
setequal(Segment_2_rooted$tip.label, Segment_9_rooted$tip.label)
setequal(Segment_2_rooted$tip.label, Segment_10_rooted$tip.label)

# Step 5: Check for polytomies
# Return value of FALSE, good to move forward
any(Segment_1_rooted$edge[,2] %in% which(table(Segment_1_rooted$edge[,1]) > 2))
any(Segment_2_rooted$edge[,2] %in% which(table(Segment_2_rooted$edge[,1]) > 2))
any(Segment_3_rooted$edge[,2] %in% which(table(Segment_3_rooted$edge[,1]) > 2))
any(Segment_4_rooted$edge[,2] %in% which(table(Segment_4_rooted$edge[,1]) > 2))
any(Segment_5_rooted$edge[,2] %in% which(table(Segment_5_rooted$edge[,1]) > 2))
any(Segment_6_rooted$edge[,2] %in% which(table(Segment_6_rooted$edge[,1]) > 2))
any(Segment_7_rooted$edge[,2] %in% which(table(Segment_7_rooted$edge[,1]) > 2))
any(Segment_8_rooted$edge[,2] %in% which(table(Segment_8_rooted$edge[,1]) > 2))
any(Segment_9_rooted$edge[,2] %in% which(table(Segment_9_rooted$edge[,1]) > 2))
any(Segment_10_rooted$edge[,2] %in% which(table(Segment_10_rooted$edge[,1]) > 2))

# Step 6: Collapse 0 branch lengths
Segment_1_rooted <- di2multi(Segment_1_rooted, tol = 0.001)
Segment_2_rooted <- di2multi(Segment_2_rooted, tol = 0.001)
Segment_3_rooted <- di2multi(Segment_3_rooted, tol = 0.001)
Segment_4_rooted <- di2multi(Segment_4_rooted, tol = 0.001)
Segment_5_rooted <- di2multi(Segment_5_rooted, tol = 0.001)
Segment_6_rooted <- di2multi(Segment_6_rooted, tol = 0.001)
Segment_7_rooted <- di2multi(Segment_7_rooted, tol = 0.001)
Segment_8_rooted <- di2multi(Segment_8_rooted, tol = 0.001)
Segment_9_rooted <- di2multi(Segment_9_rooted, tol = 0.001)
Segment_10_rooted <- di2multi(Segment_10_rooted, tol = 0.001)

# Step 7: Make cophylogeny (this makes the object but doesn't plot it yet)
tangle_2_1 <- cophylo(Segment_2_rooted, Segment_1_rooted)
tangle_2_3 <- cophylo(Segment_2_rooted, Segment_3_rooted)
tangle_2_4 <- cophylo(Segment_2_rooted, Segment_4_rooted)
tangle_2_5 <- cophylo(Segment_2_rooted, Segment_5_rooted)
tangle_2_6 <- cophylo(Segment_2_rooted, Segment_6_rooted)
tangle_2_7 <- cophylo(Segment_2_rooted, Segment_7_rooted)
tangle_2_8 <- cophylo(Segment_2_rooted, Segment_8_rooted)
tangle_2_9 <- cophylo(Segment_2_rooted, Segment_9_rooted)
tangle_2_10 <- cophylo(Segment_2_rooted, Segment_10_rooted)

# Step 8: Plot tanglegrams

plot(tangle_2_1, 
     link.type = "curved",
     link.lwd = 2,
     link.lty = "solid",
     link.col = make.transparent("gray", 0.25),
     fsize = 0.8,
     pts = FALSE,
     scale.bar = c(0.05, 0.005))

plot(tangle_2_3, 
     link.type = "curved",
     link.lwd = 2,
     link.lty = "solid",
     link.col = make.transparent("gray", 0.25),
     fsize = 0.8,
     pts = FALSE,
     scale.bar = c(0.05, 0.005))

plot(tangle_2_4, 
     link.type = "curved",
     link.lwd = 2,
     link.lty = "solid",
     link.col = make.transparent("gray", 0.25),
     fsize = 0.8,
     pts = FALSE,
     scale.bar = c(0.05, 0.005))

plot(tangle_2_5, 
     link.type = "curved",
     link.lwd = 2,
     link.lty = "solid",
     link.col = make.transparent("gray", 0.25),
     fsize = 0.8,
     pts = FALSE,
     scale.bar = c(0.05, 0.005))

plot(tangle_2_6, 
     link.type = "curved",
     link.lwd = 2,
     link.lty = "solid",
     link.col = make.transparent("gray", 0.25),
     fsize = 0.8,
     pts = FALSE,
     scale.bar = c(0.05, 0.05))

plot(tangle_2_7, 
     link.type = "curved",
     link.lwd = 2,
     link.lty = "solid",
     link.col = make.transparent("gray", 0.25),
     fsize = 0.8,
     pts = FALSE,
     scale.bar = c(0.05, 0.05))

plot(tangle_2_8, 
     link.type = "curved",
     link.lwd = 2,
     link.lty = "solid",
     link.col = make.transparent("gray", 0.25),
     fsize = 0.8,
     pts = FALSE,
     scale.bar = c(0.05, 0.005))

plot(tangle_2_9, 
     link.type = "curved",
     link.lwd = 2,
     link.lty = "solid",
     link.col = make.transparent("gray", 0.25),
     fsize = 0.8,
     pts = FALSE,
     scale.bar = c(0.05, 0.005))

plot(tangle_2_10, 
     link.type = "curved",
     link.lwd = 2,
     link.lty = "solid",
     link.col = make.transparent("gray", 0.25),
     fsize = 0.8,
     pts = FALSE,
     scale.bar = c(0.5, 0.5))


# This code was modified/adapted from: https://github.com/LKeene/DiverseCollections/blob/main/analyses/scripts/galbut_tanglegram.R
