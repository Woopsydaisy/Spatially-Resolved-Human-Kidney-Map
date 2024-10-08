# -*- coding: utf-8 -*-
"""Circular Plot R intersection fibroblast.ipynb

Automatically generated by Colab.

Original file is located at
    https://colab.research.google.com/drive/1UTxs6JrOMi0Ue5N7shQLnZCbabo8RS_i
"""

#edges <- read.csv("/content/chordplotforitaliptintersectionfibroblastreceiver.csv")
#nodes <- read.csv('nodes_bokeh_short.csv')
#edges <- read.csv("/content/chordplotforitaliptintersectionfibroblastreceiver_abbr.csv")
edges <- read.csv("/content/chordplotforitaliptintersectionfibroblastreceiver_abbr_3.csv")

edges

rownames(edges) <- edges$Index

# Remove the 'Index' column from the data frame
edges <- edges[,-1]
mat <- data.matrix(edges)
# Now 'edges' should have its row names set to the former 'Index' column
#print(edges)

rownames(edges) <- edges$X

# Remove the first column from the data frame
edges <- edges[,-1]

# Convert the rest of the data frame to a numeric matrix
mat <- data.matrix(edges)

edges

mat

#mat <- t(mat)

mat



install.packages("circlize")

# Load the circlize library
library(circlize)
## "#fdcdac", "#cbd5e8", "#f4cae4", "#e6f5c9", "#fff2ae", "#f1e2cc", "#cccccc"
# Define colors for each specified group of columns
# Adjust the color codes as per your column groupings
state_col = c("#",  # Color for C1 (Blue)
    "#DCDC8D")
# Assuming 'mat' is the matrix representation of 'edges'
mat <- as.matrix(edges)

, # Color (Orange)
    "#4B9F34",  # Color for C6 (Green)

    "#D7263D",  # Vibrant Red
    "#F9E07F",  # Bright, Sunny Yellow
    "#6B4C9A",  # Deep Purple
    "#FF6B6B",  # Soft, Bright Pink
    "#C0D6DF",  # Light Blue
    "#F9844A",  # Muted Orange
    "#6D98BA",  # Steel Blue
    "#926AA6",  # Soft Purple
    "#40394A"    )  # Color for (C17-C19)##orange:  previously: #FDB352

#state_col_transparent <- sapply(state_col, adjustcolor, alpha.f=0.9)
#col_assignments <- c(rep(state_col[1], 11),     #B           # C1
 #                    rep(state_col[2], 13))  #cDC

# Define colors for the columns
state_col = c("#C8C7C7", "#DCDC8D")
col_assignments <- c(rep(state_col[1], 11), rep(state_col[2], 13))

actual_col_count <- ncol(mat)
col_assignments <- col_assignments[1:actual_col_count]

# Create the color matrix for the chord diagram
colmat <- matrix(col_assignments, nrow = nrow(mat), ncol = ncol(mat), byrow = TRUE)

print(length(col_assignments))
print(ncol(mat))

print(head(col_assignments))

print(unique(col_assignments))

colmat

edges

mat

circos.par(cell.padding = c(0, 0, 0, 0), points.overflow.warning = FALSE)

# Generate the chord diagram with the customized colors and additional settings
cdm_res <- chordDiagram(mat, col = colmat,
    directional = TRUE, annotationTrack = "grid",
    big.gap = 20, small.gap = 2,
    preAllocateTracks = list(track.height = 0.1),
    link.target.prop = FALSE, )

circos.clear()

circos.clear()

##    E1 E2 E3 E4 E5 E6
## S1  4 14 13 17  5  2
## S2  7  1  6  8 12 15
## S3  9 10  3 16 11 18

circos.par(gap.after = c(rep(5, nrow(mat)-1), 15, rep(5, ncol(mat)-1), 15))

nrow(mat)

colnames(mat)



circos.clear()
gap.after <- c(rep(2, 10), 1, rep(2, 9))  # 10 rows + 10 columns = 20 sectors, adjusting gaps

# General circlize parameters
circos.par(cell.padding = c(0, 0, 0, 0), points.overflow.warning = FALSE,
           track.margin = c(0.01, 0.01), gap.after = gap.after)

# Generate the chord diagram with customized gaps
cdm_res <- chordDiagram(
  mat,
  col = colmat,
  directional = TRUE,
  annotationTrack = "grid",
  preAllocateTracks = list(track.height = 0.1),
  link.target.prop = FALSE
)

# Reset circlize parameters to default
circos.clear()

n <- nrow(mat)
gap_after <- c(rep(5, 10), 15, rep(5, n-11))

gap_after

circos.clear()
circos.par(gap.after = c(rep(11, nrow(mat)-1), 10, rep(11, ncol(mat)-1), 15), points.overflow.warning = FALSE)

# Generate the chord diagram with the customized colors and additional settings
cdm_res <- chordDiagram(mat, col = colmat,
    directional = TRUE)

circos.clear()

mat <- mat * 170.6

mat

mat

colmat

# Assuming 'mat' and 'colmat' are correctly prepared with original names
cdm_res <- chordDiagram(mat, col = colmat, directional = FALSE,
                        annotationTrack = "grid", big.gap = 20, small.gap = 2,
                        preAllocateTracks = list(track.height = 0.1),
                        link.target.prop = FALSE, annotationTrackHeight = c(0.05, 0.1),
                        link.lwd = 0.01)

# Adding text labels with original column names
circos.track(track.index = 2, panel.fun = function(x, y) {
  sector_name <- CELL_META$sector.index
  # Optionally, map sector_name back to a more meaningful label if necessary
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], sector_name,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA)

# Adjust the gaps
n <- nrow(mat)
gap_after <- c(rep(5, 10), 15, rep(5, n-11))

# Set the circos parameters with the adjusted gaps
circos.par(gap.after = gap_after)

# Generate the chord diagram with the new gaps
cdm_res <- chordDiagram(mat, col = colmat, directional = FALSE,
                        annotationTrack = "grid", big.gap = 20, small.gap = 2,
                        preAllocateTracks = list(track.height = 0.1),
                        link.target.prop = FALSE, annotationTrackHeight = c(0.05, 0.1),
                        link.lwd = 0.01, reduce = 0)

# Adding text labels with original column names
circos.track(track.index = 2, panel.fun = function(x, y) {
  sector_name <- CELL_META$sector.index
  # Optionally, map sector_name back to a more meaningful label if necessary
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], sector_name,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA)

# Reset circos parameters to default if necessary
circos.clear()

gap_after

library(circlize)

# Set the gap parameters
gap_size <- 1
large_gap_size <- 15

# Create the gap.after vector with sector names
sector_names <- colnames(mat)
gap_after <- rep(gap_size, length(sector_names))
gap_after[11] <- large_gap_size  # Setting the large gap after the 10th column
names(gap_after) <- sector_names

# Set the gap parameters using circos.par
circos.par(gap.after = gap_after)

# Draw the chord diagram with the specified gaps
cdm_res <- chordDiagram(mat, col = colmat, directional = FALSE,
                        annotationTrack = "grid", reduce = 0, # Set reduce to 0
                        preAllocateTracks = list(track.height = 0.1),
                        link.target.prop = FALSE, annotationTrackHeight = c(0.05, 0.1),
                        link.lwd = 0.01)

# Adding text labels with original column names
circos.track(track.index = 2, panel.fun = function(x, y) {
  sector_name <- CELL_META$sector.index
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], sector_name,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA)

# Reset circos parameters to default
circos.clear()

mat_expanded

mat

# Load necessary library
library(circlize)

# Specify the file name and path to save the diagram
output_file <- "chord_diagram.png"

# Open a PNG device with transparency and high resolution
png(filename = output_file, width = 6, height = 6, units = 'in', res = 900, bg = "transparent")

# Create the chord diagram
chordDiagram(edges, transparency = 0.5, grid.col = entity_colors, annotationTrack = "grid", preAllocateTracks = list(track.height = 0.1))

# Customize the track plot region
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  # You may adjust labels.cex to change the font size of the axis labels
  circos.axis(h = "top", labels.cex = 0.5, major.tick.length = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)

# Turn off the device to save the file
dev.off()

# Specify the file name and path to save the diagram
output_file <- "chord_diagram.png"

# Open a PNG device with transparency and high resolution
png(filename = output_file, width = 13, height = 13, units = 'in', res = 900, bg = "transparent")

# Example with a single uniform gap
cdm_res <- chordDiagram(mat, col = colmat, directional = TRUE,
                        annotationTrack = NULL,
                        preAllocateTracks = list(track.height = 0.1),
                        link.target.prop = FALSE,
                        big.gap = 15) # Adjust as necessary

# Turn off the device to save the file
dev.off()

# Specify the file name and path to save the diagram
output_file <- "chord_diagram_outside.png"

# Open a PNG device with transparency and high resolution
png(filename = output_file, width = 13, height = 13, units = 'in', res = 900, bg = "transparent")


options(repr.plot.width=13, repr.plot.height=13)

# Assuming 'mat' and 'colmat' are correctly prepared with original names
cdm_res <- chordDiagram(mat, col = colmat, directional = FALSE,
                        annotationTrack = "grid", big.gap = 10, small.gap = 1,
                        preAllocateTracks = list(track.height = 0.1),
                        link.target.prop = FALSE, annotationTrackHeight = c(0.05, 0.1))

# Adding text labels with original column names
circos.track(track.index = 2, panel.fun = function(x, y) {
  sector_name <- CELL_META$sector.index
  # Optionally, map sector_name back to a more meaningful label if necessary
  circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + mm_y(3), sector_name,
              facing = "clockwise", labels.facing="outside", adj = c(0, 0.5))
}, bg.border = NA)
dev.off()

gaps = rep(2, (ncol(mat)+ nrow(mat)))
names(gaps) = c(colnames(mat), rownames(mat))
gaps['TNFSF10.0']=5
gaps['TGFB1_2_3']=30
gaps['TNFRSF11B']=30
gaps=gaps[names(gaps)!='ALCAM']
gaps

output_file <- "chord_diagram_combined.png"
# Open a PNG device with transparency and high resolution
png(filename = output_file, width = 10, height = 10, units = 'in', res = 450, bg = "white")
circos.par(gap.after=gaps)
# Example with a single uniform gap and no additional annotation track
cdm_res <- chordDiagram(mat, col = colmat, directional = FALSE,
                        annotationTrack = NULL,
                        preAllocateTracks = list(track.height = 0.1),
                        link.target.prop = FALSE,
                        reduce = 0, transparency=0.25)
                         # Adjust as necessary
# Add labels to the diagram without an extra track
circos.track(track.index = 1, panel.fun = function(x, y) {
  sector_name <- CELL_META$sector.index
  # Optionally, map sector_name back to a more meaningful label if necessary
  circos.text(CELL_META$xcenter, CELL_META$cell.ylim[1] + mm_y(1), sector_name,
              facing = "clockwise", labels.facing="outside", adj = c(0, 0.5),
              cex = 1.5) # Adjust text size as necessary
}, bg.border = NA)
circos.clear()
# Turn off the device to save the file
dev.off()

output_file <- "chord_diagram_combined.png"

# Open a PNG device with smaller dimensions
png(filename = output_file, width = 10, height = 10, units = 'in', res = 450, bg = "white")

# Set circos parameters to adjust the size of the circle
circos.par(gap.after=gaps,
           canvas.xlim=c(-1.5, 1.5),
           canvas.ylim=c(-1.5, 1.5),
           track.margin=c(0.02, 0.02))

# Example with a single uniform gap and no additional annotation track
cdm_res <- chordDiagram(mat, col = colmat, directional = FALSE,
                        annotationTrack = NULL,
                        preAllocateTracks = list(track.height = 0.1),
                        link.target.prop = FALSE,
                        reduce = 0, transparency=0.25)

# Add labels to the diagram without an extra track
circos.track(track.index = 1, panel.fun = function(x, y) {
  sector_name <- CELL_META$sector.index
  # Optionally, map sector_name back to a more meaningful label if necessary
  circos.text(CELL_META$xcenter, CELL_META$cell.ylim[1] + mm_y(1), sector_name,
              facing = "clockwise", labels.facing="outside", adj = c(0, 0.5),
              cex = 1.3) # Adjust text size as necessary
}, bg.border = NA)

# Clear the circos plot
circos.clear()

# Turn off the device to save the file
dev.off()

circos.clear()
# Specify the file name and path to save the diagram
output_file <- "chord_diagram_combined.png"

# Open a PNG device with transparency and high resolution
png(filename = output_file, width = 15, height = 15, units = 'in', res = 450, bg = "white")
circos.par(gap.after = c(rep(2, nrow(mat)-1), 25, rep(2, 11),25,rep(2, (ncol(mat)-13)), 25))
# Example with a single uniform gap and no additional annotation track
cdm_res <- chordDiagram(mat, col = colmat, directional = FALSE,
                        annotationTrack = NULL,
                        preAllocateTracks = list(track.height = 0.1),
                        link.target.prop = FALSE,
                        reduce=0) # Adjust as necessary

# Add labels to the diagram without an extra track
circos.track(track.index = 1, panel.fun = function(x, y) {
  sector_name <- CELL_META$sector.index
  # Optionally, map sector_name back to a more meaningful label if necessary
  circos.text(CELL_META$xcenter, CELL_META$cell.ylim[1] + mm_y(1), sector_name,
              facing = "clockwise", labels.facing="outside", adj = c(0, 0.5),
              cex = 2.2) # Adjust text size as necessary
}, bg.border = NA)

# Turn off the device to save the file
dev.off()

circos.clear()



# Specify the file name and path to save the diagram
output_file <- "chord_diagram_combined.png"

# Open a PNG device with transparency and high resolution
png(filename = output_file, width = 15, height = 15, units = 'in', res = 450, bg = "white")

# Example with a single uniform gap and no additional annotation track
cdm_res <- chordDiagram(mat, col = colmat, directional = FALSE,
                        annotationTrack = NULL,
                        preAllocateTracks = list(track.height = 0.1),
                        link.target.prop = FALSE,
                        small.gap = 3,
                        big.gap = 20) # Adjust as necessary

# Add labels to the diagram without an extra track
circos.track(track.index = 1, panel.fun = function(x, y) {
  sector_name <- CELL_META$sector.index
  # Optionally, map sector_name back to a more meaningful label if necessary
  circos.text(CELL_META$xcenter, CELL_META$cell.ylim[1] + mm_y(1), sector_name,
              facing = "clockwise", labels.facing="outside", adj = c(0, 0.5),
              cex = 2.2) # Adjust text size as necessary
}, bg.border = NA)

# Turn off the device to save the file
dev.off()

# Specify the file name and path to save the diagram
output_file <- "chord_diagram_combined_black.png"

# Open a PNG device with transparency and high resolution
png(filename = output_file, width = 15, height = 15, units = 'in', res = 600, bg = "black")

# Example with a single uniform gap and no additional annotation track
cdm_res <- chordDiagram(mat, col = colmat, directional = FALSE,
                        annotationTrack = NULL,
                        preAllocateTracks = list(track.height = 0.1),
                        link.target.prop = FALSE,
                        small.gap = 1.5,
                        big.gap = 15) # Adjust as necessary

# Add labels to the diagram without an extra track
circos.track(track.index = 1, panel.fun = function(x, y) {
  sector_name <- CELL_META$sector.index
  # Optionally, map sector_name back to a more meaningful label if necessary
  circos.text(CELL_META$xcenter, CELL_META$cell.ylim[1] + mm_y(1), sector_name,
              facing = "clockwise", labels.facing="outside", adj = c(0, 0.5),
              cex = 1.6, col = "white", fontfamily = "Arial") # Adjust text size as necessary
}, bg.border = NA)

# Turn off the device to save the file
dev.off()