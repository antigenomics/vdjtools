#
# Copyright (c) 2014-2024, OOO «MiLaboratory»
#
# IN NO EVENT SHALL THE INVENTORS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
# SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
# ARISING OUT OF THE USE OF THIS SOFTWARE, EVEN IF THE INVENTORS HAS BEEN
# ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# THE SOFTWARE PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND THE LICENSOR HAS NO
# OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
# MODIFICATIONS. THE LICENSOR MAKES NO REPRESENTATIONS AND EXTENDS NO
# WARRANTIES OF ANY KIND, EITHER IMPLIED OR EXPRESS, INCLUDING, BUT NOT LIMITED
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR
# PURPOSE, OR THAT THE USE OF THE SOFTWARE WILL NOT INFRINGE ANY PATENT,
# TRADEMARK OR OTHER RIGHTS.
#
# data input

args <- commandArgs(TRUE)

require(gplots)

label    <- args[1] #"time since HSCT, months"
points   <- args[2] #"-48,-0.5,4,10,25,37,72"
file_in  <- args[3] #"tc_luc_table_collapsed.txt"
file_out <- args[4] #"out.pdf"

# load time points, create some auxillary variables
t <- apply(as.vector(read.table(text=points, sep=",")), 1, as.numeric)
n <- length(t)

# load data and remove collapse summary columns
df <- read.table(file_in, sep="\t", header = T, comment ="")
df <- df[1:(nrow(df)-2), ]

# convert abundance columns to numeric
tcols <- (ncol(df) - n + 1):ncol(df)
df[, tcols] <- apply(df[, tcols], 2, as.numeric)

# format (log-transform and scale)
y <- df[,tcols]
min <- min(y[y > 0])
max <- max(y)
x <- log10(as.matrix(y) + min / 10)
breaks = c(log10(min) - 1, seq(log10(min), log10(max), length.out = 20))
pal <- colorRampPalette(c("#2b8cbe", "#e0f3db", "#fdbb84"))

# plot

if (grepl("\\.pdf$",file_out)){
   pdf(file_out)
} else if (grepl("\\.png$",file_out)) {
   png(file_out, width     = 3.25,
                 height    = 3.25,
                 units     = "in",
                 res       = 1200,
                 pointsize = 4)
} else {
   stop('Unknown plotting format')
}

heatmap.2(x, labRow = df$cdr3aa, labCol = t, Colv = FALSE, dendrogram = "row",
          scale = "none", cexRow = 0.2 + 0.5 / log10(nrow(x)), keysize = 1.0, symkey = FALSE,
          col = c("grey50", pal(length(breaks) - 2)), breaks = breaks,
          density.info="none", trace="none",
          srtCol=90, adjCol = c(1,0), xlab = label)

dev.off()
