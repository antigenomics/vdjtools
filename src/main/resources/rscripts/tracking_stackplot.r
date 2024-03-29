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

require(ggplot2); require(reshape); require(gridExtra); require(grid)

label    <- args[1] #"time since HSCT, months"
points   <- args[2] #"-48,-0.5,4,10,25,37,72"
file_in  <- args[3] #"luc_table_collapsed.txt"
file_out <- args[4] #"out.pdf"

#

# load time points, create some auxillary variables
x <- apply(as.vector(read.table(text=points, sep=",")), 1, as.numeric)
n <- length(x)

# load data
#df <- data.frame(read.delim(file_in))
df    <- read.table(file_in, comment="", sep = "\t", header = TRUE)
xcols <- (ncol(df) - n + 1):ncol(df)
fcols <- 1:(ncol(df) - n)
xlbls <- colnames(df)[xcols]

# convert abundance columns to numeric
df[, xcols] <- apply(df[, xcols], 2, as.numeric)

# set up Non-overlapping and Not-shown
df$peak[nrow(df)-1] <- -1
df$peak[nrow(df)] <- -2

# for label placement
cs <- cumsum(df[,xcols])
peak <- df$peak + 1
peak[peak < 1] <- NA
lbl.peak <- as.vector(df[!is.na(peak), "peak"])
lbl.txt <- as.vector(df[!is.na(peak), "cdr3aa"])
lbl.x <- as.vector(na.omit(x[peak]))
lbl.y <- sapply(as.vector(na.omit(cs[cbind(1:nrow(cs), peak)])), as.numeric)
lbl.max <- sapply(as.vector(na.omit(df[cbind(1:nrow(df), xcols[peak])])), as.numeric)
lbl.p <- as.numeric(as.vector(df[!is.na(peak), "sampling.p"]))+1e-20

lbl <- data.frame(txt = lbl.txt, x = lbl.x, y = lbl.y, max = lbl.max, p.adj = lbl.p)

# reshape data
df.m <- melt(df, id = fcols)

# replace sample ids (factor) by time (numeric)
ind <- match(df.m$variable, xlbls)
df.m$variable <- x[ind]

df.m$sign <- paste(df.m$cdr3nt, df.m$v, df.m$d, df.m$j, sep="_")

pal <- colorRampPalette(c("#2b8cbe", "#e0f3db", "#fdbb84"))

#draw

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

# prepare plot
g<-ggplot() +
   geom_area(data = df.m,
             aes(variable, value, group = sign, fill = factor(peak)),
             colour = "gray25", position = 'stack', size = 0) +
   scale_fill_manual(
      name   = "Peak position",
      breaks = c(-2, -1, 0:(n-1)),
      labels = c("Non-overlapping", "Not-shown", x),
      values = c("grey50", "grey70", pal(n))
      ) +
   scale_x_continuous(expand = c(0,0), limit = c(min(x), max(x)), breaks = x) +
   scale_y_continuous(expand = c(0,0)) +
   theme_bw() +
   xlab(label) +
   ylab("abundance") +
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5), panel.border = element_blank()) +
   geom_text(data = subset(lbl, peak == 1), aes(x, y, label = txt, size = max, alpha = -log10(p.adj)), color = "gray10", hjust = 0, vjust = 0.5) +
   geom_text(data = subset(lbl, peak > 1 & peak < n), aes(x, y, label = txt, size = max, alpha = -log10(p.adj)), color = "gray10", hjust = 0.5, vjust = 0.5) +
   geom_text(data = subset(lbl, peak == n), aes(x, y, label = txt, size = max, alpha = -log10(p.adj)), color = "gray10", hjust = 1, vjust = 0.5) +
   guides(size = F)

# disable label cropping
gg_table <- ggplot_gtable(ggplot_build(g))
gg_table$layout$clip[gg_table$layout$name=="panel"] <- "off"

grid.draw(gg_table)

dev.off()
