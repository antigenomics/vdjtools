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
#
# Prints a stacked area plot for two samples
#

args<-commandArgs(TRUE)

require(ggplot2); require(reshape); require(RColorBrewer); require(FField)

# data input

sample1_id <- args[1]
sample2_id <- args[2]
file_in  <- args[3]
file_out <- args[4]

# transform data

table <- read.table(file_in, header=T, comment="", sep="\t", quote="", stringsAsFactors=F)
n <- nrow(table)
m <- ncol(table)
table <- table[n:1, ] # reorder, bring non-overlapping / not-shown to the top

sample_ids <- c(sample1_id, sample2_id)

# here we do change sample labels and further refer to them as 0 and 1
# otherwise R will re-order samples alphabetically by name ...

df <- data.frame(cdr3nt = rep(table$cdr3nt, 2),
                 cdr3aa = rep(table$cdr3aa, 2),
                 v = rep(table$v, 2),
                 d = rep(table$d, 2),
                 j = rep(table$j, 2),
                 sample = rep(c(0, 1), each = n),
                 expr   = as.numeric(c(table[,m-1], table[,m])),
                 peak   = rep(as.integer(table$peak), 2))
                 
df$sign <- paste(df$cdr3nt, df$v, df$d, df$j, sep="_")

# trick to add clonotype labels with geom_text, as suggested at StackExchange

df.0 <- subset(df, sample == 0)
df.1 <- subset(df, sample == 1)

#dd <- dd[with(dd, order(cdr3nt, levels(dd$cdr3aa))), ]
df.0$cum <- cumsum(df.0$expr)
df.0$cum <- (c(0, df.0$cum[0:(n-1)]) + df.0$cum) / 2
df.1$cum <- cumsum(df.1$expr)
df.1$cum <- (c(0, df.1$cum[0:(n-1)]) + df.1$cum) / 2

# for label placement
jitter_coord = function(y) {
    if (length(y) >= 2) {
       y.fact <- 100 / max(y)
       coords <- FFieldPtRep(coords = cbind(y * y.fact, y * y.fact), iter.max = 1000)
       y <- coords$y / y.fact
       y <- (y - min(y)) / (max(y) - min(y))
    } else {
       y
    }
}

# placing
df.00 <- subset(df.0, peak == 0)
df.00$cum2 <- jitter_coord(df.00$cum)
df.11 <- subset(df.1, peak == 1)
df.11$cum2 <- jitter_coord(df.11$cum)

# custom palette (color blind)

gs.pal <- colorRampPalette(c("#2b8cbe", "#e0f3db", "#fdbb84"))

# plotting
emax <- max(df$expr)

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

ggplot() +
    geom_area(data = df, aes(x = sample, y = expr, fill = cdr3nt,
    group = sign), colour = "grey25", size = 0.01, position="stack") +
    ylab("cumulative abundance") +
    xlab("") +
    scale_x_continuous(expand = c(0, 0), labels=sample_ids, breaks=c(0,1), limits=c(-1,2)) +
    scale_y_continuous(expand = c(0, 0), limits=c(-0.2,1.2), breaks=c(0,0.25,0.5,0.75,1)) +
    scale_fill_manual(values = c("grey50", "grey25", gs.pal(n - 2))) +
    theme_bw() +
    theme(legend.position = "none") +
    geom_text(data = df.00, 
              aes(x = -0.25, y = cum2, label = cdr3aa), hjust = 1, size=3) +
    geom_segment(data = df.00, x=0,xend=-0.23, aes(y=cum, yend=cum2), size=0.3) +
    geom_text(data = df.11, 
              aes(x = 1.25, y = cum2, label = cdr3aa), hjust = 0, size=3) +  
    geom_segment(data = df.11, x=1,xend=1.23, aes(y=cum, yend=cum2), size=0.3) +  
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

dev.off()
