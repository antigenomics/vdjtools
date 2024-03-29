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

require(ggplot2); require(RColorBrewer)

# data input

sample1_id <- args[1]
sample2_id <- args[2]
file_xy    <- args[3] # should contain cdr3nt, cdr3aa cols; two last columns should correspond to clonotype freq in sample pair
file_out   <- args[4]

# transform data

getcol_n = function(lst, col){
	apply(lst[col], 2, as.numeric)
}

getcol_c = function(lst, col){
	apply(lst[col], 2, as.character)
}

table <- read.delim(file_xy, sep="\t", quote="")
n <- nrow(table)
m <- ncol(table)
table <- table[n:1, ] # reorder, bring non-overlapping / not-shown to the top

x1 <- getcol_n(table, m - 1)
x2 <- getcol_n(table, m)

# here we do change sample labels and further refer to them as s1 and s2
# otherwise R will re-order samples alphabetically by name ...

df <- data.frame(cdr3nt = rep(getcol_c(table, "cdr3nt"), 2),
                 cdr3aa = rep(getcol_c(table, "cdr3aa"), 2),
                 v = rep(getcol_c(table, "v"), 2),
                 d = rep(getcol_c(table, "d"), 2),
                 j = rep(getcol_c(table, "j"), 2),
                 sample = rep(c("s1", "s2"), each = n),
                 expr   = c(x1, x2),
                 peak   = rep(getcol_n(table, "peak"), 2))
                 
df$sign <- paste(df$cdr3nt, df$v, df$d, df$j, sep="_")

# trick to add clonotype labels with geom_text, as suggested at StackExchange

df.0 <- subset(df, sample == "s1")
df.1 <- subset(df, sample == "s2")

#dd <- dd[with(dd, order(cdr3nt, levels(dd$cdr3aa))), ]
df.0$cum <- cumsum(df.0$expr)
df.0$cum <- (c(0, df.0$cum[0:(n-1)]) + df.0$cum) / 2
df.1$cum <- cumsum(df.1$expr)
df.1$cum <- (c(0, df.1$cum[0:(n-1)]) + df.1$cum) / 2

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
    group = sign), colour = "grey25", size = 0, position="stack") +
    ylab("cumulative abundance") +
    xlab("") +
    scale_x_discrete(expand = c(0, 0), labels=c(sample1_id, sample2_id)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_manual(values = c("grey50", "grey25", gs.pal(n - 2))) +
    theme_bw() +
    theme(legend.position = "none") + #, axis.text.x  = element_text(angle=90, vjust=0.5)) +
    geom_text(data = subset(df.0, peak == 0), aes(x = "s1", y = cum, label = cdr3aa, size = expr, alpha = rank(expr)), hjust = 0, vjust = 1) +
    geom_text(data = subset(df.1, peak == 1), aes(x = "s2", y = cum, label = cdr3aa, size = expr, alpha = rank(expr)), hjust = 1, vjust = 1) +
    scale_size_continuous(guide = "none", range = c(3, 10))

dev.off()
