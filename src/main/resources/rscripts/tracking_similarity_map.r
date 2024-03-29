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
args <- commandArgs(TRUE)

require(ggplot2); require(reshape); require(gridExtra); require(grid)

file_in  <- args[1]
file_out <- args[2]

# preprocess data

df <- read.table(file_in, header=TRUE, comment.char="", sep = "\t")

df$value <- as.numeric(as.character(df$value)) # yet again don't ask me why

get_plot <- function(metric_value) {
   ggplot(subset(df, metric == metric_value), aes(factor(X1_time), factor(X2_time))) +
       geom_tile(aes(fill = value), colour = "white") +
       scale_fill_gradient(low = "white", high = "steelblue") +
       coord_fixed(ratio = 1) +
       facet_wrap( ~ metric) +
       scale_x_discrete(expand = c(0, 0)) +
       scale_y_discrete(expand = c(0, 0)) +
       xlab("") + ylab("") +
       #guides(fill = guide_legend(title = "")) +
       theme(
                legend.title = element_blank(),
                legend.text = element_text(angle = 45, hjust = 0),
                legend.position = "bottom",
                axis.ticks = element_blank(),
                axis.text.x = element_text(angle = 90, hjust = 0)
                )
}

g1 <- get_plot("count")
g2 <- get_plot("diversity")
g3 <- get_plot("frequency")

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
             
grid.arrange(g1, g2, g3, ncol=3, nrow=1)

dev.off()
