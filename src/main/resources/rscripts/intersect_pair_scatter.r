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
# Prints a scatterplot with marginal distributions
# Note that while scatterplot is based only on vector of species present in sample intersection,
# marginal plots are built based on complete datasets
#
# Adapted from R-bloggers (http://www.r-bloggers.com/ggplot2-cheatsheet-for-visualizing-distributions/)
# With minor additions based on various StackOverflow replies
#

#options(error=traceback)

require(ggplot2); require(grid); require(gridExtra); require(reshape); require(grid)

# data input

args<-commandArgs(TRUE)

sample1_id <- args[1]
sample2_id <- args[2]
file_xy    <- args[3]
file_xx    <- args[4]
file_yy    <- args[5]
file_out   <- args[6]

# transform data
to_double = function(x) {
  log10(as.numeric(as.character(x))+1e-7)
}

xy <- read.table(file_xy, header=T)
xy$x <- to_double(xy$x)
xy$y <- to_double(xy$y)

xx <- read.table(file_xx, header=T)
xx$xx <- to_double(xx$xx)
#xx <- rbind(subset(melt(xy),variable=="x"), melt(xx))

yy <- read.table(file_yy, header=T)
yy$yy <- to_double(yy$yy)
#yy <- rbind(subset(melt(xy),variable=="y"), melt(yy))

xmin <- min(xx$xx, yy$yy)

# For regression info plotting. source: http://goo.gl/K4yh
lm_eqn = function(df){
  m = lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(coef(m)[1], digits = 2), 
                        b = format(coef(m)[2], digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq))               
}

# plotting
 
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

# placeholder (top right)

empty <- ggplot() +
     geom_point(aes(1,1), colour="white") +
     theme(                              
       plot.background  = element_blank(), 
       panel.grid.major = element_blank(), 
       panel.grid.minor = element_blank(), 
       panel.border     = element_blank(), 
       panel.background = element_blank(),
       axis.title.x     = element_blank(),
       axis.title.y     = element_blank(),
       axis.text.x      = element_blank(),
       axis.text.y      = element_blank(),
       axis.ticks       = element_blank(),
       plot.margin      = unit(c(3, -5.5, 4, 3), "mm")
     )

# scatterplot

scatter <- ggplot() + 
  theme_bw() +
  geom_point(data = xy, aes(x, y, size = (x + y) / 2),
     fill   = "red", 
     colour = "black",
     alpha  = 0.4,
     pch    = 21
     ) +
  geom_text(data = data.frame(), aes(x = xmin, y = 0, label = lm_eqn(xy)), hjust = 0, parse = TRUE) +
  stat_smooth(data = xy, aes(x, y, weight = 10^((x + y) / 2)), color= "gray25", method = "lm", fullrange = T) +
  scale_x_continuous(limit = c(xmin, 0), expand = c(0, 0.1)) +
  scale_y_continuous(limit = c(xmin, 0), expand = c(0, 0.1)) +
  scale_size_continuous(guide = "none", range = c(1, 10)) +
  xlab(sample1_id) +
  ylab(sample2_id) +  
  theme(legend.position = c(1, 1), legend.justification = c(1, 1))

# marginal density of x

plot_top <- ggplot() +  
  stat_density(data=xx, aes(x=xx, weight=10^xx/sum(10^xx), y = ..scaled..),
               fill = "grey50", colour = "gray25", size = 0.1, alpha = 0.4, adjust = 1) +
  stat_density(data=xy, aes(x=x, weight=10^x/sum(10^x), y = ..scaled..), 
               fill = "red", colour = "gray25", size = 0.1, alpha = 0.4, adjust = 1) +
  scale_x_continuous(limit = c(xmin, 0), expand = c(0, 0.25)) +
  ylab("") + theme_bw() +
  theme(legend.position = "none", axis.title.x = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(), 
        panel.grid  = element_blank())


# marginal density of y

plot_right <- ggplot() +  
  stat_density(data=yy, aes(x=yy, weight=10^yy/sum(10^yy), y = ..scaled..),
               fill = "grey50", colour = "gray25", size = 0.1, alpha = 0.4, adjust = 1) +
  stat_density(data=xy, aes(x=y, weight=10^y/sum(10^y), y = ..scaled..), 
               fill = "red", colour = "gray25", size = 0.1, alpha = 0.4, adjust = 1) +
  scale_x_continuous(limit = c(xmin, 0), expand = c(0, 0.2)) +
  coord_flip() + 
  ylab("") + theme_bw() +
  theme(legend.position = "none", axis.title.y = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid  = element_blank()) 

# arrange the plots together, with appropriate height and width for each row and column

grid.arrange(plot_top, empty, scatter, plot_right, ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))

dev.off()
