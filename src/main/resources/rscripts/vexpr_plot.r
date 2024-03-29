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
require(gplots); require(RColorBrewer); require(ggplot2); require(plotrix)

args<-commandArgs(TRUE)

input_file  <- args[1] #"ms.segments.V2.txt"
v_col_count <- as.numeric(args[2])
lbl_col     <- as.numeric(args[3])
fac_col     <- as.numeric(args[4])
cont_factor <- as.logical(args[5]) # continuous factor?
file_out <- args[6]

if (lbl_col < 1) {
   lbl_col = 1 # use sample id if not specified
}

color_by_factor <- fac_col > 0 # any color legend for samples?

# read-in data and pre-process
df<-read.table(input_file, header=T, sep="\t", comment="")
vcols<-(ncol(df)-v_col_count + 1):ncol(df)
df[, vcols] <- apply(df[, vcols], 2, as.numeric)
x <- as.matrix(df[, vcols])

# prepare palletes and factor levels
if (color_by_factor) {
    factor_name <- colnames(df)[fac_col]
    if (cont_factor) {
       df[, fac_col] <- as.numeric(df[,fac_col])
       # base pallete, to be used in legend
       scol <- colorRampPalette(c("#feb24c", "#31a354", "#2b8cbe"))(10)
       # color vector for plot
       pcol <- scol[cut(df[,fac_col], 10)]
    } else {
       df[, fac_col] <- as.factor(df[, fac_col])
       fu   <- levels(df[, fac_col])
       nlevels <- length(fu)
       # base pallete, to be used in legend
       scol <- colorRampPalette(brewer.pal(nlevels, "Set2"))(nlevels)
       # color vector for plot
       pcol <- scol[df[,fac_col]]
    }
}

# plotting function for simplicity
my.plot <- function(...) {
   heatmap.2(t(x),  na.rm=F, labCol = df[,lbl_col],
      na.col="grey50",
      col=colorRampPalette(c("#2f98ce", "#e0f3db", "#f47104")),
      density.info="none", trace="none",scale="column",
      ...)
}

custom.dev <- function(fname) {
   if (grepl("\\.pdf$",fname)){
      pdf(fname)
   } else if (grepl("\\.png$",fname)) {
      png(fname, width     = 3.25,
                 height    = 3.25,
                 units     = "in",
                 res       = 1200,
                 pointsize = 4)
   } else {
      stop('Unknown plotting format')
   }
}

if (color_by_factor) {
   custom.dev(file_out)

   # layout plot
   fig <- c(0.05, 0.75, 0, 1.0)
   mar <- c(0, 0, 0, 0)

   layout(matrix(1:2, ncol=2), width = c(2, 0.1), height = c(1, 1))

   par(fig = fig, mar = mar, xpd = NA) # this ensures labels are not cut

   # plot with coloring columns
   my.plot(ColSideColors=pcol)

   # draw a separate legend
   fig <- c(0.85, 1.00, 0, 1.0)

   par(fig = fig, mar = c(0, 0, 0, 0), xpd = NA, new=TRUE)

   if (cont_factor) {
      # set colorbar pos
      f <- 2
      ux <- grconvertX(c(0.48 - 0.05 * f, 0.48 + 0.05 * f), from = "npc", to = "user")
      uy <- grconvertY(c(0.87 - 0.07, 0.87 + 0.07), from = "npc", to = "user")

      # gradient legend
      fmin <- min(df[, fac_col])
      fmax <- max(df[, fac_col])
      color.legend(ux[1], uy[1], ux[2], uy[2],
                   legend = c(fmin, (fmin+fmax) / 2, fmax),
                   rect.col = scol, gradient = "y", align="rb")

      # position title
      uy <- grconvertY(0.96, from = "npc", to = "user")
      text(0.5, uy, factor_name, adj = c(0.5, 0.0))
   } else {
      # vanilla legend for discrete factor
      legend("top", "(x,y)", fill = scol, pch = '', box.lwd = 0, title = factor_name, legend = fu)
   }

   dev.off()
} else {
   custom.dev(file_out)

   # layout plot
   fig <- c(0.05, 0.95, 0, 1.0)
   mar <- c(0, 0, 0, 0)

   layout(matrix(1:2, ncol=2), width = c(2, 0.1), height = c(1, 1))

   par(fig = fig, mar = mar, xpd = NA) # this ensures labels are not cut

   my.plot()

   dev.off()
}
