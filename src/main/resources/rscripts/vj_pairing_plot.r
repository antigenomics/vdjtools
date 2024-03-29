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
args<-commandArgs(TRUE)

file_in  <- args[1]
file_out <- args[2]

require(circlize); require(RColorBrewer)


# load data and preproc to fit formats

temp <- read.table(file_in, sep="\t", comment="")
n <- nrow(temp)
m <- ncol(temp)
rn = as.character(temp[2:n,1])
cn = apply(temp[1,2:m], 2 , as.character)
mat <- matrix(apply(temp[2:n, 2:m], 1:2, as.numeric), n - 1, m-1) * 100

n <- nrow(temp)
m <- ncol(temp)

# Here columns and rows correspond to V and J segments respectively
# Also replace possible duplicates (undef, '.', ...)

duplicates <- intersect(rn, cn)

rownames(mat) <- replace(rn, rn==duplicates, paste("V", duplicates, sep=""))
colnames(mat) <- replace(cn, cn==duplicates, paste("J", duplicates, sep=""))

# sort

col_sum = apply(mat, 2, sum)
row_sum = apply(mat, 1, sum)

mat <- mat[order(row_sum), order(col_sum)]

# equal number of characters for visualizaiton

rn <- rownames(mat)
cn <- colnames(mat)

maxrn <- max(nchar(rn))
maxcn <- max(nchar(cn))

for(i in seq_len(length(rn))) {
      rn[i] <- paste(rn[i], paste(rep(" ", maxrn - nchar(rn[i])), collapse = ''))
}

for(i in seq_len(length(cn))) {
      cn[i] <- paste(cn[i], paste(rep(" ", maxcn - nchar(cn[i])), collapse = ''))
}

rownames(mat) <- rn
colnames(mat) <- cn

# viz using circlize

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

circos.par(gap.degree = c(rep(1, nrow(mat)-1), 10, rep(1, ncol(mat)-1), 15), start.degree = 5)

rcols <- rep(brewer.pal(12, "Paired"), nrow(mat)/12 + 1)[1:nrow(mat)]
ccols <- rep(brewer.pal(12, "Paired"), ncol(mat)/12 + 1)[1:ncol(mat)]

names(rcols) <- sort(rownames(mat))
names(ccols) <- sort(colnames(mat))

chordDiagram(mat, annotationTrack = "grid",
             grid.col = c(rcols, ccols),
             preAllocateTracks = list(track.height = 0.2), transparency = 0.5)

circos.trackPlotRegion(track.index = 1, bg.border = NA,
       panel.fun = function(x, y) {
                   sector.name = get.cell.meta.data("sector.index")
                   xlim = get.cell.meta.data("xlim")
                   ylim = get.cell.meta.data("ylim")
                   circos.text(mean(xlim), ylim[1], cex = 0.5, sector.name, facing = "clockwise", adj = c(0, 0.5))
                   }
       )

circos.clear()

dev.off()
