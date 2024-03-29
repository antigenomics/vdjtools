/*
 * Copyright (c) 2014-2024, OOO «MiLaboratory»
 *
 * IN NO EVENT SHALL THE INVENTORS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
 * SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
 * ARISING OUT OF THE USE OF THIS SOFTWARE, EVEN IF THE INVENTORS HAS BEEN
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * THE SOFTWARE PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND THE LICENSOR HAS NO
 * OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
 * MODIFICATIONS. THE LICENSOR MAKES NO REPRESENTATIONS AND EXTENDS NO
 * WARRANTIES OF ANY KIND, EITHER IMPLIED OR EXPRESS, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR
 * PURPOSE, OR THAT THE USE OF THE SOFTWARE WILL NOT INFRINGE ANY PATENT,
 * TRADEMARK OR OTHER RIGHTS.
 */

package com.antigenomics.vdjtools.overlap.permutations

class DiscreteFactorClusterStats {
    private final List<String> factorList = new ArrayList<>()
    private final List<Point> pointList = new ArrayList<>()
    private final Silhouette observedSilhouette

    public DiscreteFactorClusterStats(String fileName) {
        def reader = new File(fileName).newReader()
        reader.readLine() // id\tlabel\tfactor\tx\ty

        def line
        while ((line = reader.readLine()) != null) {
            def splitLine = line.split("\t")
            def factor = splitLine[2]
            def x = splitLine[3].toDouble(), y = splitLine[4].toDouble()

            factorList.add(factor)
            pointList.add(new Point(x, y))
        }

        this.observedSilhouette = new Silhouette(factorList, pointList)
    }

    public HashMap<String, Summary> performPermutations(int nPerms) {
        def summaryByFactor = new HashMap<String, Summary>()

        factorList.each { String factor ->
            def obsDist = observedSilhouette[factor]
            summaryByFactor.put(factor, new Summary(nPerms, obsDist.within, obsDist.between))
        }

        if (summaryByFactor.keySet().size() < 2)
            return null

        for (int i = 0; i < nPerms; i++) {
            def factorListCopy = new ArrayList<String>(factorList)
            Collections.shuffle(factorListCopy)
            def permutedSilhouette = new Silhouette(factorListCopy, pointList)

            factorList.unique(false).each { String factor ->
                def permDist = permutedSilhouette[factor]
                summaryByFactor[factor].add(permDist)
            }
        }

        summaryByFactor
    }

    public static void writeSummary(HashMap<String, Summary> summaryByFactor,
                                    String cacheFileName) {
        new File(cacheFileName).withPrintWriter { pw ->
            pw.println("factor\ttype\tperm\tobs\tp")

            summaryByFactor.each {
                def factor = it.key
                def summary = it.value
                summary.nPerms.times { int i ->
                    // todo: finish with real silhouette index
                    //pw.println([factor,
                    //            summary.getWithinPerm(i), summary.getBetweenPerm(i),
                    //            summary.withinObs, summary.betweenObs,
                    //            summary.withinP, summary.betweenP].join("\t"))

                    pw.println([factor, "within",
                                summary.getWithinPerm(i),
                                summary.withinObs,
                                summary.withinP].join("\t"))

                    pw.println([factor, "between",
                                summary.getBetweenPerm(i),
                                summary.betweenObs,
                                summary.betweenP].join("\t"))
                }
            }
        }

    }

    List<String> getFactorList() {
        Collections.unmodifiableList(factorList)
    }

    private static class Silhouette {
        private final HashMap<String, Distances> distancesMap = new HashMap<>()

        Silhouette(List<String> factorList, List<Point> pointList) {
            int n = factorList.size()

            factorList.each { distancesMap.put(it, new Distances()) }

            for (int i = 0; i < n; i++) {
                for (int j = i + 1; j < n; j++) {
                    double dist = pointList[i].dist(pointList[j])

                    if (factorList[i] != factorList[j]) {
                        distancesMap[factorList[i]].addBetween(dist)
                        distancesMap[factorList[j]].addBetween(dist)
                    } else {
                        distancesMap[factorList[i]].addWithin(dist)
                    }
                }
            }
        }

        public Distances getAt(String factor) {
            distancesMap[factor]
        }
    }

    public static class Summary {
        private int nWithin, nBetween, n
        private final int nPerms
        private final double withinObs, betweenObs
        private final double[] withinPerm, betweenPerm

        Summary(int nPerms, double withinObs, double betweenObs) {
            this.nPerms = nPerms
            this.withinObs = withinObs
            this.betweenObs = betweenObs
            this.withinPerm = new double[nPerms]
            this.betweenPerm = new double[nPerms]
        }

        public void add(Distances distancesPerm) {
            if (distancesPerm.between > betweenObs)
                nBetween++
            if (distancesPerm.within < withinObs)
                nWithin++

            betweenPerm[n] = distancesPerm.between
            withinPerm[n] = distancesPerm.within

            n++
        }

        public int getnPerms() {
            n
        }

        public double getBetweenPerm(int i) {
            betweenPerm[i]
        }

        public double getWithinPerm(int i) {
            withinPerm[i]
        }

        public double getBetweenObs() {
            betweenObs
        }

        public double getWithinObs() {
            withinObs
        }

        public double getBetweenP() {
            nBetween / (double) n
        }

        public double getWithinP() {
            nWithin / (double) n
        }
    }

    private static class Distances {
        private double withinDistSum = 0, betweenDistSum = 0
        private int kWithin = 0, kBetween = 0

        public void addWithin(double dist) {
            withinDistSum += dist
            kWithin++
        }

        public void addBetween(double dist) {
            betweenDistSum += dist
            kBetween++
        }

        public double getWithin() {
            withinDistSum / kWithin
        }

        public double getBetween() {
            betweenDistSum / kBetween
        }
    }

    private static class Point {
        private final double x, y

        Point(double x, double y) {
            this.x = x
            this.y = y
        }

        public double dist(Point other) {
            double dx = x - other.x, dy = y - other.y
            Math.sqrt(dx * dx + dy * dy)
        }
    }
}
NA
NA
NA
NA
NA
NA
NA
NA
NA
NA
NA
NA
NA
NA
NA
NA
