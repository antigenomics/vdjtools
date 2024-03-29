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

package com.antigenomics.vdjtools.misc

import java.util.zip.ZipInputStream

println "[RInstall] Opening resources stream"
def src = RInstall.class.protectionDomain.codeSource,
    jar = src.location,
    zip = new ZipInputStream(jar.openStream())
def entry
def dependencies = new HashSet<String>()

println "[RInstall] Scanning for dependencies"

while ((entry = zip.nextEntry)) {
    if (entry.name.toUpperCase().endsWith(".R")) {
        println "[RInstall] Scanning $entry.name"
        CommonUtil.resourceStreamReader(entry.name).readLines().each { String line ->
            if (line =~ /require\(.+\)/)
                line.split("require\\(").each { String token ->
                    if (token.contains(")")) {
                        def dependency = token.split("\\)")[0]
                        println "$dependency"
                        dependencies.add(dependency)
                    }
                }
        }
    }
}

println "[RInstall] Full list of dependencies to be installed:\n${dependencies.join(" ")}"

RUtil.install(dependencies as String[])

println "[RInstall] Testing"

RUtil.test(dependencies as String[])

println "[RInstall] Finished"
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
