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

/**
 * Class containing static utils for executing R scripts
 */
public class RUtil {
    public static final String PACKAGES_PATH = "$ExecUtil.MY_PATH/Rpackages/" // Local R library path
    public static final String NA = "NA", NULL = "NULL" // NaN and null in R
    public static boolean REMOVE_R_SCRIPTS = false

    private static String separatorsToSystem(String path) {
        File.separatorChar as String == '\\' ? path.replaceAll("/", "\\\\")
                .replaceAll("\\\\+", "\\\\\\\\") : path // Don't even try to ask why
    }

    /**
     * Converts a given object to numeric variable
     * @param smth object to convert
     * @return a numeric string or NA if object couldn't be converted
     */
    public static String asNumeric(smth) {
        if (smth == null)
            return NULL
        def smthStr = smth.toString()
        if (!smthStr.isDouble())
            return NA
        def value = smthStr.toDouble()
        (value.isNaN() || value.isInfinite()) ? NA : value.toString()
    }

    /**
     * Converts a given object to logical variable
     * @param smth object to convert
     * @return "T" for true or "F" for false
     */
    public static String logical(smth) {
        smth ? "T" : "F"
    }

    /**
     * Execute a given R script
     * @param scriptName R script name
     * @param params script parameters
     */
    public static void execute(String scriptName, String... params) {
        // Create a temp file to store the script
        def scriptRes = CommonUtil.resourceStreamReader("rscripts/$scriptName")
        //scriptName = UUID.randomUUID().toString() + "_" + scriptName

        def scriptFile = new File(scriptName)

        scriptFile.withPrintWriter { pw ->
            // Set up library path correctly
            // Don't do anything if packages are not installed
            // as this would misguide R not to use /usr/ library
            pw.println("#args <- c(${params.collect { '"' + it + '"' }.join(", ")})")

            if (new File(PACKAGES_PATH).exists())
                pw.println(separatorsToSystem(".libPaths(\"$PACKAGES_PATH\")"))

            // Write the rest of script to temp file
            scriptRes.readLines().each {
                pw.println(it)
            }
        }

        if (REMOVE_R_SCRIPTS) {
            scriptFile.deleteOnExit()
        }

        // Run script
        try {
            runScript(scriptName, params)
        } catch (IOException e) {
            if (e.toString().contains("error=2")) {
                println "[ERROR] Rscript not found, looks like R is not installed"
            }
        }
    }

    static void runScript(String scriptName, String... params) {
        def cmd = ["Rscript", scriptName, params]

        println "[RUtil] Executing ${cmd.flatten().join(" ")}"

        def proc = cmd.flatten().execute()

        proc.in.eachLine {
            println(it)
        }

        proc.out.close()
        proc.waitFor()

        if (proc.exitValue()) {
            println "[ERROR] ${proc.getErrorStream()}"
        }
    }

    /**
     * Install specified R packages
     * @param dependencies names of packages to install (case-sensitive)
     */
    public static void install(String... dependencies) {
        new File(PACKAGES_PATH).mkdirs()
        execute("install.r", [PACKAGES_PATH, dependencies].flatten() as String[])
    }

    /**
     * Test if specified R packages are installed
     * @param dependencies names of packages that will be checked (case-sensitive)
     */
    public static void test(String... dependencies) {
        execute("test.r", dependencies)
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
