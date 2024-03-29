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

package com.antigenomics.vdjtools.io

import com.antigenomics.vdjtools.sample.metadata.MetadataUtil
import com.antigenomics.vdjtools.misc.CommonUtil

/**
 * A file input stream factory. This factory creates a new file connection each time.
 */
public class FileInputStreamFactory implements InputStreamFactory {
    private final String fileName

    /**
     * Creates a new instance of file input stream factory associated with a given file name
     * @param fileName path to underlying file
     */
    public FileInputStreamFactory(String fileName) {
        this.fileName = fileName
    }

    /**
     * @inheritDoc
     */
    @Override
    public InputStream create() {
        CommonUtil.getFileStream(fileName)
    }

    /**
     * @inheritDoc
     */
    @Override
    public String getId() {
        MetadataUtil.fileName2id(fileName)
    }

    String getFileName() {
        return fileName
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
