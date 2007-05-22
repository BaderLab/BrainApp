package org.baderlab.csplugins.brainplugin.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

/**
 * Copyright (c) 2005 Memorial Sloan-Kettering Cancer Center
 * *
 * * Code written by: Gary Bader
 * * Authors: Gary Bader, Chris Sander
 * *
 * * This library is free software; you can redistribute it and/or modify it
 * * under the terms of the GNU Lesser General Public License as published
 * * by the Free Software Foundation; either version 2.1 of the License, or
 * * any later version.
 * *
 * * This library is distributed in the hope that it will be useful, but
 * * WITHOUT ANY WARRANTY, WITHOUT EVEN THE IMPLIED WARRANTY OF
 * * MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  The software and
 * * documentation provided hereunder is on an "as is" basis, and
 * * Memorial Sloan-Kettering Cancer Center
 * * has no obligations to provide maintenance, support,
 * * updates, enhancements or modifications.  In no event shall the
 * * Memorial Sloan-Kettering Cancer Center
 * * be liable to any party for direct, indirect, special,
 * * incidental or consequential damages, including lost profits, arising
 * * out of the use of this software and its documentation, even if
 * * Memorial Sloan-Kettering Cancer Center
 * * has been advised of the possibility of such damage.  See
 * * the GNU Lesser General Public License for more details.
 * *
 * * You should have received a copy of the GNU Lesser General Public License
 * * along with this library; if not, write to the Free Software Foundation,
 * * Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
 * *
 * * User: GaryBader
 * * Date: Oct 28, 2005
 * * Time: 2:14:45 PM
 */
public class FileReaderUtil {
    /**
     * Utility class to read a file as a list of lines in the file.
     * Best for small files, otherwise you should use the BufferedReader directly
     *
     * @param inputFile The file to read from
     * @return A list containing each line in the file as a separate list element
     * @throws java.io.IOException If there is a problem reading the file.
     */
    public static ArrayList readFileAsLineList(File inputFile) throws IOException {
        ArrayList list = new ArrayList();
        BufferedReader br = new BufferedReader(new FileReader(inputFile));
        String fileLine = null;
        while ((fileLine = br.readLine()) != null) {
            list.add(fileLine);
        }
        return list;
    }
}
