package org.baderlab.csplugins.brain;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;

/**
 * Copyright (c) 2004 Memorial Sloan-Kettering Cancer Center
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
 * * Date: Feb 5, 2005
 * * Time: 9:49:01 PM
 * * Description Handles I/O for a profile database search project file
 */

/**
 * Handles I/O for a profile database search project file
 * <p/>
 * File format:
 * 1st column: Filename of peptide list
 * 2nd column: Name of protein that peptides putatively bind to
 * 3rd column: Specify the database ID of the protein that peptides putatively bind to (e.g. swissprot:P1234).
 * <p/>
 * Comments lines start with a # and will be ignored
 */
public class ProfileSearchProjectFile {
    //In-memory representation of the project file
    private ArrayList fileNameList;

    /**
     * Main constructor for a project file
     */
    public ProfileSearchProjectFile() {
        fileNameList = new ArrayList();
    }

    /**
     * Get an iterator over all peptide files in the project file (represented as Strings)
     */
    public Iterator getFileNames() {
        return fileNameList.iterator();
    }

    /**
     * Reads a peptide list project file
     *
     * @param fileName The name of the project file
     * @throws IOException if there is any error in reading the file
     */
    public void read(String fileName) throws IOException {
        if (!isProjectFile(fileName)) {
            throw new IOException("Tried to read an invalid project file.");
        }
        String projectFileLine = null;
        String peptideListFileName = null;
        int lineCount = 0;
        BufferedReader br = new BufferedReader(new FileReader(fileName));
        File file = new File(fileName);
        String projectFileDir = file.getAbsoluteFile().getParent();
        while ((projectFileLine = br.readLine()) != null) {
            lineCount++;
            //check for a comment
            if ((projectFileLine.equals("")) || (projectFileLine.charAt(0) == '#')) {
                continue;
            }
            peptideListFileName = projectFileLine;
            File pepFile = new File(peptideListFileName);
            if (!pepFile.isAbsolute()) {
                //if peptide file is not absolute, add the project file's directory to make it absolute
                //This assumes that the peptide files will be in the same directory as the project file or relative to that directory.
                peptideListFileName = new String(projectFileDir + File.separator + peptideListFileName);
            } //else the filename is absolute, use it as is
            fileNameList.add(peptideListFileName);
        }
        br.close();
    }

    /**
     * Tests if a file is a project file format (checks to make sure the project file header is present)
     *
     * @param fileName The name of the project file.
     * @return true if fileName is a project file, false otherwise
     * @throws IOException if anything goes wrong with reading fileName
     */
    public boolean isProjectFile(String fileName) throws IOException {
        if (fileName == null) {
            return false;
        }
        BufferedReader br = new BufferedReader(new FileReader(fileName));
        String fileLine = br.readLine();
        fileLine = fileLine.trim();
        if ((fileLine != null) && (fileLine.equals("#ProjectFile"))) {
            return true;
        }
        return false;
    }
}
