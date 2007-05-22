package org.baderlab.csplugins.brainplugin;

import org.biojava.bio.BioException;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

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
 * * Date: Nov 4, 2005
 * * Time: 11:56:15 AM
 */

/**
 * Reads profiles from a single file or a project file
 */
public class PeptideToProfileReader {
    /**
     * Reads a profile file from disk - either a single profile or a project file
     * Uniques peptides by default
     *
     * @return A List of profiles
     */
    public static List readPeptidesAsProfiles(File projectFile, double fuzzFactor) {
        return (readPeptidesAsProfiles(projectFile, -1, null, fuzzFactor, null, true));
    }

    /**
     * Reads a profile file from disk - either a single profile or a project file
     * Uniques peptides by default
     *
     * @return A List of profiles
     */
    public static List readPeptidesAsProfiles(File projectFile, double fuzzFactor, File codonBiasFile, boolean uniquePeptides) {
        return readPeptidesAsProfiles(projectFile, -1, null, fuzzFactor, codonBiasFile, uniquePeptides);
    }

    //TODO: read peptide files as a list of BindingPeptideFiles

    /**
     * Reads a profile file from disk - either a single profile or a project file
     *
     * @param profileLength The length of the profiles to return
     * @param terminus      Which end of the profile to return if length is less than the profile length
     * @return A List of profiles
     */
    public static List readPeptidesAsProfiles(File projectFile, int profileLength, ProteinTerminus terminus, double fuzzFactor, File codonBiasFile, boolean uniquePeptides) {
        ProfileSearchProjectFile profileProjectFile = new ProfileSearchProjectFile();
        BindingPeptideList peptideList;
        ProteinProfile profile;
        DatabaseReference proteinReference;
        ArrayList profileList = new ArrayList();

        if (projectFile == null) {
            throw new RuntimeException("Project file must be specified.");
        }

        //test the file to see what type it is: profile or project?
        Iterator files = null;
        try {
            if (profileProjectFile.isProjectFile(projectFile.toString())) {
                //read project file
                profileProjectFile.read(projectFile.toString());
                files = profileProjectFile.getFileNames();
            } else {
                //read a single profile (list of peptides)
                ArrayList fileList = new ArrayList(1);
                fileList.add(projectFile.toString());
                files = fileList.iterator();
            }
            while (files.hasNext()) {
                String peptideListFileName = (String) files.next();
                peptideList = new BindingPeptideList();
                System.out.println("Loading peptide sequences from " + peptideListFileName);
                peptideList.read(peptideListFileName);
                proteinReference = peptideList.getProteinXref();
                String profileName = new String(peptideList.getProteinName() + "-" + Integer.toString(peptideList.getDomainNumber()));
                if (profileLength > 0) {
                    profile = new ProteinProfile(peptideList.getSequenceIteratorByLength(profileLength, terminus, false, uniquePeptides), fuzzFactor, profileName);
                } else {
                    profile = new ProteinProfile(peptideList.getSequenceIterator(false, uniquePeptides), fuzzFactor, profileName);
                }
                if (codonBiasFile != null) {
                    profile.reWeightByCodonBias(codonBiasFile);
                }                 
                profile.setProteinReference(proteinReference);
                profile.setDomainReference(peptideList.getDomainXref());
                profile.setDomainNumber(peptideList.getDomainNumber());
                profile.setDomainSequence(peptideList.getDomainSequence());
                profile.setExperimentalMethod(peptideList.getExperimentalMethod());
                profile.setDomainSequenceStart(peptideList.getDomainRangeStart());
                profile.setDomainSequenceStop(peptideList.getDomainRangeStop());
                profile.setComment(peptideList.getComment());
                profile.setProteinName(peptideList.getProteinName());
                profileList.add(profile);
            }
        } catch (IOException e) {
            e.printStackTrace();
        } catch (BioException e) {
            e.printStackTrace();
        }

        return (profileList);
    }
}
