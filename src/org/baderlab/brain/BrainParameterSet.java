package org.baderlab.brain;

import java.io.File;

/**
 * * Copyright (c) 2004 Memorial Sloan-Kettering Cancer Center
 * *
 * * Code written by: Gary Bader
 * * Authors: Gary Bader, Ethan Cerami, Chris Sander
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
 ** User: Gary Bader
 ** Date: Jan 26, 2004
 ** Time: 2:44:30 PM
 ** Description: Stores a BRAIN parameter set
 **/

/**
 * Stores a BRAIN parameter set
 */
public class BrainParameterSet {
    //parameters
    //database selection
    private File databaseFileName = null;  //database to search
    private String databaseFormat = null;  //database format e.g. FASTA, etc.
    //general search parameters
    private ProteinDatabaseSearchParams searchParams = null; //e.g. length, multiple hits, terminus
    //profile search (normalized parameters is in searchParams)
    private File profileFileName = null;  //file containing profile or project file pointing to a set of profiles
    private double scoreThreshold;
    private double scorePercentageThreshold;
    private int numberTopHits;
    private double fuzzFactor;
    private boolean profileSearch;
    //for co-expression processing
    private File expressionDataSetFile = null;
    private File samplesOfInterestFile = null;
    private File nodeXrefToAffyIDMappingFile = null;
    // codon bias specification
    private File codonBiasFileName = null; // file containing codon bias file specification
    private boolean uniquePeptides;
    private boolean uniqueQueryProteinNodes;

    /**
     * Create a new profile search parameter set
     *
     * @param databaseFileName         The filename of the database to search
     * @param databaseFormat           The format of the database records (e.g. FASTA)
     * @param searchParams             The search parameters
     * @param scoreThreshold           The score threshold for the search
     * @param scorePercentageThreshold The top score percentage threshold
     * @param numberTopHits            The number of top hits to return
     * @param fuzzFactor               The fuzz factor for the profile search
     */
    public BrainParameterSet(File databaseFileName,
                             String databaseFormat,
                             ProteinDatabaseSearchParams searchParams,
                             File profileFileName,
                             File codonBiasFileName,
                             double scoreThreshold,
                             double scorePercentageThreshold,
                             int numberTopHits,
                             double fuzzFactor,
                             boolean uniquePeptides,
                             boolean uniqueQueryProteinNodes) {

        this.databaseFileName = databaseFileName;
        this.databaseFormat = databaseFormat;
        this.searchParams = searchParams;
        this.scoreThreshold = scoreThreshold;
        this.profileFileName = profileFileName;
        this.codonBiasFileName = codonBiasFileName;
        this.scorePercentageThreshold = scorePercentageThreshold;
        this.numberTopHits = numberTopHits;
        this.fuzzFactor = fuzzFactor;
        profileSearch = true;
        this.uniquePeptides = uniquePeptides;
        this.uniqueQueryProteinNodes = uniqueQueryProteinNodes;
    }

    /**
     * Create a new search parameter set (no parameters set)
     */
    public BrainParameterSet() {
        this.databaseFileName = null;
        this.databaseFormat = null;
        this.searchParams = new ProteinDatabaseSearchParams(true);
        this.scoreThreshold = 12.0;
        this.scorePercentageThreshold = 100.0;
        this.numberTopHits = 20;
        this.fuzzFactor = 1.0;
        this.uniquePeptides = true;
        this.uniqueQueryProteinNodes = false;

    }

    /**
     * Copies a parameter set object
     *
     * @return A copy of the parameter set
     */
    public BrainParameterSet copy() {
        BrainParameterSet newParam = new BrainParameterSet(this.getDatabaseFileName(), this.getDatabaseFormat(),
                this.getSearchParams().copy(), this.getProfileFile(), this.getCodonBiasFile(), this.getScoreThreshold(),
                this.getScorePercentageThreshold(), this.getNumberTopHits(), this.getFuzzFactor(), this.getUniquePeptides(),
                this.getUniqueQueryProteinNodes());
        return newParam;
    }

    //parameter getting and setting
    public File getDatabaseFileName() {
        return databaseFileName;
    }

    public String getDatabaseFormat() {
        return databaseFormat;
    }

    public ProteinDatabaseSearchParams getSearchParams() {
        return searchParams;
    }

    public double getScoreThreshold() {
        return scoreThreshold;
    }

    public File getProfileFile() {
        return profileFileName;
    }

    public File getCodonBiasFile() {
        return codonBiasFileName;
    }

    public boolean getUniquePeptides() {
        return uniquePeptides;
    }

    public boolean getUniqueQueryProteinNodes() {
        return uniqueQueryProteinNodes;
    }

    public double getScorePercentageThreshold() {
        return scorePercentageThreshold;
    }

    public int getNumberTopHits() {
        return numberTopHits;
    }

    public double getFuzzFactor() {
        return fuzzFactor;
    }

    public void setDatabaseFileName(File databaseFileName) {
        this.databaseFileName = databaseFileName;
    }

    public void setDatabaseFormat(String databaseFormat) {
        this.databaseFormat = databaseFormat;
    }

    public void setSearchParams(ProteinDatabaseSearchParams searchParams) {
        this.searchParams = searchParams;
    }

    public void setScoreThreshold(double scoreThreshold) {
        this.scoreThreshold = scoreThreshold;
    }

    public void setProfileFileName(File profileFileName) {
        this.profileFileName = profileFileName;
    }

    public void setCodonBiasFileName(File codonBiasFileName) {
        this.codonBiasFileName = codonBiasFileName;
    }

    public void setUniquePeptides(boolean value) {
        this.uniquePeptides = value;
    }

    public void setUniqueQueryProteinNodes(boolean value) {
        this.uniqueQueryProteinNodes = value;
    }

    public void setScorePercentageThreshold(double scorePercentageThreshold) {
        this.scorePercentageThreshold = scorePercentageThreshold;
    }

    public void setNumberTopHits(int numberTopHits) {
        this.numberTopHits = numberTopHits;
    }

    public void setFuzzFactor(double fuzzFactor) {
        this.fuzzFactor = fuzzFactor;
    }

    public File getExpressionDataSetFile() {
        return expressionDataSetFile;
    }

    public void setExpressionDataSetFile(File expressionDataSetFile) {
        this.expressionDataSetFile = expressionDataSetFile;
    }

    public File getSamplesOfInterestFile() {
        return samplesOfInterestFile;
    }

    public void setSamplesOfInterestFile(File samplesOfInterestFile) {
        this.samplesOfInterestFile = samplesOfInterestFile;
    }

    public File getNodeXrefToAffyIDMappingFile() {
        return nodeXrefToAffyIDMappingFile;
    }

    public void setNodeXrefToAffyIDMappingFile(File nodeXrefToAffyIDMappingFile) {
        this.nodeXrefToAffyIDMappingFile = nodeXrefToAffyIDMappingFile;
    }

    /**
     * Checks to make sure that the minimum information required is available
     *
     * @return true if the parameter set is valid, false if not
     */
    public boolean validateDBOptions() {
        if (databaseFileName == null) {
            return false;
        }
        if (databaseFormat == null) {
            return false;
        }
        if (searchParams == null) {
            return false;
        } else {
            if (searchParams.getTerminus() == null) {
                return false;
            }
        }
        return true;
    }

    public boolean validateProfileOptions() {
        if (profileFileName == null) {
            return false;
        }
        return true;
    }

    public String toString() {
        String lineSep = System.getProperty("line.separator");
        StringBuffer sb = new StringBuffer();
        if (profileSearch) {
            sb.append("Profile search threshold " + scoreThreshold +
                    ((searchParams.isNormalized()) ? " (Normalized)" : " (Not Normalized)") + lineSep);
            sb.append("Profile Fuzz Factor: " + fuzzFactor + lineSep);
            sb.append("Top " + scorePercentageThreshold + " percent of hits." + lineSep);
            sb.append("Codon Bias File: " + ((this.codonBiasFileName == null) ? "None" : this.codonBiasFileName));
            sb.append("Unique Profile Peptides: " + ((this.uniquePeptides) ? "Yes" : "No"));
            sb.append("Node Representation: " + ((this.uniqueQueryProteinNodes) ? "Proteins" : "Domains"));
        }
        sb.append(searchParams.toString());
        if ((databaseFormat != null) && (databaseFileName != null)) {
            sb.append(databaseFormat + " database: " + databaseFileName.toString() + lineSep);
        } else {
            if (databaseFormat == null) {
                sb.append("Database format not defined." + lineSep);
            }
            if (databaseFileName == null) {
                sb.append("Database not defined." + lineSep);
            }
        }
        return sb.toString();
    }
}
