package org.baderlab.csplugins.brain;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.db.HashSequenceDB;
import org.biojava.bio.seq.db.IllegalIDException;
import org.biojava.bio.seq.db.SequenceDB;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.utils.ChangeVetoException;
import static org.baderlab.csplugins.brain.BindingPeptideList.fileFormat.*;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;
import java.util.ArrayList;

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
 * * Time: 9:49:34 PM
 * * Description Stores a list of binding peptides, as may be derived from phage display
 */

/**
 * Stores a list of binding peptides, as may be derived from phage display, and associated annotation
 */
public class BindingPeptideList {
    private HashMap nameToPeptide = null; //key is peptide name, value is peptide
    private HashMap nameToFreqMap = null; //key is peptide name, value is frequency
    private HashMap nameToQuantMap = null;//key is peptide name, value is quantitation data
    private HashMap nameToExtIdMap = null;//key is peptide name, value is external database identifier
    private DatabaseReference protXref = null;
    private ArrayList protXrefList = new ArrayList();    // list of protein identifiers from Accession field
    private String proteinName = null;
    private int domainNumber = 0;
    private DatabaseReference domainXref = null;
    private String domainName = null;
    private String experimentalMethod = null;
    private String domainSequence = null;
    private int domainRangeStart = 0;
    private int domainRangeStop = 0;
    private String comment = null;
    private int minPeptideLength = 0;
    private int maxPeptideLength = 0;
    private String organism;
    private int taxid;
    private HashMap uniquePeptide = null; //key is peptide, value is peptide - used for determining unique peptides

    // peptide file formats
    enum fileFormat { one_dot_zero, one_dot_one }

    /**
     * Create a new binding peptide list object
     */
    public BindingPeptideList() {
        nameToPeptide = new HashMap();
        nameToFreqMap = new HashMap();
        nameToQuantMap = new HashMap();
        nameToExtIdMap = new HashMap();
        minPeptideLength = Integer.MAX_VALUE;
        maxPeptideLength = 0;
    }

    /**
     * Get the minimum peptide length in the list
     */
    public int getMinPeptideLength() {
        return minPeptideLength;
    }

    /**
     * Get the maximum peptide length in the list
     */
    public int getMaxPeptideLength() {
        return maxPeptideLength;
    }

    public DatabaseReference getProteinXref() {
        return protXref;
    }

    public ArrayList getProteinXrefList() {
        return protXrefList;
    }

    public String getProteinName() {
        return proteinName;
    }

    public int getDomainNumber() {
        return domainNumber;
    }

    public DatabaseReference getDomainXref() {
        return domainXref;
    }

    public String getDomainName() {
        return domainName;
    }

    public String getExperimentalMethod() {
        return experimentalMethod;
    }

    public String getDomainSequence() {
        return domainSequence;
    }

    public int getDomainRangeStart() {
        return domainRangeStart;
    }

    public int getDomainRangeStop() {
        return domainRangeStop;
    }

    public String getOrganism() {
        return organism;
    }

    public int getTaxid() {
        return taxid;
    }

    public String getComment() {
        return comment;
    }

    /**
     * Add peptideSequence to map
     *
     * @param peptideName     The unique name of the peptide
     * @param peptideSequence The peptide
     * @throws IllegalSymbolException if the peptide contains an illegal character
     */
    private void addPeptideToMap(String peptideName, String peptideSequence) throws IllegalSymbolException {
        if ((peptideSequence == null) || (peptideName == null)) {
            throw new IllegalArgumentException("addPeptideToMap failed");
        }
        //check if peptides contain an X
        //currently commented out because X and - are useful.  May want to have this as an option in the future.
        //if ((peptideSequence.indexOf("X") >= 0) || (peptideSequence.indexOf("-") >= 0)) {
        //    System.out.println("Warning: peptide contains an X or a - character. It will not be added to the profile. (" + peptideSequence + ")");
        //    return;
        //}
        if (nameToPeptide.containsKey(peptideName)) {
            //peptideName is unique, so we're not supposed to see another peptide with the same name
            System.out.println("Warning: duplicate peptide name found. It will not be added to the profile. (" + peptideName + ")");
        } else {
            nameToPeptide.put(peptideName, ProteinTools.createProteinSequence(peptideSequence, peptideName));
        }
    }

    /**
     * Add peptideSequence to map
     *
     * @param peptideName     The unique name of the peptide
     * @param peptideSequence The peptide
     * @param frequency       The frequency of this peptide
     * @throws IllegalSymbolException if the peptide contains an illegal character
     */
    private void addPeptideToMap(String peptideName, String peptideSequence, int frequency) throws IllegalSymbolException {
        if ((peptideSequence == null) || (frequency < 1) || (peptideName == null)) {
            throw new IllegalArgumentException("addPeptideToMap failed");
        }
        //check if peptides contain an X
        //currently commented out because X and - are useful.  May want to have this as an option in the future.
        //if ((peptideSequence.indexOf("X") >= 0) || (peptideSequence.indexOf("-") >= 0)) {
        //    System.out.println("Warning: peptide contains an X or a - character. It will not be added to the profile. (" + peptideSequence + ")");
        //    return;
        //}
        if (nameToPeptide.containsKey(peptideName)) {
            //peptideName is unique, so we're not supposed to see another peptide with the same name
            System.out.println("Warning: duplicate peptide name found. It will not be added to the profile. (" + peptideName + ")");
        } else {
            nameToPeptide.put(peptideName, ProteinTools.createProteinSequence(peptideSequence, peptideName));
            nameToFreqMap.put(peptideName, new Integer(frequency));
        }
    }

    /**
     * Add peptideSequence to map
     *
     * @param peptideName     The unique name of the peptide
     * @param peptideSequence The peptide
     * @param quantData       The quantitation data value of this peptide
     * @throws IllegalSymbolException if the peptide contains an illegal character
     */
    private void addPeptideToMap(String peptideName, String peptideSequence, double quantData) throws IllegalSymbolException {
        if ((peptideSequence == null) || (peptideName == null)) {
            throw new IllegalArgumentException("addPeptideToMap failed");
        }
        //check if peptides contain an X
        //currently commented out because X and - are useful.  May want to have this as an option in the future.
        //if ((peptideSequence.indexOf("X") >= 0) || (peptideSequence.indexOf("-") >= 0)) {
        //    System.out.println("Warning: peptide contains an X or a - character. It will not be added to the profile. (" + peptideSequence + ")");
        //    return;
        //}
        if (nameToPeptide.containsKey(peptideName)) {
            //peptideName is unique, so we're not supposed to see another peptide with the same name
            System.out.println("Warning: duplicate peptide name found. It will not be added to the profile. (" + peptideName + ")");
        } else {
            nameToPeptide.put(peptideName, ProteinTools.createProteinSequence(peptideSequence, peptideName));
            nameToQuantMap.put(peptideName, new Double(quantData));
        }
    }

    /**
     * Add peptideSequence to map
     *
     * @param peptideName     The unique name of the peptide
     * @param peptideSequence The peptide
     * @param externalId      The external database identifier of this peptide
     * @throws IllegalSymbolException if the peptide contains an illegal character
     */
    private void addPeptideToMap(String peptideName, String peptideSequence, String externalId) throws IllegalSymbolException {
        if ((peptideSequence == null) || (peptideName == null)) {
            throw new IllegalArgumentException("addPeptideToMap failed");
        }
        //check if peptides contain an X
        //currently commented out because X and - are useful.  May want to have this as an option in the future.
        //if ((peptideSequence.indexOf("X") >= 0) || (peptideSequence.indexOf("-") >= 0)) {
        //    System.out.println("Warning: peptide contains an X or a - character. It will not be added to the profile. (" + peptideSequence + ")");
        //    return;
        //}
        if (nameToPeptide.containsKey(peptideName)) {
            //peptideName is unique, so we're not supposed to see another peptide with the same name
            System.out.println("Warning: duplicate peptide name found. It will not be added to the profile. (" + peptideName + ")");
        } else {
            nameToPeptide.put(peptideName, ProteinTools.createProteinSequence(peptideSequence, peptideName));
            nameToExtIdMap.put(peptideName, externalId);
        }
    }

    /**
     * Add peptideSequence to map
     *
     * @param peptideName     The unique name of the peptide
     * @param peptideSequence The peptide
     * @param frequency       The frequency of this peptide
     * @param quantData       The quantitation data value of this peptide
     * @throws IllegalSymbolException if the peptide contains an illegal character
     */
    private void addPeptideToMap(String peptideName, String peptideSequence, int frequency, double quantData) throws IllegalSymbolException {
        if ((peptideSequence == null) || (frequency < 1) || (peptideName == null)) {
            throw new IllegalArgumentException("addPeptideToMap failed");
        }
        //check if peptides contain an X
        //currently commented out because X and - are useful.  May want to have this as an option in the future.
        //if ((peptideSequence.indexOf("X") >= 0) || (peptideSequence.indexOf("-") >= 0)) {
        //    System.out.println("Warning: peptide contains an X or a - character. It will not be added to the profile. (" + peptideSequence + ")");
        //    return;
        //}
        if (nameToPeptide.containsKey(peptideName)) {
            //peptideName is unique, so we're not supposed to see another peptide with the same name
            System.out.println("Warning: duplicate peptide name found. It will not be added to the profile. (" + peptideName + ")");
        } else {
            nameToPeptide.put(peptideName, ProteinTools.createProteinSequence(peptideSequence, peptideName));
            nameToFreqMap.put(peptideName, new Integer(frequency));
            nameToQuantMap.put(peptideName, new Double(quantData));
        }
    }

    /**
     * Add peptideSequence to map
     *
     * @param peptideName     The unique name of the peptide
     * @param peptideSequence The peptide
     * @param frequency       The frequency of this peptide
     * @param externalId      The external database identifier of this peptide
     * @throws IllegalSymbolException if the peptide contains an illegal character
     */
    private void addPeptideToMap(String peptideName, String peptideSequence, int frequency, String externalId) throws IllegalSymbolException {
        if ((peptideSequence == null) || (frequency < 1) || (peptideName == null)) {
            throw new IllegalArgumentException("addPeptideToMap failed");
        }
        //check if peptides contain an X
        //currently commented out because X and - are useful.  May want to have this as an option in the future.
        //if ((peptideSequence.indexOf("X") >= 0) || (peptideSequence.indexOf("-") >= 0)) {
        //    System.out.println("Warning: peptide contains an X or a - character. It will not be added to the profile. (" + peptideSequence + ")");
        //    return;
        //}
        if (nameToPeptide.containsKey(peptideName)) {
            //peptideName is unique, so we're not supposed to see another peptide with the same name
            System.out.println("Warning: duplicate peptide name found. It will not be added to the profile. (" + peptideName + ")");
        } else {
            nameToPeptide.put(peptideName, ProteinTools.createProteinSequence(peptideSequence, peptideName));
            nameToFreqMap.put(peptideName, new Integer(frequency));
            nameToExtIdMap.put(peptideName, new String(externalId));
        }
    }

    /**
     * Add peptideSequence to map
     *
     * @param peptideName     The unique name of the peptide
     * @param peptideSequence The peptide
     * @param quantData       The quantitation data value of this peptide
     * @param externalId      The external database identifier of this peptide
     * @throws IllegalSymbolException if the peptide contains an illegal character
     */
    private void addPeptideToMap(String peptideName, String peptideSequence, double quantData, String externalId) throws IllegalSymbolException {
        if ((peptideSequence == null) || (peptideName == null)) {
            throw new IllegalArgumentException("addPeptideToMap failed");
        }
        //check if peptides contain an X
        //currently commented out because X and - are useful.  May want to have this as an option in the future.
        //if ((peptideSequence.indexOf("X") >= 0) || (peptideSequence.indexOf("-") >= 0)) {
        //    System.out.println("Warning: peptide contains an X or a - character. It will not be added to the profile. (" + peptideSequence + ")");
        //    return;
        //}
        if (nameToPeptide.containsKey(peptideName)) {
            //peptideName is unique, so we're not supposed to see another peptide with the same name
            System.out.println("Warning: duplicate peptide name found. It will not be added to the profile. (" + peptideName + ")");
        } else {
            nameToPeptide.put(peptideName, ProteinTools.createProteinSequence(peptideSequence, peptideName));
            nameToQuantMap.put(peptideName, new Double(quantData));
            nameToExtIdMap.put(peptideName, new String(externalId));
        }
    }

    /**
     * Add peptideSequence to map
     *
     * @param peptideName     The unique name of the peptide
     * @param peptideSequence The peptide
     * @param frequency       The frequency of this peptide
     * @param quantData       The quantitation data value of this peptide
     * @param externalId      The external database identifier of this peptide
     * @throws IllegalSymbolException if the peptide contains an illegal character
     */
    private void addPeptideToMap(String peptideName, String peptideSequence, int frequency, double quantData, String externalId) throws IllegalSymbolException {
        if ((peptideSequence == null) || (frequency < 1) || (peptideName == null)) {
            throw new IllegalArgumentException("addPeptideToMap failed");
        }
        //check if peptides contain an X
        //currently commented out because X and - are useful.  May want to have this as an option in the future.
        //if ((peptideSequence.indexOf("X") >= 0) || (peptideSequence.indexOf("-") >= 0)) {
        //    System.out.println("Warning: peptide contains an X or a - character. It will not be added to the profile. (" + peptideSequence + ")");
        //    return;
        //}
        if (nameToPeptide.containsKey(peptideName)) {
            //peptideName is unique, so we're not supposed to see another peptide with the same name
            System.out.println("Warning: duplicate peptide name found. It will not be added to the profile. (" + peptideName + ")");
        } else {
            nameToPeptide.put(peptideName, ProteinTools.createProteinSequence(peptideSequence, peptideName));
            nameToFreqMap.put(peptideName, new Integer(frequency));
            nameToQuantMap.put(peptideName, new Double(quantData));
            nameToExtIdMap.put(peptideName, new String(externalId));
        }
    }

    /**
     * Detect the format of a peptide file. Accepted formats are represented in the 'fileFormat' enum.
     * The format is determined by parsing the file and evaluating the size of the header section
     * and the number of columns in the peptide section, as each format is distinct in these attributes.
     *
     * @param fileName The filename of the peptide list
     * @return true if the read is successful
     */
    private fileFormat detectFileFormat(String fileName) throws IOException, IllegalSymbolException {
        int lineCount = 0, headerSize = 0, sequenceColumns = 0;
        BufferedReader br = new BufferedReader(new FileReader(fileName));
        String line = null;
        String[] lineSplit;

        while ((line = br.readLine()) != null) {

            lineCount++;

            //ignore comment or blank line
            if ((line.equals("")) || (line.charAt(0) == '#')) {
                continue;
            }

            //all lines in file (except comments/blanks) must be tab-delimited, so split the columns
            lineSplit = line.split("[\t]");

            //if <2 columns, format is invalid
            if (lineSplit.length < 2) {
                throw new RuntimeException("Expected two tab-separated columns. Unable to parse header at line " + lineCount +
                        " (Found: " + line + "). Please check the input file: " + fileName);
            }

            //if 2 columns, we're in the header
            if (lineSplit.length == 2) {
                ++headerSize;
            }
            //if >2 columns, assume we've reached sequence section
            else {
                sequenceColumns = lineSplit.length;
                break;
            }
        }
        //determine format
        if ( (headerSize == 9) && (sequenceColumns == 3) ) {
            return one_dot_zero;
        }
        else if ( (headerSize == 11) && (sequenceColumns == 5) ) {
            return one_dot_one;            
        }
        else {
            throw new RuntimeException("Unable to determine peptide file format version. File format may be invalid." +
                    " Please check the input file: " + fileName);
        }
    }

    /**
     * Reads a formatted text file (peptide file format version 1.1)
     * containing a peptide list and associated annotation.
     * <p/>
     * The current format is:
     * <pre>
     * Gene Name       NP_008944_D00001
     * Accession       Refseq:NP_008944 Ensembl:ENSG00000123124
     * Organism        Homo Sapiens (Human)
     * NCBITaxonomyID  9606
     * Domain Number   0
     * Domain Type     WW
     * Interpro ID     IPR001202
     * Technique       Protein Chip
     * Domain Sequence ETLPSGWEQRKDPHGRTYYVDHNTRTTTWERPQPL
     * Domain Range    255-300
     * Comment This is a sample peptide file in the new format (file format version 1.1)
     * PeptideName     Peptide CloneFrequency  QuantData       ExternalIdentifer
     * 1       XXXXXEFCSPPAYATLTXX     1
     * 2       XXXXXLYGCPPPYHTFEXX     1
     * 3       XXXXXLCGYPPFYEETEXX             0.5
     * </pre>
     * <p/>
     * Peptides are uniquely numbered starting at 1.  Peptides that should not be considered are labeled
     * with single letters, alphabetically ordered starting at A. Tabs separate all fields.
     *
     * Columns <i>CloneFrequency</i>, <i>QuantData</i>, and <i>ExternalIdentifier</i> are optional and
     * their values can be left empty.
     *
     * @param fileName The filename of the peptide list
     * @return true if read is successful
     * @throws IOException
     * @throws IllegalSymbolException
     */
    public boolean read(String fileName) throws IOException, IllegalSymbolException {

        //call the appropriate reader based on the detected file format
        fileFormat format = detectFileFormat(fileName);
        switch (format) {
            case one_dot_one:
                return this.readFormat1_1(fileName);

            case one_dot_zero:
                return this.readFormat1_0(fileName);
        }

        //this should not happen - format should be one of the fileFormat enum values
        return false;
    }

    /**
     * Reads a peptide file formatted according to the 1.1 format version.
     * <p/>
     * For more information, see <a>BindingPeptideList.read</a>
     *
     * @param fileName
     * @return
     * @throws IOException
     * @throws IllegalSymbolException
     */
    private boolean readFormat1_1(String fileName) throws IOException, IllegalSymbolException {
        int lineCount = 0;
        String col1 = null;
        String col2 = null;
        String accession = null;
        boolean header = true;
        String peptideName = null;
        String peptideSequence = null;
        int peptideSequenceLength = 0;
        int frequency = 0;
        double quantData = 0;
        String externalId = null;
        BufferedReader br = new BufferedReader(new FileReader(fileName));
        String line = null;
        String[] lineSplit;
        boolean haveFreq = false;
        boolean haveQuant = false;
        boolean haveExtId = false;

        while ((line = br.readLine()) != null) {
            haveFreq = false;
            haveQuant = false;
            haveExtId = false;
            
            lineCount++;
            //check for a comment or a blank line
            if ((line.equals("")) || (line.charAt(0) == '#')) {
                continue;
            }
            lineSplit = line.split("[\t]",5);
            if (header) {
                if (lineSplit.length != 2) {
                    throw new RuntimeException("Expected two tab-separated columns. Unable to parse header at line " + lineCount +
                            " (Found: " + line + "). Please check the input file: " + fileName);
                }
                col1 = lineSplit[0];
                col2 = lineSplit[1];
                //assume proper ordering of header lines
                if (col1.equalsIgnoreCase("Gene Name")) {
                    proteinName = col2;
                } else if (col1.equalsIgnoreCase("Accession")) {
                    if (proteinName == null) {
                        throw new RuntimeException("Expected gene name before accession at line " + lineCount +
                                ". Please check the input file: " + fileName);
                    }
                    // split accession string to read multiple space-separated identifiers
                    String[] accessions;
                    accessions = col2.split(" ");
                    for (int i = 0; i < accessions.length; i++) {
                        accession = accessions[i];
                        if (accession.indexOf(":") < 0) {
                            throw new RuntimeException("Expected colon separated accession at line " + lineCount +
                                    " (Found: " + col2 + "). Please check the input file: " + fileName);
                        }
                        // store the first identifier in protXref; add it to protXrefList
                        if (protXref == null) {
                            protXref = new DatabaseReference(accession);
                            protXrefList.add(protXref);
                        }
                        else {
                            protXrefList.add(new DatabaseReference(accession));
                        }
                    }
                } else if (col1.equalsIgnoreCase("Organism")) {
                    organism = col2;
                } else if (col1.equalsIgnoreCase("NCBITaxonomyID")) {
                    try {
                        taxid = Integer.valueOf(col2).intValue();
                    } catch (NumberFormatException e) {
                        throw new RuntimeException("Unable to find a valid NCBI Taxonomy ID at line " + lineCount +
                                " (Found: " + lineSplit[1] + "). Please check the input file: " + fileName, e);
                    }
                } else if (col1.equalsIgnoreCase("Domain Number")) {
                    try {
                        domainNumber = Integer.valueOf(col2).intValue();
                    } catch (NumberFormatException e) {
                        throw new RuntimeException("Unable to find a valid domain number at line " + lineCount +
                                " (Found: " + lineSplit[1] + "). Please check the input file: " + fileName, e);
                    }
                } else if (col1.equalsIgnoreCase("Domain Type")) {
                    domainName = col2;
                } else if (col1.equalsIgnoreCase("Interpro ID")) {
                    if (proteinName == null) {
                        throw new RuntimeException("Expected domain type before Interpro ID at line " + lineCount +
                                ". Please check the input file: " + fileName);
                    }
                    domainXref = new DatabaseReference("Interpro", col2);
                } else if (col1.equalsIgnoreCase("Technique")) {
                    experimentalMethod = col2;
                } else if (col1.equalsIgnoreCase("Domain sequence")) {
                    domainSequence = col2;
                } else if (col1.equalsIgnoreCase("Domain Range")) {
                    // range is given as "start-stop" where start=starting aa position, stop=ending aa position
                    if (col2.indexOf("-") < 0) {
                        throw new RuntimeException("Expected dash-separated Domain Range (e.g. 255-300) " +
                                "at line " + lineCount + " (Found: " + lineSplit[1] +
                                "). Please check the input file: " + fileName);
                    }
                    String[] range = col2.split("-");
                    if (range.length != 2) {
                        throw new RuntimeException("Expected domain range format 'start-stop' at line " + lineCount +
                                " (Found: " + lineSplit[1] + "). Please check the input file: " + fileName);
                    }
                    try {
                        domainRangeStart = Integer.valueOf(range[0]);
                        domainRangeStop = Integer.valueOf(range[1]);
                    } catch (NumberFormatException e) {
                        throw new RuntimeException("Unable to find a valid domain range value at line " + lineCount +
                        " (Found: " + lineSplit[1] + "). Please check the input file: " + fileName, e);
                    }                    
                } else if (col1.equalsIgnoreCase("Comment")) {
                    // comment is optional so only assign if it's there
                    if (col2.length() >  0) {
                        comment = col2;
                    }

                    //we now have all of the header information we need, we can start parsing peptides now
                    header = false;

                    //skip past the peptide header line
                    line = br.readLine();
                    ++lineCount;

                    if (line == null) {
                        throw new RuntimeException("Expected peptide section header at line " + lineCount +
                                ". Please check the input file: " + fileName);
                    }
                }

            } else {
                if (lineSplit.length >= 2) {
                    //Assume that file is peptideName[tab]peptide[tab]frequency

                    // PeptideName and Peptide are required in first 2 columns, remaining values are optional
                    if (lineSplit[0].length() < 1) {
                        throw new RuntimeException("Expected a value in PeptideName column at line " + lineCount +
                        " (Found: " + lineSplit + "). Please check the input file: " + fileName);
                    }
                    if (lineSplit[1].length() < 1) {
                        throw new RuntimeException("Expected a value in Peptide column at line " + lineCount +
                        " (Found: " + lineSplit + "). Please check the input file: " + fileName);
                    }
                    peptideName = lineSplit[0];
                    peptideSequence = lineSplit[1];

                    // CloneFrequency - optional
                    if ((lineSplit.length > 2) && (lineSplit[2].length() > 0)) {
                        try {
                            frequency = Integer.parseInt(lineSplit[2]);
                        } catch (NumberFormatException e) {
                            throw new RuntimeException("Unable to find a valid frequency value at line " + lineCount +
                                    " (Found: " + lineSplit[2] + "). Please check the input file: " + fileName, e);
                        }
                        haveFreq = true;
                    }

                    // QuantData - optional
                    if ((lineSplit.length > 3) && (lineSplit[3].length() > 0)) {
                        try {
                            quantData = Double.parseDouble(lineSplit[3]);
                        } catch (NumberFormatException e) {
                            throw new RuntimeException("Unable to read a valid quantitation data value at line " +
                            lineCount + "(Found: " + lineSplit[3] + "). Please check the input file: " + fileName, e);
                        }
                        haveQuant = true;
                    }

                    // ExternalIdentifier - optional
                    if ((lineSplit.length > 4) && (lineSplit[4].length() > 0)) {
                        externalId = lineSplit[4];
                        haveExtId = true;                        
                    }

                    //Only add peptide if it is a number
                    try {
                        Integer.parseInt(peptideName);

                        // add peptide and associated optional parameters
                        if (haveFreq && haveQuant && haveExtId) {
                            addPeptideToMap(peptideName, peptideSequence, frequency, quantData, externalId);
                        }
                        else if (haveFreq && haveQuant) {
                            addPeptideToMap(peptideName, peptideSequence, frequency, quantData);
                        }
                        else if (haveFreq && haveExtId) {
                            addPeptideToMap(peptideName, peptideSequence, frequency, externalId);
                        }
                        else if (haveQuant && haveExtId) {
                            addPeptideToMap(peptideName, peptideSequence, quantData, externalId);
                        }
                        else if (haveFreq) {
                            addPeptideToMap(peptideName, peptideSequence, frequency);
                        }
                        else if (haveQuant) {
                            addPeptideToMap(peptideName, peptideSequence, quantData);
                        }
                        else if (haveExtId) {
                            addPeptideToMap(peptideName, peptideSequence, externalId);
                        }
                        else {
                            addPeptideToMap(peptideName, peptideSequence);
                        }


                        //addPeptideToMap(peptideName, peptideSequence, frequency);
                        //keep track of peptide length
                        peptideSequenceLength = peptideSequence.length();
                        if (minPeptideLength > peptideSequenceLength) {
                            minPeptideLength = peptideSequenceLength;
                        }
                        if (maxPeptideLength < peptideSequenceLength) {
                            maxPeptideLength = peptideSequenceLength;
                        }
                    } catch (NumberFormatException e) {
                        //don't add if a number is not detected
                        continue;
                    }
                } else {
                    //error with this line, it doesn't have 3 columns
                    throw new RuntimeException("Line number: " + lineCount + " is not valid. (" + line + ")");
                }
            }
        }
        br.close();
        return true;
    }

    /**
     * Reads a formatted text file containing a peptide list and associated annotation.
     * </p>
     * The current format is:
     * <pre>
     * Gene Name	BAIAP1
     * Accession	Refseq:NP_004733
     * Organism	Homo Sapiens (Human)
     * NCBITaxonomyID	9606
     * Domain Number	2
     * Domain Type	PDZ
     * Interpro ID	IPR001478
     * Technique	Phage Display High Valency
     * Domain sequence	TVHIVKGPMGFGFTIADSPGGGGQRVKQIVDSPRCRGLKEGDLIVEVNKKNVQALTHNQVVDMLVECPKGSEVTLLV
     * PeptideName	Peptide	CloneFrequency
     * 1	LKHVRQTWDL	1
     * 2	ILRRGRETLL	1
     * A	XXXXSDFYFC	1
     * </pre>
     * <p/>
     * Peptides are uniquely numbered starting at 1.  Peptides that should not be considered are labeled
     * with single letters, alphabetically ordered starting at A. Tabs separate all fields.
     *
     * @param fileName The filename of the peptide list
     * @return true if read is successful
     */
    private boolean readFormat1_0(String fileName) throws IOException, IllegalSymbolException {
        int lineCount = 0;
        String col1 = null;
        String col2 = null;
        String accession = null;
        boolean header = true;
        String peptideName = null;
        String peptideSequence = null;
        int peptideSequenceLength = 0;
        int frequency = 0;
        BufferedReader br = new BufferedReader(new FileReader(fileName));
        String line = null;
        String[] lineSplit;
        while ((line = br.readLine()) != null) {
            lineCount++;
            //check for a comment or a blank line
            if ((line.equals("")) || (line.charAt(0) == '#')) {
                continue;
            }
            lineSplit = line.split("[\t]");
            if (header) {
                if (lineSplit.length != 2) {
                    throw new RuntimeException("Expected two tab-separated columns. Unable to parse header at line " + lineCount +
                            " (Found: " + line + "). Please check the input file: " + fileName);
                }
                col1 = lineSplit[0];
                col2 = lineSplit[1];
                //assume proper ordering of header lines
                if (col1.equalsIgnoreCase("Gene Name")) {
                    proteinName = col2;
                } else if (col1.equalsIgnoreCase("Accession")) {
                    if (proteinName == null) {
                        throw new RuntimeException("Expected gene name before accession at line " + lineCount +
                                ". Please check the input file: " + fileName);
                    }
                    accession = col2;
                    protXref = new DatabaseReference(accession);
                } else if (col1.equalsIgnoreCase("Organism")) {
                    organism = col2;
                } else if (col1.equalsIgnoreCase("NCBITaxonomyID")) {
                    try {
                        taxid = Integer.valueOf(col2).intValue();
                    } catch (NumberFormatException e) {
                        throw new RuntimeException("Unable to find a valid NCBI Taxonomy ID at line " + lineCount +
                                " (Found: " + lineSplit[1] + "). Please check the input file: " + fileName, e);
                    }
                } else if (col1.equalsIgnoreCase("Domain Number")) {
                    try {
                        domainNumber = Integer.valueOf(col2).intValue();
                    } catch (NumberFormatException e) {
                        throw new RuntimeException("Unable to find a valid domain number at line " + lineCount +
                                " (Found: " + lineSplit[1] + "). Please check the input file: " + fileName, e);
                    }
                } else if (col1.equalsIgnoreCase("Domain Type")) {
                    domainName = col2;
                } else if (col1.equalsIgnoreCase("Interpro ID")) {
                    if (proteinName == null) {
                        throw new RuntimeException("Expected domain type before Interpro ID at line " + lineCount +
                                ". Please check the input file: " + fileName);
                    }
                    domainXref = new DatabaseReference("Interpro", col2);
                } else if (col1.equalsIgnoreCase("Technique")) {
                    experimentalMethod = col2;
                } else if (col1.equalsIgnoreCase("Domain sequence")) {
                    domainSequence = col2;
                    //we now have all of the header information we need, we can start parsing peptides now
                    header = false;

                    //skip past the peptide header line
                    line = br.readLine();
                    ++lineCount;
                    
                    if (line == null) {
                        throw new RuntimeException("Expected peptide section header at line " + lineCount +
                                ". Please check the input file: " + fileName);
                    }
                }
            } else {
                if (lineSplit.length == 3) {
                    //Assume that file is peptideName[tab]peptide[tab]frequency
                    peptideName = lineSplit[0];
                    peptideSequence = lineSplit[1];
                    try {
                        frequency = Integer.parseInt(lineSplit[2]);
                    } catch (NumberFormatException e) {
                        throw new RuntimeException("Unable to find a valid frequency value at line " + lineCount +
                                " (Found: " + lineSplit[2] + "). Please check the input file: " + fileName, e);
                    }
                    //Only add peptide if it is a number
                    try {
                        Integer.parseInt(peptideName);
                        addPeptideToMap(peptideName, peptideSequence, frequency);
                        //keep track of peptide length
                        peptideSequenceLength = peptideSequence.length();
                        if (minPeptideLength > peptideSequenceLength) {
                            minPeptideLength = peptideSequenceLength;
                        }
                        if (maxPeptideLength < peptideSequenceLength) {
                            maxPeptideLength = peptideSequenceLength;
                        }
                    } catch (NumberFormatException e) {
                        //don't add if a number is not detected
                        continue;
                    }
                } else {
                    //error with this line, it doesn't have 3 columns
                    throw new RuntimeException("Line number: " + lineCount + " is not valid. (" + line + ")");
                }
            }
        }
        br.close();
        return true;
    }

    /**
     * Tests whether a given peptide is unique (i.e. has not already been seen)
     *
     * @return true if the peptide is unique
     */
    private boolean isUniquePeptide(String peptide) {
        if (uniquePeptide == null) {
            uniquePeptide = new HashMap();
        }
        if (uniquePeptide.containsKey(peptide)) {
            return false;
        }
        uniquePeptide.put(peptide, peptide);
        return true;
    }

    /**
     * Get sequence iterator over peptides of a specific length
     *
     * @param length   Length of returned peptides. If peptides are shorter than this length, discard, if longer, truncate
     * @param terminus The terminus to select from (N or C).
     * @return an iterator over the selected peptides from this list
     */
    public SequenceIterator getSequenceIteratorByLength(int length, ProteinTerminus terminus) {
        return (getSequenceIteratorByLength(length, terminus, false, false));
    }

    /**
     * Get sequence iterator over peptides of a specific length
     *
     * @param length             Length of returned peptides. If peptides are shorter than this length, discard, if longer, truncate
     * @param terminus           The terminus to select from (N or C).
     * @param frequencyExpansion If true, and a frequency number exists for the peptides in the list,
     *                           then the frequency number of peptides will be added to the iterator
     * @param unique             If true, only returns unique peptides (even if cut to a shorter length).
     *                           It does not make sense to set both unique and frequencyExpansion to true.
     * @return an iterator over the selected peptides from this list
     */
    public SequenceIterator getSequenceIteratorByLength(int length, ProteinTerminus terminus, boolean frequencyExpansion, boolean unique) {
        //check if length makes sense for this peptide list
        if (length < 1 || length > maxPeptideLength) {
            throw new IllegalArgumentException("Length " + length + " is invalid. Must be between 1 and " + maxPeptideLength + ".");
        }
        //check if terminus is appropriate - must be given if there is a need to truncate
        if (length < maxPeptideLength && !(terminus == ProteinTerminus.C || terminus == ProteinTerminus.N)) {
            throw new IllegalArgumentException("You must provide a terminus if peptides will be truncated.");
        }
        SequenceDB searchDB = new HashSequenceDB();
        Set peptideNames = nameToPeptide.keySet();
        SymbolList peptideToAdd = null;
        for (Iterator iterator = peptideNames.iterator(); iterator.hasNext();) {
            String peptideName = (String) iterator.next();
            Sequence peptide = (Sequence) nameToPeptide.get(peptideName);
            //check if this passes the length criteria
            int peptideLength = peptide.length();
            if (peptideLength > length) {
                //truncate based on terminus
                if (terminus == ProteinTerminus.C) {
                    //get end of string
                    peptideToAdd = peptide.subList(peptideLength - length + 1, peptideLength);
                } else if (terminus == ProteinTerminus.N) {
                    //get start of string
                    peptideToAdd = peptide.subList(1, length);
                }
            } else if (peptideLength == length) {
                peptideToAdd = peptide;
            } else {
                //don't add smaller peptides
                peptideToAdd = null;
            }
            try {
                if (peptideToAdd != null) {
                    if (frequencyExpansion) {
                        int frequency = ((Integer) nameToFreqMap.get(peptideName)).intValue();
                        for (int i = 0; i < frequency; i++) {
                            searchDB.addSequence(ProteinTools.createProteinSequence(peptideToAdd.seqString(), peptide.getName() + "_" + i));
                        }
                    } else {
                        if (unique) {
                            if (isUniquePeptide(peptideToAdd.seqString())) {
                                searchDB.addSequence(ProteinTools.createProteinSequence(peptideToAdd.seqString(), peptide.getName()));
                            } else {
                                //uncomment to print out filtered peptides
                                //System.out.println(proteinName  + "-" + domainNumber + " " +  peptide.getName() + " " + peptideToAdd.seqString());
                            }
                        } else {
                            searchDB.addSequence(ProteinTools.createProteinSequence(peptideToAdd.seqString(), peptide.getName()));
                        }
                    }
                }
            } catch (ChangeVetoException e) {
                e.printStackTrace();  //This should never happen
            } catch (IllegalIDException e) {
                e.printStackTrace();
            } catch (BioException e) {
                e.printStackTrace();
            }
        }
        return (searchDB.sequenceIterator());
    }

    /**
     * Get sequence iterator over all peptides in the list. No frequency expansion or length cutting performed.
     *
     * @return an iterator over the selected peptides from this list
     */
    public SequenceIterator getSequenceIterator() {
        return getSequenceIteratorByLength(getMaxPeptideLength(), ProteinTerminus.NONE, false, false);
    }

    /**
     * Get sequence iterator over all peptides in the list. Specified frequency expansion performed. No length cutting performed.
     *
     * @return an iterator over the selected peptides from this list
     */
    public SequenceIterator getSequenceIterator(boolean frequencyExpansion) {
        return getSequenceIteratorByLength(getMaxPeptideLength(), ProteinTerminus.NONE, frequencyExpansion, false);
    }

    /**
     * Get sequence iterator over all peptides in the list. Specified frequency expansion performed. No length cutting performed.
     * Peptides are uniqued, if specified.
     *
     * @return an iterator over the selected peptides from this list
     */
    public SequenceIterator getSequenceIterator(boolean frequencyExpansion, boolean unique) {
        return getSequenceIteratorByLength(getMaxPeptideLength(), ProteinTerminus.NONE, frequencyExpansion, unique);
    }

}
