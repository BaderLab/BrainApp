package org.baderlab.csplugins.brainplugin;

import cytoscape.CyEdge;
import cytoscape.CyNetwork;
import cytoscape.CyNode;
import cytoscape.data.CyAttributes;
import cytoscape.data.CyAttributesImpl;
import org.biojava.bio.seq.Sequence;
import org.baderlab.csplugins.brainplugin.inparanoid.InparanoidDB;
import org.baderlab.brain.BrainAlgorithm;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

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
 * * Date: Sep 1, 2005
 * * Time: 7:37:48 PM
 */

/**
 * Analysis of domain interactions that can be used for protein-protein interaction prediction
 */
public class DomainInteractionAnalysis {

    /**
     * Store all data for a single GEO sample
     */
    private static class GEODataSet {
        private String sampleID;
        private String tissueName;
        private HashMap affyID2APCall;

        public GEODataSet(String sampleID, String tissueName) {
            this.sampleID = sampleID;
            this.tissueName = tissueName;
            affyID2APCall = new HashMap(5000);
        }

        public String getSampleID() {
            return sampleID;
        }

        public void setSampleID(String sampleID) {
            this.sampleID = sampleID;
        }

        public String getTissueName() {
            return tissueName;
        }

        public void setTissueName(String tissueName) {
            this.tissueName = tissueName;
        }

        public void addAffyID2APCallMapping(String affyID, String APCall) {
            affyID2APCall.put(affyID, APCall);
        }

        public String getAffyID2APCallMapping(String affyID) {
            return ((String) affyID2APCall.get(affyID));
        }
    }

    /**
     * Scores a network by coexpression. After this method is called, each edge will contain
     * a count of the number of times out of all samples of interest the two nodes of the edge
     * are both present, given an Affymetrix absent/present call.
     *
     * @param expressionDataSetFile       The file of SOFT Sample formatted data for the dataset (usually
     *                                    multiple samples per file)
     * @param samplesOfInterestFile       A tab or space delimited list of the samples of interest
     *                                    1st column: GSM ID (GEO Sample ID)
     *                                    2nd column: Sample description string - usually a tissue name
     * @param nodeXrefToAffyIDMappingFile A tab delimited file mapping the node IDs to Affymetrix IDs
     *                                    (or more generally, the IDs in the GSM file)
     * @param net                         The network to score
     * @throws IOException if one of the files cannot be read for some reason.
     */
    public static void scoreCoExpressionByAPCall(File expressionDataSetFile, File samplesOfInterestFile,
                                                 File nodeXrefToAffyIDMappingFile, CyNetwork net) throws IOException {
        //read samples of interest so we know which samples to process
        /*File of the form:
        GSM2820\tBrain, cerebral cortex
        GSM2821\tUmbilical vein endothelial cells
        GSM2822\tLung
        */
        int samplesOfInterestCount = 0;
        HashMap samplesOfInterest2TissueNameMap = new HashMap();
        BufferedReader brInterest = new BufferedReader(new FileReader(samplesOfInterestFile));
        String samplesOfInterestFileLine;
        while ((samplesOfInterestFileLine = brInterest.readLine()) != null) {
            //check for a comment or empty line
            if ((samplesOfInterestFileLine.equals("")) || (samplesOfInterestFileLine.charAt(0) == '#')) {
                continue;
            }
            String[] gsmString = samplesOfInterestFileLine.split("\t", 2);
            if (gsmString.length == 2) {
                samplesOfInterest2TissueNameMap.put(new String(gsmString[0]), new String(gsmString[1]));
                samplesOfInterestCount++;
            } else {
                throw new IllegalArgumentException("Sample selection list line was malformed. " +
                        "Expected two tab delimited columns, but got " + samplesOfInterestFileLine + ".");
            }
        }

        //parse GNF soft file for absent/present call - move into HashMap mapping affy ID to A/M/P call
        // - one hashmap for each sample
        GEODataSet[] samples = new GEODataSet[samplesOfInterestCount];
        BufferedReader brSamples = new BufferedReader(new FileReader(expressionDataSetFile));
        String currentSample = null;
        int currentSampleIndex = -1;
        String sampleFileLine;
        boolean interestingSample = false;
        while ((sampleFileLine = brSamples.readLine()) != null) {
            /* SOFT rules:
            ^	caret lines	entity indicator line
            !	bang lines	entity attribute line
            #	hash lines	data table header description line
            n/a	data lines	data table row
            */
            char firstChar = sampleFileLine.charAt(0);
            if ((interestingSample) && (firstChar != '!') && (firstChar != '#') && (firstChar != '^')) {
                //this is a data line, parse it.  Assume that currentSample != null, as per SOFT format order
                //data line looks like this: 1001_at	4	4	16	16	14	0.25	0.23	0	0	1.00	21.46	A
                String[] affyIDString = sampleFileLine.split("\\s", 2);
                String affyID = new String(affyIDString[0]);
                if (affyID.equals("ID_REF")) {
                    continue; //skip the data header line
                }
                //TODO: generalize this by passing the column number of the AP call - not necessarily the last character
                char APCall = sampleFileLine.charAt(sampleFileLine.length() - 1);
                samples[currentSampleIndex].addAffyID2APCallMapping(affyID, String.valueOf(APCall));
            } else if (sampleFileLine.charAt(0) == '^') {
                if (sampleFileLine.indexOf("^SAMPLE") >= 0) {
                    //Parse: "^SAMPLE = GSM2819"
                    String[] gsmString = sampleFileLine.split(" = ", 2);
                    currentSample = new String(gsmString[1]);
                    if (samplesOfInterest2TissueNameMap.containsKey(currentSample)) {
                        currentSampleIndex++;
                        samples[currentSampleIndex] = new GEODataSet(currentSample, (String) samplesOfInterest2TissueNameMap.get(currentSample));
                        interestingSample = true;
                    } else {
                        interestingSample = false;
                    }
                }
            }
        }
        if ((currentSampleIndex + 1) != samplesOfInterestCount) {
            throw new IllegalArgumentException("Error detected in samples of interest file. It specifies " +
                    samplesOfInterestCount + " samples of interest, but only found " + (currentSampleIndex + 1) + ".");
        }

        //parse mapping file of IDs of interest to affy IDs
        HashMap nodeIDToAffyIDList = new HashMap();
        BufferedReader brIDMappingFile = new BufferedReader(new FileReader(nodeXrefToAffyIDMappingFile));
        String IDMappingFileLine;
        //line format is expected to be "DBNAME:DBID\tAffyID"
        //tolerate duplicate node or affy IDs
        while ((IDMappingFileLine = brIDMappingFile.readLine()) != null) {
            String[] mappingString = IDMappingFileLine.split("\\s", 2);
            if (mappingString.length == 2) {
                DatabaseReference dbref = new DatabaseReference(mappingString[0]);
                if (nodeIDToAffyIDList.containsKey(dbref)) {
                    ArrayList affyList = (ArrayList) nodeIDToAffyIDList.get(dbref);
                    affyList.add(new String(mappingString[1]));
                } else {
                    ArrayList affyList = new ArrayList();
                    affyList.add(new String(mappingString[1]));
                    nodeIDToAffyIDList.put(dbref, affyList);
                }
            }
        }

        //iterate through interactions - add 1 to an edge attribute each time both genes are present in the same tissue
        Iterator edges = net.edgesIterator();
        while (edges.hasNext()) {
            CyEdge cyEdge = (CyEdge) edges.next();
            CountNodeCoExpression(net, cyEdge, nodeIDToAffyIDList, samples);
        }
    }

    /**
     * Counts the number of times two nodes in an edge are co expressed according to and AP call
     */
    private static void CountNodeCoExpression(CyNetwork net, CyEdge cyEdge, HashMap nodeIDToAffyID, GEODataSet[] samples) {
        CyNode source = (CyNode) cyEdge.getSource();
        CyNode target = (CyNode) cyEdge.getTarget();
        DatabaseReference[] dbRefsSource = CytoscapeUtil.getDatabaseXrefsForNode(net, source);
        DatabaseReference[] dbRefsTarget = CytoscapeUtil.getDatabaseXrefsForNode(net, target);
        //try to find all affy IDs for source and target
        ArrayList affyIDsSource = new ArrayList();
        for (int i = 0; i < dbRefsSource.length; i++) {
            DatabaseReference databaseReference = dbRefsSource[i];
            ArrayList affyIDs = (ArrayList) nodeIDToAffyID.get(databaseReference);
            if (affyIDs != null) {
                affyIDsSource.addAll(affyIDs);
            }
        }
        ArrayList affyIDsTarget = new ArrayList();
        for (int i = 0; i < dbRefsTarget.length; i++) {
            DatabaseReference databaseReference = dbRefsTarget[i];
            ArrayList affyIDs = (ArrayList) nodeIDToAffyID.get(databaseReference);
            if (affyIDs != null) {
                affyIDsTarget.addAll(affyIDs);
            }
        }
        if (affyIDsSource.size() == 0 || affyIDsTarget.size() == 0) {
            //no way to link node to affy ID
            return;
        }
        int numPairPCalls = 0;
        int sourcePresentCount = 0;
        int targetPresentCount = 0;
        StringBuffer sb = null;
        for (int i = 0; i < samples.length; i++) {
            boolean sourcePresent = false;
            boolean targetPresent = false;
            GEODataSet sample = samples[i];
            //are any of the source or target affy IDs present in this sample?
            for (int j = 0; j < affyIDsSource.size(); j++) {
                String affyID = (String) affyIDsSource.get(j);
                String affyID2APCallMapping = sample.getAffyID2APCallMapping(affyID);
                if (affyID2APCallMapping != null && affyID2APCallMapping.equals("P")) {
                    sourcePresent = true;
                    sourcePresentCount++;
                    break;
                }
            }
            for (int j = 0; j < affyIDsTarget.size(); j++) {
                String affyID = (String) affyIDsTarget.get(j);
                String affyID2APCallMapping = sample.getAffyID2APCallMapping(affyID);
                if (affyID2APCallMapping != null && affyID2APCallMapping.equals("P")) {
                    targetPresent = true;
                    targetPresentCount++;
                    break;
                }
            }
            if (sourcePresent && targetPresent) {
                numPairPCalls++;
                if (sb != null) {
                    sb.append(sample.getTissueName() + "|");
                } else {
                    sb = new StringBuffer();
                    sb.append(sample.getTissueName() + "|");
                }
            }
        }
        CyAttributes cyAttributes = new CyAttributesImpl();
        cyAttributes.setAttribute(source.getIdentifier(), "ExpressionPresentCount", new Integer(sourcePresentCount));
        cyAttributes.setAttribute(target.getIdentifier(), "ExpressionPresentCount", new Integer(targetPresentCount));
        cyAttributes.setAttribute(cyEdge.getIdentifier(), "CoExpressionScore", new Integer(numPairPCalls));
        if (sb != null) {
            cyAttributes.setAttribute(cyEdge.getIdentifier(), "CoExpressionTissues", sb.toString());
        }
    }

    /**
     * Stores parameters used for the conserved hit analysis. Can also read a project file containing all
     * parameters.
     */
    private static class ConservedHitAnalysisParams {
        private ArrayList referenceParams; //One BrainParameterSet for each reference proteome search
        private ArrayList comparisonDatabaseParams; //One BrainParameterSet for each comparison proteome
        private File inparanoidSpeciesData; //The file containing a description of each species in the InparanoidDB
        private File inparanoidIDMappingFile; //The file containing the mapping of other IDs to the IDs used in the InparanoidDB
        private String inparanoidReferenceSpecies; //The reference species for the Inparanoid DB (all other proteomes compared to this one)
        private ArrayList inparanoidDataFiles; //The list of File objects pointing to all the InparanoidDB data files
        //one constant for each section of project file
        final int REF_PROTEOME = 1;
        final int REF_PROFILE = 2;
        final int COMPARISON_PROTEOME = 3;
        final int INPARANOID = 4;
        private File refProteomeFile; //The file containing the reference proteome
        private String refProteomeFormat; //The database format, e.g. fasta, if the reference proteome
        private HashMap proteomeFileName2DatabaseName; //key: (String) file name, value: (String) database name (for DatabaseReference dbname)
        private HashMap proteomeFileName2Species; //key: (String) file name, value: (String) species (e.g. ensHS)

        /**
         * Creates and reads the project file into memory
         *
         * @param projectFile The project file location
         * @throws IOException If there is an error reading the project file
         */
        public ConservedHitAnalysisParams(File projectFile) throws IOException {
            proteomeFileName2DatabaseName = new HashMap();
            proteomeFileName2Species = new HashMap();
            readProjectFile(projectFile);
        }

        public File getInparanoidSpeciesData() {
            return inparanoidSpeciesData;
        }

        public File getInparanoidIDMappingFile() {
            return inparanoidIDMappingFile;
        }

        public String getInparanoidReferenceSpecies() {
            return inparanoidReferenceSpecies;
        }

        /**
         * Get the database name (e.g. ensembl) given a database File
         *
         * @param databaseFileName The name of the proteome database file
         */
        public String getDatabaseName(File databaseFileName) {
            return ((String) proteomeFileName2DatabaseName.get(databaseFileName.toString()));
        }

        /**
         * Get the species (e.g. ensHS) given a database File
         *
         * @param databaseFileName The name of the proteome database file
         */
        public String getSpecies(File databaseFileName) {
            return ((String) proteomeFileName2Species.get(databaseFileName.toString()));
        }

        public File[] getInparanoidDataFiles() {
            File[] files = new File[inparanoidDataFiles.size()];
            inparanoidDataFiles.toArray(files);
            return files;
        }

        /**
         * Gets a list of other profiles to search against a reference proteome - stored in a BrainParameterSet object
         *
         * @return An ArrayList of BrainParameterSet objects, one for each profile set to search
         */
        public ArrayList getReferenceParams() {
            return referenceParams;
        }

        /**
         * Gets a list of comparison proteomes to search - stored in a BrainParameterSet object
         *
         * @return An ArrayList of BrainParameterSet objects, one for each comparison proteome to search
         */
        public ArrayList getComparisonDatabaseParams() {
            return comparisonDatabaseParams;
        }

        private void readProjectFile(File projectFile) throws IOException {
            BufferedReader br = new BufferedReader(new FileReader(projectFile));
            String line = null;
            String[] lineSplit;
            int section = 0;
            while ((line = br.readLine()) != null) {
                line = line.trim();
                //detect the section of the project file we're in
                if (line.equalsIgnoreCase("#ProjectFile")) {
                    continue;
                } else if (line.equalsIgnoreCase("#ReferenceProteome")) {
                    section = REF_PROTEOME;
                    if (referenceParams == null) {
                        referenceParams = new ArrayList();
                    }
                    continue;
                } else if (line.equalsIgnoreCase("#ReferenceProfiles")) {
                    section = REF_PROFILE;
                    continue;
                } else if (line.equalsIgnoreCase("#ComparisonProteomes")) {
                    section = COMPARISON_PROTEOME;
                    if (comparisonDatabaseParams == null) {
                        comparisonDatabaseParams = new ArrayList();
                    }
                    continue;
                } else if (line.equalsIgnoreCase("#Inparanoid")) {
                    section = INPARANOID;
                    continue;
                }
                if ((line.equals("")) || (line.charAt(0) == '#')) {
                    continue;
                }
                switch (section) {
                    case REF_PROTEOME:
                        //Expect e.g. D:\Gbader\Code\PDZ\data\DBs\refseq_May16_2005_human.protein.gpff\tHuman\tgenpept\trefseq(required for fasta)\n
                        lineSplit = line.split("\t");
                        if (lineSplit.length >= 3) {
                            refProteomeFile = new File(lineSplit[0]);
                            String refProteomeSpecies = lineSplit[1];
                            proteomeFileName2Species.put(refProteomeFile.toString(), refProteomeSpecies);
                            refProteomeFormat = lineSplit[2];
                            if (lineSplit.length == 4) {
                                String refProteomeDatabaseName = lineSplit[3];
                                proteomeFileName2DatabaseName.put(refProteomeFile.toString(), refProteomeDatabaseName);
                            } else if (refProteomeFormat.equalsIgnoreCase("fasta")) {
                                throw new IllegalArgumentException("#ReferenceProteome section: database format set to fasta, but no database name provided.");
                            }
                        } else {
                            throw new IllegalArgumentException("#ReferenceProteome section: expected 3 or more fields, but found " + lineSplit.length + ": " + line + ".");
                        }
                        break;
                    case REF_PROFILE:
                        //Expect: D:\Gbader\Code\PDZ\data\PDZ\BindingProfiles\Human\Final Formatted Data\TestProjectFile.txt\t11\t1.0\t20\tC\tT\n
                        lineSplit = line.split("\t");
                        if (lineSplit.length == 6) {
                            File profileFile = new File(lineSplit[0]);
                            float threshold = Float.parseFloat(lineSplit[1]);
                            float fuzz = Float.parseFloat(lineSplit[2]);
                            int top = Integer.parseInt(lineSplit[3]);
                            ProteinTerminus terminus = ProteinTerminus.parseTerminus(lineSplit[4]);
                            boolean normalized = (lineSplit[5].equalsIgnoreCase("T")) ? true : false;
                            ProteinDatabaseSearchParams dbparams = new ProteinDatabaseSearchParams(terminus);
                            dbparams.setNormalized(normalized);
                            BrainParameterSet params = new BrainParameterSet(refProteomeFile, refProteomeFormat,
                                    dbparams, profileFile, null, threshold, 100.0, top, fuzz, true, false);
                            referenceParams.add(params);
                        } else {
                            throw new IllegalArgumentException("#ReferenceProfiles section: expected 6 fields, but found " + lineSplit.length + ": " + line + ".");
                        }
                        break;
                    case COMPARISON_PROTEOME:
                        //Expect: D:\Gbader\Code\PDZ\data\PDZ\PPI Prediction\Inparanoid\downloadedInparanoid\ensAG.fa\tensAG\tfasta\tensembl(required for fasta)\n
                        lineSplit = line.split("\t");
                        if (lineSplit.length >= 3) {
                            File comparisonProteomeFile = new File(lineSplit[0]);
                            String comparisonProteomeSpecies = lineSplit[1];
                            proteomeFileName2Species.put(comparisonProteomeFile.toString(), comparisonProteomeSpecies);
                            String comparisonProteomeFormat = lineSplit[2];
                            if (lineSplit.length == 4) {
                                String comparisonProteomeDatabaseName = lineSplit[3];
                                proteomeFileName2DatabaseName.put(comparisonProteomeFile.toString(), comparisonProteomeDatabaseName);
                            } else if (comparisonProteomeFormat.equalsIgnoreCase("fasta")) {
                                throw new IllegalArgumentException("#ComparisonProteomes section: database format set to fasta, but no database name provided: " + line + ".");
                            }
                            BrainParameterSet params = new BrainParameterSet();
                            params.setDatabaseFileName(comparisonProteomeFile);
                            params.setDatabaseFormat(comparisonProteomeFormat);
                            comparisonDatabaseParams.add(params);
                        } else {
                            throw new IllegalArgumentException("#ComparisonProteomes section: expected 3 or more fields, but found " + lineSplit.length + ": " + line + ".");
                        }
                        break;
                    case INPARANOID:
                        //expect 1 line - 3 tab delimited fields
                        lineSplit = line.split("\t");
                        if (lineSplit.length == 3) {
                            inparanoidSpeciesData = new File(lineSplit[0]);
                            inparanoidReferenceSpecies = lineSplit[1];
                            inparanoidIDMappingFile = new File(lineSplit[2]);
                        } else if (lineSplit.length == 1) {
                            File inparanoidDataFile = new File(lineSplit[0]);
                            if (inparanoidDataFiles == null) {
                                inparanoidDataFiles = new ArrayList();
                            }
                            inparanoidDataFiles.add(inparanoidDataFile);
                        } else {
                            throw new IllegalArgumentException("#Inparanoid section: expected 1 or 3 fields, but found " + lineSplit.length + ": " + line + ".");
                        }
                        break;
                }
            }
        }
    }

    /**
     * Saves a single ortholog search result
     */
    private static class OrthologResult {
        public String species;
        public DatabaseReference proteinHitXref;
        public ArrayList motifHits; //List of Hit objects

        public OrthologResult(String species, DatabaseReference proteinHitXref) {
            this.species = species;
            this.proteinHitXref = proteinHitXref;
        }
    }

    /**
     * Finds orthologs between two sequence sets
     *
     * @param referenceResults       The first sequence set
     * @param referenceDatabaseName  The name of the database used to generate the first sequence set. This may be used to create
     *                               a database reference, depending on the type of sequence database i.e. GenPept contains the database reference, while
     *                               fasta format only has the database ID, no name.
     * @param comparisonResults      The second sequence set
     * @param comparisonDatabaseName The second database name
     * @param orthologyDB            The orthology database
     * @param orthologResults        Stores the results of the ortholog search - map keyed to the reference hit protein xref
     *                               Each key maps to an ArrayList of OrthologResult objects
     */
    public static void findCommonOrthologs(SequenceSearchResultSet referenceResults, String referenceDatabaseName, SequenceSearchResultSet comparisonResults, String comparisonDatabaseName, InparanoidDB orthologyDB, HashMap orthologResults) {
        Set referenceSequences = referenceResults.getSequences();
        Set comparisonSequences = comparisonResults.getSequences();
        //must check all vs all, since homology is a many:many relationship
        for (Iterator referenceIterator = referenceSequences.iterator(); referenceIterator.hasNext();) {
            Sequence referenceSequence = (Sequence) referenceIterator.next();
            Sequence referenceOriginalDBSequence = referenceResults.getOriginalSequence(referenceSequence);
            //look for a genpept sequence
            DatabaseReference referenceXref = GenPeptUtil.getAccession(referenceOriginalDBSequence);
            if (referenceXref == null) {
                //look for a fasta sequence
                referenceXref = new DatabaseReference(referenceDatabaseName, referenceOriginalDBSequence.getName());
            }
            if (referenceXref == null) {
                continue;
            }
            for (Iterator comparisonIterator = comparisonSequences.iterator(); comparisonIterator.hasNext();) {
                Sequence comparisonSequence = (Sequence) comparisonIterator.next();
                Sequence comparisonOriginalDBSequence = comparisonResults.getOriginalSequence(comparisonSequence);
                DatabaseReference comparisonXref = GenPeptUtil.getAccession(comparisonOriginalDBSequence);
                if (comparisonXref == null) {
                    //look for a fasta sequence
                    comparisonXref = new DatabaseReference(comparisonDatabaseName, comparisonOriginalDBSequence.getName());
                }
                if (comparisonXref == null) {
                    continue;
                }
                if (orthologyDB.isOrthologByAccession(referenceXref, comparisonXref)) {
                    System.out.println("Ortholog found: " + referenceXref.getDbid() + "\t" + comparisonXref.getDbid());
                    ArrayList orthologResultList = null;
                    OrthologResult result = new OrthologResult(comparisonResults.getSpecies(), comparisonXref);
                    result.motifHits = (ArrayList) comparisonResults.getHits(comparisonSequence);
                    if (orthologResults.containsKey(referenceXref)) {
                        orthologResultList = (ArrayList) orthologResults.get(referenceXref);
                        orthologResultList.add(result);
                        orthologResults.put(referenceXref, orthologResultList);
                    } else {
                        orthologResultList = new ArrayList();
                        orthologResultList.add(result);
                        orthologResults.put(referenceXref, orthologResultList);
                    }
                }
            }
        }
    }

    /**
     * Runs a conserved hit analysis. Given a list of profiles, a reference proteome, a list of comparison
     * proteomes and an ortholog database, search all proteomes with all profiles in the list and find
     * high scoring motif hits to the profiles that are orthologs
     *
     * @param projectFileName A project file which contains all parameters for the analysis
     * @param orthologDB      Optionally pass an ortholog database that may have already been initialized. If this is null,
     *                        the ortholog database will be initialized from information in the project file.
     * @throws IOException If an error occurs with any of the input files.
     */
    /*
    TODO: implement a new conserved hit analysis using a BioSQL database to store all sequences.
    Look up an ortholog group from the inparanoidDB object and then get all orthologs, examine the
    C-terminus for conservation. Determine if the reference sequence is conserved (how conserved is it?).
    Store all conserved C-termini for reporting.
    A more general motif conservation algorithm would align sequences first (external program, like Muscle),
    then see if the motif was conserved across the ortholog group.  Or you could see if there existed a
    motif (not necessarily in the same place) in one of the ortholog sequences.
    */
    public static void runConservedHitAnalysis(String projectFileName, InparanoidDB orthologDB) throws IOException {
        ConservedHitAnalysisParams analysisParams = new ConservedHitAnalysisParams(new File(projectFileName));

        if (orthologDB == null) {
            //open ortholog database if not already open
            orthologDB = new InparanoidDB();
            orthologDB.readInparanoidDataFile(analysisParams.getInparanoidDataFiles(),
                    analysisParams.getInparanoidSpeciesData(),
                    analysisParams.getInparanoidReferenceSpecies());
            orthologDB.readIDMappingFile(analysisParams.getInparanoidIDMappingFile());
        }

        ArrayList referenceParams = analysisParams.getReferenceParams();
        BrainAlgorithm alg = new BrainAlgorithm();
        for (int i = 0; i < referenceParams.size(); i++) {
            //one loop for each profile or profile set
            BrainParameterSet params = (BrainParameterSet) referenceParams.get(i);
            alg.setParams(params);
            List profileList = PeptideToProfileReader.readPeptidesAsProfiles(params.getProfileFile(), params.getFuzzFactor());
            MultiSequenceSearchResultSet referenceProteomeSearchResults = alg.runProfileSearch(profileList, null, params);
            String refDatabaseName = analysisParams.getDatabaseName(params.getDatabaseFileName());
            //set the species on all result sets
            setSpeciesOnResultSet(referenceProteomeSearchResults, analysisParams, params);

            HashMap profile2commonOrthologResultMap = new HashMap();  //key=profile, value=homolog result HashMap

            ArrayList comparisonParams = analysisParams.getComparisonDatabaseParams();
            for (int j = 0; j < comparisonParams.size(); j++) {
                //one loop for each comparison proteome
                BrainParameterSet comparisonParamSet = (BrainParameterSet) comparisonParams.get(j);
                params.setDatabaseFileName(comparisonParamSet.getDatabaseFileName());
                params.setDatabaseFormat(comparisonParamSet.getDatabaseFormat());
                MultiSequenceSearchResultSet comparisonProteomeSearchResults = alg.runProfileSearch(profileList, null, params);
                //set the species on all result sets
                setSpeciesOnResultSet(comparisonProteomeSearchResults, analysisParams, params);
                Collection referenceProteomeResults = referenceProteomeSearchResults.getAllResultSets();
                for (Iterator iterator = referenceProteomeResults.iterator(); iterator.hasNext();) {
                    //compare each reference profile search against the profile search on the comparison proteome
                    SequenceSearchResultSet referenceProteomeResultSet = (SequenceSearchResultSet) iterator.next();
                    SequenceSearchResultSet topReferenceProteomeResultSet = referenceProteomeResultSet.getTopResults(params.getNumberTopHits());
                    ProteinProfile refProfile = referenceProteomeResultSet.getProfile();
                    SequenceSearchResultSet comparisonProteomeSearchResultSet = comparisonProteomeSearchResults.getResultSet(refProfile);
                    SequenceSearchResultSet topComparisonProteomeResultSet = comparisonProteomeSearchResultSet.getTopResults(params.getNumberTopHits());
                    //create a new homolog count hashmap if required
                    HashMap homologCount;
                    if (profile2commonOrthologResultMap.containsKey(refProfile)) {
                        homologCount = (HashMap) profile2commonOrthologResultMap.get(refProfile);
                    } else {
                        homologCount = new HashMap();
                    }
                    String comparisonDatabaseName = analysisParams.getDatabaseName(params.getDatabaseFileName());
                    findCommonOrthologs(topReferenceProteomeResultSet, refDatabaseName,
                            topComparisonProteomeResultSet, comparisonDatabaseName, orthologDB, homologCount);
                    profile2commonOrthologResultMap.put(refProfile, homologCount);
                }
            }

            //print out results
            Set profiles = profile2commonOrthologResultMap.keySet();
            for (Iterator iterator = profiles.iterator(); iterator.hasNext();) {
                ProteinProfile proteinProfile = (ProteinProfile) iterator.next();
                HashMap orthologResultSet = (HashMap) profile2commonOrthologResultMap.get(proteinProfile);
                Set referenceProteinXrefs = orthologResultSet.keySet();
                for (Iterator iterator1 = referenceProteinXrefs.iterator(); iterator1.hasNext();) {
                    DatabaseReference referenceProteinXref = (DatabaseReference) iterator1.next();
                    ArrayList orthologResultList = (ArrayList) orthologResultSet.get(referenceProteinXref);
                    int count = orthologResultList.size();
                    System.out.print(proteinProfile.getName() + "\t" + referenceProteinXref.getDbid() + "\t" + count + "\t");
                    for (int j = 0; j < orthologResultList.size(); j++) {
                        OrthologResult orthologResult = (OrthologResult) orthologResultList.get(j);
                        String speciesFullName = orthologDB.getSpeciesFullName(orthologResult.species);
                        System.out.print(speciesFullName + "\t" + orthologResult.proteinHitXref.getDbid() + "\t");
                        ArrayList hits = orthologResult.motifHits;
                        for (int k = 0; k < hits.size(); k++) {
                            Hit hit = (Hit) hits.get(k);
                            System.out.print(hit + "\t");
                        }
                    }
                    System.out.println("");
                }
                System.out.println("\n");
            }
        }
    }

    /**
     * Helper method to set species on a MultiSequenceSearchResultSet
     */
    private static void setSpeciesOnResultSet(MultiSequenceSearchResultSet proteomeSearchResults, ConservedHitAnalysisParams analysisParams, BrainParameterSet params) {
        Collection referenceProteomeResults = proteomeSearchResults.getAllResultSets();
        for (Iterator iterator = referenceProteomeResults.iterator(); iterator.hasNext();) {
            SequenceSearchResultSet sequenceSearchResultSet = (SequenceSearchResultSet) iterator.next();
            sequenceSearchResultSet.setSpecies(analysisParams.getSpecies(params.getDatabaseFileName()));
        }
    }

    /* step 1: run human PDZ vs human genome (top 20)
    step 2: run worm PDZ vs worm genome (support arbitrary phage/genome pairs) (top 20)
    step 3: find homologs
    */
    public static void runConservedLinkAnalysis(String projectFileName, String comparisonProfileListFileName, String comparisonProteomeFileName, String comparisonProteomeFileFormat, InparanoidDB orthologDB) throws IOException {
        ConservedHitAnalysisParams analysisParams = new ConservedHitAnalysisParams(new File(projectFileName));

        if (orthologDB == null) {
            //open ortholog database if not already open
            orthologDB = new InparanoidDB();
            orthologDB.readInparanoidDataFile(analysisParams.getInparanoidDataFiles(),
                    analysisParams.getInparanoidSpeciesData(),
                    analysisParams.getInparanoidReferenceSpecies());
            orthologDB.readIDMappingFile(analysisParams.getInparanoidIDMappingFile());
        }

        ArrayList referenceParams = analysisParams.getReferenceParams();
        BrainAlgorithm alg = new BrainAlgorithm();
        for (int i = 0; i < referenceParams.size(); i++) {
            //one loop for each profile or profile set
            BrainParameterSet params = (BrainParameterSet) referenceParams.get(i);
            alg.setParams(params);
            MultiSequenceSearchResultSet referenceProteomeSearchResults = alg.runProfileSearch();
            String refDatabaseName = analysisParams.getDatabaseName(params.getDatabaseFileName());
            //set the species on all result sets
            setSpeciesOnResultSet(referenceProteomeSearchResults, analysisParams, params);

            //run comparison DB search with comparison profiles
            params.setProfileFileName(new File(comparisonProfileListFileName));
            params.setDatabaseFileName(new File(comparisonProteomeFileName));
            params.setDatabaseFormat(comparisonProteomeFileFormat);
            MultiSequenceSearchResultSet comparisonProteomeSearchResults = alg.runProfileSearch();
            //set the species on all result sets
            setSpeciesOnResultSet(comparisonProteomeSearchResults, analysisParams, params);

            HashMap profile2commonOrthologResultMap = new HashMap();  //key=profile, value=homolog count HashMap

            //convert both result sets to top hit lists
            Collection referenceProteomeResults = referenceProteomeSearchResults.getAllResultSets();
            MultiSequenceSearchResultSet topReferenceHitResultSet = new MultiSequenceSearchResultSet();
            for (Iterator iterator = referenceProteomeResults.iterator(); iterator.hasNext();) {
                SequenceSearchResultSet referenceProteomeResultSet = (SequenceSearchResultSet) iterator.next();
                SequenceSearchResultSet topReferenceProteomeResultSet = referenceProteomeResultSet.getTopResults(params.getNumberTopHits());
                topReferenceHitResultSet.add(topReferenceProteomeResultSet);
            }
            Collection comparisonProteomeResults = comparisonProteomeSearchResults.getAllResultSets();
            MultiSequenceSearchResultSet topComparisonHitResultSet = new MultiSequenceSearchResultSet();
            for (Iterator iterator1 = comparisonProteomeResults.iterator(); iterator1.hasNext();) {
                SequenceSearchResultSet comparisonProteomeResultSet = (SequenceSearchResultSet) iterator1.next();
                SequenceSearchResultSet topComparisonProteomeResultSet = comparisonProteomeResultSet.getTopResults(params.getNumberTopHits());
                topComparisonHitResultSet.add(topComparisonProteomeResultSet);
            }

            //match links from both sets
            referenceProteomeResults = topReferenceHitResultSet.getAllResultSets();
            for (Iterator iterator = referenceProteomeResults.iterator(); iterator.hasNext();) {
                SequenceSearchResultSet topReferenceResultSet = (SequenceSearchResultSet) iterator.next();
                ProteinProfile referenceProfile = topReferenceResultSet.getProfile();
                DatabaseReference referenceProfileProtein = referenceProfile.getProteinReference();
                comparisonProteomeResults = topComparisonHitResultSet.getAllResultSets();
                for (Iterator iterator1 = comparisonProteomeResults.iterator(); iterator1.hasNext();) {
                    SequenceSearchResultSet topComparisonResultSet = (SequenceSearchResultSet) iterator1.next();
                    ProteinProfile comparisonProfile = topComparisonResultSet.getProfile();
                    DatabaseReference comparisonProfileProtein = comparisonProfile.getProteinReference();
                    //check if reference and comparison proteins are homologs
                    if (orthologDB.isOrthologByAccession(referenceProfileProtein, comparisonProfileProtein)) {
                        System.out.println("Ortholog PDZ found: " + referenceProfile.getName() + " " + referenceProfileProtein.getDbid() + "\t" + comparisonProfile.getName() + " " + comparisonProfileProtein.getDbid());
                        //now ok to test result sets for common homologs
                        //create a new homolog count hashmap if required
                        HashMap orthologResults;
                        if (profile2commonOrthologResultMap.containsKey(referenceProfile)) {
                            orthologResults = (HashMap) profile2commonOrthologResultMap.get(referenceProfile);
                        } else {
                            orthologResults = new HashMap();
                        }
                        String comparisonDatabaseName = analysisParams.getDatabaseName(params.getDatabaseFileName());
                        findCommonOrthologs(topReferenceResultSet, refDatabaseName, topComparisonResultSet, comparisonDatabaseName, orthologDB, orthologResults);
                        profile2commonOrthologResultMap.put(referenceProfile, orthologResults);
                    }
                }
            }

            //print out results
            Set profiles = profile2commonOrthologResultMap.keySet();
            for (Iterator iterator = profiles.iterator(); iterator.hasNext();) {
                ProteinProfile proteinProfile = (ProteinProfile) iterator.next();
                HashMap orthologResultSet = (HashMap) profile2commonOrthologResultMap.get(proteinProfile);
                Set referenceProteinXrefs = orthologResultSet.keySet();
                for (Iterator iterator1 = referenceProteinXrefs.iterator(); iterator1.hasNext();) {
                    DatabaseReference referenceProteinXref = (DatabaseReference) iterator1.next();
                    ArrayList orthologResultList = (ArrayList) orthologResultSet.get(referenceProteinXref);
                    int count = orthologResultList.size();
                    System.out.print(proteinProfile.getName() + "\t" + referenceProteinXref.getDbid() + "\t" + count + "\t");
                    for (int j = 0; j < orthologResultList.size(); j++) {
                        //TODO: make output better - output species + gene name and motif from species 1.
                        OrthologResult orthologResult = (OrthologResult) orthologResultList.get(j);
                        String speciesFullName = orthologDB.getSpeciesFullName(orthologResult.species);
                        System.out.print(speciesFullName + "\t" + orthologResult.proteinHitXref.getDbid() + "\t");
                        ArrayList hits = orthologResult.motifHits;
                        for (int k = 0; k < hits.size(); k++) {
                            Hit hit = (Hit) hits.get(k);
                            System.out.print(hit + "\t");
                        }
                    }
                    System.out.println("");
                }
                System.out.println("\n");
            }
        }
    }

}
