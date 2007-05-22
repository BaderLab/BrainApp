package org.baderlab.csplugins.brainplugin;

import BiNGO.*;
import cytoscape.CyNetwork;
import cytoscape.data.annotation.Annotation;
import cytoscape.data.annotation.Ontology;
import cytoscape.task.TaskMonitor;
import org.biojava.bio.BioException;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.io.SeqIOTools;
import org.baderlab.csplugins.brainplugin.util.FileReaderUtil;

import javax.imageio.ImageIO;
import java.awt.*;
import java.io.*;
import java.math.BigDecimal;
import java.text.DateFormat;
import java.util.*;
import java.util.List;

/**
 * * Copyright (c) 2004 Memorial Sloan-Kettering Cancer Center
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
 *
 ** User: Gary Bader
 ** Date: Aug 11, 2004
 ** Time: 11:55:03 AM
 ** Description: An implementation of the BRAIN algorithm
 **/

/**
 * An implementation of the BRAIN algorithm
 */
public class BrainAlgorithm {

    private BrainParameterSet params;
    private ProteinDatabaseSearch search = null;
    private TaskMonitor taskMonitor = null;

    /**
     * The constructor.  Use this to get an instance of BRAIN to run.
     */
    public BrainAlgorithm() {
        //get current parameters
        params = BrainCurrentParameters.getInstance().getParamsCopy();
    }

    /**
     * Run a profile search from a project file or protein (list of peptides)
     *
     * @return a result set for each profile defined in the input file
     */
    public MultiSequenceSearchResultSet runProfileSearch() {
        return (runProfileSearch(null, null, null));
    }

    /**
     * Run a profile search from a project file or protein (list of peptides)
     *
     * @return a result set for each profile defined in the input file
     */
    public MultiSequenceSearchResultSet runProfileSearch(List profileList, List scoreThresholdList, BrainParameterSet internalParams) {
        MultiSequenceSearchResultSet searchResults = null;

        if (internalParams == null) {
            internalParams = params;
        }

        //get codon bias file (expects null if no bias file set)
        File codonBiasFile = params.getCodonBiasFile();

        //get unique peptides flag (expects default to be set in 'params')
        boolean uniquePeptides = params.getUniquePeptides();

        //read profile file - could be a project or single profile (list of peptides)
        if (profileList == null) {
            profileList = PeptideToProfileReader.readPeptidesAsProfiles(internalParams.getProfileFile(), internalParams.getFuzzFactor(), codonBiasFile, uniquePeptides);
        }
        //if no score list set, use the one in params for all profiles
        if (scoreThresholdList == null) {
            scoreThresholdList = new ArrayList(profileList.size());
            for (int i = 0; i < profileList.size(); i++) {
                scoreThresholdList.add(new Double(internalParams.getScoreThreshold()));
            }
        }

        //set-up a database search
        try {
            if ((internalParams.getDatabaseFileName() != null) && (internalParams.getDatabaseFormat() != null)) {
                search = new ProteinDatabaseSearch(params.getDatabaseFileName().toString(), internalParams.getDatabaseFormat());
                search.setTaskMonitor(taskMonitor);
            } else {
                System.err.println("Database filename or format not specified. Can't continue.");
                return (null);
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (BioException e) {
            e.printStackTrace();
        }
        //search the database
        try {
            searchResults = search.multiProfileSearchDB(profileList, scoreThresholdList, internalParams.getSearchParams());
        } catch (BioException e) {
            e.printStackTrace();
        }
        try {
            search.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

        return searchResults;
    }

    //run NxN profile distance calculation
    public void runAllVsAllProfileDistance(String profileListFileName) {
        final String separatorString = "\t";
        params = new BrainParameterSet();
        params.setProfileFileName(new File(profileListFileName));
        params.setFuzzFactor(0.0);
        //read profile file - could be a project or single profile (list of peptides)
        List profileList = PeptideToProfileReader.readPeptidesAsProfiles(params.getProfileFile(), params.getFuzzFactor());
        for (int i = 0; i < profileList.size(); i++) {
            ProteinProfile proteinProfile1 = (ProteinProfile) profileList.get(i);
            System.out.print(proteinProfile1.getName() + separatorString);
            for (int j = 0; j < profileList.size(); j++) {
                ProteinProfile proteinProfile2 = (ProteinProfile) profileList.get(j);
                System.out.print(ProteinProfileDistance.calculateDistributionDistance(proteinProfile1, proteinProfile2));
                if (j < profileList.size()) {
                    System.out.print(separatorString);
                }
            }
            System.out.println("");
        }
    }

    /**
     * Outputs the amino acid frequencies for all peptides contained in the profile
     * Assumes that all peptides are the same length.
     *
     * @param profileListFileName
     */
    //TODO: move this to ProteinDatabaseDistribution after that class has been refactored (and renamed to be more general)
    public void runPeptideAAByPositionHistogramAnalysis(String profileListFileName, int peptideLength) {
        //this doesn't depend on the BioJava alphabet, but is so much easier to use, and likely won't change
        final String aaList = "ACDEFGHIKLMNPQRSTVWYX-";
        long[][] histogram = new long[aaList.length()][peptideLength];

        params = new BrainParameterSet();
        params.setProfileFileName(new File(profileListFileName));
        params.setFuzzFactor(0.0);
        //read profile file - could be a project or single profile (list of peptides) - unique peptides to avoid double counting
        List profileList = PeptideToProfileReader.readPeptidesAsProfiles(params.getProfileFile(), -1, null, 0.0, null, true);
        for (int i = 0; i < profileList.size(); i++) {
            ProteinProfile proteinProfile = (ProteinProfile) profileList.get(i);
            Collection peptideList = proteinProfile.getSequenceMap();
            for (Iterator iterator = peptideList.iterator(); iterator.hasNext();) {
                Sequence sequence = (Sequence) iterator.next();
                String seqString = sequence.seqString();
                for (int k = 0; k < seqString.length(); k++) {
                    histogram[aaList.indexOf(seqString.charAt(k))][k]++;
                }
            }
        }
        System.out.print("\n");
        for (int i = 0; i < histogram.length; i++) {
            System.out.print(aaList.charAt(i) + "\t");
            for (int j = 0; j < histogram[i].length; j++) {
                System.out.print(histogram[i][j]);
                if (j < (histogram.length - 1)) {
                    System.out.print("\t");
                }
            }
            System.out.print("\n");
        }
    }

    /**
     * Outputs the amino acid frequencies for pairs of residues as position -1 and -3 of PDZ peptide ligands
     * This was created specifically for PDZ domains, but should be generalized.
     * Assumes that all peptides are the same length.
     *
     * @param profileListFileName
     */
    //TODO: change this into a more general method that finds highly correlated pairs (or higher order motifs) at given positions in a set of sequences
    public void runPeptideAAByPositionPairHistogramAnalysis(String profileListFileName, int peptideLength, int position1, int position2) {
        //this doesn't depend on the BioJava alphabet, but is so much easier to use, and likely won't change
        final String aaList = "ACDEFGHIKLMNPQRSTVWYX-";
        long[][] histogram = new long[aaList.length()][aaList.length()];

        params = new BrainParameterSet();
        params.setProfileFileName(new File(profileListFileName));
        params.setFuzzFactor(0.0);
        //read profile file - could be a project or single profile (list of peptides)
        List profileList = PeptideToProfileReader.readPeptidesAsProfiles(params.getProfileFile(), params.getFuzzFactor());
        for (int i = 0; i < profileList.size(); i++) {
            ProteinProfile proteinProfile = (ProteinProfile) profileList.get(i);
            Collection peptideList = proteinProfile.getSequenceMap();
            for (Iterator iterator = peptideList.iterator(); iterator.hasNext();) {
                Sequence sequence = (Sequence) iterator.next();
                String seqString = sequence.seqString();
                histogram[aaList.indexOf(seqString.charAt(position1))][aaList.indexOf(seqString.charAt(position2))]++;
            }
        }
        System.out.print("\n");
        for (int i = 0; i < histogram.length; i++) {
            System.out.print(aaList.charAt(i) + "\t");
            for (int j = 0; j < histogram[i].length; j++) {
                System.out.print(histogram[i][j]);
                if (j < (histogram.length - 1)) {
                    System.out.print("\t");
                }
            }
            System.out.print("\n");
        }
    }

    //cluster profiles
    //TODO: move all profile clustering code to its own profile clustering class
    public void runProfileCluster(String profileListFileName, String logoTreeTitle, File logoTreeOutput) {
        runProfileCluster(profileListFileName, -1, null, logoTreeTitle, logoTreeOutput, null, null);
    }

    /**
     * Cluster profiles and output a LogoTree
     *
     * @param profileListFileName
     * @param profileLength
     * @param terminus
     */
    public void runProfileCluster(String profileListFileName, int profileLength, ProteinTerminus terminus,
                                  String logoTreeTitle, File logoTreeOutput, File leafLabelHighlightFile, File codonBiasFile) {
        params = new BrainParameterSet();
        params.setProfileFileName(new File(profileListFileName));
        params.setFuzzFactor(0.0);
        //read profile file - could be a project or single profile (list of peptides)
        List clusterProfileList = null;
        List cutProfileList = null;
        //full profiles are always stored so we can draw a logo tree
        List originalProfileList = PeptideToProfileReader.readPeptidesAsProfiles(params.getProfileFile(), 0, null, params.getFuzzFactor(), codonBiasFile, true);
        clusterProfileList = originalProfileList;
        if (profileLength > 0) {
            //profile length specified, cluster profiles of this length (not full length ones)
            cutProfileList = new ArrayList();
            for (int i = 0; i < originalProfileList.size(); i++) {
                ProteinProfile proteinProfile = (ProteinProfile) originalProfileList.get(i);
                cutProfileList.add(proteinProfile.getTruncatedProfileCopy(profileLength, terminus));
            }
            clusterProfileList = cutProfileList;
        }
        DistanceMatrix distanceMatrix = new DistanceMatrix(clusterProfileList.size());
        distanceMatrix.calcDistances(clusterProfileList, new DistanceMetric() {
            public double calc(Object object1, Object object2) {
                ProteinProfile proteinProfile1 = (ProteinProfile) object1;
                ProteinProfile proteinProfile2 = (ProteinProfile) object2;
                String groupingByPosition = AminoAcidGrouping.getPolarChargedHydrophobeGrouping();
                return (ProteinProfileDistance.calculateAAGroupedDistributionDistance(proteinProfile1, proteinProfile2, groupingByPosition));
                /*
                //create a grouping for last 4 PDZ peptide residues
                String[] groupingByPosition = AminoAcidGrouping.getPositionSpecificPDZGrouping();
                return (ProteinProfileDistance.calculateAAGroupedByPositionDistributionDistance(proteinProfile1, proteinProfile2, groupingByPosition));
                */
            }
        });
        ArrayList al = new ArrayList(clusterProfileList.size());
        for (int i = 0; i < clusterProfileList.size(); i++) {
            ProteinProfile proteinProfile = (ProteinProfile) clusterProfileList.get(i);
            al.add(i, proteinProfile.getName());
        }
        distanceMatrix.setLabels(al);
        //cluster
        AvgLinkHierarchicalClustering cluster = new AvgLinkHierarchicalClustering(distanceMatrix);
        cluster.setOptimalLeafOrdering(true);
        cluster.run();

        /*
        try {
            cluster.writeResultsToCytoscapeFormat(new File("D:\\Gbader\\Code\\PDZ\\data\\PDZ\\PeptideProfileClustering\\GroupedByPositionDistributionDistance\\test.sif"),
                    new File("D:\\Gbader\\Code\\PDZ\\data\\PDZ\\PeptideProfileClustering\\GroupedByPositionDistributionDistance\\test.ea"), 0.25);
            writeProfileListToCytoscapeNodeAttributesFile(originalProfileList, new File("D:\\Gbader\\Code\\PDZ\\data\\PDZ\\PeptideProfileClustering\\GroupedByPositionDistributionDistance\\test.na"));

        } catch (IOException e) {
            e.printStackTrace();
        }
        */

        //summarize the clusters
        //summarizeClusters(cluster, profileList, "D:\\Gbader\\Code\\PDZ\\data\\PDZ\\PeptideClustering\\GroupedDistributionDistance\\logos");
        LogoTreeDraw ltd = new LogoTreeDraw(cluster, originalProfileList);
        ltd.setTrimNodeLogo(true, 0.2);
        ltd.setNodeLogoSubset(profileLength, terminus);
        ltd.setTitle(logoTreeTitle);
        ltd.setSymbolStyle(new PDZProteinLogoStyle());
        if (leafLabelHighlightFile != null) { // leaf label file is optional
            ArrayList list = null;
            try {
                list = FileReaderUtil.readFileAsLineList(leafLabelHighlightFile);
            } catch (IOException e) {
                e.printStackTrace();
            }
            ltd.setLeafLabelHighlightColor(list, Color.YELLOW);
        }
        ltd.setSequenceLogoStartIndex(-9);
        ltd.outputToPDF(logoTreeOutput);
    }

    /**
     * Cluster profiles and output a LogoTree for each position in a range
     */
    public void runProfileClusterPositionRange(String profileListFileName, int profileStartColumn, int profileEndColumn,
                                               int profileLength, ProteinTerminus terminus, String logoTreeTitle,
                                               File logoTreeOutput, File leafLabelHighlightFile, File codonBiasFile) {
        params = new BrainParameterSet();
        params.setProfileFileName(new File(profileListFileName));
        params.setFuzzFactor(0.0);
        //read profile file - could be a project or single profile (list of peptides)
        List clusterProfileList = null;
        List singleColumnProfileList = null;
        int maxProfileLength = 0;
        //full profiles are always stored to we can draw a logo tree
        List originalProfileList = PeptideToProfileReader.readPeptidesAsProfiles(params.getProfileFile(), 0, null, params.getFuzzFactor(), codonBiasFile, true);
        //iterate over profileLength range
        for (int i = profileStartColumn; i <= profileEndColumn; i++) {
            singleColumnProfileList = new ArrayList();
            for (int j = 0; j < originalProfileList.size(); j++) {
                ProteinProfile proteinProfile = (ProteinProfile) originalProfileList.get(j);
                singleColumnProfileList.add(proteinProfile.getProfileSubsetCopy(Integer.toString(i)));
                maxProfileLength = Math.max(maxProfileLength, proteinProfile.getNumColumns());
            }
            clusterProfileList = singleColumnProfileList;
            DistanceMatrix distanceMatrix = new DistanceMatrix(clusterProfileList.size());
            distanceMatrix.calcDistances(clusterProfileList, new DistanceMetric() {
                public double calc(Object object1, Object object2) {
                    ProteinProfile proteinProfile1 = (ProteinProfile) object1;
                    ProteinProfile proteinProfile2 = (ProteinProfile) object2;
                    //don't use a grouping for per position logo trees
                    return (ProteinProfileDistance.calculateDistributionDistance(proteinProfile1, proteinProfile2));
                }
            });
            ArrayList al = new ArrayList(clusterProfileList.size());
            for (int j = 0; j < clusterProfileList.size(); j++) {
                ProteinProfile proteinProfile = (ProteinProfile) clusterProfileList.get(j);
                al.add(j, proteinProfile.getName());
            }
            distanceMatrix.setLabels(al);
            //cluster
            AvgLinkHierarchicalClustering cluster = new AvgLinkHierarchicalClustering(distanceMatrix);
            cluster.setOptimalLeafOrdering(true);
            cluster.run();

            LogoTreeDraw ltd = new LogoTreeDraw(cluster, originalProfileList);
            ltd.setTrimNodeLogo(true, 0.2);
            ltd.setNodeLogoSubset(profileLength, terminus);
            ltd.setTitle(logoTreeTitle + "- Position " + (i - maxProfileLength + 1));
            ltd.setSymbolStyle(new PDZProteinLogoStyle());
            if (leafLabelHighlightFile != null) { // leaf label file is optional
                ArrayList list = null;
                try {
                    list = FileReaderUtil.readFileAsLineList(leafLabelHighlightFile);
                } catch (IOException e) {
                    e.printStackTrace();
                }
                ltd.setLeafLabelHighlightColor(list, Color.YELLOW);
            }
            ltd.setSequenceLogoStartIndex(-9);
            ltd.outputToPDF(new File(logoTreeOutput.toString() + "_pos" + (i - maxProfileLength + 1) + ".pdf")); //new file for each position
        }
    }

    /**
     * Cluster profiles and output a LogoTree for each position in a range
     */
    public void runProfileClusterPosition(String profileListFileName, String filter,
                                          int profileLength, ProteinTerminus terminus, String logoTreeTitle,
                                          File logoTreeOutput, File leafLabelHighlightFile, File codonBiasFile) {
        params = new BrainParameterSet();
        params.setProfileFileName(new File(profileListFileName));
        params.setFuzzFactor(0.0);
        //read profile file - could be a project or single profile (list of peptides)
        List clusterProfileList = null;
        List filteredColumnProfileList = null;
        //full profiles are always stored to we can draw a logo tree
        List originalProfileList = PeptideToProfileReader.readPeptidesAsProfiles(params.getProfileFile(), 0, null, params.getFuzzFactor(), codonBiasFile, true);
        //iterate over profileLength range
        filteredColumnProfileList = new ArrayList();
        for (int j = 0; j < originalProfileList.size(); j++) {
            ProteinProfile proteinProfile = (ProteinProfile) originalProfileList.get(j);
            filteredColumnProfileList.add(proteinProfile.getProfileSubsetCopy(filter));
        }
        clusterProfileList = filteredColumnProfileList;
        DistanceMatrix distanceMatrix = new DistanceMatrix(clusterProfileList.size());
        distanceMatrix.calcDistances(clusterProfileList, new DistanceMetric() {
            public double calc(Object object1, Object object2) {
                ProteinProfile proteinProfile1 = (ProteinProfile) object1;
                ProteinProfile proteinProfile2 = (ProteinProfile) object2;
                //                    String aaGrouping = AminoAcidGrouping.getGroupingHydrophobeBySize();
                //                   return (ProteinProfileDistance.calculateAAGroupedDistributionDistance(proteinProfile1, proteinProfile2, aaGrouping));
                return (ProteinProfileDistance.calculateDistributionDistance(proteinProfile1, proteinProfile2));
            }
        });
        ArrayList al = new ArrayList(clusterProfileList.size());
        for (int j = 0; j < clusterProfileList.size(); j++) {
            ProteinProfile proteinProfile = (ProteinProfile) clusterProfileList.get(j);
            al.add(j, proteinProfile.getName());
        }
        distanceMatrix.setLabels(al);
        //cluster
        AvgLinkHierarchicalClustering cluster = new AvgLinkHierarchicalClustering(distanceMatrix);
        cluster.setOptimalLeafOrdering(true);
        cluster.run();

        LogoTreeDraw ltd = new LogoTreeDraw(cluster, originalProfileList);
        ltd.setNodeLogoSubset(profileLength, terminus);
        ltd.setTitle(logoTreeTitle + "- Position " + filter);
        ltd.setSymbolStyle(new PDZProteinLogoStyle());
        if (leafLabelHighlightFile != null) { // leaf label file is optional
            ArrayList list = null;
            try {
                list = FileReaderUtil.readFileAsLineList(leafLabelHighlightFile);
            } catch (IOException e) {
                e.printStackTrace();
            }
            ltd.setLeafLabelHighlightColor(list, Color.YELLOW);
        }
        ltd.setSequenceLogoStartIndex(-9);
        ltd.outputToPDF(new File(logoTreeOutput.toString() + "_pos" + filter));
    }

    private void writeProfileListToCytoscapeNodeAttributesFile(List originalProfileList, File nodeAttributeFileName) throws IOException {
        BufferedWriter brNA = new BufferedWriter(new FileWriter(nodeAttributeFileName));
        brNA.write("NumberOfPeptides");
        brNA.newLine();
        for (int i = 0; i < originalProfileList.size(); i++) {
            ProteinProfile proteinProfile = (ProteinProfile) originalProfileList.get(i);
            brNA.write(proteinProfile.getName() + " = " + proteinProfile.getNumSequences());
            brNA.newLine();
        }
        brNA.close();
    }

    //run NxN profile distance calculation, output to pairs for comparison with domain sequence distance
    public void runAllVsAllProfileDistancePairOutput(String profileListFileName) {
        params = new BrainParameterSet();
        params.setProfileFileName(new File(profileListFileName));
        params.setFuzzFactor(0.0);
        //read profile file - could be a project or single profile (list of peptides)
        List profileList = PeptideToProfileReader.readPeptidesAsProfiles(params.getProfileFile(), params.getFuzzFactor());
        for (int i = 0; i < profileList.size(); i++) {
            ProteinProfile proteinProfile1 = (ProteinProfile) profileList.get(i);
            for (int j = 0; j < profileList.size(); j++) {
                ProteinProfile proteinProfile2 = (ProteinProfile) profileList.get(j);
                System.out.print(proteinProfile1.getName() + "_" + proteinProfile2.getName() + "\t");
                System.out.println(ProteinProfileDistance.calculateAAGroupedDistributionDistance(proteinProfile1, proteinProfile2, "STQN, KRH, DE, FLAMPWIVCY, G"));
            }
        }
    }

    //print out histogram of profile search score analysis
    //A graph could be output to allow people to choose a nice cutoff
    public void runHistogramScoreAnalysis(String profileListFileName, String databaseName, String dbFileFormat, ProteinTerminus term) {
        params = new BrainParameterSet();
        params.setProfileFileName(new File(profileListFileName));
        params.setFuzzFactor(1.0);
        params.setDatabaseFileName(new File(databaseName));
        params.setDatabaseFormat(dbFileFormat);
        params.setScoreThreshold(25);
        ProteinDatabaseSearchParams dbparams = new ProteinDatabaseSearchParams(term);
        dbparams.setNormalized(true);
        dbparams.setDontSaveSequences(true);
        params.setSearchParams(dbparams);
        MultiSequenceSearchResultSet searchResults = this.runProfileSearch();
        Collection results = searchResults.getAllResultSets();
        HashMap histMap = new HashMap();
        for (Iterator iterator = results.iterator(); iterator.hasNext();) {
            SequenceSearchResultSet sequenceSearchResultSet = (SequenceSearchResultSet) iterator.next();
            int histogram[] = sequenceSearchResultSet.getScoreHistogram(25);
            histMap.put(sequenceSearchResultSet.getProfile().getName(), histogram);
        }
        Set keys = histMap.keySet();
        for (Iterator iterator = keys.iterator(); iterator.hasNext();) {
            String name = (String) iterator.next();
            System.out.print(name);
            int histogram[] = (int[]) histMap.get(name);
            for (int i = 0; i < histogram.length; i++) {
                int val = histogram[i];
                System.out.print("\t" + val);
            }
            System.out.println("");
        }
    }

    //find a good score threshold for the sequence search
    //The cutoff rule is to stop when the number of hits at the current score threshold is greater than the cumulative number of hits
    //not-inclusive (does not add the last batch of hits)
    //created for phage display results of PDZ domains, where exponential score increase is observed
    public List findAutoScoreThreshold(BrainParameterSet inputParams) {
        final int maxScoreThreshold = 100;
        BrainParameterSet internalParams = new BrainParameterSet();
        internalParams.setDatabaseFileName(inputParams.getDatabaseFileName());
        internalParams.setDatabaseFormat(inputParams.getDatabaseFormat());
        internalParams.setScoreThreshold(maxScoreThreshold);
        internalParams.setFuzzFactor(inputParams.getFuzzFactor());
        ProteinDatabaseSearchParams dbparams = new ProteinDatabaseSearchParams(inputParams.getSearchParams().getTerminus());
        dbparams.setNormalized(inputParams.getSearchParams().isNormalized());
        dbparams.setDontSaveSequences(true);
        internalParams.setSearchParams(dbparams);

        List profileList = PeptideToProfileReader.readPeptidesAsProfiles(inputParams.getProfileFile(), inputParams.getFuzzFactor());

        Double[] scoreThresholdArray = new Double[profileList.size()];

        MultiSequenceSearchResultSet searchResults = runProfileSearch(profileList, null, internalParams);
        Collection results = searchResults.getAllResultSets();
        for (Iterator iterator = results.iterator(); iterator.hasNext();) {
            SequenceSearchResultSet sequenceSearchResultSet = (SequenceSearchResultSet) iterator.next();
            int histogram[] = sequenceSearchResultSet.getScoreHistogram(maxScoreThreshold);
            int cumulativeHitCount = 0;
            for (int i = 0; i < histogram.length; i++) {
                int j = histogram[i];
                if ((cumulativeHitCount > 1) && (j > cumulativeHitCount)) {
                    scoreThresholdArray[profileList.indexOf(sequenceSearchResultSet.getProfile())] = new Double(i - 1);
                    break;
                }
                cumulativeHitCount += j;
            }
            if (scoreThresholdArray[profileList.indexOf(sequenceSearchResultSet.getProfile())] == null) {
                //this means that we couldn't find any hits even when searching for maxScoreThreshold
                scoreThresholdArray[profileList.indexOf(sequenceSearchResultSet.getProfile())] = new Double(0.0);
            }
        }

        return (Arrays.asList(scoreThresholdArray));
    }

    //Answers "How many clustered protein sequences are covered by a given set of proteins?"
    //Input: set of proteins to cluster (using sequence %ID as metric) - Aligned FASTA (using e.g. Muscle)
    //Input: set of proteins of interest - how many clusters (in sequences) does this set touch? (FASTA, doesn't have to be aligned)
    //proteins from each list are matched up by name only
    //limit to columns e.g. 19,26-40,50-54,85-89,122,125-127,130-133 (index starts at 1)
    //return two columns: col1 - number of clusters; col2 - % sequence space covered by proteinsOfInterest
    public void runProteinSequenceSpaceCoverageAnalysis(File proteinsToCluster, File proteinsOfInterest, File failedProteins, String limitToColumns, int numberOfTrials) throws FileNotFoundException, BioException {
        boolean report = true; //TODO: generalize - this switches between modes - report vs. random trial
        //read proteins to cluster from an aligned fasta file and save as a List
        BufferedReader brCluster = new BufferedReader(new FileReader(proteinsToCluster));
        SequenceIterator clusterProteins = (SequenceIterator) SeqIOTools.fileToBiojava("fasta", "PROTEIN", brCluster);
        ArrayList clusterProteinSequenceStringList = new ArrayList();
        ArrayList clusterProteinNameList = new ArrayList();
        ArrayList clusterProteinSequenceList = new ArrayList();
        //create lists of sequence strings and names required for next step
        while (clusterProteins.hasNext()) {
            Sequence sequence = (Sequence) clusterProteins.nextSequence();
            clusterProteinSequenceList.add(sequence);
            if (limitToColumns != null) {
                clusterProteinSequenceStringList.add(ProteinSequenceUtil.filterSequenceByColumns(sequence, limitToColumns));
            } else {
                clusterProteinSequenceStringList.add(sequence.seqString());
            }
            clusterProteinNameList.add(sequence.getName());
        }

        //read proteins of interest from a fasta file and save as a List
        BufferedReader brInterest = new BufferedReader(new FileReader(proteinsOfInterest));
        SequenceIterator interestProteins = (SequenceIterator) SeqIOTools.fileToBiojava("fasta", "PROTEIN", brInterest);
        ArrayList interestProteinNameList = new ArrayList();
        while (interestProteins.hasNext()) {
            Sequence sequence = (Sequence) interestProteins.nextSequence();
            interestProteinNameList.add(sequence.getName());
        }

        //read failed proteins from a fasta file and save as a List
        BufferedReader brFailed = new BufferedReader(new FileReader(failedProteins));
        SequenceIterator failedProteinSequences = (SequenceIterator) SeqIOTools.fileToBiojava("fasta", "PROTEIN", brFailed);
        ArrayList failedProteinNameList = new ArrayList();
        while (failedProteinSequences.hasNext()) {
            Sequence sequence = (Sequence) failedProteinSequences.nextSequence();
            failedProteinNameList.add(sequence.getName());
        }

        //create a distance matrix
        DistanceMatrix distanceMatrix = new DistanceMatrix(clusterProteinSequenceStringList.size());
        distanceMatrix.setLabels(clusterProteinNameList);
        distanceMatrix.calcDistances(clusterProteinSequenceStringList, new AlignedProteinSequenceIdentityDistance());

        //cluster the sequences
        AvgLinkHierarchicalClustering cluster = new AvgLinkHierarchicalClustering(distanceMatrix);
        cluster.run();
        if (report) {
            System.out.print(cluster.writeResultsToGTRFormat());
            System.out.println();
            cluster.setLabelHighlightInCDTOutput(interestProteinNameList, "#FFFF00");
            cluster.setLabelHighlightInCDTOutput(failedProteinNameList, "#FF0000");
            System.out.print(cluster.toCDTString());
            System.out.println();
        }

        //find sequence space coverage of the proteins of interest sampling
        if (report) {
            System.out.println("Number of Clusters\t% Sequences Covered\tLargest Cluster Size\tSmallest Cluster Size");
        }
        double[] sequenceSpaceCoverage = new double[distanceMatrix.getMatrixDimension()];
        int[] mostBiasedClusterAssignment = null; //saves the cluster assignment that has the largest cluster size close to the size of the proteins of interest list
        int mostBiasedClusterAssignmentLargestCluster = -1;
        for (int i = 1; i <= distanceMatrix.getMatrixDimension(); i++) {
            int[] clusterAssignment = cluster.cutTree(i);
            int[] clusterSize = new int[i];
            int[] numberInterestingElementsInCluster = new int[i];
            int largestClusterSize = 0;
            int smallestClusterSize = Integer.MAX_VALUE;
            //calculate cluster stats
            for (int j = 0; j < clusterAssignment.length; j++) {
                //iterate through all clustered elements
                int clusterID = clusterAssignment[j];
                //keep track of cluster size
                clusterSize[clusterID]++;
                //keep track how many interesting genes are in each cluster
                if (interestProteinNameList.contains(distanceMatrix.getLabels().get(j))) {
                    numberInterestingElementsInCluster[clusterID]++;
                }
            }
            //find number of sequences covered
            int totalSequencesCovered = 0;
            for (int j = 0; j < numberInterestingElementsInCluster.length; j++) {
                int numElements = numberInterestingElementsInCluster[j];
                if (numElements > 0) {
                    totalSequencesCovered += clusterSize[j];
                }
            }
            //find largest and smallest cluster sizes
            int indexOfLargestCluster = -1;
            for (int j = 0; j < clusterSize.length; j++) {
                int size = clusterSize[j];
                if (largestClusterSize < size) {
                    largestClusterSize = size;
                    indexOfLargestCluster = j;
                }
                smallestClusterSize = Math.min(smallestClusterSize, size);
            }
            if ((largestClusterSize <= interestProteinNameList.size()) && (mostBiasedClusterAssignment == null)) {
                //save the cluster assignment as soon as the largest cluster is close to the size of the
                //proteins of interest list. This will be used to determine a most biased possible proteins
                //of interest selection (there could be more than one most biased selection)
                mostBiasedClusterAssignment = new int[clusterAssignment.length];
                System.arraycopy(clusterAssignment, 0, mostBiasedClusterAssignment, 0, clusterAssignment.length);
                mostBiasedClusterAssignmentLargestCluster = indexOfLargestCluster;
            }
            //save result
            sequenceSpaceCoverage[i - 1] = ((double) totalSequencesCovered / (double) distanceMatrix.getMatrixDimension());
            if (report) {
                System.out.println(i + "\t" + ((double) totalSequencesCovered / (double) distanceMatrix.getMatrixDimension()) + "\t" + largestClusterSize + "\t" + smallestClusterSize);
            }
        }
        if (report) {
            System.out.println();
        }

        //find the worst case sampling scenario (most biased sampling) - when closest sequences to each other are picked
        //create a name list for the most biased sampling
        ArrayList mostBiasedNameList = new ArrayList();
        for (int i = 0; i < mostBiasedClusterAssignment.length; i++) {
            int clusterIndex = mostBiasedClusterAssignment[i];
            if (clusterIndex == mostBiasedClusterAssignmentLargestCluster) {
                mostBiasedNameList.add(clusterProteinNameList.get(i));
            }
        }
        //calculate sequence space coverage over different cluster tree cuts
        //TODO: move this and repeated code in other places into a general method
        for (int k = 1; k <= distanceMatrix.getMatrixDimension(); k++) {
            int[] clusterAssignment = cluster.cutTree(k);
            int[] clusterSize = new int[k];
            int[] numberInterestingElementsInCluster = new int[k];
            //calculate cluster stats
            for (int j = 0; j < clusterAssignment.length; j++) {
                //iterate through all clustered elements
                int clusterID = clusterAssignment[j];
                //keep track of cluster size
                clusterSize[clusterID]++;
                //keep track how many interesting genes are in each cluster
                if (mostBiasedNameList.contains(distanceMatrix.getLabels().get(j))) {
                    numberInterestingElementsInCluster[clusterID]++;
                }
            }
            //find number of sequences covered
            int totalSequencesCovered = 0;
            for (int j = 0; j < numberInterestingElementsInCluster.length; j++) {
                int numElements = numberInterestingElementsInCluster[j];
                if (numElements > 0) {
                    totalSequencesCovered += clusterSize[j];
                }
            }
            //save result
            sequenceSpaceCoverage[k - 1] = ((double) totalSequencesCovered / (double) distanceMatrix.getMatrixDimension());
        }
        if (report) {
            for (int v = 0; v < sequenceSpaceCoverage.length; v++) {
                System.out.println((v + 1) + "\t" + sequenceSpaceCoverage[v]);
            }
        }

        //run random model of sequence space coverage (what is the coverage of a random list?)
        double[] sequenceSpaceCoverageSum = null;
        //calculate numberOfTrials number of random proteinsOfInterest lists
        for (int i = 0; i < numberOfTrials; i++) {
            //create a random interestProteinsList (sampled from the clusterProteinsList)
            ArrayList randomSequencelist = createRandomSamplingWithoutReplacement(clusterProteinSequenceList, interestProteinNameList.size());
            ArrayList randomSequenceNameList = new ArrayList();
            for (int j = 0; j < randomSequencelist.size(); j++) {
                Sequence sequence = (Sequence) randomSequencelist.get(j);
                randomSequenceNameList.add(sequence.getName());
            }
            //calculate sequence space coverage over different cluster tree cuts
            for (int k = 1; k <= distanceMatrix.getMatrixDimension(); k++) {
                int[] clusterAssignment = cluster.cutTree(k);
                int[] clusterSize = new int[k];
                int[] numberInterestingElementsInCluster = new int[k];
                //calculate cluster stats
                for (int j = 0; j < clusterAssignment.length; j++) {
                    //iterate through all clustered elements
                    int clusterID = clusterAssignment[j];
                    //keep track of cluster size
                    clusterSize[clusterID]++;
                    //keep track how many interesting genes are in each cluster
                    if (randomSequenceNameList.contains(distanceMatrix.getLabels().get(j))) {
                        numberInterestingElementsInCluster[clusterID]++;
                    }
                }
                //find number of sequences covered
                int totalSequencesCovered = 0;
                for (int j = 0; j < numberInterestingElementsInCluster.length; j++) {
                    int numElements = numberInterestingElementsInCluster[j];
                    if (numElements > 0) {
                        totalSequencesCovered += clusterSize[j];
                    }
                }
                //save result
                sequenceSpaceCoverage[k - 1] = ((double) totalSequencesCovered / (double) distanceMatrix.getMatrixDimension());
            }
            if (sequenceSpaceCoverageSum == null) {
                sequenceSpaceCoverageSum = new double[sequenceSpaceCoverage.length];
            }
            if (!report) {
                for (int v = 0; v < sequenceSpaceCoverage.length; v++) {
                    System.out.println((v + 1) + "\t" + sequenceSpaceCoverage[v]);
                }
            }
            //keep track of sum of sequence space coverage values for later average determination
            for (int j = 0; j < sequenceSpaceCoverage.length; j++) {
                sequenceSpaceCoverageSum[j] += sequenceSpaceCoverage[j];
            }
        }
        if (report) {
            double[] sequenceSpaceCoverageAverage = new double[sequenceSpaceCoverageSum.length];
            for (int i = 0; i < sequenceSpaceCoverageSum.length; i++) {
                sequenceSpaceCoverageAverage[i] = sequenceSpaceCoverageSum[i] / numberOfTrials;
                System.out.println((i + 1) + "\t" + sequenceSpaceCoverageAverage[i]);
            }
        }
    }

    private ArrayList createRandomSamplingWithoutReplacement(ArrayList itemsToSample, int numberOfSamples) {
        if (numberOfSamples > itemsToSample.size()) {
            throw new RuntimeException("Number of samples must be smaller than the size of the list given to sample.");
        }
        if (numberOfSamples == itemsToSample.size()) {
            return itemsToSample;
        }
        ArrayList randomSampling = new ArrayList();
        ArrayList sampleCopy = (ArrayList) itemsToSample.clone();
        for (int i = 0; i < numberOfSamples; i++) {
            int randomIndex = (int) Math.round((sampleCopy.size() - 1) * Math.random());
            randomSampling.add(sampleCopy.remove(randomIndex));
        }
        return randomSampling;
    }

    public void clusterProteinSequences(File proteinsToCluster) throws FileNotFoundException, BioException {
        //read proteins to cluster from an aligned fasta file and save as a List
        BufferedReader brCluster = new BufferedReader(new FileReader(proteinsToCluster));
        SequenceIterator clusterProteins = (SequenceIterator) SeqIOTools.fileToBiojava("fasta", "PROTEIN", brCluster);
        ArrayList clusterProteinSequenceStringList = new ArrayList();
        ArrayList clusterProteinNameList = new ArrayList();
        ArrayList clusterProteinSequenceList = new ArrayList();
        //create lists of sequence strings and names required for next step
        while (clusterProteins.hasNext()) {
            Sequence sequence = (Sequence) clusterProteins.nextSequence();
            clusterProteinSequenceList.add(sequence);
            clusterProteinSequenceStringList.add(sequence.seqString());
            clusterProteinNameList.add(sequence.getName());
        }

        //create a distance matrix
        DistanceMatrix distanceMatrix = new DistanceMatrix(clusterProteinSequenceStringList.size());
        distanceMatrix.setLabels(clusterProteinNameList);
        distanceMatrix.calcDistances(clusterProteinSequenceStringList, new AlignedProteinSequenceIdentityDistance());

        //cluster the sequences
        AvgLinkHierarchicalClustering cluster = new AvgLinkHierarchicalClustering(distanceMatrix);
        cluster.run();
        System.out.print(cluster.writeResultsToGTRFormat());
        System.out.println();
        System.out.print(cluster.toCDTString());
        System.out.println();
    }

    //run a full search and save results to Cytoscape
    public void runFullRegexSearchToCytoscape(String proteinName, String queryDatabaseFileName,
                                              ProteinTerminus searchTerminus, int motifLength, String[] regexList, CyNetwork net) {
        ProteinDatabaseSearch search = null;
        try {
            search = new ProteinDatabaseSearch(queryDatabaseFileName, "FASTA") {
                //TODO: move this into a more general place
                public String getIdentifier(Sequence sequence) {
                    return CytoscapeUtil.getIPI_SPIdentifier(sequence);
                }
            };
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (BioException e) {
            e.printStackTrace();
        }
        //run the search
        SequenceSearchResultSet results = null;
        for (int i = 0; i < regexList.length; i++) {
            try {
                ProteinDatabaseSearchParams params = new ProteinDatabaseSearchParams(searchTerminus, motifLength);
                results = search.regexSearchDB(regexList[i], params);
            } catch (BioException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }
            DatabaseReference protein = new DatabaseReference("SwissProt", proteinName);
            //CytoscapeUtil.addSequenceSearchResultSetToCytoscape(net, protein, results, params);
            System.out.println("Found " + results.getNumberSequencesHit() + " matching sequences out of " +
                    results.getNumberOfSequencesSearched() + " for regular expression " + regexList[i] + ".");
            try {
                search.reset();
            } catch (IOException e) {
                e.printStackTrace();
            } catch (BioException e) {
                e.printStackTrace();
            }
        }
        //System.out.println(results);
        SequenceSearchResultSet uniqueResults = results.getUniqueResultsBySequence();
        System.out.println("Found " + uniqueResults.getNumberSequencesHit() + " unique results.");
        System.out.println(uniqueResults);
        System.out.println("\n\n\n");
        PeptideChipDesign chip = null;
        try {
            chip = new PeptideChipDesign("D:\\Gbader\\Code\\PDZ\\data\\PDZ\\Protein chip design\\ecoli_preferredCodons.txt", uniqueResults);
        } catch (IOException e) {
            e.printStackTrace();
        }
        chip.setDescription("Worm chip");
        System.out.println(chip);

    }

    /**
     * Get the parameter set used for this instance of BrainAlgorithm
     *
     * @return The parameter set used
     */
    public BrainParameterSet getParams() {
        return params;
    }

    /**
     * Set the parameter set used for this instance of BrainAlgorithm
     */
    public void setParams(BrainParameterSet params) {
        this.params = params;
    }

    /**
     * If called, will cancel the task at the next convenient moment.
     */
    public void cancel() {
        search.setCancelled(true);
    }

    /**
     * Sets the task monitor for this task (if it is being run inside of a Task)
     */
    public void setTaskMonitor(TaskMonitor taskMonitor) {
        this.taskMonitor = taskMonitor;
    }

    //write sequence logos for all profiles
    public void writeSequenceLogos(String profileListFileName, String outputDirectory, String codonBiasFileName, double fuzzFactor) {
        File codonBiasFile = null;
        if (codonBiasFileName != null) {
            codonBiasFile = new File(codonBiasFileName);
        }
        List profileList = PeptideToProfileReader.readPeptidesAsProfiles(new File(profileListFileName), 0, null,
                fuzzFactor, codonBiasFile, true);
        String outFileName = null;
        for (int i = 0; i < profileList.size(); i++) {
            ProteinProfile proteinProfile = (ProteinProfile) profileList.get(i);
            outFileName = new String(outputDirectory + File.separator + proteinProfile.getName() + ".png");
            ProteinSequenceLogo logo = new ProteinSequenceLogo(proteinProfile, 240);
            try {
                logo.sequenceLogoSetStartIndex(-9);
                ImageIO.write(logo.drawSequenceLogo(), "png", new File(outFileName));
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }

    //GO term overrepresentation analysis on predicted PPIs on a per profile cluster basis
    public void findOverrepresentedGOTermsForAllPredictedPPIsForAllPossibleProfileClusters(AvgLinkHierarchicalClustering cluster, MultiSequenceSearchResultSet inputSearchResults) {
        //TODO: there is a huge amount of unnecessary code in here that is just copied from BiNGO. All of BiNGO needs to be refactored to make it modular
        //so that you can use it as a library, then this entire class has to be refactored to use the modular bingo classes
        //initialization
        //for each profile, store a list of Entrez Gene IDs of predicted interacting proteins
        HashMap profileLabelStringToEntrezGeneIDStringListMap = new HashMap();
        HashMap entrezGeneIDToGeneName = new HashMap();
        HashMap entrezGeneIDToMotifAndScoreString = new HashMap();
        Collection searchResults = inputSearchResults.getAllResultSets();
        //there is one sequence search result per profile.  Go through each result and extract the Gene IDs
        for (Iterator iterator = searchResults.iterator(); iterator.hasNext();) {
            SequenceSearchResultSet sequenceSearchResultSet = (SequenceSearchResultSet) iterator.next();
            //convert sequence search results to Entrez Gene IDs
            Set sequences = sequenceSearchResultSet.getSequences();
            ArrayList entrezGeneIDs = new ArrayList();
            for (Iterator iterator1 = sequences.iterator(); iterator1.hasNext();) {
                Sequence sequence = (Sequence) iterator1.next();
                Sequence originalDBSequence = sequenceSearchResultSet.getOriginalSequence(sequence);
                DatabaseReference dbRef = GenPeptUtil.getEntrezGeneID(originalDBSequence);
                if (dbRef != null) {
                    String entrezGeneID = dbRef.getDbid();
                    entrezGeneIDs.add(entrezGeneID);
                    if (!entrezGeneIDToGeneName.containsKey(entrezGeneID)) {
                        entrezGeneIDToGeneName.put(entrezGeneID, GenPeptUtil.getGeneName(originalDBSequence));
                    }
                    if (!entrezGeneIDToMotifAndScoreString.containsKey(entrezGeneID)) {
                        double bestScore = Double.MAX_VALUE;
                        String bestMotif = null;
                        List hits = sequenceSearchResultSet.getHits(sequence);
                        for (int i = 0; i < hits.size(); i++) {
                            Hit hit = (Hit) hits.get(i);
                            //keep track of best score
                            if (hit.getScore().doubleValue() < bestScore) {
                                bestScore = hit.getScore().doubleValue();
                                bestMotif = hit.getMatchString();
                            }
                        }
                        bestScore = truncateDouble(bestScore, 3);
                        entrezGeneIDToMotifAndScoreString.put(entrezGeneID, bestMotif + ":" + bestScore);
                    }
                } else {
                    System.out.println("Couldn't find a database reference for " + sequence.getName());
                }
            }
            //save in the map and move to the next profile search result set
            profileLabelStringToEntrezGeneIDStringListMap.put(sequenceSearchResultSet.getProfile().getName(), entrezGeneIDs);
        }

        //iterate through all tree cuts and do GO ORA (over-representation analysis) on each cluster
        for (int i = 1; i <= cluster.getNelements(); i++) {
            int[] clusterAssignment = cluster.cutTree(i);
            //collect profile labels for all leaves in a cluster
            ArrayList[] clusterIDToProfileLabelStringListArray = new ArrayList[i];
            //initialize the array
            for (int j = 0; j < i; j++) {
                clusterIDToProfileLabelStringListArray[j] = new ArrayList();
            }
            for (int j = 0; j < clusterAssignment.length; j++) {
                //iterate through all clustered elements
                int clusterID = clusterAssignment[j];
                //save the list of profiles that are in each cluster (by their unique string label)
                clusterIDToProfileLabelStringListArray[clusterID].add(cluster.getLabel(j));
            }
            //for each cluster, collect all entrez gene IDs into one large unique list for ORA
            Vector[] clusterIDToEntrezGeneIDStringVectorArray = new Vector[i];
            //initialize the array
            for (int j = 0; j < i; j++) {
                clusterIDToEntrezGeneIDStringVectorArray[j] = new Vector();
            }
            //collect the vectors of entrez gene IDs for each cluster into one vector of unique IDs per cluster
            for (int j = 0; j < i; j++) {
                //get all profile labels in the cluster
                ArrayList profileLabelList = clusterIDToProfileLabelStringListArray[j];
                for (int k = 0; k < profileLabelList.size(); k++) {
                    String profileLabel = (String) profileLabelList.get(k);
                    //get all entrez Gene IDs for each profile label
                    ArrayList entrezGeneIDs = (ArrayList) profileLabelStringToEntrezGeneIDStringListMap.get(profileLabel);
                    //add all entrez Gene IDs to the cluster, but maintain a unique vector
                    for (int l = 0; l < entrezGeneIDs.size(); l++) {
                        String entrezGeneID = (String) entrezGeneIDs.get(l);
                        if (!clusterIDToEntrezGeneIDStringVectorArray[j].contains(entrezGeneID)) {
                            clusterIDToEntrezGeneIDStringVectorArray[j].add(entrezGeneID);
                        }
                    }
                }
            }
            String clusterCutDirName = "D:\\Gbader\\Code\\PDZ\\data\\PDZ\\PPI Prediction\\ORA\\ClusterBiNGO\\ClusterCutInto" + i + "\\";
            File newDir = new File(clusterCutDirName);
            newDir.mkdirs();

            //run GO ORA on each of these entrez gene ID lists
            for (int j = 0; j < i; j++) {
                Vector entrezGeneIDList = clusterIDToEntrezGeneIDStringVectorArray[j];
                //call BiNGO to perform the ORA
                String bingoDir = "D:\\Gbader\\Code\\PDZ\\data\\DBs\\BiNGO-data-v1_1_nov2005";
                SettingsPanel settingsPanel = new SettingsPanel(bingoDir);
                BiNGOOntologyFlatFileReader readerOntology = null;
                String ontologyFile = "D:\\Gbader\\Code\\PDZ\\data\\DBs\\BiNGO-data-v1_1_nov2005\\GO_Biological_Process";
                try {
                    readerOntology = new BiNGOOntologyFlatFileReader(new File(ontologyFile));
                } catch (Exception e) {
                    e.printStackTrace();
                }
                Ontology ontology = readerOntology.getOntology();
                HashMap synonymHash = readerOntology.getSynonymHash();
                BiNGOAnnotationDefaultReader readerAnnotation = null;
                String annotationFile = "D:\\Gbader\\Code\\PDZ\\data\\DBs\\BiNGO-data-v1_1_nov2005\\H_sapiens_default";
                try {
                    String idString = "Entrez GeneID";
                    readerAnnotation = new BiNGOAnnotationDefaultReader(new File(annotationFile), synonymHash, settingsPanel, idString, "Homo Sapiens", "GO Biological Process", "GO");
                    // mdharsee 20070212 - was trying new BiNGO code (v2.0) requiring a different call to BiNGOAnnotationDefaultReader and use of BingoParameters instead of SettingsPanel
                    //  - but it seems the MultipleTestingCorrection code used below is not supported by the latest version
                    //  - reverting back to the previous working version of BiNGO referenced in original code from GB
                    //readerAnnotation = new BiNGOAnnotationDefaultReader(new File(annotationFile), synonymHash, bingoParams, idString, "Homo Sapiens");
                } catch (Exception e) {
                    e.printStackTrace();
                }
                Annotation annotation = readerAnnotation.getAnnotation();
                Vector selectedNodes = entrezGeneIDList;
                String[] nodes = annotation.getNames();
                // vector for storing the canonical names
                Vector allNodes = new Vector();
                // iterate over every node view to get the canonical names.
                for (int n = 0; n < nodes.length; n++) {
                    if (nodes[n] != null && nodes[n].length() != 0)
                        allNodes.add(nodes[n].toUpperCase());
                }
                //calculate over-representation using hypergeometric
                HypergeometricTestCalculate test = new HypergeometricTestCalculate(selectedNodes, allNodes, annotation, ontology);
                test.calculate();
                HashMap testMap = test.getTestMap();
                HashMap mapSmallX = test.getMapSmallX();
                HashMap mapSmallN = test.getMapSmallN();
                int bigX = test.getBigX();
                int bigN = test.getBigN();
                String alphaString = "0.05";
                String correctionString = "Benjamini & Hochberg False Discovery Rate (FDR) correction";
                MultipleTestingCorrection mtc = new MultipleTestingCorrection(alphaString, testMap, correctionString);
                mtc.calculate();
                HashMap correctionMap = mtc.getCorrectionMap();
                try {
                    String testString = "Hypergeometric";
                    String overUnderString = "Over-representation";
                    //create a new directory for each cluster set
                    String dirName = clusterCutDirName;
                    String fileName = "Node" + Integer.toString(j) + ".txt";
                    String clusterVsString = "Vs. genome";
                    String catString = "Over-represented categories after correction";
                    Vector selectedCanonicalNameVector = selectedNodes;
                    HashMap annotatedGenes = new HashMap();
                    String dateString = DateFormat.getDateInstance().format(new Date());
                    String timeString = DateFormat.getTimeInstance().format(new Date());
                    String NONE = "---";
                    /** constant string for the checking of numbers of categories, before correction.*/
                    String CATEGORY_BEFORE_CORRECTION = "Overrepresented categories before correction";
                    /** constant string for the checking of numbers of categories, after correction.*/
                    String CATEGORY_CORRECTION = "Overrepresented categories after correction";

                    File results = new File(dirName, fileName);
                    BufferedWriter output = new BufferedWriter(new FileWriter(results));
                    System.out.println("BiNGO results file: " + results.getPath());
                    output.write("File created with BiNGO (c) on " + dateString + " at " + timeString);
                    output.newLine();
                    output.write(ontology.toString());
                    output.newLine();
                    output.write("Selected ontology file: " + ontologyFile);
                    output.newLine();
                    output.write("Selected annotation file: " + annotationFile);
                    output.newLine();
                    output.write(overUnderString);
                    output.newLine();
                    output.write("Selected statistical test: " + testString);
                    output.newLine();
                    output.write("Selected correction: " + correctionString);
                    output.newLine();
                    output.write("Selected significance level: " + alphaString);
                    output.newLine();
                    output.write("Testing option: " + clusterVsString);
                    output.newLine();
                    output.write("The selected profiles: ");
                    ArrayList arrayList = clusterIDToProfileLabelStringListArray[j];
                    for (int l = 0; l < arrayList.size(); l++) {
                        String s = (String) arrayList.get(l);
                        output.write(s);
                        if (l < (arrayList.size() - 1)) {
                            output.write(", ");
                        }
                    }
                    output.newLine();
                    output.write("The selected predicted binders: ");
                    for (int n = 0; n < selectedCanonicalNameVector.size(); n++) {
                        int[] nodeClassifications = annotation.getClassifications(selectedCanonicalNameVector.get(n).toString());
                        for (int k = 0; k < nodeClassifications.length; k++) {
                            String cat = new Integer(nodeClassifications[k]).toString();
                            if (!annotatedGenes.containsKey(cat)) {
                                HashSet catset = new HashSet();
                                annotatedGenes.put(cat, catset);
                            }
                            ((HashSet) annotatedGenes.get(cat)).add(selectedCanonicalNameVector.get(n).toString());
                        }
                        output.write(entrezGeneIDToGeneName.get(selectedCanonicalNameVector.get(n).toString()) + " (" + entrezGeneIDToMotifAndScoreString.get(selectedCanonicalNameVector.get(n).toString()) + ")");
                        if (n < (selectedCanonicalNameVector.size() - 1)) {
                            output.write(", ");
                        }
                    }
                    output.newLine();
                    output.newLine();
                    output.write("Number of genes selected: " + bigX);
                    output.newLine();
                    output.write("Total number of genes in annotation: " + bigN);
                    output.newLine();
                    output.newLine();
                    if (testString.equals(NONE)) {
                        output.write("GO-ID" + "\t" + "# selected" + "\t" + "# total" + "\t" + "Description" + "\t" + "Genes in test set");
                        output.newLine();
                    } else if (correctionString.equals(NONE)) {
                        output.write("GO-ID" + "\t" + "p-value" + "\t" + "# selected" + "\t" + "# total" + "\t" + "Description" + "\t" + "Genes in test set");
                        output.newLine();
                    } else {
                        output.write("GO-ID" + "\t" + "p-value" + "\t" + "corr p-value" + "# selected" + "\t" + "# total" + "\t" + "Description" + "\t" + "Genes in test set");
                        output.newLine();
                    }

                    //order GO labels by increasing corrected p-value or increasing smallX

                    HashSet keySet;
                    if (!testString.equals(NONE)) {
                        keySet = new HashSet(testMap.keySet());
                    } else {
                        keySet = new HashSet(mapSmallX.keySet());
                    }
                    Iterator it = keySet.iterator();
                    String[] keyLabels = new String[keySet.size()];
                    for (int n = 0; it.hasNext(); n++) {
                        keyLabels[n] = it.next().toString();
                    }
                    String[] ordenedKeySet;
                    if (!testString.equals(NONE)) {
                        ordenedKeySet = orderKeysByPvalues(keyLabels, testMap);
                    } else {
                        ordenedKeySet = orderKeysBySmallX(keyLabels, mapSmallX);
                    }
                    boolean ok = true;

                    for (int n = 0; (n < ordenedKeySet.length) && (ok == true); n++) {

                        String termID = ordenedKeySet[n];
                        String pvalue = "";
                        String correctedPvalue = "";
                        String smallX;
                        String smallN;
                        String description;
                        // pvalue
                        if (!testString.equals(NONE)) {
                            try {
                                pvalue = SignificantFigures.sci_format(testMap.get(new Integer(termID)).toString(), 5);
                            } catch (Exception e) {
                                pvalue = "N/A";
                            }
                        } else {
                            pvalue = "N/A";
                        }
                        // corrected pvalue
                        if (!correctionString.equals(NONE)) {
                            try {
                                correctedPvalue = SignificantFigures.sci_format(correctionMap.get(termID).toString(), 5);
                            } catch (Exception e) {
                                correctedPvalue = "N/A";
                            }
                        } else {
                            correctedPvalue = "N/A";
                        }
                        // x
                        try {
                            smallX = mapSmallX.get(new Integer(termID)).toString();
                        } catch (Exception e) {
                            smallX = "N/A";
                        }
                        // n
                        try {
                            smallN = mapSmallN.get(new Integer(termID)).toString();
                        } catch (Exception e) {
                            smallN = "N/A";
                        }
                        // name
                        try {
                            description = ontology.getTerm(Integer.parseInt(termID)).getName();
                        } catch (Exception e) {
                            description = "?";
                        }

                        if (testString.equals(NONE)) {
                            output.write(termID + "\t" + smallX + "\t" + smallN + "\t" + description + "\t");
                            if (annotatedGenes.containsKey(termID)) {
                                Iterator k = ((HashSet) annotatedGenes.get(termID)).iterator();
                                while (k.hasNext()) {
                                    output.write(k.next().toString());
                                    if (k.hasNext()) {
                                        output.write('|');
                                    }
                                }
                            }
                            output.write("\n");
                        } else if (correctionString.equals(NONE)) {
                            if (catString.equals(CATEGORY_BEFORE_CORRECTION)) {
                                if ((new BigDecimal(testMap.get(new Integer(ordenedKeySet[n])).toString())).compareTo(new BigDecimal(alphaString)) < 0) {
                                    output.write(termID + "\t" + pvalue + "\t" + smallX + "\t" + smallN + "\t" + description + "\t");
                                    if (annotatedGenes.containsKey(termID)) {
                                        Iterator k = ((HashSet) annotatedGenes.get(termID)).iterator();
                                        while (k.hasNext()) {
                                            output.write(k.next().toString());
                                            if (k.hasNext()) {
                                                output.write('|');
                                            }
                                        }
                                    }
                                    output.newLine();
                                } else {
                                    ok = false;
                                }
                            } else {
                                output.write(termID + "\t" + pvalue + "\t" + smallX + "\t" + smallN + "\t" + description + "\t");
                                if (annotatedGenes.containsKey(termID)) {
                                    Iterator k = ((HashSet) annotatedGenes.get(termID)).iterator();
                                    while (k.hasNext()) {
                                        output.write(k.next().toString());
                                        if (k.hasNext()) {
                                            output.write('|');
                                        }
                                    }
                                }
                                output.write("\n");
                            }
                        } else {
                            if (catString.equals(CATEGORY_CORRECTION)) {
                                if ((new BigDecimal(correctionMap.get(ordenedKeySet[n]).toString())).compareTo(new BigDecimal(alphaString)) < 0) {
                                    output.write(termID + "\t" + pvalue + "\t" + correctedPvalue + "\t" + smallX + "\t" + smallN + "\t" + description + "\t");
                                    if (annotatedGenes.containsKey(termID)) {
                                        Iterator k = ((HashSet) annotatedGenes.get(termID)).iterator();
                                        while (k.hasNext()) {
                                            output.write((String) entrezGeneIDToGeneName.get(k.next().toString()));
                                            if (k.hasNext()) {
                                                output.write('|');
                                            }
                                        }
                                    }
                                    output.newLine();
                                } else {
                                    ok = false;
                                }
                            } else if (catString.equals(CATEGORY_BEFORE_CORRECTION)) {
                                if ((new BigDecimal(testMap.get(new Integer(ordenedKeySet[n])).toString())).compareTo(new BigDecimal(alphaString)) < 0) {
                                    output.write(termID + "\t" + pvalue + "\t" + correctedPvalue + "\t" + smallX + "\t" + smallN + "\t" + description + "\t");
                                    if (annotatedGenes.containsKey(termID)) {
                                        Iterator k = ((HashSet) annotatedGenes.get(termID)).iterator();
                                        while (k.hasNext()) {
                                            output.write(k.next().toString());
                                            if (k.hasNext()) {
                                                output.write('|');
                                            }
                                        }
                                    }
                                    output.newLine();
                                } else {
                                    ok = false;
                                }
                            } else {
                                output.write(termID + "\t" + pvalue + "\t" + correctedPvalue + "\t" + smallX + "\t" + smallN + "\t" + description + "\t");
                                if (annotatedGenes.containsKey(termID)) {
                                    Iterator k = ((HashSet) annotatedGenes.get(termID)).iterator();
                                    while (k.hasNext()) {
                                        output.write((String) entrezGeneIDToGeneName.get(k.next().toString()));
                                        if (k.hasNext()) {
                                            output.write('|');
                                        }
                                    }
                                }
                                output.newLine();
                            }
                        }
                    }

                    output.close();
                } catch (Exception e) {
                    System.out.println("Error: " + e);

                }
            }
        }
    }

    /**
     * Truncate a double-precision floating point number to a specific number of significant digits
     *
     * @return The new truncated double number (make sure to catch it)
     */
    private static double truncateDouble(double inputDouble, int numberSignificantDigits) {
        double returnDouble = inputDouble * Math.pow(10, numberSignificantDigits);
        returnDouble = Math.rint(returnDouble);
        return (returnDouble / Math.pow(10, numberSignificantDigits));
    }

    public String[] orderKeysByPvalues(String[] labels, HashMap testMap) {

        for (int i = 1; i < labels.length; i++) {
            int j = i;
            // get the first unsorted value ...
            String insert_label = labels[i];
            BigDecimal val = new BigDecimal(testMap.get(new Integer(labels[i])).toString());
            // ... and insert it among the sorted
            while ((j > 0) && (val.compareTo(new BigDecimal(testMap.get(new Integer(labels[j - 1])).toString())) < 0)) {
                labels[j] = labels[j - 1];
                j--;
            }
            // reinsert value
            labels[j] = insert_label;
        }
        return labels;
    }

    public String[] orderKeysBySmallX(String[] labels, HashMap mapSmallX) {

        for (int i = 1; i < labels.length; i++) {
            int j = i;
            // get the first unsorted value ...
            String insert_label = labels[i];
            BigDecimal val = new BigDecimal(mapSmallX.get(new Integer(labels[i])).toString());
            // ... and insert it among the sorted
            while ((j > 0) && (val.compareTo(new BigDecimal(mapSmallX.get(new Integer(labels[j - 1])).toString())) > 0)) {
                labels[j] = labels[j - 1];
                j--;
            }
            // reinsert value
            labels[j] = insert_label;
        }
        return labels;
    }

    public double calculateSpecificityScoreCorrectionFactor(String codonBiasFileName) {
        File codonBiasFile = null;
        if (codonBiasFileName != null) {
            codonBiasFile = new File(codonBiasFileName); //only used to correct the specificity score by the -9 position
        }

        final String aaList = "ACDEFGHIKLMNPQRSTVWYX";
        double[] codonBias = new double[aaList.length()];

        //read the codon bias file
        try {
            BufferedReader br = new BufferedReader(new FileReader(codonBiasFile));
            String fileLine = null;
            while ((fileLine = br.readLine()) != null) {
                String[] aaBiasString = fileLine.split("\\t");
                if (aaBiasString.length == 2) {
                    String residue = aaBiasString[0];
                    codonBias[aaList.indexOf(residue)] = Double.parseDouble(aaBiasString[1]);
                    //NOTE: Only the relative proportion of the codon numbers means anything.  You can reweight
                    //linearly as much as you want.  E.g. it doesn't matter if you reweight the codon bias
                    //percentage so that it sums to 100 without the stop codon, which is not counted in the
                    //protein profile.
                } else {
                    System.out.println("Codon bias file is malformed.  Expecting two tab-delimited columns, but found " + fileLine);
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        //normalize the codon bias file
        double sum = 0;
        for (int i = 0; i < codonBias.length; i++) {
            double bias = codonBias[i];
            sum += bias;
        }
        for (int i = 0; i < codonBias.length; i++) {
            codonBias[i] = codonBias[i] / sum;
        }

        //calculate the shannon entropy of the distribution
        //the following entropy methods modified from DistributionLogo class in BioJava
        double information = Math.log(20) / Math.log(2.0);
        for (int i = 0; i < codonBias.length; i++) {
            double p = codonBias[i];
            if (p != 0.0) {
                double lp = Math.log(p);
                information -= -p * lp / Math.log(2.0);
            }
        }

        return (Math.pow(2,information));
    }

    /**
     * Output the specificity score for all profiles in a project file
     *
     * @param profileListFileName The list of profiles to analyze
     * @param codonBiasFileName   The code bias file to apply to the profiles
     */
    public void calculateSpecificityScoreForAllProfiles(String profileListFileName, String codonBiasFileName, int profileLength, ProteinTerminus terminus) {
        //calculate the specificity score correction factor from the codon bias (for PDZ, the -9 position, which is assumed to be random)
        double specificityScoreCorrectionFactor = calculateSpecificityScoreCorrectionFactor(codonBiasFileName);
        //calculate the specificity score for each domain at each position
        List profileList = PeptideToProfileReader.readPeptidesAsProfiles(new File(profileListFileName), profileLength, terminus,
                0.0, null, true);   //use the raw peptide data - no correction
        System.out.println(specificityScoreCorrectionFactor);
        System.out.println("PDZ\t-9\t-8\t-7\t-6\t-5\t-4\t-3\t-2\t-1\tFull Profile\tEffective Specificity Residues\tNumber of Peptides");
        for (int i = 0; i < profileList.size(); i++) {
            ProteinProfile proteinProfile = (ProteinProfile) profileList.get(i);
            System.out.print(proteinProfile.getName() + "\t");
            double profileSpecificityScore=1;
            //calculate the specificity score for each position
            for (int j = 0; j < proteinProfile.getNumColumns(); j++) {
                ProteinProfile subsetPosition = proteinProfile.getProfileSubsetCopy(Integer.toString(j));
                double correctedPositionSpecificityScore = subsetPosition.calculateSpecificityScore() - specificityScoreCorrectionFactor;
                if(correctedPositionSpecificityScore<1.0) {
                    correctedPositionSpecificityScore=1.0;
                }
                profileSpecificityScore *= correctedPositionSpecificityScore;
                System.out.print(correctedPositionSpecificityScore + "\t");
            }
            //print profile specificity, effective number of specificity residues and number of sequences
            System.out.println(profileSpecificityScore + "\t" + Math.log(profileSpecificityScore) / Math.log(20) + "\t" + proteinProfile.getNumSequences());  //log base 20
        }
    }

    /**
     * Currently, this is only here for testing purposes.
     * You can run the code without loading Cytoscape using this method.
     *
     * @param args
     */
    public static void main(String[] args) {
        BrainAlgorithm alg = new BrainAlgorithm();

        //1. write sequence logos - human
/*
        alg.writeSequenceLogos("/Users/bader/Gbader/Projects/Domains/PDZ-SidhuBoone/Analysis/BindingProfiles/Human/ProfileDataSets-ProjectFiles/ProjectFile.txt",
               "/Users/bader/Gbader/Projects/Domains/PDZ-SidhuBoone/Analysis/PeptideLogos/Human-9",
                "/Users/bader/Gbader/Projects/Domains/PDZ-SidhuBoone/Analysis/BindingProfiles/PhageCodonBias/phageLibraryPosition-9CodonBias.txt",
                0.0);
*/

        //2. write sequence logos - worm
/*
        alg.writeSequenceLogos("/Users/bader/Gbader/Projects/Domains/PDZ-SidhuBoone/Analysis/BindingProfiles/Worm/ProfileDataSets-ProjectFiles/ProjectFile.txt",
               "/Users/bader/Gbader/Projects/Domains/PDZ-SidhuBoone/Analysis/PeptideLogos/Worm",
                "/Users/bader/Gbader/Projects/Domains/PDZ-SidhuBoone/Analysis/BindingProfiles/PhageCodonBias/phageLibraryNNKTheoreticalCodonBias.txt",
                0.0);
*/

        //3. LogoTree
/*
        alg.runProfileCluster("/Users/bader/Gbader/Projects/Domains/PDZ-SidhuBoone/Analysis/BindingProfiles/ProjectFileWormAndHuman.txt", 10, ProteinTerminus.C,
                "Human and Worm PDZ Binding Profile Clustering - All residues, NNK correction, (STQN, KRH, DE, LAMIVFWY, C, P, G)",
                new File("/Users/bader/Gbader/Projects/Domains/PDZ-SidhuBoone/Analysis/LogoTree/GroupedDistributionDistance/LogoTree-GroupedDistributionDistance-WormAndHuman-AllResidues-NNKCorrected.pdf"),
                new File("/Users/bader/Gbader/Projects/Domains/PDZ-SidhuBoone/Analysis/BindingProfiles/Worm/WormPDZProteinNamesWithPhageData.txt"),
                new File("/Users/bader/Gbader/Projects/Domains/PDZ-SidhuBoone/Analysis/BindingProfiles/PhageCodonBias/phageLibraryNNKTheoreticalCodonBias.txt")
                );
*/

        //Logotree per position

/*
        alg.runProfileClusterPositionRange("/Users/bader/Gbader/Projects/Domains/PDZ-SidhuBoone/Analysis/BindingProfiles/ProjectFileWormAndHuman.txt", 0,9, 10, ProteinTerminus.C,
                "Human+Worm PDZ Binding Profile Clustering - NNK bias correction",
                new File("/Users/bader/Gbader/Projects/Domains/PDZ-SidhuBoone/Analysis/LogoTree/GroupedDistributionDistanceByPosition/LogoTree-GroupedDistributionDistance-WormAndHuman-AllResidues-NNKCorrected"),
                new File("/Users/bader/Gbader/Projects/Domains/PDZ-SidhuBoone/Analysis/BindingProfiles/Worm/WormPDZProteinNamesWithPhageData.txt"),
                new File("/Users/bader/Gbader/Projects/Domains/PDZ-SidhuBoone/Analysis/BindingProfiles/PhageCodonBias/phageLibraryNNKTheoreticalCodonBias.txt"));
*/


        //4. find heavily used AAs for classification aa grouping choice
//        alg.runPeptideAAByPositionHistogramAnalysis("/Users/bader/Gbader/Projects/Domains/PDZ-SidhuBoone/Analysis/BindingProfiles/ProjectFileWormAndHuman.txt", 10);

        //5. bootstrap analysis on profile clustering
/*
        HierarchicalClusteringBootstrapAnalysis.runInputOrderRobustnessTest("/Users/bader/Gbader/Projects/Domains/PDZ-SidhuBoone/Analysis/BindingProfiles/ProjectFileWormAndHuman-filtered.txt", 10,
                ProteinTerminus.C, 10000, "10000 samplings, All Residues (STQN, KRH, DE, FLAMPWIVCY, G) - NNK Correction, human+worm",
                new File("/Users/bader/Gbader/Projects/Domains/PDZ-SidhuBoone/Analysis/LogoTree/GroupedDistributionDistance/inputOrderBootstrap/LogoTree-GroupedDistributionDistance-WormAndHuman-AllResidues_NNKCorrection-filtered-bootstrap.pdf"),
                new File("/Users/bader/Gbader/Projects/Domains/PDZ-SidhuBoone/Analysis/BindingProfiles/Worm/WormPDZProteinNamesWithPhageData.txt"),
                new File("/Users/bader/Gbader/Projects/Domains/PDZ-SidhuBoone/Analysis/BindingProfiles/PhageCodonBias/phageLibraryNNKTheoreticalCodonBias.txt"));
*/

        //6. TODO Run protein sequence space coverage analysis
/*
        try {
            //human
            alg.runProteinSequenceSpaceCoverageAnalysis(new File("D:\\Gbader\\Code\\PDZ\\data\\PDZ\\DomainDefinitionAndAlignments\\Human\\pdz_human_unique_filtered.muscle.fa"),
                    new File("D:\\Gbader\\Code\\PDZ\\data\\PDZ\\DomainDefinitionAndAlignments\\Human\\pdz_human_unique_filtered_phageResultsOnly.txt"),
                    new File("D:\\Gbader\\Code\\PDZ\\data\\PDZ\\DomainDefinitionAndAlignments\\Human\\pdz_human_unique_filtered_phageFailuresOnly.txt"),
                    null, 100);
            //"19,26-40,50-54,85-89,122,125-127,130-133", 100); //4 Angstrom
            //"27,31-34,50,53-54,85-88,122,127,130-133", 100); //2 Angstrom

            //worm
            //alg.runProteinSequenceSpaceCoverageAnalysis(new File("D:\\Gbader\\Code\\PDZ\\data\\PDZ\\DomainDefinitionAndAlignments\\Worm\\pdz_worm_unique.muscle.fa"),
            //        new File("D:\\Gbader\\Code\\PDZ\\data\\PDZ\\DomainDefinitionAndAlignments\\Worm\\pdz_worm_unique_phageResultsOnly.txt"),
            //        null, 100);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (BioException e) {
            e.printStackTrace();
        }
*/

        //7. Calculate specificity scores for all profiles

        /*
        alg.calculateSpecificityScoreForAllProfiles("/Users/bader/Gbader/Projects/Domains/PDZ-SidhuBoone/Analysis/BindingProfiles/ProjectFileWormAndHuman.txt",
                "/Users/bader/Gbader/Projects/Domains/PDZ-SidhuBoone/Analysis/BindingProfiles/PhageCodonBias/phageLibraryPosition-9CodonBias.txt",
                10, ProteinTerminus.C);
        */

        /*
        //residue-residue correlation analysis
        ResidueResidueCorrelationMatrix rrcm = null;
        try {
            //rrcm = new ResidueResidueCorrelationMatrix(new File("D:\\Gbader\\Code\\PDZ\\data\\PDZ\\DomainDefinitionAndAlignments\\HumanWorm\\pdz_human_worm_unique.muscle.fa"),
            rrcm = new ResidueResidueCorrelationMatrix(new File("D:\\Gbader\\Code\\PDZ\\data\\PDZ\\DomainDefinitionAndAlignments\\HumanWorm\\pdz_human_worm_unique.superfamily.fa"),
                    new File("D:\\Gbader\\Code\\PDZ\\data\\PDZ\\BindingProfiles\\ProjectFileWormAndHuman.txt"),
                    //new File("D:\\Gbader\\Code\\PDZ\\data\\PDZ\\BindingProfiles\\Human\\Final Formatted Data\\ERBB2IP-1-hi.pep.txt"),
                    10, ProteinTerminus.C);
            //rrcm.setDomainSequenceFilter("21,26,29,33,60,61,64,67,72,78,82,89,98,102,123"); //15 specificity residues from proteinkeys.org (index starts at 0 on human-worm alignment)
            long counts = rrcm.learn(1, 1);
            rrcm.printMostInformativeFeatures("ERBB2IP-1", -0.08);
        } catch (BioException e) {
            e.printStackTrace();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
        */

        //bootstrap analysis on profile clustering
        /*
        HierarchicalClusteringBootstrapAnalysis.runInputOrderRobustnessTest("D:\\Gbader\\Code\\PDZ\\data\\PDZ\\BindingProfiles\\ProjectFileWormAndHuman.txt", 4,
                ProteinTerminus.C, 10000, "10000 samplings, Last 4 Residues PSG - No Bias, human+worm",
                new File("D:\\Gbader\\Code\\PDZ\\data\\PDZ\\PeptideProfileClustering\\GroupedByPositionDistributionDistance\\inputOrderBootstrap\\test.pdf"));
        */

        //Conserved hit analysis
        /*
        try {
            DomainInteractionAnalysis.runConservedHitAnalysis("D:\\Gbader\\Code\\PDZ\\data\\PDZ\\PPI Prediction\\Inparanoid\\downloadedInparanoid\\projectFile.txt", null);
        } catch (IOException e) {
            e.printStackTrace();
        }
        */

        //conserved link analysis

        /*
        try {
            DomainInteractionAnalysis.runConservedLinkAnalysis("D:\\Gbader\\Code\\PDZ\\data\\PDZ\\PPI Prediction\\Inparanoid\\hs-vs-ce-Inparanoid\\projectFile.txt",
                    "D:\\Gbader\\Code\\PDZ\\data\\PDZ\\BindingProfiles\\Worm\\Final Formatted Data\\ProjectFile.txt",
                    "D:\\Gbader\\Code\\PDZ\\data\\PDZ\\PPI Prediction\\Inparanoid\\hs-vs-ce-Inparanoid\\wp147CE",
                    "fasta",
                    null);
        } catch (IOException e) {
            e.printStackTrace();
        }
        */

        /*
        try {
        DomainInteractionAnalysis.scoreCoExpressionByAPCall(new File("D:\\Gbader\\Code\\PDZ\\data\\PDZ\\PPI Prediction\\Co-expression\\Gene atlas\\GSE1133_U133A with AP calls.soft"),
        new File("D:\\Gbader\\Code\\PDZ\\data\\PDZ\\PPI Prediction\\Co-expression\\Gene atlas\\GSE1133_U133A with AP calls_filteredSampleList.TXT"),
        new File("D:\\Gbader\\Code\\PDZ\\data\\PDZ\\PPI Prediction\\Co-expression\\Gene atlas\\refseq2affyU133A_Aug23_2005_filtered.txt"), net);
        } catch (IOException e) {
        e.printStackTrace();
        }
        */

        //expression data test

        /*        alg.runConservedLinkAnalysis(
                        "D:\\Gbader\\Code\\PDZ\\data\\PDZ\\BindingProfiles\\Human\\Final Formatted Data\\TestProjectFile.txt",
        //                "D:\\Gbader\\Code\\PDZ\\data\\PDZ\\BindingProfiles\\Human\\Final Formatted Data\\ERBB2IP-1-hi.pep.txt",
                        "D:\\Gbader\\Code\\PDZ\\data\\DBs\\refseq_May16_2005_human.protein.gpff",
                        "D:\\Gbader\\Code\\PDZ\\data\\PDZ\\BindingProfiles\\Worm\\Final Formatted Data\\ProjectFile.txt",
                        "D:\\Gbader\\Code\\PDZ\\data\\DBs\\refseq_Jun9_2005_worm.protein.gpff",
                        "D:\\Gbader\\Code\\PDZ\\data\\DBs\\homologeneMay26_2005.data");*/

        //Profile distance analysis
        //alg.runAllVsAllProfileDistance("/Users/bader/Gbader/Projects/Domains/PDZ-SidhuBoone/Analysis/BindingProfiles/ProjectFileWormAndHuman.txt");
        //alg.runProfileCluster("D:\\Gbader\\Code\\PDZ\\data\\PDZ\\BindingProfiles\\Human\\ProfileDataSets-ProjectFiles\\ProjectFile.txt", 4, ProteinTerminus.C);
        //alg.runProfileCluster("D:\\Gbader\\Code\\PDZ\\data\\PDZ\\BindingProfiles\\Worm\\Final Formatted Data\\ProjectFile.txt");

        /*
        alg.runProfileClusterPositionRange("D:\\Gbader\\Code\\PDZ\\data\\PDZ\\BindingProfiles\\ProjectFileWormAndHuman.txt", 6,9, 4, ProteinTerminus.C,
                "Human+Worm PDZ Binding Profile Clustering - NNK bias correction",
                new File("D:\\Gbader\\Code\\PDZ\\data\\PDZ\\PeptideProfileClustering\\GroupedDistributionDistance\\last4positions\\GroupedHydrophobicBySize-NNK\\LogoTree-GroupedDistributionDistance-WormAndHuman-NNKBiasCorrection.pdf"),
                new File("D:\\Gbader\\Code\\PDZ\\data\\PDZ\\BindingProfiles\\Worm\\WormPDZProteinNamesWithPhageData.txt"),
                new File("D:\\Gbader\\Code\\PDZ\\data\\PDZ\\BindingProfiles\\Human\\PhageCodonBias\\phageLibraryNNKTheoreticalCodonBias.txt"));
        alg.runProfileClusterPosition("D:\\Gbader\\Code\\PDZ\\data\\PDZ\\BindingProfiles\\ProjectFileWormAndHuman.txt", "6,8", 4, ProteinTerminus.C,
                "Human+Worm PDZ Binding Profile Clustering - NNK bias correction",
                new File("D:\\Gbader\\Code\\PDZ\\data\\PDZ\\PeptideProfileClustering\\GroupedDistributionDistance\\last4positions\\GroupedHydrophobicBySize-NNK\\LogoTree-GroupedDistributionDistance-WormAndHuman-NNKBiasCorrection.pdf"),
                new File("D:\\Gbader\\Code\\PDZ\\data\\PDZ\\BindingProfiles\\Worm\\WormPDZProteinNamesWithPhageData.txt"),
                new File("D:\\Gbader\\Code\\PDZ\\data\\PDZ\\BindingProfiles\\Human\\PhageCodonBias\\phageLibraryNNKTheoreticalCodonBias.txt"));
        */
        //Erbin point mutant analysis - Sep.13.2006
/*        alg.runProfileCluster("/Users/bader/Gbader/Projects/Domains/PDZ-SidhuBoone/Analysis/BindingProfiles/Human/2006Sep-Erbin-point-mutants03/Mutants_NEW-Textfiles/projectfile.txt", 4, ProteinTerminus.C,
                "Erbin mutant panel (New) - NNK bias correction - Cluster last 4 positions (STQN, KRH, DE, FLAMPWIVCY, G)",
                new File("/Users/bader/Gbader/Projects/Domains/PDZ-SidhuBoone/Analysis/BindingProfiles/Human/2006Sep-Erbin-point-mutants03/Mutants_NEW-Textfiles/LogoTree-GroupedDistributionDistance-ErbinMutantNew-NNK-last4.pdf"),
                null, //no leaf label highlights
                new File("/Users/bader/Gbader/Projects/Domains/PDZ-SidhuBoone/Analysis/BindingProfiles/PhageCodonBias/phageLibraryNNKTheoreticalCodonBias.txt"));*/
/*        alg.runProfileClusterPositionRange("/Users/bader/Gbader/Projects/Domains/PDZ-SidhuBoone/Analysis/BindingProfiles/Human/2006Sep-Erbin-point-mutants03/Mutants_NEW-Textfiles/projectfile.txt", 0,3, 4, ProteinTerminus.C,
                "Erbin mutant panel (New) - NNK bias correction - Cluster last 4 positions",
                new File("/Users/bader/Gbader/Projects/Domains/PDZ-SidhuBoone/Analysis/BindingProfiles/Human/2006Sep-Erbin-point-mutants03/Mutants_NEW-Textfiles/LogoTree-GroupedDistributionDistance-ErbinMutantNew-NNK-last4.pdf"),
                null, //no leaf label highlights
                new File("/Users/bader/Gbader/Projects/Domains/PDZ-SidhuBoone/Analysis/BindingProfiles/PhageCodonBias/phageLibraryNNKTheoreticalCodonBias.txt"));*/
/*        alg.runProfileClusterPosition("/Users/bader/Gbader/Projects/Domains/PDZ-SidhuBoone/Analysis/BindingProfiles/Human/2006Sep-Erbin-point-mutants03/Mutants_NEW-Textfiles/projectfile.txt", "0", 4, ProteinTerminus.C,
                "Erbin mutant panel (New) - NNK bias correction - Cluster position -3",
                new File("/Users/bader/Gbader/Projects/Domains/PDZ-SidhuBoone/Analysis/BindingProfiles/Human/2006Sep-Erbin-point-mutants03/Mutants_NEW-Textfiles/LogoTree-GroupedDistributionDistance-ErbinMutantNew-NNK-position-3.pdf"),
                null, //no leaf label highlights
                new File("/Users/bader/Gbader/Projects/Domains/PDZ-SidhuBoone/Analysis/BindingProfiles/PhageCodonBias/phageLibraryNNKTheoreticalCodonBias.txt"));
        alg.runProfileClusterPosition("/Users/bader/Gbader/Projects/Domains/PDZ-SidhuBoone/Analysis/BindingProfiles/Human/2006Sep-Erbin-point-mutants03/Mutants_NEW-Textfiles/projectfile.txt", "1", 4, ProteinTerminus.C,
                "Erbin mutant panel (New) - NNK bias correction - Cluster position -2",
                new File("/Users/bader/Gbader/Projects/Domains/PDZ-SidhuBoone/Analysis/BindingProfiles/Human/2006Sep-Erbin-point-mutants03/Mutants_NEW-Textfiles/LogoTree-GroupedDistributionDistance-ErbinMutantNew-NNK-position-2.pdf"),
                null, //no leaf label highlights
                new File("/Users/bader/Gbader/Projects/Domains/PDZ-SidhuBoone/Analysis/BindingProfiles/PhageCodonBias/phageLibraryNNKTheoreticalCodonBias.txt"));
        alg.runProfileClusterPosition("/Users/bader/Gbader/Projects/Domains/PDZ-SidhuBoone/Analysis/BindingProfiles/Human/2006Sep-Erbin-point-mutants03/Mutants_NEW-Textfiles/projectfile.txt", "2", 4, ProteinTerminus.C,
                "Erbin mutant panel (New) - NNK bias correction - Cluster position -1",
                new File("/Users/bader/Gbader/Projects/Domains/PDZ-SidhuBoone/Analysis/BindingProfiles/Human/2006Sep-Erbin-point-mutants03/Mutants_NEW-Textfiles/LogoTree-GroupedDistributionDistance-ErbinMutantNew-NNK-position-1.pdf"),
                null, //no leaf label highlights
                new File("/Users/bader/Gbader/Projects/Domains/PDZ-SidhuBoone/Analysis/BindingProfiles/PhageCodonBias/phageLibraryNNKTheoreticalCodonBias.txt"));
        alg.runProfileClusterPosition("/Users/bader/Gbader/Projects/Domains/PDZ-SidhuBoone/Analysis/BindingProfiles/Human/2006Sep-Erbin-point-mutants03/Mutants_NEW-Textfiles/projectfile.txt", "3", 4, ProteinTerminus.C,
                "Erbin mutant panel (New) - NNK bias correction - Cluster position 0",
                new File("/Users/bader/Gbader/Projects/Domains/PDZ-SidhuBoone/Analysis/BindingProfiles/Human/2006Sep-Erbin-point-mutants03/Mutants_NEW-Textfiles/LogoTree-GroupedDistributionDistance-ErbinMutantNew-NNK-position-0.pdf"),
                null, //no leaf label highlights
                new File("/Users/bader/Gbader/Projects/Domains/PDZ-SidhuBoone/Analysis/BindingProfiles/PhageCodonBias/phageLibraryNNKTheoreticalCodonBias.txt"));*/

/*
        alg.runProfileCluster("/Users/bader/Gbader/Projects/Domains/PDZ-SidhuBoone/Analysis/BindingProfiles/Human/2006Sep-Erbin-point-mutants03/Mutants_OLD-Textfiles/projectfile.txt", 4, ProteinTerminus.C,
                "Erbin mutant panel (Old) - NNK bias correction - Cluster last 4 positions (STQN, KRH, DE, FLAMPWIVCY, G)",
                new File("/Users/bader/Gbader/Projects/Domains/PDZ-SidhuBoone/Analysis/BindingProfiles/Human/2006Sep-Erbin-point-mutants03/Mutants_OLD-Textfiles/LogoTree-GroupedDistributionDistance-ErbinMutantOld-NNK-last4.pdf"),
                null, //no leaf label highlights
                new File("/Users/bader/Gbader/Projects/Domains/PDZ-SidhuBoone/Analysis/BindingProfiles/PhageCodonBias/phageLibraryNNKTheoreticalCodonBias.txt"));
        alg.runProfileClusterPositionRange("/Users/bader/Gbader/Projects/Domains/PDZ-SidhuBoone/Analysis/BindingProfiles/Human/2006Sep-Erbin-point-mutants03/Mutants_OLD-Textfiles/projectfile.txt", 0,3, 4, ProteinTerminus.C,
                "Erbin mutant panel (Old) - NNK bias correction",
                new File("/Users/bader/Gbader/Projects/Domains/PDZ-SidhuBoone/Analysis/BindingProfiles/Human/2006Sep-Erbin-point-mutants03/Mutants_OLD-Textfiles/LogoTree-GroupedDistributionDistance-ErbinMutantOld-NNK.pdf"),
                null, //no leaf label highlights
                new File("/Users/bader/Gbader/Projects/Domains/PDZ-SidhuBoone/Analysis/BindingProfiles/PhageCodonBias/phageLibraryNNKTheoreticalCodonBias.txt"));
*/

        //alg.runProfileCluster("D:\\Gbader\\Code\\PDZ\\data\\PDZ\\BindingProfiles\\Human\\ProfileDataSets-ProjectFiles\\TempTestProjectFile.txt", 4, ProteinTerminus.C);

//        alg.runProfileCluster("", "SH3 Domain Logo Tree", new File(""));

        //find correlations at position -1,-3 to choose classification aa groupings
        //alg.runPeptideAAByPositionPairHistogramAnalysis("D:\\Gbader\\Code\\PDZ\\data\\PDZ\\BindingProfiles\\ProjectFileWormAndHuman.txt", 10, 6, 8);

        //run GO ORA analysis
        //GoORAAnalysis(alg);

        /*
        //output sequence search results to the screen
        Collection results = searchResults.getAllResultSets();
        for (Iterator iterator = results.iterator(); iterator.hasNext();) {
            SequenceSearchResultSet sequenceSearchResultSet = (SequenceSearchResultSet) iterator.next();
            System.out.println(sequenceSearchResultSet);
        }
        */
        /*
        CyNetwork net = CytoscapeUtil.addProfileSearchResultsToCytoscape(searchResults, alg.getParams());
        try {
            DomainInteractionAnalysis.scoreCoExpressionByAPCall(new File("D:\\Gbader\\Code\\PDZ\\data\\PDZ\\PPI Prediction\\Co-expression\\Gene atlas\\GSE1133_samples_full.txt"),
                    new File("D:\\Gbader\\Code\\PDZ\\data\\PDZ\\PPI Prediction\\Co-expression\\Gene atlas\\GSE1133 tissue list-filtered146.txtt"),
                    new File("D:\\Gbader\\Code\\PDZ\\data\\PDZ\\PPI Prediction\\Co-expression\\Gene atlas\\refseq2affyU133A_Aug23_2005_filtered.txt"), net);
        } catch (IOException e) {
            e.printStackTrace();
        }
        */

        //Profile distance analysis
        //alg.runAllVsAllProfileDistancePairOutput("D:\\Gbader\\Code\\PDZ\\data\\PDZ\\BindingProfiles\\Human\\Final Formatted Data\\ProjectFile.txt");

        //aligned sequence similarity
        /*                AlignedProteinSequenceSimilarity apss=null;
                        try {
                            apss = new AlignedProteinSequenceSimilarity("D:\\Gbader\\Code\\PDZ\\data\\PDZ\\DomainDefinitionAndAlignments\\Human\\pdz_human_unique_filtered_phageResultsOnly.fa",
                                    "15,17,21-23,29,32,33,47-50,82,86,89,90"); //2 A
                        } catch (IOException e) {
                            e.printStackTrace();
                        }
                        apss.calculateSimilarityMatrix();*/

        //histogram analysis
        //worm
        //alg.runHistogramScoreAnalysis("D:\\Gbader\\Code\\PDZ\\data\\PDZ\\BindingProfiles\\Worm\\Final Formatted Data\\ProjectFile.txt",
        //        "D:\\Gbader\\Code\\PDZ\\data\\DBs\\2005Jun10_yeast_nrpep.fasta", "FASTA", ProteinTerminus.C);
        //human
        //alg.runHistogramScoreAnalysis("D:\\Gbader\\Code\\PDZ\\data\\PDZ\\BindingProfiles\\Human\\Final Formatted Data\\ProjectFile.txt",
        //        "D:\\Gbader\\Code\\PDZ\\data\\DBs\\refseq_May16_2005_human.protein.gpff", "genpept", ProteinTerminus.C);

    }

    private static void GoORAAnalysis(BrainAlgorithm alg) {
        BrainParameterSet params = new BrainParameterSet();
        params.setProfileFileName(new File("D:\\Gbader\\Code\\PDZ\\data\\PDZ\\BindingProfiles\\Human\\ProfileDataSets-ProjectFiles\\ProjectFile.txt"));
        params.setFuzzFactor(1.0);
        params.setDatabaseFileName(new File("D:\\Gbader\\Code\\PDZ\\data\\DBs\\refseq_May16_2005_human.protein.gpff"));
        params.setDatabaseFormat("genpept");
        ProteinDatabaseSearchParams dbparams = new ProteinDatabaseSearchParams(ProteinTerminus.C);
        dbparams.setNormalized(true);
        params.setSearchParams(dbparams);
        alg.setParams(params);
        System.out.println("Searching...");
        MultiSequenceSearchResultSet searchResults = alg.runProfileSearch(null, alg.findAutoScoreThreshold(params), params);

        //read profile file - could be a project or single profile (list of peptides)
        System.out.println("Clustering...");
        List clusterProfileList = null;
        List cutProfileList = null;
        //full profiles are always stored so we can draw a logo tree
        params.setFuzzFactor(0.0);
        List originalProfileList = PeptideToProfileReader.readPeptidesAsProfiles(params.getProfileFile(), 0, null, params.getFuzzFactor(),
                new File("D:\\Gbader\\Code\\PDZ\\data\\PDZ\\BindingProfiles\\Human\\PhageCodonBias\\phageLibraryNNKTheoreticalCodonBias.txt"), true);
        clusterProfileList = originalProfileList;
        int profileLength = 4;
        if (profileLength > 0) {
            //profile length specified, cluster profiles of this length (not full length ones)
            cutProfileList = new ArrayList();
            for (int i = 0; i < originalProfileList.size(); i++) {
                ProteinProfile proteinProfile = (ProteinProfile) originalProfileList.get(i);
                cutProfileList.add(proteinProfile.getTruncatedProfileCopy(profileLength, ProteinTerminus.C));
            }
            clusterProfileList = cutProfileList;
        }
        DistanceMatrix distanceMatrix = new DistanceMatrix(clusterProfileList.size());
        distanceMatrix.calcDistances(clusterProfileList, new DistanceMetric() {
            public double calc(Object object1, Object object2) {
                ProteinProfile proteinProfile1 = (ProteinProfile) object1;
                ProteinProfile proteinProfile2 = (ProteinProfile) object2;
                String groupingByPosition = AminoAcidGrouping.getPolarChargedHydrophobeGrouping();
                return (ProteinProfileDistance.calculateAAGroupedDistributionDistance(proteinProfile1, proteinProfile2, groupingByPosition));
            }
        });
        ArrayList al = new ArrayList(clusterProfileList.size());
        for (int i = 0; i < clusterProfileList.size(); i++) {
            ProteinProfile proteinProfile = (ProteinProfile) clusterProfileList.get(i);
            al.add(i, proteinProfile.getName());
        }
        distanceMatrix.setLabels(al);
        //cluster
        AvgLinkHierarchicalClustering cluster = new AvgLinkHierarchicalClustering(distanceMatrix);
        cluster.setOptimalLeafOrdering(false);
        cluster.run();
        System.out.println("Running ORA...");
        alg.findOverrepresentedGOTermsForAllPredictedPPIsForAllPossibleProfileClusters(cluster, searchResults);
    }
}
