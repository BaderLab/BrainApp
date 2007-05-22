package org.baderlab.brain;

import mt.MatrixEntry;
import org.biojava.bio.BioException;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.DistributionTools;
import org.biojava.bio.dist.IndexedCount;
import org.biojava.bio.dp.SimpleWeightMatrix;
import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.io.SeqIOTools;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;
import org.biojava.utils.ChangeVetoException;
import smt.FlexCompRowMatrix;
import smt.SparseVector;

import javax.imageio.ImageIO;
import java.io.*;
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
 * * Date: Nov 1, 2005
 * * Time: 6:07:33 PM
 */

/**
 * Handles building and storing a protein-protein residue correlation matrix. Input is a set of aligned proteins
 * and a set of peptides/proteins known to interact with each protein in the multiple sequence alignment. Assumes
 * that one or more peptides/proteins are known to interact.
 */
public class ResidueResidueCorrelationMatrix {
    private FlexCompRowMatrix correlationMatrix = null; //the actual correlation matrix
    private int totalDomainSequenceLength = -1;  //the length of the domain MSA (all domain sequences must be the same length)
    private int totalPeptideSequenceLength = -1; //the length of the peptides (all peptides must be the same length)
    private HashMap profileNameToProfile = null; //in memory database of profiles to learn from
    private HashMap sequenceNameToSequence = null; //note: sequence name must equal profile name for matches to be made
    private boolean domainPositionsAsRows = false; //optimization for sparsematrix
    private int numDomainPositionsPerFeature = 0; //set when correlations are learned
    private int numPeptidePositionsPerFeature = 0; //set when correlations are learned
    private int multipleSequenceAlignmentLength = 0; //the length of the effective multiple sequence alignment that
    //includes aligned proteins and bound peptides/proteins (basically, the number of peptides, uniqued by domain)
    private SparseVector domainFeatureFrequencyVector = null; //stores the frequency of the domain features in the MSA
    private SparseVector peptideFeatureFrequencyVector = null; //stores the frequency of the peptide features in the MSA

    //used for sequence filtering functionality
    private String domainSequenceFilter = null; //a string representing the columns to consider from a domain sequence e.g. 1, 2..5
    private String peptideSequenceFilter = null; //a string representing the columns to consider from a peptide sequence

    private final String aaList = "ACDEFGHIKLMNPQRSTVWY"; //reference 20 letter amino acid alphabet

    /**
     * Creates a new residue-residue correlation matrix
     *
     * @param multipleSequenceAlignmentFile The multiple sequence alignment (MSA) of the domains for this matrix
     *                                      All domain sequences must be the same length.
     * @param peptideOrProjectFile          The peptide or project file specifying the peptides for this matrix
     *                                      Important note: The MSA and peptide files must use corresponding protein names e.g. proteinX in the MSA
     *                                      must also be called proteinX in the peptide file.
     * @param peptideLength                 The length of the peptides to read in.
     *                                      If <0 e.g. -1, then use the maximum peptide length found in the peptide files.
     * @param terminus                      Terminus to count from if peptide length is specified (>0)
     * @throws BioException          If the multiple sequence alignment can't be read
     * @throws FileNotFoundException If there is a problem with the peptide or project file
     */
    public ResidueResidueCorrelationMatrix(File multipleSequenceAlignmentFile, File peptideOrProjectFile, int peptideLength, ProteinTerminus terminus) throws BioException, IOException {
        List proteinProfileList = null;

        readAlignmentAndDetermineSequenceAlignmentWidth(multipleSequenceAlignmentFile);

        //read in peptide files - only consider unique peptides for the correlation matrix.
        //Fuzz factor and codon bias have no meaning, since we're just using the profile reader to get at the peptides
        proteinProfileList = PeptideToProfileReader.readPeptidesAsProfiles(peptideOrProjectFile, peptideLength, terminus, 0.0, null, true);
        profileNameToProfile = new HashMap();
        multipleSequenceAlignmentLength = 0;
        for (int i = 0; i < proteinProfileList.size(); i++) {
            ProteinProfile proteinProfile = (ProteinProfile) proteinProfileList.get(i);
            //check that all peptides are same length
            if (totalPeptideSequenceLength < 0) {
                totalPeptideSequenceLength = proteinProfile.getNumColumns();
            } else if (totalPeptideSequenceLength != proteinProfile.getNumColumns()) {
                throw new IllegalArgumentException("All peptide sequences must be the same length across all files." +
                        "Found a profile of length " + proteinProfile.getNumColumns() + " in " + proteinProfile.getName() +
                        " but was expecting length " + totalPeptideSequenceLength + " (based on the length of the" +
                        " first sequence seen).");
            }
            //save the profiles for easy lookup by name
            profileNameToProfile.put(proteinProfile.getName(), proteinProfile);
            //count the number of peptides (the length of the total effective MSA)
            multipleSequenceAlignmentLength += proteinProfile.getNumSequences();
        }
    }

    /**
     * Reads the alignment into memory and validates it.
     * Also, determines the width of the protein multiple sequence alignment (MSA). All sequences in the MSA must be
     * the same width (a.k.a. length). This method will throw an exception if that is not the case. Sets a global variable for this
     * class to store the MSA width
     *
     * @param multipleSequenceAlignmentFile The MSA file to cycle through, check and store in memory
     * @throws IOException  If there is an error with MSA file I/O
     * @throws BioException If there is an error with the MSA file format
     */
    private void readAlignmentAndDetermineSequenceAlignmentWidth(File multipleSequenceAlignmentFile) throws IOException, BioException {
        sequenceNameToSequence = new HashMap();
        //read protein domain multiple sequence alignment
        BufferedReader brMSA = new BufferedReader(new FileReader(multipleSequenceAlignmentFile));
        SequenceIterator sequenceAlignment = (SequenceIterator) SeqIOTools.fileToBiojava("fasta", "PROTEIN", brMSA);
        //determine the width of the MSA
        while (sequenceAlignment.hasNext()) {
            Sequence seq = sequenceAlignment.nextSequence();
            //check that all MSA sequences are the same length
            if (totalDomainSequenceLength < 0) {
                totalDomainSequenceLength = seq.length();
            } else if (totalDomainSequenceLength != seq.length()) {
                throw new IllegalArgumentException("All domain sequences must be the same length." +
                        "Found a domain sequence of length " + seq.length() + " called " + seq.getName() +
                        " but was expecting length " + totalDomainSequenceLength + " (based on the length of the" +
                        " first sequence seen).");
            }
            sequenceNameToSequence.put(seq.getName(), seq);
        }
        brMSA.close();
    }

    /**
     * Stores a single feature that we learn from (position, residue pair)
     */
    private class ResiduePositionPair {
        public char residue;
        public int position;  //positions start at 0

        public ResiduePositionPair(int position, char residue) {
            this.position = position;
            this.residue = residue;
        }

        public String toString() {
            return new String("(" + residue + "@" + position + ")");
        }
    }

    /**
     * Calculates the number of combinations of unordered k out of a set of n numbers (n choose k)
     * n! / (k! (n-k)!)
     *
     * @param n The size of the set
     * @param k The number of things to choose from the set
     * @return n choose k
     */
    private int getCombinations(int n, int k) {
        if (k == 1) {
            return n; //trivial case
        }
        int k1 = k;
        int k2 = n - k;
        //computes (n! / k1! k2!) for combinations
        //assure k1 >= k2
        if (k1 < k2) {
            k1 = n - k;
            k2 = k;
        }
        int t = 0;
        if (k1 > n)
            return t;
        else {
            // accumulate the factors for k2 factorial (k2 is smallest)
            t = 1;
            while (k2 > 1)
                t *= k2--;

            // accumulate the factors from n downto k1 (k1 is largest)
            int t2 = 1;
            while (n > k1)
                t2 *= n--;

            t = t2 / t;
        }
        return t;
    }

    /**
     * Finds the max matrix index for a given sequence length and 'number of positions chosen' value
     *
     * @param numberOfPositions The number of positions in the sequence considered
     * @param totalPositions    The total number of positions in the sequence
     * @return The maximum matrix index
     */
    private int getMaxIndex(int numberOfPositions, int totalPositions) {
        int maxIndex = 0;

        maxIndex = getCombinations(totalPositions, numberOfPositions);
        maxIndex *= Math.pow(20, numberOfPositions);

        return maxIndex;
    }

    /**
     * Finds the index in the matrix for this feature
     *
     * @param feature The feature to look up in the matrix. This can be one or more residue/position pairs
     * @return The matrix index for this feature
     */
    private int getFeatureIndex(ResiduePositionPair[] feature) {
        int index = 0;
        //partitions of length aaList.length()^feature.length
        int partitionNumber = 0;
        int partitionOffset = 0;
        int residueOffset = 0;

        //matrix is organized in (n choose k) partitions of (size of aa alphabet ^ number of positions)
        //where n = max number of positions and k is the number of positions in the feature
        //IMPORTANT: the partition formula assumes that pair(i).position < pair(i+1).position
        //note: partition formula is a perfect hashing function based on Pascal's triangle
        for (int i = 0; i < feature.length; i++) {
            //partition number - maxes out at the number of total positions choose feature.length
            partitionNumber += getCombinations(feature[i].position, i + 1);
        }
        //partitions of length aaList.length()^feature.length
        partitionOffset = partitionNumber * (int) Math.pow(aaList.length(), feature.length);
        for (int i = 0; i < feature.length; i++) {
            residueOffset += aaList.indexOf(feature[i].residue) * Math.pow(aaList.length(), i);
        }
        index = partitionOffset + residueOffset;

        return index;
    }

    /**
     * Converts an index generated with getFeatureIndex back to a feature
     *
     * @param index          The index to convert
     * @param feature        The allocated feature array to fill
     * @param sequenceLength The total sequence length (maximum feature position+1)
     * @return the calculated feature
     */
    private ResiduePositionPair[] indexToFeature(int index, ResiduePositionPair[] feature, int sequenceLength) {
        //partitions of length aaList.length()^feature.length
        int partitionNumber = 0;
        int partitionOffset = 0;

        //dial up positions
        //determine partition number
        partitionNumber = (int) Math.floor(index / Math.pow(aaList.length(), feature.length));
        int tempPartitionNumber = partitionNumber;
        for (int featurePosition = feature.length; featurePosition > 0; featurePosition--)
        {  //featurePosition starts at 1
            for (int i = (featurePosition - 1); i < sequenceLength; i++) {
                if (tempPartitionNumber >= getCombinations(i + 1, featurePosition)) {
                    continue;
                }
                feature[featurePosition - 1].position = i;
                tempPartitionNumber -= getCombinations(i, featurePosition);
                break;
            }
        }

        //dial up residues
        //partitions of length aaList.length()^feature.length
        partitionOffset = partitionNumber * (int) Math.pow(aaList.length(), feature.length);
        int residueOffset = index - partitionOffset;
        for (int i = feature.length - 1; i >= 0; i--) {
            feature[i].residue = aaList.charAt(residueOffset / (int) Math.pow(aaList.length(), i));
            residueOffset %= (int) Math.pow(aaList.length(), i);
        }
        return feature;
    }

    /**
     * Traverse through the sequence and all binding peptides and add features to the correlation matrix.
     * Note: the sequence filter, if set, is evaluated in this method.
     *
     * @param domainSequence                The sequence of the domain
     * @param bindingPeptides               The protein profile containing all of the binding sequences
     * @param numDomainPositionsPerFeature  The number of domain positions to consider as a feature
     * @param numPeptidePositionsPerFeature The number of peptide positions to consider as a feature
     * @return The number of counts added to the matrix
     */
    private int learnSequenceToPeptides(Sequence domainSequence, ProteinProfile bindingPeptides, int numDomainPositionsPerFeature, int numPeptidePositionsPerFeature) {
        int correlationCounts = 0;
        String domainSequenceString = null;
        if (domainSequenceFilter == null) {
            domainSequenceString = domainSequence.seqString();
        } else {
            domainSequenceString = ProteinSequenceUtil.filterSequenceByColumns(domainSequence, domainSequenceFilter);
        }
        Collection peptides = bindingPeptides.getSequenceMap();
        String peptideSequenceString = null;
        for (Iterator iterator = peptides.iterator(); iterator.hasNext();) {
            Sequence peptideSequence = (Sequence) iterator.next();
            if (peptideSequenceFilter == null) {
                peptideSequenceString = peptideSequence.seqString();
            } else {
                peptideSequenceString = ProteinSequenceUtil.filterSequenceByColumns(peptideSequence, peptideSequenceFilter);
            }
            correlationCounts += learnSequenceToPeptide(domainSequenceString, peptideSequenceString, numDomainPositionsPerFeature, numPeptidePositionsPerFeature);
        }
        return correlationCounts;
    }

    /**
     * Generates the next feature position
     *
     * @param positionArray  An already initialized array of positions, where positionArray.length is equal to the
     *                       number of positions in the feature
     * @param totalPositions The total number of positions possible (i.e. the length of the sequence)
     *                       Note: totalPositions must be larger than the length of the position array
     * @param initialize     True only for the first iteration and must be false afterwards
     * @return The positions for the next feature (an updated version of the input position array)
     */
    private int[] generateFeaturePositions(int[] positionArray, int totalPositions, boolean initialize) {
        //this basically counts upwards in base totalPositions
        if (!initialize) {
            for (int i = positionArray.length - 1; i >= 0; i--) {
                if (positionArray[i] < (totalPositions - (positionArray.length - i))) {
                    //add one to the current position
                    positionArray[i] += 1;
                    //reset all positions to the right
                    for (int j = i + 1; j < positionArray.length; j++) {
                        positionArray[j] = positionArray[j - 1] + 1;
                    }
                    break;
                }
            }
        } else {
            //true for the first call to this method
            for (int i = 0; i < positionArray.length; i++) {
                positionArray[i] = i;
            }
        }
        return positionArray;
    }

    /**
     * Traverse through the sequence and a single binding peptide and add features to the correlation matrix
     *
     * @param domainSequenceString          The sequence of the domain
     * @param peptideSequenceString         The sequence of the binding peptide
     * @param numDomainPositionsPerFeature  The number of domain positions to consider as a feature
     * @param numPeptidePositionsPerFeature The number of peptide positions to consider as a feature
     * @return The number of counts added to the matrix
     */
    private int learnSequenceToPeptide(String domainSequenceString, String peptideSequenceString, int numDomainPositionsPerFeature, int numPeptidePositionsPerFeature) {
        //initialize domain features
        int numberOfDomainFeatures = getCombinations(totalDomainSequenceLength, numDomainPositionsPerFeature);
        int[] domainPositionArray = new int[numDomainPositionsPerFeature];
        boolean initializeDomain = true;
        ResiduePositionPair[] domainFeature = allocateFeature(numDomainPositionsPerFeature);
        //initialize peptide feature
        int numberOfPeptideFeatures = getCombinations(totalPeptideSequenceLength, numPeptidePositionsPerFeature);
        int[] peptidePositionArray = new int[numPeptidePositionsPerFeature];
        boolean initializePeptide = true;
        ResiduePositionPair[] peptideFeature = allocateFeature(numPeptidePositionsPerFeature);

        int correlationCounts = 0;

        //get the ith feature of the domain sequence and the jth feature of the peptide sequence
        for (int i = 0; i < numberOfDomainFeatures; i++) {
            domainPositionArray = generateFeaturePositions(domainPositionArray, domainSequenceString.length(), initializeDomain);
            if (initializeDomain) {
                initializeDomain = false;
            }
            domainFeature = createFeature(domainFeature, domainPositionArray, domainSequenceString);
            if (!featureValid(domainFeature)) {
                continue;
            }
            for (int j = 0; j < numberOfPeptideFeatures; j++) {
                //generate next feature given the previous feature as an array of length equal to how long you want the feature to be
                peptidePositionArray = generateFeaturePositions(peptidePositionArray, peptideSequenceString.length(), initializePeptide);
                if (initializePeptide) {
                    initializePeptide = false;
                }
                peptideFeature = createFeature(peptideFeature, peptidePositionArray, peptideSequenceString);
                if (!featureValid(peptideFeature)) {
                    continue;
                }
                //add counts to the correlation matrix
                addCorrelationCount(domainFeature, peptideFeature);
                correlationCounts += 1;
            }
            initializePeptide = true;
        }
        return correlationCounts;
    }

    /**
     * Checks if the residues stored in the feature are in the aaList - accounts for gap characters or bad letters
     *
     * @param feature The feature to validate
     * @return true if the feature is valid
     */
    private boolean featureValid(ResiduePositionPair[] feature) {
        for (int i = 0; i < feature.length; i++) {
            ResiduePositionPair residuePositionPair = feature[i];
            if (aaList.indexOf(residuePositionPair.residue) < 0) {
                return false;
            }
        }
        return true;
    }

    /**
     * Helper method to allocate memory for a feature
     */
    private ResiduePositionPair[] allocateFeature(int numPositionsPerFeature) {
        ResiduePositionPair[] feature = new ResiduePositionPair[numPositionsPerFeature];
        for (int i = 0; i < feature.length; i++) {
            feature[i] = new ResiduePositionPair(-1, 'X');
        }
        return feature;
    }

    /**
     * Add a single correlation count to the correlation matrix
     *
     * @param domainFeature  The domain feature to address in the matrix
     * @param peptideFeature The peptide feature to address in the matrix
     */
    private void addCorrelationCount(ResiduePositionPair[] domainFeature, ResiduePositionPair[] peptideFeature) {
        int domainIndex = getFeatureIndex(domainFeature);
        int peptideIndex = getFeatureIndex(peptideFeature);
        if (domainPositionsAsRows) {
            correlationMatrix.add(domainIndex, peptideIndex, 1);
        } else {
            correlationMatrix.add(peptideIndex, domainIndex, 1);
        }
        //correlationMatrix.set(domainIndex, peptideIndex, correlationMatrix.get(domainIndex, peptideIndex)+1);
    }

    /**
     * Populates a feature with values
     *
     * @param feature        The feature to populate (previously allocated)
     * @param positionArray  The list of positions to use
     * @param sequenceString The protein sequence
     * @return The populated feature as an array of residue/position pairs
     */
    private ResiduePositionPair[] createFeature(ResiduePositionPair[] feature, int[] positionArray, String sequenceString) {
        for (int i = 0; i < feature.length; i++) {
            ResiduePositionPair residuePositionPair = feature[i];
            residuePositionPair.position = positionArray[i];
            residuePositionPair.residue = sequenceString.charAt(positionArray[i]);
        }
        return feature;
    }

    /**
     * Learn a correlation matrix from the MSA and lists of peptides
     *
     * @param numDomainPositionsPerFeature  The number of positions per feature on the domain sequence
     * @param numPeptidePositionsPerFeature The number of positions per feature on the peptide sequence
     * @return The number of observed correlation counts
     * @throws BioException if there is an error with the sequence alignment object
     */
    public long learn(int numDomainPositionsPerFeature, int numPeptidePositionsPerFeature) throws BioException {
        long numberCorrelationCounts = 0;

        //set the size of the features for later use
        this.numDomainPositionsPerFeature = numDomainPositionsPerFeature;
        this.numPeptidePositionsPerFeature = numPeptidePositionsPerFeature;

        //make sure filters are taken into account for sequence length calculations
        if (domainSequenceFilter != null) {
            totalDomainSequenceLength = ProteinSequenceUtil.countLengthOfFilteredStringResult(domainSequenceFilter);
        }
        if (peptideSequenceFilter != null) {
            totalPeptideSequenceLength = ProteinSequenceUtil.countLengthOfFilteredStringResult(peptideSequenceFilter);
        }

        //initialize the correlation matrix
        int numRows = 0;
        int numColumns = 0;

        int maxDomainIndex = getMaxIndex(numDomainPositionsPerFeature, totalDomainSequenceLength);
        int maxPeptideIndex = getMaxIndex(numPeptidePositionsPerFeature, totalPeptideSequenceLength);
        //we want to keep the number of rows low, because each one is a SparseVector
        //note - for some sparse matrix data structures, like RowCompressed, you may have to switch these based on size
        if (maxDomainIndex > maxPeptideIndex) {
            domainPositionsAsRows = false;
            numRows = maxPeptideIndex;
            numColumns = maxDomainIndex;
        } else {
            domainPositionsAsRows = true;
            numRows = maxDomainIndex;
            numColumns = maxPeptideIndex;
        }
        //create a sparse 2D matrix for storing the correlations.  The matrix is conceptually 2 + the total number of
        //positions considered dimensions, but this is dealt with by computing an index to map from higher D to 2D.
        correlationMatrix = new FlexCompRowMatrix(numRows, numColumns);

        //learn the correlations
        long start = System.currentTimeMillis();
        Collection profiles = profileNameToProfile.values();
        for (Iterator iterator = profiles.iterator(); iterator.hasNext();) {
            ProteinProfile proteinProfile = (ProteinProfile) iterator.next();
            if (sequenceNameToSequence.containsKey(proteinProfile.getName())) {
                //we only care about sequences from the MSA that we have binding peptide data for
                numberCorrelationCounts += learnSequenceToPeptides((Sequence) sequenceNameToSequence.get(proteinProfile.getName()), proteinProfile,
                        numDomainPositionsPerFeature, numPeptidePositionsPerFeature);
            } else {
                System.out.println("No aligned sequence was found for profile: " + proteinProfile.getName());
            }
        }

        //calculate aa frequencies at each position for domain MSA and peptide MSA
        //domain frequency calculation
        calculateDomainMSAFrequencies(numDomainPositionsPerFeature);
        //peptide frequency calculation
        calculatePeptideMSAFrequencies(numPeptidePositionsPerFeature);

        long end = System.currentTimeMillis();
        System.out.println(numberCorrelationCounts + " correlations counted in " + (end - start) / 1000 + " seconds.");

        return numberCorrelationCounts;
    }

    /**
     * Add the given frequency count to the given frequency vector for a sequence feature
     *
     * @param feature         The sequence feature to add a count to
     * @param frequencyVector The frequency vector
     * @param count           The count to add
     */
    private void addFrequencyCount(ResiduePositionPair[] feature, SparseVector frequencyVector, double count) {
        frequencyVector.add(getFeatureIndex(feature), count);
    }

    /**
     * Get the frequency count of the given sequence feature in the given vector
     *
     * @param feature         The sequence feature to address
     * @param frequencyVector The frequency vector that stores the frequency information
     * @return The frequency count
     */
    private double getFrequencyCount(ResiduePositionPair[] feature, SparseVector frequencyVector) {
        return frequencyVector.get(getFeatureIndex(feature));
    }

    /**
     * Go through all domain sequences and calculate the feature frequencies as if the domain sequence was
     * in a large multiple sequence alignment where each line was a unique domain-peptide pair
     *
     * @param numDomainPositionsPerFeature The number of domain sequence positions per feature (residues per feature)
     */
    private void calculateDomainMSAFrequencies(int numDomainPositionsPerFeature) {
        //initialize domain feature and frequency vector data
        int numberOfDomainFeatures = getCombinations(totalDomainSequenceLength, numDomainPositionsPerFeature);
        int[] domainPositionArray = new int[numDomainPositionsPerFeature];
        boolean initializeDomain = true;
        ResiduePositionPair[] domainFeature = allocateFeature(numDomainPositionsPerFeature);
        domainFeatureFrequencyVector = new SparseVector(getMaxIndex(numDomainPositionsPerFeature, totalDomainSequenceLength));

        //calculate the domain feature frequency
        Collection alignedSequences = sequenceNameToSequence.values();
        for (Iterator iterator = alignedSequences.iterator(); iterator.hasNext();) {
            Sequence domainSequence = (Sequence) iterator.next();
            if (profileNameToProfile.containsKey(domainSequence.getName())) {
                //only calculate frequency counts for this sequence if we also have the corresponding profile
                String domainSequenceString = null;
                //take into account sequence filtering, if necessary
                if (domainSequenceFilter == null) {
                    domainSequenceString = domainSequence.seqString();
                } else {
                    domainSequenceString = ProteinSequenceUtil.filterSequenceByColumns(domainSequence, domainSequenceFilter);
                }
                initializeDomain = true;
                ProteinProfile domainProfile = (ProteinProfile) profileNameToProfile.get(domainSequence.getName());
                for (int i = 0; i < numberOfDomainFeatures; i++) {
                    domainPositionArray = generateFeaturePositions(domainPositionArray, domainSequenceString.length(), initializeDomain);
                    if (initializeDomain) {
                        initializeDomain = false;
                    }
                    domainFeature = createFeature(domainFeature, domainPositionArray, domainSequenceString);
                    if (!featureValid(domainFeature)) {
                        continue;
                    }
                    //multiply the domain frequency count by the number of peptides in a profile, since we are considering
                    //a virtual multiple sequence aligment composed of unique domain-peptide pairs
                    addFrequencyCount(domainFeature, domainFeatureFrequencyVector, domainProfile.getNumSequences());
                }
            }
        }
    }

    /**
     * Go through all peptide sequences and calculate the feature frequencies as if the peptide sequence was
     * in a large multiple sequence alignment where each line was a unique domain-peptide pair
     *
     * @param numPeptidePositionsPerFeature The number of peptide sequence positions per feature (residues per feature)
     */
    private void calculatePeptideMSAFrequencies(int numPeptidePositionsPerFeature) {
        //initialize peptide feature and frequency vector data
        int numberOfPeptideFeatures = getCombinations(totalPeptideSequenceLength, numPeptidePositionsPerFeature);
        int[] peptidePositionArray = new int[numPeptidePositionsPerFeature];
        boolean initializePeptide = true;
        ResiduePositionPair[] peptideFeature = allocateFeature(numPeptidePositionsPerFeature);
        peptideFeatureFrequencyVector = new SparseVector(getMaxIndex(numPeptidePositionsPerFeature, totalPeptideSequenceLength));

        Collection proteinProfileList = profileNameToProfile.values();
        for (Iterator iterator = proteinProfileList.iterator(); iterator.hasNext();) {
            ProteinProfile proteinProfile = (ProteinProfile) iterator.next();
            Collection peptides = proteinProfile.getSequenceMap();
            for (Iterator peptideIterator = peptides.iterator(); peptideIterator.hasNext();) {
                Sequence peptideSequence = (Sequence) peptideIterator.next();
                String peptideSequenceString = null;
                //take into account sequence filtering, if necessary
                if (peptideSequenceFilter == null) {
                    peptideSequenceString = peptideSequence.seqString();
                } else {
                    peptideSequenceString = ProteinSequenceUtil.filterSequenceByColumns(peptideSequence, peptideSequenceFilter);
                }
                for (int j = 0; j < numberOfPeptideFeatures; j++) {
                    peptidePositionArray = generateFeaturePositions(peptidePositionArray, peptideSequenceString.length(), initializePeptide);
                    if (initializePeptide) {
                        initializePeptide = false;
                    }
                    peptideFeature = createFeature(peptideFeature, peptidePositionArray, peptideSequenceString);
                    if (!featureValid(peptideFeature)) {
                        continue;
                    }
                    addFrequencyCount(peptideFeature, peptideFeatureFrequencyVector, 1);
                }
                initializePeptide = true;
            }
        }
    }

    /**
     * Helper class to store a pair of integer hashes corresponding to a domain and a peptide feature
     */
    private class CorrelationResult {
        public int domainFeature;
        public int peptideFeature;
        public int correlationCount;

        public CorrelationResult(int domainFeature, int peptideFeature, int correlationCount) {
            this.domainFeature = domainFeature;
            this.peptideFeature = peptideFeature;
            this.correlationCount = correlationCount;
        }
    }

    /**
     * Fetches the most informative residue feature correlations (above a threshold)
     *
     * @param scoreThreshold The threshold to filter
     * @return A sorted map with the score as a key and an ArrayList of CorrelationResult objects as a value
     */
    private TreeMap getMostInformativeFeatures(double scoreThreshold) {
        //calculate the conditional entropy for all counted correlations
        //initialization
        TreeMap sortedResultMap = new TreeMap(); //stores a list of all features, sorted by score
        Iterator iterator = correlationMatrix.iterator();
        double domainFrequency = 0.0;
        double peptideFrequency = 0.0;
        double score = 0.0;
        ResiduePositionPair[] domainFeature = allocateFeature(numDomainPositionsPerFeature);
        ResiduePositionPair[] peptideFeature = allocateFeature(numPeptidePositionsPerFeature);
        //score each counted feature
        while (iterator.hasNext()) {
            MatrixEntry matrixEntry = (MatrixEntry) iterator.next();
            domainFeature = getDomainFeatureFromSparseMatrixEntry(matrixEntry, domainFeature);
            peptideFeature = getPeptideFeatureFromSparseMatrixEntry(matrixEntry, peptideFeature);
            domainFrequency = getFrequencyCount(domainFeature, domainFeatureFrequencyVector);
            peptideFrequency = getFrequencyCount(peptideFeature, peptideFeatureFrequencyVector);
            score = scoreFeature(matrixEntry, domainFrequency, peptideFrequency);
            if (score < scoreThreshold) {   //lower score is better - TODO: generalize
                addResultToSortedResultMap(score, matrixEntry, sortedResultMap);
                //System.out.print(featureToString(domainFeature) + "\t" + featureToString(peptideFeature) + "\t" + matrixEntry.get() + "\t");
                //System.out.println(domainFrequency + "\t" + peptideFrequency + "\t" + score);
            }
        }
        return sortedResultMap;
    }

    /**
     * Prints the most informative residue feature correlations (above a threshold) to the screen
     *
     * @param chosenSequenceName The sequence to put the domain residue feature in context of
     * @param scoreThreshold     The score threshold to filter at
     */
    public void printMostInformativeFeatures(String chosenSequenceName, double scoreThreshold) {
        ResiduePositionPair[] domainFeature = allocateFeature(numDomainPositionsPerFeature);
        ResiduePositionPair[] peptideFeature = allocateFeature(numPeptidePositionsPerFeature);
        TreeMap sortedResultMap = getMostInformativeFeatures(scoreThreshold);
        System.out.println("Features better than score " + scoreThreshold + " (lower is better) shown with " + chosenSequenceName);
        System.out.println("Total domain length: " + totalDomainSequenceLength);
        System.out.println("Total peptide length: " + totalPeptideSequenceLength);
        Set scores = sortedResultMap.keySet();
        for (Iterator iterator = scores.iterator(); iterator.hasNext();) {
            Double score = (Double) iterator.next();
            ArrayList resultList = (ArrayList) sortedResultMap.get(score);
            for (int i = 0; i < resultList.size(); i++) {
                CorrelationResult correlationResult = (CorrelationResult) resultList.get(i);
                printResult(correlationResult, chosenSequenceName, domainFeature, peptideFeature, score);
            }
        }
    }

    /**
     * Prints a single result
     *
     * @param correlationResult  The result to print
     * @param chosenSequenceName The sequence name to print the result in context of
     * @param domainFeature      The hashed domain feature
     * @param peptideFeature     The hashed peptide feature
     * @param score              The score of this result
     */
    private void printResult(CorrelationResult correlationResult, String chosenSequenceName, ResiduePositionPair[] domainFeature, ResiduePositionPair[] peptideFeature, Double score) {
        domainFeature = indexToFeature(correlationResult.domainFeature, domainFeature, totalDomainSequenceLength);
        peptideFeature = indexToFeature(correlationResult.peptideFeature, peptideFeature, totalPeptideSequenceLength);
        //print feature in short form
        int domainFrequency = (int) getFrequencyCount(domainFeature, domainFeatureFrequencyVector);
        int peptideFrequency = (int) getFrequencyCount(peptideFeature, peptideFeatureFrequencyVector);
        System.out.println(featureToString(domainFeature) + " " + featureToString(peptideFeature) + " "
                + score + " (" + correlationResult.correlationCount + "," + domainFrequency + "," + peptideFrequency + ")");
        //print domain sequence with feature highlighted
        Sequence chosenSequence = (Sequence) sequenceNameToSequence.get(chosenSequenceName);
        String sequenceString = chosenSequence.seqString();
        StringBuffer sb = new StringBuffer(sequenceString);
        for (int i = domainFeature.length - 1; i >= 0; i--) {
            //iterate in reverse order - assumes the feature will be sorted by position
            ResiduePositionPair residuePositionPair = domainFeature[i];
            sb.replace(residuePositionPair.position, residuePositionPair.position + 1, "[" + sequenceString.charAt(residuePositionPair.position) + "]");
        }
        System.out.println(sb);
        //print peptide sequence with feature highlighted
        boolean printedAResidue = false;
        for (int i = 0; i < totalPeptideSequenceLength; i++) {
            //iterate through possible peptide positions
            for (int j = 0; j < peptideFeature.length; j++) {
                //search our feature (array of residue-position pairs)
                ResiduePositionPair residuePositionPair = peptideFeature[j];
                if (residuePositionPair.position == i) {
                    System.out.print(residuePositionPair.residue);
                    printedAResidue = true;
                    break;
                }
            }
            if (!printedAResidue) {
                //if no match in the search
                System.out.print("-");
            } else {
                //reset the search
                printedAResidue = false;
            }
        }
        System.out.print("\n");
    }

    /**
     * Helper method to add a score and a matrixEntry result to the sorted result map
     */
    private void addResultToSortedResultMap(double score, MatrixEntry matrixEntry, TreeMap sortedResultMap) {
        //check if we have already seen this score in the score map
        Double scoreDouble = new Double(score);
        ArrayList resultList = null;
        //many features could have the same score, so maintain a list
        if (!sortedResultMap.containsKey(scoreDouble)) {
            //new result
            resultList = new ArrayList();
        } else {
            //existing result
            resultList = (ArrayList) sortedResultMap.get(scoreDouble);
        }
        //don't actually add the matrixEntry, since it is not allocated on each access (singleton)
        CorrelationResult result = null;
        if (domainPositionsAsRows) {
            result = new CorrelationResult(matrixEntry.row(), matrixEntry.column(), (int) matrixEntry.get());
        } else {
            result = new CorrelationResult(matrixEntry.column(), matrixEntry.row(), (int) matrixEntry.get());
        }
        resultList.add(result);
        sortedResultMap.put(scoreDouble, resultList);
    }

    /**
     * Helper method to access a domain feature from a sparse matrix entry
     */
    private ResiduePositionPair[] getDomainFeatureFromSparseMatrixEntry(MatrixEntry matrixEntry, ResiduePositionPair[] preAllocatedDomainFeature) {
        int columnIndex = 0;
        int rowIndex = 0;
        ResiduePositionPair[] domainFeature = null;
        columnIndex = matrixEntry.column();
        rowIndex = matrixEntry.row();
        if (domainPositionsAsRows) {
            domainFeature = indexToFeature(rowIndex, preAllocatedDomainFeature, totalDomainSequenceLength);
        } else {
            domainFeature = indexToFeature(columnIndex, preAllocatedDomainFeature, totalDomainSequenceLength);
        }
        return domainFeature;
    }

    /**
     * Helper method to access a peptide feature from a sparse matrix entry
     */
    private ResiduePositionPair[] getPeptideFeatureFromSparseMatrixEntry(MatrixEntry matrixEntry, ResiduePositionPair[] preAllocatedPeptideFeature) {
        int columnIndex = 0;
        int rowIndex = 0;
        ResiduePositionPair[] peptideFeature = null;
        columnIndex = matrixEntry.column();
        rowIndex = matrixEntry.row();
        if (domainPositionsAsRows) {
            peptideFeature = indexToFeature(columnIndex, preAllocatedPeptideFeature, totalPeptideSequenceLength);
        } else {
            peptideFeature = indexToFeature(rowIndex, preAllocatedPeptideFeature, totalPeptideSequenceLength);
        }
        return peptideFeature;
    }

    /**
     * Returns the score of the feature
     *
     * @param matrixEntry      contains the feature (encoded) and the feature frequency
     * @param domainFrequency  Contains the frequency of the part of the feature in the domain alignment
     * @param peptideFrequency Contains the frequency of the part of the feature in the peptide alignment
     * @return The score
     */
    private double scoreFeature(MatrixEntry matrixEntry, double domainFrequency, double peptideFrequency) {
        //currently the score is the conditional entropy
        double conditionalEntropy = (-Math.log((matrixEntry.get() * multipleSequenceAlignmentLength) / (domainFrequency * peptideFrequency)) *
                (matrixEntry.get() / multipleSequenceAlignmentLength));
        return conditionalEntropy;
    }

    /**
     * Given a sequence from the original alignment that was used in the learning step, predict
     * a protein profile based on all correlations learned
     *
     * @param alignedDomainSequence
     * @return
     */
    public ProteinProfile predictProfileFromSequence(Sequence alignedDomainSequence) {
        ProteinProfile proteinProfile = null;
        String seqString = alignedDomainSequence.seqString();
        double domainFrequency = 0.0;
        double peptideFrequency = 0.0;
        double conditionalEntropy = 0.0;
        Distribution currentDist = null;
        HashMap alphabetMap = ProteinSequenceUtil.get20aaAlphabet();  //an easy lookup table for Symbols from the protein alphabet
        Symbol residue = null;

        //create a set of empty distributions
        Distribution[] weightMatrixColumns = null;
        SimpleWeightMatrix weightMatrix = null;
        weightMatrixColumns = new Distribution[totalPeptideSequenceLength];
        IndexedCount c = new IndexedCount(ProteinTools.getAlphabet());
        for (int i = 0; i < weightMatrixColumns.length; i++) {
            weightMatrixColumns[i] = DistributionTools.countToDistribution(c);
            //set all weights to zero
            FiniteAlphabet fa = ProteinTools.getAlphabet();
            Iterator symbols = fa.iterator();
            while (symbols.hasNext()) {
                Symbol symbol = (Symbol) symbols.next();
                try {
                    weightMatrixColumns[i].setWeight(symbol, 1E-10); //set a low weight to avoid zeros (bad for drawing logos)
                } catch (IllegalSymbolException e) {
                    e.printStackTrace();  //We're using the standard alphabet
                } catch (ChangeVetoException e) {
                    e.printStackTrace();  //Not setting any change vetos
                }
            }
        }

        //iterate through the correlation matrix and filter based on the given domain sequence
        //(Note: this is more efficient because the matrix is expected to be much more sparse than
        //the possible number of features in a peptide)
        Iterator iterator = correlationMatrix.iterator();
        ResiduePositionPair[] domainFeature = allocateFeature(numDomainPositionsPerFeature);
        ResiduePositionPair[] peptideFeature = allocateFeature(numPeptidePositionsPerFeature);
        while (iterator.hasNext()) {
            MatrixEntry matrixEntry = (MatrixEntry) iterator.next();
            domainFeature = getDomainFeatureFromSparseMatrixEntry(matrixEntry, domainFeature);
            peptideFeature = getPeptideFeatureFromSparseMatrixEntry(matrixEntry, peptideFeature);
            //check if the feature is in the sequence (filter step)
            if (isFeatureInSequence(domainFeature, seqString)) {
                //now set the expected peptide frequency distribution
                domainFrequency = getFrequencyCount(domainFeature, domainFeatureFrequencyVector);
                peptideFrequency = getFrequencyCount(peptideFeature, peptideFeatureFrequencyVector);
                conditionalEntropy = scoreFeature(matrixEntry, domainFrequency, peptideFrequency);
                if (conditionalEntropy < -0.03) {
                    for (int i = 0; i < peptideFeature.length; i++) {
                        ResiduePositionPair residuePositionPair = peptideFeature[i];
                        //find the symbol for the residue at this position
                        residue = (Symbol) alphabetMap.get(String.valueOf(residuePositionPair.residue));
                        //get the distribution for the column at this position
                        currentDist = weightMatrixColumns[residuePositionPair.position];
                        //set the new weight
                        try {
                            currentDist.setWeight(residue, currentDist.getWeight(residue) + Math.abs(conditionalEntropy));
                        } catch (IllegalSymbolException e) {
                            e.printStackTrace(); //should never happen, since symbol is from the standard alphabet
                        } catch (ChangeVetoException e) {
                            e.printStackTrace(); //we are never setting a change veto
                        }
                    }
                    //System.out.print(featureToString(domainFeature) + "\t" + featureToString(peptideFeature) + "\t" + matrixEntry.get() + "\t");
                    //System.out.println(domainFrequency + "\t" + peptideFrequency + "\t" + conditionalEntropy);
                }
                //TODO: count up all peptide frequencies and divide by the max possible number, which
                //TODO: is dependent on the number of peptide features
            }
        }

        //re-weight the columns based on entropy counted per column
        for (int i = 0; i < weightMatrixColumns.length; i++) {
            Distribution weightMatrixColumn = weightMatrixColumns[i];
        }

        //create a weighmatrix and use that to create a proteinprofile
        try {
            weightMatrix = new SimpleWeightMatrix(weightMatrixColumns);
            proteinProfile = new ProteinProfile(weightMatrix, alignedDomainSequence.getName());
        } catch (BioException e) {
            e.printStackTrace(); //Should never happen, since we're using the standard protein alphabet
        }
        //proteinProfile.normalize();
        //TODO: scale down so that column with max sum of information is no more that 4.32 bits (for logo output)

        return proteinProfile;
    }

    public void outputAllLogos() {
        ProteinProfile proteinProfile = null;
        String outFileName = null;
        String outputDirectory = "D:\\Gbader\\Code\\PDZ\\data\\PDZ\\SpecificityPrediction\\Logos";
        Collection alignedSequences = sequenceNameToSequence.values();
        for (Iterator iterator = alignedSequences.iterator(); iterator.hasNext();) {
            Sequence sequence = (Sequence) iterator.next();
            //if(!sequence.getName().equals("MPDZ-2")) {continue;}
            proteinProfile = predictProfileFromSequence(sequence);
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

    /**
     * Checks if a feature is part of a given sequence
     *
     * @return true if the feature is in the sequence
     */
    private boolean isFeatureInSequence(ResiduePositionPair[] sequenceFeature, String seqString) {
        for (int i = 0; i < sequenceFeature.length; i++) {
            ResiduePositionPair residuePositionPair = sequenceFeature[i];
            if (seqString.charAt(residuePositionPair.position) != residuePositionPair.residue) {
                return false;
            }
        }
        return true;
    }

    /**
     * Converts a feature to a simple string representation
     */
    private String featureToString(ResiduePositionPair[] feature) {
        StringBuffer sb = new StringBuffer();
        for (int i = 0; i < feature.length; i++) {
            sb.append(feature[i].toString());
        }
        return sb.toString();
    }

    /**
     * Set the domain sequence filter string
     *
     * @param domainSequenceFilter A filter for domain positions e.g. "1-4,5,8".
     *                             May be set to null to turn off the filter.
     */
    public void setDomainSequenceFilter(String domainSequenceFilter) {
        this.domainSequenceFilter = domainSequenceFilter;
    }

    /**
     * Set the peptide sequence filter string
     *
     * @param peptideSequenceFilter A filter for peptide positions e.g. "1-4,5,8"
     *                              May be set to null to turn off the filter.
     */
    public void setPeptideSequenceFilter(String peptideSequenceFilter) {
        this.peptideSequenceFilter = peptideSequenceFilter;
    }
}
