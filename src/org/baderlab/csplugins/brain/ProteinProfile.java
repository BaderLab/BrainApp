package org.baderlab.csplugins.brain;

import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.DistributionTools;
import org.biojava.bio.dp.SimpleWeightMatrix;
import org.biojava.bio.dp.WeightMatrix;
import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.*;
import org.biojava.utils.ChangeVetoException;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

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
 * * Time: 4:05:58 PM
 * * Description A protein profile optimized for motif (short patterns in a protein) searching
 */

/**
 * A protein profile optimized for motif (short patterns in a protein) searching
 */
public class ProteinProfile {
    private SimpleWeightMatrix weightMatrix = null; //the actual weight matrix
    private double fuzzFactor = 0.0;  //the pseudocount weight for the profile
    private String name = null;  //the name of the profile
    private int numSequences = 0; //the number of sequences used to build the profile
    private double profileMinPValue = 0.0; //the minumum possible p-value for this profile (used in normalization)
    private double profileMaxPValue = 0.0; //the maximum possible p-value for this profile (used in normalization)
    private String experimentalMethod = null;   //The experimental method used to derive this profile
    //The following fields describe the domain that binds to the binding profile represented by this profile
    private DatabaseReference proteinReference = null;  //The database reference of the protein containing the domain
    private int domainNumber = 0;   //The domain number in the protein (counted starting at 1 from N to C)
    private DatabaseReference domainReference = null;    //The database reference of this domain
    private String domainName = null; //The domain name corresponding to domainReference
    private String domainSequence = null;   //The amino acid sequence of this domain
    private int domainSequenceStart = 0;    //The starting position of the domain sequence in the domain containing protein
    private int domainSequenceStop = 0; // The ending position of the domain sequence in the domain containing protein
    private String comment = null;  // The user-provided text comment on this profile
    private String proteinName = null; // the name of the protein

    private Map sequenceMap = null; //a map of original sequences in the profile (maps sequence name to sequence)

    private ArrayList alphabet = null;  //save the 20 aa alphabet for fast use

    /**
     * Create a new protein motif profile with fuzz factor=0.0
     *
     * @param trainingSequences A set of uniform length protein sequences that form the profile
     * @throws BioException if there is an error in the list of protein training sequences
     */
    public ProteinProfile(SequenceIterator trainingSequences) throws BioException {
        this(trainingSequences, 0.0, null);
    }

    /**
     * Create a new protein motif profile with a given fuzz factor
     *
     * @param trainingSequences A set of uniform length protein sequences that form the profile
     * @param fuzzFactor        A number to add to every entry in the protein profile matrix. The higher the number,
     *                          the fuzzier the profile will be (like adding layers of sand on a chess board, the more sand layers added,
     *                          the less you will be able to make out the shapes of the chess pieces). Fuzzier profiles will result in more
     *                          matches. This is sometimes called the pseudo-count number.
     *                          Note: setting a fuzz factor to 0.0 is not recommended if there are lots of zeros in
     *                          the protein profile, since most resulting p-values will then be zero.
     * @param name              The name of the profile
     * @throws BioException if there is an error in the list of protein training sequences
     */
    public ProteinProfile(SequenceIterator trainingSequences, double fuzzFactor, String name) throws BioException {
        this.fuzzFactor = fuzzFactor;
        alphabet = get20aaAlphabet();
        buildProfile(trainingSequences);
        this.name = name;
    }

    /**
     * Create a new protein motif profile from a list existing profiles (a combination of sequences from all profiles)
     * The combination is done by averaging over the distributions to avoid sequence number bias, but the numbers
     * of sequences used to derive the profiles are maintained (added up over all profiles)
     * Note: the sequences from the original profiles are not maintained
     *
     * @param profileList A List of ProteinProfile objects to combine
     * @param name        The name of the profile
     */
    public ProteinProfile(List profileList, String name) {
        this.fuzzFactor = 0.0;
        alphabet = get20aaAlphabet();
        this.name = name;
        numSequences = 0;
        //create a new weightmatrix so we don't overwrite the one passed
        try {
            //recreate profile #1 from its sequences (easy way to copy the weightmatrix)
            ProteinProfile proteinProfile = (ProteinProfile) profileList.get(0);
            Alignment align = new SimpleAlignment(proteinProfile.sequenceMap);
            Distribution[] dists = DistributionTools.distOverAlignment(align, false, fuzzFactor);
            weightMatrix = new SimpleWeightMatrix(dists);
            this.normalize();
            numSequences += proteinProfile.getNumSequences();
        } catch (IllegalAlphabetException e) {
            e.printStackTrace();
        }
        //start from the 2nd profile
        for (int i = 1; i < profileList.size(); i++) {
            ProteinProfile proteinProfile = (ProteinProfile) profileList.get(i);
            //average the weightmatrices
            WeightMatrix wm2 = proteinProfile.getWeightMatrix();
            double weight = 0.0;
            for (int j = 0; j < weightMatrix.columns(); j++) {
                Distribution column1 = weightMatrix.getColumn(j);
                Distribution column2 = wm2.getColumn(j);
                Iterator symbolIterator = alphabet.iterator();
                while (symbolIterator.hasNext()) {
                    Symbol symbol = (Symbol) symbolIterator.next();
                    try {
                        weight = column1.getWeight(symbol) + column2.getWeight(symbol);
                        //if we are at the end of the profile list, take the average
                        if (i != (profileList.size() - 1)) {
                            //not at end yet
                            column1.setWeight(symbol, weight); //save the sum in weightMatrix
                        } else {
                            //we're at the end of the profileList
                            column1.setWeight(symbol, weight / profileList.size()); //save the average
                        }
                    } catch (IllegalSymbolException e) {
                        //this should never happen, since we're using the standard protein alphabet
                        e.printStackTrace();
                    } catch (ChangeVetoException e) {
                        e.printStackTrace();
                    }
                }
            }
            numSequences += proteinProfile.getNumSequences();
        }
        calculateMinMaxPValues();
    }

    /**
     * Create a new protein profile based on an existing weight matrix
     *
     * @param weightMatrix The existing weight matrix
     * @param name         The profile name
     */
    public ProteinProfile(SimpleWeightMatrix weightMatrix, String name) {
        this.weightMatrix = weightMatrix;
        this.name = name;
        this.fuzzFactor = 0.0;
        alphabet = get20aaAlphabet();
        numSequences = 0;
    }

    /**
     * Saves a 20aa biojava alphabet for ease of use (ignore the U symbol)
     */
    private ArrayList get20aaAlphabet() {
        ArrayList alphabet = new ArrayList(20);
        HashMap alphabetMap = ProteinSequenceUtil.get20aaAlphabet();
        Collection symbols = alphabetMap.values();
        for (Iterator iterator = symbols.iterator(); iterator.hasNext();) {
            Symbol symbol = (Symbol) iterator.next();
            alphabet.add(symbol);
        }
        return alphabet;
    }

    /**
     * Internal method used by the constructors to build the profile from the training sequences
     */
    private void buildProfile(SequenceIterator trainingSequences) throws BioException {
        //Create an alignment using the training sequences
        numSequences = 0;
        sequenceMap = new HashMap();
        while (trainingSequences.hasNext()) {
            Sequence sequence = null;
            sequence = (Sequence) trainingSequences.nextSequence();
            sequenceMap.put(sequence.getName(), sequence);
            numSequences++;
        }
        Alignment align = new SimpleAlignment(sequenceMap);
        Distribution[] dists = DistributionTools.distOverAlignment(align, false, fuzzFactor);
        weightMatrix = new SimpleWeightMatrix(dists);
        //re-normalize
        this.normalize();
        calculateMinMaxPValues();
    }

    private void calculateMinMaxPValues() {
        if (alphabet == null) {
            throw new IllegalStateException("Alphabet not set in ProteinProfile.");
        }
        double weight = 0.0;
        profileMinPValue = 1.0;
        profileMaxPValue = 1.0;
        //find the min and max possible p-values for this profile
        for (int i = 0; i < weightMatrix.columns(); i++) {
            Distribution column = weightMatrix.getColumn(i);
            Iterator symbolIterator = alphabet.iterator();
            //find the min and max p-value for this column
            double columnMinPValue = Double.MAX_VALUE;
            double columnMaxPValue = Double.MIN_VALUE;
            while (symbolIterator.hasNext()) {
                Symbol symbol = (Symbol) symbolIterator.next();
                try {
                    weight = column.getWeight(symbol);
                } catch (IllegalSymbolException e) {
                    //this should never happen, since we're using the standard protein alphabet
                    e.printStackTrace();
                }
                columnMinPValue = Math.min(columnMinPValue, weight);
                columnMaxPValue = Math.max(columnMaxPValue, weight);
            }
            profileMinPValue = profileMinPValue * columnMinPValue;
            profileMaxPValue = profileMaxPValue * columnMaxPValue;
        }
    }

    /**
     * Get the BioJava weigth matrix describing the profile
     *
     * @return
     */
    public SimpleWeightMatrix getWeightMatrix() {
        return weightMatrix;
    }

    /**
     * Gets the name of the profile, if set.
     */
    public String getName() {
        return name;
    }

    /**
     * Gets the number of columns (length) of the profile
     */
    public int getNumColumns() {
        return weightMatrix.columns();
    }

    /**
     * Gets the number of sequences used to build the profile
     */
    public int getNumSequences() {
        return numSequences;
    }

    /**
     * Return a normalized p-value given a p-value that was returned from a search
     * A profile has a max and min possible p-value for the most optimal and non-optimal matching
     * sequence.  This method normalizes the search p-value to the interval 0..1
     * This method doesn't check the input p-value for valid input (within the min..max range of
     * the profile) and will return values outside of the 0..1 interval for these inputs. It is up to
     * the caller to check for this.
     *
     * @param pvalue The search p-value to normalize
     * @return The normalized p-value
     */
    public double getNormalizedPValue(double pvalue) {
        //normalize the p-value
        double normalizedPValue = (pvalue - profileMinPValue) / (profileMaxPValue - profileMinPValue);
        return (normalizedPValue);
    }

    /**
     * Normalize each column in the profile to sum to 1.0
     */
    public void normalize() {
        if (alphabet == null) {
            throw new IllegalStateException("Alphabet not set in ProteinProfile.");
        }
        WeightMatrix wm = this.getWeightMatrix();

        try {
            for (int i = 0; i < wm.columns(); i++) {
                //each column is normalized to sum to 1.0
                Distribution d = wm.getColumn(i);
                //get sum of column weights
                double weightSum = 0.0;
                Iterator l = alphabet.iterator();
                while (l.hasNext()) {
                    Symbol symbol = (Symbol) l.next();
                    weightSum += d.getWeight(symbol);
                }
                //now go back and normalize column to 1.0
                l = alphabet.iterator();
                while (l.hasNext()) {
                    Symbol symbol = (Symbol) l.next();
                    d.setWeight(symbol, d.getWeight(symbol) / weightSum);
                }
            }
        } catch (BioException e) {
            e.printStackTrace();
        } catch (ChangeVetoException e) {
            e.printStackTrace();
        }
    }

    /**
     * Re-weight the profile to correct for an NNK codon bias (relevant to e.g. biased phage display libraries)
     */
    public void reWeightByNNKCodonBias() {
        if (alphabet == null) {
            throw new IllegalStateException("Alphabet not set in ProteinProfile.");
        }
        final String aaList = "ACDEFGHIKLMNPQRSTVWY";
        int[] codonBias = {2, 1, 1, 1, 1, 2, 1, 1, 1, 3, 1, 1, 2, 1, 3, 3, 2, 2, 1, 1};

        WeightMatrix wm = this.getWeightMatrix();

        try {
            for (int i = 0; i < wm.columns(); i++) {
                Distribution d = wm.getColumn(i);
                double weight = 0.0;
                //divide all weights by the number of codons, keeping track of the sum of the new weights
                Iterator l = alphabet.iterator();
                while (l.hasNext()) {
                    Symbol symbol = (Symbol) l.next();
                    String symbolString = ProteinTools.getAlphabet().getTokenization("token").tokenizeSymbol(symbol);
                    if (!symbolString.equals("U")) {
                        weight = d.getWeight(symbol) / codonBias[aaList.indexOf(symbolString)];
                        d.setWeight(symbol, weight);
                    }
                }
            }
        } catch (BioException e) {
            e.printStackTrace();
        } catch (ChangeVetoException e) {
            e.printStackTrace();
        }
        //now go back and normalize columns to 1.0
        this.normalize();
        this.calculateMinMaxPValues();
    }

    /**
     * Re-weight the profile to correct for a codon bias defined in a file (relevant to e.g. biased phage display libraries)
     *
     * @param codonBiasFile A file describing the codon bias in the following format:
     *                      AminoAcidResidueOneLetterCode\tpercentageOfLibrary
     *                      Where a stop codon is represented by the letter X
     *                      For example:
     *                      A	6.25
     *                      C	3.125
     *                      D	3.125
     *                      ...
     *                      X	3.125
     *                      <p/>
     *                      Note: it doesn't matter what the percentage column sums to, but the codons must maintain their relative weight
     *                      Also, the stop codon is optional, since it is not used in the protein profile and it doesn't matter if we ignore
     *                      it, since only the relative codon difference matters.
     */
    public void reWeightByCodonBias(File codonBiasFile) {
        if (alphabet == null) {
            throw new IllegalStateException("Alphabet not set in ProteinProfile.");
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

        WeightMatrix wm = this.getWeightMatrix();

        try {
            for (int i = 0; i < wm.columns(); i++) {
                Distribution d = wm.getColumn(i);
                double weight = 0.0;
                //divide all weights by the number of codons, keeping track of the sum of the new weights
                Iterator l = alphabet.iterator();
                while (l.hasNext()) {
                    Symbol symbol = (Symbol) l.next();
                    String symbolString = ProteinTools.getAlphabet().getTokenization("token").tokenizeSymbol(symbol);
                    if (!symbolString.equals("U")) {
                        weight = d.getWeight(symbol) / codonBias[aaList.indexOf(symbolString)];
                        d.setWeight(symbol, weight);
                    }
                }
            }
        } catch (BioException e) {
            e.printStackTrace();
        } catch (ChangeVetoException e) {
            e.printStackTrace();
        }
        //now go back and normalize columns to 1.0
        this.normalize();
        this.calculateMinMaxPValues();
    }

    /**
     * Return a shortened version of the profile.
     *
     * @param profileLength The length of the profile to return. If profileLength is greater than the length of the
     *                      original profile, the original profile will be returned.
     * @param terminus      The terminus to cut from
     * @return a copy of the profile that is shorter or equal to the original profile.
     */
    public ProteinProfile getTruncatedProfileCopy(int profileLength, ProteinTerminus terminus) {
        if (profileLength >= this.getNumColumns()) {
            return this;
        }
        int loopStart = 0;
        int loopEnd = 0;
        if (terminus.equals(ProteinTerminus.C)) {
            loopStart = this.getNumColumns() - profileLength;
            loopEnd = this.getNumColumns();
        } else if (terminus.equals(ProteinTerminus.N)) {
            loopStart = 0;
            loopEnd = profileLength;
        } else {
            throw new RuntimeException("Protein terminus must be either C or N to decide on a profile cut.");
        }
        ProteinProfile cutProfile = getProfileSubsetCopy(loopStart + "-" + loopEnd);

        return cutProfile;
    }

    /**
     * Filter a profile, choosing only specific columns.
     *
     * @param filter A filter for profile columns e.g. "1-4,5,8" will return a profile containing columns 1 to 4
     *               followed by column 5 followed by column 8.
     * @return A profile containing the specified subset of columns
     */
    public ProteinProfile getProfileSubsetCopy(String filter) {
        if (filter == null || filter.equals("")) {
            return null;
        }
        //convert the string filter to a searchable list containing an expanded version of the filter
        ArrayList columnsToKeep = new ArrayList();
        String positions[] = filter.split(",");
        for (int i = 0; i < positions.length; i++) {
            String position = positions[i];
            if (position.indexOf("-") > 0) {
                //interval
                String[] startEnd = position.split("-");
                int start = Integer.parseInt(startEnd[0]);
                int end = Integer.parseInt(startEnd[1]);
                for (int j = start; j < end; j++) {
                    columnsToKeep.add(new Integer(j));
                }
            } else {
                //point
                int point = Integer.parseInt(position);
                columnsToKeep.add(new Integer(point));
            }
        }

        //change the weightmatrix directly
        //create a new weightmatrix so we don't overwrite the one passed
        Distribution[] distArray = new Distribution[columnsToKeep.size()];
        int j = 0;
        //go through all profile positions in order and test if they are to be copied
        for (int i = 0; i < this.getNumColumns(); i++) {
            if (columnsToKeep.contains(new Integer(i))) {
                distArray[j] = this.weightMatrix.getColumn(i);
                j++;
            }
        }

        ProteinProfile cutProfile = null;
        try {
            SimpleWeightMatrix newWM = new SimpleWeightMatrix(distArray);
            cutProfile = new ProteinProfile(newWM, this.getName());
        } catch (IllegalAlphabetException e) {
            e.printStackTrace();
        } catch (BioException e) {
            e.printStackTrace();
        }

        return cutProfile;
    }

    //the following entropy methods modified from DistributionLogo class in BioJava
    private static double bits = Math.log(2.0);

    /**
     * Calculate the information content of a symbol in bits.
     *
     * @param s the symbol to calculate for
     * @throws IllegalSymbolException if s is not within the dist.
     */
    private double entropy(Distribution dist, Symbol s) throws IllegalSymbolException {
        double p = dist.getWeight(s);
        if (p == 0.0) {
            return 0;
        }
        double lp = Math.log(p);

        return -p * lp / bits;
    }

    /**
     * Retrieve the maximal number of bits possible for this type of dist.
     *
     * @return maximum bits as a double
     */
    private double totalBits(Distribution dist) {
        return Math.log(20) / bits;
    }

    /**
     * <p/>
     * Calculates the total information of the dist in bits.
     * </p>
     * <p/>
     * <p/>
     * This calculates <code>totalBits - sum_r(entropy(r))</code>
     * </p>
     *
     * @return the total information in the dist
     */
    private double totalInformation(Distribution dist) {
        double inf = totalBits(dist);

        for (
                Iterator i = ((FiniteAlphabet) dist.getAlphabet()).iterator();
                i.hasNext();
                ) {
            Symbol s = (Symbol) i.next();
            try {
                inf -= entropy(dist, s);
            } catch (IllegalSymbolException ire) {
                throw new BioError(
                        "Symbol evaporated while calculating information", ire);
            }
        }

        return inf;
    }

    /**
     * Calculates the specificity score of this profile.
     * This is the product of the Shannon entropy at all positions. Each position is
     * normalized from 0 to 20, where 20 means most specific and 0 means least specific.
     *
     * @return the specificity score
     */
    public double calculateSpecificityScore() {
        double score = 1.0;
        WeightMatrix wm = this.getWeightMatrix();
        for (int pos = 0; pos < wm.columns(); pos++) {
            Distribution dist = wm.getColumn(pos);
            double informationAtPosition = totalInformation(dist);
            score *= Math.pow(2, informationAtPosition);
        }
        return score;
    }


    public String getExperimentalMethod() {
        return experimentalMethod;
    }

    public void setExperimentalMethod(String experimentalMethod) {
        this.experimentalMethod = experimentalMethod;
    }

    public DatabaseReference getProteinReference() {
        return proteinReference;
    }

    public void setProteinReference(DatabaseReference proteinReference) {
        this.proteinReference = proteinReference;
    }

    public int getDomainNumber() {
        return domainNumber;
    }

    /**
     * Get a collection of sequences used to create this profile
     *
     * @return A collection of BioJava Sequence objects or null if they are not stored in this profile
     */
    public Collection getSequenceMap() {
        if (sequenceMap != null) {
            return sequenceMap.values();
        }
        return null;
    }

    public void setDomainNumber(int domainNumber) {
        this.domainNumber = domainNumber;
    }

    public DatabaseReference getDomainReference() {
        return domainReference;
    }

    public void setDomainReference(DatabaseReference domainReference) {
        this.domainReference = domainReference;
    }

    public String getDomainName() {
        return domainName;
    }

    public void setDomainName(String domainName) {
        this.domainName = domainName;
    }

    public String getDomainSequence() {
        return domainSequence;
    }

    public void setDomainSequence(String domainSequence) {
        this.domainSequence = domainSequence;
    }

    public int getDomainSequenceStart() {
        return domainSequenceStart;
    }

    public void setDomainSequenceStart(int position) {
        this.domainSequenceStart = position;
    }

    public int getDomainSequenceStop() {
        return domainSequenceStop;
    }

    public void setDomainSequenceStop(int position) {
        this.domainSequenceStop = position;
    }

    public String getComment() {
        return comment;
    }

    public void setComment(String text) {
        this.comment = text;
    }

    public String getProteinName() {
        return proteinName;
    }

    public void setProteinName(String text) {
        this.proteinName = text;
    }

    /**
     * Returns a string representation of the protein profile (may be large!)
     */
    public String toString() {
        if (weightMatrix == null) {
            System.err.println("Matrix was null in ProteinProfile object. This should never happen.");
            return null;
        }
        if (alphabet == null) {
            throw new IllegalStateException("Alphabet not set in ProteinProfile.");
        }
        double simpleMatrix[][] = new double[weightMatrix.columns()][alphabet.size()];
        for (int i = 0; i < weightMatrix.columns(); i++) {
            Distribution d = weightMatrix.getColumn(i);
            Iterator it = alphabet.iterator();
            int j = 0;
            while (it.hasNext()) {
                Symbol symbol = (Symbol) it.next();
                try {
                    simpleMatrix[i][j] = d.getWeight(symbol);
                } catch (IllegalSymbolException e) {
                    System.err.println("Illegal symbol found in the weight weightMatrix. This should never happen.");
                }
                j++;
            }
        }
        Iterator it = alphabet.iterator();
        SymbolTokenization st = null;
        try {
            st = ProteinTools.getAlphabet().getTokenization("token");
        } catch (BioException e) {
            System.err.println("Unable to get symboltokenization. This should never happen.");
        }
        StringBuffer sb = new StringBuffer();
        String lineSep = System.getProperty("line.separator");
        for (int i = 0; i < alphabet.size(); i++) {
            String token = null;
            try {
                token = st.tokenizeSymbol((Symbol) it.next());
            } catch (IllegalSymbolException e) {
                System.err.println("Unable to convert symbol to token. This should never happen.");
            }
            sb.append(token + "\t");
            for (int j = 0; j < simpleMatrix.length; j++) {
                sb.append(simpleMatrix[j][i]);
                if (j < (simpleMatrix.length - 1)) {
                    sb.append("\t");
                }
            }
            sb.append(lineSep);
        }
        return sb.toString();
    }

    //TODO: research: fuzz smartly using aa substitution tables
}
