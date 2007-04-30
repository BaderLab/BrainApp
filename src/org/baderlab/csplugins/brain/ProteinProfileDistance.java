package org.baderlab.csplugins.brain;

import org.biojava.bio.BioException;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dp.SimpleWeightMatrix;
import org.biojava.bio.dp.WeightMatrix;
import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;

import java.util.ArrayList;
import java.util.Iterator;

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
 * * Date: Aug 22, 2005
 * * Time: 4:14:52 PM
 */

/**
 * Distance metrics for protein profiles
 */
public class ProteinProfileDistance {
    /**
     * Saves a 20aa biojava alphabet for ease of use (ignore the U symbol)
     *
     * @return a list of amino acid Symbols
     */
    static private ArrayList get20aaAlphabet() {
        ArrayList alphabet = new ArrayList(20);
        FiniteAlphabet symbols = ProteinTools.getAlphabet();
        Iterator symbolIterator = symbols.iterator();
        while (symbolIterator.hasNext()) {
            Symbol symbol = (Symbol) symbolIterator.next();
            String symbolString = null;
            try {
                symbolString = symbols.getTokenization("token").tokenizeSymbol(symbol);
            } catch (BioException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }
            if (!symbolString.equals("U")) {
                alphabet.add(symbol);
            }
        }
        return alphabet;
    }

    /**
     * Calculates the distance between profiles using only the amino acid
     * distribution statistics (not taking into account amino acid similarity)
     * <p/>
     * The following distance metric is used:
     * <p/>
     * D(a,b) = 1/w (sum     (sum      0.5* (Ai,l - Bi,l)^2) )
     * i=1..w  l=1..20
     * <p/>
     * i=number of columns in profile
     * w=length of profile
     * l=amino acid letters in alphabet
     * <p/>
     * Modified slightly from the original paper so that the distance is really fractional.
     * Transcriptional regulatory code of a eukaryotic genome.
     * Nature. 2004 Sep 2;431(7004):99-104.
     *
     * @param profile1 Profile 1 (distance from)
     * @param profile2 Profile 2 (distance to)
     * @return The fractional distance from profile 1 to 2 (0 = identical, 1 = farthest apart)
     */
    static public double calculateDistributionDistance(ProteinProfile profile1, ProteinProfile profile2) {
        if (profile1.getNumColumns() != profile2.getNumColumns()) {
            throw new IllegalArgumentException("Profiles to be compared must be the same length.");
        }

        ArrayList alphabet = get20aaAlphabet();  //save the 20 aa alphabet for fast use

        double distance = 0.0;

        WeightMatrix wm1 = profile1.getWeightMatrix();
        WeightMatrix wm2 = profile2.getWeightMatrix();

        try {
            for (int i = 0; i < wm1.columns(); i++) {
                Distribution d1 = wm1.getColumn(i);
                Distribution d2 = wm2.getColumn(i);
                Iterator l = alphabet.iterator();
                double weightSum = 0.0;
                while (l.hasNext()) {
                    Symbol symbol = (Symbol) l.next();
                    weightSum += Math.pow(d1.getWeight(symbol) - d2.getWeight(symbol), 2);
                }
                weightSum = weightSum / 2.0;
                distance += weightSum;
            }
        } catch (IllegalSymbolException e) {
            e.printStackTrace();
        }

        distance = distance / profile1.getNumColumns();
        return (distance);
    }

    /**
     * Finds the index of the group containing the symbol s
     *
     * @return an index from 0..(the number of groups-1)
     */
    static private int getGroupIndexForSymbol(String[] grouping, Symbol s) {
        int index = 0;
        String symbol = null;
        try {
            symbol = ProteinTools.getAlphabet().getTokenization("token").tokenizeSymbol(s);
        } catch (BioException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
        for (int i = 0; i < grouping.length; i++) {
            String group = grouping[i];
            if (group.indexOf(symbol) >= 0) {
                index = i;
                break;
            }
        }
        return index;
    }

    /**
     * converts a BioJava WeightMatrix to a grouped weight matrix
     *
     * @param grouping Defines the aa grouping
     * @param profile  The weight matrix to convert
     * @param alphabet The alphabet as an ArrayList, for convenience
     * @return The grouped weight matrix of size (number of weight matrix columns) x (number of groups)
     */
    static private double[][] convertWeightMatrixToGroupWeightMatrix(String[] grouping, WeightMatrix profile, ArrayList alphabet) {
        double[][] groupedWeightMatrix = new double[profile.columns()][grouping.length];
        try {
            for (int j = 0; j < profile.columns(); j++) {
                //for each column
                Distribution d1 = profile.getColumn(j);
                Iterator l = alphabet.iterator();
                while (l.hasNext()) {
                    Symbol symbol = (Symbol) l.next();
                    groupedWeightMatrix[j][getGroupIndexForSymbol(grouping, symbol)] += d1.getWeight(symbol);
                }
            }
        } catch (IllegalSymbolException e) {
            e.printStackTrace();
        }
        return groupedWeightMatrix;
    }

    /**
     * Calculates the distance between profiles using only the amino acid
     * distribution statistics grouped by amino acid similarity
     * <p/>
     * The following distance metric is used:
     * <p/>
     * D(a,b) = 1/w (sum     (sum      0.5* (Ai,l - Bi,l)^2) )
     * i=1..w  l=1..number of groups
     * <p/>
     * i=number of columns in profile
     * w=length of profile
     * l=amino acid groups in alphabet
     * <p/>
     * Modified slightly from the original paper so that the distance is really fractional.
     * Transcriptional regulatory code of a eukaryotic genome.
     * Nature. 2004 Sep 2;431(7004):99-104.
     *
     * @param profile1 Profile 1 (distance from)
     * @param profile2 Profile 2 (distance to)
     * @param grouping Defines a grouping of amino acids by physicochemical property. Groups are
     *                 separated by a comma and/or space.
     *                 Examples are:
     *                 "GSTQCYN, KRH, DE, FLAMPWIV"  (the default from weblogo)
     *                 polar amino acids (G,S,T,Y,C,Q,N), basic (K,R,H), acidic (D,E) and hydrophobic (A,V,L,I,P,W,F,M)
     *                 "STQN, KRH, DE, FLAMPWIVCY, G" (modified weblogo to emphasize hydrophobicity and size differences)
     *                 "YFW, KRH, DE, VLIMC, ASTNQ, GP" (from a 1979 paper by Chris Sander)
     *                 polar: "RKDENQ" hydrophobic: "VILMFWC" small "GASP" aliphatic: "ILVAP" basic "KRH" acidic "DE" aromatic: "FWYH"
     * @return The fractional distance from profile 1 to 2 (0 = identical, 1 = farthest apart)
     */
    static public double calculateAAGroupedDistributionDistance(ProteinProfile profile1, ProteinProfile profile2, String grouping) {
        if (profile1.getNumColumns() != profile2.getNumColumns()) {
            throw new IllegalArgumentException("Profiles to be compared must be the same length.");
        }

        ArrayList alphabet = get20aaAlphabet();  //save the 20 aa alphabet for fast use

        //convert grouping to an array
        String[] groupingArray = grouping.split("[, ]");

        //create a 'grouped' weight matrix from each profile
        double[][] groupedWeightMatrix1 = convertWeightMatrixToGroupWeightMatrix(groupingArray, profile1.getWeightMatrix(), alphabet);
        double[][] groupedWeightMatrix2 = convertWeightMatrixToGroupWeightMatrix(groupingArray, profile2.getWeightMatrix(), alphabet);
        //calculate distance between grouped matrices
        //iterate through columns
        double distance = 0.0;
        for (int j = 0; j < groupedWeightMatrix1.length; j++) {
            double[] column1 = groupedWeightMatrix1[j];
            double[] column2 = groupedWeightMatrix2[j];
            double weightSum = 0.0;
            for (int k = 0; k < column1.length; k++) { //length == number of groups in grouping
                weightSum += Math.pow(column1[k] - column2[k], 2);
            }
            weightSum = weightSum / 2.0;
            distance += weightSum;
        }

        distance = distance / profile1.getNumColumns();

        return (distance);
    }

    /**
     * converts a BioJava Distribution to a grouped column
     *
     * @param grouping     Defines the aa grouping
     * @param distribution The distribution of a single column of a weight matrix to convert
     * @param alphabet     The alphabet as an ArrayList, for convenience
     * @return The grouped column of length (number of groups)
     */
    static private double[] convertDistributionToGroupColumn(String[] grouping, Distribution distribution, ArrayList alphabet) {
        double[] groupedColumn = new double[grouping.length];
        try {
            Iterator l = alphabet.iterator();
            while (l.hasNext()) {
                Symbol symbol = (Symbol) l.next();
                groupedColumn[getGroupIndexForSymbol(grouping, symbol)] += distribution.getWeight(symbol);
            }
        } catch (IllegalSymbolException e) {
            e.printStackTrace();
        }
        return groupedColumn;
    }

    /**
     * Calculates the distance between profiles using only the amino acid
     * distribution statistics grouped by amino acid similarity, where there can be different groups at each position
     * <p/>
     * The following distance metric is used:
     * <p/>
     * D(a,b) = 1/w (sum     (sum      0.5* (Ai,l - Bi,l)^2) )
     * i=1..w  l=1..number of groups
     * <p/>
     * i=number of columns in profile
     * w=length of profile
     * l=amino acid groups in alphabet
     * <p/>
     * Modified slightly from the original paper so that the distance is really fractional.
     * Transcriptional regulatory code of a eukaryotic genome.
     * Nature. 2004 Sep 2;431(7004):99-104.
     *
     * @param profile1           Profile 1 (distance from)
     * @param profile2           Profile 2 (distance to)
     * @param groupingByPosition Defines a grouping of amino acids by physicochemical property. Groups are
     *                           separated by a comma and/or space. A grouping for each position in the profile must be
     *                           specified as elements of a String array, where the position in the array corresponds to the
     *                           position in the profile (numbered from N=0 to C=lengthOfProtein-1).
     *                           Examples are:
     *                           "GSTQCYN, KRH, DE, FLAMPWIV"  (the default from weblogo)
     *                           polar amino acids (G,S,T,Y,C,Q,N), basic (K,R,H), acidic (D,E) and hydrophobic (A,V,L,I,P,W,F,M)
     *                           "STQN, KRH, DE, FLAMPWIVCY, G" (modified weblogo to emphasize hydrophobicity and size differences)
     *                           "YFW, KRH, DE, VLIMC, ASTNQ, GP" (from a 1979 paper by Chris Sander)
     *                           polar: "RKDENQ" hydrophobic: "VILMFWC" small "GASP" aliphatic: "ILVAP" basic "KRH" acidic "DE" aromatic: "FWYH"
     * @return The fractional distance from profile 1 to 2 (0 = identical, 1 = farthest apart)
     */
    static public double calculateAAGroupedByPositionDistributionDistance(ProteinProfile profile1, ProteinProfile profile2, String[] groupingByPosition) {
        if (profile1.getNumColumns() != profile2.getNumColumns()) {
            throw new IllegalArgumentException("Profiles to be compared must be the same length.");
        }

        ArrayList alphabet = get20aaAlphabet();  //save the 20 aa alphabet for fast use

        //create a 'grouped' weight matrix from each profile
        SimpleWeightMatrix weightMatrix1 = profile1.getWeightMatrix();
        SimpleWeightMatrix weightMatrix2 = profile2.getWeightMatrix();
        //calculate distance between grouped matrices
        //iterate through columns
        double distance = 0.0;
        for (int j = 0; j < weightMatrix1.columns(); j++) {
            Distribution distribution1 = weightMatrix1.getColumn(j);
            //convert grouping to an array
            String[] groupingArray = groupingByPosition[j].split("[, ]");
            double[] column1 = convertDistributionToGroupColumn(groupingArray, distribution1, alphabet);
            Distribution distribution2 = weightMatrix2.getColumn(j);
            double[] column2 = convertDistributionToGroupColumn(groupingArray, distribution2, alphabet);
            double weightSum = 0.0;
            for (int k = 0; k < column1.length; k++) { //length == number of groups in grouping
                weightSum += Math.pow(column1[k] - column2[k], 2);
            }
            weightSum = weightSum / 2.0;
            distance += weightSum;
        }

        distance = distance / profile1.getNumColumns();

        return (distance);
    }

    /**
     * Calculates the distance between profiles using the following distance metric:
     * <p/>
     * <p/>
     * From: von Ohsen N, Zimmer R. Improving profile–profile alignments via
     * log average scoring. In Workshop on Algorithmic Bioinformatics,
     * 2001;11–26.
     *
     * @param profile1 Profile 1 (distance from)
     * @param profile2 Profile 2 (distance to)
     * @return The fractional distance from profile 1 to 2 (0 = identical, 1 = farthest apart)
     */
    static public double calculateLogAverageDistance(ProteinProfile profile1, ProteinProfile profile2) {
        if (profile1.getNumColumns() != profile2.getNumColumns()) {
            throw new IllegalArgumentException("Profiles to be compared must be the same length.");
        }

        ArrayList alphabet = get20aaAlphabet();  //save the 20 aa alphabet for fast use
        double distance = 0.0;

        WeightMatrix wm1 = profile1.getWeightMatrix();
        WeightMatrix wm2 = profile2.getWeightMatrix();

        SimilarityMatrix simMatrix = new SimilarityMatrix(SimilarityMatrix.BLOSUM62);

        try {
            double weightSum = 0.0;
            for (int i = 0; i < wm1.columns(); i++) {
                Distribution d1 = wm1.getColumn(i);
                Distribution d2 = wm2.getColumn(i);
                for (int j = 0; j < alphabet.size(); j++) {
                    for (int k = 0; k < alphabet.size(); k++) {
                        double weight = d1.getWeight((Symbol) alphabet.get(j)) * d2.getWeight((Symbol) alphabet.get(k));
                        weightSum += weight * simMatrix.getRawSimilarityScore((Symbol) alphabet.get(j), (Symbol) alphabet.get(k));
                    }
                }
                distance += Math.log(weightSum);
            }
        } catch (IllegalSymbolException e) {
            e.printStackTrace();
        }

        return (distance);
        //return (distance/pssm1.columns());  //average over all columns
    }

    /**
     * Calculate the euclidean distance for two vectors
     */
    static private double calcEuclidDistance(double[] vectorA, double[] vectorB) {
        //euclidean distance: Sqrt(Sum( (x[i]-y[i])^2 ))
        double distance = 0.0;
        for (int j = 0; j < vectorA.length; j++) {
            double a = vectorA[j];
            double b = vectorB[j];
            distance += Math.pow(a - b, 2);
        }
        return (Math.sqrt(distance));
    }

    /**
     * Calculates the distance between profiles using the following distance metric:
     * <p/>
     * <p/>
     * From: Atchley WR, Zhao J, Fernandes AD, Druke T.
     * Solving the protein sequence metric problem.
     * Proc Natl Acad Sci U S A. 2005 May 3;102(18):6395-400. Epub 2005 Apr 25.
     *
     * @param profile1 Profile 1 (distance from)
     * @param profile2 Profile 2 (distance to)
     * @return The fractional distance from profile 1 to 2 (0 = identical, 1 = farthest apart)
     */
    static public double calculateDrukeDistance(ProteinProfile profile1, ProteinProfile profile2) {
        if (profile1.getNumColumns() != profile2.getNumColumns()) {
            throw new IllegalArgumentException("Profiles to be compared must be the same length.");
        }

        ArrayList alphabet = get20aaAlphabet();  //save the 20 aa alphabet for fast use
        double distance = 0.0;

        WeightMatrix wm1 = profile1.getWeightMatrix();
        WeightMatrix wm2 = profile2.getWeightMatrix();

        SimilarityMatrix simMatrix = new SimilarityMatrix(SimilarityMatrix.DRUKE);

        try {
            for (int i = 0; i < wm1.columns(); i++) {
                Distribution d1 = wm1.getColumn(i);
                Distribution d2 = wm2.getColumn(i);
                for (int j = 0; j < alphabet.size(); j++) {
                    double weight1 = d1.getWeight((Symbol) alphabet.get(j));
                    double weight2 = d2.getWeight((Symbol) alphabet.get(j));
                    //multiply the Druke factors by the profile weight
                    double[] factor = simMatrix.getDrukeFactor((Symbol) alphabet.get(j));
                    double[] factor1 = new double[5];
                    System.arraycopy(factor, 0, factor1, 0, factor.length);
                    for (int l = 0; l < factor1.length; l++) {
                        factor1[l] = factor1[l] * weight1;
                    }
                    double[] factor2 = new double[5];
                    System.arraycopy(factor, 0, factor2, 0, factor.length);
                    for (int l = 0; l < factor2.length; l++) {
                        factor2[l] = factor2[l] * weight2;
                    }
                    distance += calcEuclidDistance(factor1, factor2);
                }
            }
        } catch (IllegalSymbolException e) {
            e.printStackTrace();
        }

        return (distance);
        //return (distance/pssm1.columns());  //average over all columns
    }


}
