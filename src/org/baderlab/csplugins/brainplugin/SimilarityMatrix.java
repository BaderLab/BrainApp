package org.baderlab.csplugins.brainplugin;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.symbol.Symbol;

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
 * * Date: Jun 19, 2005
 * * Time: 4:11:03 PM
 */

/**
 * Encapsulates an amino acid similarity matrix
 */
public class SimilarityMatrix {
    final static int BLOSUM62 = 1;
    final static int DRUKE = 100;
    private double[][] matrix = null; //matrix of scaled log odds aa similarity scores
    private double scaleFactor = 0.0; //scale factor to multiply the matrix values by to get the actual log odds scores
    //lookup index for all matrices
    final String aaList = "ACDEFGHIKLMNPQRSTVWY";

    public SimilarityMatrix(int type) {
        switch (type) {
            case BLOSUM62:
                matrix = getBlosum62();
                scaleFactor = Math.log(2) / 2;
                break;
            case DRUKE:
                matrix = getDrukeFactors();
                scaleFactor = 1;  //scale factor is not important here
                break;
        }
    }

    private String symbol2string(Symbol s) {
        String symbolToken = null;
        try {
            symbolToken = ProteinTools.getAlphabet().getTokenization("token").tokenizeSymbol(s);
        } catch (BioException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
        return (symbolToken);
    }

    /**
     * Get the similarity score for two amino acids in the similarity matrix
     *
     * @param aSymbol amino acid A
     * @param bSymbol amino acid B
     * @return The similarity score
     */
    public double getSimilarityScore(Symbol aSymbol, Symbol bSymbol) {
        String a = symbol2string(aSymbol);
        String b = symbol2string(bSymbol);
        return ((double) matrix[aaList.indexOf(a)][aaList.indexOf(b)]);
    }

    /**
     * Get the similarity score for two amino acids in the similarity matrix
     *
     * @param a amino acid A
     * @param b amino acid B
     * @return The similarity score
     */
    public double getSimilarityScore(String a, String b) {
        return ((double) matrix[aaList.indexOf(a)][aaList.indexOf(b)]);
    }

    /**
     * Get the raw similarity score for two amino acids in the similarity matrix (Pij/PiPj)
     *
     * @param aSymbol amino acid A
     * @param bSymbol amino acid B
     * @return The raw similarity score (observed over expected = odds)
     */
    public double getRawSimilarityScore(Symbol aSymbol, Symbol bSymbol) {
        double score = getSimilarityScore(aSymbol, bSymbol);
        score = score * scaleFactor;
        double rawScore = Math.exp(score);
        return (rawScore);
    }

    /**
     * Internal representation of the BLOSUM 62 matrix
     * (Values copied from NCBI BLAST data directory)
     */
    private double[][] getBlosum62() {
        //values sorted in Excel to match aaList index
        double[][] matrix = {
                {4, 0, -2, -1, -2, 0, -2, -1, -1, -1, -1, -2, -1, -1, -1, 1, 0, 0, -3, -2},
                {0, 9, -3, -4, -2, -3, -3, -1, -3, -1, -1, -3, -3, -3, -3, -1, -1, -1, -2, -2},
                {-2, -3, 6, 2, -3, -1, -1, -3, -1, -4, -3, 1, -1, 0, -2, 0, -1, -3, -4, -3},
                {-1, -4, 2, 5, -3, -2, 0, -3, 1, -3, -2, 0, -1, 2, 0, 0, -1, -2, -3, -2},
                {-2, -2, -3, -3, 6, -3, -1, 0, -3, 0, 0, -3, -4, -3, -3, -2, -2, -1, 1, 3},
                {0, -3, -1, -2, -3, 6, -2, -4, -2, -4, -3, 0, -2, -2, -2, 0, -2, -3, -2, -3},
                {-2, -3, -1, 0, -1, -2, 8, -3, -1, -3, -2, 1, -2, 0, 0, -1, -2, -3, -2, 2},
                {-1, -1, -3, -3, 0, -4, -3, 4, -3, 2, 1, -3, -3, -3, -3, -2, -1, 3, -3, -1},
                {-1, -3, -1, 1, -3, -2, -1, -3, 5, -2, -1, 0, -1, 1, 2, 0, -1, -2, -3, -2},
                {-1, -1, -4, -3, 0, -4, -3, 2, -2, 4, 2, -3, -3, -2, -2, -2, -1, 1, -2, -1},
                {-1, -1, -3, -2, 0, -3, -2, 1, -1, 2, 5, -2, -2, 0, -1, -1, -1, 1, -1, -1},
                {-2, -3, 1, 0, -3, 0, 1, -3, 0, -3, -2, 6, -2, 0, 0, 1, 0, -3, -4, -2},
                {-1, -3, -1, -1, -4, -2, -2, -3, -1, -3, -2, -2, 7, -1, -2, -1, -1, -2, -4, -3},
                {-1, -3, 0, 2, -3, -2, 0, -3, 1, -2, 0, 0, -1, 5, 1, 0, -1, -2, -2, -1},
                {-1, -3, -2, 0, -3, -2, 0, -3, 2, -2, -1, 0, -2, 1, 5, -1, -1, -3, -3, -2},
                {1, -1, 0, 0, -2, 0, -1, -2, 0, -2, -1, 1, -1, 0, -1, 4, 1, -2, -3, -2},
                {0, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1, 0, -1, -1, -1, 1, 5, 0, -2, -2},
                {0, -1, -3, -2, -1, -3, -3, 3, -2, 1, 1, -3, -2, -2, -3, -2, 0, 4, -3, -1},
                {-3, -2, -4, -3, 1, -2, -2, -3, -3, -2, -1, -4, -4, -2, -3, -3, -2, -3, 11, 2},
                {-2, -2, -3, -2, 3, -3, 2, -1, -2, -1, -1, -2, -3, -1, -2, -2, -2, -1, 2, 7}
        };
        return matrix;
    }

    public double[] getDrukeFactor(String aa) {
        return matrix[aaList.indexOf(aa)];
    }

    public double[] getDrukeFactor(Symbol aa) {
        String aaString = null;
        try {
            aaString = ProteinTools.getAlphabet().getTokenization("token").tokenizeSymbol(aa);
        } catch (BioException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
        if (aaString == null) {
            return null;
        }
        return matrix[aaList.indexOf(aaString)];
    }

    /**
     * Factor analysis of amino acid physical properties indices from:
     * Solving the protein sequence metric problem.
     * Atchley WR, Zhao J, Fernandes AD, Druke T.
     * Proc Natl Acad Sci U S A. 2005 May 3;102(18):6395-400. Epub 2005 Apr 25.
     * PMID: 15851683
     *
     * @return a matrix of 5-dimensional vectors, one for each amino acid
     */
    private double[][] getDrukeFactors() {
        double[][] matrix = {
                //Factor I (polarity)   Factor II (secondary structure)
                //Factor III (molecular size)   Factor IV (relative aa composition)
                //Factor V (electrostatic charge)
                /*A*/ {-0.591, -1.302, -0.733, 1.570, -0.146},
                /*C*/ {-1.343, 0.465, -0.862, -1.020, -0.255},
                /*D*/ {1.050, 0.302, -3.656, -0.259, -3.242},
                /*E*/ {1.357, -1.453, 1.477, 0.113, -0.837},
                /*F*/ {-1.006, -0.590, 1.891, -0.397, 0.412},
                /*G*/ {-0.384, 1.652, 1.330, 1.045, 2.064},
                /*H*/ {0.336, -0.417, -1.673, -1.474, -0.078},
                /*I*/ {-1.239, -0.547, 2.131, 0.393, 0.816},
                /*K*/ {1.831, -0.561, 0.533, -0.277, 1.648},
                /*L*/ {-1.019, -0.987, -1.505, 1.266, -0.912},
                /*M*/ {-0.663, -1.524, 2.219, -1.005, 1.212},
                /*N*/ {0.945, 0.828, 1.299, -0.169, 0.933},
                /*P*/ {0.189, 2.081, -1.628, 0.421, -1.392},
                /*Q*/ {0.931, -0.179, -3.005, -0.503, -1.853},
                /*R*/ {1.538, -0.055, 1.502, 0.440, 2.897},
                /*S*/ {-0.228, 1.399, -4.760, 0.670, -2.647},
                /*T*/ {-0.032, 0.326, 2.213, 0.908, 1.313},
                /*V*/ {-1.337, -0.279, -0.544, 1.242, -1.262},
                /*W*/ {-0.595, 0.009, 0.672, -2.128, -0.184},
                /*Y*/ {0.260, 0.830, 3.097, -0.838, 1.512}
        };

        return matrix;
    }
}
