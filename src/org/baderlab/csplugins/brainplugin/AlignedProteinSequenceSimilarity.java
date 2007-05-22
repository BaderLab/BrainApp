package org.baderlab.csplugins.brainplugin;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.io.SeqIOTools;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

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
 * * Time: 9:57:57 PM
 */

/**
 * Calculates the similarity of aligned protein sequences
 */
public class AlignedProteinSequenceSimilarity {
    private ArrayList sequences = null;

    /**
     * Reads a file of aligned sequences, which will be used to calculate an NxN similarity matrix
     *
     * @param alignedSequencesFileName A file of sequences in the following format:
     *                                 Name\tsequence
     *                                 Note: sequences must be perfectly aligned by some method and all the same length
     */
    public AlignedProteinSequenceSimilarity(String alignedSequencesFileName, String limitToColumns) throws IOException {
        //read proteins of interest from an aligned fasta file and save as a List
        BufferedReader brInterest = new BufferedReader(new FileReader(alignedSequencesFileName));
        try {
            SequenceIterator interestProteins = (SequenceIterator) SeqIOTools.fileToBiojava("fasta", "PROTEIN", brInterest);
            sequences = new ArrayList();
            while (interestProteins.hasNext()) {
                Sequence sequence = (Sequence) interestProteins.nextSequence();
                if (limitToColumns != null) {
                    //only consider the columns specified (e.g. a binding site)
                    String sequenceString = sequence.seqString();
                    StringBuffer sb = new StringBuffer();
                    //TODO: move this into a separate method
                    //create a new sequence string using only the specified positions of the original sequence
                    String positions[] = limitToColumns.split(",");
                    for (int i = 0; i < positions.length; i++) {
                        String position = positions[i];
                        if (position.indexOf("-") > 0) {
                            //interval
                            String[] startEnd = position.split("-");
                            int start = Integer.parseInt(startEnd[0]);
                            int end = Integer.parseInt(startEnd[1]);
                            sb.append(sequenceString.substring(start - 1, end)); //minus one to account for sequence index starting from 1
                        } else {
                            //point
                            int point = Integer.parseInt(position);
                            sb.append(sequenceString.charAt(point - 1)); //minus one to account for sequence index starting from 1
                        }
                    }
                    sequences.add(ProteinTools.createProteinSequence(sb.toString(), sequence.getName()));
                } else {
                    sequences.add(sequence);
                }
            }
        } catch (BioException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
    }

    private int calculateSimilarity(String sequenceA, String sequenceB) {
        int score = 0;
        SimilarityMatrix simMatrix = new SimilarityMatrix(SimilarityMatrix.BLOSUM62);
        for (int i = 0; i < sequenceA.length(); i++) {
            String aaA = sequenceA.substring(i, i + 1);
            String aaB = sequenceB.substring(i, i + 1);
            if (aaA.equals("-") || (aaB.equals("-"))) {
                continue;
            }
            score += simMatrix.getSimilarityScore(aaA, aaB);
        }
        return (score);
    }

    private double calculateSimilarityFactors(String sequenceA, String sequenceB) {
        double score = 0.0;
        SimilarityMatrix simMatrix = new SimilarityMatrix(SimilarityMatrix.DRUKE);
        for (int i = 0; i < sequenceA.length(); i++) {
            String aaA = sequenceA.substring(i, i + 1);
            String aaB = sequenceB.substring(i, i + 1);
            if (aaA.equals("-") || (aaB.equals("-"))) {
                continue;
            }
            double[] factorA = simMatrix.getDrukeFactor(aaA);
            double[] factorB = simMatrix.getDrukeFactor(aaB);
            //euclidean distance: Sqrt(Sum( (x[i]-y[i])^2 ))
            double distance = 0.0;
            for (int j = 0; j < factorA.length; j++) {
                double a = factorA[j];
                double b = factorB[j];
                distance += Math.pow(a - b, 2);
            }
            score += Math.sqrt(distance);
        }
        return (score / sequenceA.length());
    }

    public void calculateSimilarityMatrix() {
        AlignedProteinSequenceIdentityDistance apsid = new AlignedProteinSequenceIdentityDistance();
        for (int i = 0; i < sequences.size(); i++) {
            Sequence protein1 = (Sequence) sequences.get(i);
            for (int j = 0; j < sequences.size(); j++) {
                Sequence protein2 = (Sequence) sequences.get(j);
                System.out.print(protein1.getName() + "_" + protein2.getName() + "\t");
                System.out.println(apsid.calc(protein1.seqString(), protein2.seqString()));
                //System.out.println(calculateSimilarity(protein1.seqString(), protein2.seqString()));
            }
        }
    }


    public void calculateSimilarityFactorMatrix() {
        for (int i = 0; i < sequences.size(); i++) {
            Sequence protein1 = (Sequence) sequences.get(i);
            for (int j = 0; j < sequences.size(); j++) {
                Sequence protein2 = (Sequence) sequences.get(j);
                System.out.print(protein1.getName() + "_" + protein2.getName() + "\t");
                System.out.println(calculateSimilarityFactors(protein1.seqString(), protein2.seqString()));
            }
        }
    }
}
