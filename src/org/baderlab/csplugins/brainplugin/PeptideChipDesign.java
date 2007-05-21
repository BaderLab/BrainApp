package org.baderlab.csplugins.brain;

import org.biojava.bio.seq.Sequence;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;

/**
 * Copyright (c) 2004 Memorial Sloan-Kettering Cancer Center
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
 * * User: Gary Bader
 * * Date: Oct 15, 2004
 * * Time: 10:48:54 AM
 * * Description Stores peptides for a peptide chip design and handles common chip design operations
 */

/**
 * Handles peptide chip design operations
 * This class handles uniquing the peptide list, so it does not need to be done by calling code.
 */
public class PeptideChipDesign {
    HashMap aaToCodonMap; //a map of amino acid letters to preferred codon (in organism to be used to create peptides)
    SequenceSearchResultSet peptides; //a list of peptides on the chip
    String description; //a free text description of the chip

    /**
     * Create a new peptide chip design
     *
     * @param codonUsageTableFileName The preferred codon bias table in a file (tab delimited)
     *                                amino acid letter(tab)codon e.g. A(tab)GCU (one per line)
     * @throws IOException if there is an IO problem with the codon bias file
     */
    public PeptideChipDesign(String codonUsageTableFileName) throws IOException {
        aaToCodonMap = new HashMap();
        readAaToCodonMap(codonUsageTableFileName);
        this.description = null;
    }

    /**
     * Create a new peptide chip design with an initial set of peptides
     *
     * @param codonUsageTableFileName The preferred codon bias table in a file (tab delimited)
     *                                amino acid letter(tab)codon e.g. A(tab)GCU (one per line)
     * @param peptides                The initial peptide set
     * @throws IOException if there is an IO problem with the codon bias file
     */
    public PeptideChipDesign(String codonUsageTableFileName, SequenceSearchResultSet peptides) throws IOException {
        this(codonUsageTableFileName);
        //Make sure that the peptides are unique. We want to guarantee not to waste space on the chip
        this.peptides = peptides.getUniqueResultsBySequence();
    }

    public void setDescription(String description) {
        this.description = description;
    }

    /**
     * Read the amino acid to codon map
     *
     * @param fileName file name of the aa to codon map
     * @throws IOException if there is an IO problem with the codon bias file
     */
    private void readAaToCodonMap(String fileName) throws IOException {
        final String aaList = "ACDEFGHIKLMNPQRSTUVWY";
        BufferedReader br = null;
        br = new BufferedReader(new FileReader(fileName));
        String line = null;
        String[] lineSplit;
        while ((line = br.readLine()) != null) {
            lineSplit = line.split("\t");
            if (lineSplit.length == 2) {
                String aa = lineSplit[0];
                String codon = lineSplit[1];
                //do some basic error checking
                //test if the aa is recognized
                if (aaList.indexOf(aa) < 0) {
                    throw new IOException("Non standard amino acid found in codon bias table (" + aa + ").");
                }
                if (codon.length() != 3) {
                    throw new IOException("Codon found not of length 3 (" + codon + ").");
                }
                //ignore all lines that don't split into 2
                aaToCodonMap.put(aa, codon);
            }
        }
        br.close();
    }

    /**
     * Add a set of sequence search results to the chip. This method handles uniquing the peptide list, so
     * it does not need to be done by calling code.
     *
     * @param peptideResults The list of peptides to add
     */
    public void addToChip(SequenceSearchResultSet peptideResults) {
        if (peptides == null) {
            peptides = peptideResults;
        } else {
            //already have some results, add new results to existing ones
            peptides.mergeResultsFrom(peptideResults);
        }
        //Make sure that the peptides are unique. We want to guarantee not to waste space on the chip
        peptides = peptides.getUniqueResultsBySequence();
    }

    /**
     * Get the size of the chip (the number of peptides currently on the chip)
     * Peptides are guaranteed unique.
     */
    public int size() {
        return peptides.getNumberSequencesHit();
    }

    /**
     * Converts the peptide chip design to a user friendly string (a table decribing the chip)
     *
     * @return a string representation of the peptide chip
     */
    public String toString() {
        StringBuffer sb = new StringBuffer();
        if (description != null) {
            sb.append("Description: " + description + "\n");
        }
        sb.append("Peptide\tName\tDNA\n");
        Set peptideSet = peptides.getSequences();
        for (Iterator iterator = peptideSet.iterator(); iterator.hasNext();) {
            Sequence sequence = (Sequence) iterator.next();
            sb.append(sequence.seqString() + "\t" + sequence.getName() + "\t" + convertPeptideToDNA(sequence.seqString()) + "\n");
        }
        return (sb.toString());
    }

    /**
     * Converts a peptide sequence to a DNA sequence based on the registered codon bias map in this class.
     *
     * @param peptide The peptide to convert
     * @return a string of the DNA sequence
     */
    private String convertPeptideToDNA(String peptide) {
        //take peptide and replace all aa's with codons from the preferred codon table.
        StringBuffer DNA = new StringBuffer();
        char[] peptideCharArray = peptide.toCharArray();
        for (int i = 0; i < peptideCharArray.length; i++) {
            char c = peptideCharArray[i];
            String codon = (String) aaToCodonMap.get(String.valueOf(c));
            DNA = DNA.append(codon);
        }
        return (DNA.toString());
    }
}
