package org.baderlab.csplugins.brainplugin;

import org.biojava.bio.BioException;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.DistributionFactory;
import org.biojava.bio.dist.DistributionTrainerContext;
import org.biojava.bio.dist.SimpleDistributionTrainerContext;
import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.io.SeqIOTools;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Iterator;

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
 * * User: GaryBader
 * * Date: Feb 6, 2005
 * * Time: 2:22:18 PM
 * * Description Calculates various distributions over a protein database
 */
public class ProteinDatabaseDistribution {
    SequenceIterator database = null;

    //TODO: seriously refactor: create 1st, 2nd and 3rd order distribution functions, none using BioJava Distribution methods
    //TODO: figure out how to use these to get a p-value for database search methods (might want to save the distribution to avoid recalculating each time)

    /**
     * Creates a new protein database object (opens the database)
     *
     * @param dbFileName The filename of the protein database to search
     * @param dbFormat   The record format of the protein database. Uses BioJava format strings, which are
     *                   (probably in the order you would want to use them): fasta, embl, swissprot, genpept, refseq, pdb, gcg, gff, nbrf, ig
     * @throws FileNotFoundException Your sequence database was not found.
     * @throws BioException          Your database format was not recognized, the format choice does not match the actual file
     *                               format or the database does not contain protein sequences in the chosen format.
     */
    public ProteinDatabaseDistribution(String dbFileName, String dbFormat) throws FileNotFoundException, BioException {
        BufferedReader br = new BufferedReader(new FileReader(dbFileName));
        database = (SequenceIterator) SeqIOTools.fileToBiojava(dbFormat, "PROTEIN", br);
    }

    private double[][] countsTableToDistribution(long[][] countsTable, long totalPairs) {
        double[][] distribution = new double[countsTable.length][countsTable[0].length];

        for (int i = 0; i < countsTable.length; i++) {
            for (int j = 0; j < countsTable[0].length; j++) {
                distribution[i][j] = (double) countsTable[i][j] / (double) totalPairs;
            }
        }

        return (distribution);
    }

    private Object[] calcOrder2PairCount(SequenceIterator searchDB, int length, ProteinTerminus terminus, int numGapsInPair) {
        //this doesn't depend on the BioJava alphabet, but is so much easier to use, and likely won't change
        final String aaList = "ACDEFGHIKLMNPQRSTUVWY";
        long[][] countsTable = new long[aaList.length()][aaList.length()];
        long totalPairs = 0;

        Sequence sequenceFromDB = null;
        while (searchDB.hasNext()) {
            try {
                sequenceFromDB = (Sequence) searchDB.nextSequence();
            } catch (BioException e) {
                System.err.println("Can't read next sequence from database.");
            }

            //optionally filter the sequence to get a specific subset
            Sequence sequenceToSearch = ProteinTerminus.getSequenceTerminus(sequenceFromDB, length, terminus);

            String sequenceToSearchString = sequenceToSearch.seqString();
            if (sequenceToSearchString.indexOf('X') < 0) { //ignore sequences with ambiguity symbols
                for (int i = 0; i < length; i++) {
                    for (int j = (i + 1); j < length; j++) {
                        if ((j - i - 1) == numGapsInPair) {
                            countsTable[aaList.indexOf(sequenceToSearchString.charAt(i))][aaList.indexOf(sequenceToSearchString.charAt(j))]++;
                            totalPairs++;
                        }
                    }
                }
            }
        }

        Object[] returnValue = new Object[3];
        returnValue[0] = countsTable;
        returnValue[1] = new Long(totalPairs);
        returnValue[2] = aaList;

        return (returnValue);
    }

    public void printDistribution(double[][] distribution, String aaList) {
        for (int i = 0; i < aaList.length(); i++) {
            if (!(aaList.charAt(i) == 'U')) {
                System.out.print("\t" + aaList.charAt(i));
            }
        }
        System.out.print("\n");
        for (int i = 0; i < distribution.length; i++) {
            if (!(aaList.charAt(i) == 'U')) {
                System.out.print(aaList.charAt(i) + "\t");
                for (int j = 0; j < distribution[i].length; j++) {
                    if (!(aaList.charAt(j) == 'U')) {
                        System.out.print(distribution[i][j]);
                        if (j < (distribution.length - 1)) {
                            System.out.print("\t");
                        }
                    }
                }
                System.out.print("\n");
            }
        }
    }

    //calculates the distribution of an orderN alphabet for protein sequences of the given dimension over a given length
    public void calcPairDistributionSearchDB(String fastaDatabaseFileName, int length, ProteinTerminus terminus) {
        SequenceIterator searchDB = null;
        BufferedReader br = null;
        long[][] countsTable;
        double[][] distribution;
        long totalPairs = 0;
        String aaList;
        Object[] returnValue;
        //cycle through all possible gap lengths between 2 aa residues
        for (int i = 0; i <= (length - 2); i++) {
            try {
                br = new BufferedReader(new FileReader(fastaDatabaseFileName));
            } catch (FileNotFoundException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }
            try {
                searchDB = (SequenceIterator) SeqIOTools.fileToBiojava("FASTA", "PROTEIN", br);
            } catch (BioException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }
            returnValue = calcOrder2PairCount(searchDB, length, terminus, i);
            countsTable = (long[][]) returnValue[0];
            totalPairs = ((Long) returnValue[1]).longValue();
            aaList = (String) returnValue[2];
            //take countsTable and make a distribution out of it
            distribution = countsTableToDistribution(countsTable, totalPairs);
            System.out.println("Gap:" + i);
            printDistribution(distribution, aaList);
            try {
                br.close();
            } catch (IOException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }
        }
    }

    //calculate the aa distribution
    //TODO: rewrite to not use BioJava Distribution methods
    public void calcAADistributionSearchDB(int length, ProteinTerminus terminus) {
        try {
            //get a DistributionTrainerContext
            DistributionTrainerContext dtc = new SimpleDistributionTrainerContext();

            Distribution proteinDist =
                    DistributionFactory.DEFAULT.createDistribution(ProteinTools.getAlphabet());

            //register the Distribution with the trainer
            dtc.registerDistribution(proteinDist);

            //for each Sequence
            Sequence sequenceFromDB = null;
            while (database.hasNext()) {
                try {
                    sequenceFromDB = (Sequence) database.nextSequence();
                } catch (BioException e) {
                    System.err.println("Can't read next sequence from database.");
                }

                //optionally filter the sequence to get a specific subset
                Sequence sequenceToSearch = ProteinTerminus.getSequenceTerminus(sequenceFromDB, length, terminus);

                //count each Symbol to the appropriate Distribution
                for (int j = 1; j <= sequenceToSearch.length(); j++) {
                    dtc.addCount(proteinDist, sequenceToSearch.symbolAt(j), 1.0);
                }
            }

            //train the Distributions
            dtc.train();

            //print the weights of each Distribution
            SymbolTokenization st = null;
            try {
                st = ProteinTools.getAlphabet().getTokenization("token");
            } catch (BioException e) {
                System.err.println("Unable to get symboltokenization");
            }
            String token = null;
            for (Iterator iter = ((FiniteAlphabet) proteinDist.getAlphabet()).iterator(); iter.hasNext();) {
                Symbol sym = (Symbol) iter.next();
                try {
                    token = st.tokenizeSymbol((AtomicSymbol) sym);
                } catch (IllegalSymbolException e) {
                    System.err.println("Unable to convert symbol to token.");
                }
                if (!token.equalsIgnoreCase("U")) {
                    System.out.println(token + "\t" + proteinDist.getWeight(sym));
                }
            }
            System.out.println("\n");

        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }
}
