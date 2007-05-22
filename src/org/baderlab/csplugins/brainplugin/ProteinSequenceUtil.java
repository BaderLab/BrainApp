package org.baderlab.csplugins.brainplugin;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.Symbol;

import java.util.Arrays;
import java.util.HashMap;
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
 * * Date: Dec 7, 2005
 * * Time: 4:36:34 PM
 */

/**
 * Provides some utility methods for working with generic protein sequences
 */
public class ProteinSequenceUtil {

    /**
     * Filter a sequence, choosing only specific positions. E.g. select the binding site from a full sequence.
     *
     * @param sequence The sequence to filter
     * @param filter   A filter for sequence positions e.g. "1-4,5,8" will return a string containing positions 1 to 4
     *                 followed by position 5 followed by position 8.
     * @return A string containing the result of the sequence filtering
     */
    public static String filterSequenceByColumns(Sequence sequence, String filter) {
        //only consider the columns specified (e.g. a binding site)
        String sequenceString = sequence.seqString();
        return filterSequenceStringByColumns(sequenceString, filter);
    }

    /**
     * Internal helper method for filterSequenceByColumns - same arguments.
     */
    private static String filterSequenceStringByColumns(String sequenceString, String filter) {
        StringBuffer sb = new StringBuffer();
        //create a new sequence string using only the specified positions of the original sequence
        String positions[] = filter.split(",");
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
        return sb.toString();
    }

    /**
     * Returns the length of the string that would result when using the filterSequenceByColumns method
     * using the given filter
     *
     * @param filter The filter that will be used to filter a string
     * @return The length of the filtered string
     */
    public static int countLengthOfFilteredStringResult(String filter) {
        //determine last position
        String positions[] = filter.split("[,-]");
        int lastPosition = Integer.parseInt(positions[positions.length - 1]);
        //create a string containing at least 'lastPosition+1' letters
        char[] letters = new char[lastPosition + 1];
        Arrays.fill(letters, 'a');
        String dummyString = new String(letters);
        //filter it
        String newDummyString = filterSequenceStringByColumns(dummyString, filter);
        //return the length of the resulting string
        return newDummyString.length();
    }

    /**
     * Saves a 20aa biojava alphabet for ease of use (ignore the U symbol)
     *
     * @return A HashMap where the key is the String of the single letter IUPAC AA symbol and the
     *         value is the corresponding BioJava Symbol object from the standard protein alphabet
     */
    public static HashMap get20aaAlphabet() {
        HashMap alphabet = new HashMap(20);
        FiniteAlphabet symbols = ProteinTools.getAlphabet();
        Iterator symbolIterator = symbols.iterator();
        while (symbolIterator.hasNext()) {
            Symbol symbol = (Symbol) symbolIterator.next();
            String symbolString = null;
            try {
                symbolString = symbols.getTokenization("token").tokenizeSymbol(symbol);
            } catch (BioException e) {
                e.printStackTrace(); //this should never happen, since we are using a standard alphabet
            }
            if (!symbolString.equals("U")) {
                alphabet.put(symbolString, symbol);
            }
        }
        return alphabet;
    }
}
