package org.baderlab.csplugins.brainplugin;

import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.impl.SimpleSequence;
import org.biojava.bio.symbol.SymbolList;

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
 * * Time: 10:00:56 PM
 * * Description The type of protein terminus (N, C or None)
 */
public final class ProteinTerminus implements Comparable {
    public final static int TERMINUS_NONE = 0;
    public final static int TERMINUS_N = 1;
    public final static int TERMINUS_C = 2;

    public static final ProteinTerminus NONE = new ProteinTerminus(TERMINUS_NONE);
    public static final ProteinTerminus N = new ProteinTerminus(TERMINUS_N);
    public static final ProteinTerminus C = new ProteinTerminus(TERMINUS_C);

    private static final ProteinTerminus Termini[] = {
            NONE, N, C
    };

    private static final String[] TerminiStrings = {
            "NONE", "N", "C"
    };

    private int terminus;

    public static final ProteinTerminus getProteinTerminus(int terminus) {
        if (terminus >= TERMINUS_NONE && terminus <= TERMINUS_C)
            return Termini[terminus];
        return null;
    }

    private ProteinTerminus(int terminus) {
        this.terminus = terminus;
    }

    public boolean equals(Object obj) {
        // Assume proper type was given
        return (terminus == ((ProteinTerminus) obj).terminus);
    }

    public int compareTo(Object obj) {
        // Assume proper type was given
        int other = ((ProteinTerminus) obj).terminus;
        if (terminus == other)
            return 0;
        if (terminus < other)
            return -1;
        return 1;
    }

    public int getTerminus() {
        return terminus;
    }

    /**
     * parse a string describing the terminus and return the appropriate ProteinTerminus instance
     *
     * @param terminus A String e.g. C, N, NONE.
     * @return Default is NONE
     */
    public static ProteinTerminus parseTerminus(String terminus) {
        if (terminus.equalsIgnoreCase("C")) {
            return C;
        }
        if (terminus.equalsIgnoreCase("N")) {
            return N;
        }
        if (terminus.equalsIgnoreCase("NONE")) {
            return NONE;
        }
        return NONE; //default
    }

    /**
     * Returns a sequence subset corresponding to the first or last 'length' residues
     *
     * @param sequence           The input sequence to filter
     * @param lengthFromTerminus The size of the subset in amino acid residues
     * @param terminus           One of either TERMINUS_N or TERMINUS_C constants. If TERMINUS_NONE is passed, the sequence will be returned.
     */
    public static Sequence getSequenceTerminus(Sequence sequence, int lengthFromTerminus, ProteinTerminus terminus) {
        Sequence sequenceSubSet = null;
        if ((terminus.equals(N)) && (sequence.length() >= lengthFromTerminus)) {
            //only do first n residues
            SymbolList firstNresidues = sequence.subList(1, lengthFromTerminus);
            //just make a basic copy, don't maintain the annotation and features
            sequenceSubSet = new SimpleSequence(firstNresidues, "", sequence.getName(), null);
        } else if ((terminus.equals(C)) && (sequence.length() >= lengthFromTerminus)) {
            //only do last n residues
            SymbolList lastNresidues = sequence.subList(sequence.length() - (lengthFromTerminus - 1), sequence.length());
            //just make a basic copy, don't maintain the annotation and features
            sequenceSubSet = new SimpleSequence(lastNresidues, "", sequence.getName(), null);
        } else {
            //use the full sequence
            sequenceSubSet = sequence;
        }
        return (sequenceSubSet);
    }

    /**
     * Returns a short string representation of the terminus
     */
    public String toString() {
        return new String(TerminiStrings[terminus]);
    }
}
