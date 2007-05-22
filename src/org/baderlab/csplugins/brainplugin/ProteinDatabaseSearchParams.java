package org.baderlab.csplugins.brainplugin;

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
 * * Date: Feb 17, 2005
 * * Time: 10:21:32 PM
 * * Description Stores the parameters for a protein database search.
 */

/**
 * Stores the parameters for a protein database search.
 */
public class ProteinDatabaseSearchParams {
    private ProteinTerminus terminus = null;
    private Integer length = null;
    private boolean multipleHits = false; //true=multiple hits allowed
    private Boolean normalized = null;  //true=threshold is on the profile normalized p-value
    //note: normalized is an object so we can determine if it was set (null = not set)
    private int scoreType = -1; //the type of score
    public final static int SCORE_TYPE_PROBABILITY = 1;
    public final static int SCORE_TYPE_NORM_PROBABILITY = 2;
    public final static int SCORE_TYPE_DELTAG = 3;
    private boolean dontSaveSequences = false; //true = don't save sequences (useful for score statistics)

    /**
     * Use for a whole sequence search.
     *
     * @param multipleHits If true, then multiple hits are allowed, otherwise only the first hit is saved.
     */
    public ProteinDatabaseSearchParams(boolean multipleHits) {
        this.multipleHits = multipleHits;
        terminus = ProteinTerminus.NONE;
    }

    /**
     * Use for a multiple hit search over a length of sequence that starts from the N or C terminus
     *
     * @param terminus The terminus to search
     * @param length   The length of sequence from terminus to search
     */
    public ProteinDatabaseSearchParams(ProteinTerminus terminus, int length) {
        this.terminus = terminus;
        this.length = new Integer(length);
        multipleHits = true;
    }

    /**
     * Use for a single hit search at the C or N terminus based on the length of the profile or regex
     *
     * @param terminus The terminus to search from
     */
    public ProteinDatabaseSearchParams(ProteinTerminus terminus) {
        this.terminus = terminus;
        multipleHits = false;
    }

    /**
     * Set whether the score threshold will be normalized to the profile or not.
     * Default is false. This is only used for profile searching. Regex searching will
     * ignore this parameter.
     *
     * @see ProteinProfile#getNormalizedPValue(double)
     */
    public void setNormalized(boolean normalized) {
        this.normalized = new Boolean(normalized);
    }

    /**
     * Sets the type of score to be selected. Score types are defined in this class.
     */
    public void setScoreType(int scoreType) {
        this.scoreType = scoreType;
    }

    /**
     * Gets the type of score selected
     */
    public int getScoreType() {
        return scoreType;
    }

    /**
     * Gets the protein terminus choice
     */
    public ProteinTerminus getTerminus() {
        return terminus;
    }

    /**
     * Gets the length of sequence to search if set
     *
     * @return the sequence length.  If not set, returns -1.
     */
    public int getLength() {
        if (length == null) {
            return -1;
        }
        return length.intValue();
    }

    /**
     * Gets the multiple hit choice
     */
    public boolean isMultipleHits() {
        return multipleHits;
    }

    /**
     * Gets the normalized threshold choice
     */
    public boolean isNormalized() {
        if (normalized == null) {
            return false;
        }
        return normalized.booleanValue();
    }

    /**
     * Useful for GUI code, which might want to set a default value depending on whether normalized has been set
     */
    public boolean isNormalizedSet() {
        return normalized != null;
    }

    /**
     * Set the terminus to search
     */
    public void setTerminus(ProteinTerminus terminus) {
        this.terminus = terminus;
    }

    /**
     * Set the length to search - only makes sense for whole sequence and multiple hit terminal searches
     * If length==null, indicates that length is not set.
     */
    public void setLength(Integer length) {
        this.length = length;
    }

    /**
     * Set whether the search should return multiple hits
     */
    public void setMultipleHits(boolean multipleHits) {
        this.multipleHits = multipleHits;
    }

    /**
     * Gets the save sequences choice (useful for score statistics)
     */
    public boolean isDontSaveSequences() {
        return dontSaveSequences;
    }

    /**
     * Set whether the search should not save sequences
     * true = don't save sequences (useful for score statistics)
     */
    public void setDontSaveSequences(boolean dontSaveSequences) {
        this.dontSaveSequences = dontSaveSequences;
    }

    /**
     * Create a copy of this parameter object
     */
    public ProteinDatabaseSearchParams copy() {
        ProteinDatabaseSearchParams newParams = new ProteinDatabaseSearchParams(this.multipleHits);
        newParams.normalized = this.normalized;
        newParams.length = this.length;
        newParams.terminus = this.terminus;
        newParams.scoreType = this.scoreType;
        return (newParams);
    }

    public String toString() {
        String lineSep = System.getProperty("line.separator");
        StringBuffer sb = new StringBuffer();
        sb.append("Sequence search: " + ((terminus.equals(ProteinTerminus.NONE)) ? " Whole sequence" : " Sequence subset ("
                + terminus + " terminus " + length + " residues)") + lineSep);
        sb.append("Multiple hits allowed: " + ((multipleHits) ? "True" : "False") + lineSep);
        sb.append("Score type: ");
        switch (scoreType) {
            case SCORE_TYPE_PROBABILITY:
                sb.append("Probability");
            case SCORE_TYPE_NORM_PROBABILITY:
                sb.append("Normalized Probability");
            case SCORE_TYPE_DELTAG:
                sb.append("Delta G Style");
            default:
                sb.append("Unknown score type");
        }
        sb.append(lineSep);
        return sb.toString();
    }
}
