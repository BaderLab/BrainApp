package org.baderlab.csplugins.brainplugin;

import org.biojava.bio.seq.Sequence;

/**
 * class to store motif hits
 */
public class Hit {
    private int start; //motif start site, in BioJava SymbolList coordinates (1..length)
    private int end;   //motif end
    private Double score; //optional score (relevant for profile search hits), null if not used
    private Sequence sequence; //the sequence that the hit came from

    /**
     * Creates a basic motif hit (start and end are in BioJava SymbolList coordinate space 1..length)
     */
    public Hit(int start, int end, Sequence sequence) {
        this.start = start;
        this.end = end;
        this.sequence = sequence;
        score = null;
    }

    /**
     * Creates an empty Hit object. Useful for calculating score statistics without
     * saving any information.
     */
    public Hit() {
        this.sequence = null;
        score = null;
    }

    /**
     * Set the optional score parameter (relevant for profile searching)
     */
    public void setScore(Double score) {
        this.score = score;
    }

    /**
     * Get the start location of this hit on the sequence
     */
    public int getStart() {
        return start;
    }

    /**
     * Get the end location of this hit on the sequence
     */
    public int getEnd() {
        return end;
    }

    /**
     * Get the score for this hit, if derived from a profile method
     *
     * @return null if no score defined
     */
    public Double getScore() {
        return score;
    }

    /**
     * Get the sequence that was hit
     */
    public Sequence getSequence() {
        return sequence;
    }

    /**
     * Get a string representation of the pattern match
     * Does not include the start and end location, just the amino acid sequence
     */
    public String getMatchString() {
        return (sequence.subStr(start, end));
    }

    /**
     * Returns a basic string representation of this Hit
     */
    public String toString() {
        StringBuffer sb = new StringBuffer();
        sb.append(this.getScore() + "\t" + this.getSequence().seqString() + "\t(" + this.getStart() + "-" + this.getEnd() + ")");
        return sb.toString();
    }
}
