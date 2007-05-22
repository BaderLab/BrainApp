package org.baderlab.brain;

import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.baderlab.brain.util.CRC64;

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
 * * Time: 5:37:42 PM
 * * Description Stores the results of a sequence search
 */

/**
 * Stores the results of a sequence search
 */
public class SequenceSearchResultSet {
    private ProteinDatabaseSearchParams params = null; //the general parameters used to get sequence length
    private long numberOfSequencesSearched = 0;
    private long numberOfHits = 0;    //could be multiple hits per sequence
    private HashMap results = null; //key=Sequence, value=ArrayList of Hits
    private HashMap hashResults = null; //key = sequenceDigest, value = Sequence object
    private HashMap annotatedSequenceMap = null; //key = Sequence hit, value = original annotated Sequence object
    private TreeMap sortedResultMap = null; //key=score (Double), value=ArrayList of Hits
    private String species; //an organism name to optionally store.  Hits are from this organism.
    //profile search relevant
    private ProteinProfile profile = null; //exists if the sequence search used a profile
    private double scoreThreshold = 0.0;
    //regular expression search relevant
    private String regex = null; //exists if the sequence search used a regular expressions (only one regex per result set to avoid confusion)

    /**
     * Create an empty sequence search result set (profile).
     *
     * @param scoreThreshold The score threshold for the profile search
     * @param params         General search parameters
     */
    public SequenceSearchResultSet(double scoreThreshold, ProteinDatabaseSearchParams params) {
        this.scoreThreshold = scoreThreshold;
        this.params = params;
        results = new HashMap();
        hashResults = new HashMap();
    }

    /**
     * Create an empty sequence search result set (regex).
     *
     * @param params General search parameters
     */
    public SequenceSearchResultSet(ProteinDatabaseSearchParams params) {
        this.params = params;
        results = new HashMap();
        hashResults = new HashMap();
    }

    /**
     * Create a hash of the sequence for the results hashmap.
     * This is required because we can't guarantee that the sequences being saved are the
     * same object - they may be new objects created from an existing sequence e.g. for sequence
     * filtering reasons.  The hash is computed based on the actual sequence data + name, so filtered
     * sequences should be recognized as identical if they have the same name and sequence
     *
     * @param sequence Sequence to compute a hash for
     * @return The sequence hash
     */
    private String uniqueSequenceHash(Sequence sequence) {
        CRC64 crc = new CRC64();
        crc.update(sequence.getName());
        crc.update(sequence.seqString());
        long digest = crc.getValue();
        return String.valueOf(digest);
    }

    /**
     * Check if we have seen this sequence before and return the previously seen sequence
     * if we have.
     *
     * @param sequenceHit sequence to check if we've already seen
     * @param remember    If true, will remember the sequence hit
     * @return The previously seen sequence. If the sequence has not been seen and remember=true, the given
     *         sequence will be returned, otherwise null will be returned.
     */
    private Sequence findLookupSequence(Sequence sequenceHit, boolean remember) {
        Sequence lookupSequence = null;
        if (!hashResults.containsKey(uniqueSequenceHash(sequenceHit))) {
            //new sequence
            if (remember) {
                hashResults.put(uniqueSequenceHash(sequenceHit), sequenceHit);
                lookupSequence = sequenceHit;
            }
        } else {
            lookupSequence = (Sequence) hashResults.get(uniqueSequenceHash(sequenceHit));
        }
        return lookupSequence;
    }

    /**
     * Add a sequence hit (pattern match) to the result set. Multiple hits per sequence are supported.
     *
     * @param sequenceHit The sequence found with a hit
     * @param hitStart    The start of the hit
     * @param hitEnd      The end of the hit
     * @param score       The score of this hit (according to the profile)
     */
    public void addSequenceHit(Sequence sequenceHit, int hitStart, int hitEnd, Double score) {
        ArrayList hits = null;
        //check if we have already seen this sequence
        Sequence lookupSequence = findLookupSequence(sequenceHit, true);
        //check if we need to create a new hits array (yes, if multiple hits arrive for the same sequence)
        if (!results.containsKey(lookupSequence)) {
            //new sequence
            hits = new ArrayList();
        } else {
            //existing sequence
            hits = (ArrayList) results.get(lookupSequence);
        }
        //check if we have already seen this score in the score map
        ArrayList scoreHits = null;
        if (sortedResultMap != null) {
            if (!sortedResultMap.containsKey(score)) {
                //new sequence
                scoreHits = new ArrayList();
            } else {
                //existing sequence
                scoreHits = (ArrayList) sortedResultMap.get(score);
            }
        }
        Hit hit = new Hit(hitStart, hitEnd, lookupSequence);
        if (score != null) {
            hit.setScore(score);
        }
        hits.add(hit);
        //maintain the results map
        results.put(lookupSequence, hits);
        //also maintain the score map
        if (scoreHits != null) {
            scoreHits.add(hit);
            sortedResultMap.put(score, scoreHits);
        }
        numberOfHits++;
    }

    /**
     * Add a sequence hit (pattern match) to the result set. Multiple hits per sequence are supported.
     *
     * @param sequenceHit The sequence found with a hit
     * @param hitStart    The start of the hit
     * @param hitEnd      The end of the hit
     */
    public void addSequenceHit(Sequence sequenceHit, int hitStart, int hitEnd) {
        addSequenceHit(sequenceHit, hitStart, hitEnd, null);
    }

    /**
     * Add a sequence hit (pattern match) to the result set. Multiple hits per sequence are supported.
     * Note: This method is only useful for calculating score statistics. It does not save important
     * information required for search reporting - just the score.
     *
     * @param score The score of this hit (according to the profile)
     */
    public void addSequenceHit(Double score) {
        if (score == null) {
            throw new IllegalArgumentException("Score was null.");
        }
        Hit hit = new Hit();
        hit.setScore(score);
        //check if we have already seen this score in the score map
        ArrayList scoreHits = null;
        if (sortedResultMap != null) {
            if (!sortedResultMap.containsKey(score)) {
                //new sequence
                scoreHits = new ArrayList();
            } else {
                //existing sequence
                scoreHits = (ArrayList) sortedResultMap.get(score);
            }
            scoreHits.add(hit);
            sortedResultMap.put(score, scoreHits);
        }
        numberOfHits++;
    }

    /**
     * Add a reference to the original database sequence that generated this hit.
     * This is used to only keep one copy of the database sequence around to save memory
     * If the sequenceHit is not from this result set, nothing will be added
     *
     * @param sequenceHit      The hit
     * @param originalSequence The original sequence
     */
    public void addOriginalSequenceToHit(Sequence sequenceHit, Sequence originalSequence) {
        Sequence lookupSequence = findLookupSequence(sequenceHit, false);
        if (lookupSequence == null) {
            return;
        }
        if (annotatedSequenceMap == null) {
            annotatedSequenceMap = new HashMap();
        }
        annotatedSequenceMap.put(lookupSequence, originalSequence);
    }

    /**
     * Set the number of sequences that have been searched that resulted in this result set.
     */
    public void setNumberOfSequencesSearched(long numberOfSequencesSearched) {
        this.numberOfSequencesSearched = numberOfSequencesSearched;
    }

    /**
     * If the search was a profile search, set the profile that was used.
     */
    public void setProfile(ProteinProfile profile) {
        this.profile = profile;
        //this is only needed if we have scores from a profile
        //scores (-log10(p-value)) are sorted in ascending order
        sortedResultMap = new TreeMap();
    }

    /**
     * If the search was a regular expression search, set the regex that was used
     */
    public void setRegex(String regex) {
        this.regex = regex;
    }

    /**
     * Get the sequences in the result set (the results of the database search)
     *
     * @return a set of BioJava Sequence objects representing the search results
     */
    public Set getSequences() {
        return results.keySet();
    }

    /**
     * Get the species name for this result set, if it exists
     *
     * @return The species name or null if it is not set.
     */
    public String getSpecies() {
        return species;
    }

    /**
     * Sets the species name for this result set (optional)
     *
     * @param species The species (organism) name
     */
    public void setSpecies(String species) {
        this.species = species;
    }

    /**
     * Gets the original database sequence tied to this hit
     * This needs to be used to get the original annotation for the hit sequence
     *
     * @param sequenceHit The sequence hit (from this result set) to lookup the original sequence
     * @return The original database sequence, or null if it is not tied to this hit.
     */
    public Sequence getOriginalSequence(Sequence sequenceHit) {
        return (Sequence) annotatedSequenceMap.get(sequenceHit);
    }

    /**
     * Get the number of sequences in the result set
     *
     * @return the number of database sequences hit
     */
    public int getNumberSequencesHit() {
        return results.size();
    }

    /**
     * Get a list of Hit objects corresponding to the sequence
     *
     * @param sequence Get hits for this sequence
     * @return a list of Hit objects
     */
    public List getHits(Sequence sequence) {
        return ((ArrayList) results.get(sequence));
    }

    /**
     * Gets the number of sequences searched (usually the database size)
     */
    public long getNumberOfSequencesSearched() {
        return numberOfSequencesSearched;
    }

    /**
     * Gets the number of hits in this result set
     */
    public long getNumberOfHits() {
        return numberOfHits;
    }

    /**
     * Gets the profile that was used to generate this result set (if present)
     *
     * @return The profile, otherwise null
     */
    public ProteinProfile getProfile() {
        return profile;
    }

    /**
     * Gets the regex that was used to generate this result set (if present)
     *
     * @return The regular expression, otherwise null
     */
    public String getRegex() {
        return regex;
    }

    /**
     * Gets the score threshold that was used in the search resulting in this result set
     *
     * @return The double value originally set for this result set
     */
    public double getScoreThreshold() {
        return scoreThreshold;
    }

    /**
     * Gets a unique set of resulting hits by sequence. Duplicate sequences could exist in the database
     * that generate exactly the same hits. This method finds identical sequences by exact sequence match.
     *
     * @return a set of search results, unique by sequence match
     */
    public SequenceSearchResultSet getUniqueResultsBySequence() {
        SequenceSearchResultSet uniqueResults = new SequenceSearchResultSet(scoreThreshold, params);
        uniqueResults.setProfile(profile);
        uniqueResults.setRegex(regex);
        uniqueResults.setNumberOfSequencesSearched(this.numberOfSequencesSearched);
        //collect all sequence strings in a hash and then recreate a unique result set
        //this code assumes that two identical sequence strings will have the same set of hits (which should be the case)
        HashMap uniqueSequenceStrings = new HashMap();
        Set sequences = results.keySet();
        for (Iterator iterator = sequences.iterator(); iterator.hasNext();) {
            Sequence sequence = (Sequence) iterator.next();
            String seqString = sequence.seqString();
            //array[2] of old key and value
            Object[] oldResultObject = {sequence, results.get(sequence)};
            //rename sequence to keep track of all data
            if (uniqueSequenceStrings.containsKey(seqString)) {
                //update the existing sequence to add new name
                Sequence oldSequence = (Sequence) ((Object[]) uniqueSequenceStrings.get(seqString))[0];
                try {
                    oldResultObject[0] = ProteinTools.createProteinSequence(seqString, oldSequence.getName() + " " + sequence.getName());
                } catch (IllegalSymbolException e) {
                    //this should never happen, since it should have been caught much earlier on
                    e.printStackTrace();
                }
            }
            uniqueSequenceStrings.put(seqString, oldResultObject); //replaces if already present
        }
        //now put all unique hits back into result set
        Set uniqueSequences = uniqueSequenceStrings.keySet();
        for (Iterator iterator = uniqueSequences.iterator(); iterator.hasNext();) {
            String seqString = (String) iterator.next();
            Object[] oldResultObject = (Object[]) uniqueSequenceStrings.get(seqString);
            //directly insert into private fields within this class
            uniqueResults.results.put((Sequence) oldResultObject[0], (ArrayList) oldResultObject[1]);
        }
        return uniqueResults;
    }

    /**
     * Keeps track if we've seen this sequence before
     *
     * @param sequenceTracker The HashMap used to do the tracking
     * @param sequence        The sequence to test
     * @return true if the sequence has been seen before
     */
    private boolean sequenceSeenBefore(HashMap sequenceTracker, Sequence sequence) {
        if (!sequenceTracker.containsKey(sequence)) {
            //new sequence
            sequenceTracker.put(sequence, sequence);
            return false;
        }
        //we've already seen this sequence
        return true;
    }

    /**
     * Gets the top 'percentage' percent number of results (using the number of database
     * sequences stored in this result set). This only makes sense for profile search results
     * which have an associated score to threshold by.
     *
     * @param percentage The percentage to get. Specified from 0.0 to 100.0
     *                   Any percentage under 0.0 will be interpreted at as 0.0
     *                   Any percentage over 100.0 will be interpreted at as 100.0
     * @return A new result set containing only the top percentage sequence results, null if
     *         no sequences are found (or if search was not done with a profile)
     */
    public SequenceSearchResultSet getTopPercentileResults(double percentage) {
        if (profile == null) {
            return null;
        }
        if (percentage <= 0.0) {
            percentage = 0.0;
            return null;
        }
        if (percentage >= 100.0) {
            return this;
        }
        SequenceSearchResultSet newResultSet = new SequenceSearchResultSet(scoreThreshold, params);
        newResultSet.setProfile(profile);
        newResultSet.setSpecies(species);
        newResultSet.setNumberOfSequencesSearched(this.numberOfSequencesSearched);
        //figure out how many sequences to get
        long numberSequences = (long) (Math.ceil(this.getNumberSequencesHit() * (percentage / 100.0)));
        //create a map of sequences to best score per sequence
        HashMap sequenceTracker = new HashMap();
        //get the sequences
        Set scores = sortedResultMap.keySet();
        int i = 0;
        for (Iterator iterator = scores.iterator(); iterator.hasNext();) {
            Double score = (Double) iterator.next();
            ArrayList hits = (ArrayList) sortedResultMap.get(score);
            for (int j = 0; j < hits.size(); j++) {
                Hit hit = (Hit) hits.get(j);
                newResultSet.addSequenceHit(hit.getSequence(), hit.getStart(), hit.getEnd(), hit.getScore());
                newResultSet.addOriginalSequenceToHit(hit.getSequence(), this.getOriginalSequence(hit.getSequence()));
                if (!sequenceSeenBefore(sequenceTracker, hit.getSequence())) {
                    //keep track of how many sequences we've seen
                    i++;
                }
            }
            //break if we've hit the number of sequences we want to keep after getting through all
            //hits for a given score. It is possible that more sequences will be returned than we want,
            //but we stop after getting through a given score.
            if (i >= numberSequences) {
                break;
            }
        }
        return (newResultSet);
    }

    /**
     * Gets the number of top results specified from this result set
     *
     * @param numberResults The number of results to return
     * @return The number of top results to return. The actual number of results may be more than the specified
     *         number because many hits might have the same score.
     */
    public SequenceSearchResultSet getTopResults(int numberResults) {
        if (profile == null) {
            return null;
        }
        if (numberResults <= 0.0) {
            return null;
        }
        if (numberResults >= this.getNumberOfHits()) {
            return this;
        }
        //convert number of results to a percentage
        double percentage = ((double) numberResults / (double) this.getNumberSequencesHit()) * 100;
        return (this.getTopPercentileResults(percentage));
    }

    /**
     * Gets the best Hit in the result set. There are potentially many Hits with the same score
     * so this method returns a List
     *
     * @return List of best scoring Hit objects, null if no hits
     * @see Hit
     */
    public List getBestHits() {
        //only makes sense for profile hits
        if (profile == null) {
            return null;
        }
        //sortedResultMap contains scores sorted in descending order
        ArrayList hits = null;
        Set resultSet = sortedResultMap.keySet();
        for (Iterator iterator = resultSet.iterator(); iterator.hasNext();) {
            Double score = (Double) iterator.next();
            hits = (ArrayList) sortedResultMap.get(score);
            break; //only interested in the best scoring hits
        }
        return hits;
    }

    /**
     * Merge two SequenceSearchResultSet objects. Merged results are in the object method was called from.
     *
     * @param resultSetToMerge the SequenceSearchResultSet to merge with the one that the method is called from.
     */
    public void mergeResultsFrom(SequenceSearchResultSet resultSetToMerge) {
        Set toMerge = resultSetToMerge.getSequences();
        for (Iterator iterator = toMerge.iterator(); iterator.hasNext();) {
            Sequence sequence = (Sequence) iterator.next();
            List hits = resultSetToMerge.getHits(sequence);
            for (int i = 0; i < hits.size(); i++) {
                Hit hit = (Hit) hits.get(i);
                this.addSequenceHit(sequence, hit.getStart(), hit.getEnd(), hit.getScore());
            }
        }
    }

    /**
     * Cleans up a result set - should call if you are generating a lot of results to prevent memory bloat
     */
    public void clear() {
        results.clear();
        hashResults.clear();
        if (annotatedSequenceMap != null) {
            annotatedSequenceMap.clear();
        }
        sortedResultMap.clear();
    }

    /**
     * Get a histogram of the scores in this result set. E.g. bin 1 contains the number of hits with 0 < score <= 1
     * bin 2 contains the number of hits with 1 < score <= 2
     *
     * @param highestScoreThreshold Calculate the histogram up until this value
     * @return An array corresponding to the histogram from 0 to highestScoreThreshold (size highestScoreThreshold+1)
     */
    public int[] getScoreHistogram(int highestScoreThreshold) {
        //only makes sense for profiles
        if (profile == null) {
            return null;
        }
        int histogram[] = new int[highestScoreThreshold + 1];
        ArrayList hits = null;
        Set resultSet = sortedResultMap.keySet();
        for (Iterator iterator = resultSet.iterator(); iterator.hasNext();) {
            Double score = (Double) iterator.next();
            hits = (ArrayList) sortedResultMap.get(score);
            histogram[(int) Math.ceil(score.doubleValue())] += hits.size();
        }
        return histogram;
    }

    /**
     * Returns a nice string representation of this result set
     * The output is sorted by score if available (profile search results)
     * Note: output might be large
     */
    public String toString() {
        if (results == null) {
            return null;
        }
        String lineSep = System.getProperty("line.separator");
        if (results.size() == 0) {
            if (profile != null) {
                return new String("No results for " + profile.getName() + lineSep);
            } else {
                return new String("No results for " + regex + lineSep);
            }
        }
        StringBuffer sb = new StringBuffer();
        Set resultSet = null;
        Map generalResultMap = null;
        //output header
        if (profile != null) {
            //profile type search
            sb.append("Profile description: " + profile.getName() + lineSep);
            sb.append("Peptides in profile: " + profile.getNumSequences() + lineSep);
            sb.append("Sequences searched: " + this.getNumberOfSequencesSearched() + lineSep);
            sb.append("Sequences hit: " + this.getNumberSequencesHit() + lineSep);
            sb.append("Hits: " + this.getNumberOfHits() + lineSep);
            sb.append("MatchPosition\tScore\tName\tMatch\tSequence" + lineSep);
            generalResultMap = sortedResultMap;
        } else if (regex != null) {
            //regex type search
            sb.append("Pattern: " + regex + lineSep);
            sb.append("Number sequences searched: " + this.getNumberOfSequencesSearched() + lineSep);
            sb.append("Number sequences hit: " + this.getNumberSequencesHit() + lineSep);
            sb.append("Number hits: " + this.getNumberOfHits() + lineSep);
            sb.append("MatchPosition\tName\tMatch\tSequence" + lineSep);
            generalResultMap = results;
        } else {
            //this should not happen, since we're maintaining these variables internally
            throw new RuntimeException("SequenceSearchResultSet does not contain a regex or profile.");
        }
        ArrayList hits = null;
        resultSet = generalResultMap.keySet();
        for (Iterator iterator = resultSet.iterator(); iterator.hasNext();) {
            Object keyObject = iterator.next();
            hits = (ArrayList) generalResultMap.get(keyObject);
            //get sequence pattern matches
            for (int i = 0; i < hits.size(); i++) {
                Hit hit = (Hit) hits.get(i);
                sb.append(hit.getStart() + "-" + hit.getEnd());
                if (profile != null) {
                    sb.append("\t" + hit.getScore());
                }
                sb.append("\t" + hit.getSequence().getName() + "\t" + hit.getMatchString() +
                        "\t" + hit.getSequence().seqString() + lineSep);
            }
        }
        return (sb.toString());
    }
}
