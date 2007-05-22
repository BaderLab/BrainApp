package org.baderlab.brain;

import cytoscape.task.TaskMonitor;
import org.biojava.bio.BioException;
import org.biojava.bio.dp.ScoreType;
import org.biojava.bio.dp.WeightMatrixAnnotator;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceAnnotator;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.impl.SimpleSequence;
import org.biojava.bio.seq.io.SeqIOTools;
import org.biojava.bio.symbol.Location;
import org.biojava.utils.ChangeVetoException;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

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
 * * Time: 4:38:09 PM
 * * Description Runs a sequence search against a protein database using a profile or regular expression
 */

/**
 * Implements profile and regular expression searching of a protein sequence flatfile database (in many formats)
 */
public class ProteinDatabaseSearch implements SequenceIdentifierAware {
    private SequenceIterator searchDatabase = null;
    private BufferedReader br = null;
    private String dbFileName = null;
    private String dbFormat = null;
    //stats
    private long lastSearchTime;
    //task monitoring
    private boolean cancelled = false;
    private TaskMonitor taskMonitor = null;

    /**
     * Creates a new protein database search object (opens the database)
     *
     * @param dbFileName The filename of the protein database to search
     * @param dbFormat   The record format of the protein database. Uses BioJava format strings, which are:
     *                   fasta, embl, swissprot, genpept
     * @throws FileNotFoundException Your sequence database was not found.
     * @throws BioException          Your database format was not recognized, the format choice does not match the actual file
     *                               format or the database does not contain protein sequences in the chosen format.
     */
    public ProteinDatabaseSearch(String dbFileName, String dbFormat) throws FileNotFoundException, BioException {
        br = new BufferedReader(new FileReader(dbFileName));
        searchDatabase = (SequenceIterator) SeqIOTools.fileToBiojava(dbFormat, "PROTEIN", br);
        this.dbFileName = dbFileName;
        this.dbFormat = dbFormat;
    }

    /**
     * Gets the last search time in milliseconds
     */
    public long getLastSearchTime() {
        return lastSearchTime;
    }

    /**
     * Sets the task monitor to update (if search is done within a Task)
     *
     * @see cytoscape.task.Task
     */
    public void setTaskMonitor(TaskMonitor taskMonitor) {
        this.taskMonitor = taskMonitor;
    }

    /**
     * If set, will schedule the algorithm to cancelled at the next convenient opportunity
     *
     * @param cancelled Set to true if the algorithm should be cancelled
     */
    public void setCancelled(boolean cancelled) {
        this.cancelled = cancelled;
    }

    /**
     * Runs a profile sequence search against the sequence database set in the constructor
     *
     * @param profile        The protein profile to search using
     * @param scoreThreshold The probability score threshold
     *                       Probability of a match given the distribution of the weight matrix
     *                       E.g. 10 means only return matches to the profile higher than p=1E-10
     *                       The lower this number, the less matches will result
     *                       This number must be 0 or higher where 0 (p=1.0) is the most stringent threshold
     *                       where the score of the match must be equal to p=1.0 (perfect match).
     * @param params         General search parameters
     * @return The results as a sequence search result set (contains parameters and sequence hits)
     * @throws BioException if search database contains errors
     */
    public SequenceSearchResultSet profileSearchDB(ProteinProfile profile, double scoreThreshold, ProteinDatabaseSearchParams params) throws BioException {
        ArrayList profileList = new ArrayList(1);
        profileList.add(profile);
        ArrayList scoreThresholdList = new ArrayList(1);
        scoreThresholdList.add(new Double(scoreThreshold));
        MultiSequenceSearchResultSet resultSets = multiProfileSearchDB(profileList, scoreThresholdList, params);
        return (resultSets.getResultSet(profile));
    }

    /**
     * Runs a profile sequence search against the sequence database set in the constructor
     *
     * @param profileList        The list of protein profiles to search using
     * @param scoreThresholdList A list of score thresholds corresponding to the protein profile list
     *                           Each protein profile can have its own score threshold.
     *                           Negative log of probability score threshold
     *                           Probability of a match given the distribution of the weight matrix
     *                           E.g. 10 means only return matches to the profile higher than p=1E-10
     *                           The lower this number, the less matches will result
     *                           This number must be 0 or higher where 0 (p=1.0) is the most stringent threshold
     *                           where the score of the match must be equal to p=1.0 (perfect match).
     * @param params             General search parameters
     * @return The results as a sequence search result set (contains parameters and sequence hits)
     * @throws BioException if search database contains errors
     */
    public MultiSequenceSearchResultSet multiProfileSearchDB(List profileList, List scoreThresholdList, ProteinDatabaseSearchParams params) throws BioException {
        Sequence sequenceFromDB;
        SequenceAnnotator annotator;
        SequenceSearchResultSet resultSet;
        long count = 0;
        HashMap profile2annotator; //key=profile, value=annotator (either wma or npa)

        //for time stats
        long msTimeBefore = System.currentTimeMillis();

        //initialize the multi set and weight matrix annotators
        MultiSequenceSearchResultSet multiResultSet = new MultiSequenceSearchResultSet();
        profile2annotator = new HashMap();
        for (int i = 0; i < profileList.size(); i++) {
            ProteinProfile profile = (ProteinProfile) profileList.get(i);
            //set up the result set for this profile
            //convert from negative log likelihood to p-value
            double scoreThreshold = ((Double) scoreThresholdList.get(i)).doubleValue();
            double pValueScoreThreshold = Math.pow(10, -scoreThreshold);
            resultSet = new SequenceSearchResultSet(pValueScoreThreshold, params);
            resultSet.setProfile(profile);
            multiResultSet.add(resultSet);
            //set up the annotator for this profile
            if (params.getScoreType() == ProteinDatabaseSearchParams.SCORE_TYPE_NORM_PROBABILITY) {
                annotator = new NormalizedProfileAnnotator(profile, pValueScoreThreshold);
            } else if (params.getScoreType() == ProteinDatabaseSearchParams.SCORE_TYPE_PROBABILITY) {
                annotator = new WeightMatrixAnnotator(profile.getWeightMatrix(), ScoreType.PROBABILITY, pValueScoreThreshold);
            } else if (params.getScoreType() == ProteinDatabaseSearchParams.SCORE_TYPE_DELTAG) {
                annotator = new DeltaGProfileAnnotator(profile, scoreThreshold);
            } else {
                throw new IllegalArgumentException("Illegal score type passed. (" + params.getScoreType() + ")");
            }
            profile2annotator.put(profile, annotator);
            if (cancelled) {
                break;
            }
        }

        //Do the actual search
        while (searchDatabase.hasNext()) {
            if (cancelled) {
                break;
            }
            sequenceFromDB = searchDatabase.nextSequence();
            Sequence filteredSequenceFromDB = filterTrailingAsterisk(sequenceFromDB);
            //optionally filter the sequence
            Sequence sequenceToSearch = null;
            if (params.getLength() > 0) {
                //if length is set, use it, otherwise get length from profile later
                sequenceToSearch = ProteinTerminus.getSequenceTerminus(filteredSequenceFromDB, params.getLength(), params.getTerminus());
            }

            //iterate through the profiles to search the sequence
            for (int i = 0; i < profileList.size(); i++) {
                if (cancelled) {
                    break;
                }
                ProteinProfile profile = (ProteinProfile) profileList.get(i);
                resultSet = multiResultSet.getResultSet(profile);
                //check if the sequence is ready
                if (sequenceToSearch == null) {
                    //case of single hit from protein terminus, get the sequence length based on the profile length
                    //TODO: detect if profile length has changed and adjust the sequence to search
                    //TODO: if we do this, the profiles should be sorted by length, longest first, as an optimization
                    //and the old sequence should be truncated further at each shorter profile
                    sequenceToSearch = ProteinTerminus.getSequenceTerminus(filteredSequenceFromDB, profile.getNumColumns(), params.getTerminus());
                }
                //Use the weight matrix to search the sequence
                annotator = (SequenceAnnotator) profile2annotator.get(profile);
                Sequence annotatedSequenceToSearch;
                try {
                    annotatedSequenceToSearch = annotator.annotate(sequenceToSearch);
                } catch (ChangeVetoException e) {
                    //this shouldn't happen, since we're creating our own sequences internally
                    e.printStackTrace();
                    return null;
                }
                //if it was annotated, add it to the result set (should be one feature per pattern match)
                if (annotatedSequenceToSearch.countFeatures() > 0) {
                    //check if match contains X or - (these characters interfere with scoring - they
                    //generate very high scores because the Distribution class returns a weight of 1.0 for X)
                    boolean containsInvalidSymbols = false;
                    if ((annotatedSequenceToSearch.seqString().indexOf("X") >= 0) || (annotatedSequenceToSearch.seqString().indexOf("-") >= 0)) {
                        containsInvalidSymbols = true;
                    }
                    //add the hit to the result set
                    for (Iterator it = annotatedSequenceToSearch.features(); it.hasNext();) {
                        //find the score annotation in the feature
                        Feature f = (Feature) it.next();
                        if ((f.getAnnotation() != null) && (f.getAnnotation().containsProperty("score"))) {
                            Location loc = f.getLocation();
                            Double score = (Double) f.getAnnotation().getProperty("score");
                            double convertedScore = 0.0;
                            if ((params.getScoreType() == ProteinDatabaseSearchParams.SCORE_TYPE_NORM_PROBABILITY) ||
                                    (params.getScoreType() == ProteinDatabaseSearchParams.SCORE_TYPE_PROBABILITY)) {
                                //convert to a negative log of the p-value for better display (and to match the threshold parameter)
                                double neglogscore;
                                if (score.doubleValue() == 0.0) {
                                    neglogscore = Double.MAX_VALUE;
                                } else if (score.doubleValue() == 1.0) {
                                    neglogscore = 0.0;
                                } else {
                                    neglogscore = -Math.log(score.doubleValue()) / Math.log(10);
                                }
                                //check if profile match contains X or - residues, if so give it the min possible score
                                if (containsInvalidSymbols) {
                                    neglogscore = -Math.log(resultSet.getScoreThreshold()) / Math.log(10);
                                }
                                convertedScore = neglogscore;
                            } else if (params.getScoreType() == ProteinDatabaseSearchParams.SCORE_TYPE_DELTAG) {
                                convertedScore = score.doubleValue();
                            }
                            //Note: copy the sequence so it keeps its features without interfering with the annotation
                            //of the features on the next loop, making sure to keep all annotation and features
                            //Note: Sequence coordinates start at 1
                            if (params.isDontSaveSequences()) {
                                resultSet.addSequenceHit(new Double(convertedScore));
                            } else {
                                resultSet.addSequenceHit(sequenceToSearch, loc.getMin(), loc.getMax(), new Double(convertedScore));
                                resultSet.addOriginalSequenceToHit(sequenceToSearch, sequenceFromDB);
                            }
                            if (!params.isMultipleHits()) {
                                break; //don't continue if we're only interested in the first hit
                            }
                        }
                    }
                }
            }
            count++;
            if ((taskMonitor != null) && (count % 500 == 0)) {
                taskMonitor.setStatus("Searching sequence " + count);
            }
        }

        //update all of the result sets with the final sequence count
        Collection resultSets = multiResultSet.getAllResultSets();
        for (Iterator iterator = resultSets.iterator(); iterator.hasNext();) {
            resultSet = (SequenceSearchResultSet) iterator.next();
            resultSet.setNumberOfSequencesSearched(count);
        }

        //update time stats
        long msTimeAfter = System.currentTimeMillis();
        lastSearchTime = msTimeAfter - msTimeBefore;

        //update task monitor
        if (taskMonitor != null) {
            taskMonitor.setStatus("Searching sequence " + count);
        }

        return (multiResultSet);
    }

    /**
     * Runs a profile sequence search against the sequence database set in the constructor
     *
     * @param regex  The protein profile to search using
     * @param params General search parameters
     * @return The results as a sequence search result set (contains parameters and sequence hits)
     * @throws BioException if search database contains errors
     */
    public SequenceSearchResultSet regexSearchDB(String regex, ProteinDatabaseSearchParams params) throws BioException {
        ArrayList regexList = new ArrayList(1);
        regexList.add(regex);
        MultiSequenceSearchResultSet resultSets = multiRegexSearchDB(regexList, params);
        return (resultSets.getResultSet(regex));
    }

    /**
     * Runs a regex sequence search against the sequence database set in the constructor
     *
     * @param regexList The list of protein regex patterns to search using
     * @param params    General search parameters
     * @return The results as a sequence search result set (contains parameters and sequence hits)
     * @throws BioException if search database contains errors
     */
    public MultiSequenceSearchResultSet multiRegexSearchDB(List regexList, ProteinDatabaseSearchParams params) throws BioException {
        Sequence sequenceFromDB;
        SequenceSearchResultSet resultSet;
        long count = 0;

        //for time stats
        long msTimeBefore = System.currentTimeMillis();

        //initialize the multi set
        MultiSequenceSearchResultSet multiResultSet = new MultiSequenceSearchResultSet();
        for (int i = 0; i < regexList.size(); i++) {
            String regex = (String) regexList.get(i);
            resultSet = new SequenceSearchResultSet(params);
            resultSet.setRegex(regex);
            multiResultSet.add(resultSet);
        }

        while (searchDatabase.hasNext()) {
            if (cancelled) {
                break;
            }
            sequenceFromDB = searchDatabase.nextSequence();
            sequenceFromDB = filterTrailingAsterisk(sequenceFromDB);
            //optionally filter the sequence
            Sequence sequenceToSearch = null;
            if (params.getLength() > 0) {
                //if length is set, use it, otherwise get length from regex later
                sequenceToSearch = ProteinTerminus.getSequenceTerminus(sequenceFromDB, params.getLength(), params.getTerminus());
            }

            //iterate through the regex's to search the sequence
            for (int i = 0; i < regexList.size(); i++) {
                if (cancelled) {
                    break;
                }
                String regex = (String) regexList.get(i);
                resultSet = multiResultSet.getResultSet(regex);
                //check if the sequence is ready
                if (sequenceToSearch == null) {
                    sequenceToSearch = ProteinTerminus.getSequenceTerminus(sequenceFromDB, regex.length(), params.getTerminus());
                }
                Pattern p = Pattern.compile(regex);
                Matcher m = p.matcher(sequenceToSearch.seqString());
                //search for potentially overlapping hits (which is why we have to reset find each iteration)
                int start = 0;
                while (m.find(start)) {
                    start = m.start();
                    //add the hit to the result set
                    //Note: Sequence coordinates start at 1 (must translate from normal string coordinates here
                    resultSet.addSequenceHit(sequenceToSearch, start + 1, m.end());
                    start++; //move ahead to search for more hits
                    if (!params.isMultipleHits()) {
                        break; //don't continue if we're only interested in the first hit
                    }
                }
            }
            count++;
            if ((taskMonitor != null) && (count % 500 == 0)) {
                taskMonitor.setStatus("Searching sequence " + count);
            }
        }

        //update all of the result sets with the final sequence count
        Collection resultSets = multiResultSet.getAllResultSets();
        for (Iterator iterator = resultSets.iterator(); iterator.hasNext();) {
            resultSet = (SequenceSearchResultSet) iterator.next();
            resultSet.setNumberOfSequencesSearched(count);
        }

        //update time stats
        long msTimeAfter = System.currentTimeMillis();
        lastSearchTime = msTimeAfter - msTimeBefore;

        //update task monitor
        if (taskMonitor != null) {
            taskMonitor.setStatus("Searching sequence " + count);
        }

        return (multiResultSet);
    }

    /**
     * Closes the database. This should be called after searching to prevent wasting file handles.
     *
     * @throws IOException if the database can't be closed.
     */
    public void close() throws IOException {
        br.close();
    }

    /**
     * Reset a search so it can be searched again from the start. You must call this before re-searching a database.
     *
     * @throws IOException  If there is an IO problem with the file (e.g. it can't be closed or opened)
     * @throws BioException If there is a problem reading the file format which shouldn't happen unless
     *                      file was changed during the previous search
     */
    public void reset() throws IOException, BioException {
        br.close();
        br = new BufferedReader(new FileReader(dbFileName));
        searchDatabase = (SequenceIterator) SeqIOTools.fileToBiojava(dbFormat, "PROTEIN", br);
    }

    /**
     * Gets the default BioJava identifier from a Sequence. User must override to get the identifier of choice.
     * E.g. for FASTA file, this is normally the first token in the description line
     *
     * @param sequence The sequence containing the identifier
     * @return the identifier
     */
    public String getIdentifier(Sequence sequence) {
        return sequence.getName();
    }

    /**
     * Some protein sequence databases append a trailing asterisk '*' to the end of each sequence.
     * This needs to be removed here, because it can interfere with C-terminal searches.
     */
    private Sequence filterTrailingAsterisk(Sequence sequence) {
        if (!sequence.subStr(sequence.length(), sequence.length()).equals("*")) {
            return sequence;
        }
        //else sequence contains a * as the last character
        //TODO: just add a sequence feature here and check it before searching(?)
        SimpleSequence sequenceSubSet = new SimpleSequence(sequence.subList(1, sequence.length() - 1), "", sequence.getName(), null);
        return (sequenceSubSet);
    }
}
