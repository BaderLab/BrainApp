package org.baderlab.csplugins.brain;

import cytoscape.CyEdge;
import cytoscape.CyNetwork;
import cytoscape.CyNode;
import cytoscape.Cytoscape;
import org.biojava.bio.Annotation;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.Sequence;
import org.mskcc.dataservices.bio.vocab.CommonVocab;
import org.mskcc.dataservices.bio.vocab.InteractionVocab;
import org.mskcc.dataservices.bio.vocab.InteractorVocab;

import java.util.*;

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
 * * Date: Mar 28, 2005
 * * Time: 5:57:09 PM
 * * Description Utility functions to convert protein database search results to Cytoscape
 */

// mdharsee 20070212: CyNet.setEdgeAttributeValue is deprecated as of v2.2 - using cyAttributes instead

/**
 * Utility functions to convert protein database search results to Cytoscape
 */
public class CytoscapeUtil {
    static int runNumber = 1;

    /**
     * Add a result set to Cytoscape
     *
     * @param net       The Cytoscape network to add to
     * @param profile   The name (reference) of the profile to link to all members of the result set
     * @param resultSet The result set to add
     */
    public static void addSequenceSearchResultSetToCytoscape(CyNetwork net, ProteinProfile profile,
                                                             SequenceSearchResultSet resultSet, BrainParameterSet params) {
        //create the node to connect to all search results
        //(e.g. the protein that contains the domain that has the binding profile used for the search)
        CyNode nodeA = Cytoscape.getCyNode(profile.getName());
        if (nodeA == null || net.getIndex(nodeA) == 0) {
            nodeA = Cytoscape.getCyNode(profile.getName(), true);
            net.addNode(nodeA);
        }
        addPSIFeaturesToNode(net, nodeA, profile);
        Set results = resultSet.getSequences();
        String nodeBName = null;
        for (Iterator iterator = results.iterator(); iterator.hasNext();) {
            Sequence sequence = (Sequence) iterator.next();
            Sequence originalDBSequence = resultSet.getOriginalSequence(sequence);
            //NODE
            //add the edge to the network
            if (params.getDatabaseFormat().equalsIgnoreCase("genpept")) {
                nodeBName = findNameInGenPeptSequence(originalDBSequence);
                if (nodeBName == null) {
                    nodeBName = originalDBSequence.getName();
                } else {
                    nodeBName = nodeBName.concat("_" + originalDBSequence.getName());
                }
            } else {
                nodeBName = originalDBSequence.getName();
            }
            CyNode nodeB = Cytoscape.getCyNode(nodeBName);
            if (nodeB == null || net.getIndex(nodeB) == 0) {
                nodeB = Cytoscape.getCyNode(nodeBName, true);
                net.addNode(nodeB);
            }
            //EDGE
            CyEdge edge = Cytoscape.getCyEdge(nodeA.getIdentifier(), nodeA.getIdentifier() +
                    "_pp_" + nodeB.getIdentifier(), nodeB.getIdentifier(), "pp" +
                    /*to force uniqueness of edges so scores can be unique to each network*/
                    net.getIdentifier());
            net.addEdge(edge);
            //now add all hits as attributes
            double bestScore = Double.MAX_VALUE;
            String bestMotif = null;
            Hit bestHit = null;
            List hits = resultSet.getHits(sequence);
            for (int i = 0; i < hits.size(); i++) {
                Hit hit = (Hit) hits.get(i);
                //keep track of best score
                if (hit.getScore().doubleValue() < bestScore) {
                    bestHit = hit;
                    bestScore = hit.getScore().doubleValue();
                    bestMotif = hit.getMatchString();
                }
                //net.setNodeAttributeValue(nodeB, "Motif Hit", bestMotif);
                if (bestMotif != null) {
                    Cytoscape.getNodeAttributes().setAttribute(nodeB.getIdentifier(), "Motif Hit", bestMotif);
                    Cytoscape.getNodeAttributes().setAttribute(nodeB.getIdentifier(), "Motif Start", bestHit.getStart());
                    Cytoscape.getNodeAttributes().setAttribute(nodeB.getIdentifier(), "Motif End", bestHit.getEnd());
                }
            }
            //add best score as an edge attribute
            bestScore = truncateDouble(bestScore, 3);
            //net.setEdgeAttributeValue(edge, "HighestScore", new Double(bestScore));
            Cytoscape.getEdgeAttributes().setAttribute(edge.getIdentifier(), "HighestScore", new Double(bestScore));
            if (params.getDatabaseFormat().equalsIgnoreCase("genpept")) {
                addPSIFeaturesToEdge(net, edge, bestScore, bestMotif, nodeA.getIdentifier());
            }

            //add annotation to the node, based on what database format it came from
            if (params.getDatabaseFormat().equalsIgnoreCase("fasta")) {
                addFastaFeaturesToNode(net, nodeB, originalDBSequence);
            } else if (params.getDatabaseFormat().equalsIgnoreCase("genpept")) {
                addGenPeptFeaturesToNode(net, nodeB, originalDBSequence, params, sequence, bestHit);
            }

        }
    }

    /**
     * Top-level method to run a profile search and view the results in Cytoscape
     *
     * @param searchResults The search results to add to a Cytoscape network
     * @return The network the search results were added to
     */
    public static CyNetwork addProfileSearchResultsToCytoscape(MultiSequenceSearchResultSet searchResults, BrainParameterSet params) {
        CyNetwork net = Cytoscape.createNetwork("Profile Network " + runNumber++);
        if (searchResults == null) {
            System.err.println("Search error. Can't continue.");
            return null;
        }
        Collection results = searchResults.getAllResultSets();
        for (Iterator iterator = results.iterator(); iterator.hasNext();) {
            SequenceSearchResultSet sequenceSearchResultSet = (SequenceSearchResultSet) iterator.next();
            //evaluate top percentage threshold if set
            if (params.getScorePercentageThreshold() < 100.0) {
                sequenceSearchResultSet = sequenceSearchResultSet.getTopPercentileResults(params.getScorePercentageThreshold());
            }
            //evaluate number of top hits threshold if set
            if (params.getNumberTopHits() > 0) {
                sequenceSearchResultSet = sequenceSearchResultSet.getTopResults(params.getNumberTopHits());
            }
            System.out.println(sequenceSearchResultSet.getProfile().getName() + "\t" +
                    sequenceSearchResultSet.getNumberSequencesHit() + System.getProperty("line.separator"));
            CytoscapeUtil.addSequenceSearchResultSetToCytoscape(net, sequenceSearchResultSet.getProfile(), sequenceSearchResultSet, params);
            System.out.println(sequenceSearchResultSet.getNumberOfHits());
        }
        return net;
    }

    /**
     * Add Cytoscape attributes for saving to PSI-MI to a node, given a DatabaseReference object
     */
    private static void addPSIFeaturesToNode(CyNetwork net, CyNode node, ProteinProfile profile) {
        DatabaseReference protein = profile.getProteinReference();
        String[] dbName = new String[1];
        dbName[0] = protein.getDbname();
        //net.setNodeAttributeValue(node, CommonVocab.XREF_DB_NAME, dbName);
        if (dbName[0] != null) {
            Cytoscape.getNodeAttributes().setAttribute(node.getIdentifier(), (String)CommonVocab.XREF_DB_NAME, dbName[0]);
        }

        String[] dbID = new String[1];
        dbID[0] = protein.getDbid();
        //net.setNodeAttributeValue(node, CommonVocab.XREF_DB_ID, dbID);
        if (dbID[0] != null) {
            Cytoscape.getNodeAttributes().setAttribute(node.getIdentifier(), CommonVocab.XREF_DB_ID, dbID[0]);
        }

        //add profile features to node
        //net.setNodeAttributeValue(node, "Domain Number", new Integer(profile.getDomainNumber()));
        Cytoscape.getNodeAttributes().setAttribute(node.getIdentifier(), "Domain Number", new Integer(profile.getDomainNumber()));
        //net.setNodeAttributeValue(node, "Domain Sequence", profile.getDomainSequence());
        if (profile.getDomainSequence() != null) {
            Cytoscape.getNodeAttributes().setAttribute(node.getIdentifier(), "Domain Sequence", profile.getDomainSequence());
        }
        //net.setNodeAttributeValue(node, "Domain Type", profile.getDomainName());
        if (profile.getDomainName() != null) {
            Cytoscape.getNodeAttributes().setAttribute(node.getIdentifier(), "Domain Type", profile.getDomainName());
        }
        //net.setNodeAttributeValue(node, "Experimental Method", profile.getExperimentalMethod());
        if (profile.getExperimentalMethod() != null) {
            Cytoscape.getNodeAttributes().setAttribute(node.getIdentifier(), "Experimental Method", profile.getExperimentalMethod());
        }
        //net.setNodeAttributeValue(node, "ProfileNumSequences", new Integer(profile.getNumSequences()));
        Cytoscape.getNodeAttributes().setAttribute(node.getIdentifier(), "ProfileNumSequences", new Integer(profile.getNumSequences()));

        if (profile.getDomainSequenceStart() != 0) {
            Cytoscape.getNodeAttributes().setAttribute(node.getIdentifier(), "Sequence Start", profile.getDomainSequenceStart());
        }
        if (profile.getDomainSequenceStop() != 0) {
            Cytoscape.getNodeAttributes().setAttribute(node.getIdentifier(), "Sequence Stop", profile.getDomainSequenceStop());
        }
        if (profile.getComment() != null) {
            Cytoscape.getNodeAttributes().setAttribute(node.getIdentifier(), "Comment", profile.getComment());
        }
    }

    /**
     * Add features from a FASTA file to a node (these are very basic)
     *
     * @param net      The network containing the attributes
     * @param node     The node to be annotated
     * @param sequence The sequence containing the FASTA derived annotations
     */
    private static void addFastaFeaturesToNode(CyNetwork net, CyNode node, Sequence sequence) {
        String description = null;
        //retrieve the sequence description
        Annotation seqAnn = sequence.getAnnotation();
        if (seqAnn.containsProperty("description")) {
            description = (String) seqAnn.getProperty("description");
        } else {
            //if no description, just use the name
            description = sequence.getName();
        }
        //net.setNodeAttributeValue(node, "Description", description);
        Cytoscape.getNodeAttributes().setAttribute(node.getIdentifier(), "Description", description);
    }

    /**
     * Extract the value of an annotation from a GenPept file
     *
     * @param name The name of the annotation
     * @param f    The feature containing the annotation
     * @return A string containing the value of the annotation
     */
    private static ArrayList extractGenericAnnotation(String name, Feature f) {
        String xref = null;
        ArrayList xrefList = null; //the annotation is an ArrayList when there are more than one
        if (f.getAnnotation().containsProperty(name)) {
            Object xrefObj = f.getAnnotation().getProperty(name);
            if (xrefObj.getClass().equals(String.class)) {
                xref = (String) xrefObj;
                xrefList = new ArrayList(1);
                xrefList.add(xref);
            } else if (xrefObj.getClass().equals(ArrayList.class)) {
                xrefList = (ArrayList) xrefObj;
            }
        }
        return xrefList;
    }

    /**
     * Extract a db_xref annotation from a GenPept file. This needs its own method, since there are
     * some exceptions to deal with.
     *
     * @param f The feature to extract the information from
     * @return A string containing the db_xref value
     */
    private static String extractTaxDbXref(Feature f) {
        String xref = null;
        ArrayList xrefList = extractGenericAnnotation("db_xref", f);
        if (xrefList.size() == 2) {
            xref = (String) xrefList.get(1);
        } else if (xrefList.size() == 1) {
            xref = (String) xrefList.get(0);
        } else {
            return null;
        }
        return xref;
    }

    /**
     * Extract the value of an annotation from a GenPept file
     *
     * @param name The name of the annotation
     * @param f    The feature containing the annotation
     * @return A string containing the value of the annotation
     */
    private static String extractAnnotationStringByName(String name, Feature f) {
        ArrayList retVal = extractGenericAnnotation(name, f);
        if (retVal == null) {
            return null;
        }
        return (String) retVal.get(0);
    }

    /**
     * extract annotation and features from a Sequence that originated from a GenPept flat file
     *
     * @param net                The Cytoscape network that contains the Cytoscape attributes object
     * @param node               The node to annotate
     * @param originalDBSequence The originalDBSequence containing the information
     * @param params             The parameters used for this search
     * @param searchedSequence   The actual sequence searched (might be a different length than the original DB sequence)
     * @param bestHit            The best scoring motif hit for the searchedSequence
     */
    private static void addGenPeptFeaturesToNode(CyNetwork net, CyNode node, Sequence originalDBSequence, BrainParameterSet params, Sequence searchedSequence, Hit bestHit) {
        Annotation seqAnn = originalDBSequence.getAnnotation();
        if (seqAnn.containsProperty("DEFINITION")) {
            Cytoscape.getNodeAttributes().setAttribute(node.getIdentifier(), "FULL_NAME", (String) seqAnn.getProperty("DEFINITION"));
        }
        if (seqAnn.containsProperty("ACCESSION")) {
            String accession = (String) seqAnn.getProperty("ACCESSION");
            String[] accessionArray = accession.split("\\s+");
            String[] firstAccession = new String[1];
            firstAccession[0] = accessionArray[0];
            if (accessionArray.length > 0) {
                Cytoscape.getNodeAttributes().setAttribute(node.getIdentifier(), CommonVocab.XREF_DB_ID, firstAccession[0]);
            }
        }
        if (seqAnn.containsProperty("DBSOURCE")) {
            String dbSource = (String) seqAnn.getProperty("DBSOURCE");
            //extract DB source, assume format "REFSEQ: accession NM_053517.1"
            int endSource = dbSource.indexOf(":");
            String[] dbArray = new String[1];
            dbArray[0] = dbSource.substring(0, endSource);
            if (dbArray[0] != null) {
                Cytoscape.getNodeAttributes().setAttribute(node.getIdentifier(), CommonVocab.XREF_DB_NAME, dbArray[0]);
            }
        }
        if (seqAnn.containsProperty("SOURCE")) {
            String source = (String) seqAnn.getProperty("SOURCE");
            //extract common name from source, assume format e.g. "Rattus norvegicus (Norway rat)"
            int startCommonName = source.indexOf("(");
            int endCommonName = source.indexOf(")");
            if ((startCommonName >= 0) || (endCommonName >= 0)) {
                Cytoscape.getNodeAttributes().setAttribute(node.getIdentifier(), InteractorVocab.ORGANISM_COMMON_NAME, source.substring(startCommonName + 1, endCommonName));
                Cytoscape.getNodeAttributes().setAttribute(node.getIdentifier(), InteractorVocab.ORGANISM_SPECIES_NAME, source.substring(0, startCommonName - 1));
            } else {
                if (source != null) {
                    //some records don't follow the above convention
                    Cytoscape.getNodeAttributes().setAttribute(node.getIdentifier(), InteractorVocab.ORGANISM_COMMON_NAME, source);
                    Cytoscape.getNodeAttributes().setAttribute(node.getIdentifier(), InteractorVocab.ORGANISM_SPECIES_NAME, source);
                }
            }
        }
        // loop over all features in a sequence
        int protLength = 0;
        TreeMap snpMap = null;
        for (Iterator fi = originalDBSequence.features(); fi.hasNext();) {
            Feature f = (Feature) fi.next();
            if ((f.getAnnotation() != null) && (f.getType().equals("source"))) {
                String xref = extractTaxDbXref(f);
                if (xref != null) {
                    int startTaxid = xref.indexOf(":");
                    if (startTaxid > 0) {
                        Cytoscape.getNodeAttributes().setAttribute(node.getIdentifier(), InteractorVocab.ORGANISM_NCBI_TAXONOMY_ID, xref.substring(startTaxid + 1, xref.length()));
                    }
                }
            }
            if ((f.getAnnotation() != null) && (f.getType().equals("CDS"))) {
                String name = extractAnnotationStringByName("gene", f);
                if (name != null) {
                    Cytoscape.getNodeAttributes().setAttribute(node.getIdentifier(), "GENE_NAME", name);
                }
                //extract MIM ref
                ArrayList dbRefs = extractGenericAnnotation("db_xref", f);
                if (dbRefs != null) {
                    for (int i = 0; i < dbRefs.size(); i++) {
                        String s = (String) dbRefs.get(i);
                        if ((s != null) && (s.startsWith("MIM"))) {
                            if (s.indexOf(":") > 0) {
                                Cytoscape.getNodeAttributes().setAttribute(node.getIdentifier(), "OMIM", s.substring(s.indexOf(":") + 1, s.length()));
                            }
                        }
                    }
                }
            }
            //extract SNP information (may be many features)
            if ((f.getAnnotation() != null) && (f.getType().equals("variation"))) {
                ArrayList snp = extractGenericAnnotation("replace", f);
                if (snp != null) {
                    if (snpMap == null) {
                        snpMap = new TreeMap();
                    }
                    snpMap.put(new Integer(f.getLocation().getMax()), snp);
                }
            }
            //extract SNP information (may be many features)
            if ((f.getAnnotation() != null) && (f.getType().equals("Protein"))) {
                protLength = f.getLocation().getMax();
            }
            /* Information from RefSeq is formatted like this:
            Protein         1..96
            /product="keratin associated protein 12-1"
            variation       34
            /replace="m"
            /replace="v"
            /db_xref="dbSNP:9984476"*/
        }
        if (snpMap != null) {
            //save the SNPs as a string representation
            StringBuffer sbAll = null;
            StringBuffer sbMotif = null;
            Set locations = snpMap.keySet();
            //find out how long searched sequence was
            int seqSearchLength = searchedSequence.length();
            //iterate through SNPs
            for (Iterator iterator = locations.iterator(); iterator.hasNext();) {
                Integer location = (Integer) iterator.next();
                ArrayList snp = (ArrayList) snpMap.get(location);
                //filter SNPs based on actual sequence searched
                boolean addAllSNP = false;
                boolean addMotifSNP = false;
                ProteinTerminus term = params.getSearchParams().getTerminus();
                if (term.equals(ProteinTerminus.C)) {
                    if (location.intValue() > (protLength - seqSearchLength)) {
                        addAllSNP = true;
                    }
                } else if (term.equals(ProteinTerminus.N)) {
                    if (location.intValue() < seqSearchLength) {
                        addAllSNP = true;
                    }
                } else if (term.equals(ProteinTerminus.NONE)) {
                    addAllSNP = true;
                }

                // if we got a best hit (which we always should) and the SNP is on the hit subsequence,
                // add the SNP to the Motif SNP list
                if (bestHit != null) {
                    if ((location.intValue() >= bestHit.getStart()) && (location.intValue() <= bestHit.getEnd()) ) {
                        addMotifSNP = true;
                    }
                }
                if (addAllSNP) {
                    if (sbAll == null) {
                        sbAll = new StringBuffer();
                    }
                    sbAll.append(((String) snp.get(0)).toUpperCase() + location.toString() + ((String) snp.get(1)).toUpperCase() + " ");
                }
                if (addMotifSNP) {
                    if (sbMotif == null) {
                        sbMotif = new StringBuffer();
                    }
                    sbMotif.append(((String) snp.get(0)).toUpperCase() + location.toString() + ((String) snp.get(1)).toUpperCase() + " ");
                }
            }
            if (sbAll != null) {
                sbAll.append("of " + protLength);
                Cytoscape.getNodeAttributes().setAttribute(node.getIdentifier(), "SNPs (All)", sbAll.toString());
            }
            if (sbMotif != null) {
                Cytoscape.getNodeAttributes().setAttribute(node.getIdentifier(), "SNPs (Motif Hit)", sbMotif.toString());
            }
        }
    }

    /**
     * Gets a sequence name (gene name) in the CDS feature in GenPept files
     *
     * @param sequence The sequence derived from a GenPept record
     * @return The CDS gene name, or null if not found;
     */
    private static String findNameInGenPeptSequence(Sequence sequence) {
        String name = null;
        // loop over all features in a sequence
        for (Iterator fi = sequence.features(); fi.hasNext();) {
            Feature f = (Feature) fi.next();
            if ((f.getAnnotation() != null) && (f.getType().equals("CDS"))) {
                name = extractAnnotationStringByName("gene", f);
            }
        }
        return (name);
    }

    /**
     * Add PSI-MI features to an edge, so the PSI-MI writer can save them
     *
     * @param net   The network containing the annotation object
     * @param edge  The edge to annotate
     * @param score The score attached to this edge
     * @param motif The motif attached to the score on this edge
     */
    private static void addPSIFeaturesToEdge(CyNetwork net, CyEdge edge, double score, String motif, String sourceName) {
        Cytoscape.getEdgeAttributes().setAttribute(edge.getIdentifier(), InteractionVocab.EXPERIMENTAL_SYSTEM_NAME, "Sequence-based Prediction (score:" + score + ";motif:" + motif + ";from:" + sourceName + ")");
        Cytoscape.getEdgeAttributes().setAttribute(edge.getIdentifier(), InteractionVocab.EXPERIMENTAL_SYSTEM_XREF_DB, "MI");
        Cytoscape.getEdgeAttributes().setAttribute(edge.getIdentifier(), InteractionVocab.EXPERIMENTAL_SYSTEM_XREF_ID, "MI:0101");
    }

    /**
     * Truncate a double-precision floating point number to a specific number of significant digits
     *
     * @return The new truncated double number (make sure to catch it)
     */
    private static double truncateDouble(double inputDouble, int numberSignificantDigits) {
        double returnDouble = inputDouble * Math.pow(10, numberSignificantDigits);
        returnDouble = Math.rint(returnDouble);
        return (returnDouble / Math.pow(10, numberSignificantDigits));
    }

    //TODO: figure out how to generically extract the name from the sequence description
    //ensGeneName = description.substring(description.indexOf("gene:")+5, description.indexOf("gene:") + 5 +15);
    //sgdGeneName = (description.split(" ", 3))[1];
    public static String getIPI_SPIdentifier(Sequence sequence) {
        String description = sequence.getName();
        //parse an IPI FASTA file description line, that looks like this, for the swiss-prot id:
        //>IPI:IPI00000001.1|SWISS-PROT:O95793-1|REFSEQ_NP:NP_059347|ENSEMBL:ENSP00000346163|H-INV:HIT000025197;HIT000014346 Tax_Id=9606
        //TODO: parse based on "SWISS-PROT: string" because order matters
        String splitArray[] = description.split("SWISS-PROT:|TREMBL:", 2);
        String proteinId = null;
        if (splitArray.length == 2) {
            proteinId = splitArray[1];
            proteinId = proteinId.split("\\||;|-", 2)[0];
            //proteinId = proteinId.split("-", 2)[0];
        } else {
            //assume == 1
            proteinId = splitArray[0];
        }
        return (proteinId);
    }

    /**
     * Check if a node contains a given database xref.  Uses CommonVocab to find the xref
     *
     * @return true if the node contains an xref
     */
    public static boolean doesNodeContainXref(CyNode node, DatabaseReference dbref) {
        //CyNetwork net = (CyNetwork) node.getGraphPerspective();

        //String[] dbNames = (String[]) net.getNodeAttributeValue(node, CommonVocab.XREF_DB_NAME);
        //String[] dbIDs = (String[]) net.getNodeAttributeValue(node, CommonVocab.XREF_DB_ID);

        List dbNames = Cytoscape.getNodeAttributes().getListAttribute(node.getIdentifier(), CommonVocab.XREF_DB_NAME);
        List dbIDs = Cytoscape.getNodeAttributes().getListAttribute(node.getIdentifier(), CommonVocab.XREF_DB_ID);


        for (int i = 0; i < dbNames.size(); i++) {
            String dbName = (String) dbNames.get(i);
            String dbID = (String) dbIDs.get(i);
            DatabaseReference nodeXref = new DatabaseReference(dbName, dbID);
            if (nodeXref.equals(dbref)) {
                return true;
            }
        }
        return false;
    }

    /**
     * Returns the database references that have been stored on a node using the CommonVocab framework
     *
     * @param net  The network containing the node and associated attributes
     * @param node The node
     * @return An array of DatabaseReference objects containing all database references found on the node
     */
    public static DatabaseReference[] getDatabaseXrefsForNode(CyNetwork net, CyNode node) {

        //String[] dbNames = (String[]) net.getNodeAttributeValue(node, CommonVocab.XREF_DB_NAME);
        //String[] dbIDs = (String[]) net.getNodeAttributeValue(node, CommonVocab.XREF_DB_ID);

        List dbNames = Cytoscape.getNodeAttributes().getListAttribute(node.getIdentifier(), CommonVocab.XREF_DB_NAME);
        List dbIDs = Cytoscape.getNodeAttributes().getListAttribute(node.getIdentifier(), CommonVocab.XREF_DB_ID);

        if (dbNames.size() != dbIDs.size()) {
            throw new IllegalStateException("Error: " + CommonVocab.XREF_DB_NAME + " and " + CommonVocab.XREF_DB_ID + " were not the same length for node " + node.getIdentifier() + ".");
        }
        DatabaseReference[] dbrefs = new DatabaseReference[dbNames.size()];
        for (int i = 0; i < dbNames.size(); i++) {
            String dbName = (String) dbNames.get(i);
            String dbID = (String) dbIDs.get(i);
            dbrefs[i] = new DatabaseReference(dbName, dbID);
        }
        return dbrefs;
    }
}
