package org.baderlab.csplugins.brain;

import org.biojava.bio.Annotation;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.Sequence;

import java.util.ArrayList;
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
 * * Date: Sep 1, 2005
 * * Time: 8:05:31 PM
 */

/**
 * Utility classes for extracting information from GenPept files read in by BioJava
 */
public class GenPeptUtil {
    /**
     * Gets the accession number from a GenPept formatted sequence
     *
     * @param genPeptSequence A Sequence object that was read in from a GenPept formatted file
     * @return The database reference (GenPept) to this sequence record
     */
    static public DatabaseReference getAccession(Sequence genPeptSequence) {
        Annotation seqAnn = genPeptSequence.getAnnotation();
        String dbname = null;
        String dbid = null;
        if (seqAnn.containsProperty("ACCESSION")) {
            String accession = (String) seqAnn.getProperty("ACCESSION");
            String[] accessionArray = accession.split("\\s+");
            String firstAccession = null;
            firstAccession = accessionArray[0];
            if (accessionArray.length > 0) {
                dbid = firstAccession;
            }
        }
        if (seqAnn.containsProperty("DBSOURCE")) {
            String dbSource = (String) seqAnn.getProperty("DBSOURCE");
            //extract DB source, assume format "REFSEQ: accession NM_053517.1"
            int endSource = dbSource.indexOf(":");
            if (endSource < 0) {
                return null;
            }
            String dbArray = null;
            dbArray = dbSource.substring(0, endSource);
            dbname = dbArray;
        }
        DatabaseReference xref = null;
        if (dbname != null && dbid != null) {
            xref = new DatabaseReference(dbname, dbid);
        }
        return (xref);
    }

    /**
     * Gets the accession number from a GenPept formatted sequence
     *
     * @param genPeptSequence A Sequence object that was read in from a GenPept formatted file
     * @return The database reference (GenPept) to this sequence record
     */
    static public DatabaseReference getEntrezGeneID(Sequence genPeptSequence) {
        String dbname = new String("Entrez Gene");
        String dbid = null;
        for (Iterator fi = genPeptSequence.features(); fi.hasNext();) {
            Feature f = (Feature) fi.next();
            if ((f.getAnnotation() != null) && (f.getType().equals("CDS"))) {
                //extract database xref of choice
                ArrayList dbRefs = extractGenericAnnotation("db_xref", f);
                if (dbRefs != null) {
                    for (int i = 0; i < dbRefs.size(); i++) {
                        String s = (String) dbRefs.get(i);
                        if ((s != null) && (s.startsWith("GeneID"))) {
                            if (s.indexOf(":") > 0) {
                                dbid = s.substring(s.indexOf(":") + 1, s.length());
                            }
                        }
                    }
                }
            }
        }
        DatabaseReference xref = null;
        if (dbname != null && dbid != null) {
            xref = new DatabaseReference(dbname, dbid);
        }
        return (xref);
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
     * Gets a sequence name (gene name) in the CDS feature in GenPept files
     *
     * @param sequence The sequence derived from a GenPept record
     * @return The CDS gene name, or null if not found;
     */
    public static String getGeneName(Sequence sequence) {
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

}
