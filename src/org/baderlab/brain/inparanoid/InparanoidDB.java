package org.baderlab.brain.inparanoid;

import org.baderlab.brain.DatabaseReference;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

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
 * * Date: Jun 8, 2005
 * * Time: 7:04:02 PM
 */

/**
 * An in memory representation/storage of the entire Inparanoid database
 */
//TODO: make this an interface, so that it can be implemented by multiple orthology DBs
public class InparanoidDB {
    HashMap species2taxid; //key is String species, value is Integer taxonomy ID
    HashMap species2database; //key is String species, value is String database name
    HashMap species2fullName; //key is String species, value is String organism full name
    HashMap accession2orthologGroup; //key = ENSP accession dbref, value = ArrayList of Ortholog objects (not including the key)
    HashMap otherAccession2inparanoidAccession;    //stores a mapping from any other accession numer to the accession numbers used by the InparanoidDB

    /**
     * Constructor for a InparanoidDB
     */
    public InparanoidDB() {
        species2taxid = new HashMap();
        species2database = new HashMap();
        species2fullName = new HashMap();
        accession2orthologGroup = new HashMap();
    }

    /**
     * Given an InparanoidDB species short name (e.g. ensHS), return the species full name
     * as found in the species information file
     */
    public String getSpeciesFullName(String species) {
        return ((String) species2fullName.get(species));
    }

    /**
     * Internal method for adding an ortholog to the database
     */
    private void addOrtholog(DatabaseReference sourceOrtholog, Ortholog ortholog) {
        ArrayList hg = null;
        if (sourceOrtholog.equals(ortholog.getProteinID())) {
            //don't store the key in the ArrayList, otherwise it will count as an ortholog
            return;
        }
        if (accession2orthologGroup.containsKey(sourceOrtholog)) {
            hg = (ArrayList) accession2orthologGroup.get(sourceOrtholog);
        } else {
            hg = new ArrayList();
        }
        hg.add(ortholog);
        accession2orthologGroup.put(sourceOrtholog, hg);
    }

    /**
     * Reads the Inparanoid data files, as defined here:
     * http://inparanoid.cgb.ki.se/download/current
     */
    public void readInparanoidDataFile(File[] inparanoidDataFile, File speciesInfoFile, String speciesOfInterest) throws IOException {
        //parse the species information file
        //Format: ensAG\tAnopheles gambiae\t7165\tEnsembl\n
        BufferedReader br = new BufferedReader(new FileReader(speciesInfoFile));
        String line = null;
        String[] lineSplit;
        while ((line = br.readLine()) != null) {
            lineSplit = line.split("\t");
            if (lineSplit.length == 4) {
                String species = lineSplit[0];
                String speciesFullName = lineSplit[1];
                int taxid = Integer.parseInt(lineSplit[2]);
                String datasource = lineSplit[3];
                species2taxid.put(species, new Integer(taxid));
                species2database.put(species, datasource);
                species2fullName.put(species, speciesFullName);
            }
        }
        //go through all Inparanoid data files and construct in memory database of the chosen species
        for (int i = 0; i < inparanoidDataFile.length; i++) {
            File file = inparanoidDataFile[i];
            BufferedReader brData = new BufferedReader(new FileReader(file));
            ArrayList orthologsInCluster = new ArrayList(); //two or more orthologs per cluster
            ArrayList speciesOfInterestProteinsInCluster = new ArrayList(); //one or more species of interest proteins per cluster
            int count = 0;
            int prevClusterID = -1;
            while ((line = brData.readLine()) != null) {
                count++;
                lineSplit = line.split("\t");
                if (lineSplit.length < 5) {
                    throw new IOException("Inparanoid file" + file + " line " + count + " is too short. Expected at least 5 fields, but only found " + lineSplit.length + ".");
                }
                int clusterID = Integer.parseInt(lineSplit[0]);
                String species = lineSplit[2];
                String accession = lineSplit[4];
                Object taxidLookup = species2taxid.get(species);
                if (taxidLookup == null) {
                    throw new IllegalArgumentException("Wasn't able to find an NCBI taxonomy ID for species: " + species +
                            ". Please check the species info file (" + speciesInfoFile.toString() + ") to make sure it is specified.");
                }
                int ncbiTaxID = ((Integer) taxidLookup).intValue();
                DatabaseReference dbref = new DatabaseReference((String) species2database.get(species), accession);
                Ortholog ortholog = null;
                if (!species.equalsIgnoreCase(speciesOfInterest)) {
                    //only create a new ortholog for non-species of interest proteins
                    ortholog = new Ortholog(clusterID, ncbiTaxID, dbref);
                }
                if ((count > 1) && (clusterID != prevClusterID)) { //except for line 1
                    //collect orthologs until we hit the next cluster (clusters are coded in multiple lines)
                    //save the previous cluster and start the next one
                    if (speciesOfInterestProteinsInCluster.size() == 0) {
                        throw new RuntimeException("Wasn't able to find a species of interest in datafiles. Currently specified: " + speciesOfInterest + ".");
                    }
                    for (int j = 0; j < orthologsInCluster.size(); j++) {
                        Ortholog clusterOrtholog = (Ortholog) orthologsInCluster.get(j);
                        for (int k = 0; k < speciesOfInterestProteinsInCluster.size(); k++) {
                            //add this ortholog to all species of interest proteins found in the cluster
                            DatabaseReference dbrefSpeciesOfInterest = (DatabaseReference) speciesOfInterestProteinsInCluster.get(k);
                            addOrtholog(dbrefSpeciesOfInterest, clusterOrtholog);
                        }
                    }
                    //reset for the next cluster
                    speciesOfInterestProteinsInCluster.clear();
                    orthologsInCluster.clear();
                }
                if (!species.equalsIgnoreCase(speciesOfInterest)) {
                    //only save the new ortholog for non-species of interest proteins
                    orthologsInCluster.add(ortholog);
                }
                if (species.equals(speciesOfInterest)) {
                    speciesOfInterestProteinsInCluster.add(dbref);
                }
                prevClusterID = clusterID;
                /*
                As of September 2005, the Inparanoid format is defined as a tab-delimited list with the following fields:
                sqltable.?????-?????		Output files containing all Inparanoid clusters for
                each species pair in table format. See species fasta
                file list above for species abbreviations.
                e.g. sqltable.ensHS-ensCE; All Inparanoid
                clusters between Homo sapiens and Caenorhabditis
                elegans
                Each column represents:
                1 - Cluster number
                2 - Seed ortholog-pair blast score in bits
                3 - Species abbreviation
                4 - Inparanoid score
                5 - Protein identifier
                6 - Bootstrap value for seed inparalog/ortholog (present if user performed bootstrap analysis)
                */
            }
        }
    }

    //read in a mapping file from ids of interest to IDs that inparanoid knows about
    /**
     * Read an id mapping file from an inparanoid accession type of interest to another accession
     * type that might be used to query this inparanoid database object
     *
     * @param idMappingFile The file should be formatted to map inparanoid accessions to the accession
     *                      of interest.  Each ID must be formatted "databasename:id".  IDs are separated by a tab and at most
     *                      two IDs are allowed per line. E.g.
     *                      ENSEMBL:ENSP00000236850\tREFSEQ:NP_000030\n
     *                      ENSEMBL:ENSP00000352471\tREFSEQ:NP_000030\n
     *                      ENSEMBL:ENSP00000227667\tREFSEQ:NP_000031\n
     *                      ENSEMBL:ENSP00000353540\tREFSEQ:NP_000031\n
     *                      ...
     * @throws IOException if an error occurs while reading the mapping file
     */
    public void readIDMappingFile(File idMappingFile) throws IOException {
        if (otherAccession2inparanoidAccession == null) {
            otherAccession2inparanoidAccession = new HashMap();
        }
        BufferedReader br = new BufferedReader(new FileReader(idMappingFile));
        String line = null;
        String[] lineSplit;
        while ((line = br.readLine()) != null) {
            lineSplit = line.split("\t");
            if (lineSplit.length == 2) {
                String id1 = lineSplit[0];
                DatabaseReference dbref1 = new DatabaseReference(id1);
                String id2 = lineSplit[1];
                DatabaseReference dbref2 = new DatabaseReference(id2);
                otherAccession2inparanoidAccession.put(dbref2, dbref1);
            }
        }
    }

    /**
     * Tests if two accession numbers are orthologs (according to the Inparanoid database)
     *
     * @param referenceAccession  Reference accession number - is B and ortholog of A?
     *                            This is the accession number used to query the ortholog database, which must support
     *                            lookups on this accession number.
     * @param comparisonAccession Accession number B
     * @return true if the two accession numbers are orthologs
     */
    public boolean isOrthologByAccession(DatabaseReference referenceAccession, DatabaseReference comparisonAccession) {
        //try the accession numbers given directly
        ArrayList accAOrthoList = findOrthologList(referenceAccession);
        if (accAOrthoList == null) {
            return false; //can't find any ortholog information
        }
        for (int i = 0; i < accAOrthoList.size(); i++) {
            Ortholog ortholog = (Ortholog) accAOrthoList.get(i);
            if (ortholog.getProteinID().equals(comparisonAccession)) {
                return true;
            }
        }
        return false;
    }

    /**
     * Internal helper method
     */
    private ArrayList findOrthologList(DatabaseReference accession) {
        ArrayList accessionOrthoList = null;
        if (accession2orthologGroup.containsKey(accession)) {
            accessionOrthoList = (ArrayList) accession2orthologGroup.get(accession);
        } else {
            if (otherAccession2inparanoidAccession != null) {
                //try to map from another accession type of interest, if available
                if (otherAccession2inparanoidAccession.containsKey(accession)) {
                    accession = (DatabaseReference) otherAccession2inparanoidAccession.get(accession);
                    if (accession2orthologGroup.containsKey(accession)) {
                        accessionOrthoList = (ArrayList) accession2orthologGroup.get(accession);
                    }
                }
            }
        }
        return (accessionOrthoList);
    }

}