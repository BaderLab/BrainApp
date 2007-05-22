package org.baderlab.csplugins.brainplugin.test;

import junit.framework.TestCase;
import org.biojava.bio.BioException;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.io.SeqIOTools;
import org.baderlab.csplugins.brainplugin.AlignedProteinSequenceIdentityDistance;
import org.baderlab.csplugins.brainplugin.AvgLinkHierarchicalClustering;
import org.baderlab.csplugins.brainplugin.DistanceMatrix;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;

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
 * * Date: Aug 21, 2005
 * * Time: 12:13:48 PM
 */

/**
 * Tests for AvgLinkHierarchicalClustering class
 */
public class AvgLinkHierarchicalClusteringTest extends TestCase {
    DistanceMatrix distanceMatrix = null;

    public void setUp() throws FileNotFoundException, BioException {
        //read proteins to cluster from an aligned fasta file and save as a List
        BufferedReader brCluster = new BufferedReader(new FileReader("testData" + File.separator + "optimalClusterLeafOrderingTest.txt"));
        SequenceIterator clusterProteins = (SequenceIterator) SeqIOTools.fileToBiojava("fasta", "PROTEIN", brCluster);
        ArrayList clusterProteinSequenceStringList = new ArrayList();
        ArrayList clusterProteinNameList = new ArrayList();
        ArrayList clusterProteinSequenceList = new ArrayList();
        //create lists of sequence strings and names required for next step
        while (clusterProteins.hasNext()) {
            Sequence sequence = (Sequence) clusterProteins.nextSequence();
            clusterProteinSequenceList.add(sequence);
            clusterProteinSequenceStringList.add(sequence.seqString());
            clusterProteinNameList.add(sequence.getName());
        }

        //create a distance matrix
        distanceMatrix = new DistanceMatrix(clusterProteinSequenceStringList.size());
        distanceMatrix.setLabels(clusterProteinNameList);
        distanceMatrix.calcDistances(clusterProteinSequenceStringList, new AlignedProteinSequenceIdentityDistance());
    }

    public void testOptimalLeafOrdering() {
        //cluster the sequences
        AvgLinkHierarchicalClustering cluster = new AvgLinkHierarchicalClustering(distanceMatrix);
        cluster.setOptimalLeafOrdering(true);
        cluster.run();
        System.out.println("Optimal and heuristic GTR format:");
        System.out.print(cluster.writeResultsToGTRFormat());
        System.out.println();
        System.out.println("Optimal CDT format:");
        System.out.print(cluster.toCDTString());
        System.out.println();

        cluster.setOptimalLeafOrdering(false);
        cluster.run();
        System.out.println("Heuristic CDT format:");
        System.out.print(cluster.toCDTString());
        System.out.println();

    }
}
