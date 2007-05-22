package org.baderlab.csplugins.brainplugin;

import org.baderlab.csplugins.brainplugin.util.FileReaderUtil;

import java.awt.*;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

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
 * * Date: Feb 20, 2006
 * * Time: 6:06:41 PM
 */

/**
 * Perform hierarchical clustering bootstrap analysis
 */
public class HierarchicalClusteringBootstrapAnalysis {
    private static BrainParameterSet params = null;

    /**
     * Cluster profiles, run boostrap analysis of profile input order and output a LogoTree
     *
     * @param profileListFileName
     * @param profileLength
     * @param terminus
     */
    public static void runInputOrderRobustnessTest(String profileListFileName, int profileLength, ProteinTerminus terminus, int numberSamplings,
                                                   String logoTreeTitle, File logoTreeOutput, File leafLabelHighlightFile, File codonBiasFile) {
        params = new BrainParameterSet();
        params.setProfileFileName(new File(profileListFileName));
        params.setFuzzFactor(0.0);
        //read profile file - could be a project or single profile (list of peptides)
        List clusterProfileList = null;
        List cutProfileList = null;
        //full profiles are always stored to we can draw a logo tree
        List originalProfileList = PeptideToProfileReader.readPeptidesAsProfiles(params.getProfileFile(), 0, null, params.getFuzzFactor(), codonBiasFile, true);
        clusterProfileList = originalProfileList;
        if (profileLength > 0) {
            //profile length specified, cluster profiles of this length (not full length ones)
            cutProfileList = new ArrayList();
            for (int i = 0; i < originalProfileList.size(); i++) {
                ProteinProfile proteinProfile = (ProteinProfile) originalProfileList.get(i);
                cutProfileList.add(proteinProfile.getTruncatedProfileCopy(profileLength, terminus));
            }
            clusterProfileList = cutProfileList;
        }
        //do the initial clustering
        AvgLinkHierarchicalClustering originalCluster = clusterProfileList(clusterProfileList, true);

        //Do many clusterings for bootstrap analysis
        //Maintain information for the bootstrap results
        HierarchicalClusteringBootstrapResult bootstrapResults = new HierarchicalClusteringBootstrapResult();
        //For bootstrap analysis, randomize the input list bootstrapNumber many times to test for clustering robustness against input order
        AvgLinkHierarchicalClustering cluster = null;
        for (int i = 0; i < numberSamplings; i++) {
            //create a random ordered input list of profiles
            clusterProfileList = randomizeList(clusterProfileList);
            cluster = clusterProfileList(clusterProfileList, false);
            //capture the cluster results as a bootstrap result
            bootstrapResults.addClusterResults(cluster);
        }

        //Draw the bootstrap results on a single instance of the clustering results
        LogoTreeDraw ltd = new LogoTreeDraw(originalCluster, originalProfileList);
        ltd.setTrimNodeLogo(true, 0.15);
        ltd.setTitle(logoTreeTitle);
        ltd.setBootstrapResults(bootstrapResults);
        if (leafLabelHighlightFile != null) { // leaf label file is optional
            ArrayList list = null;
            try {
                list = FileReaderUtil.readFileAsLineList(leafLabelHighlightFile);
            } catch (IOException e) {
                e.printStackTrace();
            }
            ltd.setLeafLabelHighlightColor(list, Color.YELLOW);
        }
        ltd.setSequenceLogoStartIndex(-9);
        ltd.setSymbolStyle(new PDZProteinLogoStyle());
        ltd.outputToPDF(logoTreeOutput);
    }

    /**
     * Cluster a profile list (utility method)
     *
     * @param clusterProfileList  The protein profile list to cluster (list elements are ProteinProfile objects)
     * @param optimalLeafOrdering If true, will optimize the leaf ordering (more computationally intensive)
     * @return The clustering result
     */
    private static AvgLinkHierarchicalClustering clusterProfileList(List clusterProfileList, boolean optimalLeafOrdering) {
        AvgLinkHierarchicalClustering cluster;
        //calculate distance matrix
        DistanceMatrix distanceMatrix = new DistanceMatrix(clusterProfileList.size());
        distanceMatrix.calcDistances(clusterProfileList, new DistanceMetric() {
            public double calc(Object object1, Object object2) {
                ProteinProfile proteinProfile1 = (ProteinProfile) object1;
                ProteinProfile proteinProfile2 = (ProteinProfile) object2;
                //create a grouping for last 4 PDZ peptide residues
                String groupingByPosition = AminoAcidGrouping.getPolarChargedHydrophobeGrouping();
                return (ProteinProfileDistance.calculateAAGroupedDistributionDistance(proteinProfile1, proteinProfile2, groupingByPosition));
            }
        });
        ArrayList al = new ArrayList(clusterProfileList.size());
        for (int j = 0; j < clusterProfileList.size(); j++) {
            ProteinProfile proteinProfile = (ProteinProfile) clusterProfileList.get(j);
            al.add(j, proteinProfile.getName());
        }
        distanceMatrix.setLabels(al);
        //cluster
        cluster = new AvgLinkHierarchicalClustering(distanceMatrix);
        cluster.setOptimalLeafOrdering(optimalLeafOrdering);
        cluster.run();
        return cluster;
    }

    /**
     * Randomizes the order of an input list.
     *
     * @param inputList The list to randomize the order of
     * @return The exact objects in the input list, but in a random order.
     */
    private static List randomizeList(List inputList) {
        int originalInputListSize = inputList.size();
        ArrayList randomizedInputList = new ArrayList(originalInputListSize);
        for (int i = 0; i < originalInputListSize; i++) {
            randomizedInputList.add(inputList.remove((int) (Math.random() * inputList.size())));
        }

        return randomizedInputList;
    }

}
