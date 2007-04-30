package org.baderlab.csplugins.brain;

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
 * * Date: Feb 20, 2006
 * * Time: 10:27:42 PM
 */

/**
 * Stores the result of a hierarchical clustering bootstrap analysis, the counts of which
 * nodes in the clustering result are neighbors in the tree.
 */
public class HierarchicalClusteringBootstrapResult {
    private HashMap bootstrapResults = null; //key = ClusterNode, value = integer count of number of times this cluster node
    //appears in the clustering results

    /**
     * Constructor for an object that stores the result of a hierarchical clustering bootstrap analysis
     */
    public HierarchicalClusteringBootstrapResult() {
        bootstrapResults = new HashMap();
    }

    /**
     * Adds a cluster result to the bootstrap result
     *
     * @param cluster An individual cluster result
     */
    public void addClusterResults(AvgLinkHierarchicalClustering cluster) {
        //go through the cluster result and add it to the HashMap
        HierarchicalClusteringResultTree clusterResult = cluster.getResult();
        addClusterResultsRecursive(clusterResult);
    }

    /**
     * Helper method for addClusterResults recursive tree traversal
     *
     * @param node The node to recurse from
     */
    private void addClusterResultsRecursive(HierarchicalClusteringResultTree node) {
        //add all internal nodes to the bootstrap results
        if (node.leaf == false) {
            if (bootstrapResults.containsKey(node)) {
                Integer count = (Integer) bootstrapResults.get(node);
                bootstrapResults.put(node, new Integer(count.intValue() + 1));
            } else {
                bootstrapResults.put(node, new Integer(1));
            }
            addClusterResultsRecursive(node.left);
            addClusterResultsRecursive(node.right);
        }
        //ignore leaf nodes
    }

    /**
     * Gets the number of times this node has been seen in the bootstrap analysis
     *
     * @param node The node to get the count for
     * @return the number of times this node has been seen in the bootstrap analysis or -1 if
     *         the node has not been found (generally indicates an error)
     */
    public int getBootstrapCount(HierarchicalClusteringResultTree node) {
        if (bootstrapResults.containsKey(node)) {
            Integer count = (Integer) bootstrapResults.get(node);
            return count.intValue();
        }
        return -1;
    }
}
