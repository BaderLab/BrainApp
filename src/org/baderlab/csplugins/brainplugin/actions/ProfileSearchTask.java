package org.baderlab.csplugins.brainplugin.actions;

import cytoscape.CyNetwork;
import cytoscape.Cytoscape;
import cytoscape.CytoscapeInit;
import cytoscape.task.Task;
import cytoscape.task.TaskMonitor;
import org.baderlab.csplugins.brainplugin.CytoscapeUtil;
import org.baderlab.csplugins.brainplugin.MultiSequenceSearchResultSet;
import org.baderlab.csplugins.brainplugin.BrainAlgorithm;

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
 * * User: Gary Bader
 * * Date: Apr 6, 2005
 * * Time: 6:22:02 PM
 * * Description Runs a profile search task using the task framework
 */

/**
 * Runs a profile search task using the task framework
 */
public class ProfileSearchTask implements Task {
    private TaskMonitor taskMonitor = null;
    BrainAlgorithm alg = null;

    /**
     * 'Run a profile search' task
     */
    public ProfileSearchTask() {
        alg = new BrainAlgorithm();
    }

    /**
     * Run a profile search
     */
    public void run() {
        if (taskMonitor == null) {
            throw new IllegalStateException("Task Monitor is not set.");
        }
        //run profile search
        alg.setTaskMonitor(taskMonitor);
        MultiSequenceSearchResultSet searchResults = alg.runProfileSearch();
        //disable view temporarily
        int oldViewThreshold = CytoscapeInit.getViewThreshold();
        CytoscapeInit.setViewThreshold(0);
        CyNetwork net = CytoscapeUtil.addProfileSearchResultsToCytoscape(searchResults, alg.getParams());
        CytoscapeInit.setViewThreshold(oldViewThreshold);
        Cytoscape.createNetworkView(net);
        if (searchResults != null) {
            searchResults.clear();
        }
    }

    /**
     * Non-blocking call to interrupt the task.
     */
    public void halt() {
        alg.cancel();
    }

    /**
     * Sets the Task Monitor.
     *
     * @param taskMonitor TaskMonitor Object.
     */
    public void setTaskMonitor(TaskMonitor taskMonitor) throws IllegalThreadStateException {
        if (this.taskMonitor != null) {
            throw new IllegalStateException("Task Monitor is already set.");
        }
        this.taskMonitor = taskMonitor;
    }

    /**
     * Gets the Task Title.
     *
     * @return human readable task title.
     */
    public String getTitle() {
        return new String("Profile Search");
    }
}
