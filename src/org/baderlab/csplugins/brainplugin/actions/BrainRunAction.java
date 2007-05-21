package org.baderlab.csplugins.brain.actions;

import cytoscape.Cytoscape;
import cytoscape.task.ui.JTaskConfig;
import cytoscape.task.util.TaskManager;
import org.baderlab.csplugins.brain.BrainParameterSet;
import org.baderlab.csplugins.brain.BrainAlgorithm;
import org.baderlab.csplugins.brain.dialogs.BrainParameterChangeDialog;

import javax.swing.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

/**
 * Copyright (c) 2004 Memorial Sloan-Kettering Cancer Center
 * *
 * * Code written by: Gary Bader
 * * Authors: Gary Bader, Ethan Cerami, Chris Sander
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
 * * Date: Aug 10, 2004
 * * Time: 11:02:38 PM
 * * Description
 */

/**
 * The action to run BRAIN
 */
public class BrainRunAction implements ActionListener {
    /**
     * Invoked when the run action occurs.
     */
    public void actionPerformed(ActionEvent e) {
        BrainAlgorithm pmm = new BrainAlgorithm();
        BrainParameterSet params = pmm.getParams();
        if ((!params.validateDBOptions()) || (!params.validateProfileOptions())) {
            //parameters are not valid - tell the user
            JOptionPane.showMessageDialog(Cytoscape.getDesktop(), "Please set the search parameters before running the search.\n" +
                    "You must at least set the database name and profile file name.");
            BrainParameterChangeDialog paramChangeDialog = new BrainParameterChangeDialog(Cytoscape.getDesktop());
            paramChangeDialog.pack();
            paramChangeDialog.setVisible(true);
            return;
        }
        // debug - mdharsee
        //JOptionPane.showMessageDialog(Cytoscape.getDesktop(), "[DEBUG MESSAGE] Ready to run...");

        //run algorithm
        ProfileSearchTask profileSearchTask = new ProfileSearchTask();
        //Configure JTask
        JTaskConfig config = new JTaskConfig();

        //Show Cancel/Close Buttons
        config.displayCancelButton(true);
        config.displayStatus(true);

        //Execute Task via TaskManager
        //This automatically pops-open a JTask Dialog Box
        TaskManager.executeTask(profileSearchTask, config);
    }
}
