package org.baderlab.csplugins.brain;

import cytoscape.Cytoscape;
import cytoscape.plugin.CytoscapePlugin;
import org.baderlab.csplugins.brain.actions.BrainAboutAction;
import org.baderlab.csplugins.brain.actions.BrainParameterChangeAction;
import org.baderlab.csplugins.brain.actions.BrainRunAction;

import javax.swing.*;

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
 * * Time: 10:44:16 PM
 * * Description
 */


/**
 * Biologically Relevant Analysis of Interaction Networks (BRAIN).
 * Predicts potential protein-protein interactions from protein interaction domain motifs
 *
 * @author Gary Bader
 */
public class BrainPlugin extends CytoscapePlugin {
    /**
     * Constructor for the BRAIN plugin.
     */
    public BrainPlugin() {
        //set-up menu options in plugins menu
        JMenu menu = Cytoscape.getDesktop().getCyMenus().getOperationsMenu();
        JMenuItem item;
        //BRAIN submenu
        JMenu submenu = new JMenu("BRAIN");
        item = new JMenuItem("Run BRAIN");
        item.addActionListener(new BrainRunAction());
        submenu.add(item);
        //parameters dialog
        item = new JMenuItem("Set Parameters");
        item.addActionListener(new BrainParameterChangeAction());
        submenu.add(item);
        //About box
        item = new JMenuItem("About BRAIN");
        item.addActionListener(new BrainAboutAction());
        submenu.add(item);
        menu.add(submenu);
    }

    /**
     * Describes the plug in.
     *
     * @return short plug in description.
     */
    public String describe() {
        return new String("BRAIN plugin.");
    }
}
