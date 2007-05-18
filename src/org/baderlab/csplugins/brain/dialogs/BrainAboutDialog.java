package org.baderlab.csplugins.brain.dialogs;

import javax.swing.*;
import javax.swing.text.html.HTMLEditorKit;
import java.awt.*;

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
 * * Time: 10:59:03 PM
 * * Description
 */

/**
 * An about dialog box for BRAIN
 */
public class BrainAboutDialog extends JDialog {

    public BrainAboutDialog(Frame parentFrame) {
        super(parentFrame, "About BRAIN", false);
        setResizable(false);

        //main panel for dialog box
        JEditorPane editorPane = new JEditorPane();
        editorPane.setEditable(false);
        editorPane.setEditorKit(new HTMLEditorKit());
        //TODO: add a link to Cytoscape (see http://java.sun.com/j2se/1.4.2/docs/api/javax/swing/JEditorPane.html)
        /*
        editorPane.setText("<html><body><P align=center>BRAIN (Biologically Relevant Analysis of Interaction Networks)<BR>" +
                "v1.0 (April 2005)<BR>" +
                "Author: Gary Bader<BR><BR>\n" +
                "Sander Group<BR>\n" +
                "Computational Biology Center<BR>\n" +
                "Memorial Sloan-Kettering Cancer Center<BR>\n" +
                "</P></body></html>");
        */
        editorPane.setText("<html><body><P align=center>BRAIN (Biologically Relevant Analysis of Interaction Networks)<BR>" +
                "v1.0.4 alpha (May 15, 2007)<BR>" +
                "Author: Gary Bader<BR><BR>\n" +
                "Bader Lab<BR>\n" +
                "Terrence Donnelly Centre for Cellular and Biomolecular Research (CCBR)<BR>\n" +
                "University of Toronto<BR>\n" +
                "</P></body></html>");
        setContentPane(editorPane);
    }
}
