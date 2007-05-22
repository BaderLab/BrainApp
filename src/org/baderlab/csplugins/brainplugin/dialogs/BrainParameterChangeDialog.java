package org.baderlab.csplugins.brainplugin.dialogs;

import cytoscape.util.FileUtil;
import org.baderlab.brain.BrainCurrentParameters;
import org.baderlab.brain.BrainParameterSet;
import org.baderlab.brain.ProteinDatabaseSearchParams;
import org.baderlab.brain.ProteinTerminus;
import org.baderlab.csplugins.brainplugin.util.JMultiLineToolTip;

import javax.swing.*;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;
import java.text.DecimalFormat;

/**
 * * Copyright (c) 2004 Memorial Sloan-Kettering Cancer Center
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
 ** User: Gary Bader
 ** Date: Feb 6, 2005
 ** Time: 5:00:00 PM
 ** Description: The parameter change dialog which the user can use to change the parameters
 **/

/**
 * The parameter change dialog which the user can use to change the parameters
 */
public class BrainParameterChangeDialog extends JDialog {
    /*
    * WARNING: there is a little bit of messy logic in this class to make sure that the user can only enter
    * useful parameter combinations.  Before modifying 'if' statements, check side-effects.
    */

    //Parameters for BRAIN
    BrainParameterSet currentParamsCopy;    //stores current parameters - populates dialog box fields

    //resetable UI elements
    //Tab 1: database selection
    JTextField databaseFileNameTextField;
    JComboBox databaseFormatComboBox;
    //general search params
    JComboBox terminusComboBox;
    boolean wholeSequenceSearchSet = false;
    JFormattedTextField lengthFormattedTextField;
    JCheckBox multipleHitsCheckBox;
    //Tab 2: Profile search parameters
    //where to find the profiles, threshold, normalized
    JTextField profileFileNameTextField;
    JCheckBox normalizedCheckBox;
    JFormattedTextField scoreThresholdFormattedTextField;
    JFormattedTextField topHitsFormattedTextField;
    JFormattedTextField fuzzFactorFormattedTextField;
    JTextField biasFileNameTextField;
    JCheckBox uniquePeptidesCheckBox;
    //Tab 3: Advanced options
    JCheckBox groupProteinDomains;
    JComboBox nodeRepComboBox;  // what a node represents (e.g. domain, protein, etc)

    /**
     * The actual parameter change dialog that builds the UI
     *
     * @param parentFrame The parent frame for this dialog
     */
    public BrainParameterChangeDialog(Frame parentFrame) {
        super(parentFrame, "Parameters", false);
        setResizable(true);

        //get the current parameters (get a copy to avoid setting the parameters until user presses OK)
        currentParamsCopy = BrainCurrentParameters.getInstance().getParamsCopy();

        //main panel for dialog box
        JPanel panel = new JPanel(new BorderLayout());
        panel.setSize(500, 800);

        //Tab 1: Database selection
        Box databaseBox = Box.createVerticalBox();
        //databaseChooserSubPanel
        JPanel databaseChooserSubPanel = new JPanel(new BorderLayout());
        databaseChooserSubPanel.setBorder(BorderFactory.createTitledBorder("Database to Search"));

        //Database file chooser
        JPanel chooserPanel = new JPanel();
        databaseFileNameTextField = new JTextField();
        databaseFileNameTextField.setToolTipText("The file name of the database to search.");
        if (currentParamsCopy.getDatabaseFileName() != null) {
            databaseFileNameTextField.setText(currentParamsCopy.getDatabaseFileName().toString());
        }
        databaseFileNameTextField.setColumns(20);
        //this only handles VK_ENTER key events
        databaseFileNameTextField.addActionListener(new BrainParameterChangeDialog.databaseFileNameAction());
        //this handles other text change events
        databaseFileNameTextField.getDocument().addDocumentListener(new databaseFileNameAction());
        chooserPanel.add(databaseFileNameTextField);
        JButton chooseFileButton = new JButton("Browse...");
        chooseFileButton.addActionListener(new BrainParameterChangeDialog.chooseDatabaseAction());
        chooserPanel.add(chooseFileButton);
        databaseChooserSubPanel.add(chooserPanel, BorderLayout.NORTH);

        //File format combo box
        JLabel formatLabel = new JLabel("Format");
        JPanel formatLabelFieldPanel = new JPanel();
        String tipText = "Select the format of the database file.";
        formatLabelFieldPanel.setToolTipText(tipText);
        formatLabelFieldPanel.add(formatLabel);
        String[] formatChoices = {"FASTA", "SwissProt", "GenPept", "EMBL"};
        databaseFormatComboBox = new JComboBox(formatChoices);
        databaseFormatComboBox.setToolTipText(tipText);
        //if there is a format already set, maintain it (otherwise, default is used)
        if (currentParamsCopy.getDatabaseFormat() != null) {
            String currentFormat = currentParamsCopy.getDatabaseFormat();
            for (int i = 0; i < formatChoices.length; i++) {
                String formatChoice = formatChoices[i];
                if (formatChoice.equalsIgnoreCase(currentFormat)) {
                    databaseFormatComboBox.setSelectedIndex(i);
                    break;
                }
            }
        } else {
            //use the default (defined as the first format in the formatChoices list
            currentParamsCopy.setDatabaseFormat(formatChoices[0]);
        }
        formatLabelFieldPanel.add(databaseFormatComboBox);
        databaseFormatComboBox.addActionListener(new selectFormatAction());
        databaseChooserSubPanel.add(formatLabelFieldPanel, BorderLayout.SOUTH);

        //searchOptionsSubPanel
        JPanel searchOptionsSubPanel = new JPanel(new BorderLayout());
        searchOptionsSubPanel.setBorder(BorderFactory.createTitledBorder("Sequence Filtering"));

        //length (put this before the terminus set-up so we can set the enable status based on the terminus status)
        DecimalFormat decFormatInteger = new DecimalFormat();
        decFormatInteger.setParseIntegerOnly(true);
        lengthFormattedTextField = new JFormattedTextField(decFormatInteger) {
            public JToolTip createToolTip() {
                return new JMultiLineToolTip();
            }
        };
        lengthFormattedTextField.setColumns(3);
        lengthFormattedTextField.addPropertyChangeListener("value", new BrainParameterChangeDialog.formattedTextFieldAction());
        String lengthTipText = "Search only this many residues from either\nN or C terminal end of each protein in the database.\n" +
                "If not set, will use the length of the search pattern.";
        lengthFormattedTextField.setToolTipText(lengthTipText);
        if (currentParamsCopy.getSearchParams() != null) {
            int length = currentParamsCopy.getSearchParams().getLength();
            if (length > 0) {
                lengthFormattedTextField.setText((Integer.toString(length)));
            } else {
                lengthFormattedTextField.setText("5");
                currentParamsCopy.getSearchParams().setLength(5);
            }
        }
        JLabel lengthLabel = new JLabel("Search Length");
        JPanel labelFieldPanel = new JPanel() {
            public JToolTip createToolTip() {
                return new JMultiLineToolTip();
            }
        };
        labelFieldPanel.setToolTipText(lengthTipText);
        labelFieldPanel.add(lengthLabel);
        labelFieldPanel.add(lengthFormattedTextField);
        searchOptionsSubPanel.add(labelFieldPanel, BorderLayout.WEST);

        //Sequence filter option combo box
        JLabel searchLabel = new JLabel("Search");
        JPanel searchLabelFieldPanel = new JPanel();
        String searchTipText = "Optionally choose a subset of the sequence to search.";
        searchLabelFieldPanel.setToolTipText(searchTipText);
        searchLabelFieldPanel.add(searchLabel);
        String[] terminiChoices = {"Whole sequence", "N-terminus", "C-terminus"};
        ProteinTerminus[] terminiList = {ProteinTerminus.NONE, ProteinTerminus.N, ProteinTerminus.C}; //assume that terminiList matches terminiChoices
        terminusComboBox = new JComboBox(terminiChoices);
        terminusComboBox.setToolTipText(searchTipText);
        //update the length format field
        lengthFormattedTextField.setEnabled(false);
        //if there is a terminus already set, maintain it (otherwise, default is used)
        if (currentParamsCopy.getSearchParams() != null) {
            ProteinTerminus currentTerminus = currentParamsCopy.getSearchParams().getTerminus();
            //current terminus is guaranteed to return a valid terminus
            for (int i = 0; i < terminiList.length; i++) {
                ProteinTerminus terminus = terminiList[i];
                if (terminus.equals(currentTerminus)) { //assume that terminiList matches terminiChoices
                    terminusComboBox.setSelectedIndex(i);
                    break;
                }
            }
            if (!currentTerminus.equals(ProteinTerminus.NONE)) {
                if ((currentParamsCopy.getSearchParams() != null) && (currentParamsCopy.getSearchParams().isMultipleHits())) {
                    lengthFormattedTextField.setEnabled(true);
                }
                wholeSequenceSearchSet = false;
            } else {
                wholeSequenceSearchSet = true;
            }
        }
        terminusComboBox.addActionListener(new selectTerminiAction());
        searchLabelFieldPanel.add(terminusComboBox);
        searchOptionsSubPanel.add(searchLabelFieldPanel, BorderLayout.NORTH);

        //multiple hits allowed checkbox
        multipleHitsCheckBox = new JCheckBox("Multiple Hits?", false) {
            public JToolTip createToolTip() {
                return new JMultiLineToolTip();
            }
        };
        multipleHitsCheckBox.addItemListener(new multipleHitsCheckBoxAction());
        multipleHitsCheckBox.setToolTipText("If checked, search will return multiple hits,\notherwise only the first hit will be returned.");
        if (currentParamsCopy.getSearchParams() != null) {
            multipleHitsCheckBox.setSelected(currentParamsCopy.getSearchParams().isMultipleHits());
        }
        searchOptionsSubPanel.add(multipleHitsCheckBox, BorderLayout.EAST);

        //put tab 1 together
        databaseBox.add(databaseChooserSubPanel);
        databaseBox.add(searchOptionsSubPanel);
        databaseBox.add(Box.createRigidArea(new Dimension(0, 50)));

        //Tab 2: Profile search parameters
        JPanel profileParamsPanel = new JPanel(new BorderLayout());
        //profileChooserSubPanel
        JPanel profileChooserSubPanel = new JPanel(new BorderLayout());
        profileChooserSubPanel.setBorder(BorderFactory.createTitledBorder("Profile selection"));

        JPanel profileChooserPanel = new JPanel();
        profileFileNameTextField = new JTextField();
        profileFileNameTextField.setToolTipText("The file name of the profile or list of profiles to search using.");
        if (currentParamsCopy.getProfileFile() != null) {
            profileFileNameTextField.setText(currentParamsCopy.getProfileFile().toString());
        }
        profileFileNameTextField.setColumns(20);
        //this only handles VK_ENTER key events
        profileFileNameTextField.addActionListener(new BrainParameterChangeDialog.profileFileNameAction());
        //this handles other text change events
        profileFileNameTextField.getDocument().addDocumentListener(new profileFileNameAction());
        profileChooserPanel.add(profileFileNameTextField);
        JButton profileChooseFileButton = new JButton("Browse...");
        profileChooseFileButton.addActionListener(new BrainParameterChangeDialog.chooseProfileAction());
        profileChooserPanel.add(profileChooseFileButton);
        profileChooserSubPanel.add(profileChooserPanel, BorderLayout.NORTH);

        //subpanel for unique peptides checkbox
        JPanel profileChooserSubSubPanel = new JPanel();
        uniquePeptidesCheckBox = new JCheckBox("Unique peptides", false) {
            public JToolTip createToolTip() {
                return new JMultiLineToolTip();
            }
        };
        uniquePeptidesCheckBox.addItemListener(new uniquePeptidesCheckBoxAction());
        uniquePeptidesCheckBox.setToolTipText("If checked, only unique peptide sequences will be loaded from the profile.");
        uniquePeptidesCheckBox.setSelected(currentParamsCopy.getUniquePeptides());
        profileChooserSubSubPanel.add(uniquePeptidesCheckBox);
        profileChooserSubPanel.add(profileChooserSubSubPanel, BorderLayout.SOUTH);


        //profileOptionsSubPanel
        JPanel profileOptionsSubPanel = new JPanel(new BorderLayout());
        profileOptionsSubPanel.setBorder(BorderFactory.createTitledBorder("Profile search options"));

        //profileOptionsSubSubPanel
        JPanel profileOptionsSubSubPanel = new JPanel(new FlowLayout());

        //scoreTypeSubPanel
        JPanel scoreTypeSubPanel = new JPanel(new BorderLayout());
        scoreTypeSubPanel.setBorder(BorderFactory.createTitledBorder("Score Threshold"));
        profileOptionsSubSubPanel.add(scoreTypeSubPanel);
        profileOptionsSubPanel.add(profileOptionsSubSubPanel, BorderLayout.NORTH);

        //profileOptionsSubSubPanel
        Box scoreTypeInputBox = new Box(BoxLayout.Y_AXIS);

        //score threshold field
        DecimalFormat decFormatDouble = new DecimalFormat();
        decFormatDouble.setMaximumIntegerDigits(3);
        scoreThresholdFormattedTextField = new JFormattedTextField(decFormatDouble) {
            public JToolTip createToolTip() {
                return new JMultiLineToolTip();
            }
        };
        scoreThresholdFormattedTextField.setColumns(5);
        scoreThresholdFormattedTextField.addPropertyChangeListener("value", new formattedTextFieldAction());
        String scoreTipText = "Sets the score threshold for the search.\n" +
                "Matches below this score will be added to Cytoscape.\n" +
                "E.g. 10 means only return matches to the profile lower than p=1E-10\n" +
                "The lower this number, the less matches will result\n" +
                "This number must be 0 or higher where 0 (p=1.0) is the most stringent threshold\n" +
                "where the score of the match must be equal to p=1.0 (perfect match).";
        scoreThresholdFormattedTextField.setToolTipText(scoreTipText);
        scoreThresholdFormattedTextField.setText(Double.toString(currentParamsCopy.getScoreThreshold()));
        JLabel scoreThresholdLabel = new JLabel("Score threshold");
        JPanel labelFieldPanel1 = new JPanel() {
            public JToolTip createToolTip() {
                return new JMultiLineToolTip();
            }
        };
        labelFieldPanel1.setToolTipText(scoreTipText);
        labelFieldPanel1.add(scoreThresholdLabel);
        labelFieldPanel1.add(scoreThresholdFormattedTextField);
        scoreTypeInputBox.add(labelFieldPanel1);

        //score threshold field
        topHitsFormattedTextField = new JFormattedTextField(decFormatInteger) {
            public JToolTip createToolTip() {
                return new JMultiLineToolTip();
            }
        };
        topHitsFormattedTextField.setColumns(5);
        topHitsFormattedTextField.addPropertyChangeListener("value", new formattedTextFieldAction());
        String topHistTipText = "Sets the score threshold for the search\n" +
                "in terms of number of top hits. More hits than this\n" +
                "may result if many hits have the same score.\n";
        topHitsFormattedTextField.setToolTipText(topHistTipText);
        topHitsFormattedTextField.setText(Double.toString(currentParamsCopy.getNumberTopHits()));
        JLabel topHistThresholdLabel = new JLabel("Number top hits");
        JPanel labelFieldPanel2 = new JPanel() {
            public JToolTip createToolTip() {
                return new JMultiLineToolTip();
            }
        };
        labelFieldPanel2.setToolTipText(topHistTipText);
        labelFieldPanel2.add(topHistThresholdLabel);
        labelFieldPanel2.add(topHitsFormattedTextField);
        scoreTypeInputBox.add(labelFieldPanel2);

        scoreTypeSubPanel.add(scoreTypeInputBox, BorderLayout.SOUTH);

        //bottom flow panel
        JPanel profileOptionsSubPanelBottomPanel = new JPanel(new FlowLayout());

        //score threshold field
        fuzzFactorFormattedTextField = new JFormattedTextField(decFormatDouble) {
            public JToolTip createToolTip() {
                return new JMultiLineToolTip();
            }
        };
        fuzzFactorFormattedTextField.setColumns(5);
        fuzzFactorFormattedTextField.addPropertyChangeListener("value", new formattedTextFieldAction());
        String fuzzFactorTipText = "Sets the fuzz factor for the profile.\n" +
                "The higher this number, the fuzzier the profile\n" +
                "will be and the more hits will result.\n" +
                "(This is the pseudocount number for the profile)";
        fuzzFactorFormattedTextField.setToolTipText(fuzzFactorTipText);
        fuzzFactorFormattedTextField.setText(Double.toString(currentParamsCopy.getFuzzFactor()));
        JLabel fuzzFactorLabel = new JLabel("Fuzz Factor");
        JPanel fuzzFactorLabelFieldPanel = new JPanel() {
            public JToolTip createToolTip() {
                return new JMultiLineToolTip();
            }
        };
        fuzzFactorLabelFieldPanel.setToolTipText(fuzzFactorTipText);
        fuzzFactorLabelFieldPanel.add(fuzzFactorLabel);
        fuzzFactorLabelFieldPanel.add(fuzzFactorFormattedTextField);
        profileOptionsSubPanelBottomPanel.add(fuzzFactorLabelFieldPanel, BorderLayout.SOUTH);

        //normalized check box
        normalizedCheckBox = new JCheckBox("Normalized", false) {
            public JToolTip createToolTip() {
                return new JMultiLineToolTip();
            }
        };
        normalizedCheckBox.addItemListener(new normalizedCheckBoxAction());
        normalizedCheckBox.setToolTipText("If checked, match scores will be normalized to the profile.\n" +
                "(To the range of possible p-values that the profile can produce)");
        if (!currentParamsCopy.getSearchParams().isNormalizedSet()) {
            //if normalized flag is not set, we are creating the dialog for the first time - set the
            //default to be 'normalized'
            //TODO: remove the normalized check box and switch to a drop down menu of score types
            currentParamsCopy.getSearchParams().setNormalized(true);
            currentParamsCopy.getSearchParams().setScoreType(ProteinDatabaseSearchParams.SCORE_TYPE_NORM_PROBABILITY);
        }
        normalizedCheckBox.setSelected(currentParamsCopy.getSearchParams().isNormalized());
        profileOptionsSubPanelBottomPanel.add(normalizedCheckBox);

        profileOptionsSubPanel.add(profileOptionsSubPanelBottomPanel, BorderLayout.SOUTH);

        // Codon bias options subpanel
        JPanel biasFileSubPanel = new JPanel(new BorderLayout());
        biasFileSubPanel.setBorder(BorderFactory.createTitledBorder("Codon Bias File"));
        JPanel biasFileChooserPanel = new JPanel();
        biasFileNameTextField = new JTextField();
        biasFileNameTextField.setToolTipText("The file name of the codon bias file specification.");
        if (currentParamsCopy.getCodonBiasFile() != null) {
            biasFileNameTextField.setText(currentParamsCopy.getCodonBiasFile().toString());
        }
        biasFileNameTextField.setColumns(20);
        //this only handles VK_ENTER key events
        biasFileNameTextField.addActionListener(new BrainParameterChangeDialog.codonBiasFileNameAction());
        biasFileChooserPanel.add(biasFileNameTextField);
        JButton biasFileChooseFileButton = new JButton("Browse...");
        biasFileChooseFileButton.addActionListener(new chooseCodonBiasAction());
        biasFileChooserPanel.add(biasFileChooseFileButton);
        biasFileSubPanel.add(biasFileChooserPanel, BorderLayout.NORTH);


        //put tab 2 together
        JPanel tab2Panel = new JPanel(new BorderLayout());
        tab2Panel.add(profileChooserSubPanel, BorderLayout.NORTH);
        tab2Panel.add(profileOptionsSubPanel, BorderLayout.SOUTH);
        tab2Panel.add(biasFileSubPanel);
        profileParamsPanel.add(tab2Panel, BorderLayout.NORTH);


        //Tab 3: Advanced options
        Box advancedBox = Box.createVerticalBox();

        //network options
        JPanel networkSubPanel = new JPanel(new BorderLayout());
        networkSubPanel.setBorder(BorderFactory.createTitledBorder("Network Options"));

        //Node representation combo box
        JLabel nodeRepLabel = new JLabel("Nodes represent");
        JPanel nodeRepLabelFieldPanel = new JPanel();
        tipText = "A node in the network can represent a single domain, or a protein with one or more instances of a domain";
        nodeRepLabelFieldPanel.setToolTipText(tipText);
        nodeRepLabelFieldPanel.add(nodeRepLabel);
        String[] nodeRepChoices = {"Proteins","Domains"};
        nodeRepComboBox = new JComboBox(nodeRepChoices);
        nodeRepComboBox.setToolTipText(tipText);
        //set default value
        if (currentParamsCopy.getUniqueQueryProteinNodes()) {
            nodeRepComboBox.setSelectedIndex(0);
        } else {
            nodeRepComboBox.setSelectedIndex(1);
        }
        nodeRepLabelFieldPanel.add(nodeRepComboBox);
        nodeRepComboBox.addActionListener(new selectNodeRepAction());
        networkSubPanel.add(nodeRepLabelFieldPanel, BorderLayout.CENTER);

        //put Tabl 3 together
        advancedBox.add(networkSubPanel);
        advancedBox.add(Box.createRigidArea(new Dimension(0, 225)));


        //create the tabbed pane
        JTabbedPane tabbedPane = new JTabbedPane();
        //TODO: add a first, default tab that just has basic instructions for the use of the rest of the tabs
        tabbedPane.addTab("Sequence Database", null, databaseBox, "Choose the protein sequence database to search.");
        tabbedPane.addTab("Profile Search", null, profileParamsPanel, "Set parameters for profile search.");
        tabbedPane.addTab("Advanced Options", null, advancedBox, "Set advanced options.");
        panel.add(tabbedPane, BorderLayout.CENTER);

        //create the bottom panel
        JPanel bottomPanel = new JPanel(new FlowLayout());

        JButton OKButton = new JButton("OK");
        OKButton.addActionListener(new BrainParameterChangeDialog.OKAction(this));
        bottomPanel.add(OKButton);

        JButton cancelButton = new JButton("Cancel");
        cancelButton.addActionListener(new BrainParameterChangeDialog.cancelAction(this));
        bottomPanel.add(cancelButton);

        panel.add(bottomPanel, BorderLayout.SOUTH);

        setContentPane(panel);
    }

    /**
     * Saves the currently set parameters
     */
    private void saveParams() {
        BrainCurrentParameters.getInstance().setParams(currentParamsCopy);
    }

    /**
     * Action for the OK button (saves parameters)
     */
    private class OKAction extends AbstractAction {
        private JDialog dialog;

        OKAction(JDialog popup) {
            super();
            this.dialog = popup;
        }

        public void actionPerformed(ActionEvent e) {
            //validate parameters
            if (!currentParamsCopy.validateDBOptions()) {
                JOptionPane.showMessageDialog(dialog, "Please set the database name.");
                return;
            }
            //make sure that disabled fields are not set in the parameter object
            if (!lengthFormattedTextField.isEnabled()) {
                if (currentParamsCopy.getSearchParams() != null) {
                    currentParamsCopy.getSearchParams().setLength(null);
                }
            }
            saveParams();
            dialog.dispose();
        }
    }

    /**
     * Action for the cancel button (does not save parameters)
     */
    private class cancelAction extends AbstractAction {
        private JDialog dialog;

        cancelAction(JDialog popup) {
            super();
            this.dialog = popup;
        }

        public void actionPerformed(ActionEvent e) {
            dialog.dispose();
        }
    }

    /**
     * Handles choosing a database
     */
    private class chooseDatabaseAction extends AbstractAction {
        chooseDatabaseAction() {
            super();
        }

        public void actionPerformed(ActionEvent e) {
            File chosenFile = FileUtil.getFile("Select protein database...", FileUtil.LOAD);
            if (chosenFile != null) {
                currentParamsCopy.setDatabaseFileName(chosenFile);
                databaseFileNameTextField.setText(chosenFile.toString());
            }
        }
    }

    /**
     * Handles choosing a database format
     */
    private class selectFormatAction extends AbstractAction {
        selectFormatAction() {
            super();
        }

        public void actionPerformed(ActionEvent event) {
            if ("comboBoxChanged".equals(event.getActionCommand())) {
                JComboBox cb = (JComboBox) event.getSource();
                String formatName = (String) cb.getSelectedItem();
                currentParamsCopy.setDatabaseFormat(formatName);
            }
        }
    }

    /**
     * Handles choosing a sequence filter option
     */
    private class selectTerminiAction extends AbstractAction {
        String[] terminiChoices = {"Whole sequence", "N-terminus", "C-terminus"};
        ProteinTerminus[] terminiList = {ProteinTerminus.NONE, ProteinTerminus.N, ProteinTerminus.C}; //assume that terminiList matches terminiChoices

        selectTerminiAction() {
            super();
        }

        public void actionPerformed(ActionEvent event) {
            if ("comboBoxChanged".equals(event.getActionCommand())) {
                JComboBox cb = (JComboBox) event.getSource();
                String terminiSelection = (String) cb.getSelectedItem();
                if (!terminiSelection.equals(terminiChoices[0])) {  //N or C selected
                    if ((currentParamsCopy.getSearchParams() != null) && (currentParamsCopy.getSearchParams().isMultipleHits())) {
                        lengthFormattedTextField.setEnabled(true);
                    }
                    wholeSequenceSearchSet = false;
                } else { //whole sequence selected
                    lengthFormattedTextField.setEnabled(false);
                    wholeSequenceSearchSet = true;
                }
                //find termini choice and set it in the parameter object
                for (int i = 0; i < terminiChoices.length; i++) {
                    String terminiChoice = terminiChoices[i];
                    if (terminiChoice.equals(terminiSelection)) {
                        if (currentParamsCopy.getSearchParams() != null) {
                            currentParamsCopy.getSearchParams().setTerminus(terminiList[i]);
                        }
                        break;
                    }
                }
            }
        }
    }

    /**
     * Handles setting of the multiple hits parameter
     */
    private class multipleHitsCheckBoxAction implements ItemListener {
        public void itemStateChanged(ItemEvent e) {
            ProteinDatabaseSearchParams params = null;
            if (currentParamsCopy.getSearchParams() != null) {
                params = currentParamsCopy.getSearchParams();
            }
            if (e.getStateChange() == ItemEvent.DESELECTED) {
                if (params != null) {
                    params.setMultipleHits(false);
                    //If N or C terminus is set and multiple hits is false, then length format field should be disabled
                    //Length of profile will be used
                    if (lengthFormattedTextField.isEnabled()) {
                        //N or C terminus is set
                        lengthFormattedTextField.setEnabled(false);
                    }
                }
            } else {
                if (params != null) {
                    params.setMultipleHits(true);
                    if (!wholeSequenceSearchSet) {
                        lengthFormattedTextField.setEnabled(true);
                    }
                }
            }
        }
    }

    /**
     * Handles choosing a profile
     */
    private class chooseProfileAction extends AbstractAction {

        chooseProfileAction() {
            super();
        }

        public void actionPerformed(ActionEvent e) {
            File chosenFile = FileUtil.getFile("Select profiles...", FileUtil.LOAD);
            if (chosenFile != null) {
                currentParamsCopy.setProfileFileName(chosenFile);
                profileFileNameTextField.setText(chosenFile.toString());
            }
        }
    }

    /**
     * Handles choosing a codon bias file
     */
    private class chooseCodonBiasAction extends AbstractAction {

        chooseCodonBiasAction() {
            super();
        }

        public void actionPerformed(ActionEvent e) {
            File chosenFile = FileUtil.getFile("Select Codon Bias Specification...", FileUtil.LOAD);
            if (chosenFile != null) {
                currentParamsCopy.setCodonBiasFileName(chosenFile);
                biasFileNameTextField.setText(chosenFile.toString());
            }
        }
    }

    /**
     * Handles setting of the 'unique peptides' parameter
     */
    private class uniquePeptidesCheckBoxAction implements ItemListener {

        public void itemStateChanged(ItemEvent e) {
            if (e.getStateChange() == ItemEvent.DESELECTED) {
                currentParamsCopy.setUniquePeptides(false);
            }
            else if (e.getStateChange() == ItemEvent.SELECTED) {
                currentParamsCopy.setUniquePeptides(true);
            }
        }
    }
    /**
     * Handles setting of the normalized search parameter
     */
    private class normalizedCheckBoxAction implements ItemListener {
        public void itemStateChanged(ItemEvent e) {
            ProteinDatabaseSearchParams params = null;
            if (currentParamsCopy.getSearchParams() != null) {
                params = currentParamsCopy.getSearchParams();
            }
            if (e.getStateChange() == ItemEvent.DESELECTED) {
                if (params != null) {
                    params.setNormalized(false);
                    currentParamsCopy.getSearchParams().setScoreType(ProteinDatabaseSearchParams.SCORE_TYPE_PROBABILITY);
                }
            } else {
                if (params != null) {
                    params.setNormalized(true);
                    currentParamsCopy.getSearchParams().setScoreType(ProteinDatabaseSearchParams.SCORE_TYPE_NORM_PROBABILITY);
                }
            }
        }
    }

    /**
     * Handles setting database filename if the user edits the text field directly without using the file chooser
     */
    private class databaseFileNameAction extends AbstractAction implements DocumentListener {

        databaseFileNameAction() {
            super();
        }

        //any change that is made is handled the same - parameter is updated
        private void handleChange() {
            String value = databaseFileNameTextField.getText();
            if (value != null) {
                File file = new File(value);
                currentParamsCopy.setDatabaseFileName(file);
            }
        }

        public void actionPerformed(ActionEvent e) {
            handleChange();
        }

        public void changedUpdate(DocumentEvent e) {
            handleChange();
        }

        public void insertUpdate(DocumentEvent e) {
            handleChange();
        }

        public void removeUpdate(DocumentEvent e) {
            handleChange();
        }
    }

    /**
     * Handles setting profile filename if the user edits the text field directly without using the file chooser
     */
    private class profileFileNameAction extends AbstractAction implements DocumentListener {

        profileFileNameAction() {
            super();
        }

        //any change that is made is handled the same - parameter is updated
        private void handleChange() {
            String value = profileFileNameTextField.getText();
            if (value != null) {
                File file = new File(value);
                currentParamsCopy.setProfileFileName(file);
            }
        }

        public void actionPerformed(ActionEvent e) {
            handleChange();
        }

        public void changedUpdate(DocumentEvent e) {
            handleChange();
        }

        public void insertUpdate(DocumentEvent e) {
            handleChange();
        }

        public void removeUpdate(DocumentEvent e) {
            handleChange();
        }
    }

    /**
     * Handles setting codon bias specification filename if the user edits the text field directly without using the file chooser
     */
    private class codonBiasFileNameAction extends AbstractAction implements DocumentListener {

        codonBiasFileNameAction() {
            super();
        }

        //any change that is made is handled the same - parameter is updated
        private void handleChange() {
            String value = biasFileNameTextField.getText();
            if (value != null) {
                File file = new File(value);
                currentParamsCopy.setCodonBiasFileName(file);
            }
        }

        public void actionPerformed(ActionEvent e) {
            handleChange();
        }

        public void changedUpdate(DocumentEvent e) {
            handleChange();
        }

        public void insertUpdate(DocumentEvent e) {
            handleChange();
        }

        public void removeUpdate(DocumentEvent e) {
            handleChange();
        }
    }

    /**
     * Handles setting for the text field parameters that are numbers (FormattedTextField's only).
     * Makes sure that the numbers make sense.
     */
    private class formattedTextFieldAction implements PropertyChangeListener {
        public void propertyChange(PropertyChangeEvent e) {
            Object source = e.getSource();
            if (source == lengthFormattedTextField) {
                Number value = (Number) lengthFormattedTextField.getValue();
                if ((value != null) && (value.intValue() > 0)) {
                    ProteinDatabaseSearchParams params = currentParamsCopy.getSearchParams();
                    if (params != null) {
                        params.setLength(value.intValue());
                    }
                }
            } else if (source == scoreThresholdFormattedTextField) {
                Number value = (Number) scoreThresholdFormattedTextField.getValue();
                if ((value != null) && (value.doubleValue() >= 0.0)) {
                    currentParamsCopy.setScoreThreshold(value.doubleValue());
                }
            } else if (source == topHitsFormattedTextField) {
                Number value = (Number) topHitsFormattedTextField.getValue();
                if ((value != null) && (value.intValue() > 0)) {
                    currentParamsCopy.setNumberTopHits(value.intValue());
                }
            } else if (source == fuzzFactorFormattedTextField) {
                Number value = (Number) fuzzFactorFormattedTextField.getValue();
                if ((value != null) && (value.doubleValue() >= 0.0)) {
                    currentParamsCopy.setFuzzFactor(value.doubleValue());
                }
            }
        }
    }

    /**
     * Handles choosing a node representation setting
     */
    private class selectNodeRepAction extends AbstractAction {
        selectNodeRepAction() {
            super();
        }

        public void actionPerformed(ActionEvent event) {
            if ("comboBoxChanged".equals(event.getActionCommand())) {
                JComboBox cb = (JComboBox) event.getSource();
                String nodeRepName = (String) cb.getSelectedItem();
                if ("proteins".equalsIgnoreCase(nodeRepName)) {
                    currentParamsCopy.setUniqueQueryProteinNodes(true);
                }
                else if ("domains".equalsIgnoreCase(nodeRepName)) {
                    currentParamsCopy.setUniqueQueryProteinNodes(false);
                }
            }
        }
    }

    

}
