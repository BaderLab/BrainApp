package org.baderlab.brain;

import org.biojava.bio.gui.SymbolStyle;
import org.freehep.graphics2d.VectorGraphics;
import org.freehep.graphicsio.PageConstants;
import org.freehep.graphicsio.pdf.PDFGraphics2D;
import org.freehep.graphicsio.ps.PSGraphics2D;

import java.awt.*;
import java.awt.font.LineMetrics;
import java.awt.geom.AffineTransform;
import java.awt.geom.Rectangle2D;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.List;
import java.util.Properties;

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
 * * Date: Aug 24, 2005
 * * Time: 3:00:24 PM
 */

/**
 * Draw a LogoTree from hierarchical clustering and profile information
 */
public class LogoTreeDraw {
    private HierarchicalClusteringResultTree root = null; //the root of the tree
    private int totalNumberOfLeaves = 0;  //total leaves in the clustering result
    private double maxDistance = 0.0; //max distance in the clustering result
    private List profileList = null;
    private List nodeProfileList = null;
    private AvgLinkHierarchicalClustering cluster = null;  //the original clustering result
    private int[] leafOrder;
    private int sequenceLogoStartIndex = 1;
    private ArrayList labelHighlight; //a list of LabelColorPair objects for highlighting leaf labels
    //optional field to store any cluster bootstrap results. if this object is set, the bootstrap results will be drawn on the tree
    private HierarchicalClusteringBootstrapResult bootstrapResults = null;
    private SymbolStyle symbolStyle; //optional field for logo symbol style (letter coloring)

    //page sizing parameters
    public static int LETTER = 1;
    public static int LEGAL = 2;
    //width and height of the page
    private int width;
    private int height;
    //pad the left and right side
    private final double leftPadFactor = 0.13; //% padding on left
    private final double rightPadFactor = 0.13; //% padding on right
    //sizing parameters for logos
    private final double logoScaleFactorPercentage = 0.13;   //node logos scale to hit a % of the height at the root
    private double logoScaleFactor;  //the scale to maximally size logos on the left side of the tree
    //other variables calculated in constructor based on the above parameters
    private int baseLogoHeight;  //actual height of the smallest logo
    private int leftPad;  //actual pixel number to pad left
    private int rightPad;  //actual pixel number to pad right
    private int topPad = 0; //pixel padding for top - only set if there is a title
    private String title = null; //title of the output

    //trim logo fields
    private boolean trimLeafLogo = false; //required to trim leaf logo to max informative region
    private boolean trimNodeLogo = false; //required to trim leaf logo to max informative region
    private double trimLeafLogoPercentage = 0.1; //the default value to trim leaf logos at (percentage of column height)
    private double trimNodeLogoPercentage = 0.1; //the default value to trim node logos at (percentage of column height)

    //cut node profile fields (these are used for options that determine the length of logos on tree nodes only)
    private int nodeLogoProfileLength = 0;
    private ProteinTerminus nodeLogoProfileTerminus = null;

    //the maximum right margin x coordinate touched by any internal node tree logo drawn
    private int rightMarginXCoord;

    /**
     * Constructor for a LogoTree
     *
     * @param cluster     The clustering result to draw
     * @param profileList The list of profiles that were clustered
     */
    public LogoTreeDraw(AvgLinkHierarchicalClustering cluster, List profileList) {
        root = cluster.getResult();
        totalNumberOfLeaves = cluster.getNelements();
        maxDistance = cluster.getMaxDistance();
        leafOrder = cluster.getLeafOrder();
        this.profileList = profileList;
        this.cluster = cluster;
        nodeProfileList = createClusterProfileList(profileList);

        //page size - any page size can be added here - just make sure the aspect ratio is correct
        //and keep it at a high resolution, like 1250dpi
        //these numbers should be pretty high to have high resolution
        setSize(10625, 13750);   //this equals 1250 dpi for letter size paper (8.5x11 inches)
    }

    /**
     * Internally sets the width and height and makes sure other variables that are
     * based on width and height are correctly updated.
     *
     * @param width  Width of page
     * @param height Height of page
     */
    private void setSize(int width, int height) {
        this.width = width;
        this.height = height;
        //logo sizing calculations
        //determine size for placing logos on right side
        //(leaf logos all in a column on the right side + some padding between them - based on % of height)
        baseLogoHeight = (int) ((height + (height * 0.05)) / totalNumberOfLeaves);

        leftPad = (int) (width * leftPadFactor);
        rightPad = (int) (width * rightPadFactor);
        int maxLogoWidth = getMaxLogoWidth();
        if (maxLogoWidth > rightPad) {
            rightPad = (int) (maxLogoWidth * 2); //make the right margin larger if the logos won't fit
            //this still doesn't guarantee the logos will fit, since the margin can be modified by
            //the node logos.
        }

        logoScaleFactor = height * logoScaleFactorPercentage / baseLogoHeight;

        rightMarginXCoord = width - rightPad; //set the normal right margin - may be updated during tree drawing
    }

    /**
     * Draw the logo tree to a graphics object
     *
     * @param g The graphics object to drawn on
     */
    public void draw(Graphics2D g) {
        //clear the area
        g.setBackground(Color.WHITE);
        g.clearRect(0, 0, width, height);
        g.setColor(Color.BLACK);
        if (title != null) {
            //draw the title
            g.setColor(Color.BLACK);
            //calculate center of string
            //set the upper limit on string size, then scale down if required
            g.setFont(new Font("Dialog", Font.BOLD, baseLogoHeight));
            FontMetrics f = g.getFontMetrics();
            double titleWidth = f.getStringBounds(title, g).getWidth();
            LineMetrics l = f.getLineMetrics(title, g);

            //we have width space to draw the string
            double fontScaleFactor = 1.0;
            int availableStringSpace = width;
            if (titleWidth > availableStringSpace) {
                //scale the font if required to fill up width
                fontScaleFactor = availableStringSpace / titleWidth;
                g.setFont(new Font("Dialog", Font.BOLD, (int) (baseLogoHeight * fontScaleFactor)));
            }
            g.drawString(title, (int) ((availableStringSpace - titleWidth * fontScaleFactor) / 2), (int) l.getAscent());

            //set the top padding, so rest of code knows how high title was
            topPad = (int) f.getStringBounds(title, g).getHeight();
            //adjust height accordingly
            height -= topPad;
            //recalculate from above
            baseLogoHeight = (int) ((height + (height * 0.05)) / totalNumberOfLeaves);
        }
        //draw nodes recursively
        drawNode(root, g);
        //draw leaves iteratively
        drawLeaves(g);
    }

    /**
     * Output the logo tree to a PDF file
     *
     * @param outputFile The filename to output to e.g. "mylogotree.pdf"
     */
    public void outputToPDF(File outputFile) {
        Properties p = new Properties();
        //this doesn't seem to actually do anything for PDF - TODO: look into this further sometime
        p.setProperty(PageConstants.PAGE_SIZE, PageConstants.LETTER);
        VectorGraphics g = null;
        try {
            g = new PDFGraphics2D(outputFile, new Dimension(width, height));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        g.setProperties(p);
        g.startExport();
        this.draw(g);
        g.endExport();
    }

    /**
     * Output the logo tree to a PostScript file
     *
     * @param outputFile The filename to output to e.g. "mylogotree.pdf"
     * @param pageSize   Currently one of LogoTreeDraw.LETTER or LogoTreeDraw.LEGAL
     */
    public void outputToPS(File outputFile, int pageSize) {
        //page size - any page size can be added here - just make sure the aspect ratio is correct
        //and keep it at a high resolution, like 1250dpi
        if (pageSize == LogoTreeDraw.LETTER) {
            //these numbers should be pretty high to have high resolution - the page size will always be letter
            width = 10625;   //this equals 1250 dpi for letter size paper (8.5x11 inches)
            height = 13750;
        } else if (pageSize == LogoTreeDraw.LEGAL) {
            width = 10625;   //this equals 1250 dpi for legal size paper (8.5X14 inches)
            height = 17500;
        }
        Properties p = new Properties();
        p.setProperty(PageConstants.PAGE_SIZE, PageConstants.LETTER);
        VectorGraphics g = null;
        try {
            g = new PSGraphics2D(outputFile, new Dimension(width, height));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        g.setProperties(p);
        g.startExport();
        this.draw(g);
        g.endExport();
    }

    /**
     * If true, will attempt to output only the region of the leaf logo where leftmost and rightmost columns
     * contain more than the specified percentage of total logo height
     *
     * @param trimLeafLogo
     */
    public void setTrimLeafLogo(boolean trimLeafLogo, double percentageLogoHeight) {
        this.trimLeafLogo = trimLeafLogo;
        this.trimLeafLogoPercentage = percentageLogoHeight;
    }

    /**
     * If true, will attempt to output only the region of the node logo where leftmost and rightmost columns
     * contain more than the specified percentage of total logo height
     *
     * @param trimNodeLogo
     */
    public void setTrimNodeLogo(boolean trimNodeLogo, double percentageLogoHeight) {
        this.trimNodeLogo = trimNodeLogo;
        this.trimNodeLogoPercentage = percentageLogoHeight;
    }

    /**
     * Set a title of the output document. If null, no title will be drawn.
     * Default is null.
     */
    public void setTitle(String title) {
        this.title = title;
    }

    /**
     * Set the number of the start index for the sequence logo. This is the number that is drawn
     * under the leftmost sequence logo column position. The value can be positive or negative.
     */
    public void setSequenceLogoStartIndex(int sequenceLogoStartIndex) {
        this.sequenceLogoStartIndex = sequenceLogoStartIndex;
    }

    /**
     * Optional field to store any cluster bootstrap results.
     * If this object is set, the bootstrap results will be drawn on the tree
     *
     * @param bootstrapResults The boostrap results to draw on the tree
     */
    public void setBootstrapResults(HierarchicalClusteringBootstrapResult bootstrapResults) {
        this.bootstrapResults = bootstrapResults;
    }

    /**
     * Optional field to store the symbol style for logo coloring. If not set, the default is to use the
     * WebLogo symbol style.
     *
     * @param symbolStyle The symbol style to set
     */
    public void setSymbolStyle(SymbolStyle symbolStyle) {
        this.symbolStyle = symbolStyle;
    }

    /**
     * Option to create a subset of the logo only on nodes. You may want to use this to show that the logo
     * tree clustering was only done on a subset of the profile, but you want to show the full profile on the
     * leaves.
     *
     * @param profileLength The length to cut the profile to
     * @param terminus      The terminus to cut from
     */
    public void setNodeLogoSubset(int profileLength, ProteinTerminus terminus) {
        this.nodeLogoProfileLength = profileLength;
        this.nodeLogoProfileTerminus = terminus;
    }

    /**
     * Used to store label list / color pairs for highlighting leaf labels
     */
    private class LabelColorPair {
        private Color color;
        private ArrayList labels;

        public LabelColorPair(Color color, ArrayList labels) {
            this.color = color;
            this.labels = labels;
        }

        public Color getColor() {
            return color;
        }

        public void setColor(Color color) {
            this.color = color;
        }

        public ArrayList getLabels() {
            return labels;
        }

        public void setLabels(ArrayList labels) {
            this.labels = labels;
        }
    }

    /**
     * Get the leaf highlight color, if defined
     *
     * @param label The leaf label
     * @return A color if the highlight color is defined, otherwise null.
     */
    private Color getLeafLabelHighlightColor(String label) {
        for (int j = 0; j < labelHighlight.size(); j++) {
            LabelColorPair labelColorPair = (LabelColorPair) labelHighlight.get(j);
            ArrayList labelsToHighlight = labelColorPair.getLabels();
            if (labelsToHighlight.contains(label)) {
                return (labelColorPair.getColor());
            }
        }
        //not found
        return (null);
    }

    /**
     * Set a list of leaf labels to highlight
     * Labels should not overlap previously set labels
     * If a single label name is set to be highlighted in two or more colors, the behavior is undefined.
     *
     * @param labelsToHighlight The list of label Strings to highlight
     * @param color             The color to highlight
     *                          Background color will always be white.
     */
    public void setLeafLabelHighlightColor(ArrayList labelsToHighlight, Color color) {
        if (labelHighlight == null) {
            labelHighlight = new ArrayList();
        }
        LabelColorPair lcp = new LabelColorPair(color, labelsToHighlight);
        labelHighlight.add(lcp);
    }

    /**
     * Recursively draw the tree
     *
     * @param node Draw recursively down from this node
     * @param g    Draw onto this graphics object
     */
    private void drawNode(HierarchicalClusteringResultTree node, Graphics2D g) {
        //base case - don't do anything right now
        if (node.leaf) {
            return;
        }
        //each call draws one branch point - looks like '['
        int leftY = (int) ((node.left.leafOrderIndex + 0.5) / totalNumberOfLeaves * height);
        leftY += topPad;
        int rightY = (int) ((node.right.leafOrderIndex + 0.5) / totalNumberOfLeaves * height);
        rightY += topPad;
        int topX = (int) (getDistanceScale(node.distance) * (width - leftPad - rightPad) + leftPad);
        int bottomLeftX = (int) (getDistanceScale(node.left.distance) * (width - leftPad - rightPad) + leftPad);
        int bottomRightX = (int) (getDistanceScale(node.right.distance) * (width - leftPad - rightPad) + leftPad);
        Stroke s = g.getStroke(); //remember old stroke
        g.setStroke(new BasicStroke(height / 1000));
        g.drawPolyline(new int[]{bottomLeftX, topX, topX, bottomRightX}, new int[]{leftY, leftY, rightY, rightY}, 4);
        g.setStroke(s); //reset old stroke

        //draw the logo at a specific size and position
        ProteinProfile nodeProfile = (ProteinProfile) nodeProfileList.get(node.nodeIndex - 1);
        Point center = new Point(topX, (leftY + rightY) / 2);
        //scale logo size to be larger near tree root (based on cluster size)
        //scaleFactor ranges from 1.0 to logoScaleFactor
        double scaleFactor = logoScaleFactor - (getClusterSizeScale(node.numberOfLeaves) * (logoScaleFactor - 1));
        ProteinSequenceLogo logo = null;
        //first, evaluate node logo subset options
        if (nodeLogoProfileLength > 0) {
            int originalProfileLegth = nodeProfile.getNumColumns();
            nodeProfile = nodeProfile.getTruncatedProfileCopy(nodeLogoProfileLength, nodeLogoProfileTerminus);
            //correct sequenceLogoStartIndex, if required
            if (nodeLogoProfileTerminus.equals(ProteinTerminus.C)) {
                sequenceLogoStartIndex = sequenceLogoStartIndex + (originalProfileLegth - nodeLogoProfileLength);
            }
        }
        //second, evaluate node logo trim options (trim could still be used with subset)
        if (trimNodeLogo) {
            logo = new ProteinSequenceLogo(nodeProfile, trimNodeLogoPercentage, (int) (baseLogoHeight * scaleFactor));
        } else {
            logo = new ProteinSequenceLogo(nodeProfile, (int) (baseLogoHeight * scaleFactor));
        }
        logo.sequenceLogoSetStartIndex(sequenceLogoStartIndex);
        if (symbolStyle != null) {
            logo.drawSequenceLogo(g, center, symbolStyle);
        } else {
            logo.drawSequenceLogo(g, center);
        }
        int logoWidth = logo.getLogoWidth();
        //check if the logo extended past the right margin
        if ((topX + (logoWidth / 2)) > (width - rightPad)) {
            //move the right margin over to the right to fit this logo
            rightMarginXCoord = Math.max(rightMarginXCoord, topX + (logoWidth / 2));
        }

        //optionally draw the bootstrap results
        if (bootstrapResults != null) {
            int bootstrapCount = bootstrapResults.getBootstrapCount(node);
            Color origColor = g.getColor();
            g.setColor(Color.red);
            Font origFont = g.getFont();
            //calculate a reasonable font height
            int fontHeight = logo.getLogoHeight() / 2;
            if (fontHeight > (height / 80)) { //80 = heuristic
                //make sure font height is not too large
                fontHeight = (int) (height / 80);
            }
            g.setFont(new Font("Dialog", Font.BOLD, fontHeight));
            g.drawString(Integer.toString(bootstrapCount), center.x - (int) (logoWidth / 2), center.y - (int) (logo.getLogoHeight() / 2));
            g.setColor(origColor);
            g.setFont(origFont);
        }

        //recursive call
        drawNode(node.left, g);
        drawNode(node.right, g);
    }

    private void drawHighlightRectangle(String label, int x, int y, Graphics2D g) {
        Color c = getLeafLabelHighlightColor(label);
        if (c != null) {
            Color oldColor = g.getColor();
            g.setColor(c);
            FontMetrics fm = g.getFontMetrics();
            Rectangle2D r = fm.getStringBounds(label, g);
            g.fillRect(x, y - ((int) fm.getMaxAscent()), (int) r.getWidth() + 1, (int) fm.getMaxAscent() + 1);
            g.setColor(oldColor);
        }
    }

    /**
     * Draw the leaves on the right side of the tree
     *
     * @param g The graphics object to draw on
     */
    private void drawLeaves(Graphics2D g) {
        //find the widest logo - used to align text
        int maxLogoWidth = getMaxLogoWidth();
        maxLogoWidth += maxLogoWidth * 0.05; //padding for labels

        for (int i = 0; i < leafOrder.length; i++) {
            int leafIndex = leafOrder[i];
            //draw the logo of the actual profile
            ProteinProfile leafProfile = (ProteinProfile) profileList.get(leafIndex);
            ProteinSequenceLogo logo = null;
            if (trimLeafLogo) {
                logo = new ProteinSequenceLogo(leafProfile, trimLeafLogoPercentage, baseLogoHeight);
            } else {
                logo = new ProteinSequenceLogo(leafProfile, baseLogoHeight);
            }
            logo.sequenceLogoSetStartIndex(sequenceLogoStartIndex);
            int logoWidth = logo.getLogoWidth();
            //draw the logo justified as it would be if it were full length
            //take into account the number of columns that were filtered from left and right side (if any) because of setTrimLeafLogo
            LogoSizeSpecification logoSize = logo.getDetailedLogoSizeSpecification();
            int pushRight = (int) (logoSize.minColumnIndex * logoSize.columnWidth);
            Point center = new Point(rightMarginXCoord + (logoWidth) / 2 + pushRight, (int) ((i + 0.5) / totalNumberOfLeaves * height) + topPad);
            if (symbolStyle != null) {
                logo.drawSequenceLogo(g, center, symbolStyle);
            } else {
                logo.drawSequenceLogo(g, center);
            }

            //draw the profile name
            String profileName = leafProfile.getName();
            g.setColor(Color.BLACK);
            //calculate center of string
            //set the lower limit on string size, then scale up later
            g.setFont(new Font("Dialog", Font.BOLD, 12));
            FontMetrics fm = g.getFontMetrics();
            Rectangle2D stringBounds = fm.getStringBounds(profileName, g);
            //we have width - rightMarginXCoord space to draw the string
            double fontScaleFactor = 1.0;
            //add in a padding width between logo and string
            //int padBoundsWidth = (int) (logoWidth*0.1); //heuristically, 10% of logo width
            int availableStringWidth = width - rightMarginXCoord - maxLogoWidth;
            if (stringBounds.getWidth() < availableStringWidth) {
                //scale the font up if required to fill up availableStringWidth
                fontScaleFactor = availableStringWidth / stringBounds.getWidth();
                //make sure the font is not too high
                if (stringBounds.getHeight() * fontScaleFactor > logoSize.logoHeight) {
                    //reset the scaling factor by height
                    fontScaleFactor = logoSize.logoHeight / stringBounds.getHeight();
                }
            } //assume that we don't have to scale down, since the font starts off very small
            //remember the current transform
            AffineTransform startMatrix = g.getTransform();
            g.scale(fontScaleFactor, fontScaleFactor); //scale up so string fits in space
            int xRightJustified = rightMarginXCoord + maxLogoWidth;
            int yCenter = (int) ((i + 0.5) / totalNumberOfLeaves * height) + topPad;
            LineMetrics l = fm.getLineMetrics(profileName, g);
            yCenter += (l.getAscent() / 2) * fontScaleFactor;
            int x = (int) (xRightJustified / fontScaleFactor);
            int y = (int) (yCenter / fontScaleFactor);
            if (labelHighlight != null) {
                drawHighlightRectangle(profileName, x, y, g);
            }
            g.drawString(profileName, x, y);
            //restore the transform
            g.setTransform(startMatrix);
        }
    }

    /**
     * Helper method to get the maximum logo width from profiles only (not nodes)
     * Depends on baseLogoHeight being set
     */
    private int getMaxLogoWidth() {
        int maxLogoWidth = 0;
        for (int i = 0; i < profileList.size(); i++) {
            ProteinProfile proteinProfile = (ProteinProfile) profileList.get(i);
            ProteinSequenceLogo logo = null;
            if (trimLeafLogo) {
                logo = new ProteinSequenceLogo(proteinProfile, trimLeafLogoPercentage, baseLogoHeight);
            } else {
                logo = new ProteinSequenceLogo(proteinProfile, baseLogoHeight);
            }
            maxLogoWidth = Math.max(maxLogoWidth, logo.getLogoWidth());
        }
        return maxLogoWidth;
    }

    //helper method to determine scale by node joining distance
    private double getDistanceScale(double distance) {
        return ((maxDistance - distance) / maxDistance);
    }

    //helper method to determine logo scale by cluster size
    private double getClusterSizeScale(int numberOfLeavesAtThisNode) {
        return ((double) (totalNumberOfLeaves - numberOfLeavesAtThisNode) / (double) totalNumberOfLeaves);
    }

    /**
     * Create the profiles describing each internal node of the cluster tree
     *
     * @param profileList The original profile list that was clustered
     * @return The list of comboProfiles, one comboProfile per internal node. The list is
     *         indexed by the node number-1 i.e. from zero to (max node number -1). The cluster results number
     *         nodes from 1, so you must do the conversion (node number -1) to access the proper combo profile
     *         from this list.
     */
    private ArrayList createClusterProfileList(List profileList) {
        ArrayList comboProfileList = new ArrayList(profileList.size());
        //initialize comboProfileList to be a bunch of empty lists
        for (int i = 0; i < profileList.size(); i++) {
            comboProfileList.add(new String("null")); //dummy object
        }
        HierarchicalClusteringResultTree clusterResult = cluster.getResult();
        //recursively traverse tree to get the combo profiles
        collectComboProfiles(clusterResult, profileList, comboProfileList);
        return (comboProfileList);
    }

    /**
     * Recursive helper method for createClusterProfileList
     *
     * @param node             The node to act on
     * @param profileList      The original list of profiles used for the clustering
     * @param comboProfileList The list to be populated with combo profiles
     */
    private void collectComboProfiles(HierarchicalClusteringResultTree node, List profileList, ArrayList comboProfileList) {
        if (node.leaf == false) {
            //only work on internal nodes
            //at each non-leaf node, collect the leaves
            List childrenLeaves = new ArrayList(node.flatLeftChildrenLeaves);
            childrenLeaves.addAll(node.flatRightChildrenLeaves);
            //convert to profiles
            ArrayList profileListForNode = new ArrayList();
            for (int i = 0; i < childrenLeaves.size(); i++) {
                HierarchicalClusteringResultTree leaf = (HierarchicalClusteringResultTree) childrenLeaves.get(i);
                profileListForNode.add(profileList.get(leaf.nodeIndex));
            }
            //create a combo profile and add to the combo profile list
            ProteinProfile comboProfile = new ProteinProfile(profileListForNode, "Node" + node.nodeIndex);
            comboProfileList.set(node.nodeIndex - 1, comboProfile); //-1 because other logo code expects comboProfileList to start at 0
            collectComboProfiles(node.left, profileList, comboProfileList);
            collectComboProfiles(node.right, profileList, comboProfileList);
        }
        //ignore leaf nodes
    }
}
