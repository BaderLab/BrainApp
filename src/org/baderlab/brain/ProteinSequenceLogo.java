package org.baderlab.brain;

import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dp.WeightMatrix;
import org.biojava.bio.gui.DistributionLogo;
import org.biojava.bio.gui.SymbolStyle;
import org.biojava.bio.gui.TextLogoPainter;

import javax.swing.*;
import java.awt.*;
import java.awt.font.LineMetrics;
import java.awt.font.TextLayout;
import java.awt.geom.AffineTransform;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;

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
 * * Date: Oct 31, 2005
 * * Time: 9:57:40 PM
 */

/**
 * The main class to handle drawing of sequence logos for a ProteinProfile
 */
public class ProteinSequenceLogo {
    private ProteinProfile profile = null;

    private int sequenceLogoStartIndex = 1; //the start index for the sequence logo
    private boolean trimLogo = false; //required to trim logo to max informative region
    private double trimLogoPercentage = 0.1; //the default value to trim logos at (percentage of column height)
    private int logoHeight = 240; //the default logo height
    //maintains the detailed size of this logo
    private LogoSizeSpecification logoSpec = null;

    /**
     * Default constructor for a protein sequence logo. No profile trimming is performed.
     *
     * @param profile    The protein profile to draw
     * @param logoHeight The height, in pixels, to draw the logo
     */
    public ProteinSequenceLogo(ProteinProfile profile, int logoHeight) {
        this.profile = profile;
        this.trimLogo = false;
        logoSpec = new LogoSizeSpecification(this, profile, logoHeight);
    }

    /**
     * Constructor for a protein sequence logo that will trim columns on the end
     * of the profile.
     *
     * @param profile            The protein profile to draw
     * @param trimLogoPercentage If columns on either end of the profile are not above this height
     *                           (in terms of percentage of total possible column height), they will be trimmed.
     *                           This has the effect of highlighting the most informative contiguous region
     * @param logoHeight         The height, in pixels, to draw the logo
     */
    public ProteinSequenceLogo(ProteinProfile profile, double trimLogoPercentage, int logoHeight) {
        this.profile = profile;
        this.trimLogoPercentage = trimLogoPercentage;
        this.trimLogo = true;
        this.logoHeight = logoHeight;
        logoSpec = new LogoSizeSpecification(this, profile, this.logoHeight);
    }

    /**
     * Set the number of the start index for the sequence logo. This is the number that is drawn
     * under the leftmost sequence logo column position. The value can be positive or negative.
     */
    public void sequenceLogoSetStartIndex(int startIndex) {
        sequenceLogoStartIndex = startIndex;
    }

    /**
     * If true, will attempt to output only the region of the logo where leftmost and rightmost columns
     * contain more than the getTrimLogoPercentage percentage of total logo height
     */
    public boolean isTrimLogo() {
        return trimLogo;
    }

    /**
     * Only valid if isTrimLogo() returns true
     */
    public double getTrimLogoPercentage() {
        return trimLogoPercentage;
    }

    /**
     * Gets a detailed logo size specification of this logo
     */
    public LogoSizeSpecification getDetailedLogoSizeSpecification() {
        return logoSpec;
    }

    /**
     * Get the height of the logo
     */
    public int getLogoHeight() {
        return logoHeight;
    }

    /**
     * Get the width of the logo
     */
    public int getLogoWidth() {
        return ((int) logoSpec.logoWidth);
    }

    /**
     * Create a sequence logo for this profile. Residue coloring follows WebLogo (http://weblogo.berkeley.edu/)
     * This method is intended for bitmap graphics output.
     *
     * @return An image of the logo that can be saved e.g. ImageIO.write(image, "png", new File(outFileName));
     */
    public BufferedImage drawSequenceLogo() {
        BufferedImage bi = new BufferedImage((int) logoSpec.logoWidth, (int) logoSpec.logoHeight, BufferedImage.TYPE_INT_RGB);
        Graphics2D graphics = bi.createGraphics();
        graphics.setBackground(Color.WHITE);
        graphics.setColor(Color.WHITE);
        graphics.clearRect(0, 0, (int) logoSpec.logoWidth, (int) logoSpec.logoHeight);
        drawSequenceLogo(graphics, null);
        return (bi);
    }

    /**
     * Handles drawing of the logo axes
     */
    private class LogoAxisPanel extends JComponent {
        //size variables - set to defaults - should be overridden with preferred sizes
        private LogoSizeSpecification logoSpec = null;
        private String title = null;

        public void setLogoSizeSpecification(LogoSizeSpecification logoSpec) {
            this.logoSpec = logoSpec;
        }

        public void setTitle(String title) {
            this.title = title;
        }

        protected void paintComponent(Graphics gOrig) {
            Graphics2D g = (Graphics2D) gOrig;

            int logoWidth = (int) logoSpec.logoWidth;
            int logoHeight = (int) logoSpec.logoHeight;
            int padLeft = (int) logoSpec.padLeft;
            int padBottom = (int) logoSpec.padBottom;
            int columnWidth = (int) logoSpec.columnWidth;
            int columnHeight = (int) logoSpec.columnHeight;

            //add extra information around logo
            Dimension panelSize = new Dimension(logoWidth, logoHeight);
            //draw Y axis
            g.setColor(Color.BLACK);
            g.setFont(new Font("Dialog", Font.BOLD, (int) logoSpec.fontSize));
            int yAxisOffset = 3;
            int yAxisTickLength = (int) logoSpec.yAxisTickLength;
            int tickLabelPadding = 2;
            int labelWidth = 0;
            g.drawLine(padLeft - yAxisOffset, 0, padLeft - yAxisOffset, columnHeight);
            double maxBits = (Math.log(logoSpec.sequenceAlphabetSize) / Math.log(2.0));
            for (int i = 0; i < maxBits; i++) {
                int y = (int) (columnHeight - (i * columnHeight / maxBits));
                g.drawLine(padLeft - yAxisTickLength, y, padLeft - yAxisOffset, y);
                FontMetrics f = g.getFontMetrics();
                String tickLabel = Integer.toString(i);
                int tickLabelXOffset = (int) (f.getStringBounds(tickLabel, g)).getWidth();
                labelWidth = tickLabelXOffset;
                LineMetrics l = f.getLineMetrics(tickLabel, g);
                int tickLabelYOffset = (int) Math.floor(l.getAscent() / 2);
                g.drawString(tickLabel, padLeft - yAxisTickLength - tickLabelXOffset - tickLabelPadding, y + tickLabelYOffset);
            }
            //draw y axis label
            int columnLabelPadding = 2;
            String yAxisLabel = "bits";
            int yAxisLabelXPos = padLeft - yAxisTickLength - labelWidth - tickLabelPadding * 2;
            int yAxisLabelYPos = columnHeight / 2;
            drawRotatedString(yAxisLabel, g, yAxisLabelXPos, yAxisLabelYPos);
            //draw x axis numbering
            for (int i = logoSpec.minColumnIndex; i <= logoSpec.maxColumnIndex; i++) {
                String columnLabel = Integer.toString(sequenceLogoStartIndex + i);
                int x = padLeft + (columnWidth / 2) + ((i - logoSpec.minColumnIndex) * columnWidth); //x coord of the center of the string
                int y = columnHeight + columnLabelPadding; //y coord of the top of the string
                FontMetrics f = g.getFontMetrics();
                Rectangle2D stringBounds = f.getStringBounds(columnLabel, g);
                LineMetrics l = f.getLineMetrics(columnLabel, g);
                if (stringBounds.getWidth() >= columnWidth) {
                    //string is too wide, so draw it sideways
                    //see if we need to scale the text
                    if (stringBounds.getWidth() >= padBottom) {
                        //scale to make it fit
                        AffineTransform startMatrix = g.getTransform();
                        double scaleFactor = ((double) padBottom) / stringBounds.getWidth();
                        g.scale(scaleFactor, scaleFactor);
                        //center coords
                        y += (stringBounds.getWidth() / 2) * scaleFactor;
                        x += (l.getAscent() / 2) * scaleFactor;
                        drawRotatedString(columnLabel, g, (int) (x / scaleFactor), (int) (y / scaleFactor));
                        //restore the transform
                        g.setTransform(startMatrix);
                    } else {
                        //center coords
                        y += stringBounds.getWidth() / 2;
                        x += (l.getAscent() / 2);
                        drawRotatedString(columnLabel, g, x, y);
                    }
                } else {
                    //string can be drawn horizontally (normally)
                    //center coords
                    y += l.getAscent();
                    x -= (stringBounds.getWidth() / 2);
                    g.drawString(columnLabel, x, y);
                }
            }

            //optionally add information around logo
            if (title != null) {
                int padTitle = 3;
                FontMetrics f = g.getFontMetrics();
                Rectangle2D stringBounds = f.getStringBounds(title, g);
                g.drawString(title, (int) ((panelSize.getWidth() / 2) - (stringBounds.getWidth() / 2)), (int) panelSize.getHeight() - padTitle);
            }
        }
    }

    /**
     * Draws a sequence logo on a graphics object centered at a specific point. This method is
     * intended for vector graphics output. Uses the default logo color style (from Berkeley Weblogo)
     *
     * @param g      The (preferably vector) graphics object to paint to.
     * @param center Where to center the logo drawing
     */
    public void drawSequenceLogo(Graphics2D g, Point center) {
        drawSequenceLogo(g, center, new WebLogoProteinStyle());
    }

    /**
     * Draws a sequence logo on a graphics object centered at a specific point. This method is
     * intended for vector graphics output.
     *
     * @param g                The (preferably vector) graphics object to paint to.
     * @param center           Where to center the logo drawing
     * @param symbolColorStyle Specify the logo color style to use. You must implement the BioJava SymbolStyle
     *                         to create your own color scheme.
     * @see WebLogoProteinStyle
     */
    public void drawSequenceLogo(Graphics2D g, Point center, SymbolStyle symbolColorStyle) {
        int logoWidth = (int) logoSpec.logoWidth;
        int logoHeight = (int) logoSpec.logoHeight;
        int columnWidth = (int) logoSpec.columnWidth;
        int columnHeight = (int) logoSpec.columnHeight;

        WeightMatrix wm = profile.getWeightMatrix();

        JPanel logoPanel = new JPanel(new GridLayout(1, logoSpec.numberOfColumns));
        Dimension logoPanelSize = new Dimension((int) columnWidth * logoSpec.numberOfColumns, (int) columnHeight);
        logoPanel.setPreferredSize(logoPanelSize);
        logoPanel.setOpaque(false);
        RenderingHints hints = new RenderingHints(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

        try {
            //build sequence logo for all columns on a panel
            for (int pos = logoSpec.minColumnIndex; pos <= logoSpec.maxColumnIndex; pos++) {
                Distribution dist = wm.getColumn(pos);
                DistributionLogo dl = new DistributionLogo();
                dl.setRenderingHints(hints);
                dl.setOpaque(false);
                dl.setDistribution(dist);
                dl.setPreferredSize(new Dimension((int) columnWidth, (int) columnHeight));
                dl.setLogoPainter(new TextLogoPainter());
                dl.setStyle(symbolColorStyle);
                logoPanel.add(dl);
            }
        } catch (BioException ex) {
            throw new BioError(ex);
        }

        //Create axis picture
        LogoAxisPanel axisPanel = new LogoAxisPanel();
        axisPanel.setLogoSizeSpecification(logoSpec);
        axisPanel.setTitle(profile.getName() + " (" + profile.getNumSequences() + " peptides)");

        Dimension axisPanelSize = new Dimension(logoWidth, logoHeight);
        axisPanel.setPreferredSize(axisPanelSize);
        axisPanel.setOpaque(false);

        JPanel panel = new JPanel();
        panel.setOpaque(false);
        panel.setPreferredSize(axisPanelSize);
        panel.setLayout(null); //we'll handle the layout details
        axisPanel.setBounds(0, 0, (int) axisPanelSize.getWidth(), (int) axisPanelSize.getHeight());
        panel.add(axisPanel);
        logoPanel.setBounds((int) (axisPanelSize.getWidth() - logoPanelSize.getWidth()), 0, (int) logoPanelSize.getWidth(), (int) logoPanelSize.getHeight());
        panel.add(logoPanel);

        //position the graphic in the right place
        AffineTransform at = null;
        if (center != null) {
            at = g.getTransform(); //remember the transform
            g.translate(center.getX() - (logoWidth / 2), center.getY() - (logoHeight / 2));
        }

        //make sure all components are packed and print to the Graphics object
        JFrame frame = new JFrame();
        frame.getContentPane().add(panel);
        frame.pack();
        panel.print(g);
        frame.dispose();

        if (at != null) {
            g.setTransform(at);
        }
    }

    /**
     * Draws rotated text.
     *
     * @param string The string to draw
     * @param g      Draw the string on this graphics object
     * @param x      The x coordinate of the center of the string
     * @param y      The y coordinate of the center of the string
     */
    private void drawRotatedString(String string, Graphics2D g, int x, int y) {
        TextLayout layout = new TextLayout(string, g.getFont(), g.getFontRenderContext());
        //remember the current transform
        AffineTransform startMatrix = g.getTransform();
        g.translate(x, y);
        g.rotate(-Math.PI / 2);
        layout.draw(g, -layout.getAdvance() / 2, 0);
        //restore the transform
        g.setTransform(startMatrix);
    }

}
