package org.baderlab.csplugins.brain;

import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.DistributionTools;
import org.biojava.bio.dp.WeightMatrix;

/**
 * Stores the size specification for a logo
 */
public class LogoSizeSpecification {
    //these are doubles just in case we need to use them later e.g. for vector graphics drawing
    public double columnWidth;
    public double columnHeight;
    public double padLeft;
    public double padBottom;
    public double logoWidth;
    public double logoHeight;
    public double yAxisTickLength;
    public double fontSize;

    //Logo column counts
    //indices are from 0..profileWidth-1
    public int minColumnIndex; //may be different than full profile if we are trimming the logo
    public int maxColumnIndex; //may be different than full profile if we are trimming the logo
    public int numberOfColumns;

    public int sequenceAlphabetSize;

    /**
     * Calculates the dimension of a logo, given the currently set options e.g. trimLogo, logoHeight
     */
    public LogoSizeSpecification(ProteinSequenceLogo proteinSequenceLogo, ProteinProfile profile, int logoHeight) {
        //logo column is 3 wide by 10 high aspect ratio
        double columnAspectHeight = 10.0;
        double columnAspectWidth = 3.0;
        //bottom padding is 2, left padding is 2 in these units
        double leftPadAspect = 2.0;
        double bottomPadAspect = 2.0;
        //other
        double yAxisTickLengthAspect = 0.5;
        double fontAspect = 0.75;

        double unit = 0.0;

        //see if we need to trim the logo
        WeightMatrix wm = profile.getWeightMatrix();
        minColumnIndex = Integer.MAX_VALUE;
        maxColumnIndex = Integer.MIN_VALUE;
        double bits = Math.log(20.0) / Math.log(2.0);
        if (proteinSequenceLogo.isTrimLogo()) {
            //calculate the columns to draw (interval where boundary columns
            //have more information than trimLogoPercentage)
            for (int pos = 0; pos < wm.columns(); pos++) {
                Distribution dist = wm.getColumn(pos);
                double informationBits = DistributionTools.bitsOfInformation(dist);
                double informationPercentage = informationBits / bits;
                if (proteinSequenceLogo.getTrimLogoPercentage() < informationPercentage) {
                    minColumnIndex = Math.min(minColumnIndex, pos);
                    maxColumnIndex = Math.max(maxColumnIndex, pos);
                }
            }
        } else {
            //draw the whole logo
            minColumnIndex = 0;
            maxColumnIndex = wm.columns() - 1;  //minus one to convert to array index
        }
        numberOfColumns = maxColumnIndex - minColumnIndex + 1;

        columnHeight = ((double) logoHeight) / (columnAspectHeight + bottomPadAspect) * columnAspectHeight;
        unit = columnHeight / columnAspectHeight;
        columnWidth = unit * columnAspectWidth;

        padLeft = unit * leftPadAspect;
        padBottom = unit * bottomPadAspect;

        logoWidth = columnWidth * numberOfColumns + padLeft;
        this.logoHeight = logoHeight;

        yAxisTickLength = unit * yAxisTickLengthAspect;
        fontSize = unit * fontAspect;

        sequenceAlphabetSize = 20; //amino acid number
    }
}
