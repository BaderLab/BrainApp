package org.baderlab.csplugins.brain;

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
 * * Date: Feb 5, 2005
 * * Time: 4:38:09 PM
 * * Description: Annotates a sequence with hits to a protein profile (uses a normalized p-value threshold)
 */

import org.biojava.bio.BioException;
import org.biojava.bio.SimpleAnnotation;
import org.biojava.bio.dp.DP;
import org.biojava.bio.dp.ScoreType;
import org.biojava.bio.dp.WeightMatrix;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceAnnotator;
import org.biojava.bio.seq.impl.ViewSequence;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.RangeLocation;
import org.biojava.utils.ChangeVetoException;

import java.io.Serializable;

/**
 * Annotates a sequence with hits to a ProteinProfile. Uses a normalized p-value threshold.
 * Modified from biojava.org.bio.dp.WeightMatrixAnnotator.java
 * <p/>
 * This SequenceAnnotator implementation returns a new
 * ViewSequence wrapping the underlying Sequence
 * </p>
 *
 * @see ProteinProfile#getNormalizedPValue(double)
 */
public class NormalizedProfileAnnotator implements SequenceAnnotator,
        Serializable {
    private ProteinProfile profile;
    private WeightMatrix matrix;
    private double threshold;
    private final ScoreType scoreType;

    public Sequence annotate(Sequence seq) throws IllegalAlphabetException,
            BioException, ChangeVetoException {
        seq = new ViewSequence(seq);

        int cols = matrix.columns();
        Feature.Template template = new Feature.Template();
        template.source = "NormalizedProfileAnnotator";
        for (int offset = 1;
             offset <= seq.length() - cols + 1;
             offset++) {
            /*Note: sequences containing X will generate all possible symbols and will get an additive score
            It should be a score equal to the maximum in the column, not additive over all symbols according to
            our normalization model.*/
            double score = DP.scoreWeightMatrix(matrix, seq, scoreType, offset);
            double qraw = Math.exp(score);
            double q = profile.getNormalizedPValue(qraw);
            if (q > 1.0) {
                q = 1.0; //deal with rounding errors that might slightly push q over the limit
            }
            if (q >= threshold) {
                template.location = new RangeLocation(offset, offset + cols - 1);
                SimpleAnnotation ann = new SimpleAnnotation();
                ann.setProperty("score", new Double(q));
                //TODO: we might not need the raw score for anything and getting rid of it may save memory
                ann.setProperty("raw_score", new Double(qraw));
                template.annotation = ann;
                seq.createFeature(template);
            }
        }
        return seq;
    }

    /**
     * Create a new annotator that uses PROBABILITY score type.
     *
     * @param profile   The protein profile used for the annotation (search)
     * @param threshold The normalized p-value threshold
     */
    public NormalizedProfileAnnotator(ProteinProfile profile, double threshold) {
        this.profile = profile;
        this.matrix = profile.getWeightMatrix();
        this.threshold = threshold;
        this.scoreType = ScoreType.PROBABILITY;
    }
}