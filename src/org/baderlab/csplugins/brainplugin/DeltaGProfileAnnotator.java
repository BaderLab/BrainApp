package org.baderlab.csplugins.brainplugin;

import org.biojava.bio.BioException;
import org.biojava.bio.SimpleAnnotation;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dp.WeightMatrix;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceAnnotator;
import org.biojava.bio.seq.impl.ViewSequence;
import org.biojava.bio.symbol.*;
import org.biojava.utils.ChangeVetoException;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;

/**
 * Created by IntelliJ IDEA.
 * User: bader
 * Date: Jul 20, 2006
 * Time: 6:35:30 PM
 * To change this template use File | Settings | File Templates.
 */

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

public class DeltaGProfileAnnotator implements SequenceAnnotator,
        Serializable {
    private WeightMatrix matrix;
    private double threshold;
    private ArrayList alphabet;
    private double[] optimalWeightsInWeightMatrix;

    public Sequence annotate(Sequence seq) throws IllegalAlphabetException,
            BioException, ChangeVetoException {
        seq = new ViewSequence(seq);

        int cols = matrix.columns();
        Feature.Template template = new Feature.Template();
        template.source = "DeltaGProfileAnnotator";
        for (int offset = 1;
             offset <= seq.length() - cols + 1;
             offset++) {
            double score = scoreWeightMatrix(matrix, seq, offset);
            if (score <= threshold) {
                template.location = new RangeLocation(offset, offset + cols - 1);
                SimpleAnnotation ann = new SimpleAnnotation();
                ann.setProperty("score", new Double(score));
                template.annotation = ann;
                seq.createFeature(template);
            }
        }
        return seq;
    }

    /**
     * Create a new annotator that uses DeltaG score type.
     *
     * @param profile   The protein profile used for the annotation (search)
     * @param threshold The delta G threshold
     */
    public DeltaGProfileAnnotator(ProteinProfile profile, double threshold) {
        this.matrix = profile.getWeightMatrix();
        this.threshold = threshold;
        this.alphabet = get20aaAlphabet();
        try {
            this.optimalWeightsInWeightMatrix = setBestWeights(this.matrix);
        } catch (IllegalSymbolException e) {
            e.printStackTrace(); //this should never happen because we're using a standard alphabet.
        }
    }

    /**
     * Finds the symbols with the highest weights in each column of the weight matrix and stores them in
     * an array.
     *
     * @param wm The weight matrix to analyze
     * @return An array of highest weighted symbols in each position (each position in the array corresponds
     *         to a position in the weight matrix
     */
    private double[] setBestWeights(WeightMatrix wm) throws IllegalSymbolException {
        double[] bestWeights = new double[wm.columns()];
        for (int i = 0; i < wm.columns(); i++) {
            Distribution d = wm.getColumn(i);
            Iterator l = alphabet.iterator();
            double bestWeight = 0.0;
            //find the symbol with the best weight
            while (l.hasNext()) {
                Symbol symbol = (Symbol) l.next();
                if (bestWeight < d.getWeight(symbol)) {
                    bestWeight = d.getWeight(symbol);
                }
            }
            //store the symbol with the best weight
            bestWeights[i] = bestWeight;
        }
        return (bestWeights);
    }

    /**
     * Saves a 20aa biojava alphabet for ease of use (ignore the U symbol)
     */
    private ArrayList get20aaAlphabet() {
        ArrayList alphabet = new ArrayList(20);
        HashMap alphabetMap = ProteinSequenceUtil.get20aaAlphabet();
        Collection symbols = alphabetMap.values();
        for (Iterator iterator = symbols.iterator(); iterator.hasNext();) {
            Symbol symbol = (Symbol) iterator.next();
            alphabet.add(symbol);
        }
        return alphabet;
    }

    /**
     * Scores the SymbolList from symbol start to symbol (start+columns) with a
     * weight matrix using a score similar to the delta G calculation = -RTln(K).
     * <p/>
     * <p/>
     * This method allows you to use score types such as ScoreType.ODDS. The other
     * scoreWeightMatrix methods gives a result similar or identical to
     * ScoreType.PROBABILITY.
     * </p>
     *
     * @param matrix  the weight matrix used to evaluate the sequences
     * @param symList the SymbolList to assess
     * @param start   the index of the first symbol in the window to evaluate
     * @return the sum of log scores of this weight matrix
     *         having generated symbols start to (start + columns) of symList
     */
    public double scoreWeightMatrix(
            WeightMatrix matrix,
            SymbolList symList,
            int start)
            throws IllegalSymbolException {
        double score = 0;
        int cols = matrix.columns();
        Distribution d;

        for (int c = 0; c < cols; c++) {
            d = matrix.getColumn(c);
            score += Math.log(d.getWeight(symList.symbolAt(c + start)) / optimalWeightsInWeightMatrix[c]);
        }

        return -score; //negate the score
    }
}
