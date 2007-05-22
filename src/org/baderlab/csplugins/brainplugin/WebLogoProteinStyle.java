package org.baderlab.csplugins.brainplugin;

import org.biojava.bio.gui.SymbolStyle;
import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;

import java.awt.*;

/**
 * Paints colors for the sequence logo according to Berkeley WebLogo
 */
public class WebLogoProteinStyle implements SymbolStyle {

    public WebLogoProteinStyle() {
    }

    public Paint outlinePaint(Symbol s) throws IllegalSymbolException {
        return fillPaint(s);
    }

    public Paint fillPaint(Symbol s) throws IllegalSymbolException {
        Paint fillColor = null;
        //polar amino acids (G,S,T,Y,C,Q,N) are green
        if (s.equals(ProteinTools.g()) || s.equals(ProteinTools.s()) || s.equals(ProteinTools.t()) ||
                s.equals(ProteinTools.y()) || s.equals(ProteinTools.c()) || s.equals(ProteinTools.q()) ||
                s.equals(ProteinTools.n())) {
            float[] hsb = Color.RGBtoHSB(4, 206, 4, null);
            fillColor = Color.getHSBColor(hsb[0], hsb[1], hsb[2]);
        } else if (s.equals(ProteinTools.d()) || s.equals(ProteinTools.e())) {
            //acidic (D,E) red
            float[] hsb = Color.RGBtoHSB(204, 2, 4, null);
            fillColor = Color.getHSBColor(hsb[0], hsb[1], hsb[2]);
        } else if (s.equals(ProteinTools.k()) || s.equals(ProteinTools.r()) || s.equals(ProteinTools.h())) {
            //basic (K,R,H) blue
            float[] hsb = Color.RGBtoHSB(4, 2, 204, null);
            fillColor = Color.getHSBColor(hsb[0], hsb[1], hsb[2]);
        } else {
            //hydrophobic (A,V,L,I,P,W,F,M) amino acids are black
            fillColor = Color.black;
        }
        return fillColor;
    }
}
