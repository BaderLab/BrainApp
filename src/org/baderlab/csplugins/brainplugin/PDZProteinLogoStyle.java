package org.baderlab.csplugins.brainplugin;

import org.biojava.bio.gui.SymbolStyle;
import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;

import java.awt.*;

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
 * * Date: Feb 23, 2006
 * * Time: 10:54:02 PM
 */

/**
 * Paints colors for the sequence logo according to PDZ 2007
 * This coloring was developed by Sachdev Sidhu
 */
public class PDZProteinLogoStyle implements SymbolStyle {

    public PDZProteinLogoStyle() {
    }

    public Paint outlinePaint(Symbol s) throws IllegalSymbolException {
        return fillPaint(s);
    }

    public Paint fillPaint(Symbol s) throws IllegalSymbolException {
        Paint fillColor = null;
        //(STQN) are green
        if (s.equals(ProteinTools.s()) || s.equals(ProteinTools.t()) || s.equals(ProteinTools.q()) ||
                s.equals(ProteinTools.n())) {
            float[] hsb = Color.RGBtoHSB(4, 206, 4, null);
            fillColor = Color.getHSBColor(hsb[0], hsb[1], hsb[2]);
        } else if (s.equals(ProteinTools.d()) || s.equals(ProteinTools.e())) {
            //acidic/amide (DE) red
            float[] hsb = Color.RGBtoHSB(204, 2, 4, null);
            fillColor = Color.getHSBColor(hsb[0], hsb[1], hsb[2]);
        } else if (s.equals(ProteinTools.k()) || s.equals(ProteinTools.r()) || s.equals(ProteinTools.h())) {
            //basic (KRH) blue
            float[] hsb = Color.RGBtoHSB(4, 2, 204, null);
            fillColor = Color.getHSBColor(hsb[0], hsb[1], hsb[2]);
        } else if (s.equals(ProteinTools.g())) {
            float[] hsb = Color.RGBtoHSB(102, 0, 102, null);  //purple
            fillColor = Color.getHSBColor(hsb[0], hsb[1], hsb[2]);
        } else if (s.equals(ProteinTools.c())) {
            float[] hsb = Color.RGBtoHSB(255, 170, 0, null);   //orange
            fillColor = Color.getHSBColor(hsb[0], hsb[1], hsb[2]);
        } else if (s.equals(ProteinTools.p())) {
            float[] hsb = Color.RGBtoHSB(102, 51, 0, null);  //brown
            fillColor = Color.getHSBColor(hsb[0], hsb[1], hsb[2]);
        } else {
            //hydrophobic (A,V,L,I,M,FYW) amino acids are black
            fillColor = Color.black;
        }
        return fillColor;
    }
}
