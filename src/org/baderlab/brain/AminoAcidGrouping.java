package org.baderlab.brain;

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
 * * Date: Feb 20, 2006
 * * Time: 7:23:44 PM
 */

/**
 * Defines a number of amino acid groupings
 */
public class AminoAcidGrouping {

    public static String[] getPositionSpecificPDZGrouping() {
        String[] groupingByPosition = {
                "ST, KRH, DEQN, FLAMPWIVCY, G", //position -3
                "ST, KRH, DEQN, MPVCYFWLIA, G", //position -2 - ST, large hydrophobes, G
                "ST, KRH, DEQN, FLAMPIVCY, G, W", //position -1 - W is separate because of backbone binding capability
                "ST, KRH, DEQN, G, P, VAC, ILM, FWY" //position 0 - hydrophobe by size
        };
        return groupingByPosition;
    }

    public static String getHydroxylBasicAcidicpolarHydrophobeGrouping() {
        return "ST, KRH, DEQN, FLAMPWIVCY, G";
    }

    public static String getPolarChargedHydrophobeGrouping() {
        return "STQN, KRH, DE, LAMIVFWY, C, P, G";
    }

    public static String getGroupingHydrophobeBySize() {
        return "ST, KRH, DEQN, G, P, VAC, ILM, FWY";
    }

}
