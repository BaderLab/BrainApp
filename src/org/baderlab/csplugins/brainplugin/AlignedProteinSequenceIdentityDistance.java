package org.baderlab.csplugins.brainplugin;

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
 * * Date: Aug 17, 2005
 * * Time: 5:25:22 PM
 */

/**
 * Calculated the sequence identity as a distance metric (one minus % sequence identity)
 */
public class AlignedProteinSequenceIdentityDistance extends DistanceMetric {
    public double calc(Object object1, Object object2) {
        if ((!(object1 instanceof String)) || (!(object2 instanceof String))) {
            throw new RuntimeException("Non string passed to AlignedProteinSequenceIdentityDistance");
        }

        double distance = 0.0;
        String sequenceA = (String) object1, sequenceB = (String) object2;
        int count = 0;
        int numGaps = 0;
        for (int i = 0; i < sequenceA.length(); i++) {
            String aaA = sequenceA.substring(i, i + 1);
            String aaB = sequenceB.substring(i, i + 1);
            if (aaA.equals("-") || (aaB.equals("-"))) {
                numGaps++;
                continue;
            }
            if (aaA.equals(aaB)) {
                count++;
            }
        }
        distance = ((double) count / (double) (sequenceA.length() - numGaps));
        distance = 1.0 - distance; //convert from similarity to distance
        return distance;
    }
}
