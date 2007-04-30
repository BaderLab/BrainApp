package org.baderlab.csplugins.brain.test;

import junit.framework.TestCase;
import org.baderlab.csplugins.brain.*;

import java.io.File;

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
 * * Date: May 12, 2005
 * * Time: 7:54:54 PM
 */

/**
 * Test for the ProteinProfile class
 */
public class ProteinProfileTest extends TestCase {

    ProteinProfile profile1 = null;
    ProteinProfile profile2 = null;
    ProteinDatabaseSearchParams params = null;
    BindingPeptideList peptideList = null;

    /**
     * Set up a few things for this test set
     *
     * @throws Exception
     */
    public void setUp() throws Exception {
        peptideList = new BindingPeptideList();
        peptideList.read("testData" + File.separator + "SH3_pep8.txt");
        profile1 = new ProteinProfile(peptideList.getSequenceIteratorByLength(4, ProteinTerminus.C), 1, new String("testProfile"));
        peptideList = new BindingPeptideList();
        peptideList.read("testData" + File.separator + "PDZ_pep8.txt");
        profile2 = new ProteinProfile(peptideList.getSequenceIteratorByLength(4, ProteinTerminus.C), 1, new String("testProfile"));
    }

    /**
     * Test the profile searching code
     */
    public void testDistanceMetric() {
        //test isolation property
        assertEquals(0.0, ProteinProfileDistance.calculateDistributionDistance(profile1, profile1), 0.000001);

        //test distance > 0 for distant profiles
        double d1 = ProteinProfileDistance.calculateDistributionDistance(profile1, profile2);
        assertTrue(d1 > 0);

        //test symmetry
        double d2 = ProteinProfileDistance.calculateDistributionDistance(profile2, profile1);
        assertTrue(d1 == d2);
    }

    /**
     * Test the NNK reweighting code
     */
    public void testNNKCodonBiasReweight() {
        //just makes sure there are no exceptions and prints to output so user can test full calculation in Excel
        System.out.println(profile1.toString());
        profile1.reWeightByNNKCodonBias();
        System.out.println();
        System.out.println(profile1.toString());
    }

}
