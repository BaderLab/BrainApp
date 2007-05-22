package org.baderlab.brain.inparanoid.test;

import junit.framework.TestCase;
import org.biojava.bio.BioException;
import org.baderlab.brain.inparanoid.InparanoidDB;

import java.io.IOException;

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
 * * Date: Jun 8, 2005
 * * Time: 9:00:53 PM
 */

/**
 * Tests for the InparanoidDB class
 */
public class HomologeneDBTest extends TestCase {
    InparanoidDB hdb;

    /**
     * Set up a few things for this test set
     *
     * @throws Exception
     */
    public void setUp() throws Exception {
        hdb = new InparanoidDB();
        //TODO: create test
    }

    /**
     * Test the profile searching code
     */
    public void testProfileSingleHit() throws BioException, IOException {
        //assertEquals(true, hdb.isHomologByAccession("NP_031408", "NP_058682"));
        //assertEquals(false, hdb.isHomologByAccession("NP_000009", "NP_058682"));
    }

}
