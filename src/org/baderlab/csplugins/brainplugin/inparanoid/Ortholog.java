package org.baderlab.csplugins.brain.inparanoid;

import org.baderlab.csplugins.brain.DatabaseReference;

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
 * * Time: 7:08:03 PM
 */

/**
 * A single homolog (ortholog + paralog) that is part of an Inparanoid ortholog cluster
 */
public class Ortholog {
    private int clusterID;  //Inparanoid cluster ID
    private int taxid;  //NCBI taxonomy ID
    private DatabaseReference proteinID;

    public Ortholog(int clusterID, int ncbiTaxID, DatabaseReference proteinID) {
        this.clusterID = clusterID;
        this.taxid = ncbiTaxID;
        this.proteinID = proteinID;
    }

    public int getClusterID() {
        return clusterID;
    }

    public int getTaxid() {
        return taxid;
    }

    public DatabaseReference getProteinID() {
        return proteinID;
    }
}
