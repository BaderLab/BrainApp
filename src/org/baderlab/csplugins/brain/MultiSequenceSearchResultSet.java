package org.baderlab.csplugins.brain;

import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;

/**
 * Copyright (c) 2004 Memorial Sloan-Kettering Cancer Center
 * *
 * * Code written by: Gary Bader
 * * Authors: Gary Bader, Ethan Cerami, Chris Sander
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
 * * Date: Feb 16, 2005
 * * Time: 9:29:24 PM
 * * Description
 */

/**
 * A convenience container for multiple sequence search result sets (either profile or regex)
 */
public class MultiSequenceSearchResultSet {
    private HashMap resultSets;

    /**
     * Default constructor
     */
    public MultiSequenceSearchResultSet() {
        resultSets = new HashMap();
    }

    /**
     * Add a result set to the multi set
     */
    public void add(SequenceSearchResultSet results) {
        if (results.getProfile() != null) {
            resultSets.put(results.getProfile(), results);
        } else if (results.getRegex() != null) {
            resultSets.put(results.getRegex(), results);
        } else {
            throw new RuntimeException("No profile or regular expression set in the result set.");
        }
    }

    /**
     * Get all of the result sets in this multi set
     *
     * @return A collection of result sets
     */
    public Collection getAllResultSets() {
        return resultSets.values();
    }

    /**
     * Get a specific result set based on a profile of interest
     *
     * @param profile Get the result set generated with this profile
     */
    public SequenceSearchResultSet getResultSet(ProteinProfile profile) {
        return (SequenceSearchResultSet) resultSets.get(profile);
    }

    /**
     * Get a specific result set based on a regex of interest
     *
     * @param regex Get the result set generated with this regular expression
     */
    public SequenceSearchResultSet getResultSet(String regex) {
        return (SequenceSearchResultSet) resultSets.get(regex);
    }

    /**
     * Cleans up a result set - should call if you are generating a lot of results to prevent memory bloat
     */
    public void clear() {
        Collection results = this.getAllResultSets();
        for (Iterator iterator = results.iterator(); iterator.hasNext();) {
            SequenceSearchResultSet sequenceSearchResultSet = (SequenceSearchResultSet) iterator.next();
            sequenceSearchResultSet.clear();
        }
        resultSets.clear();
    }
}
