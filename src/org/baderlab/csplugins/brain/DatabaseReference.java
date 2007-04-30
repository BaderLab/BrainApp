package org.baderlab.csplugins.brain;

/**
 * Copyright (c) 2004 Memorial Sloan-Kettering Cancer Center
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
 * * Date: Mar 28, 2005
 * * Time: 6:06:31 PM
 * * Description A simple database reference
 */

/**
 * A simple database reference
 */
public class DatabaseReference {
    private String dbname = null;   //the name of the database (required)
    private String dbid = null;     //the ID of the record in the database (required)

    /**
     * Creates a DatabaseReference
     *
     * @param dbname The name of the database referenced (required)
     * @param dbid   The identifier of the referenced database record (required)
     */
    public DatabaseReference(String dbname, String dbid) {
        if (dbname == null || dbid == null) {
            throw new IllegalArgumentException("DatabaseReference must contain a database name and ID (dbname=" + dbname + "; dbid=" + dbid + ").");
        }
        this.dbname = dbname.toUpperCase(); //convert to uppercase to normalize - useful for hash maps
        this.dbid = dbid.toUpperCase();
    }

    /**
     * Creates a DatabaseReference from a colon separated dbName:ID string
     * e.g. GenBank:AA12345.1
     *
     * @param dbNameID
     */
    public DatabaseReference(String dbNameID) {
        if (dbNameID == null) {
            throw new IllegalArgumentException("DatabaseReference must contain a database name and ID.");
        }
        String dbName = null;
        String dbID = null;
        String dbXrefSplit[] = dbNameID.split(":");
        if (dbXrefSplit.length == 2) {
            dbName = dbXrefSplit[0];
            dbID = dbXrefSplit[1];
        } else {
            throw new IllegalArgumentException("dbNameID parameter was malformed.  Expected a colon separated string, but found " + dbNameID + ".");
        }
        this.dbname = dbName.toUpperCase();
        this.dbid = dbID.toUpperCase();
    }

    public String getDbname() {
        return dbname;
    }

    public String getDbid() {
        return dbid;
    }

    public void setDbname(String dbname) {
        this.dbname = dbname.toUpperCase();
    }

    public void setDbid(String dbid) {
        this.dbid = dbid.toUpperCase();
    }

    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof DatabaseReference)) return false;

        final DatabaseReference databaseReference = (DatabaseReference) o;

        if (dbid != null ? !dbid.equals(databaseReference.dbid) : databaseReference.dbid != null) return false;
        if (dbname != null ? !dbname.equals(databaseReference.dbname) : databaseReference.dbname != null) return false;

        return true;
    }

    public int hashCode() {
        int result;
        result = (dbname != null ? dbname.hashCode() : 0);
        result = 29 * result + (dbid != null ? dbid.hashCode() : 0);
        return result;
    }

    /**
     * Returns a basic string representation of this DatabaseReference
     */
    public String toString() {
        StringBuffer sb = new StringBuffer();
        sb.append(this.getDbname() + ":" + this.getDbid());
        return sb.toString();
    }
}
