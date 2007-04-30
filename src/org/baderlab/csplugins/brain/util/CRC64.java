package org.baderlab.csplugins.brain.util;

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
 * * Date: Mar 28, 2005
 * * Time: 11:02:29 PM
 * * Description CRC64 class stolen from expasy4j (http://www.isb-sib.ch/~ejain/expasy4j/)
 */

/**
 * Computes checksums.
 * <p/>
 * <pre> Checksum checksum = new CRC64();
 * checksum.update("ILIKEMATH");
 * checksum.getValue(); // 6331126313412407383
 * checksum.getString(); // '57DCADD69B133057'</pre>
 */

public class CRC64
        implements java.util.zip.Checksum {
    private static final long POLY64 = 0xD800000000000000L;
    private static final long[] crcTable = new long[256];
    private long crc;


    static {
        for (int i = 0; i < 256; ++i) {
            long part = i;
            for (int j = 0; j < 8; ++j)
                part = ((part & 1) != 0) ?
                        (part >>> 1) ^ POLY64 : (part >>> 1);
            crcTable[i] = part;
        }
    }


    public void update(int b) {
        long low = crc >>> 8;
        long high = crcTable[(int) ((crc ^ b) & 0xFF)];
        crc = low ^ high;
    }


    public void update(byte[] b, int offset, int length) {
        for (int i = offset; i < length; ++i)
            update(b[i]);
    }


    public void update(String s) {
        // update(s.getBytes(), 0, s.length());
        int size = s.length();
        for (int i = 0; i < size; ++i)
            update(s.charAt(i));

    }


    public long getValue() {
        return crc;
    }


    /**
     * Returns a zero-padded 16 character wide string
     * containing the current value of this checksum in
     * uppercase hexadecimal format.
     */

    public String toString() {

        StringBuffer buffer = new StringBuffer();
        buffer.append(Long.toHexString(crc >>> 4));
        buffer.append(Long.toHexString(crc & 0xF));
        for (int i = 16 - buffer.length(); i > 0; --i)
            buffer.insert(0, '0');
        return buffer.toString().toUpperCase();
    }


    public void reset() {
        crc = 0;
    }
}

