package org.baderlab.csplugins.brainplugin.analyses;

import org.baderlab.csplugins.brainplugin.*;
import org.baderlab.brain.BrainAlgorithm;

import javax.imageio.ImageIO;
import java.io.File;
import java.io.IOException;
import java.util.List;

/**
 * * Copyright (c) 2004 Memorial Sloan-Kettering Cancer Center
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
 *
 * User: moyez
 * Date: Feb 27, 2007
 * Time: 7:32:50 PM
 * Description: This class contains data analysis code for generating interactions, sequence logos, etc.
 * To change this template use File | Settings | File Templates.
 */
public class WWAnalysis {

    public static void main(String[] args) {
        runMultiSequenceSearch();
        //generateLogos();
    }

    private static void runMultiSequenceSearch() {
        BrainAlgorithm alg = new BrainAlgorithm();
        System.out.println("Testing [mdharsee]...");

        BrainParameterSet params = new BrainParameterSet();
        //params.setProfileFileName(new File("/Users/moyez/Code/wwdata/profiles_2007Feb21/wwprojectfile_od5.txt"));
        params.setProfileFileName(new File("/Users/moyez/Code/wwdata/profiles_2007Feb21/D00003.od5.pep.txt"));
        params.setFuzzFactor(1.0);
        params.setDatabaseFileName(new File("/Users/moyez/Code/seqdata/human.protein.gpff"));
        params.setDatabaseFormat("genpept");
        params.setScoreThreshold(0.0);
        ProteinDatabaseSearchParams dbparams = new ProteinDatabaseSearchParams(true);
        dbparams.setNormalized(true);
        dbparams.setScoreType(ProteinDatabaseSearchParams.SCORE_TYPE_NORM_PROBABILITY);
        params.setSearchParams(dbparams);
        alg.setParams(params);
        MultiSequenceSearchResultSet resultSet = alg.runProfileSearch();

    }

    // generate sequence logos
    private static void generateLogos() {
        String inFile;
        String outPath;
        String fileSuffix;

        System.out.println("brain.WWAnalysis: Generating sequence logos...");

        inFile = "/Users/moyez/Code/wwdata/analysis1/wwprojectfile_od1.txt";
        outPath = "/Users/moyez/Code/wwdata/analysis1/logos/logos_od1";
        fileSuffix = "_od1";
        System.out.println(inFile);
        writeSequenceLogos(inFile, outPath, null, 0.0, fileSuffix);

        inFile = "/Users/moyez/Code/wwdata/analysis1/wwprojectfile_od2.txt";
        outPath = "/Users/moyez/Code/wwdata/analysis1/logos/logos_od2";
        fileSuffix = "_od2";
        System.out.println(inFile);
        writeSequenceLogos(inFile, outPath, null, 0.0, fileSuffix);

        inFile = "/Users/moyez/Code/wwdata/analysis1/wwprojectfile_od3.txt";
        outPath = "/Users/moyez/Code/wwdata/analysis1/logos/logos_od3";
        fileSuffix = "_od3";
        System.out.println(inFile);
        writeSequenceLogos(inFile, outPath, null, 0.0, fileSuffix);

        inFile = "/Users/moyez/Code/wwdata/analysis1/wwprojectfile_od4.txt";
        outPath = "/Users/moyez/Code/wwdata/analysis1/logos/logos_od4";
        fileSuffix = "_od4";
        System.out.println(inFile);
        writeSequenceLogos(inFile, outPath, null, 0.0, fileSuffix);

        inFile = "/Users/moyez/Code/wwdata/analysis1/wwprojectfile_od5.txt";
        outPath = "/Users/moyez/Code/wwdata/analysis1/logos/logos_od5";
        fileSuffix = "_od5";
        System.out.println(inFile);
        writeSequenceLogos(inFile, outPath, null, 0.0, fileSuffix);
    }


    //write sequence logos for given profiles (adapted from BrainAlgorithm.writeSequenceLogos)
    public static void writeSequenceLogos(String profileListFileName,
                                          String outputDirectory,
                                          String codonBiasFileName,
                                          double fuzzFactor,
                                          String fileSuffix) {

        // filename suffix is optional - handle null value
        String suffix = fileSuffix;
        if (suffix == null) {
            suffix = "";
        }

        // codon bias file is optional - handle null value
        File codonBiasFile = null;
        if (codonBiasFileName != null) {
            codonBiasFile = new File(codonBiasFileName);
        }

        // read profiles
        List profileList = PeptideToProfileReader.readPeptidesAsProfiles(new File(profileListFileName), 0, null,
                fuzzFactor, codonBiasFile, true);
        String outFileName = null;

        // generate a sequence logo for each profile
        for (int i = 0; i < profileList.size(); i++) {
            ProteinProfile proteinProfile = (ProteinProfile) profileList.get(i);
            outFileName = new String(outputDirectory + File.separator + proteinProfile.getName() + suffix + ".png");
            ProteinSequenceLogo logo = new ProteinSequenceLogo(proteinProfile, 180);
            try {
                ImageIO.write(logo.drawSequenceLogo(), "png", new File(outFileName));
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }

}
