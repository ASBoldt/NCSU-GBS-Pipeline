/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.genome.GBS;

import net.maizegenetics.genome.solexa.LDbyBitsQuiet;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.ImportUtils;

/**
 *
 * @author glaubitz
 */
public class LoadHapMapRunLDbyBitsJeff {

    public static void main(String[] args) {
        String inFile, outFile;
        Alignment genosOnOneChr;
        if (args.length != 2) {
            System.out.println("ERROR: Two arguments expected: inFile outFile");
        }
        else {
            inFile = args[0];
            outFile = args[1];
            System.out.println("Loading hapmap file...");
            genosOnOneChr = ImportUtils.createPack1AlignmentFromFile(inFile, null, null);
            String poly = genosOnOneChr.isAllPolymorphic() ? "true" : "false";
            System.out.println("\tLocus: " + genosOnOneChr.getLocus(0) + " Sites: " + genosOnOneChr.getSiteCount()
                    + " Taxa: " + genosOnOneChr.getSequenceCount() + " allPoly: " + poly);
            LDbyBitsQuiet lbb=new LDbyBitsQuiet(genosOnOneChr);
            lbb.writeResultsToFile(outFile);
            System.out.println("...Finished");
        }
    }
}
