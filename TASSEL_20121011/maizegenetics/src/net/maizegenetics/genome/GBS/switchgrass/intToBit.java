/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.genome.GBS.switchgrass;
import java.io.*;
import net.maizegenetics.gbs.util.OpenBitSet;
import net.maizegenetics.genome.GBS.ReadsByTaxaPlasmid;
/**
 *
 * @author fl262
 */
public class intToBit {
	ReadsByTaxaPlasmid rbt;
	public intToBit (String infile, String outfile) {
		rbt = new ReadsByTaxaPlasmid(infile, true);
		output (outfile);
	}
	public void output (String outfile) {
		OpenBitSet obs;
		try {
			DataOutputStream dos  = new DataOutputStream(new BufferedOutputStream (new FileOutputStream(outfile), 65536*16));
			dos.writeInt(rbt.haplotypeNum);
			dos.writeInt(2);
			dos.writeInt(rbt.taxaNum);
			for (int i = 0; i < rbt.taxaNum; i++) {
				dos.writeUTF(rbt.taxaNames[i]);
			}
			for (int i = 0; i < rbt.haplotypeNum; i++) {
				for (int j = 0; j < 2; j++) {
					dos.writeLong(rbt.haplotype[j][i]);
				}
				dos.writeByte(64);
				obs=new OpenBitSet(rbt.taxaNum);
				for (int j = 0; j < rbt.taxaNum; j++) {
                    if(rbt.hapDist[i][j]>0) obs.set(j);
                }
                long[] obsInLong=obs.getBits();
                for (int t = 0; t < obsInLong.length; t++) {
                    dos.writeLong(obsInLong[t]);
                }
			}
			dos.flush();
			dos.close();
		}
		catch (Exception e) {

		}
	}
}
