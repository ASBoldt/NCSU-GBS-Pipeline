/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.genome.GBS;

import java.io.File;

/**
 *
 * @author jcg233
 */
public class TestVirtualDigesterJeff {
  public static void main(String[] args) {
       VirtualDigester be=new VirtualDigester(new File("C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/434LFAAXX/test/fake_genome.txt"),
               new File("C:/Documents and Settings/jcg233/My Documents/Bioinformatics/NextGen/434LFAAXX/test/fake_genome.cut.bin"));
  }

}
