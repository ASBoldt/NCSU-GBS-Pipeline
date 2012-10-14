/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package Test;

import java.io.*;
/**
 *
 * @author fl262
 */
public class readwrite {
    public static void main (String[] args)  {
        FileInputStream fis;
        FileOutputStream fos;
        BufferedInputStream bis;
        BufferedOutputStream bos;
        try {
            fis = new FileInputStream ("E:/a/b/b.txt");
            fos = new FileOutputStream ("E:/a/b/d.txt");
            bis = new BufferedInputStream (fis, 65536);
            bos = new BufferedOutputStream (fos, 65536);

            FileReader fr = new FileReader ("E:/a/b/a.txt");
            FileWriter fw = new FileWriter ("E:/a/b/c.txt");
            BufferedReader br = new BufferedReader (fr, 65536);
            BufferedWriter bw = new BufferedWriter (fw, 65536);

            DataInputStream dis = new DataInputStream (bis);
            DataOutputStream dos = new DataOutputStream (bos);

            int c;
            long ab = 234;
            while ((c = bis.read()) != -1) {
                System.out.println (c);

                dos.write(c);



            }
            dis.close();
            dos.close();
        }
        catch (Exception e) {

        }
    }
}
