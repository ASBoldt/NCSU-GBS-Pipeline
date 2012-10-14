/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.gbs.pipeline;
import java.util.HashMap;
/**
 *
 * @author qs24
 */
public class genotype_score 
{
    static double error = 0.001;
    static double v1 = Math.log10(1.0 - error * 3.0 /4.0);
    static double v2 = Math.log10(error/4);
    static double v3 = Math.log10(0.5 - (error/4.0));
    public final static HashMap <String, int[]> AllFrq2GenotypeScore = new HashMap();
    public final static HashMap <String, Byte> genotype2code = new HashMap<String, Byte>();
    public final static HashMap <Byte, String> code2genotype = new HashMap<Byte, String>();
    static
    {
        for (int t=0; t<=127; t++)
        {
            for (int k=0; k<=127; k++)
            {
                int[] scores = calcscore(t, k);
                AllFrq2GenotypeScore.put(Integer.toString(t) + "," + Integer.toString(k), scores);
            }
        }
        String alleles = "ACGT+-";
        Byte code =0;
        genotype2code.put("NN", code);
        code2genotype.put(code, "NN");
        for (int t=0; t<alleles.length(); t++)
        {
            for (int k=0; k<alleles.length(); k++)
            {
                char[] charArray = new char[] {alleles.charAt(t), alleles.charAt(k)};
                String genotype = new String(charArray);
                code ++;
                genotype2code.put(genotype, code);
                code2genotype.put(code, genotype);
                
            }
        }
        
    }

    public static byte getcodefromgenotype (char a, char b)
    {
        char[] charArray = new char[] {a, b};
        String genotype = new String(charArray);
        //System.out.println (genotype);
        return genotype2code.get(genotype);
    }
    public static String getgenotypefromcode (byte a)
    {
        return code2genotype.get(a);
    }
        
    public static int[] getscore(int a, int b)
    {
        if (a>127) a=127;
        if (b>127) b=127;
        return AllFrq2GenotypeScore.get(Integer.toString(a) + "," + Integer.toString(b));
    }
    private static int[] calcscore (int a, int b)
    {   
        int[] results= new int[4];
        int n = a+b;
        int m =a;
        if (b>m) m=b;

        double fact =0;
        if (n>m)
        {
            for (int i=n; i>m; i--)
            {
               fact += Math.log10(i);
            }
            for (int i=1; i<=(n-m); i++)
            {
               fact -= Math.log10(i);
            }
        }
        double aad = Math.pow(10, fact + (double)a * v1 + (double) b * v2);
        double abd = Math.pow(10, fact + (double)n * v3);
        double bbd = Math.pow(10, fact + (double)b * v1 + (double) a * v2);
        double md = aad;
        if (md<abd) {md=abd;}
        if (md<bbd) {md=bbd;}
        int gq = 0;
        if ((aad+abd+bbd)>0)
        { gq = (int) (md/(aad+abd+bbd) * 100);}
        
        
        
        
        int aa =(int) (-10 * (fact + (double)a * v1 + (double) b * v2));
        int ab =(int) (-10 * (fact + (double)n * v3));
        int bb =(int) ( -10 * (fact + (double)b * v1 + (double) a * v2));
        
        m = aa;
        if (m>ab) {m=ab;}
        if (m>bb) {m=bb;}
        aa-=m; ab-=m; bb-=m;
        results[0]= aa>255 ? 255 : aa;
        results[1]= ab>255 ? 255 : ab;
        results[2]= bb>255 ? 255 : bb;
        results[3]= gq;
        //System.out.println ("a" +a + "\tb" + b + "\t" + aa +  "\t" + ab + "\t" +bb + "\t" + gq);
        
        
        return results;
        
    }
    
}
