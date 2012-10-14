/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.genome.solexa;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

/**
 *
 * @author edbuckler
 */
public class SimpleTextFile {
    private ArrayList<Object> elements=new ArrayList<Object>();
    private Object[] colType;
    private int rows;


    public SimpleTextFile(String infile, String[] labels, Object[] colType, boolean header) {
        rows=countRows(infile);
        if(header==true) rows--;
        System.out.println(infile+" has "+rows+" rows");
        this.colType=colType;
        
        for(Object c: colType) {
            if(c==String.class) {
                String[] x=new String[rows];
                elements.add(x);
            } else if(c==Integer.class) {
                int[] x=new int[rows];
                elements.add(x);
            } if(c==Double.class) {
                double[] x=new double[rows];
                elements.add(x);
            } 
        }
        try{
            BufferedReader fileIn = new BufferedReader(new FileReader(infile), 100000);
            if(header) fileIn.readLine();
            for (int i = 0; i < rows; i++) {
                String[] s=fileIn.readLine().split("\\s");
                for(int j=0; j<colType.length; j++) {
                    if(colType[j]==String.class) {
                        setElement(j,i,s[j]);
                    } else if(colType[j]==Integer.class) {
                        if(s[j].equals("NaN")) {setElement(j,i,Integer.MIN_VALUE);}
                        else {setElement(j,i,Integer.parseInt(s[j]));}
                    } if(colType[j]==Double.class) {
                        setElement(j,i,Double.parseDouble(s[j]));
                    }
                }
            }
            fileIn.close();
        } catch (IOException e) {
            System.out.println("countRows error at row");
            e.printStackTrace();
        }
    }

    private void setElement(int col, int row, int value) {
        ((int[])elements.get(col))[row]=value;
    }

    private void setElement(int col, int row, double value) {
        ((double[])elements.get(col))[row]=value;
    }

    private void setElement(int col, int row, String value) {
        ((String[])elements.get(col))[row]=value;
    }

    public String getStringElement(int col, int row) {
         return ((String[])elements.get(col))[row];
    }

    public double getDoubleElement(int col, int row) {
         return ((double[])elements.get(col))[row];
    }

    public int getIntElement(int col, int row) {
         return ((int[])elements.get(col))[row];
    }

    public int getRows() {
        return rows;
    }

    public static int countRows(String infileName) {
        int row=0;
        try{
            BufferedReader fileIn = new BufferedReader(new FileReader(infileName), 100000);
            while(fileIn.readLine()!=null) {
                row++;
            }
            fileIn.close();
        } catch (IOException e) {
            System.out.println("countRows error at row"+row);
            e.printStackTrace();
        }
        return row;
    }

}
