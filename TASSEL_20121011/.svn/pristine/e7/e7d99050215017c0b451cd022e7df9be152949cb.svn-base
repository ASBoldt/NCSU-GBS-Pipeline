/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.genome.GBS.plasmid;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

/**
 *
 * @author fl262
 */
public class parseBarcode {
	public static final String[] overhang={"CAGC","CTGC"};
	public static final String[] likelyReadEnd={"GCTGGATC","GCAGGATC","GCTGAGAT", "GCAGAGAT","GCAGC","GCTGC"};
	static String nullS="AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
	barcode[] barcodeArray;
	int[] barcodeInt;
	private HashMap<Integer,Integer> quickMap;
	static int maxBarcodeLength=8;

	public parseBarcode (String qseqDirS, String keyFileS, String parsedTaxaDirS, int minQual) {
		processQseq (new File(qseqDirS), new File(keyFileS), new File(parsedTaxaDirS), minQual);
	}
	public void setupBarcodeFiles (File keyFile, File parsedTaxaDir, String flowcell, String lane) {
		try{
            BufferedReader br=new BufferedReader(new FileReader(keyFile),65536);
            ArrayList<barcode> barcodeList=new ArrayList();
            String line;
            while (((line = br.readLine()) != null)) {
                String[] temp = line.split("\\s");  //split by whitespace
                if(temp[0].equals(flowcell) && temp[1].equals(lane)) {
                    File nf=new File(parsedTaxaDir.getPath()+"/"+temp[3]+"_"+flowcell+"_"+lane+"_"+temp[2]+".cnt");
                    DataOutputStream dos=new DataOutputStream(new BufferedOutputStream(new FileOutputStream(nf),65536));
                    barcode theBC = new barcode(temp[2], overhang);
                    theBC.setTheDOS(dos);
                    barcodeList.add(theBC);
                    System.out.println(theBC.barcodeS+" "+temp[3]);
                }
            }
			barcodeArray = barcodeList.toArray(new barcode[barcodeList.size()]);
			barcodeInt = new int[barcodeArray.length*2];
			quickMap = new HashMap();
			for (int i = 0; i< barcodeArray.length; i++) {
				barcodeInt[i*2]=barcodeArray[i].barOverInt[0];
				barcodeInt[i*2+1]=barcodeArray[i].barOverInt[1];
				quickMap.put(barcodeArray[i].barOverInt[0], i);
				quickMap.put(barcodeArray[i].barOverInt[1], i);
			}
			Arrays.sort(barcodeInt);
        } catch(Exception e){
            System.out.println("Error with setupBarcodeFiles: "+e);
        }

    }
	public void parseQseq (File qseqFile, int minQual) {
		int goodBarcodedReads=0, highQual=0, totalReads=0;
		int usedLen = baseEncoder.seqLen+maxBarcodeLength;
        int max=(int)100e7;
        try{
			System.out.println("File = " + qseqFile);
			BufferedReader br = new BufferedReader (new FileReader(qseqFile), 65536);
            String line, sl, qualS;
            while (((line = br.readLine()) != null)&&(goodBarcodedReads<max)) {
                String[] temp = line.split("\t");
                sl = temp[8];
                qualS = temp[9];
                totalReads++;
                int firstBadBase=baseEncoder.getFirstLowQualityPos(qualS, minQual);
                if(firstBadBase > usedLen) {highQual++;} else {continue;}
                int miss = -1;
                miss=sl.indexOf('.');
                if((miss!=-1)&&(miss < usedLen)) continue;  //bad sequence so skip
                barcode bestBarcode=findBestBarcode(sl);
                if(bestBarcode==null) continue;  //overhang missing so skip
                String genomicSeq=sl.substring(bestBarcode.barLength, sl.length());
                String hap=removeSeqAfterSecondCutSite(genomicSeq);  //this is slow 20% of total time

				byte[] byteArray = baseEncoder.getByteArrayFromSeq(hap);
				bestBarcode.getTheDOS().write(byteArray);
                bestBarcode.getTheDOS().writeInt(1);
                goodBarcodedReads++;
                if(goodBarcodedReads%100000==0)  System.out.println("Total Reads="+totalReads+
                         "  HighQuality="+highQual+ " barcoded= " + goodBarcodedReads);
            }
			br.close();
        }
        catch(Exception e) {
			System.out.println("Total reads = "+totalReads);
            System.out.println("Catch c="+goodBarcodedReads+" e="+e);
        }
		System.out.println("Count of Tags="+goodBarcodedReads);
		try {
			for (int i = 0; i < barcodeArray.length; i++) {
				barcodeArray[i].getTheDOS().flush();
				barcodeArray[i].getTheDOS().close();
			}
		}
		catch (Exception e) {
			System.out.println("Error happened while closing theDOS "+ e);
		}
	}

	private barcode findBestBarcode(String queryS) {
        int query=baseEncoder.getIntFromSeq(queryS.substring(0, baseEncoder.intSize));
        int closestHit=Arrays.binarySearch(barcodeInt, query);
		if (closestHit < -1) {
			int index = quickMap.get(barcodeInt[-(closestHit+2)]);
			if (barcodeArray[index].compareSequence(query)) {
				return barcodeArray[index];
			}
			else {
				return null;
			}
		}
		else if (closestHit == -1) {
			return null;
		}
		else if (closestHit < barcodeInt.length) {
			int index = quickMap.get(barcodeInt[closestHit]);
			return barcodeArray[index];
		}
		else {
			return null;
		}
    }
	public void processQseq(File qseqDir, File keyFile, File parsedTaxaDir, int minQual) {
		File[] qseq = qseqDir.listFiles();
		for (int i = 0; i < qseq.length; i++) {
			String[] temp = qseq[i].getName().split("_");
			if (temp.length == 3) {
				setupBarcodeFiles(keyFile, parsedTaxaDir, temp[0], temp[1]);
			}
			else if (temp.length == 4) {
				setupBarcodeFiles(keyFile, parsedTaxaDir, temp[0], temp[2]);
			}
			else if (temp.length == 5) {
				setupBarcodeFiles(keyFile, parsedTaxaDir, temp[1], temp[3]);
			} else {
				System.out.println("Error in parsing file name (e.g. flowcell_lane_qseq.txt):" + qseq[i].toString());
				continue; // skip over files that don't match the naming convention
			}
			parseQseq (qseq[i], minQual);
		}
	}
	public String removeSeqAfterSecondCutSite(String seq) {
         //this looks for a second restriction site, and then turns the remaining sequence to AAAA
        int pos=9999;
        int startSearchBase=3;
        for(String fcs: likelyReadEnd) {
            int p=seq.indexOf(fcs, startSearchBase);
            if((p>startSearchBase)&&(p<pos)) {
                pos=p;
            }
        }
        if(pos<baseEncoder.seqLen) {
            int slen=seq.length();
              //  System.out.print("seq:"+seq);
            seq=seq.substring(0,pos+4)+nullS.substring(0, slen-pos-4);//this makes all bases after the second site AAAAAAA
                //we cut after base 4 as it could ligated to an adapter or another sequence
                //this is the length of the overhang+1

              //   System.out.println(" rseq:"+seq);
        }
        return seq;
    }
	public static void main (String[] args) {
		String qseqDirS = "M:/Plastid_GBS/qseq/";
		String keyFileS = "M:/Plastid_GBS/key/61VBRAAXX_key.txt";
		String parsedTaxaDirS = "M:/Plastid_GBS/parsedTaxa/";
		new parseBarcode (qseqDirS, keyFileS, parsedTaxaDirS, 10);
	}
}
class barcode {
	String barcodeS;
    int[] barOverInt;
    int barOverLength, barLength;
    DataOutputStream theDOS=null;

    public barcode(String barcodeS, String[] overhangS) {
        this.barcodeS=barcodeS;
        Arrays.sort(overhangS);
        barOverInt=new int[overhangS.length];
        for(int i=0; i<overhangS.length; i++) {
            barOverInt[i]=baseEncoder.getIntFromSeq(barcodeS+overhangS[i]);
        }
        barOverLength=barcodeS.length()+overhangS[0].length();
        barLength=barcodeS.length();
    }

    public DataOutputStream getTheDOS() {
        return theDOS;
    }

    public void setTheDOS(DataOutputStream theDOS) {
        this.theDOS = theDOS;
    }
    public int compareSequence(int queryInt, int maxDivCheck) {
        int div=barOverLength;
        for(int targetInt: barOverInt) {
            int c=baseEncoder.seqDifferencesForSubset(targetInt, queryInt, barOverLength, maxDivCheck);
            if(c<div) div=c;
        }
        return div;
    }
	public boolean compareSequence (int queryInt) {
		return baseEncoder.ifSeqSame(queryInt, barOverInt[0], barOverLength)|baseEncoder.ifSeqSame(queryInt, barOverInt[1], barOverLength);
	}
}