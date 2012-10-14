/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.gbs.pav;

import cern.colt.Swapper;
import cern.colt.function.IntComparator;
import cern.colt.GenericSorting;
import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.TreeSet;

/**
 *
 * @author Fei Lu
 */
public class Pedigree {
	String[] sampleFamily;
	String[] sampleName;
	String[] p1;
	String[] p2;
	float[] con1;
	float[] con2;
	String[] familyName;
	int[] familySize;
	int[] familyStartIndex;

	public Pedigree (String infileS) {
		this.readPedigree(infileS);
		this.sortByFamilyAndName();
		this.generateFamilyInfo();
	}

	public int getAllSampleSize () {
		return sampleName.length;
	}

	public int getFamilyNum () {
		return familyName.length;
	}

	public int getFamilySize (String family) {
		int hit = Arrays.binarySearch(familyName, family);
		return familySize[hit];
	}

	public String[] getFamilyStartWith (String stS, int minFamilySize) {
		ArrayList<String> nameList = new ArrayList();
		for (int i = 0; i < familyName.length; i++) {
			if (familyName[i].startsWith(stS) && familySize[i] >= minFamilySize) {
				nameList.add(familyName[i]);
			}
		}
		String[] nameArray = nameList.toArray(new String[nameList.size()]);
		return nameArray;
	}

	public String[][] getSampleOfFamilies (String[] families, int minFamilySize) {
		String[][] names = new String[families.length][];
		ArrayList<String[]> namesList = new ArrayList();
		for (int i = 0; i < families.length; i++) {
			String[] temp = this.getSampleOfFamily(families[i], minFamilySize);
			if (temp == null) continue;
			namesList.add(temp);
		}
		names = namesList.toArray(new String[namesList.size()][]);
		return names;
	}

	public String[] getSampleOfFamily (String family, int minFamilySize) {
		int hit = Arrays.binarySearch(familyName, family);
		if (familySize[hit] < minFamilySize) return null;
		int start = familyStartIndex[hit];
		String[] samples = new String[familySize[hit]];
		for (int i = start; i < start + samples.length; i++) {
			samples[i-start] = sampleName[i];
		}
		return samples;
	}

	public String[] getSampleContains (String stS) {
		ArrayList<Integer> fList = new ArrayList();
		for (int i = 0; i < familyName.length; i++) {
			if (familyName[i].contains(stS)) fList.add(i);
		}
		Integer[] fArray = fList.toArray(new Integer[fList.size()]);
		int size = 0;
		for (int i = 0; i < fArray.length; i++) {
			size += familySize[fArray[i]];
		}
		String[] samples = new String[size];
		int cnt = 0;
		for (int i = 0; i < fArray.length; i++) {
			for (int j = familyStartIndex[fArray[i]]; j < familySize[fArray[i]]; j++) {
				samples[cnt] = sampleFamily[j];
				cnt++;
			}
		}
		return samples;
	}

	public String[] getSampleStartWith (String stS) {
		ArrayList<Integer> fList = new ArrayList();
		for (int i = 0; i < familyName.length; i++) {
			if (familyName[i].startsWith(stS)) fList.add(i);
		}
		Integer[] fArray = fList.toArray(new Integer[fList.size()]);
		int size = 0;
		for (int i = 0; i < fArray.length; i++) {
			size += familySize[fArray[i]];
		}
		String[] samples = new String[size];
		int cnt = 0;
		for (int i = 0; i < fArray.length; i++) {
			for (int j = familyStartIndex[fArray[i]]; j < familySize[fArray[i]]; j++) {
				samples[cnt] = sampleFamily[j];
				cnt++;
			}
		}
		return samples;
	}

	private void generateFamilyInfo () {
		TreeSet<String> familySet = new TreeSet();
		for (int i = 0; i < sampleFamily.length; i++) {
			familySet.add(sampleFamily[i]);
		}
		familyName = familySet.toArray(new String[familySet.size()]);
		familySize = new int[familyName.length];
		familyStartIndex = new int[familyName.length];
		for (int i = 0; i < familyName.length; i++) {
			int hit = Arrays.binarySearch(sampleFamily, familyName[i]);
			while (hit > -1 && sampleFamily[hit].equals(familyName[i])) hit--;
			hit++;
			familyStartIndex[i] = hit;
			familySize[i] = 0;
			for (int j = familyStartIndex[i]; j < sampleFamily.length; j++) {
				if (familyName[i].equals(sampleFamily[j])) {
					familySize[i]++;
				}
				else {
					break;
				}
			}
		}
	}

	private void readPedigree (String infileS) {
		try {
			BufferedReader br = new BufferedReader (new FileReader(infileS), 65536);
			br.readLine();
			ArrayList<String> lineList = new ArrayList();
			String temp;
			while ((temp = br.readLine()) != null) {
				lineList.add(temp);
			}
			String[] line = lineList.toArray(new String[lineList.size()]);
			this.initialize(line.length);
			for (int i = 0; i < line.length; i++) {
				String[] tem = line[i].split("\\s+");
				sampleFamily[i] = tem[0];
				tem[1] = tem[1].replaceFirst(":\\d\\d\\d\\d\\d\\d\\d+", ""); //because of name convertion, this will be changed after the new build is done
				sampleName[i] = tem[1];
				p1[i] = tem[2];
				p2[i] = tem[3];
				con1[i] = Float.valueOf(tem[4]);
				con2[i] = Float.valueOf(tem[5]);
			}
			br.close();
		}
		catch (Exception e) {
			System.out.println(e.toString());
			System.exit(1);
		}
	}

	private void initialize (int size) {
		sampleFamily = new String[size];
		sampleName = new String[size];
		p1 = new String[size];
		p2 = new String[size];
		con1 = new float[size];
		con2 = new float[size];
	}

	public void sortByFamilyAndName () {
		GenericSorting.quickSort(0, sampleFamily.length, compFamilyAndName, swapper);
	}

	public void sortByName() {
		GenericSorting.quickSort(0, sampleFamily.length, compName, swapper);
	}

	IntComparator compName = new IntComparator() {
		public int compare(int a, int b) {
			return sampleName[a].compareTo(sampleName[b]);
		}
	};

	IntComparator compFamilyAndName = new IntComparator() {
		public int compare(int a, int b) {
			if (sampleFamily[a].equals(sampleFamily[b])) {
				return sampleName[a].compareTo(sampleName[b]);
			}
			else {
				return sampleFamily[a].compareTo(sampleFamily[b]);
			}
		}
	};

	Swapper swapper = new Swapper() {
		public void swap(int a, int b) {
			String tempS;
			tempS = sampleFamily[a]; sampleFamily[a] = sampleFamily[b]; sampleFamily[b] = tempS;
			tempS = sampleName[a]; sampleName[a] = sampleName[b]; sampleName[b] = tempS;
			tempS = p1[a]; p1[a] = p1[b]; p1[b] = tempS;
			tempS = p2[a]; p2[a] = p2[b]; p2[b] = tempS;
			float tempf;
			tempf = con1[a]; con1[a] = con1[b]; con1[b] = tempf;
			tempf = con2[a]; con2[a] = con2[b]; con2[b] = tempf;
		}
	};


}
