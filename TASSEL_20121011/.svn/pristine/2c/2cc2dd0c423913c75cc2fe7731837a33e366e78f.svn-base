/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.genome.GBS.switchgrass;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;

/**
 *
 * @author fl262
 */

//To do: debugging

public class GreedyLinkageGroupFinder {
	CorrRecord[] sigDoubleMatrix;
	CorrRecord[] bcpMatrix;
	Integer[][] linkageGroup;
	
	public GreedyLinkageGroupFinder(CorrRecord[] matrix) {
		this.sigDoubleMatrix = matrix;
		Arrays.sort(this.sigDoubleMatrix, new sortByIndex());
	}

	public void greedyFinder (int minEdge) {
		bcpMatrix = new CorrRecord[this.sigDoubleMatrix.length];
		System.arraycopy(this.sigDoubleMatrix, 0, bcpMatrix, 0, this.sigDoubleMatrix.length);
		int nullcount = 0;
		ArrayList<Integer[]> linkageGroupList = new ArrayList();
		while (bcpMatrix.length > minEdge && nullcount < 100) {
			Integer[] firstnode = new Integer[1];
			firstnode[0] = this.bcpMatrix[(int)Math.floor(Math.random() * bcpMatrix.length)].queryIndex;
			ArrayList<Integer> markerList = this.getMoreNode(firstnode, minEdge);
			if (markerList.isEmpty()) {
				nullcount++;
			}
			else {
				nullcount = 0;
				markerList.add(firstnode[0]);
				Integer[] onechrome = markerList.toArray(new Integer[markerList.size()]);
				System.out.println(onechrome.length);
				linkageGroupList.add(onechrome);
			}
		}
		linkageGroup = linkageGroupList.toArray(new Integer[linkageGroupList.size()][]);
	}

	public ArrayList<Integer> getMoreNode (Integer[] node, int minEdge) {
		ArrayList<Integer> allNextNodeList = new ArrayList();
		for (int i = 0; i < node.length; i++) {
			Integer[] nextEdge = findNextEdge(node[i], minEdge);
			if (nextEdge != null) {
				for (int j = 0; j < nextEdge.length; j++) {
					allNextNodeList.add(this.bcpMatrix[nextEdge[j]].hitIndex);
				}
				this.updateBcpMatrix(nextEdge);
			}
		}
		if (allNextNodeList.size() > minEdge) {
			Integer[] allNextNodeArray = allNextNodeList.toArray(new Integer[allNextNodeList.size()]);
			allNextNodeList.addAll(this.getMoreNode(allNextNodeArray, minEdge));
			return allNextNodeList;
		}
		else {
			return allNextNodeList;
		}
	}

	public Integer[] findNextEdge (int onenode, int minEdges) {
		CorrRecord key = new CorrRecord(onenode, -1, (float)0);
		int hit = Arrays.binarySearch(bcpMatrix, key, new sortByIndex());
		hit = -hit -1;
		if (hit == bcpMatrix.length) return null;
		ArrayList<Integer> nextEdgeList = new ArrayList();
		nextEdgeList.add(hit);
		while (hit < bcpMatrix.length - 1 && bcpMatrix[hit].queryIndex == bcpMatrix[hit + 1].queryIndex) {
			nextEdgeList.add(hit + 1);
			hit++;
		}
		Integer[] nextEdge = nextEdgeList.toArray(new Integer[nextEdgeList.size()]);
		if (nextEdge.length > minEdges) {
			return nextEdge;
		}
		return null;
	}

	public void updateBcpMatrix (Integer[] nextEdge) {
		ArrayList<Integer> deletedList = new ArrayList();
		for (int i = 0; i < nextEdge.length; i++) {
			deletedList.add(nextEdge[i]);
			CorrRecord key = new CorrRecord(this.bcpMatrix[nextEdge[i]].hitIndex, this.bcpMatrix[nextEdge[i]].queryIndex, (float)0);
			int hit = Arrays.binarySearch(this.bcpMatrix, key, new sortByIndex());
			deletedList.add(hit);
		}
		Integer[] deletedArray = deletedList.toArray(new Integer[deletedList.size()]);
		Arrays.sort(deletedArray);
		CorrRecord[] newArray = new CorrRecord[this.bcpMatrix.length - 2 * nextEdge.length];
		int count = 0;
		for (int i = 0; i < this.bcpMatrix.length; i++) {
			int hit = Arrays.binarySearch(deletedArray, i);
			if (hit < 0) {
				newArray[count] = this.bcpMatrix[i];
				count++;
			}
		}
		Arrays.sort(newArray, new sortByIndex());
		this.bcpMatrix = newArray;
		newArray = null;
	}

	public void writeLinkageGroup (String outfileS) {
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(outfileS), 65536);
			for (int i = 0; i < linkageGroup.length; i++) {
				bw.write(String.valueOf(i+1) + "\t" + String.valueOf(linkageGroup[i].length));
			}
		}
		catch (Exception e) {
			System.out.println("Error occurred while writing " + outfileS + " " + e.toString());
		}

	}
}
