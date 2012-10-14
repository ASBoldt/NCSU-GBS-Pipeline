/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.genome.GBS.switchgrass;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import net.maizegenetics.pal.tree.NeighborJoiningTree;
import net.maizegenetics.pal.distance.DistanceMatrix;
import net.maizegenetics.pal.ids.SimpleIdGroup;
import net.maizegenetics.pal.tree.Node;
import net.maizegenetics.pal.tree.NodeUtils;
import net.maizegenetics.pal.tree.ReadTree;
import net.maizegenetics.pal.tree.SimpleTree;
import net.maizegenetics.pal.tree.TreeParseException;
/**
 *
 * @author fl262
 */
public class LinkageGroupUtils {
	LinkageGroup[] lGroup;
	SimpleTree theNJT;

	public LinkageGroupUtils () {};

	public LinkageGroupUtils (Matrix dmatrixFloat) {
		this.getLinkageGroupByMatrix(dmatrixFloat);
	}

	public LinkageGroupUtils (String infileS) {
		this.readLinkageGroup(infileS);
	}

	public void getLinkageGroupByMatrix (Matrix dmatrixFloat) {
		DistanceMatrix dmatrix = getDistanceMatrixDouble(dmatrixFloat);
		System.out.println("Calculating NJ tree");
		theNJT = new NeighborJoiningTree(dmatrix);
		getLinkageGroupByTree(theNJT);
	}

	public void getLinkageGroupByTree (String infileS) {
		//this function doesn't work, because the pal.tree.ReadTree does not function well
		try {
			this.readNJTree(infileS);
			this.getLinkageGroupByTree(theNJT);
		}
		catch (TreeParseException e) {
			System.out.println("Error occurred while reading NJtree " + infileS + " " + e.toString());
			System.exit(-1);
		}
	}

	public void getLinkageGroupByTree (SimpleTree theNJT) {
		ArrayList<Node> rootNodeList = new ArrayList();
		ArrayList<Node> leafNodeList = new ArrayList();
		ArrayList<Node> topLeafNodeList = new ArrayList();
		for (int i = 0; i < theNJT.getInternalNodeCount(); i++) {
            Node inode = theNJT.getInternalNode(i);
            if(inode.getNodeHeight() < 0.48) leafNodeList.add(inode);
			else rootNodeList.add(inode);
        }
		Node[] leafNodeArray = leafNodeList.toArray(new Node[leafNodeList.size()]);
		Node[] rootNodeArray = rootNodeList.toArray(new Node[rootNodeList.size()]);
		for (int i = 0; i < leafNodeArray.length; i++) {
			Node parent = leafNodeArray[i].getParent();
			for (int j = 0; j < rootNodeArray.length; j++) {
				if (parent == rootNodeArray[j]) {
					topLeafNodeList.add(leafNodeArray[i]);
					break;
				}
			}
		}
		Node[] topLeafNodeArray = topLeafNodeList.toArray(new Node[topLeafNodeList.size()]);
		lGroup = new LinkageGroup[topLeafNodeArray.length];
		for (int i = 0; i < topLeafNodeArray.length; i++) {
			Node[] leafNode = NodeUtils.getExternalNodes(topLeafNodeArray[i]);
			int[] marker = new int[leafNode.length];
			for (int j = 0; j < leafNode.length; j++) {
				marker[j] = Integer.valueOf(leafNode[j].getIdentifier().getName());
			}
			lGroup[i] = new LinkageGroup(marker);
		}
		System.out.println("Unordered linkage groups are set up by NJ tree");
	}

	public int getLinkageGroupNumber () {
		return this.lGroup.length;
	}

	public int[] getOrderedLinkageGroupSize () {
		int[] size = new int[this.lGroup.length];
		for (int i = 0; i < size.length; i++) {
			size[i] = this.lGroup[i].getSize();
		}
		Arrays.sort(size);
		return size;
	}

	public DistanceMatrix getDistanceMatrixDouble (Matrix dmatrixFloat) {
		double[][] dmatrixDouble = new double[dmatrixFloat.value.length][dmatrixFloat.value.length];
		for (int i = 0; i < dmatrixFloat.value.length; i++) {
			for (int j = 0; j < dmatrixFloat.value.length; j++) {
				dmatrixDouble[i][j] = (double)dmatrixFloat.value[i][j];
			}
		}
		SimpleIdGroup idg = new SimpleIdGroup(dmatrixFloat.getMatrixSize(), true);
		DistanceMatrix dmatrix = new DistanceMatrix(dmatrixDouble, idg);
		return dmatrix;
	}

	public void orderLinkageGroup (Matrix r2Matrix, Boolean ifseed) {
		for (int i = 0; i < lGroup.length; i++) {
			lGroup[i].orderMarkers(r2Matrix, ifseed);
		}
		System.out.println("Markers on linkage groups are ordered");
	}

	public void deleteSmallSizeGroups (int size) {
		ArrayList<LinkageGroup> newList = new ArrayList();
		for (int i = 0; i < this.lGroup.length; i++) {
			if (lGroup[i].getSize() > size) {
				newList.add(lGroup[i]);
			}
		}
		LinkageGroup[] newArray = newList.toArray(new LinkageGroup[newList.size()]);
		this.lGroup = newArray;
	}

	public void sortLinkageGroupBySize () {
		Arrays.sort(lGroup);
	}

	public void mergeLinkageGroup () {


	}

	public void writeNJTree (String outfileS) {
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(outfileS), 65536);
			bw.write(theNJT.toString());
			bw.flush();
			bw.close();
		}
		catch (Exception e) {
			System.out.println("Error occurred while writing NJtree " + outfileS + " " + e.toString());
		}
	}

	public void writeLinkageGroup (String outfileS) {
		try {
			DataOutputStream dos =new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outfileS),65536));
			dos.writeInt(this.lGroup.length);
			for (int i = 0; i < this.lGroup.length; i++) {
				lGroup[i].writeLinkageGroup(dos);
			}
			dos.flush();
			dos.close();
		}
		catch (Exception e) {
			System.out.println("Error occurred while writing LinkageGroup " + outfileS + " " + e.toString());
		}
		System.out.println(String.valueOf(this.lGroup.length) + " linkage groups were witten");
	}

	public void readNJTree (String infileS) throws TreeParseException {
		try {
			theNJT = new ReadTree(infileS);
		}
		catch (Exception e) {
			System.out.println("Error occurred while reading NJtree " + infileS + " " + e.toString());
		}
	}

	public void readLinkageGroup (String infileS) {
		try {
			DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(infileS), 65536));
			this.lGroup = new LinkageGroup[dis.readInt()];
			for (int i = 0; i < this.lGroup.length; i++) {
				lGroup[i] = new LinkageGroup();
				lGroup[i].readLinkageGroup(dis);
			}
		}
		catch (Exception e) {
			System.out.println("Error occurred while reading LinkageGroup " + infileS + " " + e.toString());
		}
	}
}
