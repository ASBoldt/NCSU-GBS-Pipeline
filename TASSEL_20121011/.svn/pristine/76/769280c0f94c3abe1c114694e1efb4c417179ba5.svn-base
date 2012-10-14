/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.genome.GBS.repeat;

import java.awt.Color;
import java.awt.Container;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseMotionAdapter;
import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JTextField;

/**
 *
 * @author fl262
 */
public class GUIGlm extends JFrame {

	static File sourceGenoFileMerge = setQuery.sourceGenoFileMerge;
	//static File SASResultF = processGlmSAS.desFile;
	static File SASResultF = new File("E:/repeat/4QTL/R_4SAS/annotation/chloroplast_Glm.txt");
	static File blastF = new File(Blastn2Table.desFileS);
	chromInfoArray cif;
	GeneticInfoArray gif;
	GlmRecSASArray grsa, OriGrsa;
	blastRecArray bra;
	Integer[] indexDifMarker;
	Integer[] countDifMarker;
	Integer[] PositionsArray;
	int fWidth, fHeight, topMargin, leftMargin, rightMargin, labelWidth;
	chromL[] chromLs;
	String[] sortedSeqName, GlmSeqName;
	JTextField[] infoText = new JTextField[5];
	JLabel[] seqLabel;
	int clickY = -1, clickIndex = -1, moveY, moveIndex;
	HashMap<Integer, Integer>[] Y2MarkerNum;

	public GUIGlm() {
		getData();
		guiDesign();
	}

	public void guiDesign() {
		fWidth = 1650;
		fHeight = 1000;
		topMargin = 50;
		leftMargin = 30;
		rightMargin = 50;
		labelWidth = (fWidth - leftMargin - rightMargin) / cif.chromNum;
		this.setSize(fWidth, fHeight);
		this.setTitle("Result of Glmselect");
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		Container JFC = this.getContentPane();
		//JFC.setBackground(Color.WHITE);
		JFC.setLayout(null);
		float scale = fHeight * 0.80f / cif.chroms[0].chromLength;
		scaleData(scale);

		chromLs = new chromL[cif.chromNum];
		for (int i = 0; i < cif.chromNum; i++) {
			chromLs[i] = new chromL(leftMargin + i * labelWidth, topMargin, labelWidth, i);
			this.add(chromLs[i]);
		}
		for (int i = 0; i < infoText.length; i++) {
			int tWidth = 100, tHeight = 20;
			infoText[i] = new JTextField();
			infoText[i].setSize(tWidth - 10, tHeight);
			infoText[i].setLocation(leftMargin + 6 * labelWidth + tWidth * i, cif.chroms[6].chromLength + 2 * topMargin);
			this.add(infoText[i]);
		}
		infoText[1].addMouseListener(new mouseClickListener());
	}

	class chromL extends JLabel {

		int x, topMargin, width, chrIndex, diameter = 12;
		int barWidth = 10;

		chromL(int x, int topMargin, int width, int chrIndex) {
			this.x = x;
			this.topMargin = topMargin;
			this.width = width;
			this.chrIndex = chrIndex;
			this.addMouseListener(new mouseClickListener());
			this.addMouseMotionListener(new mouseMotionListener());
			//setBackground(Color.CYAN);
			setSize(width, cif.chroms[chrIndex].chromLength + topMargin);
			setLocation(x, 0);
			setOpaque(true);
		}

		@Override
		public void paint(Graphics g) {
			super.paint(g);
			Font font = new Font("markerName", Font.BOLD, 25);
			g.setFont(font);
			g.drawString(String.valueOf(chrIndex + 1), width / 2 - barWidth, topMargin - 10);
			int barX = (width - barWidth) / 2;
			g.setColor(Color.red);
			g.fillRect(barX, 0 + topMargin, barWidth, cif.chroms[chrIndex].chromLength);
			g.setColor(Color.darkGray);
			g.fillRect(barX, cif.chroms[chrIndex].centB + topMargin, barWidth, cif.chroms[chrIndex].centE - cif.chroms[chrIndex].centB + 1);
			byte[] tempChr = new byte[gif.markerNum];
			for (int i = 0; i < gif.markerNum; i++) {
				tempChr[i] = gif.markers[i].chrom;
			}
			g.setColor(Color.BLACK);
			for (int i = 0; i < gif.markerNum; i++) {
				if (gif.markers[i].chrom == chrIndex + 1) {
					g.drawLine(barX - 5, gif.markers[i].position + topMargin, barX - 1, gif.markers[i].position + topMargin);
				}
			}
			g.setColor(Color.blue);
			for (int i = 0; i < indexDifMarker.length; i++) {
				if (grsa.sigMarkers[indexDifMarker[i]].chrom == chrIndex + 1) {
					g.drawLine(barX + barWidth, grsa.sigMarkers[indexDifMarker[i]].position + topMargin, barX + barWidth + countDifMarker[i] * 2, grsa.sigMarkers[indexDifMarker[i]].position + topMargin);
				}
			}
			g.setColor(Color.YELLOW);
			if (PositionsArray != null) {
				for (int i = 0; i < PositionsArray.length; i++) {
					if (OriGrsa.sigMarkers[PositionsArray[i]].chrom == chrIndex + 1) {
						g.drawOval(barX - 18, OriGrsa.sigMarkers[PositionsArray[i]].position + topMargin - diameter / 2, diameter, diameter);
						g.fillOval(barX - 18, OriGrsa.sigMarkers[PositionsArray[i]].position + topMargin - diameter / 2, diameter, diameter);
					}
				}
			}
			g.setColor(Color.GRAY);
			if (clickIndex == chrIndex) {
				g.drawOval(barX - 18, clickY + topMargin - diameter / 2, diameter, diameter);
				g.fillOval(barX - 18, clickY + topMargin - diameter / 2, diameter, diameter);
			}
		}
	}

	public void showSeqLabel(int index) {
		if (seqLabel != null) {
			for (int i = 0; i < seqLabel.length; i++) {
				this.remove(seqLabel[i]);
				seqLabel[i] = null;
			}
			seqLabel = null;
			System.gc();
		}
		int iniX = infoText[0].getX(), iniY = infoText[0].getY() + 30, lWidth = 50, lHeight = 20;
		seqLabel = new JLabel[countDifMarker[index]];
		for (int i = 0; i < seqLabel.length; i++) {
			int items = 10;
			int nx = i % items;
			int ny = (i - nx) / items;
			seqLabel[i] = new JLabel(grsa.sigMarkers[indexDifMarker[index] + i].traitName);
			seqLabel[i].setOpaque(true);
			seqLabel[i].setSize(lWidth - 10, lHeight);
			seqLabel[i].setLocation(iniX + nx * lWidth, iniY + ny * lHeight);
			this.remove(seqLabel[i]);
			this.add(seqLabel[i]);
			seqLabel[i].addMouseListener(new mouseClickListener());
		}
		this.repaint();
	}

	public void showTextClick() {
		for (int i = 0; i < indexDifMarker.length; i++) {
			if (grsa.sigMarkers[indexDifMarker[i]].position == clickY && grsa.sigMarkers[indexDifMarker[i]].chrom == clickIndex + 1) {
				infoText[1].setText("1");
				infoText[2].setText(String.valueOf(grsa.sigMarkers[indexDifMarker[i]].chrom));
				infoText[3].setText(grsa.sigMarkers[indexDifMarker[i]].markerName);
				infoText[4].setText(String.valueOf(countDifMarker[i]));
				showSeqLabel(i);
				break;
			}
		}
	}

	public void showTextMotion() {
		infoText[0].setText(String.valueOf(Y2MarkerNum[moveIndex].get(moveY)));
	}

	public void showPosition(int index) {
		ArrayList<Integer> Positions = new ArrayList();
		int hit = Arrays.binarySearch(GlmSeqName, seqLabel[index].getText());
		while (GlmSeqName[hit].equals(seqLabel[index].getText())) {
			hit--;
		}
		hit++;
		for (int i = hit; i < GlmSeqName.length; i++) {
			if (GlmSeqName[i].equals(seqLabel[index].getText())) {
				Positions.add(i);
			}
			else {
				break;
			}
		}
		PositionsArray = Positions.toArray(new Integer[Positions.size()]);
		System.out.println(PositionsArray.length);
		repaint();
	}

	public void changeMarkerClick() {
		int total;
		if (Y2MarkerNum[clickIndex].containsKey(clickY)) {
			total = Y2MarkerNum[clickIndex].get(clickY);
			int status = Integer.valueOf(infoText[1].getText());
			if (status == total) {
				status = 0;
			}
			for (int i = 0; i < indexDifMarker.length; i++) {
				if (grsa.sigMarkers[indexDifMarker[i]].position == clickY && grsa.sigMarkers[indexDifMarker[i]].chrom == clickIndex + 1) {
					infoText[2].setText(String.valueOf(grsa.sigMarkers[indexDifMarker[i + status]].chrom));
					infoText[3].setText(grsa.sigMarkers[indexDifMarker[i + status]].markerName);
					infoText[4].setText(String.valueOf(countDifMarker[i + status]));
					showSeqLabel(i + status);
					break;
				}
			}
			status++;
			infoText[1].setText(String.valueOf(status));
		}
	}

	class mouseClickListener extends MouseAdapter {

		@Override
		public void mouseClicked(MouseEvent e) {
			for (int i = 0; i < chromLs.length; i++) {
				if (e.getSource() == chromLs[i]) {
					clickY = e.getY() - topMargin;
					clickIndex = i;
					PositionsArray = null;
					repaint();
					showTextClick();
					break;
				}
			}
			if (seqLabel != null) {
				for (int i = 0; i < seqLabel.length; i++) {
					if (e.getSource() == seqLabel[i]) {
						showPosition(i);
						new annotation(i);
					}
				}
			}
			if (e.getSource() == infoText[1]) {
				changeMarkerClick();
			}
		}
	}

	class mouseMotionListener extends MouseMotionAdapter {

		@Override
		public void mouseMoved(MouseEvent e) {
			for (int i = 0; i < chromLs.length; i++) {
				if (e.getSource() == chromLs[i]) {
					moveY = e.getY() - topMargin;
					moveIndex = i;
					showTextMotion();
					break;
				}
			}
		}
	}

	class annotation extends JFrame {

		int AFWidth = 1000;
		int AFHeight = 600;
		int AFtopMargin = 50;
		int AFleftMargin = 30;
		int index;
		int maxLine = 20;
		int interval = 20;

		annotation(int index) {
			this.setSize(AFWidth, AFHeight);
			this.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
			this.setVisible(true);
			this.index = index;
		}

		@Override
		public void paint(Graphics g) {
			super.paint(g);
			int hit = Arrays.binarySearch(sortedSeqName, seqLabel[index].getText());
			if (hit > 0) {
				while (sortedSeqName[hit].equals(seqLabel[index].getText())) {
					hit--;
				}
				hit++;
				this.setTitle("Annotation of " + bra.recs[hit].item[0]);
				for (int i = hit; i < bra.recNum; i++) {
					if (!sortedSeqName[i].equals(seqLabel[index].getText()) || i > hit + maxLine) {
						break;
					}
					g.drawString(bra.recs[i].item[2], AFleftMargin, AFtopMargin + (i - hit) * interval);
				}
			} else {
				this.setTitle("Annotation of seq " + seqLabel[index].getText() + " is not found");
			}
		}
	}

	public void getData() {
		cif = new chromInfoArray();
		gif = new GeneticInfoArray(sourceGenoFileMerge);
		Arrays.sort(gif.markers, new sortByChrom());
		grsa = new GlmRecSASArray(SASResultF, true);
		OriGrsa = new GlmRecSASArray(SASResultF, true);
		Arrays.sort(grsa.sigMarkers, new sortByMarker());
		Arrays.sort(OriGrsa.sigMarkers);
		bra = new blastRecArray(blastF);
		Arrays.sort(bra.recs);
		sortedSeqName = new String[bra.recNum];
		for (int i = 0; i < bra.recNum; i++) {
			sortedSeqName[i] = bra.recs[i].item[0].replaceAll("\\|.+", "");
		}
		ArrayList<Integer> alC = new ArrayList();
		ArrayList<Integer> alI = new ArrayList();
		alI.add(0);
		String temp = grsa.sigMarkers[0].markerName;
		int count = 0;
		for (int i = 1; i < grsa.sigMarkerNum; i++) {
			if (temp.equals(grsa.sigMarkers[i].markerName)) {
				count++;
			} else {
				alC.add(++count);
				count = 0;
				temp = grsa.sigMarkers[i].markerName;
				alI.add(i);
			}
		}
		alC.add(count);
		indexDifMarker = alI.toArray(new Integer[alI.size()]);
		countDifMarker = alC.toArray(new Integer[alC.size()]);
	}

	public void scaleData(float scale) {
		for (int i = 0; i < cif.chromNum; i++) {
			cif.chroms[i].chromLength = Math.round(scale * cif.chroms[i].chromLength);
			cif.chroms[i].centB = Math.round(scale * cif.chroms[i].centB);
			cif.chroms[i].centE = Math.round(scale * cif.chroms[i].centE);
		}
		for (int i = 0; i < gif.markerNum; i++) {
			gif.markers[i].position = Math.round(scale * gif.markers[i].position);
		}
		for (int i = 0; i < grsa.sigMarkerNum; i++) {
			grsa.sigMarkers[i].position = Math.round(scale * grsa.sigMarkers[i].position);
			OriGrsa.sigMarkers[i].position = Math.round(scale * OriGrsa.sigMarkers[i].position);
		}
		GlmSeqName = new String[OriGrsa.sigMarkerNum];
		for (int i = 0; i < OriGrsa.sigMarkerNum; i++) {
			GlmSeqName[i] = OriGrsa.sigMarkers[i].traitName;
		}
		
		Y2MarkerNum = new HashMap[cif.chromNum];
		for (int i = 0; i < cif.chromNum; i++) {
			Y2MarkerNum[i] = new HashMap();
		}
		int temp = grsa.sigMarkers[indexDifMarker[0]].position;
		int count = 1;
		for (int i = 1; i < indexDifMarker.length; i++) {
			if (temp == grsa.sigMarkers[indexDifMarker[i]].position) {
				count++;
			} else {
				Y2MarkerNum[grsa.sigMarkers[indexDifMarker[i - 1]].chrom - 1].put(temp, count);
				temp = grsa.sigMarkers[indexDifMarker[i]].position;
				count = 1;
			}
		}
		Y2MarkerNum[grsa.sigMarkers[grsa.sigMarkerNum - 1].chrom - 1].put(temp, count);
	}

	public static void main(String args[]) {
		new GUIGlm().setVisible(true);
	}
}
