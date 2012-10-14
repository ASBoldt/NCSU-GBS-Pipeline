/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.genome.GBS.repeat;

import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.image.BufferedImage;
import java.io.*;
import java.util.Arrays;
import javax.imageio.ImageIO;

/**
 *
 * @author fl262
 */
public class visualizeGlm {

	static File sourceGenoFileMerge = setQuery.sourceGenoFileMerge;
	static File SASResultF = processGlmSAS.desFile;
	static File figureFile = new File(processGlmSAS.SASResultDir, "figure.png");
	chromInfoArray cif;
	GeneticInfoArray gif;
	GlmRecSASArray grsa;

	public visualizeGlm() {
		getData();
		drawPic();
	}
	public visualizeGlm (File SASResultF, File figureFile) {
		this.SASResultF = SASResultF;
		this.figureFile = figureFile;
		getData();
		drawPic();
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
		}
	}

	public void drawPic() {
		int scaleX = 3, scaleY = 3;
		int imgWidth = 1600 * scaleX, imgHeight = 1200 * scaleY;
		int leftMargin = 100, topMargin = 50, intervalX = 150, barWidth = 10, lineLength = 7;
		float scale = imgHeight * 0.80f / cif.chroms[0].chromLength /scaleY;
		scaleData(scale);
		BufferedImage bi = new BufferedImage(imgWidth, imgHeight, BufferedImage.TYPE_4BYTE_ABGR);
		Graphics2D g2 = bi.createGraphics();
		g2.scale(scaleX, scaleY);
		Font font = new Font("chromName", Font.BOLD, 30);
		g2.setFont(font);
		for (int i = 0; i < cif.chromNum; i++) {
			int j = i + 1;
			g2.setPaint(Color.BLACK);
			g2.drawString(String.valueOf(j), leftMargin + i * intervalX, topMargin - 20);
			g2.setPaint(Color.RED);
			Rectangle chrom = new Rectangle(leftMargin + i * intervalX, topMargin, barWidth, cif.chroms[i].chromLength);
			g2.fill(chrom);
			g2.setPaint(Color.BLUE);
			Rectangle cent = new Rectangle(leftMargin + i * intervalX, topMargin + cif.chroms[i].centB, barWidth, cif.chroms[i].centE - cif.chroms[i].centB + 1);
			g2.fill(cent);
		}
		for (int i = 0; i <gif.markerNum; i++) {
			g2.setPaint(Color.BLACK);
			g2.drawLine(leftMargin + (gif.markers[i].chrom - 1) * intervalX - lineLength, topMargin + gif.markers[i].position, leftMargin + (gif.markers[i].chrom - 1) * intervalX - 1, topMargin + gif.markers[i].position);
		}
		String markerName = grsa.sigMarkers[0].markerName;
		font = new Font("markerName", Font.BOLD, 10);
		g2.setFont(font);
		int flag = 1, timesFlag = 10;
		int outMarker = 10, leftMove = 70;
		for (int i = 0; i < grsa.sigMarkerNum; i++) {
			g2.setPaint(Color.BLUE);
			if (grsa.sigMarkers[i].markerName.equals(markerName)) {
				flag++;
			}
			else {
				g2.drawLine(leftMargin + (grsa.sigMarkers[i-1].chrom - 1) * intervalX + barWidth + 1, topMargin + grsa.sigMarkers[i-1].position, leftMargin + (grsa.sigMarkers[i-1].chrom - 1) * intervalX + barWidth + 1 + timesFlag * flag, topMargin + grsa.sigMarkers[i-1].position);
				if (flag >= outMarker) {
					g2.setPaint(Color.black);
					g2.drawString(grsa.sigMarkers[i-1].markerName, leftMargin + (grsa.sigMarkers[i-1].chrom - 1) * intervalX - leftMove, topMargin + grsa.sigMarkers[i-1].position);
				}
				markerName = grsa.sigMarkers[i].markerName;
				flag = 1;
			}
			if (i == grsa.sigMarkerNum - 1) {
				g2.drawLine(leftMargin + (grsa.sigMarkers[i].chrom - 1) * intervalX + barWidth + 1, topMargin + grsa.sigMarkers[i].position, leftMargin + (grsa.sigMarkers[i].chrom - 1) * intervalX + barWidth + 1 + timesFlag * flag, topMargin + grsa.sigMarkers[i].position);
				if (flag >= outMarker) {
					g2.setPaint(Color.black);
					g2.drawString(grsa.sigMarkers[i-1].markerName, leftMargin + (grsa.sigMarkers[i-1].chrom - 1) * intervalX - leftMove, topMargin + grsa.sigMarkers[i-1].position);
				}
			}
		}
		g2.drawLine(leftMargin + (cif.chromNum - 1) * intervalX, 4 * topMargin + cif.chroms[cif.chromNum-1].chromLength, leftMargin + (cif.chromNum - 1) * intervalX + 20 * timesFlag - 1,  4 * topMargin + cif.chroms[cif.chromNum-1].chromLength);
		font = new Font("couting", Font.BOLD, 15);
		g2.setFont(font);
		g2.setPaint(Color.black);
		g2.drawString("20", leftMargin + (cif.chromNum - 1) * intervalX, 4 * topMargin + cif.chroms[cif.chromNum-1].chromLength - 5);
		try {
			bi.flush();
			ImageIO.write(bi, "png", figureFile);
		} catch (Exception e) {
			System.out.println(figureFile.toString());
		}
	}

	public void getData() {
		cif = new chromInfoArray();
		gif = new GeneticInfoArray(sourceGenoFileMerge);
		Arrays.sort(gif.markers, new sortByChrom());
		grsa = new GlmRecSASArray(SASResultF, true);
		Arrays.sort(grsa.sigMarkers, new sortByMarker());
	}

	public static void main(String[] agrs) {
		visualizeGlm vg = new visualizeGlm();
	}
}

