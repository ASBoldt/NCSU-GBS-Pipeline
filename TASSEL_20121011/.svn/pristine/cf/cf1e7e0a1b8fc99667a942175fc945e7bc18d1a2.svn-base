/*
 * LinkageDisequilibriumPlugin.java
 *
 * Created on December 22, 2006, 5:02 PM
 *
 */
package net.maizegenetics.baseplugins;

import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginEvent;
import net.maizegenetics.pal.popgen.LinkageDisequilibrium;

import javax.swing.*;

import java.awt.*;
import java.awt.event.ActionEvent;

import java.net.URL;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 *
 * @author Ed Buckler
 */
public class LinkageDisequilibriumPlugin extends AbstractPlugin {

//    boolean isRapidAnalysis = true;
//    int permutationNumber = 1000;
    int windowSize = 50;
    LinkageDisequilibrium.testDesign LDType = LinkageDisequilibrium.testDesign.SlidingWindow;
    int testSite = -1;

    /** Creates a new instance of LinkageDisequilibriumPlugin */
    public LinkageDisequilibriumPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    public DataSet performFunction(DataSet input) {

        try {
            List<Datum> alignInList = input.getDataOfType(Alignment.class);
            if (alignInList.size() < 1) {
                JOptionPane.showMessageDialog(getParentFrame(), "Invalid selection.  Please select sequence or marker alignment.");
                return null;
            }

            Datum align = alignInList.get(0);

            if (isInteractive()) {
//                LinkageDiseqDialog myDialog = new LinkageDiseqDialog(isRapidAnalysis, permutationNumber);
                LinkageDiseqDialog myDialog = new LinkageDiseqDialog(((Alignment)align.getData()).getSiteCount());
                myDialog.setLocationRelativeTo(getParentFrame());
                myDialog.setVisible(true);
                if (myDialog.isCancel()) {
                    return null;
                }
//                isRapidAnalysis = myDialog.getRapidLDAnalysis();
                windowSize = myDialog.getWindowSize();
                LDType = myDialog.getLDType();
            }

            List result = new ArrayList();
            Iterator<Datum> itr = alignInList.iterator();
            while (itr.hasNext()) {
                DataSet tds = null;
                Datum current = itr.next();
                tds = processDatum(current);
                if (tds != null) {
                    result.add(tds);
                    fireDataSetReturned(new PluginEvent(tds, LinkageDisequilibriumPlugin.class));
                }
            }

            return DataSet.getDataSet(result, this);
        } finally {
            fireProgress(100);
        }

    }

    private DataSet processDatum(Datum input) {
        Alignment aa = (Alignment) input.getData();
//        net.maizegenetics.pal.popgen.LinkageDisequilibrium theLD = new net.maizegenetics.pal.popgen.LinkageDisequilibrium(aa, permutationNumber, windowSize, isRapidAnalysis, LDType, testSite, this);
        net.maizegenetics.pal.popgen.LinkageDisequilibrium theLD = new net.maizegenetics.pal.popgen.LinkageDisequilibrium(aa, windowSize, LDType, testSite, this);
        try {
            theLD.run();
            Datum td = new Datum("LD:" + input.getName(), theLD, "LD analysis");
            DataSet tds = new DataSet(td, this);
            return tds;
        } catch (Exception e) {
            e.printStackTrace();
            StringBuilder builder = new StringBuilder();
            builder.append("Unable to run Linkage Disequilibrium analysis ");
            builder.append(e.getMessage());
            String str = builder.toString();
            JOptionPane.showMessageDialog(getParentFrame(), str);
        }
        return null;
    }

//    public boolean isRapidAnalysis() {
//        return isRapidAnalysis;
//    }

//    public void setRapidAnalysis(boolean rapidAnalysis) {
//        isRapidAnalysis = rapidAnalysis;
//    }

//    public int isPermutationNumber() {
//        return permutationNumber;
//    }

//    public void setPermutationNumber(int permutationNumber) {
//        this.permutationNumber = permutationNumber;
//    }

    public void setLDType(LinkageDisequilibrium.testDesign type) {
        LDType = type;
    }

    public LinkageDisequilibrium.testDesign getLDType() {
        return LDType;
    }

    public void setWinSize(int winSize) {
        windowSize = winSize;
    }

    public int getWinSize() {
        return windowSize;
    }

    public void setTestSite(int site) {
        testSite = site;
    }

    public int getTestSite() {
        return testSite;
    }

    /**
     * Icon for this plugin to be used in buttons, etc.
     *
     * @return ImageIcon
     */
    public ImageIcon getIcon() {
        URL imageURL = LinkageDisequilibriumPlugin.class.getResource("images/LDPlot.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    /**
     * Button name for this plugin to be used in buttons, etc.
     *
     * @return String
     */
    public String getButtonName() {
        return "Link. Diseq.";
    }

    /**
     * Tool Tip Text for this plugin
     *
     * @return String
     */
    public String getToolTipText() {
        return "Linkage Disequilibrium";
    }
}

/**
 * Title:        TASSEL
 * Description:  A java program to deal with diversity
 * Copyright:    Copyright (c) 2000
 * Company:      USDA-ARS/NCSU
 * @author Ed Buckler
 * @version 1.0
 */
class LinkageDiseqDialog extends JDialog {

//    private boolean rapidPermutations = true;
    // private int numberPermutations = 1000;
    private boolean myRunAnalysis = false;
    // private JTextField permNumberTextField = new JTextField();
    // private JLabel jLabel1 = new JLabel();
    private int myTestSite = 0;
//    private String myAlignmentForSiteList;

    private int myNumSites;

    private JPanel myPanel = new JPanel();

    private JPanel myLDSelectionPanel = new JPanel();
    private JLabel myLDTypeLabel = new JLabel("Select LD type: ");
    private JComboBox myLDType;
//    private JCheckBox myAccumulativeResultsBox = new JCheckBox("Accumulate R2 Results");

    private JPanel myLDOptionsPanel = new JPanel();
    private JLabel myFullMatrixLabel = new JLabel();
    private JTextField myWindowSizeTextField = new JTextField();
    private JLabel myWindowSizeLabel = new JLabel("LD Window Size: ");
    private int myWindowSize = 50;

    private JPanel myButtonsPanel = new JPanel();
    private JButton myRunButton = new JButton("Run");
    private JButton myCloseButton = new JButton("Close");

//    public LinkageDiseqDialog(int numSites, int numberPermutations) {
    public LinkageDiseqDialog(int numSites) {
        super((Frame) null, "Linkage Disequilibrium", true);
        myNumSites = numSites;
        // numberPermutations = numberPermutations;
        try {
            jbInit();
            pack();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    private void jbInit() throws Exception {

        String[] ldTypes = {"Full Matrix", "Sliding Window"};

        myLDType = new JComboBox(ldTypes);
        myLDType.setSelectedIndex(1);

        myFullMatrixLabel.setText("Full LD with " + (myNumSites*(myNumSites/2-1)) + " comparisons.");

        myLDSelectionPanel.setLayout(new GridBagLayout());
        myLDSelectionPanel.add(myLDTypeLabel, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0, GridBagConstraints.LINE_END, GridBagConstraints.NONE, new Insets(10, 7, 0, 0), 0, 0));
        myLDSelectionPanel.add(myLDType, new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0, GridBagConstraints.LINE_START, GridBagConstraints.NONE, new Insets(10, 0, 0, 0), 0, 0));

        myLDType.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                ldType_actionPerformed(e);
            }
        });

        myLDOptionsPanel.setLayout(new GridBagLayout());
        myLDOptionsPanel.add(myFullMatrixLabel, new GridBagConstraints(0, 0, 2, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(10, 0, 0, 0), 0, 0));
        myFullMatrixLabel.setVisible(false);

        myWindowSizeTextField.setText("" + myWindowSize);
        myLDOptionsPanel.add(myWindowSizeLabel, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0, GridBagConstraints.LINE_END, GridBagConstraints.NONE, new Insets(10, 7, 0, 0), 0, 0));
        myLDOptionsPanel.add(myWindowSizeTextField, new GridBagConstraints(1, 0, 1, 1, 1.0, 0.0, GridBagConstraints.LINE_START, GridBagConstraints.HORIZONTAL, new Insets(10, 0, 0, 7), 0, 0));

        myButtonsPanel.setLayout(new GridBagLayout());
        myButtonsPanel.add(myRunButton, new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0, GridBagConstraints.PAGE_END, GridBagConstraints.NONE, new Insets(0, 150, -5, 0), 0, 0));
        myButtonsPanel.add(myCloseButton, new GridBagConstraints(3, 0, 1, 1, 0.0, 0.0, GridBagConstraints.PAGE_END, GridBagConstraints.NONE, new Insets(0, 0, -5, 0), 0, 0));

        myPanel.setLayout(new GridLayout(3, 1));
        myPanel.add(myLDSelectionPanel);
        myPanel.add(myLDOptionsPanel);
        myPanel.add(myButtonsPanel);

        myRunButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                runButton_actionPerformed(e);
            }
        });

        myCloseButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                closeButton_actionPerformed(e);
            }
        });
        getContentPane().add(myPanel);
        getContentPane().setPreferredSize(new Dimension(300, 130));
    }

    private void hideOptions() {
        myWindowSizeTextField.setVisible(false);
        myWindowSizeLabel.setVisible(false);
    }

    public LinkageDisequilibrium.testDesign getLDType() {
        if (myLDType.getSelectedIndex() == 0) {
            return LinkageDisequilibrium.testDesign.All;
        } else if (myLDType.getSelectedIndex() == 1) {
            return LinkageDisequilibrium.testDesign.SlidingWindow;
        } else {
            throw new IllegalStateException("LinkageDisequilibriumPlugin: getLDType: No known LD Type selected.");
        }
    }

    /** Returns whether the run button was chosen*/
    public boolean isRunAnalysis() {
        return myRunAnalysis;
    }

    /** Returns whether the run button was chosen*/
    public boolean isCancel() {
        return !myRunAnalysis;
    }

    /** Return the window size */
    public int getWindowSize() {
        return myWindowSize;
    }

    public int getTestSite() {
        return myTestSite;
    }

//    public int getNumAccumulateIntervals() {
//        return myNumAccumulativeInterval;
//    }

    void ldType_actionPerformed(ActionEvent e) {
        if ( myLDType.getSelectedIndex() == 0 ) {
            myFullMatrixLabel.setVisible(true);
            myWindowSizeLabel.setVisible(false);
            myWindowSizeTextField.setVisible(false);
        } else if ( myLDType.getSelectedIndex() == 1 ) {
            myFullMatrixLabel.setVisible(false);
            myWindowSizeLabel.setVisible(true);
            myWindowSizeTextField.setVisible(true);
        }
    }

    void runButton_actionPerformed(ActionEvent e) {

        //        try {
        //            numberPermutations = Integer.parseInt(permNumberTextField.getText());
        //        } catch (Exception ee) {
        //            permNumberTextField.setText("Set Integer");
        //            return;
        //        }

        if (getLDType() == LinkageDisequilibrium.testDesign.SlidingWindow) {
            try {
                myWindowSize = Integer.parseInt(myWindowSizeTextField.getText());
            } catch (Exception ee) {
                myWindowSizeTextField.setText("Set Integer");
                return;
            }
        }

        myRunAnalysis = true;
        setVisible(false);
    }

    void closeButton_actionPerformed(ActionEvent e) {
        myRunAnalysis = false;
        setVisible(false);
    }
}
