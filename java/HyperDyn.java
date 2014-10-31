import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;

public class HyperDyn extends Basic implements ActionListener {
        /*
*********************************************************************

dl_poly/java Bias Potential Dynamics class

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
    private static MakeControl home;
    private static HyperDyn job;
    private static JTextField tebias,tvmin,target,block,black,tdelt;
    private static JTextField ttlow,track,tcatch,spring,opttol,nnebs;
    private static JTextField basin1,basin2;
    private static JComboBox<String> hypopt,tunits,optkey;
    private static JCheckBox tgoneb,tpath;
    private static JButton close;


    // Define the Graphical User Interface

    public HyperDyn() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        super();
        setTitle("Hyperdynamics");
        int n=0;

        getContentPane().setForeground(art.fore);
        getContentPane().setBackground(art.back);
        setDefaultCloseOperation(DO_NOTHING_ON_CLOSE);
        setFont(fontMain);
        GridBagLayout grd = new GridBagLayout();
        GridBagConstraints gbc = new GridBagConstraints();
        getContentPane().setLayout(grd);

        gbc.fill=GridBagConstraints.BOTH;

        // Panel label

        fix(new JLabel("Select options:",JLabel.LEFT),grd,gbc,0,n++,3,1);

        // Hyperdynamics options

        fix(new JLabel("Hyperdynamics options:",JLabel.LEFT),grd,gbc,0,n,2,1);
        hypopt = new JComboBox<String>();
        hypopt.setForeground(art.scrf);
        hypopt.setBackground(art.scrn);
        hypopt.addItem("BPD");
        hypopt.addItem("TAD");
        hypopt.addItem("NEB");
        fix(hypopt,grd,gbc,2,n++,1,1);

        fix(new JLabel(" ",JLabel.LEFT),grd,gbc,0,n++,3,1);
        fix(new JLabel("Controls common to BPD, TAD & NEB:",JLabel.LEFT),grd,gbc,0,n++,3,1);

        // Energy units

        fix(new JLabel("Energy Units:",JLabel.LEFT),grd,gbc,0,n,2,1);
        tunits = new JComboBox<String>();
        tunits.setForeground(art.scrf);
        tunits.setBackground(art.scrn);
        tunits.addItem("DL_POLY");
        tunits.addItem("e_volt");
        tunits.addItem("k_cal");
        tunits.addItem("k_joule");
        tunits.addItem("Kelvin");
        fix(tunits,grd,gbc,2,n++,1,1);

        // NEB spring force constant

        fix(new JLabel("Spring force Constant",JLabel.LEFT),grd,gbc,0,n,2,1);
        spring = new JTextField(8);
        spring.setForeground(art.scrf);
        spring.setBackground(art.scrn);
        fix(spring,grd,gbc,2,n++,1,1);

        // Minimisation control option

        fix(new JLabel("Minimise Option:",JLabel.LEFT),grd,gbc,0,n,2,1);
        optkey = new JComboBox<String>();
        optkey.setForeground(art.scrf);
        optkey.setBackground(art.scrn);
        optkey.addItem("Force");
        optkey.addItem("Energy");
        optkey.addItem("Position");
        fix(optkey,grd,gbc,2,n++,1,1);

        //Minimisation control tolerance

        fix(new JLabel("Minimise Tolerance:",JLabel.LEFT),grd,gbc,0,n,2,1);
        opttol = new JTextField(8);
        opttol.setForeground(art.scrf);
        opttol.setBackground(art.scrn);
        fix(opttol,grd,gbc,2,n++,1,1);

        fix(new JLabel(" ",JLabel.LEFT),grd,gbc,0,n++,3,1);
        fix(new JLabel("Controls common to BPD & TAD only:",JLabel.LEFT),grd,gbc,0,n++,3,1);


        // Block interval between minimisations

        fix(new JLabel("Block Interval:",JLabel.LEFT),grd,gbc,0,n,2,1);
        block = new JTextField(8);
        block.setForeground(art.scrf);
        block.setBackground(art.scrn);
        fix(block,grd,gbc,2,n++,1,1);

        // Interval between tracking configurations

        fix(new JLabel("Tracking Interval:",JLabel.LEFT),grd,gbc,0,n,2,1);
        track = new JTextField(8);
        track.setForeground(art.scrf);
        track.setBackground(art.scrn);
        fix(track,grd,gbc,2,n++,1,1);

        // Catch radius for transitions

        fix(new JLabel("Catch Radius:",JLabel.LEFT),grd,gbc,0,n,2,1);
        tcatch = new JTextField(8);
        tcatch.setForeground(art.scrf);
        tcatch.setBackground(art.scrn);
        fix(tcatch,grd,gbc,2,n++,1,1);

        fix(new JLabel(" ",JLabel.LEFT),grd,gbc,0,n++,3,1);
        fix(new JLabel("Controls specific to BPD only:",JLabel.LEFT),grd,gbc,0,n++,3,1);

        // BPD path or just dynamics

        tpath=new JCheckBox("BPD Path required?");
        tpath.setForeground(art.fore);
        tpath.setBackground(art.back);
        fix(tpath,grd,gbc,0,n++,3,1);

        // Name of target atom type

        fix(new JLabel("Target Atom Name:",JLabel.LEFT),grd,gbc,0,n,2,1);
        target = new JTextField(8);
        target.setForeground(art.scrf);
        target.setBackground(art.scrn);
        fix(target,grd,gbc,2,n++,1,1);

        // Bias potential

        fix(new JLabel("Potential Bias:",JLabel.LEFT),grd,gbc,0,n,2,1);
        tebias = new JTextField(8);
        tebias.setForeground(art.scrf);
        tebias.setBackground(art.scrn);
        fix(tebias,grd,gbc,2,n++,1,1);

        // Bias potential minimum

        fix(new JLabel("Minimum Bias:",JLabel.LEFT),grd,gbc,0,n,2,1);
        tvmin = new JTextField(8);
        tvmin.setForeground(art.scrf);
        tvmin.setBackground(art.scrn);
        fix(tvmin,grd,gbc,2,n++,1,1);

        // NEB calculation selection

        tgoneb=new JCheckBox("NEB Calculation required?");
        tgoneb.setForeground(art.fore);
        tgoneb.setBackground(art.back);
        fix(tgoneb,grd,gbc,0,n++,3,1);

        fix(new JLabel(" ",JLabel.LEFT),grd,gbc,0,n++,3,1);
        fix(new JLabel("Controls specific to TAD only:",JLabel.LEFT),grd,gbc,0,n++,3,1);

        // Blackout period

        fix(new JLabel("Blackout Period:",JLabel.LEFT),grd,gbc,0,n,2,1);
        black = new JTextField(8);
        black.setForeground(art.scrf);
        black.setBackground(art.scrn);
        fix(black,grd,gbc,2,n++,1,1);

        // Deltad control parameter

        fix(new JLabel("Deltad Parameter:",JLabel.LEFT),grd,gbc,0,n,2,1);
        tdelt = new JTextField(8);
        tdelt.setForeground(art.scrf);
        tdelt.setBackground(art.scrn);
        fix(tdelt,grd,gbc,2,n++,1,1);

        // Low temperature target

        fix(new JLabel("Target Low Temp.:",JLabel.LEFT),grd,gbc,0,n,2,1);
        ttlow = new JTextField(8);
        ttlow.setForeground(art.scrf);
        ttlow.setBackground(art.scrn);
        fix(ttlow,grd,gbc,2,n++,1,1);


        fix(new JLabel(" ",JLabel.LEFT),grd,gbc,0,n++,3,1);
        fix(new JLabel("Controls specific to NEB only:",JLabel.LEFT),grd,gbc,0,n++,3,1);

        // Required number of NEB calculations

        fix(new JLabel("Number of NEB Calcs.:",JLabel.LEFT),grd,gbc,0,n,2,1);
        nnebs = new JTextField(8);
        nnebs.setForeground(art.scrf);
        nnebs.setBackground(art.scrn);
        fix(nnebs,grd,gbc,2,n++,1,1);

        // NEB start basins

        fix(new JLabel("Basins 1:",JLabel.LEFT),grd,gbc,0,n,2,1);
        basin1 = new JTextField(8);
        basin1.setForeground(art.scrf);
        basin1.setBackground(art.scrn);
        fix(basin1,grd,gbc,2,n++,1,1);

        // NEB end basins

        fix(new JLabel("Basins 2:",JLabel.LEFT),grd,gbc,0,n,2,1);
        basin2 = new JTextField(8);
        basin2.setForeground(art.scrf);
        basin2.setBackground(art.scrn);
        fix(basin2,grd,gbc,2,n++,1,1);

        fix(new JLabel(" ",JLabel.LEFT),grd,gbc,0,n++,3,1);

        // Define the Close button

        close = new JButton("Close");
        close.setBackground(art.butn);
        close.setForeground(art.butf);
        fix(close,grd,gbc,0,n++,1,1);

        // Register action buttons

        close.addActionListener(this);

    }

    public HyperDyn(MakeControl here) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        home=here;
        job=new HyperDyn();
        job.pack();
        job.setVisible(true);
        setParams();
    }

    void setParams() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */

        // set panel

        hypopt.setSelectedIndex(hyp_key);
        tunits.setSelectedIndex(hyp_units_key);
        spring.setText(String.valueOf(neb_spring));
        optkey.setSelectedIndex(hyp_opt_key);
        opttol.setText(String.valueOf(hyp_opt_tol));
        block.setText(String.valueOf(num_block));
        track.setText(String.valueOf(num_track));
        tcatch.setText(String.valueOf(catch_radius));
        tpath.setSelected(bpdpath);
        target.setText(String.valueOf(hyp_target));
        tebias.setText(String.valueOf(ebias));
        tvmin.setText(String.valueOf(vmin));
        tgoneb.setSelected(goneb);
        black.setText(String.valueOf(num_black));
        tdelt.setText(String.valueOf(deltad));
        ttlow.setText(String.valueOf(low_temp));
        nnebs.setText(String.valueOf(num_neb));
        basin1.setText(basin_1);
        basin2.setText(basin_2);

    }

    void getParams(){
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */

        lhyp=true;
        hyp_key=hypopt.getSelectedIndex();
        hyp_units_key=tunits.getSelectedIndex();
        neb_spring=BML.giveDouble(spring.getText(),1);
        hyp_opt_key=optkey.getSelectedIndex();
        hyp_opt_tol=BML.giveDouble(opttol.getText(),1);
        num_block=BML.giveInteger(block.getText(),1);
        num_track=BML.giveInteger(track.getText(),1);
        catch_radius=BML.giveDouble(tcatch.getText(),1);
        bpdpath=tpath.isSelected();
        hyp_target=target.getText();
        ebias=BML.giveDouble(tebias.getText(),1);
        vmin=BML.giveDouble(tvmin.getText(),1);
        goneb=tgoneb.isSelected();
        num_black=BML.giveInteger(black.getText(),1);
        deltad=BML.giveDouble(tdelt.getText(),1);
        low_temp=BML.giveDouble(ttlow.getText(),1);
        num_neb=BML.giveInteger(nnebs.getText(),1);
        basin_1=basin1.getText();
        basin_2=basin2.getText();

    }

    public void actionPerformed(ActionEvent e) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */

        String arg = (String)e.getActionCommand();
        if (arg.equals("Close")) {
            byebye();
        }
    }

    void byebye() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        getParams();
        home.hyp=null;
        job.setVisible(false);
    }
}
