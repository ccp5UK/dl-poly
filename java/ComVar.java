import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;

public class ComVar extends Basic implements ActionListener {
        /*
*********************************************************************

dl_poly/java Common Variables class

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
    private static MakeControl home;
    private static ComVar job;
    private static JTextField tfcap,tistrdf,tnstraj,tistraj,tlevcon,ttscal;
    private static JTextField tistzden,tnstbts;
    private static JButton close;
    private static JCheckBox ballpairs,blcap,blzeql,blvdw,blrdf,blprdf;
    private static JCheckBox bltraj,bltscal,blzden,blpzden,blzero;

    // Define the Graphical User Interface

    public ComVar() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        super();
        setTitle("Common Variables");
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

        // All pairs

        ballpairs = new JCheckBox("All pairs");
        ballpairs.setForeground(art.fore);
        ballpairs.setBackground(art.back);
        fix(ballpairs,grd,gbc,0,n++,1,1);

        // Cap forces

        blcap = new JCheckBox("Cap Forces");
        blcap.setForeground(art.fore);
        blcap.setBackground(art.back);
        fix(blcap,grd,gbc,0,n++,1,1);

        // Force cap

        tfcap = new JTextField(8);
        tfcap.setBackground(art.scrn);
        tfcap.setForeground(art.scrf);
        fix(tfcap,grd,gbc,0,n,1,1);
        fix(new JLabel("Force cap",JLabel.LEFT),grd,gbc,1,n++,2,1);

        // Collect stats during equilibration

        blzeql = new JCheckBox("Collect");
        blzeql.setForeground(art.fore);
        blzeql.setBackground(art.back);
        fix(blzeql,grd,gbc,0,n++,1,1);

        // Disable VDW forces

        blvdw = new JCheckBox("Disable VDW forces");
        blvdw.setForeground(art.fore);
        blvdw.setBackground(art.back);
        fix(blvdw,grd,gbc,0,n++,2,1);

        // Calculate RDF and optional print

        blrdf = new JCheckBox("Calculate RDF");
        blrdf.setForeground(art.fore);
        blrdf.setBackground(art.back);
        fix(blrdf,grd,gbc,0,n,1,1);
        blprdf = new JCheckBox("Print RDF");
        blprdf.setForeground(art.fore);
        blprdf.setBackground(art.back);
        fix(blprdf,grd,gbc,2,n++,1,1);

        // RDF interval

        tistrdf = new JTextField(8);
        tistrdf.setBackground(art.scrn);
        tistrdf.setForeground(art.scrf);
        fix(tistrdf,grd,gbc,0,n,1,1);
        fix(new JLabel("RDF interval",JLabel.LEFT),grd,gbc,1,n++,2,1);

        // Produce History file

        bltraj = new JCheckBox("Produce History file");
        bltraj.setForeground(art.fore);
        bltraj.setBackground(art.back);
        fix(bltraj,grd,gbc,0,n++,2,1);

        // History file start,increment and data level

        tnstraj = new JTextField(8);
        tnstraj.setBackground(art.scrn);
        tnstraj.setForeground(art.scrf);
        fix(tnstraj,grd,gbc,0,n,1,1);
        fix(new JLabel("Start time step",JLabel.LEFT),grd,gbc,1,n++,2,1);
        tistraj = new JTextField(8);
        tistraj.setBackground(art.scrn);
        tistraj.setForeground(art.scrf);
        fix(tistraj,grd,gbc,0,n,1,1);
        fix(new JLabel("Time step  interval",JLabel.LEFT),grd,gbc,1,n++,2,1);
        tlevcon = new JTextField(8);
        tlevcon.setBackground(art.scrn);
        tlevcon.setForeground(art.scrf);
        fix(tlevcon,grd,gbc,0,n,1,1);
        fix(new JLabel("Data level key",JLabel.LEFT),grd,gbc,1,n++,2,1);

        // Enable Temp scaling

        bltscal = new JCheckBox("Enable T scaling");
        bltscal.setForeground(art.fore);
        bltscal.setBackground(art.back);
        fix(bltscal,grd,gbc,0,n++,2,1);

        // Temp scaling interval

        tnstbts = new JTextField(8);
        tnstbts.setBackground(art.scrn);
        tnstbts.setForeground(art.scrf);
        fix(tnstbts,grd,gbc,0,n,1,1);
        fix(new JLabel("T scaling interval",JLabel.LEFT),grd,gbc,1,n++,2,1);

        // Z density calculation

        blzden = new JCheckBox("Z density");
        blzden.setForeground(art.fore);
        blzden.setBackground(art.back);
        fix(blzden,grd,gbc,0,n,1,1);
        blpzden = new JCheckBox("Print Z-dens");
        blpzden.setForeground(art.fore);
        blpzden.setBackground(art.back);
        fix(blpzden,grd,gbc,2,n++,1,1);


        // Z density interval

        tistzden = new JTextField(8);
        tistzden.setBackground(art.scrn);
        tistzden.setForeground(art.scrf);
        fix(tistzden,grd,gbc,0,n,1,1);
        fix(new JLabel("Zdens interval",JLabel.LEFT),grd,gbc,1,n++,2,1);

        // Zero K MD option

        blzero = new JCheckBox("Zero K MD");
        blzero.setForeground(art.fore);
        blzero.setBackground(art.back);
        fix(blzero,grd,gbc,0,n++,1,1);

        // Define the Close button

        close = new JButton("Close");
        close.setBackground(art.butn);
        close.setForeground(art.butf);
        fix(close,grd,gbc,0,n++,1,1);

        // Register action buttons

        close.addActionListener(this);

    }

    public ComVar(MakeControl here) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        home=here;
        job=new ComVar();
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

        ballpairs.setSelected(allpairs);
        blcap.setSelected(lcap);
        tfcap.setText(String.valueOf(fcap));
        blzeql.setSelected(lzeql);
        blvdw.setSelected(lvdw);
        blprdf.setSelected(lprdf);
        blrdf.setSelected(lrdf);
        tistrdf.setText(String.valueOf(istrdf));
        bltraj.setSelected(ltraj);
        tnstraj.setText(String.valueOf(nstraj));
        tistraj.setText(String.valueOf(istraj));
        tlevcon.setText(String.valueOf(levcon));
        bltscal.setSelected(ltscal);
        tnstbts.setText(String.valueOf(nstbts));
        blzden.setSelected(lzden);
        blpzden.setSelected(lpzden);
        tistzden.setText(String.valueOf(istzden));
        blzero.setSelected(lzero);

    }

    void getParams() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */

        istrdf=BML.giveInteger(tistrdf.getText(),1);
        istzden=BML.giveInteger(tistzden.getText(),1);
        nstbts=BML.giveInteger(tnstbts.getText(),1);
        nstraj=BML.giveInteger(tnstraj.getText(),1);
        istraj=BML.giveInteger(tistraj.getText(),1);
        levcon=BML.giveInteger(tlevcon.getText(),1);
        fcap=BML.giveDouble(tfcap.getText(),1);
        allpairs=ballpairs.isSelected();
        lcap=blcap.isSelected();
        lvdw=blvdw.isSelected();
        lzeql=blzeql.isSelected();
        lrdf=blrdf.isSelected();
        lprdf=blprdf.isSelected();
        ltraj=bltraj.isSelected();
        ltscal=bltscal.isSelected();
        lzden=blzden.isSelected();
        lpzden=blpzden.isSelected();
        lzero=blzero.isSelected();

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
        home.cmv=null;
        job.setVisible(false);
    }
}
