import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;

public class SysVar extends Basic implements ActionListener {
        /*
*********************************************************************

dl_poly/java System Variables class

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
    private static MakeControl home=null;
    private static SysVar job=null;
    private static JComboBox ensemble,electro;
    private static JTextField ttemp,tpress,ttstep,trcut,tdelr,trvdw,trprim,tepsq,ttaut,ttaup,tgamt;
    private static JButton close;

    // Define the Graphical User Interface

    public SysVar() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        super();
        setTitle("System Variables");
        int n=0;

        getContentPane().setBackground(art.back);
        getContentPane().setForeground(art.fore);
        setDefaultCloseOperation(DO_NOTHING_ON_CLOSE);
        setFont(fontMain);
        GridBagLayout grd = new GridBagLayout();
        GridBagConstraints gbc = new GridBagConstraints();
        getContentPane().setLayout(grd);

        gbc.fill=GridBagConstraints.BOTH;

        // Panel label

        fix(new JLabel("Select options:",JLabel.LEFT),grd,gbc,0,n++,2,1);

        // Temperature

        ttemp = new JTextField(8);
        ttemp.setBackground(art.scrn);
        ttemp.setForeground(art.scrf);
        fix(ttemp,grd,gbc,0,n,1,1);
        fix(new JLabel("Temperature (K)",JLabel.LEFT),grd,gbc,1,n++,1,1);

        // Pressure

        tpress = new JTextField(8);
        tpress.setBackground(art.scrn);
        tpress.setForeground(art.scrf);
        fix(tpress,grd,gbc,0,n,1,1);
        fix(new JLabel("Pressure (kBars)",JLabel.LEFT),grd,gbc,1,n++,1,1);

        // Time step

        ttstep = new JTextField(8);
        ttstep.setBackground(art.scrn);
        ttstep.setForeground(art.scrf);
        fix(ttstep,grd,gbc,0,n,1,1);
        fix(new JLabel("Time step (ps)",JLabel.LEFT),grd,gbc,1,n++,1,1);

        // Forces cut off

        trcut = new JTextField(8);
        trcut.setBackground(art.scrn);
        trcut.setForeground(art.scrf);
        fix(trcut,grd,gbc,0,n,1,1);
        fix(new JLabel("Cut off (A)",JLabel.LEFT),grd,gbc,1,n++,1,1);

        // Verlet shell width

        tdelr = new JTextField(8);
        tdelr.setBackground(art.scrn);
        tdelr.setForeground(art.scrf);
        fix(tdelr,grd,gbc,0,n,1,1);
        fix(new JLabel("Verlet shell width (A)",JLabel.LEFT),grd,gbc,1,n++,1,1);

        // VDW cut off

        trvdw = new JTextField(8);
        trvdw.setBackground(art.scrn);
        trvdw.setForeground(art.scrf);
        fix(trvdw,grd,gbc,0,n,1,1);
        fix(new JLabel("VDW cut off (A)",JLabel.LEFT),grd,gbc,1,n++,1,1);

        // Primary cut off

        trprim = new JTextField(8);
        trprim.setBackground(art.scrn);
        trprim.setForeground(art.scrf);
        fix(trprim,grd,gbc,0,n,1,1);
        fix(new JLabel("Primary cut off (A)",JLabel.LEFT),grd,gbc,1,n++,1,1);

        // Dielectric Constant

        tepsq = new JTextField(8);
        tepsq.setBackground(art.scrn);
        tepsq.setForeground(art.scrf);
        fix(tepsq,grd,gbc,0,n,1,1);
        fix(new JLabel("Dielectric constant",JLabel.LEFT),grd,gbc,1,n++,1,1);

        // Temperature relaxation constant

        ttaut = new JTextField(8);
        ttaut.setBackground(art.scrn);
        ttaut.setForeground(art.scrf);
        fix(ttaut,grd,gbc,0,n,1,1);
        fix(new JLabel("Temp. relaxation (ps)",JLabel.LEFT),grd,gbc,1,n++,1,1);

        // Pressure relaxation constant

        ttaup = new JTextField(8);
        ttaup.setBackground(art.scrn);
        ttaup.setForeground(art.scrf);
        fix(ttaup,grd,gbc,0,n,1,1);
        fix(new JLabel("Press. relaxation (ps)",JLabel.LEFT),grd,gbc,1,n++,1,1);

        // Surface tension

        tgamt = new JTextField(8);
        tgamt.setBackground(art.scrn);
        tgamt.setForeground(art.scrf);
        fix(tgamt,grd,gbc,0,n,1,1);
        fix(new JLabel("Surface tension (dyn/cm)",JLabel.LEFT),grd,gbc,1,n++,1,1);

        // Choice of ensemble

        fix(new JLabel("Ensemble.............",JLabel.LEFT),grd,gbc,0,n,1,1);
        ensemble = new JComboBox();
        ensemble.setBackground(art.scrn);
        ensemble.setForeground(art.scrf);
        ensemble.addItem("NVE");
        ensemble.addItem("NVT Andersen");
        ensemble.addItem("NVT Berendsen");
        ensemble.addItem("NVT Evans");
        ensemble.addItem("NVT Hoover");
        ensemble.addItem("NVT Langevin");
        ensemble.addItem("NPT Berendsen");
        ensemble.addItem("NPT Hoover");
        ensemble.addItem("NPT Langevin");
        ensemble.addItem("NPT MTK");
        ensemble.addItem("NST Berendsen");
        ensemble.addItem("NST Hoover");
        ensemble.addItem("NST Langevin");
        ensemble.addItem("NST MTK");
        ensemble.addItem("NST-Area Berendsen");
        ensemble.addItem("NST-Area Hoover");
        ensemble.addItem("NST-Area Langevin");
        ensemble.addItem("NST-Area MTK");
        ensemble.addItem("NST-Tens Berendsen");
        ensemble.addItem("NST-Tens Hoover");
        ensemble.addItem("NST-Tens Langevin");
        ensemble.addItem("NST-Tens MTK");
        ensemble.addItem("PMF");
        fix(ensemble,grd,gbc,2,n++,1,1);

        // Choice of electrostatics

        fix(new JLabel("Electrostatics.......",JLabel.LEFT),grd,gbc,0,n,1,1);
        electro = new JComboBox();
        electro.setBackground(art.scrn);
        electro.setForeground(art.scrf);
        electro.addItem("NONE");
        electro.addItem("EWALD");
        electro.addItem("DISTAN");
        electro.addItem("T-COUL");
        electro.addItem("S-COUL");
        electro.addItem("R-FIELD");
        electro.addItem("SPME");
        electro.addItem("HK-EWALD");
        fix(electro,grd,gbc,2,n++,1,1);

        // Define the Close button

        close = new JButton("Close");
        close.setBackground(art.butn);
        close.setForeground(art.butf);
        fix(close,grd,gbc,0,n++,1,1);

        // Register action buttons

        close.addActionListener(this);

    }

    public SysVar(MakeControl here) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        home=here;
        job=new SysVar();
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

        ttemp.setText(String.valueOf(temp));
        tpress.setText(String.valueOf(press));
        ttstep.setText(String.valueOf(tstep));
        trcut.setText(String.valueOf(rcut));
        tdelr.setText(String.valueOf(delr));
        trvdw.setText(String.valueOf(rvdw));
        trprim.setText(String.valueOf(rprim));
        tepsq.setText(String.valueOf(epsq));
        ttaut.setText(String.valueOf(taut));
        ttaup.setText(String.valueOf(taup));
        tgamt.setText(String.valueOf(gamt));
        ensemble.setSelectedIndex(keyens);
        electro.setSelectedIndex(keyfce/2);
    }

    void getParams(){
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */

        temp=BML.giveDouble(ttemp.getText(),1);
        press=BML.giveDouble(tpress.getText(),1);
        tstep=BML.giveDouble(ttstep.getText(),1);
        rcut=BML.giveDouble(trcut.getText(),1);
        delr=BML.giveDouble(tdelr.getText(),1);
        rvdw=BML.giveDouble(trvdw.getText(),1);
        rprim=BML.giveDouble(trprim.getText(),1);
        epsq=BML.giveDouble(tepsq.getText(),1);
        taut=BML.giveDouble(ttaut.getText(),1);
        taup=BML.giveDouble(ttaup.getText(),1);
        keyens=ensemble.getSelectedIndex();
        keyfce=2*electro.getSelectedIndex();
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
        home.sys=null;
        job.setVisible(false);
    }

}

