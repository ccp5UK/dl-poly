import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;

public class RadDam extends Basic implements ActionListener {
        /*
*********************************************************************

dl_poly/java Radiation Damage class

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
    private static MakeControl home;
    private static RadDam job;
    private static JTextField tdstart,tdintval,tdefcut,tatom,timpstp,tenergy,tvect1;
    private static JTextField tvect2,tvect3,tvarstp,tmindis,tmaxdis,tnstmsdtmp;
    private static JTextField timsdtmp,tthick,tptemp;
    private static JComboBox pseudtyp;
    private static JButton close;
    private static JCheckBox bldefects,blvarstp,blmsdtmp,blpseudo;


    // Define the Graphical User Interface

    public RadDam() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        super();
        setTitle("Radiation Damage");
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

        // Impact description

        fix(new JLabel("Impact description:",JLabel.LEFT),grd,gbc,0,n++,3,1);
        tatom = new JTextField(8);
        tatom.setForeground(art.scrf);
        tatom.setBackground(art.scrn);
        fix(tatom,grd,gbc,0,n,1,1);
        fix(new JLabel("Atom index",JLabel.LEFT),grd,gbc,1,n++,2,1);
        timpstp = new JTextField(8);
        timpstp.setForeground(art.scrf);
        timpstp.setBackground(art.scrn);
        fix(timpstp,grd,gbc,0,n,1,1);
        fix(new JLabel("Impact time step no.",JLabel.LEFT),grd,gbc,1,n++,2,1);
        tenergy = new JTextField(8);
        tenergy.setForeground(art.scrf);
        tenergy.setBackground(art.scrn);
        fix(tenergy,grd,gbc,0,n,1,1);
        fix(new JLabel("Impact energy (keV)",JLabel.LEFT),grd,gbc,1,n++,2,1);
        fix(new JLabel("Impact vector:",JLabel.LEFT),grd,gbc,0,n++,3,1);
        tvect1 = new JTextField(6);
        tvect1.setForeground(art.scrf);
        tvect1.setBackground(art.scrn);
        fix(tvect1,grd,gbc,0,n,1,1);
        tvect2 = new JTextField(6);
        tvect2.setForeground(art.scrf);
        tvect2.setBackground(art.scrn);
        fix(tvect2,grd,gbc,1,n,1,1);
	tvect3 = new JTextField(6);
        tvect3.setForeground(art.scrf);
        tvect3.setBackground(art.scrn);
        fix(tvect3,grd,gbc,2,n++,1,1);
        fix(new JLabel(" "),grd,gbc,0,n++,1,1);

        // Variable time step option

        blvarstp = new JCheckBox("Use Variable time step");
        blvarstp.setForeground(art.fore);
        blvarstp.setBackground(art.back);
        fix(blvarstp,grd,gbc,0,n++,2,1);
        tvarstp = new JTextField(8);
        tvarstp.setBackground(art.scrn);
        tvarstp.setForeground(art.scrf);
        fix(tvarstp,grd,gbc,0,n,1,1);
        fix(new JLabel("Start time step value (ps)",JLabel.LEFT),grd,gbc,1,n++,2,1);
        tmindis = new JTextField(8);
        tmindis.setForeground(art.scrf);
        tmindis.setBackground(art.scrn);
        fix(tmindis,grd,gbc,0,n,1,1);
        fix(new JLabel("Min. allowed distance (A)",JLabel.LEFT),grd,gbc,1,n++,2,1);
        tmaxdis = new JTextField(8);
        tmaxdis.setForeground(art.scrf);
        tmaxdis.setBackground(art.scrn);
        fix(tmaxdis,grd,gbc,0,n,1,1);
        fix(new JLabel("Max. allowed distance (A)",JLabel.LEFT),grd,gbc,1,n++,2,1);
        fix(new JLabel(" "),grd,gbc,0,n++,1,1);

        // Defects

        bldefects=new JCheckBox("Produce DEFECTS file");
        bldefects.setForeground(art.fore);
        bldefects.setBackground(art.back);
        fix(bldefects,grd,gbc,0,n++,2,1);
        tdstart = new JTextField(8);
        tdstart.setForeground(art.scrf);
        tdstart.setBackground(art.scrn);
        fix(tdstart,grd,gbc,0,n,1,1);
        fix(new JLabel("Start time step number",JLabel.LEFT),grd,gbc,1,n++,2,1);
        tdintval = new JTextField(8);
        tdintval.setForeground(art.scrf);
        tdintval.setBackground(art.scrn);
        fix(tdintval,grd,gbc,0,n,1,1);
        fix(new JLabel("Time step interval",JLabel.LEFT),grd,gbc,1,n++,2,1);
        tdefcut = new JTextField(8);
        tdefcut.setForeground(art.scrf);
        tdefcut.setBackground(art.scrn);
        fix(tdefcut,grd,gbc,0,n,1,1);
        fix(new JLabel("Cut off (A):",JLabel.LEFT),grd,gbc,1,n++,2,1);
        fix(new JLabel(" "),grd,gbc,0,n++,1,1);

        // MSDTMP file option

        blmsdtmp = new JCheckBox("Produce MSDTMP file");
        blmsdtmp.setForeground(art.fore);
        blmsdtmp.setBackground(art.back);
        fix(blmsdtmp,grd,gbc,0,n++,2,1);
        tnstmsdtmp = new JTextField(8);
        tnstmsdtmp.setBackground(art.scrn);
        tnstmsdtmp.setForeground(art.scrf);
        fix(tnstmsdtmp,grd,gbc,0,n,1,1);
        fix(new JLabel("Start time step number",JLabel.LEFT),grd,gbc,1,n++,2,1);
        timsdtmp = new JTextField(8);
        timsdtmp.setBackground(art.scrn);
        timsdtmp.setForeground(art.scrf);
        fix(timsdtmp,grd,gbc,0,n,1,1);
        fix(new JLabel("Time step interval",JLabel.LEFT),grd,gbc,1,n++,2,1);
        fix(new JLabel(" "),grd,gbc,0,n++,1,1);

        // Pseudo thermal bath option

        blpseudo = new JCheckBox("Pseudo thermal bath");
        blpseudo.setForeground(art.fore);
        blpseudo.setBackground(art.back);
        fix(blpseudo,grd,gbc,0,n++,2,1);
        pseudtyp = new JComboBox();
        pseudtyp.setBackground(art.scrn);
        pseudtyp.setForeground(art.scrf);
        pseudtyp.addItem(" ");
        pseudtyp.addItem("Langevin");
        pseudtyp.addItem("Direct");
        fix(pseudtyp,grd,gbc,0,n,1,1);
        fix(new JLabel("Thermal bath type",JLabel.LEFT),grd,gbc,1,n++,2,1);
        tthick = new JTextField(8);
        tthick.setBackground(art.scrn);
        tthick.setForeground(art.scrf);
        fix(tthick,grd,gbc,0,n,1,1);
        fix(new JLabel("Layer thickness (A)",JLabel.LEFT),grd,gbc,1,n++,2,1);
	tptemp = new JTextField(8);
        tptemp.setBackground(art.scrn);
        tptemp.setForeground(art.scrf);
        fix(tptemp,grd,gbc,0,n,1,1);
        fix(new JLabel("Thermostat temp. (K)",JLabel.LEFT),grd,gbc,1,n++,2,1);
        fix(new JLabel(" "),grd,gbc,0,n++,1,1);

        // Define the Close button

        close = new JButton("Close");
        close.setBackground(art.butn);
        close.setForeground(art.butf);
        fix(close,grd,gbc,0,n++,1,1);

        // Register action buttons

        close.addActionListener(this);

    }

    public RadDam(MakeControl here) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        home=here;
        job=new RadDam();
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

        bldefects.setSelected(ldefects);
        tdstart.setText(String.valueOf(dstart));
        tdintval.setText(String.valueOf(dintval));
        tdefcut.setText(String.valueOf(defcut));
        tatom.setText(String.valueOf(atom));
        timpstp.setText(String.valueOf(impstp));
        tenergy.setText(String.valueOf(energy));
        tvect1.setText(String.valueOf(vect1));
        tvect2.setText(String.valueOf(vect2));
        tvect3.setText(String.valueOf(vect3));
        blvarstp.setSelected(lvarstp);
        tvarstp.setText(String.valueOf(varstp));
        tmindis.setText(String.valueOf(mindis));
        tmaxdis.setText(String.valueOf(maxdis));
        blmsdtmp.setSelected(lmsdtmp);
        tnstmsdtmp.setText(String.valueOf(nstmsdtmp));
        timsdtmp.setText(String.valueOf(imsdtmp));
        blpseudo.setSelected(lpseudo);
        pseudtyp.setSelectedIndex(ipseudtyp);
        tthick.setText(String.valueOf(thick));
        tptemp.setText(String.valueOf(ptemp));

    }

    void getParams(){
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */

        ldefects=bldefects.isSelected();
        dstart=BML.giveInteger(tdstart.getText(),1);
        dintval=BML.giveInteger(tdintval.getText(),1);
        defcut=BML.giveDouble(tdefcut.getText(),1);
        atom=BML.giveInteger(tatom.getText(),1);
        impstp=BML.giveInteger(timpstp.getText(),1);
        energy=BML.giveDouble(tenergy.getText(),1);
        vect1=BML.giveDouble(tvect1.getText(),1);
        vect2=BML.giveDouble(tvect2.getText(),1);
        vect3=BML.giveDouble(tvect3.getText(),1);
        lvarstp=blvarstp.isSelected();
        varstp=BML.giveDouble(tvarstp.getText(),1);
        mindis=BML.giveDouble(tmindis.getText(),1);
        maxdis=BML.giveDouble(tmaxdis.getText(),1);
        lmsdtmp=blmsdtmp.isSelected();
        nstmsdtmp=BML.giveInteger(tnstmsdtmp.getText(),1);
        imsdtmp=BML.giveInteger(timsdtmp.getText(),1);
        lpseudo=blpseudo.isSelected();
        ipseudtyp=pseudtyp.getSelectedIndex();
        thick=BML.giveDouble(tthick.getText(),1);
        ptemp=BML.giveDouble(tptemp.getText(),1);

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
        home.rdm=null;
        job.setVisible(false);
    }
}
