import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;

public class MetaDyn extends Basic implements ActionListener {
        /*
*********************************************************************

dl_poly/java Bias Potential Dynamics class

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
    private static MakeControl home;
    private static MetaDyn job;
    private static JTextField tncolv,tglbscl,tlocscl,tnq4,tnq6;
    private static JTextField tntet,metstp,gheight,gwidth,ghkey;
    private static JTextField twtdt;
    private static JCheckBox tstein,ttet,tglob,tloc;
    private static JButton close;


    // Define the Graphical User Interface

    public MetaDyn() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        super();
        setTitle("Metadynamics");
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

        // Number of order parameters

        fix(new JLabel("Number of order params:",JLabel.LEFT),grd,gbc,0,n,2,1);
        tncolv = new JTextField(8);
        tncolv.setForeground(art.scrf);
        tncolv.setBackground(art.scrn);
        fix(tncolv,grd,gbc,2,n++,1,1);

        // Steinhardt parameters selection

        tstein=new JCheckBox("Steinhardt parameters required?");
        tstein.setForeground(art.fore);
        tstein.setBackground(art.back);
        fix(tstein,grd,gbc,0,n++,3,1);

        // Number of Steinhardt Q4 parameters required

        fix(new JLabel("No. of Steinhardt Q4 params:",JLabel.LEFT),grd,gbc,0,n,2,1);
        tnq4 = new JTextField(8);
        tnq4.setForeground(art.scrf);
        tnq4.setBackground(art.scrn);
        fix(tnq4,grd,gbc,2,n++,1,1);

        // Number of Steinhardt Q6 parameters required

        fix(new JLabel("No. of Steinhardt Q6 params:",JLabel.LEFT),grd,gbc,0,n,2,1);
        tnq6 = new JTextField(8);
        tnq6.setForeground(art.scrf);
        tnq6.setBackground(art.scrn);
        fix(tnq6,grd,gbc,2,n++,1,1);

       // Tetrahedral parameters selection

        ttet=new JCheckBox("Tetrahedral parameters required?");
        ttet.setForeground(art.fore);
        ttet.setBackground(art.back);
        fix(ttet,grd,gbc,0,n++,3,1);

        // Number of Tetrahedral parameters required

        fix(new JLabel("No. of tetrahedral params:",JLabel.LEFT),grd,gbc,0,n,2,1);
        tntet = new JTextField(8);
        tntet.setForeground(art.scrf);
        tntet.setBackground(art.scrn);
        fix(tntet,grd,gbc,2,n++,1,1);

        // Global potential energy parameters selection

        tglob=new JCheckBox("Gobal PE parameters required?");
        tglob.setForeground(art.fore);
        tglob.setBackground(art.back);
        fix(tglob,grd,gbc,0,n++,3,1);

        // Global potential energy scale factor

        fix(new JLabel("Global PE scale factor:",JLabel.LEFT),grd,gbc,0,n,2,1);
        tglbscl = new JTextField(8);
        tglbscl.setForeground(art.scrf);
        tglbscl.setBackground(art.scrn);
        fix(tglbscl,grd,gbc,2,n++,1,1);

        // Local potential energy parameters selection

        tloc=new JCheckBox("Local PE parameters required?");
        tloc.setForeground(art.fore);
        tloc.setBackground(art.back);
        fix(tloc,grd,gbc,0,n++,3,1);

        // Local potential energy scale factor

        fix(new JLabel("Local PE scale factor:",JLabel.LEFT),grd,gbc,0,n,2,1);
        tlocscl = new JTextField(8);
        tlocscl.setForeground(art.scrf);
        tlocscl.setBackground(art.scrn);
        fix(tlocscl,grd,gbc,2,n++,1,1);

        fix(new JLabel(" ",JLabel.LEFT),grd,gbc,0,n++,1,1);

        // Gaussian potential deposition interval

        fix(new JLabel("Gaussian deposition interval:",JLabel.LEFT),grd,gbc,0,n,2,1);
        metstp = new JTextField(8);
        metstp.setForeground(art.scrf);
        metstp.setBackground(art.scrn);
        fix(metstp,grd,gbc,2,n++,1,1);

        // Height of Gaussian potentials

        fix(new JLabel("Gaussian height:",JLabel.LEFT),grd,gbc,0,n,2,1);
        gheight = new JTextField(8);
        gheight.setForeground(art.scrf);
        gheight.setBackground(art.scrn);
        fix(gheight,grd,gbc,2,n++,1,1);

        // Width of Gaussian potentials

        fix(new JLabel("Gaussian width",JLabel.LEFT),grd,gbc,0,n,2,1);
        gwidth = new JTextField(8);
        gwidth.setForeground(art.scrf);
        gwidth.setBackground(art.scrn);
        fix(gwidth,grd,gbc,2,n++,1,1);

        //Gaussian control parameter

        fix(new JLabel("Gaussian control key:",JLabel.LEFT),grd,gbc,0,n,2,1);
        ghkey = new JTextField(8);
        ghkey.setForeground(art.scrf);
        ghkey.setBackground(art.scrn);
        fix(ghkey,grd,gbc,2,n++,1,1);

        // Well tempered dynamics control key

        fix(new JLabel("WTD control parameter",JLabel.LEFT),grd,gbc,0,n,2,1);
        twtdt = new JTextField();
        twtdt.setForeground(art.scrf);
        twtdt.setBackground(art.scrn);
        fix(twtdt,grd,gbc,2,n++,1,1);

        // Define the Close button

        close = new JButton("Close");
        close.setBackground(art.butn);
        close.setForeground(art.butf);
        fix(close,grd,gbc,0,n++,1,1);

        // Register action buttons

        close.addActionListener(this);

    }

    public MetaDyn(MakeControl here) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        home=here;
        job=new MetaDyn();
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

        tncolv.setText(String.valueOf(ncolvar));
        tstein.setSelected(lstein);
        tnq4.setText(String.valueOf(nq4));
        tnq6.setText(String.valueOf(nq6));
        ttet.setSelected(ltet);
        tntet.setText(String.valueOf(ntet));
        tglob.setSelected(lglobpe);
        tglbscl.setText(String.valueOf(globpe_scale));
        tloc.setSelected(llocpe);
        tlocscl.setText(String.valueOf(locpe_scale));
        metstp.setText(String.valueOf(meta_step_int));
        gheight.setText(String.valueOf(ref_w_aug));
        gwidth.setText(String.valueOf(h_aug));
        ghkey.setText(String.valueOf(hkey));
        twtdt.setText(String.valueOf(wt_dt));

    }

    void getParams(){
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */

        lmetd=true;
        ncolvar=BML.giveInteger(tncolv.getText(),1);
        lstein=tstein.isSelected();
        nq4=BML.giveInteger(tnq4.getText(),1);
        nq6=BML.giveInteger(tnq6.getText(),1);
        ltet=ttet.isSelected();
        ntet=BML.giveInteger(tntet.getText(),1);
        lglobpe=tglob.isSelected();
        globpe_scale=BML.giveDouble(tglbscl.getText(),1);
        llocpe=tloc.isSelected();
        locpe_scale=BML.giveDouble(tlocscl.getText(),1);
        meta_step_int=BML.giveInteger(metstp.getText(),1);
        ref_w_aug=BML.giveDouble(gheight.getText(),1);
        h_aug=BML.giveDouble(gwidth.getText(),1);
        hkey=BML.giveInteger(ghkey.getText(),1);
        wt_dt=BML.giveDouble(twtdt.getText(),1);

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
        home.metd=null;
        job.setVisible(false);
    }
}
