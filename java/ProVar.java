import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;

public class ProVar extends Basic implements ActionListener {
        /*
*********************************************************************

dl_poly/java Program Controls class

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
    private static MakeControl home;
    private static ProVar job;
    private static JComboBox restopt,algorithm;
    private static JTextField tnstrun,tnsteql,tmult,tnstbpo,tnstack,tintsta,tewltol,tshktol,tqtntol;
    private static JTextField tjobtim,ttclose;
    private static JButton close;

    // Define the Graphical User Interface

    public ProVar() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        super();
        setTitle("Program Variables");
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

        fix(new JLabel("Select options:",JLabel.LEFT),grd,gbc,0,n++,3,1);

        // Restart option

        fix(new JLabel("Algorithm",JLabel.LEFT),grd,gbc,0,n,2,1);
        algorithm = new JComboBox();
        algorithm.setBackground(art.scrn);
        algorithm.setForeground(art.scrf);
        algorithm.addItem("Leapfrog Verlet");
        algorithm.addItem("Velocity Verlet");
        fix(algorithm,grd,gbc,2,n++,1,1);

        // Number of steps

        tnstrun = new JTextField(8);
        tnstrun.setBackground(art.scrn);
        tnstrun.setForeground(art.scrf);
        fix(tnstrun,grd,gbc,0,n,1,1);
        fix(new JLabel("Number of steps",JLabel.LEFT),grd,gbc,1,n++,2,1);

        // Equilibration steps

        tnsteql = new JTextField(8);
        tnsteql.setBackground(art.scrn);
        tnsteql.setForeground(art.scrf);
        fix(tnsteql,grd,gbc,0,n,1,1);
        fix(new JLabel("Equilibration steps",JLabel.LEFT),grd,gbc,1,n++,2,1);

        // Multiple time step

        tmult = new JTextField(8);
        tmult.setBackground(art.scrn);
        tmult.setForeground(art.scrf);
        fix(tmult,grd,gbc,0,n,1,1);
        fix(new JLabel("Multiple time step",JLabel.LEFT),grd,gbc,1,n++,2,1);

        // Print interval

        tnstbpo = new JTextField(8);
        tnstbpo.setBackground(art.scrn);
        tnstbpo.setForeground(art.scrf);
        fix(tnstbpo,grd,gbc,0,n,1,1);
        fix(new JLabel("Print interval",JLabel.LEFT),grd,gbc,1,n++,2,1);

        // Stack interval

        tnstack = new JTextField(8);
        tnstack.setBackground(art.scrn);
        tnstack.setForeground(art.scrf);
        fix(tnstack,grd,gbc,0,n,1,1);
        fix(new JLabel("Stack interval",JLabel.LEFT),grd,gbc,1,n++,2,1);

        // Stats interval

        tintsta = new JTextField(8);
        tintsta.setBackground(art.scrn);
        tintsta.setForeground(art.scrf);
        fix(tintsta,grd,gbc,0,n,1,1);
        fix(new JLabel("Stats interval",JLabel.LEFT),grd,gbc,1,n++,2,1);

        // Ewald precision

        tewltol = new JTextField(8);
        tewltol.setBackground(art.scrn);
        tewltol.setForeground(art.scrf);
        fix(tewltol,grd,gbc,0,n,1,1);
        fix(new JLabel("Ewald precision",JLabel.LEFT),grd,gbc,1,n++,2,1);

        // SHAKE tolerance

        tshktol = new JTextField(8);
        tshktol.setBackground(art.scrn);
        tshktol.setForeground(art.scrf);
        fix(tshktol,grd,gbc,0,n,1,1);
        fix(new JLabel("SHAKE tolerance",JLabel.LEFT),grd,gbc,1,n++,2,1);

        // Quaternion tolerance

        tqtntol = new JTextField(8);
        tqtntol.setBackground(art.scrn);
        tqtntol.setForeground(art.scrf);
        fix(tqtntol,grd,gbc,0,n,1,1);
        fix(new JLabel("Quaternion tolerance",JLabel.LEFT),grd,gbc,1,n++,2,1);

        // Job time

        tjobtim = new JTextField(8);
        tjobtim.setBackground(art.scrn);
        tjobtim.setForeground(art.scrf);
        fix(tjobtim,grd,gbc,0,n,1,1);
        fix(new JLabel("Job time (s)",JLabel.LEFT),grd,gbc,1,n++,2,1);

        // Close time

        ttclose = new JTextField(8);
        ttclose.setBackground(art.scrn);
        ttclose.setForeground(art.scrf);
        fix(ttclose,grd,gbc,0,n,1,1);
        fix(new JLabel("Close time (s)",JLabel.LEFT),grd,gbc,1,n++,2,1);

        // Restart option

        fix(new JLabel("Restart option",JLabel.LEFT),grd,gbc,0,n,2,1);
        restopt = new JComboBox();
        restopt.setBackground(art.scrn);
        restopt.setForeground(art.scrf);
        restopt.addItem("NONE");
        restopt.addItem("RESTART");
        restopt.addItem("RESTART SCALE");
        restopt.addItem("RESTART NOSCALE");
        fix(restopt,grd,gbc,2,n++,1,1);

        // Define the Close button

        close = new JButton("Close");
        close.setBackground(art.butn);
        close.setForeground(art.butf);
        fix(close,grd,gbc,0,n++,1,1);

        // Register action buttons

        close.addActionListener(this);

    }

    public ProVar(MakeControl here) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        home=here;
        job=new ProVar();
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

        // set panel contents

        algorithm.setSelectedIndex(keyalg);
        tmult.setText(String.valueOf(mult));
        tnstrun.setText(String.valueOf(nstrun));
        tnsteql.setText(String.valueOf(nsteql));
        tnstbpo.setText(String.valueOf(nstbpo));
        tnstack.setText(String.valueOf(nstack));
        tintsta.setText(String.valueOf(intsta));
        tewltol.setText(String.valueOf(ewltol));
        tshktol.setText(String.valueOf(shktol));
        tqtntol.setText(String.valueOf(qtntol));
        tjobtim.setText(String.valueOf(jobtim));
        ttclose.setText(String.valueOf(tclose));
        restopt.setSelectedIndex(keyres);
    }

    void getParams(){
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */

        keyalg=algorithm.getSelectedIndex();
        ewltol=BML.giveDouble(tewltol.getText(),1);
        qtntol=BML.giveDouble(tqtntol.getText(),1);
        shktol=BML.giveDouble(tshktol.getText(),1);
        jobtim=BML.giveDouble(tjobtim.getText(),1);
        tclose=BML.giveDouble(ttclose.getText(),1);
        nstrun=BML.giveInteger(tnstrun.getText(),1);
        nsteql=BML.giveInteger(tnsteql.getText(),1);
        mult=BML.giveInteger(tmult.getText(),1);
        nstbpo=BML.giveInteger(tnstbpo.getText(),1);
        nstack=BML.giveInteger(tnstack.getText(),1);
        intsta=BML.giveInteger(tintsta.getText(),1);
        keyres=restopt.getSelectedIndex();

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
        home.prv=null;
        job.setVisible(false);
    }

}
