import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;
public class Solvat extends Basic implements ActionListener {
        /*
*********************************************************************

dl_poly/java Temperature Accelerated Dynamics class

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
    private static MakeControl home;
    private static Solvat job;
    private static JTextField tstart,tintvl,tlambda,tnonlin;
    private static JTextField systema,systemb,tswtch;
    private static JCheckBox remass;
    private static JComboBox solkey,mixkey;
    private static JButton close;


    // Define the Graphical User Interface

    public Solvat() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        super();
        setTitle("Solvation Options");
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

        // Select solvation option

        fix(new JLabel("Choose solvation option:",JLabel.LEFT),grd,gbc,0,n,2,1);
        solkey = new JComboBox();
        solkey.setForeground(art.scrf);
        solkey.setBackground(art.scrn);
        solkey.addItem("Decompose");
        solkey.addItem("Free Energy");
        solkey.addItem("Excitation");
        solkey.addItem("Switch");
        fix(solkey,grd,gbc,2,n++,1,1);

	// Start of decomposition

        fix(new JLabel("Start timestep:",JLabel.LEFT),grd,gbc,0,n,2,1);
        tstart = new JTextField(8);
        tstart.setForeground(art.scrf);
        tstart.setBackground(art.scrn);
        fix(tstart,grd,gbc,2,n++,1,1);

        // Data sampling interval

        fix(new JLabel("Sampling Interval:",JLabel.LEFT),grd,gbc,0,n,2,1);
        tintvl = new JTextField(8);
        tintvl.setForeground(art.scrf);
        tintvl.setBackground(art.scrn);
        fix(tintvl,grd,gbc,2,n++,1,1);

        // Solvent relaxation switching interval

        fix(new JLabel("Switching Interval:",JLabel.LEFT),grd,gbc,0,n,2,1);
        tswtch = new JTextField(8);
        tswtch.setForeground(art.scrf);
        tswtch.setBackground(art.scrn);
        fix(tswtch,grd,gbc,2,n++,1,1);

        // Free energy mixing parameter

        fix(new JLabel("Free Energy Lambda:",JLabel.LEFT),grd,gbc,0,n,2,1);
        tlambda = new JTextField(8);
        tlambda.setForeground(art.scrf);
        tlambda.setBackground(art.scrn);
        fix(tlambda,grd,gbc,2,n++,1,1);

        // Mixing function selction

        fix(new JLabel("Free Energy Mixing Function:",JLabel.LEFT),grd,gbc,0,n,2,1);
        mixkey = new JComboBox();
        mixkey.setForeground(art.scrf);
        mixkey.setBackground(art.scrn);
        mixkey.addItem("Linear");
        mixkey.addItem("Nonlinear");
        mixkey.addItem("Trigonometric");
        mixkey.addItem("Error function");
        mixkey.addItem("Polynomial");
        mixkey.addItem("Spline kernel");
        fix(mixkey,grd,gbc,2,n++,1,1);

        // Exponent for nonlinear mixing

        fix(new JLabel("Free Eng. Nonlinear exponent:",JLabel.LEFT),grd,gbc,0,n,2,1);
        tnonlin = new JTextField(8);
        tnonlin.setForeground(art.scrf);
        tnonlin.setBackground(art.scrn);
        fix(tnonlin,grd,gbc,2,n++,1,1);

        // Reset mass option

        remass = new JCheckBox("Reset mass for Free Energy?");
        remass.setForeground(art.fore);
        remass.setBackground(art.back);
        fix(remass,grd,gbc,0,n++,3,1);

        // Define System A

        fix(new JLabel("Define System A:",JLabel.LEFT),grd,gbc,0,n,2,1);
        systema = new JTextField(8);
        systema.setForeground(art.scrf);
        systema.setBackground(art.scrn);
        fix(systema,grd,gbc,2,n++,1,1);

        // Define System B

        fix(new JLabel("Define System B:",JLabel.LEFT),grd,gbc,0,n,2,1);
        systemb = new JTextField(8);
        systemb.setForeground(art.scrf);
        systemb.setBackground(art.scrn);
        fix(systemb,grd,gbc,2,n++,1,1);

        // Define the Close button

        close = new JButton("Close");
        close.setBackground(art.butn);
        close.setForeground(art.butf);
        fix(close,grd,gbc,0,n++,1,1);

        // Register action buttons

        close.addActionListener(this);

    }

    public Solvat(MakeControl here) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        home=here;
        job=new Solvat();
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

        solkey.setSelectedIndex(sol_key);
        tstart.setText(String.valueOf(num_start));
        tintvl.setText(String.valueOf(num_intvl));
        tswtch.setText(String.valueOf(num_swtch));
        tlambda.setText(String.valueOf(lambda));
        mixkey.setSelectedIndex(mix_key-1);
        tnonlin.setText(String.valueOf(non_lin_exp));
        remass.setSelected(lremass);
        systema.setText(system_a);
        systemb.setText(system_b);

    }

    void getParams(){
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */

	lsol=true;
        sol_key=solkey.getSelectedIndex();
        num_start=BML.giveInteger(tstart.getText(),1);
        num_intvl=BML.giveInteger(tintvl.getText(),1);
        num_swtch=BML.giveInteger(tswtch.getText(),1);
        lambda=BML.giveDouble(tlambda.getText(),1);
        mix_key=mixkey.getSelectedIndex()+1;
        non_lin_exp=BML.giveInteger(tnonlin.getText(),1);
        lremass=remass.isSelected();
	system_a=systema.getText();
	system_b=systemb.getText();

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
        home.sol=null;
        job.setVisible(false);
    }
}
