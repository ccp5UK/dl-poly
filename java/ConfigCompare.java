import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;

// Define the Graphical User Interface

public class ConfigCompare extends Basic implements ActionListener {
        /*
*********************************************************************

dl_poly/java GUI class to take compare two CONFIG files

copyright - daresbury laboratory
author    - w.smith August 2011

*********************************************************************
         */
    public static ConfigCompare job;
    private static GUI home=null;
    private static ButtonGroup options;
    private static JButton compare,close;
    private static JCheckBox first,second;
    private static JTextField cut,file1,file2;
    private static String fname1,fname2;
    private static Config cfgone,cfgtwo;
    private static double cutoff;

    public ConfigCompare() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        setTitle("Cfg. Compare");

        getContentPane().setBackground(art.back);
        getContentPane().setForeground(art.fore);
        setDefaultCloseOperation(DISPOSE_ON_CLOSE);
        setFont(fontMain);
        GridBagLayout grd = new GridBagLayout();
        GridBagConstraints gbc = new GridBagConstraints();
        getContentPane().setLayout(grd);

        gbc.fill=GridBagConstraints.BOTH;

        int n=0;

        // Panel header

        fix(new JLabel("Compare CONFIG files:",JLabel.LEFT),grd,gbc,0,n++,3,1);

        // Define the Compare button

        compare = new JButton("Compare");
        compare.setBackground(art.butn);
        compare.setForeground(art.butf);
        fix(compare,grd,gbc,0,n,1,1);
        fix(new JLabel("      ",JLabel.LEFT),grd,gbc,1,n,1,1);

        // Define the Close button

        close = new JButton("  Close  ");
        close.setBackground(art.butn);
        close.setForeground(art.butf);
        fix(close,grd,gbc,2,n++,1,1);

        // Define the first configuration file

        fix(new JLabel("First configuration file:",JLabel.LEFT),grd,gbc,0,n++,3,1);
        file1 = new JTextField(12);
        file1.setBackground(art.scrn);
        file1.setForeground(art.scrf);
        fix(file1,grd,gbc,0,n++,3,1);

        // Define the second configuration file

        fix(new JLabel("Second configuration file:",JLabel.LEFT),grd,gbc,0,n++,3,1);
        file2 = new JTextField(10);
        file2.setBackground(art.scrn);
        file2.setForeground(art.scrf);
        fix(file2,grd,gbc,0,n++,3,1);

        // Choose which configuration to display

        fix(new JLabel("Show first or second configuration?",JLabel.LEFT),grd,gbc,0,n++,3,1);
        options=new ButtonGroup();
        first=new JCheckBox("Config 1",true);
        first.setForeground(art.fore);
        first.setBackground(art.back);
        options.add(first);
        fix(first,grd,gbc,0,n,1,1);
        second=new JCheckBox("Config 2",false);
        second.setForeground(art.fore);
        second.setBackground(art.back);
        options.add(second);
        fix(second,grd,gbc,2,n++,1,1);

        // Set displacement criterion

        cut = new JTextField(8);
        cut.setBackground(art.scrn);
        cut.setForeground(art.scrf);
        fix(cut,grd,gbc,0,n,1,1);
        fix(new JLabel("Max. displacement (A)  ",JLabel.LEFT),grd,gbc,1,n++,2,1);

        // Register action buttons

        compare.addActionListener(this);
        close.addActionListener(this);

    }

    // Constructor method

    public ConfigCompare(GUI here) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        println("Activated panel for Comparing CONFIG files");
        home=here;

        // Set up Graphical User interface

        job = new ConfigCompare();
        job.pack();

        // Set default values

        cutoff=1.0;
        cut.setText(String.valueOf(cutoff));
        file1.setText(new String("CONFIG"));
        file2.setText(new String("REVCON"));
        job.setVisible(true);

    }

    public void actionPerformed(ActionEvent e) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        String arg = (String)e.getActionCommand();
        if (arg.equals("Compare")) {
            cutoff=Double.parseDouble(cut.getText());
            if(first.isSelected()) {
                fname1=file1.getText();
                fname2=file2.getText();
            }
            else if(second.isSelected()) {
                fname2=file1.getText();
                fname1=file2.getText();
            }
            compare();
        }
        else if (arg.equals("  Close  ")) {
            job.dispose();
        }
    }
    void compare() {
        /*
**********************************************************************

dl_poly/java utility to compare two CONFIG files

copyright daresbury laboratory

author w. smith august 2011

**********************************************************************
         */
        boolean safe;
        double cut2,rsq;

        safe=true;
        cut2=cutoff*cutoff;

        // Read first configuration

        cfgone=new Config();
        if(!cfgone.rdCFG(fname1)) {
            println("Problem reading CONFIG file: "+fname1);
            cfgone=null;
            return;
        }

        // Read second configuration

        cfgtwo=new Config();
        if(!cfgtwo.rdCFG(fname2)) {
            println("Problem reading CONFIG file: "+fname2);
            cfgtwo=null;
            return;
        }

        // Check CONFIG files have same contents

        if(cfgone.natms != cfgtwo.natms)
            safe=false;
        else if(cfgone.pbc.imcon != cfgtwo.pbc.imcon)
            safe=false;
        else {
            for (int i=0;i<cfgone.natms;i++) {
                if(!cfgone.atoms[i].zsym.equals(cfgtwo.atoms[i].zsym))
                    safe=false;
            }
        }

        if(!safe) {
            println("Error - config files "+fname1+" and "+fname2+" not same system");
            return;
        }

        // Copy first configuration

        config=copyConfig(cfgone);

        // Calculate atomic displacement

        for (int i=0;i<cfgone.natms;i++) {
            cfgtwo.xyz[0][i]-=cfgone.xyz[0][i];
            cfgtwo.xyz[1][i]-=cfgone.xyz[1][i];
            cfgtwo.xyz[2][i]-=cfgone.xyz[2][i];
        }

        // Minimum image correction

        cfgtwo.pbc.images(cfgtwo.natms,cfgtwo.xyz);

        // Identify dispaced atoms

        for (int i=0;i<cfgtwo.natms;i++) {
            rsq=Math.pow(cfgtwo.xyz[0][i],2)+Math.pow(cfgtwo.xyz[1][i],2)+Math.pow(cfgtwo.xyz[2][i],2);
            if(rsq > cut2) {
                println("Atom "+BML.fmt(i+1,8)+" has moved "+BML.fmt(Math.sqrt(rsq),10)+" A");
            }
            else {
                config.atoms[i].dotify=true;
            }
        }

        // Backup copy of new config

        cfgsav=copyConfig(config);

        // Show displaced atoms

	if(!editor.isVisible())
	    editor.showEditor();
        editor.pane.restore();
    }
}
