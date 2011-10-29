import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;

// Define the Graphical User Interface

public class Slice extends Basic implements ActionListener {
        /*
*********************************************************************

dl_poly/java GUI class to take a slice from a CONFIG/REVCON file

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
    public static Slice job;
    private static GUI home=null;
    private static double slx,sly,slz,top,bot;
    private static JButton load,make,close;
    private static JTextField dlx,dly,dlz,ubd,lbd;

    public Slice() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        setTitle("Slice a CONFIG file");

        getContentPane().setBackground(art.back);
        getContentPane().setForeground(art.fore);
        setDefaultCloseOperation(DISPOSE_ON_CLOSE);
        setFont(fontMain);
        GridBagLayout grd = new GridBagLayout();
        GridBagConstraints gbc = new GridBagConstraints();
        getContentPane().setLayout(grd);

        gbc.fill=GridBagConstraints.BOTH;

        // Define the Load button

        load = new JButton("Load");
        load.setBackground(art.butn);
        load.setForeground(art.butf);
        fix(load,grd,gbc,0,0,1,1);

        // Define the Make button

        make = new JButton("Make");
        make.setBackground(art.butn);
        make.setForeground(art.butf);
        fix(make,grd,gbc,2,0,1,1);

        // Slice direction vector

        JLabel lab1 = new JLabel("Slice direction vector:",JLabel.LEFT);
        fix(lab1,grd,gbc,0,1,3,1);
        dlx = new JTextField(8);
        dlx.setBackground(art.scrn);
        dlx.setForeground(art.scrf);
        fix(dlx,grd,gbc,0,2,1,1);

        dly = new JTextField(8);
        dly.setBackground(art.scrn);
        dly.setForeground(art.scrf);
        fix(dly,grd,gbc,1,2,1,1);

        dlz = new JTextField(8);
        dlz.setBackground(art.scrn);
        dlz.setForeground(art.scrf);
        fix(dlz,grd,gbc,2,2,1,1);
        fix(new JLabel("        "),grd,gbc,1,3,1,1);

        // Upper bound of slice

        JLabel lab2 = new JLabel("Upper bound:",JLabel.RIGHT);
        fix(lab2,grd,gbc,0,4,2,1);
        ubd = new JTextField(8);
        ubd.setBackground(art.scrn);
        ubd.setForeground(art.scrf);
        fix(ubd,grd,gbc,2,4,1,1);

        // Lower bound of slice

        JLabel lab3 = new JLabel("Lower bound:",JLabel.RIGHT);
        fix(lab3,grd,gbc,0,5,2,1);
        lbd = new JTextField(8);
        lbd.setBackground(art.scrn);
        lbd.setForeground(art.scrf);
        fix(lbd,grd,gbc,2,5,1,1);

        // Define the Close button

        close = new JButton("Close");
        close.setBackground(art.butn);
        close.setForeground(art.butf);
        fix(close,grd,gbc,0,6,1,1);

        // Register action buttons

        load.addActionListener(this);
        make.addActionListener(this);
        close.addActionListener(this);

    }

    // Constructor method

    public Slice(GUI here) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        println("Activated panel for slicing a CONFIG file");
        home=here;

        // Set up Graphical User interface

        job = new Slice();
        job.pack();
        job.setVisible(true);
        setValues();
    }

    // Set default values

    static void setValues() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        slx=0.0;
        sly=0.0;
        slz=1.0;
        top=3.0;
        bot=-3.0;
        dlx.setText(String.valueOf(slx));
        dly.setText(String.valueOf(sly));
        dlz.setText(String.valueOf(slz));
        ubd.setText(String.valueOf(top));
        lbd.setText(String.valueOf(bot));
    }

    public void actionPerformed(ActionEvent e) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        String arg = (String)e.getActionCommand();
        if (arg.equals("Load")) {
            slx=BML.giveDouble(dlx.getText(),1);
            sly=BML.giveDouble(dly.getText(),1);
            slz=BML.giveDouble(dlz.getText(),1);
            top=BML.giveDouble(ubd.getText(),1);
            bot=BML.giveDouble(lbd.getText(),1);
            config=getConfig(home,"CFG");
            sliceFile();
        }
        if (arg.equals("Make")) {
            slx=BML.giveDouble(dlx.getText(),1);
            sly=BML.giveDouble(dly.getText(),1);
            slz=BML.giveDouble(dlz.getText(),1);
            top=BML.giveDouble(ubd.getText(),1);
            bot=BML.giveDouble(lbd.getText(),1);
            sliceFile();
        }
        else if (arg.equals("Close")) {
            job.dispose();
        }
    }
    void sliceFile() {
        /*
**********************************************************************

dl_poly/java utility to cut a slice from a CONFIG file

copyright daresbury laboratory

author w. smith march 2001

**********************************************************************
         */
        int n;
        double ddd,sss;

        if(config==null || config.natms==0)return;

        n=0;
        sss=Math.sqrt(slx*slx+sly*sly+slz*slz);
        for(int i=0;i<config.natms;i++) {
            ddd=(slx*config.xyz[0][i]+sly*config.xyz[1][i]+slz*config.xyz[2][i])/sss;
            if(ddd > bot && ddd < top) {
                config.atoms[n].znum=config.atoms[i].znum;
                config.atoms[n].zmas=config.atoms[i].zmas;
                config.atoms[n].zchg=config.atoms[i].zchg;
                config.atoms[n].zrad=config.atoms[i].zrad;
                config.atoms[n].zsym=new String(config.atoms[i].zsym);
                config.atoms[n].zcol=new Color(config.atoms[i].zcol.getRGB());
                config.atoms[n].covalent=config.atoms[i].covalent;
                config.atoms[n].dotify=config.atoms[i].dotify;
                config.xyz[0][n]=config.xyz[0][i];
                config.xyz[1][n]=config.xyz[1][i];
                config.xyz[2][n]=config.xyz[2][i];
                n++;
            }
        }
        config.natms=n;

        // write new CONFIG file

        fname="CFGSLC."+String.valueOf(numslc);
        if(config.configWrite(fname)){
            println("File "+fname+" created");
            println("Number of atoms in "+fname+" : "+config.natms);
            numslc++;
        }

        // Draw sliced structure

        config.structure=new Structure(config);
        if(!editor.isVisible())
            editor.showEditor();
        editor.pane.restore();
    }
}
