import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;

// Define the Graphical User Interface

public class InsertMolecule extends Basic implements ActionListener {
    /*
*********************************************************************

dl_poly/java GUI class to insert a molecule into a CONFIG file

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
     */
    public static InsertMolecule job;
    private static JButton load,insert,close;
    private static JTextField dsm;
    private static GUI home=null;
    private static double dis;

    public InsertMolecule() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        setTitle("Insert Molecule Panel");

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
        fix(new JLabel("        "),grd,gbc,1,0,1,1);

        // Define the Insert button

        insert = new JButton("Insert");
        insert.setBackground(art.butn);
        insert.setForeground(art.butf);
        fix(insert,grd,gbc,2,0,1,1);

        // molecule_substrate minimum distance

        JLabel lab1 = new JLabel("Min molecule-substrate distance:",JLabel.LEFT);
        fix(lab1,grd,gbc,0,1,2,1);
        dsm = new JTextField(8);
        dsm.setBackground(art.scrn);
        dsm.setForeground(art.scrf);
        fix(dsm,grd,gbc,2,1,1,1);

        // Define the Close button

        close = new JButton("Close");
        close.setBackground(art.butn);
        close.setForeground(art.butf);
        fix(close,grd,gbc,0,8,1,1);

        // Register action buttons

        load.addActionListener(this);
        insert.addActionListener(this);
        close.addActionListener(this);

    }

    // Constructor method

    public InsertMolecule(GUI here) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        home=here;
        println("Activated panel for inserting molecule into a CONFIG file");

        // Set up Graphical User interface

        job = new InsertMolecule();
        job.pack();
        job.setVisible(true);
        dis=2.5;
        dsm.setText(String.valueOf(dis));

    }

    public void actionPerformed(ActionEvent e) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        int call;
        String arg = (String)e.getActionCommand();
        if (arg.equals("Load")){
            println("Please supply the substrate CONFIG file");
            config=getConfig(home,"CFG");
            if(!editor.isVisible())
                editor.showEditor();
	    editor.pane.restore();
        }
       else if (arg.equals("Insert")) {
           if(config == null || config.natms == 0){
               println("Please supply the substrate CONFIG file");
               config=getConfig(home,"CFG");
           }
           dis=BML.giveDouble(dsm.getText(),1);
           call=addmolecule();
           if(!editor.isVisible())
               editor.showEditor();
	    editor.pane.restore();
       }
       else if (arg.equals("Close")) {
           job.dispose();
       }
    }

    int addmolecule() {
        /*
***********************************************************************

dl_poly/java utility to add a molecule to a substrate with removal of
substrate molecules to accommodate the insert.
Assumes atomic positions are in a form compatible
with the CONFIG file used in DL_POLY with periodic boundary
conditions.

Substrate is current configuration. Insert is in an external config file

copyright - daresbury laboratory
author    - w.smith 2011

***********************************************************************
         */

        int k;
        boolean zapp;
        Config newcfg;
        double ddd,dis2;
        double dsp[]=new double[3];
        int hitlist[]=new int[config.natms];

        // Get molecule to be inserted

        println("Please supply a molecule for insertion");
        newcfg=getConfig(home,"CFG");

        dis2=dis*dis;

        for(int i=0;i<config.natms;i++)
            hitlist[i]=0;

        // identify overlapping atoms

        for(int j=0;j<newcfg.natms;j++) {

            if(newcfg.atoms[j].znum > 1) {

                for(int i=0;i<config.natms;i++) {

                    if(config.atoms[i].znum > 1) {

                        // atomic separation

                        dsp[0]=config.xyz[0][i]-newcfg.xyz[0][j];
                        dsp[1]=config.xyz[1][i]-newcfg.xyz[1][j];
                        dsp[2]=config.xyz[2][i]-newcfg.xyz[2][j];

                        // minimum image

                        config.pbc.images(dsp);

                        // interatomic distance

                        ddd=dsp[0]*dsp[0]+dsp[1]*dsp[1]+dsp[2]*dsp[2];

                        // overlap check

                        if(ddd < dis2) hitlist[i]=1;

                    }
                }
            }
        }

        // identify redundant molecules

        if(config.structure.nbnds > 0) {

            for(int n=0;n<config.structure.nmols;n++) {

                zapp=false;

                for(int m=0;m<config.structure.isz[n];m++) {

                    if(hitlist[config.structure.ist[n]+m] > 0) zapp=true;

                }

                if(zapp) {
                    for(int m=0;m<config.structure.isz[n];m++) {

                        hitlist[config.structure.ist[n]+m]=1;

                    }
                }
            }
        }

        // remove redundant atoms

        k=0;

        for(int i=0;i<config.natms;i++) {

            if(hitlist[i] < 1) {

                config.atoms[k].znum=config.atoms[i].znum;
                config.atoms[k].zmas=config.atoms[i].zmas;
                config.atoms[k].zchg=config.atoms[i].zchg;
                config.atoms[k].zrad=config.atoms[i].zrad;
                config.atoms[k].zsym=new String(config.atoms[i].zsym);
                config.atoms[k].zcol=new Color(config.atoms[i].zcol.getRGB());
                config.atoms[k].covalent=config.atoms[i].covalent;
                config.atoms[k].dotify=config.atoms[i].dotify;
                config.xyz[0][k]=config.xyz[0][i];
                config.xyz[1][k]=config.xyz[1][i];
                config.xyz[2][k]=config.xyz[2][i];
                k++;

            }

        }

        // Add inserted molecule

        for(int i=0;i<newcfg.natms;i++) {

            if(k == config.atoms.length) {
                config.natms=k;
                config.resizeArrays();
            }

            config.atoms[k]=new Element();
            config.atoms[k].znum=newcfg.atoms[i].znum;
            config.atoms[k].zmas=newcfg.atoms[i].zmas;
            config.atoms[k].zchg=newcfg.atoms[i].zchg;
            config.atoms[k].zrad=newcfg.atoms[i].zrad;
            config.atoms[k].zsym=new String(newcfg.atoms[i].zsym);
            config.atoms[k].zcol=new Color(newcfg.atoms[i].zcol.getRGB());
            config.atoms[k].covalent=newcfg.atoms[i].covalent;
            config.atoms[k].dotify=newcfg.atoms[i].dotify;
            config.xyz[0][k]=newcfg.xyz[0][i];
            config.xyz[1][k]=newcfg.xyz[1][i];
            config.xyz[2][k]=newcfg.xyz[2][i];
            k++;

        }
        config.natms=k;

        // write new CONFIG file

        fname="CFGINS."+String.valueOf(numins);
        if(config.configWrite(fname)){
            println("File "+fname+" created");
            println("Number of atoms in "+fname+" : "+config.natms);
            numins++;
        }

        // rebuild molecular structure

        config.structure=new Structure(config);

        // Draw modified structure

        editor.pane.restore();
        return k;
    }

}
