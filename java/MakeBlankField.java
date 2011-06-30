import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;

public class MakeBlankField extends Basic implements ActionListener {
    /*
*********************************************************************

dl_poly/java GUI class to construct a blank FIELD file

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
     */
    public static MakeBlankField job;
    private static GUI home;
    private static JButton make,load,close;
    private static boolean lshow,lgetf;

    // Define the Graphical User Interface

    public MakeBlankField() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        setTitle("Blank FIELD Maker");

        getContentPane().setBackground(art.back);
        getContentPane().setForeground(art.fore);
        setFont(fontMain);
        GridBagLayout grd = new GridBagLayout();
        GridBagConstraints gbc = new GridBagConstraints();
        getContentPane().setLayout(grd);

        gbc.fill=GridBagConstraints.BOTH;

        // Define the Make button

        make = new JButton("Make");
        make.setBackground(art.butn);
        make.setForeground(art.butf);
        fix(make,grd,gbc,0,0,1,1);

        // Pad out

        JLabel pad1 = new JLabel(" ");
        fix(pad1,grd,gbc,0,1,1,1);

        // Define the Load button

        load = new JButton("Load");
        load.setBackground(art.butn);
        load.setForeground(art.butf);
        fix(load,grd,gbc,0,2,1,1);

        // Pad out

        JLabel pad2 = new JLabel(" ");
        fix(pad2,grd,gbc,0,3,1,1);


        // Define the Close button

        close = new JButton("Close");
        close.setBackground(art.butn);
        close.setForeground(art.butf);
        fix(close,grd,gbc,0,4,1,1);

        // Register action buttons

        make.addActionListener(this);
        load.addActionListener(this);
        close.addActionListener(this);

    }

    public MakeBlankField(GUI here) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        // initial values

        println("Activated panel for making blank FIELD file");
        home=here;
        lshow=false;
        lgetf=false;

        // instantiate application

        job=new MakeBlankField();
        job.pack();
        job.setVisible(true);

    }

    public void actionPerformed(ActionEvent e) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        int call,key;
        String arg = (String)e.getActionCommand();
        if (arg.equals("Make")) {
            lgetf=false;
            call=blankField();
        }
        else if (arg.equals("Load")) {
            lgetf=true;
            lshow=true;
            call=blankField();
        }
        else if (arg.equals("Close")) {
            job.setVisible(false);
        }
    }

    int blankField() {
        /*
*********************************************************************

dl_poly/java routine to read a DL_POLY CONFIG file
and construct a  corresponding  blank DL_POLY FIELD file

copyright - daresbury laboratory
author    - w.smith december 2000

*********************************************************************
         */
        String fldfile,cfgfile;
	Structure struct;
        Molecule mol;

        // read the CONFIG file

        if(lgetf) config=getConfig(home,"CFG");
        if(config==null || config.natms==0)return 0;

	// check if structure has been defined

	if(config.structure == null)
	    config.structure=new Structure(config);

        if(config.structure.nmols <= 0) return -1;
	struct=config.structure;

        // write the contiguized CONFIG file

        cfgfile="CFGBLK."+numblk;
        config.configWrite(cfgfile);

        // display the CONFIG file

        if(lshow){
	    if(home.editor != null)
		home.editor.job.setVisible(false);
	    home.editor=new Editor(home);
	    lshow=false;
	}

        // write the FIELD file

        fldfile="FLDBLK."+numblk;
        numblk++;
        try {
            DataOutputStream outStream = new DataOutputStream(new FileOutputStream(fldfile));

            //	FIELD file header

            outStream.writeBytes(config.title+"\n");
            outStream.writeBytes("UNITS "+"\n");
            outStream.writeBytes("MOLECULES "+BML.fmt(struct.nmoltp,8)+"\n");

            for(int m=0;m<struct.nmoltp;m++) {
                mol=struct.molecules[m];
                outStream.writeBytes(mol.molname +"\n");
                outStream.writeBytes("NUMMOLS "+BML.fmt(struct.msz[m],8)+"\n");

                // atomic data

                outStream.writeBytes("ATOMS"+BML.fmt(mol.natm,8)+"\n");

                for(int i=0;i<mol.natm;i++) {
                    outStream.writeBytes(BML.fmt(mol.atom[i].zsym,8)+BML.fmt(mol.atom[i].zmas,12)+
                    BML.fmt(0.0,12)+"\n");
                }

                // bond data

                if(mol.nbnd > 0) {
                    outStream.writeBytes("BONDS"+BML.fmt(mol.nbnd,8)+"\n");

                    for(int i=0;i<mol.nbnd;i++) {
                        outStream.writeBytes("TYPE"+BML.fmt(mol.ibnd[0][i]+1,5)+BML.fmt(mol.ibnd[1][i]+1,5)+
                        BML.fmt(0.0,12)+BML.fmt(0.0,12)+BML.fmt(0.0,12)+"\n");
                    }
                }

                // angle data

                if(mol.nang > 0) {
                    outStream.writeBytes("ANGLES"+BML.fmt(mol.nang,8)+"\n");

                    for (int i=0;i<mol.nang;i++) {
                        outStream.writeBytes("TYPE"+BML.fmt(mol.iang[0][i]+1,5)+BML.fmt(mol.iang[1][i]+1,5)+
                        BML.fmt(mol.iang[2][i]+1,5)+BML.fmt(0.0,12)+BML.fmt(0.0,12)+
                        BML.fmt(0.0,12)+"\n");
                    }
                }

                // dihedral data

                if(mol.ndih>0) {
                    outStream.writeBytes("DIHEDRALS"+BML.fmt(mol.ndih,8)+"\n");

                    for (int i=0;i<mol.ndih;i++) {
                        outStream.writeBytes("TYPE"+BML.fmt(mol.idih[0][i]+1,5)+BML.fmt(mol.idih[1][i]+1,5)+
                        BML.fmt(mol.idih[2][i]+1,5)+BML.fmt(mol.idih[3][i]+1,5)+
                        BML.fmt(0.0,12)+BML.fmt(0.0,12)+BML.fmt(0.0,12)+
                        BML.fmt(0.0,12)+BML.fmt(0.0,12)+"\n");
                    }
                }

                // inversion data

                if(mol.ninv>0) {
                    outStream.writeBytes("INVERSIONS"+BML.fmt(mol.ninv,8)+"\n");

                    for (int i=0;i<mol.ninv;i++) {
                        outStream.writeBytes("TYPE"+BML.fmt(mol.invr[0][i]+1,5)+BML.fmt(mol.invr[1][i]+1,5)+
                        BML.fmt(mol.invr[2][i]+1,5)+BML.fmt(mol.invr[3][i]+1,5)+
                        BML.fmt(0.0,12)+BML.fmt(0.0,12)+"\n");
                    }
                }

                outStream.writeBytes("FINISH"+"\n");

            }

            // VDW data

            outStream.writeBytes("VDW"+BML.fmt((struct.nunq*(struct.nunq+1))/2,8)+"\n");

            for (int i=0;i<struct.nunq;i++) {
                for(int j=i;j<struct.nunq;j++) {
                    outStream.writeBytes(BML.fmt(struct.unqatm[i],8)+BML.fmt(struct.unqatm[j],8)+"TYPE"+
                    BML.fmt(0.0,12)+BML.fmt(1.0,12)+BML.fmt(0.0,12)+"\n");
                }
            }

            // Three body forces

            outStream.writeBytes("TBP"+BML.fmt(struct.nunq*(struct.nunq*(struct.nunq+1))/2,8)+"\n");

            for(int j=0;j<struct.nunq;j++) {
                for(int i=0;i<struct.nunq;i++) {
                    for(int k=i;k<struct.nunq;k++) {
                        outStream.writeBytes(BML.fmt(struct.unqatm[i],8)+BML.fmt(struct.unqatm[j],8)+
                        BML.fmt(struct.unqatm[k],8)+"TYPE"+BML.fmt(0.0,12)+
                        BML.fmt(0.0,12)+BML.fmt(0.0,12)+BML.fmt(0.0,12)+
                        BML.fmt(0.0,12)+"\n");
                    }
                }
            }

            outStream.writeBytes("CLOSE"+"\n");
            println("FIELD file "+fldfile+" created");
            outStream.close();
        }
        catch(Exception e) {
            println("error - writing file: "+fldfile);
        }
        return 0;
    }
}

