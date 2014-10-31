import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;

public class MakeCeramField extends Basic implements ActionListener {
        /*
*********************************************************************

dl_poly/java GUI class to make a ceramics FIELD file

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
    public static MakeCeramField job;
    private static GUI home;
    private static String ffield;
    private static JComboBox<String> cpot;
    private static JButton make,load,close;
    private static JCheckBox display,tetra;
    private static JLabel lab1,pad;
    private static boolean ltetr,lshow,lgetf;
    private static String[] unqatm,name;
    private static double[] mass,aaa,rho,ccc,cga,cgb,spa,spb,rad;
    private static int natms,imcon,nunq,nmols,nrept,nshl,ntatm;

    // Define the Graphical User Interface

    public MakeCeramField() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        setTitle("Ceramic FIELD Maker");

        getContentPane().setBackground(art.back);
        getContentPane().setForeground(art.fore);
        setDefaultCloseOperation(DISPOSE_ON_CLOSE);
        setFont(fontMain);
        GridBagLayout grd = new GridBagLayout();
        GridBagConstraints gbc = new GridBagConstraints();
        getContentPane().setLayout(grd);

        gbc.fill=GridBagConstraints.BOTH;

        //        Define the Make button

        make = new JButton("Make");
        make.setBackground(art.butn);
        make.setForeground(art.butf);
        fix(make,grd,gbc,0,0,1,1);

        //        Define the Load button

        load = new JButton("Load");
        load.setBackground(art.butn);
        load.setForeground(art.butf);
        fix(load,grd,gbc,1,0,1,1);
        pad=new JLabel("xxxxxxxxxxxxxxxxxx");
        pad.setForeground(art.back);
        fix(pad,grd,gbc,1,1,1,1);

        //        Bond type specification

        lab1 = new JLabel("Force field :",JLabel.LEFT);
        fix(lab1,grd,gbc,0,2,1,1);
        cpot = new JComboBox<String>();
        cpot.setBackground(art.scrn);
        cpot.setForeground(art.scrf);
        cpot.addItem("LC_a");
        cpot.addItem("LC_b");
        cpot.addItem("LC_c1");
        cpot.addItem("LC_c2");
        cpot.addItem("GULP");
        fix(cpot,grd,gbc,1,2,1,1);
        fix(new JLabel(" "),grd,gbc,1,3,1,1);

        //        Tetragonal option

        tetra = new JCheckBox("Tetragonal ion");
        tetra.setBackground(art.back);
        tetra.setForeground(art.fore);
        fix(tetra,grd,gbc,0,4,1,1);

        //        Display option

        display = new JCheckBox("Display Cfg.");
        display.setBackground(art.back);
        display.setForeground(art.fore);
        fix(display,grd,gbc,0,5,1,1);
        fix(new JLabel(" "),grd,gbc,0,6,1,1);

        //        Define the Close button

        close = new JButton("Close");
        close.setBackground(art.butn);
        close.setForeground(art.butf);
        fix(close,grd,gbc,0,7,1,1);

        // Register action buttons

        make.addActionListener(this);
        load.addActionListener(this);
        close.addActionListener(this);

    }

    public MakeCeramField(GUI here) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        // initial values

        println("Activated panel for making ceramic FIELD files");
        home=here;
        ltetr=false;
        lshow=true;
        if(config==null)
            natms=0;
        else
            natms=config.natms;

        // instantiate application

        job=new MakeCeramField();
        job.pack();
        job.setVisible(true);

        // set GUI parameters

        cpot.setSelectedItem("LC_a");
        display.setSelected(lshow);
        tetra.setSelected(ltetr);
    }

    public void actionPerformed(ActionEvent e) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        int call,key;
        String arg = (String)e.getActionCommand();
        if (arg.equals("Make")) {
            lshow=false;
            lgetf=false;
            if(natms==0)lgetf=true;
            ffield=cpot.getSelectedItem().toString();
            ltetr=tetra.isSelected();
            call=ceramField();
        }
        else if (arg.equals("Load")) {
            lgetf=true;
            ffield=cpot.getSelectedItem().toString();
            lshow=display.isSelected();
            ltetr=tetra.isSelected();
            call=ceramField();
        }
        else if (arg.equals("Close")) {
            job.dispose();
        }
    }

    int ceramField() {
        /*
**********************************************************************

dl_poly/java routine to construct a DL_POLY FIELD file for ceramics based
on the contents of the CONFIG file and the CERAMICS force field file

copyright - daresbury laboratory
author    - w.smith january 2001

**********************************************************************
         */

        Structure struct;
        String fldfile,cfgfile;
        boolean lnew,kill,lmols,lsha,lshb;
        double det,xd,yd,zd,ssx,ssy,ssz,xss,yss,zss,rsq;
        double rcell[]=new double[9];
        int k;

        // read the CONFIG file

        if(lgetf) config=getConfig(home,"CFG");
        if(config==null || config.natms==0)return 0;
        natms=config.natms;
        imcon=config.pbc.imcon;

        // check if structure has been defined

        if(config.structure == null)
            config.structure=new Structure(config);

        if(config.structure.nmols <= 0) return -1;
        struct=config.structure;

        // write the contiguized CONFIG file

        cfgfile="CFGCER."+numcer;
        config.configWrite(cfgfile);

        // display the CONFIG file

        if(lshow) {
            if(!editor.isVisible())
                editor.showEditor();
            editor.pane.restore();
        }

        // create arrays

        nmols=struct.nmols;
        nunq=struct.nunq;
        name=struct.name;
        nrept=struct.nrept;
        nshl=struct.nshl;
        ntatm=struct.ntatm;
        unqatm=struct.unqatm;
        aaa=new double[nunq];
        rho=new double[nunq];
        ccc=new double[nunq];
        cga=new double[nunq];
        cgb=new double[nunq];
        spa=new double[nunq];
        spb=new double[nunq];
        mass=new double[nunq];
        rad=new double[natms];

        // import Ceramics field parameters

        if(ceramics() != 0) {
            println("Error - failed in processing CERAMICS file");
            return -1;
        }

        // construct the FIELD file

        fldfile="FLDCER."+numcer;
        numcer++;

        try {
            DataOutputStream outStream = new DataOutputStream(new FileOutputStream(fldfile));

            //  FIELD file header

            outStream.writeBytes(config.title+"\n");
            outStream.writeBytes("UNITS eV"+"\n");
            outStream.writeBytes("MOLECULES 1\n");
            outStream.writeBytes("FORCEFIELD "+ffield+"\n");
            outStream.writeBytes("NUMMOLS"+BML.fmt(nmols,6)+"\n");

            // atomic data

            outStream.writeBytes("ATOMS"+BML.fmt(nrept,8)+"\n");

            for(int i=0;i<nrept;i++) {
                k=i;
                for(int j=0;j<nunq;j++) {
                    if(name[i].equals(unqatm[j])) k=j;
                }
                outStream.writeBytes(name[i]+BML.fmt(mass[k],12)+BML.fmt(cga[k],12)+
                BML.fmt(1,5)+BML.fmt(0,5)+BML.fmt(0,5)+"\n");
            }

            // shell data

            if(nshl > 0) {
                outStream.writeBytes("SHELL"+BML.fmt(nshl,8)+"\n");
                rcell=AML.invert(config.pbc.cell);

                for(int i=0;i<nrept-1;i++) {
                    if(!name[i].substring(4,8).equals("    ")) {
                        k=i;
                        for(int j=0;j<nunq;j++) {
                            if(unqatm[j].equals(name[i]))k=j;
                        }
                        for(int j=i+1;j<nrept;j++) {
                            if(name[i].substring(0,4).equals(name[j].substring(0,4))) {
                                xd=config.xyz[0][i]-config.xyz[0][j];
                                yd=config.xyz[1][i]-config.xyz[1][j];
                                zd=config.xyz[2][i]-config.xyz[2][j];
                                if(imcon > 0) {
                                    ssx=(rcell[0]*xd+rcell[3]*yd+rcell[6]*zd);
                                    ssy=(rcell[1]*xd+rcell[4]*yd+rcell[7]*zd);
                                    ssz=(rcell[2]*xd+rcell[5]*yd+rcell[8]*zd);
                                    xss=ssx-BML.nint(ssx);
                                    yss=ssy-BML.nint(ssy);
                                    zss=ssz-BML.nint(ssz);
                                    xd=(config.pbc.cell[0]*xss+config.pbc.cell[3]*yss+config.pbc.cell[6]*zss);
                                    yd=(config.pbc.cell[1]*xss+config.pbc.cell[4]*yss+config.pbc.cell[7]*zss);
                                    zd=(config.pbc.cell[2]*xss+config.pbc.cell[5]*yss+config.pbc.cell[8]*zss);
                                }
                                rsq=xd*xd+yd*yd+zd*zd;
                                if(rsq < 1.0) {
                                    outStream.writeBytes(BML.fmt(i+1,5)+BML.fmt(j+1,5)+BML.fmt(spa[k],12));
                                }
                            }
                        }
                    }
                }
            }
            outStream.writeBytes("FINISH\n");

            // VDW data

            outStream.writeBytes("VDW"+BML.fmt((nunq*(nunq+1))/2,8)+"\n");

            for(int i=0;i<nunq;i++) {
                lsha=(unqatm[i].toLowerCase().substring(4,6).equals("_s")) ||
                (unqatm[i].substring(4).equals("    "));
                for(int j=i;j<nunq;j++) {
                    lshb=(unqatm[j].toLowerCase().substring(4,6).equals("_s")) ||
                    (unqatm[j].substring(4).equals("    "));
                    if(lsha && lshb) {
                        if((i == j) && (unqatm[i].substring(0,4).equals("O_2-"))) {

                            outStream.writeBytes(BML.fmt(unqatm[i],8)+BML.fmt(unqatm[i],8)+"buck"+
                            BML.fmt(aaa[i],12)+BML.fmt(rho[i],12)+BML.fmt(ccc[i],12)+"\n");
                        }
                        else if(unqatm[i].substring(0,4).equals("O_2-")) {
                            outStream.writeBytes(BML.fmt(unqatm[j],8)+BML.fmt(unqatm[i],8)+"buck"+
                            BML.fmt(aaa[j],12)+BML.fmt(rho[j],12)+BML.fmt(ccc[j],12)+"\n");
                        }
                        else if(unqatm[j].substring(0,4).equals("O_2-")) {
                            outStream.writeBytes(BML.fmt(unqatm[i],8)+BML.fmt(unqatm[j],8)+"buck"+
                            BML.fmt(aaa[i],12)+BML.fmt(rho[i],12)+BML.fmt(ccc[i],12)+"\n");
                        }
                        else {
                            outStream.writeBytes(BML.fmt(unqatm[i],8)+BML.fmt(unqatm[j],8)+"buck"+
                            BML.fmt(0.0,12)+BML.fmt(1.0,12)+BML.fmt(0.0,12)+"\n");
                        }
                    }
                    else {
                        outStream.writeBytes(BML.fmt(unqatm[i],8)+BML.fmt(unqatm[j],8)+"buck"+
                        BML.fmt(0.0,12)+BML.fmt(1.0,12)+BML.fmt(0.0,12)+"\n");
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

    int ceramics() {
        /*
*********************************************************************

dl_poly/java utility for reading the ceramics force field

copyright daresbury laboratory
author w.smith january 2001

*********************************************************************
         */
        int status,last,i,nas;
        boolean trip;

        status=0;

        // open CERAMICS file and skip to required force field

        try {
            String record="";
            InputStream instream = this.getClass().getResourceAsStream("CERAMICS");
            InputStreamReader isr = new InputStreamReader(instream);
            BufferedReader reader = new BufferedReader(isr);
            println("Reading file: CERAMICS");
            OUT2:
                while((record=reader.readLine()) != null) {
                    if(record.indexOf(ffield) == 1)break OUT2;
                }

                // read CERAMICS contents

                i=0;
                nas=0;
                while((record=reader.readLine()) != null && record.charAt(0) != '#') {
                    trip=!record.substring(4,6).equals("_t");
                    if(trip || ltetr) {
                        for(int j=0;j<nunq;j++) {
                            if(record.substring(0,4).equals(unqatm[j].substring(0,4))) {
                                if(trip)nas++;
                                if(unqatm[j].toLowerCase().substring(4,6).equals("_c")) {
                                    mass[j]=0.9*BML.giveDouble(record,2);
                                }
                                else if(unqatm[j].toLowerCase().substring(4,6).equals("_s")) {
                                    mass[j]=0.1*BML.giveDouble(record,2);
                                }
                                else {
                                    mass[j]=BML.giveDouble(record,2);
                                }
                                aaa[j]=BML.giveDouble(record,3);
                                rho[j]=BML.giveDouble(record,4);
                                ccc[j]=BML.giveDouble(record,5);
                                cga[j]=BML.giveDouble(record,6);
                                spa[j]=BML.giveDouble(record,7);
                                cgb[j]=BML.giveDouble(record,8);
                                spb[j]=BML.giveDouble(record,9);
                                i++;
                            }
                        }
                    }
                }
                last=i;
                println("CERAMICS File successfully loaded");
                reader.close();
        }
        catch(FileNotFoundException e) {
            println("Error - file not found: CERAMICS");
            return -1;
        }
        catch(Exception e) {
            println("Error reading file: CERAMICS" + " "+e);
            return -2;
        }
        if(nas != nunq) {
            println("Error - required atom type not in CERAMICS file");
            return -3;
        }

        // now sort out details of different force fields

        if(ffield.equals("LC_a") || ffield.equals("LC_b")) {
            for(i=0;i<nunq;i++) {

                if(!unqatm[i].substring(4).equals("    ")) {
                    println("Error - unrecognised forcefield atom: "+unqatm[i]);
                    return -4;
                }
            }
        }
        else {
            for(i=0;i<nunq;i++) {
                if(unqatm[i].substring(4).equals("    ")) {
                    cga[i]=strchg(unqatm[i]);
                }
                else if(BML.nint(spa[i]) == 99999) {
                    println("Error - unrecognised forcefield atom: "+unqatm[i]);
                    return -5;
                }
                else if(unqatm[i].toLowerCase().substring(4,6).equals("_c")) {
                    cga[i]=strchg(unqatm[i])-cga[i];
                }
            }
            if(ffield.equals("LC_c1") || ffield.equals("LC_c2")) {
                if(ntatm > 2) {
                    println("Force field not suitable for mixed oxides");
                    return -6;
                }
                for(i=0;i<nunq;i++) {

                    if((!unqatm[i].substring(0,4).equals("O_2-")) &&
                    (unqatm[i].substring(4,8).equals("    ") ||
                    unqatm[i].toLowerCase().substring(4,6).equals("_s"))) {
                        for(int j=0;j<nunq;j++) {
                            if(unqatm[j].substring(0,4).equals("O_2-")) {
                                if(unqatm[j].toLowerCase().substring(4,6).equals("_s")) {
                                    cga[j]=cgb[i];
                                    spa[j]=spb[i];
                                }
                                else {
                                    cga[j]=strchg(unqatm[j])-cgb[i];
                                    spa[j]=spb[i];
                                }
                            }
                        }
                    }
                }
            }
        }
        return 0;
    }

    double strchg(String sss) {
        /*
*********************************************************************

dl_poly/java GUI routine  to determine ionic charge from ion name

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        if(sss.substring(2,4).equals("2-")) {
            return -2.0;
        }
        else if(sss.substring(2,4).equals("1-")) {
            return -1.0;
        }
        else if(sss.substring(2,4).equals("_+")) {
            return 1.0;
        }
        else if(sss.substring(2,4).equals("2+")) {
            return 2.0;
        }
        else if(sss.substring(2,4).equals("3+")) {
            return 3.0;
        }
        else if(sss.substring(2,4).equals("4+")) {
            return 4.0;
        }
        else if(sss.substring(2,4).equals("5+")) {
            return 5.0;
        }
        return 0.0;
    }
}
