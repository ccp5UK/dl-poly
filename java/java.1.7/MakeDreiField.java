import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;

public class MakeDreiField extends Basic implements ActionListener {
    /*
*********************************************************************

dl_poly/java GUI class to make Dreiding Field files

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
     */
    public static MakeDreiField job;
    private static GUI home;
    private static JComboBox<String> bond,vdwt;
    private static JButton make,load,close;
    private static JCheckBox charges;
    private static JLabel lab1,lab2,pad1,pad2,pad3;
    private static boolean lbuck,lmorse,lchg,lshow,lgetf,lfixed;
    private static String[] dsym,kbnd,kang,kdih,kinv,kvdw,khbd;
    private static String[][] idvdw,idhbd;
    private static double[] dmss,dchg,drad,dang,deps,dzet,dsig;
    private static double[] bord,eig;
    private static double[][] pang,pbnd,pdih,pinv,phbd,pvdw,con,vec;
    private static int mxpbnd,mxpang,mxpdih,mxpinv,mxphbd,mxpvdw,mxvdw,lhbnd;
    private static int[] ident,jdent,ida,lab;
    private static int nvdw,nhbd,ndrei;
    private static Structure struct;

    // Define the Graphical User Interface

    public MakeDreiField() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        setTitle("Dreiding FIELD Maker");

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

        //        Pad out

        pad1 = new JLabel(" ");
        fix(pad1,grd,gbc,0,1,1,1);

        //        Bond type specification

        lab1 = new JLabel("Bond type :",JLabel.LEFT);
        fix(lab1,grd,gbc,0,2,1,1);
        bond = new JComboBox<String>();
        bond.setBackground(art.scrn);
        bond.setForeground(art.scrf);
        bond.addItem("Harmonic");
        bond.addItem("Morse");
        bond.addItem("Fixed");
        fix(bond,grd,gbc,1,2,1,1);

        //        VDW type specification

        lab2 = new JLabel("VDW type :",JLabel.LEFT);
        fix(lab2,grd,gbc,0,3,1,1);
        vdwt = new JComboBox<String>();
        vdwt.setBackground(art.scrn);
        vdwt.setForeground(art.scrf);
        vdwt.addItem("Len-Jones");
        vdwt.addItem("Buckingham");
        fix(vdwt,grd,gbc,1,3,1,1);

        //        Pad out

        pad2 = new JLabel(" ");
        fix(pad2,grd,gbc,0,4,1,1);

        //        Charges option

        charges = new JCheckBox("Use Charges");
        charges.setBackground(art.back);
        charges.setForeground(art.fore);
        fix(charges,grd,gbc,0,5,1,1);

        //        Pad out

        pad3 = new JLabel(" ");
        fix(pad3,grd,gbc,0,6,1,1);

        //        Define the Close button

        close = new JButton("Close");
        close.setBackground(art.butn);
        close.setForeground(art.butf);
        fix(close,grd,gbc,1,8,1,1);

        // Register action buttons

        make.addActionListener(this);
        load.addActionListener(this);
        close.addActionListener(this);

    }

    public MakeDreiField(GUI here) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        // initial values

        println("Activated panel for making Dreiding FIELD files");
        home=here;
        lfixed=false;
        lbuck=false;
        lmorse=false;
        lchg=false;
        lshow=false;
        lgetf=false;
        mxpbnd=3;
        mxpang=3;
        mxpdih=5;
        mxpinv=2;
        mxpvdw=3;
        mxphbd=5;

        // instantiate application

        job=new MakeDreiField();
        job.pack();
        job.setVisible(true);

        // set GUI parameters

        bond.setSelectedItem("Harmonic");
        vdwt.setSelectedItem("Len-Jones");
        charges.setSelected(lchg);
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
            lbuck=false;
            lmorse=false;
            lfixed=false;
            key=bond.getSelectedIndex();
            if(key == 1) lmorse=true;
            if(key == 2) lfixed=true;
            key=vdwt.getSelectedIndex();
            if(key > 0) lbuck=true;
            lchg=charges.isSelected();
            call=dreiField();
        }
        else if (arg.equals("Load")) {
            lshow=true;
            lgetf=true;
            lbuck=false;
            lmorse=false;
            lfixed=false;
            key=bond.getSelectedIndex();
            if(key == 1) lmorse=true;
            if(key == 2) lfixed=true;
            key=vdwt.getSelectedIndex();
            if(key > 0) lbuck=true;
            lchg=charges.isSelected();
            call=dreiField();
        }
        else if (arg.equals("Close")) {
            job.dispose();
        }
    }
    int dreiField() {
        /*
*********************************************************************

dl_poly/java routine to read a DL_POLY CONFIG file
and construct a  corresponding  DL_POLY FIELD file
using the Dreiding force model

copyright - daresbury laboratory
author    - w.smith december 2000

*********************************************************************
         */
        int newinv;
        String fldfile,cfgfile;
        Molecule mol;

        // read the CONFIG file

        if(lgetf) config=getConfig(home,"CFG");
        if(config==null || config.natms==0)return 0;

        if(config.structure == null)
            config.structure=new Structure(config);

        if(config.structure.nmols <= 0) return -1;
        struct=config.structure;

        // write the contiguized CONFIG file

        cfgfile="CFGDRE."+numdre;
        config.configWrite(cfgfile);

        // import Dreiding field parameters

        if(dreiding() != 0) {
            println("Error - failed in processing MINIDREI file");
            return -1;
        }

        // define system structure

        struct=new Structure(config);
        if(struct.nmols <= 0) return -2;

        // attach Dreiding identity to dl_poly atoms

        if(attach() < 0) return -3;

        // display the CONFIG file

        if(lshow){
            if(!editor.isVisible())
                editor.showEditor();
            editor.pane.restore();
            lshow=false;
        }

        // construct the FIELD file

        fldfile="FLDDRE."+numdre;
        numdre++;

        OUT:
            try {
                DataOutputStream outStream = new DataOutputStream(new FileOutputStream(fldfile));

                //      FIELD file header

                outStream.writeBytes(config.title+"\n");
                outStream.writeBytes("UNITS kcal"+"\n");
                outStream.writeBytes("MOLECULES "+BML.fmt(struct.nmoltp,8)+"\n");

                for(int m=0;m<struct.nmoltp;m++) {
                    mol=struct.molecules[m];
                    outStream.writeBytes(mol.molname +"\n");
                    outStream.writeBytes("NUMMOLS "+BML.fmt(struct.msz[m],8)+"\n");
                    if(identify(mol)<0) break OUT;

                    // atomic data

                    outStream.writeBytes("ATOMS"+BML.fmt(mol.natm,8)+"\n");

                    for(int i=0;i<mol.natm;i++) {
                        outStream.writeBytes(BML.fmt(mol.atom[i].zsym,8)+
                        BML.fmt(mol.atom[i].zmas,12)+
                        BML.fmt(mol.chge[i],12)+"\n");
                    }

                    // estimate bond orders

                    if(border(mol)!= 0) {
                        println("Warning - conjugation analysis failed");
                        println("Using default bond orders (1.0)");
                    }

                    // assign Dreiding bond parameters

                    if(mol.nbnd > 0) {
                        bndass(mol);
                        if(lfixed) {
                            outStream.writeBytes("CONSTRAINTS"+BML.fmt(mol.nbnd,8)+"\n");
                            for(int i=0;i<mol.nbnd;i++) {
                                outStream.writeBytes(BML.fmt(mol.ibnd[0][i]+1,5)+
                                BML.fmt(mol.ibnd[1][i]+1,5)+
                                BML.fmt(pbnd[1][i],12)+"\n");
                            }
                        }
                        else {
                            outStream.writeBytes("BONDS"+BML.fmt(mol.nbnd,8)+"\n");

                            for(int i=0;i<mol.nbnd;i++) {
                                outStream.writeBytes(kbnd[i]+BML.fmt(mol.ibnd[0][i]+1,5)+
                                BML.fmt(mol.ibnd[1][i]+1,5)+
                                BML.fmt(pbnd[0][i],12)+BML.fmt(pbnd[1][i],12)+
                                BML.fmt(pbnd[2][i],12)+"\n");
                            }
                        }
                    }

                    // assign Dreiding angle parameters

                    if(mol.nang > 0) {
                        angass(mol);
                        outStream.writeBytes("ANGLES"+BML.fmt(mol.nang,8)+"\n");

                        for (int i=0;i<mol.nang;i++) {
                            outStream.writeBytes(kang[i]+BML.fmt(mol.iang[0][i]+1,5)+
                            BML.fmt(mol.iang[1][i]+1,5)+
                            BML.fmt(mol.iang[2][i]+1,5)+
                            BML.fmt(pang[0][i],12)+
                            BML.fmt(pang[1][i],12)+
                            BML.fmt(pang[2][i],12)+"\n");
                        }
                    }

                    // assign Dreiding dihedral angles

                    if(mol.ndih>0) {
                        dihass(mol);
                        outStream.writeBytes("DIHEDRALS"+BML.fmt(mol.ndih,8)+"\n");

                        for (int i=0;i<mol.ndih;i++) {
                            outStream.writeBytes(kdih[i]+BML.fmt(mol.idih[0][i]+1,5)+
                            BML.fmt(mol.idih[1][i]+1,5)+
                            BML.fmt(mol.idih[2][i]+1,5)+
                            BML.fmt(mol.idih[3][i]+1,5)+
                            BML.fmt(pdih[0][i],12)+
                            BML.fmt(pdih[1][i],12)+
                            BML.fmt(pdih[2][i],12)+
                            BML.fmt(pdih[3][i],12)+
                            BML.fmt(pdih[4][i],12)+"\n");
                        }
                    }

                    // assign Dreiding inversion angles

                    newinv=0;
                    if(mol.ninv>0) newinv=invass(mol);
                    if(newinv>0) {
                        outStream.writeBytes("INVERSIONS"+BML.fmt(newinv,8)+"\n");

                        for (int i=0;i<mol.ninv;i++) {
                            if(!kinv[i].equals("null"))
                                outStream.writeBytes(kinv[i]+
                                BML.fmt(mol.invr[0][i]+1,5)+
                                BML.fmt(mol.invr[1][i]+1,5)+
                                BML.fmt(mol.invr[2][i]+1,5)+
                                BML.fmt(mol.invr[3][i]+1,5)+
                                BML.fmt(pinv[0][i],12)+
                                BML.fmt(pinv[1][i],12)+"\n");
                        }
                    }

                    outStream.writeBytes("FINISH"+"\n");
                }

                // assign Dreiding van der waals potentials

                nvdw=vdwass();
                outStream.writeBytes("VDW"+BML.fmt(nvdw,8)+"\n");

                for (int i=0;i<nvdw;i++) {
                    outStream.writeBytes(BML.fmt(idvdw[0][i],8)+BML.fmt(idvdw[1][i],8)+kvdw[i]+
                    BML.fmt(pvdw[0][i],12)+BML.fmt(pvdw[1][i],12)+
                    BML.fmt(pvdw[2][i],12)+"\n");
                }

                // assign Dreiding hydrogen bonds

                nhbd=hbdass();

                if(nhbd>0) {
                    outStream.writeBytes("TBP"+BML.fmt(nhbd,8)+"\n");

                    for(int i=0;i<nhbd;i++) {
                        outStream.writeBytes(BML.fmt(idhbd[0][i],8)+BML.fmt(idhbd[1][i],8)+
                        BML.fmt(idhbd[2][i],8)+
                        khbd[i]+BML.fmt(phbd[0][i],12)+BML.fmt(phbd[1][i],12)+
                        BML.fmt(phbd[2][i],12)+BML.fmt(phbd[3][i],12)+
                        BML.fmt(phbd[4][i],12)+"\n");
                    }
                }

                outStream.writeBytes("CLOSE"+"\n");
                println("FIELD file "+fldfile+" created");
                outStream.close();
            }
            catch(Exception e) {
                println("error - writing file: "+fldfile+" "+e);
            }
            return 0;
    }

    int dreiding() {
        /*
*********************************************************************

dl_poly/java utility for reading the dreiding force field

copyright daresbury laboratory
author w.smith december 2000

*********************************************************************
         */
        int mxd=100;
        dsym=new String[mxd];
        dmss=new double[mxd];
        dchg=new double[mxd];
        drad=new double[mxd];
        dang=new double[mxd];
        deps=new double[mxd];
        dsig=new double[mxd];
        dzet=new double[mxd];

        // open MINIDREI file and skip header records

        try {
            String record="";
            InputStream instream = this.getClass().getResourceAsStream("MINIDREI");
            InputStreamReader isr = new InputStreamReader(instream);
            BufferedReader reader = new BufferedReader(isr);
            println("Reading file: MINIDREI");

            for(int i=0;i<6;i++) {
                record = reader.readLine();
            }


            // read MINIDREI contents

            int i=0;
            while((record=reader.readLine()) != null) {
                dsym[i]=BML.giveWord(record,1);
                dmss[i]=BML.giveDouble(record,2);
                dchg[i]=BML.giveDouble(record,3);
                drad[i]=BML.giveDouble(record,4);
                dang[i]=BML.giveDouble(record,5);
                deps[i]=BML.giveDouble(record,6);
                dsig[i]=BML.giveDouble(record,7);
                dzet[i]=BML.giveDouble(record,8);
                i++;
            }
            ndrei=i;
            println("MINIDREI File successfully loaded");
            reader.close();
        }
        catch(FileNotFoundException e) {
            println("Error - file not found: MINIDREI");
            return -1;
        }
        catch(Exception e) {
            println("Error reading file: MINIDREI" + " "+e);
            return -2;
        }
        return 0;
    }

    int attach() {
        /*
*********************************************************************

dl_poly/java utility for attaching the dreiding force field identity
to a dl_poly atom

copyright daresbury laboratory
author w.smith october 2001

*********************************************************************
         */
        int status=0;
        ident=new int[struct.nunq];

        // attach link to dl_poly identity

        lhbnd=-1;
        for(int i=0;i<struct.nunq;i++) {
            ident[i]=-1;
            for(int j=0;j<ndrei;j++) {
                if(struct.unqatm[i].indexOf(dsym[j].trim())>=0) {
                    ident[i]=j;
                }
            }
            if(ident[i]<0) {
                println("Error - unidentified atom: "+struct.unqatm[i]);
                status=-1;
            }
            if(struct.unqatm[i].indexOf("H__HB")>=0) lhbnd=i;
        }
        return status;
    }

    int identify(Molecule mol) {
        /*
*********************************************************************

dl_poly/java utility for attaching the dreiding force field identity
to a dl_poly atom within a molecule

copyright daresbury laboratory
author w.smith october 2001

*********************************************************************
         */
        int status=0;
        jdent=new int[mol.natm];

        // attach link to dl_poly identity

        for(int i=0;i<mol.natm;i++) {
            jdent[i]=-1;
            for(int j=0;j<ndrei;j++) {
                if(mol.atom[i].zsym.indexOf(dsym[j].trim())>=0) {
                    jdent[i]=j;
                }
            }
            if(jdent[i]<0) {
                println("Error - unidentified atom: "+mol.atom[i].zsym);
                status=-1;
            }
            else {
                // assign Dreiding masses, charges and radii

                mol.atom[i].zrad=drad[jdent[i]];
                mol.atom[i].zmas=dmss[jdent[i]];
                if(lchg) {
                    mol.chge[i]=dchg[jdent[i]];
                }
                else {
                    mol.chge[i]=0.0;
                }
            }
        }
        return status;
    }

    int border(Molecule mol) {
        /*
*********************************************************************

dl_poly/java calculation of bond order in conjugated systems

copyright daresbury laboratory
author w.smith october 2000

*********************************************************************
         */
        String name;
        int ii,jj,jbnd,nmpi,status;

        ida=new int[mol.natm];
        bord=new double[mol.nbnd];

        nmpi=0;
        status=0;

        // set default bond orders

        for (int i=0;i<mol.nbnd;i++) {
            bord[i]=1.0;
        }

        // identify SP2 atoms

        jbnd=0;
        for (int i=0;i<mol.natm;i++) {
            ida[i]=-1;
            name=mol.atom[i].zsym;

            if(name.charAt(2)=='2' || name.charAt(2)=='R') {
                ida[i]=jbnd;
                jbnd++;
                if(name.charAt(0)=='C') {
                    nmpi++;
                }
                if(name.charAt(0)=='N' && mol.list[i]==2) {
                    nmpi++;
                }
                if(name.charAt(0)=='N' && mol.list[i]==3) {
                    nmpi+=2;
                }
                if(name.charAt(0)=='O' && mol.list[i]==1) {
                    nmpi++;
                }
                if(name.charAt(0)=='O' && mol.list[i]==2) {
                    nmpi+=2;
                }
            }
        }
        println("Number of SP2 atoms found:" + BML.fmt(jbnd,6));
        println("Number PI electrons found:" + BML.fmt(nmpi,6));

        // initialize huckel matrix

        println("Starting conjugated bond analysis.....");

        lab=new int[jbnd];
        eig=new double[jbnd];
        con=new double[jbnd][jbnd];
        vec=new double[jbnd][jbnd];

        for (int i=0;i<jbnd;i++) {
            for (int j=0;j<jbnd;j++) {
                con[i][j]=0.0;
                vec[i][j]=0.0;
            }
            con[i][i]=-2.0;
            vec[i][i]=1.0;
        }

        // construct huckel matrix

        for (int k=0;k<mol.nbnd;k++) {
            ii=ida[mol.ibnd[0][k]];
            jj=ida[mol.ibnd[1][k]];
            if(ii >= 0 && jj >= 0) {
                con[ii][jj]=-1.0;
                con[jj][ii]=-1.0;
            }
        }

        // diagonalise huckel matrix

        status=AML.Jacobi(jbnd,con,vec);
        if(status != 0) {
            println("Error - problem diagonalising bond order matrix");
            return -1;
        }

        // load eigenvalues

        for (int i=0;i<jbnd;i++) {
            eig[i]=con[i][i];
            lab[i]=i;
        }

        // sort eigenvalues

        AML.ShellSort(jbnd,lab,eig);

        // construct bond orders

        if(jbnd%2 > 0)
            println("Warning - possible uneven conjugation found");
        for (int k=0;k<mol.nbnd;k++) {
            ii=ida[mol.ibnd[0][k]];
            jj=ida[mol.ibnd[1][k]];

            if(ii >= 0 && jj >= 0) {
                for (int i=0;i<nmpi/2;i++) {
                    bord[k]+=2.0*vec[ii][lab[i]]*vec[jj][lab[i]];
                }
                if(nmpi%2==1)
                    bord[k]+=vec[ii][lab[nmpi/2+1]]*vec[jj][lab[nmpi/2+1]];

            }
        }

        // rationalise bond orders

        for (int k=0;k<mol.nbnd;k++) {
            ii=mol.ibnd[0][k];
            jj=mol.ibnd[1][k];

            if(bord[k] < 1.25)
                bord[k]=1.0;
            else if(bord[k] < 1.75)
                bord[k]=1.5;
            else
                bord[k]=2.0;
        }
        println("Finished conjugated bond analysis");
        return 0;
    }

    void bndass(Molecule mol) {
        /*
*********************************************************************

dl_poly/java routine for assigning bond parameters in the
dreiding force field

copyright daresbury laboratory
author w.smith october 2000

*********************************************************************
         */
        String nami,namj;
        int i,j;
        kbnd=new String[mol.nbnd];
        pbnd=new double[mxpbnd][mol.nbnd];

        for (int k=0;k<mol.nbnd;k++) {
            i=mol.ibnd[0][k];
            j=mol.ibnd[1][k];
            nami=mol.atom[i].zsym;
            namj=mol.atom[j].zsym;
            if(nami.charAt(2)=='1' && namj.charAt(2)=='1')bord[k]=3.0;

            if(lmorse) {
                kbnd[k]="mors";
                pbnd[0][k]=70.0*bord[k];
                pbnd[1][k]=drad[jdent[i]]+drad[jdent[j]]-0.01;
                pbnd[2][k]=Math.sqrt(5.0);
            }
            else {
                kbnd[k]="harm";
                pbnd[0][k]=700.0*bord[k];
                pbnd[1][k]=drad[jdent[i]]+drad[jdent[j]]-0.01;
                pbnd[2][k]=0.0;
            }
        }
        println("Bond potentials assigned");
    }

    void angass(Molecule mol) {
        /*
*********************************************************************

dl_poly/java routine for assigning valence angle parameters in the
dreiding force field

copyright daresbury laboratory
author w.smith october 2000

*********************************************************************
         */
        double th0;
        kang=new String[mol.nang];
        pang=new double[mxpang][mol.nang];

        for (int k=0;k<mol.nang;k++) {
            th0=dang[jdent[mol.iang[1][k]]];

            if(Math.abs(th0-180.0) < 0.000001) {
                kang[k]="cos ";
                pang[0][k]=100.0;
                pang[1][k]=0.0;
                pang[2][k]=1.0;
            }
            else {
                kang[k]="hcos";
                pang[0][k]=100.0/Math.pow(Math.sin(Math.PI*th0/180.0),2);
                pang[1][k]=th0;
                pang[2][k]=0.0;
            }
        }
        println("Valence angle potentials assigned");
    }

    void dihass(Molecule mol) {
        /*
*********************************************************************

dl-poly/java routine for assigning dihedral parameters in the
dreiding force field

copyright daresbury laboratory
author w.smith december 2000

*********************************************************************
         */
        String nami,namj,namk,naml;
        double degen;
        boolean lnull=true;
        int i,j,k,l;
        kdih=new String[mol.ndih];
        pdih=new double[mxpdih][mol.ndih];

        for(int m=0;m<mol.ndih;m++) {
            i=mol.idih[0][m];
            j=mol.idih[1][m];
            k=mol.idih[2][m];
            l=mol.idih[3][m];
            nami=mol.atom[i].zsym;
            namj=mol.atom[j].zsym;
            namk=mol.atom[k].zsym;
            naml=mol.atom[l].zsym;
            degen=(double)((mol.list[j]-1)*(mol.list[k]-1));

            if(namj.charAt(2)=='3' && namk.charAt(2)=='3') {
                sp3sp3(m,namj,namk,degen);
            }
            else if((namj.charAt(2)=='2' || namj.charAt(2) =='R') &&
            (namk.charAt(2)=='3')) {
                sp2sp3(m,nami,namj,namk,naml,degen);
            }
            else if((namk.charAt(2)=='2' || namk.charAt(2)=='R') &&
            (namj.charAt(2)=='3')) {
                sp2sp3(m,naml,namk,namj,nami,degen);
            }
            else if(namj.charAt(2)=='2' && namk.charAt(2)=='2') {
                sp2sp2(m,namj,namk,mol.mid[m],degen);
            }
            else if(namj.charAt(2)=='2' && namk.charAt(2)=='R') {
                sp2sp2(m,namj,namk,mol.mid[m],degen);
            }
            else if(namj.charAt(2)=='R' && namk.charAt(2)=='2') {
                sp2sp2(m,namj,namk,mol.mid[m],degen);
            }
            else if(namj.charAt(2)=='R' && namk.charAt(2)=='R') {
                sprspr(m,namj,namk,mol.mid[m],degen);
            }
            else {
                if(lnull) {
                    println("Warning - null dihedral parameters specified");
                    lnull=false;
                }
                kdih[m]="null";
                pdih[0][m]=0.0;
                pdih[1][m]=0.0;
                pdih[2][m]=0.0;
            }
            pdih[3][m]=1.0;
            pdih[4][m]=1.0;

            // check for multiple accounting of 1-4 scaling

            OUT1:

                for(int n=0;n<m-1;n++) {
                    if(i==mol.idih[0][n] && l==mol.idih[3][n]) {
                        if(!kdih[n].equals("null")) {
                            pdih[3][m]=0.0;
                            pdih[4][m]=0.0;
                        }
                        break OUT1;
                    }
                    else if(i==mol.idih[3][n] && l==mol.idih[0][n]) {
                        if(!kdih[n].equals("null")) {
                            pdih[3][m]=0.0;
                            pdih[4][m]=0.0;
                        }
                        break OUT1;
                    }
                }

                OUT2:

                    for(int n=0;n<mol.nang;n++) {

                        if(i==mol.iang[0][n] && l==mol.iang[2][n]) {
                            pdih[3][m]=0.0;
                            pdih[4][m]=0.0;
                            break OUT2;
                        }
                        else if(i==mol.iang[2][n] && l==mol.iang[0][n]) {
                            pdih[3][m]=0.0;
                            pdih[4][m]=0.0;
                            break OUT2;
                        }
                    }
        }
        println("Dihedral angles assigned");
    }

    int invass(Molecule mol) {
        /*
*********************************************************************

dl-poly routine for assigning inversion parameters in the
dreiding force field

copyright daresbury laboratory
author w.smith december 2000

*********************************************************************
         */
        String nami;
        int i,j,k,l,newinv;
        kinv=new String[mol.ninv];
        pinv=new double[mxpinv][mol.ninv];

        newinv=0;
        for(int m=0;m<mol.ninv;m++) {
            i=mol.invr[0][m];
            j=mol.invr[1][m];
            k=mol.invr[2][m];
            l=mol.invr[3][m];
            nami=mol.atom[i].zsym;

            if(nami.charAt(2)=='2' || nami.charAt(2)=='R') {
                kinv[m]="plan";
                pinv[0][m]=40.0;
                pinv[1][m]=0.0;
                newinv++;
            }
            //else if(namj.charAt(2)=='3' && namj.charAt(3)=='1')
            else if(nami.charAt(2)=='3' && nami.charAt(3)=='1') {
                kinv[m]="hcos";
                pinv[1][m]=54.73561033;
                pinv[0][m]=40.0/Math.pow(Math.sin(Math.PI*pinv[1][m]/180.0),2);
                newinv++;
            }
            else {
                kinv[m]="null";
                pinv[0][m]=0.0;
                pinv[1][m]=0.0;
            }
        }
        if(newinv > 0) println("Inversion angles assigned");
        return newinv;
    }
    void sp2sp2(int m,String namj,String namk,int n,double degen) {
        /*
**********************************************************************

dl-poly/java routine for assigning sp2-sp2 dihedral parameters
in the dreiding force field

copyright daresbury laboratory
author w.smith december 2000

*********************************************************************
         */

        kdih[m]="cos ";

        if((namj.indexOf("C_")==0) && (namk.indexOf("O_")==0)) {
            pdih[0][m]=0.5*2.0/degen;
        }
        else if((namj.indexOf("O_")==0) && (namk.indexOf("C_")==0)) {
            pdih[0][m]=0.5*2.0/degen;
        }
        else if(bord[n] < 1.3) {
            pdih[0][m]=0.5*5.0/degen;
        }
        else if(bord[n] < 1.6) {
            pdih[0][m]=0.5*10.0/degen;
        }
        else {
            pdih[0][m]=0.5*45.0/degen;
        }
        pdih[1][m]=180.0;
        pdih[2][m]=2.0;
    }
    void sprspr(int m,String namj,String namk,int n,double degen) {
        /*
*********************************************************************

dl-poly/java routine for assigning spr-spr dihedral parameters
in the dreiding force field

copyright daresbury laboratory
author w.smith december 2000

*********************************************************************
         */


        kdih[m]="cos ";
        if(bord[n] < 1.3) {
            pdih[0][m]=0.5*5.0/degen;
        }
        else if(bord[n] < 1.6) {
            pdih[0][m]=0.5*10.0/degen;
        }
        else {
            pdih[0][m]=0.5*25.0/degen;
        }
        pdih[1][m]=180.0;
        pdih[2][m]=2.0;
    }
    void sp3sp3(int m,String namj,String namk,double degen) {
        /*
*********************************************************************

dl-poly/java routine for assigning sp3-sp3 dihedral parameters
in the dreiding force field

copyright daresbury laboratory
author w.smith december 2000

*********************************************************************
         */

        if((namj.indexOf("C_")==0) && (namk.indexOf("C_")==0)) {
            kdih[m]="cos ";
            pdih[0][m]=0.5*2.0/degen;
            pdih[1][m]=0.0;
            pdih[2][m]=3.0;
        }
        else if((namj.indexOf("C_")==0) && (namk.indexOf("O_")==0)) {
            kdih[m]="cos ";
            pdih[0][m]=0.5*2.0/degen;
            pdih[1][m]=0.0;
            pdih[2][m]=3.0;
        }
        else if((namj.indexOf("O_")==0) && (namk.indexOf("C_")==0)) {
            kdih[m]="cos ";
            pdih[0][m]=0.5*2.0/degen;
            pdih[1][m]=0.0;
            pdih[2][m]=3.0;
        }
        else if((namj.indexOf("C_")==0) && (namk.indexOf("S_")==0)) {
            kdih[m]="cos ";
            pdih[0][m]=0.5*2.0/degen;
            pdih[1][m]=0.0;
            pdih[2][m]=3.0;
        }
        else if((namj.indexOf("S_")==0) && (namk.indexOf("C_")==0)) {
            kdih[m]="cos ";
            pdih[0][m]=0.5*2.0/degen;
            pdih[1][m]=0.0;
            pdih[2][m]=3.0;
        }
        else if((namj.indexOf("C_")==0) && (namk.indexOf("N_")==0)) {
            kdih[m]="cos ";
            pdih[0][m]=0.5*2.0/degen;
            pdih[1][m]=0.0;
            pdih[2][m]=3.0;
        }
        else if((namj.indexOf("N_")==0) && (namk.indexOf("C_")==0)) {
            kdih[m]="cos ";
            pdih[0][m]=0.5*2.0/degen;
            pdih[1][m]=0.0;
            pdih[2][m]=3.0;
        }
        else if((namj.indexOf("O_")==0) && (namk.indexOf("O_")==0)) {
            kdih[m]="cos ";
            pdih[0][m]=0.5*2.0;
            pdih[1][m]=0.0;
            pdih[2][m]=2.0;
        }
        else if((namj.indexOf("S_")==0) && (namk.indexOf("S_")==0)) {
            kdih[m]="cos ";
            pdih[0][m]=0.5*2.0;
            pdih[1][m]=0.0;
            pdih[2][m]=2.0;
        }
        else if((namj.indexOf("P_")==0) && (namk.indexOf("O_")==0)) {
            kdih[m]="cos ";
            pdih[0][m]=0.5*2.0/degen;
            pdih[1][m]=0.0;
            pdih[2][m]=3.0;
        }
        else if((namj.indexOf("O_")==0) && (namk.indexOf("P_")==0)) {
            kdih[m]="cos ";
            pdih[0][m]=0.5*2.0/degen;
            pdih[1][m]=0.0;
            pdih[2][m]=3.0;
        }
        else {
            println("Warning - unknown potential in sp3-sp3 dihedral");
            kdih[m]="unkn";
            pdih[0][m]=0.0;
            pdih[1][m]=0.0;
            pdih[2][m]=0.0;
        }
    }
    void sp2sp3(int m,String nami,String namj,String namk,String naml,double degen) {
        /*
*********************************************************************

dl-poly/java routine for assigning sp2-sp3 and spr-sp3 dihedral
parameters in the dreiding force field

copyright daresbury laboratory
author w.smith december 2000

*********************************************************************
         */

        if(namk.charAt(2)=='3') {
            if((namj.indexOf("C_")==0) && (namk.indexOf("C_")==0)) {
                if((nami.charAt(2)=='2') || (nami.charAt(2)=='R')) {
                    kdih[m]="cos ";
                    pdih[0][m]=0.5*1.0/degen;
                    pdih[1][m]=-180.0;
                    pdih[2][m]=6.0;
                }
                else {
                    kdih[m]="cos ";
                    pdih[0][m]=0.5*2.0/degen;
                    pdih[1][m]=0.0;
                    pdih[2][m]=3.0;
                }
            }
            else if(((namj.indexOf("C_")==0) && (namk.indexOf("O_")==0)) ||
            ((namj.indexOf("O_")==0) && (namk.indexOf("N_")==0))) {
                if((nami.charAt(2)=='2') || (nami.charAt(2)=='R')) {
                    kdih[m]="cos ";
                    pdih[0][m]=0.5*1.0/degen;
                    pdih[1][m]=-180.0;
                    pdih[2][m]=6.0;
                }
                else {
                    kdih[m]="cos ";
                    pdih[0][m]=0.5*2.0/degen;
                    pdih[1][m]=180.0;
                    pdih[2][m]=2.0;
                }
            }
            else if(((namj.indexOf("C_")==0) && (namk.indexOf("N_")==0)) ||
            ((namj.indexOf("N_")==0) && (namk.indexOf("C_")==0))) {
                if((nami.charAt(2)=='2') || (nami.charAt(2)=='R')) {
                    kdih[m]="cos ";
                    pdih[0][m]=0.5*1.0/degen;
                    pdih[1][m]=-180.0;
                    pdih[2][m]=6.0;
                }
                else {
                    kdih[m]="cos ";
                    pdih[0][m]=0.5*2.0/degen;
                    pdih[1][m]=0.0;
                    pdih[2][m]=3.0;
                }
            }
            else if((namj.indexOf("N_")==0) && (namk.indexOf("N_")==0)) {
                if((nami.charAt(2)=='2') || (nami.charAt(2)=='R')) {
                    kdih[m]="cos ";
                    pdih[0][m]=0.5*1.0/degen;
                    pdih[1][m]=-180.0;
                    pdih[2][m]=6.0;
                }
                else {
                    kdih[m]="cos ";
                    pdih[0][m]=0.5*2.0/degen;
                    pdih[1][m]=0.0;
                    pdih[2][m]=3.0;
                }
            }
            else {
                println("Warning - unknown potential in sp2-sp3 dihedral");
                kdih[m]="unkn";
                pdih[0][m]=0.0;
                pdih[1][m]=0.0;
                pdih[2][m]=0.0;
            }
        }
        else {
            if((namj.indexOf("C_")==0) && (namk.indexOf("C_")==0)) {
                if((naml.charAt(2)=='2') || (naml.charAt(2)=='R')) {
                    kdih[m]="cos ";
                    pdih[0][m]=0.5*1.0/degen;
                    pdih[1][m]=-180.0;
                    pdih[2][m]=6.0;
                }
                else {
                    kdih[m]="cos ";
                    pdih[0][m]=0.5*2.0/degen;
                    pdih[1][m]=0.0;
                    pdih[2][m]=3.0;
                }
            }
            else if(((namj.indexOf("O_")==0) && (namk.indexOf("C_")==0)) ||
            ((namj.indexOf("C_")==0) && (namk.indexOf("O_")==0))) {
                if((naml.charAt(2)=='2') || (naml.charAt(2)=='R')) {
                    kdih[m]="cos ";
                    pdih[0][m]=0.5*1.0/degen;
                    pdih[1][m]=-180.0;
                    pdih[2][m]=6.0;
                }
                else {
                    kdih[m]="cos ";
                    pdih[0][m]=0.5*2.0/degen;
                    pdih[1][m]=180.0;
                    pdih[2][m]=2.0;
                }
            }
            else if(((namj.indexOf("N_")==0) && (namk.indexOf("C_")==0)) ||
            ((namj.indexOf("C_")==0) && (namk.indexOf("N_")==0))) {
                if((naml.charAt(2)=='2') || (naml.charAt(2)=='R')) {
                    kdih[m]="cos ";
                    pdih[0][m]=0.5*1.0/degen;
                    pdih[1][m]=-180.0;
                    pdih[2][m]=6.0;
                }
                else {
                    kdih[m]="cos ";
                    pdih[0][m]=0.5*2.0/degen;
                    pdih[1][m]=0.0;
                    pdih[2][m]=3.0;
                }
            }
            else {
                println("Warning - unknown potential in sp2-sp3 dihedral");
                kdih[m]="unkn";
                pdih[0][m]=0.0;
                pdih[1][m]=0.0;
                pdih[2][m]=0.0;
            }
        }
    }

    int vdwass() {
        /*
*********************************************************************

dl-poly/java routine for assigning van der Waals parameters
in the dreiding force field

copyright daresbury laboratory
author w.smith october 2000

*********************************************************************
         */
        double aaa,bbb,ccc,dd0,rr0;
        int mxvdw=(struct.nunq*(struct.nunq+1))/2;
        int ii,jj,mvdw;
        kvdw=new String[mxvdw];
        idvdw=new String[2][mxvdw];
        pvdw=new double[mxpvdw][mxvdw];

        // determine number and identities of pair forces

        mvdw=0;

        for(int i=0;i<struct.nunq;i++) {
            ii=ident[i];
            for(int j=i;j<struct.nunq;j++) {
                jj=ident[j];
                if(mvdw >= mxvdw) {
                    println("Error - too many VDW potentials"+
                    BML.fmt(mvdw,8)+BML.fmt(mxvdw,8));
                    return -1;
                }
                idvdw[0][mvdw]=struct.unqatm[i];
                idvdw[1][mvdw]=struct.unqatm[j];
                if(lbuck) {
                    kvdw[mvdw]="buck";
                    aaa=Math.sqrt(36.0*deps[ii]*deps[jj]*Math.exp(dzet[ii]+dzet[jj])/
                    ((dzet[ii]-6.0)*(dzet[jj]-6.0)));
                    bbb=Math.sqrt(dzet[ii]*dzet[jj]*Math.pow((dsig[ii]*dsig[jj]),6)/
                    ((dzet[ii]-6.0)*(dzet[jj]-6.0)));
                    ccc=0.5*(dzet[ii]/dsig[ii]+dzet[jj]/dsig[jj]);
                    pvdw[0][mvdw]=aaa;
                    pvdw[1][mvdw]=1.0/ccc;
                    pvdw[2][mvdw]=bbb;
                }
                else {
                    kvdw[mvdw]="12-6";
                    dd0=Math.sqrt(deps[ii]*deps[jj]);
                    rr0=0.5*(dsig[ii]+dsig[jj]);
                    pvdw[0][mvdw]=dd0*Math.pow(rr0,12);
                    pvdw[1][mvdw]=2.0*dd0*Math.pow(rr0,6);
                    pvdw[2][mvdw]=0.0;
                }
                mvdw++;
            }
        }
        if(mvdw > 0)println("Van der Waals terms assigned");
        return mvdw;
    }

    int hbdass() {
        /*
*********************************************************************

dl-poly routine for assigning H-bond parameters
in the dreiding force field

copyright daresbury laboratory
author w.smith december 2000

*********************************************************************
         */
        int kdk,ii,mhbd;
        int jdk[]=new int[struct.nunq];

        // determine number and identities of bonds

        kdk=0;

        if(lhbnd>=0) {
            for(int i=0;i<struct.nunq;i++) {
                jdk[i]=0;
            }

            for(int i=0;i<struct.nunq;i++) {
                ii=ident[i];
                if(dsym[ii].indexOf("O_") == 0) {
                    kdk++;
                    jdk[i]=1;
                }
                else if(dsym[ii].indexOf("N_") == 0) {
                    kdk++;
                    jdk[i]=1;
                }
                else if(dsym[ii].indexOf("F_") == 0) {
                    kdk++;
                    jdk[i]=1;
                }
            }
        }

        mhbd=(kdk*(kdk+1))/2;
        khbd=new String[mhbd];
        idhbd=new String[3][mhbd];
        phbd=new double[mxphbd][mhbd];
        mhbd=0;
        for(int i=0;i<struct.nunq;i++) {
            if(jdk[i] > 0) {
                for(int j=i;j<struct.nunq;j++) {
                    if(jdk[j] > 0) {
                        khbd[mhbd]="hbnd";
                        idhbd[0][mhbd]=struct.unqatm[i];
                        idhbd[1][mhbd]=struct.unqatm[lhbnd];
                        idhbd[2][mhbd]=struct.unqatm[j];
                        phbd[0][mhbd]=7.0;
                        phbd[1][mhbd]=2.75;
                        phbd[2][mhbd]=0.0;
                        phbd[3][mhbd]=0.0;
                        phbd[4][mhbd]=7.5;
                        mhbd++;
                    }
                }
            }
        }
        if(mhbd > 0) println("H-bond terms assigned");
        return mhbd;
    }
}

