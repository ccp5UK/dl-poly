import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;

public class MakeOPLSField extends Basic implements ActionListener {
    /*
*********************************************************************

dl_poly/java GUI class to make OPLS Field files

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
     */
    public static MakeOPLSField job;
    private static GUI home;
    private static JButton make,load,close;
    private static boolean lshow,lgetf;
    private static String[] dsym;
    private static String[][] idvdw,kbnd,kang,kdih;
    private static double[] dmss,dchg,deps,dsig,pbnd,ppbnd;
    private static double[][] pang,pdih,pvdw,ppang,ppdih,ppvdw;
    private static int mxvdw,numh2o;
    private static int nvdw,ntpatm,ntpbnd,ntpang,ntpdih;
    private static int[] ident,jdent;
    private static Structure struct;

    // Define the Graphical User Interface

    public MakeOPLSField() {
        /*
*********************************************************************
dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
         */
        setTitle("OPLS FIELD Maker");

        getContentPane().setBackground(art.back);
        getContentPane().setForeground(art.fore);
        setFont(fontMain);
        GridBagLayout grd = new GridBagLayout();
        GridBagConstraints gbc = new GridBagConstraints();
        getContentPane().setLayout(grd);

        gbc.fill=GridBagConstraints.BOTH;

        // Pad out

        JLabel pad1 = new JLabel(" ");
        fix(pad1,grd,gbc,0,0,1,1);

        // Define the Make button

        make = new JButton("Make");
        make.setBackground(art.butn);
        make.setForeground(art.butf);
        fix(make,grd,gbc,0,1,1,1);

        // Pad out

        JLabel pad2 = new JLabel(" ");
        fix(pad2,grd,gbc,0,2,1,1);

        // Define the Load button

        load = new JButton("Load");
        load.setBackground(art.butn);
        load.setForeground(art.butf);
        fix(load,grd,gbc,0,3,1,1);

        // Pad out

        JLabel pad3 = new JLabel(" ");
        fix(pad3,grd,gbc,0,4,1,1);

        // Define the Close button

        close = new JButton("Close");
        close.setBackground(art.butn);
        close.setForeground(art.butf);
        fix(close,grd,gbc,0,5,1,1);

        // Pad out

        JLabel pad4 = new JLabel(" ");
        fix(pad4,grd,gbc,0,6,1,1);

        // Register action buttons

        make.addActionListener(this);
        load.addActionListener(this);
        close.addActionListener(this);

    }

    public MakeOPLSField(GUI here) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
         */
        // initial values

        println("Activated panel for making OPLS FIELD files");
        home=here;
        lgetf=false;

        // instantiate application

        job=new MakeOPLSField();
        job.pack();
        job.setVisible(true);

    }
    public void actionPerformed(ActionEvent e) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
         */
        int call,key;
        String arg = (String)e.getActionCommand();
        if (arg.equals("Make")) {
            lshow=false;
            lgetf=false;
            call=OPLSField();
        }
        else if (arg.equals("Load")) {
            lgetf=true;
            lshow=true;
            call=OPLSField();
        }
        else if (arg.trim().equals("Close")) {
            job.setVisible(false);
        }
    }
    int OPLSField() {
        /*
*********************************************************************

dl_poly/java routine to read a DL_POLY CONFIG file
and construct a  corresponding  DL_POLY FIELD file
using the OPLS force model

copyright - daresbury laboratory
author    - w.smith january 2002

*********************************************************************
         */
        String fldfile,cfgfile;
        Molecule mol;

        // import OPLS field parameters

        if(OPLS() != 0) {
            println("Error - failed in processing MINIOPLS file");
            return -1;
        }

        // read the CONFIG file

        if(lgetf) config=getConfig(home,"CFG");
        if(config==null || config.natms==0)return 0;
        if(config.structure.nmols <= 0) return -1;
        struct=config.structure;

        // remove redundant H atoms for OPLS

        removeHyd();

        // convert the configuration to OPLS atoms

        fixOPLSIdentity();

        // remake water molecules in OPLS style

        remakeWater();

        // define system structure

        config.structure=new Structure(config);
	struct=config.structure;
        if(struct.nmols <= 0) return -2;

        // write a new CONFIG file

        cfgfile="CFGOPL."+numopl;
        config.configWrite(cfgfile);

        // attach OPLS identity to dl_poly atoms

        if(attach() < 0) return -3;

        // display the new CONFIG file

        if(lshow){
	    if(home.editor != null)
		home.editor.job.setVisible(false);
	    home.editor=new Editor(home);
	    lshow=false;
	}

        // construct the FIELD file

        fldfile="FLDOPL."+numopl;
        numopl++;

        OUT:
            try {
                DataOutputStream outStream = new DataOutputStream(new FileOutputStream(fldfile));

                //	FIELD file header

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

                    // assign OPLS bond parameters

                    if(mol.atom[0].zsym.indexOf("OW") == 0) {
                        outStream.writeBytes("RIGID  1\n");
                        outStream.writeBytes("    4    1    2    3    4\n");
                    }
                    else {
                        if(mol.nbnd > 0) {
                            bndass(mol);
                            outStream.writeBytes("CONSTRAINTS"+BML.fmt(mol.nbnd,8)+"\n");

                            for(int i=0;i<mol.nbnd;i++) {
                                outStream.writeBytes(BML.fmt(mol.ibnd[0][i]+1,5)+
                                BML.fmt(mol.ibnd[1][i]+1,5)+
                                BML.fmt(ppbnd[i],12)+"\n");
                            }
                        }

                        // assign OPLS angle parameters

                        if(mol.nang > 0) {
                            angass(mol);
                            outStream.writeBytes("ANGLES"+BML.fmt(mol.nang,8)+"\n");

                            for (int i=0;i<mol.nang;i++) {
                                outStream.writeBytes(BML.fmt("HARM",4)+
                                BML.fmt(mol.iang[0][i]+1,5)+
                                BML.fmt(mol.iang[1][i]+1,5)+
                                BML.fmt(mol.iang[2][i]+1,5)+
                                BML.fmt(ppang[0][i],12)+
                                BML.fmt(ppang[1][i],12)+"\n");
                            }
                        }

                        // assign OPLS dihedral angles

                        if(mol.ndih>0) {
                            dihass(mol);
                            outStream.writeBytes("DIHEDRALS"+BML.fmt(mol.ndih,8)+"\n");

                            for (int i=0;i<mol.ndih;i++) {
                                outStream.writeBytes(BML.fmt("OPLS",4)+
                                BML.fmt(mol.idih[0][i]+1,5)+
                                BML.fmt(mol.idih[1][i]+1,5)+
                                BML.fmt(mol.idih[2][i]+1,5)+
                                BML.fmt(mol.idih[3][i]+1,5)+
                                BML.fmt(ppdih[0][i],12)+
                                BML.fmt(ppdih[1][i],12)+
                                BML.fmt(ppdih[2][i],12)+
                                BML.fmt(ppdih[3][i],12)+
                                BML.fmt(ppdih[4][i],12)+"\n");
                            }
                        }
                    }
                    outStream.writeBytes("FINISH"+"\n");
                }

                // assign OPLS van der waals potentials

                nvdw=vdwass();
                outStream.writeBytes("VDW"+BML.fmt(nvdw,8)+"\n");

                for (int i=0;i<nvdw;i++) {
                    outStream.writeBytes(BML.fmt(idvdw[0][i],8)+BML.fmt(idvdw[1][i],8)+
                    BML.fmt("LJ   ",5)+BML.fmt(ppvdw[0][i],12)+
                    BML.fmt(ppvdw[1][i],12)+"\n");
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

    int OPLS() {
        /*
*********************************************************************

dl_poly/java utility for reading the OPLS force field

copyright daresbury laboratory
author w.smith january 2002

*********************************************************************
         */
        // open MINIOPLS file and skip header records

        try {
            String record="";
            InputStream instream = this.getClass().getResourceAsStream("MINIOPLS");
            InputStreamReader isr = new InputStreamReader(instream);
            BufferedReader reader = new BufferedReader(isr);
            println("Reading file: MINIOPLS");
            for(int i=0;i<8;i++)
                record = reader.readLine();

            // read MINIOPLS contents

            while((record=reader.readLine()) != null) {
                if(BML.giveWord(record,1).equals("ATOMS")) {
                    ntpatm=BML.giveInteger(record,2);
                    dsym=new String[ntpatm];
                    dmss=new double[ntpatm];
                    dchg=new double[ntpatm];
                    deps=new double[ntpatm];
                    dsig=new double[ntpatm];
                    for(int i=0;i<ntpatm;i++) {
                        record=reader.readLine();
                        dsym[i]=BML.giveWord(record,1);
                        dmss[i]=BML.giveDouble(record,2);
                        deps[i]=BML.giveDouble(record,3);
                        dsig[i]=BML.giveDouble(record,4);
                        dchg[i]=BML.giveDouble(record,5);
                    }
                }
                else if(BML.giveWord(record,1).equals("CONSTR")) {
                    ntpbnd=BML.giveInteger(record,2);
                    kbnd=new String[2][ntpbnd];
                    pbnd=new double[ntpbnd];
                    for(int i=0;i<ntpbnd;i++) {
                        record=reader.readLine();
                        kbnd[0][i]=BML.giveWord(record,1);
                        kbnd[1][i]=BML.giveWord(record,2);
                        pbnd[i]=BML.giveDouble(record,3);
                    }
                }
                else if(BML.giveWord(record,1).equals("ANGLES")) {
                    ntpang=BML.giveInteger(record,2);
                    kang=new String[3][ntpang];
                    pang=new double[2][ntpang];
                    for(int i=0;i<ntpang;i++) {
                        record=reader.readLine();
                        kang[0][i]=BML.giveWord(record,1);
                        kang[1][i]=BML.giveWord(record,2);
                        kang[2][i]=BML.giveWord(record,3);
                        pang[0][i]=BML.giveDouble(record,4);
                        pang[1][i]=BML.giveDouble(record,5);
                    }
                }
                else if(BML.giveWord(record,1).equals("DIHEDRALS")) {
                    ntpdih=BML.giveInteger(record,2);
                    kdih=new String[4][ntpdih];
                    pdih=new double[5][ntpdih];
                    for(int i=0;i<ntpdih;i++) {
                        record=reader.readLine();
                        kdih[0][i]=BML.giveWord(record,1);
                        kdih[1][i]=BML.giveWord(record,2);
                        kdih[2][i]=BML.giveWord(record,3);
                        kdih[3][i]=BML.giveWord(record,4);
                        pdih[0][i]=BML.giveDouble(record,5);
                        pdih[1][i]=BML.giveDouble(record,6);
                        pdih[2][i]=BML.giveDouble(record,7);
                        pdih[3][i]=BML.giveDouble(record,8);
                        pdih[4][i]=BML.giveDouble(record,9);
                    }
                }
            }
            println("MINIOPLS File successfully loaded");
            reader.close();
        }
        catch(FileNotFoundException e) {
            println("Error - file not found: MINIOPLS");
            return -1;
        }
        catch(Exception e) {
            println("Error reading file: MINIOPLS" + " "+e);
            return -2;
        }
        return 0;
    }

    int attach() {
        /*
*********************************************************************

dl_poly/java utility for attaching the OPLS force field identity
to a dl_poly atom

copyright daresbury laboratory
author w.smith october 2001

*********************************************************************
         */
        int status=0;
        ident=new int[struct.nunq];

        numh2o=0;

        // attach link to dl_poly identity

        for(int i=0;i<struct.nunq;i++) {
            ident[i]=-1;
            for(int j=0;j<ntpatm;j++) {
                if(struct.unqatm[i].trim().equals(dsym[j].trim())) {
                    ident[i]=j;
                }
            }
            if(ident[i]<0) {
                println("Error - unidentified atom: "+struct.unqatm[i]);
                status=-1;
            }
        }
        return status;
    }

    int identify(Molecule mol) {
        /*
*********************************************************************

dl_poly/java utility for attaching the OPLS force field identity
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
            for(int j=0;j<ntpatm;j++) {
                if(mol.atom[i].zsym.trim().equals(dsym[j].trim())) {
                    jdent[i]=j;
                }
            }
            if(jdent[i]<0) {
                println("Error - unidentified atom: "+mol.atom[i].zsym);
                status=-1;
            }
            else {
                // assign OPLS masses and charges

                mol.atom[i].zmas=dmss[jdent[i]];
                mol.chge[i]=dchg[jdent[i]];
            }
        }
        return status;
    }

    void bndass(Molecule mol) {
        /*
*********************************************************************

dl_poly/java routine for assigning bond constraints in the
OPLS force field

copyright daresbury laboratory
author w.smith january 2002

*********************************************************************
         */
        int i,j;
        boolean safe;
        ppbnd=new double[mol.nbnd];

        for (int k=0;k<mol.nbnd;k++) {
            safe=false;
            i=mol.ibnd[0][k];
            j=mol.ibnd[1][k];
            for(int m=0;m<ntpbnd;m++) {
                if((dsym[jdent[i]].indexOf(kbnd[0][m]) == 0 && dsym[jdent[j]].indexOf(kbnd[1][m]) == 0)
                ||(dsym[jdent[i]].indexOf(kbnd[1][m]) == 0 && dsym[jdent[j]].indexOf(kbnd[0][m]) == 0)) {
                    ppbnd[k]=pbnd[m];
                    safe=true;
                }
            }
            if(!safe) {
                println("Warning - bond constraint unassigned: "+
                BML.fmt(mol.atom[mol.ibnd[0][k]].zsym,10)+
                BML.fmt(mol.atom[mol.ibnd[1][k]].zsym,10));
            }
        }
        println("Bond lengths assigned");
    }

    void angass(Molecule mol) {
        /*
*********************************************************************

dl_poly/java routine for assigning valence angle parameters in the
OPLS force field

copyright daresbury laboratory
author w.smith january 2002

*********************************************************************
         */
        int i,j,k;
        boolean safe;
        ppang=new double[2][mol.nang];

        for (int n=0;n<mol.nang;n++) {
            safe=false;
            i=mol.iang[1][n];
            OUT:
                for(int m=0;m<ntpang;m++) {
                    if(dsym[jdent[i]].indexOf(kang[1][m]) == 0) {
                        j=mol.iang[0][n];
                        k=mol.iang[2][n];

                        if((dsym[jdent[j]].indexOf(kang[0][m]) == 0 && dsym[jdent[k]].indexOf(kang[2][m]) == 0)
                        || (dsym[jdent[j]].indexOf(kang[2][m]) == 0 && dsym[jdent[k]].indexOf(kang[0][m]) == 0)) {
                            ppang[0][n]=pang[0][m];
                            ppang[1][n]=pang[1][m];
                            safe=true;
                        }
                    }
                }
                if(!safe) {
                    println("Warning - valence angle unassigned: "+
                    BML.fmt(mol.atom[mol.iang[0][n]].zsym,10)+
                    BML.fmt(mol.atom[mol.iang[1][n]].zsym,10)+
                    BML.fmt(mol.atom[mol.iang[2][n]].zsym,10));
                }
        }
        println("Valence angle potentials assigned");
    }

    void dihass(Molecule mol) {
        /*
*********************************************************************

dl-poly/java routine for assigning dihedral parameters in the
OPLS force field

copyright daresbury laboratory
author w.smith january 2002

*********************************************************************
         */
        String nami,namj,namk,naml;
        int i,j,k,l;
        boolean safe;
        ppdih=new double[5][mol.ndih];

        for(int m=0;m<mol.ndih;m++) {
            safe=false;
            i=mol.idih[0][m];
            j=mol.idih[1][m];
            k=mol.idih[2][m];
            l=mol.idih[3][m];
            nami=mol.atom[i].zsym.trim();
            namj=mol.atom[j].zsym.trim();
            namk=mol.atom[k].zsym.trim();
            naml=mol.atom[l].zsym.trim();

            for(int n=0;n<ntpdih;n++) {
                if((nami.indexOf(kdih[0][n]) == 0 && namj.indexOf(kdih[1][n]) == 0)
                && (namk.indexOf(kdih[2][n]) == 0 && naml.indexOf(kdih[3][n]) == 0)) {
                    ppdih[0][m]=pdih[0][n];
                    ppdih[1][m]=pdih[1][n];
                    ppdih[2][m]=pdih[2][n];
                    ppdih[3][m]=pdih[3][n];
                    ppdih[4][m]=pdih[4][n];
                    safe=true;
                }
                else if((nami.indexOf(kdih[3][n]) == 0 && namj.indexOf(kdih[2][n]) ==0)
                && (namk.indexOf(kdih[1][n]) == 0 && naml.indexOf(kdih[0][n]) == 0)) {
                    ppdih[0][m]=pdih[0][n];
                    ppdih[1][m]=pdih[1][n];
                    ppdih[2][m]=pdih[2][n];
                    ppdih[3][m]=pdih[3][n];
                    ppdih[4][m]=pdih[4][n];
                    safe=true;
                }
            }
            if(!safe) {
                println("Warning - dihedral unassigned: "+BML.fmt(nami,10)
                +BML.fmt(namj,10)+BML.fmt(namk,10)+BML.fmt(naml,10));
            }
        }
        println("Dihedral angles assigned");
    }

    int vdwass() {
        /*
*********************************************************************

dl-poly/java routine for assigning van der Waals parameters
in the OPLS force field

copyright daresbury laboratory
author w.smith january 2002

*********************************************************************
         */
        int mxvdw=(struct.nunq*(struct.nunq+1))/2;
        int ii,jj,mvdw;
        idvdw=new String[2][mxvdw];
        ppvdw=new double[2][mxvdw];

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
                ppvdw[0][mvdw]=Math.sqrt(deps[ii]*deps[jj]);
                ppvdw[1][mvdw]=Math.sqrt(dsig[ii]*dsig[jj]);
                mvdw++;
            }
        }
        if(mvdw > 0)println("Van der Waals terms assigned");
        return mvdw;
    }

    void removeHyd() {
        /*
*********************************************************************

dl-poly/java routine for removing redundant hydrogen atoms
as required by OPLS force field

copyright daresbury laboratory
author w.smith march 2003

*********************************************************************
         */
        int k=0;

        // strip away redundant hydrogen atoms

        for(int i=0;i<config.natms;i++) {
            if(config.atoms[i].znum == 1 &&
            config.atoms[struct.bond[0][i]].znum == 6)
                struct.lbnd[i]=-1;
        }

        for(int i=0;i<config.natms;i++) {
            if(struct.lbnd[i] >= 0) {
                config.atoms[k]=new Element(config.atoms[i].zsym);
                config.xyz[0][k]=config.xyz[0][i];
                config.xyz[1][k]=config.xyz[1][i];
                config.xyz[2][k]=config.xyz[2][i];
                k++;
            }
        }
        config.natms=k;
        struct=new Structure(config);

    }

    void fixOPLSIdentity() {
        /*
*********************************************************************

dl-poly/java routine for identifying and labelling OPLS atom types

copyright daresbury laboratory
author w.smith march 2003

*********************************************************************
         */
        int k,m,c,o,n,h,q;

        // identify and label atom types

        for(int i=0;i<config.natms;i++) {
            c=0;
            o=0;
            n=0;
            h=0;
            q=-1;
            if(struct.lbnd[i] == 1)
                q=struct.bond[0][i];
            else if(struct.lbnd[i] > 1)
                q=i;

            if(q >= 0) {
                for (int j=0;j<struct.lbnd[q];j++) {
                    k=struct.bond[j][q];
                    if(config.atoms[k].zsym.indexOf("C_") >= 0)
                        c++;
                    else if(config.atoms[k].zsym.indexOf("O_") >= 0)
                        o++;
                    else if(config.atoms[k].zsym.indexOf("N_") >= 0)
                        n++;
                    else if(config.atoms[k].zsym.indexOf("H_") >= 0)
                        h++;
                    else if(config.atoms[k].zsym.indexOf("HW") >= 0)
                        h++;
                }
            }
            if(config.atoms[i].zsym.indexOf("C_3") >= 0) {
                if(struct.lbnd[i] == 0) {
                    config.atoms[i].zsym=BML.fmt("C_34",8); // methane
                }
                else if(struct.lbnd[i] == 1) // methyl groups
                {
                    k=struct.bond[0][i];
                    if(struct.lbnd[k] == 1) {
                        if(config.atoms[k].zsym.indexOf("C_3") >= 0) {
                            config.atoms[i].zsym=BML.fmt("C_331",8); // ethane
                        }
                        else {
                            unknown(i);
                        }
                    }
                    else if(struct.lbnd[k] == 2) {
                        if(c == 2 && config.atoms[k].zsym.indexOf("C_3") >= 0) {
                            config.atoms[i].zsym=BML.fmt("C_332",8); // methyl on primary CH2 group
                        }
                        else if(c == 1 && o == 1 && config.atoms[k].zsym.indexOf("C_3") >= 0) {
                            config.atoms[i].zsym=BML.fmt("C_331",8); // methyl on CH(OH) group
                        }
                        else if(c == 1 && h == 1 && config.atoms[k].zsym.indexOf("O_3") >= 0) {
                            config.atoms[i].zsym=BML.fmt("C_33A",8); // methanol
                        }
                        else {
                            unknown(i);
                        }
                    }
                    else if(struct.lbnd[k] == 3) {
                        if(c == 3 && config.atoms[k].zsym.indexOf("C_3") >= 0) {
                            config.atoms[i].zsym=BML.fmt("C_333",8); // methyl on secondary CH group
                        }
                        else if(c == 2 && o == 1 && config.atoms[k].zsym.indexOf("C_3") >= 0) {
                            config.atoms[i].zsym=BML.fmt("C_332",8); // methyl in 2 propanol
                        }
                        else if(c == 3 && config.atoms[k].zsym.indexOf("C_R") >= 0) {
                            config.atoms[i].zsym=BML.fmt("C_333",8); // methyl on aromatic group
                        }
                        else if(config.atoms[k].zsym.indexOf("N_3") >= 0) {
                            if(h == 1 && c == 1) {
                                config.atoms[i].zsym=BML.fmt("C_33N2",8); // methyl on secondary amide
                            }
                            else if(h == 0 && c == 2) {
                                config.atoms[i].zsym=BML.fmt("C_33N3",8); // methyl on tertiary amide
                            }
                        }
                        else if(config.atoms[k].zsym.indexOf("C_2") >= 0) {
                            if(c == 1 && o == 1 && n == 1) {
                                config.atoms[i].zsym=BML.fmt("C_33CON",8); // methyl in acetamide
                            }
                            else {
                                unknown(i);
                            }
                        }
                        else {
                            unknown(i);
                        }
                    }
                    else if(struct.lbnd[k] == 4 && config.atoms[k].zsym.indexOf("C_3") >= 0) {
                        config.atoms[i].zsym=BML.fmt("C_334",8); // methyl on tertiary C group
                    }
                    else {
                        unknown(i);
                    }
                }
                else if(struct.lbnd[i] == 2) // CH2 groups
                {
                    if(c == 2) {
                        config.atoms[i].zsym=BML.fmt("C_32",8); // primary CH2 group
                    }
                    else if(c == 1 && o == 1) {
                        config.atoms[i].zsym=BML.fmt("C_32A",8); // primary CH2 alcohol group
                    }
                    else if(c == 1 && n == 1) {
                        k=struct.bond[0][i];
                        if(config.atoms[struct.bond[1][i]].zsym.indexOf("N_3") >= 0)
                            k=struct.bond[1][i];
                        m=0;
                        for(int j=0;j<struct.lbnd[k];j++) {
                            if(config.atoms[struct.bond[j][k]].zsym.indexOf("C_") >= 0)
                                m++;
                        }
                        if(m == 1) {
                            config.atoms[i].zsym=BML.fmt("C_32N2",8); // CH2 on secondary amide
                        }
                        else if(m == 2) {
                            config.atoms[i].zsym=BML.fmt("C_32N3",8); // CH2 on tertiary amide
                        }
                        else {
                            unknown(i);
                        }
                    }
                    else {
                        unknown(i);
                    }
                }
                else if(struct.lbnd[i] == 3) // CH groups
                {
                    if(c == 3) {
                        config.atoms[i].zsym=BML.fmt("C_31",8); // secondary CH group
                    }
                    else if(c == 2 && o == 1) {
                        config.atoms[i].zsym=BML.fmt("C_31A",8); // secondary CH alcohol group
                    }
                    else if(c == 2 && n == 1) {
                        k=struct.bond[0][i];
                        if(config.atoms[struct.bond[1][i]].zsym.indexOf("N_3") >= 0)
                            k=struct.bond[1][i];
                        else if(config.atoms[struct.bond[2][i]].zsym.indexOf("N_3") >= 0)
                            k=struct.bond[2][i];
                        m=0;
                        for(int j=0;j<struct.lbnd[k];j++) {
                            if(config.atoms[struct.bond[j][k]].zsym.indexOf("C_") >= 0)
                                m++;
                        }
                        if(m == 1) {
                            config.atoms[i].zsym=BML.fmt("C_31N2",8); // alpha CH in alanine
                        }
                        else if(m == 2) {
                            config.atoms[i].zsym=BML.fmt("C_31N3",8); // alpha CH in proline
                        }
                        else {
                            unknown(i);
                        }
                    }
                    else {
                        unknown(i);
                    }
                }
                else if(struct.lbnd[i] == 4) // tertiary C groups
                {
                    if(c == 4) {
                        config.atoms[i].zsym=BML.fmt("C_30",8); // teriary C group
                    }
                    else if(c == 3 && o == 1) {
                        config.atoms[i].zsym=BML.fmt("C_30A",8); // tertiary C alcohol group
                    }
                    else {
                        unknown(i);
                    }
                }
            }
            else if(config.atoms[i].zsym.indexOf("C_2") >= 0) // sp2 carbons
            {
                if(struct.lbnd[i] == 1) {
                    config.atoms[i].zsym=BML.fmt("C_22",8); // terminal CH2 1-ene group
                }
                else if(struct.lbnd[i] == 2) {
                    if(c == 1) {
                        config.atoms[i].zsym=BML.fmt("C_21",8); // CH 2-ene group
                    }
                    else if(o == 1 && n == 1) {
                        config.atoms[i].zsym=BML.fmt("C_21O",8); // CH in formic acid
                    }
                    else {
                        unknown(i);
                    }
                }
                else if(struct.lbnd[i] == 3) {
                    if(c == 2) {
                        config.atoms[i].zsym=BML.fmt("C_20",8); // C iso-ene group
                    }
                    else if(c == 1 &&  o == 2) {
                        config.atoms[i].zsym=BML.fmt("C_20O",8); // C in acid CO
                    }
                    else if(c == 1 &&  o == 1 && n == 1) {
                        config.atoms[i].zsym=BML.fmt("C_20O",8); // C in amide CO
                    }
                    else {
                        unknown(i);
                    }
                }
                else {
                    unknown(i);
                }
            }
            else if(config.atoms[i].zsym.indexOf("C_R") >= 0) {
                if(struct.lbnd[i] == 2) {
                    config.atoms[i].zsym=BML.fmt("C_R1",8); // aromatic CH
                }
                else if(struct.lbnd[i] == 3) {
                    config.atoms[i].zsym=BML.fmt("C_R",8); // aromatic C
                }
                else {
                    unknown(i);
                }
            }
            else if(config.atoms[i].zsym.indexOf("N_3") >= 0) {
                if(h == 0 && c == 3) {
                    config.atoms[i].zsym=BML.fmt("N_30",8); // N in tertiary amide
                }
                else if(h == 1 && c == 2) {
                    config.atoms[i].zsym=BML.fmt("N_31",8); // N in secondary amide
                }
                else if(h == 2 && c == 1) {
                    config.atoms[i].zsym=BML.fmt("N_32",8); // N in primary amide
                }
                else {
                    unknown(i);
                }
            }
            else if(config.atoms[i].zsym.indexOf("O_3") >= 0) {
                if(c == 1 && h == 1) {
                    config.atoms[i].zsym=BML.fmt("O_3A",8); // O in alcohol
                }
                else if(h == 2) {
                    numh2o++;
                    config.atoms[i].zsym=BML.fmt("OW",8); // O in water
                    config.atoms[struct.bond[0][i]].zsym=BML.fmt("HW",8); // H in water
                    config.atoms[struct.bond[1][i]].zsym=BML.fmt("HW",8); // H in water
                }
                else {
                    unknown(i);
                }
            }
            else if(config.atoms[i].zsym.indexOf("O_2") >= 0) {
                if(config.atoms[struct.bond[0][i]].zsym.indexOf("C_2") >= 0) {
                    config.atoms[i].zsym=BML.fmt("O_2C",8);
                }
                else {
                    unknown(i);
                }
            }
            else if(config.atoms[i].zsym.indexOf("OW") >= 0) {
                numh2o++;
                config.atoms[i].zsym=BML.fmt("OW",8); // O in water
                config.atoms[struct.bond[0][i]].zsym=BML.fmt("HW",8); // H in water
                config.atoms[struct.bond[1][i]].zsym=BML.fmt("HW",8); // H in water
            }
            else if(config.atoms[i].zsym.indexOf("H_") >= 0) {
                k=struct.bond[0][i];
                if(config.atoms[k].zsym.indexOf("O_3") >= 0) {
                    if(h == 2) {
                        config.atoms[i].zsym=BML.fmt("HW",8); // H in water
                    }
                    else if(h == 1 && c == 1) {
                        config.atoms[i].zsym=BML.fmt("H_A",8); // H on O in alcohol
                    }
                }
                else if(config.atoms[k].zsym.indexOf("N_3") >= 0) {
                    if(h == 1 && c == 2) {
                        config.atoms[i].zsym=BML.fmt("H_N2",8); // H on N in secondary amide
                    }
                    else if(h == 2 && c == 1) {
                        config.atoms[i].zsym=BML.fmt("H_N1",8); // H on N in primary amide
                    }
                }
                else {
                    unknown(i);
                }
            }
            else if(config.atoms[i].zsym.indexOf("HW") >= 0) {
                config.atoms[i].zsym=BML.fmt("HW",8); // H on water
            }
            else if(config.atoms[i].zsym.indexOf("QW") >= 0) {
                config.atoms[i].zsym=BML.fmt("QW",8); // charge on TIPs4P water
            }
            else {
                unknown(i);
            }
        }

    }
    void unknown(int i) {
        /*
*********************************************************************

dl-poly/java routine to write warning message

copyright daresbury laboratory
author w.smith march 2003

*********************************************************************
         */
	println("Warning - OPLS atom identity not assigned to atom"+BML.fmt(i+1,10)
        +"  "+config.atoms[i].zsym);
    }

    void remakeWater() {
        /*
*********************************************************************

dl-poly/java routine for remaking water molecules in OPLS form
as required by OPLS force field

copyright daresbury laboratory
author w.smith jan 2004

*********************************************************************
         */
        int k=0,i1,i2;

        double aaa[]=new double[3];
        double bbb[]=new double[3];
        double ccc[]=new double[3];
        double rst[]=new double[3];
        double uvw[]=new double[3];

        // create new configuration

        Config cfgnew=new Config();
        cfgnew.natms=config.natms+numh2o;
        cfgnew.title=config.title;
        cfgnew.pbc=config.pbc;
        cfgnew.pbc.cell=config.pbc.cell;
        cfgnew.atoms=new Element[cfgnew.natms];
	cfgnew.xyz=new double[3][cfgnew.natms];

        // identify water molecules

        for(int i=0;i<config.natms;i++) {
            if(config.atoms[i].zsym.equals(BML.fmt("OW",8))) {
                // first O-H bond
                i1=struct.bond[0][i];
                rst[0]=config.xyz[0][i1]-config.xyz[0][i];
                rst[1]=config.xyz[1][i1]-config.xyz[1][i];
                rst[2]=config.xyz[2][i1]-config.xyz[2][i];
                config.pbc.images(rst);
                BML.vnorm(rst);
                // second O-H bond
                i2=struct.bond[1][i];
                uvw[0]=config.xyz[0][i2]-config.xyz[0][i];
                uvw[1]=config.xyz[1][i2]-config.xyz[1][i];
                uvw[2]=config.xyz[2][i2]-config.xyz[2][i];
                config.pbc.images(uvw);
                BML.vnorm(uvw);
                // first basis vector
                ccc=BML.cross(rst,uvw);
                BML.vnorm(ccc);
                // second basis vector
                aaa[0]=rst[0]+uvw[0];
                aaa[1]=rst[1]+uvw[1];
                aaa[2]=rst[2]+uvw[2];
                BML.vnorm(aaa);
                // third basis vector
                bbb=BML.cross(aaa,ccc);
                BML.vnorm(bbb);
                // allocate oxygen atom coordinates
                cfgnew.atoms[k]=config.atoms[i];
                cfgnew.xyz[0][k]=config.xyz[0][i];
                cfgnew.xyz[1][k]=config.xyz[1][i];
                cfgnew.xyz[2][k]=config.xyz[2][i];
                // first hydrogen coordinates
                cfgnew.atoms[k+1]=config.atoms[i1];
                cfgnew.xyz[0][k+1]=config.xyz[0][i]+0.585882*aaa[0]+0.75695*bbb[0];
                cfgnew.xyz[1][k+1]=config.xyz[1][i]+0.585882*aaa[1]+0.75695*bbb[1];
                cfgnew.xyz[2][k+1]=config.xyz[2][i]+0.585882*aaa[2]+0.75695*bbb[2];
                // second hydrogen coordinates
                cfgnew.atoms[k+2]=config.atoms[i2];
                cfgnew.xyz[0][k+2]=config.xyz[0][i]+0.585882*aaa[0]-0.75695*bbb[0];
                cfgnew.xyz[1][k+2]=config.xyz[1][i]+0.585882*aaa[1]-0.75695*bbb[1];
                cfgnew.xyz[2][k+2]=config.xyz[2][i]+0.585882*aaa[2]-0.75695*bbb[2];
                // charge centre
                cfgnew.atoms[k+3]=new Element("QW");
                cfgnew.xyz[0][k+3]=config.xyz[0][i]+0.15*aaa[0];
                cfgnew.xyz[1][k+3]=config.xyz[1][i]+0.15*aaa[1];
                cfgnew.xyz[2][k+3]=config.xyz[2][i]+0.15*aaa[2];
                k+=4;
            }
            else if(!config.atoms[i].zsym.equals(BML.fmt("HW",8))) {
                cfgnew.atoms[k]=config.atoms[i];
                cfgnew.xyz[0][k]=config.xyz[0][i];
                cfgnew.xyz[1][k]=config.xyz[1][i];
                cfgnew.xyz[2][k]=config.xyz[2][i];
                k++;
            }
        }

        config=cfgnew;
        config.structure=new Structure(config);

    }

}

