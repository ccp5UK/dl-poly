import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;

// Define the Graphical User Interface

public class MakePoly extends Basic implements ActionListener {
        /*
*********************************************************************

dl_poly/java GUI class to make amorphous polymer CONFIG files

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
    public static GUI home;
    public static MakePoly job;
    private static double ccbond,chbond,ebond,dihscl,heng,volm,temp,size;
    private static double[] eps,sig,edih;
    private static int ncarbons,natms;
    private static JButton make,close;
    private static JCheckBox pbc;
    private static JTextField nc,svol,stem,ccb,chb,cce,che,hhe,ccs,chs,hhs;
    private static JTextField ccbe,de1,de2,de3,de4,dhs;
    private static boolean lrej,lpbc;
    private static double[] cell;
    private static double[][] xyz;
    private static Element[] atoms;

    public MakePoly() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        super();
        setTitle("Make Polymer");

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

        //        Define the Close button

        close = new JButton("Close");
        close.setBackground(art.butn);
        close.setForeground(art.butf);
        fix(close,grd,gbc,3,0,1,1);

        //        Spacer

        JLabel pass = new JLabel(" ");
        fix(pass,grd,gbc,1,1,3,1);

        //        Number of C atoms

        nc = new JTextField(8);
        nc.setBackground(art.scrn);
        nc.setForeground(art.scrf);
        fix(nc,grd,gbc,0,2,1,1);

        //        Instruction label 1

        JLabel lab1 = new JLabel("Number of C atoms",JLabel.LEFT);
        fix(lab1,grd,gbc,1,2,3,1);

        //        Instruction label 2

        JLabel lab2 = new JLabel("System Volume (A^3) and Temperature (K):",JLabel.LEFT);
        fix(lab2,grd,gbc,0,3,4,1);

        //        System Volume

        svol = new JTextField(8);
        svol.setBackground(art.scrn);
        svol.setForeground(art.scrf);
        fix(svol,grd,gbc,0,4,1,1);

        //       System Temperature

        stem = new JTextField(8);
        stem.setBackground(art.scrn);
        stem.setForeground(art.scrf);
        fix(stem,grd,gbc,2,4,1,1);

        //        Instruction label 3

        JLabel lab3 = new JLabel("C-C and C-H Bondlengths (A):",JLabel.LEFT);
        fix(lab3,grd,gbc,0,5,4,1);

        //        C-C Bondlength

        ccb = new JTextField(8);
        ccb.setBackground(art.scrn);
        ccb.setForeground(art.scrf);
        fix(ccb,grd,gbc,0,6,1,1);

        //        C-H Bondlength

        chb = new JTextField(8);
        chb.setBackground(art.scrn);
        chb.setForeground(art.scrf);
        fix(chb,grd,gbc,2,6,1,1);

        //        Instruction label 4

        JLabel lab4 = new JLabel("C-C, C-H and H-H LJ Epsilon (kJ/mol):",JLabel.LEFT);
        fix(lab4,grd,gbc,0,7,4,1);

        //        C-C Epsilon

        cce = new JTextField(8);
        cce.setBackground(art.scrn);
        cce.setForeground(art.scrf);
        fix(cce,grd,gbc,0,8,1,1);

        //        C-H Epsilon

        che = new JTextField(8);
        che.setBackground(art.scrn);
        che.setForeground(art.scrf);
        fix(che,grd,gbc,1,8,1,1);

        //        H-H Epsilon

        hhe = new JTextField(8);
        hhe.setBackground(art.scrn);
        hhe.setForeground(art.scrf);
        fix(hhe,grd,gbc,2,8,1,1);

        //        Instruction label 5

        JLabel lab5 = new JLabel("C-C, C-H and H-H LJ Sigma (A):",JLabel.LEFT);
        fix(lab5,grd,gbc,0,9,4,1);

        //        C-C Sigma

        ccs = new JTextField(5);
        ccs.setBackground(art.scrn);
        ccs.setForeground(art.scrf);
        fix(ccs,grd,gbc,0,10,1,1);

        //        C-H Sigma

        chs = new JTextField(5);
        chs.setBackground(art.scrn);
        chs.setForeground(art.scrf);
        fix(chs,grd,gbc,1,10,1,1);

        //        H-H Sigma

        hhs = new JTextField(5);
        hhs.setBackground(art.scrn);
        hhs.setForeground(art.scrf);
        fix(hhs,grd,gbc,2,10,1,1);

        //        C-C Bond Energy

        ccbe = new JTextField(5);
        ccbe.setBackground(art.scrn);
        ccbe.setForeground(art.scrf);
        fix(ccbe,grd,gbc,0,11,1,1);

        //        Instruction label 6

        JLabel lab6 = new JLabel("C-C Bond Energy (kJ/mol)",JLabel.LEFT);
        fix(lab6,grd,gbc,1,11,3,1);

        //        Instruction label 7

        JLabel lab7 = new JLabel("Dihedral Energy Constants (kJ/mol):",JLabel.LEFT);
        fix(lab7,grd,gbc,0,12,4,1);

        //        Dihedral energy constant 1

        de1 = new JTextField(8);
        de1.setBackground(art.scrn);
        de1.setForeground(art.scrf);
        fix(de1,grd,gbc,0,13,1,1);

        //        Dihedral energy constant 2

        de2 = new JTextField(8);
        de2.setBackground(art.scrn);
        de2.setForeground(art.scrf);
        fix(de2,grd,gbc,1,13,1,1);

        //        Dihedral energy constant 3

        de3 = new JTextField(8);
        de3.setBackground(art.scrn);
        de3.setForeground(art.scrf);
        fix(de3,grd,gbc,2,13,1,1);

        //        Dihedral energy constant 4

        de4 = new JTextField(8);
        de4.setBackground(art.scrn);
        de4.setForeground(art.scrf);
        fix(de4,grd,gbc,3,13,1,1);

        //        Dihedral LJ Scaling Factor

        dhs = new JTextField(8);
        dhs.setBackground(art.scrn);
        dhs.setForeground(art.scrf);
        fix(dhs,grd,gbc,0,14,1,1);

        //        Instruction label 9

        JLabel lab9 = new JLabel("Dihedral LJ Scaling Factor",JLabel.LEFT);
        fix(lab9,grd,gbc,1,14,3,1);

        //        Period boundary option

        pbc = new JCheckBox("PBC",true);
        pbc.setBackground(art.back);
        pbc.setForeground(art.fore);
        fix(pbc,grd,gbc,3,15,1,1);

        // Register action buttons

        make.addActionListener(this);
        close.addActionListener(this);

    }

    // Constructor method

    public MakePoly(GUI here) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */

        println("Activated panel for making amorphous polymer CONFIG files");
        home=here;

        // Define arrays

        eps = new double[3];
        sig = new double[3];
        edih = new double[4];

        // Set up Graphical User interface

        job = new MakePoly();
        job.pack();
        job.setVisible(true);
        setValues();

    }

    // Set initial values

    static void setValues() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        ncarbons=1;
        volm=1000.0;
        temp=100.0;
        ccbond=cc1b;
        chbond=ch1b;
        eps[0]=0.2977;
        eps[1]=0.1573;
        eps[2]=0.0831;
        sig[0]=3.4;
        sig[1]=3.1;
        sig[2]=2.8;
        ebond=348.0;
        edih[0]=8.832;
        edih[1]=18.087;
        edih[2]=4.88;
        edih[3]=-31.8;
        dihscl=0.5;
        nc.setText(String.valueOf(ncarbons));
        svol.setText(String.valueOf(volm));
        stem.setText(String.valueOf(temp));
        ccb.setText(String.valueOf(ccbond));
        chb.setText(String.valueOf(chbond));
        cce.setText(String.valueOf(eps[0]));
        che.setText(String.valueOf(eps[1]));
        hhe.setText(String.valueOf(eps[2]));
        ccs.setText(String.valueOf(sig[0]));
        chs.setText(String.valueOf(sig[1]));
        hhs.setText(String.valueOf(sig[2]));
        ccbe.setText(String.valueOf(ebond));
        de1.setText(String.valueOf(edih[0]));
        de2.setText(String.valueOf(edih[1]));
        de3.setText(String.valueOf(edih[2]));
        de4.setText(String.valueOf(edih[3]));
        dhs.setText(String.valueOf(dihscl));

    }
    static int polyGrow() {
        /*
*********************************************************************

dl_poly/java utility for generating an amorphous polymer chain

copyright daresbury laboratory
author    w.smith november 2000

*********************************************************************
         */
        int mxtry=25,mxfail=1000;
        int numc,nstep,lhist,nfail,call,imcon;
        boolean opr;
        double[] xbs,ybs,zbs,rot,uuu,vvv;
        double rho,xx0,yy0,zz0,alp,bet,gam;

        opr=false;
        int it[]=new int[2];
        uuu=new double[3];
        vvv=new double[3];
        xbs=new double[4];
        ybs=new double[4];
        zbs=new double[4];
        rot=new double[9];
        natms=3*ncarbons+2;
        cell=new double[9];
        atoms=new Element[natms];
        xyz=new double[3][natms];

        //     calculate mass density (amu/A^3)

        rho=(14.0268*ncarbons+2.0158)/volm;

        //     linear width of cell

        size=Math.pow(volm,(1.0/3.0));

        //     define tetrahedral groups

        xbs[0]=0.0;
        ybs[0]=0.0;
        zbs[0]=1.0;

        xbs[1]=2.0*Math.sqrt(2.0)/3.0;
        ybs[1]=0.0;
        zbs[1]=-1.0/3.0;

        xbs[2]=-Math.sqrt(2.0)/3.0;
        ybs[2]=-Math.sqrt(2.0/3.0);
        zbs[2]=-1.0/3.0;

        xbs[3]=-Math.sqrt(2.0)/3.0;
        ybs[3]=Math.sqrt(2.0/3.0);
        zbs[3]=-1.0/3.0;

        //     position first tetrahedron

        xx0=0.0;
        yy0=0.0;
        zz0=0.0;
        natms=1;
        alp=2.0*Math.PI*Math.random();
        bet=Math.PI*Math.random();
        gam=2.0*Math.PI*Math.random();
        AML.euler(alp,bet,gam,rot);
        uuu[0]=xbs[0];
        uuu[1]=ybs[0];
        uuu[2]=zbs[0];
        AML.rotate(opr,uuu,vvv,rot);
        xyz[0][0]=chbond*vvv[0]+xx0;
        xyz[1][0]=chbond*vvv[1]+yy0;
        xyz[2][0]=chbond*vvv[2]+zz0;

        numc=0;
        nstep=0;
        while(numc<ncarbons) {
            numc++;
            lhist=0;
            nfail=0;
            lrej=true;
            while(lrej) {
                lrej=false;
                lhist++;
                if(lhist>mxtry) {
                    lhist=0;
                    nfail++;
                    if(nfail>mxfail) {
                        println("Error - too many trial failures");
                        println("Number of successful numcs :"+BML.fmt(numc,8));
                        println("polyGrow terminated unsuccessfully");
                        return -2;
                    }
                    natms=Math.max(4,natms-3*(int)Math.round(3.0*Math.random()));
                    numc=(natms-1)/3+1;
                    xx0=xyz[0][natms];
                    yy0=xyz[1][natms];
                    zz0=xyz[2][natms];
                    vvv[0]=xyz[0][natms-3]-xx0;
                    vvv[1]=xyz[1][natms-3]-yy0;
                    vvv[2]=xyz[2][natms-3]-zz0;
                    vvv[0]=vvv[0]-size*Math.round(vvv[0]/size);
                    vvv[1]=vvv[1]-size*Math.round(vvv[1]/size);
                    vvv[2]=vvv[2]-size*Math.round(vvv[2]/size);
                    alp=Math.atan2(vvv[1],vvv[0]);
                    bet=Math.acos(vvv[2]/Math.sqrt(Math.pow(vvv[0],2)+Math.pow(vvv[1],2)+
                    Math.pow(vvv[2],2)));
                    gam=2.0*Math.PI*Math.random();
                }
                AML.euler(alp,bet,gam,rot);

                xyz[0][natms]=xx0;
                xyz[1][natms]=yy0;
                xyz[2][natms]=zz0;
                uuu[0]=xbs[1];
                uuu[1]=ybs[1];
                uuu[2]=zbs[1];
                AML.rotate(opr,uuu,vvv,rot);
                xyz[0][natms+1]=chbond*vvv[0]+xx0;
                xyz[1][natms+1]=chbond*vvv[1]+yy0;
                xyz[2][natms+1]=chbond*vvv[2]+zz0;
                uuu[0]=xbs[2];
                uuu[1]=ybs[2];
                uuu[2]=zbs[2];
                AML.rotate(opr,uuu,vvv,rot);
                xyz[0][natms+2]=chbond*vvv[0]+xx0;
                xyz[1][natms+2]=chbond*vvv[1]+yy0;
                xyz[2][natms+2]=chbond*vvv[2]+zz0;
                uuu[0]=xbs[3];
                uuu[1]=ybs[3];
                uuu[2]=zbs[3];
                AML.rotate(opr,uuu,vvv,rot);
                xyz[0][natms+3]=chbond*vvv[0]+xx0;
                xyz[1][natms+3]=chbond*vvv[1]+yy0;
                xyz[2][natms+3]=chbond*vvv[2]+zz0;

                //      check if added unit is energetically acceptable

                nstep++;
                if(numc>1)lrej=select();
                gam=2.0*Math.PI*Math.random();

            }

            //     growth direction vector and rotation matrix

            xx0=ccbond*vvv[0]+xx0;
            yy0=ccbond*vvv[1]+yy0;
            zz0=ccbond*vvv[2]+zz0;
            alp=Math.atan2(-vvv[1],-vvv[0]);
            bet=Math.acos(-vvv[2]);

            natms+=3;

        }
        if(numc<ncarbons) {
            println("Error - unable to propagate growth");
            println("Number of successful links :"+BML.fmt(numc,8));
            println("polyGrow terminated unsuccessfully");
            return-3;
        }
        natms++;

        //     determine simulation cell dimensions

        for (int i=0;i<9;i++) {
            cell[i]=0.0;
        }
        cell[0]=size;
        cell[4]=size;
        cell[8]=size;

        //     set atomic numbers and names

        atoms[0]=new Element("H_      ");
        atoms[natms-1]=new Element("H_      ");
        for(int i=0;i<natms-2;i+=3) {
            atoms[i+1]=new Element("C_3     ");
            atoms[i+2]=new Element("H_      ");
            atoms[i+3]=new Element("H_      ");
        }
        if(lpbc) {
            for(int i=0;i<natms;i++) {
                xyz[0][i]=xyz[0][i]-size*Math.round(xyz[0][i]/size);
                xyz[1][i]=xyz[1][i]-size*Math.round(xyz[1][i]/size);
                xyz[2][i]=xyz[2][i]-size*Math.round(xyz[2][i]/size);
            }
        }
        println("polyGrow terminated successfully");
        println("Number of trial moves:"+BML.fmt(nstep,8));
        println("Number of atoms generated:"+BML.fmt(natms,8));

        imcon=0;
        if(lpbc)imcon=1;

        // Create Config object

        config =new Config();
        config.natms=natms;
        config.pbc.imcon=imcon;
        config.atoms=atoms;
        config.pbc.cell=cell;
        config.xyz=xyz;
        config.title="Amorphous polymer with"+BML.fmt(ncarbons,4)+" units";
        config.pbc.buildBoundary(config.pbc.imcon);
        config.structure=new Structure(config);

        // write CONFIG file

        fname="CFGPOL."+String.valueOf(numpol);
        if(!config.configWrite(fname)) return -4;
        numpol++;

        // Draw structure

        if(!editor.isVisible())
            editor.showEditor();
        editor.pane.restore();

        return 0;
    }
    static boolean select() {
        /*
*********************************************************************

dl_poly routine to select or reject a new CH2 monomer
added to a chain using the boltzman factor as the
selection criterion

copyright daresbury laboratory
author w.smith march 1995

*********************************************************************
         */

        int i,j,k;

        double b1x,b1y,b1z,b2x,b2y,b2z,b3x,b3y,b3z,b12x,b12y,b12z,rat;
        double b23x,b23y,b23z,bb23,cost,ddx,ddy,ddz,fact,teng,eterm;
        double boltz,eng,bb12;


        heng=0.0;
        lrej=false;
        boltz=1.0/(temp*RGAS);
        eng=edih[0]-ebond-heng;

        //     bond vectors

        if(natms==4) {
            b1x=xyz[0][1]-xyz[0][0];
            b1y=xyz[1][1]-xyz[1][0];
            b1z=xyz[2][1]-xyz[2][0];
            b1x=b1x-size*Math.round(b1x/size);
            b1y=b1y-size*Math.round(b1y/size);
            b1z=b1z-size*Math.round(b1z/size);
        }
        else {
            b1x=xyz[0][natms-3]-xyz[0][natms-6];
            b1y=xyz[1][natms-3]-xyz[1][natms-6];
            b1z=xyz[2][natms-3]-xyz[2][natms-6];
            b1x=b1x-size*Math.round(b1x/size);
            b1y=b1y-size*Math.round(b1y/size);
            b1z=b1z-size*Math.round(b1z/size);
        }

        b2x=xyz[0][natms]-xyz[0][natms-3];
        b2y=xyz[1][natms]-xyz[1][natms-3];
        b2z=xyz[2][natms]-xyz[2][natms-3];
        b2x=b2x-size*Math.round(b2x/size);
        b2y=b2y-size*Math.round(b2y/size);
        b2z=b2z-size*Math.round(b2z/size);

        b3x=xyz[0][natms+3]-xyz[0][natms];
        b3y=xyz[1][natms+3]-xyz[1][natms];
        b3z=xyz[2][natms+3]-xyz[2][natms];
        b3x=b3x-size*Math.round(b3x/size);
        b3y=b3y-size*Math.round(b3y/size);
        b3z=b3z-size*Math.round(b3z/size);

        //     calculate dihedral angle energy

        b12x=b1y*b2z-b2y*b1z;
        b12y=b1z*b2x-b2z*b1x;
        b12z=b1x*b2y-b2x*b1y;
        bb12=Math.sqrt(Math.pow(b12x,2)+Math.pow(b12y,2)+Math.pow(b12z,2));

        b23x=b2y*b3z-b3y*b2z;
        b23y=b2z*b3x-b3z*b2x;
        b23z=b2x*b3y-b3x*b2y;
        bb23=Math.sqrt(Math.pow(b23x,2)+Math.pow(b23y,2)+Math.pow(b23z,2));

        cost=(b12x*b23x+b12y*b23y+b12z*b23z)/(bb12*bb23);
        eng=eng+edih[0]+cost*(edih[1]+cost*(edih[2]+cost*edih[3]));

        //     C-C backbone LJ interactions

        i=natms;
        for(j=1;j<natms-8;j+=3) {
            ddx=xyz[0][i]-xyz[0][j];
            ddy=xyz[1][i]-xyz[1][j];
            ddz=xyz[2][i]-xyz[2][j];
            ddx=ddx-size*Math.round(ddx/size);
            ddy=ddy-size*Math.round(ddy/size);
            ddz=ddz-size*Math.round(ddz/size);
            rat=Math.pow((Math.pow(sig[0],2)/(ddx*ddx+ddy*ddy+ddz*ddz)),3);
            eng=eng+4.0*eps[0]*(Math.pow(rat,2)-rat);
        }

        //     C-H interactions

        k=1;
        for(j=0;j<natms-3;j++) {
            if(j==k) {
                k+=3;
            }
            else {
                fact=1.0;
                if((j>i-6)||(i==7))fact=dihscl;
                ddx=xyz[0][i]-xyz[0][j];
                ddy=xyz[1][i]-xyz[1][j];
                ddz=xyz[2][i]-xyz[2][j];
                ddx=ddx-size*Math.round(ddx/size);
                ddy=ddy-size*Math.round(ddy/size);
                ddz=ddz-size*Math.round(ddz/size);
                rat=Math.pow((Math.pow(sig[1],2)/(ddx*ddx+ddy*ddy+ddz*ddz)),3);
                eng=eng+4.0*fact*eps[1]*(Math.pow(rat,2)-rat);
            }
        }

        //     initialise energy for replacement Hydrogen atom

        teng=0.0;

        //     H-H 1-4 LJ interactions

        for(i=natms+1;i<natms+4;i++) {
            k=1;
            for(j=0;j<natms;j++) {
                if(j==k) {
                    k+=3;
                }
                else {
                    fact=1.0;
                    if((j>=natms-2)||(natms==4))fact=dihscl;
                    ddx=xyz[0][i]-xyz[0][j];
                    ddy=xyz[1][i]-xyz[1][j];
                    ddz=xyz[2][i]-xyz[2][j];
                    ddx=ddx-size*Math.round(ddx/size);
                    ddy=ddy-size*Math.round(ddy/size);
                    ddz=ddz-size*Math.round(ddz/size);
                    rat=Math.pow((Math.pow(sig[2],2)/(ddx*ddx+ddy*ddy+ddz*ddz)),3);
                    eterm=4.0*fact*eps[2]*(Math.pow(rat,2)-rat);
                    eng=eng+eterm;
                    if(i==natms+3)teng=teng+eterm;
                }
            }
        }

        //     apply selection criterion

        if(eng>0.0) {
            if(Math.random()>Math.exp(-eng*boltz))lrej=true;
        }
        if(!lrej)heng=teng;

        return lrej;
    }

    // Interpret the button and textfield actions

    public void actionPerformed(ActionEvent e) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        int status;
        String arg = (String)e.getActionCommand();
        if (arg.equals("Make")) {
            ncarbons=BML.giveInteger(nc.getText(),1);
            volm=BML.giveDouble(svol.getText(),1);
            temp=BML.giveDouble(stem.getText(),1);
            ccbond=BML.giveDouble(ccb.getText(),1);
            chbond=BML.giveDouble(chb.getText(),1);
            eps[0]=BML.giveDouble(cce.getText(),1);
            eps[1]=BML.giveDouble(che.getText(),1);
            eps[2]=BML.giveDouble(hhe.getText(),1);
            sig[0]=BML.giveDouble(ccs.getText(),1);
            sig[1]=BML.giveDouble(chs.getText(),1);
            sig[2]=BML.giveDouble(hhs.getText(),1);
            ebond=BML.giveDouble(ccbe.getText(),1);
            edih[0]=BML.giveDouble(de1.getText(),1);
            edih[1]=BML.giveDouble(de2.getText(),1);
            edih[2]=BML.giveDouble(de3.getText(),1);
            edih[3]=BML.giveDouble(de4.getText(),1);
            dihscl=BML.giveDouble(dhs.getText(),1);
            lpbc=pbc.isSelected();
            status=polyGrow();
        }
        else if (arg.equals("Close")) {
            job.dispose();
        }
    }
}
