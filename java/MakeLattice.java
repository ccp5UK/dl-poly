import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;

// Define the Graphical User Interface

public class MakeLattice extends Basic implements ActionListener {
        /*
*********************************************************************

dl_poly/java GUI class to construct lattice CONFIG files

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
    public static GUI home=null;
    public static MakeLattice job;
    private static int nbas,nxnum,nynum,nznum,atnum;
    private static int status;
    private static double sxnum,synum,sznum;
    private static String atnam,aname[];
    private static int no[];
    private static double unit[],uvw[][];
    private static JButton make,clear,close,enter;
    private static JTextField ax,ay,az,bx,by,bz,cx,cy,cz,sx,sy,sz;
    private static JTextField nx,ny,nz,atna,atco;
    private static boolean kill=true;

    public MakeLattice() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        super();
        setTitle("Make Lattice");

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

        //        Define the Clear button

        clear = new JButton("Clear");
        clear.setBackground(art.butn);
        clear.setForeground(art.butf);
        fix(clear,grd,gbc,1,0,1,1);

        //        Define the Close button

        close = new JButton("Close");
        close.setBackground(art.butn);
        close.setForeground(art.butf);
        fix(close,grd,gbc,2,0,1,1);

        //        Instruction label 1

        JLabel lab1 = new JLabel("Enter unit cell vectors:",JLabel.LEFT);
        fix(lab1,grd,gbc,0,1,3,1);

        //        Instruction label 2

        JLabel lab2 = new JLabel("A vector:",JLabel.LEFT);
        fix(lab2,grd,gbc,0,2,1,1);

        //        Unit cell component Ax

        ax = new JTextField(8);
        ax.setBackground(art.scrn);
        ax.setForeground(art.scrf);
        fix(ax,grd,gbc,0,3,1,1);

        //        Unit cell component Ay

        ay = new JTextField(8);
        ay.setBackground(art.scrn);
        ay.setForeground(art.scrf);
        fix(ay,grd,gbc,1,3,1,1);

        //        Unit cell component Az

        az = new JTextField(8);
        az.setBackground(art.scrn);
        az.setForeground(art.scrf);
        fix(az,grd,gbc,2,3,1,1);

        //        Instruction label 3

        JLabel lab3 = new JLabel("B vector:",JLabel.LEFT);
        fix(lab3,grd,gbc,0,4,1,1);

        //        Unit cell component Bx

        bx = new JTextField(8);
        bx.setBackground(art.scrn);
        bx.setForeground(art.scrf);
        fix(bx,grd,gbc,0,5,1,1);

        //        Unit cell component By

        by = new JTextField(8);
        by.setBackground(art.scrn);
        by.setForeground(art.scrf);
        fix(by,grd,gbc,1,5,1,1);

        //        Unit cell component Bz

        bz = new JTextField(8);
        bz.setBackground(art.scrn);
        bz.setForeground(art.scrf);
        fix(bz,grd,gbc,2,5,1,1);

        //        Instruction label 4

        JLabel lab4 = new JLabel("C vector:",JLabel.LEFT);
        fix(lab4,grd,gbc,0,6,1,1);

        //        Unit cell component Cx

        cx = new JTextField(8);
        cx.setBackground(art.scrn);
        cx.setForeground(art.scrf);
        fix(cx,grd,gbc,0,7,1,1);

        //        Unit cell component Cy

        cy = new JTextField(8);
        cy.setBackground(art.scrn);
        cy.setForeground(art.scrf);
        fix(cy,grd,gbc,1,7,1,1);

        //        Unit cell component Cz

        cz = new JTextField(8);
        cz.setBackground(art.scrn);
        cz.setForeground(art.scrf);
        fix(cz,grd,gbc,2,7,1,1);

        //        Instruction label 5

        JLabel lab5 = new JLabel("Replication in A B C directions:",JLabel.LEFT);
        fix(lab5,grd,gbc,0,8,3,1);

        //        Replication in direction A

        nx = new JTextField(5);
        nx.setBackground(art.scrn);
        nx.setForeground(art.scrf);
        fix(nx,grd,gbc,0,9,1,1);

        //        Replication in direction B

        ny = new JTextField(5);
        ny.setBackground(art.scrn);
        ny.setForeground(art.scrf);
        fix(ny,grd,gbc,1,9,1,1);

        //        Replication in direction C

        nz = new JTextField(5);
        nz.setBackground(art.scrn);
        nz.setForeground(art.scrf);
        fix(nz,grd,gbc,2,9,1,1);

        //        Instruction label 6

        JLabel lab6 = new JLabel("Enter unit cell contents:",JLabel.LEFT);
        fix(lab6,grd,gbc,0,10,3,1);

        //        Instruction label 7

        JLabel lab7 = new JLabel("Atom name:",JLabel.RIGHT);
        fix(lab7,grd,gbc,0,11,2,1);

        //        Atom name field

        atna = new JTextField(8);
        atna.setBackground(art.scrn);
        atna.setForeground(art.scrf);
        fix(atna,grd,gbc,2,11,1,1);

        //        Instruction label 9

        JLabel lab9 = new JLabel("Fractional coordinates:",JLabel.LEFT);
        fix(lab9,grd,gbc,0,12,3,1);

        //        Fractional coordinate Sx

        sx = new JTextField(8);
        sx.setBackground(art.scrn);
        sx.setForeground(art.scrf);
        fix(sx,grd,gbc,0,13,1,1);

        //        Fractional coordinate Sy

        sy = new JTextField(8);
        sy.setBackground(art.scrn);
        sy.setForeground(art.scrf);
        fix(sy,grd,gbc,1,13,1,1);

        //        Fractional coordinate Sz

        sz = new JTextField(8);
        sz.setBackground(art.scrn);
        sz.setForeground(art.scrf);
        fix(sz,grd,gbc,2,13,1,1);

        //      Enter atomic data button

        enter = new JButton("Enter");
        enter.setBackground(art.butn);
        enter.setForeground(art.butf);
        fix(enter,grd,gbc,0,14,1,1);

        //        Instruction label 10

        JLabel lab10 = new JLabel("Atom count:",JLabel.RIGHT);
        fix(lab10,grd,gbc,1,14,1,1);

        //        Atom count field

        atco = new JTextField(8);
        atco.setBackground(art.back);
        atco.setForeground(art.fore);
        fix(atco,grd,gbc,2,14,1,1);

        // Register action buttons

        make.addActionListener(this);
        close.addActionListener(this);
        enter.addActionListener(this);
        clear.addActionListener(this);

    }

    // Constructor method

    public MakeLattice(GUI here) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */

        println("Activated panel for making lattice CONFIG files");
        kill=false;
        home=here;

        // Define arrays

        no = new int[MAXBAS];
        unit = new double[9];
        aname = new String[MAXBAS];
        uvw = new double[3][MAXBAS];

        // Set up Graphical User interface

        job = new MakeLattice();
        job.pack();
        job.setVisible(true);
        setValues();

    }

    // Set initial values

    void setValues() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        nbas=0;
        nxnum=1;
        nynum=1;
        nznum=1;
        atnum=0;
        atnam="";
        status=0;
        sxnum=0.0;
        synum=0.0;
        sznum=0.0;
        for (int i=0;i<9;i++)
            unit[i]=0.0;
        atco.setText(String.valueOf(nbas));
        ax.setText(String.valueOf(unit[0]));
        ay.setText(String.valueOf(unit[1]));
        az.setText(String.valueOf(unit[2]));
        bx.setText(String.valueOf(unit[3]));
        by.setText(String.valueOf(unit[4]));
        bz.setText(String.valueOf(unit[5]));
        cx.setText(String.valueOf(unit[6]));
        cy.setText(String.valueOf(unit[7]));
        cz.setText(String.valueOf(unit[8]));
        nx.setText(String.valueOf(nxnum));
        ny.setText(String.valueOf(nynum));
        nz.setText(String.valueOf(nznum));
        sx.setText(String.valueOf(sxnum));
        sy.setText(String.valueOf(synum));
        sz.setText(String.valueOf(sznum));
        atna.setText(atnam);
    }

    // Add atom to unit cell

    int addatom() {
        /*
**********************************************************************

dl_poly/java utility to add atoms to unit cell

copyright - daresbury laboratory
author    -  w.smith 2011

*********************************************************************
         */
        boolean safe;
        safe=true;

        if(nbas < MAXBAS) {
            if(atnam.equals("")) {
                println("Error - atom type unrecognised");
                safe=false;
            }
            else {
                atnum=AML.atmnum(atnam);
                safe=true;
            }
            if (safe) {
                nbas++;
                no[nbas-1]=atnum;
                aname[nbas-1]=BML.fmt(atnam,8);
                uvw[0][nbas-1]=sxnum;
                uvw[1][nbas-1]=synum;
                uvw[2][nbas-1]=sznum;
            }
        }
        else {
            println("Error - Maximum number of basis atoms reached");
        }
        return nbas;
    }

    // Generate the super cell

    int genlat() {
        /*
***********************************************************************

dl_poly/java routine to generate a perfect lattice with a
general unit cell

copyright - daresbury laboratory
author    - w.smith 2011

***********************************************************************
         */
        int kkk;
        fname="CFGLAT.";
        double xx,yy,zz,xs,ys,zs,wdth;
        double cprp[];
        int status=0;

        config=new Config();
        config.natms=nbas*nxnum*nynum*nznum;
        config.xyz=new double[3][config.natms];
        config.atoms=new Element[config.natms];
        config.title="Lattice file generated by FILE_MAKER utility";
        config.pbc.imcon=3;
        config.levcfg=0;
        config.pbc.cell[0]=nxnum*unit[0];
        config.pbc.cell[1]=nxnum*unit[1];
        config.pbc.cell[2]=nxnum*unit[2];
        config.pbc.cell[3]=nynum*unit[3];
        config.pbc.cell[4]=nynum*unit[4];
        config.pbc.cell[5]=nynum*unit[5];
        config.pbc.cell[6]=nznum*unit[6];
        config.pbc.cell[7]=nznum*unit[7];
        config.pbc.cell[8]=nznum*unit[8];
        config.pbc.buildBoundary(config.pbc.imcon);

        // set up CONFIG file for writing

        try {
            fname+=String.valueOf(numlat);
            DataOutputStream outStream = new DataOutputStream(new FileOutputStream(fname));
            outStream.writeBytes("Lattice file generated by FILE_MAKER utility\n");
            outStream.writeBytes(BML.fmt(0,10)+BML.fmt(3,10)+"\n");
            outStream.writeBytes(BML.fmt(config.pbc.cell[0],20)+BML.fmt(config.pbc.cell[1],20)+BML.fmt(config.pbc.cell[2],20)+"\n");
            outStream.writeBytes(BML.fmt(config.pbc.cell[3],20)+BML.fmt(config.pbc.cell[4],20)+BML.fmt(config.pbc.cell[5],20)+"\n");
            outStream.writeBytes(BML.fmt(config.pbc.cell[6],20)+BML.fmt(config.pbc.cell[7],20)+BML.fmt(config.pbc.cell[8],20)+"\n");

            // set up lattice

            kkk=0;
            for (int k=0;k<nznum;k++) {
                for (int j=0;j<nynum;j++) {
                    for (int i=0;i<nxnum;i++) {
                        for (int n=0;n<nbas;n++) {
                            xs=i+uvw[0][n]-0.5*nxnum;
                            ys=j+uvw[1][n]-0.5*nynum;
                            zs=k+uvw[2][n]-0.5*nznum;
                            xx=unit[0]*xs+unit[3]*ys+unit[6]*zs;
                            yy=unit[1]*xs+unit[4]*ys+unit[7]*zs;
                            zz=unit[2]*xs+unit[5]*ys+unit[8]*zs;
                            outStream.writeBytes(aname[n]+BML.fmt(kkk+1,10)+BML.fmt(no[n],10)+"\n");
                            outStream.writeBytes(BML.fmt(xx,20)+BML.fmt(yy,20)+BML.fmt(zz,20)+"\n");
                            config.atoms[kkk]=new Element(aname[n]);
                            config.xyz[0][kkk]=xx;
                            config.xyz[1][kkk]=yy;
                            config.xyz[2][kkk]=zz;
                            kkk++;
                        }
                    }
                }
            }
            outStream.close();
        }
        catch(Exception e) {
            println("Error - writing file: "+fname);
            status=-3;
            return status;
        }
        config.structure=new Structure(config);
	cfgsav=copyConfig(config);
        cprp=AML.dcell(unit);
        wdth=0.5*BML.min(nxnum*cprp[7],nynum*cprp[8],nznum*cprp[9]);
        println("New CONFIG file "+fname+" created");
        println("Number of atoms in system = "+BML.fmt(nbas*nxnum*nynum*nznum,10));
        println("Maximum radius of cutoff  = "+BML.fmt(wdth,10));
        numlat++;

        // Draw structure

        if(!editor.isVisible())
            editor.showEditor();
        editor.pane.restore();

        return status;
    }

    // Interpret the button and textfield actions

    public void actionPerformed(ActionEvent e) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        String arg = (String)e.getActionCommand();
        if (arg.equals("Make")) {
            unit[0]=BML.giveDouble(ax.getText(),1);
            unit[1]=BML.giveDouble(ay.getText(),1);
            unit[2]=BML.giveDouble(az.getText(),1);
            unit[3]=BML.giveDouble(bx.getText(),1);
            unit[4]=BML.giveDouble(by.getText(),1);
            unit[5]=BML.giveDouble(bz.getText(),1);
            unit[6]=BML.giveDouble(cx.getText(),1);
            unit[7]=BML.giveDouble(cy.getText(),1);
            unit[8]=BML.giveDouble(cz.getText(),1);
            nxnum=BML.giveInteger(nx.getText(),1);
            nynum=BML.giveInteger(ny.getText(),1);
            nznum=BML.giveInteger(nz.getText(),1);
            status=genlat();
        }
        else if (arg.equals("Enter")) {
            atnam=atna.getText();
            sxnum=BML.giveDouble(sx.getText(),1);
            synum=BML.giveDouble(sy.getText(),1);
            sznum=BML.giveDouble(sz.getText(),1);
            nbas=addatom();
            atnam="";
            sxnum=0.0;
            synum=0.0;
            sznum=0.0;
            atna.setText(atnam);
            atco.setText(String.valueOf(nbas));
            sx.setText(String.valueOf(sxnum));
            sy.setText(String.valueOf(synum));
            sz.setText(String.valueOf(sznum));
        }
        else if (arg.equals("Clear")) {
            setValues();
        }
        else if (arg.equals("Close")) {
            if(kill)
                System.exit(0);
            else
                job.dispose();
        }
    }
}
