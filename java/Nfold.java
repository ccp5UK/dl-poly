import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;

// Define the Graphical User Interface

public class Nfold extends Basic implements ActionListener {
    /*
*********************************************************************

dl_poly/java GUI class to multiply CONFIG file contents

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
     */
    public static GUI home=null;
    public static Nfold job;
    private static int nxnum,nynum,nznum;
    private static double zgap;
    private static JCheckBox twin,show;
    private static JButton make,load,close;
    private static JTextField nx,ny,nz,gap;
    private static boolean ltwin,lshow,lgetf;

    public Nfold() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        super();
        setTitle("N-Fold Expansion");
        getContentPane().setBackground(art.back);
        getContentPane().setForeground(art.fore);
        setDefaultCloseOperation(DISPOSE_ON_CLOSE);
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

        // Define the Load button

        load = new JButton("Load");
        load.setBackground(art.butn);
        load.setForeground(art.butf);
        fix(load,grd,gbc,2,0,1,1);

        // Instruction label 1

        JLabel lab1 = new JLabel("Replication in A,B,C directions",JLabel.LEFT);
        fix(lab1,grd,gbc,0,1,3,1);

        // Replication in direction A

        nx = new JTextField(5);
        nx.setBackground(art.scrn);
        nx.setForeground(art.scrf);
        fix(nx,grd,gbc,0,2,1,1);

        // Replication in direction B

        ny = new JTextField(5);
        ny.setBackground(art.scrn);
        ny.setForeground(art.scrf);
        fix(ny,grd,gbc,1,2,1,1);

        // Replication in direction C

        nz = new JTextField(5);
        nz.setBackground(art.scrn);
        nz.setForeground(art.scrf);
        fix(nz,grd,gbc,2,2,1,1);

        // Bilayer option

        twin = new JCheckBox("Z-Bilayer");
        twin.setBackground(art.back);
        twin.setForeground(art.fore);
        fix(twin,grd,gbc,0,3,1,1);

        // Instruction label 2

        JLabel lab2 = new JLabel("Z Gap:",JLabel.RIGHT);
        fix(lab2,grd,gbc,1,3,1,1);

        // Z gap

        gap = new JTextField(8);
        gap.setBackground(art.scrn);
        gap.setForeground(art.scrf);
        fix(gap,grd,gbc,2,3,1,1);

        // Define the Close button

        close = new JButton("Close");
        close.setBackground(art.butn);
        close.setForeground(art.butf);
        fix(close,grd,gbc,0,4,1,1);

        // Drawing option

        show = new JCheckBox("Show CFG.");
        show.setBackground(art.back);
        show.setForeground(art.fore);
        fix(show,grd,gbc,1,4,1,1);

        // Register action buttons

        make.addActionListener(this);
        load.addActionListener(this);
        close.addActionListener(this);

    }

    // Constructor method

    public Nfold(GUI here) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        home=here;
        println("Activated panel for N-fold expansion of a CONFIG file");

        // Set up Graphical User interface

        job = new Nfold();
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
        nxnum=1;
        nynum=1;
        nznum=1;
        zgap=25.0;
        ltwin=false;
        lshow=false;
        twin.setSelected(ltwin);
        show.setSelected(lshow);
        nx.setText(String.valueOf(nxnum));
        ny.setText(String.valueOf(nynum));
        nz.setText(String.valueOf(nznum));
        gap.setText(String.valueOf(zgap));
    }

    public void actionPerformed(ActionEvent e) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        int call;
        String arg = (String)e.getActionCommand();
        if (arg.equals("Make")) {
            lgetf=false;
            ltwin=twin.isSelected();
            lshow=show.isSelected();
            zgap=BML.giveDouble(gap.getText(),1);
            nxnum=BML.giveInteger(nx.getText(),1);
            nynum=BML.giveInteger(ny.getText(),1);
            nznum=BML.giveInteger(nz.getText(),1);
            call=n_fold();
        }
        else if (arg.equals("Load")) {
            lgetf=true;
            ltwin=twin.isSelected();
            lshow=show.isSelected();
            zgap=BML.giveDouble(gap.getText(),1);
            nxnum=BML.giveInteger(nx.getText(),1);
            nynum=BML.giveInteger(ny.getText(),1);
            nznum=BML.giveInteger(nz.getText(),1);
            call=n_fold();
        }
        else if (arg.equals("Close")) {
            job.dispose();
        }
    }
    int n_fold() {
        /*
**********************************************************************

dl_poly/java utility to expand a simulation cell by a
nx*ny*nz multiplication

copyright daresbury laboratory

author w. smith january 2001

**********************************************************************
         */
        Config cfgbig;
        String aname="";
        String record="";
        int kkk,ii,i,j,k,n,npt,newnum,natms,imcon;
        double fx,fy,fz,xs,ys,zs,xt,yt,zt,vxx,vyy,vzz,fxx,fyy,fzz;
        double xx,yy,zz,wdth,rrr,ttt,ca,sa;
        double cell[]=new double[9];
        double rcell[]=new double[9];
        double cprp[];

        // read the CONFIG file

        if(lgetf) config=getConfig(home,"CFG");
        if(config==null || config.natms==0)return -1;
        natms=config.natms;
        imcon=config.pbc.imcon;
        cell=config.pbc.cell;

        if(lgetf) {
            npt=fname.indexOf(".");
            fname="CFGBIG"+fname.substring(npt);
        }
        else
            fname="CFGBIG";

        println("Number of atoms in CONFIG file: "+natms);

        newnum=nxnum*nynum*nznum*natms;
        if(ltwin)newnum*=2;
        Element atoms[]=new Element[newnum];
        double xyz[][]=new double[3][newnum];

        // twin chains in the Z direction

        if(ltwin) {
            if(imcon==0 || imcon==4 || imcon==5) {
                println("Error - bilayer not possible for this PBC");
                return -2;
            }
            rrr=Math.sqrt(Math.pow(cell[6],2)+Math.pow(cell[7],2));

            if(rrr<1.e-8) {
                ca=1.0;
                sa=0.0;
            }
            else {
                ca=cell[6]/rrr;
                sa=cell[7]/rrr;
            }

            rrr=0.5*zgap/cell[8];

            for(i=0;i<natms;i++) {
                atoms[i]=config.atoms[i];
                xyz[0][i]=config.xyz[0][i]+cell[6]*rrr;
                xyz[1][i]=config.xyz[1][i]+cell[7]*rrr;
                xyz[2][i]=config.xyz[2][i]+cell[8]*rrr;

                atoms[i+natms]=config.atoms[i];
                xyz[0][i+natms]=-2.0*ca*sa*xyz[0][i]+(ca*ca-sa*sa)*xyz[1][i];
                xyz[1][i+natms]=(ca*ca-sa*sa)*xyz[0][i]+2.0*ca*sa*xyz[1][i];
                xyz[2][i+natms]=-xyz[2][i];
            }

            natms*=2;
            cell[6]=2.0*(1.0+rrr)*cell[6];
            cell[7]=2.0*(1.0+rrr)*cell[7];
            cell[8]=2.0*(1.0+rrr)*cell[8];

        }
        else {
            for (i=0;i<natms;i++) {
                atoms[i]=config.atoms[i];
                xyz[0][i]=config.xyz[0][i];
                xyz[1][i]=config.xyz[1][i];
                xyz[2][i]=config.xyz[2][i];
            }
        }

        // convert CONFIG coordinates to cell fractional coordinates

        rcell=AML.invert(cell);

        for(n=0;n<natms;n++) {
            xs=xyz[0][n];
            ys=xyz[1][n];
            zs=xyz[2][n];
            xyz[0][n]=rcell[0]*xs+rcell[3]*ys+rcell[6]*zs;
            xyz[1][n]=rcell[1]*xs+rcell[4]*ys+rcell[7]*zs;
            xyz[2][n]=rcell[2]*xs+rcell[5]*ys+rcell[8]*zs;
        }

        // new cell contents

        kkk=natms*nxnum*nynum*nznum-1;
        for(k=nznum-1;k>=0;k--) {
            for(j=nynum-1;j>=0;j--) {
                for(i=nxnum-1;i>=0;i--) {
                    for(n=natms-1;n>=0;n--) {
                        xt=xyz[0][n]+0.5;
                        yt=xyz[1][n]+0.5;
                        zt=xyz[2][n]+0.5;
                        xs=i+xt-0.5*nxnum;
                        ys=j+yt-0.5*nynum;
                        zs=k+zt-0.5*nznum;
                        xyz[0][kkk]=cell[0]*xs+cell[3]*ys+cell[6]*zs;
                        xyz[1][kkk]=cell[1]*xs+cell[4]*ys+cell[7]*zs;
                        xyz[2][kkk]=cell[2]*xs+cell[5]*ys+cell[8]*zs;
                        atoms[kkk]=new Element(atoms[n].zsym);
                        kkk--;
                    }
                }
            }
        }

        // new cell vectors

        cell[0]*=nxnum;
        cell[1]*=nxnum;
        cell[2]*=nxnum;
        cell[3]*=nynum;
        cell[4]*=nynum;
        cell[5]*=nynum;
        cell[6]*=nznum;
        cell[7]*=nznum;
        cell[8]*=nznum;

        // new CONFIG object

        cfgbig=new Config();
        cfgbig.atoms=atoms;
        cfgbig.natms=newnum;
        cfgbig.pbc.imcon=imcon;
        cfgbig.pbc.cell=cell;
        cfgbig.xyz=xyz;
        cfgbig.title=config.title;
        cfgbig.pbc.buildBoundary(cfgbig.pbc.imcon);

        // write new CONFIG file

        if(!cfgbig.configWrite(fname))return -3;
        println("File "+fname+" created");
        cprp=AML.dcell(cell);
        wdth=0.5*BML.min(cprp[6],cprp[7],cprp[8]);
        println("Number of atoms in "+fname+" : "+cfgbig.natms);
        println("Maximum cutoff radius = "+BML.fmt(wdth,10));

        // Draw expanded structure

        if(lshow) {
            config=cfgbig;
            if(!editor.isVisible())
                editor.showEditor();
            editor.pane.restore();
        }

        return 0;
    }
}
