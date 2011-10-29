import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;

public class MakeBucky extends Basic implements ActionListener {
        /*
*********************************************************************

dl_poly/java GUI class to construct fullerene CONFIG files

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
    public static GUI home;
    public static MakeBucky job;
    private static double ccbond;
    private static int numx,numz;
    private static JTextField bond,nrx,nrz;
    private static JButton make1,make2,close;

    // Define the Graphical User Interface

    public MakeBucky() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        super();
        setTitle("Make Fullerene");

        getContentPane().setBackground(art.back);
        getContentPane().setForeground(art.fore);
        setDefaultCloseOperation(DISPOSE_ON_CLOSE);
        setFont(fontMain);
        GridBagLayout grd = new GridBagLayout();
        GridBagConstraints gbc = new GridBagConstraints();
        getContentPane().setLayout(grd);

        gbc.fill=GridBagConstraints.BOTH;

        //        Define the Make C60 button

        make1 = new JButton("C60");
        make1.setBackground(art.butn);
        make1.setForeground(art.butf);
        fix(make1,grd,gbc,0,0,1,1);

        //        Define the Make Tube button

        make2 = new JButton("Tube");
        make2.setBackground(art.butn);
        make2.setForeground(art.butf);
        fix(make2,grd,gbc,2,0,1,1);
        fix(new JLabel(" "),grd,gbc,1,1,1,1);

        //        Bond length

        JLabel lab1 = new JLabel("C-C Bond (A):",JLabel.LEFT);
        fix(lab1,grd,gbc,0,2,1,1);
        bond = new JTextField(8);
        bond.setBackground(art.scrn);
        bond.setForeground(art.scrf);
        fix(bond,grd,gbc,2,2,1,1);

        //        Tube size - rings in x direction

        JLabel lab2 = new JLabel("Tube size : rings X Y",JLabel.LEFT);
        fix(lab2,grd,gbc,0,3,2,1);
        nrx = new JTextField(8);
        nrx.setBackground(art.scrn);
        nrx.setForeground(art.scrf);
        fix(nrx,grd,gbc,0,4,1,1);

        //        Tube size - rings in y direction

        nrz = new JTextField(8);
        nrz.setBackground(art.scrn);
        nrz.setForeground(art.scrf);
        fix(nrz,grd,gbc,2,4,1,1);

        //        Define the Close button

        fix(new JLabel("  "),grd,gbc,1,5,1,1);
        close = new JButton("Close");
        close.setBackground(art.butn);
        close.setForeground(art.butf);
        fix(close,grd,gbc,2,6,1,1);

        // Register action buttons

        make1.addActionListener(this);
        make2.addActionListener(this);
        close.addActionListener(this);

    }

    public MakeBucky(GUI here) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        home=here;
        println("Activated panel for making fullerene CONFIG files");
        job=new MakeBucky();
        job.pack();
        job.setVisible(true);
        ccbond=ccab;
        numx=8;
        numz=12;
        bond.setText(String.valueOf(ccbond));
        nrx.setText(String.valueOf(numx));
        nrz.setText(String.valueOf(numz));
    }

    int bucky() {
        /*
*********************************************************************

construction of buckminster fullerene molecule

copyright - daresbury laboratory
author    - w.smith november 2000

*********************************************************************
         */

        boolean same;
        int call,i,j,k,n,m;
        double dr,gam,xhi,xlo,base;
        double[][] o=new double[9][2];

        config=new Config();
        config.natms=60;
        config.xyz=new double[3][config.natms+1];
        config.atoms=new Element[config.natms+1];
        config.title="Buckminster Fullerene C60";
        config.pbc.imcon=2;
        config.levcfg=0;

        xhi=0.0;
        xlo=0.0;
        dr=Math.PI/180.0;
        config.xyz[0][0]=0.0;
        config.xyz[1][0]=0.850650808*ccbond;
        config.xyz[2][0]=2.327438436*ccbond;
        o[0][0]=Math.cos(dr*72.0);
        o[1][0]=Math.sin(dr*72.0);
        o[2][0]=0.0;
        o[3][0]=-Math.sin(dr*72.0);
        o[4][0]=Math.cos(dr*72.0);
        o[5][0]=0.0;
        o[6][0]=0.0;
        o[7][0]=0.0;
        o[8][0]=1.0;
        gam=dr*37.37736815;
        o[0][1]=Math.cos(dr*120.0);
        o[1][1]=Math.sin(dr*120.0)*Math.cos(gam);
        o[2][1]=Math.sin(dr*120.0)*Math.sin(gam);
        o[3][1]=-Math.sin(dr*120.0)*Math.cos(gam);
        o[4][1]=Math.cos(dr*120.0)*Math.pow(Math.cos(gam),2)+Math.pow(Math.sin(gam),2);
        o[5][1]=Math.cos(gam)*Math.sin(gam)*(Math.cos(dr*120.0)-1.0);
        o[6][1]=-Math.sin(dr*120.0)*Math.sin(gam);
        o[7][1]=Math.cos(gam)*Math.sin(gam)*(Math.cos(dr*120.0)-1.0);
        o[8][1]=Math.cos(dr*120.0)*Math.pow(Math.sin(gam),2)+Math.pow(Math.cos(gam),2);
        config.natms=1;
        config.atoms[0]=new Element("C_R     ");
        for(m=0;m<100;m++) {
            if(config.natms<60) {
                for(i=0;i<2;i++) {
                    for (j=0;j<60;j++) {
                        if(j < config.natms) {
                            n=config.natms;
                            config.xyz[0][n]=o[0][i]*config.xyz[0][j]+o[3][i]*config.xyz[1][j]+o[6][i]*config.xyz[2][j];
                            config.xyz[1][n]=o[1][i]*config.xyz[0][j]+o[4][i]*config.xyz[1][j]+o[7][i]*config.xyz[2][j];
                            config.xyz[2][n]=o[2][i]*config.xyz[0][j]+o[5][i]*config.xyz[1][j]+o[8][i]*config.xyz[2][j];
                            same=false;
                            for (k=0;k<config.natms;k++) {
                                if(!same) {
                                    same=true;
                                    if(Math.abs(config.xyz[0][n]-config.xyz[0][k]) > 1.e-4)
                                        same=false;
                                    if(Math.abs(config.xyz[1][n]-config.xyz[1][k]) > 1.e-4)
                                        same=false;
                                    if(Math.abs(config.xyz[2][n]-config.xyz[2][k]) > 1.e-4)
                                        same=false;
                                }
                            }
                            if(!same) {
                                config.natms++;
                                config.atoms[n]=new Element("C_R     ");
                                xhi=Math.max(xhi,config.xyz[0][n]);
                                xlo=Math.min(xlo,config.xyz[0][n]);
                            }
                        }
                    }
                }
            }
        }
        println("Number of atoms created : "+config.natms);

        // define cell vectors

        base=xhi-xlo+2.0*ccbond;
        for(i=0;i<9;i++)
            config.pbc.cell[i]=0.0;
        config.pbc.cell[0]=base;
        config.pbc.cell[4]=base;
        config.pbc.cell[8]=base;
        config.pbc.buildBoundary(config.pbc.imcon);
        config.structure=new Structure(config);

        // write CONFIG file

        fname="CFGBUK."+String.valueOf(numbuk);
        if(!config.configWrite(fname)) return -1;
        numbuk++;

        // Draw structure

        if(!editor.isVisible())
            editor.showEditor();
        editor.pane.restore();

        return 0;
    }

    int tube() {
        /*
*********************************************************************

construction of carbon nanotube

copyright - daresbury laboratory
author    - w.smith november 2000

*********************************************************************
         */

        int n,call,levels;
        double height,alp0,alp2,alp4,rad,xxx,yyy,zzz,ang;

        levels=2*(numz+1);

        config=new Config();
        config.natms=numx*levels;
        config.xyz=new double[3][config.natms];
        config.atoms=new Element[config.natms];
        config.title="Carbon Nanotube "+BML.fmt(numx,5)+"  x"+BML.fmt(numz,5);
        config.pbc.imcon=2;
        config.levcfg=0;

        n=0;
        height=ccbond*(1.5*numz+0.5);
        alp0=2.0*Math.PI/numx;
        alp2=alp0/2.0;
        alp4=alp0/4.0;
        rad=0.25*Math.sqrt(3.0)*ccbond/Math.sin(alp4);

        for(int i=0;i<9;i++)
            config.pbc.cell[i]=0.0;
        config.pbc.cell[0]=4.0*rad;
        config.pbc.cell[4]=4.0*rad;
        config.pbc.cell[8]=height+ccbond;
        config.pbc.buildBoundary(config.pbc.imcon);

        zzz=0.5*(ccbond-height-ccbond);
        ang=0.0;
        for(int k=0;k<levels;k++) {
            for(int i=0;i<numx;i++) {
                config.atoms[n]=new Element("C_R     ");
                config.xyz[0][n]=rad*Math.cos(ang+i*alp0);
                config.xyz[1][n]=rad*Math.sin(ang+i*alp0);
                config.xyz[2][n]=zzz;
                n++;
            }
            if(k%2 == 0) {
                ang+=alp2;
                zzz+=0.5*ccbond;
            }
            else {
                zzz+=ccbond;
            }
        }
        config.natms=n;
        config.structure=new Structure(config);
        println("Number of atoms created : "+config.natms);

        // write CONFIG file

        fname="CFGBUK."+String.valueOf(numbuk);
        if(!config.configWrite(fname)) return -1;
        println("Radius of tube (A) is: "+BML.fmt(rad,20));
        numbuk++;

        // Draw structure

        if(!editor.isVisible())
            editor.showEditor();
        editor.pane.restore();

        return 0;
    }

    public void actionPerformed(ActionEvent e) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        int call;
        String arg = (String)e.getActionCommand();
        if (arg.equals("C60")) {
            ccbond=BML.giveDouble(bond.getText(),1);
            call=bucky();
        }
        else if (arg.equals("Tube")) {
            ccbond=BML.giveDouble(bond.getText(),1);
            numx=BML.giveInteger(nrx.getText(),1);
            numz=BML.giveInteger(nrz.getText(),1);
            call=tube();
        }
        else if (arg.equals("Close")) {
            job.dispose();
        }
    }
}
