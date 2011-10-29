import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;

public class MakeTable extends Basic implements ActionListener {
        /*
*********************************************************************

dl_poly/java GUI class to make a TABLE file

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
    static GUI home;
    static MakeTable job;
    static String title,tname,name1,name2,eunit;
    static double rcut,xdat,ydat,ecnvrt;
    static double[] xx,yy,dd,aa,gg,zz,xi,yi;
    static int npnts,numtab,ngrid,keytab,keyfit,mxpnts;
    static JTextField atom1,atom2,gridp,cutoff,xval,yval,points;
    static JButton make,load,close,enter,clear;
    static JComboBox units,fitopt,special;
    static String[] header;

    // Define the Graphical User Interface

    public MakeTable() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        setTitle("Make TABLE File");

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
        fix(load,grd,gbc,1,0,1,1);

        // Instruction label 1

        JLabel lab1 = new JLabel("Atomic names:",JLabel.LEFT);
        fix(lab1,grd,gbc,0,1,2,1);

        // Name of first atom type

        atom1 = new JTextField(8);
        atom1.setBackground(art.scrn);
        atom1.setForeground(art.scrf);
        fix(atom1,grd,gbc,0,2,1,1);

        // Name of second atom type

        atom2 = new JTextField(8);
        atom2.setBackground(art.scrn);
        atom2.setForeground(art.scrf);
        fix(atom2,grd,gbc,1,2,1,1);

        // Instruction label 2

        JLabel lab2 = new JLabel("Grid points:",JLabel.LEFT);
        fix(lab2,grd,gbc,0,3,1,1);

        // Number of grid points

        gridp = new JTextField(8);
        gridp.setBackground(art.scrn);
        gridp.setForeground(art.scrf);
        fix(gridp,grd,gbc,1,3,1,1);

        // Instruction label 3

        JLabel lab3 = new JLabel("Cut off (A):",JLabel.LEFT);
        fix(lab3,grd,gbc,0,4,1,1);

        // Cutoff radius

        cutoff = new JTextField(8);
        cutoff.setBackground(art.scrn);
        cutoff.setForeground(art.scrf);
        fix(cutoff,grd,gbc,1,4,1,1);

        // Instruction label 4

        JLabel lab4 = new JLabel("Enter data points: r ; V(r)");
        fix(lab4,grd,gbc,0,5,2,1);

        // R coordinate

        xval = new JTextField(8);
        xval.setBackground(art.scrn);
        xval.setForeground(art.scrf);
        fix(xval,grd,gbc,0,6,1,1);

        // Function value

        yval = new JTextField(8);
        yval.setBackground(art.scrn);
        yval.setForeground(art.scrf);
        fix(yval,grd,gbc,1,6,1,1);

        // Define the Enter button

        enter = new JButton("Enter");
        enter.setBackground(art.butn);
        enter.setForeground(art.butf);
        fix(enter,grd,gbc,0,7,1,1);

        // Function value

        points = new JTextField(8);
        points.setBackground(art.back);
        points.setForeground(art.fore);
        fix(points,grd,gbc,1,7,1,1);

        // Instruction label 6

        JLabel lab6 = new JLabel("Energy units:",JLabel.LEFT);
        fix(lab6,grd,gbc,0,8,1,1);

        // Energy units required

        units = new JComboBox();
        units.setBackground(art.scrn);
        units.setForeground(art.scrf);
        units.addItem("DL_POLY");
        units.addItem("E_VOLT");
        units.addItem("K_CAL");
        units.addItem("K_JOULE");
        fix(units,grd,gbc,1,8,1,1);


        // Instruction label 7

        JLabel lab7 = new JLabel("Fitting option:",JLabel.LEFT);
        fix(lab7,grd,gbc,0,9,1,1);

        // Energy units required

        fitopt = new JComboBox();
        fitopt.setBackground(art.scrn);
        fitopt.setForeground(art.scrf);
        fitopt.addItem("SPLINE");
        fitopt.addItem("GAUSSIAN");
        fix(fitopt,grd,gbc,1,9,1,1);

        // Instruction label 8

        JLabel lab8 = new JLabel("Special options:",JLabel.LEFT);
        fix(lab8,grd,gbc,0,10,1,1);

        // Special potential functions

        special = new JComboBox();
        special.setBackground(art.scrn);
        special.setForeground(art.scrf);
        special.addItem("None");
        special.addItem("SiO2");
        special.addItem("AgI");
        fix(special,grd,gbc,1,10,1,1);

        // Define the Clear button

        clear = new JButton("Clear");
        clear.setBackground(art.butn);
        clear.setForeground(art.butf);
        fix(clear,grd,gbc,0,11,1,1);

        // Define the Close button

        close = new JButton("Close");
        close.setBackground(art.butn);
        close.setForeground(art.butf);
        fix(close,grd,gbc,1,11,1,1);

        // Register action buttons

        make.addActionListener(this);
        load.addActionListener(this);
        enter.addActionListener(this);
        clear.addActionListener(this);
        close.addActionListener(this);

    }

    public MakeTable(GUI here) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        home=here;
        println("Activated panel for making TABLE files");
        job=new MakeTable();
        job.pack();
        job.setVisible(true);
        header=new String[4];
        startUp();
    }
    void startUp() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        mxpnts=100;
        xi=new double[mxpnts];
        yi=new double[mxpnts];
        ngrid=500;
        npnts=0;
        rcut=10.0;
        xdat=0.0;
        ydat=0.0;
        keytab=0;
        keyfit=0;
        ecnvrt=1.0;
        name1="Name1";
        name2="Name2";
        eunit="DL_POLY";
        atom1.setText(name1);
        atom2.setText(name2);
        special.setSelectedItem("None");
        units.setSelectedItem(eunit);
        fitopt.setSelectedItem("SPLINE");
        points.setText(String.valueOf(npnts));
        gridp.setText(String.valueOf(ngrid));
        cutoff.setText(String.valueOf(rcut));
        xval.setText(String.valueOf(xdat));
        yval.setText(String.valueOf(ydat));
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
        ngrid=BML.giveInteger(gridp.getText(),1);

        if (arg.equals("Make")) {
            xx=new double[npnts];
            yy=new double[npnts];
            System.arraycopy(xi,0,xx,0,npnts);
            System.arraycopy(yi,0,yy,0,npnts);
            rcut=BML.giveDouble(cutoff.getText(),1);
            name1=BML.fmt(atom1.getText(),8);
            name2=BML.fmt(atom2.getText(),8);
            keytab=special.getSelectedIndex();
            keyfit=fitopt.getSelectedIndex();
            eunit=units.getSelectedItem().toString();
            ecnvrt=chkunits(eunit);
            maketab();
        }
        else if (arg.equals("Load")) {
            if((tname=selectFileNameEnds(home,"XY"))!=null) {
                npnts=getXY(tname,0);
                rcut=BML.giveDouble(cutoff.getText(),1);
                name1=BML.fmt(atom1.getText(),8);
                name2=BML.fmt(atom2.getText(),8);
                keytab=special.getSelectedIndex();
                keyfit=fitopt.getSelectedIndex();
                eunit=units.getSelectedItem().toString();
                ecnvrt=chkunits(eunit);
                maketab();
            }
            else {
                println("File selection cancelled");
            }
        }
        else if (arg.equals("Enter")) {
            xi[npnts]=BML.giveDouble(xval.getText(),1);
            yi[npnts]=BML.giveDouble(yval.getText(),1);
            npnts++;
            points.setText(String.valueOf(npnts));
        }
        else if (arg.equals("Clear")) {
            startUp();
        }
        else if (arg.equals("Close")) {
            job.dispose();
        }
    }

    void maketab() {
        /*
*********************************************************************

dl_poly.java utility to make required force tables

author    -  w.smith january 2001
copyright - daresbury laboratory

*********************************************************************
         */
        int call=0;

        println("Making TABLE file .......");

        if(keytab==0) {
            if(npnts > 0) {
                if(keyfit==0)
                    call=tabspline();
                else if(keyfit==1)
                    call=tabgauss();
            }
            else {
                println("Error - no points to fit!");
            }
        }
        else if(keytab==1) {
            call=potSiO2();
        }
        else if(keytab==2) {
            call=potAgI();
        }
    }

    double chkunits(String eunit) {
        /*
*********************************************************************

dl_poly/java utility to determine required energy units

author    -  w.smith january 2001
copyright - daresbury laboratory

*********************************************************************
         */
        double ecnvrt=1.0;

        if(eunit.equals("DL_POLY"))
            ecnvrt=1.0;
        else if(eunit.equals("E_VOLT"))
            ecnvrt=9648.5308210;
        else if(eunit.equals("K_CAL"))
            ecnvrt=418.40;
        else if(eunit.equals("K_JOULE"))
            ecnvrt=100.0;

        return ecnvrt;
    }

    static int getXY(String fname,int np) {
        /*
*********************************************************************

dl_poly/java routine to read a simple XY file

copyright - daresbury laboratory
author    - w.smith january 2001

*********************************************************************
         */
        int n,k,m,mxpnts;
        String record;

        n=0;
        k=0;
        m=0;
        mxpnts=250;
        xx=new double[mxpnts];
        yy=new double[mxpnts];

        try {
            LineNumberReader lnr = new LineNumberReader(new FileReader(fname));
            println("Reading file: "+fname);
            while((record=lnr.readLine()) != null) {
                if(record.charAt(0) == '#') {
                    header[m]=new String(record.substring(1));
                    m=(m+1)%4;
                }
                else {
                    if(record.charAt(0)== '&') {
                        k++;
                    }
                    else if(k==np) {
                        if(n==mxpnts) {
                            double xt[]=new double[2*mxpnts];
                            double yt[]=new double[2*mxpnts];
                            System.arraycopy(xx,0,xt,0,mxpnts);
                            System.arraycopy(yy,0,yt,0,mxpnts);
                            mxpnts*=2;
                            xx=xt;
                            yy=yt;
                        }
                        xx[n]=BML.giveDouble(record,1);
                        yy[n]=BML.giveDouble(record,2);
                        n++;
                    }
                }
            }
            lnr.close();
        }
        catch(FileNotFoundException e) {
            println("Error - file not found: " + fname);
            return -1;
        }
        catch(Exception e) {
            println("Error reading file: " + fname + " "+e);
            return -2;
        }
        println("Number of points loaded:"+BML.fmt(n,6));
        return n;
    }

    int tabspline() {
        /*
*********************************************************************

dl_poly/java routine for constructing a TABLE file from a small
number of data points sampled from a potential function

author - w.smith january 2001
copyright - daresbury laboratory

*********************************************************************
         */
        int np;
        double rpd,sl0,sln,cc0,ccn,x,di,dj;
        double[] uuu,ggg;

        uuu=new double[4];
        ggg=new double[4];
        zz=new double[npnts];
        aa=new double[npnts];
        dd=new double[npnts];
        gg=new double[npnts];

        // check on input variables

        if(ngrid<6) {
            println("Error - incorrect grid number specified for TABLE file");
            return -1;
        }
        else if(rcut<=0.0) {
            println("Error - incorrect cutoff specified for TABLE file");
            return -2;
        }

        // print current list of points and convert to dl_poly units

        println("Current list of spline points:");

        for(int i=0;i<npnts;i++) {
            println(BML.fmt(i+1,6)+BML.fmt(xx[i],8)+BML.fmt(yy[i],8));
            yy[i]=ecnvrt*yy[i];
        }

        // calculate spline terms

        AML.spline(npnts,xx,yy,zz,aa,dd,gg);

        // grid increment

        rpd=rcut/(ngrid-4);

        // Write TABLE file

        fname="TABSPL."+numtab;
        title=name1+" "+name2+" "+"Spline Potential - DL_POLY units";

        try {
            DataOutputStream outStream = new DataOutputStream(new FileOutputStream(fname));
            outStream.writeBytes(title+"\n");
            outStream.writeBytes(BML.fmt(rpd,15)+BML.fmt(rcut,15)+BML.fmt(ngrid,10)+"\n");
            outStream.writeBytes(BML.fmt(name1,8)+BML.fmt(name2,8)+"\n");
            np=1;
            sl0=(yy[1]-yy[0])/dd[0]-gg[1]*dd[0]/6.0;
            sln=(yy[npnts-1]-yy[npnts-2])/dd[npnts-2]+gg[npnts-2]*dd[npnts-2]/6.0;
            cc0=yy[0]-sl0*xx[0];
            ccn=yy[npnts-1]-sln*xx[npnts-1];
            for(int i=1;i<=ngrid;i+=4) {
                for(int j=0;j<4;j++) {
                    x=rpd*(i+j);
                    if(x<xx[0]) {
                        uuu[j]=sl0*x+cc0;;
                    }
                    else if(x>xx[npnts-1]) {
                        uuu[j]=sln*x+ccn;
                    }
                    else {
                        if(x>xx[np])np++;
                        di=x-xx[np-1];
                        dj=xx[np]-x;
                        uuu[j]=(di*yy[np]+dj*yy[np-1]-di*dj*
                        ((dd[np-1]+dj)*gg[np-1]+(dd[np-1]+di)*gg[np])/6.0)/dd[np-1];
                    }
                }
                outStream.writeBytes(BML.fmt(uuu[0],15)+BML.fmt(uuu[1],15)+BML.fmt(uuu[2],15)
                +BML.fmt(uuu[3],15)+"\n");
            }
            for(int i=1;i<=ngrid;i+=4) {
                for(int j=0;j<4;j++) {
                    x=rpd*(i+j);
                    if(x<xx[0])
                        ggg[j]=-sl0*x;
                    else if(x>xx[npnts-1])
                        ggg[j]=-sln*x;
                    else {
                        if(x>xx[np])np++;
                        di=x-xx[np-1];
                        dj=xx[np]-x;
                        ggg[j]=-x*(6.0*(yy[np]-yy[np-1])+(di-dj)*
                        ((dd[np-1]+dj)* gg[np-1]+(dd[np-1]+di)*gg[np])-di*dj*
                        (gg[np]-gg[np-1]))/(6.0*dd[np-1]);
                    }
                }
                outStream.writeBytes(BML.fmt(ggg[0],15)+BML.fmt(ggg[1],15)+BML.fmt(ggg[2],15)
                +BML.fmt(ggg[3],15)+"\n");
            }
            outStream.close();
        }
        catch(Exception e) {
            println("error - writing file: "+fname);
            return -3;
        }
        numtab++;
        println(fname+" file created");
        return 0;
    }

    int tabgauss() {
        /*
********************************************************************

dl_poly/java  routine for constructing a TABLE file from a small
number of data points sampled from a potential function
using a gaussian fitting procedure

author - w.smith january 2001
copyright - daresbury laboratory

*********************************************************************
         */
        double[] ccc,eee,uuu,ggg;
        double rpd,x;
        int call;

        ccc=new double[3];
        eee=new double[3];
        uuu=new double[4];
        ggg=new double[4];
        zz=new double[npnts];
        aa=new double[npnts];
        dd=new double[npnts];
        gg=new double[npnts];

        // check on input variables

        if(ngrid<6) {
            println("Error - incorrect grid number specified for TABLE file");
            return -1;
        }
        else if(rcut<=0.0) {
            println("Error - incorrect cutoff specified for TABLE file");
            return -2;
        }

        // print current list of points and convert to dl_poly units

        println("Current list of spline points:");

        for(int i=0;i<npnts;i++) {
            println(BML.fmt(i+1,6)+BML.fmt(xx[i],8)+BML.fmt(yy[i],8));
            yy[i]=ecnvrt*yy[i];
        }

        // calculate gaussian parameters

        call=AML.gaussfit(npnts,ccc,eee,xx,yy,zz,aa,dd,gg);
        if(call<0)return -3;

        // print out potential arrays for comparison

        println("       rad        y(i)        z(i)");

        for(int i=0;i<npnts;i++) {
            x=xx[i];
            zz[i]=ccc[0]*Math.exp(-eee[0]*x*x)+ccc[1]*Math.exp(-eee[1]*x*x)+ccc[2]*Math.exp(-eee[2]*x*x);
            println(BML.fmt(x,12)+BML.fmt(yy[i],12)+BML.fmt(zz[i],12));
        }

        // grid increment

        rpd=rcut/(ngrid-4);

        // Write TABLE file

        fname="TABGSS."+numtab;
        title=name1+" "+name2+" "+"Gaussian Potential - DL_POLY units";

        try {
            DataOutputStream outStream = new DataOutputStream(new FileOutputStream(fname));
            outStream.writeBytes(title+"\n");
            outStream.writeBytes(BML.fmt(rpd,15)+BML.fmt(rcut,15)+BML.fmt(ngrid,10)+"\n");
            outStream.writeBytes(BML.fmt(name1,8)+BML.fmt(name2,8)+"\n");

            for(int i=1;i<=ngrid;i+=4) {
                for(int j=0;j<4;j++) {
                    x=rpd*(i+j);
                    uuu[j]=ccc[0]*Math.exp(-eee[0]*x*x)+ccc[1]*Math.exp(-eee[1]*x*x)+
                    ccc[2]*Math.exp(-eee[2]*x*x);
                }
                outStream.writeBytes(BML.fmt(uuu[0],15)+BML.fmt(uuu[1],15)+BML.fmt(uuu[2],15)
                +BML.fmt(uuu[3],15)+"\n");
            }
            for(int i=1;i<=ngrid;i+=4) {
                for(int j=0;j<4;j++) {
                    x=rpd*(i+j);
                    ggg[j]=2.0*x*x*(eee[0]*ccc[0]*Math.exp(-eee[0]*x*x)+
                    eee[1]*ccc[1]*Math.exp(-eee[1]*x*x)+
                    eee[2]*ccc[2]*Math.exp(-eee[2]*x*x));
                }
                outStream.writeBytes(BML.fmt(ggg[0],15)+BML.fmt(ggg[1],15)+BML.fmt(ggg[2],15)
                +BML.fmt(ggg[3],15)+"\n");
            }
            outStream.close();
        }
        catch(Exception e) {
            println("error - writing file: "+fname);
            return -3;
        }
        numtab++;
        println(fname+" file created");
        return 0;
    }

    int potSiO2() {
        /*
**********************************************************************

dl_poly/java utility

construct potential arrays for silica and related systems
based on Vessal's model for constant volume systems:
Phil. Mag. B 60 (1989) 753-775

copyright - daresbury laboratory
author    - w. smith january 2001

**********************************************************************
         */
        double[] uuu,ggg;
        double rpd,fac,elrc,vlrc,rrr,rsq,www;

        uuu=new double[4];
        ggg=new double[4];

        // check on input variables

        if(ngrid<6) {
            println("Error - incorrect grid number specified for TABLE file");
            return -1;
        }
        else if(rcut<=0.0) {
            println("Error - incorrect cutoff specified for TABLE file");
            return -2;
        }

        name1="Si4+    ";
        name2="O2-     ";
        rpd=rcut/(ngrid-4);
        fac=1.60217733e-19/1.6605402e-23;

        // Write TABLE file

        fname="TABSiO2."+numtab;
        title=name1+" "+name2+" "+"Silica Potential - DL_POLY units";

        try {
            DataOutputStream outStream = new DataOutputStream(new FileOutputStream(fname));
            outStream.writeBytes(title+"\n");
            outStream.writeBytes(BML.fmt(rpd,15)+BML.fmt(rcut,15)+BML.fmt(ngrid,10)+"\n");

            // silicon-oxygen potential

            elrc=-fac*25.0/(3.0*Math.pow(rcut,3));
            vlrc=fac*6.0*25.0/(3.0*Math.pow(rcut,3));

            outStream.writeBytes(BML.fmt(name1,8)+BML.fmt(name2,8)+BML.fmt(elrc,15)+
            BML.fmt(vlrc,15)+"\n");

            for(int i=1;i<=ngrid;i+=4) {
                for(int j=0;j<4;j++) {
                    rrr=rpd*(i+j);
                    rsq=rrr*rrr;
                    if(rrr<1.5)
                        uuu[j]=fac*990.617*Math.exp(-rrr/0.3297);
                    else if(rrr<2.5) {
                        uuu[j]=fac*(728.7931332+rrr*(-1623.5166500+
                        rrr*(1496.7899792+rrr*(-704.0045699+
                        rrr*(167.0744823-rrr*15.8848138)))));
                    }
                    else if(rrr<3.5) {
                        uuu[j]=fac*(0.6187898+rrr*(-0.6702747+rrr*(0.2214821-rrr*0.0233139)));
                    }
                    else
                        uuu[j]=-fac*25.0/Math.pow(rrr,6);
                }
                outStream.writeBytes(BML.fmt(uuu[0],15)+BML.fmt(uuu[1],15)+BML.fmt(uuu[2],15)
                +BML.fmt(uuu[3],15)+"\n");
            }
            for(int i=1;i<=ngrid;i+=4) {
                for(int j=0;j<4;j++) {
                    rrr=rpd*(i+j);
                    rsq=rrr*rrr;
                    if(rrr<1.5) {
                        www=fac*990.617*Math.exp(-rrr/0.3297);
                        ggg[j]=fac*rsq*www/(rrr*0.3297);
                    }
                    else if(rrr<2.5) {
                        ggg[j]=fac*(-rsq*(-1623.5166500+rrr*(1496.7899792*
                        2.0+rrr*(-704.0045699*3.0+rrr*(167.0744823*
                        4.0-rrr*15.8848138*5.0))))/rrr);
                    }
                    else if(rrr<3.5) {
                        ggg[j]=fac*(-rsq*(-0.67027470+rrr*(0.2214821*2.0-
                        rrr*0.0233139*3.0))/rrr);
                    }
                    else
                        ggg[j]=-fac*6.0*25.0*rsq/Math.pow(rrr,8);
                }
                outStream.writeBytes(BML.fmt(ggg[0],15)+BML.fmt(ggg[1],15)+BML.fmt(ggg[2],15)
                +BML.fmt(ggg[3],15)+"\n");
            }

            // oxygen-oxygen potential

            elrc=-fac*52.12/(3.0*Math.pow(rcut,3));
            vlrc=fac*6.0*52.12/(3.0*Math.pow(rcut,3));
            outStream.writeBytes(BML.fmt(name2,8)+BML.fmt(name2,8)+BML.fmt(elrc,15)+
            BML.fmt(vlrc,15)+"\n");

            for(int i=1;i<=ngrid;i+=4) {
                for(int j=0;j<4;j++) {
                    rrr=rpd*(i+j);
                    rsq=rrr*rrr;
                    if(rrr<2.9)
                        uuu[j]=fac*4511887.2*Math.exp(-rrr/0.149);
                    else if(rrr<3.6) {
                        uuu[j]=fac*(298.0818367+rrr*(-453.1251863+rrr*(275.1200462+
                        rrr*(-83.3698778+rrr*(12.6055881-rrr*0.7606781)))));
                    }
                    else if(rrr<4.2)
                        uuu[j]=fac*(1.5952103+rrr*(-1.2208242+rrr*(0.3052061-rrr*0.0251198)));
                    else
                        uuu[j]=-fac*52.12/Math.pow(rrr,6);
                }
                outStream.writeBytes(BML.fmt(uuu[0],15)+BML.fmt(uuu[1],15)+BML.fmt(uuu[2],15)
                +BML.fmt(uuu[3],15)+"\n");
            }
            for(int i=1;i<=ngrid;i+=4) {
                for(int j=0;j<4;j++) {
                    rrr=rpd*(i+j);
                    rsq=rrr*rrr;
                    if(rrr<2.9) {
                        www=fac*4511887.2*Math.exp(-rrr/0.149);
                        ggg[j]=fac*rsq*www/(rrr*0.149);
                    }
                    else if(rrr<3.6) {
                        ggg[j]=fac*(-rsq*(-453.1251863+rrr*(275.1200462*2.+
                        rrr*(-83.3698778*3.+rrr*(12.6055881*4.-
                        rrr*0.7606781*5.))))/rrr);
                    }
                    else if(rrr<4.2) {
                        ggg[j]=fac*(-rsq*(-1.2208242+rrr*(0.3052061*2.-
                        rrr*0.0251198*3.))/rrr);
                    }
                    else
                        ggg[j]=-fac*6.0*52.12*rsq/Math.pow(rrr,8);
                }
                outStream.writeBytes(BML.fmt(ggg[0],15)+BML.fmt(ggg[1],15)+BML.fmt(ggg[2],15)
                +BML.fmt(ggg[3],15)+"\n");
            }
            outStream.close();
        }
        catch(Exception e) {
            println("error - writing file: "+fname);
            return -3;
        }
        numtab++;
        println(fname+" file created");
        return 0;
    }

    int potAgI() {
        /*
**********************************************************************

dl_poly/java utility

construct potential arrays for silver iodide
based on model by Ray, Rahman and Vashishta:
Superionics and Solid Electrolytes, Academic Press 1989

copyright - daresbury laboratory
author    - w. smith january 2001

**********************************************************************
         */
        double[] uuu,ggg;
        double rpd,fac,elrc,vlrc,rrr,rsq;

        uuu=new double[4];
        ggg=new double[4];

        // check on input variables

        if(ngrid<6) {
            println("Error - incorrect grid number specified for TABLE file");
            return -1;
        }
        else if(rcut<=0.0) {
            println("Error - incorrect cutoff specified for TABLE file");
            return -2;
        }

        name1="Ag+     ";
        name2="I-      ";
        rpd=rcut/(ngrid-4);
        fac=1.60217733e-19/1.6605402e-23;

        // Write TABLE file

        fname="TABAgI."+numtab;
        title=name1+" "+name2+" "+"Silver Iodide Potential - DL_POLY units";

        try {
            DataOutputStream outStream = new DataOutputStream(new FileOutputStream(fname));
            outStream.writeBytes(title+"\n");
            outStream.writeBytes(BML.fmt(rpd,15)+BML.fmt(rcut,15)+BML.fmt(ngrid,10)+"\n");

            // silver-iodine potential

            elrc=fac*14.39*(114.48/(6.0*Math.pow(rcut,6))-1.1736/rcut);
            vlrc=-fac*14.39*(9.*114.48/(6.0*Math.pow(rcut,6))-4.*1.1736/rcut);
            outStream.writeBytes(BML.fmt(name1,8)+BML.fmt(name2,8)+BML.fmt(elrc,15)+
            BML.fmt(vlrc,15)+"\n");

            for(int i=1;i<=ngrid;i+=4) {
                for(int j=0;j<4;j++) {
                    rrr=rpd*(i+j);
                    rsq=rrr*rrr;
                    uuu[j]=fac*(14.39*(114.48/Math.pow(rrr,5)-1.1736)/(rsq*rsq));
                }
                outStream.writeBytes(BML.fmt(uuu[0],15)+BML.fmt(uuu[1],15)+BML.fmt(uuu[2],15)
                +BML.fmt(uuu[3],15)+"\n");
            }
            for(int i=1;i<=ngrid;i+=4) {
                for(int j=0;j<4;j++) {
                    rrr=rpd*(i+j);
                    rsq=rrr*rrr;
                    ggg[j]=fac*(14.39*(9.0*114.48/Math.pow(rrr,5)-4.0*1.1736)/(rsq*rsq));
                }
                outStream.writeBytes(BML.fmt(ggg[0],15)+BML.fmt(ggg[1],15)+BML.fmt(ggg[2],15)
                +BML.fmt(ggg[3],15)+"\n");
            }

            // iodine-iodine potential

            elrc=fac*14.39*(446.64/(4.0*Math.pow(rcut,4))-6.9331/(3.0*Math.pow(rcut,3))-2.3472/rcut);
            vlrc=-fac*14.39*(7.*446.64/(4.0*Math.pow(rcut,4))-6.*6.9331/(3.0*Math.pow(rcut,3))-4.0*2.3472/rcut);
            outStream.writeBytes(BML.fmt(name2,8)+BML.fmt(name2,8)+BML.fmt(elrc,15)+
            BML.fmt(vlrc,15)+"\n");

            for(int i=1;i<=ngrid;i+=4) {
                for(int j=0;j<4;j++) {
                    rrr=rpd*(i+j);
                    rsq=rrr*rrr;
                    uuu[j]=fac*(14.39*(446.64/Math.pow(rrr,3)-6.9331/rsq-2.3472)/(rsq*rsq));
                }
                outStream.writeBytes(BML.fmt(uuu[0],15)+BML.fmt(uuu[1],15)+BML.fmt(uuu[2],15)
                +BML.fmt(uuu[3],15)+"\n");
            }
            for(int i=1;i<=ngrid;i+=4) {
                for(int j=0;j<4;j++) {
                    rrr=rpd*(i+j);
                    rsq=rrr*rrr;
                    ggg[j]=fac*(14.39*(7.*446.64/Math.pow(rrr,3)-6.*6.9331/rsq-4.*2.3472)/(rsq*rsq));
                }
                outStream.writeBytes(BML.fmt(ggg[0],15)+BML.fmt(ggg[1],15)+BML.fmt(ggg[2],15)
                +BML.fmt(ggg[3],15)+"\n");
            }

            // silver-silver potential

            elrc=fac*14.39*(0.014804/(8.0*Math.pow(rcut,8)));
            vlrc=-fac*11.0*14.39*(0.014804/(8.0*Math.pow(rcut,8)));
            outStream.writeBytes(BML.fmt(name1,8)+BML.fmt(name1,10)+BML.fmt(elrc,15)+
            BML.fmt(vlrc,15)+"\n");

            for(int i=1;i<=ngrid;i+=4) {
                for(int j=0;j<4;j++) {
                    rrr=rpd*(i+j);
                    rsq=rrr*rrr;
                    uuu[j]=fac*(14.39*(0.014804/Math.pow(rrr,11)));
                }
                outStream.writeBytes(BML.fmt(uuu[0],15)+BML.fmt(uuu[1],15)+BML.fmt(uuu[2],15)
                +BML.fmt(uuu[3],15)+"\n");
            }
            for(int i=1;i<=ngrid;i+=4) {
                for(int j=0;j<4;j++) {
                    rrr=rpd*(i+j);
                    rsq=rrr*rrr;
                    ggg[j]=fac*(14.39*(11.0*0.014804/Math.pow(rrr,11)));
                }
                outStream.writeBytes(BML.fmt(ggg[0],15)+BML.fmt(ggg[1],15)+BML.fmt(ggg[2],15)
                +BML.fmt(ggg[3],15)+"\n");
            }
            outStream.close();
        }
        catch(Exception e) {
            println("error - writing file: "+fname);
            return -3;
        }
        numtab++;
        println(fname+" file created");
        return 0;
    }
}
