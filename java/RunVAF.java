import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;

public class RunVAF extends Basic implements ActionListener {
        /*
*********************************************************************

dl_poly/java GUI class to calculate velocity autocorrelation

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
    public static RunVAF job;
    private static GUI home;
    private static String atname;
    private static boolean form;
    private static double[] xx,yy;
    private static int[] imd,msm;
    private static String[] name;
    private static double[] chge,weight;
    private static double[][] xyz,vel,frc,vaf;
    private static double[][][] vaf0;
    private static int npnts,nconf,nvaf,isampl,iovaf;
    private static JTextField atom1,history,configs,length,sample,origin;
    private static JCheckBox format;
    private static JButton run,close;

    // Define the Graphical User Interface

    public RunVAF() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        setTitle("VAF Panel");

        getContentPane().setBackground(art.back);
        getContentPane().setForeground(art.fore);
        setDefaultCloseOperation(DISPOSE_ON_CLOSE);
        setFont(fontMain);
        GridBagLayout grd = new GridBagLayout();
        GridBagConstraints gbc = new GridBagConstraints();
        getContentPane().setLayout(grd);

        gbc.fill=GridBagConstraints.BOTH;

        // Define the Run button

        run = new JButton("Run");
        run.setBackground(art.butn);
        run.setForeground(art.butf);
        fix(run,grd,gbc,0,0,1,1);

        fix(new JLabel("  "),grd,gbc,1,0,1,1);

        // Define the Close button

        close = new JButton("Close");
        close.setBackground(art.butn);
        close.setForeground(art.butf);
        fix(close,grd,gbc,2,0,1,1);


        // Name of HISTORY file

        JLabel lab1 = new JLabel("Required HISTORY file:",JLabel.LEFT);
        fix(lab1,grd,gbc,0,1,3,1);
        history = new JTextField(18);
        history.setBackground(art.scrn);
        history.setForeground(art.scrf);
        fix(history,grd,gbc,0,2,3,1);

        // History file format

        format=new JCheckBox("Formatted");
        format.setBackground(art.back);
        format.setForeground(art.fore);
        //fix(format,grd,gbc,0,3,1,1);
        JLabel lab2 = new JLabel("file?",JLabel.LEFT);
        //fix(lab2,grd,gbc,1,3,2,1);

        // Name of atomic species

        JLabel lab3 = new JLabel("Atom name:",JLabel.LEFT);
        fix(lab3,grd,gbc,0,4,2,1);
        atom1 = new JTextField(8);
        atom1.setBackground(art.scrn);
        atom1.setForeground(art.scrf);
        fix(atom1,grd,gbc,2,4,1,1);

        // Number of configurations

        JLabel lab4 = new JLabel("No. configurations:",JLabel.LEFT);
        fix(lab4,grd,gbc,0,5,2,1);
        configs = new JTextField(8);
        configs.setBackground(art.scrn);
        configs.setForeground(art.scrf);
        fix(configs,grd,gbc,2,5,1,1);

        // VAF array length

        JLabel lab5 = new JLabel("VAF array length:",JLabel.LEFT);
        fix(lab5,grd,gbc,0,6,2,1);
        length = new JTextField(8);
        length.setBackground(art.scrn);
        length.setForeground(art.scrf);
        fix(length,grd,gbc,2,6,1,1);

        // Sampling interval

        JLabel lab6 = new JLabel("Sampling interval:",JLabel.LEFT);
        fix(lab6,grd,gbc,0,7,2,1);
        sample = new JTextField(8);
        sample.setBackground(art.scrn);
        sample.setForeground(art.scrf);
        fix(sample,grd,gbc,2,7,1,1);

        // Origin interval

        JLabel lab7 = new JLabel("Origin interval:",JLabel.LEFT);
        fix(lab7,grd,gbc,0,8,2,1);
        origin = new JTextField(8);
        origin.setBackground(art.scrn);
        origin.setForeground(art.scrf);
        fix(origin,grd,gbc,2,8,1,1);

        // Register action buttons

        run.addActionListener(this);
        close.addActionListener(this);

    }

    public RunVAF(GUI here) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        home=here;
        println("Activated VAF panel");
        job=new RunVAF();
        job.pack();
        job.setVisible(true);
        npnts=0;
        form=true;
        atname="ALL";
        fname="HISTORY";
        nconf=1000;
        nvaf=512;
        isampl=1;
        iovaf=1;
        format.setSelected(form);
        atom1.setText(atname);
        history.setText(fname);
        configs.setText(String.valueOf(nconf));
        length.setText(String.valueOf(nvaf));
        sample.setText(String.valueOf(isampl));
        origin.setText(String.valueOf(iovaf));
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
        if (arg.equals("Run")) {
            form=format.isSelected();
            atname=atom1.getText();
            fname=history.getText();
            nconf=BML.giveInteger(configs.getText(),1);
            nvaf=BML.giveInteger(length.getText(),1);
            isampl=BML.giveInteger(sample.getText(),1);
            iovaf=BML.giveInteger(origin.getText(),1);
            println("Started VAF calculation .....");
            npnts=calcVAF();
            if(npnts>0) {
                vafXY(npnts,atname);
                if(graf != null)
                    graf.job.dispose();
                graf=new GraphDraw(home);
                graf.xlabel.setText("Time (ps)");
                graf.ylabel.setText("VAF");
                graf.plabel.setText("VAF of "+atname.trim());
                graf.extraPlot(npnts,xx,yy);
            }
        }
        else if (arg.equals("Close")) {
            job.dispose();
        }
    }
    int calcVAF() {
        /*
*********************************************************************

dl_poly/java routine to calculate velocity autocorrelation function
for selected atoms from dl_poly HISTORY file

copyright - daresbury laboratory
author    - w.smith march 2001

*********************************************************************
         */
        boolean all;
        int m,n,nat,natms,novaf,lsr,msr,nsvaf,imcon,iconf;
        double rnorm,vsum,tstep;
        LineNumberReader lnr=null;
        double cell[]=new double[9];
        double info[]=new double[10];

        nat=0;
        npnts=0;
        all=false;
        tstep=0.0;
        if(atname.toUpperCase().equals("ALL"))all=true;

        if(nvaf%iovaf != 0) {

            nvaf=iovaf*(nvaf/iovaf);
            println("Warning - vaf array dimension reset to "+BML.fmt(nvaf,8));
        }
        novaf=nvaf/iovaf;
        imd=new int[nvaf];
        msm=new int[nvaf];
        xx=new double[nvaf];
        yy=new double[nvaf];

        // write control variables

        println("Name of target HISTORY file   : "+fname);
        println("Label  of atom  of interest   : "+atname);
        println("Length of correlation arrays  : "+BML.fmt(nvaf,8));
        println("Number of configurations      : "+BML.fmt(nconf,8));
        println("Sampling interval             : "+BML.fmt(isampl,8));
        println("Interval between origins      : "+BML.fmt(iovaf,8));

        // initialise vaf variables

        lsr=0;
        msr=-1;
        nsvaf=0;

        // initialise control parameters for HISTORY file reader

        info[0]=0.0;
        info[1]=999999.;
        info[2]=1.0;
        info[3]=0.0;
        info[4]=0.0;
        info[5]=0.0;
        if(form) {
            lnr=hread(fname,name,lnr,info,cell,chge,weight,xyz,vel,frc);
            if(BML.nint(info[3])<0 && BML.nint(info[3])!=-1) {
                println("Error - cannot open HISTORY file");
                return -1;
            }
        }
        else {
            println("Error - unformatted read option not active");
            return -2;
        }

        natms=BML.nint(info[7]);

        //initialise vaf arrays

        name=new String[natms];
        chge=new double[natms];
        weight=new double[natms];
        xyz=new double[3][natms];
        vel=new double[3][natms];
        vaf=new double[nvaf][natms];
        vaf0=new double[nvaf][natms][3];
        for(int j=0;j<nvaf;j++) {
            msm[j]=0;
            for(int i=0;i<natms;i++) {
                vaf[j][i]=0.0;
                vaf0[j][i][0]=0.0;
                vaf0[j][i][1]=0.0;
                vaf0[j][i][2]=0.0;
            }
        }
        OUT:
            for(iconf=0;iconf<nconf;iconf++) {
                lnr=hread(fname,name,lnr,info,cell,chge,weight,xyz,vel,frc);
                if(BML.nint(info[3])<0 && BML.nint(info[3])!=-1) {
                    println("Error - HISTORY file data error");
                    info[0]=-1.0;
                    return -3;
                }
                if(lnr == null)break OUT;
                if(iconf==0)info[9]=info[6];
                if(iconf==1)tstep=info[8]*(info[6]-info[9]);
                if(BML.nint(info[3])==-1)break OUT;

                if(iconf%isampl==0) {
                    n=0;
                    for(int i=0;i<natms;i++) {
                        if(all || name[i].equals(atname)) {
                            if(n==natms) {
                                println("Error - too many atoms of specified type");
                                info[0]=-1.0;
                                lnr=hread(fname,name,lnr,info,cell,chge,weight,xyz,vel,frc);
                                return -4;
                            }
                            vel[0][n]=vel[0][i];
                            vel[1][n]=vel[1][i];
                            vel[2][n]=vel[2][i];
                            n++;
                        }
                    }
                    nat=n;
                    if(iconf==0) {
                        println("Number of atoms of selected type : "+BML.fmt(nat,8));
                    }
                    if(nat == 0) {
                        println("Error - zero atoms of specified type");
                        info[0]=-1.0;
                        lnr=hread(fname,name,lnr,info,cell,chge,weight,xyz,vel,frc);
                        return -5;
                    }

                    // calculate velocity autocorrelation

                    if(nsvaf%iovaf==0) {
                        lsr=Math.min(lsr+1,novaf);
                        msr=(msr+1)%novaf;
                        imd[msr]=0;
                        for(int i=0;i<nat;i++) {
                            vaf0[msr][i][0]=vel[0][i];
                            vaf0[msr][i][1]=vel[1][i];
                            vaf0[msr][i][2]=vel[2][i];
                        }
                    }

                    nsvaf++;

                    for(int j=0;j<lsr;j++) {
                        m=imd[j];
                        imd[j]=m+1;
                        msm[m]++;
                        for(int i=0;i<nat;i++) {
                            vaf[m][i]+=(vel[0][i]*vaf0[j][i][0]+vel[1][i]*vaf0[j][i][1]+
                            vel[2][i]*vaf0[j][i][2]);
                        }
                    }
                }
            }
            if(iconf==nconf-1) {
                info[0]=-1.0;
                lnr=hread(fname,name,lnr,info,cell,chge,weight,xyz,vel,frc);
            }

            if(BML.nint(info[3])==-1) iconf--;

            npnts=Math.min(nsvaf,nvaf);
            println("Number of configurations read: "+BML.fmt(iconf,8));

            // normalise mean square displacement

            for(int i=0;i<nat;i++) {
                rnorm=((double)msm[0])/vaf[0][i];
                for(int j=0;j<npnts;j++)
                    vaf[j][i]*=(rnorm/msm[j]);
            }

            // average velocity autocorrelation function

            for(int j=0;j<npnts;j++) {
                vsum=0.0;
                for(int i=0;i<nat;i++)
                    vsum+=vaf[j][i];
                xx[j]=tstep*j*isampl;
                yy[j]=vsum/nat;
            }
            return npnts;

    }
    void vafXY(int npts,String anam1) {
        /*
*********************************************************************

dl_poly/java GUI routine to create a VAF XY file

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        String fname;
        String[] header;
        int nhead=4,call;

        header=new String[nhead];
        fname="VAF"+String.valueOf(numvaf)+".XY";
        numvaf++;
        header[0]=" VAF Plotting Program";
        header[1]=" VAF Plot: "+anam1.trim();
        header[2]=" Time (ps)";
        header[3]=" VAF";
        call=putXY(fname,header,nhead,npts,xx,yy);
        if(call==0)
            println("PLOT file "+fname+" created");
    }
}
