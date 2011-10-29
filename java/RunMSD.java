import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;

public class RunMSD extends Basic implements ActionListener {
        /*
*********************************************************************

dl_poly/java GUI class to calculate mean square displacement

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
    public static RunMSD job;
    private static GUI home;
    private static String atname;
    private static boolean form;
    private static double[] xx,yy;
    private static String[] name;
    private static int[] imd,msm;
    private static double[] chge,weight;
    private static double[][] xyz,vel,frc,xy0,xy1,acm,msd;
    private static double[][][] msd0;
    private static int npnts,nconf,nmsd,isampl,iomsd;
    private static JTextField atom1,history,configs,length,sample,origin;
    private static JCheckBox format;
    private static JButton run,close;

    // Define the Graphical User Interface

    public RunMSD() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        setTitle("MSD Panel");

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

        // MSD array length

        JLabel lab5 = new JLabel("MSD array length:",JLabel.LEFT);
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

    public RunMSD(GUI here) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        home=here;
        println("Activated MSD panel");
        job=new RunMSD();
        job.pack();
        job.setVisible(true);
        npnts=0;
        form=true;
        atname="ALL";
        fname="HISTORY";
        nconf=1000;
        nmsd=512;
        isampl=1;
        iomsd=1;
        format.setSelected(form);
        atom1.setText(atname);
        history.setText(fname);
        configs.setText(String.valueOf(nconf));
        length.setText(String.valueOf(nmsd));
        sample.setText(String.valueOf(isampl));
        origin.setText(String.valueOf(iomsd));
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
            nmsd=BML.giveInteger(length.getText(),1);
            isampl=BML.giveInteger(sample.getText(),1);
            iomsd=BML.giveInteger(origin.getText(),1);
            println("Started MSD calculation .....");
            npnts=calcMSD();
            if(npnts>0) {
                msdXY(npnts,atname);
                if(graf != null)
                    graf.job.dispose();
                graf=new GraphDraw(home);
                graf.xlabel.setText("Time (ps)");
                graf.ylabel.setText("MSD (A^2)");
                graf.plabel.setText("MSD of "+atname.trim());
                graf.extraPlot(npnts,xx,yy);
            }
        }
        else if (arg.equals("Close")) {
            job.dispose();
        }
    }
    int calcMSD() {
        /*
*********************************************************************

dl_poly/java routine to calculate mean square displacement
for selected atoms from dl_poly HISTORY file

copyright - daresbury laboratory
author    - w.smith march 2001

*********************************************************************
         */
        boolean all;
        int m,n,nat,natms,nomsd,lsr,msr,nsmsd,imcon,iconf;
        double f1,f2,uuu,vvv,www,rmsx,rmsy,rmsz,rr2,tstep;
        LineNumberReader lnr=null;
        double cell[]=new double[9];
        double info[]=new double[10];
        double rcell[]=new double[9];
        double avcell[]=new double[9];

        nat=0;
        npnts=0;
        all=false;
        tstep=0.0;
        if(atname.toUpperCase().equals("ALL"))all=true;

        if(nmsd%iomsd != 0) {

            nmsd=iomsd*(nmsd/iomsd);
            println("Warning - msd array dimension reset to "+BML.fmt(nmsd,8));
        }
        nomsd=nmsd/iomsd;
        xx=new double[nmsd];
        yy=new double[nmsd];
        imd=new int[nmsd];
        msm=new int[nmsd];

        // write control variables

        println("Name of target HISTORY file   : "+fname);
        println("Label  of atom  of interest   : "+atname);
        println("Length of correlation arrays  : "+BML.fmt(nmsd,8));
        println("Number of configurations      : "+BML.fmt(nconf,8));
        println("Sampling interval             : "+BML.fmt(isampl,8));
        println("Interval between origins      : "+BML.fmt(iomsd,8));

        // initialize cell vector arrays

        for(int i=0;i<9;i++) {
            cell[i]=0.0;
            avcell[i]=0.0;
        }
        cell[0]=1.0;
        cell[4]=1.0;
        cell[8]=1.0;

        // initialise msd variables

        lsr=0;
        msr=-1;
        nsmsd=0;

        // initialise control parameters for HISTORY file reader

        info[0]=0.0;
        info[1]=999999.;
        info[2]=0.0;
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

        // initialise msd arrays

        natms=BML.nint(info[7]);

        name=new String[natms];
        chge=new double[natms];
        weight=new double[natms];
        xyz=new double[3][natms];
        xy0=new double[3][natms];
        xy1=new double[3][natms];
        acm=new double[3][natms];
        msd=new double[nmsd][natms];
        msd0=new double[nmsd][natms][3];
        for(int i=0;i<natms;i++) {
            acm[0][i]=0.0;
            acm[1][i]=0.0;
            acm[2][i]=0.0;
        }

        for(int j=0;j<nmsd;j++) {
            msm[j]=0;
            for(int i=0;i<natms;i++) {
                msd[j][i]=0.0;
                msd0[j][i][0]=0.0;
                msd0[j][i][1]=0.0;
                msd0[j][i][2]=0.0;
            }
        }

        OUT:
            for(iconf=0;iconf<nconf;iconf++) {
                lnr=hread(fname,name,lnr,info,cell,chge,weight,xyz,vel,frc);
                if(BML.nint(info[3])<0 && BML.nint(info[3])!=-1) {
                    println("Error - HISTORY file data error");
                    info[0]=-1.0;
                    lnr=hread(fname,name,lnr,info,cell,chge,weight,xyz,vel,frc);
                    return -3;
                }
                if(lnr == null)break OUT;
                if(iconf==0)info[9]=info[6];
                if(iconf==1)tstep=info[8]*(info[6]-info[9]);
                if(BML.nint(info[3])==-1)break OUT;

                imcon=BML.nint(info[5]);
                if(imcon>=1 && imcon<=3) {
                    rcell=AML.invert(cell);
                }
                else {
                    println("Error - incorrect periodic boundary condition");
                    info[0]=-1.0;
                    lnr=hread(fname,name,lnr,info,cell,chge,weight,xyz,vel,frc);
                    return -4;
                }

                n=0;
                for(int i=0;i<natms;i++) {
                    if(all || name[i].equals(atname)) {
                        if(n==natms) {
                            println("Error - too many atoms of specified type");
                            info[0]=-1.0;
                            lnr=hread(fname,name,lnr,info,cell,chge,weight,xyz,vel,frc);
                            return -5;
                        }
                        xy1[0][n]=xyz[0][i]*rcell[0]+xyz[1][i]*rcell[3]+xyz[2][i]*rcell[6];
                        xy1[1][n]=xyz[0][i]*rcell[1]+xyz[1][i]*rcell[4]+xyz[2][i]*rcell[7];
                        xy1[2][n]=xyz[0][i]*rcell[2]+xyz[1][i]*rcell[5]+xyz[2][i]*rcell[8];
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
                    return -6;
                }

                // running average of cell vectors

                f1=((double)iconf)/(iconf+1);
                f2=1.0/(iconf+1);
                for(int i=0;i<9;i++) {
                    avcell[i]=f1*avcell[i]+f2*cell[i];
                }

                if(iconf>0) {

                    // accumulate incremental distances

                    for(int i=0;i<nat;i++) {
                        uuu=xy1[0][i]-xy0[0][i];
                        vvv=xy1[1][i]-xy0[1][i];
                        www=xy1[2][i]-xy0[2][i];

                        uuu=uuu-BML.nint(uuu);
                        vvv=vvv-BML.nint(vvv);
                        www=www-BML.nint(www);

                        acm[0][i]+=(uuu*avcell[0]+vvv*avcell[3]+www*avcell[6]);
                        acm[1][i]+=(uuu*avcell[1]+vvv*avcell[4]+www*avcell[7]);
                        acm[2][i]+=(uuu*avcell[2]+vvv*avcell[5]+www*avcell[8]);
                    }
                }

                for(int i=0;i<nat;i++) {
                    xy0[0][i]=xy1[0][i];
                    xy0[1][i]=xy1[1][i];
                    xy0[2][i]=xy1[2][i];
                }


                // calculate mean square displacement

                if(iconf > 0) {
                    if(iconf%isampl==0) {
                        if(nsmsd%iomsd==0) {
                            lsr=Math.min(lsr+1,nomsd);
                            msr=(msr+1)%nomsd;
                            imd[msr]=0;
                            for(int i=0;i<nat;i++) {
                                msd0[msr][i][0]=0.0;
                                msd0[msr][i][1]=0.0;
                                msd0[msr][i][2]=0.0;
                            }
                        }

                        nsmsd++;
                        for(int j=0;j<lsr;j++) {
                            m=imd[j];
                            imd[j]=m+1;
                            msm[m]++;

                            for(int i=0;i<nat;i++) {
                                rmsx=msd0[j][i][0]+acm[0][i];
                                rmsy=msd0[j][i][1]+acm[1][i];
                                rmsz=msd0[j][i][2]+acm[2][i];
                                msd[m][i]+=(rmsx*rmsx+rmsy*rmsy+rmsz*rmsz);
                                msd0[j][i][0]=rmsx;
                                msd0[j][i][1]=rmsy;
                                msd0[j][i][2]=rmsz;
                            }
                        }

                        for(int i=0;i<natms;i++) {
                            acm[0][i]=0.0;
                            acm[1][i]=0.0;
                            acm[2][i]=0.0;
                        }
                    }
                }
                if(iconf==nconf-1) {
                    info[0]=-1.0;
                    lnr=hread(fname,name,lnr,info,cell,chge,weight,xyz,vel,frc);
                }
            }

            if(BML.nint(info[3])==-1) iconf--;

            npnts=Math.min(nsmsd,nmsd);
            println("Number of configurations read: "+BML.fmt(iconf,8));

            // normalise mean square displacement

            for(int j=0;j<npnts;j++) {
                rr2=0.0;
                for(int i=0;i<nat;i++) {
                    rr2+=msd[j][i];
                }
                rr2=rr2/msm[j];
                xx[j]=tstep*isampl*j;
                yy[j]=rr2/nat;
            }
            return npnts;

    }
    void msdXY(int npts,String anam1) {
        /*
*********************************************************************

dl_poly/java GUI routine to create a MSD XY file

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        String fname;
        String[] header;
        int nhead=4,call;

        header=new String[nhead];
        fname="MSD"+String.valueOf(nummsd)+".XY";
        nummsd++;
        header[0]=" MSD Plotting Program";
        header[1]=" MSD Plot: "+anam1.trim();
        header[2]=" Time (ps)";
        header[3]=" MSD";
        call=putXY(fname,header,nhead,npts,xx,yy);
        if(call==0)
            println("PLOT file "+fname+" created");
    }
}
