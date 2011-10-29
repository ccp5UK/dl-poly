import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;

public class SkwCalc extends Basic implements ActionListener {
        /*
*********************************************************************

dl_poly/java GUI class to calculate dynamic structure factors

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
    public static SkwCalc job;
    private static GUI home;
    private static SkwPlot skwplt=null;
    private static int nconf,lencor,isampl,kmax,klim,last;
    private static JTextField history,configs,length,sample,kvect;
    private static JCheckBox charge;
    private static boolean lchg,safe;
    private static JButton run,close,plot;
    private static double TPI=2*Math.PI;
    private static String[] name;
    private static double[] chge,weight;
    private static double[][] xyz,vel,frc;
    private static String sname=null;

    // Define the Graphical User Interface

    public SkwCalc() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        setTitle("S(k,w) Calculator");

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

        // Instruction label 1

        JLabel lab1 = new JLabel("Required HISTORY file:",JLabel.LEFT);
        fix(lab1,grd,gbc,0,1,3,1);

        // Name of HISTORY file

        history = new JTextField(18);
        history.setBackground(art.scrn);
        history.setForeground(art.scrf);
        fix(history,grd,gbc,0,2,3,1);

        // Number of configurations

        JLabel lab2 = new JLabel("No. configurations:",JLabel.LEFT);
        fix(lab2,grd,gbc,0,3,2,1);
        configs = new JTextField(8);
        configs.setBackground(art.scrn);
        configs.setForeground(art.scrf);
        fix(configs,grd,gbc,2,3,1,1);

        // SKW array lengths

        JLabel lab3 = new JLabel("SKW array lengths:",JLabel.LEFT);
        fix(lab3,grd,gbc,0,4,2,1);
        length = new JTextField(8);
        length.setBackground(art.scrn);
        length.setForeground(art.scrf);
        fix(length,grd,gbc,2,4,1,1);

        // Sampling interval

        JLabel lab4 = new JLabel("Sampling interval:",JLabel.LEFT);
        fix(lab4,grd,gbc,0,5,2,1);
        sample = new JTextField(8);
        sample.setBackground(art.scrn);
        sample.setForeground(art.scrf);
        fix(sample,grd,gbc,2,5,1,1);

        // Maximum k-vector index

        JLabel lab5 = new JLabel("Max k-vector index:",JLabel.LEFT);
        fix(lab5,grd,gbc,0,6,2,1);
        kvect = new JTextField(8);
        kvect.setBackground(art.scrn);
        kvect.setForeground(art.scrf);
        fix(kvect,grd,gbc,2,6,1,1);

        // Use charges option

        JLabel lab6 = new JLabel("Use charges?",JLabel.LEFT);
        fix(lab6,grd,gbc,0,7,2,1);
        charge=new JCheckBox("    ");
        charge.setBackground(art.back);
        charge.setForeground(art.fore);
        fix(charge,grd,gbc,2,7,1,1);


        // Define the Close button

        close = new JButton("Close");
        close.setBackground(art.butn);
        close.setForeground(art.butf);
        fix(close,grd,gbc,0,8,1,1);

        // Define the Plot button

        plot = new JButton("Plot");
        plot.setBackground(art.butn);
        plot.setForeground(art.butf);
        fix(plot,grd,gbc,2,8,1,1);

        // Register action buttons

        run.addActionListener(this);
        close.addActionListener(this);
        plot.addActionListener(this);

    }

    public SkwCalc(GUI here) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        home=here;
        println("Activated panel for calculating S(k,w)");
        job=new SkwCalc();
        job.pack();
        job.setVisible(true);
        fname="HISTORY";
        nconf=1000;
        lencor=512;
        isampl=1;
        kmax=1;
        safe=true;
        lchg=false;
        history.setText(fname);
        charge.setSelected(lchg);
        configs.setText(String.valueOf(nconf));
        length.setText(String.valueOf(lencor));
        sample.setText(String.valueOf(isampl));
        kvect.setText(String.valueOf(kmax));
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
            fname=history.getText();
            lchg=charge.isSelected();
            nconf=BML.giveInteger(configs.getText(),1);
            lencor=BML.giveInteger(length.getText(),1);
            isampl=BML.giveInteger(sample.getText(),1);
            kmax=BML.giveInteger(kvect.getText(),1);
            if(isampl<1)isampl=1;
            skwrun();
            numskw++;
        }
        else if (arg.equals("Close")) {
            job.dispose();
        }
        else if (arg.equals("Plot")) {
            if(skwplt != null)
                skwplt.job.dispose();
            skwplt=new SkwPlot(home,sname);
        }
    }
    void skwrun(){
        /*
*********************************************************************

dl_poly/java program for calculating the dynamic structure factor

copyright - daresbury laboratory
author    - w smith 2007

*********************************************************************
         */
        println("DL_POLY Density Correlation Program");

        // write control variables

        println("Name of target HISTORY file   : "+fname);
        println("Length of correlation array   : "+BML.fmt(lencor,8));
        println("Number of configurations      : "+BML.fmt(nconf,8));
        println("Interval between samples      : "+BML.fmt(isampl,8));
        println("Maximum k-vector index        : "+BML.fmt(kmax,8));
        if(lchg){
            println("Charge density option selected");
        }
        else{
            println("Particle density option selected");
        }
        klim=(int)((Math.pow(2*kmax+1,3)-1)/2);
        safe=forden();
        if(safe)safe=correl();
        if(safe)safe=denfft();
    }
    boolean forden(){
        /*
*********************************************************************

dl_poly/java routine to calculate fourier transform of
paricle or charge density from dl_poly HISTORY file

copyright - daresbury laboratory
author    - w.smith nov 2007

*********************************************************************
         */
        LineNumberReader lnr=null;
        DataOutputStream spc=null;
        boolean newfile=true;
        int iconf,natms,k,imcon,lmin,mmin,kkk,nn;
        double tx,ty,tz,rnrm,time0,time;
        double[] work;
        double cell[]=new double[9];
        double info[]=new double[10];
        double dxyz[][]=new double[2][klim];

        // process the HISTORY file data

        for(int i=0;i<9;i++)
            cell[i]=0.0;
        cell[0]=1.0;
        cell[4]=1.0;
        cell[8]=1.0;

        // initialise control parameters for HISTORY file reader

        info[0]=0.0;
        info[1]=999999;
        info[2]=0.0;
        info[3]=0.0;
        info[4]=0.0;
        info[5]=0.0;
        lnr=hread(fname,name,lnr,info,cell,chge,weight,xyz,vel,frc);
        if(BML.nint(info[3])<0 && BML.nint(info[3])!=-1) {
            println("Error - HISTORY file data error");
            return false;
        }
        natms=BML.nint(info[7]);
        imcon=BML.nint(info[5]);
        if(imcon < 1 || imcon > 3) {
            println("Error - incorrect periodic boundary condition");
            info[0]=-1.0;
            lnr=hread(fname,name,lnr,info,cell,chge,weight,xyz,vel,frc);
            return false;
        }

        nn=0;
        rnrm=1/Math.sqrt(natms);
        name=new String[natms];
        chge=new double[natms];
        weight=new double[natms];
        xyz=new double[3][natms];
        double[] elc=new double[natms];
        double[] els=new double[natms];
        double[] emc=new double[natms];
        double[] ems=new double[natms];
        double[] enc=new double[natms];
        double[] ens=new double[natms];
        double[][] ecx=new double[natms][kmax+1];
        double[][] esx=new double[natms][kmax+1];
        double[][] ecy=new double[natms][kmax+1];
        double[][] esy=new double[natms][kmax+1];
        double[][] ecz=new double[natms][kmax+1];
        double[][] esz=new double[natms][kmax+1];

        OUT:
            for(iconf=0;iconf<nconf;iconf++) {
                lnr=hread(fname,name,lnr,info,cell,chge,weight,xyz,vel,frc);
                if(BML.nint(info[3])<0 && BML.nint(info[3])!=-1) {
                    println("Error - HISTORY file data error");
                    info[0]=-1.0;
                    lnr=hread(fname,name,lnr,info,cell,chge,weight,xyz,vel,frc);
                    return false;
                }
                if(lnr == null) break OUT;
                if(iconf==0)info[9]=info[6];
                nn++;

                //process data at required sampling interval

                time=info[8]*(info[6]-info[9]);
                if(iconf%isampl==0){

                    work=AML.invert(cell);
                    for(int i=0;i<natms;i++) {
                        tx=xyz[0][i];
                        ty=xyz[1][i];
                        tz=xyz[2][i];
                        xyz[0][i]=work[0]*tx+work[3]*ty+work[6]*tz;
                        xyz[1][i]=work[1]*tx+work[4]*ty+work[7]*tz;
                        xyz[2][i]=work[2]*tx+work[5]*ty+work[8]*tz;
                    }
                    if(BML.nint(info[3])==-1)break OUT;

                    //calculate fourier exponential terms

                    for(int i=0;i<natms;i++){
                        ecx[i][0]=1;
                        esx[i][0]=0;
                        ecy[i][0]=1;
                        esy[i][0]=0;
                        ecz[i][0]=1;
                        esz[i][0]=0;
                        elc[i]=Math.cos(TPI*xyz[0][i]);
                        els[i]=Math.sin(TPI*xyz[0][i]);
                        emc[i]=Math.cos(TPI*xyz[1][i]);
                        ems[i]=Math.sin(TPI*xyz[1][i]);
                        enc[i]=Math.cos(TPI*xyz[2][i]);
                        ens[i]=Math.sin(TPI*xyz[2][i]);
                    }

                    for(int l=0;l<kmax;l++){
                        for(int i=0;i<natms;i++){
                            ecx[i][l+1]=ecx[i][l]*elc[i]-esx[i][l]*els[i];
                            esx[i][l+1]=ecx[i][l]*els[i]+esx[i][l]*elc[i];
                            ecy[i][l+1]=ecy[i][l]*emc[i]-esy[i][l]*ems[i];
                            esy[i][l+1]=ecy[i][l]*ems[i]+esy[i][l]*emc[i];
                            ecz[i][l+1]=ecz[i][l]*enc[i]-esz[i][l]*ens[i];
                            esz[i][l+1]=ecz[i][l]*ens[i]+esz[i][l]*enc[i];
                        }
                    }

                    //loop over k vectors

                    kkk=0;
                    mmin=0;
                    lmin=1;
                    for(int n=0;n<=kmax;n++){
                        for(int i=0;i<natms;i++){
                            enc[i]=ecz[i][n];
                            ens[i]=esz[i][n];
                        }
                        if(lchg){
                            for(int i=0;i<natms;i++){
                                enc[i]*=chge[i];
                                ens[i]*=chge[i];
                            }
                        }
                        for(int m=mmin;m<=kmax;m++){
                            if(m>=0){
                                for(int i=0;i<natms;i++){
                                    emc[i]=ecy[i][m]*enc[i]-esy[i][m]*ens[i];
                                    ems[i]=ecy[i][m]*ens[i]+esy[i][m]*enc[i];
                                }
                            }
                            else{
                                for(int i=0;i<natms;i++){
                                    emc[i]=ecy[i][-m]*enc[i]+esy[i][-m]*ens[i];
                                    ems[i]=ecy[i][-m]*ens[i]-esy[i][-m]*enc[i];
                                }
                            }
                            for(int l=lmin;l<=kmax;l++){
                                if(l>=0){
                                    for(int i=0;i<natms;i++){
                                        elc[i]=ecx[i][l]*emc[i]-esx[i][l]*ems[i];
                                        els[i]=ecx[i][l]*ems[i]+esx[i][l]*emc[i];
                                    }
                                }
                                else{
                                    for(int i=0;i<natms;i++){
                                        elc[i]=ecx[i][-l]*emc[i]+esx[i][-l]*ems[i];
                                        els[i]=ecx[i][-l]*ems[i]-esx[i][-l]*emc[i];
                                    }
                                }
                                dxyz[0][kkk]=0;
                                dxyz[1][kkk]=0;
                                for(int i=0;i<natms;i++){
                                    dxyz[0][kkk]=dxyz[0][kkk]+elc[i];
                                    dxyz[1][kkk]=dxyz[1][kkk]+els[i];
                                }
                                kkk++;
                            }
                            lmin=-kmax;
                        }
                        mmin=-kmax;
                    }

                    //normalise density fourier transform

                    for(int i=0;i<klim;i++){
                        dxyz[0][i]=rnrm*dxyz[0][i];
                        dxyz[1][i]=-rnrm*dxyz[1][i];
                    }

                    //store fourier transform of density

                    try {
                        if(newfile){
                            spc = new DataOutputStream(new FileOutputStream("SPCDEN."+String.valueOf(numskw)));
                            println("File SPCDEN."+String.valueOf(numskw)+" created");
                            if(lchg)
                                spc.writeBytes("Charge density\n");
                            else
                                spc.writeBytes("Particle density\n");
                            newfile=false;
                        }
                        spc.writeBytes(BML.fmt(time,16));
                        for(int i=0;i<klim;i++)
                            spc.writeBytes(BML.fmt(dxyz[0][i],16)+BML.fmt(dxyz[1][i],16));
                        spc.writeBytes("\n");
                    }
                    catch(Exception e) {
                        println("Error - writing file SPCDEN."+String.valueOf(numskw));
                        return false;
                    }
                }
            }
            if(iconf==nconf-1) {
                info[0]=-1.0;
                lnr=hread(fname,name,lnr,info,cell,chge,weight,xyz,vel,frc);
                try{
                    spc.close();
                }
                catch(Exception e){
                    println("Error - writing file SPCDEN."+String.valueOf(numskw));
                    return false;
                }
            }
            println("Number of configurations processed: "+BML.fmt(nn,8));
            nconf=nn;
            return true;
    }
    boolean correl(){
        /*
***********************************************************************

dl_poly java routine to calculate density correlation function

copyright - daresbury laboratory
author    - w smith nov 2007

***********************************************************************
         */
        String record="";
        double rnrm,time,tstep;
        int mm,nn,lor,mor,kkk,lmin,mmin;
        double ckr[][]=new double[2][klim];
        double cfkt[][][]=new double[2][lencor][klim];
        double ckr0[][][]=new double[2][lencor][klim];
        int num[]=new int[lencor];
        int ind[]=new int[lencor];

        nn=0;
        lor=0;
        mor=0;
        time=0;

        //initialise arrays

        for(int i=0;i<lencor;i++)
            num[i]=0;
        for(int i=0;i<klim;i++){
            for(int j=0;j<lencor;j++){
                cfkt[0][j][i]=0;
                cfkt[1][j][i]=0;
            }
        }

        //start of loop over configurations

        try{
            DataOutputStream corr = new DataOutputStream(new FileOutputStream("DENFKT."+String.valueOf(numskw)));
            File spcdata = new File("SPCDEN."+String.valueOf(numskw));
            FileInputStream instream = new FileInputStream(spcdata);
            InputStreamReader isr = new InputStreamReader(instream);
            BufferedReader reader = new BufferedReader(isr);
            println("File SPCDEN."+String.valueOf(numskw)+" opened");
            println("File DENFKT."+String.valueOf(numskw)+" created");
            if(lchg)
                corr.writeBytes("Charge density correlation functions\n");
            else
                corr.writeBytes("Particle density correlation functions\n");
            record=reader.readLine();
            OUT:
                for(int n=0;n<nconf;n++){
                    record=reader.readLine();
                    if(record.equals(null)) break OUT;
                    time=BML.giveDouble(record,1);
                    for(int i=0;i<klim;i++){
                        ckr[0][i]=BML.giveDouble(record,2*i+2);
                        ckr[1][i]=BML.giveDouble(record,2*i+3);
                    }
                    ind[mor]=0;
                    for(int k=0;k<klim;k++){
                        ckr0[0][mor][k]=ckr[0][k];
                        ckr0[1][mor][k]=-ckr[1][k];
                    }
                    lor=Math.min(lor+1,lencor);
                    mor=(mor+1)%lencor;
                    for(int l=0;l<lor;l++){
                        mm=ind[l];
                        ind[l]=mm+1;
                        num[mm]++;
                        for(int k=0;k<klim;k++){
                            cfkt[0][mm][k]+=(ckr0[0][l][k]*ckr[0][k]-ckr0[1][l][k]*ckr[1][k]);
                            cfkt[1][mm][k]+=(ckr0[0][l][k]*ckr[1][k]+ckr0[1][l][k]*ckr[0][k]);
                        }
                    }
                    nn++;
                }
                tstep=time/(nn-1);
                last=Math.min(nn,lencor);
                corr.writeBytes(BML.fmt(klim,10)+BML.fmt(last,10)+"\n");

                //normalise correlation functions

                for(int k=0;k<klim;k++){
                    rnrm=((double)num[0])/cfkt[0][0][k];
                    for(int l=0;l<last;l++){
                        cfkt[0][l][k]*=(rnrm/num[l]);
                        cfkt[1][l][k]*=(rnrm/num[l]);
                    }
                }

                //store correlation functions in DENFKT file

                kkk=0;
                mmin=0;
                lmin=1;
                for(int n=0;n<=kmax;n++){
                    for(int m=mmin;m<=kmax;m++){
                        for(int l=lmin;l<=kmax;l++){
                            corr.writeBytes("k-vector"+BML.fmt(l,5)+BML.fmt(m,5)+BML.fmt(n,5)+"\n");
                            for(int k=0;k<last;k++)
                                corr.writeBytes(BML.fmt(tstep*k,14)+BML.fmt(cfkt[0][k][kkk],14)+BML.fmt(cfkt[1][k][kkk],14)+"\n");
                            kkk++;
                        }
                        lmin=-kmax;
                    }
                    mmin=-kmax;
                }
                println("Number of correlation functions:"+BML.fmt(klim,10));
                corr.close();
                reader.close();
        }
        catch(Exception e){
            println("Error - problem calculating density correlations");
            return false;
        }
        return true;
    }
    boolean denfft(){
        /*
***********************************************************************

dl_poly java routine to calculate fourier transform of density
correlation function

copyright - daresbury laboratory
author    - w smith nov 2007

***********************************************************************
         */

        String record="";
        boolean safe=true;
        int klim=(int)((Math.pow(2*kmax+1,3)-1)/2);
        int nn,lencor2,lencor4,kkk,lmin,mmin;
        double a0,a1,a2,arg,ccc,time,tstep,omega;

        //make sure FFT is a power of 2

        lencor2=2;
        while(lencor2<=lencor){
            lencor2*=2;
        }
        lencor4=2*lencor2;

        //allocate arrays

        int key[]=new int[lencor4];
        double wind[]=new double[lencor];
        double fta[][]=new double[2][lencor4];
        double ftb[][]=new double[2][lencor4];
        double work[][]=new double[2][lencor4];
        double cfkt[][][]=new double[2][lencor][klim];

        //initialise complex fast fourier transform routine

        fft(safe,1,-1,lencor4,key,fta,work,ftb);

        //set up window function (blackman function)

        a0=0.42;
        a1=0.5;
        a2=0.08;
        arg=Math.PI/last;
        for(int i=0;i<last;i++){
            ccc=Math.cos(arg*(i+last));
            wind[i]=a0-a1*ccc+a2*(2*ccc*ccc-1);
        }

        //read the correlation functions from DENFKT file

        time=0;
        tstep=0;

        try{
            File corr = new File("DENFKT."+String.valueOf(numskw));
            FileInputStream instream = new FileInputStream(corr);
            InputStreamReader isr = new InputStreamReader(instream);
            BufferedReader reader = new BufferedReader(isr);
            println("File DENFKT."+String.valueOf(numskw)+" opened");
            record=reader.readLine();
            record=reader.readLine();
            OUT:
                for(int k=0;k<klim;k++){
                    record=reader.readLine();
                    if(record.equals(null)) break OUT;
                    nn=0;
                    for(int i=0;i<last;i++){
                        record=reader.readLine();
                        if(record.equals(null)) break OUT;
                        time=BML.giveDouble(record,1);
                        cfkt[0][i][k]=BML.giveDouble(record,2);
                        cfkt[1][i][k]=BML.giveDouble(record,3);
                        nn++;
                    }
                    tstep=time/(nn-1);
                }
                reader.close();
        }
        catch(Exception e){
            println("Error - problem reading DENFKT file");
            return false;
        }

        //Calculate dynamic structure factors

        for(int k=0;k<klim;k++){

            //apply window function

            for(int i=0;i<last;i++){
                fta[0][i]=wind[i]*cfkt[0][i][k];
                fta[1][i]=wind[i]*cfkt[1][i][k];
            }
            for(int i=last;i<lencor4;i++){
                fta[0][i]=0;
                fta[1][i]=0;
            }

            fta[0][0]/=2;
            fta[1][0]/=2;

            //calculate complex fourier transform

            fft(safe,0,-1,lencor4,key,fta,work,ftb);

            //store dynamic structure factors

            for(int i=0;i<lencor;i++){
                cfkt[0][i][k]=ftb[0][i];
                cfkt[1][i][k]=ftb[1][i];
            }
        }

        //write dynamic structure factors in DENSKW file

        omega=2*Math.PI/(tstep*lencor2);

        try{
            sname="DENSKW."+String.valueOf(numskw);
            DataOutputStream skw = new DataOutputStream(new FileOutputStream(sname));
            println("File DENSKW."+String.valueOf(numskw)+" created");
            if(lchg)
                skw.writeBytes("Charge dynamic structure factors\n");
            else
                skw.writeBytes("Particle dynamic structure factors\n");
            skw.writeBytes(BML.fmt(klim,10)+BML.fmt(lencor/4,10)+"\n");
            kkk=0;
            mmin=0;
            lmin=1;
            for(int n=0;n<=kmax;n++){
                for(int m=mmin;m<=kmax;m++){
                    for(int l=lmin;l<=kmax;l++){
                        skw.writeBytes("k-vector"+BML.fmt(l,5)+BML.fmt(m,5)+BML.fmt(n,5)+"\n");
                        for(int i=0;i<lencor/4;i++){
                            skw.writeBytes(BML.fmt(i*omega,14)+BML.fmt(cfkt[0][i][kkk],14)+BML.fmt(cfkt[1][i][kkk],14)+"\n");
                        }
                        kkk++;
                    }
                    lmin=-kmax;
                }
                mmin=-kmax;
            }
            skw.close();
        }
        catch(Exception e){
            println("Error - problem writing DENSKW file");
            return false;
        }
        println("Number of correlation functions:"+BML.fmt(klim,10));
        return true;
    }
}
