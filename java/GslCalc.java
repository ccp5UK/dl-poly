import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;

public class GslCalc extends Basic implements ActionListener {
        /*
*********************************************************************

dl_poly/java GUI class to calculate van Hove self correlation

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
    public static GslCalc job;
    private static GUI home;
    private static HovePlot hovplt=null;
    private static double tstep,rcut,delr;
    private static String atname;
    private static int nconf,lencor,isampl,iorig,npnts,mxrad;
    private static JTextField atom,history,configs,length,sample,origin,cutoff;
    private static JCheckBox format;
    private static boolean form;
    private static JButton run,close,plot;
    private static int[] imd,msm;
    private static String[] name;
    private static double[] cell,chge,weight;
    private static double[][] acm,xy0,xy1,xyz,vel,frc,gslf;
    private static double[][][] gslf0;
    private static String hname=null;

    // Define the Graphical User Interface

    public GslCalc() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        setTitle("Gs(r,t) Calculator");

        getContentPane().setBackground(art.back);
        getContentPane().setForeground(art.fore);
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

        // History file format

        JLabel lab2 = new JLabel("File is formatted?",JLabel.LEFT);
        //fix(lab2,grd,gbc,0,3,2,1);
        format=new JCheckBox("    ");
        format.setBackground(art.back);
        format.setForeground(art.fore);
        //fix(format,grd,gbc,2,3,1,1);

        // Name of first atom type

        JLabel lab3 = new JLabel("Atom name:",JLabel.LEFT);
        fix(lab3,grd,gbc,0,4,2,1);
        atom = new JTextField(8);
        atom.setBackground(art.scrn);
        atom.setForeground(art.scrf);
        fix(atom,grd,gbc,2,4,1,1);

        // Number of configurations

        JLabel lab4 = new JLabel("No. configurations:",JLabel.LEFT);
        fix(lab4,grd,gbc,0,5,2,1);
        configs = new JTextField(8);
        configs.setBackground(art.scrn);
        configs.setForeground(art.scrf);
        fix(configs,grd,gbc,2,5,1,1);

        // GSL array lengths

        JLabel lab5 = new JLabel("GSL array lengths:",JLabel.LEFT);
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

        // Cutoff radius

        JLabel lab8 = new JLabel("Cutoff radius (A):",JLabel.LEFT);
        fix(lab8,grd,gbc,0,9,2,1);
        cutoff = new JTextField(8);
        cutoff.setBackground(art.scrn);
        cutoff.setForeground(art.scrf);
        fix(cutoff,grd,gbc,2,9,1,1);


        // Define the Close button

        close = new JButton("Close");
        close.setBackground(art.butn);
        close.setForeground(art.butf);
        fix(close,grd,gbc,0,10,1,1);
        fix(new JLabel("  "),grd,gbc,1,10,1,1);

        // Define the Plot button

        plot = new JButton("Plot");
        plot.setBackground(art.butn);
        plot.setForeground(art.butf);
        fix(plot,grd,gbc,2,10,1,1);

        // Register action buttons

        run.addActionListener(this);
        close.addActionListener(this);
        plot.addActionListener(this);

    }

    public GslCalc(GUI here) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        home=here;
        println("Activated panel for calculating GSLs");
        job=new GslCalc();
        job.pack();
        job.setVisible(true);
        atname="ALL";
        fname="HISTORY";
        nconf=1000;
        lencor=512;
        isampl=1;
        iorig=1;
        rcut=3.5;
        form=true;
        atom.setText(atname);
        history.setText(fname);
        format.setSelected(form);
        configs.setText(String.valueOf(nconf));
        length.setText(String.valueOf(lencor));
        sample.setText(String.valueOf(isampl));
        origin.setText(String.valueOf(iorig));
        cutoff.setText(String.valueOf(rcut));
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
            atname=atom.getText();
            fname=history.getText();
            form=format.isSelected();
            nconf=BML.giveInteger(configs.getText(),1);
            lencor=BML.giveInteger(length.getText(),1);
            isampl=BML.giveInteger(sample.getText(),1);
            iorig=BML.giveInteger(origin.getText(),1);
            rcut=BML.giveDouble(cutoff.getText(),1);
            println("Started Gs(r,t) calculation");
            npnts=gself();
            if(npnts>0) gslFile();
        }
        else if (arg.equals("Close")) {
            job.setVisible(false);
        }
        else if (arg.equals("Plot")) {
            if(hovplt != null)
                hovplt.job.setVisible(false);
            hovplt=new HovePlot(home,hname);
        }
    }
    int gself() {
        /*
*********************************************************************

dl_poly/java routine to calculate van Hove self correlation
function for selected atoms from dl_poly HISTORY file

copyright - daresbury laboratory
author    - w.smith march 2001

*********************************************************************
         */

        boolean all;
        int nsgslf,nat,npts,nogslf,lsr,msr,iconf,imcon,natms,k,n,m;
        double rcut2,f1,f2,uuu,vvv,www,rmsx,rmsy,rmsz,rsq,tcut;
	double bcell[];
        LineNumberReader lnr=null;
        double cell[]=new double[9];
        double rcell[]=new double[9];
        double avcell[]=new double[9];
        double info[]=new double[10];

        nat=0;
        npts=0;
        all=false;
        tstep=0.0;
        if(atname.toUpperCase().equals("ALL"))all=true;

        // check on specified control variables

        if(lencor%iorig != 0) {
            lencor=iorig*(lencor/iorig);
            println("Warning - 1st dimension of Gself array reset to "+BML.fmt(lencor,8));
        }

        mxrad=Math.max(64,lencor/4);
        rcut2=rcut*rcut;
        delr=rcut/mxrad;
        nogslf=lencor/iorig;
        imd=new int[lencor];
        msm=new int[lencor];

        // write control variables

        println("Name of target HISTORY file   : "+fname);
        println("Label  of atom  of interest   : "+atname);
        println("Length of correlation arrays  : "+BML.fmt(lencor,8));
        println("Number of configurations      : "+BML.fmt(nconf,8));
        println("Sampling interval             : "+BML.fmt(isampl,8));
        println("Interval between origins      : "+BML.fmt(iorig,8));
        println("Required correlation radius   : "+BML.fmt(rcut,10));
        println("Radial correlation bin width  : "+BML.fmt(delr,10));

        // initialize average cell vectors

        for(int i=0;i<9;i++)
            avcell[i]=0.0;

        // initialise correlation arrays

        gslf=new double[lencor][mxrad];

        for(int j=0;j<lencor;j++) {
            msm[j]=0;
            for(int i=0;i<mxrad;i++)
                gslf[j][i]=0.0;
        }

        // process the HISTORY file data

        for(int ipass=0;ipass<isampl;ipass++) {
            lsr=0;
            msr=-1;
            nsgslf=0;

            // initialise control parameters for HISTORY file reader

            info[0]=0.0;
            info[1]=999999;
            info[2]=0.0;
            info[3]=0.0;
            info[4]=0.0;
            info[5]=0.0;
            if(form) {
                lnr=hread(fname,name,lnr,info,cell,chge,weight,xyz,vel,frc);
                if(BML.nint(info[3])<0 && BML.nint(info[3])!=-1) {
                    println("Error - HISTORY file data error");
                    return -1;
                }
            }
            else {
                println("Error - unformatted read option not active");
                return -2;
            }
            natms=BML.nint(info[7]);
	    imcon=BML.nint(info[5]);
	    if(imcon < 1 || imcon > 3){
		println("Error - incorrect periodic boundary condition");
		info[0]=-1.0;
		lnr=hread(fname,name,lnr,info,cell,chge,weight,xyz,vel,frc);
		return -4;
	    }

            // initialise gself arrays

            if(ipass==0) {
                name=new String[natms];
                chge=new double[natms];
                weight=new double[natms];
                acm=new double[3][natms];
                xy0=new double[3][natms];
                xy1=new double[3][natms];
                xyz=new double[3][natms];
                gslf0=new double[lencor][natms][3];
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
		if(lnr == null) break OUT;
		if(iconf==0)info[9]=info[6];
		if(iconf==1)tstep=info[8]*(info[6]-info[9]);
		if(BML.nint(info[3])==-1)break OUT;

		if(iconf == 0){
		    bcell=AML.dcell(cell);
		    tcut=0.5*BML.min(bcell[6],bcell[7],bcell[8]);
		    if(rcut > tcut){
			println("Warning - cut off reset to "+BML.fmt(tcut,10));
			rcut=tcut;
			rcut2=rcut*rcut;
			delr=rcut/mxrad;
		    }
		}

		// select relevant atoms for calculation

		n=0;
		rcell=AML.invert(cell);
		if(iconf==0) {
		    for(int i=0;i<natms;i++) {
			acm[0][i]=0.0;
			acm[1][i]=0.0;
			acm[2][i]=0.0;
		    }
		}
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
		if(ipass==0 && iconf==0)
		    println("Number of atoms of selected type : "+BML.fmt(nat,8));
		if(nat == 0) {
		    println("Error - zero atoms of specified type");
		    info[0]=-1.0;
		    lnr=hread(fname,name,lnr,info,cell,chge,weight,xyz,vel,frc);
		    return -6;
		}

		// running average of cell vectors

		f1=((double)iconf)/(iconf+1);
		f2=1.0/(iconf+1);
		for(int i=0;i<9;i++)
		    avcell[i]=f1*avcell[i]+f2*cell[i];

		if(iconf>ipass) {

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

		// calculate self correlation function

		if(iconf>ipass) {
		    if(iconf%isampl == ipass) {
			if(nsgslf%iorig==0) {
			    lsr=Math.min(lsr+1,nogslf);
			    msr=(msr+1)%nogslf;
			    imd[msr]=0;
			    for(int i=0;i<nat;i++) {
				gslf0[msr][i][0]=0.0;
				gslf0[msr][i][1]=0.0;
				gslf0[msr][i][2]=0.0;
			    }
			}
			nsgslf++;
			for(int j=0;j<lsr;j++) {
			    m=imd[j];
			    imd[j]=m+1;
			    msm[m]++;

			    for(int i=0;i<nat;i++) {
				rmsx=gslf0[j][i][0]+acm[0][i];
				rmsy=gslf0[j][i][1]+acm[1][i];
				rmsz=gslf0[j][i][2]+acm[2][i];
				rsq=rmsx*rmsx+rmsy*rmsy+rmsz*rmsz;
				gslf0[j][i][0]=rmsx;
				gslf0[j][i][1]=rmsy;
				gslf0[j][i][2]=rmsz;
				if(rsq < rcut2) {
				    k=(int)(Math.sqrt(rsq)/delr);
				    gslf[m][k]+=1.0;
				}
			    }
			}
			for(int i=0;i<nat;i++) {
			    acm[0][i]=0.0;
			    acm[1][i]=0.0;
			    acm[2][i]=0.0;
			}
		    }
		}
	    }
	    if(iconf==nconf-1) {
		info[0]=-1.0;
		lnr=hread(fname,name,lnr,info,cell,chge,weight,xyz,vel,frc);
	    }
	    if(BML.nint(info[3])==-1) iconf--;
	    npts=Math.min(nogslf,nsgslf);
	    if(ipass==0)
		println("Number of configurations read: "+BML.fmt(iconf,8));
        }

        // normalise self correlation function

        for(int j=0;j<npts;j++) {
            for(int i=0;i<mxrad;i++) {
                gslf[j][i]=(gslf[j][i]/(msm[j]*nat))/(4.0*Math.PI*Math.pow(delr,3)*(Math.pow((i-0.5),2)+1.0/12.0));
            }
        }
        return npts;
    }
    void gslFile() {
        /*
*********************************************************************

dl_poly/java GUI routine to create a Gself correlation data file

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        double time,rrr;
        hname="HOVGSL."+String.valueOf(numgsl);
        try {
            DataOutputStream out = new DataOutputStream(new FileOutputStream(hname));
            out.writeBytes("Gs(r,t) functions for atom "+atname+"\n");
            out.writeBytes(BML.fmt(atname,8)+BML.fmt(rcut,6)+"A\n");
            out.writeBytes(BML.fmt(npnts,10)+BML.fmt(mxrad,10)+"\n");

            // print out final self correlation functions

            tstep*=isampl;

            for(int j=0;j<npnts;j++) {
                time=tstep*(j+1);
                out.writeBytes(BML.fmt(time,8)+BML.fmt(time,14)+"\n");
                for(int i=0;i<mxrad;i++) {
                    rrr=delr*(i+0.5);
                    out.writeBytes(BML.fmt(rrr,14)+BML.fmt(gslf[j][i],14)+"\n");
                }
            }
            out.close();
        }
        catch(Exception e) {
            println("Error - file HOVGSL." + numgsl);
        }
        println("Number of functions created :"+BML.fmt(npnts,8));
        println("Gself file HOVGSL."+numgsl+" created");
        numgsl++;
    }

}
