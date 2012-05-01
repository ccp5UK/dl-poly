import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;

public class StatProp extends Basic implements ActionListener {
    /*
*********************************************************************

dl_poly/java GUI class to calculate statistical properties

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
     */
    public static StatProp job;
    private static GUI home;
    private static String prop;
    private static JComboBox<String> property;
    private static int npnts,keyprp,nstart;
    private static JTextField statis,start,selection;
    private static JButton run,close;
    private static ButtonGroup options;
    private static JCheckBox normal,autocorr,ftrans;
    private static int[] key={0,1,2,3,4,18,26};
    private static double[] xx,yy;
    private static double statmean;

    // Define the Graphical User Interface

    public StatProp() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        super();
        setTitle("Statistics Panel");
        int n=0;

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
        fix(run,grd,gbc,0,n,1,1);

        fix(new JLabel("  "),grd,gbc,1,n,1,1);

        // Define the Close button

        close = new JButton("Close");
        close.setBackground(art.butn);
        close.setForeground(art.butf);
        fix(close,grd,gbc,2,n++,1,1);

        // Name of STATIS file

        fix(new JLabel("Required STATIS file:",JLabel.LEFT),grd,gbc,0,n++,3,1);
        statis = new JTextField(12);
        statis.setBackground(art.scrn);
        statis.setForeground(art.scrf);
        fix(statis,grd,gbc,0,n++,3,1);

        // Choice of ensemble

        fix(new JLabel(" Property:   ",JLabel.RIGHT),grd,gbc,0,n,1,1);
        property = new JComboBox<String>();
        property.setBackground(art.scrn);
        property.setForeground(art.scrf);
        property.addItem("E-TOT");
        property.addItem("TEMP");
        property.addItem("E-CFG");
        property.addItem("E-VDW");
        property.addItem("E-COUL");
        property.addItem("VOLUME");
        property.addItem("PRESS");
        property.addItem("SELECT");
        fix(property,grd,gbc,2,n++,1,1);

        // User property selection

        selection = new JTextField(6);
        selection.setBackground(art.scrn);
        selection.setForeground(art.scrf);
        fix(selection,grd,gbc,2,n,1,1);
        fix(new JLabel("Selection",JLabel.LEFT),grd,gbc,0,n++,2,1);

        // Type of plot required

        options=new ButtonGroup();
        normal=new JCheckBox("Normal Plot",true);
        normal.setForeground(art.fore);
        normal.setBackground(art.back);
        options.add(normal);
        fix(normal,grd,gbc,0,n++,3,1);
        autocorr=new JCheckBox("Autocorrelate",false);
        autocorr.setForeground(art.fore);
        autocorr.setBackground(art.back);
        options.add(autocorr);
        fix(autocorr,grd,gbc,0,n++,3,1);
        ftrans=new JCheckBox("Fourier Transform",false);
        ftrans.setForeground(art.fore);
        ftrans.setBackground(art.back);
        options.add(ftrans);
        fix(ftrans,grd,gbc,0,n++,3,1);

        // Start time step

        start = new JTextField(6);
        start.setBackground(art.scrn);
        start.setForeground(art.scrf);
        fix(start,grd,gbc,0,n,1,1);
        fix(new JLabel("Start time step",JLabel.LEFT),grd,gbc,1,n++,2,1);


        // Register action buttons

        run.addActionListener(this);
        close.addActionListener(this);
        property.addActionListener(this);
    }

    public StatProp(GUI here) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        home=here;
        println("Activated statistics panel");
        job=new StatProp();
        job.pack();
        job.setVisible(true);
        npnts=0;
        keyprp=0;
        nstart=0;
        fname="STATIS";
        statis.setText(fname);
        start.setText(String.valueOf(nstart));
        property.setSelectedIndex(keyprp);
        selection.setText(String.valueOf(key[keyprp]));
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
        String[] header=new String[4];
        boolean safe=true;
        double a0,a1,a2,del,ccc;
        int nfft,select;

        statmean=0.0;
        header[0]=" Statistics Plotting Program";

        if (e.getSource() instanceof JComboBox) {
            select=property.getSelectedIndex();
            if(select < 7)
                selection.setText(String.valueOf(key[select]));
            else
                selection.setText("Enter: ");
        }
        else {
            if (arg.equals("Run")) {
                fname=statis.getText();
                prop=property.getSelectedItem().toString();
                select=property.getSelectedIndex();
                if(select == 7) {
                    keyprp=BML.giveInteger(selection.getText(),2);
                    selection.setText(String.valueOf(keyprp));
                }
                else
                    keyprp=key[select];
                nstart=BML.giveInteger(start.getText(),1);
                println("Started statistics calculation .....");
                // calcStats resets variable statmean
                npnts=calcStats();
                if(npnts>0) {
                    if(normal.isSelected()) {

                        // Normal plot of statistical data

                        header[1]=" Plot of "+prop;
                        header[2]=" Time (ps)";
                        header[3]=" "+prop;
                        statsXY(npnts,header);
                        graf.xlabel.setText(header[2].substring(1));
                        graf.ylabel.setText(header[3].substring(1));
                        graf.plabel.setText(header[1].substring(1));
                        graf.extraPlot(npnts,xx,yy);
                    }
                    else {
                        // determine order of required FFT
                        nfft=1;
                        while(nfft < npnts)
                            nfft*=2;

                        if(autocorr.isSelected()) {

                            // Plot of autocorrelation of statistical data

                            nfft*=2;
                            int key[]=new int[nfft];
                            double fta[][]=new double[2][nfft];
                            double ftb[][]=new double[2][nfft];
                            double work[][]=new double[2][nfft];
                            // FFT initialisation
                            fft(safe,1,-1,nfft,key,fta,work,ftb);
                            // load function for FFT
                            for(int i=0;i<nfft;i++){
                                fta[1][i]=0.0;
                                if(i<npnts)
                                    fta[0][i]=yy[i]-statmean;
                                else
                                    fta[0][i]=0.0;
                            }
                            // calculate Fourier transform
                            fft(safe,0,-1,nfft,key,fta,work,ftb);
                            // calculate Fourier transform of correlation function
                            for(int i=0;i<nfft;i++) {
                                fta[0][i]=(Math.pow(ftb[0][i],2)+Math.pow(ftb[1][i],2))/((double)nfft);
                                fta[1][i]=0.0;
                            }
                            // calculate inverse Fourier transform
                            fft(safe,0,1,nfft,key,fta,work,ftb);
                            // calculate correlation function
                            del=(xx[npnts-1]-xx[0])/(double)(npnts-1);
                            for(int i=0;i<npnts;i++) {
                                xx[i]=del*(double)i;
                                yy[i]=(Math.pow(ftb[0][i],2)+Math.pow(ftb[1][i],2))*del/((double)(npnts-i));
                            }
                            header[1]=" Autocorr. of "+prop+" fluctuations";
                            header[2]=" Time (ps)";
                            header[3]=" C(t)";
                            statsXY(npnts,header);
                            graf.xlabel.setText(header[2].substring(1));
                            graf.ylabel.setText(header[3].substring(1));
                            graf.plabel.setText(header[1].substring(1));
                            graf.extraPlot(npnts,xx,yy);
                        }
                        else if(ftrans.isSelected()) {

                            // Plot of Fourier transform of statistical data

                            int key[]=new int[nfft];
                            double fta[][]=new double[2][nfft];
                            double ftb[][]=new double[2][nfft];
                            double work[][]=new double[2][nfft];
                            // FFT initialisation
                            fft(safe,1,-1,nfft,key,fta,work,ftb);
                            // load function for FFT (use Blackman-Harris window)
                            a0=0.42;
                            a1=0.5;
                            a2=0.08;
                            del=2.0*Math.PI/npnts;
                            for(int i=0;i<nfft;i++){
                                fta[1][i]=0.0;
                                if(i<npnts) {
                                    ccc=Math.cos(del*i);
                                    fta[0][i]=(yy[i]-statmean)*(a0-a1*ccc+a2*(2*ccc*ccc-1));
                                }
                                else
                                    fta[0][i]=0.0;
                            }
                            // calculate Fourier transform
                            fft(safe,0,-1,nfft,key,fta,work,ftb);
                            // calculate real amplitudes
                            del=2.0*Math.PI*((double)(npnts-1)/(double)nfft)/(xx[npnts-1]-xx[0]);
                            nfft/=2;
                            for(int i=0;i<nfft;i++) {
                                xx[i]=del*(double)i;
                                yy[i]=Math.sqrt(Math.pow(ftb[0][i],2)+Math.pow(ftb[1][i],2));
                            }
                            header[1]=" F. Trans. of "+prop+" fluctuations";
                            header[2]=" Omega (THz)";
                            header[3]=" Amplitude";
                            statsXY(nfft,header);
                            graf.xlabel.setText(header[2].substring(1));
                            graf.ylabel.setText(header[3].substring(1));
                            graf.plabel.setText(header[1].substring(1));
                            graf.extraPlot(nfft,xx,yy);
                        }
                    }
                }
            }
            else if (arg.equals("Close")) {
                job.dispose();
            }
        }
    }

    int calcStats() {
        /*
**********************************************************************

dl_poly/java routine for blocking method for estimating standard
error of mean (Flyvbjerg and Petersen  JCP 91 (1989) 461)

copyright - daresbury laboratory
author    - w.smith march 2001


**********************************************************************
         */

        LineNumberReader lnr=null;
        int ndeg=128;
        int i,j,n,mmm,mcols,mxpnts;
        double aver;
        String record;

        mxpnts=250;
        xx=new double[mxpnts];
        yy=new double[mxpnts];
        double alpha[]=new double[mxpnts];
        double sigma[]=new double[mxpnts];
        double enigma[]=new double[mxpnts];

        // open the statistical data file

        try {
            lnr = new LineNumberReader(new FileReader(fname));
            println("Reading file: "+fname);
            record = lnr.readLine();
            println("File header record: "+record);
            record = lnr.readLine();
            println(record);

            // read in statistical data

            i=0;
            while((record=lnr.readLine())!= null) {
                mcols=BML.giveInteger(record,3);
                if(i==mxpnts) {
                    mxpnts*=2;
                    double aaa[]=new double[mxpnts];
                    double bbb[]=new double[mxpnts];
                    double ccc[]=new double[mxpnts];
                    double ddd[]=new double[mxpnts];
                    double eee[]=new double[mxpnts];
                    System.arraycopy(xx,0,aaa,0,i);
                    System.arraycopy(yy,0,bbb,0,i);
                    System.arraycopy(alpha,0,ccc,0,i);
                    System.arraycopy(sigma,0,ddd,0,i);
                    System.arraycopy(enigma,0,eee,0,i);
                    enigma=eee;
                    sigma=ddd;
                    alpha=ccc;
                    xx=aaa;
                    yy=bbb;
                }
                if(BML.giveInteger(record,1) >= nstart) {
                    xx[i]=BML.giveDouble(record,2);
                    for(j=0;j<=(mcols-1)/5;j++) {
                        record=lnr.readLine();
                        if(j==keyprp/5) {
                            yy[i]=BML.giveDouble(record,keyprp-5*j+1);
                            alpha[i]=BML.giveDouble(record,keyprp-5*j+1);
                        }
                    }
                    i++;
                }
                else {
                    for(j=0;j<=(mcols-1)/5;j++)
                        record=lnr.readLine();
                }
            }
            lnr.close();
        }
        catch(FileNotFoundException e) {
            println("Error - STATIS file not found");
            return -2;
        }
        catch(IOException e) {
            println("Error - STATIS file processing error");
            try{lnr.close();}catch(IOException ee){}
            return -3;
        }
        npnts=i;
        println("Number of data points processed: "+BML.fmt(npnts,8));

        // calculate average

        aver=0.0;
        for(i=0;i<npnts;i++)
            aver+=alpha[i];
        aver=aver/npnts;
        statmean=aver;

        println("Calculated average: "+BML.fmt(aver,15));

        // calculate error using blocking method

        println("Standard error and error uncertainty");

        sigma[0]=0.0;
        for(i=0;i<npnts;i++) {
            alpha[i]-=aver;
            sigma[0]+=(alpha[i]*alpha[i]);
        }
        sigma[0]=Math.sqrt(sigma[0]/(npnts*(npnts-1)));
        enigma[0]=sigma[0]/Math.sqrt(2*(npnts-1));
        println(BML.fmt(sigma[0],15)+"+/-"+BML.fmt(enigma[0],15));

        n=0;
        mmm=npnts;
        while(mmm>3 && n<ndeg) {
            j=0;
            n++;
            mmm=mmm/2;
            sigma[n]=0.0;
            for(i=0;i<mmm;i++) {
                alpha[i]=0.5*(alpha[j]+alpha[j+1]);
                sigma[n]+=(alpha[i]*alpha[i]);
                j+=2;
            }
            sigma[n]=Math.sqrt(sigma[n]/(mmm*(mmm-1)));
            enigma[n]=sigma[n]/Math.sqrt((2*(mmm-1)));
            println(BML.fmt(sigma[n],15)+"+/-"+BML.fmt(enigma[n],15));
        }

        return npnts;
    }
    void statsXY(int npts, String[] header) {
        /*
*********************************************************************

dl_poly/java GUI routine to create a STATS XY file

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        String fname;

        fname="STATS"+String.valueOf(numstat++)+".XY";
        if(putXY(fname,header,4,npts,xx,yy)==0)
            println("PLOT file "+fname+" created");
    }
}
