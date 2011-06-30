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
    private static JComboBox property;
    private static int npnts,keyprp;
    private static JTextField statis;
    private static JButton run,close;
    private static int[] key={0,1,2,3,4,18,26,25};
    private static double[] xx,yy;

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

        // Define the Close button

        close = new JButton("Close");
        close.setBackground(art.butn);
        close.setForeground(art.butf);
        fix(close,grd,gbc,2,0,1,1);

        // Name of STATIS file

        JLabel lab1 = new JLabel("Required STATIS file:",JLabel.LEFT);
        fix(lab1,grd,gbc,0,1,3,1);
        statis = new JTextField(12);
        statis.setBackground(art.scrn);
        statis.setForeground(art.scrf);
        fix(statis,grd,gbc,0,2,3,1);

        // Choice of ensemble

        JLabel lab2 = new JLabel(" Property:   ",JLabel.RIGHT);
        fix(lab2,grd,gbc,0,3,1,1);
        property = new JComboBox();
        property.setBackground(art.scrn);
        property.setForeground(art.scrf);
        property.addItem("E-TOT");
        property.addItem("TEMP");
        property.addItem("E-CFG");
        property.addItem("E-VDW");
        property.addItem("E-COUL");
        property.addItem("VOLUME");
        property.addItem("PRESS");
        property.addItem("PMF");
        fix(property,grd,gbc,2,3,1,1);

        // Register action buttons

        run.addActionListener(this);
        close.addActionListener(this);

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
        fname="STATIS";
        statis.setText(fname);
        property.setSelectedIndex(keyprp);
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
            fname=statis.getText();
            prop=property.getSelectedItem().toString();
            keyprp=key[property.getSelectedIndex()];
            println("Started statistics calculation .....");
            npnts=calcStats();
            if(npnts>0) {
                statsXY(npnts,prop);
                graf.xlabel.setText("Time (ps)");
                graf.ylabel.setText(prop);
                graf.plabel.setText("Plot of "+prop);
                graf.extraPlot(npnts,xx,yy);
            }
        }
        else if (arg.equals("Close")) {
            job.setVisible(false);
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
    void statsXY(int npts,String prop) {
        /*
*********************************************************************

dl_poly/java GUI routine to create a STATS XY file

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        String fname;
        String[] header;
        int nhead=4,call;

        header=new String[nhead];
        fname="STATS"+String.valueOf(numstat)+".XY";
        numstat++;
        header[0]=" Statistics Plotting Program";
        header[1]=" Plot of "+prop;
        header[2]=" Time (ps)";
        header[3]=" "+prop;
        call=putXY(fname,header,nhead,npts,xx,yy);
        if(call==0)
            println("PLOT file "+fname+" created");
    }
}
