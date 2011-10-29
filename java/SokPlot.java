import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;

public class SokPlot extends Basic implements ActionListener {
        /*
*********************************************************************

dl_poly/java GUI class to plot structure factots

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
    public static SokPlot job;
    private static GUI home;
    private static String name1,name2,title;
    private static int npnts;
    private static double[] xx,yy;
    private static JTextField atom1,atom2,rdfdat;
    private static JButton plot,close;

    // Define the Graphical User Interface

    public SokPlot() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        setTitle("S(k) Plotter");

        getContentPane().setBackground(art.back);
        getContentPane().setForeground(art.fore);
        setDefaultCloseOperation(DISPOSE_ON_CLOSE);
        setFont(fontMain);
        GridBagLayout grd = new GridBagLayout();
        GridBagConstraints gbc = new GridBagConstraints();
        getContentPane().setLayout(grd);

        gbc.fill=GridBagConstraints.BOTH;

        // Define the Plot button

        plot = new JButton("Plot");
        plot.setBackground(art.butn);
        plot.setForeground(art.butf);
        fix(plot,grd,gbc,0,0,1,1);

        fix(new JLabel("  "),grd,gbc,1,0,1,1);

        // Define the Close button

        close = new JButton("Close");
        close.setBackground(art.butn);
        close.setForeground(art.butf);
        fix(close,grd,gbc,2,0,1,1);

        // Instruction label 1

        JLabel lab1 = new JLabel("Input file:",JLabel.LEFT);
        fix(lab1,grd,gbc,0,1,3,1);

        // Name of first atom type

        rdfdat = new JTextField(18);
        rdfdat.setBackground(art.scrn);
        rdfdat.setForeground(art.scrf);
        fix(rdfdat,grd,gbc,0,2,3,1);

        // Instruction label 2

        JLabel lab2 = new JLabel("Atomic names:",JLabel.LEFT);
        fix(lab2,grd,gbc,0,3,3,1);

        // Name of first atom type

        atom1 = new JTextField(8);
        atom1.setBackground(art.scrn);
        atom1.setForeground(art.scrf);
        fix(atom1,grd,gbc,0,4,1,1);

        // Name of second atom type

        atom2 = new JTextField(8);
        atom2.setBackground(art.scrn);
        atom2.setForeground(art.scrf);
        fix(atom2,grd,gbc,2,4,1,1);

        // Register action buttons

        plot.addActionListener(this);
        close.addActionListener(this);

    }

    public SokPlot(GUI here) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        home=here;
        println("Activated panel for plotting S(k)");
        job=new SokPlot();
        job.pack();
        job.setVisible(true);
        npnts=0;
        name1="Name1";
        name2="Name2";
        fname="RDFDAT";
        atom1.setText(name1);
        atom2.setText(name2);
        rdfdat.setText(fname);
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
        if (arg.equals("Plot")) {
            name1=atom1.getText();
            name2=atom2.getText();
            fname=rdfdat.getText();
            npnts=rdsok(fname,name1,name2);
            if(npnts>0) {
                if(graf != null)
                    graf.job.dispose();
                graf=new GraphDraw(home);
                sokXY(npnts,name1,name2);
                graf.xlabel.setText("k (1/A)");
                graf.ylabel.setText("S(k)");
                graf.plabel.setText("S(k) of "+name1.trim()+" - "+name2.trim());
                graf.extraPlot(npnts,xx,yy);
            }
        }
        else if (arg.equals("Close")) {
            job.dispose();
        }
    }
    int rdsok(String fname,String atnam1,String atnam2) {
        /*
*********************************************************************

dl_poly/java routine to read a DL_POLY RDFDAT file

copyright - daresbury laboratory
author    - w.smith february 2001

*********************************************************************
         */
        int nrdfs,npts;
        double delr;
        boolean found=false;
        String record,anam1,anam2;

        try {
            LineNumberReader lnr = new LineNumberReader(new FileReader(fname));
            println("Reading file: "+fname);
            title = lnr.readLine();
            println("File header record: "+title);
            record = lnr.readLine();
            nrdfs=BML.giveInteger(record,1);
            npts=BML.giveInteger(record,2);
            xx=new double[4*npts];
            yy=new double[4*npts];
            OUT:
                for(int n=0;n<nrdfs;n++) {
                    record=lnr.readLine();
                    anam1=BML.giveWord(record,1);
                    anam2=BML.giveWord(record,2);

                    for(int i=0;i<npts;i++) {
                        record=lnr.readLine();
                        xx[i]=BML.giveDouble(record,1);
                        yy[i]=BML.giveDouble(record,2);
                    }

                    if((atnam1.equals(anam1) && atnam2.equals(anam2)) ||
                    (atnam1.equals(anam2) && atnam2.equals(anam1))) {
                        found=true;
                        break OUT;
                    }
                }
                if(!found) {
                    println("Error - required RDF not found in file "+fname);
                    lnr.close();
                    return -1;
                }
                lnr.close();
        }
        catch(FileNotFoundException e) {
            println("Error - file not found: " + fname);
            return -2;
        }
        catch(Exception e) {
            println("Error reading file: " + fname + " "+e);
            return -3;
        }
        println("Number of points loaded:"+BML.fmt(npts,6));

        // Calculate radial Fourier transform of RDF

        npts--;
        for(int i=0;i<npts;i++) {
            xx[i]=0.5*(xx[i]+xx[i+1]);
            yy[i]=0.5*(yy[i]+yy[i+1])-1.0;
        }
        yy[npts-1]*=0.5;
        for(int i=npts;i<4*npts;i++) {
            yy[i]=0.0;
        }
        delr=(xx[npts-1]-xx[0])/(npts-1);
        npts*=4;
        radfft(1,npts,delr,xx,yy);
        return npts;
    }
    void sokXY(int npts,String anam1,String anam2) {
        /*
*********************************************************************

dl_poly/java GUI routine to create a S(k) XY file

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        String fname;
        String[] header;
        int nhead=4,call;

        header=new String[nhead];
        fname="SOK"+String.valueOf(numsok)+".XY";
        numsok++;
        header[0]=" S(k) Plotting Program";
        header[1]=" S(k) Plot: "+anam1.trim()+" + "+anam2.trim();
        header[2]=" Radius (A)";
        header[3]=" S(k)";
        call=putXY(fname,header,nhead,npts,xx,yy);
        if(call==0)
            println("PLOT file "+fname+" created");
    }

    void radfft(int isw,int nnn,double delr,double aaa[],double bbb[]) {
        /*
***********************************************************************

dl_poly/java 3D radial fourier transform routine using lado's method
reference: j. comput. phys. 8 (1971) 417

copyright - daresbury laboratory
author    - w.smith february 2001

note: first data point is delr not 0

***********************************************************************
         */
        double sw,delk;

        // perform fourier transform

        sw=isw*Math.PI/nnn;

        for(int j=0;j<nnn;j++) {
            aaa[j]=0.0;
            for(int i=0;i<nnn;i++) {
                aaa[j]=(i+1)*bbb[i]*Math.sin(sw*(j+1)*(i+1))+aaa[j];
            }
            aaa[j]=(4.0*nnn*Math.pow(delr,3)/(j+1))*aaa[j];
        }
        delk=Math.PI/(delr*nnn);
        for(int i=0;i<nnn;i++) {
            bbb[i]=aaa[i];
            aaa[i]=delk*(i+1);
        }
    }
}

