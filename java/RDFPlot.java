import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;

public class RDFPlot extends Basic implements ActionListener {
        /*
*********************************************************************

dl_poly/java GUI class to plot radial distribution functions

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
    public static RDFPlot job;
    private static GUI home;
    private static String name1,name2,title;
    private static double[] xx,yy;
    private static int npnts,numrdf;
    private static JTextField atom1,atom2,rdfdat;
    private static JButton plot,close;

    public RDFPlot() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        setTitle("RDF Plotter");

        getContentPane().setBackground(art.back);
        getContentPane().setForeground(art.fore);
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

        JLabel lab1 = new JLabel("Required RDFDAT file:",JLabel.LEFT);
        fix(lab1,grd,gbc,0,1,3,1);

        // Name of RDFDAT file

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

    public RDFPlot(GUI here) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        home=here;
        println("Activated panel for plotting RDFs");
        job=new RDFPlot();
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
            npnts=rdrdf(fname,name1,name2);
            if(npnts>0) {
                if(graf != null)
                    graf.job.setVisible(false);
                graf=new GraphDraw(home);
                rdfXY(npnts,name1,name2);
                graf.xlabel.setText("Radius (A)");
                graf.ylabel.setText("G(r)");
                graf.plabel.setText("RDF of "+name1.trim()+" - "+name2.trim());
                graf.extraPlot(npnts,xx,yy);
            }
        }
        else if (arg.equals("Close")) {
            job.setVisible(false);
        }
    }
    int rdrdf(String fname,String atnam1,String atnam2) {
        /*
*********************************************************************

dl_poly/java routine to read a DL_POLY RDFDAT file

copyright - daresbury laboratory
author    - w.smith february 2001

*********************************************************************
         */
        int nrdfs,npts;
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
            xx=new double[npts];
            yy=new double[npts];
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
            return -1;
        }
        catch(Exception e) {
            println("Error reading file: " + fname + " "+e);
            return -3;
        }
        println("Number of points loaded:"+BML.fmt(npts,6));
        return npts;

    }
    void rdfXY(int npts,String anam1,String anam2) {
        /*
*********************************************************************

dl_poly/java GUI routine  to create a RDF XY file

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        String fname;
        String[] header;
        int nhead=4,call;

        header=new String[nhead];
        fname="RDF"+String.valueOf(numrdf)+".XY";
        numrdf++;
        header[0]=" RDF Plotting Program";
        header[1]=" RDF Plot: "+anam1.trim()+" + "+anam2.trim();
        header[2]=" Radius (A)";
        header[3]=" RDF";
        call=putXY(fname,header,nhead,npts,xx,yy);
        if(call==0)
            println("PLOT file "+fname+" created");
    }
}

