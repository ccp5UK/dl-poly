import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;

public class ZdenPlot extends Basic implements ActionListener {
        /*
*********************************************************************

dl_poly/java GUI class to plot z-densities

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
    public static ZdenPlot job;
    private static GUI home;
    private static String name,title;
    private static int npnts;
    private static JTextField atom,zdndat;
    private static JButton plot,close;
    private static double[] xx,yy;

    // Define the Graphical User Interface

    public ZdenPlot() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        setTitle("Z-Density Plotter");

        getContentPane().setBackground(art.back);
        getContentPane().setForeground(art.fore);
	setDefaultCloseOperation(DISPOSE_ON_CLOSE);
        setFont(fontMain);
        GridBagLayout grd = new GridBagLayout();
        GridBagConstraints gbc = new GridBagConstraints();
        getContentPane().setLayout(grd);

        gbc.fill=GridBagConstraints.BOTH;

        // Define the Plot button

        plot = new JButton(" Plot ");
        plot.setBackground(art.butn);
        plot.setForeground(art.butf);
        fix(plot,grd,gbc,0,0,1,1);

        fix(new JLabel("  "),grd,gbc,1,0,1,1);

        // Define the Close button

        close = new JButton("Close");
        close.setBackground(art.butn);
        close.setForeground(art.butf);
        fix(close,grd,gbc,2,0,1,1);

        // Name of ZDNDAT file

        JLabel lab1 = new JLabel("Required ZDNDAT file:",JLabel.LEFT);
        fix(lab1,grd,gbc,0,1,3,1);
        zdndat = new JTextField(12);
        zdndat.setBackground(art.scrn);
        zdndat.setForeground(art.scrf);
        fix(zdndat,grd,gbc,0,2,3,1);

        // Name of atom

        JLabel lab2 = new JLabel("Atom name:  ",JLabel.LEFT);
        fix(lab2,grd,gbc,0,3,1,1);
        atom = new JTextField(8);
        atom.setBackground(art.scrn);
        atom.setForeground(art.scrf);
        fix(atom,grd,gbc,2,3,1,1);

        // Register action buttons

        plot.addActionListener(this);
        close.addActionListener(this);

    }

    public ZdenPlot(GUI here) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        home=here;
        println("Activated panel for plotting Z-density");
        job=new ZdenPlot();
        job.pack();
        job.setVisible(true);
        npnts=0;
        name="Name";
        fname="ZDNDAT";
        atom.setText(name);
        zdndat.setText(fname);
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
            name=atom.getText();
            fname=zdndat.getText();
            npnts=rdzden(fname,name);
            if(npnts>0) {
                zdenXY(npnts,name);
                if(graf != null)
                    graf.job.dispose();
                graf=new GraphDraw(home);
                graf.xlabel.setText("Z-Distance (A)");
                graf.ylabel.setText("D(z)");
                graf.plabel.setText("Z-Density of "+name.trim());
                graf.extraPlot(npnts,xx,yy);
            }
        }
        else if (arg.equals("Close")) {
            job.dispose();
        }
    }
    int rdzden(String fname,String atnam1) {
        /*
*********************************************************************

dl_poly/java routine to read a DL_POLY ZDNDAT file

copyright - daresbury laboratory
author    - w.smith march 2001

*********************************************************************
         */
        int nzdens,npts;
        boolean found=false;
        String record,anam1="";

        try {
            LineNumberReader lnr = new LineNumberReader(new FileReader(fname));
            println("Reading file: "+fname);
            title = lnr.readLine();
            println("File header record: "+title);
            record = lnr.readLine();
            nzdens=BML.giveInteger(record,1);
            npts=BML.giveInteger(record,2);
            xx=new double[npts];
            yy=new double[npts];
            OUT:
                for(int n=0;n<nzdens;n++) {
                    record=lnr.readLine();
                    anam1=BML.giveWord(record,1);

                    for(int i=0;i<npts;i++) {
                        record=lnr.readLine();
                        xx[i]=BML.giveDouble(record,1);
                        yy[i]=BML.giveDouble(record,2);
                    }

                    if(atnam1.equals(anam1)) {
                        found=true;
                        break OUT;
                    }
                }
                if(!found) {
                    println("Error - required Z-Density not found in file "+fname);
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
        return npts;

    }
    void zdenXY(int npts,String aname) {
        /*
*********************************************************************

dl_poly/java GUI routine to create a Z-Density XY file

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        String fname;
        String[] header;
        int nhead=4,call;

        header=new String[nhead];
        fname="Zden"+String.valueOf(numzdn)+".XY";
        numzdn++;
        header[0]=" Z-Density Plotting Program";
        header[1]=" Z-Density Plot: "+aname.trim();
        header[2]=" Z-Distance (A)";
        header[3]=" D(z)";
        call=putXY(fname,header,nhead,npts,xx,yy);
        if(call==0)
            println("PLOT file "+fname+" created");
    }
}

