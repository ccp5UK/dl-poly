import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;

public class HovePlot extends Basic implements ActionListener {
        /*
*********************************************************************

dl_poly/java GUI class to plot van Hove correlation functions

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
    public static HovePlot job;
    private static GUI home;
    private static String name1,name2,title,hname;
    private static double[] xx,yy;
    private static int npnts,nseq;
    private static JTextField seqno,hovefile;
    private static JCheckBox dens;
    private static JButton plot,close;
    private static boolean ldens;
    private static double time;

    // Define the Graphical User Interface

    public HovePlot() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        setTitle("van Hove Plotter");

        getContentPane().setBackground(art.back);
        getContentPane().setForeground(art.fore);
        setFont(fontMain);
        GridBagLayout grd = new GridBagLayout();
        GridBagConstraints gbc = new GridBagConstraints();
        getContentPane().setLayout(grd);

        gbc.fill=GridBagConstraints.BOTH;

        // Define the Plot button

        plot = new JButton("  Plot  ");
        plot.setBackground(art.butn);
        plot.setForeground(art.butf);
        fix(plot,grd,gbc,0,0,1,1);
        fix(new JLabel("  "),grd,gbc,1,0,1,1);

        // Define the Close button

        close = new JButton("Close");
        close.setBackground(art.butn);
        close.setForeground(art.butf);
        fix(close,grd,gbc,2,0,1,1);

        // Name of van Hove data file

        JLabel lab0 = new JLabel("van Hove plot file:",JLabel.LEFT);
        fix(lab0,grd,gbc,0,1,2,1);
        hovefile = new JTextField(14);
        hovefile.setBackground(art.scrn);
        hovefile.setForeground(art.scrf);
        fix(hovefile,grd,gbc,0,2,3,1);

        // Function sequence number

        JLabel lab1 = new JLabel("Plot sequence No.:",JLabel.LEFT);
        fix(lab1,grd,gbc,0,3,2,1);
        seqno = new JTextField(5);
        seqno.setBackground(art.scrn);
        seqno.setForeground(art.scrf);
        fix(seqno,grd,gbc,2,3,1,1);

        // Radial density option

        JLabel lab2 = new JLabel("Radial density opt.:",JLabel.LEFT);
        fix(lab2,grd,gbc,0,4,2,1);
        dens = new JCheckBox("    ");
        dens.setBackground(art.back);
        dens.setForeground(art.fore);
        fix(dens,grd,gbc,2,4,1,1);

        // Register action buttons

        plot.addActionListener(this);
        close.addActionListener(this);

    }

    public HovePlot(GUI here){
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        home=here;
        println("Activated panel for plotting van Hove functions");
        job=new HovePlot();
        job.pack();
        job.setVisible(true);
        npnts=0;
        nseq=0;
        ldens=false;
	hovefile.setText(null);
        dens.setSelected(ldens);
        seqno.setText(String.valueOf(nseq));
    }

    public HovePlot(GUI here,String hovename) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        home=here;
        println("Activated panel for plotting van Hove functions");
        job=new HovePlot();
        job.pack();
        job.setVisible(true);
        npnts=0;
        nseq=0;
        ldens=false;
	hovefile.setText(hovename);
        dens.setSelected(ldens);
        seqno.setText(String.valueOf(nseq));
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
        if (arg.equals("  Plot  ")) {
	    hname=hovefile.getText();
	    if(hname.length()==0)hname=null;
            nseq=BML.giveInteger(seqno.getText(),1);
            ldens=dens.isSelected();
	    if(hname == null){
		hname=selectFileNameBegins(home,"HOVG");
		println("File selected is "+hname+hname.toUpperCase().indexOf("HOVGSL"));
	    }
            if(hname != null) {
                npnts=rdhove(ldens,nseq);
                if(npnts>0) {
                    hovXY(npnts,name1,name2);
                    if(graf != null)
                        graf.job.setVisible(false);
                    graf=new GraphDraw(home);
                    graf.xlabel.setText("Radius (A)");
                    if(hname.toUpperCase().indexOf("HOVGSL")>=0) {
                        graf.ylabel.setText("Gs(r,t)");
                        graf.plabel.setText("Gs(r,t) of "+name1.trim()+" - "+name2.trim()+" at"+BML.fmt(time,8)+" ps");
                        graf.extraPlot(npnts,xx,yy);
                    }
                    else if(hname.toUpperCase().indexOf("HOVGDF")>=0) {
                        graf.ylabel.setText("Gd(r,t)");
                        graf.plabel.setText("Gd(r,t) of "+name1.trim()+" - "+name2.trim()+" at"+BML.fmt(time,8)+" ps");
                        graf.extraPlot(npnts,xx,yy);
                    }
                    else
                        println("Error nominated file not van Hove data file");
                }
            }
        }
        else if (arg.equals("Close")) {
            job.setVisible(false);
        }
    }
    int rdhove(boolean lden,int nseq) {
        /*
*********************************************************************

dl_poly/java routine to read a DL_POLY van Hove Data file

copyright - daresbury laboratory
author    - w.smith march 2001

*********************************************************************
         */
        int nhovs,npts;
        String record;

        try {
            LineNumberReader lnr = new LineNumberReader(new FileReader(hname));
            println("Reading file: "+hname);
            title = lnr.readLine();
            println("File header record: "+title);
            record = lnr.readLine();
            name1=BML.giveWord(record,1);
            name2=BML.giveWord(record,2);
            record = lnr.readLine();
            nhovs=BML.giveInteger(record,1);
            npts=BML.giveInteger(record,2);
            xx=new double[npts];
            yy=new double[npts];
            if(nseq>nhovs) {
                println("Error - requested function not present in file");
                lnr.close();
                return -1;
            }
            OUT:
                for(int n=0;n<nhovs;n++) {
                    record=lnr.readLine();
                    time=BML.giveDouble(record,2);
                    for(int i=0;i<npts;i++) {
                        record=lnr.readLine();
                        xx[i]=BML.giveDouble(record,1);
                        yy[i]=BML.giveDouble(record,2);
                    }

                    if(n==nseq) {
                        break OUT;
                    }
                }
                lnr.close();
        }
        catch(FileNotFoundException e) {
            println("Error - file not found: " + hname);
            return -2;
        }
        catch(Exception e) {
            println("Error reading file: " + hname + " "+e);
            return -3;
        }
        println("Number of points ploted:"+BML.fmt(npts,6));

        // activate radial density option

        if(lden) {
            for(int i=0;i<npts;i++)
                yy[i]*=(4.0*Math.PI*Math.pow(xx[i],2));
        }

        return npts;
    }
    void hovXY(int npts,String anam1,String anam2) {
        /*
*********************************************************************

dl_poly/java GUI routine to create a van Hove XY file

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        String gname;
        String[] header;
        int nhead=4,call;

        header=new String[nhead];
        gname="HOV"+String.valueOf(numhov)+".XY";
        numhov++;
        header[0]=" van Hove Plotting Program";
        header[2]=" Radius (A)";
        if(hname.toUpperCase().indexOf("HOVGSL")>=0) {
            header[1]=" Gs(r,t) Plot: "+anam1.trim()+" at"+BML.fmt(time,8)+" ps";
            header[3]=" Gs(r,t)";
        }
        else {
            header[1]=" Gd(r,t) Plot: "+anam1.trim()+"+"+anam2.trim()+" at"+BML.fmt(time,8)+" ps";
            header[3]=" Gd(r,t)";
        }
        call=putXY(gname,header,nhead,npts,xx,yy);
        if(call==0)
            println("PLOT file "+gname+" created");
    }
}

