import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;

public class EopVar extends Basic implements ActionListener {
        /*
*********************************************************************

dl_poly/java  Extra Options class

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
    public static MakeControl home;
    public static EopVar job;
    public static JTextField tbinsize,tdensvar,tdump,tewldev,tmxquat,tmxshak;
    public static JTextField tnfold1,tnfold2,tnfold3,tspmeev,tregauss;
    public static JButton close;
    public static JCheckBox blexclude,blmetdir,blvdwdir,blnoindex,blnostrict,blreplay;
    public static JCheckBox blpslab,blvdwshift,blnotopo;
    public static boolean lexclude,lmetdir,lvdwdir,lnoindex,lnostrict,lreplay,lpslab;
    public static boolean lvdwshift,lnotopo;
    public static int ndump,ewldev,spmeev,mxquat,mxshak,nfold1,nfold2,nfold3,nregauss;
    public static double binsize,densvar;

    // Define the Graphical User Interface

    public EopVar() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        super();
        setTitle("Extra Controls");
        int n=0;

        getContentPane().setForeground(art.fore);
        getContentPane().setBackground(art.back);
        setDefaultCloseOperation(DO_NOTHING_ON_CLOSE);
        setFont(fontMain);
        GridBagLayout grd = new GridBagLayout();
        GridBagConstraints gbc = new GridBagConstraints();
        getContentPane().setLayout(grd);

        gbc.fill=GridBagConstraints.BOTH;

        // Panel label

        fix(new JLabel("Select options:",JLabel.LEFT),grd,gbc,0,n++,3,1);

        // Binsize

        tbinsize = new JTextField(8);
        tbinsize.setForeground(art.scrf);
        tbinsize.setBackground(art.scrn);
        fix(tbinsize,grd,gbc,0,n,1,1);
        fix(new JLabel("RDF/ZDEN binsize",JLabel.LEFT),grd,gbc,1,n++,2,1);

        // Densvar

        tdensvar = new JTextField(8);
        tdensvar.setForeground(art.scrf);
        tdensvar.setBackground(art.scrn);
        fix(tdensvar,grd,gbc,0,n,1,1);
        fix(new JLabel("DensVar",JLabel.LEFT),grd,gbc,1,n++,2,1);

        // Dump Interval

        tdump = new JTextField(8);
        tdump.setBackground(art.scrn);
        tdump.setForeground(art.scrf);
        fix(tdump,grd,gbc,0,n,1,1);
        fix(new JLabel("Dump Interval",JLabel.LEFT),grd,gbc,1,n++,2,1);

        // Ewald Evaluation Interval

        tewldev = new JTextField(8);
        tewldev.setForeground(art.scrf);
        tewldev.setBackground(art.scrn);
        fix(tewldev,grd,gbc,0,n,1,1);
        fix(new JLabel("Ewald Interval",JLabel.LEFT),grd,gbc,1,n++,2,1);

        // SPME Evaluation Interval

        tspmeev = new JTextField(8);
        tspmeev.setForeground(art.scrf);
        tspmeev.setBackground(art.scrn);
        fix(tspmeev,grd,gbc,0,n,1,1);
        fix(new JLabel("SPME Interval",JLabel.LEFT),grd,gbc,1,n++,2,1);

        // Extended Coulomb Exclusions

        blexclude = new JCheckBox("Extend Coulomb Exclusions");
        blexclude.setForeground(art.fore);
        blexclude.setBackground(art.back);
        fix(blexclude,grd,gbc,0,n++,2,1);

        // Maximum Quaternion Iterations

        tmxquat = new JTextField(8);
        tmxquat.setForeground(art.scrf);
        tmxquat.setBackground(art.scrn);
        fix(tmxquat,grd,gbc,0,n,1,1);
        fix(new JLabel("Max Quatn. Iterations",JLabel.LEFT),grd,gbc,1,n++,2,1);

        // Maximum SHAKE Iterations

        tmxshak = new JTextField(8);
        tmxshak.setForeground(art.scrf);
        tmxshak.setBackground(art.scrn);
        fix(tmxshak,grd,gbc,0,n,1,1);
        fix(new JLabel("Max SHAKE Iterations",JLabel.LEFT),grd,gbc,1,n++,2,1);

        // NFold Expand CONFIG and FIELD

        fix(new JLabel("Nfold Expand",JLabel.LEFT),grd,gbc,0,n++,2,1);
        tnfold1 = new JTextField(6);
        tnfold1.setBackground(art.scrn);
        tnfold1.setForeground(art.scrf);
        fix(tnfold1,grd,gbc,0,n,1,1);
        tnfold2 = new JTextField(6);
        tnfold2.setBackground(art.scrn);
        tnfold2.setForeground(art.scrf);
        fix(tnfold2,grd,gbc,1,n,1,1);
        tnfold3 = new JTextField(6);
        tnfold3.setBackground(art.scrn);
        tnfold3.setForeground(art.scrf);
        fix(tnfold3,grd,gbc,2,n++,1,1);

        // Direct Metal Evaluation

        blmetdir = new JCheckBox("Metal Direct");
        blmetdir.setForeground(art.fore);
        blmetdir.setBackground(art.back);
        fix(blmetdir,grd,gbc,0,n++,1,1);

        // Direct VDW Evaluation

        blvdwdir = new JCheckBox("VDW Direct");
        blvdwdir.setForeground(art.fore);
        blvdwdir.setBackground(art.back);
        fix(blvdwdir,grd,gbc,0,n++,2,1);

        // Ignore CONFIG file Indices

        blnoindex = new JCheckBox("No Index");
        blnoindex.setForeground(art.fore);
        blnoindex.setBackground(art.back);
        fix(blnoindex,grd,gbc,0,n++,2,1);

        // Abort strict data checks

        blnostrict = new JCheckBox("No Strict");
        blnostrict.setForeground(art.fore);
        blnostrict.setBackground(art.back);
        fix(blnostrict,grd,gbc,0,n++,2,1);

        // No topology print option

        blnotopo = new JCheckBox("No Topology");
        blnotopo.setForeground(art.fore);
        blnotopo.setBackground(art.back);
        fix(blnotopo,grd,gbc,0,n++,2,1);

        // Regauss interval

        tregauss = new JTextField(8);
        tregauss.setBackground(art.scrn);
        tregauss.setForeground(art.scrf);
        fix(tregauss,grd,gbc,0,n,1,1);
        fix(new JLabel("Regauss Interval",JLabel.LEFT),grd,gbc,1,n++,2,1);

        // Enable replay option

        blreplay = new JCheckBox("Enable Replay");
        blreplay.setForeground(art.fore);
        blreplay.setBackground(art.back);
        fix(blreplay,grd,gbc,0,n++,2,1);

        // Processor slab option

        blpslab = new JCheckBox("Processor Slab");
        blpslab.setForeground(art.fore);
        blpslab.setBackground(art.back);
        fix(blpslab,grd,gbc,0,n++,2,1);

        // Shift VDW forces option

        blvdwshift = new JCheckBox("VDW Shift");
        blvdwshift.setForeground(art.fore);
        blvdwshift.setBackground(art.back);
        fix(blvdwshift,grd,gbc,0,n++,2,1);

        // Define the Close button

        close = new JButton("Close");
        close.setBackground(art.butn);
        close.setForeground(art.butf);
        fix(close,grd,gbc,0,n++,1,1);

        // Register action buttons

        close.addActionListener(this);

    }

    public EopVar(MakeControl here) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        home=here;
        job=new EopVar();
        job.pack();
        job.setVisible(true);
        setParams();
    }

    void setParams() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */

        // set parameters values

        home.seteopValues();

        // set panel

        tbinsize.setText(String.valueOf(binsize));
        tdensvar.setText(String.valueOf(densvar));
        tdump.setText(String.valueOf(ndump));
        tewldev.setText(String.valueOf(ewldev));
        tspmeev.setText(String.valueOf(spmeev));
        blexclude.setSelected(lexclude);
        tmxquat.setText(String.valueOf(mxquat));
        tmxshak.setText(String.valueOf(mxshak));
        tnfold1.setText(String.valueOf(nfold1));
        tnfold2.setText(String.valueOf(nfold2));
        tnfold3.setText(String.valueOf(nfold3));
        blmetdir.setSelected(lmetdir);
        blvdwdir.setSelected(lvdwdir);
        blnoindex.setSelected(lnoindex);
        blnostrict.setSelected(lnostrict);
        blnotopo.setSelected(lnotopo);
        tregauss.setText(String.valueOf(nregauss));
        blreplay.setSelected(lreplay);
        blpslab.setSelected(lpslab);
        blvdwshift.setSelected(lvdwshift);

    }

    public void actionPerformed(ActionEvent e) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */

        String arg = (String)e.getActionCommand();
        if (arg.equals("Close")) {
            byebye();
        }
    }

    void byebye() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        home.geteopValues();
        home.eop=null;
        job.setVisible(false);
    }
}
