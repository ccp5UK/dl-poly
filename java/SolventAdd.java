import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;

// Define the Graphical User Interface

public class SolventAdd extends Basic implements ActionListener {
    /*
*********************************************************************

dl_poly/java GUI class to add solvent to a CONFIG file

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
     */
    public static SolventAdd job;
    private static Config solute,solvent;
    private static String solventfile;
    private static double slx,sly,slz,top,bot,dws,dww;
    private static JButton load,make,close;
    private static JCheckBox slab;
    private static JTextField dlx,dly,dlz,ubd,lbd,wsd,wwd,sfil;
    private static boolean lslab,lfetch;
    private static int[] hitlist;
    private static GUI home=null;

    public SolventAdd() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        setTitle("Add Solvent Panel");

        getContentPane().setBackground(art.back);
        getContentPane().setForeground(art.fore);
        setFont(fontMain);
        GridBagLayout grd = new GridBagLayout();
        GridBagConstraints gbc = new GridBagConstraints();
        getContentPane().setLayout(grd);

        gbc.fill=GridBagConstraints.BOTH;

        // Define the Load button

        load = new JButton("Load");
        load.setBackground(art.butn);
        load.setForeground(art.butf);
        fix(load,grd,gbc,0,0,1,1);
        fix(new JLabel("        "),grd,gbc,1,0,1,1);

        // Define the Make button

        make = new JButton("Make");
        make.setBackground(art.butn);
        make.setForeground(art.butf);
        fix(make,grd,gbc,2,0,1,1);

        // Solvent-solute minimum distance

        JLabel lab1 = new JLabel("Min solvent-solute distance:",JLabel.LEFT);
        fix(lab1,grd,gbc,0,1,2,1);
        wsd = new JTextField(8);
        wsd.setBackground(art.scrn);
        wsd.setForeground(art.scrf);
        fix(wsd,grd,gbc,2,1,1,1);

        // Solvent-solvent minimum distance

        JLabel lab2 = new JLabel("Min solvent-solvent distance:",JLabel.LEFT);
        fix(lab2,grd,gbc,0,2,2,1);
        wwd = new JTextField(8);
        wwd.setBackground(art.scrn);
        wwd.setForeground(art.scrf);
        fix(wwd,grd,gbc,2,2,1,1);

        // Solvent file nomination

        JLabel lab3 = new JLabel("Solvent file:",JLabel.LEFT);
        fix(lab3,grd,gbc,0,3,2,1);
        sfil = new JTextField(8);
        sfil.setBackground(art.scrn);
        sfil.setForeground(art.scrf);
        fix(sfil,grd,gbc,2,3,1,1);

        // Solvent slab option

        slab=new JCheckBox("Slab?");
        slab.setBackground(art.back);
        slab.setForeground(art.fore);
        fix(slab,grd,gbc,0,4,1,1);

        // Slab direction vector

        JLabel lab4 = new JLabel("Slab direction vector:",JLabel.LEFT);
        fix(lab4,grd,gbc,0,5,3,1);
        dlx = new JTextField(8);
        dlx.setBackground(art.scrn);
        dlx.setForeground(art.scrf);
        fix(dlx,grd,gbc,0,6,1,1);

        dly = new JTextField(8);
        dly.setBackground(art.scrn);
        dly.setForeground(art.scrf);
        fix(dly,grd,gbc,1,6,1,1);

        dlz = new JTextField(8);
        dlz.setBackground(art.scrn);
        dlz.setForeground(art.scrf);
        fix(dlz,grd,gbc,2,6,1,1);

        // Upper bound of slab

        JLabel lab5 = new JLabel("Upper bound:",JLabel.RIGHT);
        fix(lab5,grd,gbc,0,7,2,1);
        ubd = new JTextField(8);
        ubd.setBackground(art.scrn);
        ubd.setForeground(art.scrf);
        fix(ubd,grd,gbc,2,7,1,1);

        // Lower bound of slab

        JLabel lab6 = new JLabel("Lower bound:",JLabel.RIGHT);
        fix(lab6,grd,gbc,0,8,2,1);
        lbd = new JTextField(8);
        lbd.setBackground(art.scrn);
        lbd.setForeground(art.scrf);
        fix(lbd,grd,gbc,2,8,1,1);

        // Define the Close button

        close = new JButton("Close");
        close.setBackground(art.butn);
        close.setForeground(art.butf);
        fix(close,grd,gbc,0,9,1,1);

        // Register action buttons

        load.addActionListener(this);
        make.addActionListener(this);
        close.addActionListener(this);

    }

    // Constructor method

    public SolventAdd(GUI here) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        home=here;
        println("Activated panel for adding solvent to a CONFIG file");

        // Set up Graphical User interface

        job = new SolventAdd();
        job.pack();
        job.setVisible(true);
        setValues();
    }

    // Set default values

    static void setValues() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        dww=2.5;
        dws=2.5;
        slx=0.0;
        sly=0.0;
        slz=1.0;
        top=3.0;
        bot=-3.0;
        lslab=false;
        lfetch=false;
        solventfile="WATER300K";
        sfil.setText(solventfile);
        slab.setSelected(lslab);
        wsd.setText(String.valueOf(dws));
        wwd.setText(String.valueOf(dww));
        dlx.setText(String.valueOf(slx));
        dly.setText(String.valueOf(sly));
        dlz.setText(String.valueOf(slz));
        ubd.setText(String.valueOf(top));
        lbd.setText(String.valueOf(bot));
    }

    public void actionPerformed(ActionEvent e) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        int call;
        String arg = (String)e.getActionCommand();
        if (arg.equals("Load")) {
            lfetch=true;
            lslab=slab.isSelected();
            dws=BML.giveDouble(wsd.getText(),1);
            dww=BML.giveDouble(wwd.getText(),1);
            slx=BML.giveDouble(dlx.getText(),1);
            sly=BML.giveDouble(dly.getText(),1);
            slz=BML.giveDouble(dlz.getText(),1);
            top=BML.giveDouble(ubd.getText(),1);
            bot=BML.giveDouble(lbd.getText(),1);
            if(!addsolvent())
                println("Problem encountered adding solvent");
        }
        else if (arg.equals("Make")) {
            lfetch=false;
            if(config == null)lfetch=true;
            lslab=slab.isSelected();
            dws=BML.giveDouble(wsd.getText(),1);
            dww=BML.giveDouble(wwd.getText(),1);
            slx=BML.giveDouble(dlx.getText(),1);
            sly=BML.giveDouble(dly.getText(),1);
            slz=BML.giveDouble(dlz.getText(),1);
            top=BML.giveDouble(ubd.getText(),1);
            bot=BML.giveDouble(lbd.getText(),1);
            if(!addsolvent())
                println("Problem encountered adding solvent");
        }
        else if (arg.equals("Close")) {
            job.setVisible(false);
        }
    }

    boolean addsolvent() {
        /*
***********************************************************************

dl_poly/java utility to add SPC solvent molecules to a structure to fill
out the MD cell.
Assumes the periodic boundary is cubic, orthorhombic, trunctated octahedral
or rhombic dodecahedral

copyright - daresbury laboratory
author    - w. smith 2011

***********************************************************************
         */

        // Load solute if required

        if(lfetch || config == null)
            solute=getConfig(home,"CFG");
        else
            solute=copyConfig(config);

        // Check periodic boundary

        if(solute.pbc.imcon == 0 || solute.pbc.imcon == 3 ||
           solute.pbc.imcon == 6 || solute.pbc.imcon == 7) {
            println("Error - unsuitable CONFIG PBC for solvent addition");
            return false;
        }

        // load and check solvent data

        if(!getSolvent()) {
            println("Error - problem processing solvent data");
            return false;
        }

        // Load solvent molecules into CONFIG file

        if(fillCell())
            println("Solvent added to cell");
        else {
            println("Problem encountered adding solvent to cell");
            return false;
        }

        // Remove overlapping solvent molecules

        if(killOverlaps())
            println("Overlapping solvent molecules removed");
        else {
            println("Problem encountered removing overlapping solvent molecules");
            return false;
        }

        // Cut solvent slab if required

        if(lslab)cutSlab();

        // Put solute back into solvent

        if(insertSolute())
            println("Sovent addition completed");
        else {
            println("Problem encountered re-inserting solute");
            return false;
        }

        // write new CONFIG file

        fname="CFGSOL."+String.valueOf(numsol);
        if(config.configWrite(fname)){
            println("File "+fname+" created");
            println("Number of atoms in "+fname+" : "+config.natms);
            numsol++;
        }

        // Draw modified structure

	config.pbc.images(config.natms,config.xyz);
        config.structure=new Structure(config);
        if(home.editor != null)
            home.editor.job.setVisible(false);
        home.editor=new Editor(home);

        return true;
    }

    boolean getSolvent() {
        /*
***********************************************************************

dl_poly/java utility to process a solvent file for SolventAdd program

copyright - daresbury laboratory
author    - w. smith 2011

***********************************************************************
*/

        int j;
        double uvw[]=new double[3];

        // Load solvent

        solvent=new Config();

        if(!solvent.rdCFG(solventfile)) {
            println("Error - reading solvent file:" + solventfile);
            return false;
        }

        // Check periodic boundary of solvent file

        if(solvent.pbc.imcon != 1 && solvent.pbc.imcon != 2) {
            println("Error - solvent file PBC must be 1 or 2");
            return false;
        }

        // Make sure solvent molecules not split by PBC

        for(int i=0;i<solvent.natms;i++) {
            for(int m=0;m<solvent.structure.lbnd[i];m++) {
                j=solvent.structure.bond[m][i];
                if(j > i) {
                    uvw[0]=solvent.xyz[0][j]-solvent.xyz[0][i];
                    uvw[1]=solvent.xyz[1][j]-solvent.xyz[1][i];
                    uvw[2]=solvent.xyz[2][j]-solvent.xyz[2][i];
                    solvent.pbc.images(uvw);
                    solvent.xyz[0][j]=solvent.xyz[0][i]+uvw[0];
                    solvent.xyz[1][j]=solvent.xyz[1][i]+uvw[1];
                    solvent.xyz[2][j]=solvent.xyz[2][i]+uvw[2];
                }
            }
        }

        return true;
    }

    boolean fillCell() {
        /*
***********************************************************************

dl_poly/java utility to fill a CONFIG simulation cell with solvent
molecules

copyright - daresbury laboratory
author    - w. smith 2011

***********************************************************************
         */
        boolean keep;
        int n,ncx,ncy,ncz,mol;
        double xxx,yyy,zzz,basx,basy,basz,cellx,celly,cellz,octa,rhom;
        double RT2;

        // Create new CONFIG file

        config=new Config();
        config.pbc.imcon=solute.pbc.imcon;
        for(int i=0;i<9;i++)
            config.pbc.cell[i]=solute.pbc.cell[i];
        config.pbc.buildBoundary(config.pbc.imcon);

        // insert solvent molecules into cell

        n=0;
        ncx=(int)(config.pbc.cell[0]/solvent.pbc.cell[0])+1;
        ncy=(int)(config.pbc.cell[4]/solvent.pbc.cell[4])+1;
        ncz=(int)(config.pbc.cell[8]/solvent.pbc.cell[8])+1;
        basz=-0.5*solvent.pbc.cell[8]*(double)(ncz-1);
        for(int k=0;k<ncz;k++) {
            basy=-0.5*solvent.pbc.cell[4]*(double)(ncy-1);
            for(int j=0;j<ncy;j++) {
                basx=-0.5*solvent.pbc.cell[0]*(double)(ncx-1);
                for(int i=0;i<ncx;i++) {
                    for(int m=0;m<solvent.natms;m++) {
                        if(n == config.atoms.length) {
                            config.natms=n;
                            config.resizeArrays();
                        }
                        config.atoms[n]=new Element();
                        config.atoms[n].znum=solvent.atoms[m].znum;
                        config.atoms[n].zmas=solvent.atoms[m].zmas;
                        config.atoms[n].zchg=solvent.atoms[m].zchg;
                        config.atoms[n].zrad=solvent.atoms[m].zrad;
                        config.atoms[n].zsym=new String(solvent.atoms[m].zsym);
                        config.atoms[n].zcol=new Color(solvent.atoms[m].zcol.getRGB());
                        config.atoms[n].covalent=solvent.atoms[m].covalent;
                        config.xyz[0][n]=basx+solvent.xyz[0][m];
                        config.xyz[1][n]=basy+solvent.xyz[1][m];
                        config.xyz[2][n]=basz+solvent.xyz[2][m];
                        n++;
                    }
                    basx+=solvent.pbc.cell[0];
                }
                basy+=solvent.pbc.cell[4];
            }
            basz+=solvent.pbc.cell[8];
        }
        config.natms=n;

        // flag atoms outside the PBC

        hitlist=new int[config.natms];

        RT2=Math.sqrt(2.0);
        cellx=config.pbc.cell[0]/2.0;
        celly=config.pbc.cell[4]/2.0;
        cellz=config.pbc.cell[8]/2.0;
        octa=1.5*cellx+0.5*dww;
        rhom=2.0*cellx+0.5*dww;
        for(int i=0;i<config.natms;i++) {
            hitlist[i]=1;
            if(config.atoms[i].znum > 1) {
                xxx=Math.abs(config.xyz[0][i]);
                yyy=Math.abs(config.xyz[1][i]);
                zzz=Math.abs(config.xyz[2][i]);
                if(xxx > cellx+0.5*dww)
                    hitlist[i]=-1;
                else if(yyy > celly+0.5*dww)
                    hitlist[i]=-1;
                else if(zzz > cellz+0.5*dww)
                    hitlist[i]=-1;
                else if(config.pbc.imcon == 4) {
                    if(xxx+yyy+zzz > octa)
                        hitlist[i]=-1;
                }
                else if(config.pbc.imcon == 5) {
                    if(xxx+yyy+RT2*zzz >= rhom)
                        hitlist[i]=-1;
                }
            }
        }

        // flag molecules for removal

        n=0;
        mol=1;
        for(int k=0;k<ncx*ncy*ncz;k++) {
            for(int m=0;m<solvent.structure.nmols;m++) {
                keep=true;
                for(int j=0;j<solvent.structure.isz[m];j++) {
                    if(hitlist[n] < 0) keep=false;
                    n++;
                }
                n-=solvent.structure.isz[m];
                for(int j=0;j<solvent.structure.isz[m];j++) {
                    if(keep)
                        hitlist[n]=mol;
                    else
                        hitlist[n]=-mol;
                    n++;
                }
                if(keep) mol++;
            }
        }

        // remove molecules outside the PBC

        n=0;
        for(int i=0;i<config.natms;i++) {
            if(hitlist[i] > 0) {
                config.atoms[n].znum=config.atoms[i].znum;
                config.atoms[n].zmas=config.atoms[i].zmas;
                config.atoms[n].zchg=config.atoms[i].zchg;
                config.atoms[n].zrad=config.atoms[i].zrad;
                config.atoms[n].zsym=new String(config.atoms[i].zsym);
                config.atoms[n].zcol=new Color(config.atoms[i].zcol.getRGB());
                config.atoms[n].covalent=config.atoms[i].covalent;
                config.xyz[0][n]=config.xyz[0][i];
                config.xyz[1][n]=config.xyz[1][i];
                config.xyz[2][n]=config.xyz[2][i];
                hitlist[n]=hitlist[i];
                n++;
            }
        }
        config.natms=n;

        return true;
    }

    boolean killOverlaps() {
        /*
***********************************************************************

dl_poly/java utility to remove molecules that are too close across a
periodic boundary

copyright - daresbury laboratory
author    - w. smith 2011

***********************************************************************
         */

        boolean zapp;
        int n,k0,mol;
        double uvw[]=new double[3];
        double dsq,rsq;

        rsq=0.0;
        dsq=dww*dww;

        // identify atoms close to cell boundary

        k0=0;
        zapp=false;
        mol=hitlist[0];

        for(int i=0;i<config.natms-1;i++) {

            if(config.atoms[i].znum > 1) {
                if(!zapp) {
                    for(int j=i+1;j<config.natms;j++) {
                        if(hitlist[j] > mol && config.atoms[j].znum > 1) {
                            uvw[0]=config.xyz[0][i]-config.xyz[0][j];
                            uvw[1]=config.xyz[1][i]-config.xyz[1][j];
                            uvw[2]=config.xyz[2][i]-config.xyz[2][j];
                            config.pbc.images(uvw);
                            rsq=Math.pow(uvw[0],2)+Math.pow(uvw[1],2)+Math.pow(uvw[2],2);
                            if(dsq > rsq) {
                                zapp=true;
                                break;
                            }
                        }
                    }
                }
            }
            if(hitlist[i+1] > mol) {
                if(zapp) {
                    for(int k=k0;k<=i;k++)
                        hitlist[k]=-mol;
                }
                mol=hitlist[i+1];
                zapp=false;
                k0=i+1;
            }

        }

        // remove marked molecules

        n=0;
        for(int i=0;i<config.natms;i++) {
            if(hitlist[i] > 0) {
                config.atoms[n].znum=config.atoms[i].znum;
                config.atoms[n].zmas=config.atoms[i].zmas;
                config.atoms[n].zchg=config.atoms[i].zchg;
                config.atoms[n].zrad=config.atoms[i].zrad;
                config.atoms[n].zsym=new String(config.atoms[i].zsym);
                config.atoms[n].zcol=new Color(config.atoms[i].zcol.getRGB());
                config.atoms[n].covalent=config.atoms[i].covalent;
                config.xyz[0][n]=config.xyz[0][i];
                config.xyz[1][n]=config.xyz[1][i];
                config.xyz[2][n]=config.xyz[2][i];
                hitlist[n]=hitlist[i];
                n++;
            }
        }
        config.natms=n;

        return true;
    }

    void cutSlab() {
        /*
***********************************************************************

dl_poly/java utility to cut a slab out of a solvent configuration

copyright - daresbury laboratory
author    - w.smith 2011

***********************************************************************
         */

        int k,j0,mol;
        boolean zapp;
        double ddd,sss;

	config.pbc.images(config.natms,config.xyz);

        // identify solvent atoms outside of the slab

	j0=0;
	zapp=false;
	mol=hitlist[0];
	sss=Math.sqrt(slx*slx+sly*sly+slz*slz);
	for(int i=0;i<config.natms;i++) {
	    if(config.atoms[i].znum > 1) {
                ddd=(slx*config.xyz[0][i]+sly*config.xyz[1][i]+slz*config.xyz[2][i])/sss;
                if(ddd < bot || ddd > top) zapp=true;
	    }
	    if(i+1 == config.natms || hitlist[i+1] > mol) {
		if(zapp) {
		    for(int j=j0;j<=i;j++)
			hitlist[j]=-mol;
		}
		if(i+1 < config.natms)mol=hitlist[i+1];
		zapp=false;
		j0=i+1;
	    }
	}

        // Remove redundant molecules

        k=0;
        for(int i=0;i<config.natms;i++) {
            if(hitlist[i] > 0) {
                config.atoms[k].znum=config.atoms[i].znum;
                config.atoms[k].zmas=config.atoms[i].zmas;
                config.atoms[k].zchg=config.atoms[i].zchg;
                config.atoms[k].zrad=config.atoms[i].zrad;
                config.atoms[k].zsym=new String(config.atoms[i].zsym);
                config.atoms[k].zcol=new Color(config.atoms[i].zcol.getRGB());
                config.atoms[k].covalent=config.atoms[i].covalent;
                config.xyz[0][k]=config.xyz[0][i];
                config.xyz[1][k]=config.xyz[1][i];
                config.xyz[2][k]=config.xyz[2][i];
		hitlist[k]=hitlist[i];
                k++;
            }
        }
        config.natms=k;

    }

    boolean insertSolute() {
        /*
***********************************************************************

dl_poly/java utility to add a solute to a solvent with removal of
solven molecules to accommodate the insert.
Assumes atomic positions are in a form compatible
with the CONFIG file used in DL_POLY with periodic boundary
conditions.

copyright - daresbury laboratory
author    - w.smith 2011

***********************************************************************
         */

        int k,j0,mol;
        boolean zapp;
        double ddd,dis2;
        double dsp[]=new double[3];

        dis2=dws*dws;

        // identify solvent atoms overlapping with solute atoms

        for(int j=0;j<solute.natms;j++) {
            if(solute.atoms[j].znum > 1) {
                for(int i=0;i<config.natms;i++) {
                    if(config.atoms[i].znum > 1) {
                        dsp[0]=config.xyz[0][i]-solute.xyz[0][j];
                        dsp[1]=config.xyz[1][i]-solute.xyz[1][j];
                        dsp[2]=config.xyz[2][i]-solute.xyz[2][j];
                        config.pbc.images(dsp);
                        ddd=dsp[0]*dsp[0]+dsp[1]*dsp[1]+dsp[2]*dsp[2];
                        if(ddd < dis2) hitlist[i]=-Math.abs(hitlist[i]);
                    }
                }
            }
        }

        // flag overlapping molecules for removal

	j0=0;
	zapp=false;
	mol=Math.abs(hitlist[0]);
	for(int i=0;i<config.natms;i++) {
	    if(config.atoms[i].znum > 1) {
		if(hitlist[i] < 0) zapp=true;
	    }
	    if(i+1 == config.natms || Math.abs(hitlist[i+1]) > mol) {
		if(zapp) {
		    for(int j=j0;j<=i;j++)
			hitlist[j]=-mol;
		}
		if(i+1 < config.natms)
		    mol=Math.abs(hitlist[i+1]);
		zapp=false;
		j0=i+1;
	    }
	}

	// Remove redundant molecules

        k=0;
        for(int i=0;i<config.natms;i++) {
            if(hitlist[i] > 0) {
                config.atoms[k].znum=config.atoms[i].znum;
                config.atoms[k].zmas=config.atoms[i].zmas;
                config.atoms[k].zchg=config.atoms[i].zchg;
                config.atoms[k].zrad=config.atoms[i].zrad;
                config.atoms[k].zsym=new String(config.atoms[i].zsym);
                config.atoms[k].zcol=new Color(config.atoms[i].zcol.getRGB());
                config.atoms[k].covalent=config.atoms[i].covalent;
                config.xyz[0][k]=config.xyz[0][i];
                config.xyz[1][k]=config.xyz[1][i];
                config.xyz[2][k]=config.xyz[2][i];
                k++;
            }
        }

        // Add solute to config

        for(int i=0;i<solute.natms;i++) {
            if(k == config.atoms.length) {
                config.natms=k;
                config.resizeArrays();
            }
            config.atoms[k]=new Element();
            config.atoms[k].znum=solute.atoms[i].znum;
            config.atoms[k].zmas=solute.atoms[i].zmas;
            config.atoms[k].zchg=solute.atoms[i].zchg;
            config.atoms[k].zrad=solute.atoms[i].zrad;
            config.atoms[k].zsym=new String(solute.atoms[i].zsym);
            config.atoms[k].zcol=new Color(solute.atoms[i].zcol.getRGB());
            config.atoms[k].covalent=solute.atoms[i].covalent;
            config.xyz[0][k]=solute.xyz[0][i];
            config.xyz[1][k]=solute.xyz[1][i];
            config.xyz[2][k]=solute.xyz[2][i];
            k++;
        }
        config.natms=k;

        return true;
    }

}