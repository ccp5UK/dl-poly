import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;

public class MakeControl extends Basic implements ActionListener {
        /*
*********************************************************************

dl_poly/java GUI class to make CONTROL files

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
    public static GUI home=null;
    public static MakeControl job=null;
    public static SysVar sys=null;
    public static ProVar prv=null;
    public static ComVar cmv=null;
    public static EopVar eop=null;
    public static RadDam rdm=null;
    private static JTextField ttitle;
    private static String title,fname;
    private static JButton make,edit,close,sysbut,prvbut,cmvbut,rdmbut,eopbut;
    private static int numctr,nstrun,nsteql,mult,nstbpo,nstack,intsta,keyens,keyres,istrdf;
    private static int nstraj,istraj,levcon,nstbts,keyfce,keyalg,ndump,ewldev,spmeev,mxquat;
    private static int mxshak,nfold1,nfold2,nfold3,nregauss,istzden;
    private static int dstart,impstp,dintval,atom,nstmsdtmp,imsdtmp,ipseudtyp;
    private static boolean allpairs,lcap,lzeql,lrdf,lprdf,ltraj,ltscal,lzden,lzero,lvdw;
    private static boolean lpzden,lexclude,ldefects,lvarstp,lmsdtmp,lpseudo,lnotopo;
    private static boolean lmetdir,lvdwdir,lnoindex,lnostrict,lreplay,lpslab,lvdwshift;
    private static double temp,press,tstep,rcut,delr,rvdw,rprim,epsq,taut,taup,ewltol;
    private static double jobtim,tclose,fcap,gamt,binsize,densvar,shktol,qtntol;
    public static double defcut,energy,vect1,vect2,vect3,varstp;
    public static double mindis,maxdis,thick,ptemp;

    // Define the Graphical User Interface

    public MakeControl() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        super();
        setTitle("Make CONTROL File");
        int n=0;

        getContentPane().setBackground(art.back);
        getContentPane().setForeground(art.fore);
        setDefaultCloseOperation(DO_NOTHING_ON_CLOSE);
        setFont(fontMain);
        GridBagLayout grd = new GridBagLayout();
        GridBagConstraints gbc = new GridBagConstraints();
        getContentPane().setLayout(grd);

        gbc.fill=GridBagConstraints.BOTH;

        // Define the Make button

        make = new JButton("Make");
        make.setBackground(art.butn);
        make.setForeground(art.butf);
        fix(make,grd,gbc,0,n,1,1);

        // Define the Edit button

        edit = new JButton("Edit");
        edit.setBackground(art.butn);
        edit.setForeground(art.butf);
        fix(edit,grd,gbc,2,n++,1,1);

        // Blank line

        fix(new JLabel(" "),grd,gbc,0,n++,1,1);

        // File header

        fix(new JLabel("File Header",JLabel.LEFT),grd,gbc,0,n++,1,1);
        ttitle = new JTextField(22);
        ttitle.setBackground(art.scrn);
        ttitle.setForeground(art.scrf);
        fix(ttitle,grd,gbc,0,n++,3,1);

        // Blank line

        fix(new JLabel(" "),grd,gbc,0,n++,1,1);

        // Define the system variables button

        sysbut = new JButton("System Controls");
        sysbut.setBackground(art.butn);
        sysbut.setForeground(art.butf);
        fix(sysbut,grd,gbc,1,n++,1,1);

        // Define the program variables button

        prvbut = new JButton("Program Controls");
        prvbut.setBackground(art.butn);
        prvbut.setForeground(art.butf);
        fix(prvbut,grd,gbc,1,n++,1,1);

        // Define the common options button

        cmvbut = new JButton("Common Options");
        cmvbut.setBackground(art.butn);
        cmvbut.setForeground(art.butf);
        fix(cmvbut,grd,gbc,1,n++,1,1);

        // Define the Extra Controls button

        eopbut = new JButton("Extra Controls");
        eopbut.setBackground(art.butn);
        eopbut.setForeground(art.butf);
        fix(eopbut,grd,gbc,1,n++,1,1);

        // Define the Radiation Damage button

        rdmbut = new JButton("Radiation Damage Controls");
        rdmbut.setBackground(art.butn);
        rdmbut.setForeground(art.butf);
        fix(rdmbut,grd,gbc,1,n++,1,1);

        // Blank line

        fix(new JLabel(" "),grd,gbc,0,n++,1,1);

        // Define the Close button

        close = new JButton("Close");
        close.setBackground(art.butn);
        close.setForeground(art.butf);
        fix(close,grd,gbc,2,n++,1,1);

        // Register action buttons

        make.addActionListener(this);
        edit.addActionListener(this);
        sysbut.addActionListener(this);
        prvbut.addActionListener(this);
        cmvbut.addActionListener(this);
        rdmbut.addActionListener(this);
        eopbut.addActionListener(this);
        close.addActionListener(this);
    }

    public MakeControl(GUI here) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        home=here;
        println("Activated panel for making CONTROL files");
        job=new MakeControl();
        job.pack();
        job.setVisible(true);
        setParams();
        ttitle.setText(title);
    }

    void setParams() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        // set default values

        title="CONTROL file generated by DL_POLY/java utility";
        temp=295.0;
        press=0.0;
        tstep=0.001;
        rcut=0.0;
        delr=0.0;
        rvdw=0.0;
        rprim=0.0;
        epsq=1.0;
        gamt=0.0;
        taut=1.0;
        taup=1.0;
        fcap=10000.0;
        ewltol=1.E-5;
        qtntol=1.E-5;
        shktol=1.E-5;
        jobtim=100.0;
        tclose=10.0;
        nstrun=10;
        nsteql=10;
        mult=1;
        nstbts=10;
        nstbpo=10;
        nstack=100;
        intsta=10;
        keyens=0;
        keyres=0;
        keyfce=0;
        keyalg=0;
        nstraj=100;
        istraj=10;
        levcon=0;
        allpairs=false;
        lcap=false;
        lzeql=false;
        lrdf=false;
        lprdf=false;
        istrdf=10;
        ltraj=false;
        ltscal=false;
        lzden=false;
        lpzden=false;
        istzden=10;
        lzero=false;
        binsize=0.05;
        densvar=0.0;
        ndump=1000;
        ewldev=1;
        spmeev=1;
        lexclude=false;
        mxquat=100;
        mxshak=250;
        nfold1=1;
        nfold2=1;
        nfold3=1;
        lmetdir=false;
        lvdwdir=false;
        lnoindex=false;
        lnostrict=false;
        lnotopo=false;
        nregauss=0;
        lreplay=false;
        lpslab=false;
        lvdwshift=false;
        ldefects=false;
        lvarstp=false;
        lmsdtmp=false;
        lpseudo=false;
        dstart=0;
        impstp=0;
        dintval=0;
        atom=0;
        nstmsdtmp=0;
        imsdtmp=0;
        ipseudtyp=0;
        defcut=0.0;
        energy=0.0;
        vect1=0.0;
        vect2=0.0;
        vect3=0.0;
        varstp=0.0;
        mindis=0.0;
        maxdis=0.0;
        thick=0.0;
        ptemp=0.0;

    }

    void getParams() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        title=ttitle.getText();
        if(sys != null) getsysValues();
        if(prv != null) getprvValues();
        if(cmv != null) getcmvValues();
        if(eop != null) geteopValues();
        if(rdm != null) getrdmValues();
    }

    int ctrMake() {
        /*
***********************************************************************

dl_poly/java program for constructing the dl_poly CONTROL file
defining the simulation control parameters

copyright - daresbury laboratory
author    - w. smith december 2000

***********************************************************************
         */
        boolean kill,lcut,ldelr,lprim,lrvdw;

        // intitialize system variables

        kill=false;
        lcut   = (rcut>0.0);
        ldelr  = (delr>0.0);
        lprim  = (rprim>0.0);
        lrvdw  = (rvdw>0.0);

        // check internal consistency of data

        if(lvdw)
            keyfce = 2*(keyfce/2);
        else
            keyfce = 2*(keyfce/2)+1;

        if(!lcut) {
            kill = true;
            println("Error - no cutoff specified");
        }

        if(!lrvdw && keyfce%2==1) {
            if(lcut) {
                rvdw=rcut;
            }
            else {
                kill = true;
                println("Error - no vdw cutoff set");
            }
        }

        if(!ldelr) {
            kill = true;
            println("Error - no verlet shell width set");
        }

        if(mult>1) {
            if(!lprim) {
                kill = true;
                println("Error - no primary cutoff set");
            }
        }
        else if(rprim>rcut) {
            kill = true;
            println("Error - primary cutoff too large");
        }

        if(keyens >= 1 && keyens <= 5 && keyens != 3) {
            if(taut <= 0.0) {
                kill=true;
                println("Error - temperature relaxation time incorrect");
            }
        }

        if(keyens >= 6 && keyens <= 21) {
            if(taut <= 0.0) {
                kill=true;
                println("Error - temperature relaxation time incorrect");
            }
            if(taup <= 0.0) {
                kill=true;
                println("Error - pressure relaxation time incorrect");
            }
        }

        if(mult>1) {
            if(rcut-rprim < delr) {
                kill = true;
                println("Error - primary and secondary cutoffs incorrect");
            }
        }

        if(rcut < rvdw) {
            kill = true;
            println("Error - VdW cutoff exceeds general cutoff");
        }

        if(allpairs) {
            if(mult == 1) {
                kill = true;
                println("All pairs must use multiple timestep");
            }
            if(keyfce/2 < 2 || keyfce/2>3) {
                kill = true;
                println("Error - electrostatics incorrect for all pairs");
            }
        }

        if (kill) {
            println("CONTROL file NOT created - see above errors");
            return -1;
        }

        //  open the CONTROL file for output

        fname="CNTROL."+numctr;

        try {
            DataOutputStream outStream = new DataOutputStream(new FileOutputStream(fname));

            outStream.writeBytes(title+"\n"+"\n");

            if(lreplay)
                outStream.writeBytes("replay\n");
            else if(keyalg == 0)
                outStream.writeBytes("integration verlet leapfrog \n");
            else if(keyalg == 1)
                outStream.writeBytes("integration velocity verlet \n");
            if(lzero)
                outStream.writeBytes("zero temperature optimisation");
            outStream.writeBytes("\n");

            outStream.writeBytes("# System state information\n");
            outStream.writeBytes("temperature    "+BML.fmt(temp,8)+"\n");
            outStream.writeBytes("pressure       "+BML.fmt(press,8)+"\n");
            outStream.writeBytes("\n");

            outStream.writeBytes("# Ensemble specification\n");
            if(keyens ==0)
                outStream.writeBytes("ensemble nve \n");
            else if(keyens ==1)
                outStream.writeBytes("ensemble nvt andersen "+BML.fmt(taut,8)+BML.fmt(taup,8)+"\n");
            else if(keyens ==2)
                outStream.writeBytes("ensemble nvt berendsen "+BML.fmt(taut,8)+"\n");
            else if(keyens ==3)
                outStream.writeBytes("ensemble nvt evans \n");
            else if(keyens ==4)
                outStream.writeBytes("ensemble nvt hoover "+BML.fmt(taut,8)+"\n");
            else if(keyens ==5)
                outStream.writeBytes("ensemble nvt langevin "+BML.fmt(taut,8)+"\n");
            else if(keyens ==6)
                outStream.writeBytes("ensemble npt berendsen "+BML.fmt(taut,8)+BML.fmt(taup,8)+"\n");
            else if(keyens ==7)
                outStream.writeBytes("ensemble npt hoover "+BML.fmt(taut,8)+BML.fmt(taup,8)+"\n");
            else if(keyens ==8)
                outStream.writeBytes("ensemble npt langevin "+BML.fmt(taut,8)+BML.fmt(taup,8)+"\n");
            else if(keyens ==9)
                outStream.writeBytes("ensemble npt mtk "+BML.fmt(taut,8)+BML.fmt(taup,8)+"\n");
            else if(keyens ==10)
                outStream.writeBytes("ensemble nst berendsen "+BML.fmt(taut,8)+BML.fmt(taup,8)+"\n");
            else if(keyens ==11)
                outStream.writeBytes("ensemble nst hoover "+BML.fmt(taut,8)+BML.fmt(taup,8)+"\n");
            else if(keyens ==12)
                outStream.writeBytes("ensemble nst langevin "+BML.fmt(taut,8)+BML.fmt(taup,8)+"\n");
            else if(keyens ==13)
                outStream.writeBytes("ensemble nst mtk "+BML.fmt(taut,8)+BML.fmt(taup,8)+"\n");
            else if(keyens ==14)
                outStream.writeBytes("ensemble nst berendsen "+BML.fmt(taut,8)+BML.fmt(taup,8)+" area\n");
            else if(keyens ==15)
                outStream.writeBytes("ensemble nst hoover "+BML.fmt(taut,8)+BML.fmt(taup,8)+" area\n");
            else if(keyens ==16)
                outStream.writeBytes("ensemble nst langevin "+BML.fmt(taut,8)+BML.fmt(taup,8)+" area\n");
            else if(keyens ==17)
                outStream.writeBytes("ensemble nst mtk "+BML.fmt(taut,8)+BML.fmt(taup,8)+" area\n");
            else if(keyens ==18)
                outStream.writeBytes("ensemble nst berendsen "+BML.fmt(taut,8)+BML.fmt(taup,8)+" tension "+BML.fmt(gamt,8)+"\n");
            else if(keyens ==19)
                outStream.writeBytes("ensemble nst hoover "+BML.fmt(taut,8)+BML.fmt(taup,8)+" tension "+BML.fmt(gamt,8)+"\n");
            else if(keyens ==20)
                outStream.writeBytes("ensemble nst langevin "+BML.fmt(taut,8)+BML.fmt(taup,8)+" tension "+BML.fmt(gamt,8)+"\n");
            else if(keyens ==21)
                outStream.writeBytes("ensemble nst mtk "+BML.fmt(taut,8)+BML.fmt(taup,8)+" tension "+BML.fmt(gamt,8)+"\n");
            else if(keyens ==22)
                outStream.writeBytes("ensemble pmf \n");
            outStream.writeBytes("\n");

            if(keyres > 0) outStream.writeBytes("# Restart controls\n");
            if(keyres ==1)
                outStream.writeBytes("restart \n");
            else if(keyres ==2)
                outStream.writeBytes("restart scale \n");
            else if(keyres == 3)
                outStream.writeBytes("restart noscale \n");
            outStream.writeBytes("\n");

            outStream.writeBytes("# Startup options\n");
            if(densvar > 0.0)
                outStream.writeBytes("densvar        "+BML.fmt(densvar,8)+"\n");
            if(nfold1*nfold2*nfold3 > 1) {
                outStream.writeBytes("nfold          "+BML.fmt(nfold1,8)+BML.fmt(nfold2,8)+BML.fmt(nfold3,8)+"\n");
            }
            if(lnoindex)
                outStream.writeBytes("no index\n");
            if(lnostrict)
                outStream.writeBytes("no strict\n");
            if(lnotopo)
                outStream.writeBytes("no topology\n");
            if(lpslab)
                outStream.writeBytes("slab\n");
            outStream.writeBytes("\n");

            outStream.writeBytes("# Simulation controls\n");
            outStream.writeBytes("timestep       "+BML.fmt(tstep,8)+"\n");
            outStream.writeBytes("steps          "+BML.fmt(nstrun,8)+"\n");
            outStream.writeBytes("equilibration  "+BML.fmt(nsteql,8)+"\n");
            outStream.writeBytes("multiple step  "+BML.fmt(mult,8)+"\n");
            if(nregauss > 0)
                outStream.writeBytes("regauss        "+BML.fmt(nregauss,8)+"\n");
            if(nstbts>0)
                outStream.writeBytes("scale          "+BML.fmt(nstbts,8)+"\n");
            if(lzeql)
                outStream.writeBytes("collect \n");
            outStream.writeBytes("print          "+BML.fmt(nstbpo,8)+"\n");
            outStream.writeBytes("stack          "+BML.fmt(nstack,8)+"\n");
            outStream.writeBytes("stats          "+BML.fmt(intsta,8)+"\n");
            if(ndump != 1000)
                outStream.writeBytes("dump           "+BML.fmt(ndump,8)+"\n");
            if(ltraj)
                outStream.writeBytes("trajectory     "+BML.fmt(nstraj,8)+BML.fmt(istraj,8)+BML.fmt(levcon,8)+"\n");
            outStream.writeBytes("\n");

            if(lrdf || lzden)outStream.writeBytes("# RDF/Z-density controls\n");
            if(lrdf)
                outStream.writeBytes("rdf            "+BML.fmt(istrdf,8)+"\n");
            if(lrdf && lprdf)
                outStream.writeBytes("print rdf \n");
            if(lzden)
                outStream.writeBytes("zden           "+BML.fmt(istzden,8)+"\n");
            if(lzden && lpzden)
                outStream.writeBytes("print zden \n");
            if(Math.abs(binsize-0.05) > 1.E-10)
                outStream.writeBytes("binsize        "+BML.fmt(binsize,8)+"\n");
            outStream.writeBytes("\n");

            outStream.writeBytes("# Force cut-offs and controls\n");
            if(rprim>0.0)
                outStream.writeBytes("primary cutoff "+BML.fmt(rprim,8)+"\n");
            outStream.writeBytes("cutoff         "+BML.fmt(rcut,8)+"\n");
            outStream.writeBytes("delr width     "+BML.fmt(delr,8)+"\n");
            outStream.writeBytes("rvdw cutoff    "+BML.fmt(rvdw,8)+"\n");
            if(lvdw)
                outStream.writeBytes("no vdw forces\n");
            if(lexclude)
                outStream.writeBytes("exclude\n");
            if(lmetdir)
                outStream.writeBytes("metal direct\n");
            if(lvdwdir)
                outStream.writeBytes("vdw direct\n");
            if(lvdwshift)
                outStream.writeBytes("vdw shift\n");
            if(lcap)
                outStream.writeBytes("cap forces     "+BML.fmt(fcap,8)+"\n");
            if(allpairs)
                outStream.writeBytes("all pairs \n");
            outStream.writeBytes("\n");

            outStream.writeBytes("# Electrostatics controls\n");
            if(keyfce/2 ==0)
                outStream.writeBytes("no electrostatics \n");
            else if(keyfce/2 ==1) {
                outStream.writeBytes("ewald precision"+BML.fmt(ewltol,8)+"\n");
                if(ewldev > 1)
                    outStream.writeBytes("ewald evaluate "+BML.fmt(ewldev,8)+"\n");
            }
            else if(keyfce/2 ==2)
                outStream.writeBytes("distan \n");
            else if(keyfce/2 ==3)
                outStream.writeBytes("coulomb \n");
            else if(keyfce/2 ==4)
                outStream.writeBytes("shift \n");
            else if(keyfce/2 == 5) {
                outStream.writeBytes("reaction field \n");
                outStream.writeBytes("eps constant   "+BML.fmt(epsq,8)+"\n");
            }
            else if(keyfce/2 ==6) {
                outStream.writeBytes("spme precision "+BML.fmt(ewltol,8)+"\n");
                if(spmeev > 1)
                    outStream.writeBytes("spme evaluate  "+BML.fmt(spmeev,8)+"\n");
            }
            else if(keyfce/2 ==7)
                outStream.writeBytes("hke precision  "+BML.fmt(ewltol,8)+"\n");
            outStream.writeBytes("\n");

            if(Math.abs(shktol-1.E-5) > 1.E-10 || Math.abs(shktol-1.E-5) > 1.E-10 || mxshak != 250 || mxquat != 100)
                outStream.writeBytes("# Quaternion and SHAKE controls\n");
            if(Math.abs(shktol-1.E-5) > 1.E-10)
                outStream.writeBytes("shake tolerance"+BML.fmt(shktol,8)+"\n");
            if(mxshak != 250)
                outStream.writeBytes("mxshak         "+BML.fmt(mxshak,8)+"\n");
            if(Math.abs(qtntol-1.E-5) > 1.E-10)
                outStream.writeBytes("quaternion tolerance"+BML.fmt(qtntol,8)+"\n");
            if(mxquat != 100)
                outStream.writeBytes("mxquat         "+BML.fmt(mxquat,8)+"\n");
            outStream.writeBytes("\n");

            if(atom > 0 || ldefects || lvarstp || lmsdtmp || lpseudo)
               outStream.writeBytes("# Radiation damage controls\n");
            if(atom > 0)
                outStream.writeBytes("impact         "+BML.fmt(atom,8)+BML.fmt(impstp,8)+BML.fmt(energy,8)+
                   BML.fmt(vect1,8)+BML.fmt(vect2,8)+BML.fmt(vect3,8)+"\n");
            if(lvarstp) {
                outStream.writeBytes("variable timestep"+BML.fmt(varstp,8)+"\n");
                outStream.writeBytes("mindis           "+BML.fmt(mindis,8)+"\n");
                outStream.writeBytes("maxdis           "+BML.fmt(maxdis,8)+"\n");
            }
            if(ldefects)
                outStream.writeBytes("defects        "+BML.fmt(dstart,8)+BML.fmt(dintval,8)+BML.fmt(defcut,8)+"\n");
            if(lmsdtmp)
                outStream.writeBytes("msdtmp         "+BML.fmt(nstmsdtmp,8)+BML.fmt(imsdtmp,8)+"\n");
            if(lpseudo) {
                outStream.writeBytes("pseudo ");
                if(ipseudtyp == 1)
                    outStream.writeBytes("langevin ");
                else if(ipseudtyp == 2)
                    outStream.writeBytes("direct   ");
                outStream.writeBytes(BML.fmt(thick,8)+BML.fmt(ptemp,8)+"\n");
            }
            outStream.writeBytes("\n");

            outStream.writeBytes("# Job running time and closedown time\n");
            outStream.writeBytes("job time       "+BML.fmt(jobtim,8)+"\n");
            outStream.writeBytes("close time     "+BML.fmt(tclose,8)+"\n");
            outStream.writeBytes("\n");

            outStream.writeBytes("finish \n");

            // close CONTROL file

            outStream.close();
        }
        catch(Exception e) {
            println("error - writing file: "+fname);
            return -2;
        }
        println(fname+" file created");
        numctr++;
        return 0;
    }

    int ctrLoad() {
        /*
***********************************************************************

dl_poly/java program for loading the dl_poly CONTROL file

copyright - daresbury laboratory
author    - w. smith december 2000

**********************************************************************
         */

        int pass=0,n=0;
        String word="",record="",word2="";
        LineNumberReader lnr;

        title="";

        // open CONTROL file

        try {
            lnr = new LineNumberReader(new FileReader(fname));
            println("Reading file: "+fname);
            title=lnr.readLine();
            println("File header record: "+title);
            pass++;

            // scan through the CONTROL file

            while((record=lnr.readLine()) != null) {
                n=BML.countWords(record);
                if(n > 0 && record.charAt(0) != '#') {
                    word=BML.giveWord(record,1).toLowerCase();

                    if(word.indexOf("integ")>=0) {
                        word=BML.giveWord(record,2).toLowerCase();
                        word2=BML.giveWord(record,3).toLowerCase();
                        if(word.indexOf("velocity")>=0 || word2.indexOf("velocity")>=0)
                            keyalg=1;
                        else
                            keyalg=0;
                    }
                    else if(word.indexOf("steps")>=0) {
                        nstrun=BML.giveInteger(record,n);
                    }
                    else if(word.indexOf("equil")>=0) {
                        nsteql=BML.giveInteger(record,n);
                    }
                    else if(word.indexOf("restart")>=0) {
                        if(record.indexOf("noscale")>=0)
                            keyres = 3;
                        else if(record.indexOf("scale")>=0)
                            keyres = 2;
                        else
                            keyres = 1;
                    }
                    else if(word.indexOf("ensemble")>=0) {
                        word=BML.giveWord(record,2).toLowerCase();
                        if(word.indexOf("nve")>=0)
                            keyens = 0;
                        else if(word.indexOf("nvt")>=0) {
                            word=BML.giveWord(record,3).toLowerCase();
                            if(word.indexOf("andersen")>=0)
                                keyens = 1;
                            else if(word.indexOf("berendsen")>=0) {
                                keyens = 2;
                                taut=BML.giveDouble(record,n);
                            }
                            else if(word.indexOf("evans")>=0)
                                keyens = 3;
                            else if(word.indexOf("hoover")>=0) {
                                keyens = 4;
                                taut=BML.giveDouble(record,n);
                            }
                            else if(word.indexOf("langevin")>=0) {
                                keyens = 5;
                                taut=BML.giveDouble(record,n);
                            }
                        }
                        else if(word.indexOf("npt")>=0) {
                            word=BML.giveWord(record,3).toLowerCase();
                            if(word.indexOf("berendsen")>=0) {
                                keyens = 6;
                                taut = BML.giveDouble(record,n-1);
                                taup = BML.giveDouble(record,n);
                            }
                            else if(word.indexOf("hoover")>=0) {
                                keyens = 7;
                                taut = BML.giveDouble(record,n-1);
                                taup = BML.giveDouble(record,n);
                            }
                            else if(word.indexOf("langevin")>=0) {
                                keyens = 8;
                                taut = BML.giveDouble(record,n-1);
                                taup = BML.giveDouble(record,n);
                            }
                            else if(word.indexOf("mtk")>=0) {
                                keyens = 9;
                                taut = BML.giveDouble(record,n-1);
                                taup = BML.giveDouble(record,n);
                            }
                        }
                        else if(word.indexOf("nst")>=0) {
                            word=BML.giveWord(record,3).toLowerCase();
                            taut = BML.giveDouble(record,4);
                            taup = BML.giveDouble(record,5);
                            if(word.indexOf("berendsen")>=0) {
                                keyens = 10;
                                word=BML.giveWord(record,6).toLowerCase();
                                if(word.indexOf("area")>=0)
                                    keyens = 14;
                                else if(word.indexOf("tens")>=0) {
                                    keyens = 18;
                                    gamt = BML.giveDouble(record,n);
                                }
                            }
                            else if(word.indexOf("hoover")>=0) {
                                keyens = 11;
                                word=BML.giveWord(record,6).toLowerCase();
                                if(word.indexOf("area")>=0)
                                    keyens = 15;
                                else if(word.indexOf("tens")>=0) {
                                    keyens = 19;
                                    gamt = BML.giveDouble(record,n);
                                }
                            }
                            else if(word.indexOf("langevin")>=0) {
                                keyens = 12;
                                word=BML.giveWord(record,6).toLowerCase();
                                if(word.indexOf("area")>=0)
                                    keyens = 16;
                                else if(word.indexOf("tens")>=0) {
                                    keyens = 20;
                                    gamt = BML.giveDouble(record,n);
                                }
                            }
                            else if(word.indexOf("mtk")>=0) {
                                keyens = 13;
                                word=BML.giveWord(record,6).toLowerCase();
                                if(word.indexOf("area")>=0)
                                    keyens = 17;
                                else if(word.indexOf("tens")>=0) {
                                    keyens = 21;
                                    gamt = BML.giveDouble(record,n);
                                }
                            }
                        }
                        else if(word.indexOf("pmf")>=0) {
                            keyens = 22;
                        }
                    }
                    else if(word.indexOf("scale")>=0) {
                        ltscal=true;
                        nstbts=BML.giveInteger(record,n);
                    }
                    else if(word.indexOf("rdf")>=0) {
                        lrdf=true;
                        istrdf=BML.giveInteger(record,n);
                    }
                    else if(word.indexOf("zden")>=0) {
                        lzden=true;
                        istzden=BML.giveInteger(record,n);
                    }
                    else if(word.indexOf("collect")>=0) {
                        lzeql=true;
                    }
                    else if(word.indexOf("zero")>=0) {
                        lzero=true;
                    }
                    else if(word.indexOf("print")>=0) {
                        word=BML.giveWord(record,2).toLowerCase();
                        if(word.indexOf("rdf")>=0) {
                            lprdf=true;
                        }
                        else if(word.indexOf("zden")>=0) {
                            lpzden=true;
                        }
                        else {
                            nstbpo=BML.giveInteger(record,n);
                            nstbpo = Math.max(nstbpo,1);
                        }
                    }
                    else if(word.indexOf("regauss")>=0) {
                        nregauss=BML.giveInteger(record,n);
                    }
                    else if(word.indexOf("stack")>=0) {
                        nstack=BML.giveInteger(record,n);
                    }
                    else if(word.indexOf("stats")>=0) {
                        intsta=BML.giveInteger(record,n);
                    }
                    else if(word.indexOf("traj")>=0) {
                        ltraj=true;
                        nstraj=BML.giveInteger(record,2);
                        istraj=BML.giveInteger(record,3);
                        levcon=BML.giveInteger(record,4);
                    }
                    else if(word.indexOf("ewald")>=0) {
                        word=BML.giveWord(record,2).toLowerCase();
                        if(word.indexOf("eval")>=0) {
                            ewldev=BML.giveInteger(record,n);
                        }
                        else {
                            keyfce = 2;
                            ewltol=BML.giveDouble(record,n);
                        }
                    }
                    else if(word.indexOf("distan")>=0) {
                        keyfce = 4;
                    }
                    else if(word.indexOf("coul")>=0) {
                        keyfce = 6;
                    }
                    else if(word.indexOf("shift")>=0) {
                        keyfce = 8;
                    }
                    else if(word.indexOf("reaction")>=0) {
                        keyfce = 10;
                    }
                    else if(word.indexOf("spme")>=0) {
                        word=BML.giveWord(record,2).toLowerCase();
                        if(word.indexOf("eval")>=0) {
                            spmeev=BML.giveInteger(record,n);
                        }
                        else {
                            keyfce = 12;
                            ewltol=BML.giveDouble(record,n);
                        }
                    }
                    else if(word.indexOf("hke")>=0) {
                        keyfce = 14;
                        ewltol=BML.giveDouble(record,n);
                    }
                    else if(word.indexOf("cap")>=0) {
                        lcap = true;
                        fcap=BML.giveDouble(record,n);
                    }
                    else if(word.indexOf("exclu")>=0) {
                        lexclude=true;
                    }
                    else if(word.indexOf("no")>=0) {
                        word=BML.giveWord(record,2).toLowerCase();
                        if (word.indexOf("vdw")>=0)
                            lvdw=true;
                        else if(word.indexOf("elec")>=0)
                            keyfce=0;
                        else if(word.indexOf("str")>=0)
                            lnostrict=true;
                        else if(word.indexOf("topo")>=0)
                            lnotopo=true;
                        else if(word.indexOf("ind")>=0)
                            lnoindex=true;
                    }
                    else if(word.indexOf("nfold")>=0) {
                        nfold1=BML.giveInteger(record,n-2);
                        nfold2=BML.giveInteger(record,n-1);
                        nfold3=BML.giveInteger(record,n);
                    }
                    else if(word.indexOf("mult")>=0) {
                        mult=BML.giveInteger(record,n);
                    }
                    else if(word.indexOf("mxquat")>=0) {
                        mxquat=BML.giveInteger(record,n);
                    }
                    else if(word.indexOf("mxshak")>=0) {
                        mxshak=BML.giveInteger(record,n);
                    }
                    else if(word.indexOf("timestep")>=0) {
                        tstep=BML.giveDouble(record,n);
                    }
                    else if(word.indexOf("temp")>=0) {
                        temp=BML.giveDouble(record,n);
                    }
                    else if(word.indexOf("pres")>=0) {
                        press=BML.giveDouble(record,n);
                    }
                    else if(word.indexOf("prim")>=0) {
                        rprim=BML.giveDouble(record,n);
                    }
                    else if(word.indexOf("cut")>=0) {
                        rcut=BML.giveDouble(record,n);
                    }
                    else if(word.indexOf("rvdw")>=0) {
                        rvdw=BML.giveDouble(record,n);
                    }
                    else if(word.indexOf("metal")>=0) {
                        word=BML.giveWord(record,2).toLowerCase();
                        if(word.indexOf("direct")>=0)
                            lmetdir=true;
                    }
                    else if(word.indexOf("vdw")>=0) {
                        word=BML.giveWord(record,2).toLowerCase();
                        if(word.indexOf("shift")>=0)
                            lvdwshift=true;
                        else if(word.indexOf("direct")>=0)
                            lvdwdir=true;
                    }
                    else if(word.indexOf("delr")>=0) {
                        delr=BML.giveDouble(record,n);
                    }
                    else if(word.indexOf("eps")>=0) {
                        epsq=BML.giveDouble(record,n);
                    }
                    else if(word.indexOf("shake")>=0) {
                        shktol=BML.giveDouble(record,n);
                    }
                    else if(word.indexOf("quat")>=0) {
                        qtntol=BML.giveDouble(record,n);
                    }
                    else if(word.indexOf("job")>=0) {
                        jobtim=BML.giveDouble(record,n);
                    }
                    else if(word.indexOf("close")>=0) {
                        tclose=BML.giveDouble(record,n);
                    }
                    else if(word.indexOf("replay")>=0) {
                        lreplay=true;
                    }
                    else if(word.indexOf("slab")>=0) {
                        lpslab=true;
                    }
                    else if(word.indexOf("all")>=0) {
                        allpairs = true;
                    }
                    else if(word.indexOf("binsize")>=0) {
                        binsize=BML.giveDouble(record,n);
                    }
                    else if(word.indexOf("densvar")>=0) {
                        densvar=BML.giveDouble(record,n);
                    }
                    else if(word.indexOf("dump")>=0) {
                        ndump=BML.giveInteger(record,n);
                    }
                    else if(word.indexOf("impact")>=0) {
                        atom=BML.giveInteger(record,n-5);
                        impstp=BML.giveInteger(record,n-4);
                        energy=BML.giveDouble(record,n-3);
                        vect1=BML.giveDouble(record,n-2);
                        vect2=BML.giveDouble(record,n-1);
                        vect3=BML.giveDouble(record,n);
                    }
                    else if(word.indexOf("variable")>=0) {
                        word=BML.giveWord(record,2).toLowerCase();
                        if(word.indexOf("timestep")>=0) {
                           varstp=BML.giveDouble(record,n);
                           lvarstp=true;
                        }
                    }
                    else if(word.indexOf("mindis")>=0) {
                        mindis=BML.giveDouble(record,n);
                    }
                    else if(word.indexOf("maxdis")>=0) {
                        maxdis=BML.giveDouble(record,n);
                    }
                    else if(word.indexOf("defects")>=0) {
                        dstart=BML.giveInteger(record,n-2);
                        dintval=BML.giveInteger(record,n-1);
                        defcut=BML.giveDouble(record,n);
                        ldefects=true;
                    }
                    else if(word.indexOf("msdtmp")>=0) {
                        nstmsdtmp=BML.giveInteger(record,n-1);
                        imsdtmp=BML.giveInteger(record,n);
                        lmsdtmp=true;
                    }
                    else if(word.indexOf("pseudo")>=0) {
                        word=BML.giveWord(record,2).toLowerCase();
                        if(word.indexOf("langevin")>=0)
                           ipseudtyp=1;
                        else if(word.indexOf("direct")>=0)
                            ipseudtyp=2;
                        else
                            ipseudtyp=0;
                        thick=BML.giveDouble(record,n-1);
                        ptemp=BML.giveDouble(record,n);
                        lpseudo=true;
                    }
                    else if(word.indexOf("finish")>=0) {
                        pass++;
                        println("File "+fname+" processed lines "+pass);
                    }
                }
                pass++;
            }
            lnr.close();
        }
        catch(FileNotFoundException e) {
            println("Error - file not found: " + fname);
            return -1;
        }
        catch(Exception e) {
            println("Error reading file: " + fname + " "+e);
            return -2;
        }
        return 0;
    }

    void setsysValues(){
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */

        sys.temp=temp;
        sys.press=press;
        sys.tstep=tstep;
        sys.rcut=rcut;
        sys.delr=delr;
        sys.rvdw=rvdw;
        sys.rprim=rprim;
        sys.gamt=gamt;
        sys.epsq=epsq;
        sys.taut=taut;
        sys.taup=taup;
        sys.keyens=keyens;
        sys.keyfce=keyfce;
    }

    void setprvValues(){
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */

        prv.keyalg=keyalg;
        prv.mult=mult;
        prv.nstbpo=nstbpo;
        prv.nstrun=nstrun;
        prv.nsteql=nsteql;
        prv.nstack=nstack;
        prv.intsta=intsta;
        prv.ewltol=ewltol;
        prv.shktol=shktol;
        prv.qtntol=qtntol;
        prv.jobtim=jobtim;
        prv.tclose=tclose;
        prv.keyres=keyres;

    }

    void setcmvValues(){
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */

        cmv.allpairs=allpairs;
        cmv.lcap=lcap;
        cmv.fcap=fcap;
        cmv.lzeql=lzeql;
        cmv.lvdw=lvdw;
        cmv.lprdf=lprdf;
        cmv.lrdf=lrdf;
        cmv.istrdf=istrdf;
        cmv.ltraj=ltraj;
        cmv.nstraj=nstraj;
        cmv.istraj=istraj;
        cmv.levcon=levcon;
        cmv.ltscal=ltscal;
        cmv.nstbts=nstbts;
        cmv.lpzden=lpzden;
        cmv.lzden=lzden;
        cmv.istzden=istzden;
        cmv.lzero=lzero;

    }

    void seteopValues(){
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */

        eop.binsize=binsize;
        eop.densvar=densvar;
        eop.ndump=ndump;
        eop.ewldev=ewldev;
        eop.spmeev=spmeev;
        eop.lexclude=lexclude;
        eop.mxquat=mxquat;
        eop.mxshak=mxshak;
        eop.nfold1=nfold1;
        eop.nfold2=nfold2;
        eop.nfold3=nfold3;
        eop.lmetdir=lmetdir;
        eop.lvdwdir=lvdwdir;
        eop.lnoindex=lnoindex;
        eop.lnostrict=lnostrict;
        eop.lnotopo=lnotopo;
        eop.nregauss=nregauss;
        eop.lreplay=lreplay;
        eop.lpslab=lpslab;
        eop.lvdwshift=lvdwshift;

    }

    void setrdmValues(){
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */

        rdm.ldefects=ldefects;
        rdm.lvarstp=lvarstp;
        rdm.lmsdtmp=lmsdtmp;
        rdm.lpseudo=lpseudo;
        rdm.dstart=dstart;
        rdm.impstp=impstp;
        rdm.dintval=dintval;
        rdm.atom=atom;
        rdm.nstmsdtmp=nstmsdtmp;
        rdm.imsdtmp=imsdtmp;
        rdm.ipseudtyp=ipseudtyp;
        rdm.defcut=defcut;
        rdm.energy=energy;
        rdm.vect1=vect1;
        rdm.vect2=vect2;
        rdm.vect3=vect3;
        rdm.varstp=varstp;
        rdm.mindis=mindis;
        rdm.maxdis=maxdis;
        rdm.thick=thick;
        rdm.ptemp=ptemp;

    }

    void getsysValues(){
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */

        temp=BML.giveDouble(sys.ttemp.getText(),1);
        press=BML.giveDouble(sys.tpress.getText(),1);
        tstep=BML.giveDouble(sys.ttstep.getText(),1);
        rcut=BML.giveDouble(sys.trcut.getText(),1);
        delr=BML.giveDouble(sys.tdelr.getText(),1);
        rvdw=BML.giveDouble(sys.trvdw.getText(),1);
        rprim=BML.giveDouble(sys.trprim.getText(),1);
        epsq=BML.giveDouble(sys.tepsq.getText(),1);
        taut=BML.giveDouble(sys.ttaut.getText(),1);
        taup=BML.giveDouble(sys.ttaup.getText(),1);
        keyens=sys.ensemble.getSelectedIndex();
        keyfce=2*sys.electro.getSelectedIndex();
    }

    void getprvValues(){
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */

        keyalg=prv.algorithm.getSelectedIndex();
        ewltol=BML.giveDouble(prv.tewltol.getText(),1);
        qtntol=BML.giveDouble(prv.tqtntol.getText(),1);
        shktol=BML.giveDouble(prv.tshktol.getText(),1);
        jobtim=BML.giveDouble(prv.tjobtim.getText(),1);
        tclose=BML.giveDouble(prv.ttclose.getText(),1);
        nstrun=BML.giveInteger(prv.tnstrun.getText(),1);
        nsteql=BML.giveInteger(prv.tnsteql.getText(),1);
        mult=BML.giveInteger(prv.tmult.getText(),1);
        nstbpo=BML.giveInteger(prv.tnstbpo.getText(),1);
        nstack=BML.giveInteger(prv.tnstack.getText(),1);
        intsta=BML.giveInteger(prv.tintsta.getText(),1);
        keyres=prv.restopt.getSelectedIndex();

    }

    void getcmvValues(){
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */

        istrdf=BML.giveInteger(cmv.tistrdf.getText(),1);
        istzden=BML.giveInteger(cmv.tistzden.getText(),1);
        nstbts=BML.giveInteger(cmv.tnstbts.getText(),1);
        nstraj=BML.giveInteger(cmv.tnstraj.getText(),1);
        istraj=BML.giveInteger(cmv.tistraj.getText(),1);
        levcon=BML.giveInteger(cmv.tlevcon.getText(),1);
        fcap=BML.giveDouble(cmv.tfcap.getText(),1);
        allpairs=cmv.ballpairs.isSelected();
        lcap=cmv.blcap.isSelected();
        lvdw=cmv.blvdw.isSelected();
        lzeql=cmv.blzeql.isSelected();
        lrdf=cmv.blrdf.isSelected();
        lprdf=cmv.blprdf.isSelected();
        ltraj=cmv.bltraj.isSelected();
        ltscal=cmv.bltscal.isSelected();
        lzden=cmv.blzden.isSelected();
        lpzden=cmv.blpzden.isSelected();
        lzero=cmv.blzero.isSelected();

    }

    void geteopValues(){
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */

        binsize=BML.giveDouble(eop.tbinsize.getText(),1);
        densvar=BML.giveDouble(eop.tdensvar.getText(),1);
        ndump=BML.giveInteger(eop.tdump.getText(),1);
        ewldev=BML.giveInteger(eop.tewldev.getText(),1);
        spmeev=BML.giveInteger(eop.tspmeev.getText(),1);
        lexclude=eop.blexclude.isSelected();
        mxquat=BML.giveInteger(eop.tmxquat.getText(),1);
        mxshak=BML.giveInteger(eop.tmxshak.getText(),1);
        nfold1=BML.giveInteger(eop.tnfold1.getText(),1);
        nfold2=BML.giveInteger(eop.tnfold2.getText(),1);
        nfold3=BML.giveInteger(eop.tnfold3.getText(),1);
        lmetdir=eop.blmetdir.isSelected();
        lvdwdir=eop.blvdwdir.isSelected();
        lnoindex=eop.blnoindex.isSelected();
        lnostrict=eop.blnostrict.isSelected();
        lnotopo=eop.blnotopo.isSelected();
        nregauss=BML.giveInteger(eop.tregauss.getText(),1);
        lreplay=eop.blreplay.isSelected();
        lpslab=eop.blpslab.isSelected();
        lvdwshift=eop.blvdwshift.isSelected();

    }

    void getrdmValues(){
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */

        ldefects=rdm.bldefects.isSelected();
        dstart=BML.giveInteger(rdm.tdstart.getText(),1);
        dintval=BML.giveInteger(rdm.tdintval.getText(),1);
        defcut=BML.giveDouble(rdm.tdefcut.getText(),1);
        atom=BML.giveInteger(rdm.tatom.getText(),1);
        impstp=BML.giveInteger(rdm.timpstp.getText(),1);
        energy=BML.giveDouble(rdm.tenergy.getText(),1);
        vect1=BML.giveDouble(rdm.tvect1.getText(),1);
        vect2=BML.giveDouble(rdm.tvect2.getText(),1);
        vect3=BML.giveDouble(rdm.tvect3.getText(),1);
        lvarstp=rdm.blvarstp.isSelected();
        varstp=BML.giveDouble(rdm.tvarstp.getText(),1);
        mindis=BML.giveDouble(rdm.tmindis.getText(),1);
        maxdis=BML.giveDouble(rdm.tmaxdis.getText(),1);
        lmsdtmp=rdm.blmsdtmp.isSelected();
        nstmsdtmp=BML.giveInteger(rdm.tnstmsdtmp.getText(),1);
        imsdtmp=BML.giveInteger(rdm.timsdtmp.getText(),1);
        lpseudo=rdm.blpseudo.isSelected();
        ipseudtyp=rdm.pseudtyp.getSelectedIndex();
        thick=BML.giveDouble(rdm.tthick.getText(),1);
        ptemp=BML.giveDouble(rdm.tptemp.getText(),1);

    }

    public void actionPerformed(ActionEvent evnt) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        int call;
        Object origin=evnt.getSource();

        if (origin == make) {
            getParams();
            call=ctrMake();
        }
        else if (origin == edit) {
            if((fname=selectFileNameBegins(home,"CNT"))!=null) {
                call=ctrLoad();
                if(sys != null) sys.setParams();
                if(prv != null) prv.setParams();
                if(cmv != null) cmv.setParams();
                if(eop != null) eop.setParams();
            }
            else {
                println("File selection cancelled");
            }
        }
        else if (origin == close) {
            byebye();
        }
        else if (origin == sysbut) {
            if(sys == null)
                sys=new SysVar(job);
        }
        else if (origin == prvbut) {
            if(prv == null)
                prv=new ProVar(job);
        }
        else if (origin == cmvbut) {
            if(cmv == null)
                cmv=new ComVar(job);
        }
        else if (origin == eopbut) {
            if(eop == null)
                eop=new EopVar(job);
        }
        else if (origin == rdmbut) {
            if(rdm == null)
                rdm=new RadDam(job);
        }
    }

    static void byebye() {
        if(sys != null) sys.byebye();
        if(prv != null) prv.byebye();
        if(cmv != null) cmv.byebye();
        if(eop != null) eop.byebye();
        if(rdm != null) rdm.byebye();
        job.setVisible(false);
    }

}

