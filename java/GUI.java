import java.awt.*;
import java.io.*;
import java.awt.event.*;
import java.awt.print.*;
import java.awt.geom.*;
import javax.swing.*;

// Define the DL_POLY Graphical User Interface

public class GUI extends Basic implements Printable{
        /*
*********************************************************************

dl_poly/java GUI class

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
    private static GUI job;
    private static String page="Welcome to the DL_POLY GUI (c) W. Smith 2011.";

    public GUI() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */

        // Set up Frame

        super();
        setTitle("DL_POLY GUI");
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        fontMain=new Font("Verdana",Font.BOLD,14);
        setFont(fontMain);
        grd = new GridLayout(1,1);
        getContentPane().setLayout(grd);
        board = new JTextArea();
        board.setBackground(art.scrn);
        board.setForeground(art.scrf);
        board.setFont(fontMain);
        board.setEditable(false);
        board.setLineWrap(true);
        board.setWrapStyleWord(true);
        scroll = new JScrollPane(board);
        scroll.setPreferredSize(new Dimension(SCREENX,SCREENY));
        getContentPane().add(scroll);

        // Define menu bar

        defMenuBar();

    }

    // Define the menu bar

    void defMenuBar() {
        /*
*********************************************************************

dl_poly/java GUI routine to define the graphics buttons

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        JMenuBar top = new JMenuBar();
        top.setForeground(art.fore);
        top.setBackground(art.back);
        top.setBorderPainted(false);

        // File menu A

        JMenuItem itemAAa,itemAAb,itemAAc,itemAAd;
        JMenuItem itemAB,itemAC,itemAD,itemAE,itemAF,itemAG,itemAH,itemAI;
        JMenu File = MyMenu("File");
        JMenu subFileAA = MyMenu("Load Config");
        subFileAA.add(itemAAa = MyMenuItem("CFG"));
        subFileAA.add(itemAAb = MyMenuItem("XYZ"));
        subFileAA.add(itemAAc = MyMenuItem("PDB"));
        subFileAA.add(itemAAd = MyMenuItem("MSI"));
        subFileAA.setForeground(Color.black);
        File.add(subFileAA);
        File.add(itemAB = MyMenuItem("Delete Config"));
        File.add(itemAC = MyMenuItem("View File"));
        File.add(itemAD = MyMenuItem("Delete File"));
        File.add(itemAE = MyMenuItem("Defaults"));
        File.add(itemAF = MyMenuItem("Print Text"));
        File.add(itemAG = MyMenuItem("Print Graphics"));
        File.add(itemAH = MyMenuItem("Reset"));
        File.add(itemAI = MyMenuItem("Quit"));
        File.setFont(fontMain);
        File.setForeground(art.fore);
        File.setBackground(art.back);
        top.add(File);
        itemAAa.addActionListener(new GUIMenuHandler());
        itemAAb.addActionListener(new GUIMenuHandler());
        itemAAc.addActionListener(new GUIMenuHandler());
        itemAAd.addActionListener(new GUIMenuHandler());
        itemAB.addActionListener(new GUIMenuHandler());
        itemAC.addActionListener(new GUIMenuHandler());
        itemAD.addActionListener(new GUIMenuHandler());
        itemAE.addActionListener(new GUIMenuHandler());
        itemAF.addActionListener(new GUIMenuHandler());
        itemAG.addActionListener(new GUIMenuHandler());
        itemAH.addActionListener(new GUIMenuHandler());
        itemAI.addActionListener(new GUIMenuHandler());

        // FileMaker menu C

        JMenuItem itemFA,itemCA;
        JMenuItem itemCBa,itemCBb,itemCBc,itemCBd;
        JMenuItem itemCCa,itemCCb,itemCCc,itemCCd,itemCCe;
        JMenuItem itemCDa,itemCDb,itemCDc,itemCDd,itemCDe;
        JMenu FileMaker  = MyMenu("FileMaker");
        FileMaker.add(itemFA = MyMenuItem("Show Editor"));
        FileMaker.add(itemCA = MyMenuItem("CONTROL"));
        JMenu subFileMakerCB = MyMenu("CONFIG");
        subFileMakerCB.add(itemCBa = MyMenuItem("Lattice"));
        subFileMakerCB.add(itemCBb = MyMenuItem("Chain"));
        subFileMakerCB.add(itemCBc = MyMenuItem("Polymer"));
        subFileMakerCB.add(itemCBd = MyMenuItem("Bucky"));
        subFileMakerCB.setForeground(Color.black);
        FileMaker.add(subFileMakerCB);
        JMenu subFileMakerCC = MyMenu("FIELD");
        subFileMakerCC.add(itemCCa = MyMenuItem("Blank"));
        subFileMakerCC.add(itemCCb = MyMenuItem("Dreiding"));
        subFileMakerCC.add(itemCCc = MyMenuItem("OPLS"));
        subFileMakerCC.add(itemCCd = MyMenuItem("Ceramics"));
        subFileMakerCC.add(itemCCe = MyMenuItem("Table"));
        subFileMakerCC.setForeground(Color.black);
        FileMaker.add(subFileMakerCC);
        JMenu subFileMakerCD = MyMenu("Tools");
        subFileMakerCD.add(itemCDa = MyMenuItem("Add Solvent"));
        subFileMakerCD.add(itemCDb = MyMenuItem("BondLengths"));
        subFileMakerCD.add(itemCDc = MyMenuItem("Insert Molecule"));
        subFileMakerCD.add(itemCDd = MyMenuItem("N_fold"));
        subFileMakerCD.add(itemCDe = MyMenuItem("Slice"));
        subFileMakerCD.setForeground(Color.black);
        FileMaker.add(subFileMakerCD);
        FileMaker.setForeground(art.fore);
        FileMaker.setBackground(art.back);
        FileMaker.setFont(fontMain);
        top.add(FileMaker);
        itemFA.addActionListener(new GUIMenuHandler());
        itemCA.addActionListener(new GUIMenuHandler());
        itemCBa.addActionListener(new GUIMenuHandler());
        itemCBb.addActionListener(new GUIMenuHandler());
        itemCBc.addActionListener(new GUIMenuHandler());
        itemCBd.addActionListener(new GUIMenuHandler());
        itemCCa.addActionListener(new GUIMenuHandler());
        itemCCb.addActionListener(new GUIMenuHandler());
        itemCCc.addActionListener(new GUIMenuHandler());
        itemCCd.addActionListener(new GUIMenuHandler());
        itemCCe.addActionListener(new GUIMenuHandler());
        itemCDa.addActionListener(new GUIMenuHandler());
        itemCDb.addActionListener(new GUIMenuHandler());
        itemCDc.addActionListener(new GUIMenuHandler());
        itemCDd.addActionListener(new GUIMenuHandler());
        itemCDe.addActionListener(new GUIMenuHandler());

        // Execute menu D

        JMenuItem itemDA,itemDB;
        JMenu Execute  = MyMenu("Execute");
        Execute.add(itemDA = MyMenuItem("Run DL_POLY"));
        Execute.add(itemDB = MyMenuItem("Store/Fetch Data"));
        Execute.setForeground(art.fore);
        Execute.setBackground(art.back);
        Execute.setFont(fontMain);
        top.add(Execute);
        itemDA.addActionListener(new GUIMenuHandler());
        itemDB.addActionListener(new GUIMenuHandler());

        // Analysis menu E

        JMenuItem itemE0;
        JMenuItem itemEAa,itemEAb,itemEAc,itemEAd;
        JMenuItem itemEBa,itemEBb,itemEBc;
        JMenuItem itemECa,itemECb,itemECc;
        JMenuItem itemEDa,itemEDb,itemEDc;
        JMenuItem itemEEa,itemEEb;
        JMenu Analysis  = MyMenu("Analysis");
        Analysis.add(itemE0 = MyMenuItem("Statistics"));
        JMenu subAnalysisEA = MyMenu("Structure");
        subAnalysisEA.add(itemEAa = MyMenuItem("RDF_Plot"));
        subAnalysisEA.add(itemEAb = MyMenuItem("RDF_Calc"));
        subAnalysisEA.add(itemEAc = MyMenuItem("S(k)"));
        subAnalysisEA.add(itemEAd = MyMenuItem("Z_Density"));
        subAnalysisEA.setForeground(Color.black);
        Analysis.add(subAnalysisEA);
        JMenu subAnalysisEB = MyMenu("Dynamics");
        subAnalysisEB.add(itemEBa = MyMenuItem("MSD"));
        subAnalysisEB.add(itemEBb = MyMenuItem("VAF"));
        subAnalysisEB.add(itemEBc = MyMenuItem("FAF"));
        subAnalysisEB.setForeground(Color.black);
        Analysis.add(subAnalysisEB);
        JMenu subAnalysisEC = MyMenu("van Hove");
        subAnalysisEC.add(itemECa = MyMenuItem("Gs(r,t)"));
        subAnalysisEC.add(itemECb = MyMenuItem("Gd(r,t)"));
        subAnalysisEC.add(itemECc = MyMenuItem("S(k,w)"));
        subAnalysisEC.setForeground(Color.black);
        Analysis.add(subAnalysisEC);
        JMenu subAnalysisED = MyMenu("Display");
        subAnalysisED.add(itemEDa = MyMenuItem("CONFIG"));
        subAnalysisED.add(itemEDb = MyMenuItem("REVCON"));
        subAnalysisED.setForeground(Color.black);
        Analysis.add(subAnalysisED);
        JMenu subAnalysisEE = MyMenu("Tools");
        subAnalysisEE.add(itemEEa = MyMenuItem("What Atoms?"));
        subAnalysisEE.add(itemEEb = MyMenuItem("Plot"));
        subAnalysisEE.setForeground(Color.black);
        Analysis.add(subAnalysisEE);
        Analysis.setForeground(art.fore);
        Analysis.setBackground(art.back);
        Analysis.setFont(fontMain);
        top.add(Analysis);
        itemE0.addActionListener(new GUIMenuHandler());
        itemEAa.addActionListener(new GUIMenuHandler());
        itemEAb.addActionListener(new GUIMenuHandler());
        itemEAc.addActionListener(new GUIMenuHandler());
        itemEAd.addActionListener(new GUIMenuHandler());
        itemEBa.addActionListener(new GUIMenuHandler());
        itemEBb.addActionListener(new GUIMenuHandler());
        itemEBc.addActionListener(new GUIMenuHandler());
        itemECa.addActionListener(new GUIMenuHandler());
        itemECb.addActionListener(new GUIMenuHandler());
        itemECc.addActionListener(new GUIMenuHandler());
        itemEDa.addActionListener(new GUIMenuHandler());
        itemEDb.addActionListener(new GUIMenuHandler());
        itemEEa.addActionListener(new GUIMenuHandler());
        itemEEb.addActionListener(new GUIMenuHandler());

        // Information menu G

        JMenuItem itemGA,itemGB,itemGC,itemGD,itemGE,itemGF;
        JMenuItem itemGG,itemGH,itemGI;
        JMenu Help  = MyMenu("Information");
        Help.add(itemGA = MyMenuItem("About DL_POLY"));
        Help.add(itemGB = MyMenuItem("Disclaimer"));
        Help.add(itemGC = MyMenuItem("Licence"));
        Help.add(itemGD = MyMenuItem("Acknowledgements"));
        Help.add(itemGE = MyMenuItem("MINIDREI"));
        Help.add(itemGF = MyMenuItem("MINIOPLS"));
        Help.add(itemGG = MyMenuItem("CERAMICS"));
        Help.add(itemGH = MyMenuItem("Test Cases"));
        Help.add(itemGI = MyMenuItem("Clear Text"));
        Help.setForeground(art.fore);
        Help.setBackground(art.back);
        Help.setFont(fontMain);
        top.add(Help);
        itemGA.addActionListener(new GUIMenuHandler());
        itemGB.addActionListener(new GUIMenuHandler());
        itemGC.addActionListener(new GUIMenuHandler());
        itemGD.addActionListener(new GUIMenuHandler());
        itemGE.addActionListener(new GUIMenuHandler());
        itemGF.addActionListener(new GUIMenuHandler());
        itemGG.addActionListener(new GUIMenuHandler());
        itemGH.addActionListener(new GUIMenuHandler());
        itemGI.addActionListener(new GUIMenuHandler());

        // Invoke menu bar

        setJMenuBar(top);
    }

    class GUIMenuHandler implements ActionListener {
        /*
*********************************************************************

dl_poly/java GUI class to handle menu events

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        public void actionPerformed(ActionEvent e) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
            String arg = e.getActionCommand();

            if (e.getSource() instanceof JMenuItem) {
                if(arg.equals("Quit")) {
                    System.exit(0);
                }
                else if(arg.equals("View File")) {
                    getViewFile(job);
                }
                else if(arg.equals("Delete Config")) {
                    config=null;
		    cfgsav=null;
		    if(editor != null)
			editor.closeEditor();
                }
                else if(arg.equals("Delete File")) {
                    zappFile(job);
                }
                else if(arg.equals("Defaults")) {
                    if(newdefs != null)
                        newdefs.job.setVisible(false);
                    newdefs=new ChgDefaults(job);
                }
                else if(arg.equals("Reset")) {
                    startGUI();
                    clearScreen();
		    config=null;
		    if(editor != null)
			editor.closeEditor();
                    println("GUI variables reinitialized");
                }
                else if(arg.equals("Print Text")) {
                    printOut();
                }
                else if(arg.equals("Print Graphics")) {
		    if(config == null)
			println("Error - No configuration to print");
		    else
			editor.pane.printOut();
                }
                else if(arg.equals("Show Editor")) {
                    if(editor != null)
                        editor.closeEditor();
		    editor=new Editor(job);
                }
                else if(arg.equals("CONTROL")){
                    if(makctr != null)
                        makctr.byebye();
                    makctr=new MakeControl(job);
                }
                else if(arg.equals("CFG")) {
                    ftype="CFG";
                    println("CFG file input selected");
		    config=getConfig(job,ftype);
		    if(config != null) {
			if(editor != null)
			    editor.closeEditor();
			editor=new Editor(job);
		    }
                }
                else if(arg.equals("XYZ")) {
                    ftype="XYZ";
                    println("XYZ file input selected");
		    config=getConfig(job,ftype);
		    if(config != null) {
			if(editor != null)
			    editor.closeEditor();
			editor=new Editor(job);
		    }
                }
                else if(arg.equals("PDB")) {
                    ftype="PDB";
                    println("PDB file input selected");
		    config=getConfig(job,ftype);
		    if(config != null) {
			if(editor != null)
			    editor.closeEditor();
			editor=new Editor(job);
		    }
                }
                else if(arg.equals("MSI")) {
                    ftype="MSI";
                    println("MSI file input selected");
		    config=getConfig(job,ftype);
		    if(config != null) {
			if(editor != null)
			    editor.closeEditor();
			editor=new Editor(job);
		    }
                }
                else if(arg.equals("Lattice")){
                    if(maklat != null)
                        maklat.job.setVisible(false);
                    maklat=new MakeLattice(job);
                }
                else if(arg.equals("Bucky")){
                    if(makbuk != null)
                        makbuk.job.setVisible(false);
                    makbuk=new MakeBucky(job);
                }
                else if(arg.equals("Polymer")){
                    if(makpol != null)
                        makpol.job.setVisible(false);
                    makpol=new MakePoly(job);
                }
                else if(arg.equals("Chain")){
                    if(makchain != null)
                        makchain.job.setVisible(false);
                    makchain=new MakeChain(job);
                }
                else if(arg.equals("N_fold")){
                    if(enfold != null)
                        enfold.job.setVisible(false);
                    enfold=new Nfold(job);
                }
                else if(arg.equals("BondLengths")){
                    if(bondlen != null)
                        bondlen.job.setVisible(false);
                    bondlen=new BondLengths(job);
                }
                else if(arg.equals("Add Solvent")){
                    if(addh2o != null)
                        addh2o.job.setVisible(false);
                    addh2o=new SolventAdd(job);
                }
                else if(arg.equals("Slice")){
                    if(slcrev != null)
                        slcrev.job.setVisible(false);
                    slcrev=new Slice(job);
                }
                else if(arg.equals("Insert Molecule")){
                    if(insmol != null)
                        insmol.job.setVisible(false);
                    insmol=new InsertMolecule(job);
                }
                else if(arg.equals("Blank")){
                    if(makblank != null)
                        makblank.job.setVisible(false);
                    makblank=new MakeBlankField(job);
                }
                else if(arg.equals("Dreiding")){
                    if(makdrei != null)
                        makdrei.job.setVisible(false);
                    makdrei=new MakeDreiField(job);
                }
                else if(arg.equals("OPLS")){
                    if(makopls != null)
                        makopls.job.setVisible(false);
                    makopls=new MakeOPLSField(job);
                }
                else if(arg.equals("Ceramics")){
                    if(makceram != null)
                        makceram.job.setVisible(false);
                    makceram=new MakeCeramField(job);
                }
                else if(arg.equals("Table")){
                    if(maktable != null)
                        maktable.job.setVisible(false);
                    maktable=new MakeTable(job);
                }
                else if(arg.equals("Run DL_POLY")){
                    if(runjob != null)
                        runjob.job.setVisible(false);
                    runjob=new Execute(job);
                }
                else if(arg.equals("Store/Fetch Data")){
                    if(datarc != null)
                        datarc.job.setVisible(false);
                    datarc=new DataArchiver(job);
                }
                else if(arg.equals("RDF_Plot")){
                    if(rdfplt != null)
                        rdfplt.job.setVisible(false);
                    rdfplt=new RDFPlot(job);
                }
                else if(arg.equals("RDF_Calc")){
                    if(rdfcal != null)
                        rdfcal.job.setVisible(false);
                    rdfcal=new RDFCalc(job);
                }
                else if(arg.equals("S(k)")){
                    if(sokplt != null)
                        sokplt.job.setVisible(false);
                    sokplt=new SokPlot(job);
                }
                else if(arg.equals("Z_Density")){
                    if(zdnplt != null)
                        zdnplt.job.setVisible(false);
                    zdnplt=new ZdenPlot(job);
                }
                else if(arg.equals("MSD")){
                    if(msdrun != null)
                        msdrun.job.setVisible(false);
                    msdrun=new RunMSD(job);
                }
                else if(arg.equals("VAF")){
                    if(vafrun != null)
                        vafrun.job.setVisible(false);
                    vafrun=new RunVAF(job);
                }
                else if(arg.equals("FAF")){
                    if(fafrun != null)
                        fafrun.job.setVisible(false);
                    fafrun=new RunFAF(job);
                }
                else if(arg.equals("Gs(r,t)")){
                    if(gslcal != null)
                        gslcal.job.setVisible(false);
                    gslcal=new GslCalc(job);
                }
                else if(arg.equals("Gd(r,t)")){
                    if(gdfcal != null)
                        gdfcal.job.setVisible(false);
                    gdfcal=new GdfCalc(job);
                }
                else if(arg.equals("S(k,w)")) {
                    if(skwcal != null)
                        skwcal.job.setVisible(false);
                    skwcal=new SkwCalc(job);
                }
                else if(arg.equals("CONFIG")) {
                    config=getConfig(job,"CONFIG");
                }
                else if(arg.equals("REVCON")) {
                    config=getConfig(job,"REVCON");
                }
                else if(arg.equals("Plot")) {
                    if(graf != null)
                        graf.job.setVisible(false);
                    graf=new GraphDraw(job);
                }
                else if(arg.equals("What Atoms?")) {
                    println("Select required CONFIG file for input");
                    if((fname=selectFileNameContains(job,"CFG"))!=null) {
                        AML.whatAtoms(job,fname);
                    }
                    else {
                        println("File selection cancelled");
                    }
                }
                else if(arg.equals("Statistics")){
                    if(graf != null)
                        graf.job.setVisible(false);
                    graf=new GraphDraw(job);
                    if(staprp != null)
                        staprp.job.setVisible(false);
                    staprp=new StatProp(job);
                }
                else if(arg.equals("About DL_POLY")) {
                    viewResource("About_DL_POLY");
                }
                else if(arg.equals("Disclaimer")) {
                    viewResource("Disclaimer");
                }
                else if(arg.equals("Licence")) {
                    viewResource("Licence");
                }
                else if(arg.equals("Acknowledgements")) {
                    viewResource("Acknowledge");
                }
                else if(arg.equals("MINIDREI")) {
                    viewResource("MINIDREI");
                }
                else if(arg.equals("MINIOPLS")) {
                    viewResource("MINIOPLS");
                }
                else if(arg.equals("CERAMICS")) {
                    viewResource("CERAMICS");
                }
                else if(arg.equals("Test Cases")) {
                    viewResource("TestInfo");
                }
                else if(arg.equals("Clear Text")) {
                    clearScreen();
                }
            }
        }
    }

    // Main method

    public static void main(String args[]) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        // Select colour scheme

        String scheme="picasso";
        if(args.length>0)scheme=args[0].toLowerCase();
        art=new ColorScheme(scheme);

	// Set starting values

        startGUI();

        // Set up Graphical User interface

        job = new GUI();
	job.pack();
        job.setVisible(true);
        println(page);

    }

    static void startGUI() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */

        // Initialise GUI variables

        ftype="CFG";

        // File numbers

        numlat=0;
        numpol=0;
        numbuk=0;
        numblk=0;
        numchn=0;
        numctr=0;
        numtab=0;
        numxyz=0;
        nummsi=0;
        numpdb=0;
        numrdf=0;
        numsok=0;
        numsol=0;
        numzdn=0;
        nummsd=0;
        numvaf=0;
        numfaf=0;
        numstat=0;
        numgdf=0;
        numhov=0;
        numgsl=0;
        numskw=0;
        numsko=0;
        numh2o=0;
        numdre=0;
        numopl=0;
        numcer=0;
        numsav=0;
        numslc=0;

        // Define default bondlengths

        cc1b=1.54;
        cc2b=1.34;
        cc3b=1.20;
        ccab=1.39;
        ch1b=1.09;
        cn1b=1.47;
        cn2b=1.27;
        cn3b=1.16;
        cnab=1.35;
        co1b=1.43;
        co2b=1.22;
        coab=1.36;
        nh1b=1.01;
        oh1b=0.96;
        bondpc=BONDPERCENT;
        rotdef=DEFINE_ROTATION;
        tradef=DEFINE_TRANSLATION;
        incx=tradef;
        incy=tradef;
        incz=tradef;
        rotcos=Math.cos(Math.PI*rotdef/180.0);
        rotsin=Math.sin(Math.PI*rotdef/180.0);

        // nullify processes
        //
        //	proc=null;
        //	prdf=null;
        //	pgdf=null;
        //	pskw=null;
    }

    void draggedResponse(int a,int b){}

    JMenu MyMenu(String s) {
        JMenu  mine = new JMenu(s);
        mine.setFont(fontMain);
        mine.setForeground(art.fore);
        mine.setBackground(art.back);
        return mine;
    }

    JMenuItem MyMenuItem(String s) {
        JMenuItem  mine = new JMenuItem(s);
        mine.setFont(fontMain);
        mine.setForeground(art.fore);
        mine.setBackground(art.back);
        return mine;
    }

    public void printOut(){
        /*
*********************************************************************

dl_poly/java GUI routine to print text from the GUI

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        PrinterJob prnjob=PrinterJob.getPrinterJob();
        PageFormat pf=prnjob.defaultPage();
        pf.setOrientation(PageFormat.PORTRAIT);
        prnjob.setPrintable(job,pf);
        try{
            if(prnjob.printDialog()) {
                println("Initiating print .......");
                prnjob.print();
                println("Print complete");
            }
            else {
                println("Print cancelled");
            }
        }
        catch(PrinterException pe){
            System.out.println("Error - problem with print "+pe.toString());
        }

    }

    public int print(Graphics g,PageFormat pf,int pageIndex) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        if(pageIndex > 0) {
            return Printable.NO_SUCH_PAGE;
        }
        Graphics2D g2d=(Graphics2D)g;
        g2d.translate(pf.getImageableX(),pf.getImageableY());
        paint(g2d);
        return Printable.PAGE_EXISTS;
    }
}
