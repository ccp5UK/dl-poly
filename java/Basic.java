import java.awt.*;
import java.io.*;
import java.awt.event.*;
import java.io.File;
import javax.swing.*;
import javax.swing.filechooser.*;
import java.awt.print.*;
import java.awt.geom.*;
import java.util.*;

// Define the DL_POLY Graphical User Interface

abstract class Basic extends JFrame {
        /*
*********************************************************************

main dl_poly/java GUI super class

copyright - daresbury laboratory
author    - w.smith 2006

*********************************************************************
         */

    public static final int TILEX=500;
    public static final int TILEY=500;
    public static final int SCREENX=600;
    public static final int SCREENY=150;
    public static final int SCREENZ=525;
    public static final int FILENAME_HEAD=0;
    public static final int FILENAME_TAIL=1;
    public static final int FILENAME_CONTAINS=2;
    public static final int APPROVE_OPTION=0;
    public static final int CANCEL_OPTION=1;
    public static final double BONDPERCENT=10.0;
    public static final double DEFINE_ROTATION=30.0;
    public static final double DEFINE_TRANSLATION=1.0;
    public static final double RGAS=8.3145e-3;
    public static final int MXCONNECT=10;
    public static final int MXUNIQUE=100;
    public static final int MAXBAS=20;
    public static final int MXATMS=100;
    public static final int MXJOIN=100;
    public static String ftype=null;
    public static Font fontMain;
    public static boolean perspective=true;
    public static boolean showbonds=true;
    public static boolean showwater=false;
    public static boolean showcell=true;
    public static ColorScheme art;
    public static BondLengths bondlen=null;
    public static ChgDefaults newdefs=null;
    public static MakeControl makctr=null;
    public static MakeLattice maklat=null;
    public static MakeBucky makbuk=null;
    public static MakePoly makpol=null;
    public static MakeChain makchain=null;
    public static Nfold enfold=null;
    public static SolventAdd addh2o=null;
    public static Execute runjob=null;
    public static DataArchiver datarc=null;
    public static RDFPlot rdfplt=null;
    public static RDFCalc rdfcal=null;
    public static SokPlot sokplt=null;
    public static ZdenPlot zdnplt=null;
    public static Slice slcrev=null;
    public static InsertMolecule insmol=null;
    public static RunMSD msdrun=null;
    public static RunVAF vafrun=null;
    public static RunFAF fafrun=null;
    public static GslCalc gslcal=null;
    public static GdfCalc gdfcal=null;
    public static SkwCalc skwcal=null;
    public static StatProp staprp=null;
    public static MakeBlankField makblank=null;
    public static MakeDreiField makdrei=null;
    public static MakeOPLSField makopls=null;
    public static MakeTable maktable=null;
    public static MakeCeramField makceram=null;
    public static GraphDraw graf=null;
    public GridLayout grd;
    public static JTextArea board;
    public static Builder pane;
    public static JScrollPane scroll;
    public static Config config=null;
    public static Config cfgsav=null;
    public static Editor editor=null;
    public static boolean edit=false;
    public static boolean alert=false;
    public static int numbuk,numlat,numpol,numchn,numctr,numtab,numxyz;
    public static int nummsi,numpdb,numrdf,numsok,nummsd,numvaf,numfaf;
    public static int numstat,numzdn,numgdf,numhov,numgsl,numskw,numh2o;
    public static int numblk,numdre,numopl,numcer,numsav,numslc,numsko;
    public static int numsol;
    public static double cc1b,cc2b,cc3b,ccab,ch1b,cn1b,cn2b,cn3b,cnab;
    public static double coab,nh1b,oh1b,co1b,co2b,bondpc,tradef,rotdef;
    public static double rotcos,rotsin,incx,incy,incz;
    public static Process proc=null;
    public static String fname=null;
    public GuiFileFilter mf;
    public JFileChooser fc;

    void fix(Component cmp,GridBagLayout grd,GridBagConstraints gbc,
    int gx,int gy,int gw,int gh,int wx,int wy) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        gbc.gridx=gx;
        gbc.gridy=gy;
        gbc.gridwidth=gw;
        gbc.gridheight=gh;
        gbc.weightx=wx;
        gbc.weighty=wy;
        gbc.ipadx=2;
        gbc.ipady=2;
        grd.setConstraints(cmp,gbc);
        getContentPane().add(cmp);
    }

    void fix(Component cmp,GridBagLayout grd,GridBagConstraints gbc,int gx, int gy,
    int gw, int gh) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        gbc.gridx=gx;
        gbc.gridy=gy;
        gbc.gridwidth=gw;
        gbc.gridheight=gh;
        grd.setConstraints(cmp,gbc);
        gbc.ipadx=2;
        gbc.ipady=2;
        getContentPane().add(cmp);
    }

    boolean copyFile(String aname,String bname) {
        /*
*********************************************************************

dl_poly/java routine to copy formatted DL_POLY files

copyright - daresbury laboratory
author    - w.smith february 2001

*********************************************************************
         */
        String record="";

        try {
            LineNumberReader lnr = new LineNumberReader(new FileReader(aname));
            DataOutputStream out = new DataOutputStream(new FileOutputStream(bname));
            println("Copying file: "+aname+" to "+bname);
            record = lnr.readLine();
            out.writeBytes(record+"\n");
            println("File header record: "+record);
            while((record=lnr.readLine()) != null) {
                out.writeBytes(record+"\n");
            }
            lnr.close();
            out.close();
        }
        catch(FileNotFoundException e) {
            println("Error - file not found: " + " "+e);
            return false;
        }
        catch(Exception e) {
            println("Error copying file: " +" "+e);
            return false;
        }
        return true;
    }

    boolean copyData(String aname,String bname) {
        /*
*********************************************************************

dl_poly/java routine to copy unformatted DL_POLY files

copyright - daresbury laboratory
author    - w.smith february 2001

*********************************************************************
         */
        DataInputStream fdi=null;
        DataOutputStream fdo=null;
        try {
            fdi = new DataInputStream(new FileInputStream(aname));
            fdo = new DataOutputStream(new FileOutputStream(bname));
            println("Copying file: "+aname+" to "+bname);
            while(true) {
                fdo.writeByte(fdi.readByte());
            }
        }
        catch(FileNotFoundException e) {
            println("Error - file not found: " + " "+e);
            return false;
        }
        catch(IOException e) {
            println("File copy complete. Final report:" + " "+e);
        }
        finally {
            try{fdi.close();fdo.close();}catch(IOException e){}
        }
        return true;
    }

    int putXY(String fname,String header[],int nhead,int npnts,double xx[],double yy[]) {
        /*
*********************************************************************

dl_poly/java routine to write a simple XY file

copyright - daresbury laboratory
author    - w.smith february 2001

*********************************************************************
         */
        try {
            DataOutputStream dos = new DataOutputStream(new FileOutputStream(fname));
            for(int i=0;i<nhead;i++)
                dos.writeBytes("#"+header[i]+"\n");
            for(int i=0;i<npnts;i++) {
                dos.writeBytes(BML.fmt(xx[i],20)+BML.fmt(yy[i],20)+"\n");
            }
            dos.close();
        }
        catch(Exception e) {
            println("Error - writing file: "+fname);
            return -1;
        }
        return 0;
    }

    void clearScreen() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        board.setText("");
    }

    static void println(String s) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        try {
            board.append(s+"\n");
            board.scrollRectToVisible(new Rectangle(0,board.getHeight()-2,1,1));
        }
        catch(Exception e) {
            System.out.println(s);
        }
    }

    static void print(String s) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        try {
            board.append(s);
        }
        catch(Exception e) {
            System.out.println(s);
        }
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

    void getViewFile(GUI home) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */

        // View a text file

        JFileChooser fc = new JFileChooser(new File("./"));
        fc.setDialogTitle("Select file for viewing");
        fc.setForeground(art.fore);
        fc.setBackground(art.back);
        int key=fc.showDialog(home,"View");
        if(key == JFileChooser.APPROVE_OPTION) {
            fname=fc.getSelectedFile().getPath();
            viewFile(fname);
        }
        else if(key == JFileChooser.CANCEL_OPTION) {
            println("File selection cancelled");
        }
        else {
            println("Error - file viewer failure");
        }
    }

    static String selectFileNameBegins(GUI home,String filter) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2003

*********************************************************************
         */
        String filename=null;
        JFileChooser fc = new JFileChooser(new File("./"));
        fc.setDialogTitle("Select required file");
        GuiFileFilter mf = new GuiFileFilter(FILENAME_HEAD);
        mf.filter=filter;
        fc.addChoosableFileFilter(mf);
        int key=fc.showOpenDialog(home);
        if(key == JFileChooser.APPROVE_OPTION) {
            filename=fc.getSelectedFile().getPath();
        }
        else if(key == CANCEL_OPTION) {
            filename=null;
        }
        else {
            println("Error - file selection failure");
            filename=null;
        }
        return filename;
    }

    static String selectFileNameContains(GUI home,String filter) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2003

*********************************************************************
         */
        String filename=null;
        JFileChooser fc = new JFileChooser(new File("./"));
        fc.setDialogTitle("Select required file");
        GuiFileFilter mf = new GuiFileFilter(FILENAME_CONTAINS);
        mf.filter=filter;
        fc.addChoosableFileFilter(mf);
        int key=fc.showOpenDialog(home);
        if(key == JFileChooser.APPROVE_OPTION) {
            filename=fc.getSelectedFile().getPath();
        }
        else if(key == CANCEL_OPTION) {
            filename=null;
        }
        else {
            println("Error - file selection failure");
            filename=null;
        }
        return filename;
    }

    static String selectFileNameEnds(GUI home,String filter) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2003

*********************************************************************
         */
        String filename=null;
        JFileChooser fc = new JFileChooser(new File("./"));
        fc.setDialogTitle("Select required file");
        GuiFileFilter mf = new GuiFileFilter(FILENAME_TAIL);
        mf.filter=filter;
        fc.addChoosableFileFilter(mf);
        int key=fc.showOpenDialog(home);
        if(key == JFileChooser.APPROVE_OPTION) {
            filename=fc.getSelectedFile().getPath();
        }
        else if(key == CANCEL_OPTION) {
            filename=null;
        }
        else {
            println("Error - file selection failure");
            filename=null;
        }
        return filename;
    }

    void viewResource(String fileName) {
        /*
*********************************************************************

dl_poly/java routine to display a Jar resource file

copyright - daresbury laboratory
author    - w.smith may 2005

*********************************************************************
         */
        clearScreen();

        try {
            String text="";
            InputStream instream = this.getClass().getResourceAsStream(fileName);
            InputStreamReader isr = new InputStreamReader(instream);
            BufferedReader reader = new BufferedReader(isr);
            while((text=reader.readLine()) != null) {
                println(text);
            }
            reader.close();
        }
        catch(Exception e) {
            println("Error: reading file " + fileName);
        }
    }

    String dateToday() {
        /*
*********************************************************************

dl_poly/java routine to display current date as a text string

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        String today;
        int mnth;

        Calendar datum=Calendar.getInstance();

        mnth=datum.get(Calendar.MONTH)+1;
        today=datum.get(Calendar.YEAR)+"/"+mnth+
            "/"+datum.get(Calendar.DAY_OF_MONTH)+"/"
            +datum.get(Calendar.HOUR_OF_DAY)+":"
            +datum.get(Calendar.MINUTE);

        return today;
    }

    void viewFile(String fileName) {
        /*
*********************************************************************

dl_poly/java routine to display a text file

copyright - daresbury laboratory
author    - w.smith october 2000

*********************************************************************
         */
        String text = "";
        clearScreen();
        LineNumberReader info=null;
        try {
            info = new LineNumberReader(new FileReader(fileName));
            do {
                text = info.readLine();
                if (text != null)
                    println(text);
            }
            while (text != null);
            info.close();
        }
        catch(Exception e) {
            println("Error: reading file " + fileName);
        }
    }

    void zappFile(GUI home) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */

        // Delete a file selected from browser

        JFileChooser fc = new JFileChooser(new File("./"));
        fc.setDialogTitle("Select file for deletion");
        fc.setForeground(art.fore);
        fc.setBackground(art.back);
        int key=fc.showDialog(home,"Delete");
        if(key == JFileChooser.APPROVE_OPTION) {
            WarningBox danger=new WarningBox(home,"Warning!",true);
            danger.setVisible(true);
            if(alert) {
                fname=fc.getSelectedFile().getPath();
                (new File(fname)).delete();
                println("File "+fname+" has been deleted");
            }
            else {
                println("File deletion aborted");
            }
        }
        else if(key == JFileChooser.CANCEL_OPTION) {
            println("File selection cancelled");
        }
        else {
            println("Error - file deletion failure");
        }
    }

    LineNumberReader hread(String fname,String name[],LineNumberReader lnr,
    double info[],double cell[],double chge[],double weight[],
    double xyz[][],double vel[][],double frc[][]) {
        /*
***********************************************************************

dl_poly/java routine for reading the formatted history file

copyright - daresbury laboratory
author    - w. smith february 2001

Notes:
info[0] - 0 to open file, 1 to continue, -1 to close file
info[1] - mxatms - maximum number of atoms in configuration
info[2] - keytrj - trajectory key required by calling routine
info[3] - iflg - error flag
          0  - no error
         -1  - end of file reached
         -2  - File not found
         -3  - no velocity data in file
         -4  - no force data in file
         -5  - too many atoms in configuration
info[4] - ktrj - trajectory key of HISTORY file
info[5] - imcon - image convention (PBC of cell)
info[6] - nstep - current time step number
info[7] - natms - actual number of atoms in configuration
info[8] - tstep - time between HISTORY configurations (ps)
info[9] - nstep0 - first time step in HISTORY file

***********************************************************************
         */
        String record,step;
        double tstep;
        int begin,mxatms,keytrj,ktrj,imcon,nstep,natms;

        begin=BML.nint(info[0]);
        mxatms=BML.nint(info[1]);
        keytrj=BML.nint(info[2]);
        ktrj=BML.nint(info[4]);
        imcon=BML.nint(info[5]);

        // open the HISTORY file

        try {
            if(begin==0) {
                info[0]=1.0;
                lnr = new LineNumberReader(new FileReader(fname));
                println("Reading file: "+fname);
                record = lnr.readLine();
                println("HISTORY file header: "+record.trim());
                record=lnr.readLine();
                ktrj=BML.giveInteger(record,1);
                imcon=BML.giveInteger(record,2);
                natms=BML.giveInteger(record,3);
                info[4]=(double)ktrj;
                info[5]=(double)imcon;
                info[7]=(double)natms;
                if(keytrj>ktrj) {
                    if(ktrj==0) {
                        println("Error - no velocities in file");
                        lnr.close();
                        info[3]=-3.0;
                        return null;
                    }
                    if(keytrj>1) {
                        println("Error - no forces in file");
                        lnr.close();
                        info[3]=-4.0;
                        return null;
                    }
                }
                return lnr;
            }
            if(begin < 0) {
                System.out.println(".");
                lnr.close();
                return null;
            }

            // Configuration header

            record=lnr.readLine();
            if(record==null) {
                info[3]=-1.0;
                lnr.close();
                return null;
            }
            nstep=BML.giveInteger(record,2);
            natms=BML.giveInteger(record,3);
            tstep=BML.giveDouble(record,6);
            info[6]=(double)nstep;
            info[7]=(double)natms;
            info[8]=tstep;
            if(natms>mxatms) {
                println("Error - too many atoms in MD cell");
                println("File contains "+BML.fmt(natms,8)+" atoms");
                lnr.close();
                info[3]=-5.0;
                return null;
            }
            if(imcon>0) {
                record=lnr.readLine();
                cell[0]=BML.giveDouble(record,1);
                cell[1]=BML.giveDouble(record,2);
                cell[2]=BML.giveDouble(record,3);
                record=lnr.readLine();
                cell[3]=BML.giveDouble(record,1);
                cell[4]=BML.giveDouble(record,2);
                cell[5]=BML.giveDouble(record,3);
                record=lnr.readLine();
                cell[6]=BML.giveDouble(record,1);
                cell[7]=BML.giveDouble(record,2);
                cell[8]=BML.giveDouble(record,3);
            }
            System.out.print(".");
            for(int i=0;i<natms;i++) {
                record=lnr.readLine();
                name[i]=BML.giveWord(record,1);
                weight[i]=BML.giveDouble(record,3);
                chge[i]=BML.giveDouble(record,4);
                record=lnr.readLine();
                xyz[0][i]=BML.giveDouble(record,1);
                xyz[1][i]=BML.giveDouble(record,2);
                xyz[2][i]=BML.giveDouble(record,3);
                if(keytrj>0) {
                    record=lnr.readLine();
                    vel[0][i]=BML.giveDouble(record,1);
                    vel[1][i]=BML.giveDouble(record,2);
                    vel[2][i]=BML.giveDouble(record,3);
                }
                else if(ktrj>0) {
                    record=lnr.readLine();
                }
                if(keytrj>1) {
                    record=lnr.readLine();
                    frc[0][i]=BML.giveDouble(record,1);
                    frc[1][i]=BML.giveDouble(record,2);
                    frc[2][i]=BML.giveDouble(record,3);
                }
                else if(ktrj>1) {
                    record=lnr.readLine();
                }
            }
        }
        catch(FileNotFoundException e) {
            println("Error - file not found: " + fname);
            info[3]=-2.0;
            return lnr;
        }
        catch(Exception e) {
            println("Error reading file: " + fname + " "+e);
            try{lnr.close();}catch(IOException t){}
            return lnr;
        }
        return lnr;
    }

    static void fft(boolean safe,int ind,int isw,int ndiv,int key[],
             double aaa[][],double wfft[][],double bbb[][]){
        /*
***********************************************************************

DL_POLY java fast fourier transform routine (for radix 2 only)

copyright - daresbury laboratory
author    - w smith  Nov 2007

Note the following (complex) array dimensions:

double aaa[][]=new double[2][ndiv];
double bbb[][]=new double[2][ndiv];
double wfft[][]=new double[2][ndiv];

***********************************************************************
*/
        boolean check;
        int nt,nu,iii,jjj,kkk,jj2,np2,k12,nu1,kd2;
        double tpn,arg,tt0,tt1;

        safe=true;

        //check that radix 2 rule applies

        nu=0;
        nt=1;
        check=true;
        for(int i=1;i<21;i++){
            nt=2*nt;
            if(nt==ndiv){
                check=false;
                nu=i;
            }
        }

        if(check){
            println("Error - number of points not a power of two");
            safe=false;
            return;
        }

        if(ind > 0){

            //set reverse bit address array

            for(kkk=0;kkk<ndiv;kkk++){
                iii=0;
                jjj=kkk;
                for(int j=1;j<=nu;j++){
                    jj2=jjj/2;
                    iii=2*(iii-jj2)+jjj;
                    jjj=jj2;
                }
                key[kkk]=iii;
            }

            //initialise complex exponential factors

            arg=0;
            np2=ndiv/2;
            tpn=2*Math.PI/ndiv;
            wfft[0][0]=1;
            wfft[1][0]=0;
            for(int i=1;i<=np2;i++){
                arg=tpn*i;
                wfft[0][i]=Math.cos(arg);
                wfft[1][i]=Math.sin(arg);
                wfft[0][ndiv-i]=wfft[0][i];
                wfft[1][ndiv-i]=wfft[1][i];
            }
            return;
        }

        //take conjugate of exponentials if required

        if(isw < 0){
            for(int i=1;i<ndiv;i++)
                wfft[1][i]=-wfft[1][i];
        }

        //make copy of input array

        for(int i=0;i<ndiv;i++){
            bbb[0][i]=aaa[0][i];
            bbb[1][i]=aaa[1][i];
        }

        //perform fourier transform

        kkk=0;
        nu1=nu-1;
        np2=ndiv/2;
        for(int l=0;l<nu;l++){
            kd2=BML.nint(Math.pow(2,nu1));
            while(kkk < ndiv){
                for(int i=0;i<np2;i++){
                    iii=key[kkk/kd2];
                    k12=kkk+np2;
                    tt0=bbb[0][k12]*wfft[0][iii]-bbb[1][k12]*wfft[1][iii];
                    tt1=bbb[0][k12]*wfft[1][iii]+bbb[1][k12]*wfft[0][iii];
                    bbb[0][k12]=bbb[0][kkk]-tt0;
                    bbb[1][k12]=bbb[1][kkk]-tt1;
                    bbb[0][kkk]=bbb[0][kkk]+tt0;
                    bbb[1][kkk]=bbb[1][kkk]+tt1;
                    kkk++;
                }
                kkk+=np2;
            }
            kkk=0;
            np2/=2;
            nu1--;
        }

        //unscramble the fft using bit address array

        for(kkk=0;kkk<ndiv;kkk++){
            iii=key[kkk];
            if(iii > kkk){
                tt0=bbb[0][kkk];
                tt1=bbb[1][kkk];
                bbb[0][kkk]=bbb[0][iii];
                bbb[1][kkk]=bbb[1][iii];
                bbb[0][iii]=tt0;
                bbb[1][iii]=tt1;
            }
        }

        //restore exponentials to unconjugated values if necessary

        if(isw < 0){
            for(int i=1;i<ndiv;i++)
                wfft[1][i]=-wfft[1][i];
        }

        return;
    }

    static Config getConfig(GUI here,String ftype) {
        /*
*********************************************************************

dl_poly/java method for reading CONFIG files

copyright - daresbury laboratory
author    - w.smith june 2001

*********************************************************************
         */

        Config newcfg=new Config();
        newcfg.home=here;

        if(ftype.equals("CFG")) {
            println("Select required CONFIG file for input");
            if((newcfg.fname=selectFileNameBegins(here,"CFG"))!=null) {
                if(!newcfg.rdCFG(newcfg.fname)) {
                    println("Problem reading CONFIG file");
                    newcfg=null;
                }
            }
            else {
                println("File selection cancelled");
                newcfg=null;
            }
        }
        else if(ftype.equals("XYZ")) {
            println("Select required XYZ file for input");
            if((newcfg.fname=selectFileNameEnds(here,"XYZ"))!=null) {
                if(newcfg.rdXYZ(newcfg.fname)) {
                    newcfg.fname="CFGXYZ."+numxyz;
                    if(newcfg.configWrite(newcfg.fname)) numxyz++;
                }
                else {
                    println("Problem reading XYZ file");
                    newcfg=null;
                }
            }
            else {
                println("File selection cancelled");
                newcfg=null;
            }
        }
        else if(ftype.equals("PDB")) {
            println("Select required PDB file for input");
            if((newcfg.fname=selectFileNameEnds(here,"PDB"))!=null) {
                if(newcfg.rdPDB(newcfg.fname)) {
                    newcfg.fname="CFGPDB."+numpdb;
                    if(newcfg.configWrite(newcfg.fname)) numpdb++;
                }
                else {
                    println("Problem reading PDB file");
                    newcfg=null;
                }
            }
            else {
                println("File selection cancelled");
                newcfg=null;
            }
        }
        else if(ftype.equals("MSI")) {
            println("Select required CERIUS/MSI file for input");
            if((newcfg.fname=selectFileNameEnds(here,"MSI"))!=null) {
                if(newcfg.rdMSI(newcfg.fname)) {
                    newcfg.fname="CFGMSI."+nummsi;
                    if(newcfg.configWrite(newcfg.fname)) nummsi++;
                }
                else {
                    println("Problem reading MSI file");
                    newcfg=null;
                }
            }
            else {
                println("File selection cancelled");
                newcfg=null;
            }
        }
        else if(ftype.equals("CONFIG")) {
            newcfg.fname="CONFIG";
            if(!newcfg.rdCFG(newcfg.fname)) {
                println("Problem reading CONFIG file");
                newcfg=null;
            }
        }
        else if(ftype.equals("REVCON")) {
            newcfg.fname="REVCON";
            if(!newcfg.rdCFG(newcfg.fname)) {
                println("Problem reading REVCON file");
                newcfg=null;
            }
        }
        return newcfg;
    }

    static Config copyConfig(Config cfg) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        Config cfgcpy=new Config();
        cfgcpy.home=cfg.home;
        cfgcpy.fname=cfg.fname;
        cfgcpy.title=cfg.title;
        cfgcpy.natms=cfg.natms;
        cfgcpy.pbc.imcon=cfg.pbc.imcon;
        cfgcpy.atoms=new Element[cfg.natms];
        cfgcpy.xyz=new double[3][cfg.natms];
        for (int i=0;i<9;i++)
            cfgcpy.pbc.cell[i]=cfg.pbc.cell[i];
        cfgcpy.pbc.buildBoundary(cfg.pbc.imcon);
        for(int i=0;i<cfg.natms;i++){
            cfgcpy.atoms[i]=new Element(cfg.atoms[i].zsym);
            cfgcpy.xyz[0][i]=cfg.xyz[0][i];
            cfgcpy.xyz[1][i]=cfg.xyz[1][i];
            cfgcpy.xyz[2][i]=cfg.xyz[2][i];
        }
        cfgcpy.structure=new Structure(cfgcpy);

        return cfgcpy;
    }

}
