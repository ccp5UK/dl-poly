import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;

public class DataArchiver extends Basic implements ActionListener {
        /*
*********************************************************************

dl_poly/java GUI class for data archive methods

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
    private static GUI home;
    public static DataArchiver job;
    private static JButton select,store,fetch,info,close;
    private static JTextField dirold,dirnew;
    private static JComboBox test;
    private static String dname="DEFAULT";
    private static WarningBox danger;

    // Define the Graphical User Interface

    public DataArchiver() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        super();
        setTitle("Data Archive");

        getContentPane().setBackground(art.back);
        getContentPane().setForeground(art.fore);
	setDefaultCloseOperation(DISPOSE_ON_CLOSE);
        setFont(fontMain);
        GridBagLayout grd = new GridBagLayout();
        GridBagConstraints gbc = new GridBagConstraints();
        getContentPane().setLayout(grd);

        gbc.fill=GridBagConstraints.BOTH;

        // Input file selection

        JLabel lab1=new JLabel("Select standard test case.....");
        fix(lab1,grd,gbc,0,0,3,1);

        // Select button

        select = new JButton("Select");
        select.setBackground(art.butn);
        select.setForeground(art.butf);
        fix(select,grd,gbc,0,1,1,1);

        // Test case choice

        test = new JComboBox();
        test.setBackground(art.scrn);
        test.setForeground(art.scrf);
        test.addItem("TEST1");
        test.addItem("TEST2");
        test.addItem("TEST3");
        test.addItem("TEST4");
        test.addItem("TEST5");
        test.addItem("TEST6");
        test.addItem("TEST7");
        test.addItem("TEST8");
        test.addItem("TEST9");
        test.addItem("TEST10");
        test.addItem("TEST11");
        test.addItem("TEST12");
        test.addItem("TEST13");
        test.addItem("TEST14");
        test.addItem("TEST15");
        test.addItem("TEST16");
        test.addItem("TEST17");
        test.addItem("TEST18");
        test.addItem("TEST19");
        test.addItem("TEST20");
        fix(test,grd,gbc,2,1,1,1);

        // Copy files from archive

        JLabel lab2=new JLabel("Data retrieval:");
        fix(lab2,grd,gbc,0,2,1,1);
        JLabel lab4=new JLabel("Directory:    ",JLabel.RIGHT);
        fix(lab4,grd,gbc,2,2,1,1);

        // Fetch files from archive

        fetch = new JButton("Fetch");
        fetch.setBackground(art.butn);
        fetch.setForeground(art.butf);
        fix(fetch,grd,gbc,0,3,1,1);

        dirold = new JTextField(dname);
        dirold.setBackground(art.scrn);
        dirold.setForeground(art.scrf);
        fix(dirold,grd,gbc,2,3,1,1);

        // Move files to archive

        JLabel lab3=new JLabel("Data storage:");
        fix(lab3,grd,gbc,0,4,1,1);
        JLabel lab5=new JLabel("Directory:    ",JLabel.RIGHT);
        fix(lab5,grd,gbc,2,4,1,1);

        store = new JButton("Store");
        store.setBackground(art.butn);
        store.setForeground(art.butf);
        fix(store,grd,gbc,0,5,1,1);

        dirnew = new JTextField(dname);
        dirnew.setBackground(art.scrn);
        dirnew.setForeground(art.scrf);
        fix(dirnew,grd,gbc,2,5,1,1);

        // pad out

        JLabel lab6=new JLabel("    ");
        fix(lab6,grd,gbc,1,6,1,1);

        // Information button

        info = new JButton("Info");
        info.setBackground(art.butn);
        info.setForeground(art.butf);
        fix(info,grd,gbc,0,7,1,1);

        // Close button

        close = new JButton("Close");
        close.setBackground(art.butn);
        close.setForeground(art.butf);
        fix(close,grd,gbc,2,7,1,1);

        // Register action buttons

        select.addActionListener(this);
        info.addActionListener(this);
        fetch.addActionListener(this);
        store.addActionListener(this);
        close.addActionListener(this);

    }

    public DataArchiver(GUI here) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        home=here;
        println("Activated DL_POLY DataArchiver");
        job=new DataArchiver();
        job.pack();
        job.setVisible(true);
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
        if (arg.equals("Select")) {
            dname=test.getSelectedItem().toString()+"/VV";
            println("About to overwrite current DL_POLY I/O files");
            danger=new WarningBox(home,"Warning!",true);
            danger.setVisible(true);
            if(alert)
                fetchFiles();
            else
                println("Operation cancelled");
        }
        else if (arg.equals("Fetch")) {
            dname=dirold.getText();
            println("About to overwrite current DL_POLY I/O files");
            danger=new WarningBox(home,"Warning!",true);
            danger.setVisible(true);
            if(alert)
                fetchFiles();
            else
                println("Operation cancelled");
        }
        else if (arg.equals("Store")) {
            dname=dirnew.getText();
            storeFiles();
        }
        else if (arg.equals("Info")) {
            viewResource("TestInfo");
        }
        else if (arg.equals("Close")) {
	    job.dispose();
        }
    }
    void fetchFiles() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        File f=null;
        String fname="";
        println("Fetching DL_POLY files from directory "+dname+".....");

        if((f=new File(fname="../data/"+dname.trim()+"/CONTROL")).exists()) {
            if(copyFile(fname,"CONTROL"))
                println("File CONTROL copied");
        }
        else {
            println("CONTROL file not found");
        }
        if((f=new File(fname="../data/"+dname.trim()+"/CONFIG")).exists()) {
            if(copyFile(fname,"CONFIG"))
                println("File CONFIG copied");
        }
        else {
            println("CONFIG file not found");
        }
        if((f=new File(fname="../data/"+dname.trim()+"/FIELD")).exists()) {
            if(copyFile(fname,"FIELD"))
                println("File FIELD copied");
        }
        else {
            println("FIELD file not found");
        }
        if((f=new File(fname="../data/"+dname.trim()+"/TABLE")).exists()) {
            if(copyFile(fname,"TABLE"))
                println("File TABLE copied");
        }
        else {
            println("TABLE file not found");
        }
        if((f=new File(fname="../data/"+dname.trim()+"/REVIVE")).exists()) {
            if(copyFile(fname,"REVIVE"))
                println("File REVIVE copied");
        }
        else {
            println("REVIVE file not found");
        }
        if((f=new File(fname="../data/"+dname.trim()+"/REVCON")).exists()) {
            if(copyFile(fname,"REVCON"))
                println("File REVCON copied");
        }
        else {
            println("REVCON file not found");
        }
        if((f=new File(fname="../data/"+dname.trim()+"/OUTPUT")).exists()) {
            if(copyFile(fname,"OUTPUT"))
                println("File OUTPUT copied");
        }
        else {
            println("OUTPUT file not found");
        }
        if((f=new File(fname="../data/"+dname.trim()+"/STATIS")).exists()) {
            if(copyFile(fname,"STATIS"))
                println("File STATIS copied");
        }
        else {
            println("STATIS file not found");
        }
        if((f=new File(fname="../data/"+dname.trim()+"/HISTORY")).exists()) {
            if(copyFile(fname,"HISTORY"))
                println("File HISTORY copied");
        }
        else {
            println("HISTORY file not found");
        }
        if((f=new File(fname="../data/"+dname.trim()+"/REVOLD")).exists()) {
            if(copyFile(fname,"REVOLD"))
                println("File REVOLD copied");
        }
        else {
            println("REVOLD file not found");
        }
        println("DL_POLY file copy completed");
    }
    void storeFiles() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        File f=null;
        String fname="";
        println("Storing DL_POLY files in directory "+dname+".....");

        if(!(f=new File("../data/"+dname.trim())).exists()) {

            (new File("../data/"+dname.trim())).mkdir();

            if((f=new File("CONTROL")).exists()) {
                fname="../data/"+dname.trim()+"/CONTROL";
                copyFile("CONTROL",fname);
                println("File CONTROL stored");
                f.delete();
            }
            if((f=new File("CONFIG")).exists()) {
                fname="../data/"+dname.trim()+"/CONFIG";
                copyFile("CONFIG",fname);
                println("File CONFIG stored");
                f.delete();
            }
            if((f=new File("FIELD")).exists()) {
                fname="../data/"+dname.trim()+"/FIELD";
                copyFile("FIELD",fname);
                println("File FIELD stored");
                f.delete();
            }
            if((f=new File("TABLE")).exists()) {
                fname="../data/"+dname.trim()+"/TABLE";
                copyFile("TABLE",fname);
                println("File TABLE stored");
                f.delete();
            }
            if((f=new File("OUTPUT")).exists()) {
                fname="../data/"+dname.trim()+"/OUTPUT";
                copyFile("OUTPUT",fname);
                println("File OUTPUT stored");
                f.delete();
            }
            if((f=new File("REVIVE")).exists()) {
                fname="../data/"+dname.trim()+"/REVIVE";
                copyFile("REVIVE",fname);
                println("File REVIVE stored");
                f.delete();
            }
            if((f=new File("REVOLD")).exists()) {
                fname="../data/"+dname.trim()+"/REVOLD";
                copyFile("REVOLD",fname);
                println("File REVOLD stored");
                f.delete();
            }
            if((f=new File("REVCON")).exists()) {
                fname="../data/"+dname.trim()+"/REVCON";
                copyFile("REVCON",fname);
                println("File REVCON stored");
                f.delete();
            }
            if((f=new File("STATIS")).exists()) {
                fname="../data/"+dname.trim()+"/STATIS";
                copyFile("STATIS",fname);
                println("File STATIS stored");
                f.delete();
            }
            if((f=new File("HISTORY")).exists()) {
                fname="../data/"+dname.trim()+"/HISTORY";
                copyFile("HISTORY",fname);
                println("File HISTORY stored");
                f.delete();
            }
            println("DL_POLY files stored");
        }
        else {
            println("Error - nominated storage directory already exists");
        }
    }
}
