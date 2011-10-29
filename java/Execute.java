import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;

public class Execute extends Basic implements ActionListener {
        /*
*********************************************************************

dl_poly/java GUI class to control DL_POLY execution

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
    private static GUI home;
    private static double time0,time1;
    public static Execute job;
    private static JTextField task;
    private static JButton exec,stat,zapp,updt,dlte,ctrl,fild,cnfg,tabl,close;
    private static WarningBox danger=null;

    // Define the Graphical User Interface

    public Execute() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        super();
        setTitle("Run DL_POLY");

        getContentPane().setBackground(art.back);
        getContentPane().setForeground(art.fore);
        setDefaultCloseOperation(DISPOSE_ON_CLOSE);
        setFont(fontMain);
        GridBagLayout grd = new GridBagLayout();
        GridBagConstraints gbc = new GridBagConstraints();
        getContentPane().setLayout(grd);

        gbc.fill=GridBagConstraints.BOTH;

        int n=0;

        // Run button

        exec = new JButton("Run");
        exec.setBackground(art.butn);
        exec.setForeground(art.butf);
        fix(exec,grd,gbc,0,n,1,1);

        // Define the Close button

        close = new JButton("Close");
        close.setBackground(art.butn);
        close.setForeground(art.butf);
        fix(close,grd,gbc,1,n++,1,1);

        // Describe the run command

        fix(new JLabel("Current run command:"),grd,gbc,0,n++,2,1);
        task=new JTextField();
        task.setBackground(art.scrn);
        task.setForeground(art.scrf);
        fix(task,grd,gbc,0,n,2,1);

        // Input file selection

        fix(new JLabel("Select required input files:"),grd,gbc,0,n++,2,1);

        // Select the CONTROL file

        ctrl = new JButton("CONTROL");
        ctrl.setBackground(art.butn);
        ctrl.setForeground(art.butf);
        fix(ctrl,grd,gbc,0,n,1,1);

        // Select the CONFIG file

        cnfg = new JButton("CONFIG");
        cnfg.setBackground(art.butn);
        cnfg.setForeground(art.butf);
        fix(cnfg,grd,gbc,1,n++,1,1);

        // Select the FIELD file

        fild = new JButton("FIELD");
        fild.setBackground(art.butn);
        fild.setForeground(art.butf);
        fix(fild,grd,gbc,0,n,1,1);

        // Select the TABLE file

        tabl = new JButton("TABLE");
        tabl.setBackground(art.butn);
        tabl.setForeground(art.butf);
        fix(tabl,grd,gbc,1,n++,1,1);

        // Job monitoring options

        fix(new JLabel("Job monitoring options:"),grd,gbc,0,n++,2,1);

        // Kill job

        zapp = new JButton("Kill");
        zapp.setBackground(art.butn);
        zapp.setForeground(art.butf);
        fix(zapp,grd,gbc,0,n,1,1);

        // Job status

        stat = new JButton("Status");
        stat.setBackground(art.butn);
        stat.setForeground(art.butf);
        fix(stat,grd,gbc,1,n++,1,1);

        // File handling options

        fix(new JLabel("File handling options:"),grd,gbc,0,n++,2,1);

        // Clear data files

        dlte = new JButton("Clear");
        dlte.setBackground(art.butn);
        dlte.setForeground(art.butf);
        fix(dlte,grd,gbc,0,n,1,1);

        // Update data files

        updt = new JButton("Update");
        updt.setBackground(art.butn);
        updt.setForeground(art.butf);
        fix(updt,grd,gbc,1,n++,1,1);

        // Register action buttons

        exec.addActionListener(this);
        ctrl.addActionListener(this);
        cnfg.addActionListener(this);
        fild.addActionListener(this);
        tabl.addActionListener(this);
        updt.addActionListener(this);
        zapp.addActionListener(this);
        dlte.addActionListener(this);
        stat.addActionListener(this);
        close.addActionListener(this);

    }

    public Execute(GUI here) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        home=here;
        println("Activated DL_POLY Execute panel");
        job=new Execute();
        job.pack();
        if(executable==null)
            executable=new String(defaultexecutable);
        task.setText(executable);
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
        String fname;
        String arg = (String)e.getActionCommand();
        if (arg.equals("Run")) {
            try {
                if(proc == null) {
                    executable=task.getText();
                    if(!(new File(executable)).exists()) {
                        println("Error - "+executable+" program not available");
                    }
                    else {
                        time0=(double)System.currentTimeMillis();
                        proc=Runtime.getRuntime().exec(executable+"&");
                        println(executable+" job submitted: "+String.valueOf(proc));
                    }
                }
                else {
                    println("Error - "+executable+" job already running");
                }
            }
            catch(IOException ee) {
                println("Error - "+executable+" job submission failure");
            }
        }
        else if (arg.equals("Close")) {
            job.dispose();
        }
        else if (arg.equals("CONTROL")) {
            println("Select required CONTROL file");
            if((fname=selectFileNameBegins(home,"CNT"))!=null) {
                if(copyFile(fname,"CONTROL"))
                    println("Control file selected");
            }
            else
                println("No file selected");
        }
        else if (arg.equals("CONFIG")) {
            println("Select required CONFIG file");
            if((fname=selectFileNameBegins(home,"CFG"))!=null) {
                if(copyFile(fname,"CONFIG"))
                    println("CONFIG file selected");
            }
            else
                println("No file selected");
        }
        else if (arg.equals("FIELD")) {
            println("Select required FIELD file");
            if((fname=selectFileNameBegins(home,"FLD"))!=null) {
                if(copyFile(fname,"FIELD"))
                    println("FIELD file selected");
            }
            else
                println("No file selected");
        }
        else if (arg.equals("TABLE")) {
            println("Select required TABLE file");
            if((fname=selectFileNameBegins(home,"TAB"))!=null) {
                if(copyFile(fname,"TABLE"))
                    println("TABLE file selected");
            }
            else
                println("No file selected");
        }
        else if (arg.equals("Kill") && proc != null) {
            println("Cancelling "+executable+" job: "+String.valueOf(proc));
            proc.destroy();
            proc=null;
            println(executable+" job cancelled");
        }
        else if (arg.equals("Status") && proc != null) {
            try {
                int state=proc.exitValue();
                if(state==0)
                    println(executable+" has terminated normally");
                else
                    println(executable+" has terminated abnormally");
                proc=null;
            }
            catch(IllegalThreadStateException ee) {
                time1=((double)System.currentTimeMillis()-time0)/1000.;

                println(executable+" job is running. Elapsed time (s)= "+BML.fmt(time1,9));
            }
        }
        else if (arg.equals("Clear")) {
            println("About to delete all current DL_POLY I/O files");
            danger=new WarningBox(home,"Warning!",true);
            danger.setVisible(true);
            if(alert)
                wipeOutFiles();
            else
                println("Operation cancelled");
        }
        else if (arg.equals("Update")) {
            println("About to overwrite some current DL_POLY I/O files");
            danger=new WarningBox(home,"Warning!",true);
            danger.setVisible(true);
            if(alert)
                updateFiles();
            else
                println("Operation cancelled");
        }
    }
    void wipeOutFiles() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        File f=null;
        println("Deleting all current DL_POLY I/O files .....");
        if((f=new File("CONTROL")).exists()) {
            f.delete();
            println("File CONTROL deleted");
        }
        if((f=new File("CONFIG")).exists()) {
            f.delete();
            println("File CONFIG deleted");
        }
        if((f=new File("FIELD")).exists()) {
            f.delete();
            println("File FIELD deleted");
        }
        if((f=new File("TABLE")).exists()) {
            f.delete();
            println("File TABLE deleted");
        }
        if((f=new File("OUTPUT")).exists()) {
            f.delete();
            println("File OUTPUT deleted");
        }
        if((f=new File("REVCON")).exists()) {
            f.delete();
            println("File REVCON deleted");
        }
        if((f=new File("REVIVE")).exists()) {
            f.delete();
            println("File REVIVE deleted");
        }
        if((f=new File("REVOLD")).exists()) {
            f.delete();
            println("File REVOLD deleted");
        }
        if((f=new File("STATIS")).exists()) {
            f.delete();
            println("File STATIS deleted");
        }
        if((f=new File("HISTORY")).exists()) {
            f.delete();
            println("File HISTORY deleted");
        }
        println("All current DL_POLY I/O files deleted");
    }
    void updateFiles() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        File f=null;
        println("Updating DL_POLY I/O files .....");
        if((f=new File("CONFIG")).exists())
            if(copyFile("CONFIG","CONFIG.BAK"))
                println("CONFIG file backup taken: CONFIG.BAK");
            else
                println("Error CONFIG file not found");
        if((f=new File("REVCON")).exists())
            if(copyFile("REVCON","CONFIG"))
                println("REVCON file renamed as CONFIG file");
            else
                println("Error REVCON file not found");
        if((f=new File("REVOLD")).exists())
            if(copyData("REVOLD","REVOLD.BAK"))
                println("REVOLD file backup taken: REVOLD.BAK");
            else
                println("Error REVOLD file not found");
        if((f=new File("REVIVE")).exists())
            if(copyData("REVIVE","REVOLD"))
                println("REVIVE file renamed as REVOLD file");
            else
                println("Error REVIVE file not found");
        println("DL_POLY I/O files updated.");
    }
}
