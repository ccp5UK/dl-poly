import java.awt.*;
import java.awt.event.*;
import javax.swing.*;

public class ChgDefaults extends Basic implements ActionListener {
    /*
**********************************************************************

dl_poly/java routine to change the defaults for the GUI

copyright - daresbury laboratory
author    - w. smith january 2001

**********************************************************************
     */
    public static GUI home;
    public static ChgDefaults job;
    private static JTextField rotang,tradis,bondtol;
    private static JButton set,close;

    public ChgDefaults() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        super();
        setTitle("Change Defaults");

        getContentPane().setForeground(art.fore);
        getContentPane().setBackground(art.back);
        setFont(fontMain);
        getContentPane().setLayout(new GridLayout(4,2));

        // define the buttons

        set=new JButton("Set");
        set.setForeground(art.butf);
        set.setBackground(art.butn);
        close=new JButton("Close");
        close.setForeground(art.butf);
        close.setBackground(art.butn);
        getContentPane().add(set);
        getContentPane().add(close);

        // Fixed rotation angle

        JLabel l07=new JLabel("Rotation angle (deg)");
        getContentPane().add(l07);
        rotang=new JTextField(8);
        rotang.setBackground(art.scrn);
        rotang.setForeground(art.scrf);
        getContentPane().add(rotang);

        // Fixed translation distance

        JLabel l08=new JLabel("Translation dist (A)");
        getContentPane().add(l08);
        tradis=new JTextField(8);
        tradis.setBackground(art.scrn);
        tradis.setForeground(art.scrf);
        getContentPane().add(tradis);

        // Bond acceptance tolerance (percent)

        JLabel l09=new JLabel("Bond Tolerance (%)");
        getContentPane().add(l09);
        bondtol=new JTextField(8);
        bondtol.setBackground(art.scrn);
        bondtol.setForeground(art.scrf);
        getContentPane().add(bondtol);

        // Register action buttons

        set.addActionListener(this);
        close.addActionListener(this);


    }

    public ChgDefaults(GUI here) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        home=here;
        println("Activated panel for changing defaults");
        job=new ChgDefaults();
        job.pack();
        job.setVisible(true);
        defaultSet();
    }

    public void defaultSet() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        rotang.setText(String.valueOf(rotdef));
        tradis.setText(String.valueOf(tradef));
        bondtol.setText(String.valueOf(bondpc));
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
        if (arg.equals("Set")) {
            rotdef=BML.giveDouble(rotang.getText(),1);
            tradef=BML.giveDouble(tradis.getText(),1);
            bondpc=BML.giveDouble(bondtol.getText(),1);
            rotcos=Math.cos(Math.PI*rotdef/180.0);
            rotsin=Math.sin(Math.PI*rotdef/180.0);
            incx=tradef;
            incy=tradef;
            incz=tradef;
            println("GUI defaults now changed");
        }
        else if (arg.equals("Close")) {
            job.setVisible(false);
        }
    }
}



