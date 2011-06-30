import java.awt.*;
import java.awt.event.*;
import javax.swing.*;

public class BondLengths extends Basic implements ActionListener {
            /*
**********************************************************************

dl_poly/java class to define the default bondlengths for the GUI

copyright - daresbury laboratory
author    - w. smith november 2000

**********************************************************************
             */
    static GUI home;
    static BondLengths job;
    static JTextField b01,b02,b03,b04,b05,b06,b07,b08,b09,b10,b11,b12,b13,b14;
    static JLabel l01,l02,l03,l04,l05,l06,l07,l08,l09,l10,l11,l12,l13,l14;
    static JButton set,close;
    static Font fontMain;

    public BondLengths() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */

        super();
	setTitle("Set Bond Lengths");
        getContentPane().setBackground(art.back);
        getContentPane().setForeground(art.fore);
        setFont(fontMain);
        getContentPane().setLayout(new GridLayout(15,2));

        // define the buttons

        set=new JButton("Set");
        set.setForeground(art.butf);
        set.setBackground(art.butn);
        close=new JButton("Close");
        close.setForeground(art.butf);
        close.setBackground(art.butn);
        getContentPane().add(set);
        getContentPane().add(close);

        // C-C single bondlength

        l01=new JLabel("C-C 1-Bond");
        getContentPane().add(l01);
        b01=new JTextField(8);
        b01.setBackground(art.scrn);
        b01.setForeground(art.scrf);
        getContentPane().add(b01);

        // C-C double bondlength

        l02=new JLabel("C-C 2-Bond");
        getContentPane().add(l02);
        b02=new JTextField(8);
        b02.setBackground(art.scrn);
        b02.setForeground(art.scrf);
        getContentPane().add(b02);

        // C-C triple bondlength

        l03=new JLabel("C-C 3-Bond");
        getContentPane().add(l03);
        b03=new JTextField(8);
        b03.setBackground(art.scrn);
        b03.setForeground(art.scrf);
        getContentPane().add(b03);

        // C-C aromatic bondlength

        l04=new JLabel("C-C R-Bond");
        getContentPane().add(l04);
        b04=new JTextField(8);
        b04.setBackground(art.scrn);
        b04.setForeground(art.scrf);
        getContentPane().add(b04);

        // C-H bondlength

        l05=new JLabel("C-H Bond");
        getContentPane().add(l05);
        b05=new JTextField(8);
        b05.setBackground(art.scrn);
        b05.setForeground(art.scrf);
        getContentPane().add(b05);

        // C-N single bondlength

        l06=new JLabel("C-N 1-Bond");
        getContentPane().add(l06);
        b06=new JTextField(8);
        b06.setBackground(art.scrn);
        b06.setForeground(art.scrf);
        getContentPane().add(b06);

        // C-N double bondlength

        l07=new JLabel("C-N 2-Bond");
        getContentPane().add(l07);
        b07=new JTextField(8);
        b07.setBackground(art.scrn);
        b07.setForeground(art.scrf);
        getContentPane().add(b07);

        // C-N triple bondlength

        l08=new JLabel("C-N 3-Bond");
        getContentPane().add(l08);
        b08=new JTextField(8);
        b08.setBackground(art.scrn);
        b08.setForeground(art.scrf);
        getContentPane().add(b08);

        // C-N aromatic bondlength

        l09=new JLabel("C-N R-Bond");
        getContentPane().add(l09);
        b09=new JTextField(8);
        b09.setBackground(art.scrn);
        b09.setForeground(art.scrf);
        getContentPane().add(b09);

        // C-O single bondlength

        l10=new JLabel("C-O 1-Bond");
        getContentPane().add(l10);
        b10=new JTextField(8);
        b10.setBackground(art.scrn);
        b10.setForeground(art.scrf);
        getContentPane().add(b10);

        // C-O double bondlength

        l11=new JLabel("C-O 2-Bond");
        getContentPane().add(l11);
        b11=new JTextField(8);
        b11.setBackground(art.scrn);
        b11.setForeground(art.scrf);
        getContentPane().add(b11);

        // C-O aromatic bondlength

        l12=new JLabel("C-O R-Bond");
        getContentPane().add(l12);
        b12=new JTextField(8);
        b12.setBackground(art.scrn);
        b12.setForeground(art.scrf);
        getContentPane().add(b12);

        // N-H bondlength

        l13=new JLabel("N-H Bond");
        getContentPane().add(l13);
        b13=new JTextField(8);
        b13.setBackground(art.scrn);
        b13.setForeground(art.scrf);
        getContentPane().add(b13);

        // O-H bondlength

        l14=new JLabel("O-H Bond");
        getContentPane().add(l14);
        b14=new JTextField(8);
        b14.setBackground(art.scrn);
        b14.setForeground(art.scrf);
        getContentPane().add(b14);

        // Register action buttons

        set.addActionListener(this);
        close.addActionListener(this);

    }
    public BondLengths(GUI here) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        home=here;
        println("Activated panel for bondlength adjustment");
        job=new BondLengths();
        job.pack();
        job.setVisible(true);
        bondLengthSet();
    }
    public void bondLengthSet() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        b01.setText(String.valueOf(cc1b));
        b02.setText(String.valueOf(cc2b));
        b03.setText(String.valueOf(cc3b));
        b04.setText(String.valueOf(ccab));
        b05.setText(String.valueOf(ch1b));
        b06.setText(String.valueOf(cn1b));
        b07.setText(String.valueOf(cn2b));
        b08.setText(String.valueOf(cn3b));
        b09.setText(String.valueOf(cnab));
        b10.setText(String.valueOf(co1b));
        b11.setText(String.valueOf(co2b));
        b12.setText(String.valueOf(coab));
        b13.setText(String.valueOf(nh1b));
        b14.setText(String.valueOf(oh1b));
    }
    public void actionPerformed(ActionEvent e) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        int call;
        String arg = (String)e.getActionCommand();
        if (arg.equals("Set")) {
            cc1b=BML.giveDouble(b01.getText(),1);
            cc2b=BML.giveDouble(b02.getText(),1);
            cc3b=BML.giveDouble(b03.getText(),1);
            ccab=BML.giveDouble(b04.getText(),1);
            ch1b=BML.giveDouble(b05.getText(),1);
            cn1b=BML.giveDouble(b06.getText(),1);
            cn2b=BML.giveDouble(b07.getText(),1);
            cn3b=BML.giveDouble(b08.getText(),1);
            cnab=BML.giveDouble(b09.getText(),1);
            co1b=BML.giveDouble(b10.getText(),1);
            co2b=BML.giveDouble(b11.getText(),1);
            coab=BML.giveDouble(b12.getText(),1);
            nh1b=BML.giveDouble(b13.getText(),1);
            oh1b=BML.giveDouble(b14.getText(),1);
            println("Bondlengths now reset");
        }
        else if (arg.equals("Close")) {
            job.setVisible(false);
        }
    }
}

