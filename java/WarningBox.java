import java.awt.*;
import java.awt.event.*;

public class WarningBox extends Dialog implements ActionListener {
            /*
**********************************************************************

dl_poly/java class to generate a warning box

copyright - daresbury laboratory
author    - w. smith march 2001

**********************************************************************
             */
    static GUI home;
    static Color back,fore,butn,butf;
    static Button pos,neg;
    static Font fontMain;

    public WarningBox(GUI here,String header,boolean modal) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        super(here,"Warning!",modal);
        setSize(400,70);
        home=here;
        back=new Color(255,0,0);
        fore=new Color(0,0,0);
        butn=new Color(255,255,0);
        butf=new Color(0,0,0);
        fontMain=home.fontMain;
        setBackground(back);
        setForeground(fore);
        setFont(fontMain);
        setLayout(new FlowLayout());

        // Warning label

        Label holdit = new Label("Warning - Do you wish to proceed?");
        add(holdit);

        // define the buttons

        pos=new Button("Yes");
        pos.setForeground(butf);
        pos.setBackground(butn);
        neg=new Button(" No ");
        neg.setForeground(butf);
        neg.setBackground(butn);
        add(pos);
        add(neg);

        // Register action buttons

        pos.addActionListener(this);
        neg.addActionListener(this);

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
        if (arg.equals("Yes")) {
            home.alert=true;
        }
        else if (arg.equals(" No ")) {
            home.alert=false;
        }
	dispose();
    }
}
