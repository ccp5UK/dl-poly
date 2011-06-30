import java.awt.*;
import java.io.*;
import java.awt.event.*;
import java.awt.print.*;
import java.awt.geom.*;
import javax.swing.*;

// Define the DL_POLY GUI Editor

public class Editor extends Basic {
        /*
*********************************************************************

main dl_poly/java Editor class

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */

    public static GUI home;
    public static Editor job;
    private GridBagLayout grd;
    private GridBagConstraints gbc;
    public static boolean safe;
    private static JButton AC,BC,CC,DC,EC,FC,GC,HC,IC,JC;
    private static JMenu Options;

    public Editor() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        super();

        // Set up Frame

        setTitle("DL_POLY Editor");

        setFont(fontMain);
        grd = new GridBagLayout();
        gbc = new GridBagConstraints();
        getContentPane().setLayout(grd);
        gbc.fill=GridBagConstraints.BOTH;

        // Define the Molecular Builder

        pane=new Builder(this);
        pane.setPreferredSize(new Dimension(TILEX,TILEY));
        pane.setBackground(art.scrn);
        pane.setForeground(art.scrf);
        fix(pane,grd,gbc,0,0,1,14,TILEX,TILEY);

        // Define the Editor buttons

        defEditorButtons();

        // Define menu bar

        defMenuBar();

    }

    // Main Constructor

    public Editor(GUI here) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        home=here;
        safe=true;

        if(config == null){
            if(cfgsav == null){
                config=new Config();
            }
            else{
                config=copyConfig(cfgsav);
            }
        }
        else{
            cfgsav=copyConfig(config);
        }

        // Set up Editor interface

        job = new Editor();
        job.pack();
        job.setVisible(true);
    }

    // define the editor buttons

    void defEditorButtons() {
        /*
*********************************************************************

dl_poly/java GUI routine to define the graphics buttons

copyrigh2t - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        JButton AA,AB,BA,BB,CA,CB,DA,DB,EA,EB;
        JButton FA,FB,GA,GB,HA,HB,IA,IB,JA,JB;

        AA=MyButton("New");
        fix(AA,grd,gbc,2,0,1,1,20,10);
        AB=MyButton("Edt");
        fix(AB,grd,gbc,3,0,1,1,20,10);
        AC=MyButton("Drw");
        fix(AC,grd,gbc,4,0,1,1,20,10);
        BA=MyButton("Rst");
        fix(BA,grd,gbc,2,1,1,1,20,10);
        BB=MyButton("Clr");
        fix(BB,grd,gbc,3,1,1,1,20,10);
        BC=MyButton("Lnk");
        fix(BC,grd,gbc,4,1,1,1,20,10);

        // second block of buttons

        CA=MyButton("Tx-");
        fix(CA,grd,gbc,2,3,1,1,20,10);
        CB=MyButton("Tx+");
        fix(CB,grd,gbc,3,3,1,1,20,10);
        CC=MyButton("Del");
        fix(CC,grd,gbc,4,3,1,1,20,10);
        DA=MyButton("Ty-");
        fix(DA,grd,gbc,2,4,1,1,20,10);
        DB=MyButton("Ty+");
        fix(DB,grd,gbc,3,4,1,1,20,10);
        DC=MyButton("ADH");
        fix(DC,grd,gbc,4,4,1,1,20,10);

        // third block of buttons

        EA=MyButton("Tz-");
        fix(EA,grd,gbc,2,6,1,1,20,10);
        EB=MyButton("Tz+");
        fix(EB,grd,gbc,3,6,1,1,20,10);
        EC=MyButton("Grp");
        fix(EC,grd,gbc,4,6,1,1,20,10);
        FA=MyButton("Rot");
        fix(FA,grd,gbc,2,7,1,1,20,10);
        FB=MyButton("Tra");
        fix(FB,grd,gbc,3,7,1,1,20,10);
        FC=MyButton("Opt");
        fix(FC,grd,gbc,4,7,1,1,20,10);

        // fourth block of buttons

        GA=MyButton("Rx-");
        fix(GA,grd,gbc,2,9,1,1,20,10);
        GB=MyButton("Rx+");
        fix(GB,grd,gbc,3,9,1,1,20,10);
        GC=MyButton("Sav");
        fix(GC,grd,gbc,4,9,1,1,20,10);
        HA=MyButton("Ry-");
        fix(HA,grd,gbc,2,10,1,1,20,10);
        HB=MyButton("Ry+");
        fix(HB,grd,gbc,3,10,1,1,20,10);
        HC=MyButton("Dup");
        fix(HC,grd,gbc,4,10,1,1,20,10);

        // fifth block of buttons

        IA=MyButton("Rz-");
        fix(IA,grd,gbc,2,12,1,1,20,10);
        IB=MyButton("Rz+");
        fix(IB,grd,gbc,3,12,1,1,20,10);
        IC=MyButton("Box");
        fix(IC,grd,gbc,4,12,1,1,20,10);
        JA=MyButton("H2O");
        fix(JA,grd,gbc,2,13,1,1,20,10);
        JB=MyButton("Bnd");
        fix(JB,grd,gbc,3,13,1,1,20,10);
        JC=MyButton("Frg");
        fix(JC,grd,gbc,4,13,1,1,20,10);

        // Enable/Disable active edit buttons

        AC.setEnabled(edit);
        BC.setEnabled(edit);
        CC.setEnabled(edit);
        DC.setEnabled(edit);
        EC.setEnabled(edit);
        FC.setEnabled(edit);
        GC.setEnabled(edit);
        HC.setEnabled(edit);
        IC.setEnabled(edit);
        JC.setEnabled(edit);

        // Register button events listeners

        AA.addActionListener(new EditorButtonHandler());
        AB.addActionListener(new EditorButtonHandler());
        AC.addActionListener(new EditorButtonHandler());
        BA.addActionListener(new EditorButtonHandler());
        BB.addActionListener(new EditorButtonHandler());
        BC.addActionListener(new EditorButtonHandler());
        CA.addActionListener(new EditorButtonHandler());
        CB.addActionListener(new EditorButtonHandler());
        CC.addActionListener(new EditorButtonHandler());
        DA.addActionListener(new EditorButtonHandler());
        DB.addActionListener(new EditorButtonHandler());
        DC.addActionListener(new EditorButtonHandler());
        EA.addActionListener(new EditorButtonHandler());
        EB.addActionListener(new EditorButtonHandler());
        EC.addActionListener(new EditorButtonHandler());
        FA.addActionListener(new EditorButtonHandler());
        FB.addActionListener(new EditorButtonHandler());
        FC.addActionListener(new EditorButtonHandler());
        GA.addActionListener(new EditorButtonHandler());
        GB.addActionListener(new EditorButtonHandler());
        GC.addActionListener(new EditorButtonHandler());
        HA.addActionListener(new EditorButtonHandler());
        HB.addActionListener(new EditorButtonHandler());
        HC.addActionListener(new EditorButtonHandler());
        IA.addActionListener(new EditorButtonHandler());
        IB.addActionListener(new EditorButtonHandler());
        IC.addActionListener(new EditorButtonHandler());
        JA.addActionListener(new EditorButtonHandler());
        JB.addActionListener(new EditorButtonHandler());
        JC.addActionListener(new EditorButtonHandler());

    }

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

        // Edit menu A

        JMenuItem itemAAa,itemAAb,itemAAc,itemAAd,itemAAe;
        JMenuItem itemAAf,itemAAg,itemAAh,itemAAi,itemAAj;
        JMenuItem itemAAk,itemAAl,itemAAm,itemAAn,itemAAo;
        JMenuItem itemAAp,itemAAq;
        JMenuItem itemABa,itemABb,itemABc,itemABd,itemABe;
        JMenuItem itemACa,itemACb,itemACc,itemACd,itemACe;
        JMenuItem itemACf,itemACg,itemACh,itemACi,itemACj;
        JMenuItem itemACk,itemABf;
        Options = MyMenu("Options");
        JMenu subOptionsAA = MyMenu("Atoms");
        subOptionsAA.add(itemAAa = MyMenuItem("H_"));
        subOptionsAA.add(itemAAb = MyMenuItem("C_3"));
        subOptionsAA.add(itemAAc = MyMenuItem("C_2"));
        subOptionsAA.add(itemAAd = MyMenuItem("C_1"));
        subOptionsAA.add(itemAAe = MyMenuItem("C_R"));
        subOptionsAA.add(itemAAf = MyMenuItem("O_3"));
        subOptionsAA.add(itemAAg = MyMenuItem("O_2"));
        subOptionsAA.add(itemAAh = MyMenuItem("N_3"));
        subOptionsAA.add(itemAAi = MyMenuItem("N_2"));
        subOptionsAA.add(itemAAj = MyMenuItem("N_1"));
        subOptionsAA.add(itemAAk = MyMenuItem("P_3"));
        subOptionsAA.add(itemAAl = MyMenuItem("P_2"));
        subOptionsAA.add(itemAAm = MyMenuItem("S_3"));
        subOptionsAA.add(itemAAn = MyMenuItem("S_2"));
        subOptionsAA.add(itemAAo = MyMenuItem("OW"));
        subOptionsAA.add(itemAAp = MyMenuItem("HW"));
        subOptionsAA.add(itemAAq = MyMenuItem("H__HB"));
        subOptionsAA.setForeground(Color.black);
        Options.add(subOptionsAA);
        JMenu subOptionsAB = MyMenu("Box");
        subOptionsAB.add(itemABa = MyMenuItem("None"));
        subOptionsAB.add(itemABb = MyMenuItem("Cubic"));
        subOptionsAB.add(itemABc = MyMenuItem("O-Rhombic"));
        subOptionsAB.add(itemABd = MyMenuItem("T-Octahedral"));
        subOptionsAB.add(itemABe = MyMenuItem("R-Dodecahedral"));
        subOptionsAB.add(itemABf = MyMenuItem("Hexagonal"));
        subOptionsAB.setForeground(Color.black);
        Options.add(subOptionsAB);
        JMenu subOptionsAC = MyMenu("Fragment");
        subOptionsAC.add(itemACa = MyMenuItem("Search"));
        subOptionsAC.add(itemACb = MyMenuItem("Alanine"));
        subOptionsAC.add(itemACc = MyMenuItem("Benzene"));
        subOptionsAC.add(itemACd = MyMenuItem("Glucose"));
        subOptionsAC.add(itemACe = MyMenuItem("i-Butane"));
        subOptionsAC.add(itemACf = MyMenuItem("Naphthalene"));
        subOptionsAC.add(itemACg = MyMenuItem("Styrene"));
        subOptionsAC.add(itemACh = MyMenuItem("c-Hexane"));
        subOptionsAC.add(itemACi = MyMenuItem("n-Butane"));
        subOptionsAC.add(itemACj = MyMenuItem("n-Decane"));
        subOptionsAC.add(itemACk = MyMenuItem("n-Hexane"));
        subOptionsAC.setForeground(Color.black);
        Options.add(subOptionsAC);
        Options.setForeground(art.fore);
        Options.setBackground(art.back);
        Options.setFont(fontMain);
        top.add(Options);
        itemAAa.addActionListener(new EditMenuHandler());
        itemAAb.addActionListener(new EditMenuHandler());
        itemAAc.addActionListener(new EditMenuHandler());
        itemAAd.addActionListener(new EditMenuHandler());
        itemAAe.addActionListener(new EditMenuHandler());
        itemAAf.addActionListener(new EditMenuHandler());
        itemAAg.addActionListener(new EditMenuHandler());
        itemAAh.addActionListener(new EditMenuHandler());
        itemAAi.addActionListener(new EditMenuHandler());
        itemAAj.addActionListener(new EditMenuHandler());
        itemAAk.addActionListener(new EditMenuHandler());
        itemAAl.addActionListener(new EditMenuHandler());
        itemAAm.addActionListener(new EditMenuHandler());
        itemAAn.addActionListener(new EditMenuHandler());
        itemAAp.addActionListener(new EditMenuHandler());
        itemAAq.addActionListener(new EditMenuHandler());
        itemAAo.addActionListener(new EditMenuHandler());
        itemABa.addActionListener(new EditMenuHandler());
        itemABb.addActionListener(new EditMenuHandler());
        itemABc.addActionListener(new EditMenuHandler());
        itemABd.addActionListener(new EditMenuHandler());
        itemABe.addActionListener(new EditMenuHandler());
        itemABf.addActionListener(new EditMenuHandler());
        itemACa.addActionListener(new EditMenuHandler());
        itemACb.addActionListener(new EditMenuHandler());
        itemACc.addActionListener(new EditMenuHandler());
        itemACd.addActionListener(new EditMenuHandler());
        itemACe.addActionListener(new EditMenuHandler());
        itemACf.addActionListener(new EditMenuHandler());
        itemACg.addActionListener(new EditMenuHandler());
        itemACh.addActionListener(new EditMenuHandler());
        itemACi.addActionListener(new EditMenuHandler());
        itemACj.addActionListener(new EditMenuHandler());
        itemACk.addActionListener(new EditMenuHandler());


        // Invoke menu bar

        setJMenuBar(top);

        // Enable editor menu

        Options.setEnabled(edit);

    }

    JButton MyButton(String s) {
        JButton mine=new JButton(s);
        mine.setBackground(art.butn);
        mine.setForeground(art.butf);
        return mine;
    }

    // Handle Button events

    class EditorButtonHandler implements ActionListener {
        /*
*********************************************************************

dl_poly/java GUI class to handle button events

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

            if (e.getSource() instanceof JButton) {
                if(arg.equals("New")) {
                    config=getConfig(home,ftype);
                    if(config != null)
                        cfgsav=copyConfig(config);
                    pane.restore();
                }
                else if(arg.equals("Clr")) {
                    if(edit && !safe) saveEdit(0);
                    config=null;
                    cfgsav=null;
                    pane.newBuild();
                    pane.restore();
                }
                else if(arg.equals("Tx+"))
                    switchTxp();
                else if(arg.equals("Tx-"))
                    switchTxm();
                else if(arg.equals("Ty+"))
                    switchTyp();
                else if(arg.equals("Ty-"))
                    switchTym();
                else if(arg.equals("Tz+"))
                    switchTzp();
                else if(arg.equals("Tz-"))
                    switchTzm();
                else if(arg.equals("Tra"))
                    switchTranslate();
                else if(arg.equals("Rx+"))
                    switchRxp();
                else if(arg.equals("Rx-"))
                    switchRxm();
                else if(arg.equals("Ry+"))
                    switchRyp();
                else if(arg.equals("Ry-"))
                    switchRym();
                else if(arg.equals("Rz+"))
                    switchRzp();
                else if(arg.equals("Rz-"))
                    switchRzm();
                else if(arg.equals("Rot"))
                    switchRotate();
                else if(arg.equals("Rst"))
                    switchRestore();
                else if(arg.equals("Edt"))
                    switchEdit();
                else if(arg.equals("Drw"))
                    switchDraw();
                else if(arg.equals("Bnd"))
                    switchBonds();
                else if(arg.equals("Grp"))
                    switchGroup();
                else if(arg.equals("Opt"))
                    switchOpt();
                else if(arg.equals("Sav"))
                    saveEdit(1);
                else if(arg.equals("Del"))
                    switchDelete();
                else if(arg.equals("ADH"))
                    switchAddHydrogen();
                else if(arg.equals("Lnk"))
                    switchAddLinks();
                else if(arg.equals("H2O"))
                    showWater();
                else if(arg.equals("Dup"))
                    switchDuplicate();
                else if(arg.equals("Box"))
                    switchBoxEditor();
                else if(arg.equals("Frg"))
                    switchAddFragment();
            }
        }
    }

    class EditMenuHandler implements ActionListener {
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
                if(arg.equals("C_1")) {
                    setDrawAtom("C_1");
                }
                else if(arg.equals("C_2")) {
                    setDrawAtom("C_2");
                }
                else if(arg.equals("C_R")) {
                    setDrawAtom("C_R");
                }
                else if(arg.equals("C_3")) {
                    setDrawAtom("C_3");
                }
                else if(arg.equals("O_2")) {
                    setDrawAtom("O_2");
                }
                else if(arg.equals("O_3")) {
                    setDrawAtom("O_3");
                }
                else if(arg.equals("H_")) {
                    setDrawAtom("H_");
                }
                else if(arg.equals("N_1")) {
                    setDrawAtom("N_1");
                }
                else if(arg.equals("N_2")) {
                    setDrawAtom("N_2");
                }
                else if(arg.equals("N_3")) {
                    setDrawAtom("N_3");
                }
                else if(arg.equals("S_2")) {
                    setDrawAtom("S_2");
                }
                else if(arg.equals("S_3")) {
                    setDrawAtom("S_3");
                }
                else if(arg.equals("P_2")) {
                    setDrawAtom("P_2");
                }
                else if(arg.equals("P_3")) {
                    setDrawAtom("P_3");
                }
                else if(arg.equals("OW")) {
                    setDrawAtom("OW");
                }
                else if(arg.equals("HW")) {
                    setDrawAtom("HW");
                }
                else if(arg.equals("H__HB")) {
                    setDrawAtom("H__HB");
                }
                else if(arg.equals("None")) {
                    setBoxType("NON");
                }
                else if(arg.equals("Cubic")) {
                    setBoxType("CUB");
                }
                else if(arg.equals("O-Rhombic")) {
                    setBoxType("ORH");
                }
                else if(arg.equals("T-Octahedral")) {
                    setBoxType("OCT");
                }
                else if(arg.equals("R-Dodecahedral")) {
                    setBoxType("DEC");
                }
                else if(arg.equals("Hexagonal")) {
                    setBoxType("HEX");
                }
                else if(arg.equals("Search")) {
                    setFragmentType("SELECTED");
                }
                else if(arg.equals("Alanine")) {
                    setFragmentType("ALANINE");
                }
                else if(arg.equals("Benzene")) {
                    setFragmentType("BENZENE");
                }
                else if(arg.equals("Naphthalene")) {
                    setFragmentType("NAPHTHALENE");
                }
                else if(arg.equals("n-Hexane")) {
                    setFragmentType("n_HEXANE");
                }
                else if(arg.equals("n-Decane")) {
                    setFragmentType("n_DECANE");
                }
                else if(arg.equals("n-Butane")) {
                    setFragmentType("n_BUTANE");
                }
                else if(arg.equals("c-Hexane")) {
                    setFragmentType("c_HEXANE");
                }
                else if(arg.equals("i-Butane")) {
                    setFragmentType("i_BUTANE");
                }
                else if(arg.equals("Styrene")) {
                    setFragmentType("STYRENE");
                }
                else if(arg.equals("Glucose")) {
                    setFragmentType("GLUCOSE");
                }
            }
        }
    }

    void editRestore() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        pane.fac=0.3;
        pane.oper=0;
        pane.ngroup=0;
        pane.mark0=-1;
        pane.mark1=-1;
        pane.mark2=-1;
        //pane.scale=37.5;
        if(cfgsav == null) {
            println("No edit backup available");
        }
        else {
            config=copyConfig(cfgsav);
            println("Edit backup restored");
            safe=true;
            pane.lpen=false;
            pane.link=false;
            pane.news="NULL";
            pane.repaint();
        }
    }

    void saveEdit(int n) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        String fname="CFGEDT";

        if(n > 0)safe=false;

        // Save backup data file

        if(!safe) {
            if(n > 0){
                fname+="."+numsav;
                numsav++;
            }
            config.structure=new Structure(config);
            cfgsav=copyConfig(config);
            if(config.configWrite(fname)){
                safe=true;
            }
            else{
                println("Error - failed to write Edit backup file: "+fname);
            }
        }
    }

    void showWater() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        if(pane.water) {
            println("Water visibility is OFF");
            pane.water=false;
        }
        else {
            println("Water visibility is ON");
            pane.water=true;
        }
        pane.repaint();
    }

    void switchEdit() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        if(edit) {
            edit=false;
            pane.fragment=null;
            saveEdit(0);
        }
        else {
            if(cfgsav!=null)
                config=copyConfig(cfgsav);
            edit=true;
            pane.fragment=null;
            pane.oper=0;
            pane.news="NULL";
            setDrawAtom("C_3");
        }
        Options.setEnabled(edit);
        AC.setEnabled(edit);
        BC.setEnabled(edit);
        CC.setEnabled(edit);
        DC.setEnabled(edit);
        EC.setEnabled(edit);
        FC.setEnabled(edit);
        GC.setEnabled(edit);
        HC.setEnabled(edit);
        IC.setEnabled(edit);
        JC.setEnabled(edit);

        pane.shx=0;
        pane.shy=0;
        pane.restore();
    }

    void switchDuplicate() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        if(edit){
            pane.mark0=-1;
            pane.mark1=-1;
            pane.mark2=-1;
            if(pane.oper == 11) {
                pane.oper=0;
                pane.news="NULL";
            }
            else {
                pane.oper=11;
                pane.news="DUPLICATE";
            }
            pane.repaint();
        }
    }

    void switchBoxEditor() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        if(edit){
            pane.shx=0;
            pane.shy=0;
            pane.mark0=-1;
            pane.mark1=-1;
            pane.mark2=-1;
            if(pane.oper == 12) {
                pane.oper=0;
                pane.news="NULL";
            }
            else {
                pane.oper=12;
                pane.defineBox();
                pane.news="EDIT BOX";
            }
            pane.repaint();
        }
    }

    void switchDraw() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        if(edit){
            pane.mark0=-1;
            pane.mark1=-1;
            pane.mark2=-1;
            if(pane.oper == 1) {
                pane.oper=0;
                pane.lpen=false;
                pane.news="NULL";
            }
            else {
                pane.oper=1;
                pane.lpen=true;
                pane.link=false;
                pane.news="DRAW";
            }
            pane.repaint();
        }
    }

    void switchBonds() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        if(showbonds)
            showbonds=false;
        else
            showbonds=true;

        pane.repaint();
    }

    void switchAddFragment() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        if(edit){
            if(pane.oper == 13) {
                pane.oper=0;
                pane.news="NULL";
            }
            else {
                pane.oper=13;
                pane.news="FRAGMENT";
            }
            pane.repaint();
        }
    }

    void switchAddLinks() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        if(edit){
            pane.mark0=-1;
            pane.mark1=-1;
            pane.mark2=-1;
            if(pane.oper == 2) {
                pane.oper=0;
                pane.news="NULL";
            }
            else {
                pane.oper=2;
                pane.news="LINK";
            }
            pane.repaint();
        }
    }

    void switchAddHydrogen() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        if(edit){
            pane.mark0=-1;
            pane.mark1=-1;
            pane.mark2=-1;
            if(pane.oper == 10) {
                pane.oper=0;
                pane.news="NULL";
            }
            else {
                pane.oper=10;
                pane.news="ADD/DELETE HYDROGEN";
                if(pane.ngroup != 0)
                    pane.addDeleteHydrogen();
            }
            pane.repaint();
        }
    }

    void switchGroup() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        if(edit){
            pane.mark0=-1;
            pane.mark1=-1;
            pane.mark2=-1;
            if(pane.oper == 5) {
                pane.oper=0;
                pane.news="NULL";
            }
            else {
                pane.oper=5;
                pane.news="GROUP";
            }
            pane.repaint();
        }
    }

    void switchOpt() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        if(edit){
            if(pane.oper == 9 && pane.keyopt == 999) {
                pane.oper=0;
                pane.news="NULL";
            }
            else {
                pane.oper=9;
                pane.news="OPTIMIZATION";
                pane.Optimize();
            }
            pane.repaint();
        }
    }

    void switchDelete() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        if(edit){
            if(pane.oper == 3) {
                pane.oper=0;
                pane.news="NULL";
            }
            else {
                pane.oper=3;
                pane.mark0=-1;
                pane.mark1=-1;
                pane.mark2=-1;
                pane.news="DELETE";
                if(pane.ngroup != 0) pane.deleteGroup();
            }
            pane.repaint();
        }
    }

    void switchTzm() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        if(edit) {
            if(pane.oper == 6) {
                if(pane.activity.equals("Tz-")) {
                    pane.activity="none";
                    pane.news="SHIFT";
                }
                else {
                    pane.activity="Tz-";
                    pane.news="SHIFT Z";
                }
            }
            else {
                if(pane.oper == 40) {
                    pane.oper=0;
                    pane.activity="none";
                    pane.news="NULL";
                }
                else {
                    pane.oper=40;
                    pane.activity="Tz-";
                    pane.news="SHIFT Z-";
                }
            }
        }
        else {
            pane.displace(2,-incz);
        }
        pane.repaint();
    }

    void switchTzp() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        if(edit) {
            if(pane.oper == 6) {
                if(pane.activity.equals("Tz+")) {
                    pane.activity="none";
                    pane.news="SHIFT";
                }
                else {
                    pane.activity="Tz+";
                    pane.news="SHIFT Z";
                }
            }
            else {
                if(pane.oper == 41) {
                    pane.oper=0;
                    pane.activity="none";
                    pane.news="NULL";
                }
                else {
                    pane.oper=41;
                    pane.activity="Tz+";
                    pane.news="SHIFT Z+";
                }
            }
        }
        else {
            pane.displace(2,incz);
        }
        pane.repaint();
    }

    void switchTxm() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        if(edit) {
            if(pane.oper == 6) {
                if(pane.activity.equals("Tx-")) {
                    pane.activity="none";
                    pane.news="SHIFT";
                }
                else {
                    pane.activity="Tx-";
                    pane.news="SHIFT X";
                }
            }
            else {
                if(pane.oper == 42) {
                    pane.oper=0;
                    pane.activity="none";
                    pane.news="NULL";
                }
                else {
                    pane.oper=42;
                    pane.activity="Tx-";
                    pane.news="SHIFT X-";
                }
            }
        }
        else {
            pane.displace(0,-incx);
        }
        pane.repaint();
    }

    void switchTxp() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        if(edit) {
            if(pane.oper == 6) {
                if(pane.activity.equals("Tx+")) {
                    pane.activity="none";
                    pane.news="SHIFT";
                }
                else {
                    pane.activity="Tx+";
                    pane.news="SHIFT X";
                }
            }
            else {
                if(pane.oper == 43) {
                    pane.oper=0;
                    pane.activity="none";
                    pane.news="NULL";
                }
                else {
                    pane.oper=43;
                    pane.activity="Tx+";
                    pane.news="SHIFT X+";
                }
            }
        }
        else {
            pane.displace(0,incx);
        }
        pane.repaint();
    }

    void switchTym() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        if(edit) {
            if(pane.oper == 6) {
                if(pane.activity.equals("Ty-")) {
                    pane.activity="none";
                    pane.news="SHIFT";
                }
                else {
                    pane.activity="Ty-";
                    pane.news="ZOOM";
                }
            }
            else {
                if(pane.oper == 44) {
                    pane.oper=0;
                    pane.activity="none";
                    pane.news="NULL";
                }
                else {
                    pane.oper=44;
                    pane.activity="Ty-";
                    pane.news="SHIFT Y-";
                }
            }
        }
        else {
            pane.zoom(-1);
        }
        pane.repaint();
    }

    void switchTyp() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        if(edit) {
            if(pane.oper == 6) {
                if(pane.activity.equals("Ty+")) {
                    pane.activity="none";
                    pane.news="SHIFT";
                }
                else {
                    pane.activity="Ty+";
                    pane.news="PAN";
                }
            }
            else {
                if(pane.oper == 45) {
                    pane.oper=0;
                    pane.activity="none";
                    pane.news="NULL";
                }
                else {
                    pane.oper=45;
                    pane.activity="Ty+";
                    pane.news="SHIFT Y+";
                }
            }
        }
        else {
            pane.zoom(1);
        }
        pane.repaint();
    }

    void switchTranslate() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        pane.activity="none";
        if(pane.oper == 6) {
            pane.oper=0;
            pane.news="NULL";
        }
        else {
            pane.oper=6;
            pane.news="SHIFT";
        }
        pane.repaint();
    }

    void switchRzm() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        if(edit) {
            if(pane.oper == 7) {
                if(pane.activity.equals("Rz-")) {
                    pane.activity="none";
                    pane.news="ROTATE";
                }
                else {
                    pane.activity="Rz-";
                    pane.news="ROTATE Z";
                }
            }
            else {
                if(pane.oper == 80) {
                    pane.oper=0;
                    pane.activity="none";
                    pane.news="NULL";
                }
                else {
                    pane.oper=80;
                    pane.activity="Rz-";
                    pane.news="ROTATE Z-";
                }
            }
        }
        else {
            pane.rotate(0,1,rotcos,-rotsin);
        }
        pane.repaint();
    }

    void switchRzp() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        if(edit) {
            if(pane.oper == 7) {
                if(pane.activity.equals("Rz+")) {
                    pane.activity="none";
                    pane.news="ROTATE";
                }
                else {
                    pane.activity="Rz+";
                    pane.news="ROTATE Z";
                }
            }
            else {
                if(pane.oper == 81) {
                    pane.oper=0;
                    pane.activity="none";
                    pane.news="NULL";
                }
                else {
                    pane.oper=81;
                    pane.activity="Rz+";
                    pane.news="ROTATE Z+";
                }
            }
        }
        else {
            pane.rotate(0,1,rotcos,rotsin);
        }
        pane.repaint();
    }

    void switchRxm() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        if(edit) {
            if(pane.oper == 7) {
                if(pane.activity.equals("Rx-")) {
                    pane.activity="none";
                    pane.news="ROTATE";
                }
                else {
                    pane.activity="Rx-";
                    pane.news="ROTATE X";
                }
            }
            else {
                if(pane.oper == 82) {
                    pane.oper=0;
                    pane.activity="none";
                    pane.news="NULL";
                }
                else {
                    pane.oper=82;
                    pane.activity="Rx-";
                    pane.news="ROTATE X-";
                }
            }
        }
        else {
            pane.rotate(1,2,rotcos,-rotsin);
        }
        pane.repaint();
    }

    void switchRxp() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        if(edit) {
            if(pane.oper == 7) {
                if(pane.activity.equals("Rx+")) {
                    pane.activity="none";
                    pane.news="ROTATE";
                }
                else {
                    pane.activity="Rx+";
                    pane.news="ROTATE X";
                }
            }
            else {
                if(pane.oper == 83) {
                    pane.oper=0;
                    pane.activity="none";
                    pane.news="NULL";
                }
                else {
                    pane.oper=83;
                    pane.activity="Rx+";
                    pane.news="ROTATE X+";
                }
            }
        }
        else {
            pane.rotate(1,2,rotcos,rotsin);
        }
        pane.repaint();
    }

    void switchRym() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        if(edit) {
            if(pane.oper == 7) {
                if(pane.activity.equals("Ry-")) {
                    pane.activity="none";
                    pane.news="ROTATE";
                }
                else {
                    pane.activity="Ry-";
                    pane.news="ROTATE Y";
                }
            }
            else {
                if(pane.oper == 84) {
                    pane.oper=0;
                    pane.activity="none";
                    pane.news="NULL";
                }
                else {
                    pane.oper=84;
                    pane.activity="Ry-";
                    pane.news="ROTATE Y-";
                }
            }
        }
        else {
            pane.rotate(2,0,rotcos,-rotsin);
        }
        pane.repaint();
    }

    void switchRyp() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        if(edit) {
            if(pane.oper == 7) {
                if(pane.activity.equals("Ry+")) {
                    pane.activity="none";
                    pane.news="ROTATE";
                }
                else {
                    pane.activity="Ry+";
                    pane.news="ROTATE Y";
                }
            }
            else {
                if(pane.oper == 85) {
                    pane.oper=0;
                    pane.activity="none";
                    pane.news="NULL";
                }
                else {
                    pane.oper=85;
                    pane.activity="Ry+";
                    pane.news="ROTATE Y+";
                }
            }
        }
        else {
            pane.rotate(2,0,rotcos,rotsin);
        }
        pane.repaint();
    }

    void switchRotate() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        pane.activity="none";
        if(pane.oper == 7) {
            pane.oper=0;
            pane.news="NULL";
        }
        else {
            pane.oper=7;
            pane.news="ROTATE";
        }
        pane.repaint();
    }

    void switchRestore() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        if(edit)
            editRestore();
        else{
            config=copyConfig(cfgsav);
            pane.restore();
        }
    }

    void setDrawAtom(String atm) {
         /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        if(edit){
            pane.atom=new Element(atm);
            pane.repaint();
        }
    }

    void setFragmentType(String fra) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        boolean chk;

        if(edit){
            pane.frgname=fra;
            if(fra.equals("SELECTED")) {
                pane.fragment=getConfig(home,"CFG");
            }
            else {
                pane.fragment=new Config();
                //if(!pane.fragment.rdCFG("../java/FRGLIB/"+fra+".FRG"))
                if(!pane.fragment.rdFRG("FRGLIB/"+fra+".FRG"))
                    pane.fragment=null;
            }
            pane.repaint();
        }
    }

    void setBoxType(String box) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        if(box.equals("NON"))
            pane.boxtyp=0;
        else if(box.equals("CUB"))
            pane.boxtyp=1;
        else if(box.equals("ORH"))
            pane.boxtyp=2;
        else if(box.equals("OCT"))
            pane.boxtyp=4;
        else if(box.equals("DEC"))
            pane.boxtyp=5;
        else if(box.equals("HEX"))
            pane.boxtyp=7;
        pane.repaint();
    }

    void closeEditor() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */

            if(edit && !safe) saveEdit(0);
            edit=false;
            if(config!= null && config.natms == 0) config=null;
            job.setVisible(false);
    }

}
