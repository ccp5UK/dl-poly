import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;

// Define the graphical User Interface

class MakeChain extends Basic implements ActionListener {
        /*
*********************************************************************

dl_poly/java GUI class to make chain polymer CONFIG files

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
    public static GUI home;
    public static MakeChain job;
    private static JButton make,close;
    private static JLabel lab1,lab2,lab3,lab4,lab5,pass;
    private static JTextField nc,area,eon,gapz;
    private static JCheckBox chk1,chk2;
    private static JComboBox<String> head;
    private static boolean flip,twin;
    private static String headgroup;
    private static int keyhed,ncarbons,nethos;
    private static double zgap,harea;
    private static double beta=.615479708;
    private static double[] rot,uuu,vvv;
    private static double[] xbs,ybs,zbs;
    private static double[] cell;
    private static double[][] xyz;
    private static Element[] atoms;

    public MakeChain(){
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        super();
        setTitle("Make Chain");

        getContentPane().setBackground(art.back);
        getContentPane().setForeground(art.fore);
        setDefaultCloseOperation(DISPOSE_ON_CLOSE);
        setFont(fontMain);
        GridBagLayout grd = new GridBagLayout();
        GridBagConstraints gbc = new GridBagConstraints();
        getContentPane().setLayout(grd);

        gbc.fill=GridBagConstraints.BOTH;

        // Define the Make button

        make = new JButton("Make");
        make.setBackground(art.butn);
        make.setForeground(art.butf);
        fix(make,grd,gbc,0,0,1,1);

        // Define the Close button

        close = new JButton("Close");
        close.setBackground(art.butn);
        close.setForeground(art.butf);
        fix(close,grd,gbc,2,0,1,1);

        // Spacer label

        pass = new JLabel(" ");
        fix(pass,grd,gbc,0,1,2,1);

        // Instruction label 1

        lab1 = new JLabel("Number of C atoms:",JLabel.LEFT);
        fix(lab1,grd,gbc,0,2,2,1);

        // Number of C atoms

        nc = new JTextField(6);
        nc.setBackground(art.scrn);
        nc.setForeground(art.scrf);
        fix(nc,grd,gbc,2,2,1,1);

        // Instruction label 2

        lab2 = new JLabel("Headgroup Area (A^2)",JLabel.LEFT);
        fix(lab2,grd,gbc,0,3,2,1);

        // Head group area

        area = new JTextField(6);
        area.setBackground(art.scrn);
        area.setForeground(art.scrf);
        fix(area,grd,gbc,2,3,1,1);

        // Instruction label 3

        lab3 = new JLabel("Head Group:",JLabel.LEFT);
        fix(lab3,grd,gbc,0,4,2,1);

        // Head group choice

        head = new JComboBox<String>();
        head.setBackground(art.scrn);
        head.setForeground(art.scrf);
        head.addItem("None");
        head.addItem("Soap");
        head.addItem("Carboxy");
        head.addItem("Phenol");
        head.addItem("TAB");
        head.addItem("(EO)n");
        fix(head,grd,gbc,2,4,1,1);

        // Instruction label 4

        lab4 = new JLabel("Number of (EO)n groups:",JLabel.LEFT);
        fix(lab4,grd,gbc,0,5,2,1);

        // Number of (EO)n groups

        eon = new JTextField(6);
        eon.setBackground(art.scrn);
        eon.setForeground(art.scrf);
        fix(eon,grd,gbc,2,5,1,1);

        // Twin checkbox

        chk1 = new JCheckBox("Twin");
        chk1.setBackground(art.back);
        chk1.setForeground(art.fore);
        fix(chk1,grd,gbc,0,6,1,1);

        // Flip checkbox

        chk2 = new JCheckBox("Flip");
        chk2.setBackground(art.back);
        chk2.setForeground(art.fore);
        fix(chk2,grd,gbc,2,6,1,1);

        // Instruction label 5

        lab5 = new JLabel("Z - gap (A):",JLabel.LEFT);
        fix(lab5,grd,gbc,0,7,2,1);

        // Z - gap between twins

        gapz = new JTextField(6);
        gapz.setBackground(art.scrn);
        gapz.setForeground(art.scrf);
        fix(gapz,grd,gbc,2,7,1,1);

        // Register action buttons

        make.addActionListener(this);
        close.addActionListener(this);

    }

    // Constructor method

    public MakeChain(GUI here) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */

        println("Activated panel for making chain polymer CONFIG file");
        home=here;

        // Define arrays

        uuu=new double[3];
        vvv=new double[3];
        rot=new double[9];
        xbs=new double[4];
        ybs=new double[4];
        zbs=new double[4];

        // Set up Panel

        job = new MakeChain();
        job.pack();
        job.setVisible(true);
        setValues();
    }

    // Set initial values

    void setValues() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        // set default values

        keyhed=0;
        ncarbons=1;
        nethos=1;
        zgap=10.0;
        harea=25.0;
        flip=false;
        twin=false;

        // define tetrahedral groups

        setTetra();

        nc.setText(String.valueOf(ncarbons));
        area.setText(String.valueOf(harea));
        eon.setText(String.valueOf(nethos));
        gapz.setText(String.valueOf(zgap));
        chk1.setSelected(twin);
        chk2.setSelected(flip);

    }
    void setTetra() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        // define tetrahedral group

        xbs[0]=0.577350269;
        ybs[0]=0.0;
        zbs[0]=-0.816496581;

        xbs[1]=-0.577350269;
        ybs[1]=0.816496581;
        zbs[1]=0.0;

        xbs[2]=-0.577350269;
        ybs[2]=-0.816496581;
        zbs[2]=0.0;

        xbs[3]=0.577350269;
        ybs[3]=0.0;
        zbs[3]=0.816496581;
    }
    int chain() {
        /*
*********************************************************************

dl_poly/java utility for generating a linear polymer chain

copyright daresbury laboratory
author    w.smith 2011

*********************************************************************
         */

        double size,xx0,yy0,zz0,dens,disp,gap,base,cc23,co23;
        int imcon,natms;
        cell=new double[9];

        headgroup="";

        // start position of chain

        xx0=0.0;
        yy0=0.0;
        zz0=0.0;

        // select head group

        if(keyhed==1) {
            headgroup="Soap";
            cc23=0.5*(cc1b+cc2b);
            natms=soap();
        }
        else if(keyhed==2) {
            headgroup="Carboxy";
            cc23=0.5*(cc1b+cc2b);
            co23=0.5*(co1b+co2b);
            natms=carboxy();
        }
        else if(keyhed==3) {
            headgroup="Phenol";
            cc23=0.5*(cc1b+ccab);
            natms=phenol();
        }
        else if(keyhed==4) {
            headgroup="Trimethylammonium Bromide";
            natms=trimethyl();
        }
        else if(keyhed==5) {
            headgroup="Ethoxy(n)";
            natms=polyoxyeth();
        }
        else {
            natms=3*ncarbons+2;
            if(twin)natms*=2;
            atoms=new Element[natms];
            xyz=new double[3][natms];
            headgroup="No";
            atoms[0]=new Element("H_");
            xyz[0][0]=ch1b*xbs[0]+xx0;
            xyz[1][0]=ch1b*ybs[0]+yy0;
            xyz[2][0]=ch1b*zbs[0]+zz0;
            natms=1;
        }

        // attach chain to headgroup

        for(int i=0;i<ncarbons;i++) {
            atoms[natms]=new Element("C_3");
            xyz[0][natms]=xx0;
            xyz[1][natms]=yy0;
            xyz[2][natms]=zz0;

            atoms[natms+1]=new Element("H_");
            xyz[0][natms+1]=ch1b*xbs[1]+xx0;
            xyz[1][natms+1]=ch1b*ybs[1]+yy0;
            xyz[2][natms+1]=ch1b*zbs[1]+zz0;

            atoms[natms+2]=new Element("H_");
            xyz[0][natms+2]=ch1b*xbs[2]+xx0;
            xyz[1][natms+2]=ch1b*ybs[2]+yy0;
            xyz[2][natms+2]=ch1b*zbs[2]+zz0;

            atoms[natms+3]=new Element("H_");
            xyz[0][natms+3]=ch1b*xbs[3]+xx0;
            xyz[1][natms+3]=ch1b*ybs[3]+yy0;
            xyz[2][natms+3]=ch1b*zbs[3]+zz0;

            // growth direction vector

            xx0=cc1b*xbs[3]+xx0;
            yy0=cc1b*ybs[3]+yy0;
            zz0=cc1b*zbs[3]+zz0;
            xbs[0]=-xbs[0];
            xbs[1]=-xbs[1];
            xbs[2]=-xbs[2];
            xbs[3]=-xbs[3];

            natms+=3;
        }

        natms++;

        // determine simulation cell dimensions

        size=Math.sqrt(harea/0.866025403);
        cell[0]=size*0.866025403;
        cell[1]=-size*0.5;
        cell[2]=0.0;
        cell[3]=size*0.866025403;
        cell[4]=size*0.5;
        cell[5]=0.0;
        cell[6]=0.0;
        cell[7]=0.0;
        cell[8]=xyz[2][natms-1]-xyz[2][0]+0.5*size;

        // set first atom of chain to z=0

        base=xyz[2][0];

        for(int i=0;i<natms;i++) {
            xyz[2][i]=xyz[2][i]-base;
        }

        // flip structure if required

        if(flip) {
            base=xyz[2][natms-1];
            for(int i=0;i<natms;i++) {
                xyz[1][i]=-xyz[1][i];
                xyz[2][i]=-xyz[2][i]+base;
            }
        }

        // twin the chains in the Z direction

        if(twin) {
            gap=0.5*zgap;
            for(int i=0;i<natms;i++) {
                xyz[2][i]+=gap;
                atoms[i+natms]=new Element(atoms[i].zsym);
                xyz[0][i+natms]=-xyz[0][i];
                xyz[1][i+natms]= xyz[1][i];
                xyz[2][i+natms]=-xyz[2][i];
            }
            natms=2*natms;
            cell[8]=2.0*cell[8]+zgap;
        } else {
            for(int i=0;i<natms;i++)
                xyz[2][i]-=cell[8]/2.0;
        }
        println("Chain terminated successfully");
        println("Number of atoms generated:"+BML.fmt(natms,8));

        imcon=3;

        // Create Config object

        config =new Config();
        config.natms=natms;
        config.pbc.imcon=imcon;
        config.atoms=atoms;
        config.pbc.cell=cell;
        config.xyz=xyz;
        config.title="Chain polymer with"+BML.fmt(ncarbons,4)+" units and "+headgroup+" Head Group";
        config.pbc.buildBoundary(config.pbc.imcon);
        config.structure=new Structure(config);
	cfgsav=copyConfig(config);

        // write CONFIG file

        fname="CFGCHN."+String.valueOf(numchn);
        if(!config.configWrite(fname)) return -4;
        numchn++;

        // Draw structure

        if(!editor.isVisible())
            editor.showEditor();
        editor.pane.restore();

        return 0;
    }

    // make the chosen head group

    int soap()  {
        /*
**********************************************************************

dl_poly/java utility to make the head of a hydrocarbon chain
this one is for soap

copyright daresbury laboratory
author w.smith 2011

**********************************************************************
         */
        int natms;
        boolean op=false;

        natms=3*ncarbons+5;
        if(twin)natms*=2;
        atoms=new Element[natms];
        xyz=new double[3][natms];

        atoms[0]=new Element("Na");
        xyz[0][0]=0.0;
        xyz[1][0]=0.0;
        xyz[2][0]=-3.0*cc1b;

        atoms[1]=new Element("C_2");
        xyz[0][1]=0.0;
        xyz[1][1]=0.0;
        xyz[2][1]=-cc1b;

        atoms[2]=new Element("O_R");
        xyz[0][2]= coab*0.866025403;
        xyz[1][2]=0.0;
        xyz[2][2]=-cc1b-coab*0.5;

        atoms[3]=new Element("O_R");
        xyz[0][3]=-coab*0.866025403;
        xyz[1][3]=0.0;
        xyz[2][3]=-cc1b-coab*0.5;

        AML.euler(0.0,-beta,0.0,rot);

        for(int i=0;i<4;i++) {
            uuu[0]=xyz[0][i];
            uuu[1]=xyz[1][i];
            uuu[2]=xyz[2][i];
            AML.rotate(op,uuu,vvv,rot);
            xyz[0][i]=vvv[0];
            xyz[1][i]=vvv[1];
            xyz[2][i]=vvv[2];
        }
        return 4;
    }

    int carboxy() {
        /*
**********************************************************************

dl_poly/java utility to make the head of a hydrocarbon chain
this one is for carboxy head group

copyright daresbury laboratory
author w.smith 2011

**********************************************************************
         */
        int natms;
        boolean op=false;

        natms=3*ncarbons+5;
        if(twin)natms*=2;
        atoms=new Element[natms];
        xyz=new double[3][natms];

        atoms[0]=new Element("H_");
        xyz[0][0]= co1b*0.866025403;
        xyz[1][0]=0.0;
        xyz[2][0]=-cc1b-co1b*0.5-oh1b;

        atoms[1]=new Element("C_2");
        xyz[0][1]=0.0;
        xyz[1][1]=0.0;
        xyz[2][1]=-cc1b;

        atoms[2]=new Element("O_3");
        xyz[0][2]= co1b*0.866025403;
        xyz[1][2]=0.0;
        xyz[2][2]=-cc1b-co1b*0.5;

        atoms[3]=new Element("O_2");
        xyz[0][3]=-co2b*0.866025403;
        xyz[1][3]=0.0;
        xyz[2][3]=-cc1b-co2b*0.5;

        AML.euler(0.0,-beta,0.0,rot);

        for(int i=0;i<4;i++) {
            uuu[0]=xyz[0][i];
            uuu[1]=xyz[1][i];
            uuu[2]=xyz[2][i];
            AML.rotate(op,uuu,vvv,rot);
            xyz[0][i]=vvv[0];
            xyz[1][i]=vvv[1];
            xyz[2][i]=vvv[2];
        }
        return 4;
    }

    int trimethyl() {
        /*
*********************************************************************

dl_poly/java utility to make the head of a hydrocarbon chain
this one is for trimethylammonium

copyright daresbury laboratory
author w.smith 2011

**********************************************************************
         */
        int natms;
        boolean op=false;
        double cgam=0.942809041;
        double sgam=0.333333333;

        natms=3*ncarbons+15;
        if(twin)natms*=2;
        atoms=new Element[natms];
        xyz=new double[3][natms];

        atoms[0]=new Element("Br");
        xyz[0][0]=0.0;
        xyz[1][0]=0.0;
        xyz[2][0]=-3.0*cn1b-ch1b;

        atoms[1]=new Element("H_");
        xyz[0][1]=cgam*cn1b;
        xyz[1][1]=0.0;
        xyz[2][1]=-ch1b-cn1b-cn1b*sgam;

        atoms[2]=new Element("C_3");
        xyz[0][2]=cgam*cn1b;
        xyz[1][2]=0.0;
        xyz[2][2]=-cn1b-cn1b*sgam;

        atoms[3]=new Element("H_");
        xyz[0][3]=cn1b*cgam+ch1b*cgam*0.5;
        xyz[1][3]=ch1b*cgam*0.866025403;
        xyz[2][3]=-cn1b-cn1b*sgam+ch1b*sgam;

        atoms[4]=new Element("H_");
        xyz[0][4]=cn1b*cgam+ch1b*cgam*0.5;
        xyz[1][4]=-ch1b*cgam*0.866025403;
        xyz[2][4]=-cn1b-cn1b*sgam+ch1b*sgam;

        for(int i=1;i<5;i++) {
            atoms[i+4]=new Element(atoms[i].zsym);
            xyz[0][i+4]=-0.5*xyz[0][i]-0.866025403*xyz[1][i];
            xyz[1][i+4]=xyz[0][i]*0.866025403-0.5*xyz[1][i];
            xyz[2][i+4]=xyz[2][i];

            atoms[i+8]=new Element(atoms[i].zsym);
            xyz[0][i+8]=-0.5*xyz[0][i]+0.866025403*xyz[1][i];
            xyz[1][i+8]=-xyz[0][i]*0.866025403-0.5*xyz[1][i];
            xyz[2][i+8]=xyz[2][i];
        }

        atoms[13]=new Element("N_3");
        xyz[0][13]=0.0;
        xyz[1][13]=0.0;
        xyz[2][13]=-cn1b;

        AML.euler(0.0,-beta,0.0,rot);

        for(int i=0;i<14;i++) {
            uuu[0]=xyz[0][i];
            uuu[1]=xyz[1][i];
            uuu[2]=xyz[2][i];
            AML.rotate(op,uuu,vvv,rot);
            xyz[0][i]=vvv[0];
            xyz[1][i]=vvv[1];
            xyz[2][i]=vvv[2];
        }
        return 14;
    }

    int phenol() {
        /*
**********************************************************************

dl_poly/java utility to make the head of a hydrocarbon chain
this one is for p-phenol head group

copyright daresbury laboratory
author w.smith 2011

**********************************************************************
         */
        int natms;
        boolean op=false;
        double f;

        natms=3*ncarbons+13;
        if(twin)natms*=2;
        atoms=new Element[natms];
        xyz=new double[3][natms];

        f=1.03/1.09;

        atoms[0]=new Element("H_");
        xyz[0][0]= oh1b*0.866025403;
        xyz[1][0]=0.0;
        xyz[2][0]=-cc1b-ccab*2.0-coab-oh1b*0.5;

        atoms[1]=new Element("O_R");
        xyz[0][1]=0.0;
        xyz[1][1]=0.0;
        xyz[2][1]=-cc1b-ccab*2.0-coab;

        atoms[2]=new Element("C_R");
        xyz[0][2]=0.0;
        xyz[1][2]=0.0;
        xyz[2][2]=-cc1b-ccab*2.0;

        atoms[3]=new Element("C_R");
        xyz[0][3]=-ccab*0.866025403;
        xyz[1][3]=0.0;
        xyz[2][3]=-cc1b-ccab*1.5;

        atoms[4]=new Element("C_R");
        xyz[0][4]= ccab*0.866025403;
        xyz[1][4]=0.0;
        xyz[2][4]=-cc1b-ccab*1.5;

        atoms[5]=new Element("C_R");
        xyz[0][5]=-ccab*0.866025403;
        xyz[1][5]=0.0;
        xyz[2][5]=-cc1b-ccab*0.5;

        atoms[6]=new Element("C_R");
        xyz[0][6]= ccab*0.866025403;
        xyz[1][6]=0.0;
        xyz[2][6]=-cc1b-ccab*0.5;

        atoms[7]=new Element("C_R");
        xyz[0][7]=0.0;
        xyz[1][7]=0.0;
        xyz[2][7]=-cc1b;

        atoms[8]=new Element("H_");
        xyz[0][8]= (f*ch1b+ccab)*0.866025403;
        xyz[1][8]=0.0;
        xyz[2][8]=-cc1b-ccab*1.5-f*ch1b*0.5;

        atoms[9]=new Element("H_");
        xyz[0][9]=-(f*ch1b+ccab)*0.866025403;
        xyz[1][9]=0.0;
        xyz[2][9]=-cc1b-ccab*1.5-f*ch1b*0.5;

        atoms[10]=new Element("H_");
        xyz[0][10]=-(f*ch1b+ccab)*0.866025403;
        xyz[1][10]=0.0;
        xyz[2][10]=-cc1b+(f*ch1b-ccab)*0.5;

        atoms[11]=new Element("H_");
        xyz[0][11]= (f*ch1b+ccab)*0.866025403;
        xyz[1][11]=0.0;
        xyz[2][11]=-cc1b+(f*ch1b-ccab)*0.5;

        AML.euler(0.0,-beta,0.0,rot);

        for(int i=0;i<12;i++) {
            uuu[0]=xyz[0][i];
            uuu[1]=xyz[1][i];
            uuu[2]=xyz[2][i];
            AML.rotate(op,uuu,vvv,rot);
            xyz[0][i]=vvv[0];
            xyz[1][i]=vvv[1];
            xyz[2][i]=vvv[2];
        }
        return 12;
    }

    int polyoxyeth() {
        /*
*********************************************************************

dl_poly/java utility for generating a polyoxyethylene head group

copyright daresbury laboratory
author    w.smith 2011

*********************************************************************
         */
        int natms;
        double xx0,yy0,zz0;

        natms=3*ncarbons+7*nethos+3;
        if(twin)natms*=2;
        atoms=new Element[natms];
        xyz=new double[3][natms];

        // start position of chain

        xx0=0.0;
        yy0=0.0;
        zz0=0.0;

        natms=7*nethos+3;

        for(int i=0;i<=3*nethos;i++) {
            if(i%3==0) {
                atoms[natms-2]=new Element("O_3");
                xx0=xx0+co1b*xbs[0];
                yy0=yy0+co1b*ybs[0];
                zz0=zz0+co1b*zbs[0];
                xyz[0][natms-2]=xx0;
                xyz[1][natms-2]=yy0;
                xyz[2][natms-2]=zz0;
                natms--;
            }
            else {
                xx0=xx0+cc1b*xbs[0];
                yy0=yy0+cc1b*ybs[0];
                zz0=zz0+cc1b*zbs[0];

                atoms[natms-4]=new Element("C_3");
                xyz[0][natms-4]=xx0;
                xyz[1][natms-4]=yy0;
                xyz[2][natms-4]=zz0;

                atoms[natms-3]=new Element("H_");
                xyz[0][natms-3]=xx0-ch1b*xbs[1];
                xyz[1][natms-3]=yy0+ch1b*ybs[1];
                xyz[2][natms-3]=zz0+ch1b*zbs[1];

                atoms[natms-2]=new Element("H_");
                xyz[0][natms-2]=xx0-ch1b*xbs[2];
                xyz[1][natms-2]=yy0+ch1b*ybs[2];
                xyz[2][natms-2]=zz0+ch1b*zbs[2];

                natms-=3;
            }

            // growth direction vector

            xbs[0]=-xbs[0];
            xbs[1]=-xbs[1];
            xbs[2]=-xbs[2];
            xbs[3]=-xbs[3];

        }

        atoms[0]=new Element("H_");
        xyz[0][0]=xx0+oh1b*xbs[0];
        xyz[1][0]=yy0+oh1b*ybs[0];
        xyz[2][0]=zz0+oh1b*zbs[0];

        natms=7*nethos+2;

        if(nethos%2 == 0) {
            for(int i=0;i<natms;i++)
                xyz[0][i]=-xyz[0][i];
        }

        return natms;
    }

    // Interpret the button and textfield actions

    public void actionPerformed(ActionEvent e) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        int status;
        String arg = (String)e.getActionCommand();
        if (arg.equals("Make")) {
            setTetra();
            ncarbons=BML.giveInteger(nc.getText(),1);
            harea=BML.giveDouble(area.getText(),1);
            nethos=BML.giveInteger(eon.getText(),1);
            zgap=BML.giveDouble(gapz.getText(),1);
            keyhed=head.getSelectedIndex();
            twin=chk1.isSelected();
            flip=chk2.isSelected();
            status=chain();
        }
        else if (arg.equals("Close")) {
            job.dispose();
        }
    }
}
