import java.awt.*;

public class Element {
    /*
*********************************************************************

dl_poly/java class to determine the element object given
the corresponding atomic symbol

author    - w.smith March 2011

*********************************************************************
     */
    public int znum;
    public String zsym;
    public Color zcol;
    public double zrad,zmas,zchg;
    public boolean covalent;

    public Element(){
        /*
*********************************************************************

dl_poly/java GUI routine: Element class constructor

copyright - daresbury laboratory
author    - w.smith march 2011

*********************************************************************
         */
        znum=0;
        zmas=0.0;
        zchg=0.0;
        zrad=0.0;
        zsym=null;
        zcol=null;
	covalent=false;
    }

    public Element(String atnam) {
        /*
*********************************************************************

dl_poly/java GUI routine: Element class constructor

copyright - daresbury laboratory
author    - w.smith march 2011

*********************************************************************
         */
        znum=0;
        zmas=0.0;
        zchg=0.0;
        zrad=1.0;
        zcol=null;
        zsym=BML.fmt(atnam,8);
	covalent=false;

        if(zsym.charAt(0)=='H') {
            if(zsym.indexOf("H ") == 0) {
                znum=1;
                zmas=1.00797;
                zrad=0.37;
                zcol=Color.white;
		covalent=true;
            }
            else if(zsym.indexOf("H_") == 0) {
                znum=1;
                zmas=1.00797;
                zrad=0.37;
                zcol=Color.white;
		covalent=true;
            }
            else if(zsym.indexOf("HW") == 0) {
                znum=1;
                zmas=1.00797;
                zrad=0.37;
                zcol=Color.white;
		covalent=true;
            }
            else if(zsym.indexOf("He") == 0) {
                znum=2;
                zmas=4.0026;
                zcol=new Color(0,200,0);
            }
            else if(zsym.indexOf("Ho") == 0) {
                znum=67;
                zmas=164.930;
                zcol=Color.lightGray;
            }
            else if(zsym.indexOf("Hf") == 0) {
                znum=72;
                zmas=178.49;
                zrad=1.44;
                zcol=Color.lightGray;
            }
            else if(zsym.indexOf("Hg") == 0) {
                znum=80;
                zmas=200.59;
                zrad=1.44;
                zcol=new Color(216,216,216);
            }
        }
        else if(zsym.charAt(0)=='C') {
            if(zsym.indexOf("C_1") == 0) {
                znum=6;
                zmas=12.01115;
                zcol=Color.gray;
                zrad=0.60;
		covalent=true;
            }
            else if (zsym.indexOf("C_2") == 0) {
                znum=6;
                zmas=12.01115;
                zcol=Color.gray;
                zrad=0.67;
		covalent=true;
            }
            else if (zsym.indexOf("C_R") == 0) {
                znum=6;
                zmas=12.01115;
                zcol=Color.gray;
                zrad=0.70;
		covalent=true;
            }
            else if (zsym.indexOf("C_3") == 0) {
                znum=6;
                zmas=12.01115;
                zcol=Color.gray;
                zrad=0.77;
		covalent=true;
            }
            else if(zsym.indexOf("C ") == 0) {
                znum=6;
                zmas=12.01115;
                zrad=0.77;
                zcol=Color.gray;
		covalent=true;
            }
            else if(zsym.indexOf("C_") == 0) {
                znum=6;
                zmas=12.01115;
                zrad=0.77;
                zcol=Color.gray;
		covalent=true;
            }
            else if(zsym.indexOf("Cl") == 0) {
                znum=17;
                zmas=35.453;
                zrad=0.99;
                zcol=new Color(0,128,64);
		covalent=true;
            }
            else if(zsym.indexOf("Ca") == 0) {
                znum=20;
                zmas=40.08;
                zrad=1.74;
                zcol=Color.lightGray;
            }
            else if(zsym.indexOf("Cr") == 0) {
                znum=24;
                zmas=51.996;
                zrad=1.17;
                zcol=Color.lightGray;
            }
            else if(zsym.indexOf("Co") == 0) {
                znum=27;
                zmas=58.9332;
                zrad=1.16;
                zcol=new Color(255,128,255);
            }
            else if(zsym.indexOf("Cu") == 0) {
                znum=29;
                zmas=63.546;
                zrad=1.17;
                zcol=new Color(106,106,255);
            }
            else if(zsym.indexOf("Cd") == 0) {
                znum=48;
                zmas=112.40;
                zrad=1.41;
                zcol=Color.lightGray;
            }
            else if(zsym.indexOf("Cs") == 0) {
                znum=55;
                zmas=132.905;
                zrad=2.35;
                zcol=Color.lightGray;
            }
            else if(zsym.indexOf("Ce") == 0) {
                znum=58;
                zmas=140.12;
                zcol=Color.lightGray;
            }
            else if(zsym.indexOf("Cm") == 0) {
                znum=96;
                zmas=247.;
                zcol=Color.lightGray;
            }
            else if(zsym.indexOf("Cf") == 0) {
                znum=98;
                zmas=251.;
                zcol=Color.lightGray;
            }
        }
        else if(zsym.charAt(0)=='D') {
            if(zsym.indexOf("D ") == 0) {
                znum=1;
                zmas=2.0;
                zrad=0.37;
                zcol=new Color(225,225,225);
		covalent=true;
            }
            else if(zsym.indexOf("D_") == 0) {
                znum=1;
                zmas=2.0;
                zrad=0.37;
                zcol=new Color(225,225,225);
            }
            else if(zsym.indexOf("Dy") == 0) {
                znum=66;
                zmas=162.50;
                zcol=Color.lightGray;
            }
        }
        else if(zsym.charAt(0)=='O') {
            if (zsym.indexOf("O_2") == 0) {
                znum=8;
                zmas=15.9994;
                zcol=Color.red;
                zrad=0.61;
		covalent=true;
            }
            else if (zsym.indexOf("O_3") == 0) {
                znum=8;
                zmas=15.9994;
                zcol=Color.red;
                zrad=0.66;
		covalent=true;
            }
            else if(zsym.indexOf("O ") == 0) {
                znum=8;
                zmas=15.9994;
                zrad=0.74;
                zcol=Color.red;
		covalent=true;
            }
            else if(zsym.indexOf("O_") == 0) {
                znum=8;
                zmas=15.9994;
                zrad=0.74;
                zcol=Color.red;
		covalent=true;
            }
            else if(zsym.indexOf("OW") == 0) {
                znum=8;
                zmas=15.9994;
                zrad=0.74;
                zcol=Color.red;
		covalent=true;
            }
            else if(zsym.indexOf("Os") == 0) {
                znum=76;
                zmas=190.2;
                zrad=1.26;
                zcol=Color.lightGray;
            }
        }
        else if(zsym.charAt(0)=='N') {
            if (zsym.indexOf("N_1") == 0) {
                znum=7;
                zmas=14.0067;
                zcol=Color.blue;
                zrad=0.55;
		covalent=true;
            }
            else if (zsym.indexOf("N_2") == 0) {
                znum=7;
                zmas=14.0067;
                zcol=Color.blue;
                zrad=0.60;
		covalent=true;
            }
            else if (zsym.indexOf("N_3") == 0) {
                znum=7;
                zmas=14.0067;
                zcol=Color.blue;
                zrad=0.73;
		covalent=true;
            }
            else if(zsym.indexOf("N ") == 0) {
                znum=7;
                zmas=14.0067;
                zrad=0.74;
                zcol=Color.blue;
		covalent=true;
            }
            else if(zsym.indexOf("N_") == 0) {
                znum=7;
                zmas=14.0067;
                zrad=0.74;
                zcol=Color.blue;
		covalent=true;
            }
            else if(zsym.indexOf("Ne") == 0) {
                znum=10;
                zmas=20.179;
                zcol=new Color(0,200,0);
            }
            else if(zsym.indexOf("Na") == 0) {
                znum=11;
                zmas=22.9898;
                zrad=1.57;
                zcol=new Color(255,255,128);
            }
            else if(zsym.indexOf("Ni") == 0) {
                znum=28;
                zmas=58.71;
                zrad=1.15;
                zcol=new Color(0,198,99);
            }
            else if(zsym.indexOf("Nb") == 0) {
                znum=41;
                zmas=92.906;
                zrad=1.34;
                zcol=Color.lightGray;
            }
            else if(zsym.indexOf("Nd") == 0) {
                znum=60;
                zmas=144.24;
                zcol=Color.lightGray;
            }
            else if(zsym.indexOf("Np") == 0) {
                znum=93;
                zmas=237.;
                zcol=Color.lightGray;
            }
            else if(zsym.indexOf("No") == 0) {
                znum=102;
                zmas=254.;
                zcol=Color.lightGray;
            }
        }
        else if(zsym.charAt(0)=='P') {
            if (zsym.indexOf("P_2") == 0) {
                znum=15;
                zmas=30.9738;
                zcol=Color.orange;
                zrad=0.77;
		covalent=true;
            }
            else if (zsym.indexOf("P_3") == 0) {
                znum=15;
                zmas=30.9738;
                zrad=1.11;
                zcol=Color.orange;
		covalent=true;
            }
            else if(zsym.indexOf("P ") == 0) {
                znum=15;
                zmas=30.9738;
                zrad=1.10;
                zcol=Color.orange;
		covalent=true;
            }
            else if(zsym.indexOf("P_") == 0) {
                znum=15;
                zmas=30.9738;
                zrad=1.10;
                zcol=Color.orange;
		covalent=true;
            }
            else if(zsym.indexOf("Pd") == 0) {
                znum=46;
                zmas=106.4;
                zrad=1.28;
                zcol=Color.lightGray;
            }
            else if(zsym.indexOf("Pr") == 0) {
                znum=59;
                zmas=140.907;
                zcol=Color.lightGray;
            }
            else if(zsym.indexOf("Pm") == 0) {
                znum=61;
                zmas=147.;
                zcol=Color.lightGray;
            }
            else if(zsym.indexOf("Pt") == 0) {
                znum=78;
                zmas=195.09;
                zrad=1.29;
                zcol=Color.lightGray;
            }
            else if(zsym.indexOf("Pb") == 0) {
                znum=82;
                zmas=207.19;
                zrad=1.54;
                zcol=Color.lightGray;
            }
            else if(zsym.indexOf("Po") == 0) {
                znum=84;
                zmas=210.;
                zrad=1.46;
                zcol=Color.lightGray;
            }
            else if(zsym.indexOf("Pa") == 0) {
                znum=91;
                zmas=231.;
                zcol=Color.lightGray;
            }
            else if(zsym.indexOf("Pu") == 0) {
                znum=94;
                zmas=242.;
                zcol=Color.lightGray;
            }
        }
        else if(zsym.charAt(0)=='S') {
            if(zsym.indexOf("Si") == 0) {
                znum=14;
                zmas=28.086;
                zrad=1.17;
                zcol=new Color(128,64,0);
		covalent=true;
            }
            else if (zsym.indexOf("S_2") == 0) {
                znum=16;
                zmas=32.064;
                zcol=Color.yellow;
                zrad=0.94;
		covalent=true;
            }
            else if (zsym.indexOf("S_3") == 0) {
                znum=16;
                zmas=32.064;
                zcol=Color.yellow;
                zrad=1.04;
		covalent=true;
            }
            else if(zsym.indexOf("S ") == 0) {
                znum=16;
                zmas=32.064;
                zrad=1.04;
                zcol=new Color(255,255,0);
		covalent=true;
            }
            else if(zsym.indexOf("S_") == 0) {
                znum=16;
                zmas=32.064;
                zrad=1.04;
                zcol=new Color(255,255,0);
		covalent=true;
            }
            else if(zsym.indexOf("Sc") == 0) {
                znum=21;
                zmas=44.956;
                zrad=1.44;
                zcol=Color.lightGray;
            }
            else if(zsym.indexOf("Se") == 0) {
                znum=34;
                zmas=78.96;
                zrad=1.17;
                zcol=Color.lightGray;
            }
            else if(zsym.indexOf("Sr") == 0) {
                znum=38;
                zmas=87.62;
                zrad=1.91;
                zcol=Color.lightGray;
            }
            else if(zsym.indexOf("Sn") == 0) {
                znum=50;
                zmas=118.69;
                zrad=1.40;
                zcol=Color.lightGray;
            }
            else if(zsym.indexOf("Sb") == 0) {
                znum=51;
                zmas=121.75;
                zrad=1.41;
                zcol=Color.lightGray;
            }
            else if(zsym.indexOf("Sm") == 0) {
                znum=62;
                zmas=150.35;
                zcol=Color.lightGray;
            }
        }
        else if(zsym.charAt(0)=='B') {
            if(zsym.indexOf("Be") == 0) {
                znum=4;
                zmas=9.0122;
                zrad=0.89;
                zcol=Color.lightGray;
            }
            else if(zsym.indexOf("B ") == 0) {
                znum=5;
                zmas=10.811;
                zrad=0.80;
                zcol=new Color(0,0,225);
		covalent=true;
            }
            else if(zsym.indexOf("B_") == 0) {
                znum=5;
                zmas=10.811;
                zrad=0.80;
                zcol=new Color(0,0,225);
		covalent=true;
            }
            else if(zsym.indexOf("Br") == 0) {
                znum=35;
                zmas=79.904;
                zrad=1.14;
                zcol=new Color(145,0,72);
		covalent=true;
            }
            else if(zsym.indexOf("Ba") == 0) {
                znum=56;
                zmas=137.34;
                zrad=1.98;
                zcol=Color.lightGray;
            }
            else if(zsym.indexOf("Bi") == 0) {
                znum=83;
                zmas=208.980;
                zrad=1.52;
                zcol=Color.lightGray;
            }
            else if(zsym.indexOf("Bk") == 0) {
                znum=97;
                zmas=247.;
                zcol=Color.lightGray;
            }
        }
        else if(zsym.charAt(0)=='K') {
            if(zsym.indexOf("K ") == 0) {
                znum=19;
                zmas=39.102;
                zrad=2.03;
                zcol=new Color(255,0,255);
            }
            else if(zsym.indexOf("K_") == 0) {
                znum=19;
                zmas=39.102;
                zrad=2.03;
                zcol=new Color(255,0,255);
            }
            else if(zsym.indexOf("Kr") == 0) {
                znum=36;
                zmas=83.80;
                zcol=new Color(0,200,0);
            }
        }
        else if(zsym.charAt(0)=='F') {
            if(zsym.indexOf("F ") == 0) {
                znum=9;
                zmas=18.9984;
                zrad=0.72;
                zcol=new Color(255,128,0);
		covalent=true;
            }
            else if(zsym.indexOf("F_") == 0) {
                znum=9;
                zmas=18.9984;
                zrad=0.72;
                zcol=new Color(255,128,0);
		covalent=true;
            }
            else if(zsym.indexOf("Fe") == 0) {
                znum=26;
                zmas=55.847;
                zrad=1.16;
                zcol=new Color(0,255,64);
            }
            else if(zsym.indexOf("Fr") == 0) {
                znum=87;
                zmas=223.;
                zrad=2.50;
                zcol=Color.lightGray;
            }
            else if(zsym.indexOf("Fm") == 0) {
                znum=100;
                zmas=253.;
                zcol=Color.lightGray;
            }
        }
        else if(zsym.charAt(0)=='M') {
            if(zsym.indexOf("Mg") == 0) {
                znum=12;
                zmas=24.304;
                zrad=1.36;
                zcol=new Color(128,128,64);
            }
            else if(zsym.indexOf("Mn") == 0) {
                znum=25;
                zmas=54.9380;
                zrad=1.17;
                zcol=Color.lightGray;
            }
            else if(zsym.indexOf("Mo") == 0) {
                znum=42;
                zmas=95.94;
                zrad=1.29;
                zcol=Color.lightGray;
            }
            else if(zsym.indexOf("Md") == 0) {
                znum=101;
                zmas=256.;
                zcol=Color.lightGray;
            }
        }
        else if(zsym.charAt(0)=='A') {
            if(zsym.indexOf("Al") == 0) {
                znum=13;
                zmas=26.9815;
                zrad=1.25;
                zcol=Color.lightGray;
            }
            else if(zsym.indexOf("Ar") == 0) {
                znum=18;
                zmas=39.948;
                zcol=new Color(0,200,0);
            }
            else if(zsym.indexOf("As") == 0) {
                znum=33;
                zmas=74.9216;
                zrad=1.21;
                zcol=Color.lightGray;
            }
            else if(zsym.indexOf("Ag") == 0) {
                znum=47;
                zmas=107.868;
                zrad=1.34;
                zcol=new Color(216,216,216);
            }
            else if(zsym.indexOf("Au") == 0) {
                znum=79;
                zmas=196.967;
                zrad=1.34;
                zcol=new Color(232,226,32);
            }
            else if(zsym.indexOf("At") == 0) {
                znum=85;
                zmas=210.;
                zrad=1.40;
                zcol=Color.lightGray;
            }
            else if(zsym.indexOf("Ac") == 0) {
                znum=89;
                zmas=227.;
                zcol=Color.lightGray;
            }
            else if(zsym.indexOf("Am") == 0) {
                znum=95;
                zmas=243.;
                zcol=Color.lightGray;
            }
        }
        else if(zsym.charAt(0)=='L') {
            if(zsym.indexOf("Li") == 0) {
                znum=3;
                zmas=6.939;
                zrad=1.23;
                zcol=new Color(255,0,50);
            }
            else if(zsym.indexOf("La") == 0) {
                znum=57;
                zmas=138.91;
                zrad=1.69;
                zcol=Color.lightGray;
            }
            else if(zsym.indexOf("Lu") == 0) {
                znum=71;
                zmas=174.97;
                zcol=Color.lightGray;
            }
            else if(zsym.indexOf("Lr") == 0) {
                znum=103;
                zmas=257.;
                zcol=Color.lightGray;
            }
        }
        else if(zsym.charAt(0)=='T') {
            if(zsym.indexOf("Ti") == 0) {
                znum=22;
                zmas=47.90;
                zrad=1.32;
                zcol=Color.lightGray;
            }
            else if(zsym.indexOf("Tc") == 0) {
                znum=43;
                zmas=99.;
                zrad=1.27;
                zcol=Color.lightGray;
            }
            else if(zsym.indexOf("Te") == 0) {
                znum=52;
                zmas=127.60;
                zrad=1.37;
                zcol=Color.lightGray;
            }
            else if(zsym.indexOf("Tb") == 0) {
                znum=65;
                zmas=158.924;
                zcol=Color.lightGray;
            }
            else if(zsym.indexOf("Tm") == 0) {
                znum=69;
                zmas=168.934;
                zcol=Color.lightGray;
            }
            else if(zsym.indexOf("Ta") == 0) {
                znum=73;
                zmas=180.948;
                zrad=1.34;
                zcol=Color.lightGray;
            }
            else if(zsym.indexOf("Tl") == 0) {
                znum=81;
                zmas=204.37;
                zrad=1.55;
                zcol=Color.lightGray;
            }
            else if(zsym.indexOf("Th") == 0) {
                znum=90;
                zmas=232.038;
                zcol=Color.lightGray;
            }
        }
        else if(zsym.charAt(0)=='I') {
            if(zsym.indexOf("I ") == 0) {
                znum=53;
                zmas=126.9044;
                zrad=1.33;
                zcol=new Color(128,0,64);
		covalent=true;
            }
            else if(zsym.indexOf("I_") == 0) {
                znum=53;
                zmas=126.9044;
                zrad=1.33;
                zcol=new Color(128,0,64);
		covalent=true;
            }
            else if(zsym.indexOf("In") == 0) {
                znum=49;
                zmas=114.82;
                zrad=1.50;
                zcol=Color.lightGray;
            }
            else if(zsym.indexOf("Ir") == 0) {
                znum=77;
                zmas=192.2;
                zrad=1.26;
                zcol=Color.lightGray;
            }
        }
        else if(zsym.charAt(0)=='R') {
            if(zsym.indexOf("Rb") == 0) {
                znum=37;
                zmas=85.47;
                zrad=2.16;
                zcol=new Color(193,193,193);
            }
            else if(zsym.indexOf("Ru") == 0) {
                znum=44;
                zmas=101.07;
                zrad=1.24;
                zcol=Color.lightGray;
            }
            else if(zsym.indexOf("Rh") == 0) {
                znum=45;
                zmas=102.905;
                zrad=1.25;
                zcol=Color.lightGray;
            }
            else if(zsym.indexOf("Re") == 0) {
                znum=75;
                zmas=186.2;
                zrad=1.28;
                zcol=Color.lightGray;
            }
            else if(zsym.indexOf("Rn") == 0) {
                znum=86;
                zmas=222.;
                zcol=Color.lightGray;
            }
            else if(zsym.indexOf("Ra") == 0) {
                znum=88;
                zmas=226.;
                zcol=Color.lightGray;
            }
        }
        else if(zsym.charAt(0)=='E') {
            if(zsym.indexOf("Eu") == 0) {
                znum=63;
                zmas=151.96;
                zcol=Color.lightGray;
            }
            else if(zsym.indexOf("Er") == 0) {
                znum=68;
                zmas=167.26;
                zcol=Color.lightGray;
            }
            else if(zsym.indexOf("Es") == 0) {
                znum=99;
                zmas=254.;
                zcol=Color.lightGray;
            }
        }
        else if(zsym.charAt(0)=='Z') {
            if(zsym.indexOf("Zn") == 0) {
                znum=30;
                zmas=65.97;
                zrad=1.25;
                zcol=new Color(233,233,233);
            }
            else if(zsym.indexOf("Zr") == 0) {
                znum=40;
                zmas=91.22;
                zrad=1.45;
                zcol=Color.lightGray;
            }
        }
        else if(zsym.charAt(0)=='G') {
            if(zsym.indexOf("Ga") == 0) {
                znum=31;
                zmas=69.72;
                zrad=1.25;
                zcol=Color.lightGray;
            }
            else if(zsym.indexOf("Ge") == 0) {
                znum=32;
                zmas=72.59;
                zrad=1.22;
                zcol=Color.lightGray;
            }
            else if(zsym.indexOf("Gd") == 0) {
                znum=64;
                zmas=157.25;
                zcol=Color.lightGray;
            }
        }
        else if(zsym.charAt(0)=='Y') {
            if(zsym.indexOf("Y ") == 0) {
                znum=39;
                zmas=88.905;
                zrad=1.62;
                zcol=Color.lightGray;
            }
            else if(zsym.indexOf("Y_") == 0) {
                znum=39;
                zmas=88.905;
                zrad=1.62;
                zcol=Color.lightGray;
            }
            else if(zsym.indexOf("Yb") == 0) {
                znum=70;
                zmas=173.04;
                zcol=Color.lightGray;
            }
        }
        else if(zsym.charAt(0)=='W') {
            if(zsym.indexOf("W ") == 0) {
                znum=74;
                zmas=183.85;
                zrad=1.30;
                zcol=Color.lightGray;
            }
            else if(zsym.indexOf("W_") == 0) {
                znum=74;
                zmas=183.85;
                zrad=1.30;
                zcol=Color.lightGray;
            }
        }
        else if(zsym.charAt(0)=='V') {
            if(zsym.indexOf("V ") == 0) {
                znum=23;
                zmas=50.942;
                zrad=1.22;
                zcol=Color.lightGray;
            }
            else if(zsym.indexOf("V_") == 0) {
                znum=23;
                zmas=50.942;
                zrad=1.22;
                zcol=Color.lightGray;
            }
        }
        else if(zsym.charAt(0)=='U') {
            if(zsym.indexOf("U ") == 0) {
                znum=92;
                zmas=238.03;
                zcol=Color.lightGray;
            }
            else if(zsym.indexOf("U_") == 0) {
                znum=92;
                zmas=238.03;
                zcol=Color.lightGray;
            }
        }
        else if(zsym.indexOf("Xe") == 0) {
            znum=54;
            zmas=131.30;
            zcol=new Color(0,200,0);
        }
        else if(zsym.indexOf("QW") == 0) {
            znum=999;
            zmas=0.0;
            zrad=0.2;
            zcol=Color.pink;
	    covalent=true;
        }
    }
}
