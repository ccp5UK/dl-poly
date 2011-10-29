import java.io.*;

public abstract class AML {
    /*
*********************************************************************
     
main dl_poly/java GUI application methods library
     
copyright - daresbury laboratory
author    - w.smith 2000
     
*********************************************************************
     */    
    static double[] invert(double a[]) {
        /*
***********************************************************************
         
dl_poly subroutine to invert a 3 * 3 matrix using cofactors
         
copyright - daresbury laboratory
author    - w. smith december 2000
         
***********************************************************************
         */
        double b[]=new double[9];
        double d,r;
        
        // calculate adjoint matrix
        
        b[0]=a[4]*a[8]-a[5]*a[7];
        b[1]=a[2]*a[7]-a[1]*a[8];
        b[2]=a[1]*a[5]-a[2]*a[4];
        b[3]=a[5]*a[6]-a[3]*a[8];
        b[4]=a[0]*a[8]-a[2]*a[6];
        b[5]=a[2]*a[3]-a[0]*a[5];
        b[6]=a[3]*a[7]-a[4]*a[6];
        b[7]=a[1]*a[6]-a[0]*a[7];
        b[8]=a[0]*a[4]-a[1]*a[3];
        
        // calculate determinant
        
        d=a[0]*b[0]+a[3]*b[1]+a[6]*b[2];
        r=0.0;
        if(Math.abs(d) > 0.0)r=1.0/d;
        
        // complete inverse matrix
        
        b[0]=r*b[0];
        b[1]=r*b[1];
        b[2]=r*b[2];
        b[3]=r*b[3];
        b[4]=r*b[4];
        b[5]=r*b[5];
        b[6]=r*b[6];
        b[7]=r*b[7];
        b[8]=r*b[8];
        
        return b;
    }
    
    static int atmnum(String atnam) {
        /*
*********************************************************************
         
dl_poly/java routine to determine the atomic number given
the corresponding atomic symbol
         
copyright - daresbury laboratory
author    - w.smith october 2000
         
*********************************************************************
         */
        int atnum=0;
        if(atnam.indexOf("Ac") == 0) atnum=89;
        else if(atnam.indexOf("Ag") == 0) atnum=47;
        else if(atnam.indexOf("Al") == 0) atnum=13;
        else if(atnam.indexOf("Am") == 0) atnum=95;
        else if(atnam.indexOf("Ar") == 0) atnum=18;
        else if(atnam.indexOf("As") == 0) atnum=33;
        else if(atnam.indexOf("At") == 0) atnum=85;
        else if(atnam.indexOf("Au") == 0) atnum=79;
        else if(atnam.indexOf("Ba") == 0) atnum=56;
        else if(atnam.indexOf("Be") == 0) atnum=4;
        else if(atnam.indexOf("Bi") == 0) atnum=83;
        else if(atnam.indexOf("Bk") == 0) atnum=97;
        else if(atnam.indexOf("Br") == 0) atnum=35;
        else if(atnam.indexOf("B ") == 0) atnum=5;
        else if(atnam.indexOf("B_") == 0) atnum=5;
        else if(atnam.indexOf("Ca") == 0) atnum=20;
        else if(atnam.indexOf("Cd") == 0) atnum=48;
        else if(atnam.indexOf("Ce") == 0) atnum=58;
        else if(atnam.indexOf("Cf") == 0) atnum=98;
        else if(atnam.indexOf("Cl") == 0) atnum=17;
        else if(atnam.indexOf("Cm") == 0) atnum=96;
        else if(atnam.indexOf("Co") == 0) atnum=27;
        else if(atnam.indexOf("Cr") == 0) atnum=24;
        else if(atnam.indexOf("Cs") == 0) atnum=55;
        else if(atnam.indexOf("Cu") == 0) atnum=29;
        else if(atnam.indexOf("C ") == 0) atnum=6;
        else if(atnam.indexOf("C_") == 0) atnum=6;
        else if(atnam.indexOf("Dy") == 0) atnum=66;
        else if(atnam.indexOf("D ") == 0) atnum=1;
        else if(atnam.indexOf("D_") == 0) atnum=1;
        else if(atnam.indexOf("Er") == 0) atnum=68;
        else if(atnam.indexOf("Es") == 0) atnum=99;
        else if(atnam.indexOf("Eu") == 0) atnum=63;
        else if(atnam.indexOf("Fe") == 0) atnum=26;
        else if(atnam.indexOf("Fm") == 0) atnum=100;
        else if(atnam.indexOf("Fr") == 0) atnum=87;
        else if(atnam.indexOf("F ") == 0) atnum=9;
        else if(atnam.indexOf("F_") == 0) atnum=9;
        else if(atnam.indexOf("Ga") == 0) atnum=31;
        else if(atnam.indexOf("Gd") == 0) atnum=64;
        else if(atnam.indexOf("Ge") == 0) atnum=32;
        else if(atnam.indexOf("He") == 0) atnum=2;
        else if(atnam.indexOf("Hf") == 0) atnum=72;
        else if(atnam.indexOf("Hg") == 0) atnum=80;
        else if(atnam.indexOf("Ho") == 0) atnum=67;
        else if(atnam.indexOf("H ") == 0) atnum=1;
        else if(atnam.indexOf("H_") == 0) atnum=1;
        else if(atnam.indexOf("HW") == 0) atnum=1;
        else if(atnam.indexOf("In") == 0) atnum=49;
        else if(atnam.indexOf("Ir") == 0) atnum=77;
        else if(atnam.indexOf("I ") == 0) atnum=53;
        else if(atnam.indexOf("I_") == 0) atnum=53;
        else if(atnam.indexOf("Kr") == 0) atnum=36;
        else if(atnam.indexOf("K ") == 0) atnum=19;
        else if(atnam.indexOf("K_") == 0) atnum=19;
        else if(atnam.indexOf("La") == 0) atnum=57;
        else if(atnam.indexOf("Li") == 0) atnum=3;
        else if(atnam.indexOf("Lr") == 0) atnum=103;
        else if(atnam.indexOf("Lu") == 0) atnum=71;
        else if(atnam.indexOf("Md") == 0) atnum=101;
        else if(atnam.indexOf("Mg") == 0) atnum=12;
        else if(atnam.indexOf("Mn") == 0) atnum=25;
        else if(atnam.indexOf("Mo") == 0) atnum=42;
        else if(atnam.indexOf("Na") == 0) atnum=11;
        else if(atnam.indexOf("Nb") == 0) atnum=41;
        else if(atnam.indexOf("Nd") == 0) atnum=60;
        else if(atnam.indexOf("Ne") == 0) atnum=10;
        else if(atnam.indexOf("Ni") == 0) atnum=28;
        else if(atnam.indexOf("No") == 0) atnum=102;
        else if(atnam.indexOf("Np") == 0) atnum=93;
        else if(atnam.indexOf("N ") == 0) atnum=7;
        else if(atnam.indexOf("N_") == 0) atnum=7;
        else if(atnam.indexOf("Os") == 0) atnum=76;
        else if(atnam.indexOf("O ") == 0) atnum=8;
        else if(atnam.indexOf("O_") == 0) atnum=8;
        else if(atnam.indexOf("OW") == 0) atnum=8;
        else if(atnam.indexOf("Pa") == 0) atnum=91;
        else if(atnam.indexOf("Pb") == 0) atnum=82;
        else if(atnam.indexOf("Pd") == 0) atnum=46;
        else if(atnam.indexOf("Pm") == 0) atnum=61;
        else if(atnam.indexOf("Po") == 0) atnum=84;
        else if(atnam.indexOf("Pr") == 0) atnum=59;
        else if(atnam.indexOf("Pt") == 0) atnum=78;
        else if(atnam.indexOf("Pu") == 0) atnum=94;
        else if(atnam.indexOf("P ") == 0) atnum=15;
        else if(atnam.indexOf("P_") == 0) atnum=15;
        else if(atnam.indexOf("Ra") == 0) atnum=88;
        else if(atnam.indexOf("Rb") == 0) atnum=37;
        else if(atnam.indexOf("Re") == 0) atnum=75;
        else if(atnam.indexOf("Rh") == 0) atnum=45;
        else if(atnam.indexOf("Rn") == 0) atnum=86;
        else if(atnam.indexOf("Ru") == 0) atnum=44;
        else if(atnam.indexOf("Sb") == 0) atnum=51;
        else if(atnam.indexOf("Sc") == 0) atnum=21;
        else if(atnam.indexOf("Se") == 0) atnum=34;
        else if(atnam.indexOf("Si") == 0) atnum=14;
        else if(atnam.indexOf("Sm") == 0) atnum=62;
        else if(atnam.indexOf("Sn") == 0) atnum=50;
        else if(atnam.indexOf("Sr") == 0) atnum=38;
        else if(atnam.indexOf("S ") == 0) atnum=16;
        else if(atnam.indexOf("S_") == 0) atnum=16;
        else if(atnam.indexOf("Ta") == 0) atnum=73;
        else if(atnam.indexOf("Tb") == 0) atnum=65;
        else if(atnam.indexOf("Tc") == 0) atnum=43;
        else if(atnam.indexOf("Te") == 0) atnum=52;
        else if(atnam.indexOf("Th") == 0) atnum=90;
        else if(atnam.indexOf("Ti") == 0) atnum=22;
        else if(atnam.indexOf("Tl") == 0) atnum=81;
        else if(atnam.indexOf("Tm") == 0) atnum=69;
        else if(atnam.indexOf("U ") == 0) atnum=92;
        else if(atnam.indexOf("U_") == 0) atnum=92;
        else if(atnam.indexOf("V ") == 0) atnum=23;
        else if(atnam.indexOf("V_") == 0) atnum=23;
        else if(atnam.indexOf("W ") == 0) atnum=74;
        else if(atnam.indexOf("W_") == 0) atnum=74;
        else if(atnam.indexOf("Xe") == 0) atnum=54;
        else if(atnam.indexOf("Yb") == 0) atnum=70;
        else if(atnam.indexOf("Y ") == 0) atnum=39;
        else if(atnam.indexOf("Y_") == 0) atnum=39;
        else if(atnam.indexOf("Zn") == 0) atnum=30;
        else if(atnam.indexOf("Zr") == 0) atnum=40;
        return atnum;
    }
    
    static double[] dcell(double aaa[]){
        /*
***********************************************************************
         
dl_poly/java method to calculate the dimensional properties
of a simulation cell specified by the input matrix aaa.
the results are returned in the array bbb, with :
         
bbb(0 to 2) - lengths of cell vectors
bbb(3 to 5) - cosines of cell angles
bbb(6 to 8) - perpendicular cell widths
bbb(9)     - cell volume
         
copyright - daresbury laboratory
author    - w. smith october 2000
         
***********************************************************************
         */
        double axb1,axb2,axb3,bxc1,bxc2,bxc3,cxa1,cxa2,cxa3;
        double bbb[]=new double[10];
	
        // calculate lengths of cell vectors
        
        bbb[0]=Math.sqrt(aaa[0]*aaa[0]+aaa[1]*aaa[1]+aaa[2]*aaa[2]);
        bbb[1]=Math.sqrt(aaa[3]*aaa[3]+aaa[4]*aaa[4]+aaa[5]*aaa[5]);
        bbb[2]=Math.sqrt(aaa[6]*aaa[6]+aaa[7]*aaa[7]+aaa[8]*aaa[8]);
        
        // calculate cosines of cell angles
        
        bbb[3]=(aaa[0]*aaa[3]+aaa[1]*aaa[4]+aaa[2]*aaa[5])/(bbb[0]*bbb[1]);
        bbb[4]=(aaa[0]*aaa[6]+aaa[1]*aaa[7]+aaa[2]*aaa[8])/(bbb[0]*bbb[2]);
        bbb[5]=(aaa[3]*aaa[6]+aaa[4]*aaa[7]+aaa[5]*aaa[8])/(bbb[1]*bbb[2]);
        
        // calculate vector products of cell vectors
        
        axb1=aaa[1]*aaa[5]-aaa[2]*aaa[4];
        axb2=aaa[2]*aaa[3]-aaa[0]*aaa[5];
        axb3=aaa[0]*aaa[4]-aaa[1]*aaa[3];
        bxc1=aaa[4]*aaa[8]-aaa[5]*aaa[7];
        bxc2=aaa[5]*aaa[6]-aaa[3]*aaa[8];
        bxc3=aaa[3]*aaa[7]-aaa[4]*aaa[6];
        cxa1=aaa[7]*aaa[2]-aaa[1]*aaa[8];
        cxa2=aaa[0]*aaa[8]-aaa[2]*aaa[6];
        cxa3=aaa[1]*aaa[6]-aaa[0]*aaa[7];
        
        // calculate volume of cell
        
        bbb[9]=Math.abs(aaa[0]*bxc1+aaa[1]*bxc2+aaa[2]*bxc3);
        
        // calculate cell perpendicular widths
        
        bbb[6]=bbb[9]/Math.sqrt(bxc1*bxc1+bxc2*bxc2+bxc3*bxc3);
        bbb[7]=bbb[9]/Math.sqrt(cxa1*cxa1+cxa2*cxa2+cxa3*cxa3);
        bbb[8]=bbb[9]/Math.sqrt(axb1*axb1+axb2*axb2+axb3*axb3);
	
	return bbb;
    }
    
    static void rotate(boolean op,double uuu[],double vvv[],double rot[]) {
        /*
*********************************************************************
         
dl_poly routine to transform a vector from the local (molecule-
fixed) frame of reference (uuu[0],[1],[2]) to the global frame of
reference (vvv[0],[1],[2]) and vice-versa using the rotation matrix
defined by the matrix equation uuu = (Rot) vvv.
         
copyright - daresbury laboratory
author    - w.smith november 2000
         
note: if op=true  converts global to local
if op=false converts local to global
         
*********************************************************************
         */
        if(op) {
            uuu[0]=vvv[0]*rot[0]+vvv[1]*rot[3]+vvv[2]*rot[6];
            uuu[1]=vvv[0]*rot[1]+vvv[1]*rot[4]+vvv[2]*rot[7];
            uuu[2]=vvv[0]*rot[2]+vvv[1]*rot[5]+vvv[2]*rot[8];
        }
        else {
            vvv[0]=uuu[0]*rot[0]+uuu[1]*rot[1]+uuu[2]*rot[2];
            vvv[1]=uuu[0]*rot[3]+uuu[1]*rot[4]+uuu[2]*rot[5];
            vvv[2]=uuu[0]*rot[6]+uuu[1]*rot[7]+uuu[2]*rot[8];
        }
    }
    
    static void euler(double alp,double bet,double gam,double rot[]) {
        /*
*********************************************************************
         
dl_poly routine to construct a rotation matrix as defined by the
Z-(-)Y'-Z'' euler angle convention and the equation r'=Rr, where
r' is a vector in a molecular frame, r the corresponding vector
in the laboratory frame and R is the rotation matrix
         
copyright - daresbury laboratory
author    - w. smith november 2000
         
*********************************************************************
         */
        double ca,cb,cg,sa,sb,sg;
        
        ca=Math.cos(alp);
        cb=Math.cos(bet);
        cg=Math.cos(gam);
        sa=Math.sin(alp);
        sb=Math.sin(bet);
        sg=Math.sin(gam);
        
        rot[0]= ca*cb*cg-sa*sg;
        rot[1]=-ca*cb*sg-sa*cg;
        rot[2]= ca*sb;
        rot[3]= sa*cb*cg+ca*sg;
        rot[4]=-sa*cb*sg+ca*cg;
        rot[5]= sa*sb;
        rot[6]=-sb*cg;
        rot[7]= sb*sg;
        rot[8]= cb;
    }
    
    static void spline(int npnts,double xx[],double yy[],double zz[],
    double aa[],double dd[],double gg[]) {
        /*
***********************************************************************
         
dl_poly/java utility for spline fit of discrete function
         
copyright - daresbury laboratory
author    - w.smith january 2001
         
***********************************************************************
         */
        
        int n1=npnts-1;
        int n2=npnts-2;
        
        gg[0]=0.0;
        dd[0]=xx[1]-xx[0];
        for(int i=1;i<n1;i++) {
            dd[i]=xx[i+1]-xx[i];
            gg[i]=2.0*(xx[i+1]-xx[i-1]);
            zz[i]=6.0*((yy[i+1]-yy[i])/dd[i]-(yy[i]-yy[i-1])/dd[i-1]);
        }
        gg[n1]=0.0;
        dd[n1]=0.0;
        aa[0]=0.0;
        aa[1]=dd[1]/gg[1];
        for(int i=2;i<n2;i++) {
            gg[i]=gg[i]-dd[i-1]*aa[i-1];
            aa[i]=dd[i]/gg[i];
        }
        gg[n1-1]=gg[n1-1]-dd[n2-1]*aa[n2-1];
        gg[1]=zz[1]/gg[1];
        for(int i=2;i<n1;i++) {
            gg[i]=(zz[i]-dd[i-1]*gg[i-1])/gg[i];
        }
        for(int i=1;i<n2;i++) {
            gg[n1-i]=gg[n1-i]-aa[n1-i]*gg[npnts-i];
        }
    }
    
    static void whatAtoms(String fname) {
        /*
*********************************************************************
         
dl_poly/java method to determine what atom types are present
in a CONFIG or REVCON file
         
copyright - daresbury laboratory
author    - w.smith february 2001
         
*********************************************************************
         */
        int nunq=0,nua=10;
        int levcfg,imcon;
        String record,namstr;
        String[] name;
        int[] num;
        boolean lunq;
        
        num=new int[nua];
        name=new String[nua];
        
        // open the CONFIG file
        
        try {
            LineNumberReader lnr = new LineNumberReader(new FileReader(fname));
            Basic.println("Reading file: "+fname);
            record = lnr.readLine();
            Basic.println("File header record: "+record);
            record = lnr.readLine();
            levcfg=BML.giveInteger(record,1);
            imcon =BML.giveInteger(record,2);
            if(imcon > 0) {
                record = lnr.readLine();
                record = lnr.readLine();
                record = lnr.readLine();
            }
            
            // read configuration
            
            int i=0;
            int k=levcfg+2;
            while((record=lnr.readLine()) != null) {
                if(i%k == 0) {
                    namstr = BML.fmt(BML.giveWord(record,1),8);
                    lunq=true;
                    for(int j=0;j<nunq;j++) {
                        if(namstr.equals(name[j])) {
                            lunq=false;
                            num[j]++;
                        }
                    }
                    
                    if(lunq) {
                        if(nunq==nua) {
                            nua*=2;
                            int nun[]=new int[nua];
                            String ttt[]=new String[nua];
                            System.arraycopy(name,0,ttt,0,nua/2);
                            System.arraycopy(num,0,nun,0,nua/2);
                            name=ttt;
                            num=nun;
                        }
                        name[nunq]=namstr;
                        num[nunq]=1;
                        nunq++;
                    }
                }
                i++;
            }
            lnr.close();
        }
        catch(FileNotFoundException e) {
            Basic.println("Error - file not found: " + fname);
        }
        catch(Exception e) {
            Basic.println("Error reading file: " + fname + " "+e);
        }
        if(nunq>0) {
            Basic.println("Unique atom types in system and populations:");
            int tot=0;
            for(int i=0;i<nunq;i++) {
                Basic.println(BML.fmt(name[i],8)+BML.fmt(num[i],6));
                tot+=num[i];
            }
            Basic.println("Total number of atoms:"+BML.fmt(tot,8));
        }
    }
    
    static void whatAtoms(Config cfg) {
        /*
*********************************************************************
         
dl_poly/java method to determine what atom types are present
in a configuration loaded into the GUI
         
copyright - daresbury laboratory
author    - w.smith august 2011
         
*********************************************************************
         */
        int nunq=0,nua=10;
        String[] name;
        int[] num;
        boolean lunq;
        
        num=new int[nua];
        name=new String[nua];
        
        // scan the configuration
        
	for(int i=0;i<cfg.natms;i++) {
	    
	    lunq=true;
	    for(int j=0;j<nunq;j++) {
		if(cfg.atoms[i].zsym.equals(name[j])) {
		    lunq=false;
		    num[j]++;
		}
	    }
            
	    if(lunq) {
		if(nunq==nua) {
		    nua*=2;
		    int nun[]=new int[nua];
		    String ttt[]=new String[nua];
		    System.arraycopy(name,0,ttt,0,nua/2);
		    System.arraycopy(num,0,nun,0,nua/2);
		    name=ttt;
		    num=nun;
		}
		name[nunq]=cfg.atoms[i].zsym;
		num[nunq]=1;
		nunq++;
	    }
	}

        if(nunq>0) {
            Basic.println("Unique atom types in system and populations:");
            int tot=0;
            for(int i=0;i<nunq;i++) {
                Basic.println(BML.fmt(name[i],8)+BML.fmt(num[i],6));
                tot+=num[i];
            }
            Basic.println("Total number of atoms:"+BML.fmt(tot,8));
        }
    }
    
    static int Jacobi(int n,double a[][],double v[][]) {
        /*
*********************************************************************
         
dl_poly/java GUI routine for diagonalisation of real symmetric
matices by Jacobi method
         
copyright - daresbury laboratory
author    - w.smith 2000
         
*********************************************************************
         */
        int i,j,k,status;
        double tes,scl,rho,tem,v1,v2,v3,omg,u,c,s;
        boolean pas;
        
        tes=0;
        scl=0;
        status=0;
        rho=1e-8;
        
        // rescale matrix for optimal accuracy
        
        for (i=0;i<n;i++) {
            if (Math.abs(a[i][i])>scl) scl=Math.abs(a[i][i]);
        }
        
        for (i=0;i<n;i++) {
            for (j=0;j<=i;j++) {
                a[i][j]=a[i][j]/scl;
            }
        }
        
        // Set initial value of moving tolerance
        
        for (i=1;i<n;i++) {
            for (j=0;j<i;j++) {
                tes=tes+2*a[i][j]*a[i][j];
            }
        }
        
        tes=Math.sqrt(tes);
        
        // Jacobi diagonalisation
        
        while (tes>rho) {
            tes=tes/n;
            if (tes<rho) tes=rho;
            
            pas=true;
            
            while (pas) {
                pas=false;
                
                for (i=1;i<n;i++) {
                    for (j=0;j<i;j++)
                        
                        if (Math.abs(a[i][j])>=tes) {
                            pas=true;
                            v1=a[j][j];
                            v2=a[i][j];
                            v3=a[i][i];
                            u=.5*(v1-v3);
                            if (Math.abs(u)<rho)
                                omg=-1;
                            else
                                omg=-v2/Math.sqrt(v2*v2+u*u);
                            if (u<0)
                                omg=-omg;
                            s=omg/Math.sqrt(2*(1+Math.sqrt(1-omg*omg)));
                            c=Math.sqrt(1-s*s);
                            for (k=0;k<n;k++) {
                                if (k>=i) {
                                    tem=a[k][j]*c-a[k][i]*s;
                                    a[k][i]=a[k][j]*s+a[k][i]*c;
                                    a[k][j]=tem;
                                }
                                else if (k<j) {
                                    tem=a[j][k]*c-a[i][k]*s;
                                    a[i][k]=a[j][k]*s+a[i][k]*c;
                                    a[j][k]=tem;
                                }
                                else {
                                    tem=a[k][j]*c-a[i][k]*s;
                                    a[i][k]=a[k][j]*s+a[i][k]*c;
                                    a[k][j]=tem;
                                }
                                tem=v[k][j]*c-v[k][i]*s;
                                v[k][i]=v[k][j]*s+v[k][i]*c;
                                v[k][j]=tem;
                            }
                            a[j][j]=v1*c*c+v3*s*s-2*v2*s*c;
                            a[i][i]=v1*s*s+v3*c*c+2*v2*s*c;
                            a[i][j]=(v1-v3)*s*c+v2*(c*c-s*s);
                        }
                }
            }
        }
        
        
        //	rescale matrix
        
        for (i=0;i<n;i++) {
            for (j=0;j<=i;j++) {
                a[i][j]=scl*a[i][j];
            }
        }
        return status;
    }
    
    static void ShellSort(int n,int lst[], double aaa[]) {
        /*
**********************************************************************
         
dl_poly/java routine to sort list of real numbers into ascending order
using the shell method
         
copyright - daresbury laboratory
author    - w.smith november 2000
         
**********************************************************************
         */
        int i,j,k,l,m;
        double tmp;
        
        m=n;
        while (m > 0) {
            m=m/2;
            for (j=0;j<n-m;j++) {
                i=j;
                while (i>=0) {
                    l=i+m;
                    if (aaa[l]<aaa[i]) {
                        k=lst[i];
                        lst[i]=lst[l];
                        lst[l]=k;
                        tmp=aaa[i];
                        aaa[i]=aaa[l];
                        aaa[l]=tmp;
                        i=i-m;
                    }
                    else
                        i=-1;
                }
            }
        }
    }
    
    static int gaussfit(int npnts,double ccc[],double eee[],double x[],
    double y[],double z[],double g1[],double g2[],double g3[]) {
        /*
***********************************************************************
         
routine to fit a lennard jones potential with a sum of gaussians
author k. singer
         
java version january 2001 author w.smith
         
***********************************************************************
         */
        
        double varm,d1,d2,d3,dr1,dr2,dr3,a11,a12,a22,a21,a31,a32,a33,a13;
        double a23,b1,b2,b3,r,rr,det,absd,var;
        double c1=0.0,c2=0.0,c3=0.0;
        int k0,k,l,m,i;
        double[] ddd;
        
        k0=15;
        varm=1.0e+25;
        ddd=new double[3];
        
        // set initial guesses for exponents
        
        ddd[0]=0.443;
        ddd[1]=1.544;
        ddd[2]=8.502;
        eee[0]=0.003889;
        eee[1]=0.003889;
        eee[2]=0.007779;
        
        if(npnts<6) {
            eee[0]=0.0;
            ddd[0]=0.0;
            k0=-1;
        }
        else if(npnts<4) {
            Basic.println("Error - too few data points for gaussian fit");
            return -1;
        }
        
        // adjustment of first exponent
        
        for(k=0;k<k0;k++) {
            d1=ddd[0]+(k-7)*eee[0];
            
            // adjustment of second exponent
            
            for(l=0;l<15;l++) {
                d2=ddd[1]+(l-7)*eee[1];
                
                // adjustment of third exponent
                
                for(m=0;m<15;m++) {
                    
                    d3=ddd[2]+(m-7)*eee[2];
                    
                    // store current exponents
                    
                    dr1=d1;
                    dr2=d2;
                    dr3=d3;
                    
                    // initialise least squares parameters
                    
                    a11=0.0;
                    a12=0.0;
                    a22=0.0;
                    a21=0.0;
                    a31=0.0;
                    a32=0.0;
                    a33=0.0;
                    a13=0.0;
                    a23=0.0;
                    b1=0.0;
                    b2=0.0;
                    b3=0.0;
                    
                    // calculate least squares parameters
                    
                    if(k0<0) {
                        for(i=0;i<npnts;i++) {
                            r=x[i];
                            rr=r*r;
                            g2[i]=Math.exp(-rr*dr2);
                            g3[i]=Math.exp(-rr*dr3);
                            a22=a22+Math.pow(g2[i],2);
                            a23=a23+g2[i]*g3[i];
                            a33=a33+Math.pow(g3[i],2);
                            b2=b2+g2[i]*y[i];
                            b3=b3+g3[i]*y[i];
                        }
                        det=(a22*a33-a23*a23);
                    }
                    else {
                        for(i=0;i<npnts;i++) {
                            r=x[i];
                            rr=r*r;
                            g1[i]=Math.exp(-rr*dr1);
                            g2[i]=Math.exp(-rr*dr2);
                            g3[i]=Math.exp(-rr*dr3);
                            a11=a11+Math.pow(g1[i],2);;
                            a12=a12+g1[i]*g2[i];
                            a22=a22+Math.pow(g2[i],2);;
                            a13=a13+g1[i]*g3[i];
                            a23=a23+g2[i]*g3[i];
                            a33=a33+Math.pow(g3[i],2);
                            b1=b1+g1[i]*y[i];
                            b2=b2+g2[i]*y[i];
                            b3=b3+g3[i]*y[i];
                        }
                        det=a11*(a22*a33-a23*a23)+a12*(a23*a13-a12*a33)+a13*(a12*a23-a22*a13);
                    }
                    absd=Math.abs(det);
                    if(absd>2.e-20) {
                        // calculate coefficients
                        
                        if(k0<0) {
                            c2=b2*a33-b3*a23;
                            c3=b3*a22-b2*a23;
                        }
                        else {
                            c1=b1*(a22*a33-a23*a23)+b2*(a13*a23-a12*a33)+b3*(a12*a23-a13*a22);
                            c2=b1*(a23*a13-a12*a33)+b2*(a11*a33-a13*a13)+b3*(a12*a13-a11*a23);
                            c3=b1*(a12*a23-a22*a13)+b2*(a12*a13-a11*a23)+b3*(a11*a22-a12*a12);
                        }
                        c1=c1/det;
                        c2=c2/det;
                        c3=c3/det;
                        
                        // construct approximating potential array
                        
                        var=0.0;
                        if(k0<0) {
                            for(i=0;i<npnts;i++) {
                                z[i]=c2*g2[i]+c3*g3[i];
                                var=var+x[i]*Math.pow((y[i]-z[i]),2);
                            }
                        }
                        else {
                            for(i=0;i<npnts;i++) {
                                z[i]=c1*g1[i]+c2*g2[i]+c3*g3[i];
                                var=var+x[i]*Math.pow((y[i]-z[i]),2);
                            }
                        }
                        
                        // update potential parameters
                        
                        if(var<=varm) {
                            varm=var;
                            ccc[0]=c1;
                            ccc[1]=c2;
                            ccc[2]=c3;
                            eee[0]=d1;
                            eee[1]=d2;
                            eee[2]=d3;
                        }
                    }
                    ddd[2]=eee[2];
                }
                ddd[1]=eee[1];
            }
            ddd[0]=eee[0];
        }
        
        // print out best parameters
        
        Basic.println("Best coefficients: "+BML.fmt(ccc[0],12)+BML.fmt(ccc[1],12)+BML.fmt(ccc[2],12));
        Basic.println("Best exponents   : "+BML.fmt(eee[0],12)+BML.fmt(eee[1],12)+BML.fmt(eee[2],12));
        return 0;
    }
    
    static void ShellSort(int n,int lst[], int iii[]) {
        /*
*********************************************************************
         
dl_poly/java routine to sort list of integers into ascending order
using the shell method
         
copyright - daresbury laboratory
author    - w.smith 2000
         
*********************************************************************
         */
        int i,j,k,l,m;
        
        m=n;
        while (m > 0) {
            m=m/2;
            for (j=0;j<n-m;j++) {
                i=j;
                while (i>=0) {
                    l=i+m;
                    if (iii[lst[l]]<iii[lst[i]]) {
                        k=lst[i];
                        lst[i]=lst[l];
                        lst[l]=k;
                        i=i-m;
                    }
                    else
                        i=-1;
                }
            }
        }
    }
    static double ArcTan2(double x,double y) {
        /*
*********************************************************************
         
stars routine
author - w.smith 2002
         
calculate the arctangent with quadrant correction
         
*********************************************************************
*/
        double t,ax,ay;
        ax=Math.abs(x);
        ay=Math.abs(y);
        
        if(ax+ay < 2e-10) {
            t=0;
        }
        else if( ax < 1e-10) {
            if(y > 0)
                t=0.5*Math.PI;
            else
                t=1.5*Math.PI;
        }
        else {
            t=Math.atan(y/x);
            if(x < 0)
                t+=Math.PI;
            else if(x>=0 && y<0)
                t+=2*Math.PI;
        }
        return t;
    }
}
