import java.io.*;

public class Molecule extends Basic {
        /*
*********************************************************************

dl_poly/java class to define a molecule

copyright - daresbury laboratory
author    - w.smith may 2006

*********************************************************************
         */
    public int natm,nbnd,nang,ndih,ninv;
    public String molname;
    public Element[] atom;
    public double[] chge;
    public int[] list,mid;
    public int[][] ibnd,iang,idih,invr,link;
    public double[][] xyz;
    private Config cfg;

    Molecule(Config home,int ibeg,int itot,String molnam) {
        /*
*********************************************************************

dl_poly/java constructor to define Molecule object from a CONFIG list

copyright - daresbury laboratory
author    - w.smith june 2001

*********************************************************************
         */

	super();
	cfg=home;
        int i=0;
        nbnd=0;
        nang=0;
        ndih=0;
        ninv=0;
        natm=itot;
        molname=molnam;
	chge=new double[natm];
        atom=new Element[natm];
        xyz=new double[3][natm];
        for(int j=ibeg;j<ibeg+itot;j++) {
	    chge[i]=0;
            atom[i]=cfg.atoms[j];
            xyz[0][i]=cfg.xyz[0][j];
            xyz[1][i]=cfg.xyz[1][j];
            xyz[2][i]=cfg.xyz[2][j];
            i++;
        }
        bonds();
        angles();
    }

    void bonds() {
        /*
**********************************************************************

dl_poly/java utility to assign bonds in a molecule
based on a bond distance criterion

copyright daresbury laboratory
author w. smith  september 2001

**********************************************************************
         */

        int last,m1,m2,j,kkk,kbnd,ncl;
        double det,ff,xd,yd,zd,ssx,ssy,ssz,xss,yss,zss,rsq,bndfac;

        // initial values

        ncl=12;
        xd=0.0;
        yd=0.0;
        zd=0.0;
        ssx=0.0;
        ssy=0.0;
        ssz=0.0;
        xss=0.0;
        yss=0.0;
        zss=0.0;
        rsq=0.0;
        kbnd=natm;
        bndfac=1.0+bondpc/100.0;
        ibnd=new int[2][kbnd];
        link=new int[ncl][natm];

        // zero list array

        list=new int[natm];
        for(int i=0;i<natm;i++)
            list[i]=0;

        // determine bonds in molecule

        nbnd=0;
        last=natm;
        m1=natm/2;
        m2=(natm-1)/2;

        for(int m=1;m<=m1;m++) {
            if(m>m2)last=m1;

            for(int i=0;i<last;i++) {
                j=i+m;
                if(j>=natm)j=j-natm;
                ff=bndfac*(atom[i].zrad+atom[j].zrad);
                xd=xyz[0][i]-xyz[0][j];
                yd=xyz[1][i]-xyz[1][j];
                zd=xyz[2][i]-xyz[2][j];
                if(cfg.pbc.imcon>0) {
                    ssx=(cfg.pbc.rcell[0]*xd+cfg.pbc.rcell[3]*yd+cfg.pbc.rcell[6]*zd);
                    ssy=(cfg.pbc.rcell[1]*xd+cfg.pbc.rcell[4]*yd+cfg.pbc.rcell[7]*zd);
                    ssz=(cfg.pbc.rcell[2]*xd+cfg.pbc.rcell[5]*yd+cfg.pbc.rcell[8]*zd);
                    xss=ssx-BML.nint(ssx);
                    yss=ssy-BML.nint(ssy);
                    zss=ssz-BML.nint(ssz);
                    xd=(cfg.pbc.cell[0]*xss+cfg.pbc.cell[3]*yss+cfg.pbc.cell[6]*zss);
                    yd=(cfg.pbc.cell[1]*xss+cfg.pbc.cell[4]*yss+cfg.pbc.cell[7]*zss);
                    zd=(cfg.pbc.cell[2]*xss+cfg.pbc.cell[5]*yss+cfg.pbc.cell[8]*zss);
                }
                if((Math.abs(xd)<=ff) && (Math.abs(yd)<=ff) && (Math.abs(zd)<=ff)) {
                    rsq=xd*xd+yd*yd+zd*zd;
                    if(rsq<=ff*ff) {
                        if(nbnd==kbnd) {
                            int jbnd[][]=new int[2][2*kbnd];
                            for(int n=0;n<kbnd;n++) {
                                jbnd[0][n]=ibnd[0][n];
                                jbnd[1][n]=ibnd[1][n];
                            }
                            ibnd=jbnd;
                            kbnd*=2;
                        }
                        ibnd[0][nbnd]=Math.min(i,j);
                        ibnd[1][nbnd]=Math.max(i,j);
                        nbnd++;
                        kkk=Math.max(list[i],list[j]);
                        if(kkk==ncl) {
                            int kink[][]=new int[2*ncl][natm];
                            for (int n=0;n<natm;n++) {
                                for(int k=0;k<ncl;k++) {
                                    kink[k][n]=link[k][n];
                                }
                            }
                            link=kink;
                            ncl*=2;
                        }
                        link[list[i]][i]=j;
                        link[list[j]][j]=i;
                        list[i]++;
                        list[j]++;
                    }
                }
            }
        }
        //println("Number of possible pair bonds found: "+BML.fmt(nbnd,6));
    }

    void angles() {
        /*
**********************************************************************

dl_poly/java utility to assign bonds, valence angles,
dihedrals and inversion angles

copyright daresbury laboratory
author w. smith  sepetember 2001

**********************************************************************
         */

        int kang,kdih,kinv,jj,kk;

        // determine valence angles in system

        nang=0;
        kang=natm;
        iang=new int[3][kang];

        for(int i=0;i<natm;i++) {
            for(int j=1;j<list[i];j++) {
                for(int k=0;k<j;k++) {
                    if(nang==kang) {
                        int jang[][]=new int[3][2*kang];
                        for(int n=0;n<kang;n++) {
                            jang[0][n]=iang[0][n];
                            jang[1][n]=iang[1][n];
                            jang[2][n]=iang[2][n];
                        }
                        iang=jang;
                        kang*=2;
                    }
                    iang[0][nang]=Math.min(link[j][i],link[k][i]);
                    iang[1][nang]=i;
                    iang[2][nang]=Math.max(link[j][i],link[k][i]);
                    nang++;
                }
            }
        }
        //println("Number of possible valence angles found: "+BML.fmt(nang,6));

        // determine dihedral angles in system

        ndih=0;
        kdih=natm;
        mid=new int[kdih];
        idih=new int[4][kdih];

        for(int i=0;i<nbnd;i++) {
            jj=ibnd[0][i];
            kk=ibnd[1][i];
            for(int j=0;j<list[jj];j++) {
                if(kk!=link[j][jj]) {
                    for(int k=0;k<list[kk];k++) {
                        if(jj!=link[k][kk]) {
                            if(ndih==kdih) {
                                int jid[]=new int[2*kdih];
                                int jdih[][]=new int[4][2*kdih];
                                for(int n=0;n<kdih;n++) {
                                    jid[n]=mid[n];
                                    jdih[0][n]=idih[0][n];
                                    jdih[1][n]=idih[1][n];
                                    jdih[2][n]=idih[2][n];
                                    jdih[3][n]=idih[3][n];
                                }
                                idih=jdih;
                                mid=jid;
                                kdih*=2;
                            }
			    if(link[j][jj] != link[k][kk])
				{
				    mid[ndih]=i;
				    idih[0][ndih]=link[j][jj];
				    idih[1][ndih]=jj;
				    idih[2][ndih]=kk;
				    idih[3][ndih]=link[k][kk];
				    ndih++;
				}
                        }
                    }
                }
            }
        }
        //println("Number of possible dihedral angles found: "+BML.fmt(ndih,6));

        // determine inversion angles in system

        ninv=0;
        kinv=natm;
        invr=new int[4][kinv];
        for(int i=0;i<natm;i++) {
            if(list[i]>2) {
                for(int j=0;j<list[i]-2;j++) {
                    for(int k=j+1;k<list[i]-1;k++) {
                        for(int m=k+1;m<list[i];m++) {
                            if(ninv==kinv) {
                                int jinv[][]=new int[4][2*kinv];
                                for(int n=0;n<kinv;n++) {
                                    jinv[0][n]=invr[0][n];
                                    jinv[1][n]=invr[1][n];
                                    jinv[2][n]=invr[2][n];
                                    jinv[3][n]=invr[3][n];
                                }
                                invr=jinv;
                                kinv*=2;
                            }
                            invr[0][ninv]=i;
                            invr[1][ninv]=link[j][i];
                            invr[2][ninv]=link[k][i];
                            invr[3][ninv]=link[m][i];
                            ninv++;
                        }
                    }
                }
            }
        }
        //println("Number of possible inversion angles found: "+BML.fmt(ninv,6));
    }
}
