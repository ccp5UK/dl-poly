import java.io.*;

public class Structure extends Basic {
    /*
*********************************************************************

dl_poly/java class to define the structure of a configuration
in terms of molecules and other components

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
     */

    public int nbnds,nunq,nrept,nshl,ntatm,nmols,nmoltp;
    public int[] lbnd,msz;
    public int[][] join,bond;
    public Molecule[] molecules;
    public String[] name,unqatm;
    public int[] ist,isz;
    private double[][] xyz;
    private Config cfg;
    private int[] idc,mtp,mst;

    Structure() {
        /*
*********************************************************************

dl_poly/java constructor for Structure class

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */

        lbnd=new int[MXATMS];
        bond=new int[MXCONNECT][MXATMS];
        join=new int[2][MXJOIN];

    }

    Structure(Config home) {
        /*
*********************************************************************

dl_poly/java constructor for Structure class

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */

        super();

        cfg=home;

        // identify bonds in system

        nbnds=getBonds();

        // contiguize the structure

        contiguize();

        // identify unique atoms

        nunq=uniqueAtoms();

        // atomic repeat pattern in CONFIG file

        nrept=numRepeatAtoms();

        // number of core-shell units

        nshl=numCoreShells();

        // number of atoms types

        ntatm=numAtomTypes();
        if(nbnds > 0) nmols=molFind();
        if(nmols > 0) nmoltp=molSame();
        if(nmols == 0) nmols=cfg.natms/nrept;

    }

    int getBonds() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

********************g*************************************************
         */
        boolean safe=true;
        int i,j,k,n,mxcon;
        double xd,yd,zd,ff,rsq,radi,radj,bndfac;
        double uuu[]=new double[3];

        lbnd=new int[cfg.natms];
        bond=new int[MXCONNECT][cfg.natms];

        nbnds=0;
        bndfac=1.0+bondpc/100.0;

        for(i=0;i<cfg.natms;i++)
            lbnd[i]=0;

        OUT:
        for(i=0;i<cfg.natms-1;i++) {

	    if(cfg.atoms[i].covalent) {

		radi=cfg.atoms[i].zrad;

		for(j=i+1;j<cfg.natms;j++) {

		    if(cfg.atoms[j].covalent) {
			radj=cfg.atoms[j].zrad;
			ff=bndfac*(radi+radj);
			uuu[0]=cfg.xyz[0][i]-cfg.xyz[0][j];
			uuu[1]=cfg.xyz[1][i]-cfg.xyz[1][j];
			uuu[2]=cfg.xyz[2][i]-cfg.xyz[2][j];
			cfg.pbc.images(uuu);
			rsq=uuu[0]*uuu[0]+uuu[1]*uuu[1]+uuu[2]*uuu[2];
			if(rsq <= ff*ff) {
			    if((lbnd[i] < MXCONNECT)&&(lbnd[j] < MXCONNECT)) {
				bond[lbnd[i]][i]=j;
				bond[lbnd[j]][j]=i;
				lbnd[i]+=1;
				lbnd[j]+=1;
				nbnds++;
			    }
			    else{
				safe=false;
				break OUT;
			    }
			}
		    }
		}
	    }
	}

        if(!safe){
            println("Warning - too many bonds found in Structure");
            nbnds=0;
        }

        // construct list of bonds

        if(nbnds > 0) {
            n=0;
            join=new int[2][nbnds];
            for(i=0;i<cfg.natms;i++) {
                for(j=0;j<lbnd[i];j++) {
                    k=bond[j][i];
                    if(k > i) {
                        join[0][n]=i;
                        join[1][n]=k;
                        n++;
                    }
                }
            }
        }

        return nbnds;

    }

    void contiguize() {
        /*
*********************************************************************

dl_poly/java GUI routine to make a contiguous version of a Config

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */

        Element[] atms;
        double[][] uuu;
        int k,m,n,p,q;
        int[] key,lok,kkk,lll;
        int[][] bbb;

        if(nbnds > 0) {

            key=new int[cfg.natms];
            lok=new int[cfg.natms];

            for (int i=0;i<cfg.natms;i++) {
                lok[i]=i;
                key[i]=i;
            }

            //identify clusters

            for (int i=0;i<nbnds;i++) {

                m=Math.min(join[0][i],join[1][i]);
                n=Math.max(join[0][i],join[1][i]);

                if(key[n] == n) {
                    key[n]=key[m];
                }
                else if (key[n] != key[m]) {

                    p=Math.min(key[m],key[n]);
                    q=Math.max(key[m],key[n]);

                    for (int j=0;j<cfg.natms;j++) {
                        if(key[j] == q) key[j]=p;
                    }
                }
            }

            // sort clusters in ascending order

            AML.ShellSort(cfg.natms,lok,key);

            // construct contiguous configuration

            atms=new Element[cfg.natms];
            uuu=new double[3][cfg.natms];
            for (int i=0;i<cfg.natms;i++) {
                p=lok[i];
                atms[i]=cfg.atoms[p];
                uuu[0][i]=cfg.xyz[0][p];
                uuu[1][i]=cfg.xyz[1][p];
                uuu[2][i]=cfg.xyz[2][p];
            }

            // reconstruct bonding arrays

            join=new int[2][nbnds];
            kkk=new int[cfg.natms];
            lll=new int[cfg.natms];
            bbb=new int[MXCONNECT][cfg.natms];

            for(int i=0;i<cfg.natms;i++)
                kkk[lok[i]]=i;

            n=0;
            for(int i=0;i<cfg.natms;i++) {
                p=lok[i];
                lll[i]=lbnd[p];
                for (int j=0;j<lll[i];j++){
                    k=kkk[bond[j][p]];
                    bbb[j][i]=k;
                    if(k < i){
                        join[0][n]=k;
                        join[1][n]=i;
                        n++;
                    }
                }
            }

            cfg.atoms=atms;
            cfg.xyz=uuu;
            lbnd=lll;
            bond=bbb;
        }
    }

    int uniqueAtoms() {
        /*
*********************************************************************

dl_poly/java class to determine unique atoms in a configuration

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        int kkk,mxunq=MXUNIQUE;
        boolean lnew;

        name=new String[cfg.natms];
        unqatm=new String[mxunq];

        for(int i=0;i<cfg.natms;i++)
            name[i]=BML.fmt(cfg.atoms[i].zsym,8);

        // determine unique atom types

        kkk=1;
        unqatm[0]=name[0];

        for(int i=1;i<cfg.natms;i++) {

            lnew=true;
            for(int j=0;j<kkk;j++) {
                if(name[i].equals(unqatm[j])) lnew=false;
            }
            if(lnew) {
                if(kkk == mxunq) {
                    String unqnew[]=new String[2*mxunq];
                    System.arraycopy(unqatm,0,unqnew,0,mxunq);
                    unqatm=unqnew;
                    mxunq*=2;
                }
                unqatm[kkk]=name[i];
                kkk++;
            }

        }
        return kkk;
    }

    int numRepeatAtoms() {
        /*
*********************************************************************

dl_poly/java routine to determine the repeat pattern of atoms in a
CONFIG file

copyright - daresbury laboratory
author    - w.smith 2011

**********************************************************************
         */
        boolean lnum;
        int nnum;

        nnum=1;
        nrept=cfg.natms;

        OUT1:
            for(int j=1;j<=cfg.natms/2;j++) {

                if(cfg.natms%j == 0) {

                    lnum=true;
                    for(int i=0;i<cfg.natms-j;i++) {
                        if(!name[i].equals(name[i+j]))lnum=false;
                    }
                    if(lnum) {
                        nrept=j;
                        nnum=cfg.natms/j;
                        break OUT1;
                    }
                }
            }
            return nrept;
    }

    int numCoreShells() {
        /*
**********************************************************************

dl_poly/java routine to determine number of core-shell units in a
configuration

copyright - daresbury laboratory
author    - w.smith 2011

**********************************************************************
         */

        nshl=0;
        for(int i=0;i<nrept;i++)
            if(name[i].charAt(4)=='_' && name[i].toLowerCase().charAt(5)=='s')nshl++;

        return nshl;

    }

    int numAtomTypes() {
        /*
**********************************************************************

dl_poly/java routine to determine number of atom types in a
configuration allowing for core-shell units

copyright - daresbury laboratory
author    - w.smith 2011

**********************************************************************
         */

        ntatm=0;
        for(int i=0;i<nunq;i++) {

            if(unqatm[i].substring(4).equals("    ")) {
                ntatm++;
            }
            else if(unqatm[i].toLowerCase().substring(4,6).equals("_s")) {
                ntatm++;
            }

        }
        return ntatm;

    }

    int molFind() {
        /*
*********************************************************************

dl_poly/java routine for identifying individual molecules in a
configuration using cluster analysis

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        int nmols,jcl,kcl,iatm,jatm,idm,nmax;

        idc=new int[cfg.natms];

        // initialise cluster arrays

        for (int i=0;i<cfg.natms;i++)
            idc[i]=i;

        // search for molecules

        nmols=cfg.natms;
        for(int i=0;i<nbnds;i++) {

            iatm=Math.min(join[0][i],join[1][i]);
            jatm=Math.max(join[0][i],join[1][i]);

            if(idc[jatm] == jatm) {
                idc[jatm]=idc[iatm];
                nmols--;
            }
            else if (idc[jatm] != idc[iatm]) {
                jcl=Math.min(idc[iatm],idc[jatm]);
                kcl=Math.max(idc[iatm],idc[jatm]);

                for(int k=0;k<cfg.natms;k++) {
                    if(idc[k] == kcl) idc[k]=jcl;
                }
                nmols--;
            }

        }

        //println("Number of molecules found: "+BML.fmt(nmols,6));

        if(nmols == 0)return 0;

        // define  molecule locations in CONFIG arrays

        ist=new int[nmols];
        isz=new int[nmols];

        for(int i=0;i<nmols;i++) {
            ist[i]=0;
            isz[i]=0;
        }

        idm=0;
        for(int i=1;i<cfg.natms;i++) {

            if(idc[i] != idc[i-1]) {
                isz[idm]=i-ist[idm];
                idm++;
                if(idm == nmols) {
                    println("Error - molecule data not contiguous in CONFIG");
                    return -1;
                }
                ist[idm]=i;
            }

        }
        isz[idm]=cfg.natms-ist[idm];
        idm++;

        nmax=1;
        for(int i=0;i<nmols;i++)
            nmax=Math.max(nmax,isz[i]);

        //println("Largest molecule found: "+BML.fmt(nmax,6));

        return nmols;
    }

    int molSame() {
        /*
*********************************************************************

dl_poly/java routine for identifying sequences of identical molecules
in a configuration (NB this operation valid for dl_poly only)

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        Molecule mole;
        String molnam;
        boolean same,found;
        int ia,ib,kk,mm,idm,kkk,ja,jb;

        if(nmols == 0) return 0;

        mtp=new int[nmols];
        mst=new int[nmols];
        msz=new int[nmols];

        for(int i=0;i<nmols;i++) {
            mtp[i]=i;
            mst[i]=0;
            msz[i]=0;
        }

        // now identify equivalent molecules

        for(int i=1;i<nmols;i++) {

            // check molecule size

            same=true;
            OUT:
                if(isz[i] == isz[i-1]) {
                    ia=ist[i];
                    ib=ist[i-1];

                    // compare corresponding atoms

                    for(int j=0;j<isz[i];j++) {

                        // check atom types

                        if(!(cfg.atoms[ia+j].zsym.equals(cfg.atoms[ib+j].zsym))) {
                            same=false;
                            break OUT;
                        }

                        // check atom valencies

                        if(lbnd[ia+j] != lbnd[ib+j]) {
                            same=false;
                            break OUT;
                        }
                    }

                    // check molecular topology

                    found=false;
                    OUT1:
                        for(int j=0;j<isz[i];j++) {

                            kk=ia+j;
                            mm=ib+j;
                            for(int k=0;k<lbnd[kk];k++) {

                                kkk=bond[k][kk]-ist[i];
                                for(int m=0;m<lbnd[mm];m++) {
                                    if(kkk == bond[m][mm]-ist[i-1]) {
                                        found=true;
                                    }
                                }
                                if(!found) break OUT1;

                            }
                        }
                        same=found;
                        if(!same) break OUT;
                        mtp[i]=mtp[i-1];
                }
        }

        // make list of unique molecules

        idm=0;
        mst[0]=0;
        for(int i=1;i<nmols;i++) {

            if(mtp[i] != mtp[i-1]) {
                msz[idm]=i-mst[idm];
                mst[++idm]=i;
            }

        }
        msz[idm]=nmols-mst[idm];
        nmoltp=++idm;

        // assign structural details to molecules

        molecules=new Molecule[nmoltp];

        for(int i=0;i<nmoltp;i++) {

            ja=ist[mst[i]];
            jb=isz[mst[i]];
            molnam="Species "+BML.fmt(i,6);
            mole=new Molecule(cfg,ja,jb,molnam);
            molecules[i]=mole;

        }

        return nmoltp;
    }

    void resizeJoinArray() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        int[][] jjn=new int[2][join[0].length+MXJOIN];
        for(int i=0;i<nbnds;i++) {
            jjn[0][i]=join[0][i];
            jjn[1][i]=join[1][i];
        }
        join=jjn;
    }

    void resizeBondArrays() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        int[] lll=new int[lbnd.length+MXATMS];
        int[][] bbb=new int[MXCONNECT][bond[0].length+MXATMS];
        for(int i=0;i<config.natms;i++) {
            lll[i]=lbnd[i];
            for(int j=0;j<MXCONNECT;j++)
                bbb[j][i]=bond[j][i];
        }
        for(int i=config.natms;i<lll.length;i++)
            lll[i]=0;
        lbnd=lll;
        bond=bbb;
    }

}
