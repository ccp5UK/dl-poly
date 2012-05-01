import java.io.*;
import java.awt.*;

public class Config extends Basic {
        /*
*********************************************************************

dl_poly/java class to define the contents of a CONFIG file

copyright - daresbury laboratory
author    - w.smith may 2006

*********************************************************************
         */
    public GUI home;
    public String fname;
    public String title;
    public int natms,levcfg;
    public CellBoundary pbc;
    public double[][] xyz;
    public Element[] atoms;
    public Structure structure;

    Config(){
        /*
*********************************************************************

dl_poly/java constructor to define Config class

copyright - daresbury laboratory
author    - w.smith march 2011

*********************************************************************
         */
        super();
        fname=null;
        home=null;
        natms=0;
        levcfg=0;
        title="GUI Reference : "+dateToday();
        structure=null;
        pbc=new CellBoundary();
        xyz=new double[3][MXATMS];
        atoms=new Element[MXATMS];
    }

    boolean configWrite(String filename) {
        /*
*********************************************************************

dl_poly/java routine to write a DL_POLY CONFIG file

copyright - daresbury laboratory
author    - w.smith march 2011

*********************************************************************
         */
        try {
            DataOutputStream outStream = new DataOutputStream(new FileOutputStream(filename));
            outStream.writeBytes(title+"\n");
            outStream.writeBytes(BML.fmt(0,10)+BML.fmt(pbc.imcon,10)+"\n");
            if(pbc.imcon>0) {
                outStream.writeBytes(BML.fmt(pbc.cell[0],20)+BML.fmt(pbc.cell[1],20)+BML.fmt(pbc.cell[2],20)+"\n");
                outStream.writeBytes(BML.fmt(pbc.cell[3],20)+BML.fmt(pbc.cell[4],20)+BML.fmt(pbc.cell[5],20)+"\n");
                outStream.writeBytes(BML.fmt(pbc.cell[6],20)+BML.fmt(pbc.cell[7],20)+BML.fmt(pbc.cell[8],20)+"\n");
            }
            for (int k=0;k<natms;k++) {
                outStream.writeBytes(BML.fmt(atoms[k].zsym,8)+BML.fmt(k+1,10)+BML.fmt(atoms[k].znum,10)+"\n");
                outStream.writeBytes(BML.fmt(xyz[0][k],20)+BML.fmt(xyz[1][k],20)+BML.fmt(xyz[2][k],20)+"\n");
            }
            outStream.close();
            println("New CONFIG file created: "+filename);
        }
        catch(Exception e) {
            println("Error - writing file: "+filename);
            return false;
        }
        return true;
    }

    boolean rdCFG(String fname) {
        /*
*********************************************************************

dl_poly/java routine to read a DL_POLY CONFIG file

copyright - daresbury laboratory
author    - w.smith march 2011

*********************************************************************
         */
        int i,j,k,m;
        LineNumberReader lnr=null;
        String record="",namstr="";
        double xlo=0,xhi=0,ylo=0,yhi=0,zlo=0,zhi=0,ddd;

        // open the CONFIG file

        try {
            lnr = new LineNumberReader(new FileReader(fname));
            println("Reading file: "+fname);
            title = lnr.readLine();
            println("File header record: "+title);
            record = lnr.readLine();
            levcfg=BML.giveInteger(record,1);
            pbc.imcon =BML.giveInteger(record,2);
            if(pbc.imcon > 0) {
                record = lnr.readLine();
                pbc.cell[0]=BML.giveDouble(record,1);
                pbc.cell[1]=BML.giveDouble(record,2);
                pbc.cell[2]=BML.giveDouble(record,3);
                record = lnr.readLine();
                pbc.cell[3]=BML.giveDouble(record,1);
                pbc.cell[4]=BML.giveDouble(record,2);
                pbc.cell[5]=BML.giveDouble(record,3);
                record = lnr.readLine();
                pbc.cell[6]=BML.giveDouble(record,1);
                pbc.cell[7]=BML.giveDouble(record,2);
                pbc.cell[8]=BML.giveDouble(record,3);
            }

            // read coordinates

            i=0;
            j=0;
            natms=0;
            k=levcfg+2;
            while((record=lnr.readLine()) != null) {

                i=j/k;
                m=j-k*i;
                if(m == 0) {
                    namstr = BML.fmt(BML.giveWord(record,1),8);

                    if(natms == atoms.length)
                        resizeArrays();
                    if(!namstr.equals("")) {
                        atoms[i]=new Element(namstr);
                        natms++;
                    }
                    else {
                        println("Error  - unknown atom type in CONFIG file: "+namstr);
                        lnr.close();
                        return false;
                    }
                }
                if(m == 1) {
                    xyz[0][i]=BML.giveDouble(record,1);
                    xyz[1][i]=BML.giveDouble(record,2);
                    xyz[2][i]=BML.giveDouble(record,3);
                    if(pbc.imcon == 0) {
                        xlo=Math.min(xlo,xyz[0][i]);
                        ylo=Math.min(ylo,xyz[1][i]);
                        zlo=Math.min(zlo,xyz[2][i]);
                        xhi=Math.max(xhi,xyz[0][i]);
                        yhi=Math.max(yhi,xyz[1][i]);
                        zhi=Math.max(zhi,xyz[2][i]);
                    }
                }
                j++;
            }
            lnr.close();
        }
        catch(FileNotFoundException e) {
            println("Error - file not found: " + fname);
            return false;
        }
        catch(Exception e) {
            println("Error reading file: " + fname + " "+e);
            return false;
        }

        // construct virtual cell for imcon=0

        if(pbc.imcon == 0) {
            pbc.imcon=2;
            ddd=Math.max(xhi-xlo,yhi-ylo);
            ddd=Math.max(ddd,zhi-zlo)+10;

            pbc.cell[0]=ddd;
            pbc.cell[4]=ddd;
            pbc.cell[8]=ddd;
        }
        pbc.buildBoundary(pbc.imcon);
        if(pbc.imcon > 0)
            pbc.images(natms,xyz);
        structure = new Structure(this);
        println("Selected CONFIG file loaded successfully");
        println("Number of atoms found: "+natms);

        return true;
    }

     boolean rdFRG(String fname) {
        /*
*********************************************************************

dl_poly/java routine to read a DL_POLY FRAGMENT file

copyright - daresbury laboratory
author    - w.smith march 2011

*********************************************************************
         */
        int i,j,k,m;
        LineNumberReader lnr=null;
        String record="",namstr="";

        // open the FRAGMENT file

        try {
            InputStream instream = this.getClass().getResourceAsStream(fname);
            InputStreamReader isr = new InputStreamReader(instream);
            BufferedReader reader = new BufferedReader(isr);
            println("Reading file: "+fname);
            title = reader.readLine();
            println("File header record: "+title);
            record = reader.readLine();
            levcfg=BML.giveInteger(record,1);
            pbc.imcon =BML.giveInteger(record,2);

            // read coordinates

            i=0;
            j=0;
            natms=0;
            k=levcfg+2;
            while((record=reader.readLine()) != null) {

                i=j/k;
                m=j-k*i;
                if(m == 0) {
                    namstr = BML.fmt(BML.giveWord(record,1),8);

                    if(natms == atoms.length)
                        resizeArrays();
                    if(!namstr.equals("")) {
                        atoms[i]=new Element(namstr);
                    natms++;
                    }
                    else {
                        println("Error  - unknown atom type in CONFIG file: "+namstr);
                        reader.close();
                        return false;
                    }
                }
                if(m == 1) {
                    xyz[0][i]=BML.giveDouble(record,1);
                    xyz[1][i]=BML.giveDouble(record,2);
                    xyz[2][i]=BML.giveDouble(record,3);
                }
                j++;
            }
            reader.close();
        }
        catch(FileNotFoundException e) {
            println("Error - file not found: " + fname);
            return false;
        }
        catch(Exception e) {
            println("Error reading file: " + fname + " "+e);
            return false;
        }

        structure = new Structure(this);
        println("Selected FRAGMENT file loaded successfully");
        println("Number of atoms found: "+natms);

        return true;
    }

    boolean rdXYZ(String fname) {
        /*
*********************************************************************

dl_poly/java routine to read an XYZ configuration file

copyright - daresbury laboratory
author    - w.smith march 2011

*********************************************************************
         */
        LineNumberReader lnr=null;
        String record="",namstr="";
        double xlo,xhi,ylo,yhi,zlo,zhi,ddd;

        xlo=0.0;
        ylo=0.0;
        zlo=0.0;
        xhi=0.0;
        yhi=0.0;
        zhi=0.0;
        pbc.imcon=2;
        levcfg=0;

        // open XYZ file

        try {
            lnr = new LineNumberReader(new FileReader(fname));
            println("Reading file: "+fname);
            record = lnr.readLine();
            natms = BML.giveInteger(record,1);
            println("Number of atoms in file : "+natms);
            if(natms > atoms.length) {
                atoms=new Element[natms];
                xyz=new double[3][natms];
            }
            title = lnr.readLine();
            println("File header record: "+title);

            // read coordinates

            for(int i=0;i<natms;i++) {
                record = lnr.readLine();
                namstr = BML.fmt(BML.giveWord(record,1),8);
                if(!namstr.equals(""))
                    atoms[i]=new Element(namstr);
                else {
                    println("Error  - unknown atom type in XYZ file: "+namstr);
                    lnr.close();
                    return false;
                }
                xyz[0][i]=BML.giveDouble(record,2);
                xyz[1][i]=BML.giveDouble(record,3);
                xyz[2][i]=BML.giveDouble(record,4);
                xlo=Math.min(xlo,xyz[0][i]);
                ylo=Math.min(ylo,xyz[1][i]);
                zlo=Math.min(zlo,xyz[2][i]);
                xhi=Math.max(xhi,xyz[0][i]);
                yhi=Math.max(yhi,xyz[1][i]);
                zhi=Math.max(zhi,xyz[2][i]);
            }
            lnr.close();
        }
        catch(FileNotFoundException e) {
            println("Error - file not found: " + fname);
            return false;
        }
        catch(Exception e) {
            println("Error reading file: " + fname + " "+e);
            return false;
        }

        // centre the cell contents

        for(int i=0;i<natms;i++) {
            xyz[0][i]-=(xhi-xlo)/2.0;
            xyz[1][i]-=(yhi-ylo)/2.0;
            xyz[2][i]-=(zhi-zlo)/2.0;
        }

        // build boundary condition

        ddd=Math.pow(((xhi-xlo)*(yhi-ylo)*(zhi-zlo)/natms),(1.0/3.0));
        pbc.cell[0]=-(xlo-xhi)+ddd;
        pbc.cell[4]=-(ylo-yhi)+ddd;
        pbc.cell[8]=-(zlo-zhi)+ddd;
        pbc.buildBoundary(pbc.imcon);
        structure = new Structure(this);
        println("Selected XYZ file loaded successfully");

        return true;
    }

    boolean rdPDB(String fname) {
        /*
*********************************************************************

dl_poly/java routine to read a PDB configuration file

copyright - daresbury laboratory
author    - w.smith may 2011

*********************************************************************
         */
        LineNumberReader lnr=null;
        String record="",namstr="";
        String label,header,atom,hetatm,pdbres,remark;
        double x,y,z,xlo,xhi,ylo,yhi,zlo,zhi,ddd;

        natms=0;
        pbc.imcon=2;
        levcfg=0;
        xlo=0.0;
        ylo=0.0;
        zlo=0.0;
        xhi=0.0;
        yhi=0.0;
        zhi=0.0;
        header="HEADER";
        atom="ATOM";
        hetatm="HETATM";
        pdbres="PDBRES";
	remark="REMARK";

        // open PDB file

        try {
            lnr = new LineNumberReader(new FileReader(fname));
            println("Reading file: "+fname);

            // read configuration

            while((record=lnr.readLine()) != null) {

                if(record.indexOf(header)>=0) {
                    title=record.substring(9,79);
                    println("Header: "+title);
                }
                else if(record.indexOf(remark)>=0) {
                    println("Remark: "+record.substring(9));
                }
                else if(record.indexOf(atom)>=0 || record.indexOf(hetatm)>=0) {
                    if(natms == atoms.length)
                        resizeArrays();
		    label=BML.giveWord(record,1);
                    namstr=BML.fmt(BML.giveWord(record,3).trim(),8);
                    x=BML.giveDouble(record,6);
                    y=BML.giveDouble(record,7);
                    z=BML.giveDouble(record,8);
                    if(label.indexOf(hetatm)>=0 && namstr.equals("O       ")) namstr="OW      ";
                    atoms[natms]=new Element(namstr);
                    xyz[0][natms]=x;
                    xyz[1][natms]=y;
                    xyz[2][natms]=z;
                    if(natms == 0) {
                        xlo=xhi=x;
                        ylo=yhi=y;
                        zlo=zhi=z;
		    }
                    else {
                        xlo=Math.min(xlo,x);
                        ylo=Math.min(ylo,y);
                        zlo=Math.min(zlo,z);
                        xhi=Math.max(xhi,x);
                        yhi=Math.max(yhi,y);
                        zhi=Math.max(zhi,z);
                    }
                    natms++;
                }
            }
            lnr.close();
        }
        catch(FileNotFoundException e) {
            println("Error - file not found: " + fname);
            return false;
        }
        catch(Exception e) {
            println("Error reading file: " + fname + " "+e);
            return false;
        }

        // centre the cell contents

        for(int i=0;i<natms;i++) {
            xyz[0][i]-=(xhi+xlo)/2.0;
            xyz[1][i]-=(yhi+ylo)/2.0;
            xyz[2][i]-=(zhi+zlo)/2.0;
        }

        // build boundary condition

        ddd=Math.max(Math.max(xhi-xlo,yhi-ylo),zhi-zlo);
        ddd=Math.pow((ddd*ddd*ddd/natms),(1.0/3.0));
        pbc.cell[0]=-(xlo-xhi)+ddd;
        pbc.cell[4]=-(ylo-yhi)+ddd;
        pbc.cell[8]=-(zlo-zhi)+ddd;
        pbc.buildBoundary(pbc.imcon);
        structure = new Structure(this);
	println("Number of atoms found: "+natms);
        println("Selected PDB file loaded successfully");

        return true;
    }

    boolean rdMSI(String fname) {
        /*
*********************************************************************

dl_poly/java utility to read a CERIUS 2 configuration file

CERIUS 2 is the copyright of Molecular Simulations Inc

copyright - daresbury laboratory
author    - w.smith  may 2011

*********************************************************************
         */
        LineNumberReader lnr=null;
        String record="",namstr="";
        int keypbc,i,j,k,m;
        double xlo=0,xhi=0,ylo=0,yhi=0,zlo=0,zhi=0,ddd;

        natms=0;
        pbc.imcon=2;
        levcfg=0;

        // open the CERIUS file

        try {
            lnr = new LineNumberReader(new FileReader(fname));
            println("Reading file: "+fname);
            title=lnr.readLine();
            println("File header record: "+title);
            while((record=lnr.readLine()) != null) {
                if(record.indexOf("Model")>=0) {
                    // do nothing in current implementation
                }
                else if((m=record.indexOf("PeriodicType"))>=0) {
                    record=record.substring(m+12);
                    keypbc=BML.giveInteger(record,1);
                    if(keypbc==100)pbc.imcon=1;
                }
                else if((m=record.indexOf("A3"))>=0) {
                    record=record.substring(m+2);
                    pbc.cell[0]=BML.giveDouble(record,1);
                    pbc.cell[1]=BML.giveDouble(record,2);
                    pbc.cell[2]=BML.giveDouble(record,3);
                }
                else if((m=record.indexOf("B3"))>=0) {
                    record=record.substring(m+2);
                    pbc.cell[3]=BML.giveDouble(record,1);
                    pbc.cell[4]=BML.giveDouble(record,2);
                    pbc.cell[5]=BML.giveDouble(record,3);
                }
                else if((m=record.indexOf("C3"))>=0) {
                    record=record.substring(m+2);
                    pbc.cell[6]=BML.giveDouble(record,1);
                    pbc.cell[7]=BML.giveDouble(record,2);
                    pbc.cell[8]=BML.giveDouble(record,3);
                }
                else if(record.indexOf("Atom1")>=0) {
                    //Ignore!
                }
                else if(record.indexOf("Atom2")>=0) {
                    //Ignore!
                }
                else if(record.indexOf("Atom")>=0) {
                    if(natms == atoms.length) {
                        Element atomz[]=new Element[2*atoms.length];
                        double uvw[][]=new double[3][2*atoms.length];
                        for(int n=0;n<atoms.length;n++) {
                            atomz[n]=new Element(atoms[n].zsym);
                            uvw[0][n]=xyz[0][n];
                            uvw[1][n]=xyz[1][n];
                            uvw[2][n]=xyz[2][n];
                        }
                        atoms=atomz;
                        xyz=uvw;
                    }
                    xyz[0][natms]=0.0;
                    xyz[1][natms]=0.0;
                    xyz[2][natms]=0.0;
                    namstr="        ";

                    OUT:
                        for(k=0;k<100;k++) {
                            record=lnr.readLine();
                            record=record.trim();
                            if(record.charAt(0)==')') {
                                natms++;
                                break OUT;
                            }
                            else if((m=record.indexOf("FFType"))>=0) {
                                record=record.substring(m+6);
                                namstr=BML.giveWord(record,1);
                                if(namstr.equals("H___A")) namstr="H__HB";
                                atoms[natms]=new Element(namstr);
                            }
                            else if((m=record.indexOf("Charge"))>=0) {
                                record=record.substring(m+6);
                                atoms[natms].zchg=BML.giveDouble(record,1);
                            }
                            else if((m=record.indexOf("Mass"))>=0) {
                                record=record.substring(m+4);
                                atoms[natms].zmas=BML.giveDouble(record,1);
                            }
                            else if((m=record.indexOf("XYZ"))>=0) {
                                record=record.substring(m+3);
                                xyz[0][natms]=BML.giveDouble(record,1);
                                xyz[1][natms]=BML.giveDouble(record,2);
                                xyz[2][natms]=BML.giveDouble(record,3);
                                if(pbc.imcon == 0) {
                                    xlo=Math.min(xlo,xyz[0][natms]);
                                    ylo=Math.min(ylo,xyz[1][natms]);
                                    zlo=Math.min(zlo,xyz[2][natms]);
                                    xhi=Math.max(xhi,xyz[0][natms]);
                                    yhi=Math.max(yhi,xyz[1][natms]);
                                    zhi=Math.max(zhi,xyz[2][natms]);
                                }
                            }
                        }
                }
            }
            lnr.close();
        }
        catch(FileNotFoundException e) {
            println("Error - file not found: " + fname);
            return false;
        }
        catch(Exception e) {
            println("Error reading file: " + fname + " "+e);
            println(record);
            return false;
        }
        if(pbc.imcon == 0) {

            // centre the cell contents

            for(i=0;i<natms;i++) {
                xyz[0][i]-=(xhi-xlo)/2.0;
                xyz[1][i]-=(yhi-ylo)/2.0;
                xyz[2][i]-=(zhi-zlo)/2.0;
            }

            // build boundary condition

            pbc.imcon=2;
            ddd=Math.pow((xhi-xlo)*(yhi-ylo)*(zhi-zlo)/natms,(1.0/3.0));
            pbc.cell[0]=-(xlo-xhi)+ddd;
            pbc.cell[4]=-(ylo-yhi)+ddd;
            pbc.cell[8]=-(zlo-zhi)+ddd;
        }
        pbc.buildBoundary(pbc.imcon);
        structure = new Structure(this);
        println("CERIUS file loaded successfully");
        println("Number of atoms found: "+natms);

        return true;
    }

    public void resizeArrays(){
        /*
*********************************************************************

dl_poly/java routine to resize the arrays of a DL_POLY CONFIG file

copyright - daresbury laboratory
author    - w.smith march 2011

*********************************************************************
         */

        Element atomz[]=new Element[atoms.length+MXATMS];
        double uvw[][]=new double[3][atoms.length+MXATMS];

        for(int n=0;n<natms;n++) {
            atomz[n]=new Element();
            atomz[n].znum=atoms[n].znum;
            atomz[n].zmas=atoms[n].zmas;
            atomz[n].zchg=atoms[n].zchg;
            atomz[n].zrad=atoms[n].zrad;
            atomz[n].zsym=new String(atoms[n].zsym);
            atomz[n].zcol=new Color(atoms[n].zcol.getRGB());
            atomz[n].covalent=atoms[n].covalent;
            atomz[n].dotify=atoms[n].dotify;
            uvw[0][n]=xyz[0][n];
            uvw[1][n]=xyz[1][n];
            uvw[2][n]=xyz[2][n];
        }
        atoms=atomz;
        xyz=uvw;
    }
}
