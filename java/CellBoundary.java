import java.io.*;
import java.awt.*;

public class CellBoundary{
        /*
*********************************************************************
         
dl_poly/java class to define simulation cell boundary conditions
         
copyright - daresbury laboratory
author    - w.smith March 2011
         
*********************************************************************
         */
    public String kind;
    public int imcon;
    public int num_edges;
    public int num_vertices;
    public int[][] edge;
    public double[][] vrt;
    public double[] cell;
    public double[] rcell;
    
    public CellBoundary() {
        /*
*********************************************************************
         
dl_poly/java base constructor to define cell boundary
         
copyright - daresbury laboratory
author    - w.smith march 2011
         
*********************************************************************
         */
        kind="Unspecified";
        imcon=0;
        num_edges=0;
        num_vertices=0;
        cell=new double[9];
        for(int i=0;i<9;i++)
            cell[i]=0.0;
    }
    
    public void buildBoundary(int imcon){
        /*
*********************************************************************
         
dl_poly/java routine to build cell boundaries
         
copyright - daresbury laboratory
author    - w.smith may 2006
         
*********************************************************************
         */

        rcell=AML.invert(cell);
        
        if(imcon == 0) {
            kind="Open Boundary";
        }
        else if(imcon == 1) {
            kind="Cubic Boundary";
            boxCELL();
        }
        else if(imcon == 2) {
            kind="Orthorhombic Boundary";
            boxCELL();
        }
        else if(imcon == 3) {
            kind="Triclinic Boundary";
            boxCELL();
        }
        else if(imcon == 4) {
            kind="Truncated Octahedron";
            boxOCT();
        }
        else if(imcon == 5) {
            kind="Rhombic Dodecahedron";
            boxDEC();
        }
        else if(imcon == 7) {
            kind="Hexagonal Prism";
            boxHEX();
        }
    }
    
    public void boxCELL() {
        /*
*********************************************************************
         
dl_poly/java GUI routine - standard periodic boundaries
         
copyright - daresbury laboratory
author    - w.smith 2001
         
*********************************************************************
         */
        num_vertices=8;
        num_edges=12;
        
        vrt=new double[3][num_vertices];
        edge=new int[2][num_edges];
        
        vrt[0][0]=-0.5*(cell[0]+cell[3]+cell[6]);
        vrt[1][0]=-0.5*(cell[1]+cell[4]+cell[7]);
        vrt[2][0]=-0.5*(cell[2]+cell[5]+cell[8]);
        
        vrt[0][1]= 0.5*(cell[0]-cell[3]-cell[6]);
        vrt[1][1]= 0.5*(cell[1]-cell[4]-cell[7]);
        vrt[2][1]= 0.5*(cell[2]-cell[5]-cell[8]);
        
        vrt[0][2]= 0.5*(cell[0]+cell[3]-cell[6]);
        vrt[1][2]= 0.5*(cell[1]+cell[4]-cell[7]);
        vrt[2][2]= 0.5*(cell[2]+cell[5]-cell[8]);
        
        vrt[0][3]=-0.5*(cell[0]-cell[3]+cell[6]);
        vrt[1][3]=-0.5*(cell[1]-cell[4]+cell[7]);
        vrt[2][3]=-0.5*(cell[2]-cell[5]+cell[8]);
        
        vrt[0][4]=-0.5*(cell[0]+cell[3]-cell[6]);
        vrt[1][4]=-0.5*(cell[1]+cell[4]-cell[7]);
        vrt[2][4]=-0.5*(cell[2]+cell[5]-cell[8]);
        
        vrt[0][5]= 0.5*(cell[0]-cell[3]+cell[6]);
        vrt[1][5]= 0.5*(cell[1]-cell[4]+cell[7]);
        vrt[2][5]= 0.5*(cell[2]-cell[5]+cell[8]);
        
        vrt[0][6]= 0.5*(cell[0]+cell[3]+cell[6]);
        vrt[1][6]= 0.5*(cell[1]+cell[4]+cell[7]);
        vrt[2][6]= 0.5*(cell[2]+cell[5]+cell[8]);
        
        vrt[0][7]=-0.5*(cell[0]-cell[3]-cell[6]);
        vrt[1][7]=-0.5*(cell[1]-cell[4]-cell[7]);
        vrt[2][7]=-0.5*(cell[2]-cell[5]-cell[8]);
        
        edge[0][0]=0;
        edge[1][0]=1;
        edge[0][1]=1;
        edge[1][1]=2;
        edge[0][2]=2;
        edge[1][2]=3;
        edge[0][3]=3;
        edge[1][3]=0;
        edge[0][4]=4;
        edge[1][4]=5;
        edge[0][5]=5;
        edge[1][5]=6;
        edge[0][6]=6;
        edge[1][6]=7;
        edge[0][7]=7;
        edge[1][7]=4;
        edge[0][8]=0;
        edge[1][8]=4;
        edge[0][9]=1;
        edge[1][9]=5;
        edge[0][10]=2;
        edge[1][10]=6;
        edge[0][11]=3;
        edge[1][11]=7;
    }
    
    public void boxHEX() {
        /*
*********************************************************************
         
dl_poly/java GUI routine - hexagonal periodic boundaries
         
copyright - daresbury laboratory
author    - w.smith 2001
         
*********************************************************************
         */
        num_vertices=12;
        num_edges=18;
        double aaa=cell[0]/3.0;
        double bbb=0.5*aaa*Math.sqrt(3.0);
        double ccc=0.5*cell[8];
        
        vrt=new double[3][num_vertices];
        edge=new int[2][num_edges];
        
        vrt[0][0]= aaa;
        vrt[1][0]= 0.0;
        vrt[2][0]= ccc;
        
        vrt[0][1]= 0.5*aaa;
        vrt[1][1]= bbb;
        vrt[2][1]= ccc;
        
        vrt[0][2]=-0.5*aaa;
        vrt[1][2]= bbb;
        vrt[2][2]= ccc;
        
        vrt[0][3]=-aaa;
        vrt[1][3]= 0.0;
        vrt[2][3]= ccc;
        
        vrt[0][4]=-0.5*aaa;
        vrt[1][4]=-bbb;
        vrt[2][4]= ccc;
        
        vrt[0][5]= 0.5*aaa;
        vrt[1][5]=-bbb;
        vrt[2][5]= ccc;
        
        vrt[0][6]= aaa;
        vrt[1][6]= 0.0;
        vrt[2][6]=-ccc;
        
        vrt[0][7]= 0.5*aaa;
        vrt[1][7]= bbb;
        vrt[2][7]=-ccc;
        
        vrt[0][8]=-0.5*aaa;
        vrt[1][8]= bbb;
        vrt[2][8]=-ccc;
        
        vrt[0][9]=-aaa;
        vrt[1][9]= 0.0;
        vrt[2][9]=-ccc;
        
        vrt[0][10]=-0.5*aaa;
        vrt[1][10]=-bbb;
        vrt[2][10]=-ccc;
        
        vrt[0][11]= 0.5*aaa;
        vrt[1][11]=-bbb;
        vrt[2][11]=-ccc;
        
        edge[0][0]=0;
        edge[1][0]=1;
        edge[0][1]=1;
        edge[1][1]=2;
        edge[0][2]=2;
        edge[1][2]=3;
        edge[0][3]=3;
        edge[1][3]=4;
        edge[0][4]=4;
        edge[1][4]=5;
        edge[0][5]=5;
        edge[1][5]=0;
        edge[0][6]=6;
        edge[1][6]=7;
        edge[0][7]=7;
        edge[1][7]=8;
        edge[0][8]=8;
        edge[1][8]=9;
        edge[0][9]=9;
        edge[1][9]=10;
        edge[0][10]=10;
        edge[1][10]=11;
        edge[0][11]=11;
        edge[1][11]=6;
        edge[0][12]=0;
        edge[1][12]=6;
        edge[0][13]=1;
        edge[1][13]=7;
        edge[0][14]=2;
        edge[1][14]=8;
        edge[0][15]=3;
        edge[1][15]=9;
        edge[0][16]=4;
        edge[1][16]=10;
        edge[0][17]=5;
        edge[1][17]=11;
    }
    
    public void boxOCT() {
        /*
*********************************************************************
         
dl_poly/java GUI routine - truncated octahedral periodic boundaries
         
copyright - daresbury laboratory
author    - w.smith 2001
         
*********************************************************************
         */
        num_vertices=24;
        num_edges=36;
        
        vrt=new double[3][num_vertices];
        edge=new int[2][num_edges];
        
        vrt[0][0]= 0.5*cell[0];
        vrt[1][0]= 0.25*cell[0];
        vrt[2][0]= 0.0;
        
        vrt[0][1]= 0.5*cell[0];
        vrt[1][1]= 0.0;
        vrt[2][1]= 0.25*cell[0];
        
        vrt[0][2]= 0.5*cell[0];
        vrt[1][2]=-0.25*cell[0];
        vrt[2][2]= 0.0;
        
        vrt[0][3]= 0.5*cell[0];
        vrt[1][3]= 0.0;
        vrt[2][3]=-0.25*cell[0];
        
        vrt[0][4]=-0.5*cell[0];
        vrt[1][4]= 0.25*cell[0];
        vrt[2][4]= 0.0;
        
        vrt[0][5]=-0.5*cell[0];
        vrt[1][5]= 0.0;
        vrt[2][5]= 0.25*cell[0];
        
        vrt[0][6]=-0.5*cell[0];
        vrt[1][6]=-0.25*cell[0];
        vrt[2][6]= 0.0;
        
        vrt[0][7]=-0.5*cell[0];
        vrt[1][7]= 0.0;
        vrt[2][7]=-0.25*cell[0];
        
        vrt[0][8]= 0.25*cell[0];
        vrt[1][8]= 0.5*cell[0];
        vrt[2][8]= 0.0;
        
        vrt[0][9]= 0.0;
        vrt[1][9]= 0.5*cell[0];
        vrt[2][9]= 0.25*cell[0];
        
        vrt[0][10]=-0.25*cell[0];
        vrt[1][10]= 0.5*cell[0];
        vrt[2][10]= 0.0;
        
        vrt[0][11]= 0.0;
        vrt[1][11]= 0.5*cell[0];
        vrt[2][11]=-0.25*cell[0];
        
        vrt[0][12]= 0.25*cell[0];
        vrt[1][12]=-0.5*cell[0];
        vrt[2][12]= 0.0;
        
        vrt[0][13]= 0.0;
        vrt[1][13]=-0.5*cell[0];
        vrt[2][13]= 0.25*cell[0];
        
        vrt[0][14]=-0.25*cell[0];
        vrt[1][14]=-0.5*cell[0];
        vrt[2][14]= 0.0;
        
        vrt[0][15]= 0.0;
        vrt[1][15]=-0.5*cell[0];
        vrt[2][15]=-0.25*cell[0];
        
        vrt[0][16]= 0.25*cell[0];
        vrt[1][16]= 0.0;
        vrt[2][16]= 0.5*cell[0];
        
        vrt[0][17]= 0.0;
        vrt[1][17]= 0.25*cell[0];
        vrt[2][17]= 0.5*cell[0];
        
        vrt[0][18]=-0.25*cell[0];
        vrt[1][18]= 0.0;
        vrt[2][18]= 0.5*cell[0];
        
        vrt[0][19]= 0.0;
        vrt[1][19]=-0.25*cell[0];
        vrt[2][19]= 0.5*cell[0];
        
        vrt[0][20]= 0.25*cell[0];
        vrt[1][20]= 0.0;
        vrt[2][20]=-0.5*cell[0];
        
        vrt[0][21]= 0.0;
        vrt[1][21]= 0.25*cell[0];
        vrt[2][21]=-0.5*cell[0];
        
        vrt[0][22]=-0.25*cell[0];
        vrt[1][22]= 0.0;
        vrt[2][22]=-0.5*cell[0];
        
        vrt[0][23]= 0.0;
        vrt[1][23]=-0.25*cell[0];
        vrt[2][23]=-0.5*cell[0];
        
        edge[0][0]=0;
        edge[1][0]=1;
        edge[0][1]=1;
        edge[1][1]=2;
        edge[0][2]=2;
        edge[1][2]=3;
        edge[0][3]=3;
        edge[1][3]=0;
        edge[0][4]=4;
        edge[1][4]=5;
        edge[0][5]=5;
        edge[1][5]=6;
        edge[0][6]=6;
        edge[1][6]=7;
        edge[0][7]=7;
        edge[1][7]=4;
        edge[0][8]=8;
        edge[1][8]=9;
        edge[0][9]=9;
        edge[1][9]=10;
        edge[0][10]=10;
        edge[1][10]=11;
        edge[0][11]=11;
        edge[1][11]=8;
        edge[0][12]=12;
        edge[1][12]=13;
        edge[0][13]=13;
        edge[1][13]=14;
        edge[0][14]=14;
        edge[1][14]=15;
        edge[0][15]=15;
        edge[1][15]=12;
        edge[0][16]=16;
        edge[1][16]=17;
        edge[0][17]=17;
        edge[1][17]=18;
        edge[0][18]=18;
        edge[1][18]=19;
        edge[0][19]=19;
        edge[1][19]=16;
        edge[0][20]=20;
        edge[1][20]=21;
        edge[0][21]=21;
        edge[1][21]=22;
        edge[0][22]=22;
        edge[1][22]=23;
        edge[0][23]=23;
        edge[1][23]=20;
        edge[0][24]=0;
        edge[1][24]=8;
        edge[0][25]=10;
        edge[1][25]=4;
        edge[0][26]=6;
        edge[1][26]=14;
        edge[0][27]=12;
        edge[1][27]=2;
        edge[0][28]=1;
        edge[1][28]=16;
        edge[0][29]=18;
        edge[1][29]=5;
        edge[0][30]=7;
        edge[1][30]=22;
        edge[0][31]=20;
        edge[1][31]=3;
        edge[0][32]=9;
        edge[1][32]=17;
        edge[0][33]=19;
        edge[1][33]=13;
        edge[0][34]=15;
        edge[1][34]=23;
        edge[0][35]=21;
        edge[1][35]=11;
        
    }
    
    public void boxDEC() {
        /*
*********************************************************************
         
dl_poly/java GUI routine - rhombic docecahedral periodic boundaries
         
copyright - daresbury laboratory
author    - w.smith 2001
         
*********************************************************************
         */
        num_vertices=14;
        num_edges=24;
        double ddd=cell[0]/Math.sqrt(2.0);
        
        vrt=new double[3][num_vertices];
        edge=new int[2][num_edges];
        
        vrt[0][0]= 0.0;
        vrt[1][0]= 0.0;
        vrt[2][0]= ddd;
        
        vrt[0][1]= 0.5*cell[0];
        vrt[1][1]= 0.0;
        vrt[2][1]= 0.5*ddd;
        
        vrt[0][2]= 0.0;
        vrt[1][2]= 0.5*cell[0];
        vrt[2][2]= 0.5*ddd;
        
        vrt[0][3]=-0.5*cell[0];
        vrt[1][3]= 0.0;
        vrt[2][3]= 0.5*ddd;
        
        vrt[0][4]= 0.0;
        vrt[1][4]=-0.5*cell[0];
        vrt[2][4]= 0.5*ddd;
        
        vrt[0][5]= 0.5*cell[0];
        vrt[1][5]=-0.5*cell[0];
        vrt[2][5]= 0.0;
        
        vrt[0][6]= 0.5*cell[0];
        vrt[1][6]= 0.5*cell[0];
        vrt[2][6]= 0.0;
        
        vrt[0][7]=-0.5*cell[0];
        vrt[1][7]= 0.5*cell[0];
        vrt[2][7]= 0.0;
        
        vrt[0][8]=-0.5*cell[0];
        vrt[1][8]=-0.5*cell[0];
        vrt[2][8]= 0.0;
        
        vrt[0][9]= 0.5*cell[0];
        vrt[1][9]= 0.0;
        vrt[2][9]=-0.5*ddd;
        
        vrt[0][10]= 0.0;
        vrt[1][10]= 0.5*cell[0];
        vrt[2][10]=-0.5*ddd;
        
        vrt[0][11]=-0.5*cell[0];
        vrt[1][11]= 0.0;
        vrt[2][11]=-0.5*ddd;
        
        vrt[0][12]= 0.0;
        vrt[1][12]=-0.5*cell[0];
        vrt[2][12]=-0.5*ddd;
        
        vrt[0][13]= 0.0;
        vrt[1][13]= 0.0;
        vrt[2][13]=-ddd;
        
        edge[0][0]=0;
        edge[1][0]=1;
        edge[0][1]=0;
        edge[1][1]=2;
        edge[0][2]=0;
        edge[1][2]=3;
        edge[0][3]=0;
        edge[1][3]=4;
        edge[0][4]=1;
        edge[1][4]=5;
        edge[0][5]=1;
        edge[1][5]=6;
        edge[0][6]=2;
        edge[1][6]=6;
        edge[0][7]=2;
        edge[1][7]=7;
        edge[0][8]=3;
        edge[1][8]=7;
        edge[0][9]=3;
        edge[1][9]=8;
        edge[0][10]=4;
        edge[1][10]=5;
        edge[0][11]=4;
        edge[1][11]=8;
        edge[0][12]=9;
        edge[1][12]=5;
        edge[0][13]=9;
        edge[1][13]=6;
        edge[0][14]=10;
        edge[1][14]=6;
        edge[0][15]=10;
        edge[1][15]=7;
        edge[0][16]=11;
        edge[1][16]=7;
        edge[0][17]=11;
        edge[1][17]=8;
        edge[0][18]=12;
        edge[1][18]=5;
        edge[0][19]=12;
        edge[1][19]=8;
        edge[0][20]=13;
        edge[1][20]=9;
        edge[0][21]=13;
        edge[1][21]=10;
        edge[0][22]=13;
        edge[1][22]=11;
        edge[0][23]=13;
        edge[1][23]=12;
    }
    
    void images(double xyz[]) {
        /*
*********************************************************************
         
dl_poly/java GUI routine
         
copyright - daresbury laboratory
author    - w.smith 2001
         
*********************************************************************
         */
        
        double ssx,ssy,ssz;
        double RT2=Math.sqrt(2.0);
        
        if(imcon > 0) {
            
            ssx=(rcell[0]*xyz[0]+rcell[3]*xyz[1]+rcell[6]*xyz[2]);
            ssy=(rcell[1]*xyz[0]+rcell[4]*xyz[1]+rcell[7]*xyz[2]);
            ssz=(rcell[2]*xyz[0]+rcell[5]*xyz[1]+rcell[8]*xyz[2]);
            
            ssx-=BML.nint(ssx);
            ssy-=BML.nint(ssy);
            if(imcon != 6) {
                ssz-=BML.nint(ssz);
                
                if((imcon==4 && (Math.abs(ssx)+Math.abs(ssy)+Math.abs(ssz) >= 0.75)) ||
                (imcon==5 && (Math.abs(ssx)+Math.abs(ssy)+Math.abs(2.0*ssz) >= 1.0))) {
                    ssx-=0.5*BML.sign(ssx);
                    ssy-=0.5*BML.sign(ssy);
                    ssz-=0.5*BML.sign(ssz);
                }
            }
            
            xyz[0]=(cell[0]*ssx+cell[3]*ssy+cell[6]*ssz);
            xyz[1]=(cell[1]*ssx+cell[4]*ssy+cell[7]*ssz);
            xyz[2]=(cell[2]*ssx+cell[5]*ssy+cell[8]*ssz);
        }
    }
    
    public void images(int natms,double xyz[][]) {
        /*
*********************************************************************
         
dl_poly/java GUI routine
         
copyright - daresbury laboratory
author    - w.smith 2001
         
*********************************************************************
         */
        
        double ssx,ssy,ssz;
        double RT2=Math.sqrt(2.0);
        
        if(imcon > 0) {
            
            for(int i=0;i<natms;i++) {
                ssx=(rcell[0]*xyz[0][i]+rcell[3]*xyz[1][i]+rcell[6]*xyz[2][i]);
                ssy=(rcell[1]*xyz[0][i]+rcell[4]*xyz[1][i]+rcell[7]*xyz[2][i]);
                ssz=(rcell[2]*xyz[0][i]+rcell[5]*xyz[1][i]+rcell[8]*xyz[2][i]);
                
                ssx-=BML.nint(ssx);
                ssy-=BML.nint(ssy);
                if(imcon != 6) {
                    ssz-=BML.nint(ssz);
                    
                    if((imcon==4 && (Math.abs(ssx)+Math.abs(ssy)+Math.abs(ssz) >= 0.75)) ||
                    (imcon==5 && (Math.abs(ssx)+Math.abs(ssy)+Math.abs(2.0*ssz) >= 1.0))) {
		        ssx-=0.5*BML.sign(ssx);
                        ssy-=0.5*BML.sign(ssy);
                        ssz-=0.5*BML.sign(ssz);
                    }
                }
                
                xyz[0][i]=(cell[0]*ssx+cell[3]*ssy+cell[6]*ssz);
                xyz[1][i]=(cell[1]*ssx+cell[4]*ssy+cell[7]*ssz);
                xyz[2][i]=(cell[2]*ssx+cell[5]*ssy+cell[8]*ssz);
            }
        }
    }
    
}
