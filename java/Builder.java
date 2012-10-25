import java.io.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.print.*;
import javax.swing.*;

public class Builder extends JComponent implements Printable {
        /*
*********************************************************************

dl_poly/java GUI class for rendering molecular pictures

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
    Editor home;
    Image offscreenImg;
    Graphics offscreenG;
    Config fragment=null;
    boolean link,lpen,water=true,lmove=false;
    int oper=0,nangs,tilex,tiley,boxtyp=1,nhdel=0,pen=0,shx=0,shy=0;
    int ngroup,sx,sy,fx,fy,keyopt,mark0=-1,mark1=-1,mark2=-1;
    int[] lst,xpos,ypos,zpos,dia,xvrt,yvrt,zvrt,group;
    double scale,xx2,yy2,zz2,dd2,fac;
    boolean[] led;
    String news,activity,frgname;
    Element atom;
    private int MXGROUP=100;

    public Builder(Editor here) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */

        // Define the parent Editor

        home=here;

        // set up initial CONFIG if currently null

        if(home.config == null || home.config.natms == 0)
            newBuild();


        // Register Mouse Listeners

        addMouseListener(new MousePoints());
        addMouseMotionListener(new MouseMoves());

        // Register object for screen resize

        addComponentListener(new MyScreenSizer());

    }

    public void newBuild(){
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/

        link=false;
        activity="none";
        home.config=new Config();
        home.config.pbc.imcon=0;
        home.config.pbc.cell[0]=10;
        home.config.pbc.cell[4]=10;
        home.config.pbc.cell[8]=10;
        home.config.pbc.buildBoundary(1);
        home.config.structure=new Structure();
        led=new boolean[home.config.pbc.num_edges];
        xvrt=new int[home.config.pbc.num_vertices];
        yvrt=new int[home.config.pbc.num_vertices];
        zvrt=new int[home.config.pbc.num_vertices];
        news="NULL";
        oper=0;
    }

    public void printOut() {
        /*
*********************************************************************

dl_poly/java GUI routine to print the GUI

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        PrinterJob prnjob=PrinterJob.getPrinterJob();
        PageFormat pf=prnjob.defaultPage();
        pf.setOrientation(PageFormat.LANDSCAPE);
        prnjob.setPrintable(this,pf);
        try {
            if(prnjob.printDialog()) {
                home.println("Initiating print .......");
                prnjob.print();
                home.println("Print complete");
            }
            else {
                home.println("Print cancelled");
            }
        }
        catch(PrinterException pe) {
            home.println("Error - problem with print "+pe.toString());
        }

    }

    public int print(Graphics g,PageFormat pf,int pageIndex) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        if(pageIndex > 0) {
            return Printable.NO_SUCH_PAGE;
        }
        Graphics2D g2d=(Graphics2D)g;
        g2d.translate(pf.getImageableX(),pf.getImageableY());
        paint(g2d);
        return Printable.PAGE_EXISTS;
    }

    public void paint(Graphics g) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        boolean showbond;
        int j,k,u,w,rr,xd,yd,zd,ka,kb,ja,jb;

        tilex=getSize().width;
        tiley=getSize().height;
        offscreenImg = createImage(tilex,tiley);
        offscreenG = offscreenImg.getGraphics();
        offscreenG.setColor(home.art.scrn);
        offscreenG.fillRect(0,0,tilex,tiley);
        offscreenG.setColor(Color.black);
        g.setFont(new Font("Monaco",Font.PLAIN,10));
        if(home.edit) {
            offscreenG.drawString("EDIT: "+news,10,10);
            offscreenG.drawString(atom.zsym,tilex-50,10);
            offscreenG.drawLine(tilex/2+shx,tiley/2+10+shy,tilex/2+shx,tiley/2-10+shy);
            offscreenG.drawLine(tilex/2+10+shx,tiley/2+shy,tilex/2-10+shx,tiley/2+shy);
            if(boxtyp > 0 && home.config != null) {
                g.setFont(new Font("Monaco",Font.PLAIN,8));
                if(home.config.pbc.imcon > 0)
                    if(oper == 12)
                        offscreenG.drawString("Box size:"+BML.fmt(home.config.pbc.cell[0],10)+
                BML.fmt(home.config.pbc.cell[4],10)+
                BML.fmt(home.config.pbc.cell[8],10),10,tiley-10);
                g.setFont(new Font("Monaco",Font.PLAIN,10));
            }
            if(oper == 13 && fragment != null)
                offscreenG.drawString("FRAG: "+frgname,10,tiley-10);
        }
        else if(oper == 6 || oper == 7)
            offscreenG.drawString("VIEW: "+news,10,10);
        else
            offscreenG.drawString("VIEW",10,10);

        if(home.config != null) {

            //      order in z direction

            zOrder();

            //      Mark clicked atom

            if(mark0>=0) {
                w=dia[mark0];
                offscreenG.setColor(Color.red);
                offscreenG.fillOval(xpos[mark0]-4,ypos[mark0]-4,w+8,w+8);
            }
            offscreenG.setColor(Color.black);

            if(home.config.pbc.imcon > 0) {

                //      check visibility of cell edges

                checkEdges();

                //      Draw hidden box edges

                for (int i=0;i<home.config.pbc.num_edges;i++) {
                    if(!led[i]) {
                        j=home.config.pbc.edge[0][i];
                        k=home.config.pbc.edge[1][i];
                        offscreenG.drawLine(xvrt[j],yvrt[j],xvrt[k],yvrt[k]);
                    }
                }

            }

            //      Draw atoms

            for (int i=0;i<home.config.natms;i++) {
                if((water ||!(home.config.atoms[lst[i]].zsym.indexOf("OW") == 0 ||
                home.config.atoms[lst[i]].zsym.indexOf("HW") == 0)) &&
                !(home.config.atoms[lst[i]].zsym.indexOf("QW") == 0)) {
                    j=lst[i];
                    w=dia[j];
                    if(j == mark0 || j == mark1 || j == mark2 || (ngroup > 0 && group[j] > -1)) {
                        offscreenG.setColor(Color.red);
                        offscreenG.fillOval(xpos[j]-4,ypos[j]-4,w+8,w+8);
                    }
                    offscreenG.fillOval(xpos[j]-2,ypos[j]-2,w+4,w+4);
                    offscreenG.setColor(home.config.atoms[j].zcol);
                    offscreenG.fillOval(xpos[j],ypos[j],w,w);
                    offscreenG.setColor(Color.black);
                    w=w/2;
                    if(home.config.structure.nbnds>0) {
                        for(int m=0;m<home.config.structure.lbnd[j];m++) {
                            showbond=true;
                            k=home.config.structure.bond[m][j];
                            if(home.config.pbc.imcon > 0) {
                                if(Math.abs((double)(xpos[k]-xpos[j])) > scale*4)showbond=false;
                                if(Math.abs((double)(ypos[k]-ypos[j])) > scale*4)showbond=false;
                                if(Math.abs((double)(zpos[k]-zpos[j])) > scale*4)showbond=false;
                                if(home.edit)showbond=true;
                            }
                            if(!home.showbonds)showbond=false;
                            if(home.config.atoms[k].zsym.indexOf("QW") == 0)showbond=false;
                            if(showbond) {
                                u=dia[k]/2;
                                if(j>k && zpos[k]==zpos[j]) {
                                    xd=xpos[k]-xpos[j]+u-w;
                                    yd=ypos[k]-ypos[j]+u-w;
                                    rr=(int)Math.sqrt(xd*xd+yd*yd);
                                    ja=xpos[j]+w+(int)((double)(xd*w)/rr);
                                    jb=ypos[j]+w+(int)((double)(yd*w)/rr);
                                    ka=xpos[k]+u-(int)((double)(xd*u)/rr);
                                    kb=ypos[k]+u-(int)((double)(yd*u)/rr);
                                    offscreenG.drawLine(ja,jb+1,ka,kb+1);
                                    offscreenG.drawLine(ja+1,jb,ka+1,kb);
                                    offscreenG.drawLine(ja,jb,ka,kb);
                                    offscreenG.drawLine(ja-1,jb,ka-1,kb);
                                    offscreenG.drawLine(ja,jb-1,ka,kb-1);
                                }
                                else if(zpos[k]>zpos[j]) {
                                    xd=xpos[k]-xpos[j]+u-w;
                                    yd=ypos[k]-ypos[j]+u-w;
                                    zd=zpos[k]-zpos[j]+u-w;
                                    rr=(int)Math.sqrt(xd*xd+yd*yd+zd*zd);
                                    ja=xpos[j]+w+(int)((double)(xd*w)/rr);
                                    jb=ypos[j]+w+(int)((double)(yd*w)/rr);
                                    ka=xpos[k]+u;
                                    kb=ypos[k]+u;
                                    offscreenG.drawLine(ja,jb+1,ka,kb+1);
                                    offscreenG.drawLine(ja+1,jb,ka+1,kb);
                                    offscreenG.drawLine(ja,jb,ka,kb);
                                    offscreenG.drawLine(ja-1,jb,ka-1,kb);
                                    offscreenG.drawLine(ja,jb-1,ka,kb-1);
                                }
                            }
                        }
                    }
                }
            }

            if(home.config.pbc.imcon > 0) {

                //      Draw non hidden box edges

                for (int i=0;i<home.config.pbc.num_edges;i++) {
                    if(led[i]) {
                        j=home.config.pbc.edge[0][i];
                        k=home.config.pbc.edge[1][i];
                        offscreenG.drawLine(xvrt[j],yvrt[j],xvrt[k],yvrt[k]);
                    }
                }
            }
            if(lmove && oper == 5) {
                offscreenG.drawLine(sx,sy,fx,sy);
                offscreenG.drawLine(sx,sy,sx,fy);
                offscreenG.drawLine(sx,fy,fx,fy);
                offscreenG.drawLine(fx,sy,fx,fy);
            }
        }
        g.drawImage(offscreenImg,0,0,this);
    }

    void zOrder() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2003

*********************************************************************
        */
        double rad;

        if(home.config.natms > 0){
            if(home.edit)
                fac=0.3;
            else
                fac=0.5;
            lst=new int[home.config.natms];
            if(home.edit && home.config.pbc.imcon > 0 && home.config.pbc==null)
                home.config.pbc.buildBoundary(home.config.pbc.imcon);
            xvrt=new int[home.config.pbc.num_vertices];
            yvrt=new int[home.config.pbc.num_vertices];
            zvrt=new int[home.config.pbc.num_vertices];
            led=new boolean[home.config.pbc.num_edges];
            xpos=new int[home.config.natms];
            ypos=new int[home.config.natms];
            zpos=new int[home.config.natms];
            dia=new int[home.config.natms];
            for (int i=0;i<home.config.natms;i++) {
                lst[i]=i;
                rad=fac*home.config.atoms[i].zrad;
                if(home.config.atoms[i].dotify)
                    rad/=3.0;;
                dia[i]=(int)(2*scale*rad);
                xpos[i]=(int)(scale*(home.config.xyz[0][i]-rad))+tilex/2+shx;
                ypos[i]=-(int)(scale*(home.config.xyz[2][i]+rad))+tiley/2+shy;
                zpos[i]=-(int)(scale*home.config.xyz[1][i]);
            }
            for (int i=0;i<home.config.pbc.num_vertices;i++) {
                xvrt[i]=(int)(scale*home.config.pbc.vrt[0][i])+tilex/2+shx;
                yvrt[i]=-(int)(scale*home.config.pbc.vrt[2][i])+tiley/2+shy;
                zvrt[i]=-(int)(scale*home.config.pbc.vrt[1][i]);
            }
        }
        AML.ShellSort0(home.config.natms,lst,zpos);
    }

    void restore() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        oper=0;
        fac=0.5;
        ngroup=0;
        mark0=-1;
        mark1=-1;
        mark2=-1;
        home.safe=true;
        if(home.config != null && home.config.pbc.imcon > 0)
            scale=setBoxScale();
        repaint();
    }

    void checkEdges() {
        /*
*********************************************************************

dl_poly/java GUI routine to check visibility of box edges

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        int w,n1,n2,xa,xb,ya,yb;
        double px,py,pp,sx,sy,ss,uu,dd,ww,rr,zz,tt;

        if(home.config != null) {
            for (int n=0;n<home.config.pbc.num_edges;n++) {
                led[n]=true;
                n1=home.config.pbc.edge[0][n];
                n2=home.config.pbc.edge[1][n];
                xa=Math.min(xvrt[n1],xvrt[n2]);
                xb=Math.max(xvrt[n1],xvrt[n2]);
                ya=Math.min(yvrt[n1],yvrt[n2]);
                yb=Math.max(yvrt[n1],yvrt[n2]);
                OUT:
                    for(int i=0;i<home.config.natms;i++) {
                        if(xpos[i]>xa-dia[i] && xpos[i]<xb && ypos[i]>ya-dia[i] && ypos[i]<yb) {
                            w=dia[i]/2;
                            ww=0.25*(double)(dia[i]*dia[i]);
                            px=(double)(xvrt[n2]-xvrt[n1]);
                            py=(double)(yvrt[n2]-yvrt[n1]);
                            pp=Math.sqrt(px*px+py*py);
                            sx=(double)(xpos[i]+w-xvrt[n1]);
                            sy=(double)(ypos[i]+w-yvrt[n1]);
                            uu=(sx*px+sy*py)/pp;
                            ss=sx*sx+sy*sy-uu*uu;
                            if(ss < ww) {
                                dd=(double)(zpos[i]-zvrt[n1]);
                                zz=(double)(zvrt[n2]-zvrt[n1]);
                                tt=uu*zz/pp;
                                if (tt < dd) {
                                    led[n]=false;
                                    break OUT;
                                }
                            }
                        }
                    }
            }
        }
    }

    double setBoxScale() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2007

*********************************************************************
         */
        double size=BML.max(home.config.pbc.cell[0]+home.config.pbc.cell[3]+home.config.pbc.cell[6],
        home.config.pbc.cell[1]+home.config.pbc.cell[4]+home.config.pbc.cell[7],
        home.config.pbc.cell[2]+home.config.pbc.cell[5]+home.config.pbc.cell[8]);
        return (0.75*Math.min(tilex,tiley)/Math.max(size,1.0));
    }

    void zoom(int k) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        double rad;

        //      if(home.edit)return;

        if (home.config != null) {
            if(scale==0.0) {
                if(home.config.pbc.imcon > 0)
                    scale=setBoxScale();
                else
                    scale=37.5;
            }
            else {
                if(k > 0)
                    scale/=1.5;
                else
                    scale*=1.5;
            }
            repaint();
        }
    }

    public void displace(int k,double inc) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        double rad;

        if(home.edit)return;

        if(home.config != null) {
            for (int i=0;i<home.config.natms;i++)
                home.config.xyz[k][i]+=inc;
            for (int i=0;i<home.config.pbc.num_vertices;i++)
                home.config.pbc.vrt[k][i]+=inc;
            if (k == 0) {
                for (int i=0;i<home.config.natms;i++) {
                    rad=fac*home.config.atoms[i].zrad;
                    xpos[i]=(int)(scale*(home.config.xyz[0][i]-rad))+tilex/2+shx;
                }
                for (int i=0;i<home.config.pbc.num_vertices;i++) {
                    xvrt[i]=(int)(scale*home.config.pbc.vrt[0][i])+tilex/2+shx;
                }
            }
            else if (k == 1) {
                for (int i=0;i<home.config.natms;i++) {
                    rad=fac*home.config.atoms[i].zrad;
                    ypos[i]=-(int)(scale*(home.config.xyz[2][i]+rad))+tiley/2+shy;
                }
                for (int i=0;i<home.config.pbc.num_vertices;i++) {
                    yvrt[i]=-(int)(scale*home.config.pbc.vrt[2][i])+tiley/2+shy;
                }
            }
            else if (k == 2) {
                for (int i=0;i<home.config.natms;i++) {
                    rad=fac*home.config.atoms[i].zrad;
                    zpos[i]=-(int)(scale*home.config.xyz[1][i]);
                }
                for (int i=0;i<home.config.pbc.num_vertices;i++) {
                    zvrt[i]=-(int)(scale*home.config.pbc.vrt[1][i]);
                }
            }
            repaint();
        }
    }

    public void rotate(int j,int k,double c,double s) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        double u,v,rad;

        if(home.config != null) {
            for (int i=0;i<home.config.natms;i++) {
                u=home.config.xyz[j][i];
                v=home.config.xyz[k][i];
                home.config.xyz[j][i]=(c*u-s*v);
                home.config.xyz[k][i]=(s*u+c*v);
            }
            for (int i=0;i<home.config.pbc.num_vertices;i++) {
                u=home.config.pbc.vrt[j][i];
                v=home.config.pbc.vrt[k][i];
                home.config.pbc.vrt[j][i]=(c*u-s*v);
                home.config.pbc.vrt[k][i]=(s*u+c*v);
            }
            if ((j == 0)&&(k == 1)) {
                for (int i=0;i<home.config.natms;i++) {
                    rad=fac*home.config.atoms[i].zrad;
                    xpos[i]=(int)(scale*(home.config.xyz[0][i]-rad))+tilex/2+shx;
                    ypos[i]=-(int)(scale*(home.config.xyz[2][i]+rad))+tiley/2+shy;
                }
                for (int i=0;i<home.config.pbc.num_vertices;i++) {
                    xvrt[i]=(int)(scale*home.config.pbc.vrt[0][i])+tilex/2+shx;
                    yvrt[i]=-(int)(scale*home.config.pbc.vrt[2][i])+tiley/2+shy;
                }
            }
            else if ((j == 2)&&(k == 0)) {
                for (int i=0;i<home.config.natms;i++) {
                    rad=fac*home.config.atoms[i].zrad;
                    xpos[i]=(int)(scale*(home.config.xyz[0][i]-rad))+tilex/2+shx;
                    zpos[i]=-(int)(scale*home.config.xyz[1][i]);
                }
                for (int i=0;i<home.config.pbc.num_vertices;i++) {
                    xvrt[i]=(int)(scale*home.config.pbc.vrt[0][i])+tilex/2+shx;
                    zvrt[i]=-(int)(scale*home.config.pbc.vrt[1][i]);
                }
            }
            else if ((j == 1)&&(k == 2)) {
                for (int i=0;i<home.config.natms;i++) {
                    rad=fac*home.config.atoms[i].zrad;
                    ypos[i]=-(int)(scale*(home.config.xyz[2][i]+rad))+tiley/2+shx;
                    zpos[i]=-(int)(scale*home.config.xyz[1][i]);
                }
                for (int i=0;i<home.config.pbc.num_vertices;i++) {
                    yvrt[i]=-(int)(scale*home.config.pbc.vrt[2][i])+tiley/2+shy;
                    zvrt[i]=-(int)(scale*home.config.pbc.vrt[1][i]);
                }
            }
            repaint();
        }
    }

    void shiftAtoms() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        double rad;

        if (home.config != null) {

            for (int i=0;i<home.config.natms;i++) {
                home.config.xyz[0][i]+=(fx-sx)/scale;
                home.config.xyz[2][i]-=(fy-sy)/scale;
            }
            for (int i=0;i<home.config.pbc.num_vertices;i++) {
                home.config.pbc.vrt[0][i]+=(fx-sx)/scale;
                home.config.pbc.vrt[2][i]-=(fy-sy)/scale;
            }
            repaint();
        }
    }

    void turnAtoms() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        double rad,ca,cb,sa,sb,tx,ty,tz;

        if(home.edit)return;

        if (home.config != null) {
            ca=Math.cos((sx-fx)*Math.PI/tilex);
            sa=Math.sin((sx-fx)*Math.PI/tilex);
            cb=Math.cos((sy-fy)*Math.PI/tiley);
            sb=Math.sin((sy-fy)*Math.PI/tiley);
            for (int i=0;i<home.config.natms;i++) {
                tx=home.config.xyz[0][i];
                ty=home.config.xyz[2][i];
                tz=home.config.xyz[1][i];
                home.config.xyz[0][i]=ca*tx+sa*tz;
                home.config.xyz[1][i]=-sa*cb*tx+sb*ty+ca*cb*tz;
                home.config.xyz[2][i]=sa*sb*tx+cb*ty-ca*sb*tz;
            }
            for (int i=0;i<home.config.pbc.num_vertices;i++) {
                tx=home.config.pbc.vrt[0][i];
                ty=home.config.pbc.vrt[2][i];
                tz=home.config.pbc.vrt[1][i];
                home.config.pbc.vrt[0][i]=ca*tx+sa*tz;
                home.config.pbc.vrt[1][i]=-sa*cb*tx+sb*ty+ca*cb*tz;
                home.config.pbc.vrt[2][i]=sa*sb*tx+cb*ty-ca*sb*tz;
            }
            repaint();
        }
    }

    void identifyAtom() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
         */
        int k;
        double aaa,xx1,yy1,zz1,dd1;

        if(oper == 0) {

            k=getAtom();

            if(k < 0) {
                mark0=-1;
                mark1=-1;
                mark2=-1;
                repaint();
                return;
            }

            // Identify clicked atoms

            if(mark0==k) {
                if(home.edit) {

                    // Swap atom identities if different

                    if(!atom.zsym.equals(home.config.atoms[k].zsym)) {
                        home.config.atoms[k]=new Element(atom.zsym);
                        home.println("Atom substitution:");
                        home.println("ATOM: "+BML.fmt(k+1,10)+"     "+BML.fmt(home.config.atoms[k].zsym,8));
                        home.safe=false;
                    }
                }
                else {
                    mark0=-1;
                    mark1=-1;
                    mark2=-1;
                }
                repaint();
            }
            else {
                mark2=mark1;
                mark1=mark0;
                mark0=k;
                if(mark1 < 0)
                    home.println("ATOM: "+BML.fmt(k+1,10)+"     "+BML.fmt(home.config.atoms[k].zsym,8));
                else {
                    xx1=home.config.xyz[0][mark1]-home.config.xyz[0][mark0];
                    yy1=home.config.xyz[1][mark1]-home.config.xyz[1][mark0];
                    zz1=home.config.xyz[2][mark1]-home.config.xyz[2][mark0];
                    dd1=Math.sqrt(Math.pow(xx1,2)+Math.pow(yy1,2)+Math.pow(zz1,2));
                    if(mark2 < 0) {
                        home.println("ATOM: "+BML.fmt(k+1,10)+"     "+
                        BML.fmt(home.config.atoms[k].zsym,8)+
                        "     "+BML.fmt(dd1,10));
                    }
                    else {
                        aaa=-(xx1*xx2+yy1*yy2+zz1*zz2)/(dd1*dd2);
                        if(aaa > 1)aaa=1;
                        if(aaa < -1)aaa=-1;
                        aaa=(180/Math.PI)*Math.acos(aaa);
                        home.println("ATOM: "+BML.fmt(k+1,10)+"     "+
                        BML.fmt(home.config.atoms[k].zsym,8)+
                        "     "+BML.fmt(dd1,10)+"     "+
                        BML.fmt(aaa,10));
                    }
                    xx2=xx1;
                    yy2=yy1;
                    zz2=zz1;
                    dd2=dd1;
                }
                repaint();
            }
        }
    }

    int getAtom() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
         */
        int k;
        double ddd,width,rscale;

        // Get identity of clicked atom

        k=-1;
        if(home.config != null) {
            rscale=1.0/scale;
            for(int i=0;i<home.config.natms;i++) {
                ddd=Math.sqrt(Math.pow(home.config.xyz[0][i]-(sx-shx)*rscale,2)+
                Math.pow(home.config.xyz[2][i]-(sy+shy)*rscale,2));
                width=fac*home.config.atoms[i].zrad;
                if(home.config.atoms[i].dotify)
                    width=fac*home.config.atoms[i].zrad/3.0;
                if(ddd < width) {
                    if(k < 0)
                        k=i;
                    else if(home.config.xyz[1][k] > home.config.xyz[1][i])
                        k=i;
                }
            }
        }
        return k;
    }

    int cgropt(int natms,double step,int ida[],double hnrm[],
    double fff[],double grad[],double hhh[][],double dxyz[][]) {
        /*
*********************************************************************

dl-poly/java GUI conjugate gradient routine

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
         */
        int j;
        double ggg,stride,gam2;

        // Magnitude of current gradient vector

        ggg=0.0;
        for(int i=0;i<natms;i++) {
            ggg+=(dxyz[0][i]*dxyz[0][i]+dxyz[1][i]*dxyz[1][i]+
            dxyz[2][i]*dxyz[2][i]);
        }
        ggg=Math.sqrt(ggg);

        if(keyopt == 0) {

            // Set original search direction (vector hhh)

            keyopt=1;
            hnrm[0]=ggg;
            grad[0]=ggg;
            grad[2]=ggg;
            fff[2]=fff[0];

            for(int i=0;i<natms;i++) {
                j=ida[i];
                hhh[0][i]=dxyz[0][i];
                hhh[1][i]=dxyz[1][i];
                hhh[2][i]=dxyz[2][i];
                home.config.xyz[0][j]+=(step*hhh[0][i]);
                home.config.xyz[1][j]+=(step*hhh[1][i]);
                home.config.xyz[2][j]+=(step*hhh[2][i]);
            }
        }
        else if(keyopt == 1) {
            // Line search along chosen direction

            stride=step;
            fff[1]=fff[2];
            fff[2]=fff[0];
            grad[1]=grad[2];

            grad[2]=0.0;
            for(int i=0;i<natms;i++) {
                grad[2]+=(hhh[0][i]*dxyz[0][i]+hhh[1][i]*dxyz[1][i]+
                hhh[2][i]*dxyz[2][i]);
            }
            grad[2]=grad[2]/hnrm[0];

            // Linear extrapolation to minimum

            if(grad[2] < 0) {
                stride=step*grad[2]/(grad[1]-grad[2]);
                keyopt=2;
            }
            for(int i=0;i<natms;i++) {
                j=ida[i];
                home.config.xyz[0][j]+=(stride*hhh[0][i]);
                home.config.xyz[1][j]+=(stride*hhh[1][i]);
                home.config.xyz[2][j]+=(stride*hhh[2][i]);
            }
        }
        else if(keyopt == 2) {
            fff[1]=fff[2];
            fff[2]=fff[0];

            // Check for global convergence

            if(Math.abs(ggg/natms) < 0.0000001) {
                return 999;
            }

            // Construct conjugate search vector

            gam2=Math.pow((ggg/grad[0]),2);

            grad[0]=ggg;
            hnrm[0]=0.0;
            grad[2]=0.0;
            for(int i=0;i<natms;i++) {
                hhh[0][i]=dxyz[0][i]+gam2*hhh[0][i];
                hhh[1][i]=dxyz[1][i]+gam2*hhh[1][i];
                hhh[2][i]=dxyz[2][i]+gam2*hhh[2][i];
                hnrm[0]+=(hhh[0][i]*hhh[0][i]+hhh[1][i]*hhh[1][i]+
                hhh[2][i]*hhh[2][i]);
                grad[2]+=(hhh[0][i]*dxyz[0][i]+hhh[1][i]*dxyz[1][i]+
                hhh[2][i]*dxyz[2][i]);
            }
            hnrm[0]=Math.sqrt(hnrm[0]);
            grad[2]/=hnrm[0];
            for(int i=0;i<natms;i++) {
                j=ida[i];
                home.config.xyz[0][j]+=(step*hhh[0][i]);
                home.config.xyz[1][j]+=(step*hhh[1][i]);
                home.config.xyz[2][j]+=(step*hhh[2][i]);
            }
            keyopt=1;
        }
        return keyopt;
    }

    void moleculeBuilder() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
         */
        int j,k,m;
        double ddd,rscale;

        rscale=1.0/scale;

        if(oper == 1) {

            // Build molecular structure

            OUT:
                if(lpen) {

                    // extend arrays if atom count equals current array sizes

                    if((home.config.atoms.length - home.config.natms) < 2) {
                        home.config.resizeArrays();
                        home.config.structure.resizeBondArrays();
                    }
                    if((home.config.structure.join[0].length - home.config.structure.nbnds) < 2)
                        home.config.structure.resizeJoinArray();

                    // Check identity of new atom

                    k=-1;

                    for(int i=0;i<home.config.natms;i++) {
                        ddd=Math.sqrt(Math.pow(home.config.xyz[0][i]-(sx-shx)*rscale,2)+
                        Math.pow(home.config.xyz[2][i]-(sy+shy)*rscale,2));
                        if(ddd < fac*home.config.atoms[i].zrad) {
                            if(k == -1)
                                k=i+1;
                            else if(home.config.xyz[1][i] < home.config.xyz[1][k])
                                k=i+1;
                        }
                    }

                    // Same atom as last time - finish current fragment

                    if(k == home.config.natms) {
                        link=false;
                        lpen=false;
                        break OUT;
                    }

                    if(k < 0)k=home.config.natms;

                    // Add bond to structure

                    if(link) {
                        j=Math.min(home.config.natms,k)-1;
                        if(k == home.config.natms) {

                            // New atom - bonds to previous entered atom
                            repaint();

                            m=home.config.natms;
                            home.config.structure.bond[home.config.structure.lbnd[m]][m]=j;
                            home.config.structure.bond[home.config.structure.lbnd[j]][j]=m;
                            home.config.structure.lbnd[m]++;
                            home.config.structure.lbnd[j]++;
                            home.config.structure.join[0][home.config.structure.nbnds]=Math.min(j,m);
                            home.config.structure.join[1][home.config.structure.nbnds]=Math.max(j,m);
                            home.config.structure.nbnds++;
                        }
                        else {

                            // Old atom - connect to existing fragment

                            m=home.config.natms-1;
                            home.config.structure.bond[home.config.structure.lbnd[m]][m]=j;
                            home.config.structure.bond[home.config.structure.lbnd[j]][j]=m;
                            home.config.structure.lbnd[m]++;
                            home.config.structure.lbnd[j]++;
                            home.config.structure.join[0][home.config.structure.nbnds]=Math.min(j,m);
                            home.config.structure.join[1][home.config.structure.nbnds]=Math.max(j,m);
                            home.config.structure.nbnds++;
                        }
                    }

                    if(k == home.config.natms) {

                        // Add new atom to fragment

                        home.config.atoms[home.config.natms]=new Element(atom.zsym);
                        home.config.xyz[0][home.config.natms]=(sx-shx)*rscale;
                        home.config.xyz[2][home.config.natms]=(sy+shy)*rscale;
                        home.config.xyz[1][home.config.natms]=0.0;
                        link=true;
                        home.safe=false;
                        home.config.natms++;
                    }
                    else {

                        // Old atom - finish current fragment

                        link=false;
                        lpen=false;
                        home.config.xyz[1][home.config.natms-1]=home.config.xyz[1][k-1];
                        break OUT;
                    }
                }
                else
                    lpen=true;

                repaint();
        }
    }

    void linkMaker() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
         */
        if(oper == 2) {
            if(mark0 < 0) {
                mark0=getAtom();
                repaint();
            }
            else if(mark1 < 0) {
                mark1=getAtom();
                if(mark0 == mark1) {
                    mark1=-1;
                    repaint();
                }
                else if(mark1 >= 0) {
                    if(!duplicateBond()) {
                        home.config.structure.bond[home.config.structure.lbnd[mark0]][mark0]=mark1;
                        home.config.structure.bond[home.config.structure.lbnd[mark1]][mark1]=mark0;
                        home.config.structure.lbnd[mark0]++;
                        home.config.structure.lbnd[mark1]++;
                        home.config.structure.join[0][home.config.structure.nbnds]=Math.min(mark0,mark1);
                        home.config.structure.join[1][home.config.structure.nbnds]=Math.max(mark0,mark1);
                        home.config.structure.nbnds++;
                        mark0=-1;
                        mark1=-1;
                        repaint();
                    }
                }
            }
            if(home.config.structure.nbnds == home.config.structure.join[0].length)
                    home.config.structure.resizeJoinArray();

        }
    }

    void makeGroup() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
         */
        int k;

        k=getAtom();
        if(k < 0) return;

        if(ngroup == 0) {
            group=new int[MXGROUP];
            while(k >= group.length)
                resizeGroupArray();
            for(int i=0;i<MXGROUP;i++)
                group[i]=-1;
            group[k]=1;
            ngroup++;
        }
        else if(group[k] > 0) {
            group=null;
            ngroup=0;
        }
        else {
            while(k >= group.length)
                resizeGroupArray();
            group[k]=1;
            ngroup++;
        }

        repaint();
    }

    void frameMolecule() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
         */
        double rscale;
        int bx,ex,by,ey;

        sx=sx-tilex/2;
        sy=tiley/2-sy;
        fx=fx-tilex/2;
        fy=tiley/2-fy;
        rscale=1.0/scale;
        bx=Math.min(sx,fx)-shx;
        ex=Math.max(sx,fx)-shx;
        by=Math.min(sy,fy)+shy;
        ey=Math.max(sy,fy)+shy;

        if(ngroup == 0) {
            group=new int[MXGROUP];
            for(int i=0;i<MXGROUP;i++)
                group[i]=-1;
        }

        while(group.length < home.config.atoms.length)
            resizeGroupArray();

        for(int i=0;i<home.config.natms;i++) {
            if((home.config.xyz[0][i] >= bx*rscale) && (home.config.xyz[0][i] <= ex*rscale) &&
            (home.config.xyz[2][i] >= by*rscale) && (home.config.xyz[2][i] <= ey*rscale)) {
                group[i]=1;
                ngroup++;
            }
        }
        repaint();
    }

    void addDeleteHydrogen() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
       */
        int k,n;

        if(oper == 10) {
            if(ngroup == 0) {

                // Get identity of clicked atom

                k=getAtom();
                if(k >= 0) {
                    nhdel=0;
                    for (int i=0;i<home.config.structure.lbnd[k];i++)
                        if(home.config.atoms[home.config.structure.bond[i][k]].znum == 1) nhdel++;
                    if(nhdel > 0)
                        delHAtoms(k);
                    else
                        hydrogenAtom(k);
                }
                else {
                    nhdel=0;
                    n=home.config.natms;
                    for(int i=0;i<n;i++)
                        if(home.config.atoms[i].znum == 1) nhdel++;
                    if(nhdel > 0) {
                        delHAtoms(-1);
                    }
                    else {
                        for(int i=0;i<n;i++) {
                            hydrogenAtom(i);
                        }
                    }
                }
            }
            else {
                nhdel=0;
                n=home.config.natms;
                for(int i=0;i<n;i++)
                    if(group[i] > 0 && home.config.atoms[i].znum == 1) nhdel++;
                if(nhdel > 0) {
                    delHAtoms(-1);
                }
                else {
                    for (int i=0;i<n;i++) {
                        if(group[i] > 0) {
                            hydrogenAtom(i);
                        }
                    }
                }
                ngroup=0;
            }
            home.safe=false;
            repaint();
        }
    }

    void hydrogenAtom(int k) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
       */
        String str;
        int nhy,hyb,val;
        double stick;

        nhy=0;
        hyb=0;
        val=0;
        str=home.config.atoms[k].zsym;
        if( str.indexOf("OW") >= 0)
            stick=home.config.atoms[k].zrad+(new Element("HW")).zrad;
        else
            stick=home.config.atoms[k].zrad+(new Element("H_")).zrad;
        if(str.indexOf("C_1") >= 0) {
            hyb=1;
            val=2;
            nhy=2-home.config.structure.lbnd[k];
        }
        else if (str.indexOf("C_2") >= 0) {
            hyb=2;
            val=3;
            nhy=3-home.config.structure.lbnd[k];
        }
        else if (str.indexOf("C_R") >= 0) {
            hyb=2;
            val=3;
            nhy=3-home.config.structure.lbnd[k];
        }
        else if (str.indexOf("C_3") >= 0) {
            hyb=3;
            val=4;
            nhy=4-home.config.structure.lbnd[k];
        }
        else if (str.indexOf("O_3") >= 0) {
            hyb=3;
            val=2;
            nhy=2-home.config.structure.lbnd[k];
        }
        else if (str.indexOf("OW") >= 0) {
            hyb=3;
            val=2;
            nhy=2-home.config.structure.lbnd[k];
        }
        else if (str.indexOf("N_2") >= 0) {
            hyb=2;
            val=2;
            nhy=2-home.config.structure.lbnd[k];
        }
        else if (str.indexOf("N_3") >= 0) {
            hyb=3;
            val=3;
            nhy=3-home.config.structure.lbnd[k];
        }
        else if (str.indexOf("S_3") >= 0) {
            hyb=3;
            val=2;
            nhy=2-home.config.structure.lbnd[k];
        }
        else if (str.indexOf("P_2") >= 0) {
            hyb=2;
            val=3;
            nhy=3-home.config.structure.lbnd[k];
        }
        else if (str.indexOf("P_3") >= 0) {
            hyb=3;
            val=4;
            nhy=4-home.config.structure.lbnd[k];
        }
        if(nhy > 0)

            addHAtom(k,nhy,hyb,val,stick);
    }

    void addHAtom(int k,int nhy,int hyb,int val,double stick) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
         */
        int a,b,c;
        double[] vv,va,vb,vc,vd;
        double xx,yy,zz,ct,st,cp,sp,atn;

        // extend arrays if required

        if((home.config.atoms.length - home.config.natms) <= nhy) {
            home.config.resizeArrays();
            home.config.structure.resizeBondArrays();
        }
        if((home.config.structure.join[0].length - home.config.structure.nbnds) <= nhy)
            home.config.structure.resizeJoinArray();

        if(hyb == 1 && nhy == 1 && val == 2) {
            a=home.config.structure.bond[0][k];
            va=BML.ndxyz(k,a,home.config.xyz);
            addVH(k,stick,va);
        }
        else if(hyb == 2 && nhy == 1 && val == 2) {
            vv=new double[3];
            xx=0.5;
            yy=0.86602540375;
            a=home.config.structure.bond[0][k];
            va=BML.ndxyz(k,a,home.config.xyz);
            atn=Math.atan2(va[1],va[0]);
            cp=Math.cos(atn);
            sp=Math.sin(atn);
            vv[0]=va[0]*xx-sp*yy;
            vv[1]=va[1]*xx+cp*yy;
            vv[2]=va[2]*xx;
            addVH(k,stick,vv);
        }
        else if(hyb == 2 && nhy == 1 && val == 3) {
            vv=new double[3];
            a=home.config.structure.bond[0][k];
            b=home.config.structure.bond[1][k];
            va=BML.ndxyz(k,a,home.config.xyz);
            vb=BML.ndxyz(k,b,home.config.xyz);
            vv[0]=va[0]+vb[0];
            vv[1]=va[1]+vb[1];
            vv[2]=va[2]+vb[2];
            BML.vnorm(vv);
            addVH(k,stick,vv);
        }
        else if(hyb == 2 && nhy == 2 && val == 2) {
            vv=new double[3];
            vv[0]=0.5;
            vv[1]=0.86602540384;
            vv[2]=0.0;
            addVH(k,stick,vv);
            vv[1]=-0.86602540384;
            addVH(k,stick,vv);
        }
        else if(hyb == 2 && nhy == 2 && val == 3) {
            vv=new double[3];
            xx=0.5;
            yy=0.86602540375;
            a=home.config.structure.bond[0][k];
            va=BML.ndxyz(k,a,home.config.xyz);
            atn=Math.atan2(va[1],va[0]);
            cp=Math.cos(atn);
            sp=Math.sin(atn);
            vv[0]=va[0]*xx-sp*yy;
            vv[1]=va[1]*xx+cp*yy;
            vv[2]=va[2]*xx;
            addVH(k,stick,vv);
            vv[0]=va[0]*xx+sp*yy;
            vv[1]=va[1]*xx-cp*yy;
            vv[2]=va[2]*xx;
            addVH(k,stick,vv);
        }
        else if(hyb == 2 && nhy == 3 && val == 3) {
            vv=new double[3];
            vv[0]=1.0;
            vv[1]=0.0;
            vv[2]=0.0;
            addVH(k,stick,vv);
            vv[0]=0.5;
            vv[1]=0.86602540384;
            addVH(k,stick,vv);
            vv[2]=-0.86602540384;
            addVH(k,stick,vv);
        }
        else if(hyb == 3 && nhy == 1 && val == 2) {
            vv=new double[3];
            a=home.config.structure.bond[0][k];
            va=BML.ndxyz(k,a,home.config.xyz);
            atn=Math.atan2(va[1],va[0]);
            cp=Math.cos(atn);
            sp=Math.sin(atn);
            xx=1.0/3.0;
            yy=0.94280904165;
            vv[0]=xx*va[0]-yy*sp;
            vv[1]=xx*va[1]+yy*cp;
            vv[2]=xx*va[2];
            addVH(k,stick,vv);
        }
        else if(hyb == 3 && nhy == 1 && val == 3) {
            vv=new double[3];
            vc=new double[3];
            a=home.config.structure.bond[0][k];
            b=home.config.structure.bond[1][k];
            va=BML.ndxyz(k,a,home.config.xyz);
            vb=BML.ndxyz(k,b,home.config.xyz);
            vc[0]=va[0]+vb[0];
            vc[1]=va[1]+vb[1];
            vc[2]=va[2]+vb[2];
            BML.vnorm(vc);
            vd=BML.cross(va,vb);
            BML.vnorm(vd);
            xx=0.57735026924;
            yy=0.81649658092;
            vv[0]=xx*vc[0]+yy*vd[0];
            vv[1]=xx*vc[1]+yy*vd[1];
            vv[2]=xx*vc[2]+yy*vd[2];
            addVH(k,stick,vv);
        }
        else if(hyb == 3 && nhy == 1 &&  val == 4) {
            vv=new double[3];
            a=home.config.structure.bond[0][k];
            b=home.config.structure.bond[1][k];
            c=home.config.structure.bond[2][k];
            va=BML.ndxyz(k,a,home.config.xyz);
            vb=BML.ndxyz(k,b,home.config.xyz);
            vc=BML.ndxyz(k,c,home.config.xyz);
            vv[0]=va[0]+vb[0]+vc[0];
            vv[1]=va[1]+vb[1]+vc[1];
            vv[2]=va[2]+vb[2]+vc[2];
            BML.vnorm(vv);
            addVH(k,stick,vv);
        }
        else if(hyb == 3 && nhy == 2 && val == 2) {
            vv=new double[3];
            vv[0]=1.0/3.0;
            vv[1]=0.81649658092;
            vv[2]=-0.47140452082;
            addVH(k,stick,vv);
            vv[1]=-0.81649658092;
            addVH(k,stick,vv);
        }
        else if(hyb == 3 && nhy == 2 && val == 3) {
            vv=new double[3];
            a=home.config.structure.bond[0][k];
            va=BML.ndxyz(k,a,home.config.xyz);
            atn=Math.atan2(va[1],va[0]);
            cp=Math.cos(atn);
            sp=Math.sin(atn);
            st=va[2];
            ct=Math.sqrt(va[0]*va[0]+va[1]*va[1]);
            xx=1.0/3.0;
            yy=0.81649658092;
            zz=-0.47140452082;
            vv[0]=xx*va[0]-yy*sp-zz*cp*st;
            vv[1]=xx*va[1]+yy*cp-zz*sp*st;
            vv[2]=xx*st+zz*ct;
            addVH(k,stick,vv);
            vv[0]=xx*va[0]+yy*sp-zz*cp*st;
            vv[1]=xx*va[1]-yy*cp-zz*sp*st;
            vv[2]=xx*st+zz*ct;
            addVH(k,stick,vv);
        }
        else if(hyb == 3 && nhy == 2 && val == 4) {
            vv=new double[3];
            vc=new double[3];
            a=home.config.structure.bond[0][k];
            b=home.config.structure.bond[1][k];
            va=BML.ndxyz(k,a,home.config.xyz);
            vb=BML.ndxyz(k,b,home.config.xyz);
            vc[0]=va[0]+vb[0];
            vc[1]=va[1]+vb[1];
            vc[2]=va[2]+vb[2];
            BML.vnorm(vc);
            vd=BML.cross(va,vb);
            BML.vnorm(vd);
            xx=0.57735026924;
            yy=0.81649658092;
            vv[0]=xx*vc[0]+yy*vd[0];
            vv[1]=xx*vc[1]+yy*vd[1];
            vv[2]=xx*vc[2]+yy*vd[2];
            addVH(k,stick,vv);
            vv[0]=xx*vc[0]-yy*vd[0];
            vv[1]=xx*vc[1]-yy*vd[1];
            vv[2]=xx*vc[2]-yy*vd[2];
            addVH(k,stick,vv);
        }
        else if(hyb == 3 && nhy == 3 && val == 3) {
            vv=new double[3];
            vv[0]=1.0/3.0;
            vv[1]=0.0;
            vv[2]=0.94280904165;
            addVH(k,stick,vv);
            vv[2]=-0.47140452076;
            vv[1]=0.81649658092;
            addVH(k,stick,vv);
            vv[1]=-0.81649658092;
            addVH(k,stick,vv);
        }
        else if(hyb == 3 && nhy == 3 && val == 4) {
            vv=new double[3];
            a=home.config.structure.bond[0][k];
            va=BML.ndxyz(k,a,home.config.xyz);
            atn=Math.atan2(va[1],va[0]);
            cp=Math.cos(atn);
            sp=Math.sin(atn);
            ct=va[2];
            st=Math.sqrt(va[0]*va[0]+va[1]*va[1]);
            xx=0.94280904165;
            yy=0.0;
            zz=1.0/3.0;
            vv[0]=ct*cp*xx-sp*yy+cp*st*zz;
            vv[1]=ct*sp*xx+cp*yy+sp*st*zz;
            vv[2]=ct*zz-st*xx;
            addVH(k,stick,vv);
            xx=-0.47140452076;
            yy=0.81649658092;
            vv[0]=ct*cp*xx-sp*yy+cp*st*zz;
            vv[1]=ct*sp*xx+cp*yy+sp*st*zz;
            vv[2]=ct*zz-st*xx;
            addVH(k,stick,vv);
            yy=-0.81649658092;
            vv[0]=ct*cp*xx-sp*yy+cp*st*zz;
            vv[1]=ct*sp*xx+cp*yy+sp*st*zz;
            vv[2]=ct*zz-st*xx;
            addVH(k,stick,vv);
        }
        else if(hyb == 3 && nhy == 4 && val == 4) {
            vv=new double[3];
            vv[0]=0.57735026918;
            vv[1]=0.57735026918;
            vv[2]=0.57735026918;
            addVH(k,stick,vv);
            vv[0]=-0.57735026918;
            vv[1]=-0.57735026918;
            vv[2]=0.57735026918;
            addVH(k,stick,vv);
            vv[0]=0.57735026918;
            vv[1]=-0.57735026918;
            vv[2]=-0.57735026918;
            addVH(k,stick,vv);
            vv[0]=-0.57735026918;
            vv[1]=0.57735026918;
            vv[2]=-0.57735026918;
            addVH(k,stick,vv);
        }
    }

    void addVH(int k,double aaa,double vvv[]) {
        /*
*********************************************************************

dl_poly/java GUI routine to add hydrogen atom and bond vector

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
         */

        if (home.config.atoms[k].zsym.indexOf("OW") >=0 )
            home.config.atoms[home.config.natms]=new Element("HW");
        else
            home.config.atoms[home.config.natms]=new Element("H_");
        home.config.xyz[0][home.config.natms]=home.config.xyz[0][k]+aaa*vvv[0];
        home.config.xyz[1][home.config.natms]=home.config.xyz[1][k]+aaa*vvv[1];
        home.config.xyz[2][home.config.natms]=home.config.xyz[2][k]+aaa*vvv[2];
        home.config.structure.join[0][home.config.structure.nbnds]=Math.min(k,home.config.natms);
        home.config.structure.join[1][home.config.structure.nbnds]=Math.max(k,home.config.natms);
        home.config.structure.bond[home.config.structure.lbnd[k]][k]=home.config.natms;
        home.config.structure.bond[0][home.config.natms]=k;
        home.config.structure.lbnd[home.config.natms]=1;
        home.config.structure.lbnd[k]++;
        home.config.structure.nbnds++;
        home.config.natms++;

    }

    void delHAtoms(int k) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
         */
        int n,m=0;

        if(ngroup == 0) {
            if(k >= 0) {
                m=0;
                int zapp[]=new int[nhdel];
                for (int i=0;i<home.config.structure.lbnd[k];i++) {
                    if(home.config.atoms[home.config.structure.bond[i][k]].znum == 1) {
                        zapp[m]=home.config.structure.bond[i][k];
                        m++;
                    }
                }
                for(int i=0;i<nhdel;i++) {
                    for(int j=i+1;j<nhdel;j++) {
                        n=zapp[j];
                        if(zapp[j] < zapp[i]) {
                            n=zapp[i];
                            zapp[i]=zapp[j];
                        }
                        zapp[j]=n-1;
                    }
                    deleteAtom(zapp[i]);
                }
            }
            else {
                for (int i=0;i<nhdel;i++) {
                    OUT:
                    for(int j=m;j<home.config.natms;j++) {
                        if(home.config.atoms[j].znum == 1) {
                            deleteAtom(j);
                            m=j;
                            break OUT;
                        }
                    }
                }
            }
        }
        else {
            for (int i=0;i<nhdel;i++) {
                OUT:
                    for(int j=m;j<home.config.natms;j++) {
                        if(group[j] > 0 && home.config.atoms[j].znum == 1) {
                            deleteAtom(j);
                            m=j;
                            break OUT;
                        }
                    }
            }
        }
        home.safe=false;
        repaint();
    }

    void resizeGroupArray() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        int[] ggg=new int[group.length+MXGROUP];
        for(int i=0;i<ggg.length;i++) {
            if(i < group.length)
                ggg[i]=group[i];
            else
                ggg[i]=-1;

        }
        group=ggg;
    }

    boolean duplicateBond() {
        /*
*********************************************************************

dl_poly/java GUI routine to remove duplicated bonds

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
         */
        int n,m,k;

        k=-1;
        n=Math.min(mark0,mark1);
        m=Math.max(mark0,mark1);

        OUT:
            for(int i=0;i<home.config.structure.nbnds;i++) {
                if(n == home.config.structure.join[0][i] && m == home.config.structure.join[1][i]) {
                    k=i;
                    home.config.structure.join[0][i]=-1;
                    break OUT;
                }
            }
            if(k < 0) return false;

            // Delete duplicated bond

            home.config.structure.nbnds--;
            for(int i=k;i<home.config.structure.nbnds;i++) {
                home.config.structure.join[0][i]=home.config.structure.join[0][i+1];
                home.config.structure.join[1][i]=home.config.structure.join[1][i+1];
            }
            mark0=-1;
            mark1=-1;

            // Rebuild connection table

            rebuildConnections();
            repaint();

            return true;
    }

    void deleteAtom(int k) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
         */
        int j,m,n;

        // Remove atom

        home.config.natms--;
        home.safe=false;
        for(int i=k;i<home.config.natms;i++) {
            home.config.atoms[i]=new Element(home.config.atoms[i+1].zsym);
            home.config.xyz[0][i]=home.config.xyz[0][i+1];
            home.config.xyz[1][i]=home.config.xyz[1][i+1];
            home.config.xyz[2][i]=home.config.xyz[2][i+1];
            if(ngroup > 0)group[i]=group[i+1];
        }

        // Delete cancelled bonds

        j=0;
        for(int i=0;i<home.config.structure.nbnds;i++) {
            m=home.config.structure.join[0][i];
            n=home.config.structure.join[1][i];
            if(m != k && n != k) {
                if(m > k) m--;
                if(n > k) n--;
                home.config.structure.join[0][j]=m;
                home.config.structure.join[1][j]=n;
                j++;
            }
        }
        home.config.structure.nbnds=j;

        // Rebuild connection table

        rebuildConnections();
    }

    void deleteGroup() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
         */
        int k;
        int n=0;

        if(oper == 3) {
            if(ngroup == 0) {
                // Get identity of clicked atom

                k=getAtom();
                if(k >= 0) deleteAtom(k);
            }
            else {
                for (int i=0;i<ngroup;i++) {
                    OUT:
                        for(int j=n;j<home.config.natms;j++) {
                            if(group[j] > 0) {
                                deleteAtom(j);
                                n=j;
                                break OUT;
                            }
                        }
                }
                ngroup=0;
            }
            home.safe=false;
            repaint();
        }
    }

    void duplicateGroup() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
         */
        int k,n,m;
        int iab[];

        if(oper == 11 && ngroup > 0) {
            m=home.config.natms+ngroup;
            while(m > group.length)
                resizeGroupArray();
            while(m > home.config.atoms.length) {
                home.config.resizeArrays();
                home.config.structure.resizeBondArrays();
            }
            iab=new int[home.config.atoms.length];
            k=home.config.natms;
            for (int i=0;i<home.config.natms;i++) {
                iab[i]=-1;
                if(group[i] > 0) {
                    iab[i]=k;
                    home.config.atoms[k]=new Element(home.config.atoms[i].zsym);
                    home.config.xyz[0][k]=home.config.xyz[0][i]+home.incx;
                    home.config.xyz[1][k]=home.config.xyz[1][i];
                    home.config.xyz[2][k]=home.config.xyz[2][i]+home.incz;
                    group[k]=1;
                    k++;
                }
            }
            k=home.config.structure.nbnds;
            for (int i=0;i<home.config.structure.nbnds;i++) {
                n=home.config.structure.join[0][i];
                m=home.config.structure.join[1][i];
                if(group[n] > 0 && group[m] > 0) {
                    n=iab[n];
                    m=iab[m];
                    if(k == home.config.structure.join[0].length)
                        home.config.structure.resizeJoinArray();
                    home.config.structure.join[0][k]=n;
                    home.config.structure.join[1][k]=m;
                    home.config.structure.bond[home.config.structure.lbnd[n]][n]=m;
                    home.config.structure.bond[home.config.structure.lbnd[m]][m]=n;
                    home.config.structure.lbnd[n]++;
                    home.config.structure.lbnd[m]++;
                    k++;
                }
            }
            home.config.structure.nbnds=k;
            for(int i=0;i<home.config.natms;i++)
                group[i]=-1;
            home.config.natms+=ngroup;
            home.safe=false;
            oper=6;
            news="SHIFT";
        }
        repaint();
    }

    void insertFragment() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
         */
        int k,n,m,base;

        base=home.config.natms;

        if(fragment != null) {
            if(home.config == null) home.config=fragment;

            // ensure arrays are big enough

            if(ngroup == 0)
                group=new int[MXGROUP];
            m=home.config.natms+fragment.natms;
            while(m > group.length)
                resizeGroupArray();
            while(m > home.config.atoms.length) {
                home.config.resizeArrays();
                home.config.structure.resizeBondArrays();
            }
            m=home.config.structure.nbnds+fragment.structure.nbnds;
            while(m > home.config.structure.join[0].length)
                home.config.structure.resizeJoinArray();

            // add fragment atoms

            k=home.config.natms;
            for (int i=0;i<home.config.natms;i++)
                group[i]=-1;
            for (int i=0;i<fragment.natms;i++) {
                home.config.atoms[k]=new Element(fragment.atoms[i].zsym);
                home.config.xyz[0][k]=fragment.xyz[0][i];
                home.config.xyz[1][k]=fragment.xyz[1][i];
                home.config.xyz[2][k]=fragment.xyz[2][i];
                group[k]=1;
                k++;
            }
            home.config.natms=k;

            // add fragment bonds

            k=home.config.structure.nbnds;
            for (int i=0;i<fragment.structure.nbnds;i++) {
                n=fragment.structure.join[0][i]+base;
                m=fragment.structure.join[1][i]+base;
                home.config.structure.join[0][k]=n;
                home.config.structure.join[1][k]=m;
                home.config.structure.bond[home.config.structure.lbnd[n]][n]=m;
                home.config.structure.bond[home.config.structure.lbnd[m]][m]=n;
                home.config.structure.lbnd[n]++;
                home.config.structure.lbnd[m]++;
                k++;
            }
            home.config.structure.nbnds=k;
            ngroup=fragment.natms;
            home.safe=false;
            oper=6;
            news="SHIFT";
        }
        repaint();
    }

    void rebuildConnections() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
         */
        int j,k;

        for(int i=0;i<home.config.structure.lbnd.length;i++)
            home.config.structure.lbnd[i]=0;

        //if(home.config.structure.nbnds > 0) {
            for(int i=0;i<home.config.structure.nbnds;i++) {
                j=home.config.structure.join[0][i];
                k=home.config.structure.join[1][i];
                home.config.structure.bond[home.config.structure.lbnd[j]][j]=k;
                home.config.structure.bond[home.config.structure.lbnd[k]][k]=j;
                home.config.structure.lbnd[j]++;
                home.config.structure.lbnd[k]++;
            }
            //}
    }

    void shiftMolecule(String activity) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
         */
        int k;
        double rad;

        k=getAtom();

        if(ngroup == 0) {
            if(k >= 0) {
                if(activity.equals("Tz+"))
                    home.config.xyz[2][k]+=0.1;
                else if(activity.equals("Tz-"))
                    home.config.xyz[2][k]-=0.1;
                else if(activity.equals("Tx+"))
                    home.config.xyz[0][k]+=0.1;
                else if(activity.equals("Tx-"))
                    home.config.xyz[0][k]-=0.1;
                else if(activity.equals("Ty+"))
                    home.config.xyz[1][k]+=0.1;
                else if(activity.equals("Ty-"))
                    home.config.xyz[1][k]-=0.1;
            }
            else {
                for(int i=0;i<home.config.natms;i++) {
                    if(activity.equals("Tz+"))
                        home.config.xyz[2][i]+=0.1;
                    else if(activity.equals("Tz-"))
                        home.config.xyz[2][i]-=0.1;
                    else if(activity.equals("Tx+"))
                        home.config.xyz[0][i]+=0.1;
                    else if(activity.equals("Tx-"))
                        home.config.xyz[0][i]-=0.1;
                    else if(activity.equals("Ty+"))
                        home.config.xyz[1][i]+=0.1;
                    else if(activity.equals("Ty-"))
                        home.config.xyz[1][i]-=0.1;
                }
            }
        }
        else {
            if(k >= 0 && group[k] < 0) {
                if(activity.equals("Tz+"))
                    home.config.xyz[2][k]+=0.1;
                else if(activity.equals("Tz-"))
                    home.config.xyz[2][k]-=0.1;
                else if(activity.equals("Tx+"))
                    home.config.xyz[0][k]+=0.1;
                else if(activity.equals("Tx-"))
                    home.config.xyz[0][k]-=0.1;
                else if(activity.equals("Ty+"))
                    home.config.xyz[1][k]+=0.1;
                else if(activity.equals("Ty-"))
                    home.config.xyz[1][k]-=0.1;
            }
            else {
                for (int i=0;i<home.config.natms;i++) {
                    if(group[i] > 0) {
                        if(activity.equals("Tz+"))
                            home.config.xyz[2][i]+=0.1;
                        else if(activity.equals("Tz-"))
                            home.config.xyz[2][i]-=0.1;
                        else if(activity.equals("Tx+"))
                            home.config.xyz[0][i]+=0.1;
                        else if(activity.equals("Tx-"))
                            home.config.xyz[0][i]-=0.1;
                        else if(activity.equals("Ty+"))
                            home.config.xyz[1][i]+=0.1;
                        else if(activity.equals("Ty-"))
                            home.config.xyz[1][i]-=0.1;
                    }
                }
            }
        }
        home.safe=false;
        repaint();
    }

    void moveMolecule() {
        /*
*********************************************************************

dl_Poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
         */
        int k;
        double[] dd=new double[2];
        double dx,dz,rad;

        dx=(fx-sx)/scale;
        dz=(sy-fy)/scale;

        if(activity.equals("Tx+") || activity.equals("Tx-")) {
            dd[0]=dx;
            dd[1]=0;
        }
        else if(activity.equals("Tz+") || activity.equals("Tz-")) {
            dd[0]=0;
            dd[1]=dz;
        }
        else {
            dd[0]=dx;
            dd[1]=dz;
        }

        if(ngroup == 0) {
            sx=sx-tilex/2;
            sy=tiley/2-sy;
            k=getAtom();
            if(k >= 0) {
                home.config.xyz[0][k]+=dd[0];
                home.config.xyz[2][k]+=dd[1];
            }
            else {
                for(int i=0;i<home.config.natms;i++) {
                    home.config.xyz[0][i]+=dd[0];
                    home.config.xyz[2][i]+=dd[1];
                }
            }
        }
        else {
            for(int i=0;i<home.config.natms;i++) {
                if(group[i] > 0) {
                    home.config.xyz[0][i]+=dd[0];
                    home.config.xyz[2][i]+=dd[1];
                }
            }
        }
        home.safe=false;
        repaint();
    }

    void turnMolecule(String activity) {
        /*
*********************************************************************

dl_Poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
         */
        int n,j,k;
        double u,v,ca,sa,rad;

        j=-1;
        k=-1;
        sa=home.rotsin;
        ca=home.rotcos;

        if(activity.equals("Rx+")) {
            j=1;
            k=2;
            sa=home.rotsin;
        }
        else if(activity.equals("Rx-")) {
            j=1;
            k=2;
            sa=-home.rotsin;
        }
        else if(activity.equals("Ry+")) {
            j=2;
            k=0;
            sa=home.rotsin;
        }
        else if(activity.equals("Ry-")) {
            j=2;
            k=0;
            sa=-home.rotsin;
        }
        else if(activity.equals("Rz+")) {
            j=0;
            k=1;
            sa=home.rotsin;
        }
        else if(activity.equals("Rz-")) {
            j=0;
            k=1;
            sa=-home.rotsin;
        }

        if(ngroup == 0) {
            n=getAtom();
            if(n >= 0) {
                u=home.config.xyz[j][n];
                v=home.config.xyz[k][n];
                home.config.xyz[j][n]=(ca*u-sa*v);
                home.config.xyz[k][n]=(sa*u+ca*v);
            }
            else {
                for (int i=0;i<home.config.natms;i++) {
                    u=home.config.xyz[j][i];
                    v=home.config.xyz[k][i];
                    home.config.xyz[j][i]=(ca*u-sa*v);
                    home.config.xyz[k][i]=(sa*u+ca*v);
                }
            }
        }
        else {
            for (int i=0;i<home.config.natms;i++) {
                if(group[i] > 0) {
                    u=home.config.xyz[j][i];
                    v=home.config.xyz[k][i];
                    home.config.xyz[j][i]=(ca*u-sa*v);
                    home.config.xyz[k][i]=(sa*u+ca*v);
                }
            }
        }
        home.safe=false;
        repaint();
    }

    void rotateMolecule() {
        /*
*********************************************************************

dl_Poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
         */
        double[] r={1,0,0,0,1,0,0,0,1};
        double rad,r1,r2,ca,cb,sa,sb,tx,ty,tz;

        if(activity.equals("Rx+") || activity.equals("Rx-")) {
            ca=Math.cos((fy-sy)*Math.PI/tiley);
            sa=Math.sin((fy-sy)*Math.PI/tiley);
            r[4]= ca;
            r[5]= sa;
            r[7]=-sa;
            r[8]= ca;
        }
        else if(activity.equals("Rz+") || activity.equals("Rz-")) {
            cb=Math.cos((fx-sx)*Math.PI/tilex);
            sb=Math.sin((fx-sx)*Math.PI/tilex);
            r[0]= cb;
            r[1]= sb;
            r[3]=-sb;
            r[4]= cb;
        }
        else if(activity.equals("Ry+") || activity.equals("Ry-")) {
            if(sx > tilex/2)
                rad=(fy-sy)*Math.PI/tiley;
            else
                rad=-(fy-sy)*Math.PI/tiley;
            ca=Math.cos(rad);
            sa=Math.sin(rad);
            r[0]= ca;
            r[2]=-sa;
            r[6]= sa;
            r[8]= ca;
        }
        else {
            ca=Math.cos((sx-fx)*Math.PI/tilex);
            sa=Math.sin((sx-fx)*Math.PI/tilex);
            cb=Math.cos((sy-fy)*Math.PI/tiley);
            sb=Math.sin((sy-fy)*Math.PI/tiley);
            r[0]= ca;
            r[2]= sa*sb;
            r[1]=-sa*cb;
            r[6]= 0;
            r[8]= cb;
            r[7]= sb;
            r[3]= sa;
            r[5]=-ca*sb;
            r[4]= ca*cb;
        }
        if(ngroup > 0) {
            for (int i=0;i<home.config.natms;i++) {
                if(group[i] > 0) {
                    tx=home.config.xyz[0][i];
                    ty=home.config.xyz[1][i];
                    tz=home.config.xyz[2][i];
                    home.config.xyz[0][i]=r[0]*tx+r[3]*ty+r[6]*tz;
                    home.config.xyz[1][i]=r[1]*tx+r[4]*ty+r[7]*tz;
                    home.config.xyz[2][i]=r[2]*tx+r[5]*ty+r[8]*tz;
                }
            }
        }
        else {
            for (int i=0;i<home.config.natms;i++) {
                tx=home.config.xyz[0][i];
                ty=home.config.xyz[1][i];
                tz=home.config.xyz[2][i];
                home.config.xyz[0][i]=r[0]*tx+r[3]*ty+r[6]*tz;
                home.config.xyz[1][i]=r[1]*tx+r[4]*ty+r[7]*tz;
                home.config.xyz[2][i]=r[2]*tx+r[5]*ty+r[8]*tz;
            }
        }
        home.safe=false;
        repaint();
    }

    void Optimize() {
        /*
*********************************************************************

dl_Poly/java GUI routine to optimize structures

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
         */

        double step;
        int k,matms,mbnds,mxangs;
        double[] fff,grad,hnrm,length,angcon;
        double[][] hhh,dxyz;
        int[][] jopt,tbond,angle;
        int[] ida,idb,tlbnd;

        if(home.config.natms < 2) return;

        if(oper == 9) {
            if(ngroup == 0) {
                jopt=home.config.structure.join;
                tbond=home.config.structure.bond;
                tlbnd=home.config.structure.lbnd;
                mbnds=home.config.structure.nbnds;
                matms=home.config.natms;
                ida=new int[home.config.natms];
                idb=new int[home.config.natms];
                for(int i=0;i<home.config.natms;i++) {
                    ida[i]=i;
                    idb[i]=i;
                }
            }
            else {
                mbnds=0;
                matms=ngroup;
                jopt=new int[2][home.config.structure.nbnds];
                for (int i=0;i<home.config.structure.nbnds;i++) {
                    if(group[home.config.structure.join[0][i]] > 0 || group[home.config.structure.join[1][i]] > 0) {
                        jopt[0][mbnds]=home.config.structure.join[0][i];
                        jopt[1][mbnds]=home.config.structure.join[1][i];
                        mbnds++;
                    }
                }
                k=0;
                ida=new int[matms];
                idb=new int[home.config.natms];
                tlbnd=new int[matms];
                tbond=new int[Basic.MXCONNECT][matms];
                for (int i=0;i<home.config.natms;i++) {
                    idb[i]=-1;
                    if(group[i] > 0) {
                        ida[k]=i;
                        idb[i]=k;
                        tlbnd[k]=home.config.structure.lbnd[i];
                        for (int j=0;j<home.config.structure.lbnd[i];j++) {
                            tbond[j][k]=home.config.structure.bond[j][i];
                        }
                        k++;
                    }
                }
            }
            keyopt=0;
            step=0.01;
            fff=new double[3];
            grad=new double[3];
            hnrm=new double[1];
            hhh=new double[3][matms];
            dxyz=new double[3][matms];
            length=new double[mbnds];
            mxangs=matms*12;
            angle=new int[3][mxangs];
            angcon=new double[mxangs];


            // Define bond properties

            bondProperties(matms,mbnds,mxangs,ida,jopt,angle,tlbnd,tbond,length,angcon);

            k=0;
            while (k < 1000 && keyopt != 999) {
                fff[0]=energyGradient(matms,mbnds,idb,jopt,tlbnd,tbond,angle,
                length,angcon,dxyz);
                keyopt=cgropt(matms,step,ida,hnrm,fff,grad,hhh,dxyz);
                k++;
            }
            if(keyopt == 999)
                news="OPTIMIZATION: CONVERGED";
            else
                news="OPTIMIZATION: ENERGY ="+BML.fmt(fff[0],12);
            home.safe=false;
            repaint();
        }
    }

    void bondProperties(int natms,int nbnds,int mxangs,int ida[],int join[][],int angle[][],
    int lbnd[],int bond[][],double length[],double angcon[]) {
        /*
*********************************************************************

dl_Poly/java GUI routine to define bond and angle properties

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
         */
        int n,m;

        // Define target bond lengths

        for(int i=0;i<nbnds;i++) {
            n=join[0][i];
            m=join[1][i];
            length[i]=home.config.atoms[n].zrad+home.config.atoms[m].zrad;
        }

        // Define target bond angles

        n=0;
        for(int i=0;i<natms;i++) {
            if(lbnd[i] > 1) {
                for(int j=1;j<lbnd[i];j++) {
                    for(int k=0;k<j;k++) {
                        m=ida[i];
                        angle[0][n]=bond[j][i];
                        angle[1][n]=m;
                        angle[2][n]=bond[k][i];
                        if(home.config.atoms[m].zsym.charAt(1) == 'W')
                            angcon[n]=-1.0/3.0;
                        else if(home.config.atoms[m].zsym.charAt(2) == '3')
                            angcon[n]=-1.0/3.0;
                        else if(home.config.atoms[m].zsym.charAt(2) == '2')
                            angcon[n]=-1.0/2.0;
                        else if(home.config.atoms[m].zsym.charAt(2) == 'R')
                            angcon[n]=-1.0/2.0;
                        else if(home.config.atoms[m].zsym.charAt(2) == '1')
                            angcon[n]=-1.0;
                        n++;
                        if(n == mxangs) {
                            mxangs*=2;
                            int aaa[][]=new int[3][mxangs];
                            double bbb[]=new double[mxangs];
                            for(int o=0;o<n;o++) {
                                bbb[o]=angcon[o];
                                aaa[0][o]=angle[0][o];
                                aaa[1][o]=angle[1][o];
                                aaa[2][o]=angle[2][o];
                            }
                            angle=aaa;
                            angcon=bbb;
                        }
                    }
                }
            }
        }
        nangs=n;
    }

    double energyGradient(int natms,int nbnds,int idb[],int join[][],int lbnd[],int bond[][],
    int angle[][],double length[],double angcon[],double dxyz[][]) {
        /*
*********************************************************************

dl_Poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
         */
        int k,l,n,m,q,b,c;
        double energy,harm,hcos,dx,dy,dz,dxa,dya,dza,dxb,dyb,dzb;
        double rr,ra,rb,coz,gamma,gxa,gya,gza,gxb,gyb,gzb,rt6,sig;
        double sig2,rsig,rs3,rs6,rrr,rsq,gmax;

        harm=1.0;
        hcos=0.5;
        gmax=0.0;
        energy=0.0;

        // Initialise the gradient array

        for(int i=0;i<natms;i++) {
            dxyz[0][i]=0.0;
            dxyz[1][i]=0.0;
            dxyz[2][i]=0.0;
        }

        // Contribution from harmonic bonds

        for(int i=0;i<nbnds;i++) {
            n=join[0][i];
            m=join[1][i];
            dx=home.config.xyz[0][m]-home.config.xyz[0][n];
            dy=home.config.xyz[1][m]-home.config.xyz[1][n];
            dz=home.config.xyz[2][m]-home.config.xyz[2][n];
            rr=Math.sqrt(dx*dx+dy*dy+dz*dz);
            energy+=0.5*harm*Math.pow((rr-length[i]),2);
            gamma=-harm*(rr-length[i])/rr;
            gmax=Math.max(gmax,Math.abs(gamma));
            m=idb[m];
            if(m >= 0) {
                dxyz[0][m]+=(dx*gamma);
                dxyz[1][m]+=(dy*gamma);
                dxyz[2][m]+=(dz*gamma);
            }
            n=idb[n];
            if(n >= 0) {
                dxyz[0][n]-=(dx*gamma);
                dxyz[1][n]-=(dy*gamma);
                dxyz[2][n]-=(dz*gamma);
            }
        }

        // Contribution from cosine angle potentials

        for(int i=0;i<nangs;i++) {
            k=angle[0][i];
            n=angle[1][i];
            m=angle[2][i];
            dxa=home.config.xyz[0][k]-home.config.xyz[0][n];
            dya=home.config.xyz[1][k]-home.config.xyz[1][n];
            dza=home.config.xyz[2][k]-home.config.xyz[2][n];
            ra=Math.sqrt(dxa*dxa+dya*dya+dza*dza);
            dxa/=ra;
            dya/=ra;
            dza/=ra;
            dxb=home.config.xyz[0][m]-home.config.xyz[0][n];
            dyb=home.config.xyz[1][m]-home.config.xyz[1][n];
            dzb=home.config.xyz[2][m]-home.config.xyz[2][n];
            rb=Math.sqrt(dxb*dxb+dyb*dyb+dzb*dzb);
            dxb/=rb;
            dyb/=rb;
            dzb/=rb;
            coz=dxa*dxb+dya*dyb+dza*dzb;
            energy+=(0.5*hcos*Math.pow((coz-angcon[i]),2));
            gamma=-hcos*(coz-angcon[i]);
            gmax=Math.max(gmax,Math.abs(gamma));
            gxa=gamma*(dxb-dxa*coz)/ra;
            gya=gamma*(dyb-dya*coz)/ra;
            gza=gamma*(dzb-dza*coz)/ra;
            gxb=gamma*(dxa-dxb*coz)/rb;
            gyb=gamma*(dya-dyb*coz)/rb;
            gzb=gamma*(dza-dzb*coz)/rb;
            k=idb[k];
            if(k >= 0) {
                dxyz[0][k]+=gxa;
                dxyz[1][k]+=gya;
                dxyz[2][k]+=gza;
            }
            n=idb[n];
            if(n >= 0) {
                dxyz[0][n]-=(gxa+gxb);
                dxyz[1][n]-=(gya+gyb);
                dxyz[2][n]-=(gza+gzb);
            }
            m=idb[m];
            if(m >= 0) {
                dxyz[0][m]+=gxb;
                dxyz[1][m]+=gyb;
                dxyz[2][m]+=gzb;
            }

        }

        return energy;
    }

    void defineBox() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        if(home.config != null) {
            if(home.config.pbc.imcon != boxtyp) {
                home.config.pbc.imcon=boxtyp;
                home.config.pbc.cell[0]=10;
                home.config.pbc.cell[4]=10;
                home.config.pbc.cell[8]=10;
                if(boxtyp == 5)
                    home.config.pbc.cell[8]=10*Math.sqrt(2.0);
                home.config.pbc.buildBoundary(home.config.pbc.imcon);
                home.safe=false;
            }
            scale=setBoxScale();
        }
    }

    void boxResize() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        double scy=Math.max(0.1,1+(fy-sy)/(double)tilex);
        if(home.config != null) {
            if(home.config.pbc.imcon == 2 || home.config.pbc.imcon == 7) {
                double scx=Math.max(0.1,1+(fx-sx)/(double)tilex);
                home.config.pbc.cell[0]*=scx;
                home.config.pbc.cell[4]*=scx;
                home.config.pbc.cell[8]*=scy;
                for (int i=0;i<home.config.pbc.num_vertices;i++) {
                    home.config.pbc.vrt[0][i]*=scx;
                    home.config.pbc.vrt[1][i]*=scx;
                    home.config.pbc.vrt[2][i]*=scy;
                }
            }
            else {
                home.config.pbc.cell[0]*=scy;
                home.config.pbc.cell[4]*=scy;
                home.config.pbc.cell[8]*=scy;
                for (int i=0;i<home.config.pbc.num_vertices;i++) {
                    home.config.pbc.vrt[0][i]*=scy;
                    home.config.pbc.vrt[1][i]*=scy;
                    home.config.pbc.vrt[2][i]*=scy;
                }
            }
            pen++;
            home.safe=false;
            if(pen == 5) {
                pen=0;
                sx=fx;
                sy=fy;
                scale=setBoxScale();
                repaint();
            }
        }
    }

    class MyScreenSizer implements ComponentListener {
        /*
*********************************************************************

dl_poly/java GUI class to handle screen events

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        public void componentResized(ComponentEvent e) {
            /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
             */
            Dimension arg = ((Component)e.getSource()).getSize();

            tilex = arg.width;
            tiley = arg.height;
            if(home.config!=null) {
                if(home.config.pbc.imcon > 0)
                    scale=setBoxScale();
                else
                    scale=37.5;
                repaint();
            }
        }
        public void componentMoved(ComponentEvent e){}
        public void componentShown(ComponentEvent e){}
        public void componentHidden(ComponentEvent e){}
    }

    class MousePoints implements MouseListener {
        /*
*********************************************************************

dl_poly/java GUI class to handle mouse events

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        public void mouseClicked(MouseEvent e) {
            /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
             */
            sx=e.getX()-tilex/2;
            sy=tiley/2-e.getY();

            if(home.edit) {
                if(oper == 0) {
                    identifyAtom();
                }
                else if(oper == 1) {
                    moleculeBuilder();
                }
                else if(oper == 2) {
                    linkMaker();
                }
                else if(oper == 3) {
                    deleteGroup();
                }
                else if(oper/10 == 4) {
                    shiftMolecule(activity);
                }
                else if(oper == 5) {
                    makeGroup();
                }
                else if(oper/10 == 8) {
                    turnMolecule(activity);
                }
                else if(oper == 9) {
                    Optimize();
                }
                else if(oper == 10) {
                    addDeleteHydrogen();
                }
                else if(oper == 11) {
                    duplicateGroup();
                }
                else if(oper == 13) {
                    insertFragment();
                }
            }
            else {
                oper=0;
                identifyAtom();
            }
        }

        public void mousePressed(MouseEvent e) {
            /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
             */
            sx=e.getX();
            sy=e.getY();
            lmove=false;
        }

        public void mouseReleased(MouseEvent e) {
            /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
             */
            int j;
            double aaa,dd1,xx1,yy1,zz1,w;

            fx=e.getX();
            fy=e.getY();

            if(sx != fx && sy != fy) {
                if(home.edit) {
                    lmove=false;
                    if(oper == 5) {
                        frameMolecule();
                    }
                    else if(oper == 6) {
                        if(activity.equals("Ty-")) {
                            scale*=Math.max(0.1,1+(fx-sx)/(double)tilex);
                            repaint();
                        }
                        else if(activity.equals("Ty+")) {
                            shx+=(fx-sx);
                            shy+=(fy-sy);
                            repaint();
                        }
                        else
                            moveMolecule();
                    }
                    else if(oper == 7) {
                        rotateMolecule();
                    }
                }
                else {
                    if(oper == 6)
                        shiftAtoms();
                    else if(oper == 7)
                        turnAtoms();
                }
            }
        }
        public void mouseEntered(MouseEvent e){}
        public void mouseExited(MouseEvent e){}
    }

    class MouseMoves implements MouseMotionListener {
        public void mouseDragged(MouseEvent e) {
            /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
             */
            fx=e.getX();
            fy=e.getY();

            if(home.edit) {
                if(oper == 5) {
                    lmove=true;
                    repaint();
                }
                else if(oper == 6) {
                    if(activity.equals("Ty-")) {
                        scale*=Math.max(0.1,1+(fx-sx)/(double)tilex);
                        repaint();
                    }
                    else if(activity.equals("Ty+")) {
                        shx+=(fx-sx);
                        shy+=(fy-sy);
                        repaint();
                    }
                    else
                        moveMolecule();
                    sx=fx;
                    sy=fy;
                }
                else if(oper == 7) {
                    rotateMolecule();
                    sx=fx;
                    sy=fy;
                }
                else if(oper == 12) {
                    boxResize();
                }
            }
        }
        public void mouseMoved(MouseEvent e){}
    }
}
