import java.io.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.print.*;
import javax.swing.*;

public class GraphDraw extends Basic implements ActionListener {
    /*
*********************************************************************

dl_poly/java GUI class to draw graph plots

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
     */
    public static GUIPlot plot;
    public static GraphDraw job;
    public static JTextField xlabel,ylabel,plabel;
    private static GUI home;
    private static ColorScheme art=null;
    private static String plab,xlab,ylab,fname;
    private static int keyplot,ngrid,line,npnts,sx,fx,sy,fy;
    private static JButton load,spline,close,dots,lines,prin,zoom,clear;
    private static Font fontMain=null;
    private static double[] uu,vv,xx,yy;
    private static boolean ldots,lzoom,lmove,lscale;
    private static String[] header;
    private static boolean[] lshow;
    private static final int SCREENPX=600;
    private static final int SCREENPY=525;
    private static int ptilex=SCREENPX,ptiley=SCREENPY;

    // Define the Graphical User Interface

    public GraphDraw() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
	super();
        setTitle("Graph Plotter");

        getContentPane().setBackground(art.back);
        getContentPane().setForeground(art.fore);
        setFont(fontMain);
        GridBagLayout grd = new GridBagLayout();
        GridBagConstraints gbc = new GridBagConstraints();
        getContentPane().setLayout(grd);
        gbc.fill=GridBagConstraints.BOTH;

        // Define the Plot canvas

        plot=new GUIPlot(this);
        plot.setSize(ptilex,ptiley);
        plot.setBackground(art.scrn);
        plot.setForeground(art.scrf);
        fix(plot,grd,gbc,0,0,4,8,ptilex,ptiley);

        // Load button

        load = new JButton("Load");
        load.setBackground(art.butn);
        load.setForeground(art.butf);
        fix(load,grd,gbc,4,0,1,1,20,10);

        // Spline button

        spline = new JButton("Spline");
        spline.setBackground(art.butn);
        spline.setForeground(art.butf);
        fix(spline,grd,gbc,4,1,1,1,20,10);

        // Blank button 4

        dots = new JButton("Dots");
        dots.setBackground(art.butn);
        dots.setForeground(art.butf);
        fix(dots,grd,gbc,4,2,1,1,20,10);

        // Blank button 5

        lines = new JButton("Lines");
        lines.setBackground(art.butn);
        lines.setForeground(art.butf);
        fix(lines,grd,gbc,4,3,1,1,20,10);

        // Blank button 6

        prin = new JButton("Print");
        prin.setBackground(art.butn);
        prin.setForeground(art.butf);
        fix(prin,grd,gbc,4,4,1,1,20,10);

        // Blank button 7

        zoom = new JButton("Zoom");
        zoom.setBackground(art.butn);
        zoom.setForeground(art.butf);
        fix(zoom,grd,gbc,4,5,1,1,20,10);

        // Blank button 8

        clear = new JButton("Clear");
        clear.setBackground(art.butn);
        clear.setForeground(art.butf);
        fix(clear,grd,gbc,4,6,1,1,20,10);

        // Close button

        close = new JButton("Close");
        close.setBackground(art.butn);
        close.setForeground(art.butf);
        fix(close,grd,gbc,4,7,1,1,20,10);

        // X label textfield

        xlabel=new JTextField("X axis");
        xlabel.setBackground(art.scrn);
        xlabel.setForeground(art.scrf);
        fix(xlabel,grd,gbc,0,8,1,1,5,10);

        // Y label textfield

        ylabel=new JTextField("Y axis");
        ylabel.setBackground(art.scrn);
        ylabel.setForeground(art.scrf);
        fix(ylabel,grd,gbc,1,8,1,1,5,10);

        // Plot label textfield

        plabel=new JTextField("Title");
        plabel.setBackground(art.scrn);
        plabel.setForeground(art.scrf);
        fix(plabel,grd,gbc,2,8,2,1,20,10);
        fix(new JLabel(""),grd,gbc,4,8,1,1,20,10);

        // Register action buttons

        load.addActionListener(this);
        spline.addActionListener(this);
        dots.addActionListener(this);
        prin.addActionListener(this);
        zoom.addActionListener(this);
        lines.addActionListener(this);
        clear.addActionListener(this);
        close.addActionListener(this);
        xlabel.addActionListener(this);
        ylabel.addActionListener(this);
        plabel.addActionListener(this);

    }

    public GraphDraw(GUI here) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        home=here;
        ngrid=1000;
        println("Activated Graph Drawing panel");
	art=home.art;
        line=0;
        npnts=0;
        ldots=true;
        lzoom=false;
        lmove=false;
        lscale=false;
        header= new String[4];
        job=new GraphDraw();
        job.pack();
        job.setVisible(true);
        xlabel.setText("");
        ylabel.setText("");
        plabel.setText("");

    }

    void fix(Component cmp,GridBagLayout grd,GridBagConstraints gbc,
	     int gx,int gy,int gw,int gh,int wx,int wy){
	/*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
	gbc.gridx=gx;
	gbc.gridy=gy;
	gbc.gridwidth=gw;
	gbc.gridheight=gh;
	gbc.weightx=wx;
	gbc.weighty=wy;
	gbc.ipadx=2;
	gbc.ipady=2;
	grd.setConstraints(cmp,gbc);
	getContentPane().add(cmp);
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
        if(e.getSource() instanceof JTextField) {
            xlab=xlabel.getText();
            ylab=ylabel.getText();
            plab=plabel.getText();
            plot.repaint();
        }
        else if (arg.equals("Load")) {
            if((fname=selectFileNameEnds(home,"XY"))!=null) {
                keyplot=0;
                lzoom=false;
                lscale=false;
                npnts=getXY(fname,0);
                lshow=new boolean[npnts];
                plab=header[1];
                xlab=header[2];
                ylab=header[3];
                xlabel.setText(xlab);
                ylabel.setText(ylab);
                plabel.setText(plab);
                plot.repaint();
            }
            else {
                println("File selection cancelled");
            }
        }
        else if (arg.equals("Spline")) {
            if(npnts>0 && npnts<ngrid) {
                if(keyplot==0) {
                    keyplot=1;
                    uu=new double[ngrid];
                    vv=new double[ngrid];
                    lshow=new boolean[ngrid];
                    interpolate(npnts,ngrid,xx,yy,uu,vv);
                }
                else {
                    keyplot=0;
                }
                plot.repaint();
            }
        }
        else if (arg.equals("Print")) {
            plot.printOut();
        }
        else if (arg.equals("Dots")) {
            if(ldots) {
                ldots=false;
            }
            else {
                ldots=true;
            }
            plot.repaint();
        }
        else if (arg.equals("Lines")) {
            line++;
            if(line>2)line=0;
            plot.repaint();
        }
        else if (arg.equals("Clear")) {
            npnts=0;
            plot.repaint();
        }
        else if (arg.equals("Zoom")) {
            if(lzoom) {
                lzoom=false;
                lscale=false;
                plot.repaint();
            }
            else {
                lzoom=true;
                lscale=false;
            }
        }
        else if (arg.equals("Close")) {
            job.setVisible(false);
        }
    }

    void extraPlot(int npts,double xi[],double yi[]) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        // Entry for other classes

        line=0;
        keyplot=0;
        npnts=npts;
        ldots=false;
        lzoom=false;
        lmove=false;
        lscale=false;
        xx=new double[npnts];
        yy=new double[npnts];
        lshow=new boolean[npnts];
        System.arraycopy(xi,0,xx,0,npnts);
        System.arraycopy(yi,0,yy,0,npnts);
        plab=plabel.getText();
        xlab=xlabel.getText();
        ylab=ylabel.getText();
        plot.repaint();
    }

    int getXY(String fname,int np) {
        /*
*********************************************************************

dl_poly/java routine to read a simple XY file

copyright - daresbury laboratory
author    - w.smith january 2001

*********************************************************************
         */
        int n,k,m,mxpnts;
        String record;

        n=0;
        k=0;
        m=0;
        mxpnts=250;
        xx=new double[mxpnts];
        yy=new double[mxpnts];
        header[0]="DATA SOURCE";
        header[1]="DATA TITLE";
        header[2]="X-Coordinate";
        header[3]="Y-Coordinate";

        try {
            LineNumberReader lnr = new LineNumberReader(new FileReader(fname));
            println("Reading file: "+fname);
            while((record=lnr.readLine()) != null) {
                if(record.charAt(0) == '#') {
                    header[m]=new String(record.substring(1));
                    m=(m+1)%4;
                }
                else {
                    if(record.charAt(0)== '&') {
                        k++;
                    }
                    else if(k==np) {
                        if(n==mxpnts) {
                            double xt[]=new double[2*mxpnts];
                            double yt[]=new double[2*mxpnts];
                            System.arraycopy(xx,0,xt,0,mxpnts);
                            System.arraycopy(yy,0,yt,0,mxpnts);
                            mxpnts*=2;
                            xx=xt;
                            yy=yt;
                        }
                        xx[n]=BML.giveDouble(record,1);
                        yy[n]=BML.giveDouble(record,2);
                        n++;
                    }
                }
            }
            lnr.close();
        }
        catch(FileNotFoundException e) {
            println("Error - file not found: " + fname);
            return -1;
        }
        catch(Exception e) {
            println("Error reading file: " + fname + " "+e);
            return -2;
        }
        println("Number of points loaded:"+BML.fmt(n,6));
        return n;
    }

    void interpolate(int npnts,int ngrid, double xx[],double yy[],double uu[],
    double vv[]) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        // Interpolation using spline coefficients

        int np=1;
        double di,dj,sl0,sln,cc0,ccn,rpd;

        double zz[]=new double[npnts];
        double aa[]=new double[npnts];
        double dd[]=new double[npnts];
        double gg[]=new double[npnts];

        AML.spline(npnts,xx,yy,zz,aa,dd,gg);

        sl0=(yy[1]-yy[0])/dd[0]-gg[1]*dd[0]/6.0;
        sln=(yy[npnts-1]-yy[npnts-2])/dd[npnts-2]+gg[npnts-2]*dd[npnts-2]/6.0;
        cc0=yy[0]-sl0*xx[0];
        ccn=yy[npnts-1]-sln*xx[npnts-1];
        rpd=(xx[npnts-1]-xx[0])/ngrid;

        for(int i=0;i<ngrid;i++) {
            uu[i]=rpd*i+xx[0];
            while(np<npnts && uu[i]>xx[np]){np++;}
            di=uu[i]-xx[np-1];
            dj=xx[np]-uu[i];
            vv[i]=(di*yy[np]+dj*yy[np-1]-di*dj*
            ((dd[np-1]+dj)*gg[np-1]+(dd[np-1]+di)*gg[np])/6.0)/dd[np-1];
        }
    }

    class GUIScreenSizer extends ComponentAdapter {
        /*
*********************************************************************

dl_poly/java GUI class to adjust screen size

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        public void componentResized(ComponentEvent e) {
            Dimension arg = ((Component)e.getSource()).getSize();
            ptilex = arg.width;
            ptiley = arg.height;
            lzoom=false;
            lscale=false;
            plot.repaint();
        }
    }

    class GUIPlot extends Canvas implements Printable {
        /*
*********************************************************************

dl_poly/java GUI class to plot arrays

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        GraphDraw home;
        Image offscreenImg;
        Graphics offscreenG;
        int[] xpos,ypos;
        double[] scale=new double[6];

        public GUIPlot(GraphDraw here) {
            /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
             */
            home=here;
            xpos=new int[ngrid];
            ypos=new int[ngrid];

            // register object for screen resize

            addComponentListener(new GUIScreenSizer());

            // Register Mouse Listeners

            addMouseListener(new MousePoints());
            addMouseMotionListener(new MouseMoves());

        }

        void printOut() {
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
                    println("Initiating print .......");
                    prnjob.print();
                    println("Print complete");
                }
                else {
                    println("No print initiated");
                }
            }
            catch(PrinterException pe) {
                println("Error - problem with print "+pe.toString());
            }

        }

        public int print(Graphics g,PageFormat pf,int pageIndex) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

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

        public void update(Graphics g) {
            if(!lmove) {
                g.setColor(getBackground());
                g.fillRect(0,0,ptilex,ptiley);
                g.setColor(getForeground());
            }
            paint(g);
        }

        public void paint(Graphics g) {
            /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
             */
            int j,k,u,w,rr,xd,yd,zd,ka,kb,ja,jb,dx,dy,px,py,tx,ty;

            dx=BML.nint(0.1*ptilex);
            dy=BML.nint(0.1*ptiley);
            tx=BML.nint(0.9*ptilex);
            ty=BML.nint(0.9*ptiley);

            if(npnts>0) {
                ptilex=getSize().width;
                ptiley=getSize().height;
                offscreenImg = createImage(ptilex,ptiley);
                offscreenG = offscreenImg.getGraphics();
                offscreenG.setColor(art.scrn);
                offscreenG.fillRect(0,0,ptilex,ptiley);
                offscreenG.setColor(Color.black);

                // Plot Data

                if(lmove) {
                    offscreenG.drawLine(sx,sy,sx,fy);
                    offscreenG.drawLine(sx,sy,fx,sy);
                    offscreenG.drawLine(sx,fy,fx,fy);
                    offscreenG.drawLine(fx,sy,fx,fy);
                }
                if(keyplot==0) {
                    plotScale(npnts,scale,lshow,xx,yy);
                    graphFrame(g,scale);
                    if(ldots && lshow[0])offscreenG.fillRect(xpos[0]-3,ypos[0]-3,6,6);
                    for (int i=1;i<npnts;i++) {
                        if(ldots && lshow[i])offscreenG.fillRect(xpos[i]-3,ypos[i]-3,6,6);
                        if(lshow[i-1] && lshow[i]) {
                            offscreenG.drawLine(xpos[i-1],ypos[i-1],xpos[i],ypos[i]);
                            if(line>0)offscreenG.drawLine(xpos[i-1],ypos[i-1]+1,xpos[i],ypos[i]+1);
                            if(line>1)offscreenG.drawLine(xpos[i-1],ypos[i-1]-1,xpos[i],ypos[i]-1);
                        }
                    }
                }
                else {
                    plotScale(ngrid,scale,lshow,uu,vv);
                    graphFrame(g,scale);
                    for(int i=1;i<ngrid;i++) {
                        if(lshow[i-1] && lshow[i]) {
                            offscreenG.drawLine(xpos[i-1],ypos[i-1],xpos[i],ypos[i]);
                            if(line>0)offscreenG.drawLine(xpos[i-1],ypos[i-1]+1,xpos[i],ypos[i]+1);
                            if(line>1)offscreenG.drawLine(xpos[i-1],ypos[i-1]-1,xpos[i],ypos[i]-1);
                        }
                    }
                    if(ldots) {
                        for (int i=0;i<npnts;i++) {
                            px=BML.nint(scale[4]*(xx[i]-scale[0])+dx);
                            py=ptiley-BML.nint(scale[5]*(yy[i]-scale[2])+dy);
                            if(px>=dx && px<=tx && py>=dy && py<=ty) offscreenG.fillRect(px-3,py-3,6,6);
                        }
                    }
                }
                g.drawImage(offscreenImg,0,0,this);
            }
        }

        void plotScale(int nnn,double scale[],boolean lshow[],double xx[],double yy[]) {
            /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
             */
            int dx,dy,px,py,tx,ty;
            double x0,x1,y0,y1;

            px=BML.nint(0.8*ptilex);
            py=BML.nint(0.8*ptiley);
            dx=BML.nint(0.1*ptilex);
            dy=BML.nint(0.1*ptiley);
            tx=BML.nint(0.9*ptilex);
            ty=BML.nint(0.9*ptiley);

            // find graph scaling factors

            x0=xx[0];
            x1=xx[nnn-1];
            y0=yy[0];
            y1=yy[0];
            for(int i=1;i<nnn;i++) {
                y0=Math.min(yy[i],y0);
                y1=Math.max(yy[i],y1);
            }
            if(y0 == y1) {
                if(y0==0.0) {
                    y0=-1.0;
                    y1= 1.0;
                }
                else {
                    y0=y0-0.1*(Math.abs(y0)+Math.abs(y1));
                    y1=y1+0.1*(Math.abs(y0)+Math.abs(y1));
                }
            }
            scale[0]=xx[0];
            scale[1]=xx[nnn-1];
            scale[2]=y0;
            scale[3]=y1;
            if(lscale) {
                scale[0]=x0+(Math.min(sx,fx)-dx)*(x1-x0)/px;
                scale[1]=x0+(Math.max(sx,fx)-dx)*(x1-x0)/px;
                scale[2]=y0+(ptiley-Math.max(sy,fy)-dy)*(y1-y0)/py;
                scale[3]=y0+(ptiley-Math.min(sy,fy)-dy)*(y1-y0)/py;
            }
            scale[5]=py/(scale[3]-scale[2]);
            if(BML.nint(Math.abs(scale[2]*scale[5]))>0)
                scale[2]-=(10.0/scale[5]);
            scale[3]+=(10.0/scale[5]);
            scale[5]=py/(scale[3]-scale[2]);
            scale[4]=px/(scale[1]-scale[0]);

            // calculate plot positions

            if(nnn > ngrid) {
                ngrid=nnn;
                xpos=new int[ngrid];
                ypos=new int[ngrid];
            }
            for(int i=0;i<nnn;i++) {
                lshow[i]=true;
                xpos[i]=BML.nint(scale[4]*(xx[i]-scale[0])+dx);
                ypos[i]=ptiley-BML.nint(scale[5]*(yy[i]-scale[2])+dy);
                if(xpos[i]<dx || xpos[i]>tx || ypos[i]<dy || ypos[i]>ty)lshow[i]=false;
            }
        }

        void graphFrame(Graphics g,double scale[]) {
            /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
             */
            Font f=new Font("Helvetica",Font.PLAIN,12);
            g.setFont(f);
            int ax,ay,xp,yp,dx,dy,px,py,nticks=5;
            double[] xs,ys;
            String[] xl,yl;

            px=BML.nint(0.8*ptilex);
            dx=BML.nint(0.1*ptilex);
            py=BML.nint(0.8*ptiley);
            dy=BML.nint(0.1*ptiley);
            ax=dx;
            ay=dy+py;

            // draw box

            offscreenG.drawLine(dx,dy,(dx+px),dy);
            if(line>0)offscreenG.drawLine(dx,(dy-1),(dx+px),(dy-1));
            if(line>1)offscreenG.drawLine(dx,(dy+1),(dx+px),(dy+1));
            offscreenG.drawLine(dx,(dy+py),(dx+px),(dy+py));
            if(line>0)offscreenG.drawLine(dx,(dy+py-1),(dx+px),(dy+py-1));
            if(line>1)offscreenG.drawLine(dx,(dy+py+1),(dx+px),(dy+py+1));
            offscreenG.drawLine(dx,dy,dx,(dy+py));
            if(line>0)offscreenG.drawLine((dx-1),dy,(dx-1),(dy+py));
            if(line>1)offscreenG.drawLine((dx+1),dy,(dx+1),(dy+py));
            offscreenG.drawLine((dx+px),dy,(dx+px),(dy+py));
            if(line>0)offscreenG.drawLine((dx+px-1),dy,(dx+px-1),(dy+py));
            if(line>1)offscreenG.drawLine((dx+px+1),dy,(dx+px+1),(dy+py));

            // x axis

            if(scale[2]<0.0 && scale[3]>0.0) {
                ay=BML.nint(dy+scale[5]*scale[3]);
                offscreenG.drawLine(dx,ay,(dx+px),ay);
                if(line>0)offscreenG.drawLine(dx,(ay-1),(dx+px),(ay-1));
                if(line>1)offscreenG.drawLine(dx,(ay+1),(dx+px),(ay+1));
            }
            ay=dy+py;
            xs=new double[nticks];
            xl=new String[nticks];
            axisValues(nticks,scale[0],scale[1],xs,xl);
            for(int i=0;i<nticks;i++) {
                if(xs[i]>=scale[0] && xs[i]<=scale[1]) {
                    xp=BML.nint(scale[4]*(xs[i]-scale[0])+dx);
                    offscreenG.drawLine(xp,ay,xp,ay-7);
                    if(line>0)offscreenG.drawLine(xp-1,ay,xp-1,ay-7);
                    if(line>1)offscreenG.drawLine(xp+1,ay,xp+1,ay-7);
                    offscreenG.drawLine(xp,dy,xp,dy+7);
                    if(line>0)offscreenG.drawLine(xp-1,dy,xp-1,dy+7);
                    if(line>1)offscreenG.drawLine(xp+1,dy,xp+1,dy+7);
                    offscreenG.drawString(xl[i],xp-24,ay+18);
                }
            }

            // y axis

            if(scale[0]<0.0 && scale[1]>0.0) {
                ax=BML.nint(dx-scale[4]*scale[0]);
                offscreenG.drawLine(ax,dy,ax,(dy+py));
                if(line>0)offscreenG.drawLine((ax-1),dy,(ax-1),(dy+py));
                if(line>0)offscreenG.drawLine((ax+1),dy,(ax+1),(dy+py));
            }
            ys=new double[nticks];
            yl=new String[nticks];
            axisValues(nticks,scale[2],scale[3],ys,yl);
            for(int i=0;i<nticks;i++) {
                if(ys[i]>=scale[2] && ys[i]<=scale[3]) {
                    yp=ptiley-BML.nint(scale[5]*(ys[i]-scale[2])+dy);
                    offscreenG.drawLine(ax,yp,ax+7,yp);
                    if(line>0)offscreenG.drawLine(ax,yp-1,ax+7,yp-1);
                    if(line>1)offscreenG.drawLine(ax,yp+1,ax+7,yp+1);
                    offscreenG.drawLine(ax+px,yp,ax+px-7,yp);
                    if(line>0)offscreenG.drawLine(ax+px,yp-1,ax+px-7,yp-1);
                    if(line>1)offscreenG.drawLine(ax+px,yp+1,ax+px-7,yp+1);
                    offscreenG.drawString(yl[i],(int)(dx*1.25)-10*yl[i].length(),yp+5);
                }
            }

            // Plot labels

            Font h=new Font("Helvetica",Font.BOLD,14);
            g.setFont(h);
            offscreenG.drawString(plab,ptilex/2-5*plab.length(),dy/2);
            xp=BML.nint(scale[4]*(0.5*(xs[nticks-2]+xs[nticks-3])-scale[0])+dx);
            offscreenG.drawString(xlab,xp,dy+py+dy/2);
            yp=ptiley-BML.nint(scale[5]*(0.5*(ys[nticks-1]+ys[nticks-2])-scale[2])+dy);
            offscreenG.drawString(ylab,dx/4,yp);

        }

        void axisValues(int nticks,double xmin,double xmax,double xs[],String xl[]) {
            /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
             */
            int u,d;
            double t0,t1,xstp;

            u=decade(Math.max(Math.abs(xmin),Math.abs(xmax)));
            xstp=(xmax-xmin)/(nticks-1);
            d=decade(xstp);
            t0=(double)BML.nint(xmin*Math.pow(10.0,-d))*Math.pow(10.0,d);
            t1=(double)BML.nint(xstp*Math.pow(10.0,-d))*Math.pow(10.0,d);
            for(int i=0;i<nticks;i++) {
                xs[i]=t0+t1*i;
                xl[i]=BML.fmt(xs[i]*Math.pow(10.0,-u),u-d+4)+"E"+String.valueOf(u);
            }
        }
        int decade(double r) {
            int d=0;
            while(r<1.0){r*=10.0;d--;}
            while(r>=10.0){r/=10.0;d++;}
            return d;
        }

        //support classes

        class MousePoints extends MouseAdapter {
            public void mousePressed(MouseEvent e) {
                /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
                 */
                if(lzoom) {
                    sx=e.getX();
                    sy=e.getY();
                }
            }
            public void mouseReleased(MouseEvent e) {
                /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
                 */
                if(lzoom && !lscale) {
                    fx=e.getX();
                    fy=e.getY();
                    if(sx != fx && sy != fy) {
                        lmove=false;
                        lscale=true;
                        repaint();
                    }
                }
            }
        }

        class MouseMoves extends MouseMotionAdapter {
            public void mouseDragged(MouseEvent e) {
                /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
                 */
                if(lzoom && !lscale) {
                    lmove=true;
                    fx=e.getX();
                    fy=e.getY();
                    repaint();
                }
            }
            public void mouseMoved(MouseEvent e){}
        }
    }
}
