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
    public static JTextField xlabel,ylabel,plabel,info;
    private static GUI home;
    private static ColorScheme art=null;
    private static String plab,xlab,ylab,fname;
    private static int keyplot,line,npnts,ngrid,sx,fx,sy,fy;
    private static JButton load,clear,prin,close,dots,lines,zoom,spline;
    private static JButton ddx,sdx,dft,auto,lsq,limit;
    private static Font fontMain=null;
    private static double[] uu,vv,xx,yy,sss;
    private static double[] scale=new double[6];
    private static int[] xpos,ypos,ssx,ssy;
    private static boolean ldots,lzoom,lmove,llimit,llsq;
    private static String[] header;
    private static double xmin,xmax,ymin,ymax,ymid;
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
        setDefaultCloseOperation(DISPOSE_ON_CLOSE);
        setFont(fontMain);
        GridBagLayout grd = new GridBagLayout();
        GridBagConstraints gbc = new GridBagConstraints();
        getContentPane().setLayout(grd);
        gbc.fill=GridBagConstraints.BOTH;
        int n=0;

        // Define the Plot canvas

        plot=new GUIPlot(this);
        plot.setSize(ptilex,ptiley);
        plot.setBackground(art.scrn);
        plot.setForeground(art.scrf);
        fix(plot,grd,gbc,0,0,4,7,ptilex,ptiley);

        // Load button

        load = new JButton("Load");
        load.setBackground(art.butn);
        load.setForeground(art.butf);
        fix(load,grd,gbc,4,n,1,1,20,10);

        // Clear button

        clear = new JButton("Clear");
        clear.setBackground(art.butn);
        clear.setForeground(art.butf);
        fix(clear,grd,gbc,5,n++,1,1,20,10);

        // Print button

        prin = new JButton("Print");
        prin.setBackground(art.butn);
        prin.setForeground(art.butf);
        fix(prin,grd,gbc,4,n,1,1,20,10);

        // Apply limits to data set

        limit = new JButton("Limit");
        limit.setBackground(art.butn);
        limit.setForeground(art.butf);
        fix(limit,grd,gbc,5,n++,1,1,20,10);

        // Dots button

        dots = new JButton("Dots");
        dots.setBackground(art.butn);
        dots.setForeground(art.butf);
        fix(dots,grd,gbc,4,n,1,1,20,10);

        // Lines button

        lines = new JButton("Lines");
        lines.setBackground(art.butn);
        lines.setForeground(art.butf);
        fix(lines,grd,gbc,5,n++,1,1,20,10);

        // Zoom button

        zoom = new JButton("Zoom");
        zoom.setBackground(art.butn);
        zoom.setForeground(art.butf);
        fix(zoom,grd,gbc,4,n,1,1,20,10);

        // Spline button

        spline = new JButton("Spline");
        spline.setBackground(art.butn);
        spline.setForeground(art.butf);
        fix(spline,grd,gbc,5,n++,1,1,20,10);

        // Integrate button

        sdx = new JButton("S_dx");
        sdx.setBackground(art.butn);
        sdx.setForeground(art.butf);
        fix(sdx,grd,gbc,4,n,1,1,20,10);

        // Radial integration button

        ddx = new JButton("d/dx");
        ddx.setBackground(art.butn);
        ddx.setForeground(art.butf);
        fix(ddx,grd,gbc,5,n++,1,1,20,10);

        // Direct Fourier Transform button

        dft = new JButton("DFT");
        dft.setBackground(art.butn);
        dft.setForeground(art.butf);
        fix(dft,grd,gbc,4,n,1,1,20,10);

        // Radial Fourier Transform button

        auto = new JButton("Auto");
        auto.setBackground(art.butn);
        auto.setForeground(art.butf);
        fix(auto,grd,gbc,5,n++,1,1,20,10);

        // Linear least squares button

        lsq = new JButton("LSQ");
        lsq.setBackground(art.butn);
        lsq.setForeground(art.butf);
        fix(lsq,grd,gbc,4,n,1,1,20,10);

        // Close button

        close = new JButton("Close");
        close.setBackground(art.butn);
        close.setForeground(art.butf);
        fix(close,grd,gbc,5,n++,1,1,20,10);

        // X label textfield

        xlabel=new JTextField("X axis");
        xlabel.setBackground(art.scrn);
        xlabel.setForeground(art.scrf);
        fix(xlabel,grd,gbc,0,n,1,1,5,10);

        // Y label textfield

        ylabel=new JTextField("Y axis");
        ylabel.setBackground(art.scrn);
        ylabel.setForeground(art.scrf);
        fix(ylabel,grd,gbc,1,n,1,1,5,10);

        // Plot label textfield

        plabel=new JTextField("Plot title");
        plabel.setBackground(art.scrn);
        plabel.setForeground(art.scrf);
        fix(plabel,grd,gbc,2,n,2,1,20,10);

        // Information textfield

        info=new JTextField(" ");
        info.setBackground(art.scrn);
        info.setForeground(art.scrf);
        fix(info,grd,gbc,4,n,2,1,5,10);

        // Register action buttons

        load.addActionListener(this);
        clear.addActionListener(this);
        prin.addActionListener(this);
        close.addActionListener(this);
        dots.addActionListener(this);
        lines.addActionListener(this);
        zoom.addActionListener(this);
        spline.addActionListener(this);
        ddx.addActionListener(this);
        sdx.addActionListener(this);
        dft.addActionListener(this);
        auto.addActionListener(this);
        lsq.addActionListener(this);
        limit.addActionListener(this);
        xlabel.addActionListener(this);
        ylabel.addActionListener(this);
        plabel.addActionListener(this);
        info.addActionListener(this);

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
        println("Activated Graph Drawing panel");
        art=home.art;
        line=0;
        npnts=0;
        ldots=true;
        lzoom=false;
        llimit=false;
        llsq=false;
        lmove=false;
        header= new String[4];
        job=new GraphDraw();
        job.pack();
        job.setVisible(true);

    }

    public void actionPerformed(ActionEvent e) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        double sum;

        String arg = (String)e.getActionCommand();
        if(e.getSource() instanceof JTextField) {
            xlab=xlabel.getText();
            ylab=ylabel.getText();
            plab=plabel.getText();
            plot.repaint();
        }
        else if (arg.equals("Load")) {
            lzoom=false;
            llimit=false;
            if((fname=selectFileNameEnds(home,"XY"))!=null) {
                keyplot=0;
                npnts=getXY(fname,0);
                xmin=xx[0];
                xmax=xx[0];
                ymin=yy[0];
                ymax=yy[0];
                for(int i=1;i<npnts;i++) {
                    xmin=Math.min(xmin,xx[i]);
                    xmax=Math.max(xmax,xx[i]);
                    ymin=Math.min(ymin,yy[i]);
                    ymax=Math.max(ymax,yy[i]);
                }
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
            llimit=false;
            if(npnts>0) {
                if(keyplot==0) {
                    keyplot=1;
                    ngrid=interpolate();
                }
                else
                    keyplot=0;
                plot.repaint();
            }
        }
        else if (arg.equals("S_dx")) {
            llimit=false;
            if(npnts>0) {
                sum=integrate();
                if(keyplot > 0)
                    ngrid=interpolate();
                info.setText("S_dx = "+BML.fmt(sum,8));
                plab="Integral of "+ylab;
                plabel.setText(plab);
                plot.repaint();
            }
        }
        else if (arg.equals("d/dx")) {
            llimit=false;
            if(npnts>0) {
                differentiate();
                if(keyplot > 0)
                    ngrid=interpolate();
                info.setText(" ");
                plab="d/dx of "+ylab;
                plabel.setText(plab);
                plot.repaint();
            }
        }
        else if (arg.equals("LSQ")) {
            llimit=false;
            if(npnts>0) {
                llsq=true;
                keyplot=0;
                sss=leastsquares();
                info.setText("C="+BML.fmt(sss[0],9)+" "+"M="+BML.fmt(sss[1],9));
                println("GraphDraw: Least squares fit for: "+plab);
                println("Intercept = "+BML.fmt(sss[0],10));
                println("Gradient  = "+BML.fmt(sss[1],10));
                plab="LSQ fit of "+ylab;
                plabel.setText(plab);
                plot.repaint();
            }
        }
        else if (arg.equals("Print")) {
            llimit=false;
            plot.printOut();
        }
        else if (arg.equals("Dots")) {
            llimit=false;
            if(ldots) {
                ldots=false;
            }
            else {
                ldots=true;
            }
            plot.repaint();
        }
        else if (arg.equals("Lines")) {
            llimit=false;
            line++;
            if(line>2)line=0;
            plot.repaint();
        }
        else if (arg.equals("Clear")) {
            npnts=0;
            lzoom=false;
            llimit=false;
            xlabel.setText("X axis");
            ylabel.setText("Y axis");
            plabel.setText("Plot title");
            plot.repaint();
        }
        else if (arg.equals("DFT")) {
            lzoom=false;
            llimit=false;
            dftPlot();
            plot.repaint();
        }
        else if (arg.equals("Auto")) {
            lzoom=false;
            llimit=false;
            autoPlot();
            plot.repaint();
        }
        else if (arg.equals("Limit")) {
            lzoom=false;
            if(llimit) {
                llimit=false;
                xmin=xx[0];
                xmax=xx[0];
                ymin=yy[0];
                ymax=yy[0];
                for(int i=0;i<npnts;i++) {
                    xmin=Math.min(xmin,xx[i]);
                    xmax=Math.max(xmax,xx[i]);
                    ymin=Math.min(ymin,yy[i]);
                    ymax=Math.max(ymax,yy[i]);
                }
                plot.repaint();
            }
            else {
                llimit=true;
            }
        }
        else if (arg.equals("Zoom")) {
            llimit=false;
            if(lzoom) {
                lzoom=false;
                xmin=xx[0];
                xmax=xx[0];
                ymin=yy[0];
                ymax=yy[0];
                for(int i=1;i<npnts;i++) {
                    xmin=Math.min(xmin,xx[i]);
                    xmax=Math.max(xmax,xx[i]);
                    ymin=Math.min(ymin,yy[i]);
                    ymax=Math.max(ymax,yy[i]);
                }
                plot.repaint();
            }
            else {
                lzoom=true;
            }
        }
        else if (arg.equals("Close")) {
            job.dispose();
        }
    }

    void dftPlot() {
        /*
*********************************************************************

dl_poly/java GUI routine for Discrete Fourier Transform Plot

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */

        double a0,a1,a2,avg,ccc,del;
        boolean safe=true;
        int nfft=1;

        while(nfft < npnts)
            nfft*=2;
        int key[]=new int[nfft];
        double fta[][]=new double[2][nfft];
        double ftb[][]=new double[2][nfft];
        double work[][]=new double[2][nfft];

        // FFT initialisation

        fft(safe,1,-1,nfft,key,fta,work,ftb);

        // Average of initial function

        avg=0.0;
        for(int i=0;i<npnts;i++)
            avg+=(yy[i]/npnts);

        // load function for FFT (use Blackman-Harris window)

        a0=0.42;
        a1=0.5;
        a2=0.08;
        del=2.0*Math.PI/npnts;
        for(int i=0;i<nfft;i++){
            fta[1][i]=0.0;
            if(i<npnts) {
                ccc=Math.cos(del*i);
                fta[0][i]=(yy[i]-avg)*(a0-a1*ccc+a2*(2*ccc*ccc-1));
            }
            else
                fta[0][i]=0.0;
        }

        // calculate Fourier transform

        fft(safe,0,-1,nfft,key,fta,work,ftb);

        // calculate real amplitudes

        del=2.0*Math.PI*((double)(npnts-1)/(double)nfft)/(xx[npnts-1]-xx[0]);
        npnts=nfft/2;
        xx=new double[npnts];
        yy=new double[npnts];
        xmin=0.0;
        xmax=del*(double)(npnts-1);
        ymin=Math.sqrt(Math.pow(ftb[0][0],2)+Math.pow(ftb[1][0],2));
        ymax=ymin;
        for(int i=0;i<npnts;i++) {
            xx[i]=del*(double)i;
            yy[i]=Math.sqrt(Math.pow(ftb[0][i],2)+Math.pow(ftb[1][i],2));
            ymin=Math.min(ymin,yy[i]);
            ymax=Math.max(ymax,yy[i]);
        }

        // Plotting information

        xlab="Omega";
        ylab="Amplitude";
        plab="Fourier Transform of "+ylab;
        xlabel.setText(xlab);
        ylabel.setText(ylab);
        plabel.setText(plab);

    }


    void autoPlot() {
        /*
*********************************************************************

dl_poly/java GUI routine for Autocorrelation Plot

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */

        double avg,del;
        boolean safe=true;
        int nfft=1;

        while(nfft < npnts)
            nfft*=2;

        nfft*=2;
        int key[]=new int[nfft];
        double fta[][]=new double[2][nfft];
        double ftb[][]=new double[2][nfft];
        double work[][]=new double[2][nfft];

        // FFT initialisation

        fft(safe,1,-1,nfft,key,fta,work,ftb);

        // Average of initial function

        avg=0.0;
        for(int i=0;i<npnts;i++)
            avg+=(yy[i]/npnts);

        // load function for FFT

        for(int i=0;i<nfft;i++){
            fta[1][i]=0.0;
            if(i<npnts)
                fta[0][i]=yy[i]-avg;
            else
                fta[0][i]=0.0;
        }

        // calculate Fourier transform

        fft(safe,0,-1,nfft,key,fta,work,ftb);

        // calculate Fourier transform of correlation function

        for(int i=0;i<nfft;i++) {
            fta[0][i]=(Math.pow(ftb[0][i],2)+Math.pow(ftb[1][i],2))/((double)nfft);
            fta[1][i]=0.0;
        }

        // calculate inverse Fourier transform

        fft(safe,0,1,nfft,key,fta,work,ftb);

        // calculate correlation function

        del=(xx[npnts-1]-xx[0])/(double)(npnts-1);
        xmin=0.0;
        xmax=del*(double)(npnts-1);
        ymin=(Math.pow(ftb[0][0],2)+Math.pow(ftb[1][0],2))*del/((double)(npnts));
        ymax=ymin;
        for(int i=0;i<npnts;i++) {
            xx[i]=del*(double)i;
            yy[i]=(Math.pow(ftb[0][i],2)+Math.pow(ftb[1][i],2))*del/((double)(npnts-i));
            ymin=Math.min(ymin,yy[i]);
            ymax=Math.max(ymax,yy[i]);
        }

        // Plotting information

        plab=" Autocorrelation of "+ylab;
        ylab=" C(t)";
        xlabel.setText(xlab);
        ylabel.setText(ylab);
        plabel.setText(plab);

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
        llimit=false;
        lmove=false;
        xx=new double[npnts];
        yy=new double[npnts];
        xmin=xi[0];
        xmax=xi[0];
        ymin=yi[0];
        ymax=yi[0];
        for(int i=0;i<npnts;i++) {
            xx[i]=xi[i];
            yy[i]=yi[i];
            xmin=Math.min(xmin,xx[i]);
            xmax=Math.max(xmax,xx[i]);
            ymin=Math.min(ymin,yy[i]);
            ymax=Math.max(ymax,yy[i]);
        }
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
        mxpnts=100;
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
                            double xt[]=new double[mxpnts+100];
                            double yt[]=new double[mxpnts+100];
                            System.arraycopy(xx,0,xt,0,mxpnts);
                            System.arraycopy(yy,0,yt,0,mxpnts);
                            mxpnts+=100;
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

    int interpolate() {
        /*
*********************************************************************

dl_poly/java GUI routine: Interpolation using spline coefficients

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */

        int np=1;
        double di,dj,rpd;

        double zz[]=new double[npnts];
        double aa[]=new double[npnts];
        double dd[]=new double[npnts];
        double gg[]=new double[npnts];

        ngrid=Math.max(npnts,5000);
        uu=new double[ngrid];
        vv=new double[ngrid];

        AML.spline(npnts,xx,yy,zz,aa,dd,gg);

        rpd=(xx[npnts-1]-xx[0])/ngrid;
        for(int i=0;i<ngrid;i++) {
            uu[i]=rpd*i+xx[0];
            while(np<npnts && uu[i]>xx[np]){np++;}
            di=uu[i]-xx[np-1];
            dj=xx[np]-uu[i];
            vv[i]=(di*yy[np]+dj*yy[np-1]-di*dj*
            ((dd[np-1]+dj)*gg[np-1]+(dd[np-1]+di)*gg[np])/6.0)/dd[np-1];
        }
        return ngrid;
    }

    double integrate() {
        /*
*********************************************************************

dl_poly/java GUI routine: Integration using spline coefficients

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */

        double aaa,sum;

        double zz[]=new double[npnts];
        double aa[]=new double[npnts];
        double dd[]=new double[npnts];
        double gg[]=new double[npnts];

        AML.spline(npnts,xx,yy,zz,aa,dd,gg);

        sum=0.0;
        aaa=yy[0];
        yy[0]=0.0;
        ymin=0.0;
        ymax=0.0;

        for(int i=1;i<npnts;i++) {
            sum+=0.5*Math.pow(xx[i]-xx[i-1],2)/dd[i-1]*(aaa+yy[i]-(xx[i]-xx[i-1])*(gg[i]+gg[i-1])/18.0*(dd[i-1]+0.5*(xx[i]-xx[i-1])));
            ymin=Math.min(ymin,sum);
            ymax=Math.max(ymax,sum);
            aaa=yy[i];
            yy[i]=sum;
        }
        return sum;
    }

    void differentiate() {
        /*
*********************************************************************

dl_poly/java GUI routine: Differentiation using spline coefficients

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */

        double aaa,diff;

        double zz[]=new double[npnts];
        double aa[]=new double[npnts];
        double dd[]=new double[npnts];
        double gg[]=new double[npnts];

        AML.spline(npnts,xx,yy,zz,aa,dd,gg);

        aaa=yy[0];
        diff=(yy[1]-aaa-(xx[1]-xx[0])/6.0*((gg[1]+gg[0])*dd[0]+(xx[1]-xx[0])*gg[0]))/dd[0];
        yy[0]=diff;
        ymin=diff;
        ymax=diff;

        for(int i=1;i<npnts;i++) {
            diff=(yy[i]-aaa+(xx[i]-xx[i-1])/6.0*((gg[i]+gg[i-1])*dd[i-1]+(xx[i]-xx[i-1])*gg[i]))/dd[i-1];
            ymin=Math.min(ymin,diff);
            ymax=Math.max(ymax,diff);
            aaa=yy[i];
            yy[i]=diff;
        }
    }

    double[] leastsquares() {
        /*
*********************************************************************

dl_poly/java GUI routine: Linear least squares fit

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
*/
        double det;
        double[] uuu= new double[2];
        double[] sss= new double[2];
        double[][] mat= new double[2][2];

        uuu[0]=0.0;
        uuu[1]=0.0;
        mat[0][1]=0.0;
        mat[1][0]=0.0;
        mat[1][1]=0.0;
        mat[0][0]=(double)npnts;
        for(int i=0;i<npnts;i++){
            uuu[0]=uuu[0]+yy[i];
            uuu[1]=uuu[1]+yy[i]*xx[i];
            mat[0][1]=mat[0][1]+xx[i];
            mat[1][0]=mat[1][0]+xx[i];
            mat[1][1]=mat[1][1]+xx[i]*xx[i];
        }
        det=mat[0][0]*mat[1][1]-mat[0][1]*mat[1][0];
        sss[0]=(mat[1][1]*uuu[0]-mat[0][1]*uuu[1])/det;
        sss[1]=(mat[0][0]*uuu[1]-mat[1][0]*uuu[0])/det;

        return sss;
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
            llimit=false;
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

        public GUIPlot(GraphDraw here) {
            /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
             */
            home=here;

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
            boolean lgo;
            int dx,dy,px,py;

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
                    plotScale(npnts,scale,xx,yy);
                    graphFrame(g,scale);
                    lgo=false;
                    for (int i=0;i<npnts;i++) {
                        if(xx[i] >= scale[0] && xx[i] <= scale[1] && yy[i] >= scale[2] && yy[i] <= scale[3]) {
                            if(ldots)
                                offscreenG.fillRect(xpos[i]-3,ypos[i]-3,6,6);
                            if(lgo && xx[i-1] >= scale[0] && xx[i-1] <= scale[1] && yy[i-1] >= scale[2] && yy[i-1] <= scale[3]) {
                                offscreenG.drawLine(xpos[i-1],ypos[i-1],xpos[i],ypos[i]);
                                if(line>0)offscreenG.drawLine(xpos[i-1],ypos[i-1]+1,xpos[i],ypos[i]+1);
                                if(line>1)offscreenG.drawLine(xpos[i-1],ypos[i-1]-1,xpos[i],ypos[i]-1);
                            }
                            lgo=true;
                        }
                    }
                    if(llsq) {
                        offscreenG.drawLine(ssx[0],ssy[0],ssx[1],ssy[1]);
                        offscreenG.drawLine(ssx[0],ssy[0]+1,ssx[1],ssy[1]+1);
                        offscreenG.drawLine(ssx[0],ssy[0]-1,ssx[1],ssy[1]-1);
                        llsq=false;
                    }
                }
                else {
                    plotScale(ngrid,scale,uu,vv);
                    graphFrame(g,scale);
                    lgo=false;
                    for(int i=1;i<ngrid;i++) {
                        if(uu[i] >= scale[0] && uu[i] <= scale[1] && vv[i] >= scale[2] && vv[i] <= scale[3]) {
                            if(lgo && uu[i-1] >= scale[0] && uu[i-1] <= scale[1] && vv[i-1] >= scale[2] && vv[i-1] <= scale[3]) {
                                offscreenG.drawLine(xpos[i-1],ypos[i-1],xpos[i],ypos[i]);
                                if(line>0)offscreenG.drawLine(xpos[i-1],ypos[i-1]+1,xpos[i],ypos[i]+1);
                                if(line>1)offscreenG.drawLine(xpos[i-1],ypos[i-1]-1,xpos[i],ypos[i]-1);
                            }
                            lgo=true;
                        }
                    }
                    if(ldots) {
                        dx=BML.nint(0.1*ptilex);
                        dy=BML.nint(0.1*ptiley);
                        for (int i=0;i<npnts;i++) {
                            if(xx[i] >= scale[0] && xx[i] <= scale[1] && yy[i] >= scale[2] && yy[i] <= scale[3]) {
                                px=BML.nint(scale[4]*(xx[i]-scale[0])+dx);
                                py=ptiley-BML.nint(scale[5]*(yy[i]-scale[2])+dy);
                                offscreenG.fillRect(px-3,py-3,6,6);
                            }
                        }
                    }
                }
                g.drawImage(offscreenImg,0,0,this);
            }
        }

        void plotScale(int nnn,double scale[],double xx[],double yy[]) {
            /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
             */
            int dx,dy,px,py;
            double ds;

            px=BML.nint(0.8*ptilex);
            py=BML.nint(0.8*ptiley);
            dx=BML.nint(0.1*ptilex);
            dy=BML.nint(0.1*ptiley);

            // find graph scaling factors

            scale[0]=xmin;
            scale[1]=xmax;
            if(ymin == ymax) {
                if(ymin==0.0) {
                    scale[2]=-10.0;
                    scale[3]= 10.0;
                }
                else {
                    scale[2]=0;
                    scale[3]=2*ymin;
                }
            }
            else {
                if(ymin < 0 || ymin > (ymax-ymin)/20.0)
                    scale[2]=ymin-(ymax-ymin)/20.0;
                else
                    scale[2]=ymin;
                if(ymax > 0 || ymax < -(ymax-ymin)/20.0)
                    scale[3]=ymax+(ymax-ymin)/20.0;
                else
                    scale[3]=ymax;
            }
            scale[4]=px/(scale[1]-scale[0]);
            scale[5]=py/(scale[3]-scale[2]);

            // calculate plot positions

            xpos=new int[nnn];
            ypos=new int[nnn];
            for(int i=0;i<nnn;i++) {
                xpos[i]=BML.nint(scale[4]*(xx[i]-scale[0])+dx);
                ypos[i]=ptiley-BML.nint(scale[5]*(yy[i]-scale[2])+dy);
            }
            if(llsq) {
                ssx=new int[2];
                ssy=new int[2];
                ssx[0]=BML.nint(scale[4]*(xx[0]-scale[0])+dx);
                ssx[1]=BML.nint(scale[4]*(xx[nnn-1]-scale[0])+dx);
                ssy[0]=ptiley-BML.nint(scale[5]*(sss[0]+sss[1]*xx[0]-scale[2])+dy);
                ssy[1]=ptiley-BML.nint(scale[5]*(sss[0]+sss[1]*xx[nnn-1]-scale[2])+dy);
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
            int ax,ay,xp,yp,dx,dy,px,py,ntkx,ntky;
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
            xs=new double[8];
            xl=new String[8];
            ntkx=axisValues(scale[0],scale[1],xs,xl);
            for(int i=0;i<ntkx;i++) {
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
            ys=new double[8];
            yl=new String[8];
            ntky=axisValues(scale[2],scale[3],ys,yl);
            for(int i=0;i<ntky;i++) {
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
            xp=BML.nint(scale[4]*(0.5*(xs[ntkx-2]+xs[ntkx-3])-scale[0])+dx);
            offscreenG.drawString(xlab,xp,dy+py+dy/2);
            yp=ptiley-BML.nint(scale[5]*(0.5*(ys[ntky-1]+ys[ntky-2])-scale[2])+dy);
            offscreenG.drawString(ylab,dx/4,yp);

        }

        int axisValues(double xlo,double xhi,double xs[],String xl[]) {
            /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
             */
            int i,u,d,nticks;
            double t0,t1,xstp;

            nticks=5;
            u=decade(Math.max(Math.abs(xlo),Math.abs(xhi)));
            xstp=(xhi-xlo)/(nticks-1);
            d=decade(xstp);
            t0=(double)BML.nint(xlo*Math.pow(10.0,-d))*Math.pow(10.0,d);
            t1=(double)BML.nint(xstp*Math.pow(10.0,-d))*Math.pow(10.0,d);
            i=0;
            while((t0+t1*i)<=xhi){i++;}
            nticks=i;
            for(i=0;i<nticks;i++) {
                xs[i]=t0+t1*i;
                xl[i]=BML.fmt(xs[i]*Math.pow(10.0,-u),u-d+4)+"E"+String.valueOf(u);
            }
            return nticks;
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
                int dx,dy;
                double xr,yr;
                dx=BML.nint(0.1*ptilex);
                dy=BML.nint(0.1*ptiley);

                sx=e.getX();
                sy=e.getY();
                if(npnts>0) {
                    xr=(double)(sx-dx)/scale[4]+scale[0];
                    yr=(double)(ptiley-sy-dy)/scale[5]+scale[2];
                    info.setText("(X,Y)=("+BML.fmt(xr,8)+","+BML.fmt(yr,8)+")");
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
                int n,dx,dy,xp,yp,bx,tx,by,ty;

                fx=e.getX();
                fy=e.getY();
                bx=Math.min(sx,fx);
                tx=Math.max(sx,fx);
                by=Math.min(sy,fy);
                ty=Math.max(sy,fy);
                dx=BML.nint(0.1*ptilex);
                dy=BML.nint(0.1*ptiley);
                info.setText(" ");
                if(lzoom) {
                    if(sx != fx && sy != fy) {
                        lmove=false;
                        n=0;
                        for(int i=0;i<npnts;i++) {
                            xp=BML.nint(scale[4]*(xx[i]-scale[0])+dx);
                            yp=ptiley-BML.nint(scale[5]*(yy[i]-scale[2])+dy);
                            if(xp >= bx && xp <= tx && yp >= by && yp <= ty) {n++;}
                        }
                        if(n >= 5) {
                            xmin=(double)(bx-dx)/scale[4]+scale[0];
                            xmax=(double)(tx-dx)/scale[4]+scale[0];
                            ymin=(double)(ptiley-ty-dy)/scale[5]+scale[2];
                            ymax=(double)(ptiley-by-dy)/scale[5]+scale[2];
                        }
                        repaint();
                    }
                }
                else if(llimit) {
                    if(sx != fx && sy != fy) {
                        keyplot=0;
                        lmove=false;
                        n=0;
                        for(int i=0;i<npnts;i++) {
                            xp=BML.nint(scale[4]*(xx[i]-scale[0])+dx);
                            yp=ptiley-BML.nint(scale[5]*(yy[i]-scale[2])+dy);
                            if(xp >= bx && xp <= tx && yp >= by && yp <= ty) {
                                xx[n]=xx[i];
                                yy[n]=yy[i];
                                n++;
                            }
                        }
                        xmin=(double)(bx-dx)/scale[4]+scale[0];
                        xmax=(double)(tx-dx)/scale[4]+scale[0];
                        ymin=(double)(ptiley-ty-dy)/scale[5]+scale[2];
                        ymax=(double)(ptiley-by-dy)/scale[5]+scale[2];
                        llimit=false;
                        println("Number of points is now:"+BML.fmt(n,6));
                        npnts=n;
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
                if(lzoom || llimit) {
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
