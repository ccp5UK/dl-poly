import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;

// Define the Graphical User Interface

public class ShowClusters extends Basic implements ActionListener {
        /*
*********************************************************************

dl_poly/java GUI class to show a atomic clusters in a configuration

copyright - daresbury laboratory
author    - w.smith August 2011

*********************************************************************
         */
    public static ShowClusters job;
    private static GUI home=null;
    private static JButton show,load,close;
    private static JTextField cut,target;
    private static String targets;
    private static double cutoff;
    private static Color[] shade;
    public static int num_clusters,num_selected,max_size,min_size;;
    public static int[] size,range,distribution;
    public static int[] knd,key,loc,lok;

    public ShowClusters() {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        setTitle("Show Clusters");

        getContentPane().setBackground(art.back);
        getContentPane().setForeground(art.fore);
        setDefaultCloseOperation(DISPOSE_ON_CLOSE);
        setFont(fontMain);
        GridBagLayout grd = new GridBagLayout();
        GridBagConstraints gbc = new GridBagConstraints();
        getContentPane().setLayout(grd);

        gbc.fill=GridBagConstraints.BOTH;

        int n=0;

        // Panel header

        fix(new JLabel("Show Clusters",JLabel.LEFT),grd,gbc,0,n++,3,1);

        // Define the Show button

        show = new JButton("Show");
        show.setBackground(art.butn);
        show.setForeground(art.butf);
        fix(show,grd,gbc,0,n,1,1);
        fix(new JLabel("      ",JLabel.LEFT),grd,gbc,1,n,1,1);

        // Define the load button

        load = new JButton("Load");
        load.setBackground(art.butn);
        load.setForeground(art.butf);
        fix(load,grd,gbc,2,n++,1,1);

        // Define the atoms in target clusters

        fix(new JLabel("Target Atoms:",JLabel.LEFT),grd,gbc,0,n++,3,1);
        target = new JTextField(10);
        target.setBackground(art.scrn);
        target.setForeground(art.scrf);
        fix(target,grd,gbc,0,n++,3,1);

        // Set displacement criterion

        cut = new JTextField(8);
        cut.setBackground(art.scrn);
        cut.setForeground(art.scrf);
        fix(cut,grd,gbc,0,n,1,1);
        fix(new JLabel("Max. displacement (A)  ",JLabel.LEFT),grd,gbc,1,n++,2,1);

        // Define the Close button

        close = new JButton("Close");
        close.setBackground(art.butn);
        close.setForeground(art.butf);
        fix(close,grd,gbc,2,n++,1,1);

        // Register action buttons

        show.addActionListener(this);
        load.addActionListener(this);
        close.addActionListener(this);

    }

    // Constructor method

    public ShowClusters(GUI here) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        println("Activated panel for Showing Radiation Clusters");
        home=here;

        // Set up Graphical User interface

        job = new ShowClusters();
        job.pack();

        // Set default values

        cutoff=2.5;
        cut.setText(String.valueOf(cutoff));
        target.setText("");
        job.setVisible(true);
        if(config != null && config.natms > 0)
            AML.whatAtoms(config);
    }

    public void actionPerformed(ActionEvent e) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */
        String arg = (String)e.getActionCommand();
        if (arg.equals("Show")) {
            cutoff=Double.parseDouble(cut.getText());
            targets=target.getText();
            showclus();
        }
        else if (arg.equals("Load")) {
            config=getConfig(home,"CFG");
            if(config != null && config.natms > 0)
                AML.whatAtoms(config);
            if(!editor.isVisible())
                editor.showEditor();
            editor.pane.restore();
        }
        else if (arg.equals("Close")) {
            job.dispose();
        }
    }
    void showclus() {
        /*
**********************************************************************

dl_poly/java utility to show atomic clusters

copyright daresbury laboratory

author w. smith august 2011

**********************************************************************
         */
        boolean safe;
        String newtarget;
        int n,k,r,g,b;

        safe=false;

        if(config == null || config.natms == 0) {
            println("Error - no configuration loaded");
            return;
        }

        // Identify target atoms

        knd=new int[config.natms];

        for(int i=0;i<config.natms;i++)
            knd[i]=0;

        num_selected=0;
        for(int j=0;j<BML.countWords(targets);j++) {
            newtarget=BML.giveWord(targets,j+1);
            for(int i=0;i<config.natms;i++) {
                if(newtarget.equals(BML.giveWord(config.atoms[i].zsym,1))) {
                   knd[i]=1;
                   safe=true;
                   num_selected++;
                }
            }
        }

        // Exit if no required atoms found

        if(!safe) {
            println("Error - no atoms of requested type(s) found");
            return;
        }

        // Identify clusters in configuration

        clusters();

        // Exit of no clusters found

        if(max_size < 2) {
            println("Error - no clusters larger than one atom found");
            return;
        }

        // Colour code atoms not selected for cluster analysis

        for (int i=0;i<config.natms;i++) {
            if(knd[i]==0) {
                config.atoms[i].zcol=Color.lightGray;
                config.atoms[i].dotify=true;
            }
            else
                config.atoms[i]=new Element(config.atoms[i].zsym);
        }

        // Set cluster colours

        shade=new Color[num_clusters];
        for(int i=0;i<num_clusters;i++) {
            r=(int)(256*Math.random());
            g=(int)(256*Math.random());
            b=(int)(256*Math.random());
            shade[i]=new Color(r,g,b);
        }

        // Colour code atoms selected for cluster analysis

        n=0;
        config.atoms[loc[lok[0]]].zcol=shade[n];
        for(int i=1;i<num_selected;i++) {
            if(key[i] != key[i-1]) n++;
            config.atoms[loc[lok[i]]].zcol=shade[n];
        }

        // Backup copy of new config

        cfgsav=copyConfig(config);

        // Show clusters

        if(!editor.isVisible())
            editor.showEditor();
        editor.pane.restore();
    }

    void clusters() {
        /*
*********************************************************************

dl_poly/java GUI routine to determine clusters in a configuration

copyright - daresbury laboratory
author    - w.smith 2011

*********************************************************************
         */

        double[] vvv;
        double cut2;
        int j,m,n,p,q,ii,jj,lim,npm1,npm2;
        key=new int[num_selected];
        lok=new int[num_selected];
        loc=new int[num_selected];

        cut2=cutoff*cutoff;

        // Initialise arrays

        n=0;
        for (int i=0;i<config.natms;i++) {
            if(knd[i] > 0) {
                loc[n]=i;
                lok[n]=n;
                key[n]=n;
                n++;
            }
        }

        // Identify clusters

        lim=num_selected;
        npm1=lim/2;
        npm2=(lim-1)/2;
        vvv=new double[3];
        for(int k=1;k<=npm1;k++) {
            if(k>npm2)lim=npm1;
            for (int i=0;i<lim;i++) {
                j=i+k;
                if(j>=num_selected) j-=num_selected;
                ii=loc[i];
                jj=loc[j];
                vvv[0]=config.xyz[0][ii]-config.xyz[0][jj];
                vvv[1]=config.xyz[1][ii]-config.xyz[1][jj];
                vvv[2]=config.xyz[2][ii]-config.xyz[2][jj];
                config.pbc.images(vvv);
                if(cut2 >= vvv[0]*vvv[0]+vvv[1]*vvv[1]+vvv[2]*vvv[2]){

                    m=Math.min(i,j);
                    n=Math.max(i,j);

                    if(key[n] == n) {
                        key[n]=key[m];
                    }
                    else if (key[n] != key[m]) {

                        p=Math.min(key[m],key[n]);
                        q=Math.max(key[m],key[n]);

                        for (int u=0;u<num_selected;u++) {
                            if(key[u] == q) key[u]=p;
                        }
                    }
                }

            }
        }

        // Sort clusters in ascending order

        AML.ShellSort(num_selected,lok,key);

        // Find number of clusters

        n=1;
        for(int i=1;i<num_selected;i++) {
            if(key[i] != key[i-1])n++;
        }
        num_clusters=n;
        println("The number of clusters found is "+BML.fmt(n,6));

        // Find cluster sizes

        size=new int[num_clusters];
        m=0;
        n=1;
        min_size=num_selected;
        max_size=1;
        for(int i=1;i<num_selected;i++) {
            if(key[i] != key[i-1]) {
                min_size=Math.min(min_size,n);
                max_size=Math.max(max_size,n);
                size[m]=n;
                m++;
                n=1;
            }
            else
                n++;
        }
        min_size=Math.min(min_size,n);
        max_size=Math.max(max_size,n);
        size[m]=n;
        println("The  largest cluster found is "+BML.fmt(max_size,6));
        println("The smallest cluster found is "+BML.fmt(min_size,6));

        // Calculate population distribution of sizes

        range=new int[max_size-min_size+1];
        distribution=new int[max_size-min_size+1];

        println("Size Distribution"+"\n"+"  size     pop");

        for(int i=0;i<max_size-min_size+1;i++) {
            range[i]=min_size+i;
            distribution[i]=0;
        }

        for(int i=0;i<num_clusters;i++)
            distribution[size[i]-min_size]+=1;

        for(int i=0;i<max_size-min_size+1;i++) {
            if(distribution[i] > 0)
                println(BML.fmt(range[i],6)+BML.fmt(distribution[i],6));
        }
    }

}
