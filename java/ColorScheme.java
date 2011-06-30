import java.awt.*;

public class ColorScheme {
        /*
**********************************************************************
         
dl_poly/java class to define GUI colour scheme
         
copyright - daresbury laboratory
author    - w. smith december 2000
         
**********************************************************************
         */
    public Color back,fore,scrn,butn,scrf,butf;
    
    public ColorScheme(String scheme) {
        /*
*********************************************************************
         
dl_poly/java GUI routine
         
copyright - daresbury laboratory
author    - w.smith 2000
         
*********************************************************************
         */
        
        if(scheme.equals("picasso")) {
            back=new Color(80,80,220);
            fore=new Color(255,255,255);
            scrn=new Color(240,240,255);
            butn=new Color(128,128,230);
            scrf=new Color(0,0,0);
            butf=new Color(0,0,0);
        }
        else if(scheme.equals("monet")) {
            back=new Color(176,109,216);
            fore=new Color(255,255,255);
            scrn=new Color(222,194,239);
            butn=new Color(0,64,128);
            scrf=new Color(0,0,0);
            butf=new Color(255,255,255);
        }
        else if(scheme.equals("vangogh")) {
            back=new Color(255,255,0);
            fore=new Color(0,0,0);
            scrn=new Color(255,255,191);
            butn=new Color(255,128,64);
            scrf=new Color(0,0,0);
            butf=new Color(0,0,0);
        }
        else if(scheme.equals("cezanne")) {
            back=new Color(0,128,0);
            fore=new Color(255,255,255);
            scrn=new Color(193,255,193);
            butn=new Color(119,60,0);
            scrf=new Color(0,0,0);
            butf=new Color(255,255,255);
        }
        else if(scheme.equals("mondrian")) {
            back=new Color(0,0,255);
            fore=new Color(0,0,0);
            scrn=new Color(255,255,0);
            butn=new Color(255,0,0);
            scrf=new Color(0,0,0);
            butf=new Color(0,0,0);
        }
        else {
            // Picasso as default
            
            back=new Color(80,80,220);
            fore=new Color(0,0,0);
            scrn=new Color(205,205,245);
            butn=new Color(128,128,230);
            scrf=new Color(0,0,0);
            butf=new Color(0,0,0);
        }
    }
}
