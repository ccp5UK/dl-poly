public abstract class BML {
        /*
*********************************************************************

dl_poly/java GUI basic methods library

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */

    static boolean testWord(String word) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        boolean b=true;
        for(int i=0;i<word.length();i++) {
            byte j=(byte)word.charAt(i);
            if(!((j > 32 && j < 59) || (j > 96 && j < 123)))b=false ;
        }
        return b;
    }
    static int countWords(String text) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        int n=0;
        boolean gap=true;

        for (int i=0;i<text.length();i++) {
            if((text.charAt(i) != ' ') && (text.charAt(i) != ',')) {
                if(gap) {
                    gap=false;
                    n++;
                }
            }
            else {
                gap=true;
            }
        }
        return n;
    }
    static String giveWord(String text, int n) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */

        int k=0,a=0,b=0,c=0;
        String word="";
        char u;
        while (k < n && c < text.length()) {
            u=text.charAt(c);
            if (u==',' || u==' ' || u=='"' || u=='(' || u==')') {
                if (b>a) k++;
                if (k<n) a=c+1;
            }
            else {
                b=c+1;
            }
            c++;
        }
        if (k == n) {
            word=text.substring(a,b);
        }
        else {
            word=text.substring(a);
        }
        return word;
    }
    static int giveInteger(String text, int n) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        char u;
        int k=0,a=0,b=0,c=0,d=0;
        while (k < n && c < text.length()) {
            u=text.charAt(c);
            if (u==',' || u==' ' || u=='"' || u=='(' || u==')') {
                if (b>a) k++;
                if (k<n) a=c+1;
            }
            else {
                b=c+1;
            }
            c++;
        }
        if (k == n) {
            d=Integer.parseInt(text.substring(a,b));
        }
        else {
            d=Integer.parseInt(text.substring(a));
        }
        return d;
    }
    static double giveDouble(String text, int n) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        int k=0,a=0,b=0,c=0;
        double d;
        char u;
        Double dub;
        while (k < n && c < text.length()) {
            u=text.charAt(c);
            if (u==',' || u==' ' || u=='"' || u=='(' || u==')') {
                if (b>a) k++;
                if (k<n) a=c+1;
            }
            else {
                b=c+1;
            }
            c++;
        }
        if (k == n) {
            dub = new Double(text.substring(a,b));
            d=dub.doubleValue();
        }
        else {
            dub = new Double(text.substring(a));
            d=dub.doubleValue();
        }
        return d;
    }

    // Format a double into a fixed length string

    static String fmt(double a, int n) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        String f="00000000000000000000";
        String t="  ******************";
        String s=Double.toString(a);
        int e=Math.max(s.indexOf("e"),s.indexOf("E"));
        if(e > 0) {
            String u=s.substring(e);
            int k=u.length();
            if(a >= 0)
                s="  "+s.substring(0,e)+f;
            else
                s=" "+s.substring(0,e)+f;
            return  s.substring(0,n-k)+u;
        }
        else {
            if(s.indexOf(".") < 0)s=s+".";
            if(s.indexOf(".") > n-2)
                s=t;
            else {
                if(a >= 0)
                    s="  "+s+f;
                else
                    s=" "+s+f;
            }
            return s.substring(0,n);
        }
    }

    // Format an int into a fixed length string

    static String fmt(int j, int n) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        String f="                    ";
        String t=String.valueOf(j);
        int k=Math.min(n,t.length());
        return  (f.substring(0,n-k)+t.substring(0,k));
    }

    // Format a string in a specified number of characters

    static String fmt(String a, int n) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        String pad = "                    ";
        String word = a+pad;
        return word.substring(0,n);
    }

    // Minimum of three real numbers

    static double min(double a, double b, double c) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        double z=a;
        z=Math.min(z,b);
        z=Math.min(z,c);
        return z;
    }

    // Maximum of three real numbers

    static double max(double a, double b, double c) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
        double z=a;
        z=Math.max(z,b);
        z=Math.max(z,c);
        return z;
    }

    // Nearest integer function

    static int nint(double a) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        int k=(int)a;
        if(a >= 0) {
            if(a-k >= 0.5)k++;
        }
        else {
            if(k-a >= 0.5)k--;
        }
        return k;
    }

    // Sign of number

    static double sign(double a) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
         */
        if(a>0.0)
            return 1.0;
        else if(a<0.0)
            return -1.0;
        return 0.0;
    }

    // calculate a normalised difference vector

    static double[] ndxyz(int i,int j,double xyz[][]) {
        /*
*********************************************************************

dl_poly/java GUI routine to calculate normalised xyz difference vector

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
         */
        double rrr;
        double uuu[]=new double[3];

        uuu[0]=xyz[0][i]-xyz[0][j];
        uuu[1]=xyz[1][i]-xyz[1][j];
        uuu[2]=xyz[2][i]-xyz[2][j];
        rrr=1.0/Math.sqrt(uuu[0]*uuu[0]+uuu[1]*uuu[1]+uuu[2]*uuu[2]);
        uuu[0]*=rrr;
        uuu[1]*=rrr;
        uuu[2]*=rrr;

        return uuu;
    }

    // normalise a vector

    static void vnorm(double vvv[]) {
        /*
*********************************************************************

dl_poly/java GUI routine

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
         */
        double rrr;

        rrr=1.0/Math.sqrt(vvv[0]*vvv[0]+vvv[1]*vvv[1]+vvv[2]*vvv[2]);
        vvv[0]*=rrr;
        vvv[1]*=rrr;
        vvv[2]*=rrr;

    }

    // vector cross product

    static double[] cross(double uuu[],double vvv[]) {
        /*
*********************************************************************

dl_poly/java GUI routine vector cross product

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
         */

        double www[]=new double[3];

        www[0]=uuu[1]*vvv[2]-uuu[2]*vvv[1];
        www[1]=uuu[2]*vvv[0]-uuu[0]*vvv[2];
        www[2]=uuu[0]*vvv[1]-uuu[1]*vvv[0];

        return www;
    }
}


