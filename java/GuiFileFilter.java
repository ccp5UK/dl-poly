import java.io.File;
import javax.swing.filechooser.*;

class GuiFileFilter extends FileFilter {
        /*
*********************************************************************

dl_poly/java GUI class for defining file filters

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
         */
    int action;
    String filter;
    String filename;
    String description;

    final int FILENAME_HEAD=0;
    final int FILENAME_TAIL=1;
    final int FILENAME_CONTAINS=2;

    public GuiFileFilter(int act) {
	super();
        action=act;
    }
    public boolean accept(File f) {
        if(f.isDirectory()) return true;

        filename = f.getName();

        if(action == FILENAME_CONTAINS) {
            if(filename.indexOf(filter.toUpperCase())>=0 || filename.indexOf(filter.toLowerCase())>=0) {
                return true;
            }
        }
        else if(action == FILENAME_HEAD) {
            if(filename.startsWith(filter.toUpperCase()) || filename.startsWith(filter.toLowerCase())) {
                return true;
            }
        }
        else if(action == FILENAME_TAIL) {
            if(filename.endsWith(filter.toUpperCase()) || filename.endsWith(filter.toLowerCase())) {
                return true;
            }
        }
        return false;
    }
    public String getDescription() {
        if(filter.equals("CFG"))
            description="CONFIG Files";
        else if(filter.equals("FLD"))
            description="FIELD Files";
        else if(filter.equals("XYZ"))
            description="XYZ Files";
        else if(filter.equals("PDB"))
            description="PDB Files";
        else if(filter.equals("MSI"))
            description="CERIUS_2 Files";
        else if(filter.equals("CNT"))
            description="CONTROL Files";
        else if(filter.equals("TAB"))
            description="TABLE Files";
        else if(filter.equals("XY"))
            description="XY Plot Files";
        else if(filter.equals("HOVG"))
            description="van Hove Plot Files";
        else if(filter.equals("DEN"))
            description="S(k,w) Plot Files";
        return description;
    }
}
