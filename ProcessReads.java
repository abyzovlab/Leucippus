// package packlh000;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Vector;

/**
 * 
 * @author :
 * test svn 04/25/2016
 * ProcessReads is a simple class that contains 'getNameCoNum', 'reverseComplement', and 'reverse' methods
 * The 'getNameId' method is not currently used.
 * It actually retrieves read name-type, it complement-reverses a read, and it reverses the base quality String  
 * 
 */
public class ProcessReads 
{

    /**
     * Method that returns a unique part of the actual FASTQ read name
     * It was created to reduce the amount of final calculations. (less String length when comparison
     * in search is implemented	
     **/
    public String getNameId(String name)
    {
	String nms="";
	String[] r1wrdss = name.split(" " );
	String[] r1wrdsp = r1wrdss[0].split(":");
	nms = r1wrdsp[r1wrdsp.length-3] + ":" +  r1wrdsp[r1wrdsp.length-2] + ":" + r1wrdsp[r1wrdsp.length-1];
	return nms;
    }
    /**
     *
     * Method that accepts a read name and returns only the first unique part of it. Info at the end of the read and 
     * the type of read is missing. The method takes care of two different formats of read names.
     * @param ln : String the entire read name in the field 
     * @return : auxi : Simple read name.
     */
    public String getSimpleNameBack(String ln)
    {
	String auxi="";
	ln = ln.replaceAll("[/]", " ");
	
	if (ln.length()>1)
	    {
		String[] aux=ln.split(" ");
		if (aux.length>1)
		    auxi = aux[0].substring(1, aux[0].length());
		else
		    //auxi = aux[0].substring(1, aux[0].length()-2);
		    auxi = aux[0].substring(1, aux[0].length());
	    }
	return auxi;
    }

    public String getTypeBack(String ln)
    {
	String type="";
	ln = ln.replaceAll("[/]", " ");
	String[] aux=ln.split(" ");
	if (aux.length > 1)
	    type= Character.toString(aux[1].charAt(0));
	if (aux.length == 1) {
	    //type=Character.toString(aux[0].charAt(aux[0].length()-1));
	    type="";
	}
	// System.out.println(type);
	return type;
    }
    

    /**
     *
     * Method that accepts a read name and returns only the first unique part of it. Info at the end of the read and 
     * the type of read is missing. The method takes care of two different formats of read names.
     * The 'getSimpleName' and getType methods work well for the following read name and type formats
     * 	String name1 = "@MISEQ:56:000000000-A34:1:2134345:9234435:24156";  type = no type
     * 	String name2 = "@MISEQ:45:000000000-ATR:1:12346:1084534:18356331 2:N:0:GTTTTCCGC";  type = 2
     * 	String name3 = "@MISEQ:61:000000000-AHJK:1:23431Y:9741:18133461 2:N:0:GTCCGCCC";  type = 2
     * 	String name4 = "@MISEQ:35:000000000-AZAX:1:23145:923530:2462353156 1:N:0:ACAGAATG"; type = 1
     * 	String name5 = "@MISEQ:89:000000000-A45ER:1:232345:92341:27:2"; type = 2
     * 	String name6 = "@MISEQ:32:000000000-A123421:1:2115:93543:2413458/2"; type = 2
     * 
     * @param ln : String the entire read name in the field 
     * @return : auxi : Simple read name.
     */
    public String getSimpleName(String ln, String separator)
    {
	String auxi="";
	char sep='l';
	if (separator.isEmpty()) {
	    if (ln.length() > 1) {
		String[] aux=ln.split(" ");
		if (aux.length>1)
		    auxi = aux[0].substring(1, aux[0].length());
		else {
		    if ( (ln.charAt(ln.length()-2)=='/') || (ln.charAt(ln.length()-2)==':') )
			auxi = ln.substring(1, ln.length()-2);
		    else
			auxi = ln.substring(1, ln.length());
		}
	    }
	} else {
	    if (separator.length() == 1) {
		sep=separator.charAt(0);
		if (ln.charAt(ln.length()-2) == sep)
		    auxi = ln.substring(1, ln.length()-2);
		else
		    auxi = ln.substring(1, ln.length());
	    }
	}
	return auxi;
    }
    /**
     * See getSimpleName info
     * @param ln
     * @return
     */
    public String getType(String ln, String separator)
    {
	String type="";
	char sep='l';
	if (separator.isEmpty()) {
	    if (ln.length() > 1) {
		String[] aux=ln.split(" ");
		if (aux.length > 1)
		    type= Character.toString(aux[1].charAt(0));
		else
		    if (ln.charAt(ln.length()-2) == '/' || ln.charAt(ln.length() - 2) == ':')
			type = Character.toString(ln.charAt(ln.length() - 1));
		    else
			type="";
	    }
	} else {
	    if (separator.length() == 1) {
		sep=separator.charAt(0);
		if (ln.charAt(ln.length()-2) == sep)
		    type = Character.toString(ln.charAt(ln.length() - 1));
	    }
	}
	return type;
    }


    /**
     * 'getNameCoNum' method receives the entire String that is the read's name 
     * Next it identifies its format, and it retrieves the actual read name and the
     * type of the read. It creates and returns a comma separated String that contains
     * the actual read name and the  read type.
     * 08/21/2015 : Reads 2 were replaced with Reads3 in a particular data
     * 08/21/2015 : the method where modified to replace 3 by 2
     * @param ln : initial mixed read name
     * @return auxi : actual-read-name + "," + read-type
     */
    public String getNameCoNum(String ln)
    {
	String auxi="";
	String[] aux=ln.split(" ");
	if (aux.length>1)
	    {
		if (aux[1].charAt(0)=='3')
		    auxi = aux[0].substring(1, aux[0].length()) + "," + "2";
		else
		    auxi = aux[0].substring(1, aux[0].length()) + "," + aux[1].charAt(0);
		
	    }
	else
	    {
		if (aux[0].charAt(aux[0].length()-1)=='3')
		    auxi = aux[0].substring(1, aux[0].length()) + "," + "2";
		else
		    auxi = aux[0].substring(1, aux[0].length()-2) + "," + aux[0].charAt(aux[0].length()-1);
	    }	
	return auxi;
    }

    /**
     * 'reverseComplement' creates a reverse complement of the particular read.
     * 'fastq' format provides reads with the actual sequence. This read can be derived from reverse strand.
     *  The demand is : all reads must have a format like they are read from left to right.
     *  So the method accepts the read from the reverse strand reverses it properly and returns it.
     *  
     * @param rrd
     * @return
     */
    public static String reverseComplement(String rrd)
    {
	String resu="";
	char[] charrread = rrd.toCharArray();
	char[] res = new char[charrread.length];
	int count=-1;
	
	for (int i=charrread.length-1; i>-1; i--)
	    {
		count = count+1;
		res[count]=charrread[i];
	    }
	
	for(int i=0; i< res.length; i++)
	    {
		if (res[i]=='A')
		    res[i]='t';
		if (res[i]=='T')
		    res[i]='a';
		if (res[i]=='C')
		    res[i]='g';
		if (res[i]=='G')
		    res[i]='c';
	    }
	
	for(int i=0; i< res.length; i++)
	    {
		if (res[i]=='t')
		    res[i]='T';
		if (res[i]=='a')
		    res[i]='A';
		if (res[i]=='c')
		    res[i]='C';
		if (res[i]=='g')
		    res[i]='G';
	    }		
	
	resu = new String(res);
	return resu;
    }
    
    public static String reverse(String str)
    {
	String res="";
	char[] charstr = str.toCharArray();
	char[] aux = new char[charstr.length];
	int count=-1;
	for (int i=charstr.length-1; i>-1; i--)
	    {
		count = count+1;
		aux[count]=charstr[i];
	    }
	res = new String(aux);
	return res;
    }
}
