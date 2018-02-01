// package packlh000;
// Remainder change hash tables to put as key the entire read name
public class Readob
{
    String namefull, name, seq, qual;

    public Readob(String namefull, String name, String seq, String qual)
    {
	super();
	
	int l = seq.length();
	while (l > 0) {
	    char c = qual.charAt(l - 1);
//  	    System.out.println(c + "  " + (Integer.valueOf(c) - 33));
	    if (Integer.valueOf(c) - 33 > 2) break;
	    l--;
	}
	
	this.namefull = namefull;
	this.name     = name;
	this.seq      = seq.substring(0,l);
	this.qual     = qual.substring(0,l);
    }

    public Readob() 
    {
	super();
	this.name = "";
	this.seq = "";
	this.qual = "";
	this.namefull = "";
    }
	
    public String getNamefull()
    {
	return namefull;
    }

    public void setNamefull(String namefull)
    {
	this.namefull = namefull;
    }

    public String getName()
    {
	return name;
    }

    public void setName(String name)
    {
	this.name = name;
    }

    public String getSeq() {
	return seq;
    }
    
    public void setSeq(String seq) {
	this.seq = seq;
    }
    
    public String getQual() {
	return qual;
    }
    
    public void setQual(String qual) {
	this.qual = qual;
    }
    
    @Override public String toString()
    {
	return "\nread object:\ncomplete name=" + namefull + "\nname=" + name + "\nseq=" + seq + "\nqual=" + qual+"\n";
    }
}
