// package packlh000;
// Remainder change hash tables to put as key the entire read name
public class Readob
{
	String namefull, name, seq, qual;
//, nm;

	public Readob(String namefull, String name, String seq, String qual)
//, String nm) 
	{
		super();
		this.namefull = namefull;
		this.name = name;
		this.seq = seq;
		this.qual = qual;
		;
		//this.nm = nm;
	}
	public Readob() 
	{
		super();
		this.name = "";
		this.seq = "";
		this.qual = "";
		this.namefull = "";
		//this.nm = "";
	}
	
	public String getNamefull() {
		return namefull;
	}

	public void setNamefull(String namefull) {
		this.namefull = namefull;
	}

	public String getName() {
		return name;
	}

	//public String getNm() {
		//return nm;
	//}
	//public void setNm(String nm) {
	//	this.nm = nm;
	//}
	public void setName(String name) {
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

	@Override
	public String toString() {
		return "\nread object:\ncomplete name=" + namefull + "\nname=" + name + "\nseq=" + seq + "\nqual=" + qual+"\n";
				// + ", nm=" + nm + "]";
	}
}
