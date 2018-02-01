# source("/home/m026918/Leucippus/Pvalue_Mutrate.r")
# p-value-and-mutation-rate-rscript

# arguments
source = "";		# input table(s) one table or two for comparison
			# case : one table : graph or posgraph
			# case two tables  : graph
			# if graph then use regular noise calculation
			# if pos graph use only some restrictions
			# there is not input for germline border in posgraph or 
			# graph inclusion for noise : #reference/total)>

destination = ""; 	# output graph file
graphortable = "";	# remnand from previousp value option(used particular loop)
#overlaps = "";		# restricts xval, yval to be greater than overlap (string)
overlap = "";		# restricts xval, yval to be greater than overlap (string)
coverage = "";		# user customizable restriction (seems it is working)
ranges="";
range = 0;		# numeric variable for error cutoff(#alternative/total)
			# reforms the shape(x inclusion and extension)
type = "";		# pvalue or mutation rate
curarg="";		# in followng loop : stores curent argument as a string 
curargvc=c("");		# current argument array : is used when agument contains
			# two variables (split when "=" is present in name and value)
cover=0;		# user customizable numeric variable for coverage-restiction.
  
prf = 0;
args <- commandArgs(trailingOnly = TRUE)


# source = "";
# destination = "";
# graphortable = "";
# overlap = "";
# coverage = "";
# prf = 0;
# type = "";
# curarg="";
# curargvc=c("");
# cover=0;
# ranges=""
# range=0

#args <- commandArgs(trailingOnly = TRUE)
#print("here");
#print(args[1]);
#print(args[2]);
#print(args[3]);
#print(args[4]);

for (i in 1:length(args))
{
	curarg = args[i];
	curargvc=strsplit(curarg, "=")
	curargvc=curargvc[[1]];
	if(length(curargvc)==2)
	{

		if(curargvc[1]=="src")
		{
			source=curargvc[2];
			print(source);
		}

		if(curargvc[1]=="dst")
		{
			destination=curargvc[2];
			print(destination);
		}

		if(curargvc[1]=="tp")
		{
			type=curargvc[2];
			print(type);
		}

		if(curargvc[1]=="ov")
		{
			overlap=curargvc[2];
			print(overlap);
		}

		if(curargvc[1]=="gt")
		{
			graphortable=curargvc[2];
			print(graphortable);
		}

		if(curargvc[1]=="cov")
		{
			coverage=curargvc[2];
			cover = strtoi(coverage);
			print(coverage);
		}

		if(curargvc[1]=="gla")
		{
			germlineAFs=curargvc[2];
			# prf = as.numeric(germlineAFs);
			print(prf);
		}

		if(curargvc[1]=="rng")
		{
			ranges=curargvc[2];
			range = as.numeric(ranges);
			prf = as.numeric(ranges);

			print(paste("range : ",  range, sep=""))
		}
	}
}

# range = prf
#prf=0.95
if (range==0)
range = 0.005
cover=cover-1
oneortwo=0;
# print(graphortable)

# put all table names to a vector

splat <- strsplit(source, ",")[[1]]
tblnms=c();
nmOfiles=0;

for(i in 1:length(splat))
{
	tblnms = c(tblnms, basename(splat[i]))
}

nmOfiles=length(splat);

if(nmOfiles == 2)
{
	print("Table files are two!");
	oneortwo=2;
}

if(nmOfiles == 1)
{
	print("Only one Table file.");
	oneortwo=1;
}

#create a list vector to hold the data frames of each table as its elements
print(length(splat));
vectorOfTables <- vector(mode = "list", length = length(splat));

# populate list-data-frames
# Error generated because R converts String Values of Single Capital Ts to logic TRUE
# Now all data frame values initial are converted to characters and later the numeric 
# columns are converted to numeric values  
for(i in 1:length(splat))
{
	dfrm <- read.table(splat[i],sep="\t",na.strings = "NA",header = TRUE, colClasses = "character")
	vectorOfTables[[i]] <- dfrm
}

i=0
xl=c(-0.0001,range)
yl=c(0,10000)
#brk <- seq(0,0.005, 6.25e-5)
brk <- seq(0,range, 6.25e-5)
xval=overlap
yval=overlap
j=0;
i=1;
lf=1;
# exclude site
# Reference  > 95%
# Quality 0 to 50 AC AT AG, TG TA TC, CA CT CG, GT GA GC

if ((graphortable=="groutfpval") | (graphortable=="grout") | (graphortable=="posgrout"))
{
	pdf(paste(destination, ".pdf", sep=""), width=11, height=8.5)
}
#setwd("C:\\Users\\m026918\\Desktop");
#prf=0.95 # percentage of reference
#prf=0.55 # percentage of reference
# Ctreate Comparison empirical pvalue graphs
if( (oneortwo==2) & (type=="pvalue") & ( (graphortable=="groutfpval") | (graphortable=="grout") ) )
{
#----------------
	prf=0.95
#----------------
	old.par <- par(mfrow=c(1, 1))
	print("comparison pvalue graph 2");
	dfrm <- vectorOfTables[[1]]
	nbe  <- as.numeric(dfrm$number_base_expected);
	dfrm <- subset(dfrm, nbe>0);
	chr  <- dfrm$chromosome;
	be   <- dfrm$base_expected;
	alt  <- dfrm$alternative;
	site <- as.numeric(dfrm$site);
	tA   <- as.numeric(dfrm$number_of_As);
	tC   <- as.numeric(dfrm$number_of_Cs);
	tT   <- as.numeric(dfrm$number_of_Ts);
	tG   <- as.numeric(dfrm$number_of_Gs);
	nbe  <- as.numeric(dfrm$number_base_expected);
	xv   <- as.numeric(dfrm$x);
	yv   <- as.numeric(dfrm$y);
	dv = nbe;
	print("table 1");

	sbt_ac <- subset(tC/dv, (be=='A' & alt=='A' & tA/dv>prf & xv>xval & yv>yval & tC/dv < range & dv>cover)) # red
	sbt_at <- subset(tT/dv, (be=='A' & alt=='A' & tA/dv>prf & xv>xval & yv>yval & tT/dv < range & dv>cover)) # blue 
	sbt_ag <- subset(tG/dv, (be=='A' & alt=='A' & tA/dv>prf & xv>xval & yv>yval & tG/dv < range & dv>cover)) # green

	sbt_tg <- subset(tG/dv, (be=='T' & alt=='T' & tT/dv>prf & xv>xval & yv>yval & tG/dv < range & dv>cover)) # purple
	sbt_ta <- subset(tA/dv, (be=='T' & alt=='T' & tT/dv>prf & xv>xval & yv>yval & tA/dv < range & dv>cover)) # cyan
	sbt_tc <- subset(tC/dv, (be=='T' & alt=='T' & tT/dv>prf & xv>xval & yv>yval & tC/dv < range & dv>cover)) # orange

	sbt_ca <- subset(tA/dv, (be=='C' & alt=='C' & tC/dv>prf & xv>xval & yv>yval & tA/dv < range & dv>cover)) # red 
	sbt_ct <- subset(tT/dv, (be=='C' & alt=='C' & tC/dv>prf & xv>xval & yv>yval & tT/dv < range & dv>cover)) # blue 
	sbt_cg <- subset(tG/dv, (be=='C' & alt=='C' & tC/dv>prf & xv>xval & yv>yval & tG/dv < range & dv>cover)) # green

	sbt_gt <- subset(tT/dv, (be=='G' & alt=='G' & tG/dv>prf & xv>xval & yv>yval & tT/dv < range & dv>cover)) # purple
	sbt_ga <- subset(tA/dv, (be=='G' & alt=='G' & tG/dv>prf & xv>xval & yv>yval & tA/dv < range & dv>cover)) # cyan 
	sbt_gc <- subset(tC/dv, (be=='G' & alt=='G' & tG/dv>prf & xv>xval & yv>yval & tC/dv < range & dv>cover)) # orange

	sbt_ac.cut = cut(sbt_ac, brk, right=FALSE)
	sbt_ac.freq = table(sbt_ac.cut)
	cumfreq0 = c(0, cumsum(sbt_ac.freq))
	sumfreq0=sum(sbt_ac.freq)
	sbt_ac1 = cumfreq0/sumfreq0
	sbt_ac1_p = 1-cumfreq0/sumfreq0
	
	sbt_at.cut = cut(sbt_at, brk, right=FALSE)
	sbt_at.freq = table(sbt_at.cut)
	cumfreq0 = c(0, cumsum(sbt_at.freq))
	sumfreq0=sum(sbt_at.freq)
	sbt_at1 = cumfreq0/sumfreq0
	sbt_at1_p = 1-cumfreq0/sumfreq0
	
	sbt_ag.cut = cut(sbt_ag, brk, right=FALSE)
	sbt_ag.freq = table(sbt_ag.cut)
	cumfreq0 = c(0, cumsum(sbt_ag.freq))
	sumfreq0=sum(sbt_ag.freq)
	sbt_ag1 = cumfreq0/sumfreq0
	sbt_ag1_p = 1-cumfreq0/sumfreq0

	sbt_tg.cut = cut(sbt_tg, brk, right=FALSE)
	sbt_tg.freq = table(sbt_tg.cut)
	cumfreq0 = c(0, cumsum(sbt_tg.freq))
	sumfreq0=sum(sbt_tg.freq)
	sbt_tg1 = cumfreq0/sumfreq0
	sbt_tg1_p = 1-cumfreq0/sumfreq0

	sbt_ta.cut = cut(sbt_ta, brk, right=FALSE)
	sbt_ta.freq = table(sbt_ta.cut)
	cumfreq0 = c(0, cumsum(sbt_ta.freq))
	sumfreq0=sum(sbt_ta.freq)
	sbt_ta1 = cumfreq0/sumfreq0
	sbt_ta1_p = 1-cumfreq0/sumfreq0

	sbt_tc.cut = cut(sbt_tc, brk, right=FALSE)
	sbt_tc.freq = table(sbt_tc.cut)
	cumfreq0 = c(0, cumsum(sbt_tc.freq))
	sumfreq0=sum(sbt_tc.freq)
	sbt_tc1 = cumfreq0/sumfreq0
	sbt_tc1_p = 1-cumfreq0/sumfreq0

	sbt_ca.cut = cut(sbt_ca, brk, right=FALSE)
	sbt_ca.freq = table(sbt_ca.cut)
	cumfreq0 = c(0, cumsum(sbt_ca.freq))
	sumfreq0=sum(sbt_ca.freq)
	sbt_ca1 = cumfreq0/sumfreq0
	sbt_ca1_p = 1-cumfreq0/sumfreq0

	sbt_ct.cut = cut(sbt_ct, brk, right=FALSE)
	sbt_ct.freq = table(sbt_ct.cut)
	cumfreq0 = c(0, cumsum(sbt_ct.freq))
	sumfreq0=sum(sbt_ct.freq)
	sbt_ct1 = cumfreq0/sumfreq0
	sbt_ct1_p = 1-cumfreq0/sumfreq0

	sbt_cg.cut = cut(sbt_cg, brk, right=FALSE)
	sbt_cg.freq = table(sbt_cg.cut)
	cumfreq0 = c(0, cumsum(sbt_cg.freq))
	sumfreq0=sum(sbt_cg.freq)
	sbt_cg1 = cumfreq0/sumfreq0
	sbt_cg1_p = 1-cumfreq0/sumfreq0

	sbt_gt.cut = cut(sbt_gt, brk, right=FALSE)
	sbt_gt.freq = table(sbt_gt.cut)
	cumfreq0 = c(0, cumsum(sbt_gt.freq))
	sumfreq0=sum(sbt_gt.freq)
	sbt_gt1 = cumfreq0/sumfreq0
	sbt_gt1_p = 1-cumfreq0/sumfreq0

	sbt_ga.cut = cut(sbt_ga, brk, right=FALSE)
	sbt_ga.freq = table(sbt_ga.cut)
	cumfreq0 = c(0, cumsum(sbt_ga.freq))
	sumfreq0=sum(sbt_ga.freq)
	sbt_ga1 = cumfreq0/sumfreq0
	sbt_ga1_p = 1-cumfreq0/sumfreq0

	sbt_gc.cut = cut(sbt_gc, brk, right=FALSE)
	sbt_gc.freq = table(sbt_gc.cut)
	cumfreq0 = c(0, cumsum(sbt_gc.freq))
	sumfreq0=sum(sbt_gc.freq)
	sbt_gc1 = cumfreq0/sumfreq0
	sbt_gc1_p = 1-cumfreq0/sumfreq0

	dfrm <- vectorOfTables[[2]]
	nbe  <- as.numeric(dfrm$number_base_expected);
	dfrm <- subset(dfrm, nbe>0);
	chr  <- dfrm$chromosome;
	be   <- dfrm$base_expected;
	alt  <- dfrm$alternative;
	site <- as.numeric(dfrm$site);
	tA   <- as.numeric(dfrm$number_of_As);
	tC   <- as.numeric(dfrm$number_of_Cs);
	tT   <- as.numeric(dfrm$number_of_Ts);
	tG   <- as.numeric(dfrm$number_of_Gs);
	nbe  <- as.numeric(dfrm$number_base_expected);
	xv   <- as.numeric(dfrm$x);
	yv   <- as.numeric(dfrm$y);
	dv = nbe;

	print("table 2");

	sbt_ac <- subset(tC/dv, (be=='A' & alt=='A' & tA/dv>prf & xv>xval & yv>yval & tC/dv < range & dv>cover)) # red
	sbt_at <- subset(tT/dv, (be=='A' & alt=='A' & tA/dv>prf & xv>xval & yv>yval & tT/dv < range & dv>cover)) # blue 
	sbt_ag <- subset(tG/dv, (be=='A' & alt=='A' & tA/dv>prf & xv>xval & yv>yval & tG/dv < range & dv>cover)) # green

	sbt_tg <- subset(tG/dv, (be=='T' & alt=='T' & tT/dv>prf & xv>xval & yv>yval & tG/dv < range & dv>cover)) # purple
	sbt_ta <- subset(tA/dv, (be=='T' & alt=='T' & tT/dv>prf & xv>xval & yv>yval & tA/dv < range & dv>cover)) # cyan
	sbt_tc <- subset(tC/dv, (be=='T' & alt=='T' & tT/dv>prf & xv>xval & yv>yval & tC/dv < range & dv>cover)) # orange

	sbt_ca <- subset(tA/dv, (be=='C' & alt=='C' & tC/dv>prf & xv>xval & yv>yval & tA/dv < range & dv>cover)) # red 
	sbt_ct <- subset(tT/dv, (be=='C' & alt=='C' & tC/dv>prf & xv>xval & yv>yval & tT/dv < range & dv>cover)) # blue 
	sbt_cg <- subset(tG/dv, (be=='C' & alt=='C' & tC/dv>prf & xv>xval & yv>yval & tG/dv < range & dv>cover)) # green

	sbt_gt <- subset(tT/dv, (be=='G' & alt=='G' & tG/dv>prf & xv>xval & yv>yval & tT/dv < range & dv>cover)) # purple
	sbt_ga <- subset(tA/dv, (be=='G' & alt=='G' & tG/dv>prf & xv>xval & yv>yval & tA/dv < range & dv>cover)) # cyan 
	sbt_gc <- subset(tC/dv, (be=='G' & alt=='G' & tG/dv>prf & xv>xval & yv>yval & tC/dv < range & dv>cover)) # orange

	sbt_ac.cut = cut(sbt_ac, brk, right=FALSE)
	sbt_ac.freq = table(sbt_ac.cut)
	cumfreq0 = c(0, cumsum(sbt_ac.freq))
	sumfreq0=sum(sbt_ac.freq)
	sbt_ac2 = cumfreq0/sumfreq0
	sbt_ac2_p = 1-cumfreq0/sumfreq0
	
	sbt_at.cut = cut(sbt_at, brk, right=FALSE)
	sbt_at.freq = table(sbt_at.cut)
	cumfreq0 = c(0, cumsum(sbt_at.freq))
	sumfreq0=sum(sbt_at.freq)
	sbt_at2 = cumfreq0/sumfreq0
	sbt_at2_p = 1-cumfreq0/sumfreq0
	
	sbt_ag.cut = cut(sbt_ag, brk, right=FALSE)
	sbt_ag.freq = table(sbt_ag.cut)
	cumfreq0 = c(0, cumsum(sbt_ag.freq))
	sumfreq0=sum(sbt_ag.freq)
	sbt_ag2 = cumfreq0/sumfreq0
	sbt_ag2_p = 1-cumfreq0/sumfreq0

	sbt_tg.cut = cut(sbt_tg, brk, right=FALSE)
	sbt_tg.freq = table(sbt_tg.cut)
	cumfreq0 = c(0, cumsum(sbt_tg.freq))
	sumfreq0=sum(sbt_tg.freq)
	sbt_tg2 = cumfreq0/sumfreq0
	sbt_tg2_p = 1-cumfreq0/sumfreq0

	sbt_ta.cut = cut(sbt_ta, brk, right=FALSE)
	sbt_ta.freq = table(sbt_ta.cut)
	cumfreq0 = c(0, cumsum(sbt_ta.freq))
	sumfreq0=sum(sbt_ta.freq)
	sbt_ta2 = cumfreq0/sumfreq0
	sbt_ta2_p = 1-cumfreq0/sumfreq0

	sbt_tc.cut = cut(sbt_tc, brk, right=FALSE)
	sbt_tc.freq = table(sbt_tc.cut)
	cumfreq0 = c(0, cumsum(sbt_tc.freq))
	sumfreq0=sum(sbt_tc.freq)
	sbt_tc2 = cumfreq0/sumfreq0
	sbt_tc2_p = 1-cumfreq0/sumfreq0

	sbt_ca.cut = cut(sbt_ca, brk, right=FALSE)
	sbt_ca.freq = table(sbt_ca.cut)
	cumfreq0 = c(0, cumsum(sbt_ca.freq))
	sumfreq0=sum(sbt_ca.freq)
	sbt_ca2 = cumfreq0/sumfreq0
	sbt_ca2_p = 1-cumfreq0/sumfreq0

	sbt_ct.cut = cut(sbt_ct, brk, right=FALSE)
	sbt_ct.freq = table(sbt_ct.cut)
	cumfreq0 = c(0, cumsum(sbt_ct.freq))
	sumfreq0=sum(sbt_ct.freq)
	sbt_ct2 = cumfreq0/sumfreq0
	sbt_ct2_p = 1-cumfreq0/sumfreq0

	sbt_cg.cut = cut(sbt_cg, brk, right=FALSE)
	sbt_cg.freq = table(sbt_cg.cut)
	cumfreq0 = c(0, cumsum(sbt_cg.freq))
	sumfreq0=sum(sbt_cg.freq)
	sbt_cg2 = cumfreq0/sumfreq0
	sbt_cg2_p = 1-cumfreq0/sumfreq0

	sbt_gt.cut = cut(sbt_gt, brk, right=FALSE)
	sbt_gt.freq = table(sbt_gt.cut)
	cumfreq0 = c(0, cumsum(sbt_gt.freq))
	sumfreq0=sum(sbt_gt.freq)
	sbt_gt2 = cumfreq0/sumfreq0
	sbt_gt2_p = 1-cumfreq0/sumfreq0

	sbt_ga.cut = cut(sbt_ga, brk, right=FALSE)
	sbt_ga.freq = table(sbt_ga.cut)
	cumfreq0 = c(0, cumsum(sbt_ga.freq))
	sumfreq0=sum(sbt_ga.freq)
	sbt_ga2 = cumfreq0/sumfreq0
	sbt_ga2_p = 1-cumfreq0/sumfreq0

	sbt_gc.cut = cut(sbt_gc, brk, right=FALSE)
	sbt_gc.freq = table(sbt_gc.cut)
	cumfreq0 = c(0, cumsum(sbt_gc.freq))
	sumfreq0=sum(sbt_gc.freq)
	sbt_gc2 = cumfreq0/sumfreq0
	sbt_gc2_p = 1-cumfreq0/sumfreq0


#ymaxi = max((h11$counts+h21$counts)/2, h12$counts, h22$counts, h13$counts, h23$counts, ho11$counts, ho21$counts, ho12$counts, ho22$counts, ho13$counts, ho23$counts)
	#ym = ymaxi + 0.1*ymaxi
	#yl=c(0,ym)
	
	#ymaxi = max(log2((sbt_ac1_p + sbt_tg1_p)/2), log2((sbt_ac1_po + sbt_ta1_po)/2), y=log2((sbt_at1_p + sbt_ta1_p)/2), log2((sbt_at1_po + sbt_ta1_po)/2), log2((sbt_ag1_p + sbt_tc1_p)/2), log2((sbt_ag1_po + sbt_tc1_po)/2), log2((sbt_ca1_p + sbt_gt1_p)/2), log2((sbt_ca1_po + sbt_gt1_po)/2), log2((sbt_ct1_p + sbt_ga1_p)/2), log2((sbt_ct1_po + sbt_ga1_po)/2), log2((sbt_cg1_p + sbt_gc1_p)/2), log2((sbt_cg1_po + sbt_gc1_po)/2))
	
	ymaxi = max(log10((sbt_ac1_p + sbt_tg1_p)/2), y=log10((sbt_at1_p + sbt_ta1_p)/2), log10((sbt_ag1_p + sbt_tc1_p)/2), log10((sbt_ca1_p + sbt_gt1_p)/2), log10((sbt_ct1_p + sbt_ga1_p)/2), log10((sbt_cg1_p + sbt_gc1_p)/2))
	
	
	ym = ymaxi + 0.1*ymaxi
	yl=c(-3,0)
	subi = paste("file1:  A->C & T->G = blue,   A->T & T->A = red   A->G & T->C = green   C->A & G->T = cyan,   C->T & G->A = purple,   C->G & G->C = orange \n", sep="")
	subia = paste(subi, "file2:  A->C & T->G = dashblue,   A->T & T->A = dashred   A->G & T->C = dashgreen   C->A & G->T = dashcyan,   C->T & G->A = dashpurple,   C->G & G->C = dashorange \n", sep="")
	
#subia = paste(subi, "i6_sep: A -> C & T -> G = d.blue,  A -> T & T -> A = d.red,   A -> G & T -> C = d.green \n C -> A & G -> T = d.cyan,  C -> T & G -> A = d.purple,  C -> G & G -> C = d.orange \n", sep="");
	
	#subi = paste(subia, qual[i], " < A -> C, T -> G = dot_blue,    A -> T, T -> A = dot_red,    A -> G, T -> C = dot_green \n", qual[i], " < C -> A, G -> T = dot_cyan,    C -> T, G -> A = dot_purple,    C -> G, G -> C = dot_orange", sep="")
	#ttl = "Normalized_Cumulative"
if (is.null(cover))
	cover=0;
	ttl  = paste("p_value: sites size>=", (cover+1), ", reference > ", (prf*100), "%,  y=log(p-value)",sep="");

	# p value graphs
	# AC TG >= 0 sbt_ac  jjjjjj sbt_tg blue
	plot(x=brk, y=(lf)*log10((sbt_ac1_p + sbt_tg1_p)/2), col="blue", lwd = 1, type="l", ylim=yl, xlim=xl, yaxt='n', xaxt="n", sub = subia, cex.sub = 0.70, font.sub = 1, xlab=NA, main=paste(ttl, "\nfile1: ", tblnms[1], "\nfile2: ", tblnms[2], sep=""), cex.main=0.70, ylab=paste("file1: ", tblnms[1], "  file2: ",  tblnms[2], sep="") )

	lines(x=brk, y=(lf)*log10((sbt_ac2_p + sbt_tg2_p)/2), col="blue", lty="twodash", lwd = 1)
	# AT TA >= 0 sbt_at sbt_ta red
	lines(x=brk, y=(lf)*log10((sbt_at1_p + sbt_ta1_p)/2), col="red", lwd = 1)
	lines(x=brk, y=(lf)*log10((sbt_at2_p + sbt_ta2_p)/2), col="red", lty="twodash", lwd = 1)

	# AG TC >= 0 sbt_ag sbt_tc green
	lines(x=brk, y=(lf)*log10((sbt_ag1_p + sbt_tc1_p)/2), col="green", lwd = 1)
	lines(x=brk, y=(lf)*log10((sbt_ag2_p + sbt_tc2_p)/2), col="green", lty="twodash", lwd = 1)
	# CA GT >= 0 sbt_ca sbt_gt cyan
	lines(x=brk, y=(lf)*log10((sbt_ca1_p + sbt_gt1_p)/2), col="cyan",lwd = 1)
	lines(x=brk, y=(lf)*log10((sbt_ca2_p + sbt_gt2_p)/2), col="cyan", lty="twodash", lwd = 1)

	# CT GA >= 0 sbt_ct sbt_ga purple
	lines(x=brk, y=(lf)*log10((sbt_ct1_p + sbt_ga1_p)/2), col="purple", lwd = 1)
	lines(x=brk, y=(lf)*log10((sbt_ct2_p + sbt_ga2_p)/2), col="purple", lty="twodash", lwd = 1)

	# CG GC >= 0 sbt_cg sbt_gc orange
	lines(x=brk, y=(lf)*log10((sbt_cg1_p + sbt_gc1_p)/2), col="orange", lwd = 1) 
	lines(x=brk, y=(lf)*log10((sbt_cg2_p + sbt_gc2_p)/2), col="orange", lty="twodash", lwd = 1) 

	axis(1, cex.axis = .75)
	#axis(2, pos=c(0,0), cex.axis = .5, labels = FALSE)
	axis(2, pos=c(0,-10), cex.axis = .5, yaxp = c(-10, 0, 10), las=2)

	par(old.par)
	dev.off()
}


if( (oneortwo==1)  & (type=="pvalue") & ( (graphortable=="groutfpval") | (graphortable=="grout") ))
{
	old.par <- par(mfrow=c(1, 1))
	
#----------------
	prf=0.95
#----------------
	dfrm <- vectorOfTables[[1]]
	nbe  <- as.numeric(dfrm$number_base_expected);
	dfrm <- subset(dfrm, nbe>0);
	chr  <- dfrm$chromosome;
	be   <- dfrm$base_expected;
	alt  <- dfrm$alternative;
	site <- as.numeric(dfrm$site);
	tA   <- as.numeric(dfrm$number_of_As);
	tC   <- as.numeric(dfrm$number_of_Cs);
	tT   <- as.numeric(dfrm$number_of_Ts);
	tG   <- as.numeric(dfrm$number_of_Gs);
	nbe  <- as.numeric(dfrm$number_base_expected);
	xv   <- as.numeric(dfrm$x);
	yv   <- as.numeric(dfrm$y);
	dv = nbe

	print("table 1");

	sbt_ac <- subset(tC/dv, (be=='A' & alt=='A' & tA/dv>prf & xv>xval & yv>yval & tC/dv < range & dv>cover)) # red
	sbt_at <- subset(tT/dv, (be=='A' & alt=='A' & tA/dv>prf & xv>xval & yv>yval & tT/dv < range & dv>cover)) # blue 
	sbt_ag <- subset(tG/dv, (be=='A' & alt=='A' & tA/dv>prf & xv>xval & yv>yval & tG/dv < range & dv>cover)) # green

	sbt_tg <- subset(tG/dv, (be=='T' & alt=='T' & tT/dv>prf & xv>xval & yv>yval & tG/dv < range & dv>cover)) # purple
	sbt_ta <- subset(tA/dv, (be=='T' & alt=='T' & tT/dv>prf & xv>xval & yv>yval & tA/dv < range & dv>cover)) # cyan
	sbt_tc <- subset(tC/dv, (be=='T' & alt=='T' & tT/dv>prf & xv>xval & yv>yval & tC/dv < range & dv>cover)) # orange

	sbt_ca <- subset(tA/dv, (be=='C' & alt=='C' & tC/dv>prf & xv>xval & yv>yval & tA/dv < range & dv>cover)) # red 
	sbt_ct <- subset(tT/dv, (be=='C' & alt=='C' & tC/dv>prf & xv>xval & yv>yval & tT/dv < range & dv>cover)) # blue 
	sbt_cg <- subset(tG/dv, (be=='C' & alt=='C' & tC/dv>prf & xv>xval & yv>yval & tG/dv < range & dv>cover)) # green

	sbt_gt <- subset(tT/dv, (be=='G' & alt=='G' & tG/dv>prf & xv>xval & yv>yval & tT/dv < range & dv>cover)) # purple
	sbt_ga <- subset(tA/dv, (be=='G' & alt=='G' & tG/dv>prf & xv>xval & yv>yval & tA/dv < range & dv>cover)) # cyan 
	sbt_gc <- subset(tC/dv, (be=='G' & alt=='G' & tG/dv>prf & xv>xval & yv>yval & tC/dv < range & dv>cover)) # orange

	sbt_ac.cut = cut(sbt_ac, brk, right=FALSE)
	sbt_ac.freq = table(sbt_ac.cut)
	cumfreq0 = c(0, cumsum(sbt_ac.freq))
	sumfreq0=sum(sbt_ac.freq)
	sbt_ac1 = cumfreq0/sumfreq0
	sbt_ac1_p = 1-cumfreq0/sumfreq0
	
	sbt_at.cut = cut(sbt_at, brk, right=FALSE)
	sbt_at.freq = table(sbt_at.cut)
	cumfreq0 = c(0, cumsum(sbt_at.freq))
	sumfreq0=sum(sbt_at.freq)
	sbt_at1 = cumfreq0/sumfreq0
	sbt_at1_p = 1-cumfreq0/sumfreq0
	
	sbt_ag.cut = cut(sbt_ag, brk, right=FALSE)
	sbt_ag.freq = table(sbt_ag.cut)
	cumfreq0 = c(0, cumsum(sbt_ag.freq))
	sumfreq0=sum(sbt_ag.freq)
	sbt_ag1 = cumfreq0/sumfreq0
	sbt_ag1_p = 1-cumfreq0/sumfreq0

	sbt_tg.cut = cut(sbt_tg, brk, right=FALSE)
	sbt_tg.freq = table(sbt_tg.cut)
	cumfreq0 = c(0, cumsum(sbt_tg.freq))
	sumfreq0=sum(sbt_tg.freq)
	sbt_tg1 = cumfreq0/sumfreq0
	sbt_tg1_p = 1-cumfreq0/sumfreq0

	sbt_ta.cut = cut(sbt_ta, brk, right=FALSE)
	sbt_ta.freq = table(sbt_ta.cut)
	cumfreq0 = c(0, cumsum(sbt_ta.freq))
	sumfreq0=sum(sbt_ta.freq)
	sbt_ta1 = cumfreq0/sumfreq0
	sbt_ta1_p = 1-cumfreq0/sumfreq0

	sbt_tc.cut = cut(sbt_tc, brk, right=FALSE)
	sbt_tc.freq = table(sbt_tc.cut)
	cumfreq0 = c(0, cumsum(sbt_tc.freq))
	sumfreq0=sum(sbt_tc.freq)
	sbt_tc1 = cumfreq0/sumfreq0
	sbt_tc1_p = 1-cumfreq0/sumfreq0

	sbt_ca.cut = cut(sbt_ca, brk, right=FALSE)
	sbt_ca.freq = table(sbt_ca.cut)
	cumfreq0 = c(0, cumsum(sbt_ca.freq))
	sumfreq0=sum(sbt_ca.freq)
	sbt_ca1 = cumfreq0/sumfreq0
	sbt_ca1_p = 1-cumfreq0/sumfreq0

	sbt_ct.cut = cut(sbt_ct, brk, right=FALSE)
	sbt_ct.freq = table(sbt_ct.cut)
	cumfreq0 = c(0, cumsum(sbt_ct.freq))
	sumfreq0=sum(sbt_ct.freq)
	sbt_ct1 = cumfreq0/sumfreq0
	sbt_ct1_p = 1-cumfreq0/sumfreq0

	sbt_cg.cut = cut(sbt_cg, brk, right=FALSE)
	sbt_cg.freq = table(sbt_cg.cut)
	cumfreq0 = c(0, cumsum(sbt_cg.freq))
	sumfreq0=sum(sbt_cg.freq)
	sbt_cg1 = cumfreq0/sumfreq0
	sbt_cg1_p = 1-cumfreq0/sumfreq0

	sbt_gt.cut = cut(sbt_gt, brk, right=FALSE)
	sbt_gt.freq = table(sbt_gt.cut)
	cumfreq0 = c(0, cumsum(sbt_gt.freq))
	sumfreq0=sum(sbt_gt.freq)
	sbt_gt1 = cumfreq0/sumfreq0
	sbt_gt1_p = 1-cumfreq0/sumfreq0

	sbt_ga.cut = cut(sbt_ga, brk, right=FALSE)
	sbt_ga.freq = table(sbt_ga.cut)
	cumfreq0 = c(0, cumsum(sbt_ga.freq))
	sumfreq0=sum(sbt_ga.freq)
	sbt_ga1 = cumfreq0/sumfreq0
	sbt_ga1_p = 1-cumfreq0/sumfreq0

	sbt_gc.cut = cut(sbt_gc, brk, right=FALSE)
	sbt_gc.freq = table(sbt_gc.cut)
	cumfreq0 = c(0, cumsum(sbt_gc.freq))
	sumfreq0=sum(sbt_gc.freq)
	sbt_gc1 = cumfreq0/sumfreq0
	sbt_gc1_p = 1-cumfreq0/sumfreq0


#ymaxi = max((h11$counts+h21$counts)/2, h12$counts, h22$counts, h13$counts, h23$counts, ho11$counts, ho21$counts, ho12$counts, ho22$counts, ho13$counts, ho23$counts)
	#ym = ymaxi + 0.1*ymaxi
	#yl=c(0,ym)
	
	#ymaxi = max(log2((sbt_ac1_p + sbt_tg1_p)/2), log2((sbt_ac1_po + sbt_ta1_po)/2), y=log2((sbt_at1_p + sbt_ta1_p)/2), log2((sbt_at1_po + sbt_ta1_po)/2), log2((sbt_ag1_p + sbt_tc1_p)/2), log2((sbt_ag1_po + sbt_tc1_po)/2), log2((sbt_ca1_p + sbt_gt1_p)/2), log2((sbt_ca1_po + sbt_gt1_po)/2), log2((sbt_ct1_p + sbt_ga1_p)/2), log2((sbt_ct1_po + sbt_ga1_po)/2), log2((sbt_cg1_p + sbt_gc1_p)/2), log2((sbt_cg1_po + sbt_gc1_po)/2))
	
	ymaxi = max(log10((sbt_ac1_p + sbt_tg1_p)/2), y=log10((sbt_at1_p + sbt_ta1_p)/2), log10((sbt_ag1_p + sbt_tc1_p)/2), log10((sbt_ca1_p + sbt_gt1_p)/2), log10((sbt_ct1_p + sbt_ga1_p)/2), log10((sbt_cg1_p + sbt_gc1_p)/2))
	
	ym = ymaxi + 0.1*ymaxi
	yl=c(-3,0)
	subia = paste("A->C & T->G = blue,   A->T & T->A = red   A->G & T->C = green   C->A & G->T = cyan,   C->T & G->A = purple,   C->G & G->C = orange \n", sep="")
	#subia = paste(subi, "i6_sep: A -> C & T -> G = d.blue,  A -> T & T -> A = d.red,   A -> G & T -> C = d.green \n C -> A & G -> T = d.cyan,  C -> T & G -> A = d.purple,  C -> G & G -> C = d.orange \n", sep="");
	
	#subi = paste(subia, qual[i], " < A -> C, T -> G = dot_blue,    A -> T, T -> A = dot_red,    A -> G, T -> C = dot_green \n", qual[i], " < C -> A, G -> T = dot_cyan,    C -> T, G -> A = dot_purple,    C -> G, G -> C = dot_orange", sep="")
	#ttl = "Normalized_Cumulative"
	#ttl  = "p_value, sites size>=100 \ny=log(p-value)"
	#ttl  = paste("Mutation Rate: sites size>=", (cover+1), ",  y=average counts\n", sep="")
	#ttl  = paste("p_value: sites size>=", (cover+1), ",  y=log(p-value), reference>95%\n", sep="")
	ttl  = paste("p_value: sites size>=", (cover+1), ", reference > ", (prf*100), "%,  y=log(p-value)\n",sep="");
	# p value graphs
	# AC TG >= 0 sbt_ac sbt_tg blue
	plot(x=brk, y=(lf)*log10((sbt_ac1_p + sbt_tg1_p)/2), col="blue", lwd = 1, type="l", ylim=yl, xlim=xl, yaxt='n', xaxt="n", sub = subia, cex.sub = 0.70, font.sub = 1, xlab=NA, main=paste(ttl,"file: ", tblnms[1], sep=""), cex.main=0.70, ylab=tblnms[1])
	
	# AT TA >= 0 sbt_at sbt_ta red
	lines(x=brk, y=(lf)*log10((sbt_at1_p + sbt_ta1_p)/2), col="red", lwd = 1)

	# AG TC >= 0 sbt_ag sbt_tc green
	lines(x=brk, y=(lf)*log10((sbt_ag1_p + sbt_tc1_p)/2), col="green", lwd = 1)

	# CA GT >= 0 sbt_ca sbt_gt cyan
	lines(x=brk, y=(lf)*log10((sbt_ca1_p + sbt_gt1_p)/2), col="cyan",lwd = 1)

	# CT GA >= 0 sbt_ct sbt_ga purple
	lines(x=brk, y=(lf)*log10((sbt_ct1_p + sbt_ga1_p)/2), col="purple", lwd = 1)

	# CG GC >= 0 sbt_cg sbt_gc orange
	lines(x=brk, y=(lf)*log10((sbt_cg1_p + sbt_gc1_p)/2), col="orange", lwd = 1) 
	
	axis(1, cex.axis = .75)
	#axis(2, pos=c(0,0), cex.axis = .5, labels = FALSE)
	axis(2, pos=c(0,-10), cex.axis = .5, yaxp = c(-10, 0, 10), las=2)
	par(old.par)
	dev.off()
}


if( (oneortwo==2) & (type=="mutrate") & ( (graphortable=="groutfpval") | (graphortable=="grout") ))
{
	old.par <- par(mfrow=c(1, 1))
	
#----------------
	prf=0.95
#----------------
	dfrm <- vectorOfTables[[1]]
	nbe  <- as.numeric(dfrm$number_base_expected);
	dfrm <- subset(dfrm, nbe>0);
	chr  <- dfrm$chromosome;
	be   <- dfrm$base_expected;
	alt  <- dfrm$alternative;
	site <- as.numeric(dfrm$site);
	tA   <- as.numeric(dfrm$number_of_As);
	tC   <- as.numeric(dfrm$number_of_Cs);
	tT   <- as.numeric(dfrm$number_of_Ts);
	tG   <- as.numeric(dfrm$number_of_Gs);
	nbe  <- as.numeric(dfrm$number_base_expected);
	xv   <- as.numeric(dfrm$x);
	yv   <- as.numeric(dfrm$y);
	dv = nbe

	print("table 1");

	sbt_ac <- subset(tC/dv, (be=='A' & alt=='A' & tA/dv>prf & xv>xval & yv>yval & tC/dv < 0.005 & dv>cover)) # red
	sbt_at <- subset(tT/dv, (be=='A' & alt=='A' & tA/dv>prf & xv>xval & yv>yval & tT/dv < 0.005 & dv>cover)) # blue 
	sbt_ag <- subset(tG/dv, (be=='A' & alt=='A' & tA/dv>prf & xv>xval & yv>yval & tG/dv < 0.005 & dv>cover)) # green

	sbt_tg <- subset(tG/dv, (be=='T' & alt=='T' & tT/dv>prf & xv>xval & yv>yval & tG/dv < 0.005 & dv>cover)) # purple
	sbt_ta <- subset(tA/dv, (be=='T' & alt=='T' & tT/dv>prf & xv>xval & yv>yval & tA/dv < 0.005 & dv>cover)) # cyan
	sbt_tc <- subset(tC/dv, (be=='T' & alt=='T' & tT/dv>prf & xv>xval & yv>yval & tC/dv < 0.005 & dv>cover)) # orange

	sbt_ca <- subset(tA/dv, (be=='C' & alt=='C' & tC/dv>prf & xv>xval & yv>yval & tA/dv < 0.005 & dv>cover)) # red 
	sbt_ct <- subset(tT/dv, (be=='C' & alt=='C' & tC/dv>prf & xv>xval & yv>yval & tT/dv < 0.005 & dv>cover)) # blue 
	sbt_cg <- subset(tG/dv, (be=='C' & alt=='C' & tC/dv>prf & xv>xval & yv>yval & tG/dv < 0.005 & dv>cover)) # green

	sbt_gt <- subset(tT/dv, (be=='G' & alt=='G' & tG/dv>prf & xv>xval & yv>yval & tT/dv < 0.005 & dv>cover)) # purple
	sbt_ga <- subset(tA/dv, (be=='G' & alt=='G' & tG/dv>prf & xv>xval & yv>yval & tA/dv < 0.005 & dv>cover)) # cyan 
	sbt_gc <- subset(tC/dv, (be=='G' & alt=='G' & tG/dv>prf & xv>xval & yv>yval & tC/dv < 0.005 & dv>cover)) # orange
	

	h11 <- hist(sbt_ac, breaks=brk, plot=FALSE)		# ac-tg
	h12 <- hist(sbt_at, breaks=brk, plot=FALSE)		# at-ta
	h13 <- hist(sbt_ag, breaks=brk, plot=FALSE)		# ag-tc
	
	h21 <- hist(sbt_tg, breaks=brk, plot=FALSE)		# tg-ac
	h22 <- hist(sbt_ta, breaks=brk, plot=FALSE)		# ta-at
	h23 <- hist(sbt_tc, breaks=brk, plot=FALSE)		# tc-ag
	
	h31 <- hist(sbt_ca, breaks=brk, plot=FALSE)		# ca-gt
	h32 <- hist(sbt_ct, breaks=brk, plot=FALSE)		# ct-ga
	h33 <- hist(sbt_cg, breaks=brk, plot=FALSE)		# cg-gc

	h41 <- hist(sbt_gt, breaks=brk, plot=FALSE)		# gt-ac
	h42 <- hist(sbt_ga, breaks=brk, plot=FALSE)		# ga-ct
	h43 <- hist(sbt_gc, breaks=brk, plot=FALSE)		# gc-cg

	dfrm <- vectorOfTables[[2]]
	nbe  <- as.numeric(dfrm$number_base_expected);
	dfrm <- subset(dfrm, nbe>0);
	chr  <- dfrm$chromosome;
	be   <- dfrm$base_expected;
	alt  <- dfrm$alternative;
	site <- as.numeric(dfrm$site);
	tA   <- as.numeric(dfrm$number_of_As);
	tC   <- as.numeric(dfrm$number_of_Cs);
	tT   <- as.numeric(dfrm$number_of_Ts);
	tG   <- as.numeric(dfrm$number_of_Gs);
	nbe  <- as.numeric(dfrm$number_base_expected);
	xv   <- as.numeric(dfrm$x);
	yv   <- as.numeric(dfrm$y);
	dv = nbe

	print("table 2");

	sbt_ac <- subset(tC/dv, (be=='A' & alt=='A' & tA/dv>prf & xv>xval & yv>yval & tC/dv < range & dv>cover)) # red
	sbt_at <- subset(tT/dv, (be=='A' & alt=='A' & tA/dv>prf & xv>xval & yv>yval & tT/dv < range & dv>cover)) # blue 
	sbt_ag <- subset(tG/dv, (be=='A' & alt=='A' & tA/dv>prf & xv>xval & yv>yval & tG/dv < range & dv>cover)) # green

	sbt_tg <- subset(tG/dv, (be=='T' & alt=='T' & tT/dv>prf & xv>xval & yv>yval & tG/dv < range & dv>cover)) # purple
	sbt_ta <- subset(tA/dv, (be=='T' & alt=='T' & tT/dv>prf & xv>xval & yv>yval & tA/dv < range & dv>cover)) # cyan
	sbt_tc <- subset(tC/dv, (be=='T' & alt=='T' & tT/dv>prf & xv>xval & yv>yval & tC/dv < range & dv>cover)) # orange

	sbt_ca <- subset(tA/dv, (be=='C' & alt=='C' & tC/dv>prf & xv>xval & yv>yval & tA/dv < range & dv>cover)) # red 
	sbt_ct <- subset(tT/dv, (be=='C' & alt=='C' & tC/dv>prf & xv>xval & yv>yval & tT/dv < range & dv>cover)) # blue 
	sbt_cg <- subset(tG/dv, (be=='C' & alt=='C' & tC/dv>prf & xv>xval & yv>yval & tG/dv < range & dv>cover)) # green

	sbt_gt <- subset(tT/dv, (be=='G' & alt=='G' & tG/dv>prf & xv>xval & yv>yval & tT/dv < range & dv>cover)) # purple
	sbt_ga <- subset(tA/dv, (be=='G' & alt=='G' & tG/dv>prf & xv>xval & yv>yval & tA/dv < range & dv>cover)) # cyan 
	sbt_gc <- subset(tC/dv, (be=='G' & alt=='G' & tG/dv>prf & xv>xval & yv>yval & tC/dv < range & dv>cover)) # orange

	ho11 <- hist(sbt_ac, breaks=brk, plot=FALSE)		# ac-tg
	ho12 <- hist(sbt_at, breaks=brk, plot=FALSE)		# at-ta
	ho13 <- hist(sbt_ag, breaks=brk, plot=FALSE)		# ag-tc
	
	ho21 <- hist(sbt_tg, breaks=brk, plot=FALSE)		# tg-ac
	ho22 <- hist(sbt_ta, breaks=brk, plot=FALSE)		# ta-at
	ho23 <- hist(sbt_tc, breaks=brk, plot=FALSE)		# tc-ag
	
	ho31 <- hist(sbt_ca, breaks=brk, plot=FALSE)		# ca-gt
	ho32 <- hist(sbt_ct, breaks=brk, plot=FALSE)		# ct-ga
	ho33 <- hist(sbt_cg, breaks=brk, plot=FALSE)		# cg-gc

	ho41 <- hist(sbt_gt, breaks=brk, plot=FALSE)		# gt-ac
	ho42 <- hist(sbt_ga, breaks=brk, plot=FALSE)		# ga-ct
	ho43 <- hist(sbt_gc, breaks=brk, plot=FALSE)		# gc-cg


	ymaxi = max(h11$counts, h21$counts, h12$counts, h22$counts, h13$counts, h23$counts, h31$counts, h41$counts, h32$counts, h42$counts, h33$counts, h43$counts, ho11$counts, ho21$counts, ho12$counts, ho22$counts, ho13$counts, ho23$counts, ho31$counts, ho41$counts, ho32$counts, ho42$counts, ho33$counts, ho43$counts)
	ym = ymaxi + 0.1*ymaxi
	yl=c(0,ym)
	ttl  = paste("Mutation Rate: sites size>=", (cover+1), ",  y=combined average counts\n", sep="")
	subi = paste("file1:  A->C & T->G = blue,   A->T & T->A = red   A->G & T->C = green   C->A & G->T = cyan,   C->T & G->A = purple,   C->G & G->C = orange \n", sep="")
	subia = paste(subi, "file2:  A->C & T->G = dashblue,   A->T & T->A = dashred   A->G & T->C = dashgreen   C->A & G->T = dashcyan,   C->T & G->A = dashpurple,   C->G & G->C = dashorange \n", sep="")

	#subi = paste("A -> C = blue, T -> G = dot_cyan\nA -> T = red, T -> A = dot_purple\nA -> G = green, T -> C = dot_orange \n ", ym, sep="")
	
	plot(x=(h11$mids+h21$mids)/2, y=(h11$counts+h21$counts)/2, lwd = 1, type="l", yaxt="n", xlim=xl, ylim = yl, xaxt="n", col="blue", sub = subia, cex.sub = 0.75, font.sub = 1, xlab=NA, cex.main = 0.75, main=paste( ttl, "file1: ", tblnms[1], "\nfile2: ", tblnms[2], sep=""), ylab=paste(tblnms[1], "  ", tblnms[2], sep="") )
	
	lines((x=h12$mids+h22$mids)/2, y=(h12$counts+h22$counts)/2, lwd = 1, col="red")
	lines((x=h13$mids+h23$mids)/2, y=(h13$counts+h23$counts)/2, lwd = 1, col="green")
	lines((x=h31$mids+h41$mids)/2, y=(h31$counts+h41$counts)/2, lwd = 1, col="cyan")
	lines((x=h32$mids+h42$mids)/2, y=(h32$counts+h42$counts)/2, lwd = 1, col="purple")
	lines((x=h33$mids+h43$mids)/2, y=(h33$counts+h43$counts)/2, lwd = 1, col="orange")

	lines((x=ho11$mids+ho21$mids)/2, y=(ho11$counts+ho21$counts)/2, lwd = 1,  lty="twodash", col="blue")	
	lines((x=ho12$mids+ho22$mids)/2, y=(ho12$counts+ho22$counts)/2, lwd = 1,  lty="twodash", col="red")
	lines((x=ho13$mids+ho23$mids)/2, y=(ho13$counts+ho23$counts)/2, lwd = 1,  lty="twodash", col="green")
	lines((x=ho31$mids+ho41$mids)/2, y=(ho31$counts+ho41$counts)/2, lwd = 1,  lty="twodash", col="cyan")
	lines((x=ho32$mids+ho42$mids)/2, y=(ho32$counts+ho42$counts)/2, lwd = 1,  lty="twodash", col="purple")
	lines((x=ho33$mids+ho43$mids)/2, y=(ho33$counts+ho43$counts)/2, lwd = 1,  lty="twodash", col="orange")


	axis(1, cex.axis = .75)
	axis(2, pos=c(0,0), cex.axis = .75)

	par(old.par)
	dev.off()
}


if( (oneortwo==1)  & (type=="mutrate")  & ( (graphortable=="groutfpval") | (graphortable=="grout") ) )
{
	old.par <- par(mfrow=c(1, 1))

	
#----------------
	prf=0.95
#----------------	
	dfrm0 <- vectorOfTables[[1]]
	nbe  <- as.numeric(dfrm$number_base_expected);
	dfrm <- subset(dfrm, nbe>0);
	chr  <- dfrm$chromosome;
	be   <- dfrm$base_expected;
	alt  <- dfrm$alternative;
	site <- as.numeric(dfrm$site);
	tA   <- as.numeric(dfrm$number_of_As);
	tC   <- as.numeric(dfrm$number_of_Cs);
	tT   <- as.numeric(dfrm$number_of_Ts);
	tG   <- as.numeric(dfrm$number_of_Gs);
	nbe  <- as.numeric(dfrm$number_base_expected);
	xv   <- as.numeric(dfrm$x);
	yv   <- as.numeric(dfrm$y);
	dv = nbe

	print(1)

	sbt_ac <- subset(tC/dv, (be=='A' & alt=='A' & tA/dv>prf & xv>xval & yv>yval & tC/dv < range & dv>cover)) # red
	sbt_at <- subset(tT/dv, (be=='A' & alt=='A' & tA/dv>prf & xv>xval & yv>yval & tT/dv < range & dv>cover)) # blue 
	sbt_ag <- subset(tG/dv, (be=='A' & alt=='A' & tA/dv>prf & xv>xval & yv>yval & tG/dv < range & dv>cover)) # green

	sbt_tg <- subset(tG/dv, (be=='T' & alt=='T' & tT/dv>prf & xv>xval & yv>yval & tG/dv < range & dv>cover)) # purple
	sbt_ta <- subset(tA/dv, (be=='T' & alt=='T' & tT/dv>prf & xv>xval & yv>yval & tA/dv < range & dv>cover)) # cyan
	sbt_tc <- subset(tC/dv, (be=='T' & alt=='T' & tT/dv>prf & xv>xval & yv>yval & tC/dv < range & dv>cover)) # orange

	sbt_ca <- subset(tA/dv, (be=='C' & alt=='C' & tC/dv>prf & xv>xval & yv>yval & tA/dv < range & dv>cover)) # red 
	sbt_ct <- subset(tT/dv, (be=='C' & alt=='C' & tC/dv>prf & xv>xval & yv>yval & tT/dv < range & dv>cover)) # blue 
	sbt_cg <- subset(tG/dv, (be=='C' & alt=='C' & tC/dv>prf & xv>xval & yv>yval & tG/dv < range & dv>cover)) # green

	sbt_gt <- subset(tT/dv, (be=='G' & alt=='G' & tG/dv>prf & xv>xval & yv>yval & tT/dv < range & dv>cover)) # purple
	sbt_ga <- subset(tA/dv, (be=='G' & alt=='G' & tG/dv>prf & xv>xval & yv>yval & tA/dv < range & dv>cover)) # cyan 
	sbt_gc <- subset(tC/dv, (be=='G' & alt=='G' & tG/dv>prf & xv>xval & yv>yval & tC/dv < range & dv>cover)) # orange

	h11 <- hist(sbt_ac, breaks=brk, plot=FALSE)		# ac-tg
	h12 <- hist(sbt_at, breaks=brk, plot=FALSE)		# at-ta
	h13 <- hist(sbt_ag, breaks=brk, plot=FALSE)		# ag-tc
	
	h21 <- hist(sbt_tg, breaks=brk, plot=FALSE)		# tg-ac
	h22 <- hist(sbt_ta, breaks=brk, plot=FALSE)		# ta-at
	h23 <- hist(sbt_tc, breaks=brk, plot=FALSE)		# tc-ag
	
	h31 <- hist(sbt_ca, breaks=brk, plot=FALSE)		# ca-gt
	h32 <- hist(sbt_ct, breaks=brk, plot=FALSE)		# ct-ga
	h33 <- hist(sbt_cg, breaks=brk, plot=FALSE)		# cg-gc

	h41 <- hist(sbt_gt, breaks=brk, plot=FALSE)		# gt-ac
	h42 <- hist(sbt_ga, breaks=brk, plot=FALSE)		# ga-ct
	h43 <- hist(sbt_gc, breaks=brk, plot=FALSE)		# gc-cg	

	ymaxi = max(h11$counts, h21$counts, h12$counts, h22$counts, h13$counts, h23$counts, h31$counts, h41$counts, h32$counts, h42$counts, h33$counts, h43$counts)
	ym = ymaxi + 0.1*ymaxi
	yl=c(0,ym)
	ttl  = paste("Mutation Rate: sites size>=", (cover+1), ", y=counts\n", sep="")
	subi = paste("file:  A->C & T->G = blue,   A->T & T->A = red   A->G & T->C = green   C->A & G->T = cyan,   C->T & G->A = purple,   C->G & G->C = orange \n", sep="")
	#subia = paste(subi, "file2:  A->C & T->G = dashblue,   A->T & T->A = dashred   A->G & T->C = dashgreen   C->A & G->T = dashcyan,   C->T & G->A = dashpurple,   C->G & G->C = dashorange \n", sep="")

	#subi = paste("A -> C = blue, T -> G = dot_cyan\nA -> T = red, T -> A = dot_purple\nA -> G = green, T -> C = dot_orange \n ", ym, sep="")
	
	plot(x=(h11$mids+h21$mids)/2, y=(h11$counts+h21$counts)/2, lwd = 1, type="l", yaxt="n", xlim=xl, ylim = yl, xaxt="n", col="blue", sub = subi, cex.sub = 0.75, font.sub = 1, xlab=NA,  cex.main = 0.75, main=paste(ttl, "file: ", tblnms[1], sep=""), ylab=tblnms[1])
	lines((x=h12$mids+h22$mids)/2, y=(h12$counts+h22$counts)/2, lwd = 1, col="red")
	lines((x=h13$mids+h23$mids)/2, y=(h13$counts+h23$counts)/2, lwd = 1, col="green")
	lines((x=h31$mids+h41$mids)/2, y=(h31$counts+h41$counts)/2, lwd = 1, col="cyan")
	lines((x=h32$mids+h42$mids)/2, y=(h32$counts+h42$counts)/2, lwd = 1, col="purple")
	lines((x=h33$mids+h43$mids)/2, y=(h33$counts+h43$counts)/2, lwd = 1, col="orange")

	axis(1, cex.axis = .75)
	axis(2, pos=c(0,0), cex.axis = .75)

	par(old.par)
	dev.off()

}

if((oneortwo==1)  &  ((graphortable=="groutfpval") | (graphortable=="fpval")) )
{
	options(scipen=999)
	#source = "C:\\Users\\m026918\\Desktop\\MosNew\\MosaicsNewDataCorrected\\Mosaics_i6\\Tables\\ge_0\\MosaicNoiseTable_BQ_gt_00.tsv"
	#destination = "C:\\Users\\m026918\\Desktop\\dest\\table\\"

	#dfrm <- read.table(source,sep="\t",na.strings = "NA",header = TRUE)
	
	yl=c(0,10000)
	#	brk <- seq(0,0.005, 6.25e-5)
	#	brk <- seq(0,1.0, 10e-5)
	right_points <- seq(0,1.0, 10e-5)
	brk <- seq(-0.0001,1.0, 10e-5)
	xl=c(-0.0,1.00)
##	xval=0
##	yval=0
	j=0

	##chr  <- dfrm$chromosome;
	##site <- dfrm$site;
	#chrsite=paste(chr, " ", site, sep="")
	# dfrm=subset(dfrm, !chrsite%in%excl) # excludes the sites of interest excl is a vector generated from the primers file(element : chr, space, site of interest. exemple : "2 2345673" means chromosome: 2 site: 2345673
	dfrm <- vectorOfTables[[1]]
	nbe  <- as.numeric(dfrm$number_base_expected);
	dfrm <- subset(dfrm, nbe>0);
	chr  <- dfrm$chromosome;
	be   <- dfrm$base_expected;
	alt  <- dfrm$alternative;
	site <- as.numeric(dfrm$site);
	tA   <- as.numeric(dfrm$number_of_As);
	tC   <- as.numeric(dfrm$number_of_Cs);
	tT   <- as.numeric(dfrm$number_of_Ts);
	tG   <- as.numeric(dfrm$number_of_Gs);
	nbe  <- as.numeric(dfrm$number_base_expected);
	xv   <- as.numeric(dfrm$x);
	yv   <- as.numeric(dfrm$y);
	dv = nbe;
	crlimit=1.5

#	xl=c(-0.0001,0.005)
	yl=c(0,10000)
#	brk <- seq(0,0.005, 6.25e-5)
#	brk <- seq(0,1.0, 10e-5)
	right_points <- seq(0,1.0, 10e-5)
	brk <- seq(-0.0001,1.0, 10e-5)
	#setwd("C:\\Users\\m026918\\Desktop");
	#prf=0.95 # percentage of reference
	xl=c(-0.0,1.00)
##	xval=0
##	yval=0
	j=0
	i=1

sbt_ac <- subset(tC/dv, (be=='A' & alt=='A' & tA/dv>prf & xv>xval & yv>yval & tC/dv < range & dv>cover)) # red
sbt_at <- subset(tT/dv, (be=='A' & alt=='A' & tA/dv>prf & xv>xval & yv>yval & tT/dv < range & dv>cover)) # blue 
sbt_ag <- subset(tG/dv, (be=='A' & alt=='A' & tA/dv>prf & xv>xval & yv>yval & tG/dv < range & dv>cover)) # green

sbt_tg <- subset(tG/dv, (be=='T' & alt=='T' & tT/dv>prf & xv>xval & yv>yval & tG/dv < range & dv>cover)) # purple
sbt_ta <- subset(tA/dv, (be=='T' & alt=='T' & tT/dv>prf & xv>xval & yv>yval & tA/dv < range & dv>cover)) # cyan
sbt_tc <- subset(tC/dv, (be=='T' & alt=='T' & tT/dv>prf & xv>xval & yv>yval & tC/dv < range & dv>cover)) # orange

sbt_ca <- subset(tA/dv, (be=='C' & alt=='C' & tC/dv>prf & xv>xval & yv>yval & tA/dv < range & dv>cover)) # red 
sbt_ct <- subset(tT/dv, (be=='C' & alt=='C' & tC/dv>prf & xv>xval & yv>yval & tT/dv < range & dv>cover)) # blue 
sbt_cg <- subset(tG/dv, (be=='C' & alt=='C' & tC/dv>prf & xv>xval & yv>yval & tG/dv < range & dv>cover)) # green

sbt_gt <- subset(tT/dv, (be=='G' & alt=='G' & tG/dv>prf & xv>xval & yv>yval & tT/dv < range & dv>cover)) # purple
sbt_ga <- subset(tA/dv, (be=='G' & alt=='G' & tG/dv>prf & xv>xval & yv>yval & tA/dv < range & dv>cover)) # cyan 
sbt_gc <- subset(tC/dv, (be=='G' & alt=='G' & tG/dv>prf & xv>xval & yv>yval & tC/dv < range & dv>cover)) # orange

	sbt_ac.cut = cut(sbt_ac, brk)
#	, right=FALSE) 
	sbt_ac.freq = table(sbt_ac.cut)
	
	sbt_tg.cut = cut(sbt_tg, brk)
#	, right=FALSE)
	sbt_tg.freq = table(sbt_tg.cut)

	sbt_ac_tg.freq = sbt_ac.freq + sbt_tg.freq	
	
	sumfreq0=sum(sbt_ac_tg.freq)
	cumfreq0 = c(0, cumsum(sbt_ac_tg.freq))
	numerator=sumfreq0-cumfreq0
	
	sbt_ac_tg_p = numerator/sumfreq0
	sbt_ac_tg_p <- sbt_ac_tg_p[1:length(sbt_ac_tg_p)-1] #
	
	ac_tg_frq_pval_table1 <- cbind(right_points, sbt_ac_tg.freq, sbt_ac_tg_p)
#	---------------------------------------------------------------------------------------------------	#	

#	---------------------------------------------------------------------------------------------------	#	
	sbt_at.cut = cut(sbt_at, brk)
#	, right=FALSE) 
	sbt_at.freq = table(sbt_at.cut)
	
	sbt_ta.cut = cut(sbt_ta, brk)
#	, right=FALSE)
	sbt_ta.freq = table(sbt_ta.cut)
	
	sbt_at_ta.freq = sbt_at.freq + sbt_ta.freq	

	sumfreq0=sum(sbt_at_ta.freq)
	cumfreq0 = c(0, cumsum(sbt_at_ta.freq))
	numerator=sumfreq0-cumfreq0
	
	sbt_at_ta_p = numerator/sumfreq0
	sbt_at_ta_p <- sbt_at_ta_p[1:length(sbt_at_ta_p)-1] #

	at_ta_frq_pval_table1 <- cbind(right_points, sbt_at_ta.freq, sbt_at_ta_p)
#	---------------------------------------------------------------------------------------------------	#	

#	---------------------------------------------------------------------------------------------------	#	
	sbt_ag.cut = cut(sbt_ag, brk)
#	, right=FALSE)
	sbt_ag.freq = table(sbt_ag.cut)

	sbt_tc.cut = cut(sbt_tc, brk)
#	, right=FALSE)
	sbt_tc.freq = table(sbt_tc.cut)

	sbt_ag_tc.freq = sbt_ag.freq + sbt_tc.freq	

	sumfreq0=sum(sbt_ag_tc.freq)
	cumfreq0 = c(0, cumsum(sbt_ag_tc.freq))
	numerator=sumfreq0-cumfreq0
	
	sbt_ag_tc_p = numerator/sumfreq0
	sbt_ag_tc_p <- sbt_ag_tc_p[1:length(sbt_ag_tc_p)-1] #

	ag_tc_frq_pval_table1 <- cbind(right_points, sbt_ag_tc.freq, sbt_ag_tc_p)
#	---------------------------------------------------------------------------------------------------	#	

#	---------------------------------------------------------------------------------------------------	#	
	sbt_ca.cut = cut(sbt_ca, brk)
#	, right=FALSE) 
	sbt_ca.freq = table(sbt_ca.cut)
	
	sbt_gt.cut = cut(sbt_gt, brk)
#	, right=FALSE)
	sbt_gt.freq = table(sbt_gt.cut)

	sbt_ca_gt.freq = sbt_ca.freq + sbt_gt.freq
	
	sumfreq0=sum(sbt_ca_gt.freq)
	cumfreq0 = c(0, cumsum(sbt_ca_gt.freq))
	numerator=sumfreq0-cumfreq0
	
	sbt_ca_gt_p = numerator/sumfreq0
	sbt_ca_gt_p <- sbt_ca_gt_p[1:length(sbt_ca_gt_p)-1] #

	ca_gt_frq_pval_table1 <- cbind(right_points, sbt_ca_gt.freq, sbt_ca_gt_p)
#	---------------------------------------------------------------------------------------------------	#	

#	---------------------------------------------------------------------------------------------------	#	
	sbt_ct.cut = cut(sbt_ct, brk)
#	, right=FALSE) 
	sbt_ct.freq = table(sbt_ct.cut)
	
	sbt_ga.cut = cut(sbt_ga, brk)
#	, right=FALSE)
	sbt_ga.freq = table(sbt_ga.cut)

	sbt_ct_ga.freq = sbt_ct.freq + sbt_ga.freq	

	sumfreq0=sum(sbt_ct_ga.freq)
	cumfreq0 = c(0, cumsum(sbt_ct_ga.freq))
	numerator=sumfreq0-cumfreq0
	
	sbt_ct_ga_p = numerator/sumfreq0
	sbt_ct_ga_p <- sbt_ct_ga_p[1:length(sbt_ct_ga_p)-1] #

	ct_ga_frq_pval_table1 <- cbind(right_points, sbt_ct_ga.freq, sbt_ct_ga_p)
#	---------------------------------------------------------------------------------------------------	#	

#	---------------------------------------------------------------------------------------------------	#	
	sbt_cg.cut = cut(sbt_cg, brk)
#	, right=FALSE) 
	sbt_cg.freq = table(sbt_cg.cut)
	
	sbt_gc.cut = cut(sbt_gc, brk)
#	, right=FALSE)
	sbt_gc.freq = table(sbt_gc.cut)

	sbt_cg_gc.freq = sbt_cg.freq + sbt_gc.freq	
	
	sumfreq0=sum(sbt_cg_gc.freq)
	cumfreq0 = c(0, cumsum(sbt_cg_gc.freq))
	numerator=sumfreq0-cumfreq0

	sbt_cg_gc_p = numerator/sumfreq0
	sbt_cg_gc_p <- sbt_cg_gc_p[1:length(sbt_cg_gc_p)-1] #

	cg_gc_frq_pval_table1 <- cbind(right_points, sbt_cg_gc.freq, sbt_cg_gc_p)

	frq_pval_table <- cbind(right_points, sbt_ac_tg.freq, sbt_ac_tg_p, sbt_at_ta.freq, sbt_at_ta_p, sbt_ag_tc.freq, sbt_ag_tc_p, sbt_ca_gt.freq, sbt_ca_gt_p, sbt_ct_ga.freq, sbt_ct_ga_p, sbt_cg_gc.freq, sbt_cg_gc_p)
	write.table(frq_pval_table, file= paste(destination, ".fpvtb.tsv", sep=""), sep="\t", row.names=FALSE)

}


if((oneortwo==2)  &  ((graphortable=="groutfpval") | (graphortable=="fpval")) )
{
	options(scipen=999)
	#source = "C:\\Users\\m026918\\Desktop\\MosNew\\MosaicsNewDataCorrected\\Mosaics_i6\\Tables\\ge_0\\MosaicNoiseTable_BQ_gt_00.tsv"
	#destination = "C:\\Users\\m026918\\Desktop\\dest\\table\\"

	#dfrm <- read.table(source,sep="\t",na.strings = "NA",header = TRUE)

	yl=c(0,10000)
	#	brk <- seq(0,0.005, 6.25e-5)
	#	brk <- seq(0,1.0, 10e-5)
	right_points <- seq(0,1.0, 10e-5)
	brk <- seq(-0.0001,1.0, 10e-5)
	xl=c(-0.0,1.00)
##	xval=0
##	yval=0
	j=0

	#chr  <- dfrm$chromosome;
	#site <- dfrm$site;
	#chrsite=paste(chr, " ", site, sep="")
	# dfrm=subset(dfrm, !chrsite%in%excl) # excludes the sites of interest excl is a vector generated from the primers file(element : chr, space, site of interest. exemple : "2 2345673" means chromosome: 2 site: 2345673
	dfrm <- vectorOfTables[[1]]
	nbe  <- as.numeric(dfrm$number_base_expected);
	dfrm <- subset(dfrm, nbe>0);
	chr  <- dfrm$chromosome;
	be   <- dfrm$base_expected;
	alt  <- dfrm$alternative;
	site <- as.numeric(dfrm$site);
	tA   <- as.numeric(dfrm$number_of_As);
	tC   <- as.numeric(dfrm$number_of_Cs);
	tT   <- as.numeric(dfrm$number_of_Ts);
	tG   <- as.numeric(dfrm$number_of_Gs);
	nbe  <- as.numeric(dfrm$number_base_expected);
	xv   <- as.numeric(dfrm$x);
	yv   <- as.numeric(dfrm$y);
	dv = nbe;
	crlimit=1.5

#	xl=c(-0.0001,0.005)
	yl=c(0,10000)
#	brk <- seq(0,0.005, 6.25e-5)
#	brk <- seq(0,1.0, 10e-5)
	right_points <- seq(0,1.0, 10e-5)
	brk <- seq(-0.0001,1.0, 10e-5)
	#setwd("C:\\Users\\m026918\\Desktop");
	prf=0.95 # percentage of reference
	xl=c(-0.0,1.00)
##	xval=0
##	yval=0
	j=0
	i=1

sbt_ac <- subset(tC/dv, (be=='A' & alt=='A' & tA/dv>prf & xv>xval & yv>yval & tC/dv < range & dv>cover)) # red
sbt_at <- subset(tT/dv, (be=='A' & alt=='A' & tA/dv>prf & xv>xval & yv>yval & tT/dv < range & dv>cover)) # blue 
sbt_ag <- subset(tG/dv, (be=='A' & alt=='A' & tA/dv>prf & xv>xval & yv>yval & tG/dv < range & dv>cover)) # green

sbt_tg <- subset(tG/dv, (be=='T' & alt=='T' & tT/dv>prf & xv>xval & yv>yval & tG/dv < range & dv>cover)) # purple
sbt_ta <- subset(tA/dv, (be=='T' & alt=='T' & tT/dv>prf & xv>xval & yv>yval & tA/dv < range & dv>cover)) # cyan
sbt_tc <- subset(tC/dv, (be=='T' & alt=='T' & tT/dv>prf & xv>xval & yv>yval & tC/dv < range & dv>cover)) # orange

sbt_ca <- subset(tA/dv, (be=='C' & alt=='C' & tC/dv>prf & xv>xval & yv>yval & tA/dv < range & dv>cover)) # red 
sbt_ct <- subset(tT/dv, (be=='C' & alt=='C' & tC/dv>prf & xv>xval & yv>yval & tT/dv < range & dv>cover)) # blue 
sbt_cg <- subset(tG/dv, (be=='C' & alt=='C' & tC/dv>prf & xv>xval & yv>yval & tG/dv < range & dv>cover)) # green

sbt_gt <- subset(tT/dv, (be=='G' & alt=='G' & tG/dv>prf & xv>xval & yv>yval & tT/dv < range & dv>cover)) # purple
sbt_ga <- subset(tA/dv, (be=='G' & alt=='G' & tG/dv>prf & xv>xval & yv>yval & tA/dv < range & dv>cover)) # cyan 
sbt_gc <- subset(tC/dv, (be=='G' & alt=='G' & tG/dv>prf & xv>xval & yv>yval & tC/dv < range & dv>cover)) # orange

	sbt_ac.cut = cut(sbt_ac, brk)
#	, right=FALSE) 
	sbt_ac.freq = table(sbt_ac.cut)
	
	sbt_tg.cut = cut(sbt_tg, brk)
#	, right=FALSE)
	sbt_tg.freq = table(sbt_tg.cut)

	sbt_ac_tg.freq = sbt_ac.freq + sbt_tg.freq	
	
	sumfreq0=sum(sbt_ac_tg.freq)
	cumfreq0 = c(0, cumsum(sbt_ac_tg.freq))
	numerator=sumfreq0-cumfreq0
	
	sbt_ac_tg_p = numerator/sumfreq0
	sbt_ac_tg_p <- sbt_ac_tg_p[1:length(sbt_ac_tg_p)-1] #
	
	ac_tg_frq_pval_table1 <- cbind(right_points, sbt_ac_tg.freq, sbt_ac_tg_p)
#	---------------------------------------------------------------------------------------------------	#	

#	---------------------------------------------------------------------------------------------------	#	
	sbt_at.cut = cut(sbt_at, brk)
#	, right=FALSE) 
	sbt_at.freq = table(sbt_at.cut)
	
	sbt_ta.cut = cut(sbt_ta, brk)
#	, right=FALSE)
	sbt_ta.freq = table(sbt_ta.cut)
	
	sbt_at_ta.freq = sbt_at.freq + sbt_ta.freq	

	sumfreq0=sum(sbt_at_ta.freq)
	cumfreq0 = c(0, cumsum(sbt_at_ta.freq))
	numerator=sumfreq0-cumfreq0
	
	sbt_at_ta_p = numerator/sumfreq0
	sbt_at_ta_p <- sbt_at_ta_p[1:length(sbt_at_ta_p)-1] #

	at_ta_frq_pval_table1 <- cbind(right_points, sbt_at_ta.freq, sbt_at_ta_p)
#	---------------------------------------------------------------------------------------------------	#	

#	---------------------------------------------------------------------------------------------------	#	
	sbt_ag.cut = cut(sbt_ag, brk)
#	, right=FALSE)
	sbt_ag.freq = table(sbt_ag.cut)

	sbt_tc.cut = cut(sbt_tc, brk)
#	, right=FALSE)
	sbt_tc.freq = table(sbt_tc.cut)

	sbt_ag_tc.freq = sbt_ag.freq + sbt_tc.freq	

	sumfreq0=sum(sbt_ag_tc.freq)
	cumfreq0 = c(0, cumsum(sbt_ag_tc.freq))
	numerator=sumfreq0-cumfreq0
	
	sbt_ag_tc_p = numerator/sumfreq0
	sbt_ag_tc_p <- sbt_ag_tc_p[1:length(sbt_ag_tc_p)-1] #

	ag_tc_frq_pval_table1 <- cbind(right_points, sbt_ag_tc.freq, sbt_ag_tc_p)
#	---------------------------------------------------------------------------------------------------	#	

#	---------------------------------------------------------------------------------------------------	#	
	sbt_ca.cut = cut(sbt_ca, brk)
#	, right=FALSE) 
	sbt_ca.freq = table(sbt_ca.cut)
	
	sbt_gt.cut = cut(sbt_gt, brk)
#	, right=FALSE)
	sbt_gt.freq = table(sbt_gt.cut)

	sbt_ca_gt.freq = sbt_ca.freq + sbt_gt.freq
	
	sumfreq0=sum(sbt_ca_gt.freq)
	cumfreq0 = c(0, cumsum(sbt_ca_gt.freq))
	numerator=sumfreq0-cumfreq0
	
	sbt_ca_gt_p = numerator/sumfreq0
	sbt_ca_gt_p <- sbt_ca_gt_p[1:length(sbt_ca_gt_p)-1] #

	ca_gt_frq_pval_table1 <- cbind(right_points, sbt_ca_gt.freq, sbt_ca_gt_p)
#	---------------------------------------------------------------------------------------------------	#	

#	---------------------------------------------------------------------------------------------------	#	
	sbt_ct.cut = cut(sbt_ct, brk)
#	, right=FALSE) 
	sbt_ct.freq = table(sbt_ct.cut)
	
	sbt_ga.cut = cut(sbt_ga, brk)
#	, right=FALSE)
	sbt_ga.freq = table(sbt_ga.cut)

	sbt_ct_ga.freq = sbt_ct.freq + sbt_ga.freq	

	sumfreq0=sum(sbt_ct_ga.freq)
	cumfreq0 = c(0, cumsum(sbt_ct_ga.freq))
	numerator=sumfreq0-cumfreq0
	
	sbt_ct_ga_p = numerator/sumfreq0
	sbt_ct_ga_p <- sbt_ct_ga_p[1:length(sbt_ct_ga_p)-1] #

	ct_ga_frq_pval_table1 <- cbind(right_points, sbt_ct_ga.freq, sbt_ct_ga_p)
#	---------------------------------------------------------------------------------------------------	#	

#	---------------------------------------------------------------------------------------------------	#	
	sbt_cg.cut = cut(sbt_cg, brk)
#	, right=FALSE) 
	sbt_cg.freq = table(sbt_cg.cut)
	
	sbt_gc.cut = cut(sbt_gc, brk)
#	, right=FALSE)
	sbt_gc.freq = table(sbt_gc.cut)

	sbt_cg_gc.freq = sbt_cg.freq + sbt_gc.freq	
	
	sumfreq0=sum(sbt_cg_gc.freq)
	cumfreq0 = c(0, cumsum(sbt_cg_gc.freq))
	numerator=sumfreq0-cumfreq0

	sbt_cg_gc_p = numerator/sumfreq0
	sbt_cg_gc_p <- sbt_cg_gc_p[1:length(sbt_cg_gc_p)-1] #

	cg_gc_frq_pval_table1 <- cbind(right_points, sbt_cg_gc.freq, sbt_cg_gc_p)

	frq_pval_table <- cbind(right_points, sbt_ac_tg.freq, sbt_ac_tg_p, sbt_at_ta.freq, sbt_at_ta_p, sbt_ag_tc.freq, sbt_ag_tc_p, sbt_ca_gt.freq, sbt_ca_gt_p, sbt_ct_ga.freq, sbt_ct_ga_p, sbt_cg_gc.freq, sbt_cg_gc_p)
	write.table(frq_pval_table, file= paste(destination, ".fpvtb1.tsv", sep=""), sep="\t", row.names=FALSE)





	
##	xval=0
##	yval=0
	j=0
	#chr  <- dfrm$chromosome;
	#site <- dfrm$site;
	#chrsite=paste(chr, " ", site, sep="")
	# dfrm=subset(dfrm, !chrsite%in%excl) # excludes the sites of interest excl is a vector generated from the primers file(element : chr, space, site of interest. exemple : "2 2345673" means chromosome: 2 site: 2345673
	dfrm <- vectorOfTables[[2]]
	nbe  <- as.numeric(dfrm$number_base_expected);
	dfrm <- subset(dfrm, nbe>0);
	chr  <- dfrm$chromosome;
	be   <- dfrm$base_expected;
	alt  <- dfrm$alternative;
	site <- as.numeric(dfrm$site);
	tA   <- as.numeric(dfrm$number_of_As);
	tC   <- as.numeric(dfrm$number_of_Cs);
	tT   <- as.numeric(dfrm$number_of_Ts);
	tG   <- as.numeric(dfrm$number_of_Gs);
	nbe  <- as.numeric(dfrm$number_base_expected);
	xv   <- as.numeric(dfrm$x);
	yv   <- as.numeric(dfrm$y);
	dv = nbe;
	crlimit=1.5

	# prf=0.95 # percentage of reference
	xl=c(-0.0,1.00)
##	xval=0
##	yval=0
	j=0
	i=1

sbt_ac <- subset(tC/dv, (be=='A' & alt=='A' & tA/dv>prf & xv>xval & yv>yval & tC/dv < range & dv>cover)) # red
sbt_at <- subset(tT/dv, (be=='A' & alt=='A' & tA/dv>prf & xv>xval & yv>yval & tT/dv < range & dv>cover)) # blue 
sbt_ag <- subset(tG/dv, (be=='A' & alt=='A' & tA/dv>prf & xv>xval & yv>yval & tG/dv < range & dv>cover)) # green

sbt_tg <- subset(tG/dv, (be=='T' & alt=='T' & tT/dv>prf & xv>xval & yv>yval & tG/dv < range & dv>cover)) # purple
sbt_ta <- subset(tA/dv, (be=='T' & alt=='T' & tT/dv>prf & xv>xval & yv>yval & tA/dv < range & dv>cover)) # cyan
sbt_tc <- subset(tC/dv, (be=='T' & alt=='T' & tT/dv>prf & xv>xval & yv>yval & tC/dv < range & dv>cover)) # orange

sbt_ca <- subset(tA/dv, (be=='C' & alt=='C' & tC/dv>prf & xv>xval & yv>yval & tA/dv < range & dv>cover)) # red 
sbt_ct <- subset(tT/dv, (be=='C' & alt=='C' & tC/dv>prf & xv>xval & yv>yval & tT/dv < range & dv>cover)) # blue 
sbt_cg <- subset(tG/dv, (be=='C' & alt=='C' & tC/dv>prf & xv>xval & yv>yval & tG/dv < range & dv>cover)) # green

sbt_gt <- subset(tT/dv, (be=='G' & alt=='G' & tG/dv>prf & xv>xval & yv>yval & tT/dv < range & dv>cover)) # purple
sbt_ga <- subset(tA/dv, (be=='G' & alt=='G' & tG/dv>prf & xv>xval & yv>yval & tA/dv < range & dv>cover)) # cyan 
sbt_gc <- subset(tC/dv, (be=='G' & alt=='G' & tG/dv>prf & xv>xval & yv>yval & tC/dv < range & dv>cover)) # orange

	sbt_ac.cut = cut(sbt_ac, brk)
#	, right=FALSE) 
	sbt_ac.freq = table(sbt_ac.cut)
	
	sbt_tg.cut = cut(sbt_tg, brk)
#	, right=FALSE)
	sbt_tg.freq = table(sbt_tg.cut)

	sbt_ac_tg.freq = sbt_ac.freq + sbt_tg.freq	
	
	sumfreq0=sum(sbt_ac_tg.freq)
	cumfreq0 = c(0, cumsum(sbt_ac_tg.freq))
	numerator=sumfreq0-cumfreq0
	
	sbt_ac_tg_p = numerator/sumfreq0
	sbt_ac_tg_p <- sbt_ac_tg_p[1:length(sbt_ac_tg_p)-1] #
	
	ac_tg_frq_pval_table1 <- cbind(right_points, sbt_ac_tg.freq, sbt_ac_tg_p)
#	---------------------------------------------------------------------------------------------------	#	

#	---------------------------------------------------------------------------------------------------	#	
	sbt_at.cut = cut(sbt_at, brk)
#	, right=FALSE) 
	sbt_at.freq = table(sbt_at.cut)
	
	sbt_ta.cut = cut(sbt_ta, brk)
#	, right=FALSE)
	sbt_ta.freq = table(sbt_ta.cut)
	
	sbt_at_ta.freq = sbt_at.freq + sbt_ta.freq	

	sumfreq0=sum(sbt_at_ta.freq)
	cumfreq0 = c(0, cumsum(sbt_at_ta.freq))
	numerator=sumfreq0-cumfreq0
	
	sbt_at_ta_p = numerator/sumfreq0
	sbt_at_ta_p <- sbt_at_ta_p[1:length(sbt_at_ta_p)-1] #

	at_ta_frq_pval_table1 <- cbind(right_points, sbt_at_ta.freq, sbt_at_ta_p)
#	---------------------------------------------------------------------------------------------------	#	

#	---------------------------------------------------------------------------------------------------	#	
	sbt_ag.cut = cut(sbt_ag, brk)
#	, right=FALSE)
	sbt_ag.freq = table(sbt_ag.cut)

	sbt_tc.cut = cut(sbt_tc, brk)
#	, right=FALSE)
	sbt_tc.freq = table(sbt_tc.cut)

	sbt_ag_tc.freq = sbt_ag.freq + sbt_tc.freq	

	sumfreq0=sum(sbt_ag_tc.freq)
	cumfreq0 = c(0, cumsum(sbt_ag_tc.freq))
	numerator=sumfreq0-cumfreq0
	
	sbt_ag_tc_p = numerator/sumfreq0
	sbt_ag_tc_p <- sbt_ag_tc_p[1:length(sbt_ag_tc_p)-1] #

	ag_tc_frq_pval_table1 <- cbind(right_points, sbt_ag_tc.freq, sbt_ag_tc_p)
#	---------------------------------------------------------------------------------------------------	#	

#	---------------------------------------------------------------------------------------------------	#	
	sbt_ca.cut = cut(sbt_ca, brk)
#	, right=FALSE) 
	sbt_ca.freq = table(sbt_ca.cut)
	
	sbt_gt.cut = cut(sbt_gt, brk)
#	, right=FALSE)
	sbt_gt.freq = table(sbt_gt.cut)

	sbt_ca_gt.freq = sbt_ca.freq + sbt_gt.freq
	
	sumfreq0=sum(sbt_ca_gt.freq)
	cumfreq0 = c(0, cumsum(sbt_ca_gt.freq))
	numerator=sumfreq0-cumfreq0
	
	sbt_ca_gt_p = numerator/sumfreq0
	sbt_ca_gt_p <- sbt_ca_gt_p[1:length(sbt_ca_gt_p)-1] #

	ca_gt_frq_pval_table1 <- cbind(right_points, sbt_ca_gt.freq, sbt_ca_gt_p)
#	---------------------------------------------------------------------------------------------------	#	

#	---------------------------------------------------------------------------------------------------	#	
	sbt_ct.cut = cut(sbt_ct, brk)
#	, right=FALSE) 
	sbt_ct.freq = table(sbt_ct.cut)
	
	sbt_ga.cut = cut(sbt_ga, brk)
#	, right=FALSE)
	sbt_ga.freq = table(sbt_ga.cut)

	sbt_ct_ga.freq = sbt_ct.freq + sbt_ga.freq	

	sumfreq0=sum(sbt_ct_ga.freq)
	cumfreq0 = c(0, cumsum(sbt_ct_ga.freq))
	numerator=sumfreq0-cumfreq0
	
	sbt_ct_ga_p = numerator/sumfreq0
	sbt_ct_ga_p <- sbt_ct_ga_p[1:length(sbt_ct_ga_p)-1] #

	ct_ga_frq_pval_table1 <- cbind(right_points, sbt_ct_ga.freq, sbt_ct_ga_p)
#	---------------------------------------------------------------------------------------------------	#	

#	---------------------------------------------------------------------------------------------------	#	
	sbt_cg.cut = cut(sbt_cg, brk)
#	, right=FALSE) 
	sbt_cg.freq = table(sbt_cg.cut)
	
	sbt_gc.cut = cut(sbt_gc, brk)
#	, right=FALSE)
	sbt_gc.freq = table(sbt_gc.cut)

	sbt_cg_gc.freq = sbt_cg.freq + sbt_gc.freq	
	
	sumfreq0=sum(sbt_cg_gc.freq)
	cumfreq0 = c(0, cumsum(sbt_cg_gc.freq))
	numerator=sumfreq0-cumfreq0

	sbt_cg_gc_p = numerator/sumfreq0
	sbt_cg_gc_p <- sbt_cg_gc_p[1:length(sbt_cg_gc_p)-1] #

	cg_gc_frq_pval_table1 <- cbind(right_points, sbt_cg_gc.freq, sbt_cg_gc_p)

	frq_pval_table <- cbind(right_points, sbt_ac_tg.freq, sbt_ac_tg_p, sbt_at_ta.freq, sbt_at_ta_p, sbt_ag_tc.freq, sbt_ag_tc_p, sbt_ca_gt.freq, sbt_ca_gt_p, sbt_ct_ga.freq, sbt_ct_ga_p, sbt_cg_gc.freq, sbt_cg_gc_p)
	write.table(frq_pval_table, file= paste(destination, ".fpvtb2.tsv", sep=""), sep="\t", row.names=FALSE)

}

# Create position substitution rate  graph
if( (oneortwo==1)  & (type=="mutrate")  & (graphortable=="posgrout") )
{
xl=c(-0.0001,range)
# yl=c(0,10000)
brk <- seq(0,range, 6.25e-5)
#brk <- seq(0,range, 6.25e-5/range)
#brk <- seq(0,0.05, 6.25e-5)
# xval=overlap
# yval=overlap
	prf=0
	old.par <- par(mfrow=c(1, 1))

	print("mutrate posgrout 1");
	dfrm <- vectorOfTables[[1]]
	nbe  <- as.numeric(dfrm$number_base_expected);
	dfrm <- subset(dfrm, nbe>0);
	chr  <- dfrm$chromosome;
	be   <- dfrm$base_expected;
	alt  <- dfrm$alternative;
	site <- as.numeric(dfrm$site);
	tA   <- as.numeric(dfrm$number_of_As);
	tC   <- as.numeric(dfrm$number_of_Cs);
	tT   <- as.numeric(dfrm$number_of_Ts);
	tG   <- as.numeric(dfrm$number_of_Gs);
	nbe  <- as.numeric(dfrm$number_base_expected);
	xv   <- as.numeric(dfrm$x);
	yv   <- as.numeric(dfrm$y);

	dv = nbe
	#print(range)

	sbt_ac <- subset(tC/dv, (be=='A' & tA/dv>prf & xv>xval & yv>yval & tC/dv < range & dv>cover)) # blue
	sbt_at <- subset(tT/dv, (be=='A' & tA/dv>prf & xv>xval & yv>yval & tT/dv < range & dv>cover)) # cyan 
	sbt_ag <- subset(tG/dv, (be=='A' & tA/dv>prf & xv>xval & yv>yval & tG/dv < range & dv>cover)) # violet
	sbt_tg <- subset(tG/dv, (be=='T' & tT/dv>prf & xv>xval & yv>yval & tG/dv < range & dv>cover)) # yellow
	sbt_ta <- subset(tA/dv, (be=='T' & tT/dv>prf & xv>xval & yv>yval & tA/dv < range & dv>cover)) # green
	sbt_tc <- subset(tC/dv, (be=='T' & tT/dv>prf & xv>xval & yv>yval & tC/dv < range & dv>cover)) # red
	sbt_ca <- subset(tA/dv, (be=='C' & tC/dv>prf & xv>xval & yv>yval & tA/dv < range & dv>cover)) # purple 
	sbt_ct <- subset(tT/dv, (be=='C' & tC/dv>prf & xv>xval & yv>yval & tT/dv < range & dv>cover)) # pink
	sbt_cg <- subset(tG/dv, (be=='C' & tC/dv>prf & xv>xval & yv>yval & tG/dv < range & dv>cover)) # brown
	sbt_gt <- subset(tT/dv, (be=='G' & tG/dv>prf & xv>xval & yv>yval & tT/dv < range & dv>cover)) # black
	sbt_ga <- subset(tA/dv, (be=='G' & tG/dv>prf & xv>xval & yv>yval & tA/dv < range & dv>cover)) # orange
	sbt_gc <- subset(tC/dv, (be=='G' & tG/dv>prf & xv>xval & yv>yval & tC/dv < range & dv>cover)) # grey


	h11s <- hist(sbt_ac, breaks=brk, plot=FALSE)		# ac-tg
	h12s <- hist(sbt_at, breaks=brk, plot=FALSE)		# at-ta
	h13s <- hist(sbt_ag, breaks=brk, plot=FALSE)		# ag-tc
	
	h21s <- hist(sbt_tg, breaks=brk, plot=FALSE)		# tg-ac
	h22s <- hist(sbt_ta, breaks=brk, plot=FALSE)		# ta-at
	h23s <- hist(sbt_tc, breaks=brk, plot=FALSE)		# tc-ag
	
	h31s <- hist(sbt_ca, breaks=brk, plot=FALSE)		# ca-gt
	h32s <- hist(sbt_ct, breaks=brk, plot=FALSE)		# ct-ga
	h33s <- hist(sbt_cg, breaks=brk, plot=FALSE)		# cg-gc

	h41s <- hist(sbt_gt, breaks=brk, plot=FALSE)		# gt-ac
	h42s <- hist(sbt_ga, breaks=brk, plot=FALSE)		# ga-ct
	h43s <- hist(sbt_gc, breaks=brk, plot=FALSE)		# gc-cg	

	ysmaxi = max(h11s$counts, h21s$counts, h12s$counts, h22s$counts, h13s$counts, h23s$counts, h31s$counts, h41s$counts, h32s$counts, h42s$counts, h33s$counts, h43s$counts)
	ysm = ysmaxi + 0.1*ysmaxi
	ysl=c(0,ysm)
	ttl  = paste("Mutation Rate: sites size>=", (cover+1), ", reference > ", (prf), ", y=counts\n", sep="")
	subi = paste("A->C = blue, T->G = yellow, A->T = cyan, T->A = green, A->G =violet, T->C = red\nC->A = purple, G->T = black, C->T = pink, G->A = orange, C->G = brown, G->C = grey\n", sep="")

	if(be[[1]][1]=='A')
	{
	plot(x=h11s$mids, y=h11s$counts, lwd = 1, col="blue", type="l", yaxt="n", xlim=xl, ylim = ysl, xaxt="n", sub = subi, cex.sub = 0.75, font.sub = 1, xlab=NA,  cex.main = 0.75, main=paste(ttl, "file: ", tblnms[1], sep=""), ylab=tblnms[1])
	lines(x=h12s$mids, y=h12s$counts, lwd = 1, col="cyan")
	lines(x=h13s$mids, y=h13s$counts, lwd = 1, col="violet")
	}

	if(be[[1]][1]=='T')
	{
	plot(x=h21s$mids, y=h21s$counts, lwd = 1, col="yellow", type="l", yaxt="n", xlim=xl, ylim = ysl, xaxt="n", sub = subi, cex.sub = 0.75, font.sub = 1, xlab=NA,  cex.main = 0.75, main=paste(ttl, "file: ", tblnms[1], sep=""), ylab=tblnms[1])
	lines(x=h22s$mids, y=h22s$counts, lwd = 1, col="green")
	lines(x=h23s$mids, y=h23s$counts, lwd = 1, col="red")
	}

	if(be[[1]][1]=='C')
	{
	plot(x=h31s$mids, y=h31s$counts, lwd = 1, col="purple", type="l", yaxt="n", xlim=xl, ylim = ysl, xaxt="n", sub = subi, cex.sub = 0.75, font.sub = 1, xlab=NA,  cex.main = 0.75, main=paste(ttl, "file: ", tblnms[1], sep=""), ylab=tblnms[1])
	lines(x=h32s$mids, y=h32s$counts, lwd = 1, col="pink")
	lines(x=h33s$mids, y=h33s$counts, lwd = 1, col="brown")
	}

	if(be[[1]][1]=='G')
	{
	plot(x=h41s$mids, y=h41s$counts, lwd = 1, col="black", type="l", yaxt="n", xlim=xl, ylim = ysl, xaxt="n", sub = subi, cex.sub = 0.75, font.sub = 1, xlab=NA,  cex.main = 0.75, main=paste(ttl, "file: ", tblnms[1], sep=""), ylab=tblnms[1])
	lines(x=h42s$mids, y=h42s$counts, lwd = 1, col="orange")
	lines(x=h43s$mids, y=h43s$counts, lwd = 1, col="grey")
	}

	axis(1, cex.axis = .75)
	axis(2, pos=c(0,0), cex.axis = .75)

	par(old.par)
	dev.off()
}

# Create position empirical pvalue graph 
if( (oneortwo==1)  & (type=="pvalue")  & (graphortable=="posgrout") )
{
xl=c(-0.0001,range)

	print("pvalue posgrout 1");
	dfrm <- vectorOfTables[[1]]
	nbe  <- as.numeric(dfrm$number_base_expected);
	dfrm <- subset(dfrm, nbe>0);
	chr  <- dfrm$chromosome;
	be   <- dfrm$base_expected;
	alt  <- dfrm$alternative;
	site <- as.numeric(dfrm$site);
	tA   <- as.numeric(dfrm$number_of_As);
	tC   <- as.numeric(dfrm$number_of_Cs);
	tT   <- as.numeric(dfrm$number_of_Ts);
	tG   <- as.numeric(dfrm$number_of_Gs);
	nbe  <- as.numeric(dfrm$number_base_expected);
	xv   <- as.numeric(dfrm$x);
	yv   <- as.numeric(dfrm$y);

	dv = nbe;
	#print(range);

	crlimit=1.5
	xl=c(-0.0001,range)
	#xl=c(-0.0001,0.005)
	yl=c(0,10000)
#	brk <- seq(0,0.005, 6.25e-5)
#	brk <- seq(0,1.0, 10e-5)
	right_points <- seq(0,1.0, 10e-5)
#	brk <- seq(-0.0001,1.0, 10e-5)
#	xl=c(-0.0,0.03)
##	xval=0
##	yval=0
	j=0
	i=1
	prf=0
#	sbt_ac <- subset(tC/dv, (be=='A' & alt=='C' & tA/dv>prf & xv>xval & yv>yval & tC/dv < 0.005 & dv>cover)) # blue
#	sbt_at <- subset(tT/dv, (be=='A' & alt=='T' & tA/dv>prf & xv>xval & yv>yval & tT/dv < 0.005 & dv>cover)) # cyan 
#	sbt_ag <- subset(tG/dv, (be=='A' & alt=='G' & tA/dv>prf & xv>xval & yv>yval & tG/dv < 0.005 & dv>cover)) # violet

#	sbt_tg <- subset(tG/dv, (be=='T' & alt=='G' & tT/dv>prf & xv>xval & yv>yval & tG/dv < 0.005 & dv>cover)) # yellow
#	sbt_ta <- subset(tA/dv, (be=='T' & alt=='A' & tT/dv>prf & xv>xval & yv>yval & tA/dv < 0.005 & dv>cover)) # green
#	sbt_tc <- subset(tC/dv, (be=='T' & alt=='C' & tT/dv>prf & xv>xval & yv>yval & tC/dv < 0.005 & dv>cover)) # red

#	sbt_ca <- subset(tA/dv, (be=='C' & alt=='A' & tC/dv>prf & xv>xval & yv>yval & tA/dv < 0.005 & dv>cover)) # purple 
#	sbt_ct <- subset(tT/dv, (be=='C' & alt=='T' & tC/dv>prf & xv>xval & yv>yval & tT/dv < 0.005 & dv>cover)) # pink
#	sbt_cg <- subset(tG/dv, (be=='C' & alt=='G' & tC/dv>prf & xv>xval & yv>yval & tG/dv < 0.005 & dv>cover)) # brown

#	sbt_gt <- subset(tT/dv, (be=='G' & alt=='T' & tG/dv>prf & xv>xval & yv>yval & tT/dv < 0.005 & dv>cover)) # black
#	sbt_ga <- subset(tA/dv, (be=='G' & alt=='A' & tG/dv>prf & xv>xval & yv>yval & tA/dv < 0.005 & dv>cover)) # orange
#	sbt_gc <- subset(tC/dv, (be=='G' & alt=='C' & tG/dv>prf & xv>xval & yv>yval & tC/dv < 0.005 & dv>cover)) # grey

#	sbt_ac <- subset(tC/dv, (be=='A' & alt=='C' & tA/dv>prf & xv>xval & yv>yval & tC/dv < range & dv>cover)) # blue
#	sbt_at <- subset(tT/dv, (be=='A' & alt=='T' & tA/dv>prf & xv>xval & yv>yval & tT/dv < range & dv>cover)) # cyan 
#	sbt_ag <- subset(tG/dv, (be=='A' & alt=='G' & tA/dv>prf & xv>xval & yv>yval & tG/dv < range & dv>cover)) # violet
#	sbt_tg <- subset(tG/dv, (be=='T' & alt=='G' & tT/dv>prf & xv>xval & yv>yval & tG/dv < range & dv>cover)) # yellow
#	sbt_ta <- subset(tA/dv, (be=='T' & alt=='A' & tT/dv>prf & xv>xval & yv>yval & tA/dv < range & dv>cover)) # green
#	sbt_tc <- subset(tC/dv, (be=='T' & alt=='C' & tT/dv>prf & xv>xval & yv>yval & tC/dv < range & dv>cover)) # red
#	sbt_ca <- subset(tA/dv, (be=='C' & alt=='A' & tC/dv>prf & xv>xval & yv>yval & tA/dv < range & dv>cover)) # purple 
#	sbt_ct <- subset(tT/dv, (be=='C' & alt=='T' & tC/dv>prf & xv>xval & yv>yval & tT/dv < range & dv>cover)) # pink
#	sbt_cg <- subset(tG/dv, (be=='C' & alt=='G' & tC/dv>prf & xv>xval & yv>yval & tG/dv < range & dv>cover)) # brown
#	sbt_gt <- subset(tT/dv, (be=='G' & alt=='T' & tG/dv>prf & xv>xval & yv>yval & tT/dv < range & dv>cover)) # black
#	sbt_ga <- subset(tA/dv, (be=='G' & alt=='A' & tG/dv>prf & xv>xval & yv>yval & tA/dv < range & dv>cover)) # orange
#	sbt_gc <- subset(tC/dv, (be=='G' & alt=='C' & tG/dv>prf & xv>xval & yv>yval & tC/dv < range & dv>cover)) # grey


	sbt_ac <- subset(tC/dv, (be=='A' & tA/dv>prf & xv>xval & yv>yval & tC/dv < range & dv>cover)) # blue
	sbt_at <- subset(tT/dv, (be=='A' & tA/dv>prf & xv>xval & yv>yval & tT/dv < range & dv>cover)) # cyan 
	sbt_ag <- subset(tG/dv, (be=='A' & tA/dv>prf & xv>xval & yv>yval & tG/dv < range & dv>cover)) # violet
	sbt_tg <- subset(tG/dv, (be=='T' & tT/dv>prf & xv>xval & yv>yval & tG/dv < range & dv>cover)) # yellow
	sbt_ta <- subset(tA/dv, (be=='T' & tT/dv>prf & xv>xval & yv>yval & tA/dv < range & dv>cover)) # green
	sbt_tc <- subset(tC/dv, (be=='T' & tT/dv>prf & xv>xval & yv>yval & tC/dv < range & dv>cover)) # red
	sbt_ca <- subset(tA/dv, (be=='C' & tC/dv>prf & xv>xval & yv>yval & tA/dv < range & dv>cover)) # purple 
	sbt_ct <- subset(tT/dv, (be=='C' & tC/dv>prf & xv>xval & yv>yval & tT/dv < range & dv>cover)) # pink
	sbt_cg <- subset(tG/dv, (be=='C' & tC/dv>prf & xv>xval & yv>yval & tG/dv < range & dv>cover)) # brown
	sbt_gt <- subset(tT/dv, (be=='G' & tG/dv>prf & xv>xval & yv>yval & tT/dv < range & dv>cover)) # black
	sbt_ga <- subset(tA/dv, (be=='G' & tG/dv>prf & xv>xval & yv>yval & tA/dv < range & dv>cover)) # orange
	sbt_gc <- subset(tC/dv, (be=='G' & tG/dv>prf & xv>xval & yv>yval & tC/dv < range & dv>cover)) # grey


#sbt_ga <- subset(tA/dv, (be=='G' & alt=='A' & tG/dv>prf & xv>xval & yv>yval & dv>cover)) # orange 

	sbt_ac.cut = cut(sbt_ac, brk, right=FALSE)
	sbt_ac.freq = table(sbt_ac.cut)
	cumfreq0 = c(0, cumsum(sbt_ac.freq))
	sumfreq0=sum(sbt_ac.freq)
	sbt_ac1 = cumfreq0/sumfreq0
	sbt_ac1_p = 1-cumfreq0/sumfreq0
	
	sbt_at.cut = cut(sbt_at, brk, right=FALSE)
	sbt_at.freq = table(sbt_at.cut)
	cumfreq0 = c(0, cumsum(sbt_at.freq))
	sumfreq0=sum(sbt_at.freq)
	sbt_at1 = cumfreq0/sumfreq0
	sbt_at1_p = 1-cumfreq0/sumfreq0
	
	sbt_ag.cut = cut(sbt_ag, brk, right=FALSE)
	sbt_ag.freq = table(sbt_ag.cut)
	cumfreq0 = c(0, cumsum(sbt_ag.freq))
	sumfreq0=sum(sbt_ag.freq)
	sbt_ag1 = cumfreq0/sumfreq0
	sbt_ag1_p = 1-cumfreq0/sumfreq0

	sbt_tg.cut = cut(sbt_tg, brk, right=FALSE)
	sbt_tg.freq = table(sbt_tg.cut)
	cumfreq0 = c(0, cumsum(sbt_tg.freq))
	sumfreq0=sum(sbt_tg.freq)
	sbt_tg1 = cumfreq0/sumfreq0
	sbt_tg1_p = 1-cumfreq0/sumfreq0

	sbt_ta.cut = cut(sbt_ta, brk, right=FALSE)
	sbt_ta.freq = table(sbt_ta.cut)
	cumfreq0 = c(0, cumsum(sbt_ta.freq))
	sumfreq0=sum(sbt_ta.freq)
	sbt_ta1 = cumfreq0/sumfreq0
	sbt_ta1_p = 1-cumfreq0/sumfreq0

	sbt_tc.cut = cut(sbt_tc, brk, right=FALSE)
	sbt_tc.freq = table(sbt_tc.cut)
	cumfreq0 = c(0, cumsum(sbt_tc.freq))
	sumfreq0=sum(sbt_tc.freq)
	sbt_tc1 = cumfreq0/sumfreq0
	sbt_tc1_p = 1-cumfreq0/sumfreq0

	sbt_ca.cut = cut(sbt_ca, brk, right=FALSE)
	sbt_ca.freq = table(sbt_ca.cut)
	cumfreq0 = c(0, cumsum(sbt_ca.freq))
	sumfreq0=sum(sbt_ca.freq)
	sbt_ca1 = cumfreq0/sumfreq0
	sbt_ca1_p = 1-cumfreq0/sumfreq0

	sbt_ct.cut = cut(sbt_ct, brk, right=FALSE)
	sbt_ct.freq = table(sbt_ct.cut)
	cumfreq0 = c(0, cumsum(sbt_ct.freq))
	sumfreq0=sum(sbt_ct.freq)
	sbt_ct1 = cumfreq0/sumfreq0
	sbt_ct1_p = 1-cumfreq0/sumfreq0

	sbt_cg.cut = cut(sbt_cg, brk, right=FALSE)
	sbt_cg.freq = table(sbt_cg.cut)
	cumfreq0 = c(0, cumsum(sbt_cg.freq))
	sumfreq0=sum(sbt_cg.freq)
	sbt_cg1 = cumfreq0/sumfreq0
	sbt_cg1_p = 1-cumfreq0/sumfreq0

	sbt_gt.cut = cut(sbt_gt, brk, right=FALSE)
	sbt_gt.freq = table(sbt_gt.cut)
	cumfreq0 = c(0, cumsum(sbt_gt.freq))
	sumfreq0=sum(sbt_gt.freq)
	sbt_gt1 = cumfreq0/sumfreq0
	sbt_gt1_p = 1-cumfreq0/sumfreq0

	sbt_ga.cut = cut(sbt_ga, brk, right=FALSE)
	sbt_ga.freq = table(sbt_ga.cut)
	cumfreq0 = c(0, cumsum(sbt_ga.freq))
	sumfreq0=sum(sbt_ga.freq)
	sbt_ga1 = cumfreq0/sumfreq0
	sbt_ga1_p = 1-cumfreq0/sumfreq0

	sbt_gc.cut = cut(sbt_gc, brk, right=FALSE)
	sbt_gc.freq = table(sbt_gc.cut)
	cumfreq0 = c(0, cumsum(sbt_gc.freq))
	sumfreq0=sum(sbt_gc.freq)
	sbt_gc1 = cumfreq0/sumfreq0
	sbt_gc1_p = 1-cumfreq0/sumfreq0

	ymaxi = max(log10(sbt_ac),log10(sbt_at),log10(sbt_ag),log10(sbt_tg),log10(sbt_ta),log10(sbt_tc),log10(sbt_ca),log10(sbt_ct),log10(sbt_cg),log10(sbt_gt),log10(sbt_ga),log10(sbt_gc))
	ym = ymaxi + 0.1*ymaxi
	
	yl=c(-3,0)

	ttl  = paste("p_value: sites size>=", (cover+1), ", reference > ", (prf), "%,  y=log(p-value)\n",sep="");
	subia = paste("A->C = blue, T->G = yellow, A->T = cyan, T->A = green, A->G =violet, T->C = red\nC->A = purple, G->T = black, C->T = pink, G->A = orange, C->G = brown, G->C = grey\n", sep="")
#	tblnms = c("table01.tsv")

	if(be[[1]][1]=='A')
	{
	 plot(x=brk, y=(lf)*log10(sbt_ac1_p), col="blue", lwd = 1, type="l", ylim=yl, xlim=xl, yaxt='n', xaxt="n", sub = subia, cex.sub = 0.70, font.sub = 1, xlab=NA, main=paste(ttl,"file: ", tblnms[1], sep=""), cex.main=0.70, ylab=tblnms[1])
	lines(x=brk, y=(lf)*log10(sbt_at1_p), col="cyan", lwd = 1)
	lines(x=brk, y=(lf)*log10(sbt_ag1_p), col="violet", lwd = 1)
	}

	if(be[[1]][1]=='T')
	{
	 plot(x=brk, y=(lf)*log10(sbt_tg1_p), col="yellow", lwd = 1, type="l", ylim=yl, xlim=xl, yaxt='n', xaxt="n", sub = subia, cex.sub = 0.70, font.sub = 1, xlab=NA, main=paste(ttl,"file: ", tblnms[1], sep=""), cex.main=0.70, ylab=tblnms[1])
	lines(x=brk, y=(lf)*log10(sbt_ta1_p), col="green", lwd = 1)
	lines(x=brk, y=(lf)*log10(sbt_tc1_p), col="red", lwd = 1)
	}

	if(be[[1]][1]=='C')
	{
	 plot(x=brk, y=(lf)*log10(sbt_ca1_p), col="purple", lwd = 1, type="l", ylim=yl, xlim=xl, yaxt='n', xaxt="n", sub = subia, cex.sub = 0.70, font.sub = 1, xlab=NA, main=paste(ttl,"file: ", tblnms[1], sep=""), cex.main=0.70, ylab=tblnms[1])
	lines(x=brk, y=(lf)*log10(sbt_ct1_p), col="pink", lwd = 1)
	lines(x=brk, y=(lf)*log10(sbt_cg1_p), col="brown", lwd = 1)
	}

	if(be[[1]][1]=='G')
	{
	 plot(x=brk, y=(lf)*log10(sbt_gt1_p), col="black", lwd = 1, type="l", ylim=yl, xlim=xl, yaxt='n', xaxt="n", sub = subia, cex.sub = 0.70, font.sub = 1, xlab=NA, main=paste(ttl,"file: ", tblnms[1], sep=""), cex.main=0.70, ylab=tblnms[1])
	lines(x=brk, y=(lf)*log10(sbt_ga1_p), col="orange", lwd = 1)
	lines(x=brk, y=(lf)*log10(sbt_gc1_p), col="grey", lwd = 1)
	}
	axis(1, cex.axis = .75)
	#axis(2, pos=c(0,0), cex.axis = .5, labels = FALSE)
	axis(2, pos=c(0,-10), cex.axis = .5, yaxp = c(-10, 0, 10), las=2)
	dev.off()
}
