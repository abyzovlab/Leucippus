/** Description of Leucippus 
 * @author  1: Alexej Abyzov
 * @author  2: Nikolaos Vasmatzis
 * @authors 3: Others
 * @version 1.0 Build 08/06/2015.
 * Description : Software for Noise Estimation
 */
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Hashtable;
import java.util.Set;
import java.util.List;
import java.util.Random;
import java.util.Vector;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import java.io.RandomAccessFile;
import java.util.UUID;
import java.net.URL;
import java.util.Enumeration;
public class Leucippus
{
	private static String[][] noisetable = new String[0][0];
// 	Two dimensional String array used to hold table results
	private static List<Double[]> dbtb = new ArrayList<Double[]>();
// 	Data structure that is used to host arrays with various size.
// 	that contain coefficients from Pascal's Triangle.
// 	It is used in calculating the statistics for choosing the window 
// 	with the best overlap.
	private static String q = "!\"#$%&\'()*+,-./0123456789:;<=>?@ABCDE"
	        +"FGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~";
	private static char[] qar = q.toCharArray();
// 	the qar is a base quality grading 'char' array;
// 	the index of qar array provides with the actual base quality. 
	private static int[] percmis = new int[11];
// 	integer array that holds the number of reads that present within groups of 
// 	percentage of mismatches.
	private static Vector<String> readswithindels = new Vector<String>();
//	Vector that holds all reads that present indels. The first element of a
//  particular group contains information about the interval that was used to
// 	(coordinates, sequence, and locus) to retrieve the reads from 'bam' files.
//	Each read line of the vector contains the origin (bam file name) the read name
//  the modified read(according to cigar and position, the cigar, and the position.
	private static Vector<String> missed_sites = new Vector<String>();
	/**
	 * @param args : String array that holds user arguments
	 * @throws InterruptedException
	 * @throws IOException
	 */
	public static void main(String[] args) throws InterruptedException, IOException
	{ 
// 		Main method (START)
//		System.out.println("main method start!");	
// 		Information printed on screen(START)
// 		General Information (START)
		String leu_gen_nfo = "\n\nProgram: Leucippus (Mosaics: Site Noise-Estimation, Validation, Identification)\n"
				+ "Version: 0.1.0\n\n"
				+ "Usage:   Leucippus <command> [options]\n\n"
				+ "Command: frag		create long reads\n"
				+ "	 noisetab	create noise table\n"
				+ "	 graph		create graph\n"
				+ "	 decide		decide [somatic | germline  | unknown | omitted]\n"
				+ "	 extract 	extract reads that present alternative alelle\n\n"
				+ "Config options: " + "	-cf+suffix <value> [none]\n\n";
// 		General Information (END)
// 		Specific Information (START)
		String fragsinfo = "\nUsage: Leucippus frag -o ofile [options] [ file1.fq.gz file2.fq.gz ]\n\n"
				+ "Options: -o  FILE        output file in gzip (e.g., longreads.fq.gz)\n"
                                + "         -o1 FILE        output file in gzip (e.g., remainedreads1.fq.gz)\n"
                                + "         -o2 FILE        output file in gzip (e.g., remainedreads2.fq.gz)\n"
                                + "         -l              maximum long reads length [500]\n"
                                + "         -minov          minimum overlapping region [50]\n"
				+ "Input:   file1.fq.gz     input gzipped file with reads in fastq format\n"
				+ "                         (if not provided, input is taken from standard input)\n"
				+ "         file2.fq.gz     input gzipped file with reads in fastq format\n";

		String tableinfo = "\nUsage: Leucippus noisetab [options] file1.bam [file2.bam, ...]\n\n"
				+ "Options:     -interval FILE  file with intervals\n"
				+ "             -ref      FILE  FASTA file with reference sequence (could be .gz file)\n"
				+ "             -o        FILE  output file\n"
				+ "             -q              base quality cut off 0-100 [20]\n"
				+ "             -pad            padding range at read ends [1]\n"
				+ "             -d              distance cut off [-1]\n\n";

		String graphinfo = "\nUsage: Leucippus graph [options] -o <prefix> table1.file [table2.file]\n\n"
				+ "Options:    -type            graph type: pvalue|mutrate [pvalue]\n"
				+ "	    -coverage INT     minimum site/position coverage [100]\n"
				+ "            -range    DOUBLE  maximum range for error [0.05]\n"
				+ "            -overlap  INT     use positions in overlapping 3'-ends of reads\n"
				+ "                                (the number specifies read length and it is\n"
				+ "                                 only useful for analysis of amplicon-seq data)\n\n"
//				+ "            -fpval           generate frequency p-value table\n"
				+ "            -o               prefix for output files: prefix.pdf\n\n";
//				+ "                                                      prefix.fpvtb[1,2].tsv\n\n";

		String decideinfo = "\nUsage: Leucippus decide [options] table.file\n\n"
				+ "Options: -o           FILE    results\n"
				+ "         -coverage            minimum site/position coverage [100]\n"
				+ "         -germlineAF  DOUBLE  minimum AF to call variant as germline [0.35]\n"
				+ "         -pvalue      DOUBLE  p-value for somatic call[0.05]\n\n";

		String extractinfo = "\nUsage: Leucippus extract [options] bam.file\n\n"
				+ "Options: -o           FILE   results\n"
				+ "         -a                  alelle_X [X could be A,C,T,G]\n"
				+ "         -p                  chr:position\n\n";
		
		String confinfo = "\nUsage: Not specified. Instead Options are used.\n\n"
				+ "Options:   -cfsam  <path to samtools executable [samtools]>\n"
				+ "           -cfR    <path to R executable [R]>\n"
//				+ "  -cfbwa  <path to bwa (none)>\n"
				+ "           -cftmp  <path to temporary folder [./Temp]>\n"
				+ "           -cfbash <path to bash shell [/bin/bash]>\n"
				+ "           -cfseparator character\n\n";
//				+ "  -cfref  <path to ref (none)>\n\n";
// 		Specific Information(END)
// 		Information printed on screen (END)

// 		final argument arrays (START)
		String[][] OverlapArgs = new String[8][2]; // frag
		String[][] TablesArgs =  new String[7][2]; // noisetab
		String[][] GraphsArgs =  new String[6][2]; // graph
		String[][] DecideArgs =  new String[5][2]; // decide
		String[][] ExtractArgs = new String[4][2]; // extract
// 		final argument arrays (END)

// 		Leucippus States(States 0 = no state, 1-4 = task state
		int state = 0;
		boolean doubleoccur = false;
// 		Identify state (START)
// 		For each command (frag, noisetab, graph, decide, and extract) there is a
// 			corresponding state : 1, 2, 3, 4, 5
		for (int i = 0; i < args.length; i++)
		{
			if ((args[i].equals("frag")) && ((state == 0) || (state == 1)))
				state = 1;
			if ((args[i].equals("noisetab")) && ((state == 0) || (state == 2)))
				state = 2;
			if ((args[i].equals("graph")) && ((state == 0) || (state == 3)))
				state = 3;
			if ((args[i].equals("decide")) && ((state == 0) || (state == 4)))
				state = 4;
			if ((args[i].equals("extract")) && ((state == 0) || (state == 5)))
				state = 5;
			if ((args[i].equals("frag"))
					&& ((state == 2) || (state == 3) || (state == 4) || (state == 5)))
				doubleoccur = true;
			if ((args[i].equals("noisetab"))
					&& ((state == 1) || (state == 3) || (state == 4) || (state == 5)))
				doubleoccur = true;
			if ((args[i].equals("graph"))
					&& ((state == 1) || (state == 2) || (state == 4) || (state == 5)))
				doubleoccur = true;
			if ((args[i].equals("decide"))
					&& ((state == 1) || (state == 2) || (state == 3) || (state == 5)))
				doubleoccur = true;
			if ((args[i].equals("extract"))
					&& ((state == 1) || (state == 2) || (state == 3) || (state == 4)))
				doubleoccur = true;
		}
// 		Identify state (END)
		if (doubleoccur == true) {
			System.out.println("The main command must be 'frag' or 'noisetab' " +
					"or 'graph' or 'decide' or 'extract' (only one of them).");
			Exit();
		}
// 		End States and one instance restriction

		if (args.length == 0)
		{
			System.out.println(leu_gen_nfo);
			Exit();
		}

		if (args.length == 1)
			if ((args[0].equals("-h")) || (args[0].equals("-help"))
					|| (args[0].equals("help"))) {
				System.out.println(leu_gen_nfo);
				Exit();
			}

// 		Provide information and Exit(START) if "args" length is <3
		if (args.length < 3) {
			for (int i = 0; i < args.length; i++) {
				if ((args[i].equals("frag")) || (args[i].equals("-frag")))

				{
					System.out.println(fragsinfo);
					Exit();
				}

				if ((args[i].equals("graph")) || (args[i].equals("-graph"))) {
					System.out.println(graphinfo);
					Exit();
				}

				if ((args[i].equals("noisetab"))
						|| (args[i].equals("-noisetab"))) {
					System.out.println(tableinfo);
					Exit();
				}

				if ((args[i].equals("decide")) || (args[i].equals("-decide"))) {
					System.out.println(decideinfo);
					Exit();
				}
				if ((args[i].equals("extract")) || (args[i].equals("-extract"))) {
					System.out.println(extractinfo);
					Exit();
				}

				if ((args[i].equals("config"))
						|| (args[i].equals("configuration"))
						|| (args[i].equals("-config"))
						|| (args[i].equals("-configuration"))) {
					System.out.println(confinfo);
					Exit();
				}
			}
		}

// 		Provide information and Exit (END)

// 		Declare Variables (START)

		String // OS = "",
				// slash = "/",
		root_pth = "", fast1 = "", fast2 = "", fast1_cm_fast2 = "", 
		fragout = "", shortread1out="", shortread2out="", shrts1c2out="", grorpvtb = "", 
		lrleng = "", separator="", gref_pth = "", pad="", 
		intable = "", outdecide = "", coverages = "", germlineAFs = "", pvalues="", 
		prmrs_pth = "", fragpipenp = "", minvrlps="", tables_path = "", 
		graphs_path = "", graph_type = "", graphs_overlap = "", dcut = "",
		outextract="", bamsextract="", positionextract="", alellextract="";
		double pvalued=0.0;
		int k = 1, p = 0, dcuti = 500;
		int lrlengi = 0, minvrlpi=0;
//		 Declare Variables(END)
		
		String cfpath = System.getProperty("user.dir");
		String[][] cfargs = ArgsConfigCheck(args);
		String bsh_pth = cfargs[0][1];
		String smtls_pth = cfargs[1][1];
		// String bwa_pth = cfargs[2][1];
		// String ref_pth = cfargs[3][1];
		String rscrp_path = cfargs[4][1];
		String tmppath = cfargs[5][1];
		separator = cfargs[6][1];
		for (int i = 0; i < 6; i++)
			System.out.println(cfargs[i][0] + " " + cfargs[i][1]);

		Vector<String> bms = new Vector<String>();
		Vector<String> prrfsqv = new Vector<String>();

// 		STATES (START)
// 		correspondence states : frag 1, noisetab 2, graph 3, decide 4, extract 5
		if (state != 0) {

			// Four "if" statements for the four states
			// Each statement has a test "args" method call at the beginning
			// If the arguments have errors then the testing method will
			// provide an informative error message and then it will terminate
			// the program.
			// If the arguments are correct then the testing method will
			// populate
			// the corresponding two dimensional array (it holds modified
			// arguments) and will returned here.
			// Next the if statement will generate the final command method
			// call. Once the final command method finishes its task 
			// the main method terminates the program.

// 			frag (START) state 1
			if (state == 1) {
				System.out.println("Make Long Reads!  " + args.length);

				BufferedReader br1 = null, br2 = null;
				BufferedWriter bw = null;
				GZIPOutputStream gzipout = null;
				// gzipoutr1 = null, gzipoutr2 = null;

				GZIPInputStream gzip1 = null, gzip2 = null;
				
				OverlapArgs = testArgsforFrags(args);

				for (int po = 0; po < OverlapArgs.length; po++)
					System.out.println(OverlapArgs[po][0] + " : "
							+ OverlapArgs[po][1]);

				if (OverlapArgs[2][0].equals("bamfile"))
					OverlapArgs[2][1] = "piped input bam file expeced.";

				fragpipenp = OverlapArgs[0][1];
				fast1_cm_fast2 = OverlapArgs[1][1];
				fragout = OverlapArgs[2][1]; // output path
				shortread1out = OverlapArgs[3][1];
				shortread2out = OverlapArgs[4][1];
				if( !(shortread1out==null) && !(shortread2out==null) )
					shrts1c2out=shortread1out + "," + shortread2out;

				lrleng = OverlapArgs[5][1]; // long read length
				lrlengi = Integer.parseInt(lrleng);
				minvrlps = OverlapArgs[6][1];
				minvrlpi = Integer.parseInt(minvrlps);

				if (OverlapArgs[0][1].equals("pipe")) {
					System.out.println("main method 2 'pipe' no input file.\n"
							+ "Proceed with Make Long Reads!");
					gzipout = new GZIPOutputStream(
							new FileOutputStream(fragout));
					bw = new BufferedWriter(new OutputStreamWriter(gzipout));
					// if then .., else null(START)
					// gzipoutr1 = new GZIPOutputStream(
					//		new FileOutputStream(shortread1out));
//					bw1 = new BufferedWriter(new OutputStreamWriter(gzipoutr1));
//
//					gzipoutr2 = new GZIPOutputStream(
//							new FileOutputStream(shortread2out));
//					bw2 = new BufferedWriter(new OutputStreamWriter(gzipoutr2));
					// if then .., else null(END)
					br1 = new BufferedReader(new InputStreamReader(System.in));
					br2 = null;
					MakeLongReads(bw, shrts1c2out, br1, br2, lrlengi, p, minvrlpi, separator);
				}

				if ((OverlapArgs[0][1].equals("notpipe_gz"))
						|| (OverlapArgs[0][1].equals("notpipe"))) {
					if (OverlapArgs[1][0].equals("InTwofastq")) {
						System.out
								.println("main method 2 'not pipe' two input " +
								"fastq.gz files.\n"
										+ "Proceed with Make Long Reads!");
						fast1 = "";
						fast2 = "";
						String[] wrdes = fast1_cm_fast2.split(",");
						fast1 = wrdes[0];
						fast2 = wrdes[1];
						if (!(fast1.isEmpty())) {
							gzip1 = new GZIPInputStream(new FileInputStream(
									fast1));
							br1 = new BufferedReader(new InputStreamReader(
									gzip1));
						}
						if (!(fast2.isEmpty())) {
							gzip2 = new GZIPInputStream(new FileInputStream(
									fast2));
							br2 = new BufferedReader(new InputStreamReader(
									gzip2));
						}
						gzipout = new GZIPOutputStream(new FileOutputStream(
								fragout));
						bw = new BufferedWriter(new OutputStreamWriter(gzipout));
						// if then .., else null(START)
//						System.out.println(shortread1out + "  " + shortread2out);

						// if then .., else null(END)
						// MakeLongReads(bw, br1, br2, lrlengi, p);
						MakeLongReads(bw, shrts1c2out, br1, br2, lrlengi, p, minvrlpi, separator);

					}
					if (OverlapArgs[1][0].equals("InOnefastq")) {
						System.out
								.println("main method 2 'not pipe' one input " +
								"fastq.gz file.\n"
										+ "Proceed with Make Long Reads!");
						fast1 = fast1_cm_fast2;
						fast2 = "";
						if (!(fast1.isEmpty())) {
							gzip1 = new GZIPInputStream(new FileInputStream(
									fast1));
							br1 = new BufferedReader(new InputStreamReader(
									gzip1));
						}
						if (!(fast2.isEmpty())) {
							gzip2 = new GZIPInputStream(new FileInputStream(
									fast2));
							br2 = new BufferedReader(new InputStreamReader(
									gzip2));
						}
						gzipout = new GZIPOutputStream(new FileOutputStream(
								fragout));
						bw = new BufferedWriter(new OutputStreamWriter(gzipout));


//						gzipoutr1 = new GZIPOutputStream(
//							new FileOutputStream(shortread1out));
//						bw1 = new BufferedWriter(new OutputStreamWriter(gzipoutr1));
//						gzipoutr2 = new GZIPOutputStream(
//							new FileOutputStream(shortread2out));
//						bw2 = new BufferedWriter(new OutputStreamWriter(gzipoutr2));

						// if then .., else null(END)

						// MakeLongReads(bw, br1, br2, lrlengi, p);
						MakeLongReads(bw, shrts1c2out, br1, br2, lrlengi, p, minvrlpi, separator);
					}
				}
			}
// 			frag (END)

// 			noisetab (START) state 2
			if (state == 2) {
				System.out.println("Make Tables! " + args.length);
				TablesArgs = testArgsforTables(args);

				for (int o1 = 0; o1 < TablesArgs.length; o1++) {
					System.out.println(TablesArgs[o1][0] + "  "
							+ TablesArgs[o1][1]);
				}

				System.out.println();
				String perc = TablesArgs[3][1];
				String input = "";
				String output = TablesArgs[2][1];
				dcut = TablesArgs[4][1];
				dcuti = Integer.parseInt(dcut);
				prmrs_pth = TablesArgs[1][1];
				gref_pth = TablesArgs[5][1];
				pad=TablesArgs[6][1];
				// prrfsqv = GenomeRefParser(prmrs_pth, gref_pth);
				System.out.println("Request Genome reference sequences!");
				prrfsqv = IndependentGenomeRefParser(prmrs_pth, gref_pth);
				if (prrfsqv.size() == 0) {
					System.out
							.println("Reference parser returned zero referense"+
							"sequences\nExit.");
					Exit();
				}

				if (TablesArgs[0][0].equals("bams")) {
					String[] pwrds = TablesArgs[0][1].split(",");
					for (int bm = 0; bm < pwrds.length; bm++)
						bms.add(pwrds[bm]);
					input = "none";
				}

				output = TablesArgs[2][1];
				root_pth = new java.io.File(".").getCanonicalPath();
				root_pth = root_pth + "/";
				System.out.println("root path = " + root_pth);

				for (int jk = 0; jk < bms.size(); jk++)
					System.out.println(bms.get(jk));
				/*
				 * TablesArgs[0][1]="bam,bam,...";
				 * TablesArgs[1][1]="interval"; TablesArgs[2][1]="output";
				 * TablesArgs[3][1]="quality"; TablesArgs[4][1]="distance";
				 * TablesArgs[5][1]="reference"; TablesArgs[6][1]="pad";
				 */
				System.out.println("main method 2 Proceed with Table!");
				makeTables(bsh_pth, smtls_pth, input, output, perc, prrfsqv,
						bms, tmppath, dcuti, pad);
				k = 0;
			}
// 			noisetab (END)
			
// 			graph (START) state 3
			if (state == 3) {
				// type, coverage, range, overlap, fpval(optional) o input
				// 1 5 6 4 3 2
				System.out.println("Make Graphs! " + args.length);
				// Thread.sleep(1000);
				GraphsArgs = testArgsforGraphs(args);
				graph_type = GraphsArgs[0][1];
				tables_path = GraphsArgs[1][1];
				graphs_path = GraphsArgs[2][1];
				graphs_overlap = GraphsArgs[3][1];
				coverages = GraphsArgs[4][1];
				germlineAFs = GraphsArgs[5][1];

				System.out.println("Graph Type " + GraphsArgs[0][0] + "  "
						+ GraphsArgs[0][1]);
				System.out.println("Input Table(s) " + GraphsArgs[1][0] + "  "
						+ GraphsArgs[1][1]);
				System.out.println("Output " + GraphsArgs[2][0] + "  "
						+ GraphsArgs[2][1]);
				System.out.println("Overlap " + GraphsArgs[3][0] + "  "
						+ GraphsArgs[3][1]);
				System.out.println("Coverage " + GraphsArgs[4][0] + "  "
						+ GraphsArgs[4][1]);
				System.out.println("Range " + GraphsArgs[5][0] + "  "
						+ (1 - Double.parseDouble(GraphsArgs[5][1])));

				if(rscrp_path.isEmpty())
				{
					System.out.println("R script executable has not been provided.\nExit.");
					Exit();
				}

				// added to provide option for creation of frequency-p-value
				// table
				grorpvtb = GraphsArgs[2][0];

				System.out.println("graph or pvalue table : " + grorpvtb);
				System.out.println("main method 2 Proceed with Graph!\n\n");
				CreateGraphs(graph_type, tables_path, graphs_path, rscrp_path,
						bsh_pth, cfpath, grorpvtb, graphs_overlap, coverages,
						germlineAFs);
			}
// 			graph(END)

// 			decide (START) state 4
			if (state == 4) {
				System.out.println("Decide!");
				DecideArgs = testArgsforDecide(args);
				intable = DecideArgs[0][1];
				outdecide = DecideArgs[1][1];
				coverages = DecideArgs[2][1];
				germlineAFs = DecideArgs[3][1];
				pvalues = DecideArgs[4][1];
				pvalued = Double.parseDouble(pvalues);
				System.out.println("input table : " + intable);
				System.out.println("output results : " + outdecide);
				System.out.println(DecideArgs[2][0] + " : " + coverages);
				System.out.println(DecideArgs[3][0] + " : " + germlineAFs);
				System.out.println(DecideArgs[4][0] + " : " + pvalued + "\n");
				System.out.println("main method 2 Proceed with Decide!");
				Decide(intable, outdecide, coverages, germlineAFs, pvalued);
			}
// 			decide(END)

// 			extract (START) state 5
			if (state == 5) {
				System.out.println("Extract!");
				ExtractArgs = testArgsforExtract(args);
				bamsextract=ExtractArgs[0][1];
				if (ExtractArgs[0][0].equals("bams")) {
					String[] pwrds = bamsextract.split(",");
					for (int bm = 0; bm < pwrds.length; bm++)
						bms.add(pwrds[bm]);
					//input = "none";
				}
				else if (ExtractArgs[0][0].equals("nobams"))
				{
					System.out.println("No bams to proceed!");
					Exit();
				}
				//outextract="", bamsextract="", positionextract="", alellextract="";
				outextract = ExtractArgs[3][1];
				root_pth = new java.io.File(".").getCanonicalPath();
				root_pth = root_pth + "/";
				System.out.println("root path = " + root_pth);

				//for (int jk = 0; jk < bms.size(); jk++)
					//System.out.println(bms.get(jk));
				alellextract = ExtractArgs[1][1];
				positionextract = ExtractArgs[2][1];
				
				System.out.println("bam files : " + bamsextract);
				System.out.println("output results : " + outextract);
				System.out.println("positionextract : " + positionextract);
				System.out.println("alelextract : " + alellextract);
				
				System.out.println("main method 2 Proceed with Extract!");
				// Extract(intable, outdecide, coverages, germlineAFs, pvalued);
				Extract(bsh_pth, smtls_pth, outextract, positionextract, alellextract,
						bms, tmppath);
			}
// 			decide(END)

		}
// 		STATES(END)
		System.out.println("End of Main method!");
	} 
// 	Main method(END)

	/**
	 * Method that checks for validity of user entered arguments. The method
	 * performs three tasks. It tests the paths/files if they exist. In returned
	 * two dimensional array usually it populates the [x][0] element with an
	 * informative value, and the [x][1] with the actual value. The method has
	 * no capability of termination.
	 * 
	 * @param argmns
	 * @return resargs String[][] a two dimensional String array is returned
	 * @throws IOException
	 *             (the method was designed to read a configuration file; now it
	 *             tests for path existence)
	 * 
	 */
	public static String[][] ArgsConfigCheck(String[] argmns)
			throws IOException {
	    // ArgsConfigCheck method (START)
		System.out.println("\nArgsConfigCheck method Start! ----------");
		String confinfo = "configuration [main command doesnt exist)]\n"
				+ "insread a simple prefix '-cf' followed by the suffix is used.\n"
				+ "-cfsam       [-cfsam]  [path for samtools]\n"
				+ "-cfR     [-cfR]    [path for R]\n"
				//+ "-cfbwa       [-cfbwa]  [path for bwa]\n"
				+ "-cftmp       [-cftmp]  [path for temporary folder]\n"
				+ "-cfbash      [-cfbash] [path for bash]\n";
				//+ "-cfref       [-cfref]  [path for ref]";

		String OS = "", root_pth = "", bsh_pth = "", smtls_pth = "", 
		bwa_pth = "", ref_bam_pth = "", rscrp_path = "", tmp_pth = "", inpth="", 
		slash = "/", separator="";
		int corcnt = 0;
		File fl;
		String jcur="";
		String[][] resargs = new String[7][2];
// Indices [x][0] are for info; [x][1] are for values
		
		resargs[0][0] = "bsh_pth";
		resargs[0][1] = "";
		resargs[1][0] = "smtls_pth";
		resargs[1][1] = "";		
		resargs[2][0] = "bwa_pth";
		resargs[2][1] = "";
		resargs[3][0] = "ref_bam_pth";
		resargs[3][1] = "";
		resargs[4][0] = "rscrp_path";
		resargs[4][1] = "";
		resargs[5][0] = "tmp_path";
		resargs[5][0] = "";
		resargs[6][0] = "separatordef";
		resargs[6][1] = "";
		for (int i = 0; i < 6; i++)
			resargs[i][1] = "";

		for (int i = 0; i < argmns.length; i++) 
		{
			if ((argmns[i].equals("-cfbash")) && (i < argmns.length - 1)) 
			{
				
				if (testFileExistence(argmns[i + 1]))
				{
					bsh_pth = argmns[i + 1];
					resargs[0][1] = bsh_pth;
					corcnt = corcnt + 1;
					System.out.println("bash executable exists!");
				}
				else
				{
					System.out.println("bash executable not found.\nExit.");
					Exit();
				}

//				if (testParentDirectoryExistence(argmns[i + 1])) {
//					bsh_pth = argmns[i + 1];
//					resargs[0][1] = bsh_pth;
//					corcnt = corcnt + 1;
//					System.out.println("bash path directory exists!");
//				} else {
//					System.out.println("bash path directory not found.");
//					resargs[0][1] = "";
//				}
			}
			if ((argmns[i].equals("-cfR")) && (i < argmns.length - 1)) 
			{
				if (testFileExistence(argmns[i + 1]))
				{
					rscrp_path = argmns[i + 1];
					resargs[4][1] = rscrp_path;
					corcnt = corcnt + 1;
					System.out.println("R executable exists!");
				}
				else
				{
					System.out.println("R executable not found.\nExit.");
					Exit();
				}

				//if (testParentDirectoryExistence(argmns[i + 1])) {
				//	rscrp_path = argmns[i + 1];
									
				//	resargs[4][1] = rscrp_path;
				//	corcnt = corcnt + 1;
				//	System.out.println("R path directory exists!");
				//} else {
			//		System.out.println("R path directory not found.");
			//		resargs[4][1] = "";
			//	}
			}



			if ((argmns[i].equals("-cfsam")) && (i < argmns.length - 1)) 
			{

				if (testFileExistence(argmns[i + 1]))
				{
					smtls_pth = argmns[i + 1];
					resargs[1][1] = smtls_pth;
					corcnt = corcnt + 1;
					System.out.println("samtools executable exists!");
				}
				else
				{
					System.out.println("samtools executable not found.\nExit.");
					Exit();
				}


				//if (testParentDirectoryExistence(argmns[i + 1])) {
				//	smtls_pth = argmns[i + 1];
				//	resargs[1][1] = smtls_pth;
				//	corcnt = corcnt + 1;
				//	System.out.println("samtools path directory exists!");
				//} else {
				//	System.out.println("samtools path directory not found.");
				//	resargs[1][1] = "";
				//}
			}

			if ((argmns[i].equals("-cfbwa")) && (i < argmns.length - 1)) {
				if (testParentDirectoryExistence(argmns[i + 1])) {
					bwa_pth = argmns[i + 1];
					resargs[2][1] = bwa_pth;
					corcnt = corcnt + 1;
					System.out.println("bwa path directory exists!");
				} else {
					System.out.println("bwa path directory not found.");
					resargs[2][1] = "";
				}
			}
			if ((argmns[i].equals("-cfseparator")) && (i < argmns.length - 1))
			{
					separator = argmns[i + 1];
					resargs[6][0] = "separator";
					resargs[6][1] = separator;
					corcnt = corcnt + 1;
					System.out.println("separator has been stated!");
			}


			if ((argmns[i].equals("-cftmp")) && (i < argmns.length - 1)) 
			{
				jcur= argmns[i + 1];
				fl=new File(jcur);
				inpth = fl.getName();
				if (jcur.equals(inpth)) 
				{
					root_pth = new java.io.File(".").getCanonicalPath();
					root_pth = root_pth + "/";
					jcur = root_pth + jcur;
					fl=new File(jcur);
					if(!fl.exists())
						fl.mkdir();
					if(!fl.exists()) 
					{
  						System.out.println("Temporary path directory cannot be created.\nExit.");
						Exit();
					}

					tmp_pth = jcur;
					resargs[5][1] = tmp_pth;
					corcnt = corcnt + 1;
				}
				else if (testParentDirectoryExistence(jcur)) 
				{
					fl=new File(jcur);
					if(!fl.exists())
						fl.mkdir();
					if(!fl.exists())
					{
  						System.out.println("Temporary path directory cannot be created.\nExit.");
						Exit();
					}
					else
					{					
						tmp_pth = jcur;
						resargs[5][1] = tmp_pth;
						corcnt = corcnt + 1;
						System.out.println("Temporary path directory exists!");
					}
				} 
				else if (!(testParentDirectoryExistence(jcur)))
				{
					System.out.println("Temporary parent path directory not found.\nExit.");
					Exit();
				} 
					//a" +
					//		" new 'Temp' directory will be created inside the " +
					//		"current directory.");
					// resargs[5][1] = "none";
					
					//root_pth = new java.io.File(".").getCanonicalPath();
					//root_pth = root_pth + "/Temp";
					//jcur = root_pth;
					//fl=new File(jcur);
					//if(!fl.exists())
					//	fl.mkdir();
					//if(!fl.exists())
					//{
  				//		System.out.println("Temporary path directory cannot be created.\nExit.");
				//		Exit();
				//	}
				//	else
				//	{
				//		resargs[5][1] = jcur;
						
				//	}
				//	}//
			}
		}
		if( (argmns[0].equals("noisetab")) || (argmns[0].equals("extract")) )
		{
			System.out.println("args noisetab ----" + argmns[0]);
			if (resargs[5][1].equals(""))
			{
				// tmp_pth = new java.io.File( "." ).getCanonicalPath();
				tmp_pth = System.getProperty("user.dir");
				tmp_pth = tmp_pth + slash + "Temp";
				fl = new File(tmp_pth);
				if (!(fl.exists()))
				fl.mkdir();
				// fl.createNewFile();
				resargs[5][1] = tmp_pth;
			}
		}
		if (resargs[4][1].equals("")) 
		{
			resargs[4][1] = "";
		}
		if (resargs[3][1].equals("")) {
			resargs[3][1] = "none";
		}
		if (resargs[2][1].equals("")) {
			resargs[2][1] = "";
		}
		if (resargs[1][1].equals("")) {
			resargs[1][1] = "none";
		}
		if (resargs[0][1].equals("")) {
			resargs[0][1] = "bin/bash";
		}

		System.out.println("-------ArgsConfigCheck method End.\n");
		return resargs;
	      // ArgsConfigCheck method (END)

	}


	/**
	 * Method that checks for validity of user entered arguments. The method
	 * performs three tasks. It tests the paths/files if they exist. In returned
	 * two dimensional array usually it populates the [x][0] element with an
	 * informative value, and the [x][1] with the actual value. The method has
	 * no capability of termination.
	 * 
	 * @param argmns
	 * @return resargs String[][] a two dimensional String array is returned
	 * @throws IOException
	 *             (the method was designed to read a configuration file; now it
	 *             tests for path existence)
	 * 
	 */
	public static String[][] ArgsConfigCheckBack2(String[] argmns)
			throws IOException {
	    // ArgsConfigCheck method (START)
		System.out.println("\nArgsConfigCheck method Start! ----------");
		String confinfo = "configuration [main command doesnt exist)]\n"
				+ "insread a simple prefix '-cf' followed by the suffix is used.\n"
				+ "-cfsam       [-cfsam]  [path for samtools]\n"
				+ "-cfR     [-cfR]    [path for R]\n"
				//+ "-cfbwa       [-cfbwa]  [path for bwa]\n"
				+ "-cftmp       [-cftmp]  [path for temporary folder]\n"
				+ "-cfbash      [-cfbash] [path for bash]\n";
				//+ "-cfref       [-cfref]  [path for ref]";

		String OS = "", root_pth = "", bsh_pth = "", smtls_pth = "", 
		bwa_pth = "", ref_bam_pth = "", rscrp_path = "", tmp_pth = "", 
		slash = "/";
		int corcnt = 0;
		File fl;
		String[][] resargs = new String[6][2];
// Indices [x][0] are for info; [x][1] are for values
		
		resargs[0][0] = "bsh_pth";
		resargs[1][0] = "smtls_pth";
		resargs[2][0] = "bwa_pth";
		resargs[3][0] = "ref_bam_pth";
		resargs[4][0] = "rscrp_path";
		resargs[5][0] = "tmp_path";
		for (int i = 0; i < 6; i++)
			resargs[i][1] = "";

		for (int i = 0; i < argmns.length; i++) {
			if ((argmns[i].equals("-cfbash")) && (i < argmns.length - 1)) {
				if (testParentDirectoryExistence(argmns[i + 1])) {
					bsh_pth = argmns[i + 1];
					resargs[0][1] = bsh_pth;
					corcnt = corcnt + 1;
					System.out.println("bash path directory exists!");
				} else {
					System.out.println("bash path directory not found.");
					resargs[0][1] = "";
				}
			}

			if ((argmns[i].equals("-cfsam")) && (i < argmns.length - 1)) {
				if (testParentDirectoryExistence(argmns[i + 1])) {
					smtls_pth = argmns[i + 1];
					resargs[1][1] = smtls_pth;
					corcnt = corcnt + 1;
					System.out.println("samtools path directory exists!");
				} else {
					System.out.println("samtools path directory not found.");
					resargs[1][1] = "";
				}
			}

			if ((argmns[i].equals("-cfbwa")) && (i < argmns.length - 1)) {
				if (testParentDirectoryExistence(argmns[i + 1])) {
					bwa_pth = argmns[i + 1];
					resargs[2][1] = bwa_pth;
					corcnt = corcnt + 1;
					System.out.println("bwa path directory exists!");
				} else {
					System.out.println("bwa path directory not found.");
					resargs[2][1] = "";
				}
			}

			if ((argmns[i].equals("-cfR")) && (i < argmns.length - 1)) {
				if (testParentDirectoryExistence(argmns[i + 1])) {
					rscrp_path = argmns[i + 1];
					resargs[4][1] = rscrp_path;
					corcnt = corcnt + 1;
					System.out.println("R path directory exists!");
				} else {
					System.out.println("R path directory not found.");
					resargs[4][1] = "";
				}
			}

			if ((argmns[i].equals("-cftmp")) && (i < argmns.length - 1)) {
				if (testParentDirectoryExistence(argmns[i + 1])) {
					tmp_pth = argmns[i + 1];
					resargs[5][1] = tmp_pth;
					corcnt = corcnt + 1;
					System.out.println("Temporary path directory exists!");
				} else {
					System.out
							.println("Temporary path directory not found.\n a" +
							" new 'Temp' directory will be created inside the " +
							"current directory.");
					resargs[5][1] = "none";
				}
			}
		}
		if(argmns[0].equals("noisetab"))
		{
			System.out.println("args noisetab ----" + argmns[0]);
			if (resargs[5][1].equals(""))
			{
				// tmp_pth = new java.io.File( "." ).getCanonicalPath();
				tmp_pth = System.getProperty("user.dir");
				tmp_pth = tmp_pth + slash + "Temp";
				fl = new File(tmp_pth);
				if (!(fl.exists()))
				fl.mkdir();
				// fl.createNewFile();
				resargs[5][1] = tmp_pth;
			}
		}
		if (resargs[4][1].equals("")) {
				resargs[4][1] = "";
		}
		if (resargs[3][1].equals("")) {
			resargs[4][1] = "none";
		}
		if (resargs[2][1].equals("")) {
			resargs[4][1] = "";
		}
		if (resargs[1][1].equals("")) {
			resargs[4][1] = "";
		}
		if (resargs[0][1].equals("")) {
			resargs[4][1] = "bin/bash";
		}

		System.out.println("-------ArgsConfigCheck method End.\n");
		return resargs;
	      // ArgsConfigCheck method (END)

	}

	// test frag(START)
	/**
	 * Method that checks for validity of user entered arguments. The method
	 * performs three tasks. It tests the paths/files if they exist. In returned
	 * two dimensional array usually it populates the [x][0] element with an
	 * informative value, and the [x][1] with the actual value. The method
	 * terminates the program if it finds errors.
	 * 
	 * @param argmns
	 *            : Strin[]
	 * @return resargs String[][] a two dimensional String array is returned
	 * @throws IOException
	 *             (the method was designed to read a configuration file; now it
	 *             tests for path existence)
	 * @param argmns
	 * @return
	 * @throws IOException
	 */
	public static String[][] testArgsforFragsBackres(String[] argmns)
			throws IOException {

		System.out.println("Start of testArgsforFrags method!");
		String fragsinfo = "\nUsage: Leucippus frag -o ofile [options] [ file1.fq.gz file2.fq.gz ]\n\n"
				+ "Options: -o  FILE        output file in gzip (e.g., longreads.fq.gz)\n"
                                + "         -o1 FILE        output file in gzip (e.g., remainedreads1.fq.gz)\n"
                                + "         -o2 FILE        output file in gzip (e.g., remainedreads2.fq.gz)\n"
                                + "         -l              maximum long reads length [500]\n"
                                + "         -minov          minimum overlapping region [50]\n"
				+ "Input:   file1.fq.gz     input gzipped file with reads in fastq format\n"
				+ "                         (if not provided, input is taken from standard input)\n"
				+ "         file2.fq.gz     input gzipped file with reads in fastq format\n";


		boolean missinparrd12out = false;
		int corcnt = 0;
		String root_pth = "", inpth = "", flpth = "";

		String maxlrleg = "";
		int maxlrlegi = 0;
		double db = 0.0;

		String ftnm = "";
		String correct_report = "", missed_report = "";
		String cucrinf = "", cumsinf = "";
		File outfl, inpfl;
		String lread = "", rread = "", longread = "", shortread1="", shortread2="";
		// test for pipe
		int pibal = 0;
		String minovs = "";
		int minovi=0;

		String[][] resargs = new String[7][2];
		resargs[0][0] = "pipeornot";
		resargs[1][0] = "input";
		resargs[2][0] = "output";
		resargs[3][0] = "outputread1path";
		resargs[4][0] = "outputread2path";
		resargs[5][0] = "longreadmaxlength";
		resargs[5][1] = "500";
		resargs[6][0] = "minoveralpdef";
		resargs[6][1] = "50";
		// default value for long max length = 500
		Vector<String> infiles = new Vector<String>();
		String cur = "", jcur = "", inflss = "";

		String drct = "", testvalue = "";
		if ((argmns[argmns.length - 2].equals("-o"))
				|| (argmns[argmns.length - 2].equals("-o1"))
				|| (argmns[argmns.length - 2].equals("-o2"))
				|| (argmns[argmns.length - 2].equals("-l"))
				|| (argmns[argmns.length - 2].equals("-minov")))
			drct = "pipe";
		else
			drct = "not_pipe";
		System.out.println(drct);
// 		Test for input file(s) existence.
		if (drct.equals("not_pipe")) {
			for (int j = argmns.length - 1; j > 0; j--) {
				if (!(argmns[j - 1].charAt(0) == '-')) {
					jcur = argmns[j];
					inpfl = new File(jcur);
					inpth = inpfl.getName();
					System.out.println(inpth);
					if (jcur.equals(inpth)) {
						root_pth = new java.io.File(".").getCanonicalPath();
						root_pth = root_pth + "/";
						jcur = root_pth + jcur;
					}
					if (testFileExistence(jcur)) {
						infiles.add(jcur);
						System.out.println("input file exists : " + jcur);
					} else {
						System.out.println("input file not found : " + jcur
								+ "\nExit.");
						Exit();
					}
				} else
					j = 0;
			}
		}

// 		Test(END) for input file(s) existence.
		if (infiles.size() == 1) {
			System.out.println("One fastq file as input!");
			resargs[1][0] = "InOnefastq";
		}
		if (infiles.size() == 2) {
			System.out.println("Two fastq files as input!");
			resargs[1][0] = "InTwofastq";
		}
		if (infiles.size() > 2) {
			System.out
					.println("More than two valid files as an input found. Exit.");
			Exit();
		}
		String testq = "";
		int infcnt = 0;
		for (int u = 0; u < infiles.size(); u++) {
			infcnt = 0;
			testq = infiles.get(u);
			for (int uj = 0; uj < infiles.size(); uj++) {
				if (testq.equals(infiles.get(uj)))
					infcnt = infcnt + 1;
				if (infcnt > 1) {
					System.out
							.println("input file was typed two times in " +
							"arguments.\nExit.");
					Exit();
				}
			}
			// System.out.println(infiles.get(u));
		}

		for (int i = 0; i < argmns.length; i++) {
// 		Test for outputs file root directory existence.
			if (argmns[i].equals("-o")) {
				if (i < argmns.length - 1) {
					longread = argmns[i + 1];
					outfl = new File(longread);
					ftnm = outfl.getName();
					if (longread.equals(ftnm)) {
						root_pth = new java.io.File(".").getCanonicalPath();
						root_pth = root_pth + "/";
						flpth = root_pth + ftnm;
						resargs[2][0] = "outputpath";
						resargs[2][1] = flpth;
						corcnt = corcnt + 1;
						cucrinf = "Output is OK!";
						correct_report = correct_report + cucrinf + "\n";
						System.out.println("Output is OK!");
						//if (pibal == 1)
						//	resargs[3][1] = flpth;

					} else if (!(longread.equals(ftnm))) {
						if (testParentDirectoryExistence(argmns[i + 1])) {
							resargs[2][0] = "outputpath";
							longread = argmns[i + 1];
							resargs[2][1] = longread;
							corcnt = corcnt + 1;
							cucrinf = "Output is OK!";
							correct_report = correct_report + cucrinf + "\n";
							System.out.println("Output is OK!");
						//	if (pibal == 1)
						//		resargs[3][1] = longread;

						} else {
							cumsinf = "Output root folder not found.";
				
							missed_report = missed_report + cumsinf + "\n";
							// System.out.println("Output root folder not found.");
						}
					}
				}
			}

			if (argmns[i].equals("-o1")) {
				if (i < argmns.length - 1) {
					shortread1 = argmns[i + 1];
					outfl = new File(shortread1);
					ftnm = outfl.getName();
					if (shortread1.equals(ftnm)) {
						root_pth = new java.io.File(".").getCanonicalPath();
						root_pth = root_pth + "/";
						flpth = root_pth + ftnm;
						resargs[3][0] = "outputRead1Path";
						resargs[3][1] = flpth;
						corcnt = corcnt + 1;
						cucrinf = "OutputRead1Path is OK!";
						correct_report = correct_report + cucrinf + "\n";
						System.out.println(cucrinf);
					//	if (pibal == 1)
					//		resargs[3][1] = flpth;

					} else if (!(shortread1.equals(ftnm))) {
						if (testParentDirectoryExistence(argmns[i + 1])) {
							resargs[3][0] = "outputRead1Path";
							shortread1 = argmns[i + 1];
							resargs[3][1] = shortread1;
							corcnt = corcnt + 1;
							cucrinf = "OutputRead1Path is OK!";
							correct_report = correct_report + cucrinf + "\n";
							System.out.println(cucrinf);
					//		if (pibal == 1)
					//			resargs[3][1] = shoread1;

						} else {
							cumsinf = "Output read1-root folder not found.";
							missed_report = missed_report + cumsinf + "\n";
							missinparrd12out=true;
						// System.out.println("Output root folder not found.");
						}
					}
				}
			}

			if (argmns[i].equals("-o2")) {
				if (i < argmns.length - 1) {
					shortread2 = argmns[i + 1];
					outfl = new File(shortread2);
					ftnm = outfl.getName();
					if (shortread2.equals(ftnm)) {
						root_pth = new java.io.File(".").getCanonicalPath();
						root_pth = root_pth + "/";
						flpth = root_pth + ftnm;
						resargs[4][0] = "outputRead2Path";
						resargs[4][1] = flpth;
						corcnt = corcnt + 1;
						cucrinf = "OutputRead2Path is OK!";
						correct_report = correct_report + cucrinf + "\n";
						System.out.println(cucrinf);
					//	if (pibal == 1)
					//		resargs[3][1] = flpth;

					} else if (!(shortread2.equals(ftnm))) {
						if (testParentDirectoryExistence(argmns[i + 1])) {
							resargs[4][0] = "outputRead2Path";
							shortread2 = argmns[i + 1];
							resargs[4][1] = shortread2;
							corcnt = corcnt + 1;
							cucrinf = "OutputRead2Path is OK!";
							correct_report = correct_report + cucrinf + "\n";
							System.out.println(cucrinf);
					//		if (pibal == 1)
					//			resargs[3][1] = shoread2;

						} else {
							cumsinf = "Output read2-root folder not found.";
							missed_report = missed_report + cumsinf + "\n";
							missinparrd12out=true;
						// System.out.println("Output root folder not found.");
						}
					}
				}
			}

			if (argmns[i].equals("-l")) {
				if (i < argmns.length - 1) {
					maxlrleg = argmns[i + 1];
					if (isNumeric(maxlrleg)) {
						db = Double.parseDouble(maxlrleg);
						db = Math.ceil(db);
						maxlrlegi = (int) db;
						if (maxlrlegi > 3) {
							cumsinf = "Maximum long read length OK!";
							resargs[5][1] = Integer.toString(maxlrlegi);
							corcnt = corcnt + 1;
						} else {
							System.out
									.println("Maximum long read length too " +
									"small.");
							Exit();
						}
					} else {
						System.out
								.println("Maximum long read length Not " + 
								"Numeric.");
						Exit();
					}

				} else {
					System.out.println("Maximum long read length is missing.");
					Exit();
				}
			}

			if (argmns[i].equals("-minov")) {
				if (i < argmns.length - 1) {
					minovs = argmns[i + 1];
					if (isNumeric(minovs)) {
						db = Double.parseDouble(minovs);
						db = Math.ceil(db);
						minovi = (int) db;
						if (minovi > 2) {
							cumsinf = "Minimum overlapping region is OK!";
							resargs[6][0] = "minoverlap";
							resargs[6][1] = Integer.toString(minovi);
							corcnt = corcnt + 1;
						} else {
							System.out
									.println("Minimum overlapping region is too " +
									"small.\nExit.");
							Exit();
						}
					} else {
						System.out
								.println("Minimum overlapping region is Not " + 
								"Numeric.\nExit.");
						Exit();
					}

				} else {
					System.out.println("Minimum overlapping region error.\nExit.");
					Exit();
				}
			}
		}

		System.out.println(missed_report);

//		if( !( resargs[3][1].isEmpty()) && !( resargs[4][1].isEmpty()) )
		if( ( (resargs[3][1]==null) && (resargs[4][1]!=null) ) || ( (resargs[3][1]!=null) && (resargs[4][1]==null) ) )
		{		
			System.out.println("Both output for failed to overlap reads are required\n when optional '-o1' and '-o2' are used.\nOne output for failed to overlap reads is missing\nExit.");
			Exit();
		}


		if( (resargs[3][1]!=null) && (resargs[4][1]!=null) )		
		{
			if( resargs[3][1].equals(resargs[4][1]))
			{
				System.out.println("outputs for failed to overlap reads are same.\nExit.");
				Exit();
			}
			if( ( resargs[2][1].equals(resargs[3][1])) || ( resargs[2][1].equals(resargs[4][1])) )
			{
				System.out.println("outputs for Long reads and failed to overlap reads are same.\nExit.");
				Exit();
			}
		}
		//for (int i=2; i<5; i++)
		//	System.out.println(resargs[i][1]);
		// System.out.println("before exit");
		if(missinparrd12out==true)
			Exit();
		// for(int i=0; i<infiles.size(); i++)
		// {
		// if(i<infiles.size()-1)
		// inflss=inflss + infiles.get(i) + ",";
		// if(i==infiles.size()-1)
		// inflss=inflss + infiles.get(i);
		// }

		if (infiles.size() > 1) {
			for (int i = infiles.size() - 1; i > -1; i--) {
				if (i > 0)
					inflss = inflss + infiles.get(i) + ",";
				if (i == 0)
					inflss = inflss + infiles.get(i);
			}
		}

		// Error fixed on 05/14/2015
		// if only one fastq file was as an input the path was missed
		// fixed with the following if statement
		if (infiles.size() == 1)
			inflss = inflss + infiles.get(0);

		File tszpd;
		if (infiles.size() > 0) {
			tszpd = new File(infiles.get(0));
			if (isFileGZipped(tszpd)) {
				// resargs[1][0] = "leftrightreaddfqgz";
				resargs[0][1] = "notpipe_gz";
			}
			if (!(isFileGZipped(tszpd))) {
				// resargs[1][0] = "leftrightreaddfqgz";
				resargs[0][1] = "notpipe";
			}
		}
		// Test if the input file is a bam file that contains left and right
		// read (needs implementation)
		if (infiles.size() < 1) {
			resargs[1][0] = "bamfile";
			resargs[0][1] = "pipe";
			System.out.println("Piped bam file expected.");
			// resargs[1][1] = "none";
		}

		if (infiles.size() > 0) // string >0
		{
			resargs[1][1] = inflss;
		}
		System.out.println("End of testArgsforFrags method!");

		return resargs;
	}
	//test frag (END)

	public static String[][] testArgsforFrags(String[] argmns)
			throws IOException {

		System.out.println("Start of testArgsforFrags method!");
		String fragsinfo = "\nUsage: Leucippus frag -o ofile [options] [ file1.fq.gz file2.fq.gz ]\n\n"
				+ "Options: -o  FILE        output file in gzip (e.g., longreads.fq.gz)\n"
                                + "         -o1 FILE        output file in gzip (e.g., remainedreads1.fq.gz)\n"
                                + "         -o2 FILE        output file in gzip (e.g., remainedreads2.fq.gz)\n"
                                + "         -l              maximum long reads length [500]\n"
                                + "         -minov          minimum overlapping region [50]\n"
				+ "Input:   file1.fq.gz     input gzipped file with reads in fastq format\n"
				+ "                         (if not provided, input is taken from standard input)\n"
				+ "         file2.fq.gz     input gzipped file with reads in fastq format\n";


		boolean missinparrd12out = false;
		int corcnt = 0;
		String root_pth = "", inpth = "", flpth = "";

		String maxlrleg = "";
		int maxlrlegi = 0;
		double db = 0.0;

		String ftnm = "";
		
		// cumulative
		String correct_report = "", missed_report = "", general_report=""; 
		
		// cumulative correct report, missing report, general report for 
		// informative variables
		String cucrinf = "", cumsinf = ""; // current added to cumulative
		File outfl, inpfl;
		String lread = "", rread = "", longread = "", shortread1="", shortread2="";
		// test for pipe
		int pibal = 0;
		String minovs = "";
		int minovi=0;

		String[][] resargs = new String[7][2];
		resargs[0][0] = "pipeornot";
		resargs[1][0] = "input";
		resargs[2][0] = "output";
		resargs[3][0] = "outputread1path";
		resargs[4][0] = "outputread2path";
		resargs[5][0] = "longreadmaxlengthdef";
		resargs[5][1] = "500";
		resargs[6][0] = "minoveralpdef";
		resargs[6][1] = "50";
		// default value for long max length = 500
		Vector<String> infiles = new Vector<String>();
		String cur = "", jcur = "", inflss = "";

		String drct = "", testvalue = "";
		System.out.println(argmns.length);
		// int argindex = 0;
		int actualargsize = 0;
		for (int i=0; i<argmns.length; i++)
		{
		    if (argmns[i]==null)
		    {
		        actualargsize=i; 
// ( the size is always 1 + last index, and last index = size -1)
		        break;    
		    }
			actualargsize=i; 
		}
		if(actualargsize!=0)
		{
//		    if ((argmns[argmns.length - 2].equals("-o"))
//		            || (argmns[argmns.length - 2].equals("-o1"))
//		            || (argmns[argmns.length - 2].equals("-o2"))
//		            || (argmns[argmns.length - 2].equals("-l"))
//		            || (argmns[argmns.length - 2].equals("-minov")))
//		        drct = "pipe";
//		    else
//		        drct = "not_pipe";
		  System.out.println(argmns[actualargsize - 2] + "  size = " + actualargsize);
          if((argmns[actualargsize - 2].equals("-o"))
          || (argmns[actualargsize - 2].equals("-o1"))
          || (argmns[actualargsize - 2].equals("-o2"))
          || (argmns[actualargsize - 2].equals("-l"))
          || (argmns[actualargsize - 2].equals("-minov")))
                  drct = "pipe";
          else
                  drct = "not_pipe";
		    
		}
		general_report = general_report + drct + "\n";
	
		System.out.println(drct);
// 		Test for input file(s) existence.
		if (drct.equals("not_pipe")) {
			for (int j = argmns.length - 1; j > 0; j--) 
			{
				if(argmns[j]!=null)
				{
				if (!(argmns[j - 1].charAt(0) == '-')) {
					jcur = argmns[j];
					inpfl = new File(jcur);
					inpth = inpfl.getName();
					System.out.println(inpth);
					if (jcur.equals(inpth)) {
						root_pth = new java.io.File(".").getCanonicalPath();
						root_pth = root_pth + "/";
						jcur = root_pth + jcur;
					}
					if (testFileExistence(jcur)) {
						infiles.add(jcur);
						System.out.println("input file exists : " + jcur);
					} else {
						System.out.println("input file not found : " + jcur
								+ "\nExit.");
						Exit();
					}
				} else
					j = 0;
			}
		  }
		}

// 		Test(END) for input file(s) existence.
		if (infiles.size() == 1) {
			System.out.println("One fastq file as input!");
			resargs[1][0] = "InOnefastq";
		}
		if (infiles.size() == 2) 
		{
			if(infiles.get(0).equals(infiles.get(1)))
			{
				System.out.println("One input file entered twise.\nExit.");
				Exit();
			}
			System.out.println("Two fastq files as input!");
			resargs[1][0] = "InTwofastq";
		}
		if (infiles.size() > 2) 
		{
			if( (infiles.get(0).equals(infiles.get(1))) && (infiles.get(0).equals(infiles.get(2))))
			{
				System.out.println("More than two input valid files have been found.");
				System.out.println("One input file entered three times.\nExit.");
				Exit();
			}
			
			if( (infiles.get(0).equals(infiles.get(1))) || (infiles.get(0).equals(infiles.get(2))) || (infiles.get(1).equals(infiles.get(2))))
			{
				System.out.println("More than two input valid files have been found.");
				System.out.println("One input file entered twise.\nExit.");
				Exit();
			}

			System.out
					.println("More than two valid files as an input found.\nExit.");
			Exit();
		}
		String testq = "";
		int infcnt = 0;
		for (int u = 0; u < infiles.size(); u++) {
			infcnt = 0;
			testq = infiles.get(u);
			for (int uj = 0; uj < infiles.size(); uj++) {
				if (testq.equals(infiles.get(uj)))
					infcnt = infcnt + 1;
				if (infcnt > 1) {
					System.out
							.println("input file was typed two times in " +
							"arguments.\nExit.");
					Exit();
				}
			}
			// System.out.println(infiles.get(u));
		}

		for (int i = 0; i < argmns.length; i++) {
			if(argmns[i]!=null)
			{
// 		Test for outputs file root directory existence.
			if (argmns[i].equals("-o")) {
				if (i < argmns.length - 1) {
					longread = argmns[i + 1];
					outfl = new File(longread);
					ftnm = outfl.getName();
					if (longread.equals(ftnm)) {
						root_pth = new java.io.File(".").getCanonicalPath();
						root_pth = root_pth + "/";
						flpth = root_pth + ftnm;
						resargs[2][0] = "outputpath";
						resargs[2][1] = flpth;
						corcnt = corcnt + 1;
						cucrinf = "Output is OK!";
						correct_report = correct_report + cucrinf + "\n";
						System.out.println("Output is OK!");
						//if (pibal == 1)
						//	resargs[3][1] = flpth;

					} else if (!(longread.equals(ftnm))) {
						if (testParentDirectoryExistence(argmns[i + 1])) {
							resargs[2][0] = "outputpath";
							longread = argmns[i + 1];
							resargs[2][1] = longread;
							corcnt = corcnt + 1;
							cucrinf = "Output is OK!";
							correct_report = correct_report + cucrinf + "\n";
							System.out.println("Output is OK!");
						//	if (pibal == 1)
						//		resargs[3][1] = longread;

						} else {
							cumsinf = "Output root folder not found.";
				
							missed_report = missed_report + cumsinf + "\n";
							// System.out.println("Output root folder not found.");
						}
					}
				}
			}

			if (argmns[i].equals("-o1")) {
				if (i < argmns.length - 1) {
					shortread1 = argmns[i + 1];
					outfl = new File(shortread1);
					ftnm = outfl.getName();
					if (shortread1.equals(ftnm)) {
						root_pth = new java.io.File(".").getCanonicalPath();
						root_pth = root_pth + "/";
						flpth = root_pth + ftnm;
						resargs[3][0] = "outputRead1Path";
						resargs[3][1] = flpth;
						corcnt = corcnt + 1;
						cucrinf = "OutputRead1Path is OK!";
						correct_report = correct_report + cucrinf + "\n";
						System.out.println(cucrinf);
					//	if (pibal == 1)
					//		resargs[3][1] = flpth;

					} else if (!(shortread1.equals(ftnm))) {
						if (testParentDirectoryExistence(argmns[i + 1])) {
							resargs[3][0] = "outputRead1Path";
							shortread1 = argmns[i + 1];
							resargs[3][1] = shortread1;
							corcnt = corcnt + 1;
							cucrinf = "OutputRead1Path is OK!";
							correct_report = correct_report + cucrinf + "\n";
							System.out.println(cucrinf);
					//		if (pibal == 1)
					//			resargs[3][1] = shoread1;

						} else {
							cumsinf = "Output read1-root folder not found.";
							missed_report = missed_report + cumsinf + "\n";
							missinparrd12out=true;
						// System.out.println("Output root folder not found.");
						}
					}
				}
			}

			if (argmns[i].equals("-o2")) {
				if (i < argmns.length - 1) {
					shortread2 = argmns[i + 1];
					outfl = new File(shortread2);
					ftnm = outfl.getName();
					if (shortread2.equals(ftnm)) {
						root_pth = new java.io.File(".").getCanonicalPath();
						root_pth = root_pth + "/";
						flpth = root_pth + ftnm;
						resargs[4][0] = "outputRead2Path";
						resargs[4][1] = flpth;
						corcnt = corcnt + 1;
						cucrinf = "OutputRead2Path is OK!";
						correct_report = correct_report + cucrinf + "\n";
						System.out.println(cucrinf);
					//	if (pibal == 1)
					//		resargs[3][1] = flpth;

					} else if (!(shortread2.equals(ftnm))) {
						if (testParentDirectoryExistence(argmns[i + 1])) {
							resargs[4][0] = "outputRead2Path";
							shortread2 = argmns[i + 1];
							resargs[4][1] = shortread2;
							corcnt = corcnt + 1;
							cucrinf = "OutputRead2Path is OK!";
							correct_report = correct_report + cucrinf + "\n";
							System.out.println(cucrinf);
					//		if (pibal == 1)
					//			resargs[3][1] = shoread2;

						} else {
							cumsinf = "Output read2-root folder not found.";
							missed_report = missed_report + cumsinf + "\n";
							missinparrd12out=true;
						// System.out.println("Output root folder not found.");
						}
					}
				}
			}

			if (argmns[i].equals("-l")) {
				if (i < argmns.length - 1) {
					maxlrleg = argmns[i + 1];
					if (isNumeric(maxlrleg)) {
						db = Double.parseDouble(maxlrleg);
						db = Math.ceil(db);
						maxlrlegi = (int) db;
						if (maxlrlegi > 3) {
							cumsinf = "Maximum long read length OK!";
							resargs[5][1] = Integer.toString(maxlrlegi);
							corcnt = corcnt + 1;
						} else {
							System.out
									.println("Maximum long read length too " +
									"small.");
							Exit();
						}
					} else {
						System.out
								.println("Maximum long read length Not " + 
								"Numeric.");
						Exit();
					}

				} else {
					System.out.println("Maximum long read length is missing.");
					Exit();
				}
			}

			if (argmns[i].equals("-minov")) {
				if (i < argmns.length - 1) {
					minovs = argmns[i + 1];
					if (isNumeric(minovs)) {
						db = Double.parseDouble(minovs);
						db = Math.ceil(db);
						minovi = (int) db;
						if (minovi > 2) {
							cumsinf = "Minimum overlapping region is OK!";
							resargs[6][0] = "minoverlap";
							resargs[6][1] = Integer.toString(minovi);
							corcnt = corcnt + 1;
						} else {
							System.out
									.println("Minimum overlapping region is too " +
									"small.\nExit.");
							Exit();
						}
					} else {
						System.out
								.println("Minimum overlapping region is Not " + 
								"Numeric.\nExit.");
						Exit();
					}

				} else {
					System.out.println("Minimum overlapping region error.\nExit.");
					Exit();
				}
			}
		}
		}
		System.out.println(missed_report);

//		if( !( resargs[3][1].isEmpty()) && !( resargs[4][1].isEmpty()) )
		if( ( (resargs[3][1]==null) && (resargs[4][1]!=null) ) || ( (resargs[3][1]!=null) && (resargs[4][1]==null) ) )
		{		
			System.out.println("Both output for failed to overlap reads are required\n when optional '-o1' and '-o2' are used.\nOne output for failed to overlap reads is missing\nExit.");
			Exit();
		}


		if( (resargs[3][1]!=null) && (resargs[4][1]!=null) &&  (resargs[2][1]!=null) )		
		{
			if( resargs[3][1].equals(resargs[4][1]))
			{
				System.out.println("outputs for failed to overlap reads are same.\nExit.");
				Exit();
			}
			if( ( resargs[2][1].equals(resargs[3][1])) || ( resargs[2][1].equals(resargs[4][1])) )
			{
				System.out.println("outputs for Long reads and failed to overlap reads are same.\nExit.");
				Exit();
			}
		}

		//for (int i=2; i<5; i++)
		//	System.out.println(resargs[i][1]);
		// System.out.println("before exit");
		if(missinparrd12out==true)
			Exit();
		// for(int i=0; i<infiles.size(); i++)
		// {
		// if(i<infiles.size()-1)
		// inflss=inflss + infiles.get(i) + ",";
		// if(i==infiles.size()-1)
		// inflss=inflss + infiles.get(i);
		// }

		if (infiles.size() > 1) {
			for (int i = infiles.size() - 1; i > -1; i--) {
				if (i > 0)
					inflss = inflss + infiles.get(i) + ",";
				if (i == 0)
					inflss = inflss + infiles.get(i);
			}
		}

		// Error fixed on 05/14/2015
		// if only one fastq file was as an input the path was missed
		// fixed with the following if statement
		if (infiles.size() == 1)
			inflss = inflss + infiles.get(0);

		File tszpd;
		if (infiles.size() > 0) {
			tszpd = new File(infiles.get(0));
			if (isFileGZipped(tszpd)) {
				// resargs[1][0] = "leftrightreaddfqgz";
				resargs[0][1] = "notpipe_gz";
			}
			if (!(isFileGZipped(tszpd))) {
				// resargs[1][0] = "leftrightreaddfqgz";
				resargs[0][1] = "notpipe";
			}
		}
		// Test if the input file is a bam file that contains left and right
		// read (needs implementation)
		if (infiles.size() < 1) {
			resargs[1][0] = "bamfile";
			resargs[0][1] = "pipe";
			System.out.println("Piped bam file expected.");
			// resargs[1][1] = "none";
		}
		if(resargs[2][1]==null) 		
		{
			System.out.println("output for long reads is missing.\nExit.");
			Exit();
		}

		if (infiles.size() > 0) // string >0
		{
			resargs[1][1] = inflss;
		}
		System.out.println("End of testArgsforFrags method!");

		return resargs;
	}
	//test frag (END)





	public static String[][] testArgsforFragsBack08(String[] argmns)
			throws IOException {

		System.out.println("Start of testArgsforFrags method!");
		int corcnt = 0;
		String root_pth = "", inpth = "", flpth = "";

		String maxlrleg = "";
		int maxlrlegi = 0;
		double db = 0.0;

		String ftnm = "";
		String correct_report = "", missed_report = "";
		String cucrinf = "", cumsinf = "";
		File outfl, inpfl;
		String lread = "", rread = "", longread = "";
		// test for pipe
		int pibal = 0;

		String[][] resargs = new String[4][2];
		resargs[0][0] = "pipeornot";
		resargs[1][0] = "input";
		resargs[2][0] = "output";
		resargs[3][0] = "longreadmaxlength";
		resargs[3][1] = "500";
		// default value for long max length = 500

		Vector<String> infiles = new Vector<String>();
		String cur = "", jcur = "", inflss = "";

		String drct = "", testvalue = "";
		if ((argmns[argmns.length - 2].equals("-o"))
				|| (argmns[argmns.length - 2].equals("-l")))
			drct = "pipe";
		else
			drct = "not_pipe";
		System.out.println(drct);

		if (drct.equals("not_pipe")) {
			for (int j = argmns.length - 1; j > 0; j--) {
				if (!(argmns[j - 1].charAt(0) == '-')) {
					jcur = argmns[j];
					inpfl = new File(jcur);
					inpth = inpfl.getName();
					System.out.println(inpth);
					if (jcur.equals(inpth)) {
						root_pth = new java.io.File(".").getCanonicalPath();
						root_pth = root_pth + "/";
						jcur = root_pth + jcur;
					}
					if (testFileExistence(jcur)) {
						infiles.add(jcur);
						System.out.println("input file exists : " + jcur);
					} else {
						System.out.println("input file not found : " + jcur
								+ "\nExit.");
						Exit();
					}
				} else
					j = 0;
			}
		}

		if (infiles.size() == 1) {
			System.out.println("One fastq file as input!");
			resargs[1][0] = "InOnefastq";
		}
		if (infiles.size() == 2) {
			System.out.println("Two fastq files as input!");
			resargs[1][0] = "InTwofastq";
		}
		if (infiles.size() > 2) {
			System.out
					.println("More than two valid files as an input found. Exit.");
			Exit();
		}
		String testq = "";
		int infcnt = 0;
		for (int u = 0; u < infiles.size(); u++) {
			infcnt = 0;
			testq = infiles.get(u);
			for (int uj = 0; uj < infiles.size(); uj++) {
				if (testq.equals(infiles.get(uj)))
					infcnt = infcnt + 1;
				if (infcnt > 1) {
					System.out
							.println("input file was typed two times in " +
							"arguments.\nExit.");
					Exit();
				}
			}
			// System.out.println(infiles.get(u));
		}

		for (int i = 0; i < argmns.length; i++) {
			if (argmns[i].equals("-o")) {
				if (i < argmns.length - 1) {
					longread = argmns[i + 1];
					outfl = new File(longread);
					ftnm = outfl.getName();
					if (longread.equals(ftnm)) {
						root_pth = new java.io.File(".").getCanonicalPath();
						root_pth = root_pth + "/";
						flpth = root_pth + ftnm;
						resargs[2][0] = "outputpath";
						resargs[2][1] = flpth;
						corcnt = corcnt + 1;
						cucrinf = "Output is OK!";
						correct_report = correct_report + cucrinf + "\n";
						System.out.println("Output is OK!");
						if (pibal == 1)
							resargs[3][1] = flpth;

					} else if (!(longread.equals(ftnm))) {
						if (testParentDirectoryExistence(argmns[i + 1])) {
							resargs[2][0] = "outputpath";
							longread = argmns[i + 1];
							resargs[2][1] = longread;
							corcnt = corcnt + 1;
							cucrinf = "Output is OK!";
							correct_report = correct_report + cucrinf + "\n";
							System.out.println("Output is OK!");
							if (pibal == 1)
								resargs[3][1] = longread;

						} else {
							cumsinf = "Output root folder not found.";
							missed_report = missed_report + cumsinf + "\n";
						// System.out.println("Output root folder not found.");
						}
					}
				}
			}
			if (argmns[i].equals("-l")) {
				if (i < argmns.length - 1) {
					maxlrleg = argmns[i + 1];
					if (isNumeric(maxlrleg)) {
						db = Double.parseDouble(maxlrleg);
						db = Math.ceil(db);
						maxlrlegi = (int) db;
						if (maxlrlegi > 3) {
							cumsinf = "Maximum long read length OK!";
							resargs[3][1] = Integer.toString(maxlrlegi);
							corcnt = corcnt + 1;
						} else {
							System.out
									.println("Maximum long read length too " +
									"small.");
							Exit();
						}
					} else {
						System.out
								.println("Maximum long read length Not " + 
								"Numeric.");
						Exit();
					}

				} else {
					System.out.println("Maximum long read length is missing.");
					Exit();
				}
			}

		}

		// for(int i=0; i<infiles.size(); i++)
		// {
		// if(i<infiles.size()-1)
		// inflss=inflss + infiles.get(i) + ",";
		// if(i==infiles.size()-1)
		// inflss=inflss + infiles.get(i);
		// }

		if (infiles.size() > 1) {
			for (int i = infiles.size() - 1; i > -1; i--) {
				if (i > 0)
					inflss = inflss + infiles.get(i) + ",";
				if (i == 0)
					inflss = inflss + infiles.get(i);
			}
		}

		// Error fixed on 05/14/2015
		// if only one fastq file was as an input the path was missed
		// fixed with the following if statement
		if (infiles.size() == 1)
			inflss = inflss + infiles.get(0);

		File tszpd;
		if (infiles.size() > 0) {
			tszpd = new File(infiles.get(0));
			if (isFileGZipped(tszpd)) {
				// resargs[1][0] = "leftrightreaddfqgz";
				resargs[0][1] = "notpipe_gz";
			}
			if (!(isFileGZipped(tszpd))) {
				// resargs[1][0] = "leftrightreaddfqgz";
				resargs[0][1] = "notpipe";
			}
		}
		// Test if the input file is a bam file that contains left and right
		// read (needs implementation)
		if (infiles.size() < 1) {
			resargs[1][0] = "bamfile";
			resargs[0][1] = "pipe";
			System.out.println("Piped bam file expected.");
			// resargs[1][1] = "none";
		}

		if (infiles.size() > 0) // string >0
		{
			resargs[1][1] = inflss;
		}
		System.out.println("End of testArgsforFrags method!");

		return resargs;
	}
	// test frag (END)
	// test noisetab (START)

	
	
	
	//extract -a T -p chr2:12307676 my.bam
	//output in fastq
	public static String[][] testArgsforExtract(String[] argmns)
			throws IOException
	{
		System.out
		.println("\n\nStart of test arguments for Extract Reads!");
		
		String info = "\nUsage: Leucippus extract [options] bam.file\n\n"
				+ "Options: -o           FILE   results\n"
				+ "         -a                  alelle_X [X could be A,C,T,G]"
				+ "         -p                  chr:position\n\n";

		String[][] resargs = new String[4][2];
		boolean donotproceed=false;
		String missed_report="", cumsinf="";
		String reason="";
		String bmcorrect="", correct_report="", cucrinf="";
		String alle="", position="";
		
		resargs[0][0] = "bamsq";
		resargs[1][0] = "alternative";
		resargs[2][0] = "position";
		resargs[3][0] = "output";

		if ((argmns[argmns.length - 2].equals("-o"))
				|| (argmns[argmns.length - 2].equals("-a"))
				|| (argmns[argmns.length - 2].equals("-p"))) 
		{
			bmcorrect = "bammissing";
			resargs[0][0] = bmcorrect;
			cumsinf = "Bam file(s) is(are) missing.";
			missed_report = missed_report + cumsinf + "\n";
			donotproceed=true;
			reason=reason+cumsinf+"\n";
		}
		if ((argmns[argmns.length - 1].equals("-o"))
				|| (argmns[argmns.length - 1].equals("-a"))
				|| (argmns[argmns.length - 1].equals("-p"))) {
			bmcorrect = "bammissing";
			resargs[0][0] = bmcorrect;
		} else {
			bmcorrect = "bamexpect";
			resargs[0][0] = bmcorrect;
		}
		String outtest = "", ftnm = "", root_pth = "", flpth = "";
		// resargs[0][0] = "bams";
		int bmcnt = 0, corcnt=0;

		File bmfl;
		String crarg = "", bamfls="";
		String outpth="";
		// System.out.println(" --- > " + resargs[0][0]);
		// System.out.println("bc > " + bmcorrect);
		if (resargs[0][0].equals("bamexpect")) {
			for (int i = argmns.length - 1; i >= 2; i--) {
				if (!(argmns[i - 1].charAt(0) == '-')) {
					crarg = argmns[i];
					bmfl = new File(crarg);
					ftnm = bmfl.getName();
					if (crarg.equals(ftnm)) {
						root_pth = new java.io.File(".").getCanonicalPath();
						root_pth = root_pth + "/";
						crarg = root_pth + ftnm;
					}
					if (testFileExistence(crarg)) {
						bmcnt = bmcnt + 1;
						if (bmcnt == 1)
							bamfls = crarg;
						if (bmcnt > 1)
							bamfls = crarg + "," + bamfls;
					} else if (!(testFileExistence(argmns[i])))
						System.out.println("bam file '" + argmns[i]
								+ "' not found.");
				} else if (argmns[i - 1].charAt(0) == '-')
					i = 1;
			}
		}
		if (!(bamfls.equals("")))
		{
			resargs[0][0] = "bams";
			resargs[0][1] = bamfls;
			corcnt = corcnt + 1;
		}
		else if (bamfls.equals(""))
		{
			cumsinf="There aren't bamfiles to process.";
			donotproceed=true;
			reason=reason+cumsinf+"\n";
		}
		
		for (int i = 0; i < argmns.length; i++) 
		{
			if (argmns[i].equals("-o"))
			{
				if (i < argmns.length - 1)
				{
					outtest = argmns[i + 1];
					File outfl = new File(outtest);
					ftnm = outfl.getName();
					if (outtest.equals(ftnm))
					{
						root_pth = new java.io.File(".").getCanonicalPath();
						root_pth = root_pth + "/";
						flpth = root_pth + ftnm;
						resargs[3][0] = "outputpath";
						resargs[3][1] = flpth;
						corcnt = corcnt + 1;
						cucrinf = "Output is OK!";
					}
					else if (!(outtest.equals(ftnm))) 
					{
					if (testParentDirectoryExistence(argmns[i + 1]))
					{
						resargs[3][0] = "outputpath";
						outpth = argmns[i + 1];
						resargs[3][1] = outpth;
						corcnt = corcnt + 1;
						cucrinf = "Output is OK!";
						correct_report = correct_report + cucrinf + "\n";
							// System.out.println("Output is OK!");
					}
					else 
					{
						cumsinf = "Output root folder not found.";
						missed_report = missed_report + cumsinf + "\n";
						donotproceed=true;
						reason=reason+cumsinf+"\n";
							// System.out.println("Output root folder not found.");
					}
				}
			  }
			}
			if ((argmns[i].equals("-a")) && (i < argmns.length - 1))
			{
				alle=argmns[i + 1];
				if( (alle.equals("A")) || 
					(alle.equals("C")) || 
					(alle.equals("T")) || 
					(alle.equals("G")))
				{
					resargs[1][1] = alle;
					corcnt = corcnt + 1;
					cucrinf = "Alelle is OK!";
				}
				else
				{
					cumsinf = "Alelle not stated correctly.";
					missed_report = missed_report + cumsinf + "\n";
					donotproceed=true;
					reason=reason+cumsinf+"\n";
				}
			}
			
			if ((argmns[i].equals("-p")) && (i < argmns.length - 1))
			{
				position=argmns[i + 1];
				//chr2:12307676
				String[] pswds = position.split(":");
				if(pswds.length==2)
				{
					if( (isPosUnsignInteger(pswds[1])) && isChromosome(pswds[0]))
					{
						resargs[2][1]=position + "-" + pswds[1];
						corcnt = corcnt + 1;
						cucrinf = "Position is OK!";					
					}
					else
					{
						cumsinf = "Position has not been stated correct.";
						missed_report = missed_report + cumsinf + "\n";
						donotproceed=true;
						reason=reason+cumsinf+"\n";
					}
				}
				else
				{
					cumsinf = "Position must have the following format :\nchr2:12307676";
					missed_report = missed_report + cumsinf + "\n";
					donotproceed=true;
					reason=reason+cumsinf+"\n";
				}
			}
		}

		System.out
		.println("\n....End of test arguments for Extract Reads!\n\n");
		
		if(donotproceed==true)
		{
			System.out.println("Error/Missing :\n" + reason+"Exit.\n");
			Exit();
		}
		
		if(resargs[3][0].equals("output"))
		{
			resargs[3][1]="p";
			corcnt = corcnt + 1;
		}
		return resargs;
	}

	public static boolean isChromosome(String chr)
	{
		boolean result = true;
		int chnm = 0;
		String chlt="";
		// potential number
		if(chr.length()>3)
		{
			if(chr.substring(0,3).equals("chr"))
			{
				if((chr.charAt(3)=='X') || (chr.charAt(3)=='Y'))
					result = true;
				else if(isPosUnsignInteger(chr.substring(3,chr.length())))
				{
					if ( (Integer.parseInt(chr.substring(3,chr.length()))<23) &&
					    (Integer.parseInt(chr.substring(3,chr.length()))>0) )
						result = true;
					else if ( (Integer.parseInt(chr.substring(3,chr.length()))>=23) &&
					    (Integer.parseInt(chr.substring(3,chr.length()))<=0) )
						result = false;
				}
				else
					result = false;
			}
		}
		else if(chr.length()<=2)
		{
			if(isPosUnsignInteger(chr))
			{
				chnm=Integer.parseInt(chr);
				if((chnm>0) && (chnm<23))
					result = true;
				else
					result=false;	
			}
			else if( (chr.equals("X")) || (chr.equals("Y")) )
				result = true;
			else		
				result=false;
		}
		return result; 
	}

	/**
	 * @param argmns : String[]
	 * @return resargs : Strjng[][]
	 * @throws IOException
	 **/
	public static String[][] testArgsforTables(String[] argmns)
			throws IOException {

		System.out
				.println("\n\nStart of test arguments for creating noise " +
				"tables .....\n");

		String info = "\nUsage: leucippus noisetab [options] file1.bam [file2.bam, ...]\n\n"
				+ "Options:     -interval FILE  file with intervals\n"
				+ "             -ref      FILE  FASTA file with reference sequence (could be .gz file)\n"
				+ "             -o        FILE  output file\n"
				+ "             -q              quality cut off 0-100 [20]\n"
				+ "             -d              distance cut off [-1]\n\n";

		// Default values for base quality and distance qual = 20 and distance =
		// -1
		boolean donotproceed=false;
		String reason="";
		int qbal = 0, dbal = 0, pbal=0; // indicators for quality and distance
		// if ==0 then the default value will be used.
		int dis = 0;
		// Variables that were created in order to allow output file name to be
		// in the application current path
		String outtest = "", ftnm = "", root_pth = "", flpth = "";
		// (fixed size) Array table that holds two elements in each row : the
		// description[x][0] and the value[x][1]
		// The size must be changed here if a new variable needs to be
		// used.
		int padi=1;
		String pad="";
		String[][] resargs = new String[7][2];
		String bmcorrect = "";
		String missed_report = "", cumsinf = "";
		String correct_report = "", cucrinf = "";
		double distanced = 0.0;
		resargs[0][0] = "bamsq";
		resargs[1][0] = "interval";
		resargs[2][0] = "output";
		resargs[3][0] = "quality";
		resargs[4][0] = "distance";
		resargs[5][0] = "reference";
		resargs[6][0] = "defpad";
		resargs[6][1] = "1";

		String bamfls = "", bampth = "", prmpath = "", outpth = "", 
		dist = "", qual = "", refpath = "";
		for (int i1 = 0; i1 < 5; i1++)
			for (int i2 = 0; i2 < 2; i2++)
				resargs[i1][i2] = "";
		if ((argmns[argmns.length - 2].equals("-interval"))
				|| (argmns[argmns.length - 2].equals("-q"))
				|| (argmns[argmns.length - 2].equals("-o"))
				|| (argmns[argmns.length - 2].equals("-d"))
				|| (argmns[argmns.length - 2].equals("-pad"))
				|| (argmns[argmns.length - 2].equals("-ref"))) {
			bmcorrect = "bammissing";
			resargs[0][0] = bmcorrect;
		}
		if ((argmns[argmns.length - 1].equals("-interval"))
				|| (argmns[argmns.length - 1].equals("-q"))
				|| (argmns[argmns.length - 1].equals("-o"))
				|| (argmns[argmns.length - 1].equals("-d"))
				|| (argmns[argmns.length - 2].equals("-pad"))
				|| (argmns[argmns.length - 1].equals("-ref"))) {
			bmcorrect = "bammissing";
			resargs[0][0] = bmcorrect;
		} else {
			bmcorrect = "bamexpect";
			resargs[0][0] = bmcorrect;
		}
		// creat counters to test the uniquenes of entered variable-values
		int bcn = 0, icn = 0, qcn = 0, dcn = 0, ocn = 0, pdcn=0, refcnt = 0;
		// ------If there is double instance then Exit();

		for (int i = 0; i < argmns.length; i++) {
			// if( (argmns[i].equals("-bam")) || (argmns[i].equals("-dir")))
			// if (argmns[i].equals("-bam"))
			// bcn = bcn + 1;
			if (argmns[i].equals("-interval"))
				icn = icn + 1;
			if (argmns[i].equals("-q"))
				qcn = qcn + 1;
			if (argmns[i].equals("-o"))
				ocn = ocn + 1;
			if (argmns[i].equals("-d"))
				dcn = dcn + 1;
			if (argmns[i].equals("-pad"))
				pdcn = pdcn + 1;
			if (argmns[i].equals("-ref"))
				refcnt = refcnt + 1;
		}
// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
// ___________________________________________________________________________________________________
		if ((bcn > 1) || (icn > 1) || (qcn > 1) || (ocn > 1) || (dcn > 1)
				|| (pdcn > 1) || (refcnt > 1)) {
			System.out
					.println("Found more than one instace of operators, \nExit.");
			Exit();
		}

		if ((icn != 1) || (refcnt != 1)) {
			System.out
					.println("leucippus needs one reference file and one interval file to proceed with 'noisetab' command.\nExit.");
			Exit();
		}
		if (dcn == 0)
			missed_report = "Distance value is missing. Distance default value(not applicable) will be used.\n";

		if (qcn == 0)
			missed_report = missed_report
					+ "Base Quality value is missing. Base Quality default value will be used. q = 20\n";
		if (pdcn == 0)
			missed_report = missed_report
					+ "Pading value is missing. Padding default value will be used. -pad = 1\n";
		if (ocn == 0)
			missed_report = missed_report + "Output file is missing.\n";

		if (resargs[0][0].equals("bammissing")) {
			System.out
					.println("bam files must provided as arguments (missing bam file(s)).\nExit.");
			Exit();
		}
		// ------- If there is double instance then Exit();

		int corcnt = 0; // counter that terminates the program if the entered
						// input
						// arguments are not enough
		// -q and -d are set initialized (by default) to 20 and -1
		// correspondingly;
		// and the 'corcnt' counter will be incremented for them anyway.
		int quali = 20, disti = -1;

		// resargs[0][0] = "bams";
		int bmcnt = 0;

		File bmfl;
		String crarg = "";
		// System.out.println(" --- > " + resargs[0][0]);
		// System.out.println("bc > " + bmcorrect);
		if (resargs[0][0].equals("bamexpect")) {
			for (int i = argmns.length - 1; i >= 2; i--) {
				if (!(argmns[i - 1].charAt(0) == '-')) {
					crarg = argmns[i];
					bmfl = new File(crarg);
					ftnm = bmfl.getName();
					if (crarg.equals(ftnm)) {
						root_pth = new java.io.File(".").getCanonicalPath();
						root_pth = root_pth + "/";
						crarg = root_pth + ftnm;
					}
					if (testFileExistence(crarg)) {
						bmcnt = bmcnt + 1;
						if (bmcnt == 1)
							bamfls = crarg;
						if (bmcnt > 1)
							bamfls = crarg + "," + bamfls;
					} else if (!(testFileExistence(argmns[i])))
						System.out.println("bam file '" + argmns[i]
								+ "' not found.");
				} else if (argmns[i - 1].charAt(0) == '-')
					i = 1;
			}
		}
		if (!(bamfls.equals("")))
		{
			resargs[0][0] = "bams";
			resargs[0][1] = bamfls;
			corcnt = corcnt + 1;
		}
		else if (bamfls.equals(""))
		{
			cumsinf="There aren't bamfiles to process.";
			donotproceed=true;
			reason=reason+cumsinf+"\n";
		}
		
		resargs[2][0] = "output";
		for (int i = 0; i < argmns.length; i++) {
			if ((argmns[i].equals("-interval")) && (i < argmns.length - 1)) {
				if (testFileExistence(argmns[i + 1])) {

					prmpath = argmns[i + 1];
					resargs[1][0] = "interval";
					resargs[1][1] = prmpath;
					corcnt = corcnt + 1;
					cucrinf = "Interval file is OK!";
					correct_report = correct_report + cucrinf + "\n";
					// System.out.println("Interval file is OK!");
				} else if (!(testFileExistence(argmns[i + 1]))) {
					cumsinf = "Interval File '" + argmns[i + 1]
							+ "' not found.";
					missed_report = missed_report + cumsinf + "\n";
					donotproceed=true;
					reason=reason+cumsinf+"\n";
					// System.out.println("Interval file '" + argmns[i+1] +
					// "' not found.");
				}
			}
			if ((argmns[i].equals("-ref")) && (i < argmns.length - 1)) {
				if (testFileExistence(argmns[i + 1])) {

					refpath = argmns[i + 1];
					resargs[5][0] = "refpath";
					resargs[5][1] = refpath;
					corcnt = corcnt + 1;
					cucrinf = "Genome Reference File is OK!";
					correct_report = correct_report + cucrinf + "\n";
					// System.out.println("Genome Reference File is OK!");
				} else if (!(testFileExistence(argmns[i + 1]))) {
					cumsinf = "Genome Reference File '" + argmns[i + 1]
							+ "' not found.";
					missed_report = missed_report + cumsinf + "\n";
					donotproceed=true;
					reason=reason+cumsinf+"\n";
					// System.out.println("Genome Reference File '" +
					// argmns[i+1] + "' not found.");
				}
			}

			if (argmns[i].equals("-o")) {
				if (i < argmns.length - 1) {
					outtest = argmns[i + 1];
					File outfl = new File(outtest);
					ftnm = outfl.getName();
					if (outtest.equals(ftnm)) {
						root_pth = new java.io.File(".").getCanonicalPath();
						root_pth = root_pth + "/";
						flpth = root_pth + ftnm;
						resargs[2][0] = "outputpath";
						resargs[2][1] = flpth;
						corcnt = corcnt + 1;
						cucrinf = "Output is OK!";
						correct_report = correct_report + cucrinf + "\n";
						// System.out.println("Output is OK!");
					} else if (!(outtest.equals(ftnm))) {
						if (testParentDirectoryExistence(argmns[i + 1])) {
							resargs[2][0] = "outputpath";
							outpth = argmns[i + 1];
							resargs[2][1] = outpth;
							corcnt = corcnt + 1;
							cucrinf = "Output is OK!";
							correct_report = correct_report + cucrinf + "\n";
							// System.out.println("Output is OK!");
						} else {
							cumsinf = "Output root folder not found.";
							missed_report = missed_report + cumsinf + "\n";
							donotproceed=true;
							reason=reason+cumsinf+"\n";
							// System.out.println("Output root folder not found.");
						}
					}
				}
			}
			if (argmns[i].equals("-q")) {
				if (i < argmns.length - 1)
					// restriction for integer
					qual = argmns[i + 1];
				if (isPosUnsignInteger(qual)) {
					resargs[3][0] = "quality";
					resargs[3][1] = qual;
					corcnt = corcnt + 1;
					qbal = 1;
					cucrinf = "Quality is OK!";
					correct_report = correct_report + cucrinf + "\n";
					// System.out.println("Quality is OK!");
				} else {
					resargs[3][0] = "quality";
					cumsinf = "Quality is not valid.";
					missed_report = missed_report + cumsinf + "\n";
					// System.out.println("Quality is not valid.");
					qbal = 3;
				}
			}

			if (argmns[i].equals("-pad")) {
				if (i < argmns.length - 1)
					// restriction for integer
					pad = argmns[i + 1];
				if (isPosUnsignInteger(pad)) {
					resargs[6][0] = "pad";
					resargs[6][1] = pad;
					corcnt = corcnt + 1;
					cucrinf = "Pading is OK!";
					correct_report = correct_report + cucrinf + "\n";
					pbal=1;
					// System.out.println("Quality is OK!");
				} else {
					resargs[6][0] = "paderror";
					cumsinf = "Pading is not valid.";
					missed_report = missed_report + cumsinf + "\n";
					// System.out.println("Quality is not valid.");
					pbal = 3;
				}
			}


			if (argmns[i].equals("-d")) {
				if (i < argmns.length - 1) {
					dist = argmns[i + 1];

					if (isNumeric(dist)) // if it is numeric is OK(in the future
					{
						distanced = Double.parseDouble(dist);
						dis = (int) Math.ceil(distanced);
						// if (isPosUnsignInteger(dist))
						if (distanced > 0) {

							resargs[4][0] = "distance";
							// dist=argmns[i+1];
							resargs[4][1] = Integer.toString(dis);
							corcnt = corcnt + 1;
							dbal = 1;
							cucrinf = "Distance is OK!";
							correct_report = correct_report + cucrinf + "\n";
							// System.out.println("Distance is OK!");

						} else if (distanced < 0) {
							resargs[4][0] = "distance";
							cucrinf = "Distance not applicable;  OK!";
							dist = "-1";
							resargs[4][1] = dist;
							dbal = 1;
							corcnt = corcnt + 1;
						}
					} else // it is not numeric
					// therefore it must be corrected
					{
						resargs[4][0] = "distance";
						cumsinf = "Distance is not valid.";
						missed_report = missed_report + cumsinf + "\n";
						// System.out.println("Distance is not valid.");
						donotproceed=true;
						reason=reason+cumsinf+"\n";
						dbal = 3;
					}
				}
			}
		}

		if (bamfls.length() >= 1) {
			corcnt = corcnt + 1;
			resargs[0][1] = bamfls;
			cucrinf = "Bam files is OK!";
			correct_report = correct_report + cucrinf + "\n"; // System.out.println("Bam files is OK!");
		} else if (bamfls.length() < 1) {
			// System.out.println("No Bam files to process.");
			cumsinf = "No Bam files to process.";
			missed_report = missed_report + cumsinf + "\n";
		}
		// System.out.println(dbal);
		if (qbal == 0) {
			resargs[3][0] = "defquality";
			resargs[3][1] = "20";
			corcnt = corcnt + 1;
			cucrinf = "Default base quality cut off = 20  OK!";
			correct_report = correct_report + cucrinf + "\n";
			// System.out.println("Default base quality cut off = 20  OK!");
		}
		if (dbal == 0) {
			resargs[4][0] = "defdistance";
			resargs[4][1] = "-1";
			corcnt = corcnt + 1;
			cucrinf = "Default Distance cut off = -1  OK!";
			correct_report = correct_report + cucrinf + "\n";
			// System.out.println("Default Distance cut off = -1  OK!");
		}
		if (pbal == 0)
		{
			resargs[6][0] = "defpad";
			resargs[6][1] = "1";
			corcnt = corcnt + 1;
			cucrinf = "Default pading value = 1  OK!";
			correct_report = correct_report + cucrinf + "\n";

		}
		System.out.println("Missing|Errors Values :\n" + missed_report + "\n");
		System.out
				.println("Correct|Default Values :\n" + correct_report + "\n");

		System.out.println("corcnt = " + corcnt);

		for (int h = 0; h < resargs.length; h++)
			System.out.println(resargs[h][0] + " : " + resargs[h][1]);

		System.out
				.println("\n....End of test arguments for creating noise table.\n\n");
		if(donotproceed==true)
		{
			System.out.println("Error/Missing :\n" + reason+"Exit.\n");
			Exit();
		}
		if (corcnt < 7) {
			Exit();
		}
		return resargs;
	}

    	// test graph (START)
	/**
	 * 'testArgsforGraphs' method --> to be completed
	 * 
	 * @param argmns
	 * @return
	 * @throws IOException
	 * @throws InterruptedException
	 **/
	public static String[][] testArgsforGraphs(String[] argmns)
			throws IOException, InterruptedException {
		String graphinfo = "\nUsage: leucippus graph [options] -o <prefix> table1.file [table2.file]\n\n"
				+ "Options:    -type            graph type: pvalue|mutrate [pvalue]\n"
				+ "	    -coverage        	minimum site/position coverage [100]\n"
				+ "            -range  DOUBLE  minimum cutoff to exclude germline from noise graph [0.95]\n"
				+ "                                found/expected when base found is same as reference(expected) base\n"
				+ "            -overlap INT     use positions in overlapping 3'-ends of reads.\n"
				+ "                             The number specifies read length.\n"
				+ "                             Only useful for analysis of amplicon-seq data.\n"
				+ "            -fpval           generate frequency p-value table\n"
				+ "            -o               prefix for output files: prefix.pdf,\n"
				+ "                                                      prefix.fpvtb[1,2].tsv\n\n";
		System.out.println("\ntestArgsforGraphs method Start ------------");

		String tblpthstr = "";
		String[][] resargs = new String[6][2];
		resargs[0][0] = "type";
		resargs[0][1] = "";
		resargs[1][0] = "checkTable(s)"; // String that could contain one or two
											// // table(absolutepaths)
		resargs[2][0] = "";
		resargs[3][0] = "dfault";
		resargs[3][1] = "0"; // Default overlap value = 0
		resargs[4][0] = "";
		resargs[4][1] = "0";
		resargs[5][0] = "";
		resargs[5][1] = "0.95";
		double rangedb = 0.0, fnrgdb = 0.0;

		String grtst = "", crarg = "", ranges = "";
		String type = "", tbsrce = "", grphout = "", overlaps = "", coverages = "", germlineAFs = "";
		// gprim="";
		int overlapi = 0, corcnt = 0;
		;
		String root_pth = "", inpth = "", flpth = "";
		String ftnm = "", jcur = "";
		File outfl, inpfl;
		String correct_report = "", missed_report = "", tbinf = "";
		String cucrinf = "", cumsinf = "", ftbstr = "";
		int typecnt = 0, tbsrcecnt = 0, grphoutcnt = 0, fpvalcnt = 0, ovrlapcnt = 0, cvcnt = 0, grmlcnt = 0, tbfcnt = 0;
		;
		// , gprmcnt=0;
		// constuct a common (type + fpval) string and keep only one version
		// ("groutfpval")

		for (int i = 0; i < argmns.length; i++) {
			if (argmns[i] != null) {
				if (argmns[i].equals("-type"))
					typecnt = typecnt + 1;

				if (argmns[i].equals("-coverage"))
					cvcnt = cvcnt + 1;

				if (argmns[i].equals("-o")) {
					grphoutcnt = grphoutcnt + 1;
					resargs[2][0] = resargs[2][0] + "grout";
				}
				if (argmns[i].equals("-overlap"))
					ovrlapcnt = ovrlapcnt + 1;
				if (argmns[i].equals("-range"))
					grmlcnt = grmlcnt + 1;
				if (argmns[i].equals("-fpval")) {
					fpvalcnt = fpvalcnt + 1;
					resargs[2][0] = resargs[2][0] + "fpval";
				}
			}
		}

		if ((argmns[argmns.length - 1].equals("-overlap"))
				|| (argmns[argmns.length - 1].equals("-range"))
				|| (argmns[argmns.length - 1].equals("-coverage"))
				|| (argmns[argmns.length - 1].equals("-type"))
				|| (argmns[argmns.length - 1].equals("-o"))
				|| (argmns[argmns.length - 2].equals("-overlap"))
				|| (argmns[argmns.length - 2].equals("-range"))
				|| (argmns[argmns.length - 2].equals("-coverage"))
				|| (argmns[argmns.length - 2].equals("-type"))
				|| (argmns[argmns.length - 2].equals("-o"))) {
			tbinf = "tablemiss";
			resargs[1][0] = tbinf;
			tbsrcecnt = 0;
			System.out
					.println("Missing Table(s) as Arguments at the end. \n Exit");
			Exit();
		} else {
			tbsrcecnt = tbsrcecnt + 1;
			resargs[1][0] = "TableExpect";
		}
		//

		if ((resargs[2][0].equals("groutfpval"))
				|| (resargs[2][0].equals("fpvalgrout")))
			resargs[2][0] = "groutfpval";
		System.out.println("graph  : " + resargs[2][0]);

		// System.out.println("testArgsforGraphs method  2 ");

		if ((typecnt > 1) || (tbsrcecnt > 1) || (grphoutcnt > 1)
				|| (ovrlapcnt > 1) || (cvcnt > 1) || (grmlcnt > 1)) {
			System.out.println("Arguments must be stated once. \n Exit");
			Exit();
		}

		if ((tbsrcecnt == 0) || ((grphoutcnt == 0) && (fpvalcnt == 0))) {
			System.out
					.println("Table(s) input, (graph output or table output (path-prefix) are required.\n Exit");
			System.out.println("typecnt = " + typecnt + " tablecnt = "
					+ tbsrcecnt + " outcnt = " + grphoutcnt
					+ " optional coverage " + cvcnt);
			Exit();
		}

		if (cvcnt == 0) {
			resargs[4][1] = "100";
			System.out
					.println("Coverage default value will be used : Coverage = "
							+ resargs[4][1]);
			resargs[4][0] = "coveragedef";
			corcnt = corcnt + 1;
		}
		// System.out.println("grmlcnt  : " + grmlcnt);
		if (grmlcnt == 0) {
			System.out
					.println("range default valuew will be used : range = 0.95");
			resargs[5][0] = "range";
			resargs[5][1] = "0.95";
			corcnt = corcnt + 1;
		}

		if (resargs[1][0].equals("TableExpect")) {
			for (int i = argmns.length - 1; i >= 2; i--) {
				if ((argmns[i - 1].equals("-fpval"))
						|| ((!(isNumeric(argmns[i]))) && (!(argmns[i - 1]
								.charAt(0) == '-')))
						|| (isNumeric(argmns[i - 1]))) {
					jcur = argmns[i];
					inpfl = new File(crarg);
					ftnm = inpfl.getName();
					if (jcur.equals(ftnm)) {
						root_pth = new java.io.File(".").getCanonicalPath();
						root_pth = root_pth + "/";
						jcur = root_pth + ftnm;
					}
					if (testFileExistence(jcur)) {
						System.out.println("input table file exists : " + jcur);
						tbfcnt = tbfcnt + 1;
						if (tbfcnt == 1)
							tbsrce = jcur;
						if (tbfcnt > 1)
							tbsrce = tbsrce + "," + jcur;
					} else if (!(testFileExistence(jcur)))
						System.out.println("table file '" + jcur
								+ "' not found.");
				} else if (argmns[i - 1].charAt(0) == '-')
					i = 1;
			}
		}

		if (!(tbsrce.equals(""))) {
			String[] wds = tbsrce.split(",");
			if (wds.length > 2) {
				System.out
						.println("Exceeded maximum of two input tables.\nExit.");
				Exit();
			}
			if (wds.length == 2) {
				ftbstr = wds[1] + "," + wds[0];
				resargs[1][0] = "tbftwo";
				resargs[1][1] = ftbstr;
				corcnt = corcnt + 1;
			}
			if (wds.length == 1) {
				ftbstr = wds[0];
				resargs[1][0] = "tbfone";
				resargs[1][1] = ftbstr;
				corcnt = corcnt + 1;
			}
		}

		for (int i = 0; i < argmns.length; i++) {
			if (argmns[i].equals("-type")) {
				if (i < argmns.length - 1) {
					type = argmns[i + 1];
					if ((type.equals("pvalue")) || (type.equals("mutrate"))) {
						resargs[0][1] = type;
						System.out.println("graph type OK!");
						corcnt = corcnt + 1;
					} else {
						System.out
								.println("graph type not found.\ngraph type default = pvalue  OK");
						resargs[0][1] = "pvalue";
						corcnt = corcnt + 1;
					}
				}
			}

			if (argmns[i].equals("-overlap")) {
				if (i < argmns.length - 1) {
					overlaps = argmns[i + 1];
					if (isPosUnsignInteger(overlaps)) {
						resargs[3][1] = overlaps;
						resargs[3][0] = "overlap";
						System.out.println("overlap OK!");
						corcnt = corcnt + 1;
					} else {
						System.out
								.println("overlap cutoff value error, default = 0  OK");
						resargs[3][0] = "overlap";
						corcnt = corcnt + 1;
					}
				}
			}

			if (argmns[i].equals("-coverage")) {
				if (i < argmns.length - 1) {
					coverages = argmns[i + 1];

					if (isPosUnsignInteger(coverages)) {
						resargs[4][1] = coverages;
						resargs[4][0] = "coverage";
						System.out.println("coverages OK!");
						corcnt = corcnt + 1;
					} else {
						resargs[4][1] = "100";
						System.out.println("coverage value error, default = "
								+ resargs[4][1] + " OK");
						resargs[4][0] = "coveragedef";
						corcnt = corcnt + 1;
					}
				}
			}

			if (argmns[i].equals("-range")) {
				if (i < argmns.length - 1) {
					ranges = argmns[i + 1];
					// System.out.println(ranges);
					if (isNumeric(ranges)) {
						rangedb = Double.parseDouble(ranges);

						if (rangedb >= 0 && rangedb <= 1) {

							fnrgdb = 1 - rangedb;
							resargs[5][1] = Double.toString(fnrgdb);
							resargs[5][0] = "range";
							System.out.println("range : " + rangedb + " OK!");
							corcnt = corcnt + 1;
						} else {
							System.out
									.println("range value error, default = 0.05  OK");
							resargs[5][0] = "range";
							resargs[5][1] = "0.95";
							corcnt = corcnt + 1;
						}
					}
				}
			}

			if (argmns[i].equals("-o")) {
				if (i < argmns.length - 1) {
					if (argmns[i + 1].charAt(0) != '-') {
						grtst = "";
						grtst = argmns[i + 1];
						outfl = new File(grtst);
						flpth = outfl.getName();
						System.out
								.println("Graph |& frq-pv-tab out : " + flpth);
						// System.out.println("Out!");
						if (grtst.equals(flpth)) {
							root_pth = new java.io.File(".").getCanonicalPath();
							jcur = root_pth + "/" + flpth;
							// System.out.println("In!");
						} else
							jcur = grtst;

						if (testParentDirectoryExistence(jcur)) {
							grphout = jcur;
							resargs[2][1] = grphout;
							corcnt = corcnt + 1;
							System.out.println("graph directory exists OK!");
							System.out.println("graph path = " + grphout);
						} else
							System.out.println("graph directory not found.");
					}
				}
			}
		}

		if (resargs[0][1].equals("")) {
			resargs[0][1] = "pvalue";
			System.out
					.println("graph type not found.\ngraph type default = pvalue  OK");
			corcnt = corcnt + 1;
		}

		if (resargs[2][0].equals("fpval")) {
			System.out
					.println("graph type default = pvalue  OK.\nOnly table(s) will be provided.");
			resargs[0][1] = "pvalue";
		}

		if (resargs[3][0].equals("dfault"))

		{

			System.out.println("graph -overlap default = 0 cutoff OK.");
			resargs[3][1] = "0";
			resargs[3][0] = "overlap";
			corcnt = corcnt + 1;
		}

		if (resargs[5][0].isEmpty()) {
			resargs[5][0] = "range";
			resargs[5][1] = "0.95";
			System.out.println("range value error, default = 0.05  OK");
			corcnt = corcnt + 1;
		}

		System.out.println("Graph args count = " + corcnt);
		System.out.println("-----testArgsforGraphs method End.");
		if (corcnt < 6) {
			System.out.println("corcnt = " + corcnt);
			System.out.println("Graph arguments Error. Exit.");
			Exit();
		}
		System.out.println();
		return resargs;

		// System.out.println("testArgsforGraphs method!");
		//
		// String graphinfo = "graph 		[main command]\n"+
		// "[-type]		[pvalue or mutrate]\n"+ "[-tbdir]	[tables path]\n"+
		// "[-grout]	[graph destination]";
		//
		// "mutrate"; return resargs;
		//
	}
    // test graph (END)
 
	// decide Test (START)
	/**
	 * 
	 * @param argmns
	 * @return
	 * @throws IOException
	 */
	public static String[][] testArgsforDecide(String[] argmns)
			throws IOException {

		// Replace -tbf
		String decideinfo = "\nUsage: Leucippus decide [options] table.file\n\n"
				+ "Options: -o           FILE    results\n"
				+ "         -coverage            minimum site/position coverage [100]\n"
				+ "         -germlineAF  DOUBLE  minimum AF to call variant as germline [0.35]\n"
				+ "         -pvalue      DOUBLE  p-value for somatic call[0.05]\n\n";
		
		//String decideinfo = "\nUsage: leucippus decide [options] table.file\n\n"
		//		+ "Options: -o           FILE    results\n"
		//		+ "         -coverage            minimum site/position coverage [100]\n"
		//		+ "         -germlineAF  DOUBLE  minimum AF to call variant as germline [0.35]\n\n";

		//
		// System.out.println("Start TestArgsforDecide method!");
		// String decideinfo
		// ="\nUsage:	leucippus decide [[option] <value> [option] <value>]\n\n"
		// + "Options:	-tbf			FILE	table\n"
		// + "		-o			FILE	results\n\n";

//		java Leucippus decide -coverage 100 -pvalue 40 -o decidew.tsv /data5/experpath/vasm/vasm/NextGen/niko/MosaicsTablesGraphs021115mdy/i1/Tables/d7/00_i1_d7_0_.tsv

		System.out.println("\ntestArgsforDecide method Start ------------");
		String[][] resargs = new String[5][2];
		resargs[0][0] = "";
		resargs[0][1] = "";
		resargs[1][0] = "outres";
		resargs[1][1] = "";
		resargs[2][0] = "coverage";
		resargs[2][1] = ""; // default =100
		resargs[3][0] = "";
		resargs[3][1] = ""; // default =100
		resargs[4][0] = "pvalue"; // pvalue
		resargs[4][1] = "0.0"; // default =0.05

		String correct_report = "", missed_report = "";
		String cucrinf = "", cumsinf = "";

		String grtst = "";
		String intab = "", outres = "", jcur = "", cvs = "";
		int cvi = 0;
		double gmlnd = 0.0, pvld=0.0;
		String root_pth = "", inpth = "", flpth = "", tbsrce = "", tblpthstr = "", outst = "", output = "";
		String ftnm = "";
		File outfl, inpfl;
		int cn = 0, incnt = 0, outcnt = 0, cvcnt = 0, gmlncnt = 0, pvaluecnt=0;
		// , gprmcnt=0;
		for (int i = 0; i < argmns.length; i++) {
			if (argmns[i] != null) {
				// if (argmns[i].equals("-tbf"))
				// incnt = incnt + 1;
				if (argmns[i].equals("-o"))
					outcnt = outcnt + 1;
				if (argmns[i].equals("-coverage"))
					cvcnt = cvcnt + 1;
				if (argmns[i].equals("-germlineAF"))
					gmlncnt = gmlncnt + 1;
				if (argmns[i].equals("-pvalue"))
					pvaluecnt = pvaluecnt + 1;
				cn = cn + 1;
			}
		}

		if (argmns[argmns.length - 2].charAt(0) == '-') {
			System.out
					.println("The input table must be entered as the last argument.\nInstead identifier for the last argument found.\nExit.");
			Exit();
		} else {
			resargs[0][0] = "intabexp";
		}

		// System.out.println("testArgsforGraphs method  2 ");
		if (outcnt > 1 || cvcnt > 1 || gmlncnt > 1 || pvaluecnt >1) {
			System.out.println("Arguments must be stated once. \n Exit");
			Exit();
		}

		if (outcnt == 0) {
			// System.out.println("arguments for source table file name and output file name are required.\n(coverage is optional)\n Exit");
			System.out
					.println("arguments for output file name is required.\n(coverage and germlineAF are optional)\n Exit");
			Exit();
		}
		int corcnt = 0;

		if (resargs[0][0].equals("intabexp")) {
			jcur = argmns[argmns.length - 1];
			inpfl = new File(jcur);
			inpth = inpfl.getName();
			System.out.println(inpth);
			if (jcur.equals(inpth)) {
				root_pth = new java.io.File(".").getCanonicalPath();
				root_pth = root_pth + "/";
				jcur = root_pth + jcur;
			}
			if (testFileExistence(jcur)) {
				System.out.println("input table file exists : " + jcur);
				resargs[0][1] = jcur;
				corcnt = corcnt + 1;
			} else if (!(testFileExistence(jcur))) {
				System.out
						.println("Table file not found : " + jcur + "\nExit.");
				Exit();
			}
		}

		for (int i = 0; i < cn; i++) {

			jcur = "";
			if (argmns[i].equals("-germlineAF")) {
				gmlnd = 0.0;
				System.out.println(argmns[i]);
				if (i < argmns.length - 1) {
					jcur = argmns[i + 1];
					System.out.println(jcur);
					if (isNumeric(jcur))
						gmlnd = Double.parseDouble(jcur);
					if (gmlnd >= 0) {
						resargs[3][0] = "germlineAF";
						resargs[3][1] = jcur;
						System.out.println("germlineAF = " + jcur + " is OK.");
						corcnt = corcnt + 1;
					} else {
						System.out.println("germlineAF must be double >= 0.");
						resargs[3][0] = "DefgermlineAF";
						resargs[3][1] = "0.35";
						System.out.println("Default germlineAF = 0.35.");
						corcnt = corcnt + 1;
					}
				}
			} // if (argmns[i].equals("-germlineAF");

			jcur = "";
			if (argmns[i].equals("-pvalue")) 
			{
				pvld = 0.0;
				System.out.println(argmns[i]);
				if (i < argmns.length - 1) 
				{
					jcur = argmns[i + 1];
					System.out.println(jcur);
					if (isNumeric(jcur))
					{
						pvld = Double.parseDouble(jcur);
						if((pvld>1) || (pvld<0))
						{
							System.out.println("pvalue = " + jcur + " is not within (0-1).\n");
							System.out.println("Default pvalue = 0.05 will be used");
							resargs[4][0] = "pvaluedf";
							resargs[4][1] = "0.05";
							System.out.println(resargs[4][0] + " = " + resargs[4][1] + " is OK.");
							corcnt = corcnt + 1;
						}
						if ( (pvld > 0) && (pvld < 1))
						{
							resargs[4][0] = "pvalue";
							resargs[4][1] = jcur;
							System.out.println("pvalue = " + jcur + " is OK.");
							corcnt = corcnt + 1;
						}
					}
					else 
					{
						System.out.println("pvalue is not numeric.\n");
						System.out.println("Default pvalue = 0.05 will be used");
						resargs[4][0] = "pvaluedf";
						resargs[4][1] = "0.05";
						System.out.println(resargs[4][0] + " = " + resargs[4][1] + " is OK.");
						corcnt = corcnt + 1;
					}
				}
			} // if (argmns[i].equals("-germlineAF");

			jcur = "";
			if (argmns[i].equals("-o")) {
				if (i < argmns.length - 1) {
					if (argmns[i + 1].charAt(0) != '-') {
						outst = "";
						outst = argmns[i + 1];
						outfl = new File(outst);
						flpth = outfl.getName();
						System.out.println("out file for decide : " + flpth);
						if (outst.equals(flpth)) {
							root_pth = new java.io.File(".").getCanonicalPath();
							jcur = root_pth + "/" + flpth;
							System.out.println("In!");
						} else
							jcur = outst;

						if (testParentDirectoryExistence(jcur)) {
							output = jcur;
							resargs[1][1] = output;
							resargs[1][0] = "outres";
							corcnt = corcnt + 1;
							System.out
									.println("output parent directory exists OK!");
							System.out.println("output path = " + output);
						} else
							System.out
									.println("output parent directory not found.");
						// the counter will not be increased if
						// the previous else statement is correct, and it will cause
						// termination of the program (line close to the end of this method).

					}
				}
			}

			jcur = "";
			if (argmns[i].equals("-coverage")) {
				cvi = -1;
				System.out.println(argmns[i]);
				if (i < argmns.length - 1) {
					jcur = argmns[i + 1];
					System.out.println(jcur);
					if (isNumeric(jcur))
						cvi = Integer.parseInt(jcur);
					if (cvi >= 0) {
						resargs[2][1] = jcur;
						resargs[2][0] = "coverage";
						System.out.println("Coverage = " + jcur + " is OK.");
						corcnt = corcnt + 1;
					} else {
						System.out.println("Coverage must be integer >= 0");
						resargs[2][0] = "defcoverage";
						resargs[2][1] = "100";
						System.out
								.println("Default Coverage will be used.\ncoverage >= 100");
						corcnt = corcnt + 1;
					}
				}
			}
		}

		if (resargs[2][1].length() == 0) {
			resargs[2][0] = "defcoverage";
			resargs[2][1] = "100";
			System.out
					.println("Default value for coverage will be used.\ncoverage = 100");
			corcnt = corcnt + 1;
		}

		if (resargs[3][0].length() == 0) {
			System.out.println("germlineAF not entered.");
			resargs[3][0] = "defgermlineAF";
			resargs[3][1] = "0.35";
			System.out
					.println("Default value for germlinAf will be used.\ngermlineAF >= 0.35");
			System.out.println(resargs[3][0] + " = " + resargs[3][1] + " is OK.");
			corcnt = corcnt + 1;
		}
		
		if (resargs[4][1].equals("0.0"))
 		{
			System.out.println("pvalue not entered.");
			resargs[4][0] = "pvaluedf";
			resargs[4][1] = "0.05";
			System.out.println("Default value for pvalue will be used.\npvalue = 0.05");
			System.out.println(resargs[4][0] + " = " + resargs[4][1] + " is OK.");
			corcnt = corcnt + 1;
		}



		if (corcnt != 5) {
			System.out.println("TestArgsforDecide method found errors.");
			Exit();
		}
		System.out.println("-----------testArgsforDecide method End!\n");
		return resargs;
	}
    // decide Test (END)

// --------------- TEST ARGUMENTS-----------------------------------------\\
// -----------------------------------------------------------------------\\
// ---------------------END-----------------------------------------------\\

	/**
	 * 'CreateGraphs' method accepts all necessary parameters in order to
	 * construct an R script call that is sent to bash command line. It uses
	 * java Process builder to do so.
	 * 
	 * @param mutratepvalue
	 *            : String that holds the type of the graph to be created
	 *            mutation rate or p-value
	 * @param tables_path
	 *            : String for input table(s)
	 * @param graphs_path
	 *            : String for output(prefix) graph (and/or) pvalue-frquency
	 *            taqble
	 * @param rscript_path
	 *            : String path where R is installed
	 * @param bash_path
	 *            : String path for bash
	 * 
	 * @param cfpath
	 *            : String path for configuration file
	 * @param grortb
	 *            : String graph or table or both
	 * @param overlap
	 *            : String that holds a number for the lowest acceptable value
	 *            for x and y in input table file. All excluded positions(in
	 *            table file) are not in the overlaping region.
	 * @param coverages
	 *            : String minimum value for number of expected bases
	 * @param germlineAFs
	 * @throws IOException
	 * @throws InterruptedException
	 **/
	protected static void CreateGraphs(String mutratepvalue,
			String tables_path, String graphs_path, String rscript_path,
			String bash_path, String cfpath, String grortb, String overlap,
			String coverages, String germlineAFs) throws IOException,
			InterruptedException {
	    // CreateGraphs method (START)
		// grortb : graph or table or both
		System.out.println("\nStart Create Graph method!\n");

		System.out.println("Graph type : " + mutratepvalue + "\ntables_path : "
				+ tables_path + "\ngraphs_path : " + graphs_path
				+ "\nrscript_path : " + rscript_path + "\nbash_path : "
				+ bash_path + "\ncfpath : " + cfpath + "\nOverlap : " + overlap
				+ "\ncoverage : " + coverages + "\nrange : "
				+ (1 - Double.parseDouble(germlineAFs)));

//		ClassLoader loader = Leucippus082815.class.getClassLoader();
		ClassLoader loader = Leucippus.class.getClassLoader();

		URL crpathU = loader.getResource("");
		String crpaths = crpathU.toString();
		String path = crpaths.substring(5, crpaths.length());

		String rscgrpth = "";
		// if(rscript_path.equals("none"))
		//System.out.println("Rscript Path : " + rscript_path);
		//if (rscript_path.isEmpty())
			//rscript_path = "/usr/local/biotools/r/R-3.0.1/bin/Rscript";
		System.out.println("Rscript Path : " + rscript_path);

		if (((mutratepvalue.equals("pvalue")) || (mutratepvalue
				.equals("mutrate"))) && (cfpath.length() > 0))
			rscgrpth = path + "Pvalue_Mutrate.r";

		else if (((mutratepvalue.equals("pvalue")) || (mutratepvalue
				.equals("mutrate"))) && (cfpath.length() < 1))
			rscgrpth = path + "/" + "Pvalue_Mutrate.r";
		// rscgrpth = "Pvalue_Mutrate.r";
		else {
			System.out.println("Cannot find '" + mutratepvalue
					+ "' type graph.\nExit.");
			Exit();
		}

		String rnrscrpt = "src=" + tables_path + "; " + "dst=" + graphs_path
				+ "; tp=" + mutratepvalue + "; gt=" + grortb + "; ov="
				+ overlap + "; cov=" + coverages + "; gla=" + germlineAFs
				+ "; " + rscript_path + " " + rscgrpth
				+ " -s $src -d $dst -t $tp -h $gt -o $ov -c $cov -l $gla";
		System.out.println(rnrscrpt);
		try {
			ProcessBuilder pb = new ProcessBuilder("/bin/bash", "-c", rnrscrpt);
			// Map<String, String> environ = pb.environment();
			Process process = pb.start();
			// Vector<String> vct = new Vector<String>();
			// java.io.OutputStream stdin = process.getOutputStream() ;
			// java.io.InputStream stderr = process.getErrorStream();
			// java.io.InputStream stdout = process.getInputStream();
			// BufferedReader reader = new BufferedReader (new
			// InputStreamReader(stdout));
			// while ((line = reader.readLine()) != null)
			// {
			// vct.add(line+"\n");
			// System.out.println(line);
			// }
			// reader.close();
			process.waitFor();
		} finally {
			System.out
					.println("Process builder end. 'Create Graph(10)' method.\n");
		}		
		System.out.println("End 'Create Graph' method.\n");
	      // CreateGraphs method (END)
	}
    public static void MakeLongReads(BufferedWriter bw, String srd1c2out, BufferedReader br1,
            BufferedReader br2, int lrmlength, int p, int mvvrlp, String separator)
                    // throws IOException, InterruptedException
        {
        //String st = "/data5/experpath/vasm/vasm/NextGen/niko/S1123-02_F_LongReads/R1<R2_1123-02_RR_005.fastq.gz";
        // Create try catch statement : if first input file creates 
        // ending error (not properly closed in its creation)then go 
        // in finaly and if normal passing variable is not updated then
        // go through while loop through all remained reeds in the other 
        // input file(create variable that identifies the current while
        // loop (esxternal (first; input file) or internal(second; 
        // input file.
        // If the second file creates the ending error then in finally 
        // statement go through while loop through all remained reads 
        // of the first file. After all loops are finished in finally
        // or if normal termination occured then close all open files
        // in finally.      
        System.out.println("            Start  of  TestBufferedReaderPass method!");

        boolean storefailed=false;
        String shortread1out=null, shortread2out=null;
        BufferedWriter bw1=null, bw2=null;
        int br1rcnti=0, br2rcnti=0;
        if(!((srd1c2out==null) || srd1c2out.isEmpty()) )
        {
            String[] srd12out = srd1c2out.split(","); 
            storefailed=true;
            if(srd12out.length>1)
            {
                shortread1out=srd12out[0];
                shortread2out=srd12out[1];
            }
        }


        if(!( (shortread1out==null) || (shortread2out==null) ))
        {
            GZIPOutputStream gzipoutr1;
            try {
                gzipoutr1 = new GZIPOutputStream(
                new FileOutputStream(shortread1out));
                bw1 = new BufferedWriter(new OutputStreamWriter(gzipoutr1));
                
            } catch (IOException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
            
            GZIPOutputStream gzipoutr2;
            try {
                gzipoutr2 = new GZIPOutputStream(
                new FileOutputStream(shortread2out));
                bw2 = new BufferedWriter(new OutputStreamWriter(gzipoutr2));
            } catch (IOException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
            
        }
        else
        {
            storefailed=false;
        }


/*
        if( (bw1!=null) && (bw2!=null))
            storefailed=true;
        System.out.println(" --->  " + storefailed);
        if(bw1==null)
            System.out.println("bw1 = null");       
        if(bw2==null)
            System.out.println("bw2 = null");
*/


        long time1 = 0, time2 = 0, timeSpent = 0;
        String line_1 = "", line_2 = "", flname="", name = "", read = "", qual = "", namecnum; 
                                // line_1 for buffered reader 1, line_2 for 2
                                // Current name, read, quality
        String lread_1 = "", lread_2 = "", lr = ""; // long read from buffered
                                // reader 1, long read from
                                // buffered reader 2, type 1
                                // or 2
        // each buffered reader can produce only one long read during one
        // iteration in the loop
        // therefore if participate two buffered readers it is possible for two
        // long reads to be created after only one iteration
        // thus there are two lread_1 and lread_2 variables.
        int counterout_1 = 0, counterin_1 = 0; // read counters; fastq format
        // provides 4 lines for each read counterout_1 and counterout_2 counters 
        // count the number of the reads
        // for a particular bufferreader because they increase their value
        // once in the outer loop when the name of the read is provided.
        // counterin_1 and 2 : they are used in inner loop to retrieve read
        // sequence and quality( only one and three value is used (2 contains
        // '+')

        int counterout_2 = 0, counterin_2 = 0;
        int lrcnt = 0;

        ProcessReads pr = new ProcessReads();
        Readob ro = new Readob();
        Readob robq = new Readob();
        Hashtable<String, Readob> reads = new Hashtable<String, Readob>();
        try {
            calculateBinCoef(lrmlength);
        } catch (IOException | InterruptedException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        System.out.println("Binomial Coefficients have been Restored!");
        if (br2 == null)
            System.out.println("Make Long Reads One input file or Pipe!");

        if (!(br2 == null))
            System.out.println("Make Long Reads Two input files!");

        // The main loop iterates through the first buffered reader that it for
        // surte contains the left reads It is also possible to contain the
        // right reads too
    
        String nmprog="", type="", type1="", type2="";
        time1 = System.nanoTime();
        String name1="", name2="", shortrd_1="", shortrd_2="";
        String inout="";

        try {
            while ((line_1 = br1.readLine()) != null)
            // && (lrcnt< 15501) )
            {
                inout="outout";
                type="";
                ro = new Readob();
                lread_1 = "none";
                lread_2 = "none";
                flname="";
                flname=line_1;
                name = pr.getSimpleName(line_1, separator);
                type = pr.getType(line_1, separator);
                
                        if(type1.equals(""))
                        {
                            type1=type;
                        }
                        if((type2.equals("")) && (!type1.equals(type)))
                        {
                            type2=type;
                        }

                counterin_1 = 1;
                // inout="outout";
                while ((counterin_1 <= 3) && ((line_1 = br1.readLine()) != null)) 
                {
                    inout="outin";
                    if (counterin_1 == 1)
                        read = line_1;
                    if (counterin_1 == 3)
                        qual = line_1;
                    counterin_1 = counterin_1 + 1;
                    // System.out.println(" counterin_1 :  " + counterin_1);
                }
                robq = null;
                br1rcnti=br1rcnti+1;
                ro = buildRead02(flname, name, read, qual);
                robq = reads.get(name);
//          in scenario with two input files :
//              It is possible that the following 'if' never happens if the positions of reads in the 
//          files have direct sequential correspondence. For example the first file contains read
//          name1 as first read. The second file contains read name1 as paired read to previous read
//          the external while loop will constract the read and then it will query the hash table for 
//          paired read. The paired read is in the second file thus the has table will return none.
//          the read from the first file will be pushed in the hash table. In the next while loop
//          when the paired read is eventually constructed the hash table will return the read that 
//          has been pushed to it from the first external while loop.
//          If this senarion is happening all the time the following if statement will be always skipped. 
                if (!(robq == null))
                {
                    reads.remove(name);
                    robq=rvcomp(robq);
                    lread_1 = makeLongRead(ro, robq, mvvrlp);
                    // System.out.println(lread_1);
                    if (!(lread_1.equals("none")))
                    {
                        // System.out.println("Long Read type = " + type);
                        lrcnt = lrcnt + 1;
    
                        if (lrcnt % 500 == 0)
                            System.out.println("LR nm = " + lrcnt);
                        
                        bw.write(lread_1);
                        bw.newLine();
//                      if (lrcnt>10000)
//                      {
//                        bw.close();
//                        Exit();
//                      }
                        
                        
                    }
 
                    else if ( (lread_1.equals("none")) && (storefailed==true) ) 
                    {
                        if (type.equals(type1))
                        {
                            shortrd_1 = ro.getNamefull() + "\n" + ro.getSeq() + "\n+\n" + ro.getQual();
                            // System.out.println(ro.getName() + " " + type);
                            bw1.write(shortrd_1);
                            bw1.newLine();
                            
                            shortrd_2 = robq.getNamefull() + "\n" + robq.getSeq() + "\n+\n" + robq.getQual();
                            bw2.write(shortrd_2);
                            bw2.newLine();
                        }
                        else if (!(type.equals(type1)))
                        {
                            shortrd_2 = ro.getNamefull() + "\n" + ro.getSeq() + "\n+\n" + ro.getQual();
                            // System.out.println(ro.getName() + " " + type);
                            bw2.write(shortrd_2);
                            bw2.newLine();
                            
                            shortrd_1 = robq.getNamefull() + "\n" + robq.getSeq() + "\n+\n" + robq.getQual();
                            bw1.write(shortrd_1);
                            bw1.newLine();
                        }

                    }
                    // System.out.println("!Here ! ++++++++++++++++++ ! Here!");
                }
                // it is possible that the following else statement is happening all times
                    else if (robq == null)
                    reads.put(name, ro);
                //  System.out.println(robq);
                // still inout="outout";
                if (!(br2 == null)) // if second buffered reader is not null(=
                {
                    // System.out.println("!Here ! &&&&&&&&&&&&&&&&& ! Here!");
                    counterin_2 = 1;
                    inout="outinin";
                    while ((counterin_2 <= 3)
                            && ((line_2 = br2.readLine()) != null)) // Every time // inout="outin";
                    {
                        inout="in";
                        if (!(line_2.isEmpty()))
                        {
                            type="";
                            ro = null;
                            robq = null;
                            lread_2 = "none";
                            flname="";
                            flname=line_2;
                            name = pr.getSimpleName(line_2, separator);
                            type = pr.getType(line_2, separator);
                                    if((type2.equals("")) && (!type1.equals(type)))
                                    {
                                        type2=type;
                                    }

                            while ((counterin_2 <= 3)
                                     && ((line_2 = br2.readLine()) != null))
                            {
                                if (counterin_2 == 1)
                                    read = line_2;
                                if (counterin_2 == 3)
                                    qual = line_2;
                                counterin_2 = counterin_2 + 1;
                            }
                            br2rcnti=br2rcnti+1;
                            ro = buildRead02(flname, name, read, qual);
                            robq = reads.get(name);
                                if (robq == null)
                                reads.put(name, ro);
//                      it is possible that the following else statement happens all times
                                else if (!(robq == null))
                            {
                                reads.remove(name);
                                ro=rvcomp(ro);
                                lread_2 = makeLongRead(robq, ro, mvvrlp);
                                if (!(lread_2.equals("none")))
                                {
                                    // System.out.println("Long Read type = " + type);
                                    lrcnt = lrcnt + 1;
                                    if (lrcnt % 500 == 0)
                                        System.out.println(lrcnt);
                                    bw.write(lread_2);
                                    bw.newLine();
                                    //if (lrcnt>10000)
                                    //{
                                    //	bw.close();
                                    //	Exit();
                                    //}
                                }

                                else if ((lread_2.equals("none")) && (storefailed==true) )
                                {

                                    // System.out.println(storefailed);
                                    // shortrd_1 = robq.getName() + " " + type;
                                    if (type.equals(type2))
                                    {
                                        shortrd_2 = ro.getNamefull() + "\n" + ro.getSeq() + "\n+\n" + ro.getQual();
                                        // System.out.println(ro.getName() + " not long  " + type);
                                        bw2.write(shortrd_2);   
                                        bw2.newLine();
                                        shortrd_1 = robq.getNamefull() + "\n" + robq.getSeq() + "\n+\n" + robq.getQual();
                                        bw1.write(shortrd_1);
                                        bw1.newLine();
                                    }
                                    else if (!(type.equals(type2)))
                                    {
                                        shortrd_1 = ro.getNamefull() + "\n" + ro.getSeq() + "\n+\n" + ro.getQual();
                                        // System.out.println(ro.getName() + " not long  " + type);
                                        bw1.write(shortrd_1);   
                                        bw1.newLine();
                                        shortrd_2 = robq.getNamefull() + "\n" + robq.getSeq() + "\n+\n" + robq.getQual();
                                        bw2.write(shortrd_2);
                                        bw2.newLine();
                                    }
                                }
                            }
                        }
                    }
                    inout="outinout";
                }

                //if (lrcnt % 2000 == 0)
                //{
                    //System.out.println(lrcnt);
                    //break;
                //}
            }
        } 
        catch (IOException e) 
        {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        finally
        {
            System.out.println("Exception catched." + inout);
            System.out.println("Number of Unpaired Reads remained in Hash Table : "
                    + reads.size());
            try
            {
            if ( (inout.equals("outout")) || (inout.equals("outin")) || (inout.equals("outinout")) )
            {
            	System.out.println(inout);
                if (!(br2 == null)) // if second buffered reader is not null(=
                {
                    //System.out.println("!Here ! &&&&&&&&&&&&&&&&& ! Here!");
                    counterin_2 = 1;

                        if(!(br2.readLine() == null))
                            System.out.println("NOT NULL");
                        if(reads.isEmpty())
                            System.out.println("Hash Empty!");
                        else
                        	System.out.println("Hash Not Empty!");
                        
                        if(!(reads.isEmpty()))
                        while ((counterin_2 <= 3)
                                && ((line_2 = br2.readLine()) != null) && (!(reads.isEmpty()))) // Every time
                        {
                           // System.out.println("COUNT !");
                            if (!(line_2.isEmpty()))
                            {
                                type="";
                                ro = null;
                                robq = null;
                                lread_2 = "none";
                                flname="";
                                flname=line_2;
                                name = pr.getSimpleName(line_2, separator);
                                type = pr.getType(line_2, separator);
                                if( (type2.equals("")) && (!type1.equals(type)) )
                                    type2=type;

                                while ((counterin_2 <= 3)
                                        && ((line_2 = br2.readLine()) != null))
                                {
                                    if (counterin_2 == 1)
                                        read = line_2;
                                    if (counterin_2 == 3)
                                        qual = line_2;
                                    counterin_2 = counterin_2 + 1;
                                }
                                br2rcnti=br2rcnti+1;
                                ro = buildRead02(flname, name, read, qual);
                                robq = reads.get(name);
                                //if (robq == null)
                                  //  reads.put(name, ro);
//                              it is possible that the following else statement happens all times
                                   // else if (!(robq == null))
                               
                                if (robq != null)
                                {
                                    reads.remove(name);
                                    ro=rvcomp(ro);
                                
                                    lread_2 = makeLongRead(robq, ro, mvvrlp);
                                    if (!(lread_2.equals("none")))
                                    {
                                        // System.out.println("Long Read type = " + type);
                                        lrcnt = lrcnt + 1;
                                       // if (lrcnt % 100 == 0)
                                        System.out.println("LR out 1 = " + lrcnt);
                                        bw.write(lread_2);
                                        bw.newLine();
                                    }

                                    else if ((lread_2.equals("none")) && (storefailed==true) )
                                    {

                                        // System.out.println(storefailed);
                                        // shortrd_1 = robq.getName() + " " + type;
                                        if (type.equals(type2))
                                        {
                                            shortrd_2 = ro.getNamefull() + "\n" + ro.getSeq() + "\n+\n" + ro.getQual();
                                            // System.out.println(ro.getName() + " not long  " + type);
                                            bw2.write(shortrd_2);   
                                            bw2.newLine();
                                            shortrd_1 = robq.getNamefull() + "\n" + robq.getSeq() + "\n+\n" + robq.getQual();
                                            bw1.write(shortrd_1);
                                            bw1.newLine();
                                        }
                                        else if (!(type.equals(type2)))
                                        {
                                            shortrd_1 = ro.getNamefull() + "\n" + ro.getSeq() + "\n+\n" + ro.getQual();
                                            // System.out.println(ro.getName() + " not long  " + type);
                                            bw1.write(shortrd_1);   
                                            bw1.newLine();
                                            shortrd_2 = robq.getNamefull() + "\n" + robq.getSeq() + "\n+\n" + robq.getQual();
                                            bw2.write(shortrd_2);
                                            bw2.newLine();
                                        }
                                    }
                                }
                            }
                            counterin_2 = 1;
                        }
                      }  
                   } 
            
            if (inout.equals("in"))
            {
            	System.out.println("In!");
            	if(br1 != null)
            	{	 
                while ((line_1 = br1.readLine()) != null)
                    // && (lrcnt< 15501) )
                {
                        inout="outout";
                        type="";
                        ro = new Readob();
                        lread_1 = "none";
                        lread_2 = "none";
                        flname="";
                        flname=line_1;
                        name = pr.getSimpleName(line_1, separator);
                        type = pr.getType(line_1, separator);
                        
                        if(type1.equals(""))
                        {
                            type1=type;
                        }
                        if((type2.equals("")) && (!type1.equals(type)))
                        {
                           type2=type;
                        }
                        
                        counterin_1 = 1;
                        // inout="outout";
                        while ((counterin_1 <= 3) && ((line_1 = br1.readLine()) != null)) 
                        {
                            inout="outin";
                            if (counterin_1 == 1)
                                read = line_1;
                            if (counterin_1 == 3)
                                qual = line_1;
                            counterin_1 = counterin_1 + 1;
                            // System.out.println(" counterin_1 :  " + counterin_1);
                        }
                        robq = null;
                        br1rcnti=br1rcnti+1;
                        ro = buildRead02(flname, name, read, qual);
                        robq = reads.get(name);
                        
                        if (robq != null)
                        {
                            reads.remove(name);
                            robq=rvcomp(robq);
                            lread_1 = makeLongRead(ro, robq, mvvrlp);
                            // System.out.println(lread_1);
                            if (!(lread_1.equals("none")))
                            {
                                // System.out.println("Long Read type = " + type);
                                lrcnt = lrcnt + 1;
            
                                if (lrcnt % 100 == 0)
                                    System.out.println("LR nm = " + lrcnt);
                                
                                bw.write(lread_1);
                                bw.newLine();      
                            }
         
                            else if ( (lread_1.equals("none")) && (storefailed==true) ) 
                            {
                                if (type.equals(type1))
                                {
                                    shortrd_1 = ro.getNamefull() + "\n" + ro.getSeq() + "\n+\n" + ro.getQual();
                                    // System.out.println(ro.getName() + " " + type);
                                    bw1.write(shortrd_1);
                                    bw1.newLine();
                                    
                                    shortrd_2 = robq.getNamefull() + "\n" + robq.getSeq() + "\n+\n" + robq.getQual();
                                    bw2.write(shortrd_2);
                                    bw2.newLine();
                                }
                                else if (!(type.equals(type1)))
                                {
                                    shortrd_2 = ro.getNamefull() + "\n" + ro.getSeq() + "\n+\n" + ro.getQual();
                                    // System.out.println(ro.getName() + " " + type);
                                    bw2.write(shortrd_2);
                                    bw2.newLine();
                                    
                                    shortrd_1 = robq.getNamefull() + "\n" + robq.getSeq() + "\n+\n" + robq.getQual();
                                    bw1.write(shortrd_1);
                                    bw1.newLine();
                                }
                            }
                            // System.out.println("!Here ! ++++++++++++++++++ ! Here!");
                        }
                        // it is possible that the following else statement is happening all times
                        counterin_1 = 1;
                }
            }	
            }
            }
            catch (IOException e) 
            {
                        // TODO Auto-generated catch block
                        e.printStackTrace();
            }
            finally
            {
                System.out.println("Second Finally works!");
                System.out.println("total reads 1 = " + br1rcnti);
                System.out.println("total reads 2 = " + br2rcnti);
            }
                
            if (!(br1 == null))
                try {
                    br1.close();
                } catch (IOException e) {
                    // TODO Auto-generated catch block
                    System.out.println("2");
                    e.printStackTrace();
                }
            if (!(br2 == null))
                try {
                    br2.close();
                } catch (IOException e1) {
                    // TODO Auto-generated catch block
                    System.out.println("3");
                    e1.printStackTrace();
                }
            try {
                bw.close();
            } catch (IOException e) {
                // TODO Auto-generated catch block
                System.out.println("4");
                e.printStackTrace();
            }

            try {
                if(bw1!=null)
                bw1.close();
            } catch (IOException e) {
                // TODO Auto-generated catch block
                System.out.println("5");
                e.printStackTrace();
            }
            try {
                if(bw2!=null)
                bw2.close();
            } catch (IOException e) {
                // TODO Auto-generated catch block
                System.out.println("6");
                e.printStackTrace();
            }
            System.out.println("At least files were closed well.");
        }


        if(storefailed==false)
        {
            System.out.println(" Delete both short reads failed to create long reads.?");
        }

        System.out.println("Number of Unpaired Reads remained in Hash Table : "
                + reads.size());

        time2 = System.nanoTime();
        timeSpent = time2 - time1;
        double db = 0.0;
        db = timeSpent / 1000000000;
        db = db / 60;
        System.out.println("Time Spent : " + db + " minutes.");
        System.out.println("            End  of  MakeLongReads method.");
    }



    public static void MakeLongReadsBack110915(BufferedWriter bw, String srd1c2out, BufferedReader br1,
            BufferedReader br2, int lrmlength, int p, int mvvrlp, String separator)throws InterruptedException, IOException
                    // throws IOException, InterruptedException
        {
        // Create try catch statement : if first input file creates 
        // ending error (not properly closed in its creation)then go 
        // in finaly and if normal passing variable is not updated then
        // go through while loop through all remained reeds in the other 
        // input file(create variable that identifies the current while
        // loop (esxternal (first; input file) or internal(second; 
        // input file.
        // If the second file creates the ending error then in finally 
        // statement go through while loop through all remained reads 
        // of the first file. After all loops are finished in finally
        // or if normal termination occured then close all open files
        // in finally.      
        System.out.println("            Start  of  TestBufferedReaderPass method!");
	Vector<String> remainedinhash = new Vector<String>();
	String rmanpth = "path";
        

	boolean storefailed=false;
        String shortread1out=null, shortread2out=null;
        BufferedWriter bw1=null, bw2=null;
        int br1rcnti=0, br2rcnti=0;
        if(!((srd1c2out==null) || srd1c2out.isEmpty()) )
        {
            String[] srd12out = srd1c2out.split(","); 
            storefailed=true;
            if(srd12out.length>1)
            {
                shortread1out=srd12out[0];
                shortread2out=srd12out[1];
            }
        }


        if(!( (shortread1out==null) || (shortread2out==null) ))
        {
            GZIPOutputStream gzipoutr1;
            try {
                gzipoutr1 = new GZIPOutputStream(
                new FileOutputStream(shortread1out));
                bw1 = new BufferedWriter(new OutputStreamWriter(gzipoutr1));
                
            } catch (IOException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
            
            GZIPOutputStream gzipoutr2;
            try {
                gzipoutr2 = new GZIPOutputStream(
                new FileOutputStream(shortread2out));
                bw2 = new BufferedWriter(new OutputStreamWriter(gzipoutr2));
            } catch (IOException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
            
        }
        else
        {
            storefailed=false;
        }


/*
        if( (bw1!=null) && (bw2!=null))
            storefailed=true;
        System.out.println(" --->  " + storefailed);
        if(bw1==null)
            System.out.println("bw1 = null");       
        if(bw2==null)
            System.out.println("bw2 = null");
*/


        long time1 = 0, time2 = 0, timeSpent = 0;
        String line="", line_1 = "", line_2 = "", flname="", name = "", read = "", qual = "", namecnum; 
                                // line_1 for buffered reader 1, line_2 for 2
                                // Current name, read, quality
        String lread_1 = "", lread_2 = "", lr = ""; // long read from buffered
                                // reader 1, long read from
                                // buffered reader 2, type 1
                                // or 2
        // each buffered reader can produce only one long read during one
        // iteration in the loop
        // therefore if participate two buffered readers it is possible for two
        // long reads to be created after only one iteration
        // thus there are two lread_1 and lread_2 variables.
        int counterout_1 = 0, counterin_1 = 0; // read counters; fastq format
        // provides 4 lines for each read counterout_1 and counterout_2 counters 
        // count the number of the reads
        // for a particular bufferreader because they increase their value
        // once in the outer loop when the name of the read is provided.
        // counterin_1 and 2 : they are used in inner loop to retrieve read
        // sequence and quality( only one and three value is used (2 contains
        // '+')
	String remanpath="";
        int counterout_2 = 0, counterin_2 = 0;
        int lrcnt = 0;

        ProcessReads pr = new ProcessReads();
        Readob ro = new Readob();
        Readob robq = new Readob();
	Readob rdb = new Readob();
        Hashtable<String, Readob> reads = new Hashtable<String, Readob>();
        try {
            calculateBinCoef(lrmlength);
        } catch (IOException | InterruptedException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        System.out.println("Binomial Coefficients have been Restored!");
        if (br2 == null)
            System.out.println("Make Long Reads One input file or Pipe!");

        if (!(br2 == null))
            System.out.println("Make Long Reads Two input files!");

        // The main loop iterates through the first buffered reader that it for
        // surte contains the left reads It is also possible to contain the
        // right reads too
    
        String nmprog="", type="", type1="", type2="";
        time1 = System.nanoTime();
        String name1="", name2="", shortrd_1="", shortrd_2="";
        String inout="";

        try {
            while ((line_1 = br1.readLine()) != null)
            // && (lrcnt< 15501) )
            {
                inout="outout";
                type="";
                ro = new Readob();
                lread_1 = "none";
                lread_2 = "none";
                flname="";
                flname=line_1;
                name = pr.getSimpleName(line_1, separator);
                type = pr.getType(line_1, separator);
                
                        if(type1.equals(""))
                        {
                            type1=type;
                        }
                        if((type2.equals("")) && (!type1.equals(type)))
                        {
                            type2=type;
                        }

                counterin_1 = 1;
                // inout="outout";
                while ((counterin_1 <= 3) && ((line_1 = br1.readLine()) != null)) 
                {
                    inout="outin";
                    if (counterin_1 == 1)
                        read = line_1;
                    if (counterin_1 == 3)
                        qual = line_1;
                    counterin_1 = counterin_1 + 1;
                    // System.out.println(" counterin_1 :  " + counterin_1);
                }
                robq = null;
                br1rcnti=br1rcnti+1;
                ro = buildRead02(flname, name, read, qual);
                robq = reads.get(name);
//          in scenario with two input files :
//              It is possible that the following 'if' never happens if the positions of reads in the 
//          files have direct sequential correspondence. For example the first file contains read
//          name1 as first read. The second file contains read name1 as paired read to previous read
//          the external while loop will constract the read and then it will query the hash table for 
//          paired read. The paired read is in the second file thus the has table will return none.
//          the read from the first file will be pushed in the hash table. In the next while loop
//          when the paired read is eventually constructed the hash table will return the read that 
//          has been pushed to it from the first external while loop.
//          If this senarion is happening all the time the following if statement will be always skipped. 
                if (!(robq == null))
                {
                    reads.remove(name);
                    robq=rvcomp(robq);
                    lread_1 = makeLongRead(ro, robq, mvvrlp);
                    // System.out.println(lread_1);
                    if (!(lread_1.equals("none")))
                    {
                        // System.out.println("Long Read type = " + type);
                        lrcnt = lrcnt + 1;
    
                        if (lrcnt % 100 == 0)
                            System.out.println("LR nm = " + lrcnt);
                        
                        bw.write(lread_1);
                        bw.newLine();
                        
                        
                        
                    }
 
                    else if ( (lread_1.equals("none")) && (storefailed==true) ) 
                    {
                        if (type.equals(type1))
                        {
                            shortrd_1 = ro.getNamefull() + "\n" + ro.getSeq() + "\n+\n" + ro.getQual();
                            // System.out.println(ro.getName() + " " + type);
                            bw1.write(shortrd_1);
                            bw1.newLine();
                            
                            shortrd_2 = robq.getNamefull() + "\n" + robq.getSeq() + "\n+\n" + robq.getQual();
                            bw2.write(shortrd_2);
                            bw2.newLine();
                        }
                        else if (!(type.equals(type1)))
                        {
                            shortrd_2 = ro.getNamefull() + "\n" + ro.getSeq() + "\n+\n" + ro.getQual();
                            // System.out.println(ro.getName() + " " + type);
                            bw2.write(shortrd_2);
                            bw2.newLine();
                            
                            shortrd_1 = robq.getNamefull() + "\n" + robq.getSeq() + "\n+\n" + robq.getQual();
                            bw1.write(shortrd_1);
                            bw1.newLine();
                        }

                    }
                    // System.out.println("!Here ! ++++++++++++++++++ ! Here!");
                }
                // it is possible that the following else statement is happening all times
                    else if (robq == null)
                    reads.put(name, ro);
                //  System.out.println(robq);
                // still inout="outout";
                if (!(br2 == null)) // if second buffered reader is not null(=
                {
                    // System.out.println("!Here ! &&&&&&&&&&&&&&&&& ! Here!");
                    counterin_2 = 1;
                    inout="outinin";
                    while ((counterin_2 <= 3)
                            && ((line_2 = br2.readLine()) != null)) // Every time // inout="outin";
                    {
                        inout="in";
                        if (!(line_2.isEmpty()))
                        {
                            type="";
                            ro = null;
                            robq = null;
                            lread_2 = "none";
                            flname="";
                            flname=line_2;
                            name = pr.getSimpleName(line_2, separator);
                            type = pr.getType(line_2, separator);
                                    if((type2.equals("")) && (!type1.equals(type)))
                                    {
                                        type2=type;
                                    }

                            while ((counterin_2 <= 3)
                                     && ((line_2 = br2.readLine()) != null))
                            {
                                if (counterin_2 == 1)
                                    read = line_2;
                                if (counterin_2 == 3)
                                    qual = line_2;
                                counterin_2 = counterin_2 + 1;
                            }
                            br2rcnti=br2rcnti+1;
                            ro = buildRead02(flname, name, read, qual);
                            robq = reads.get(name);
                                if (robq == null)
                                reads.put(name, ro);
//                      it is possible that the following else statement happens all times
                                else if (!(robq == null))
                            {
                                reads.remove(name);
                                ro=rvcomp(ro);
                                lread_2 = makeLongRead(robq, ro, mvvrlp);
                                if (!(lread_2.equals("none")))
                                {
                                    // System.out.println("Long Read type = " + type);
                                    lrcnt = lrcnt + 1;
                                    if (lrcnt % 100 == 0)
                                        System.out.println(lrcnt);
                                    bw.write(lread_2);
                                    bw.newLine();
                                }

                                else if ((lread_2.equals("none")) && (storefailed==true) )
                                {

                                    // System.out.println(storefailed);
                                    // shortrd_1 = robq.getName() + " " + type;
                                    if (type.equals(type2))
                                    {
                                        shortrd_2 = ro.getNamefull() + "\n" + ro.getSeq() + "\n+\n" + ro.getQual();
                                        // System.out.println(ro.getName() + " not long  " + type);
                                        bw2.write(shortrd_2);   
                                        bw2.newLine();
                                        shortrd_1 = robq.getNamefull() + "\n" + robq.getSeq() + "\n+\n" + robq.getQual();
                                        bw1.write(shortrd_1);
                                        bw1.newLine();
                                    }
                                    else if (!(type.equals(type2)))
                                    {
                                        shortrd_1 = ro.getNamefull() + "\n" + ro.getSeq() + "\n+\n" + ro.getQual();
                                        // System.out.println(ro.getName() + " not long  " + type);
                                        bw1.write(shortrd_1);   
                                        bw1.newLine();
                                        shortrd_2 = robq.getNamefull() + "\n" + robq.getSeq() + "\n+\n" + robq.getQual();
                                        bw2.write(shortrd_2);
                                        bw2.newLine();
                                    }
                                }
                            }
                        }
                    }
                    inout="outinout";
                }

                //if (lrcnt % 2000 == 0)
                //{
                    //System.out.println(lrcnt);
                    //break;
                //}
            }
        } 
        catch (IOException e) 
        {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        finally
        {
            System.out.println("Exception catched." + inout);
            System.out.println("Number of Unpaired Reads remained in Hash Table : "
                    + reads.size());
            try
            {
            if ( (inout.equals("outout")) || (inout.equals("outin")) || (inout.equals("outinout")) )
            {
            	System.out.println(inout);
                if (!(br2 == null)) // if second buffered reader is not null(=
                {
                    System.out.println("!Here ! &&&&&&&&&&&&&&&&& ! Here!");
                    counterin_2 = 1;

                        if(!(br2.readLine() == null))
                            System.out.println("NOT NULL");
                        if(reads.isEmpty())
                            System.out.println("Hash Empty!");
                        else
                        	System.out.println("Hash Not Empty!");
                        
                        if(!(reads.isEmpty()))
                        while ((counterin_2 <= 3)
                                && ((line_2 = br2.readLine()) != null) && (!(reads.isEmpty()))) // Every time
                        {
                           // System.out.println("COUNT !");
                            if (!(line_2.isEmpty()))
                            {
                                type="";
                                ro = null;
                                robq = null;
                                lread_2 = "none";
                                flname="";
                                flname=line_2;
                                name = pr.getSimpleName(line_2, separator);
                                type = pr.getType(line_2, separator);
                                if( (type2.equals("")) && (!type1.equals(type)) )
                                    type2=type;

                                while ((counterin_2 <= 3)
                                        && ((line_2 = br2.readLine()) != null))
                                {
                                    if (counterin_2 == 1)
                                        read = line_2;
                                    if (counterin_2 == 3)
                                        qual = line_2;
                                    counterin_2 = counterin_2 + 1;
                                }
                                br2rcnti=br2rcnti+1;
                                ro = buildRead02(flname, name, read, qual);
                                robq = reads.get(name);
                                //if (robq == null)
                                  //  reads.put(name, ro);
//                              it is possible that the following else statement happens all times
                                   // else if (!(robq == null))
                               
                                if (robq != null)
                                {
                                    reads.remove(name);
                                    ro=rvcomp(ro);
                                
                                    lread_2 = makeLongRead(robq, ro, mvvrlp);
                                    if (!(lread_2.equals("none")))
                                    {
                                        // System.out.println("Long Read type = " + type);
                                        lrcnt = lrcnt + 1;
                                       // if (lrcnt % 100 == 0)
                                        System.out.println("LR out 1 = " + lrcnt);
                                        bw.write(lread_2);
                                        bw.newLine();
                                    }

                                    else if ((lread_2.equals("none")) && (storefailed==true) )
                                    {

                                        // System.out.println(storefailed);
                                        // shortrd_1 = robq.getName() + " " + type;
                                        if (type.equals(type2))
                                        {
                                            shortrd_2 = ro.getNamefull() + "\n" + ro.getSeq() + "\n+\n" + ro.getQual();
                                            // System.out.println(ro.getName() + " not long  " + type);
                                            bw2.write(shortrd_2);   
                                            bw2.newLine();
                                            shortrd_1 = robq.getNamefull() + "\n" + robq.getSeq() + "\n+\n" + robq.getQual();
                                            bw1.write(shortrd_1);
                                            bw1.newLine();
                                        }
                                        else if (!(type.equals(type2)))
                                        {
                                            shortrd_1 = ro.getNamefull() + "\n" + ro.getSeq() + "\n+\n" + ro.getQual();
                                            // System.out.println(ro.getName() + " not long  " + type);
                                            bw1.write(shortrd_1);   
                                            bw1.newLine();
                                            shortrd_2 = robq.getNamefull() + "\n" + robq.getSeq() + "\n+\n" + robq.getQual();
                                            bw2.write(shortrd_2);
                                            bw2.newLine();
                                        }
                                    }
                                }
                            }
                            counterin_2 = 1;
                        }
                      }  
                   } 
            
            if (inout.equals("in"))
            {
            	System.out.println("In!");
            	if(br1 != null)
            	{	 
                while ((line_1 = br1.readLine()) != null)
                    // && (lrcnt< 15501) )
                {
                        inout="outout";
                        type="";
                        ro = new Readob();
                        lread_1 = "none";
                        lread_2 = "none";
                        flname="";
                        flname=line_1;
                        name = pr.getSimpleName(line_1, separator);
                        type = pr.getType(line_1, separator);
                        
                        if(type1.equals(""))
                        {
                            type1=type;
                        }
                        if((type2.equals("")) && (!type1.equals(type)))
                        {
                           type2=type;
                        }
                        
                        counterin_1 = 1;
                        // inout="outout";
                        while ((counterin_1 <= 3) && ((line_1 = br1.readLine()) != null)) 
                        {
                            inout="outin";
                            if (counterin_1 == 1)
                                read = line_1;
                            if (counterin_1 == 3)
                                qual = line_1;
                            counterin_1 = counterin_1 + 1;
                            // System.out.println(" counterin_1 :  " + counterin_1);
                        }
                        robq = null;
                        br1rcnti=br1rcnti+1;
                        ro = buildRead02(flname, name, read, qual);
                        robq = reads.get(name);
                        
                        if (robq != null)
                        {
                            reads.remove(name);
                            robq=rvcomp(robq);
                            lread_1 = makeLongRead(ro, robq, mvvrlp);
                            // System.out.println(lread_1);
                            if (!(lread_1.equals("none")))
                            {
                                // System.out.println("Long Read type = " + type);
                                lrcnt = lrcnt + 1;
            
                                if (lrcnt % 100 == 0)
                                    System.out.println("LR nm = " + lrcnt);
                                
                                bw.write(lread_1);
                                bw.newLine();      
                            }
         
                            else if ( (lread_1.equals("none")) && (storefailed==true) ) 
                            {
                                if (type.equals(type1))
                                {
                                    shortrd_1 = ro.getNamefull() + "\n" + ro.getSeq() + "\n+\n" + ro.getQual();
                                    // System.out.println(ro.getName() + " " + type);
                                    bw1.write(shortrd_1);
                                    bw1.newLine();
                                    
                                    shortrd_2 = robq.getNamefull() + "\n" + robq.getSeq() + "\n+\n" + robq.getQual();
                                    bw2.write(shortrd_2);
                                    bw2.newLine();
                                }
                                else if (!(type.equals(type1)))
                                {
                                    shortrd_2 = ro.getNamefull() + "\n" + ro.getSeq() + "\n+\n" + ro.getQual();
                                    // System.out.println(ro.getName() + " " + type);
                                    bw2.write(shortrd_2);
                                    bw2.newLine();
                                    
                                    shortrd_1 = robq.getNamefull() + "\n" + robq.getSeq() + "\n+\n" + robq.getQual();
                                    bw1.write(shortrd_1);
                                    bw1.newLine();
                                }
                            }
                            // System.out.println("!Here ! ++++++++++++++++++ ! Here!");
                        }
                        // it is possible that the following else statement is happening all times
                        counterin_1 = 1;
                }
            }	
            }
            }
            catch (IOException e) 
            {
                        // TODO Auto-generated catch block
                        e.printStackTrace();
            }
            finally
            {
                System.out.println("Second Finally works!");
                System.out.println("total reads 1 = " + br1rcnti);
                System.out.println("total reads 2 = " + br2rcnti);


		//  name = "", read = "", qual
    		Enumeration dn;
		if(reads.size()>1)
		{
    			dn = reads.elements();
			while (dn.hasMoreElements())
			{
				rdb=(Readob)dn.nextElement();
				name=rdb.getNamefull();
				read=rdb.getSeq();
				qual=rdb.getQual();
				line = name + "\n"+ read +"\n+\n" + qual;
				remainedinhash.add(line);
			}
		}
		if(remainedinhash.size()>0)
		{
			writeToFile(rmanpth, remainedinhash);
		}


            }
                
            if (!(br1 == null))
                try {
                    br1.close();
                } catch (IOException e) {
                    // TODO Auto-generated catch block
                    System.out.println("2");
                    e.printStackTrace();
                }
            if (!(br2 == null))
                try {
                    br2.close();
                } catch (IOException e1) {
                    // TODO Auto-generated catch block
                    System.out.println("3");
                    e1.printStackTrace();
                }
            try {
                bw.close();
            } catch (IOException e) {
                // TODO Auto-generated catch block
                System.out.println("4");
                e.printStackTrace();
            }

            try {
                if(bw1!=null)
                bw1.close();
            } catch (IOException e) {
                // TODO Auto-generated catch block
                System.out.println("5");
                e.printStackTrace();
            }
            try {
                if(bw2!=null)
                bw2.close();
            } catch (IOException e) {
                // TODO Auto-generated catch block
                System.out.println("6");
                e.printStackTrace();
            }
            System.out.println("At least files were closed well.");
        }


        if(storefailed==false)
        {
            System.out.println(" Delete both short reads failed to create long reads.?");
        }

        System.out.println("Number of Unpaired Reads remained in Hash Table : "
                + reads.size());

        time2 = System.nanoTime();
        timeSpent = time2 - time1;
        double db = 0.0;
        db = timeSpent / 1000000000;
        db = db / 60;
        System.out.println("Time Spent : " + db + " minutes.");
        System.out.println("            End  of  MakeLongReads method.");
    }

	public static void MakeLongReadsBack082(BufferedWriter bw, String srd1c2out, BufferedReader br1,
			BufferedReader br2, int lrmlength, int p, int mvvrlp, String separator)
			        //throws IOException, InterruptedException{
		// Create try catch statement : if first input file creates 
		// ending error (not properly closed in its creation)then go 
		// in finaly and if normal passing variable is not updated then
		// go through while loop through all remained reeds in the other 
		// input file(create variable that identifies the current while
		// loop (esxternal (first; input file) or internal(second; 
		// input file.
		// If the second file creates the ending error then in finally 
		// statement go through while loop through all remained reads 
		// of the first file. After all loops are finished in finally
		// or if normal termination occured then close all open files
		// in finally.	
	    {
		System.out.println("			Start  of  TestBufferedReaderPass method!");

		boolean storefailed=false;
		String shortread1out=null, shortread2out=null;
		BufferedWriter bw1=null, bw2=null;
		
		if(!((srd1c2out==null) || srd1c2out.isEmpty()) )
		{
			String[] srd12out = srd1c2out.split(","); 
			storefailed=true;
			if(srd12out.length>1)
			{
				shortread1out=srd12out[0];
				shortread2out=srd12out[1];
			}
		}

		try 
		{
		    if(!( (shortread1out==null) || (shortread2out==null) ))
		    {
			GZIPOutputStream gzipoutr1;
            
                gzipoutr1 = new GZIPOutputStream(
                new FileOutputStream(shortread1out));
  
			bw1 = new BufferedWriter(new OutputStreamWriter(gzipoutr1));
			GZIPOutputStream gzipoutr2 = new GZIPOutputStream(
			new FileOutputStream(shortread2out));
			bw2 = new BufferedWriter(new OutputStreamWriter(gzipoutr2));
		    }
		    else
		    {
		        storefailed=false;
		    }
        } 
		catch (IOException e)
		{
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
/*
		if( (bw1!=null) && (bw2!=null))
			storefailed=true;
		System.out.println(" --->  " + storefailed);
		if(bw1==null)
			System.out.println("bw1 = null");		
		if(bw2==null)
			System.out.println("bw2 = null");

*/
		long time1 = 0, time2 = 0, timeSpent = 0;
		String line_1 = "", line_2 = "", flname="", name = "", read = "", qual = "", namecnum; 
								// line_1 for buffered reader 1, line_2 for 2
								// Current name, read, quality
		String lread_1 = "", lread_2 = "", lr = "";	// long read from buffered
								// reader 1, long read from
								// buffered reader 2, type 1
								// or 2
		// each buffered reader can produce only one long read during one
		// iteration in the loop
		// therefore if participate two buffered readers it is possible for two
		// long reads to be created after only one iteration
		// thus there are two lread_1 and lread_2 variables.
		int counterout_1 = 0, counterin_1 = 0; // read counters; fastq format
		// provides 4 lines for each read counterout_1 and counterout_2 counters 
		// count the number of the reads
		// for a particular bufferreader because they increase their value
		// once in the outer loop when the name of the read is provided.
		// counterin_1 and 2 : they are used in inner loop to retrieve read
		// sequence and quality( only one and three value is used (2 contains
		// '+')

		int counterout_2 = 0, counterin_2 = 0;
		int lrcnt = 0;

		ProcessReads pr = new ProcessReads();
		Readob ro = new Readob();
		Readob robq = new Readob();
		Hashtable<String, Readob> reads = new Hashtable<String, Readob>();
		
		try {
            calculateBinCoef(lrmlength);
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (InterruptedException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
		System.out.println("Binomial Coefficients have been Restored!");
		if (br2 == null)
			System.out.println("Make Long Reads One input file or Pipe!");

		if (!(br2 == null))
			System.out.println("Make Long Reads Two input files!");

		// The main loop iterates through the first buffered reader that it for
		// surte contains the left reads It is also possible to contain the
		// right reads too
	
		String nmprog="", type="", type1="", type2="";
		time1 = System.nanoTime();
		String name1="", name2="", shortrd_1="", shortrd_2="";
		String whileinout="";

		try {
            while ((line_1 = br1.readLine()) != null)
            // && (lrcnt< 15501) )
            {
                whileinout="out";
            	type="";
            	ro = new Readob();
            	lread_1 = "none";
            	lread_2 = "none";
            	flname="";
            	flname=line_1;
            	name = pr.getSimpleName(line_1, separator);
            	type = pr.getType(line_1, separator);
            	
                		if(type1.equals(""))
                		{
                    		type1=type;
                		}
                		if((type2.equals("")) && (!type1.equals(type)))
                		{
                    		type2=type;
                		}

            	counterin_1 = 1;

            	while ((counterin_1 <= 3) && ((line_1 = br1.readLine()) != null)) 
            	{
            		if (counterin_1 == 1)
            			read = line_1;
            		if (counterin_1 == 3)
            			qual = line_1;
            		counterin_1 = counterin_1 + 1;
            		// System.out.println(" counterin_1 :  " + counterin_1);
            	}
            	robq = null;
            	ro = buildRead02(flname, name, read, qual);
            	robq = reads.get(name);
//			in scenario with two input files :
//		        It is possible that the following 'if' never happens if the positions of reads in the 
//			files have direct sequential correspondence. For example the first file contains read
//			name1 as first read. The second file contains read name1 as paired read to previous read
//			the external while loop will constract the read and then it will query the hash table for 
//			paired read. The paired read is in the second file thus the has table will return none.
//			the read from the first file will be pushed in the hash table. In the next while loop
//			when the paired read is eventually constructed the hash table will return the read that 
//			has been pushed to it from the first external while loop.
//			If this senarion is happening all the time the following if statement will be always skipped. 
            	if (!(robq == null))
            	{
            		reads.remove(name);
            		robq=rvcomp(robq);
            		lread_1 = makeLongRead(ro, robq, mvvrlp);
            		// System.out.println(lread_1);
            		if (!(lread_1.equals("none")))
            		{
            			// System.out.println("Long Read type = " + type);
            			lrcnt = lrcnt + 1;
            			bw.write(lread_1);
            			bw.newLine();
            			// System.out.println("!Here ! ---- ! Here!");
            		}
 
            		else if ( (lread_1.equals("none")) && (storefailed==true) ) 
            		{
            			if (type.equals(type1))
            			{
            				shortrd_1 = ro.getNamefull() + "\n" + ro.getSeq() + "\n+\n" + ro.getQual();
            				// System.out.println(ro.getName() + " " + type);
            				bw1.write(shortrd_1);
            				bw1.newLine();
            				
            				shortrd_2 = robq.getNamefull() + "\n" + robq.getSeq() + "\n+\n" + robq.getQual();
            				bw2.write(shortrd_2);
            				bw2.newLine();
            			}
            			else if (!(type.equals(type1)))
            			{
            				shortrd_2 = ro.getNamefull() + "\n" + ro.getSeq() + "\n+\n" + ro.getQual();
            				// System.out.println(ro.getName() + " " + type);
            				bw2.write(shortrd_2);
            				bw2.newLine();
            				
            				shortrd_1 = robq.getNamefull() + "\n" + robq.getSeq() + "\n+\n" + robq.getQual();
            				bw1.write(shortrd_1);
            				bw1.newLine();
            			}

            		}
            		// System.out.println("!Here ! ++++++++++++++++++ ! Here!");
            	}
            	// it is possible that the following else statement is happening all times
                	else if (robq == null)
            		reads.put(name, ro);
            	//	System.out.println(robq);
            	if (!(br2 == null)) // if second buffered reader is not null(=
            	{
            	   
            		// System.out.println("!Here ! &&&&&&&&&&&&&&&&& ! Here!");
            		counterin_2 = 1;
            		while ((counterin_2 <= 3)
            				&& ((line_2 = br2.readLine()) != null)) // Every time
            		{
            		    whileinout="in";
            			if (!(line_2.isEmpty()))
            			{
            				type="";
            				ro = null;
            				robq = null;
            				lread_2 = "none";
            				flname="";
            				flname=line_2;
            				name = pr.getSimpleName(line_2, separator);
            				type = pr.getType(line_2, separator);
                					if((type2.equals("")) && (!type1.equals(type)))
                					{
                    					type2=type;
                					}

            				while ((counterin_2 <= 3)
            						&& ((line_2 = br2.readLine()) != null))
            				{
            					if (counterin_2 == 1)
            						read = line_2;
            					if (counterin_2 == 3)
            						qual = line_2;
            					counterin_2 = counterin_2 + 1;
            				}

            				ro = buildRead02(flname, name, read, qual);
            				robq = reads.get(name);
            		        	if (robq == null)
            					reads.put(name, ro);
// 						it is possible that the following else statement happens all times
            		        	else if (!(robq == null))
            				{
            					reads.remove(name);
            					ro=rvcomp(ro);
            					lread_2 = makeLongRead(robq, ro, mvvrlp);
            					if (!(lread_2.equals("none")))
            					{
            						// System.out.println("Long Read type = " + type);
            						lrcnt = lrcnt + 1;
            						if (lrcnt % 1000 == 0)
            							System.out.println(lrcnt);
            						bw.write(lread_2);
            						bw.newLine();
            					}

            					else if ((lread_2.equals("none")) && (storefailed==true) )
            					{

            						// System.out.println(storefailed);
            						// shortrd_1 = robq.getName() + " " + type;
            						if (type.equals(type2))
            						{
            							shortrd_2 = ro.getNamefull() + "\n" + ro.getSeq() + "\n+\n" + ro.getQual();
            							// System.out.println(ro.getName() + " not long  " + type);
            							bw2.write(shortrd_2);	
            							bw2.newLine();
            							shortrd_1 = robq.getNamefull() + "\n" + robq.getSeq() + "\n+\n" + robq.getQual();
            							bw1.write(shortrd_1);
            							bw1.newLine();
            						}
            						else if (!(type.equals(type2)))
            						{
            							shortrd_1 = ro.getNamefull() + "\n" + ro.getSeq() + "\n+\n" + ro.getQual();
            							// System.out.println(ro.getName() + " not long  " + type);
            							bw1.write(shortrd_1);	
            							bw1.newLine();
            							shortrd_2 = robq.getNamefull() + "\n" + robq.getSeq() + "\n+\n" + robq.getQual();
            							bw2.write(shortrd_2);
            							bw2.newLine();
            						}
            					}
            				}
            			}
            		}
            	}

            	//if (lrcnt % 2000 == 0)
            	//{
            		//System.out.println(lrcnt);
            	    //break;
            	//}
            }
        }
	catch (IOException e) 
	{
 //           if(something) then do something
            e.printStackTrace();
        }
	finally
	{
		if(whileinout.equals("in"))
		{ 
	            System.out.println("Exception catched." + whileinout);
	            {
	                if (!(br2 == null)) // if second buffered reader is not null(=
	                {
	                    // System.out.println("!Here ! &&&&&&&&&&&&&&&&& ! Here!");
	                    counterin_2 = 1;
	                    try {
	                        while ((counterin_2 <= 3)
	                                && ((line_2 = br2.readLine()) != null) && (!(reads.isEmpty()))) // Every time
	                        {
	                            if (!(line_2.isEmpty()))
	                            {
	                                type="";
	                                ro = null;
	                                robq = null;
	                                lread_2 = "none";
	                                flname="";
	                                flname=line_2;
	                                name = pr.getSimpleName(line_2, separator);
	                                type = pr.getType(line_2, separator);
	                                        if((type2.equals("")) && (!type1.equals(type)))
	                                        {
	                                            type2=type;
	                                        }

	                                while ((counterin_2 <= 3)
	                                        && ((line_2 = br2.readLine()) != null))
	                                {
	                                    if (counterin_2 == 1)
	                                        read = line_2;
	                                    if (counterin_2 == 3)
	                                        qual = line_2;
	                                    counterin_2 = counterin_2 + 1;
	                                }

	                                ro = buildRead02(flname, name, read, qual);
	                                robq = reads.get(name);
	                                    if (robq == null)
	                                    reads.put(name, ro);
//	                      it is possible that the following else statement happens all times
	                                    else if (!(robq == null))
	                                {
	                                    reads.remove(name);
	                                    ro=rvcomp(ro);
	                                    lread_2 = makeLongRead(robq, ro, mvvrlp);
	                                    if (!(lread_2.equals("none")))
	                                    {
	                                        // System.out.println("Long Read type = " + type);
	                                        lrcnt = lrcnt + 1;
	                                        if (lrcnt % 100 == 0)
	                                            System.out.println(lrcnt);
	                                        bw.write(lread_2);
	                                        bw.newLine();
	                                    }

	                                    else if ((lread_2.equals("none")) && (storefailed==true) )
	                                    {

	                                        // System.out.println(storefailed);
	                                        // shortrd_1 = robq.getName() + " " + type;
	                                        if (type.equals(type2))
	                                        {
	                                            shortrd_2 = ro.getNamefull() + "\n" + ro.getSeq() + "\n+\n" + ro.getQual();
	                                            // System.out.println(ro.getName() + " not long  " + type);
	                                            bw2.write(shortrd_2);   
	                                            bw2.newLine();
	                                            shortrd_1 = robq.getNamefull() + "\n" + robq.getSeq() + "\n+\n" + robq.getQual();
	                                            bw1.write(shortrd_1);
	                                            bw1.newLine();
	                                        }
	                                        else if (!(type.equals(type2)))
	                                        {
	                                            shortrd_1 = ro.getNamefull() + "\n" + ro.getSeq() + "\n+\n" + ro.getQual();
	                                            // System.out.println(ro.getName() + " not long  " + type);
	                                            bw1.write(shortrd_1);   
	                                            bw1.newLine();
	                                            shortrd_2 = robq.getNamefull() + "\n" + robq.getSeq() + "\n+\n" + robq.getQual();
	                                            bw2.write(shortrd_2);
	                                            bw2.newLine();
	                                        }
	                                    }
	                                }
	                            }
	            
	                        }
	                    } 
	                    catch (IOException e)
	                    {
	                        e.printStackTrace();
	                    }
	                    finally
	                    {
	                        System.out.println("Second Finally works!");
	                        try
	                        {
	                            
	                            if (!(br1 == null))
	                                br1.close();
	                            if (!(br2 == null))
	                                br2.close();
	                            bw.close();
	                            if (!(bw1 == null))
	                                bw1.close();
                                if (!(bw2 == null))
                                    bw2.close();
                                
                                System.out.println("Files were closed correctly!");
	                        }
	                        catch (IOException e)
	                        {
	                            e.printStackTrace();
	                        }
	                        
	                    }
	                }
	                
	            }
		        
		  }
		  else if(whileinout.equals("out"))
		  {
		        System.out.println("out was the event here!");
		  }
		}
		
		System.out.println("It went through here!");
		try
		{
		    if (!(br1 == null))
               		 br1.close();
		    if (!(br2 == null))
		        br2.close();

		    if (!(bw == null))
	            	bw.close();    
		    if (!(bw1 == null))
		        bw1.close();
		    if (!(bw2 == null))
		        bw2.close();
        	} 
		catch (IOException e) 
		{
            		e.printStackTrace();
        	}
	
		if(storefailed==false)
		{
			System.out.println(" Delete both short reads failed to create long reads.");
		}

		System.out.println("Number of Unpaired Reads remained in Hash Table : "
				+ reads.size());

		time2 = System.nanoTime();
		timeSpent = time2 - time1;
		double db = 0.0;
		db = timeSpent / 1000000000;
		db = db / 60;
		System.out.println("Time Spent : " + db + " minutes.");
		System.out.println("			End  of  MakeLongReads method.");
	}


	
	/*
	 * public static void example_of_two_BR(BufferedWriter bw, BufferedReader
	 * br1, BufferedReader br2, int lrmlength, int p) throws IOException { while
	 * ((line1 = br1.readLine()) != null) { // Extract values from line1 and
	 * extra lines ro1=new Readob(name1,read1,qual1); robq=reads.get(name1);
	 * if(robq != null) { reads.remove(name1); lreadlr = makeLongRead(ro1,
	 * robq); } else // if (robq == null) reads.put(name1, ro1);
	 * 
	 * if (br2 != null && (line2 = br2.readLine()) != null) { // Extract values
	 * from line2 ro2=new Readob(name2,read2,qual2); robq=reads.get(name2); if
	 * (robq != null) { reads.remove(name2); lreadrl = makeLongRead(robq, ro2);
	 * } else // if (robq == null) reads.put(name2, ro2); } } // End of While
	 * loop
	 * 
	 * while (br2 != null && (line2 = br2.readLine()) != null) { // Extract
	 * values from line2 and extra lines ro2=new Readob(name2,read2,qual2);
	 * robq=reads.get(name2); if (robq != null) { reads.remove(name2); lreadrl =
	 * makeLongRead(robq, ro2); } else // if (robq == null) reads.put(name2,
	 * ro2); } }
	 * 
	 * public static void example_of_two_BR(BufferedWriter bw, BufferedReader
	 * br1, BufferedReader br2, int lrmlength, int p) throws IOException {
	 * String line="", name="", read="", qual="", namecnum; String line1="",
	 * name1="", read1="", qual1=""; String line2="", name2="", read2="",
	 * qual2=""; String lreadlr="", lreadrl="", lr=""; int counterout=0,
	 * counterin=0, lhcnt=0; int longreascnt=0;
	 * 
	 * Readob ro = new Readob(); Readob robq = new Readob();
	 * 
	 * while ((line1 = br1.readLine()) != null) { //Extract values from line1
	 * and extra lines counterout=counterout+1; if(counterout%500==0)
	 * System.out.println(counterout); ro=null; robq=null; lreadrl="none"; //
	 * System.out.println("READ1 counterout = " + counterout + "  " + line1);
	 * namecnum = pr.getNameCoNum(line); String[] ws1 = namecnum.split(",");
	 * name=ws1[0]; lr=ws1[1]; counterin=1;
	 * 
	 * while((counterin<=3) && ((line = in.readLine())!= null)) { if
	 * (counterin==1) { read=line; } if (counterin==3) { qual=line; }
	 * counterin=counterin+1; } ro=new Readob(name,read,qual);
	 * robq=reads.get(name); if(lr.equals("1")) lreadrl = makeLongRead(ro,
	 * robq);
	 * 
	 * if (br2 != null) {
	 */


	// MAKE FRAG (START)
	/**
	 * "MakeLongReads" method creates long reads from short reads. To do so it
	 * constructs read objects and inserts/deletes them to/from a hash table
	 * using read name as keys. In order to create a long read from two sorted
	 * paired reads it calls "makelongread" method. "makelongread" method slides
	 * the reads against each other; then it finds the best matcing overlapping window
	 * using binomial propability distribution. The process is described in "mekelongread"
	 * method.
	 * 
	 * "MakeLongReads" method accepts a buff writer, two buff readers, and the maximmum 
	 * long read length. "br2" can be null (in the case where both reads(1 and 2 are in "br1"). 
	 * If "br2" is not null then two input files have been provited. Piped input also will work.
	 * 
	 * @param bw : BufferedWriter : Each time a succesfull long read have been constructed, its
	 * information is written to bw
	 *                             
	 * @param br1 : BufferedReader holds read-1 lines in fastq forlmat 
	 * @param br2 : BufferedReader holds read-2 lines in fastq forlmat 
	 * @param lrmlength : integer for long read length (user entered or default(500)) 
	 * @param p integer as boolean (print long read to screen or store it regularly)
	 * @throws IOException
	 * @throws InterruptedException
	 * 
	 * prerequisites : at least one Buffered reader must not be null
	 * read length must be at least 51 base paires (default value = 500)
	 */

	public static void MakeLongReadsBB(BufferedWriter bw, String srd1c2out, BufferedReader br1,
			BufferedReader br2, int lrmlength, int p, int mvvrlp, String separator) throws IOException, InterruptedException{
		// Create try catch statement : if first input file creates 
		// ending error (not properly closed in its creation)then go 
		// in finaly and if normal passing variable is not updated then
		// go through while loop through all remained reeds in the other 
		// input file(create variable that identifies the current while
		// loop (esxternal (first; input file) or internal(second; 
		// input file.
		// If the second file creates the ending error then in finally 
		// statement go through while loop through all remained reads 
		// of the first file. After all loops are finished in finally
		// or if normal termination occured then close all open files
		// in finally.		
		System.out.println("			Start  of  TestBufferedReaderPass method!");

		boolean storefailed=false;
		String shortread1out=null, shortread2out=null;
		BufferedWriter bw1=null, bw2=null;
		
		if(!((srd1c2out==null) || srd1c2out.isEmpty()) )
		{
			String[] srd12out = srd1c2out.split(","); 
			storefailed=true;
			if(srd12out.length>1)
			{
				shortread1out=srd12out[0];
				shortread2out=srd12out[1];
			}
		}


		if(!( (shortread1out==null) || (shortread2out==null) ))
		{
			GZIPOutputStream gzipoutr1 = new GZIPOutputStream(
			new FileOutputStream(shortread1out));
			bw1 = new BufferedWriter(new OutputStreamWriter(gzipoutr1));
			GZIPOutputStream gzipoutr2 = new GZIPOutputStream(
			new FileOutputStream(shortread2out));
			bw2 = new BufferedWriter(new OutputStreamWriter(gzipoutr2));
		}
		else
		{
			storefailed=false;
		}

/*
		if( (bw1!=null) && (bw2!=null))
			storefailed=true;
		System.out.println(" --->  " + storefailed);
		if(bw1==null)
			System.out.println("bw1 = null");		
		if(bw2==null)
			System.out.println("bw2 = null");

*/
		long time1 = 0, time2 = 0, timeSpent = 0;
		String line_1 = "", line_2 = "", flname="", name = "", read = "", qual = "", namecnum; 
								// line_1 for buffered reader 1, line_2 for 2
								// Current name, read, quality
		String lread_1 = "", lread_2 = "", lr = "";	// long read from buffered
								// reader 1, long read from
								// buffered reader 2, type 1
								// or 2
		// each buffered reader can produce only one long read during one
		// iteration in the loop
		// therefore if participate two buffered readers it is possible for two
		// long reads to be created after only one iteration
		// thus there are two lread_1 and lread_2 variables.
		int counterout_1 = 0, counterin_1 = 0; // read counters; fastq format
		// provides 4 lines for each read counterout_1 and counterout_2 counters 
		// count the number of the reads
		// for a particular bufferreader because they increase their value
		// once in the outer loop when the name of the read is provided.
		// counterin_1 and 2 : they are used in inner loop to retrieve read
		// sequence and quality( only one and three value is used (2 contains
		// '+')

		int counterout_2 = 0, counterin_2 = 0;
		int lrcnt = 0;

		ProcessReads pr = new ProcessReads();
		Readob ro = new Readob();
		Readob robq = new Readob();
		Hashtable<String, Readob> reads = new Hashtable<String, Readob>();
		calculateBinCoef(lrmlength);
		System.out.println("Binomial Coefficients have been Restored!");
		if (br2 == null)
			System.out.println("Make Long Reads One input file or Pipe!");

		if (!(br2 == null))
			System.out.println("Make Long Reads Two input files!");

		// The main loop iterates through the first buffered reader that it for
		// surte contains the left reads It is also possible to contain the
		// right reads too
	
		String nmprog="", type="", type1="", type2="";
		time1 = System.nanoTime();
		String name1="", name2="", shortrd_1="", shortrd_2="";


		while ((line_1 = br1.readLine()) != null)
		// && (lrcnt< 15501) )
		{
			type="";
			ro = new Readob();
			lread_1 = "none";
			lread_2 = "none";
			flname="";
			flname=line_1;
			name = pr.getSimpleName(line_1, separator);
			type = pr.getType(line_1, separator);
			
            		if(type1.equals(""))
            		{
                		type1=type;
            		}
            		if((type2.equals("")) && (!type1.equals(type)))
            		{
                		type2=type;
            		}

			counterin_1 = 1;

			while ((counterin_1 <= 3) && ((line_1 = br1.readLine()) != null)) 
			{
				if (counterin_1 == 1)
					read = line_1;
				if (counterin_1 == 3)
					qual = line_1;
				counterin_1 = counterin_1 + 1;
				// System.out.println(" counterin_1 :  " + counterin_1);
			}
			robq = null;
			ro = buildRead02(flname, name, read, qual);
			robq = reads.get(name);
//			in scenario with two input files :
//		        It is possible that the following 'if' never happens if the positions of reads in the 
//			files have direct sequential correspondence. For example the first file contains read
//			name1 as first read. The second file contains read name1 as paired read to previous read
//			the external while loop will constract the read and then it will query the hash table for 
//			paired read. The paired read is in the second file thus the has table will return none.
//			the read from the first file will be pushed in the hash table. In the next while loop
//			when the paired read is eventually constructed the hash table will return the read that 
//			has been pushed to it from the first external while loop.
//			If this senarion is happening all the time the following if statement will be always skipped. 
			if (!(robq == null))
			{
				reads.remove(name);
				robq=rvcomp(robq);
				lread_1 = makeLongRead(ro, robq, mvvrlp);
				// System.out.println(lread_1);
				if (!(lread_1.equals("none")))
				{
					// System.out.println("Long Read type = " + type);
					lrcnt = lrcnt + 1;
					bw.write(lread_1);
					bw.newLine();
					// System.out.println("!Here ! ---- ! Here!");
				}
 
				else if ( (lread_1.equals("none")) && (storefailed==true) ) 
				{
					if (type.equals(type1))
					{
						shortrd_1 = ro.getNamefull() + "\n" + ro.getSeq() + "\n+\n" + ro.getQual();
						// System.out.println(ro.getName() + " " + type);
						bw1.write(shortrd_1);
						bw1.newLine();
						
						shortrd_2 = robq.getNamefull() + "\n" + robq.getSeq() + "\n+\n" + robq.getQual();
						bw2.write(shortrd_2);
						bw2.newLine();
					}
					else if (!(type.equals(type1)))
					{
						shortrd_2 = ro.getNamefull() + "\n" + ro.getSeq() + "\n+\n" + ro.getQual();
						// System.out.println(ro.getName() + " " + type);
						bw2.write(shortrd_2);
						bw2.newLine();
						
						shortrd_1 = robq.getNamefull() + "\n" + robq.getSeq() + "\n+\n" + robq.getQual();
						bw1.write(shortrd_1);
						bw1.newLine();
					}

				}
				// System.out.println("!Here ! ++++++++++++++++++ ! Here!");
			}
			// it is possible that the following else statement is happening all times
	        	else if (robq == null)
				reads.put(name, ro);
			//	System.out.println(robq);
			if (!(br2 == null)) // if second buffered reader is not null(=
			{
				// System.out.println("!Here ! &&&&&&&&&&&&&&&&& ! Here!");
				counterin_2 = 1;
				while ((counterin_2 <= 3)
						&& ((line_2 = br2.readLine()) != null)) // Every time
				{
					if (!(line_2.isEmpty()))
					{
						type="";
						ro = null;
						robq = null;
						lread_2 = "none";
						flname="";
						flname=line_2;
						name = pr.getSimpleName(line_2, separator);
						type = pr.getType(line_2, separator);
            					if((type2.equals("")) && (!type1.equals(type)))
            					{
                					type2=type;
            					}

						while ((counterin_2 <= 3)
								&& ((line_2 = br2.readLine()) != null))
						{
							if (counterin_2 == 1)
								read = line_2;
							if (counterin_2 == 3)
								qual = line_2;
							counterin_2 = counterin_2 + 1;
						}

						ro = buildRead02(flname, name, read, qual);
						robq = reads.get(name);
				        	if (robq == null)
							reads.put(name, ro);
// 						it is possible that the following else statement happens all times
				        	else if (!(robq == null))
						{
							reads.remove(name);
							ro=rvcomp(ro);
							lread_2 = makeLongRead(robq, ro, mvvrlp);
							if (!(lread_2.equals("none")))
							{
								// System.out.println("Long Read type = " + type);
								lrcnt = lrcnt + 1;
								if (lrcnt % 1000 == 0)
									System.out.println(lrcnt);
								bw.write(lread_2);
								bw.newLine();
							}

							else if ((lread_2.equals("none")) && (storefailed==true) )
							{

								// System.out.println(storefailed);
								// shortrd_1 = robq.getName() + " " + type;
								if (type.equals(type2))
								{
									shortrd_2 = ro.getNamefull() + "\n" + ro.getSeq() + "\n+\n" + ro.getQual();
									// System.out.println(ro.getName() + " not long  " + type);
									bw2.write(shortrd_2);	
									bw2.newLine();
									shortrd_1 = robq.getNamefull() + "\n" + robq.getSeq() + "\n+\n" + robq.getQual();
									bw1.write(shortrd_1);
									bw1.newLine();
								}
								else if (!(type.equals(type2)))
								{
									shortrd_1 = ro.getNamefull() + "\n" + ro.getSeq() + "\n+\n" + ro.getQual();
									// System.out.println(ro.getName() + " not long  " + type);
									bw1.write(shortrd_1);	
									bw1.newLine();
									shortrd_2 = robq.getNamefull() + "\n" + robq.getSeq() + "\n+\n" + robq.getQual();
									bw2.write(shortrd_2);
									bw2.newLine();
								}
							}
						}
					}
				}
			}
			//if (lrcnt % 2000 == 0)
			//{
				//System.out.println(lrcnt);
			    //break;
			//}
		}

		if (!(br1 == null))
			br1.close();
		if (!(br2 == null))
			br2.close();
		bw.close();

		bw1.close();
		bw2.close();

		if(storefailed==false)
		{
			System.out.println(" Delete both short reads failed to create long reads.");
		}

		System.out.println("Number of Unpaired Reads remained in Hash Table : "
				+ reads.size());

		time2 = System.nanoTime();
		timeSpent = time2 - time1;
		double db = 0.0;
		db = timeSpent / 1000000000;
		db = db / 60;
		System.out.println("Time Spent : " + db + " minutes.");
		System.out.println("			End  of  MakeLongReads method.");
	}

	public static Readob rvcomp(Readob rorc)
    {
    	ProcessReads pr = null;
    	rorc.setSeq(pr.reverseComplement(rorc.getSeq()));
        rorc.setQual(pr.reverse( rorc.getQual()));
     return rorc;
    }
    	
  
    public static Readob  buildRead02(String fname, String name, String read, String qual)
    {
        ProcessReads pr = new ProcessReads();
        Readob result = new Readob();
        result = new Readob(fname, name, read, qual);
        return result;
    }


/**
 * "makeLongRead" method accepts two parameters Readob (Read object). 
 * The reads passed are paired reads. The sequence in one of them have been reverse 
 * complemented The sequences of the reads can be directly slided one to another
 * from both sides kepping the direction unchenged(a. slide from first read end on
 * second read start, b. slide from second read end on first read start)
 *  
 * @param read1 Readob Object : Read (read name read sequence, read base quality)
 * @param read2 Readob Object : Read (read name read sequence, read base quality)
 * @return
 */
	private static String makeLongRead(Readob read1, Readob read2, int minoverlap) {
	    // makeLongRead method (START)
		// System.out.println("'frags' Make Long Read Method!");
//		09/28/2015 In frag specify minimum overlap (user defined, 
//		default >=50 
	    int  minovrlp = 50;
	    if (minoverlap > 0) minovrlp = minoverlap;
		
		Readob lrd = new Readob();
		String leftread = "", rightCIread = "", 
		        leftqual = "", rightIqual = "", 
		        name1 = "", name2="";
		leftread = read1.getSeq();
		rightCIread = read2.getSeq();
		leftqual = read1.getQual();
		rightIqual = read2.getQual();
		name1 = read1.getName();
		name2 = read2.getName();

		char[] lrchar = leftread.toCharArray();
		char[] rrchar = rightCIread.toCharArray(); // right already complemented
													// and inverted read
		char[] lqchar = leftqual.toCharArray();
		char[] rqchar = rightIqual.toCharArray(); // right already inverted
													// quality
		int minn = Math.min(lrchar.length, rrchar.length);
		int n = 0, k = 0, m = 0, l = 0, r = 0, nreal = 0;
		Vector<String> cmpreads = new Vector<String>();
		Vector<Integer> matcheslr = new Vector<Integer>();
		Vector<Integer> matchesrl = new Vector<Integer>();

		matcheslr.add(0); // mutches start from 1
		int fncnt = -1, bal = 0, indexlri = 0, indexrli = 0, rindexi = 0, lrres_n = 0, rlres_n = 0, tstq = 0;
		cmpreads = new Vector<String>();
		String lefread = "", rigread = "", lefqual = "", riqual = "";
		double db = 0, minlr = 0, minrl = 0, gd = 0, quotlr = 0.0, quotrl = 0.0;
		String qualine = "", namline = "", prdname = "";
		int longreadln = 0;
		char ot = 's'; // overlapping type : s = straight, i = inverse;
		// Vector<String> binresults = new Vector<String>();
		// Line 34
		Vector<Double> rplr = new Vector<Double>(); // initialize to zero :
													// straight
		// l----------l
		// <- r--------r

		Vector<Double> rprl = new Vector<Double>(); // initialize to zero :
													// inverse
		// r-------------r
		// l---------l ->
		String index = "", qindex = "";
		String longread = "", longqual = "";

		String lgread = "";
		Vector<String> fstq = new Vector<String>();
		Vector<Character> bqres = new Vector<Character>();
		for (n = 0; n < minn; n++) {
			nreal = n + 1;
			lefread = "";
			rigread = "";
			m = 0;
			r = -1; // r = right counter
			// l = left counter
			for (l = lrchar.length - n - 1; l < lrchar.length; l++) {
				r = r + 1; // first itteration r = 0;
				// lefread=lefread + lrchar[l];
				// rigread=rigread+rrchar[r];
				if (lrchar[l] == rrchar[r])
// && lrchar[l]) != 'N' && lrchar[l] != 'n') 
						m++;
			}
			matcheslr.add(m);
			// cmpreads.add("(" + (n+1) + ", " + m + ") " + lefread + "_" +
			// rigread);
		}

		n = 0;
		k = 0;
		m = 0;
		l = 0;
		r = 0;
		nreal = 0;
		matchesrl.add(0);
		for (n = 0; n < minn; n++) {
			nreal = n + 1;
			lefread = "";
			rigread = "";
			m = 0;
			l = -1; // l = left counter
			// r = right counter
			for (r = rrchar.length - n - 1; r < rrchar.length; r++) {
				l = l + 1; // first itteration r = 0;
				// lefread=lefread + lrchar[l];
				// rigread=rigread+rrchar[r];
				// System.out.println("l = " + l + "    r = " + r);
				if (lrchar[l] == rrchar[r])
// && lrchar[l]))!= 'N' && lrchar[l] != 'n') 
m++;
			}
			matchesrl.add(m);
			// cmpreads.add("(" + (n+1) + ", " + m + ") " + lefread + "_" +
			// rigread);
			// System.out.println("n = " + (n+1) + ",  m = " + m);
		}

		rplr = new Vector<Double>(); // initialize to empty
		db = 0;
		rplr.add(db);

		// Use the index and the index value of matces to retrieve the
		// coefficient.
		// by both cocmparing matching calculatuions
		for (int o = 1; o < matcheslr.size(); o++) // index zero element is not
													// used
		// for(int o=1; o<20; o++)
		{
			gd = 0;
			// System.out.println(matcheslr.get(o));
			gd = binomCDF(o, o - matcheslr.get(o));
			rplr.add(gd);
			// System.out.println("bin(" + (o) + ", " + matches.get(o) +
			// ")  results = " + g + "\n" + cmpreads.get(o-1));
		}

		rprl = new Vector<Double>();
		db = 0;
		rprl.add(db);
		//
		for (int o = 1; o < matchesrl.size(); o++) // index zero element is not
													// used
		{
			gd = 0;
			// System.out.println(matches.get(o));
			gd = binomCDF(o, o - matchesrl.get(o));
			rprl.add(gd);
			// System.out.println("bin(" + (o) + ", " + matches.get(o) +
			// ")  results = " + g + "\n" + cmpreads.get(o-1));
		}

		// Find the minimum coefficient for standard comparison
		// left-end-right-start
		minlr = 1;
		indexlri = 0;
		// for(int i=1; i<rplr.size(); i++)
		for (int il = 1; il < minn + 1; il++) {
			// -----------------------------------------------------------------------------------------------------------------
			// ---------------------12/10/2014 corrected <1/3 to
			// <1/4--------------------------------------------------------------------------------------------
			// -----------------------------------------------------------------------------------------------------------------
			// if((rplr.get(il)<minlr) && (
			// (il-matcheslr.get(il)<((double)(il))/3)) && (il>=50) )
			if ((rplr.get(il) < minlr)
					&& ((il - matcheslr.get(il) <= ((double) (il)) / 4))
					&& (il >= minovrlp))
			{
				minlr = rplr.get(il);
				indexlri = il;
			}
		}

		lrres_n = matcheslr.get(indexlri);

		quotlr = 0.0;

		// Find the minimum coefficient for the other comparison right-end
		// left-start
		minrl = 1;
		indexrli = 0;

		for (int ir = 1; ir < minn; ir++) {
			// -----------------------------------------------------------------------------------------------------------------
			// ---------------------12/10/2014 corrected <1/3 to
			// <1/4--------------------------------------------------------------------------------------------
			// -----------------------------------------------------------------------------------------------------------------
			// if((rprl.get(ir)<minrl) && (
			// (ir-matchesrl.get(ir)<((double)(ir))/3)) && (ir>=50) )
			if ((rprl.get(ir) < minrl)
					&& ((ir - matchesrl.get(ir) <= ((double) (ir)) / 4))
					&& (ir >= minovrlp))
			// -----------------------------------------------------------------------------------------------------------------
			// ---------------------12/10/2014 corrected <1/3 to
			// <1/4--------------------------------------------------------------------------------------------
			// -----------------------------------------------------------------------------------------------------------------

			{
				minrl = rprl.get(ir);
				indexrli = ir;
			}
		}

		rlres_n = matchesrl.get(indexrli);
		if (minlr < minrl)
			ot = 's';
		if (minrl < minlr)
			ot = 'i';
		if ((minrl == minlr) && !(minrl == 1))
			ot = 's';
		if ((minrl == minlr) && (minrl == 1))
			ot = 'n';

		if (matcheslr.get(indexlri) > 0) // first element is zero and not used,
											// therefore the index corresponds
											// to the actual n
			quotlr = (double) ((double) (matcheslr.get(indexlri)))
					/ ((double) (indexlri));// = m/n (as double)
		else
			quotlr = 0;

		if (matchesrl.get(indexrli) > 0) // first element is zero and not used,
											// therefore the index corresponds
											// to the actual n
			quotrl = (double) ((double) (matchesrl.get(indexrli)))
					/ ((double) (indexrli));// = m/n (as double)
		else
			quotrl = 0;
		// CalcBetterQuality(char lb, char rb, char leftq, char rightq)
		// Standard creation of long read
		if (ot == 's') {
			longread = "";
			longqual = "";
			longread = leftread.substring(0, leftread.length() - indexlri); // read
			longqual = leftqual.substring(0, leftqual.length() - indexlri); // quality
			rindexi = -1;
			for (int b = leftread.length() - indexlri; b < leftread.length(); b++) {
				rindexi = rindexi + 1;
				// bqres = CalcBetterQualityEf(lrchar[b], rrchar[rindexi],
				// lqchar[b], rqchar[rindexi]);
				bqres = ChooseBaseCalculateQuality(lrchar[b], rrchar[rindexi],
						lqchar[b], rqchar[rindexi]);

				longread = longread + bqres.get(0);
				longqual = longqual + bqres.get(1);
			}
			for (int c = indexlri; c < rightCIread.length(); c++) {
				longread = longread + rrchar[c];
				longqual = longqual + rqchar[c];
			}
			longreadln = longread.length();
		} else if (ot == 'i') {
			longread = "";
			longqual = "";
			r = rrchar.length - indexrli - 1; // test this
			for (l = 0; l < indexrli; l++) {
				r = r + 1;
				// bqres = CalcBetterQualityEf(lrchar[l], rrchar[r], lqchar[l],
				// rqchar[r]);
				bqres = ChooseBaseCalculateQuality(lrchar[l], rrchar[r],
						lqchar[l], rqchar[r]);

				longread = longread + bqres.get(0);
				longqual = longqual + bqres.get(1);
			}
			longreadln = longread.length();
		}
		// No suitable mutch found
		// else if(ot == 'n')
		// {
		// longread="_";
		// longqual="_";
		// }

		if (ot == 's') {
			// resultsline = index +"\tmin = " + minlr + "\t" + "bin(" +
			// (indexlri) + ", " + matcheslr.get(indexlri) + ")\t" + quotlr +
			// "\t" + ot + "\t" + prdname + "\t" + longreadln +"\t" + longread +
			// "\t" + longqual;
			// lgread fstq
			String[] wds = name1.split(" ");
			longread = ProcessReads.reverseComplement(longread);
			longqual = ProcessReads.reverse(longqual);
			lgread = "@" + wds[0] + "\n" + longread + "\n+\n" + longqual;

		}
		if (ot == 'i') {
			// resultsline = index +"\tmin = " + minrl + "\t" + "bin(" +
			// (indexrli) + ", " + matchesrl.get(indexrli) + ")\t" + quotrl +
			// "\t" + ot + "\t" + prdname + "\t" + longreadln +"\t" + longread +
			// "\t" + longqual;
			// lgread fstq
			String[] wds = name1.split(" ");
			longread = ProcessReads.reverseComplement(longread);
			longqual = ProcessReads.reverse(longqual);
			lgread = "@" + wds[0] + "\n" + longread + "\n+\n" + longqual;
		}
		// if(ot == 'n')
		// resultsline = index +"\tmin = " + minlr + "\t" + "bin(" + (indexlri)
		// + ", " + matcheslr.get(indexlri) + ")\t" + quotlr + "\t" + ot + "\t"
		// + prdname + "\t" + longreadln +"\t" + longread + "\t" + longqual;

		if (!(ot == 's') && !(ot == 'i')) {
			lgread = "none";
		}
		longread = "";
		longqual = "";
		longreadln = 0;
		return lgread;
		// makeLongRead method (END)
	}

	// /error if((lb=='N') || (lb=='N')) corected on 05/14/2015 statment :
	// if((lb!='N') && (rb=='N')) was refused
	public static Vector<Character> ChooseBaseCalculateQuality(char lb,
			char rb, char lq, char rq) {
		Vector<Character> resc = new Vector<Character>();
		int lqi = 0, rqi = 0;
		int bld = 0, bal = 0, randomNum = 0;
		if ((lb == 'N') || (rb == 'N'))
		// error if((lb=='N') || (lb=='N')) corected on 05/14/2015 statment :
		// if((lb!='N') && (rb=='N')) was refused
		{
			if ((lb == 'N') && (rb == 'N')) {

				resc.add('N');
				resc.add(qar[0]);
				bld = 1;
			}
			if ((lb == 'N') && (rb != 'N')) {

				resc.add(rb);
				resc.add(rq);
				bld = 1;
			}

			if ((lb != 'N') && (rb == 'N')) {
				resc.add(lb);
				resc.add(lq);
				bld = 1;
			}
			// the previous error was generating ocasionally cases where the bld
			// remained = 0 in [if (lb!='N') && (rb=='N')] statement since
			// (if((lb=='N') || (lb=='N'))
			// in this case if only rb=='N' couldn't be able to be captured and
			// in the following if statement it was going in the option of
			// choosing radomly the base that could involve the righrt
			// uncuptured N given lb!=N.
		}
		if (bld == 0) {
			lqi = 0;
			rqi = 0;
			bal = 0;
			for (int j = 0; j < qar.length; j++) {
				// calculate qualities first
				if (lq == qar[j]) {
					lqi = j;
					bal = bal + 1; // if both qualities have been calculated
									// then exit that loop
				}
				if (rq == qar[j]) {
					rqi = j;
					bal = bal + 1; // if both qualities have been calculated
									// then exit that loop
				}
				if (bal == 2) // if both qualities have been calculated then
								// exit that loop
					j = qar.length;
			}
			if (lb == rb) {
				resc.add(lb);
				resc.add(qar[(lqi + rqi)]);
			}
			if (lb != rb) {
				if (lqi == rqi) {
					Random rand = new Random();
					randomNum = rand.nextInt((2 - 1) + 1) + 1;
					// System.out.println(randomNum);
					if (randomNum == 1)
						resc.add(lb);
					if (randomNum == 2)
						resc.add(rb);
					resc.add(qar[0]);
				} else if (lqi != rqi) {
					if (lqi > rqi)
						resc.add(lb);
					if (lqi < rqi)
						resc.add(rb);
					resc.add(qar[Math.abs(lqi - rqi)]);
				}
			}// End of if(lb!=rb)
		}// End of if(bld==0)
		return resc;
	}

	public static boolean baseQuality(int gt, char ql) {
		boolean bl = false;
		double qi = 0;
		int qsum = 0;
		for (int j = 0; j < qar.length; j++) {
			if (ql == qar[j]) {
				qi = j;
				// System.out.println("Leter : " + seqAr[i] + "  Quality : " +
				// quaAr[i] + "  Score = " + qar[j] + " index = " + j +
				// "  Sum = " + qsum);
			}
		}
		// qi=qi-33;
		if (qi > gt)
			bl = true;
		return bl;
	}

	/*
	 * public static void RestoreBiCoef(String root) throws IOException,
	 * InterruptedException {
	 * 
	 * // Populate the standard binomial coeficient Long List of Long arrays //
	 * from a file Vector<String> bcv = new Vector<String>(); String binpath =
	 * root + "binocoefs.ssv"; bcv = readTheFileIncludeFirstLine(binpath,
	 * "binocoefs.ssv");
	 * 
	 * String str = ""; long lg = 0;
	 * 
	 * // add the first element to correspond to zero which w2il no use it. //
	 * The triangle will start from n = 1
	 * 
	 * Vector<Long> aux = new Vector<Long>(); aux.add(lg); Long[] var =
	 * aux.toArray(new Long[aux.size()]); dbtb.add(var);
	 * 
	 * // 0 //not used // 0 1 // 0 1 1 // 0 1 2 1 // 0 1 3 3 1 // 0 1 4 6 4 1 //
	 * 0 1 5 10 10 5 1 for (int i = 0; i < bcv.size(); i++) { lg = 0; aux = new
	 * Vector<Long>(); aux.add(lg); // add first element to vector (index 0 is
	 * not used)
	 * 
	 * str = bcv.get(i); String[] strwrds = str.split(" "); for (int y = 0; y <
	 * strwrds.length; y++) { lg = Long.parseLong(strwrds[y]); aux.add(lg); }
	 * Long[] nvkarr = aux.toArray(new Long[aux.size()]); dbtb.add(nvkarr); } }
	 * 
	 * public static void calculateBinCoef(int lrlength) { double val = 0.0;
	 * long vall = 0; long lg = 0; Vector<Long> temp = new Vector<Long>(); for
	 * (int i = 0; i <= lrlength; i++) {
	 * 
	 * // System.out.println(i); temp = new Vector<Long>(); temp.add(lg); // 0
	 * //not used // 0 1 // 0 1 1 // 0 1 2 1 // 0 1 3 3 1 // 0 1 4 6 4 1 // 0 1
	 * 5 10 10 5 1 for (int j = 0; j <= i; j++) { val = 1; for (int count = 0;
	 * count < j; count++) { val = val * (i - count) / (j - count); } vall =
	 * Math.round(val);
	 * 
	 * temp.add(vall); } Long[] nvkarr = temp.toArray(new Long[temp.size()]);
	 * dbtb.add(nvkarr); } }
	 */

	// Changed to use double 05/22/2015
	public static void calculateBinCoef(int lrlength) throws IOException,
			InterruptedException {
		Vector<Double> dbv = new Vector<Double>();
		for (int i = 0; i <= lrlength; i++) {
			dbv = new Vector<Double>();
			dbv.add(0.0);
			// 0 //not used
			// 0 1
			// 0 1 1
			// 0 1 2 1
			// 0 1 3 3 1
			// 0 1 4 6 4 1
			// 0 1 5 10 10 5 1
			for (int j = 0; j <= i; j++) {
				double val = 1;
				for (int count = 0; count < j; count++) {
					val *= (i - count) / (j - count);
				}
				dbv.add(val);
			}
			Double[] ndvkarr = dbv.toArray(new Double[dbv.size()]);
			dbtb.add(ndvkarr);
		}
	}

	public static double binomCDF(int n, int m) {
		double bval = 0;
		double sum = 0;
		for (int i = 0; i <= m; i++) {
			bval = dbtb.get(n)[i];
			sum += bval * Math.pow(0.75, i) * Math.pow(0.25, (n - i));
		}
		return sum;
	}
	// CREATE FRAGS(END)

	// CREATE TABLES(START)


	// public static void RetrieveDataFromBam( String bash_path,
	// String samtools_path,
	// String results3_path,
	// Vector<String> dtflnmspth,
	// Vector<String> dtflnms,
	// Vector<String> prms) throws IOException, InterruptedException
	// {

	public static Vector<String> RetrieveDataFromBam(String bash_path,
			String samtools_path, String temp_path, String locus,
			Vector<String> dtflnmspth) throws IOException, InterruptedException {
		temp_path = temp_path + "/";
		String samtoolsview = "";
		if ((samtools_path.equals("none")) || (samtools_path.equals("")))
			samtoolsview = " samtools view ";
		else
			samtoolsview = " " + samtools_path + " view ";
		//	samtoolsview = " " + samtools_path + "/samtools view ";
		System.out.println("samtoolsview : " + samtoolsview);
		Vector<String> currentsample = new Vector<String>();
		Vector<String> totaldatacol = new Vector<String>();
		Vector<String> bmnames = new Vector<String>();
		File fil;
		String bmdtnm = "";
		for (int o = 0; o < dtflnmspth.size(); o++) {
			fil = new File(dtflnmspth.get(o));
			bmdtnm = fil.getName();
			bmnames.add(bmdtnm);
		}

		ProcessBuilder pb;
		try {

			System.out.println(locus);
			for (int i = 0; i < dtflnmspth.size(); i++) // go through each bam
														// file
			{

				pb = new ProcessBuilder("/bin/bash", "-c", samtoolsview
						+ dtflnmspth.get(i) + " " + locus + " " + "-o "
						+ temp_path + i + ".txt");
				// pb = new ProcessBuilder("/bin/bash", "-c",
				// " /projects/bsi/bictools/bin/samtools view "
				// + dtflnmspth.get(i) + " " + locus + " "
				// + "-o " + results1_path + i + ".txt");

				Process p = pb.start();
				System.out.println("Process Started ....");
				p.waitFor();

				System.out.println("......Process Ended.");
				System.out.println("File ready.");

				currentsample = readTheFileIncludeFirstLine(temp_path + i
						+ ".txt", i + ".txt");
				System.out.println("Current Sample size = " + currentsample.size());
				File f = new File(temp_path + i + ".txt");
				f.delete();

				// System.out.println("Vector is done for mynewfile" + i +
				// ".txt");
				for (int k = 0; k < currentsample.size(); k++)
					totaldatacol.add(bmnames.get(i) + "\t"
							+ currentsample.get(k));
				// It looks like the new data line includes the origine (bam
				// file name) and the actual line (starting from the read name)

			}
		} finally {
			System.out.println("ProcessEnded");
		}
		// writeToFile(temp_path + "/testerror.tsv", totaldatacol);
		return totaldatacol;
	}

	/**
	 * Method to be documented
	 * 
	 * @param vlns
	 * @param refseq
	 * @return
	 */
	public static int gettenpercentorgtmimath(Vector<String> vlns, String refseq) {
		char[] refsech = refseq.toCharArray();
		int res = 0;
		String line = "", seq = "";
		int counter;
		int readactualsize = 0;
		for (int i = 0; i < vlns.size(); i++) {
			counter = 0;
			readactualsize = 0;
			line = vlns.get(i);
			String[] words = line.split("\t");
			seq = words[3];

			char[] sech = seq.toCharArray();
			for (int p = 0; p < sech.length; p++)
				if (sech[p] != 'P')
					readactualsize = readactualsize + 1;
			if (refsech.length != sech.length) {
				System.out.println("Primer length = " + refsech.length
						+ "  Read Length = " + sech.length);
				System.out.println("refseq     = " + refseq);
				System.out.println("Paded Read = " + seq);
				// endOfTheProgram();

			}
			for (int j = 0; j < Math.min((double) refsech.length,
					(double) sech.length); j++) {
				if (sech[j] != 'P')
					if (sech[j] == refsech[j])
						counter = counter + 1;
			}
			double resu = (double) counter / (double) readactualsize;
			// System.out.println(resu);
			if (resu < 0.9)
				res = res + 1;
		}

		return res;
	}

	// QC quality control implemented to give additional information
	// that will help evaluate and minimize procedural error
	// The methog gradualy increment counters hold as element of
	// the global integer array 'percmis'. For each read the percentage of
	// mismatches is calculated
	// and if it falls in the [1-10] area the coresponding counter is
	// incremented.
	// Before the program ends a file with these results is created.
	public static void QC(Vector<String> vlns, String refseq) {
		// percmis
		char[] refsech = refseq.toCharArray();
		int res = 0;
		String line = "", seq = "";
		int counter;
		int readactualsize = 0;
		for (int i = 0; i < vlns.size(); i++) {
			counter = 0;
			readactualsize = 0;
			line = vlns.get(i);
			String[] words = line.split("\t");
			seq = words[3];

			char[] sech = seq.toCharArray();

// 			QUALITY DEFINE No -> Ps Ds ds
// 			Calculate read actual size
			for (int p = 0; p < sech.length; p++)
			{
				if ( (sech[p] != 'P') && (sech[p] != 'D') && (sech[p] != 'd') && (sech[p] != 'N') )
					readactualsize = readactualsize + 1;
			}

// 			Print error message if pading didn't work and paded read size is not equal to reference size 	
			if (refsech.length != sech.length)
			{
				System.out.println("Primer length = " + refsech.length
						+ "  Read Length = " + sech.length);
				System.out.println("refseq     = " + refseq);
				System.out.println("Paded Read = " + seq);
				// endOfTheProgram();

			}
			
// 			From the rest read nucleotides calculate the number of matches
			for (int j = 0; j < Math.min((double) refsech.length,
					(double) sech.length); j++) 
			{
				if ( (sech[j] != 'P') && (sech[j] != 'D') && (sech[j] != 'd') && (sech[j] != 'N') )
					if (sech[j] == refsech[j])
						counter = counter + 1;
			}
			
			double resu = 1 - (double) counter / (double) readactualsize;
			// System.out.println(resu);
			if ((resu >= 0.0) && (resu <= 0.1))
				percmis[0] = percmis[0] + 1;
			if ((resu > 0.01) && (resu <= 0.02))
				percmis[1] = percmis[1] + 1;
			if ((resu > 0.02) && (resu <= 0.03))
				percmis[2] = percmis[2] + 1;
			if ((resu > 0.03) && (resu <= 0.04))
				percmis[3] = percmis[3] + 1;
			if ((resu > 0.04) && (resu <= 0.05))
				percmis[4] = percmis[4] + 1;
			if ((resu > 0.05) && (resu <= 0.06))
				percmis[5] = percmis[5] + 1;
			if ((resu > 0.06) && (resu <= 0.07))
				percmis[6] = percmis[6] + 1;
			if ((resu > 0.07) && (resu <= 0.08))
				percmis[7] = percmis[7] + 1;
			if ((resu > 0.08) && (resu <= 0.09))
				percmis[8] = percmis[8] + 1;
			if ((resu > 0.09) && (resu <= 0.10))
				percmis[9] = percmis[9] + 1;
			if ((resu > 0.10))
				percmis[10] = percmis[10] + 1;
		}
	}

	/**
	 * 'maketables' method creates tables from bam files that contain reads
	 * generated from amplicon sequencing or capture. The method accepts sites
	 * -information contained in elements of the 'prms' String Vector (passed 
	 * as a parameter). Each site information is organized into words separated by
	 * tabs in a sentence(String Vector element). The method uses an outer loop 
	 * that goes through each site. At the beginning of the loop the site 
	 * (Vector String element) is splited into words. Next the chromosomal 
	 * coordinates for the locus are constructed and used to query all bam files.
	 * The data lines returned from all bam files are filtered out using 
	 * mapping quality >= 30. The remaining data will undergo a sequence(read) and
	 * quality transformation. In sort by looking at the read position, and the cigar
	 * information the sequence and its quality will be transformed in the 
	 * following way : The read will start and end in the site range start and end 
	 * points(reads that exced those points will be trimmed out; reads that are within
	 * those points Ps will be added to the side(s). Eventully all reads after 
	 * transformation will have the same start and end points. Before the transformation
	 * process the reads are transformed acording to cigar( D(s) are added to the 
	 * sequence; insertion(s) are removed from sequence; Ss and Hs are handled properly.
	 * If there is N(s) information in cigar then lower case d(s) are added to the the
	 * corresponding read positions in both read(sequence) and quality fields.
	 * After this primary preparation all Read Nucleotides can be mapped directly to
	 * the reference. In sort there is a primary treansformation that regains the length 
	 * of read(expected); information for deletions and insertions also are retrieved 
	 * along with the read. In the secondary transformation there is padding or trimming 
	 * of the reads in such a way that all reads start and end at the same start and end
	 * points as indicated from site information. This was the alignment process. Now the 
	 * alligned reads can be put in a imaginary table where each column represents a 
	 * reference position and each row represents the (paded or trimited) read sequence 
	 * itself. The sites information contain also the reference sequence. This means 
	 * that for each position there is only one Nucleotide that is expected. By looking
	 * at the column(particular position) someone can count all expected nucleotides 
	 * found. In similar way all other three Nucleotides have their own count number for 
	 * a particular position. This was a simple way to describe the function of this 
	 * method. The method populates a two dimensional String array. This can be imagined as
	 * a table where the rows represent a particular position, and the columns represent 
	 * information about that position. Column headers are Chromosome, Position, Counts Total 
	 * Nucleotide found, Expected Nucleotide, Alternative Nucleotide, Counts of As, 
	 * Counts of Cs, Counts of Ts, Counts of Gs. Other positional information has been
	 * added to the array in order to provide information about next position deletion(s)
	 * (one(only one deletion in the next position found) or more deletions(consequtive 
	 * deletions found starting from next position). About insertion(s) the next position 
	 * insertion(s) is recorded. The number of next position insertions and the sequence
	 * and its counts is also recorded. More information about this method will be added
	 * later.
	 * INDELS Explaining
	 * The indels are reported as followng:
	 * counting of counting of number of deletions + informative of 
	 * for example A,A,A,A means that for the next position (A) have been found
     * four valid reads that have the following position deleted The deletions
	 * are all single nucleotide deletions thas the del position one will
     * be incremented by four
     * similarly ACCT means that starting from next position there is a deletion
     * found of the folowing four nucleotides. The del position four will be
     * icreased by one(ACCT are the sequence of nucleotides after this position
	 * The insertion follow the similar rules.
	 * Insertions and deletions informative fields are null if there 
     * is zero coverage 
	 * It calls samtools command
	 * to query bam files from command line.
	 * Parameters Explain   
	 * @param bash_path : path to bash  
	 * @param samtools_path : path to samtools
	 * @param input : 
	 * @param output : 
	 * @param perc : 
 	 * @param prms : String Vector with sites information
	 * @param bms : comma separated String with bampaths
	 * @param tmppath : path to temporary folder
	 * @param dcuti : int varible that is used to filter out all reads that 
	 * start and end in a distance more than dcuti ( abs(read start - site range start)>dcuti  ) 
	 * @throws IOException
	 * @throws InterruptedException
	 */
	
// 	add example of returned String Vector element in each step
// 	makeTables(bsh_pth, smtls_pth, input, output, perc, prmsv, bms, tmppath);

	public static void makeTables(String bash_path, String samtools_path,
			String input, String output, String perc, Vector<String> prms,
			Vector<String> bms, String tmppath, int dcuti, String pad) throws IOException,
			InterruptedException {

// **************************************************************************************************************************\\
		int padi=0;
		if(isPosUnsignInteger(pad))
			padi=Integer.parseInt(pad);
		String rdnm = "", rdnmln = "";
		String outFilename = output + ".insdels.gz";
		GZIPOutputStream gzipout = new GZIPOutputStream(new FileOutputStream(
				outFilename));
		BufferedWriter bw1 = new BufferedWriter(new OutputStreamWriter(gzipout));
// **************************************************************************************************************************\\

		Vector<String> pmis = new Vector<String>();
		for (int i = 0; i < percmis.length; i++)
			percmis[i] = 0;
//		1. Start
// 		2. Go throug each site and retrieve the corresponding data from all bam files
// 		3. Populate the pareticular part of the table array.
// 		4. Continue until no sites have been left
// 		5. Write to table file
//		6. End
		System.out.println("Distance = " + dcuti);

		String pline = "";
		int psz = 0, tblsz = 0; // table size = maximum pair primers sequence
		// prim seq = left primer start - sequence - right primer start(end)
		int tpogtmm = 0;
		/*
		 * if(prms.size()==1) { pline=prms.get(0); String[] pwrds =
		 * pline.split("\t"); tblsz = pwrds[3].length(); }
		 */
		// for(int i=0; i<prms.size()-1; i++)
		for (int i = 0; i < prms.size(); i++) {
			pline = prms.get(i);
			String[] pwrds = pline.split("\t");
			// System.out.println(pline);
			psz = pwrds[3].length();
			// System.out.println(psz);
			// calculate the maximum to update the results array
			if (tblsz < psz)
				tblsz = psz;
		}
		noisetable = new String[tblsz][35];
		int mqi = 30;
		System.out.println("Maping Quality > " + mqi);
		Vector<String> reslts = new Vector<String>();
		String cursmline = "", fnline = "", bflnm = "", chromosome = "", positions = "", MAPQ = "", cigar = "", sequence = "", flg = "", qual = "", truefalse = "", strand = "", readname = "";
		String primline = "";
		int b = 0;
		int position = 0;
		int refchrn = 0, refstart = 0, refend = 0, refsize = 0, refactualsize = 0;
		String refseq = "", alter = "", altern="", smplnmi = "", smplnmi1 = "";
		
		String dataline = "", padedseq = "", datseq = "", datqual = "", newline = "";
		int refid = 0, actualend = 0, posii = 0, posiii=0;
		String dseq = ""; // holds deletion sequences
		String vcline = "", vcseqpaded = "";
		int vcposition = 0, vcrefid = 0, vcend = 0;
// 05/03/2015 Changed to accept Chromosome X and Y	
		String vcrefidS="";
// 05/03/2015 Changed to accept Chromosome X and Y	
		String vcseq = "", vcqual = "";
		int bqlcount = 0;
		char expected_base = 'N';
		String vDeletions = "", vInsertions = "";
		int d1 = 0, d2 = 0;
		int x = 0, y = 0;

		String resline = ""; // results line
// 		Prepare table Header String 
		resline = "chromosome" + "\t" + "site" + "\t" + "x" + "\t"
				+ "y"
				+ "\t"
				+ "base_expected"
				+ "\t"
				+ "alternative"
				+ // 5th column has added(alternative
				"\t" + "number_base_expected" + "\t" + "number_of_As" + "\t"
				+ "number_of_Cs" + "\t" + "number_of_Ts" + "\t"
				+ "number_of_Gs" + "\t" + "number_of_Ds" + "\t"
				+ "number_of_Ps" + "\tD1\tD2\tD3\tD4\tD5\tD6\tD7\tD8\tD9\tD10"
				+ "\tI1\tI2\tI3\tI4\tI5\tI6\tI7\tI8\tI9\tI10"
				+ "\tDelInfo\tInsInfo";

		// + "\t" + "N";
		// System.out.println(resline);
		reslts.add(resline);
		resline = "";

		// public static void CreateNoiseTables(String tables_path, String
		// results2_path, Vector<String> prim, int qsc) throws IOException,
		// InterruptedException
		// {
		int qsc = Integer.parseInt(perc);

		Vector<String> currentsample = new Vector<String>();
		Vector<String> cursmlq = new Vector<String>();
		// Vector<String> results = new Vector<String>();
		Vector<String> vlns = new Vector<String>();

		Vector<String> insposnumseq = new Vector<String>();
		int curposition = 0;
		Vector<Integer> delposv = new Vector<Integer>();
		Vector<Integer> delnumv = new Vector<Integer>();
		String dsinf = "";
		Vector<Integer> insposv = new Vector<Integer>();
		Vector<Integer> insnumv = new Vector<Integer>();
		Vector<String> insseqv = new Vector<String>();
		String isinf = "";

		int dlbfpos = 0, dlnum = 0, inbfpos = 0, innum = 0;
		// position before deletion insertion and number of insertion or deletion

		String inseq = "";

		String headers = "";
		Vector<String> totaldatacol = new Vector<String>();
		String prmline = "", CHROM = "", Left_start = "", Right_start = "", seq = "", site = "", REF = "", ALT = "", POSs="", REFs="", ALTs="";
		String Left_primer = "", Right_primer = "", locus = "";
		String strs = "", ends = "";
		Integer Left_starti = 0, Right_starti = 0, MAPQi = 0;
		// Vector<String> PrimCumula =new Vector<String>();

		// Create a folder in the temporary folder with a unique name to hold
		// temporary data files
		String tmpuiid = "";
		tmpuiid = tmppath + "/"
				+ UUID.randomUUID().toString().replaceAll("-", "i");
		File dutmp = new File(tmpuiid);
		if (!dutmp.exists())
			dutmp.mkdir();
		// Create a folder in the temporary folder with a unique name to hold
		// temporary data files
		
//		06/08/2015 added in the case more sites are included when investigating one range	
		Vector<Integer> posv = new Vector<Integer>();	
		Vector<String> refv = new Vector<String>();
		Vector<String> altv = new Vector<String>();

//		all three fields will have same 
		int ercnt=0;
		Vector<String> testError = new Vector<String>();

// 05/03/2015 Changed to accept Chromosome X and Y		
		String refidS="";
// 05/03/2015 Changed to accept Chromosome X and Y
		for (int h = 0; h < prms.size(); h++)
		//for(int h=0; h<4; h++)
		// for(int h=64; h<66; h++)
		{
			cursmlq = new Vector<String>();
			vlns = new Vector<String>();
			prmline = prms.get(h);
			String[] prwords = prmline.split("\t");
			System.out.println("-------------------> " + prmline);
			if (prwords.length > 1) // *************************************************************---------------------------------------------
			{
				// one part of the comparison <----
				// cr startc endc ResSeq POS REF ALT
				//  0   1     2     3     4   5   6
				CHROM = prwords[0];
				Left_start = prwords[1];
				Right_start = prwords[2];
				// site = prwords[4];
				POSs = prwords[4];
				REFs = prwords[5];
				ALTs = prwords[6];

// 05/03/2015 Changed to accept Chromosome X and Y
				// refchrn = Integer.parseInt(CHROM); // CHROM=prwords[0];
				refstart = Integer.parseInt(Left_start); // start of left primer
				refend = Integer.parseInt(Right_start); // start of right primer
				posv = new Vector<Integer>();	
				refv = new Vector<String>();
				altv = new Vector<String>();

				String[] pwds = POSs.split(",");
				site = pwds[0];
				String[] rwds = REFs.split(",");
				String[] awds = ALTs.split(",");
				for(int mls=0; mls<pwds.length; mls++)			// multiple sites
				{
				   posv.add(Integer.parseInt(pwds[mls]));
				   refv.add(rwds[mls]);
				   altv.add(awds[mls]);
				}

//	**************************************************************************				
//				posii = Integer.parseInt(site);
//				alter = prwords[6];
//	**************************************************************************
// Take the first elements of the vectors in order to keep the program running

				REF=refv.get(0);
				ALT=altv.get(0);
				posii=posv.get(0);

				refseq = prwords[3]; // The reference sequence itself
				char[] carrs = refseq.toCharArray();
				refsize = refseq.length();

				System.out.println("Chr : chr" + CHROM + "\nrefstart: "
						+ refstart + "\nrefend : " + refend + "\nrefsize : "
						+ refsize);
				refactualsize = refend - refstart;

				Left_starti = Integer.parseInt(Left_start);
				Right_starti = Integer.parseInt(Right_start);
// 	06/03/2015 Change made in the extension of the border when construct locus(coordinates to query bam files
//	06/03/2015 In new data the sites could be close enough. So the extension have been reduced to 10 bases in each side(it was 150)
				locus = CHROM + ":" + (Left_starti - 10) + "-"
						+ (Right_starti + 10);
				System.out.println(locus);
				rdnmln = refseq;
				bw1.write(rdnmln);
				bw1.newLine();

				// System.out.println("bam[0] = " + bms.get(0));

				totaldatacol = RetrieveDataFromBam(bash_path, samtools_path,
						tmpuiid, locus, bms);
				System.out.println("Elaborate Data for primer : " + h);

				//if (h==2)
				//writeToFile("/data5/experpath/vasm/vasm/NextGen/niko/pancreatic_capture110615/samples/1329805207/tables/data_Error_Test.tsv", totaldatacol);
				System.out.println("Current Vector Size : " + totaldatacol.size());
				testError = new Vector<String>();
				for (int i = 0; i < totaldatacol.size(); i++)
				{
					cursmline = totaldatacol.get(i);
					String[] cslwords = cursmline.split("\t");
					// System.out.println(i + "<-    ->" + cslwords.length);
					if (cslwords.length > 9)
					{
						bflnm = cslwords[0];
						readname = cslwords[1];

						// inserted on 03/16/2015 to be
						// farther elaborated from
						// PrepareTheDataAndInDls method

						strand = cslwords[2];
						chromosome = cslwords[3];
						// System.out.println(chromosome);
						positions = cslwords[4];
						MAPQ = cslwords[5];
						MAPQi = Integer.parseInt(MAPQ);
						cigar = cslwords[6];
						flg = cslwords[10];
						sequence = cslwords[10];
						qual = cslwords[11];
						// System.out.println("MAPQi : " + MAPQi +
						// "  Read Length = " + sequence.length() +
						// "  Qual Length = " + qual.length());
						// truefalse=cslwords[12];
						// fnline="", chromosome="", position="", Ncigar=,
						// cigar="", sequence="", flg="", qual="", truefalse="";
						if (MAPQi > mqi)
						// if(MAPQi>0)
						{
							fnline = bflnm + "\t" + site + "\t" + CHROM + "\t"
									+ REF + "\t" + ALT + "\t" + chromosome
									+ "\t" + positions + "\t" + MAPQ + "\t"
									+ cigar + "\t" + sequence + "\t" + qual
									+ "\t" + strand + "\t" + readname;
							cursmlq.add(fnline);
							//if(h==2)
								//testError.add(sequence);
						}
					}
					// if(cslwords.length<9)
					// System.out.println(cslwords[3]);
				}

				// if(h==2)
//					writeToFile("/data5/experpath/vasm/vasm/NextGen/niko/pancreatic_capture110615/samples/1329805207/tables/data_Error_seq_Test.tsv", testError);
//				System.out.println("Remained Vector Size : " + cursmlq.size());
				
				int totallines = cursmlq.size();
				// // for (int k=0; k<totallines; k++)
				// //cursmlq.set(k, cursmlq.get(k)+ "\t" + totallines);
				// System.out.println("JJJJJ----> Size = " + cursmlq.size());
				cursmlq = PrepareTheDataAndInDls(cursmlq);
//				System.out.println("Size after Prepare The data: " + cursmlq.size());
				// 03/16/2015 cursmlq contains prepared reads with additional
				// information
				// each cursmlq's line contains read name at index 10 starting
				// from 0
				// for common reads(no indels) read name = cmn
				// for reads that have insertions or/and deletions the read name
				// is the actual name
				// for (int k=0; k<cursmlq.size(); k++)
				// PrimCumula.add(cursmlq.get(k));
				// writeToFile(tmppath + "/_" + h + "_" + h + "_test007.tsv",
				// cursmlq);
				// writeToFile(tmppath + "/_" + h + "_" + h +
				// "_test00hfgh7.tsv", PrimCumula);
				// PrimCumula = new Vector<String>();
				for (int i2 = 0; i2 < cursmlq.size(); i2++) // go through each
															// read(line)
				{
					rdnm = "";
					dataline = cursmlq.get(i2);
					String[] dawords = dataline.split("\t");
//					refid = Integer.parseInt(dawords[0]);
					refidS = dawords[0];
					position = Integer.parseInt(dawords[1]);
					actualend = Integer.parseInt(dawords[2]);
					datseq = dawords[3];
					// System.out.println(" datseq : " + datseq);
					datqual = dawords[4];
					// System.out.println(" datqual : " + datqual);
					vDeletions = "";
					vInsertions = "";
					cigar = "";
					vDeletions = dawords[5];
					vInsertions = dawords[6];
					cigar = dawords[7];
					rdnm = dawords[10]; // read name "cnm" if the read is
										// normal.
					smplnmi = dawords[11];
//  Limits
					if (//(position >= refstart - 100)
							//&& (position <= refstart + 100)
							//&& 
							(position + datseq.length() > refstart)
							&& (position <= refend)) // position matches
					{
						// ercnt=ercnt+1;
						padedseq = padDataSequence(datseq, datqual, position,
								refstart, refsize, padi);
						// results.add(padedseq); //
						// LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
						// after modifications are done for quality control the
						// padedseq
						// contains the paded sequence and its quality separated
						// by tab.

						// Distance / No Distance Control
						if (dcuti != -1) {
							d1 = Math.abs(refstart - position);
							d2 = Math.abs(actualend - refend);
							// System.out.println("d1 = " + d1 + "  d2 = " +
							// d2);
							if ((padedseq.length() > 2)
									&& ((d1 < dcuti) && (d2 < dcuti))) {
							//	newline = refid + "\t" + refstart + "\t"
								newline = refidS + "\t" + refstart + "\t"
										+ actualend + "\t" + padedseq + "\t"
										+ vDeletions + "\t" + vInsertions
										+ "\t" + cigar + "\t" + rdnm + "\t"
										+ smplnmi;
								vlns.add(newline);
							}
						} else if (dcuti == -1) {
							if (padedseq.length() > 2) {
								newline = refidS + "\t" + refstart + "\t"
										+ actualend + "\t" + padedseq + "\t"
										+ vDeletions + "\t" + vInsertions
										+ "\t" + cigar + "\t" + rdnm + "\t"
										+ smplnmi;
								vlns.add(newline);
								if(h==2)
								testError.add(padedseq);
							}
						}

					}

					// Distance / No Distance Contro

				}
//				System.out.println("Size after pading the data: " + vlns.size() + "  counter : " + ercnt);
//				Retrieve paded data fo visualization
//				if(h==2)
//					writeToFile("/data5/experpath/vasm/vasm/NextGen/niko/pancreatic_capture110615/samples/1329805207/tables/padeddata_Error_seq_Test.tsv", testError);
// 				System.out.println("Size after pading : " + cursmlq.size())
// 				QUALITY CONTROL Calculate the percenta ge of mismatches for				
// 				each read, group those reads and then count them
// 				creatte a report log file
				QC(vlns, refseq);
// 				QUALITY CONTROL Calculate the percentge of mismatches

				// tpogtmm = tpogtmm + gettenpercentorgtmimath(vlns,refseq);
				// writeToFile(tmppath + "/_" + h + "_" + h +
				// "HereError_vlns.tsv", vlns);

				// ------------- Add component to the table

				if (vlns.size() > 0) 
				{
					System.out.println("lines size = " + vlns.size());
					// prepare the array numbers that hold string numbers with
					// zeros
					for (int nt = 0; nt < noisetable.length; nt++)
						// for(int nt1=5; nt1<13; nt1++)
						for (int nt1 = 0; nt1 < 33; nt1++)
							noisetable[nt][nt1] = "0";
					for (int nt = 0; nt < noisetable.length; nt++)
						// for(int nt1=5; nt1<13; nt1++)
						for (int nt1 = 33; nt1 < 35; nt1++)
							noisetable[nt][nt1] = "";

					// Go trough each data line take the read and create a
					// character arrey that
					// holds its bases as elements
					// The collective array has designed in such a way that can
					// update all its elements
					// with each data line passesd. One can imagine it as a
					// table that has # (noisetablelength) rows
					// and 14 columns. Each row contains information about one
					// referense position
					// and how this position has found in all particular
					// primer-pair resulting reads.
					// The update procedure is recorded line by line

					for (int i3 = 0; i3 < vlns.size(); i3++) 
					{
						vcline = vlns.get(i3);
						// System.out.println(" here - > " + vcline);
						String[] vcwords = vcline.split("\t");
				// 05/03/2015 Changed to accept Chromosome X and Y
				//		vcrefid = Integer.parseInt(vcwords[0]);
						vcrefidS = vcwords[0];
						// 05/03/2015 Changed to accept Chromosome X and Y
						vcposition = Integer.parseInt(vcwords[1]);
						vcend = Integer.parseInt(vcwords[2]);
						vcseq = vcwords[3]; // sequence was returned paded in
											// combination with the quality
											// paded(tab separated one String)
						vcqual = vcwords[4];
						vDeletions = vcwords[5];
						vInsertions = vcwords[6];
						cigar = vcwords[7];
						rdnm = vcwords[8];
						smplnmi1 = vcwords[9];
						// parallel integer vectors for deletion position and
						// size
						delposv = new Vector<Integer>();
						delnumv = new Vector<Integer>();
						// parallel integer vectors for insertion position,
						// size, sequence
						insposv = new Vector<Integer>();
						insnumv = new Vector<Integer>();
						insseqv = new Vector<String>();
						insposnumseq = new Vector<String>();

						// ****************************************************************************************************************************

						if (!(rdnm.equals("cmn"))) {
							rdnmln = rdnm + "\t" + smplnmi1 + "\t" + vcrefidS
									+ "\t" + position + "\t" + cigar + "\t"
									+ vcseq + "\t" + vDeletions + "\t"
									+ vInsertions;
							// System.out.println(" rdnmln - > " + rdnmln);
							bw1.write(rdnmln);
							bw1.newLine();
						}

						if (!(vDeletions.equals("none"))) {
							String[] vdlswrds = vDeletions.split(",");
							for (int ds = 0; ds < vdlswrds.length; ds++) {
								String[] singdl = vdlswrds[ds].split("_");
								dlbfpos = Integer.parseInt(singdl[0]);
								delposv.add(dlbfpos);
								dlnum = Integer.parseInt(singdl[1]);
								delnumv.add(dlnum);
							}
						}

						if (!(vInsertions.equals("none"))) {
							String[] vinswrds = vInsertions.split(",");
							for (int is = 0; is < vinswrds.length; is++) {
								String[] singins = vinswrds[is].split("_"); // single
																			// insertion
								inbfpos = Integer.parseInt(singins[0]); // position
																		// before
																		// insertion(s)
								insposv.add(inbfpos);
								innum = Integer.parseInt(singins[1]);
								insnumv.add(innum);
								inseq = singins[2];
								insseqv.add(inseq);
							}
						}
						// System.out.println(vcseq);
						// System.out.println(vcqual);

						// vcseqpaded = padDataSequence(vcseq, vcposition,
						// refstart);

						char[] carsq = vcseq.toCharArray(); // second part of
															// the comparison
															// <-----
						char[] carql = vcqual.toCharArray();
						// for(int s=0; s<carql.length; s++)
						// System.out.println(carsq[s] + " " + carql[s]);

						// first index 0 0 will hold the position
						// example noisetable[0][0] : chromosome
						// example noisetable[0][1] : site
						// example noisetable[0][2] : x
						// example noisetable[0][3] : y
						// example noisetable[0][4] : base expected
						// example noisetable[0][5] : alternative
						// examole noisetable[0][6] : #base expected
						// example noisetable[0][7] : #A
						// example noisetable[0][8] : #C
						// example noisetable[0][9] : #T
						// example noisetable[0][10] : #G
						// example noisetable[0][11] : #D
						// example noisetable[0][12] : #P
						// example noisetable[0][13] : D1
						// example noisetable[0][14] : D2
						// example noisetable[0][15] : D3
						// example noisetable[0][16] : D4
						// example noisetable[0][17] : D5
						// example noisetable[0][18] : D6
						// example noisetable[0][19] : D7
						// example noisetable[0][20] : D8
						// example noisetable[0][21] : D9
						// example noisetable[0][22] : D10
						// example noisetable[0][23] : I1
						// example noisetable[0][24] : I2
						// example noisetable[0][25] : I3
						// example noisetable[0][26] : I4
						// example noisetable[0][27] : I5
						// example noisetable[0][28] : I6
						// example noisetable[0][29] : I7
						// example noisetable[0][30] : I8
						// example noisetable[0][31] : I9
						// example noisetable[0][32] : I10
						// example noisetable[0][33] : DelsInfo
						// example noisetable[0][34] : InssInfo

						x = 1;
						y = 1;
						// Imagine that for each read the following loop will
						// update all
						// the same elements (in (#) = (noise-table length)
						// lines and 14 columns) the number of data lines.
						// example: Pass through all (#) read bases (that
						// correspond to position references)
						// and update. Imagine: create a big sheet with all
						// particular data alligned reads
						// the one under the other and look at each column((#)
						// total) verticaly to compare them
						// ACTCTCCTGGTA
						// ACTCTCCTGGTA
						// ACTCTCCTGGTA
						// ACTCTCCTGGTA
						// CCTCTCCTGGTA //This has a substitution in first
						// position
						// ACTCTCCTGGTA // align them and look verticaly for
						// changes
						// ACTCTCCTGGTA
						// ACTCTCCTGGTA

						// System.out.println(carsq.length);
						// System.out.println(carrs.length);
						// for (int i4=0; i4<#tlength; i4++)
						// System.out.println("vcseq  = " + vcseq);
						// System.out.println("vcqual = " + vcqual);

						// Create minimum upper bounter for the next loop
						int minub = Math.min(carsq.length, carrs.length);
						// System.out.println("carsq.length = " + carsq.length);
						// System.out.println("prim.length = " + carrs.length);
						// System.out.println("min = " + minub);

						// for (int i4=0; i4<carsq.length; i4++)
						for (int i4 = 0; i4 < minub; i4++) 
						{
							x = i4 + 1;
							y = carsq.length - x + 1;
							// if (x>300)
							// x=0;
							// if(y>300)
							// y=0;
							// System.out.println(" i4 : " + i4);
							expected_base = carrs[i4];
							curposition = vcposition + i4;
//	Deletions	Deletions	Deletions	Deletions	Deletions	Deletions	Deletions
//	Deletions	Deletions	Deletions	Deletions	Deletions	Deletions	Deletions
							for (int dcnt = 0; dcnt < delposv.size(); dcnt++) 
							{
								if (curposition == delposv.get(dcnt))
								{
									if (delnumv.get(dcnt) <= 10)
										noisetable[i4][12 + delnumv.get(dcnt)] = Integer
												.toString(Integer
														.parseInt(noisetable[i4][12 + delnumv
																.get(dcnt)]) + 1);
									else
										noisetable[i4][12 + 10] = Integer
												.toString(Integer
														.parseInt(noisetable[i4][12 + 10]) + 1);
									// dsinf="D_"+delposv.get(dcnt) + "_" +
									// delnumv.get(dcnt) + ",";
									// System.out.println("----> " +
									// noisetable[i4][13] + " " +
									// noisetable[i4][12+delnumv.get(dcnt)] +
									// " " + dsinf);
									// if(delnumv.get(dcnt)<=10)
									// noisetable[i4][33] = noisetable[i4][33] +
									// "_" +
									// noisetable[i4][12+delnumv.get(dcnt)] +
									// dsinf;
									// else
									// noisetable[i4][33] = noisetable[i4][33] +
									// "_" + noisetable[i4][12+10] + dsinf;

									if (delnumv.get(dcnt) <= 10) {
										for (int qe = 0; qe < delnumv.get(dcnt); qe++) {
											if (i4 + qe + 1 < carrs.length)
												dseq = dseq
														+ carrs[i4 + qe + 1];
										}
										if (dseq.length() > 0)
											noisetable[i4][33] = noisetable[i4][33]
													+ dseq + ",";
										dseq = "";
									}
								}
							}
//	Deletions	Deletions	Deletions	Deletions	Deletions	Deletions	Deletions
//	Deletions	Deletions	Deletions	Deletions	Deletions	Deletions	Deletions


//	Insertions	Insertions	Insertions	Insertions	Insertions	Insertions	Insertions
//	Insertions	Insertions	Insertions	Insertions	Insertions	Insertions	Insertions
							// Ins start with 22
							for (int icnt = 0; icnt < insposv.size(); icnt++) {
								if (curposition == insposv.get(icnt)) {
									// System.out.println(cigar);
									// System.out.println("----> " +
									// noisetable[i4][22+insnumv.get(icnt)] +
									// " index  " + insnumv.get(icnt));
									// The total columns for recording
									// insertions are 10
									// Insertions> 10 were creating errors
									// If insertions > 10 the record insertion =
									// 10 (02/02/2015)
									if (insnumv.get(icnt) <= 10)
										noisetable[i4][22 + insnumv.get(icnt)] = Integer
												.toString(Integer
														.parseInt(noisetable[i4][22 + insnumv
																.get(icnt)]) + 1);
									else
										noisetable[i4][22 + 10] = Integer
												.toString(Integer
														.parseInt(noisetable[i4][22 + 10]) + 1);
									// isinf="I_"+insposv.get(icnt) + "_" +
									// insnumv.get(icnt) + "_" +
									// insseqv.get(icnt) +",";
									isinf = insseqv.get(icnt) + ",";

									// System.out.println("----> " +
									// noisetable[i4][23] + " " +
									// noisetable[i4][22+insnumv.get(icnt)] +
									// " " + isinf);
									if (insnumv.get(icnt) <= 10) {
										// noisetable[i4][34] =
										// noisetable[i4][34] + "_" +
										// noisetable[i4][22+insnumv.get(icnt)]
										// + isinf;
										noisetable[i4][34] = noisetable[i4][34]
												+ isinf;
									} else {
										noisetable[i4][34] = noisetable[i4][34]
												+ isinf;
										// noisetable[i4][34] =
										// noisetable[i4][34] + "_" +
										// noisetable[i4][22+10] + isinf;
									}
								}
							}
//	Insertions	Insertions	Insertions	Insertions	Insertions	Insertions	Insertions
//	Insertions	Insertions	Insertions	Insertions	Insertions	Insertions	Insertions

							// Ins start with 22
							//noisetable[i4][0] = Integer.toString(vcrefid); // chromosome

							noisetable[i4][0] = vcrefidS;
							noisetable[i4][1] = Integer.toString(vcposition
									+ i4); // site
							noisetable[i4][2] = Integer.toString(x); // x
							noisetable[i4][3] = Integer.toString(y); // y
							noisetable[i4][4] = Character
									.toString(expected_base); // expected
							noisetable[i4][5] = Character
									.toString(expected_base); // alternative if
																// position for
							for(int mls=0; mls<posv.size(); mls++)			// multiple sites
							{
				   				posiii = posv.get(mls);
								//System.out.println(posiii + " <---->  " + curposition); 
				   				altern = altv.get(mls);
								if (posiii == curposition) 
								{
									noisetable[i4][5] = altern;
								//	System.out.println("In");
								//	System.out.println(noisetable[i4][5] + "  --  " + noisetable[i4][4] );
								}
							}

							// if (expected_base=='N')
							// noisetable[i4][13]="1";
							// if (!(expected_base=='N'))
							// noisetable[i4][13]="0";
							// if (carsq[i4]==expected_base)
							// {
							// noisetable[i4][5]=Integer.toString(Integer.parseInt(noisetable[i4][5])+1);
							// }
							// for(int s=0; s<carql.length; s++)
							// System.out.println(carql[s]);
							// taken out of the next loop in order to keep track
							// of deletions without modifying the expected
							// number of bases.
							if (carsq[i4] == 'D') {
								noisetable[i4][11] = Integer.toString(Integer
										.parseInt(noisetable[i4][11]) + 1);
							}
							if (baseQuality(qsc, carql[i4])
									&& (!(carsq[i4] == 'P'))
									&& (!(carsq[i4] == 'D'))
									&& (!(carsq[i4] == 'N'))
									&& (!(carsq[i4] == 'd'))) {
								bqlcount = bqlcount + 1;
								// noisetable[i4][6]=Integer.toString(vlns.size());
								noisetable[i4][6] = Integer.toString(Integer
										.parseInt(noisetable[i4][6]) + 1);
								if (carsq[i4] == 'A') {
									noisetable[i4][7] = Integer
											.toString(Integer
													.parseInt(noisetable[i4][7]) + 1);
								}

								if (carsq[i4] == 'C') {
									noisetable[i4][8] = Integer
											.toString(Integer
													.parseInt(noisetable[i4][8]) + 1);
								}

								if (carsq[i4] == 'T') {
									noisetable[i4][9] = Integer
											.toString(Integer
													.parseInt(noisetable[i4][9]) + 1);
								}

								if (carsq[i4] == 'G') {
									noisetable[i4][10] = Integer
											.toString(Integer
													.parseInt(noisetable[i4][10]) + 1);
								}
								if (carsq[i4] == 'D') {
									noisetable[i4][11] = Integer
											.toString(Integer
													.parseInt(noisetable[i4][11]) + 1);
								}
								if (carsq[i4] == 'P') {
									noisetable[i4][12] = Integer
											.toString(Integer
													.parseInt(noisetable[i4][12]) + 1);
								}
							}
						}
					}
					// vcposition = Integer.parseInt(vcwords[5]);
					// vcrefid = Integer.parseInt(vcwords[4]);
					// vcseq =
					for (int i5 = 0; i5 < carrs.length; i5++) 
					{
						resline = noisetable[i5][0] + "\t" + noisetable[i5][1]
								+ "\t" + noisetable[i5][2] + "\t"
								+ noisetable[i5][3] + "\t" + noisetable[i5][4]
								+ "\t" + noisetable[i5][5] + "\t"
								+ noisetable[i5][6] + "\t" + noisetable[i5][7]
								+ "\t" + noisetable[i5][8] + "\t"
								+ noisetable[i5][9] + "\t" + noisetable[i5][10]
								+ "\t" + noisetable[i5][11] + "\t"
								+ noisetable[i5][12] + "\t"
								+ noisetable[i5][13] + "\t"
								+ noisetable[i5][14] + "\t"
								+ noisetable[i5][15] + "\t"
								+ noisetable[i5][16] + "\t"
								+ noisetable[i5][17] + "\t"
								+ noisetable[i5][18] + "\t"
								+ noisetable[i5][19] + "\t"
								+ noisetable[i5][20] + "\t"
								+ noisetable[i5][21] + "\t"
								+ noisetable[i5][22] + "\t"
								+ noisetable[i5][23] + "\t"
								+ noisetable[i5][24] + "\t"
								+ noisetable[i5][25] + "\t"
								+ noisetable[i5][26] + "\t"
								+ noisetable[i5][27] + "\t"
								+ noisetable[i5][28] + "\t"
								+ noisetable[i5][29] + "\t"
								+ noisetable[i5][30] + "\t"
								+ noisetable[i5][31] + "\t"
								+ noisetable[i5][32] + "\t"
								+ noisetable[i5][33] + "\t"
								+ noisetable[i5][34];

						// + "\t" + noisetable[i5][13];
						reslts.add(resline);
					}
				}


				else if (vlns.size()==0) 
				{

					// zeros
					for (int nt = 0; nt < noisetable.length; nt++)
						// for(int nt1=5; nt1<13; nt1++)
						for (int nt1 = 0; nt1 < 33; nt1++)
							noisetable[nt][nt1] = "0";
					for (int nt = 0; nt < noisetable.length; nt++)
						// for(int nt1=5; nt1<13; nt1++)
						for (int nt1 = 33; nt1 < 35; nt1++)

					x=0;
					y=carrs.length;
					for (int i5 = 0; i5 < carrs.length; i5++) 
					{
						
						expected_base = carrs[i5];
						x = i5 + 1;
						y = carrs.length - x + 1;
						curposition = refstart + i5;
						noisetable[i5][0] = CHROM;
						noisetable[i5][1] = Integer.toString(curposition); // site
						noisetable[i5][2] = Integer.toString(x); // x
						noisetable[i5][3] = Integer.toString(y); // y
						noisetable[i5][4] = Character
								.toString(expected_base); // expected
						noisetable[i5][5] = Character
								.toString(expected_base); // alternative if
												// position for
						for(int mls=0; mls<posv.size(); mls++)			// multiple sites
						{
			   				posiii = posv.get(mls);
							//System.out.println(posiii + " <---->  " + curposition); 
			   				altern = altv.get(mls);
							if (posiii == curposition) 
							{
								noisetable[i5][5] = altern;
							//	System.out.println("In");
							//	System.out.println(noisetable[i4][5] + "  --  " + noisetable[i4][4] );
							}
						}
						noisetable[i5][33]=null;
						noisetable[i5][34]=null;			
						resline = noisetable[i5][0] + "\t" + noisetable[i5][1]
								+ "\t" + noisetable[i5][2] + "\t"
								+ noisetable[i5][3] + "\t" + noisetable[i5][4]
								+ "\t" + noisetable[i5][5] + "\t"
								+ noisetable[i5][6] + "\t" + noisetable[i5][7]
								+ "\t" + noisetable[i5][8] + "\t"
								+ noisetable[i5][9] + "\t" + noisetable[i5][10]
								+ "\t" + noisetable[i5][11] + "\t"
								+ noisetable[i5][12] + "\t"
								+ noisetable[i5][13] + "\t"
								+ noisetable[i5][14] + "\t"
								+ noisetable[i5][15] + "\t"
								+ noisetable[i5][16] + "\t"
								+ noisetable[i5][17] + "\t"
								+ noisetable[i5][18] + "\t"
								+ noisetable[i5][19] + "\t"
								+ noisetable[i5][20] + "\t"
								+ noisetable[i5][21] + "\t"
								+ noisetable[i5][22] + "\t"
								+ noisetable[i5][23] + "\t"
								+ noisetable[i5][24] + "\t"
								+ noisetable[i5][25] + "\t"
								+ noisetable[i5][26] + "\t"
								+ noisetable[i5][27] + "\t"
								+ noisetable[i5][28] + "\t"
								+ noisetable[i5][29] + "\t"
								+ noisetable[i5][30] + "\t"
								+ noisetable[i5][31] + "\t"
								+ noisetable[i5][32] + "\t"
								+ noisetable[i5][33] + "\t"
								+ noisetable[i5][34];
						reslts.add(resline);
					}
				}


			}
		}

		File fotnm = new File(output);
		String nm = fotnm.getName();
		// System.out.println("Output ----. " + output);
		// output=output.substring(0, output.length()-4)+ tpogtmm + "_" +
		// output.substring(output.length()-4, output.length());

		writeToFile(output, reslts);

		// System.out.println("[0-1]%, (1-2]%, (2-3]%, (3-4]%, (4-5]%, (5-6]%, (6-7]%, (7-8]%, (8-9]%, (9-10]% >10%");

		pmis.add("[0-1]%	(1-2]%	(2-3]%	(3-4]%	(4-5]%	(5-6]%	(6-7]%	(7-8]%	(8-9]%	(9-10]%	>10%");
		pmis.add(percmis[0] + "\t" + percmis[1] + "\t" + percmis[2] + "\t"
				+ percmis[3] + "\t" + percmis[4] + "\t" + percmis[5] + "\t"
				+ percmis[6] + "\t" + percmis[7] + "\t" + percmis[8] + "\t"
				+ percmis[9] + "\t" + percmis[10]);
		tpogtmm = percmis[10];
		File ofli = new File(output);

		File diri = ofli.getParentFile();
		String oprpth = diri.getName();
		String outlog = oprpth + ".qc";

		writeToFile(output + ".qc", pmis);
		bw1.close();

		// Read Counts Percentage of Mismatches.
        	System.out.println("\nRead Counts Percentage of Mismatches");
        	System.out.println("====================================");
		for(int qci=0; qci<pmis.size(); qci++)
	           System.out.println(pmis.get(qci));

		System.out.println(nm
				+ "\nTotal Number of Reads with > 10% mismatches = " + tpogtmm);
		// Delete the folder in the temporary folder with a unique name to hold
		// temporary data files
		// tmpuiid
		if (dutmp.exists()) {
			dutmp.delete();
			// System.out.println("By Now the folder '" + dutmp.getName() +
			// "'  should have been deleted!");
		}
		/*
		 * } } } } //System.exit(0); //System.out.println("lines size = " +
		 * vlns.size()); // prepare the array numbers that hold string numbers
		 * with zeros for(int nt=0; nt<noisetable.length; nt++) // for(int
		 * nt1=5; nt1<13; nt1++) for(int nt1=0; nt1<35; nt1++)
		 * noisetable[nt][nt1]="0"; if(vlns.size()>0)
		 */
	}

	/**
	 * Method to be documented
	 * 
	 * @param vcseq
	 * @param datql
	 * @param vcposition
	 * @param refstart
	 * @param refsize
	 * @return
	 */
	public static String padDataSequence(String vcseq, String datql,
			int vcposition, int refstart, int refsize, int padi) {
		// constant 450 length : substitute refsize with 450 everywhere in this
		// method
		String ipsq="", ipql="", inpad="", resultdtseq = "", resultedqual = "", reseqqual = "";
		int diff = 0, size = 0, dfgt = 0, dflt = 0, ngtsize = 0;
		size = vcseq.length();
		resultdtseq = vcseq;
		resultedqual = datql;
		if(padi>0)
		{
			for(int p=0; p<padi; p++)
				inpad=inpad+"P";
		}

		if((padi>0) && (2*padi<size))
		{
		 	resultdtseq=inpad + vcseq.substring(padi, vcseq.length()-padi) + inpad;
			resultedqual=inpad + datql.substring(padi, datql.length()-padi) + inpad;
		}

		// System.out.println("fragment  position = " + vcposition);
		// System.out.println("reference position = " + refstart);
		// System.out.println("difference = " + diff);

		if (vcposition == refstart) {
			// if (size<450) // line 455
			if (size < refsize) {
				diff = refsize - size;
				for (int i = 0; i < diff; i++) {
					resultdtseq = resultdtseq + "P";
					resultedqual = resultedqual + "P";
				}
			}
			if (size > refsize) {
				resultdtseq = resultdtseq.substring(0, refsize);
				resultedqual = resultedqual.substring(0, refsize);
			}
		}

		if (vcposition > refstart) {
			// the data starting position is after the reference position : add
			// Ps in front
			diff = vcposition - refstart;
			// System.out.println("fragment  position = " + vcposition);
			// System.out.println("reference position = " + refstart);
			// System.out.println("difference = " + diff);
			for (int i = 0; i < diff; i++) {
				resultdtseq = "P" + resultdtseq;
				resultedqual = "P" + resultedqual;
			}
			size = resultdtseq.length();

			if (size < refsize) {
				diff = refsize - size;
				for (int i = 0; i < diff; i++) {
					resultdtseq = resultdtseq + "P";
					resultedqual = resultedqual + "P";
				}
			}

			if (size > refsize) {
				resultdtseq = resultdtseq.substring(0, refsize);
				resultedqual = resultedqual.substring(0, refsize);
			}
		}

		if (vcposition < refstart) {
			diff = refstart - vcposition;
			// System.out.println("fragment  position = " + vcposition);
			// System.out.println("reference position = " + refstart);
			// System.out.println("difference = " + diff);
			resultdtseq = resultdtseq.substring(diff, resultdtseq.length());
			resultedqual = resultedqual.substring(diff, resultedqual.length());

			size = resultdtseq.length();

			if (resultdtseq.length() < refsize) {
				diff = refsize - resultdtseq.length();
				for (int i = 0; i < diff; i++) {
					resultdtseq = resultdtseq + "P";
					resultedqual = resultedqual + "P";
				}
			}
			if (resultdtseq.length() > refsize) {
				resultdtseq = resultdtseq.substring(0, refsize);
				resultedqual = resultedqual.substring(0, refsize);
			}
		}
		reseqqual = resultdtseq + "\t" + resultedqual;

		return reseqqual;
	}

	public static Vector<String> retrieveBamflNmspths(String dt_pth,
			String slash, int choice) {
		Vector<String> rslts = new Vector<String>();
		// Go through data path and collect all files with 'bam' extension;
		// generate the complete data path
		// in order to be used in process guider by samtools
		Vector<String> fls = new Vector<String>();
		File folder = new File(dt_pth);
		// System.out.println(data_path);
		fls = listFilesForFolder(folder);
		String flnm = "";
		dt_pth = dt_pth + slash;
		if (choice == 1) {
			for (int i = 0; i < fls.size(); i++)
				if (fls.get(i)
						.substring(fls.get(i).length() - 3, fls.get(i).length())
						.equals("bam"))
					rslts.add(dt_pth + fls.get(i));
		}
		if (choice == 0) {
			for (int i = 0; i < fls.size(); i++)
				if (fls.get(i)
						.substring(fls.get(i).length() - 3, fls.get(i).length())
						.equals("bam"))
					rslts.add(fls.get(i));
		}
		rslts = sortStrings(rslts);
		return rslts;
	}

	// Handle Cigar String and adapt read and quality Strings
	public static String HandleCigSs(String cigar, String read, String qual) {
		String result = "";
		String testcigar = "";
		testcigar = retrieveCigarInfoExHs(cigar); // retrieve space separated
													// cigar string with no Hs
													// and their numbers in it

		// System.out.println("testcigar = " + testcigar );
		String[] wcs = testcigar.split(" "); // retrieve words with letters and
												// numbers from transformed(no
												// Hs) and space separated cigar
												// string
		Vector<Integer> cnv = new Vector<Integer>(); // # Matches and Mismatches
														// or # Insertions or #
														// Deletions or Ses all
														// together numeric
														// vector
		Vector<String> clv = new Vector<String>(); // the corresponding to
													// previous(numeric) letter
													// vector
		char[] rch = read.toCharArray(); // transform read to character array(in
											// order to deal with it(delete or
											// insert particular letters)
		char[] qch = qual.toCharArray(); // transform quality string in same way
		// separate and put the elements of cigar space separated string into
		// the vector
		// populate the vectors
		// System.out.println("testcigar = " + testcigar);
		for (int i1 = 0; i1 < wcs.length; i1++) {
			if ((i1 % 2 == 0) && (!wcs[i1].equals(""))) {
				// System.out.println(wcs[i1]);
				cnv.add(Integer.parseInt(wcs[i1].replace("\"", "")));
			} else {
				if (!wcs[i1].equals(""))
					clv.add(wcs[i1]);
			}
		}
		/*
		 * int asize=0; // find the actual size in order to utilize
		 * correspondence with mapping positions
		 * 
		 * for (int i1=0; i1<cnv.size(); i1++) { if ( (clv.get(i1).equals("M"))
		 * || (clv.get(i1).equals("D"))) asize=asize+cnv.get(i1); }
		 * System.out.println(asize);
		 */
		// Loop that creates the final adapted read (no Ss) and final quality
		// string.
		int crps = -1, cp = 0;
		String finalread = "", finalquality = "", finalCigar = "", posi = "";
		int j = 0, bal = 0;// used to increase (positions string in order to
							// keep truck of the positions

		for (int i = 0; i < clv.size(); i++) {
			// System.out.println("----------->  " + clv.get(i));
			// System.out.println("----------->  " + cnv.get(i));

			if (clv.get(i).equals("S")) {
				// System.out.println(clv.get(i) + " : " + cp);

				for (int k = 0; k < cnv.get(i); k++) {
					crps = crps + 1;
				}
				crps = crps + 1;
				bal = 0;
			}

			else if (!(clv.get(i).equals("S")) && !(clv.get(i).equals("D")) && !(clv.get(i).equals("N"))) {
				if ((crps == -1))
					crps = crps + 1;
				cp = crps;
				posi = posi + " ";
				// posi=posi+crps;
				// System.out.println(clv.get(i) + " --:-- " + cp);
				if (cp + cnv.get(i) < rch.length + 1) {
					for (j = cp; j < cp + cnv.get(i); j++) {
						posi = posi + "," + crps;
						finalread = finalread + rch[j];
						//System.out.println("Here 1");
						//System.out.println(cigar);
						//System.out.println(read);
						//System.out.println(qual);

						finalquality = finalquality + qch[j];
						//System.out.println("Here 2");
						crps = crps + 1;
						bal = 1;
					}
					// crps=crps-1;
				}
			}
			// else if((clv.get(i).equals("S")) || (clv.get(i).equals("H")))
		}
		Vector<String> auxclv = new Vector<String>();
		Vector<Integer> auxcnv = new Vector<Integer>();
		Vector<String> aux2clv = new Vector<String>(); 
		// these vectors work with previous auxiliary and are used to 
		// generate the final cigar string (no Ss)
		Vector<Integer> aux2cnv = new Vector<Integer>();
		// Loop that populates auxiliary vectors for numbers and Letters of the
		// cigar string.
		String finalcigar = "";
		for (int i = 0; i < clv.size(); i++) {
			if (!(clv.get(i).equals("S")) && !(clv.get(i).equals("H"))) {
				auxclv.add(clv.get(i));
				auxcnv.add(cnv.get(i));
				finalcigar = finalcigar + cnv.get(i) + clv.get(i);
			}
		}
		//
		// System.out.println(auxclv);
		// System.out.println(auxcnv);
		// Loop that generates the final cigar string vector components
		String prv = "", cur = "";
		int prvnm = 0, curnm = 0, sum = 0;
		//
		for (int i = 0; i < auxclv.size(); i++) {
			cur = auxclv.get(i);
			curnm = auxcnv.get(i);
			if (i > 0)// means there must be previous element
			{
				if (cur.equals(prv)) {
					sum = sum + curnm;
				}
				if (!(cur.equals(prv))) {
					aux2clv.add(prv);
					aux2cnv.add(sum);
					sum = 0;
				}
				if (!(cur.equals(prv)) && (i == auxclv.size() - 1)) {
					aux2clv.add(cur);
					aux2cnv.add(curnm);
				}
				if ((cur.equals(prv)) && (i == auxclv.size() - 1)) {
					aux2clv.add(cur);
					aux2cnv.add(sum);
				}
			}
			prv = auxclv.get(i);
			prvnm = auxcnv.get(i);
			// System.out.println("Sum = " + sum);
			if (sum == 0) {
				sum = prvnm;
				// System.out.println("after Sum 0 =  " + sum);
			}
		}

		String finalcigars = "";
		// Loop that generates the final cigar string from aux2 vector.
		for (int i = 0; i < aux2clv.size(); i++)
			finalcigars = finalcigars + aux2cnv.get(i) + aux2clv.get(i);

		// System.out.println("finalcigar = " + finalcigar);
		// System.out.println("finalcigars = " + finalcigars);
		result = finalcigar + "\n" + finalread + "\n" + finalquality;
		/*
		 * System.out.println(cigar); System.out.println(finalcigar);
		 * System.out.println(finalcigars); System.out.println("R " + read);
		 * System.out.println("f " + finalread); // System.out.println("r " +
		 * reff); System.out.println(posi);
		 */
		return result;
	}

	public static String retrieveCigarInfo(String cigar) {
		// System.out.println( Start of 'retrieveCigarInfo' method ..........");
		String res = "";
		// System.out.println(cigar);
		char[] cigchArray = cigar.toCharArray();
		// System.out.println("----hi " + cigchArray[0]);
		String numst = "";
		Vector<Integer> cignumbs = new Vector<Integer>();
		Vector<String> cigRes = new Vector<String>();

		for (int i = 0; i < cigchArray.length; i++) {
			if (Character.isDigit(cigchArray[i])) {
				numst = numst + cigchArray[i];
				// System.out.println(" " + Character.toString(cigchArray[i]));
			}
			if (!Character.isDigit(cigchArray[i])) {
				cignumbs.add(Integer.parseInt(numst));
				numst = "";
				cigRes.add(Character.toString(cigchArray[i]));
			}
		}

		// System.out.println(cigar);
		// System.out.println(cignumbs);
		// System.out.println(cigRes);
		for (int i = 0; i < cignumbs.size(); i++) {
			if (i == 0)
				res = cignumbs.get(i) + " " + cigRes.get(i);
			if (i > 0)
				res = res + " " + cignumbs.get(i) + " " + cigRes.get(i);
		}
		// System.out.println(........End of 'retrieveCigarInfo' method!");
		return res;
	}

	/**
	 *
	 * return res
	 */
	public static String retrieveCigarInfoExHs(String cigar) {
		// System.out.println( Start of 'retrieveCigarInfo' method ..........");
		String res = "";
		// System.out.println(cigar);
		// System.out.println("cigar : " + cigar);
		char[] cigchArray = cigar.toCharArray();
		// System.out.println("----hi " + cigchArray[0]);
		String numst = "";
		Vector<Integer> cignumbs = new Vector<Integer>();
		Vector<String> cigRes = new Vector<String>();

		for (int i = 0; i < cigchArray.length; i++) {
			if (Character.isDigit(cigchArray[i])) {
				numst = numst + cigchArray[i];
				// System.out.println(" " + Character.toString(cigchArray[i]));
			}
			// exclude Hs and their number
			if (((!Character.isDigit(cigchArray[i])) && (cigchArray[i]) != 'H')) {
				cignumbs.add(Integer.parseInt(numst));
				numst = "";
				cigRes.add(Character.toString(cigchArray[i]));
			}

			if (((!Character.isDigit(cigchArray[i])) && (cigchArray[i]) == 'H')) {
				numst = "";
			}
		}

		// the modified loop for handling Hs
		res = "";
		for (int i = 0; i < cignumbs.size(); i++) {
			// if(!(cigRes.get(i).equals("H")))
			// {
			if (i == 0)
				res = cignumbs.get(i) + " " + cigRes.get(i);
			if (i > 0)
				res = res + " " + cignumbs.get(i) + " " + cigRes.get(i);
			// }
		}

		// System.out.println(........End of 'retrieveCigarInfo' method!");
		return res;
	}

	// Handle Cigar String and adapt read and quality Strings
	// Handle Cigar String and adapt read and quality Strings
	// Handle Cigar String and adapt read and quality Strings
	// Handle Cigar String and adapt read and quality Strings
	// Handle Cigar String and adapt read and quality Strings
	/**
	 * Active method that needs to be doccumented when the last check for
	 * noisetab command will be performed Pay attention to the commented
	 * examples in the method.
	 * 
	 * 
	 * @param res1data
	 * @return
	 **/



	public static Vector<String> PrepareTheDataAndInDlsBack111215(Vector<String> res1data) {
		System.out.println("	Start of 'PrepareTheData' Method!");
		System.out.println("		Processing Start..................");
		System.out.println("----> Size = " + res1data.size());
		int regectread = 0;
		/*
		 * String example_line = + origflnm + "\t" // 0 + site +"\t" // 1 +
		 * CHROM + "\t" // 2 + REF + "\t" // 3 + ALT + "\t" // 4 + chromosome +
		 * "\t" // 5 + position + "\t" // 6 + MAPQ + "\t" // 7 + cigar + "\t" //
		 * 8 + sequence + "\t" // 9 + qual + "\t" // 10 + strand + "\t" // 11 +
		 * readname + "\t" // 12 + numbelines; // 13
		 */
		String line = "", cigar = "", fragment = "", site_ref = "", filenm = "", quality = "", strand = "", readname = "", readnmreturn = "", numberlines = "";
		int chrom = 0, site = 0, position = 0, MAPQ = 0, expval = 0, frsize = 0;
		int asize = 0;
		String auxstr = "";
		String ciginfo = "";
		String fifrag = "";
		String fiqual = "";
		Vector<Integer> cnv = new Vector<Integer>(); // populate 'C'igar
														// 'N'umbers 'V'ector
														// with numbers
		Vector<String> clv = new Vector<String>(); // populate 'C'igar 'L'etters
													// 'V'ector with letters
		Vector<String> ffrag = new Vector<String>(); // final transformed
														// fragment
		Vector<String> fqual = new Vector<String>();
		Vector<Integer> delpoints = new Vector<Integer>();
		Vector<Integer> inspoints = new Vector<Integer>();
		int Ds = 0;

		Vector<String> res1fnms = new Vector<String>(); // results 1 file names
		Vector<String> res2fnms = new Vector<String>(); // results 2 file names
		
		Vector<String> secres = new Vector<String>();
		String prvcig = "";
		String rd = "", rd1 = "";
		String gic = "", gic1 = "";
		String chromS="";
		for (int i = 0; i < res1data.size(); i++)
		// for(int i=0; i<10; i++)
		// for(int i=1890; i<1914; i++)
		{
			regectread = 0;
			line = res1data.get(i);
			// System.out.println("\n line = " + line);
			String mwords[] = line.split("\t");
			// for(int y=0; y<mwords.length; y++)
			for (int y = 0; y < 10; y++) {
				mwords[y] = mwords[y].replace("\"", "");
			}

			/*
			 * + origflnm + "\t" // 0 + site +"\t" // 1 + CHROM + "\t" // 2 +
			 * REF + "\t" // 3 + ALT + "\t" // 4 + chromosome + "\t" // 5 +
			 * position + "\t" // 6 + MAPQ + "\t" // 7 + cigar + "\t" // 8 +
			 * sequence + "\t" // 9 + qual + "\t" // 10 + strand + "\t" // 11 +
			 * readname + "\t" // 12 + numbelines; // 13
			 */
			filenm = mwords[0]; // filename
			readname = mwords[12];
			site = Integer.parseInt(mwords[1].replace("\"", ""));
			site_ref = mwords[2].replace("\"", "");
			// expval=Integer.parseInt(mwords[12].replace("\"", ""));
			
			// chrom = Integer.parseInt(mwords[5].replace("\"", ""));
			chromS = (mwords[5].replace("\"", ""));
			position = Integer.parseInt(mwords[6].replace("\"", ""));
			MAPQ = Integer.parseInt(mwords[7].replace("\"", ""));
			cigar = mwords[8];
			gic = mwords[8];
			cigar = cigar.replace("\"", "");
			rd = mwords[9];
			if (cigar.length() <= 1)
				regectread = 1;
			if (regectread == 0) {
				prvcig = cigar;
				fragment = mwords[9];
				fragment = fragment.replace("\"", "");

				frsize = fragment.length();
				quality = mwords[10];
				// base quality field contains characters that are quality
				// indicators and always are important (thus this field must not
				// be cleaned
				// quality = quality.replace("\"", "");
				strand = mwords[11];
				// System.out.println("File Name : " + filenm +
				// "\nChromosome : chr" + chrom + "\nPosition : " + position +
				// "\nSeqSize : " + frsize + "\nSite : " + site +
				// "\nSite Reference : " + site_ref + "\nExpected Value : " +
				// expval + "\nnCigar : " + nCigar + "\nCigar : " + cigar);
				// System.out.println(cigar);
				// for(int w=0; w<mwords.length; w++)
				// System.out.println(mwords[w]);

				String rescireadqual = "";
				rescireadqual = HandleCigSs(cigar, fragment, quality);
				String[] cfq = rescireadqual.split("\n");
				String curcig = cfq[0];
				cigar = cfq[0];
				gic1 = cfq[0];
				fragment = cfq[1];
				rd1 = cfq[1];
				quality = cfq[2];
				ciginfo = retrieveCigarInfo(cigar); // retrieve a string that
													// contains digestible cigar
													// info
				// example : "123 M 1 I 25 M 2 D 2 M" (elements separated by
				// space; odds are letters; evens are numbers.)
				String[] wcs = ciginfo.split(" "); // retrieve words with
													// letters and numbers
													// separately

				// System.out.println(ciginfo);
				// System.out.println("Chromosome : chr" + chrom);
				cnv = new Vector<Integer>(); // # matches or # mismatches or #
												// Insertions or # deletions all
												// together numeric vector
				clv = new Vector<String>(); // the corresponding to
											// previous(numeric) letter vector

				// separate and put the elements of cigar space separated string
				// into the vector
				for (int i1 = 0; i1 < wcs.length; i1++) {
					if ((i1 % 2 == 0) && (!wcs[i1].equals(""))) {
						// System.out.println(wcs[i1]);
						cnv.add(Integer.parseInt(wcs[i1].replace("\"", "")));
					} else {
						if (!wcs[i1].equals(""))
							clv.add(wcs[i1]);
					}
				}

				// find the actual size
				asize = 0;
				// System.out.println(cnv.size());
				// for loop that calculates the length of the read that
				// finally will be created by adding the deletion points too
				// This length might be greater than the initial length if
				// insertion weren't found but deletions exist.
				// Example 1 Cigar = 20M2D39M, initial length = 59, final
				// length=61 (final cigar = 20M2M39M=61M)

				// Example 2 Cigar = 20M1D1I39M, initial length = 60, final
				// length = 60
				// (one insertion was removed and one D have been added)
				// (final cigar = 20M1D39M=60M since information about Ds
				// included in the read now)

				for (int i1 = 0; i1 < cnv.size(); i1++) {
					if ((clv.get(i1).equals("M")) || (clv.get(i1).equals("D")))
						asize = asize + cnv.get(i1);
				}

				String[] fragpl = fragment.split(""); // fragment string
														// character array
				String[] qualpl = quality.split(""); // quality string character
														// array

				ffrag = new Vector<String>(); // final transformed fragment
				ffrag.addAll(Arrays.asList(fragpl));
				ffrag.remove(0);

				fqual = new Vector<String>(); // final transformed fragment
				fqual.addAll(Arrays.asList(qualpl));
				fqual.remove(0);

				// System.out.println(" Hello  - > : " + ffrag.get(0));
				// System.out.println(" Hello  - >   : " +
				// ffrag.get(ffrag.size()-1));
				// find and store all insertion and deletion points in
				// corresponding vectors

				String delsstr = "", inssstr = "";
				delpoints = new Vector<Integer>();
				Ds = 0; // Ds counter
				inspoints = new Vector<Integer>();
				int fsize = 0; // fragment size (number of total bases Ms and Ds
				int dis = 0;

				Vector<String> dels = new Vector<String>();
				String delstr = "";

				Vector<String> inss = new Vector<String>();
				String insstr = "", insbs = "";
				int currpos = 0; // Variable used to deal with the actual read
									// string char array]
									// it is used to target the info about
									// insertions and retrieve the involved
									// bases
									// it is incremented only with presence of
									// base. Deletion cannot increment that
									// variable.

				for (int i2 = 0; i2 < clv.size(); i2++) {
					if (clv.get(i2).equals("M")) {
						fsize = fsize + cnv.get(i2);
						currpos = currpos + cnv.get(i2); // this is very useful
															// for Insertions
						// currpos is an indicator of before insertion position
						// that depends on the
						// initial size of the read. the final (and mapping
						// position are retrieved form 'fize' variable
					}
					if (clv.get(i2).equals("D")) // if cigar indicates D then
													// retrieve its number and
					// go through a loop in it
					{
						delstr = (position + fsize - 1) + "_" + cnv.get(i2);
						dels.add(delstr);
						for (dis = 1; dis <= cnv.get(i2); dis++) {
							delpoints.add(fsize + dis - 1 + Ds);
							// System.out.print("\n->" + (fsize+dis-1 + Ds)
							// +"<-\n");
						}
						Ds = Ds + dis - 1;
						fsize = fsize + cnv.get(i2);

						// fsize=fsize+cnv.get(i);
						// System.out.println("Deletion point will not count in fragment size for now.");
						// System.out.println(delpoints);
					}
					// System.out.println(fsize);
					if (clv.get(i2).equals("I")) {

						// for(int ins=1; ins<=cnv.get(i2); ins++)
						// inspoints.add(fsize+ins-1);
						insbs = "";
						for (int j = currpos; j < currpos + cnv.get(i2); j++) {
							insbs = insbs + ffrag.get(j);
							inspoints.add(j);
							// ffrag.remove(j);
							// j=j-1;
						}
						insstr = (position + fsize - 1) + "_" + cnv.get(i2)
								+ "_" + insbs;
						inss.add(insstr);
						// fsize=fsize+cnv.get(i2);
						currpos = currpos + cnv.get(i2);
					}

				}

				for (int i3 = 0; i3 < inspoints.size(); i3++) {
					ffrag.set(inspoints.get(i3), "I");
					// System.out.print(inspoints.get(i3)+ " ");
				}

				int ICounter = 0;
				// System.out.println();
				for (int i3 = 0; i3 < ffrag.size(); i3++) {
					if (ffrag.get(i3).equals("I")) {
						ffrag.remove(i3);
						fqual.remove(i3);
						// System.out.println("'I' detected!  : " + i3);
						if (i3 > 1)
							i3 = i3 - 1; // Vectors change dynamically their
											// size
						// in 'element removal' the next element will acquire
						// the current index
						// the current index must be then decreased in order for
						// that element to
						// be able to be checked with the next iteration.
						ICounter = ICounter + 1;
					}
				}
				/*
				 * System.out.println(ffrag); System.out.println(fqual);
				 * System.out.println(delpoints); System.out.println(fragment);
				 * System.out.println(rd); System.out.println(rd1);
				 * System.out.println(gic); System.out.println(gic1);
				 */
				for (int i3 = 0; i3 < delpoints.size(); i3++) {
					// ffrag.add(delpoints.get(i), "D" + i);
					ffrag.add(delpoints.get(i3) - i3, "D");
					fqual.add(delpoints.get(i3) - i3, "D");
				}
				fifrag = "";
				for (int i3 = 0; i3 < ffrag.size(); i3++)
					fifrag = fifrag + ffrag.get(i3);

				fiqual = "";
				for (int i3 = 0; i3 < fqual.size(); i3++)
					fiqual = fiqual + fqual.get(i3);

				delsstr = "";
				// returned read name will be cmn = common if no deletions and
				// no insertions have been detected.
				readnmreturn = "cmn";
				if (dels.size() > 0) {
					for (int in = 0; in < dels.size(); in++) {
						if (in < dels.size() - 1)
							delsstr = delsstr + dels.get(in) + ",";
						if (in == dels.size() - 1)
							delsstr = delsstr + dels.get(in);
					}
					// returned read name will be cmn = common if no deletions
					// and no insertions have been detected.
					readnmreturn = readname;
				}
				if (dels.size() == 0)
					delsstr = "none";

				inssstr = "";
				if (inss.size() > 0) {
					for (int in = 0; in < inss.size(); in++) {
						if (in < inss.size() - 1)
							inssstr = inssstr + inss.get(in) + ",";
						if (in == inss.size() - 1)
							inssstr = inssstr + inss.get(in);

					}
					// returned read name will be cmn = common if no deletions
					// and no insertions have been detected.
					readnmreturn = readname;
				}
				if (inss.size() == 0)
					inssstr = "none";
// 03/06/2015 Changed to accept X and Y chromosome (chrom to chromS)
				secres.add(chromS + "\t" + (position) + "\t"
// 03/06/2015 Changed to accept X and Y chromosome (chrom to chromS)
						+ (position + fifrag.length() - 1) + "\t" + fifrag
						+ "\t" + fiqual + "\t" + delsstr + "\t" + inssstr
						+ "\t" + cigar + "\t" + curcig + "\t" + prvcig + "\t"
						+ readnmreturn + "\t" + filenm);
			}
			// Position already includes the first place of length therefore to
			// find read end actual position add length and subtract one
		}

		System.out.println("................Preparation  Ended.");
		System.out.println("	End of 'PrepareTheData' Method!");
		return secres;
	}


	public static Vector<String> PrepareTheDataAndInDls(Vector<String> res1data) {
		System.out.println("	Start of 'PrepareTheData' Method!");
		System.out.println("		Processing Start..................");
		System.out.println("    Size = " + res1data.size());
		int regectread = 0;
		//
		// String example_line = + origflnm + "\t" // 0 + site +"\t" // 1 +
		// CHROM + "\t" // 2 + REF + "\t" // 3 + ALT + "\t" // 4 + chromosome +
		// "\t" // 5 + position + "\t" // 6 + MAPQ + "\t" // 7 + cigar + "\t" //
		// 8 + sequence + "\t" // 9 + qual + "\t" // 10 + strand + "\t" // 11 +
		// readname + "\t" // 12 + numbelines; // 13
		//
		String line = "", cigar = "", fragment = "", 
				site_ref = "", filenm = "", quality = "", 
				strand = "", readname = "", readnmreturn = "", 
				numberlines = "", delenstr="", prvcig = "", 
				rd = "", rd1 = "", gic = "", gic1 = "", chromS = "",
				rescireadqual = "", curcig = "",
				auxstr = "", ciginfo = "", fifrag = "", fiqual = "",
				insstr = "", insbs = "";
		
		int chrom = 0, site = 0, position = 0, MAPQ = 0, 
				expval = 0, frsize = 0, asize = 0, DNs=0, currpos=0;
	
		Vector<Integer> cnv = new Vector<Integer>(); // populate 'C'igar// 'N'umbers 'V'ector// with numbers
		Vector<String>  clv = new Vector<String>(); // populate 'C'igar 'L'etter// 'V'ector with letters
		Vector<String> ffrag = new Vector<String>(); // finaltransformed// fragment
		Vector<String> fqual = new Vector<String>();
		Vector<Integer> delpoints = new Vector<Integer>();
		Vector<Integer> delsenspoints = new Vector<Integer>();
		Vector<String> delsensletters = new Vector<String>();
		Vector<Integer> inspoints = new Vector<Integer>();
		Vector<String> secres = new Vector<String>();
		Vector<String> delsEns = new Vector<String>();
		Vector<String> dels = new Vector<String>();
		Vector<String> inss = new Vector<String>();
		
		
		for (int i = 0; i < res1data.size(); i++)
		{
			regectread = 0;
			line = res1data.get(i);
			// System.out.println("\n line = " + line);
			String mwords[] = line.split("\t");
			// for(int y=0; y<mwords.length; y++)
			for (int y = 0; y < 10; y++) {
				mwords[y] = mwords[y].replace("\"", "");
			}

			//
			// + origflnm + "\t" // 0 + site +"\t" // 1 + CHROM + "\t" // 2 +
			// REF + "\t" // 3 + ALT + "\t" // 4 + chromosome + "\t" // 5 +
			// position + "\t" // 6 + MAPQ + "\t" // 7 + cigar + "\t" // 8 +
			// sequence + "\t" // 9 + qual + "\t" // 10 + strand + "\t" // 11 +
			// readname + "\t" // 12 + numbelines; // 13
			//
			filenm = mwords[0]; // filename
			readname = mwords[12];
			site = Integer.parseInt(mwords[1].replace("\"", ""));
			site_ref = mwords[2].replace("\"", "");
			// expval=Integer.parseInt(mwords[12].replace("\"", ""));

			// chrom = Integer.parseInt(mwords[5].replace("\"", ""));
			chromS = (mwords[5].replace("\"", ""));
			position = Integer.parseInt(mwords[6].replace("\"", ""));
			MAPQ = Integer.parseInt(mwords[7].replace("\"", ""));
			cigar = mwords[8];
			gic = mwords[8];
			cigar = cigar.replace("\"", "");
			for(int yo=0;yo<cigar.length(); yo++)
				if(cigar.charAt(yo)=='N')
					System.out.println("Found N in cigar : " + cigar);
			rd = mwords[9];
			quality = mwords[10];
			if ((cigar.length() <= 1) || (quality.length() <= 1) )
				regectread = 1;

			if (regectread == 0) 
			{
				prvcig = cigar;
				fragment = mwords[9];
				fragment = fragment.replace("\"", "");

				frsize = fragment.length();
				quality = mwords[10];
				// base quality field contains characters that are quality
				// indicators and always are important (thus this field must not
				// be cleaned
				// quality = quality.replace("\"", "");
				strand = mwords[11];
				// System.out.println("File Name : " + filenm +
				// "\nChromosome : chr" + chrom + "\nPosition : " + position +
				// "\nSeqSize : " + frsize + "\nSite : " + site +
				// "\nSite Reference : " + site_ref + "\nExpected Value : " +
				// expval + "\nnCigar : " + nCigar + "\nCigar : " + cigar);
				// System.out.println(cigar);
				// for(int w=0; w<mwords.length; w++)
				// System.out.println(mwords[w]);
				rescireadqual = HandleCigSs(cigar, fragment, quality);
				String[] cfq = rescireadqual.split("\n");
				curcig = cfq[0];
				cigar = cfq[0];
				gic1 = cfq[0];
				fragment = cfq[1];
				rd1 = cfq[1];
				quality = cfq[2];
				ciginfo = retrieveCigarInfo(cigar); // retrieve a string that
													// contains digestible cigar
													// info
				// example : "123 M 1 I 25 M 2 D 2 M" (elements separated by
				// space; odds are letters; evens are numbers.)
				String[] wcs = ciginfo.split(" "); // retrieve words with
													// letters and numbers
				dels = new Vector<String>();								// separately
				delsenspoints = new Vector<Integer>();
				// System.out.println(ciginfo);
				// System.out.println("Chromosome : chr" + chrom);
				cnv = new Vector<Integer>(); // # matches or # mismatches or #
				// Insertions or # deletions all
				// together numeric vector
				clv = new Vector<String>(); // the corresponding to
				// previous(numeric) letter vector
				// separate and put the elements of cigar space separated string
				// into the vector
				for (int i1 = 0; i1 < wcs.length; i1++)
				{
					if ((i1 % 2 == 0) && (!wcs[i1].equals("")))
					{
						// System.out.println(wcs[i1]);
						cnv.add(Integer.parseInt(wcs[i1].replace("\"", "")));
					}
					else
					{
						if (!wcs[i1].equals(""))
							clv.add(wcs[i1]);
					}
				}

				// find the actual size
				asize = 0;
				// System.out.println(cnv.size());
				// for loop that calculates the length of the read that
				// finally will be created by adding the deletion points too
				// This length might be greater than the initial length if
				// insertion weren't found but deletions exist.
				// Example 1 Cigar = 20M2D39M, initial length = 59, final
				// length=61 (final cigar = 20M2M39M=61M)
				// Example 2 Cigar = 20M1D1I39M, initial length = 60, final
				// length = 60
				// (one insertion was removed and one D have been added)
				// (final cigar = 20M1D39M=60M since information about Ds
				// included in the read now)

//				Here the Ns must be added
				for (int i1 = 0; i1 < cnv.size(); i1++) 
				{
					if ((clv.get(i1).equals("M")) 
					|| (clv.get(i1).equals("D")) 
					|| (clv.get(i1).equals("N")))
						asize = asize + cnv.get(i1);
				}

				String[] fragpl = fragment.split(""); // fragment string
														// character array
				String[] qualpl = quality.split(""); // quality string character
														// array

				ffrag = new Vector<String>(); // final transformed fragment
				ffrag.addAll(Arrays.asList(fragpl));
				ffrag.remove(0);

				fqual = new Vector<String>(); // final transformed fragment
				fqual.addAll(Arrays.asList(qualpl));
				fqual.remove(0);

				String delsstr = "", ensstr="", inssstr = "";
				delpoints = new Vector<Integer>();
				
				DNs = 0; // Ds counter
				inspoints = new Vector<Integer>();
				int fsize = 0; // fragment size (number of total bases Ms and Ds
				int dnis = 0;
				String delstr = "";
				String enstr = "";
				//enspoints = new Vector<Integer>();
				
				inss = new Vector<String>();
				insstr = ""; insbs = "";
				currpos = 0; // Variable used to deal with the actual read
									// string char array]
									// it is used to target the info about
									// insertions and retrieve the involved
									// bases
									// it is incremented only with presence of
									// base. Deletion cannot increment that
									// variable.

				for (int i2 = 0; i2 < clv.size(); i2++)
				{
					if (clv.get(i2).equals("M"))
					{
						fsize = fsize + cnv.get(i2);
						currpos = currpos + cnv.get(i2); // this is very useful
															// for Insertions
						// currpos is an indicator of before insertion position
						// that depends on the
						// initial size of the read. the final (and mapping
						// position are retrieved form 'fize' variable
					}
					if ( (clv.get(i2).equals("D")) ||  (clv.get(i2).equals("N")) )// if cigar indicates D then
													// retrieve its number and
					// go through a loop in it
					{

						if((clv.get(i2).equals("D")) || clv.get(i2).equals("N"))
						{
							delenstr = (position + fsize - 1) + "_" + cnv.get(i2);
							if(clv.get(i2).equals("D"))
								dels.add(delenstr);
							delsEns.add(delenstr);
						}
				
						for (dnis = 1; dnis <= cnv.get(i2); dnis++) 
						{
							//delpoints.add(fsize + dnis - 1 + Ds);
							delsenspoints.add(fsize + dnis - 1 + DNs);
							if((clv.get(i2).equals("D")))
							{
								//System.out.println("DDD " + (fsize + dnis - 1 + DNs));
								delsensletters.add("D");
								delpoints.add(fsize + dnis - 1 + DNs);
							}
							if((clv.get(i2).equals("N")))
								delsensletters.add("d");
								// System.out.print("\n->" + (fsize+dis-1 + Ds)
								// +"<-\n");
						}
						DNs = DNs + dnis - 1;
						fsize = fsize + cnv.get(i2);

						// fsize=fsize+cnv.get(i);
						// System.out.println("Deletion point will not count in fragment size for now.");
						// System.out.println(delpoints);
					}
					// System.out.println(fsize);
					
					if (clv.get(i2).equals("I")) 
					{
						// for(int ins=1; ins<=cnv.get(i2); ins++)
						// inspoints.add(fsize+ins-1);
						insbs = "";
						for (int j = currpos; j < currpos + cnv.get(i2); j++)
						{
							insbs = insbs + ffrag.get(j);
							inspoints.add(j);
							// ffrag.remove(j);
							// j=j-1;
						}
						insstr = (position + fsize - 1) + "_" + cnv.get(i2)
								+ "_" + insbs;
						inss.add(insstr);
						// fsize=fsize+cnv.get(i2);
						currpos = currpos + cnv.get(i2);
					}
				}
				
//				put insertions(Is) on read using the insertion stored points(indeces)
//				and then remove those Is and their mirrored qualities using the inserted Ises' indeces
				for (int i3 = 0; i3 < inspoints.size(); i3++)
				{
					ffrag.set(inspoints.get(i3), "I");
					// System.out.print(inspoints.get(i3)+ " ");
				}

				int ICounter = 0;
				// System.out.println();

				for (int i3 = 0; i3 < ffrag.size(); i3++)
				{
					if (ffrag.get(i3).equals("I"))
					{
						ffrag.remove(i3);
						fqual.remove(i3);
						// System.out.println("'I' detected!  : " + i3);
						if (i3 > 1)
							i3 = i3 - 1; // Vectors change dynamically their size
						// in 'element removal' the next element will acquire
						// the current index
						// the current index must be then decreased in order for
						// that element to
						// be able to be checked with the next iteration.
						ICounter = ICounter + 1;
					}
				}
//				Now quality and frag strings have no insertions on them 
				
				//
				// System.out.println(ffrag); System.out.println(fqual);
				// System.out.println(delpoints); System.out.println(fragment);
				// System.out.println(rd); System.out.println(rd1);
				// System.out.println(gic); System.out.println(gic1);
				//

//				
//				for (int i3 = 0; i3 < delpoints.size(); i3++)
//				{
//					// ffrag.add(delpoints.get(i), "D" + i);
//					ffrag.add(delpoints.get(i3) - i3, "D");
//					fqual.add(delpoints.get(i3) - i3, "D");
//				}
//
				// System.out.println(delsenspoints.size());
				for (int i3 = 0; i3 < delsenspoints.size(); i3++)
				{
					// ffrag.add(delpoints.get(i), "D" + i);
					if(delsensletters.get(i3).equals("D"))
					{
						//System.out.println(ffrag.size() + "  " + (delsenspoints.get(i3) -i3));
						//System.out.println(cigar);
						ffrag.add(delsenspoints.get(i3) - i3, "D");
						fqual.add(delsenspoints.get(i3) - i3, "D");
					}
					else if(delsensletters.get(i3).equals("d"))
					{
						ffrag.add(delsenspoints.get(i3) - i3, "d");
						fqual.add(delsenspoints.get(i3) - i3, "d");
					}
				}
				// System.out.println(ffrag);
				// System.out.println(fqual);
		
		
				fifrag = "";
				for (int i3 = 0; i3 < ffrag.size(); i3++)
					fifrag = fifrag + ffrag.get(i3);

				fiqual = "";			
				for (int i3 = 0; i3 < fqual.size(); i3++)
					fiqual = fiqual + fqual.get(i3);

				delsstr = "";
				// returned read name will be cmn = common if no deletions and
				// no insertions have been detected.
				readnmreturn = "cmn";
				if (dels.size() > 0)
				{
					for (int in = 0; in < dels.size(); in++)
					{
						if (in < dels.size() - 1)
							delsstr = delsstr + dels.get(in) + ",";
						if (in == dels.size() - 1)
							delsstr = delsstr + dels.get(in);
					}
					// returned read name will be cmn = common if no deletions
					// and no insertions have been detected.
					readnmreturn = readname;
				}
				if (dels.size() == 0)
					delsstr = "none";

				inssstr = "";
				if (inss.size() > 0)
				{
					for (int in = 0; in < inss.size(); in++)
					{
						if (in < inss.size() - 1)
							inssstr = inssstr + inss.get(in) + ",";
						if (in == inss.size() - 1)
							inssstr = inssstr + inss.get(in);
					}
					// returned read name will be cmn = common if no deletions
					// and no insertions have been detected.
					readnmreturn = readname;
				}
				if (inss.size() == 0)
					inssstr = "none";
				// 03/06/2015 Changed to accept X and Y chromosome (chrom to
				// chromS)
				secres.add(chromS + "\t"
						+ (position)
						+ "\t"
						// 03/06/2015 Changed to accept X and Y chromosome
						// (chrom to chromS)
						+ (position + fifrag.length() - 1) + "\t" + fifrag
						+ "\t" + fiqual + "\t" + delsstr + "\t" + inssstr
						+ "\t" + cigar + "\t" + curcig + "\t" + prvcig + "\t"
						+ readnmreturn + "\t" + filenm);
			}
			// Position already includes the first place of length therefore to
			// find read end actual position add length and subtract one
		}

		System.out.println("................Preparation  Ended.");
		System.out.println("	End of 'PrepareTheData' Method!");
		return secres;
	}

	
	public static Vector<String> PrepareTheDataAndInDlsEx(Vector<String> res1data) {
		System.out.println("	Start of 'PrepareTheData' Method!");
		System.out.println("		Processing Start..................");
		System.out.println("    Size = " + res1data.size());
		int regectread = 0;
		//
		// String example_line = + origflnm + "\t" // 0 + site +"\t" // 1 +
		// CHROM + "\t" // 2 + REF + "\t" // 3 + ALT + "\t" // 4 + chromosome +
		// "\t" // 5 + position + "\t" // 6 + MAPQ + "\t" // 7 + cigar + "\t" //
		// 8 + sequence + "\t" // 9 + qual + "\t" // 10 + strand + "\t" // 11 +
		// readname + "\t" // 12 + numbelines; // 13
		//
		String line = "", cigar = "", fragment = "", 
				site_ref = "", filenm = "", quality = "", 
				strand = "", readname = "", readnmreturn = "", 
				numberlines = "", delenstr="", prvcig = "", 
				rd = "", rd1 = "", gic = "", gic1 = "", chromS = "",
				rescireadqual = "", curcig = "",
				auxstr = "", ciginfo = "", fifrag = "", fiqual = "",
				insstr = "", insbs = "";
		
		int chrom = 0, site = 0, position = 0, MAPQ = 0, 
				expval = 0, frsize = 0, asize = 0, DNs=0, currpos=0;
	
		Vector<Integer> cnv = new Vector<Integer>(); // populate 'C'igar// 'N'umbers 'V'ector// with numbers
		Vector<String>  clv = new Vector<String>(); // populate 'C'igar 'L'etter// 'V'ector with letters
		Vector<String> ffrag = new Vector<String>(); // finaltransformed// fragment
		Vector<String> fqual = new Vector<String>();
		Vector<Integer> delpoints = new Vector<Integer>();
		Vector<Integer> delsenspoints = new Vector<Integer>();
		Vector<String> delsensletters = new Vector<String>();
		Vector<Integer> inspoints = new Vector<Integer>();
		Vector<String> secres = new Vector<String>();
		Vector<String> delsEns = new Vector<String>();
		Vector<String> dels = new Vector<String>();
		Vector<String> inss = new Vector<String>();
		
		String index="";
		
		for (int i = 0; i < res1data.size(); i++)
		{
			regectread = 0;
			line = res1data.get(i);
			// System.out.println("\n line = " + line);
			String mwords[] = line.split("\t");
			// for(int y=0; y<mwords.length; y++)
			for (int y = 0; y < 10; y++) {
				mwords[y] = mwords[y].replace("\"", "");
			}

			//
			// + origflnm + "\t" // 0 + site +"\t" // 1 + CHROM + "\t" // 2 +
			// REF + "\t" // 3 + ALT + "\t" // 4 + chromosome + "\t" // 5 +
			// position + "\t" // 6 + MAPQ + "\t" // 7 + cigar + "\t" // 8 +
			// sequence + "\t" // 9 + qual + "\t" // 10 + strand + "\t" // 11 +
			// readname + "\t" // 12 + numbelines; // 13 // index
			//
			filenm = mwords[0]; // filename
			readname = mwords[12];
			site = Integer.parseInt(mwords[1].replace("\"", ""));
			site_ref = mwords[2].replace("\"", "");
			// expval=Integer.parseInt(mwords[12].replace("\"", ""));

			// chrom = Integer.parseInt(mwords[5].replace("\"", ""));
			chromS = (mwords[5].replace("\"", ""));
			position = Integer.parseInt(mwords[6].replace("\"", ""));
			MAPQ = Integer.parseInt(mwords[7].replace("\"", ""));
			cigar = mwords[8];
			gic = mwords[8];
			cigar = cigar.replace("\"", "");
			for(int yo=0;yo<cigar.length(); yo++)
				if(cigar.charAt(yo)=='N')
					System.out.println("Found N in cigar : " + cigar);
			rd = mwords[9];
			quality = mwords[10];
			index = mwords[13];
			if ((cigar.length() <= 1) || (quality.length() <= 1) )
				regectread = 1;

			if (regectread == 0) 
			{
				prvcig = cigar;
				fragment = mwords[9];
				fragment = fragment.replace("\"", "");

				frsize = fragment.length();
				quality = mwords[10];
				// base quality field contains characters that are quality
				// indicators and always are important (thus this field must not
				// be cleaned
				// quality = quality.replace("\"", "");
				strand = mwords[11];
				// System.out.println("File Name : " + filenm +
				// "\nChromosome : chr" + chrom + "\nPosition : " + position +
				// "\nSeqSize : " + frsize + "\nSite : " + site +
				// "\nSite Reference : " + site_ref + "\nExpected Value : " +
				// expval + "\nnCigar : " + nCigar + "\nCigar : " + cigar);
				// System.out.println(cigar);
				// for(int w=0; w<mwords.length; w++)
				// System.out.println(mwords[w]);
				rescireadqual = HandleCigSs(cigar, fragment, quality);
				String[] cfq = rescireadqual.split("\n");
				curcig = cfq[0];
				cigar = cfq[0];
				gic1 = cfq[0];
				fragment = cfq[1];
				rd1 = cfq[1];
				quality = cfq[2];
				ciginfo = retrieveCigarInfo(cigar); // retrieve a string that
													// contains digestible cigar
													// info
				// example : "123 M 1 I 25 M 2 D 2 M" (elements separated by
				// space; odds are letters; evens are numbers.)
				String[] wcs = ciginfo.split(" "); // retrieve words with
													// letters and numbers
				dels = new Vector<String>();								// separately
				delsenspoints = new Vector<Integer>();
				// System.out.println(ciginfo);
				// System.out.println("Chromosome : chr" + chrom);
				cnv = new Vector<Integer>(); // # matches or # mismatches or #
				// Insertions or # deletions all
				// together numeric vector
				clv = new Vector<String>(); // the corresponding to
				// previous(numeric) letter vector
				// separate and put the elements of cigar space separated string
				// into the vector
				for (int i1 = 0; i1 < wcs.length; i1++)
				{
					if ((i1 % 2 == 0) && (!wcs[i1].equals("")))
					{
						// System.out.println(wcs[i1]);
						cnv.add(Integer.parseInt(wcs[i1].replace("\"", "")));
					}
					else
					{
						if (!wcs[i1].equals(""))
							clv.add(wcs[i1]);
					}
				}

				// find the actual size
				asize = 0;
				// System.out.println(cnv.size());
				// for loop that calculates the length of the read that
				// finally will be created by adding the deletion points too
				// This length might be greater than the initial length if
				// insertion weren't found but deletions exist.
				// Example 1 Cigar = 20M2D39M, initial length = 59, final
				// length=61 (final cigar = 20M2M39M=61M)
				// Example 2 Cigar = 20M1D1I39M, initial length = 60, final
				// length = 60
				// (one insertion was removed and one D have been added)
				// (final cigar = 20M1D39M=60M since information about Ds
				// included in the read now)

//				Here the Ns must be added
				for (int i1 = 0; i1 < cnv.size(); i1++) 
				{
					if ((clv.get(i1).equals("M")) 
					|| (clv.get(i1).equals("D")) 
					|| (clv.get(i1).equals("N")))
						asize = asize + cnv.get(i1);
				}

				String[] fragpl = fragment.split(""); // fragment string
														// character array
				String[] qualpl = quality.split(""); // quality string character
														// array

				ffrag = new Vector<String>(); // final transformed fragment
				ffrag.addAll(Arrays.asList(fragpl));
				ffrag.remove(0);

				fqual = new Vector<String>(); // final transformed fragment
				fqual.addAll(Arrays.asList(qualpl));
				fqual.remove(0);

				String delsstr = "", ensstr="", inssstr = "";
				delpoints = new Vector<Integer>();
				
				DNs = 0; // Ds counter
				inspoints = new Vector<Integer>();
				int fsize = 0; // fragment size (number of total bases Ms and Ds
				int dnis = 0;
				String delstr = "";
				String enstr = "";
				//enspoints = new Vector<Integer>();
				
				inss = new Vector<String>();
				insstr = ""; insbs = "";
				currpos = 0; // Variable used to deal with the actual read
									// string char array]
									// it is used to target the info about
									// insertions and retrieve the involved
									// bases
									// it is incremented only with presence of
									// base. Deletion cannot increment that
									// variable.

				for (int i2 = 0; i2 < clv.size(); i2++)
				{
					if (clv.get(i2).equals("M"))
					{
						fsize = fsize + cnv.get(i2);
						currpos = currpos + cnv.get(i2); // this is very useful
															// for Insertions
						// currpos is an indicator of before insertion position
						// that depends on the
						// initial size of the read. the final (and mapping
						// position are retrieved form 'fize' variable
					}
					if ( (clv.get(i2).equals("D")) ||  (clv.get(i2).equals("N")) )// if cigar indicates D then
													// retrieve its number and
					// go through a loop in it
					{

						if((clv.get(i2).equals("D")) || clv.get(i2).equals("N"))
						{
							delenstr = (position + fsize - 1) + "_" + cnv.get(i2);
							if(clv.get(i2).equals("D"))
								dels.add(delenstr);
							delsEns.add(delenstr);
						}
				
						for (dnis = 1; dnis <= cnv.get(i2); dnis++) 
						{
							//delpoints.add(fsize + dnis - 1 + Ds);
							delsenspoints.add(fsize + dnis - 1 + DNs);
							if((clv.get(i2).equals("D")))
							{
								//System.out.println("DDD " + (fsize + dnis - 1 + DNs));
								delsensletters.add("D");
								delpoints.add(fsize + dnis - 1 + DNs);
							}
							if((clv.get(i2).equals("N")))
								delsensletters.add("d");
								// System.out.print("\n->" + (fsize+dis-1 + Ds)
								// +"<-\n");
						}
						DNs = DNs + dnis - 1;
						fsize = fsize + cnv.get(i2);

						// fsize=fsize+cnv.get(i);
						// System.out.println("Deletion point will not count in fragment size for now.");
						// System.out.println(delpoints);
					}
					// System.out.println(fsize);
					
					if (clv.get(i2).equals("I")) 
					{
						// for(int ins=1; ins<=cnv.get(i2); ins++)
						// inspoints.add(fsize+ins-1);
						insbs = "";
						for (int j = currpos; j < currpos + cnv.get(i2); j++)
						{
							insbs = insbs + ffrag.get(j);
							inspoints.add(j);
							// ffrag.remove(j);
							// j=j-1;
						}
						insstr = (position + fsize - 1) + "_" + cnv.get(i2)
								+ "_" + insbs;
						inss.add(insstr);
						// fsize=fsize+cnv.get(i2);
						currpos = currpos + cnv.get(i2);
					}
				}
				
//				put insertions(Is) on read using the insertion stored points(indeces)
//				and then remove those Is and their mirrored qualities using the inserted Ises' indeces
				for (int i3 = 0; i3 < inspoints.size(); i3++)
				{
					ffrag.set(inspoints.get(i3), "I");
					// System.out.print(inspoints.get(i3)+ " ");
				}

				int ICounter = 0;
				// System.out.println();

				for (int i3 = 0; i3 < ffrag.size(); i3++)
				{
					if (ffrag.get(i3).equals("I"))
					{
						ffrag.remove(i3);
						fqual.remove(i3);
						// System.out.println("'I' detected!  : " + i3);
						if (i3 > 1)
							i3 = i3 - 1; // Vectors change dynamically their size
						// in 'element removal' the next element will acquire
						// the current index
						// the current index must be then decreased in order for
						// that element to
						// be able to be checked with the next iteration.
						ICounter = ICounter + 1;
					}
				}
//				Now quality and frag strings have no insertions on them 
				
				//
				// System.out.println(ffrag); System.out.println(fqual);
				// System.out.println(delpoints); System.out.println(fragment);
				// System.out.println(rd); System.out.println(rd1);
				// System.out.println(gic); System.out.println(gic1);
				//

//				
//				for (int i3 = 0; i3 < delpoints.size(); i3++)
//				{
//					// ffrag.add(delpoints.get(i), "D" + i);
//					ffrag.add(delpoints.get(i3) - i3, "D");
//					fqual.add(delpoints.get(i3) - i3, "D");
//				}
//
				// System.out.println(delsenspoints.size());
				for (int i3 = 0; i3 < delsenspoints.size(); i3++)
				{
					// ffrag.add(delpoints.get(i), "D" + i);
					if(delsensletters.get(i3).equals("D"))
					{
						//System.out.println(ffrag.size() + "  " + (delsenspoints.get(i3) -i3));
						//System.out.println(cigar);
						ffrag.add(delsenspoints.get(i3) - i3, "D");
						fqual.add(delsenspoints.get(i3) - i3, "D");
					}
					else if(delsensletters.get(i3).equals("d"))
					{
						ffrag.add(delsenspoints.get(i3) - i3, "d");
						fqual.add(delsenspoints.get(i3) - i3, "d");
					}
				}
				// System.out.println(ffrag);
				// System.out.println(fqual);
		
		
				fifrag = "";
				for (int i3 = 0; i3 < ffrag.size(); i3++)
					fifrag = fifrag + ffrag.get(i3);

				fiqual = "";			
				for (int i3 = 0; i3 < fqual.size(); i3++)
					fiqual = fiqual + fqual.get(i3);

				delsstr = "";
				// returned read name will be cmn = common if no deletions and
				// no insertions have been detected.
				readnmreturn = "cmn";
				if (dels.size() > 0)
				{
					for (int in = 0; in < dels.size(); in++)
					{
						if (in < dels.size() - 1)
							delsstr = delsstr + dels.get(in) + ",";
						if (in == dels.size() - 1)
							delsstr = delsstr + dels.get(in);
					}
					// returned read name will be cmn = common if no deletions
					// and no insertions have been detected.
					readnmreturn = readname;
				}
				if (dels.size() == 0)
					delsstr = "none";

				inssstr = "";
				if (inss.size() > 0)
				{
					for (int in = 0; in < inss.size(); in++)
					{
						if (in < inss.size() - 1)
							inssstr = inssstr + inss.get(in) + ",";
						if (in == inss.size() - 1)
							inssstr = inssstr + inss.get(in);
					}
					// returned read name will be cmn = common if no deletions
					// and no insertions have been detected.
					readnmreturn = readname;
				}
				if (inss.size() == 0)
					inssstr = "none";
				// 03/06/2015 Changed to accept X and Y chromosome (chrom to
				// chromS)
				secres.add(chromS + "\t"
						+ (position)
						+ "\t"
						// 03/06/2015 Changed to accept X and Y chromosome
						// (chrom to chromS)
						+ (position + fifrag.length() - 1) + "\t" + fifrag
						+ "\t" + fiqual + "\t" + delsstr + "\t" + inssstr
						+ "\t" + cigar + "\t" + curcig + "\t" + prvcig + "\t"
						+ readnmreturn + "\t" + filenm + "\t" +index);
			}
			// Position already includes the first place of length therefore to
			// find read end actual position add length and subtract one
		}

		System.out.println("................Preparation  Ended.");
		System.out.println("	End of 'PrepareTheData' Method!");
		return secres;
	}	
	
	// ---------------------------------------------------------------------------------------------------\\
	// --START--------------------------------------------------------------------------------------------\\
	// -------------------------------TASK 4 -- DECIDE
	// ---------------------------------------------------\\
	// DECIDE Process (START)
	/**
	 * // Method that accepts the table path the output path and the coverage //
	 * It reads the table identifies and stores the sites of interest // Next it
	 * calculates for each site the range of its frequency towards the // range
	 * of the data frequency // The method became significant faster than the
	 * original. // In short : it populates column parallel arrays from the data
	 * table. // For each site it uses the particular set of these column arrays
	 * to calculate the // data frequency through 'retriveDataFreq' method (the
	 * arrays are sent to that method // as parameters, along with reference and
	 * alternative bases and the coverage limit with germlines exclusion
	 * percentage(prf)). // // 06/18/2015 // Updated : Decide (and
	 * 'retriveDataFreq' method) on 06/18/2015 The method is really fast // The
	 * method creates six combined mirrored Double data(table) frequency
	 * arrays(Vectors) // Next for each site of interest calculates the p value
	 * by positioning sites frequency // in the corresponding previously created
	 * Double array. // End 06/18/2015
	 * 
	 * 
	 * @param intab
	 *            : String input table
	 * @param outres
	 *            : String output file
	 * @param coverage
	 *            : String lower limit for # reads found that contain a valid
	 *            position for this data line
	 * @param germlineAFs
	 *            : String lower limit for germline
	 * @throws IOException
	 *             : Open and read file is involved
	 **/
	public static void Decide(String intab, String outres, String coverage,
			String germlineAFs, Double pvd) throws IOException
	// , InterruptedException
	{
		/*
		 * Somatic   : p-value < 0.05 && AF <= 0.35 
		 * Germline  : AF > 0.35
		 * Omitted   : Coverage < 100 or user entered value 
		 * Undefined : the rest 
         * prf = 0.0 to 0.95 
         * (	xv positional distance from left, 
         * 		yv positional distance from	right, 
         * 		prf lower bound for (number of base expected / coverage) 
		 * 		dv coverage   )
		 */
		int covr = 0;
		double germlineAF = 0.0;
		germlineAF = Double.parseDouble(germlineAFs);

		if (pvd<=0.0)
			pvd=0.05;

		covr = Integer.parseInt(coverage);
		System.out.println("Start Decide method!");
		System.out.println("Table in : " + intab);
		System.out.println("Results out : " + outres);
		// a) Retrieve table data and store each line to a String Vector
		File tbfl = new File(intab);
		String flnm = tbfl.getName();
		Vector<String> tbdt = new Vector<String>(); // vector for
		tbdt = readTheFile(intab, flnm);
		int tbsz = tbdt.size();
		
		// data table size to bee used in loops
		// Make it faster by populating the table String column arrays only one
		// time
		// Declare Variables

		System.out.println("Declare and populate the data arrays!");
		Vector<Double> sbst = new Vector<Double>();
		String[] chromosome = new String[tbdt.size()], site = new String[tbdt
				.size()], x = new String[tbdt.size()], y = new String[tbdt
				.size()], base_expected = new String[tbdt.size()], alternative = new String[tbdt
				.size()], number_base_expected = new String[tbdt.size()], number_of_As = new String[tbdt
				.size()], number_of_Cs = new String[tbdt.size()], number_of_Ts = new String[tbdt
				.size()], number_of_Gs = new String[tbdt.size()], number_of_Ds = new String[tbdt
				.size()], number_of_Ps = new String[tbdt.size()];

		// Some future changes maybe will change the sequence of the columns
		// Please test the output if the order is not correct.

		// Access row element of a particular column by index(same index = same
		// row = same position
		// chromo site x y base altern number number number number number
		// some expec ative base _of_As _of_Cs _of_Ts _of_Gs
		// ted expec
		// index ted

		// ===============================================================================================\\
		// [0] 1 2074785 1 101 A A 546 545 0 1 0
		// [1] 1 2074786 2 100 C C 550 0 550 0 0
		// [2] 1 2074787 3 99 T T 551 0 0 551 0
		// [3] 1 2074788 4 98 C C 557 0 556 0 1
		// [4] 1 2074789 5 97 A A 566 565 0 1 0
		// [5] 1 2074790 6 96 G G 571 0 0 0 571
		// [6] 1 2074791 7 95 G G 577 1 0 0 576
		// [7] 1 2074792 8 94 T T 583 1 0 580 2
		// [8] 1 2074793 9 93 G G 590 0 0 0 590
		// ===============================================================================================\\

		double prf = 0;
		String chrm = "", be = "", ref = "", alt = "";
		int nbe = 0, dv = 0, tC = 0, tA = 0, tRef = 0, tMut = 0;
		int xv = 0, yv = 0;
		double qt = 0.0;
		// String site="";
		String line = "";
		Vector<String> stofint = new Vector<String>(); // Vector for sites of interest
		// populate the parallel arrays using the table
		// and Retrieve from the data table Vector all table site lines
		String mutline="";
		double sitefreq = 0.0;
		int tref = 0, talt = 0, stdv = 0;

		for (int i = 0; i < tbdt.size(); i++) 
		{
			line = tbdt.get(i);
			String[] wrs = line.split("\t");
			chromosome[i] = wrs[0];
			site[i] = wrs[1];
			x[i] = wrs[2];
			y[i] = wrs[3];
			ref  = wrs[4];
			alt  = wrs[5];
			base_expected[i] = ref;
			alternative[i] = alt;
			number_base_expected[i] = wrs[6];
			stdv = Integer.parseInt(wrs[6]);
			number_of_As[i] = wrs[7];
			number_of_Cs[i] = wrs[8];
			number_of_Ts[i] = wrs[9];
			number_of_Gs[i] = wrs[10];
			number_of_Ds[i] = wrs[11];
			number_of_Ps[i] = wrs[12];
			// if( (!(ref.equals(alt))) && (stdv>=covr) )
//		11/03/15 adapted to handle more alternative options for same position
//		if alternative field contains two or more letters(Nucleotides), then
//		the next enclosed if and for loop statments will create separate
//		site line for each alternative Nucleotide for this particular position.
			if (!ref.equals(alt))
			{
				if(alt.length()>1)
				{
					for(int j=0; j<alt.length(); j++)
					{
						mutline=wrs[0] + "\t" + wrs[1] + 
						"\t" + wrs[2] + "\t" + wrs[3] + 
						"\t" + wrs[4] + "\t" + alt.charAt(j) + 
						"\t" + wrs[6] + "\t" + wrs[7] + 
						"\t" + wrs[8] + "\t" + wrs[9] + 
						"\t" + wrs[10] + "\t" + wrs[11] + 
						"\t" + wrs[12]; 
						stofint.add(mutline);
					}
				}
				else
					stofint.add(line);
			}
		}

//		11/03/15 adapted to handle more alternative options for same position
//		if alternative field contains two or more letters(Nucleotides), then
//		the previous enclosed if and for loop statment will create separate
//		site line for each alternative Nucleotide for this particular position.

		// End Make it faster by populating the table String column arrays only
		// one time

		// coverage for sites of interest
		System.out.println("Sites size = " + stofint.size());

		String decide = "", lineres = "";

		// Combined
		String STATEMENT = "";
		// String Header = "CHR\tPOS\tREF\tEXP\tCOV\tFREQ(#ALT/COV)\tP-VALUE\tCONCLUSION";
		// (somatic|germline|unknown|omitted)";

		Vector<String> results = new Vector<String>();

		Vector<Double> AtoCandTtoG = new Vector<Double>();
		Vector<Double> AtoGandTtoC = new Vector<Double>();
		Vector<Double> AtoTandTtoA = new Vector<Double>();
		Vector<Double> CtoTandGtoA = new Vector<Double>();
		Vector<Double> CtoAandGtoT = new Vector<Double>();
		Vector<Double> CtoGandGtoC = new Vector<Double>();
		Vector<Double> mirror = new Vector<Double>();

		// int tbsz, int cvri, double prf, String ref, String mu
		AtoCandTtoG = retrieveDataFreq(tbsz, covr, prf, "A", "C", chromosome,
				site, x, y, base_expected, alternative, number_base_expected,
				number_of_As, number_of_Cs, number_of_Ts, number_of_Gs);
		mirror = retrieveDataFreq(tbsz, covr, prf, "T", "G", chromosome, site,
				x, y, base_expected, alternative, number_base_expected,
				number_of_As, number_of_Cs, number_of_Ts, number_of_Gs);
		for (int p = 0; p < mirror.size(); p++)
			AtoCandTtoG.add(mirror.get(p));
		AtoCandTtoG = sortVector(AtoCandTtoG);

		AtoGandTtoC = retrieveDataFreq(tbsz, covr, prf, "A", "G", chromosome,
				site, x, y, base_expected, alternative, number_base_expected,
				number_of_As, number_of_Cs, number_of_Ts, number_of_Gs);
		mirror = retrieveDataFreq(tbsz, covr, prf, "T", "C", chromosome, site,
				x, y, base_expected, alternative, number_base_expected,
				number_of_As, number_of_Cs, number_of_Ts, number_of_Gs);
		for (int p = 0; p < mirror.size(); p++)
			AtoGandTtoC.add(mirror.get(p));
		AtoGandTtoC = sortVector(AtoGandTtoC);

		AtoTandTtoA = retrieveDataFreq(tbsz, covr, prf, "A", "T", chromosome,
				site, x, y, base_expected, alternative, number_base_expected,
				number_of_As, number_of_Cs, number_of_Ts, number_of_Gs);
		mirror = retrieveDataFreq(tbsz, covr, prf, "T", "A", chromosome, site,
				x, y, base_expected, alternative, number_base_expected,
				number_of_As, number_of_Cs, number_of_Ts, number_of_Gs);
		for (int p = 0; p < mirror.size(); p++)
			AtoTandTtoA.add(mirror.get(p));
		AtoTandTtoA = sortVector(AtoTandTtoA);

		CtoTandGtoA = retrieveDataFreq(tbsz, covr, prf, "C", "T", chromosome,
				site, x, y, base_expected, alternative, number_base_expected,
				number_of_As, number_of_Cs, number_of_Ts, number_of_Gs);
		mirror = retrieveDataFreq(tbsz, covr, prf, "G", "A", chromosome, site,
				x, y, base_expected, alternative, number_base_expected,
				number_of_As, number_of_Cs, number_of_Ts, number_of_Gs);
		for (int p = 0; p < mirror.size(); p++)
			CtoTandGtoA.add(mirror.get(p));
		CtoTandGtoA = sortVector(CtoTandGtoA);

		CtoAandGtoT = retrieveDataFreq(tbsz, covr, prf, "C", "A", chromosome,
				site, x, y, base_expected, alternative, number_base_expected,
				number_of_As, number_of_Cs, number_of_Ts, number_of_Gs);
		mirror = retrieveDataFreq(tbsz, covr, prf, "G", "T", chromosome, site,
				x, y, base_expected, alternative, number_base_expected,
				number_of_As, number_of_Cs, number_of_Ts, number_of_Gs);
		for (int p = 0; p < mirror.size(); p++)
			CtoAandGtoT.add(mirror.get(p));
		CtoAandGtoT = sortVector(CtoAandGtoT);

		CtoGandGtoC = retrieveDataFreq(tbsz, covr, prf, "C", "G", chromosome,
				site, x, y, base_expected, alternative, number_base_expected,
				number_of_As, number_of_Cs, number_of_Ts, number_of_Gs);
		mirror = retrieveDataFreq(tbsz, covr, prf, "G", "C", chromosome, site,
				x, y, base_expected, alternative, number_base_expected,
				number_of_As, number_of_Cs, number_of_Ts, number_of_Gs);
		for (int p = 0; p < mirror.size(); p++)
			CtoGandGtoC.add(mirror.get(p));
		CtoGandGtoC = sortVector(CtoGandGtoC);

		System.out.println("Combined");

		// c) 1. Go through each site (String line element of 'stofint' Vector)
		// and calculate its frequency :
		// using its reference and alternative field (calculate: alt found /
		// total) and proceed in 2.
		// 2. Go to the corresponding data frequency vector and through a for
		// loop Statement find the index
		// where the value of data is greater than or equal to that of the
		// frequency of the site.
		// 4. calculate p-value by p = (size-index) / size
		// size = curfrq.size()
		// index = as stated in the previous sentence
		// chr pos id(.) REF ALT COV FREQ(#ALT/COV) P-VALUE
		// STATEMENT(somatic|germline|unknown)

		
		for (int i = 0; i < stofint.size(); i++) // 1.
		// for(int i=0; i<20; i++) // 1
		{
			line = stofint.get(i);
			String[] wrds = line.split("\t");
			ref = wrds[4];
			alt = wrds[5];
			if (ref.equals("A"))
				tref = Integer.parseInt(wrds[7]);
			if (ref.equals("C"))
				tref = Integer.parseInt(wrds[8]);
			if (ref.equals("T"))
				tref = Integer.parseInt(wrds[9]);
			if (ref.equals("G"))
				tref = Integer.parseInt(wrds[10]);
			if (alt.equals("A"))
				talt = Integer.parseInt(wrds[7]);
			if (alt.equals("C"))
				talt = Integer.parseInt(wrds[8]);
			if (alt.equals("T"))
				talt = Integer.parseInt(wrds[9]);
			if (alt.equals("G"))
				talt = Integer.parseInt(wrds[10]);
			dv = Integer.parseInt(wrds[6]);
			if (dv > 0)
				sitefreq = (double) talt / (double) dv; // it could be zero
			else
				sitefreq = 0.0;

			Vector<Double> curfrq = new Vector<Double>();
			if ((ref.equals("A") && alt.equals("C"))
					|| (ref.equals("T") && alt.equals("G")))
				curfrq = AtoCandTtoG;
			else if ((ref.equals("A") && alt.equals("G"))
					|| (ref.equals("T") && alt.equals("C")))
				curfrq = AtoGandTtoC;
			else if ((ref.equals("A") && alt.equals("T"))
					|| (ref.equals("T") && alt.equals("A")))
				curfrq = AtoTandTtoA;
			else if ((ref.equals("C") && alt.equals("T"))
					|| (ref.equals("G") && alt.equals("A")))
				curfrq = CtoTandGtoA;
			else if ((ref.equals("C") && alt.equals("A"))
					|| (ref.equals("G") && alt.equals("T")))
				curfrq = CtoAandGtoT;
			else if ((ref.equals("C") && alt.equals("G"))
					|| (ref.equals("G") && alt.equals("C")))
				curfrq = CtoGandGtoC;

			double pval = 1.0;
			int cur_size = curfrq.size();
			for (int o = 0; o < cur_size; o++) {
				pval = 1 - (double) o / (double) cur_size;
				if (curfrq.get(o) >= sitefreq)
					break;
				if (o == cur_size - 1)
					pval = 0.0;
			}

			// germlineAF = 0.35
			if (dv < covr)
				STATEMENT = "omitted";
			else if (sitefreq > germlineAF)
				STATEMENT = "germline";
			//else if (pval < 0.05)
			else if (pval < pvd)
				STATEMENT = "somatic";
			else
				STATEMENT = "undefined";
			//String Header = "CHR\tPOS\tREF\tEXP\tCOV\tFREQ(#EXP/COV)\tP-VALUE\tCONCLUSION";	
			lineres = wrds[0] + "\t" + // chromosome
					wrds[1] + "\t" + // site
					// "." + "\t" +
					// wrds[2] + "\t" + // x
					// wrds[3] + "\t" + // y
					wrds[4] + "\t" + // base_expected (same as reference) 
					wrds[5] + "\t" + // alternative (not same as reference) now it is expected (we expect it as somatic) 
					wrds[6] + "\t" + // number_base_expected 
					// wrds[7] + "\t" + // number_of_As
					// wrds[8] + "\t" + // number_of_Cs
					// wrds[9] + "\t" + // number_of_Ts
					// wrds[10]; // number_of_Gs
					// + "\t" +
					// wrds[11] + "\t" + // number_of_Ds
					// wrds[12]; // number_of_Ps
					+sitefreq + "\t" + pval + "\t" + STATEMENT;

			// /(double)curfrq.size();
			// System.out.println(decide + "\t" + pval + "\t" + frac + "\t" +
			// lineres);
			// if(i%20==0)
			// System.out.println(i);
			// results.add(decide + "\t" + pval + "\t" + frac + "\t" + lineres);
			results.add(lineres);
		}
		
		String Header = "CHR\tPOS\tREF\tEXP\tCOV\tFREQ(#EXP/COV)\tP-VALUE\tCONCLUSION";	
//		find minimum pvalue greater than zero
		Vector<String>finalresults = new Vector<String>();
		finalresults.add(Header);
		String chromosomes="", sites="", nmbexps="", 
				sitefreqs="", base_expecteds="", 
				alternatives="", pvals="", statements="";
		int cov=0;
		
		double pvl=0.0, pvmin=1.0; // pvalue, minimum pvalue
		double pvrplacezero=0.0;  // minimum pvalue greater than zero divided by 2
		
		for(int i=0; i<results.size(); i++)
		{
			line=results.get(i);
			String[] wrds = line.split("\t");
			pvals =  wrds[6];
			pvl=Double.parseDouble(pvals);
			if( (pvl>0) && (pvl<pvmin) )
				pvmin=pvl;
		}
		pvrplacezero=pvmin/2;
		
		String firesln="";
		for(int i=0; i<results.size(); i++)
		{
			firesln="";
			line=results.get(i);
			String[] wrds = line.split("\t");
			chromosomes = wrds[0];		//	chromosome
			sites = wrds[1];		//	site
			base_expecteds = wrds[2];	//	base_expected
			alternatives = wrds[3];		//	alternative
			nmbexps = wrds[4];		//	number_base_expected
			cov = Integer.parseInt(nmbexps);
			sitefreqs = wrds[5];
			pvals =  wrds[6];
			pvl=Double.parseDouble(pvals);
			statements = wrds[7];
			
			//if(pvl==0.0)
			//{
			//	pvl=pvrplacezero;
			//	pvals=Double.toString(pvl);
			//	statements="SensSomatic";
			//}

			if(pvl==0.0)
			{
				pvl=pvrplacezero;
				pvals=Double.toString(pvl);
				if(cov>=covr)
					statements="Sens_"+ statements;
			}
			firesln=chromosomes  + "\t" + 
					sites + "\t" + 
					base_expecteds + "\t" + 
					alternatives + "\t" + 
					nmbexps + "\t" + 
					sitefreqs + "\t" + 
					pvals + "\t" + 
					statements;
			finalresults.add(firesln);
		}
		writeToFile(outres, finalresults);
		System.out.println("End   Decide method.");
	}

	
	
	
	

	// Method that accepts the table path the output path and the coverage
	// It reads the table identifies and stores the sites of interest
	// Next it calculates for each site the range of its frequency towards the
	// range of the data frequency
	// The method became significant faster than the original.
	// In short : it populates column parallel arrays from the data table.
	// For each site it uses the particular set of these column arrays to
	// calculate the
	// data frequency through 'retriveDataFreq' method (the arrays are sent to
	// that method
	// as parameters, along with reference and alternative bases and the
	// coverage limit with germlines exclusion percentage(prf)).

	// The following method generates a positional vector (# specific quotients
	// = mutations/coverage) after
	// filtering the passed table. It also accepts two more parameters: the
	// reference and the
	// alternative base. The metghod derived from a conversion of the original R
	// script
	// which is used to generate pvalue graphs. It uses many parallel Arrays to
	// store values for each table row word
	// then when it retrieves those values it uses the index that corresponds to
	// a particular column.
	// Variables that could be set up:
	// xv positional distance from left
	// yv positional distance from right,
	// prf lower bound of (number of base expected / coverage)
	// dv coverage

	/**
	 * 'retrieveDataFreq' method accepts eleven columns of data table as String
	 * arrays. Data columns populated the arrays(in 'Decide' method) that are
	 * passed as parameters here. elements of different arrays with same index
	 * belong to the same table row. Each data row holds information for a
	 * particular position. The descriptive names(table column names) of the
	 * arrays are : chromosome, site, x, y, base_expected, alternative,
	 * number_base_expected, number_of_As, number_of_Cs, number_of_Ts, and
	 * number_of_Gs. The method generates a data array of frequencies of all
	 * data positions that have expected base = 'ref'. The frequency for each
	 * included position is given by the formula:
	 * number_of_(mut)/number_base_expected. Both 'ref' and 'mut' are passed as
	 * parameters in the method.
	 * 
	 * Main Idea of returned data array: include all positions that have: 1)
	 * expected base = ref, 2) number of expected bases>minimum coverage(cvri), 3)
	 * (number of mut found/number of expected bases)<prf then insert to the
	 * returned array(Vector) the frequency = (number of mut found/number of
	 * expected baases) Parameters
	 * 
	 * @param tbsz
	 *            : integer table total rows size
	 * @param cvri
	 *            : integer minimum position coverage
	 * @param prf
	 *            : double that represents the lower frequency border of data
	 *            inclusion example : prf <
	 *            number_base_expected_found/number_base_expected If for example
	 *            prf = 0.95 then any position with germline characteristics is
	 *            excluded.
	 * @param ref
	 *            : String reference base
	 * @param mut
	 *            : String
	 * 
	 * @param chromosome
	 * @param site
	 * @param x
	 * @param y
	 * @param base_expected
	 * @param alternative
	 * @param number_base_expected
	 * @param number_of_As
	 * @param number_of_Cs
	 * @param number_of_Ts
	 * @param number_of_Gs
	 * @return
	 * 
	 */
	public static Vector<Double> retrieveDataFreq(int tbsz, int cvri,
			double prf, String ref, String mut, String[] chromosome,
			String[] site, String[] x, String[] y, String[] base_expected,
			String[] alternative, String[] number_base_expected,
			String[] number_of_As, String[] number_of_Cs,
			String[] number_of_Ts, String[] number_of_Gs) {
		// System.out.println("Method Start");
		// Retrieve the total number of (Ref) Base Expected and store it in a
		// index driven Array
		// Retrieve also the total number of the mutation in question(
		// The Two arrays derived from the following loop will be filtered and
		// used to
		//
		Vector<Double> sbst = new Vector<Double>();
		int xv = 0, yv = 0, dv = 0, tRef = 0, tMut = 0;
		String be = "", alt = "";
		double qt = 0.0;
		for (int i = 0; i < tbsz; i++) {
			// System.out.println("--- " + i + "    " + tbdt.size());
			// chromosome[i];
			// site[i];
			xv = Integer.parseInt(x[i]);
			yv = Integer.parseInt(y[i]);
			be = base_expected[i];
			alt = alternative[i];
			if (ref.equals("A"))
				tRef = Integer.parseInt(number_of_As[i]);
			if (ref.equals("C"))
				tRef = Integer.parseInt(number_of_Cs[i]);
			if (ref.equals("T"))
				tRef = Integer.parseInt(number_of_Ts[i]);
			if (ref.equals("G"))
				tRef = Integer.parseInt(number_of_Gs[i]);
			if (mut.equals("A"))
				tMut = Integer.parseInt(number_of_As[i]);
			if (mut.equals("C"))
				tMut = Integer.parseInt(number_of_Cs[i]);
			if (mut.equals("T"))
				tMut = Integer.parseInt(number_of_Ts[i]);
			if (mut.equals("G"))
				tMut = Integer.parseInt(number_of_Gs[i]);
			dv = Integer.parseInt(number_base_expected[i]);
			// dv=nbe;
			// add tC/dv if :
			// 1 (be=='A' & alt=='A')
			// 2 &

			// 3 & xv>xval
			// 4 & yv>yval
			// 5 & tC/dv < 0.005
			// 6 & dv>100))
			if ((xv > 0) && (yv > 0)) // 3, 4
			{
				if (dv >= cvri)// 6 (dv =coverage) include only positions with
								// coverage > 1, 100
				{
					if ((double) tRef / (double) dv > prf) // 2
					{
// System.out.println("Here1");
						if ((be.equals(ref)) && (alt.equals(ref))) 
// 1 this restricts that the sites of interest are out(in sites ref not equals to alt)
						{
							// System.out.println("Here1,2    " + tC + "  " + i
							// + "   " + number_of_Cs[i]);
							// if(tMut!=0)
							{
								// System.out.println("Here2");
								qt = (double) tMut / (double) dv; 
// 0 divide the number of total specific mutations found in a position
								{ // by the position coverage
									if ((double) tMut / (double) dv < 0.005) 
// 5 include only the range of cut off
									{
										sbst.add(qt); // Then add the qt in the
														// the vector
									}
								}
							}
						}
					}
				}
			}

			// Integer.parseInt(number_of_Gs[i]);
			// Integer.parseInt(number_of_Ds[i]);
			// Integer.parseInt(number_of_Ps[i]);
			// tC/dv, (be=='A' & alt=='A' & tA/dv>prf & xv>xval & yv>yval &
			// tC/dv < 0.005 & dv>100))
		}
		sbst = sortVector(sbst);
		// System.out.println("Method End");

		return sbst;
	}

	/**
	 * 'sortVector' method sorts in ascending order a vector of Doubles it
	 * converts the Vector<Double> to a Double[] Array; it uses Arrays.sort to
	 * sort the array; it then converts the array to vector and returns it.
	 * 
	 * @param vec:
	 *            is a Double Vector(unsorted
	 * @return res: Double Vector sorted
	 */
	public static Vector<Double> sortVector(Vector<Double> vec) {
		Vector<Double> res = new Vector<Double>();
		Double[] arvc = new Double[vec.size()];
		vec.toArray(arvc);
		Arrays.sort(arvc);
		res = new Vector<Double>(Arrays.asList(arvc));
		return res;
	}

	// DECIDE Process (END)
	// -------------------------------TASK 4 -- DECIDE
	// ---------------------------------------------------\\
	// ------------------------------------------------------------------------------END------------------\\
	// ---------------------------------------------------------------------------------------------------\\

	
	// ----START-----------------------------------------------------------------\\
	// --------------- TEST ARGUMENTS--------------------------------------------\\
	// --------------------------------------------------------------------------\\	
//		EXTRACT EXTRACT (START)--------------------------------------------------------------------------\\	
	// -------------------------------------------------------------------------------------------\\

		public static void Extract(String bash_pth, String samtools_path, String outextract, 
				String positionextract, String alellextract,
				Vector<String> bms, String tmppath) throws IOException, InterruptedException
		{
			char alt=alellextract.charAt(0);
			
			int refstart=0, refend=0;
			String site = "";
			String[] wds = positionextract.split("-");
			if(wds.length==2)
				site = wds[1];
			refstart = Integer.parseInt(site)-100;
			refend = Integer.parseInt(site)+100;
			String[] wdch = positionextract.split(":");
			String chrp="";
			chrp=wdch[0];
			//chrp=wdch[0].substring(3,wdch[0].length());
			String locus= chrp + ":" + refstart + "-" + refend;
			System.out.println(locus);
			// positionextract = locus
			System.out.println("Start Extract Method!");
			System.out.println(bash_pth + "\n" + samtools_path + 
						"\n" + outextract + "\n" + locus + 
						"\n" + alellextract + "\n" +  tmppath);
			for(int i=0; i<bms.size(); i++)
				System.out.println(bms.get(i));
			
			
			Vector<String> crdvct = new Vector<String>();
			Vector<String> totaldatacol = new Vector<String>();
			Vector<String> cursmlq = new Vector<String>();
			String tmpuiid = "";
			tmpuiid = tmppath + "/"
					+ UUID.randomUUID().toString().replaceAll("-", "i");
			File dutmp = new File(tmpuiid);
			if (!dutmp.exists())
				dutmp.mkdir();
			String cursmline="", fnline="";
			String bflnm ="", readname ="", strand ="", chromosome ="", 
					positions ="", MAPQ ="", cigar ="", flg ="",
					sequence ="", qual ="", reference="", alternative="";
			alternative = alellextract;
			int MAPQi=0, mqi=0; 
			
			
			
			totaldatacol = RetrieveDataFromBam(bash_pth, samtools_path,
						tmpuiid, locus, bms);
			//System.out.println("Elaborate Data for primer : " + h);
			System.out.println("Current Vector Size : " + totaldatacol.size());
			Vector<String> testError = new Vector<String>();
			for (int i = 0; i < totaldatacol.size(); i++) 
			{
				cursmline = totaldatacol.get(i);
				String[] cslwords = cursmline.split("\t");
				// System.out.println(i + "<-    ->" + cslwords.length);
				if (cslwords.length > 9)
				{
					bflnm = cslwords[0];
					readname = cslwords[1]; // inserted on 03/16/2015 to be
											// farther elaborated from
											// PrepareTheDataAndInDls method
					strand = cslwords[2];
					chromosome = cslwords[3];
					// System.out.println(chromosome);
					positions = cslwords[4];
					MAPQ = cslwords[5];
					MAPQi = Integer.parseInt(MAPQ);
					cigar = cslwords[6];
					flg = cslwords[10];
					sequence = cslwords[10];
					qual = cslwords[11];
					// System.out.println("MAPQi : " + MAPQi +
					// "  Read Length = " + sequence.length() +
					// "  Qual Length = " + qual.length());
					// truefalse=cslwords[12];
					// fnline="", chromosome="", position="", Ncigar=,
					// cigar="", sequence="", flg="", qual="", truefalse="";
					if (MAPQi > mqi)
					// if(MAPQi>0)
					{
						fnline = bflnm + "\t" + site + "\t" + chromosome + "\t"
								+ reference + "\t" + alternative + "\t" + chromosome
								+ "\t" + positions + "\t" + MAPQ + "\t"
								+ cigar + "\t" + sequence + "\t" + qual
								+ "\t" + strand + "\t" + readname +"\t" + i;
						cursmlq.add(fnline);
						//System.out.println("fnline : " + fnline);
						//if(h==2)
							//testError.add(sequence);
					}
				}

				// if(cslwords.length<9)
				// System.out.println(cslwords[3]);
			}
			
			int totallines = cursmlq.size();
			String redlines="", readseq="", readqual="";
			Vector<String> results = new Vector<String>();
			Vector<String> results2 = new Vector<String>();
			//
			//
			//
			// // for (int k=0; k<totallines; k++)
			// //cursmlq.set(k, cursmlq.get(k)+ "\t" + totallines);
			//

			System.out.println("Size = " + cursmlq.size());
			cursmlq = PrepareTheDataAndInDlsEx(cursmlq);
//			System.out.println("Size after Prepare The data: " + cursmlq.size());
			// 03/16/2015 cursmlq contains prepared reads with additional
			// information
			// each cursmlq's line contains read name at index 10 starting
			// from 0
			// for common reads(no indels) read name = cmn
			// for reads that have insertions or/and deletions the read name
			// is the actual name
			// for (int k=0; k<cursmlq.size(); k++)
			// PrimCumula.add(cursmlq.get(k));
			// writeToFile(tmppath + "/_" + h + "_" + h + "_test007.tsv",
			// cursmlq);
			// writeToFile(tmppath + "/_" + h + "_" + h +
			// "_test00hfgh7.tsv", PrimCumula);
			// PrimCumula = new Vector<String>();
			String rdnm="",
			dataline="",
			datseq="",
			datqual="",
			vDeletions="",
			vInsertions="",
			refidS="",
			smplnmi = "",
			padedseq="",
			index="",
			newline="";
			int refsize=200;
			int position = 0;
			int actualend =0;
			int indexi=0;
			Vector<String> vlns = new Vector<String>();
								// normal.
			
			for (int i2 = 0; i2 < cursmlq.size(); i2++) // go through each
														// read(line)
			{
				rdnm = "";
				dataline = cursmlq.get(i2);
				String[] dawords = dataline.split("\t");
//				refid = Integer.parseInt(dawords[0]);
				refidS = dawords[0];
				position = Integer.parseInt(dawords[1]);
				actualend = Integer.parseInt(dawords[2]);
				datseq = dawords[3];
				// System.out.println(" datseq : " + datseq);
				datqual = dawords[4];
				// System.out.println(" datqual : " + datqual);
				vDeletions = "";
				vInsertions = "";
				cigar = "";
				vDeletions = dawords[5];
				vInsertions = dawords[6];
				cigar = dawords[7];
				rdnm = dawords[10]; // read name "cnm" if the read is
									// normal.
				smplnmi = dawords[11];
				index = dawords[12];
	//Limits
				//if (//(position >= refstart - 100)
						//&& (position <= refstart + 100)
						//&& 
					//	(position + datseq.length() > refstart)
					//	&& (position <= refend)) // position matches
				// System.out.println(dataline);
				{
					// ercnt=ercnt+1;
					padedseq = padDataSequenceExtract(datseq, datqual, position,
							refstart, refsize);
					//System.out.println(padedseq);
					// results.add(padedseq); //
					// LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
					// after modifications are done for quality control the
					// padedseq
					// contains the paded sequence and its quality separated
					// by tab.

					// Distance / No Distance Control
					//if (dcuti != -1) {
						//d1 = Math.abs(refstart - position);
						//d2 = Math.abs(actualend - refend);
						// System.out.println("d1 = " + d1 + "  d2 = " +
						// d2);
					//System.out.println("padedseq.charAt(100) : " + padedseq.charAt(100));
						if ((padedseq.length() > 2))
								//&& ((d1 < dcuti) && (d2 < dcuti))) 
						{
						//	newline = refid + "\t" + refstart + "\t"
							//System.out.println("padedseq.charAt(100) : " + padedseq.charAt(100));
							newline = refidS + "\t" + refstart + "\t"
									+ actualend + "\t" + padedseq + "\t"
									+ vDeletions + "\t" + vInsertions
									+ "\t" + cigar + "\t" + rdnm + "\t"
									+ smplnmi + "\t" + index;
							vlns.add(newline);
							indexi=Integer.parseInt(index);
							//String[] stri = padedseq.split(" ");
							//System.out.println("padedseq.charAt(100) : " + padedseq.charAt(100) + "  " + stri[1]);
							if(padedseq.charAt(100)==alt)
									//alternative,toString())
							{
								
								cursmline = totaldatacol.get(indexi);
								String[] cslwords = cursmline.split("\t");
								// System.out.println(i + "<-    ->" + cslwords.length);
								if (cslwords.length > 9)
								{
									readname = "@"+cslwords[1];
									sequence = cslwords[10];
									qual = cslwords[11];
									redlines=readname+"\n"+sequence+"\n+\n"+qual;
									results.add(redlines);
									//String[] str = padedseq.split("\t");
									//if(str.length==2)
									//{
										//readseq=str[0];
									//readqual=str[1];
									//redlines=rdnm+"\n"+readseq+"\n+\n"+readqual;
									//results2.add(redlines);
								}
							}
							
						}
					}
				/*
				else if (dcuti == -1) {
						if (padedseq.length() > 2) {
							newline = refidS + "\t" + refstart + "\t"
									+ actualend + "\t" + padedseq + "\t"
									+ vDeletions + "\t" + vInsertions
									+ "\t" + cigar + "\t" + rdnm + "\t"
									+ smplnmi;
							vlns.add(newline);
							if(h==2)
							testError.add(padedseq);
						}
					}
				*/

				}

				// Distance / No Distance Contro
			if(outextract.equals("p"))
			{
				for(int i=0; i<results.size(); i++)
					System.out.println(results.get(i));
			}
			else
				writeToFile(outextract, results);
			//writeToFile(outextract+"vlns.tsv", vlns);
			System.out.println("End Extract Method!");
		}
		
		public static String padDataSequenceExtract(String vcseq, String datql,
				int vcposition, int refstart, int refsize) {
			// constant 450 length : substitute refsize with 450 everywhere in this
			// method
			String resultdtseq = "", resultedqual = "", reseqqual = "";
			int diff = 0, size = 0, dfgt = 0, dflt = 0, ngtsize = 0;
			size = vcseq.length();
			resultdtseq = vcseq;
			resultedqual = datql;
			// System.out.println("fragment  position = " + vcposition);
			// System.out.println("reference position = " + refstart);
			// System.out.println("difference = " + diff);

			if (vcposition == refstart) {
				// if (size<450) // line 455
				if (size < refsize) {
					diff = refsize - size;
					for (int i = 0; i < diff; i++) {
						resultdtseq = resultdtseq + "P";
						resultedqual = resultedqual + "P";
					}
				}
				if (size > refsize) {
					resultdtseq = resultdtseq.substring(0, refsize);
					resultedqual = resultedqual.substring(0, refsize);
				}
			}

			if (vcposition > refstart) {
				// the data starting position is after the reference position : add
				// Ps in front
				diff = vcposition - refstart;
				// System.out.println("fragment  position = " + vcposition);
				// System.out.println("reference position = " + refstart);
				// System.out.println("difference = " + diff);
				for (int i = 0; i < diff; i++) {
					resultdtseq = "P" + resultdtseq;
					resultedqual = "P" + resultedqual;
				}
				size = resultdtseq.length();

				if (size < refsize) {
					diff = refsize - size;
					for (int i = 0; i < diff; i++) {
						resultdtseq = resultdtseq + "P";
						resultedqual = resultedqual + "P";
					}
				}

				if (size > refsize) {
					resultdtseq = resultdtseq.substring(0, refsize);
					resultedqual = resultedqual.substring(0, refsize);
				}
			}

			if (vcposition < refstart) {
				diff = refstart - vcposition;
				// System.out.println("fragment  position = " + vcposition);
				// System.out.println("reference position = " + refstart);
				// System.out.println("difference = " + diff);
				resultdtseq = resultdtseq.substring(diff, resultdtseq.length());
				resultedqual = resultedqual.substring(diff, resultedqual.length());

				size = resultdtseq.length();

				if (resultdtseq.length() < refsize) {
					diff = refsize - resultdtseq.length();
					for (int i = 0; i < diff; i++) {
						resultdtseq = resultdtseq + "P";
						resultedqual = resultedqual + "P";
					}
				}
				if (resultdtseq.length() > refsize) {
					resultdtseq = resultdtseq.substring(0, refsize);
					resultedqual = resultedqual.substring(0, refsize);
				}
			}
			reseqqual = resultdtseq + "\t" + resultedqual;

			return reseqqual;
		}
		
//		EXTRACT EXTRACT (END)--------------------------------------------------------------------------\\	
	// -------------------------------------------------------------------------------------------\\
		
	
	
	
	
	
	
	
	// Parse Reference Human Genome

    /**
     * IndependentGenomeRefParser method.
     * The method is used in 'noisetab' command in the first interval-sites preparative 
     * step This tskes place before the creation of noise tables. Leuicippus needs the 
     * reference sequence(normal) retrieved from Genome reference in order to know what
     * nucleotide is expected. While creating table this information is used to fill
     * the reference field in each table row that corresponds to a particular position. 
     * The interval(file) path, and the Genome reference (file) path  are passed 
     * as parameter here. 'IndependentGenomeRefParser' method reterieves all sites information 
     * using the 'ReadTheFile' method; next it performs sites-sort and merges overlapping sites using 
     * 'MergeOverlappinSites' and other methods call from MergeOverlappinSites. When information
     * about sites is ready(merged and sorted) the array of included chromosomes is retrieved 
     * and is used to generate an array list the contains site arrays of same chromosome.
     * The parser reads the Genome file and gradually reconstracts the current chromosome
     * (using a String vector. Each Vector element represent a Genome file line; The length of
     * the line is to the chromosomal comulative length. Using this information and the index of
     * Vector element the reference sequence is retrieved for each prepared site for this chromosome.
     * If there are not sites for this chromosome then the next chromosome is reconstracted in the 
     * same way. A new line for each site(merged) with included reference sequence is generated  
     * and is added to the returned String Vector(String elements contain tab separated words)  
     * @param prm_pth
     *            :String that holds the path for input file with sites
     * @param genref_pth
     *            :String that holds the path for Genome reference file The
     *            genome reference file must hold all chromosomes included in
     *            sites in ascending order first line is the lowest. Last
     *            chromosome is the highest.
     * @return Vector<String> : A modified sites Vector that holds also
     *         reference sequences from Left_start to Right_start.
     *
     * @throws FileNotFoundException
     * @throws IOException
     * @throws InterruptedException
     *
     *             sites format(tab separated file):
     *             ================================= 
     *             CHROM POS     REF ALT Left_start Right_start 
     *             2     5627740  C   G   5627478    5627921
     *             2     6912245  T   C   6912108    6912479 
     *             2     16266496 G   T   16266417   16266808
     *             
     *             
     *             Merged Sites format :
     *     CHROM POS                       REF     ALT     Left_start            Right_start             
     *     2     5627740,5627743,5627747   C,T,A   G,G,C   Minimum-of-left-start maximmun-of-Right_start
     *
     **/


    public static Vector<String> IndependentGenomeRefParser(String prm_pth,
            String genref_pth) throws FileNotFoundException, IOException,
            InterruptedException {

        System.out.println("Reference Parser Start");

        List<String[]> sites = new ArrayList<String[]>();

        File pr = new File(prm_pth);
        String prfnm = pr.getName();
        Vector<String> prmsits = new Vector<String>(); 
        // initial next sorted and next merge sites vector
        Vector<String> prmsfnl = new Vector<String>(); 
        // final returned sites vector (sorted merged and with sequences)
        Vector<String> chromosomes = new Vector<String>(); 
        // vector that controls the chromosomal flow of the loop
        // it gives sign for what list of sites(grouped in 
        // respect of chromosome) will be processed retrieve sort and merge 
        prmsits = readTheFileIncludeFirstLine(prm_pth, prfnm);
        prmsits = SortSites(prmsits);
	//for(int u=0; u<prmsits.size(); u++)
		//System.out.println(prmsits.get(u));
        //  prmsits = MergeSamePosSites(prmsits);


        prmsits = MergeOverlappinSites(prmsits);
        //for (int i = 0; i < prmsits.size(); i++)
            //System.out.println(prmsits.get(i));
	// Exit();
        // Create a sorted set of primers chromosome names in order to be used
        // as a guide for processing the sites in a chromosomal packets
        // (next use example: go through each chromosome found in sites
        // and retrieve the chromosome reference; process all sites of that
        // chromosome
        // go to the next chromosome etc.
        String sln = "", chrm = "", prchrom = "", crchrom = "";
        for (int k = 0; k < prmsits.size(); k++) {
            sln = prmsits.get(k);
            String[] ws = sln.split("\t");
            crchrom = ws[0];
            if (!(crchrom.equals(prchrom)))
                chromosomes.add(crchrom);
            prchrom = crchrom;
        }
        for(int i=0; i<chromosomes.size(); i++)
            System.out.println(chromosomes.get(i));
        
        // Populate the sorted primer/sites list(Each list is for a particular
        // chromosome
        Vector<String> crsit = new Vector<String>();
        String cm = "", stln = "", stcm = "";
        for (int i = 0; i < chromosomes.size(); i++) {
            cm = chromosomes.get(i);

            for (int k = 0; k < prmsits.size(); k++) {
                stln = prmsits.get(k);
                String[] wds = stln.split("\t");
                stcm = wds[0];
                if (stcm.equals(cm))
                    crsit.add(stln);
            }

            String[] arvc = new String[crsit.size()];
            crsit.toArray(arvc);
            sites.add(arvc);
            crsit = new Vector<String>();
        }
        // End Populate the sorted primer/sites list(Each list is for a
        // particular chromosome

        File gfl = new File(genref_pth);
        // >1 dna:chromosome chromosome:GRCh37:1:1:249250621:1
        String line = "", chromline = "", chrmnums = "", dnachck = "", 
                chmnmszsts = "", cnms = "", cstrss = "", cends = "", cactst;
        int counter = 0, start = 0, end = 0, ext = 0;
        String cr = "";
        String ResSeq = "";
        String CHROM = "", POS = "", REF = "", ALT = "", Left_start = "",
                Right_start = "", curchrom = "";
        String chomnmline = "", firstln = "", test = "";
        int crst = 0, crrfstrt = 0, stractstart = 0, 
                endactstart = 0, crcnt = 0;
        int chromcounter = -1, bil = 0, Chromosomelength = 0;

        Vector<String> chrtst = new Vector<String>();

        String[] curentprims;
        String hit = "";
        String ht = "no"; // hit= the position is not found yet
                            // yes the start position is found
        boolean chromlinefound = false;

        BufferedReader br = null;
        if (isFileGZipped(gfl)) {
            GZIPInputStream gzip = new GZIPInputStream(new FileInputStream(
                    genref_pth));
            br = new BufferedReader(new InputStreamReader(gzip));
        } else if (!(isFileGZipped(gfl))) {
            InputStreamReader reffl = new InputStreamReader(
                    new FileInputStream(genref_pth));
            br = new BufferedReader(reffl);
        }
        boolean foundom = false;
        while ((line = br.readLine()) != null) 
        {
            if ((counter == chromosomes.size()) && (foundom == true))
                break;
            cr = chromosomes.get(counter);
            curentprims = sites.get(counter); // get all sites for a particular chromosome
            //ext = 1; // exit the main loop
            foundom = false;
 
            if (!(line == null))
                if (line.length() > 0)
                {
                    if ((line.charAt(0) == '>') || (chromlinefound == true))
                    {
                        if (line.charAt(0) == '>') 
                        {
                            System.out.println(line);
                            chomnmline = line; // question for chromosome one
                                                // and answer chromosome one
                            chromlinefound = true;
                        }
                        String[] words = chomnmline.split(" ");
                        if ( (chromlinefound == true) && (line.charAt(0) == '>'))
                            firstln = line;
                        
                        if (words.length > 1)
                        {
                            chrmnums = "";
                            dnachck = "";
                            chmnmszsts = "";
                            // Retrieve the chromosome(1,2,..,22,X,Y,..)
                            chrmnums = words[0].substring(1, words[0].length());
                            dnachck = words[1];
                            chmnmszsts = words[2];
                            String[] chinf = chmnmszsts.split(":");
                            cnms = chinf[2];
                            cstrss = chinf[3];
                            cends = chinf[4];
                            cactst = chinf[5];
                            if (chrmnums.equals(cr))
                            {
                                foundom = true;
                                // System.out.println("Here chrmnums = " +
                                // chrmnums);
                                chrtst = new Vector<String>();
                                //if(chromosomes.size()>1)
                                counter = counter + 1;
                                curchrom = "";
                                chromcounter = chromcounter + 1;
                                chromline = "";

                                // the following loop ends when next chromosome
                                // informative line found
                                // and the if informative line completes its
                                // single participation
                                // since after that the counter is always > 0
                                // the next informatyive chrom line ia stored in
                                // the chomnmline String
                                // this is used in the next big if (if counter >
                                // 0)
                                 chromlinefound = false;
                                if (!(firstln.isEmpty()))
                                {
                                    // chrtst.add(firstln);// dont add
                                    // chrom-info in order to keeep position =
                                    // prev count + index + 1
                                    System.out.println("cr  " + chomnmline);
                                    System.out.println("ppp " + firstln);
                                    if (counter != 0)
                                        Chromosomelength = Chromosomelength
                                                + firstln.length();
                                }
                                firstln = "";
                                bil = 0;
                                if (!(line.charAt(0) == '>'))
                                {
                                    chrtst.add(line);
                                    Chromosomelength = Chromosomelength
                                            + line.length();   
                                }
                                while ((bil != 1)
                                        && ((line = br.readLine()) != null))
                                {
                                    // curchrom = curchrom + line;
                                    if (!(line.charAt(0) == '>')) {
                                        chrtst.add(line);
                                        Chromosomelength = Chromosomelength
                                                + line.length();
                                    } else if (line.charAt(0) == '>') {
                                        chomnmline = line; // store the chro
                                                            // info line
                                        System.out.println("--- " + line);
                                        bil = 1;
                                        chromlinefound = true;
                                    }
                                }
                                // chromosome vector is Ready
                                for (int sr = 0; sr < curentprims.length; sr++) {
                                    String[] prwords = curentprims[sr]
                                            .split("\t");
                                    CHROM = prwords[0];
                                    POS = prwords[1];
                                    REF = prwords[2];
                                    ALT = prwords[3];
                                    Left_start = prwords[4];
                                    Right_start = prwords[5];
                                    start = Integer.parseInt(Left_start);
                                    end = Integer.parseInt(Right_start);
                                    crcnt = 0;
                                    ht = "no";
                                    hit = "";
                                    for (int f = 0; f < chrtst.size(); f++) {
                                        crst = crcnt; // start point
                                        crcnt = crcnt + chrtst.get(f).length();

                                        if (start <= crcnt) {
                                            hit = hit + chrtst.get(f);
                                            if (ht.equals("no")) {
                                                crrfstrt = crst;
                                                stractstart = start
                                                        - (crrfstrt + 1); // find local start index
                                                endactstart = stractstart
                                                        + (end - start) + 1; // find local end index
                                            }
                                            ht = "yes";
                                        }
                                        if (end <= crcnt) {
                                            hit = hit + chrtst.get(f);
                                            break;
                                        }
                                    }
                                    
                                    // System.out.println(hit.length());
                                    ResSeq = hit.substring(stractstart,
                                            endactstart);
                                    if (ResSeq.length() > 0) {
                                        test = cr + "\t" + start + "\t" + end
                                                + "\t" + ResSeq + "\t" + POS
                                                + "\t" + REF + "\t" + ALT;
                                        prmsfnl.add(test);
                                        System.out.println(test);
                                    } else {
                                        System.out
                                                .println(" sequence : "
                                                        + cr
                                                        + "\t"
                                                        + start
                                                        + "\t"
                                                        + end
                                                        + " site = "
                                                        + POS
                                                        + "\t"
                                                        + REF
                                                        + "\t"
                                                        + ALT
                                                        + "\n has zero sequence return.");
                                        missed_sites.add(curentprims[sr]);
                                    }
                                }
                            }
                        }// if(words.length==3)
                        if (words.length == 1)
                        {
                            System.out.println("length=1");
                            
              
                        }
                        
                    }
                }
            
            
            
            
            
        }
        if (missed_sites.size() > 0) {
            System.out.println("The following sites will not be processed : ");
            for (int i = 0; i < missed_sites.size(); i++)
                System.out.println(missed_sites.get(i));
        }
        return prmsfnl;
    }

/**
 * Method that identifies overlapping sites referred to same chromosome.
 * It creates a common information about overlapping sites and updates the
 * sites Vector accordingly. The new information for merged sites includes
 * the original in comma separated form in the general tab separated words
 * Added on 10/29/2015 :Sites in same position with different alternatives
 * will be represented with only one word with more than one letters in the 
 * alternative field.
 * Example :alternative A and alternative C = AC
 * @param tbdt : Vector<String> It holds all sites; some of them present overlaps 
 * @return remained : Vector<String> Non overlapping sites. Some of the sites
 * have been merged.
 **/
	public static Vector<String> MergeOverlappinSites(Vector<String> tbdt) {

		String line = "", lstrt = "", rstrt = "", lincpc = "", lincpp = "", lincp = "";
		String crms = "", crm = "", crmq;
		String CHROM = "", POS = "", REF = "", ALT = "", Left_start = "", Right_start = "";
		String CurPOS = "", CurREF = "", CurALT = "";
		String FinPOS = "", FinREF = "", FinALT = "";

		// Create String line to hold all chromosomes separated by space
		// This will generate an array of chromosome names from 1 to Y
		for (int i = 1; i < 23; i++)
			crms = crms + i + " ";
		crms = crms + "X Y";
		Vector<String> remained = new Vector<String>();
		Vector<String> remained2 = new Vector<String>();
		Vector<String> chrmsite = new Vector<String>();
		int lstrti = 0, rstrti = 0, diffi = 0;
		int prst = 0, pren = 0, crst = 0, cren = 0;

		String[] chroms = crms.split(" "); // Chromosome array

		remained = new Vector<String>();
		for (int j = 0; j < chroms.length; j++) // Go through each chromosome
		{
			crmq = chroms[j];
			chrmsite = new Vector<String>();
			for (int i = 0; i < tbdt.size(); i++) 
			{
				line = tbdt.get(i);
				String[] ws = line.split("\t");
				crm = ws[0];
				if (crmq.equals(crm))
				{
					chrmsite.add(line); // populate this vector with all sites
										// for a particular chromosome
					// System.out.println("---> " + line);
				}
			}
			for (int k = 1; k < chrmsite.size(); k++) 
			{

				lincpp = chrmsite.get(k - 1); // previous site
				lincpc = chrmsite.get(k); // current site
				
				String[] wordsp = lincpp.split("\t");
				CHROM = wordsp[0];
				POS = wordsp[1];
				REF = wordsp[2];
				ALT = wordsp[3];
				Left_start = wordsp[4];
				Right_start = wordsp[5];

				String[] wordsc = lincpc.split("\t");

				CurPOS = wordsc[1];
				CurREF = wordsc[2];
				CurALT = wordsc[3];

				prst = Integer.parseInt(wordsp[4]);
				pren = Integer.parseInt(wordsp[5]);

				crst = Integer.parseInt(wordsc[4]);
				cren = Integer.parseInt(wordsc[5]);

				// if((crst<prst) || (cren<pren))
				// System.out.println("not sorted");
				// System.out.println("Prev POS =  " + POS + "    Cur POS" + CurPOS);	
				if( (CurPOS.equals(POS)) && (!(ALT.equals(CurALT))))
				{
					//System.out.println("same position found");
					FinPOS=POS;
					FinALT = ALT + CurALT;
					FinALT = removeDuplicates(FinALT);
					FinREF = REF;
					lincp = CHROM + "\t" + FinPOS + "\t" + FinREF + "\t"
							+ FinALT + "\t" + prst + "\t" + cren;
					chrmsite.set(k - 1, lincp);
					chrmsite.removeElementAt(k);
					//System.out.println("New : " + lincp);
					k = 0;	
                 }
			}
			for (int o = 0; o < chrmsite.size(); o++)
				remained.add(chrmsite.get(o));
		}

		remained2 = new Vector<String>();

		for (int j = 0; j < chroms.length; j++) // Go through each chromosome
		{
			crmq = chroms[j];
			chrmsite = new Vector<String>();
			for (int i = 0; i < remained.size(); i++) 
			{
				line = remained.get(i);
				String[] ws = line.split("\t");
				crm = ws[0];
				if (crmq.equals(crm))
				{
					chrmsite.add(line); // populate this vector with all sites
										// for a particular chromosome
					// System.out.println("---> " + line);
				}
			}

			for (int k = 1; k < chrmsite.size(); k++) 
			{

				lincpp = chrmsite.get(k - 1); // previous site
				lincpc = chrmsite.get(k); // current site

				String[] wordsp = lincpp.split("\t");
				CHROM = wordsp[0];
				POS = wordsp[1];
				REF = wordsp[2];
				ALT = wordsp[3];
				Left_start = wordsp[4];
				Right_start = wordsp[5];

				String[] wordsc = lincpc.split("\t");

				CurPOS = wordsc[1];
				CurREF = wordsc[2];
				CurALT = wordsc[3];

				prst = Integer.parseInt(wordsp[4]);
				pren = Integer.parseInt(wordsp[5]);

				crst = Integer.parseInt(wordsc[4]);
				cren = Integer.parseInt(wordsc[5]);

				// if((crst<prst) || (cren<pren))
				// System.out.println("not sorted");

				if (  (crst < pren) && (!(CurPOS.equals(POS))))
				{
					//System.out.println("overlap found");
					//System.out.println(lincpp);
					//System.out.println(lincpc);
					FinPOS = POS + "," + CurPOS;
					FinREF = REF + "," + CurREF;
					FinALT = ALT + "," + CurALT;
					lincp = CHROM + "\t" + FinPOS + "\t" + FinREF + "\t"
							+ FinALT + "\t" + prst + "\t" + cren;
					chrmsite.set(k - 1, lincp);
					chrmsite.removeElementAt(k);
					//System.out.println("New : " + lincp);
					k = 0;
				}

			}

			for (int o = 0; o < chrmsite.size(); o++)
				remained2.add(chrmsite.get(o));
			// System.out.println("\n\n");
		}

		// int min=2000;
		// Vector<Double> db = new Vector<Double>();
		return remained2;
	}

	/**

 	 * The following method is useful when alternative fields contain more than 
	 * one nucleotides. It is called from 'MergeOverlappinSites' method.
	 * The combination of these methods is part of IndependentGenomeReferenceParser
	 * process(retrieve reference sequences (noisetab command).   
	 * The removeDuplicates method accepts a String with Letters
	 * (four possible nucleotides ('A', 'C', 'T', 'G')). The method removes
	 * the more than one nucleotide same instances from the passed String.
	 * The returned String contains only one instance of included nucleotides
	 * in a particular order
	 * @param str
	 * @return results 
	 * Prerequisite: The passed String must not be nu
	 */
	public static String removeDuplicates(String str)
	{
		String results="";
		for(int i=0; i<str.length(); i++)
		{
			if(str.charAt(i)=='A')
			{
				results = results+"A";
				break;
			}
		}
		for(int i=0; i<str.length(); i++)
		{
			if(str.charAt(i)=='T')
			{
				results = results+"T";
				break;
			}
		}
		for(int i=0; i<str.length(); i++)
		{
			if(str.charAt(i)=='C')
			{
				results = results+"C";
				break;
			}
		}
		for(int i=0; i<str.length(); i++)
		{
			if(str.charAt(i)=='G')
			{
				results = results+"G";
				break;
			}
		}
		return results;
	}

/**
 * Method that identifies overlapping sites refered to same chromosome.
 * It creates a common information about overlapping sites and updates the
 * sites Vector accordingly. The new information for merged sites includes
 * the original in comma separated form in the general tab separated words 
 * @param tbdt : Vector<String> It holds all sites some of them present overlaps 
 * @return remained : Vector<String> Non overlapping sites. Some of the sites
 *  have been merged.
 */
	public static Vector<String> MergeOverlappinSitesBack(Vector<String> tbdt) {
		String line = "", lstrt = "", rstrt = "", lincpc = "", lincpp = "", lincp = "";
		String crms = "", crm = "", crmq;
		String CHROM = "", POS = "", REF = "", ALT = "", Left_start = "", Right_start = "";
		String CurPOS = "", CurREF = "", CurALT = "";
		String FinPOS = "", FinREF = "", FinALT = "";

		// Create String line to hold all chromosomes separated by space
		// This will generate an array of chromosome names from 1 to Y
		for (int i = 1; i < 23; i++)
			crms = crms + i + " ";
		crms = crms + "X Y";
		Vector<String> remained = new Vector<String>();
		Vector<String> chrmsite = new Vector<String>();
		int lstrti = 0, rstrti = 0, diffi = 0;
		int prst = 0, pren = 0, crst = 0, cren = 0;

		String[] chroms = crms.split(" "); // Chromosome array
		for (int j = 0; j < chroms.length; j++) // Go through each chromosome
												// and
												// generate a
		{
			crmq = chroms[j];

			chrmsite = new Vector<String>();
			for (int i = 0; i < tbdt.size(); i++) {
				line = tbdt.get(i);

				String[] ws = line.split("\t");
				crm = ws[0];
				if (crmq.equals(crm)) {
					chrmsite.add(line); // populate this vector with all sites
										// for a particular chromosome
					// System.out.println("---> " + line);
				}
			}
			for (int k = 1; k < chrmsite.size(); k++) {

				lincpp = chrmsite.get(k - 1); // previous site
				lincpc = chrmsite.get(k); // current site
				String[] wordsp = lincpp.split("\t");
				CHROM = wordsp[0];
				POS = wordsp[1];
				REF = wordsp[2];
				ALT = wordsp[3];
				Left_start = wordsp[4];
				Right_start = wordsp[5];

				String[] wordsc = lincpc.split("\t");

				CurPOS = wordsc[1];
				CurREF = wordsc[2];
				CurALT = wordsc[3];

				prst = Integer.parseInt(wordsp[4]);
				pren = Integer.parseInt(wordsp[5]);

				crst = Integer.parseInt(wordsc[4]);
				cren = Integer.parseInt(wordsc[5]);

				// if((crst<prst) || (cren<pren))
				// System.out.println("not sorted");
				if ((crst < pren)) {
					System.out.println("overlap found");
					System.out.println(lincpp);
					System.out.println(lincpc);
					FinPOS = POS + "," + CurPOS;
					FinREF = REF + "," + CurREF;
					FinALT = ALT + "," + CurALT;
					lincp = CHROM + "\t" + FinPOS + "\t" + FinREF + "\t"
							+ FinALT + "\t" + prst + "\t" + cren;
					chrmsite.set(k - 1, lincp);
					chrmsite.removeElementAt(k);
					System.out.println("New : " + lincp);
					k = 0;
				}
			}
			for (int o = 0; o < chrmsite.size(); o++)
				remained.add(chrmsite.get(o));
			// System.out.println("\n\n");
		}
		// int min=2000;
		// Vector<Double> db = new Vector<Double>();
		return remained;
	}

/**
	 * 'SortSites' method accepts a String Vector that holds all unsorted sites.
     * It sorts them in ascending order.
	 * @param chnsites Vector<String> 
	 *            : String Vector with unsorted sites
	 * @return results : Vector<String> with sorted sites
*/
	
	public static Vector<String> SortSites(Vector<String> chnsites) {

		Vector<String> results = new Vector<String>();
		Vector<String> chmsts = new Vector<String>();
		// build a sorted array of chromosomes. ( 1, 2, 3, ... 22, X, Y )
		String crms = "", crm = "", crmq;
		for (int i = 1; i < 23; i++)
			crms = crms + i + " ";
		crms = crms + "X Y";
		String[] chroms = crms.split(" "); // Chromosome array
		// String header = "CHROM	POS	REF	ALT	Left_start	Right_start";
		String header = "";
		header = chnsites.get(0);
		String stlne = "", CHROM = "", POS = "", REF = "", ALT = "", Left_start = "", Right_start = "";
		if (chnsites.size() > 0) {
			header = chnsites.get(0);
			results.add(header);
//      loop that goes through all cromosomes in ascending order 
//		retrieves all sites of current chromosome
//      sorts those sites according to 'Left Start' word(using 'SortChromSites' method
//	    inserts the sorted sited in the results Vector
//      When the loop ends results Vector has all sites sorted according to chromosom
//      and to Left Start.			
			for (int j = 0; j < chroms.length; j++) // Go through each
													// chromosome and
			{
				chmsts = new Vector<String>();
				crm = chroms[j];
				for (int i = 0; i < chnsites.size(); i++) {
					stlne = chnsites.get(i);
					String[] ws = stlne.split("\t");
					CHROM = ws[0];
					if (crm.equals(CHROM))
						chmsts.add(stlne);
				}
				
				chmsts = SortChromSites(chmsts);
				for (int r = 0; r < chmsts.size(); r++) {
					results.add(chmsts.get(r));
				}
			}
		}
		System.out.println("Sites Sorted (Chromosome, Left Start).");
		return results;
	}
	
/**
 * The method accepts a String Vector. Its elements are information for 
 * sites of the same chromosome.
 * A site information is a tab delimited String that holds the following
 * information : Chromosome, Position(site), Reference Base, 
 * Alternative Base, Left Start, Right Start.
 * The method sorts the vector elements in ascending order according Left Start, 
 * Example sorted: 
 *  Chromosome 1 Left Start 1000
 *  Chromosome 1 Left Start 3020
 *  Chromosome 1 Left Start 4992
 *  ...
 *  
 *  The method uses a paralel vector that holds all Left Start words(Strings) converted
 *  to integers. Next a bubble sort attempts to sort that Vector along
 *  with the original Vector.
 *  The method is called from the 'SortSites' method which is called from
 *  'IndependentGenomeReferenceParser' method.
 *  The "SortSites" method sorts the sites according to cromosomes. 
 * @param chnsites : Vector<String> unsorted
 * @return results Vector<String> sorted
 */
	public static Vector<String> SortChromSites(Vector<String> chnsites) 
	{

		Vector<Integer> lstts = new Vector<Integer>();
		String CHROM = "", POS = "", REF = "", ALT = "", Left_start = "";
		String line = "";
		int lstrti = 0;
		for (int i = 0; i < chnsites.size(); i++) {
			line = chnsites.get(i);
			String[] wds = line.split("\t");
			Left_start = wds[4];
			lstrti = Integer.parseInt(Left_start);
			lstts.add(lstrti);
		}
		int prvlsti = 0, curlstsi = 0;
		String prvln = "", curln = "";

		for (int i = 1; i < chnsites.size(); i++) {

			prvln = chnsites.get(i - 1);
			curln = chnsites.get(i);

			prvlsti = lstts.get(i - 1);
			curlstsi = lstts.get(i);

			if (prvlsti > curlstsi) {
				lstts.set(i - 1, curlstsi);
				lstts.set(i, prvlsti);
				chnsites.set(i - 1, curln);
				chnsites.set(i, prvln);
				i = 0;
			}
		}
		return chnsites;
	}

	/**
	 * Method that tests if a file is gzipped
	 * @param tstfl : File
	 * @return boolean : true if the file is gzipped; false otherwise
	 * 
	 */
	public static boolean isFileGZipped(File tstfl) {
		int mg = 0;
		try {
			RandomAccessFile rnd = new RandomAccessFile(tstfl, "r");
			mg = rnd.read() & 0xff | ((rnd.read() << 8) & 0xff00);
			rnd.close();
		} catch (Throwable e) {
			e.printStackTrace(System.err);
		}
		return mg == GZIPInputStream.GZIP_MAGIC;
	}

	/**
	 * Method used to test if the collection directory for the output file
	 * exists
	 * @param flpth
	 *            : The path with the file name
	 * @return res : true if parent directory exists; false otherwise
	 */
	public static boolean testParentDirectoryExistence(String flpth) {
		boolean res = false;
		File fl = new File(flpth);
		File parent = fl.getParentFile();
		System.out.println("DIR = " + parent.getName());
		if (parent.exists())
			res = true;
		return res;
	}

	/**
	 * Method that tests the existence of a particular file
	 * @param flpth
	 *            : String : file path including the name
	 * @return boolean res : boolean (it is true when indeed the file exists; otherwise
	 *         it is false)
	 */
	public static boolean testFileExistence(String flpth) {
		boolean res = false;
		File fl = new File(flpth);
		if (fl.exists() && !fl.isDirectory())
			res = true;
		return res;
	}
	
	/**
	 * 'isPosUnsignInteger' method accepts a String as a parameter and checks if
	 * it is a positive integer. If yes then it returns true; otherwise it
	 * returns false
	 * @param qual
	 *            : String
	 * @return res : boolean requirements : none
	 */
	public static boolean isPosUnsignInteger(String qual)
	{
		boolean res = false;
		if (qual.matches("\\d+"))
			res = true;
		return res;
	}

	public static boolean isNumeric(String inputData) {
		return inputData.matches("[-+]?\\d+(\\.\\d+)?");
	}

	/**
	 * Method that exits wen its called : Mainly it is used when there is an
	 * error in input values
	 * 
	 */
	public static void Exit() {
		System.exit(1);
	}

	/**
	 * Method that returns all files contained in a particular directory
	 * @param folder
	 *            : path to the directory
	 * @return filenames : String vector that contains all file names in that
	 *         directory (only file names not directories)
	 */
	public static Vector<String> listFilesForFolder(final File folder) {
		int nothing = 0;
		Vector<String> filenames = new Vector<String>();
		for (final File fileEntry : folder.listFiles()) {
			nothing = 0;
			if (fileEntry.isDirectory()) {
				// listFilesForFolder(fileEntry);
				nothing = 1;
			} else {
				filenames.add(fileEntry.getName());
				// System.out.println(fileEntry.getName());
			}
		}
		return filenames;
	}

	/**
	 * Method that sorts a String Vector using language background methods. It
	 * first converts the Vector to a String array; next it sorts the Array;
	 * finally it converts the resulted Array to a String Vector
	 * 
	 * @param results
	 *            : input String Vector
	 * @return res : output sorted String Vector
	 */
	public static Vector<String> sortStrings(Vector<String> results) {
		Vector<String> res = new Vector<String>();
		ArrayList<String> list = new ArrayList<String>(results);
		Collections.sort(list);

		String[] ar = list.toArray(new String[0]);
		// for(int i=0; i<ar.length; i++)
		// System.out.println(ar[i]);
		res = new Vector<String>(Arrays.asList(ar));
		return res;
	}

	/**
	 * Method that Reads a file line by line and stores each line as a String
	 * Vector element. Then it returns the String Vector. The method requires
	 * that the file exists, otherwise it throws an exception
	 * 
	 * @param fnpath
	 *            : The file name with the path
	 * @param filename
	 *            : The file name
	 * @return vct : String Vector that holds the file lines.
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public static Vector<String> readTheFileIncludeFirstLine(String fnpath,
			String filename) throws IOException, InterruptedException {
		System.out.println("	    Start ReadTheFile method");
		System.out.println("	    Start Read File : " + filename);
		// Thread.sleep(4000);
		File fl = new File(fnpath);
		Vector<String> vct = new Vector<String>();
		BufferedReader br;
		String line = "";
		int count = -1;
		br = new BufferedReader(new FileReader(fnpath));
		while ((line = br.readLine()) != null) {
			count = count + 1;
			if ((count >= 0) && (line.length() > 0))
				vct.add(line);
		}
		System.out.println("		file size = " + fl.length() + " Bytes");
		System.out.println("		file size = "
				+ ((double) fl.length() / (double) 1000000) + " MB");

		br.close();
		// System.out.println("First line :" + vct.get(0));
		// System.out.println("Last line :" + vct.get(vct.size()-1));
		System.out.println("	    End  Read  File : " + filename);
		System.out.println("       ....End Read The File Method!\n");
		return vct;
	}

	/**
	 * 'readTheFile' method reads a file line by line and stores each line as a
	 * String Vector element. Then it returns the String Vector(The first file
	 * line is skipped. The method requires that the file exists, otherwise it
	 * throws an exception
	 * 
	 * @param fnpath
	 *            : The file name with the path
	 * @param filename
	 *            : The file name
	 * @return vct : String Vector that holds the file lines.
	 * @throws IOException
	 * @throws InterruptedException
	 */

	public static Vector<String> readTheFile(String fnpath, String filename)
			throws IOException
	// InterruptedException
	{

		System.out.println("	    Start ReadTheFile method");
		System.out.println("	    Start Read File : " + filename);
		// Thread.sleep(4000);
		File fl = new File(fnpath);
		Vector<String> vct = new Vector<String>();
		BufferedReader br;
		String line = "";
		int count = -1;
		br = new BufferedReader(new FileReader(fnpath));

		while ((line = br.readLine()) != null) {
			count = count + 1;
			if ((count > 0) && (line.length() > 0))
				vct.add(line);
		}

		System.out.println("		file size = " + fl.length() + " Bytes");
		System.out.println("		file size = "
				+ ((double) fl.length() / (double) 1000000) + " MB");

		br.close();
		// System.out.println("First line :" + vct.get(0));
		// System.out.println("Last line :" + vct.get(vct.size()-1));
		System.out.println("	    End  Read  File : " + filename);
		System.out.println("       ....End Read The File Method!\n");
		return vct;
	}

	/**
	 * Method that writes all elements of a String vector to a file.
	 * @param fnpath : String : destination file path name.
	 * @param vc : Vector<String> it holds results data String lines. 
	 * @throws IOException
	 * Prerequisites : Permissions must be writeable in the destination directory 
	 */
	public static void writeToFile(String fnpath, Vector<String> vc)
			throws IOException {
		System.out.println("writeToFile method.");
		File result = new File(fnpath);
		FileWriter fw = new FileWriter(fnpath);
		BufferedWriter bw = new BufferedWriter(fw);
		for (int i = 0; i < vc.size(); i++)
			bw.write(vc.elementAt(i) + "\n");
		bw.close();
	}
}
