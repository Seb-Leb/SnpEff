package org.snpeff.snpEffect.commandLine;

import java.util.HashMap;
import java.util.List;

import org.snpeff.SnpEff;
import org.snpeff.align.SmithWaterman;
import org.snpeff.codons.CodonTable;
import org.snpeff.codons.CodonTables;
import org.snpeff.fileIterator.FastaFileIterator;
import org.snpeff.interval.Gene;
import org.snpeff.interval.Transcript;
import org.snpeff.snpEffect.Config;
import org.snpeff.util.Gpr;
import org.snpeff.util.Timer;

/**
 * Command line: Calculate coding sequences from a file and compare them to the ones calculated from our data structures
 *
 * Note: This is done in order to see potential incompatibility
 *       errors between genome sequence and annotation.
 *
 * @author pcingola
 */
public class SnpEffCmdCds extends SnpEff {

	public static boolean onlyOneError = false; // This is used in some test-cases
	public static double maxErrorPercentage = 0.01; // Maximum allowed error is 1% (otherwise test fails)
	public static int MAX_ALIGN_LENGTH = 33000;

	boolean storeAlignments; // Store alignments (used for some test cases)
	boolean checkNumOk = true;
	int totalErrors = 0;
	int totalOk = 0;
	int totalWarnings = 0;
	int totalNotFound = 0;
	String cdsFile = "";
	HashMap<String, String> cdsByTrId;

	public SnpEffCmdCds() {
	}

	public SnpEffCmdCds(Config config) {
		this.config = config;
		cdsFile = config.getFileNameCds();
	}

	public SnpEffCmdCds(String genomeVer, String configFile, String cdsFile) {
		this.configFile = configFile;
		this.genomeVer = genomeVer;
		this.cdsFile = cdsFile;
	}

	void add(String trId, String seq, int lineNum, boolean check) {
		// Repeated transcript Id? => Check that Protein is the same
		if (check && (cdsByTrId.get(trId) != null) && (!cdsByTrId.get(trId).equals(seq))) //
			System.err.println("ERROR: Different CDS for the same transcript ID. This should never happen!!!"//
					+ "\n\tLine number   : " + lineNum //
					+ "\n\tTranscript ID : '" + trId + "'"//
					+ "\n\tProtein       : " + cdsByTrId.get(trId) //
					+ "\n\tProtein (new) : " + seq //
			);

		// Use whole trId
		cdsByTrId.put(trId, seq); // Add it to the hash
		if (debug) Gpr.debug("Adding cdsByTrId{'" + trId + "'} :\t" + seq);
	}

	/**
	 * Compare all CDS
	 */
	double cdsCompare() {
		int i = 0;

		if (verbose) {
			// Show labels
			System.err.println("\tLabels:");
			System.err.println("\t\t'+' : OK");
			System.err.println("\t\t'.' : Missing");
			System.err.println("\t\t'*' : Error");
			System.err.print("\t");
		}

		// Compare all genes
		for (Gene gene : config.getGenome().getGenes())
			for (Transcript tr : gene) {
				char status = ' ';

				boolean ok = false;
				String cds = tr.cds().toUpperCase();
				String mRna = tr.mRna().toUpperCase();
				String cdsReference = cdsByTrId.get(tr.getId());

				if (cdsReference != null) cdsReference = cdsReference.toUpperCase();

				if (cdsReference == null) {
					status = '.';
					if (debug) System.err.println("\nWARNING:Cannot find reference CDS for transcript '" + tr.getId() + "'");
				} else if (cds.isEmpty()) {
					status = '.';
					if (debug) System.err.println("\nWARNING:Empty CDS for transcript '" + tr.getId() + "'");
				} else if (cds.equals(cdsReference)) {
					status = '+';

					// Sanity check: Start and stop codons
					if ((cds != null) && (cds.length() >= 3)) {
						CodonTable ctable = CodonTables.getInstance().getTable(config.getGenome(), tr.getChromosomeName());

						// Check start codon
						String startCodon = cds.substring(0, 3);
						if (!ctable.isStart(startCodon)) {
							if (debug) System.err.println("\nWARNING: CDS for transcript '" + tr.getId() + "' does not start with a start codon:\t" + startCodon + "\t" + cds);
							totalWarnings++;
						}

						// Check stop codon
						String stopCodon = cds.substring(cds.length() - 3, cds.length());
						if (!ctable.isStop(stopCodon)) {
							if (debug) System.err.println("\nWARNING: CDS for transcript '" + tr.getId() + "' does not end with a stop codon:\t" + stopCodon + "\t" + cds);
							totalWarnings++;
						}
					}
				} else if (mRna.equals(cdsReference)) { // May be the file has mRNA instead of CDS?
					status = '+';
				} else if ((mRna.length() < cdsReference.length()) // CDS longer than mRNA? May be it is actually an mRNA + poly-A tail (instead of a CDS)
						&& cdsReference.substring(mRna.length()).replace('A', ' ').trim().isEmpty() // May be it is an mRNA and it has a ploy-A tail added
						&& cdsReference.substring(0, mRna.length()).equals(mRna) // Compare cutting poly-A tail
				) {
					// OK, it was a mRNA +  polyA
					status = '+';
				} else if ((mRna.length() > cdsReference.length()) // PolyA in the reference?
						&& mRna.substring(cdsReference.length()).replace('A', ' ').trim().isEmpty() //
						&& mRna.substring(0, cdsReference.length()).equals(mRna) //
				) {
					// OK, it was a mRNA +  polyA
					status = '+';
				} else if (cdsReference.indexOf(cds) >= 0) { // CDS fully included in reference?
					status = '+';
				} else {
					status = '*';

					if (debug || storeAlignments || onlyOneError) {
						// Create a string indicating differences
						SmithWaterman sw = new SmithWaterman(cds, cdsReference);
						if (Math.max(cds.length(), cdsReference.length()) < MAX_ALIGN_LENGTH) sw.align();

						int maxScore = Math.min(cds.length(), cdsReference.length());
						int score = sw.getAlignmentScore();

						if (debug || onlyOneError) {
							System.err.println("\nERROR: CDSs do not match for transcript " + tr.getId() //
									+ "\tStrand:" + tr.isStrandMinus()//
									+ "\tExons: " + tr.numChilds() //
									+ "\n" //
									+ String.format("\tSnpEff CDS  (%6d) : '%s'\n", cds.length(), cds.toLowerCase()) //
									+ String.format("\tReference   (%6d) : '%s'\n", cdsReference.length(), cdsReference.toLowerCase()) //
									+ "\tAlignment (Snpeff CDS vs Reference CDS)." //
									+ "\tScore: " + score //
									+ "\tMax. possible score: " + maxScore //
									+ "\tDiff: " + (maxScore - score) //
									+ "\n" + sw //
							);
							System.err.println("Transcript details:\n" + tr);
						}
						if (onlyOneError) {
							System.err.println("Transcript details:\n" + tr);
							throw new RuntimeException("Showing only one error!");
						}

					}
				}

				// Update counters
				switch (status) {
				case '.':
					totalNotFound++;
					break;

				case '+':
					totalOk++;
					ok = true;
					break;

				case '*':
					totalErrors++;
					break;

				case ' ':
					break;

				default:
					throw new RuntimeException("Unknown status '" + status + "'");
				}

				// Update transcript DnaCheck status
				if (ok) tr.setDnaCheck(true);

				// Show a mark
				if (verbose && (status != ' ')) {
					if (status == '.' && !tr.isProteinCoding()) {
						// OK, may be this is a non-protein coding transcript and that's why we could not find the CDS
					} else {
						System.out.print(status);
						i++;
						if (i % 100 == 0) System.out.print("\n\t");
					}

				}
			}

		double perc = ((double) totalErrors) / ((double) (totalErrors + totalOk));
		System.out.println("\n\tCDS check:\t" //
				+ config.getGenome().getVersion() //
				+ "\tOK: " + totalOk //
				+ "\tWarnings: " + totalWarnings //
				+ "\tNot found: " + totalNotFound //
				+ "\tErrors: " + totalErrors //
				+ "\tError percentage: " + (100 * perc) + "%" //
		);

		// Sanity check
		if (checkNumOk && totalOk <= 0) fatalErrorNoTranscriptsChecked();

		return perc;
	}

	/**
	 * Show an error message that actually helps to solve the problem
	 */
	void fatalErrorNoTranscriptsChecked() {
		StringBuilder sb = new StringBuilder();

		// Show some transcript IDs
		int maxTrIds = 20;
		sb.append("Transcript IDs from database (sample):\n" + sampleTrIds(maxTrIds));
		sb.append("Transcript IDs from database (fasta file):\n" + sampleTrIdsFasta(maxTrIds));
		fatalError("No CDS checked. This is might be caused by differences in FASTA file transcript IDs respect to database's transcript's IDs.\n" + sb);
	}

	/**
	 * Parse command line arguments
	 */
	@Override
	public void parseArgs(String[] args) {
		for (int i = 0; i < args.length; i++) {

			// Argument starts with '-'?
			if (args[i].startsWith("-")) {
				if ((args[i].equals("-c") || args[i].equalsIgnoreCase("-config"))) {
					if ((i + 1) < args.length) configFile = args[++i];
					else usage("Option '-c' without config file argument");
				} else if (args[i].equals("-v") || args[i].equalsIgnoreCase("-debug")) {
					debug = true;
				} else usage("Unknown option '" + args[i] + "'");
			} else if (genomeVer.isEmpty()) genomeVer = args[i];
			else if (cdsFile.isEmpty()) cdsFile = args[i];
			else usage("Unknown parameter '" + args[i] + "'");
		}

		// Check: Do we have all required parameters?
		if (genomeVer.isEmpty()) usage("Missing genomer_version parameter");
		if (cdsFile.isEmpty()) usage("Missing cds_file parameter");
	}

	/**
	 * Read a file that has all CDS
	 */
	void readCdsFile() {
		cdsByTrId = new HashMap<>();

		if (cdsFile.endsWith("txt") || cdsFile.endsWith("txt.gz")) readCdsFileTxt();
		else readCdsFileFasta();

		if (cdsByTrId.isEmpty()) fatalError("CDS file is empty!");
	}

	/**
	 * Read CDSs from a file
	 * Format: Tab-separated format, containing "sequence \t transcriptId"
	 */
	void readCdsFileFasta() {
		// Load file
		FastaFileIterator ffi = new FastaFileIterator(cdsFile);
		for (String seq : ffi) {
			String trId = ffi.getName();
			add(trId, seq, ffi.getLineNum(), true);

			// Also try processing header line using different separators
			List<String> ids = ffi.fastaHeader2Ids();
			for (String id : ids) {
				// We don't check for uniqueness here since many items in this
				// list are tokens that are expected to be repeated
				add(id, seq, ffi.getLineNum(), false);
			}
		}
	}

	/**
	 * Read CDSs from a file
	 * Format: Tab-separated format, containing "sequence \t transcriptId"
	 */
	void readCdsFileTxt() {
		// Load file
		String cdsData = Gpr.readFile(cdsFile);
		String cdsLines[] = cdsData.split("\n");

		// Parse each line
		int lineNum = 1;
		for (String cdsLine : cdsLines) {
			// Split tab separated fields
			String field[] = cdsLine.split("\\s+");

			// Parse fields
			if (field.length >= 2) {
				// OK Parse fields
				String seq = field[1].trim();
				String trId = field[0].trim();

				// Repeated transcript Id? => Check that CDS is the same
				if ((cdsByTrId.get(trId) != null) && (!cdsByTrId.get(trId).equals(seq))) System.err.println("ERROR: Different CDS for the same transcript ID. This should never happen!!!\n\tLine number: " + lineNum + "\n\tTranscript ID:\t" + trId + "\n\tCDS:\t\t" + cdsByTrId.get(trId) + "\n\tCDS (new):\t" + seq);

				cdsByTrId.put(trId, seq); // Add it to the hash
			}

			lineNum++;
		}
	}

	/**
	 * Run command
	 */
	@Override
	public boolean run() {
		if (verbose) Timer.showStdErr("Checking database using CDS sequences");

		// Load config
		if (config == null) loadConfig();

		// Read CDS form file
		if (verbose) Timer.showStdErr("Reading CDSs from file '" + cdsFile + "'...");
		readCdsFile(); // Load CDS
		if (verbose) Timer.showStdErr("done (" + cdsByTrId.size() + " CDSs).");

		// Load predictor
		if (config.getSnpEffectPredictor() == null) loadDb();

		// Compare CDS
		if (verbose) Timer.showStdErr("Comparing CDS...");
		cdsCompare();
		if (verbose) Timer.showStdErr("done");

		return true;
	}

	/**
	 * Show same Transcript IDs
	 */
	String sampleTrIds(int maxTrIds) {
		StringBuilder sb = new StringBuilder();
		int count = 0;
		for (Gene g : config.getGenome().getGenes())
			for (Transcript tr : g) {
				sb.append("\t'" + tr.getId() + "'\n");
				if (count++ > maxTrIds) return sb.toString();
			}
		return sb.toString();
	}

	/**
	 * Show same Transcript IDs from FASTA file
	 */
	String sampleTrIdsFasta(int maxTrIds) {
		StringBuilder sb = new StringBuilder();
		int count = 0;
		for (String trid : cdsByTrId.keySet()) {
			sb.append("\t'" + trid + "'\n");
			if (count++ > maxTrIds) return sb.toString();
		}
		return sb.toString();
	}

	public void setCheckNumOk(boolean checkNumOk) {
		this.checkNumOk = checkNumOk;
	}

	public void setStoreAlignments(boolean storeAlignments) {
		this.storeAlignments = storeAlignments;
	}

	/**
	 * Show usage and exit
	 */
	@Override
	public void usage(String message) {
		if (message != null) System.err.println("Error: " + message + "\n");
		System.err.println("snpEff version " + SnpEff.VERSION);
		System.err.println("Usage: snpEff cds [options] genome_version cds_file");
		System.err.println("\nOptions:");
		System.err.println("\t-c , -config <file> : Specify config file");
		System.err.println("\t-noLog              : Do not report usage statistics to server");
		System.err.println("\t-v                  : Verbose mode");
		System.err.println("\t-d                  : Debug mode");
		System.exit(-1);
	}
}
