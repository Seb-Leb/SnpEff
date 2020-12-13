package org.snpeff.interval;

import java.util.HashMap;
import java.util.concurrent.ArrayBlockingQueue;

import org.snpeff.interval.Exon.ExonSpliceType;
import org.snpeff.snpEffect.Config;
import org.snpeff.snpEffect.SnpEffectPredictor;
import org.snpeff.stats.CountByType;
import org.snpeff.util.Gpr;
import org.snpeff.util.Timer;

/**
 * Characterize exons based on alternative splicing
 *
 * References: "Alternative splicing and evolution - diversification, exon definition and function"  (see Box 1)
 *
 * @author pablocingolani
 */
public class ExonSpliceCharacterizer {

	public static final int MAX_EXONS = 1000; // Do not characterize transcripts having more than this number of exons
	public static final int SHOW_EVERY = 100;

	private volatile int threadsCompleted;
	private WorkerThread[] workers;
	private ArrayBlockingQueue<Runnable> taskQueue;
	private int taskCount;
	private volatile int tasksCompleted;
	private volatile boolean jobInProgress;

	boolean verbose = true;
	Genome genome;
	HashMap<Exon, Exon.ExonSpliceType> typeByExon;
	CountByType countByType = new CountByType();

	public ExonSpliceCharacterizer(Genome genome) {
		this.genome = genome;
		typeByExon = new HashMap<Exon, Exon.ExonSpliceType>();
	}

	public ExonSpliceCharacterizer(String genomeVer) {
		Config config = new Config(genomeVer);
		SnpEffectPredictor snpEffectPredictor = config.loadSnpEffectPredictor();
		genome = snpEffectPredictor.getGenome();
		typeByExon = new HashMap<Exon, Exon.ExonSpliceType>();
	}

	/**
	 * Characterize all exons
	 */
	public CountByType characterize() {
		type();
		return countByType;
	}

	/**
	 * Count number of exons
	 */
	int countExons() {
		int count = 0;
		for (Gene g : genome.getGenes())
			for (Transcript tr : g)
				count += tr.numChilds();
		return count;
	}

	/**
	 * Does the marker intersect any exon in 'tr'?
	 * @param m
	 * @param tr
	 * @return
	 */
	boolean intersectsAnyExon(Marker m, Transcript tr) {
		for (Exon e : tr)
			if (m.intersects(e)) return true;
		return false;
	}

	/**
	 * Is thins an ALTERNATIVE_3SS exon?
	 * @param exon
	 * @param gene
	 * @return
	 */
	boolean isAlt3ss(Exon exon, Gene gene) {
		for (Transcript tr : gene)
			for (Exon e : tr) {
				if (exon.intersects(e)) {
					if (exon.isStrandPlus()) {
						// Same exon end, different exon start?
						if ((exon.getStart() != e.getStart()) && (exon.getEnd() == e.getEnd())) return true;
					} else {
						// Same exon end, different exon start? (negative strand)
						if ((exon.getStart() == e.getStart()) && (exon.getEnd() != e.getEnd())) return true;
					}
				}
			}

		return false;
	}

	/**
	 * Is thins an ALTERNATIVE_5SS exon?
	 * @param exon
	 * @param gene
	 * @return
	 */
	boolean isAlt5ss(Exon exon, Gene gene) {
		for (Transcript tr : gene)
			for (Exon e : tr) {
				if (exon.intersects(e)) {
					if (exon.isStrandPlus()) {
						// Same exon start, different exon end?
						if ((exon.getStart() == e.getStart()) && (exon.getEnd() != e.getEnd())) return true;
					} else {
						// Same exon start, different exon end? (negative strand)
						if ((exon.getStart() != e.getStart()) && (exon.getEnd() == e.getEnd())) return true;
					}
				}
			}

		return false;
	}

	/**
	 * Is this exon mutually exclusive with another exon?
	 * @param exon
	 * @param gene
	 * @return
	 */
	boolean isMutEx(Exon exon, Gene gene) {
		if (gene.numChilds() <= 1) return false;

		//---
		// Make a list of all unique exons
		//---
		String exonKey = key(exon);
		HashMap<String, Exon> uniqEx = new HashMap<String, Exon>();
		for (Transcript tr : gene)
			for (Exon e : tr) {
				String ekey = key(e);
				if (!exonKey.equals(ekey)) uniqEx.put(ekey, e);
			}

		//---
		// For each unique exon, compare if it is mutually exclusive with 'exon'
		//---
		Transcript exonTr = (Transcript) exon.getParent();
		for (Exon e : uniqEx.values()) {
			ExonSpliceType type = typeByExon.get(e);;
			if (type == Exon.ExonSpliceType.SKIPPED) { // Only check for these exons (avoid all ALT*)
				boolean xor = true;
				for (Transcript tr : gene) {
					if (exonTr.intersects(tr) && exon.intersects(tr) && e.intersects(exonTr)) { // Make sure both exons intersect both transcripts (otherwise they cannot be mutually exclusive)
						xor &= intersectsAnyExon(e, tr) ^ intersectsAnyExon(exon, tr);
					} else xor = false;
				}

				// XOR is true? => Mutually exclusive
				if (xor) return true;
			}
		}

		return false;
	}

	/**
	 * Create a simple hash key based on choromosomal position
	 * @param m
	 * @return
	 */
	String key(Marker m) {
		return m.getChromosomeName() + ":" + m.getStart() + "-" + m.getEnd();
	}

	public void setVerbose(boolean verbose) {
		this.verbose = verbose;
	}

	
	/**
	 * Mark exons types
	 */
	
	void type() {
		if (verbose) {
			Timer.showStdErr("Caracterizing exons by splicing (stage 1) : ");
			System.out.print("\t");
		}

		jobInProgress = true;
		taskCount = 0;
		tasksCompleted = 0;

		// create the queue and workers
		taskQueue = new ArrayBlockingQueue<Runnable>(10);
		int processors = Runtime.getRuntime().availableProcessors();
		workers = new WorkerThread[processors];
		for (int i = 0; i < processors; i++) {
			workers[i] = new WorkerThread();
		}

		// wait for the worker to start...
		try {
			Thread.sleep(5000);
		}
		catch(InterruptedException e){
			e.printStackTrace();
		}

		// Find retained exons
		int numTr = 1;
		for (Gene g : genome.getGenes()) {
			// Count exons
			CountByType count = new CountByType();
			for (Transcript tr : g)
				for (Exon e : tr)
					count.inc(key(e));

			// Label exons
			int countTr = g.numChilds();

			Timer.showStdErr("number of trxps: ");
			System.out.print(countTr);
			System.out.print("\t");

			for (Transcript tr : g) {
				if (verbose) Gpr.showMark(numTr++, SHOW_EVERY, "\t");
				labelExonsTask task = new labelExonsTask(g, tr, countTr, count);
				try {
					taskCount++;
					taskQueue.put(task);
					Thread.sleep(10);
				}
				catch(InterruptedException e){
					e.printStackTrace();
				}
			}
		}
		try {
			while (jobInProgress) {Thread.sleep(1000);}
		}
		catch(InterruptedException e){
		}

		// Now analyze if there are mutually exclusive exons
		if (verbose) {
			System.err.println("");
			Timer.showStdErr("Caracterizing exons by splicing (stage 2) : ");
			System.out.print("\t");
		}

		jobInProgress = true;
		taskCount = 0;
		tasksCompleted = 0;
		numTr = 1;
		for (Gene g : genome.getGenes()) {
			for (Transcript tr : g) {
					if (verbose) Gpr.showMark(numTr++, SHOW_EVERY, "\t");
					exonMutExTask task = new exonMutExTask(g, tr);
				try {
					taskCount++;
					taskQueue.put(task);
				}
				catch(InterruptedException e){
					e.printStackTrace();
				}
			}
		}
		try {
			while (jobInProgress) {Thread.sleep(1000);}
		}
		catch(InterruptedException e){
		}
		if (verbose) Timer.showStdErr("done.");
	}

	/**
	 * Mark this exons as 'type'
	 * @param e
	 * @param type
	 */
	synchronized void type(Exon e, Exon.ExonSpliceType type) {
		e.spliceType = type;
		countByType.inc(type.toString());
		typeByExon.put(e, type);
	}

	synchronized private void jobFinished() {
		jobInProgress = false;
	}

   synchronized private void taskFinished() {
      tasksCompleted++;
      if (tasksCompleted == taskCount) { // all threads have finished
         jobFinished();
      }
   }

   	private class labelExonsTask implements Runnable{
   		Gene g;
   		Transcript tr;
   		Exon e;
   		int countTr;
   		CountByType count;
   		labelExonsTask(Gene g, Transcript tr, int countTr, CountByType count){
   			this.g = g;
   			this.tr = tr;
   			this.countTr = countTr;
   			this.count = count;
   		}
   		public void run() {
			for (Exon e : tr) {
				String eKey = key(e);
				int countEx = (int) count.get(eKey);

				// Is this exon maintained in all transcripts?
				if (countEx == countTr) type(e, Exon.ExonSpliceType.RETAINED);
				else {
					if (isAlt3ss(e, g)) type(e, Exon.ExonSpliceType.ALTTENATIVE_3SS);
					else if (isAlt5ss(e, g)) type(e, Exon.ExonSpliceType.ALTTENATIVE_5SS);
					else if (tr.numChilds() > 1) {
						if (e.getRank() == 1) type(e, Exon.ExonSpliceType.ALTTENATIVE_PROMOMOTER);
						else if (e.getRank() == tr.numChilds()) type(e, Exon.ExonSpliceType.ALTTENATIVE_POLY_A);
						else type(e, Exon.ExonSpliceType.SKIPPED);
					}
				}
			}
			taskFinished();
   		}
   	}

   	private class exonMutExTask implements Runnable{
   		Gene g;
   		Transcript tr;
   		Exon e;
   		exonMutExTask(Gene g, Transcript tr){
   			this.g = g;
   			this.tr = tr;
   		}
   		public void run() {
			if (tr.numChilds() < MAX_EXONS) {
				for (Exon e : tr) {
					ExonSpliceType type = typeByExon.get(e);
					if (type == ExonSpliceType.SKIPPED) { // Try to re-annotate only these
						if (isMutEx(e, g)) type(e, Exon.ExonSpliceType.MUTUALLY_EXCLUSIVE);
					}
				}
			} else {
				System.err.println("");
				Gpr.debug("WARNING: Gene '" + g.getId() + "', transcript '" + tr.getId() + "' has too many exons (" + tr.numChilds() + " exons). Skipped");
			}
			taskFinished();
   		}
   	}

   private class WorkerThread extends Thread {
      WorkerThread() {
         try {
            setDaemon(true);
         }
         catch (Exception e) {
         }
         start();
      }
      public void run() {
         while (true) {
            try {
               Runnable task = taskQueue.take();
               task.run();
            }
            catch (InterruptedException e) {
            }
         }
      }
   }
}
