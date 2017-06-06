package org.snpeff.snpEffect.testCases.integration;

import java.util.List;

import org.junit.Test;
import org.snpeff.SnpEff;
import org.snpeff.interval.Gene;
import org.snpeff.interval.Transcript;
import org.snpeff.snpEffect.Config;
import org.snpeff.snpEffect.EffectType;
import org.snpeff.snpEffect.SnpEffectPredictor;
import org.snpeff.snpEffect.commandLine.SnpEffCmdEff;
import org.snpeff.util.Gpr;
import org.snpeff.vcf.VcfEffect;
import org.snpeff.vcf.VcfEntry;

import junit.framework.Assert;

/**
 * Test case for exon frames
 *
 * @author pcingola
 */
public class TestCasesIntegrationExonFrame {

	boolean verbose = false;

	public TestCasesIntegrationExonFrame() {
		super();
	}

	/**
	 * Test database: Build, check and annotate
	 */
	@Test
	public void test_01() {
		Gpr.debug("Test");

		//---
		// Build database
		//---
		String genomeName = "testLukas";
		String args[] = { "build", "-noLog", "-gff3", genomeName };

		SnpEff snpEff = new SnpEff(args);
		snpEff.setVerbose(verbose);
		snpEff.setSupressOutput(!verbose);
		boolean ok = snpEff.run();
		Assert.assertTrue(ok);

		//---
		// Load database and check some numbers
		//---
		String configFile = Config.DEFAULT_CONFIG_FILE;
		Config config = new Config(genomeName, configFile);
		if (verbose) System.out.println("Loading database");
		SnpEffectPredictor snpEffectPredictor = config.loadSnpEffectPredictor();

		// Find transcript (there is only one)
		Transcript transcript = null;
		for (Gene gene : snpEffectPredictor.getGenome().getGenes())
			for (Transcript tr : gene)
				transcript = tr;

		// Check parameters
		Assert.assertEquals(454126, transcript.getCdsStart());
		Assert.assertEquals(450599, transcript.getCdsEnd());

		//---
		// Check annotations
		//---
		String vcfFileName = "tests/integration/exonFrame/testLukas.vcf";
		String argsEff[] = { "-classic", "-noHgvs", "-ud", "0", genomeName, vcfFileName };

		// Annotate
		SnpEff cmd = new SnpEff(argsEff);
		SnpEffCmdEff cmdEff = (SnpEffCmdEff) cmd.cmd();
		cmdEff.setVerbose(verbose);
		cmdEff.setSupressOutput(!verbose);
		List<VcfEntry> vcfEntries = cmdEff.run(true);
		Assert.assertTrue("Errors while executing SnpEff", cmdEff.getTotalErrs() <= 0);

		// Analyze annotations
		for (VcfEntry ve : vcfEntries) {
			if (verbose) System.out.println(ve.toStringNoGt());

			EffectType expectedEffect = EffectType.valueOf(ve.getInfo("EXP_EFF"));
			String expectedAa = ve.getInfo("EXP_AA");
			String expectedCodon = ve.getInfo("EXP_CODON");

			boolean found = false;
			for (VcfEffect veff : ve.getVcfEffects()) {
				String eff = veff.getEffectType().toString();

				if (verbose) {
					System.out.println("\t" + veff);
					System.out.println("\t\tExpecing: '" + expectedEffect + "'\tFound: '" + eff + "'");
					System.out.println("\t\tExpecing: '" + expectedAa + "'\tFound: '" + veff.getAa() + "'");
					System.out.println("\t\tExpecing: '" + expectedCodon + "'\tFound: '" + veff.getCodon() + "'");
				}

				// Effect matches expected?
				if (veff.hasEffectType(expectedEffect) //
						&& ((veff.getAa() == null) || veff.getAa().isEmpty() || expectedAa.equals(veff.getAa())) //
						&& ((veff.getCodon() == null) || veff.getCodon().isEmpty() || expectedCodon.equals(veff.getCodon())) //
				) //
					found = true;
			}

			if (!found) throw new RuntimeException("Cannot find expected effect '" + expectedEffect + "', amino acid change '" + expectedAa + "' and codon change '" + expectedCodon + "'");
		}
	}

	/**
	 * Build genome (no exceptions should be thrown)
	 */
	@Test
	public void test_02() {
		Gpr.debug("Test");

		// Build database
		String genomeName = "testMacuminata";
		String args[] = { "build", "-noLog", genomeName };

		SnpEff snpEff = new SnpEff(args);
		snpEff.setVerbose(verbose);
		snpEff.setSupressOutput(!verbose);
		boolean ok = snpEff.run();
		Assert.assertTrue(ok);
	}
}
