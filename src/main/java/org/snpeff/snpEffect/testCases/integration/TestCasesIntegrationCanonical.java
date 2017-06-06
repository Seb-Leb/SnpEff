package org.snpeff.snpEffect.testCases.integration;

import org.junit.Test;
import org.snpeff.SnpEff;
import org.snpeff.interval.Gene;
import org.snpeff.interval.Transcript;
import org.snpeff.snpEffect.commandLine.SnpEffCmdEff;
import org.snpeff.util.Gpr;

import junit.framework.Assert;

/**
 *
 * Test cases for canonical transcript selection
 *
 * @author pcingola
 */
public class TestCasesIntegrationCanonical extends TestCasesIntegrationBase {

	public TestCasesIntegrationCanonical() {
		super();
	}

	/**
	 * Test canonical transcripts
	 */
	@Test
	public void test_01() {
		Gpr.debug("Test");
		String geneId = "APOBEC3H";
		String trId = "NM_001166003.2";
		String args[] = { "-canon", "testHg19Chr22", "tests/integration/canonical/empty.vcf" };

		SnpEff cmd = new SnpEff(args);
		SnpEffCmdEff cmdEff = (SnpEffCmdEff) cmd.cmd();
		cmdEff.setVerbose(verbose);
		cmdEff.setSupressOutput(!verbose);
		cmdEff.load();

		Gene gene = cmdEff.getConfig().getSnpEffectPredictor().getGene(geneId);
		Transcript tr = gene.subIntervals().iterator().next();
		Assert.assertEquals("Expecting transcript ID does not match", trId, tr.getId());
	}

	/**
	 * Test Somatic vs Germline
	 */
	@Test
	public void test_02() {
		Gpr.debug("Test");
		String geneId = "APOBEC3H";
		String trId = "NM_181773.4";
		String args[] = { "-canonList", "tests/integration/canonical/canon_geneId2trId_test02.txt", "testHg19Chr22", "tests/integration/canonical/empty.vcf" };

		SnpEff cmd = new SnpEff(args);
		SnpEffCmdEff cmdEff = (SnpEffCmdEff) cmd.cmd();
		cmdEff.setVerbose(verbose);
		cmdEff.setSupressOutput(!verbose);
		cmdEff.load();

		Gene gene = cmdEff.getConfig().getSnpEffectPredictor().getGene(geneId);
		Transcript tr = gene.subIntervals().iterator().next();
		Assert.assertEquals("Expecting transcript ID does not match", trId, tr.getId());
	}

}
