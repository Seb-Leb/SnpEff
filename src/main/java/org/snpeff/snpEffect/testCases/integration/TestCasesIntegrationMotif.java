package org.snpeff.snpEffect.testCases.integration;

import org.junit.Test;
import org.snpeff.SnpEff;
import org.snpeff.snpEffect.VariantEffect.EffectImpact;
import org.snpeff.util.Gpr;

/**
 * Test Motif databases
 *
 * @author pcingola
 */
public class TestCasesIntegrationMotif extends TestCasesIntegrationBase {

	public TestCasesIntegrationMotif() {
		super();
	}

	@Test
	public void test_01() {
		Gpr.debug("Test");
		checkMotif("testHg3770Chr22", "tests/integration/motif/test_motif_01.vcf", "MA0099.2:AP1", EffectImpact.LOW, false);
	}

	@Test
	public void test_01_ann() {
		Gpr.debug("Test");
		checkMotif("testHg3770Chr22", "tests/integration/motif/test_motif_01.vcf", "MA0099.2", EffectImpact.LOW, true);
	}

	@Test
	public void test_02() {
		Gpr.debug("Test");
		checkMotif("testHg3770Chr22", "tests/integration/motif/test_motif_02.vcf", "MA0099.2:AP1", EffectImpact.MODIFIER, false);
	}

	@Test
	public void test_02_ann() {
		Gpr.debug("Test");
		checkMotif("testHg3770Chr22", "tests/integration/motif/test_motif_02.vcf", "MA0099.2", EffectImpact.MODIFIER, true);
	}

	@Test
	public void test_03() {
		Gpr.debug("Test");
		checkMotif("testHg3770Chr22", "tests/integration/motif/test_motif_03.vcf", "MA0099.2:AP1", EffectImpact.LOW, false);
	}

	@Test
	public void test_03_ann() {
		Gpr.debug("Test");
		checkMotif("testHg3770Chr22", "tests/integration/motif/test_motif_03.vcf", "MA0099.2", EffectImpact.LOW, true);
	}

	/**
	 * MNP outside Motif: Should not throw any exception
	 */
	@Test
	public void test_04() {
		Gpr.debug("Test");
		String genome = "testHg3775Chr11";
		String vcf = "tests/integration/motif/craig_chr11.vcf";

		String args[] = { "-noLog", genome, vcf };
		SnpEff snpEff = new SnpEff(args);
		snpEff.setVerbose(verbose);
		snpEff.setSupressOutput(!verbose);
		snpEff.setDebug(debug);
		snpEff.run();
	}

	/**
	 * Motif has 9 bases but ENSEMBL file marks it as a 10 base interval
	 * SNP affect last base (as marked by ENSEMBL), since there is no sequence for that base position, an exception is thrown.
	 */
	@Test
	public void test_05() {
		Gpr.debug("Test");
		String genome = "testHg3775Chr14";
		String vcf = "tests/integration/motif/craig_chr14.vcf";

		String args[] = { "-noLog", genome, vcf };
		SnpEff snpEff = new SnpEff(args);
		snpEff.setVerbose(verbose);
		snpEff.setSupressOutput(!verbose);
		snpEff.setDebug(debug);
		snpEff.run();
	}

}
