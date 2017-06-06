package org.snpeff.snpEffect.testCases.integration;

import org.junit.Test;
import org.snpeff.util.Gpr;

/**
 *
 * Test case
 */
public class TestCasesIntegrationHgvsHard extends TestCasesIntegrationBase {

	public TestCasesIntegrationHgvsHard() {
		super();
		shiftHgvs = true;
	}

	/**
	 * Test some 'dup' form chr17
	 */
	@Test
	public void test_hgvs_dup() {
		Gpr.debug("Test");

		String genome = "testHg19Chr17";
		String vcf = "tests/integration/hgvsHard/hgvs_dup.vcf";

		compareHgvs(genome, vcf);
	}

	@Test
	public void test_hgvs_md_chr17() {
		Gpr.debug("Test");
		String genome = "testHg19Chr17";
		String vcf = "tests/integration/hgvsHard/hgvs_md.chr17.vcf";
		compareHgvs(genome, vcf, false);
	}

	@Test
	public void test_hgvs_walk_and_roll_1() {
		Gpr.debug("Test");

		String genome = "testHg19Chr1";
		String vcf = "tests/integration/hgvsHard/hgvs_jeremy_1.vcf";

		compareHgvs(genome, vcf);
	}

	@Test
	public void test_hgvs_walk_and_roll_2() {
		Gpr.debug("Test");

		String genome = "testHg19Chr17";
		String vcf = "tests/integration/hgvsHard/hgvs_walk_and_roll.1.vcf";

		compareHgvs(genome, vcf, true);
	}

	@Test
	public void test_hgvs_walk_and_roll_3() {
		Gpr.debug("Test");
		String genome = "testHg19Chr13";
		String vcf = "tests/integration/hgvsHard/hgvs_savant.vcf";
		compareHgvs(genome, vcf, true);
	}

	/**
	 * HGVS (protein) was reporting a wrong deletion length
	 */
	@Test
	public void test_HgvsP_deleteion_length() {
		Gpr.debug("Test");
		String genome = "testHg19Chr1";
		String vcf = "tests/integration/hgvsHard/hgvs_protein_deleteion_length.vcf";
		compareHgvs(genome, vcf);
	}

	/**
	 * Some difficult HGVS test case
	 */
	@Test
	public void test_HgvsP_deleteion_length_2() {
		Gpr.debug("Test");
		String genome = "testHg19Chr2";
		String vcf = "tests/integration/hgvsHard/hgvs_protein_deleteion_length_2.vcf";
		compareHgvs(genome, vcf);
	}
}
