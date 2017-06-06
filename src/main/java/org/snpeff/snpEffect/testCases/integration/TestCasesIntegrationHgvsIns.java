package org.snpeff.snpEffect.testCases.integration;

import org.junit.Test;
import org.snpeff.util.Gpr;

/**
 * Test cases for HGVS notation on insertions
 */
public class TestCasesIntegrationHgvsIns extends TestCasesIntegrationBase {

	public TestCasesIntegrationHgvsIns() {
		super();
	}

	/**
	 * Insertion / duplication issues
	 */
	@Test
	public void test_02_hgvs_insertions_chr1() {
		Gpr.debug("Test");
		checkHgvs("testHg19Chr1", "tests/integration/hgvsIns/hgvs_ins_dups_chr1.vcf", 4);
	}

	/**
	 * Insertion / duplication issues
	 */
	@Test
	public void test_03_hgvs_insertions_chr3() {
		Gpr.debug("Test");
		checkHgvs("testHg19Chr3", "tests/integration/hgvsIns/hgvs_ins_dups_chr3.vcf", 2);
	}

	/**
	 * Insertion / duplication issues
	 */
	@Test
	public void test_04_hgvs_insertions_chr4() {
		Gpr.debug("Test");
		checkHgvs("testHg19Chr4", "tests/integration/hgvsIns/hgvs_ins_dups_chr4.vcf", 2);
	}

	/**
	 * Insertion / duplication issues
	 */
	@Test
	public void test_05_hgvs_insertions_chr19() {
		Gpr.debug("Test");
		checkHgvs("testHg19Chr19", "tests/integration/hgvsIns/hgvs_ins_dups_chr19.vcf", 2);
	}

	/**
	 * Insertion / duplication issues
	 */
	@Test
	public void test_06_hgvs_insertions_chr7() {
		Gpr.debug("Test");
		checkHgvs("testHg19Chr7", "tests/integration/hgvsIns/hgvs_ins_chr7.vcf", 2);
	}

	/**
	 * This frameshift caused an exception while processing HGVS protein notation
	 */
	@Test
	public void test_07_hgvs_insertions1() {
		Gpr.debug("Test");

		String genomeName = "testHg3775Chr1";
		String vcf = "tests/integration/hgvsIns/hgvs_ins_07.vcf";

		snpEffect(genomeName, vcf, null);
	}

}
