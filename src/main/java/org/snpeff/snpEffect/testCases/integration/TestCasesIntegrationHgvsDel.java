package org.snpeff.snpEffect.testCases.integration;

import org.junit.Test;
import org.snpeff.util.Gpr;

/**
 * Test cases for HGVS notation on insertions
 */
public class TestCasesIntegrationHgvsDel extends TestCasesIntegrationBase {

	public TestCasesIntegrationHgvsDel() {
		super();
	}

	/**
	 * This frameshift caused an exception while processing HGVS protein notation
	 */
	@Test
	public void test_01_hgvs_deletions_chr11() {
		Gpr.debug("Test");

		String genomeName = "testHg19Chr11";
		String vcf = "tests/integration/hgvsDel/test_01_hgvs_deletions_chr11.vcf";

		snpEffect(genomeName, vcf, null);

	}

}
