package org.snpeff.snpEffect.testCases.integration;

import java.util.List;

import org.junit.Assert;
import org.junit.Test;
import org.snpeff.SnpEff;
import org.snpeff.snpEffect.commandLine.SnpEffCmdEff;
import org.snpeff.util.Gpr;
import org.snpeff.vcf.VcfEffect;
import org.snpeff.vcf.VcfEntry;

/**
 * Test random SNP changes
 *
 * @author pcingola
 */
public class TestCasesIntegrationInsEtc {

	public static boolean debug = false;
	public static boolean verbose = false || debug;

	/**
	 * Insertion on minus strand
	 */
	@Test
	public void test_01_InsOffByOne() {
		Gpr.debug("Test");
		String args[] = { "-classic", "-noHgvs", "testENST00000268124", "tests/integration/insEtc/ins_off_by_one.vcf" };

		SnpEff cmd = new SnpEff(args);
		SnpEffCmdEff snpeff = (SnpEffCmdEff) cmd.cmd();
		snpeff.setSupressOutput(!verbose);
		snpeff.setVerbose(verbose);

		List<VcfEntry> vcfEnties = snpeff.run(true);
		int count = 0;
		for (VcfEntry ve : vcfEnties) {

			// Get first effect (there should be only one)
			List<VcfEffect> veffs = ve.getVcfEffects();
			VcfEffect veff = veffs.get(0);

			Assert.assertEquals("Q53QQ", veff.getAa());
			count++;
		}

		Assert.assertTrue(count > 0);
	}

}
