package org.snpeff.snpEffect.testCases.integration;

import org.junit.Test;
import org.snpeff.SnpEff;
import org.snpeff.snpEffect.commandLine.SnpEffCmdEff;
import org.snpeff.util.Gpr;

import junit.framework.Assert;

/**
 * Calculate missense over silent ratio
 *
 * @author pcingola
 */
public class TestCasesIntegrationMissenseSilentRatio {

	boolean verbose = false;

	public TestCasesIntegrationMissenseSilentRatio() {
		super();
	}

	@Test
	public void test_01() {
		Gpr.debug("Test");
		String args[] = { //
				"-classic" //
				, "-useLocalTemplate" //
				, "testHg3765Chr22" //
				, "./tests/integration/missenseSilentRatio/missenseSilent.chr22.vcf.gz" //
		};

		SnpEff cmd = new SnpEff(args);
		cmd.setVerbose(verbose);
		cmd.setSupressOutput(!verbose);
		SnpEffCmdEff snpeff = (SnpEffCmdEff) cmd.cmd();

		snpeff.run();

		double silentRatio = snpeff.getChangeEffectResutStats().getSilentRatio();
		if (verbose) System.err.println("Missense / Silent ratio: " + silentRatio);

		Assert.assertEquals(1.19, silentRatio, 0.1);
	}
}
