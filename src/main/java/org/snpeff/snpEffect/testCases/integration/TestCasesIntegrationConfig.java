package org.snpeff.snpEffect.testCases.integration;

import org.junit.Test;
import org.snpeff.SnpEff;
import org.snpeff.snpEffect.Config;
import org.snpeff.util.Gpr;

import junit.framework.Assert;

/**
 * Test case
 */
public class TestCasesIntegrationConfig {

	boolean debug = false;
	boolean verbose = false;

	/**
	 * Check that config file can be overriden by command line options
	 */
	@Test
	public void test_01_ConfigOverride() {
		Gpr.debug("Test");

		// Create command
		String repo = "http://nonsense.url/test/zzz";
		String args[] = { //
				"-configOption" //
				, Config.KEY_DATABASE_REPOSITORY + "=" + repo //
				, "testHg3775Chr22" //
				, "tests/integration/config/test_ann_01.vcf" //
		};

		// Create command and run
		SnpEff cmd = new SnpEff(args);
		cmd.setSupressOutput(!verbose);
		cmd.setVerbose(verbose);
		cmd.setDebug(debug);
		cmd.run();

		// Check that config option really changed
		if (verbose) System.out.println("Repository: " + cmd.getConfig().getDatabaseRepository());
		Assert.assertEquals(repo, cmd.getConfig().getDatabaseRepository());
	}

}
