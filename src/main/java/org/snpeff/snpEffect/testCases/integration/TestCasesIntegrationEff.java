package org.snpeff.snpEffect.testCases.integration;

import java.util.List;

import org.junit.Test;
import org.snpeff.snpEffect.VariantEffect.EffectImpact;
import org.snpeff.util.Gpr;
import org.snpeff.vcf.VcfEffect;
import org.snpeff.vcf.VcfEntry;

import junit.framework.Assert;

/**
 *
 * Test cases for other 'effect' issues
 *
 * @author pcingola
 */
public class TestCasesIntegrationEff extends TestCasesIntegrationBase {

	public TestCasesIntegrationEff() {
		super();
	}

	/**
	 * Test output order
	 */
	@Test
	public void test_01() {
		Gpr.debug("Test");
		List<VcfEntry> vcfEntries = snpEffect("testHg3770Chr22", "tests/integration/eff/eff_sort.vcf", null);

		for (VcfEntry ve : vcfEntries) {
			if (verbose) System.out.println(ve);

			EffectImpact impPrev = EffectImpact.HIGH;
			for (VcfEffect veff : ve.getVcfEffects()) {
				EffectImpact imp = veff.getImpact();

				if (verbose) System.out.println("\t" + imp + "\t" + impPrev + "\t" + imp.compareTo(impPrev) + "\t" + veff);
				Assert.assertTrue(impPrev.compareTo(imp) <= 0); // Higher impact go first
				impPrev = imp;
			}
		}
	}

	/**
	 * Test output order: Canonical first
	 */
	@Test
	public void test_01_canonical() {
		Gpr.debug("Test");
		List<VcfEntry> vcfEntries = snpEffect("testHg3775Chr8", "tests/integration/eff/eff_sort_canon.vcf", null);

		// Only one entry in this file
		Assert.assertEquals(1, vcfEntries.size());

		VcfEntry ve = vcfEntries.get(0);
		VcfEffect veff = ve.getVcfEffects().get(0);

		Assert.assertEquals("ENST00000456015", veff.getTranscriptId());
	}

	/**
	 * Test GATK option: At most one effect per VCF entry
	 */
	@Test
	public void test_02() {
		Gpr.debug("Test");
		String args[] = { "-gatk" };
		List<VcfEntry> vcfEntries = snpEffect("testHg3770Chr22", "tests/integration/eff/eff_sort.vcf", args);

		for (VcfEntry ve : vcfEntries) {
			int numEffs = ve.getVcfEffects().size();
			if (verbose) System.out.println("Num effects:" + numEffs + "\t" + ve);
			Assert.assertTrue(numEffs <= 1);
		}
	}

	/**
	 * Make sure that empty VCF does not trigger an exception when creating the summary
	 */
	@Test
	public void test_03_EmptyVcf() {
		Gpr.debug("Test");
		String args[] = { "eff", "-noLog" };
		snpEffect("testHg3770Chr22", "tests/integration/eff/empty_only_header.vcf", args);
	}

	/**
	 * Test that CSV summary does not throw any error
	 */
	@Test
	public void test_04() {
		Gpr.debug("Test");
		String args[] = { "-csvStats", "test_04_TestCasesEff.csv" };
		snpEffect("testHg3770Chr22", "tests/integration/eff/eff_sort.vcf", args);
	}

	/**
	 * GATK mode should not have SPLICE_REGION (it is currently not supported)
	 */
	@Test
	public void test_05() {
		Gpr.debug("Test");
		String genomeName = "testHg3775Chr1";
		String vcf = "tests/integration/eff/gatk_NO_splice_regions.vcf";
		String args[] = { "eff", "-noLog", "-gatk" };
		List<VcfEntry> vcfEntries = snpEffect(genomeName, vcf, args);

		for (VcfEntry ve : vcfEntries) {
			if (verbose) System.out.println(ve);

			for (VcfEffect veff : ve.getVcfEffects()) {
				if (verbose) System.out.println("\t'" + veff.getEffectsStr() + "'\t" + veff);
				if (veff.getEffectsStr().indexOf("SPLICE_SITE_REGION") >= 0) throw new RuntimeException("Splice region effects should not present in GATK compatible mode");
			}
		}
	}

	/**
	 * Test an MNP at the end of the transcript: We should be able to annotate without throwing any error
	 */
	@Test
	public void test_06() {
		Gpr.debug("Test");
		String args[] = {};
		List<VcfEntry> list = snpEffect("testHg3775Chr15", "tests/integration/eff/mnp_insertion_at_transcript_end.vcf", args);

		// We should be able to annotate this entry (if INFO is empty, something went wrong)
		Assert.assertFalse(list.get(0).getInfoStr().isEmpty());
	}

	/**
	 * Test an MNP at the end of the transcript: We should be able to annotate without throwing any error
	 */
	@Test
	public void test_07() {
		Gpr.debug("Test");
		String args[] = {};
		List<VcfEntry> list = snpEffect("testHg3775Chr10", "tests/integration/eff/mnp_deletion.vcf", args);

		// We should be able to annotate this entry (if INFO is empty, something went wrong)
		Assert.assertFalse(list.get(0).getInfoStr().isEmpty());
	}

	/**
	 * Fixing bug: GATK does not annotate all VCF entries
	 */
	@Test
	public void test_08_gatk_missing_annotations() {
		Gpr.debug("Test");

		String genomeName = "testMycobacterium_tuberculosis_CCDC5079_uid203790";
		String vcf = "tests/integration/eff/test_gatk_no_annotations.vcf";
		String args[] = { "-noLog", "-gatk" };
		List<VcfEntry> vcfEntries = snpEffect(genomeName, vcf, args);

		for (VcfEntry ve : vcfEntries) {
			int count = 0;
			if (verbose) System.out.println(ve);

			for (VcfEffect veff : ve.getVcfEffects()) {
				if (verbose) System.out.println("\t'" + veff.getEffectsStr() + "'\t" + veff);
				count++;
			}

			// Check that there is one and only one annotation
			Assert.assertEquals(1, count);
		}
	}

}
