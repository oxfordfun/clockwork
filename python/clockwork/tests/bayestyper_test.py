import filecmp
import os
import shutil
import unittest

from clockwork import bayestyper

modules_dir = os.path.dirname(os.path.abspath(bayestyper.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data', 'bayestyper')


def vcfs_the_same(vcf1, vcf2):
    '''Checks VCFs are the "same". Ignores the header line added
    by bcftools ##bcftools_normCommand=norm ... because the paths to files
    in that line are unknown when we run the tests'''
    with open(vcf1) as f:
        lines1 = [x for x in f if not x.startswith('##bcftools_normCommand=norm')]
    with open(vcf2) as f:
        lines2 = [x for x in f if not x.startswith('##bcftools_normCommand=norm')]

    return lines1 == lines2


class TestBayestyper(unittest.TestCase):
    def test_normalise_vcf(self):
        '''test normalise_vcf'''
        vcf_in = os.path.join(data_dir, 'normalise_vcf.in.1.vcf')
        expect_vcf = os.path.join(data_dir, 'normalise_vcf.out.1.vcf')
        ref_fa = os.path.join(data_dir, 'normalise_vcf.ref.fa')
        tmp_out = 'tmp.normalise_vcf.vcf'
        bayestyper.normalise_vcf(vcf_in, ref_fa, tmp_out, 'bcftools')
        self.assertTrue(vcfs_the_same(expect_vcf, tmp_out))
        os.unlink(tmp_out)


    def test_normalise_vcfs(self):
        '''test normalise_vcfs'''
        vcf_files = {
            'tool1': os.path.join(data_dir, 'normalise_vcf.in.1.vcf'),
            'tool2': os.path.join(data_dir, 'normalise_vcf.in.2.vcf'),
        }
        ref_fa = os.path.join(data_dir, 'normalise_vcf.ref.fa')
        outprefix = 'tmp.normalise_vcfs'
        normalised_vcfs = bayestyper.normalise_vcfs(vcf_files, ref_fa, 'bcftools', outprefix)
        self.assertTrue(vcfs_the_same(normalised_vcfs['tool1'], os.path.join(data_dir, 'normalise_vcf.out.1.vcf')))
        self.assertTrue(vcfs_the_same(normalised_vcfs['tool2'], os.path.join(data_dir, 'normalise_vcf.out.2.vcf')))
        os.unlink(normalised_vcfs['tool1'])
        os.unlink(normalised_vcfs['tool2'])


    def test_vcf_add_contigs_to_header_and_get_sample(self):
        '''test vcf_add_contigs_to_header_and_get_sample'''
        vcf_in = os.path.join(data_dir, 'vcf_add_contigs_to_header_and_get_sample.in.vcf')
        tmp_out = 'tmp.vcf_add_contigs_to_header_and_get_sample.out.vcf'
        ref_lengths = {'ref.1': 1000, 'ref.2': 42, 'ref.4': 100}
        got_sample = bayestyper.vcf_add_contigs_to_header_and_get_sample(vcf_in, tmp_out, ref_lengths)
        self.assertEqual('samplemcsampleface', got_sample)
        expect_vcf = os.path.join(data_dir, 'vcf_add_contigs_to_header_and_get_sample.expect.vcf')
        self.assertTrue(filecmp.cmp(expect_vcf, tmp_out, shallow=False))
        os.unlink(tmp_out)


    def test_fix_vcf_headers_and_get_sample_name(self):
        '''test fix_vcf_headers_and_get_sample_name'''
        vcf_files = {
            'tool1': os.path.join(data_dir, 'fix_vcf_headers_and_get_sample_name.in.1.vcf'),
            'tool2': os.path.join(data_dir, 'fix_vcf_headers_and_get_sample_name.in.2.vcf'),
        }
        ref_lengths = {'ref.1': 1000, 'ref.2': 42, 'ref.4': 100}
        outprefix = 'tmp.fix_vcf_headers_and_get_sample_name'
        got_vcfs, got_sample = bayestyper.fix_vcf_headers_and_get_sample_name(vcf_files, outprefix, ref_lengths)
        self.assertEqual('samplemcsampleface', got_sample)
        self.assertTrue(vcfs_the_same(got_vcfs['tool1'], os.path.join(data_dir, 'fix_vcf_headers_and_get_sample_name.out.1.vcf')))
        self.assertTrue(vcfs_the_same(got_vcfs['tool2'], os.path.join(data_dir, 'fix_vcf_headers_and_get_sample_name.out.2.vcf')))
        os.unlink(got_vcfs['tool1'])
        os.unlink(got_vcfs['tool2'])


    def test_get_reads_type_option_for_kmc(self):
        '''test get_reads_type_option_for_kmc'''
        self.assertEqual('-fbam', bayestyper.get_reads_type_option_for_kmc('foo.bam'))
        self.assertEqual('-fq', bayestyper.get_reads_type_option_for_kmc('foo.fq'))
        self.assertEqual('-fq', bayestyper.get_reads_type_option_for_kmc('foo.fastq'))
        self.assertEqual('-fq', bayestyper.get_reads_type_option_for_kmc('foo.fq.gz'))
        self.assertEqual('-fq', bayestyper.get_reads_type_option_for_kmc('foo.fastq.gz'))
        with self.assertRaises(Exception):
            bayestyper.get_reads_type_option_for_kmc('foo.bar')


    def test_remove_starred_alts_from_vcf(self):
        '''test remove_starred_alts_from_vcf'''
        vcf_in = os.path.join(data_dir, 'remove_starred_alts_from_vcf.in.vcf')
        tmp_out = 'tmp.remove_starred_alts_from_vcf.vcf'
        bayestyper.remove_starred_alts_from_vcf(vcf_in, tmp_out)
        expect_vcf = os.path.join(data_dir, 'remove_starred_alts_from_vcf.expect.vcf')
        self.assertTrue(filecmp.cmp(expect_vcf, tmp_out, shallow=False))
        os.unlink(tmp_out)


    def test_run(self):
        '''test run'''
        ref_fasta = os.path.join(data_dir, 'run.ref.fa')
        reads_file = os.path.join(data_dir, 'run.reads.fq')
        vcf_files = {
            'tool1': os.path.join(data_dir, 'run.calls.1.vcf'),
            'tool2': os.path.join(data_dir, 'run.calls.2.vcf'),
        }
        # bayesTyper crashes on small test set, presumably because it
        # really needs lots of variants to gather whatever stats it's using.
        # So use the testing option, which runs as much as possible of the
        # bayestyper workflow, then just copies one of the input VCF files
        # to be the final output file
        outdir = 'tmp.bayestyper.run.out'
        bayestyper.run(ref_fasta, reads_file, vcf_files, outdir, testing=True)
        self.assertTrue(os.path.exists(os.path.join(outdir, '05.final.vcf')))
        shutil.rmtree(outdir)

