import logging
import os
import subprocess

import cluster_vcf_records
import pyfastaq

from clockwork import utils


def normalise_vcf(vcf_in, ref_fasta, vcf_out, bcftools):
    command = f'{bcftools} norm -cx -f {ref_fasta} {vcf_in} > {vcf_out}'
    logging.info(f'Normalise VCF file: {command}')
    utils.syscall(command)


def normalise_vcfs(vcf_files, ref_fasta, bcftools, outprefix):
    normalised_vcfs = {}

    for tool, vcf_file in vcf_files.items():
        normalised_vcf = outprefix + '.' + str(len(normalised_vcfs)) + '.vcf'
        normalise_vcf(vcf_file, ref_fasta, normalised_vcf, bcftools)
        normalised_vcfs[tool] = normalised_vcf

    return normalised_vcfs

# From bayesTyper documentation:
# Important: bayesTyperTools combine requires the vcf header to contain contig entries (e.g.##contig=<ID=8,length=146364022>) for all reference sequences containing variants in the vcf; the contigs further need to appear in the same order in the header and for the variant entries.
#
# This function adds contig names to the header.
# It also gets the sample name, which is used elsewhere
def vcf_add_contigs_to_header_and_get_sample(vcf_in, vcf_out, ref_lengths):
    header_lines, vcf_records = cluster_vcf_records.vcf_file_read.vcf_file_to_list(vcf_in)
    sample = cluster_vcf_records.vcf_file_read.get_sample_name_from_vcf_header_lines(header_lines)
    header_lines = [x for x in header_lines if not x.startswith('##contig=<ID=')]
    seen_contigs = set()
    contig_header_lines = []

    for record in vcf_records:
        if record.CHROM not in seen_contigs:
            seen_contigs.add(record.CHROM)
            contig_header_lines.append(f'##contig=<ID={record.CHROM},length={ref_lengths[record.CHROM]}>')

    with open(vcf_out, 'w') as f:
        print(*header_lines[:-1], sep='\n', file=f)
        print(*contig_header_lines, sep='\n', file=f)
        print(header_lines[-1], file=f)
        print(*vcf_records, sep='\n', file=f)

    return sample


def fix_vcf_headers_and_get_sample_name(vcf_files, outprefix, ref_lengths):
    samples = set()
    fixed_vcfs = {}

    for tool, vcf_file in vcf_files.items():
        fixed_vcf = outprefix + '.' + str(len(fixed_vcfs)) + '.vcf'
        sample = vcf_add_contigs_to_header_and_get_sample(vcf_file, fixed_vcf, ref_lengths)
        samples.add(sample)
        fixed_vcfs[tool] = fixed_vcf

    if len(samples) > 1:
        logging.warning('More than one sample name found in VCF files. Going to use one of them at random')

    sample = samples.pop()
    return fixed_vcfs, sample


def get_reads_type_option_for_kmc(filename):
    extensions = {'.bam': '-fbam',
        '.fq': '-fq',
        '.fq.gz': '-fq',
        '.fastq': '-fq',
        '.fastq.gz': '-fq',
    }

    for extension in extensions:
        if filename.endswith(extension):
            return extensions[extension]

    raise Exception(f'Reads file format not recognised: {filename}')


def remove_starred_alts_from_vcf(vcf_in, vcf_out):
    with open(vcf_in) as f_in, open(vcf_out, 'w') as f_out:
        for line in f_in:
            if line.startswith('#') or '*' not in line.split('\t')[4]:
                print(line, end='', file=f_out)


def run(ref_fasta, reads_file, vcf_files, outdir, bayestyper=None, kmc=None, bcftools=None, kmc_mem=2):
    if bayestyper is None:
        bayestyper = os.environ.get('CLOCKWORK_BAYESTYPER', 'bayesTyper')
    if kmc is None:
        kmc = os.environ.get('CLOCKWORK_KMC', 'kmc')
    if bcftools is None:
        bcftools = os.environ.get('CLOCKWORK_BCFTOOLS', 'bcftools')

    bayestyper_tools = f'{bayestyper}Tools'

    reads_file = os.path.abspath(reads_file)
    if not os.path.exists(reads_file):
        raise FileNotFoundError(f'Error finding reads file {reads_file}')

    ref_fasta = os.path.abspath(ref_fasta)
    if not os.path.exists(ref_fasta):
        raise FileNotFoundError(f'Error finding reference fasta file {ref_fasta}')


    ref_fasta_fai = f'{ref_fasta}.fai'
    if not os.path.exists(ref_fasta_fai):
        raise FileNotFoundError(f'Error finding reference fasta fai file {ref_fasta_fai}. Please run samtools faidx {ref_fasta}')

    for tool in vcf_files:
        if not os.path.exists(vcf_files[tool]):
            raise FileNotFoundError(f'Error finding VCF file {vcf_files[tool]}')
        vcf_files[tool] = os.path.abspath(vcf_files[tool])

    cwd = os.getcwd()
    outdir = os.path.abspath(outdir)
    os.mkdir(outdir)
    os.chdir(outdir)
    f_log = open('log.txt', 'w')
    ref_lengths = {}
    pyfastaq.tasks.lengths_from_fai(ref_fasta_fai, ref_lengths)
    print('Adding relevant info to headers of input VCF files', file=f_log)
    fixed_vcfs, sample_name = fix_vcf_headers_and_get_sample_name(vcf_files, '00.fix_header', ref_lengths)
    print('Normalising VCF files using bcftools', file=f_log)
    normalised_vcfs = normalise_vcfs(fixed_vcfs, ref_fasta, bcftools, '01.normalised')
    kmc_reads_option = get_reads_type_option_for_kmc(reads_file)
    kmc_out = '02.kmc'
    command = f'kmc -m{kmc_mem} -k55 -ci1 {kmc_reads_option} {reads_file} {kmc_out} {outdir}'
    print('Running kmc:', command, file=f_log)
    utils.syscall(command)
    command = f'{bayestyper_tools} makeBloom -k {kmc_out}'
    print('Running bayesTyperTools makeBloom:', command, file=f_log)
    utils.syscall(command)
    vcfs_string = ','.join([f'{x}:{normalised_vcfs[x]}' for x in normalised_vcfs])
    typertools_combine_out = '03.typertools.combine'
    utils.syscall(f'{bayestyper_tools} combine -v {vcfs_string} -o {typertools_combine_out} -z')
    ploidy_file = '04.ploidy.tsv'
    with open(ploidy_file, 'w') as f:
        for ref_name in ref_lengths:
            print(ref_name, 2, 2, sep='\t', file=f)
    samples_tsv = '04.samples.tsv'
    with open(samples_tsv, 'w') as f:
        print(sample_name, 'M', os.path.join(outdir, kmc_out), sep='\t', file=f)

    command = f'{bayestyper} cluster -v {typertools_combine_out}.vcf.gz -s {samples_tsv} -g {ref_fasta}'
    print('Running bayesTyper cluster:', command, file=f_log)
    utils.syscall(command)
    # This commented out command is for bayesTyper v1.4.1. But we can't use that version
    # because it crashes:
    # https://github.com/bioinformatics-centre/BayesTyper/issues/4
    # So we're using 1.3.1 instead. 1.3.1 doesn't have the -y option, but
    # leave the command here as a note for the fuutre.
    #utils.syscall(f'{bayestyper} genotype -y {ploidy_file} -v bayestyper_unit_1/variant_clusters.bin -c bayestyper_cluster_data -s {samples_tsv} -g {ref_fasta} -o bayestyper_unit_1/bayestyper')
    command = f'{bayestyper} genotype -v bayestyper_unit_1/variant_clusters.bin -c bayestyper_cluster_data -s {samples_tsv} -g {ref_fasta} -o bayestyper_unit_1/bayestyper'
    print('Running bayesTyper genotype:', command, file=f_log)
    utils.syscall(command)
    bayestyper_vcf = 'bayestyper_unit_1/bayestyper.vcf'
    filtered_bayestyper_vcf = '05.final.vcf'
    remove_starred_alts_from_vcf(bayestyper_vcf, filtered_bayestyper_vcf)
    f.close()
    os.chdir(cwd)

