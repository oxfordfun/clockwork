from clockwork import bayestyper

def run(options):
    vcf_files = {}
    with open(options.vcf_file_list) as f:
        for line in f:
            tool, filename = line.rstrip().split('\t')
            vcf_files[tool] = filename

    bayestyper.run(
        options.ref_fasta,
        options.reads_file,
        vcf_files,
        options.output_dir,
        kmc_mem=options.kmc_mem,
        testing=options.testing,
    )

