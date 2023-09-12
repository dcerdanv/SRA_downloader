import pandas as pd
from snakemake.utils import validate

#### GLOBAL PARAMETERS ####

configfile: "config.yaml"
validate(config, schema="schemas/config.schema.yaml")

OUTDIR = config['outdir']
LOGDIR = config['logdir']
TMPDIR = config['tmpdir']


compress = config['parameters']['compress']
remove_sra = config['parameters']['remove_sra']

#### LOAD SAMPLES TABLE ###
samples = pd.read_table(config["samples"]).set_index("sample_ID", drop=False)
validate(samples, schema="schemas/samples.schema.yaml")


#### GLOBAL scope functions ####
def get_resource(rule,resource) -> int:
	'''
	Attempt to parse config.yaml to retrieve resources available for a given
	rule. It will revert to default if a key error is found. Returns an int.
	with the allocated resources available for said rule. Ex: "threads": 1
	'''

	try:
		return config['resources'][rule][resource]
	except KeyError: # TODO: LOG THIS
		print(f'Failed to resolve resource for {rule}/{resource}: using default parameters')
		return config["resources"]['default'][resource]


def get_all_input():

    all_input = []

    if compress == True:
        for sample in samples['sample_ID']:
            all_input += [f"{OUTDIR}/{sample}.fastq.gz"]
            all_input += [f"{LOGDIR}/{sample}/{sample}_validate.log"]

            if remove_sra == True:
                all_input += [f"{LOGDIR}/{sample}/{sample}_remove.log"]
    else:
        for sample in samples['sample_ID']:
            all_input += [f"{OUTDIR}/{sample}.fastq"]
            all_input += [f"{LOGDIR}/{sample}/{sample}_validate.log"]

            if remove_sra == True:
                all_input += [f"{LOGDIR}/{sample}/{sample}_remove.log"]

    return all_input

def get_remove_sra_input():
    
    if compress == True:
        remove_sra = [f"{OUTDIR}/{{sample}}.fastq.gz"]
    else:
        remove_sra = [f"{OUTDIR}/{{sample}}.fastq"]

    return remove_sra



### RULES ###
rule all:
	input:
		get_all_input()


rule prefetch_sra:
    output:
        prefetch_sra=f"{OUTDIR}/{{sample}}/{{sample}}.sra"
    threads:
        get_resource('prefetch_sra', 'threads')
    resources:
        mem_mb=get_resource('prefetch_sra', 'mem_mb'),
        runtime=get_resource('prefetch_sra', 'runtime')
    params:
        sample=f"{{sample}}",
        # optional parameters
        extra=config['parameters']['prefetch_options']
    log:
        f"{LOGDIR}/{{sample}}/{{sample}}_prefetch.log"
    conda:
        'envs/sra_downloader.yaml'
    shell:
        "echo Downloading {params.sample} SRA object. > {log} 2>&1 && date >> {log} 2>&1 && prefetch --output-directory {OUTDIR} {params.extra} {params.sample} >> {log} 2>&1 && date >> {log} 2>&1"


rule validate_sra:
    input:
        prefetch_sra=f"{OUTDIR}/{{sample}}/{{sample}}.sra"
    output:
        f"{LOGDIR}/{{sample}}/{{sample}}_validate.log"
    threads:
        get_resource('validate_sra', 'threads')
    resources:
        mem_mb=get_resource('validate_sra', 'mem_mb'),
        runtime=get_resource('validate_sra', 'runtime')
    params:
        sample=f"{{sample}}",
        # optional parameters
        extra=config['parameters']['validate_options']
    log:
        f"{LOGDIR}/{{sample}}/{{sample}}_validate.log"
    conda:
        'envs/sra_downloader.yaml'
    shell:
        "echo Validating {params.sample} SRA object. > {log} 2>&1 && date >> {log} 2>&1 && vdb-validate {params.extra} {input.prefetch_sra} >> {log} 2>&1 && date >> {log} 2>&1"


rule dump_fastq:
    input:
        validation_log = f"{LOGDIR}/{{sample}}/{{sample}}_validate.log",
        prefetch_sra=f"{OUTDIR}/{{sample}}/{{sample}}.sra"
    output:
        fastq=f"{OUTDIR}/{{sample}}.fastq"
    threads:
        get_resource('dump_fastq', 'threads')
    resources:
        mem_mb=get_resource('dump_fastq', 'mem_mb'),
        runtime=get_resource('dump_fastq', 'runtime')
    params:
        sample=f"{OUTDIR}/{{sample}}",
        # optional parameters
        extra=config['parameters']['dump_fastq_options']
    log:
        f"{LOGDIR}/{{sample}}/{{sample}}_extract.log"
    conda:
        'envs/sra_downloader.yaml'
    shell:
        "echo Extracting fastq from {params.sample}. > {log} && 2>&1 && date >> {log} 2>&1 && fasterq-dump -v -L info {params.extra} --outdir {OUTDIR} --temp {TMPDIR} -e {threads} {params.sample} && echo Extraction finished. >> {log} 2>&1 && date >> {log} 2>&1"


rule compress_fastq:
    input:
        fastq = f"{OUTDIR}/{{sample}}.fastq"
    output:
        fastq_gz=f"{OUTDIR}/{{sample}}.fastq.gz"
    threads:
        get_resource('compress_fastq', 'threads')
    resources:
        mem_mb=get_resource('compress_fastq', 'mem_mb'),
        runtime=get_resource('compress_fastq', 'runtime')
    params:
        # optional parameters
        extra=config['parameters']['pigz_options']
    log:
        f"{LOGDIR}/{{sample}}/{{sample}}_dump.log"
    conda:
        'envs/sra_downloader.yaml'
    shell:
        "echo Compressing {input.fastq} to GZ format. > {log} && 2>&1  date >> {log} 2>&1 && pigz --processes {threads} {params.extra} {input.fastq} >> {log} 2>&1 && date >> {log} 2>&1"


rule remove_sra:
    input:
        get_remove_sra_input()
    output:
        f"{LOGDIR}/{{sample}}/{{sample}}_remove.log"
    threads:
        get_resource('default', 'threads')
    resources:
        mem_mb=get_resource('default', 'mem_mb'),
        runtime=get_resource('default', 'runtime')
    params:
        sample=f"{{sample}}",
        sample_dir = [f"{OUTDIR}/{{sample}}"]

    log:
        f"{LOGDIR}/{{sample}}/{{sample}}_remove.log"
    shell:
        "rm -r {params.sample_dir} && echo SRA folder for sample {params.sample} removed successfully > {log} 2>&1 && date >> {log} 2>&1"