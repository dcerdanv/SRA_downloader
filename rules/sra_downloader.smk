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
        '../envs/sra_downloader.yaml'
    shell:
        "echo Downloading {params.sample} SRA object. > {log} 2>&1 && date >> {log} 2>&1 && prefetch --output-directory {OUTDIR} {params.extra} {params.sample} >> {log} 2>&1 && date >> {log} 2>&1"


rule validate_sra:
    input:
        prefetch_sra=f"{OUTDIR}/{{sample}}/{{sample}}.sra"
    output:
        touch(temp(f"{LOGDIR}/{{sample}}/flags/{{sample}}_validate.done"))
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
        '../envs/sra_downloader.yaml'
    shell:
        "echo Validating {params.sample} SRA object. > {log} 2>&1 && date >> {log} 2>&1 && vdb-validate {params.extra} {input.prefetch_sra} >> {log} 2>&1 && date >> {log} 2>&1"


checkpoint dump_fastq:
    input:
        validation_log = f"{LOGDIR}/{{sample}}/flags/{{sample}}_validate.done",
        prefetch_sra=f"{OUTDIR}/{{sample}}/{{sample}}.sra"
    output:
        directory(f"{OUTDIR}/{{sample}}/fastq")
    threads:
        get_resource('dump_fastq', 'threads')
    resources:
        mem_mb=get_resource('dump_fastq', 'mem_mb'),
        runtime=get_resource('dump_fastq', 'runtime')
    params:
        sample=f"{OUTDIR}/{{sample}}",
        # optional parameters
        extra=config['parameters']['fasterq-dump_options']
    log:
        f"{LOGDIR}/{{sample}}/{{sample}}_dump_fastq.log"
    conda:
        '../envs/sra_downloader.yaml'
    shell:
        "echo Extracting fastq from {params.sample}. > {log} && 2>&1 && date >> {log} 2>&1 && fasterq-dump -v -L info {params.extra} --outdir {output} --temp {TMPDIR} -e {threads} {params.sample} && echo Extraction finished. >> {log} 2>&1 && date >> {log} 2>&1"


rule move_fastq:
    input:
        fastq = f"{OUTDIR}/{{sample}}/fastq/{{sample_part}}.fastq"
    output:
        fastq_gz = f"{OUTDIR}/{{sample}}/{{sample_part}}.fastq"
    threads:
        get_resource('default', 'threads')
    resources:
        mem_mb=get_resource('default', 'mem_mb'),
        runtime=get_resource('default', 'runtime')
    params:
        sample_dir = f"{OUTDIR}/{{sample}}",
    shell:
        "mv {input.fastq} {params.sample_dir}"


rule compress_fastq:
    input:
        fastq = f"{OUTDIR}/{{sample}}/{{sample_part}}.fastq"
    output:
        fastq_gz = f"{OUTDIR}/{{sample}}/{{sample_part}}.fastq.gz"
    threads:
        get_resource('compress_fastq', 'threads')
    resources:
        mem_mb=get_resource('compress_fastq', 'mem_mb'),
        runtime=get_resource('compress_fastq', 'runtime')
    params:
        # optional parameters
        extra=config['parameters']['pigz_options'],
        sample_dir = f"{OUTDIR}/{{sample}}",
    log:
        f"{LOGDIR}/{{sample}}/{{sample_part}}_compress.log"
    conda:
        '../envs/sra_downloader.yaml'
    shell:
        "echo Compressing {input.fastq} to GZ format. > {log} && 2>&1  date >> {log} 2>&1 && pigz --processes {threads} {params.extra} {input.fastq} >> {log} 2>&1 && date >> {log} 2>&1"


def get_remove_sra_input(wildcards):
        
    checkpoint_output = checkpoints.dump_fastq.get(**wildcards).output[0]

    if compress == True:
        remove_sra_input = expand(f"{OUTDIR}/{wildcards.sample}/{{sample_part}}.fastq.gz", \
                    sample_part=glob_wildcards(os.path.join(checkpoint_output, f"{{sample_part}}.fastq")).sample_part)
    else:
        remove_sra_input = expand(f"{OUTDIR}/{wildcards.sample}/{{sample_part}}.fastq", \
                    sample_part=glob_wildcards(os.path.join(checkpoint_output, f"{{sample_part}}.fastq")).sample_part)
    
    return remove_sra_input


rule remove_sra:
    input:
        get_remove_sra_input
    output:
        touch(f"{LOGDIR}/{{sample}}/flags/{{sample}}_remove.done")
    threads:
        get_resource('default', 'threads')
    resources:
        mem_mb=get_resource('default', 'mem_mb'),
        runtime=get_resource('default', 'runtime')
    params:
        sample=f"{{sample}}",
        sample_sra = f"{OUTDIR}/{{sample}}/{{sample}}.sra",
        fastq_dir = f"{OUTDIR}/{{sample}}/fastq",
    log:
        f"{LOGDIR}/{{sample}}/{{sample}}_remove.log"
    shell:
        "rm -r {params.fastq_dir} && rm -r {params.sample_sra} && echo SRA folder for sample {params.sample} removed successfully > {log} 2>&1 && date >> {log} 2>&1"


rule not_remove_sra:
    input:
        get_remove_sra_input
    output:
        touch(f"{LOGDIR}/{{sample}}/flags/{{sample}}_not_remove.done")
    threads:
        get_resource('default', 'threads')
    resources:
        mem_mb=get_resource('default', 'mem_mb'),
        runtime=get_resource('default', 'runtime')
    params:
        fastq_dir = f"{OUTDIR}/{{sample}}/fastq",
    shell:
        "rm -r {params.fastq_dir}"


def get_all_input(wildcards):

    all_input = []

    for sample in samples['sample_ID']:

        if remove_sra == True:
            all_input += [f"{LOGDIR}/{sample}/flags/{sample}_remove.done"]
        else:
            all_input += [f"{LOGDIR}/{sample}/flags/{sample}_not_remove.done"]

    return all_input

