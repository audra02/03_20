import glob
import os

# List your sample names (without file extension)
samples = "chr17"

input_path = "/home/mif/03_20"
  
output_path = "results"

thread_no = 2

reference_file = "/home/shared/raw_data/raw_data_exercise_20250320/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta"

path_to_vep_cache = "/home/shared/raw_data/raw_data_exercise_20250320"

rule all:
    input:
        vcf = expand("{output}/{sample}.vcf.gz", output=output_path, sample=samples),
        annotation = expand("{output}/{sample}_annotated.vcf", output=output_path, sample=samples),
        filtered_vcf = expand("{output}/{sample}_annotated_filtered.vcf", output=output_path, sample=samples),


rule variant_calling:
    input:
        cram = expand("{input}/{{sample}}.cram", input=input_path),
        cram_index = expand("{input}/{{sample}}.cram.crai", input=input_path)
    output:
        vcf = "{output}/{sample}.vcf.gz",
        vcf_index = "{output}/{sample}.vcf.gz.tbi",
    params:
        ref = reference_file,
        gatk_opts = "--native-pair-hmm-threads 4 ",
        quality = 20
    log:
        "{output}/{sample}.haplotypecaller.log"
    threads: thread_no
    conda: "envs/variant_calling.yml"
    shell:
        """
        gatk HaplotypeCaller \
            -R {params.ref} \
            -I {input.cram} \
            -O {output.vcf} \
            --min-base-quality-score {params.quality} \
            {params.gatk_opts} \
            > {log} 2>&1
        """

rule vep:
    input: 
        vcf = "{output}/{sample}.vcf.gz",
        vcf_index = "{output}/{sample}.vcf.gz.tbi",
    output:
        annotated_vcf = "{output}/{sample}_annotated.vcf",
    conda: "envs/vep.yml"
    params:
        vcf = "{output}/{sample}.vcf",
        cache = path_to_vep_cache
    threads: thread_no
    shell:
        """
            gzip -dkc {input.vcf} > {params.vcf}
            vep -i {params.vcf} \
            -o {output.annotated_vcf} \
            --vcf \
            --cache  --dir_cache {params.cache} \
            --af --af_gnomade --af_gnomadg --sift b \
            --polyphen b --symbol --protein --biotype \
            --af --variant_class  \
            --force_overwrite 
        """

rule filter_vcf:
    input:
        annotated_vcf = "{output}/{sample}_annotated.vcf",
    output:
        filtered_vcf = "{output}/{sample}_annotated_filtered.vcf",
    conda: "envs/vep.yml"   
    shell:
        """
            filter_vep \
            -i {input.annotated_vcf} \
            -o {output.filtered_vcf} \
            --vcf \
            --filter "Consequence match missense_variant or Consequence match stop_gained or Consequence match frameshift_variant" \
            --filter "gnomADg_AF < 0.05" \
            --filter "CLIN_SIG match uncertain_significance or CLIN_SIG match pathogenic or CLIN_SIG match likely_pathogenic or CLIN_SIG match risk_factor" \
            --force_overwrite
        """