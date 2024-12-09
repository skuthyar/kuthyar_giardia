# Define paths and parameters
SINGULARITY = "/projects/b1057/qiime2-core-2019-4.simg"
INPUT_PATH = "/projects/p30050/skuthyar/AcarayaProjectNewVersion"
OUTPUT_PATH = INPUT_PATH
METADATA_FILE = "/projects/p30050/skuthyar/AcarayaProjectNewVersion/metadata-allhowlers.txt"
GIARDIA_METADATA_FILE = "/projects/p30050/skuthyar/AcarayaProjectNewVersion/giardia-metadata.txt"
GG_CLASSIFIER = "/projects/b1057/gg-13-8-99-nb-classifier-qiime2019-4.qza"

rule all:
    input:
        expand("{OUTPUT_PATH}/giardia-taxonomy-table-from-biom.txt", OUTPUT_PATH=OUTPUT_PATH)

rule import_data:
    input:
        manifest = "/projects/p30050/skuthyar/AcarayaProjectNewVersion/howlermanifest.csv"
    output:
        paired_end_demux = f"{OUTPUT_PATH}/paired-end-demux.qza"
    shell:
        "singularity exec -B /projects/p30050 -B /projects/p30050/skuthyar {SINGULARITY} "
        "qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' "
        "--input-path {input.manifest} --output-path {output.paired_end_demux} "
        "--input-format PairedEndFastqManifestPhred33"

rule dada2_denoise:
    input:
        demux = f"{OUTPUT_PATH}/paired-end-demux.qza"
    output:
        table = f"{OUTPUT_PATH}/table-ee5.qza",
        rep_seqs = f"{OUTPUT_PATH}/rep-seqs-ee5.qza",
        stats = f"{OUTPUT_PATH}/dada2denoising-stats-ee5.qza"
    shell:
        "singularity exec -B /projects/p30050 -B /projects/p30050/skuth
