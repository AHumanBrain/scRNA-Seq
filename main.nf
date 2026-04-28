nextflow.enable.dsl=2

// =========================================================================
// PARAMETERS — Auto-detect entry point based on what the user provides
// =========================================================================
params.bcl_dir      = null      // Stage 1: Raw BCL run folder
params.fastq_dir    = null      // Stage 2: Directory of FASTQ files
params.matrix_dir   = null      // Stage 4: Pre-computed count matrix (MEX format)
params.sample_id    = 'sample'  // Sample identifier for outputs
params.star_index   = null      // Pre-built STAR genome index (optional)
params.genome_fasta = null      // Genome FASTA for building STAR index
params.genome_gtf   = null      // Gene annotation GTF for building STAR index
params.chemistry    = 'auto'    // 10x chemistry version (auto, v2, v3)
params.sample_sheet = null      // Illumina SampleSheet.csv for BCL demultiplexing

// =========================================================================
// STAGE 0: Demo Mode — Fetch 10x Genomics example datasets
// =========================================================================
process FETCH_DATA {
    publishDir "results/data", mode: 'copy'
    container 'ubuntu:22.04'

    output:
    path "pbmc3k_data", emit: data3k
    path "pbmc1k_data", emit: data1k

    script:
    """
    apt-get update && apt-get install -y wget tar
    mkdir -p pbmc3k_data pbmc1k_data

    wget https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
    tar -xzf pbmc3k_filtered_gene_bc_matrices.tar.gz -C pbmc3k_data

    wget https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_filtered_feature_bc_matrix.tar.gz
    tar -xzf pbmc_1k_v3_filtered_feature_bc_matrix.tar.gz -C pbmc1k_data
    """
}

// =========================================================================
// STAGE 1: BCL to FASTQ — Demultiplex raw Illumina sequencing output
// =========================================================================
process BCL2FASTQ {
    publishDir "results/fastqs", mode: 'copy'

    input:
    path bcl_run_dir
    path sample_sheet

    output:
    path "fastqs_out/*_R{1,2}_*.fastq.gz", emit: fastqs
    path "fastqs_out/Reports/", emit: reports
    path "fastqs_out/Stats/", emit: stats

    script:
    def sheet_arg = sample_sheet.name != 'NO_SHEET' ? "--sample-sheet ${sample_sheet}" : ''
    """
    mkdir -p fastqs_out
    bcl2fastq \\
        --runfolder-dir ${bcl_run_dir} \\
        --output-dir fastqs_out \\
        ${sheet_arg} \\
        --no-lane-splitting \\
        --minimum-trimmed-read-length 8 \\
        -r 4 -p 8 -w 4
    """
}

// =========================================================================
// STAGE 2: Read QC — Quality check raw FASTQ reads
// =========================================================================
process FASTQC {
    publishDir "results/fastqc", mode: 'copy'

    input:
    path fastq_files

    output:
    path "*.html", emit: html
    path "*.zip", emit: zip

    script:
    """
    fastqc --threads ${task.cpus} --outdir . ${fastq_files}
    """
}

// =========================================================================
// STAGE 3a: Build STAR Genome Index (only if no pre-built index provided)
// =========================================================================
process STAR_INDEX {
    publishDir "results/star_index", mode: 'copy'
    storeDir "star_index_cache"

    input:
    path genome_fasta
    path genome_gtf

    output:
    path "star_index", emit: index

    script:
    """
    mkdir -p star_index
    STAR \\
        --runMode genomeGenerate \\
        --runThreadN ${task.cpus} \\
        --genomeDir star_index \\
        --genomeFastaFiles ${genome_fasta} \\
        --sjdbGTFfile ${genome_gtf} \\
        --sjdbOverhang 100
    """
}

// =========================================================================
// STAGE 3b: Alignment — Map reads to genome and quantify with STARsolo
// =========================================================================
process STARSOLO {
    publishDir "results/starsolo", mode: 'copy'

    input:
    path star_index
    path reads       // paired FASTQs: R1 (barcode+UMI), R2 (cDNA)
    path whitelist

    output:
    path "Solo.out/Gene/filtered/", emit: matrix
    path "Aligned.sortedByCoord.out.bam", emit: bam, optional: true
    path "Log.final.out", emit: log

    script:
    """
    # Auto-download whitelist if placeholder provided
    WL="${whitelist}"
    if [ ! -s "\${WL}" ] || grep -q "NO_WHITELIST" "\${WL}" 2>/dev/null; then
        echo "Downloading 10x Chromium v3 barcode whitelist..."
        wget -q https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/3M-february-2018.txt.gz
        gunzip 3M-february-2018.txt.gz
        WL="3M-february-2018.txt"
    fi

    STAR \\
        --runThreadN ${task.cpus} \\
        --genomeDir ${star_index} \\
        --readFilesIn ${reads} \\
        --readFilesCommand zcat \\
        --soloType CB_UMI_Simple \\
        --soloCBwhitelist \${WL} \\
        --soloCBstart 1 --soloCBlen 16 \\
        --soloUMIstart 17 --soloUMIlen 12 \\
        --outSAMtype BAM SortedByCoordinate \\
        --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \\
        --soloFeatures Gene \\
        --soloCellFilter EmptyDrops_CR \\
        --clipAdapterType CellRanger4
    """
}

// =========================================================================
// STAGE 4: Integrated Scanpy Analysis (Our unique value-add)
// =========================================================================
process SCANPY_ANALYSIS {
    publishDir "results", mode: 'copy'

    input:
    path matrix_dirs

    output:
    path "figures/", emit: figures
    path "pbmc_integrated.h5ad", emit: h5ad
    path "cell_type_assignments.csv", emit: csv
    path "scanpy_analysis.log", emit: log

    script:
    """
    python3 ${projectDir}/bin/analyze.py ${matrix_dirs} > scanpy_analysis.log 2>&1
    """
}

// =========================================================================
// STAGE 5: MultiQC — Aggregate all QC reports
// =========================================================================
process MULTIQC {
    publishDir "results", mode: 'copy'

    input:
    path "inputs/*"

    output:
    path "multiqc_report.html", emit: report, optional: true
    path "multiqc_data/", emit: data, optional: true

    script:
    """
    echo "# plot_type: 'table'" > scanpy_stats_mqc.tsv
    echo "Sample\tPipeline_Status\tEntry_Point" >> scanpy_stats_mqc.tsv
    echo "${params.sample_id}\tComplete\t${params.bcl_dir ? 'BCL' : params.fastq_dir ? 'FASTQ' : params.matrix_dir ? 'Matrix' : 'Demo'}" >> scanpy_stats_mqc.tsv

    multiqc . --title "scRNA-Seq Pipeline Report" -n "multiqc_report.html"
    """
}

// =========================================================================
// WORKFLOW — Conditional execution based on entry point
// =========================================================================
workflow {

    // Collect QC channels for MultiQC aggregation
    mqc_inputs = Channel.empty()

    // -------------------------------------------------------------------
    // PATH 1: User provides a BCL run folder → full pipeline
    // -------------------------------------------------------------------
    if (params.bcl_dir) {
        bcl_ch     = Channel.fromPath(params.bcl_dir, checkIfExists: true)
        sheet_ch   = params.sample_sheet
                        ? Channel.fromPath(params.sample_sheet, checkIfExists: true)
                        : Channel.fromPath("${projectDir}/assets/NO_SHEET")

        BCL2FASTQ(bcl_ch, sheet_ch)
        FASTQC(BCL2FASTQ.out.fastqs.flatten())

        // Build or use pre-built STAR index
        if (params.star_index) {
            star_idx = Channel.fromPath(params.star_index, checkIfExists: true)
        } else {
            STAR_INDEX(
                Channel.fromPath(params.genome_fasta, checkIfExists: true),
                Channel.fromPath(params.genome_gtf, checkIfExists: true)
            )
            star_idx = STAR_INDEX.out.index
        }

        whitelist = Channel.fromPath("${projectDir}/assets/whitelists/3M-february-2018.txt")
        STARSOLO(star_idx, BCL2FASTQ.out.fastqs.collect(), whitelist)

        SCANPY_ANALYSIS(STARSOLO.out.matrix.collect())

        mqc_inputs = FASTQC.out.zip
            .mix(STARSOLO.out.log)
            .mix(BCL2FASTQ.out.stats.ifEmpty([]))
            .mix(SCANPY_ANALYSIS.out.log)
            .collect()
    }

    // -------------------------------------------------------------------
    // PATH 2: User provides FASTQ files → skip demultiplexing
    // -------------------------------------------------------------------
    else if (params.fastq_dir) {
        fastq_ch = Channel.fromPath("${params.fastq_dir}/*.fastq.gz", checkIfExists: true)
        FASTQC(fastq_ch)

        if (params.star_index) {
            star_idx = Channel.fromPath(params.star_index, checkIfExists: true)
        } else {
            STAR_INDEX(
                Channel.fromPath(params.genome_fasta, checkIfExists: true),
                Channel.fromPath(params.genome_gtf, checkIfExists: true)
            )
            star_idx = STAR_INDEX.out.index
        }

        whitelist = Channel.fromPath("${projectDir}/assets/whitelists/3M-february-2018.txt")
        STARSOLO(star_idx, fastq_ch.collect(), whitelist)

        SCANPY_ANALYSIS(STARSOLO.out.matrix.collect())

        mqc_inputs = FASTQC.out.zip
            .mix(STARSOLO.out.log)
            .mix(SCANPY_ANALYSIS.out.log)
            .collect()
    }

    // -------------------------------------------------------------------
    // PATH 3: User provides a count matrix → skip upstream entirely
    // -------------------------------------------------------------------
    else if (params.matrix_dir) {
        matrix_ch = Channel.fromPath(params.matrix_dir, checkIfExists: true)
        SCANPY_ANALYSIS(matrix_ch.collect())

        mqc_inputs = SCANPY_ANALYSIS.out.log.collect()
    }

    // -------------------------------------------------------------------
    // PATH 4: Demo mode — fetch 10x example datasets (default)
    // -------------------------------------------------------------------
    else {
        FETCH_DATA()
        SCANPY_ANALYSIS(FETCH_DATA.out.data3k.mix(FETCH_DATA.out.data1k).collect())

        mqc_inputs = SCANPY_ANALYSIS.out.figures
            .mix(SCANPY_ANALYSIS.out.log)
            .collect()
    }

    // Aggregate all QC into a single report
    MULTIQC(mqc_inputs)
}
