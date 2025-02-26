process PREPARE_LIGAND {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(sdf_file, stageAs: 'input.sdf'), val(smiles)

    output:
    tuple val(meta), path("*_prepared.sdf"), emit: ligand
    tuple val(meta), path("*_comparison.png"), optional: true, emit: image
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def image_out = task.ext.save_image ? "--image ${prefix}_comparison.png" : ""

    """
    prepare_sdf.py \\
        ${sdf_file} \\
        --smiles '${smiles}' \\
        --output ${prefix}_prepared.sdf \\
        ${image_out} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rdkit: \$(python -c "import rdkit; print(rdkit.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_prepared.sdf
    touch ${prefix}_comparison.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rdkit: "2024.03.5"
    END_VERSIONS
    """
}
