process MERGE_COMPLEX {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(protein, stageAs: 'protein.pdb'), path(ligand, stageAs: 'ligand.pdb')

    output:
    tuple val(meta), path("*_complex.pdb"), emit: complex
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    merge_complex.py \\
        ${protein} \\
        ${ligand} \\
        --output ${prefix}_complex.pdb \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mdtraj: \$(pip list | grep mdtr | cut -d ' ' -f 3- | tr -d ' ')
        openmm: \$(python -c "import openmm; print(openmm.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_complex.pdb

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mdtraj: "1.10.0"
        openmm: "8.2"
    END_VERSIONS
    """
}
