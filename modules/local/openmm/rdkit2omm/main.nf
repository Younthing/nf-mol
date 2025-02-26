process RDKIT2OMM {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(sdf_file, stageAs: 'input.sdf')

    output:
    tuple val(meta), path("*_omm.pdb"), emit: pdb
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def name = task.ext.name ?: "LIG"

    """
    rdkit2omm.py \\
        ${sdf_file} \\
        --name ${name} \\
        --output ${prefix}_omm.pdb \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rdkit: \$(python -c "import rdkit; print(rdkit.__version__)")
        openff: \$(python -c "import openff.toolkit; print(openff.toolkit.__version__)")
        openmm: \$(python -c "import openmm; print(openmm.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_omm.pdb

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rdkit: "2024.03.5"
        openff: "0.16.7"
        openmm: "8.2"
    END_VERSIONS
    """
}
