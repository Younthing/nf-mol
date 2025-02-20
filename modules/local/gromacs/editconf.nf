process EDITCONF_GROMACS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/gromacs:2021.3--mpi_openmpi_h9969a6a_0'
        : 'docker.io/gromacs/gromacs:2022.2'}"

    publishDir [
        path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}/${meta.id}" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename },
    ]

    input:
    tuple val(meta), path(structure)

    output:
    tuple val(meta), path("*_newbox.gro"), emit: structure
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${structure.baseName}"

    // 从params读取参数,设置默认值
    def distance = params.editconf_distance ?: '1.0'
    def box_type = params.editconf_box_type ?: 'dodecahedron'
    def center = params.editconf_center ? '-c' : ''

    """
    gmx editconf \\
        -f ${structure} \\
        -o ${prefix}_newbox.gro \\
        ${center} \\
        -d ${distance} \\
        -bt ${box_type} \\
        ${args} &> editconf.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gromacs: \$(gmx --version | grep 'GROMACS version' | grep -oE '[0-9]+[0-9]+')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_newbox.gro

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gromacs: \$(gmx --version | grep 'GROMACS version' | grep -oE '[0-9]+[0-9]+')
    END_VERSIONS
    """
}
