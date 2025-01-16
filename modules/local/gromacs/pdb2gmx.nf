process PDB2GMX_GROMACS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/gromacs:2021.3--mpi_openmpi_h9969a6a_0'
        : 'biocontainers/gromacs:2021.3--mpi_openmpi_h9969a6a_0'}"

    input:
    tuple val(meta), path(pdb)

    output:
    tuple val(meta), path("*.gro"), emit: structure
    tuple val(meta), path("*.top"), emit: topology
    tuple val(meta), path("*.itp"), emit: include_topology
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gmx pdb2gmx \
        -f ${pdb} \
        -o ${prefix}_processed.gro \
        -p ${prefix}_topol.top \
        -i ${prefix}_posre.itp \
        -water tip3p \
        -ff charmm27 \
        -ignh  # 忽略氢原子\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gromacs: \$(gmx --version | grep 'GROMACS version' | grep -oE '[0-9]+[0-9]+')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_processed.gro
    touch ${prefix}_topol.top
    touch ${prefix}_posre.itp

    cat <<-END_VERSIONS > versions.yml
    
    "${task.process}":
        gromacs: \$(gmx --version | gmx --version | grep 'GROMACS version' | grep -oE '[0-9]+[0-9]+')
    END_VERSIONS
    """
}
