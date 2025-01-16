process PDB2GMX_GROMACS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/gromacs:2021.3--mpi_openmpi_h9969a6a_0'
        : 'biocontainers/gromacs:2021.3--mpi_openmpi_h9969a6a_0'}"

    // 按样本保存目录
    publishDir [
        path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}/${meta.id}" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename },
    ]

    input:
    tuple val(meta), path(pdb)

    output:
    tuple val(meta), path("*_processed.gro"), emit: structure
    tuple val(meta), path("*_topol.top"), emit: topology
    tuple val(meta), path("*.itp"), emit: include_topology
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${pdb.baseName}"

    // 从params读取参数,设置默认值
    def water = params.pdb2gmx_water_model ?: 'tip3p'
    def ff = params.pdb2gmx_forcefield ?: 'charmm27'
    def h_flag = params.pdb2gmx_ignore_h ? '-ignh' : ''

    """
    # 首先检查PDB文件
    if ! gmx check -f ${pdb} &> check.log; then
        echo "PDB文件需要修复,退出代码1"
        exit 1
    fi
    gmx pdb2gmx \\
        -f ${pdb} \\
        -o ${prefix}_processed.gro \\
        -p ${prefix}_topol.top \\
        -i ${prefix}_posre.itp \\
        -water ${water} \\
        -ff ${ff} \\
        ${h_flag} \\
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
