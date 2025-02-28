process OPENMM_SIMULATE {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(complex, stageAs: 'input.pdb'), path(ligand, stageAs: 'ligand.sdf')

    output:
    tuple val(meta), path("${meta.id}"), emit: results_dir
    tuple val(meta), path("${meta.id}/topology.pdb"), emit: topology
    tuple val(meta), path("${meta.id}/trajectory.xtc"), emit: trajectory
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // 定义参数变量
    def prefix = task.ext.prefix ?: "${meta.id}"
    def md_steps = params.md_steps ?: 5000
    def temperature = params.temperature ?: 300.0
    def write_interval = params.write_interval ?: 10
    def log_interval = params.log_interval ?: 250
    def protein_ff = params.protein_ff ?: 'amber14-all.xml'
    def solvent_ff = params.solvent_ff ?: 'amber14/tip3pfb.xml'
    def padding = params.padding ?: 1.0
    def ionic_strength = params.ionic_strength ?: 0.15
    def friction = params.friction ?: 1.0
    def stepsize = params.stepsize ?: 2.0
    def args = task.ext.args ?: ''
    def ligand_arg = ligand.name != 'NO_FILE' ? "--ligand ${ligand}" : ''

    """
    mkdir -p ${prefix}
    
    run_simulation.py \\
        --pdb ${complex} \\
        ${ligand_arg} \\
        --output ${prefix} \\
        --steps ${md_steps} \\
        --temperature ${temperature} \\
        --write-interval ${write_interval} \\
        --log-interval ${log_interval} \\
        --protein-ff "${protein_ff}" \\
        --solvent-ff "${solvent_ff}" \\
        --padding ${padding} \\
        --ionic-strength ${ionic_strength} \\
        --friction ${friction} \\
        --stepsize ${stepsize} \\
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        openmm: \$(python -c "import openmm; print(openmm.__version__)")
        rdkit: \$(python -c "import rdkit; print(rdkit.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/topology.pdb
    touch ${prefix}/trajectory.xtc
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        openmm: "8.0.0"
        rdkit: "2022.09.5"
    END_VERSIONS
    """
}
