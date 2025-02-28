process OPENMM_MDANALYSIS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(topology), path(trajectory)

    output:
    tuple val(meta), path("${meta.id}"), emit: analysis_dir
    tuple val(meta), path("${meta.id}/aligned_md.html"), emit: aligned_trajectory
    tuple val(meta), path("${meta.id}/rmsd_results.csv"), emit: rmsd_data
    tuple val(meta), path("${meta.id}/rmsd_plot.png"), emit: rmsd_plot
    tuple val(meta), path("${meta.id}/frame_rmsd_plot.png"), emit: frame_rmsd
    tuple val(meta), path("${meta.id}/pocket_view.html"), emit: pocket_view
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ligand_name = params.ligand_name ?: 'LIG'
    """
    visualization.py \\
        ${topology} \\
        ${trajectory} \\
        ${ligand_name} \\
        --output_dir ${prefix} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mdanalysis: \$(python -c "import MDAnalysis; print(MDAnalysis.__version__)")
        openmm: \$(python -c "import openmm; print(openmm.__version__)")
        nglview: \$(python -c "import nglview; print(nglview.__version__)")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${meta.id}
    touch ${meta.id}/aligned_md.html
    touch ${meta.id}/rmsd_results.csv
    touch ${meta.id}/rmsd_plot.png
    touch ${meta.id}/frame_rmsd_plot.png
    touch ${meta.id}/pocket_view.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "3.8.0"
        mdanalysis: "2.0.0"
        openmm: "8.0.0"
        nglview: "3.0.3"
    END_VERSIONS
    """
}
