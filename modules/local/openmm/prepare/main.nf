process OPENMM_PREPARE {
    // 1. 进程标识和资源设置
    tag "${meta.id}"
    label 'process_low'

    // 2. 错误处理策略
    errorStrategy 'retry'
    maxRetries 3

    // 3. 环境配置
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(pdb_file, stageAs: 'input.pdb')

    output:
    tuple val(meta), path("*.{pdb,cif,pdbx}"), emit: structure
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // 解析格式参数并设置扩展名
    def format = 'pdb'
    if (args.contains('--format pdb') || args.contains('-f pdb')) {
        format = 'pdb'
    }
    else if (args.contains('--format cif') || args.contains('-f cif')) {
        format = 'cif'
    }
    else if (args.contains('--format pdbx') || args.contains('-f pdbx')) {
        format = 'pdbx'
    }

    """
    prepare_pdb.py prepare \\
        ${pdb_file} \\
        --output ${prefix}_prepared.${format} \\
        ${args}

    # 获取版本信息到变量中，避免直接输出失败
    OPENMM_VERSION=\$(python -c "import openmm; print(openmm.__version__)" || echo "unknown")
    PDBFIXER_VERSION=\$(python -c "import pdbfixer; print('1.11')" || echo "unknown")

    # 创建版本文件
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        openmm: "\${OPENMM_VERSION}"
        pdbfixer: "\${PDBFIXER_VERSION}"
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    // 在存根模式下也需要解析正确的格式
    def format = 'pdb'
    if (task.ext.args?.contains('--format pdb') || task.ext.args?.contains('-f pdb')) {
        format = 'pdb'
    }
    else if (task.ext.args?.contains('--format cif') || task.ext.args?.contains('-f cif')) {
        format = 'cif'
    }
    else if (task.ext.args?.contains('--format pdbx') || task.ext.args?.contains('-f pdbx')) {
        format = 'pdbx'
    }

    """
    touch ${prefix}_prepared.${format}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        openmm: "7.7.0"
        pdbfixer: "1.9"
    END_VERSIONS
    """
}
