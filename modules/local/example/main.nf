process TEST_BASE {

    // 1. 进程标识和资源管理
    tag "${meta.id}"
    // 进程标签，用于日志和报告唯一标识
    label 'process_low'
    // 资源标签, 在 conf/base.config 中定义

    errorStrategy 'retry'
    // 错误处理策略
    maxRetries 3
    // 最大重试次数

    // 2. 环境和容器管理
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(pdb_file, stageAs: "pdb/*"), path(sdf_file, stageAs: "sdf/*")

    output:
    // tuple val(meta), path("*.html"), emit: html 
    // tuple val(meta), path("*.zip") , emit: zip
    // path  "versions.yml"           , emit: versions
    tuple val(meta), path("${meta.id}_module.txt"), emit: module
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // 可选参数
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo ${meta.id}
    echo ${pdb_file}
    echo ${sdf_file}
    
    example.py -f  ${prefix}_module.txt ${args}    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version| grep -oP '(?<=Python ).*')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_module.txt

     cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version| grep -oP '(?<=Python ).*')
    END_VERSIONS
    """
}
