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
    
    // *直接使用 conda 安装环境
    // conda "conda-forge::python=3.12 conda-forge::plotly=5.13.1 conda-forge::pandas=1.5.3"
    // *Environment Modules
    // module 'ncbi-blast/2.2.27'

    // 使用 Singularity 且不直接拉取 Docker 镜像时，需要使用 docker:// 前缀
    // container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'docker://august777/radiomics-plastimatch:1.4.7' :
    //     'august777/radiomics-plastimatch:1.4.7'}"

    // 3. 输出文件发布路径
    // 需要放在输入输出之前，conf/modules.config 中定义了输出文件的默认发布路径
    // publishDir "${params.outdir}", mode: 'copy',pattern: "**.txt"
   
    // 4. 进程输入输出定义
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
