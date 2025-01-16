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
    
    // 2.1 conda yaml 配置文件，好处是可以优雅包含pip包
    conda "${moduleDir}/environment.yml"
    
    // 2.2 直接使用 conda 安装环境
    // conda "conda-forge::python=3.12 conda-forge::plotly=5.13.1 conda-forge::pandas=1.5.3"
    
    // 2.3 Environment Modules
    // module 'ncbi-blast/2.2.27'

    // 2.4 docker和sif容器
    // 使用 Singularity 且不直接拉取 Docker 镜像时，需要使用 docker:// 前缀
    // container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'docker://august777/radiomics-plastimatch:1.4.7' :
    //     'august777/radiomics-plastimatch:1.4.7'}"

    // 3. 输出文件发布路径
    // 定义需要放在输入输出之前，conf/modules.config 中定义了输出文件的默认发布路径：切分进程名取第一节
    // publishDir "${params.outdir}", mode: 'copy',pattern: "**.txt"
   
    // 4. 进程输入输出定义
    // 注意：显式说明的路径才会在容器中自动挂载
    input:
    tuple val(meta), path(pdb_file, stageAs: "pdb/*"), path(sdf_file, stageAs: "sdf/*")

    output:
    // tuple val(meta), path("*.html"), emit: html 
    // tuple val(meta), path("*.zip") , emit: zip
    // path  "versions.yml"           , emit: versions
    tuple val(meta), path("${meta.id}_module.txt"), emit: module
    path "versions.yml", emit: versions
    
    // 5. 执行条件控制
    when:
    //task.ext 在 conf/modules.config 中定义 例：
    //withName: FASTQC { // 如果模块名称是 "FASTQC"
    //     ext.args = '--quiet' // 添加一个额外参数 `--quiet`
    // }
    task.ext.when == null || task.ext.when 
    
    // 6. 主要脚本逻辑
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
    // task 就是指当前运行的进程实例的配置对象
    // task.process 为进程名称，可用于版本信息的记录
    // task.memory 为进程内存，task.cpus 为进程 CPU 数量

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
