/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
/* IMPORT 语句 */
// 导入格式: include { 组件名 } from '路径' 路径可以到main.nf也可以到它的父目录

include { TEST_BASE } from '../modules/local/example/main'
include { TEST_PLOTLY } from '../modules/local/example_plotly'
include { MULTIQC } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap } from 'plugin/nf-schema'
include { paramsSummaryMultiqc } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_molflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MOLFLOW {
    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    //
    // MODULE: Run TEST_MODULE
    //
    TEST_BASE(
        ch_samplesheet
    )

    TEST_PLOTLY()

    // mix() 合并通道示例:
    // 假设 TEST_MODULE.out.module 输出为: ['test1', 'result1.txt']
    //                                    ['test2', 'result2.txt']
    // collect { it[1] } 提取第二个元素: ['result1.txt', 'result2.txt']
    // mix 操作后: ch_multiqc_files 现在包含: ['result1.txt', 'result2.txt']
    ch_multiqc_files = ch_multiqc_files.mix(
        TEST_BASE.out.module.collect { it[1] },
        TEST_PLOTLY.out.plot.collect(),
    )

    // first() 获取第一个元素
    // 假设 TEST_MODULE.out.versions 输出: ['v1.0', 'v2.0']
    // first() 后获得: 'v1.0'
    ch_versions = ch_versions.mix(
        TEST_BASE.out.versions,
        TEST_PLOTLY.out.versions,
    )

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: '' + 'pipeline_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )
        .set { ch_collated_versions }
    // 生成的YAML文件内容:
    // test_module: v1.2.3
    // xxx: v 
    // xxx: v 

    //
    // MODULE: MultiQC
    //
    // 读取默认的MultiQC配置文件，定义报告的基本布局和设置
    ch_multiqc_config = Channel.fromPath(
        "${projectDir}/assets/multiqc_config.yml",
        checkIfExists: true
    )

    // 允许用户提供自定义配置文件，覆盖默认设置
    ch_multiqc_custom_config = params.multiqc_config
        ? Channel.fromPath(params.multiqc_config, checkIfExists: true)
        : Channel.empty()

    // 允许用户提供自定义log，覆盖默认设置
    ch_multiqc_logo = params.multiqc_logo
        ? Channel.fromPath(params.multiqc_logo, checkIfExists: true)
        : Channel.empty()

    // 使用 paramsSummaryMap 函数生成工作流参数摘要
    // 参数包括 workflow 对象和 schema 文件路径
    summary_params = paramsSummaryMap(
        workflow,
        parameters_schema: "nextflow_schema.json"
    )

    // 创建一个值为参数摘要的通道,用于 MultiQC 报告
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

    // 将工作流摘要文件混入 multiqc 文件通道
    // collectFile 操作符将内容收集到指定名称的文件中
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml')
    )

    // 三元运算符检查是否提供了自定义方法描述文件
    // 如果没有则使用默认模板文件
    ch_multiqc_custom_methods_description = params.multiqc_methods_description
        ? file(params.multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)

    // 创建包含方法描述文本的通道
    ch_methods_description = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description)
    )

    // 将版本信息文件混入 multiqc 文件通道
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)

    // 将方法描述文件混入 multiqc 文件通道
    // sort:true 表示对内容进行排序
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true,
        )
    )

    // 调用 MULTIQC 流程
    // collect() 和 toList() 操作符将通道内容转换为列表
    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        [],
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // 输出 MultiQC HTML 报告路径
    versions = ch_versions // 输出软件版本信息文件路径
}
