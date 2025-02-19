//
// 具有特定于 open-bio/molflow 管道的功能的子工作流
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    导入功能/模块/子工作流程
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap } from 'plugin/nf-schema'
include { samplesheetToList } from 'plugin/nf-schema'
include { completionEmail } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    初始化管道的子工作流
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {
    take:
    version // boolean: 显示版本并退出
    validate_params // boolean：布尔值是否在运行时根据模式验证参数
    monochrome_logs // 布尔值：不使用彩色日志输出
    nextflow_cli_args //   array：位置 nextflow CLI 参数列表
    outdir //  string：保存结果的输出目录
    input //  string：输入样本表的路径

    main:

    ch_versions = Channel.empty()

    //
    // 打印版本并退出（如果需要）并将管道参数转储到 JSON 文件
    //
    UTILS_NEXTFLOW_PIPELINE(
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1,
    )

    //
    // 验证参数并生成参数摘要到标准输出
    //
    UTILS_NFSCHEMA_PLUGIN(
        workflow,
        validate_params,
        null,
    )

    //
    // 检查提供给管道的配置
    //
    UTILS_NFCORE_PIPELINE(
        nextflow_cli_args
    )


    //
    // 从通过 params.input 提供的输入文件创建通道
    //

    ch_samplesheet = Channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        .map { meta, pdb, mol ->
            // 如果是地址，将其转换为文件
            pdb = pdb instanceof File ? pdb : file(pdb)
            mol = mol instanceof File ? mol : file(mol)

            return [meta, pdb, mol]
        }

    emit:
    samplesheet = ch_samplesheet
    versions = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    管道完成的子工作流程
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {
    take:
    email //  字符串：电子邮件地址
    email_on_fail //  字符串：管道故障时发送的电子邮件地址
    plaintext_email // 布尔值：发送纯文本电子邮件而不是 HTML
    outdir //    路径：将发布结果的输出目录的路径
    monochrome_logs // 布尔值：在日志输出中禁用 ANSI 颜色代码
    hook_url //  string：通知的钩子 URL
    multiqc_report //  字符串：MultiQC 报告的路径

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

    //
    // 完成电子邮件和摘要
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_report.toList(),
            )
        }

        completionSummary(monochrome_logs)
        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error("Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting")
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    功能
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//

//

// 生成 MultiQC 的方法描述
//
def toolCitationText() {
    // TODO nf-core：可以选择将文本引用工具添加到此列表中。
    // 可以使用三元运算符动态构造基础条件，例如参数[“run_xyz”]？ “工具（Foo 等人，2023）”：“”，
    // 取消注释 methodDescriptionText 中的函数以在 MultiQC 报告中呈现
    def citation_text = [
        "Tools used in the workflow included:",
        "FastQC (Andrews 2010),",
        "MultiQC (Ewels et al. 2016)",
        ".",
    ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core：可以选择将书目条目添加到此列表中。
    // 可以使用三元运算符动态构造基础条件，例如参数[“run_xyz”]？ "<li>作者 (2023) 出版商名称、期刊、DOI</li>" : "",
    // 取消注释 methodDescriptionText 中的函数以在 MultiQC 报告中呈现
    def reference_text = [
        "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
        "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>",
    ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // 转换为命名映射，以便可以与 MultiQC YML 文件中熟悉的 NXF ${workflow} 变量语法一起使用
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // 管道 DOI
    if (meta.manifest_map.doi) {
        // 使用循环处理多个 DOI
        // 删除“https://doi.org/”以使用 DOI 与 DOI 解析器处理管道
        // 删除``，因为manifest.doi是一个字符串而不是一个正确的列表
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    }
    else {
        meta["doi_text"] = ""
    }
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // 工具参考
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core：仅当 toolCitationText/toolBibliographyText 中的逻辑已填充时才取消下面的注释！
    // meta["tool_itations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", “。”）
    // 元["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine = new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
