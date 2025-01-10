include { TEST_MODULE } from './modules/local/example'

workflow {
    // 打印所有参数
    println "\nWorkflow parameters:"
    println "===================="
    params.each { key, value ->
        println "Parameter: $key -> $value"
    }

    // 从CSV读取或创建默认测试数据
    def input_data = params.input ? 
        Channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> 
                def meta = [id: row.id]
                def pdb = file(row.pdb)
                [meta, pdb]
            } :
        Channel.of(
            [[id: 'test_sample'], file("${projectDir}/test-data")]
        )
    
    // 添加日志记录
    input_data.subscribe { meta, pdb -> 
        println "Processing: ID=${meta.id}, PDB=${pdb}"
    }
    
    def result = TEST_MODULE(input_data)

    result.module.subscribe { println "Module output file: ${it}" }
    result.versions.subscribe { println "Versions output file: ${it}" }
}