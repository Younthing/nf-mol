process TEST_PLOTLY {
    tag "Plotly plot"

    conda "conda-forge::python=3.10 conda-forge::plotly=5.13.1 conda-forge::pandas=1.5.3"

    output:
    path "scatter_plot_mqc.html", emit: plot
    // 改为 _mqc.html 后缀
    path "versions.yml", emit: versions

    script:
    """
    example_plotly.py  # 生成 scatter_plot_mqc.html 文件

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plotly: \$(python -c "import plotly; print(plotly.__version__)")
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    
    """
}
