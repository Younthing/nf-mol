#!/usr/bin/env python3
import plotly.express as px
import numpy as np
import pandas as pd

# 生成随机数据
np.random.seed(42)
n_points = 100
df = pd.DataFrame(
    {"x": np.random.normal(0, 1, n_points), "y": np.random.normal(0, 1, n_points)}
)

# 创建散点图
fig = px.scatter(df, x="x", y="y", title="Random Scatter Plot")

# 设置图表的大小和边距
fig.update_layout(width=800, height=600, margin=dict(l=50, r=50, t=50, b=50))

# 保存为 HTML，设置全屏按钮和 ModeBar
fig.write_html(
    "scatter_plot_mqc.html",
    include_plotlyjs="cdn",
    config={"displayModeBar": True, "responsive": True},
)
