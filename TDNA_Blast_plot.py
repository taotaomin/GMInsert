import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math
import argparse

def draw_figure(ref, result, o):
    # 初始化变量，用于存储碱基数量
    tdna_length = 0
    # 读取FASTA文件
    with open(ref, 'r') as file:
        # 跳过注释行
        header = file.readline()
        # 读取碱基数量
        for line in file:
            tdna_length += len(line.strip())

    # 读取txt文件的1、9、10行
    df = pd.read_csv(result, header=None, usecols=[0, 8, 9], sep="\t")
    df.columns = ['id', 'start', 'end']
    df['start_position'] = df[['start', 'end']].min(axis=1)
    df['length'] = (df['end'] - df['start']).abs()
    df['color'] = np.where(df['end'] > df['start'], 'green', 'red')

    # 提取需要用到的id
    df['id'] = df['id'].apply(lambda x: "_".join(x.split('_')[1:4]))

    # 设置图形大小
    plt.figure(figsize=(10, len(df) * 0.6))

    # 向上取整到最接近的 1000 的倍数
    rounded_length = math.ceil(tdna_length / 1000) * 1000
    plt.xlim(0, rounded_length)

    # 根据起始位置和长度绘制条形图
    plt.barh(df.index, df['length'], color=df['color'], left=df['start_position'])
    # 在每个条形图的坐标添加ID
    for index, row in df.iterrows():
        plt.text(row['start_position'] - 500, index, row['id'], va='center', ha='right')

    # 根据颜色添加位置信息
    for index, row in df.iterrows():
        position_info = f"{row['start']} - {row['end']}" if row['color'] == 'green' else f"{row['end']} - {row['start']}"
        plt.text(row['start_position'] + row['length'] + 500, index, position_info, va='center', ha='left')

    # 画TDNA图
    plt.barh(-1, tdna_length, color='lightblue')
    plt.text(tdna_length / 2, -1, f'TDNA Length: {tdna_length}', va='center', ha='center', color='black', label='TDNA')
    plt.title('assemble_TDNA_blast', y=1.1)
    plt.gca().invert_yaxis()
    # 设置标签
    for color in df['color'].unique():
        subset = df[df['color'] == color]
        # 基于颜色选择标签
        if color == 'green':
            label = 'Forward alignment'  # 绿色代表正向
        elif color == 'red':
            label = 'Reverse alignment'  # 红色代表反向
        else:
            label = 'Other Strand'  # 其他颜色可以自行定义标签

        # Plot only the first bar for each color with the label
        plt.barh(subset.index[0], subset['length'].values[0], color=color, left=subset['start_position'].values[0],
                 label=label)
    plt.legend()

    # plt.xlabel('Position') 未添加
    # 隐藏y轴坐标
    plt.yticks([])
    # 去除框线
    plt.box(False)
    # 将X轴放到图形的最上方
    plt.gca().xaxis.tick_top()
    # 显示图形
    plt.savefig(o)


def parse_args():
    parser = argparse.ArgumentParser(description='Argument Parser for file paths and save paths')
    parser.add_argument('--ref', help='File path for reference file')
    parser.add_argument('--result', help='File path for result file')
    parser.add_argument('--o', help='File path to save the generated plot')
    args = parser.parse_args()

    # Check if the specified output file has a PDF extension
    if not args.o.lower().endswith('.pdf'):
        parser.error("The output file must have a '.pdf' extension.")

def main():
    args=parse_args()
    draw_figure(args.ref, args.result, args.o)

if __name__ == "__main__":
    main()