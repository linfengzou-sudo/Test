import pandas as pd
import matplotlib.pyplot as plt

def read_fai(fai_file):
    """
    读取fai文件，返回染色体及长度的DataFrame
    """
    df = pd.read_csv(fai_file, sep='\t', header=None, usecols=[0,1], names=['chr','length'])
    return df

def calc_cumstart(chrom_sizes):
    """
    计算染色体累积起始位置
    """
    chrom_sizes = chrom_sizes.sort_values('chr').reset_index(drop=True)
    chrom_sizes['cum_start'] = chrom_sizes['length'].cumsum() - chrom_sizes['length']
    return chrom_sizes

def add_cumpos(df_snp, chrom_sizes):
    """
    根据染色体累积起始位置计算SNP累积坐标
    """
    chr_to_cumstart = dict(zip(chrom_sizes['chr'], chrom_sizes['cum_start']))
    df_snp['cum_pos'] = df_snp.apply(lambda r: chr_to_cumstart.get(r['chr'], 0) + r['pos'], axis=1)
    return df_snp

def plot_mapq(dfA, chrom_sizes_A, dfB, chrom_sizes_B, output_png):
    """
    绘制两个参考基因组的SNP平均MAPQ散点图
    """
    plt.figure(figsize=(16,6))

    plt.scatter(dfA['cum_pos'], dfA['avg_mapq'], c='blue', s=10, alpha=0.6, label='Ref V2')
    plt.scatter(dfB['cum_pos'], dfB['avg_mapq'], c='red', s=10, alpha=0.6, label='Ref V5')

    plt.xlabel('Genomic position (cumulative by chromosome)')
    plt.ylabel('Average MAPQ')
    plt.title('Average MAPQ per SNP (cumulative genomic coordinate)')

    # 画染色体边界竖线 - 参考基因组A
    for start in chrom_sizes_A['cum_start']:
        plt.axvline(x=start, color='blue', linestyle='--', linewidth=0.5, alpha=0.3)
    # 参考基因组B边界竖线
    for start in chrom_sizes_B['cum_start']:
        plt.axvline(x=start, color='red', linestyle='--', linewidth=0.5, alpha=0.3)

    # 染色体标签位置：分别在两个基因组的中点标注
    def add_chr_labels(chrom_sizes, y_pos, color):
        chrom_mid = chrom_sizes['cum_start'] + chrom_sizes['length'] / 2
        for x, chrname in zip(chrom_mid, chrom_sizes['chr']):
            plt.text(x, y_pos, chrname, color=color, horizontalalignment='center', fontsize=9, alpha=0.7)

    ylim = plt.ylim()
    add_chr_labels(chrom_sizes_A, ylim[1]*0.95, 'blue')
    add_chr_labels(chrom_sizes_B, ylim[1]*0.90, 'red')

    plt.legend()
    plt.tight_layout()
    plt.savefig(output_png, dpi=300)
    print(f"Plot saved as {output_png}")

def main():
    # 文件路径，改成你的文件名
    refA_snp_file = "SNP_MAPQ_avg_V2.txt"
    refB_snp_file = "SNP_MAPQ_avg_V5.txt"
    refA_fai = "v2.fa.fai"
    refB_fai = "v5.fa.fai"
    output_png = "SNP_avg_MAPQ_cumulative_scatter.png"

    # 读入数据，指定列名（无表头）
    dfA = pd.read_csv(refA_snp_file, sep="\t", header=None, names=['chr','pos','avg_mapq'])
    dfB = pd.read_csv(refB_snp_file, sep="\t", header=None, names=['chr','pos','avg_mapq'])

    chrom_sizes_A = read_fai(refA_fai)
    chrom_sizes_B = read_fai(refB_fai)

    chrom_sizes_A = calc_cumstart(chrom_sizes_A)
    chrom_sizes_B = calc_cumstart(chrom_sizes_B)

    dfA = add_cumpos(dfA, chrom_sizes_A)
    dfB = add_cumpos(dfB, chrom_sizes_B)

    plot_mapq(dfA, chrom_sizes_A, dfB, chrom_sizes_B, output_png)

if __name__ == "__main__":
    main()

