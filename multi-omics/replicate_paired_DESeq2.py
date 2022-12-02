import pandas as pd
from rpy2 import robjects
from rpy2.robjects import pandas2ri, Formula
from rpy2.robjects.packages import importr
import click
import numpy as np

def cal_BaseMeanAB(counts,group):
    df1 =  np.log(counts)
    df2 = df1.mean(1).to_frame()
    df2.columns = ["Average of log values"]
    df3 = df2.replace(-np.inf, np.nan).dropna()
    df4 = df1.sub(df1.mean(1), axis=0).dropna()
    df5 = df4.median().to_frame().T
    df6 = np.e**df5
    df7 = counts.div(df6.iloc[0, :])

    case_name = group["Group"].tolist()[0]
    control_name = group["Group"].tolist()[-1]

    dic = group["Group"].to_dict()
    df8 = df7.groupby(dic, axis=1).mean()
    df8 = df8[[control_name,case_name]]
    df8.rename(columns={df8.columns[0]:f"BaseMean_control_{control_name}"},inplace=True)
    df8.rename(columns={df8.columns[1]:f"BaseMean_case_{case_name}"},inplace=True)
    df8.index.name="gene_id"
    return df8,f"BaseMean_case_{case_name}",f"BaseMean_control_{control_name}"


@click.command()
@click.option("-c", "--counts", help="result/diff/{group}/{group}_gene_counts.xls")
@click.option("-s", "--sample_group", help="result/diff/{group}/{group}_sample_group.txt")
@click.option("-o", "--outfile", help="result/diff/{group}/{group}_diff_unfiltered.xls")
def run_replicate_paired_DESeq2(counts, sample_group,outfile):
    # 必须加，否则转换数据类型会报错
    pandas2ri.activate()
    to_dataframe = robjects.r('function(x) data.frame(x)')

    # pandas读入文件
    countData_py = pd.read_csv(counts, sep='\t', index_col=0)
    colData_py = pd.read_csv(sample_group, sep='\t', index_col=0)
    BaseMeanAB,case_name,control_name = cal_BaseMeanAB(countData_py, colData_py)
    countData_py.index.name = None
    colData_py.index.name = None

    # 为了DESeq2.results(contrast=)
    casename = colData_py["Group"].tolist()[0]
    controlname = colData_py["Group"].tolist()[-1]

    # 转换
    countData_r = pandas2ri.py2rpy(countData_py)
    colData_r = pandas2ri.py2rpy(colData_py)

    #py转r的过程中会把-自动替换成.,需要转换回去
    countData_r.colnames = [i.replace(".", "-") for i in countData_r.colnames]
    colData_r.rownames = [i.replace(".", "-") for i in colData_r.rownames]

    DESeq2 = importr('DESeq2')

    dds = DESeq2.DESeqDataSetFromMatrix(countData=countData_r, colData=colData_r, design=Formula("~ Type + Group"))
    dds = DESeq2.DESeq(dds)
    res_r = DESeq2.results(dds, independentFiltering=False, cooksCutoff=False,contrast=robjects.StrVector(["Group",casename,controlname]))

    res_r = to_dataframe(res_r)
    res_py = robjects.conversion.rpy2py(res_r)
    res_py.index.name = "gene_id"

    result = pd.merge(res_py, BaseMeanAB, how="left", left_index=True, right_index=True)
    result.insert(1,control_name,result.pop(control_name))
    result.insert(2,case_name,result.pop(case_name))

    result.to_csv(outfile,sep='\t')


if __name__ == '__main__':
    run_replicate_paired_DESeq2()