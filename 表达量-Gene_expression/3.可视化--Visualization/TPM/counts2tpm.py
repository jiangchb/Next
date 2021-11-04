import pandas as pd

import click

@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option("-c","--counts",type=str,help='counts file')
@click.option("-l","--lens", type=str,help='gene(all-exon) length file')
@click.option("-o","--output", type=str,help='fpkm file name')
def counts2tpm(counts, lens, output):
    counts = pd.read_csv(counts, sep='\t', index_col=0)
    anno = pd.read_csv(lens, sep='\t', index_col=0)
    anno.columns = ['len']
    d = anno.to_dict()['len']
    lens = pd.DataFrame(d, index=counts.columns).T
    f1 = counts*1000/lens
    sum_reads = f1.sum()
    result = f1*1000000/sum_reads
    result.index.name='id'
    result.to_csv(output, sep='\t')

if __name__=='__main__':
    counts2tpm()
