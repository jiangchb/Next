import sys
if sys.platform.startswith("linux"):
    import matplotlib as mt

    mt.use('agg')
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.patches as patches

import pandas as pd
import numpy as np

##准备文件：富集top文件和对应的差异文件；目前还没有效支持miRNA转录本富集

def circos_plot(df,diff,typ,p,category,figname):
    ###数据准备
    dic = {'biological_process': 'darkorange', 'cellular_component': 'deepskyblue', 'molecular_function': 'yellowgreen', 'KEGG': 'darkorange'}
    karyotype_df = df.copy()
    karyotype_df['-log10(p-value)'] = -np.log10(karyotype_df[p])
    print (karyotype_df['-log10(p-value)'].head)
#    karyotype_df['-log10(p-value)'] = karyotype_df['-log10(p-value)'].replace(np.inf,sorted(set(karyotype_df['-log10(p-value)']))[-2]*1.1)
    karyotype_df['-log10(p-value)'] = round(karyotype_df['-log10(p-value)'],2)
    #karyotype_df['P'] = karyotype_df['-log10(p-value)']*1000/karyotype_df['-log10(p-value)'].max()
    karyotype_df['length'] = 1000
    all_width = 2 * np.pi - karyotype_df.shape[0] * 0.02
    #dic = {k:v for k,v in zip(sorted(list(set(df.category)),key=list(df.category).index),list(mcolors.cnames.keys())[:catogry])}
    if typ == 'GO':
        karyotype_df['color'] = karyotype_df[category].apply(lambda x:dic[x])
    else:
        karyotype_df['color'] = 'darkorange'
    
    karyotype_df['width'] = (karyotype_df['length']/ karyotype_df['length'].sum()) * all_width
    
    start_i = 0
    starts = []
    for length in karyotype_df['width']:
        starts.append(start_i)
        start_i = start_i + length + 0.02
    karyotype_df['P'] = karyotype_df['-log10(p-value)']*karyotype_df['width'][0]/karyotype_df['-log10(p-value)'].max()
    karyotype_df['num'] = karyotype_df['ListHits']*karyotype_df['width'][0]/karyotype_df['ListHits'].max()

    karyotype_df['start'] = starts
    karyotype_df['center'] = karyotype_df.apply(lambda x: x.start + x.width / 2, axis=1)
    karyotype_df['end'] = karyotype_df.apply(lambda x: x.start + x.width, axis=1)
    ###数据准备
# =============================================================================
#     scale = []
#     for index, row in karyotype_df.iterrows():
#         for i in range(row.length // 20000000 + 1):
#             scale.append([row.name, i * 20, row.width * i * 20000000 / row.length + row.start])
#     scale_df = pd.DataFrame(scale, columns=['chrom', 'scale_name', 'start'])
# =============================================================================
    #绘图调色
    map_vir = cm.get_cmap(name='Reds')
    
    #定义画布和坐标并去掉坐标的坐标轴信息
    fig = plt.figure(figsize=[16, 16])
    ax = fig.add_subplot(1, 1, 1, projection='polar')
    ax.axis('off')

    #通路圆圈
    ax.bar(x=karyotype_df['start'],
           height=0.36, bottom=6,
           width=karyotype_df['width'],
           color=karyotype_df['color'],
           align='edge',
           edgecolor='#A9A9A9',
           alpha=0.8)
    
    
    #通路加名称
    for index, row in karyotype_df.iterrows():
        ax.text(row.center, 6.2, row.name,
                ha='center', va='center', fontsize=11,
                rotation=np.rad2deg(
                    row.center if row.center < np.pi else row.center + np.pi) - 90)
    #配置p值颜色    
    print(karyotype_df)
    norm = plt.Normalize(min(karyotype_df['-log10(p-value)'])*0.6,max(karyotype_df['-log10(p-value)'])*1.2)
    norm_values = norm(karyotype_df['-log10(p-value)'].values)
    colors = map_vir(norm_values)
    
    #通路的p值
    ax.bar(x=karyotype_df['center']-karyotype_df['num']/2,
       height=0.25, bottom=5.6,
       width=karyotype_df['num'],
       color=colors,
       align='edge',
       edgecolor='',
       alpha=0.6)
    #通路对应的差异基因数
# =============================================================================
#     for index, row in karyotype_df.iterrows():
#         ax.text(row.center, 5.4, row['ListHits'],
#                 ha='center', va='center', fontsize=11,
#                 rotation=np.rad2deg(
#                     row.center if row.center < np.pi else row.center + np.pi) - 90)
# =============================================================================
    #ax.legend()
    #添加上下调信息
    karyotype_df['Up'] = karyotype_df.Gene.apply(lambda x:[diff[n] for n in x.split('; ')].count('Up'))
    karyotype_df['Down'] = karyotype_df.Gene.apply(lambda x:[diff[n] for n in x.split('; ')].count('Down'))
    karyotype_df['all'] = karyotype_df['Up'] + karyotype_df['Down'] 
    karyotype_df['Up1'] = karyotype_df['Up']*karyotype_df['width'][0]/karyotype_df['all']
    karyotype_df['Down1'] = karyotype_df['Down']*karyotype_df['width'][0]/karyotype_df['all']
    karyotype_df['Down2'] = karyotype_df['start']+karyotype_df['Up1']
    ax.bar(x=karyotype_df['start'],
           height=0.2, bottom=5.1,
           width=karyotype_df['Up1'],
           color='orangered',
           align='edge',
           edgecolor='#A9A9A9')
    for index, row in karyotype_df.iterrows():
        ax.text(row.start+row.Up1/2, 4.9, row.Up,
                ha='center', va='center', fontsize=15,color='orangered',
                rotation=np.rad2deg(
                    row.center if row.center < np.pi else row.center + np.pi) - 90)
        
    ax.bar(x=karyotype_df['Down2'],
           height=0.2, bottom=5.1,
           width=karyotype_df['Down1'],
           color='mediumseagreen',
           align='edge',
           edgecolor='#A9A9A9')  
    for index, row in karyotype_df.iterrows():
        ax.text(row.Down2+row.Down1/2, 4.9, row.Down,
                ha='center', va='center', fontsize=15,color='mediumseagreen',
                rotation=np.rad2deg(
                    row.center if row.center < np.pi else row.center + np.pi) - 90)
        
    #富集分数圆圈
    i,j=2.7,1
    while i >0:
        ax.bar(x=karyotype_df['start'],
           height=i, bottom=2,
           width=karyotype_df['width'],
           color='gray',
           align='edge',
           edgecolor='white',
           alpha=0.2*j)
        i = i-0.3
        j = j*0.5
    
    ax.bar(x=karyotype_df['start'],
           height=np.log10(karyotype_df['Enrichment_score']+1)*2.5/max(np.log10(karyotype_df['Enrichment_score']+1)), bottom=2,
           width=karyotype_df['width'],
           color=karyotype_df['color'],
           align='edge',
           edgecolor='#A9A9A9',
           alpha=0.8)
    

# =============================================================================
#     fig = plt.figure(figsize=[16, 16])
#     ax = fig.add_subplot(1, 1, 1, projection='polar')
# =============================================================================
    #添加右侧图例并分成上下两块区域
    legend_area = fig.add_axes([0.9, 0.15, 0.1, 0.75])
    legend_area.axis('off')
    legends_loc = [[0, 0.6, 0.4, 0.28],
                   [0, 0.3, 0.5, 0.28]]
    
    #图例1：GO分类
    if typ == 'GO':
        markers = sorted(list(set(df[category])),key=list(df[category]).index)
        legend_1 = legend_area.inset_axes(legends_loc[0])
        legend_1.axis('off')  # turn off all spines and ticks
        
        #marker_category = {v: k for k, v in category_marker.items()}
        for index, marker in enumerate(markers):
            legend_1.scatter(x=0, y=index*0.4, marker='s', color=dic[marker], edgecolors='face', s=160)
            legend_1.text(x=0.1, y=index*0.4, s=marker, va='center', ha='left', color=dic[marker],fontsize=24)
        legend_1.set_ylim(-0.5, index if index else 0.5)
        legend_1.set_xlim(-0.1, 0.1)
    
    #图例2：p值
    legend_2 = legend_area.inset_axes(legends_loc[1])
    legend_2.axis('off')  # turn off all spines and ticks    
    unique_p = karyotype_df['-log10(p-value)'].drop_duplicates().tolist() 
    if len(unique_p) >= 6:
        x2 = [0] * 6
        num2 = list(np.percentile(unique_p, (100, 80, 60, 40, 20, 0))) 
    else:
        x2 = [0] * len(unique_p)
        num2 = unique_p
    #s2 = [unique_p.index(n) for n in num2 ]
    y2 = [2 * i for i in range(0, len(x2))]
    map_vir = cm.get_cmap(name='Reds')
    norm = plt.Normalize(min(num2)*0.6,max(num2)*1.2)
    norm_values = norm(num2)
    colors = map_vir(norm_values)

    legend_2.scatter(x=x2, y=y2, s=500, c=colors, edgecolors='face')
    for i, value in enumerate(num2):
        legend_2.text(x=0.1, y=y2[i], s=round(value,2),
                      va='center', ha='left', fontsize=18)
    legend_2.set_ylim([-1, y2[-1]*1.1  if y2[-1] else 0.5])
    legend_2.set_xlim(-0.1, 0.1)
    legend_2.set_title('-log$_{10}$(p-value)', fontdict={'size': 20}, loc='left')
 
# =============================================================================
#     fig = plt.figure(figsize=[16, 16])
#     ax = fig.add_subplot(1, 1, 1, projection='polar')
# =============================================================================
    #添加圈内图例，分成四小块区域
    legend_area1 = fig.add_axes([0.42, 0.4, 0.24, 0.24])
    legend_area1.axis('off')
    legends_loc1 = [[0, 0.6, 0.2, 0.06],
                    [0, 0.45, 0.2, 0.06],
                    [0, 0.3, 0.2, 0.06],
                    [0, 0.15, 0.2, 0.06]]
    #添加基因数图例
    legend1_1 = legend_area1.inset_axes(legends_loc1[0])
    legend1_1.axis('off')
    legend1_1.bar([1,2,3],height=0.1, width=1.2,color=colors[3])
    legend1_1.text(x=4, y=0.05, s='Number of Genes', va='center', ha='left',fontsize=15)
    #添加上下调图例
    legend1_2 = legend_area1.inset_axes(legends_loc1[1])
    legend1_2.axis('off')
    legend1_2.bar([1,2,3],height=0.1, width=1.2,color='orangered')
    legend1_2.text(x=4, y=0.05, s='Up-Regulated', va='center', ha='left',fontsize=15)
    legend1_3 = legend_area1.inset_axes(legends_loc1[2])
    legend1_3.axis('off')
    legend1_3.bar([1,2,3],height=0.1, width=1.2,color='mediumseagreen')
    legend1_3.text(x=4, y=0.05, s='Down-Regulated', va='center', ha='left',fontsize=15)
    #添加富集分数图例
    legend1_4 = legend_area1.inset_axes(legends_loc1[3])
    legend1_4.axis('off')
    if typ == 'KEGG' or (typ == 'GO' and len(set(karyotype_df[category])) ==1):
        x0=[[0.1,0.95,0.95,0.1]]
        y0=[[0.25,0.05,0.95,0.75]]
        if typ == 'KEGG':
            z0=['darkorange']
        else:
            z0=[dic.get(karyotype_df[category][0],'darkorange')]
    elif len(set(karyotype_df[category])) ==2:
        x0=[[0.1,0.1+0.85/2,0.1+0.85/2,0.1],[0.1+0.85/2,0.95,0.95,0.1+0.85/2]]
        y0=[[0.25,0.15,0.85,0.75],[0.15,0.05,0.95,0.85]]
        z0=sorted([dic[n] for n in set(karyotype_df[category])])
    else:
        x0=[[0.1,0.1+0.85/3,0.1+0.85/3,0.1],[0.1+0.85/3,0.1+1.7/3,0.1+1.7/3,0.1+0.85/3],[0.1+1.7/3,0.95,0.95,0.1+1.7/3]]
        y0=[[0.25,0.25-0.2/3,0.75+0.2/3,0.75],[0.25-0.2/3,0.25-0.4/3,0.75+0.4/3,0.75+0.2/3],[0.25-0.4/3,0.05,0.95,0.75+0.4/3]]
        z0=['darkorange','deepskyblue','yellowgreen']
    for x,y,z in zip(x0,y0,z0):
        legend1_4.add_patch(patches.Polygon(xy=list(zip(x,y)), color=z))
    legend1_4.text(x=1.08, y=0.48, s='Enrichment-Score', va='center', ha='left',fontsize=15)
    plt.savefig('%s_circos.pdf' %figname,bbox_inches='tight')
    plt.savefig('%s_circos.png' %figname,bbox_inches='tight',dpi=300)

diff = pd.read_csv('BA-VS-CC-Protein-diff-expression.xls',sep='\t',index_col=0)['Regulation'].to_dict()
df = pd.read_csv('GO.qValue_order.enrichment.xls',sep='\t',index_col=0)
circos_plot(df,diff,'GO','pValue','Category','GO')  

#df = pd.read_csv('KEGG.top.Total.xls',sep='\t',index_col=0) 
#circos_plot(df,diff,'KEGG','pValue','Category','KEGG')  
    
    
    
    


