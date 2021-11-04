import argparse
import os

def GSVA_analysis(i, k, g, F, p, o, c, h, w, n, R, B,  G):
    path = os.path.dirname(os.path.abspath(__file__))
    comparison = ""
    if c:
        group = c.split(",")
        comparison = "-vs-".join([group[0],group[2]])
    #Complete GSVA analysis
    command = R + " " + path+"/GSVA_diffexp.R"
    if not i:
        return "Please input a file"
    else:
        command = " ".join(str(s) for s in [command,"-i",i, '-o', o])
    #GSVA
    if g:
        if c:
            print(" ".join(str(s) for s in [command,'-r GO_GSVA_enrichment.xls', '-g', g, '-c', c, '-A GO_allexp_GSVA_' + comparison + '.xls', '-D GO_diffexp_GSVA_' + comparison + '.xls', "-F", F, '-p', p,'-B GO_GSVA_score_' + comparison, '-H', h, '-W', w, '-N', n, '-G', G]))
            os.system(" ".join(str(s) for s in [command,'-r GO_GSVA_enrichment.xls', '-g', g, '-c', c, '-A GO_allexp_GSVA_' + comparison + '.xls', '-D GO_diffexp_GSVA_' + comparison + '.xls', "-F", F, '-p', p,'-B GO_GSVA_score_' + comparison, '-H', h, '-W', w, '-N', n, '-G', G]))
        else:
            print(" ".join(str(s) for s in [command, '-g', g, '-r GO_GSVA_enrichment.xls']))
            os.system(" ".join(str(s) for s in [command, '-g', g, '-r GO_GSVA_enrichment.xls']))
    if k:
        if c:
            print(" ".join(str(s) for s in [command, '-r KEGG_GSVA_enrichment.xls', '-g', k, '-c', c, '-A KEGG_allexp_GSVA_' + comparison + '.xls', '-D KEGG_diffexp_GSVA_' + comparison + '.xls', "-F", F, '-p', p, '-B KEGG_GSVA_score_' + comparison, '-H', h, '-W', w, '-N', n, '-G', G]))
            os.system(" ".join(str(s) for s in [command, '-r KEGG_GSVA_enrichment.xls', '-g', k, '-c', c, '-A KEGG_allexp_GSVA_' + comparison + '.xls', '-D KEGG_diffexp_GSVA_' + comparison + '.xls', "-F", F, '-p', p, '-B KEGG_GSVA_score_' + comparison, '-H', h, '-W', w, '-N', n, '-G', G]))
        else:
            print(" ".join(str(s) for s in [command, '-g', k, '-r KEGG_GSVA_enrichment.xls']))
            os.system(" ".join(str(s) for s in [command, '-g', k, '-r KEGG_GSVA_enrichment.xls']))

    #Differential analysis
    if not g and not k and c:
        print(" ".join(str(s) for s in [command, '-c', c, '-A allexp_GSVA_' + comparison + '.xls', '-D diffexp_GSVA_' + comparison + '.xls','-B diffexp_GSVA_score_' + comparison, '-R reg', '-H', h, '-W', w, '-N', n, '-G', G]))
        os.system(" ".join(str(s) for s in [command, '-c', c,'-A allexp_GSVA_' + comparison + '.xls', '-D diffexp_GSVA_' + comparison + '.xls', "-F", F, '-p', p, '-B diffexp_GSVA_score_' + comparison, '-R reg', '-H', h, '-W', w, '-N', n, '-G', G]))

    #Bar plot:
    if not g and not k and not c:
        print(" ".join(str(s) for s in [command, '-B', B, '-H', h, '-W', w, '-N', n, '-G', G]))
        os.system(" ".join(str(s) for s in [command, '-B', B, "-F", F, '-p', p, '-H', h, '-W', w, '-N', n, '-G', G]))


def run():
    parser = argparse.ArgumentParser("Run gsva for epxression set. Caution: ")
    parser.add_argument('-i', '--input', help = "Input an expression set such as fpkm.xls; or gsva enrichement set to do a differentiation; or differential result to draw bar plot. The former samples should belong to case group, and the latters should belong to control group")
    parser.add_argument('-k', '--kegg_gmt', help = "Input a kegg gmt file")
    parser.add_argument('-g', '--go_gmt', help = "Input a go gmt file")
    parser.add_argument('-F', '--FDR', default = '1', help = "The thresholds of FDR, default None")
    parser.add_argument('-p', '--pval', default ='0.05', help = "Thresholds of p value, default 0.05")
    parser.add_argument('-o', '--OUTDIR', default = os.getcwd(), help = "Output directory, default currect directory")
    parser.add_argument('-c', '--group', help = "The group names and sample number. For example A,4,B,4, namely the former four samples belong to group A, and the latters belong to group B")
    parser.add_argument('-H','--height', default = '40', help = "The height of plot")
    parser.add_argument('-W','--width', default = '20', help="The width of plot")
    parser.add_argument('-R', '--Rscript', default='/home/shenyitian/miniconda3/bin/Rscript', help="The interpreter for R script")
    parser.add_argument('-B', '--bar_plot', default= "diffexp_genesets_GSVA_score", help = "The name of bar plot, only work when merely darwing a plot")
    parser.add_argument('-N', '--topn', default= "100", help = "The number of bars, ranked by FDR if filtering by FDR, or similarily by pval")
    parser.add_argument('-G', '--sig_reg', default= "reg", help = "Color by significance(pval or FDR) or regulation(Up down), sig or reg")
    args = parser.parse_args()
    GSVA_analysis(i = args.input, k = args.kegg_gmt, g = args.go_gmt, F = args.FDR, p = args.pval, o = args.OUTDIR, c = args.group, h = args.height, w = args.width, R = args.Rscript, B = args.bar_plot, n = args.topn, G = args.sig_reg)
if __name__ == '__main__':
    run()
