import os

RESULTS_FOLDERNAME = "/home/shpc_100839/PaperCD8/data/Tex/scenic/res/"

SAMPLE_ID = 'Tex'

EXP_FILTER_LOOM_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.h5ad.filter.loom'.format(SAMPLE_ID))
MOTIFS_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.motifs.csv'.format(SAMPLE_ID))
AUCELL_LOOM_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.aucell.csv'.format(SAMPLE_ID))

os.system("pyscenic aucell %s %s \
    --output %s \
    --num_workers 8" %(EXP_FILTER_LOOM_FNAME,MOTIFS_FNAME,AUCELL_LOOM_FNAME))
