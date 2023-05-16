import pandas as pd
import numpy as np
import scipy as sp
from rpy2.robjects import r
from rpy2.robjects.conversion import localconverter
from scipy.stats import ranksums
from sklearn.neighbors import kneighbors_graph
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.metrics import roc_curve, auc
from sklearn.decomposition import PCA
from sklearn.metrics import average_precision_score
import matplotlib.pyplot as plt
import scanpy as sc
import signature_utils as su

def random_siganture(exp, graph, size=300):
    genelist = exp.sample(size).index
    sig = su.signature_values(exp, genelist)
    return su.propagation(sig, graph), list(genelist)


def random_genes(exp, size=30):
    return exp.sample(size).index


def calculate_roc_auc(idents, predict):
    fpr, tpr, _ = roc_curve(idents, predict, pos_label=1)
    return auc(fpr, tpr)

def calculate_aupr(idents, predict):
    return average_precision_score(idents, predict)

def plot_aucs(auc, min_s, step, max_s):
    plt.plot(range(min_s, max_s, step), auc)
    plt.title("ROC AUC for sub signatures")
    plt.xlabel("size of signature")
    plt.ylabel("ROC AUC")
    plt.show()

def build_snn_graph(pca,k=20):
    #not good enogh 
    graph = kneighbors_graph(pca, k, mode='connectivity', include_self=True).toarray()
    graph = pd.DataFrame(graph)
    return graph  


def reduce_signature(exp, sig, graph, idents, down=None, min_s=10, step=10, max_s=None):
    aucs = []
    if max_s is None:
        if down is not None:
            max_s = min(len(sig), len(down))
        max_s = len(sig)

    sig = pd.DataFrame(sig).sample(frac=1, random_state=42)
    for size in range(min_s, max_s, step):
        sig_scores = su.signature_values(exp, sig[:size], down_sig=down[:size])
        pred = su.propagation(sig_scores, graph)
        aucs.append(calculate_roc_auc(idents, pred))
    plot_aucs(aucs, min_s, step, max_s)
    return aucs


def reduce_genes(exp, sig, idents, down=None, min_s=1000, step=1000):
    aucs = []
    exp = exp.sample(frac=1, random_state=42)
    for size in range(min_s, exp.shape[0], step):
        exp_reduce = exp.iloc[:size]
        exp_reduce = exp_reduce.apply(lambda x: np.log1p(x / x.sum()) * 1000, axis=0)
        pca_func = PCA(n_components=20)
        pca = pca_func.fit_transform(exp_reduce.T)
        graph = build_snn_graph(pca)
        sig_scores = su.signature_values(exp_reduce, sig, down_sig=down)
        pred = su.propagation(sig_scores, graph)
        aucs.append(calculate_roc_auc(idents, pred))
    plt.plot(range(min_s, exp.shape[0], step), aucs)
    plt.title("ROC AUC for subset of genes")
    plt.xlabel("Number of genes")
    plt.ylabel("ROC AUC")
    plt.show()
    return aucs


def compare_multy_signature(exp, graph, sigs, idents, names, downs=None, min_s=10, step=10, max_s=None):
    aucs = []
    if max_s is None:
        if downs is not None:
            max_s = min(list(map(lambda x: len(x), sigs)) + list(map(lambda x: len(x), downs)))
        else:
            max_s = min(map(lambda x: len(x), sigs))

    for i in range(len(sigs)):
        if downs is not None:
            aucs.append(reduce_signature(exp, sigs[i], graph, idents[i], downs[i], min_s, step, max_s))
        else:
            aucs.append(reduce_signature(exp, sigs[i], graph, idents[i], min_s=min_s, step=step, max_s=max_s))

    for i in range(len(aucs)):
        plt.plot(range(min_s, max_s, step), aucs[i], label=names[i])
        plt.title("ROC AUC for sub signatures")
        plt.xlabel("size of signature")
        plt.ylabel("ROC AUC")
        plt.legend()
    plt.show()
    return aucs

def make_result_barplot(d,p,cluster,ROC_flag):
    plt.bar(["Scanpy","Propagation"],[d,p],color=["coral","aquamarine"])    
    if ROC_flag:
        plt.ylabel('ROCAUC')
    else:
        plt.ylabel('AUPR')
    plt.title(f"Preformance of Signature Classifaction in {cluster}")

plt.show()


def test_scanpy_obj_roc(obj,sig_up, ident_name,ident_col,cluster_name,conn_flag = True,sig_down=None, ROC_flag=True,umap_flag = False):
    targets = obj.obs[ident_col] == ident_name
    if sig_down is not None:
        sc.tl.score_genes(obj,sig_up,gene_pool=list(sig_down)+list(sig_up),use_raw=False)
    sc.tl.score_genes(obj,sig_up,use_raw=False)
    if umap_flag:
        sc.pl.umap(obj, color=["score"],color_map="magma")
    else:
        sc.pl.tsne(obj, color=["score"],color_map="magma")
    obj = su.run_signature_on_obj(obj,sig_up,sig_down ,conn_flag=conn_flag,umap_flag=umap_flag)
   # return calculate_roc_auc(targets,obj.obs.score) , calculate_roc_auc(targets,obj.obs.SigScore) 
    if not ROC_flag:
        d,p =  calculate_aupr(targets,obj.obs.score) , calculate_aupr(targets,obj.obs.SigScore) 
    else:
        d,p = calculate_roc_auc(targets,obj.obs.score) , calculate_roc_auc(targets,obj.obs.SigScore)
    make_result_barplot(d,p,cluster_name,ROC_flag)
    return d,p

def import_seurat_to_scanpy(path=r"./data/pbmc_SELP.RData", temp_path = 'C:/Users/ronsh/Downloads'):
        r("library(Seurat)")
        r("library(SeuratDisk)")
        r("obj = readRDS(r'(" + path +")')")
        r("obj = UpdateSeuratObject(obj)")
      #  r("SaveH5Seurat(obj, 'data/obj.h5Seurat',overwrite = TRUE)")
        r("SaveH5Seurat(obj, '" + temp_path+  "/obj.h5Seurat' ,overwrite = TRUE)")
        r("Convert('"+  temp_path + "/obj.h5Seurat', dest = 'h5ad', overwrite = TRUE,verbose = TRUE)")
        sc_adata = sc.read_h5ad(temp_path + "/obj.h5ad")
        sc.pp.neighbors(sc_adata, n_neighbors=15, n_pcs=20)
        return sc_adata

def add_noise(obj,alpha = 0.80, drop_out = False):
    obj_noise = obj.copy()
    #obj_noise.X = (1-alpha) *obj_noise.X + alpha*np.random.randn(*obj.X.shape)
    if drop_out:
        obj_noise.X = obj_noise.X * np.random.binomial(1,(1-alpha),obj.X.shape)
    else:
        obj_noise.X = (1-alpha) *obj_noise.X + alpha*np.random.randn(*obj.X.shape)
    sc.tl.pca(obj_noise, svd_solver='arpack',use_highly_variable = False)
    sc.pp.neighbors(obj_noise,n_pcs=15, n_neighbors=100)

    return obj_noise

if __name__ == "__main__":
    obj = sc.read_h5ad("./data/scRNA-LN-compressed.h5ad")
    up_sig = pd.read_csv(r".\signatures\treg_up.csv",index_col=0) 
    down_sig =  pd.read_csv(r".\signatures\treg_down.csv",index_col=0) 
    d_score, p_score = test_scanpy_obj_roc(obj,up_sig.squeeze(),"Tregs",'cell_types',"Tregs",sig_down=down_sig.squeeze(),ROC_flag=False,umap_flag=True)
