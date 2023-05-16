import pandas as pd
import numpy as np 
import scipy as sp
from scipy.stats import ranksums
import signature_utils as su
#mport pingouin as pg
from sklearn.neighbors import kneighbors_graph
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.metrics import roc_curve, auc
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt 
import networkx as nx 
import test_utils as tu
import scanpy as sc
from scipy.spatial import distance_matrix
import joblib
import seaborn as sns
import os

# testing pathway - not yet
# def test_learning_signature():
#     exp = pd.read_csv(r".\data\exp_tregs.csv",index_col=0)
#     idents_tregs = pd.read_csv(r".\targets\\cluster2_idents.csv",index_col=0) 
#     pca =  pd.read_csv(r".\graphs\pca_tregs.csv",index_col=0)
#     tu.build_signature(exp,pca, idents_tregs)
    # up = pd.read_csv(r"C:\Users\Ron\Desktop\MHCII\tregs_dataset_up.csv",index_col=0) 
    # down = pd.read_csv(r"C:\Users\Ron\Desktop\MHCII\tregs_down_up.csv",index_col=0) 
    # graph = pd.read_csv(r"C:\Users\Ron\Desktop\MHCII\graph_tregs.csv",index_col=0) 
    # sig_scores = signature_values(exp,up,down=down)
    # pred = propagition(sig_scores,graph)
    # calculate_roc_auc(idents_tregs, pred)


def test_reducing_features():
    exp = pd.read_csv(r".\data\exp_counts.csv",index_col=0)
    idents_naive = pd.read_csv(r".\targets\idents_naive.csv",index_col=0) 
    naive_up =  pd.read_csv(r".\signatures\naive_up.csv",index_col=0) 
    naive_down =  pd.read_csv(r".\signatures\naive_down.csv",index_col=0) 
    aucs_naive = tu.reduce_genes(exp,naive_up,idents_naive,naive_down)
    
    idents_tregs = pd.read_csv(r".\targets\idents_tregs.csv",index_col=0) 
    tregs_up =  pd.read_csv(r".\signatures\treg_up.csv",index_col=0) 
    treg_down =  pd.read_csv(r".\signatures\treg_down.csv",index_col=0) 
    aucs_tregs = tu.reduce_genes(exp,tregs_up,idents_tregs,treg_down)
    
    idents_prog = pd.read_csv(r".\data\idents_prog.csv",index_col=0) 
    prog_up =  pd.read_csv(r".\signatures\prog_up.csv",index_col=0) 
    prog_down =  pd.read_csv(r".\signatures\prog_down.csv",index_col=0) 
    aucs_prog = tu.reduce_genes(exp,prog_up,idents_prog,prog_down)
    return aucs_naive, aucs_tregs, aucs_prog


def test_reducing_sigs():
    # load data
    exp = pd.read_csv(r".\data\exp.csv",index_col=0)   
    graph =pd.read_csv(r".\graphs\graph.csv",index_col=0)
    
    # naive t cells
    idents_naive = pd.read_csv(r".\targets\idents_naive.csv",index_col=0) 
    naive_up =  pd.read_csv(r".\signatures\naive_up.csv",index_col=0) 
    naive_down =  pd.read_csv(r".\signatures\naive_down.csv",index_col=0) 
    aucs_naive = tu.reduce_signature(exp,naive_up,graph,idents_naive,down=naive_down)
    
    # tregs
    idents_tregs = pd.read_csv(r".\targets\idents_tregs.csv",index_col=0) 
    tregs_up =  pd.read_csv(r".\signatures\treg_up.csv",index_col=0) 
    treg_down =  pd.read_csv(r".\signatures\treg_down.csv",index_col=0) 
    aucs_tregs = tu.reduce_signature(exp,tregs_up,graph,idents_tregs,down=treg_down)
    
    # proliferationg
    idents_prog = pd.read_csv(r".\targets\idents_prog.csv",index_col=0) 
    prog_up =  pd.read_csv(r".\signatures\prog_up.csv",index_col=0) 
    prog_down =  pd.read_csv(r".\signatures\prog_down.csv",index_col=0) 
    aucs_prog = tu.reduce_signature(exp,prog_up,graph,idents_prog,down=prog_down)

    return aucs_naive, aucs_tregs, aucs_prog


def test_reducing_sigs():
    # load data
    exp = pd.read_csv(r".\data\exp.csv",index_col=0)
    graph =pd.read_csv(r".\graphs\graph.csv",index_col=0)

    # naive t cells
    idents_naive = pd.read_csv(r".\targets\idents_naive.csv",index_col=0)
    naive_up =  pd.read_csv(r".\signatures\naive_up.csv",index_col=0)
    naive_down =  pd.read_csv(r".\signatures\naive_down.csv",index_col=0)
    aucs_naive = tu.reduce_signature(exp,naive_up,graph,idents_naive,down=naive_down)

    # tregs
    idents_tregs = pd.read_csv(r".\targets\idents_tregs.csv",index_col=0)
    tregs_up =  pd.read_csv(r".\signatures\treg_up.csv",index_col=0)
    treg_down =  pd.read_csv(r".\signatures\treg_down.csv",index_col=0)
    aucs_tregs = tu.reduce_signature(exp,tregs_up,graph,idents_tregs,down=treg_down)

    # proliferationg
    idents_prog = pd.read_csv(r".\targets\idents_prog.csv",index_col=0)
    prog_up =  pd.read_csv(r".\signatures\prog_up.csv",index_col=0)
    prog_down =  pd.read_csv(r".\signatures\prog_down.csv",index_col=0)
    aucs_prog = tu.reduce_signature(exp,prog_up,graph,idents_prog,down=prog_down)

    return aucs_naive, aucs_tregs, aucs_prog

def test_scanpy_with_noise():
    tregs_up= pd.read_csv(r".\signatures\treg_up.csv",index_col=0)
    treg_down =  pd.read_csv(r".\signatures\treg_down.csv",index_col=0)
    naive_up =  pd.read_csv(r".\signatures\naive_up.csv",index_col=0)
    naive_down =  pd.read_csv(r".\signatures\naive_down.csv",index_col=0)
    prog_up =  pd.read_csv(r".\signatures\prog_up.csv",index_col=0)
    prog_down =  pd.read_csv(r".\signatures\prog_down.csv",index_col=0)

    obj = tu.import_seurat_to_scanpy()
    obj_noise  = tu.add_noise(obj)

    tu.test_scanpy_obj_roc(obj_noise,tregs_up.squeeze(),4,"seurat_clusters","Tregs Noise",False,treg_down.squeeze(),ROC_flag=False)
    tu.test_scanpy_obj_roc(obj_noise,naive_up.squeeze(),3,"seurat_clusters","Naive Noise",False,naive_down.squeeze(),ROC_flag=False)
    tu.test_scanpy_obj_roc(obj_noise,prog_up.squeeze(),8,"seurat_clusters","Prolifrating Noise",False,prog_down.squeeze(),ROC_flag=False)


def test_sig_python(path="./data/scRNA-LN-compressed.h5ad"):
    obj = sc.read_h5ad(path)
    sc.pp.neighbors(obj, n_neighbors=10, n_pcs=20)
    naive_up =  pd.read_csv(r".\signatures\naive_up.csv",index_col=0)
    naive_down =  pd.read_csv(r".\signatures\naive_down.csv",index_col=0)
    su.run_signature_propagation_first(obj, naive_up, naive_down)
   # su.run(obj, naive_up, naive_down)


    #naive_up =  pd.read_csv(r".\signatures\naive_up.csv",index_col=0)
    #naive_down =  pd.read_csv(r".\signatures\naive_down.csv",index_col=0)
    #su.run_signature_on_obj(obj, naive_up, naive_down)

    #b_cell_up = pd.read_csv(r".\signatures\b_cell_up.csv",index_col=0)
    #su.run_signature_on_obj(obj, b_cell_up, None)

    #DC1_up = pd.read_csv(r".\signatures\cdc1_up.csv",index_col=0)
    #su.run_signature_on_obj(obj, DC1_up, None)

    #monocytes_up = pd.read_csv(r".\signatures\Monocytes_up.csv",index_col=0)
    #su.run_signature_on_obj(obj, monocytes_up, None)

def test_gene_propagation_signature(path):
    obj = sc.read_h5ad(path)
    sc.pp.neighbors(obj, n_neighbors=15, n_pcs=20)
    #naive_up =  pd.read_csv(r".\signatures\naive_up.csv",index_col=0)
    #naive_down =  pd.read_csv(r".\signatures\naive_down.csv",index_col=0)
    # naive t cells
    idents_naive = pd.read_csv(r".\targets\idents_naive.csv",index_col=0)
    naive_up =  pd.read_csv(r".\signatures\naive_up.csv",index_col=0)
    naive_down =  pd.read_csv(r".\signatures\naive_down.csv",index_col=0)
    su.run_signature_on_obj(obj, naive_up, naive_down,umap_flag=False)
    print(f"AUC Score of orig propSig: {tu.calculate_roc_auc(idents_naive,obj.obs.SigScore)}")
    print(f"AUPR Score of orig propSig: {tu.calculate_aupr(idents_naive,obj.obs.SigScore)}")

    su.run_signature_propagation_first(obj, naive_up, naive_down,umap_flag=False)
    print(f"AUC Score of gene_propSig: {tu.calculate_roc_auc(idents_naive,obj.obs.SigScore)}")
    print(f"AUPR Score of orig gene_propSig: {tu.calculate_aupr(idents_naive,obj.obs.SigScore)}")

    # tregs
    idents_tregs = pd.read_csv(r".\targets\idents_tregs.csv",index_col=0)
    tregs_up =  pd.read_csv(r".\signatures\treg_up.csv",index_col=0)
    treg_down =  pd.read_csv(r".\signatures\treg_down.csv",index_col=0)
    su.run_signature_on_obj(obj, tregs_up, treg_down,umap_flag=False)
    print(f"AUC Score of orig propSig: {tu.calculate_roc_auc(idents_tregs,obj.obs.SigScore)}")
    print(f"AUPR Score of orig propSig: {tu.calculate_aupr(idents_tregs,obj.obs.SigScore)}")

    su.run_signature_propagation_first(obj, tregs_up, treg_down,umap_flag=False)
    print(f"AUC Score of gene_propSig: {tu.calculate_roc_auc(idents_tregs,obj.obs.SigScore)}")
    print(f"AUPR Score of orig gene_propSig: {tu.calculate_aupr(idents_tregs,obj.obs.SigScore)}")

    # proliferationg
    idents_prog = pd.read_csv(r".\targets\idents_prog.csv",index_col=0)
    prog_up =  pd.read_csv(r".\signatures\prog_up.csv",index_col=0)
    prog_down =  pd.read_csv(r".\signatures\prog_down.csv",index_col=0)

    su.run_signature_on_obj(obj, prog_up, prog_down,umap_flag=False)
    print(f"AUC Score of orig propSig: {tu.calculate_roc_auc(idents_prog,obj.obs.SigScore)}")
    print(f"AUPR Score of orig propSig: {tu.calculate_aupr(idents_prog,obj.obs.SigScore)}")

    su.run_signature_propagation_first(obj, prog_up, prog_down,umap_flag=False)
    print(f"AUC Score of gene_propSig: {tu.calculate_roc_auc(idents_prog,obj.obs.SigScore)}")
    print(f"AUPR Score of orig gene_propSig: {tu.calculate_aupr(idents_prog,obj.obs.SigScore)}")

    #su.run_signature_on_obj(obj, monocytes_up, None)

    
def test_propagation_genes(path="./data/scRNA-LN-compressed.h5ad"):
    obj = sc.read_h5ad(path)
    sc.pp.neighbors(obj, n_neighbors=10, n_pcs=20)
    graph = obj.obsp["connectivities"].toarray()
    exp = obj.raw.to_adata().to_df().T
    return su.propagate_all_genes(graph,exp)


def get_neighbor_graph(st_adata):
    spots_location = pd.DataFrame(st_adata.obsm["location"], index=st_adata._obs.index, columns=["x", "y"])
    distance = distance_matrix(spots_location, spots_location)
    distance_df = pd.DataFrame(distance, index=spots_location.index, columns=spots_location.index)
    np.fill_diagonal(distance_df.values, float('inf'))
    distance_nei = distance_df < 200
    neighbors_list_df = distance_nei.apply(lambda row: (distance_nei.index[row].tolist()))

    neighbor_graph = pd.DataFrame(np.zeros((len(spots_location.index), len(spots_location.index))),index=spots_location.index, columns=spots_location.index)
    for cell in spots_location.index:
        neighbor_graph[cell].loc[neighbors_list_df[cell]] = 1

    return(neighbor_graph)

def test_sig_STD():
    PARENT_DIR = "../../Users/Ella/PycharmProjects/Spatial/"
    st_adata = sc.read_h5ad(f"{PARENT_DIR}ST-LN-compressed.h5ad")
    tregs_up = pd.read_csv(r".\signatures\treg_up.csv",index_col=0)
    treg_down = pd.read_csv(r".\signatures\treg_down.csv",index_col=0)
    spots_location = pd.DataFrame(st_adata.obsm["location"], index=st_adata._obs.index, columns=["x", "y"])

    # neighbor_graph
    neighbor_graph = get_neighbor_graph(st_adata)
    #su.run_signature_on_STD(st_adata, neighbor_graph, tregs_up, treg_down)


    # proportion types
    st_annot = joblib.load("../../Users/Ella/Spatial/st_data.joblib")
    proportions_type = st_annot.obsm["proportions"]
    proportions_type_normed = proportions_type / abs(proportions_type).sum(axis=0)

    A = kneighbors_graph(st_adata.to_df(), n_neighbors = 5, mode='connectivity', include_self=True)
    knn_neighbor_graph = pd.DataFrame(A.toarray(),index=spots_location.index, columns=spots_location.index)
    #new_neighbor = (neighbor_graph[neighbor_graph.columns]==1) & (knn_neighbor_graph[knn_neighbor_graph.columns]==1)

    # both knn and phisical
    new_neighbor = (neighbor_graph[neighbor_graph.columns]==1) & (knn_neighbor_graph[knn_neighbor_graph.columns]==1)
    su.run_signature_on_STD(st_adata, new_neighbor.astype(int), "cDC1 cells signature - spatial nei", "location", tregs_up)

    # tregs
    su.run_signature_on_STD(st_adata, knn_neighbor_graph, "Tregs signature", "location", tregs_up, treg_down)
    plot_proportion_sig(proportions_type, st_adata, "Tregs")

    # b cells #
    b_cell_up = pd.read_csv(r".\signatures\b_cell_up.csv",index_col=0)
    su.run_signature_on_STD(st_adata, knn_neighbor_graph, "B cells signature", "location", b_cell_up)
    plot_proportion_sig(proportions_type, st_adata, "B cells")
    su.run_signature_on_STD(st_adata, neighbor_graph, "B cells signature - spatial nei", "location", b_cell_up)

    DC1_up = pd.read_csv(r".\signatures\cdc1_up.csv",index_col=0)
    su.run_signature_on_STD(st_adata, knn_neighbor_graph, "cDC1 signature", "location", DC1_up)
    plot_proportion_sig(proportions_type, st_adata, "cDC1s")
    su.run_signature_on_STD(st_adata, neighbor_graph, "cDC1 cells signature - spatial nei", "location", DC1_up)

    monocytes_up = pd.read_csv(r".\signatures\Monocytes_up.csv",index_col=0)
    su.run_signature_on_STD(st_adata, knn_neighbor_graph, "monocytes signature", "location", monocytes_up)
    plot_proportion_sig(proportions_type, st_adata, "Monocytes")
    su.run_signature_on_STD(st_adata, neighbor_graph, "Monocytes signature - spatial nei", "location", monocytes_up)



def test_sig_HUM_STD():
    PARENT_DIR = "../../Users/Ella/PycharmProjects/Spatial/data/"
    adata_vis = sc.datasets.visium_sge(sample_id="V1_Human_Lymph_Node")

    #################################### gene signature for Malaria #####################################
    Atypical_signature = pd.read_csv(r".\signatures\Malaria signature\Atypical_signature.csv")
    Atypical_signature_sig = Atypical_signature[Atypical_signature["adj.P.Val"]<0.05]
    Atypical_signature_up = Atypical_signature_sig[Atypical_signature_sig["logFC"]>0]["Genes"]
    # length is 220

    Naive_signature = pd.read_csv(r".\signatures\Malaria signature\Naive_signature.csv")
    Naive_signature_sig = Naive_signature[Naive_signature["adj.P.Val"]<0.02]
    Naive_signature_up = Naive_signature_sig[Naive_signature_sig["logFC"]>0]["Genes"]
    # length is 226

    Memory_signature = pd.read_csv(r".\signatures\Malaria signature\CMBC_signature.csv")
    Memory_signature_sig = Memory_signature[Memory_signature["adj.P.Val"]<0.05]
    Memory_signature_up = Memory_signature_sig[Memory_signature_sig["logFC"]>0]["Genes"]
    # length is 50

    AcMemory_signature = pd.read_csv(r".\signatures\Malaria signature\AcMB_signature.csv")
    ACMemory_signature_sig = AcMemory_signature[AcMemory_signature["adj.P.Val"]<0.01]
    AcMemory_signature_up = ACMemory_signature_sig[ACMemory_signature_sig["logFC"]>0]["Genes"]
    # length is 143

    A = kneighbors_graph(adata_vis.to_df(), n_neighbors = 5, mode='connectivity', include_self=True)
    knn_neighbor_graph = pd.DataFrame(A.toarray(),index=adata_vis._obs.index, columns=adata_vis._obs.index)

    b_cell_up = pd.read_csv(r".\signatures\b_cell_up.csv",index_col=0)
    b_cell_up = b_cell_up["x"].str.upper()
    #su.run_signature_on_STD(adata_vis, knn_neighbor_graph, "B cells signature", "spatial", b_cell_up)

    su.run_signature_on_STD(adata_vis, knn_neighbor_graph, "Naive B cells signature", "spatial", Naive_signature_up)
    su.run_signature_on_STD(adata_vis, knn_neighbor_graph, "Atypical B cells signature", "spatial", Atypical_signature_up)
    su.run_signature_on_STD(adata_vis, knn_neighbor_graph, "Memory B cells signature", "spatial", Memory_signature_up)
    su.run_signature_on_STD(adata_vis, knn_neighbor_graph, "AcMemory B cells signature", "spatial", AcMemory_signature_up)


def plot_proportion_sig(proportions_type_normed, st_adata, type_cell):
    d = {f"proportion {type_cell}":proportions_type_normed[type_cell], f"signature score {type_cell}": st_adata.obs["SigScore"]}
    data = pd.DataFrame(data = d)
    data.sort_values(by = [f"proportion {type_cell}"],ascending = True, inplace=True)
    plt.scatter(data[f"signature score {type_cell}"], data[f"proportion {type_cell}"], alpha=0.5)
    plt.title(f"proportion to sig score {type_cell}")
    plt.ylabel(f'proportion {type_cell}')
    plt.xlabel(f'sig score {type_cell}')
    plt.show()


def test_alphas(path, min_alpha, max_alpha, jump):
    obj = sc.read_h5ad(path)
   # sc.pp.neighbors(obj, n_neighbors=15, n_pcs=20)
    idents_naive = pd.read_csv(r".\targets\idents_naive.csv",index_col=0)
    naive_up =  pd.read_csv(r".\signatures\naive_up.csv",index_col=0)
    naive_down =  pd.read_csv(r".\signatures\naive_down.csv",index_col=0)

    for alpha in range(min_alpha,max_alpha,jump):
        alpha /= 10
        print(alpha)
        su.run_signature_on_obj(obj, naive_up, naive_down,umap_flag=False)
        print(f"AUC Score of orig propSig: {tu.calculate_roc_auc(idents_naive,obj.obs.SigScore)}")
        print(f"AUPR Score of orig propSig: {tu.calculate_aupr(idents_naive,obj.obs.SigScore)}")


if __name__ == "__main__":
    os.chdir(r"C:\Users\ronsh\OneDrive\שולחן העבודה\Signature")
    #test_scanpy_with_noise()
    #matrix = test_propagation_genes(r"C:\Users\ronsh\OneDrive\שולחן העבודה\Signature\data\obj.h5ad")
    #print(matrix.shape)
    #test_gene_propagation_signature(r"C:\Users\ronsh\OneDrive\שולחן העבודה\Signature\data\obj.h5ad")
    test_alphas(r"C:\Users\ronsh\OneDrive\שולחן העבודה\Signature\data\obj_sig_test.h5ad", 5, 10, 1)
    # test_reducing_sigs()
    # test_sig_python()
    #test_sig_STD()
    #test_sig_HUM_STD()
