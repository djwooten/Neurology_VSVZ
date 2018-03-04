import pandas as pd
import numpy as np
from sklearn.manifold import TSNE, Isomap
from sklearn.decomposition import PCA
from matplotlib import pyplot as plt

def clean_sample_name(name):
    return name[:12].replace(".","-")

def plot(data,clusters,labels,title,xlab,ylab,fname=None, xmax=None):

    colors = ["#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6"]
    
    cluster_dict = dict()
    for x,y,c in zip(data[:,0],data[:,1],clusters):
        if cluster_dict.has_key(c):
            cluster_dict[c][0].append(x)
            cluster_dict[c][1].append(y)
        else: cluster_dict[c] = ([x,],[y,])
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i,c in enumerate(cluster_dict.keys()):
        ax.scatter(cluster_dict[c][0], cluster_dict[c][1], c=colors[i], s=40, label=labels[c])
    ax.legend(scatterpoints=1,numpoints=1)
    ax.grid()
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    plt.title(title)
    
    if xmax is not None: ax.set_xlim(ax.get_xlim()[0], xmax)
    
    if fname is None: plt.show()
    else: plt.savefig(fname)
    plt.cla()
    plt.close()
    plt.clf()

# Read in raw data
data = pd.read_csv("../expression.csv",index_col=0)
data.columns = [clean_sample_name(i) for i in data.columns]
for i in data.columns:
    data[i] = pd.to_numeric(data[i],errors='coerce')    
data = data.dropna().transpose()

# Read in SVZ status
svz_status = pd.read_csv("../data/patient_status.csv",index_col=0)

# Find which samples have known SVZ status
shared_samples = [i for i in data.index if i in svz_status.index]
shared_sample_indices = []
for i in data.index:
    shared_sample_indices.append(i in svz_status.index)
shared_sample_indices = np.where(shared_sample_indices)[0]

# Read SVZ status and subtype
svz = svz_status.loc[shared_samples,"SVZ"]
subtype = svz_status.loc[shared_samples,"Subtype"]
svz_labels = {0:"VSVZ-",1:"VSVZ+"}
subtype_labels = dict()
for i in subtype.unique(): subtype_labels[i] = i

# Do dimensionality reduction
y_pca = PCA(n_components=30).fit_transform(data)
y_tsne = TSNE(n_components=2).fit_transform(y_pca)
y_isomap = Isomap(n_components=2).fit_transform(y_pca)


plot(y_pca[shared_sample_indices,:],svz,svz_labels,"PCA: VSVZ Status","PC1","PC2",fname="svz_pca.pdf")
plot(y_tsne[shared_sample_indices,:],svz,svz_labels,"t-SNE: VSVZ Status","tsne1","tsne2",xmax=60,fname="svz_tsne.pdf")
plot(y_isomap[shared_sample_indices,:],svz,svz_labels,"Isomap: VSVZ Status","Isomap1","Isomap2",xmax=400,fname="svz_isomap.pdf")

plot(y_pca[shared_sample_indices,:],subtype,subtype_labels,"PCA: Subtype","PC1","PC2",fname="subtype_pca.pdf")
plot(y_tsne[shared_sample_indices,:],subtype,subtype_labels,"t-SNE: Subtype","tsne1","tsne2",xmax=60,fname="subtype_tsne.pdf")
plot(y_isomap[shared_sample_indices,:],subtype,subtype_labels,"Isomap: Subtype","Isomap1","Isomap2",xmax=400,fname="subtype_isomap.pdf")
