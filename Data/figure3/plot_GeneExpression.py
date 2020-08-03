#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2020/7/30 20:02
# @Author  : YinLei Hu
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
def plot_TSNE_Heatmap(TSNE=None,x = None,
                      s = 1, vmax = None, vmin = None,
                      cmap="coolwarm",title = None,save_path = None):
    if save_path:
        pdf = PdfPages(save_path+title+".pdf")
        plt.figure()
        plt.tight_layout()
        p = plt.scatter(TSNE[:, 0], TSNE[:, 1],
                        s=s, c=x, vmax=vmax,vmin=vmin,
                        cmap=cmap)
        plt.colorbar(p)
        plt.title(title)
        print("savefig...")
        pdf.savefig()
        plt.close()
        pdf.close()
    else:
        p = plt.scatter(TSNE[:, 0], TSNE[:, 1],
                        s=s, c=x, vmax=vmax, vmin=vmin,
                        cmap=cmap)
        plt.colorbar(p)
        plt.title(title)
        plt.show()

path = "G:/program1/10x_data/code/figure3/"

Data_raw = pd.read_csv(path + "/Raw_fig3_Data.csv", header = 0, index_col=0)
Data_WEDGE = pd.read_csv(path + "WEDGE_fig3_Data.csv", header = 0, index_col=0)

TSNE_raw = pd.read_csv(path + "Raw_fig3_TSNE.csv", header = 0, index_col=0)
TSNE_WEDGE = pd.read_csv(path + "WEDGE_fig3_TSNE.csv", header = 0, index_col=0)

##########################################################################################
plot_TSNE_Heatmap(TSNE_WEDGE.values,x=Data_raw.values[0,:],title="Cr2"
                  ,vmin= 0.5,vmax=2.75,save_path = path+"WEDGE_TSNE_RawValue_")
plot_TSNE_Heatmap(TSNE_WEDGE.values,x=Data_WEDGE.values[0,:], title="Cr2",
                  vmin= 1.5,vmax=2, save_path = path+"WEDGE_TSNE_WEDGEValue_")

##########################################################################################
plot_TSNE_Heatmap(TSNE_WEDGE.values,x=Data_raw.values[1,:],title="Fcer2a"
                  ,vmin= 1,vmax=2.7,save_path = path+"WEDGE_TSNE_RawValue_")
plot_TSNE_Heatmap(TSNE_WEDGE.values,x=Data_WEDGE.values[1,:], title="Fcer2a",
                  vmin= 1.4,vmax=1.9, save_path = path+"WEDGE_TSNE_WEDGEValue_")

##########################################################################################
plot_TSNE_Heatmap(TSNE_raw.values,x=Data_raw.values[0,:],title="Cr2"
                  ,vmin= 0.5,vmax=2.75,save_path = path+"Raw_TSNE_RawValue_")
plot_TSNE_Heatmap(TSNE_raw.values,x=Data_raw.values[1,:],title="Fcer2a"
                  ,vmin= 1,vmax=2.7,save_path = path+"Raw_TSNE_RawValue_")

##########################################################################################
plot_TSNE_Heatmap(TSNE_WEDGE.values,x=Data_WEDGE.values[2,:],title="Cd93",
                  vmin= 0.0,vmax=2.0, save_path = path+"WEDGE_TSNE_WEDGEValue_")
plot_TSNE_Heatmap(TSNE_WEDGE.values,x=Data_raw.values[2,:],title="Cd93",
                  vmin= 0.0,vmax=2.0,save_path = path+"WEDGE_TSNE_RawValue_")

##########################################################################################
plot_TSNE_Heatmap(TSNE_WEDGE.values,x=Data_raw.values[3,:],title="nGene",
                  vmin= 600,vmax=1000,save_path = path+"WEDGE_TSNE_nGene_cut")


#DCA
##########################################################################################
Data = pd.read_csv(path + "DCA_Data.csv", header=0, index_col=0)
TSNE = pd.read_csv(path + "DCA_TSNE.csv", header=0, index_col=0)
plot_TSNE_Heatmap(TSNE.values,x=Data.values[0,:] ,title="Cr2",
                  vmin= 0.75,vmax=2,save_path = path+"DCA_TSNE_DCA_Value_")
plot_TSNE_Heatmap(TSNE.values,x=Data.values[1,:],title="Fcer2a",
                  vmin= 0.6,vmax=1.8,save_path = path+"DCA_TSNE_DCA_Value_")

#SAVERX
##########################################################################################
Data = pd.read_csv(path + "SAVERX_Data.csv", header=0, index_col=0)
TSNE = pd.read_csv(path + "SAVERX_TSNE.csv", header=0, index_col=0)
plot_TSNE_Heatmap(TSNE.values,x=Data.values[0,:],title="Cr2",
                  vmin= 0,vmax=2.5,save_path = path+"SAVERX_TSNE_SAVERX_Value_")
plot_TSNE_Heatmap(TSNE.values,x=Data.values[1,:],title="Fcer2a",
                  vmin= 0,vmax=2.5,save_path = path+"SAVERX_TSNE_SAVERX_Value_")

#MAGIC
##########################################################################################
Data = pd.read_csv(path + "MAGIC_Data.csv", header=0, index_col=0)
TSNE = pd.read_csv(path + "MAGIC_TSNE.csv", header=0, index_col=0)
plot_TSNE_Heatmap(TSNE.values,x=Data.values[0,:],title="Cr2",
                  vmin= 3,vmax=6,save_path = path+"MAGIC_TSNE_MAGIC_Value_")
plot_TSNE_Heatmap(TSNE.values,x=Data.values[1,:],title="Fcer2a",
                  vmin= 1.5,vmax=5.5,save_path = path+"MAGIC_TSNE_MAGIC_Value_")

#ENHANCE
##########################################################################################
Data = pd.read_csv(path + "ENHANCE_Data.csv", header=0, index_col=0)
TSNE = pd.read_csv(path + "ENHANCE_TSNE.csv", header=0, index_col=0)
plot_TSNE_Heatmap(TSNE.values,x=Data.values[0,:],title="Cr2",
                  vmin= 1.0,vmax=1.4,save_path = path+"ENHANCE_TSNE_ENHANCE_Value_")
plot_TSNE_Heatmap(TSNE.values,x=Data.values[1,:],title="Fcer2a",
                  vmin= 1.0,vmax=1.5,save_path = path+"ENHANCE_TSNE_ENHANCE_Value_")
#ALRA
##########################################################################################
Data = pd.read_csv(path + "ALRA_Data.csv", header=0, index_col=0)
TSNE = pd.read_csv(path + "ALRA_TSNE.csv", header=0, index_col=0)
plot_TSNE_Heatmap(TSNE.values,x=Data.values[0,:],title="Cr2",
                  vmin= 1.0,vmax=2.4,save_path = path+"ALRA_TSNE_ALRA_Value_")
plot_TSNE_Heatmap(TSNE.values,x=Data.values[1,:],title="Fcer2a",
                  vmin= 1.0,vmax=2.4,save_path = path+"ALRA_TSNE_ALRA_Value_")