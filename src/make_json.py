import pandas as pd
import numpy as np
import json
import sys,os
#from clustergrammer2 import net

#from clustergrammer_widget import *
#net = Network(clustergrammer_widget)

from clustergrammer import Network
net = Network()

import seaborn as sns
import warnings
import multiprocessing as mp
warnings.filterwarnings('ignore')



class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        else:
            return super(NpEncoder, self).default(obj)
        
def save_to_json(inst_dict, filename, indent=True):
        fw = open(filename, 'w')
        if indent == True:
                fw.write( json.dumps(inst_dict, indent=2,cls=NpEncoder) )
        else:
                fw.write( json.dumps(inst_dict,cls=NpEncoder) )

def get_heatmap_json_mat(dataframe,celltype):
        cols = []
        for i in range(0,celltype.shape[0]):
                try:
                        float(celltype.iloc[0,0])
                        cell_type = 'cellType: '+str(celltype.iloc[i,0])
                except:
                         cell_type = 'cellType: '+ celltype.iloc[i,0]
                type_tuple = (cell_type, cell_type)
                #type_tuple = (celltype.index[i], cell_type)
                cols.append(type_tuple)
        dataframe.columns = cols
        return dataframe

def RGB_to_Hex(rgb):
    color = '#'
    for i in rgb:
        color+=str(hex(int(i*255)))[-2:].replace('x','0').upper()
    return color

def set_cat_colors(colors, cellTypes,axis, cat_index, cat_title=False):
    i = 0
    for inst_ct in range(len(cellTypes)):
        inst_ct = cellTypes[i]
        if cat_title != False:
            cat_name = cat_title + ': ' + inst_ct
        else:
            cat_name = inst_ct           
        inst_color = colors[i]
        net.set_cat_color(axis=axis, cat_index=cat_index, cat_name=cat_name, inst_color=inst_color)
        i+=1

def set_cat_colors2(colors, axis, cat_index, cat_title=False):
    i = 0
    for inst_ct in range(colors.shape[0]):
        inst_ct = colors.index[i]     ##cell name
        if cat_title != False:
            cat_name = cat_title + ': ' + str(inst_ct)
        else:
            cat_name = inst_ct
        inst_color = colors.iloc[i,0]
        net.set_cat_color(axis=axis, cat_index=cat_index, cat_name=cat_name, inst_color=inst_color)
        i+=1


def clustergrammer_json(mat_path,
                        celltype_path,
                        out_path,
                        colors):
    #print("############ clustergrammer_json")
    #print(mat_path)
    #print(celltype_path)
    mat = pd.read_csv(mat_path, index_col=0)
    cellType = pd.read_csv(celltype_path, index_col=0, header=None)
    cellType = cellType.loc[mat.columns.values.tolist()]
    #print(mat.shape)
    #print(cellType.shape)
    pas = get_heatmap_json_mat(mat, cellType)

    net.load_df(pas)

    set_cat_colors2(colors,'col',1,'cellType')

    
    net.cluster()
    #net.widget()
    

    net.write_json_to_file(net_type = 'viz',filename = out_path)


def scatter_2d(dataframe, title, x_lab, y_lab, colors):
    body = {'chart': {
        'type': 'scatter',
        'zoomType': 'xy',
         "marginRight": 100
    },
    'accessibility': {
        'description': 'A scatter plot compares the dimensional reduction of cells by cell type. '
    },
    'title': {
        'text': title
    },
    'xAxis': {
        'title': {
            'enabled': True,
            'text': x_lab
        },
        'startOnTick': True,
        'endOnTick': True,
        'showLastLabel': True
    },
    'yAxis': {
        'title': {
            'text': y_lab
        }
    },
    "legend": {
    "layout": "vertical",
    "align": "right",
    "verticalAlign": "middle",
    "floating": True,
    "backgroundColor": "white",
    "borderWidth": 1
    },
    'plotOptions': {
        'scatter': {
            'marker': {
                #"symbol": "circle",
                'radius': 4,
                'states': {
                    'hover': {
                        'enabled': True,
                        'lineColor': 'rgb(100,100,100)'
                    }
                }
            },
            'states': {
                'hover': {
                    'marker': {
                        'enabled': False
                    }
                }
            },
            'tooltip': {
                'headerFormat': '<b>{series.name}</b><br>',
                'pointFormat': '{point.x}, {point.y}'
            }
        }
    },
     'series': data_series(dataframe,colors)     
          }
    return body

def scatter_3d(dataframe, title, x_lab, y_lab,z_lab,colors):
    body = {
    'chart': {
        "marginRight": 100,
        'renderTo': 'container',
        'margin': 100,
        'type': 'scatter3d',
        'animation': False,
        'options3d': {
            'enabled': True,
            'alpha': 10,
            'beta': 30,
            'depth': 250,
            'viewDistance': 5,
            'fitToPlot': False,
            'frame': {
                'bottom': { 'size': 1, 'color': 'rgba(0,0,0,0.02)' },
                'back': { 'size': 1, 'color': 'rgba(0,0,0,0.04)' },
                'side': { 'size': 1, 'color': 'rgba(0,0,0,0.06)' }
            }
        }
    },
    'title': {
        'text': title
    },
    'xAxis': {
        'title': {
            'enabled': True,
            'text': x_lab
        },
    },
    'yAxis': {
        'title': {
            'text': y_lab
        }
    },
    'zAxis': {
        'title': {
            'enabled': True,
            'text': z_lab
        },
    },
    'plotOptions': {
        'scatter': {
            'width': 10,
            'height': 10,
            'depth': 10,
            'marker':{
                'symbol':'circle',
            }
        },
    },
    "legend": {
        "layout": "vertical",
        "align": "right",
        "verticalAlign": "middle",
        "floating": True,
        "backgroundColor": "white",
        "borderWidth": 1
    },
    'series': data_series(dataframe,colors)
}
    return body

def getColor(color_turple, hue=0.8):
    r = int(color_turple[0]*255)
    g = int(color_turple[1]*255)
    b = int(color_turple[2]*255)
    return f'rgba({r}, {g}, {b}, {hue})'
def getColor2(hex, hue = 0.8):
    r = int(hex[1:3],16)
    g = int(hex[3:5],16)
    b = int(hex[5:7], 16)
    return f'rgba({r}, {g}, {b}, {hue})'


def getData(dataframe):
    data = []
    for i in range(dataframe.shape[0]):
        data_i = []
        for j in range(dataframe.shape[1]-1):
            data_i.append(dataframe.iloc[i,j])
        data.append(data_i)
    return data

def data_series(dataframe, colors):
    #### 设定最后一列是细胞类型
    cell_types = list(set(dataframe.iloc[:,-1]))
    cell_types.sort()
    print("from make_json............................")
    print(cell_types)
    #colors = sns.hls_palette(len(cell_types), l=.65, s=1)
    series = []
    for i in range(len(cell_types)):
        #print(cell_types[i])
        series.append({'name':cell_types[i],
                      # 'color':'getColor(colors[i])',
                      'color':getColor2(colors.loc[cell_types[i]][0]),
                      'data':getData(dataframe[dataframe.iloc[:,-1] == cell_types[i]])
                      })
    return series  

def scatter_json(input_path,
                 out_path,
                 title, 
                 x_lab,
                 y_lab,
                 z_lab,
                 colors):
    dataframe = pd.read_csv(input_path,index_col=0)
    if dataframe.shape[1] == 3:
        ret_json = scatter_2d(dataframe, title, x_lab, y_lab, colors)
    elif dataframe.shape[1] == 4:
        ret_json = scatter_3d(dataframe, title, x_lab, y_lab, z_lab, colors)
    else:
        print(f'error for dim of {input_path}')
        return 0
    save_to_json(ret_json, filename = out_path)
    
def marker_json(input_path, output_path, files):
    ret = {}
    dataframe = pd.read_csv(input_path,dtype='str')
    cell_types = list(set(dataframe['cluster']))
    for i in range(len(cell_types)):
        cell = cell_types[i]
        ret[cell] = []
        cell_frame = dataframe[dataframe['cluster'] == cell]
        for j in range(cell_frame.shape[0]):
            if cell_frame.iloc[j,0] in files:
                ret[cell].append(cell_frame.iloc[j,0])
    save_to_json(ret, filename = output_path)

def paserFileMat(path_f,f, celltype_path, colors):
    #print(f)
    #if os.path.isdir(path_f):
    #print(path_f)
    input_path = os.path.join(path_f,f+'.csv')
    output_path = os.path.join(path_f,f+'.json')
    #print(output_path)
    if os.path.exists(input_path):
        #print(input_path)
        #print(output_path)
        
        clustergrammer_json(input_path, celltype_path, output_path, colors)
    else:
        print("no such file:",input_path)

def parseFolderMat(folder, celltype_path, colors):
    files = os.listdir(folder)
    #print(celltype_path)
    #print(colors)
    if 'uncompress' in files:
        files.remove('uncompress')
    #print(files)
    pool = mp.Pool()
    #print(files)
    #i = 0
    #print(files[1:5])
    for f in files:
        path_f = os.path.join(folder,f)
        if os.path.isdir(path_f):
            #print(f)
            #paserFileMat(path_f,f, celltype_path, colors)
            pool.apply_async(paserFileMat, (path_f,f, celltype_path, colors, ))
    #       if os.path.isdir(path_f):
    #        input_path = os.path.join(path_f,f+'.csv')
    #        output_path = os.path.join(path_f,f+'.json')
            #print(output_path)
            #clustergrammer_json(input_path, celltype_path, output_path, colors)
    #       if i == 0:
    #           output_path = os.path.join(folder,'first.json')
    #           clustergrammer_json(input_path, celltype_path, output_path)
    #           i = 1
    pool.close()
    pool.join()

def excel_data_json(dataframe):
    colnames = dataframe.columns
    ret = []
    for i in range(dataframe.shape[0]):
        i_dic = {}
        i_dic['id'] = i+1
        for j in range(len(colnames)):
            i_dic[colnames[j]] = dataframe.iloc[i,j]
        ret.append(i_dic)
    return ret

def excel_json(input_path, out_path):
    dataframe = pd.read_csv(input_path,dtype='str')
    dataframe = dataframe.drop(['geneList'],axis=1)
    n_rows = dataframe.shape[0]
    series_obj = {
        "total": n_rows,
        "totalNotFiltered": n_rows,
        "rows": excel_data_json(dataframe)}
    save_to_json(series_obj, filename = out_path)

if __name__ == '__main__':
    out_dir = sys.argv[1]
    
    ###### pas_heatmap_json
    #pas_path = os.path.join(out_dir, 'pas.csv')
    colors = pd.read_csv(os.path.join(out_dir, 'pas_color.csv'), index_col=0)
     
    celltype_path = os.path.join(out_dir, 'cell_type.csv')
    #pas_json_path = os.path.join(out_dir, 'pas.json')
    #clustergrammer_json(pas_path, celltype_path, pas_json_path)
    
    ###### tsne scatter
    tsne2d_path = os.path.join(out_dir, 'tsne_2D.csv')
    tsne3d_path = os.path.join(out_dir, 'tsne_3D.csv')
    tsne2d_json_path = os.path.join(out_dir, 'T-SNE-2D.json')
    tsne3d_json_path = os.path.join(out_dir, 'T-SNE-3D.json')
    scatter_json(tsne2d_path, tsne2d_json_path, 't-SNE 2D', 'tSNE_1', 'tSNE_2', 'none', colors)
    scatter_json(tsne3d_path, tsne3d_json_path, 't-SNE 3D', 'tSNE_1', 'tSNE_2', 'tSNE_3', colors)
    
    ##### umap scatter
    umap2d_path = os.path.join(out_dir, 'umap_2D.csv')
    umap3d_path = os.path.join(out_dir, 'umap_3D.csv')
    umap2d_json_path = os.path.join(out_dir, 'UMAP-2D.json')
    umap3d_json_path = os.path.join(out_dir, 'UMAP-3D.json')
    scatter_json(umap2d_path, umap2d_json_path, 'UMAP 2D', 'UMAP_1', 'UMAP_2', 'none', colors)
    scatter_json(umap3d_path, umap3d_json_path, 'UMAP 3D', 'UMAP_1', 'UMAP_2', 'UMAP_3', colors)
 
    ##### cell type markers
    files = os.listdir(out_dir)
    marker_path = os.path.join(out_dir, 'pas_markers.csv')
    marker_json_path = os.path.join(out_dir, 'cellType_pathways.json')
    marker_json(marker_path, marker_json_path, files)
    print("cell type markers success")
    
    ##### pas_json
    marker_json_path = os.path.join(out_dir, 'pas_markers_excel.json')
    excel_json(marker_path, marker_json_path)

    #### single pathway clustergrammer
    parseFolderMat(out_dir, celltype_path, colors)
    print("all folder success")

    #### excel json
    #tsne2d_json_path = os.path.join(out_dir, 'T-SNE-2D_excel.json')
    #excel_json(tsne2d_path, tsne2d_json_path)
    #tsne3d_json_path = os.path.join(out_dir, 'T-SNE-3D_excel.json')
    #excel_json(tsne3d_path, tsne3d_json_path)
    #umap2d_json_path = os.path.join(out_dir, 'UMAP-2D_excel.json')
    #excel_json(umap2d_path, umap2d_json_path)
    #umap3d_json_path = os.path.join(out_dir, 'UMAP-3D_excel.json')
    #excel_json(umap3d_path, umap3d_json_path)
    


