import pandas as pd
import os,sys



def writeTable(file):
    #html = '''<table id="table1" style="font-size: smaller; margin-top: -44.8px;" class="table table-bordered table-hover">'''
    html = ''
    f_table = pd.read_csv(file,header=None)
    ncol = f_table.shape[1]
    nrow = min(50, f_table.shape[0])
    for i in range(nrow):
        if i == 0:
               
            html += '''<tr data-index="0" style="font-weight: bold">'''
        else:
            i_row = str(i)
            html += f'''<tr data-index="{i_row}">'''
        for j in range(ncol):
            j_col = str(f_table.iloc[i,j])
            html += f'''<td style="">{j_col}</td>'''
        html += '</tr>'
    #html += '</tbody></table><div class="fixed-table-border" style="height: 399px;"></div></div>'
    return html

def writePathways(data_folder):
    paths_html = ''
    allFiles = os.listdir(data_folder)
    i = 0
    for f in allFiles:
        if os.path.isdir(os.path.join(data_folder, f)):
            if i == 0:
                first = f
            paths_html += f'''<option value="{f}">{f}</option>'''
            i+=1
    return paths_html,first



def pasteHTML(f_path, jobID, email, out_file, data_folder, str_isCelltype): 
    ############### xiong
    #out_folder = 'sctpa.bio-data.cn/sctpa/results/'
    dirCurrent=os.path.split(os.path.realpath(__file__))[0]
    print("#######jonResult_html.py:")
    print(dirCurrent)
    out_folder=dirCurrent.split("/WEB-INF")[0]+"/results/"
    ############### xiong end
    ret_lines = ''
    with open(f_path) as f:
        lines = f.readlines()
    del_line = 0
    for i in range(len(lines)):
        f_i = lines[i].strip()
        if f_i[0:14] == '''<td id="jobId"''':
            ret_lines += f_i.replace('example1', jobID)+'\n'
        #    print('aaaaaaaaaaaaaaaaaaaaaaaaa')
        elif f_i[0:14] == '''<td id="email"''':
            ret_lines += f_i.replace('N/A',email)+'\n'
        #    print('bbbbbbbbbbbb')
        elif f_i[0:24] == '''<a id="downloadCellType"''':
            if str_isCelltype == 'yes':
                ret_lines += '\n'
            else:
                ret_lines += f_i.replace('example1', jobID)+'\n'
        elif f_i[0:19] == '''<a id="downloadPDF"''':
            ret_lines += f_i.replace('example1', jobID)+'\n'
        #    print('ddddddddd')
        elif f_i[0:26] == '''<a id="downloadHeatmapPDF"''':
            ret_lines += f_i.replace('example1', jobID)+'\n'
        #    print('ddddddddd')
        elif f_i[0:25] == '''<a id="downloadBubblePDF"''':
            ret_lines += f_i.replace('example1', jobID)+'\n'
        elif f_i[0:36] == '''<a id="downloadPathwayActivityScore"''':
            ret_lines += f_i.replace('example1', jobID)+'\n'
            print("pas")
        elif f_i[0:29] == '''<a id="downloadMarkerPathway"''':
            ret_lines += f_i.replace('example1', jobID)+'\n'
        #    print('ddddddddd')
        #elif f_i[0:28] == '''<div id="nonLinearReduction"''' and str_isCelltype == 'yes':
        #    ret_lines += '''<div id="nonLinearReduction" class="container my-4" style="display: none"'''
        #    ret_lines += '\n'
        #    print('ddddddd')
        else:
            ret_lines += lines[i]
    
    with open(out_file, 'w') as fw:
        fw.write(ret_lines)
    print(out_file)    
        
if __name__ == '__main__':
    
    ############### xiong
    #temp_html_path = '/var/www/sctpa.bio-data.cn/sctpa/index2221.html'
    #temp_html_path_noDiff = '/var/www/sctpa.bio-data.cn/sctpa/index_noMarkers.html'
    #temp_html_path_noDime = '/var/www/sctpa.bio-data.cn/sctpa/index_noDimen.html'
    #temp_html_path_noClust = '/var/www/sctpa.bio-data.cn/sctpa/index_noClust.html'

    dirCurrent=os.path.split(os.path.realpath(__file__))[0]
    dirResults=dirCurrent+"/../data/result_template/"
    temp_html_path = dirResults+'index2221.html'
    temp_html_path_noDiff = dirResults+'index_noMarkers.html'
    temp_html_path_noDime = dirResults+'index_noDimen.html'
    temp_html_path_noClust = dirResults+'index_noClust.html'
    ############### xiong end

    data_folder = sys.argv[1]
    jobID = sys.argv[2]
    email = sys.argv[3]
    str_isCelltype = sys.argv[4]
    has_diffGene = sys.argv[5]
    print(has_diffGene)    
    if has_diffGene == 'yes':
        f_path = temp_html_path
    elif has_diffGene == 'noDime':
        f_path = temp_html_path_noDime
    elif has_diffGene == 'noClust':
        f_path = temp_html_path_noClust
    else:
        f_path = temp_html_path_noDiff
    
    out_file = os.path.join(data_folder,'index.html')
    #print(out_file)
    if email == 'unacquainted':
        email = ''
    pasteHTML(f_path, jobID, email, out_file, data_folder, str_isCelltype)
