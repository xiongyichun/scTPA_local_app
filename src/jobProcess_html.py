import sys,os


def pasteJobProcess(cur_num,message):
    if cur_num == 1:
        return '''<div id="output"><p>Your job is running.</p><p>(1/5) Start running.</p>'''
    elif cur_num == 2:
        return '''<div id="output"><p>Your job is running.</p><p>(1/5) Start running.</p><p>(2/5) Data is being pre-processed.</p>'''
    elif cur_num == 3:
        return '''<div id="output"><p>Your job is running.</p><p>(1/5) Start running.</p><p>(2/5) Data is being pre-processed.</p><p>(3/5)  Computing pathway activity score successfully.</p>'''
    elif cur_num == 4:
        return '''<div id="output"><p>Your job is running.</p><p>(1/5) Start running.</p><p>(2/5) Data is being pre-processed.</p><p>(3/5) Computing pathway activity score successfully.</p><p>(4/5) Clustering cells successfully.</p>'''
    elif cur_num == 5:
        return  '''<div id="output"><p>Your job is running.</p><p>(1/5) Start running.</p><p>(2/5) Data is being pre-processed.</p><p>(3/5) Computing pathway activity score successfully.</p><p>(4/5) Clustering cells successfully.</p><p>(5/5) Identifying pathway signatures successfully.</p><p>Job finished.</p></div>'''
    elif cur_num == 6:
        return '''<div id="output"><p> Error: The file format didn't match the required one.</p><p> Please read the <a target="_blank" href="../../help.html#input">Help</a> for more information and check your file before re-uploading.</p>'''
    elif cur_num == 7:
        return '''<div id="output"><p> Error: Some cell names in the expression matrix are not included in the first column of the cell type label file.</p><p> Please check files and reupload.</p>'''
    elif cur_num == 8:
        return '''<div id="output"><p> Errors in loading expression matrix. </p><p> Please check your file before re-uploading.</p>'''
    elif cur_num == 9:
        return '''<div id="output"><p> Errors in choosing the pathway database. </p><p> Please contact author.</p>'''
    elif cur_num == 10:
        return '''<div id="output"><p> Errors in loading the user-defined pathway. </p><p> Please check your file before re-uploading.</p>'''
    elif cur_num == 11:
        return '''<div id="output"><p> Errors in calculating the pathway activity score matrix. </p><p> Please check your file before re-uploading.</p>'''
    elif cur_num == 12:
        return '''<div id="output"><p> Errors in clustering cells. </p><p> Please check your file before re-uploading.</p>'''
    elif cur_num == 13:
        return '''<div id="output"><p> Errors in identifying pathway signatures (posibble reason: file is too large; gene names dosn't appear in any pathways). </p><p> Please check your file before re-uploading.</p>'''
    elif cur_num == 14:
        return '''<div id="output"><p> Errors unexcepted. </p><p> Please check your file before re-uploading.</p>'''
    elif cur_num == 15:
        return '''<div id="output"><p> Errors in the missing values imputation. </p><p> Please check your file before re-uploading.</p>'''
    elif cur_num == 16:
        return '''<div id="output"><p> Error: Cannot allocate memory. </p>'''
    elif cur_num == 17:
        return '''<div id="output"><p> The uploaded expression profile is too big.</p> <p>Please obtain the results using result retrieve function in the Job Retrieve page via the corresponding IDs. </p><p> This procedure may cost more than 2 hours. </p>'''
    elif cur_num == 18:
        return f'''<div id="output"><p> {message} </p>'''
    
def pasteProcessBar(cur_num):
    if cur_num == 1:
        return '''<div id="progressBar" class="progress-bar  progress-bar-striped progress-bar-animated" role="progressbar" style="width: 10%;" aria-valuenow="25" aria-valuemin="0" aria-valuemax="10">10%</div>'''
    elif cur_num == 2:
        return '''<div id="progressBar" class="progress-bar  progress-bar-striped progress-bar-animated" role="progressbar" style="width: 30%;" aria-valuenow="25" aria-valuemin="0" aria-valuemax="30">30%</div>'''
    elif cur_num == 3:
        return '''<div id="progressBar" class="progress-bar  progress-bar-striped progress-bar-animated" role="progressbar" style="width: 70%;" aria-valuenow="25" aria-valuemin="0" aria-valuemax="70">70%</div>'''
    elif cur_num == 4:
        return '''<div id="progressBar" class="progress-bar  progress-bar-striped progress-bar-animated" role="progressbar" style="width: 90%;" aria-valuenow="25" aria-valuemin="0" aria-valuemax="90">90%</div>'''
    elif cur_num == 5:
        return  '''<div id="progressBar" class="progress-bar  progress-bar-striped progress-bar-animated" role="progressbar" style="width: 100%;" aria-valuenow="25" aria-valuemin="0" aria-valuemax="100">100%</div>'''
    
def pasteHTML(f_path, cur_num, out_file, message=None):
    #print(cur_num)
    ori_bar = '''<div id="progressBar" class="progress-bar  progress-bar-striped progress-bar-animated" role="progressbar" style="width: 0%;" aria-valuenow="25" aria-valuemin="0" aria-valuemax="0">0%</div>''' 
    ret_lines = ''
    with open(f_path) as f:
        lines = f.readlines()
    for i in range(len(lines)):
        #if i == 89:
        #    ret_lines = ret_lines + '    '+jobID+'                        </td>'+'\n'
        #elif i == 98:
        #    ret_lines += email+'                        </td>'+'\n'
        line_i = lines[i].strip()
        #if cur_num == 5 and i == 5:
        #    print(line_i)
        #    ret_lines = ret_lines + '<meta http-equiv="Refresh" content="5">' + '\n'
        if line_i[0:17] =='''<div id="output">''':
            process_line = pasteJobProcess(cur_num,message)
            #print(process_line)
            ret_lines += process_line + '\n'
        elif line_i[0:21] == '''<div id="progressBar"''':
            #print('aaaaaaaaaaaaabbbbbbbbbbbbbbbbbbbbbbbbbbbbcccccccccccccccccccccccccccc')
            if cur_num < 6:
                process_line = pasteProcessBar(cur_num)
                ret_lines += process_line + '\n'
            else:
                ret_lines += '\n'
        #elif i == 122 and len(line_i) == 0:
        #    ret_lines += '''<div class="lds-gear" style="width: 100%;height:100%">'''
        #    ret_lines += '\n'    
        else:
            ret_lines += lines[i]
    if not cur_num == 5:
        with open(out_file, 'w') as fw:
            fw.write(ret_lines)
   
    with open(f_path, 'w') as fw:
        fw.write(ret_lines)
        
if __name__ == '__main__':
    #f_path = '/var/www/sctpa.bio-data.cn/Web1203/result_frame.html'
    #out_file = f_path
    #out_file = '/var/www/sctpa.bio-data.cn/Web1203/jobProcess_temp_res.html'
    #print(sys.argv)
    #print('bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb------------------------------------------')
    folder = sys.argv[1]
    f_path = os.path.join(folder, 'processing.html')
    print(f_path)
    #out_file = f_path.replace('_frame', '')
    out_file = os.path.join(folder, 'result.html')
    print(out_file)
    
    #email = sys.argv[2]
    #print(email)
    cur_num = int(sys.argv[2])
    try:
        message = sys.argv[3]
    except:
        message = None
    #print(cur_num)
    #out_file = '/var/www/sctpa.bio-data.cn/Web1203/jobProcess_temp_res'+str(cur_num)+'.html'
    pasteHTML(f_path, cur_num, out_file, message)
