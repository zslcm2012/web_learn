#!/usr/bin/python3
# -*- coding:utf-8 -*-

## 2019/02/19
## modify IB_PGS'TYPE to IB_PGS
## modify IB_NICS'TYPE to IB_NICS
## added Project_YKCS_PGS_*_LC TYPE to PGS_LC

## 2019/02/25
## modify _PGS_ and NO logo resolove 10M
## when the parents' karyotype both are empty ,modify its resolove to Whole_arm

## 2019/02/27
## debug when --flag is "Y", new_csv has only title

## 2019/03/01
## added PGDC resolove 4M

## 2019/03/07 version V1.4.2
## modify the module of parse karyotype to match when the chromosomes joined by (,|;)
## added parse karyotype to match when karyotype is descriped:dic(15;21)
## added sequence_system_paltform argument:(I/L),default = I
## when project_name re.search _NICS1_ ,the TYPE equal NICS1,the resolove equal Whole_Arm
## when project_name re.search _NICS2_ ,the TYPE equal NICS,the resolove equal 10M
## modify the type to PGS when project_name re.search Project_YKCS_PGS_*_LC_*

## 2019/03/08 version V1.4.3
## when the resolove is not conventional,print info to confirm
## when the karyotype_str re.search AZF,target chr is chrY
## when the karyotype_str re.search XXX,target chr is chrX
## when the karyotype_str re.search XYY,target chr is chrY
## when the karyotype_str as t(Y;3),target chr is sucessfully drawed

## 2019/04/16 version V1.4.4
## when the _CNV_ projects' resolove is not 10M, lib is not MALBAC_1_step_lib,the Project_warning log is made to result/local_report
## modify parse karyotype to chr_one,chr_another:inv(\d+|\w+)
## 2019/04/29 version V1.4.4
## test git

import os
import sys
import csv
import re
import argparse


### make dir ###
def safe_make_dir(path):
    try:
        os.makedirs(path)
    except OSError:
        if os.path.exists(path):
            pass
        else:
            raise

### read SJD to get the infomation ###
def parse_sjd(lims_dir,run_id):
    pro_sjd_info = {}
    run_flag = "_" + run_id
    sjd_head_dict = {"男方染色体核型": "ManKaryotype", "女方染色体核型": "WomanKaryotype", "报告是否带logo和章": "logo", "基因": "gene",
                     "分辨率": "resolove", "其他需要关注的染色体/核型": "others"}
    dict_key = ["ManKaryotype","WomanKaryotype","logo","gene","resolove","others","Merge","Parent_karyo","target_karyo"]
    dict_value = ["NA","NA","NA","NA","NA","NA","NA","NA","NA"]
    files = os.listdir(lims_dir)
    for file in files:
        file = str(file)
        if file.endswith(run_flag):
            project_sjd = "Project_" + file
            if project_sjd not in pro_sjd_info:
                pro_sjd_info[project_sjd] =dict(zip(dict_key,dict_value))

                sjdfile = file + ".txt"
                sjdfile_path = lims_dir + "/" + file + "/" + sjdfile
                sjd_info = open(sjdfile_path, 'r', encoding='gbk') # newline = ""
                head_sjd = sjd_info.readline().strip('\n')
                detail_sjd = sjd_info.readline().strip('\n')
                head_sjd_list = head_sjd.split('\t')
                detail_sjd_list = detail_sjd.split('\t')
                for num in range(0, len(head_sjd_list)):
                    for zip_key in pro_sjd_info[project_sjd].keys():
                        if (detail_sjd_list[num]):
                            if (sjd_head_dict.get(head_sjd_list[num]) == zip_key):
                                pro_sjd_info[project_sjd][zip_key] = detail_sjd_list[num]
                                if (sjd_head_dict[head_sjd_list[num]] == "resolove"):
                                    if (len(detail_sjd_list[num].split('+')) > 1):
                                        pro_sjd_info[project_sjd]['resolove'] = detail_sjd_list[-1].split('+')[0]
                                        print(project_sjd + ' 需要分析两种分辨率:' + detail_sjd_list[-1].split('+')[0] + '和' + detail_sjd_list[-1].split('+')[1])
                sjd_info.close()

    return pro_sjd_info


### write new csv for CNV analysis project ###
def write_new_csv(head_01,head_02,csv_cnv_list,new_csv_path):
    new_csv_info = csv.writer(open(new_csv_path, 'w'), delimiter=',', lineterminator='\n', quoting=csv.QUOTE_ALL)
    new_csv_info.writerow(head_01)
    new_csv_info.writerow(head_02)
    for csv_line in csv_cnv_list:
        new_csv_info.writerow(csv_line)

### from csv get two dict:  project_line  and   FILE_DICT['lib']  and get the new csv of  cnv_analysis  ###

def read_samplesheet(samplesheet_path,flag):
    samplesheet_info = samplesheet_path.split('/')
    lims_dir = "/".join(samplesheet_info[0:-1])
    samplesheet_name = samplesheet_info[-1]
    samplesheet_name_info = samplesheet_name.split('_')
    run_id = samplesheet_name_info[0]
    FILE_DICT = {}
    project_line = {}
    ONCPGD_dict = {}

    csv_cnv_list = []
    csv_ONCPGD_cnv_list = []
    list1 = ['PGS_', 'NICS_', '_NICS1_','_NICS2_','_CNV_', '_ONCPGD_', '_CPGD_']
    list2 = ['_PGD', '_MaReCs_']
    run_dir = os.getcwd()
    new_csv_path = run_dir + '/' + run_id + '_CNV_samplesheet.csv'
    ONCPGD_new_csv_path = run_dir + '/' + run_id + '_ONCPGD_CNV_samplesheet.csv'

    try:
        fh = open(samplesheet_path,'r')
    except IOError:
        print('Could not read file:' + samplesheet_path)
        sys.exit()

    reader = csv.reader(fh)
    head_01 = next(reader)
    head_02 = next(reader)
    for line in reader:
        lane_number = line[1]
        project = line[10]
        project_run = project + '_' + run_id
        sampleid = line[3]
        wga_lib = line[6]
        if project_run not in FILE_DICT:
            FILE_DICT[project_run] = {}
        if not 'LIB' in FILE_DICT[project_run]:
            FILE_DICT[project_run]['LIB'] = wga_lib
        if not 'CNV' in FILE_DICT[project_run]:
            FILE_DICT[project_run]['CNV'] = 'NA'
        if not 'SNP' in FILE_DICT[project_run]:  ### 为SNP分析预留的接口
            FILE_DICT[project_run]['SNP'] = 'NA'
        if not 'TYPE' in FILE_DICT[project_run]:
            FILE_DICT[project_run]['TYPE'] = 'NA'

        for type1 in list1:
            if re.search(type1,project_run):
                csv_cnv_list.append(line)
                if re.search('_CNV_',project_run):
                    FILE_DICT[project_run]['TYPE'] = "PGS"
                else:
                    FILE_DICT[project_run]['TYPE'] = type1.replace('_','')
                FILE_DICT[project_run]['CNV'] = 'CNV'

                if (flag == 'Y' or flag == 'y'):
                    if re.search('_ONCPGD_',project_run):
                        csv_ONCPGD_cnv_list.append(line)
                    if re.search('_PGS_',project_run) and re.search('_Control_',project_run):
                        csv_ONCPGD_cnv_list.append(line)

        for type2 in list2:
            if re.search(type2,project_run) and not re.search('SNP',sampleid):
                csv_cnv_list.append(line)
                FILE_DICT[project_run]['TYPE'] = type2.replace('_','')
                FILE_DICT[project_run]['CNV'] = 'CNV'
            if re.search(type2, project_run) and  re.search('SNP', sampleid):
                FILE_DICT[project_run]['SNP'] = 'SNP'

###------- the special project_type---------------------------###
        if re.search('_IB', project_run) and re.search('PGS', project_run):
            FILE_DICT[project_run]['TYPE'] = "IB_PGS"
        if re.search('_IB_NICS_', project_run) or re.search('_IBNICS_', project_run):
            FILE_DICT[project_run]['TYPE'] = "IB_NICS"
        if re.search('Project_YKCS_PGS_', project_run) and re.search('_LC_', project_run):
            FILE_DICT[project_run]['TYPE'] = "PGS_LC"
        if re.search('_NICS2_', project_run):
            FILE_DICT[project_run]['TYPE'] = "NICS"
###-----------------------------------------------------------###

    if (flag == 'Y' or flag == 'y'):
        write_new_csv(head_01,head_02,csv_ONCPGD_cnv_list,ONCPGD_new_csv_path)
    else:
        write_new_csv(head_01,head_02,csv_cnv_list,new_csv_path)

    return(lims_dir,run_id,project_line,FILE_DICT,ONCPGD_dict)

### parse karyotype to chr_one,chr_another ###
def parse_karyotype(karyotye_str):
    karyo_list = []
    chr_list = []
    if re.search('rob', karyotye_str) or re.search('t', karyotye_str) or re.search('der', karyotye_str) or re.search('dic',karyotye_str):
        pattern_list = re.findall("(?:['rob','der','t','dic']\((\w+|\d+)[';',','](\w+|\d+)\)?[';']?(\w+|\d+)?\)?)",karyotye_str)
        #pattern_list = re.findall("(?:['rob','der','t','dic']\((\d+)[';',','](\d+)\)?[';']?(\d+)?\)?)",karyotye_str)
        #pattern_list = re.findall("(?:['rob','der','t','dic']\((\d+)[';',','](\d+)\)?[';',',']?(\d+)?\)?)", karyotye_str)
        #pattern_list = re.findall("(?:['rob','der','t']\((\d+);(\d+)\)?;?(\d+)?\)?)", karyotye_str)
        # pattern_list = re.findall("(?:['rob','der','t']\((\d+);(\d+)\))", karyotye_str)
        if pattern_list:
            for num in range(0, len(pattern_list[0])):
                if pattern_list[0][num]:
                    karyo_list.append(pattern_list[0][num])
    if re.search('inv\(', karyotye_str):
        pattern_list_01 = re.findall("inv\((\d+|\w+)\)", karyotye_str)
        karyo_list.append(pattern_list_01[0])
    if re.search('AZF',karyotye_str):
        karyo_list.append('Y')
    if re.search('XYY',karyotye_str):
        karyo_list.append('Y')
    if re.search('XXX',karyotye_str):
        karyo_list.append('X')
    karyo_list = list(set(karyo_list))
    for chr_num in karyo_list:
        chr_num = 'chr' + chr_num
        chr_list.append(chr_num)

    chr_karyo = ','.join(chr_list)
    return chr_karyo


### parse_resolove when sjd's resolove is empty ###
def parse_resolove(FILE_CNV_DICT):
    list4 = ['_IB_PGS_','_IBPGS_','_CNV_','_IB_NICS_','_IBNICS_']                              #10M
    list5 = ['_PGD_','_ONCPGD_','_CPGD_','_PGDM_','_PGDA_','_PGDC_','_PGS_']                  #4M
    list6 = ['_MaReCs_']                                                                      #1M
    list_resolove = ['Whole_Chromosome','Whole_Arm','10M','4M','1M']
    for proname in FILE_CNV_DICT.keys():
        if FILE_CNV_DICT[proname]['resolove'] == "NA":
            for pro_type_04 in list4:
                if re.search(pro_type_04,proname):
                    FILE_CNV_DICT[proname]['resolove'] = '10M'
            for pro_type_05 in list5:
                if re.search(pro_type_05,proname):
                    FILE_CNV_DICT[proname]['resolove'] = '4M'
            for pro_type_06 in list6:
                if re.search(pro_type_06,proname):
                    FILE_CNV_DICT[proname]['resolove'] = '1M'
            if re.search('_NICS_', proname):
                if ( FILE_CNV_DICT[proname]['Parent_karyo'] != "NA,NA"):
                    if re.fullmatch('46,XY,46,XX',FILE_CNV_DICT[proname]['Parent_karyo']):
                        FILE_CNV_DICT[proname]['resolove'] = 'Whole_Arm'
                    else:
                        FILE_CNV_DICT[proname]['resolove'] = '10M'
                else:
                    FILE_CNV_DICT[proname]['resolove'] = 'Whole_Arm'

###-------------------NICS1 & NICS2 resolove------------------------------------------------------------###
            if re.search('_NICS1_', proname):
                FILE_CNV_DICT[proname]['resolove'] = 'Whole_Arm'
            if re.search('_NICS2_',proname):
                FILE_CNV_DICT[proname]['resolove'] = '10M'
        else:
            if (re.search('_NICS1_', proname) and FILE_CNV_DICT[proname]['resolove'] != 'Whole_Arm'):
                print(proname + '  -----please to confirm the resolove from inspection')
                FILE_CNV_DICT[proname]['resolove'] = 'Whole_Arm'
            if (re.search('_NICS2_', proname) and FILE_CNV_DICT[proname]['resolove'] != '10M'):
                print(proname + '  -----please to confirm the resolove from inspection')
                FILE_CNV_DICT[proname]['resolove'] = '10M'
###-----------------------------------------------------------------------------------------------------###
###----if the resolove from is not in['Whole_Chromosome','Whole_Arm','10M','4M','1M'],print info to confirm----###
    for pro_key in FILE_CNV_DICT.keys():
        if (FILE_CNV_DICT[pro_key]['resolove'] not in list_resolove):
            print(pro_key + ' *****the resolove is special,please confirm')
    return FILE_CNV_DICT

### according to FILE_CNV_DICT to write_project_CNV_list ###
def write_project_CNV_list(FILE_CNV_DICT,project_CNV_list_file):
    project_CNV_list_info = open(project_CNV_list_file,'w')
    for project_name in FILE_CNV_DICT.keys():
        project_CNV_list_info.writelines(project_name + '\t' + FILE_CNV_DICT[project_name]['Final_CNV_info'] + '\n')
    project_CNV_list_info.close()

### according to FILE_CNV_DICT to write_project_config_json ###
def write_project_config_json(FILE_CNV_DICT,project_config_json_dir,sequence_system_platform):
    safe_make_dir(project_config_json_dir)
    class_type = {'Whole_Chromosome':'1', 'Whole_Arm':'2', '10M':'3','4M':'4','1M':'5'}
    diruuid = '    "diruuid": "b83d20c7-656a-4edf-a1f9-606abb1d257f",\n'
    report_mosaic = '    "report_mosaic": true,\n'
    filetype = '    "filetype": ".fastq.gz",\n'
    report_gender = '    "report_gender": true,\n'
    clientid = '    "clientid": "NysECxCCTf4rXAxLAAAE",\n'
    if (sequence_system_platform == "I" or sequence_system_platform == "i"):
        platform = '    "platform": "Illumina",\n'
    if (sequence_system_platform == "L" or sequence_system_platform == "l"):
        platform = '    "platform": "Ion_Torrent",\n'
    for key_pro in FILE_CNV_DICT.keys():
        project_config_json_file = project_config_json_dir + '/' + key_pro + "_config.json"
        project_id = key_pro
        project_id_str = '    "proj_id": "' + project_id + '",\n'
        reagent = FILE_CNV_DICT[key_pro]['LIB']
        reagent_str = '{\n    "reagent": "'  + reagent + '",\n'
        report_level  = str(class_type.get(FILE_CNV_DICT[key_pro]['resolove']))
        report_level_str = '    "report_level": "' + report_level + '"\n   }'
        if (FILE_CNV_DICT[key_pro]['TYPE'] == 'NICS' or FILE_CNV_DICT[key_pro]['TYPE'] == 'NICS1'):
            type  = "NICS"
        elif (FILE_CNV_DICT[key_pro]['TYPE'] == 'MaReCs'):
            type = "MaReCs"
        elif (FILE_CNV_DICT[key_pro]['TYPE'] == 'PGS_LC'):
            type = "PGS_LC"
        else:
            type = "PGS"
        type_str = '    "type": "' + type + '",\n'
        project_config_json_info = open(project_config_json_file,'w')
        project_config_json_info.writelines(reagent_str)
        project_config_json_info.writelines(diruuid)
        project_config_json_info.writelines(report_mosaic)
        project_config_json_info.writelines(project_id_str)
        project_config_json_info.writelines(type_str)
        project_config_json_info.writelines(filetype)
        project_config_json_info.writelines(report_gender)
        project_config_json_info.writelines(clientid)
        project_config_json_info.writelines(platform)
        project_config_json_info.writelines(report_level_str)
        project_config_json_info.close()
###when the _CNV_ projects' resolove is not 10M, lib is not MALBAC_1_step_lib,the Project_warning log is made to result/local_report ###
def warning_unresonable_info(FILE_CNV_DICT):
    for pro in FILE_CNV_DICT.keys():
        local_report_warning_dir  = "/data/Project/" + pro + "/result/local_report/"
        safe_make_dir(local_report_warning_dir)
        local_report_warning_log = "/data/Project/" + pro + "/result/local_report/" + pro + "_warning.log"
        if re.search('_CNV_',pro):
            for all_keys in FILE_CNV_DICT[pro].keys():
                if ( FILE_CNV_DICT[pro]['resolove'] != "10M") or (FILE_CNV_DICT[pro]['LIB'] != "MALBAC_1_step_lib"):
                    resolove_info = open(local_report_warning_log,'w',encoding='utf-8')
                    if ( FILE_CNV_DICT[pro]['resolove'] != "10M" and FILE_CNV_DICT[pro]['LIB'] == "MALBAC_1_step_lib"):
                        resolove_info.writelines("请确认该项目分辨率是不是：" + FILE_CNV_DICT[pro]['resolove'] + "?\n")
                    elif (FILE_CNV_DICT[pro]['resolove'] == "10M" and FILE_CNV_DICT[pro]['LIB'] != "MALBAC_1_step_lib"):
                        resolove_info.writelines("请确认该项目的扩增建库方式是不是："  + FILE_CNV_DICT[pro]['LIB'] + "?\n")
                    elif (FILE_CNV_DICT[pro]['resolove'] != "10M" and FILE_CNV_DICT[pro]['LIB'] != "MALBAC_1_step_lib"):
                        resolove_info.writelines("请确认该项目分辨率是不是：" + FILE_CNV_DICT[pro]['resolove'] +  "?\n")
                        resolove_info.writelines("请确认该项目的扩增建库方式是不是：" + FILE_CNV_DICT[pro]['LIB'] + "?\n")
                    resolove_info.close()


### parse argument ###
def main():
    parser = argparse.ArgumentParser(description="make project_CNV_info.list for a run.")
    parser.add_argument('--samplesheet_path',action="store", required=True,help="The path of samplesheet of Lims.")
    parser.add_argument('--sequence_system_platform',action="store", required=False,default="I",choices=["I", "i", "L", "l","default=I"], help="[default=I];The platform type of sysenquence system :I/i:Iluminna;L/l:Life.")
    parser.add_argument('--flag', action="store", required=False, default="N",choices=["Y", "y", "N", "n","default=N"],help="[default=N];whether do ONCPGD analysis.")
    args = parser.parse_args()
    samplesheet_path = args.samplesheet_path
    sequence_system_platform = args.sequence_system_platform
    flag = args.flag

### start to  call function ###
    lims_dir,run_id,project_line,FILE_DICT,ONCPGD_dict = read_samplesheet(samplesheet_path,flag)
    pro_sjd_info = parse_sjd(lims_dir,run_id,)
    NEW_FILE_DICT = {}
    for pro in pro_sjd_info.keys():
        if not pro in NEW_FILE_DICT:
            NEW_FILE_DICT[pro] = {}
            key_merge = '{},{},{},{}'.format(pro_sjd_info[pro].get('ManKaryotype'), pro_sjd_info[pro].get('WomanKaryotype'),pro_sjd_info[pro].get('others'),pro_sjd_info[pro].get('gene'))
            key_parent_karyo = '{},{}'.format(pro_sjd_info[pro].get('ManKaryotype'), pro_sjd_info[pro].get('WomanKaryotype'))
            pro_sjd_info[pro]['Merge'] = str(key_merge)
            pro_sjd_info[pro]['Parent_karyo'] = str(key_parent_karyo)
            if (pro_sjd_info[pro]['Parent_karyo'] != "NA,NA"):
                karyotype_str = str(pro_sjd_info[pro]['Parent_karyo'])
                chr_karyo = parse_karyotype(karyotype_str)
                pro_sjd_info[pro]['target_karyo'] = chr_karyo
                if re.search('arr',karyotype_str):
                    print(pro + ' ' + karyotype_str)

        NEW_FILE_DICT[pro] = dict(FILE_DICT[pro], ** pro_sjd_info[pro])

### extract the project of only CNV analysis  ###

    FILE_CNV_DICT = {}

    for pro_name in NEW_FILE_DICT.keys():
        if ( NEW_FILE_DICT[pro_name]['CNV'] == 'CNV'):
            if not pro_name in FILE_CNV_DICT:
                FILE_CNV_DICT[pro_name] = {}
            FILE_CNV_DICT[pro_name].update(NEW_FILE_DICT[pro_name])


    FILE_CNV_DICT = parse_resolove(FILE_CNV_DICT)

    for project_name in FILE_CNV_DICT.keys():
        #last_merge = '{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(FILE_CNV_DICT[project_name].get('CNV'),FILE_CNV_DICT[project_name].get('TYPE'),FILE_CNV_DICT[project_name].get('resolove'),FILE_CNV_DICT[project_name].get('target_karyo'),FILE_CNV_DICT[project_name].get('Merge'),FILE_CNV_DICT[project_name].get('logo'),FILE_CNV_DICT[project_name].get('LIB'))
        last_merge = '{}\t{}\t{}\t{}\t{}\t{}'.format(FILE_CNV_DICT[project_name].get('CNV'),FILE_CNV_DICT[project_name].get('TYPE'),FILE_CNV_DICT[project_name].get('resolove'),FILE_CNV_DICT[project_name].get('target_karyo'),FILE_CNV_DICT[project_name].get('Merge'),FILE_CNV_DICT[project_name].get('logo'))
        last_merge = last_merge.replace('NA','')
        last_merge = last_merge.replace('None','')
        last_merge = '{}\t{}'.format(last_merge,FILE_CNV_DICT[project_name].get('LIB'))
        last_merge = last_merge.strip('\n')
        if not 'Final_CNV_info' in FILE_CNV_DICT[project_name]:
            FILE_CNV_DICT[project_name]['Final_CNV_info'] = last_merge

### judeg _CNV_ projects resolove and LIB whether is right ###
    #warning_unresonable_info(FILE_CNV_DICT)
### extract ONCPGD and PGS Control project to make a new dict:FILE_ONCPGD_CNV_DICT ###

    ONCPGD_CNV_DICT = {}
    for pro_on in FILE_CNV_DICT.keys():
        if re.findall("PGS(.*)Control", pro_on) or re.search('_ONCPGD_', pro_on):
            if not pro_on in ONCPGD_CNV_DICT:
                ONCPGD_CNV_DICT[pro_on] = {}

            ONCPGD_CNV_DICT[pro_on].update(FILE_CNV_DICT[pro_on])

### write_project_CNV_info.list ###

    pwd_dir = os.getcwd()
    project_CNV_list = pwd_dir + '/' + run_id + '_project_CNV_info.list'
    project_config_json_dir = pwd_dir + '/' + run_id + '_project_config'
    ONCPGD_CNV_list = pwd_dir + '/' + run_id + '_ONCPGD_project_CNV_info.list'
    ONCPGD_config_json_dir = pwd_dir + '/' + run_id + '_ONCPGD_project_config'


    if flag == 'N' or flag == 'n':
        write_project_CNV_list(FILE_CNV_DICT,project_CNV_list)
        write_project_config_json(FILE_CNV_DICT,project_config_json_dir,sequence_system_platform)
    if flag == 'Y' or flag == 'y':
        write_project_CNV_list(ONCPGD_CNV_DICT,ONCPGD_CNV_list)
        write_project_config_json(ONCPGD_CNV_DICT,ONCPGD_config_json_dir,sequence_system_platform)

###---------------if projects  not in FILE_CNV_DICT and not in SNP_dict,print them----------###
    for proname_01 in pro_sjd_info.keys():
        if proname_01 not in FILE_CNV_DICT.keys():
            print(proname_01 + " 这个项目需要一起分析吗?")
        if proname_01 in FILE_CNV_DICT.keys():
            if re.search('_RE_',proname_01) or re.search('_Re_',proname_01) or re.search('_re_',proname_01):
                print(proname_01,"  please confirm the resolove ")

if __name__ == '__main__':
    main()





