#!usr/bin/python3
# -*- coding:utf-8 -*-
# @Time     :2019/3/6 11:09
# @Author   : zhousai@yikongenomics.com
# @File     : md5sum_for_Miseq_GW000059.py
# @Function : GW000059 raw_data_dir include MD5.txt,after copy these raw_data to our service,need to make new MD5SUM

import os
import re
import sys
import subprocess
import argparse

def read_old_MD5_value(old_value_dir):
    line_number_MD5 = {}
    MD5_raw_name = {}
    MD5_new_name = {}
    MD5_old_file = old_value_dir + '/MD5.txt'
    line_num = 0
    with open(MD5_old_file,'r',encoding='gbk') as OLD:
        for line in OLD.readlines():
            line = line.strip()
            line_num += 1
            line_number_MD5[line_num] = line

    for number_MD5 in line_number_MD5.keys():
        if re.search('.fastq.gz',line_number_MD5[number_MD5]):
            #print(number_MD5)
            temp_list = line_number_MD5[number_MD5].split('\\')
            fastq_name = temp_list[-1]
            fastq_name = fastq_name.replace('):','')
            temp_list_Md5 = line_number_MD5[number_MD5 + 1].split(' ')
            Md5_raw_value = "".join(temp_list_Md5)
            MD5_raw_name[fastq_name] = Md5_raw_value
