#-*- coding:gbk -*-
import os
import json
from bin import getSampleInfo
from bin import getQC
from bin import getVar

#---------------------
# 本程序用来处理中山肿瘤BRCA扩增子+BRCA组织出一份报告的需求
# 1. 输入sample_list，解析各个样本的json文件，输出结构化数据（gbrca_v1_data 和 tbrca_data）
# 2. 根据模板展示内容对结构化数据进行合并
#  需要展示的内容如下：
#  1）订单信息-以brcav1的内容为主，增加一个tbrca
#  2）数据质控-分为brcav1和tbrca
#  3）胚系变异结果-从brcav1的结果中获取
#  4）体细胞变异-从tbrca数据中提取，并且跟brcav1作比较，区分胚系和体细胞
#----------------------

config = os.path.join(os.path.abspath(os.path.dirname(os.path.dirname(__file__))), "config")

# 结构化数据生成
def get_data(json_name, outfile, config):
	data = {}
	# 解析json文件
	with open(outfile+"/"+json_name+".json", "r", encoding='utf-8') as file_json:
		jsonDict = json.load(file_json)
		jsonDict["sample_info"]["report_module_type"] = "rummage" if jsonDict["sample_info"]["report_module_type"] == "clinical" else jsonDict["sample_info"]["report_module_type"]
		jsonDict["sample_info"]["tumor_names_cn"] = jsonDict["sample_info"]["tumor_type"]
	data["sample"] = getSampleInfo.getSample(jsonDict)
	data["qc"], data["lib_quality_control"] = getQC.getJsonQC(jsonDict)
	data["var"], data["var_brca"], data["var_hrr_shsy"] = getVar.getVar(jsonDict, config, "")
	return data

# 过滤出体细胞变异
def get_somatic(brca_v1_data, tbrca_data):
    # 1. brca_v1汇总变异为列表供参考
    # 2. 遍历tBRCA中的变异，如果在brca_v1中存在的，var_origin对应改为germline，否则为somatic
    # 3. 额外处理：tbrca中的胚系5/4类、3类；brca_v1中的体细胞I/II类
    # 不考虑brca_v1检出，tbrca未检出的情况
    brca_v1_var_list = {}
    for i in ["B1_L5", "B2_L5", "B1_L4", "B2_L4", "B1_L3", "B2_L3", "B1_L2", "B2_L2", "B1_L1", "B2_L1"]:
        for var in brca_v1_data[i]:
            var_info = var["gene_symbol"]+var["hgvs_c"]+var["hgvs_p"]
            #print ("1111111111", var_info)
            brca_v1_var_list[var_info] = var["evi_sum"]
    for i in ["B1_L5", "B2_L5", "B1_L4", "B2_L4", "B1_L3", "B2_L3", "B1_L2", "B2_L2", "B1_L1", "B2_L1"]:
        for var in tbrca_data[i]:
            var_info = var["gene_symbol"]+var["hgvs_c"]+var["hgvs_p"]
            var["var_origin"] = "germline" if var_info in brca_v1_var_list.keys() else "somatic"
            var["evi_sum"] = brca_v1_var_list[var_info] if var_info in brca_v1_var_list.keys() else var["evi_sum"]
    somatic_345_brca1 = [var for var in tbrca_data["B1_L5"] + tbrca_data["B1_L4"] + tbrca_data["B1_L3"] if var["var_origin"] == "somatic"]
    somatic_345_brca2 = [var for var in tbrca_data["B2_L5"] + tbrca_data["B2_L4"] + tbrca_data["B2_L3"] if var["var_origin"] == "somatic"]
    somatic_12 = [var for var in tbrca_data["B1_L2"] + tbrca_data["B2_L2"] + tbrca_data["B1_L1"] + tbrca_data["B2_L1"] if var["var_origin"] == "somatic"] 
    #germline_45 = [var for var in tbrca_data["B1_L5"] + tbrca_data["B2_L5"] + tbrca_data["B1_L4"] + tbrca_data["B2_L4"] if var["var_origin"] == "germline"]
    #germline_3 = [var for var in tbrca_data["B1_L3"] + tbrca_data["B2_L3"] if var["var_origin"] == "germline"]
    somatic_I_II = [var for var in tbrca_data["B1_L5"] + tbrca_data["B2_L5"] + tbrca_data["B1_L4"] + tbrca_data["B2_L4"] if var["var_origin"] == "somatic"]
    # 2026.01.07-体系增加III类变异
    somatic_III = [var for var in tbrca_data["B1_L3"] + tbrca_data["B2_L3"] if var["var_origin"] == "somatic"]
    # 2026.01.07-增加完成  
    return somatic_345_brca1, somatic_345_brca2, somatic_12, somatic_I_II, somatic_III

# 1. 遍历样本信息，逐个样本进行处理，获取get_data函数转化过的数据
def stran_json_data(sample_list, json_dir):
    brca_v1_data = {}
    tbrca_data = {}
    for sample in sample_list:
        data = get_data(sample["order_code"], json_dir, config)
        if data["sample"]["prod_names"] == "BRCA1/2（扩增子）":
            brca_v1_data = data
        elif data["sample"]["prod_names"] == "BRCA1/BRCA2（组织）":
            tbrca_data = data
    return brca_v1_data, tbrca_data

# 2. 合并转化数据
def merge_data(brca_v1_data, tbrca_data):
    render_data = {
        "sample" : {},
        "var" : {"brca_v1" : {}, "tbrca" : {}},
        "qc" : {"brca_v1" : {}, "tbrca" : {}},
        "lib_qc" : {"brca_v1" : {}, "tbrca" : {}}
    }
    render_data["sample"] = brca_v1_data["sample"]
    render_data["sample"]["tbrca"] = tbrca_data["sample"]
    render_data["qc"]["brca_v1"] = brca_v1_data["qc"]["dna_data_qc"]
    render_data["qc"]["tbrca"] = tbrca_data["qc"]["dna_data_qc"]
    render_data["lib_qc"]["brca_v1"]  = brca_v1_data["lib_quality_control"]["lib_dna_qc"]
    render_data["lib_qc"]["tbrca"]  = tbrca_data["lib_quality_control"]["lib_dna_qc"]
    render_data["var"]["brca_v1"] = brca_v1_data["var_brca"]["snv_s"]
    somatic_345_brca1, somatic_345_brca2, somatic_12, somatic_I_II, somatic_III = get_somatic(brca_v1_data["var_brca"]["snv_s"], tbrca_data["var_brca"]["snv_s"])
    render_data["var"]["tbrca"]["somatic_345_brca1"] = somatic_345_brca1
    render_data["var"]["tbrca"]["somatic_345_brca2"] = somatic_345_brca2
    render_data["var"]["tbrca"]["somatic_12"] = somatic_12
    render_data["var"]["tbrca"]["somatic_I_II"] = somatic_I_II
    render_data["var"]["tbrca"]["somatic_III"] = somatic_III
    return render_data

# 3. 主函数-返回合并后的数据
def MergeReport_main(sample_list, outfile):
    brca_v1_data, tbrca_data = stran_json_data(sample_list, outfile)
    render_data = merge_data(brca_v1_data, tbrca_data)
    return render_data