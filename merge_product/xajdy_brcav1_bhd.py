#-*- coding:gbk -*-
import os
import json
from bin import getSampleInfo
from bin import getQC
from bin import getVar

#---------------------
# 本程序用来处理西安交大一BRCA扩增子+BHD出一份报告的需求
# 1. 输入sample_list，解析各个样本的json文件，输出结构化数据（gbrca_v1_data 和 bhd_data）
# 2. 根据模板展示内容对结构化数据进行合并
#  需要展示的内容如下：
#  1）订单信息-以bhd的内容为主，增加一个brcav1
#  2）数据质控-分为brcav1和bhd
#  3）胚系变异结果-从brcav1的结果中获取
#  4）体细胞变异-从bhd数据中提取，并且跟brcav1作比较，区分胚系和体细胞
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
def get_somatic(brca_v1_data, bhd_data):
    # 1. brca_v1汇总变异为列表供参考
    # 2. 遍历bhd中的变异，如果在brca_v1中存在的，var_origin对应改为germline，否则为somatic
    # 3. 额外处理：bhd体细胞I/II/III类
    # 不考虑brca_v1检出，bhd未检出的情况
    brca_v1_var_list = {}
    for i in ["B1_L5", "B2_L5", "B1_L4", "B2_L4", "B1_L3", "B2_L3", "B1_L2", "B2_L2", "B1_L1", "B2_L1"]:
        for var in brca_v1_data[i]:
            var_info = var["gene_symbol"]+var["hgvs_c"]+var["hgvs_p"]
            brca_v1_var_list[var_info] = var["evi_sum"]
    for i in ["level_I", "level_II", "level_onco_nodrug", "level_III"]:
        for var in bhd_data[i]:
            var_info = var["gene_symbol"]+var["hgvs_c"]+var["hgvs_p"] if var["bio_category"] == "Snvindel" else "HD_var"
            var["var_origin"] = "germline" if var_info in brca_v1_var_list.keys() else "somatic"
    somatic_I = [var for var in bhd_data["level_I"] if var["var_origin"] == "somatic"]
    somatic_II = [var for var in bhd_data["level_II"] if var["var_origin"] == "somatic"]
    somatic_III = [var for var in bhd_data["level_onco_nodrug"] + bhd_data["level_III"] if var["var_origin"] == "somatic"]
    return somatic_I, somatic_II, somatic_III

# 1. 遍历样本信息，逐个样本进行处理，获取get_data函数转化过的数据
def stran_json_data(sample_list, json_dir):
    brca_v1_data = {}
    bhd_data = {}
    for sample in sample_list:
        data = get_data(sample["order_code"], json_dir, config)
        if data["sample"]["prod_names"] == "BRCA1/BRCA2（全血）":
            brca_v1_data = data
        elif data["sample"]["prod_names"] == "BRCA1/BRCA2（组织 全血）":
            bhd_data = data
    return brca_v1_data, bhd_data

# 2. 合并转化数据
def merge_data(brca_v1_data, bhd_data):
    render_data = {
        "sample" : {},
        "var" : {"brca_v1" : {}, "bhd" : {}},
        "qc" : {"brca_v1" : {}, "bhd" : {}},
        "lib_qc" : {"brca_v1" : {}, "bhd" : {}}
    }
    render_data["sample"] = bhd_data["sample"]
    render_data["sample"]["brca_v1"] = brca_v1_data["sample"]
    render_data["qc"]["brca_v1"] = brca_v1_data["qc"]["dna_data_qc"]
    render_data["qc"]["bhd"] = bhd_data["qc"]["dna_data_qc"]
    render_data["lib_qc"]["brca_v1"]  = brca_v1_data["lib_quality_control"]["lib_dna_qc"]
    render_data["lib_qc"]["bhd"]  = bhd_data["lib_quality_control"]["lib_dna_qc"]
    render_data["var"]["brca_v1"] = brca_v1_data["var_brca"]["snv_s"]
    somatic_I, somatic_II, somatic_III = get_somatic(brca_v1_data["var_brca"]["snv_s"], bhd_data["var"]["var_somatic"])
    render_data["var"]["bhd"]["somatic_I"] = somatic_I
    render_data["var"]["bhd"]["somatic_II"] = somatic_II
    render_data["var"]["bhd"]["somatic_III"] = somatic_III
    return render_data

# 3. 主函数-返回合并后的数据
def MergeReport_main(sample_list, outfile):
    brca_v1_data, bhd_data = stran_json_data(sample_list, outfile)
    render_data = merge_data(brca_v1_data, bhd_data)
    return render_data