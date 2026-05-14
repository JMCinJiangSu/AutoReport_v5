#-*- coding:gbk -*-
import os
import json
from bin import getSampleInfo
from bin import getQC
from bin import getVar

#---------------------
# 本程序用来处理中山肿瘤BRCA扩增子+HRD出一份报告的需求
# 1. 输入sample_list，解析各个样本的json文件，输出结构化数据（gbrca_v1_data 和 hrd_data）
# 2. 根据模板展示内容对结构化数据进行合并
#  需要展示的内容如下：
#  1）订单信息-以brcav1的内容为主，增加一个hrd
#  2）数据质控-分为brcav1和hrd
#  3）胚系变异结果-从brcav1的结果中获取
#  4）体细胞变异-从hrd数据中提取，并且跟brcav1作比较，区分胚系和体细胞
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

def parp_summary(var_list):
	gene_list = ["ATM","BARD1","BRCA1","BRCA2","BRIP1","CDK12","CHEK1","CHEK2","FANCL","PALB2","RAD51B","RAD51C","RAD51D","RAD54L", "FANCA"]
	detect_gene = [var["gene_symbol"] for var in var_list if var["gene_symbol"] in gene_list]
	result = []
	for var in var_list:
		result.append(var)
	for gene in set(gene_list) - set(detect_gene):
		result.append({
			"gene_symbol" : gene
		})	
	return sorted(result, key=lambda i:i["gene_symbol"])

# 过滤出体细胞变异
def get_somatic(brca_v1_data, hrd_data):
    # 1. brca_v1汇总变异为列表供参考
    # 2. 遍历tBRCA中的变异，如果在brca_v1中存在的，var_origin对应改为germline，否则为somatic
    # 3. 额外处理：tbrca中的胚系5/4类、3类；brca_v1中的体细胞I/II类
    # 不考虑brca_v1检出，tbrca未检出的情况
    brca_v1_var_list = {}
    for i in ["B1_L5", "B2_L5", "B1_L4", "B2_L4", "B1_L3", "B2_L3", "B1_L2", "B2_L2", "B1_L1", "B2_L1"]:
        for var in brca_v1_data[i]:
            var_info = var["gene_symbol"]+var["hgvs_c"]+var["hgvs_p"]
            brca_v1_var_list[var_info] = var["evi_sum"]
    for i in ["level_I", "level_II", "level_onco_nodrug", "level_III", "level_IV"]:
        for var in hrd_data[i]:
            var_info = var["gene_symbol"]+var["hgvs_c"]+var["hgvs_p"]
            var["var_origin"] = "germline" if var_info in brca_v1_var_list.keys() else "somatic"
            var["evi_sum"] = brca_v1_var_list[var_info] if var_info in brca_v1_var_list.keys() else var["evi_sum"]
    
    somatic_345_brca1_notinbrcav1 = [var for var in hrd_data["level_I"] + \
                                                    hrd_data["level_II"] + \
                                                    hrd_data["level_onco_nodrug"] + \
                                                    hrd_data["level_III"] if var["var_origin"] == "somatic" and var["gene_symbol"] == "BRCA1"]
    
    somatic_345_brca2_notinbrcav1 = [var for var in hrd_data["level_I"] + \
                                                    hrd_data["level_II"] + \
                                                    hrd_data["level_onco_nodrug"] + \
                                                    hrd_data["level_III"] if var["var_origin"] == "somatic" and var["gene_symbol"] == "BRCA2"]
    
    somatic_I_II_other = [var for var in hrd_data["level_I"] + \
                                       hrd_data["level_II"] if var["var_origin"] == "somatic" and var["gene_symbol"] not in  ["BRCA1", "BRCA2"]]
    
    somatic_III_other = [var for var in hrd_data["level_onco_nodrug"] + \
                                      hrd_data["level_III"] if var["var_origin"] == "somatic" and var["gene_symbol"] not in  ["BRCA1", "BRCA2"]]
    somatic_I_II_brca = [var for var in hrd_data["level_I"] + \
                                       hrd_data["level_II"] if var["var_origin"] == "somatic" and var["gene_symbol"] in  ["BRCA1", "BRCA2"]]
    somatic_I_II_brca1 = [var for var in hrd_data["level_I"] + \
                                       hrd_data["level_II"] if var["var_origin"] == "somatic" and var["gene_symbol"] == "BRCA1"]
    somatic_I_II_brca2 = [var for var in hrd_data["level_I"] + \
                                       hrd_data["level_II"] if var["var_origin"] == "somatic" and var["gene_symbol"] == "BRCA2"]
    parp_result = parp_summary(brca_v1_data["B1_L5"] + brca_v1_data["B2_L5"] + brca_v1_data["B1_L4"] + brca_v1_data["B2_L4"] + \
                               brca_v1_data["B1_L3"] + brca_v1_data["B2_L3"] + somatic_345_brca1_notinbrcav1 + somatic_345_brca2_notinbrcav1 + \
                               somatic_I_II_other + somatic_III_other)
    # 2026.03.06-新增体细胞BRCA III类变异
    somatic_III_brca = [var for var in hrd_data["level_onco_nodrug"] + \
                                       hrd_data["level_III"] if var["var_origin"] == "somatic" and var["gene_symbol"] in  ["BRCA1", "BRCA2"]]
    # 2026.03.06-新增结束
    
    result = {
         "somatic_345_brca1_notinbrcav1" : somatic_345_brca1_notinbrcav1,
         "somatic_345_brca2_notinbrcav1" : somatic_345_brca2_notinbrcav1,
         "somatic_I_II_other" : somatic_I_II_other,
         "somatic_III_other" : somatic_III_other,
         "somatic_I_II_brca" : somatic_I_II_brca,
         "somatic_I_II_brca1" : somatic_I_II_brca1,
         "somatic_I_II_brca2" : somatic_I_II_brca2,
         "parp_result" : parp_result,
         "somatic_III_brca" : somatic_III_brca
    }
        
    return result

# 1. 遍历样本信息，逐个样本进行处理，获取get_data函数转化过的数据
def stran_json_data(sample_list, json_dir):
    brca_v1_data = {}
    tbrca_data = {}
    for sample in sample_list:
        data = get_data(sample["order_code"], json_dir, config)
        if data["sample"]["prod_names"] == "BRCA1/2（扩增子）":
            brca_v1_data = data
        elif data["sample"]["prod_names"] == "HRD Complete（组织）":
            tbrca_data = data
    return brca_v1_data, tbrca_data

# 2. 合并转化数据
def merge_data(brca_v1_data, hrd_data):
    render_data = {
        "sample" : {},
        "var" : {"brca_v1" : {}, "hrd" : {}},
        "qc" : {"brca_v1" : {}, "hrd" : {}},
        "lib_qc" : {"brca_v1" : {}, "hrd" : {}}
    }
    render_data["sample"] = brca_v1_data["sample"]
    render_data["sample"]["hrd"] = hrd_data["sample"]
    render_data["qc"]["brca_v1"] = brca_v1_data["qc"]["dna_data_qc"]
    render_data["qc"]["hrd"] = hrd_data["qc"]["dna_data_qc"]
    render_data["lib_qc"]["brca_v1"]  = brca_v1_data["lib_quality_control"]["lib_dna_qc"]
    render_data["lib_qc"]["hrd"]  = hrd_data["lib_quality_control"]["lib_dna_qc"]
    render_data["var"]["brca_v1"] = brca_v1_data["var_brca"]["snv_s"]
    render_data["var"]["hrd"] = get_somatic(brca_v1_data["var_brca"]["snv_s"], hrd_data["var"]["var_somatic"])
    render_data["var"]["hrd"]["gss"] = hrd_data["var"]["gss"]["gss"]
    # PARP抑制剂相关标志物检测结果

    return render_data

# 3. 主函数-返回合并后的数据
def MergeReport_main(sample_list, outfile):
    brca_v1_data, hrd_data = stran_json_data(sample_list, outfile)
    render_data = merge_data(brca_v1_data, hrd_data)
    dataJson = json.dumps(render_data, ensure_ascii = False)
    with open(outfile+"/test_to_word.json", "w", encoding = "utf-8") as outFile:
         outFile.write(dataJson)
    return render_data