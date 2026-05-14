#-*- coding:gbk -*-
import os
import json
from bin import getSampleInfo
from bin import getQC
from bin import getVar
import copy

#---------------------
# 本程序用来处理山东齐鲁HRD+gBRCA出一份报告的需求
# 1. 输入sample_list，解析各个样本的json文件，输出结构化数据（gbrca_data 和 hrd_data）
# 2. 根据模板展示内容对结构化数据进行合并
#  需要展示的内容如下：
#  1）订单信息-仅展示hrd的内容，有需要再增加gbrca
#  2）数据质控-分为gbrca和hrd
#  3）胚系变异结果-从gbrca的结果中获取
#  4）体细胞变异-从hrd数据中提取，并且跟gbrca作比较，区分胚系和体细胞
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

# 2. 判断是否匹配到治疗/辅助诊断/预后的证据
def judgeRegimen(var):
	if var["evi_sum"]["evi_split"] and set(["Diagnostic","Predictive","Prognostic"]) & set(var["evi_sum"]["evi_split"].keys()):
		return 1
	else:
		return 0
      
# 致病/疑似致病但无用药的变异，在结果汇总表中归为III类
def S_level(var):
	s_level = var["clinic_num_s"] if judgeRegimen(var) else 3
	return s_level

def getSum_forGermlinePanel(var_list_raw, gene_list_all, gene_list_appr):
	var_list = copy.deepcopy(var_list_raw)
	result = []
	detect_gene = [var["gene_symbol"] for var in var_list]
	for var in var_list:
		# 致病/疑似致病但无用药的归为III类
		var["clinic_num_s"] = S_level(var)
		result.append(var)
	for gene in set(gene_list_all) - set(detect_gene):
		result.append({"gene_symbol" : gene})
	# 标记获批基因
	for i in result:
		if i["gene_symbol"] in  gene_list_appr:
			i["appr"] = "T"

	return sorted(result, key = lambda i:i["gene_symbol"])

def parp_summary(var_list):
    gene_list_all = ["BRCA1", "BRCA2", "ATM", "BARD1", "BRIP1", "CDH1", "CDK12", "CHEK1", "CHEK2", "FANCA", "FANCL", "HDAC2", "PALB2", "PPP2R2A", \
					"PTEN", "RAD51B", "RAD51C", "RAD51D", "RAD54L", "TP53"]
    gene_list_appr = ["BRCA1", "BRCA2", "ATM", "BARD1", "BRIP1", "CDK12", "CHEK1", "CHEK2", "FANCL", "PALB2", "RAD51B", "RAD51C", "RAD51D", "RAD54L", "FANCA"]
    return getSum_forGermlinePanel(var_list, gene_list_all, gene_list_appr)

# 过滤出体细胞变异
def get_somatic(gbrca_data, hrd_data):
    # 1. gbrca汇总变异为列表供参考
    # 2. 遍历hrd中的BRCA变异，如果在gbrca中存在的，var_origin对应改为germline，否则为somatic
    # 3. germline对应的丰度需要用gbrca的
    # 不考虑gbrca检出，hrd未检出的情况
    gbrca_var_list = {}
    for i in ["B1_L5", "B2_L5", "B1_L4", "B2_L4", "B1_L3", "B2_L3", "B1_L2", "B2_L2", "B1_L1", "B2_L1"]:
        for var in gbrca_data[i]:
            var_info = var["gene_symbol"]+var["hgvs_c"]+var["hgvs_p"]
            gbrca_var_list[var_info] = var["freq"]
    for i in ["level_I", "level_II", "level_onco_nodrug", "level_III", "level_IV"]:
        for var in hrd_data[i]:
            var_info = var["gene_symbol"]+var["hgvs_c"]+var["hgvs_p"]
            var["var_origin"] = "germline" if var_info in gbrca_var_list.keys() else "somatic" if var["gene_symbol"] in ["BRCA1", "BRCA2"] else ""
            var["freq"] = gbrca_var_list[var_info] if var_info in gbrca_var_list.keys() else var["freq"]

    brca1_var = [var for var in hrd_data["level_I"] + hrd_data["level_II"] if var["gene_symbol"] == "BRCA1"]
    brca2_var = [var for var in hrd_data["level_I"] + hrd_data["level_II"] if var["gene_symbol"] == "BRCA2"]

    somatic_I = hrd_data["level_I"]
    somatic_II = hrd_data["level_II"]
    somatic_onco_nodrug = hrd_data["level_onco_nodrug"]
    somatic_III = hrd_data["level_III"]

    parp_result = parp_summary(somatic_I + somatic_II + somatic_onco_nodrug + somatic_III)
    
    result = {
         "brca1_var" : brca1_var,
         "brca2_var" : brca2_var,
         "somatic_I" : somatic_I,
         "somatic_II" : somatic_II,
         "somatic_onco_nodrug" : somatic_onco_nodrug,
         "somatic_III" : somatic_III,
         "parp_result" : parp_result
    }
        
    return result

# 1. 遍历样本信息，逐个样本进行处理，获取get_data函数转化过的数据
def stran_json_data(sample_list, json_dir):
    gbrca_data = {}
    hrd_data = {}
    for sample in sample_list:
        data = get_data(sample["order_code"], json_dir, config)
        if data["sample"]["prod_names"] == "BRCA1/BRCA2（全血）":
            gbrca_data = data
        elif data["sample"]["prod_names"] == "HRD Complete（组织）":
            hrd_data = data
    return gbrca_data, hrd_data

# 2. 合并转化数据
def merge_data(gbrca_data, hrd_data):
    render_data = {
        "sample" : {},
        "var" : {"gbrca" : {}, "hrd" : {}},
        "qc" : {"gbrca" : {}, "hrd" : {}},
        "lib_qc" : {"gbrca" : {}, "hrd" : {}}
    }
    render_data["sample"] = hrd_data["sample"]
    #render_data["sample"]["hrd"] = hrd_data["sample"]
    render_data["qc"]["gbrca"] = gbrca_data["qc"]["dna_data_qc"]
    render_data["qc"]["hrd"] = hrd_data["qc"]["dna_data_qc"]
    render_data["lib_qc"]["gbrca"]  = gbrca_data["lib_quality_control"]["lib_dna_qc"]
    render_data["lib_qc"]["hrd"]  = hrd_data["lib_quality_control"]["lib_dna_qc"]
    render_data["var"]["gbrca"] = gbrca_data["var_brca"]["snv_s"]
    render_data["var"]["hrd"] = get_somatic(gbrca_data["var_brca"]["snv_s"], hrd_data["var"]["var_somatic"])
    render_data["var"]["hrd"]["gss"] = hrd_data["var"]["gss"]["gss"]
    # PARP抑制剂相关标志物检测结果

    return render_data

# 3. 主函数-返回合并后的数据
def MergeReport_main(sample_list, outfile):
    gbrca_data, hrd_data = stran_json_data(sample_list, outfile)
    render_data = merge_data(gbrca_data, hrd_data)
    #dataJson = json.dumps(render_data, ensure_ascii = False)
    #with open(outfile+"/test_to_word.json", "w", encoding = "utf-8") as outFile:
     #    outFile.write(dataJson)
    return render_data