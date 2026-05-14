#-*- coding:gbk -*-
import re
from libs.getEvi import varRegimen
import copy
from libs.specialRequest import varInfo_FDZS
from libs.getConfig import clinicalNumStran, functionNumStran

'''
Discription
	
	处理MLPA格式。 
	返回格式：
	mlpa : {
		"B1_LOSS" : [],
		"B1_Gain" : [],
		"B2_LOSS" : [],
		"B2_Gain" : []
	}

'''

def process_mlpa(jsonDict, config):
	mlpa = copy.deepcopy(jsonDict["mlpa"])
	result = {}
	for var in mlpa:
		# evi_sum转换
		var["evi_sum"] = varRegimen(jsonDict, var["evi_sum"], config, var) if "evi_sum" in var.keys() and var["evi_sum"] else []
		# 复旦中山变异特殊需求
		var["var_info_forFDZS"] = varInfo_FDZS(var["gene_symbol"], "", "", config)

	# 返回格式
	for gene in ["BRCA1", "BRCA2"]:
		for i in ["Loss", "Gain"]:
			result["B"+gene[-1]+"_"+i] = [var for var in mlpa if var["gene_symbol"] == gene and re.search(i, var["type"])]

	# 返回MLPA图-2022.11.15
	mlpa_image = []
	for var in mlpa:
		if re.search("Loss|Gain", var["type"]):
			if var["file_path"] and var["file_path"] not in mlpa_image:
				mlpa_image.append(var["file_path"])

	# 返回MLPA del图-2023.03.30
	mlpa_image_del = []
	for var in mlpa:
		if re.search("Loss", var["type"]):
			if var["file_path"] and var["file_path"] not in mlpa_image_del:
				mlpa_image_del.append(var["file_path"])
	
	return result, mlpa_image,  mlpa_image_del
		
# 2024.10.31-增加一个规则-MLPA不再以del为疑似致病、dup为意义不明，看具体的clinical_significance来判断等级
def process_mlpa_v2(jsonDict, config):
	mlpa = copy.deepcopy(jsonDict["mlpa"])
	result = {}
	for var in mlpa:
		# evi_sum转换
		var["evi_sum"] = varRegimen(jsonDict, var["evi_sum"], config, var) if "evi_sum" in var.keys() and var["evi_sum"] else []
		# 复旦中山变异特殊需求
		var["var_info_forFDZS"] = varInfo_FDZS(var["gene_symbol"], "", "", config)
		var["clinic_num_g"] = clinicalNumStran(config).get(var["clinical_significance"], 3) if var["clinical_significance"] and var["clinical_significance"] != "-" else \
							  clinicalNumStran(config).get(var["function_classification"], 3)
		var["clinic_num_s"] = functionNumStran(config).get(var["function_classification"], 3) if var["function_classification"] and var["function_classification"] != "-" else \
							  functionNumStran(config).get(var["clinical_significance"], 3)
	# 返回格式
	'''
	{
	"B1_mlpa_1" : [], "B2_mlpa_1" : [],
	"B1_mlpa_2" : [], "B2_mlpa_2" : [],
	"B1_mlpa_3" : [], "B2_mlpa_3" : [],
	"B1_mlpa_4" : [], "B2_mlpa_4" : [],
	"B1_mlpa_5" : [], "B2_mlpa_5" : []
	}
	'''
	for gene in ["BRCA1", "BRCA2"]:
		for level in [1,2,3,4,5]:
			result["B"+gene[-1]+"_mlpa_L"+str(level)] = [var for var in mlpa if var["gene_symbol"] == gene and var["clinic_num_g"] == level]
	return result