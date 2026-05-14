#-*- coding:gbk -*-
import re
from libs.getEvi import varRegimen
import copy
from libs.specialRequest import varInfo_FDZS
from libs.getConfig import clinicalNumStran, functionNumStran
from libs.rule import S_function

'''
Discription
	
	处理cnv格式。 (胚系变异)
	适用BRCA cnv Loss/Gain，用来代替mlpa-2024.08.30新增

'''

def process_gcnv(jsonDict, config):
	cnv = [var for var in copy.deepcopy(jsonDict["cnv"]) if var["var_origin"] == "germline" and var["gene_symbol"] in ["BRCA1", "BRCA2"]]
	result = {}
	for var in cnv:
		var["evi_sum"] = varRegimen(jsonDict, var["evi_sum"], config, var)
		var["var_info_forFDZS"] = varInfo_FDZS(var["gene_symbol"], "", "", config)
		# 字段名跟mlpa保持一致，这样模板少改一些
		var["type"] = var["cnv_type"]
		var["value"] = var["region_exon"]
		# 新增字段-type_raw-用于区分del中的Loss、Hetedel、homodel
		var["type_raw"] = var["type"]
		if var["type"] in ["HeteDel", "homoDel"]:
			var["type"] == "Loss"

	for gene in ["BRCA1", "BRCA2"]:
		for i in ["Loss", "Gain"]:
			result["B"+gene[-1]+"_"+i] = [var for var in cnv if var["gene_symbol"] == gene and re.search(i, var["cnv_type"])]
	return result

# 2024.10.31-增加一个规则-不再以del为疑似致病、dup为意义不明，看具体的clinical_significance来判断等级
def process_gcnv_v2(jsonDict, config):
	cnv = [var for var in copy.deepcopy(jsonDict["cnv"]) if var["var_origin"] == "germline" and var["gene_symbol"] in ["BRCA1", "BRCA2"]]
	result = {}
	for var in cnv:
		var["evi_sum"] = varRegimen(jsonDict, var["evi_sum"], config, var)
		var["var_info_forFDZS"] = varInfo_FDZS(var["gene_symbol"], "", "", config)
		# 字段名跟mlpa保持一致，这样模板少改一些
		#var["type"] = var["cnv_type"]
		# cnv_type用于区分del中的Loss、Hetedel、homodel、dup(Gain)，type保持跟mlpa一致
		var["type"] = "Gain" if var["cnv_type"] == "Gain" else "Loss"
		var["value"] = var["region_exon"]
		var["clinic_num_g"] = clinicalNumStran(config).get(var["clinical_significance"], 3) if var["clinical_significance"] and var["clinical_significance"] != "-" else \
							  clinicalNumStran(config).get(var["function_classification"], 3)
		var["clinic_num_s"] = functionNumStran(config).get(var["function_classification"], 3) if var["function_classification"] and var["function_classification"] != "-" else \
							  functionNumStran(config).get(var["clinical_significance"], 3)
		
		# 2025.12.10-安徽省立BRCA需要展示I/II类，clinic_num_s根据治疗方案等级重新设置下
		var["clinic_num_s"], var["top_level"] = S_function(var)

	for gene in ["BRCA1", "BRCA2"]:
		for level in [1,2,3,4,5]:
			result["B"+gene[-1]+"_gcnv_L"+str(level)] = [var for var in cnv if var["gene_symbol"] == gene and var["clinic_num_g"] == level]
	return result

# 2025.03.11-增加规则-适用150-BRCA1/2、林奇5个基因和VHL大片段缺失/重复突变
def process_gcnv_allgene(jsonDict, config):
	result = {}
	cnv = [var for var in copy.deepcopy(jsonDict["cnv"]) if var["var_origin"] == "germline"]
	cnv = sorted(cnv, key = lambda i:i["gene_symbol"])
	for var in cnv:
		var["evi_sum"] = varRegimen(jsonDict, var["evi_sum"], config, var)
		var["var_info_forFDZS"] = varInfo_FDZS(var["gene_symbol"], "", "", config)
		# 字段名跟mlpa保持一致，这样模板少改一些
		var["type"] = "Gain" if var["cnv_type"] == "Gain" else "Loss"
		var["value"] = var["region_exon"]
		var["clinic_num_g"] = clinicalNumStran(config).get(var["clinical_significance"], 3) if var["clinical_significance"] and var["clinical_significance"] != "-" else \
							  clinicalNumStran(config).get(var["function_classification"], 3)
		var["clinic_num_s"] = functionNumStran(config).get(var["function_classification"], 3) if var["function_classification"] and var["function_classification"] != "-" else \
							  functionNumStran(config).get(var["clinical_significance"], 3)
	for level in [1,2,3,4,5]:
		result["level_"+str(level)] = [var for var in cnv if var["clinic_num_g"] == level]

	return result