#-*- coding:gbk -*-
import re
from libs.getConfig import clinicalNumStran, functionNumStran
from libs.getEvi import varRegimen
from libs.specialRequest import varInfo_FJZL
from libs.rule import S_function
import copy
from libs.rule import decimal_float, decimal_percen
from libs.getConfig import get_gene_class

'''
Discription
	
	处理cnv格式。 (非胚系变异)

'''

def process_cnv(jsonDict, config):
	# 2024.08.30-CNV过滤掉germline的变异
	#cnv = copy.deepcopy(jsonDict["cnv"])
	cnv = [var for var in copy.deepcopy(jsonDict["cnv"]) if var["var_origin"] != "germline"]
	# 2024.08.30-更新完成
	for var in cnv:
		#var["cn_mean"] = format(round(float(var["cn_mean"]), 2)+0.00, ".2f") if "cn_mean" in var.keys() and var["cn_mean"] else 0
		var["cn_mean"] = decimal_float(var["cn_mean"]) if "cn_mean" in var.keys() and var["cn_mean"] else 0
		# 变异分类加个临时功能，后续极元开发完成后再禁用（只改SNVindel和CNV，SV的默认按体细胞的方案进行判定，CNV由于有BRCA MLPA，也加一下吧）-2023.05.23
		#if var["gene_symbol"] in get_gene_class(config)["function"]:
		#	var["clinic_num_g"] = clinicalNumStran(config).get(var["function_classification"], 3)
		#	var["clinic_num_s"] = functionNumStran(config).get(var["function_classification"], 3)
		#elif var["gene_symbol"] in get_gene_class(config)["clinical"]:
		#	var["clinic_num_g"] = clinicalNumStran(config).get(var["clinical_significance"], 3)
		#	var["clinic_num_s"] = functionNumStran(config).get(var["clinical_significance"], 3)
		#else:
			# 变异分类更新，致癌性和致病性只返回1个-2023.05.22
			#var["clinic_num_g"] = clinicalNumStran().get(var["clinical_significance"], 3) 
			#var["clinic_num_s"] = functionNumStran().get(var["function_classification"], 3)
		#	var["clinic_num_g"] = clinicalNumStran(config).get(var["clinical_significance"], 3) if var["clinical_significance"] and var["clinical_significance"] != "-" else \
		#						  clinicalNumStran(config).get(var["function_classification"], 3)
		#	var["clinic_num_s"] = functionNumStran(config).get(var["function_classification"], 3) if var["function_classification"] and var["function_classification"] != "-" else \
		#						  functionNumStran(config).get(var["clinical_significance"], 3)
			# 2023.05.22 更新完成

		# 临床意义更新-2023.11.21
		var["clinic_num_g"] = clinicalNumStran(config).get(var["clinical_significance"], 3) if var["clinical_significance"] and var["clinical_significance"] != "-" else \
							  clinicalNumStran(config).get(var["function_classification"], 3)
		var["clinic_num_s"] = functionNumStran(config).get(var["function_classification"], 3) if var["function_classification"] and var["function_classification"] != "-" else \
							  functionNumStran(config).get(var["clinical_significance"], 3)
		# 临床意义更新完成-2023.11.21
		
		# 2023.05.23 更新完成
		var["evi_sum"] = varRegimen(jsonDict, var["evi_sum"], config, var)
		var["clinic_num_s"], var["top_level"] = S_function(var)
		# 福建肿瘤：返回变异频率相关信息（来源配置表）和治疗方案汇总
		var["var_info_forFJZL"], var["var_regimen_forFJZL"] = varInfo_FJZL(var, jsonDict["sample_info"]["tumor_list"], config)
	# 按cn_mean排序下
	cnv = sorted(cnv, key=lambda i:float(i["cn_mean"]), reverse=True)

	return cnv