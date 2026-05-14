#-*- coding:gbk -*-
import re
from libs.getConfig import clinicalNumStran, functionNumStran
from libs.getEvi import varRegimen
from libs.rule import S_function
import copy
from libs.specialRequest import varInfo_FJZL

'''
Discription
	
	处理hd格式。 (非胚系变异)

'''

def process_hd(jsonDict, config):
	hd = [var for var in copy.deepcopy(jsonDict["hd"])] if "hd" in jsonDict.keys() else []
	# 判断HD是否要展示
	judge_hd_inter = False
	for var in hd:
		if "evi_sum" in var.keys():
			judge_hd_inter = True
			break
	# 2025.08.19-进院北三CP200要展示HD
	# 2025.11.14-新增Master Panel（组织_北医）
	# 2026.01.12-新增BHD组织
	# 2026.01.28-BHD改名BRCA（前列腺癌）
	# 2026.02.10-西安交大一BHD产品名改为BRCA1/BRCA2（组织）、BRCA1/BRCA2（组织 全血）
	# 2026.03.04-温附一-JY和进院CP200要展示HD
	# 2026.04.07-德阳-JY和进院CP200要展示HD
	# 2026.04.21-南充市中心医院-JYCP200要展示HD
	if judge_hd_inter and (jsonDict["sample_info"]["prod_names"] in ["Master Panel（组织）", "Master Panel", "HANDLE OncoPro", "HANDLE OncoPro Compact", "Master Panel（组织_北医）", "BHD（组织）", "BRCA（前列腺癌）"] or \
						(jsonDict["sample_info"]["prod_names"] in ["OncoPro（组织）", "Classic Panel 200（组织）"] and jsonDict["sample_info"]["company"] == "北京大学第三医院" and jsonDict["sample_info"]["report_module_type"] == "hospital") or \
						(jsonDict["sample_info"]["prod_names"] in ["BRCA1/BRCA2（组织）", "BRCA1/BRCA2（组织 全血）"] and jsonDict["sample_info"]["company"] in ["西安交通大学医学院第一附属医院", "西安交通大学第一附属医院"]) or \
						(jsonDict["sample_info"]["prod_names"] in ["OncoPro（组织）", "Classic Panel 200（组织）"] and jsonDict["sample_info"]["company"] == "温州医科大学附属第一医院" and jsonDict["sample_info"]["report_module_type"] == "hospital") or \
						(jsonDict["sample_info"]["prod_names"] in ["OncoPro（组织）", "Classic Panel 200（组织）"] and jsonDict["sample_info"]["origin_company"] in ["温州医科大学附属第一医院-JY", "广东省胸部疾病学会"] and jsonDict["sample_info"]["report_module_type"] == "rummage") or \
						(jsonDict["sample_info"]["prod_names"] in ["OncoPro（组织）", "Classic Panel 200（组织）"] and jsonDict["sample_info"]["company"] == "温州医科大学附属第一医院" and jsonDict["sample_info"]["report_module_type"] == "rummage" and jsonDict["sample_info"]["order_type"] == "汇总进院") or \
						(jsonDict["sample_info"]["prod_names"] in ["OncoPro（组织）", "Classic Panel 200（组织）"] and jsonDict["sample_info"]["company"] == "德阳市人民医院" and jsonDict["sample_info"]["report_module_type"] == "hospital") or \
						(jsonDict["sample_info"]["prod_names"] in ["OncoPro（组织）", "Classic Panel 200（组织）"] and jsonDict["sample_info"]["origin_company"] in ["德阳市人民医院-JY"] and jsonDict["sample_info"]["report_module_type"] == "rummage") or \
						(jsonDict["sample_info"]["prod_names"] in ["OncoPro（组织）", "Classic Panel 200（组织）"] and jsonDict["sample_info"]["origin_company"] in ["南充市中心医院-JY"] and jsonDict["sample_info"]["report_module_type"] == "rummage")):
		for var in hd:
			var["clinic_num_g"] = clinicalNumStran(config).get(var["clinical_significance"], 3) if var["clinical_significance"] and var["clinical_significance"] != "-" else \
								  clinicalNumStran(config).get(var["function_classification"], 3)
			var["clinic_num_s"] = functionNumStran(config).get(var["function_classification"], 3) if var["function_classification"] and var["function_classification"] != "-" else \
								  functionNumStran(config).get(var["clinical_significance"], 3)
			var["evi_sum"] = varRegimen(jsonDict, var["evi_sum"], config, var)
			var["clinic_num_s"], var["top_level"] = S_function(var)
			# 福建肿瘤：返回变异频率相关信息（来源配置表）和治疗方案汇总
			var["var_info_forFJZL"], var["var_regimen_forFJZL"] = varInfo_FJZL(var, jsonDict["sample_info"]["tumor_list"], config)
		return hd
	else:
		return []
	