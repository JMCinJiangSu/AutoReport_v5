#-*- coding:gbk -*-
import copy
import re
from libs import listResultToDict
from libs.getPinyin import topinyin
'''
Discription
	
	治疗方案先后顺序制定成排序规则（包含单药和联合治疗） 

'''

def sort_rule(jsonDict):
	#gss_copy = copy.deepcopy(jsonDict["gss"]) if "gss" in jsonDict.keys() and jsonDict["gss"] else {}
	# gss应该返回个字典，新版本变成了列表，这边加个兼容-2023.09.26
	gss_copy = listResultToDict.ListToDict(copy.deepcopy(jsonDict["gss"])) if "gss" in jsonDict.keys() and jsonDict["gss"] else {}
	
	# 新增HD-2025.09.16
	# 判断HD是否要展示
	hd = [var for var in copy.deepcopy(jsonDict["hd"])] if "hd" in jsonDict.keys() else []
	judge_hd_inter = False
	for var in hd:
		if "evi_sum" in var.keys():
			judge_hd_inter = True
			break
	# 2025.08.19-进院北三CP200要展示HD
	# 2025.11.14-新增Master Panel（组织_北医）
	# 2026.01.12-新增BHD组织
	# 2026.01.28-BHD改名BRCA（前列腺癌）
	hd_copy = []
	if judge_hd_inter and (jsonDict["sample_info"]["prod_names"] in ["Master Panel（组织）", "Master Panel", "HANDLE OncoPro", "HANDLE OncoPro Compact", "Master Panel（组织_北医）", "BHD（组织）", "BRCA（前列腺癌）"] or 
						(jsonDict["sample_info"]["prod_names"] in ["OncoPro（组织）", "Classic Panel 200（组织）"] and jsonDict["sample_info"]["company"] == "北京大学第三医院" and jsonDict["sample_info"]["report_module_type"] == "hospital")):
		hd_copy = copy.deepcopy(jsonDict["hd"])
	# 新增完成

	# 兼容完成-2023.09.26
	var_list = copy.deepcopy(jsonDict["snvindel"]+\
							 jsonDict["cnv"]+\
							 jsonDict["sv"]+\
							 jsonDict["rna_sv"]+\
							 jsonDict["knb"]+\
							 jsonDict["msi"]+\
							 jsonDict["tmb"]+\
							 jsonDict["pdl1"]+\
							 jsonDict["mlpa"]+\
							 jsonDict["hrd"]+\
							 [gss_copy]+\
							 hd_copy)
	regimen_list = []
	for var in var_list:
		if "evi_sum" in var.keys():
			# 2025.05.14-新增排序规则-v4返回的治疗方案是随机顺序，程序会重新排序，这边按通用规则也排一下
			# 通用规则：等级A>B>C>D，按拼音顺序
			for evi in var["evi_sum"]:
				evi["evi_conclusion_simple"] = evi["evi_conclusion"][0] if evi["evi_conclusion"] else ""
				evi["regimen_name_py"] = topinyin(evi["regimen_name"]) if evi["regimen_name"] else "0"
			var["evi_sum"] = sorted(var["evi_sum"], key = lambda i : (i["evi_conclusion_simple"], i["regimen_name_py"].upper()))
			# 2025.05.14-新增完成
			regimen_list += [evi["regimen_name"] for evi in var["evi_sum"] if evi["regimen_name"] and re.search("Sensitive", evi["clinical_significance"])]

	return regimen_list