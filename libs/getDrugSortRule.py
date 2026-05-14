#-*- coding:gbk -*-
import copy
import re
from libs import listResultToDict
from libs.getPinyin import topinyin
'''
Discription
	
	药物先后顺序制定成排序规则（治疗方案中若有多种药物，则拆分处理） 

'''

def sort_rule(jsonDict):
	rule = []
	#gss_copy = copy.deepcopy(jsonDict["gss"]) if "gss" in jsonDict.keys() and jsonDict["gss"] else {}
	# gss应该返回个字典，新版本变成了列表，这边加个兼容-2023.09.26
	gss_copy = listResultToDict.ListToDict(copy.deepcopy(jsonDict["gss"])) if "gss" in jsonDict.keys() and jsonDict["gss"] else {}
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
							 [gss_copy])
	drug_list = []
	for var in var_list:
		if "evi_sum" in var.keys():
			# 2025.05.14-新增排序规则-v4返回的治疗方案是随机顺序，程序会重新排序，这边按通用规则也排一下
			# 通用规则：等级A>B>C>D，按拼音顺序
			for evi in var["evi_sum"]:
				evi["evi_conclusion_simple"] = evi["evi_conclusion"][0] if evi["evi_conclusion"] else ""
				evi["regimen_name_py"] = topinyin(evi["regimen_name"]) if evi["regimen_name"] else "0"
			var["evi_sum"] = sorted(var["evi_sum"], key = lambda i : (i["evi_conclusion_simple"], i["regimen_name_py"].upper()))
			# 2025.05.14-新增完成
			drug_list += [evi["regimen_name"] for evi in var["evi_sum"] if evi["regimen_name"] and re.search("Sensitive", evi["clinical_significance"])]
	for drug in drug_list:
		for i in re.split("\+", drug):
			if i.strip() not in rule:
				rule.append(i.strip())
	return rule