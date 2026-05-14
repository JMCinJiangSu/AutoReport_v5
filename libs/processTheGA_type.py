#-*- coding:gbk -*-
import re
from libs.getEvi import varRegimen
from libs import listResultToDict
from libs.rule import judge_var, get_varsimpleinfo
import copy

'''
Discription
	
	处理胃癌分子分型格式。
	需要ebv_type、MSI和部分基因检测结果 

'''

def process_ga_type(jsonDict, var_data):
	ga_result = {}
	ebv_gene_list = ["PIK3CA", "ARID1A", "BCOR"]
	gs_gene_list = ["CDH1", "RHOA", "CLDN18"]
	cin_gene_list = ["TP53", "ERBB2", "KRAS", "EGFR", "CDK6", "CCNE1", "APC", "CTNNB1", "SMAD2", "SMAD4", "PTEN"]
	ebv_type_dict = listResultToDict.ListToDict(copy.deepcopy(jsonDict["ebv_type"])) if "ebv_type" in jsonDict.keys() and jsonDict["ebv_type"] else {}
	# 胃癌分子分型有：EBV、MSI、GS、CIN四种，胃癌分子分型相关标志物部分检测到的分型都展示
	# 检测小结部分
	# 检测到EBV、MSI则展示EBV、MSI（同时存在则都展示）
	# 未检测到EBV、MSI时，若存在GS、CIN分型，则展示（同时存在则都展示）
	### 变异未限定变异类型
	# 展示体细胞I/II/肿瘤发生发展相关变异和胚系致病/疑似致病变异
	level_12_var = [var for var in var_data if judge_var(var, [4,5], [4,5])]
	# 输入基因列表，返回检测结果
	def get_var(gene_list):
		result_dict = {}
		for var in level_12_var:
			for gene in re.split(",", var["gene_symbol"]):
				if gene in gene_list:
					if gene not in result_dict.keys():
						result_dict.setdefault(gene, [])
					var_info = get_varsimpleinfo(var)
					result_dict[gene].append(var_info)
		#return result_dict, ", ".join(["{0} {1}".format(k, v) for k, v in result_dict.items()])
		return result_dict, ", ".join((["{0} {1}".format(k, a) for k,v in result_dict.items() for a in v]))

	ebv_gene_result, ebv_gene_str = get_var(ebv_gene_list)
	gs_gene_list, gs_gene_str = get_var(gs_gene_list)
	cin_gene_list, cin_gene_str = get_var(cin_gene_list)
	ga_result = {
		"ebv_type" : ebv_type_dict,
		"ebv_gene" : ebv_gene_result,
		"ebv_sum" : ebv_gene_str,
		"gs_gene" : gs_gene_list,
		"gs_sum" : gs_gene_str,
		"cin_gene" : cin_gene_list,
		"cin_sum" : cin_gene_str
	}
	
	return ga_result

# 2025.04.02-保留胃癌相关变异的所有信息
def cp200_process_ga_type(jsonDict, var_data):
	ga_result = {}
	ebv_gene_list = ["PIK3CA", "ARID1A"]
	gs_gene_list = ["CDH1", "CLDN18"]
	cin_gene_list = ["TP53", "ERBB2", "KRAS", "EGFR", "CDK6", "CCNE1", "APC", "CTNNB1", "SMAD4", "PTEN"]
	ebv_type_dict = listResultToDict.ListToDict(copy.deepcopy(jsonDict["ebv_type"])) if "ebv_type" in jsonDict.keys() and jsonDict["ebv_type"] else {}
	# 胃癌分子分型有：EBV、MSI、GS、CIN四种，胃癌分子分型相关标志物部分检测到的分型都展示
	# 检测小结部分
	# 检测到EBV、MSI则展示EBV、MSI（同时存在则都展示）
	# 未检测到EBV、MSI时，若存在GS、CIN分型，则展示（同时存在则都展示）
	### 变异未限定变异类型
	# 展示体细胞I/II/肿瘤发生发展相关变异和胚系致病/疑似致病变异
	level_12_var = [var for var in var_data if judge_var(var, [4,5], [4,5])]
	# 输入基因列表，返回检测结果
	def get_var(gene_list):
		result_dict = {}
		result_dict_all = {}
		for var in level_12_var:
			for gene in re.split(",", var["gene_symbol"]):
				if gene in gene_list:
					if gene not in result_dict.keys():
						result_dict.setdefault(gene, [])
					if gene not in result_dict_all.keys():
						result_dict_all.setdefault(gene, [])
					var_info = get_varsimpleinfo(var)
					result_dict[gene].append(var_info)
					result_dict_all[gene].append(var)
		#return result_dict, ", ".join(["{0} {1}".format(k, v) for k, v in result_dict.items()])
		return result_dict_all, ", ".join((["{0} {1}".format(k, a) for k,v in result_dict.items() for a in v]))

	ebv_gene_result, ebv_gene_str = get_var(ebv_gene_list)
	gs_gene_list, gs_gene_str = get_var(gs_gene_list)
	cin_gene_list, cin_gene_str = get_var(cin_gene_list)
	ga_result = {
		"ebv_type" : ebv_type_dict,
		"ebv_gene" : ebv_gene_result,
		"ebv_sum" : ebv_gene_str,
		"gs_gene" : gs_gene_list,
		"gs_sum" : gs_gene_str,
		"cin_gene" : cin_gene_list,
		"cin_sum" : cin_gene_str
	}
	
	return ga_result

# 2026.02.04-氨基酸改为三字母，适配武汉协和MP
def process_ga_type_abbr(jsonDict, var_data):
	ga_result = {}
	ebv_gene_list = ["PIK3CA", "ARID1A", "BCOR"]
	gs_gene_list = ["CDH1", "RHOA", "CLDN18"]
	cin_gene_list = ["TP53", "ERBB2", "KRAS", "EGFR", "CDK6", "CCNE1", "APC", "CTNNB1", "SMAD2", "SMAD4", "PTEN"]
	ebv_type_dict = listResultToDict.ListToDict(copy.deepcopy(jsonDict["ebv_type"])) if "ebv_type" in jsonDict.keys() and jsonDict["ebv_type"] else {}
	# 胃癌分子分型有：EBV、MSI、GS、CIN四种，胃癌分子分型相关标志物部分检测到的分型都展示
	# 检测小结部分
	# 检测到EBV、MSI则展示EBV、MSI（同时存在则都展示）
	# 未检测到EBV、MSI时，若存在GS、CIN分型，则展示（同时存在则都展示）
	### 变异未限定变异类型
	# 展示体细胞I/II/肿瘤发生发展相关变异和胚系致病/疑似致病变异
	level_12_var = [var for var in var_data if judge_var(var, [4,5], [4,5])]
	# 输入基因列表，返回检测结果
	def get_var(gene_list):
		result_dict = {}
		for var in level_12_var:
			for gene in re.split(",", var["gene_symbol"]):
				if gene in gene_list:
					if gene not in result_dict.keys():
						result_dict.setdefault(gene, [])
					var_info = get_varsimpleinfo_abbr(var)
					result_dict[gene].append(var_info)
		#return result_dict, ", ".join(["{0} {1}".format(k, v) for k, v in result_dict.items()])
		return result_dict, ", ".join((["{0} {1}".format(k, a) for k,v in result_dict.items() for a in v]))

	ebv_gene_result, ebv_gene_str = get_var(ebv_gene_list)
	gs_gene_list, gs_gene_str = get_var(gs_gene_list)
	cin_gene_list, cin_gene_str = get_var(cin_gene_list)
	ga_result = {
		"ebv_type" : ebv_type_dict,
		"ebv_gene" : ebv_gene_result,
		"ebv_sum" : ebv_gene_str,
		"gs_gene" : gs_gene_list,
		"gs_sum" : gs_gene_str,
		"cin_gene" : cin_gene_list,
		"cin_sum" : cin_gene_str
	}
	
	return ga_result

# 返回变异结果，如Snvindel返回hgvs_p_abbr或hgvs_c
def get_varsimpleinfo_abbr(var):
	if var["bio_category"] == "Snvindel":
		return var["hgvs_p_abbr"] if var["hgvs_p"] != "p.?" else var["hgvs_c"]
	elif var["bio_category"] == "Cnv":
		return "扩增"
	elif var["bio_category"] in ["Sv", "PSeqRnaSv"]:
		return "{0}-{1}融合".format(var["five_prime_gene"], var["three_prime_gene"])
	# 20225.06.17-新增HD
	elif var["bio_category"] == "PHd":
		if var["type"] == "HomoDel":
			return "纯合缺失"
		elif var["type"] == "HeteDel":
			return "杂合缺失"
		else:
			return "未知变异类型！"