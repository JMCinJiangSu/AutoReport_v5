#-*- coding:gbk -*-
import jinja2
import re
from docxtpl import RichText
from collections import Counter

import jinja2.filters
from libs.rule import S_level
import copy
from pypinyin import pinyin, Style
from itertools import chain
from datetime import datetime, timedelta
from itertools import groupby
from collections import Counter

'''
自定义jinja2过滤器
'''

### 定制类 ###
# 1. 复旦中山厦门医院CP40，药物“、”.join()展示为一行
def regimen_sum(a):
	result = []
	evi_sum = a["evi_sum"]
	if "evi_split" in evi_sum.keys():
		if "Predictive" in evi_sum["evi_split"].keys():
			for i in evi_sum["evi_split"]["Predictive"]:
				result.append("{0}（{1}，{2}级）".format(i["regimen_name"], i["clinical_significance_cn"], i["evi_conclusion_simple"]))
		if "Prognostic" in evi_sum["evi_split"].keys():
			for i in evi_sum["evi_split"]["Prognostic"]:
				result.append("{0}（{1}，{2}级）".format("预后"+i["clinical_significance_cn"], " / ", i["evi_conclusion_simple"]))
		if "Diagnostic" in evi_sum["evi_split"].keys():
			for i in evi_sum["evi_split"]["Diagnostic"]:
				result.append("{0}（{1}，{2}级）".format("辅助诊断", " / ", i["evi_conclusion_simple"]))
	return "、".join(result) if result else "-"
jinja2.filters.FILTERS["regimen_sum_filter"] = regimen_sum 

# 2. 北大人民检验科HRR-出BRCA报告-检测结论
# gene_symbol gene_region hgvs_c hgvs_p, clinic_num_g.strans
def var_brca_sum(var_list):
	result = []
	clinic_trans = {5 : "致病性变异", 4 : "疑似致病性变异"}
	for var in var_list:
		if var["type"] == "Loss":
			result.append("{0} {1} del, {2}".format(var["gene_symbol"], var["value"], "疑似致病性变异"))

		else:
			if var["hgvs_p"] != "p.?":
				result.append("{0} {1} {2} {3}, {4}".format(var["gene_symbol"], var["gene_region"], var["hgvs_c"], var["hgvs_p"], clinic_trans[var["clinic_num_g"]]))
			else:
				result.append("{0} {1} {2}, {3}".format(var["gene_symbol"], var["gene_region"], var["hgvs_c"], clinic_trans[var["clinic_num_g"]]))
	return "；".join(result)
jinja2.filters.FILTERS["var_brca_sum"] = var_brca_sum

# 2. 北大人民检验科HRR-出BRCA报告-检测结论-兼容CNV/MLPA
# gene_symbol gene_region hgvs_c hgvs_p, clinic_num_g.strans
def var_brca_sum_v2(var_list):
	result = []
	clinic_trans = {5 : "致病性变异", 4 : "疑似致病性变异"}
	for var in var_list:
		if var["type"] == "Loss":
			result.append("{0} {1} del, {2}".format(var["gene_symbol"], var["value"], clinic_trans[var["clinic_num_g"]]))
		elif var["type"] == "Gain":
			result.append("{0} {1} dup, {2}".format(var["gene_symbol"], var["value"], clinic_trans[var["clinic_num_g"]]))

		else:
			if var["hgvs_p"] != "p.?":
				result.append("{0} {1} {2} {3}, {4}".format(var["gene_symbol"], var["gene_region"], var["hgvs_c"], var["hgvs_p"], clinic_trans[var["clinic_num_g"]]))
			else:
				result.append("{0} {1} {2}, {3}".format(var["gene_symbol"], var["gene_region"], var["hgvs_c"], clinic_trans[var["clinic_num_g"]]))
	return "；".join(result)
jinja2.filters.FILTERS["var_brca_sum_v2"] = var_brca_sum_v2

# 3. 北大人民检验科HRR-出BRCA报告-变异解读结论
def inter_brca_sum(var_list):
	result_dict = {}
	result = []
	for var in var_list:
		if var["type"] == "Loss":
			var["clinic_num_g"] = 4
		elif var["type"] == "Gain":
			var["clinic_num_g"] = 3

		if var["clinic_num_g"] not in result_dict.keys():
			result_dict.setdefault(var["clinic_num_g"], [])
		result_dict[var["clinic_num_g"]].append(var)
	clinic_trans = {5:"致病性变异", 4:"疑似致病性变异", 3:"意义不明确变异", 2:"疑似良性变异"}
	for i in [5,4,3,2]:
		if i in result_dict.keys():
			result.append(str(len(result_dict[i]))+"个"+clinic_trans[i])
	return "、".join(result)
jinja2.filters.FILTERS["inter_brca_sum"] = inter_brca_sum

# 3. 北大人民检验科HRR-出BRCA报告-变异解读结论-2024.12.27
def inter_brca_sum_v2(var_list):
	result_dict = {}
	result = []
	for var in var_list:
		#if var["type"] == "Loss":
		#	var["clinic_num_g"] = 4
		#elif var["type"] == "Gain":
		#	var["clinic_num_g"] = 3

		if var["clinic_num_g"] not in result_dict.keys():
			result_dict.setdefault(var["clinic_num_g"], [])
		result_dict[var["clinic_num_g"]].append(var)
	clinic_trans = {5:"致病性变异", 4:"疑似致病性变异", 3:"意义不明确变异", 2:"疑似良性变异"}
	for i in [5,4,3,2]:
		if i in result_dict.keys():
			result.append(str(len(result_dict[i]))+"个"+clinic_trans[i])
	return "、".join(result)
jinja2.filters.FILTERS["inter_brca_sum_v2"] = inter_brca_sum_v2

# 孙逸仙116-结果小结
def summary_SYX_116(var_list):
	result = []
	for var in var_list:
		if var["bio_category"] == "Snvindel":
			if var["hgvs_p"] != "p.?":
				result.append(var["gene_symbol"]+" "+var["hgvs_p"])
			else:
				result.append(var["gene_symbol"]+" "+var["hgvs_c"])
		elif var["bio_category"] == "Cnv":
			result.append(var["gene_symbol"]+" 扩增")
		elif var["bio_category"] == "Sv":
			result.append(var["five_prime_gene"]+"-"+var["three_prime_gene"]+" 融合")
	return ", ".join(result)
jinja2.filters.FILTERS["summary_SYX_116"] = summary_SYX_116

# 中山人民HRDC检测结果小结-2023.06.02
def summary_ZSRM_hrd(info):
	# 阳性：基于本次送检样本，检出BRCA1基因BRCA1 gene_region hgvs_c hgvs_p I类变异，HRD状态阳性以及TP53基因TP53 gene_region hgvs_c hgvs_p II类变异。
	# 阴性：基于本次送检样本，未检出致病性或可能致病性变异。
	var_dict = info[0]
	hrd = info[1]
	class_stran = {5 : "I类", 4 : "II类"}
	
	brca_result = []
	for var in var_dict["ec_type"]["BRCA1_level12"] + var_dict["ec_type"]["BRCA2_level12"]:
		if var["hgvs_p"] != "p.?":
			brca_result.append("{0}基因{0} {1} {2} {3} {4}变异".format(var["gene_symbol"], var["gene_region"], var["hgvs_c"], var["hgvs_p"], class_stran.get(var["clinic_num_s"])))
		else:
			brca_result.append("{0}基因{0} {1} {2} {3}变异".format(var["gene_symbol"], var["gene_region"], var["hgvs_c"], class_stran.get(var["clinic_num_s"])))
	brca_result_str = "，".join(brca_result)
	
	hrd_result = "HRD状态阳性" if hrd["var_id"] == "HRD+" else "HRD状态阴性"

	other_result = []
	for var in var_dict["var_somatic"]["level_I"] + var_dict["var_somatic"]["level_II"]:
		if var["gene_symbol"] not in ["BRCA1", "BRCA2"]:
			if var["hgvs_p"] != "p.?":
				other_result.append("{0}基因{0} {1} {2} {3} {4}变异".format(var["gene_symbol"], var["gene_region"], var["hgvs_c"], var["hgvs_p"], class_stran.get(var["clinic_num_s"])))
			else:
				other_result.append("{0}基因{0} {1} {2} {3}变异".format(var["gene_symbol"], var["gene_region"], var["hgvs_c"], class_stran.get(var["clinic_num_s"])))
	other_result_str = ",".join(other_result)
	hrd_result_str = "以及".join([hrd_result, other_result_str])

	result = []
	if brca_result_str:
		result.append(brca_result_str)
	if hrd_result_str:
		result.append(hrd_result_str)

	# 组装-阳性/阴性
	if var_dict["var_somatic"]["level_I"] + var_dict["var_somatic"]["level_II"]:
		return "基于本次送检样本，检出{0}。".format("，".join(result))
	else:
		return "基于本次送检样本，未检出致病性或疑似致病性变异。"
jinja2.filters.FILTERS["summary_ZSRM_hrd"] = summary_ZSRM_hrd

##############################################################################################################

### 通用类 ###

### 变异列表类 ###
# 1. 体细胞：返回肿瘤发生发展相关变异列表
def return_somatic_onconodrug_list(info):
	var = info[0]
	sample = info[1]
	if re.search("Master", sample["prod_names"]):
		if not sample["control_sample_id"]:
			return var["var_somatic"]["level_onco_nodrug"] + var["var_germline_nodrug"]
		else:
			return var["var_somatic"]["level_onco_nodrug"]
	else:
		return var["var_somatic"]["level_onco_nodrug"]
jinja2.filters.FILTERS["somatic_onconodrug_list"] = return_somatic_onconodrug_list

# 2. 体细胞：返回体细胞/未知来源变异，III类变异列表
def return_somatic_level3_list(info):
	var = info[0]
	sample = info[1]
	if sample["prod_names"] == "BPTM（5基因）":
		return var["ec_type"]["POLE_level3"]+var["ec_type"]["TP53_level3"]+var["ec_type"]["BRCA1_level3"]+var["ec_type"]["BRCA2_level3"]
	elif sample["prod_names"] == "PTM（3基因）":
		return var["ec_type"]["POLE_level3"]+var["ec_type"]["TP53_level3"]
	else:
		return var["var_somatic"]["level_III"]
jinja2.filters.FILTERS["somatic_level3_list"] = return_somatic_level3_list

# 3. 返回林奇变异列表，未检测到变异的基因也要展示
def return_gLS5_list(gLS5):
	result = []
	gene_list = ["EPCAM", "MLH1", "MSH2", "MSH6", "PMS2"]
	for gene in gene_list:
		if gene in gLS5.keys() and gLS5[gene]:
			result.extend(gLS5[gene])
		else:
			result.append({"gene_symbol" : gene})
	return result
jinja2.filters.FILTERS["gLS5_list"] = return_gLS5_list
	
# 4. 返回BRCA基因检测结果，未检测到变异的基因也要展示
def return_brca_list(brca):
	result = []
	brca1_list = brca["snv_s"]["B1_L5"] + brca["snv_s"]["B1_L4"]+brca["mlpa"]["B1_Loss"]+brca["mlpa"]["B1_Gain"]+brca["snv_s"]["B1_L3"]
	brca2_list = brca["snv_s"]["B2_L5"] + brca["snv_s"]["B2_L4"]+brca["mlpa"]["B2_Loss"]+brca["mlpa"]["B2_Gain"]+brca["snv_s"]["B2_L3"]
	if brca1_list:
		result.extend(brca1_list)
	else:
		result.append({"gene_symbol" : "BRCA1"})
	
	if brca2_list:
		result.extend(brca2_list)
	else:
		result.append({"gene_symbol" : "BRCA2"})
	return result
jinja2.filters.FILTERS["brca_list"] = return_brca_list

# 4. 返回BRCA基因检测结果，未检测到变异的基因也要展示-2024.12.10-区分MLPA和CNV，并且变异等级不再del为4类，dup为3类
def return_brca_list_v2(info):
	brca = info[0]
	var_type = info[1]
	result = []
	if var_type == "mlpa":
		brca1_list = brca["snv_s"]["B1_L5"] + brca["mlpa_v2"]["B1_mlpa_L5"] + \
					 brca["snv_s"]["B1_L4"] + brca["mlpa_v2"]["B1_mlpa_L4"] + \
					 brca["snv_s"]["B1_L3"] + brca["mlpa_v2"]["B1_mlpa_L3"]
		brca2_list = brca["snv_s"]["B2_L5"] + brca["mlpa_v2"]["B2_mlpa_L5"] + \
					 brca["snv_s"]["B2_L4"] + brca["mlpa_v2"]["B2_mlpa_L4"] + \
					 brca["snv_s"]["B2_L3"] + brca["mlpa_v2"]["B2_mlpa_L3"]
	else:
		brca1_list = brca["snv_s"]["B1_L5"] + brca["gcnv_v2"]["B1_gcnv_L5"] + \
					 brca["snv_s"]["B1_L4"] + brca["gcnv_v2"]["B1_gcnv_L4"] + \
					 brca["snv_s"]["B1_L3"] + brca["gcnv_v2"]["B1_gcnv_L3"]
		brca2_list = brca["snv_s"]["B2_L5"] + brca["gcnv_v2"]["B2_gcnv_L5"] + \
					 brca["snv_s"]["B2_L4"] + brca["gcnv_v2"]["B2_gcnv_L4"] + \
					 brca["snv_s"]["B2_L3"] + brca["gcnv_v2"]["B2_gcnv_L3"]

	if brca1_list:
		result.extend(brca1_list)
	else:
		result.append({"gene_symbol" : "BRCA1"})
	
	if brca2_list:
		result.extend(brca2_list)
	else:
		result.append({"gene_symbol" : "BRCA2"})
	return result
jinja2.filters.FILTERS["brca_list_v2"] = return_brca_list_v2
	

#----------------------------------------------------------------------------------------------------------------

### 变异解读类 ###
# 1. 体细胞：返回体细胞/未知来源变异 I、II类解读变异列表
def return_somatic_level12_inter(info):
	var_brca = info[0]
	var = info[1]
	sample = info[2]
	if re.search("BRCA", sample["prod_names"]):
		# var_origin != "germline"
		brca_result = []
		for var in var_brca["snv_s"]["B1_L5"] + var_brca["snv_s"]["B2_L5"] + var_brca["snv_s"]["B1_L4"] + var_brca["snv_s"]["B2_L4"]:
			if var["var_origin"] != "germline":
				brca_result.append(var)
		return brca_result
	elif re.search("BPTM（5基因）", sample["prod_names"]):
		return var["ec_type"]["POLE_level12"] + var["ec_type"]["TP53_level12"] + var["ec_type"]["BRCA1_level12"] + var["ec_type"]["BRCA2_level12"]
	elif re.search("PTM（3基因）", sample["prod_names"]):
		return var["ec_type"]["POLE_level12"] + var["ec_type"]["TP53_level12"]
	# 如果是MP且无配对时，var_origin为germline的且有用药的也要展示
	elif re.search("Master Panel（组织）|Master Panel（组织-综合）|Master Panel V1（组织）", sample["prod_names"]) and not sample["control_sample_id"]:
		return var["var_somatic"]["level_I"] + var["var_germline"]["regimen_level_I"] + var["var_somatic"]["level_II"] + var["var_germline"]["regimen_level_II"]
	# HRD 完整通用模板，过滤掉BRCA突变-2023.12.06-暂时不用
	#elif sample["prod_names"] == "HRD Complete（组织）" and "complete" in sample["report_name"]:
	#	return [var for var in var["var_somatic"]["level_I"] + var["var_somatic"]["level_II"] if var["gene_symbol"] not in ["BRCA1", "BRCA2"]]
	# HRD 添加完成-2023.12.06
	else:
		return var["var_somatic"]["level_I"] + var["var_somatic"]["level_II"]
jinja2.filters.FILTERS["somatic_level12_inter"] = return_somatic_level12_inter

# 2. 胚系：返回胚系变异4/5类列表，用于解读
def return_germline_level45_inter(info):
	var_brca = info[0]
	var = info[1]
	sample = info[2]
	if re.search("BRCA", sample["prod_names"]):
		return var_brca["snv_m"]["B1_G_L5"]+var_brca["snv_m"]["B2_G_L5"]+var_brca["snv_m"]["B1_G_L4"]+var_brca["snv_m"]["B2_G_L4"]+\
			   var_brca["mlpa"]["B1_Loss"] + var_brca["mlpa"]["B2_Loss"]
	elif re.search("HRR", sample["prod_names"]):
		return var["var_germline"]["level_5"]+var["var_germline"]["level_4"]+var_brca["mlpa"]["B1_Loss"]+var_brca["mlpa"]["B2_Loss"]
	else:
		return var["var_germline"]["level_5"]+var["var_germline"]["level_4"]
jinja2.filters.FILTERS["germline_level45_inter"] = return_germline_level45_inter

# 2. 胚系：返回胚系变异4/5类列表，用于解读-v2
def return_germline_level45_inter_v2(info):
	var_brca = info[0]
	var = info[1]
	sample = info[2]
	var_type = info[3]
	if re.search("BRCA", sample["prod_names"]):
		if var_type == "mlpa":
			return var_brca["snv_m"]["B1_G_L5"] + var_brca["snv_m"]["B2_G_L5"] + \
				   var_brca["mlpa_v2"]["B1_mlpa_L5"] + var_brca["mlpa_v2"]["B2_mlpa_L5"] + \
				   var_brca["snv_m"]["B1_G_L4"]+var_brca["snv_m"]["B2_G_L4"] + \
				   var_brca["mlpa_v2"]["B1_mlpa_L4"] + var_brca["mlpa_v2"]["B2_mlpa_L4"]
		else:
			return var_brca["snv_m"]["B1_G_L5"] + var_brca["snv_m"]["B2_G_L5"] + \
				   var_brca["gcnv_v2"]["B1_gcnv_L5"] + var_brca["gcnv_v2"]["B2_gcnv_L5"] + \
				   var_brca["snv_m"]["B1_G_L4"]+var_brca["snv_m"]["B2_G_L4"] + \
				   var_brca["gcnv_v2"]["B1_gcnv_L4"] + var_brca["gcnv_v2"]["B2_gcnv_L4"]
		
	elif re.search("HRR", sample["prod_names"]):
		if var_type == "mlpa":
			return var["var_germline"]["level_5"] + var_brca["mlpa_v2"]["B1_mlpa_L5"] + var_brca["mlpa_v2"]["B2_mlpa_L5"] + \
				   var["var_germline"]["level_4"] + var_brca["mlpa_v2"]["B1_mlpa_L4"] + var_brca["mlpa_v2"]["B2_mlpa_L4"]
		else:
			return var["var_germline"]["level_5"] + var_brca["gcnv_v2"]["B1_gcnv_L5"] + var_brca["gcnv_v2"]["B2_gcnv_L5"] + \
				   var["var_germline"]["level_4"] + var_brca["gcnv_v2"]["B1_gcnv_L4"] + var_brca["gcnv_v2"]["B2_gcnv_L4"]
		
	else:
		return var["var_germline"]["level_5"]+var["var_germline"]["level_4"]
jinja2.filters.FILTERS["germline_level45_inter_v2"] = return_germline_level45_inter_v2

#3. 胚系：返回胚系变异3类列表，用于解读
def return_germline_level3_inter(info):
	var_brca = info[0]
	var = info[1]
	sample = info[2]
	if re.search("BRCA", sample["prod_names"]):
		return var_brca["mlpa"]["B1_Gain"]+var_brca["mlpa"]["B2_Gain"]+var_brca["snv_m"]["B1_G_L3"]+var_brca["snv_m"]["B2_G_L3"]
	elif re.search("HRR", sample["prod_names"]):
		return var_brca["mlpa"]["B1_Gain"]+var_brca["mlpa"]["B2_Gain"]+var["var_germline"]["level_3"]
	else:
		return var["var_germline"]["level_3"]
jinja2.filters.FILTERS["germline_level3_inter"] = return_germline_level3_inter

#3. 胚系：返回胚系变异3类列表，用于解读-v2
def return_germline_level3_inter_v2(info):
	var_brca = info[0]
	var = info[1]
	sample = info[2]
	var_type = info[3]
	if re.search("BRCA", sample["prod_names"]):
		if var_type == "mlpa":
			return var_brca["snv_m"]["B1_G_L3"] + var_brca["snv_m"]["B2_G_L3"] + \
				   var_brca["mlpa_v2"]["B1_mlpa_L3"] + var_brca["mlpa_v2"]["B2_mlpa_L3"]
		else:
			return var_brca["snv_m"]["B1_G_L3"] + var_brca["snv_m"]["B2_G_L3"] + \
				   var_brca["gcnv_v2"]["B1_gcnv_L3"] + var_brca["gcnv_v2"]["B2_gcnv_L3"]
		
	elif re.search("HRR", sample["prod_names"]):
		if var_type == "mlpa":
			return var["var_germline"]["level_3"] + var_brca["mlpa_v2"]["B1_mlpa_L3"] + var_brca["mlpa_v2"]["B2_mlpa_L3"]
		else:
			return var["var_germline"]["level_3"] + var_brca["gcnv_v2"]["B1_gcnv_L3"] + var_brca["gcnv_v2"]["B2_gcnv_L3"]
	else:
		return var["var_germline"]["level_3"]
jinja2.filters.FILTERS["germline_level3_inter_v2"] = return_germline_level3_inter_v2

#----------------------------------------------------------------------------------------------------------------

### 变异结果中的内容 ###
# 1. 返回基因，单个或多个
def return_gene_symbol(var):
	result = []
	if "bio_category" in var.keys() and var["bio_category"]:
		if var["bio_category"] in ["Snvindel", "Cnv"]:
			result.append(var["gene_symbol"])
		elif var["bio_category"] in ["Sv", "PSeqRnaSv"]:
			if "," in var["gene_symbol"]:
				if var["five_prime_gene"] != var["three_prime_gene"]:
					result.append(var["five_prime_gene"])
					result.append(var["three_prime_gene"])
				else:
					result.append(var["five_prime_gene"])
			else:
				result.append(var["gene_symbol"])
		else:
			result.append(var["gene_symbol"])
	else:
		# gene_symbol 可能不存在，这边加个兼容-2023.09.22
		if "gene_symbol" in var.keys():
			result.append(var["gene_symbol"])
		# 兼容完成-2023.09.22
	rt = RichText()
	rt.add("\n".join(result), italic=True)
	return rt
jinja2.filters.FILTERS["gene_symbol"] = return_gene_symbol

# 1. 返回基因，单个或多个-2025.03.24
def return_gene_symbol_v2(var):
	result = []
	if "bio_category" in var.keys() and var["bio_category"]:
		if var["bio_category"] in ["Snvindel", "Cnv"]:
			result.append(var["gene_symbol"])
		elif var["bio_category"] in ["Sv", "PSeqRnaSv"]:
			if "," in var["gene_symbol"]:
				if var["five_prime_gene"] != var["three_prime_gene"]:
					result.append(var["five_prime_gene"])
					result.append(var["three_prime_gene"])
				else:
					result.append(var["five_prime_gene"])
			else:
				result.append(var["gene_symbol"])
		else:
			result.append(var["gene_symbol"])
	else:
		# gene_symbol 可能不存在，这边加个兼容-2023.09.22
		if "gene_symbol" in var.keys():
			result.append(var["gene_symbol"])
		# 兼容完成-2023.09.22
	#rt = RichText()
	#rt.add("\n".join(result), italic=True)
	return "\n".join(result)
jinja2.filters.FILTERS["gene_symbol_v2"] = return_gene_symbol_v2

# 2. 变异检测结果
def var_info(var):
	result = []
	if "bio_category" in var.keys():
		if var["bio_category"] == "Snvindel":
			if var["hgvs_p"] != "p.?":
				result.append(var["gene_region"]+" "+var["hgvs_c"]+" "+var["hgvs_p"])
			else:
				result.append(var["gene_region"]+" "+var["hgvs_c"])
			result.append(var["transcript_primary"])
		elif var["bio_category"] == "Cnv":
			result.append("扩增") 
		elif var["bio_category"] in ["Sv", "PSeqRnaSv"]:
			if var["five_prime_gene"] == "MET" and var["three_prime_gene"] == "MET":
				result.append("MET exon14 跳跃")
				result.append(var["five_prime_transcript"])
			else:
				result.append("{0}:{1}-{2}:{3}融合".format(var["five_prime_gene"], var["five_prime_cds"], var["three_prime_gene"], var["three_prime_cds"]))
				result.append("{0}/{1}".format(var["five_prime_transcript"], var["three_prime_transcript"]))
		elif var["bio_category"] == "PMLPA":
			if var["type"] == "Loss":
				result.append(var["value"]+" del")
			elif var["type"] == "Gain":
				result.append(var["value"]+" dup")
		# 2025.07.01-新增PHd
		elif var["bio_category"] == "PHd":
			if var["type"] == "HomoDel":
				result.append("纯合缺失")
			elif var["type"] == "HeteDel":
				result.append("杂合缺失")
			else:
				result.append("未知变异类型！")
			result.append("（{0}）".format(var["region"]))
		# 2025.07.01-新增完成
	rt = RichText()
	rt.add("\n".join(result))
	return rt
	#return result
jinja2.filters.FILTERS["var_info"] = var_info

# 2. 变异检测结果-v2
def var_info_v2(var):
	result = []
	if "bio_category" in var.keys():
		if var["bio_category"] == "Snvindel":
			if var["hgvs_p"] != "p.?":
				result.append(var["gene_region"]+" "+var["hgvs_c"]+" "+var["hgvs_p"])
			else:
				result.append(var["gene_region"]+" "+var["hgvs_c"])
			result.append(var["transcript_primary"])
		elif var["bio_category"] == "Cnv":
			if var["var_origin"] == "germline" and var["gene_symbol"] in ["BRCA1", "BRCA2"]:
				if var["type"] == "Loss":
					result.append(var["value"] + " del")
				else:
					result.append(var["value"] + " dup")
			else:
				result.append("扩增") 
		elif var["bio_category"] in ["Sv", "PSeqRnaSv"]:
			if var["five_prime_gene"] == "MET" and var["three_prime_gene"] == "MET":
				result.append("MET exon14 跳跃")
				result.append(var["five_prime_transcript"])
			else:
				result.append("{0}:{1}-{2}:{3}融合".format(var["five_prime_gene"], var["five_prime_cds"], var["three_prime_gene"], var["three_prime_cds"]))
				result.append("{0}/{1}".format(var["five_prime_transcript"], var["three_prime_transcript"]))
		elif var["bio_category"] == "PMLPA":
			if var["type"] == "Loss":
				result.append(var["value"]+" del")
			elif var["type"] == "Gain":
				result.append(var["value"]+" dup")
	rt = RichText()
	rt.add("\n".join(result))
	return rt
	#return result
jinja2.filters.FILTERS["var_info_v2"] = var_info_v2

# 2. 变异检测结果-v3,2025.03.24
def var_info_v3(var):
	result = []
	if "bio_category" in var.keys():
		if var["bio_category"] == "Snvindel":
			if var["hgvs_p"] != "p.?":
				result.append(var["gene_region"]+" "+var["hgvs_c"]+" "+var["hgvs_p"])
			else:
				result.append(var["gene_region"]+" "+var["hgvs_c"])
			result.append(var["transcript_primary"])
		elif var["bio_category"] == "Cnv":
			if var["var_origin"] == "germline" and var["gene_symbol"] in ["BRCA1", "BRCA2"]:
				if var["type"] == "Loss":
					result.append(var["value"] + " del")
				else:
					result.append(var["value"] + " dup")
			else:
				result.append("扩增") 
		elif var["bio_category"] in ["Sv", "PSeqRnaSv"]:
			if var["five_prime_gene"] == "MET" and var["three_prime_gene"] == "MET":
				result.append("MET exon14 跳跃")
				result.append(var["five_prime_transcript"])
			else:
				result.append("{0}:{1}-{2}:{3}融合".format(var["five_prime_gene"], var["five_prime_cds"], var["three_prime_gene"], var["three_prime_cds"]))
				result.append("{0}/{1}".format(var["five_prime_transcript"], var["three_prime_transcript"]))
		elif var["bio_category"] == "PMLPA":
			if var["type"] == "Loss":
				result.append(var["value"]+" del")
			elif var["type"] == "Gain":
				result.append(var["value"]+" dup")
		# 2025.07.17-新增PHd
		elif var["bio_category"] == "PHd":
			if var["type"] == "HomoDel":
				result.append("纯合缺失")
			elif var["type"] == "HeteDel":
				result.append("杂合缺失")
			else:
				result.append("未知变异类型！")
			result.append("（{0}）".format(var["region"]))
		# 2025.07.17-新增完成
	#rt = RichText()
	#rt.add("\n".join(result))
	return "\n".join(result)
	#return result
jinja2.filters.FILTERS["var_info_v3"] = var_info_v3


# 3. 返回变异丰度，适用BRCA等，胚系展示纯合/杂合，体细胞展示频率，MLPA用-
def freq_stran(info):
	'''
	体细胞snvindel：freq_str
	胚系snvindel：厦门项目freq， 上海项目freq_rc，注意MP还要考虑是否有对照样本，无对照样本，来源为germline的直接展示freq_str
	sv：CP40 copies，其他freq_str
	PSeqRnaSv：freq
	'''
	var = info[0]
	sample = info[1]
	if "var_origin" in var.keys():
		if var["bio_category"] == "Snvindel":
			if var["var_origin"] == "germline":
				if "Master" in sample["prod_names"]:
					if sample["control_sample_id"]:
						return "纯合" if var["freq_rc"] and float(var["freq_rc"]) >= 0.85 else "杂合"
					else:
						return var["freq_str"]
				elif re.search("116|76|25|21|18", sample["prod_names"]):
					return "纯合" if var["freq_rc"] and float(var["freq_rc"]) >= 0.85 else "杂合" if var["freq_rc"] and float(var["freq_rc"]) < 0.85 else "未提取到freq_rc！"
				else:
					return "纯合" if float(var["freq"]) >= 0.85 else "杂合"
			else:
				return var["freq_str"]
		elif var["bio_category"] == "Cnv":
			return var["cn_mean"]
		elif var["bio_category"] == "Sv":
			if re.search("Classic|CRC12", sample["prod_names"]):
				return str(var["copies"])+" copies"
			else:
				return var["freq_str"]
		elif var["bio_category"] == "PSeqRnaSv":
			return str(var["freq"])+" copies"
		# 2025.07.01-新增PHd
		elif var["bio_category"] == "PHd":
			return "-"
		# 2025.07.01-新增完成
	else:
		if var["type"] in ["Loss", "Gain"]:
			return "/"
jinja2.filters.FILTERS["freq_stran"] = freq_stran

# 3. 返回变异丰度，适用BRCA等，胚系展示纯合/杂合，体细胞展示频率，MLPA用-，新增gCNV-2024.12.26
def freq_stran_v2(info):
	'''
	体细胞snvindel：freq_str
	胚系snvindel：厦门项目freq， 上海项目freq_rc，注意MP还要考虑是否有对照样本，无对照样本，来源为germline的直接展示freq_str
	sv：CP40 copies，其他freq_str
	PSeqRnaSv：freq
	'''
	var = info[0]
	sample = info[1]
	if "var_origin" in var.keys():
		if var["bio_category"] == "Snvindel":
			if var["var_origin"] == "germline":
				if "Master" in sample["prod_names"]:
					if sample["control_sample_id"]:
						return "纯合" if float(var["freq_rc"]) >= 0.85 else "杂合"
					else:
						return var["freq_str"]
				elif re.search("116|76|25|21|18", sample["prod_names"]):
					return "纯合" if var["freq_rc"] and float(var["freq_rc"]) >= 0.85 else "杂合" if var["freq_rc"] and float(var["freq_rc"]) < 0.85 else "未提取到freq_rc！"
				else:
					return "纯合" if float(var["freq"]) >= 0.85 else "杂合"
			else:
				return var["freq_str"]
		elif var["bio_category"] == "Cnv":
			if var["var_origin"] == "germline" and var["gene_symbol"] in ["BRCA1", "BRCA2"]:
				return "杂合" if var["cnv_type"] == "HeteDel" else "纯合" if var["cnv_type"] == "HomoDel" else "/"
			else:
				return var["cn_mean"]
		elif var["bio_category"] == "Sv":
			if re.search("Classic|CRC12", sample["prod_names"]):
				return str(var["copies"])+" copies"
			else:
				return var["freq_str"]
		elif var["bio_category"] == "PSeqRnaSv":
			return str(var["freq"])+" copies"
	else:
		if var["type"] in ["Loss", "Gain"]:
			return "/"
jinja2.filters.FILTERS["freq_stran_v2"] = freq_stran_v2

# 4. 临床意义-药物
def significance_regimen(var):
	result = []
	#rt = RichText()
	evi_sum = var["evi_sum"]
	if evi_sum and "evi_split" in evi_sum.keys():
		if "Predictive" in evi_sum["evi_split"].keys():
			for i in evi_sum["evi_split"]["Predictive"]:
				result.append("{0}（{1}，{2}级）".format(i["regimen_name"], i["clinical_significance_cn"], i["evi_conclusion_simple"]))
		if "Prognostic" in evi_sum["evi_split"].keys():
			for i in evi_sum["evi_split"]["Prognostic"]:
				result.append("{0}（{1}，{2}级）".format("预后"+i["clinical_significance_cn"], " / ", i["evi_conclusion_simple"]))
		if "Diagnostic" in evi_sum["evi_split"].keys():
			for i in evi_sum["evi_split"]["Diagnostic"]:
				result.append("{0}（{1}，{2}级）".format("辅助诊断", " / ", i["evi_conclusion_simple"]))
		#rt.add("\n".join(result)) 
	if not result:
		result = ["-"]
	rt = RichText()
	rt.add("\n".join(result))
	return rt
	#return result
jinja2.filters.FILTERS["significance_regimen"] = significance_regimen

# 4. 临床意义-药物-v2-取消超文本-2025.12.24
def significance_regimen_v2(var):
	result = []
	evi_sum = var["evi_sum"]
	if evi_sum and "evi_split" in evi_sum.keys():
		if "Predictive" in evi_sum["evi_split"].keys():
			for i in evi_sum["evi_split"]["Predictive"]:
				result.append("{0}（{1}，{2}级）".format(i["regimen_name"], i["clinical_significance_cn"], i["evi_conclusion_simple"]))
		if "Prognostic" in evi_sum["evi_split"].keys():
			for i in evi_sum["evi_split"]["Prognostic"]:
				result.append("{0}（{1}，{2}级）".format("预后"+i["clinical_significance_cn"], " / ", i["evi_conclusion_simple"]))
		if "Diagnostic" in evi_sum["evi_split"].keys():
			for i in evi_sum["evi_split"]["Diagnostic"]:
				result.append("{0}（{1}，{2}级）".format("辅助诊断", " / ", i["evi_conclusion_simple"]))
		#rt.add("\n".join(result)) 
	if not result:
		result = ["-"]
	return "\n".join(result)
	#return result
jinja2.filters.FILTERS["significance_regimen_v2"] = significance_regimen_v2

# 5. 返回临床意义，适用BRCA等
def clinic_stran(info):
	var = info[0]
	s_type = info[1]
	g_dict = {
		5 : "致病性",
		4 : "疑似致病性",
		3 : "意义不明确",
		2 : "疑似良性",
		1 : "良性"
	}
	s_dict1 = {
		5 : "I类-强临床意义",
		4 : "II类-潜在临床意义",
		3 : "III类-临床意义不明",
		2 : "IV类-良性/可能良性",
		1 : "IV类-良性/可能良性"
	}
	s_dict2 = {
		5 : "I类",
		4 : "II类",
		3 : "III类",
		2 : "IV类",
		1 : "IV类"
	}
	s_dict = s_dict1 if s_type == "s1" else s_dict2 if s_type == "s2" else {}
	if "var_origin" in var.keys():
		if var["var_origin"] == "germline":
			return g_dict.get(var["clinic_num_g"], "")
		else:
			# 肿瘤发生发展相关归为III类
			if var["clinic_num_s"] in [4, 5] and (not var["evi_sum"]["evi_split"] or \
				(var["evi_sum"]["evi_split"] and not set(["Diagnostic","Predictive","Prognostic"]) & set(var["evi_sum"]["evi_split"].keys()))) :
				return s_dict.get(3)
			else:
				return s_dict.get(var["clinic_num_s"], "")
	else:
		if var["type"] == "Loss":
			return g_dict.get(4)
		elif var["type"] == "Gain":
			return g_dict.get(3)
jinja2.filters.FILTERS["clinic_stran"] = clinic_stran

# 5. 返回临床意义，适用BRCA等-大片段重复/缺失等级改为动态的-2024.12.26
def clinic_stran_v2(info):
	var = info[0]
	s_type = info[1]
	g_dict = {
		5 : "致病性",
		4 : "疑似致病性",
		3 : "意义不明确",
		2 : "疑似良性",
		1 : "良性"
	}
	s_dict1 = {
		5 : "I类-强临床意义",
		4 : "II类-潜在临床意义",
		3 : "III类-临床意义不明",
		2 : "IV类-良性/可能良性",
		1 : "IV类-良性/可能良性"
	}
	s_dict2 = {
		5 : "I类",
		4 : "II类",
		3 : "III类",
		2 : "IV类",
		1 : "IV类"
	}
	s_dict = s_dict1 if s_type == "s1" else s_dict2 if s_type == "s2" else {}
	if "var_origin" in var.keys():
		if var["var_origin"] == "germline":
			return g_dict.get(var["clinic_num_g"], "")
		else:
			# 肿瘤发生发展相关归为III类
			if var["clinic_num_s"] in [4, 5] and (not var["evi_sum"]["evi_split"] or \
				(var["evi_sum"]["evi_split"] and not set(["Diagnostic","Predictive","Prognostic"]) & set(var["evi_sum"]["evi_split"].keys()))) :
				return s_dict.get(3)
			else:
				return s_dict.get(var["clinic_num_s"], "")
	else:
		if var["type"] == "Loss":
			return g_dict.get(var["clinic_num_g"], "")
		elif var["type"] == "Gain":
			return g_dict.get(var["clinic_num_g"], "")
jinja2.filters.FILTERS["clinic_stran_v2"] = clinic_stran_v2

# 5.1 返回临床意义，适用MP
def clinic_stran_MP(info):
	var = info[0]
	s_type = info[1]
	sample = info[2]
	g_dict = {
		5 : "致病性",
		4 : "疑似致病性",
		3 : "意义不明确",
		2 : "疑似良性",
		1 : "良性"
	}
	s_dict1 = {
		5 : "I类-强临床意义",
		4 : "II类-潜在临床意义",
		3 : "III类-临床意义不明",
		2 : "IV类-良性/可能良性",
		1 : "IV类-良性/可能良性"
	}
	s_dict2 = {
		5 : "I类",
		4 : "II类",
		3 : "III类",
		2 : "IV类",
		1 : "IV类"
	}
	s_dict = s_dict1 if s_type == "s1" else s_dict2 if s_type == "s2" else {}
	if "var_origin" in var.keys():
		if var["var_origin"] == "germline":
			if sample["control_sample_id"]:
				return g_dict.get(var["clinic_num_g"], "")
			else:
				return s_dict.get(var["clinic_num_s"], "")
		else:
			# 肿瘤发生发展相关归为III类
			if var["clinic_num_s"] in [4, 5] and (not var["evi_sum"]["evi_split"] or \
				(var["evi_sum"]["evi_split"] and not set(["Diagnostic","Predictive","Prognostic"]) & set(var["evi_sum"]["evi_split"].keys()))) :
				return s_dict.get(var[3])
			else:
				return s_dict.get(var["clinic_num_s"], "")
	else:
		if var["type"] == "Loss":
			return g_dict.get(4)
		elif var["type"] == "Gain":
			return g_dict.get(3)
jinja2.filters.FILTERS["clinic_stran_MP"] = clinic_stran_MP


# 6. 变异来源转化
def var_origin_stran(var):
	if "var_origin" in var.keys():
		if var["var_origin"] == "germline":
			return "胚系"
		elif var["var_origin"] == "somatic":
			return "体细胞"
		else:
			return "待定"
	else:
		if var["type"] in ["Loss", "Gain"]:
			return "胚系"
jinja2.filters.FILTERS["var_origin_stran"] = var_origin_stran

# 7. 结果小结表格表头-丰度-适用于BRCA和林奇
def title_freq_brca(prod_names):
	if prod_names in ["BRCA1/BRCA2（全血）", "BRCA1/2（扩增子）", "HRR（全血）", "林奇综合征"]:
		return  "基因型"
	elif prod_names in ["BRCA1/BRCA2（组织）"]:
		return "丰度"
	elif prod_names in ["BRCA1/BRCA2（组织 全血）"]:
		return "基因型/丰度"
jinja2.filters.FILTERS["title_freq_brca"] = title_freq_brca

# 8. III 类变异表头-丰度-体细胞
# 2024.05.07新增BPTM Plus（组织）
def title_freq_III(sample):
	if sample["prod_names"] in ["BRCA1/BRCA2（全血）", "BRCA1/2（扩增子）", "HRR（全血）", "林奇综合征"]:
		return  "基因型"
	elif sample["prod_names"] in ["BRCA1/BRCA2（组织）","BPTM（5基因）", "PTM（3基因）", "HRR（组织）", "HRD Complete（组织）", "10基因（血液）", \
								  "Pan116（血液）", "LC76（血液）", "CRC25（血液）", "TC21（血液）", "GA18（血液）", "HRR（组织）-前列腺", \
								  "Pan116（血液）-胰腺", "Pan116（血液）-乳腺", "BPTM Plus（组织）"]:
		return "丰度"
	elif sample["prod_names"] in ["BRCA1/BRCA2（组织 全血）", "HRR（组织 全血）"]:
		return "丰度"
	elif sample["prod_names"] in ["10基因（组织）", "Classic Panel", "CRC12-MSI",\
								  "Pan116（组织）", "LC76（组织）", "CRC25（组织）", "TC21（组织）", "GA18（组织）",\
								  "Master Panel（组织）", "Master Panel（血液）",\
								  "Classic Panel-胃", "Classic Panel-胃肠", "Classic Panel-甲状腺", "Classic Panel-肝胆", \
								  "Classic Panel-骨", "Classic Panel-黑色素", "Classic Panel-头颈",\
								  "Pan116（组织）-胰腺", "Pan116（组织）-乳腺", "Master Panel（血液-综合）", "Master Panel（组织-综合）", "Master Panel V1（组织）"]:
		return "丰度/拷贝数"
jinja2.filters.FILTERS["title_freq_III"] = title_freq_III

# 9. 肿瘤发生发展相关变异-丰度
# 目前涉及的项目有HRR组织/HRR配对/HRD/LC10组织/LC10血液/PAN116组织/PAN116血液/LC76组织/LC76血液/CRC25组织/CRC25血液/TC21组织/TC21血液/GA18组织/GA18血液
# Master组织/Master血液/CP40/CRC12
# 丰度/拷贝数
# 2024.05.07新增BPTM Plus（组织）
# 2026.03.11新增WES（组织）
def title_freq_onconodrug(sample):
	if sample["prod_names"] in ["HRR（组织）", "HRR（组织 全血）", "HRD Complete（组织）", "10基因（血液）", \
							   "Pan116（血液）", "LC76（血液）", "CRC25（血液）", "TC21（血液）", "GA18（血液）", "HRR（组织）-前列腺", \
							   "Pan116（血液）-胰腺", "Pan116（血液）-乳腺", "BPTM Plus（组织）"]:
		return "丰度"
	elif sample["prod_names"] in ["Pan116（组织）", "LC76（组织）", "CRC25（组织）", "TC21（组织）", "GA18（组织）",\
								 "Master Panel（组织）", "Master Panel（血液）", "Classic Panel", "CRC12-MSI", "10基因（组织）", \
								 "Classic Panel-胃", "Classic Panel-胃肠", "Classic Panel-甲状腺", "Classic Panel-肝胆", \
								 "Classic Panel-骨", "Classic Panel-黑色素", "Classic Panel-头颈",\
								 "Pan116（组织）-胰腺", "Pan116（组织）-乳腺", "Master Panel（血液-综合）", "Master Panel（组织-综合）", "Master Panel V1（组织）", \
								 "WES（组织）"]:
		return "丰度/拷贝数"
jinja2.filters.FILTERS["title_freq_onconodrug"] = title_freq_onconodrug

# 10. 靶向治疗表-丰度
# 适用于MP/116/LC10/CP40/CRC12
def title_freq_targetRegimen(sample):
	freq = ["丰度"]
	if re.search("组织", sample["prod_names"]):
		freq.append("拷贝数")
	if sample["prod_names"] in ["Classic Panel", "CRC12-MSI"]:
		freq.append("拷贝数")
	if sample["prod_names"] == "Master Panel（血液）" or sample["prod_names"] == "Master Panel（血液-综合）":
		freq.append("拷贝数")
	if sample["control_sample_id"]:
		freq.append("基因型")
	return "/".join(freq)
jinja2.filters.FILTERS["title_freq_targetRegimen"] = title_freq_targetRegimen

# 11 PARP结果汇总表-丰度表头
def title_parp(sample):
	if sample["prod_names"] in ["BRCA1/BRCA2（全血）", "HRR（全血）"]:
		return "基因型"
	elif sample["prod_names"] in ["HRR（组织）", "HRD Complete（组织）"]:
		return "丰度"
	elif sample["prod_names"] in ["HRR（组织 全血）"]:
		return "基因型/丰度"
jinja2.filters.FILTERS["title_parp"] = title_parp
#----------------------------------------------------------------------------------------------------------------

### 结果小结类 ###
# 1. PD-L1
def pdl1_summary(pdl1):
	return "阴性。" if pdl1["result"] == "阴性" else "阳性，{0}为{1}。".format(pdl1["type"], pdl1["value"])
jinja2.filters.FILTERS["pdl1_summary"] = pdl1_summary

# 2. MSI
def msi_summary(msi):
	return "微卫星稳定（MSS）。" if msi["var_id"] == "MSS" else "微卫星不稳定（MSI-H）。" if msi["var_id"] == "MSI-H" else "未返回MSI结果！"
jinja2.filters.FILTERS["msi_summary"] = msi_summary

# 3. TMB
def tmb_summary(tmb):
	TMB_result = "低" if tmb["var_id"] == "TMB-L" else "高"
	return "{0} Muts/Mb, 肿瘤突变负荷较{1}（{2}）。".format(tmb["TMB_value"], TMB_result, "TMB-H" if tmb["var_id"] == "TMB-H" else "TMB-L")
jinja2.filters.FILTERS["tmb_summary"] = tmb_summary

# 4. GEP
def gep_summary(gep):
	return "GEP分值为{0}分".format(gep["gep_score"])
jinja2.filters.FILTERS["gep_summary"] = gep_summary

# 5. TME
def tme_summary(tme):
	tme_dict = {
		"IE/F" : "免疫富集/纤维化亚型(IE/F)",
		"IE" : "免疫富集/非纤维化亚型(IE)",
		"F" : "纤维化亚型(F)",
		"D" : "免疫荒漠型(D)"
	}
	return "TME分型为{0}".format(tme_dict.get(tme["tme_type"], tme["tme_type"]))
jinja2.filters.FILTERS["tme_summary"] = tme_summary

# 6. 免疫检查点抑制剂疗效相关基因
def io_summary(io):
	result = []
	if io["io_p_summary"]:
		result.append(io["io_p_summary"]+"（疗效正相关）")
	if io["io_n_summary"]:
		result.append(io["io_n_summary"]+"（疗效负相关）")
	if not result:
		result = ["未检出相关变异"]
	rt = RichText()
	rt.add("；\n".join(result)+"。")
	return rt
jinja2.filters.FILTERS["io_summary"] = io_summary

# 6. 免疫检查点抑制剂疗效相关基因-2024.02.20更新，目前只有MP和116用到
def io_summary_v2(info):
	io = info[0]
	sample = info[1]
	result = []
	if sample["prod_names"] in ["Master Panel（组织）", "Master Panel（血液）", "Master Panel V1（组织）"]:
		if io["io_p_summary_new_cnvlist"]:
			result.append(io["io_p_summary_new_cnvlist"]+"（疗效正相关）")
		if io["io_n_summary_new_cnvlist"]:
			result.append(io["io_n_summary_new_cnvlist"]+"（疗效负相关）")
	else:
		if io["io_p_summary_116"]:
			result.append(io["io_p_summary_116"]+"（疗效正相关）")
		if io["io_n_summary_116"]:
			result.append(io["io_n_summary_116"]+"（疗效负相关）")
	if not result:
		result = ["未检出相关变异"]
	rt = RichText()
	rt.add("；\n".join(result)+"。")
	return rt
jinja2.filters.FILTERS["io_summary_v2"] = io_summary_v2

# 6. 免疫检查点抑制剂疗效相关基因-2025.03.24更新，目前只有MP和116用到
# 不使用超文本了，模板中使用replace("\n", "<w:br/>")进行换行符替换
def io_summary_v3(info):
	io = info[0]
	sample = info[1]
	result = []
	if sample["prod_names"] in ["Master Panel（组织）", "Master Panel（血液）", "Master Panel V1（组织）"]:
		if io["io_p_summary_new_cnvlist"]:
			result.append(io["io_p_summary_new_cnvlist"]+"（疗效正相关）")
		if io["io_n_summary_new_cnvlist"]:
			result.append(io["io_n_summary_new_cnvlist"]+"（疗效负相关）")
	else:
		if io["io_p_summary_116"]:
			result.append(io["io_p_summary_116"]+"（疗效正相关）")
		if io["io_n_summary_116"]:
			result.append(io["io_n_summary_116"]+"（疗效负相关）")
	if not result:
		result = ["未检出相关变异"]
	return "；\n".join(result)+"。"
jinja2.filters.FILTERS["io_summary_v3"] = io_summary_v3

# 6. 免疫检查点抑制剂疗效相关基因-2025.07.01更新，目前只有MP和116用到-MP新增HD
def io_summary_v4(info):
	io = info[0]
	sample = info[1]
	result = []
	if sample["prod_names"] in ["Master Panel（组织）", "Master Panel（血液）", "Master Panel V1（组织）"]:
		if io["hd_io_p_summary"]:
			result.append(io["hd_io_p_summary"]+"（疗效正相关）")
		if io["hd_io_n_summary"]:
			result.append(io["hd_io_n_summary"]+"（疗效负相关）")
	else:
		if io["io_p_summary_116"]:
			result.append(io["io_p_summary_116"]+"（疗效正相关）")
		if io["io_n_summary_116"]:
			result.append(io["io_n_summary_116"]+"（疗效负相关）")
	if not result:
		result = ["未检出相关变异"]
	return "；\n".join(result)+"。"
jinja2.filters.FILTERS["io_summary_v4"] = io_summary_v4

# 识别是否为数值
def is_number(i):
	try:
		float(i)
		return True
	except:
		pass 
	if i.isnumeric():
		return True
	return False

# 7. HRD检测结果
def hrd_summary(info):
	gss = info[0]
	sample = info[1]
	result = []
	hrd_result = "HRD阳性" if gss["BRCA1"] or gss["BRCA2"] or float(gss["gss"]["gsscore"]) >= 50 else "HRD阴性"
	# 加结果仅供参考的一些情况-2023.08.29
	# 未检测到BRCA致病/疑似致病变异+以下任意一个条件
	# 1）baf_noise > 0.055
	# 2) depth_noise > 0.35
	# 3) 无肿瘤细胞含量
	# 4）有肿瘤细胞含量但是格式不正确的
	# 5）有肿瘤细胞含量，格式正确，但是数值<30
	note = "（结果仅供参考）" if not gss["BRCA1"] and not gss["BRCA2"] and \
							(float(gss["gss"]["baf_noise"]) > 0.055 or \
							 float(gss["gss"]["depth_noise"]) > 0.35 or \
							 not sample["tumor_content"] or \
							 (sample["tumor_content"] and not is_number(sample["tumor_content_num"])) or \
							 (sample["tumor_content"] and is_number(sample["tumor_content_num"]) and float(sample["tumor_content_num"]) < 30)) \
							else ""
	result.append(hrd_result + note)
	if gss["summary"]:
		result.append("HRR通路相关基因突变：{0}".format(gss["summary"]))
	rt = RichText()
	rt.add("；\n".join(result)+"。")
	return rt
jinja2.filters.FILTERS["hrd_summary"] = hrd_summary

# 7. HRD检测结果-v2,20240402
# baf_noise/depth_noise 改为判定var_auto_result，gss阈值改为45
def hrd_summary_v2(info):
	gss = info[0]
	sample = info[1]
	result = []
	hrd_result = "HRD阳性" if gss["BRCA1"] or gss["BRCA2"] or float(gss["gss"]["gsscore"]) >= 45 else "HRD阴性"
	# 加结果仅供参考的一些情况-2023.08.29
	# 未检测到BRCA致病/疑似致病变异+以下任意一个条件
	# 1）baf_noise > 0.055（取消，2024.04.02）
	# 2) depth_noise > 0.35（取消，2024.04.02）
	# 3) 无肿瘤细胞含量
	# 4）有肿瘤细胞含量但是格式不正确的
	# 5）有肿瘤细胞含量，格式正确，但是数值<30
	# 6）var_auto_result为F（新增，2024.04.02）
	note = "（结果仅供参考）" if not gss["BRCA1"] and not gss["BRCA2"] and \
							(("var_auto_result" in gss["gss"].keys() and gss["gss"]["var_auto_result"] and gss["gss"]["var_auto_result"] == "F") or \
							 not sample["tumor_content"] or \
							 (sample["tumor_content"] and not is_number(sample["tumor_content_num"])) or \
							 (sample["tumor_content"] and is_number(sample["tumor_content_num"]) and float(sample["tumor_content_num"]) < 30)) \
							else ""
	result.append(hrd_result + note)
	if gss["summary"]:
		result.append("HRR通路相关基因突变：{0}".format(gss["summary"]))
	rt = RichText()
	rt.add("；\n".join(result)+"。")
	return rt
jinja2.filters.FILTERS["hrd_summary_v2"] = hrd_summary_v2

# 7. HRD检测结果-v3,2025.09.22
# baf_noise/depth_noise 改为判定var_auto_result，gss阈值改为45
# 不使用富文本
def hrd_summary_v3(info):
	gss = info[0]
	sample = info[1]
	result = []
	hrd_result = "HRD阳性" if gss["BRCA1"] or gss["BRCA2"] or float(gss["gss"]["gsscore"]) >= 45 else "HRD阴性"
	# 加结果仅供参考的一些情况-2023.08.29
	# 未检测到BRCA致病/疑似致病变异+以下任意一个条件
	# 1）baf_noise > 0.055（取消，2024.04.02）
	# 2) depth_noise > 0.35（取消，2024.04.02）
	# 3) 无肿瘤细胞含量
	# 4）有肿瘤细胞含量但是格式不正确的
	# 5）有肿瘤细胞含量，格式正确，但是数值<30
	# 6）var_auto_result为F（新增，2024.04.02）
	note = "（结果仅供参考）" if not gss["BRCA1"] and not gss["BRCA2"] and \
							(("var_auto_result" in gss["gss"].keys() and gss["gss"]["var_auto_result"] and gss["gss"]["var_auto_result"] == "F") or \
							 not sample["tumor_content"] or \
							 (sample["tumor_content"] and not is_number(sample["tumor_content_num"])) or \
							 (sample["tumor_content"] and is_number(sample["tumor_content_num"]) and float(sample["tumor_content_num"]) < 30)) \
							else ""
	result.append(hrd_result + note)
	if gss["summary"]:
		result.append("HRR通路相关基因突变：{0}".format(gss["summary"]))
	return ("；\n".join(result)+"。")
jinja2.filters.FILTERS["hrd_summary_v3"] = hrd_summary_v3

# 8. GA检测结果
def ga_summary(info):
	ga = info[0]
	msi = info[1]
	result = []
	if ga["ebv_type"]["ebv_type"] == "P" and msi["var_id"] == "MSI-H":
		if ga["ebv_sum"]:
			result.append("EB病毒感染型（EBV），{0}".format(ga["ebv_sum"]))
		else:
			result.append("EB病毒感染型（EBV）")
		result.append("微卫星不稳定型（MSI）")
	elif ga["ebv_type"]["ebv_type"] == "P":
		result.append("EB病毒感染型（EBV）")
	elif msi["var_id"] == "MSI-H":
		result.append("微卫星不稳定型（MSI）")
	else:
		if ga["gs_sum"]:
			result.append("基因组稳定型（GS）")
		if ga["cin_sum"]:
			result.append("染色体不稳定型（CIN）")
	if not result:
		result = ["-"]
	rt = RichText()
	rt.add("\n".join(result))
	return rt
jinja2.filters.FILTERS["ga_summary"] = ga_summary

# 8. GA检测结果-2025.03.24-更新
def ga_summary_v2(info):
	ga = info[0]
	msi = info[1]
	result = []
	if ga["ebv_type"]["ebv_type"] == "P" and msi["var_id"] == "MSI-H":
		if ga["ebv_sum"]:
			result.append("EB病毒感染型（EBV），{0}".format(ga["ebv_sum"]))
		else:
			result.append("EB病毒感染型（EBV）")
		result.append("微卫星不稳定型（MSI）")
	elif ga["ebv_type"]["ebv_type"] == "P":
		result.append("EB病毒感染型（EBV）")
	elif msi["var_id"] == "MSI-H":
		result.append("微卫星不稳定型（MSI）")
	else:
		if ga["gs_sum"]:
			result.append("基因组稳定型（GS）")
		if ga["cin_sum"]:
			result.append("染色体不稳定型（CIN）")
	if not result:
		result = ["-"]
	#rt = RichText()
	#rt.add("\n".join(result))
	return "\n".join(result)
jinja2.filters.FILTERS["ga_summary_v2"] = ga_summary_v2

# 9. 子宫内膜癌分子分型结果
def ec_summary(ec_type):
	ec_dict = {
		"POLE-ultramutated type EC" : "POLE突变型（POLE mutation，POLE mut）",
		"MSI-H type EC" : "错配修复功能缺陷（Mismatch repair deficiency，MMRd）",
		"CNH type EC" : "TP53基因突变（p53 abnormality，p53 abn）",
		"CNL type EC" : "非特异性分子谱（Non-specific molecular profile，NSMP）"
	}
	return ec_dict.get(ec_type, ec_type)
jinja2.filters.FILTERS["ec_summary"] = ec_summary
	
# 10. 胚系变异（需要展示具体变异）
def var_g_summary(var_list):
	g_var_list = var_list["var_germline"]["level_5"] + var_list["var_germline"]["level_4"]
	result = []
	if g_var_list:
		for var in g_var_list:
			hgvs = var["hgvs_c"] if var["hgvs_p"] == "p.?" else var["hgvs_p"]
			clinic_g_cn = "致病性变异" if var["clinic_num_g"] == 5 else "疑似致病性变异"
			result.append("检出{0} {1}，为{2}".format(var["gene_symbol"], hgvs, clinic_g_cn))
	else:
		result = ["在检测范围内，未检出致病性/疑似致病性变异"]
	return "；".join(result)+"。"
jinja2.filters.FILTERS["var_g_summary"] = var_g_summary

# 10. 体细胞/来源不明变异
def var_s_summary(info):
	var_list = info[0]
	sample = info[1]
	result = []
	def sum_var(vlist):
		v_result = []
		for var in vlist:
			if var["bio_category"] == "Snvindel":
				if var["hgvs_p"] != "p.?":
					v_result.append(var["gene_symbol"]+" "+var["hgvs_p"])
				else:
					v_result.append(var["gene_symbol"]+" "+var["hgvs_c"])
			elif var["bio_category"] == "Cnv":
				v_result.append(var["gene_symbol"]+" 扩增")
			elif var["bio_category"] in ["Sv", "PSeqRnaSv"]:
				if var["five_prime_gene"] == "MET" and var["three_prime_gene"] == "MET":
					v_result.append("MET exon14 跳跃")
				else:
					if var["five_prime_gene"]+"-"+var["three_prime_gene"]+" 融合" not in v_result:
						v_result.append(var["five_prime_gene"]+"-"+var["three_prime_gene"]+" 融合")
			# 2025.07.01-新增PHd
			elif var["bio_category"] == "PHd":
				if var["type"] == "HomoDel":
					if var["gene_symbol"]+" 纯合缺失" not in v_result:
						v_result.append(var["gene_symbol"]+" 纯合缺失")
				elif var["type"] == "HeteDel":
					if var["gene_symbol"]+" 杂合缺失" not in v_result:
						v_result.append(var["gene_symbol"]+" 杂合缺失")
				else:
					v_result.append(var["gene_symbol"]+" 未知变异类型！")
			# 2025.07.01-新增完成
		return ", ".join(v_result)
	# 有对照样本时
	c_var_all_num = len(var_list["var_somatic"]["level_I"] + var_list["var_somatic"]["level_II"] + \
				 var_list["var_somatic"]["level_onco_nodrug"] + var_list["var_somatic"]["level_III"])
	c_var_onco_drug_num = len(var_list["var_somatic"]["level_I"] + var_list["var_somatic"]["level_II"])
	c_var_onco_nodrug_num = len(var_list["var_somatic"]["level_onco_nodrug"])
	c_var_onco_drug_str = sum_var(var_list["var_somatic"]["level_I"] + var_list["var_somatic"]["level_II"])
	# 无对照样本时
	# 2024.05.20-总变异数需要加上疑似胚系非致病但是有用药的位点
	#var_all_num = len(var_list["var_somatic"]["level_I"] + var_list["var_somatic"]["level_II"] + \
	#			 var_list["var_somatic"]["level_onco_nodrug"] + var_list["var_somatic"]["level_III"] +\
	#			 var_list["var_germline"]["level_5"] + var_list["var_germline"]["level_4"])
	var_all_num = len(var_list["var_somatic"]["level_I"] + var_list["var_somatic"]["level_II"] + \
				 var_list["var_somatic"]["level_onco_nodrug"] + var_list["var_somatic"]["level_III"] +\
				 var_list["var_germline"]["level_5"] + var_list["var_germline"]["level_4"] + var_list["var_germline"]["regimen_noonco_var"])
	var_onco_drug_num = len(var_list["var_somatic"]["level_I"] + var_list["var_somatic"]["level_II"]+\
							var_list["var_germline"]["regimen_level_I"] + var_list["var_germline"]["regimen_level_II"])
	var_onco_nodrug_num = len(var_list["var_somatic"]["level_onco_nodrug"] + var_list["var_germline_nodrug"])
	var_onco_drug_str = sum_var(var_list["var_somatic"]["level_I"] + var_list["var_somatic"]["level_II"]+\
								  var_list["var_germline"]["regimen_level_I"] + var_list["var_germline"]["regimen_level_II"])
	if sample["control_sample_id"]:
		str_1 = "检出{0}个体细胞变异，其中具有临床意义的变异有{1}个，肿瘤发生发展相关变异有{2}个。".format(str(c_var_all_num), str(c_var_onco_drug_num), str(c_var_onco_nodrug_num))
		str_2 = "具有临床意义的变异有{0}。".format(c_var_onco_drug_str)
		if c_var_all_num != 0:
			if c_var_onco_drug_num != 0:
				return str_1 + str_2
			else:
				return str_1
		else:
			return "在检测范围内，未检出体细胞变异。"
	else:
		str_1 = "检出{0}个变异，其中具有临床意义的变异有{1}个，肿瘤发生发展相关变异有{2}个。".format(str(var_all_num), str(var_onco_drug_num), str(var_onco_nodrug_num))
		str_2 = "具有临床意义的变异有{0}。".format(var_onco_drug_str)
		if var_all_num != 0:
			if var_onco_drug_num != 0:
				return str_1 + str_2
			else:
				return str_1
		else:
			return "在检测范围内，未检出变异。"
jinja2.filters.FILTERS["var_s_summary"] = var_s_summary

#----------------------------------------------------------------------------------------------------------------

### HRD ###
# 1. BRCA检测结果
def hrd_brca_result(var_list):
	result = []
	if var_list:
		for var in var_list:
			# 2025.07.01-新增PHd
			if var["bio_category"] == "PHd":
				if var["type"] == "HomoDel":
					result.append("纯合缺失 " + str(var["region"]))
				elif var["type"] == "HeteDel":
					result.append("杂合缺失 " + str(var["region"]))
				else:
					result.append("未知变异类型！" + str(var["region"]))
			elif var["bio_category"] == "Snvindel":
				if var["hgvs_p"] != "p.?":
					result.append(var["gene_region"]+" "+var["hgvs_c"]+" "+var["hgvs_p"])
				else:
					result.append(var["gene_region"]+" "+var["hgvs_c"])
	else:
		result = ["未检出致病性或疑似致病性变异"]
	rt = RichText()
	rt.add("\n".join(result))
	return rt
jinja2.filters.FILTERS["hrd_brca_result"] = hrd_brca_result

### 治疗策略-解读 ###
def get_evi_sum(evi_sum):
	result = []
	if evi_sum and "evi_split" in evi_sum.keys():
		if "Predictive_merge" in evi_sum["evi_split"].keys() and evi_sum["evi_split"]["Predictive_merge"]:
			result.extend([{"regimen" : a["regimen_name"], "inter" : a["evi_interpretation"]} for a in evi_sum["evi_split"]["Predictive_merge"]])
		if "Prognostic" in evi_sum["evi_split"].keys() and evi_sum["evi_split"]["Prognostic"]:
			result.extend([{"regimen" : "预后相关", "inter" : a["evi_interpretation"]} for a in evi_sum["evi_split"]["Prognostic"]])
		if "Diagnostic" in evi_sum["evi_split"].keys() and evi_sum["evi_split"]["Diagnostic"]:
			result.extend([{"regimen" : "辅助诊断相关", "inter" : a["evi_interpretation"]} for a in evi_sum["evi_split"]["Diagnostic"]])
	if not result:
		result = [{"regimen" : "", "inter": "目前关于该变异的临床治疗实践尚不明确。"}]
	return result
jinja2.filters.FILTERS["evi_sum"] = get_evi_sum

### 其他类 ###
# 1. 治疗方案介绍-生物标志物
def approval_regimen_biomarker(info):
	biomaker_list = info[0]
	judge_mergeMET = info[1]
	result = []
	for i in biomaker_list:
		if "hgvs_c" in i.keys() and i["hgvs_c"]:
			if i["hgvs_p"] != "p.?":
				result.append("{0} {1} {2} {3}".format(i["gene_symbol"], i["gene_region"], i["hgvs_c"], i["hgvs_p"]))
			else:
				result.append("{0} {1} {2}".format(i["gene_symbol"], i["gene_region"], i["hgvs_c"]))
		elif "cnv_type" in i.keys() and i["cnv_type"]:
			# 2024.08.30-CNV 区分Loss，其他的写扩增
			if i["cnv_type"] == "Loss" and i["gene_symbol"] in ["BRCA1", "BRCA2"]:
				result.append("{0} 大片段缺失".format(i["gene_symbol"]))
			else:
				# 2026.01.12-新增oncopro血液-loss（小写）为缺失，其他为扩增
				if i["cnv_type"] in ["loss", "Loss"]:
					result.append("{0} 缺失".format(i["gene_symbol"]))
				else:
					result.append("{0} 扩增".format(i["gene_symbol"]))
		elif "five_prime_gene" in i.keys() and i["five_prime_gene"]:
			if i["five_prime_gene"] == "MET" and i["three_prime_gene"] == "MET":
				result.append("MET exon14 跳跃")
			else:
				# 重新拆分hgvs，CP40的region少了ins
				if "hgvs" in i.keys() and i["hgvs"]:
					i["five_prime_region"] = "-".join(re.split("-", (re.split(":", i["hgvs"])[2]))[:-1]) \
						  					 if not re.search("--", i["hgvs"]) \
											 else re.split("_", (re.split("--", i["hgvs"])[0]))[-1]
					i["three_prime_region"] = re.split(":", i["hgvs"])[-1] \
											  if not re.search("--", i["hgvs"]) \
											  else re.split("_", (re.split("--", i["hgvs"])[1]))[-1]
					# 加一个兼容-2023.10.19
					# var_hgvs新格式，gene1:NM_xxx:exon1--gene2:NM_xxx:exon2, 旧的为gene1:NM_xxx_exon1--gene2:NM_xxx_exon2
					# cds会变成xxx:exon1和xxx:exon2
					i["five_prime_region"] = re.split(":", i["five_prime_region"])[-1] if re.search(":", i["five_prime_region"]) else i["five_prime_region"]
					i["three_prime_region"] = re.split(":", i["three_prime_region"])[-1] if re.search(":", i["three_prime_region"]) else i["three_prime_region"]
					# 兼容完成-2023.10.19
					# 加一个兼容-2024.01.25
					# 4. gene:转录本_exon-gene:转录本_exon 重新提取后，five_prime_cds为空，以此做为重新拆分的判定依据
					if not i["five_prime_region"]:
						i["five_prime_region"] = re.split("_", re.split("-", i["hgvs"])[0])[-1]
						i["three_prime_region"] = re.split("_", i["hgvs"])[-1]
					# 兼容完成-2024.01.25
				result.append("{0}:{1}:{2}-{3}:{4}:{5}".format(i["five_prime_gene"], i["five_prime_transcript"], i["five_prime_region"],\
															   i["three_prime_gene"], i["three_prime_transcript"], i["three_prime_region"]))
		elif "biomarker_type" in i.keys() and i["biomarker_type"]:
			if i["biomarker_type"] == "KRAS/NRAS/BRAF WT":
				result.append("KRAS/NRAS/BRAF 野生型")
			elif i["biomarker_type"] == "HRD-":
				result.append("HRD阴性")
			elif i["biomarker_type"] == "HRD+":
				result.append("HRD阳性")
			else:
				result.append(i["biomarker_type"])
		else:
			result.append("无法分辨的分子标志物！")
	# 若存在MET 14跳跃DNA/RNA共检的话，则删除RNA里的结果，仅保留DNA
	result_redup = []
	for i in result:
		if i not in result_redup:
			result_redup.append(i)
	# 2024.02.27更新-MET DNA/RNA共检时都要展示
	#if judge_mergeMET:
	#	if "MET exon14 跳跃" in result_redup:
	#		result_redup.remove("MET exon14 跳跃")
	# 2024.02.27更新完成
	if not result_redup:
		result_redup = ["-"]
	rt = RichText()
	rt.add("\n".join(result_redup))
	return rt
jinja2.filters.FILTERS["approval_regimen_biomarker"] = approval_regimen_biomarker

# 1. 治疗方案介绍-生物标志物-v2-2024.12.04-CNV新增Loss、Gain、HeteDel和HomoDel
def approval_regimen_biomarker_v2(info):
	biomaker_list = info[0]
	judge_mergeMET = info[1]
	result = []
	for i in biomaker_list:
		if "hgvs_c" in i.keys() and i["hgvs_c"]:
			if i["hgvs_p"] != "p.?":
				result.append("{0} {1} {2} {3}".format(i["gene_symbol"], i["gene_region"], i["hgvs_c"], i["hgvs_p"]))
			else:
				result.append("{0} {1} {2}".format(i["gene_symbol"], i["gene_region"], i["hgvs_c"]))
		elif "cnv_type" in i.keys() and i["cnv_type"]:
			# 2024.08.30-CNV 区分Loss，其他的写扩增
			if i["cnv_type"] == "Loss" and i["gene_symbol"] in ["BRCA1", "BRCA2"]:
				result.append("{0} 大片段缺失".format(i["gene_symbol"]))
			elif i["cnv_type"] == "Gain" and i["gene_symbol"] in ["BRCA1", "BRCA2"]:
				result.append("{0} 大片段重复".format(i["gene_symbol"]))
			elif i["cnv_type"] == "HeteDel" and i["gene_symbol"] in ["BRCA1", "BRCA2"]:
				result.append("{0} 杂合大片段缺失".format(i["gene_symbol"]))
			elif i["cnv_type"] == "HomoDel" and i["gene_symbol"] in ["BRCA1", "BRCA2"]:
				result.append("{0} 纯合大片段缺失".format(i["gene_symbol"]))
			else:
				result.append("{0} 扩增".format(i["gene_symbol"]))
		elif "five_prime_gene" in i.keys() and i["five_prime_gene"]:
			if i["five_prime_gene"] == "MET" and i["three_prime_gene"] == "MET":
				result.append("MET exon14 跳跃")
			else:
				# 重新拆分hgvs，CP40的region少了ins
				if "hgvs" in i.keys() and i["hgvs"]:
					i["five_prime_region"] = "-".join(re.split("-", (re.split(":", i["hgvs"])[2]))[:-1]) \
						  					 if not re.search("--", i["hgvs"]) \
											 else re.split("_", (re.split("--", i["hgvs"])[0]))[-1]
					i["three_prime_region"] = re.split(":", i["hgvs"])[-1] \
											  if not re.search("--", i["hgvs"]) \
											  else re.split("_", (re.split("--", i["hgvs"])[1]))[-1]
					# 加一个兼容-2023.10.19
					# var_hgvs新格式，gene1:NM_xxx:exon1--gene2:NM_xxx:exon2, 旧的为gene1:NM_xxx_exon1--gene2:NM_xxx_exon2
					# cds会变成xxx:exon1和xxx:exon2
					i["five_prime_region"] = re.split(":", i["five_prime_region"])[-1] if re.search(":", i["five_prime_region"]) else i["five_prime_region"]
					i["three_prime_region"] = re.split(":", i["three_prime_region"])[-1] if re.search(":", i["three_prime_region"]) else i["three_prime_region"]
					# 兼容完成-2023.10.19
					# 加一个兼容-2024.01.25
					# 4. gene:转录本_exon-gene:转录本_exon 重新提取后，five_prime_cds为空，以此做为重新拆分的判定依据
					if not i["five_prime_region"]:
						i["five_prime_region"] = re.split("_", re.split("-", i["hgvs"])[0])[-1]
						i["three_prime_region"] = re.split("_", i["hgvs"])[-1]
					# 兼容完成-2024.01.25
				result.append("{0}:{1}:{2}-{3}:{4}:{5}".format(i["five_prime_gene"], i["five_prime_transcript"], i["five_prime_region"],\
															   i["three_prime_gene"], i["three_prime_transcript"], i["three_prime_region"]))
		elif "biomarker_type" in i.keys() and i["biomarker_type"]:
			if i["biomarker_type"] == "KRAS/NRAS/BRAF WT":
				result.append("KRAS/NRAS/BRAF 野生型")
			elif i["biomarker_type"] == "HRD-":
				result.append("HRD阴性")
			elif i["biomarker_type"] == "HRD+":
				result.append("HRD阳性")
			# 2025.07.09-新增HD
			elif "HomoDel" in i["biomarker_type"] or "HeteDel" in i["biomarker_type"]:
				split_str = re.split(":", i["biomarker_type"])
				if len(split_str) == 4:
					hd_type = "纯合缺失" if split_str[-1] == "HomoDel" else "杂合缺失" if split_str[-1] == "HeteDel" else "未知变异类型！"
					result.append(split_str[1] + " " + hd_type + " " + split_str[2])
					#result.append(" ".join(split_str[1:4]))
				else:
					result.append(i["biomarker_type"])
			else:
				result.append(i["biomarker_type"])
		else:
			result.append("无法分辨的分子标志物！")
	# 若存在MET 14跳跃DNA/RNA共检的话，则删除RNA里的结果，仅保留DNA
	result_redup = []
	for i in result:
		if i not in result_redup:
			result_redup.append(i)
	# 2024.02.27更新-MET DNA/RNA共检时都要展示
	#if judge_mergeMET:
	#	if "MET exon14 跳跃" in result_redup:
	#		result_redup.remove("MET exon14 跳跃")
	# 2024.02.27更新完成
	if not result_redup:
		result_redup = ["-"]
	rt = RichText()
	rt.add("\n".join(result_redup))
	return rt
jinja2.filters.FILTERS["approval_regimen_biomarker_v2"] = approval_regimen_biomarker_v2

# 1. 治疗方案介绍-生物标志物-v3-2025.03.24-取消超文本
def approval_regimen_biomarker_v3(info):
	biomaker_list = info[0]
	judge_mergeMET = info[1]
	result = []
	for i in biomaker_list:
		if "hgvs_c" in i.keys() and i["hgvs_c"]:
			if i["hgvs_p"] != "p.?":
				result.append("{0} {1} {2} {3}".format(i["gene_symbol"], i["gene_region"], i["hgvs_c"], i["hgvs_p"]))
			else:
				result.append("{0} {1} {2}".format(i["gene_symbol"], i["gene_region"], i["hgvs_c"]))
		elif "cnv_type" in i.keys() and i["cnv_type"]:
			# 2024.08.30-CNV 区分Loss，其他的写扩增
			if i["cnv_type"] == "Loss" and i["gene_symbol"] in ["BRCA1", "BRCA2"]:
				result.append("{0} 大片段缺失".format(i["gene_symbol"]))
			elif i["cnv_type"] == "Gain" and i["gene_symbol"] in ["BRCA1", "BRCA2"]:
				result.append("{0} 大片段重复".format(i["gene_symbol"]))
			elif i["cnv_type"] == "HeteDel" and i["gene_symbol"] in ["BRCA1", "BRCA2"]:
				result.append("{0} 杂合大片段缺失".format(i["gene_symbol"]))
			elif i["cnv_type"] == "HomoDel" and i["gene_symbol"] in ["BRCA1", "BRCA2"]:
				result.append("{0} 纯合大片段缺失".format(i["gene_symbol"]))
			else:
				result.append("{0} 扩增".format(i["gene_symbol"]))
		elif "five_prime_gene" in i.keys() and i["five_prime_gene"]:
			if i["five_prime_gene"] == "MET" and i["three_prime_gene"] == "MET":
				result.append("MET exon14 跳跃")
			else:
				# 重新拆分hgvs，CP40的region少了ins
				if "hgvs" in i.keys() and i["hgvs"]:
					i["five_prime_region"] = "-".join(re.split("-", (re.split(":", i["hgvs"])[2]))[:-1]) \
						  					 if not re.search("--", i["hgvs"]) \
											 else re.split("_", (re.split("--", i["hgvs"])[0]))[-1]
					i["three_prime_region"] = re.split(":", i["hgvs"])[-1] \
											  if not re.search("--", i["hgvs"]) \
											  else re.split("_", (re.split("--", i["hgvs"])[1]))[-1]
					# 加一个兼容-2023.10.19
					# var_hgvs新格式，gene1:NM_xxx:exon1--gene2:NM_xxx:exon2, 旧的为gene1:NM_xxx_exon1--gene2:NM_xxx_exon2
					# cds会变成xxx:exon1和xxx:exon2
					i["five_prime_region"] = re.split(":", i["five_prime_region"])[-1] if re.search(":", i["five_prime_region"]) else i["five_prime_region"]
					i["three_prime_region"] = re.split(":", i["three_prime_region"])[-1] if re.search(":", i["three_prime_region"]) else i["three_prime_region"]
					# 兼容完成-2023.10.19
					# 加一个兼容-2024.01.25
					# 4. gene:转录本_exon-gene:转录本_exon 重新提取后，five_prime_cds为空，以此做为重新拆分的判定依据
					if not i["five_prime_region"]:
						i["five_prime_region"] = re.split("_", re.split("-", i["hgvs"])[0])[-1]
						i["three_prime_region"] = re.split("_", i["hgvs"])[-1]
					# 兼容完成-2024.01.25
				result.append("{0}:{1}:{2}-{3}:{4}:{5}".format(i["five_prime_gene"], i["five_prime_transcript"], i["five_prime_region"],\
															   i["three_prime_gene"], i["three_prime_transcript"], i["three_prime_region"]))
		elif "biomarker_type" in i.keys() and i["biomarker_type"]:
			if i["biomarker_type"] == "KRAS/NRAS/BRAF WT":
				result.append("KRAS/NRAS/BRAF 野生型")
			elif i["biomarker_type"] == "HRD-":
				result.append("HRD阴性")
			elif i["biomarker_type"] == "HRD+":
				result.append("HRD阳性")
			# 2025.07.17-新增HD
			elif "HomoDel" in i["biomarker_type"] or "HeteDel" in i["biomarker_type"]:
				split_str = re.split(":", i["biomarker_type"])
				if len(split_str) == 4:
					hd_type = "纯合缺失" if split_str[-1] == "HomoDel" else "杂合缺失" if split_str[-1] == "HeteDel" else "未知变异类型！"
					result.append(split_str[1] + " " + hd_type + " " + split_str[2])
					#result.append(" ".join(split_str[1:4]))
				else:
					result.append(i["biomarker_type"])
			else:
				result.append(i["biomarker_type"])
		else:
			result.append("无法分辨的分子标志物！")
	# 若存在MET 14跳跃DNA/RNA共检的话，则删除RNA里的结果，仅保留DNA
	result_redup = []
	for i in result:
		if i not in result_redup:
			result_redup.append(i)
	# 2024.02.27更新-MET DNA/RNA共检时都要展示
	#if judge_mergeMET:
	#	if "MET exon14 跳跃" in result_redup:
	#		result_redup.remove("MET exon14 跳跃")
	# 2024.02.27更新完成
	if not result_redup:
		result_redup = ["-"]
	#rt = RichText()
	#rt.add("\n".join(result_redup))
	return "\n".join(result_redup)
jinja2.filters.FILTERS["approval_regimen_biomarker_v3"] = approval_regimen_biomarker_v3

# 2. PARP结果汇总表，根据产品（HRR和HRD）选择对应展示列表
def choose_parp_list(info):
	hrr_list = info[0]["cdx"]["format5_forHRR_for_new_vesion"]
	hrd_list = info[0]["cdx"]["format4_forHRDC_for_new_vesion"]
	sample = info[1]
	if sample["prod_names"] in ["HRR（全血）", "HRR（组织）", "HRR（组织 全血）"]:
		return hrr_list
	elif sample["prod_names"] in ["HRD Complete（组织）"]:
		return hrd_list
jinja2.filters.FILTERS["choose_parp_list"] = choose_parp_list

# 体细胞等级转化
def somatic_class_stran(clinic_num_s):
	if clinic_num_s == 5:
		return "I类"
	elif clinic_num_s == 4:
		return "II类"
jinja2.filters.FILTERS["somatic_class_stran"] = somatic_class_stran


# 5. 检测详细结果-检测结果
def detect_result(var):
	result = []
	if var["bio_category"] == "Snvindel":
		if var["hgvs_p"] != "p.?":
			result.append(var["gene_region"]+" "+var["hgvs_c"]+" "+var["hgvs_p"])
		else:
			result.append(var["gene_region"]+" "+var["hgvs_c"])
		result.append(var["transcript_primary"])
	elif var["bio_category"] == "Cnv":
		result.append("扩增")
	elif var["bio_category"] in ["Sv", "PSeqRnaSv"]:
		if var["five_prime_gene"] == var["three_prime_gene"] and var["three_prime_gene"] == "MET":
			result.append("MET exon14 跳跃")
			result.append(var["three_prime_transcript"])
		else:
			result.append(var["five_prime_gene"]+":"+var["five_prime_cds"]+"-"+var["three_prime_gene"]+":"+var["three_prime_cds"]+"融合")
			result.append(var["five_prime_transcript"]+"/"+var["three_prime_transcript"])
	rt = RichText()
	rt.add("\n".join(result))
	return rt
jinja2.filters.FILTERS["detect_result"] = detect_result

# 6. 检测详细结果-频率-体细胞
def detect_freq_s(var):
	freq = ""
	if var["bio_category"] == "Snvindel":
		freq = var["freq_str"]
	elif var["bio_category"] == "Cnv":
		freq = var["cn_mean"]
	elif var["bio_category"] == "Sv":
		freq = var["freq_str"]
	elif var["bio_category"] == "PSeqRnaSv":
		freq = str(var["freq"])+" copies"
	return freq
jinja2.filters.FILTERS["freq_s"] = detect_freq_s

# 7. 检测详细结果-频率-胚系
# 适用上海项目，有对照样本的情况
def detect_freq_g_SH(var):
	return "未提取到数据！" if not var["freq_sc"] else "纯合" if float(var["freq_sc"]) >= 0.85 else "杂合"
jinja2.filters.FILTERS["freq_g_SH"] = detect_freq_g_SH

# io展示整理为方法
def io_stran(io_list):
	if not io_list:
		io_list = ["-"]
	rt = RichText()
	rt.add("\n".join(io_list))
	return rt
jinja2.filters.FILTERS["io_stran"] = io_stran

# io展示整理为方法-2025.03.24
def io_stran_v2(io_list):
	if not io_list:
		io_list = ["-"]
	#rt = RichText()
	#rt.add("\n".join(io_list))
	return "\n".join(io_list)
jinja2.filters.FILTERS["io_stran_v2"] = io_stran_v2

# NCCN指南推荐基因检测结果，检测结果展示-MP
def cdx_type1(info):
	cdx_list = info[0]
	sample = info[1]
	rt = RichText()
	result = []
	if cdx_list:
		for var in cdx_list:
			if var["bio_category"] in ["Sv", "PSeqRnaSv"]:
				if var["five_prime_gene"] == "MET" and var["three_prime_gene"] == "MET":
					result.append("MET exon14 跳跃，{0}".format(var["freq"]))
				else:
					if "var_info_M" in var.keys() and var["var_info_M"]:
						result.append(var["var_info_M"]+"，"+var["freq"])
					else:
						result.append(var["var_info"]+"，"+var["freq"])
			else:
				if var["gene_symbol"] in ["BRCA1", "BRCA2"]:
					if sample["control_sample_id"]:
						if var["var_origin"] == "germline":
							result.append("{0}，{1}".format(var["var_info"], var["freq_rc_str"]))
							result.append("胚系变异")
						else:
							result.append("{0}，{1}".format(var["var_info"], var["freq"]))
							result.append("体细胞变异")
					else:
						result.append("{0}，{1}".format(var["var_info"], var["freq"]))
				else:
					result.append("{0}，{1}".format(var["var_info"], var["freq"]))
	rt.add("\n".join(result))
	return rt
jinja2.filters.FILTERS["cdx_type1"] = cdx_type1

# NCCN指南推荐基因检测结果，检测结果展示-MP-2025.03.24
def cdx_type1_v2(info):
	cdx_list = info[0]
	sample = info[1]
	rt = RichText()
	result = []
	if cdx_list:
		for var in cdx_list:
			if var["bio_category"] in ["Sv", "PSeqRnaSv"]:
				if var["five_prime_gene"] == "MET" and var["three_prime_gene"] == "MET":
					result.append("MET exon14 跳跃，{0}".format(var["freq"]))
				else:
					if "var_info_M" in var.keys() and var["var_info_M"]:
						result.append(var["var_info_M"]+"，"+var["freq"])
					else:
						result.append(var["var_info"]+"，"+var["freq"])
			else:
				if var["gene_symbol"] in ["BRCA1", "BRCA2"]:
					if sample["control_sample_id"]:
						if var["var_origin"] == "germline":
							result.append("{0}，{1}".format(var["var_info"], var["freq_rc_str"]))
							result.append("胚系变异")
						else:
							result.append("{0}，{1}".format(var["var_info"], var["freq"]))
							result.append("体细胞变异")
					else:
						result.append("{0}，{1}".format(var["var_info"], var["freq"]))
				else:
					result.append("{0}，{1}".format(var["var_info"], var["freq"]))
	#rt.add("\n".join(result))
	return "\n".join(result)
jinja2.filters.FILTERS["cdx_type1_v2"] = cdx_type1_v2

# NCCN指南推荐基因检测结果，检测结果展示-116变异检测结果
def cdx_type2_var_info(var):
	if "bio_category" in var.keys():
		if var["bio_category"] == "Snvindel":
			if var["hgvs_p"] != "p.?":
				return var["hgvs_p"]
			else:
				return var["hgvs_c"]
		elif var["bio_category"] == "Cnv":
			return "扩增"
		elif var["bio_category"] in ["Sv", "PSeqRnaSv"]:
			if var["five_prime_gene"] == "MET" and var["three_prime_gene"] == "MET":
				return "MET exon14 跳跃"
			else:
				return "{0}-{1} 融合".format(var["five_prime_gene"], var["three_prime_gene"])
	else:
		return "未检测到"
jinja2.filters.FILTERS["cdx_type2_var_info"] = cdx_type2_var_info

# NCCN指南推荐基因检测结果，检测结果展示-116丰度
def cdx_type2_freq(info):
	var = info[0]
	sample = info[1]
	rt = RichText()
	result = []
	if "bio_category" in var.keys():
		if var["bio_category"] == "Snvindel":
			if var["var_origin"] == "germline":
				result.append("未提取到freq_rc！" if not var["freq_rc"] else "纯合" if float(var["freq_rc"]) >= 0.85 else "杂合")
				#rt.add("未提取到freq_rc！" if not var["freq_rc"] else "纯合" if float(var["freq_rc"]) >= 0.85 else "杂合", \
				#	   size = 17, color = "#000000", font = "Source Han Sans Normal")
			else:
				result.append(var["freq_str"])
				#rt.add(var["freq_str"], size = 17, color = "#000000", font = "Source Han Sans Normal")
			if var["gene_symbol"] in ["BRCA1", "BRCA2"] and sample["control_sample_id"]:
				result.append("胚系变异" if var["var_origin"] == "germline" else "体细胞变异")
				#rt.add("\n胚系变异" if var["var_origin"] == "germline" else "\n体细胞变异", size = 17, color = "#000000", font = "Source Han Sans Normal")
		elif var["bio_category"] == "Cnv":
			result.append(var["cn_mean"])
			#rt.add(var["cn_mean"], size = 17, color = "#000000", font = "Source Han Sans Normal")
		elif var["bio_category"] == "Sv":
			result.append(var["freq_str"])
			#rt.add(var["freq_str"], size = 17, color = "#000000", font = "Source Han Sans Normal")
	else:
		result.append("-")
		#rt.add("-", size = 17, color = "#000000", font = "Source Han Sans Normal")
	rt.add("\n".join(result))
	return rt
jinja2.filters.FILTERS["cdx_type2_freq"] = cdx_type2_freq

# NCCN指南推荐基因检测结果，检测结果展示-116临床意义
def cdx_type2_class(var):
	stran = {
		"germline" : {
			5 : "致病性",
			4 : "疑似致病性"
		},
		"somatic" : {
			5 : "I类",
			4 : "II类",
			3 : "III类"
		}
	}
	if "bio_category" in var.keys():
		if var["var_origin"] == "germline":
			return stran["germline"].get(var["clinic_num_g"])
		else:
			return stran["somatic"].get(var["clinic_num_s"])
	else:
		return "-"
jinja2.filters.FILTERS["cdx_type2_class"] = cdx_type2_class

# NCCN指南推荐基因检测结果-CP40
def cdx_type3_var_info(var):
	result = []
	if "var_info" in var.keys() and var["var_info"]:
		if "MET-MET" in var["var_info"]:
			result.append("MET exon14 skipping")
		elif var["merge_sv_list"]:
			for a in var["merge_sv_list"]:
				result.append(a+"融合")
		else:
			result.append(var["var_info"])
		if var["note"]:
			result.append("(MET exon14 skipping)")
	else:
		result = ["未检测到"]
	rt = RichText()
	rt.add("\n".join(result))
	return rt
jinja2.filters.FILTERS["cdx_type3_var_info"] = cdx_type3_var_info

# NCCN指南推荐基因检测结果-丰度
def cdx_type3_freq(var):
	result = []
	if "freq" in var.keys() and var["freq"]:
		if "融合" in var["var_info"]:
			result.append(str(var["freq"])+" copies")
		else:
			result.append(var["freq"])
		if var["note_freq"]:
			result.append("("+var["note_freq"]+")")
	else:
		result.append("-")
	rt = RichText()
	rt.add("\n".join(result))
	return rt
jinja2.filters.FILTERS["cdx_type3_freq"] = cdx_type3_freq

# NCCN指南推荐基因检测结果-临床意义
def cdx_type3_class(var):
	if "level" in var.keys() and var["level"]:
		if var["level"] == "5":
			return "I类"
		elif var["level"] == "4":
			return "II类"
		else:
			return "III类"
	else:
		return "-"
jinja2.filters.FILTERS["cdx_type3_class"] = cdx_type3_class

# 选择CRC25、TC21、GA18检测结果汇总
def choose_116_list(info):
	CRC25_list = info[0]["CRC25"]
	GA18_list = info[0]["GA18"]
	TC21_list = info[0]["TC21"]
	prod_names = info[1]
	if prod_names in ["CRC25（组织）", "CRC25（血液）"]:
		return CRC25_list
	elif prod_names in ["GA18（组织）", "GA18（血液）"]:
		return GA18_list
	elif prod_names in ["TC21（组织）", "TC21（血液）"]:
		return TC21_list
jinja2.filters.FILTERS["choose_116_list"] = choose_116_list

### 产品声明 ###
# 1. 序号
def product_state_index(info):
	var_brca = info[0]
	sample = info[1]
	# BRCA和HRD完整版的先替换成旧版本，这段代码先注释掉-2023.12.18
	#if sample["prod_names"] in ["BRCA1/BRCA2（全血）", "BRCA1/BRCA2（组织）", "BRCA1/BRCA2（组织 全血）", "HRD Complete（组织）"] and "complete" in sample["report_name"]:
	#	if sample["prod_names"] in ["BRCA1/BRCA2（全血）", "BRCA1/BRCA2（组织 全血）"]:
	#		if var_brca["snv_m"]["B1_G_L5"] or var_brca["snv_m"]["B2_G_L5"] or var_brca["snv_m"]["B1_G_L4"] or var_brca["snv_m"]["B2_G_L4"] or\
	#																	   var_brca["mlpa"]["B1_Loss"] or var_brca["mlpa"]["B2_Loss"]:
	#			return "5"
	#		else:
	#			return "4"
	#	elif sample["prod_names"] in ["BRCA1/BRCA2（组织）"]:
	#		if var_brca["snv_s"]["B1_L5"] or var_brca["snv_s"]["B2_L5"] or var_brca["snv_s"]["B1_L4"] or var_brca["snv_s"]["B2_L4"]:
	#			return "5"
	#		else:
	#			return "4"
	#	elif sample["prod_names"] in ["HRD Complete（组织）"]:
	#		return "4"
	#elif "hospital" in sample["report_name"]:
	if "hospital" in sample["report_name"]:
		if sample["prod_names"] in ["BRCA1/BRCA2（全血）", "BRCA1/BRCA2（组织 全血）"]:
			if var_brca["snv_m"]["B1_G_L5"] or var_brca["snv_m"]["B2_G_L5"] or var_brca["snv_m"]["B1_G_L4"] or var_brca["snv_m"]["B2_G_L4"] or\
																		   var_brca["mlpa"]["B1_Loss"] or var_brca["mlpa"]["B2_Loss"]:
				return "5"
			else:
				return "4"
		elif sample["prod_names"] in ["BRCA1/BRCA2（组织）"]:
			if var_brca["snv_s"]["B1_L5"] or var_brca["snv_s"]["B2_L5"] or var_brca["snv_s"]["B1_L4"] or var_brca["snv_s"]["B2_L4"]:
				return "5"
			else:
				return "4"
		else:
			return "4"
	else:
		if sample["prod_names"] in ["BRCA1/BRCA2（全血）", "BRCA1/2（扩增子）"]:
			if var_brca["snv_s"]["B1_L5"] or var_brca["snv_s"]["B2_L5"] or var_brca["snv_s"]["B1_L4"] or var_brca["snv_s"]["B2_L4"] or\
																		   var_brca["mlpa"]["B1_Loss"] or var_brca["mlpa"]["B2_Loss"]:
				if re.search("健康", sample["tumor_names_cn"]):
					return "5"
				else:
					return "7"
			else:
				if re.search("健康", sample["tumor_names_cn"]):
					return "4"
				else:
					return "6"
		elif sample["prod_names"] in ["BRCA1/BRCA2（组织）"]:
			if var_brca["snv_s"]["B1_L5"] or var_brca["snv_s"]["B2_L5"] or var_brca["snv_s"]["B1_L4"] or var_brca["snv_s"]["B2_L4"]:
				return "7"
			else:
				return "6"
		elif sample["prod_names"] in ["BRCA1/BRCA2（组织 全血）"]:
			if var_brca["snv_m"]["B1_G_L5"] or var_brca["snv_m"]["B2_G_L5"] or var_brca["snv_m"]["B1_G_L4"] or var_brca["snv_m"]["B2_G_L4"] or\
																		   var_brca["mlpa"]["B1_Loss"] or var_brca["mlpa"]["B2_Loss"]:
				return "7"
			else:
				return "6"
		elif sample["prod_names"] in ["林奇综合征"]:
			return "4"
		# 加一个HRR全血，健康人4，患者6
		elif sample["prod_names"] in ["HRR（全血）"]:
			if re.search("健康", sample["tumor_names_cn"]):
				return "4"
			else:
				return "6"
		else:
			return "6"
jinja2.filters.FILTERS["product_state_index"] = product_state_index

# 1. 序号
def product_state_index_v2(info):
	var_brca = info[0]
	sample = info[1]
	var_type = info[2]
	if "hospital" in sample["report_name"]:
		if sample["prod_names"] in ["BRCA1/BRCA2（全血）", "BRCA1/BRCA2（组织 全血）"]:
			if var_type == "mlpa":
				if var_brca["snv_m"]["B1_G_L5"] + var_brca["snv_m"]["B2_G_L5"] + \
				   var_brca["snv_m"]["B1_G_L4"] + var_brca["snv_m"]["B2_G_L4"] + \
				   var_brca["mlpa_v2"]["B1_mlpa_L5"] + var_brca["mlpa_v2"]["B2_mlpa_L5"] + \
				   var_brca["mlpa_v2"]["B1_mlpa_L4"] + var_brca["mlpa_v2"]["B2_mlpa_L4"]:
					return "5"
				else:
					return "4"
			else:
				if var_brca["snv_m"]["B1_G_L5"] + var_brca["snv_m"]["B2_G_L5"] + \
				   var_brca["snv_m"]["B1_G_L4"] + var_brca["snv_m"]["B2_G_L4"] + \
				   var_brca["gcnv_v2"]["B1_gcnv_L5"] + var_brca["gcnv_v2"]["B2_gcnv_L5"] + \
				   var_brca["gcnv_v2"]["B1_gcnv_L4"] + var_brca["gcnv_v2"]["B2_gcnv_L4"]:
					return "5"
				else:
					return "4"
				
		elif sample["prod_names"] in ["BRCA1/BRCA2（组织）"]:
			if var_brca["snv_s"]["B1_L5"] or var_brca["snv_s"]["B2_L5"] or var_brca["snv_s"]["B1_L4"] or var_brca["snv_s"]["B2_L4"]:
				return "5"
			else:
				return "4"
		else:
			return "4"
		
	else:
		if sample["prod_names"] in ["BRCA1/BRCA2（全血）", "BRCA1/2（扩增子）"]:
			if var_type == "mlpa":
				if var_brca["snv_m"]["B1_G_L5"] + var_brca["snv_m"]["B2_G_L5"] + \
				   var_brca["snv_m"]["B1_G_L4"] + var_brca["snv_m"]["B2_G_L4"] + \
				   var_brca["mlpa_v2"]["B1_mlpa_L5"] + var_brca["mlpa_v2"]["B2_mlpa_L5"] + \
				   var_brca["mlpa_v2"]["B1_mlpa_L4"] + var_brca["mlpa_v2"]["B2_mlpa_L4"]:
					if re.search("健康", sample["tumor_names_cn"]):
						return "5"
					else:
						return "7"
				else:
					if re.search("健康", sample["tumor_names_cn"]):
						return "4"
					else:
						return "6"
			else:
				if var_brca["snv_m"]["B1_G_L5"] + var_brca["snv_m"]["B2_G_L5"] + \
				   var_brca["snv_m"]["B1_G_L4"] + var_brca["snv_m"]["B2_G_L4"] + \
				   var_brca["gcnv_v2"]["B1_gcnv_L5"] + var_brca["gcnv_v2"]["B2_gcnv_L5"] + \
				   var_brca["gcnv_v2"]["B1_gcnv_L4"] + var_brca["gcnv_v2"]["B2_gcnv_L4"]:
					if re.search("健康", sample["tumor_names_cn"]):
						return "5"
					else:
						return "7"
				else:
					if re.search("健康", sample["tumor_names_cn"]):
						return "4"
					else:
						return "6"
		elif sample["prod_names"] in ["BRCA1/BRCA2（组织）"]:
			if var_brca["snv_s"]["B1_L5"] or var_brca["snv_s"]["B2_L5"] or var_brca["snv_s"]["B1_L4"] or var_brca["snv_s"]["B2_L4"]:
				return "7"
			else:
				return "6"
		elif sample["prod_names"] in ["BRCA1/BRCA2（组织 全血）"]:
			if var_type == "mlpa":
				if var_brca["snv_m"]["B1_G_L5"] + var_brca["snv_m"]["B2_G_L5"] + \
				   var_brca["snv_m"]["B1_G_L4"] + var_brca["snv_m"]["B2_G_L4"] + \
				   var_brca["mlpa_v2"]["B1_mlpa_L5"] + var_brca["mlpa_v2"]["B2_mlpa_L5"] + \
				   var_brca["mlpa_v2"]["B1_mlpa_L4"] + var_brca["mlpa_v2"]["B2_mlpa_L4"]:
					return "7"
				else:
					return "6"
			else:
				if var_brca["snv_m"]["B1_G_L5"] + var_brca["snv_m"]["B2_G_L5"] + \
				   var_brca["snv_m"]["B1_G_L4"] + var_brca["snv_m"]["B2_G_L4"] + \
				   var_brca["gcnv_v2"]["B1_gcnv_L5"] + var_brca["gcnv_v2"]["B2_gcnv_L5"] + \
				   var_brca["gcnv_v2"]["B1_gcnv_L4"] + var_brca["gcnv_v2"]["B2_gcnv_L4"]:
					return "7"
				else:
					return "6"
		elif sample["prod_names"] in ["林奇综合征"]:
			return "4"
		# 加一个HRR全血，健康人4，患者6
		elif sample["prod_names"] in ["HRR（全血）"]:
			if re.search("健康", sample["tumor_names_cn"]):
				return "4"
			else:
				return "6"
		else:
			return "6"
jinja2.filters.FILTERS["product_state_index_v2"] = product_state_index_v2

### sample ###
# 1. 样本类型
# 116血液项目，样本类型会将样本和对照的写一起
def sample_type(sample):
	result = ""
	if sample["prod_names"] in ["Pan116（血液）", "TC21（血液）", "GA18（血液）", "LC76（血液）", "CRC25（血液）", "Master Panel（血液）"] and sample["control_sample_id"]:
		if sample["sample_type"] and sample["control_sample_type"]:
			if sample["sample_type"] != sample["control_sample_type"]:
				result = sample["sample_type"]+","+sample["control_sample_type"]
			else:
				result = sample["sample_type"]
		elif sample["sample_type"] and not sample["control_sample_type"]:
			result = sample["sample_type"]
		elif not sample["sample_type"] and sample["control_sample_type"]:
			result = sample["control_sample_type"]
	else:
		result = sample["sample_type"]
	if not result:
		result = "-"
	return result
jinja2.filters.FILTERS["sample_type"] = sample_type

# 2. 样本数量
# 116血液项目，样本数量会将样本和对照的写一起
def sample_amount(sample):
	result = ""
	if sample["prod_names"] in ["Pan116（血液）", "TC21（血液）", "GA18（血液）", "LC76（血液）", "CRC25（血液）", "Master Panel（血液）"] and sample["control_sample_id"]:
		if sample["sample_type"] and sample["control_sample_type"]:
			if sample["sample_type"] != sample["control_sample_type"]:
				result = sample["sample_amount"]+","+sample["control_sample_amount"]
			else:
				result = sample["sample_amount"]
		elif sample["sample_type"] and not sample["control_sample_type"]:
			result = sample["sample_amount"]
		elif not sample["sample_type"] and sample["control_sample_type"]:
			result = sample["control_sample_amount"]
	else:
		result = sample["sample_amount"]
	if not result:
		result = "-"
	return result
jinja2.filters.FILTERS["sample_amount"] = sample_amount

# 3. 样本采集日期
# 血液项目采集日期从blood_collection_date中提取，其他从gather_data中提取
def gather_data(sample):
	result = ""
	if sample["prod_names"] in ["BRCA1/BRCA2（全血）", "HRR（全血）", "10基因（血液）", "61遗传基因", "Pan116（血液）", "林奇综合征",\
							    "TC21（血液）", "GA18（血液）", "LC76（血液）", "CRC25（血液）", "BPTM（全血）", "Master Panel（血液）", "BRCA1/2（扩增子）"]:
		result = sample["blood_collection_date"]
	else:
		result = sample["gather_data"]
	if not result:
		result = "-"
	return result
jinja2.filters.FILTERS["gather_data"] = gather_data


### 治疗方案介绍-过滤 ###
# MP组织，仅卵巢癌、乳腺癌、前列腺癌时展示HRD解读，故其他癌种在治疗方案部分，需要过滤掉HRD+和HRD-的内容
def approval_regimen_filter_hrd(info):
	regimen_list = info[0]
	sample = info[1]
	hrd_p = {"biomarker_type" : "HRD+"}
	hrd_n = {"biomarker_type" : "HRD-"}
	result = []
	if "Master" in sample["prod_names"] and not set(["乳腺癌", "卵巢癌", "前列腺癌"])&set(sample["tumor_list"]) and "实体瘤" not in sample["tumor_names_cn"]:
		for regimen in regimen_list:
			for a in regimen["var"]:
				if a == hrd_p or a == hrd_n:
					regimen["var"].remove(a)
			if regimen["var"]:
				result.append(regimen)
	else:
		result = regimen_list
	return result
jinja2.filters.FILTERS["approval_regimen_filter_hrd"] = approval_regimen_filter_hrd


### QC ###
# 1. HRR全血
def qc_gHRR(qc, lib_quality_control):
	ngs_qc = qc["dna_data_qc"] if "dna_data_qc" in qc.keys() else {}
	lib_qc = lib_quality_control["lib_dna_qc"] if "lib_dna_qc" in lib_quality_control.keys() else {}
	if float(ngs_qc["cleandata_q30_num"]) >= 0.75 and float(ngs_qc["depth_ssbc_num"]) >= 100:
		if lib_qc and "dna_qty" in lib_qc.keys() and lib_qc["dna_qty"] and "library_qty" in lib_qc.keys() and lib_qc["library_qty"]:
			if float(lib_qc["dna_qty"]) >= 20 and float(lib_qc["library_qty"]) >= 200:
				return "合格"
			else:
				return "风险"
		else:
			return "合格（质控项有缺失，请补齐数据后自行评估）"
	else:
		return "风险"
jinja2.filters.FILTERS["qc_gHRR"] = qc_gHRR

# 2. HRR组织

qc_stand = {
	"gHRR" : {},
	"tHRR" : {
		"qc" : {
			"dna_data_qc" : {
				"cleandata_q30_num" : 0.75,
				"depth_ssbc_num" : 300
			}
		},
		"lib_quality_control" : {
			"lib_dna_qc" : {
				"dna_qty" : 30,
				"library_qty" : 200
			}
		},
		"tumor_content" : 20
	}
}

def qc_tHRR(qc, lib_quality_control, sample, prod):
	stand = qc_stand.get(prod)
	# qc结果
	qc_result = []
	if "dna_data_qc" in qc.keys() and qc["dna_data_qc"]:
		for item in stand["qc"]["dna_data_qc"].keys():
			if item in qc["dna_data_qc"].keys() and float(qc["dna_data_qc"][item]) >= stand["dna_data_qc"][item]:
				qc_result.append("合格")
			else:
				qc_result.append("风险")
	# 湿实验结果
	if "lib_dna_qc" in lib_quality_control.keys() and lib_quality_control["lib_dna_qc"]:
		for item in stand["lib_quality_control"]["dna_data_qc"].keys():
			pass

# v4临检通用模板使用
# MP小结展示具体变异-体细胞变异
def var_sum_s(var_list):
	result = []
	for var in var_list:
		if var["bio_category"] == "Snvindel":
			if var["hgvs_p"] != "p.?":
				result.append(var["gene_symbol"]+" "+var["hgvs_p"])
			else:
				result.append(var["gene_symbol"]+" "+var["hgvs_c"])
		elif var["bio_category"] == "Cnv":
			# 2026.01.12-若cnv_type为Loss/loss，返回缺失
			if "cnv_type" in var.keys() and var["cnv_type"] and var["cnv_type"] in ["Loss", "loss"]:
				result.append(var["gene_symbol"]+" 缺失")
			else:
				result.append(var["gene_symbol"]+" 扩增")
		elif var["bio_category"] in ["Sv", "PSeqRnaSv"]:
			if var["five_prime_gene"] == "MET" and var["three_prime_gene"] == "MET":
				result.append("MET exon14 跳跃")
			else:
				# 融合可能会有重复（rna exon相同，断点不同的情况）
				if var["five_prime_gene"]+"-"+var["three_prime_gene"]+" 融合" not in result:
					result.append(var["five_prime_gene"]+"-"+var["three_prime_gene"]+" 融合")
		# 2025.07.15-新增PHd
		elif var["bio_category"] == "PHd":
			if var["type"] == "HomoDel":
				if var["gene_symbol"]+" 纯合缺失" not in result:
					result.append(var["gene_symbol"]+" 纯合缺失")
			elif var["type"] == "HeteDel":
				if var["gene_symbol"]+" 杂合缺失" not in result:
					result.append(var["gene_symbol"]+" 杂合缺失")
			else:
				if var["gene_symbol"]+" 未知变异类型！" not in result:
					result.append(var["gene_symbol"]+" 未知变异类型！")
		# 2025.07.15-新增完成
	return ", ".join(result)
jinja2.filters.FILTERS["var_sum_s_filter"] = var_sum_s

# MP小结展示具体变异-胚系变异
def var_sum_g(var_list):
	result = []
	for var in var_list:
		hgvs = var["hgvs_c"] if var["hgvs_p"] == "p.?" else var["hgvs_p"]
		clinic_g_cn = "致病变异" if var["clinic_num_g"] == 5 else "疑似致病变异"
		result.append("检出{0} {1}，为{2}".format(var["gene_symbol"], hgvs, clinic_g_cn))
	return "；".join(result)
jinja2.filters.FILTERS["var_sum_g_filter"] = var_sum_g


######################################## 定制模板过滤器#################################################
# 重庆西南CP40，检测结果中要列出【研究证据】参考文献
# 规则如下：
# 1. A或C3等级，如果匹配到NCCN的，展示“NCCN临床实践指南”
# 2. 有PMID的展示PMID号，没有的展示NCT号，都没有的展示斜杠，都有的展示PMID号
# 3. NCCN和（PMID or NCT）可能会同时出现，NCCN要展示，（PMID or NCT）放其中1个
def get_refer_CQXN(evi_sum):
	# 获取PMID
	def getPMID_from_inter(inter):
		pmid_list = []
		mat = re.compile(r"PMID.\s?\d+")
		for i in mat.findall(str(inter)):
			if re.search(":|: |：|： ", i):
				pmid = (re.split(":|: |：|： ", i))[1].replace(" ", "")
			else:
				pmid = (re.split("PMID", i))[1]
			pmid_list.append(pmid)
		return pmid_list
	# 获取NCT
	def getNCT_from_inter(inter):
		nct_list = []
		mat = re.compile(r"NCT\d+")
		for i in mat.findall(str(inter)):
			nct_list.append(i)
		return nct_list
	# 汇总研究证据，包含治疗、辅助诊断和预后
	evi_list = []
	for i in ["Predictive", "Prognostic", "Diagnostic"]:
		if i in evi_sum["evi_split"].keys() and evi_sum["evi_split"][i]:
			evi_list.extend(evi_sum["evi_split"][i])
	# 提取参考文献
	refer = []
	for evi in evi_list:
		# 判断是否有nccn
		if "A" in evi["evi_conclusion"] or "C3" in evi["evi_conclusion"]:
			if re.search("NCCN", evi["evi_interpretation"]):
				refer.append("NCCN临床实践指南")
		# 判断pmid和NCT
		pmid_list = getPMID_from_inter(evi["evi_interpretation"])
		nct_list = getNCT_from_inter(evi["evi_interpretation"])
		if pmid_list:
			for i in pmid_list:
				refer.append("PMID: "+str(i))
		else:
			if nct_list:
				refer.extend(nct_list)
	# 去重
	refer_redup = []
	for i in refer:
		if i not in refer_redup:
			refer_redup.append(i)
	if not refer_redup:
		refer_redup = ["/"]
	rt = RichText()
	rt.add("\n".join(refer_redup))
	return rt
jinja2.filters.FILTERS["refer_CQXN"] = get_refer_CQXN

def ZDY_HRR_germline_summary(var_list):
	var_dict = {}
	for var in var_list:
		if var["type"] == "Loss":
			var_info = "{0} {1} del".format(var["gene_symbol"], var["value"])
			var["clinic_num_g"] = 4
		else:
			if var["hgvs_p"] != "p.?":
				var_info = "{0} {1}:{2}:{3}".format(var["gene_symbol"], var["transcript_primary"], var["hgvs_c"], var["hgvs_p"])
			else:
				var_info = "{0} {1}:{2}".format(var["gene_symbol"], var["transcript_primary"], var["hgvs_c"])
		if var["clinic_num_g"] not in var_dict.keys():
			var_dict.setdefault(var["clinic_num_g"], [])
		var_dict[var["clinic_num_g"]].append(var_info)
	result = []
	if 5 in var_dict.keys():
		result.append("检出{0}个致病性变异，{1}".format(str(len(var_dict[5])), ", ".join(var_dict[5])))
	if 4 in var_dict.keys():
		result.append("检出{0}个疑似致病性变异，{1}".format(str(len(var_dict[4])), ", ".join(var_dict[4])))
	if not var_dict:
		result.append("未检出致病或疑似致病性胚系变异")
	return "；".join(result)
jinja2.filters.FILTERS["ZDY_HRR_germline_summary"] = ZDY_HRR_germline_summary

# 中山人民HRDC检测结果小结-更新-2023.07.12
def summary_ZSRM_hrd_v2(info):
	# 阳性：基于本次送检样本，检出HRD状态【阳性/阴性】、BRCA gene_region hgvs_c hgvs_p I类变异、TP53 gene_region hgvs_c hgvs_p。
	# 阴性：基于本次送检样本，检出HRD【阳性/阴性】，未检出I类、II类变异。
	var_dict = info[0]
	hrd = info[1]
	class_stran = {5 : "I类", 4 : "II类"}

	var_result = []
	for var in var_dict["var_somatic"]["level_I"] + var_dict["var_somatic"]["level_II"]:
		if var["hgvs_p"] != "p.?":
			var_result.append("{0} {1} {2} {3} {4}变异".format(var["gene_symbol"], var["gene_region"], var["hgvs_c"], var["hgvs_p"], class_stran.get(var["clinic_num_s"])))
		else:
			var_result.append("{0} {1} {2} {3}变异".format(var["gene_symbol"], var["gene_region"], var["hgvs_c"], class_stran.get(var["clinic_num_s"])))
	
	hrd_result = "HRD状态阳性" if hrd["var_id"] == "HRD+" else "HRD状态阴性"

	if var_result:
		return "基于本次送检样本，检出"+hrd_result+"、"+"、".join(var_result)+"。"
	else:
		return "基于本次送检样本，检出"+hrd_result+"，未检出I类、II类变异。"
jinja2.filters.FILTERS["summary_ZSRM_hrd_v2"] = summary_ZSRM_hrd_v2

# 重庆西南-三类变异删除解读为3的非编码区变异（剪接变异除外）
def filter_III_var_CQXN(var_list):
	result = []
	for var in var_list:
		if var["bio_category"] in ["Sv", "Cnv"]:
			result.append(var)
		elif var["bio_category"] == "Snvindel":
			if var["clinic_num_s"] == 3 and var["type"] in ["Intronic", "3'UTR", "5'UTR", "FlankingRegion3", "FlankingRegion5"]:
				pass
			else:
				result.append(var)
	return result
jinja2.filters.FILTERS["filter_III_var_CQXN"] = filter_III_var_CQXN

# 深圳二院HRR，可能获益临床试验仅展示BRCA基因-2023.07.17
def filter_clinicalTrial_SZEY(clinical_list):
	return [a for a in clinical_list if a["gene_symbol"] in ["BRCA1", "BRCA2"]]
jinja2.filters.FILTERS["filter_clinicalTrial_SZEY"] = filter_clinicalTrial_SZEY

# 深圳二院HRR，药物介绍仅展示BRCA基因相关的-2023.07.17
def filter_drug_SZEY(drug_list):
	for drug in drug_list:
		drug["var_filter"] = [a for a in drug["var"] if ("gene_symbol" in a.keys() and a["gene_symbol"] and a["gene_symbol"] in ["BRCA1", "BRCA2"]) or\
							  "biomarker_type" in a.keys() and "BRCA1" in a["biomarker_type"] or \
							  "biomarker_type" in a.keys() and "BRCA2" in a["biomarker_type"]]
	return [a for a in drug_list if a["var_filter"]]
jinja2.filters.FILTERS["filter_drug_SZEY"] = filter_drug_SZEY

# 河北省人民CP40，I类变异仅展示AB药物，II类变异仅展示CD药物
# info格式为[msi, knb, var_level_I, var_level_II, regimen_approval_list]
def filter_drug_HNRM(info):
	# 处理变异列表，返回指定列表、指定等级药物的结果，格式为[("阿法替尼", var1), ("阿法替尼", var2), ("厄洛替尼", var1)]
	def get_list(var_list, level_list):
		result = []
		for var in var_list:
			bio = []
			if "bio_category" in var.keys():
				if var["bio_category"] == "Snvindel":
					bio.append({
						"bio_category" : var["bio_category"],
						"gene_symbol" : var["gene_symbol"],
						"hgvs_c" : var["hgvs_c"],
						"hgvs_p" : var["hgvs_p"]
					})
				elif var["bio_category"] == "Cnv":
					bio.append({"bio_category" : var["bio_category"], "gene_symbol" : var["gene_symbol"]})
				elif var["bio_category"] == "Sv":
					if var["five_prime_gene"] == "MET" and var["three_prime_gene"] == "MET":
						bio.append({"bio_category" : var["bio_category"], "hgvs" : "MET exon14 skipping"})
					else:
						for i in var["merge_sv_list"]:
							bio.append({"bio_category" : var["bio_category"], "hgvs" : i+"融合"})
			elif "var_id" in var.keys():
				if var["var_id"] == "KRAS/NRAS/BRAF WT":
					bio.append({"var_id" : "KRAS/NRAS/BRAF 野生型"})
				elif var["var_id"] == "MSI-H":
					bio.append({"var_id" : "MSI-H"})
			if "evi_sum" in var.keys() and var["evi_sum"] and "Predictive" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Predictive"]:
				for evi in var["evi_sum"]["evi_split"]["Predictive"]:
					if evi["evi_conclusion_simple"] in level_list:
						result.append((evi["regimen_name"], bio))
		return result
	
	# 1. 汇总MSI-H AB、KNB AB、I类AB和II类CD药物
	regimen_list_filter = []
	regimen_list_filter.extend(get_list([info[0]], ["A", "B"]))
	regimen_list_filter.extend(get_list([info[1]], ["A", "B"]))
	regimen_list_filter.extend(get_list(info[2], ["A", "B"]))
	regimen_list_filter.extend(get_list(info[3], ["C", "D"]))

	# 2. 格式转化{"阿法替尼" : [var1, var2], "厄洛替尼" : [var1]}
	regimen_dict_filter = {}
	for i in regimen_list_filter:
		regimen_dict_filter.setdefault(i[0], [])
		for j in i[1]:
			regimen_dict_filter[i[0]].append(j)

	# 对治疗方案介绍列表进行过滤，并新增过滤后的变异列表var_hnrm
	regimen_approval_list = info[4]
	regimen_approval_filter = []
	for regimen in regimen_approval_list:
		regimen_name = regimen["regimen_cn"] if regimen["regimen_cn"] else regimen["regimen_en"]
		if regimen_name in regimen_dict_filter.keys():
			regimen["var_hnrm"] = regimen_dict_filter[regimen_name]
			regimen_approval_filter.append(regimen)

	return regimen_approval_filter
jinja2.filters.FILTERS["filter_drug_HNRM"] = filter_drug_HNRM

# 杭州区域肺癌需要判断EGFR基因检出情况
def judge_EGFR(info):
	var_list = info[0]
	sample = info[1]
	gene_list = [var["gene_symbol"] for var in var_list]
	if "EGFR" not in gene_list and "肺癌" in sample["tumor_list"] and "实体瘤" not in sample["tumor_names_cn"]:
		return "EGFRnone"
jinja2.filters.FILTERS["judge_EGFR"] = judge_EGFR

# 杭州区域肺癌需要判断EGFR T790M突变情况
def sort_var_forhz_egfr(info):
	var_list = info[0]
	sample = info[1]
	for var in var_list:
		var["egfr_sort"] = 0
	t790m = {"gene_symbol" : "EGFR", "t790m_var_info" : "未检出T790M变异", "egfr_sort" : 1}
	egfr_hgvs_p_list = [var["hgvs_p"].replace("p.","").replace("(","").replace(")","") for var in var_list if var["gene_symbol"] == "EGFR" and var["bio_category"] == "Snvindel"]
	
	for var in var_list:
		if var["bio_category"] == "Cnv" and var["gene_symbol"] == "EGFR":
			egfr_hgvs_p_list.append("扩增")

	if "T790M" not in egfr_hgvs_p_list and egfr_hgvs_p_list and "肺癌" in sample["tumor_list"] and "实体瘤" not in sample["tumor_names_cn"]:
		# 计算EGFR最后一个变异所在index
		egfr_last_index = 0
		for var in var_list:
			if var["gene_symbol"] == "EGFR":
				egfr_last_index = var_list.index(var)
		# 插入T790M检测结果
		var_list.insert(egfr_last_index+1, t790m)

	return var_list
jinja2.filters.FILTERS["sort_var_forhz_egfr"] = sort_var_forhz_egfr

# 复旦中山CP40-小结部分展示化疗位点结果，格式上基因名需要合并单元格-2023.07.31
def chemo_stran_fdzs(chemo_list):
	result = {}
	for i in chemo_list:
		result.setdefault(i["gene_symbol"], [])
		result[i["gene_symbol"]].append({
			"dbsnp" : i["dbsnp"],
			"genotype" : i["genotype"]
		})
	chemo_result = []
	for k, v in result.items():
		v_sort = sorted(v, key=lambda i:i["dbsnp"])
		chemo_result.append({
			"gene_symbol" : k,
			"info" : v_sort
		})
	return sorted(chemo_result, key=lambda i:i["gene_symbol"])
jinja2.filters.FILTERS["chemo_stran_fdzs"] = chemo_stran_fdzs

# 福建附一新增参考文献-2023.07.31
def refer_nccn_fjfy(sample):
	refer = {
		"非小细胞肺癌" : [
			"中国非小细胞肺癌患者EGFR+T790M基因突变检测专家共识（2018年版）",
			"二代测序技术在NSCLC中的临床应用中国专家共识（2020版）",
			"非小细胞肺癌分子病理检测临床实践指南（2021版）",
			"中国非小细胞肺癌RET基因融合临床检测专家共识（2021年版）",
			"非小细胞肺癌MET临床检测中国专家共识（2022年版）",
			"非小细胞肺癌细针穿刺细胞学标本基因检测专家共识（2022年版）",
			"非小细胞肺癌恶性浆膜腔积液分子病理检测中国专家共识（2022年版）",
			"非小细胞肺癌融合基因检测临床实践中国专家共识（2023年版）",
			"中国晚期非小细胞肺癌BRAF突变诊疗专家共识（2023年版）"
		],
		"结直肠癌" : [
			"结直肠癌及其他相关实体瘤微卫星不稳定性检测中国专家共识（2019年版）",
			"结直肠癌分子标志物临床检测中国专家共识（2021年版）",
			"结直肠癌分子检测高通量测序中国专家共识（2021年版）"
		],
		"卵巢癌" : [
			"卵巢上皮性癌BRCA基因检测的中国专家讨论（2017年版）",
			"基于下一代测序技术的BRCA1_2基因检测指南(2019年版)",
			"上皮性卵巢癌PARP抑制剂相关生物标志物检测的中国专家共识（2020年版）",
			"BRCA1_2数据解读中国专家共识(2021年版)"
		],
		"前列腺癌" : [
			"中国前列腺癌患者基因检测专家共识（2020年版）",
			"前列腺癌同源重组修复基因检测及变异解读专家共识（2022年版）"
		],
		"乳腺癌" : [
			"中国乳腺癌患者BRCA基因检测与临床应用专家共识（2018年版）",
			"基于下一代测序技术的BRCA1_2基因检测指南(2019年版)",
			"复发/转移性乳腺癌标志物临床应用专家共识（2019年版）",
			"晚期乳腺癌基因检测热点问题中国专家共识（2021年版）",
			"基于靶标指导乳腺癌精准治疗标志物临床应用专家共识（2022年版）"
		],
		"胃癌" : [
			"胃癌HER2检测指南（2016年版）",
			"HER2阳性晚期胃癌分子靶向治疗的中国专家共识（2016年版）",
			"胃癌高通量测序临床应用中国专家共识（2022年版）"
		]
	}
	result = []
	if "实体瘤" not in sample["tumor_names_cn"]:
		set_list = set(sample["tumor_list"]) & set([i for i in refer.keys()])
		for tumor in set_list:
			result.extend(refer[tumor])
	return result
jinja2.filters.FILTERS["refer_nccn_fjfy"] = refer_nccn_fjfy

# 删除变异中的BCL2L11变异-2023.08.02
def remove_bcl2l11(var_list):
	for var in var_list:
		if var["gene_symbol"] == "BCL2L11" and "hgvs_c" in var.keys() and var["hgvs_c"] and var["hgvs_c"] == "c.394+1479_394+4381del":
			var_list.remove(var)
	return var_list
jinja2.filters.FILTERS["remove_bcl2l11"] = remove_bcl2l11

# 西安交大一LC10/PAN116/CP40，删除KRAS、NRAS、BRAF、KNB中含有“CSCO”的证据
# 无需考虑删除证据后对变异等级的影响-2023.08.15
def filter_csco(evi_list):
	return [i for i in evi_list if "CSCO" not in i["evi_interpretation"]]
jinja2.filters.FILTERS["filter_csco"] = filter_csco

# 安徽省立BRCA-结果小结部分，需要基因斜体，并且变异之间用、间隔-2023.08.28
# 在每个变异（除了最后一个）后面加个顿号，就能在模板中使用for循环展示变异了
def ahsl_brca_sum(info):
	var_list = info[0]
	judge_intumor = info[1]
	var_g = {5 : "致病性变异", 4 : "疑似致病性变异", 3 : "意义不明确变异"}
	for var in var_list:
		if var["type"] == "Loss":
			var["clinic_num_g_stran"] = "疑似致病性变异"
		elif var["type"] == "Gain":
			var["clinic_num_g_stran"] = "意义不明确变异"
		else:
			if var["var_origin"] == "germline":
				var["clinic_num_g_stran"] = var_g.get(var["clinic_num_g"])
			else:
				if var["clinic_num_g"] in [4, 5]:
					var["clinic_num_s_stran"] = "I类变异" if judge_intumor == "intumor" else "II类变异"
				else:
					var["clinic_num_s_stran"] = "III类变异"
	if len(var_list) > 1:
		for var in var_list[0:-1]:
			if var["type"] in ["Loss", "Gain"] or var["var_origin"] == "germline":
				var["clinic_num_g_stran"] = var["clinic_num_g_stran"]+"、"
			else:
				var["clinic_num_s_stran"] = var["clinic_num_s_stran"]+"、"
	return var_list
jinja2.filters.FILTERS["ahsl_brca_sum"] = ahsl_brca_sum

# 南方医院CP40， NCCN指南推荐标志物检测结果汇总表，需要拆分为一/二/三级变异展示-2023.09.01
# 该过滤器可以返回指定等级结果
def split_cdx_nfyy(info):
	var_list = [var for var in info[0] if "var_info" in var.keys()]
	level = info[1]
	return [var for var in var_list if var["level"] == str(int(level))]
jinja2.filters.FILTERS["split_cdx_nfyy"] = split_cdx_nfyy

# 上海仁济CP40，化疗需要展示变异/野生型-2023.09.18
# wt 为野生型，只有一个为ref为杂合变异型，两个都不为ref为纯合变异型
def judge_chemo_genotype(info):
	dbsnp = info[0]
	alt = info[1]
	chemo_dict = {
		"rs3918290" : {"ref": "C", "wt": "C/C"},
		"rs55886062" : {"ref": "A", "wt": "A/A"},
		"rs67376798" : {"ref": "T", "wt": "T/T"},
		"rs75017182" : {"ref": "G", "wt": "G/G"},
		"rs10929302" : {"ref": "G", "wt": "G/G"},
		"rs4148323" : {"ref": "G", "wt": "G/G"},
		"rs8175347" : {"ref": "(TA)6", "wt": "(TA)6/(TA)6"}
		}
	split_alt = re.split("/", alt)
	if dbsnp in chemo_dict.keys():
		if split_alt[0] == chemo_dict.get(dbsnp)["ref"] and split_alt[1] == chemo_dict.get(dbsnp)["ref"]:
			return "野生型"
		else:
			if len(set(split_alt)) == 1:
				return "纯合变异型"
			else:
				return "杂合变异型"
		#elif split_alt[0] != chemo_dict.get(dbsnp)["ref"] and split_alt[1] != chemo_dict.get(dbsnp)["ref"]:
		#	return "纯合变异型"
		#else:
		#	return "杂合变异型"
	else:
		return "位点不在配置中"
jinja2.filters.FILTERS["judge_chemo_genotype"] = judge_chemo_genotype

# 孙逸仙116，子宫内膜癌时，结果小结和变异解读要把POLE、TP53变异升为I类，且变异解读第一段要加相关证据。
def syx_ec_var_class(info):
	var_list = info[0]
	class_num = info[1]
	level_I = [var for var in var_list if var["clinic_num_s"] == 5 or (var["clinic_num_s"] == 4 and var["gene_symbol"] in ["POLE", "TP53"])]
	level_II = [var for var in var_list if var["clinic_num_s"] == 4 and var["gene_symbol"] not in ["POLE", "TP53"]]
	if class_num == 1:
		return level_I
	elif class_num == 2:
		return level_II
	else:
		return []
jinja2.filters.FILTERS["syx_ec_var_class"] = syx_ec_var_class

# 中山人民HRDC检测结果小结-更新-2023.10.30
def summary_ZSRM_hrd_v3(info):
	# 阳性：基于本次送检样本，检出HRD状态【阳性/阴性】、BRCA gene_region hgvs_c hgvs_p I类变异、TP53 gene_region hgvs_c hgvs_p。
	# 阴性：基于本次送检样本，检出HRD【阳性/阴性】，未检出I类、II类变异。
	# 2023.10.30更新内容：兼容生信v1.3.0版本，hrd结果在gss字段中的情况。
	var_dict = info[0]
	hrd = info[1]
	class_stran = {5 : "I类", 4 : "II类"}
	gss = info[2]

	var_result = []
	for var in var_dict["var_somatic"]["level_I"] + var_dict["var_somatic"]["level_II"]:
		if var["hgvs_p"] != "p.?":
			var_result.append("{0} {1} {2} {3} {4}变异".format(var["gene_symbol"], var["gene_region"], var["hgvs_c"], var["hgvs_p"], class_stran.get(var["clinic_num_s"])))
		else:
			var_result.append("{0} {1} {2} {3}变异".format(var["gene_symbol"], var["gene_region"], var["hgvs_c"], class_stran.get(var["clinic_num_s"])))
	
	#hrd_result = "HRD状态阳性" if hrd["var_id"] == "HRD+" else "HRD状态阴性"
	if hrd and "var_id" in hrd.keys():
		hrd_result = "HRD状态阳性" if hrd["var_id"] == "HRD+" else "HRD状态阴性"
	elif gss and "var_id" in gss.keys():
		hrd_result = "HRD状态阳性" if gss["var_id"] == "HRD+" else "HRD状态阴性"
	else:
		hrd_result = "未获取到HRD结果！"

	if var_result:
		return "基于本次送检样本，检出"+hrd_result+"、"+"、".join(var_result)+"。"
	else:
		return "基于本次送检样本，检出"+hrd_result+"，未检出I类、II类变异。"
jinja2.filters.FILTERS["summary_ZSRM_hrd_v3"] = summary_ZSRM_hrd_v3

# 中山人民HRDC检测结果小结-更新-2024.12.05
def summary_ZSRM_hrd_v4(info):
	# 阳性：基于本次送检样本，检出HRD状态【阳性/阴性】、BRCA gene_region hgvs_c hgvs_p I类变异、TP53 gene_region hgvs_c hgvs_p。
	# 阴性：基于本次送检样本，检出HRD【阳性/阴性】，未检出I类、II类变异。
	# 2023.10.30更新内容：兼容生信v1.3.0版本，hrd结果在gss字段中的情况。
	# 2024.12.05更新内容：兼容生信v1.4.0版本，1.0使用hrd字段，其他版本使用gss字段。
	var_dict = info[0]
	hrd = info[1]
	class_stran = {5 : "I类", 4 : "II类"}
	gss = info[2]
	json_batch_name = info[3]

	var_result = []
	for var in var_dict["var_somatic"]["level_I"] + var_dict["var_somatic"]["level_II"]:
		if var["hgvs_p"] != "p.?":
			var_result.append("{0} {1} {2} {3} {4}变异".format(var["gene_symbol"], var["gene_region"], var["hgvs_c"], var["hgvs_p"], class_stran.get(var["clinic_num_s"])))
		else:
			var_result.append("{0} {1} {2} {3}变异".format(var["gene_symbol"], var["gene_region"], var["hgvs_c"], class_stran.get(var["clinic_num_s"])))
	
	#hrd_result = "HRD状态阳性" if hrd["var_id"] == "HRD+" else "HRD状态阴性"
	#if hrd and "var_id" in hrd.keys():
	#	hrd_result = "HRD状态阳性" if hrd["var_id"] == "HRD+" else "HRD状态阴性"
	#elif gss and "var_id" in gss.keys():
	#	hrd_result = "HRD状态阳性" if gss["var_id"] == "HRD+" else "HRD状态阴性"
	#else:
	#	hrd_result = "未获取到HRD结果！"
	
	# v1.1.0使用hrd，其他版本使用gss-2024.12.05
	if "v1.1.0" in json_batch_name:
		if hrd and "var_id" in hrd.keys():
			hrd_result = "HRD状态阳性" if hrd["var_id"] == "HRD+" else "HRD状态阴性"
		else:
			hrd_result = "未获取到HRD结果！"
	else:
		if gss and "var_id" in gss.keys():
			hrd_result = "HRD状态阳性" if gss["var_id"] == "HRD+" else "HRD状态阴性"
		else:
			hrd_result = "未获取到HRD结果！"
	# 2024.12.05-更新完成

	if var_result:
		return "基于本次送检样本，检出"+hrd_result+"、"+"、".join(var_result)+"。"
	else:
		return "基于本次送检样本，检出"+hrd_result+"，未检出I类、II类变异。"
jinja2.filters.FILTERS["summary_ZSRM_hrd_v4"] = summary_ZSRM_hrd_v4

# 获取变异列表中的Snvindel-2023.11.16

# 获取变异列表中的RNA SV-2023.11.16

# 北京医院治疗方案介绍-生物标志物
# 2025.11.12-新增HD
def approval_regimen_biomarker_BJYY(biomaker_list):
	result = []
	for i in biomaker_list:
		if "hgvs_c" in i.keys() and i["hgvs_c"]:
			if i["hgvs_p"] != "p.?":
				result.append("{0} {1}".format(i["gene_symbol"], i["hgvs_p"].replace("(", "").replace(")", "")))
			else:
				result.append("{0} {1}".format(i["gene_symbol"], i["hgvs_c"]))
		elif "cnv_type" in i.keys() and i["cnv_type"]:
			result.append("{0} 扩增".format(i["gene_symbol"]))
		elif "five_prime_gene" in i.keys() and i["five_prime_gene"]:
			if i["five_prime_gene"] == "MET" and i["three_prime_gene"] == "MET":
				result.append("MET exon14 skipping")
			else:
				# 重新拆分hgvs，CP40的region少了ins
				if "hgvs" in i.keys() and i["hgvs"]:
					i["five_prime_region"] = "-".join(re.split("-", (re.split(":", i["hgvs"])[2]))[:-1]) \
						  					 if not re.search("--", i["hgvs"]) \
											 else re.split("_", (re.split("--", i["hgvs"])[0]))[-1]
					i["three_prime_region"] = re.split(":", i["hgvs"])[-1] \
											  if not re.search("--", i["hgvs"]) \
											  else re.split("_", (re.split("--", i["hgvs"])[1]))[-1]
					# 加一个兼容-2023.10.19
					# var_hgvs新格式，gene1:NM_xxx:exon1--gene2:NM_xxx:exon2, 旧的为gene1:NM_xxx_exon1--gene2:NM_xxx_exon2
					# cds会变成xxx:exon1和xxx:exon2
					i["five_prime_region"] = re.split(":", i["five_prime_region"])[-1] if re.search(":", i["five_prime_region"]) else i["five_prime_region"]
					i["three_prime_region"] = re.split(":", i["three_prime_region"])[-1] if re.search(":", i["three_prime_region"]) else i["three_prime_region"]
					# 兼容完成-2023.10.19
					# 加一个兼容-2024.01.25
					# 4. gene:转录本_exon-gene:转录本_exon 重新提取后，five_prime_cds为空，以此做为重新拆分的判定依据
					if not i["five_prime_region"]:
						i["five_prime_region"] = re.split("_", re.split("-", i["hgvs"])[0])[-1]
						i["three_prime_region"] = re.split("_", i["hgvs"])[-1]
					# 兼容完成-2024.01.25
				result.append("{0}:{1}-{2}:{3} 融合".format(i["five_prime_gene"], i["five_prime_region"], i["three_prime_gene"], i["three_prime_region"]))
		elif "biomarker_type" in i.keys() and i["biomarker_type"]:
			if i["biomarker_type"] == "KRAS/NRAS/BRAF WT":
				result.append("KRAS/NRAS/BRAF WT")
			elif i["biomarker_type"] == "HRD-":
				result.append("HRD-")
			elif i["biomarker_type"] == "HRD+":
				result.append("HRD+")
			# 2025.11.12-新增HD
			elif "HomoDel" in i["biomarker_type"] or "HeteDel" in i["biomarker_type"]:
				split_str = re.split(":", i["biomarker_type"])
				if len(split_str) == 4:
					hd_type = "纯合缺失" if split_str[-1] == "HomoDel" else "杂合缺失" if split_str[-1] == "HeteDel" else "未知变异类型！"
					result.append(split_str[1] + " " + hd_type + " " + split_str[2])
					#result.append(" ".join(split_str[1:4]))
				else:
					result.append(i["biomarker_type"])
			else:
				result.append(i["biomarker_type"])
		else:
			result.append("无法分辨的分子标志物！")
	# 若存在MET 14跳跃DNA/RNA共检的话，则删除RNA里的结果，仅保留DNA
	result_redup = []
	for i in result:
		if i not in result_redup:
			result_redup.append(i)
	#print (result_redup)
	#print (judge_mergeMET)
	if not result_redup:
		result_redup = ["-"]
	rt = RichText()
	rt.add("\n".join(result_redup))
	return rt
	#return result_redup
jinja2.filters.FILTERS["approval_regimen_biomarker_BJYY"] = approval_regimen_biomarker_BJYY

# 获取变异列表中的Snvindel-2023.11.16-适用药企模板
def filter_snvindel(var_list):
	return [var for var in var_list if var["bio_category"] == "Snvindel"]
jinja2.filters.FILTERS["filter_snvindel"] = filter_snvindel

# 获取变异列表中的RNA SV-2023.11.16-适用药企模板
def filter_sv(var_list):
	return [var for var in var_list if var["bio_category"] in ["Sv", "PSeqRnaSv"]]
jinja2.filters.FILTERS["filter_sv"] = filter_sv

# 过滤出CP40基因变异结果-适用云肿单独116血液出CP40报告-2023.11.23
def filter_cp40(var_list):
	gene_list = ["ALK", "FGFR1", "FGFR2", "FGFR3", "RET", "ROS1", "ERBB2", "AKT1", "BRAF", \
				 "CTNNB1", "DDR2", "DPYD", "EGFR", "ESR1", "FGFR4", "HRAS", "IDH1", "IDH2", \
				 "KEAP1", "KIT", "KRAS", "NFE2L2", "NRAS", "PDGFRA", "PIK3CA", "POLE", "PTEN",\
				 "RB1", "STK11", "TP53", "UGT1A1", "MAP2K1", "NRG1", "NTRK1", "NTRK2", "NTRK3", \
				 "CDK4", "MYC", "NKX2-1", "MET"]
	return [var for var in var_list if set(re.split(",", var["gene_symbol"])) & set(gene_list)]
jinja2.filters.FILTERS["filter_cp40"] = filter_cp40

# 过滤出CP40基因变异结果-适用云肿单独116血液出CP40报告-2025.04.01
def filter_cp40_v2(var_list):
	gene_list = ["ALK", "FGFR1", "FGFR2", "FGFR3", "RET", "ROS1", "ERBB2", "AKT1", "BRAF", \
				 "CTNNB1", "DDR2", "DPYD", "EGFR", "ESR1", "FGFR4", "HRAS", "IDH1", "IDH2", \
				 "CDKN2A", "KIT", "AKT3", "KRAS", "NF1", "NRAS", "PDGFRA", "PIK3CA", "POLE", "PTEN",\
				 "RB1", "STK11", "TP53", "UGT1A1", "MAP2K1", "NRG1", "NTRK1", "NTRK2", "NTRK3", \
				 "CDK4", "MYC", "MET"]
	return [var for var in var_list if set(re.split(",", var["gene_symbol"])) & set(gene_list)]
jinja2.filters.FILTERS["filter_cp40_v2"] = filter_cp40_v2

# 贵州肿瘤BRCA检测结果-2023.11.30
# 阳性时：检出BRCA1基因p.Q858*突变
def gzzl_brca_sum(var_list):
	result = []
	for var in var_list:
		if var["type"] == "Loss":
			result.append("{0}基因{1}大片段缺失".format(var["gene_symbol"], var["value"])) 
		# 2024.11.22-增加一个Gain
		elif var["type"] == "Gain":
			result.append("{0}基因{1}大片段重复".format(var["gene_symbol"], var["value"])) 
		# 2024.11.22-增加完成
		else:
			if var["hgvs_p"] != "p.?":
				result.append("{0}基因{1}突变".format(var["gene_symbol"], var["hgvs_p"]))
			else:
				result.append("{0}基因{1}突变".format(var["gene_symbol"], var["hgvs_c"]))
	return "检出"+"、".join(result)
jinja2.filters.FILTERS["gzzl_brca_sum"] = gzzl_brca_sum

# HRD 完整版的过滤掉BRCA变异-2023.12.06
def filter_brca_forHRD(info):
	var_list = info[0]
	sample = info[1]
	result = []
	for var in var_list:
		if sample["prod_names"] == "HRD Complete（组织）" and "complete" in sample["report_name"]:
			if var["gene_symbol"] not in ["BRCA1", "BRCA2"] and var["bio_category"] == "Snvindel":
				result.append(var)
		else:
			result.append(var)
	return result
jinja2.filters.FILTERS["filter_brca_forHRD"] = filter_brca_forHRD

# 处理药物介绍中的药理机制-2023.12.19
# 若治疗方案中各个药物均缺少药理机制，则报告中添加提示，若>=1个药物有药理机制，则不加提示。
# 化疗仍旧不展示
def get_drug_mechanism_cn(drug_details):
	result = []
	for i in drug_details:
		if i["drug_mechanism_cn"] and i["drug_mechanism_cn"] not in result and i["drug_name"] != "化疗;Chemotherapy":
			result.append(i["drug_mechanism_cn"])
	return result
jinja2.filters.FILTERS["get_drug_mechanism_cn"] = get_drug_mechanism_cn

# 福建附一CP40，肠癌的KRAS/NRAS变异，结果小结中要加备注-2023.12.22
# 描述格式为“XX基因检测到XX，丰度为X，对XXX、XXX敏感（A级）；XX基因检测到XX，拷贝数为X”，插入在丰度后面。
def fjfy_cp40_add_crcnote(info):
	inter_sum = info[0]
	sample = info[1]
	note = "（中华胃肠外科杂志2021年3月第24卷第3期发布《结直肠癌分子标志物临床检测中国专家共识》，推荐以5%作为组织学的RAS基因检测的突变丰度截断值）"
	if "肠癌" in sample["tumor_list"] and "实体瘤" not in sample["tumor_names_cn"]:
		result = []
		for var in re.split("；", inter_sum):
			gene_symbol = re.split("基因", var)[0]
			if gene_symbol in ["KRAS", "NRAS"]:
				tmp = re.split("，", var)
				tmp.insert(2, note)
				result.append("，".join(tmp))
			else:
				result.append(var)
		return "；".join(result)
	else:
		return inter_sum
jinja2.filters.FILTERS["fjfy_cp40_add_crcnote"] = fjfy_cp40_add_crcnote

# 上海仁济-CP40，判断I/II类变异中是否存在ERBB2扩增-2024.01.02
def judge_RJ_ERBB2(var_list):
	judge_erbb2 = 0
	for var in var_list:
		if var["gene_symbol"] == "ERBB2" and var["bio_category"] == "Cnv":
			judge_erbb2 = 1
	return judge_erbb2
jinja2.filters.FILTERS["judge_RJ_ERBB2"] = judge_RJ_ERBB2

## 药企项目合并-2024.01.12
# 获取变异列表中的Snvindel-2023.11.16
def filter_snvindel(var_list):
	return [var for var in var_list if var["bio_category"] == "Snvindel"]
jinja2.filters.FILTERS["filter_snvindel"] = filter_snvindel

# 获取变异列表中的RNA SV-2023.11.16
def filter_sv(var_list):
	return [var for var in var_list if var["bio_category"] in ["Sv", "PSeqRnaSv"]]
jinja2.filters.FILTERS["filter_sv"] = filter_sv

# 获取变异列表中的EGFR -2023.11.22
def filter_egfr(var_list):
	return [var for var in var_list if var["bio_category"] == "Snvindel" and var["gene_symbol"] == "EGFR"]
jinja2.filters.FILTERS["filter_egfr"] = filter_egfr

#  XW5301-汇总基因检测结果，包含未检测到变异的基因-2023.11.27
def getSum_forXW5301(var_list, gene_rule, gene_transcript):
	result = []
	detect_gene = []
	for var in var_list:
		for gene in set(re.split(",", var["gene_symbol"])):
			detect_gene.append(gene)
			if gene in gene_rule:
				var["transcript_primary_xw5031"] = gene_transcript.get(gene)
				result.append(var)
	for gene in set(gene_rule) - set(detect_gene):
		result.append({
			"gene_symbol" : gene,
			"transcript_primary_xw5031" : gene_transcript.get(gene)
		})
	return sorted(result, key = lambda i:gene_rule.index(i["gene_symbol"]))

# XW5301-snvindel排序-2023.11.27
def sort_for5301_snvindel(var_list):
	gene_rule = ["MET", "EGFR", "ALK", "KRAS", "ROS1", "RET", "ERBB2", "BRAF", "NRAS", "PIK3CA"]
	gene_transcript = {
		"MET" : "NM_000245",
		"EGFR" : "NM_005228",
		"ALK" : "NM_004304",
		"KRAS" : "NM_033360",
		"ROS1" : "NM_002944",
		"RET" : "NM_020975",
		"ERBB2" : "NM_004448",
		"BRAF" : "NM_004333",
		"NRAS" : "NM_002524",
		"PIK3CA" : "NM_006218"
	}
	return getSum_forXW5301(var_list, gene_rule, gene_transcript)
jinja2.filters.FILTERS["sort_for5301_snvindel"] = sort_for5301_snvindel

# XW5301-sv排序-2023.11.27
def sort_for5301_sv(var_list):
	gene_rule = ["ALK", "ROS1", "RET"]
	gene_transcript = {
		"ALK" : "NM_004304",
		"ROS1" : "NM_002944",
		"RET" : "NM_020975"
	}
	return getSum_forXW5301(var_list, gene_rule, gene_transcript)
jinja2.filters.FILTERS["sort_for5301_sv"] = sort_for5301_sv

# AD4701-tBRCA 检测结果小结-2023.11.30
def ad4701_tbrca_sum(var_list):
	clinic_5_count = len([var for var in var_list if var["clinic_num_g"] == 5])
	clinic_4_count = len([var for var in var_list if var["clinic_num_g"] == 4])
	if clinic_5_count and clinic_4_count:
		return "本次实验检出{0}个致病性变异和{1}个疑似致病性变异。".format(str(clinic_5_count), str(clinic_4_count))
	elif clinic_5_count and not clinic_4_count:
		return "本次实验检出{0}个致病性变异。".format(str(clinic_5_count))
	elif not clinic_5_count and clinic_4_count:
		return "本次实验检出{0}个疑似致病性变异。".format(str(clinic_4_count))
	else:
		return "本次实验未检出致病性或疑似致病性变异。"
jinja2.filters.FILTERS["ad4701_tbrca_sum"] = ad4701_tbrca_sum
## 药企项目合并结束-2024.01.12

# 判断肿瘤细胞含量数值是否合格-2024.01.22
# 输入格式[sample, lib_quality_control, 30]|judge_tumor_content
### 更新-2024.02.01-南昌附一MP不考虑肿瘤细胞含量了，这边统一返回True
def judge_tumor_content(info):
	sample = info[0]
	lib_qc = info[1]
	threshold = info[2]
	tumor_content = sample["tumor_content"] if "tumor_content" in sample.keys() and sample["tumor_content"] else \
					lib_qc["lib_dna_qc"]["tumor_content"] if \
						"lib_dna_qc" in lib_qc.keys() and lib_qc["lib_dna_qc"] and \
						"tumor_content" in lib_qc["lib_dna_qc"].keys() and lib_qc["lib_dna_qc"]["tumor_content"] else \
					0
	tumor_content = str(tumor_content).replace("%", "")
	threshold = float(threshold) if is_number(threshold) else 0
	#if is_number(tumor_content):
	#	if float(tumor_content) >= threshold:
	#		return True
	#	else:
	#		return False
	#else:
	#	return False
	return True
jinja2.filters.FILTERS["judge_tumor_content"] = judge_tumor_content

# 南昌附一结果小结-需要基因斜体，并且变异之间用、间隔-2024.01.23
# 在每个变异（除了最后一个）后面加个顿号，就能在模板中使用for循环展示变异了
# 仅考虑snvindel和MLPA
def ncfy_var_stran(var_list):
	for var in var_list:
		if var["type"] == "Loss":
			var["ncfy_var_info"] = var["value"]+" del"
		elif var["type"] == "Gain":
			var["ncfy_var_info"] = var["value"]+" dup"
		else:
			var["ncfy_var_info"] = var["hgvs_p"] if var["hgvs_p"] != "p.?" else var["hgvs_c"]
	if len(var_list) > 1:
		for var in var_list[0:-1]:
			var["ncfy_var_info"] = var["ncfy_var_info"] + "、"
	return var_list
jinja2.filters.FILTERS["ncfy_var_stran"] = ncfy_var_stran

# 南昌附一HRD结果小结处其他基因检测结果需要过滤掉BRCA变异-2024.01.24
def filter_brca_forncfy(var_list):
	return [var for var in var_list if var["gene_symbol"] not in ["BRCA1", "BRCA2"]]
jinja2.filters.FILTERS["filter_brca_forncfy"] = filter_brca_forncfy

# 温附二BRCA-胚系-小结处输出“致病性变异、疑似致病性变异和意义不明确变异”-2024.02.04
def wzfe_summary_germline(var_list):
	result = []
	germline_dict = {
		5 : "致病性变异",
		4 : "疑似致病性变异",
		3 : "意义不明确变异"
	}
	for var in var_list:
		if germline_dict.get(var["clinic_num_g"]) not in result:
			result.append(germline_dict.get(var["clinic_num_g"]))
	if len(result) <= 2:
		return "和".join(result)
	else:
		return "致病性变异、疑似致病性变异和意义不明确变异"
jinja2.filters.FILTERS["wzfe_summary_germline"] = wzfe_summary_germline

# 温附二BRCA-体细胞-小结处输出“临床意义明确、有潜在临床意义和临床意义不明确变异”-2024.02.05
def wzfe_summary_somatic(var_list):
	result = []
	germline_dict = {
		5 : "临床意义明确",
		4 : "有潜在临床意义",
		3 : "临床意义不明确"
	}
	for var in var_list:
		if germline_dict.get(var["clinic_num_s"]) not in result:
			result.append(germline_dict.get(var["clinic_num_s"]))
	if len(result) <= 2:
		return "和".join(result)
	else:
		return "临床意义明确、有潜在临床意义和临床意义不明确"
jinja2.filters.FILTERS["wzfe_summary_somatic"] = wzfe_summary_somatic

# 上海仁济-CP40，判断I/II类变异中是否存在KRAS G12D-2024.02.21
def judge_RJ_KRAS_G12D(var_list):
	judge_kras_g12d = 0
	for var in var_list:
		if var["gene_symbol"] == "KRAS" and var["bio_category"] == "Snvindel" and var["hgvs_p"] == "p.(G12D)":
			judge_kras_g12d = 1
	return judge_kras_g12d
jinja2.filters.FILTERS["judge_RJ_KRAS_G12D"] = judge_RJ_KRAS_G12D

# 上海仁济-CP40，判断I/II类变异中是否存在KRAS G12C-2024.02.21
def judge_RJ_KRAS_G12C(var_list):
	judge_kras_g12c = 0
	for var in var_list:
		if var["gene_symbol"] == "KRAS" and var["bio_category"] == "Snvindel" and var["hgvs_p"] == "p.(G12C)":
			judge_kras_g12c = 1
	return judge_kras_g12c
jinja2.filters.FILTERS["judge_RJ_KRAS_G12C"] = judge_RJ_KRAS_G12C

# 上海仁济-CP40，判断I/II类变异中是否存在KRAS G12D/C之外的KRAS变异-2024.02.21
def judge_RJ_KRAS_other(var_list):
	judge_kras_other = 0
	for var in var_list:
		if var["gene_symbol"] == "KRAS" and var["bio_category"] == "Snvindel" and var["hgvs_p"] not in ["p.(G12C)", "p.(G12D)"]:
			judge_kras_other = 1
	return judge_kras_other
jinja2.filters.FILTERS["judge_RJ_KRAS_other"] = judge_RJ_KRAS_other

# 上海仁济-CP40，判断I/II类变异中是否存在BRAF V600E-2024.02.21
def judge_RJ_BRAF_V600E(var_list):
	judge_braf_v600e = 0
	for var in var_list:
		if var["gene_symbol"] == "BRAF" and var["bio_category"] == "Snvindel" and var["hgvs_p"] == "p.(V600E)":
			judge_braf_v600e = 1
	return judge_braf_v600e
jinja2.filters.FILTERS["judge_RJ_BRAF_V600E"] = judge_RJ_BRAF_V600E

# 上海仁济-CP40，判断I/II类变异中是否存在BRAF 除V600E以外的变异-2024.02.21
def judge_RJ_BRAF_other(var_list):
	judge_braf_other = 0
	for var in var_list:
		if var["gene_symbol"] == "BRAF" and var["bio_category"] == "Snvindel" and var["hgvs_p"] != "p.(V600E)":
			judge_braf_other = 1
	return judge_braf_other
jinja2.filters.FILTERS["judge_RJ_BRAF_other"] = judge_RJ_BRAF_other

# XW2402 CP 获取G12C结果 2023年12月27日
def filter_G12C(var_list):
	return [var for var in var_list if var['gene_symbol'] == 'KRAS' and var['hgvs_p'] == 'p.(G12C)']
jinja2.filters.FILTERS['filter_G12C'] = filter_G12C

# 获取变异列表中的CNV-2023年12月27日
def filter_cnv(var_list):
	return [var for var in var_list if var["bio_category"] == "Cnv"]
jinja2.filters.FILTERS["filter_cnv"] = filter_cnv

# XW5101 汇总结果
def getSum_forXW5101(var_list, gene_rule):
	result = []
	detect_gene = []
	for var in var_list:
		for gene in set(re.split(",", var["gene_symbol"])):
			detect_gene.append(gene)
			if gene in gene_rule:
				result.append(var)
	for gene in set(gene_rule) - set(detect_gene):
		result.append({'gene_symbol' : gene})
	return sorted(result, key = lambda i:gene_rule.index(i["gene_symbol"]))

def sort_5101_snv(var_list):
    gene_rule = ['FLT3', 'KMT2A', 'NPM1', 'NUP98', 'TP53']
    return getSum_forXW5101(var_list, gene_rule)
jinja2.filters.FILTERS['sort_5101_snv'] = sort_5101_snv

def sort_5101_sv(var_list):
    gene_rule = ['KMT2A', 'NPM1', 'NUP98']
    return getSum_forXW5101(var_list, gene_rule)
jinja2.filters.FILTERS['sort_5101_sv'] = sort_5101_sv

# 判断治疗方案是否被NMPA批准-2024.03.12
def judge_nmpa_regimen(info):
	regimen_name = info[0]
	therapeutic_regimen = info[1]
	regimen_dict = {}
	for regimen in therapeutic_regimen:
		if regimen["regimen_cn"] and regimen["regimen_cn"] not in regimen_dict:
			regimen_dict[regimen["regimen_cn"]] = regimen["approval_organization"]
		if regimen["regimen_en"] and regimen["regimen_en"] not in regimen_dict:
			regimen_dict[regimen["regimen_en"]] = regimen["approval_organization"]
	if regimen_name in regimen_dict.keys() and regimen_dict[regimen_name] and "NMPA" in regimen_dict[regimen_name]:
		return True
	else:
		return False
jinja2.filters.FILTERS["judge_nmpa_regimen"] = judge_nmpa_regimen

# 子宫内膜癌时，删除POLE、TP53变异-适用湘雅116-2024.03.26
def remove_ecvar(info):
	var_list = info[0]
	sample = info[1]
	if "子宫内膜癌" in sample["tumor_list"]:
		return [var for var in var_list if var["gene_symbol"] not in ["POLE", "TP53"]]
	else:
		return var_list
jinja2.filters.FILTERS["remove_ecvar"] = remove_ecvar

# 判断变异是否有ABC等级证据-2024.03.29
def judge_abc_evi(var):
	regimen_level = [i["evi_conclusion_simple"] for i in var["evi_sum"]["regimen_evi_sum"] if i["evidence_type"] in ["Predictive", "Prognostic", "Diagnostic"]] \
					 if var["evi_sum"] and "regimen_evi_sum" in var["evi_sum"].keys() \
					 else []
	if set(["A", "B", "C"]) & set(regimen_level):
		return True
	else:
		return False
jinja2.filters.FILTERS["judge_abc_evi"] = judge_abc_evi

# gene_region转化为中文
def gene_region_strn(gene_region):
	region_dict = {
		"exon" : "外显子",
		"intron" : "内含子",
		"3'UTR" : "3'UTR",
		"5'UTR" : "5'UTR",
		"3'FLANKING" : "非编码区",
		"5'FLANKING" : "非编码区"
		}

	region_list_en = re.split("_", gene_region)
	region_list_cn = []
	for i in region_list_en:
		if re.search("exon", i):
			region_list_cn.append(i.replace("exon", "")+"号外显子")
		elif re.search("intron", i):
			region_list_cn.append(i.replace("intron", "")+"号内含子")
		else:
			region_cn = region_dict[i] if i in region_dict.keys() else i
			region_list_cn.append(region_cn)
	return "到".join(region_list_cn)

# 上海肺科CP新增“本次检测结果”-2024.04.16
def shfk_cp_sum(var_list):
	type_stran = {
		"3'UTR" : "3'UTR区突变",
		"5'UTR" : "5'UTR区突变",
		"Intronic" : "内含子区突变",
		"FlankingRegion3" : "侧翼区突变",
		"FlankingRegion5" : "侧翼区突变"
	}
	result = []
	for var in var_list:
		var_str = ""
		if var["bio_category"] == "Snvindel":
			gene_region_cn = gene_region_strn(var["gene_region"])
			type_cn = var["type_cn"] if var["type_cn"] != "--" else type_stran.get(var["type"], var["type"])
			if "judge_mergeMET" in var.keys() and var["judge_mergeMET"]:
				if var["hgvs_p"] != "p.?":
					var_str = "MET基因14号外显子跳跃MET {0}:{1}:{2}:{3}（MET exon14 skipping），突变丰度/拷贝数{4}（{5}）".format(
						var["transcript_primary"], var["gene_region"], var["hgvs_c"], var["hgvs_p"], var["freq_str"], var["freq_2"]
						)
				else:
					var_str = "MET基因14号外显子跳跃MET {0}:{1}:{2}（MET exon14 skipping），突变丰度/拷贝数{3}（{4}）".format(
						var["transcript_primary"], var["gene_region"], var["hgvs_c"], var["freq_str"], var["freq_2"]
						)
			else:
				if var["hgvs_p"] != "p.?":
					var_str = "{0}基因{1}{2}{3}:{4}，突变丰度{5}".format(
						var["gene_symbol"], gene_region_cn, type_cn, var["hgvs_c"], var["hgvs_p"], var["freq_str"]
						)
				else:
					var_str = "{0}基因{1}{2}{3}，突变丰度{4}".format(
						var["gene_symbol"], gene_region_cn, type_cn, var["hgvs_c"], var["freq_str"]
						)
		elif var["bio_category"] == "Cnv":
			var_str = "{0}基因扩增，拷贝数{1}".format(
				var["gene_symbol"], var["cn_mean"]
				)
		elif var["bio_category"] == "Sv":
			if var["five_prime_gene"] == "MET" and var["three_prime_gene"] == "MET":
				var_str = "MET基因14号外显子跳跃（MET exon14 skipping），拷贝数{0} copies".format(var["copies"]) 
			else:
				# EML4-ALK融合需要展示具体亚型-2024.07.18
				if var["five_prime_gene"] == "EML4" and var["three_prime_gene"] == "ALK":
					var_desc_merge = shfk_cp_alksv(var["merge_sv_list"])
				else:
					var_desc_merge = var["var_desc_merge"]
				#var_str = "{0}-{1}融合，具体的融合型为{2}，拷贝数{3} copies".format(
				#	var["five_prime_gene"], var["three_prime_gene"], var["var_desc_merge"], var["copies"]
				#	)
				var_str = "{0}-{1}融合，具体的融合型为{2}，拷贝数{3} copies".format(
					var["five_prime_gene"], var["three_prime_gene"], var_desc_merge, var["copies"]
					)
		result.append(var_str)
	# 最后一个结尾加句号，其他结尾加分号
	result_ap = []
	if len(result) > 1:
		for var in result[0:-1]:
			var = var + "；"
			result_ap.append(var)
	if len(result) >= 1:
		result_ap.append(result[-1] + "。")
	return result_ap
jinja2.filters.FILTERS["shfk_cp_sum"] = shfk_cp_sum

# 上海肺科CP40 EML4-ALK融合需要展示具体的亚型-2024.07.18
def shfk_cp_alksv(merge_sv_list):
	# gene1:exon1-gene2:exon2
	alk_sv_dict = {
		("exon13", "exon20") : "v1",
		("exon20", "exon20") : "v2",
		("exon6", "exon20") : "v3a/b",
		("exon15", "exon20") : "v4'",
		("exon2", "exon20") : "v5a/b",
		("exon18", "exon20") : "v5'",
		("exon14", "exon20") : "v7",
		("exon17", "exon20") : "v8a/b"
	}
	result = []
	for var in merge_sv_list:
		five_prime_cds = re.split("-", re.split(":", var)[1])[0]
		three_prime_cds = re.split("-", re.split(":", var)[-1])[0]
		key = (five_prime_cds, three_prime_cds)
		sv_type = alk_sv_dict.get(key, "")
		if sv_type:
			result.append(var+"（{0}型）".format(sv_type))
		else:
			result.append(var)
	return "、".join(result)


# 北京医院MP，点突变，小片段插入缺失检测结果中，预测胚系的放在了体系后面，排序需要进行调整-2024.04.24
def sort_bjyy_mp_snvindel(info):
	sample = info[0]
	var_list = info[1]
	result_var = []
	if "judge_mergeMET" in var_list["special"].keys() and var_list["special"]["judge_mergeMET"]:
		for var in var_list["special"]["BJYY_Sv"]:
			if var["bio_category"] == "Snvindel" and var["gene_symbol"] == "MET":
				result_var.append(var)
	for var in var_list["special"]["BJYY_SNV"]:
		result_var.append(var)
	# 单样本的需要展示预测胚系变异
	if not sample["control_sample_id"]:
		for var in var_list["var_germline"]["level_5"] + var_list["var_germline"]["level_4"] + var_list["var_germline"]["level_3"]:
			result_var.append(var)
	
	# 重新排序下，跟临检的规则一样
	# I/II/III类
	# 相同的等级的看证据最高等级，按A/B/C/D排
	# 证据最高等级一样的，按致病性/致癌性排
	# 基因名升序
	# 频率降序（在前面处理过了，这边不再添加排序条件）
	top_level_rule = {"A" : 0, "B" : 1, "C" : 2, "D" : 3, "N" : 4}
	clinic_num_s_rule = {5 : 0, 4 : 1, 3 : 2, 2 : 3, 1 : 4}
	result_var = sorted(result_var, key=lambda i : (top_level_rule.get(i["top_level"]), clinic_num_s_rule.get(i["clinic_num_s"]), i["gene_symbol"]))
	return result_var
jinja2.filters.FILTERS["sort_bjyy_mp_snvindel"] = sort_bjyy_mp_snvindel

# 北京医院MP，融合检测结果中，RNASV排在了SV前面，排序需要进行调整-2024.04.24
def sort_bjyy_mp_sv(var_list):
	# 重新排序下，跟临检的规则一样
	# I/II/III类
	# 相同的等级的看证据最高等级，按A/B/C/D排
	# 证据最高等级一样的，按致病性/致癌性排
	# 基因名升序
	# 频率降序（在前面处理过了，这边不再添加排序条件）
	top_level_rule = {"A" : 0, "B" : 1, "C" : 2, "D" : 3, "N" : 4}
	clinic_num_s_rule = {5 : 0, 4 : 1, 3 : 2, 2 : 3, 1 : 4}
	var_list = sorted(var_list, key=lambda i : (top_level_rule.get(i["top_level"]), clinic_num_s_rule.get(i["clinic_num_s"]), i["gene_symbol"]))
	return var_list
jinja2.filters.FILTERS["sort_bjyy_mp_sv"] = sort_bjyy_mp_sv

# 福建省立LC10/116需要判定MET exon14 skipping - 2024.04.26
# （后面处理）如果返回数据中有变异分组（var_category_names），则识别“MET Exon14 Skipping”，有的判定为14跳跃
# 无变异分组则Snvindel MET变异，证据中有A的判定为14跳跃
def judge_fjsl_met14(var):
	if var["bio_category"] == "Snvindel" and var["gene_symbol"] == "MET":
		#if "var_category_names" in var.keys():
		#	if var["var_category_names"] and "MET Exon14 Skipping" in var["var_category_names"]:
		#		return True
		#	else:
		#		return False
		#else:
		evi_level_list = []
		for evi_type in ["Diagnostic","Predictive","Prognostic"]:
			if evi_type in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"][evi_type]:
				evi_level_list.extend([evi["evi_conclusion_simple"] for evi in var["evi_sum"]["evi_split"][evi_type]])
		if "A" in evi_level_list:
			return True
		else:
			return False
	else:
		return False
jinja2.filters.FILTERS["judge_fjsl_met14"] = judge_fjsl_met14

# 江苏省中医院BPTM，结果小结需要拼接-2024.04.29
def jszy_bptm_summary(info):
	region_dict = {
		"POLE" : "的外切酶结构域（Exon3-14）以及Exon19部分区域",
		"TP53" : "全部编码区域",
		"BRCA1" : "全编码区、外显子-内含子连接区、部分内含子和UTR区",
		"BRCA2" : "全编码区、外显子-内含子连接区、部分内含子和UTR区"
	}
	ec_type = info[0]
	msi =info[1]
	result = []
	for gene in ["POLE", "MSI", "TP53", "BRCA1", "BRCA2"]:
		if gene == "MSI":
			if msi["var_id"] == "MSI-H":
				result.append("MSI状态为微卫星不稳定（MSI-H）")
			else:
				result.append("MSI状态为微卫星稳定（MSS）")
		else:
			for var in ec_type[gene+"_level12"]:
				if var["hgvs_p"] != "p.?":
					result.append("在{0}基因{1}中检测到{2}基因{3}，{4} {5} {6} {7}".format(gene, region_dict.get(gene), gene, var["varInfo_XAJDY"],\
																					   var["transcript_primary"], var["gene_region"], var["hgvs_c"], var["hgvs_p"]))
				else:
					result.append("在{0}基因{1}中检测到{2}基因{3}，{4} {5} {6}".format(gene, region_dict.get(gene), gene, var["varInfo_XAJDY"],\
																					   var["transcript_primary"], var["gene_region"], var["hgvs_c"]))
	if ec_type["POLE_level12"]:
		return "患者"+"，".join(result)+"。"
	else:
		return "患者的"+"，".join(result)+"。"
jinja2.filters.FILTERS["jszy_bptm_summary"] = jszy_bptm_summary

# 吉林省肿瘤医院gBRCA，结果小结需要拼接-2024.05.07
def jlzl_gbrca_summary(var_brca):
	result = []
	level_5_count = len(var_brca["snv_s"]["B1_L5"] + var_brca["snv_s"]["B2_L5"])
	level_4_count = len(var_brca["snv_s"]["B1_L4"] + var_brca["snv_s"]["B2_L4"] + var_brca["mlpa"]["B1_Loss"] + var_brca["mlpa"]["B2_Loss"])
	level_3_count = len(var_brca["snv_s"]["B1_L3"] + var_brca["snv_s"]["B2_L3"] + var_brca["mlpa"]["B1_Gain"] + var_brca["mlpa"]["B2_Gain"])
	if level_5_count:
		result.append("{0}项致病性变异".format(str(level_5_count)))
	if level_4_count:
		result.append("{0}项疑似致病性变异".format(str(level_4_count)))
	if level_3_count:
		result.append("{0}项意义不明确变异".format(str(level_3_count)))
	return "、".join(result)
jinja2.filters.FILTERS["jlzl_gbrca_summary"] = jlzl_gbrca_summary

# 吉林省肿瘤医院gBRCA，结果小结需要拼接(CNV代替MLPA)-2024.09.13
def jlzl_gbrca_summary_cnv(var_brca):
	result = []
	level_5_count = len(var_brca["snv_s"]["B1_L5"] + var_brca["snv_s"]["B2_L5"])
	level_4_count = len(var_brca["snv_s"]["B1_L4"] + var_brca["snv_s"]["B2_L4"] + var_brca["gcnv"]["B1_Loss"] + var_brca["gcnv"]["B2_Loss"])
	level_3_count = len(var_brca["snv_s"]["B1_L3"] + var_brca["snv_s"]["B2_L3"] + var_brca["gcnv"]["B1_Gain"] + var_brca["gcnv"]["B2_Gain"])
	if level_5_count:
		result.append("{0}项致病性变异".format(str(level_5_count)))
	if level_4_count:
		result.append("{0}项疑似致病性变异".format(str(level_4_count)))
	if level_3_count:
		result.append("{0}项意义不明确变异".format(str(level_3_count)))
	return "、".join(result)
jinja2.filters.FILTERS["jlzl_gbrca_summary_cnv"] = jlzl_gbrca_summary_cnv

# BPTM plus展示变异小结-2024.05.07
# 仅包含snvindel变异
def bptm_plus_var_sum(var_list):
	result = []
	for var in var_list:
		if var["hgvs_p"] != "p.?":
			result.append(var["gene_symbol"]+" "+var["hgvs_c"]+" "+var["hgvs_p"])
		else:
			result.append(var["gene_symbol"]+" "+var["hgvs_c"])
	return ", ".join(result)
jinja2.filters.FILTERS["bptm_plus_var_sum"] = bptm_plus_var_sum

# 判断是否检出CNV变异-适用福建附一CP-2024.05.15
def judge_cnv(var_list):
	cnv_list = [var for var in var_list if var["bio_category"] == "Cnv"]
	if cnv_list:
		return True 
	else:
		return False
jinja2.filters.FILTERS["judge_cnv"] = judge_cnv

# 孙逸仙HRD需要把exon10转化为10号外显子-2024.05.15
jinja2.filters.FILTERS["stran_gene_region_syx"] = gene_region_strn

# 判断重庆西南ptBRCA首页样品总体质量评估-2024.05.23
# 血液和组织分别：质控均合格->良好；1项超出标准->合格；>=2项超出标准->风险
# 总体评估结合血液和组织的结果，if其中有风险->判定为风险，elif其中有合格->合格，else ->良好
def judge_cqxn_ptbrca_qc(info):
	lib_qc = info[0]
	ngs_qc = info[1]
	order_tumor_content = info[2]
	# gBRCA质控标准
	gbrca_stand = {
		"dna_qty" : 30,
		"library_qty" : 150,
		"cleandata_q30_num" : 0.75,
		"mapping_ratio_num" : 0.85,
		"cover_ratio_num" : 1,
		"uni20_num" : 0.9,
		"depth_mean_num" : 500,
		"depth_min_num" : 50
	}
	# tBRCA质控标准
	tbrca_stand = {
		"dna_qty" : 30,
		"library_qty" : 150,
		"cleandata_q30_num" : 0.75,
		"mapping_ratio_num" : 0.85,
		"cover_ratio_num" : 1,
		"uni20_num" : 0.9,
		"depth_mean_num" : 500,
		"depth_ssbc_num" : 300
	}

	# 分别判断gBRCA和tBRCA各项质控
	# 输出列表["pass", "fail", "deletion"]，代表通过、不通过和质控结果缺失
	def judge_qc_items(lib_name, ngs_name, stand):
		ngs_item = ["cleandata_q30_num", "mapping_ratio_num", "cover_ratio_num", "uni20_num", "depth_mean_num", "depth_min_num", "depth_ssbc_num"]
		lib_item = ["dna_qty", "library_qty"]
		result = []
		for i in ngs_item:
			if i in stand.keys():
				if is_number(ngs_qc[ngs_name][i]) and float(ngs_qc[ngs_name][i]) >= stand.get(i):
					result.append("pass")
				else:
					result.append("fail")
		for i in lib_item:
			if lib_name in lib_qc.keys() and lib_qc[lib_name] and i in lib_qc[lib_name].keys() and lib_qc[lib_name][i]:
				if is_number(lib_qc[lib_name][i]) and float(lib_qc[lib_name][i]) >= stand.get(i):
					result.append("pass")
				else:
					result.append("fail")
			else:
				result.append("deletion")
		return result

	# 分别判断gBRCA和tBRCA的综合评估，输出良好、合格和风险
	def judge_qc(qc_judge_list):
		result = ""
		qc_counter = Counter(qc_judge_list)
		if "fail" in qc_counter.keys() and qc_counter["fail"] >= 2:
			return "risk"
		elif "fail" in qc_counter.keys() and qc_counter["fail"] == 1:
			return "pass"
		else:
			return "good"

	gbrca_qc_item_result_list = judge_qc_items("control_lib_dna_qc", "control_data_qc", gbrca_stand)
	tbrca_qc_item_result_list = judge_qc_items("lib_dna_qc", "dna_data_qc", tbrca_stand)
	# tBRCA还需判断肿瘤细胞含量
	tumor_content = order_tumor_content if order_tumor_content else \
					lib_qc["lib_dna_qc"]["tumor_content"] if "lib_dna_qc" in lib_qc.keys() and \
															 lib_qc["lib_dna_qc"] and \
															 "tumor_content" in lib_qc["lib_dna_qc"].keys() and \
															 lib_qc["lib_dna_qc"]["tumor_content"] else \
					""
	if tumor_content:
		if is_number(tumor_content.replace("%", "")) and float(tumor_content.replace("%", "")) >= 20:
			tbrca_qc_item_result_list.append("pass")
		else:
			tbrca_qc_item_result_list.append("fail")

	gbrca_sum_result = judge_qc(gbrca_qc_item_result_list)
	tbrca_sum_result = judge_qc(tbrca_qc_item_result_list)
	
	ptbrca_sum_result = ""
	if "risk" in [gbrca_sum_result, tbrca_sum_result]:
		ptbrca_sum_result = "风险"
	elif "pass" in [gbrca_sum_result, tbrca_sum_result]:
		ptbrca_sum_result = "合格"
		if "deletion" in gbrca_qc_item_result_list or "deletion" in tbrca_qc_item_result_list:
			ptbrca_sum_result += "（质控项有缺失，请补齐数据后自行评估）"
	else:
		ptbrca_sum_result = "良好"
		if "deletion" in gbrca_qc_item_result_list or "deletion" in tbrca_qc_item_result_list:
			ptbrca_sum_result += "（质控项有缺失，请补齐数据后自行评估）"
	#print ("gBRCA各项比对结果：", gbrca_qc_item_result_list)
	#print ("gBRCA总体判断结果：", gbrca_sum_result)
	#print ("tBRCA各项比对结果：", tbrca_qc_item_result_list)
	#print ("tBRCA总体判断结果：", tbrca_sum_result)
	#print ("ptBRCA判定结果：",ptbrca_sum_result)
	return ptbrca_sum_result
jinja2.filters.FILTERS["judge_cqxn_ptbrca_qc"] = judge_cqxn_ptbrca_qc

# 厦门市一116附录5-1本次检测还发现以下结果，需要对列表进行过滤并且排序-2024.05.24
# 展示非10基因体系I/II类和胚系4/5类变异
# 2024.06.03-体系新增III类变异
def xmsy_116_filter_sort(info):
	# 胚系变异及排序，致病>疑似致病
	germline_var = [var for var in info if "var_info" in var.keys() and var["var_info"] and var["var_origin"] == "germline" and var["clinic_num_g"] in [4, 5]]
	germline_var = sorted(germline_var, key=lambda i:i["clinic_num_g"], reverse = True)
	# 体细胞变异及排序,I>II类，同等级比较最高等级证据A>B>C>D，基因名，变异类型snvindel>cnv>sv
	# 体细胞新增III类变异-2024.06.03
	#somatic_var = [var for var in info if "var_info" in var.keys() and var["var_info"] and var["var_origin"] != "germline" and var["level"] in ["4", "5"]]
	somatic_var = [var for var in info if "var_info" in var.keys() and var["var_info"] and var["var_origin"] != "germline" and var["level"] in ["3", "4", "5"]]
	var_type_rule = {"Snvindel" : 0, "Cnv" : 1, "Sv" : 2, "PSeqRnaSv" : 3}
	top_level_rule = {"A" : 0, "B" : 1, "C" : 2, "D" : 3, "N" : 4}
	clinic_num_s_rule = {5 : 0, 4 : 1, 3 : 2, 2 : 3, 1 : 4}
	somatic_var = sorted(somatic_var, key=lambda i:(top_level_rule.get(i["top_level"]),clinic_num_s_rule.get(i["clinic_num_s"]), i["gene_symbol"],var_type_rule.get(i["bio_category"])))

	return germline_var + somatic_var
jinja2.filters.FILTERS["xmsy_116_filter_sort"] = xmsy_116_filter_sort

# 中南大学附属湘雅医院ptHRR-其他拓展基因变异列表需要过滤掉BRCA突变-2024.05.27
def znxy_pthrr_filter_brca(var_list):
	return [var for var in var_list if var["gene_symbol"] not in ["BRCA1", "BRCA2"]]
jinja2.filters.FILTERS["znxy_pthrr_filter_brca"] = znxy_pthrr_filter_brca

# 孙逸仙HRD 附录部分展示基因介绍（除BRCA）-2024.05.29
def syx_hrd_gene_info(var_list):
	result = []
	for var in var_list:
		if var["gene_symbol"] not in ["BRCA1", "BRCA2"]:
			tmp_dict = {
				"gene_symbol" : var["gene_symbol"],
				"gene_function" : var["gene_function"]
			}
			if tmp_dict not in result:
				result.append(tmp_dict)
	return result
jinja2.filters.FILTERS["syx_hrd_gene_info"] = syx_hrd_gene_info

# 孙逸仙检HRD 测结果小结-2024.05.29
def syx_var_sum(var_list):
	type_stran = {
		"3'UTR" : "3'UTR区突变",
		"5'UTR" : "5'UTR区突变",
		"Intronic" : "内含子区突变",
		"FlankingRegion3" : "非编码区突变",
		"FlankingRegion5" : "非编码区突变"
	}
	result = []
	for var in var_list:
		hgvs = var["hgvs_p_ZJZL"] if var["hgvs_p_ZJZL"] != "p.?" else var["hgvs_c"]
		type_cn = var["type_cn"] if var["type_cn"] != "--" else type_stran.get(var["type"], var["type"])
		result.append({"gene_symbol" : var["gene_symbol"], "info" : "基因" + hgvs + type_cn})
	# 每个变异（除了最后一个）info后面加个顿号，在模板里使用for循环展示变异
	if len(var_list) > 1:
		for var in result[0:-1]:
			var["info"] = var["info"] + "、"
	return result
jinja2.filters.FILTERS["syx_var_sum"] = syx_var_sum

# 温附一BPTM Plus组织-所有基因检测结果汇总-2024.05.30
def get_tbptm_plus_summary(info):
	var_list = info[0]
	msi = info[1]
	gene_list = ["POLE", "MSI", "TP53", "BRCA1", "BRCA2", "CTNNB1", "EPCAM", "MLH1", "MSH2", "MSH6", "PMS2"]
	detect_gene_list = [var["gene_symbol"] for var in var_list]
	for gene in gene_list:
		if gene not in detect_gene_list and gene != "MSI":
			var_list.append({
				"gene_symbol" : gene
			})
	msi["gene_symbol"] = "MSI"
	var_list.append(msi)
	return sorted(var_list, key=lambda i:gene_list.index(i["gene_symbol"]))
jinja2.filters.FILTERS["get_tbptm_plus_summary"] = get_tbptm_plus_summary

# 温附一BPTM Plus血液-所有基因检测结果汇总-2024.05.30
def get_gbptm_plus_summary(var_list):
	gene_list = ["BRCA1", "BRCA2", "POLE", "TP53", "CTNNB1", "EPCAM", "MLH1", "MSH2", "MSH6", "PMS2"]
	detect_gene_list = [var["gene_symbol"] for var in var_list]
	for gene in gene_list:
		if gene not in detect_gene_list:
			var_list.append({
				"gene_symbol" : gene
			})
	return sorted(var_list, key = lambda i:gene_list.index(i["gene_symbol"]))
jinja2.filters.FILTERS["get_gbptm_plus_summary"] = get_gbptm_plus_summary

# 北大人民BPTM Plus组织-检测结果展示I/II类+肿瘤发生发展相关变异-2024.06.04
# 单个基因结果
def bdrm_bptm_plus_generesult(info):
	var_list = info[0]
	gene = info[1]
	result = [var for var in var_list if var["gene_symbol"] == gene]
	if result:
		return result
	else:
		return [{"gene_symbol" : gene}]
jinja2.filters.FILTERS["bdrm_bptm_plus_generesult"] = bdrm_bptm_plus_generesult

# 北大人民BPTM Plus组织，具有临床意义详细解读和临床意义不明变异排序-2024.06.04
def bdrm_var_sort(var_list):
	gene_list = ["BRCA1", "BRCA2", "POLE", "TP53", "CTNNB1", "EPCAM", "MLH1", "MSH2", "MSH6", "PMS2"]
	return sorted(var_list, key = lambda i:gene_list.index(i["gene_symbol"]))
jinja2.filters.FILTERS["bdrm_var_sort"] = bdrm_var_sort

# 华东医院HRR，未检测到变异的基因列表需要展示-2024.06.26
def filter_hrr_fdhd_noresult_gene(var_list):
	hrr_gene = ["AR", "ATM", "ATR", "BARD1", "BRAF", "BRCA1", "BRCA2", "BRIP1", "CDH1", "CDK12", \
				"CHEK1", "CHEK2", "ERBB2", "ESR1", "FANCA", "FANCL", "HDAC2", "HOXB13", "KRAS", \
				"MRE11", "NBN", "NRAS", "PALB2", "PIK3CA", "PPP2R2A", "PTEN", "RAD51B", "RAD51C", \
				"RAD51D", "RAD54L", "STK11", "TP53"]
	var_gene = [var["gene_symbol"] for var in var_list]
	result = []
	for gene in hrr_gene:
		if gene not in var_gene:
			result.append(gene)
	return ", ".join(sorted(result))
jinja2.filters.FILTERS["filter_hrr_fdhd_noresult_gene"] = filter_hrr_fdhd_noresult_gene

# 华东医院HRR，变异解释需要处理下-2024.06.26
def fdhd_variant_desc_cn(variant_desc_cn):
	desc_list = re.split("，", variant_desc_cn)
	if len(desc_list) <= 1:
		return ""
	else:
		return "，".join(desc_list[1:])
jinja2.filters.FILTERS["fdhd_variant_desc_cn"] = fdhd_variant_desc_cn

# 烟台毓璜顶BPTM Plus，展示每个基因的检出结果-2024.07.01
def ytyhd_gene_detect_bptmplus(info):
	gene = info[3]
	var_list_I = [var for var in info[0] if var["gene_symbol"] == gene]
	var_list_II = [var for var in info[1] if var["gene_symbol"] == gene]
	var_list_onconodrug = [var for var in info[2] if var["gene_symbol"] == gene]
	result = []
	def get_class_result(level, v_list):
		tmp = []
		tmp.append("检出{0}个{1}变异".format(str(len(v_list)), level))
		for var in v_list:
			if var["hgvs_p"] != "p.?":
				tmp.append(var["gene_region"]+" "+var["hgvs_c"]+" "+var["hgvs_p"])
			else:
				tmp.append(var["gene_region"]+" "+var["hgvs_c"])
		return "，".join(tmp)
	if var_list_I:
		result.append(get_class_result("I类", var_list_I))
	if var_list_II:
		result.append(get_class_result("II类", var_list_II))
	if var_list_onconodrug:
		result.append(get_class_result("可能与肿瘤发生发展相关", var_list_onconodrug))
	return result
jinja2.filters.FILTERS["ytyhd_gene_detect_bptmplus"] = ytyhd_gene_detect_bptmplus


# 江苏省中医BPTM plus调整变异顺序，按基因列表顺序排-2024.07.02
def jszy_bptm_plus_sort(var_list):
	gene_rule = ["POLE", "TP53", "BRCA1", "BRCA2", "CTNNB1", "EPCAM", "MLH1", "MSH2", "MSH6", "PMS2"]
	return sorted(var_list, key = lambda i:gene_rule.index(i["gene_symbol"]))
jinja2.filters.FILTERS["jszy_bptm_plus_sort"] = jszy_bptm_plus_sort


# 浙二CP小结-2024.07.03
def zjfe_cp40_var_sum(var_list):
	result = []
	for var in var_list:
		if var["bio_category"] == "Snvindel":
			if var["hgvs_p"] != "p.?":
				result.append("{0} {1}（{2}）".format(var["gene_symbol"], var["hgvs_p"], var["freq_str"]))
			else:
				result.append("{0} {1}（{2}）".format(var["gene_symbol"], var["hgvs_c"], var["freq_str"]))
		elif var["bio_category"] == "Cnv":
			result.append("{0} {1}（{2}）".format(var["gene_symbol"], "扩增", var["cn_mean"]))
		elif var["bio_category"] in ["Sv", "PSeqRnaSv"]:
			if var["five_prime_gene"] == "MET" and var["three_prime_gene"] == "MET":
				result.append("MET exon14 跳跃（{0} copies）".format(var["copies"]))
			else:
				result.append("{0}-{1} 融合（{2} copies）".format(var["five_prime_gene"], var["three_prime_gene"], var["copies"]))
	return ", ".join(result)
jinja2.filters.FILTERS["zjfe_cp40_var_sum"] = zjfe_cp40_var_sum

# 福建附一HRR和CP40检测结果第一段数据质控总结-2024.07.03
# 2025.01.26-新增CP200
def judge_fjfy_summary(info):
	prod_names = info[2]
	# 文库质量-文库浓度和文库总量
	library_concn = info[0]["library_concn"] if info[0] and "library_concn" in info[0].keys() and info[0]["library_concn"] else ""
	library_qty = info[0]["library_qty"] if info[0] and "library_qty" in info[0].keys() and info[0]["library_qty"] else ""
	if prod_names == "Classic Panel":
		lib_judge = "未填写" if not library_concn else "合格" if float(library_concn) >= 10 else "警戒"
	# oncopro组织改名Classic Panel 200（组织）-2025.05.22
	elif prod_names == "OncoPro（组织）" or prod_names == "Classic Panel 200（组织）":
		lib_judge = "未填写" if not library_concn else "合格" if float(library_concn) >= 5 else "警戒"
	else:
		lib_judge = "未填写" if not library_qty else "合格" if float(library_qty) >= 200 else "警戒"
	# 测序质量
	if prod_names == "Classic Panel":
		ngs_judge = "合格" if float(info[1]["cleandata_q30_num"]) >= 0.75 and float(info[1]["depth_ssbc_num"]) >= 400 and float(info[1]["depth_rna_ctrl_num"]) >= 20 else \
					"警戒"
	# oncopro组织改名Classic Panel 200（组织）-2025.05.22
	elif prod_names == "OncoPro（组织）" or prod_names == "Classic Panel 200（组织）":
		ngs_judge = "合格" if float(info[1]["cleandata_q30_num"]) >= 0.75 and float(info[1]["depth_ssbc_num"]) >= 400 and float(info[1]["depth_rna_ctrl_num"]) >= 10 else \
					"警戒"
	elif prod_names == "HRR（全血）":
		ngs_judge = "合格" if float(info[1]["cleandata_q30_num"]) >= 0.75 and float(info[1]["depth_ssbc_num"]) >= 100 else "警戒"
	else:
		ngs_judge = "合格" if float(info[1]["cleandata_q30_num"]) >= 0.75 and float(info[1]["depth_ssbc_num"]) >= 300 else "警戒"
	result = [
		{"type" : "文库质量", "info" : lib_judge},
		{"type" : "测序质量", "info" : ngs_judge}
	]
	# 排序，1 合格>未填写>警戒，2 文库质量>测序质量
	type_rule = ["文库质量", "测序质量"]
	info_rule = ["合格", "未填写", "警戒"]
	result = sorted(result, key = lambda i:(info_rule.index(i["info"]), type_rule.index(i["type"])))
	true_item = [i["type"] for i in result if i["info"] == "合格"]
	false_item = [i["type"] for i in result if i["info"] == "警戒"]
	bank_item = [i["type"] for i in result if i["info"] == "未填写"]
	# 文字组装
	tmp = []
	if true_item:
		if len(true_item) == 1:
			tmp.append(true_item[0]+"合格")
		else:
			tmp.append("、".join(true_item)+"均合格")
	if bank_item:
		tmp.append("、".join(bank_item)+"（数据未填写，无法判断）")
	if false_item:
		if len(false_item) == 1:
			tmp.append(false_item[0]+"警戒")
		else:
			tmp.append("、".join(false_item)+"均为警戒")
	return "，".join(tmp)
jinja2.filters.FILTERS["judge_fjfy_summary"] = judge_fjfy_summary

# 甘肃武威CP40，需要分成10基因和其他基因，展示I/II/III类变异-2024.07.05
def gsww_cp40_var(info):
	var_raw_list = info[0]
	gene_type = info[1]
	#---------------------------------------------------------------------------------------------------------
	def get_var(var_raw_list, gene_list):
		var_list = [var for var in var_raw_list if set(re.split(",", var["gene_symbol"])) & set(gene_list)]
		result = []
		detect_gene = []
		for var in var_list:
			result.append(var)
			for gene in re.split(",", var["gene_symbol"]):
				detect_gene.append(gene)
		for gene in sorted(set(gene_list) - set(detect_gene)):
			result.append({"gene_symbol" : gene})
		return result
	#----------------------------------------------------------------------------------------------------------
	gene10_list = ["ALK", "BRAF", "EGFR", "ERBB2", "KRAS", "MET", "NRAS", "PIK3CA", "RET", "ROS1"]	
	gene30_list = ["AKT1", "CDK4", "CTNNB1", "DDR2", "DPYD", "ESR1", "FGFR1", "FGFR2", "FGFR3", "FGFR4", \
				   "HRAS", "IDH1", "IDH2", "KEAP1", "KIT", "MAP2K1", "MYC", "NFE2L2", "NKX2-1", "NRG1", \
				   "NTRK1", "NTRK2", "NTRK3", "PDGFRA", "POLE", "PTEN", "RB1", "STK11", "TP53", "UGT1A1"]
	result_gene10 = get_var(var_raw_list, gene10_list)
	result_gene30 = get_var(var_raw_list, gene30_list)
	
	if gene_type == "10gene":
		return result_gene10
	elif gene_type == "30gene":
		return result_gene30
		#return sorted(result_gene30, key=lambda i:i["gene_symbol"])
jinja2.filters.FILTERS["gsww_cp40_var"] = gsww_cp40_var

# 甘肃武威HRR-BRCA和其他基因检测结果-2024.07.08
def gsww_hrr_var(info):
	var_raw_list = info[0]
	gene_type = info[1]
	# brca基因结果
	result_brca = [var for var in var_raw_list if var["gene_symbol"] in ["BRCA1", "BRCA2"]]
	for gene in set(["BRCA1", "BRCA2"]) - set([var["gene_symbol"] for var in result_brca]):
		result_brca.append({"gene_symbol" : gene})
	result_brca = sorted(result_brca, key = lambda i:i["gene_symbol"])
	# 其他基因结果
	other_gene_list = ["AR", "ATM", "ATR", "BARD1", "BRAF", "BRIP1", "CDH1", "CDK12", "CHEK1", "CHEK2", \
					   "ERBB2", "ESR1", "FANCA", "FANCL", "HDAC2", "HOXB13", "KRAS", "MRE11", "NBN", "NRAS", \
					   "PALB2", "PIK3CA", "PPP2R2A", "PTEN", "RAD51B", "RAD51C", "RAD51D", "RAD54L", "STK11", "TP53"]
	result_other = [var for var in var_raw_list if var["gene_symbol"] not in ["BRCA1", "BRCA2"]]
	for gene in set(other_gene_list) - set([var["gene_symbol"] for var in result_other]):
		result_other.append({"gene_symbol" : gene})
	result_other = sorted(result_other, key = lambda i:i["gene_symbol"])
	# 根据需求返回对应的结果
	if gene_type == "brca":
		return result_brca
	elif gene_type == "other":
		return result_other
jinja2.filters.FILTERS["gsww_hrr_var"] = gsww_hrr_var

# 甘肃武威gHRR小结，需要展示致病/疑似致病变异-2024.07.09
def gsww_ghrr_var_sum_g(var_list):
	result = []
	for var in var_list:
		if var["bio_category"] == "PMLPA":
			if var["type"] == "Loss":
				result.append(var["gene_symbol"]+" "+var["value"]+" del")
			elif var["type"] == "Gain":
				result.append(var["gene_symbol"]+" "+var["value"]+" dup")
		elif var["bio_category"] == "Snvindel":
			if var["hgvs_p"] != "p.?":
				result.append(var["gene_symbol"]+" "+var["hgvs_p"])
			else:
				result.append(var["gene_symbol"]+" "+var["hgvs_c"])
	return ", ".join(result)
jinja2.filters.FILTERS["gsww_ghrr_var_sum_g"] = gsww_ghrr_var_sum_g

# 相同证据描述的合并治疗方案展示
def merge_Predictive_evi(datainfo):
	tmp_dict = {}
	for evi in datainfo:
		tmp_dict.setdefault(evi["evi_interpretation"], [])
		tmp_dict[evi["evi_interpretation"]].append(
			{
				"regimen_name" : evi["regimen_name"],
				"evi_conclusion_simple" : evi["evi_conclusion_simple"],
				"clinical_significance_cn" : evi["clinical_significance_cn"],
				"regimen_name_py" : evi["regimen_name_py"],
				"evi_conclusion" : evi["evi_conclusion"]
			}
		)

	merge_result = []
	for k, v  in tmp_dict.items():
		merge_result.append(
			{
				"regimen_name" : "、".join([i["regimen_name"] for i in v]),
				"evi_conclusion_simple" : "/".join([i["evi_conclusion_simple"] for i in v]),
				"clinical_significance_cn" : "/".join([i["clinical_significance_cn"] for i in v]),
				"regimen_name_py" : "/".join([i["regimen_name_py"] for i in v]),
				"evi_interpretation" : k,
				"evi_conclusion" : "/".join([i["evi_conclusion"] for i in v])
			}
		)

	return merge_result

# 湘雅二医院-HRD+BRCA结果解读，治疗策略需要加上BRCA证据
# 规则：HRD治疗方案+BRCA治疗方案，若治疗方案有重复，展示HRD的即可
def xyey_evi_sum(info):
	hrd_evi_sum = info[0] if info[0] else {}
	brca_var = info[1]
	result_evi = []
	regimen_name = []
	if "evi_split"  in hrd_evi_sum.keys() and hrd_evi_sum["evi_split"] and \
					"Predictive" in hrd_evi_sum["evi_split"].keys() and hrd_evi_sum["evi_split"]["Predictive"]:
		for evi in hrd_evi_sum["evi_split"]["Predictive"]:
			if evi["regimen_name"] not in regimen_name:
				regimen_name.append(evi["regimen_name"])
				result_evi.append(evi)
	for var in brca_var:
		if "evi_split" in var["evi_sum"].keys() and var["evi_sum"]["evi_split"] and \
					"Predictive" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Predictive"]:
			for evi in var["evi_sum"]["evi_split"]["Predictive"]:
				if evi["regimen_name"] not in regimen_name:
					regimen_name.append(evi["regimen_name"])
					result_evi.append(evi)
	# 合并相同描述的证据
	merge_evi_sum = merge_Predictive_evi(result_evi)
	return merge_evi_sum
jinja2.filters.FILTERS["xyey_evi_sum"] = xyey_evi_sum

# 获取分析包版本，如ADXHS-tBPTMplus_v1.0.1 -2024.07.12
def get_analysis_version(json_batch_name):
	result = ""
	if json_batch_name:
		json_batch_list = re.split("_", json_batch_name)
		# zip上传包：日期_测序仪_测序芯片_生信分析包_版本_系统日期流水号
		if len(json_batch_list) == 6:
			result = json_batch_list[3]+"_"+json_batch_list[4]
		# export上传包：日期_测序仪_测序芯片_生信分析包_版本.export系统日期流水号
		# 2025.09.24-v4上zip上传包格式和export一样
		elif len(json_batch_list) == 5:
			version = re.split("\.export|\.zip", json_batch_list[-1])
			result = json_batch_list[3]+"_"+version[0]
	return result
jinja2.filters.FILTERS["get_analysis_version"] = get_analysis_version

# 重庆肿瘤116，变异解读后半段的描述需要做合并-2024.07.15
# 该变异为XX变异，导致XXX。【该（致病/可能致病性）变异可能与（药物疗效、XX癌预后、XX癌辅助诊断相关）。】
def cqxn_t116_evi_sum(var):
	clinic_inter = ""
	# snvindel-致病性和致癌性仅返回1个，这边可以用clinic_num_g判断是否为致病/疑似致病
	if var["bio_category"] == "Snvindel":
		clinic_inter = "致病/可能致病性" if var["clinic_num_g"] in [4, 5] else ""
	# cnv和sv返回两个，cnv和sv目前看function_classification
	elif var["bio_category"] in ["Cnv", "Sv"]:
		clinic_inter = "致病/可能致病性" if var["function_classification"] in ["Oncogenic", "Likely oncogenic"] else ""
	
	evi_type = []
	if "Predictive" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Predictive"]:
		evi_type.append("药物疗效")

	if "Prognostic" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Prognostic"]:
		tumor_lst = []
		for evi in var["evi_sum"]["evi_split"]["Prognostic"]:
			if "tumor_name_cn" in evi.keys() and evi["tumor_name_cn"] not in tumor_lst:
				tumor_lst.append(evi["tumor_name_cn"])
		evi_type.append("、".join(tumor_lst)+"预后")
	
	if "Diagnostic" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Diagnostic"]:
		tumor_lst = []
		for evi in var["evi_sum"]["evi_split"]["Diagnostic"]:
			if "tumor_name_cn" in evi.keys() and evi["tumor_name_cn"] not in tumor_lst:
				tumor_lst.append(evi["tumor_name_cn"])
		evi_type.append("、".join(tumor_lst)+"辅助诊断")

	return "该{0}变异可能与{1}相关。".format(clinic_inter, "、".join(evi_type))
jinja2.filters.FILTERS["cqxn_t116_evi_sum"] = cqxn_t116_evi_sum

# 重庆肿瘤116，NCCN指南推荐基因需要加上POLD1和POLE基因-2024.07.15
#def cqxn_116_nccn(info):
#	cdx_for_116 = info[0]
#	var_list = info[1]
#	if not set(["POLD1", "POLE"]) & set([i["gene_symbol"] for i in cdx_for_116]):
#		gene_dict = {
#			"POLD1" : "结直肠癌、小肠腺癌",
#			"POLE" : "结直肠癌、小肠腺癌"
#		}
#		detect_gene = []
#		for var in var_list:
#			for gene in set(re.split(",", var["gene_symbol"])):
#				if gene in ["POLD1", "POLE"]:
#					detect_gene.append(gene)
#					var["gene_symbol"] = gene
#					var["disease"] = gene_dict[gene]
#					var["clinic_num_s"] = S_level(var)
#					cdx_for_116.append(var)
#		for gene in set(["POLD1", "POLE"]) - set(detect_gene):
#			cdx_for_116.append({
#				"gene_symbol" : gene,
#				"disease" : gene_dict[gene]
#			})
#	return sorted(cdx_for_116, key=lambda i:i["gene_symbol"])
#jinja2.filters.FILTERS["cqxn_116_nccn"] = cqxn_116_nccn

# 吉林大学第一医院LC10-治疗方案重新排序下，根据A1>A2>A3……排序-2024.07.19
def jlyy_lc10_sort_evi(evi_sum):
	return sorted(evi_sum, key=lambda i:(i["evi_conclusion"], i["sense_rule"], i["regimen_name_py"].upper()))
jinja2.filters.FILTERS["jlyy_lc10_sort_evi"] = jlyy_lc10_sort_evi

# 吉林大学第一院LC10-检测结果汇总表-2024.07.17
def jlyy_lc10_summary(info):
	sum_type = info[2]
	if info[0]:
		var_list = [info[0]] + info[1]
	else:
		var_list = info[1]
	# 提示敏感/耐药
	def sense_resist(var_list, sense_or_resist):
		result = []
		for var in var_list:
			if "Predictive" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Predictive"]:
				regimen_list = ["{0}（{1}）".format(evi["regimen_name"], evi["evi_conclusion"][::-1]) for \
								evi in jlyy_lc10_sort_evi(var["evi_sum"]["evi_split"]["Predictive"]) if evi["evi_conclusion_simple"] != "D" and \
								evi["clinical_significance_cn"] == sense_or_resist]
				if regimen_list:
					var["regimen_sum_for_jdyy"] = "，".join(regimen_list)
					result.append(var)
		return result
	# 提示A级敏感
	def sense_A(var_list):
		A_result = []
		for var in var_list:
			# 变异信息----------------------------------------------------
			var_info = ""
			if var["var_id"] == "KRAS/NRAS/BRAF WT":
				var_info = "KRAS/NRAS/BRAF V600E 野生型"
			else:
				if var["bio_category"] == "Snvindel":
					if var["hgvs_p"] != "p.?":
						var_info = var["gene_symbol"] + " " + var["hgvs_p_ZJZL"]
					else:
						var_info = var["gene_symbol"] + " " + var["hgvs_c"]
				elif var["bio_category"] == "Cnv":
					var_info = var["gene_symbol"] + " 扩增"
				elif var["bio_category"] == "Sv":
					var_info = var["five_prime_gene"]+"-"+var["three_prime_gene"]+" 融合"
			# 变异信息-----------------------------------------------------
			sense_A_regimen = []
			resis_A_regimen = []
			A_sum = []
			if "Predictive" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Predictive"]:
				for evi in jlyy_lc10_sort_evi(var["evi_sum"]["evi_split"]["Predictive"]):
					if evi["evi_conclusion_simple"] == "A" and evi["clinical_significance_cn"] == "敏感" and evi["regimen_name"] not in sense_A_regimen:
						sense_A_regimen.append(evi["regimen_name"])
					if evi["evi_conclusion_simple"] == "A" and evi["clinical_significance_cn"] == "耐药" and evi["regimen_name"] not in sense_A_regimen:
						resis_A_regimen.append(evi["regimen_name"])
			if sense_A_regimen:
				A_sum.append("、".join(sense_A_regimen)+"敏感")
			if resis_A_regimen:
				A_sum.append("、".join(resis_A_regimen)+"耐药")
			if A_sum:
				A_result.append(var_info+" "+"，".join(A_sum)+"，A级推荐请结合临床")
		return "；".join(A_result)

	if sum_type == "sense":
		return sense_resist(var_list, "敏感")
	elif sum_type == "resist":
		return sense_resist(var_list, "耐药")
	elif sum_type == "sense_A":
		return sense_A(var_list)
jinja2.filters.FILTERS["jlyy_lc10_summary"] = jlyy_lc10_summary

# 吉林大学第一医院LC10-返回治疗指定证据（敏感/耐药、等级）-2024.07.17
def jlyy_lc10_get_evi(info):
	evi_sense = info[0]
	evi_level = info[1]
	evi_predictive = info[2]
	result = []
	if evi_predictive:
		for evi in evi_predictive:
			if evi["clinical_significance_cn"] == evi_sense and evi["evi_conclusion_simple"] in evi_level and evi["evi_interpretation"] not in result:
				result.append(evi["evi_interpretation"])
	return result
jinja2.filters.FILTERS["jlyy_lc10_get_evi"] = jlyy_lc10_get_evi

# 吉林大学第一医院LC10-返回治疗-敏感 A级证据总结-2024.07.17
def jlyy_lc10_sense_a_evi(var):
	var_info = ""
	if var["var_id"] == "KRAS/NRAS/BRAF WT":
		var_info = "KRAS/NRAS/BRAF V600E 野生型"
	else:
		if var["bio_category"] == "Snvindel":
			if var["hgvs_p"] != "p.?":
					var_info = var["gene_symbol"] + " " + var["hgvs_p_ZJZL"]
			else:
				var_info = var["gene_symbol"] + " " + var["hgvs_c"]
		elif var["bio_category"] == "Cnv":
			var_info = var["gene_symbol"] + " 扩增"
		elif var["bio_category"] == "Sv":
			var_info = var["five_prime_gene"]+"-"+var["three_prime_gene"]+" 融合"
	result = []
	evi_dict = {}
	level_a_evi = [evi for evi in var["evi_sum"]["evi_split"]["Predictive"] if evi["clinical_significance_cn"] == "敏感" and evi["evi_conclusion_simple"] == "A"]
	# 2024.11.05-按ABCD排序
	#level_a_evi = jlyy_lc10_sort_evi(level_a_evi)
	for evi in level_a_evi:
		if evi["refer_agency"] not in evi_dict.keys():
			evi_dict.setdefault(evi["refer_agency"], {})
		if evi["tumor_name_cn"] not in evi_dict[evi["refer_agency"]].keys():
			evi_dict[evi["refer_agency"]].setdefault(evi["tumor_name_cn"], [])
		if evi["regimen_name"] not in evi_dict[evi["refer_agency"]][evi["tumor_name_cn"]]:
			evi_dict[evi["refer_agency"]][evi["tumor_name_cn"]].append(evi["regimen_name"])
	for agency, info in evi_dict.items():
		for tumor, regimen_list in info.items():
			# 语句格式更新-2024.11.04
			#if agency in ["FDA", "NMPA"]:
			#	result.append("{0}适用于{1}变异的{2}患者（{3}批准）".format("、".join(regimen_list), var_info, tumor, agency))
			#else:
			#	result.append("{0}适用于{1}变异的{2}患者（{3}指南推荐）".format("、".join(regimen_list), var_info, tumor, agency))
			if agency in ["FDA", "NMPA"]:
				result.append("{0}批准{1}适用于{2}变异的{3}患者".format(agency, "、".join(regimen_list), var_info, tumor))
			else:
				result.append("{0}指南推荐{1}适用于{2}变异的{3}患者".format(agency, "、".join(regimen_list), var_info, tumor))
	return "；".join(result)
jinja2.filters.FILTERS["jlyy_lc10_sense_a_evi"] = jlyy_lc10_sense_a_evi

# 重庆肿瘤116基因-检测结果小结-具有临床意义的变异需要区分为10基因和其他基因-2024.07.25
def cqzl_116_var12_summary(info):
	all_var_list = info[0]
	gene_type = info[1]
	LC10_gene_list = ["ALK", "BRAF", "EGFR", "ERBB2", "KRAS", "NRAS", "MET", "PIK3CA", "RET", "ROS1"]
	var_list = []
	if gene_type == "inLC10":
		var_list = [var for var in all_var_list if set(re.split(",", var["gene_symbol"])) & set(LC10_gene_list)]
	elif gene_type == "outLC10":
		var_list = [var for var in all_var_list if not set(re.split(",", var["gene_symbol"])) & set(LC10_gene_list)]
	# 2024.11.08-加一个全部返回的
	else:
		var_list = all_var_list
	
	for var in var_list:
		if var["bio_category"] == "Snvindel":
			var["cqzl_var_info"] = var["hgvs_p"] if var["hgvs_p"] != "p.?" else var["hgvs_c"]
		elif var["bio_category"] == "Cnv":
			var["cqzl_var_info"] = "扩增"
		elif var["bio_category"] == "Sv":
			var["cqzl_var_info"] = "融合"
		else:
			var["cqzl_var_info"] = ""
	
	# 每个变异（除了最后一个）cqzl_var_info后面加个逗号，在模板中使用for循环展示变异
	if len(var_list) > 1:
		for var in var_list[0:-1]:
			var["cqzl_var_info"] = var["cqzl_var_info"] + "，"
	return var_list
jinja2.filters.FILTERS["cqzl_116_var12_summary"] = cqzl_116_var12_summary

# 复旦中山200基因，汇总，其中变异涉及基因数190个，其中一个EBV（不展示了），化疗6个，MSI 1个，这边只汇总190-2024.07.26
def fdzs_190_gene_sum(var_list):
	gene_list = ["ABRAXAS1", "AKT1", "AKT2", "AKT3", "ALK", "APC", "AR", "ARAF", "ARID1A", "ARID1B", \
				 "ARID2", "ATM", "ATR", "ATRX", "AURKA", "AXIN1", "B2M", "BAP1", "BARD1", "BMPR1A", \
				 "BRAF", "BRCA1", "BRCA2", "BRD3", "BRD4", "BRIP1", "CCND1", "CCNE1", "CD274", "CDH1", \
				 "CDK12", "CDK4", "CDK6", "CDKN1B", "CDKN2A", "CDKN2B", "CHD1", "CHEK1", "CHEK2", "CLDN18", \
				 "CLIP1", "CTNNB1", "CYP17A1", "DDR2", "DICER1", "DNAJB1", "DNMT3A", "EGFR", "EIF1AX", \
				 "ELK4", "EPCAM", "ERBB2", "ERBB3", "ERCC2", "ERG", "ESR1", "ETV1", "ETV4", "ETV5", \
				 "EZH2", "FANCA", "FANCC", "FANCD2", "FANCL", "FAT1", "FBXW7", "FGF19", "FGFR1", "FGFR2", \
				 "FGFR3", "FGFR4", "FH", "FLT3", "FOXA1", "GATA3", "GATA6", "GEN1", "GLI1", "GLI2", \
				 "GNA11", "GNA14", "GNAQ", "GNAS", "GTF2I", "H3-3A", "H3C2", "H3C3", "HDAC2", "HOXB13", \
				 "HRAS", "HSD3B1", "IDH1", "IDH2", "IFNGR1", "IFNGR2", "JAK1", "JAK2", "KDM5C", "KDM6A", \
				 "KEAP1", "KIT", "KMT2B", "KMT2C", "KMT2D", "KRAS", "MAP2K1", "MAP2K4", "MDM2", "MDM4", \
				 "MET", "MLH1", "MLH3", "MMS22L", "MRE11", "MSH2", "MSH3", "MSH6", "MTAP", "MTOR", \
				 "MUTYH", "MYB", "MYC", "MYCN", "NBN", "NF1", "NF2", "NFE2L2", "NKX2-1", "NRAS", \
				 "NRG1", "NSD3", "NTRK1", "NTRK2", "NTRK3", "NUTM1", "PALB2", "PAX8", "PBRM1", "PDCD1LG2", \
				 "PDGFRA", "PIK3CA", "PIK3R1", "PMS2", "POLD1", "POLE", "PPP2R1A", "PRKACA", "PTCH1", "PTCH2", \
				 "PTEN", "QKI", "RAD50", "RAD51B", "RAD51C", "RAD51D", "RAD54L", "RAF1", "RASA1", "RASGRF1", \
				 "RB1", "RBM10", "RET", "RICTOR", "RIT1", "RNF43", "ROS1", "SETD2", "SF3B1", "SMAD4", \
				 "SMARCA2", "SMARCA4", "SMARCB1", "SMO", "SPOP", "STAT3", "STK11", "SUFU", "TERT", "TET2", \
				 "TGFBR2", "TP53", "TSC1", "TSC2", "VEGFA", "VHL", "YAP1", "ZFTA", "ZNF532", "ZNF592"]

	result = []
	
	detect_gene = []
	for var in var_list:
		for gene in set(re.split(",", var["gene_symbol"])):
			detect_gene.append(gene)
			if gene in gene_list:
				result.append(var)
	for gene in set(gene_list) - set(detect_gene):
		result.append({"gene_symbol" : gene})
	
	return sorted(result, key = lambda i:i["gene_symbol"])
jinja2.filters.FILTERS["fdzs_190_gene_sum"] = fdzs_190_gene_sum

# 复旦中山200基因，需要展示指定化疗位点的基因型-2024.07.30
def fdzs_chemo_result(info):
	chemo = info[0]
	dbsnp = info[1]
	chemo_dict = {}
	for i in chemo:
		chemo_dict[i["dbsnp"]] = i["genotype"]
	if dbsnp in chemo_dict.keys():
		return chemo_dict[dbsnp]
	else:
		return ""
jinja2.filters.FILTERS["fdzs_chemo_result"] = fdzs_chemo_result


# 山东齐鲁CP40-胃肠道间质瘤/黑色素瘤-变异需要过滤掉非检测范围内的-2024.08.07
def sdql_cp_mel_gist_var_filter(info):
	var_list = info[0]
	gene_list = info[1]
	return [var for var in var_list if var["gene_symbol"] in gene_list and var["bio_category"] == "Snvindel"]
jinja2.filters.FILTERS["sdql_cp_mel_gist_var_filter"] = sdql_cp_mel_gist_var_filter
	
# 山东齐鲁CP40-胃肠道间质瘤/黑色素瘤-结果汇总，没检测到变异的基因也要展示-2024.08.07
def sdql_cp_mel_gist_var_summary(info):
	var_list = info[0]
	gene_list = info[1]
	detect_gene = [var["gene_symbol"] for var in var_list]
	result = [var for var in var_list if var["gene_symbol"] in gene_list and var["bio_category"] == "Snvindel"]
	for gene in set(gene_list) - set(detect_gene):
		result.append({"gene_symbol" : gene})
	return sorted(result, key = lambda i:i["gene_symbol"])
jinja2.filters.FILTERS["sdql_cp_mel_gist_var_summary"] = sdql_cp_mel_gist_var_summary

# 福建附一116检测结果第一段数据质控总结-2024.08.08
def judge_fjfy_116_summary(info):
	prod_names = info[3]
	# 文库质量-文库总量
	# 文库总量字段改为dna_pre_library_qty-2024.08.14
	#library_qty = info[0]["library_qty"] if info[0] and "library_qty" in info[0].keys() and info[0]["library_qty"] else ""
	library_qty = info[0]["dna_pre_library_qty"] if info[0] and "dna_pre_library_qty" in info[0].keys() and info[0]["dna_pre_library_qty"] else ""
	lib_judge = "未填写" if not library_qty else "合格" if float(library_qty) > 900 else "警戒"
	# 测序质量
	dna_data_qc = info[1]
	qc_gradient = info[2]
	ngs_judge = "合格" if float(qc_gradient["coverage_ratio_uniq_hot_180"].replace("%", "")) >= 0.95 and \
						  float(dna_data_qc["cleandata_q30_num"]) >= 0.75 and \
						  float(dna_data_qc["cover_ratio_num"]) >= 0.95 else \
				"警戒"
	result = [
		{"type" : "文库质量", "info" : lib_judge},
		{"type" : "测序质量", "info" : ngs_judge}
	]
	# 排序，1 合格>未填写>警戒，2 文库质量>测序质量
	type_rule = ["文库质量", "测序质量"]
	info_rule = ["合格", "未填写", "警戒"]
	result = sorted(result, key = lambda i:(info_rule.index(i["info"]), type_rule.index(i["type"])))
	true_item = [i["type"] for i in result if i["info"] == "合格"]
	false_item = [i["type"] for i in result if i["info"] == "警戒"]
	bank_item = [i["type"] for i in result if i["info"] == "未填写"]
	# 文字组装
	tmp = []
	if true_item:
		if len(true_item) == 1:
			tmp.append(true_item[0]+"合格")
		else:
			tmp.append("、".join(true_item)+"均合格")
	if bank_item:
		tmp.append("、".join(bank_item)+"（数据未填写，无法判断）")
	if false_item:
		if len(false_item) == 1:
			tmp.append(false_item[0]+"警戒")
		else:
			tmp.append("、".join(false_item)+"均为警戒")
	return "，".join(tmp)
jinja2.filters.FILTERS["judge_fjfy_116_summary"] = judge_fjfy_116_summary

# 聊城人民BPTM Plus，展示每个基因的检出结果-2024.08.12（跟毓璜顶不同的是聊城需要展示III类变异）
def lcrm_gene_detect_bptmplus(info):
	gene = info[3]
	var_list_I = [var for var in info[0] if var["gene_symbol"] == gene]
	var_list_II = [var for var in info[1] if var["gene_symbol"] == gene]
	var_list_III = [var for var in info[2] if var["gene_symbol"] == gene]
	result = []
	def get_class_result(level, v_list):
		tmp = []
		tmp.append("检出{0}个{1}变异".format(str(len(v_list)), level))
		for var in v_list:
			if var["hgvs_p"] != "p.?":
				tmp.append(var["gene_region"]+" "+var["hgvs_c"]+" "+var["hgvs_p"])
			else:
				tmp.append(var["gene_region"]+" "+var["hgvs_c"])
		return "，".join(tmp)
	if var_list_I:
		result.append(get_class_result("I类", var_list_I))
	if var_list_II:
		result.append(get_class_result("II类", var_list_II))
	if var_list_III:
		result.append(get_class_result("III类", var_list_III))
	return result
jinja2.filters.FILTERS["lcrm_gene_detect_bptmplus"] = lcrm_gene_detect_bptmplus

# 聊城人民HRD，展示每个基因的检测情况，未检测到变异的基因也要展示-2024.08.12
def lcrm_hrd_summary(var_list):
	gene_list = ["ATM", "BARD1", "BRIP1", "CDH1", "CDK12", "CHEK1", "CHEK2", "FANCA", "FANCL", \
				 "HDAC2", "PALB2", "PPP2R2A", "PTEN", "RAD51B", "RAD51C", "RAD51D", "RAD54L", "TP53"]
	appr_gene = ["ATM", "BARD1", "BRIP1", "CDK12", "CHEK1", "CHEK2", "FANCA", \
				 "FANCL", "PALB2", "RAD51B", "RAD51C", "RAD51D", "RAD54L"]
	result = [var for var in var_list if var["gene_symbol"] in gene_list]
	detect_gene = [var["gene_symbol"] for var in result]
	for gene in set(gene_list) - set(detect_gene):
		result.append({"gene_symbol" : gene})
	for var in result:
		var["appr"] = "T" if var["gene_symbol"] in appr_gene else ""
	return sorted(result, key = lambda i:i["gene_symbol"])
jinja2.filters.FILTERS["lcrm_hrd_summary"] = lcrm_hrd_summary

# 吉林大学第一医院-LC10-本癌种重要基因要根据癌种进行区分-2024.08.21
def jlyy_lc10_important_gene(info):
	withoutvar_genelist = info[0]
	tumor_list = info[1]
	if "肺癌" in tumor_list:
		return "，".join(sorted(list(set(withoutvar_genelist) & set(["ALK", "BRAF", "EGFR", "ERBB2", "KRAS", "MET", "RET", "ROS1"]))))
	elif "肠癌" in tumor_list:
		return "，".join(sorted(list(set(withoutvar_genelist) & set(["KRAS", "NRAS", "BRAF"]))))
	else:
		return "，".join(withoutvar_genelist)
jinja2.filters.FILTERS["jlyy_lc10_important_gene"] = jlyy_lc10_important_gene

# 北大一HRR汇总-需要展示每个基因的结果，包含转录本、检测区域、检测变异类型等-2024.08.22
def bdy_hrr_summary(var_list):
	hrr_dict = {
		"AR" : ("NM_000044", "Exon 1~Exon 8、外显子-内含子连接区", "点突变、插入/缺失"),
		"ATM" : ("NM_000051", "Exon 2~Exon 63、外显子-内含子连接区", "点突变、插入/缺失"),
		"ATR" : ("NM_001184", "Exon 1~Exon 47、外显子-内含子连接区", "点突变、插入/缺失"),
		"BARD1" : ("NM_000465", "Exon 1~Exon 11、外显子-内含子连接区", "点突变、插入/缺失"),
		"BRAF" : ("NM_004333", "Exon 11/12/15/18 热点区域", "点突变、插入/缺失"),
		"BRCA1" : ("NM_007294", "Exon 2~Exon 3、 Exon 5~Exon 24、外显子-内含子连接区、部分内含子 和 UTR 区", "点突变、插入/缺失"),
		"BRCA2" : ("NM_000059", "Exon 2~Exon 27、外显子-内含子连接区、部分内含子和 UTR 区", "点突变、插入/缺失"),
		"BRIP1" : ("NM_032043", "Exon 2~Exon 20、外显子-内含子连接区", "点突变、插入/缺失"),
		"CDH1" : ("NM_004360", "Exon 1~Exon 16、外显子-内含子连接区", "点突变、插入/缺失"),
		"CDK12" : ("NM_016507", "Exon 1~Exon 14、外显子-内含子连接区", "点突变、插入/缺失"),
		"CHEK1" : ("NM_001274", "Exon 2~Exon 13、外显子-内含子连接区", "点突变、插入/缺失"),
		"CHEK2" : ("NM_007194", "Exon 2~Exon 15、外显子-内含子连接区", "点突变、插入/缺失"),
		"ERBB2" : ("NM_004448", "Exon 1~3/6~9/11~12/16~21/23/25~27 热点区域", "点突变、插入/缺失"),
		"ESR1" : ("NM_001122740", "Exon 2~Exon 9、外显子-内含子连接区", "点突变、插入/缺失"),
		"FANCA" : ("NM_000135", "Exon 1~Exon 43、外显子-内含子连接区", "点突变、插入/缺失"),
		"FANCL" : ("NM_018062", "Exon 1~Exon 14、外显子-内含子连接区", "点突变、插入/缺失"),
		"HDAC2" : ("NM_001527", "Exon 1~Exon 14、外显子-内含子连接区", "点突变、插入/缺失"),
		"HOXB13" : ("NM_006361", "Exon 1~Exon 2、外显子-内含子连接区", "点突变、插入/缺失"),
		"KRAS" : ("NM_033360", "Exon 2~4 热点区域", "点突变、插入/缺失"),
		"MRE11" : ("NM_005591", "Exon 2~Exon 20、外显子-内含子连接区", "点突变、插入/缺失"),
		"NBN" : ("NM_002485", "Exon 1~Exon 16、外显子-内含子连接区", "点突变、插入/缺失"),
		"NRAS" : ("NM_002524", "Exon 2~4 热点区域", "点突变、插入/缺失"),
		"PALB2" : ("NM_024675", "Exon 1~Exon 13、外显子-内含子连接区", "点突变、插入/缺失"),
		"PIK3CA" : ("NM_006218", "Exon 2/5/6/8/10/14/21 热点区域", "点突变、插入/缺失"),
		"PPP2R2A" : ("NM_002717", "Exon 2~Exon 10、外显子-内含子连接区", "点突变、插入/缺失"),
		"PTEN" : ("NM_000314", "Exon 1~Exon 9、外显子-内含子连接区", "点突变、插入/缺失"),
		"RAD51B" : ("NM_133509", "Exon 2~Exon 11、外显子-内含子连接区", "点突变、插入/缺失"),
		"RAD51C" : ("NM_058216", "Exon 1~Exon 9、外显子-内含子连接区", "点突变、插入/缺失"),
		"RAD51D" : ("NM_002878", "Exon 1~2/ 6~10、外显子-内含子连接区", "点突变、插入/缺失"),
		"RAD54L" : ("NM_001142548", "Exon 2~Exon 19、外显子-内含子连接区", "点突变、插入/缺失"),
		"STK11" : ("NM_000455", "Exon 1~Exon 9、外显子-内含子连接区", "点突变、插入/缺失"),
		"TP53" : ("NM_000546", "Exon 2~Exon 11、外显子-内含子连接区", "点突变、插入/缺失")
	}
	detect_gene = [var["gene_symbol"] for var in var_list]
	for var in var_list:
		var["detect_region"] = hrr_dict[var["gene_symbol"]][1] if var["gene_symbol"] in hrr_dict.keys() else ""
		var["detect_type"] = hrr_dict[var["gene_symbol"]][2] if var["gene_symbol"] in hrr_dict.keys() else ""
		if "transcript_primary_simple" not in var.keys():
			var["transcript_primary_simple"] = hrr_dict[var["gene_symbol"]][0] if var["gene_symbol"] in hrr_dict.keys() else ""
	for gene in set([k for k in hrr_dict.keys()]) - set(detect_gene):
		var_list.append({
			"gene_symbol" : gene,
			"transcript_primary_simple" : hrr_dict[gene][0],
			"detect_region" : hrr_dict[gene][1],
			"detect_type" : hrr_dict[gene][2]
		})
	return sorted(var_list, key = lambda i:i["gene_symbol"])
jinja2.filters.FILTERS["bdy_hrr_summary"] = bdy_hrr_summary

# 获取分析包版本，如v1.0.1 -2024.08.22
def get_analysis_version_ynzl(json_batch_name):
	result = ""
	if json_batch_name:
		json_batch_list = re.split("_", json_batch_name)
		# zip上传包：日期_测序仪_测序芯片_生信分析包_版本_系统日期流水号
		if len(json_batch_list) == 6:
			result = json_batch_list[4]
		# export上传包：日期_测序仪_测序芯片_生信分析包_版本.export系统日期流水号
		elif len(json_batch_list) == 5:
			result = re.split("\.export|\.zip", json_batch_list[-1])[0]
	return result
jinja2.filters.FILTERS["get_analysis_version_ynzl"] = get_analysis_version_ynzl

# 云肿遗传150-检测结果小结-2024.08.23
def ynzl_150_summary(var_list):
	# 基于本次送检样本，【RAD51D基因检测到1个疑似致病性变异】，建议至肿瘤基因门诊进行遗传咨询
	result = []
	detect_gene = []
	for var in var_list:
		if var["gene_symbol"] not in detect_gene:
			detect_gene.append(var["gene_symbol"])

	def get_single_gene(var_list, gene):
		summary = []
		gene_var_list_4 = [var for var in var_list if var["gene_symbol"] == gene and var["clinic_num_g"] == 4]
		gene_var_list_5 = [var for var in var_list if var["gene_symbol"] == gene and var["clinic_num_g"] == 5]
		if gene_var_list_5:
			summary.append("{0}个致病性变异".format(str(len(gene_var_list_5))))
		if gene_var_list_4:
			summary.append("{0}个疑似致病性变异".format(str(len(gene_var_list_4))))
		return "和".join(summary)

	for gene in detect_gene:
		single_gene_summary = get_single_gene(var_list, gene)
		if single_gene_summary:
			result.append("{0}基因检测到{1}".format(gene, single_gene_summary))
	return "、".join(result)
jinja2.filters.FILTERS["ynzl_150_summary"] = ynzl_150_summary

# 云南肿瘤判断3类是否存在具有遗传相关证据的变异-2024.08.23
def ynzl_150_judge_level3_Predisposing(var_list):
	judge_result = False
	for var in var_list:
		if "evi_split" in var["evi_sum"].keys() and var["evi_sum"]["evi_split"] and \
			"Predisposing" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Predisposing"]:
			judge_result = True
	return judge_result
jinja2.filters.FILTERS["ynzl_150_judge_level3_Predisposing"] = ynzl_150_judge_level3_Predisposing

# 云肿gBRCA-判断变异是否为获批的，最高等级A/B为获批，最高等级C/D为非获批-2024.08.23
def brca_judge_inApprovalTumor(var_brca):
	var = var_brca["snv_s"]["B1_L4"] + var_brca["snv_s"]["B1_L5"] + \
		  var_brca["snv_s"]["B2_L4"] + var_brca["snv_s"]["B2_L5"] + \
		  var_brca["mlpa"]["B1_Loss"] + var_brca["mlpa"]["B2_Loss"]
	level_list = []
	for i in var:
		# 这边若治疗方案为空，会被判定成潜在临床意义！！
		if "evi_sum" in i.keys() and i["evi_sum"] and "evi_split" in i["evi_sum"].keys() and i["evi_sum"]["evi_split"] and "Predictive" in i["evi_sum"]["evi_split"].keys():
			level_list += [j["evi_conclusion_simple"] for j in i["evi_sum"]["evi_split"]["Predictive"]]
	if level_list:
		if set(["A", "B"]) & set(level_list):
			return "intumor"
		else:
			return "outtumor"
jinja2.filters.FILTERS["brca_judge_inApprovalTumor"] = brca_judge_inApprovalTumor

# 湘雅二院HRD，小结新增表格，展示HRD状态、BRCA基因和其他基因检测情况-2024.08.27
# 列1需要合并，模板中实现会有问题
def xyey_hrd_merge(var_dict):
	brca_var = var_dict["ec_type"]["BRCA1_level12"] + var_dict["ec_type"]["BRCA2_level12"]
	all_var = var_dict["var_somatic"]["level_I"] + var_dict["var_somatic"]["level_II"]
	result = []
	if brca_var:
		result.append({
			"type" : "BRCA",
			"var_list" : brca_var
		})
	other_var = [var for var in all_var if var["gene_symbol"] not in ["BRCA1", "BRCA2"]]
	if other_var:
		result.append({
			"type" : "other",
			"var_list" : other_var
		})
	return result
jinja2.filters.FILTERS["xyey_hrd_merge"] = xyey_hrd_merge

# 复旦中山MP+肉瘤-肉瘤指南/共识推荐基因检测结果-2024.08.27
def fdzs_mp_remap_sum(info):
	var_list = info[0]["var_somatic"]["level_I"] + info[0]["var_somatic"]["level_II"]
	tumor = info[1]
	tumor_biomarker_dict = {
		"非典型脂肪瘤样肿瘤_高分化_去分化脂肪肉瘤" : ["CDK4_cnv", "MDM2_cnv"],
		"黏液样_圆形细胞_脂肪肉瘤" : ["FUS-DDIT3", "EWSR1-DDIT3"],
		"孤立性纤维性肿瘤" : ["NAB2-STAT6"],
		"炎性成肌纤维细胞瘤" : ["TPM3-ALK", "TPM4-ALK", "CLTC-ALK", "RANBP2-ALK", "ATIC-ALK", \
				 			  "CARS-ALK", "SEC31L1-ALK", "PPFIBP1-ALK", "TIMP3-ALK", "IGFBP5-ALK", \
							  "THBS1-ALK", "RRBP1-ALK", "ETV6-NTRK3", "TFG-ROS1", "YWHAE-ROS1"],
		"隆突性皮肤纤维肉瘤_巨细胞成纤维细胞瘤" : ["COL1A1-PDGFB", "COL6A3-PDGFD", "EMILIN2-PDGFD"],
		"婴儿型纤维肉瘤" : ["ETV6-NTRK3"],
		"低级别纤维黏液样肉瘤" : ["FUS-CREB3L1", "FUS-CREB3L2"],
		"硬化性上皮样纤维肉瘤" : ["EWSR1-CREB3L1", "FUS-CREB3L1", "FUS-CREB3L2", "YAP1-KMT2A"],
		"血管瘤样纤维组织细胞瘤" : ["EWSR1-ATF1", "EWSR1-CREB1", "FUS-ATF1"],
		"侵袭性纤维瘤" : ["APC_all", "CTNNB1_all"],
		"腱鞘巨细胞瘤" : ["CSF1_sv"],
		"上皮样血管内皮瘤" : ["WWTR1-CAMTA1", "YAP1-TFE3", "MBNL1-FOS", "VIM-FOS", "ZFP36-FOSB", \
							 "WWTR1-FOSB", "ACTB-FOSB", "SERPINE1-FOSB", "YAP1-MAML2", "PTBP1-MAML2"],
		"血管肉瘤" : ["MYC_cnv"],
		"腺泡状横纹肌肉瘤" : ["PAX3-FOXO1", "PAX7-FOXO1", "PAX3-FOXO4", "PAX3-NCOA1", "PAX3-NCOA2", "FOXO1-FGFR1", "PAX3-INO80D"],
		"先天性_婴儿型梭形细胞_横纹肌肉瘤" : ["SRF-NCOA2", "TEAD1-NCOA2", "VGLL2-NCOA2", "VGLL2-CITED2"],
		"成人梭形细胞_硬化性横纹肌肉瘤" : ["MYOD1_all"],
		"胚胎性横纹肌肉瘤" : ["BCOR_all", "FBXW7_all", "FGFR4_all", "HRAS_all", "KRAS_all", \
							 "MYOD1_all", "NF1_all", "NRAS_all", "PIK3CA_all", "TP53_all"],
		"间叶性软骨肉瘤" : ["HEY1-NCOA2"],
		"恶性周围神经鞘膜瘤" : ["NF1_all", "CDKN2A_all", "CDKN2B_snvidel", "EED_all", "SUZ12_all"],
		"恶性色素性神经鞘膜肿瘤" : ["PRKAR1A_all"],
		"软组织肌上皮肿瘤" : ["EWSR1-POU5F1", "EWSR1-PBX1", "FUS-KLF17", "EWSR1-PBX3", "EWSR1-ZNF444"],
		"NTRK重排梭形细胞肿瘤" : ["NTRK1_sv", "NTRK2_sv", "NTRK3_sv"],
		"滑膜肉瘤" : ["SS18-SSX1", "SS18-SSX2", "SS18-SSX4", "SS18L1-SSX1"],
		"上皮样肉瘤" : ["SMARCB1_all"],
		"腺泡状软组织肉瘤" : ["ASPSCR1-TFE3"],
		"软组织透明细胞肉瘤" : ["EWSR1-ATF1", "EWSR1-CREB1"],
		"软骨肉瘤" : ["IDH1_all", "IDH2_all"],
		"骨外黏液样软骨肉瘤" : ["EWSR1-NR4A3", "TAF15-NR4A3", "TCF12-NR4A3", "TFG-NR4A3"],
		"促结缔组织增生性小圆细胞肿瘤" : ["EWSR1-WT1"],
		"肾外横纹肌样瘤" : ["SMARCB1_all"],
		"未分化圆细胞肉瘤" : ["CIC-DUX4", "BCOR-CCNB3"],
		"内膜肉瘤" : ["CDK4_cnv", "MDM2_cnv"],
		"血管周上皮样细胞肿瘤" : ["TSC2_all", "TFE3_sv", "RAD51B_sv", "HTR4-ST3GAL1"],
		"尤因肉瘤_周围神经_外胚层肿瘤" : ["EWSR1-FLI1", "EWSR1-ERG", "EWSR1-FEV", "EWSR1-ETV1", "EWSR1-ETV4", "EWSR1-PATZ1", "FUS-ERG", "FUS-FEV"],
		"CIC重排肉瘤" : ["CIC-DUX4", "CIC-FOXO4", "CIC-NUTM1", "CIC-NUTM2B"],
		"BCOR重排肉瘤" : ["BCOR-CCNB3", "BCOR-MAML3", "ZC3H7B-BCOR", "YWHAE-NUTM2B"],
		"婴幼儿未分化圆细胞肉瘤_婴幼儿原始黏液样间叶性肿瘤" : ["YWHAE-NUTM2B"],
		"EWSR1-非ETS融合的圆细胞肉瘤" : ["EWSR1-NFATC2", "EWSR1-SP3", "EWSR1-POU5F1", "EWSR1-PATZ1", "EWSR1-SMARCA5", "FUS-NFATC2"],
		"胃肠道间质瘤_Carney-Stratakis综合征" : ["KIT_all", "PDGFRA_all", "SDHB_all", "SDHC_all", "SDHD_all"],
		"上皮样纤维组织细胞瘤" : ["SQSTM1-ALK", "VCL-ALK"],
		"骨上皮样和梭形细胞横纹肌肉瘤" : ["EWSR1-TFCP2", "FUS-TFCP2", "MEIS1-NCOA2"],
		"GLI1扩增_重排恶性上皮样肿瘤" : ["ACTB-GLI1", "MALAT1-GLI1", "PTCH1-GLI1"],
		"胃母细胞瘤_胃丛状纤维黏液瘤" : ["MALAT1-GLI1"],
		"富细胞性肌样肿瘤" : ["SRF-ICA1L"],
		"肌上皮瘤样玻璃样变肿瘤" : ["OGT-FOXO3"],
		"肾脏原始梭形细胞肉瘤" : ["MEIS1-NCOA2"],
		"部分不能分类的圆细胞肉瘤" : ["EWSR1-CREB1"],
		"血管球瘤" : ["MIR143-NOTCH1", "MIR143-NOTCH2", "MIR143-NOTCH3"],
		"其他双表达CD34和S-100的梭形细胞肿瘤" : ["PDZRN3-RAF1", "SLMAP-RAF1", "TMF1-RAF1", "MTAP-RAF1", "TFG-RET", "MYH10-RET", \
							  				   "NCOA4-RET", "VCL-RET", "CLIP2-RET", "CCDC6-RET", "KHDRBS1-RET", "SPECC1L-RET", \
											   "KIAA1217-RET", "PPP1CB-ALK", "CUX1-BRAF", "SEPTIN7-BRAF", "CDC42SE2-BRAF"],
		"上皮样平滑肌肉瘤" : ["PGR_sv"],
		"黏液样平滑肌肉瘤" : ["PLAG1_sv"],
		"动脉瘤样骨囊肿" : ["CDH11-USP6", "CNBP-USP6", "COL1A1-USP6", "MYH9-USP6", "OMD-USP6", "THRAP3-USP6"],
		"高级别子宫内膜间质肉瘤" : ["YWHAE-NUTM2A", "YWHAE-NUTM2B", "ZC3H7B-BCOR"],
		"低级别子宫内膜间质肉瘤" : ["JAZF1-SUZ12", "JAZF1-PHF1", "MEAF6-PHF1", "EPC1-PHF1", "MBTD1-CXorf67", "BRD8-PHF1", "EPC2-PHF1", "EPC1-SUZ12"],
		"类似于卵巢性索肿瘤的子宫肿瘤" : ["ESR1_sv", "GREB1_sv"],
		"Mullerian腺肉瘤" : ["NCOA2_sv", "NCOA3_sv"]
	}

	# 1. 获取gene融合的结果
	def get_sv(var_list, biomarker):
		return [biomarker for var in var_list if var["bio_category"] in ["Sv", "PSeqRnaSv"] and biomarker in list(re.split(",", var["gene_symbol"]))]
	# 2. 获取gene基因扩增的结果
	def get_cnv(var_list, biomarker):
		return [biomarker for var in var_list if var["bio_category"] in ["Cnv"] and biomarker == var["gene_symbol"]]
	# 3. 获取gene基因变异的结果
	def get_all(var_list, biomarker):
		return [biomarker for var in var_list if biomarker in list(re.split(",", var["gene_symbol"]))]
	# 4. 获取gene1-gene2融合的结果
	def get_sv_detail(var_list, biomarker):
		return [biomarker for var in var_list if var["bio_category"] in ["Sv", "PSeqRnaSv"] and \
		  									re.split("-", biomarker)[0] == var["five_prime_gene"] and \
											re.split("-", biomarker)[1] == var["three_prime_gene"]]

	# 获取肿瘤对应的分子标志物
	biomarker_list = tumor_biomarker_dict.get(tumor, [])
	# 对是否有检测到分子标志物进行判断
	biomarker_result = {
		"sv_detail" : [],
		"sv" : [],
		"cnv" : [],
		"all" : []
	}
	for biomarker in biomarker_list:
		if "-" in biomarker:
			if get_sv_detail(var_list, biomarker):
				biomarker_result["sv_detail"].extend(get_sv_detail(var_list, biomarker))
		elif re.split("_", biomarker)[1] == "all":
			if get_all(var_list, re.split("_", biomarker)[0]):
				biomarker_result["all"].extend(get_all(var_list, re.split("_", biomarker)[0]))
		elif re.split("_", biomarker)[1] == "cnv":
			if get_cnv(var_list, re.split("_", biomarker)[0]):
				biomarker_result["cnv"].extend(get_cnv(var_list, re.split("_", biomarker)[0]))
		elif re.split("_", biomarker)[1] == "sv":
			if get_sv(var_list, re.split("_", biomarker)[0]):
				biomarker_result["sv"].extend(get_sv(var_list, re.split("_", biomarker)[0]))
	# 对分子标志物检测结果进行整理
	# 返回[gene基因变异, gene融合, gene1-gene2融合, gene基因扩增]	
	biomarker_result_sum = []
	if biomarker_result["all"]:
		biomarker_result_sum.append({"gene" : "、".join(biomarker_result["all"]), "add_info" : "基因变异"})
	if biomarker_result["sv"]:
		biomarker_result_sum.append({"gene" : "、".join(biomarker_result["sv"]), "add_info" : "融合"})
	if biomarker_result["sv_detail"]:
		biomarker_result_sum.append({"gene" : "、".join(biomarker_result["sv_detail"]), "add_info" : "融合"})
	if biomarker_result["cnv"]:
		biomarker_result_sum.append({"gene" : "、".join(biomarker_result["cnv"]), "add_info" : "基因扩增"})
	
	return biomarker_result_sum
jinja2.filters.FILTERS["fdzs_mp_remap_sum"] = fdzs_mp_remap_sum

# --------------------这个功能先不使用，使用CNV时F掉MLPA结果，使用MLPA时F掉CNV结果就好--------------------------
# BRCA CNV和MLPA-治疗方案/药物介绍删除相关变异等-2024.09.02
# 默认只会存在一个del或一个dup
def brca_cnv_mlpa_delete_drug_var(info):
	# MLPA结果
	mlpa_var = info[0]
	# gCNV结果
	gcnv_var = info[1]
	# 大片段变异结果以哪个为主，mlpa/gcnv
	judge_brca_cnv = info[2]
	# 治疗方案/药物介绍
	regimen_list = info[3]

	# 汇总mlpa的分子标志物
	mlpa_biomarker = []
	for var in mlpa_var:
		mlpa_type = "del" if var["type"] == "Loss" else "dup" if var["type"] == "Gain" else ""
		key_ = {
			"biomarker_type" : var["gene_symbol"]+" "+var["value"]+" "+mlpa_type
		}
		if key_ not in mlpa_biomarker:
			mlpa_biomarker.append(key_)
	# ("mlpa_biomarker", mlpa_biomarker)
	# 汇总gcnv的分子标志物，biomarker_type/cnv_type两种都会返回（看是否检出MLPA及drug/regimen有区别）
	gcnv_biomarker_key1 = []
	gcnv_biomarker_key2 = []
	for var in gcnv_var:
		gcnv_type = "del" if var["type"] == "Loss" else "dup" if var["type"] == "Gain" else ""
		key1_ = {
			"biomarker_type" : var["gene_symbol"]+" "+var["value"]+" "+gcnv_type
		}
		if key1_ not in gcnv_biomarker_key1:
			gcnv_biomarker_key1.append(key1_)

		key2_ = {
			"gene_symbol" : var["gene_symbol"],
			"cnv_type" : var["cnv_type"]
		}
		if key2_ not in gcnv_biomarker_key2:
			gcnv_biomarker_key2.append(key2_)
	# ("gcnv_biomarker 1", gcnv_biomarker_key1)
	#print ("gcnv_biomarker 2", gcnv_biomarker_key2)
	# 结果以mlpa为准时，需要删除gcnv的分子标志物（bio1 需要跟mlpa比较，考虑要不要删，bio2直接删了吧）
	if judge_brca_cnv == "mlpa":
		delete_bio = []
		for bio in gcnv_biomarker_key2:
			delete_bio.append(bio)
		for bio in gcnv_biomarker_key1:
			if bio not in mlpa_biomarker:
				delete_bio.append(bio)
	# 结果以gcnv为准时，需要删除mlpa的分子标志物
	if judge_brca_cnv == "gcnv":
		delete_bio = []
		# 若存在gcnv
		# 若存在mlpa（共检），则删除gcnv bio2；若不存在mlpa，则删除gcnv bio1
		if gcnv_var:
			if mlpa_var:
				delete_bio.extend(gcnv_biomarker_key2)
			else:
				delete_bio.extend(gcnv_biomarker_key1)
		# 不存在gcnv，则mlpa的都删了
		else:
			delete_bio.extend(mlpa_biomarker)
	#print ("delete_bio", delete_bio)
	for regimen in regimen_list:
		for a in regimen["var"]:
			if a in delete_bio:
				regimen["var"].remove(a)
	return_regimen = [i for i in regimen_list if i["var"]]
	return return_regimen
jinja2.filters.FILTERS["brca_cnv_mlpa_delete_drug_var"] = brca_cnv_mlpa_delete_drug_var
#----------------------------先不用--------------------------------------------------------------------

# 重庆附一CP40小结展示具体变异
def cqfy_cp40_var_sum_s(var_list):
	result = []
	for var in var_list:
		if var["bio_category"] == "Snvindel":
			if "judge_mergeMET" in var.keys() and var["judge_mergeMET"]:
				if var["hgvs_p"] != "p.?":
					result.append(var["gene_symbol"]+" "+var["hgvs_p"]+"（MET exon14 skipping）")
				else:
					result.append(var["gene_symbol"]+" "+var["hgvs_c"]+"（MET exon14 skipping）")
			else:
				if var["hgvs_p"] != "p.?":
					result.append(var["gene_symbol"]+" "+var["hgvs_p"])
				else:
					result.append(var["gene_symbol"]+" "+var["hgvs_c"])
		elif var["bio_category"] == "Cnv":
			result.append(var["gene_symbol"]+"扩增")
		elif var["bio_category"] in ["Sv", "PSeqRnaSv"]:
			if var["five_prime_gene"] == "MET" and var["three_prime_gene"] == "MET":
				result.append("MET exon14 skipping")
			else:
				# 融合可能会有重复（rna exon相同，断点不同的情况）
				if var["five_prime_gene"]+"-"+var["three_prime_gene"]+"融合" not in result:
					result.append(var["five_prime_gene"]+"-"+var["three_prime_gene"]+"融合")
	return ", ".join(result)
jinja2.filters.FILTERS["cqfy_cp40_var_sum_s"] = cqfy_cp40_var_sum_s

# 重庆附一CP40，变异解读后半段的描述需要做合并-2024.09.03
# 该变异为XX变异，导致XXX。该变异可能与【药物疗效、XX癌预后、XX癌辅助诊断相关）。】
def cqfy_cp40_evi_sum(var):
	evi_type = []
	if "Predictive" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Predictive"]:
		evi_type.append("药物疗效")

	if "Prognostic" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Prognostic"]:
		tumor_lst = []
		for evi in var["evi_sum"]["evi_split"]["Prognostic"]:
			if "tumor_name_cn" in evi.keys() and evi["tumor_name_cn"] not in tumor_lst:
				tumor_lst.append(evi["tumor_name_cn"])
		evi_type.append("、".join(tumor_lst)+"预后")
	
	if "Diagnostic" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Diagnostic"]:
		tumor_lst = []
		for evi in var["evi_sum"]["evi_split"]["Diagnostic"]:
			if "tumor_name_cn" in evi.keys() and evi["tumor_name_cn"] not in tumor_lst:
				tumor_lst.append(evi["tumor_name_cn"])
		evi_type.append("、".join(tumor_lst)+"辅助诊断")

	return "、".join(evi_type)
jinja2.filters.FILTERS["cqfy_cp40_evi_sum"] = cqfy_cp40_evi_sum

# 河北省人民CP40-需要拆分成10基因和其他基因-2024-09-06
# var_type分为lc10和other,lc10_nofound
def hdrm_cp40_split_gene(info):
	var_type = info[0]
	var_list = info[1]
	LC10_gene_list = ["ALK", "BRAF", "EGFR", "ERBB2", "KRAS", "MET", "NRAS", "PIK3CA", "RET", "ROS1"]
	if var_type == "lc10":
		return [var for var in var_list if set(re.split(",", var["gene_symbol"])) & set(LC10_gene_list)]
	elif var_type == "other":
		return [var for var in var_list if not set(re.split(",", var["gene_symbol"])) & set(LC10_gene_list)]
	else:
		detect_gene = []
		for var in var_list:
			for gene in re.split(",", var["gene_symbol"]):
				if gene in LC10_gene_list:
					detect_gene.append(gene)
		return sorted(list(set(LC10_gene_list) - set(detect_gene)))
jinja2.filters.FILTERS["hdrm_cp40_split_gene"] = hdrm_cp40_split_gene

# 复旦中山RNAseq，临床意义提示部分，可能会有多个辅助诊断/临床预后，需要删除重复的-2024.09.09
# 预后：XX癌临床预后较好（A级）；辅助诊断：XX癌辅助诊断标志物（B级）
def fdzs_rnaseq_redup_evi(info):
	evi_type = info[0]
	evi_list = info[1]
	result = []
	for evi in evi_list:
		key = ""
		if evi_type == "临床预后":
			key = "{0}临床预后{1}（{2}级）".format(evi["tumor_name_cn"], evi["clinical_significance_cn"], evi["evi_conclusion_simple"])
		elif evi_type == "辅助诊断":
			key = "{0}辅助诊断标志物（{1}级）".format(evi["tumor_name_cn"], evi["evi_conclusion_simple"])
		if key not in result:
			result.append(key)
	return result
jinja2.filters.FILTERS["fdzs_rnaseq_redup_evi"] = fdzs_rnaseq_redup_evi

# 复旦中山RNAseq，详细解读，辅助诊断/临床预后多条证据时，癌种相同的需要放一起展示-2024.09.09
def fdzs_rnaseq_tumor_group(evi_list):
	tumor_list = []
	evi_dict = {}
	for evi in evi_list:
		if evi["tumor_name_cn"] not in tumor_list:
			tumor_list.append(evi["tumor_name_cn"])
			evi_dict.setdefault(evi["tumor_name_cn"], [])
		evi_dict[evi["tumor_name_cn"]].append(evi["evi_interpretation"])
	
	result = []
	for tumor, inter in evi_dict.items():
		result.append(
			{
				"tumor_name_cn" : tumor,
				"evi_interpretation_list" : inter
			}
		)
	return sorted(result, key=lambda i:tumor_list.index(i["tumor_name_cn"]))
jinja2.filters.FILTERS["fdzs_rnaseq_tumor_group"] = fdzs_rnaseq_tumor_group

# 吉大一LC10-结果小结规则更新-2024.09.11
# 1. 检出含A级证据的变异，则优先展示A级证据：EGFR基因 p.L858R突变，阿法替尼、吉非替尼敏感（耐药分开写），A级推荐请结合临床
# 2. 无A则判断含B级证据的变异，展示B级证据
# 3. 无B则判断含C级证据的变异，展示C级证据
# 4. 无C则判断D级证据变异，展示：检出与靶向药敏感性有关体细胞变异，但药物证据等级有限
def jlyy_lc10_summary_new_rule(info):
	knb = info[0]
	raw_var_list = info[1]
	var_list = []
	if knb:
		var_list.append(knb)
	var_list.extend(raw_var_list)

	# 获取变异证据最高等级
	for var in var_list:
		regimen_level = [i["evi_conclusion_simple"] for i in var["evi_sum"]["regimen_evi_sum"] if i["evidence_type"] in ["Predictive"]] \
					 	if var["evi_sum"] and "regimen_evi_sum" in var["evi_sum"].keys() \
					 	else []				
		var["top_level"] = "A" if "A" in regimen_level else \
						   "B" if "B" in regimen_level else \
						   "C" if "C" in regimen_level else \
						   "D" if "D" in regimen_level else \
						   "N"
	# 将变异列表按最高证据等级A/B/C/D进行拆分
	var_type_dict = {}
	for level in ["A", "B", "C", "D"]:
		var_type_dict["level_"+level+"_var"] = [var for var in var_list if var["top_level"] == level]

	# 处理单个变异的治疗方案
	def get_regimen(evi_list, level):
		# 对evi_list进行排序，按A0/A1/A2这样子排
		evi_list_sort = jlyy_lc10_sort_evi(evi_list)
		sense_regimen = [evi["regimen_name"] for evi in evi_list_sort if evi["clinical_significance_cn"] == "敏感" and evi["evi_conclusion_simple"] == level]
		resis_regimen = [evi["regimen_name"] for evi in evi_list_sort if evi["clinical_significance_cn"] == "耐药" and evi["evi_conclusion_simple"] == level]
		sense_str = "{0}敏感，{1}级推荐请结合临床".format("、".join(sense_regimen), level) if sense_regimen else ""
		resis_str = "{0}耐药，{1}级推荐请结合临床".format("、".join(resis_regimen), level) if resis_regimen else ""
		return sense_str, resis_str
	
	result = []
	# 最高等级A的变异
	if var_type_dict["level_A_var"]:
		for var in var_type_dict["level_A_var"]:
			sense_str, resis_str = get_regimen(var["evi_sum"]["evi_split"]["Predictive"], "A")
			var["jdyy_sense_str"] = sense_str
			var["jdyy_resis_str"] = resis_str
			result.append(var)
	# 最高等级B的变异
	elif var_type_dict["level_B_var"]:
		for var in var_type_dict["level_B_var"]:
			sense_str, resis_str = get_regimen(var["evi_sum"]["evi_split"]["Predictive"], "B")
			var["jdyy_sense_str"] = sense_str
			var["jdyy_resis_str"] = resis_str
			result.append(var)
	# 最高等级C的变异
	elif var_type_dict["level_C_var"]:
		for var in var_type_dict["level_C_var"]:
			sense_str, resis_str = get_regimen(var["evi_sum"]["evi_split"]["Predictive"], "C")
			var["jdyy_sense_str"] = sense_str
			var["jdyy_resis_str"] = resis_str
			result.append(var)
	elif var_type_dict["level_D_var"]:
		result.append(var_type_dict["level_D_var"][0])
	
	return result
jinja2.filters.FILTERS["jlyy_lc10_summary_new_rule"] = jlyy_lc10_summary_new_rule

# 北大人民检验科gHRR，过滤掉BRCA变异-2024.09.18
def bdrmjyk_hrr_filter_brca(var_list):
	return [var for var in var_list if var["gene_symbol"] not in ["BRCA1", "BRCA2"]]
jinja2.filters.FILTERS["bdrmjyk_hrr_filter_brca"] = bdrmjyk_hrr_filter_brca

# 孙逸仙116结果解读，需要把“改变异为XX变异，导致”改成具体基因-2024.09.19
# 内含子/同义突变等：该变异为内含子区突变。（改为“XXX基因内含子区突变。”）
# 其他突变：该变异为错义突变，导致基因编码蛋白……（改为“XXX基因编码蛋白……”）
# 其他情况就照原始的返回好了。
def syx_116_inter(info):
	variant_desc_cn = info[0]
	gene_symbol = info[1]
	var_inter = re.split("，", variant_desc_cn)
	if len(var_inter) == 1:
		return variant_desc_cn.replace("该变异为", gene_symbol+"基因")
	elif len(var_inter) == 2:
		return var_inter[1].replace("导致", gene_symbol)
	else:
		return variant_desc_cn
jinja2.filters.FILTERS["syx_116_inter"] = syx_116_inter

# 吉林大学第一院LC10-新增10个基因的汇总表-2024.09.19
def jlyy_10gene_result(var_list):
	# 提示敏感/耐药
	def sense_resist(gene_var_list, sense_or_resist):
		result = []
		copy_gene_var_list = copy.deepcopy(gene_var_list)
		for var in copy_gene_var_list:
			if "Predictive" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Predictive"]:
				regimen_list = ["{0}（{1}）".format(evi["regimen_name"], evi["evi_conclusion"][::-1]) for \
								evi in jlyy_lc10_sort_evi(var["evi_sum"]["evi_split"]["Predictive"]) if evi["evi_conclusion_simple"] != "D" and \
								evi["clinical_significance_cn"] == sense_or_resist]
				if regimen_list:
					var["regimen_sum_for_jdyy"] = "，".join(regimen_list)
					var["sense_or_resist"] = sense_or_resist
					result.append(var)
		return result
	gene_list = ["ALK", "BRAF", "EGFR", "ERBB2", "KRAS", "MET", "NRAS", "PIK3CA", "RET", "ROS1"]
	all_gene_result = []
	all_gene_result.append({
		"gene_symbol" : "基因",
		"result" : "title",
		"sort_rule" : 0
	})
	for gene in gene_list:
		gene_var_list = [var for var in var_list if set(re.split(",", var["gene_symbol"])) & set([gene])]
		single_gene_result = []
		if sense_resist(gene_var_list, "敏感"):
			single_gene_result.extend(sense_resist(gene_var_list, "敏感"))
		if sense_resist(gene_var_list, "耐药"):
			single_gene_result.extend(sense_resist(gene_var_list, "耐药"))
		if gene_var_list:
			if single_gene_result:
				all_gene_result.append({
					"gene_symbol" : gene,
					"var_result" : single_gene_result,
					"sort_rule" : 1
				})
			else:
				all_gene_result.append({
					"gene_symbol" : gene,
					"result" : "nofoundsignificant",
					"sort_rule" : 2
				})
		else:
			all_gene_result.append({
				"gene_symbol" : gene,
				"result" : "nofoundvar",
				"sort_rule" : 3
			})
	#print (sorted(all_gene_result, key=lambda i:i["sort_rule"]))
	return sorted(all_gene_result, key=lambda i:i["sort_rule"])
jinja2.filters.FILTERS["jlyy_10gene_result"] = jlyy_10gene_result

# BPTM plus配对展示变异小结-2024.09.20
# 包含snvindel和MLPA变异
def pt_bptm_plus_var_sum(var_list):
	result = []
	for var in var_list:
		if var["type"] == "Loss" and var["gene_symbol"] in ["BRCA1", "BRCA2"]:
			result.append(var["gene_symbol"] + " " + var["value"] + " del")
		elif var["type"] == "Gain" and var["gene_symbol"] in ["BRCA1", "BRCA2"]:
			result.append(var["gene_symbol"] + " " + var["value"] + " dup")
		else:
			if var["hgvs_p"] != "p.?":
				result.append(var["gene_symbol"]+" "+var["hgvs_c"]+" "+var["hgvs_p"])
			else:
				result.append(var["gene_symbol"]+" "+var["hgvs_c"])
	return ", ".join(result)
jinja2.filters.FILTERS["pt_bptm_plus_var_sum"] = pt_bptm_plus_var_sum

# 判断变异是否有治疗方案ABC等级证据+预后/辅助诊断ABCD等级证据-2024.09.23
def judge_abc_evi_regimen(var):
	regimen_level = [i["evi_conclusion_simple"] for i in var["evi_sum"]["regimen_evi_sum"] if i["evidence_type"] in ["Predictive"]] \
					 if var["evi_sum"] and "regimen_evi_sum" in var["evi_sum"].keys() \
					 else []
	other_level = [i["evi_conclusion_simple"] for i in var["evi_sum"]["regimen_evi_sum"] if i["evidence_type"] in ["Prognostic", "Diagnostic"]] \
					 if var["evi_sum"] and "regimen_evi_sum" in var["evi_sum"].keys() \
					 else [] 
	if set(["A", "B", "C"]) & set(regimen_level) or set(["A", "B", "C", "D"]) & set(other_level):
		return True
	else:
		return False
jinja2.filters.FILTERS["judge_abc_evi_regimen"] = judge_abc_evi_regimen

# 云南肿瘤BPMT-检出BRCA变异需要在小结备注处添加提示-2024.09.24
def ynzl_brca_note(var_list):
	# BRCA1基因检测到X个致病性变异和X个疑似致病性变异，BRCA2基因未检测到致病性或疑似致病性变异。
	result = []
	brca1_5_count = len([var for var in var_list if var["gene_symbol"] == "BRCA1" and var["clinic_num_g"] == 5])
	brca1_4_count = len([var for var in var_list if var["gene_symbol"] == "BRCA1" and var["clinic_num_g"] == 4])
	brca2_5_count = len([var for var in var_list if var["gene_symbol"] == "BRCA2" and var["clinic_num_g"] == 5])
	brca2_4_count = len([var for var in var_list if var["gene_symbol"] == "BRCA2" and var["clinic_num_g"] == 4])
	if brca1_5_count and brca1_4_count:
		result.append("BRCA1基因检测到{0}个致病性变异和{1}个疑似致病性变异".format(str(brca1_5_count), str(brca1_4_count)))
	elif brca1_5_count and not brca1_4_count:
		result.append("BRCA1基因检测到{0}个致病性变异".format(str(brca1_5_count)))
	elif not brca1_5_count and brca1_4_count:
		result.append("BRCA1基因检测到{0}个疑似致病性变异".format(str(brca1_4_count)))

	if brca2_5_count and brca2_4_count:
		result.append("BRCA2基因检测到{0}个致病性变异和{1}个疑似致病性变异".format(str(brca2_5_count), str(brca2_4_count)))
	elif brca2_5_count and not brca2_4_count:
		result.append("BRCA2基因检测到{0}个致病性变异".format(str(brca2_5_count)))
	elif not brca2_5_count and brca2_4_count:
		result.append("BRCA2基因检测到{0}个疑似致病性变异".format(str(brca2_4_count)))

	for gene in set(["BRCA1", "BRCA2"]) - set([var["gene_symbol"] for var in var_list]):
		result.append("{0}基因未检测到致病性或疑似致病性变异".format(gene))

	return "，".join(result)
jinja2.filters.FILTERS["ynzl_brca_note"] = ynzl_brca_note

# ptBPTM Plus-子宫内膜癌分子分型表格-需要结合体细胞和胚系变异-2024.09.24
# 体细胞I/II + 胚系4/5,POLE需要限定变异分组，但是系统没有返回，就都放吧
def ptbptm_plus_ec_sum(info):
	var_list = info[0]
	msi = info[1]
	result = []
	ec_gene = ["POLE", "TP53"]
	pole_var = [var for var in var_list if var["gene_symbol"] == "POLE"]
	if pole_var:
		result.extend(pole_var)
	else:
		result.append({"gene_symbol" : "POLE", "result" : "nofound"})
	result.append(msi)
	tp53_var = [var for var in var_list if var["gene_symbol"] == "TP53"]
	if tp53_var:
		result.extend(tp53_var)
	else:
		result.append({"gene_symbol" : "TP53", "result" : "nofound"})

	return result
jinja2.filters.FILTERS["ptbptm_plus_ec_sum"] = ptbptm_plus_ec_sum

# 吉林大学第一院LC10-返回单个基因结果-2024.09.26
def jlyy_10gene_result_singlegene(info):
	var_list = info[0]
	gene = info[1]
	sense_or_resist = info[2]
	gene_var_list = [var for var in var_list if set(re.split(",", var["gene_symbol"])) & set([gene])]
	result = []
	for var in gene_var_list:
		if "Predictive" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Predictive"]:
			regimen_list = ["{0}（{1}）".format(evi["regimen_name"], evi["evi_conclusion"][::-1]) for \
								evi in jlyy_lc10_sort_evi(var["evi_sum"]["evi_split"]["Predictive"]) if evi["evi_conclusion_simple"] != "D" and \
								evi["clinical_significance_cn"] == sense_or_resist]
			if regimen_list:
				var["regimen_sum_for_jdyy"] = "，".join(regimen_list)
				result.append(var)
	if gene_var_list and not result:
		result = [{"gene_symbol" : "nofoundsignificant"}]
	if not gene_var_list:
		result = [{"gene_symbol" : "nofoundvar"}]
	return result
jinja2.filters.FILTERS["jlyy_10gene_result_singlegene"] = jlyy_10gene_result_singlegene

# 吉林大学第一院LC10-返回单个基因结果-V2-2024.11.04
def jlyy_10gene_result_singlegene_v2(info):
	# 与上一版相比区别在
	# 1、删除证据等级的数字，排序按通用的来
	# 2、只提示最高级别药物（D级不要）
	var_list = info[0]
	gene = info[1]
	sense_or_resist = info[2]
	gene_var_list = [var for var in var_list if set(re.split(",", var["gene_symbol"])) & set([gene])]
	#print (gene_var_list)
	result = []
	D_sense_regimen = []
	for var in gene_var_list:
		#print (var["gene_symbol"])
		if "Predictive" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Predictive"]:
			regimen_list = []
			A_regimen = ["{0}（{1}）".format(evi["regimen_name"], evi["evi_conclusion_simple"]) for evi in var["evi_sum"]["evi_split"]["Predictive"] if \
						evi["evi_conclusion_simple"] == "A" and evi["clinical_significance_cn"] == sense_or_resist]
			B_regimen = ["{0}（{1}）".format(evi["regimen_name"], evi["evi_conclusion_simple"]) for evi in var["evi_sum"]["evi_split"]["Predictive"] if \
						evi["evi_conclusion_simple"] == "B" and evi["clinical_significance_cn"] == sense_or_resist]
			C_regimen = ["{0}（{1}）".format(evi["regimen_name"], evi["evi_conclusion_simple"]) for evi in var["evi_sum"]["evi_split"]["Predictive"] if \
						evi["evi_conclusion_simple"] == "C" and evi["clinical_significance_cn"] == sense_or_resist]
			D_regimen = ["{0}（{1}）".format(evi["regimen_name"], evi["evi_conclusion_simple"]) for evi in var["evi_sum"]["evi_split"]["Predictive"] if \
						evi["evi_conclusion_simple"] == "D" and evi["clinical_significance_cn"] == sense_or_resist]
			regimen_list = A_regimen if A_regimen else \
						   B_regimen if B_regimen else \
						   C_regimen if C_regimen else D_sense_regimen.extend(D_regimen)
			#print ("A, ", len(A_regimen))
			#print ("B, ", len(B_regimen))
			#print ("C, ", len(C_regimen))
			#print ("D, ", len(D_regimen))
			if regimen_list:
				var["regimen_sum_for_jdyy"] = "，".join(regimen_list)
				result.append(var)
	if gene_var_list and not result and D_sense_regimen:
		result = [{"gene_symbol" : "nofoundsignificant"}]
	if not result:
		result = [{"gene_symbol" : "nofoundvar"}]
	#print (result)
	return result
jinja2.filters.FILTERS["jlyy_10gene_result_singlegene_v2"] = jlyy_10gene_result_singlegene_v2

# 吉林大学第一院LC10-返回耐药结果-2024.09.26
def jlyy_10gene_result_resis_result(var_list):
	result = []
	regimen_D_list = []
	for var in var_list:
		if "Predictive" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Predictive"]:
			regimen_list = ["{0}（{1}）".format(evi["regimen_name"], evi["evi_conclusion"][::-1]) for \
								evi in jlyy_lc10_sort_evi(var["evi_sum"]["evi_split"]["Predictive"]) if evi["evi_conclusion_simple"] != "D" and \
								evi["clinical_significance_cn"] == "耐药"]
			
			regimen_D = ["{0}（{1}）".format(evi["regimen_name"], evi["evi_conclusion"][::-1]) for \
								evi in jlyy_lc10_sort_evi(var["evi_sum"]["evi_split"]["Predictive"]) if evi["evi_conclusion_simple"] == "D" and \
								evi["clinical_significance_cn"] == "耐药"]
			if regimen_list:
				var["regimen_sum_for_jdyy"] = "，".join(regimen_list)
				result.append(var)
			if regimen_D:
				regimen_D_list.append(regimen_D)
	if result:
		result = sorted(result, key = lambda i:i["gene_symbol"])
	else:
		if regimen_D_list:
			result = [{"gene_symbol" : "nofoundsignificant"}]
		else:
			result = [{"gene_symbol" : "nofoundvar"}]
	return result
jinja2.filters.FILTERS["jlyy_10gene_result_resis_result"] = jlyy_10gene_result_resis_result

# 吉林大学第一院LC10-返回耐药结果-V2-2024.11.04
def jlyy_10gene_result_resis_result_v2(var_list):
	# 与上一版相比区别在
	# 1、删除证据等级的数字，排序按通用的来
	# 2、只提示最高级别药物（D级也要）
	result = []
	for var in var_list:
		if "Predictive" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Predictive"]:
			A_regimen = ["{0}（{1}）".format(evi["regimen_name"], evi["evi_conclusion_simple"]) for evi in var["evi_sum"]["evi_split"]["Predictive"] if \
						evi["evi_conclusion_simple"] == "A" and evi["clinical_significance_cn"] == "耐药"]
			B_regimen = ["{0}（{1}）".format(evi["regimen_name"], evi["evi_conclusion_simple"]) for evi in var["evi_sum"]["evi_split"]["Predictive"] if \
						evi["evi_conclusion_simple"] == "B" and evi["clinical_significance_cn"] == "耐药"]
			C_regimen = ["{0}（{1}）".format(evi["regimen_name"], evi["evi_conclusion_simple"]) for evi in var["evi_sum"]["evi_split"]["Predictive"] if \
						evi["evi_conclusion_simple"] == "C" and evi["clinical_significance_cn"] == "耐药"]
			D_regimen = ["{0}（{1}）".format(evi["regimen_name"], evi["evi_conclusion_simple"]) for evi in var["evi_sum"]["evi_split"]["Predictive"] if \
						evi["evi_conclusion_simple"] == "D" and evi["clinical_significance_cn"] == "耐药"]
			regimen_list = A_regimen if A_regimen else \
						   B_regimen if B_regimen else \
						   C_regimen if C_regimen else \
						   D_regimen if D_regimen else []
			if regimen_list:
				var["regimen_sum_for_jdyy"] = "，".join(regimen_list)
				result.append(var)
	if result:
		result = sorted(result, key = lambda i:i["gene_symbol"])
	else:
		result = [{"gene_symbol" : "nofoundvar"}]
	return result
jinja2.filters.FILTERS["jlyy_10gene_result_resis_result_v2"] = jlyy_10gene_result_resis_result_v2

# 吉林大学第一院LC10-返回耐药结果-V3-2025.05.08
def jlyy_10gene_result_resis_result_v3(var_list):
	# 与上一版相比区别在
	# 只提示最高级别药物（仅有D的不展示，并且模板里输出“检出耐药相关基因，但证据有限”）
	result = []
	D_resis_regimen = []
	for var in var_list:
		if "Predictive" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Predictive"]:
			A_regimen = ["{0}（{1}）".format(evi["regimen_name"], evi["evi_conclusion_simple"]) for evi in var["evi_sum"]["evi_split"]["Predictive"] if \
						evi["evi_conclusion_simple"] == "A" and evi["clinical_significance_cn"] == "耐药"]
			B_regimen = ["{0}（{1}）".format(evi["regimen_name"], evi["evi_conclusion_simple"]) for evi in var["evi_sum"]["evi_split"]["Predictive"] if \
						evi["evi_conclusion_simple"] == "B" and evi["clinical_significance_cn"] == "耐药"]
			C_regimen = ["{0}（{1}）".format(evi["regimen_name"], evi["evi_conclusion_simple"]) for evi in var["evi_sum"]["evi_split"]["Predictive"] if \
						evi["evi_conclusion_simple"] == "C" and evi["clinical_significance_cn"] == "耐药"]
			D_regimen = ["{0}（{1}）".format(evi["regimen_name"], evi["evi_conclusion_simple"]) for evi in var["evi_sum"]["evi_split"]["Predictive"] if \
						evi["evi_conclusion_simple"] == "D" and evi["clinical_significance_cn"] == "耐药"]
			regimen_list = A_regimen if A_regimen else \
						   B_regimen if B_regimen else \
						   C_regimen if C_regimen else D_resis_regimen.extend(D_regimen)
			if regimen_list:
				var["regimen_sum_for_jdyy"] = "，".join(regimen_list)
				result.append(var)
	if result:
		result = sorted(result, key = lambda i:i["gene_symbol"])
	else:
		if D_resis_regimen:
			result = [{"gene_symbol" : "nofoundsignificant"}]
		else:
			result = [{"gene_symbol" : "nofoundvar"}]
	return result
jinja2.filters.FILTERS["jlyy_10gene_result_resis_result_v3"] = jlyy_10gene_result_resis_result_v3

# 吉林大学第一院LC10-返回敏感结果-2024.11.04
def jlyy_10gene_result_sense_result(var_list):
	# 与上一版相比区别在
	# 1、删除证据等级的数字，排序按通用的来
	# 2、只提示最高级别药物（D级也要）
	result = []
	#D_sense_regimen = []
	for var in var_list:
		if "Predictive" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Predictive"]:
			regimen_list = []
			A_regimen = ["{0}（{1}）".format(evi["regimen_name"], evi["evi_conclusion_simple"]) for evi in var["evi_sum"]["evi_split"]["Predictive"] if \
						evi["evi_conclusion_simple"] == "A" and evi["clinical_significance_cn"] == "敏感"]
			B_regimen = ["{0}（{1}）".format(evi["regimen_name"], evi["evi_conclusion_simple"]) for evi in var["evi_sum"]["evi_split"]["Predictive"] if \
						evi["evi_conclusion_simple"] == "B" and evi["clinical_significance_cn"] == "敏感"]
			C_regimen = ["{0}（{1}）".format(evi["regimen_name"], evi["evi_conclusion_simple"]) for evi in var["evi_sum"]["evi_split"]["Predictive"] if \
						evi["evi_conclusion_simple"] == "C" and evi["clinical_significance_cn"] == "敏感"]
			D_regimen = ["{0}（{1}）".format(evi["regimen_name"], evi["evi_conclusion_simple"]) for evi in var["evi_sum"]["evi_split"]["Predictive"] if \
						evi["evi_conclusion_simple"] == "D" and evi["clinical_significance_cn"] == "敏感"]
			regimen_list = A_regimen if A_regimen else \
						   B_regimen if B_regimen else \
						   C_regimen if C_regimen else []
			if regimen_list:
				var["regimen_sum_for_jdyy"] = "，".join(regimen_list)
				result.append(var)
	if result:
		result = sorted(result, key = lambda i:i["gene_symbol"])
	return result
jinja2.filters.FILTERS["jlyy_10gene_result_sense_result"] = jlyy_10gene_result_sense_result


# 吉大一LC10-结果小结规则更新-按基因分类-2024.09.26
# 1. 检出含A级证据的变异，则优先展示A级证据：EGFR基因 p.L858R突变，阿法替尼、吉非替尼敏感（耐药分开写），A级推荐请结合临床
# 2. 无A则判断含B级证据的变异，展示B级证据
# 3. 无B则判断含C级证据的变异，展示C级证据
# 4. 无C则判断D级证据变异，展示：检出与靶向药敏感性有关体细胞变异，但药物证据等级有限
def jlyy_lc10_summary_new_rule_single_gene(var_list):
	#knb = info[0]
	#raw_var_list = info[1]
	#var_list = []
	#if knb:
	#	var_list.append(knb)
	#var_list.extend(raw_var_list)

	# 获取变异证据最高等级
	for var in var_list:
		regimen_level = [i["evi_conclusion_simple"] for i in var["evi_sum"]["regimen_evi_sum"] if i["evidence_type"] in ["Predictive"]] \
					 	if var["evi_sum"] and "regimen_evi_sum" in var["evi_sum"].keys() \
					 	else []				
		var["top_level"] = "A" if "A" in regimen_level else \
						   "B" if "B" in regimen_level else \
						   "C" if "C" in regimen_level else \
						   "D" if "D" in regimen_level else \
						   "N"
	# 将变异列表按最高证据等级A/B/C/D进行拆分
	var_type_dict = {}
	for level in ["A", "B", "C", "D"]:
		var_type_dict["level_"+level+"_var"] = [var for var in var_list if var["top_level"] == level]

	# 处理单个变异的治疗方案
	def get_regimen(evi_list, level):
		# 对evi_list进行排序，按A0/A1/A2这样子排
		evi_list_sort = jlyy_lc10_sort_evi(evi_list)
		sense_regimen = [evi["regimen_name"] for evi in evi_list_sort if evi["clinical_significance_cn"] == "敏感" and evi["evi_conclusion_simple"] == level]
		resis_regimen = [evi["regimen_name"] for evi in evi_list_sort if evi["clinical_significance_cn"] == "耐药" and evi["evi_conclusion_simple"] == level]
		sense_str = "{0}敏感，{1}级推荐请结合临床".format("、".join(sense_regimen), level) if sense_regimen else ""
		resis_str = "{0}耐药，{1}级推荐请结合临床".format("、".join(resis_regimen), level) if resis_regimen else ""
		return sense_str, resis_str
	
	result = []
	# 最高等级A的变异
	if var_type_dict["level_A_var"]:
		for var in var_type_dict["level_A_var"]:
			sense_str, resis_str = get_regimen(var["evi_sum"]["evi_split"]["Predictive"], "A")
			var["jdyy_sense_str"] = sense_str
			var["jdyy_resis_str"] = resis_str
			result.append(var)
	# 最高等级B的变异
	elif var_type_dict["level_B_var"]:
		for var in var_type_dict["level_B_var"]:
			sense_str, resis_str = get_regimen(var["evi_sum"]["evi_split"]["Predictive"], "B")
			var["jdyy_sense_str"] = sense_str
			var["jdyy_resis_str"] = resis_str
			result.append(var)
	# 最高等级C的变异
	elif var_type_dict["level_C_var"]:
		for var in var_type_dict["level_C_var"]:
			sense_str, resis_str = get_regimen(var["evi_sum"]["evi_split"]["Predictive"], "C")
			var["jdyy_sense_str"] = sense_str
			var["jdyy_resis_str"] = resis_str
			result.append(var)
	#elif var_type_dict["level_D_var"]:
	#	result.append(var_type_dict["level_D_var"][0])
	
	return result

def jlyy_lc10_summary_new_rule_summary(info):
	knb = info[0]
	raw_var_list = info[1]
	result = []
	if knb:
		result.extend(jlyy_lc10_summary_new_rule_single_gene([knb]))
	gene_list = ["ALK", "BRAF", "EGFR", "ERBB2", "KRAS", "MET", "NRAS", "PIK3CA", "RET", "ROS1"]
	for gene in gene_list:
		var_list  = [var for var in raw_var_list if set(re.split(",", var["gene_symbol"])) & set([gene])]
		gene_regimen = jlyy_lc10_summary_new_rule_single_gene(var_list)
		if gene_regimen:
			result.extend(gene_regimen)
	return result
jinja2.filters.FILTERS["jlyy_lc10_summary_new_rule_summary"] = jlyy_lc10_summary_new_rule_summary

def jlyy_lc10_summary_top_level(info):
	knb = info[0]
	raw_var_list = info[1]
	var_list = []
	if knb:
		var_list.append(knb)
	var_list.extend(raw_var_list)

	# 获取变异证据最高等级
	top_level = ""
	for var in var_list:
		regimen_level = [i["evi_conclusion_simple"] for i in var["evi_sum"]["regimen_evi_sum"] if i["evidence_type"] in ["Predictive"]] \
					 	if var["evi_sum"] and "regimen_evi_sum" in var["evi_sum"].keys() \
					 	else []				
		top_level = "A" if "A" in regimen_level else \
						   "B" if "B" in regimen_level else \
						   "C" if "C" in regimen_level else \
						   "D" if "D" in regimen_level else \
						   "N"
	return top_level
jinja2.filters.FILTERS["jlyy_lc10_summary_top_level"] = jlyy_lc10_summary_top_level

# 删除KNB-2024.11.05
def jlyy_lc10_summary_top_level_v2(var_list):
	# 获取变异证据最高等级
	top_level = ""
	level_sum = []
	for var in var_list:
		regimen_level = [i["evi_conclusion_simple"] for i in var["evi_sum"]["regimen_evi_sum"] if i["evidence_type"] in ["Predictive"] and i["clinical_significance_cn"] == "敏感"] \
					 	if var["evi_sum"] and "regimen_evi_sum" in var["evi_sum"].keys() else []
		level_sum.extend(regimen_level)				
	top_level = "A" if "A" in level_sum else \
				"B" if "B" in level_sum else \
				"C" if "C" in level_sum else \
				"D" if "D" in level_sum else \
				"N"
	return top_level
jinja2.filters.FILTERS["jlyy_lc10_summary_top_level_v2"] = jlyy_lc10_summary_top_level_v2

# 吉林大学第一医院116-治疗方案重新排序下，根据敏感>耐药，A1>A2>A3……排序-2024.09.27
def jlyy_116_sort_evi(evi_sum):
	return sorted(evi_sum, key=lambda i:(i["sense_rule"], i["evi_conclusion"], i["regimen_name_py"].upper()))
jinja2.filters.FILTERS["jlyy_116_sort_evi"] = jlyy_116_sort_evi

# 吉林大学第一院116-返回敏感/耐药结果-2024.09.27
def jlyy_116gene_result(info):
	knb = info[0]
	raw_var_list = info[1]
	var_list = []
	if knb:
		if type(knb).__name__=="dict":
			var_list.extend([knb])
		else:
			var_list.extend(knb)
	var_list.extend(raw_var_list)
	sense_or_resis = info[2]

	result = []
	regimen_D_list = []
	for var in var_list:
		if "Predictive" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Predictive"]:
			regimen_list = ["{0}（{1}）".format(evi["regimen_name"], evi["evi_conclusion"][::-1]) for \
								evi in jlyy_lc10_sort_evi(var["evi_sum"]["evi_split"]["Predictive"]) if evi["evi_conclusion_simple"] != "D" and \
								evi["clinical_significance_cn"] == sense_or_resis]
			
			regimen_D = ["{0}（{1}）".format(evi["regimen_name"], evi["evi_conclusion"][::-1]) for \
								evi in jlyy_lc10_sort_evi(var["evi_sum"]["evi_split"]["Predictive"]) if evi["evi_conclusion_simple"] == "D" and \
								evi["clinical_significance_cn"] == sense_or_resis]
			if regimen_list:
				var["regimen_sum_for_jdyy"] = "，".join(regimen_list)
				result.append(var)
			if regimen_D:
				regimen_D_list.append(regimen_D)
	if not result:
		if regimen_D_list:
			result = [{"gene_symbol" : "nofoundsignificant"}]
		else:
			result = [{"gene_symbol" : "nofoundvar"}]
	return result
jinja2.filters.FILTERS["jlyy_116gene_result"] = jlyy_116gene_result

# 吉林大学第一院116-返回敏感结果-2024.11.05
# 敏感药物只展示等级最高的（除了C）
def jlyy_116gene_result_sense(info):
	knb = info[0]
	raw_var_list = info[1]
	var_list = []
	if knb:
		if type(knb).__name__=="dict":
			var_list.extend([knb])
		else:
			var_list.extend(knb)
	var_list.extend(raw_var_list)
	result_sense = []
	regimen_D_list_sense = []
	for var in var_list:
		if "Predictive" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Predictive"]:
			tmp_dict = {}
			regimen_sense_list = []
			for level in ["A", "B", "C", "D"]:
				tmp_dict["sense_"+level] = ["{0}（{1}）".format(evi["regimen_name"], evi["evi_conclusion_simple"]) for \
										   evi in var["evi_sum"]["evi_split"]["Predictive"] if evi["evi_conclusion_simple"] == level and \
										   evi["clinical_significance_cn"] == "敏感"]
			regimen_sense_list = tmp_dict["sense_A"] if tmp_dict["sense_A"] else \
								 tmp_dict["sense_B"] if tmp_dict["sense_B"] else \
								 tmp_dict["sense_C"] if tmp_dict["sense_C"] else \
								 regimen_D_list_sense.extend(tmp_dict["sense_D"])
		
			if regimen_sense_list:
				var["regimen_sum_for_jdyy"] = "，".join(regimen_sense_list)
				result_sense.append(var)

	if not result_sense:
		if regimen_D_list_sense:
			result_sense = [{"gene_symbol" : "nofoundsignificant"}]
		else:
			result_sense = [{"gene_symbol" : "nofoundvar"}]
	return result_sense
jinja2.filters.FILTERS["jlyy_116gene_result_sense"] = jlyy_116gene_result_sense

# 吉林大学第一院116-返回敏感结果-肺癌耐抵扣-2025.05.15
# 敏感药物只展示等级最高的（除了C）
def jlyy_116gene_result_sense_filter_egfr(info):
	knb = info[0]
	var_list = info[1]
	tumor_list = info[2]
	# 判断是否检出MET扩增
	judge_met_cnv = "F"
	for var in var_list:
		if var["bio_category"] == "Cnv" and var["gene_symbol"] == "MET":
			judge_met_cnv = "T"
			break
	result_sense = []
	regimen_D_list_sense = []
	for var in var_list:
		if "Predictive" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Predictive"]:
			if not (var["bio_category"] == "Snvindel" and var["gene_symbol"] == "EGFR" and judge_met_cnv == "T" and "肺癌" in tumor_list):
				tmp_dict = {}
				regimen_sense_list = []
				for level in ["A", "B", "C", "D"]:
					tmp_dict["sense_"+level] = ["{0}（{1}）".format(evi["regimen_name"], evi["evi_conclusion_simple"]) for \
											   evi in var["evi_sum"]["evi_split"]["Predictive"] if evi["evi_conclusion_simple"] == level and \
											   evi["clinical_significance_cn"] == "敏感"]
				regimen_sense_list = tmp_dict["sense_A"] if tmp_dict["sense_A"] else \
									 tmp_dict["sense_B"] if tmp_dict["sense_B"] else \
									 tmp_dict["sense_C"] if tmp_dict["sense_C"] else \
									 regimen_D_list_sense.extend(tmp_dict["sense_D"])
		
				if regimen_sense_list:
					var["regimen_sum_for_jdyy"] = "，".join(regimen_sense_list)
					result_sense.append(var)
	if knb:
		knb["regimen_sum_for_jdyy"] = "西妥昔单抗，帕尼单抗，西妥昔单抗β"
		result_sense.insert(0, knb)

	if not result_sense:
		if regimen_D_list_sense:
			result_sense = [{"gene_symbol" : "nofoundsignificant"}]
		else:
			result_sense = [{"gene_symbol" : "nofoundvar"}]
	return result_sense
jinja2.filters.FILTERS["jlyy_116gene_result_sense_filter_egfr"] = jlyy_116gene_result_sense_filter_egfr

def jdy_116_evi_sum(Predictive, filter_regimen):
	tmp_dict = {}
	regimen_sense_list = []
	for level in ["A", "B", "C", "D"]:
		tmp_dict["sense_"+level] = ["{0}（{1}）".format(evi["regimen_name"], evi["evi_conclusion_simple"]) for \
								   evi in Predictive if evi["evi_conclusion_simple"] == level and \
								   evi["clinical_significance_cn"] == "敏感" and evi["regimen_name"] not in filter_regimen]
	regimen_sense_list = tmp_dict["sense_A"] if tmp_dict["sense_A"] else \
						 tmp_dict["sense_B"] if tmp_dict["sense_B"] else \
						 tmp_dict["sense_C"] if tmp_dict["sense_C"] else []
	return regimen_sense_list

# 吉林大学第一院116-返回敏感结果-肺癌耐抵扣-2025.06.17
# 敏感药物只展示等级最高的（除了C）
# 2025.06.17新增规则：存在EGFR T790M，EGFR 部分变异不展示四个药物
def jlyy_116gene_result_sense_filter_egfr_v2(info):
	knb = info[0]
	var_list = info[1]
	tumor_list = info[2]
	# 判断是否检出MET扩增
	judge_met_cnv = "F"
	for var in var_list:
		if var["bio_category"] == "Cnv" and var["gene_symbol"] == "MET":
			judge_met_cnv = "T"
			break
	# 判断是否检出EGFR T790M
	judge_egfr_t790m = "F"
	for var in var_list:
		if var["bio_category"] == "Snvindel" and var["gene_symbol"] == "EGFR" and var["hgvs_p"] == "p.(T790M)":
			judge_egfr_t790m = "T"
			break
	result_sense = []
	regimen_D_list_sense = []
	for var in var_list:
		if "Predictive" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Predictive"]:
			regimen_sense_list = []
			if "肺癌" in tumor_list:
				if var["bio_category"] == "Snvindel" and var["gene_symbol"] == "EGFR":
					judge_egfr_sense_var = "F"
					if var["hgvs_p"] in ["p.(L858R)", "p.(L861Q)", "p.(S768I)"]:
						judge_egfr_sense_var = "T"
					elif "G719" in var["hgvs_p"]:
						judge_egfr_sense_var = "T"
					elif var["gene_region"] == "exon19" and "del" in var["hgvs_p"]:
						judge_egfr_sense_var = "T"
					
					if judge_met_cnv == "T":
						continue
					elif judge_egfr_t790m == "T" and judge_egfr_sense_var == "T":
						regimen_sense_list = jdy_116_evi_sum(var["evi_sum"]["evi_split"]["Predictive"], ["阿法替尼", "达可替尼", "厄洛替尼", "吉非替尼"])
					else:
						regimen_sense_list = jdy_116_evi_sum(var["evi_sum"]["evi_split"]["Predictive"], [])						
				else:
					regimen_sense_list = jdy_116_evi_sum(var["evi_sum"]["evi_split"]["Predictive"], [])	
			else:
				regimen_sense_list = jdy_116_evi_sum(var["evi_sum"]["evi_split"]["Predictive"], [])	

			if regimen_sense_list:
				var["regimen_sum_for_jdyy"] = "，".join(regimen_sense_list)
				result_sense.append(var)
	if knb:
		knb["regimen_sum_for_jdyy"] = "西妥昔单抗，帕尼单抗，西妥昔单抗β"
		result_sense.insert(0, knb)

	if not result_sense:
		if regimen_D_list_sense:
			result_sense = [{"gene_symbol" : "nofoundsignificant"}]
		else:
			result_sense = [{"gene_symbol" : "nofoundvar"}]
	return result_sense
jinja2.filters.FILTERS["jlyy_116gene_result_sense_filter_egfr_v2"] = jlyy_116gene_result_sense_filter_egfr_v2

# 吉林大学第一院116-返回耐药结果-2024.11.05
# 耐药药物展示等级最高的（仅有D的话D级证据也要展示）
def jlyy_116gene_result_resis(info):
	knb = info[0]
	raw_var_list = info[1]
	var_list = []
	if knb:
		if type(knb).__name__=="dict":
			var_list.extend([knb])
		else:
			var_list.extend(knb)
	var_list.extend(raw_var_list)
	result_resis = []
	for var in var_list:
		if "Predictive" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Predictive"]:
			tmp_dict = {}
			for level in ["A", "B", "C", "D"]:
				tmp_dict["resis_"+level] = ["{0}（{1}）".format(evi["regimen_name"], evi["evi_conclusion_simple"]) for \
								  		  evi in var["evi_sum"]["evi_split"]["Predictive"] if evi["evi_conclusion_simple"] == level and \
										  evi["clinical_significance_cn"] == "耐药"]
			regimen_resis_list = tmp_dict["resis_A"] if tmp_dict["resis_A"] else \
								 tmp_dict["resis_B"] if tmp_dict["resis_B"] else \
								 tmp_dict["resis_C"] if tmp_dict["resis_C"] else \
								 tmp_dict["resis_D"] if tmp_dict["resis_D"] else []
			if regimen_resis_list:
				var["regimen_sum_for_jdyy"] = "，".join(regimen_resis_list)
				result_resis.append(var)
	if not result_resis:
		result_resis = [{"gene_symbol" : "nofoundvar"}]
	return result_resis
jinja2.filters.FILTERS["jlyy_116gene_result_resis"] = jlyy_116gene_result_resis

# 吉林大学第一院116-返回耐药结果-V2-2025.05.15
# 耐药药物展示等级最高的（仅有D的话不报）
def jlyy_116gene_result_resis_v2(info):
	knb = info[0]
	raw_var_list = info[1]
	var_list = []
	if knb:
		if type(knb).__name__=="dict":
			var_list.extend([knb])
		else:
			var_list.extend(knb)
	var_list.extend(raw_var_list)
	result_resis = []
	regimen_D_list_resis = []
	for var in var_list:
		if "Predictive" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Predictive"]:
			tmp_dict = {}
			for level in ["A", "B", "C", "D"]:
				tmp_dict["resis_"+level] = ["{0}（{1}）".format(evi["regimen_name"], evi["evi_conclusion_simple"]) for \
								  		  evi in var["evi_sum"]["evi_split"]["Predictive"] if evi["evi_conclusion_simple"] == level and \
										  evi["clinical_significance_cn"] == "耐药"]
			regimen_resis_list = tmp_dict["resis_A"] if tmp_dict["resis_A"] else \
								 tmp_dict["resis_B"] if tmp_dict["resis_B"] else \
								 tmp_dict["resis_C"] if tmp_dict["resis_C"] else \
								 regimen_D_list_resis.extend(tmp_dict["resis_D"])
			if regimen_resis_list:
				var["regimen_sum_for_jdyy"] = "，".join(regimen_resis_list)
				result_resis.append(var)
	if not result_resis:
		if regimen_D_list_resis:
			result_resis = [{"gene_symbol" : "nofoundsignificant"}]
		else:
			result_resis = [{"gene_symbol" : "nofoundvar"}]
	return result_resis
jinja2.filters.FILTERS["jlyy_116gene_result_resis_v2"] = jlyy_116gene_result_resis_v2

def jlyy_116_summary_new_rule_summary(info):
	knb = info[0]
	raw_var_list = info[1]
	result = []
	if knb:
		result.extend(jlyy_lc10_summary_new_rule_single_gene([knb]))
	gene_list = []
	for var in raw_var_list:
		for gene in re.split(",", var["gene_symbol"]):
			if gene not in gene_list:
				gene_list.append(gene)
	for gene in gene_list:
		var_list  = [var for var in raw_var_list if set(re.split(",", var["gene_symbol"])) & set([gene])]
		gene_regimen = jlyy_lc10_summary_new_rule_single_gene(var_list)
		if gene_regimen:
			result.extend(gene_regimen)
	# 去重
	sort_result = []
	for var in result:
		if var not in sort_result:
			sort_result.append(var)
	return sort_result
jinja2.filters.FILTERS["jlyy_116_summary_new_rule_summary"] = jlyy_116_summary_new_rule_summary

# 吉林大学第一医院116-返回治疗-敏感 A级证据总结-2024.09.27
def jlyy_116_sense_a_evi(var):
	var_info = ""
	if var["var_id"] == "KRAS/NRAS/BRAF WT":
		var_info = "KRAS/NRAS/BRAF V600E 野生型"
	else:
		if var["bio_category"] == "Snvindel":
			if var["hgvs_p"] != "p.?":
					var_info = var["gene_symbol"] + " " + var["hgvs_p_ZJZL"]
			else:
				var_info = var["gene_symbol"] + " " + var["hgvs_c"]
		elif var["bio_category"] == "Cnv":
			var_info = var["gene_symbol"] + " 扩增"
		elif var["bio_category"] == "Sv":
			var_info = var["five_prime_gene"]+"-"+var["three_prime_gene"]+" 融合"
	result = []
	evi_dict = {}
	level_a_evi = [evi for evi in var["evi_sum"]["evi_split"]["Predictive"] if evi["clinical_significance_cn"] == "敏感" and evi["evi_conclusion_simple"] == "A"]
	# 2024.11.05-按ABCD排序
	#level_a_evi = jlyy_lc10_sort_evi(level_a_evi)
	for evi in level_a_evi:
		if evi["refer_agency"] not in evi_dict.keys():
			evi_dict.setdefault(evi["refer_agency"], {})
		if evi["tumor_name_cn"] not in evi_dict[evi["refer_agency"]].keys():
			evi_dict[evi["refer_agency"]].setdefault(evi["tumor_name_cn"], [])
		if evi["regimen_name"] not in evi_dict[evi["refer_agency"]][evi["tumor_name_cn"]]:
			evi_dict[evi["refer_agency"]][evi["tumor_name_cn"]].append(evi["regimen_name"])
	for agency, info in evi_dict.items():
		for tumor, regimen_list in info.items():
			# 2024.11.05-引用机构放到前面
			#if agency in ["FDA", "NMPA"]:
			#	result.append("{0}适用于{1}变异的{2}患者（{3}批准）".format("、".join(regimen_list), var_info, tumor, agency))
			#else:
			#	result.append("{0}适用于{1}变异的{2}患者（{3}指南推荐）".format("、".join(regimen_list), var_info, tumor, agency))
			if agency in ["FDA", "NMPA"]:
				result.append("{0}批准{1}适用于{2}变异的{3}患者".format(agency, "、".join(regimen_list), var_info, tumor))
			else:
				result.append("{0}指南推荐{1}适用于{2}变异的{3}患者".format(agency, "、".join(regimen_list), var_info, tumor))
			# 2024.11.05-更新完成
	return result
jinja2.filters.FILTERS["jlyy_116_sense_a_evi"] = jlyy_116_sense_a_evi

# 吉林大学第一医院-116-本癌种重要基因要根据癌种进行区分-2024.09.27
def jlyy_116_important_gene(info):
	var_list = info[0]
	tumor_list = info[1]
	lung_gene = ["ALK", "EGFR", "ROS1", "RET", "MET", "ERBB2", "KRAS", "BRAF", "NTRK1", "NTRK2", "NTRK3"]
	colo_gene = ["KRAS", "NRAS", "BRAF"]
	ga_gene = ["ERBB2", "NTRK1", "NTRK2", "NTRK3", "RET"]
	detect_gene = []
	for var in var_list:
		for gene in re.split(",", var["gene_symbol"]):
			if gene not in detect_gene:
				detect_gene.append(gene)
	if "肠癌" in tumor_list:
		return "，".join(sorted([gene for gene in colo_gene if gene not in detect_gene]))
	elif "胃癌" in tumor_list:
		return "，".join(sorted([gene for gene in ga_gene if gene not in detect_gene]))
	else:
		return "，".join(sorted([gene for gene in lung_gene if gene not in detect_gene]))
jinja2.filters.FILTERS["jlyy_116_important_gene"] = jlyy_116_important_gene

# 吉林大学第一医院-116-本癌种重要基因要根据癌种进行区分-V2-2024.11.05
def jlyy_116_important_gene_v2(info):
	var_list = info[0]
	tumor_list = info[1]
	lung_gene = ["ALK", "EGFR", "ROS1", "RET", "MET", "ERBB2", "KRAS", "BRAF", "NTRK1", "NTRK2", "NTRK3"]
	colo_gene = ["KRAS", "NRAS", "BRAF"]
	ga_gene = ["ERBB2", "NTRK1", "NTRK2", "NTRK3", "RET", "BRAF"]
	tc_gene = ["ALK", "BRAF", "RET", "NTRK1", "NTRK2", "NTRK3"]
	mel_gene = ["BRAF"]
	oc_gene = ["BRCA"]
	gist_gene = ["KIT", "PDGFRA"]
	# 2025.05.15-新增子宫内膜癌分子分型
	ec_gene = ["MLH1", "MSH2", "MSH6", "PMS2", "POLE", "TP53"]
	detect_gene = []
	for var in var_list:
		for gene in re.split(",", var["gene_symbol"]):
			if gene not in detect_gene:
				detect_gene.append(gene)
	if "肠癌" in tumor_list:
		return "，".join(sorted([gene for gene in colo_gene if gene not in detect_gene]))
	elif "胃癌" in tumor_list:
		return "，".join(sorted([gene for gene in ga_gene if gene not in detect_gene]))
	elif "甲状腺癌" in tumor_list:
		return "，".join(sorted([gene for gene in tc_gene if gene not in detect_gene]))
	elif "黑色素瘤" in tumor_list:
		return "，".join(sorted([gene for gene in mel_gene if gene not in detect_gene]))
	elif "卵巢癌" in tumor_list:
		return "，".join(sorted([gene for gene in oc_gene if gene not in detect_gene]))
	elif "胃肠道间质瘤" in tumor_list:
		return "，".join(sorted([gene for gene in gist_gene if gene not in detect_gene]))
	elif "子宫内膜癌" in tumor_list:
		return "，".join(sorted([gene for gene in ec_gene if gene not in detect_gene]))
	else:
		return "，".join(sorted([gene for gene in lung_gene if gene not in detect_gene]))
jinja2.filters.FILTERS["jlyy_116_important_gene_v2"] = jlyy_116_important_gene_v2

# 山东齐鲁ptm plus结果汇总-2024.09.29
def sdql_ptm_plus_summary(var_list):
	gene_list = ["CTNNB1", "EPCAM", "MLH1", "MSH2", "MSH6", "PMS2", "POLE", "TP53"]
	result = []
	result.extend(var_list)
	detect_gene = [var["gene_symbol"] for var in var_list]
	for gene in set(gene_list) - set(detect_gene):
		result.append({"gene_symbol" : gene})
	return sorted(result, key=lambda i:gene_list.index(i["gene_symbol"]))
jinja2.filters.FILTERS["sdql_ptm_plus_summary"] = sdql_ptm_plus_summary

# 云肿150，判断是否检测到3类有遗传风险的变异-2024.10.11
def judge_ynzl_150_level3_risk(var_list):
	result = ""
	for var in var_list:
		if var["evi_sum"] and "evi_split" in var["evi_sum"].keys() and var["evi_sum"]["evi_split"] and \
			"Predisposing" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Predisposing"]:
			result = "True"
			break
	return result
jinja2.filters.FILTERS["judge_ynzl_150_level3_risk"] = judge_ynzl_150_level3_risk

# 处理CP200 HD结果-判断是否存在阳性结果-2024.10.12
def judge_cp200_hd_positive(hd_list):
	if [var for var in hd_list if var["var_auto_result"] == "T"]:
		return True
	else:
		return False
jinja2.filters.FILTERS["judge_cp200_hd_positive"] = judge_cp200_hd_positive

# 处理CP200 HD结果-存在阳性-返回全部结果-2024.10.12
def get_cp200_hd(hd_list):
	result = []
	gene_dict = {
		"CDKN2A" : "E1-E3，NM_000077.5",
		"CDKN2B" : "E1-E2，NM_004936.4",
		"MTAP" : "E1-E8，NM_002451.4"
	}
	for var in hd_list:
		if var["var_auto_result"] == "T":
			var["homodel_result"] = "阳性"
			result.append(var)
	for gene in set(["CDKN2A", "CDKN2B", "MTAP"]) - set([var["gene_symbol"] for var in hd_list if var["var_auto_result"] == "T"]):
		result.append({
			"gene_symbol" : gene,
			"region_transcript_primary" : gene_dict[gene],
			"homodel_result" : "阴性"
		})
	return sorted(result, key=lambda i:i["gene_symbol"])
jinja2.filters.FILTERS["get_cp200_hd"] = get_cp200_hd

# 西安交大一Master DNA/RNA共检SV计算1，展示成两个-2024.10.12
def xajdy_master_var_list(var_list):
	result = []
	for var in var_list:
		result.append(var)
		if "rna_detect" in var.keys() and var["rna_detect"]:
			result.append(var["rna_detect"])
	return result
jinja2.filters.FILTERS["xajdy_master_var_list"] = xajdy_master_var_list

# 西安交大一-输出最高等级药物，AB或CD-2024.10.12
def xajdy_filter_regimen_level(regimen_list):
	ab_regimen_list = [regimen for regimen in regimen_list if regimen["evi_conclusion_simple"] in ["A", "B"]]
	cd_regimen_list = [regimen for regimen in regimen_list if regimen["evi_conclusion_simple"] in ["C", "D"]]
	if ab_regimen_list:
		return ab_regimen_list
	else:
		return cd_regimen_list
jinja2.filters.FILTERS["xajdy_filter_regimen_level"] = xajdy_filter_regimen_level

# 重庆附一CP40-新模板-指定基因检测结果汇总-2024.10.16
def cqfy_cp40_gene_var_list(info):
	var_list = info[0]
	gene = info[1]
	return [var for var in var_list if gene in re.split(",", var["gene_symbol"])]
jinja2.filters.FILTERS["cqfy_cp40_gene_var_list"] = cqfy_cp40_gene_var_list

# 重庆附一CP40-新模板-结论-I类（包含用药）-2024.10.16
def cqfy_cp40_I_sum(var_list):
	# 2024.12.24-gene_list改为10基因
	#gene_list = ["EGFR", "ALK", "BRAF", "ERBB2", "KRAS", "MET", "NTRK1", "NTRK2", "NTRK3", "RET", "ROS1"]
	gene_list = ["EGFR", "ALK", "BRAF", "ERBB2", "KRAS", "MET", "NRAS", "PIK3CA", "RET", "ROS1"]
	# 2024.12.24-更新完成
	# 排序，上面基因的变异优先展示
	for var in var_list:
		# 兼容MSI-2024.12.24
		if var["var_id"] in ["MSI-H", "MSS"]:
			var["gene_symbol"] = "MSS"
			var["bio_category"] = "MSS"
		# 2024.12.24-兼容完成
		var["cqfy_cp40_sort_rule"] = 99
		for gene in re.split(",", var["gene_symbol"]):
			if gene in gene_list:
				var["cqfy_cp40_sort_rule"] = gene_list.index(gene)
				break
	var_list = sorted(var_list, key = lambda i:i["cqfy_cp40_sort_rule"])

	result = []
	# EGFR-TKI药物-区分一二三代
	# 2025.01.21-EGFR-TKI三代新增药物利厄替尼
	egfr_tki_dict = {
		"first" : ["吉非替尼", "厄洛替尼", "埃克替尼"],
		"second" : ["阿法替尼", "达可替尼"],
		"third" : ["奥希替尼", "阿美替尼", "伏美替尼", "贝福替尼", "瑞厄替尼", "拉泽替尼", "瑞齐替尼", "利厄替尼"],
		"new" : ["佐利替尼"]
	}
	# EGFR-TKI药物-不区分几代
	egfr_tik_list = egfr_tki_dict["first"] + egfr_tki_dict["second"] + egfr_tki_dict["third"] + egfr_tki_dict["new"]
	for var in var_list:
		# 只输出A/B中最高等级药物
		sense_A = [regimen["regimen_name"] for regimen in var["evi_sum"]["regimen_S"] if var["evi_sum"]["regimen_S"] and regimen["evi_conclusion_simple"] == "A"]
		sense_B = [regimen["regimen_name"] for regimen in var["evi_sum"]["regimen_S"] if var["evi_sum"]["regimen_S"] and regimen["evi_conclusion_simple"] == "B"]
		resis_A = [regimen["regimen_name"] for regimen in var["evi_sum"]["regimen_R"] if var["evi_sum"]["regimen_R"] and regimen["evi_conclusion_simple"] == "A"]
		resis_B = [regimen["regimen_name"] for regimen in var["evi_sum"]["regimen_R"] if var["evi_sum"]["regimen_R"] and regimen["evi_conclusion_simple"] == "B"]
		sense_regimen = sense_A if sense_A else sense_B if sense_B else []
		resis_regimen = resis_A if resis_A else resis_B if resis_B else []
		# 变异类型1：EGFR 19del/L858R/L861Q/S768I/G719。
		# 变异类型2：EGFR T790M
		# 变异类型3：其他变异
		var_type = "other"
		if var["bio_category"] == "Snvindel":
			hgvs_p = var["hgvs_p_ZJZL"].replace("p.", "")
			ref = var["ref"] if var["ref"] != "-" else ""
			alt = var["alt"] if var["alt"] != "-" else ""
			if hgvs_p in ["L858R", "L861Q", "S768I"]:
				var_type = "egfr_sense"
			elif "G719" in hgvs_p:
				var_type = "egfr_sense"
			elif var["gene_region"] == "exon19" and "del" in var["hgvs_p"] and len(ref) > len(alt):
				var_type = "egfr_sense"
			elif hgvs_p == "T790M":
				var_type = "egfr_resis"
		tmp_list = []
		# 根据变异类型来组装结果
		if var_type == "egfr_sense":
			# 对EGFR-TKI（……）、……敏感，对……耐药。
			sense_egfr = [regimen for regimen in sense_regimen if regimen in egfr_tik_list]
			sense_other = [regimen for regimen in sense_regimen if regimen not in egfr_tik_list]
			if sense_egfr and sense_other:
				tmp_list.append("对EGFR-TKI（{0}）、{1}敏感".format("、".join(sense_egfr), "、".join(sense_other)))
			elif sense_egfr and not sense_other:
				tmp_list.append("对EGFR-TKI（{0}）敏感".format("、".join(sense_egfr)))
			elif not sense_egfr and sense_other:
				tmp_list.append("对{0}敏感".format("、".join(sense_other)))
			if resis_regimen:
				tmp_list.append("对{0}耐药".format("、".join(resis_regimen)))
		elif var_type == "egfr_resis":
			# 对三代EGFR-TKI（……）、……敏感，对一、二代EGFR-TKI（……）、……耐药。
			sense_egfr = [regimen for regimen in sense_regimen if regimen in egfr_tki_dict["third"]]
			sense_other = [regimen for regimen in sense_regimen if regimen not in egfr_tki_dict["third"]]
			resis_egfr = [regimen for regimen in resis_regimen if regimen in egfr_tki_dict["first"] + egfr_tki_dict["second"]]
			resis_other = [regimen for regimen in resis_regimen if regimen not in egfr_tki_dict["first"] + egfr_tki_dict["second"]]
			if sense_egfr and sense_other:
				tmp_list.append("对三代EGFR-TKI（{0}）、{1}敏感".format("、".join(sense_egfr), "、".join(sense_other)))
			elif sense_egfr and not sense_other:
				tmp_list.append("对三代EGFR-TKI（{0}）敏感".format("、".join(sense_egfr)))
			elif not sense_egfr and sense_other:
				tmp_list.append("对{0}敏感".format("、".join(sense_other)))
			
			if resis_egfr and resis_other:
				tmp_list.append("对一、二代EGFR-TKI（{0}）、{1}耐药".format("、".join(resis_egfr), "、".join(resis_other)))
			elif resis_egfr and not resis_other:
				tmp_list.append("对一、二代EGFR-TKI（{0}）耐药".format("、".join(resis_egfr)))
			elif not resis_egfr and resis_other:
				tmp_list.append("对{0}耐药".format("、".join(resis_other)))
		elif var_type == "other":
			# 对……敏感，对……耐药。
			if sense_regimen:
				tmp_list.append("对{0}敏感".format("、".join(sense_regimen)))
			if resis_regimen:
				tmp_list.append("对{0}耐药".format("、".join(resis_regimen)))
		var["cqfy_cp40_var_sum"] = "，".join(tmp_list)
		if var["cqfy_cp40_var_sum"]:
			result.append(var)
	# 每个变异（除了最后一个）cqfy_cp40_var_sum加个分号，在模板里使用for循环展示变异
	if len(result) > 1:
		for var in result[0:-1]:
			var["cqfy_cp40_var_sum"] = var["cqfy_cp40_var_sum"]+"；"
	return result
jinja2.filters.FILTERS["cqfy_cp40_I_sum"] = cqfy_cp40_I_sum

# 重庆附一CP40-新模板-结论-I类-未检测到变异的基因-2024.10.16
def cqfy_cp40_I_nofundvar_gene(var_list):
	gene_list = ["EGFR", "ALK", "BRAF", "ERBB2", "KRAS", "MET", "NTRK1", "NTRK2", "NTRK3", "RET", "ROS1"]
	detect_gene = []
	for var in var_list:
		for gene in re.split(",", var["gene_symbol"]):
			if gene not in detect_gene:
				detect_gene.append(gene)
	result = []
	for gene in gene_list:
		if gene not in detect_gene:
			if gene == "ERBB2":
				result.append("ERBB2(HER2)")
			else:
				result.append(gene)
	if set(detect_gene) & set(gene_list):
		return ""
	else:
		return "、".join(result)
jinja2.filters.FILTERS["cqfy_cp40_I_nofundvar_gene"] = cqfy_cp40_I_nofundvar_gene

# 重庆附一CP40-新模板-结论-II、III类-2024.10.16
def cqfy_cp_var_II_III_sum(var_list):
	gene_list = ["EGFR", "ALK", "BRAF", "ERBB2", "KRAS", "MET", "NTRK1", "NTRK2", "NTRK3", "RET", "ROS1"]
	# 排序，上面基因的变异优先展示
	for var in var_list:
		var["cqfy_cp40_sort_rule"] = 99
		for gene in re.split(",", var["gene_symbol"]):
			if gene in gene_list:
				var["cqfy_cp40_sort_rule"] = gene_list.index(gene)
				break
	var_list = sorted(var_list, key = lambda i:i["cqfy_cp40_sort_rule"])

	for var in var_list:
		if var["bio_category"] == "Snvindel":
			# c+p展示-2024.10.21
			#var["cqfy_var_info"] = var["hgvs_p_ZJZL"] if var["hgvs_p_ZJZL"] != "p.?" else var["hgvs_c"]
			if var["hgvs_p_ZJZL"] != "p.?":
				var["cqfy_var_info"] = var["hgvs_c"]+" "+var["hgvs_p_ZJZL"]
			else:
				var["cqfy_var_info"] = var["hgvs_c"]
			if "judge_mergeMET" in var.keys() and var["judge_mergeMET"]:
				var["cqfy_var_info"] += "（MET exon14 skipping）"
		elif var["bio_category"] == "Cnv":
			var["cqfy_var_info"] = "基因扩增"
		elif var["bio_category"] == "Sv":
			if var["five_prime_gene"] == "MET" and var["three_prime_gene"] == "MET":
				var["cqfy_var_info"] = "14号外显子跳跃突变"
			else:
				var["cqfy_var_info"] = "基因融合"
	if len(var_list) > 1:
		for var in var_list[0:-1]:
			var["cqfy_var_info"] = var["cqfy_var_info"]+"、"
	return var_list
jinja2.filters.FILTERS["cqfy_cp_var_II_III_sum"] = cqfy_cp_var_II_III_sum

# 重庆附一CP40-新模板-结论-I、II、III类-2024.10.16
def cqfy_cp_var_sum(var_list):
	gene_list = ["EGFR", "ALK", "BRAF", "ERBB2", "KRAS", "MET", "NTRK1", "NTRK2", "NTRK3", "RET", "ROS1"]
	# 排序，上面基因的变异优先展示
	for var in var_list:
		var["cqfy_cp40_sort_rule"] = 99
		for gene in re.split(",", var["gene_symbol"]):
			if gene in gene_list:
				var["cqfy_cp40_sort_rule"] = gene_list.index(gene)
				break
	var_list = sorted(var_list, key = lambda i:i["cqfy_cp40_sort_rule"])

	for var in var_list:
		if var["bio_category"] == "Snvindel":
			# c+p展示-2024.10.21
			#var["cqfy_var_info"] = var["hgvs_p_ZJZL"] if var["hgvs_p_ZJZL"] != "p.?" else var["hgvs_c"]
			if var["hgvs_p_ZJZL"] != "p.?":
				var["cqfy_var_info"] = var["hgvs_p_ZJZL"]
			else:
				var["cqfy_var_info"] = var["hgvs_c"]
			if "judge_mergeMET" in var.keys() and var["judge_mergeMET"]:
				var["cqfy_var_info"] += "（MET exon14 skipping）"
		elif var["bio_category"] == "Cnv":
			var["cqfy_var_info"] = "基因扩增"
		elif var["bio_category"] == "Sv":
			if var["five_prime_gene"] == "MET" and var["three_prime_gene"] == "MET":
				var["cqfy_var_info"] = "14号外显子跳跃突变"
			else:
				var["cqfy_var_info"] = "基因融合"
	if len(var_list) > 1:
		for var in var_list[0:-1]:
			var["cqfy_var_info"] = var["cqfy_var_info"]+"、"
	return var_list
jinja2.filters.FILTERS["cqfy_cp_var_sum"] = cqfy_cp_var_sum

# 西安交大一HRD-HRD临床提示，HRD最高等级+BRCA最高等级-2024.10.23
def xajdy_hrd_regimen(info):
	result = []
	hrd = info[0]
	brca_list = info[1]
	hrd_regimen_sum = [
		{
			"regimen_name":i["regimen_name"], 
			"evidence_type":i["evidence_type"], 
			"clinical_significance_cn":i["clinical_significance_cn"], 
			"evi_conclusion_simple":i["evi_conclusion_simple"],
			"regimen_name_py":i["regimen_name_py"]
   		} 
   		for i in hrd["evi_sum"]["regimen_evi_sum"] if i["evidence_type"] in ["Diagnostic","Predictive","Prognostic"]] \
			if hrd["evi_sum"] and "regimen_evi_sum" in hrd["evi_sum"].keys() and hrd["evi_sum"]["regimen_evi_sum"] else []
	hrd_regimen_sum_AB = [regimen for regimen in hrd_regimen_sum if regimen["evi_conclusion_simple"] in ["A", "B"]]
	hrd_regimen_sum_CD = [regimen for regimen in hrd_regimen_sum if regimen["evi_conclusion_simple"] in ["C", "D"]]
	if hrd_regimen_sum_AB:
		result.extend(hrd_regimen_sum_AB)
	else:
		result.extend(hrd_regimen_sum_CD)
	
	brca_regimen_sum = []
	if brca_list:
		for var in brca_list:
			brca_regimen_sum += [
				{
					"regimen_name":i["regimen_name"], 
					"evidence_type":i["evidence_type"], 
					"clinical_significance_cn":i["clinical_significance_cn"], 
					"evi_conclusion_simple":i["evi_conclusion_simple"],
					"regimen_name_py":i["regimen_name_py"]
				} 
				for i in var["evi_sum"]["regimen_evi_sum"] if i["evidence_type"] in ["Diagnostic","Predictive","Prognostic"]]
	brca_regimen_sum_AB = [regimen for regimen in brca_regimen_sum if regimen["evi_conclusion_simple"] in ["A", "B"]]
	brca_regimen_sum_CD = [regimen for regimen in brca_regimen_sum if regimen["evi_conclusion_simple"] in ["C", "D"]]
	if brca_regimen_sum_AB:
		result.extend(brca_regimen_sum_AB)
	else:
		result.extend(brca_regimen_sum_CD)

	result_redup = []
	for i in result:
		if i not in result_redup:
			result_redup.append(i)

	return sorted(result_redup, key=lambda i:(i["evi_conclusion_simple"], i["clinical_significance_cn"], i["regimen_name_py"]))
jinja2.filters.FILTERS["xajdy_hrd_regimen"] = xajdy_hrd_regimen

# 内蒙古人民-HRR-结果小结处分为BRCA和其他基因-2024.10.29
def nmrm_hrr_filter_var(info):
	var_list = info[0]
	gene_type = info[1]
	#result = []
	if gene_type == "brca":
		result = [var for var in var_list if var["gene_symbol"] in ["BRCA1", "BRCA2"]]
		for gene in set(["BRCA1", "BRCA2"]) - set([var["gene_symbol"] for var in var_list]):
			result.append({"gene_symbol" : gene})
		return sorted(result, key = lambda i:i["gene_symbol"])
	else:
		return [var for var in var_list if var["gene_symbol"] not in ["BRCA1", "BRCA2"]]
	#return result
jinja2.filters.FILTERS["nmrm_hrr_filter_var"] = nmrm_hrr_filter_var

# 吉林大学第一医院-临床提示耐药放到敏感后面-2024.11.07
def jlyy_regimen_sort(regimen_list):
	return sorted(regimen_list, key = lambda i:i["sense_rule"])
jinja2.filters.FILTERS["jlyy_regimen_sort"] = jlyy_regimen_sort


# 适用北大三MP-2024.11.11
def io_detect_for_BDS_MP(info):
	var_list = info[0]
	return_type = info[1]
	# 返回结果中的io_result用于填充IO表，", ".join(io_p_list), ", ".join(io_n_list)用于填充检测结果小结
	io_result = {}
	# 汇总体细胞I/II/肿瘤发生发展相关变异+胚系致病/疑似致病性变异
	io_gene_P = ["ATM","ATR","BRCA1","BRCA2","BRIP1","CHEK1","CHEK2","ERCC1","FANCA","MRE11",\
				 "PALB2","RAD50","XRCC1","MLH1","MSH2","MSH6","PMS2","POLE","POLD1","TP53",\
				 "KRAS","CD274","ARID1A","LRP1B","SETD2","PRKDC","TERT","KMT2D","FAT1","CDK12",\
				 "SERPINB3","SERPINB4"]
	io_gene_N = ["EGFR","ALK","MDM2","MDM4","CDKN2A","CDKN2B","DNMT3A","STK11","IFNGR1","IRF1",\
				 "JAK1","JAK2","APC","CTNNB1","B2M","PTEN","CCND1","FGF3","FGF19"]
	cnv_gene_list = ["CD274", "MDM2", "MDM4"]
	both_cnv_list = ["CCND1","FGF3","FGF19"]
	
	#level_12_var = [var for var in var_data if judge_var(var, [4,5], [4,5])]

	for var in var_list:
		# 仅展示扩增的基因
		if var["bio_category"] == "Cnv" and var["gene_symbol"] in cnv_gene_list:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			io_result[var["gene_symbol"]].append(var)
		# 仅展示融合的基因
		elif var["bio_category"] in ["Sv", "PSeqRnaSv"] and set(re.split(",", var["gene_symbol"])) & set(["ALK"]):
			if "ALK" not in io_result.keys():
				io_result.setdefault("ALK", [])
			io_result["ALK"].append(var)
		# 其余基因展示Snvindel
		elif var["bio_category"] == "Snvindel" and var["gene_symbol"] in io_gene_P + io_gene_N:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			io_result[var["gene_symbol"]].append(var)

	# summary展示
	io_p_list = [i for k,v in io_result.items() for i in v if k not in both_cnv_list and k in io_gene_P]
	io_n_list = [i for k,v in io_result.items() for i in v if k not in both_cnv_list and k in io_gene_N]

	# 处理CNV共突变
	judge_both_cnv = ""
	if "CCND1" in io_result.keys() and "FGF3" in io_result.keys() and "FGF19" in io_result.keys():
		io_n_list.append({"gene_symbol" : "CCND1/FGF3/FGF19", "var_type" : "扩增共突变"})
	
	if return_type == "p":
		return io_p_list
	elif return_type == "n":
		return io_n_list
jinja2.filters.FILTERS["io_detect_for_BDS_MP"] = io_detect_for_BDS_MP

# 适用北大三MP-2025.06.17
# 增加HD
def io_detect_for_BDS_MP_v2(info):
	var_list = info[0]
	return_type = info[1]
	# 返回结果中的io_result用于填充IO表，", ".join(io_p_list), ", ".join(io_n_list)用于填充检测结果小结
	io_result = {}
	# 汇总体细胞I/II/肿瘤发生发展相关变异+胚系致病/疑似致病性变异
	io_gene_P = ["ATM","ATR","BRCA1","BRCA2","BRIP1","CHEK1","CHEK2","ERCC1","FANCA","MRE11",\
				 "PALB2","RAD50","XRCC1","MLH1","MSH2","MSH6","PMS2","POLE","POLD1","TP53",\
				 "KRAS","CD274","ARID1A","LRP1B","SETD2","PRKDC","TERT","KMT2D","FAT1","CDK12",\
				 "SERPINB3","SERPINB4"]
	io_gene_N = ["EGFR","ALK","MDM2","MDM4","CDKN2A","CDKN2B","DNMT3A","STK11","IFNGR1","IRF1",\
				 "JAK1","JAK2","APC","CTNNB1","B2M","PTEN","CCND1","FGF3","FGF19"]
	cnv_gene_list = ["CD274", "MDM2", "MDM4", "CCND1","FGF3","FGF19"]
	both_cnv_list = ["CCND1","FGF3","FGF19"]
	# hd 是额外需要展示的，原有的输出类型不变
	hd_gene_list = ["ATM", "BRCA1", "BRCA2", "BRIP1", "CDK12", "CHEK1", "CHEK2", "FANCA", \
				 	"PALB2", "SETD2", "TP53", "CDKN2A", "CDKN2B", "PTEN", "STK11"]
	
	#level_12_var = [var for var in var_data if judge_var(var, [4,5], [4,5])]

	for var in var_list:
		# 仅展示扩增的基因
		if var["bio_category"] == "Cnv" and var["gene_symbol"] in cnv_gene_list:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			io_result[var["gene_symbol"]].append(var)
		# 仅展示融合的基因
		elif var["bio_category"] in ["Sv", "PSeqRnaSv"] and set(re.split(",", var["gene_symbol"])) & set(["ALK"]):
			if "ALK" not in io_result.keys():
				io_result.setdefault("ALK", [])
			io_result["ALK"].append(var)
		# HD基因经确认同时展示snvindel和hd
		elif var["bio_category"] == "Snvindel" or var["bio_category"] == "PHd" and var["gene_symbol"] in hd_gene_list:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			io_result[var["gene_symbol"]].append(var)
		# 其余基因展示Snvindel
		elif var["bio_category"] == "Snvindel" and var["gene_symbol"] in io_gene_P + io_gene_N:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			io_result[var["gene_symbol"]].append(var)

	# summary展示
	io_p_list = [i for k,v in io_result.items() for i in v if k not in both_cnv_list and k in io_gene_P]
	io_n_list = [i for k,v in io_result.items() for i in v if k not in both_cnv_list and k in io_gene_N]

	# 处理CNV共突变
	judge_both_cnv = ""
	if "CCND1" in io_result.keys() and "FGF3" in io_result.keys() and "FGF19" in io_result.keys():
		io_n_list.append({"gene_symbol" : "CCND1/FGF3/FGF19", "var_type" : "扩增共突变"})
	
	if return_type == "p":
		return io_p_list
	elif return_type == "n":
		return io_n_list
jinja2.filters.FILTERS["io_detect_for_BDS_MP_v2"] = io_detect_for_BDS_MP_v2

# 复旦中山200基因，汇总，其中变异涉及基因数190个，其中一个EBV（不展示了），化疗6个，MSI 1个，这边只汇总190-2024.07.26
# 加上HD结果（仅胶质瘤、质控合格时展示）
def fdzs_190_gene_sum_v2(info):
	var_list = info[0]
	hd = info[1]
	tumor_list = info[2]
	snp_cover_ratio_num = info[3]
	gene_list = ["ABRAXAS1", "AKT1", "AKT2", "AKT3", "ALK", "APC", "AR", "ARAF", "ARID1A", "ARID1B", \
				 "ARID2", "ATM", "ATR", "ATRX", "AURKA", "AXIN1", "B2M", "BAP1", "BARD1", "BMPR1A", \
				 "BRAF", "BRCA1", "BRCA2", "BRD3", "BRD4", "BRIP1", "CCND1", "CCNE1", "CD274", "CDH1", \
				 "CDK12", "CDK4", "CDK6", "CDKN1B", "CDKN2A", "CDKN2B", "CHD1", "CHEK1", "CHEK2", "CLDN18", \
				 "CLIP1", "CTNNB1", "CYP17A1", "DDR2", "DICER1", "DNAJB1", "DNMT3A", "EGFR", "EIF1AX", \
				 "ELK4", "EPCAM", "ERBB2", "ERBB3", "ERCC2", "ERG", "ESR1", "ETV1", "ETV4", "ETV5", \
				 "EZH2", "FANCA", "FANCC", "FANCD2", "FANCL", "FAT1", "FBXW7", "FGF19", "FGFR1", "FGFR2", \
				 "FGFR3", "FGFR4", "FH", "FLT3", "FOXA1", "GATA3", "GATA6", "GEN1", "GLI1", "GLI2", \
				 "GNA11", "GNA14", "GNAQ", "GNAS", "GTF2I", "H3-3A", "H3C2", "H3C3", "HDAC2", "HOXB13", \
				 "HRAS", "HSD3B1", "IDH1", "IDH2", "IFNGR1", "IFNGR2", "JAK1", "JAK2", "KDM5C", "KDM6A", \
				 "KEAP1", "KIT", "KMT2B", "KMT2C", "KMT2D", "KRAS", "MAP2K1", "MAP2K4", "MDM2", "MDM4", \
				 "MET", "MLH1", "MLH3", "MMS22L", "MRE11", "MSH2", "MSH3", "MSH6", "MTAP", "MTOR", \
				 "MUTYH", "MYB", "MYC", "MYCN", "NBN", "NF1", "NF2", "NFE2L2", "NKX2-1", "NRAS", \
				 "NRG1", "NSD3", "NTRK1", "NTRK2", "NTRK3", "NUTM1", "PALB2", "PAX8", "PBRM1", "PDCD1LG2", \
				 "PDGFRA", "PIK3CA", "PIK3R1", "PMS2", "POLD1", "POLE", "PPP2R1A", "PRKACA", "PTCH1", "PTCH2", \
				 "PTEN", "QKI", "RAD50", "RAD51B", "RAD51C", "RAD51D", "RAD54L", "RAF1", "RASA1", "RASGRF1", \
				 "RB1", "RBM10", "RET", "RICTOR", "RIT1", "RNF43", "ROS1", "SETD2", "SF3B1", "SMAD4", \
				 "SMARCA2", "SMARCA4", "SMARCB1", "SMO", "SPOP", "STAT3", "STK11", "SUFU", "TERT", "TET2", \
				 "TGFBR2", "TP53", "TSC1", "TSC2", "VEGFA", "VHL", "YAP1", "ZFTA", "ZNF532", "ZNF592"]
	
	if "胶质瘤" in tumor_list and float(snp_cover_ratio_num) >= 0.9:
		for var in hd:
			if var["var_auto_result"] == "T":
				var_list.append({
					"gene_symbol" : var["gene_symbol"],
					"region" : var["region"],
					"bio_category" : "hd",
					"var_id" : "hd"
				})
	
	result = []
	
	detect_gene = []
	for var in var_list:
		for gene in set(re.split(",", var["gene_symbol"])):
			detect_gene.append(gene)
			if gene in gene_list:
				result.append(var)
	for gene in set(gene_list) - set(detect_gene):
		result.append({"gene_symbol" : gene})
	
	return sorted(result, key = lambda i:i["gene_symbol"])
jinja2.filters.FILTERS["fdzs_190_gene_sum_v2"] = fdzs_190_gene_sum_v2

# 北大三-MP-基因变异和解读中，证据分为敏感和耐药，其中预后较好和辅助诊断放在敏感中，预后较差的放在耐药中-2024.11.14
# 预后中等的先放在敏感中，有需要再改
def bds_mp_prognostic(info):
	evi_type = info[0]
	evi_list = info[1]
	poor_evi = [evi for evi in evi_list if evi["clinical_significance_cn"] == "较差"]
	better_evi = [evi for evi in evi_list if evi["clinical_significance_cn"] != "较差"]
	if evi_type == "Poor":
		return poor_evi
	else:
		return better_evi
jinja2.filters.FILTERS["bds_mp_prognostic"] = bds_mp_prognostic

# 判断变异是否有治疗方案ABC等级证据+预后/辅助诊断ABC等级证据-2024.11.14
def judge_abc_evi_regimen_v2(var):
	regimen_level = [i["evi_conclusion_simple"] for i in var["evi_sum"]["regimen_evi_sum"] if i["evidence_type"] in ["Predictive","Prognostic", "Diagnostic"]] \
					 if var["evi_sum"] and "regimen_evi_sum" in var["evi_sum"].keys() \
					 else []
	if set(["A", "B", "C"]) & set(regimen_level):
		return True
	else:
		return False
jinja2.filters.FILTERS["judge_abc_evi_regimen_v2"] = judge_abc_evi_regimen_v2

# 检测结果详细解读中敏感/耐药证据-2024.11.14
def get_bds_mp_inter_evi(info):
	var = info[0]
	sense_or_resis = info[1]
	clinical_list_sense = ["敏感"]
	clinical_list_resis = ["耐药"]
	clinical_list = clinical_list_sense if sense_or_resis == "sense" else clinical_list_resis
	evi_list_ABC = [i for i in var["evi_sum"]["regimen_evi_sum"] if i["evidence_type"] in ["Predictive"] and \
		 			 i["evi_conclusion_simple"] in ["A", "B", "C"] and i["clinical_significance_cn"] in clinical_list] \
					 if var["evi_sum"] and "regimen_evi_sum" in var["evi_sum"].keys() \
					 else []

	evi_list_D = [i for i in var["evi_sum"]["regimen_evi_sum"] if i["evidence_type"] in ["Predictive"] and \
		 			 i["evi_conclusion_simple"] in ["D"] and i["clinical_significance_cn"] in clinical_list] \
					 if var["evi_sum"] and "regimen_evi_sum" in var["evi_sum"].keys() \
					 else []
	if var["top_level"] in ["A", "B", "C"]:
		return evi_list_ABC
	else:
		return evi_list_D
jinja2.filters.FILTERS["get_bds_mp_inter_evi"] = get_bds_mp_inter_evi

# 上海肺科CP200-HD-2024.11.19
def shfk_cp200_hd_filter(info):
	hd = info[0]
	hd_type = info[1]
	posi = [var["gene_symbol"] for var in hd if var["var_auto_result"] == "T"]
	nega = []
	for gene in set(["CDKN2A", "CDKN2B", "MTAP"]) - set([var["gene_symbol"] for var in hd if var["var_auto_result"] == "T"]):
		nega.append(gene)
	if hd_type == "positive":
		return "、".join(posi)
	else:
		return "、".join(nega)
jinja2.filters.FILTERS["shfk_cp200_hd_filter"] = shfk_cp200_hd_filter

# 北大三MP-变异排序按基因首字母-2024.11.21
def bds_mp_var_sort(var_list):
	return sorted(var_list, key = lambda i:i["gene_symbol"])
jinja2.filters.FILTERS["bds_mp_var_sort"] = bds_mp_var_sort

# 聊城人民BPTM Plus，展示每个基因的检出结果-2024.11.25
# 与前一个相比，表格多了丰度和变异解读列
def lcrm_gene_detect_bptmplus_v2(info):
	gene = info[1]
	return [var for var in info[0] if var["gene_symbol"] == gene]
jinja2.filters.FILTERS["lcrm_gene_detect_bptmplus_v2"] = lcrm_gene_detect_bptmplus_v2

def bds_cp_qc(info):
	lib_qc = info[0]
	ngs_qc = info[1]
	qc_standard = {
		"cleandata_q30_num" : 0.75,
		"depth_ssbc_num" : 400,
		"depth_rna_ctrl_num" : 20
	}
	qc_dict = {
		"cleandata_q30_num" : "Q30",
		"depth_ssbc_num" : "平均有效深度",
		"depth_rna_ctrl_num" : "RNA内参绝对拷贝数"
	}
	Fail_item = []
	if lib_qc and "lib_dna_qc" in lib_qc.keys() and lib_qc["lib_dna_qc"] and "library_concn" in lib_qc["lib_dna_qc"].keys() and \
		lib_qc["lib_dna_qc"]["library_concn"] and is_number(lib_qc["lib_dna_qc"]["library_concn"]) and float(lib_qc["lib_dna_qc"]["library_concn"]) < 10:
		#print ("1")
		Fail_item.append("文库DNA浓度")
	for item in qc_standard.keys():
		if item in ngs_qc.keys() and ngs_qc[item] and is_number(ngs_qc[item]) and float(ngs_qc[item]) < qc_standard.get(item):
			Fail_item.append(qc_dict.get(item))
	#print (Fail_item)
	return "，".join(Fail_item)
jinja2.filters.FILTERS["bds_cp_qc"] = bds_cp_qc

# 北大三CP40-第X外显子XX突变/第X内含子区突变-2024.11.25
def bds_cp_var_info(var):
	type_dict = {
		"Intronic" : "内含子区突变",
		"3'UTR" : "3'UTR区突变",
		"5'UTR" : "5'UTR区突变",
		"FlankingRegion3" : "侧翼区突变",
		"FlankingRegion5" : "侧翼区突变"
	}
	if re.search("exon|intron", var["gene_region"]):
		region = "第"+gene_region_strn(var["gene_region"]).replace("号", "")
	else:
		region = gene_region_strn(var["gene_region"]).replace("号", "")
	type_cn = var["type_cn"] if var["type_cn"] != "--" else type_dict.get(var["type"], var["type"])
	result = region + type_cn
	result = result.replace("经典剪接位点突变", "内含子区突变")
	result = result.replace("内含子内含子区突变", "内含子区突变")
	return result
jinja2.filters.FILTERS["bds_cp_var_info"] = bds_cp_var_info

# 吉林省肿瘤医院gBRCA，结果小结需要拼接-2024.11.26-v2
# mlpa del/dup不再固定为4/3类，改为动态的
def jlzl_gbrca_summary_v2(var_brca):
	result = []
	level_5_count = len(var_brca["snv_s"]["B1_L5"] + var_brca["snv_s"]["B2_L5"] + var_brca["mlpa_v2"]["B1_mlpa_L5"] + var_brca["mlpa_v2"]["B2_mlpa_L5"])
	level_4_count = len(var_brca["snv_s"]["B1_L4"] + var_brca["snv_s"]["B2_L4"] + var_brca["mlpa_v2"]["B1_mlpa_L4"] + var_brca["mlpa_v2"]["B2_mlpa_L4"])
	level_3_count = len(var_brca["snv_s"]["B1_L3"] + var_brca["snv_s"]["B2_L3"] + var_brca["mlpa_v2"]["B1_mlpa_L3"] + var_brca["mlpa_v2"]["B2_mlpa_L3"])
	if level_5_count:
		result.append("{0}项致病性变异".format(str(level_5_count)))
	if level_4_count:
		result.append("{0}项疑似致病性变异".format(str(level_4_count)))
	if level_3_count:
		result.append("{0}项意义不明确变异".format(str(level_3_count)))
	return "、".join(result)
jinja2.filters.FILTERS["jlzl_gbrca_summary_v2"] = jlzl_gbrca_summary_v2

# 吉林省肿瘤医院gBRCA，结果小结需要拼接(CNV代替MLPA)-2024.09.13
def jlzl_gbrca_summary_cnv_v2(var_brca):
	result = []
	level_5_count = len(var_brca["snv_s"]["B1_L5"] + var_brca["snv_s"]["B2_L5"] + var_brca["gcnv_v2"]["B1_gcnv_L5"] + var_brca["gcnv_v2"]["B2_gcnv_L5"])
	level_4_count = len(var_brca["snv_s"]["B1_L4"] + var_brca["snv_s"]["B2_L4"] + var_brca["gcnv_v2"]["B1_gcnv_L4"] + var_brca["gcnv_v2"]["B2_gcnv_L4"])
	level_3_count = len(var_brca["snv_s"]["B1_L3"] + var_brca["snv_s"]["B2_L3"] + var_brca["gcnv_v2"]["B1_gcnv_L3"] + var_brca["gcnv_v2"]["B2_gcnv_L3"])
	if level_5_count:
		result.append("{0}项致病性变异".format(str(level_5_count)))
	if level_4_count:
		result.append("{0}项疑似致病性变异".format(str(level_4_count)))
	if level_3_count:
		result.append("{0}项意义不明确变异".format(str(level_3_count)))
	return "、".join(result)
jinja2.filters.FILTERS["jlzl_gbrca_summary_cnv_v2"] = jlzl_gbrca_summary_cnv_v2

# 贵州肿瘤BRCA检测结果-2024.11.27
# 阳性时：检出BRCA1基因p.Q858*突变
# v2-CNV Loss分为HeteDel、HomeDel和Loss三种
def gzzl_brca_sum_v2(var_list):
	result = []
	for var in var_list:
		if var["type"] == "Loss":
			if "cnv_type" in var.keys() and var["cnv_type"]:
				if var["cnv_type"] == "HeteDel":
					result.append("{0}基因{1}杂合大片段缺失".format(var["gene_symbol"], var["value"])) 
				elif var["cnv_type"] == "HomoDel":
					result.append("{0}基因{1}纯合大片段缺失".format(var["gene_symbol"], var["value"])) 
				else:
					result.append("{0}基因{1}大片段缺失".format(var["gene_symbol"], var["value"]))
			else: 
				result.append("{0}基因{1}大片段缺失".format(var["gene_symbol"], var["value"])) 
		# 2024.11.22-增加一个Gain
		elif var["type"] == "Gain":
			result.append("{0}基因{1}大片段重复".format(var["gene_symbol"], var["value"])) 
		# 2024.11.22-增加完成
		else:
			if var["hgvs_p"] != "p.?":
				result.append("{0}基因{1}突变".format(var["gene_symbol"], var["hgvs_p"]))
			else:
				result.append("{0}基因{1}突变".format(var["gene_symbol"], var["hgvs_c"]))
	return "检出"+"、".join(result)
jinja2.filters.FILTERS["gzzl_brca_sum_v2"] = gzzl_brca_sum_v2

# 北大三-MP-筛选出不合格的质控项-2024.11.28
def bds_mp_qc(info):
	lib_qc = info[0]
	ngs_qc = info[1]
	prod_type = info[2]
	ngs_qc_standard = {
		"dna_v1" : {
			"cleandata_q30_num" : 0.75,
			"cover_ratio_num" : 0.95,
			"uni20_uniq_hot_num" : 0.9,
			"uni20_uniq_nonhot_num" : 0.8,
			"depth_mean_uniq_hot_num" : 1000,
			"depth_mean_uniq_nonhot_num" : 500
		},
		"rna_v1" : {
			"cleandata_q30_num" : 0.75,
			"read_ontarget_ratio_num" : 0.5,
			"mapping_ratio_num" : 0.8,
			"end2sense_ratio_num" : 0.9,
			"totaldata" : 2
		},
		"dna_v2" : {
			"cleandata_q30_num" : 0.75,
			"cover_ratio_num" : 0.95,
			"uni20_uniq_hot_num" : 0.9,
			"uni20_uniq_nonhot_num" : 0.8,
			"depth_mean_uniq_hot_num" : 800,
			"depth_mean_uniq_nonhot_num" : 400
		},
		"rna_v2" : {
			"cleandata_q30_num" : 0.75,
			"mapping_ratio_num" : 0.8,
			"end2sense_ratio_num" : 0.9,
			"effectivereads_num" : 4000000
		}
	}

	lib_qc_standard = {
		"dna" : {
			"dna_qty" : 150,
			"break_dna_qty" : 60,
			"library_qty" : 500
		},
		"rna" : {
			"rna_qty" : 200,
			"library_qty" : 500
		}
	}

	ngs_qc_key = {
		"dna" : {
			"cleandata_q30_num" : "DNA样本Q30",
			"cover_ratio_num" : "覆盖度",
			"uni20_uniq_hot_num" : "组织均一性（热点区域）",
			"uni20_uniq_nonhot_num" : "组织均一性（非热点区域）",
			"depth_mean_uniq_hot_num" : "组织平均有效深度（热点区域）",
			"depth_mean_uniq_nonhot_num" : "组织平均有效深度（非热点区域）"
		},
		"rna_v1" : {
			"cleandata_q30_num" : "RNA样本Q30",
			"read_ontarget_ratio_num" : "捕获效率",
			"mapping_ratio_num" : "比对率",
			"end2sense_ratio_num" : "RNA链特异性指标",
			"totaldata" : "下机数据量（G）"
		},
		"rna_v2" : {
			"cleandata_q30_num" : "RNA样本Q30",
			"mapping_ratio_num" : "比对率",
			"end2sense_ratio_num" : "RNA链特异性指标",
			"effectivereads_num" : "有效Reads数"
		}
	}

	lib_qc_key = {
		"dna" : {
			"dna_qty" : "DNA总量（ng）",
			"break_dna_qty" : "片段化DNA总量（ng）",
			"library_qty" : "DNA预文库总量（ng）"
		},
		"rna" : {
			"rna_qty" : "RNA总量（ng）",
			"library_qty" : "RNA预文库总量（ng）"
		}
	}

	ngs_dna_qc_stand = ngs_qc_standard.get("dna_v1") if prod_type == "v1" else ngs_qc_standard.get("dna_v2")
	ngs_rna_qc_stand = ngs_qc_standard.get("rna_v1") if prod_type == "v1" else ngs_qc_standard.get("rna_v2")
	lib_dna_qc_stand = lib_qc_standard.get("dna")
	lib_rna_qc_stand = lib_qc_standard.get("rna")
	ngs_dna_qc_key = ngs_qc_key.get("dna")
	ngs_rna_qc_key = ngs_qc_key.get("rna_v1") if prod_type == "v1" else ngs_qc_key.get("rna_v2")
	lib_dna_qc_key = lib_qc_key.get("dna")
	lib_rna_qc_key = lib_qc_key.get("rna")

	ngs_dna_qc = ngs_qc["dna_data_qc"] if "dna_data_qc" in ngs_qc.keys() and ngs_qc["dna_data_qc"] else {}
	ngs_rna_qc = ngs_qc["rna_data_qc"] if "rna_data_qc" in ngs_qc.keys() and ngs_qc["rna_data_qc"] else {}
	lib_dna_qc = lib_qc["lib_dna_qc"] if "lib_dna_qc" in lib_qc.keys() and lib_qc["lib_dna_qc"] else {}
	lib_rna_qc = lib_qc["rna_lib_qc"] if "rna_lib_qc" in lib_qc.keys() and lib_qc["rna_lib_qc"] else {}

	Fail_item = []
	# DNA 湿实验质控
	for item in lib_dna_qc_stand.keys():
		if item in lib_dna_qc.keys() and lib_dna_qc[item] and is_number(lib_dna_qc[item]) and float(lib_dna_qc[item]) < lib_dna_qc_stand.get(item):
			Fail_item.append(lib_dna_qc_key.get(item))
	# DNA NGS质控
	for item in ngs_dna_qc_stand.keys():
		if item in ngs_dna_qc.keys() and ngs_dna_qc[item] and is_number(ngs_dna_qc[item]) and float(ngs_dna_qc[item]) < ngs_dna_qc_stand.get(item):
			Fail_item.append(ngs_dna_qc_key.get(item))
	# RNA 湿实验质控
	for item in lib_rna_qc_stand.keys():
		if item in lib_rna_qc.keys() and lib_rna_qc[item] and is_number(lib_rna_qc[item]) and float(lib_rna_qc[item]) < lib_rna_qc_stand.get(item):
			Fail_item.append(lib_rna_qc_key.get(item))
	# RNA NGS质控
	for item in ngs_rna_qc_stand.keys():
		if item == "totaldata":
			if item in ngs_rna_qc.keys() and ngs_rna_qc[item]:
				total_data = ngs_rna_qc[item].replace("G", "")
				if is_number(total_data):
					if float(total_data) < ngs_rna_qc_stand.get(item):
						Fail_item.append(ngs_rna_qc_key.get(item))
				else:
					Fail_item.append(ngs_rna_qc_key.get(item))
		else:
			if item in ngs_rna_qc.keys() and ngs_rna_qc[item] and is_number(ngs_rna_qc[item]) and float(ngs_rna_qc[item]) < ngs_rna_qc_stand.get(item):
				Fail_item.append(ngs_rna_qc_key.get(item))
	return "，".join(Fail_item)
jinja2.filters.FILTERS["bds_mp_qc"] = bds_mp_qc

# 北大三-MP-筛选出不合格的质控项-2025.06.17-QC阈值有变
def bds_mp_qc_v2(info):
	lib_qc = info[0]
	ngs_qc = info[1]
	prod_type = info[2]
	ngs_qc_standard = {
		"dna_v1" : {
			"cleandata_q30_num" : 0.75,
			"cover_ratio_num" : 0.95,
			"uni20_uniq_hot_num" : 0.9,
			"uni20_uniq_nonhot_num" : 0.8,
			"depth_mean_uniq_hot_num" : 1000,
			"depth_mean_uniq_nonhot_num" : 500
		},
		"rna_v1" : {
			"cleandata_q30_num" : 0.75,
			"read_ontarget_ratio_num" : 0.5,
			"mapping_ratio_num" : 0.8,
			"end2sense_ratio_num" : 0.9,
			"totaldata" : 2
		},
		"dna_v2" : {
			"cleandata_q30_num" : 0.75,
			"cover_ratio_num" : 0.95,
			"uni20_uniq_hot_num" : 0.9,
			"uni20_uniq_nonhot_num" : 0.8,
			"depth_mean_uniq_hot_num" : 800,
			"depth_mean_uniq_nonhot_num" : 400
		},
		"rna_v2" : {
			"cleandata_q30_num" : 0.75,
			"mapping_ratio_num" : 0.8,
			"end2sense_ratio_num" : 0.9,
			"effectivereads_num" : 3000000
		}
	}

	lib_qc_standard = {
		"dna" : {
			"dna_qty" : 150,
			"break_dna_qty" : 60,
			"library_qty" : 500
		},
		"rna" : {
			"rna_qty" : 200,
			"library_qty" : 500
		}
	}

	ngs_qc_key = {
		"dna" : {
			"cleandata_q30_num" : "DNA样本Q30",
			"cover_ratio_num" : "覆盖度",
			"uni20_uniq_hot_num" : "组织均一性（热点区域）",
			"uni20_uniq_nonhot_num" : "组织均一性（非热点区域）",
			"depth_mean_uniq_hot_num" : "组织平均有效深度（热点区域）",
			"depth_mean_uniq_nonhot_num" : "组织平均有效深度（非热点区域）"
		},
		"rna_v1" : {
			"cleandata_q30_num" : "RNA样本Q30",
			"read_ontarget_ratio_num" : "捕获效率",
			"mapping_ratio_num" : "比对率",
			"end2sense_ratio_num" : "RNA链特异性指标",
			"totaldata" : "下机数据量（G）"
		},
		"rna_v2" : {
			"cleandata_q30_num" : "RNA样本Q30",
			"mapping_ratio_num" : "比对率",
			"end2sense_ratio_num" : "RNA链特异性指标",
			"effectivereads_num" : "有效Reads数"
		}
	}

	lib_qc_key = {
		"dna" : {
			"dna_qty" : "DNA总量（ng）",
			"break_dna_qty" : "片段化DNA总量（ng）",
			"library_qty" : "DNA预文库总量（ng）"
		},
		"rna" : {
			"rna_qty" : "RNA总量（ng）",
			"library_qty" : "RNA预文库总量（ng）"
		}
	}

	ngs_dna_qc_stand = ngs_qc_standard.get("dna_v1") if prod_type == "v1" else ngs_qc_standard.get("dna_v2")
	ngs_rna_qc_stand = ngs_qc_standard.get("rna_v1") if prod_type == "v1" else ngs_qc_standard.get("rna_v2")
	lib_dna_qc_stand = lib_qc_standard.get("dna")
	lib_rna_qc_stand = lib_qc_standard.get("rna")
	ngs_dna_qc_key = ngs_qc_key.get("dna")
	ngs_rna_qc_key = ngs_qc_key.get("rna_v1") if prod_type == "v1" else ngs_qc_key.get("rna_v2")
	lib_dna_qc_key = lib_qc_key.get("dna")
	lib_rna_qc_key = lib_qc_key.get("rna")

	ngs_dna_qc = ngs_qc["dna_data_qc"] if "dna_data_qc" in ngs_qc.keys() and ngs_qc["dna_data_qc"] else {}
	ngs_rna_qc = ngs_qc["rna_data_qc"] if "rna_data_qc" in ngs_qc.keys() and ngs_qc["rna_data_qc"] else {}
	lib_dna_qc = lib_qc["lib_dna_qc"] if "lib_dna_qc" in lib_qc.keys() and lib_qc["lib_dna_qc"] else {}
	lib_rna_qc = lib_qc["rna_lib_qc"] if "rna_lib_qc" in lib_qc.keys() and lib_qc["rna_lib_qc"] else {}

	Fail_item = []
	# DNA 湿实验质控
	for item in lib_dna_qc_stand.keys():
		if item in lib_dna_qc.keys() and lib_dna_qc[item] and is_number(lib_dna_qc[item]) and float(lib_dna_qc[item]) < lib_dna_qc_stand.get(item):
			Fail_item.append(lib_dna_qc_key.get(item))
	# DNA NGS质控
	for item in ngs_dna_qc_stand.keys():
		if item in ngs_dna_qc.keys() and ngs_dna_qc[item] and is_number(ngs_dna_qc[item]) and float(ngs_dna_qc[item]) < ngs_dna_qc_stand.get(item):
			Fail_item.append(ngs_dna_qc_key.get(item))
	# RNA 湿实验质控
	for item in lib_rna_qc_stand.keys():
		if item in lib_rna_qc.keys() and lib_rna_qc[item] and is_number(lib_rna_qc[item]) and float(lib_rna_qc[item]) < lib_rna_qc_stand.get(item):
			Fail_item.append(lib_rna_qc_key.get(item))
	# RNA NGS质控
	for item in ngs_rna_qc_stand.keys():
		if item == "totaldata":
			if item in ngs_rna_qc.keys() and ngs_rna_qc[item]:
				total_data = ngs_rna_qc[item].replace("G", "")
				if is_number(total_data):
					if float(total_data) < ngs_rna_qc_stand.get(item):
						Fail_item.append(ngs_rna_qc_key.get(item))
				else:
					Fail_item.append(ngs_rna_qc_key.get(item))
		else:
			if item in ngs_rna_qc.keys() and ngs_rna_qc[item] and is_number(ngs_rna_qc[item]) and float(ngs_rna_qc[item]) < ngs_rna_qc_stand.get(item):
				Fail_item.append(ngs_rna_qc_key.get(item))
	return "，".join(Fail_item)
jinja2.filters.FILTERS["bds_mp_qc_v2"] = bds_mp_qc_v2

def hrr_parp(var_list):
	gene_list = ["ATM","BARD1","BRCA1","BRCA2","BRIP1","CDK12","CHEK1","CHEK2","FANCL","PALB2",\
			  	 "RAD51B","RAD51C","RAD51D","RAD54L", "FANCA", "ATR", "MRE11", "NBN"]
	detect_gene = [var["gene_symbol"] for var in var_list if var["gene_symbol"] in gene_list]
	result = []
	# 2025.12.23:修复bug，var_list需要根据基因列表筛一下
	for var in var_list:
		if var["gene_symbol"] in gene_list:
			result.append(var)
	#result.extend(var_list)			
	for gene in set(gene_list) - set(detect_gene):
		result.append({
			"gene_symbol" : gene
		})	
	return sorted(result, key=lambda i:i["gene_symbol"])
jinja2.filters.FILTERS["hrr_parp"] = hrr_parp

# 山东区域模板，BRCA返回1-3类变异等级-2024.12.04
def brca_summary_sd_level123(var_list):
	result = []
	level_dict = {3 : "意义不明确突变", 2 : "疑似无害突变", 1 : "无害突变"}
	sort_rule = ["意义不明确突变", "疑似无害突变", "无害突变", "多态性改变"]
	for var in var_list:
		if "tag" in var.keys() and re.search("Polymorphism", var["tag"]) and var["clinic_num_g"] == 1:
			result.append("多态性改变")
		else:
			if var["clinic_num_g"] in level_dict.keys():
				result.append(level_dict.get(var["clinic_num_g"]))
	return "、".join(sorted(set(result), key=lambda i:sort_rule.index(i)))
jinja2.filters.FILTERS["brca_summary_sd_level123"] = brca_summary_sd_level123
	
# 厦门市一116临床提示排序更新-2024.12.06
# 同等级FDA/NMPA获批证据优先展示
def xmsy_evi_sort(info):
	evi_list = info[0]
	var = info[1]
	for evi in evi_list:
		if set(re.split("、", evi["refer_agency"])) & set(["FDA", "NMPA"]):
			evi["xmsy_116_agency"] = 0
		else:
			evi["xmsy_116_agency"] = 1
	# 排序-在通用的基础上增加-FDA/NMPA的证据优先展示
	# RET 融合，用药排序更新-由之前的卡博替尼、普拉替尼、赛普替尼更新为普拉提尼、赛普替尼、卡博替尼-2023.05.11
	drug_list_rule = {"普拉替尼" : 0, "塞普替尼" : 1,"卡博替尼" : 2}
	# 2024.08.09-新增ERBB2 snv/cnv药物排序，德曲妥珠单抗同等级的放在最前面
	erbb2_drug_list_rule = {"德曲妥珠单抗" : 0}
	if "bio_category" in var.keys() and var["bio_category"] in ["Sv", "PSeqRnaSv"] and (var["five_prime_gene"] == "RET" or var["three_prime_gene"] == "RET"):
		evi_list = sorted(evi_list, key = lambda i:(i["evi_conclusion_simple"], \
												 i["sense_rule"], \
												 drug_list_rule.get(i["regimen_name"], 3) , \
												 i["xmsy_116_agency"], \
												 i["regimen_name_py"].upper()))
	# 2024.08.09-新增ERBB2 snv/cnv药物排序，德曲妥珠单抗同等级的放在最前面
	elif "bio_category" in var.keys() and var["bio_category"] in ["Snvindel", "Cnv"] and var["gene_symbol"] == "ERBB2":
		evi_list = sorted(evi_list, key = lambda i:(i["evi_conclusion_simple"], \
												 i["sense_rule"], \
												 erbb2_drug_list_rule.get(i["regimen_name"], 1) , \
												 i["xmsy_116_agency"], \
												 i["regimen_name_py"].upper()))
	else:
		evi_list = sorted(evi_list, key = lambda i:(i["evi_conclusion_simple"], \
											  		i["sense_rule"], \
													i["xmsy_116_agency"], \
													i["regimen_name_py"].upper()))
		#for evi in evi_list:
		#	print (evi["regimen_name"], evi["evi_conclusion_simple"], evi["sense_rule"], evi["xmsy_116_agency"], evi["regimen_name_py"].upper())
	return evi_list
jinja2.filters.FILTERS["xmsy_evi_sort"] = xmsy_evi_sort

# 湘雅二院ptBRCA小结，需要区分体细胞/胚系、BRCA1/2-2024.12.09
def xyey_ptbrca_summary(info):
	var_list_raw = info[0]
	origin = info[1]
	var_list_g = [var for var in var_list_raw if "var_origin" in var.keys() and var["var_origin"] and var["var_origin"] == "germline" or var["type"] in ["Loss", "Gain"]]
	var_list_s = [var for var in var_list_raw if "var_origin" in var.keys() and var["var_origin"] and var["var_origin"] == "somatic"]
	for gene in set(["BRCA1", "BRCA2"]) - set([var["gene_symbol"] for var in var_list_g]):
		var_list_g.append({"gene_symbol" : gene, "xyey_detect_result" : "noFound"})
	for gene in set(["BRCA1", "BRCA2"]) - set([var["gene_symbol"] for var in var_list_s]):
		var_list_s.append({"gene_symbol" : gene, "xyey_detect_result" : "noFound"})
	if origin == "germline":
		return sorted(var_list_g, key=lambda i:i["gene_symbol"])
	else:
		return sorted(var_list_s, key=lambda i:i["gene_symbol"])
jinja2.filters.FILTERS["xyey_ptbrca_summary"] = xyey_ptbrca_summary

# 上海仁济BRCA-结果汇总-2024.12.09
def var_summary_SHRJ(var_list):
	'''
	上海仁济（单独组织、单独全血）
	检测结果需要拼接
	'''
	varInter_list = []
	for var in var_list:
		if var["type"] == "Loss":
			inter = "检测到{0}基因{1} del突变，{2}".format(var["gene_symbol"], var["value"], var["variant_desc_cn"])
		elif var["type"] == "Gain":
			inter = "检测到{0}基因{1} dup突变，{2}".format(var["gene_symbol"], var["value"], var["variant_desc_cn"])
		else:
			if var["hgvs_p"] != "p.?":
				inter = "检测到"+var["gene_symbol"]+"基因的"+var["hgvs_c"]+"（"+var["hgvs_p"]+"）突变，"+var["variant_desc_cn"]
			else:
				inter = "检测到"+var["gene_symbol"]+"基因的"+var["hgvs_c"]+"突变，"+var["variant_desc_cn"]
			if re.search("Splicing|Intronic", var["type"]):
				inter = inter[0:-1]+"，可能造成异常剪接。"
			if re.search("FrameShift|Nonsense", var["type"]):
				inter = inter[0:-1]+"，可能形成功能损伤或失活的蛋白。"
		varInter_list.append(inter)
	
	return "".join(varInter_list)
jinja2.filters.FILTERS["var_summary_SHRJ"] = var_summary_SHRJ

# 山东省立bptm plus结果汇总-2024.09.29
def sdsl_bptm_plus_summary(var_list):
	gene_list = ["BRCA1", "BRCA2", "CTNNB1", "EPCAM", "MLH1", "MSH2", "MSH6", "PMS2", "POLE", "TP53"]
	result = []
	result.extend(var_list)
	detect_gene = [var["gene_symbol"] for var in var_list]
	for gene in set(gene_list) - set(detect_gene):
		result.append({"gene_symbol" : gene})
	return sorted(result, key=lambda i:gene_list.index(i["gene_symbol"]))
jinja2.filters.FILTERS["sdsl_bptm_plus_summary"] = sdsl_bptm_plus_summary

# 吉大一HRD-列出BRCA1/2检出变异等级-2024.12.16
def jdy_hrd_brca_level_summary(var_list):
	result = []
	brca1_level_5 = [var for var in var_list if var["gene_symbol"] == "BRCA1" and var["clinic_num_g"] == 5]
	brca2_level_5 = [var for var in var_list if var["gene_symbol"] == "BRCA2" and var["clinic_num_g"] == 5]
	brca1_level_4 = [var for var in var_list if var["gene_symbol"] == "BRCA1" and var["clinic_num_g"] == 4]
	brca2_level_4 = [var for var in var_list if var["gene_symbol"] == "BRCA2" and var["clinic_num_g"] == 4]
	if brca1_level_5:
		result.append("BRCA1： 致病性变异：是")
	if brca2_level_5:
		result.append("BRCA2： 致病性变异：是")
	if brca1_level_4:
		result.append("BRCA1： 疑似致病性变异：是")
	if brca2_level_4:
		result.append("BRCA2： 疑似致病性变异：是")
	return "；".join(result)
jinja2.filters.FILTERS["jdy_hrd_brca_level_summary"] = jdy_hrd_brca_level_summary

# 吉大一HRD-输入证据，返回敏感最高等级-2024.12.16
def jdy_hrd_evi_filter(info):
	evi_list = info[0]
	evi_type = info[1]
	sense_evi = [evi for evi in evi_list if evi["clinical_significance_cn"] == "敏感"]
	sense_evi_D = [evi for evi in evi_list if evi["clinical_significance_cn"] == "敏感" and evi["evi_conclusion_simple"] == "D"]
	tmp_dict = {}
	for level in ["A", "B", "C", "D"]:
		tmp_dict["sense_"+level] = ["{0}（{1}）".format(evi["regimen_name"], evi["evi_conclusion_simple"]) for \
								   evi in sense_evi if evi["evi_conclusion_simple"] == level]
	regimen_sense_list_gss = "，".join(tmp_dict["sense_A"]) if tmp_dict["sense_A"] else \
							 "，".join(tmp_dict["sense_B"]) if tmp_dict["sense_B"] else \
							 "，".join(tmp_dict["sense_C"]) if tmp_dict["sense_C"] else \
							 []
	regimen_sense_list_hrd = "，".join(tmp_dict["sense_A"]) + "，A级推荐" if tmp_dict["sense_A"] else \
							 "，".join(tmp_dict["sense_B"]) + "，B级推荐" if tmp_dict["sense_B"] else \
							 "，".join(tmp_dict["sense_C"]) + "，C级推荐" if tmp_dict["sense_C"] else \
							 []
	if regimen_sense_list_gss:
		if evi_type == "gss":
			return regimen_sense_list_gss
		else:
			return regimen_sense_list_hrd
	else:
		if sense_evi_D:
			return "nofoundsignificant"
		else:
			return "nofoundvar"
jinja2.filters.FILTERS["jdy_hrd_evi_filter"] = jdy_hrd_evi_filter

# 吉林大学第一院HRD-返回其他基因敏感结果-2024.12.17
# 敏感药物只展示等级最高的（除了C）
def jlyy_hrd_othergene_result_sense(raw_var_list):
	result_sense = []
	regimen_D_list_sense = []
	var_list = [var for var in raw_var_list if var["gene_symbol"] not in ["BRCA1", "BRCA2"]]
	for var in var_list:
		if "Predictive" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Predictive"]:
			tmp_dict = {}
			regimen_sense_list = []
			for level in ["A", "B", "C", "D"]:
				tmp_dict["sense_"+level] = ["{0}（{1}）".format(evi["regimen_name"], evi["evi_conclusion_simple"]) for \
										   evi in var["evi_sum"]["evi_split"]["Predictive"] if evi["evi_conclusion_simple"] == level and \
										   evi["clinical_significance_cn"] == "敏感"]
			regimen_sense_list = tmp_dict["sense_A"] if tmp_dict["sense_A"] else \
								 tmp_dict["sense_B"] if tmp_dict["sense_B"] else \
								 tmp_dict["sense_C"] if tmp_dict["sense_C"] else \
								 regimen_D_list_sense.extend(tmp_dict["sense_D"])
			
			regimen_sense_list_sum = "，".join(tmp_dict["sense_A"]) + "，A级推荐" if tmp_dict["sense_A"] else \
								     "，".join(tmp_dict["sense_B"]) + "，B级推荐" if tmp_dict["sense_B"] else \
								     "，".join(tmp_dict["sense_C"]) + "，C级推荐" if tmp_dict["sense_C"] else \
								     regimen_D_list_sense.extend(tmp_dict["sense_D"])
		
			if regimen_sense_list:
				var["regimen_sum_for_jdyy"] = "，".join(regimen_sense_list)
				var["regimen_sum_for_jdyy_total"] = regimen_sense_list_sum
				result_sense.append(var)

	if not result_sense:
		if regimen_D_list_sense:
			result_sense = [{"gene_symbol" : "nofoundsignificant"}]
		else:
			result_sense = [{"gene_symbol" : "nofoundvar"}]
	return result_sense
jinja2.filters.FILTERS["jlyy_hrd_othergene_result_sense"] = jlyy_hrd_othergene_result_sense

#吉大一HRD-过滤BRCA变异-2024.12.17
def jlyy_hrd_filter_brca_var(var_list):
	return [var for var in var_list if var["gene_symbol"] not in ["BRCA1", "BRCA2"]]
jinja2.filters.FILTERS["jlyy_hrd_filter_brca_var"] = jlyy_hrd_filter_brca_var

# 北大人民MP血液-变异需要分为snvindel、cnv和sv-2024.12.18
def bdrm_gmp_split_var(info):
	var_list = info[0]
	var_type = info[1]
	return [var for var in var_list if var["bio_category"] == var_type]
jinja2.filters.FILTERS["bdrm_gmp_split_var"] = bdrm_gmp_split_var

# 广东人民tHRR结果汇总-2024.12.23
def gdrm_thrr_total_gene_sum(var_list):
	gene_list = ["BRCA1","BRCA2","AR","ATM","ATR","BARD1","BRIP1","CDH1","CDK12","CHEK1",\
				 "CHEK2","ESR1","FANCA","FANCL","HDAC2","HOXB13","MRE11","NBN","PALB2","PPP2R2A",\
				 "PTEN","RAD51B","RAD51C","RAD51D","RAD54L","STK11","TP53","BRAF","ERBB2","KRAS",\
				 "NRAS","PIK3CA"]
	detect_gene = [var["gene_symbol"] for var in var_list if var["gene_symbol"] in gene_list]
	for gene in set(gene_list) - set(detect_gene):
		var_list.append({
			"gene_symbol" : gene
		})
	return sorted(var_list, key=lambda i:gene_list.index(i["gene_symbol"]))
jinja2.filters.FILTERS["gdrm_thrr_total_gene_sum"] = gdrm_thrr_total_gene_sum

# 广东人民tHRR结果小结，仅展示I/II类变异-2024.12.24
# 基于本次送检样本，本次实验检测到X个具有明确临床意义的变异：……
def gdrm_thrr_I_II_sum(var_list):
	def get_var_info(var):
		if var["hgvs_p"] != "p.?":
			return "{0}基因({1}) {2}{3} {4}".format(var["gene_symbol"], var["transcript_primary"], \
										 			var["varInfo_XAJDY"].replace("第", ""), var["hgvs_c"], var["hgvs_p"])
		else:
			return "{0}基因({1}) {2}{3}".format(var["gene_symbol"], var["transcript_primary"], \
										 			var["varInfo_XAJDY"].replace("第", ""), var["hgvs_c"])
	level_I = [get_var_info(var) for var in var_list if var["clinic_num_s"] == 5]
	level_II = [get_var_info(var) for var in var_list if var["clinic_num_s"] == 4]
	result = []
	if level_I:
		result.append("检测到{0}个具有明确临床意义的变异：{1}。".format(len(level_I), "、".join(level_I)))
	if level_II:
		result.append("检测到{0}个具有潜在临床意义的变异：{1}。".format(len(level_II), "、".join(level_II)))
	return "基于本次送检样本，本次实验{0}".format("".join(result))
jinja2.filters.FILTERS["gdrm_thrr_I_II_sum"] = gdrm_thrr_I_II_sum

# 北大三MP-非子宫内膜癌需要展示分子分型-2024.12.26
def judge_ec_type(info):
	var = info[0]
	sample = info[1]
	msi = info[2]
	if sample["control_sample_id"]:
		var_list = var["var_somatic"]["level_I"] + var["var_somatic"]["level_II"]
	else:
		var_list = var["var_somatic"]["level_I"] + var["var_germline"]["regimen_level_I"] + \
				   var["var_somatic"]["level_II"] + var["var_germline"]["regimen_level_II"]
	pole_var = ["p.(M294R)", "p.(P286R)", "p.(Y458F)", "p.(F285_P286delinsLR)", "p.(N423_L424delinsKI)", \
			 	"p.(P286C)", "p.(P436R)", "p.(N363K)", "p.(F367S)", "p.(M295R)", \
				"p.(L424I)", "p.(L424V)", "p.(P286H)", "p.(S461P)", "p.(S297F)", \
				"p.(D368N)", "p.(V411L)", "p.(D368Y)", "p.(M444K)", "p.(A456P)", "p.(S459F)"]
	pole_list = [var for var in var_list if var["bio_category"] == "Snvindel" and var["gene_symbol"] == "POLE" and var["hgvs_p"] in pole_var]
	tp53_list = [var for var in var_list if var["bio_category"] == "Snvindel" and var["gene_symbol"] == "TP53"]
	if pole_list:
		ec_type = "POLE突变亚型（POLE mut）"
	elif msi["var_id"] == "MSI-H":
		ec_type = "错配修复功能缺陷亚型（MMRd）"
	elif tp53_list:
		ec_type = "p53异常亚型（p53 abnormality，p53abn）"
	else:
		ec_type = "非特异性分子特征亚型（Non-specific molecular profile，NSMP）"
	return ec_type
jinja2.filters.FILTERS["judge_ec_type"] = judge_ec_type

# 北大三MP-非子宫内膜癌需要展示分子分型-2025.06.17
# TP53增加PHd
def judge_ec_type_v2(info):
	var = info[0]
	sample = info[1]
	msi = info[2]
	if sample["control_sample_id"]:
		var_list = var["var_somatic"]["level_I"] + var["var_somatic"]["level_II"]
	else:
		var_list = var["var_somatic"]["level_I"] + var["var_germline"]["regimen_level_I"] + \
				   var["var_somatic"]["level_II"] + var["var_germline"]["regimen_level_II"]
	pole_var = ["p.(M294R)", "p.(P286R)", "p.(Y458F)", "p.(F285_P286delinsLR)", "p.(N423_L424delinsKI)", \
			 	"p.(P286C)", "p.(P436R)", "p.(N363K)", "p.(F367S)", "p.(M295R)", \
				"p.(L424I)", "p.(L424V)", "p.(P286H)", "p.(S461P)", "p.(S297F)", \
				"p.(D368N)", "p.(V411L)", "p.(D368Y)", "p.(M444K)", "p.(A456P)", "p.(S459F)"]
	pole_list = [var for var in var_list if var["bio_category"] == "Snvindel" and var["gene_symbol"] == "POLE" and var["hgvs_p"] in pole_var]
	tp53_list = [var for var in var_list if var["bio_category"] in ["Snvindel", "PHd"] and var["gene_symbol"] == "TP53"]
	if pole_list:
		ec_type = "POLE突变亚型（POLE mut）"
	elif msi["var_id"] == "MSI-H":
		ec_type = "错配修复功能缺陷亚型（MMRd）"
	elif tp53_list:
		ec_type = "p53异常亚型（p53 abnormality，p53abn）"
	else:
		ec_type = "非特异性分子特征亚型（Non-specific molecular profile，NSMP）"
	return ec_type
jinja2.filters.FILTERS["judge_ec_type_v2"] = judge_ec_type_v2

# 北大三MP首页和基本信息页病理号展示-2025.01.03
def bds_mp_pathological_id(sample):
	# 进院：若对照样本类型不含有“血”，则对病理号“XXX1(肿瘤)/XXX2(对照)”进行拆分处理，首页展示XXX1，第二页展示XXX1(肿瘤)/XXX2(对照)
	#      若对照样本类型含有“血”，则首页和第二页展示完整的病理号
	# 临检：mark.样本类型通过"、"进行拆分
	#      判断对照是否含有血，规则同进院
	
	# 对照样本类型
	pathological_first = ""
	if sample["control_sample_id"]:
		# 对照样本类型处理
		control_sample_type = sample["control_sample_type"] if sample["report_module_type"] == "hospital" else \
							  sample["ZJFB_mark_dict"]["样本类型"] if "样本类型" in sample["ZJFB_mark_dict"].keys() and sample["ZJFB_mark_dict"]["样本类型"] else \
							  ""
		#print (control_sample_type)
		# 首页病理号处理
		patho_list = re.split("/", sample["pathological_id"])
		#print ("patho_list", patho_list)
		#pathological_first = patho_list[0].replace("(肿瘤)", "") if len(patho_list) >= 1 else sample["pathological_id"]

		# 首页病理号展示判定
		if "血" not in control_sample_type:
			#print ("血in control")
			pathological_first = patho_list[0].replace("(肿瘤)", "") if len(patho_list) >= 1 else sample["pathological_id"]
		else:
			pathological_first = sample["pathological_id"]
	else:
		pathological_first = sample["pathological_id"]
	#print (pathological_first)
	return pathological_first
jinja2.filters.FILTERS["bds_mp_pathological_id"] = bds_mp_pathological_id

# 北大三MP样本类型去重-2025.01.03
def bds_sample_type_redup(sample):
	sample_list = []
	if sample["report_module_type"] == "hospital":
		if sample["sample_type"] and sample["sample_type"] not in sample_list:
			sample_list.append(sample["sample_type"])
		if sample["control_sample_type"] and sample["control_sample_type"] not in sample_list:
			sample_list.append(sample["control_sample_type"])
	else:
		mark_sample = sample["ZJFB_mark_dict"]["样本类型"] if "样本类型" in sample["ZJFB_mark_dict"].keys() and sample["ZJFB_mark_dict"]["样本类型"] else ""
		mark_sample_split = re.split("、", mark_sample)
		for i in mark_sample_split:
			if i not in sample_list:
				sample_list.append(i)
	return "、".join(sample_list)
jinja2.filters.FILTERS["bds_sample_type_redup"] = bds_sample_type_redup

# 青岛附属BPTM Plus组织肿瘤发生发展相关变异需要分为林奇相关和其他-2025.01.07
def qdfs_bptm_lyn(info):
	var_list = info[0]
	var_type = info[1]
	lyn_gene = ["MLH1", "MSH2", "MSH6", "PMS2", "EPCAM"]
	in_lyn_var = []
	out_lyn_var = []
	for var in var_list:
		if var["gene_symbol"] in lyn_gene and float(var["freq"]) >= 0.2:
			in_lyn_var.append(var)
		else:
			out_lyn_var.append(var)
	if var_type == "in":
		return in_lyn_var
	else:
		return out_lyn_var
jinja2.filters.FILTERS["qdfs_bptm_lyn"] = qdfs_bptm_lyn

# 同济gBRCA，靶向药物汇总-2025.01.14
def tj_gbrca_regimen_sum(var_list):
	result = []
	for var in var_list:
		if var["evi_sum"] and "evi_split" in var["evi_sum"].keys() and var["evi_sum"]["evi_split"] and \
			"Predictive" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Predictive"]:
			for evi in var["evi_sum"]["evi_split"]["Predictive"]:
				regimen_str = "{0}（{1}，{2}级）".format(evi["regimen_name"], evi["clinical_significance_cn"], evi["evi_conclusion_simple"])
				if regimen_str not in result:
					result.append(regimen_str)
	return "、".join(result)
jinja2.filters.FILTERS["tj_gbrca_regimen_sum"] = tj_gbrca_regimen_sum		

# 同济gBRCA，靶向药物汇总-2025.03.04
# 增加批准机构，仅A级展示
def tj_gbrca_regimen_sum_v2(var_list):
	result = []
	for var in var_list:
		if var["evi_sum"] and "evi_split" in var["evi_sum"].keys() and var["evi_sum"]["evi_split"] and \
			"Predictive" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Predictive"]:
			for evi in var["evi_sum"]["evi_split"]["Predictive"]:
				if "regimen_refer_agency" in evi.keys() and evi["regimen_refer_agency"] and evi["evi_conclusion_simple"] == "A":
					regimen_refer_agency = "/".join(re.split(",", evi["regimen_refer_agency"]))
					regimen_str = "{0}（{1}，{2}级，{3}）".format(evi["regimen_name"], evi["clinical_significance_cn"], evi["evi_conclusion_simple"], regimen_refer_agency)
				else:
					regimen_str = "{0}（{1}，{2}级，-）".format(evi["regimen_name"], evi["clinical_significance_cn"], evi["evi_conclusion_simple"])
				if regimen_str not in result:
					result.append(regimen_str)
	return "、".join(result)
jinja2.filters.FILTERS["tj_gbrca_regimen_sum_v2"] = tj_gbrca_regimen_sum_v2	

# 同济gBRCA，靶向药物汇总-2025.04.23
# 增加批准机构，仅A级展示
def tj_gbrca_regimen_sum_v3(var_list):
	result = []
	for var in var_list:
		if var["evi_sum"] and "evi_split" in var["evi_sum"].keys() and var["evi_sum"]["evi_split"] and \
			"Predictive" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Predictive"]:
			for evi in var["evi_sum"]["evi_split"]["Predictive"]:
				if "regimen_refer_agency" in evi.keys() and evi["regimen_refer_agency"] and evi["evi_conclusion_simple"] == "A":
					regimen_refer_agency = "/".join(re.split(",", evi["regimen_refer_agency"]))
					regimen_str = "{0}（{1}，{2}级，{3}）".format(evi["regimen_name"], evi["clinical_significance_cn"], evi["evi_conclusion_simple"], regimen_refer_agency)
				else:
					regimen_str = "{0}（{1}，{2}级）".format(evi["regimen_name"], evi["clinical_significance_cn"], evi["evi_conclusion_simple"])
				if regimen_str not in result:
					result.append(regimen_str)
	return "、".join(result)
jinja2.filters.FILTERS["tj_gbrca_regimen_sum_v3"] = tj_gbrca_regimen_sum_v3

# 南昌附一gBPTM plus展示变异小结-2025.01.14
# 包含snvindel、MLPA/CNV变异
def ncfy_gbptm_var_sum(var_list):
	result = []
	for var in var_list:
		if var["type"] == "Loss":
			result.append(var["gene_symbol"] + " " + var["value"] +" del")
		elif var["type"] == "Gain":
			result.append(var["gene_symbol"] + " " + var["value"] +" dup")
		else:
			if var["hgvs_p"] != "p.?":
				result.append(var["gene_symbol"]+" "+var["hgvs_c"]+" "+var["hgvs_p"])
			else:
				result.append(var["gene_symbol"]+" "+var["hgvs_c"])
	return ", ".join(result)
jinja2.filters.FILTERS["ncfy_gbptm_var_sum"] = ncfy_gbptm_var_sum

# 福建协和tBRCA、tLC10临床提示排序更新-2025.01.15
# 同等级NMPA获批证据优先展示，FDA其次，再展示其他（获批字段看regimen_refer_agency）
# 再按首字母排
def get_regimen_first_str(regimen_name):
    regimen_first_str = []
    for regimen in re.split("\+", regimen_name):
        if "\u4e00" <= regimen <= "\u9fff":
            _str = []
            for i in regimen:
                if "\u4e00" <= i <= "\u9fff":
                    _str.append("".join(chain.from_iterable(pinyin(i, Style.TONE3)))[0].upper())
                else:
                    _str.append(i)
            regimen_first_str.extend(_str)
        else:
            regimen_first_str.append(regimen[0].upper())
    return "".join(regimen_first_str)

def fjxh_evi_sort(evi_list):
	for evi in evi_list:
		regimen_refer_agency = evi["regimen_refer_agency"] if "regimen_refer_agency" in evi.keys() and evi["regimen_refer_agency"] else ""
		if set(re.split(",", regimen_refer_agency)) & set(["NMPA"]):
			evi["fjxh_agency"] = 0
		elif set(re.split(",", regimen_refer_agency)) & set(["FDA"]):
			evi["fjxh_agency"] = 1
		else:
			evi["fjxh_agency"] = 2
		evi["regimen_name_FJXH_str"] = get_regimen_first_str(evi["regimen_name"])
	evi_list = sorted(evi_list, key = lambda i:(i["evi_conclusion_simple"], \
											  	i["sense_rule"], \
												i["fjxh_agency"], \
												i["regimen_name_FJXH_str"]))
	#for evi in evi_list:
	#	print (evi["regimen_name"], evi["evi_conclusion_simple"], evi["sense_rule"], evi["regimen_refer_agency"] if "regimen_refer_agency" in evi.keys() else "", evi["regimen_name_FJXH_str"])
	return evi_list
jinja2.filters.FILTERS["fjxh_evi_sort"] = fjxh_evi_sort

# 福建协和-tLC10相同证据描述的合并治疗方案展示-2025.01.15
def fjxh_lc10_merge_Predictive_evi(evi_list):
	tmp_dict = {}
	for evi in evi_list:
		tmp_dict.setdefault(evi["evi_interpretation"], [])
		tmp_dict[evi["evi_interpretation"]].append({"regimen_name" : evi["regimen_name"]})
	merge_result = []
	for k, v  in tmp_dict.items():
		merge_result.append(
			{
				"regimen_name" : "、".join([i["regimen_name"] for i in v]),
				"evi_interpretation" : k,
			}
		)
	return merge_result
jinja2.filters.FILTERS["fjxh_lc10_merge_Predictive_evi"] = fjxh_lc10_merge_Predictive_evi

# 2. PARP结果汇总表，根据产品（HRR和HRD）选择对应展示列表-兼容MLPA/CNV-2025.01.15
def choose_parp_list_v2(info):
	hrr_list_mlpa = info[0]["cdx"]["format5_forHRR_for_new_vesion_mlpa"]
	hrr_list_gcnv = info[0]["cdx"]["format5_forHRR_for_new_vesion_gcnv"]
	hrd_list = info[0]["cdx"]["format4_forHRDC_for_new_vesion"]
	sample = info[1]
	mlpa_or_cnv = info[2]
	if sample["prod_names"] in ["HRR（全血）", "HRR（组织）", "HRR（组织 全血）"]:
		if mlpa_or_cnv == "mlpa":
			return hrr_list_mlpa
		else:
			return hrr_list_gcnv
	elif sample["prod_names"] in ["HRD Complete（组织）"]:
		return hrd_list
jinja2.filters.FILTERS["choose_parp_list_v2"] = choose_parp_list_v2

# 2. PARP结果汇总表，根据产品（HRR和HRD）选择对应展示列表-兼容MLPA/CNV-2025.01.15
# 2026.03.16-新增HRR 34基因，增加MLH1基因
def choose_parp_list_v3(info):
	hrr_list_mlpa = info[0]["cdx"]["format5_forHRR_for_new_vesion_mlpa_34"]
	hrr_list_gcnv = info[0]["cdx"]["format5_forHRR_for_new_vesion_gcnv_34"]
	hrd_list = info[0]["cdx"]["format4_forHRDC_for_new_vesion"]
	sample = info[1]
	mlpa_or_cnv = info[2]
	if sample["prod_names"] in ["HRR（全血）", "HRR（组织）", "HRR（组织 全血）"]:
		if mlpa_or_cnv == "mlpa":
			return hrr_list_mlpa
		else:
			return hrr_list_gcnv
	elif sample["prod_names"] in ["HRD Complete（组织）"]:
		return hrd_list
jinja2.filters.FILTERS["choose_parp_list_v3"] = choose_parp_list_v3

# 复旦中山gHRR-列出3类变异基因-2025.01.17
def fdzs_ghrr_var3_gene(var_list):
	result = []
	for var in var_list:
		if var["gene_symbol"] not in result:
			result.append(var["gene_symbol"])
	return "、".join(result)
jinja2.filters.FILTERS["fdzs_ghrr_var3_gene"] = fdzs_ghrr_var3_gene

# 甘肃武威gHRR小结，需要展示致病/疑似致病变异-兼容MLPA+CNV-2025.01.17
def gsww_ghrr_var_sum_g_v2(var_list):
	result = []
	for var in var_list:
		if var["bio_category"] == "PMLPA":
			if var["type"] == "Loss":
				result.append(var["gene_symbol"]+" "+var["value"]+" del")
			elif var["type"] == "Gain":
				result.append(var["gene_symbol"]+" "+var["value"]+" dup")
		elif var["bio_category"] == "Cnv" and var["gene_symbol"] in ["BRCA1", "BRCA2"]:
			if var["type"] == "Loss":
				result.append(var["gene_symbol"]+" "+var["value"]+" del")
			elif var["type"] == "Gain":
				result.append(var["gene_symbol"]+" "+var["value"]+" dup")
		elif var["bio_category"] == "Snvindel":
			if var["hgvs_p"] != "p.?":
				result.append(var["gene_symbol"]+" "+var["hgvs_p"])
			else:
				result.append(var["gene_symbol"]+" "+var["hgvs_c"])
	return ", ".join(result)
jinja2.filters.FILTERS["gsww_ghrr_var_sum_g_v2"] = gsww_ghrr_var_sum_g_v2

# 广东省人民gHRR结果小结-2025.01.19
# 基于本次送检样本，检测到【XXX基因gene_region hgvs_c hgvs_p致病/疑似致病性变异】。
def gdrm_ghrr_summary(var_list):
	clinic_num_g_stran = {5 : "致病性变异", 4 : "疑似致病性变异"}
	result = []
	for var in var_list:
		if var["type"] == "Loss":
			result.append("{0}基因{1} del {2}".format(var["gene_symbol"], var["value"], clinic_num_g_stran.get(var["clinic_num_g"]))) 
		elif var["type"] == "Gain":
			result.append("{0}基因{1} dup {2}".format(var["gene_symbol"], var["value"], clinic_num_g_stran.get(var["clinic_num_g"]))) 
		else:
			if var["hgvs_p"] != "p.?":
				result.append("{0}基因{1} {2} {3}{4}".format(var["gene_symbol"], var["gene_region"], var["hgvs_c"], var["hgvs_p"],clinic_num_g_stran.get(var["clinic_num_g"])))
			else:
				result.append("{0}基因{1} {2}{3}".format(var["gene_symbol"], var["gene_region"], var["hgvs_c"],clinic_num_g_stran.get(var["clinic_num_g"])))
	return ";".join(result)
jinja2.filters.FILTERS["gdrm_ghrr_summary"] = gdrm_ghrr_summary

# 适用于广东人民gHRR，展示胚系5/4/3类变异-2025.01.19
def gdrm_ghrr_summary_all(var_list):
	gene_list = ["BRCA1","BRCA2","AR","ATM","ATR","BARD1","BRIP1","CDH1","CDK12","CHEK1",\
				 "CHEK2","ESR1","FANCA","FANCL","HDAC2","HOXB13","MRE11","NBN","PALB2","PPP2R2A",\
				 "PTEN","RAD51B","RAD51C","RAD51D","RAD54L","STK11","TP53","BRAF","ERBB2","KRAS",\
				 "NRAS","PIK3CA"]
	detect_gene = [var["gene_symbol"] for var in var_list if var["gene_symbol"] in gene_list]
	for gene in set(gene_list) - set(detect_gene):
		var_list.append({
			"gene_symbol" : gene
		})	
	return sorted(var_list, key=lambda i:gene_list.index(i["gene_symbol"]))
jinja2.filters.FILTERS["gdrm_ghrr_summary_all"] = gdrm_ghrr_summary_all

# 适用于山东齐鲁gHRR，展示胚系5/4/3类变异-2025.01.19
def sdql_ghrr_summary_all(var_list):
	gene_list_all = ["BRCA1", "BRCA2", "AR", "ATM", "ATR", "BARD1", "BRIP1", "CDH1", "CDK12", "CHEK1", "CHEK2", "ESR1", \
					 "FANCA", "FANCL", "HDAC2", "HOXB13", "MRE11", "NBN", "PALB2", "PPP2R2A", "PTEN", "RAD51B", "RAD51C", \
					 "RAD51D", "RAD54L", "STK11", "TP53", "BRAF", "ERBB2", "KRAS", "NRAS", "PIK3CA"]
	gene_list_appr = ["BRCA1", "BRCA2", "ATM", "BARD1", "BRIP1", "CDK12", "CHEK1", "CHEK2", "FANCL", "PALB2", "RAD51B", "RAD51C", "RAD51D", "RAD54L"]
	detect_gene = [var["gene_symbol"] for var in var_list]
	for gene in set(gene_list_all) - set(detect_gene):
		var_list.append({"gene_symbol" : gene})
	for i in var_list:
		if i["gene_symbol"] in  gene_list_appr:
			i["appr"] = "T"
	return sorted(var_list, key = lambda i:i["gene_symbol"])
jinja2.filters.FILTERS["sdql_ghrr_summary_all"] = sdql_ghrr_summary_all

# 广东省人民ptHRR结果小结-2025.01.22
# 胚系结果从上面的函数中获取，这边主要实现体细胞变异
def gdrm_pthrr_summary(info):
	result = []
	germline_list = info[0]
	somatic_list = info[1]
	germline_str = gdrm_ghrr_summary(germline_list)
	if germline_str:
		result.append(germline_str)
	clinic_num_s_stran = {5 : "I类变异", 4 : "II类变异"}
	somatic_var_list = []
	for var in somatic_list:
		if var["hgvs_p"] != "p.?":
			somatic_var_list.append("{0}基因{1} {2} {3} {4}".format(var["gene_symbol"], var["gene_region"], var["hgvs_c"], var["hgvs_p"],clinic_num_s_stran.get(var["clinic_num_s"])))
		else:
			somatic_var_list.append("{0}基因{1} {2} {3}".format(var["gene_symbol"], var["gene_region"], var["hgvs_c"],clinic_num_s_stran.get(var["clinic_num_s"])))
	if somatic_var_list:
		result.append(";".join(somatic_var_list))

	return ";".join(result)
jinja2.filters.FILTERS["gdrm_pthrr_summary"] = gdrm_pthrr_summary

# 福建附一CP200新增参考文献-2025.01.26
def refer_nccn_fjfy_cp200(sample):
	refer = {
		"非小细胞肺癌" : [
			"中国非小细胞肺癌患者EGFR+T790M基因突变检测专家共识（2018年版）",
			"二代测序技术在NSCLC中的临床应用中国专家共识（2020版）",
			"非小细胞肺癌分子病理检测临床实践指南（2021版）",
			"中国非小细胞肺癌RET基因融合临床检测专家共识（2021年版）",
			"非小细胞肺癌MET临床检测中国专家共识（2022年版）",
			"非小细胞肺癌细针穿刺细胞学标本基因检测专家共识（2022年版）",
			"非小细胞肺癌恶性浆膜腔积液分子病理检测中国专家共识（2022年版）",
			"非小细胞肺癌融合基因检测临床实践中国专家共识（2023年版）",
			"中国晚期非小细胞肺癌BRAF突变诊疗专家共识（2023年版）",
			"非小细胞肺癌表皮生长因子受体20号外显子插入突变检测临床实践中国专家共识（2024版）"
		],
		"结直肠癌" : [
			"结直肠癌及其他相关实体瘤微卫星不稳定性检测中国专家共识（2019年版）",
			"结直肠癌分子标志物临床检测中国专家共识（2021年版）",
			"结直肠癌分子检测高通量测序中国专家共识（2021年版）",
			"二代测序技术在消化系统肿瘤临床应用的中国专家共识（2024版）"
		],
		"卵巢癌" : [
			"卵巢上皮性癌BRCA基因检测的中国专家讨论（2017年版）",
			"基于下一代测序技术的BRCA1_2基因检测指南(2019年版)",
			"上皮性卵巢癌PARP抑制剂相关生物标志物检测的中国专家共识（2020年版）"
		],
		"前列腺癌" : [
			"中国前列腺癌患者基因检测专家共识（2020年版）",
			"前列腺癌同源重组修复基因检测及变异解读专家共识（2022年版）"
		],
		"乳腺癌" : [
			"中国乳腺癌患者BRCA基因检测与临床应用专家共识（2018年版）",
			"复发/转移性乳腺癌标志物临床应用专家共识（2019年版）",
			"基于靶标指导乳腺癌精准治疗标志物临床应用专家共识（2022年版）",
			"中国晚期三阴性乳腺癌临床诊疗指南（2024版）",
			"乳腺癌HER2检测指南（2024版）"
		],
		"胃癌" : [
			"胃癌HER2检测指南（2016年版）",
			"HER2阳性晚期胃癌分子靶向治疗的中国专家共识（2016年版）",
			"胃癌高通量测序临床应用中国专家共识（2023年版）"
		],
		"胃肠间质瘤" : [
			"胃肠间质瘤基因检测与临床应用的中国专家共识（2021版）"
		],
		"子宫内膜癌" : [
			"子宫内膜癌分子检测中国专家共识（2021版）",
			"子宫内膜癌分子分型临床应用中国专家共识（2024年版）"
		],
		"甲状腺癌" : [
			"甲状腺癌RET基因检测与临床应用专家共识（2021版）",
			"甲状腺肿瘤高通量测序技术基因检测报告模板中国专家共识（2023版）",
			"RET基因变异甲状腺癌诊疗中国专家共识(2024版)"
		],
		"胰腺癌" : [
			"胰腺癌诊疗指南（2022年版)",
			"早期胰腺癌分子诊断专家共识（2023年版）"
		]
	}
	result = []
	if "实体瘤" not in sample["tumor_names_cn"]:
		set_list = set(sample["tumor_list"]) & set([i for i in refer.keys()])
		for tumor in set_list:
			result.extend(refer[tumor])
	return result
jinja2.filters.FILTERS["refer_nccn_fjfy_cp200"] = refer_nccn_fjfy_cp200

# 华西厦门116-临床提示-按厦门市一的要求排序，并且处理为一行展示-2025.02.11
def hxxm_116_regimen_sum(a):
	result = []
	evi_sum = a["evi_sum"]
	if "evi_split" in evi_sum.keys():
		if "Predictive" in evi_sum["evi_split"].keys():
			for i in xmsy_evi_sort([evi_sum["evi_split"]["Predictive"], a]):
				result.append("{0}（{1}，{2}级）".format(i["regimen_name"], i["clinical_significance_cn"], i["evi_conclusion_simple"]))
		if "Prognostic" in evi_sum["evi_split"].keys():
			for i in xmsy_evi_sort([evi_sum["evi_split"]["Prognostic"], a]):
				result.append("{0}（{1}，{2}级）".format("预后"+i["clinical_significance_cn"], " / ", i["evi_conclusion_simple"]))
		if "Diagnostic" in evi_sum["evi_split"].keys():
			for i in xmsy_evi_sort([evi_sum["evi_split"]["Diagnostic"], a]):
				result.append("{0}（{1}，{2}级）".format("辅助诊断", " / ", i["evi_conclusion_simple"]))
	return "、".join(result) if result else "-"
jinja2.filters.FILTERS["hxxm_116_regimen_sum"] = hxxm_116_regimen_sum 

# CP200 IO中融合突变（ALK）需要写为gene1-gene2的格式-2025.02.13
def cp200_alk_io(var_list):
	result = []
	for var in var_list:
		if "融合" in var:
			gene1 = re.split(":", var)[0]
			gene2 = re.split(":", re.split("-", var)[-1])[0]
			if "{0}-{1}融合".format(gene1, gene2) not in result:
				result.append("{0}-{1}融合".format(gene1, gene2))
		else:
			result.append(var)
	return result
jinja2.filters.FILTERS["cp200_alk_io"] = cp200_alk_io

# 云肿BPTM Plus组织-2025.02.19
# 判断是否检出林奇相关基因4/5类变异
def ynzl_bptmplue_judge_lyn(var_list):
	gene_list = ["EPCAM", "MLH1", "MSH2", "MSH6", "PMS2"]
	filter_var = [var for var in var_list if var["gene_symbol"] in gene_list and var["clinic_num_g"] in [4, 5]]
	if filter_var:
		return True
	else:
		return False
jinja2.filters.FILTERS["ynzl_bptmplue_judge_lyn"] = ynzl_bptmplue_judge_lyn

# 国际部HRD Focus-2025.02.20
# BRCA变异等级看致病性（胚系规则）
def hrd_focus_ID(info):
	gene = info[0]
	level = info[1]
	var = info[2]
	var_list = []
	for i in ["level_I", "level_II", "level_onco_nodrug", "level_III"]:
		if var["var_somatic"][i]:
			var_list.extend(var["var_somatic"][i])
	return [var for var in var_list if var["gene_symbol"] == gene and var["clinic_num_g"] == level]
jinja2.filters.FILTERS["hrd_focus_ID"] = hrd_focus_ID

# 统计MP免疫正负相关变异数-2025.02.21
def mp_io_count(var_dict):
	num = 0
	for gene, var_list in var_dict.items():
		for var in var_list:
			num += 1
	return num
jinja2.filters.FILTERS["mp_io_count"] = mp_io_count

# 获取描述中的PMID号
# 示例：……（PMID：111, 222, 333）……
def getPMID_from_inter(inter):
	pmid_list = []
	mat = re.compile(r"\（PMID.*?\）")
	for i in mat.findall(str(inter)):
		if re.search(":|: |：|： ", i):
			pmid = (re.split(":|: |：|： ", i))[1].replace(" ", "").replace("）", "")
			tmp_list = re.split(",", pmid)
		pmid_list.extend(tmp_list)
	return pmid_list

# 获取免疫正负相关检出变异的对应参考文献-2025.02.21
def kjdyfs_io_refer(info):
	inter_list = info[0]
	refer_dict = info[1]
	pmid_list = []
	for i in inter_list:
		pmid_list.extend(getPMID_from_inter(i["inter"]))

	refer = []
	for pmid in pmid_list:
		if pmid in refer_dict.keys() and refer_dict[pmid] not in refer:
			refer.append(refer_dict[pmid])
	return refer
jinja2.filters.FILTERS["kjdyfs_io_refer"] = kjdyfs_io_refer

# 山东省立MP-拆分出指定基因（单个）变异结果-LC10+BRCA-2025.02.25
def sdsl_mp_single_gene_var_list(info):
	var_list = info[0]
	gene = info[1]
	# 指定基因snvindel均返回、DNA融合只返回ALK/ROS1/RET基因的、扩增只返回MET基因、RNA融合返回ALK/ROS1/RET和MET 14跳跃
	# 2025.07.16-新增，BRCA1/2需要输出HD
	result = []
	for var in var_list:
		if gene in re.split(",", var["gene_symbol"]):
			if var["bio_category"] == "Snvindel":
				result.append(var)
			elif var["bio_category"] == "Cnv":
				if var["gene_symbol"] == "MET":
					result.append(var)
			elif var["bio_category"] == "Sv":
				if set(["ALK", "RET", "ROS1"]) & set(re.split(",", var["gene_symbol"])):
					result.append(var)
			# 2025.07.16-新增HD
			elif var["bio_category"] == "PHd":
				if var["gene_symbol"] in ["BRCA1", "BRCA2"]:
					result.append(var)
			# 2025.07.16-新增结束
			else:
				if set(["ALK", "RET", "ROS1"]) & set(re.split(",", var["gene_symbol"])):
					result.append(var)
				if var["five_prime_gene"] == "MET" and var["three_prime_gene"] == "MET":
					result.append(var)
	return result
jinja2.filters.FILTERS["sdsl_mp_single_gene_var_list"] = sdsl_mp_single_gene_var_list

# 山东省立MP-拆分出指定基因（固定列表）变异结果-LC10+BRCA-2025.02.25
def sdsl_mp_list_gene_var_list(var_list):
	# 指定基因snvindel均返回、融合只返回ALK/ROS1/RET基因的、扩增只返回MET基因
	# RNA MET 14 跳跃突变也返回吧
	gene_list = ["EGFR", "ALK", "ROS1", "RET", "BRAF", "MET", "KRAS", "NRAS", "ERBB2", "PIK3CA", "BRCA1", "BRCA2"]
	result = []
	for var in var_list:
		if set(gene_list) & set(re.split(",", var["gene_symbol"])):
			if var["bio_category"] == "Snvindel":
				result.append(var)
			elif var["bio_category"] == "Cnv":
				if var["gene_symbol"] == "MET":
					result.append(var)
			elif var["bio_category"] == "Sv":
				if set(["ALK", "RET", "ROS1"]) & set(re.split(",", var["gene_symbol"])):
					result.append(var)
			# 2025.07.16-新增HD
			elif var["bio_category"] == "PHd":
				if var["gene_symbol"] in ["BRCA1", "BRCA2"]:
					result.append(var)
			# 2025.07.16-新增结束
			else:
				if set(["ALK", "RET", "ROS1"]) & set(re.split(",", var["gene_symbol"])):
					result.append(var)
				if var["five_prime_gene"] == "MET" and var["three_prime_gene"] == "MET":
					result.append(var)
	return result
jinja2.filters.FILTERS["sdsl_mp_list_gene_var_list"] = sdsl_mp_list_gene_var_list

# 删除字符串中的%-2025.02.26
def remove_str(num):
	return num.replace("%", "")
jinja2.filters.FILTERS["remove_str"] = remove_str

# 云肿结果小结-相同等级变异、相同基因的要汇总到一起展示-2025.03.04
def ynzl_summary_var(var_list):
	'''
	ALK基因：p.G1202R、ALK扩增、EML4-ALK融合（EML4:exon13-ALK:exon20; EML4:exon13-ins39-ALK:exon20）
	'''
	# 获取基因列表
	gene_list = []
	for var in var_list:
		if var["gene_symbol"] not in gene_list:
			gene_list.append(var["gene_symbol"])
	
	result = {}
	for var in var_list:
		var_info = ""
		if var["bio_category"] == "Snvindel":
			var_info = var["hgvs_p_ZJZL"] if var["hgvs_p_ZJZL"] != "p.?" else var["hgvs_c"]
			if "judge_mergeMET" in var.keys() and var["judge_mergeMET"]:
				var_info += "（MET exon14 skipping）"
		elif var["bio_category"] == "Cnv":
			var_info = "扩增"
		else:
			if var["five_prime_gene"] == "MET" and var["three_prime_gene"] == "MET":
				var_info = "MET exon14 skipping"
			else:
				if "var_desc_merge" in var.keys() and var["var_desc_merge"]:
					var_desc_merge = "；".join(re.split("、", var["var_desc_merge"]))
				else:
					var_desc_merge = "{0}:{1}-{2}:{3}".format(var["five_prime_gene"], var["five_prime_cds"], var["three_prime_gene"], var["three_prime_cds"])
				var_info = "{0}-{1}融合（{2}）".format(var["five_prime_gene"], var["three_prime_gene"], var_desc_merge)
		if var["gene_symbol"] not in result.keys():
			result.setdefault(var["gene_symbol"], [])
		result[var["gene_symbol"]].append(var_info)

	result_list = []
	for gene, info in result.items():
		result_list.append(
			{
				"gene_symbol" : gene,
				"var_info" : "、".join(info)
			}
		)
	result_list = sorted(result_list, key=lambda i:gene_list.index(i["gene_symbol"]))
	return result_list
jinja2.filters.FILTERS["ynzl_summary_var"] = ynzl_summary_var

# 云南肿瘤详细解读部分-每个基因一行，若有用药相同的变异需要合并展示-2025.03.04
def ynzl_var_inter(var_list):
	gene_list = []
	for var in var_list:
		if var["gene_symbol"] not in gene_list:
			gene_list.append(gene_list)
	
	gene_dict = {}
	for var in var_list:
		if var["gene_symbol"] not in gene_dict.keys():
			gene_dict.setdefault(var["gene_symbol"], [])
		gene_dict[var["gene_symbol"]].append(var)

	# 格式转化
	# [{gene_symbol:XX,"gene_function":XX,var_info:[1,2,3]}]
	result = []
	for gene, info in gene_dict.items():
		###1. 获取基因介绍
		gene_function = []
		if "," in gene:
			for var in info:
				if var["three_prime_gene_function"] not in gene_function:
					gene_function.append(var["three_prime_gene_function"])
				if var["five_prime_gene_function"] not in gene_function:
					gene_function.append(var["five_prime_gene_function"])
		else:
			for var in info:
				if var["gene_function"] not in gene_function:
					gene_function.append(var["gene_function"])
		
		###2. 变异和证据格式转化
		# 转为格式1
		# [{药物1：变异1}, {辅助诊断1：变异2}, {预后1：变异1}]
		tmp_list1 = []
		for var in info:
			var_str = ""
			if var["bio_category"] == "Snvindel":
				var_str = var["gene_symbol"]+" "+var["hgvs_p_ZJZL"] if var["hgvs_p_ZJZL"] != "p.?" else var["gene_symbol"] + " "+ var["hgvs_c"]
				if "judge_mergeMET" in var.keys() and var["judge_mergeMET"]:
					var_str += "（MET exon14 skipping）"
			elif var["bio_category"] == "Cnv":
				var_str = var["gene_symbol"]+" 扩增"
			else:
				if var["five_prime_gene"] == "MET" and var["three_prime_gene"] == "MET":
					var_str = "MET exon14 skipping"
				else:
					if "var_desc_merge" in var.keys() and var["var_desc_merge"]:
						var_str = var["var_desc_merge"]+"融合"
					else:
						var_str = "{0}:{1}-{2}:{3}融合".format(var["five_prime_gene"], var["five_prime_cds"], var["three_prime_gene"], var["three_prime_cds"])
			if "Predictive" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Predictive"]:
				for evi in var["evi_sum"]["evi_split"]["Predictive"]:
					if evi["evi_conclusion_simple"] in ["A", "B"]:
						regimen = "{0}（{1}，{2}级）".format(evi["regimen_name"], evi["clinical_significance_cn"], evi["evi_conclusion_simple"])
						tmp_list1.append({("Predictive", regimen) : var_str})
			if "Prognostic" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Prognostic"]:
				for evi in var["evi_sum"]["evi_split"]["Prognostic"]:
					if evi["evi_conclusion_simple"] in ["A", "B"]:
						tmp_list1.append({("Prognostic", evi["evi_interpretation"]) : var_str})
			if "Diagnostic" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Diagnostic"]:
				for evi in var["evi_sum"]["evi_split"]["Diagnostic"]:
					if evi["evi_conclusion_simple"] in ["A", "B"]:
						tmp_list1.append({("Diagnostic", evi["evi_interpretation"]) : var_str})
		
		# 转为格式2
		# {药物1：变异1、变异2，药物2：变异3}
		tmp_dict2 = {}
		for i in tmp_list1:
			for k, v in i.items():
				if k not in tmp_dict2.keys():
					tmp_dict2.setdefault(k, [])
				tmp_dict2[k].append(v)

		# 转为格式3
		# {变异1、变异2：药物1、药物2、药物3，变异3：药物3}
		tmp_dict3 = {}
		for evi, var in tmp_dict2.items():
			if tuple(var) not in tmp_dict3.keys():
				tmp_dict3.setdefault(tuple(var), {})
			if evi[0] not in tmp_dict3[tuple(var)].keys():
				tmp_dict3[tuple(var)].setdefault(evi[0], [])
			tmp_dict3[tuple(var)][evi[0]].append(evi[1])
		
		# 转为格式4
		# [{var:变异1，变异2，type:治疗/辅助诊断/预后，inter:描述信息}]
		tmp_list4 = []
		for var, info in tmp_dict3.items():
			for type_, inter_ in info.items():
				#print (evi, var)
				tmp_list4.append({
					"evi_name" : "、".join(list(var)),
					"evi_type" : type_,
					"inter" : "、".join(inter_) if type_ == "Predictive" else "".join(inter_)
				})
	
		result.append({
			"gene_symbol" : gene,
			"gene_function" : gene_function,
			"var_info" : tmp_list4
		})
		
	return result	
jinja2.filters.FILTERS["ynzl_var_inter"] = ynzl_var_inter

# 湖南省肿瘤BRCA-输出指定变异-2025.03.06
def hnzl_var_sum(var_list):
	result = []
	for var in var_list:
		if var["bio_category"] == "Snvindel":
			if var["hgvs_p"] != "p.?":
				result.append("{0} {1}".format(var["gene_symbol"], var["hgvs_p"]))
			else:
				result.append("{0} {1}".format(var["gene_symbol"], var["hgvs_c"]))
		elif var["type"] == "Loss":
			result.append("{0} {1} del".format(var["gene_symbol"], var["value"]))
		elif var["type"] == "Gain":
			result.append("{0} {1} dup".format(var["gene_symbol"], var["value"]))
	return "、".join(result)
jinja2.filters.FILTERS["hnzl_var_sum"] = hnzl_var_sum

# 云肿结果小结-相同基因的要汇总到一起展示-2025.03.07
def ynzl_brca_summary_var(var_list):
	'''
	BRCA1基因：p.G1202R、exon16 del、exon16 dup
	'''
	# 获取基因列表
	gene_list = []
	for var in var_list:
		if var["gene_symbol"] not in gene_list:
			gene_list.append(var["gene_symbol"])
	
	result = {}
	for var in var_list:
		var_info = ""
		if var["bio_category"] == "Snvindel":
			var_info = var["hgvs_p_ZJZL"] if var["hgvs_p_ZJZL"] != "p.?" else var["hgvs_c"]
		elif var["type"] == "Loss":
			var_info = var["value"] + " del"
		elif var["type"] == "Gain":
			var_info = var["value"] + " dup"
		
		if var["gene_symbol"] not in result.keys():
			result.setdefault(var["gene_symbol"], [])
		result[var["gene_symbol"]].append(var_info)

	result_list = []
	for gene, info in result.items():
		result_list.append(
			{
				"gene_symbol" : gene,
				"var_info" : "、".join(info)
			}
		)
	result_list = sorted(result_list, key=lambda i:gene_list.index(i["gene_symbol"]))
	return result_list
jinja2.filters.FILTERS["ynzl_brca_summary_var"] = ynzl_brca_summary_var

# 云南肿瘤详细解读部分-BRCA-每个基因一行，若有用药相同的变异需要合并展示-2025.03.07
def ynzl_brca_var_inter(var_list):
	gene_dict = {}
	for var in var_list:
		if var["gene_symbol"] not in gene_dict.keys():
			gene_dict.setdefault(var["gene_symbol"], [])
		gene_dict[var["gene_symbol"]].append(var)

	# 格式转化
	# [{gene_symbol:XX,"gene_function":XX,var_info:[1,2,3]}]
	result = []
	for gene, info in gene_dict.items():
		###1. 获取基因介绍、遗传风险
		gene_function = []
		risk_inter = []
		for var in info:
			if var["gene_function"] not in gene_function:
				gene_function.append(var["gene_function"])
			if "evi_split" in var["evi_sum"].keys() and var["evi_sum"]["evi_split"] and \
				"Predisposing" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Predisposing"]:
				for evi in var["evi_sum"]["evi_split"]["Predisposing"]:
					if evi["evi_interpretation"] not in risk_inter:
						risk_inter.append(evi["evi_interpretation"])
		
		###2. 变异和证据格式转化
		# 转为格式1
		# [{药物1：变异1}, {辅助诊断1：变异2}, {预后1：变异1}]
		tmp_list1 = []
		for var in info:
			var_str = ""
			if var["bio_category"] == "Snvindel":
				var_str = var["gene_symbol"]+" "+var["hgvs_p_ZJZL"] if var["hgvs_p_ZJZL"] != "p.?" else var["gene_symbol"] + " "+ var["hgvs_c"]
			elif var["type"] == "Loss":
				var_str = var["gene_symbol"]+" "+var["value"] + " del"
			elif var["type"] == "Gain":
				var_str = var["gene_symbol"]+" "+var["value"] + " dup"

			if "Predictive" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Predictive"]:
				for evi in var["evi_sum"]["evi_split"]["Predictive"]:
					if evi["evi_conclusion_simple"] in ["A", "B"]:
						regimen = "{0}（{1}，{2}级）".format(evi["regimen_name"], evi["clinical_significance_cn"], evi["evi_conclusion_simple"])
						tmp_list1.append({("Predictive", regimen) : var_str})
		
		# 转为格式2
		# {药物1：变异1、变异2，药物2：变异3}
		tmp_dict2 = {}
		for i in tmp_list1:
			for k, v in i.items():
				if k not in tmp_dict2.keys():
					tmp_dict2.setdefault(k, [])
				tmp_dict2[k].append(v)

		# 转为格式3
		# {变异1、变异2：药物1、药物2、药物3，变异3：药物3}
		tmp_dict3 = {}
		for evi, var in tmp_dict2.items():
			if tuple(var) not in tmp_dict3.keys():
				tmp_dict3.setdefault(tuple(var), {})
			if evi[0] not in tmp_dict3[tuple(var)].keys():
				tmp_dict3[tuple(var)].setdefault(evi[0], [])
			tmp_dict3[tuple(var)][evi[0]].append(evi[1])
		
		# 转为格式4
		# [{var:变异1，变异2，type:治疗/辅助诊断/预后，inter:描述信息}]
		tmp_list4 = []
		for var, info in tmp_dict3.items():
			for type_, inter_ in info.items():
				#print (evi, var)
				tmp_list4.append({
					"evi_name" : "、".join(list(var)),
					"evi_type" : type_,
					"inter" : "、".join(inter_) if type_ == "Predictive" else "".join(inter_)
				})
	
		result.append({
			"gene_symbol" : gene,
			"gene_function" : gene_function,
			"risk_inter" : risk_inter,
			"var_info" : tmp_list4
		})
		
	return result	
jinja2.filters.FILTERS["ynzl_brca_var_inter"] = ynzl_brca_var_inter

# 云南肿瘤详细解读部分-每个基因一行，若有用药相同(所有用药)的变异需要合并展示-2025.03.10
def ynzl_var_inter_v2(var_list):
	#gene_list = []
	#for var in var_list:
	#	if var["gene_symbol"] not in gene_list:
	#		gene_list.append(gene_list)
	
	gene_dict = {}
	for var in var_list:
		if var["gene_symbol"] not in gene_dict.keys():
			gene_dict.setdefault(var["gene_symbol"], [])
		gene_dict[var["gene_symbol"]].append(var)

	# 格式转化
	# [{gene_symbol:XX,"gene_function":XX,var_info:[1,2,3]}]
	result = []
	for gene, info in gene_dict.items():
		###1. 获取基因介绍
		gene_function = []
		if "," in gene:
			for var in info:
				if var["three_prime_gene_function"] not in gene_function:
					gene_function.append(var["three_prime_gene_function"])
				if var["five_prime_gene_function"] not in gene_function:
					gene_function.append(var["five_prime_gene_function"])
		else:
			for var in info:
				if var["gene_function"] not in gene_function:
					gene_function.append(var["gene_function"])
		risk_inter = []
		if "evi_split" in var["evi_sum"].keys() and var["evi_sum"]["evi_split"] and \
				"Predisposing" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Predisposing"]:
				for evi in var["evi_sum"]["evi_split"]["Predisposing"]:
					if evi["evi_interpretation"] not in risk_inter:
						risk_inter.append(evi["evi_interpretation"])
		
		###2. 变异和证据格式转化
		# 转为格式1
		# [{(用药，变异1) : 证据}, {(诊断，变异2) : 证据}, {(预后，变异2) : 证据}]
		tmp_list1 = []
		for var in info:
			var_str = ""
			if var["bio_category"] == "Snvindel":
				var_str = var["gene_symbol"]+" "+var["hgvs_p_ZJZL"] if var["hgvs_p_ZJZL"] != "p.?" else var["gene_symbol"] + " "+ var["hgvs_c"]
				if "judge_mergeMET" in var.keys() and var["judge_mergeMET"]:
					var_str += "（MET exon14 skipping）"
			elif var["bio_category"] == "Cnv":
				if var["var_origin"] == "germline":
					if var["type"] == "Loss":
						var_str = var["gene_symbol"] + " " + var["value"] + " del"
					elif var["type"] == "Gain":
						var_str = var["gene_symbol"] + " " + var["value"] + " dup"
				else:
					var_str = var["gene_symbol"]+" 扩增"
			else:
				if var["five_prime_gene"] == "MET" and var["three_prime_gene"] == "MET":
					var_str = "MET exon14 skipping"
				else:
					if "var_desc_merge" in var.keys() and var["var_desc_merge"]:
						var_str = var["var_desc_merge"]+"融合"
					else:
						var_str = "{0}:{1}-{2}:{3}融合".format(var["five_prime_gene"], var["five_prime_cds"], var["three_prime_gene"], var["three_prime_cds"])

			if "Predictive" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Predictive"]:
				tmp_regimen_sum = []
				for evi in var["evi_sum"]["evi_split"]["Predictive"]:
					if evi["evi_conclusion_simple"] in ["A", "B"]:
						regimen = "{0}（{1}，{2}级）".format(evi["regimen_name"], evi["clinical_significance_cn"], evi["evi_conclusion_simple"])
						if regimen not in tmp_regimen_sum:
							tmp_regimen_sum.append(regimen)
				if tmp_regimen_sum:
					tmp_list1.append({("Predictive", var_str) : "、".join(tmp_regimen_sum) + "。"})
						
			if "Prognostic" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Prognostic"]:
				tmp_prognostiv_sum = []
				for evi in var["evi_sum"]["evi_split"]["Prognostic"]:
					if evi["evi_conclusion_simple"] in ["A", "B"]:
						if evi["evi_interpretation"] not in tmp_prognostiv_sum:
							tmp_prognostiv_sum.append(evi["evi_interpretation"])
				if 	tmp_prognostiv_sum:
					tmp_list1.append({("Prognostic", var_str) : "".join(tmp_prognostiv_sum)})
				
			if "Diagnostic" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Diagnostic"]:
				tmp_diagnostic_sum = []
				for evi in var["evi_sum"]["evi_split"]["Diagnostic"]:
					if evi["evi_conclusion_simple"] in ["A", "B"]:
						if evi["evi_interpretation"] not in tmp_diagnostic_sum:
							tmp_diagnostic_sum.append(evi["evi_interpretation"])
				if tmp_diagnostic_sum:
					tmp_list1.append({("Diagnostic", var_str) : "".join(tmp_diagnostic_sum)})
		
		# 转为格式2
		# {证据1 ： [（变异1，类型1）, （变异2，类型1）]}
		tmp_dict2 = {}
		for i in tmp_list1:
			for var_info, evi in i.items():
				if evi not in tmp_dict2.keys():
					tmp_dict2.setdefault(evi, [])
				tmp_dict2[evi].append(var_info)
		
		# 转为格式3
		tmp_list3 = []
		for evi, var_info in tmp_dict2.items():
			tmp_list3.append(
				{
					"evi_name" : "、".join([var[1] for var in var_info]),
					"evi_type" : [var[0] for var in var_info][0],
					"inter" : evi
				}
			)
	
		result.append({
			"gene_symbol" : gene,
			"gene_function" : gene_function,
			"var_info" : tmp_list3,
			"risk_inter" : risk_inter
		})
		
	return result	
jinja2.filters.FILTERS["ynzl_var_inter_v2"] = ynzl_var_inter_v2

# 云南肿瘤CP40-融合检出多个变异型时，不同变异型之间用";"隔开-2025.03.10
def ynzl_sv_stran(var_list):
	if len(var_list) > 1:
		for var in var_list[0:-1]:
			var["three_prime_cds"] = var["three_prime_cds"] + ";"
	return var_list
jinja2.filters.FILTERS["ynzl_sv_stran"] = ynzl_sv_stran

# 云南肿瘤CP40-融合检出多个变异型时，拷贝数之间用";"隔开-2025.03.10
def ynzl_sv_copies_stran(copies_list):
	for var in copies_list:
		var["copies"] += " copies"
	if len(copies_list) > 1:
		for var in copies_list[0:-1]:
			var["copies"] += ";"
	return copies_list
jinja2.filters.FILTERS["ynzl_sv_copies_stran"] = ynzl_sv_copies_stran

# 云南肿瘤详细解读部分BRCA-每个基因一行，若有用药相同(所有用药)的变异需要合并展示-2025.03.10
def ynzl_brca_var_inter_v2(var_list):
	gene_dict = {}
	for var in var_list:
		if var["gene_symbol"] not in gene_dict.keys():
			gene_dict.setdefault(var["gene_symbol"], [])
		gene_dict[var["gene_symbol"]].append(var)
		
	# 格式转化
	# [{gene_symbol:XX,"gene_function":XX,var_info:[1,2,3]}]
	result = []
	risk_inter = []
	for gene, info in gene_dict.items():
		###1. 获取基因介绍
		gene_function = []
		if var["gene_function"] not in gene_function:
			gene_function.append(var["gene_function"])
		if "evi_split" in var["evi_sum"].keys() and var["evi_sum"]["evi_split"] and \
				"Predisposing" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Predisposing"]:
				for evi in var["evi_sum"]["evi_split"]["Predisposing"]:
					if evi["evi_interpretation"] not in risk_inter:
						risk_inter.append(evi["evi_interpretation"])
		
		###2. 变异和证据格式转化
		# 转为格式1
		# [{(用药，变异1) : 证据}, {(诊断，变异2) : 证据}, {(预后，变异2) : 证据}]
		tmp_list1 = []
		for var in info:
			var_str = ""
			if var["bio_category"] == "Snvindel":
				var_str = var["gene_symbol"]+" "+var["hgvs_p_ZJZL"] if var["hgvs_p_ZJZL"] != "p.?" else var["gene_symbol"] + " "+ var["hgvs_c"]
			elif var["type"] == "Loss":
				var_str = var["gene_symbol"]+" "+var["value"] + " del"
			elif var["type"] == "Gain":
				var_str = var["gene_symbol"]+" "+var["value"] + " dup"

			if "Predictive" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Predictive"]:
				tmp_regimen_sum = []
				for evi in var["evi_sum"]["evi_split"]["Predictive"]:
					if evi["evi_conclusion_simple"] in ["A", "B"]:
						regimen = "{0}（{1}，{2}级）".format(evi["regimen_name"], evi["clinical_significance_cn"], evi["evi_conclusion_simple"])
						if regimen not in tmp_regimen_sum:
							tmp_regimen_sum.append(regimen)
				if tmp_regimen_sum:
					tmp_list1.append({("Predictive", var_str) : "、".join(tmp_regimen_sum) + "。"})
		
		# 转为格式2
		# {证据1 ： [（变异1，类型1）, （变异2，类型1）]}
		tmp_dict2 = {}
		for i in tmp_list1:
			for var_info, evi in i.items():
				if evi not in tmp_dict2.keys():
					tmp_dict2.setdefault(evi, [])
				tmp_dict2[evi].append(var_info)
		
		# 转为格式3
		tmp_list3 = []
		for evi, var_info in tmp_dict2.items():
			tmp_list3.append(
				{
					"evi_name" : "、".join([var[1] for var in var_info]),
					"evi_type" : [var[0] for var in var_info][0],
					"inter" : evi
				}
			)
	
		result.append({
			"gene_symbol" : gene,
			"gene_function" : gene_function,
			"var_info" : tmp_list3,
			"risk_inter" : risk_inter
		})
		
	return result	
jinja2.filters.FILTERS["ynzl_brca_var_inter_v2"] = ynzl_brca_var_inter_v2

# 云南肿瘤HRD-汇总BRCA变异对应参考药物-2025.03.10
def ynzl_hrd_brca_regimen(var_list):
	result = []
	for var in var_list:
		tmp_regimen_sum = []
		if "Predictive" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Predictive"]:
				if var["clinic_num_s"] == 5:
					tmp_regimen_sum = ["{0}（{1}，{2}级）".format(evi["regimen_name"], evi["clinical_significance_cn"], evi["evi_conclusion_simple"]) for \
										evi in var["evi_sum"]["evi_split"]["Predictive"] if evi["evi_conclusion_simple"] in ["A", "B"]]
					
				else:
					tmp_regimen_sum = ["{0}（{1}，{2}级）".format(evi["regimen_name"], evi["clinical_significance_cn"], evi["evi_conclusion_simple"]) for \
										evi in var["evi_sum"]["evi_split"]["Predictive"]]
		for evi in tmp_regimen_sum:
				if evi not in result:
					result.append(evi)
	return result
jinja2.filters.FILTERS["ynzl_hrd_brca_regimen"] = ynzl_hrd_brca_regimen

# 湖南省肿瘤ptBRCA-需要计算体细胞变异个数-2025.03.11
def hnzl_ptbrca_somatic_num(var_list):
	return [var for var in var_list if var["var_origin"] == "somatic"]
jinja2.filters.FILTERS["hnzl_ptbrca_somatic_num"] = hnzl_ptbrca_somatic_num

# 新疆附一HRD-需要过滤出BRCA和其他-2025.03.11
def xjfy_hrd_brca(var_list):
	return [var for var in var_list if var["gene_symbol"] in ["BRCA1", "BRCA2"]]
jinja2.filters.FILTERS["xjfy_hrd_brca"] = xjfy_hrd_brca

def xjfy_hrd_other(var_list):
	return [var for var in var_list if var["gene_symbol"] not in ["BRCA1", "BRCA2"]]
jinja2.filters.FILTERS["xjfy_hrd_other"] = xjfy_hrd_other

# 福建省立150-需要区分是10基因还是其他基因-2025.03.12
def fjsl_150_lc10(var_list):
	return [var for var in var_list if var["gene_symbol"] in ["ALK", "BRAF", "EGFR", "MET", "RET"]]
jinja2.filters.FILTERS["fjsl_150_lc10"] = fjsl_150_lc10

def fjsl_150_other(var_list):
	return [var for var in var_list if var["gene_symbol"] not in ["ALK", "BRAF", "EGFR", "MET", "RET"]]
jinja2.filters.FILTERS["fjsl_150_other"] = fjsl_150_other

# 新疆附一CP200-过滤出10个基因中未检出I/II变异的基因-2025.03.14
def xjfy_cp200_withoutPathVar_lc10(var_list):
	lc10_gene = ["ALK", "BRAF", "EGFR", "ERBB2", "KRAS", "MET", "NRAS", "PIK3CA", "RET", "ROS1"]
	detect_gene = []
	for var in var_list:
		gene = re.split(",", var["gene_symbol"])
		detect_gene.extend(gene)
	return sorted(list(set(lc10_gene) - set(detect_gene)))
jinja2.filters.FILTERS["xjfy_cp200_withoutPathVar_lc10"] = xjfy_cp200_withoutPathVar_lc10

#---------------------中山六院CP200相关 开始----------------------------------------------------#
# 中山六院CP200证据合并-2025.03.17
# 输入证据列表和类型（同癌种/其他癌种）
def zsly_evi_merge(evi_list, tumor_type):
	evi_stran = {}
	for evi in evi_list:
		evi["tumor_name_cn"] = evi["tumor_name_cn"] if evi["tumor_name_cn"] else evi["tumor_name_en"] if "tumor_name_en" in evi.keys() and evi["tumor_name_en"] else ""
		key_intumor = (evi["evidence_type"], evi["regimen_name"], evi["clinical_significance_cn"], evi["evi_conclusion_simple"])
		key_outtumor = (evi["evidence_type"], evi["regimen_name"], evi["clinical_significance_cn"], evi["evi_conclusion_simple"], evi["tumor_name_cn"])
		if tumor_type == "intumor":
			key_ = key_intumor
		else:
			key_ = key_outtumor
		if key_ not in evi_stran.keys():
			evi_stran.setdefault(key_, [])
		evi_stran[key_].append(evi)

	evi_result = []
	for k, v in evi_stran.items():
		value = v[0]
		# 2025.05.30-汇总下证据对应肿瘤列表
		tumor_name_cn_list = []
		for evi in v:
			tumor_name_cn = evi["tumor_name_cn"] if evi["tumor_name_cn"] else evi["tumor_name_en"] if "tumor_name_en" in evi.keys() and evi["tumor_name_en"] else ""
			if tumor_name_cn not in tumor_name_cn_list:
				tumor_name_cn_list.append(tumor_name_cn)
		value["tumor_name_list"] = "、".join(tumor_name_cn_list)
		#print (k, value["tumor_name_list"])
		# 2025.05.30-新增完成
		evi_origin = []
		# 2025.05.16-新增参考文献evi_refer_dict
		evi_refer_dict = {}
		# 2025.05.16-新增完成
		for evi in v:
			if evi["refer_agency"] and evi["refer_agency"]:
				#if evi["refer_agency"] not in evi_origin:
				#	evi_origin.append(evi["refer_agency"])
				# 2025.04.21-引用机构可能会多条，这边兼容
				for i in re.split(",", evi["refer_agency"]):
					if i not in evi_origin:
						evi_origin.append(i)
					# 2025.05.16-新增参考文献
					if i not in evi_refer_dict.keys():
						evi_refer_dict.setdefault(i, [])
					#print (re.split(",", evi["literature_evi_sum"]))
					#evi_refer_dict[i].extend(re.split(",", evi["literature_evi_sum"]))
					# 2025.06.03-加个兼容，原来没有参考文献是个列表，里面是空值，现在改为null
					evi["literature_evi_sum"] = evi["literature_evi_sum"] if evi["literature_evi_sum"] else []
					for refer_item in evi["literature_evi_sum"]:
						if refer_item:
							evi_refer_dict[i].append(refer_item)
					# 2025.05.16-新增完成
			elif evi["evidence_level"] and zsly_evidence_level_stran(evi["evidence_level"]):
				if evi["evidence_level"] not in evi_origin:
					evi_origin.append(zsly_evidence_level_stran(evi["evidence_level"]))
					# 2025.05.16-新增参考文献
					#evi_refer_dict[zsly_evidence_level_stran(evi["evidence_level"])] = []
					# 2025.06.13-debug
					if zsly_evidence_level_stran(evi["evidence_level"]) not in evi_refer_dict.keys():
						evi_refer_dict.setdefault(zsly_evidence_level_stran(evi["evidence_level"]), [])
					# 2025.06.13-debug完成
					#evi_refer_dict[zsly_evidence_level_stran(evi["evidence_level"])].extend(re.split(",", evi["literature_evi_sum"]))
					# 2025.06.03-加个兼容，原来没有参考文献是个列表，里面是空值，现在改为null
					evi["literature_evi_sum"] = evi["literature_evi_sum"] if evi["literature_evi_sum"] else []
					for refer_item in evi["literature_evi_sum"]:
						if refer_item:
							evi_refer_dict[zsly_evidence_level_stran(evi["evidence_level"])].append(refer_item)
					#print ("111111111111111", evi["regimen_name"], evi_refer_dict[zsly_evidence_level_stran(evi["evidence_level"])])
					# 2025.05.16-新增完成
		#print ("111111111111111111111",evi_refer_dict)
		evi_origin_set = []
		for origin in evi_origin:
			if origin not in evi_origin_set:
				evi_origin_set.append(origin)
		value["evi_origin"] = "、".join(evi_origin_set)
		# 2025.05.16-新增参考文献
		value["evi_refer_dict"] = evi_refer_dict
		# 2025.05.16-新增完成
		evi_result.append(value)
	#for evi in evi_result:
	#	print ("11111111111", evi["regimen_name"], evi["evi_refer_dict"])
	return evi_result

# 1. 同等级变异需要根据snvindel>cnv>sv进行排序 - 2025.03.14
def zsly_cp200_sort_var(var_list):
	sort_rule = ["Snvindel", "Cnv", "Sv"]
	return sorted(var_list, key = lambda i : sort_rule.index(i["bio_category"]))
jinja2.filters.FILTERS["zsly_cp200_sort_var"] = zsly_cp200_sort_var

# 2. 治疗方案没有获批机构的，将研究类型进行转化-2025.03.14
def zsly_evidence_level_stran(evi_level):
	evi_level_dict = {
		"Clinical-phase I" : "临床研究",
		"Clinical-phase II" : "临床研究",
		"Clinical-phase III" : "临床研究",
		"Clinical-phase IV" : "临床研究",
		"Clinical-retrospective" : "临床研究",
		"Clinical-unknown phase" : "临床研究",
		"Case report" : "案例报道",
		"Preclinical-in vitro" : "临床前研究",
		"Preclinical-in vivo" : "临床前研究"
		}
	return evi_level_dict.get(evi_level, evi_level) if evi_level else ""
jinja2.filters.FILTERS["zsly_evidence_level_stran"] = zsly_evidence_level_stran

# 3. 单突变-获取同癌种证据 - 2025.03.14
def zsly_cp200_filter_intumor_evi(evi_sum):
	evi_all = []
	if evi_sum and "evi_split" in evi_sum.keys() and evi_sum["evi_split"]:
		for evi_type in ["Predictive", "Prognostic", "Diagnostic"]:
			if evi_type in evi_sum["evi_split"] and evi_sum["evi_split"][evi_type]:
				evi_all.extend(evi_sum["evi_split"][evi_type])
	intumor_evi = [evi for evi in evi_all if evi["evi_conclusion_simple"] in ["A", "B"] or evi["evi_conclusion"] in ["C1", "C2", "D1", "D2", "D5", "E1", "E2"]]
	#for evi in intumor_evi:
	#	print ("1111111111", evi["regimen_name"])
	merge_intumor_evi = zsly_evi_merge(intumor_evi, "intumor")
	#for i in merge_intumor_evi:
		#print ("                    1111111111", i)
	return sorted(merge_intumor_evi, key = lambda i : (i["evi_conclusion_simple"], i["sense_rule"], i["regimen_name_py"].upper()))
jinja2.filters.FILTERS["zsly_cp200_filter_intumor_evi"] = zsly_cp200_filter_intumor_evi

# 4. 单突变-获取其他癌种证据 - 2025.03.14
def zsly_cp200_filter_outtumor_evi(info):
	var_level = info[0]
	evi_sum = info[1]
	evi_all = []
	if evi_sum and "evi_split" in evi_sum.keys() and evi_sum["evi_split"]:
		for evi_type in ["Predictive", "Prognostic", "Diagnostic"]:
			if evi_type in evi_sum["evi_split"] and evi_sum["evi_split"][evi_type]:
				evi_all.extend(evi_sum["evi_split"][evi_type])
	C_evi = [evi for evi in evi_all if evi["evi_conclusion"] in ["C3", "C4"]]
	merege_C_evi = zsly_evi_merge(C_evi, "outtumor")
	merege_C_evi = sorted(merege_C_evi, key = lambda i : (i["evi_conclusion_simple"], i["sense_rule"], i["regimen_name_py"].upper()))
	D_evi = [evi for evi in evi_all if evi["evi_conclusion"] in ["D3", "D4"]]
	merege_D_evi = zsly_evi_merge(D_evi, "outtumor")
	merege_D_evi = sorted(merege_D_evi, key = lambda i : (i["evi_conclusion_simple"], i["sense_rule"], i["regimen_name_py"].upper()))
	if var_level == "I":
		return merege_C_evi
	else:
		if merege_C_evi:
			return merege_C_evi
		else:
			return merege_D_evi
jinja2.filters.FILTERS["zsly_cp200_filter_outtumor_evi"] = zsly_cp200_filter_outtumor_evi

# *** 融合位点合并 ***
def zsly_sv_merge(var_list_raw):
	var_list = copy.deepcopy(var_list_raw)
	#print (var_list)
	# 这段待定
	#for var in var_list:
	#	if "bio_category" not in var.keys() or ("bio_category" in var.keys and not var["bio_category"]):
	#		var["bio_category"]	 = var["var_category"]
	sv_merge = []
	sv_data = [var for var in var_list if var["bio_category"] == "Sv"]
	snvindel_data = [var for var in var_list if var["bio_category"] == "Snvindel"]
	cnv_data = [var for var in var_list if var["bio_category"] == "Cnv"]
	special_markers = [var for var in var_list if var["bio_category"] == "Special_markers"]
	sv_dict = {}
	for var in sv_data:
		var["five_prime_gene"] = re.split(":", var["var_hgvs"])[0]
		var["three_prime_gene"] = re.split(":", (re.split("-", var["var_hgvs"])[-1]))[0]
		# var_hgvs兼容4种格式格式-参考sv处理代码
		var["five_prime_cds"] = "-".join(re.split("-", (re.split(":", var["var_hgvs"])[2]))[:-1]) \
								if not re.search("--", var["var_hgvs"]) else re.split("_", (re.split("--", var["var_hgvs"])[0]))[-1]
		var["three_prime_cds"] = re.split(":", var["var_hgvs"])[-1] if not re.search("--", var["var_hgvs"]) else re.split("_", (re.split("--", var["var_hgvs"])[1]))[-1]
		var["five_prime_cds"] = re.split(":", var["five_prime_cds"])[-1] if re.search(":", var["five_prime_cds"]) else var["five_prime_cds"]
		var["three_prime_cds"] = re.split(":", var["three_prime_cds"])[-1] if re.search(":", var["three_prime_cds"]) else var["three_prime_cds"]
		if not var["five_prime_cds"]:
			var["five_prime_cds"] = re.split("_", re.split("-", var["var_hgvs"])[0])[-1]
			var["three_prime_cds"] = re.split("_", var["var_hgvs"])[-1]

	for var in sv_data:
		key = (var["five_prime_gene"], var["three_prime_gene"])
		if key not in sv_dict.keys():
			sv_dict.setdefault(key, [])
		sv_dict[key].append(var)
	
	for key, value in sv_dict.items():
		sv_merge.append(
			{
				"gene_symbol" : value[0]["gene_symbol"],
				"bio_category" : "Sv",
				"five_prime_gene" : key[0],
				"three_prime_gene" : key[1],
				"sv_list" : [{"five_prime_cds" : i["five_prime_cds"], "three_prime_cds" : i["three_prime_cds"], "copies" : i["copies"]} for i in value]
			}
		)
	#print (snvindel_data + cnv_data + sv_merge + special_markers)
	return snvindel_data + cnv_data + sv_merge + special_markers
jinja2.filters.FILTERS["zsly_sv_merge"] = zsly_sv_merge

# 5. 判断是否有共突变-2025.03.14
def zsly_co_mutation(info):
	var_list = info[0]
	knb = info[1]
	msi = info[2]
	var_level = info[3]
	if knb:
		var_list.append(knb)
	if msi["var_id"] == "MSI-H":
		var_list.append(msi)
	### 汇总所有共突变证据 ###
	co_evi_sum = []
	for var in var_list:
		evi_sum = []
		#print (type(var))
		for evi_type in ["Predictive", "Prognostic", "Diagnostic"]:
			if "evi_split" in var["evi_sum"].keys() and var["evi_sum"]["evi_split"] and \
			evi_type in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"][evi_type]:
				evi_sum.extend(var["evi_sum"]["evi_split"][evi_type])
		co_evi_sum.extend([evi for evi in evi_sum if "co_mutation" in evi.keys() and evi["co_mutation"]])

	# 2025.04.11-新增-对co_mutation中的变异进行排序，
	# 1. 防止json返回中同一个共突变在不同证据里顺序不一致
	# 2. 防止不同样本返回共突变顺序不一致
	for evi in co_evi_sum:
		var_list = evi["molecular_var"]
		# snvindel，基因升、频率降，freq兼容下小数、百分数和其他
		snvindel_data = [var for var in var_list if var["bio_category"] == "Snvindel"]
		for var in snvindel_data:
			var["freq"] = float(var["freq"].replace("%", ""))/100 if re.search("%", str(var["freq"])) else float(var["freq"]) if var["freq"] and is_number(var["freq"]) else var["freq"]
		snvindel_data = sorted(snvindel_data, key = lambda i :  i["freq"], reverse=True)
		snvindel_data = sorted(snvindel_data, key = lambda i : i["gene_symbol"])
		# cnv，基因升、拷贝数降
		cnv_data = [var for var in var_list if var["bio_category"] == "Cnv"]
		cnv_data = sorted(cnv_data, key = lambda i : float(i["cn_mean"]), reverse=True)
		cnv_data = sorted(cnv_data, key = lambda i : i["gene_symbol"])
		# sv，基因升，拷贝数降
		sv_data = [var for var in var_list if var["bio_category"] == "Sv"]
		sv_data = sorted(sv_data, key = lambda i : float(i["copies"]), reverse=True)
		sv_data = sorted(sv_data, key = lambda i : i["gene_symbol"])
		# 特殊标志物，字符串升
		special_markers = [var for var in var_list if var["bio_category"] == "Special_markers"]
		special_markers = sorted(special_markers, key = lambda i : i["biomarkers_name"])
		evi["molecular_var"] = snvindel_data + cnv_data + sv_data + special_markers
	# 2025.04.11-新增完成

	### 转化为变异:[证据] ###
	co_evi_dict = {}
	for co_evi in co_evi_sum:
		var_info = []
		for var in co_evi["molecular_var"]:
			tmp = tuple([tuple([i for i in var.keys()]), tuple([i for i in var.values()])])
			var_info.append(tmp)
		var_info_tuple = tuple(var_info)
		if var_info_tuple not in co_evi_dict.keys():
			co_evi_dict.setdefault(var_info_tuple, [])
		co_evi_dict[var_info_tuple].append(co_evi)

	### 判断共突变为I类还是II类 ###
	level_I_co_mutation = []
	level_II_co_mutation = []
	for var, evi in co_evi_dict.items():
		co_mutation = evi[0]["molecular_var"]
		evi_level_sum = [i["evi_conclusion_simple"] for i in evi]
		if set(["A", "B"]) & set(evi_level_sum):
			level_I_co_mutation.append(
				{
					"var" : co_mutation,
					"evi" : evi
				}
			)
		else:
			level_II_co_mutation.append(
				{
					"var" : co_mutation,
					"evi" : evi
				}
			)
	
	### 根据模板中调用参数，返回对应结果 ###
	if var_level == "I":
		return level_I_co_mutation
	else:
		return level_II_co_mutation
jinja2.filters.FILTERS["zsly_co_mutation"] = zsly_co_mutation

# 6. 共突变-获取同癌种证据 - 2025.03.14
def zsly_cp200_co_get_intumor_evi(evi_sum):
	intumor_evi = [evi for evi in evi_sum if evi["evi_conclusion_simple"] in ["A", "B"] or evi["evi_conclusion"] in ["C1", "C2", "D1", "D2", "D5", "E1", "E2"]]
	merge_intumor_evi = zsly_evi_merge(intumor_evi, "intumor")
	return sorted(merge_intumor_evi, key = lambda i : (i["evi_conclusion_simple"], i["sense_rule"], i["regimen_name_py"].upper()))
jinja2.filters.FILTERS["zsly_cp200_co_get_intumor_evi"] = zsly_cp200_co_get_intumor_evi

# 7. 共突变-获取其他癌种证据 - 2025.03.14
def zsly_co_get_outtumor_evi(info):
	var_level = info[0]
	evi_sum = info[1]
	C_evi = [evi for evi in evi_sum if evi["evi_conclusion"] in ["C3", "C4"]]
	merge_C_evi = zsly_evi_merge(C_evi, "outtumor")
	merge_C_evi = sorted(merge_C_evi, key = lambda i : (i["evi_conclusion_simple"], i["sense_rule"], i["regimen_name_py"].upper()))
	D_evi = [evi for evi in evi_sum if evi["evi_conclusion"] in ["D3", "D4"]]
	merge_D_evi = zsly_evi_merge(D_evi, "outtumor")
	merge_D_evi = sorted(merge_D_evi, key = lambda i : (i["evi_conclusion_simple"], i["sense_rule"], i["regimen_name_py"].upper()))
	if var_level == "I":
		return merge_C_evi
	else:
		if merge_C_evi:
			return merge_C_evi
		else:
			return merge_D_evi
jinja2.filters.FILTERS["zsly_co_get_outtumor_evi"] = zsly_co_get_outtumor_evi

# 8. 化疗按要求排序 - 2025.03.14
def zsly_chemo_sort(chemo_list):
	chemo_rule = ["1A", "1B", "2A", "2B", "3", "4"]
	for var in chemo_list:
		if var["gene_symbol"] == "UGT1A1":
			var["gene_index_num"] = 0
		else:
			var["gene_index_num"] = 1
	return sorted(chemo_list, key = lambda i : (i["gene_index_num"], i["gene_symbol"], chemo_rule.index(i["evi_level"]), i["dbsnp"]))
jinja2.filters.FILTERS["zsly_chemo_sort"] = zsly_chemo_sort

# 9. 化疗需要汇总野生型/变异型等 - 2025.03.17
def zsly_chemo_stran(chemo_list):
	var_or_wt = {
		"rs3918290" : "C",
		"rs55886062" : "A",
		"rs67376798" : "T",
		"rs75017182" : "G",
		"rs11615" : "A",
		"rs1695" : "A",
		"rs10929302" : "G",
		"rs4148323" : "G",
		"rs8175347" : "(TA)6"
	}
	result = {}
	for chemo in chemo_list:
		if chemo["gene_symbol"] not in result.keys():
			result.setdefault(chemo["gene_symbol"], [])
		genotype_str = ""
		if chemo["dbsnp"] in var_or_wt.keys():
			ref = var_or_wt[chemo["dbsnp"]]
			split_alt = re.split("/", chemo["genotype"])
			if split_alt[0] == ref and split_alt[1] == ref:
				genotype_str == "野生型"
			else:
				if len(set(split_alt)) == 1:
					genotype_str = "纯合变异"
				else:
					genotype_str = "杂合变异"
				result[chemo["gene_symbol"]].append("{0} {1}（{2}）".format(chemo["dbsnp"], chemo["genotype"], genotype_str))
	result_list = []
	for k, v in result.items():
		if not v:
			v = ["野生型"]
		result_list.append({
			"gene_symbol" : k,
			"result" : "；".join(v)
		})
	for i in result_list:
		if i["gene_symbol"] == "UGT1A1":
			i["gene_index_num"] = 0
		else:
			i["gene_index_num"] = 1
	return sorted(result_list, key = lambda i : (i["gene_index_num"], i["gene_symbol"]))
jinja2.filters.FILTERS["zsly_chemo_stran"] = zsly_chemo_stran

# 10. 免疫正负相关解析 - 2025.03.24
# 输入参数为[配置信息，检测结果，tumor_list，正相关/负相关]
def zsly_io_inter(info):
	database = info[0]
	detect_result = info[1]
	tumor_list = info[2]
	tumor_list.append("实体瘤")
	io_type = info[3]
	#io_gene_P = ["ATM","ATR","BRCA1","BRCA2","BRIP1","CHEK1","CHEK2","FANCA","MRE11",\
	#			 "PALB2","RAD50","MLH1","MSH2","MSH6","PMS2","POLE","POLD1","TP53",\
	#			 "KRAS","CD274","ARID1A","SETD2","TERT","KMT2D","FAT1","CDK12"]
	#2026.05.07 增加DNMT3A 孟智悦
	io_gene_P = ["ATM","ATR","BRCA1","BRCA2","BRIP1","CHEK1","CHEK2","FANCA","MRE11",\
				 "PALB2","RAD50","MLH1","MSH2","MSH6","PMS2","POLE","POLD1","TP53",\
				 "KRAS","CD274","ARID1A","SETD2","TERT","KMT2D","FAT1","CDK12","DNMT3A"]
	
	# 2025.04.07 删除基因CDKN2B
	#io_gene_N = ["EGFR","ALK","MDM2","MDM4","CDKN2A","DNMT3A","STK11","IFNGR1",\
	#			 "JAK1","JAK2","APC","CTNNB1","B2M","PTEN","FGF19"]
	# 2026.05.07 删除基因DNMT3A 孟智悦
	io_gene_N = ["EGFR","ALK","MDM2","MDM4","CDKN2A","STK11","IFNGR1",\
				 "JAK1","JAK2","APC","CTNNB1","B2M","PTEN","FGF19"]
	
	gene_list = io_gene_P if io_type == "p" else io_gene_N
	result = []
	for gene, info in detect_result.items():
		if gene in gene_list:
			intumor = ""
			outtumor = ""
			if gene in database.keys():
				if set(tumor_list) & set(re.split("、", database.get(gene)["tumor"])):
					intumor = (database.get(gene))["inter"]
				else:
					outtumor = (database.get(gene))["inter"]
			else:
				intumor = "配置表中未有相关信息！"
				outtumor = "配置表中未有相关信息！"
			result.append(
				{
					"gene_symbol" : gene,
					"var_info" : info,
					"intumor" :  intumor,
					"outtumor" : outtumor,
					# 2025.05.16-新增参考文献
					"refer" : (database.get(gene))["refer"]
					# 2025.05.16-新增结束
				}
			)
	return result
jinja2.filters.FILTERS["zsly_io_inter"] = zsly_io_inter

# 参考文献-PMID - 2025.03.25
def zsly_pmid(refer):
	pmid_list = []
	for i in refer:
		pmid = (re.split("PMID:", i)[-1]).replace("]", "")
		if pmid not in pmid_list:
			pmid_list.append(pmid)

	pmid_all_index = []
	num = 1
	for i in pmid_list:
		pmid_all_index.append({"num" : num, "pmid" : i})
		num += 1
	
	if pmid_all_index:
		rest_num = len(pmid_all_index) % 4
		key = int(len(pmid_all_index) / 4)
		list1 = []
		list2 = []
		list3 = []
		list4 = []
		if rest_num == 1:
			list1 = pmid_all_index[0 : key+1]
			list2 = pmid_all_index[key+1 : key*2+1]
			list2.append({"num" : "", "pmid" : ""})
			list3 = pmid_all_index[key*2+1 : key*3+1]
			list3.append({"num" : "", "pmid" : ""})
			list4 = pmid_all_index[key*3+1 :]
			list4.append({"num" : "", "pmid" : ""})
		elif rest_num == 2:
			list1 = pmid_all_index[0 : key+1]
			list2 = pmid_all_index[key+1 : key*2+2]
			list3 = pmid_all_index[key*2+2 : key*3+2]
			list3.append({"num" : "", "pmid" : ""})
			list4 = pmid_all_index[key*3+2 :]
			list4.append({"num" : "", "pmid" : ""})
		elif rest_num == 3:
			list1 = pmid_all_index[0 : key+1]
			list2 = pmid_all_index[key+1 : key*2+2]
			list3 = pmid_all_index[key*2+2 : key*3+3]
			list4 = pmid_all_index[key*3+3 :]
			list4.append({"num" : "", "pmid" : ""})
		else:
			list1 = pmid_all_index[0 : key]
			list2 = pmid_all_index[key : key*2]
			list3 = pmid_all_index[key*2 : key*3]
			list4 = pmid_all_index[key*3 :]

	result = []
	if pmid_all_index:
		for i in range(0, len(list1)):
			result.append({
				"num_1" : list1[i]["num"] if list1[i]["num"] else "",
				"pmid_1" : list1[i]["pmid"] if list1[i]["pmid"] else "",
				"num_2" : list2[i]["num"] if list2[i]["num"] else "",
				"pmid_2" : list2[i]["pmid"] if list2[i]["pmid"] else "",
				"num_3" : list3[i]["num"] if list3[i]["num"] else "",
				"pmid_3" : list3[i]["pmid"] if list3[i]["pmid"] else "",
				"num_4" : list4[i]["num"] if list4[i]["num"] else "",
				"pmid_4" : list4[i]["pmid"] if list4[i]["pmid"] else ""
			})

	return result
jinja2.filters.FILTERS["zsly_pmid"] = zsly_pmid

# freq改为百分数并且保留两位小数 - 2025.03.26
def zsly_freq(freq):
	freq = float(freq.replace("%", ""))/100 if re.search("%", str(freq)) else float(freq) if freq and is_number(freq) else freq
	return "{:.2%}".format(float(int(float(freq) * 10000 + 0.5) / 10000)) if freq else ""
jinja2.filters.FILTERS["zsly_freq"] = zsly_freq

# 免疫正负相关解析表中，融合要改为gene1-gene2融合
def zsly_io_sv(var_info):
	if "融合" in var_info:
		five_gene = re.split(":", var_info)[0]
		three_gene = re.split(":", (re.split("-", var_info)[-1]))[0]
		return "{0}-{1}融合".format(five_gene, three_gene)
	else:
		return var_info
jinja2.filters.FILTERS["zsly_io_sv"] = zsly_io_sv

# 2025.04.02
# 共突变，snvindel freq >= 0.01 且freq <= 0.05, 
# cnv ERBB2/MET cn_mean >= 3.5 且 cn_mean <= 5， cnv其他cn_mean >= 6且cn_mean <= 8添加上角标
def zsly_co_mutation_add_subscript(var_list):
	result = []
	for var in var_list:
		if "bio_category" in var.keys() and var["bio_category"] and var["bio_category"] == "Snvindel":
			if float(var["freq"]) >= 0.01 and float(var["freq"]) <= 0.05:
				if "Snvindel" not in result:
					result.append("Snvindel")
			elif float(var["freq"]) < 0.01:
				if "Snvindel_2" not in result: #修改Snvidel为Snvindel 2026.05.07 孟智悦
					result.append("Snvindel_2")
		elif "bio_category" in var.keys() and var["bio_category"] and var["bio_category"] == "Cnv":
			if var["gene_symbol"] in ["ERBB2", "MET"]:
				if float(var["cn_mean"]) >= 3.5 and float(var["cn_mean"]) <= 5:
					if "Cnv" not in result:
						result.append("Cnv")
			else:
				if float(var["cn_mean"]) >= 6 and float(var["cn_mean"]) <= 8:
					if "Cnv" not in result:
						result.append("Cnv")
	return result
jinja2.filters.FILTERS["zsly_co_mutation_add_subscript"] = zsly_co_mutation_add_subscript

# 化疗需要判断各位点频率-2025.04.02
def zsly_chemo_var(info):
	var = info[0]
	all_detect_result = info[1]
	freq = 0
	for i in all_detect_result:
		if i["dbsnp"] == var["dbsnp"]:
			freq = i["freq"]
	if (float(freq) >= 0.15 and float(freq) <= 0.25) or (float(freq) >= 0.75 and float(freq) <= 0.9):
		return True
	else:
		return False
jinja2.filters.FILTERS["zsly_chemo_var"] = zsly_chemo_var

# 化疗需要判断各位点频率-2025.04.14
# 需要展示频率区间
def zsly_chemo_var_v2(info):
	var = info[0]
	all_detect_result = info[1]
	freq = 0
	for i in all_detect_result:
		if i["dbsnp"] == var["dbsnp"]:
			freq = i["freq"]
	if (float(freq) >= 0.15 and float(freq) <= 0.25) or (float(freq) >= 0.75 and float(freq) <= 0.9):
		if (float(freq) >= 0.15 and float(freq) <= 0.25):
			return "loss"
		else:
			return "more"
	else:
		return False
jinja2.filters.FILTERS["zsly_chemo_var_v2"] = zsly_chemo_var_v2

# 化疗需要判断各位点频率-2025.04.27
# 需要展示频率区间
# 野生型/杂合型灰区阈值为>=15% <=25% 和 >=75% <85%
# 杂合型/纯合型灰区阈值为>=85% <=90%
def zsly_chemo_var_v3(info):
	var = info[0]
	all_detect_result = info[1]
	freq = 0
	for i in all_detect_result:
		if i["dbsnp"] == var["dbsnp"]:
			freq = i["freq"]
		# 返回阈值兼容百分数-2025.12.31
			if "%" in freq:
				freq = float(freq.replace("%", "")) / 100
	if (float(freq) >= 0.15 and float(freq) <= 0.25) or (float(freq) >= 0.75 and float(freq) < 0.85) or (float(freq) >= 0.85 and float(freq) <= 0.90):
		if (float(freq) >= 0.15 and float(freq) <= 0.25):
			return "loss1"
		elif (float(freq) >= 0.75 and float(freq) < 0.85):
			return "loss2"
		else:
			return "more"
	else:
		return False
jinja2.filters.FILTERS["zsly_chemo_var_v3"] = zsly_chemo_var_v3

# 化疗判断所有知识库有记录位点频率-2025.04.02
def zsly_chemo_summary(info):
	chemo_list = info[0]
	all_detect_result = info[1]
	result = False
	for var in chemo_list:
		if zsly_chemo_var([var, all_detect_result]):
			result = True
			break
	return result
jinja2.filters.FILTERS["zsly_chemo_summary"] = zsly_chemo_summary

# 化疗判断所有知识库有记录位点频率-2025.04.27
def zsly_chemo_summary_v2(info):
	chemo_list = info[0]
	all_detect_result = info[1]
	result = []
	for var in chemo_list:
		if zsly_chemo_var_v3([var, all_detect_result]) and zsly_chemo_var_v3([var, all_detect_result]) in ["loss1", "loss2"]:
			if "loss1" not in result:
				result.append("loss")
		if zsly_chemo_var_v3([var, all_detect_result]) and zsly_chemo_var_v3([var, all_detect_result]) == "more":
			if "more" not in result:
				result.append("more")
	return result
jinja2.filters.FILTERS["zsly_chemo_summary_v2"] = zsly_chemo_summary_v2

# MSI治疗方案仅展示A级并去重-2025.04.03
def zsly_msi_regimen(regimen_list):
	result = []
	for regimen in regimen_list:
		if regimen["evi_conclusion_simple"] == "A" and regimen["regimen_name"] not in result:
			result.append(regimen["regimen_name"])
	return result
jinja2.filters.FILTERS["zsly_msi_regimen"] = zsly_msi_regimen

# 对治疗方案进行整合/去重处理-2025.04.08
# 相关药物、相应类型相同的，取证据等级最高的
def zsly_evi_redup(evi_list):
	# 每个治疗方案/预后/诊断汇总成字典
	predictive_result = {}
	prognostic_result = {}
	diagnostic_result = {}
	for evi in evi_list:
		if evi["evidence_type"] == "Predictive":
			if evi["regimen_name"] not in predictive_result.keys():
				predictive_result.setdefault(evi["regimen_name"], [])
			predictive_result[evi["regimen_name"]].append(evi)
		if evi["evidence_type"] == "Prognostic":
			if "Prognostic" not in prognostic_result.keys():
				prognostic_result.setdefault("Prognostic", [])
			prognostic_result["Prognostic"].append(evi)
		if evi["evidence_type"] == "Diagnostic":
			#print ("111111111111", evi)
			if "Diagnostic" not in diagnostic_result.keys():
				diagnostic_result.setdefault("Diagnostic", [])
			diagnostic_result["Diagnostic"].append(evi)
	#print ("11111111111111111", predictive_result)

	# 获取每个治疗方案/预后/诊断最高等级的证据，汇总为列表
	def get_top_evi(evi_dict):
		top_level_evi = []
		for regimen, evi in evi_dict.items():
			evi_concolusion_list = [i["evi_conclusion_simple"] for i in evi]
			top_level = "A" if "A" in evi_concolusion_list else \
						"B" if "B" in evi_concolusion_list else \
						"C" if "C" in evi_concolusion_list else \
						"D" if "D" in evi_concolusion_list else \
						"E"
			top_level_evi.extend([i for i in evi if i["evi_conclusion_simple"] == top_level])
		return top_level_evi
	#print ("111111111111111111", predictive_result)
	predictive_top_evi = get_top_evi(predictive_result) if predictive_result else []
	prognostic_top_evi = get_top_evi(prognostic_result) if prognostic_result else []
	diagnostic_top_evi = get_top_evi(diagnostic_result) if diagnostic_result else []
	#print ("1111111111111111111111111", diagnostic_top_evi)
		
	# 相关药物、相应类型、证据来源、证据等级相同的进行合并
	def get_evi_sum(evi_list):
		tmp_dict = {}
		for evi in evi_list:
			key = (evi["regimen_name"], evi["clinical_significance_cn"], evi["evi_origin"], evi["evi_conclusion_simple"])
			if key not in tmp_dict.keys():
				tmp_dict.setdefault(key, [])
			tmp_dict[key].append(evi)
		tmp_result = []
		for key, evi in tmp_dict.items():
			tumor_name_cn_list = []
			# 2025.05.16-新增参考文献
			evi_refer_type = {}
			# 2025.05.16-新增完成
			for i in evi:
				#print (i["regimen_name"],i["tumor_name_cn"])
				if i["tumor_name_cn"] not in tumor_name_cn_list:
					tumor_name_cn_list.append(i["tumor_name_cn"])	
				# 2025.05.16-新增参考文献
				#print (i["evi_refer_dict"])
				for k, v in i["evi_refer_dict"].items():
					if k not in evi_refer_type.keys():
						evi_refer_type.setdefault(k, [])
					evi_refer_type[k].extend(v)
				# 增加一个癌种列表-用来展示辅助诊断支持证据的-2025.05.30
				dia_tumor_list = []
				for j in re.split("、", i["tumor_name_list"]):
					if j not in dia_tumor_list:
						dia_tumor_list.append(j)
				# 2025.05.30-新增完成

			evi_refer_list = []
			for k, v in evi_refer_type.items():
				evi_refer_list.append({
					"appr_type" : k,
					"refer_code" : v
				})
			#print (evi_refer_list)
			# 2025.05.16-新增完成
			tmp_result.append({
				"regimen_name" : key[0],
				"clinical_significance_cn" : key[1],
				"evi_origin" : key[2],
				"evi_conclusion_simple" : key[3],
				"sense_rule" : "".join(list(set(i["sense_rule"] for i in evi))),
				"regimen_name_py" : "".join(list(set(i["regimen_name_py"] for i in evi))),
				"tumor_name_cn" : "、".join(tumor_name_cn_list),
				"evidence_type" : "".join(list(set(i["evidence_type"] for i in evi))),
				# 2025.05.16-新增参考文献
				"evi_refer_type" : evi_refer_list,
				# 2025.05.16-新增完成
				# 2025.05.30-新增诊断支持证据癌种
				"dia_tumor_list" : "、".join(dia_tumor_list)
				# 2025.05.30-新增完成
			})
		return tmp_result
	
	result = []
	if predictive_top_evi:
		#print ("1111111111111", predictive_top_evi)
		result.extend(get_evi_sum(predictive_top_evi))
	if prognostic_top_evi:
		result.extend(get_evi_sum(prognostic_top_evi))
	if diagnostic_top_evi:
		#print (get_evi_sum(diagnostic_top_evi))
		result.extend(get_evi_sum(diagnostic_top_evi))
	#print ("11111111111", result)
	return sorted(result, key=lambda i : (i["evi_conclusion_simple"], i["sense_rule"], i["regimen_name_py"].upper()))
jinja2.filters.FILTERS["zsly_evi_redup"] = zsly_evi_redup

# 非PMID参考文献序号由上个表格继续 - 2025.04.24
def zsly_refer_id(refer):
	pmid_list = []
	for i in refer:
		pmid = (re.split("PMID:", i)[-1]).replace("]", "")
		if pmid not in pmid_list:
			pmid_list.append(pmid)
	return len(pmid_list) + 1
jinja2.filters.FILTERS["zsly_refer_id"] = zsly_refer_id

#--------------------中山六院CP200相关 结束------------------------------------------------------#

# 国际部-不匹配用药的情况下，只能通过clinical_category进行变异分类-2025.03.18
def International_get_I_II(var_list):
	I_II_var = [var for var in var_list if var["clinical_category"] in ["Tier I", "Tier II"]]
	return sorted(I_II_var, key = lambda i : i["clinical_category"])
jinja2.filters.FILTERS["International_get_I_II"] = International_get_I_II

def International_get_III(var_list):
	III_var = [var for var in var_list if var["clinical_category"] in ["Tier III"]]
	return sorted(III_var, key = lambda i : i["clinical_category"])
jinja2.filters.FILTERS["International_get_III"] = International_get_III

# 保留4位数小数-2025.03.19
def xajdy_decimal_float(a):
	'''
	四舍五入，处理浮点数
	'''
	#return str(Decimal(str(a)).quantize(Decimal("0.01"), rounding="ROUND_HALF_UP"))

	# 更新为<5舍去，>=5进1，和excel数据一致，但有可能和系统展示的不同-2023.07.04
	return "{:.4f}".format(float(int(float(a) * 10000 + 0.5) / 10000))
jinja2.filters.FILTERS["xajdy_decimal_float"] = xajdy_decimal_float

# 国际部MP化疗检测结果需要拆分为两列展示-2025.03.19
# num指的是模板中要按几列展示，基因+TPM为一列
def international_chemo(chemo_data):
	num = 1
	for i in chemo_data:
		i["num"] = num
		num += 1
	chemo_list1 = []
	chemo_list2 = []
	if len(chemo_data) % 2 == 1:
		chemo_data.append({"num" : len(chemo_data) + 1})
		mid_index = len(chemo_data) // 2
		chemo_list1 = chemo_data[:mid_index]
		chemo_list2 = chemo_data[mid_index:]
	result = []
	if chemo_list1:
		for i in range(0, len(chemo_list1)):
			result.append({
				"num_1" : chemo_list1[i]["num"] if chemo_list1[i]["num"] else "",
				"gene_symbol_1" : chemo_list1[i]["gene_symbol"] if chemo_list1[i]["gene_symbol"] else "",
				"dbsnp_1" : chemo_list1[i]["dbsnp"] if chemo_list1[i]["dbsnp"] else "",
				"genotype_1" : chemo_list1[i]["genotype"] if chemo_list1[i]["genotype"] else "",
				"num_2" : chemo_list2[i]["num"] if "gene_symbol" in chemo_list2[i].keys() else "",
				"gene_symbol_2" : chemo_list2[i]["gene_symbol"] if "gene_symbol" in chemo_list2[i].keys() else "",
				"dbsnp_2" : chemo_list2[i]["dbsnp"] if "gene_symbol" in chemo_list2[i].keys() else "",
				"genotype_2" : chemo_list2[i]["genotype"] if "gene_symbol" in chemo_list2[i].keys() else "",	
			})
	return result
jinja2.filters.FILTERS["international_chemo"] = international_chemo

# 保留2位数小数-2025.03.19
def international_decimal_float(a):
	'''
	四舍五入，处理浮点数
	'''
	#return str(Decimal(str(a)).quantize(Decimal("0.01"), rounding="ROUND_HALF_UP"))

	# 更新为<5舍去，>=5进1，和excel数据一致，但有可能和系统展示的不同-2023.07.04
	return "{:.2f}".format(float(int(float(a) * 100 + 0.5) / 100))
jinja2.filters.FILTERS["international_decimal_float"] = international_decimal_float

# 保留3位数小数-2025.03.19
def xajdy_mp_decimal_float(a):
	'''
	四舍五入，处理浮点数
	'''
	#return str(Decimal(str(a)).quantize(Decimal("0.01"), rounding="ROUND_HALF_UP"))

	# 更新为<5舍去，>=5进1，和excel数据一致，但有可能和系统展示的不同-2023.07.04
	return "{:.3f}".format(float(int(float(a) * 1000 + 0.5) / 1000))
jinja2.filters.FILTERS["xajdy_mp_decimal_float"] = xajdy_mp_decimal_float

# 复旦华山116-结果总结-2025.03.30
def fdhs_116_var_sum(var_list):
	'''
	复旦华东116&76，报告中需要展示I/II/肿瘤发生发展相关/III类变异详情/胚系
	I类+II类：XX基因X号外显子/内含子检测到XX突变hgvs_p/hgvs_c，突变丰度为X
	III类：   XX基因X号外显子/内含子检测到XX突变hgvs_p/hgvs_c，突变丰度为X

	snvindel：XX基因X号外显子/内含子XX突变hgvs_c,hgvs_p,丰度为X；
	   胚系： XX基因X号外显子/内含子XX突变hgvs_c,hgvs_p,基因型为X
	cnv：XX基因扩增，拷贝数为X；
	sv：XX-XX基因融合，风度为X；
	'''
	def var_info_stran(var):
		var_info = ""
		region_dict = {
		"exon" : "外显子",
		"intron" : "内含子",
		"3'UTR" : "3'UTR",
		"5'UTR" : "5'UTR",
		"3'FLANKING" : "非编码区",
		"5'FLANKING" : "非编码区"
		}
		if var["bio_category"] == "Cnv":
			var_info = var["gene_symbol"]+"基因扩增，拷贝数为"+str(var["cn_mean"])
		elif var["bio_category"] == "Sv":
			var_info = var["five_prime_gene"]+"-"+var["three_prime_gene"]+"基因融合，丰度为"+str(var["freq_str"])
		elif var["bio_category"] == "Snvindel":
			freq = ""
			if var["var_origin"] == "germline":
				if "freq_rc" in var.keys() and var["freq_rc"] and float(var["freq_rc"]) >= 0.85:
					freq = "基因型为纯合型"
				elif "freq_rc" in var.keys() and var["freq_rc"] and float(var["freq_rc"]) < 0.85:
					freq = "基因型为杂合型"
				else:
					freq = "基因型为【手动填写！】"
			else:
				freq = "丰度为{0}".format(var["freq_str"])
			# 加一个MET exon 14跳跃突变-2023.07.13
			if var["gene_symbol"] == "MET" and "exon14 skipping" in var["variant_interpret_cn"]:
				if var["hgvs_p"] != "p.?":
					var_info  = "MET基因14号外显子跳跃突变"+var["hgvs_c"]+"，"+var["hgvs_p"]+"，" + freq
				else:
					var_info = "MET基因14号外显子跳跃突变"+var["hgvs_c"]+"，" + freq
			# MET exon14跳跃添加结束-2023.07.13
			else:
				region_list_en = re.split("_", var["gene_region"])
				region_list_cn = []
				for i in region_list_en:
					if re.search("exon", i):
						region_list_cn.append(i.replace("exon", "")+"号外显子")
					elif re.search("intron", i):
						region_list_cn.append(i.replace("intron", "")+"号内含子")
					else:
						region_cn = region_dict[i] if i in region_dict.keys() else i
						region_list_cn.append(region_cn)
				if var["hgvs_p"] != "p.?":
					var_info = var["gene_symbol"]+"基因"+"到".join(region_list_cn) + var["type_cn"] +  var["hgvs_c"] + " " + var["hgvs_p"]+"，" + freq
				else:
					if var["type_cn"] != "--":
						var_info = var["gene_symbol"]+"基因"+"到".join(region_list_cn) + var["type_cn"] + var["hgvs_c"]+"突变，" + freq
					else:
						var_info = var["gene_symbol"]+"基因"+"到".join(region_list_cn) + var["hgvs_c"]+"突变，" + freq
		return var_info
	return "；".join([var_info_stran(var) for var in var_list])
jinja2.filters.FILTERS["fdhs_116_var_sum"] = fdhs_116_var_sum

# 济宁医学院附属-报告中拆分为LC10和其他基因 - 2025.03.28
def jnfy_filter_lc10(info):
	var_list = info[0]
	var_type = info[1]
	LC10_snvindel = ["ALK", "BRAF", "EGFR", "ERBB2", "KRAS", "MET", "NRAS", "PIK3CA", "RET", "ROS1"]
	LC10_sv = ["ALK", "RET", "ROS1", "MET"]
	LC10_cnv = ["MET"]
	lc10_var = []
	other_var = []
	for var in var_list:
		if (var["bio_category"] == "Snvindel" and var["gene_symbol"] in LC10_snvindel) or \
		   (var["bio_category"] == "Sv" and set(re.split(",", var["gene_symbol"])) & set(LC10_sv)) or \
		   (var["bio_category"] == "Cnv" and var["gene_symbol"] in LC10_cnv):
			lc10_var.append(var)
		else:
			other_var.append(var)
	if var_type == "lc10":
		return lc10_var
	else:
		return other_var
jinja2.filters.FILTERS["jnfy_filter_lc10"] = jnfy_filter_lc10

# 北京同仁tHRR-筛选出BRCA相关变异-2025.04.07
def bjtr_thrr_brca_var(var_list):
	return [var for var in var_list if var["gene_symbol"] in ["BRCA1", "BRCA2"]]
jinja2.filters.FILTERS["bjtr_thrr_brca_var"] = bjtr_thrr_brca_var

# 北京同仁tHRR-可能获益药物仅展示BRCA相关-2025.04.07
def bjtr_thrr_regimen(regimen_list):
	result = []
	for regimen in regimen_list:
		for var in regimen["var"]:
			if var["gene_symbol"] in ["BRCA1", "BRCA2"] and regimen not in result:
				result.append(regimen)
	return result
jinja2.filters.FILTERS["bjtr_thrr_regimen"] = bjtr_thrr_regimen

# 福建省肿瘤HRD-需要分为BRCA、HRR和其他基因，并且按胚系变异等级来展示-2025.04.08
def fjzl_hrd_var_filter(info):
	var_list = info[0]
	var_type = info[1]
	brca_gene = ["BRCA1", "BRCA2"]
	hrr_gene = ["ATM", "BARD1", "BRIP1",  "CDK12", "CHEK1", "CHEK2", "FANCA", "FANCL", \
			 	"HDAC2", "PALB2", "PPP2R2A",  "RAD51B", "RAD51C", "RAD51D", "RAD54L"]
	other_gene = ["CDH1", "PTEN", "TP53"]
	gene_list = []
	if var_type == "brca":
		gene_list = brca_gene
	elif var_type == "hrr":
		gene_list = hrr_gene
	else:
		gene_list = other_gene
	level_5 = [var for var in var_list if var["clinic_num_g"] == 5 and var["gene_symbol"] in gene_list]
	level_4 = [var for var in var_list if var["clinic_num_g"] == 4 and var["gene_symbol"] in gene_list]
	level_5 = sorted(level_5, key=lambda i : (i["gene_symbol"], i["freq"]))
	level_4 = sorted(level_4, key=lambda i : (i["gene_symbol"], i["freq"]))
	return level_5 + level_4
jinja2.filters.FILTERS["fjzl_hrd_var_filter"] = fjzl_hrd_var_filter

# 同济BPTM MSI，药物“、”.join()展示为一行，A级还需要展示引用机构-2025.04.09
def tj_bptm_regimen_sum(a):
	result = []
	evi_sum = a["evi_sum"]
	if "evi_split" in evi_sum.keys():
		if "Predictive" in evi_sum["evi_split"].keys():
			for i in evi_sum["evi_split"]["Predictive"]:
				if "regimen_refer_agency" in i.keys() and i["regimen_refer_agency"] and i["evi_conclusion_simple"] == "A":
					regimen_refer_agency = "/".join(re.split(",", i["regimen_refer_agency"]))
					result.append("{0}（{1}，{2}级，{3}）".format(i["regimen_name"], i["clinical_significance_cn"], i["evi_conclusion_simple"], regimen_refer_agency))
				else:
					result.append("{0}（{1}，{2}级，-）".format(i["regimen_name"], i["clinical_significance_cn"], i["evi_conclusion_simple"]))
		if "Prognostic" in evi_sum["evi_split"].keys():
			for i in evi_sum["evi_split"]["Prognostic"]:
				if "regimen_refer_agency" in i.keys() and i["regimen_refer_agency"] and i["evi_conclusion_simple"] == "A":
					regimen_refer_agency = "/".join(re.split(",", i["regimen_refer_agency"]))
					result.append("{0}（{1}，{2}级，{3}}）".format("预后"+i["clinical_significance_cn"], " / ", i["evi_conclusion_simple"], regimen_refer_agency))
				else:
					result.append("{0}（{1}，{2}级，-）".format("预后"+i["clinical_significance_cn"], " / ", i["evi_conclusion_simple"]))
		if "Diagnostic" in evi_sum["evi_split"].keys():
			for i in evi_sum["evi_split"]["Diagnostic"]:
				if "regimen_refer_agency" in i.keys() and i["regimen_refer_agency"] and i["evi_conclusion_simple"] == "A":
					regimen_refer_agency = "/".join(re.split(",", i["regimen_refer_agency"]))
					result.append("{0}（{1}，{2}级，{3}）".format("辅助诊断", " / ", i["evi_conclusion_simple"], regimen_refer_agency))
				else:
					result.append("{0}（{1}，{2}级，-）".format("辅助诊断", " / ", i["evi_conclusion_simple"]))
	return "、".join(result) if result else "-"
jinja2.filters.FILTERS["tj_bptm_regimen_sum"] = tj_bptm_regimen_sum 

# 同济BPTM MSI，药物“、”.join()展示为一行，A级还需要展示引用机构-2025.04.23
# 无引用机构的不显示“-”
def tj_bptm_regimen_sum_v2(a):
	result = []
	evi_sum = a["evi_sum"]
	if "evi_split" in evi_sum.keys():
		if "Predictive" in evi_sum["evi_split"].keys():
			for i in evi_sum["evi_split"]["Predictive"]:
				if "regimen_refer_agency" in i.keys() and i["regimen_refer_agency"] and i["evi_conclusion_simple"] == "A":
					regimen_refer_agency = "/".join(re.split(",", i["regimen_refer_agency"]))
					result.append("{0}（{1}，{2}级，{3}）".format(i["regimen_name"], i["clinical_significance_cn"], i["evi_conclusion_simple"], regimen_refer_agency))
				else:
					result.append("{0}（{1}，{2}级）".format(i["regimen_name"], i["clinical_significance_cn"], i["evi_conclusion_simple"]))
		if "Prognostic" in evi_sum["evi_split"].keys():
			for i in evi_sum["evi_split"]["Prognostic"]:
				if "regimen_refer_agency" in i.keys() and i["regimen_refer_agency"] and i["evi_conclusion_simple"] == "A":
					regimen_refer_agency = "/".join(re.split(",", i["regimen_refer_agency"]))
					result.append("{0}（{1}，{2}级，{3}}）".format("预后"+i["clinical_significance_cn"], " / ", i["evi_conclusion_simple"], regimen_refer_agency))
				else:
					result.append("{0}（{1}，{2}级）".format("预后"+i["clinical_significance_cn"], " / ", i["evi_conclusion_simple"]))
		if "Diagnostic" in evi_sum["evi_split"].keys():
			for i in evi_sum["evi_split"]["Diagnostic"]:
				if "regimen_refer_agency" in i.keys() and i["regimen_refer_agency"] and i["evi_conclusion_simple"] == "A":
					regimen_refer_agency = "/".join(re.split(",", i["regimen_refer_agency"]))
					result.append("{0}（{1}，{2}级，{3}）".format("辅助诊断", " / ", i["evi_conclusion_simple"], regimen_refer_agency))
				else:
					result.append("{0}（{1}，{2}级）".format("辅助诊断", " / ", i["evi_conclusion_simple"]))
	return "、".join(result) if result else "-"
jinja2.filters.FILTERS["tj_bptm_regimen_sum_v2"] = tj_bptm_regimen_sum_v2

# 2025.04.23-贵州人民BRCA-snvindel
def gzrm_brca_var(var):
	if var["hgvs_p"] != "p.?":
		return "{0} {1} {2}".format(var["gene_region"], var["hgvs_c"], var["hgvs_p"])
	else:
		return "{0} {1}".format(var["gene_region"], var["hgvs_c"])
jinja2.filters.FILTERS["gzrm_brca_var"] = gzrm_brca_var

# 2025.05.06-华西150，变异解读需要根据。进行拆分
def hx_150_split_varinter(var_list):
	for var in var_list:
		num = 1
		var["inter_split"] = []
		if var["variant_interpret_cn"]:
			varinter_list = re.split("。", var["variant_interpret_cn"])
			for i in varinter_list:
				if i:
					var["inter_split"].append({
						"num" : str(num),
						"inter" : i
					})
				num += 1
	return var_list
jinja2.filters.FILTERS["hx_150_split_varinter"] = hx_150_split_varinter

# 2025.05.06-华西150-风险和管理-2025.05.06
# 2025.10.09-新增需求-风险管理需要根据性别做筛选
def hx_150_risk(info):
	var_list = info[0]
	hx_150_gene_risk_inter = info[1]
	hx_150_gene_risk_table = info[2]
	hx_150_gene_risk_management = info[3]
	gender = info[4]
	for var in var_list:
		var["hx_150_gene_risk_inter"] = hx_150_gene_risk_inter.get(var["gene_symbol"], [])
		var["hx_150_gene_risk_table"] = hx_150_gene_risk_table.get(var["gene_symbol"], [])
		# 2025.09.24-新增需求，DIS3L2、MBD4、MSH3、TRIM37基因，如果检出变异为杂合，就不展示风险管理
		if var["gene_symbol"] in ["DIS3L2", "MBD4", "MSH3", "TRIM37"] and float(var["freq"]) <= 0.85:
			var["hx_150_gene_risk_management"] = []
		else:
			var["hx_150_gene_risk_management"] = hx_150_gene_risk_management.get(var["gene_symbol"], [])
		#print ("hx_150_gene_risk_inter", var["hx_150_gene_risk_inter"])
		#print ("hx_150_gene_risk_table", var["hx_150_gene_risk_table"])
		#print ("hx_150_gene_risk_management", var["hx_150_gene_risk_management"])
		# 2025.10.09-新增需求-风险管理需要根据性别做筛选
		var["hx_150_gene_risk_management"] = [i for i in var["hx_150_gene_risk_management"] if i["gender"] == gender or not i["gender"]]
		# 2025.10.09-新增完成
		# 1. 基因与肿瘤的关联
		# 位置1：变异后面
		# 位置2：风险表格下面
		# 位置3：风险管理后面
		# 1.1 若一个基因只有一条描述，则直接放在位置1
		# 1.2 若一个基因有多条描述
		# 1.2.1 变异分组和配置表的变异分组一致的话用变异分组对应的数据，否则用通用数据（ABRAXAS1 R361Q需要展示特殊和通用两条）
		# 1.2.2 上一步骤过滤完剩一条，则直接展示到位置1
		# 1.2.3 剩多条，常显展示到位置1，常隐需要结合变异基因型，纯合展示到位置2，杂合展示到位置3	
		# 增加一个参考文献-2025.05.09
		reference = []	
		gene_type = ""
		if var["type"] in ["Loss", "Gain"]:
			gene_type = "纯合" if var["cnv_type"] == "HomoDel" else "杂合" if var["cnv_type"] == "HeteDel" else "" 
		else:
			gene_type = "纯合" if var["freq"] and float(var["freq"]) >= 0.85 else "杂合"
		var_category_names = var["var_category_names"] if "var_category_names" in var.keys() and var["var_category_names"] else ""
		site1_list = []
		site2_list = []
		site3_list = []
		# 2025.08.15-如果样本检出变异均为常隐+杂合时，需要删除家族其他成员风险管理，这边加个判定
		# 常隐+杂合 添加“yes”，其他添加“no”，模板中判定“no”在列表里，则展示家族其他成员风险管理
		judge_ad_hete = []
		# 2025.08.15-华西150体检没有展示ABRAXAS1基因
		if var["gene_symbol"] == "ABRAXAS1" and var["bio_category"] == "Snvindel" and var["hgvs_p"] == "p.(R361Q)":
			site1_list = [i["inter"] for i in var["hx_150_gene_risk_inter"]]
			reference.extend([i["reference"] for i in var["hx_150_gene_risk_inter"] if i["reference"]])
		else:
			category_list = [i for i in var["hx_150_gene_risk_inter"] if var_category_names and set(re.split(",", var_category_names)) & set(re.split("\|\|", i["var_category_names"]))]
			#print ("category_list", category_list)
			nocategory_list = [i for i in var["hx_150_gene_risk_inter"] if not i["var_category_names"]]
			#print ("nocategory_list", nocategory_list)
			filter_list = category_list if category_list else nocategory_list
			if len(filter_list) == 1:
				#site1_list = [i["inter"] for i in filter_list]
				reference.extend([i["reference"] for i in filter_list if i["reference"]])

				# 2025.08.15-新增-若疾病仅有常染色体隐性遗传，且为杂合(2025.08.18-非杂合也要删除)时，“若您处于育龄期，”前加一句“您仅为携带者，”
				site1_list = []
				if filter_list[0]["genetic_mode"] == "常染色体隐性遗传" and gene_type == "杂合":
					judge_ad_hete.append("yes")
					for i in filter_list:
						if "若您处于育龄期，" in i["inter"]:
							site1_split_str = re.split("若您处于育龄期，", i["inter"])
							site1_list.append(site1_split_str[0] + "您仅为携带者，若您处于育龄期，" + "".join(site1_split_str[1:]))
						else:
							site1_list.append(i["inter"])
				else:
					if filter_list[0]["genetic_mode"] == "常染色体隐性遗传":
						judge_ad_hete.append("yes")
					else:
						judge_ad_hete.append("no")
					site1_list = [i["inter"] for i in filter_list]
				# 2025.08.15-新增完成

			else:
				judge_ad_hete.append("no")
				site1_list = [i["inter"] for i in filter_list if i["genetic_mode"] == "常染色体显性遗传"]
				reference.extend([i["reference"] for i in filter_list if i["reference"] and i["genetic_mode"] == "常染色体显性遗传"])
				site2_list = [i["inter"] for i in filter_list if i["genetic_mode"] == "常染色体隐性遗传"] if gene_type == "纯合" else []
				if gene_type == "纯合":
					reference.extend([i["reference"] for i in filter_list if i["reference"] and i["genetic_mode"] == "常染色体隐性遗传"])
				site3_list = [i["inter"] for i in filter_list if i["genetic_mode"] == "常染色体隐性遗传"] if gene_type == "杂合" or not gene_type else []
		#print ("site1_list", site1_list)
		#print ("site2_list", site2_list)
		#print ("site3_list", site3_list)
		var["hx_150_gene_risk_inter_site1"] = "".join(site1_list)
		var["hx_150_gene_risk_inter_site2"] = "".join(site2_list)
		var["hx_150_gene_risk_inter_site3"] = "".join(site3_list)
		var["judge_ad_hete"] = judge_ad_hete
		#print ("hx_150_gene_risk_inter_site1", var["hx_150_gene_risk_inter_site1"] )
		#print ("hx_150_gene_risk_inter_site2", var["hx_150_gene_risk_inter_site2"] )
		#print ("hx_150_gene_risk_inter_site3", var["hx_150_gene_risk_inter_site3"] )
		reference_redup = []
		for i in reference:
			if i not in reference_redup:
				reference_redup.append(i)
		var["reference"] = ",".join(reference_redup)

		# 2. 风险表格
		# 样式1：三列，包含肿瘤、绝对风险和与疾病的相关性
		# 样式2：四列，包含肿瘤、年龄、携带者风险和普通人风险
		# 样式3：总结一句话
		risk_table_tmp_dict = {}
		if var["hx_150_gene_risk_table"]:
			hx_150_gene_risk_table_filter_category = [i for i in var["hx_150_gene_risk_table"] if var_category_names and set(re.split(",", var_category_names)) & set(re.split("\|\|", i["var_category_names"]))]
			hx_150_gene_risk_table_filter_nocategory = [i for i in var["hx_150_gene_risk_table"] if not i["var_category_names"]]
			hx_150_gene_risk_table_filter = hx_150_gene_risk_table_filter_category if hx_150_gene_risk_table_filter_category else hx_150_gene_risk_table_filter_nocategory
			for i in hx_150_gene_risk_table_filter:
				if i["absolute_risk"]:
					i["risk_table_type"] = "column_3"
				elif i["summary"]:
					i["risk_table_type"] = "summary"
				else:
					i["risk_table_type"] = "column_4"
				note = i["note"] if i["note"] else "nofound"
				if note not in risk_table_tmp_dict.keys():
					risk_table_tmp_dict.setdefault(note, [])
				risk_table_tmp_dict[note].append(i)
		risk_table_list = []
		for k, v in risk_table_tmp_dict.items():
			risk_table_type = v[0].get("risk_table_type", "nofound") if v and len(v) >= 1 else ""
			risk_table_list.append({
				"note" : k,
				"risk_table_type" : risk_table_type,
				"info" : v
			})
		#print ("risk_table_list", risk_table_list)
		var["risk_table_list"] = risk_table_list

		# 3. 风险管理表格
		# 样式1：不区分人群
		# 样式2：需要区分儿童和成年人
		hx_150_gene_risk_management_dict = {}
		if var["hx_150_gene_risk_management"]:
			hx_150_gene_risk_management_filter_category = [i for i in var["hx_150_gene_risk_management"] if var_category_names and set(re.split(",", var_category_names)) & set(re.split("\|\|", i["var_category_names"]))]
			hx_150_gene_risk_management_filter_nocategory = [i for i in var["hx_150_gene_risk_management"] if not i["var_category_names"]]
			hx_150_gene_risk_management_filter = hx_150_gene_risk_management_filter_category if hx_150_gene_risk_management_filter_category else hx_150_gene_risk_management_filter_nocategory
			#tumor_sort = []
			for i in hx_150_gene_risk_management_filter:
				if i["note3"]:
					i["management_type"] = "summary"
				else:
					i["management_type"] = "table"
			#	if i["tumor"] not in tumor_sort:
			#		tumor_sort.append(i["tumor"])
				key = i["note2"] if i["note2"] else "noinfo"
				if key not in hx_150_gene_risk_management_dict.keys():
					hx_150_gene_risk_management_dict.setdefault(key,{})
				if i["tumor"] not in hx_150_gene_risk_management_dict[key].keys():
					hx_150_gene_risk_management_dict[key].setdefault(i["tumor"],[])
				hx_150_gene_risk_management_dict[key][i["tumor"]].append(i)
		hx_150_gene_risk_management_list = []
		for agetype, v in hx_150_gene_risk_management_dict.items():
			tmp_list = []
			for tumor, info in v.items():
				tmp_list.append(
					{
						"tumor" : tumor,
						"info" : info
					}
				)
			hx_150_gene_risk_management_list.append(
				{
					"agetype" : agetype,
					"management_list" : tmp_list,
					"management_type" : tmp_list[0]["info"][0]["management_type"] if tmp_list and \
										len(tmp_list) >= 1 and \
										"info" in tmp_list[0] and \
										tmp_list[0]["info"] and \
										len(tmp_list[0]["info"]) >= 1 and \
										"management_type" in tmp_list[0]["info"][0].keys() and \
										tmp_list[0]["info"][0]["management_type"] else ""
				}
			)
		var["hx_150_gene_risk_management_list"] = hx_150_gene_risk_management_list
	return var_list
jinja2.filters.FILTERS["hx_150_risk"] = hx_150_risk

# 2025.05.07-华西150-风险管理去重
def hx_150_redup_risk_management(var_list):
	result = []
	for var in var_list:
		if var["hx_150_gene_risk_management_list"] and var["hx_150_gene_risk_management_list"] not in result:
			result.append(var["hx_150_gene_risk_management_list"])
	return result
jinja2.filters.FILTERS["hx_150_redup_risk_management"] = hx_150_redup_risk_management

# 2025.05.07-华西150-汇总有风险管理的变异
def hx_150_redup_risk_management_count(var_list):
	result = []
	for var in var_list:
		if var["hx_150_gene_risk_management_list"]:
			result.append(var)
	return result
jinja2.filters.FILTERS["hx_150_redup_risk_management_count"] = hx_150_redup_risk_management_count

# 2025.05.07-华西150-风险管理下面的描述去重
def hx_150_redup_risk_management_inter(var_list):
	result = []
	for var in var_list:
		if var["hx_150_gene_risk_inter_site3"] and var["hx_150_gene_risk_inter_site3"] not in result:
			result.append(var["hx_150_gene_risk_inter_site3"])
	return result
jinja2.filters.FILTERS["hx_150_redup_risk_management_inter"] = hx_150_redup_risk_management_inter

# BPTM plus展示变异小结-2025.05.07
# 仅包含snvindel变异
def bptm_plus_var_sum_whxh(var_list):
	result = []
	for var in var_list:
		if var["hgvs_p_abbr"] != "p.?":
			result.append(var["gene_symbol"]+" "+var["hgvs_c"]+" "+var["hgvs_p_abbr"])
		else:
			result.append(var["gene_symbol"]+" "+var["hgvs_c"])
	return ", ".join(result)
jinja2.filters.FILTERS["bptm_plus_var_sum_whxh"] = bptm_plus_var_sum_whxh

# 浙二150变异小结-2025.05.07
# 包含snvindel和CNV变异
def zjey_150_var_sum(var_list):
	result = []
	for var in var_list:
		if var["type"] == "Loss":
			result.append(var["gene_symbol"] + " " + var["value"] + " del")
		elif var["type"] == "Gain":
			result.append(var["gene_symbol"] + " " + var["value"] + " dup")
		else:
			if var["hgvs_p"] != "p.?":
				result.append(var["gene_symbol"]+" "+var["hgvs_c"]+" "+var["hgvs_p"])
			else:
				result.append(var["gene_symbol"]+" "+var["hgvs_c"])
	return ", ".join(result)
jinja2.filters.FILTERS["zjey_150_var_sum"] = zjey_150_var_sum

# 浙二变异区分林奇和其他-2025.05.27
# 2025.05.29-更新，添加指南推荐基因
def zjey_filter_var(info):
	var_list = info[0]
	var_type = info[1]
	gls5_gene = ["MLH1", "MSH2", "MSH6", "PMS2", "EPCAM"]
	nccn_gene = ["APC", "ASXL1", "ATM", "AXIN2", "BAP1", "BARD1","BLM","BMPR1A","BRCA1","BRCA2", \
				 "BRIP1","CDH1","CDK4","CDKN1B","CDKN1C","CDKN2A","CHEK2","DDB2","DDX41","DICER1", \
				 "DIS3L2","EGFR","ERCC2","ERCC3","ERCC4","ERCC5","FANCA","FANCB","FANCC", \
				 "FANCD2","FANCE","FANCF","FANCG","FANCI","FANCL","FANCM","FH","FLCN","GREM1", \
				 "HOXB13","KIT","MAD2L2","MAX","MBD4","MC1R","MEN1","MET","MITF", \
				 "MLH3","MSH3","MUTYH","NF1","NTHL1","PALB2","PDGFRA", \
				 "POLD1","POLE","POLH","PRSS1","PTEN","RAD51","RAD51C","RAD51D","RET","RHBDF2", \
				 "RNF43","RPS20","SDHA","SDHAF2","SDHB","SDHC","SDHD","SLX4","SMAD4","SMARCA4", \
				 "SPINK1","STK11","TERT","TMEM127","TP53","TRIM37","TRIP13","TSC1","TSC2","UBE2T", \
				 "VHL","WT1","XPA","XPC","XRCC2","GALNT12"]
	gls5_var = [var for var in var_list if var["gene_symbol"] in gls5_gene]
	nccn_var = [var for var in var_list if var["gene_symbol"] in nccn_gene]
	other_var = [var for var in var_list if var["gene_symbol"] not in gls5_gene and var["gene_symbol"] not in nccn_gene]
	if var_type == "gls5":
		return gls5_var
	elif var_type == "nccn":
		return nccn_var
	else:
		return other_var
jinja2.filters.FILTERS["zjey_filter_var"] = zjey_filter_var

# 吉大一10基因-结果小结出KRAS/NRAS/BRAF特殊处理-仅用药，预后和辅助诊断不变-2025.05.08
# 仅适用肠癌-2025.05.19
def jdy_10_summary_ras(info):
	var_list = info[0]
	tumor_list = info[1]
	result_var = []
	for var in var_list:
		if "肠癌" in tumor_list:
		# 第一遍过滤，过滤掉NRAS/KRAS 肠癌和胰腺癌时的其他癌种相关药物（仅敏感）--取消过滤2025.05.15
		#if var["bio_category"] == "Snvindel" and var["gene_symbol"] in ["KRAS", "NRAS"] and set(["肠癌", "胰腺癌"]) & set(tumor_list):
		#	if "Predictive" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Predictive"]:
		#		var["evi_sum"]["evi_split"]["Predictive"] = [evi for evi in var["evi_sum"]["evi_split"]["Predictive"] if \
		#										 ((set(tumor_list) & set(re.split("\/", evi["tumor_name_cn"])) or "实体瘤" in evi["tumor_name_cn"]) and \
		#	  									 evi["clinical_significance_cn"] == "敏感") or \
		#										 (evi["clinical_significance_cn"] == "耐药")]
		#		print ("1. ",len(var["evi_sum"]["evi_split"]["Predictive"]))
		# 第二遍过滤，KRAS、NRAS、BRAF过滤掉呋喹替尼、贝伐单抗、瑞戈非尼（敏感+耐药都是）（KNB也是，但KNB在模板里固定展示，这边不再处理了）
			if var["bio_category"] == "Snvindel" and var["gene_symbol"] in ["KRAS", "NRAS", "BRAF"]:
				if "Predictive" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Predictive"]:
					var["evi_sum"]["evi_split"]["Predictive"] = [evi for evi in var["evi_sum"]["evi_split"]["Predictive"] if evi["regimen_name"] not in ["呋喹替尼", "贝伐珠单抗", "瑞戈非尼"]]
		#		print ("2. ",len(var["evi_sum"]["evi_split"]["Predictive"]))
		# 第三遍过滤，KRAS 除了G12C有药，其他变异仅推荐D级用药（仅敏感）-取消过滤2025.05.15
		#if var["bio_category"] == "Snvindel" and var["gene_symbol"] == "KRAS" and var["hgvs_p"] != "p.(G12C)":
		#	if "Predictive" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Predictive"]:
		#		var["evi_sum"]["evi_split"]["Predictive"] = [evi for evi in var["evi_sum"]["evi_split"]["Predictive"] if \
		#										 (evi["evi_conclusion_simple"] == "D" and evi["clinical_significance_cn"] == "敏感") or \
		#										 (evi["clinical_significance_cn"] == "耐药")]
		#		print ("3. ",len(var["evi_sum"]["evi_split"]["Predictive"]))
		# 2025.05.15-新增过滤条件-耐药相关-KRAS/NRAS/BRAF耐药只保留A级证据
			if var["bio_category"] == "Snvindel" and var["gene_symbol"] in ["KRAS", "NRAS", "BRAF"]:
				if "Predictive" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Predictive"]:
					var["evi_sum"]["evi_split"]["Predictive"] = [evi for evi in var["evi_sum"]["evi_split"]["Predictive"] if evi["clinical_significance_cn"] == "敏感" or (evi["clinical_significance_cn"] == "耐药" and evi["evi_conclusion_simple"] == "A")]
		result_var.append(var)
	return result_var
jinja2.filters.FILTERS["jdy_10_summary_ras"] = jdy_10_summary_ras

# 吉林大学第一院LC10-返回敏感结果-2025.05.08
def jlyy_10gene_result_sense_result_v2(info):
	# 与上一版相比区别在
	# 增加KNB野生型
	# 适配结果小结，肺癌若检出MET扩增和EGFR，则EGFR敏感药物不体现
	var_list = info[0]
	knb = info[1]
	tumor_list = info[2]
	# 判断是否检出MET 扩增
	judge_met_cnv = "F"
	for var in var_list:
#		print (var["bio_category"], var["gene_symbol"])
		if var["bio_category"] == "Cnv" and var["gene_symbol"] == "MET":
			judge_met_cnv = "T"
			break
	#print ("judge_met_cnv", judge_met_cnv)
	result = []
	#D_sense_regimen = []
	for var in var_list:
		#if var["bio_category"] == "Snvindel" and var["gene_symbol"] == "EGFR" and judge_met_cnv == "T":
		#	continue
		if "Predictive" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Predictive"]:
			if not (var["bio_category"] == "Snvindel" and var["gene_symbol"] == "EGFR" and judge_met_cnv == "T" and "肺癌" in tumor_list):
				regimen_list = []
				A_regimen = ["{0}（{1}）".format(evi["regimen_name"], evi["evi_conclusion_simple"]) for evi in var["evi_sum"]["evi_split"]["Predictive"] if \
							evi["evi_conclusion_simple"] == "A" and evi["clinical_significance_cn"] == "敏感"]
				B_regimen = ["{0}（{1}）".format(evi["regimen_name"], evi["evi_conclusion_simple"]) for evi in var["evi_sum"]["evi_split"]["Predictive"] if \
							evi["evi_conclusion_simple"] == "B" and evi["clinical_significance_cn"] == "敏感"]
				C_regimen = ["{0}（{1}）".format(evi["regimen_name"], evi["evi_conclusion_simple"]) for evi in var["evi_sum"]["evi_split"]["Predictive"] if \
							evi["evi_conclusion_simple"] == "C" and evi["clinical_significance_cn"] == "敏感"]
				D_regimen = ["{0}（{1}）".format(evi["regimen_name"], evi["evi_conclusion_simple"]) for evi in var["evi_sum"]["evi_split"]["Predictive"] if \
							evi["evi_conclusion_simple"] == "D" and evi["clinical_significance_cn"] == "敏感"]
				regimen_list = A_regimen if A_regimen else \
							   B_regimen if B_regimen else \
							   C_regimen if C_regimen else []
				if regimen_list:
					var["regimen_sum_for_jdyy"] = "，".join(regimen_list)
					result.append(var)
	result = sorted(result, key = lambda i:i["gene_symbol"])
	if knb:
		knb["regimen_sum_for_jdyy"] = "帕尼单抗，西妥昔单抗，西妥昔单抗β+FOLFIRI"
		result.insert(0, knb)

	return result
jinja2.filters.FILTERS["jlyy_10gene_result_sense_result_v2"] = jlyy_10gene_result_sense_result_v2

def jdy_evi_sum(Predictive, filter_regimen):
	regimen_list = []
	A_regimen = ["{0}（{1}）".format(evi["regimen_name"], evi["evi_conclusion_simple"]) for evi in Predictive if \
				evi["evi_conclusion_simple"] == "A" and evi["clinical_significance_cn"] == "敏感" and evi["regimen_name"] not in filter_regimen]
	B_regimen = ["{0}（{1}）".format(evi["regimen_name"], evi["evi_conclusion_simple"]) for evi in Predictive if \
				evi["evi_conclusion_simple"] == "B" and evi["clinical_significance_cn"] == "敏感" and evi["regimen_name"] not in filter_regimen]
	C_regimen = ["{0}（{1}）".format(evi["regimen_name"], evi["evi_conclusion_simple"]) for evi in Predictive if \
				evi["evi_conclusion_simple"] == "C" and evi["clinical_significance_cn"] == "敏感" and evi["regimen_name"] not in filter_regimen]
	D_regimen = ["{0}（{1}）".format(evi["regimen_name"], evi["evi_conclusion_simple"]) for evi in Predictive if \
				evi["evi_conclusion_simple"] == "D" and evi["clinical_significance_cn"] == "敏感" and evi["regimen_name"] not in filter_regimen]
	regimen_list = A_regimen if A_regimen else \
				   B_regimen if B_regimen else \
				   C_regimen if C_regimen else []
	return regimen_list

# 吉林大学第一院LC10-返回敏感结果-2025.06.17
def jlyy_10gene_result_sense_result_v3(info):
	# 与上一版相比区别在
	# 增加KNB野生型
	# 适配结果小结，肺癌若检出MET扩增和EGFR，则EGFR敏感药物不体现
	# 增加耐药抵扣-2025.06.17：若检出EGFR T790M，则L858R、19del、L861Q、G719、S768I不展示阿法替尼、达可替尼、厄洛替尼、吉非替尼
	var_list = info[0]
	knb = info[1]
	tumor_list = info[2]
	# 判断是否检出MET 扩增
	judge_met_cnv = "F"
	for var in var_list:
		if var["bio_category"] == "Cnv" and var["gene_symbol"] == "MET":
			judge_met_cnv = "T"
			break
	# 判断是否检出EGFR T790M
	judge_egfr_t790m = "F"
	for var in var_list:
		if var["bio_category"] == "Snvindel" and var["gene_symbol"] == "EGFR" and var["hgvs_p"] == "p.(T790M)":
			judge_egfr_t790m = "T"
			break
	#print ("judge_egfr_t790m", judge_egfr_t790m)
	result = []
	# 肺癌，若检出MET扩增，则所有EGFR 药物不展示
	# 肺癌，若检出EGFR T790M，则部分EGFR 变异的部分药物不展示
	for var in var_list:
		if "Predictive" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Predictive"]:
			regimen_list = []
			if "肺癌" in tumor_list:
				if var["bio_category"] == "Snvindel" and var["gene_symbol"] == "EGFR":
					judge_egfr_sense_var = "F"
					if var["hgvs_p"] in ["p.(L858R)", "p.(L861Q)", "p.(S768I)"]:
						judge_egfr_sense_var = "T"
					elif "G719" in var["hgvs_p"]:
						judge_egfr_sense_var = "T"
					# 19del规则暂定如下
					elif var["gene_region"] == "exon19" and "del" in var["hgvs_p"]:
						judge_egfr_sense_var = "T"
					#print ("judge_egfr_sense_var", judge_egfr_sense_var)
					
					if judge_met_cnv == "T":
						continue
					elif judge_egfr_t790m == "T" and judge_egfr_sense_var == "T":
						regimen_list = jdy_evi_sum(var["evi_sum"]["evi_split"]["Predictive"], ["阿法替尼", "达可替尼", "厄洛替尼", "吉非替尼"])
						#print (regimen_list)
					else:
						regimen_list = jdy_evi_sum(var["evi_sum"]["evi_split"]["Predictive"], [])
				else:
					regimen_list = jdy_evi_sum(var["evi_sum"]["evi_split"]["Predictive"], [])
			else:
				regimen_list = jdy_evi_sum(var["evi_sum"]["evi_split"]["Predictive"], [])
	
			if regimen_list:
				var["regimen_sum_for_jdyy"] = "，".join(regimen_list)
				result.append(var)
	result = sorted(result, key = lambda i:i["gene_symbol"])
	if knb:
		knb["regimen_sum_for_jdyy"] = "帕尼单抗，西妥昔单抗，西妥昔单抗β"
		result.insert(0, knb)

	return result
jinja2.filters.FILTERS["jlyy_10gene_result_sense_result_v3"] = jlyy_10gene_result_sense_result_v3

# 2025.05.12-华西乳腺癌模板-获取临床意义
# 2025.05.21-EPCMA风险描述要区分变异分组，药物和预后没有先不管了-考虑直接加在模板里
def hx_breast_template(info):
	var_list = info[0]
	tumor_list = info[1]
	hx_breast_template_significance = info[2]
	for var in var_list:
		num = 1
		var_significance = hx_breast_template_significance.get(var["gene_symbol"], {})
		# 1. 预后信息-仅乳腺癌
		prognostic = var_significance["prognostic_breast"] if "乳腺癌" in tumor_list and \
															  "prognostic_breast" in var_significance and \
															  var_significance["prognostic_breast"]  and \
															  var_significance["prognostic_breast"] != "-" \
														   else ""
		prognostic_num = ""
		if prognostic:
			prognostic_num = str(num)
			num += 1
		# 2. 用药-区分乳腺癌/卵巢癌
		predictive = ""
		if "乳腺癌" in tumor_list:
			predictive = var_significance["predictive_breast"] if "predictive_breast" in var_significance and \
																  var_significance["predictive_breast"] and \
																  var_significance["predictive_breast"] != "-" \
															   else ""
		elif "卵巢癌" in tumor_list:
			predictive = var_significance["predictive_ovarian"] if "predictive_ovarian" in var_significance and \
																   var_significance["predictive_ovarian"] and \
																   var_significance["predictive_ovarian"] != "-" \
																else ""
		else:
			predictive = "nofound"
				
		predictive_num = ""
		if predictive:
			predictive_num = str(num)
			num += 1
		# 3. 风险提示
		risk_suggest = var_significance["risk_inter"] if "risk_inter" in var_significance and \
														  var_significance["risk_inter"] \
													  else ""
		# 2025.05.21-新增EPCAM只适用Deletion分组
		#var_category_names = var["var_category_names"] if "var_category_names" in var.keys() and var["var_category_names"] else ""
		#if risk_suggest:
		#	if var["gene_symbol"] == "EPCAM" and "EPCAM Deletion (include exon8-9)" not in var_category_names:
		#		risk_suggest = "EPCAM基因的双等位基因致病性或疑似致病性胚系变异（纯合或复合杂合）与先天性腹泻5型伴簇绒状肠病（常染色体隐性遗传）相关（ClinGen）。\
		#						先天性腹泻5型伴簇绒状肠病是一种罕见的婴儿期慢性腹泻疾病，通常表现为在出生头几个月出现严重的慢性水样腹泻和生长受限，诊断后需要进行全\
		#						肠外营养治疗，有效绕过肠道保证热量和液体摄入，预防死亡（PMID: 18572020）。"
		# 2025.05.21-新增完成
		risk_suggest_num = "" 
		if risk_suggest:
			risk_suggest_num = num
			num += 1
		# 4. 总结：该变异与乳腺癌/卵巢癌的发病风险、预后、药物敏感性有关。
		tumor_risk = var_significance["tumor_risk"] if "tumor_risk" in var_significance and \
													   var_significance["tumor_risk"] and \
													   var_significance["tumor_risk"] != "-" \
													else ""
		summary_list = []
		summary_inter = ""
		if tumor_risk != "-":
			summary_list.append(tumor_risk+"的发病风险")
			if prognostic:
				summary_list.append("预后")
			if predictive and predictive != "nofound":
				summary_list.append("药物敏感性")
			summary_inter = "该变异与{0}有关。".format("、".join(summary_list))
		else:
			summary_inter = "该变异外显率低，与乳腺癌/卵巢癌的关系有限。"

		result = {
			"summary_inter" : summary_inter,
			"risk_suggest" : risk_suggest,
			"risk_suggest_num" : risk_suggest_num,
			"prognostic" : prognostic,
			"prognostic_num" : prognostic_num,
			"predictive" : predictive,
			"predictive_num" : predictive_num,
			"refer" : var_significance.get("refer", "")
		}
		var["hx_breast_result"] = result
	return var_list
jinja2.filters.FILTERS["hx_breast_template"] = hx_breast_template

# 重庆附一CP40，提示用药、预后、诊断信息-2025.05.13
# 该变异可能与（药物疗效、XX癌预后、XX癌辅助诊断相关）。
def cqfy_cp40_evi_sum(var):
	evi_type = []
	if "Predictive" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Predictive"]:
		evi_type.append("药物疗效")

	if "Prognostic" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Prognostic"]:
		tumor_lst = []
		for evi in var["evi_sum"]["evi_split"]["Prognostic"]:
			if "tumor_name_cn" in evi.keys() and evi["tumor_name_cn"] not in tumor_lst:
				tumor_lst.append(evi["tumor_name_cn"])
		evi_type.append("、".join(tumor_lst)+"预后")
	
	if "Diagnostic" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Diagnostic"]:
		tumor_lst = []
		for evi in var["evi_sum"]["evi_split"]["Diagnostic"]:
			if "tumor_name_cn" in evi.keys() and evi["tumor_name_cn"] not in tumor_lst:
				tumor_lst.append(evi["tumor_name_cn"])
		evi_type.append("、".join(tumor_lst)+"辅助诊断")

	return "该变异可能与{0}相关。".format("、".join(evi_type))
jinja2.filters.FILTERS["cqfy_cp40_evi_sum"] = cqfy_cp40_evi_sum

# 吉大一116新增甲状腺癌结果汇总-2025.05.15
def jdy_116_get_tc_summary(info):
	gene = info[0]
	var_list = info[1]
	# 2025.05.27-RET不限制在SV里
	#if gene in ["RET", "NTRK1", "NTRK3"]:
	if gene in ["NTRK1", "NTRK3"]:
		return [var for var in var_list if gene in re.split(",",var["gene_symbol"]) and var["bio_category"] == "Sv"]
	else:
		return [var for var in var_list if gene in re.split(",",var["gene_symbol"])]
jinja2.filters.FILTERS["jdy_116_get_tc_summary"] = jdy_116_get_tc_summary


# ---------------- 中山六院-添加参考文献序号-2025.05.16 ------------------开始
# MSI药物展示-2025.04.29
# 增加refer_code
def zsly_msi_regimen_v2(regimen_list):
	result = []
	regimen_sort = []
	regimen_dict = {}
	for regimen in regimen_list:
		if regimen["evi_conclusion_simple"] == "A":
			if regimen["regimen_name"] not in regimen_sort:
				regimen_sort.append(regimen["regimen_name"])
			if regimen["regimen_name"] not in regimen_dict.keys():
				regimen_dict.setdefault(regimen["regimen_name"], [])
			# 2025.06.03-加个兼容，原来没有参考文献是个列表，里面是空值，现在改为null
			regimen["literature_evi_sum"] = regimen["literature_evi_sum"] if regimen["literature_evi_sum"] else []
			regimen_dict[regimen["regimen_name"]].extend(regimen["literature_evi_sum"])
	for k, v in regimen_dict.items():
		result.append({
			"regimen_name" : k,
			"refer_code" : v
		})
	return sorted(result, key = lambda i:regimen_sort.index(i["regimen_name"]))
jinja2.filters.FILTERS["zsly_msi_regimen_v2"] = zsly_msi_regimen_v2

# 获取io所有参考文献-线下配置的，无code
def zsly_get_io_refer_all(tumor_inter, zsly_result, tumor_list):
	io_p_var = zsly_io_inter([tumor_inter, zsly_result, tumor_list, "p"])
	io_n_var = zsly_io_inter([tumor_inter, zsly_result, tumor_list, "n"])
	io_refer_list = []
	for var_item in io_p_var + io_n_var:
		for i in re.split(", ", var_item["refer"]):
			j = re.split(": ", i)
			if j[0] == "PMID":
				if {"type" : "pmid", "pmid" : j[1], "title" : ""} not in io_refer_list:
					io_refer_list.append({"type" : "pmid", "pmid" : j[1], "title" : ""})
			else:
				if {"type" : "", "pmid" : "", "title" : "{0}: {1}".format(j[0], j[1])} not in io_refer_list:
					io_refer_list.append({"type" : "", "pmid" : "", "title" : "{0}: {1}".format(j[0], j[1])})
	return io_refer_list

# 中山六院参考文献-2025.04.29
def zsly_all_refer(var, msi, tumor_list, refer_original):
	# 汇总治疗方案顺序
	result = []
	num = 1
	# code并排序
	refer_sort = []
	regimen_list = []
	# KNB 本癌种+其他癌种
	if var["knb"]:
		regimen_list.extend(zsly_evi_redup(zsly_cp200_filter_intumor_evi(var["knb"]["evi_sum"])))
		regimen_list.extend(zsly_evi_redup(zsly_cp200_filter_outtumor_evi(["I", var["knb"]["evi_sum"]])))
	# I类单突变 本癌种+其他癌种
	var_I_list = zsly_cp200_sort_var(var["var_somatic"]["level_I"])
	for i in var_I_list:
		regimen_list.extend(zsly_evi_redup(zsly_cp200_filter_intumor_evi(i["evi_sum"])))
		regimen_list.extend(zsly_evi_redup(zsly_cp200_filter_outtumor_evi(["I", i["evi_sum"]])))
	# I类共突变 本癌种+其他癌种
	var_I_co_list = zsly_co_mutation([var["var_somatic"]["level_I"] + var["var_somatic"]["level_II"], var["knb"], msi, "I"])
	for i in var_I_co_list:
		regimen_list.extend(zsly_evi_redup(zsly_cp200_co_get_intumor_evi(i["evi"])))
		regimen_list.extend(zsly_evi_redup(zsly_co_get_outtumor_evi(["I", i["evi"]])))
	# II类单突变 本癌种+其他癌种
	var_II_list = zsly_cp200_sort_var(var["var_somatic"]["level_II"])
	for i in var_II_list:
		regimen_list.extend(zsly_evi_redup(zsly_cp200_filter_intumor_evi(i["evi_sum"])))
		regimen_list.extend(zsly_evi_redup(zsly_cp200_filter_outtumor_evi(["II", i["evi_sum"]])))
	# II类共突变 本癌种+其他癌种
	var_II_co_list = zsly_co_mutation([var["var_somatic"]["level_I"] + var["var_somatic"]["level_II"], var["knb"], msi, "II"])
	for i in var_II_co_list:
		regimen_list.extend(zsly_evi_redup(zsly_cp200_co_get_intumor_evi(i["evi"])))
		regimen_list.extend(zsly_evi_redup(zsly_co_get_outtumor_evi(["II", i["evi"]])))
	
	for evi in regimen_list:
		for i in evi["evi_refer_type"]:
			#print (i)
			for j in i["refer_code"]:
				if j not in refer_sort:
					refer_sort.append(j)

	# MSI
	if msi["var_id"] == "MSI-H":
		msi_regimen_list = zsly_msi_regimen_v2(msi["evi_sum"]["evi_split"]["Predictive"])
		for regimen in msi_regimen_list:
			for j in regimen["refer_code"]:
				if j not in refer_sort:
					refer_sort.append(j)
	for refer_code in refer_sort:
		result.append(
			{
				"num" : num,
				"info" : zsly_get_refer_code_info([refer_code, refer_original])
			}
		)
		num += 1

	pmid_list = [i["info"]["pmid"] for i in result if i["info"]["pmid"]]
	other_list = [i["info"]["title"] for i in result if not i["info"]["pmid"]]

	# IO
	io_refer_list = zsly_get_io_refer_all(var["io"]["tumor_inter"], var["io"]["zsly_result"], tumor_list)
	io_refer_need_add_to_result = []
	for i in io_refer_list:
		if (i["pmid"] and i["pmid"] not in pmid_list) or (not i["pmid"] and i["title"] not in other_list):
			io_refer_need_add_to_result.append(i)
	for i in io_refer_need_add_to_result:
		result.append(
			{
				"num" : num,
				"info" : i,
				"type" : ""
			}
		)
		num += 1

	# 对result进行转化
	result_stran = {}
	for i in result:
		#code = i["info"]["refer_code"] if "refer_code" in i["info"].keys() and i["info"]["refer_code"] else \
		#	   i["info"]["pmid"] if "pmid" in i["info"].keys() and i["info"]["pmid"] else \
		#	   i["info"]["title"] if "title" in i["info"].keys() else ""
		#result_stran[code] = i
		# 这部分改进一下，key包含code、pmid和title
		code = i["info"]["refer_code"] if "refer_code" in i["info"].keys() and i["info"]["refer_code"] else ""
		if code:
			result_stran[code] = i
		pmid = i["info"]["pmid"] if "pmid" in i["info"].keys() and i["info"]["pmid"] else ""
		if pmid:
			result_stran[pmid] = i
		title = i["info"]["title"] if "title" in i["info"].keys() else ""
		if title:
			result_stran[title] = i

	return result, result_stran


# 获取参考文献code对应数据-2025.04.29
def zsly_get_refer_code_info(info):
	code_info = {}
	code = info[0]
	refer_original = info[1]
	for i in refer_original:
		if "refer_code" in i.keys() and i["refer_code"] and i["refer_code"] == code:
			i["pmid"] = (re.split("PMID:", i["pmid"])[-1]).replace("]", "") if i["pmid"] else ""
			code_info = i
			break
	return code_info

# 输入参考文献，输出序号-模板调用的
def zsly_get_refer_num(info):
	var = info[0]
	msi = info[1]
	tumor_list = info[2]
	refer_code = info[3]
	refer_original = info[4]
	refer_sort_list, refer_sort_dict = zsly_all_refer(var, msi, tumor_list, refer_original)
	final_num = []
	#print ("参考文献排序结果",refer_sort_list)
	#print (refer_sort_dict)
	#print ("输入参考文献code",refer_code)
	for i in refer_code:
		if i in refer_sort_dict.keys():
			#print ("yes")
			index_num = refer_sort_dict[i]["num"]
		#index_num = refer_sort_list.index(i)+1 if i in refer_sort_list else ""
			if index_num and index_num not in final_num:
				final_num.append(index_num)

	final_num = sorted(final_num)
	final_num_2 = [str(i) for i in final_num]
	#print ("输出参考文献序号",final_num_2)
	#return final_num
	return final_num_2
	#return ",".join(final_num_2)
jinja2.filters.FILTERS["zsly_get_refer_num"] = zsly_get_refer_num

# 处理IO参考文献
# 格式为“PMID: xxx, PMID: xxx, DOI: xxx”
# 输出[xxx, xxx, xxx]
def zsly_io_refer_to_list(info):
	refer_list = []
	for i in re.split(", ", info):
			j = re.split(": ", i)
			if j[0] == "PMID":
				if j[1] not in refer_list:
					refer_list.append(j[1])
			else:
				if "{0}: {1}".format(j[0], j[1]) not in refer_list:
					refer_list.append("{0}: {1}".format(j[0], j[1]))
	return refer_list
jinja2.filters.FILTERS["zsly_io_refer_to_list"] = zsly_io_refer_to_list

# 参考文献汇总表格有/无PMID-2025.04.29
def zsly_refer_sum(info):
	var = info[0]
	msi = info[1]
	tumor_list = info[2]
	refer_original = info[3]
	refer_sort_list, refer_sort_dict = zsly_all_refer(var, msi, tumor_list, refer_original)
	return zsly_pmid_v2(refer_sort_list)
jinja2.filters.FILTERS["zsly_refer_sum"] = zsly_refer_sum

# 参考文献汇总表-无PMID-2025.04.29
def zsly_refer_without_pmid(info):
	var = info[0]
	msi = info[1]
	tumor_list = info[2]
	refer_original = info[3]
	refer_sort_list, refer_sort_dict = zsly_all_refer(var, msi, tumor_list, refer_original)
	return [i for i in refer_sort_list if not i["info"]["pmid"]]
jinja2.filters.FILTERS["zsly_refer_without_pmid"] = zsly_refer_without_pmid

def zsly_pmid_v2(refer):
	pmid_list = []
	for i in refer:
		pmid = (re.split("PMID:", i["info"]["pmid"])[-1]).replace("]", "") if i["info"]["pmid"] else ""
		if (i["num"], pmid) not in pmid_list:
			pmid_list.append((i["num"], pmid))

	pmid_all_index = []
	num = 1
	for i in pmid_list:
		pmid_all_index.append({"num" : num, "pmid" : i})
		num += 1
	
	if pmid_all_index:
		rest_num = len(pmid_all_index) % 4
		key = int(len(pmid_all_index) / 4)
		list1 = []
		list2 = []
		list3 = []
		list4 = []
		if rest_num == 1:
			list1 = pmid_all_index[0 : key+1]
			list2 = pmid_all_index[key+1 : key*2+1]
			list2.append({"num" : "", "pmid" : ""})
			list3 = pmid_all_index[key*2+1 : key*3+1]
			list3.append({"num" : "", "pmid" : ""})
			list4 = pmid_all_index[key*3+1 :]
			list4.append({"num" : "", "pmid" : ""})
		elif rest_num == 2:
			list1 = pmid_all_index[0 : key+1]
			list2 = pmid_all_index[key+1 : key*2+2]
			list3 = pmid_all_index[key*2+2 : key*3+2]
			list3.append({"num" : "", "pmid" : ""})
			list4 = pmid_all_index[key*3+2 :]
			list4.append({"num" : "", "pmid" : ""})
		elif rest_num == 3:
			list1 = pmid_all_index[0 : key+1]
			list2 = pmid_all_index[key+1 : key*2+2]
			list3 = pmid_all_index[key*2+2 : key*3+3]
			list4 = pmid_all_index[key*3+3 :]
			list4.append({"num" : "", "pmid" : ""})
		else:
			list1 = pmid_all_index[0 : key]
			list2 = pmid_all_index[key : key*2]
			list3 = pmid_all_index[key*2 : key*3]
			list4 = pmid_all_index[key*3 :]

	result = []
	if pmid_all_index:
		for i in range(0, len(list1)):
			result.append({
				"num_1" : list1[i]["num"] if list1[i]["num"] else "",
				"pmid_1" : list1[i]["pmid"] if list1[i]["pmid"] else "",
				"num_2" : list2[i]["num"] if list2[i]["num"] else "",
				"pmid_2" : list2[i]["pmid"] if list2[i]["pmid"] else "",
				"num_3" : list3[i]["num"] if list3[i]["num"] else "",
				"pmid_3" : list3[i]["pmid"] if list3[i]["pmid"] else "",
				"num_4" : list4[i]["num"] if list4[i]["num"] else "",
				"pmid_4" : list4[i]["pmid"] if list4[i]["pmid"] else ""
			})

	return result
	
# 中山六院-参考文献序号相邻的合并-2025.05.23
def zsly_merge_num(num_list):
	sort_list = sorted(list(set([int(n) for n in num_list])))
	len_list = len(sort_list)
	i = 0
	split_list = []
	tmp_list = [sort_list[i]]
	while True:
		if i + 1 == len_list:
			break
		next_n = sort_list[i+1]
		if sort_list[i] + 1 == next_n:
			tmp_list.append(next_n)
		else:
			split_list.append(tmp_list)
			tmp_list = [next_n]
		i += 1
	split_list.append(tmp_list)
		
	refer_num = []
	for i in split_list:
		if len(i) == 1:
			refer_num.append(str(i[0]))
		else:
			refer_num.append(str(i[0]) + "-" + str(i[-1]))
	refer_num_str = ", ".join(refer_num)

	return refer_num_str
jinja2.filters.FILTERS["zsly_merge_num"] = zsly_merge_num

# 化疗判断所有知识库有记录位点频率-2025.05.23
# 与上一版相比区别在：loss2归到杂合型/纯合型灰区，上一版是在野生型/杂合型灰区
def zsly_chemo_summary_v3(info):
	chemo_list = info[0]
	all_detect_result = info[1]
	result = []
	for var in chemo_list:
		if zsly_chemo_var_v3([var, all_detect_result]) and zsly_chemo_var_v3([var, all_detect_result]) in ["loss1"]:
			if "loss" not in result:
				result.append("loss")
		if zsly_chemo_var_v3([var, all_detect_result]) and zsly_chemo_var_v3([var, all_detect_result]) in ["loss2", "more"]:
			if "more" not in result:
				result.append("more")
	return result
jinja2.filters.FILTERS["zsly_chemo_summary_v3"] = zsly_chemo_summary_v3

# 化疗判断所有知识库有记录位点频率-2025.06.11
# 与上一版相比区别在：loss2归到杂合型/纯合型灰区，上一版是在野生型/杂合型灰区
# 2026.06.11-新增具体灰区
def zsly_chemo_summary_v4(info):
	chemo_list = info[0]
	all_detect_result = info[1]
	result = []
	for var in chemo_list:
		if zsly_chemo_var_v3([var, all_detect_result]) and zsly_chemo_var_v3([var, all_detect_result]) in ["loss1"]:
			if "loss1" not in result:
				result.append("loss1")
		if zsly_chemo_var_v3([var, all_detect_result]) and zsly_chemo_var_v3([var, all_detect_result]) in ["loss2"]:
			if "loss2" not in result:
				result.append("loss2")
		if zsly_chemo_var_v3([var, all_detect_result]) and zsly_chemo_var_v3([var, all_detect_result]) in ["more"]:
			if "more" not in result:
				result.append("more")
	return result
jinja2.filters.FILTERS["zsly_chemo_summary_v4"] = zsly_chemo_summary_v4

# ---------------- 中山六院-添加参考文献序号-2025.05.16 ------------------完成

# 北大人民MP-胃肠道间质瘤结果小结需要额外展示KIT、PDGFRA、SDH（包含SDHA/B/C/D）基因检测结果-2025.05.26
def dbrm_ga_filter(info):
	var_list = info[0]
	gene_list = info[1]
	var_origin = info[2]
	gene_var_list = [var for var in var_list if set(re.split(",", var["gene_symbol"])) & set(gene_list)]
	result = []
	for var in gene_var_list:
		if var_origin == "somatic":
			if var["bio_category"] == "Snvindel":
				if set(["SDHA", "SDHB", "SDHC", "SDHD"]) & set(gene_list):
					result.append("{0} {1}（{2}类）".format(var["gene_symbol"], \
										  				    var["hgvs_p_BJYY"] if var["hgvs_p_BJYY"] != "p.?" else var["hgvs_c"], \
										   					"I" if var["clinic_num_s"] == 5 else "II"))
				else:
					result.append("{0}（{1}类）".format(var["hgvs_p_BJYY"] if var["hgvs_p_BJYY"] != "p.?" else var["hgvs_c"], \
										   				"I" if var["clinic_num_s"] == 5 else "II"))
			elif var["bio_category"] == "Cnv":
				result.append("扩增（{0}类）".format("I" if var["clinic_num_s"] == 5 else "II"))
			# 2025.07.21-新增HD
			elif var["bio_category"] == "PHd":
				hd_type = var["gene_symbol"] + " 未知变异类型！"
				if var["type"] == "HomoDel":
					hd_type = var["gene_symbol"] + "纯合缺失"
				elif var["type"] == "HeteDel":
					hd_type = var["gene_symbol"] + "杂合缺失"
				if hd_type not in result:
					result.append(hd_type)
			# 2025.07.21-新增完成
			elif var["bio_category"] in ["Sv", "PSeqRnaSv"]:
				result.append("{0}:{1}-{2}:{3}融合（{4}类）".format(var["five_prime_gene"], var["five_prime_cds"], \
												   				   var["three_prime_gene"], var["three_prime_cds"], \
									   							   "I" if var["clinic_num_s"] == 5 else "II"))
		else:
			if set(["SDHA", "SDHB", "SDHC", "SDHD"]) & set(gene_list):
				result.append("{0} {1}（{2}）".format(var["gene_symbol"], \
													  var["hgvs_p_BJYY"] if var["hgvs_p_BJYY"] != "p.?" else var["hgvs_c"], \
										   			  "致病性变异" if var["clinic_num_g"] == 5 else "疑似致病性变异"))
			else:
				result.append("{0}（{1}）".format(var["hgvs_p_BJYY"] if var["hgvs_p_BJYY"] != "p.?" else var["hgvs_c"], \
										   		  "致病性变异" if var["clinic_num_g"] == 5 else "疑似致病性变异"))
	return "、".join(result)
jinja2.filters.FILTERS["dbrm_ga_filter"] = dbrm_ga_filter

# 西安交大一-用药后面加个“等”-2025.05.26
def xajdy_regimen_add_deng(regimen_list):
	result = []
	for regimen in regimen_list:
		result.append("{0}({1}级)".format(regimen["regimen_name"], regimen["evi_conclusion_simple"]))
	return "\t".join(result) + "等"
jinja2.filters.FILTERS["xajdy_regimen_add_deng"] = xajdy_regimen_add_deng


# 西部战区总医院-116-展示10个基因检测结果-2025.05.28
def xbzqzy_116_lc10(info):
	gene = info[0]
	var_list = info[1]
	gene_var_list = [var for var in var_list if gene in re.split(",", var["gene_symbol"])]
	return gene_var_list
jinja2.filters.FILTERS["xbzqzy_116_lc10"] = xbzqzy_116_lc10

# 西安交大一结果解读部分增加临床证据汇总描述（仅用药和肺癌KRAS A级预后）-2025.06.03
# 1. 存在AB则展示AB，否则展示CD
# NCCN A - 敏感和耐药：展示描述
# CSCO A - 敏感和耐药：展示描述
# 其他 A - 敏感：      XXX、XXX已获批用于携带该类变异的患者
# NCCN C3 - 敏感和耐药：展示描述
# CSCO C3 - 敏感和耐药：展示描述
# 其他 C3 - 敏感：      XXX、XXX已获批用于携带该类变异的XXX、XXX（证据癌种）患者
# 其他 B/C - 敏感：临床试验表明，XXX、XXX对携带该类变异的XXX、XXX（证据癌种）患者可能敏感
# 其他 B/C - 耐药：临床试验表明，携带该类变异的XXX（证据癌种）患者可能对XXX耐药
# D - 敏感：有研究表明，该类变异可能对XXX、XXX敏感
# D - 耐药：有研究表明，该类变异可能对XXX、XXX耐药
def xajdy_evi_summary_old(evi_list):
	# A级证据处理
	AB_evi = []
	A_appr_regimen = []
	for evi in evi_list:
		if evi["evi_conclusion_simple"] == "A":
			if "NCCN" in evi["refer_agency"]:
				if evi["evi_interpretation"] not in AB_evi:
					AB_evi.append(evi["evi_interpretation"])
			elif "CSCO" in evi["refer_agency"]:
				if evi["evi_interpretation"] not in AB_evi:
					AB_evi.append(evi["evi_interpretation"])
			else:
				if evi["regimen_name"] not in A_appr_regimen:
					A_appr_regimen.append(evi["regimen_name"])
	if A_appr_regimen:
		AB_evi.append("{0}已获批用于携带该类变异的患者。".format("、".join(A_appr_regimen)))	
	# B级证据处理
	B_snese_tumor_and_regimen = {}
	B_resis_tumor_and_regimen = {}
	for evi in evi_list:
		if evi["evi_conclusion_simple"] == "B":
			tumor_name = evi["tumor_name_cn"] if evi["tumor_name_cn"] else evi["tumor_name_en"] if "tumor_name_en" in evi.keys() and evi["tumor_name_en"] else ""
			if evi["clinical_significance_cn"] == "敏感":
				if tumor_name not in B_snese_tumor_and_regimen.keys():
					B_snese_tumor_and_regimen.setdefault(tumor_name, [])
				B_snese_tumor_and_regimen[tumor_name].append(evi["regimen_name"])
			else:
				if tumor_name not in B_resis_tumor_and_regimen.keys():
					B_resis_tumor_and_regimen.setdefault(tumor_name, [])
				B_resis_tumor_and_regimen[tumor_name].append(evi["regimen_name"])
	if B_snese_tumor_and_regimen:
		for tumor, regimen_list in B_snese_tumor_and_regimen.items():
			AB_evi.append("临床试验表明，{0}对携带该类变异的{1}患者可能敏感。".format("、".join(regimen_list), tumor))
	if B_resis_tumor_and_regimen:
		for tumor, regimen_list in B_resis_tumor_and_regimen.items():
			AB_evi.append("临床试验表明，携带该类变异的{0}患者可能对{1}耐药。".format(tumor, "、".join(regimen_list)))
	# C3级处理
	CD_evi = []
	C_appr_regimen = {}
	for evi in evi_list:
		if evi["evi_conclusion"] == "C3":
			if "NCCN" in evi["refer_agency"]:
				if evi["evi_interpretation"] not in CD_evi:
					CD_evi.append(evi["evi_interpretation"])
			elif "CSCO" in evi["refer_agency"]:
				if evi["evi_interpretation"] not in CD_evi:
					CD_evi.append(evi["evi_interpretation"])
			else:
				tumor_name = evi["tumor_name_cn"] if evi["tumor_name_cn"] else evi["tumor_name_en"] if "tumor_name_en" in evi.keys() and evi["tumor_name_en"] else ""
				if tumor_name not in C_appr_regimen.keys():
					C_appr_regimen.setdefault(tumor_name, [])
				C_appr_regimen[tumor_name].append(evi["regimen_name"])
	for tumor, regimen in C_appr_regimen.items():
		CD_evi.append("{0}已获批用于携带该类变异的{1}患者。".format("、".join(regimen), tumor))
	# 其他C级处理
	C_snese_tumor_and_regimen = {}
	C_resis_tumor_and_regimen = {}
	for evi in evi_list:
		if evi["evi_conclusion_simple"] == "C" and evi["evi_conclusion"] != "C3":
			tumor_name = evi["tumor_name_cn"] if evi["tumor_name_cn"] else evi["tumor_name_en"] if "tumor_name_en" in evi.keys() and evi["tumor_name_en"] else ""
			if evi["clinical_significance_cn"] == "敏感":
				if tumor_name not in C_snese_tumor_and_regimen.keys():
					C_snese_tumor_and_regimen.setdefault(tumor_name, [])
				C_snese_tumor_and_regimen[tumor_name].append(evi["regimen_name"])
			else:
				if tumor_name not in C_resis_tumor_and_regimen.keys():
					C_resis_tumor_and_regimen.setdefault(tumor_name, [])
				C_resis_tumor_and_regimen[tumor_name].append(evi["regimen_name"])
	if C_snese_tumor_and_regimen:
		for tumor, regimen_list in C_snese_tumor_and_regimen.items():
			CD_evi.append("临床试验表明，{0}对携带该类变异的{1}患者可能敏感。".format("、".join(regimen_list), tumor))
	if C_resis_tumor_and_regimen:
		for tumor, regimen_list in C_resis_tumor_and_regimen.items():
			CD_evi.append("临床试验表明，携带该类变异的{0}患者可能对{1}耐药。".format(tumor, "、".join(regimen_list)))
	# D级处理
	D_sense_regimen = [evi["regimen_name"] for evi in evi_list if evi["evi_conclusion_simple"] == "D" and evi["clinical_significance_cn"] == "敏感"]
	D_resis_regimen = [evi["regimen_name"] for evi in evi_list if evi["evi_conclusion_simple"] == "D" and evi["clinical_significance_cn"] == "耐药"]
	if D_sense_regimen:
		CD_evi.append("有研究表明，该类变异可能对{0}敏感。".format("、".join(D_sense_regimen)))
	if D_resis_regimen:
		CD_evi.append("该类变异可能对{0}耐药。".format("、".join(D_resis_regimen)))
	
	if AB_evi:
		return AB_evi
	else:
		return CD_evi
jinja2.filters.FILTERS["xajdy_evi_summary_old"] = xajdy_evi_summary_old

# ---2025.06.05-规则更新---------------------------------------------------------------
# 有A级就仅展示A级，没有则展示BCD
# A-敏感：针对携带该变异的患者，获批或指南推荐的药物有XXX、XXX等
# A-耐药：展示描述
# C3-敏感：针对该突变，其他癌种中已有获批或指南推荐的药物，如XXX药用于XXX癌，YYY药用于YYY癌等。
# B/其他C-敏感：临床试验表明，XXX对携带该类变异的XXX癌可能敏感、YYY对携带该类变异的YYY癌可能敏感。
# B/C-耐药：临床试验表明，XXX对携带该类变异的XXX癌可能耐药、YYY对携带该类变异的YYY癌可能耐药。
# D-敏感：有研究表明，给类变异可能对XXX、XXX敏感。
# D-耐药：该类变异可能对XXX、XXX耐药。
# -------------------------------------------------------------------------------------
def xajdy_evi_summary(evi_list):
	# A级证据处理
	A_evi = []
	A_sense_regimen = []
	A_resis_evi = []
	for evi in evi_list:
		if evi["evi_conclusion_simple"] == "A":
			if evi["clinical_significance_cn"] == "敏感":
				if evi["regimen_name"] not in A_sense_regimen:
					A_sense_regimen.append(evi["regimen_name"])
			else:
				if evi["evi_interpretation"] not in A_resis_evi:
					A_resis_evi.append(evi["evi_interpretation"])
	if A_sense_regimen:
		A_evi.append("针对携带该变异的患者，获批或指南推荐的药物有{0}等。".format("、".join(A_sense_regimen)))
	if A_resis_evi:
		A_evi.extend(A_resis_evi)
	
	# BCD级证据处理
	BCD_evi = []
	c3_sense_tumor_and_regimen = {}
	c3_resis_evi = []
	b_c_other_sense_tumor_and_regimen = {}
	b_c_resis_tumor_and_regimen = {}
	for evi in evi_list:
		tumor_name = evi["tumor_name_cn"] if "tumor_name_cn" in evi.keys() and evi["tumor_name_cn"] else evi["tumor_name_en"] if "tumor_name_en" in evi.keys() and evi["tumor_name_en"] else ""
		if evi["evi_conclusion"] == "C3" and evi["clinical_significance_cn"] == "敏感":
			if tumor_name not in c3_sense_tumor_and_regimen.keys():
				c3_sense_tumor_and_regimen.setdefault(tumor_name, [])
			c3_sense_tumor_and_regimen[tumor_name].append(evi["regimen_name"])
		elif evi["evi_conclusion"] == "C3" and evi["clinical_significance_cn"] == "耐药":
			if evi["evi_interpretation"] not in c3_resis_evi:
				c3_resis_evi.append(evi["evi_interpretation"])
		elif evi["evi_conclusion_simple"] in ["C", "B"] and evi["clinical_significance_cn"] == "敏感":
			if tumor_name not in b_c_other_sense_tumor_and_regimen.keys():
				b_c_other_sense_tumor_and_regimen.setdefault(tumor_name, [])
			b_c_other_sense_tumor_and_regimen[tumor_name].append(evi["regimen_name"])
		elif evi["evi_conclusion_simple"] in ["C", "B"] and evi["clinical_significance_cn"] == "耐药":
			if tumor_name not in b_c_resis_tumor_and_regimen.keys():
				b_c_resis_tumor_and_regimen.setdefault(tumor_name, [])
			b_c_resis_tumor_and_regimen[tumor_name].append(evi["regimen_name"])
	d_sense_regimen = [evi["regimen_name"] for evi in evi_list if evi["evi_conclusion_simple"] == "D" and evi["clinical_significance_cn"] == "敏感"]
	d_resis_regimen = [evi["regimen_name"] for evi in evi_list if evi["evi_conclusion_simple"] == "D" and evi["clinical_significance_cn"] == "耐药"]
	# 拼接
	if c3_sense_tumor_and_regimen:
		tmp_list = []
		for tumor, regimen in c3_sense_tumor_and_regimen.items():
			tmp_list.append("{0}用于{1}".format("、".join(regimen), tumor))
		BCD_evi.append("针对该突变，其他癌种中已有获批或指南推荐的药物，如{0}等。".format("、".join(tmp_list)))
	if c3_resis_evi:
		BCD_evi.extend(c3_resis_evi)
	b_c_tmp_list = []
	if b_c_other_sense_tumor_and_regimen:
		for tumor, regimen in b_c_other_sense_tumor_and_regimen.items():
			b_c_tmp_list.append("{0}对携带该类变异的{1}可能敏感".format("、".join(regimen), tumor))
	if b_c_resis_tumor_and_regimen:
		for tumor, regimen in b_c_resis_tumor_and_regimen.items():
			b_c_tmp_list.append("{0}对携带该类变异的{1}可能耐药".format("、".join(regimen), tumor))
	if b_c_tmp_list:
		BCD_evi.append("临床试验表明，{0}。".format("、".join(b_c_tmp_list)))
	d_tmp_list = []
	if d_sense_regimen:
		d_tmp_list.append("可能对{0}敏感".format("、".join(d_sense_regimen)))
	if d_resis_regimen:
		d_tmp_list.append("可能对{0}耐药".format("、".join(d_resis_regimen)))
	if d_tmp_list:
		BCD_evi.append("有研究表明，该类变异{0}。".format("，".join(d_tmp_list)))
	
	if A_evi:
		return "".join(A_evi)
	else:
		return "".join(BCD_evi)
jinja2.filters.FILTERS["xajdy_evi_summary"] = xajdy_evi_summary

# ---2025.06.11-规则更新---------------------------------------------------------------
# 有A级就仅展示A级，没有则展示BCD
# A-敏感：针对携带该变异的患者，获批或指南推荐的药物有XXX、XXX等
# A-耐药：展示描述
# C3-敏感：针对该突变，其他癌种中已有获批或指南推荐的药物，如XXX药用于XXX癌，YYY药用于YYY癌等。
# B/其他C-敏感：临床试验表明，XXX对携带该类变异的XXX癌可能敏感、YYY对携带该类变异的YYY癌可能敏感。
# B/C-耐药：临床试验表明，XXX对携带该类变异的XXX癌可能耐药、YYY对携带该类变异的YYY癌可能耐药。
# D-敏感：有研究表明，给类变异可能对XXX、XXX敏感。
# D-耐药：该类变异可能对XXX、XXX耐药。
# 2025.06.11-新增参考文献
# -------------------------------------------------------------------------------------
def get_pmid_from_inter(inter):
	pmid_list = []
	mat = re.compile(r"PMID.\s?\d+")
	for i in mat.findall(str(inter)):
		if re.search(":|: |：|： ", i):
			pmid = (re.split(":|: |：|： ", i))[1].replace(" ", "")
		else:
			pmid = (re.split("PMID", i))[1]
		pmid_list.append(pmid)
	return pmid_list

def xajdy_evi_summary_v2(evi_list):
	# A级证据处理
	A_evi = []
	A_sense_regimen = []
	A_resis_evi = []
	for evi in evi_list:
		if evi["evi_conclusion_simple"] == "A":
			if evi["clinical_significance_cn"] == "敏感":
				if evi["regimen_name"] not in A_sense_regimen:
					A_sense_regimen.append(evi["regimen_name"])
			else:
				if evi["evi_interpretation"] not in A_resis_evi:
					A_resis_evi.append(evi["evi_interpretation"])				
	if A_sense_regimen:
		A_evi.append("针对携带该变异的患者，获批或指南推荐的药物有{0}等。".format("、".join(A_sense_regimen)))
	if A_resis_evi:
		A_evi.extend(A_resis_evi)
	
	# BCD级证据处理
	BCD_evi = []
	c3_sense_tumor_and_regimen = {}
	c3_resis_evi = []
	b_c_other_sense_tumor_and_regimen = {}
	b_c_resis_tumor_and_regimen = {}
	for evi in evi_list:
		pmid_list = get_pmid_from_inter(evi["evi_interpretation"])
		tumor_name = evi["tumor_name_cn"] if "tumor_name_cn" in evi.keys() and evi["tumor_name_cn"] else evi["tumor_name_en"] if "tumor_name_en" in evi.keys() and evi["tumor_name_en"] else ""
		# C3无需PMID
		if evi["evi_conclusion"] == "C3" and evi["clinical_significance_cn"] == "敏感":
			if tumor_name not in c3_sense_tumor_and_regimen.keys():
				c3_sense_tumor_and_regimen.setdefault(tumor_name, [])
			c3_sense_tumor_and_regimen[tumor_name].append(evi["regimen_name"])
		elif evi["evi_conclusion"] == "C3" and evi["clinical_significance_cn"] == "耐药":
			if evi["evi_interpretation"] not in c3_resis_evi:
				c3_resis_evi.append(evi["evi_interpretation"])
		# 其他B/C/D需要PMID
		elif evi["evi_conclusion_simple"] in ["C", "B"] and evi["clinical_significance_cn"] == "敏感":
			if tumor_name not in b_c_other_sense_tumor_and_regimen.keys():
				b_c_other_sense_tumor_and_regimen.setdefault(tumor_name, [])
			b_c_other_sense_tumor_and_regimen[tumor_name].append((evi["regimen_name"], pmid_list))
		elif evi["evi_conclusion_simple"] in ["C", "B"] and evi["clinical_significance_cn"] == "耐药":
			if tumor_name not in b_c_resis_tumor_and_regimen.keys():
				b_c_resis_tumor_and_regimen.setdefault(tumor_name, [])
			b_c_resis_tumor_and_regimen[tumor_name].append((evi["regimen_name"], pmid_list))
	d_sense_regimen = [(evi["regimen_name"], get_pmid_from_inter(evi["evi_interpretation"])) for evi in evi_list if evi["evi_conclusion_simple"] == "D" and evi["clinical_significance_cn"] == "敏感"]
	d_resis_regimen = [(evi["regimen_name"], get_pmid_from_inter(evi["evi_interpretation"])) for evi in evi_list if evi["evi_conclusion_simple"] == "D" and evi["clinical_significance_cn"] == "耐药"]
	# 拼接
	if c3_sense_tumor_and_regimen:
		tmp_list = []
		for tumor, regimen in c3_sense_tumor_and_regimen.items():
			tmp_list.append("{0}用于{1}".format("、".join(regimen), tumor))
		BCD_evi.append("针对该突变，其他癌种中已有获批或指南推荐的药物，如{0}等。".format("，".join(tmp_list)))
	if c3_resis_evi:
		BCD_evi.extend(c3_resis_evi)

	b_c_tmp_list = []
	if b_c_other_sense_tumor_and_regimen:
		for tumor, regimen in b_c_other_sense_tumor_and_regimen.items():
			regimen_str = "、".join([i[0] for i in regimen])
			pmid_list = []
			for item in regimen:
				for pmid in item[1]:
					if pmid not in pmid_list:
						pmid_list.append(pmid)
			if pmid_list:
				b_c_tmp_list.append("{0}在携带该类变异的{1}中可能敏感（{2}）".format(regimen_str, tumor, ", ".join(["PMID:"+str(i) for i in pmid_list])))
			else:
				b_c_tmp_list.append("{0}在携带该类变异的{1}中可能敏感".format(regimen_str, tumor))
	if b_c_resis_tumor_and_regimen:
		for tumor, regimen in b_c_resis_tumor_and_regimen.items():
			regimen_str = "、".join([i[0] for i in regimen])
			pmid_list = []
			for item in regimen:
				for pmid in item[1]:
					if pmid not in pmid_list:
						pmid_list.append(pmid)
			if pmid_list:
				b_c_tmp_list.append("{0}在携带该类变异的{1}中可能耐药（{2}）".format(regimen_str, tumor, ", ".join(["PMID:"+str(i) for i in pmid_list])))
			else:
				b_c_tmp_list.append("{0}在携带该类变异的{1}中可能耐药".format(regimen_str, tumor))
	if b_c_tmp_list:
		BCD_evi.append("临床试验表明，{0}。".format("，".join(b_c_tmp_list)))
		
	d_tmp_list = []
	if d_sense_regimen:
		regimen_str = "、".join([i[0] for i in d_sense_regimen])
		pmid_list = []
		for item in d_sense_regimen:
			for pmid in item[1]:
				if pmid not in pmid_list:
					pmid_list.append(pmid)
		if pmid_list:
			d_tmp_list.append("可能对{0}敏感（{1}）".format(regimen_str, ", ".join(["PMID:"+str(i) for i in pmid_list])))
		else:
			d_tmp_list.append("可能对{0}敏感".format(regimen_str))
	if d_resis_regimen:
		regimen_str = "、".join([i[0] for i in d_resis_regimen])
		pmid_list = []
		for item in d_resis_regimen:
			for pmid in item[1]:
				if pmid not in pmid_list:
					pmid_list.append(pmid)
		if pmid_list:
			d_tmp_list.append("可能对{0}耐药（{1}）".format(regimen_str, ", ".join(["PMID:"+str(i) for i in pmid_list])))
		else:
			d_tmp_list.append("可能对{0}耐药".format(regimen_str))
	if d_tmp_list:
		BCD_evi.append("有研究表明，该类变异{0}。".format("，".join(d_tmp_list)))
	
	if A_evi:
		return "".join(A_evi)
	else:
		return "".join(BCD_evi)
jinja2.filters.FILTERS["xajdy_evi_summary_v2"] = xajdy_evi_summary_v2

# 西安交大一-新增NCCN指南推荐基因列表（区分癌种）-2025.06.03
def xajdy_nccn_result(info):
	nccn_list = info[0]
	var_list = info[1]
	#tumor_list = info[2]
	sample = info[2]
	prod_name = sample["prod_names"]
	# 1. 汇总癌种对应的基因列表
	nccn_tumor_gene = {}
	nccn_all_gene = []
	for i in nccn_list:
		if i["gene_symbol"] not in nccn_all_gene:
			nccn_all_gene.append(i["gene_symbol"])
		for tumor in re.split("、", i["disease"]):
			if tumor not in nccn_tumor_gene:
				nccn_tumor_gene.setdefault(tumor, [])
			nccn_tumor_gene[tumor].append(i["gene_symbol"])
	# 2. 筛选出需要展示的基因列表
	raw_gene_list_all = []
	# 2025.06.11-新增规则，如果癌种为 其他/实体瘤，则展示所有基因
	if sample["tumor_type"] in ["其他", "实体瘤"]:
		raw_gene_list_all = nccn_all_gene
	else:
		for tumor in sample["tumor_list"]:
			if tumor in nccn_tumor_gene.keys():
				raw_gene_list_all.extend(nccn_tumor_gene[tumor])
		raw_gene_list_all.extend(nccn_tumor_gene["实体瘤"])
	### 跟产品检测范围取交集-116/MP无需再过滤，其他产品有新增的往下再加
	lc10_gene_list = ["ALK", "BRAF", "EGFR", "ERBB2", "KRAS", "MET", "NRAS", "PIK3CA", "RET", "ROS1"]
	cp40_gene_list = ["AKT1", "ALK", "BRAF", "CDK4", "CTNNB1", "DDR2", "DPYD", "EGFR", \
				      "ERBB2", "ESR1", "FGFR1", "FGFR2", "FGFR3", "FGFR4", "HRAS", "IDH1", \
					  "IDH2", "KEAP1", "KIT", "KRAS", "MAP2K1", "MET", "MYC", "NFE2L2", \
					  "NKX2-1", "NRAS", "NRG1", "NTRK1", "NTRK2", "NTRK3", "PDGFRA", "PIK3CA", \
					  "POLE", "PTEN", "RB1", "RET", "ROS1", "STK11", "TP53", "UGT1A1"]
	tc21_gene_list = ["AKT1", "ALK", "BRAF", "CTNNB1", "EIF1AX", "GNAS", "HRAS", "KRAS", \
					  "NRAS", "NTRK1", "NTRK3", "PAX8", "PDGFRA", "PIK3CA", "PTEN", "RASAL1", \
					  "RET", "TERT", "TP53", "TSC2", "TSHR"]
	ga18_gene_list = ["AKT1", "BRAF", "EGFR", "ERBB2", "FGF19", "FGFR1", "FGFR2", "FGFR3", \
				      "HRAS", "KIT", "KRAS", "MET", "NRAS", "PDGFRA", "PIK3CA", "PTEN", \
					  "TP53"]
	brca_gene_list = ["BRCA1", "BRCA2"]
	if "10基因" in prod_name:
		raw_gene_list = list(set(lc10_gene_list) & set(raw_gene_list_all))
	elif prod_name == "Classic Panel":
		raw_gene_list = list(set(cp40_gene_list) & set(raw_gene_list_all))
	elif "TC21" in prod_name:
		raw_gene_list = list(set(tc21_gene_list) & set(raw_gene_list_all))
	elif "GA18" in prod_name:
		raw_gene_list = list(set(ga18_gene_list) & set(raw_gene_list_all))
	elif "BRCA" in prod_name:
		raw_gene_list = list(set(brca_gene_list) & set(raw_gene_list_all))
	else:
		raw_gene_list = raw_gene_list_all
	# 3. 基因列表去重下
	gene_list = []
	for gene in raw_gene_list:
		if gene not in gene_list:
			gene_list.append(gene)
	# 4. 汇总结果
	raw_result = [var for var in var_list if set(re.split(",", var["gene_symbol"])) & set(gene_list)]
	detect_gene = []
	for var in var_list:
		for gene in re.split(",", var["gene_symbol"]):
			detect_gene.append(gene)
	for gene in set(gene_list) - set(detect_gene):
		raw_result.append({"gene_symbol" : gene})
	raw_result = sorted(raw_result, key = lambda i:i["gene_symbol"])
	# 5. 结果分成两列展示
	result = []
	for i in range(0, len(raw_result)-2, 2):
		tmp_dict = {}
		for j in range(1, 3):
			tmp_dict["var"+str(j)] = raw_result[i+j-1]
		result.append(tmp_dict)

	rest_num= len(raw_result) % 2
	rest_tmp_dict = {}
	for j in range(1, 3):
		rest_tmp_dict["var"+str(j)] = ""

	num = 1
	last_row_num = len(raw_result)-rest_num if rest_num != 0 else len(raw_result)-rest_num-2
	for j in range(last_row_num, len(raw_result)):
		rest_tmp_dict["var"+str(num)] = raw_result[j]
		num += 1
	result.append(rest_tmp_dict)
	return result
jinja2.filters.FILTERS["xajdy_nccn_result"] = xajdy_nccn_result

# 西安交大一判断指定基因列表是否存在胚系4/5类变异-2025.06.04
def xajdy_germline_45(info):
	gene_type = info[0]
	gls5_gene = ["MSH2", "MSH6", "PMS2", "MLH1", "EPCAM"]
	var_list = info[1]
	gls5_45 = [var for var in var_list if var["gene_symbol"] in gls5_gene]
	other_45 = [var for var in var_list if var["gene_symbol"] not in gls5_gene]
	if gene_type == "gls5":
		return gls5_45
	else:
		return other_45
jinja2.filters.FILTERS["xajdy_germline_45"] = xajdy_germline_45

# 西安交大一-胚系项目新增NCCN指南推荐基因列表（不区分癌种，按指定基因列表返回）-2025.06.06
def xajdy_nccn_result_germline(info):
	var_list = info[0]
	prod_name = info[1]
	brca_gene_list = ["BRCA1", "BRCA2"]
	ptm_gene_list = ["POLE", "TP53"]
	hrd_gene_list = ["ATM", "BARD1", "BRCA1", "BRCA2", "BRIP1", "CDH1", \
				  	 "CDK12", "CHEK1", "CHEK2", "FANCA", "FANCL", "HDAC2", \
					 "PALB2", "PPP2R2A", "PTEN", "RAD51B", "RAD51C", "RAD51D", \
					 "RAD54L", "TP53"]
	gene61_gene_list = ["APC", "ATM", "BARD1", "BMPR1A", "BRAF", "BRCA1", \
					 	"BRCA2", "BRIP1", "CDH1", "CDK4", "CDKN2A", "CHEK2", \
						"ELAC2", "EPCAM", "FANCC", "FH", "FLCN", "GNAS", \
						"GREM1", "HOXB13", "HRAS", "KIT", "MAX", "MEN1", \
						"MET", "MLH1", "MRE11", "MSH2", "MSH6", "MUTYH", \
						"NBN", "NF1", "NTRK1", "PALB2", "PALLD", "PDGFRA", \
						"PMS2", "PRKAR1A", "PTCH1", "PTEN", "RAD50", "RAD51", \
						"RAD51C", "RAD51D", "RB1", "RET", "SDHA", "SDHAF2", \
						"SDHB", "SDHC", "SDHD", "SMAD4", "SMARCA4", "SMARCB1", \
						"STK11", "TMEM127", "TP53", "TSC1", "TSC2", "VHL", "WRN"]
	gene150_gene_list = ["ABRAXAS1","AIP","ALK","APC","ASXL1","ATM","ATR","ATRIP","AXIN2","BAP1",\
					     "BARD1","BLM","BMPR1A","BRAF","BRCA1","BRCA2","BRIP1","BUB1","BUB3","CDC73",\
						 "CDH1","CDK12","CDK4","CDKN1B","CDKN1C","CDKN2A","CFTR","CHEK1","CHEK2","CTNNA1",\
						 "DDB2","DDX41","DICER1","DIS3L2","DLST","EGFR","ELAC2","EPCAM","ERCC1","ERCC2",\
						 "ERCC3","ERCC4","ERCC5","EXT1","EXT2","FAN1","FANCA","FANCB","FANCC","FANCD2",\
						 "FANCE","FANCF","FANCG","FANCI","FANCL","FANCM","FH","FLCN","FOXE1","GALNT12",\
						 "GEN1","GNAS","GPR101","GREM1","HDAC2","HNF1A","HNF1B","HOXB13","HRAS","KIT",\
						 "LRP6","MAD2L2","MAX","MBD4","MC1R","MEN1","MET","MITF","MLH1","MLH3",\
						 "MRE11","MSH2","MSH3","MSH6","MUTYH","NBN","NF1","NF2","NTHL1","NTRK1",\
						 "PALB2","PALLD","PDGFRA","PMS2","POLD1","POLE","POLE2","POLH","POT1","PPP2R2A",\
						 "PRKAR1A","PRSS1","PTCH1","PTEN","RAD50","RAD51","RAD51B","RAD51C","RAD51D","RAD54L",\
						 "RB1","RBBP8","RECQL","RECQL4","REST","RET","RHBDF2","RNASEL","RNF43","RPS20",\
						 "SDHA","SDHAF2","SDHB","SDHC","SDHD","SLC25A11","SLX4","SMAD4","SMARCA2","SMARCA4",\
						 "SMARCB1","SPINK1","STK11","SUFU","TERT","TGFBR2","TMEM127","TP53","TRIM37","TRIP13",\
						 "TSC1","TSC2","UBE2T","VHL","WRN","WT1","XPA","XPC","XRCC2","XRCC3"]
	gene_list = []
	if "BRCA" in prod_name:
		gene_list = brca_gene_list
	elif prod_name == "PTM（3基因）":
		gene_list = ptm_gene_list
	elif prod_name == "HRD Complete（组织）" or prod_name == "M HRD（组织）":
		gene_list = hrd_gene_list
	elif prod_name == "61遗传基因":
		gene_list = gene61_gene_list
	elif prod_name == "遗传易感150基因":
		gene_list = gene150_gene_list
	# 汇总结果
	raw_result = [var for var in var_list if set(re.split(",", var["gene_symbol"])) & set(gene_list)]
	detect_gene = []
	for var in var_list:
		for gene in re.split(",", var["gene_symbol"]):
			detect_gene.append(gene)
	for gene in set(gene_list) - set(detect_gene):
		raw_result.append({"gene_symbol" : gene})
	raw_result = sorted(raw_result, key = lambda i:i["gene_symbol"])
	# 结果分成两列展示
	result = []
	for i in range(0, len(raw_result)-2, 2):
		tmp_dict = {}
		for j in range(1, 3):
			tmp_dict["var"+str(j)] = raw_result[i+j-1]
		result.append(tmp_dict)

	rest_num= len(raw_result) % 2
	rest_tmp_dict = {}
	for j in range(1, 3):
		rest_tmp_dict["var"+str(j)] = ""

	num = 1
	last_row_num = len(raw_result)-rest_num if rest_num != 0 else len(raw_result)-rest_num-2
	for j in range(last_row_num, len(raw_result)):
		rest_tmp_dict["var"+str(num)] = raw_result[j]
		num += 1
	result.append(rest_tmp_dict)
	return result
jinja2.filters.FILTERS["xajdy_nccn_result_germline"] = xajdy_nccn_result_germline

# 2025.06.09-云肿BPTM Plus组织、HRD-详细解析
# 一个基因展示一行，多突变的汇总所有证据，存在AB展示AB，否则展示CD，不展示预后和辅助诊断
def ynzl_germline_prod_var_inter(var_list):
	gene_rule = []
	gene_inter = {}
	for var in var_list:
		if var["gene_symbol"] not in gene_rule:
			gene_rule.append(var["gene_symbol"])
			gene_inter[var["gene_symbol"]] = var["gene_function"]
	gene_var_dict = {}
	for var in var_list:
		if var["gene_symbol"] not in gene_var_dict.keys():
			gene_var_dict.setdefault(var["gene_symbol"], [])
		if "Predictive" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Predictive"]:
			for evi in var["evi_sum"]["evi_split"]["Predictive"]:
				evi_tmp = {
					"regimen_name" : evi["regimen_name"],
					"evi_interpretation" : evi["evi_interpretation"],
					"evi_conclusion_simple" : evi["evi_conclusion_simple"]
				}
				if evi_tmp not in gene_var_dict[var["gene_symbol"]]:
					gene_var_dict[var["gene_symbol"]].append(evi_tmp)
	result = []
	for gene, evi_list in gene_var_dict.items():
		AB_evi = [i for i in evi_list if i["evi_conclusion_simple"] in ["A", "B"]]
		CD_evi = [i for i in evi_list if i["evi_conclusion_simple"] in ["C", "D"]]
		if AB_evi:
			result.append({
				"gene_symbol" : gene,
				"gene_function" : gene_inter.get(gene, ""),
				"evi" : ynzl_merge_evi(AB_evi)
			})
		else:
			result.append({
				"gene_symbol" : gene,
				"gene_function" : gene_inter.get(gene, ""),
				"evi" : ynzl_merge_evi(CD_evi)
			})
	return sorted(result, key=lambda i:gene_rule.index(i["gene_symbol"]))
jinja2.filters.FILTERS["ynzl_germline_prod_var_inter"] = ynzl_germline_prod_var_inter

# 相同证据描述的合并治疗方案展示--云肿
def ynzl_merge_evi(datainfo):
	tmp_dict = {}
	for evi in datainfo:
		tmp_dict.setdefault(evi["evi_interpretation"], [])
		tmp_dict[evi["evi_interpretation"]].append(evi["regimen_name"])

	merge_result = []
	for k, v  in tmp_dict.items():
		merge_result.append(
			{
				"regimen_name" : "、".join(v),
				"evi_interpretation" : k
			}
		)

	return merge_result

# 重庆附一116小结展示具体变异-hgvs_p不带括号-2025.06.10
def cqfy_116_var_sum_s(var_list):
	result = []
	for var in var_list:
		if var["bio_category"] == "Snvindel":
			if var["hgvs_p"] != "p.?":
				result.append(var["gene_symbol"]+"基因 "+var["hgvs_p_ZJZL"]+var["varInfo_XAJDY"])
			else:
				result.append(var["gene_symbol"]+" "+var["hgvs_c"]+var["varInfo_XAJDY"])
		elif var["bio_category"] == "Cnv":
			result.append(var["gene_symbol"]+"基因扩增")
		elif var["bio_category"] in ["Sv", "PSeqRnaSv"]:
			if var["five_prime_gene"] == "MET" and var["three_prime_gene"] == "MET":
				result.append("MET exon14 跳跃")
			else:
				# 融合可能会有重复（rna exon相同，断点不同的情况）
				sv_info = "{0}-{1}({2}:{3})融合".format(var["five_prime_gene"], var["three_prime_gene"], \
										  				var["five_prime_gene"][0] + var["five_prime_cds"].replace("exon", ""),\
														var["three_prime_gene"][0] + var["three_prime_cds"].replace("exon", ""))
				if sv_info not in result:
					result.append(sv_info)
	return ", ".join(result)
jinja2.filters.FILTERS["cqfy_116_var_sum_s"] = cqfy_116_var_sum_s

# 云肿结果小结-相同等级变异、相同基因的要汇总到一起展示-2025.06-10
# hgvs_p改为带括号的
def ynzl_summary_var_v2(var_list):
	'''
	ALK基因：p.(G1202R)、ALK扩增、EML4-ALK融合（EML4:exon13-ALK:exon20; EML4:exon13-ins39-ALK:exon20）
	'''
	# 获取基因列表
	gene_list = []
	for var in var_list:
		if var["gene_symbol"] not in gene_list:
			gene_list.append(var["gene_symbol"])
	
	result = {}
	for var in var_list:
		var_info = ""
		if var["bio_category"] == "Snvindel":
			var_info = var["hgvs_p"] if var["hgvs_p"] != "p.?" else var["hgvs_c"]
			if "judge_mergeMET" in var.keys() and var["judge_mergeMET"]:
				var_info += "（MET exon14 skipping）"
		elif var["bio_category"] == "Cnv":
			var_info = "扩增"
		else:
			if var["five_prime_gene"] == "MET" and var["three_prime_gene"] == "MET":
				var_info = "MET exon14 skipping"
			else:
				if "var_desc_merge" in var.keys() and var["var_desc_merge"]:
					var_desc_merge = "；".join(re.split("、", var["var_desc_merge"]))
				else:
					var_desc_merge = "{0}:{1}-{2}:{3}".format(var["five_prime_gene"], var["five_prime_cds"], var["three_prime_gene"], var["three_prime_cds"])
				var_info = "{0}-{1}融合（{2}）".format(var["five_prime_gene"], var["three_prime_gene"], var_desc_merge)
		if var["gene_symbol"] not in result.keys():
			result.setdefault(var["gene_symbol"], [])
		result[var["gene_symbol"]].append(var_info)

	result_list = []
	for gene, info in result.items():
		result_list.append(
			{
				"gene_symbol" : gene,
				"var_info" : "、".join(info)
			}
		)
	result_list = sorted(result_list, key=lambda i:gene_list.index(i["gene_symbol"]))
	return result_list
jinja2.filters.FILTERS["ynzl_summary_var_v2"] = ynzl_summary_var_v2

# 2025.06.11-新增-嵇梦晨-适用山东肿瘤
def sum_var(var_list):
	v_result = []
	for var in var_list:
		if var["bio_category"] == "Snvindel":
			if var["hgvs_p"] != "p.?":
				v_result.append(var["gene_symbol"]+" "+var["hgvs_p"])
			else:
				v_result.append(var["gene_symbol"]+" "+var["hgvs_c"])
		elif var["bio_category"] == "Cnv":
			# 2025.06.20-区分loss-嵇梦晨
			# 缺失兼容小写loss
			if var["cnv_type"] in ["Loss", "loss"]:
				v_result.append(var["gene_symbol"] + " 缺失")
			else:
				v_result.append(var["gene_symbol"]+" 扩增")
		elif var["bio_category"] in ["Sv", "PSeqRnaSv"]:
			if var["five_prime_gene"] == "MET" and var["three_prime_gene"] == "MET":
				v_result.append("MET exon14 跳跃")
			else:
				if var["five_prime_gene"]+"-"+var["three_prime_gene"]+" 融合" not in v_result:
					v_result.append(var["five_prime_gene"]+"-"+var["three_prime_gene"]+" 融合")
	return ", ".join(v_result)
jinja2.filters.FILTERS["sum_var"] = sum_var

def var_filter_40(var_list):
	'''
	OncoPro Panel
	CP40 40基因单独展示
	'''
	gene_list = ["AKT1", "ALK", "BRAF", "CDK4", "CTNNB1", "DDR2", "DPYD", "EGFR", "ERBB2", 
    "ESR1", "FGFR1", "FGFR2", "FGFR3", "FGFR4", "HRAS", "IDH1", "IDH2", 
    "KEAP1", "KIT", "KRAS", "MAP2K1", "MET", "MYC", "NFE2L2", "NKX2-1", 
    "NRAS", "NRG1", "NTRK1", "NTRK2", "NTRK3", "PDGFRA", "PIK3CA", "POLE", 
    "PTEN", "RB1", "RET", "ROS1", "STK11", "TP53", "UGT1A1"]
	# 不知道V3和V4 SV的gene_symbol怎么返回
	result = []
	for var in var_list:
		if var["bio_category"] in ["Sv", "PSeqRnaSv"]:
			# 2025.09.19-刘炜芬-融合规则改一下，gene_symbol在40基因内的报
			if set(re.split(",", var["gene_symbol"])) & set(gene_list):
				result.append(var)
			#if var["five_prime_gene"] in gene_list or var["three_prime_gene"] in gene_list:
			#	result.append(var)
			# 2025.09.19-更新完成
		else:
			if var["gene_symbol"] in gene_list:
				result.append(var)
	return result
jinja2.filters.FILTERS["var_filter_40"] = var_filter_40

def var_filter_160(var_list):
	gene_list = ["AKT1", "ALK", "BRAF", "CDK4", "CTNNB1", "DDR2", "DPYD", "EGFR", "ERBB2", 
    "ESR1", "FGFR1", "FGFR2", "FGFR3", "FGFR4", "HRAS", "IDH1", "IDH2", 
    "KEAP1", "KIT", "KRAS", "MAP2K1", "MET", "MYC", "NFE2L2", "NKX2-1", 
    "NRAS", "NRG1", "NTRK1", "NTRK2", "NTRK3", "PDGFRA", "PIK3CA", "POLE", 
    "PTEN", "RB1", "RET", "ROS1", "STK11", "TP53", "UGT1A1"]
	result = []
	for var in var_list:
		if var["bio_category"] in ["Sv", "PSeqRnaSv"]:
			# 2025.09.19-刘炜芬-融合规则改一下，gene_symbol不在40基因内的报
			if not set(re.split(",", var["gene_symbol"])) & set(gene_list):
				result.append(var)
			#if not (var["five_prime_gene"] in gene_list and var["three_prime_gene"] in gene_list):
			#	result.append(var)
			# 2025.09.19-更新完成
		else:
			if var["gene_symbol"] not in gene_list:
				result.append(var)
	return result
jinja2.filters.FILTERS["var_filter_160"] = var_filter_160
# 2025.06.11-新增完成

# 2025.06.12-重庆西南116-过滤10基因获批位点
def cqxn_10gene_appr_var(var_list):
	gene_list = ["ALK", "BRAF", "EGFR", "ERBB2", "KRAS", "MET", "NRAS", "PIK3CA", "RET", "ROS1"]
	appr_dict = {
		"EGFR" : ["p.G719A","p.E746_A750del","p.L747_T751del","p.V769_D770insASV","p.V774_C775insHV",\
				  "p.H773_V774insAH","p.H773_V774insY","p.P772_H773insPNP","p.V769_D770insGG","p.N771_P772insRHN",\
				  "p.L861Q","p.G719S","p.L747_P753delinsS","p.T790M","p.D770_N771insSVD",\
				  "p.D770delinsGY","p.N771_P772insN","p.H773_V774insTH","p.N771_P772insG","p.V769_D770insGTL",\
				  "p.P772_H773insGHP","p.G719C","p.E746_S752delinsV","p.S768I","p.A763_Y764insFQEA",\
				  "p.H773_V774insH","p.N771_P772insH","p.N771delinsGY","p.N771delinsGF","p.D770_N771insGF",\
				  "p.H773delinsYNPY","p.E746_A750del","p.L747_A750delinsP","p.D770_N771insG","p.H773_V774insNPH",\
				  "p.H773_V774insPH","p.P772_H773insYNP","p.P772_H773insQ","p.P772_H773insGNP","p.D770_N771insSLA","p.L858R"],
		"KRAS" : ["p.G12D", "p.G12A", "p.G12V", "p.G12S", "p.G12C", "p.Q61H"],
		"NRAS" : ["p.G12D", "p.Q61R", "p.Q61K"],
		"PIK3CA" : ["p.H1047R"],
		"BRAF" : ["p.V600E"],
		"ERBB2" : ["p.A775_G776insYVMA"],
		"MET" : ["c.3082+1G>T"],
		"ALK" : ["EML4:exon13-ALK:exon20", "EML4:exon6-ALK:exon20", "EML4:exon20-ALK:exon20"],
		"ROS1" : ["CD74:exon6-ROS1:exon34", "GOPC:exon8-ROS1:exon35"],
		"RET" : ["KIF5B:exon15-RET:exon12"]
	}
	appr_var = []
	for var in var_list:
		if var["bio_category"] == "Snvindel":
			snv_info = var["hgvs_p_ZJZL"] if var["hgvs_p_ZJZL"] != "p.?" else var["hgvs_c"]
			if var["gene_symbol"] in appr_dict.keys() and snv_info in appr_dict[var["gene_symbol"]]:
				appr_var.append(var)
		elif var["bio_category"] == "Sv":
			sv_info = "{0}:{1}-{2}:{3}".format(var["five_prime_gene"], var["five_prime_cds"], var["three_prime_gene"], var["three_prime_cds"])
			# ALK/ROS1/RET几个获批变异型，主基因都在3端
			if var["three_prime_gene"] in appr_dict.keys() and sv_info in appr_dict[var["three_prime_gene"]]:
				appr_var.append(var)
	detect_gene = []
	for var in appr_var:
		for gene in re.split(",", var["gene_symbol"]):
			if gene in gene_list and gene not in detect_gene:
				detect_gene.append(gene)
	result = []
	result.extend(appr_var)
	for gene in sorted(list(set(gene_list) - set(detect_gene))):
		result.append({"gene_symbol" : gene})
	return result
jinja2.filters.FILTERS["cqxn_10gene_appr_var"] = cqxn_10gene_appr_var

# 武汉协和1/2类变异需要过滤出符合下面规则的变异展示到3类中-2025.06.13
# >=1个BA的不展示到III类
# >=2个BS的不展示到III类
# 不符合上述两条并且证据含有致病性的（首字母为P）的展示到III类
def whxh_12var_filter(var_list):
	BA_list = ["BA1","BS1_Alone","BS2_Alone","BP6_Alone","BS4_Alone","BP5_Alone","BS3_Alone",\
			   "BP2_Alone","BP1_Alone","BP3_Alone","BP4_Alone","BP7_Alone","BP7_Alone(RNA)"]
	BS_list = ["BA1_Strong","BS1","BS2","BP6_Strong","BS4","BP5_Strong","BS3","BP2_Strong",\
			   "BP1_Strong","BP3_Strong","BP4_Strong","BP7_Strong","BP7_Strong(RNA)"]
	result = []
	for var in var_list:
		patho = False
		evi_category = re.split(";", var["evidence_categorys"]) if "evidence_categorys" in var.keys() and var["evidence_categorys"] else []
		BA_num = len(set(evi_category) & set(BA_list))
		BS_num = len(set(evi_category) & set(BS_list))
		if BA_num >= 1:
			continue
		elif BS_num >= 2:
			continue
		else:
			for i in evi_category:
				if i[0] == "P":
					patho = True
					break
		#print (var["gene_symbol"], var["hgvs_c"], var["hgvs_p"], "BA_num",BA_num,"BS_num",BS_num, patho)
		if patho:
			result.append(var)
	return result
jinja2.filters.FILTERS["whxh_12var_filter"] = whxh_12var_filter

# 武汉协和1/2类变异需要过滤出符合下面规则的变异不展示到良性变异表格中-2025.06.16
# >=1个BA的不展示到III类
# >=2个BS的不展示到III类
# 不符合上述两条并且证据含有致病性的（首字母为P）的展示到III类
def whxh_12var_filter_opposite(var_list):
	BA_list = ["BA1","BS1_Alone","BS2_Alone","BP6_Alone","BS4_Alone","BP5_Alone","BS3_Alone",\
			   "BP2_Alone","BP1_Alone","BP3_Alone","BP4_Alone","BP7_Alone","BP7_Alone(RNA)"]
	BS_list = ["BA1_Strong","BS1","BS2","BP6_Strong","BS4","BP5_Strong","BS3","BP2_Strong",\
			   "BP1_Strong","BP3_Strong","BP4_Strong","BP7_Strong","BP7_Strong(RNA)"]
	result = []
	for var in var_list:
		evi_category = re.split(";", var["evidence_categorys"]) if "evidence_categorys" in var.keys() and var["evidence_categorys"] else []
		BA_num = len(set(evi_category) & set(BA_list))
		BS_num = len(set(evi_category) & set(BS_list))
		if BA_num >= 1:
			result.append(var)
		elif BS_num >= 2:
			result.append(var)
		else:
			patho = False
			for i in evi_category:
				if i[0] == "P":
					patho = True
					break
		#print (var["gene_symbol"], var["hgvs_c"], var["hgvs_p"], "BA_num",BA_num,"BS_num",BS_num, patho)
			if not patho:
				result.append(var)
	return result
jinja2.filters.FILTERS["whxh_12var_filter_opposite"] = whxh_12var_filter_opposite

# 武汉协和HRR结果汇总-2025.06.13
def whxh_hrr_summary(raw_var_list):
	gene_list = ["ATM","BARD1","BRCA1","BRCA2","BRIP1","CDK12","CHEK1","CHEK2","FANCL",\
			  	 "PALB2","RAD51B","RAD51C","RAD51D","RAD54L", "FANCA", "ATR", "MRE11", "NBN"]
	var_list = [var for var in raw_var_list if var["gene_symbol"] in gene_list]
	detect_gene = [var["gene_symbol"] for var in var_list if var["gene_symbol"] in gene_list]
	result = []
	result.extend(var_list)	
	for gene in set(gene_list) - set(detect_gene):
		result.append({
			"gene_symbol" : gene
		})	
	return sorted(result, key=lambda i:i["gene_symbol"])
jinja2.filters.FILTERS["whxh_hrr_summary"] = whxh_hrr_summary

# 华西日期补充0，如“2025/5/1”改为“2025-05-01” - 2025.06.16
def schx_date(indate):
	if indate:
		indate = indate.replace("-", "/")
		date_item = re.split("/", indate)
		if len(date_item) == 3:
			m_item = "0" + date_item[1] if len(date_item[1]) == 1 else date_item[1]
			d_item = "0" + date_item[2] if len(date_item[2]) == 1 else date_item[2]
			return "{0}-{1}-{2}".format(date_item[0], m_item, d_item)
		else:
			return indate
	else:
		return ""
jinja2.filters.FILTERS["schx_date"] = schx_date

# 武汉协和HRD结果汇总-2025.06.16
def whxh_hrd_summary(raw_var_list):
	gene_list = ["ATM","BARD1","BRCA1","BRCA2","BRIP1","CDK12","CHEK1","CHEK2",\
			     "FANCL","PALB2","RAD51B","RAD51C","RAD51D","RAD54L", "FANCA"]
	var_list = [var for var in raw_var_list if var["gene_symbol"] in gene_list]
	detect_gene = [var["gene_symbol"] for var in var_list if var["gene_symbol"] in gene_list]
	result = []
	result.extend(var_list)	
	for gene in set(gene_list) - set(detect_gene):
		result.append({
			"gene_symbol" : gene
		})	
	return sorted(result, key=lambda i:i["gene_symbol"])
jinja2.filters.FILTERS["whxh_hrd_summary"] = whxh_hrd_summary

# 吉大一日期-2天-2025.06.17
def jdy_date_delta(indate):
	date = datetime.strptime(indate, "%Y-%m-%d")
	new_date = date - timedelta(days = 2)
	return new_date.strftime("%Y-%m-%d")
jinja2.filters.FILTERS["jdy_date_delta"] = jdy_date_delta

# 吉大一116新增子宫内膜癌分子分型，POLE需要具体变异-2025.06.17
def jdy_116_ec_pole(var_list):
	result = []
	if var_list:
		for var in var_list:
			if var["hgvs_p_ZJZL"] != "p.?":
				result.append(var["hgvs_p_ZJZL"])
			else:
				result.append(var["hgvs_c"])
	return "、".join(result)
jinja2.filters.FILTERS["jdy_116_ec_pole"] = jdy_116_ec_pole

# 温附一CP40-肠癌小结处输出KRAS/NRAS exon2-4 I/II/肿瘤发生发展变异+BRAF V600E变异-2025.06.20
# 均为snvindel
# 肠癌除了KNB，其他变异正常输出
def wfy_knb_var_list(info):
	var_list = info[0]
	var_type = info[1]
	gene_sort_rule = ["KRAS", "NRAS", "BRAF"]
	result_knb = []
	result_other = []
	detect_gene = []
	for var in var_list:
		if var["bio_category"] == "Snvindel":
			if var["gene_symbol"] in ["KRAS", "NRAS"] and set(["exon2", "exon3", "exon4"]) & set(re.split("_", var["gene_region"])):
				result_knb.append(var)
				detect_gene.append(var["gene_symbol"])
			elif var["gene_symbol"] == "BRAF" and var["hgvs_p"] == "p.(V600E)":
				result_knb.append(var)
				detect_gene.append(var["gene_symbol"])
			else:
				result_other.append(var)
		else:
			result_other.append(var)
	for gene in set(gene_sort_rule) - set(detect_gene):
		result_knb.append(
			{"gene_symbol" : gene}
		)
	if var_type == "knb":
		return sorted(result_knb, key = lambda i:gene_sort_rule.index(i["gene_symbol"]))
	else:
		return result_other
jinja2.filters.FILTERS["wfy_knb_var_list"] = wfy_knb_var_list

# 广东医科附属-新增NCCN指南推荐基因列表（区分癌种）-2025.06.30
def gdykfs_nccn_result(info):
	nccn_list = info[0]
	var_list = info[1]
	sample = info[2]
	# 1. 汇总癌种对应的基因列表
	nccn_tumor_gene = {}
	nccn_all_gene = []
	for i in nccn_list:
		if i["gene_symbol"] not in nccn_all_gene:
			nccn_all_gene.append(i["gene_symbol"])
		for tumor in re.split("、", i["disease"]):
			if tumor not in nccn_tumor_gene:
				nccn_tumor_gene.setdefault(tumor, [])
			nccn_tumor_gene[tumor].append(i["gene_symbol"])
	# 2. 筛选出需要展示的基因列表
	raw_gene_list_all = []
	if sample["tumor_type"] in ["其他", "实体瘤"]:
		raw_gene_list_all = nccn_all_gene
	else:
		for tumor in sample["tumor_list"]:
			if tumor in nccn_tumor_gene.keys():
				raw_gene_list_all.extend(nccn_tumor_gene[tumor])
		raw_gene_list_all.extend(nccn_tumor_gene["实体瘤"])
	# 3. 基因列表去重下
	gene_list = []
	for gene in raw_gene_list_all:
		if gene not in gene_list:
			gene_list.append(gene)
	# 4. 汇总结果
	raw_result = [var for var in var_list if set(re.split(",", var["gene_symbol"])) & set(gene_list)]
	detect_gene = []
	for var in var_list:
		for gene in re.split(",", var["gene_symbol"]):
			detect_gene.append(gene)
	for gene in set(gene_list) - set(detect_gene):
		raw_result.append({"gene_symbol" : gene})
	raw_result = sorted(raw_result, key = lambda i:i["gene_symbol"])
	return raw_result
jinja2.filters.FILTERS["gdykfs_nccn_result"] = gdykfs_nccn_result

# 南昌附一证据来源治疗方案没有获批机构的，将研究类型进行转化-2025.07.01
def ncfy_evidence_level_stran(evi_level):
	evi_level_dict = {
		"Clinical-phase I" : "临床试验",
		"Clinical-phase II" : "临床试验",
		"Clinical-phase III" : "临床试验",
		"Clinical-phase IV" : "临床试验",
		"Clinical-retrospective" : "临床试验",
		"Clinical-unknown phase" : "临床试验",
		"Case report" : "案例报道",
		"Preclinical-in vitro" : "临床前研究",
		"Preclinical-in vivo" : "临床前研究"
		}
	return evi_level_dict.get(evi_level, evi_level) if evi_level else ""
jinja2.filters.FILTERS["ncfy_evidence_level_stran"] = ncfy_evidence_level_stran

# 广东医科CP200，变异小结处变异需要分三列展示-2025.07.03-孟智悦
def split_into_three_columns(var_data):
	'''
	将变异拆分为一行三个变异展示
	'''
	# 构建变异列表（处理RNA多变异情况）
	var_list = []
	for var in var_data:
		var_list.append(var)
		if "dup_rna_detect" in var.keys() and len(var["dup_rna_detect"]) > 1:
			for rna_var in var["dup_rna_detect"]:
				var_list.append(rna_var)

	result = []
	# 每3个变异为一组进行分组
	for i in range(0, len(var_list), 3):
		row = {}
		# 当前组的变异数量
		group_size = min(3, len(var_list) - i)

		# 填充当前组的变异（1-3个）
		for j in range(1, 4):
			idx = i + j - 1
			if j <= group_size:
				row[f"var{j}"] = var_list[idx]
			else:
				row[f"var{j}"] = ""  # 不足3个时用空字符串填充

		result.append(row)

	return result
jinja2.filters.FILTERS["split_into_three_columns"] = split_into_three_columns

#广东医科附属-展示变异的解读结果-2025.07.03-孟智悦
# -刘炜芬-snvindel function_classification和clinical_significance只返回1个（两个都返回也要兼容下，都返回体系就看function，胚系看clinical）
# 可能存在3类有药的情况，这边补齐判定
def gdykfs_evi_sum(var):
	clinic_inter = ""
	level_stran = {
		"Oncogenic" : "致癌性变异（Oncogenic）",
		"Likely oncogenic" : "可能致癌性变异（Likely Oncogenic）",
		"Pathogenic" : "致病性变异（Pathogenic）",
		"Likely pathogenic" : "疑似致病性变异（Likely pathogenic）",
		"Uncertain" : "意义不明确变异（Uncertain）",
		"Likely benign" : "疑似良性变异（Likely benign）",
		"Benign" : "良性变异（Benign）"
	}
	# 处理 SNV/INDEL 类型
	if var["bio_category"] == "Snvindel":
		if var["var_origin"] == "germline":
			if var["clinical_significance"] != "-":
				clinic_inter = level_stran.get(var["clinical_significance"], "")
			elif var["function_classification"] != "-":
				clinic_inter = level_stran.get(var["function_classification"], "")
		else:
			if var["function_classification"] != "-":
				clinic_inter = level_stran.get(var["function_classification"], "")
			elif var["clinical_significance"] != "-":
				clinic_inter = level_stran.get(var["clinical_significance"], "")
	# 处理 CNV 类型
	elif var["bio_category"] == "Cnv":
		if var["var_origin"] == "germline":
			clinic_inter = level_stran.get(var["clinical_significance"], "")
		else:
			clinic_inter = level_stran.get(var["function_classification"], "")
	# 处理 SV 类型
	elif var["bio_category"] == "Sv":
		clinic_inter = level_stran.get(var["function_classification"], "")
		
	return clinic_inter
jinja2.filters.FILTERS["gdykfs_evi_sum"] = gdykfs_evi_sum

# 广东医科附属变异解读需要展示靶向/预后/辅助诊断总结-2025.07.03-孟智悦
def gdykfs_evi_summary(evi_list):
	# A级证据处理
	A_evi = []
	A_sense_regimen = []
	A_resis_evi = []
	for evi in evi_list:
		if evi["evi_conclusion_simple"] == "A":
			if evi["clinical_significance_cn"] == "敏感":
				if evi["regimen_name"] not in A_sense_regimen:
					A_sense_regimen.append(evi["regimen_name"])
			else:
				if evi["evi_interpretation"] not in A_resis_evi:
					A_resis_evi.append(evi["evi_interpretation"])
	if A_sense_regimen:
		A_evi.append("针对携带该变异的患者，获批或指南推荐的药物有{0}等。".format("、".join(A_sense_regimen)))
	if A_resis_evi:
		A_evi.extend(A_resis_evi)
	A_str = "".join(A_evi)

	#C3-敏感证据处理
	# C3无辅助诊断和预后的情况
	C3_evi = []
	C3_sense_tumor_and_regimen = {}
	C3_resis_evi = []
	for evi in evi_list:
		tumor_name = evi["tumor_name_cn"] if "tumor_name_cn" in evi.keys() and evi["tumor_name_cn"] else evi[
			"tumor_name_en"] if "tumor_name_en" in evi.keys() and evi["tumor_name_en"] else ""
		if evi["evi_conclusion"] == "C3" and evi["clinical_significance_cn"] == "敏感":
			if tumor_name not in C3_sense_tumor_and_regimen.keys():
				C3_sense_tumor_and_regimen.setdefault(tumor_name, [])
			C3_sense_tumor_and_regimen[tumor_name].append(evi["regimen_name"])
		elif evi["evi_conclusion"] == "C3" and evi["clinical_significance_cn"] == "耐药":
			if evi["evi_interpretation"] not in C3_resis_evi:
				C3_resis_evi.append(evi["evi_interpretation"])
	if C3_sense_tumor_and_regimen:
		tmp_list = []
		for tumor, regimen in C3_sense_tumor_and_regimen.items():
			tmp_list.append("{0}用于{1}".format("、".join(regimen), tumor))
		C3_evi.append("针对该突变，其他癌种中已有获批或指南推荐的药物，如{0}等。".format("，".join(tmp_list)))
	if C3_resis_evi:
		C3_evi.extend(C3_resis_evi)
	C3_str = "".join(C3_evi)

	# 其他B/C级证据处理
	BC_evi = []
	b_c_other_sense_tumor_and_regimen = {}
	b_c_resis_tumor_and_regimen = {}
	b_c_pron_tumor_and_clinical_significance_cn = {}
	b_c_diag_tumor = {}
	for evi in evi_list:
		pmid_list = get_pmid_from_inter(evi["evi_interpretation"])
		tumor_name = evi["tumor_name_cn"] if "tumor_name_cn" in evi.keys() and evi["tumor_name_cn"] else evi[
			"tumor_name_en"] if "tumor_name_en" in evi.keys() and evi["tumor_name_en"] else ""
		clinical_significance_cn = evi["clinical_significance_cn"]
		# 其他B/C/D需要PMID
		if evi["evi_conclusion_simple"] in ["C", "B"] and evi["evi_conclusion"] != "C3":
			if evi["clinical_significance_cn"] == "敏感":
				if tumor_name not in b_c_other_sense_tumor_and_regimen.keys():
					b_c_other_sense_tumor_and_regimen.setdefault(tumor_name, [])
				b_c_other_sense_tumor_and_regimen[tumor_name].append((evi["regimen_name"], pmid_list))
			elif evi["clinical_significance_cn"] == "耐药":
				if tumor_name not in b_c_resis_tumor_and_regimen.keys():
					b_c_resis_tumor_and_regimen.setdefault(tumor_name, [])
				b_c_resis_tumor_and_regimen[tumor_name].append((evi["regimen_name"], pmid_list))
			elif evi["clinical_significance_cn"] in ["较好", "较差", "中等"]:
				b_c_pron_tumor_and_clinical_significance_cn.setdefault(tumor_name, []).append((evi["clinical_significance_cn"], pmid_list))
				#print (b_c_pron_tumor_and_clinical_significance_cn)
			elif evi["clinical_significance"] in ["Positive", "Negative"]:
				b_c_diag_tumor.setdefault(tumor_name, []).append(pmid_list)

	b_c_tmp_list = []
	if b_c_other_sense_tumor_and_regimen:
		for tumor, regimen in b_c_other_sense_tumor_and_regimen.items():
			regimen_str = "、".join([i[0] for i in regimen])
			pmid_list = []
			for item in regimen:
				for pmid in item[1]:
					if pmid not in pmid_list:
						pmid_list.append(pmid)
			if pmid_list:
				b_c_tmp_list.append("{0}在携带该类变异的{1}中可能敏感（{2}）".format(regimen_str, tumor, ", ".join(
					["PMID:" + str(i) for i in pmid_list])))
			else:
				b_c_tmp_list.append("{0}在携带该类变异的{1}中可能敏感".format(regimen_str, tumor))
	if b_c_resis_tumor_and_regimen:
		for tumor, regimen in b_c_resis_tumor_and_regimen.items():
			regimen_str = "、".join([i[0] for i in regimen])
			pmid_list = []
			for item in regimen:
				for pmid in item[1]:
					if pmid not in pmid_list:
						pmid_list.append(pmid)
			if pmid_list:
				b_c_tmp_list.append("{0}在携带该类变异的{1}中可能耐药（{2}）".format(regimen_str, tumor, ", ".join(
					["PMID:" + str(i) for i in pmid_list])))
			else:
				b_c_tmp_list.append("{0}在携带该类变异的{1}中可能耐药".format(regimen_str, tumor))
	
	# 辅助诊断和预后只默认有一条，多了可能报告生成会出问题
	# 这边加个无PMID的情况-刘炜芬-2025.07.03
	if b_c_pron_tumor_and_clinical_significance_cn:
		if pmid_list:
			b_c_tmp_list.append(f"该变异提示预后{clinical_significance_cn}（{ ', '.join(['PMID:' + str(i) for i in pmid_list])}）")
		else:
			b_c_tmp_list.append(f"该变异提示预后{clinical_significance_cn}")
	
	if b_c_diag_tumor:
		if pmid_list:
			b_c_tmp_list.append(f"该变异与{tumor_name}的辅助诊断有关（{ ', '.join(['PMID:' + str(i) for i in pmid_list])}）")
		else:
			b_c_tmp_list.append(f"该变异与{tumor_name}的辅助诊断有关")


	if b_c_tmp_list:
		BC_evi.append("临床试验表明，{0}。".format("，".join(b_c_tmp_list)))

	BC_str = "".join(BC_evi)

    #D级证据处理
	D_evi = []
	d_sense_regimen = [(evi["regimen_name"], get_pmid_from_inter(evi["evi_interpretation"])) for evi in evi_list if
					   evi["evi_conclusion_simple"] == "D" and evi["clinical_significance_cn"] == "敏感"]
	d_resis_regimen = [(evi["regimen_name"], get_pmid_from_inter(evi["evi_interpretation"])) for evi in evi_list if
					   evi["evi_conclusion_simple"] == "D" and evi["clinical_significance_cn"] == "耐药"]
	d_pron_clinical_significance_cn = [(evi["clinical_significance_cn"], get_pmid_from_inter(evi["evi_interpretation"])) for evi in evi_list if
					   evi["evi_conclusion_simple"] == "D" and evi["clinical_significance_cn"] in ["较好", "较差", "中等"]]
	d_diag_tumor = [(get_pmid_from_inter(evi["evi_interpretation"])) for evi in evi_list if
					   evi["evi_conclusion_simple"] == "D" and evi["clinical_significance"] in ["Positive", "Negative"]]
	pmid_list = get_pmid_from_inter(evi["evi_interpretation"])
	tumor_name = evi["tumor_name_cn"] if "tumor_name_cn" in evi.keys() and evi["tumor_name_cn"] else evi[
		"tumor_name_en"] if "tumor_name_en" in evi.keys() and evi["tumor_name_en"] else ""
	clinical_significance_cn = evi["clinical_significance_cn"]

	d_tmp_list = []
	if d_sense_regimen:
		regimen_str = "、".join([i[0] for i in d_sense_regimen])
		pmid_list = []
		for item in d_sense_regimen:
			for pmid in item[1]:
				if pmid not in pmid_list:
					pmid_list.append(pmid)
		if pmid_list:
			d_tmp_list.append(
				"可能对{0}敏感（{1}）".format(regimen_str, ", ".join(["PMID:" + str(i) for i in pmid_list])))
		else:
			d_tmp_list.append("可能对{0}敏感".format(regimen_str))
	if d_resis_regimen:
		regimen_str = "、".join([i[0] for i in d_resis_regimen])
		pmid_list = []
		for item in d_resis_regimen:
			for pmid in item[1]:
				if pmid not in pmid_list:
					pmid_list.append(pmid)
		if pmid_list:
			d_tmp_list.append(
				"可能对{0}耐药（{1}）".format(regimen_str, ", ".join(["PMID:" + str(i) for i in pmid_list])))
		else:
			d_tmp_list.append("可能对{0}耐药".format(regimen_str))
	# 辅助诊断和预后只默认有一条，多了可能报告生成会出问题
	# 这边加个无PMID的情况-刘炜芬-2025.07.03
	if d_pron_clinical_significance_cn:
		if pmid_list:
			d_tmp_list.append(f"提示预后{clinical_significance_cn}（{ ', '.join(['PMID:' + str(i) for i in pmid_list])}）")
		else:
			d_tmp_list.append(f"提示预后{clinical_significance_cn}")
	if d_diag_tumor:
		if pmid_list:
			d_tmp_list.append(f"与{tumor_name}的辅助诊断有关（{ ', '.join(['PMID:' + str(i) for i in pmid_list])}）")
		else:
			d_tmp_list.append(f"与{tumor_name}的辅助诊断有关")

	if d_tmp_list:
		D_evi.append("有研究表明，该类变异{0}。".format("，".join(d_tmp_list)))
	D_str =  "".join(D_evi)

#将A_str, C3_str, BC_str, D_str合并展示
	# parts = []
	# if A_str: parts.append(A_str)
	# if C3_str: parts.append(C3_str)
	# if BC_str: parts.append(BC_str)
	# if D_str: parts.append(D_str)
	#
	# # 用句号连接各部分（避免重复句号）
	# result = ""
	# for p in parts:
	# 	if p.endswith('。') or p.endswith('.') or not result:
	# 		result += p
	# 	else:
	# 		result += '。' + p  # 确保段落分隔
	#
	# return result

	return A_str, C3_str, BC_str, D_str

jinja2.filters.FILTERS["gdykfs_evi_summary"] = gdykfs_evi_summary

# 济宁医学院附属-胃肠道间质瘤-报告中拆分为胃肠道间质瘤相关基因、LC10和其他基因 - 2025.07.03
def jnfy_ga_filter_lc10(info):
	var_list = info[0]
	var_type = info[1]
	ga_snvindel = ["KIT", "PDGFRA", "BRAF", "NF1", "KRAS", "PIK3CA"]
	ga_sv = ["FGFR1", "NTRK3", "BRAF"]
	# ga中展示的基因在lc10中不重复展示
	LC10_snvindel = ["ALK", "EGFR", "ERBB2", "MET", "NRAS", "RET", "ROS1"]
	LC10_sv = ["ALK", "RET", "ROS1", "MET"]
	LC10_cnv = ["MET"]
	lc10_var = []
	ga_var = []
	other_var = []
	for var in var_list:
		if (var["bio_category"] == "Snvindel" and var["gene_symbol"] in LC10_snvindel) or \
		   (var["bio_category"] == "Sv" and set(re.split(",", var["gene_symbol"])) & set(LC10_sv)) or \
		   (var["bio_category"] == "Cnv" and var["gene_symbol"] in LC10_cnv):
			lc10_var.append(var)
		elif (var["bio_category"] == "Snvindel" and var["gene_symbol"] in ga_snvindel) or \
		   (var["bio_category"] == "Sv" and set(re.split(",", var["gene_symbol"])) & set(ga_sv)):
			ga_var.append(var)
		else:
			other_var.append(var)
	if var_type == "lc10":
		return lc10_var
	elif var_type == "ga":
		return ga_var
	else:
		return other_var
jinja2.filters.FILTERS["jnfy_ga_filter_lc10"] = jnfy_ga_filter_lc10

# 华西HRD-变异汇总-2025.07.03
def schx_hrd_var_sum(var_list):
	result = []
	for var in var_list:
		if var["hgvs_p_ZJZL"] != "p.?":
			result.append("{0} {1}".format(var["gene_symbol"], var["hgvs_p_ZJZL"]))
		else:
			result.append("{0} {1}".format(var["gene_symbol"], var["hgvs_c"]))
	return "，".join(result)
jinja2.filters.FILTERS["schx_hrd_var_sum"] = schx_hrd_var_sum

# 武汉协和BPTM Plus血液-所有基因检测结果汇总-2025.07.07
def whxh_gbptm_plus_summary(var_list):
	gene_list = ["POLE", "TP53", "BRCA1", "BRCA2", "CTNNB1", "EPCAM", "MLH1", "MSH2", "MSH6", "PMS2"]
	detect_gene_list = [var["gene_symbol"] for var in var_list]
	for gene in gene_list:
		if gene not in detect_gene_list:
			var_list.append({
				"gene_symbol" : gene
			})
	return sorted(var_list, key = lambda i:gene_list.index(i["gene_symbol"]))
jinja2.filters.FILTERS["whxh_gbptm_plus_summary"] = whxh_gbptm_plus_summary

# 江门中心150-乳腺癌需要区分乳腺癌常见基因和其他基因-2025.07.07
def jmzx_150_summary(info):
	var_list = info[0]
	tumor_type = info[1]
	breast_gene = ["BRCA1", "BRCA2", "CHEK2", "PTEN", "PALB2", "TP53", "RAD51C", "RAD51D"]
	breast_var = [var for var in var_list if var["gene_symbol"] in breast_gene]
	other_var = [var for var in var_list if var["gene_symbol"] not in breast_gene]
	if tumor_type == "intumor":
		return breast_var
	else:
		return other_var
jinja2.filters.FILTERS["jmzx_150_summary"] = jmzx_150_summary

def cds_region(str):
	#安徽省立 Sv five_cds_region, three_cds_region
	#EML4:NM_019063.5:intron6:c.668-6137
	str = str.split(":")
	return ":".join([str[1], str[3]])
jinja2.filters.FILTERS["cds_region"] = cds_region

# 重庆西南116-III类变异过滤掉同义、内含子、UTR和flankingregion-2025.07.11
def cqxn_filter_III_var(var_list):
	filter_type = ["Intronic", "Synonymous_Substitution", "5'UTR", "3'UTR", "FlankingRegion5", "FlankingRegion3"]
	result = []
	for var in var_list:
		if var["bio_category"] == "Snvindel":
			if var["type"] and var["type"] not in filter_type:
				result.append(var)
		else:
			result.append(var)
	return result
jinja2.filters.FILTERS["cqxn_filter_III_var"] = cqxn_filter_III_var

# 重庆西南116-获取变异解读中的参考文献-2025.07.11
def cqxn_116_getpmidfrominter(inter):
	pmid_list = get_pmid_from_inter(inter)
	result = []
	for pmid in pmid_list:
		if pmid not in result:
			result.append(pmid)
	return result
jinja2.filters.FILTERS["cqxn_116_getpmidfrominter"] = cqxn_116_getpmidfrominter

# 复旦中山CP200 EML4-ALK融合需要展示具体的亚型-2025.07.14
# 变异解读
def fdzs_cp200_alksv_inter(merge_sv_list):
	return shfk_cp_alksv(merge_sv_list)
jinja2.filters.FILTERS["fdzs_cp200_alksv_inter"] = fdzs_cp200_alksv_inter

# 复旦中山CP200 EML4-ALK融合需要展示具体的亚型-2025.07.14
# 小结展示
# 2025.09.01-“未知”改为“其他”
def fdzs_cp200_alksv_sum(sv_info):
	# gene1:exon1-gene2:exon2
	alk_sv_dict = {
		("exon13", "exon20") : "v1",
		("exon20", "exon20") : "v2",
		("exon6", "exon20") : "v3a/b",
		("exon15", "exon20") : "v4'",
		("exon2", "exon20") : "v5a/b",
		("exon18", "exon20") : "v5'",
		("exon14", "exon20") : "v7",
		("exon17", "exon20") : "v8a/b"
	}
	five_prime_cds = re.split("-", re.split(":", sv_info)[1])[0]
	three_prime_cds = re.split("-", re.split(":", sv_info)[-1])[0]
	# 2025.09.01-“未知”改为“其他”
	#return alk_sv_dict.get((five_prime_cds, three_prime_cds), "未知")
	return alk_sv_dict.get((five_prime_cds, three_prime_cds), "其他")
	# 2025.09.01-完成
jinja2.filters.FILTERS["fdzs_cp200_alksv_sum"] = fdzs_cp200_alksv_sum

# 2025.07.14
# 共突变，snvindel freq >= 0.01 且freq <= 0.05, 
# cnv ERBB2/MET cn_mean >= 3.5 且 cn_mean <= 5， cnv其他cn_mean >= 6且cn_mean <= 8添加上角标
# 2025.07.14-新增-若snvindel freq < 0.01，也需要加角标和备注
def zsly_co_mutation_add_subscript_v2(var_list):
	result = []
	for var in var_list:
		if "bio_category" in var.keys() and var["bio_category"] and var["bio_category"] == "Snvindel":
			if float(var["freq"]) >= 0.01 and float(var["freq"]) <= 0.05:
				if "Snvindel" not in result:
					result.append("Snvindel")
			elif float(var["freq"]) < 0.01:
				if "Svnidel_2" not in result:
					result.append("Snvindel_2")
		elif "bio_category" in var.keys() and var["bio_category"] and var["bio_category"] == "Cnv":
			if var["gene_symbol"] in ["ERBB2", "MET"]:
				if float(var["cn_mean"]) >= 3.5 and float(var["cn_mean"]) <= 5:
					if "Cnv" not in result:
						result.append("Cnv")
			else:
				if float(var["cn_mean"]) >= 6 and float(var["cn_mean"]) <= 8:
					if "Cnv" not in result:
						result.append("Cnv")
	return result
jinja2.filters.FILTERS["zsly_co_mutation_add_subscript_v2"] = zsly_co_mutation_add_subscript_v2

# 2025.07.17-重庆西南116-区分10基因获批位点和其他位点
def cqxn_10gene_appr_var_v2(info):
	var_list = info[0]
	gene_type = info[1]
	gene_list = ["ALK", "BRAF", "EGFR", "ERBB2", "KRAS", "MET", "NRAS", "PIK3CA", "RET", "ROS1"]
	appr_dict = {
		"EGFR" : ["p.G719A","p.E746_A750del","p.L747_T751del","p.V769_D770insASV","p.V774_C775insHV",\
				  "p.H773_V774insAH","p.H773_V774insY","p.P772_H773insPNP","p.V769_D770insGG","p.N771_P772insRHN",\
				  "p.L861Q","p.G719S","p.L747_P753delinsS","p.T790M","p.D770_N771insSVD",\
				  "p.D770delinsGY","p.N771_P772insN","p.H773_V774insTH","p.N771_P772insG","p.V769_D770insGTL",\
				  "p.P772_H773insGHP","p.G719C","p.E746_S752delinsV","p.S768I","p.A763_Y764insFQEA",\
				  "p.H773_V774insH","p.N771_P772insH","p.N771delinsGY","p.N771delinsGF","p.D770_N771insGF",\
				  "p.H773delinsYNPY","p.E746_A750del","p.L747_A750delinsP","p.D770_N771insG","p.H773_V774insNPH",\
				  "p.H773_V774insPH","p.P772_H773insYNP","p.P772_H773insQ","p.P772_H773insGNP","p.D770_N771insSLA","p.L858R"],
		"KRAS" : ["p.G12D", "p.G12A", "p.G12V", "p.G12S", "p.G12C", "p.Q61H"],
		"NRAS" : ["p.G12D", "p.Q61R", "p.Q61K"],
		"PIK3CA" : ["p.H1047R"],
		"BRAF" : ["p.V600E"],
		"ERBB2" : ["p.A775_G776insYVMA"],
		"MET" : ["c.3082+1G>T"],
		"ALK" : ["EML4:exon13-ALK:exon20", "EML4:exon6-ALK:exon20", "EML4:exon20-ALK:exon20"],
		"ROS1" : ["CD74:exon6-ROS1:exon34", "GOPC:exon8-ROS1:exon35"],
		"RET" : ["KIF5B:exon15-RET:exon12"]
	}
	appr_var = []
	other_var = []
	for var in var_list:
		if var["bio_category"] == "Snvindel":
			snv_info = var["hgvs_p_ZJZL"] if var["hgvs_p_ZJZL"] != "p.?" else var["hgvs_c"]
			if var["gene_symbol"] in appr_dict.keys() and snv_info in appr_dict[var["gene_symbol"]]:
				appr_var.append(var)
			else:
				other_var.append(var)
		elif var["bio_category"] == "Sv":
			sv_info = "{0}:{1}-{2}:{3}".format(var["five_prime_gene"], var["five_prime_cds"], var["three_prime_gene"], var["three_prime_cds"])
			# ALK/ROS1/RET几个获批变异型，主基因都在3端
			if var["three_prime_gene"] in appr_dict.keys() and sv_info in appr_dict[var["three_prime_gene"]]:
				appr_var.append(var)
			else:
				other_var.append(var)
		else:
			other_var.append(var)
	detect_gene = []
	for var in appr_var:
		for gene in re.split(",", var["gene_symbol"]):
			if gene in gene_list and gene not in detect_gene:
				detect_gene.append(gene)
	result = []
	result.extend(appr_var)
	for gene in sorted(list(set(gene_list) - set(detect_gene))):
		result.append({"gene_symbol" : gene})
	if gene_type == "LC10":
		return result
	else:
		return other_var
jinja2.filters.FILTERS["cqxn_10gene_appr_var_v2"] = cqxn_10gene_appr_var_v2

# 适用北大三CP200-2025.07.22
def io_detect_for_BDS_CP200(info):
	var_list = info[0]
	return_type = info[1]
	io_result = {}
	# 汇总体细胞I/II/肿瘤发生发展相关变异+胚系致病/疑似致病性变异
	io_gene_P = ["ATM","ATR","BRCA1","BRCA2","BRIP1","CHEK1","CHEK2","FANCA","MRE11",\
				 "PALB2","RAD50","MLH1","MSH2","MSH6","PMS2","POLE","POLD1","TP53",\
				 "KRAS","CD274","ARID1A","SETD2","TERT","KMT2D","FAT1","CDK12"]
	io_gene_N = ["EGFR","ALK","MDM2","MDM4","CDKN2A","CDKN2B","DNMT3A","STK11","IFNGR1",\
				 "JAK1","JAK2","APC","CTNNB1","B2M","PTEN","FGF19"]
	cnv_gene_list = ["CD274", "MDM2", "MDM4", "FGF19"]
	
	for var in var_list:
		# 仅展示扩增的基因
		if var["bio_category"] == "Cnv" and var["gene_symbol"] in cnv_gene_list:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			io_result[var["gene_symbol"]].append(var)
		# 仅展示融合的基因
		elif var["bio_category"] in ["Sv", "PSeqRnaSv"] and set(re.split(",", var["gene_symbol"])) & set(["ALK"]):
			if "ALK" not in io_result.keys():
				io_result.setdefault("ALK", [])
			io_result["ALK"].append(var)
		# 其余基因展示Snvindel
		elif var["bio_category"] == "Snvindel" and var["gene_symbol"] in io_gene_P + io_gene_N:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			io_result[var["gene_symbol"]].append(var)

	# summary展示
	io_p_list = [i for k,v in io_result.items() for i in v if k in io_gene_P]
	io_n_list = [i for k,v in io_result.items() for i in v if k in io_gene_N]
	
	if return_type == "p":
		return io_p_list
	elif return_type == "n":
		return io_n_list
jinja2.filters.FILTERS["io_detect_for_BDS_CP200"] = io_detect_for_BDS_CP200

# 国际部胚系项目-2025.07.23
# 1. 通过基因、变异分组和性别匹配遗传肿瘤
# 2. 相同基因相同配置信息的放在一个表格里展示
def international_risk(info):
	var_list = info[0]
	config = info[1]
	sample = info[2]
	for var in var_list:
		category_risk = []
		var_category_names = var["var_category_names"] if var["var_category_names"] else ""
		for category in re.split(",", var_category_names):
			if var["gene_symbol"] in config.keys() and category in config[var["gene_symbol"]].keys():
				category_risk.extend(config[var["gene_symbol"]][category])
		
		var["var_tumor_risk"] = []
		var["var_disease_risk"] = []
		for cate_item in category_risk:
			tumor_str_ = "{0} ({1})".format(cate_item["tumor"], cate_item["tumor_hereditary_pattern"]) if cate_item["tumor_hereditary_pattern"] != "-" else cate_item["tumor"]+" (-)"
			disease_str_ = "{0} ({1})".format(cate_item["disease"], cate_item["disease_hereditary_pattern"]) if cate_item["disease_hereditary_pattern"] != "-" else cate_item["disease"]+" (-)"
			if cate_item["gender"] == "-" or cate_item["gender"].lower() == sample["gender"].lower():
				if tumor_str_ not in var["var_tumor_risk"]:
					var["var_tumor_risk"].append(tumor_str_)
			if disease_str_ not in var["var_disease_risk"]:
				var["var_disease_risk"].append(disease_str_)
	
	tmp_dict = {}
	for var in var_list:
		key = (var["gene_symbol"], tuple(var["var_tumor_risk"]), tuple(var["var_disease_risk"]))
		if key not in tmp_dict.keys():
			tmp_dict.setdefault(key, [])
		tmp_dict[key].append(var)
	
	result = []
	for key, var_l in tmp_dict.items():
		result.append(
			{
				"var_list" : var_l,
				"var_tumor_risk" : var_l[0]["var_tumor_risk"],
				"var_disease_risk" : var_l[0]["var_disease_risk"]
			}
		)
	return result
jinja2.filters.FILTERS["international_risk"] = international_risk

# 北大三-CP200-筛选出不合格的质控项-2025.07.25
def bds_cp200_qc(info):
	lib_qc = info[0]
	ngs_qc = info[1]
	ngs_qc_standard = {
		"cleandata_q30_num" : 0.75,
		"depth_ssbc_num" : 400,
		"depth_rna_ctrl_num" : 10
	}

	lib_qc_standard = {
		"library_concn" : 5	
	}

	ngs_qc_key = {
		"cleandata_q30_num" : "DNA样本Q30",
		"depth_ssbc_num" : "平均有效深度(x)",
		"depth_rna_ctrl_num" : "RNA内参绝对拷贝数(x)"
	}

	lib_qc_key = {
		"library_concn" : "文库DNA浓度(ng/μL)"	
	}

	ngs_dna_qc = ngs_qc["dna_data_qc"] if "dna_data_qc" in ngs_qc.keys() and ngs_qc["dna_data_qc"] else {}
	lib_dna_qc = lib_qc["lib_dna_qc"] if "lib_dna_qc" in lib_qc.keys() and lib_qc["lib_dna_qc"] else {}

	Fail_item = []
	# DNA 湿实验质控
	for item in lib_qc_standard.keys():
		if item in lib_dna_qc.keys() and lib_dna_qc[item] and is_number(lib_dna_qc[item]) and float(lib_dna_qc[item]) < lib_qc_standard.get(item):
			Fail_item.append(lib_qc_key.get(item))
	# DNA NGS质控
	for item in ngs_qc_standard.keys():
		if item in ngs_dna_qc.keys() and ngs_dna_qc[item] and is_number(ngs_dna_qc[item]) and float(ngs_dna_qc[item]) < ngs_qc_standard.get(item):
			Fail_item.append(ngs_qc_key.get(item))
		# 考虑质控结果为0的情况
		elif item in ngs_dna_qc.keys() and not ngs_dna_qc[item]:
			Fail_item.append(ngs_qc_key.get(item))
	return "，".join(Fail_item)
jinja2.filters.FILTERS["bds_cp200_qc"] = bds_cp200_qc

# 国际部胚系项目-2025.07.23
# 1. 通过基因、变异分组和性别匹配遗传肿瘤
# 2. 相同基因相同配置信息的放在一个表格里展示
# 3. 疾病不同的分多行展示
def international_risk_v2(info):
	var_list = info[0]
	config = info[1]
	sample = info[2]
	for var in var_list:
		category_risk = []
		var_category_names = var["var_category_names"] if var["var_category_names"] else ""
		for category in re.split(",", var_category_names):
			if var["gene_symbol"] in config.keys() and category in config[var["gene_symbol"]].keys():
				category_risk.extend(config[var["gene_symbol"]][category])
		
		# 筛选出无性别区分或性别相同的信息
		filter_category_risk = [i for i in category_risk if i["gender"] == "-" or i["gender"].lower() == sample["gender"].lower()]
		disease_dict = {}
		for cate_item in filter_category_risk:
			tumor_str_ = "{0} ({1})".format(cate_item["tumor"], cate_item["tumor_hereditary_pattern"]) if cate_item["tumor_hereditary_pattern"] != "-" else cate_item["tumor"]+" (-)"
			disease_str_ = "{0} ({1})".format(cate_item["disease"], cate_item["disease_hereditary_pattern"]) if cate_item["disease_hereditary_pattern"] != "-" else cate_item["disease"]+" (-)"
			if disease_str_ not in disease_dict.keys():
				disease_dict.setdefault(disease_str_, [])
			disease_dict[disease_str_].append(tumor_str_)
		var["disease_risk"] = disease_dict
	
	tmp_dict = {}
	for var in var_list:
		key = (var["gene_symbol"], tuple(var["disease_risk"]))
		if key not in tmp_dict.keys():
			tmp_dict.setdefault(key, [])
		tmp_dict[key].append(var)
	
	result = []
	for key, var_l in tmp_dict.items():
		result.append(
			{
				"var_list" : var_l,
				"disease_risk" : [{"disease" : k, "tumor" : v} for k,v in var_l[0]["disease_risk"].items()]
			}
		)
	return result
jinja2.filters.FILTERS["international_risk_v2"] = international_risk_v2

# 郑大一gBRCA小结-2025.07.27
def ZDY_gbrca_germline_summary(var_list):
	var_dict = {}
	for var in var_list:
		if var["type"] == "Loss":
			var_info = "{0} {1} del".format(var["gene_symbol"], var["value"])
		elif var["type"] == "Gain":
			var_info = "{0} {1} dup".format(var["gene_symbol"], var["value"])
		else:
			if var["hgvs_p"] != "p.?":
				var_info = "{0} {1}:{2}:{3}".format(var["gene_symbol"], var["transcript_primary"], var["hgvs_c"], var["hgvs_p"])
			else:
				var_info = "{0} {1}:{2}".format(var["gene_symbol"], var["transcript_primary"], var["hgvs_c"])
		if var["clinic_num_g"] not in var_dict.keys():
			var_dict.setdefault(var["clinic_num_g"], [])
		var_dict[var["clinic_num_g"]].append(var_info)
	result = []
	if 5 in var_dict.keys():
		result.append("检出{0}个致病性变异，{1}".format(str(len(var_dict[5])), ", ".join(var_dict[5])))
	if 4 in var_dict.keys():
		result.append("检出{0}个疑似致病性变异，{1}".format(str(len(var_dict[4])), ", ".join(var_dict[4])))
	if not var_dict:
		result.append("未检出致病或疑似致病性胚系变异")
	return "；".join(result)
jinja2.filters.FILTERS["ZDY_gbrca_germline_summary"] = ZDY_gbrca_germline_summary

# 浙二HRR，变异总结-2025.07.27
def zjey_hrr_varsum(var_list):
	result = []
	for var in var_list:
		if var["bio_category"] == "Snvindel":
			hgvs = var["hgvs_p"] if var["hgvs_p"] != "p.?" else var["hgvs_c"]
			freq_ = ""
			if var["var_origin"] == "germline":
				freq_ = "纯合" if float(var["freq"]) >= 0.85 else "杂合"
			else:
				freq_ = var["freq_str"]
			result.append("{0} {1}（{2}）".format(var["gene_symbol"], hgvs, freq_))
		elif var["type"] == "Loss":
			hete_or_homo = ""
			if "cnv_type" in var.keys() and var["cnv_type"]:
				hete_or_homo = "杂合" if var["cnv_type"] == "HeteDel" else "纯合" if var["cnv_type"] == "HomoDel" else ""
			if hete_or_homo:
				result.append("{0} {1} del（{2}）".format(var["gene_symbol"], var["value"], hete_or_homo))
			else:
				result.append("{0} {1} del".format(var["gene_symbol"], var["value"]))
		elif var["type"] == "Gain":
			result.append("{0} {1} dup".format(var["gene_symbol"], var["value"]))
	return "、".join(result)
jinja2.filters.FILTERS["zjey_hrr_varsum"] = zjey_hrr_varsum


# 获取描述中的PMID号-v2-2025.07.30
# 示例：……（PMID：111, PMID: 222, PMID: 333）……
def getPMID_from_inter_v2(inter):
	pmid_list = []
	mat = re.compile(r"\（PMID.*?\）")
	mat2 = re.compile(r"\(PMID.*?\)")
	for i in mat.findall(str(inter)) + mat2.findall(str(inter)):
		replace_blank = [" ", "PMID", "（", "）", ":", "：", "(", ")"]
		for str_ in replace_blank:
			i = i.replace(str_, "")
		i = i.replace("，", ",")
		for r in re.split(",", i):
			if r not in pmid_list:
				pmid_list.append(r)
	return pmid_list

# 中山六院150参考文献-2025.07.30
# 基因介绍、变异解读获取描述里的PMID号
# 遗传风险获取literature_evi_sum
def szly_gene150_refer(info):
	var_list = info[0]
	refer_original = info[1]
	# refer_type分为pmid 和 other
	refer_type = info[2]
	# 参考文献code对应内容
	refer_code_dict = {}
	for refer in refer_original:
		refer["pmid"] = (re.split("PMID:", refer["pmid"])[-1]).replace("]", "") if refer["pmid"] else ""
		refer["title"] = refer["title"] if refer["title"] else ""
		if "refer_code" in refer.keys() and refer["refer_code"]:
			refer_code_dict[refer["refer_code"]] = refer
	
	# 获取变异中的参考文献
	pmid_refer = []
	other_refer = []
	for var in var_list:
		# 获取基因介绍pmid
		for pmid in getPMID_from_inter_v2(var["gene_function"] if var["gene_function"] else ""):
			if pmid not in pmid_refer:
				pmid_refer.append(pmid)
		# 获取变异解读pmid
		for pmid in getPMID_from_inter_v2(var["variant_interpret_cn"] if var["variant_interpret_cn"] else ""):
			if pmid not in pmid_refer:
				pmid_refer.append(pmid)
		# 获取遗传风险证据pmid和其他参考文献
		literature_evi_sum_list = []
		if "Predisposing" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Predisposing"]:
			for evi in var["evi_sum"]["evi_split"]["Predisposing"]:
				if "literature_evi_sum" in evi.keys() and evi["literature_evi_sum"]:
					literature_evi_sum_list.extend(evi["literature_evi_sum"])
		for refer_code in literature_evi_sum_list:
			if refer_code in refer_code_dict.keys():
				refer_info = refer_code_dict.get(refer_code)
				if refer_info["pmid"] and refer_info["pmid"] not in pmid_refer:
					pmid_refer.append(refer_info["pmid"])
				elif refer_info["title"] and refer_info["title"] not in other_refer:
					other_refer.append(refer_info["title"])
	#  针对不同表格，返回PMID或非PMID参考文献
	if refer_type == "pmid":
		return zsly_pmid(pmid_refer)
	else:
		other_refer_list = []
		num = len(pmid_refer) + 1
		for i in other_refer:
			other_refer_list.append({
				"num" : num,
				"title" : i
			})
			num += 1
		return other_refer_list
jinja2.filters.FILTERS["szly_gene150_refer"] = szly_gene150_refer

# 福建附一CP200-2025.08.02
# 结果小结处结直肠癌/子宫内膜癌需要展示林奇4个基因检测情况
# 频率>=40% <=60%，致病性、疑似致病性
def fjfy_cp200_lyn_result(var_list):
	gene_list = ["MLH1", "PMS2", "MSH2", "MSH6"]
	lyn_var_list = [var for var in var_list if var["bio_category"] == "Snvindel" and \
					var["gene_symbol"] in gene_list and var["clinic_num_g"] in [4,5] and \
					float(var["freq"]) >= 0.4 and float(var["freq"]) <= 0.6]
	for var in lyn_var_list:
		var["gene_region_cn"] = gene_region_strn(var["gene_region"])

	return lyn_var_list
jinja2.filters.FILTERS["fjfy_cp200_lyn_result"] = fjfy_cp200_lyn_result

# 福建附一CP200-2025.08.02
# 结直肠癌/子宫内膜癌判断是否为胚系相关变异
def fjfy_cp200_judge_germline_gene(var):
	germline_gene = ["APC", "BRCA1", "BRCA2", "EPCAM", "FH", "MLH1", "MSH2", "MSH6", \
					 "NF1", "NF2", "PMS2", "RAD50", "RAD51B", "RAD51C", "RAD51D", "RAD54L", \
					 "STK11", "TP53", "TSC1", "TSC2", "VHL"]
	lyn_gene = ["MLH1", "PMS2", "MSH2", "MSH6"]
	result = ""
	if var["bio_category"] == "Snvindel" and float(var["freq"]) >= 0.4 and float(var["freq"]) <= 0.6 \
		and var["clinic_num_g"] in [4,5]:
		if var["gene_symbol"] in lyn_gene:
			result = "lyn"
		elif var["gene_symbol"] in germline_gene:
			result = "germline"
	return result
jinja2.filters.FILTERS["fjfy_cp200_judge_germline_gene"] = fjfy_cp200_judge_germline_gene

# 福建附一CP200-2025.08.02
# 结直肠癌/子宫内膜癌过滤指定变异列表中可能为胚系的变异
def fjfy_filter_germline(info):
	var_list = info[0]
	var_type = info[1]
	if var_type == "ec":
		return [var for var in var_list if fjfy_cp200_judge_germline_gene(var) in ["lyn", "germline"]]
	else:
		return [var for var in var_list if fjfy_cp200_judge_germline_gene(var) not in ["lyn", "germline"]]
jinja2.filters.FILTERS["fjfy_filter_germline"] = fjfy_filter_germline

# 福建附一CP200-2025.08.02
# 肺癌判断是否为EGFR L858R/19del/T790M、ALK融合
def fjfy_cp200_judge_egfr_alk(var):
	result = False
	if var["bio_category"] == "Sv":
		if "ALK" in re.split(",", var["gene_symbol"]):
			result = True
	elif var["bio_category"] == "Snvindel" and var["gene_symbol"] == "EGFR":
		if var["hgvs_p"] in ["p.(L858R)", "p.(T790M)"]:
			result = True
		elif "EGFR Exon19 del" in var["var_category_names"]:
			result = True
	return result
jinja2.filters.FILTERS["fjfy_cp200_judge_egfr_alk"] = fjfy_cp200_judge_egfr_alk

# 福建附一CP200-2025.08.02
# 筛选出有适应症的证据
def fjfy_filter_adaptation(evi_split):
	result = []
	for evi_type in ["Predictive_merge", "Prognostic", "Diagnostic"]:
		if evi_type in evi_split.keys() and evi_split[evi_type]:
			for evi in evi_split[evi_type]:
				print (evi["evi_adaptation_disease_cn"])
			result.extend([evi for evi in evi_split[evi_type] if evi["evi_adaptation_disease_cn"]])
	return result
jinja2.filters.FILTERS["fjfy_filter_adaptation"] = fjfy_filter_adaptation

# 北大三，CNV/PHd需要单独表格展示-2025.08.02
def bds_mp_filter_cnv(info):
	var_list = info[0]
	var_type = info[1]
	if var_type == "cnv":
		return [var for var in var_list if var["bio_category"] in ["Cnv", "PHd"]]
	else:
		return [var for var in var_list if var["bio_category"] not in ["Cnv", "PHd"]]
jinja2.filters.FILTERS["bds_mp_filter_cnv"] = bds_mp_filter_cnv

# 解放军总医院第七医学中心-MP组织-荧光定量PCR位点检出状态-少红-2025.08.05
def jfjzy_s_pcr_trans(info):
	# EGFR
	egfr_pcr_dict = {"E746_A750del (1)": "c.2235_2249del",
					 "E746_A750del (2)": "c.2236_2250del",
					 "E746_A750del(1)": "c.2235_2249del",
					 "E746_A750del(2)": "c.2236_2250del",
					 "L747_P753>S": "L747_P753delinsS",
					 "E746_T751>I": "E746_T751delinsI",
					 "E746_S752>V": "E746_S752delinsV",
					 "L747_T751>Q": "L747_T751delinsQ",
					 "L747_A750>P": "L747_A750delinsP",
					 "L747_P753>Q": "L747_P753delinsQ",
					 "L747_T751>P": "L747_T751delinsP"}
	var_dict = info[0]
	gene = info[1]
	exon = info[2]
	pcr_p_list = [egfr_pcr_dict[i] if i in egfr_pcr_dict else i for i in info[3].replace(' ', '').split(',')]
	
	var_list = var_dict["var_somatic"]["level_I"] + \
			   var_dict["var_somatic"]["level_II"] + \
			   var_dict["var_somatic"]["level_onco_nodrug"] + \
			   var_dict["var_somatic"]["level_III"] + \
			   var_dict["var_germline"]["level_5"] + \
			   var_dict["var_germline"]["level_4"] + \
			   var_dict["var_germline"]["level_3"] 
	for var in var_list:
		if var["bio_category"] == "Snvindel":
			if var["gene_symbol"] == gene and var["gene_region"] == exon:
				hgvs_p = var["hgvs_p"].replace("p.","").replace("(","").replace(")","")
				hgvs_c = var["hgvs_c"]
				if hgvs_c in pcr_p_list or hgvs_p in pcr_p_list:
					return "突变型"
	return "未见变异"
jinja2.filters.FILTERS["jfjzy_s_pcr_trans"] = jfjzy_s_pcr_trans

# 解放军总医院第七医学中心-MP组织-位点筛选-少红-2025.08.05
def jfjzy_s_var_filter(info):
	var_dict = info[0]
	gene = info[1]
	
	result = []
	# 20250805-预测胚系3类不展示（BRCA基因3类无用药，所以这边不考虑3类有药的情况了-刘炜芬）
	# 20250806-客户确认后又要加预测胚系3类变异-刘炜芬
	var_list = var_dict["var_somatic"]["level_I"] + \
			   var_dict["var_somatic"]["level_II"] + \
			   var_dict["var_somatic"]["level_onco_nodrug"] + \
			   var_dict["var_somatic"]["level_III"] + \
			   var_dict["var_germline"]["level_5"] + \
			   var_dict["var_germline"]["level_4"] + \
			   var_dict["var_germline"]["level_3"]
	
	for var in var_list:
		if var["gene_symbol"] == gene and var["bio_category"] == "Snvindel":
			result.append(var)
	return result
jinja2.filters.FILTERS["jfjzy_s_var_filter"] = jfjzy_s_var_filter

# 1. 治疗方案介绍-生物标志物-v3-2025.08.06-取消超文本-适用吉大一MP，氨基酸变化用hgvs_p_ZJZL
def approval_regimen_biomarker_v3_jdyy(info):
	biomaker_list = info[0]
	judge_mergeMET = info[1]
	result = []
	for i in biomaker_list:
		if "hgvs_c" in i.keys() and i["hgvs_c"]:
			if i["hgvs_p"] != "p.?":
				hgvs_p = i["hgvs_p"].replace("(", "").replace(")", "") if "=" not in i["hgvs_p"] else i["hgvs_p"]
				result.append("{0} {1} {2} {3}".format(i["gene_symbol"], i["gene_region"], i["hgvs_c"], hgvs_p))
			else:
				result.append("{0} {1} {2}".format(i["gene_symbol"], i["gene_region"], i["hgvs_c"]))
		elif "cnv_type" in i.keys() and i["cnv_type"]:
			# 2024.08.30-CNV 区分Loss，其他的写扩增
			if i["cnv_type"] == "Loss" and i["gene_symbol"] in ["BRCA1", "BRCA2"]:
				result.append("{0} 大片段缺失".format(i["gene_symbol"]))
			elif i["cnv_type"] == "Gain" and i["gene_symbol"] in ["BRCA1", "BRCA2"]:
				result.append("{0} 大片段重复".format(i["gene_symbol"]))
			elif i["cnv_type"] == "HeteDel" and i["gene_symbol"] in ["BRCA1", "BRCA2"]:
				result.append("{0} 杂合大片段缺失".format(i["gene_symbol"]))
			elif i["cnv_type"] == "HomoDel" and i["gene_symbol"] in ["BRCA1", "BRCA2"]:
				result.append("{0} 纯合大片段缺失".format(i["gene_symbol"]))
			else:
				result.append("{0} 扩增".format(i["gene_symbol"]))
		elif "five_prime_gene" in i.keys() and i["five_prime_gene"]:
			if i["five_prime_gene"] == "MET" and i["three_prime_gene"] == "MET":
				result.append("MET exon14 跳跃")
			else:
				# 重新拆分hgvs，CP40的region少了ins
				if "hgvs" in i.keys() and i["hgvs"]:
					i["five_prime_region"] = "-".join(re.split("-", (re.split(":", i["hgvs"])[2]))[:-1]) \
						  					 if not re.search("--", i["hgvs"]) \
											 else re.split("_", (re.split("--", i["hgvs"])[0]))[-1]
					i["three_prime_region"] = re.split(":", i["hgvs"])[-1] \
											  if not re.search("--", i["hgvs"]) \
											  else re.split("_", (re.split("--", i["hgvs"])[1]))[-1]
					# 加一个兼容-2023.10.19
					# var_hgvs新格式，gene1:NM_xxx:exon1--gene2:NM_xxx:exon2, 旧的为gene1:NM_xxx_exon1--gene2:NM_xxx_exon2
					# cds会变成xxx:exon1和xxx:exon2
					i["five_prime_region"] = re.split(":", i["five_prime_region"])[-1] if re.search(":", i["five_prime_region"]) else i["five_prime_region"]
					i["three_prime_region"] = re.split(":", i["three_prime_region"])[-1] if re.search(":", i["three_prime_region"]) else i["three_prime_region"]
					# 兼容完成-2023.10.19
					# 加一个兼容-2024.01.25
					# 4. gene:转录本_exon-gene:转录本_exon 重新提取后，five_prime_cds为空，以此做为重新拆分的判定依据
					if not i["five_prime_region"]:
						i["five_prime_region"] = re.split("_", re.split("-", i["hgvs"])[0])[-1]
						i["three_prime_region"] = re.split("_", i["hgvs"])[-1]
					# 兼容完成-2024.01.25
				result.append("{0}:{1}:{2}-{3}:{4}:{5}".format(i["five_prime_gene"], i["five_prime_transcript"], i["five_prime_region"],\
															   i["three_prime_gene"], i["three_prime_transcript"], i["three_prime_region"]))
		elif "biomarker_type" in i.keys() and i["biomarker_type"]:
			if i["biomarker_type"] == "KRAS/NRAS/BRAF WT":
				result.append("KRAS/NRAS/BRAF 野生型")
			elif i["biomarker_type"] == "HRD-":
				result.append("HRD阴性")
			elif i["biomarker_type"] == "HRD+":
				result.append("HRD阳性")
			# 2025.07.17-新增HD
			elif "HomoDel" in i["biomarker_type"] or "HeteDel" in i["biomarker_type"]:
				split_str = re.split(":", i["biomarker_type"])
				if len(split_str) == 4:
					hd_type = "纯合缺失" if split_str[-1] == "HomoDel" else "杂合缺失" if split_str[-1] == "HeteDel" else "未知变异类型！"
					result.append(split_str[1] + " " + hd_type + " " + split_str[2])
					#result.append(" ".join(split_str[1:4]))
				else:
					result.append(i["biomarker_type"])
			else:
				result.append(i["biomarker_type"])
		else:
			result.append("无法分辨的分子标志物！")
	# 若存在MET 14跳跃DNA/RNA共检的话，则删除RNA里的结果，仅保留DNA
	result_redup = []
	for i in result:
		if i not in result_redup:
			result_redup.append(i)
	# 2024.02.27更新-MET DNA/RNA共检时都要展示
	#if judge_mergeMET:
	#	if "MET exon14 跳跃" in result_redup:
	#		result_redup.remove("MET exon14 跳跃")
	# 2024.02.27更新完成
	if not result_redup:
		result_redup = ["-"]
	#rt = RichText()
	#rt.add("\n".join(result_redup))
	return "\n".join(result_redup)
jinja2.filters.FILTERS["approval_regimen_biomarker_v3_jdyy"] = approval_regimen_biomarker_v3_jdyy

# 10. 体细胞/来源不明变异-中国人民解放军总医院第七医学中心-（预测）胚系3类也要展示-2025.08.07
def var_s_summary_JFJZYDQ(info):
	var_list = info[0]
	sample = info[1]
	result = []
	def sum_var(vlist):
		v_result = []
		for var in vlist:
			if var["bio_category"] == "Snvindel":
				if var["hgvs_p"] != "p.?":
					v_result.append(var["gene_symbol"]+" "+var["hgvs_p"])
				else:
					v_result.append(var["gene_symbol"]+" "+var["hgvs_c"])
			elif var["bio_category"] == "Cnv":
				v_result.append(var["gene_symbol"]+" 扩增")
			elif var["bio_category"] in ["Sv", "PSeqRnaSv"]:
				if var["five_prime_gene"] == "MET" and var["three_prime_gene"] == "MET":
					v_result.append("MET exon14 跳跃")
				else:
					if var["five_prime_gene"]+"-"+var["three_prime_gene"]+" 融合" not in v_result:
						v_result.append(var["five_prime_gene"]+"-"+var["three_prime_gene"]+" 融合")
			# 2025.07.01-新增PHd
			elif var["bio_category"] == "PHd":
				if var["type"] == "HomoDel":
					if var["gene_symbol"]+" 纯合缺失" not in v_result:
						v_result.append(var["gene_symbol"]+" 纯合缺失")
				elif var["type"] == "HeteDel":
					if var["gene_symbol"]+" 杂合缺失" not in v_result:
						v_result.append(var["gene_symbol"]+" 杂合缺失")
				else:
					v_result.append(var["gene_symbol"]+" 未知变异类型！")
			# 2025.07.01-新增完成
		return ", ".join(v_result)
	# 有对照样本时
	c_var_all_num = len(var_list["var_somatic"]["level_I"] + var_list["var_somatic"]["level_II"] + \
				 var_list["var_somatic"]["level_onco_nodrug"] + var_list["var_somatic"]["level_III"])
	c_var_onco_drug_num = len(var_list["var_somatic"]["level_I"] + var_list["var_somatic"]["level_II"])
	c_var_onco_nodrug_num = len(var_list["var_somatic"]["level_onco_nodrug"])
	c_var_onco_drug_str = sum_var(var_list["var_somatic"]["level_I"] + var_list["var_somatic"]["level_II"])
	# 无对照样本时
	# 2024.05.20-总变异数需要加上疑似胚系非致病但是有用药的位点
	#var_all_num = len(var_list["var_somatic"]["level_I"] + var_list["var_somatic"]["level_II"] + \
	#			 var_list["var_somatic"]["level_onco_nodrug"] + var_list["var_somatic"]["level_III"] +\
	#			 var_list["var_germline"]["level_5"] + var_list["var_germline"]["level_4"])
	var_all_num = len(var_list["var_somatic"]["level_I"] + var_list["var_somatic"]["level_II"] + \
				 var_list["var_somatic"]["level_onco_nodrug"] + var_list["var_somatic"]["level_III"] +\
				 var_list["var_germline"]["level_5"] + var_list["var_germline"]["level_4"] + var_list["var_germline"]["level_3"])
	var_onco_drug_num = len(var_list["var_somatic"]["level_I"] + var_list["var_somatic"]["level_II"]+\
							var_list["var_germline"]["regimen_level_I"] + var_list["var_germline"]["regimen_level_II"])
	var_onco_nodrug_num = len(var_list["var_somatic"]["level_onco_nodrug"] + var_list["var_germline_nodrug"])
	var_onco_drug_str = sum_var(var_list["var_somatic"]["level_I"] + var_list["var_somatic"]["level_II"]+\
								  var_list["var_germline"]["regimen_level_I"] + var_list["var_germline"]["regimen_level_II"])
	if sample["control_sample_id"]:
		str_1 = "检出{0}个体细胞变异，其中具有临床意义的变异有{1}个，肿瘤发生发展相关变异有{2}个。".format(str(c_var_all_num), str(c_var_onco_drug_num), str(c_var_onco_nodrug_num))
		str_2 = "具有临床意义的变异有{0}。".format(c_var_onco_drug_str)
		if c_var_all_num != 0:
			if c_var_onco_drug_num != 0:
				return str_1 + str_2
			else:
				return str_1
		else:
			return "在检测范围内，未检出体细胞变异。"
	else:
		str_1 = "检出{0}个变异，其中具有临床意义的变异有{1}个，肿瘤发生发展相关变异有{2}个。".format(str(var_all_num), str(var_onco_drug_num), str(var_onco_nodrug_num))
		str_2 = "具有临床意义的变异有{0}。".format(var_onco_drug_str)
		if var_all_num != 0:
			if var_onco_drug_num != 0:
				return str_1 + str_2
			else:
				return str_1
		else:
			return "在检测范围内，未检出变异。"
jinja2.filters.FILTERS["var_s_summary_JFJZYDQ"] = var_s_summary_JFJZYDQ

# 福建附一150-2025.08.07
# 林奇5个基因检测情况需要分开展示
# 致病性、疑似致病性
def fjfy_150_lyn_result(var_list):
	gene_list = ["MLH1", "PMS2", "MSH2", "MSH6", "EPCAM"]
	lyn_var_list = [var for var in var_list if var["gene_symbol"] in gene_list]
	for var in lyn_var_list:
		if var["bio_category"] == "Snvindel":
			var["gene_region_cn"] = gene_region_strn(var["gene_region"])
		else:
			var["gene_region_cn"] = var["value"].replace("exon", "") + "号外显子" if "value" in var.keys() and var["value"] else ""

	return lyn_var_list
jinja2.filters.FILTERS["fjfy_150_lyn_result"] = fjfy_150_lyn_result

# 福建附一150-2025.08.08
# 林奇5个基因3类+其他基因3/4/5类
def fjfy_150_other_result(var_list):
	result = []
	lyn_gene_list = ["MLH1", "PMS2", "MSH2", "MSH6", "EPCAM"]
	other_gene_list = ["APC", "BRCA1", "BRCA2", "FH", "NF1", "NF2", "RAD50", "RAD51B", \
					   "RAD51C", "RAD51D", "RAD54L", "STK11", "TP53", "TSC1", "TSC2", "VHL"]
	for var in var_list:
		if var["gene_symbol"] in lyn_gene_list and var["clinic_num_g"] == 3:
			result.append(var)
		elif var["gene_symbol"] in other_gene_list and var["clinic_num_g"] in [3, 4, 5]:
			result.append(var)
	return result
jinja2.filters.FILTERS["fjfy_150_other_result"] = fjfy_150_other_result

# 福建附一150-2025.08.08
# 林奇5个基因5/4类+其他基因5/4类
def fjfy_150_45_result(info):
	var_list = info[0]
	gene_type = info[1]
	level_list = info[2]
	gene_list = ["MLH1", "PMS2", "MSH2", "MSH6", "EPCAM", \
			  	 "APC", "BRCA1", "BRCA2", "FH", "NF1", "NF2", "RAD50", "RAD51B", \
				 "RAD51C", "RAD51D", "RAD54L", "STK11", "TP53", "TSC1", "TSC2", "VHL"]
	if gene_type == "in":
		return [var for var in var_list if var["gene_symbol"] in gene_list and var["clinic_num_g"] in level_list]
	else:
		return [var for var in var_list if var["gene_symbol"] not in gene_list and var["clinic_num_g"] in level_list]
jinja2.filters.FILTERS["fjfy_150_45_result"] = fjfy_150_45_result

# 甘肃武威CP200-变异统计时不包含10基因-2025.08.11
def gsww_cp200_filter_10(var_list):
	gene_list = ["MET", "EGFR", "ALK", "KRAS", "ROS1", "RET", "ERBB2", "BRAF", "NRAS", "PIK3CA"]
	return [var for var in var_list if not set(re.split(",", var["gene_symbol"])) & set(gene_list)]
jinja2.filters.FILTERS["gsww_cp200_filter_10"] = gsww_cp200_filter_10

# 甘肃武威MP-变异统计时不包含12基因-2025.08.11
def gsww_mp_filter_12(var_list):
	gene_list = ["MET", "EGFR", "ALK", "KRAS", "ROS1", "RET", "ERBB2", "BRAF", "NRAS", "PIK3CA", "BRCA1", "BRCA2"]
	return [var for var in var_list if not set(re.split(",", var["gene_symbol"])) & set(gene_list)]
jinja2.filters.FILTERS["gsww_mp_filter_12"] = gsww_mp_filter_12

# 国际部胚系项目-2025.08.12
# 1. 通过基因、变异分组和性别匹配遗传肿瘤
# 2. 相同基因相同配置信息的放在一个表格里展示
# 3. 疾病不同的分多行展示
# 2025.08.12-新增筛选条件，根据纯合杂合-AR/AD进行筛选
def international_risk_v3(info):
	var_list = info[0]
	config = info[1]
	sample = info[2]
	for var in var_list:
		category_risk = []
		var_category_names = var["var_category_names"] if var["var_category_names"] else ""
		for category in re.split(",", var_category_names):
			if var["gene_symbol"] in config.keys() and category in config[var["gene_symbol"]].keys():
				category_risk.extend(config[var["gene_symbol"]][category])
		
		# 筛选出无性别区分或性别相同的信息
		filter_category_risk = [i for i in category_risk if i["gender"] == "-" or i["gender"].lower() == sample["gender"].lower()]
		# 2025.08.19-新增筛选规则：
		# 1. snvindel 杂合变异，筛选出disease对应为AD、-、XL、SD的信息
		# 2. snvindel 纯合变异，筛选出disease对应为AR的信息
		# 3. cnv 均按disease对应AD来筛选
		gene_type = "纯合" if var["bio_category"] == "Snvindel" and float(var["freq"]) >= 0.85 else "杂合"
		filter_category_risk_gene_type = []
		if gene_type == "纯合":
			filter_category_risk_gene_type = [i for i in filter_category_risk if i["disease_hereditary_pattern"] == "AR"]
		else:
			filter_category_risk_gene_type = [i for i in filter_category_risk if i["disease_hereditary_pattern"] in ["AD", "-", "XL", "SD"]]
		# 2025.08.19-新增完成
		disease_dict = {}
		for cate_item in filter_category_risk_gene_type:
			# 2025.08.12-不要（AD）
			#tumor_str_ = "{0} ({1})".format(cate_item["tumor"], cate_item["tumor_hereditary_pattern"]) if cate_item["tumor_hereditary_pattern"] != "-" else cate_item["tumor"]+" (-)"
			disease_str_ = "{0} ({1})".format(cate_item["disease"], cate_item["disease_hereditary_pattern"]) if cate_item["disease_hereditary_pattern"] != "-" else cate_item["disease"]+" (-)"
			tumor_str_ = "{0}".format(cate_item["tumor"])
			if disease_str_ not in disease_dict.keys():
				disease_dict.setdefault(disease_str_, [])
			disease_dict[disease_str_].append(tumor_str_)
		var["disease_risk"] = disease_dict
	
	tmp_dict = {}
	for var in var_list:
		key = (var["gene_symbol"], tuple(var["disease_risk"]))
		if key not in tmp_dict.keys():
			tmp_dict.setdefault(key, [])
		tmp_dict[key].append(var)
	
	result = []
	for key, var_l in tmp_dict.items():
		result.append(
			{
				"var_list" : var_l,
				"disease_risk" : [{"disease" : k, "tumor" : v} for k,v in var_l[0]["disease_risk"].items()]
			}
		)
	return result
jinja2.filters.FILTERS["international_risk_v3"] = international_risk_v3

# 国际部MP，HD不考虑等级，在2.1和2.2都要展示-2025.08.12
# HD变异可能涉及I/II/肿瘤发生发展相关/III
def inter_mp_filter_hd(info):
	var_list = info[0]
	judge_hd = info[1]
	if judge_hd == "inhd":
		return [var for var in var_list if var["bio_category"] == "PHd"]
	else:
		return [var for var in var_list if var["bio_category"] != "PHd"]
jinja2.filters.FILTERS["inter_mp_filter_hd"] = inter_mp_filter_hd

# 华西体检，汇总各个变异对应的judge_ad_hete
def schx_150_judge_ad_hete(var_list):
	judge_ad_hete_list = []
	for var in var_list:
		judge_ad_hete_list.extend(var["judge_ad_hete"])
	return judge_ad_hete_list
jinja2.filters.FILTERS["schx_150_judge_ad_hete"] = schx_150_judge_ad_hete

# 广东医科ptBPTM Plus-ECtype表格的POLE/TP53需要展示体细胞I/II和胚系4/5-2025.08.18
def gdyk_ptbptmplus_ec_sum(info):
	var_list = info[0]
	gene = info[1]
	if not var_list:
		var_list = [{"gene_symbol" : gene, "result" : "nofound"}]
	return var_list
jinja2.filters.FILTERS["gdyk_ptbptmplus_ec_sum"] = gdyk_ptbptmplus_ec_sum

# 温附一BPTM Plus血液-所有基因检测结果汇总-2025.08.19
# 按BRCA、林奇5基因、其他基因排
def get_gbptm_plus_summary_v2(var_list):
	gene_list = ["BRCA1", "BRCA2", "EPCAM", "MLH1", "MSH2", "MSH6", "PMS2", "POLE", "TP53", "CTNNB1"]
	detect_gene_list = [var["gene_symbol"] for var in var_list]
	for gene in gene_list:
		if gene not in detect_gene_list:
			var_list.append({
				"gene_symbol" : gene
			})
	return sorted(var_list, key = lambda i:gene_list.index(i["gene_symbol"]))
jinja2.filters.FILTERS["get_gbptm_plus_summary_v2"] = get_gbptm_plus_summary_v2

# 温附一BPTM Plus全血判断是否检出林奇相关基因变异-2025.08.19
def wfy_gbptm_judge_lyn(info):
	var_list = info[0]
	gene_type = info[1]
	gene_list = ["EPCAM", "MLH1", "MSH2", "MSH6", "PMS2"]
	lyn_var = [var for var in var_list if var["gene_symbol"] in gene_list]
	other_var = [var for var in var_list if var["gene_symbol"] not in gene_list]
	if gene_type == "lyn":
		return lyn_var
	else:
		return other_var
jinja2.filters.FILTERS["wfy_gbptm_judge_lyn"] = wfy_gbptm_judge_lyn

# 温附一BPTM Plus组织-所有基因检测结果汇总-2025.08.19
# MSI放基因前展示，这边就删掉MSI了
def get_tbptm_plus_summary_v2(var_list):
	gene_list = ["POLE", "TP53", "BRCA1", "BRCA2", "CTNNB1", "EPCAM", "MLH1", "MSH2", "MSH6", "PMS2"]
	detect_gene_list = [var["gene_symbol"] for var in var_list]
	for gene in gene_list:
		if gene not in detect_gene_list:
			var_list.append({
				"gene_symbol" : gene
			})
	return sorted(var_list, key=lambda i:gene_list.index(i["gene_symbol"]))
jinja2.filters.FILTERS["get_tbptm_plus_summary_v2"] = get_tbptm_plus_summary_v2

# 适用北大三CP200-2025.08.19
# 新增HD-CDKN2A/CDKN2B
def io_detect_for_BDS_CP200_v2(info):
	var_list = info[0]
	return_type = info[1]
	io_result = {}
	# 汇总体细胞I/II/肿瘤发生发展相关变异+胚系致病/疑似致病性变异
	io_gene_P = ["ATM","ATR","BRCA1","BRCA2","BRIP1","CHEK1","CHEK2","FANCA","MRE11",\
				 "PALB2","RAD50","MLH1","MSH2","MSH6","PMS2","POLE","POLD1","TP53",\
				 "KRAS","CD274","ARID1A","SETD2","TERT","KMT2D","FAT1","CDK12"]
	io_gene_N = ["EGFR","ALK","MDM2","MDM4","CDKN2A","CDKN2B","DNMT3A","STK11","IFNGR1",\
				 "JAK1","JAK2","APC","CTNNB1","B2M","PTEN","FGF19"]
	cnv_gene_list = ["CD274", "MDM2", "MDM4", "FGF19"]
	
	for var in var_list:
		# 仅展示扩增的基因
		if var["bio_category"] == "Cnv" and var["gene_symbol"] in cnv_gene_list:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			io_result[var["gene_symbol"]].append(var)
		# 仅展示融合的基因
		elif var["bio_category"] in ["Sv", "PSeqRnaSv"] and set(re.split(",", var["gene_symbol"])) & set(["ALK"]):
			if "ALK" not in io_result.keys():
				io_result.setdefault("ALK", [])
			io_result["ALK"].append(var)
		# 展示HD和snvindel
		elif var["bio_category"] in ["Snvindel", "PHd"] and var["gene_symbol"] in ["CDKN2A", "CDKN2B"]:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			io_result[var["gene_symbol"]].append(var)
		# 其余基因展示Snvindel
		elif var["bio_category"] == "Snvindel" and var["gene_symbol"] in io_gene_P + io_gene_N:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			io_result[var["gene_symbol"]].append(var)

	# summary展示
	io_p_list = [i for k,v in io_result.items() for i in v if k in io_gene_P]
	io_n_list = [i for k,v in io_result.items() for i in v if k in io_gene_N]
	
	if return_type == "p":
		return io_p_list
	elif return_type == "n":
		return io_n_list
jinja2.filters.FILTERS["io_detect_for_BDS_CP200_v2"] = io_detect_for_BDS_CP200_v2

# 2025.08.28-国际部胚系2.1，变异放在一个表格中展示，对应肿瘤也放同个表格里
# 注意去重
def international_risk_merge(var_list):
	result = []
	for var in var_list:
		if var["disease_risk"]:
			for i in var["disease_risk"]:
				if i not in result:
					result.append(i)
	return result
jinja2.filters.FILTERS["international_risk_merge"] = international_risk_merge

# 2025.09.02-郑大一小结-兼容CNV和MLPA，且大片段重复/缺失变异等级改为动态判定
def ZDY_HRR_germline_summary_v2(var_list):
	var_dict = {}
	for var in var_list:
		if var["type"] == "Loss":
			var_info = "{0} {1} del".format(var["gene_symbol"], var["value"])
		elif var["type"] == "Gain":
			var_info = "{0} {1} dup".format(var["gene_symbol"], var["value"])
		else:
			if var["hgvs_p"] != "p.?":
				var_info = "{0} {1}:{2}:{3}".format(var["gene_symbol"], var["transcript_primary"], var["hgvs_c"], var["hgvs_p"])
			else:
				var_info = "{0} {1}:{2}".format(var["gene_symbol"], var["transcript_primary"], var["hgvs_c"])
		if var["clinic_num_g"] not in var_dict.keys():
			var_dict.setdefault(var["clinic_num_g"], [])
		var_dict[var["clinic_num_g"]].append(var_info)
	result = []
	if 5 in var_dict.keys():
		result.append("检出{0}个致病性变异，{1}".format(str(len(var_dict[5])), ", ".join(var_dict[5])))
	if 4 in var_dict.keys():
		result.append("检出{0}个疑似致病性变异，{1}".format(str(len(var_dict[4])), ", ".join(var_dict[4])))
	if not var_dict:
		result.append("未检出致病或疑似致病性胚系变异")
	return "；".join(result)
jinja2.filters.FILTERS["ZDY_HRR_germline_summary_v2"] = ZDY_HRR_germline_summary_v2

# 福建附一CP200 小结处，10个基因结果总结-2025.09.08
def FJFY_CP200_10gene_summary(info):
	var_list = info[0]
	level = info[1]
	# 福建附一EML4-ALK融合只展示1-3型
	def eml4_alk_type(merge_sv_list):
		# gene1:exon1-gene2:exon2
		alk_sv_dict = {
			("exon13", "exon20") : "V1",
			("exon20", "exon20") : "V2",
			("exon6", "exon20") : "V3a/b"
		}
		result = []
		for var in merge_sv_list:
			five_prime_cds = re.split("-", re.split(":", var)[1])[0]
			three_prime_cds = re.split("-", re.split(":", var)[-1])[0]
			key = (five_prime_cds, three_prime_cds)
			sv_type = alk_sv_dict.get(key, "")
			if sv_type:
				result.append(var+"（{0}）".format(sv_type))
			else:
				result.append(var)
		return "、".join(result)
	
	def var_info_stran(var):
		var_info = ""
		region_dict = {
		"exon" : "外显子",
		"intron" : "内含子",
		"3'UTR" : "3'UTR",
		"5'UTR" : "5'UTR",
		"3'FLANKING" : "非编码区",
		"5'FLANKING" : "非编码区"
		}
		if var["bio_category"] == "Cnv":
			var_info = var["gene_symbol"]+"基因检测到扩增，拷贝数为"+str(var["cn_mean"])
		elif var["bio_category"] == "Sv":
			var_desc_merge = var["var_desc_merge"] if "var_desc_merge" in var.keys() and var["var_desc_merge"] else ""
			if var["five_prime_gene"] == "MET" and var["three_prime_gene"] == "MET":
				var_info = var["gene_symbol"]+"基因检测到"+var["five_prime_gene"]+"-"+var["three_prime_gene"]+"融合，拷贝数为"+str(var["copies"])+" copies"
			elif var["five_prime_gene"] == "EML4" and var["three_prime_gene"] == "ALK":
				var_info = var["gene_symbol"]+"基因检测到"+var["five_prime_gene"]+"-"+var["three_prime_gene"]+"融合，融合模式为"+ eml4_alk_type(var["merge_sv_list"]) +"，拷贝数为"+str(var["copies"])+" copies"
			else:
				var_info = var["gene_symbol"]+"基因检测到"+var["five_prime_gene"]+"-"+var["three_prime_gene"]+"融合，融合模式为"+var_desc_merge+"，拷贝数为"+str(var["copies"])+" copies"
		elif var["bio_category"] == "Snvindel":
			# 加一个MET exon 14跳跃突变-2023.07.13
			if var["gene_symbol"] == "MET" and "exon14 skipping" in var["variant_interpret_cn"]:
				if var["hgvs_p_FJFY"] != "p.?":
					var_info  = "MET基因14号外显子跳跃突变"+var["hgvs_c"]+"("+var["hgvs_p_FJFY"]+")，突变丰度为"+str(var["freq_str"])
				else:
					var_info = "MET基因14号外显子跳跃突变"+var["hgvs_c"]+"，突变丰度为"+str(var["freq_str"])
			# MET exon14跳跃添加结束-2023.07.13
			else:
				region_list_en = re.split("_", var["gene_region"])
				region_list_cn = []
				for i in region_list_en:
					if re.search("exon", i):
						region_list_cn.append(i.replace("exon", "")+"号外显子")
					elif re.search("intron", i):
						region_list_cn.append(i.replace("intron", "")+"号内含子")
					else:
						region_cn = region_dict[i] if i in region_dict.keys() else i
						region_list_cn.append(region_cn)
				if var["hgvs_p_FJFY"] != "p.?":
					var_info = var["gene_symbol"]+"基因"+"到".join(region_list_cn)+"检测到" + var["type_cn"] + var["hgvs_p_FJFY"]+"，突变丰度为"+str(var["freq_str"])
				else:
					if var["type_cn"] != "--":
						var_info = var["gene_symbol"]+"基因"+"到".join(region_list_cn)+"检测到"+ var["type_cn"] + var["hgvs_c"]+"突变，突变丰度为"+str(var["freq_str"])
					else:
						var_info = var["gene_symbol"]+"基因"+"到".join(region_list_cn)+"检测到"+ var["hgvs_c"]+"突变，突变丰度为"+str(var["freq_str"])
		return var_info

	def regimen_info(var):
		sense_regimen = [i["regimen_name"] for i in var["evi_sum"]["evi_split"]["Predictive"] if \
							re.search("Sensitive", i["clinical_significance"]) and \
							i["evi_conclusion_simple"]=="A"] if \
								"Predictive" in var["evi_sum"]["evi_split"].keys() else \
								[]
		resis_regimen = [i["regimen_name"] for i in var["evi_sum"]["evi_split"]["Predictive"] if \
							re.search("Resistant", i["clinical_significance"]) and \
							i["evi_conclusion_simple"]=="A"] if \
								"Predictive" in var["evi_sum"]["evi_split"].keys() else \
								[]

		sense_regimen_str = "对" + "、".join(sense_regimen) + "可能敏感（A级）" if sense_regimen else ""
		resis_regimen_str = "对" + "、".join(resis_regimen) + "可能耐药（A级）" if resis_regimen else ""

		return sense_regimen_str, resis_regimen_str
	
	result_list = []
	result = ""
	lc10_gene = ["BRAF", "EGFR", "ERBB2", "KRAS", "MET", "ALK", "ROS1", "RET", "NRAS", "PIK3CA"]
	gene10_var_list = [var for var in var_list if set(lc10_gene) & set(re.split(",", var["gene_symbol"]))]
	if level == "I":
		for var in gene10_var_list:
			tmp_list = []
			var_info = var_info_stran(var)
			sense_regimen_str, resis_regimen_str = regimen_info(var)
			tmp_list.append(var_info)
			if sense_regimen_str:
				tmp_list.append(sense_regimen_str)
			if resis_regimen_str:
				tmp_list.append(resis_regimen_str)
			if "，".join(tmp_list) not in result_list:
				result_list.append("，".join(tmp_list))
		result = "；".join(result_list)
	else:
		for var in gene10_var_list:
			if var_info_stran(var) and var_info_stran(var) not in result_list:
				result_list.append(var_info_stran(var))
		result = "；".join(result_list)
	
	return result
jinja2.filters.FILTERS["FJFY_CP200_10gene_summary"] = FJFY_CP200_10gene_summary

# 福建附一CP200 EML4-ALK融合需要展示具体的亚型（仅1/2/3，存在123时同时检出其他亚型，标红处理）-2025.09.08
def FJFY_CP200_10gene_var(merge_sv_list):
	# gene1:exon1-gene2:exon2
	alk_sv_dict = {
		("exon13", "exon20") : "V1",
		("exon20", "exon20") : "V2",
		("exon6", "exon20") : "V3a/b",
		("exon15", "exon20") : "V4",
		("exon2", "exon20") : "V5a/b",
		("exon18", "exon20") : "V5",
		("exon14", "exon20") : "V7",
		("exon17", "exon20") : "V8a/b"
	}
	v123 = []
	vother = []
	for sv_info in merge_sv_list:
		five_prime_cds = re.split("-", re.split(":", sv_info)[1])[0]
		three_prime_cds = re.split("-", re.split(":", sv_info)[-1])[0]
		if alk_sv_dict.get((five_prime_cds, three_prime_cds), "") in ["V1", "V2", "V3a/b"] and alk_sv_dict.get((five_prime_cds, three_prime_cds), "") not in v123:
			v123.append(alk_sv_dict.get((five_prime_cds, three_prime_cds), ""))
		elif alk_sv_dict.get((five_prime_cds, three_prime_cds), "") in ["V4", "V5", "V5a/b", "V7", "V8a/b"]:
			vother.append(alk_sv_dict.get((five_prime_cds, three_prime_cds), ""))	

	return ", ".join(v123), vother
jinja2.filters.FILTERS["FJFY_CP200_10gene_var"] = FJFY_CP200_10gene_var

# 温附一HRR，过滤出指定基因结果-20250916
def wfy_filter_gene_var(info):
	var_list = info[0]
	gene_symbol = info[1]
	return [var for var in var_list if var["gene_symbol"] == gene_symbol]
jinja2.filters.FILTERS["wfy_filter_gene_var"] = wfy_filter_gene_var

# 重庆肿瘤116基因-检测结果小结-具有临床意义的变异需要区分为10基因和其他基因-2025.09.16
# 需要增加MET 14 跳跃突变标记
def cqzl_116_var12_summary_v2(info):
	all_var_list = info[0]
	gene_type = info[1]
	LC10_gene_list = ["ALK", "BRAF", "EGFR", "ERBB2", "KRAS", "NRAS", "MET", "PIK3CA", "RET", "ROS1"]
	var_list = []
	if gene_type == "inLC10":
		var_list = [var for var in all_var_list if set(re.split(",", var["gene_symbol"])) & set(LC10_gene_list)]
	elif gene_type == "outLC10":
		var_list = [var for var in all_var_list if not set(re.split(",", var["gene_symbol"])) & set(LC10_gene_list)]
	# 2024.11.08-加一个全部返回的
	else:
		var_list = all_var_list
	
	for var in var_list:
		if var["bio_category"] == "Snvindel":
			var["cqzl_var_info"] = var["hgvs_p"] if var["hgvs_p"] != "p.?" else var["hgvs_c"]
			if "var_category_names" in var.keys() and var["var_category_names"] and "MET Exon14 Skipping" in var["var_category_names"]:
				var["cqzl_var_info"] = var["cqzl_var_info"] + "（MET 14跳跃缺失） "
		elif var["bio_category"] == "Cnv":
			var["cqzl_var_info"] = "扩增"
		elif var["bio_category"] == "Sv":
			var["cqzl_var_info"] = "融合"
		else:
			var["cqzl_var_info"] = ""
	
	# 每个变异（除了最后一个）cqzl_var_info后面加个逗号，在模板中使用for循环展示变异
	if len(var_list) > 1:
		for var in var_list[0:-1]:
			var["cqzl_var_info"] = var["cqzl_var_info"] + "，"
	return var_list
jinja2.filters.FILTERS["cqzl_116_var12_summary_v2"] = cqzl_116_var12_summary_v2

# 解放军总医院第七医学中心 master-需要展示指定POLE变异检测情况-2025.09.16
def jfjzydq_mp_filter_pole(var_list):
	pole_result = []
	for var in var_list:
		if var["bio_category"] == "Snvindel" and var["gene_symbol"] == "POLE" and var["hgvs_p"].replace("p.(", "").replace(")", "") in ["P286R", "V411L", "A456P"]:
			if var["hgvs_p"].replace("p.(", "").replace(")", "") not in pole_result:
				pole_result.append(var["hgvs_p"].replace("p.(", "").replace(")", ""))
	if pole_result:
		return "POLE {0}突变".format("、".join(pole_result))
	else:
		return ""
jinja2.filters.FILTERS["jfjzydq_mp_filter_pole"] = jfjzydq_mp_filter_pole

# 中六150-功能改变-2025.09.18
# Splicing +- 1/2 归为Splice，gene_region中同时含有exon和intron的归为Span
def zsly_150_vartype(var):
	type_stran = {
		"nonSynonymous_Substitution" : "Missense",
		"Nonsense_Mutation" : "Nonsense",
		"FrameShift_Deletion" : "Frameshift",
		"FrameShift_Duplication" : "Frameshift",
		"FrameShift_Insertion" : "Frameshift",
		"FrameShift_Substitution" : "Frameshift",
		"nonFrameShift_Insertion" : "Insertion",
		"nonFrameShift_Duplication" : "Insertion",
		"nonFrameShift_Deletion" : "Deletion",
		"nonFrameShift_Substitution" : "Delins",
		"Intronic" : "Intronic",
		"Synonymous_Substitution" : "Synonymous",
		"3'UTR" : "UTR",
		"5'UTR" : "UTR",
		"CorruptedStart" : "Alt_start",
		"FlankingRegion5" : "Flanking",
		"FlankingRegion3" : "Flanking",
		"Extension" : "Extension"
	}
	if var["type"] == "Splicing":
		if "exon" in var["gene_region"] and "intron" in var["gene_region"]:
			return "Span"
		else:
			return "Splice"
	else:
		return type_stran.get(var["type"], var["type"])
jinja2.filters.FILTERS["zsly_150_vartype"] = zsly_150_vartype

# 中六150-遗传方式-2025.09.19
# 使用国际部配置表，看disease_hereditary_pattern列
# 4/5类变异需要匹配变异分组，3类？（先都按变异分组来匹配，不管等级）
# XL和SD归哪里待定吧
def zsly_150_hereditary_pattern(info):
	var = info[0]
	config = info[1]
	# 4/5类变异
	category_risk = []
	var_category_names = var["var_category_names"] if "var_category_names" in var.keys() and var["var_category_names"] else ""
	if not var_category_names:
		var_category_names = var["var_name"] if "var_name" in var.keys() and var["var_name"] else ""
	for category in re.split(",", var_category_names):
		if var["gene_symbol"] in config.keys() and category in config[var["gene_symbol"]].keys():
			category_risk.extend(config[var["gene_symbol"]][category])

	pattern_list = []
	for i in category_risk:
		if i["disease_hereditary_pattern"] and i["disease_hereditary_pattern"] != "-" and i["disease_hereditary_pattern"] not in pattern_list:
			pattern_list.append(i["disease_hereditary_pattern"])
	
	return "+".join(pattern_list)
jinja2.filters.FILTERS["zsly_150_hereditary_pattern"] = zsly_150_hereditary_pattern

# 吉林大学第一医院-116-本癌种重要基因要根据癌种进行区分-V3-2025.09.19
# 新增癌种
def jlyy_116_important_gene_v3(info):
	var_list = info[0]
	tumor_list = info[1]
	tumor_dict = {
		"肺癌" : ["ALK", "EGFR", "ROS1", "RET", "MET", "ERBB2", "KRAS", "BRAF", "NTRK1", "NTRK2", "NTRK3"],
		"肠癌" : ["KRAS", "NRAS", "BRAF"],
		"胃癌" : ["ERBB2", "NTRK1", "NTRK2", "NTRK3", "RET", "BRAF"],
		"甲状腺癌" : ["ALK", "BRAF", "RET", "NTRK1", "NTRK2", "NTRK3"],
		"黑色素瘤" : ["BRAF", "KIT", "NRAS", "NTRK1", "NTRK2", "NTRK3", "RET"],
		"卵巢癌" : ["BRAF", "BRCA1", "BRCA2", "KRAS", "NTRK1", "NTRK2", "NTRK3", "RET"],
		"子宫内膜癌" : ["MLH1", "MSH2", "MSH6", "PMS2", "POLE", "TP53"],
		"肾癌" : ["BRAF", "NTRK1", "NTRK2", "NTRK3", "RET", "VHL"],
		"膀胱癌" : ["ERBB2", "FGFR2", "FGFR3", "NTRK1", "NTRK2", "NTRK3"],
		"尿路上皮癌" : ["ERBB2", "FGFR2", "FGFR3", "NTRK1", "NTRK2", "NTRK3", "RET"],
		"乳腺癌" : ["AKT1", "BRCA1", "BRCA2", "ERBB2", "ESR1", "NTRK1", "NTRK2", "NTRK3", "PIK3CA", "PTEN", "RET"],
		"鼻咽癌" : ["EGFR"],
		"涎腺肿瘤" : ["AR", "ERBB2", "HRAS", "NTRK1", "NTRK2", "NTRK3", "PIK3CA"],
		"胶质瘤" : ["BRAF", "NTRK1", "NTRK2", "NTRK3"],
		"前列腺癌" : ["BRAF", "BRCA1", "BRCA2", "NTRK1", "NTRK2", "NTRK3", "RET"],
		"肝癌" : ["NTRK1", "NTRK2", "NTRK3", "RET"],
		"胆管癌" : ["BRAF", "ERBB2", "FGFR2", "IDH1", "KRAS", "NRG1", "NTRK1", "NTRK2", "NTRK3", "RET"],
		"胰腺癌" : ["BRAF", "BRCA1", "BRCA2", "ERBB2", "FGFR1", "FGFR2", "FGFR3", "KRAS", "NRG1", "NTRK1", "NTRK2", "NTRK3", "PALB2", "RET"],
		"胃肠道间质瘤" : ["BRAF", "KIT", "NF1", "NTRK1", "NTRK2", "NTRK3", "PDGFRA", "RET"],
		"宫颈癌" : ["ERBB2", "NTRK1", "NTRK2", "NTRK3", "RET"]
	}
	detect_gene = []
	for var in var_list:
		for gene in re.split(",", var["gene_symbol"]):
			if gene not in detect_gene:
				detect_gene.append(gene)
	nccn_gene_list = []
	for tumor in tumor_list:
		if tumor in tumor_dict.keys():
			nccn_gene_list.extend(tumor_dict[tumor])
	nccn_gene_list = list(set(nccn_gene_list))
	return "，".join(sorted([gene for gene in nccn_gene_list if gene not in detect_gene]))
jinja2.filters.FILTERS["jlyy_116_important_gene_v3"] = jlyy_116_important_gene_v3

# 临沂肿瘤BRCA、HRR等-诊断结果-2025.09.22
def lyzl_germline_sum(var_list):
	result = []
	for var in var_list:
		if var["type"] == "Loss":
			result.append("{0}基因{1} del突变".format(var["gene_symbol"], var["value"]))
		elif var["type"] == "Gain":
			result.append("{0}基因{1} dup突变".format(var["gene_symbol"], var["value"]))
		elif var["bio_category"] == "Snvindel":
			if var["hgvs_p"] != "p.?":
				result.append("{0}基因{1}突变".format(var["gene_symbol"], var["hgvs_p"]))
			else:
				result.append("{0}基因{1}突变".format(var["gene_symbol"], var["hgvs_c"]))
	return "，".join(result)
jinja2.filters.FILTERS["lyzl_germline_sum"] = lyzl_germline_sum

# 浙江人民HRR，检测结果总结需要列出治疗方案，<=3个治疗方案按实际情况写，>3个的只取前三个-2025.09.22
def zjrm_hrr_sum(var_list):
	result = []
	for var in var_list:
		regimen_list = []
		if "evi_split" in var["evi_sum"].keys():
			if "Predictive" in var["evi_sum"]["evi_split"].keys():
				for i in var["evi_sum"]["evi_split"]["Predictive"]:
					regimen_list.append("{0}（{1}，{2}级）".format(i["regimen_name"], i["clinical_significance_cn"], i["evi_conclusion_simple"]))
		if len(regimen_list) == 0:
			var["zjrm_hrr_sum_str"] = ""
		elif len(regimen_list) <= 3:
			var["zjrm_hrr_sum_str"] = "推荐{0}。".format("、".join(regimen_list))
		else:
			var["zjrm_hrr_sum_str"] = "推荐{0}等，详见第四部分详细检测结果。".format("、".join(regimen_list[0:3]))
		if regimen_list:
			result.append(var)
	return result
jinja2.filters.FILTERS["zjrm_hrr_sum"] = zjrm_hrr_sum

# 毓璜顶oncopro41基因-NCCN指南推荐基因检测结果，需要过滤及区分癌种-2025.09.29
def yhd_oncopro41_nccn_filter(info):
	sample = info[0]
	nccn_list = info[1]
	gene_dict = {
		"AKT1" : "NM_001382430.1",
		"ALK" : "NM_004304.5",
		"BRAF" : "NM_004333.6",
		"CTNNB1" : "NM_001904.4",
		"DDR2" : "NM_006182.4",
		"DPYD" : "NM_000110.4",
		"EGFR" : "NM_005228.5",
		"ERBB2" : "NM_004448.4",
		"ESR1" : "NM_000125.4",
		"FGFR1" : "NM_023110.3",
		"FGFR2" : "NM_000141.5",
		"FGFR3" : "NM_000142.5",
		"FGFR4" : "NM_213647.3",
		"HRAS" : "NM_005343.4",
		"IDH1" : "NM_005896.4",
		"IDH2" : "NM_002168.4",
		"KEAP1" : "NM_203500.2",
		"KIT" : "NM_000222.3",
		"KRAS" : "NM_004985.5",
		"MAP2K1" : "NM_002755.4",
		"MET" : "NM_000245.4",
		"NFE2L2" : "NM_006164.5",
		"NRAS" : "NM_002524.5",
		"NTRK1" : "NM_002529.4",
		"NTRK2" : "NM_006180.6",
		"NTRK3" : "NM_001012338.3",
		"PDGFRA" : "NM_006206.6",
		"PIK3CA" : "NM_006218.4",
		"POLD1" : "NM_002691.4",
		"POLE" : "NM_006231.4",
		"PTEN" : "NM_000314.8",
		"RB1" : "NM_000321.3",
		"RET" : "NM_020975.6",
		"ROS1" : "NM_002944.2",
		"STK11" : "NM_000455.5",
		"TP53" : "NM_000546.6",
		"UGT1A1" : "NM_000463.3",
		"CDK4" : "NM_000075.4",
		"MYC" : "NM_002467.6",
		"NKX2-1" : "NM_001079668.3",
		"NRG1" : "NM_013964.5"
	}
	
	result = []
	for i in nccn_list:
		if "实体瘤" in i["disease"] or set(sample["tumor_list"]) & set(re.split("、", i["disease"])):
			if i["gene_symbol"] in gene_dict.keys():
				i["transcript_primary"] = gene_dict.get(i["gene_symbol"])
				result.append(i)
	return result
jinja2.filters.FILTERS["yhd_oncopro41_nccn_filter"] = yhd_oncopro41_nccn_filter

# 战略支援部队特色医学中心 MP组织结果小结-io-基因斜体-2025.10.09
def jfjdj_mp_io_sum(io_result):
	result = []
	for i in io_result:
		for var in i["var_info"]:
			tmp_dict = {}
			if var["var_type"] == "snvindel":
				tmp_dict = {
					"gene_symbol" : i["gene_symbol"],
					"var_info" : var["hgvs_p"] if var["hgvs_p"] != "p.?" else var["hgvs_c"],
					"var_type" : "snvindel"
				}
			elif var["var_type"] == "cnv":
				tmp_dict = {
					"gene_symbol" : i["gene_symbol"],
					"var_info" : "扩增",
					"var_type" : "cnv"
				}
			elif i["gene_symbol"] == "CCND1/FGF3/FGF19":
				tmp_dict = {
					"gene_symbol" : i["gene_symbol"],
					"var_info" : "共扩增",
					"var_type" : "cnv"
				}
			elif var["var_type"] == "sv":
				tmp_dict = {
					"gene_symbol" : i["gene_symbol"],
					"var_info" : "融合",
					"var_type" : "sv",
					"three_prime_gene" : var["three_prime_gene"],
					"three_prime_cds" : var["three_prime_cds"],
					"five_prime_gene" : var["five_prime_gene"],
					"five_prime_cds" : var["five_prime_cds"]
				}
			elif var["var_type"] == "PHd":
				tmp_dict = {
					"gene_symbol" : i["gene_symbol"],
					"var_info" : "纯合缺失" if var["type"] == "HomoDel" else "杂合缺失" if var["type"] == "HeteDel" else "未知变异类型",
					"var_type" : "PHd"
				}
			if tmp_dict not in result:
				result.append(tmp_dict)
	# 每个变异（除了最后一个）cqzl_var_info后面加个逗号，在模板中使用for循环展示变异
	if len(result) > 1:
		for var in result[0:-1]:
			var["var_info"] = var["var_info"] + ", "
	return result
jinja2.filters.FILTERS["jfjdj_mp_io_sum"] = jfjdj_mp_io_sum

# 2025.07.17-重庆西南116-区分10基因获批位点和其他位点
# 2025.10.10-新增需求-10基因不在获批位点列表的I类变异也要展示
def cqxn_10gene_appr_var_v3(info):
	# 判断是为I类变异
	def judge_level_I_var(var):
		if var["evi_sum"]["evi_split"] and set(["Diagnostic","Predictive","Prognostic"]) & set(var["evi_sum"]["evi_split"].keys()):
			if var["clinic_num_s"] == 5:
				return 1
			else:
				return 0
		else:
			return 0

	var_list = info[0]
	gene_type = info[1]
	gene_list = ["ALK", "BRAF", "EGFR", "ERBB2", "KRAS", "MET", "NRAS", "PIK3CA", "RET", "ROS1"]
	appr_dict = {
		"EGFR" : ["p.G719A","p.E746_A750del","p.L747_T751del","p.V769_D770insASV","p.V774_C775insHV",\
				  "p.H773_V774insAH","p.H773_V774insY","p.P772_H773insPNP","p.V769_D770insGG","p.N771_P772insRHN",\
				  "p.L861Q","p.G719S","p.L747_P753delinsS","p.T790M","p.D770_N771insSVD",\
				  "p.D770delinsGY","p.N771_P772insN","p.H773_V774insTH","p.N771_P772insG","p.V769_D770insGTL",\
				  "p.P772_H773insGHP","p.G719C","p.E746_S752delinsV","p.S768I","p.A763_Y764insFQEA",\
				  "p.H773_V774insH","p.N771_P772insH","p.N771delinsGY","p.N771delinsGF","p.D770_N771insGF",\
				  "p.H773delinsYNPY","p.E746_A750del","p.L747_A750delinsP","p.D770_N771insG","p.H773_V774insNPH",\
				  "p.H773_V774insPH","p.P772_H773insYNP","p.P772_H773insQ","p.P772_H773insGNP","p.D770_N771insSLA","p.L858R"],
		"KRAS" : ["p.G12D", "p.G12A", "p.G12V", "p.G12S", "p.G12C", "p.Q61H"],
		"NRAS" : ["p.G12D", "p.Q61R", "p.Q61K"],
		"PIK3CA" : ["p.H1047R"],
		"BRAF" : ["p.V600E"],
		"ERBB2" : ["p.A775_G776insYVMA"],
		"MET" : ["c.3082+1G>T"],
		"ALK" : ["EML4:exon13-ALK:exon20", "EML4:exon6-ALK:exon20", "EML4:exon20-ALK:exon20"],
		"ROS1" : ["CD74:exon6-ROS1:exon34", "GOPC:exon8-ROS1:exon35"],
		"RET" : ["KIF5B:exon15-RET:exon12"]
	}
	appr_var = []
	other_var = []
	for var in var_list:
		if var["bio_category"] == "Snvindel":
			snv_info = var["hgvs_p_ZJZL"] if var["hgvs_p_ZJZL"] != "p.?" else var["hgvs_c"]
			if var["gene_symbol"] in appr_dict.keys() and snv_info in appr_dict[var["gene_symbol"]]:
				appr_var.append(var)
			elif var["gene_symbol"] in gene_list and judge_level_I_var(var):
				appr_var.append(var)
			else:
				other_var.append(var)
		elif var["bio_category"] == "Sv":
			sv_info = "{0}:{1}-{2}:{3}".format(var["five_prime_gene"], var["five_prime_cds"], var["three_prime_gene"], var["three_prime_cds"])
			# ALK/ROS1/RET几个获批变异型，主基因都在3端
			if var["three_prime_gene"] in appr_dict.keys() and sv_info in appr_dict[var["three_prime_gene"]]:
				appr_var.append(var)
			elif set(re.split(",", var["gene_symbol"])) & set(gene_list) and judge_level_I_var(var):
				appr_var.append(var)
			else:
				other_var.append(var)
		else:
			if var["gene_symbol"] in gene_list and judge_level_I_var(var):
				appr_var.append(var)
			else:
				other_var.append(var)
	detect_gene = []
	for var in appr_var:
		for gene in re.split(",", var["gene_symbol"]):
			if gene in gene_list and gene not in detect_gene:
				detect_gene.append(gene)
	result = []
	result.extend(appr_var)
	for gene in sorted(list(set(gene_list) - set(detect_gene))):
		result.append({"gene_symbol" : gene})
	if gene_type == "LC10":
		return result
	else:
		return other_var
jinja2.filters.FILTERS["cqxn_10gene_appr_var_v3"] = cqxn_10gene_appr_var_v3

# 2025.10.11-浙江省中医院BRCA-汇总药物
def zjzy_regimen_sum(appr_regimen_list):
	result = []
	for i in appr_regimen_list:
		if "regimen_cn" in i.keys() and i["regimen_cn"]:
			result.append(i["regimen_cn"])
		elif "regimen_en" in i.keys() and i["regimen_en"]:
			result.append(i["regimen_en"])
	return "，".join(result)
jinja2.filters.FILTERS["zjzy_regimen_sum"] = zjzy_regimen_sum

# 2025.10.16-西安交大一Master-判断输入的变异列表中是否包含证据
def xajdy_mp_judge_evi(var_list):
	result = False
	for var in var_list:
		if var["evi_sum"]["evi_split"] and set(["Diagnostic","Predictive","Prognostic"]) & set(var["evi_sum"]["evi_split"].keys()):
			result = True
			break
	return result
jinja2.filters.FILTERS["xajdy_mp_judge_evi"] = xajdy_mp_judge_evi

# 北三MP HD需要加基因座信息-2025.10.23
def bds_mp_hd_locus(gene_symbol):
	locus_dict = {
		"ATM" : "11q22.3",
		"BARD1" : "2q35",
		"BRCA1" : "17q21.31",
		"BRCA2" : "13q13.1",
		"BRIP1" : "17q23.2",
		"CDH1" : "16q22.1",
		"CDK12" : "17q12",
		"CDKN2A" : "9p21.3",
		"CDKN2B" : "9p21.3",
		"CHEK1" : "11q24.2",
		"CHEK2" : "22q12.1", 
		"FANCA" : "16q24.3",
		"FANCL" : "2p16.1",
		"FH" : "1q43",
		"HDAC2" : "6q21",
		"KEAP1" : "19p13.2",
		"MTAP" : "9p21.3",
		"PALB2" : "16p12.2",
		"PPP2R2A" : "8p21.2",
		"PTEN" : "10q23.31",
		"RAD51B" : "14q24.1",
		"RAD51C" : "17q22",
		"RAD51D" : "17q12",
		"RAD54L" : "1p34.1",
		"RB1" : "13q14.2",
		"SETD2" : "3p21.31",
		"SMARCA4" : "19p13.2",
		"STK11" : "19p13.3",
		"TP53" : "17p13.1"
	}

	return locus_dict.get(gene_symbol, "配置表中未找到对应信息，请检查！")
jinja2.filters.FILTERS["bds_mp_hd_locus"] = bds_mp_hd_locus

# 2025.10.28-返回肿瘤细胞含量或富集后肿瘤细胞含量
def get_tumor_content(info):
	sample = info[0]
	lib = info[1]
	key = info[2]
	if key in sample.keys() and sample[key]:
		return sample[key]
	else:
		if "lib_dna_qc" in lib.keys() and key in lib["lib_dna_qc"].keys() and lib["lib_dna_qc"][key]:
			return lib["lib_dna_qc"][key]
		else:
			return ""
jinja2.filters.FILTERS["get_tumor_content"] = get_tumor_content

# 2025.10.31-同济BPTM、tBRCA和gBRCA需要增加一句“19号外显子移码缺失突变”
def tj_var_info(var):
	result = ""
	type_stran_before = {
		"FrameShift_Duplication" : "移码突变",
		"nonFrameShift_Duplication" : "重复突变",
		"nonFrameShift_Insertion" : "插入突变",
		"FrameShift_Insertion" : "移码突变",
		"FrameShift_Deletion" : "移码突变",
		"nonFrameShift_Deletion" : "缺失突变",
		"FrameShift_Substitution" : "移码突变",
		"nonFrameShift_Substitution" : "缺失插入突变"
	}

	type_stran_after = {
		"FrameShift_Duplication" : "移码插入突变",
		"nonFrameShift_Duplication" : "框内插入突变",
		"nonFrameShift_Insertion" : "框内插入突变",
		"FrameShift_Insertion" : "移码插入突变",
		"FrameShift_Deletion" : "移码缺失突变",
		"nonFrameShift_Deletion" : "框内缺失突变",
		"FrameShift_Substitution" : "移码缺失插入突变",
		"nonFrameShift_Substitution" : "框内缺失插入突变"
	}
	if var["type"] == "Loss":
		result = var["value"].replace("exon", "") + "号外显子大片段缺失突变"
	elif var["type"] == "Gain":
		result = var["value"].replace("exon", "") + "号外显子拷贝数扩增"
	else:
		if var["type"] in type_stran_before.keys():
			result = var["varInfo_XAJDY"].replace("第", "").replace(type_stran_before[var["type"]], type_stran_after[var["type"]])
		else:
			result = var["varInfo_XAJDY"].replace("第", "")
	return result
jinja2.filters.FILTERS["tj_var_info"] = tj_var_info

# 山东省立ptHRR，需要总结胚系和体系变异-2025.11.10
# BRCA基因的hgvs_c（hgvs_p）
def sdsl_pthrr_sum(info):
	var_list = info[0]
	var_origin = info[1]
	germline = [var for var in var_list if var["type"] in ["Loss", "Gain"] or var["var_origin"] == "germline"]
	somatic = [var for var in var_list if "var_origin" in var.keys() and var["var_origin"] == "somatic"]
	result_var = germline if var_origin == "germline" else somatic
	result = []
	for var in result_var:
		if var["type"] in ["Loss", "Gain"]:
			cnv_type = "del" if var["type"] == "Loss" else "dup"
			result.append("{0}基因的{1} {2}".format(var["gene_symbol"], var["value"], cnv_type))
		else:
			if var["hgvs_p"] != "p.?":
				result.append("{0}基因的{1}（{2}）".format(var["gene_symbol"], var["hgvs_c"], var["hgvs_p"]))
			else:
				result.append("{0}基因的{1}".format(var["gene_symbol"], var["hgvs_p"]))
	return "、".join(result)
jinja2.filters.FILTERS["sdsl_pthrr_sum"] = sdsl_pthrr_sum

# 新疆附一CP200，I-IV类变异展示-2025.11.11
def xjfy_cp200_i_iv_sum(var_list):
	for var in var_list:
		if var["bio_category"] == "Snvindel":
			var["xjfy_var_info"] = var["hgvs_p"] if var["hgvs_p"] != "p.?" else var["hgvs_c"]
			if "judge_mergeMET" in var.keys() and var["judge_mergeMET"]:
				var["xjfy_var_info"] += "（MET exon14 skipping）"
		elif var["bio_category"] == "Cnv":
			var["xjfy_var_info"] = "扩增"
		elif var["bio_category"] == "Sv":
			if var["five_prime_gene"] == "MET" and var["three_prime_gene"] == "MET":
				var["xjfy_var_info"] = "exon14 skipping"
			else:
				var["xjfy_var_info"] = "融合"
	
	if len(var_list) > 1:
		for var in var_list[0:-1]:
			var["xjfy_var_info"] = var["xjfy_var_info"] + "，"
	return var_list
jinja2.filters.FILTERS["xjfy_cp200_i_iv_sum"] = xjfy_cp200_i_iv_sum

# 北京医院MP v2-统计含辅助诊断/预后相关证据的变异数-2025.11.12
def bjyy_mpv2_evi_type(info):
	var_list = info[0]
	evi_type = info[1]
	result = [var for var in var_list if set(var["evi_sum"]["evi_split"].keys()) & set([evi_type])]
	return len(result)
jinja2.filters.FILTERS["bjyy_mpv2_evi_type"] = bjyy_mpv2_evi_type

# 北京医院MP v2-CNV+HD单独展示-2025.11.12
def bjyy_filter_cnv_hd(var_list):
	return [var for var in var_list if var["bio_category"] in ["Cnv", "PHd"]]
jinja2.filters.FILTERS["bjyy_filter_cnv_hd"] = bjyy_filter_cnv_hd

# 北大三MP-疑似胚系变异需要在小结备注中展示出来-2025.11.12
def bds_likely_germline_sum(var_list):
	result = []
	for var in var_list:
		if var["hgvs_p_FJFY"] != "p.?":
			result.append("{0}:{1} ({2})突变丰度{3}".format(var["gene_symbol"], var["hgvs_c"], var["hgvs_p_FJFY"], var["freq_str"]))
		else:
			#print ("111111", var)
			result.append("{0}:{1}突变丰度{2}".format(var["gene_symbol"], var["hgvs_c"], var["freq_str"]))
	return "、".join(result)
jinja2.filters.FILTERS["bds_likely_germline_sum"] = bds_likely_germline_sum

# 复旦大学附属华山医院150-3类变异不展示同义、内含子、UTR变异-2025.11.12
def fdhs_level3_filter(var_list):
	filter_type = ["Synonymous_Substitution", "Intronic", "3'UTR", "5'UTR"]
	result = []
	for var in var_list:
		if var["bio_category"] == "Snvindel":
			if var["type"] not in filter_type:
				result.append(var)
		else:
			result.append(var)
	return result
jinja2.filters.FILTERS["fdhs_level3_filter"] = fdhs_level3_filter

# 复旦大学附属华山医院150-报告日期格式 2025年9月08日，月份第一位不补位，日期补位-2025.11.12
# 默认输入格式2025-09-08
def fdhs_report_date(indate):
	if indate:
		date_obj = datetime.strptime(str(indate), "%Y-%m-%d")
		result = f"{date_obj.year}年{date_obj.month}月{date_obj.strftime('%d')}日"
		return result
	else:
		return ""
jinja2.filters.FILTERS["fdhs_report_date"] = fdhs_report_date

# 2025.05.12-华西乳腺癌模板-获取临床意义
# 2025.05.21-EPCMA风险描述要区分变异分组，药物和预后没有先不管了-考虑直接加在模板里
# 2025.11.27-新增前列腺癌；额外匹配变异分组信息
# 2025.12.29-新增消化道肿瘤：结直肠癌和胰腺癌
# 目前一个基因对应只有一条信息，gene:info暂不需要改为gene:[info list]
def hx_breast_template_v2(info):
	var_list = info[0]
	tumor_list = info[1]
	hx_breast_template_significance = info[2]
	for var in var_list:
		num = 1
		var_significance = hx_breast_template_significance.get(var["gene_symbol"], {})
		if not (var_significance and set(re.split("\|\|", var_significance["var_category_names"])) & set(re.split(",", var["var_category_names"]))):
			var_significance = {}
		# 1. 预后信息-仅乳腺癌
		prognostic = var_significance["prognostic_breast"] if "乳腺癌" in tumor_list and \
															  "prognostic_breast" in var_significance and \
															  var_significance["prognostic_breast"]  and \
															  var_significance["prognostic_breast"] != "-" \
														   else ""
		prognostic_num = ""
		if prognostic:
			prognostic_num = str(num)
			num += 1
		# 2. 用药-区分乳腺癌/卵巢癌
		predictive = ""
		if "乳腺癌" in tumor_list:
			predictive = var_significance["predictive_breast"] if "predictive_breast" in var_significance and \
																  var_significance["predictive_breast"] and \
																  var_significance["predictive_breast"] != "-" \
															   else ""
		elif "卵巢癌" in tumor_list:
			predictive = var_significance["predictive_ovarian"] if "predictive_ovarian" in var_significance and \
																   var_significance["predictive_ovarian"] and \
																   var_significance["predictive_ovarian"] != "-" \
																else ""
		elif "前列腺癌" in tumor_list:
			predictive = var_significance["predictive_prostate"] if "predictive_prostate" in var_significance and \
																	var_significance["predictive_prostate"] and \
																	var_significance["predictive_prostate"] != "-" \
																 else ""
		# 2025.12.29-新增结直肠癌和胰腺癌
		elif "结直肠癌" in tumor_list:
			predictive = var_significance["predictive_colorectal"] if "predictive_colorectal" in var_significance and \
																	var_significance["predictive_colorectal"] and \
																	var_significance["predictive_colorectal"] != "-" \
																 else ""
		elif "胰腺癌" in tumor_list:
			predictive = var_significance["predictive_pancreatic"] if "predictive_pancreatic" in var_significance and \
																	var_significance["predictive_pancreatic"] and \
																	var_significance["predictive_pancreatic"] != "-" \
																 else ""
		# 2025.12.29-新增完成
		else:
			predictive = "nofound"
				
		predictive_num = ""
		if predictive:
			predictive_num = str(num)
			num += 1
		# 3. 风险提示
		risk_suggest = var_significance["risk_inter"] if "risk_inter" in var_significance and \
														  var_significance["risk_inter"] \
													  else ""
		risk_suggest_num = "" 
		if risk_suggest:
			risk_suggest_num = num
			num += 1

		result = {
			"risk_suggest" : risk_suggest,
			"risk_suggest_num" : risk_suggest_num,
			"prognostic" : prognostic,
			"prognostic_num" : prognostic_num,
			"predictive" : predictive,
			"predictive_num" : predictive_num,
			"refer" : var_significance.get("refer", "")
		}
		var["hx_breast_result"] = result
	return var_list
jinja2.filters.FILTERS["hx_breast_template_v2"] = hx_breast_template_v2

# 浙二HRR药物排序更新-2025.11.27
# 伊那利塞+哌柏西利+氟维司群 同等级排到第一位
def zjey_hrr_regimen_sort(regimen_list):
	for regimen in regimen_list:
		if regimen["regimen_name"] == "伊那利塞+哌柏西利+氟维司群":
			regimen["zjey_rule"] = 0
		else:
			regimen["zjey_rule"] = 1
	return sorted(regimen_list, key = lambda i : (i["evi_conclusion_simple"], i["zjey_rule"]))
jinja2.filters.FILTERS["zjey_hrr_regimen_sort"] = zjey_hrr_regimen_sort

# 浙二HRR药物排序更新后，详细解读合并相同描述的治疗方案-2025.11.27
def zjey_hrr_merge_evi(regimen_list):
	return merge_Predictive_evi(regimen_list)
jinja2.filters.FILTERS["zjey_hrr_merge_evi"] = zjey_hrr_merge_evi

# 福建附一CP200、116对治疗方案进行总结-2025.11.28
# 该突变提示【可能对XXX、XXX（证据等级为A），XXX、XXX（证据等级为B）敏感；可能对XXX（证据等级为C）耐药。】
def fjfy_regimen_sum(regimen_list):
	regimen_dict = {}
	for regimen in regimen_list:
		if (regimen["evi_conclusion_simple"], regimen["clinical_significance_cn"]) not in regimen_dict.keys():
			regimen_dict.setdefault((regimen["evi_conclusion_simple"], regimen["clinical_significance_cn"]), [])
		regimen_dict[(regimen["evi_conclusion_simple"], regimen["clinical_significance_cn"])].append(regimen["regimen_name"])
	
	result_sense = []
	result_resis = []
	result = []
	for level in ["A", "B", "C", "D"]:
		if (level, "敏感") in regimen_dict.keys():
			result_sense.append("{0}（证据等级为{1}）".format("、".join(regimen_dict[(level, "敏感")]), level))
		if (level, "耐药") in regimen_dict.keys():
			result_resis.append("{0}（证据等级为{1}）".format("、".join(regimen_dict[(level, "耐药")]), level))
	if result_sense:
		result.append("可能对{0}敏感".format("，".join(result_sense)))
	if result_resis:
		result.append("可能对{0}耐药".format("，".join(result_resis)))

	return "；".join(result)
jinja2.filters.FILTERS["fjfy_regimen_sum"] = fjfy_regimen_sum

# 福建附一CP200-肺癌ALK融合、EGFR I类仅展示AB等级药物-2025.11.28
def fjfy_filter_alk_egfr_ab_regimen(regimen_list):
	return [i for i in regimen_list if i["evi_conclusion_simple"] in ["A", "B"]]
jinja2.filters.FILTERS["fjfy_filter_alk_egfr_ab_regimen"] = fjfy_filter_alk_egfr_ab_regimen

# 华西HRD，检测小结变异展示为一行-2025.12.02
def schx_hrd_var_summary(var_list):
	result = []
	for var in var_list:
		if var["hgvs_p_ZJZL"] != "p.?":
			result.append(var["gene_symbol"] + " " + var["hgvs_p_ZJZL"])
		else:
			result.append(var["gene_symbol"] + " " + var["hgvs_c"])
	return "，".join(result)
jinja2.filters.FILTERS["schx_hrd_var_summary"] = schx_hrd_var_summary

# 华西HRD，BRCA变异药物特殊展示规则-2025.12.02
# 第一部分只展示 奥拉帕利+贝伐珠单抗、尼拉帕利阿比特龙片+泼尼松、尼拉帕利阿比特龙片+泼尼松龙三个药物
# 第二部分不展示上述三个药物
# 有获批适应症的展示适应症描述，没有则展示详细描述
# 描述相同的治疗方案需要合并展示
def schx_hrd_evi_sum(info):
	evi_raw = info[0]
	evi_type = info[1]
	gene_symbol = info[2]
	# 处理evi_list-只保留regimen_name、clinical_significance_cn、evi_conclusion_simple和inter
	evi_list = []
	for evi in evi_raw:
		evi_list.append(
			{
				"regimen_name" : evi["regimen_name"],
				"inter" : evi["evi_adaptation_disease_cn"] if evi["evi_adaptation_disease_cn"] else evi["evi_interpretation"],
				"clinical_significance_cn" : evi["clinical_significance_cn"],
				"evi_conclusion_simple" : evi["evi_conclusion_simple"]
			}
		)
	# 拆分为特殊治疗方案和其他
	special_regimen = ["奥拉帕利+贝伐珠单抗", "尼拉帕利阿比特龙片+泼尼松", "尼拉帕利阿比特龙片+泼尼松龙"]
	in_evi_list = [i for i in evi_list if i["regimen_name"] in special_regimen]
	out_evi_list = [i for i in evi_list if i["regimen_name"] not in special_regimen]
	# 对描述相同的治疗方案进行合并
	def merge_regimen_for_SCHX(_evi_list):
		tmp_dict = {}
		for evi in _evi_list:
			if (evi["inter"], evi["clinical_significance_cn"], evi["evi_conclusion_simple"]) not in tmp_dict.keys():
				tmp_dict.setdefault((evi["inter"], evi["clinical_significance_cn"], evi["evi_conclusion_simple"]), [])
			tmp_dict[(evi["inter"], evi["clinical_significance_cn"], evi["evi_conclusion_simple"])].append(evi["regimen_name"])
		merge_result = []
		for inter_info, regimen in tmp_dict.items():
			merge_result.append(
				{
					"regimen_name" : "{0}（{1}，{2}级）".format("、".join(regimen), inter_info[1], inter_info[2]),
					"inter" : inter_info[0]
				}
			)
		return merge_result
	if gene_symbol in ["BRCA1", "BRCA2"] and evi_type == "in":
		return merge_regimen_for_SCHX(in_evi_list)
	elif gene_symbol in ["BRCA1", "BRCA2"] and evi_type == "out":
		return merge_regimen_for_SCHX(out_evi_list)
	else:
		return merge_regimen_for_SCHX(evi_list)
jinja2.filters.FILTERS["schx_hrd_evi_sum"] = schx_hrd_evi_sum

# 华西HRD，过滤出BRCA基因变异-2025.12.02
def schx_filter_brca(var_list):
	return [var for var in var_list if var["gene_symbol"] in ["BRCA1", "BRCA2"]]
jinja2.filters.FILTERS["schx_filter_brca"] = schx_filter_brca

# 安徽省立BRCA去除C4、D4证据并返回AB或CD证据-2025.12.08
def ahsl_filter_AB_or_CD(regimen_list):
	AB = [i for i in regimen_list if i["evi_conclusion_simple"] in ["A", "B"]]
	CD = [i for i in regimen_list if i["evi_conclusion_simple"] in ["C", "D"] and i["evi_conclusion"] not in ["C4", "D4"]]
	if AB:
		return AB
	else:
		return CD
jinja2.filters.FILTERS["ahsl_filter_AB_or_CD"] = ahsl_filter_AB_or_CD

# 安徽省立BRCA去除C4、D4证据并返回A或BCD证据-2025.12.08
def ahsl_filter_A_or_BCD(regimen_list):
	A = [i for i in regimen_list if i["evi_conclusion_simple"] in ["A"]]
	# 靶向药物全都展示的话，这边就不过滤掉C4、D4了吧
	#BCD = [i for i in regimen_list if i["evi_conclusion_simple"] in ["B", "C", "D"] and i["evi_conclusion"] not in ["C4", "D4"]]
	BCD = [i for i in regimen_list if i["evi_conclusion_simple"] in ["B", "C", "D"]]
	if A:
		return A
	else:
		return BCD
jinja2.filters.FILTERS["ahsl_filter_A_or_BCD"] = ahsl_filter_A_or_BCD

# 安徽省立BRCA去除C4、D4证据-2025.12.11
def ahsl_filter_C4_D4(regimen_list):
	return [i for i in regimen_list if i["evi_conclusion"] not in ["C4", "D4"]]
jinja2.filters.FILTERS["ahsl_filter_C4_D4"] = ahsl_filter_C4_D4

# 安徽省立BRCA-治疗方案返回所有证据，需要进行筛选-2025.12.08
# 治疗方案对应多条证据（等级相同）
def ahsl_brca_summary_regimen(raw_regimen_list):
	# 对证据进行排序
	appr_rule = ["NMPA_FDA", "NMPA", "FDA", "CSCO_NCCN", "CSCO", "NCCN", \
			   "Clinical-phase I", "Clinical-phase II", "Clinical-phase III", "Clinical-phase IV", "Clinical-retrospective", "Clinical-unknown phase", \
				"Case report", "Preclinical-in vitro", "Preclinical-in vivo"]
	for evi in raw_regimen_list:
		evi["ahsl_appr_sort_key"] = evi["refer_agency"] if evi["refer_agency"] else evi["evidence_level"]
		#print (evi["regimen_name"], evi["ahsl_appr_sort_key"])
	regimen_list = sorted(raw_regimen_list, key = lambda i : appr_rule.index(i["ahsl_appr_sort_key"]))

	regimen_appr_list = {}
	for evi in raw_regimen_list:
		if evi["regimen_name"] not in regimen_appr_list.keys():
			regimen_appr_list.setdefault(evi["regimen_name"], [])
		if evi["refer_agency"]:
			regimen_appr_list[evi["regimen_name"]].append(evi["refer_agency"])
	# 二次排序 NMPA/FDA > NMPA > FDA > CSCO/NCCN > CSCO > NCCN >其他
	for evi in raw_regimen_list:
		if evi["regimen_name"] in regimen_appr_list.keys() and regimen_appr_list[evi["regimen_name"]]:
			if "NMPA" in regimen_appr_list[evi["regimen_name"]] and "FDA" in regimen_appr_list[evi["regimen_name"]]:
				evi["ahsl_appr_sort_key_2"] = 0
			elif "NMPA" in regimen_appr_list[evi["regimen_name"]]:
				evi["ahsl_appr_sort_key_2"] = 1
			elif "FDA" in regimen_appr_list[evi["regimen_name"]]:
				evi["ahsl_appr_sort_key_2"] = 2
			elif "CSCO" in regimen_appr_list[evi["regimen_name"]] and "NCCN" in regimen_appr_list[evi["regimen_name"]]:
				evi["ahsl_appr_sort_key_2"] = 3
			elif "CSCO" in regimen_appr_list[evi["regimen_name"]]:
				evi["ahsl_appr_sort_key_2"] = 4
			elif "NCCN" in regimen_appr_list[evi["regimen_name"]]:
				evi["ahsl_appr_sort_key_2"] = 5
			else:
				evi["ahsl_appr_sort_key_2"] = 6
		else:
			evi["ahsl_appr_sort_key_2"] = 6
	regimen_list = sorted(raw_regimen_list, key = lambda i : i["ahsl_appr_sort_key_2"])

	# 1. 根据获批机构进行分类，NMPA/FDA、NMPA、FDA、CSCO/NCCN、CSCO、NCCN、临床试验/回顾性分析、案例报道/临床前研究
	regimen_dict = {}
	for evi in regimen_list:
		evi["inter"] = evi["evi_adaptation_disease_cn"] if evi["evi_adaptation_disease_cn"] and evi["evi_conclusion_simple"] == "A" else evi["evi_interpretation"]
		if evi["regimen_name"] not in regimen_dict.keys():
			regimen_dict.setdefault(evi["regimen_name"], {})
		key = evi["refer_agency"] if evi["refer_agency"] and evi["evi_conclusion_simple"] == "A" else evi["evidence_level"]
		if key not in regimen_dict[evi["regimen_name"]].keys():
			regimen_dict[evi["regimen_name"]].setdefault(key, [])
		regimen_dict[evi["regimen_name"]][key].append(evi)
	
	# 2. 对NMPA/FDA、CSCO/NCCN、临床试验/回顾性分析、案例报道/临床前研究进行合并
	def ahsh_stran_evi(group_evi_list):
		# 描述可能重复，这边加个去重的操作
		inter_list = []
		for i in group_evi_list:
			if i["inter"] not in inter_list:
				inter_list.append(i["inter"])
		result = {}
		result["inter"] = "\n".join(inter_list)
		result["evi_conclusion_simple"] = group_evi_list[0]["evi_conclusion_simple"] if group_evi_list else ""
		return result if group_evi_list else {}
	
	clinic_phase = ["Clinical-phase I", "Clinical-phase II", "Clinical-phase III", "Clinical-phase IV", "Clinical-retrospective", "Clinical-unknown phase"]
	other = ["Case report", "Preclinical-in vitro", "Preclinical-in vivo"]
	regimen_dict_merge = {}
	#regimen_appr_list = {}
	for regimen, info in regimen_dict.items():
		#regimen_appr_list[regimen] = []
		tmp_dict = {
			"NMPA_FDA" : {},
			"NMPA" : {},
			"FDA" : {},
			"CSCO_NCCN" : {},
			"CSCO" : {},
			"NCCN" : {},
			"clinic" : {},
			"other" : {}
		}
		if "NMPA" in info.keys() and "FDA" in info.keys():
			tmp_dict["NMPA_FDA"] = ahsh_stran_evi(info["NMPA"] + info["FDA"])
			#regimen_appr_list[regimen].append("NMPA_FDA")
		elif "NMPA" in info.keys():
			tmp_dict["NMPA"] = ahsh_stran_evi(info["NMPA"])
			#regimen_appr_list[regimen].append("NMPA")
		elif "FDA" in info.keys():
			tmp_dict["FDA"] = ahsh_stran_evi(info["FDA"])
			#regimen_appr_list[regimen].append("FDA")

		if "CSCO" in info.keys() and "NCCN" in info.keys():
			tmp_dict["CSCO_NCCN"] = ahsh_stran_evi(info["CSCO"] + info["NCCN"])
			#regimen_appr_list[regimen].append("CSCO_NCCN")
		elif "CSCO" in info.keys():
			tmp_dict["CSCO"] = ahsh_stran_evi(info["CSCO"])
			#regimen_appr_list[regimen].append("CSCO")
		elif "NCCN" in info.keys():
			tmp_dict["NCCN"] = ahsh_stran_evi(info["NCCN"])
			#regimen_appr_list[regimen].append("NCCN")

		ahsl_regimen_evi_sort_rule = ["NMPA", "FDA", "CSCO", "NCCN", "Clinical-phase I", "Clinical-phase II", "Clinical-phase III", "Clinical-phase IV", "Clinical-retrospective", "Clinical-unknown phase", "Case report", "Preclinical-in vitro", "Preclinical-in vivo"]
		tmp_clinic_list = []
		for clinic in clinic_phase:
			if clinic in info.keys():
				tmp_clinic_list.extend(info[clinic])
		if tmp_clinic_list:
			#print ("before", tmp_clinic_list)
			# 排序： NMPA > FDA > CSCO > NCCN >其他
			for i in tmp_clinic_list:	
				i["ahsl_regimen_evi_sort_key"] = i["refer_agency"] if i["refer_agency"] else i["evidence_level"]
			tmp_clinic_list = sorted(tmp_clinic_list, key = lambda i : ahsl_regimen_evi_sort_rule.index(i["ahsl_regimen_evi_sort_key"]))
			#print ("after", tmp_clinic_list)
			tmp_dict["clinic"] = ahsh_stran_evi(tmp_clinic_list)

		tmp_other_item = []
		for other_item in other:
			if other_item in info.keys():
				tmp_other_item.extend(info[other_item])
		if tmp_other_item:
			for i in tmp_other_item:
				i["ahsl_regimen_evi_sort_key"] = i["refer_agency"] if i["refer_agency"] else i["evidence_level"]
			tmp_other_item = sorted(tmp_other_item, key = lambda i : ahsl_regimen_evi_sort_rule.index(i["ahsl_regimen_evi_sort_key"]))
			tmp_dict["other"] = ahsh_stran_evi(tmp_other_item)
		regimen_dict_merge[regimen] = tmp_dict
	
	# 3. 格式调整
	regimen_dict_merge_stran = {}
	for key in ["NMPA_FDA", "NMPA", "FDA", "CSCO_NCCN", "CSCO", "NCCN", "clinic", "other"]:
		if key not in regimen_dict_merge_stran.keys():
			regimen_dict_merge_stran.setdefault(key, [])
		for regimen, info in regimen_dict_merge.items():
			if info[key]:
				regimen_dict_merge_stran[key].append({
					"regimen_name" : regimen,
					"inter" : info[key]["inter"],
					"evi_conclusion_simple" : info[key]["evi_conclusion_simple"],
					"appr_list" : regimen_appr_list.get(regimen, [])
				})
	return regimen_dict_merge_stran
jinja2.filters.FILTERS["ahsl_brca_summary_regimen"] = ahsl_brca_summary_regimen

# 安徽省立BRCA对描述+等级相同的治疗方案进行合并展示-2025.12.08
def ahsl_brca_inter_merge(ahsl_inter_list):
	tmp_dict = {}
	for evi in ahsl_inter_list:
		if (evi["inter"], evi["evi_conclusion_simple"]) not in tmp_dict.keys():
			tmp_dict.setdefault((evi["inter"], evi["evi_conclusion_simple"]), [])
		tmp_dict[(evi["inter"], evi["evi_conclusion_simple"])].append(evi["regimen_name"])
	merge_result = []
	for k, v in tmp_dict.items():
		merge_result.append(
			{
				"regimen_name" : "、".join(v),
				"inter" : k[0],
				"evi_conclusion_simple" : k[1]
			}
		)
	return merge_result
jinja2.filters.FILTERS["ahsl_brca_inter_merge"] = ahsl_brca_inter_merge	

# 中南湘雅BRCA-健康人需要展示遗传风险-2025.12.12
def znxy_brca_predisposing(var_list):
	evi_result = []
	for var in var_list:
		if var["evi_sum"] and "evi_split" in var["evi_sum"].keys() and var["evi_sum"]["evi_split"] and "Predisposing" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Predisposing"]:
			for evi in var["evi_sum"]["evi_split"]["Predisposing"]:
				if evi["evi_interpretation"] not in evi_result:
					evi_result.append(evi["evi_interpretation"])
	return evi_result
jinja2.filters.FILTERS["znxy_brca_predisposing"] = znxy_brca_predisposing

# 重庆西南116-综合评估需要结合除肿瘤细胞含量外，其他所有质控项-2025.12.16
# 全部合格-良好
# 一项不合格-合格
# 两项及以上不合格-风险
def cqxn_116_qc_judge(info):
	lib_qc = info[0]
	ngs_qc = info[1]
	false_item = 0
	# 湿实验，没有值就不管了，报告里也不提示
	if "lib_dna_qc" in lib_qc.keys() and lib_qc["lib_dna_qc"]:
		if "dna_qty" in lib_qc["lib_dna_qc"].keys() and lib_qc["lib_dna_qc"]["dna_qty"]:
			if not (is_number(lib_qc["lib_dna_qc"]["dna_qty"]) and float(lib_qc["lib_dna_qc"]["dna_qty"]) >= 150):
				false_item += 1
		if "library_qty" in lib_qc["lib_dna_qc"].keys() and lib_qc["lib_dna_qc"]["library_qty"]:
			if not (is_number(lib_qc["lib_dna_qc"]["library_qty"]) and float(lib_qc["lib_dna_qc"]["library_qty"]) >= 500):
				false_item += 1
	# ngs质控结果
	if ngs_qc["dna_data_qc"]["totaldata"] and "G" in ngs_qc["dna_data_qc"]["totaldata"]:
		if float(ngs_qc["dna_data_qc"]["totaldata"].replace("G", "")) <= 1:
			false_item += 1
	elif ngs_qc["dna_data_qc"]["totaldata"] and "M" in ngs_qc["dna_data_qc"]["totaldata"]:
		if float(ngs_qc["dna_data_qc"]["totaldata"].replace("M", "")) <= 1000:
			false_item += 1
	else:
		false_item += 1

	if ngs_qc["dna_data_qc"]["cleandata_q30_num"] <= 0.75:
		false_item += 1
	if ngs_qc["dna_data_qc"]["cover_ratio_num"] < 0.95:
		false_item += 1
	if ngs_qc["dna_data_qc"]["uni20_uniq_num"] <= 0.9:
		false_item += 1
	if float(ngs_qc["qc_gradient"]["coverage_ratio_uniq_hot_180"].replace("%", "")) < 95:
		false_item += 1
	if float(ngs_qc["dna_data_qc"]["cnv_cv"]) >= 0.4:
		false_item += 1
	if float(ngs_qc["dna_data_qc"]["cnv_uni"]) >= 1.5:
		false_item += 1
	if ngs_qc["dna_data_qc"]["depth_mean_raw_num"] < 1500:
		false_item += 1
	if ngs_qc["dna_data_qc"]["depth_mean_uniq_num"] < 500:
		false_item += 1	
	if ngs_qc["dna_data_qc"]["mapping_ratio_num"] < 0.95:
		false_item += 1

	if false_item == 0:
		return "良好"
	elif false_item == 1:
		return "合格"
	else:
		return "风险"
jinja2.filters.FILTERS["cqxn_116_qc_judge"] = cqxn_116_qc_judge

# 中山人民HRR首页BRCA结果需要合并，未检出变异的基因放空行-2025.12.19
def szrm_hrr_sum(var_list):
	result = []
	result = copy.deepcopy(var_list)
	detect_gene = [a["gene_symbol"] for a in var_list]
	for gene in ["BRCA1", "BRCA2"]:
		if gene not in detect_gene:
			result.append({
				"gene_symbol" : gene
			})
	return sorted(result, key = lambda i : i["gene_symbol"])
jinja2.filters.FILTERS["szrm_hrr_sum"] = szrm_hrr_sum

# 武汉同济Master，需要过滤KNB相关结果-2025.12.23
# 仅限制在snvindel中
def whtj_mp_filter_knb(info):
	var_list = info[0]
	var_type = info[1]
	knb_result = []
	other_result = []
	for var in var_list:
		if var["bio_category"] == "Snvindel" and (var["gene_symbol"] in ["KRAS", "NRAS"] or (var["gene_symbol"] in ["BRAF"] and var["hgvs_p"] == "p.(V600E)")):
			knb_result.append(var)
		else:
			other_result.append(var)
	knb_result_gene = [var["gene_symbol"] for var in knb_result]
	for gene in ["KRAS", "NRAS", "BRAF"]:
		if gene not in knb_result_gene:
			knb_result.append({
				"gene_symbol" : gene
			})

	if var_type == "inknb":
		return knb_result
	else:
		return other_result
jinja2.filters.FILTERS["whtj_mp_filter_knb"] = whtj_mp_filter_knb

# 2025.12.24-同济Master-免疫正负相关只展示有检出的，并且变异展示格式跟前面的一样
def whtj_mp_io(var_list):
	io_result = []
	# 汇总体细胞I/II/肿瘤发生发展相关变异+胚系致病/疑似致病性变异
	io_gene_P = ["ATM","ATR","BRCA1","BRCA2","BRIP1","CHEK1","CHEK2","ERCC1","FANCA","MRE11",\
				 "PALB2","RAD50","XRCC1","MLH1","MSH2","MSH6","PMS2","POLE","POLD1","TP53",\
				 "KRAS","CD274","ARID1A","LRP1B","SETD2","PRKDC","TERT","KMT2D","FAT1","CDK12",\
				 "SERPINB3","SERPINB4"]
	io_gene_N = ["EGFR","ALK","MDM2","MDM4","CDKN2A","CDKN2B","DNMT3A","STK11","IFNGR1","IRF1",\
				 "JAK1","JAK2","APC","CTNNB1","B2M","PTEN","CCND1","FGF3","FGF19"]
	cnv_gene_list = ["CD274", "MDM2", "MDM4", "CCND1", "FGF3", "FGF19"]
	# hd 是额外需要展示的，原有的输出类型不变
	hd_gene_list = ["ATM", "BRCA1", "BRCA2", "BRIP1", "CDK12", "CHEK1", "CHEK2", "FANCA", \
				 	"PALB2", "SETD2", "TP53", "CDKN2A", "CDKN2B", "PTEN", "STK11"]

	for var in var_list:
		# 仅展示扩增的基因
		if var["bio_category"] == "Cnv" and var["gene_symbol"] in cnv_gene_list:
			io_result.append(var)
		# 仅展示融合的基因
		elif var["bio_category"] in ["Sv", "PSeqRnaSv"] and set(re.split(",", var["gene_symbol"])) & set(["ALK"]):
			io_result.append(var)
		# HD基因经确认同时展示snvindel和hd
		elif (var["bio_category"] == "Snvindel" or var["bio_category"] == "PHd") and var["gene_symbol"] in hd_gene_list:
			io_result.append(var)
		# 其余基因展示Snvindel
		elif var["bio_category"] == "Snvindel" and var["gene_symbol"] in io_gene_P + io_gene_N:
			io_result.append(var)

	io_p_result = [var for var in io_result if var["gene_symbol"] in io_gene_P]
	for var in io_p_result:
		var["io_type"] = "P"
	io_n_result = [var for var in io_result if var["gene_symbol"] in io_gene_N and not (var["bio_category"] == "Cnv" and var["gene_symbol"] in ["CCND1","FGF3","FGF19"])]
	for var in io_n_result:
		var["io_type"] = "N"
	# 判断共突变
	io_cnv_result_gene = [var["gene_symbol"] for var in io_result if var["bio_category"] == "Cnv"]
	if "CCND1" in io_cnv_result_gene and "FGF3" in io_cnv_result_gene and "FGF19" in io_cnv_result_gene:
		io_n_result.append(
			{
				"gene_symbol" : "CCND1/FGF3/FGF19",
				"bio_category" : "Cnv",
				"io_type" : "N"
			}
		)
	return io_p_result + io_n_result
jinja2.filters.FILTERS["whtj_mp_io"] = whtj_mp_io

# 武汉同济-MP和CP200-治疗方案返回最高等级的所有证据，整合证据引用机构的信息（仅A级）-2025.12.26
def whtj_mp_regimen_list(a):
	result = []
	# 用药相关
	Predictive_list = a["evi_sum"]["evi_split"]["Predictive"] if "evi_sum" in a.keys() and a["evi_sum"] and "evi_split" in  a["evi_sum"].keys() and a["evi_sum"]["evi_split"] and "Predictive" in a["evi_sum"]["evi_split"].keys() and a["evi_sum"]["evi_split"]["Predictive"] else []
	Predictive_stran_dict = {}
	for evi in Predictive_list:
		key = (evi["regimen_name"], evi["evi_conclusion_simple"], evi["clinical_significance_cn"])
		if key not in Predictive_stran_dict.keys():
			Predictive_stran_dict.setdefault(key, [])
		if evi["refer_agency"] and evi["refer_agency"] not in Predictive_stran_dict[key]:
			Predictive_stran_dict[key].append(evi["refer_agency"])
	for k, v in Predictive_stran_dict.items():
		if k[1] == "A" and v:
			result.append("{0}（{1}，{2}级，{3}）".format(k[0], k[2], k[1], "/".join(v)))
		else:
			result.append("{0}（{1}，{2}级）".format(k[0], k[2], k[1]))
	# 预后相关
	#Prognostic_list = a["evi_sum"]["evi_split"]["Prognostic"] if "Prognostic" in a["evi_sum"]["evi_split"].keys() else []
	Prognostic_list = a["evi_sum"]["evi_split"]["Prognostic"] if "evi_sum" in a.keys() and a["evi_sum"] and "evi_split" in  a["evi_sum"].keys() and a["evi_sum"]["evi_split"] and "Prognostic" in a["evi_sum"]["evi_split"].keys() and a["evi_sum"]["evi_split"]["Prognostic"] else []
	Prognostic_stran_dict = {}
	for evi in Prognostic_list:
		key = (evi["clinical_significance_cn"], evi["evi_conclusion_simple"])
		if key not in Prognostic_stran_dict.keys():
			Prognostic_stran_dict.setdefault(key, [])
		if evi["refer_agency"] and evi["refer_agency"] not in Prognostic_stran_dict[key]:
			Prognostic_stran_dict[key].append(evi["refer_agency"])
	for k, v in Prognostic_stran_dict.items():
		if k[1] == "A" and v:
			result.append("{0}（{1}，{2}级，{3}）".format("预后"+k[0], " / ", k[1], "/".join(v)))
		else:
			result.append("{0}（{1}，{2}级）".format("预后"+k[0], " / ", k[1]))
	# 辅助诊断相关
	#Diagnostic_list = a["evi_sum"]["evi_split"]["Diagnostic"] if "Diagnostic" in a["evi_sum"]["evi_split"].keys() else []
	Diagnostic_list = a["evi_sum"]["evi_split"]["Diagnostic"] if "evi_sum" in a.keys() and a["evi_sum"] and "evi_split" in  a["evi_sum"].keys() and a["evi_sum"]["evi_split"] and "Diagnostic" in a["evi_sum"]["evi_split"].keys() and a["evi_sum"]["evi_split"]["Diagnostic"] else []
	Diagnostic_stran_dict = {}
	for evi in Diagnostic_list:
		key = evi["evi_conclusion_simple"]
		if key not in Diagnostic_stran_dict.keys():
			Diagnostic_stran_dict.setdefault(key, [])
		if evi["refer_agency"] and evi["refer_agency"] not in Diagnostic_stran_dict[key]:
			Diagnostic_stran_dict[key].append(evi["refer_agency"])
	for k, v in Diagnostic_stran_dict.items():
		if k == "A" and v:
			result.append("{0}（{1}，{2}级，{3}）".format("辅助诊断", " / ", k, "/".join(v)))
		else:
			result.append("{0}（{1}，{2}级）".format("辅助诊断", " / ", k))
	return "、".join(result) if result else "-"
jinja2.filters.FILTERS["whtj_mp_regimen_list"] = whtj_mp_regimen_list

# 同济Master-胚系遗传风险需要去重-2025.12.26
def whtj_mp_evi_inter_redup(evi_list):
	result = []
	for evi in evi_list:
		if evi["evi_interpretation"] not in result:
			result.append(evi["evi_interpretation"])
	return result
jinja2.filters.FILTERS["whtj_mp_evi_inter_redup"] = whtj_mp_evi_inter_redup

# 同济Master/CP200-参考文献-evi_sum-2025.12.29
def whtj_evi_refer(info):
	literature_evi_sum_list = []
	var_list = info[0]
	refer_original = info[1]
	# 1. refer_original处理为字段，key 为refer_code
	refer_original_dict = {}
	for refer in refer_original:
		if "refer_code" in refer.keys():
			refer_original_dict[refer["refer_code"]] = refer

	# 2. 获取变异evi_sum中的refer_code合集（这边遗传相关的一起加上）
	for var in var_list:
		if var:
			for evi_type in ["Predictive", "Prognostic", "Diagnostic", "Predisposing"]:
				if evi_type in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"][evi_type]:
					for evi in var["evi_sum"]["evi_split"][evi_type]:
						if "literature_evi_sum" in evi.keys() and evi["literature_evi_sum"]:
							literature_evi_sum_list.extend(evi["literature_evi_sum"])
	
	# 3. 变异evi_sum refer_code合集获取具体参考文献信息
	result = []
	for refer_code in literature_evi_sum_list:
		if refer_code in refer_original_dict.keys():
			refer_info = refer_original_dict[refer_code]
			# 存在PMID
			if refer_info["pmid"]:
				pmid = refer_info["pmid"] if refer_info["pmid"] else ""
				authors = refer_info["authors"] if refer_info["authors"] else ""
				date = refer_info["date"] if refer_info["date"] else ""
				title = refer_info["title"] if refer_info["title"] else ""
				journal = refer_info["journal"] if refer_info["journal"] else ""
				vol = refer_info["vol"] if refer_info["vol"] else ""
				pmid_info = " ".join([authors, date, title, journal, vol, pmid])
				if pmid_info not in result:
					result.append(pmid_info)
			# 不存在PMID
			else:
				if refer_info["title"] not in result:
					result.append(refer_info["title"])
	return result
jinja2.filters.FILTERS["whtj_evi_refer"] = whtj_evi_refer


# 浙江人民HRR，检测结果总结需要列出治疗方案，<=3个治疗方案按实际情况写，>3个的只取前三个-2026.01.04
# PIK3CA 伊那利塞+哌柏西利+氟维司群、卡匹色替+氟维司群、阿培利司+氟维司群 三个优先排序
def zjrm_hrr_sum_v2(var_list):
	pik3ca_regimen_dict = {"伊那利塞+哌柏西利+氟维司群" : 0, "卡匹色替+氟维司群" : 1, "阿培利司+氟维司群" : 2}
	result = []
	for var in var_list:
		regimen_list = []
		if "evi_split" in var["evi_sum"].keys():
			if "Predictive" in var["evi_sum"]["evi_split"].keys():
				for i in var["evi_sum"]["evi_split"]["Predictive"]:
					i["zjrm_pik3ca_sort_key"] = pik3ca_regimen_dict.get(i["regimen_name"], 999)
				var["evi_sum"]["evi_split"]["Predictive"] = sorted(var["evi_sum"]["evi_split"]["Predictive"], key = lambda i : (i["evi_conclusion_simple"], i["zjrm_pik3ca_sort_key"]))
				for i in var["evi_sum"]["evi_split"]["Predictive"]:
					regimen_list.append("{0}（{1}，{2}级）".format(i["regimen_name"], i["clinical_significance_cn"], i["evi_conclusion_simple"]))
		if len(regimen_list) == 0:
			var["zjrm_hrr_sum_str"] = ""
		elif len(regimen_list) <= 3:
			var["zjrm_hrr_sum_str"] = "推荐{0}。".format("、".join(regimen_list))
		else:
			var["zjrm_hrr_sum_str"] = "推荐{0}等，详见第四部分详细检测结果。".format("、".join(regimen_list[0:3]))
		if regimen_list:
			result.append(var)
	return result
jinja2.filters.FILTERS["zjrm_hrr_sum_v2"] = zjrm_hrr_sum_v2


# 浙江人民HRR PIK3CA治疗方案排序-2026.01.04
# 伊那利塞+哌柏西利+氟维司群、卡匹色替+氟维司群、阿培利司+氟维司群 三个优先排序
def zjrm_hrr_pik3ca_sort(info):
	regimen_list = info[0]
	gene = info[1]
	pik3ca_regimen_dict = {"伊那利塞+哌柏西利+氟维司群" : 0, "卡匹色替+氟维司群" : 1, "阿培利司+氟维司群" : 2}
	if gene == "PIK3CA":
		for evi in regimen_list:
			evi["zjrm_pik3ca_sort_key"] = pik3ca_regimen_dict.get(evi["regimen_name"], 999)
		return sorted(regimen_list, key = lambda i : (i["evi_conclusion_simple"], i["zjrm_pik3ca_sort_key"]))
	else:
		return regimen_list
jinja2.filters.FILTERS["zjrm_hrr_pik3ca_sort"] = zjrm_hrr_pik3ca_sort

# 中山人民检测结果汇总-2026.01.04
def zsrm_hrr_hrd_summary(var_list):
	result = []
	for var in var_list:
		significance_clinic = ""
		if var["type"] in ["Loss", "Gain"] or var["var_origin"] == "germline":
			significance_clinic = "致病性变异" if var["clinic_num_g"] == 5 else "疑似致病性变异"
		else:
			significance_clinic = "I类变异" if var["clinic_num_s"] == 5 else "II类变异"
		if var["type"] == "Loss":
			result.append("{0} {1} del {2}".format(var["gene_symbol"], var["value"], significance_clinic))
		elif var["type"] == "Gain":
			result.append("{0} {1} dup {2}".format(var["gene_symbol"], var["value"], significance_clinic))
		else:
			if var["hgvs_p"] != "p.?":
				result.append("{0} {1} {2} {3} {4}".format(var["gene_symbol"], var["gene_region"], var["hgvs_c"], var["hgvs_p"], significance_clinic))
			else:
				result.append("{0} {1} {2} {3}".format(var["gene_symbol"], var["gene_region"], var["hgvs_c"], significance_clinic))
	return "；".join(result)
jinja2.filters.FILTERS["zsrm_hrr_hrd_summary"] = zsrm_hrr_hrd_summary

# BHD-治疗方案需要进行分类，伴随诊断药物和其他PARP抑制剂-2026.01.12
def bhd_regimen_split(info):
	regimen_list = info[0]
	appr_type = info[1]
	appr_list = []
	other_list = []
	for regimen in regimen_list:
		if regimen["regimen_name"] in ["尼拉帕利阿比特龙片+泼尼松", "尼拉帕利阿比特龙片+泼尼松龙"]:
			appr_list.append(regimen)
		else:
			other_list.append(regimen)
	if appr_type == "cd":
		return appr_list
	else:
		return other_list
jinja2.filters.FILTERS["bhd_regimen_split"] = bhd_regimen_split

# 复旦中山MP-过滤指定基因变异结果-2026.01.16
def fdzs_filter_genelist(info):
	gene_list = info[0]
	var_list = info[1]
	return [var for var in var_list if set(re.split(",", var["gene_symbol"])) & set(gene_list)]
jinja2.filters.FILTERS["fdzs_filter_genelist"] = fdzs_filter_genelist

# 安徽省立BRCA-汇总4/5类变异-2026.01.19
# BRCAX基因的c.（p.）变异为胚系致病性变异
def ahsl_brca_summary(var_list):
	result = []
	for var in var_list:
		significance = "致病性变异" if var["clinic_num_g"] == 5 else "疑似致病性变异"
		if var["type"] == "Loss":
			result.append("{0}基因的{1} del为胚系{2}".format(var["gene_symbol"], var["value"], significance))
		elif var["type"] == "Gain":
			result.append("{0}基因的{1} dup为胚系{2}".format(var["gene_symbol"], var["value"], significance))
		else:
			if var["hgvs_p"] != "p.?":
				result.append("{0}基因的{1}（{2}）变异为胚系{3}".format(var["gene_symbol"], var["hgvs_c"], var["hgvs_p"], significance))
			else:
				result.append("{0}基因的{1}变异为胚系{2}".format(var["gene_symbol"], var["hgvs_c"], significance))
	return "、".join(result)
jinja2.filters.FILTERS["ahsl_brca_summary"] = ahsl_brca_summary

# 温附一CP43-判断变异基因是否在非10基因获批中-2026.01.23
def wfy_cp43_judge_var_inothergene(var):
	lc10_gene = ["ALK", "BRAF", "EGFR", "ERBB2", "KRAS", "MET", "NRAS", "PIK3CA", "RET", "ROS1"]
	if set(re.split(",", var["gene_symbol"])) & set(lc10_gene):
		return False
	else:
		return True
jinja2.filters.FILTERS["wfy_cp43_judge_var_inothergene"] = wfy_cp43_judge_var_inothergene

# 温附一CP43-返回非10基因变异-2026.01.23
def wfy_cp43_get_othergene(var_list):
	result = []
	for var in var_list:
		if wfy_cp43_judge_var_inothergene(var):
			result.append(var)
	return result
jinja2.filters.FILTERS["wfy_cp43_get_othergene"] = wfy_cp43_get_othergene

# 上海肺科116新增“本次检测结果”-2026.01.27
def shfk_116_sum(var_list):
	type_stran = {
		"3'UTR" : "3'UTR区突变",
		"5'UTR" : "5'UTR区突变",
		"Intronic" : "内含子区突变",
		"FlankingRegion3" : "侧翼区突变",
		"FlankingRegion5" : "侧翼区突变"
	}
	result = []
	for var in var_list:
		var_str = ""
		if var["bio_category"] == "Snvindel":
			gene_region_cn = gene_region_strn(var["gene_region"])
			type_cn = var["type_cn"] if var["type_cn"] != "--" else type_stran.get(var["type"], var["type"])
			if var["hgvs_p"] != "p.?":
				var_str = "{0}基因{1}{2}{3}:{4}，突变丰度{5}".format(var["gene_symbol"], gene_region_cn, type_cn, var["hgvs_c"], var["hgvs_p"], var["freq_str"])
			else:
				var_str = "{0}基因{1}{2}{3}，突变丰度{4}".format(var["gene_symbol"], gene_region_cn, type_cn, var["hgvs_c"], var["freq_str"])
		elif var["bio_category"] == "Cnv":
			var_str = "{0}基因扩增，拷贝数{1}".format(var["gene_symbol"], var["cn_mean"])
		elif var["bio_category"] == "Sv":
			# EML4-ALK融合需要展示具体亚型-2024.07.18
			if var["five_prime_gene"] == "EML4" and var["three_prime_gene"] == "ALK":
				var_desc_merge = shfk_cp_alksv([var["five_prime_gene"]+":"+var["five_prime_cds"]+"-"+var["three_prime_gene"]+":"+var["three_prime_cds"]])
			else:
				var_desc_merge = var["five_prime_gene"]+":"+var["five_prime_cds"]+"-"+var["three_prime_gene"]+":"+var["three_prime_cds"]
			var_str = "{0}-{1}融合，具体的融合型为{2}，突变丰度{3}".format(var["five_prime_gene"], var["three_prime_gene"], var_desc_merge, var["freq_str"])
		result.append(var_str)
	# 最后一个结尾加句号，其他结尾加分号
	result_ap = []
	if len(result) > 1:
		for var in result[0:-1]:
			var = var + "；"
			result_ap.append(var)
	if len(result) >= 1:
		result_ap.append(result[-1] + "。")
	return result_ap
jinja2.filters.FILTERS["shfk_116_sum"] = shfk_116_sum

def shfk_116_sv(var_list):
	'''
	上海肺科CP40：检测到融合时，其他说明中要写“本次检测到gene1-gene2融合。具体的融合型为gene1:exon-gene2:exon”
	'''
	sv_var = [var for var in var_list if var["bio_category"] == "Sv"]
	sv_str_list = []
	for var in sv_var:
		var_desc_merge = var["five_prime_gene"]+":"+var["five_prime_cds"]+"-"+var["three_prime_gene"]+":"+var["three_prime_cds"]
		sv_str_list.append("检测到"+var["five_prime_gene"]+"-"+var["three_prime_gene"]+"融合，具体的融合型为"+var_desc_merge)	
	sv_str = "本次"+"；".join(sv_str_list)+"。" if sv_str_list else ""
	return sv_str
jinja2.filters.FILTERS["shfk_116_sv"] = shfk_116_sv

# 德阳人民150，药物介绍过滤出存在非BRCA基因变异的结果-2026.01.27
def dyrm_filter_drug_other(drug_list):
	for drug in drug_list:
		drug["judge_other_var"] = False
		for a in drug["var"]:
			if ("gene_symbol" in a.keys() and a["gene_symbol"] and a["gene_symbol"] not in ["BRCA1", "BRCA2"] and "cnv_type" not in a.keys()) or "biomarker_type" in a.keys() or "cnv_type" in a.keys():
				drug["judge_other_var"] = True
	return [a for a in drug_list if a["judge_other_var"]]
jinja2.filters.FILTERS["dyrm_filter_drug_other"] = dyrm_filter_drug_other


# 德阳人民150，药物介绍过滤出存在BRCA snvindel基因变异的结果-2026.01.28
def dyrm_filter_drug_brca_snvindel(drug_list):
	for drug in drug_list:
		drug["judge_brca_snvindel"] = False
		for a in drug["var"]:
			if "gene_symbol" in a.keys() and a["gene_symbol"] and a["gene_symbol"] in ["BRCA1", "BRCA2"] and "cnv_type" not in a.keys():
				drug["judge_brca_snvindel"] = True
	return [a for a in drug_list if a["judge_brca_snvindel"]]
jinja2.filters.FILTERS["dyrm_filter_drug_brca_snvindel"] = dyrm_filter_drug_brca_snvindel

# 日期处理，输入2026-01-01，输出2026年01月01日-2026.01.29
def date_stran(indate):
	date_obj = datetime.strptime(indate, "%Y-%m-%d")
	date_str = date_obj.strftime("%Y年%m月%d日")
	return date_str
jinja2.filters.FILTERS["date_stran"] = date_stran

# 2026.01.29-用来解决for循环必须有空白行的问题
# 列表最后一个元素标记上laster，模板中判断不存在laster则加空白行
def mark_laster(_list):
	if _list:
		_list[-1]["laster"] = True
	return _list
jinja2.filters.FILTERS["mark_laster"] = mark_laster

# 2026.02.03-武汉协和MP小结展示具体变异-体细胞变异-氨基酸展示三字母
def whxh_var_sum_s(var_list):
	result = []
	for var in var_list:
		if var["bio_category"] == "Snvindel":
			if var["hgvs_p"] != "p.?":
				result.append(var["gene_symbol"]+" "+var["hgvs_p_abbr"])
			else:
				result.append(var["gene_symbol"]+" "+var["hgvs_c"])
		elif var["bio_category"] == "Cnv":
			# 2026.01.12-若cnv_type为Loss/loss，返回缺失
			if "cnv_type" in var.keys() and var["cnv_type"] and var["cnv_type"] in ["Loss", "loss"]:
				result.append(var["gene_symbol"]+" 缺失")
			else:
				result.append(var["gene_symbol"]+" 扩增")
		elif var["bio_category"] in ["Sv", "PSeqRnaSv"]:
			if var["five_prime_gene"] == "MET" and var["three_prime_gene"] == "MET":
				result.append("MET exon14 跳跃")
			else:
				# 融合可能会有重复（rna exon相同，断点不同的情况）
				if var["five_prime_gene"]+"-"+var["three_prime_gene"]+" 融合" not in result:
					result.append(var["five_prime_gene"]+"-"+var["three_prime_gene"]+" 融合")
		# 2025.07.15-新增PHd
		elif var["bio_category"] == "PHd":
			if var["type"] == "HomoDel":
				if var["gene_symbol"]+" 纯合缺失" not in result:
					result.append(var["gene_symbol"]+" 纯合缺失")
			elif var["type"] == "HeteDel":
				if var["gene_symbol"]+" 杂合缺失" not in result:
					result.append(var["gene_symbol"]+" 杂合缺失")
			else:
				if var["gene_symbol"]+" 未知变异类型！" not in result:
					result.append(var["gene_symbol"]+" 未知变异类型！")
		# 2025.07.15-新增完成
	return ", ".join(result)
jinja2.filters.FILTERS["whxh_var_sum_s"] = whxh_var_sum_s

# 2026.02.03-武汉协和MP小结展示具体变异-胚系变异-氨基酸展示三字母
def whxh_var_sum_g(var_list):
	result = []
	for var in var_list:
		hgvs = var["hgvs_c"] if var["hgvs_p"] == "p.?" else var["hgvs_p_abbr"]
		clinic_g_cn = "致病变异" if var["clinic_num_g"] == 5 else "疑似致病变异"
		result.append("检出{0} {1}，为{2}".format(var["gene_symbol"], hgvs, clinic_g_cn))
	return "；".join(result)
jinja2.filters.FILTERS["whxh_var_sum_g"] = whxh_var_sum_g

# 2026.02.03-武汉协和MP组织-变异检测结果-氨基酸展示三字母
def whxh_var_info(var):
	result = []
	if "bio_category" in var.keys():
		if var["bio_category"] == "Snvindel":
			if var["hgvs_p"] != "p.?":
				result.append(var["gene_region"]+" "+var["hgvs_c"]+" " + var["hgvs_p"] + " " +var["hgvs_p_abbr"])
			else:
				result.append(var["gene_region"]+" "+var["hgvs_c"])
			result.append(var["transcript_primary"])
		elif var["bio_category"] == "Cnv":
			if var["var_origin"] == "germline" and var["gene_symbol"] in ["BRCA1", "BRCA2"]:
				if var["type"] == "Loss":
					result.append(var["value"] + " del")
				else:
					result.append(var["value"] + " dup")
			else:
				result.append("扩增") 
		elif var["bio_category"] in ["Sv", "PSeqRnaSv"]:
			if var["five_prime_gene"] == "MET" and var["three_prime_gene"] == "MET":
				result.append("MET exon14 跳跃")
				result.append(var["five_prime_transcript"])
			else:
				result.append("{0}:{1}-{2}:{3}融合".format(var["five_prime_gene"], var["five_prime_cds"], var["three_prime_gene"], var["three_prime_cds"]))
				result.append("{0}/{1}".format(var["five_prime_transcript"], var["three_prime_transcript"]))
		elif var["bio_category"] == "PMLPA":
			if var["type"] == "Loss":
				result.append(var["value"]+" del")
			elif var["type"] == "Gain":
				result.append(var["value"]+" dup")
		# 2025.07.17-新增PHd
		elif var["bio_category"] == "PHd":
			if var["type"] == "HomoDel":
				result.append("纯合缺失")
			elif var["type"] == "HeteDel":
				result.append("杂合缺失")
			else:
				result.append("未知变异类型！")
			result.append("（{0}）".format(var["region"]))
		# 2025.07.17-新增完成
	return "\n".join(result)
jinja2.filters.FILTERS["whxh_var_info"] = whxh_var_info

def whxh_get_varsimpleinfo(var):
	if var["bio_category"] == "Snvindel":
		return var["hgvs_p_abbr"] if var["hgvs_p_abbr"] != "p.?" else var["hgvs_c"]
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

# 2026.02.03-武汉协和MP组织-NCCN指南推荐基因检测结果，检测结果展示-三字母氨基酸
def whxh_cdx_type(info):
	cdx_list = info[0]
	sample = info[1]
	result = []
	if cdx_list:
		for var in cdx_list:
			if var["bio_category"] in ["Sv", "PSeqRnaSv"]:
				if var["five_prime_gene"] == "MET" and var["three_prime_gene"] == "MET":
					result.append("MET exon14 跳跃，{0}".format(var["freq"]))
				else:
					if "var_info_M" in var.keys() and var["var_info_M"]:
						result.append(var["var_info_M"]+"，"+var["freq"])
					else:
						result.append(var["var_info"]+"，"+var["freq"])
			else:
				var_info = whxh_get_varsimpleinfo(var)
				if var["gene_symbol"] in ["BRCA1", "BRCA2"]:
					if sample["control_sample_id"]:
						if var["var_origin"] == "germline":
							result.append("{0}，{1}".format(var_info, var["freq_rc_str"]))
							result.append("胚系变异")
						else:
							result.append("{0}，{1}".format(var_info, var["freq"]))
							result.append("体细胞变异")
					else:
						result.append("{0}，{1}".format(var_info, var["freq"]))
				else:
					result.append("{0}，{1}".format(var_info, var["freq"]))
	return "\n".join(result)
jinja2.filters.FILTERS["whxh_cdx_type"] = whxh_cdx_type

# 2026.02.03-武汉协和MP-免疫检查点抑制剂疗效相关基因-氨基酸展示三字母
def whxh_io_summary(io):
	result = []
	if io["abbr_hd_io_p_summary"]:
		result.append(io["abbr_hd_io_p_summary"]+"（疗效正相关）")
	if io["abbr_hd_io_n_summary"]:
		result.append(io["abbr_hd_io_n_summary"]+"（疗效负相关）")
	if not result:
		result = ["未检出相关变异"]
	return "；\n".join(result)+"。"
jinja2.filters.FILTERS["whxh_io_summary"] = whxh_io_summary

# 2026.02.03-武汉协和MP治疗方案介绍-生物标志物-氨基酸三字母
def whxh_approval_regimen_biomarker(info):
	biomaker_list = info[0]
	judge_mergeMET = info[1]
	result = []
	for i in biomaker_list:
		if "hgvs_c" in i.keys() and i["hgvs_c"]:
			if i["hgvs_p"] != "p.?":
				i["hgvs_p_abbr"] = i["hgvs_p_abbr"] if "hgvs_p_abbr" in i.keys() and i["hgvs_p_abbr"] else splitAA(i["hgvs_p"])
				result.append("{0} {1} {2} {3}".format(i["gene_symbol"], i["gene_region"], i["hgvs_c"], i["hgvs_p_abbr"]))
			else:
				result.append("{0} {1} {2}".format(i["gene_symbol"], i["gene_region"], i["hgvs_c"]))
		elif "cnv_type" in i.keys() and i["cnv_type"]:
			# 2024.08.30-CNV 区分Loss，其他的写扩增
			if i["cnv_type"] == "Loss" and i["gene_symbol"] in ["BRCA1", "BRCA2"]:
				result.append("{0} 大片段缺失".format(i["gene_symbol"]))
			elif i["cnv_type"] == "Gain" and i["gene_symbol"] in ["BRCA1", "BRCA2"]:
				result.append("{0} 大片段重复".format(i["gene_symbol"]))
			elif i["cnv_type"] == "HeteDel" and i["gene_symbol"] in ["BRCA1", "BRCA2"]:
				result.append("{0} 杂合大片段缺失".format(i["gene_symbol"]))
			elif i["cnv_type"] == "HomoDel" and i["gene_symbol"] in ["BRCA1", "BRCA2"]:
				result.append("{0} 纯合大片段缺失".format(i["gene_symbol"]))
			else:
				result.append("{0} 扩增".format(i["gene_symbol"]))
		elif "five_prime_gene" in i.keys() and i["five_prime_gene"]:
			if i["five_prime_gene"] == "MET" and i["three_prime_gene"] == "MET":
				result.append("MET exon14 跳跃")
			else:
				# 重新拆分hgvs，CP40的region少了ins
				if "hgvs" in i.keys() and i["hgvs"]:
					i["five_prime_region"] = "-".join(re.split("-", (re.split(":", i["hgvs"])[2]))[:-1]) \
						  					 if not re.search("--", i["hgvs"]) \
											 else re.split("_", (re.split("--", i["hgvs"])[0]))[-1]
					i["three_prime_region"] = re.split(":", i["hgvs"])[-1] \
											  if not re.search("--", i["hgvs"]) \
											  else re.split("_", (re.split("--", i["hgvs"])[1]))[-1]
					# 加一个兼容-2023.10.19
					# var_hgvs新格式，gene1:NM_xxx:exon1--gene2:NM_xxx:exon2, 旧的为gene1:NM_xxx_exon1--gene2:NM_xxx_exon2
					# cds会变成xxx:exon1和xxx:exon2
					i["five_prime_region"] = re.split(":", i["five_prime_region"])[-1] if re.search(":", i["five_prime_region"]) else i["five_prime_region"]
					i["three_prime_region"] = re.split(":", i["three_prime_region"])[-1] if re.search(":", i["three_prime_region"]) else i["three_prime_region"]
					# 兼容完成-2023.10.19
					# 加一个兼容-2024.01.25
					# 4. gene:转录本_exon-gene:转录本_exon 重新提取后，five_prime_cds为空，以此做为重新拆分的判定依据
					if not i["five_prime_region"]:
						i["five_prime_region"] = re.split("_", re.split("-", i["hgvs"])[0])[-1]
						i["three_prime_region"] = re.split("_", i["hgvs"])[-1]
					# 兼容完成-2024.01.25
				result.append("{0}:{1}:{2}-{3}:{4}:{5}".format(i["five_prime_gene"], i["five_prime_transcript"], i["five_prime_region"],\
															   i["three_prime_gene"], i["three_prime_transcript"], i["three_prime_region"]))
		elif "biomarker_type" in i.keys() and i["biomarker_type"]:
			if i["biomarker_type"] == "KRAS/NRAS/BRAF WT":
				result.append("KRAS/NRAS/BRAF 野生型")
			elif i["biomarker_type"] == "HRD-":
				result.append("HRD阴性")
			elif i["biomarker_type"] == "HRD+":
				result.append("HRD阳性")
			# 2025.07.17-新增HD
			elif "HomoDel" in i["biomarker_type"] or "HeteDel" in i["biomarker_type"]:
				split_str = re.split(":", i["biomarker_type"])
				if len(split_str) == 4:
					hd_type = "纯合缺失" if split_str[-1] == "HomoDel" else "杂合缺失" if split_str[-1] == "HeteDel" else "未知变异类型！"
					result.append(split_str[1] + " " + hd_type + " " + split_str[2])
					#result.append(" ".join(split_str[1:4]))
				else:
					result.append(i["biomarker_type"])
			else:
				result.append(i["biomarker_type"])
		else:
			result.append("无法分辨的分子标志物！")
	# 若存在MET 14跳跃DNA/RNA共检的话，则删除RNA里的结果，仅保留DNA
	result_redup = []
	for i in result:
		if i not in result_redup:
			result_redup.append(i)
	# 2024.02.27更新-MET DNA/RNA共检时都要展示
	#if judge_mergeMET:
	#	if "MET exon14 跳跃" in result_redup:
	#		result_redup.remove("MET exon14 跳跃")
	# 2024.02.27更新完成
	if not result_redup:
		result_redup = ["-"]
	#rt = RichText()
	#rt.add("\n".join(result_redup))
	return "\n".join(result_redup)
jinja2.filters.FILTERS["whxh_approval_regimen_biomarker"] = whxh_approval_regimen_biomarker

# 2026.02.04-单字母氨基酸转化为三字母
def splitAA(hgvs_p):
	transAA = {
		"A":"Ala",
		"C":"Cys",
		"D":"Asp",
		"E":"Glu",
		"F":"Phe",
		"G":"Gly",
		"H":"His",
		"I":"Ile",
		"K":"Lys",
		"L":"Leu",
		"M":"Met",
		"N":"Asn",
		"P":"Pro",
		"Q":"Gln",
		"R":"Arg",
		"S":"Ser",
		"T":"Thr",
		"V":"Val",
		"W":"Trp",
		"Y":"Tyr",
		"*":"Ter"
		}
	AA_str = []
	for i in list(hgvs_p):
		AA_str.append(transAA.get(i, i))
	return "".join(AA_str)

# 陕西人民BPTM Plus血液-所有基因检测结果汇总-2026.02.09
def sxrm_gbptm_plus_summary(var_list):
	gene_list = ["BRCA1", "BRCA2", "CTNNB1", "EPCAM", "MLH1", "MSH2", "MSH6", "PMS2", "POLE", "TP53"]
	detect_gene_list = [var["gene_symbol"] for var in var_list]
	for gene in gene_list:
		if gene not in detect_gene_list:
			var_list.append({
				"gene_symbol" : gene
			})
	return sorted(var_list, key = lambda i:gene_list.index(i["gene_symbol"]))
jinja2.filters.FILTERS["sxrm_gbptm_plus_summary"] = sxrm_gbptm_plus_summary

# 2026.02.12-复旦中山RNAseq，five_prime_gene/five_prime_cds/three_prime_gene/three_prime_cds相同时变异合并，copies相加，报告里展示一条。
def fdzs_rnaseq_merge_sv(var_list):
	sv_data = copy.deepcopy(var_list)
	for var in sv_data:
		# 获取key 赋值给hgvs_p
		var["hgvs_p"] = var["five_prime_gene"] + ":" + var["five_prime_cds"] + "-" + var["three_prime_gene"] + ":" + var["three_prime_cds"]
		print ("1111", var["hgvs_p"], var["freq"])
		Predictive_evi = [evi["regimen_name"] for evi in var["evi_sum"]["evi_split"]["Predictive"]] if "Predictive" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Predictive"] else []
		Prognostic_evi = [evi["clinical_significance"] for evi in var["evi_sum"]["evi_split"]["Prognostic"]] if "Prognostic" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Prognostic"] else []
		Diagnostic_evi = ["Diagnostic" for evi in var["evi_sum"]["evi_split"]["Diagnostic"]] if "Diagnostic" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Diagnostic"] else []
		evi_key = []
		evi_key.extend(sorted(Predictive_evi))
		evi_key.extend(sorted(Prognostic_evi))
		evi_key.extend(sorted(Diagnostic_evi))
		var["hgvs_p"] += "-".join(evi_key)
	# 加个排序，使两端基因相同的变异排到一起，groupby对于ABBA排序的变异处理有点问题
	sv_data = sorted(sv_data, key=lambda i:(i["five_prime_gene"], var["five_prime_cds"], i["three_prime_gene"], var["three_prime_cds"], i["hgvs_p"]))
	sv_sum = groupby(sv_data, key=lambda x : x["hgvs_p"])
	sv_combination = []
	for i, j in sv_sum:
		sv_group = list(j)
		major_sv = sv_group[0]
		major_sv['freq'] = sum(int(k["freq"]) for k in sv_group if k["freq"])	
		sv_combination.append(major_sv)
	return sv_combination
jinja2.filters.FILTERS["fdzs_rnaseq_merge_sv"] = fdzs_rnaseq_merge_sv

# 2026.02.12-菏泽市立ptHRR-解读部分BRCA需要过滤出体细胞的
def hzsl_pthrr_filter_somatic(var_list):
	return [var for var in var_list if var["var_origin"] == "somatic"]
jinja2.filters.FILTERS["hzsl_pthrr_filter_somatic"] = hzsl_pthrr_filter_somatic

# 2026.02.25-北三MP过滤出VHL基因变异(包含)（snvindel）
def bds_filter_vhl(var_list):
	return [var for var in var_list if var["gene_symbol"] == "VHL" and var["bio_category"] == "Snvindel"]
jinja2.filters.FILTERS["bds_filter_vhl"] = bds_filter_vhl

# 2026.02.25-北三MP过滤出VHL基因变异（不包含）（snvindel）
def bds_filter_novhl(var_list):
	return [var for var in var_list if not (var["gene_symbol"] == "VHL" and var["bio_category"] == "Snvindel")]
jinja2.filters.FILTERS["bds_filter_novhl"] = bds_filter_novhl

# 2026.02.25-北三MP胚系有药的过滤出4/5类变异
def bds_mp_filter_germline45(var_list):
	return [var for var in var_list if var["clinic_num_g"] in [4,5]]
jinja2.filters.FILTERS["bds_mp_filter_germline45"] = bds_mp_filter_germline45

# 2026.02.26-输入两个版本号A,B，比较大小，A<B的返回false，A>=B的返回true
def judge_version(info):
	json_batch_name = info[0]
	versionA = get_analysis_version(json_batch_name)
	versionB = info[1]
	# 处理为0.1.2的格式，并且以.为隔断进行拆分
	versionA = versionA.replace("V", "").replace("v", "").split("_")[1].split(".") if versionA else ""
	versionB = versionB.replace("V", "").replace("v", "").split(".") if versionB else ""
	# 处理后的格式为[0, 7, 0]，三位版本号。未返回版本号的就按最新的展示
	if versionA and versionB:
		for i in range(len(versionA)):
			a = int(versionA[i])
			b = int(versionB[i]) if i < len(versionB) else 0
			if a == b:
				continue
			else:
				if a - b > 0:
					return True
				else:
					return False
		return True
	else:
		return True	
jinja2.filters.FILTERS["judge_version"] = judge_version

# 2026.03.04-温附一HRD-综合评估表格的临床意义，需要展示HRD和BRCA的证据
# BRCA和HRD均有的治疗方案，以BRCA的描述为准
def wfy_hrd_regimen_sum(info):
	gss_evi = info[0]
	brca_list = info[1]
	sum_evi = []
	brca_regimen = []
	for var in brca_list:
		if "Predictive" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Predictive"]:
			for evi in var["evi_sum"]["evi_split"]["Predictive"]:
				if evi not in sum_evi:
					sum_evi.append(evi)
					brca_regimen.append(evi["regimen_name"])
	if gss_evi and "evi_split" in gss_evi.keys() and gss_evi["evi_split"] and "Predictive" in gss_evi["evi_split"].keys() and gss_evi["evi_split"]["Predictive"]:
		for evi in gss_evi["evi_split"]["Predictive"]:
			if evi["regimen_name"] not in brca_regimen:
				sum_evi.append(evi)
	
	# 排序、合并
	sum_evi = sorted(sum_evi, key = lambda i:(i["evi_conclusion_simple"], i["sense_rule"], i["regimen_name_py"].upper()))

	tmp_dict = {}
	for evi in sum_evi:
		interpretation = evi["evi_adaptation_disease_cn"] if evi["evi_adaptation_disease_cn"] else evi["evi_interpretation"]
		if interpretation not in tmp_dict.keys():
			tmp_dict.setdefault(interpretation, [])
		tmp_dict[interpretation].append(evi["regimen_name"])
	
	merge_evi_sum = []
	for k, v in tmp_dict.items():
		merge_evi_sum.append(
			{
				"regimen_name" : "、".join(v),
				"interpretation" : k
			}
		)

	if info[2] == "regimen":
		return sum_evi
	else:
		return merge_evi_sum
jinja2.filters.FILTERS["wfy_hrd_regimen_sum"] = wfy_hrd_regimen_sum

# 2026-03-06-中六150统计各证据项数量
def zsly_transform_evi_codes(info):
	evidence_categorys = info[0]
	# 对带下划线的证据项进行转化，首字母+后缀转化
	suffix_dict = {
		"Moderate" : "M",
		"Supporting" : "P",
		"Strong" : "S",
		"Alone" : "A",
		"VeryStrong" : "VS"
	}
	code_list = re.split(";", evidence_categorys) if evidence_categorys else []
	code_stran = []
	for code in code_list:
		if "_" in code:
			code_stran.append(code.split("_")[0][0] + suffix_dict.get(code.split("_")[1], code.split("_")[1]))
		else:
			code_stran.append(code)

	prefixes = []
	for code in code_stran:
		match = re.match(r"^([A-Z]+)", code)
		if match:
			prefixes.append(match.group(1))
	counts = Counter(prefixes)

	# 排序
	# 4/5类变异按照PVS > PS > PM > PP > BA > BS > BM >BP；3类按B > P
	priority_P = {
		"PVS": 1, "PS": 2, "PM": 3, "PP": 4, 
		"BA": 5, "BS": 6, "BM": 7, "BP": 8
	}
	priority_N = {
		"PVS": 5, "PS": 6, "PM": 7, "PP": 8, 
		"BA": 1, "BS": 2, "BM": 3, "BP": 4
	}
	if info[1] == "P":
		sorted_items = sorted(counts.items(), key=lambda x: priority_P.get(x[0], 99))
	else:
		sorted_items = sorted(counts.items(), key=lambda x: priority_N.get(x[0], 99))
	
	return "+".join([f"{count}{prefix}" for prefix, count in sorted_items])
jinja2.filters.FILTERS["zsly_transform_evi_codes"] = zsly_transform_evi_codes

# 2026.03.06-新增MLH1基因，适用最新分析流程
def hrr_parp_v2(var_list):
	gene_list = ["ATM","BARD1","BRCA1","BRCA2","BRIP1","CDK12","CHEK1","CHEK2","FANCL","PALB2",\
			  	 "RAD51B","RAD51C","RAD51D","RAD54L", "FANCA", "ATR", "MRE11", "NBN", "MLH1"]
	detect_gene = [var["gene_symbol"] for var in var_list if var["gene_symbol"] in gene_list]
	result = []
	# 2025.12.23:修复bug，var_list需要根据基因列表筛一下
	for var in var_list:
		if var["gene_symbol"] in gene_list:
			result.append(var)
	#result.extend(var_list)			
	for gene in set(gene_list) - set(detect_gene):
		result.append({
			"gene_symbol" : gene
		})	
	return sorted(result, key=lambda i:i["gene_symbol"])
jinja2.filters.FILTERS["hrr_parp_v2"] = hrr_parp_v2

# 2026.03.11-WES通用模板，III类变异只取前200个变异
# 按频率降序重新排一下
# >200个变异，取snvindel前200个，不足200个则展示所有snvindel
# <= 200个变异，展示所有snvindel和cnv
def wes_get200var(var_list):
	snv_var_list = [var for var in var_list if var["bio_category"] == "Snvindel"]
	cnv_var_list = [var for var in var_list if var["bio_category"] == "Cnv"]
	snv_var_list = sorted(snv_var_list, key = lambda i : float(i["freq"]), reverse=True)
	cnv_var_list = sorted(cnv_var_list, key = lambda i : float(i["cn_mean"]), reverse=True)

	if len(var_list) > 200:
		return snv_var_list[:200] if len(snv_var_list) > 200 else snv_var_list
	else:
		return snv_var_list + cnv_var_list
jinja2.filters.FILTERS["wes_get200var"] = wes_get200var

# 同济MP-I类变异不展示D级证据-2026.03.18
def tj_mp_filter_D(regimen_list_raw):
	if regimen_list_raw != "-":
		regimen_list = re.split("、", regimen_list_raw)
		result = []
		for regimen in regimen_list:
			if re.search("A级|B级|C级",regimen):
				result.append(regimen)
		return "、".join(result)
	else:
		return "-"
jinja2.filters.FILTERS["tj_mp_filter_D"] = tj_mp_filter_D

# 西京医院（中国人民解放军空军军医大学第一附属医院）MP-I类变异汇总A级敏感药物
def xj_mp_a_regimen_sense(regimen_list):
	a_sense_list = [i["regimen_name"] for i in regimen_list if i["clinical_significance_cn"] == "敏感" and i["evi_conclusion_simple"] == "A"]
	if a_sense_list:
		return "针对携带该变异的患者，获批或指南推荐的药物有" + "、".join(a_sense_list) + "等。"
	else:
		return ""
jinja2.filters.FILTERS["xj_mp_a_regimen_sense"] = xj_mp_a_regimen_sense

# 2026.03.20-同济MP输出EML4-ALK融合具体的亚型
def tj_mp_alksv(var):
	alk_sv_dict = {
		("exon13", "exon20") : "V1",
		("exon20", "exon20") : "V2",
		("exon6", "exon20") : "V3a/b",
		("exon15", "exon20") : "V4'",
		("exon2", "exon20") : "V5a/b",
		("exon18", "exon20") : "V5'",
		("exon14", "exon20") : "V7",
		("exon17", "exon20") : "V8a/b"
	}
	sv_type = ""
	if var["five_prime_gene"] == "EML4" and var["three_prime_gene"] == "ALK":
		sv_type = alk_sv_dict.get((var["five_prime_cds"], var["three_prime_cds"]), "")
	return sv_type
jinja2.filters.FILTERS["tj_mp_alksv"] = tj_mp_alksv

# 重庆附一EGFR PACC位点添加标注-2026.04.02
# 仅考虑snvindel变异
def cqfy_116_judge_pacc(var):
	hgvs_p_list = ["p.A647T", "p.R675W", "p.E709A", "p.E709K", "p.E709V", "p.E709_T710delinsD", "p.L718Q", "p.L718V", \
				   "p.G719A", "p.G719C", "p.G719S", "p.G724S", "p.E736K", "p.I740_K745dup", "p.L747P", "p.L747S", \
				   "p.A750_I759delinsPN", "p.T751_I759delinsN", "p.S752_I759del", "p.A755D", "p.K757M", "p.K757R", "p.V765L", "p.S768C", \
				   "p.S768I", "p.V769L", "p.V769M", "p.N771G", "p.V774M", "p.R776H", "p.R776C", "p.G779F", \
				   "p.L792H", "p.G796S", "p.C797S", "p.T854I", "p.T854S"]
	if var["bio_category"] == "Snvindel" and var["gene_symbol"] == "EGFR" and var["hgvs_p"].replace("(", "").replace(")", "") in hgvs_p_list:
		return True
	else:
		return False
jinja2.filters.FILTERS["cqfy_116_judge_pacc"] = cqfy_116_judge_pacc

# 吉大一MP获取变异列表中的CNV+HD-2026.04.07
def filter_cnv_hd(var_list):
	return [var for var in var_list if var["bio_category"] in ["Cnv", "PHd"]]
jinja2.filters.FILTERS["filter_cnv_hd"] = filter_cnv_hd

# 获取变异列表中的DNA SV-2026.04.07
def filter_dna_sv(var_list):
	return [var for var in var_list if var["bio_category"] in ["Sv"]]
jinja2.filters.FILTERS["filter_dna_sv"] = filter_dna_sv

# 适用于山东齐鲁gHRR-v2，展示胚系5/4/3类变异-2026.04.13
# 增加基因AKT1和MLH1，其中MLH1是获批基因
def sdql_ghrr_summary_all_v2(var_list):
	gene_list_all = ["BRCA1", "BRCA2", "AR", "ATM", "ATR", "BARD1", "BRIP1", "CDH1", "CDK12", "CHEK1", "CHEK2", "ESR1", \
					 "FANCA", "FANCL", "HDAC2", "HOXB13", "MRE11", "NBN", "PALB2", "PPP2R2A", "PTEN", "RAD51B", "RAD51C", \
					 "RAD51D", "RAD54L", "STK11", "TP53", "BRAF", "ERBB2", "KRAS", "NRAS", "PIK3CA", "AKT1", "MLH1"]
	gene_list_appr = ["BRCA1", "BRCA2", "ATM", "BARD1", "BRIP1", "CDK12", "CHEK1", "CHEK2", "FANCL", "PALB2", "RAD51B", "RAD51C", "RAD51D", "RAD54L", "MLH1"]
	detect_gene = [var["gene_symbol"] for var in var_list]
	for gene in set(gene_list_all) - set(detect_gene):
		var_list.append({"gene_symbol" : gene})
	for i in var_list:
		if i["gene_symbol"] in  gene_list_appr:
			i["appr"] = "T"
	return sorted(var_list, key = lambda i:i["gene_symbol"])
jinja2.filters.FILTERS["sdql_ghrr_summary_all_v2"] = sdql_ghrr_summary_all_v2

# 2026.04.20-南充中心116，变异区分为获证基因（EGFR 20ins）和其他基因
def nczx_judge_egfr(info):
	var_list = info[0]
	egfr_var = []
	other_var = []
	for var in var_list:
		if var["bio_category"] == "Snvindel" and var["gene_symbol"] == "EGFR" and var["gene_region"] == "exon20" and "EGFR Exon20 ins" in var["var_category_names"]:
			egfr_var.append(var)
		else:
			other_var.append(var)

	if info[1] == "inegfr":
		return egfr_var
	else:
		return other_var
jinja2.filters.FILTERS["nczx_judge_egfr"] = nczx_judge_egfr	

# 南充中心HRR，药物介绍不展示BRCA基因相关的-2026.04.20
def nczx_drug_nobrca(drug_list):
	for drug in drug_list:
		drug["var_filter"] = [a for a in drug["var"] if ("gene_symbol" in a.keys() and a["gene_symbol"] and a["gene_symbol"] not in ["BRCA1", "BRCA2"]) and not \
							  ("biomarker_type" in a.keys() and "BRCA1" in a["biomarker_type"] or \
							  "biomarker_type" in a.keys() and "BRCA2" in a["biomarker_type"])]
	return [a for a in drug_list if a["var_filter"]]
jinja2.filters.FILTERS["nczx_drug_nobrca"] = nczx_drug_nobrca

# 2026.04.21-特医MP-拷贝数变异按照变异等级排序
def ty_mp_cnv_hd_sort(var_list):
	for var in var_list:
		if var["evi_sum"]["evi_split"] and set(["Diagnostic","Predictive","Prognostic"]) & set(var["evi_sum"]["evi_split"].keys()):
			var["ty_sort_num"] = 0 if var["clinic_num_s"] == 5 else 1
		else:
			var["ty_sort_num"] = 2 if var["clinic_num_s"] == 5 else 3 if var["clinic_num_s"] == 4 else 4
	return sorted(var_list, key = lambda i:i["ty_sort_num"])
jinja2.filters.FILTERS["ty_mp_cnv_hd_sort"] = ty_mp_cnv_hd_sort

# 华东医院HRR-34基因，未检测到变异的基因列表需要展示-2026.04.21
def filter_hrr_fdhd_noresult_gene_v2(var_list):
	hrr_gene = ["AR", "ATM", "ATR", "BARD1", "BRAF", "BRCA1", "BRCA2", "BRIP1", "CDH1", "CDK12", \
				"CHEK1", "CHEK2", "ERBB2", "ESR1", "FANCA", "FANCL", "HDAC2", "HOXB13", "KRAS", \
				"MRE11", "NBN", "NRAS", "PALB2", "PIK3CA", "PPP2R2A", "PTEN", "RAD51B", "RAD51C", \
				"RAD51D", "RAD54L", "STK11", "TP53", "AKT1", "MLH1"]
	var_gene = [var["gene_symbol"] for var in var_list]
	result = []
	for gene in hrr_gene:
		if gene not in var_gene:
			result.append(gene)
	return ", ".join(sorted(result))
jinja2.filters.FILTERS["filter_hrr_fdhd_noresult_gene_v2"] = filter_hrr_fdhd_noresult_gene_v2

# 适用于广东人民gHRR-34基因，展示胚系5/4/3类变异-2026.04.21
def gdrm_ghrr_summary_all_v2(var_list):
	gene_list = ["BRCA1","BRCA2","AKT1", "AR","ATM","ATR","BARD1","BRIP1","CDH1","CDK12","CHEK1",\
				 "CHEK2","ESR1","FANCA","FANCL","HDAC2","HOXB13", "MLH1", "MRE11","NBN","PALB2","PPP2R2A",\
				 "PTEN","RAD51B","RAD51C","RAD51D","RAD54L","STK11","TP53","BRAF","ERBB2","KRAS",\
				 "NRAS","PIK3CA"]
	detect_gene = [var["gene_symbol"] for var in var_list if var["gene_symbol"] in gene_list]
	for gene in set(gene_list) - set(detect_gene):
		var_list.append({
			"gene_symbol" : gene
		})	
	return sorted(var_list, key=lambda i:gene_list.index(i["gene_symbol"]))
jinja2.filters.FILTERS["gdrm_ghrr_summary_all_v2"] = gdrm_ghrr_summary_all_v2

# 甘肃武威HRR-34基因-BRCA和其他基因检测结果-2026.04.21
def gsww_hrr_var_v2(info):
	var_raw_list = info[0]
	gene_type = info[1]
	# brca基因结果
	result_brca = [var for var in var_raw_list if var["gene_symbol"] in ["BRCA1", "BRCA2"]]
	for gene in set(["BRCA1", "BRCA2"]) - set([var["gene_symbol"] for var in result_brca]):
		result_brca.append({"gene_symbol" : gene})
	result_brca = sorted(result_brca, key = lambda i:i["gene_symbol"])
	# 其他基因结果
	other_gene_list = ["AR", "ATM", "ATR", "BARD1", "BRAF", "BRIP1", "CDH1", "CDK12", "CHEK1", "CHEK2", \
					   "ERBB2", "ESR1", "FANCA", "FANCL", "HDAC2", "HOXB13", "KRAS", "MRE11", "NBN", "NRAS", \
					   "PALB2", "PIK3CA", "PPP2R2A", "PTEN", "RAD51B", "RAD51C", "RAD51D", "RAD54L", "STK11", "TP53", "AKT1", "MLH1"]
	result_other = [var for var in var_raw_list if var["gene_symbol"] not in ["BRCA1", "BRCA2"]]
	for gene in set(other_gene_list) - set([var["gene_symbol"] for var in result_other]):
		result_other.append({"gene_symbol" : gene})
	result_other = sorted(result_other, key = lambda i:i["gene_symbol"])
	# 根据需求返回对应的结果
	if gene_type == "brca":
		return result_brca
	elif gene_type == "other":
		return result_other
jinja2.filters.FILTERS["gsww_hrr_var_v2"] = gsww_hrr_var_v2

# 2026.04.21-南通附属CP200小结，需要增加PACC标签
def ntfs_cp200_var_sum_s(var_list):
	result = []
	for var in var_list:
		if var["bio_category"] == "Snvindel":
			if var["hgvs_p"] != "p.?":
				if cqfy_116_judge_pacc(var):
					result.append(var["gene_symbol"]+" "+var["hgvs_p"]+"，该突变属于EGFR基因PACC变异")
				else:
					result.append(var["gene_symbol"]+" "+var["hgvs_p"])
			else:
				if cqfy_116_judge_pacc(var):
					result.append(var["gene_symbol"]+" "+var["hgvs_c"]+"，该突变属于EGFR基因PACC变异")
				else:
					result.append(var["gene_symbol"]+" "+var["hgvs_c"])
		elif var["bio_category"] == "Cnv":
			# 2026.01.12-若cnv_type为Loss/loss，返回缺失
			if "cnv_type" in var.keys() and var["cnv_type"] and var["cnv_type"] in ["Loss", "loss"]:
				result.append(var["gene_symbol"]+" 缺失")
			else:
				result.append(var["gene_symbol"]+" 扩增")
		elif var["bio_category"] in ["Sv", "PSeqRnaSv"]:
			if var["five_prime_gene"] == "MET" and var["three_prime_gene"] == "MET":
				result.append("MET exon14 跳跃")
			else:
				# 融合可能会有重复（rna exon相同，断点不同的情况）
				if var["five_prime_gene"]+"-"+var["three_prime_gene"]+" 融合" not in result:
					result.append(var["five_prime_gene"]+"-"+var["three_prime_gene"]+" 融合")
		# 2025.07.15-新增PHd
		elif var["bio_category"] == "PHd":
			if var["type"] == "HomoDel":
				if var["gene_symbol"]+" 纯合缺失" not in result:
					result.append(var["gene_symbol"]+" 纯合缺失")
			elif var["type"] == "HeteDel":
				if var["gene_symbol"]+" 杂合缺失" not in result:
					result.append(var["gene_symbol"]+" 杂合缺失")
			else:
				if var["gene_symbol"]+" 未知变异类型！" not in result:
					result.append(var["gene_symbol"]+" 未知变异类型！")
		# 2025.07.15-新增完成
	return "；".join(result)
jinja2.filters.FILTERS["ntfs_cp200_var_sum_s"] = ntfs_cp200_var_sum_s

# 武汉协和HRR结果汇总-34基因，获批新增MLH1-2026.04.22
def whxh_hrr_summary_v2(raw_var_list):
	gene_list = ["ATM","BARD1","BRCA1","BRCA2","BRIP1","CDK12","CHEK1","CHEK2","FANCL",\
			  	 "PALB2","RAD51B","RAD51C","RAD51D","RAD54L", "FANCA", "ATR", "MRE11", "NBN", "MLH1"]
	var_list = [var for var in raw_var_list if var["gene_symbol"] in gene_list]
	detect_gene = [var["gene_symbol"] for var in var_list if var["gene_symbol"] in gene_list]
	result = []
	result.extend(var_list)	
	for gene in set(gene_list) - set(detect_gene):
		result.append({
			"gene_symbol" : gene
		})	
	return sorted(result, key=lambda i:i["gene_symbol"])
jinja2.filters.FILTERS["whxh_hrr_summary_v2"] = whxh_hrr_summary_v2

# 2026.04.22-北京朝阳tLC10-检出指定变异才展示
def bjcy_tlc10_var(var_list):
	var_dict = {
		"KRAS" : ["G12D", "G12A", "G12V", "G12S", "G12C", "Q61H"],
		"NRAS" : ["G12D", "Q61R", "Q61K"],
		"PIK3CA" : ["H1047R"],
		"BRAF" : ["V600E"]
	}
	return [var for var in var_list if var["bio_category"] == "Snvindel" and \
		                               var["gene_symbol"] in var_dict.keys() and \
									   var["hgvs_p"].replace("p.", "").replace("(", "").replace(")", "") in var_dict[var["gene_symbol"]]]
jinja2.filters.FILTERS["bjcy_tlc10_var"] = bjcy_tlc10_var

# 2026.04.27-武汉协和MP组织-变异检测结果-氨基酸展示三字母-HD和CNV增加chr+location
def whxh_var_info_v2(var):
	result = []
	if "bio_category" in var.keys():
		if var["bio_category"] == "Snvindel":
			if var["hgvs_p"] != "p.?":
				result.append(var["gene_region"]+" "+var["hgvs_c"]+" " + var["hgvs_p"] + " " +var["hgvs_p_abbr"])
			else:
				result.append(var["gene_region"]+" "+var["hgvs_c"])
			result.append(var["transcript_primary"])
		elif var["bio_category"] == "Cnv":
			if var["var_origin"] == "germline" and var["gene_symbol"] in ["BRCA1", "BRCA2"]:
				if var["type"] == "Loss":
					result.append(var["value"] + " del")
				else:
					result.append(var["value"] + " dup")
			else:
				result.append("扩增") 
				result.append(var["chr"].replace("chr", "") + var["location"])
		elif var["bio_category"] in ["Sv", "PSeqRnaSv"]:
			if var["five_prime_gene"] == "MET" and var["three_prime_gene"] == "MET":
				result.append("MET exon14 跳跃")
				result.append(var["five_prime_transcript"])
			else:
				result.append("{0}:{1}-{2}:{3}融合".format(var["five_prime_gene"], var["five_prime_cds"], var["three_prime_gene"], var["three_prime_cds"]))
				result.append("{0}/{1}".format(var["five_prime_transcript"], var["three_prime_transcript"]))
		elif var["bio_category"] == "PMLPA":
			if var["type"] == "Loss":
				result.append(var["value"]+" del")
			elif var["type"] == "Gain":
				result.append(var["value"]+" dup")
		# 2025.07.17-新增PHd
		elif var["bio_category"] == "PHd":
			if var["type"] == "HomoDel":
				result.append("纯合缺失")
			elif var["type"] == "HeteDel":
				result.append("杂合缺失")
			else:
				result.append("未知变异类型！")
			result.append("（{0}）".format(var["region"]))
			hd_dict = {
				"ATM" : "11q22.3",
				"BARD1" : "2q35",
				"BRCA1" : "17q21.31",
				"BRCA2" : "13q13.1",
				"BRIP1" : "17q23.2",
				"CDH1" : "16q22.1",
				"CDK12" : "17q12",
				"CHEK1" : "11q24.2",
				"CHEK2" : "22q12.1",
				"FANCA" : "16q24.3",
				"FANCL" : "2p16.1",
				"HDAC2" : "6q21",
				"PALB2" : "16p12.2",
				"PPP2R2A" : "8p21.2",
				"PTEN" : "10q23.31",
				"RAD51B" : "14q24.1",
				"RAD51C" : "17q22",
				"RAD51D" : "17q12",
				"RAD54L" : "1p34.1",
				"TP53" : "17p13.1",
				"SMARCA4" : "19p13.2",
				"KEAP1" : "19p13.2",
				"STK11" : "19p13.3",
				"RB1" : "13q14.2",
				"FH" : "1q43",
				"SETD2" : "3p21.31",
				"MTAP" : "9p21.3",
				"CDKN2A" : "9p21.3",
				"CDKN2B" : "9p21.3"
			}
			result.append(hd_dict.get(var["gene_symbol"], "配置表未找到对应染色体位置，请手动补充！"))
		# 2025.07.17-新增完成
	return "\n".join(result)
jinja2.filters.FILTERS["whxh_var_info_v2"] = whxh_var_info_v2

# 广东人民tHRR结果汇总-34基因-2026.04.27
# 新增AKT1和MLH1
def gdrm_thrr_total_gene_sum_v2(var_list):
	gene_list = ["BRCA1","BRCA2","AKT1","AR","ATM","ATR","BARD1","BRIP1","CDH1","CDK12","CHEK1",\
				 "CHEK2","ESR1","FANCA","FANCL","HDAC2","HOXB13","MLH1","MRE11","NBN","PALB2","PPP2R2A",\
				 "PTEN","RAD51B","RAD51C","RAD51D","RAD54L","STK11","TP53","BRAF","ERBB2","KRAS",\
				 "NRAS","PIK3CA"]
	detect_gene = [var["gene_symbol"] for var in var_list if var["gene_symbol"] in gene_list]
	for gene in set(gene_list) - set(detect_gene):
		var_list.append({
			"gene_symbol" : gene
		})
	return sorted(var_list, key=lambda i:gene_list.index(i["gene_symbol"]))
jinja2.filters.FILTERS["gdrm_thrr_total_gene_sum_v2"] = gdrm_thrr_total_gene_sum_v2

# 武汉协和返回变异丰度，适用BRCA等，胚系展示纯合/杂合，体细胞展示频率，MLPA用-
def freq_stran_whxh(info):
	'''
	体细胞snvindel：freq_str
	胚系snvindel：厦门项目freq， 上海项目freq_rc，注意MP还要考虑是否有对照样本，无对照样本，来源为germline的直接展示freq_str
	sv：CP40 copies，其他freq_str
	PSeqRnaSv：freq
	'''
	var = info[0]
	sample = info[1]
	if "var_origin" in var.keys():
		if var["bio_category"] == "Snvindel":
			if var["var_origin"] == "germline":
				if "Master" in sample["prod_names"]:
					if sample["control_sample_id"]:
						return "纯合" if var["freq_rc"] and float(var["freq_rc"]) >= 0.85 else "杂合"
					else:
						return var["freq_str"]
				elif re.search("116|76|25|21|18", sample["prod_names"]):
					return "纯合" if var["freq_rc"] and float(var["freq_rc"]) >= 0.85 else "杂合" if var["freq_rc"] and float(var["freq_rc"]) < 0.85 else "未提取到freq_rc！"
				else:
					return "纯合" if float(var["freq"]) >= 0.85 else "杂合"
			else:
				return var["freq_str"]
		elif var["bio_category"] == "Cnv":
			return "CN：" + var["cn_mean"].replace(".00", ".0")
		elif var["bio_category"] == "Sv":
			if re.search("Classic|CRC12", sample["prod_names"]):
				return str(var["copies"])+" copies"
			else:
				return var["freq_str"]
		elif var["bio_category"] == "PSeqRnaSv":
			return str(var["freq"])+" copies"
		# 2025.07.01-新增PHd
		elif var["bio_category"] == "PHd":
			return "CN：0.0"
		# 2025.07.01-新增完成
	else:
		if var["type"] in ["Loss", "Gain"]:
			return "/"
jinja2.filters.FILTERS["freq_stran_whxh"] = freq_stran_whxh

# 2026.04.30-南充中心116，变异区分为获证基因（EGFR）和其他基因
def nczx_judge_egfr_v2(info):
	var_list = info[0]
	egfr_var = []
	other_var = []
	for var in var_list:
		if var["bio_category"] == "Snvindel" and var["gene_symbol"] == "EGFR":
			egfr_var.append(var)
		else:
			other_var.append(var)

	if info[1] == "inegfr":
		return egfr_var
	else:
		return other_var
jinja2.filters.FILTERS["nczx_judge_egfr_v2"] = nczx_judge_egfr_v2	

# 2026.04.30-中六CP200小结需要添加角标，表格底部根据角标出现顺序添加备注
# 小结展示顺序：I类、II类、MSI、免疫正负相关、化疗
# snvindel freq >= 0.01 且freq <= 0.05 ==> Snvindel
# snvindel freq < 0.01 ==> Snvindel_2
# cnv ERBB2/MET cn_mean >= 3.5 且 cn_mean <= 5， cnv其他cn_mean >= 6且cn_mean <= 8添加上角标 ==> Cnv
# MSI msi_score >=0.12 且msi_score <= 0.24添加上角标 ==>MSI
# 化疗 ==> loss1/loss2/more
def zsly_cp200_summary_add_subscript(info):
	var_list = info[0]
	msi = info[1]
	io_list = info[2]
	chemo_list = info[3]
	all_detect_result = info[4]
	result = []
	for var in var_list:
		if "bio_category" in var.keys() and var["bio_category"] and var["bio_category"] == "Snvindel":
			if float(var["freq"]) >= 0.01 and float(var["freq"]) <= 0.05:
				if "Snvindel" not in result:
					result.append("Snvindel")
			elif float(var["freq"]) < 0.01:
				if "Snvindel_2" not in result: #修改Snvidel为Snvindel 2026.05.07 孟智悦
					result.append("Snvindel_2")
		elif "bio_category" in var.keys() and var["bio_category"] and var["bio_category"] == "Cnv":
			if var["gene_symbol"] in ["ERBB2", "MET"]:
				if float(var["cn_mean"]) >= 3.5 and float(var["cn_mean"]) <= 5:
					if "Cnv" not in result:
						result.append("Cnv")
			else:
				if float(var["cn_mean"]) >= 6 and float(var["cn_mean"]) <= 8:
					if "Cnv" not in result:
						result.append("Cnv")
	if float(msi["msi_score"]) >= 0.12 and float(msi["msi_score"]) <= 0.24:
		result.append("MSI")
	for i in io_list:
		for var in i["var_info"]:
			if "bio_category" in var.keys() and var["bio_category"] and var["bio_category"] == "Snvindel":
				if float(var["freq"]) >= 0.01 and float(var["freq"]) <= 0.05:
					if "Snvindel" not in result:
						result.append("Snvindel")
				elif float(var["freq"]) < 0.01:
					if "Snvindel_2" not in result: #修改Snvidel为Snvindel 2026.05.07 孟智悦
						result.append("Snvindel_2")
			elif "bio_category" in var.keys() and var["bio_category"] and var["bio_category"] == "Cnv":
				if var["gene_symbol"] in ["ERBB2", "MET"]:
					if float(var["cn_mean"]) >= 3.5 and float(var["cn_mean"]) <= 5:
						if "Cnv" not in result:
							result.append("Cnv")
				else:
					if float(var["cn_mean"]) >= 6 and float(var["cn_mean"]) <= 8:
						if "Cnv" not in result:
							result.append("Cnv")
	for i in chemo_list:
		for var in i["result"]:
			if var["info"] != "野生型":
				if zsly_chemo_var_v3([var, all_detect_result]) and zsly_chemo_var_v3([var, all_detect_result]) in ["loss1"]:
					if "loss1" not in result:
						result.append("loss1")
				if zsly_chemo_var_v3([var, all_detect_result]) and zsly_chemo_var_v3([var, all_detect_result]) in ["loss2"]:
					if "loss2" not in result:
						result.append("loss2")
				if zsly_chemo_var_v3([var, all_detect_result]) and zsly_chemo_var_v3([var, all_detect_result]) in ["more"]:
					if "more" not in result:
						result.append("more")
	return result
jinja2.filters.FILTERS["zsly_cp200_summary_add_subscript"] = zsly_cp200_summary_add_subscript

# 中六CP200对化疗结果进行变形，方便添加上标，化疗需要汇总野生型/变异型等 - 2026.04.30
# [{gene_symbol : AA, result : BB}, {gene_symbol : CC, result : DD}]
# [{gene_symbol : AA, result : [{info : BB1, dbsnp : jj1}, {info : BB2, dbsnp : jj2}]}]
def zsly_chemo_stran_v2(chemo_list):
	var_or_wt = {
		"rs3918290" : "C",
		"rs55886062" : "A",
		"rs67376798" : "T",
		"rs75017182" : "G",
		"rs11615" : "A",
		"rs1695" : "A",
		"rs10929302" : "G",
		"rs4148323" : "G",
		"rs8175347" : "(TA)6"
	}
	result = {}
	for chemo in chemo_list:
		if chemo["gene_symbol"] not in result.keys():
			result.setdefault(chemo["gene_symbol"], [])
		genotype_str = ""
		if chemo["dbsnp"] in var_or_wt.keys():
			ref = var_or_wt[chemo["dbsnp"]]
			split_alt = re.split("/", chemo["genotype"])
			if split_alt[0] == ref and split_alt[1] == ref:
				genotype_str == "野生型"
			else:
				if len(set(split_alt)) == 1:
					genotype_str = "纯合变异"
				else:
					genotype_str = "杂合变异"
				result[chemo["gene_symbol"]].append({
					"info" : "{0} {1}（{2}）".format(chemo["dbsnp"], chemo["genotype"], genotype_str),
					"dbsnp" : chemo["dbsnp"]
				})
	result_list = []
	for k, v in result.items():
		if not v:
			v = [{"info" : "野生型"}]
		result_list.append({
			"gene_symbol" : k,
			"result" : v
		})
	for i in result_list:
		if i["gene_symbol"] == "UGT1A1":
			i["gene_index_num"] = 0
		else:
			i["gene_index_num"] = 1
	# 对v进行格式调整，>=1的，除了最后一个，其他info末尾添加“；”
	for i in result_list:
		if len(i["result"]) > 1:
			for j in i["result"][0:-1]:
				j["info"] += "；"
			#for j in range(len(i["result"]) - 1):
			#	i["result"][j]["info"] += "；"
	#print (result_list)
	return sorted(result_list, key = lambda i : (i["gene_index_num"], i["gene_symbol"]))
jinja2.filters.FILTERS["zsly_chemo_stran_v2"] = zsly_chemo_stran_v2

# 复旦中山CP200和CP40的cnv，met/erbb2大于等于6，显示扩增，当3.5≤CN＜6，显示“拷贝数增加\n（结合IHC及FISH结果综合考虑)”。其他基因的CNV规则更新为：当CN≥15，显示“扩增”；当CN＜15，显示“拷贝数增加”
def fdzs_cnv_stran(var):
	if var['gene_symbol'] in ['MET', 'ERBB2', 'HER2']:
		if 'cn_mean' in var.keys() and var['cn_mean']:
			if float(var['cn_mean']) >= 6:
				return "扩增"
			elif float(var['cn_mean']) >= 3.5 and float(var['cn_mean']) < 6:
				return "拷贝数增加\n（结合IHC及FISH结果综合考虑)"
	else:
		if 'cn_mean' in var.keys() and var['cn_mean']:
			if float(var['cn_mean']) >= 15:
				return "扩增"
			else:
				return "拷贝数增加"
jinja2.filters.FILTERS["fdzs_cnv_stran"] = fdzs_cnv_stran

# 复旦中山MP的CNV，针对MET/ERBB2 5-9报拷贝数增加，≥10报扩增；其他基因5-14报拷贝数增加，≥15报扩增
def fdzs_mp_cnv_stran(var):
	if var['gene_symbol'] in ['MET', 'ERBB2', 'HER2']:
		if 'cn_mean' in var.keys() and var['cn_mean']:
			if float(var['cn_mean']) >= 10:
				return "扩增"
			elif float(var['cn_mean']) >= 5 and float(var['cn_mean']) < 10:
				return "拷贝数增加"
	else:
		if 'cn_mean' in var.keys() and var['cn_mean']:
			if float(var['cn_mean']) >= 15:
				return "扩增"
			elif float(var['cn_mean']) >= 5 and float(var['cn_mean']) < 15:
				return "拷贝数增加"
jinja2.filters.FILTERS["fdzs_mp_cnv_stran"] = fdzs_mp_cnv_stran

# 复旦中山MP检测小结处也要修改CNV描述
def var_sum_s_fdzs(var_list):
	result = []
	for var in var_list:
		if var["bio_category"] == "Snvindel":
			if var["hgvs_p"] != "p.?":
				result.append(var["gene_symbol"]+" "+var["hgvs_p"])
			else:
				result.append(var["gene_symbol"]+" "+var["hgvs_c"])
		elif var["bio_category"] == "Cnv":
			# 2026.01.12-若cnv_type为Loss/loss，返回缺失
			if "cnv_type" in var.keys() and var["cnv_type"] and var["cnv_type"] in ["Loss", "loss"]:
				result.append(var["gene_symbol"]+" 缺失")
			else:
				result.append(var["gene_symbol"]+" "+fdzs_mp_cnv_stran(var)) #修改CNV描述 2026.05.14, jmc
		elif var["bio_category"] in ["Sv", "PSeqRnaSv"]:
			if var["five_prime_gene"] == "MET" and var["three_prime_gene"] == "MET":
				result.append("MET exon14 跳跃")
			else:
				# 融合可能会有重复（rna exon相同，断点不同的情况）
				if var["five_prime_gene"]+"-"+var["three_prime_gene"]+" 融合" not in result:
					result.append(var["five_prime_gene"]+"-"+var["three_prime_gene"]+" 融合")
		# 2025.07.15-新增PHd
		elif var["bio_category"] == "PHd":
			if var["type"] == "HomoDel":
				if var["gene_symbol"]+" 纯合缺失" not in result:
					result.append(var["gene_symbol"]+" 纯合缺失")
			elif var["type"] == "HeteDel":
				if var["gene_symbol"]+" 杂合缺失" not in result:
					result.append(var["gene_symbol"]+" 杂合缺失")
			else:
				if var["gene_symbol"]+" 未知变异类型！" not in result:
					result.append(var["gene_symbol"]+" 未知变异类型！")
		# 2025.07.15-新增完成
	return ", ".join(result)
jinja2.filters.FILTERS["var_sum_s_fdzs"] = var_sum_s_fdzs

# 复旦中山MP免疫相关变异总结，嵇梦晨，2026.05.14
def io_summary_fdzs(io):
	result = []
	if io["io_p_summary"]:
		result.append(io["io_p_fdzs_mp"]+"（疗效正相关）")
	if io["io_n_summary"]:
		result.append(io["io_n_fdzs_mp"]+"（疗效负相关）")
	if not result:
		result = ["未检出相关变异"]
	rt = RichText()
	rt.add("；\n".join(result)+"。")
	return rt
jinja2.filters.FILTERS["io_summary_fdzs"] = io_summary_fdzs