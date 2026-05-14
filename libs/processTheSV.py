#-*- coding:gbk -*-
import re
from libs.getConfig import clinicalNumStran, functionNumStran
from libs.getEvi import varRegimen
import copy
from itertools import groupby
from libs.specialRequest import varInfo_FJZL
from libs.rule import S_function
from libs.rule import decimal_float, decimal_percen

'''
Discription	
	处理sv格式。
	包含sv 和rna_sv
	处理格式包含但不限于下列内容：
	1. 融合主基因判断
	   通过检测范围，5端基因在范围内，但3端不在，则主基因为5端基因，其余情况均默认以3端为主基因(取消，直接使用json返回的gene_symbol，2022.07.26)
	2. 基因合并，仅限于CP40 sv
	3. RNA SV 和DNA SV匹配 
'''

# 用于拼接融合：gene1:exon1-gene2:exon2
def var_info(var):
	return var["five_prime_gene"]+":"+var["five_prime_cds"]+"-"+var["three_prime_gene"]+":"+var["three_prime_cds"]

def process_sv(jsonDict, config):
	sv = copy.deepcopy(jsonDict["sv"])
	rna_sv = copy.deepcopy(jsonDict["rna_sv"])
	# sv处理
	for var in sv:
		# 2025.04.24-新增evi_sum_raw
		#var["evi_sum_raw"] = copy.deepcopy(var["evi_sum"]) if "evi_sum" in var.keys() and var["evi_sum"] else []
		# 2025.04.24-新增完成
		# 变异分类更新，致癌性和致病性只返回1个-2023.05.22
		#var["clinic_num_g"] = clinicalNumStran().get(var["clinical_significance"], 3) 
		#var["clinic_num_s"] = functionNumStran().get(var["function_classification"], 3) 
		var["clinic_num_g"] = clinicalNumStran(config).get(var["clinical_significance"], 3) if var["clinical_significance"] and var["clinical_significance"] != "-" else \
							  clinicalNumStran(config).get(var["function_classification"], 3)
		var["clinic_num_s"] = functionNumStran(config).get(var["function_classification"], 3) if var["function_classification"] and var["function_classification"] != "-" else \
							  functionNumStran(config).get(var["clinical_significance"], 3)
		# 2023.05.22 更新完成
		var["evi_sum"] = varRegimen(jsonDict, var["evi_sum"], config, var)
		var["clinic_num_s"], var["top_level"] = S_function(var)
		# var_hgvs兼容4种格式格式
		# three_prime_cds 和 five_prime_cds重新提取，返回json中删除了“-ins”等内容，适用格式如下
		## 1. gene:转录本_exon--gene:转录本_exon
		## 2. gene:转录本:exon--gene:转录本:exon
		var["five_prime_cds"] = "-".join(re.split("-", (re.split(":", var["var_hgvs"])[2]))[:-1]) \
								if not re.search("--", var["var_hgvs"]) else re.split("_", (re.split("--", var["var_hgvs"])[0]))[-1]
		var["three_prime_cds"] = re.split(":", var["var_hgvs"])[-1] if not re.search("--", var["var_hgvs"]) else re.split("_", (re.split("--", var["var_hgvs"])[1]))[-1]
		#print (var["five_prime_cds"], var["three_prime_cds"])

		# 加一个兼容-2023.10.19
		# var_hgvs新格式，3. gene1:NM_xxx:exon1--gene2:NM_xxx:exon2, 旧的为gene1:NM_xxx_exon1--gene2:NM_xxx_exon2
		# cds会变成xxx:exon1和xxx:exon2
		var["five_prime_cds"] = re.split(":", var["five_prime_cds"])[-1] if re.search(":", var["five_prime_cds"]) else var["five_prime_cds"]
		var["three_prime_cds"] = re.split(":", var["three_prime_cds"])[-1] if re.search(":", var["three_prime_cds"]) else var["three_prime_cds"]
		# 兼容完成-2023.10.19
		#print (var["five_prime_cds"], var["three_prime_cds"])

		# 加一个兼容-2024.01.25
		# 4. gene:转录本_exon-gene:转录本_exon 重新提取后，five_prime_cds为空，以此做为重新拆分的判定依据
		if not var["five_prime_cds"]:
			var["five_prime_cds"] = re.split("_", re.split("-", var["var_hgvs"])[0])[-1]
			var["three_prime_cds"] = re.split("_", var["var_hgvs"])[-1]
		# 兼容完成-2024.01.25

		# 主基因转录本
		tmp_dict = {var["five_prime_gene"] : var["five_prime_transcript"], var["three_prime_gene"] : var["three_prime_transcript"]}
		var["transcript_primary"] = tmp_dict.get(var["gene_symbol"], "")
		# 处理频率 CP40用copies，其余项目用freq
#		var["freq_str"] = "{:.2%}".format(float(var["freq"])) if "freq" in var.keys() and var["freq"] else ""

		# freq 兼容百分数格式，均转化为小数格式-适配v4-刘炜芬，2023.11.09
		freq_key = ["freq"]
		for key in freq_key:
			if key in var.keys():
				var[key] = float(var[key].replace("%", ""))/100 if re.search("%", str(var[key])) else var[key]
		# 兼容结束-2023.11.09

		var["freq_str"] = decimal_percen(var["freq"]) if "freq" in var.keys() and var["freq"] else ""
		var["copies"] = var["copies"] if "copies" in var.keys() and var["copies"] else ""
		# 福建肿瘤：返回变异频率相关信息（来源配置表）和治疗方案汇总
		var["var_info_forFJZL"], var["var_regimen_forFJZL"] = varInfo_FJZL(var, jsonDict["sample_info"]["tumor_list"], config)
		# 多个主基因原以“,”分隔，现有以“\n”分隔的情况，考虑代码中做兼容-2023.03.31
		var["gene_symbol"] = var["gene_symbol"].replace("\n", ",")
	# sv合并结果 仅CP40/CRC12要合并，其余项目直接使用sv的数据
	#sv_combination = combinationSV(sv) if re.search("Classic|CRC12", jsonDict["sample_info"]["prod_names"]) else sv
	# sv合并结果，新增OncoPro项目，目前仅CP40/CRC12/OncoPro要合并，其余项目直接使用sv数据-2024.07.31
	sv_combination = combinationSV(sv) if re.search("Classic|CRC12|OncoPro", jsonDict["sample_info"]["prod_names"]) else sv

	# rna_sv处理
	for var in rna_sv:
		#var["evi_sum_raw"] = copy.deepcopy(var["evi_sum"]) if "evi_sum" in var.keys() and var["evi_sum"] else []
		# 变异分类更新，致癌性和致病性只返回1个-2023.05.22
		#var["clinic_num_g"] = clinicalNumStran().get(var["clinical_significance"], 3) 
		#var["clinic_num_s"] = functionNumStran().get(var["function_classification"], 3) 
		var["clinic_num_g"] = clinicalNumStran(config).get(var["clinical_significance"], 3) if var["clinical_significance"] and var["clinical_significance"] != "-" else \
							  clinicalNumStran(config).get(var["function_classification"], 3)
		var["clinic_num_s"] = functionNumStran(config).get(var["function_classification"], 3) if var["function_classification"] and var["function_classification"] != "-" else \
							  functionNumStran(config).get(var["clinical_significance"], 3)
		# 2023.05.22 更新完成
		var["evi_sum"] = varRegimen(jsonDict, var["evi_sum"], config, var)
		var["clinic_num_s"], var["top_level"] = S_function(var)
		# 处理three_prime_cds 和 five_prime_cds
		var["five_prime_cds"] = re.split(":", var["five_prime_cds"])[0]
		var["three_prime_cds"] = re.split(":", var["three_prime_cds"])[0]
		# 频率（部分模板需要展示）
		var["freq"] = int(float(var["supp_splt_reads"]) + float(var["supp_span_reads"]))
		var["var_info_forFJZL"], var["var_regimen_forFJZL"] = varInfo_FJZL(var, jsonDict["sample_info"]["tumor_list"], config)
		# 多个主基因原以“,”分隔，现有以“\n”分隔的情况，考虑代码中做兼容-2023.03.31
		var["gene_symbol"] = var["gene_symbol"].replace("\n", ",")
		# 兼容v4，v3的reads字段v4不返回了-2023.11.06
		var["reads"] = int(float(var["supp_splt_reads"]) + float(var["supp_span_reads"])) if "reads" not in var.keys() else var["reads"]
		# reads处理为整数-2023.11.06
		var["reads"] = int(var["reads"]) if "reads" in var.keys() and var["reads"] else ""

		# 兼容完成-2023.11.06

	# rna sv 和dna sv匹配
	sv_combination_match, rna_sv_only = matchRNA_DNA_sv(sv_combination, rna_sv)

	# 按频率排序下
	# 2024.03.06-新增Myeloid
	# 2024.07.31-新增OncoPro
	if not re.search("Classic|CRC12|Myeloid|OncoPro", jsonDict["sample_info"]["prod_names"]):
#		print (sv_combination_match)
		sv_combination_match = sorted(sv_combination_match, key=lambda i:float(str(i["freq"]).replace("%","")), reverse=True)
	else:
		sv_combination_match = sorted(sv_combination_match, key=lambda i:float(i["copies"]), reverse=True)
	
	rna_sv = sorted(rna_sv, key=lambda i:float(i["freq"]), reverse=True)
	rna_sv_only = sorted(rna_sv_only, key=lambda i:float(i["freq"]), reverse=True)

	return sv_combination_match, rna_sv, rna_sv_only

def combinationSV(sv_data):
	# 2025.06.24-hgvs_p需要加到前面来，排序时hgvs_p一致的排到一起
	for var in sv_data:
		var["hgvs_p"] = var["five_prime_gene"] + "-" + var["three_prime_gene"]
		Predictive_evi = [evi["regimen_name"] for evi in var["evi_sum"]["evi_split"]["Predictive"]] if "Predictive" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Predictive"] else []
		Prognostic_evi = [evi["clinical_significance"] for evi in var["evi_sum"]["evi_split"]["Prognostic"]] if "Prognostic" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Prognostic"] else []
		Diagnostic_evi = ["Diagnostic" for evi in var["evi_sum"]["evi_split"]["Diagnostic"]] if "Diagnostic" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Diagnostic"] else []
		evi_key = []
		evi_key.extend(sorted(Predictive_evi))
		evi_key.extend(sorted(Prognostic_evi))
		evi_key.extend(sorted(Diagnostic_evi))
		var["hgvs_p"] += "-".join(evi_key)
	# 2025.06.24-新增完成
	# 加个排序，使两端基因相同的变异排到一起，groupby对于ABBA排序的变异处理有点问题-2023.05.25
	#sv_data = sorted(sv_data, key=lambda i:(i["five_prime_gene"], i["three_prime_gene"]))
	# 2025.06.24-排序新增hgvs_p
	sv_data = sorted(sv_data, key=lambda i:(i["five_prime_gene"], i["three_prime_gene"], i["hgvs_p"]))
	# 2025.06.24-更新完成
	# 2023.05.25-更新完成
	# 仅适用于CP40和CRC12
	#for var in sv_data:
		#var["hgvs_p"] = var["five_prime_gene"] + "-" + var["three_prime_gene"]

		# 2025.04.24-hgvs_p增加证据，防止不同基因型的融合证据不同但是被合并到一起了
		#Predictive_evi = [evi["regimen_name"] for evi in var["evi_sum_raw"] if evi["evidence_type"] == "Predictive"]
		#Prognostic_evi = [evi["clinical_significance"] for evi in var["evi_sum_raw"] if evi["evidence_type"] == "Prognostic"] 
		#Diagnostic_evi = ["Diagnostic" for evi in var["evi_sum_raw"] if evi["evidence_type"] == "Diagnostic"] 

		#Predictive_evi = [evi["regimen_name"] for evi in var["evi_sum"]["evi_split"]["Predictive"]] if "Predictive" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Predictive"] else []
		#Prognostic_evi = [evi["clinical_significance"] for evi in var["evi_sum"]["evi_split"]["Prognostic"]] if "Prognostic" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Prognostic"] else []
		#Diagnostic_evi = ["Diagnostic" for evi in var["evi_sum"]["evi_split"]["Prognostic"]] if "Diagnostic" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Diagnostic"] else []
		#evi_key = []
		#evi_key.extend(sorted(Predictive_evi))
		#evi_key.extend(sorted(Prognostic_evi))
		#evi_key.extend(sorted(Diagnostic_evi))
		#var["hgvs_p"] += "-".join(evi_key)
		# 2025.04.24-新增结束
	sv_sum = groupby(sv_data, key=lambda x : x["hgvs_p"])
	sv_combination = []
	for i, j in sv_sum:
		sv_group = list(j)
		major_sv = sv_group[0]
		# 添加各个融合型及对应的拷贝数，适用于云南肿瘤CP40-2022.08.23
		major_sv["sv_YNZL"] = []
		for k in sv_group:
			major_sv["sv_YNZL"].append({
				"five_prime_gene" : k["five_prime_gene"],
				"five_prime_cds" : k["five_prime_cds"],
				"three_prime_gene" : k["three_prime_gene"],
				"three_prime_cds" : k["three_prime_cds"],
				"copies" : k["copies"] if "copies" in k.keys() and k["copies"] else ""
				})
		# 添加完成-2022.08.23
		major_sv['copies'] = sum(int(k["copies"]) for k in sv_group if k["copies"])
		# 融合模式拼接
		major_sv["var_desc_merge"] = "、".join([var_info(k) for k in sv_group])
		# 各个融合型以列表的形式存，方便报告展示
		major_sv["merge_sv_list"] = [var_info(k) for k in sv_group]
	
		sv_combination.append(major_sv)

	return sv_combination

def matchRNA_DNA_sv(sv, rna_sv):
	rna_sv_dict = {}
	# 301 MP存在DNA1条+RNA多条（exon同，断点不同）的情况，这边做兼容，新增个字典-2023.04.12
	dup_rna_sv_dict = {}
	for var in rna_sv:
		rna_sv_dict.setdefault(var_info(var), {})
		rna_sv_dict[var_info(var)] = var
	# 新增字典-2023.04.12
		dup_rna_sv_dict.setdefault(var_info(var), [])
		dup_rna_sv_dict[var_info(var)].append(var)
	
	# 若存在RNA exon相同，断点不同的情况时，添加个标签用于判断是否展示具体断点-2023.04.12
	for var in rna_sv:
		key = var_info(var)
		if key in dup_rna_sv_dict.keys():
			if len(dup_rna_sv_dict[key]) > 1:
				var["dup_rna"] = "yes"
	
	rna_sv_only = copy.deepcopy(rna_sv)
	rna_sv_pop_key = []
	for var in sv:
		if var_info(var) in rna_sv_dict.keys():
			var["rna_detect"] = rna_sv_dict[var_info(var)]
			# 新增字典-2023.04.12
			var["dup_rna_detect"] = dup_rna_sv_dict[var_info(var)]
			rna_sv_pop_key.append(var_info(var))

	# 匹配完的，rna_sv_only把DNA/RNA共检变异删掉
	# 代码更新，原remove无法将重复的全部删除，改为下面的代码-2023.02.28
	rna_sv_only_reduce = []
	for var in rna_sv_only:
		if var_info(var) not in rna_sv_pop_key:
			rna_sv_only_reduce.append(var)

#	for var in rna_sv_only:
#		if var["five_prime_gene"]+":"+var["five_prime_cds"]+"-"+var["three_prime_gene"]+":"+var["three_prime_cds"] in rna_sv_pop_key:
#			rna_sv_only.remove(var)

	return sv, rna_sv_only_reduce
	