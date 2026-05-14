#-*- coding:gbk -*-
import os
import xlrd
import re
from libs.getConfig import CP40_RECOMGene, PAN116_RECOMGene, CP40Gene
import copy
from libs.rule import S_level, s_var_rule, judge_var, get_varsimpleinfo, judgeRegimen
# 2025.11.26-新增CP43配置信息
from libs.getConfig import CP43_RECOMGene
# 2025.11.26新增完成

'''
Discription
	
	伴随诊断推荐基因检测结果，格式多种。 
	
'''
def getNCCN_detect(var_data, tumor_list, mlpa_data, config, mlpa_data_v2, gcnv_data_v2):
	cdx_result = {}
	cdx_result["format1_forMP"] = format1(var_data)
	cdx_result["format2_forCP"] = format2(var_data, config)
	cdx_result["format3_forCP_splittumor"], cdx_result["format3_forCP_splittumor_simple"] = format3(var_data, tumor_list, config)
	cdx_result["format4_forHRDC"] = format4(var_data)
	cdx_result["format5_forHRR"] = format5(var_data, mlpa_data)
	# 新版报告模板-HRD和HRR新增两个基因-这边加一个版本-2023.12.04
	cdx_result["format4_forHRDC_for_new_vesion"] = format4_for_new_vesion(var_data)
	cdx_result["format5_forHRR_for_new_vesion"] = format5_for_new_vesion(var_data, mlpa_data)
	# 2023.12.04-新增完成
	# 2025.01.15-新增HRR MLPA/CNV汇总表
	cdx_result["format5_forHRR_for_new_vesion_mlpa"] = format5_for_new_vesion_mlpa(var_data, mlpa_data_v2)
	cdx_result["format5_forHRR_for_new_vesion_gcnv"] = format5_for_new_vesion_gcnv(var_data, gcnv_data_v2)
	# 2025.01.15-新增完成
	# 2026.03.16-新增HRR 34基因，增加MLH1基因
	cdx_result["format5_forHRR_for_new_vesion_mlpa_34"] = format5_for_new_vesion_mlpa_34(var_data, mlpa_data_v2)
	cdx_result["format5_forHRR_for_new_vesion_gcnv_34"] = format5_for_new_vesion_gcnv_34(var_data, gcnv_data_v2)
	# 2026.03.16-新增完成
	cdx_result["format6_for116"] = format6(var_data, config)
	cdx_result["format7_for116_allVar"] = format7(var_data, config)
	# 新增福建省立CP40，III类变异删掉同义突变-20221013
	var_data_forFJSL = copy.deepcopy(var_data)
	FJSL_var = s_var_rule(var_data_forFJSL)
	FJSL_var_for_summary = FJSL_var["level_I"]+FJSL_var["level_II"]+FJSL_var["level_onco_nodrug"]+FJSL_var["level_III_without_Syn"]
	cdx_result["format8_forCP_FJSL_without_Syn"], cdx_result["format8_forCP_FJSL_without_Syn_simple"]= format3(FJSL_var_for_summary, tumor_list, config)
	# 新增重庆西南，III类变异删掉非编码区（除剪接）且解读为3的变异-2023.07.13
	var_data_forCQXN = copy.deepcopy(var_data)
	var_CQXN_filter = []
	for var in var_data_forCQXN:
		if var["bio_category"] in ["Sv", "Cnv"]:
			var_CQXN_filter.append(var)
		elif var["bio_category"] == "Snvindel":
			if var["clinic_num_s"] == 3 and var["type"] in ["Intronic", "3'UTR", "5'UTR", "FlankingRegion3", "FlankingRegion5"] and not set(["Diagnostic","Predictive","Prognostic"]) & set(var["evi_sum"]["evi_split"].keys()):
				pass
			else:
				var_CQXN_filter.append(var)
	cdx_result["format9_forCP_CQXN"], cdx_result["format9_forCP_CQXN_simple"] = format3(var_CQXN_filter, tumor_list, config)

	# 合并梦晨代码-2024.01.11
	cdx_result["format7_forTC30FJSL"] = format7_FJSL(var_data)
	cdx_result["format8_forTC30FJXH"] = format8_FJXH(var_data)
	# 合并梦晨代码结束-2024.01.11

	# 2024.06.06-福建省立116，展示体细胞I/II类变异
	cdx_result["format9_forFJSL_116"] = format9(var_data, config)
	# 2024.06.06-新增结束

	# 2025.03.03-云南肿瘤CP40，展示I/II类变异
	cdx_result["format10_forYNZL_CP40"] = format_ynzl_cp40(var_data, tumor_list, config)

	# 2025.11.26-新增CP43
	cdx_result["format3_forCP43_splittumor"], cdx_result["format3_forCP43_splittumor_simple"] = format3_cp43(var_data, tumor_list, config)
	# 2025.11.26-新增CP43结束

	return cdx_result

# 各种规则判断######################################################################################################################################
# 1. 返回SV的第二种展示方式，gene1:exon1-gene2:exon2
def get_svinfo2(var):
	return "{0}:{1}-{2}:{3}融合".format(var["five_prime_gene"], var["five_prime_cds"], var["three_prime_gene"], var["three_prime_cds"])

# 2. 返回MP临检的变异频率
def get_freq(var):
	if var["bio_category"] == "Snvindel":
		return var["freq_str"]
	elif var["bio_category"] == "Cnv":
		return var["cn_mean"]
	elif var["bio_category"] == "Sv":
		# DNA/RNA共检时，展示RNA的频率
		return "{0} copies".format(str(var["rna_detect"]["freq"])) if "rna_detect" in var.keys() and var["rna_detect"] else var["freq_str"]
	elif var["bio_category"] == "PSeqRnaSv":
		return "{0} copies".format(str(var["freq"]))

# 3. 判断变异是否为I类
def judge_level_I(var):
	if var["evi_sum"]["evi_split"] and set(["Diagnostic","Predictive","Prognostic"]) & set(var["evi_sum"]["evi_split"].keys()) and var["clinic_num_s"] == 5:
		return 1
	else:
		return 0

# 4. 返回CP40的变异频率
def get_freq_cp(var):
	if var["bio_category"] == "Snvindel":
		return var["freq_str"]
	elif var["bio_category"] == "Cnv":
		return var["cn_mean"]
	else:
		if "freq_str" in var.keys() and var["freq_str"]:
			return var["freq_str"]
		elif "copies" in var.keys() and var["copies"]:
			return var["copies"]
		elif "reads" in var.keys() and var["reads"]:
			return var["reads"]

# 判定规则结束###############################################################################################################################################


# 格式1：用于MP临检通用、浙肿MP和CP40，NCCN指南推荐基因
# 需要判断检测内容，如点突变、小片段插入缺失、融合、MET14跳跃、KRAS G12等
def format1(var_data):
	var_data_copy = copy.deepcopy(var_data)
	# 新增snvindel
	# snvindel 新增ERBB2 - 2022.08.30
	# snvindel 新增POLD1、POLE、AKT1、PTEN和VHL - 2025.07.01
	snvindel_gene = ["EGFR", "KIT", "PDGFRA","BRCA1","BRCA2","ATM","BARD1","BRIP1","CDH1","CDK12","CHEK1","CHEK2",\
					 "FANCA","FANCL","HDAC2","PALB2","PPP2R2A","PTEN","RAD51B","RAD51C","RAD51D","RAD54L","TP53", "ERBB2", \
					 "POLD1", "POLE", "AKT1", "PTEN", "VHL"]
	# snv 新增POLE和POLD1
	# snv 新增ESR1-2025.07.01
	snv_gene = ["ALK", "ROS1", "RET", "BRAF", "ERBB2", "PIK3CA", "FGFR2", "FGFR3", "IDH1", "IDH2","BRCA1","BRCA2", "POLE", "POLD1", "ESR1"]
	fusion_gene = ["ALK", "ROS1", "RET", "FGFR2", "FGFR3", "NTRK1", "NTRK2", "NTRK3"]
	cnv_gene = ["MET", "ERBB2"]
	other_gene = ["MET", "KRAS", "NRAS"]
	result = {}
	# 只展示I/II/肿瘤发生发展变异和胚系4/5类变异
	level_12_var = [var for var in var_data_copy if judge_var(var, [4,5], [4,5])]
	# 2025.02.21- 更新展示变异等级
	# 1. 有配对：体细胞I/II/肿瘤发生发展 + 胚系4/5 + 胚系3类但有用药
	# 2. 无配对：体细胞I/II/肿瘤发生发展 + 预测胚系4/5（归为I/II/肿瘤发生发展） + 预测胚系3类但有用药
	# ==> 即新增（确认/预测）胚系3类但有用药的变异即可
	germline_level3_regimen = [var for var in var_data_copy if var["var_origin"] == "germline" and var["clinic_num_g"] == 3 and judgeRegimen(var)]
	if germline_level3_regimen:
		level_12_var.extend(germline_level3_regimen)
	# 2025.02.21-新增完成
	for i in level_12_var:
		# 新增几个字段用于报告展示
		i["var_info"] = get_varsimpleinfo(i)
		i["var_info_M"] = get_svinfo2(i) if i["bio_category"] in ["Sv", "PSeqRnaSv"] else ""
		i["freq"] = get_freq(i)
		# 下面这行好像重复了，有时间看下
		i["var"] = copy.deepcopy(i)
		# 结果处理
		# 1. 检测snvindel基因的结果
		if i["bio_category"] == "Snvindel" and i["gene_symbol"] in snvindel_gene:
			if i["gene_symbol"]+"_snvindel" not in result.keys():
				result.setdefault(i["gene_symbol"]+"_snvindel", []) 
			result[i["gene_symbol"]+"_snvindel"].append(i)
		# 2. 检测snv基因的结果
		if i["bio_category"] == "Snvindel" and i["gene_symbol"] in snv_gene and not re.search("del|ins|fs", i["hgvs_p"]):
			if i["gene_symbol"]+"_snv" not in result.keys():
				result.setdefault(i["gene_symbol"]+"_snv", [])
			result[i["gene_symbol"]+"_snv"].append(i)		
		# 3. 检测cnv基因的结果
		if i["bio_category"] == "Cnv" and i["gene_symbol"] in cnv_gene:
			if i["gene_symbol"]+"_cn" not in result.keys():
				result.setdefault(i["gene_symbol"]+"_cn", [])
			result[i["gene_symbol"]+"_cn"].append(i)		
		# 4. 检测sv基因的结果
		if i["bio_category"] in ["Sv", "PSeqRnaSv"] and set(re.split(",", i["gene_symbol"])) & set(fusion_gene):
			for sv_gene in set(re.split(",", i["gene_symbol"])):
				if sv_gene+"_sv" not in result.keys():
					result.setdefault(sv_gene+"_sv", [])
				result[sv_gene+"_sv"].append(i)		
		# 5. 特殊检测内容处理
		# 5.1 KRAS、NRAS 分为G12|G13|Q61|A146和其他变异
		if i["bio_category"] == "Snvindel" and i["gene_symbol"] in ["KRAS", "NRAS"]:
			if re.search("G12|G13|Q61|A146", i["hgvs_p"]):
				if i["gene_symbol"]+"_sp" not in result.keys():
					result.setdefault(i["gene_symbol"]+"_sp", [])
				result[i["gene_symbol"]+"_sp"].append(i)
			else:
				if i["gene_symbol"]+"_ot" not in result.keys():
					result.setdefault(i["gene_symbol"]+"_ot", [])
				result[i["gene_symbol"]+"_ot"].append(i)
		# 5.2 MET 14跳跃突变，仅存在snvindel时，无法准确判定，只能匹配I类MET变异，仅存在sv时，判定为14跳跃，同时存在snvindel和sv时，判定为14跳跃
		if (i["gene_symbol"] == "MET" and i["bio_category"] == "Snvindel" and judge_level_I(i)) or \
		   ("judge_mergeMET" in i.keys() and i["judge_mergeMET"]):
			if "MET_skip" not in result.keys():
				result.setdefault("MET_skip", [])
			result["MET_skip"].append(i)
			# MET 14跳跃-新增MET-MET-20220901
		if i["bio_category"] in ["Sv", "PSeqRnaSv"] and i["five_prime_gene"] == "MET" and i["three_prime_gene"] == "MET":
			if "MET_skip" not in result.keys():
				result.setdefault("MET_skip", [])
			result["MET_skip"].append(i)

	# 检测结果去重（正常情况下应该不会有重复的）	
	sort_result = {}
	for k, v in result.items():
		sort_result[k] = []
		for i in v:
			if i not in sort_result[k]:
				sort_result[k].append(i)
	return sort_result

# 格式2：用于CP40
# 展示基因、转录本、相关肿瘤、检测结果，不分区癌种，包含I/II/III类变异
def format2(var_data, config):
	var_data_copy = copy.deepcopy(var_data)
	recom_list = CP40_RECOMGene(config)["RECORDS"]
	result = []
	# 处理recom_list
	recom_dict = {}
	for i in recom_list:
		if i["gene_symbol"] not in recom_dict.keys():
			recom_dict.setdefault(i["gene_symbol"], {})
		disease = "、".join([j["disease"] for j in recom_list if i["gene_symbol"] == j["gene_symbol"]])
		
		recom_dict[i["gene_symbol"]] = {
			"transcript" : i["transcript_refseq"],
			"disease" : disease,
			"var_type" : i["var_type"]
		}
	recom_gene = [i for i in recom_dict.keys()]
	detect_gene = []
	# 返回检测到变异的指南推荐基因
	var_list = [var for var in var_data_copy if var["clinic_num_s"] in [3,4,5]]
	for var in var_list:
		for gene in set(re.split(",", var["gene_symbol"])):
			if gene in recom_dict.keys():
				detect_gene.append(gene)
				var["var_info"] = get_varsimpleinfo(var)
				var["note"] = var["judge_mergeMET"] if "judge_mergeMET" in var.keys() else ""
				var["note_freq"] = var["freq_2"] if "judge_mergeMET" in var.keys() else ""
				var["freq"] = get_freq_cp(var)
				var["level"] = str(int(S_level(var)))
				var["gene_symbol"] = gene
				var["transcript"] = recom_dict[gene]["transcript"]
				var["disease"] = recom_dict[gene]["disease"]
				var["var_type"] = recom_dict[gene]["var_type"]
				result.append(var)
	# 返回未检测到变异的指南推荐基因			
	for gene in set(recom_gene) - set(detect_gene):
		result.append({
			"gene_symbol" : gene,
			"transcript" : recom_dict[gene]["transcript"],
			"disease" : recom_dict[gene]["disease"],
			"var_type" : recom_dict[gene]["var_type"]
		})
	return sorted(result, key=lambda i:i["gene_symbol"])

# 格式3：用于CP40
# 展示基因、检测内容和结果，区分癌种，包含I/II/III类变异
def format3(var_data, tumor_list, config):
	var_data_copy = copy.deepcopy(var_data)
	recom_list = CP40_RECOMGene(config)["RECORDS"]
	# 将基因对应所有癌种合并起来
	gene_tumor = {}
	for i in recom_list:
		if i["gene_symbol"] not in gene_tumor.keys():
			gene_tumor.setdefault(i["gene_symbol"], [])
		gene_tumor[i["gene_symbol"]].append(i["disease"])
	# 筛选出对应癌种的NCCN推荐基因，其中均包含实体瘤
	tumor_recom_list = [i for i in recom_list if i["disease"] in tumor_list or i["disease"] == "实体瘤"]
	recom_dict = {}
	for i in tumor_recom_list:
		if i["gene_symbol"] not in recom_dict.keys():
			recom_dict[i["gene_symbol"]] = {
				"transcript" : i["transcript_refseq"],
				"disease" : i["disease"],
				"var_type" : i["var_type"],
				"gene_tumor" : "、".join(gene_tumor.get(i["gene_symbol"]))
			}
	# result用于返回肿瘤相关的全部基因，result_simple只返回检测到变异的基因，如贵州肿瘤CP40
	result = []
	result_simple = []
	recom_gene = [i for i in recom_dict.keys()]
	detect_gene = []
	# 返回检测到变异的指南推荐基因
	var_list =[var for var in var_data_copy if var["clinic_num_s"] in [3,4,5]]
	for var in var_list:
		for gene in set(re.split(",", var["gene_symbol"])):
			if gene in recom_dict.keys():
				detect_gene.append(gene)
				var["var_info"] = get_varsimpleinfo(var)
				var["note"] = var["judge_mergeMET"] if "judge_mergeMET" in var.keys() else ""
				var["note_freq"] = var["freq_2"] if "judge_mergeMET" in var.keys() else ""
				var["freq"] = get_freq_cp(var)			
				# 加一个融合的分型列表 2022.07.26
				# 加一个融合分型及拷贝数列表 2022.08.23 适用云南肿瘤CP40
				var["merge_sv_list"] = var["merge_sv_list"] if "merge_sv_list" in var.keys() else []
				var["sv_YNZL"] = var["sv_YNZL"] if "sv_YNZL" in var.keys() else []
				var["level"] = str(int(S_level(var)))
				var["gene_symbol"] = gene
				var["transcript"] = recom_dict[gene]["transcript"]
				var["disease"] = recom_dict[gene]["disease"]
				var["var_type"] = recom_dict[gene]["var_type"]
				var["gene_tumor"] = recom_dict[gene]["gene_tumor"]
				#print (var["transcript"], var["disease"], var["gene_tumor"])
				# 加一个氨基酸三字母
				var["var_info_abbr"] = get_varsimpleinfo(var) if var["bio_category"] != "Snvindel" else var["hgvs_p_abbr"] if var["hgvs_p_abbr"] != "p.?" else var["hgvs_c"]
				result.append(var)
				result_simple.append(var)
	# 返回未检测到变异的指南推荐基因-仅result
	for gene in set(recom_gene) - set(detect_gene):
		result.append({
			"gene_symbol" : gene,
			"transcript" : recom_dict[gene]["transcript"],
			"disease" : recom_dict[gene]["disease"],
			"var_type" : recom_dict[gene]["var_type"],
			"gene_tumor" : recom_dict[gene]["gene_tumor"]
		})

	return sorted(result, key=lambda i:i["gene_symbol"]), sorted(result_simple, key=lambda i:i["gene_symbol"])

# 2025.11.26-新增
# 格式3-CP43：用于CP43
# 展示基因、检测内容和结果，区分癌种，包含I/II/III类变异
def format3_cp43(var_data, tumor_list, config):
	var_data_copy = copy.deepcopy(var_data)
	recom_list = CP43_RECOMGene(config)["RECORDS"]
	# 将基因对应所有癌种合并起来
	gene_tumor = {}
	for i in recom_list:
		if i["gene_symbol"] not in gene_tumor.keys():
			gene_tumor.setdefault(i["gene_symbol"], [])
		gene_tumor[i["gene_symbol"]].append(i["disease"])
	# 筛选出对应癌种的NCCN推荐基因，其中均包含实体瘤
	tumor_recom_list = [i for i in recom_list if i["disease"] in tumor_list or i["disease"] == "实体瘤"]
	recom_dict = {}
	for i in tumor_recom_list:
		if i["gene_symbol"] not in recom_dict.keys():
			recom_dict[i["gene_symbol"]] = {
				"transcript" : i["transcript_refseq"],
				"disease" : i["disease"],
				"var_type" : i["var_type"],
				"gene_tumor" : "、".join(gene_tumor.get(i["gene_symbol"]))
			}
	# result用于返回肿瘤相关的全部基因，result_simple只返回检测到变异的基因，如贵州肿瘤CP40
	result = []
	result_simple = []
	recom_gene = [i for i in recom_dict.keys()]
	detect_gene = []
	# 返回检测到变异的指南推荐基因
	var_list =[var for var in var_data_copy if var["clinic_num_s"] in [3,4,5]]
	for var in var_list:
		for gene in set(re.split(",", var["gene_symbol"])):
			if gene in recom_dict.keys():
				detect_gene.append(gene)
				var["var_info"] = get_varsimpleinfo(var)
				var["note"] = var["judge_mergeMET"] if "judge_mergeMET" in var.keys() else ""
				var["note_freq"] = var["freq_2"] if "judge_mergeMET" in var.keys() else ""
				var["freq"] = get_freq_cp(var)			
				# 加一个融合的分型列表 2022.07.26
				# 加一个融合分型及拷贝数列表 2022.08.23 适用云南肿瘤CP40
				var["merge_sv_list"] = var["merge_sv_list"] if "merge_sv_list" in var.keys() else []
				var["sv_YNZL"] = var["sv_YNZL"] if "sv_YNZL" in var.keys() else []
				var["level"] = str(int(S_level(var)))
				var["gene_symbol"] = gene
				var["transcript"] = recom_dict[gene]["transcript"]
				var["disease"] = recom_dict[gene]["disease"]
				var["var_type"] = recom_dict[gene]["var_type"]
				var["gene_tumor"] = recom_dict[gene]["gene_tumor"]
				#print (var["transcript"], var["disease"], var["gene_tumor"])
				# 加一个氨基酸三字母
				var["var_info_abbr"] = get_varsimpleinfo(var) if var["bio_category"] != "Snvindel" else var["hgvs_p_abbr"] if var["hgvs_p_abbr"] != "p.?" else var["hgvs_c"]
				result.append(var)
				result_simple.append(var)
	# 返回未检测到变异的指南推荐基因-仅result
	for gene in set(recom_gene) - set(detect_gene):
		result.append({
			"gene_symbol" : gene,
			"transcript" : recom_dict[gene]["transcript"],
			"disease" : recom_dict[gene]["disease"],
			"var_type" : recom_dict[gene]["var_type"],
			"gene_tumor" : recom_dict[gene]["gene_tumor"]
		})

	return sorted(result, key=lambda i:i["gene_symbol"]), sorted(result_simple, key=lambda i:i["gene_symbol"])

# 格式4：用于HRD Complete，包含I/II/III类变异，无SV，所以不用考虑gene_symbol两个的情况
# 2024.02.20-所有模板均加上FANCA基因
def format4(var_data):
	var_data_copy = copy.deepcopy(var_data)
	gene_list = ["ATM","BARD1","BRCA1","BRCA2","BRIP1","CDK12","CHEK1","CHEK2","FANCL","PALB2","RAD51B","RAD51C","RAD51D","RAD54L","FANCA"]
	var_list = [var for var in var_data_copy if var["clinic_num_s"] in [3,4,5] and var["gene_symbol"] in gene_list and var["bio_category"] == "Snvindel"]
	detect_gene = [var["gene_symbol"] for var in var_list if var["gene_symbol"] in gene_list]
	result = []
	for var in var_list:
		var["clinic_num_s"] = S_level(var)
		result.append(var)
	for gene in set(gene_list) - set(detect_gene):
		result.append({
			"gene_symbol" : gene
		})	
	return sorted(result, key=lambda i:i["gene_symbol"])

# 新版报告模板汇总，需要新增两个基因结果（FANCA和HDAC2），这边添加一个更新版本-2023.12.04
# 格式4：用于HRD Complete，包含I/II/III类变异，无SV，所以不用考虑gene_symbol两个的情况
# 2024.01.16-删除HDAC2基因
def format4_for_new_vesion(var_data):
	var_data_copy = copy.deepcopy(var_data)
	gene_list = ["ATM","BARD1","BRCA1","BRCA2","BRIP1","CDK12","CHEK1","CHEK2","FANCL","PALB2","RAD51B","RAD51C","RAD51D","RAD54L", "FANCA"]
	var_list = [var for var in var_data_copy if var["clinic_num_s"] in [3,4,5] and var["gene_symbol"] in gene_list and var["bio_category"] == "Snvindel"]
	detect_gene = [var["gene_symbol"] for var in var_list if var["gene_symbol"] in gene_list]
	result = []
	for var in var_list:
		var["clinic_num_s"] = S_level(var)
		result.append(var)
	for gene in set(gene_list) - set(detect_gene):
		result.append({
			"gene_symbol" : gene
		})	
	return sorted(result, key=lambda i:i["gene_symbol"])
# 2023.12.04-新增完成

# 格式5：用于HRR全血、组织和配对，包含体细胞I/II/III类变异和胚系3/4/5类变异，还需考虑MLPA的情况
# 2024.02.20-所有模板均加上基因ATR、FANCA、MRE11、NBN
def format5(var_data, mlpa_data):
	var_data_copy = copy.deepcopy(var_data)
	mlpa_data_copy = copy.deepcopy(mlpa_data)
	gene_list = ["ATM","BARD1","BRCA1","BRCA2","BRIP1","CDK12","CHEK1","CHEK2","FANCL","PALB2","RAD51B","RAD51C","RAD51D","RAD54L", "FANCA", "ATR", "MRE11", "NBN"]
	var_list = [var for var in var_data_copy if judge_var(var, [3,4,5], [3,4,5]) and var["gene_symbol"] in gene_list and var["bio_category"] == "Snvindel"]
	detect_gene = [var["gene_symbol"] for var in var_list if var["gene_symbol"] in gene_list]
	for k, v in mlpa_data_copy.items():
		if v:
			var_list += v
			detect_gene+= [var["gene_symbol"] for var in v]
	result = []
	for var in var_list:
		# MLPA结果
		if "type" in var.keys() and var["type"] and var["bio_category"] != "Snvindel" and re.search("Loss|Gain",var["type"]):
			result.append(var)
		# Snvindel
		else:
			var["clinic_num_s"] = S_level(var)
			result.append(var)				
	for gene in set(gene_list) - set(detect_gene):
		result.append({
			"gene_symbol" : gene
		})	
	return sorted(result, key=lambda i:i["gene_symbol"])

# 新版报告模板汇总，需要新增两个基因结果（FANCA和HDAC2），这边添加一个更新版本-2023.12.04
# 格式5：用于HRR全血、组织和配对，包含体细胞I/II/III类变异和胚系3/4/5类变异，还需考虑MLPA的情况
# 2024.01.04-新增ATR、MRE11、NBN基因
# 2024.01.16-删除HDAC2基因
def format5_for_new_vesion(var_data, mlpa_data):
	var_data_copy = copy.deepcopy(var_data)
	mlpa_data_copy = copy.deepcopy(mlpa_data)
	gene_list = ["ATM","BARD1","BRCA1","BRCA2","BRIP1","CDK12","CHEK1","CHEK2","FANCL","PALB2","RAD51B","RAD51C","RAD51D","RAD54L", "FANCA", "ATR", "MRE11", "NBN"]
	var_list = [var for var in var_data_copy if judge_var(var, [3,4,5], [3,4,5]) and var["gene_symbol"] in gene_list and var["bio_category"] == "Snvindel"]
	detect_gene = [var["gene_symbol"] for var in var_list if var["gene_symbol"] in gene_list]
	for k, v in mlpa_data_copy.items():
		if v:
			var_list += v
			detect_gene+= [var["gene_symbol"] for var in v]
	result = []
	for var in var_list:
		# MLPA结果
		if "type" in var.keys() and var["type"] and var["bio_category"] != "Snvindel" and re.search("Loss|Gain",var["type"]):
			result.append(var)
		# Snvindel
		else:
			var["clinic_num_s"] = S_level(var)
			result.append(var)				
	for gene in set(gene_list) - set(detect_gene):
		result.append({
			"gene_symbol" : gene
		})	
	return sorted(result, key=lambda i:i["gene_symbol"])
# 2023.12.04-新增完成

# 格式5：用于HRR全血、组织和配对，包含体细胞I/II/III类变异和胚系3/4/5类变异，还需考虑MLPA的情况
# 2025.01.15-MLPA等级改为动态提取的
def format5_for_new_vesion_mlpa(var_data, mlpa_data_v2):
	var_data_copy = copy.deepcopy(var_data)
	mlpa_data_v2_copy = copy.deepcopy(mlpa_data_v2)
	gene_list = ["ATM","BARD1","BRCA1","BRCA2","BRIP1","CDK12","CHEK1","CHEK2","FANCL","PALB2","RAD51B","RAD51C","RAD51D","RAD54L", "FANCA", "ATR", "MRE11", "NBN"]
	var_list = [var for var in var_data_copy if judge_var(var, [3,4,5], [3,4,5]) and var["gene_symbol"] in gene_list and var["bio_category"] == "Snvindel"]
	detect_gene = [var["gene_symbol"] for var in var_list if var["gene_symbol"] in gene_list]
	for var in mlpa_data_v2_copy["B1_mlpa_L5"] + mlpa_data_v2_copy["B2_mlpa_L5"] + \
			   mlpa_data_v2_copy["B1_mlpa_L4"] + mlpa_data_v2_copy["B2_mlpa_L4"] + \
			   mlpa_data_v2_copy["B1_mlpa_L3"] + mlpa_data_v2_copy["B2_mlpa_L3"]:
		var_list.append(var)
		detect_gene.append(var["gene_symbol"])
	result = []
	for var in var_list:
		#print (var)
		# MLPA结果
		if "type" in var.keys() and var["type"] and var["bio_category"] != "Snvindel" and re.search("Loss|Gain",var["type"]):
			result.append(var)
		# Snvindel
		else:
			var["clinic_num_s"] = S_level(var)
			result.append(var)				
	for gene in set(gene_list) - set(detect_gene):
		result.append({
			"gene_symbol" : gene
		})	
	return sorted(result, key=lambda i:i["gene_symbol"])
# 2025.01.15-新增完成

# 格式5：用于HRR全血、组织和配对，包含体细胞I/II/III类变异和胚系3/4/5类变异，还需考虑MLPA的情况
# 2025.01.15-MLPA等级改为动态提取的
# 2026.03.16-新增MLH1，适用34基因
def format5_for_new_vesion_mlpa_34(var_data, mlpa_data_v2):
	var_data_copy = copy.deepcopy(var_data)
	mlpa_data_v2_copy = copy.deepcopy(mlpa_data_v2)
	gene_list = ["ATM","BARD1","BRCA1","BRCA2","BRIP1","CDK12","CHEK1","CHEK2","FANCL","PALB2","RAD51B","RAD51C","RAD51D","RAD54L", "FANCA", "ATR", "MRE11", "NBN", "MLH1"]
	var_list = [var for var in var_data_copy if judge_var(var, [3,4,5], [3,4,5]) and var["gene_symbol"] in gene_list and var["bio_category"] == "Snvindel"]
	detect_gene = [var["gene_symbol"] for var in var_list if var["gene_symbol"] in gene_list]
	for var in mlpa_data_v2_copy["B1_mlpa_L5"] + mlpa_data_v2_copy["B2_mlpa_L5"] + \
			   mlpa_data_v2_copy["B1_mlpa_L4"] + mlpa_data_v2_copy["B2_mlpa_L4"] + \
			   mlpa_data_v2_copy["B1_mlpa_L3"] + mlpa_data_v2_copy["B2_mlpa_L3"]:
		var_list.append(var)
		detect_gene.append(var["gene_symbol"])
	result = []
	for var in var_list:
		#print (var)
		# MLPA结果
		if "type" in var.keys() and var["type"] and var["bio_category"] != "Snvindel" and re.search("Loss|Gain",var["type"]):
			result.append(var)
		# Snvindel
		else:
			var["clinic_num_s"] = S_level(var)
			result.append(var)				
	for gene in set(gene_list) - set(detect_gene):
		result.append({
			"gene_symbol" : gene
		})	
	return sorted(result, key=lambda i:i["gene_symbol"])
# 2025.01.15-新增完成

# 格式5：用于HRR全血、组织和配对，包含体细胞I/II/III类变异和胚系3/4/5类变异，还需考虑CNV的情况
# 2025.01.15-CNV等级改为动态提取的
def format5_for_new_vesion_gcnv(var_data, gcnv_data_v2):
	var_data_copy = copy.deepcopy(var_data)
	gcnv_data_v2_copy = copy.deepcopy(gcnv_data_v2)
	gene_list = ["ATM","BARD1","BRCA1","BRCA2","BRIP1","CDK12","CHEK1","CHEK2","FANCL","PALB2","RAD51B","RAD51C","RAD51D","RAD54L", "FANCA", "ATR", "MRE11", "NBN"]
	var_list = [var for var in var_data_copy if judge_var(var, [3,4,5], [3,4,5]) and var["gene_symbol"] in gene_list and var["bio_category"] == "Snvindel"]
	detect_gene = [var["gene_symbol"] for var in var_list if var["gene_symbol"] in gene_list]
	for var in gcnv_data_v2_copy["B1_gcnv_L5"] + gcnv_data_v2_copy["B2_gcnv_L5"] + \
			   gcnv_data_v2_copy["B1_gcnv_L4"] + gcnv_data_v2_copy["B2_gcnv_L4"] + \
			   gcnv_data_v2_copy["B1_gcnv_L3"] + gcnv_data_v2_copy["B2_gcnv_L3"]:
		var_list.append(var)
		detect_gene.append(var["gene_symbol"])
	result = []
	for var in var_list:
		# MLPA结果
		if "type" in var.keys() and var["type"] and var["bio_category"] != "Snvindel" and re.search("Loss|Gain",var["type"]):
			result.append(var)
		# Snvindel
		else:
			var["clinic_num_s"] = S_level(var)
			result.append(var)				
	for gene in set(gene_list) - set(detect_gene):
		result.append({
			"gene_symbol" : gene
		})	
	return sorted(result, key=lambda i:i["gene_symbol"])
# 2025.01.15-新增完成

# 格式5：用于HRR全血、组织和配对，包含体细胞I/II/III类变异和胚系3/4/5类变异，还需考虑CNV的情况
# 2025.01.15-CNV等级改为动态提取的
# 2026.03.16-新增MLH1，适用34基因
def format5_for_new_vesion_gcnv_34(var_data, gcnv_data_v2):
	var_data_copy = copy.deepcopy(var_data)
	gcnv_data_v2_copy = copy.deepcopy(gcnv_data_v2)
	gene_list = ["ATM","BARD1","BRCA1","BRCA2","BRIP1","CDK12","CHEK1","CHEK2","FANCL","PALB2","RAD51B","RAD51C","RAD51D","RAD54L", "FANCA", "ATR", "MRE11", "NBN", "MHL1"]
	var_list = [var for var in var_data_copy if judge_var(var, [3,4,5], [3,4,5]) and var["gene_symbol"] in gene_list and var["bio_category"] == "Snvindel"]
	detect_gene = [var["gene_symbol"] for var in var_list if var["gene_symbol"] in gene_list]
	for var in gcnv_data_v2_copy["B1_gcnv_L5"] + gcnv_data_v2_copy["B2_gcnv_L5"] + \
			   gcnv_data_v2_copy["B1_gcnv_L4"] + gcnv_data_v2_copy["B2_gcnv_L4"] + \
			   gcnv_data_v2_copy["B1_gcnv_L3"] + gcnv_data_v2_copy["B2_gcnv_L3"]:
		var_list.append(var)
		detect_gene.append(var["gene_symbol"])
	result = []
	for var in var_list:
		# MLPA结果
		if "type" in var.keys() and var["type"] and var["bio_category"] != "Snvindel" and re.search("Loss|Gain",var["type"]):
			result.append(var)
		# Snvindel
		else:
			var["clinic_num_s"] = S_level(var)
			result.append(var)				
	for gene in set(gene_list) - set(detect_gene):
		result.append({
			"gene_symbol" : gene
		})	
	return sorted(result, key=lambda i:i["gene_symbol"])
# 2025.01.15-新增完成


# 格式6：用于PAN116单独组织、配对组织和全血，包含体细胞I/II/肿瘤发生发展相关和胚系4/5类变异
def format6(var_data, config):
	var_data_copy = copy.deepcopy(var_data)
	recom_dict = {}
	for i in PAN116_RECOMGene(config)["RECORDS"]:
		recom_dict[i["gene_symbol"]] = i["disease"]
	recom_gene = [i for i in recom_dict.keys()]
	var_list = [var for var in var_data_copy if judge_var(var, [4,5], [4,5])]
	result = []
	detect_gene = []
	for var in var_list:
		for gene in set(re.split(",", var["gene_symbol"])):
			if gene in recom_gene:
				detect_gene.append(gene)
				var["gene_symbol"] = gene
				var["disease"] = recom_dict[gene]
				var["clinic_num_s"] = S_level(var)
				result.append(var)
	for gene in set(recom_gene) - set(detect_gene):
		result.append({
			"gene_symbol" : gene,
			"disease" : recom_dict[gene]
		})
	return sorted(result, key=lambda i:i["gene_symbol"])

# 格式7：用于PAN116单独组织、配对组织和全血，包含体细胞I/II/III/肿瘤发生发展相关和胚系4/5类变异
# 适用模板：西安交大一
def format7(var_data, config):
	var_data_copy = copy.deepcopy(var_data)
	recom_dict = {}
	for i in PAN116_RECOMGene(config)["RECORDS"]:
		recom_dict[i["gene_symbol"]] = i["disease"]
	recom_gene = [i for i in recom_dict.keys()]
	var_list = [var for var in var_data_copy if judge_var(var, [3,4,5], [4,5])]
	result = []
	detect_gene = []
	for var in var_list:
		for gene in set(re.split(",", var["gene_symbol"])):
			if gene in recom_gene:
				detect_gene.append(gene)
				var["gene_symbol"] = gene
				var["disease"] = recom_dict[gene]
				var["clinic_num_s"] = S_level(var)
				result.append(var)	
	for gene in set(recom_gene) - set(detect_gene):
		result.append({
			"gene_symbol" : gene,
			"disease" : recom_dict[gene]
		})
	return sorted(result, key=lambda i:i["gene_symbol"])


# 合并梦晨代码-2024.01.11
# 格式7：用于TC30福建省立医院,2022年9月7日，嵇梦晨
def format7_FJSL(var_data):
	var_data_copy = copy.deepcopy(var_data)
	gene_list = ["BRAF", "TERT", "KRAS", "NRAS", "HRAS", "RET", "PAX8", "NTRK1", "NTRK2", "NTRK3", "ALK", "EGFR", "ERBB2", "MET", "ROS1"]
	var_list = [var for var in var_data_copy if var["clinic_num_s"] in [5, 4, 3] and var["var_origin"] != "germline"]
	result = {}
	if var_list:
		for i in var_list:
			# 2025.06.16-新增HD，代码兼容
			#var_info = "扩增" if i["bio_category"] == "Cnv" else "" \
			#	if i["bio_category"] in ["Sv", "PSeqRnaSv"] else i["hgvs_p_ZJZL"] if (i["hgvs_p"]!="p.?") else i["hgvs_c"]
			var_info = ""
			if i["bio_category"] == "PHd":
				if i["type"] == "HomoDel":
					var_info = "纯合缺失"
				elif i["type"] == "HeteDel":
					var_info = "杂合缺失"
				else:
					var_info = "未知变异类型！"
			elif i["bio_category"] == "Cnv":
				var_info = "扩增"
			elif i["bio_category"] in ["Sv", "PSeqRnaSv"]:
				var_info = ""
			elif i["bio_category"] == "Snvindel":
				if i["hgvs_p"] != "p.?":
					var_info = i["hgvs_p_ZJZL"]
				else:
					var_info = i["hgvs_c"]

			#var_info = "纯合缺失" if i["bio_category"] == "PHd" and i["type"] == "HomoDel" else "杂合缺失" if i["bio_category"] == "PHd" and i["type"] == "HeteDel" else "扩增" if i["bio_category"] == "Cnv" else "" \
			#	if i["bio_category"] in ["Sv", "PSeqRnaSv"] else i["hgvs_p_ZJZL"] if (i["hgvs_p"]!="p.?") else i["hgvs_c"]
			var_info_sv = i["five_prime_gene"]+":"+i["five_prime_cds"]+"-"+i["three_prime_gene"]+":"+i["three_prime_cds"] if i["bio_category"] == "Sv" else ""
			
			#freq = i["cn_mean"] if i["bio_category"] == "Cnv" else i["rnasv_reads"] if i["bio_category"] == "PSeqRnaSv" else i["freq_str"]
			# 2025.06.16-新增HD，代码兼容
			#freq = i["cn_mean"] if i["bio_category"] == "Cnv" else i["reads"] if i["bio_category"] == "PSeqRnaSv" else i["freq_str"]
			freq = i["cn_mean"] if i["bio_category"] == "Cnv" else i["reads"] if i["bio_category"] == "PSeqRnaSv" else i["freq_str"] if i["bio_category"] == "Snvindel" else ""
			#clinic_num_s = i["clinic_num_s"]
			clinic_num_s = S_level(i)
			if i["gene_symbol"] in gene_list:
				if i["gene_symbol"] not in result.keys():
					result.setdefault(i["gene_symbol"], [])
				result[i["gene_symbol"]].append({"var_info" : var_info, "freq" : freq, "var_info_sv" : var_info_sv, "clinic_num_s" : clinic_num_s})	
	return result
	
# 格式8：用于TC30福建协和,2022年9月7日，嵇梦晨
def format8_FJXH(var_data):
	var_data_copy = copy.deepcopy(var_data)
	gene_list = ["ALK", "BRAF", "EGFR", "FGFR2", "FGFR3", "KIT", "KRAS", "MET", "NRAS", "PDGFRA", "RET", "PIK3CA", "NTRK1", "NTRK2", "NTRK3"]
	var_list = [var for var in var_data_copy if var["clinic_num_s"] in [5, 4, 3] and var["var_origin"] != "germline"]
	result = {}
	if var_list:
		for i in var_list:
			# 2025.06.16-新增HD，代码兼容
			#var_info = "扩增" if i["bio_category"] == "Cnv" else "" \
			#	if i["bio_category"] in ["Sv", "PSeqRnaSv"] else i["hgvs_p_ZJZL"] if (i["hgvs_p"]!="p.?") else i["hgvs_c"]
			var_info = ""
			if i["bio_category"] == "PHd":
				if i["type"] == "HomoDel":
					var_info = "纯合缺失"
				elif i["type"] == "HeteDel":
					var_info = "杂合缺失"
				else:
					var_info = "未知变异类型！"
			elif i["bio_category"] == "Cnv":
				var_info = "扩增"
			elif i["bio_category"] in ["Sv", "PSeqRnaSv"]:
				var_info = ""
			elif i["bio_category"] == "Snvindel":
				if i["hgvs_p"] != "p.?":
					var_info = i["hgvs_p_ZJZL"]
				else:
					var_info = i["hgvs_c"]
			#var_info = "纯合缺失" if i["bio_category"] == "PHd" and i["type"] == "HomoDel" else "杂合缺失" if i["bio_category"] == "PHd" and i["type"] == "HeteDel" else "扩增" if i["bio_category"] == "Cnv" else "" \
			#	if i["bio_category"] in ["Sv", "PSeqRnaSv"] else i["hgvs_p_ZJZL"] if (i["hgvs_p"]!="p.?") else i["hgvs_c"]
			var_info_sv = i["five_prime_gene"]+":"+i["five_prime_cds"]+"-"+i["three_prime_gene"]+":"+i["three_prime_cds"] if i["bio_category"] == "Sv" else ""
			
			#freq = i["cn_mean"] if i["bio_category"] == "Cnv" else i["rnasv_reads"] if i["bio_category"] == "PSeqRnaSv" else i["freq_str"]
			# 2025.06.16-新增HD，代码兼容
			#freq = i["cn_mean"] if i["bio_category"] == "Cnv" else i["reads"] if i["bio_category"] == "PSeqRnaSv" else i["freq_str"]
			freq = i["cn_mean"] if i["bio_category"] == "Cnv" else i["reads"] if i["bio_category"] == "PSeqRnaSv" else i["freq_str"] if i["bio_category"] == "Snvindel" else ""
			#clinic_num_s = i["clinic_num_s"]
			clinic_num_s = S_level(i)
			# 变异类型 del ins dup delins
			# 2025.06.16-新增HD
			#bio_category = '拷贝数变异' if i['bio_category'] == 'Cnv' else '融合' if i['bio_category'] in ['Sv', 'PSeqRnaSv'] else '点突变' if '>' in i['hgvs_c'] else '插入/缺失' if 'delins' in i['hgvs_c'] \
			#	else '缺失' if 'del' in i['hgvs_c'] else '插入' if ('ins' in i['hgvs_c'] or 'dup' in i['hgvs_c']) else ''
			bio_category = ""
			if i["bio_category"] == "PHd":
				if i["type"] == "HomoDel":
					bio_category = "纯合缺失"
				elif i["type"] == "HeteDel":
					bio_category = "杂合缺失"
				else:
					bio_category = "未知变异类型！"
			elif i["bio_category"] == "Cnv":
				bio_category = "拷贝数变异"
			elif i['bio_category'] in ['Sv', 'PSeqRnaSv']:
				bio_category = "融合"
			elif i["bio_category"] == "Snvindel":
				bio_category = "点突变" if '>' in i['hgvs_c'] else '插入/缺失' if 'delins' in i['hgvs_c'] \
								else '缺失' if 'del' in i['hgvs_c'] else '插入' if ('ins' in i['hgvs_c'] or 'dup' in i['hgvs_c']) else ''

			#bio_category = '纯合缺失' if i['bio_category'] == 'PHd' and i["type"] == "HomoDel" else "杂合缺失" if i['bio_category'] == 'PHd' and i["type"] == "HeteDel" else '拷贝数变异' if i['bio_category'] == 'Cnv' else '融合' if i['bio_category'] in ['Sv', 'PSeqRnaSv'] else '点突变' if '>' in i['hgvs_c'] else '插入/缺失' if 'delins' in i['hgvs_c'] \
			#	else '缺失' if 'del' in i['hgvs_c'] else '插入' if ('ins' in i['hgvs_c'] or 'dup' in i['hgvs_c']) else ''
			if i["gene_symbol"] in gene_list:
				if i["gene_symbol"] not in result.keys():
					result.setdefault(i["gene_symbol"], [])
				result[i["gene_symbol"]].append({"var_info" : var_info, "freq" : freq, "var_info_sv" : var_info_sv, "clinic_num_s" : clinic_num_s, 'bio_category' : bio_category})

	return result
# 合并梦晨代码结束-2024.01.11

# 格式9：用于PAN116单独组织-福建省立-仅展示体细胞I/II类变异
def format9(var_data, config):
	var_data_copy = copy.deepcopy(var_data)
	recom_dict = {}
	for i in PAN116_RECOMGene(config)["RECORDS"]:
		recom_dict[i["gene_symbol"]] = i["disease"]
	recom_gene = [i for i in recom_dict.keys()]
	#var_list = [var for var in var_data_copy if judge_var(var, [4,5], [4,5])]
	var_list = [var for var in var_data_copy if var["var_origin"] != "germline" and var["clinic_num_s"] in [4,5] and judgeRegimen(var)]
	result = []
	detect_gene = []
	for var in var_list:
		for gene in set(re.split(",", var["gene_symbol"])):
			if gene in recom_gene:
				detect_gene.append(gene)
				var["gene_symbol"] = gene
				var["disease"] = recom_dict[gene]
				var["clinic_num_s"] = S_level(var)
				result.append(var)
	for gene in set(recom_gene) - set(detect_gene):
		result.append({
			"gene_symbol" : gene,
			"disease" : recom_dict[gene]
		})
	return sorted(result, key=lambda i:i["gene_symbol"])

# 2025.03.03-新增云南肿瘤CP40-与癌种/实体瘤相关、I/II类变异结果
# 展示基因、检测变异类型和结果，区分癌种，包含I/II类变异
def format_ynzl_cp40(var_data, tumor_list, config):
	var_data_copy = copy.deepcopy(var_data)
	recom_list = CP40_RECOMGene(config)["RECORDS"]
	gene_type = CP40Gene(config)
	# 将基因对应所有癌种合并起来
	gene_tumor = {}
	for i in recom_list:
		if i["gene_symbol"] not in gene_tumor.keys():
			gene_tumor.setdefault(i["gene_symbol"], [])
		gene_tumor[i["gene_symbol"]].append(i["disease"])
	# 筛选出对应癌种的NCCN推荐基因，其中均包含实体瘤
	tumor_recom_list = [i for i in recom_list if i["disease"] in tumor_list or i["disease"] == "实体瘤"]
	recom_dict = {}
	for i in tumor_recom_list:
		if i["gene_symbol"] not in recom_dict.keys():
			recom_dict[i["gene_symbol"]] = {
				"transcript" : i["transcript_refseq"],
				"disease" : i["disease"],
				"var_type" : i["var_type"],
				"gene_tumor" : "、".join(gene_tumor.get(i["gene_symbol"]))
			}

	result = []
	recom_gene = [i for i in recom_dict.keys()]
	detect_gene = []
	# 返回检测到变异的指南推荐基因
	var_list = [var for var in var_data_copy if var["var_origin"] != "germline" and var["clinic_num_s"] in [4,5] and judgeRegimen(var)]
	for var in var_list:
		for gene in set(re.split(",", var["gene_symbol"])):
			if gene in recom_dict.keys():
				detect_gene.append(gene)
				var["var_type"] = gene_type.get(gene, "-")
				result.append(var)
	# 返回未检测到变异的指南推荐基因-仅result
	for gene in set(recom_gene) - set(detect_gene):
		result.append({
			"gene_symbol" : gene,
			"var_type" : gene_type.get(gene, "-")
		})

	return sorted(result, key=lambda i:i["gene_symbol"])