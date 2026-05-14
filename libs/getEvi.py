#-*- coding:gbk -*-
import re
from libs.getConfig import senseStran
from libs.getPinyin import topinyin
from libs.getInterRef import getRef_from_inter
import copy
import itertools
from collections import defaultdict
from datetime import datetime

'''
Discription	

	将变异对应的治疗方案、证据描述等转化为适合填充的格式。 

'''

# 相同证据描述的合并治疗方案展示
def merge_Predictive_evi(datainfo):
	# 2024.12.09-结果新增refer_agency和sense_rule
	# 2025.08.02-结果新增evi_adaptation_disease_cn
	tmp_dict = {}
	for evi in datainfo:
		tmp_dict.setdefault(evi["evi_interpretation"], [])
		tmp_dict[evi["evi_interpretation"]].append(
			{
				"regimen_name" : evi["regimen_name"],
				"evi_conclusion_simple" : evi["evi_conclusion_simple"],
				"clinical_significance_cn" : evi["clinical_significance_cn"],
				"regimen_name_py" : evi["regimen_name_py"],
				"evi_conclusion" : evi["evi_conclusion"],
				"refer_agency" : evi["refer_agency"] if evi["refer_agency"] else "",
				"sense_rule" : evi["sense_rule"],
				"evi_adaptation_disease_cn" : evi["evi_adaptation_disease_cn"]
			}
		)

	merge_result = []
	for k, v  in tmp_dict.items():
		merge_result.append(
			{
				"regimen_name" : "、".join([i["regimen_name"] for i in v]),
				# 2024.12.09-evi_conclusion取最高等级
				#"evi_conclusion_simple" : "/".join([i["evi_conclusion_simple"] for i in v]),
				"evi_conclusion_simple" : "A" if "A" in [i["evi_conclusion_simple"] for i in v] else \
										  "B" if "B" in [i["evi_conclusion_simple"] for i in v] else \
										  "C" if "C" in [i["evi_conclusion_simple"] for i in v] else \
										  "D",
				"clinical_significance_cn" : "/".join([i["clinical_significance_cn"] for i in v]),
				"regimen_name_py" : "/".join([i["regimen_name_py"] for i in v]),
				"evi_interpretation" : k,
				"evi_conclusion" : "/".join([i["evi_conclusion"] for i in v]),
				"refer_agency" : "、".join([i["refer_agency"] for i in v]),
				"sense_rule" : "0" if "0" in [i["sense_rule"] for i in v] else "1",
				"evi_adaptation_disease_cn" : [i["evi_adaptation_disease_cn"] for i in v][0] if [i["evi_adaptation_disease_cn"] for i in v] else "",
				# 2025.08.25：华西HRD需要展示敏感/耐药，这边默认展示第一个值
				"clinical_significance_cn_first" : [i["clinical_significance_cn"] for i in v][0]
				# 2025.08.25：新增结束
			}
		)

	return merge_result

def varRegimen(jsonDict, evi_sum, config, var):
	data = {}
	data["refer_evi"] = []
	data["evi_split"] = {}
	data["refer_evi_risk"] = []

	# 1. 新增字段、证据描述基础处理、排序等
	for evi in evi_sum:
		## 新增字段-用于报告展示-1.临床意义转化为中文
		evi["clinical_significance_cn"] = senseStran(config).get(evi["clinical_significance"], evi["clinical_significance"])
		## 新增字段-用于报告展示-2.等级由A0、A1等转化为A
		evi["evi_conclusion_simple"] = evi["evi_conclusion"][0] if evi["evi_conclusion"] else ""
		## 新增字段-用于报告展示-3.适用浙肿，治疗方案获批机构，仅保留FDA、MNPA、EMA
		evi["regimen_refer_agency_ZJZL"] = "/".join(list(set(["FDA","NMPA","EMA"]) & set(re.split(",", evi["regimen_refer_agency"])))) \
										   if "regimen_refer_agency" in evi.keys() and evi["regimen_refer_agency"] \
										   else "-"
		
		## 新增字段-用于排序-1. 治疗方案转化为拼音，便于后续排序
		evi["regimen_name_py"] = topinyin(evi["regimen_name"]) if evi["regimen_name"] else "0"
		## 新增字段-用于排序-2. 用于对敏感和耐药进行排序，敏感前耐药后
		evi["sense_rule"] = "0" if re.search("Sensitive", evi["clinical_significance"]) else \
							"1" if re.search("Resistant", evi["clinical_significance"]) else \
							evi["clinical_significance"]
		
		# 其他功能-1.证据描述去掉末尾空格
		evi["evi_interpretation"] = evi["evi_interpretation"].strip() if evi["evi_interpretation"] else ""
		
	# 对证据进行排序，按治疗方案等级、敏感/耐药、治疗方案拼音（字母统一大写，便于排序）排序
	#evi_sum = sorted(evi_sum, key = lambda i:(i["evi_conclusion_simple"], i["sense_rule"], i["regimen_name_py"].upper()))
	# RET 融合，用药排序更新-由之前的卡博替尼、普拉替尼、赛普替尼更新为普拉提尼、赛普替尼、卡博替尼-2023.05.11
	drug_list_rule = {
		"普拉替尼" : 0, 
		"塞普替尼" : 1,
		"卡博替尼" : 2
		}
	# 2024.08.09-新增ERBB2 snv/cnv药物排序，德曲妥珠单抗同等级的放在最前面
	erbb2_drug_list_rule = {"德曲妥珠单抗" : 0}
	# 2024.08.09-新增完成
	if "bio_category" in var.keys() and var["bio_category"] in ["Sv", "PSeqRnaSv"] and (var["five_prime_gene"] == "RET" or var["three_prime_gene"] == "RET"):
		evi_sum = sorted(evi_sum, key = lambda i:(i["evi_conclusion_simple"], i["sense_rule"], drug_list_rule.get(i["regimen_name"], 3) , i["regimen_name_py"].upper()))
	# 2024.08.09-新增ERBB2 snv/cnv药物排序，德曲妥珠单抗同等级的放在最前面
	elif "bio_category" in var.keys() and var["bio_category"] in ["Snvindel", "Cnv"] and var["gene_symbol"] == "ERBB2":
		evi_sum = sorted(evi_sum, key = lambda i:(i["evi_conclusion_simple"], i["sense_rule"], erbb2_drug_list_rule.get(i["regimen_name"], 1) , i["regimen_name_py"].upper()))
	# 2024.08.09-新增完成
	else:
		evi_sum = sorted(evi_sum, key = lambda i:(i["evi_conclusion_simple"], i["sense_rule"], i["regimen_name_py"].upper()))
	# 排序更新完成-2023.05.11

	# 2. 对排序后的证据进行参考文献提取、证据分类、治疗证据合并等
	for evi in evi_sum:
		## 获取参考文献-证据描述
		data["refer_evi"].extend(getRef_from_inter(jsonDict, evi["evi_interpretation"]))

		## 获取参考文献-遗传风险
		if evi["evidence_type"] == "Predisposing":
			data["refer_evi_risk"].extend(getRef_from_inter(jsonDict, evi["evi_interpretation"]))
		
		## 创建evi_split，将证据进行分类，治疗、辅助诊断、预后、遗传风险
		if evi["evidence_type"] not in data["evi_split"].keys():
			data["evi_split"].setdefault(evi["evidence_type"], [])
		data["evi_split"][evi["evidence_type"]].append(evi)

		## 目前仅有Predictive时才要合并相同证据，额外设置“Predictive_merge”字段，根据报告需求进行选用
		if "Predictive" in data["evi_split"].keys():
			data["evi_split"]["Predictive_merge"] = merge_Predictive_evi(data["evi_split"]["Predictive"])

	# 一些特殊需求
	## 1. 汇总所有证据
	data["regimen_evi_sum"] = evi_sum
	## 2. 汇总A级敏感治疗方案
	data["regimen_FDA_S"] = [
		{
			"regimen_name" : var["regimen_name"], 
			"evi_conclusion_simple" : var["evi_conclusion_simple"]
			} 
		for var in evi_sum if re.search("Sensitive",var["clinical_significance"]) and var["evi_conclusion_simple"] == "A"
		]
	## 3. 汇总非A级敏感治疗方案
	data["regimen_noFDA_S"] = [
		{
			"regimen_name" : var["regimen_name"], 
			"evi_conclusion_simple" : var["evi_conclusion_simple"]
			} 
		for var in evi_sum if re.search("Sensitive",var["clinical_significance"]) and var["evi_conclusion_simple"] != "A"
		]

	## 4. 汇总敏感治疗方案-新增证据描述-2023.08.15
	data["regimen_S"] = [
		{
			"regimen_name" : var["regimen_name"], 
			"evi_conclusion_simple" : var["evi_conclusion_simple"],
			"evi_interpretation" : var["evi_interpretation"]
			} 
		for var in evi_sum if re.search("Sensitive",var["clinical_significance"])
		]

	## 5. 汇总耐药治疗方案-新增证据描述-2023.08.15
	data["regimen_R"] = [
		{
			"regimen_name" : var["regimen_name"], 
			"evi_conclusion_simple" : var["evi_conclusion_simple"],
			"evi_interpretation" : var["evi_interpretation"]
			} 
		for var in evi_sum if re.search("Resistant",var["clinical_significance"])
		]

	## 6. 敏感治疗方案列表转化为字符串
	data["regimen_S_str"] = "、".join([i["regimen_name"] for i in data["regimen_S"]])

	## 7. 耐药治疗方案列表转化为字符串
	data["regimen_R_str"] = "、".join([i["regimen_name"] for i in data["regimen_R"]])

	# 孙逸仙：药物拆分
	# I类 FDA/NMPA/NCCN药物：A0/A1/A2/C3 
	# I类 临床试验药物：B1/B2/B3/C1/C2/C4/D1/D2/D3/D4/D5/D6
	# II类 FDA/NMPA批准在其他癌种的药物： C3
	# II类 临床试验药物： C1/C2/C4/D1/D2/D3/D4/D5/D6
	# 可总结为两个字段，获批和临床试验药物，I/II类可根据变异等级来判定，应该不会混淆
	SYX_regimen_appr = [evi for evi in evi_sum if evi["evi_conclusion"] in ["A0","A1","A2","C3"]]
	# 2026.01.29:D等级新增D7、D8，这边加一下
	SYX_regimen_clinic = [evi for evi in evi_sum if evi["evi_conclusion"] in ["B1","B2","B3","C1","C2","C4","D1","D2","D3","D4","D5","D6", "D7", "D8"]]
	#根据治疗方案、预后、辅助诊断、风险等拆分
	data["SYX_regiman_appr"] = {}
	for evi in SYX_regimen_appr:
		if evi["evidence_type"] not in  data["SYX_regiman_appr"].keys():
			data["SYX_regiman_appr"].setdefault(evi["evidence_type"], [])
		data["SYX_regiman_appr"][evi["evidence_type"]].append(evi)
		
	data["SYX_regimen_clinic"] = {}
	for evi in SYX_regimen_clinic:
		if evi["evidence_type"] not in  data["SYX_regimen_clinic"].keys():
			data["SYX_regimen_clinic"].setdefault(evi["evidence_type"], [])
		data["SYX_regimen_clinic"][evi["evidence_type"]].append(evi)

	## 8. 安徽省立LC10 -2025.06.27-嵇梦晨
	#regimen_list = getRegimen_Approval_AHSL(evi_sum)
	#data["regimen_list_AHSL"] = Regimen_AHSL(regimen_list, evi_sum)
	#data["regimen_inter_AHSL"] = Regimen_inter_AHSL(regimen_list, evi_sum)	
	# 2025.06.27-新增完成

	return data

# 2025.06.27-嵇梦晨-安徽省立
def preprocess_AHSL(regimen_list, datainfo):
	"""
	加一个判断是否NMPA获批 “NMPA”:"yes"
	排序规则 A类：NMPA/FDA > NMPA > FDA > CSCO/NCCN > CSCO > NCCN
	无A类：B > C > D, regimen_name_py,临床试验、回顾性分析、案例报道、临床前证据
	C4 D4证据不展示，报告系统无此功能，先用脚本过滤
	B、C、D类证据，regime_name相同，evi_conclusion B1&B2 B2&B3 C1&C2 D1&D2 D5&D6同时出现，删除这条治疗方案
	"""
	datainfo = [evi for evi in datainfo if evi["evi_conclusion"] not in ["C4", "D4"]] if datainfo else []
	conflict_groups = [
        {'B1', 'B2'},  # B1和B2同时出现会产生冲突
		{'B2', 'B3'},  # B2和B3同时出现会产生冲突
		{'C1', 'C2'},  # C1和C2同时出现会产生冲突
        {'D1', 'D2'},  # D1和D2同时出现会产生冲突
        {'D5', 'D6'}   # D5和D6同时出现会产生冲突
    ]
	regimen_groups = defaultdict(list)
	for evi in datainfo:
		regimen_groups[evi['regimen_name']].append(evi)
	
	for regimen_name, items in regimen_groups.items():
		conclusions = set(evi['evi_conclusion'] for evi in items)

		for group in conflict_groups:
			if group.issubset(conclusions):
				datainfo = [evi for evi in datainfo if evi['regimen_name'] != regimen_name]
		if len(items) > 1:
			# NMPA FDA获批的药物，合并适应症evi_adaptation_disease_cn, 在datainfo中只保留一个
			if any(re.search("NMPA|FDA", evi["refer_agency"]) for evi in items):
				merged_adaptation = "\n".join({evi["evi_adaptation_disease_cn"] for evi in items if evi["evi_adaptation_disease_cn"]})
				for evi in datainfo:
					if evi["regimen_name"] == regimen_name and re.search("NMPA|FDA", evi["refer_agency"]):
						evi["evi_adaptation_disease_cn"] = merged_adaptation

	if not datainfo:
		return datainfo
	for evi in datainfo:
		evi["refer_agency_2"] = regimen_list.get(evi["regimen_name"], "-")
#		evi["refer_agency_sort"] = 0 if re.search("NMPA", evi["refer_agency"]) else 1 if re.search("FDA", evi["refer_agency"]) else 2 if re.search("CSCO", evi["refer_agency"]) \
#			else 3 if re.search("NCCN", evi["refer_agency"]) else 4 if re.search("Clinical", evi["evidence_level"]) else 5 if re.search("Retrospecific", evi["evidence_level"]) \
#			else 6 if re.search("Case", evi["evidence_level"]) else 7 if re.search("Preclinical", evi["evidence_level"]) else 8
		evi["refer_agency_sort"] = 0 if re.search("NMPA/FDA", evi["refer_agency_2"]) else 1 if re.search("NMPA", evi["refer_agency_2"]) else \
			2 if re.search("FDA", evi["refer_agency_2"]) else 3 if re.search("CSCO/NCCN", evi["refer_agency_2"]) else \
			4 if re.search("CSCO", evi["refer_agency_2"]) else 5 if re.search("NCCN", evi["refer_agency_2"]) else \
			6 if re.search("Clinical", evi["evidence_level"]) else 7 if re.search("Retrospecific", evi["evidence_level"]) \
			else 8 if re.search("Case", evi["evidence_level"]) else 9 if re.search("Preclinical", evi["evidence_level"]) else 10
		#evi["NMPA"] = "yes" if re.search("NMPA", evi["refer_agency"]) and evi["evi_conclusion_simple"] == "A" else "no"
		evi["NMPA"] = "yes" if re.search("NMPA", evi["refer_agency_2"]) else "yes" if re.search("NMPA", evi["regimen_refer_agency"]) else "no"
	datainfo = sorted(datainfo, key= lambda i:(i["evi_conclusion_simple"], i["refer_agency_sort"], i["regimen_name_py"].upper()))

	return datainfo

def Regimen_AHSL(regimen_list, datainfo):
	datainfo = preprocess_AHSL(regimen_list, datainfo)
	if not datainfo:
		return {"sensi" : [], "resist" : []}
	
	sensi_list = [evi for evi in datainfo if "Sensitive" in evi["clinical_significance"]]
	resist_list = [evi for evi in datainfo if "Resistant" in evi["clinical_significance"]]

	# 去重 2025.06.11
	def deduplicate_evidence(evidence_list):
		seen = set()
		deduped = []
		for evi in evidence_list:
			# 创建唯一标识符
			identifier = (evi["regimen_name"], evi["evi_conclusion_simple"], evi["refer_agency_2"])
			if identifier not in seen:
				seen.add(identifier)
				deduped.append(evi)
		keys = ["regimen_name", "evi_conclusion_simple", "refer_agency_2", "NMPA"]
		for evi in deduped:
			evi = {key : evi[key] for key in keys}
		return deduped
	
	sensi_deduped = deduplicate_evidence(sensi_list)
	resist_deduped = deduplicate_evidence(resist_list)

	return {"sensi": sensi_deduped, "resist": resist_deduped}

def getRegimen_Approval_AHSL(evi_sum):
	'''
	evi_sum返回了所有证据，从evi_sum的refer_agency字段获取治疗方案的获批机构
	regimen_name, refer_agency
	'''
	apprlist = ["FDA", "NMPA", "NCCN", "CSCO"]
	appdict = {"NMPA" : 0, "FDA" : 1, "CSCO" : 2, "NCCN" : 3}
	regimen_groups = {}
	if evi_sum:
		for evi in evi_sum:
			regimen = evi["regimen_name"]
			agency = evi["refer_agency"]
			if agency in apprlist:
				if regimen not in regimen_groups:
					regimen_groups[regimen] = set()
				regimen_groups[regimen].add(agency)
	
	result = {}
	for regimen, agencies in regimen_groups.items():
		sorted_agencies = sorted(agencies, key=lambda x: appdict.get(x, 4))
		result[regimen] = "/".join(sorted_agencies)
	
	return result

def process_top_evi(evi_sum):
	"""
	最高等级返回所有证据
	B/C/D为最高等级时，同一个治疗方案取最新的一条
	"""
	if not evi_sum:
		return []
	
	result = []
	conclusion_priority = {
		'B1': 9, 'B2': 8, 'B3': 7, 
		'C1': 6, 'C2': 5, 'C3': 4,
		'D1': 3, 'D2': 3,  # D1和D2视为同一等级
		'D3': 2,
		'D5': 1, 'D6': 1   # D5和D6视为同一等级
	}
	regimen_groups = defaultdict(list)
	for evi in evi_sum:
		regimen_groups[evi['regimen_name']].append(evi)
	
	for regimen_name, items in regimen_groups.items():
		highest_priority = -1
		candidate_items = []
		for item in items:
			priority = conclusion_priority.get(item['evi_conclusion'], 0)
			if priority > highest_priority:
				highest_priority = priority
				candidate_items = [item]
			elif priority == highest_priority:
				candidate_items.append(item)
		
		if len(candidate_items) > 1:
			# 如果有多条证据，按发布时间排序，取最新的一条
			print(len(candidate_items))
			candidate_items.sort(key=lambda x: datetime.strptime(x.get('publish_time', '1900-01-01'), '%Y-%m-%d') if x.get('publish_time') else datetime.min, reverse=True)

		# 保留优先级最高且最近发布的证据
		result.append(candidate_items[0])

	return result

def merge_AHSL(evi_sum):
    if not evi_sum:
        return []

    interpretation_groups = defaultdict(list)
    for evi in evi_sum:
        key = evi.get('evi_interpretation', '')
        interpretation_groups[key].append(evi)
    
    merged_evi = []
    for interpretation, items in interpretation_groups.items():
        if len(items) == 1:
            merged_evi.append(items[0])
            continue
        
        base_item = items[0].copy()
        regimen_names = [item['regimen_name'] for item in items]
        base_item['regimen_name'] = "、".join(regimen_names)
        merged_evi.append(base_item)
    
    return merged_evi

def Regimen_inter_AHSL(regimen_list, datainfo):
	"""
	有A类药物的治疗方案，分为敏感和耐药，加一个判断是否NMPA获批 “NMPA”:"yes"
	无A类药物的治疗方案，分为敏感和耐药，细分为临床试验/回顾性分析、案例报道/临床前证据两类
	相同临床证据，合并
	排序规则 A类：FDA/NMPA > NMPA > CSCO/NCCN > CSCO > NCCN
	无A类：B > C > D, regimen_name_py
	C4 D4证据不展示，报告系统无此功能，先用脚本过滤
	最高等级为A：同一治疗方案展示所有证据，NMPA/FDA展示获批适应症，NCCN/CSCO展示证据，取最新的一条
	"""
	#------------------------------------------------------------------------------------------------------------------------------
	# 规则更新 2025.06.12
	# a)知识库返回最高等级的所有证据
	#     A类证据，NCCN&CSCO-regimen_name相同，取发布时间最近的一条，报告展示evi_interpretation内容
	#     B、C、D类证据，regime_name相同，取发布时间最近的一条，报告展示evi_interpretation内容
	# b)新增发布时间字段，暂定为publish_time:"2011-12-27"
	# c)evi_adaptation_disease_cn相同的药物，regimen_name “、”连接
	# d)evi_interpretation相同的药物，regimen_name “、”连接
	#----------------------------------------------------------------------------------------------
	# 2025年7月3日 jmc
	# 分类逻辑更改，先区分敏感和耐药，再根据不同的证据等级进行分类
	#------------------------------------------------------------------------------------------------------------------------------
	tmp_dict = {
		"sensi" : {"level_A" : [], "level_A2" : [], "level_Clinic" : [], "level_D" : [], "regimen_list" : []},
		"resist" : {"level_A" : [], "level_A2" : [], "level_Clinic" : [], "level_D" : [], "regimen_list" : []}
	}
	datainfo = preprocess_AHSL(regimen_list, datainfo)

	if not datainfo:
		return tmp_dict
	
	sensi_list = [evi for evi in datainfo if "Sensitive" in evi["clinical_significance"]]
	resist_list = [evi for evi in datainfo if "Resistant" in evi["clinical_significance"]]
	tmp_dict["sensi"]["level_A"] = [evi for evi in sensi_list if re.search("NMPA|FDA", evi["refer_agency"]) and evi["evi_conclusion_simple"] == "A"]
	tmp_dict["resist"]["level_A"] = [evi for evi in resist_list if re.search("NMPA|FDA", evi["refer_agency"]) and evi["evi_conclusion_simple"] == "A"]
	tmp_dict["sensi"]["level_A2"] = [evi for evi in sensi_list if re.search("CSCO|NCCN", evi["refer_agency"]) and evi["evi_conclusion_simple"] == "A"]
	tmp_dict["resist"]["level_A2"] = [evi for evi in resist_list if re.search("CSCO|NCCN", evi["refer_agency"]) and evi["evi_conclusion_simple"] == "A"]
	tmp_dict["sensi"]["level_Clinic"] = [evi for evi in sensi_list if re.search("Clinical|Retrospecific", evi["evidence_level"]) and evi["evi_conclusion_simple"] != "A"]
	tmp_dict["resist"]["level_Clinic"] = [evi for evi in resist_list if re.search("Clinical|Retrospecific", evi["evidence_level"]) and evi["evi_conclusion_simple"] != "A"]
	tmp_dict["sensi"]["level_D"] = [evi for evi in sensi_list if re.search("Case|Preclinical", evi["evidence_level"]) and evi["evi_conclusion_simple"] != "A"]
	tmp_dict["resist"]["level_D"] = [evi for evi in resist_list if re.search("Case|Preclinical", evi["evidence_level"]) and evi["evi_conclusion_simple"] != "A"]

	regimen_list_A = [{"name" : evi["regimen_name"], "NMPA" : evi["NMPA"]} for evi in tmp_dict["sensi"]["level_A"] + tmp_dict["sensi"]["level_A2"]]
	unique_dict = {(r["name"], r["NMPA"]): r for r in regimen_list_A}
	regimen_list_A = list(unique_dict.values())
	regimen_list_B = [{"name" : evi["regimen_name"], "NMPA" : evi["NMPA"]} for evi in tmp_dict["sensi"]["level_Clinic"] + tmp_dict["sensi"]["level_D"]]
	unique_dict = {(r["name"], r["NMPA"]): r for r in regimen_list_B}
	regimen_list_B = list(unique_dict.values())
	regimen_list_C = list(dict.fromkeys([evi["regimen_name"] for evi in tmp_dict["resist"]["level_A"] + tmp_dict["resist"]["level_A2"]]))
	regimen_list_D = list(dict.fromkeys([evi["regimen_name"] for evi in tmp_dict["resist"]["level_Clinic"] + tmp_dict["resist"]["level_D"]]))
	if regimen_list_A:
		tmp_dict["sensi"]["regimen_list"] = regimen_list_A
	elif regimen_list_B:
		tmp_dict["sensi"]["regimen_list"] = regimen_list_B
	else:
		tmp_dict["sensi"]["regimen_list"] = []
	if regimen_list_C:
		tmp_dict["resist"]["regimen_list"] = regimen_list_C
	elif regimen_list_D:
		tmp_dict["resist"]["regimen_list"] = regimen_list_D
	else:
		tmp_dict["resist"]["regimen_list"] = []

	if tmp_dict["sensi"]["level_A"]:
		adaptation_groups = defaultdict(list)
		for evi in tmp_dict["sensi"]["level_A"]:
			adaptation_groups[evi["evi_adaptation_disease_cn"]].append(evi)
		tmp_dict["sensi"]["level_A"] = []
		for adaptation, items in adaptation_groups.items():
			base_item = items[0].copy()
			regimen_names = set(evi["regimen_name"] for evi in items)
			base_item['regimen_name'] = "、".join(regimen_names)
			tmp_dict["sensi"]["level_A"].append(base_item)
		tmp_dict["sensi"]["level_A"] = sorted(tmp_dict["sensi"]["level_A"], key = lambda i:(i["refer_agency_sort"], i["regimen_name_py"].upper()))
	if tmp_dict["resist"]["level_A"]:
		adaptation_groups = defaultdict(list)
		for evi in tmp_dict["resist"]["level_A"]:
			adaptation_groups[evi["evi_adaptation_disease_cn"]].append(evi)
		tmp_dict["resist"]["level_A"] = []
		for adaptation, items in adaptation_groups.items():
			base_item = items[0].copy()
			regimen_names = set(evi["regimen_name"] for evi in items)
			base_item['regimen_name'] = "、".join(regimen_names)
			tmp_dict["resist"]["level_A"].append(base_item)
		tmp_dict["resist"]["level_A"] = sorted(tmp_dict["resist"]["level_A"], key = lambda i:(i["refer_agency_sort"], i["regimen_name_py"].upper()))
	
	tmp_dict["sensi"]["level_A2"] = sorted(tmp_dict["sensi"]["level_A2"], key = lambda i:(i["evi_conclusion_simple"], i["refer_agency_sort"], i["regimen_name_py"].upper()))
	tmp_dict["resist"]["level_A2"] = sorted(tmp_dict["resist"]["level_A2"], key = lambda i:(i["evi_conclusion_simple"], i["refer_agency_sort"], i["regimen_name_py"].upper()))
	tmp_dict["sensi"]["level_A2"] = merge_AHSL(tmp_dict["sensi"]["level_A2"])
	tmp_dict["resist"]["level_A2"] = merge_AHSL(tmp_dict["resist"]["level_A2"])
	tmp_dict["sensi"]["level_Clinic"] = process_top_evi(tmp_dict["sensi"]["level_Clinic"])
	tmp_dict["resist"]["level_Clinic"] = process_top_evi(tmp_dict["resist"]["level_Clinic"])
	tmp_dict["sensi"]["level_Clinic"] = sorted(tmp_dict["sensi"]["level_Clinic"], key = lambda i:(i["evi_conclusion_simple"], i["refer_agency_sort"], i["regimen_name_py"].upper()))
	tmp_dict["resist"]["level_Clinic"] = sorted(tmp_dict["resist"]["level_Clinic"], key = lambda i:(i["evi_conclusion_simple"], i["refer_agency_sort"], i["regimen_name_py"].upper()))
	tmp_dict["sensi"]["level_Clinic"] = merge_AHSL(tmp_dict["sensi"]["level_Clinic"])
	tmp_dict["resist"]["level_Clinic"] = merge_AHSL(tmp_dict["resist"]["level_Clinic"])
	tmp_dict["sensi"]["level_D"] = process_top_evi(tmp_dict["sensi"]["level_D"])
	tmp_dict["resist"]["level_D"] = process_top_evi(tmp_dict["resist"]["level_D"])
	tmp_dict["sensi"]["level_D"] = sorted(tmp_dict["sensi"]["level_D"], key = lambda i:(i["evi_conclusion_simple"], i["refer_agency_sort"], i["regimen_name_py"].upper()))
	tmp_dict["resist"]["level_D"] = sorted(tmp_dict["resist"]["level_D"], key = lambda i:(i["evi_conclusion_simple"], i["refer_agency_sort"], i["regimen_name_py"].upper()))
	tmp_dict["sensi"]["level_D"] = merge_AHSL(tmp_dict["sensi"]["level_D"])
	tmp_dict["resist"]["level_D"] = merge_AHSL(tmp_dict["resist"]["level_D"])
	
	return tmp_dict

# 2025.06.27-新增完成