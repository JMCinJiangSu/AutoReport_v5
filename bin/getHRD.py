#-*- coding:gbk -*-
import copy
from libs import listResultToDict
from libs.getEvi import varRegimen
from datetime import datetime

'''
Discription
	
	该脚本用来获取HRD结果。
	提取内容：
	1. 结果
	2. 治疗方案
	3. 合并BRCA突变的治疗方案（不包含解释） 

'''

def getHRD(jsonDict, BRCA_data, config):
	hrd_dict = listResultToDict.ListToDict(copy.deepcopy(jsonDict["hrd"]))
	if "evi_sum" in hrd_dict.keys() and hrd_dict["evi_sum"]:
		# hrd治疗方案转化
		hrd_dict["evi_sum"] = varRegimen(jsonDict, hrd_dict["evi_sum"], config, hrd_dict)
		# hrd治疗方案汇总
		# 2025.07.02-增加refer_agency和evidence_level
		# 2025.08.01-增加日期
		regimen_sum = [
			{
			"regimen_name":i["regimen_name"], 
			"evidence_type":i["evidence_type"], 
			"clinical_significance_cn":i["clinical_significance_cn"], 
			"evi_conclusion_simple":i["evi_conclusion_simple"],
			"regimen_name_py":i["regimen_name_py"],
			"refer_agency":i["refer_agency"] if "refer_agency" in i.keys() and i["refer_agency"] else "",
			"evidence_level":i["evidence_level"] if "evidence_level" in i.keys() and i["evidence_level"] else "",
			"publish_time":i["publish_time"] if "publish_time" in i.keys() and i["publish_time"] else "1900-01-01"
			} 
			for i in hrd_dict["evi_sum"]["regimen_evi_sum"] if i["evidence_type"] in ["Diagnostic","Predictive","Prognostic"]]
		# BRCA治疗方案汇总
		# 2025.07.02-增加refer_agency和evidence_level
		# 2025.08.01-增加日期
		if BRCA_data:
			for var in BRCA_data:
				regimen_sum += [
					{
						"regimen_name":i["regimen_name"], 
						"evidence_type":i["evidence_type"], 
						"clinical_significance_cn":i["clinical_significance_cn"], 
						"evi_conclusion_simple":i["evi_conclusion_simple"],
						"regimen_name_py":i["regimen_name_py"],
						"refer_agency":i["refer_agency"] if "refer_agency" in i.keys() and i["refer_agency"] else "",
						"evidence_level":i["evidence_level"] if "evidence_level" in i.keys() and i["evidence_level"] else "",
						"publish_time":i["publish_time"] if "publish_time" in i.keys() and i["publish_time"] else "1900-01-01"
						} for i in var["evi_sum"]["regimen_evi_sum"] if i["evidence_type"] in ["Diagnostic","Predictive","Prognostic"]]
		# hrd + BRCA治疗方案去重、排序

		
		# 更新去重规则-仅根据evidence_type、regimen_name、clinical_significance_cn和evi_conclusion_simple进行过滤-2025.08.01
		# 有多条时
		# 1. 等级最高的
		# 2. 保留日期最新的
		# 3. NMPA>FDA
		# 4. Clinical-phase IV > Clinical-phase III > Clinical-phase II > Clinical-phase I > Clinical-retrospective > Clinical-unknown phase > \
		#    Case report > Preclinical-in vivo > Preclinical-in vitro
		#regimen_sum_redup = []
		#for i in regimen_sum:
		#	if i not in regimen_sum_redup:
		#		regimen_sum_redup.append(i)	

		# 按日期、引用机构、研究类型排序		
		appr_rule = ["NMPA", "FDA"]
		clinic_rule = ["Clinical-phase IV", "Clinical-phase III", "Clinical-phase II", "Clinical-phase I", "Clinical-retrospective", \
				 	   "Clinical-unknown phase", "Case report", "Preclinical-in vivo", "Preclinical-in vitro"]
		regimen_sum_sort = copy.deepcopy(regimen_sum)
		for i in regimen_sum_sort:
			i["refer_agency_num"] = appr_rule.index(i["refer_agency"]) if i["refer_agency"] in appr_rule else 999
			i["evidence_level_num"] = clinic_rule.index(i["evidence_level"]) if i["evidence_level"] in clinic_rule else 999
		# 先不考虑耐药的排序情况（默认为敏感的）
		regimen_sum_sort = sorted(regimen_sum_sort, key = lambda i:(i["refer_agency_num"], i["evidence_level_num"]))
		regimen_sum_sort = sorted(regimen_sum_sort, key = lambda i:(datetime.strptime(i["publish_time"], "%Y-%m-%d")), reverse=True)
		regimen_sum_sort = sorted(regimen_sum_sort, key = lambda i:(i["evi_conclusion_simple"]))
		# 去重
		tmp_dict = {}
		for i in regimen_sum_sort:
			key = (i["evidence_type"], i["regimen_name"], i["clinical_significance_cn"], i["evi_conclusion_simple"])
			if key not in tmp_dict:
				tmp_dict.setdefault(key, [])
			tmp_dict[key].append(i)
		
		regimen_sum_redup = []
		for k, v in tmp_dict.items():
			regimen_sum_redup.append(v[0])

	
		regimen_sum_redup = sorted(regimen_sum_redup, key=lambda i:(i["evi_conclusion_simple"], i["clinical_significance_cn"], i["regimen_name_py"]))
		# 等级判断（判断依据包含hrd和BRCA的结果）
		regimen_level = [i["evi_conclusion_simple"] for i in regimen_sum_redup]
		hrd_dict["level_num"] = 5 if set(["A", "B"]) & set(regimen_level) else 4 if set(["C", "D"]) & set(regimen_level) else 3
		# 治疗方案分类展示
		hrd_dict["regimen"] = {}
		for regimen in regimen_sum_redup:
			if regimen["evidence_type"] not in hrd_dict["regimen"]:
				hrd_dict["regimen"].setdefault(regimen["evidence_type"], [])
			hrd_dict["regimen"][regimen["evidence_type"]].append(regimen)

	return hrd_dict